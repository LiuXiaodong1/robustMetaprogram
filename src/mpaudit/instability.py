"""Pilot P1a: how much do 3CA-style MPs move when arbitrary choices change?

Perturbs one choice at a time (NMF seed, rank grid, filter/clustering thresholds, gene universe,
leave-one-study-out, cell bootstrap) and matches each perturbed MP set back to the reference run.
Fits are cached, so threshold and leave-study-out perturbations cost no NMF. Each perturbation row
(``instability_runs.csv``) carries the run's NMF fit counts (``n_fits``, ``n_capped_fits``,
``n_failed_fits``; see :mod:`mpaudit.gavish`).

With few tumours per study, leave-one-study-out can remove enough programs that an MP no longer has
more than ``min_group_size`` founder partners; low leave-one-study-out recall on small cohorts (e.g.
the 10-tumour toy of :mod:`mpaudit.sim`, where a study holds 3-4 tumours and EMT, IFN and the CNA
are planted in 5, 6 and 7 tumours) is then mechanical, not evidence of instability.

Method variants (kind ``"method_variant"``; currently ``uncentred``, NMF without gene centring)
change the method rather than an arbitrary setting, and an MP that does not survive them is not
unstable. They are reported in their own column (``recovery_method_variant``) and summary rows,
and are left out of ``recovery_all``, ``median_best_jaccard`` and the gene-frequency table.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .compare import match_sets, summarize_match
from .data import Keys, sample_map
from .gavish import GavishParams, GavishResult, NMFCache, assemble, check_keys, compute_programs

METHOD_VARIANT = "method_variant"


@dataclass
class Perturbation:
    name: str
    kind: str
    params: dict = field(default_factory=dict)
    drop_samples: tuple[str, ...] = ()
    cell_fraction: float = 1.0
    subset_seed: int = 0


def default_perturbations(
    base: GavishParams,
    *,
    study_samples: dict[str, list[str]] | None = None,
    seeds: tuple[int, ...] = (1, 2, 3),
    rank_grids: tuple[tuple[int, ...], ...] = ((3, 4, 5, 6, 7, 8), (5, 6, 7, 8, 9, 10), (4, 5, 6)),
    threshold_steps: dict | None = None,
    universes: tuple[tuple[str, int], ...] = (("top_mean", 5000), ("hvg", 2000), ("hvg", 5000)),
    method_variants: bool = True,
    n_boot: int = 3,
    boot_fraction: float = 0.8,
) -> list[Perturbation]:
    """A one-factor-at-a-time grid around ``base``; ``method_variants`` adds the method variants
    (reported separately, see the module docstring)."""
    out: list[Perturbation] = []
    for s in seeds:
        out.append(Perturbation(f"seed={s}", "seed", {"seed": s}))
    for rg in rank_grids:
        if tuple(rg) != tuple(base.ranks):
            out.append(Perturbation(f"ranks={rg[0]}-{rg[-1]}", "ranks", {"ranks": tuple(rg)}))
    steps = threshold_steps or {
        "intra_min": (30, 40), "inter_min": (8, 12), "intra_max": (5, 15),
        "min_intersect_cluster": (8, 12), "min_group_size": (3, 7),
    }
    for k, vals in steps.items():
        for v in vals:
            if getattr(base, k) != v:
                out.append(Perturbation(f"{k}={v}", "threshold", {k: v}))
    for u, n in universes:
        if (u, n) != (base.gene_universe, base.n_universe):
            out.append(Perturbation(f"universe={u}{n}", "universe", {"gene_universe": u, "n_universe": n}))
    if method_variants and base.center:
        out.append(Perturbation("uncentred", METHOD_VARIANT, {"center": False}))
    for study, samples in (study_samples or {}).items():
        out.append(Perturbation(f"leave_out={study}", "loso", drop_samples=tuple(samples)))
    for b in range(n_boot):
        out.append(Perturbation(f"cells{int(boot_fraction * 100)}%#{b + 1}", "bootstrap",
                                cell_fraction=boot_fraction, subset_seed=b + 1))
    return out


@dataclass
class InstabilityResult:
    reference: GavishResult
    runs: dict[str, dict[str, list[str]]]
    summary: pd.DataFrame
    per_mp: pd.DataFrame
    genes: pd.DataFrame


def run_instability(
    adata,
    keys: Keys,
    base: GavishParams,
    perturbations: list[Perturbation],
    *,
    counts_layer: str | None = None,
    cache: NMFCache | None = None,
    recurrence_key: str | None = None,
    min_jaccard: float = 0.3,
    n_jobs: int = 1,
    reference: GavishResult | None = None,
    log=print,
) -> InstabilityResult:
    """Run the reference (unless given) and every perturbation, and match each perturbed MP set
    to the reference.

    Every perturbed :class:`GavishParams` is built and checked (:func:`mpaudit.gavish.check_keys`)
    before the first fit, so a bad perturbation fails before any NMF rather than after the
    reference and the perturbations before it. Each row of ``summary`` carries the NMF fit
    diagnostics of its run (``n_fits``, ``n_capped_fits`` stopped at ``max_iter``,
    ``n_failed_fits`` stopped by the restart limit; fits served from the cache are counted again).
    """
    perturbed = [base.replace(**p.params) for p in perturbations]
    for params in [base, *perturbed]:
        check_keys(keys, params, adata, recurrence_key=recurrence_key)
    cache = cache or NMFCache()
    unit = sample_map(adata, keys, recurrence_key) if recurrence_key else None
    study = sample_map(adata, keys, keys.study) if keys.study else None

    def run(params, p: Perturbation | None = None) -> GavishResult:
        spg = compute_programs(adata, keys, params, counts_layer=counts_layer, cache=cache,
                               cell_fraction=p.cell_fraction if p else 1.0,
                               subset_seed=p.subset_seed if p else 0, n_jobs=n_jobs)
        return assemble(spg, params, unit_of_sample=unit, study_of_sample=study,
                        drop_samples=frozenset(p.drop_samples) if p else frozenset())

    ref = reference or run(base)
    ref_sets = ref.mp_sets()
    log(f"[instability] reference: {len(ref.robust)} robust programs, {len(ref_sets)} MPs")

    runs, rows, matches = {}, [], {}
    for p, params in zip(perturbations, perturbed):
        res = run(params, p)
        alt = res.mp_sets()
        runs[p.name] = alt
        m = match_sets(ref_sets, alt, min_jaccard)
        matches[p.name] = (m, alt)
        rows.append({"perturbation": p.name, "kind": p.kind, "n_robust": len(res.robust),
                     **summarize_match(m, len(alt)), "n_fits": res.n_fits, "n_capped_fits": res.n_capped_fits,
                     "n_failed_fits": res.n_failed_fits})
        log(f"[instability] {p.name}: {len(alt)} MPs, recall {rows[-1]['recall']:.2f}")
    summary = pd.DataFrame(rows)

    per_mp_rows, gene_rows = [], []
    kinds = sorted({p.kind for p in perturbations} | {METHOD_VARIANT})
    stability = [p.name for p in perturbations if p.kind != METHOD_VARIANT]
    for mp, genes in ref_sets.items():
        rec = {"mp": mp, "n_genes": len(genes)}
        for kind in kinds:
            names = [p.name for p in perturbations if p.kind == kind]
            hits = [bool(matches[n][0].set_index("ref").loc[mp, "matched"]) for n in names]
            rec[f"recovery_{kind}"] = float(np.mean(hits)) if hits else np.nan
        all_hits, jacs = [], []
        gene_hits = {g: 0 for g in genes}
        for name in stability:
            m, alt = matches[name]
            row = m.set_index("ref").loc[mp]
            all_hits.append(bool(row["matched"]))
            jacs.append(float(row["best_jaccard"]))
            if row["matched"]:
                for g in set(alt[row["alt"]]) & set(genes):
                    gene_hits[g] += 1
        rec["recovery_all"] = float(np.mean(all_hits)) if all_hits else np.nan
        rec["median_best_jaccard"] = float(np.median(jacs)) if jacs else np.nan
        per_mp_rows.append(rec)
        n = max(len(stability), 1)
        gene_rows += [{"mp": mp, "gene": g, "frequency": c / n} for g, c in gene_hits.items()]
    return InstabilityResult(reference=ref, runs=runs, summary=summary,
                             per_mp=pd.DataFrame(per_mp_rows), genes=pd.DataFrame(gene_rows))
