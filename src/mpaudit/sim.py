"""Toy tumour cohort with known ground truth, for tests and demonstrations.

Planted in malignant cells:

* biology: ``cycle`` and ``hypoxia`` (all tumours), ``EMT`` (half the tumours; higher
  prevalence under the fixed protocol, i.e. selection against fragile cells in fresh samples);
* TME-shared biology: ``IFN`` (active in malignant *and* TME cells of the same samples);
* private programs: two per tumour;
* nuisance: dissociation ``stress`` (immediate-early/heat-shock genes; sample-level intensity,
  shared by malignant and TME cells, absent under the fixed protocol), a recurrent ``CNA``
  (3x focal amplification of ~40 genes on one arm in a subclone of 70% of tumours;
  a 2x gain at these expression levels was too weak for 3CA-style NMF to pick up), ambient RNA from TME cells
  (per-cell fraction varies), and malignant-TME doublets.

Outcomes of an audit on this cohort are true by construction. Use it to check that code runs
and that each axis responds as designed, not as evidence about real data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .typer import IEG_HSP_CORE

T_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GZMA"]
MYELOID_MARKERS = ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "TYROBP", "FCER1G", "AIF1"]
FIB_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]
RIBO = [f"RPL{i}" for i in range(3, 13)] + [f"RPS{i}" for i in range(2, 12)]
MITO = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ND4", "MT-CYB", "MT-ATP6"]


def _counts(rate: np.ndarray, depth: np.ndarray, theta: float, rng: np.random.Generator) -> np.ndarray:
    rate = rate / rate.sum(axis=1, keepdims=True)
    mu = depth[:, None] * rate
    lam = rng.gamma(theta, mu / theta)
    return rng.poisson(lam)


def simulate_cohort(
    n_tumours: int = 10,
    n_twin: int = 4,
    n_region: int = 2,
    n_mal: int = 200,
    n_tme_per_type: int = 40,
    n_genes: int = 1500,
    program_size: int = 40,
    n_studies: int = 3,
    seed: int = 0,
):
    """Return ``(adata, truth)``. ``truth`` holds planted gene sets, their kinds and a cytoband table."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    named = list(dict.fromkeys(IEG_HSP_CORE + T_MARKERS + MYELOID_MARKERS + FIB_MARKERS + RIBO + MITO))
    generic = [f"G{i:04d}" for i in range(n_genes - len(named))]
    genes = np.array(named + generic)
    G = genes.size
    gi = {g: i for i, g in enumerate(genes)}

    # positions: 22 chromosomes of 100 Mb, centromere at 40-45 Mb; shuffle genes over the genome
    chrom = np.array([str(1 + (i % 22)) for i in rng.permutation(G)])
    start = rng.integers(1, 100_000_000, G)
    cyto = pd.DataFrame([(c, 40_000_000, 45_000_000, "cen", "acen") for c in map(str, range(1, 23))] +
                        [(c, 0, 40_000_000, "p1", "gneg") for c in map(str, range(1, 23))],
                        columns=["chrom", "start", "end", "band", "stain"])
    arm = np.where(start < 42_500_000, "p", "q")
    cna_arm_genes = np.flatnonzero((chrom == "7") & (arm == "q") & np.char.startswith(genes, "G"))

    base = rng.lognormal(-1.0, 1.5, G)
    base[[gi[g] for g in RIBO]] *= 30
    base[[gi[g] for g in MITO]] *= 20
    pool = rng.permutation([gi[g] for g in generic if gi[g] not in set(cna_arm_genes)])
    take = iter(np.array_split(pool[: program_size * (5 + 2 * n_tumours)], 5 + 2 * n_tumours))
    programs = {name: next(take) for name in ["cycle", "hypoxia", "EMT", "IFN", "stress_extra"]}
    stress_idx = np.concatenate([[gi[g] for g in IEG_HSP_CORE], programs.pop("stress_extra")[:16]])
    private = {f"T{t + 1:02d}": [next(take), next(take)] for t in range(n_tumours)}

    tme_types = {"T": T_MARKERS, "Myeloid": ["PTPRC"] + MYELOID_MARKERS, "Fibroblast": FIB_MARKERS}
    tme_profile = {}
    for ct, markers in tme_types.items():
        prof = base * rng.lognormal(0, 0.5, G)
        prof[[gi[g] for g in markers]] *= 40
        prof[rng.choice(len(generic), 30, replace=False) + len(named)] *= 10
        tme_profile[ct] = prof

    tumours = [f"T{t + 1:02d}" for t in range(n_tumours)]
    t_ifn = set(rng.choice(tumours, int(round(0.6 * n_tumours)), replace=False))
    t_emt = set(rng.choice(tumours, int(round(0.5 * n_tumours)), replace=False))
    t_cna = set(rng.choice(tumours, int(round(0.7 * n_tumours)), replace=False))
    samples = [(t, "fresh", "R1") for t in tumours]
    samples += [(t, "fixed", "R1") for t in tumours[:n_twin]]
    samples += [(t, "fresh", "R2") for t in tumours[n_twin:n_twin + n_region]]

    blocks, obs_rows = [], []
    for t, protocol, region in samples:
        ti = tumours.index(t)
        s_name = f"{t}_{protocol}_{region}"
        stress_level = 0.0 if protocol == "fixed" else rng.uniform(0.5, 2.0)
        mal_base = base * rng.lognormal(0, 0.3, G)
        rate = np.tile(mal_base, (n_mal, 1))

        def activate(idx, frac, fc, rows=None):
            rows = np.arange(n_mal) if rows is None else rows
            on = rows[rng.random(rows.size) < frac]
            rate[np.ix_(on, idx)] *= fc

        activate(programs["cycle"], 0.25, 4.0)
        activate(programs["hypoxia"], 0.30, 3.0)
        if t in t_emt:
            activate(programs["EMT"], 0.35 if protocol == "fixed" else 0.20, 3.0)
        if t in t_ifn:
            activate(programs["IFN"], 0.30, 3.0)
        for p in private[t]:
            activate(p, 0.30, 3.0)
        clone = np.where(rng.random(n_mal) < 0.5, "c2", "c1") if t in t_cna else np.full(n_mal, "c1")
        if t in t_cna:
            rate[np.ix_(clone == "c2", cna_arm_genes)] *= 3.0
        stress_act = stress_level * rng.gamma(2.0, 0.5, n_mal)
        rate[:, stress_idx] *= (1 + 3 * stress_act)[:, None]

        # doublets with TME cells
        doublet = rng.random(n_mal) < 0.04
        cts = list(tme_types)
        for r in np.flatnonzero(doublet):
            prof = tme_profile[cts[rng.integers(len(cts))]]
            rate[r] = rate[r] / rate[r].sum() + prof / prof.sum()

        # ambient RNA from TME profiles
        amb = sum(tme_profile.values())
        amb = amb / amb.sum()
        rho = rng.beta(1.0, 20.0, n_mal)
        rate = (1 - rho)[:, None] * rate / rate.sum(axis=1, keepdims=True) + rho[:, None] * amb[None, :]

        depth = rng.lognormal(np.log(4000), 0.5, n_mal)
        blocks.append(_counts(rate, depth, 10.0, rng))
        for c in range(n_mal):
            obs_rows.append({"sample": s_name, "tumour": t, "protocol": protocol, "region": region,
                             "cell_type": "Malignant", "malignant": True, "clone": clone[c],
                             "doublet_score": float(rng.beta(8, 2) if doublet[c] else rng.beta(2, 8)),
                             "contamination": float(np.clip(rho[c] + rng.normal(0, 0.01), 0, 1))})

        for ct, prof in tme_profile.items():
            r = np.tile(prof, (n_tme_per_type, 1))
            if t in t_ifn:
                on = rng.random(n_tme_per_type) < 0.3
                r[np.ix_(np.flatnonzero(on), programs["IFN"])] *= 3.0
            act = stress_level * rng.gamma(2.0, 0.5, n_tme_per_type)
            r[:, stress_idx] *= (1 + 3 * act)[:, None]
            depth = rng.lognormal(np.log(3000), 0.5, n_tme_per_type)
            blocks.append(_counts(r, depth, 10.0, rng))
            for _ in range(n_tme_per_type):
                obs_rows.append({"sample": s_name, "tumour": t, "protocol": protocol, "region": region,
                                 "cell_type": ct, "malignant": False, "clone": "none",
                                 "doublet_score": float(rng.beta(2, 8)), "contamination": 0.0})

    obs = pd.DataFrame(obs_rows)
    obs["study"] = obs["tumour"].map({t: f"S{(i % n_studies) + 1}" for i, t in enumerate(tumours)})
    obs["cancer_type"] = obs["tumour"].map({t: ("A" if i % 2 == 0 else "B") for i, t in enumerate(tumours)})
    obs.index = [f"cell{i}" for i in range(len(obs))]
    X = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    var = pd.DataFrame({"chromosome": chrom, "start": start}, index=pd.Index(genes, name="gene"))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()

    truth_sets = {name: list(genes[idx]) for name, idx in programs.items()}
    truth_sets["stress"] = list(genes[stress_idx])
    truth_sets["CNA_7q"] = list(genes[cna_arm_genes])
    truth_sets["contamination"] = list(dict.fromkeys(T_MARKERS + MYELOID_MARKERS + FIB_MARKERS))
    kinds = {"cycle": "biology", "hypoxia": "biology", "EMT": "biology (protocol-selected)",
             "IFN": "TME-shared biology", "stress": "technical (dissociation)",
             "CNA_7q": "genetic (dosage)", "contamination": "technical (ambient/doublet)"}
    truth = {"sets": truth_sets, "kinds": kinds, "cytoband": cyto,
             "private": {t: [list(genes[p]) for p in ps] for t, ps in private.items()},
             "tumours_ifn": sorted(t_ifn), "tumours_emt": sorted(t_emt), "tumours_cna": sorted(t_cna)}
    return adata, truth
