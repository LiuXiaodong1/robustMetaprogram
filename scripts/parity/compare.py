"""Compare R's ``NMF::nmf(method = "snmf/r")`` with the mpaudit port on the same exported matrices.

    python scripts/parity/compare.py --export parity/export --r-out parity/r_out --out parity/report \\
        [--h5ad cohort.h5ad | --python-bases DIR] [--python-replicate-seed 1] [--r-replicate DIR]

The port's bases come from :func:`mpaudit.gavish.compute_programs` on the exported samples with the
exported parameters (``--h5ad``, default: the export's source file; ``--cache-dir`` reuses an
``mpaudit p1`` NMF cache), or are loaded from a basis directory written by an earlier run
(``--python-bases``; every run writes ``<out>/python_bases``). R's bases are read from
``run_r_snmf.R``'s output. Three levels are compared, always on the ranks both sides have:

1. Programs. Per sample and rank, R's k programs are paired one-to-one with the port's k programs
   by the Hungarian assignment that maximises the summed top-``n_top`` overlap (dead, all-zero
   columns are left out). Statistic: matched overlap, 0-50 genes. Also reported: the cosine of the
   paired basis columns, and each basis's "profile objective", the sNMF/R objective R reports
   (``.snmf.objective``, eta = max(A)) with H re-fitted to the exported matrix for that W (a
   convex problem). Lower is a better fit; it tells "different but equally good optimum" from
   "worse optimum", whichever side stopped early.
2. Robust programs (the filter in :func:`mpaudit.gavish.assemble`, run on each side's bases): per
   sample, robust programs are paired by top-gene overlap; a robust program is reproduced when its
   partner shares >= ``intra_min`` [35] genes, the threshold the filter itself uses for "the same
   program at another rank".
3. MPs: :func:`mpaudit.gavish.assemble` on each side, MPs paired with
   :func:`mpaudit.compare.match_sets` (Jaccard >= 0.3, mpaudit's "recovered" threshold); per pair,
   shared genes, member counts and member agreement (R members mapped to their level-1 partners).

Verdict (constants below; also stated in every report):

* PASS: median matched program overlap >= 45/50, every robust program reproduced with the same
  count in every sample, and the same number of MPs, each paired with >= 45/50 shared genes and the
  same number of members. If neither side forms an MP, the MP level is "not tested".
* WITHIN NOISE: not PASS, but R versus port agrees at least as well as every supplied seed
  replicate (port versus port at another seed, R versus R) at every level: median program overlap
  within 2 genes of the replicate's, and robust-program agreement (reproduced / larger count) and
  MP agreement (identical / larger count) at least the replicate's. The solvers differ no more
  than either differs from itself.
* FAIL: otherwise. INVALID: the two sides did not fit the same matrices.

Caveats. sNMF/R is non-convex and both codes keep the best of ``nrun`` random restarts, so even R
against itself at another seed does not give 50/50; run a replicate before reading a FAIL as a port
defect. Program overlaps from one sample's ranks are not independent. With only 3-5 samples the
3CA MP threshold cannot be met (a founder needs > ``min_group_size`` [5] partners, i.e. >= 7
samples), so by default the MP stage lowers ``min_group_size`` so that an MP must recur in all but
one exported sample (at least 3), and switches the v2 post-filters and gene order off; both sides
get the same settings, and ``--mp-as-configured`` keeps the exported ones.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parity_io as pio  # noqa: E402

from mpaudit import gavish  # noqa: E402
from mpaudit.compare import match_sets  # noqa: E402
from mpaudit.data import Keys  # noqa: E402
from mpaudit.gavish import (GavishParams, GavishResult, NMFCache, SamplePrograms, assemble,  # noqa: E402
                            compute_programs, program_id)

PROGRAM_OVERLAP_FRAC = 0.9   # median matched overlap >= 0.9 * n_top (45 of 50)
MP_OVERLAP_FRAC = 0.9        # paired MPs share >= 0.9 * n_top genes
MP_MATCH_JACCARD = 0.3       # pairing threshold for MPs (mpaudit's "recovered")
NOISE_TOLERANCE = 2          # genes: WITHIN NOISE if median overlap >= replicate median - 2
OBJECTIVE_RTOL = 1e-6        # smaller relative differences in profile objective are numerical noise
CHECK_RTOL = 1e-9            # R's matrix sum/max against the export (R may parse the last bit differently)


# --------------------------------------------------------------------------- level 1: programs

def top_lists(W: np.ndarray, universe: np.ndarray, n_top: int) -> list[list[str] | None]:
    """Top genes of each basis column (``None`` for a dead, all-zero column), as in mpaudit."""
    return [None if not np.any(W[:, j] > 0) else [str(g) for g in universe[np.argsort(-W[:, j], kind="stable")[:n_top]]]
            for j in range(W.shape[1])]


def pair_by_overlap(tops_a: list, tops_b: list) -> list[tuple[int, int, int]]:
    """One-to-one pairs ``(i, j, overlap)`` maximising the summed overlap; ``None`` entries are skipped."""
    ia = [i for i, t in enumerate(tops_a) if t is not None]
    ib = [j for j, t in enumerate(tops_b) if t is not None]
    if not ia or not ib:
        return []
    O = np.array([[len(set(tops_a[i]) & set(tops_b[j])) for j in ib] for i in ia], dtype=int)
    ri, ci = linear_sum_assignment(-O)
    return [(ia[r], ib[c], int(O[r, c])) for r, c in zip(ri, ci)]


def _aligned(a: SamplePrograms, b: SamplePrograms) -> tuple[np.ndarray, np.ndarray]:
    """Row indices of the genes ``a`` and ``b`` share (identity when the universes are equal)."""
    if np.array_equal(a.universe, b.universe):
        idx = np.arange(len(a.universe))
        return idx, idx
    ib = b.gene_index()
    shared = [(i, ib[g]) for i, g in enumerate(a.universe) if g in ib]
    return np.array([i for i, _ in shared], dtype=int), np.array([j for _, j in shared], dtype=int)


def compare_programs(a: dict[str, SamplePrograms], b: dict[str, SamplePrograms], ranks, n_top: int) -> pd.DataFrame:
    """Level 1: one row per paired program (sample, rank, program_a, program_b, overlap, cosine)."""
    rows = []
    for s, sa in a.items():
        sb = b[s]
        ia, ib = _aligned(sa, sb)
        for r in ranks:
            Wa, Wb = np.asarray(sa.W[r], dtype=float), np.asarray(sb.W[r], dtype=float)
            for i, j, ov in pair_by_overlap(top_lists(Wa, sa.universe, n_top), top_lists(Wb, sb.universe, n_top)):
                x, y = Wa[ia, i], Wb[ib, j]
                cos = float(x @ y / (np.linalg.norm(x) * np.linalg.norm(y))) if x.any() and y.any() else np.nan
                rows.append({"sample": s, "rank": r, "program_a": program_id(s, r, i), "program_b": program_id(s, r, j),
                             "overlap": ov, "cosine": cos})
    return pd.DataFrame(rows, columns=["sample", "rank", "program_a", "program_b", "overlap", "cosine"])


def profile_objective(A: np.ndarray, W: np.ndarray, beta: float, eta: float, tol: float = 1e-10,
                      max_iter: int = 5000) -> float:
    """R's ``.snmf.objective`` for basis ``W`` with ``H >= 0`` re-fitted to ``A``:
    ``0.5 * (||A - WH||^2 + eta ||W||^2 + beta sum_c (sum_j H[j, c])^2)``.

    For fixed W the problem in H is a convex quadratic program (Gram matrix ``W'W + beta 11'``),
    solved here by exact coordinate descent (HALS) until H changes by less than ``tol`` relative.
    """
    W = np.asarray(W, dtype=np.float64)
    k = W.shape[1]
    WtW = W.T @ W
    G = WtW + beta
    WtA = W.T @ A
    H = np.zeros((k, A.shape[1]))
    for _ in range(max_iter):
        delta = 0.0
        for j in range(k):
            new = np.maximum(H[j] + (WtA[j] - G[j] @ H) / G[j, j], 0.0)
            delta += float(np.sum((new - H[j]) ** 2))
            H[j] = new
        if delta <= tol**2 * max(float(np.sum(H * H)), 1e-300):
            break
    resid = float(np.einsum("ij,ij->", A, A)) - 2.0 * float(np.sum(H * WtA)) + float(np.sum(WtW * (H @ H.T)))
    return 0.5 * (max(resid, 0.0) + eta * float(np.sum(W * W)) + beta * float(np.sum(H.sum(axis=0) ** 2)))


def objectives(export_dir: Path, manifest: dict, r: dict[str, SamplePrograms], py: dict[str, SamplePrograms],
               ranks, params: GavishParams, r_out: Path) -> pd.DataFrame:
    """Per sample and rank: R's reported best objective and both sides' profile objectives."""
    rows = []
    for entry in manifest["samples"]:
        s = entry["sample"]
        A = pio.read_matrix(export_dir, entry)
        eta = params.eta if params.eta is not None and params.eta >= 0 else float(A.max())
        reported = {}
        fits_path = r_out / entry["dir"] / pio.FITS
        if fits_path.exists():
            f = pd.read_csv(fits_path)
            reported = dict(zip(f["rank"].astype(int), pd.to_numeric(f["objective"], errors="coerce")))
        for k in ranks:
            obj_r = profile_objective(A, r[s].W[k], params.beta, eta)
            obj_p = profile_objective(A, py[s].W[k], params.beta, eta)
            rows.append({"sample": s, "rank": k, "r_reported_objective": float(reported.get(k, np.nan)),
                         "r_profile_objective": obj_r, "py_profile_objective": obj_p,
                         "py_over_r": obj_p / obj_r if obj_r > 0 else np.nan,
                         "py_vs_r_pct": 100.0 * (obj_p / obj_r - 1.0) if obj_r > 0 else np.nan})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- levels 2 and 3

def mp_stage_params(params: GavishParams, ranks, n_samples: int, min_group_size: int | None,
                    as_configured: bool) -> tuple[GavishParams, list[str]]:
    """Parameters for :func:`assemble` on both sides, and notes on what was changed."""
    p = params.replace(ranks=tuple(ranks))
    notes = []
    if tuple(ranks) != tuple(params.ranks):
        notes.append(f"ranks restricted to {list(ranks)} (configured {list(params.ranks)}): R output has no other ranks")
    if as_configured:
        return p, notes
    g = min_group_size if min_group_size is not None else min(params.min_group_size, max(1, n_samples - 3))
    names = {f.name for f in dataclasses.fields(GavishParams)}
    changes = {"min_group_size": g, "min_programs_per_mp": 0, "min_studies_per_mp": 1, "v2_gene_order": False}
    changes = {k: v for k, v in changes.items() if k in names and getattr(p, k) != v}
    if changes:
        notes.append("MP stage run with " + ", ".join(f"{k}={v}" for k, v in changes.items())
                     + f" ({n_samples} samples; use --mp-as-configured for the exported values)")
    return p.replace(**changes), notes


def compare_robust(ra: GavishResult, rb: GavishResult, samples: list[str], intra_min: int) -> pd.DataFrame:
    """Level 2: per sample, robust programs paired by top-gene overlap (unpaired ones get an empty partner)."""
    rows = []
    for s in samples:
        pa = [p for p in ra.robust if ra.program_sample[p] == s]
        pb = [p for p in rb.robust if rb.program_sample[p] == s]
        pairs = pair_by_overlap([ra.top[p] for p in pa], [rb.top[p] for p in pb])
        for i, j, ov in pairs:
            rows.append({"sample": s, "program_a": pa[i], "program_b": pb[j], "overlap": ov, "reproduced": ov >= intra_min})
        done_a, done_b = {i for i, _, _ in pairs}, {j for _, j, _ in pairs}
        rows += [{"sample": s, "program_a": p, "program_b": "", "overlap": 0, "reproduced": False}
                 for i, p in enumerate(pa) if i not in done_a]
        rows += [{"sample": s, "program_a": "", "program_b": p, "overlap": 0, "reproduced": False}
                 for j, p in enumerate(pb) if j not in done_b]
    return pd.DataFrame(rows, columns=["sample", "program_a", "program_b", "overlap", "reproduced"])


def compare_mps(ra: GavishResult, rb: GavishResult, pairs: dict[str, str]) -> pd.DataFrame:
    """Level 3: MPs paired by Jaccard; ``pairs`` maps side-a program ids to their level-1 partners."""
    sa, sb = ra.mp_sets(), rb.mp_sets()
    mem_a = {m.name: m.members for m in ra.mps}
    mem_b = {m.name: m.members for m in rb.mps}
    rows = []
    m = match_sets(sa, sb, MP_MATCH_JACCARD)
    for row in m.itertuples(index=False):
        alt = row.alt
        mapped = {pairs.get(p, f"unpaired:{p}") for p in mem_a[row.ref]}
        members_b = set(mem_b.get(alt, []))
        union = mapped | members_b
        rows.append({"mp_a": row.ref, "mp_b": alt, "jaccard": row.jaccard if alt else np.nan,
                     "shared_genes": len(set(sa[row.ref]) & set(sb[alt])) if alt else 0,
                     "n_members_a": len(mem_a[row.ref]), "n_members_b": len(members_b),
                     "member_jaccard": len(mapped & members_b) / len(union) if alt and union else np.nan})
    paired = set(m["alt"]) if len(m) else set()
    rows += [{"mp_a": "", "mp_b": name, "jaccard": np.nan, "shared_genes": 0, "n_members_a": 0,
              "n_members_b": len(mem_b[name]), "member_jaccard": np.nan} for name in sb if name not in paired]
    return pd.DataFrame(rows, columns=["mp_a", "mp_b", "jaccard", "shared_genes", "n_members_a", "n_members_b",
                                       "member_jaccard"])


def compare_pair(a: dict[str, SamplePrograms], b: dict[str, SamplePrograms], ranks, params: GavishParams,
                 mp_params: GavishParams, unit: dict | None, study: dict | None) -> dict:
    """All three levels for two sets of bases over the same samples."""
    prog = compare_programs(a, b, ranks, params.n_top)
    res_a = assemble(a, mp_params, unit_of_sample=unit, study_of_sample=study)
    res_b = assemble(b, mp_params, unit_of_sample=unit, study_of_sample=study)
    robust = compare_robust(res_a, res_b, list(a), mp_params.intra_min)
    mps = compare_mps(res_a, res_b, dict(zip(prog["program_a"], prog["program_b"])))
    n_top = params.n_top
    mp_tested = bool(res_a.mps or res_b.mps)
    summary = {
        "n_programs_paired": int(len(prog)),
        "median_overlap": float(prog["overlap"].median()) if len(prog) else np.nan,
        "p10_overlap": float(prog["overlap"].quantile(0.1)) if len(prog) else np.nan,
        "min_overlap": int(prog["overlap"].min()) if len(prog) else None,
        "median_cosine": float(prog["cosine"].median()) if len(prog) else np.nan,
        "n_robust_a": len(res_a.robust), "n_robust_b": len(res_b.robust),
        "n_robust_reproduced": int(robust["reproduced"].sum()),
        "n_mps_a": len(res_a.mps), "n_mps_b": len(res_b.mps),
        "n_mps_paired": int((mps["mp_a"].ne("") & mps["mp_b"].ne("")).sum()),
        "n_mps_identical": int(((mps["shared_genes"] >= MP_OVERLAP_FRAC * n_top)
                                & (mps["n_members_a"] == mps["n_members_b"]) & mps["mp_a"].ne("")
                                & mps["mp_b"].ne("")).sum()),
        "mp_tested": mp_tested,
    }
    summary["robust_agreement"] = (summary["n_robust_reproduced"] / max(summary["n_robust_a"], summary["n_robust_b"])
                                   if summary["n_robust_a"] or summary["n_robust_b"] else 1.0)
    summary["mp_agreement"] = (summary["n_mps_identical"] / max(summary["n_mps_a"], summary["n_mps_b"])
                               if mp_tested else 1.0)
    summary["program_ok"] = bool(len(prog)) and summary["median_overlap"] >= PROGRAM_OVERLAP_FRAC * n_top
    summary["robust_ok"] = summary["n_robust_a"] == summary["n_robust_b"] == summary["n_robust_reproduced"]
    summary["mp_ok"] = summary["n_mps_a"] == summary["n_mps_b"] == summary["n_mps_identical"]
    return {"programs": prog, "robust": robust, "mps": mps, "summary": summary, "res_a": res_a, "res_b": res_b}


def verdict(main: dict, baselines: dict[str, dict], issues: list[str]) -> tuple[str, str]:
    """Overall verdict and a one-line reason (rules in the module docstring)."""
    s = main["summary"]
    if issues:
        return "INVALID", "the two sides did not fit the same matrices: " + "; ".join(issues[:3])
    if s["program_ok"] and s["robust_ok"] and (s["mp_ok"] or not s["mp_tested"]):
        return "PASS", "programs, robust programs and MPs agree" if s["mp_tested"] else \
            "programs and robust programs agree; MP level not tested (no MPs on either side)"
    failed = [name for name, ok in (("programs", s["program_ok"]), ("robust programs", s["robust_ok"]),
                                    ("MPs", s["mp_ok"] or not s["mp_tested"])) if not ok]
    if baselines and all(within_noise(s, b["summary"]) for b in baselines.values()):
        return "WITHIN NOISE", (f"{', '.join(failed)} differ, but no more than between seed replicates "
                                f"(median program overlap {s['median_overlap']:g})")
    if baselines:
        return "FAIL", f"{', '.join(failed)} differ, more than between seed replicates"
    return "FAIL", f"{', '.join(failed)} differ (no seed replicate supplied)"


def within_noise(s: dict, rep: dict) -> bool:
    """R versus port agrees at least as well as a seed replicate at every level: median program
    overlap within NOISE_TOLERANCE genes, robust-program and MP agreement fractions at least the
    replicate's."""
    return (s["median_overlap"] >= rep["median_overlap"] - NOISE_TOLERANCE
            and s["robust_agreement"] >= rep["robust_agreement"] and s["mp_agreement"] >= rep["mp_agreement"])


# --------------------------------------------------------------------------- inputs

def check_inputs(export_dir: Path, manifest: dict, r_out: Path, sides: dict[str, dict[str, SamplePrograms]]) -> list[str]:
    """Problems showing that a side did not fit the exported matrix (empty when all is well)."""
    issues = []
    for entry in manifest["samples"]:
        s = entry["sample"]
        genes = np.asarray(pio.read_lines(export_dir / entry["dir"] / pio.GENES), dtype=object)
        for side, programs in sides.items():
            if not np.array_equal(np.asarray(programs[s].universe, dtype=object), genes):
                issues.append(f"{side} genes of {s} differ from the export")
        check = r_out / entry["dir"] / pio.INPUT_CHECK
        if not check.exists():
            issues.append(f"no {pio.INPUT_CHECK} for {s} in the R output")
            continue
        c = pd.read_csv(check).iloc[0]
        if int(c["n_genes"]) != entry["n_genes"] or int(c["n_cells"]) != entry["n_cells"]:
            issues.append(f"R read a {c['n_genes']} x {c['n_cells']} matrix for {s}, "
                          f"exported {entry['n_genes']} x {entry['n_cells']}")
        elif not (np.isclose(float(c["sum"]), entry["sum"], rtol=CHECK_RTOL, atol=0)
                  and np.isclose(float(c["max"]), entry["max"], rtol=CHECK_RTOL, atol=0)):
            issues.append(f"R's matrix sum/max for {s} differ from the export")
    return issues


def python_programs(manifest: dict, params: GavishParams, h5ad: str, cache_dir: str | None,
                    n_jobs: int, log=print) -> dict[str, SamplePrograms]:
    """The port's fits for the exported samples (manifest order)."""
    import anndata as ad

    samples = [e["sample"] for e in manifest["samples"]]
    log(f"[compare] fitting the port on {len(samples)} samples, ranks {list(params.ranks)}, seed {params.seed}")
    programs = compute_programs(ad.read_h5ad(h5ad), Keys.from_dict(manifest["keys"]), params,
                                counts_layer=manifest["counts_layer"], samples=samples,
                                cache=NMFCache(cache_dir), n_jobs=n_jobs)
    missing = [s for s in samples if s not in programs]
    if missing:
        raise RuntimeError(f"The port did not fit {missing}: cell selection differs from the export")
    return {s: programs[s] for s in samples}


# --------------------------------------------------------------------------- report

def guide(n_top: int) -> str:
    t = int(np.ceil(PROGRAM_OVERLAP_FRAC * n_top))
    return f"""## How to read this report

- **PASS**: median matched program overlap >= {t}/{n_top}, every robust program reproduced (partner
  shares >= intra_min genes) with the same count per sample, and the same number of MPs, each paired
  with >= {int(np.ceil(MP_OVERLAP_FRAC * n_top))}/{n_top} shared genes and the same number of members. Absolute numbers from the
  port can be quoted as 3CA-equivalent for this setting.
- **WITHIN NOISE**: not PASS, but R versus port agrees at least as well as every seed replicate
  (port versus port, or R versus R, at another seed) at every level: median program overlap within
  {NOISE_TOLERANCE} genes, robust_agreement and mp_agreement at least the replicate's. The solvers differ no
  more than a rerun of either; quote absolute numbers with that caveat and the replicate's numbers.
- **FAIL**: look at the profile objectives. If the port's basis fits worse at most sample-ranks
  (py_vs_r_pct > 0), however small the difference (early stopping costs little objective but moves
  top genes), the port stops short of R: lower `tol` (e.g. 1e-7) or raise `max_iter`, re-export,
  rerun. If the signs are mixed, the solvers reach different, equally good optima (NMF is not
  identifiable); then robust-program and MP agreement decide whether MP-level results survive the
  solver change, and absolute counts must not be quoted as 3CA-equivalent.
  Without a seed replicate a FAIL is not interpretable: rerun with `--python-replicate-seed`.
- **INVALID**: the two sides did not fit the same matrices (see Input checks); nothing else in the
  report means anything.

Statistics. Programs: per sample and rank, R's and the port's programs paired one-to-one
(Hungarian, maximising summed top-{n_top} overlap); overlaps of one sample's ranks are not
independent. Cosine: of paired basis columns. Profile objective: R's `.snmf.objective` with H
re-fitted for each W; lower is better. Robust programs and MPs: `mpaudit.gavish.assemble` on each
side's bases with the same parameters; MPs paired by `mpaudit.compare.match_sets` (Jaccard >=
{MP_MATCH_JACCARD}); member_jaccard compares member programs after mapping R's to their paired port programs.
"""


def objective_summary(obj: pd.DataFrame | None) -> dict:
    """Which side's basis fits the matrix worse (higher profile objective), counted over sample-ranks.

    The direction matters more than the size: stopping early raises the objective very little
    (on the toy cohort, tol=1e-5 instead of 1e-8 cost 0.002-0.03%) yet moves top-50 lists.
    """
    if obj is None or obj.empty:
        return {}
    ratio = obj["py_over_r"].dropna()
    return {"n": int(ratio.size), "median_py_vs_r_pct": float(100.0 * (ratio.median() - 1.0)),
            "n_port_worse": int((ratio > 1 + OBJECTIVE_RTOL).sum()), "n_r_worse": int((ratio < 1 - OBJECTIVE_RTOL).sum())}


def objective_text(o: dict) -> str:
    if not o:
        return "Profile objectives not computed."
    return (f"Profile objective, port versus R: median {o['median_py_vs_r_pct']:+.4f}%; the port's basis fits worse "
            f"at {o['n_port_worse']} of {o['n']} sample-ranks, R's at {o['n_r_worse']} (ties within a relative "
            f"{OBJECTIVE_RTOL:g} not counted).")


SUMMARY_COLS = ["median_overlap", "p10_overlap", "min_overlap", "median_cosine", "n_robust_a", "n_robust_b",
                "n_robust_reproduced", "robust_agreement", "n_mps_a", "n_mps_b", "n_mps_identical", "mp_agreement"]


def rank_table(prog: pd.DataFrame, obj: pd.DataFrame | None) -> pd.DataFrame:
    tbl = prog.groupby(["sample", "rank"], sort=False).agg(
        n_paired=("overlap", "size"), median_overlap=("overlap", "median"), min_overlap=("overlap", "min"),
        median_cosine=("cosine", "median")).reset_index()
    if obj is not None and len(obj):
        tbl = tbl.merge(obj, on=["sample", "rank"], how="left")
    return tbl


def write_report(out: Path, main: dict, baselines: dict, ranks_tbl: pd.DataFrame, info: dict) -> None:
    s = main["summary"]
    lines = ["# R versus port parity report", "",
             f"**Verdict: {info['verdict']}** - {info['reason']}", "",
             "## Inputs", "",
             f"- export: `{info['export']}` (mpaudit {info['export_version']}, code {info['export_code']}); "
             f"{len(info['samples'])} samples: {', '.join(info['samples'])}",
             f"- R output: `{info['r_out']}`"
             + "".join(f"; {k} {v}" for k, v in info["r_versions"].items()),
             f"- port: {info['python_source']} (mpaudit {info['mpaudit_version']}, code {info['code_version']})",
             f"- ranks {info['ranks']}, nrun {info['nrun']}, n_top {info['n_top']}, R seed {info['r_seed']}, "
             f"port seed {info['py_seed']}"]
    lines += [f"- note: {n}" for n in info["notes"]]
    lines += ["", "## Input checks", ""]
    lines += [f"- {i}" for i in info["issues"]] or ["- all samples: R read the exported matrix (dimensions, sum, max) "
                                                    "and both sides' gene rows equal the export"]
    lines += ["", "## Summary", "",
              pio.md_table(pd.DataFrame([{"comparison": k, **{c: v["summary"][c] for c in SUMMARY_COLS}}
                                         for k, v in {"R vs port": main, **baselines}.items()])),
              "## Programs by rank (all samples)", "",
              pio.md_table(main["programs"].groupby("rank").agg(
                  n_paired=("overlap", "size"), median_overlap=("overlap", "median"),
                  p10_overlap=("overlap", lambda x: x.quantile(0.1)), min_overlap=("overlap", "min"),
                  median_cosine=("cosine", "median")).reset_index()),
              "## Programs by sample and rank", "", objective_text(info["objective"]), "", pio.md_table(ranks_tbl),
              "## Robust programs", "",
              f"R {s['n_robust_a']}, port {s['n_robust_b']}, reproduced {s['n_robust_reproduced']}.", "",
              pio.md_table(main["robust"].groupby("sample", sort=False).agg(
                  n_r=("program_a", lambda x: int((x != "").sum())), n_port=("program_b", lambda x: int((x != "").sum())),
                  reproduced=("reproduced", "sum"), min_overlap=("overlap", "min")).reset_index()),
              "## MPs", "",
              f"R {s['n_mps_a']}, port {s['n_mps_b']}, paired {s['n_mps_paired']}, identical {s['n_mps_identical']}.", "",
              pio.md_table(main["mps"]), guide(info["n_top"])]
    (out / "report.md").write_text("\n".join(lines))


def jsonable(v):
    """``v`` with numpy scalars turned into Python ones and NaN into None (strict JSON)."""
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [jsonable(x) for x in v]
    if isinstance(v, (np.integer, np.bool_)):
        return v.item()
    if isinstance(v, (float, np.floating)):
        return None if np.isnan(v) else float(v)
    return v


# --------------------------------------------------------------------------- main

def _h5ad(args, manifest: dict) -> str:
    h5ad = args.h5ad or manifest["sources"].get("h5ad")
    if not h5ad:
        raise SystemExit("Give --h5ad (the export does not name its source file)")
    return h5ad


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", required=True, help="directory written by export_matrices.py")
    ap.add_argument("--r-out", required=True, help="directory written by run_r_snmf.R")
    ap.add_argument("--out", default="parity/report")
    ap.add_argument("--h5ad", help="cohort to fit the port on (default: the export's source file)")
    ap.add_argument("--python-bases", help="load port bases written by an earlier run (<out>/python_bases)")
    ap.add_argument("--cache-dir", help="mpaudit NMF cache directory to reuse (e.g. the one from 'mpaudit p1')")
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--python-replicate-seed", type=int, help="also fit the port at this seed (noise floor)")
    ap.add_argument("--r-replicate", help="a second run_r_snmf.R output at another --seed (noise floor)")
    ap.add_argument("--mp-min-group-size", type=int, help="min_group_size for the MP stage (default: see above)")
    ap.add_argument("--mp-as-configured", action="store_true", help="MP stage with the exported parameters unchanged")
    ap.add_argument("--no-objective", action="store_true", help="skip the profile objectives (faster)")
    args = ap.parse_args(argv)

    export_dir, r_out, out = Path(args.export), Path(args.r_out), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    manifest = pio.read_manifest(export_dir)
    params = pio.params_from_manifest(manifest)
    samples = [e["sample"] for e in manifest["samples"]]
    ranks = [k for k in params.ranks if k in set(pio.available_ranks(r_out, manifest))]
    if not ranks:
        raise SystemExit(f"No rank of {list(params.ranks)} has R bases for every sample under {r_out}")
    fit_params = params.replace(ranks=tuple(ranks))
    if args.python_replicate_seed == params.seed:
        raise SystemExit(f"--python-replicate-seed must differ from the exported seed ({params.seed})")

    r_bases = pio.read_bases(r_out, manifest, ranks)
    if args.python_bases:
        py_bases = pio.read_bases(args.python_bases, manifest, ranks)
        python_source = f"bases loaded from `{args.python_bases}`"
    else:
        py_bases = python_programs(manifest, fit_params, _h5ad(args, manifest), args.cache_dir, args.n_jobs)
        pio.write_bases(py_bases, out / "python_bases", manifest)
        python_source = f"fitted on `{_h5ad(args, manifest)}` (bases written to `{out / 'python_bases'}`)"

    issues = check_inputs(export_dir, manifest, r_out, {"R": r_bases, "port": py_bases})
    mp_params, notes = mp_stage_params(params, ranks, len(samples), args.mp_min_group_size, args.mp_as_configured)
    unit = {e["sample"]: e["unit"] for e in manifest["samples"]} if manifest.get("recurrence_key") else None
    studies = {e["sample"]: e["study"] for e in manifest["samples"]}
    study = studies if all(v is not None for v in studies.values()) else None
    main_cmp = compare_pair(r_bases, py_bases, ranks, params, mp_params, unit, study)

    baselines = {}
    if args.python_replicate_seed is not None:
        seed = args.python_replicate_seed
        rep = python_programs(manifest, fit_params.replace(seed=seed), _h5ad(args, manifest), args.cache_dir,
                              args.n_jobs)
        pio.write_bases(rep, out / f"python_bases_seed{seed}", manifest)
        baselines[f"port vs port (seed {seed})"] = compare_pair(py_bases, rep, ranks, params, mp_params, unit, study)
    if args.r_replicate:
        baselines["R vs R (replicate)"] = compare_pair(r_bases, pio.read_bases(args.r_replicate, manifest, ranks),
                                                       ranks, params, mp_params, unit, study)

    obj = None if args.no_objective or issues else objectives(export_dir, manifest, r_bases, py_bases, ranks, params, r_out)
    ranks_tbl = rank_table(main_cmp["programs"], obj)
    v, reason = verdict(main_cmp, baselines, issues)

    r_versions = {}
    if (r_out / "versions.csv").exists():
        vv = pd.read_csv(r_out / "versions.csv", keep_default_na=False).set_index("key")["value"]
        r_versions = {k: vv[k] for k in ("R", "NMF", "seed") if k in vv.index}
    code_now = gavish.code_version() if hasattr(gavish, "code_version") else None
    if manifest.get("code_version") not in (None, code_now):
        notes.append(f"export made with mpaudit code {manifest['code_version']}, port now {code_now}: "
                     "re-export if prepare_matrix changed (the input checks would show it)")
    info = {
        "verdict": v, "reason": reason, "export": str(export_dir), "r_out": str(r_out), "samples": samples,
        "export_version": manifest.get("mpaudit_version"), "export_code": manifest.get("code_version"),
        "mpaudit_version": metadata.version("mpaudit"), "code_version": code_now, "r_versions": r_versions,
        "python_source": python_source, "ranks": ranks, "nrun": params.nrun, "n_top": params.n_top,
        "r_seed": manifest.get("r_seed"), "py_seed": params.seed, "notes": notes, "issues": issues,
        "objective": objective_summary(obj),
    }
    main_cmp["programs"].to_csv(out / "programs.csv", index=False)
    ranks_tbl.to_csv(out / "ranks.csv", index=False)
    main_cmp["robust"].to_csv(out / "robust.csv", index=False)
    main_cmp["mps"].to_csv(out / "mps.csv", index=False)
    for name, b in baselines.items():
        b["programs"].to_csv(out / f"programs_{name.split(' (')[0].replace(' ', '_')}.csv", index=False)
    write_report(out, main_cmp, baselines, ranks_tbl, info)
    summary = {**info, "r_vs_port": main_cmp["summary"], "baselines": {k: b["summary"] for k, b in baselines.items()},
               "thresholds": {"program_overlap_frac": PROGRAM_OVERLAP_FRAC, "mp_overlap_frac": MP_OVERLAP_FRAC,
                              "mp_match_jaccard": MP_MATCH_JACCARD, "noise_tolerance": NOISE_TOLERANCE,
                              "objective_rtol": OBJECTIVE_RTOL, "check_rtol": CHECK_RTOL},
               "mp_stage_params": dataclasses.asdict(mp_params)}
    (out / "summary.json").write_text(json.dumps(jsonable(summary), indent=2, allow_nan=False))
    print(f"[compare] {v}: {reason}; see {out / 'report.md'}")
    return summary


if __name__ == "__main__":
    main()
