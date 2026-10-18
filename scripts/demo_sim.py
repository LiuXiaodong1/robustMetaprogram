"""Run the P1 audit on the toy cohort and compare each audited MP with the planted truth.

    python scripts/demo_sim.py --out results/demo --n-jobs 4
    python scripts/demo_sim.py --clean --out results/demo_clean --n-jobs 4
    python scripts/demo_sim.py --true-doublet --out results/demo_true_doublet --n-jobs 4

The default cohort plants biology, TME-shared biology, a CNA, and technical nuisance (dissociation
stress, ambient RNA and doublets, an abundance/%ribo axis, a long-gene axis), together with the
confounders that can fake nuisance (malignant->TME spillover, cycling cells with higher doublet
scores, protocol capture bias). The decoys condition on ``clone_est``, ``doublet_est`` and
``contamination_est``, estimated from the counts (never the planted clone, doublet score or
contamination), on ``pct_ribo``, computed from the counts, and on ``intronic_fraction``. The toy
draws ``intronic_fraction`` as the planted nuclear fraction plus N(0, 0.03): a near-true readout
standing in for a read-level measurement (STARsolo, velocyto, DropletQC), so long-gene recruitment
by the technical arm is true by construction. ``--clean`` removes every planted nuisance program and
the CNA but keeps the confounders: it is the null cohort on which the frozen rule should not say GO.

On the default cohort the estimated doublet-cycle coupling is weak (seed 0: a ``doublet_est`` arm
keeps 0.00 of the cycle program's excess kappa, the registered technical arm 0.01; see the
``mpaudit.decoy`` docstring), so the default run hardly exercises the ``kappa_retained`` >= 0.5
guard of criterion (a). On the ``--clean`` cohort, without TME ambient RNA, ``doublet_est`` tracks
cycling more (a ``doublet_est`` arm keeps 0.26). ``--true-doublet`` puts the planted
``doublet_score`` (cycling cells score 0.3 higher) into the technical arm instead of
``doublet_est``, to run a real program tightly coupled to a covariate through NMF recruitment and
that guard. Arm contents are registered as filled in but not compared by the code, so all three runs
carry the same rule label; the ``--true-doublet`` arm is not one a registration would fix (the
planted score is not observable in real data).

The toy has no curated artefacts and no curated labels, and its MPs are re-derived: under the
registered rule ``mps: rederive`` needs ``typer.labels`` and ``typer.labelled_gmt``, so every run
here is the unregistered variant ``P1-GO v1.3 (unregistered variant: typer.labels=None)``
(``unregistered_variant: true``), with no other deviation.

Every outcome here is true by construction (see mpaudit.sim). The point is to check that the code
runs, that each axis responds to the nuisance it was designed for, and how often real biology is
flagged (the false-flag rate). ``truth_summary.json`` holds the verdict and rule label, the counts
of planted biology and technical MPs that were counted, biology MPs typed dosage
(``biology_dosage``), biology MPs that the technical decoys recruit in >= ``min_hits`` of >=
``min_decoys`` decoys whatever their ``kappa_retained`` (``biology_tech_recruited``; per MP whether
the ``kappa_retained`` guard passed and whether it is typed dosage, the two things that keep such an
MP from counting), how many counted MPs meet each criterion, and the code commit.
"""
import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd

from mpaudit.compare import jaccard_matrix
from mpaudit.pipeline import run_p1
from mpaudit.report import validate_rule
from mpaudit.sim import simulate_cohort

KEYS = dict(sample="sample", malignant="malignant", tumour="tumour", study="study", cancer_type="cancer_type",
            cell_type="cell_type", protocol="protocol", clone="clone_est", doublet="doublet_est",
            intronic="intronic_fraction")
FLAG_COLS = ["tech_decoy_hits", "tech_decoy_n", "tech_kappa_retained", "dosage_decoy_hits", "spill_adjusted",
             "partial_r2", "null_q99", "slope", "partial_r2_rho", "null_q99_rho", "rho_sensitivity",
             "partial_r2_unadjusted", "n_paired_tumours", "frac_loss_beyond_null_kappa", "positional_arm",
             "positional_arm_excluding_known", "q_excluding_known_clusters", "recovery_all", "expected_tme_shared",
             "dosage", "crit_a_tech_decoy", "crit_b_tme", "crit_c_protocol", "counted"]


def code_commit() -> dict:
    """Short commit of the repository holding this script, and whether tracked files differ from it."""
    here = Path(__file__).resolve().parent
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=here, capture_output=True, text=True,
                                check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=here,
                               capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return {"code_commit": "unknown", "code_uncommitted_changes": None}
    return {"code_commit": commit or "unknown", "code_uncommitted_changes": bool(dirty)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/demo")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--nrun", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--clean", action="store_true", help="null cohort: no planted nuisance programs and no CNA")
    ap.add_argument("--no-instability", action="store_true")
    ap.add_argument("--cache-dir", help="keep NMF fits on disk so an interrupted run resumes")
    ap.add_argument("--true-doublet", action="store_true",
                    help="technical arm on the planted doublet_score (coupled to cycling) instead of doublet_est")
    args = ap.parse_args()

    flags = dict(dissociation_stress=False, tme_contamination=False, cna=False, technical_axes=False) if args.clean else {}
    adata, truth = simulate_cohort(seed=args.seed, **flags)
    cfg = {
        "keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour", "n_jobs": args.n_jobs,
        "cache_dir": args.cache_dir,
        "gavish": {"nrun": args.nrun, "min_genes_per_cell": 0},  # toy cohort has only 1,500 genes
        "instability": {"enabled": not (args.no_instability or args.clean)},
        "decoy": {"n": 5, "arms": {"technical": ["doublet_score" if args.true_doublet else "doublet_est",
                                                 "contamination_est", "pct_ribo", "intronic_fraction"],
                                   "dosage": ["clone_est"]}},
        "positional": {"cytoband": truth["cytoband"]},
        "protocol": {"reference": "fresh", "test": "fixed"},
        # pre-declared before the run: IFN is expected to co-vary with TME cells while being real. The MPs
        # are re-derived (named MP1, MP2, ... by the run), so the exclusion is a gene set, not a name.
        "decision": {"expected_tme_shared_sets": {"IFN": truth["sets"]["IFN"]}},
        # the toy has no curated labels, which mps: rederive needs under the registered rule (the only deviation)
        "unregistered_variant": True,
    }
    res = run_p1(adata, cfg, args.out)

    fl = res["tables"]["mp_flags"].set_index("mp")
    J = jaccard_matrix(res["mp_sets"], truth["sets"])
    names = list(truth["sets"])
    rows = []
    for i, mp in enumerate(res["mp_sets"]):
        j = int(J[i].argmax()) if J.shape[1] else 0
        hit = J.shape[1] and J[i, j] >= 0.2
        kind = truth["kinds"][names[j]] if hit else "private / unmatched"
        rows.append({"mp": mp, "planted": names[j] if hit else "(none)", "truth_kind": kind,
                     "jaccard": round(float(J[i, j]), 2) if J.shape[1] else 0.0,
                     **{c: fl.loc[mp, c] for c in FLAG_COLS if c in fl}})
    table = pd.DataFrame(rows)
    table.to_csv(Path(args.out) / "truth_comparison.csv", index=False)
    biology = table["truth_kind"].str.contains("biology")
    technical = table["truth_kind"].str.startswith("technical")
    rule = validate_rule(cfg["decision"])
    counted = table["counted"].astype(bool)
    recruited = biology & (table["tech_decoy_hits"] >= rule["min_hits"]) & (table["tech_decoy_n"] >= rule["min_decoys"])
    summary = {
        "verdict": res["decision"]["verdict"], "n_counted": res["decision"]["n_counted"],
        "biology_mps": int(biology.sum()), "biology_counted": int((biology & counted).sum()),
        "technical_mps": int(technical.sum()), "technical_counted": int((technical & counted).sum()),
        "dosage_mps": res["decision"]["dosage_mps"],
        # biology typed dosage (e.g. a biology MP that absorbed CNA genes): labelled, never counted
        "biology_dosage": int((biology & table["dosage"].astype(bool)).sum()),
        # biology MPs the technical decoys recruit whatever their kappa_retained: one that is not counted was
        # stopped by the kappa_retained >= min_kappa_retained guard, by a dosage label, or by both
        "biology_tech_recruited": int(recruited.sum()),
        "biology_tech_recruited_mps": [
            {"mp": r.mp, "planted": r.planted, "hits": int(r.tech_decoy_hits),
             "kappa_retained": None if pd.isna(r.tech_kappa_retained) else round(float(r.tech_kappa_retained), 3),
             "kappa_guard_passed": bool(r.tech_kappa_retained >= rule["min_kappa_retained"]),
             "dosage": bool(r.dosage), "counted": bool(r.counted)} for r in table[recruited].itertuples()],
        "counted_by_criterion": {c: int((counted & table[col].astype(bool)).sum()) for c, col in
                                 (("a", "crit_a_tech_decoy"), ("b", "crit_b_tme"), ("c", "crit_c_protocol"))},
        "technical_arm": cfg["decoy"]["arms"]["technical"],
        "rule_version": res["decision"]["rule_version"], "registered_rule": res["decision"]["registered_rule"],
        "deviations": res["decision"]["deviations"],
        **code_commit(),
    }
    (Path(args.out) / "truth_summary.json").write_text(json.dumps(summary, indent=2))
    with pd.option_context("display.width", 300, "display.max_columns", 40, "display.precision", 3):
        print(table)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
