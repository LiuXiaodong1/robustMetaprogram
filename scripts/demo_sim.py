"""Run the P1 audit on the toy cohort and compare each audited MP with the planted truth.

    python scripts/demo_sim.py --out results/demo --n-jobs 4
    python scripts/demo_sim.py --clean --out results/demo_clean --n-jobs 4

The default cohort plants biology, TME-shared biology, a CNA, and technical nuisance (dissociation
stress, ambient RNA and doublets, an abundance/%ribo axis, a long-gene axis), together with the
confounders that can fake nuisance (malignant->TME spillover, cycling cells with higher doublet
scores, protocol capture bias). Decoys condition only on covariates estimated from the counts.
``--clean`` removes every planted nuisance program and the CNA but keeps the confounders: it is the
null cohort on which the frozen rule should not say GO.

Every outcome here is true by construction (see mpaudit.sim). The point is to check that the code
runs, that each axis responds to the nuisance it was designed for, and how often real biology is
flagged (the false-flag rate).
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from mpaudit.compare import jaccard_matrix
from mpaudit.pipeline import run_p1
from mpaudit.sim import simulate_cohort

KEYS = dict(sample="sample", malignant="malignant", tumour="tumour", study="study", cancer_type="cancer_type",
            cell_type="cell_type", protocol="protocol", clone="clone_est", doublet="doublet_est",
            intronic="intronic_fraction")
FLAG_COLS = ["tech_decoy_hits", "tech_kappa_retained", "dosage_decoy_hits", "spill_adjusted", "partial_r2", "null_q99",
             "slope", "partial_r2_rho", "null_q99_rho", "rho_sensitivity", "partial_r2_unadjusted", "n_paired_tumours",
             "frac_loss_beyond_null_kappa", "positional_arm", "positional_arm_excluding_known",
             "q_excluding_known_clusters", "recovery_all", "expected_tme_shared", "dosage",
             "crit_a_tech_decoy", "crit_b_tme", "crit_c_protocol", "counted"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/demo")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--nrun", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--clean", action="store_true", help="null cohort: no planted nuisance programs and no CNA")
    ap.add_argument("--no-instability", action="store_true")
    args = ap.parse_args()

    flags = dict(dissociation_stress=False, tme_contamination=False, cna=False, technical_axes=False) if args.clean else {}
    adata, truth = simulate_cohort(seed=args.seed, **flags)
    cfg = {
        "keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour", "n_jobs": args.n_jobs,
        "gavish": {"nrun": args.nrun, "min_genes_per_cell": 0},  # toy cohort has only 1,500 genes
        "instability": {"enabled": not (args.no_instability or args.clean)},
        "decoy": {"n": 5, "arms": {"technical": ["doublet_est", "contamination_est", "pct_ribo", "intronic_fraction"],
                                   "dosage": ["clone_est"]}},
        "positional": {"cytoband": truth["cytoband"]},
        "protocol": {"reference": "fresh", "test": "fixed"},
        # pre-declared before the run: IFN is expected to co-vary with TME cells while being real. The MPs
        # are re-derived (named MP1, MP2, ... by the run), so the exclusion is a gene set, not a name.
        "decision": {"expected_tme_shared_sets": {"IFN": truth["sets"]["IFN"]}},
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
    summary = {
        "verdict": res["decision"]["verdict"], "n_counted": res["decision"]["n_counted"],
        "biology_mps": int(biology.sum()), "biology_counted": int((biology & table["counted"]).sum()),
        "technical_mps": int(technical.sum()), "technical_counted": int((technical & table["counted"]).sum()),
        "dosage_mps": res["decision"]["dosage_mps"],
    }
    (Path(args.out) / "truth_summary.json").write_text(json.dumps(summary, indent=2))
    with pd.option_context("display.width", 300, "display.max_columns", 40, "display.precision", 3):
        print(table)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
