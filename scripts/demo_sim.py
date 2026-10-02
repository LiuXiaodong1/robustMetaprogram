"""Run the P1 audit on the toy cohort and compare each audited MP with the planted truth.

    python scripts/demo_sim.py --out results/demo --n-jobs 4

Every outcome here is true by construction (see mpaudit.sim). The point is to check that the
code runs end to end and that each axis responds to the nuisance it was designed for.
"""
import argparse
from pathlib import Path

import pandas as pd

from mpaudit.compare import jaccard_matrix
from mpaudit.pipeline import run_p1
from mpaudit.sim import simulate_cohort

KEYS = dict(sample="sample", malignant="malignant", tumour="tumour", study="study", cancer_type="cancer_type",
            cell_type="cell_type", protocol="protocol", clone="clone", doublet="doublet_score")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/demo")
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--nrun", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    adata, truth = simulate_cohort(seed=args.seed)
    cfg = {
        "keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour", "n_jobs": args.n_jobs,
        "gavish": {"nrun": args.nrun},
        "decoy": {"n": 3, "strata": ["clone", "doublet_score", "contamination"]},
        "positional": {"cytoband": truth["cytoband"]},
        "protocol": {"reference": "fresh", "test": "fixed"},
    }
    res = run_p1(adata, cfg, args.out)

    flags = res["tables"]["mp_flags"].set_index("mp")
    J = jaccard_matrix(res["mp_sets"], truth["sets"])
    names = list(truth["sets"])
    rows = []
    for i, mp in enumerate(res["mp_sets"]):
        j = int(J[i].argmax())
        rows.append({"mp": mp, "planted": names[j] if J[i, j] >= 0.2 else "(none)",
                     "truth_kind": truth["kinds"][names[j]] if J[i, j] >= 0.2 else "",
                     "jaccard": round(float(J[i, j]), 2),
                     **{c: flags.loc[mp, c] for c in ["recruitment", "partial_r2", "kappa_ratio_median",
                                                      "mean_delta_intensity", "frac_beyond_null_intensity",
                                                      "recovery_all", "best_arm", "q_excluding_known_clusters",
                                                      "flagged"] if c in flags}})
    table = pd.DataFrame(rows)
    table.to_csv(Path(args.out) / "truth_comparison.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.precision", 3):
        print(table)


if __name__ == "__main__":
    main()
