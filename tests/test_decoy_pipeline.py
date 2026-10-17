import numpy as np
import pytest

from mpaudit import Keys
from mpaudit.coexpr import kappa
from mpaudit.data import get_counts, malignant_mask
from mpaudit.decoy import make_decoy
from mpaudit.pipeline import run_p1
from mpaudit.sim import simulate_cohort

KEYS = dict(sample="sample", malignant="malignant", tumour="tumour", study="study", cancer_type="cancer_type",
            cell_type="cell_type", protocol="protocol", clone="clone", doublet="doublet_score")


@pytest.fixture(scope="module")
def small():
    return simulate_cohort(n_tumours=6, n_twin=2, n_region=1, n_mal=90, n_tme_per_type=25, n_genes=700, seed=3)


def _sample_rows(adata, keys, s):
    mal = malignant_mask(adata, keys)
    return np.flatnonzero(mal & (adata.obs["sample"].astype(str).to_numpy() == s))


def test_permute_decoy_keeps_stratum_marginals_and_tme(small):
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", strata=("clone",), n_depth_bins=1, method="permute", seed=1)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    tme = ~malignant_mask(adata, keys)
    assert abs(X[tme] - D[tme]).sum() == 0
    s = adata.obs["sample"].iloc[0]
    rows = _sample_rows(adata, keys, s)
    for c in adata.obs["clone"].iloc[rows].unique():
        r = rows[adata.obs["clone"].iloc[rows].to_numpy() == c]
        assert np.allclose(np.asarray(X[r].sum(axis=0)), np.asarray(D[r].sum(axis=0)))


def test_residual_decoy_destroys_program_coexpression(small):
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", method="residual", seed=2)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    pos = {g: i for i, g in enumerate(adata.var_names)}
    idx = [pos[g] for g in truth["sets"]["cycle"]]
    rows = _sample_rows(adata, keys, adata.obs["sample"].iloc[0])

    def mean_corr(M):  # depth-robust latent co-expression; raw log counts would keep depth effects
        Y = M[rows]
        return kappa(Y[:, idx].toarray(), np.asarray(Y.sum(axis=1)).ravel())

    assert mean_corr(X) > 0.15
    assert abs(mean_corr(D)) < 0.05
    tot_x, tot_d = X[rows].sum(), D[rows].sum()
    assert abs(tot_d - tot_x) / tot_x < 0.15


def test_run_p1_end_to_end(small, tmp_path):
    adata, truth = small
    cfg = {
        "keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour",
        "gavish": {"ranks": [3, 4], "nrun": 1, "max_iter": 200, "min_genes_per_cell": 0},
        "instability": {"seeds": [1], "n_boot": 1},
        "decoy": {"n": 1, "arms": {"technical": ["doublet_est", "contamination_est"], "dosage": ["clone_est"]}},
        "tme": {"n_boot": 20, "n_null": 40},
        "positional": {"cytoband": truth["cytoband"]},
        "protocol": {"reference": "fresh", "test": "fixed", "n_null_sets": 20},
        "decision": {"expected_tme_shared_sets": {"IFN": truth["sets"]["IFN"]}},
    }
    out = run_p1(adata, cfg, tmp_path, log=lambda *a: None)
    d = out["decision"]
    assert d["rule_version"] == "P1-GO v1.0"
    assert d["verdict"] in {"GO", "PIVOT", "INCOMPLETE"}
    assert not d["powered"]["technical_decoys"]  # one decoy per arm is below the rule's minimum of 5
    for f in ["summary.md", "mp_flags.csv", "run_config.json", "typer_features.csv", "decoy_recruitment.csv",
              "protocol_shift.csv", "tme_covariation.csv"]:
        assert (tmp_path / f).exists()
    assert set(out["tables"]["decoy_recruitment"]["arm"]) == {"technical", "dosage"}


def test_residual_decoy_keeps_covariate_driven_coexpression(small):
    # estimated covariates, as an analyst would have them; the planted doublet_score is coupled to
    # cycling by design (see test_decoy_arms.py for that false-recruitment check)
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", strata=("clone_est", "doublet_est", "contamination_est"),
                     seed=4)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    pos = {g: i for i, g in enumerate(adata.var_names)}

    def k(M, name):
        idx = [pos[g] for g in truth["sets"][name]]
        vals = []
        for s in adata.obs["sample"].unique()[:4]:
            Y = M[_sample_rows(adata, keys, s)]
            vals.append(kappa(Y[:, idx].toarray(), np.asarray(Y.sum(axis=1)).ravel()))
        return np.nanmean(vals)

    assert k(D, "contamination") > 0.6 * k(X, "contamination")  # continuous covariate kept
    assert k(D, "cycle") < 0.2 * k(X, "cycle")                  # real program destroyed
