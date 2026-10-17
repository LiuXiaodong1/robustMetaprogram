"""Decoy arms, kappa_retained and the unpenalised-intercept GLM (mpaudit.decoy)."""
import numpy as np
import pandas as pd
import pytest

from mpaudit import GavishParams, Keys, NMFCache
from mpaudit.gavish import DECOY_TAG
from mpaudit.data import get_counts, malignant_mask
from mpaudit.decoy import (PER_DECOY_COLUMNS, PER_MP_COLUMNS, _covariates, _glm_mean, _strata, decoy_recruitment,
                           decoy_tag, kappa_retained, make_decoy)
from mpaudit.sim import simulate_cohort

KEYS = Keys(sample="sample", malignant="malignant", tumour="tumour", study="study")
SETS = ("cycle", "contamination", "CNA_7q", "stress")
ARMS = {"technical": ("doublet_est", "contamination_est"), "dosage": ("clone_est",)}
PARAMS = GavishParams(ranks=(3, 4), nrun=1, max_iter=200, min_genes_per_cell=0)


@pytest.fixture(scope="module")
def cohort():
    return simulate_cohort(n_tumours=6, n_twin=2, n_region=1, n_mal=120, n_tme_per_type=25, n_genes=700, seed=3)


@pytest.fixture(scope="module")
def sets(cohort):
    return {k: cohort[1]["sets"][k] for k in SETS}


def _retained(cohort, strata, sets, seed=7) -> pd.Series:
    adata, _ = cohort
    dec = make_decoy(adata, KEYS, counts_layer="counts", strata=strata, seed=seed)
    return kappa_retained(adata, dec, KEYS, sets, counts_layer="counts").set_index("mp")["kappa_retained"]


def _malignant_rows(adata):
    mal = malignant_mask(adata, KEYS)
    samples = adata.obs["sample"].astype(str).to_numpy()
    return {s: np.flatnonzero(mal & (samples == s)) for s in pd.unique(samples[mal])}


def test_glm_intercepts_are_not_shrunk(cohort):
    """Unpenalised stratum intercepts: fitted means sum to the observed counts per stratum, low-count genes too."""
    adata, _ = cohort
    X = get_counts(adata, "counts")
    for s, rows in list(_malignant_rows(adata).items())[:3]:
        Y = X[rows].toarray()
        Y = Y[:, Y.sum(axis=0) > 0]
        obs, totals = adata.obs.iloc[rows], Y.sum(axis=1)
        st = _strata(obs, totals, ("clone_est",), 1, 3)
        mu = _glm_mean(Y, totals / np.median(totals), st, _covariates(obs, ("doublet_est", "contamination_est")))
        low = Y.mean(axis=0) < 0.1
        assert low.sum() >= 10
        for h in np.unique(st):
            fit, seen = mu[st == h].sum(axis=0), Y[st == h].sum(axis=0)
            np.testing.assert_allclose(fit[seen > 0], seen[seen > 0], rtol=1e-4)
        assert np.isfinite(mu).all() and mu.max() < 10 * Y.max()  # no IRLS overshoot


def test_decoy_keeps_low_count_gene_totals(cohort):
    """Per sample, decoy totals of low-count genes stay within 10% of the original (no upward pull)."""
    adata, _ = cohort
    dec = make_decoy(adata, KEYS, counts_layer="counts", strata=("clone_est", "doublet_est", "contamination_est"),
                     seed=11)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    ratios = []
    for rows in _malignant_rows(adata).values():
        mean = np.asarray(X[rows].mean(axis=0)).ravel()
        low = np.flatnonzero((mean > 0) & (mean < 0.25))
        ratios.append(D[rows][:, low].sum() / X[rows][:, low].sum())
        assert abs(D[rows].sum() / X[rows].sum() - 1) < 0.05
    ratios = np.array(ratios)
    assert np.all(np.abs(ratios - 1) < 0.1), ratios
    assert abs(ratios.mean() - 1) < 0.03


def test_contamination_arm_keeps_contamination_coexpression(cohort, sets):
    kr = _retained(cohort, ("contamination_est",), sets)
    assert kr["contamination"] >= 0.5


def test_doublet_arms_keep_little_of_cycle(cohort, sets):
    """False-recruitment check: cycling cells get higher doublet scores (cell size). A doublet_est-only
    decoy keeps almost none of the cycle program's co-expression; one on the planted doublet_score,
    shifted up by 0.3 in cycling cells, keeps a visible share that stays below the rule's 0.5."""
    est = _retained(cohort, ("doublet_est",), sets)
    assert est["cycle"] < 0.15 and est["stress"] < 0.15
    planted = _retained(cohort, ("doublet_score",), sets)
    assert 0.2 < planted["cycle"] < 0.5


def test_kappa_retained_identity_and_nan(cohort, sets):
    adata, _ = cohort
    same = kappa_retained(adata, adata, KEYS, sets, counts_layer="counts").set_index("mp")
    np.testing.assert_allclose(same.loc[["cycle", "CNA_7q"], "kappa_retained"], 1.0)
    assert (same.loc[["cycle", "CNA_7q"], "n_samples"] > 0).all()
    none = kappa_retained(adata, adata, KEYS, sets, counts_layer="counts", min_kappa=2.0)
    assert none["kappa_retained"].isna().all() and (none["n_samples"] == 0).all()


def test_arm_output_format_and_typing(cohort, sets):
    adata, _ = cohort
    cache = NMFCache()
    per_mp, per_decoy, decoy_mps = decoy_recruitment(adata, KEYS, PARAMS, sets, arms=ARMS, n_decoys=2,
                                                     counts_layer="counts", cache=cache, recurrence_key="tumour",
                                                     log=lambda *a: None)
    assert list(per_mp.columns) == PER_MP_COLUMNS and list(per_decoy.columns) == PER_DECOY_COLUMNS
    assert len(per_mp) == len(ARMS) * len(sets) and set(per_mp["arm"]) == set(ARMS)
    assert (per_mp["n_decoys"] == 2).all() and per_mp["hits"].between(0, 2).all()
    np.testing.assert_allclose(per_mp["recruitment"], per_mp["hits"] / per_mp["n_decoys"])
    assert (per_mp["max_best_jaccard"] >= per_mp["mean_best_jaccard"]).all()
    assert len(per_decoy) == len(ARMS) * 2 and set(per_decoy["decoy"]) == set(decoy_mps)
    # NMF fit diagnostics of each decoy cohort: samples x ranks x nrun fits, of which some capped or failed
    n_samples = cohort[0].obs.loc[malignant_mask(cohort[0], KEYS), "sample"].nunique()
    assert (per_decoy["n_fits"] == n_samples * len(PARAMS.ranks) * PARAMS.nrun).all()
    assert per_decoy["n_capped_fits"].between(0, per_decoy["n_fits"]).all()
    assert per_decoy["n_failed_fits"].between(0, per_decoy["n_fits"]).all()
    tech = per_decoy.loc[per_decoy["arm"] == "technical", "decoy"]
    assert tech.str.contains("contamination_est+doublet_est", regex=False).all()  # sorted columns in the cache tag
    assert tech.str.contains("-residual-").all() and tech.str.contains("-d1-").all()
    assert per_decoy["decoy"].str.startswith(DECOY_TAG).all()  # so the cache key hashes the decoy generator
    assert set(tech.str.extract(r"-s(\d+)-[0-9a-f]{12}$")[0]) == {"1000", "1001"}

    kr = per_mp.set_index(["arm", "mp"])["kappa_retained"]
    assert kr["dosage", "CNA_7q"] >= 0.5 > kr["technical", "CNA_7q"]  # clone arm types CNA as dosage
    assert kr["technical", "contamination"] >= 0.5 > kr["dosage", "contamination"]

    # `strata` alone is one arm called "all"; same columns, method, depth bins and seed -> same cache tag
    n_cached = len(cache.mem)
    per_mp_all, per_decoy_all, _ = decoy_recruitment(adata, KEYS, PARAMS, sets, strata=("clone_est",), n_decoys=1,
                                                     counts_layer="counts", cache=cache, recurrence_key="tumour",
                                                     log=lambda *a: None)
    assert set(per_mp_all["arm"]) == {"all"} and len(per_mp_all) == len(sets)
    assert per_decoy_all["decoy"].iloc[0] == per_decoy.loc[per_decoy["arm"] == "dosage", "decoy"].iloc[0]
    assert len(cache.mem) == n_cached  # served from the cache

    with pytest.raises(ValueError):
        decoy_recruitment(adata, KEYS, PARAMS, sets, arms=ARMS, strata=("clone_est",))


def test_decoy_tag_follows_every_argument_and_the_covariate_values(cohort):
    """The cache tag is built from the make_decoy call itself, defaults included, and from the values it
    conditions on: with method="permute" changed covariate values leave the counts' totals unchanged."""
    adata, _ = cohort
    kw = dict(counts_layer="counts", strata=("clone_est",), method="permute", seed=1)
    tag = decoy_tag(adata, KEYS, **kw)
    assert tag.startswith(f"{DECOY_TAG}permute-clone_est-d1-s1-")
    assert decoy_tag(adata, KEYS, **kw, n_cont_bins=3) == tag         # a default passed explicitly
    assert decoy_tag(adata, KEYS, **kw, min_stratum_cells=25) != tag  # any other argument
    assert decoy_tag(adata, KEYS, **{**kw, "seed": 2}) != tag
    moved = adata.copy()
    assert decoy_tag(moved, KEYS, **kw) == tag
    vals = moved.obs["clone_est"].to_numpy().copy()
    rows = np.flatnonzero(malignant_mask(moved, KEYS))
    vals[rows] = vals[np.roll(rows, 1)]  # the same column name and labels, other values in malignant cells
    moved.obs["clone_est"] = pd.Series(vals, index=moved.obs.index, dtype=moved.obs["clone_est"].dtype)
    assert decoy_tag(moved, KEYS, **kw) != tag
