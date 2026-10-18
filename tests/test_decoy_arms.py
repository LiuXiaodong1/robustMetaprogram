"""Decoy arms, kappa_retained (P1-GO v1.3: the excess over matched random sets) and the unpenalised-intercept
GLM (mpaudit.decoy)."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit import GavishParams, Keys, NMFCache
from mpaudit.coexpr import kappa, matched_random_sets, set_kappas
from mpaudit.gavish import DECOY_TAG
from mpaudit.data import add_qc_obs, get_counts, malignant_mask
from mpaudit.decoy import (EXCESS_ENTRY_SD, MIN_CELLS, MIN_DETECT, MIN_FLOOR_SETS, MIN_GENES, N_NULL_FLOOR,
                           PER_DECOY_COLUMNS, PER_MP_COLUMNS, _covariates, _glm_mean, _reference_kappa, _set_rng,
                           _strata, decoy_recruitment, decoy_tag, kappa_retained, make_decoy)
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
    # the samples entering kappa_retained, decided on the real counts: the same in every arm, as kappa_retained counts
    n = per_mp.set_index(["arm", "mp"])["kappa_retained_n_samples"]
    alone = kappa_retained(adata, adata, KEYS, sets, counts_layer="counts", seed=1000).set_index("mp")["n_samples"]
    for arm in ARMS:
        assert n[arm].to_dict() == alone.to_dict()
    assert (n > 0).all()

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
    """The cache tag is built from the make_decoy call itself, defaults included, and from what it reads in
    malignant cells: the barcodes, and the dtype and values of the columns it conditions on."""
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
    # the dtype: make_decoy bins a numeric column but stratifies a categorical one by value (same hashed values)
    cont = {**kw, "strata": ("doublet_est",)}
    as_category = adata.copy()
    as_category.obs["doublet_est"] = as_category.obs["doublet_est"].astype("category")
    assert decoy_tag(as_category, KEYS, **cont) != decoy_tag(adata, KEYS, **cont)
    # only malignant cells enter: dropping a TME cell keeps the tag, renaming a malignant cell does not
    mal = malignant_mask(adata, KEYS)
    keep = np.ones(adata.n_obs, bool)
    keep[np.flatnonzero(~mal)[5]] = False
    assert decoy_tag(adata[keep].copy(), KEYS, **kw) == tag
    renamed = adata.copy()
    renamed.obs_names = [f"x{n}" if i == rows[0] else n for i, n in enumerate(renamed.obs_names)]
    assert decoy_tag(renamed, KEYS, **kw) != tag


# --------------------------------------------------------------------------- kappa_retained, P1-GO v1.3

def test_kappa_retained_constants():
    """Registered and fixed constants of kappa_retained (a change is a new rule version)."""
    assert N_NULL_FLOOR == 50 and MIN_FLOOR_SETS == 20 and EXCESS_ENTRY_SD == 3.0
    assert MIN_CELLS == 20 and MIN_DETECT == 0.05 and MIN_GENES == 5
    with pytest.raises(ValueError, match="n_null_floor must be >= 20"):
        _reference_kappa(None, KEYS, {}, None, 20, 0.05, 5, n_null_floor=19)


def test_set_kappas_equal_kappa_set_by_set(cohort, sets):
    """Every kappa of a sample from one set_kappas call: the values of kappa() set by set (here bit for bit),
    edge cases included (a one-gene set, a gene without counts, a gene listed twice)."""
    adata, _ = cohort
    X = get_counts(adata, "counts")
    rows = next(iter(_malignant_rows(adata).values()))
    Xs = X[rows]
    s = np.asarray(Xs.sum(axis=1)).ravel()
    detect = np.asarray((Xs > 0).mean(axis=0)).ravel()
    mean_expr = np.asarray(Xs.sum(axis=0)).ravel() / s.sum()
    pos = {g: i for i, g in enumerate(adata.var_names)}
    g = np.array(sorted(pos[x] for x in sets["cycle"] if detect[pos[x]] >= 0.05))
    rs = matched_random_sets(g, mean_expr, detect, 50, np.random.default_rng(0), candidates=detect >= 0.05)
    zero = int(np.flatnonzero(np.asarray(Xs.sum(axis=0)).ravel() == 0)[0])
    all_sets = [g, *rs, g[:1], np.r_[g[:5], zero], np.r_[g[:5], g[:2]]]
    one_by_one = np.array([kappa(Xs[:, r].toarray(), s) for r in all_sets])
    together = set_kappas(Xs.tocsc(), s, all_sets, chunk=97)  # a chunk size that splits the union unevenly
    assert np.isnan(together[-3]) and np.isnan(one_by_one[-3])
    np.testing.assert_allclose(together, one_by_one, rtol=1e-12, atol=1e-14)
    assert np.array_equal(set_kappas(Xs.toarray(), s, all_sets), together, equal_nan=True)  # dense input
    # the path of real-scale samples (cells x union genes above max_dense: each set read column by column, and the
    # gene moments one gene at a time when the memory budget allows no more), on CSR, float32 and integer input
    for counts in (Xs.tocsc(), Xs.tocsr(), Xs.astype(np.float32).tocsc(), Xs.toarray().astype(np.int64)):
        np.testing.assert_allclose(set_kappas(counts, s, all_sets, max_dense=0), one_by_one, rtol=1e-12, atol=1e-14)


def _dilution_sample(rng, n: int, a2_mp: float, n_mp=30, n_bg=600, n_ab=40) -> np.ndarray:
    """Counts (cells x genes) of one sample: ``n_mp`` program genes whose log latent expression shares a factor
    of variance ``a2_mp``, ``n_bg`` background genes, and ``n_ab`` abundant genes that carry about half of the
    counts and vary from cell to cell. Normalising by total counts dilutes every other gene by one per-cell
    factor: the background co-expression that every gene set of the sample shares (kappa about 0.3)."""
    G = n_mp + n_bg + n_ab
    mu = np.r_[2.0 * rng.uniform(0.8, 1.25, n_mp + n_bg), np.full(n_ab, 30.0)]
    depth, c, f = rng.lognormal(0, 0.3, n), rng.normal(0, 0.5, n), rng.standard_normal(n)
    load = np.zeros(G)
    load[:n_mp] = np.sqrt(a2_mp)
    sd = np.r_[np.sqrt(0.1 - load[:n_mp + n_bg] ** 2), np.full(n_ab, 0.1)]
    z = f[:, None] * load[None, :] + rng.standard_normal((n, G)) * sd[None, :]
    z[:, n_mp + n_bg:] += c[:, None]
    return rng.poisson(depth[:, None] * mu[None, :] * np.exp(z))


def _dilution_adata(blocks: list[np.ndarray], n_mp=30, n_bg=600, n_ab=40) -> ad.AnnData:
    obs = pd.concat([pd.DataFrame({"sample": f"S{k:02d}", "malignant": True}, index=[f"S{k}_c{i}" for i in range(len(b))])
                     for k, b in enumerate(blocks)])
    var = pd.DataFrame(index=[f"MP{i}" for i in range(n_mp)] + [f"BG{i}" for i in range(n_bg)]
                       + [f"AB{i}" for i in range(n_ab)])
    return ad.AnnData(sp.csr_matrix(np.vstack(blocks).astype(np.float32)), obs=obs, var=var)


DILUTION_KEYS = Keys(sample="sample", malignant="malignant")
DILUTION_SET = {"mp": [f"MP{i}" for i in range(30)]}


def test_kappa_retained_is_the_excess_share():
    """A decoy keeping half of the program's own covariance keeps about 0.5; one keeping only the background
    every gene set shares keeps about 0, where the v1.2 ratio kappa_decoy / kappa_real was about 0.4."""
    def sample(a2, seed):  # one sample of 1500 cells, the same draws but for the program factor's share
        return _dilution_sample(np.random.default_rng(seed), 1500, a2)

    real = _dilution_adata([sample(0.08, 0)])
    half, background = (_dilution_adata([sample(a2, 1)]) for a2 in (0.04, 0.0))
    info = _reference_kappa(real, DILUTION_KEYS, DILUTION_SET, None, 20, 0.05, 5)[0][2]["mp"]
    assert info["enter"] and info["bg_real"] > 0.2 and len(info["sets"]) == N_NULL_FLOOR
    k = {name: kappa_retained(real, dec, DILUTION_KEYS, DILUTION_SET).iloc[0]
         for name, dec in (("half", half), ("background", background))}
    assert abs(k["half"]["kappa_retained"] - 0.5) < 0.15 and abs(k["background"]["kappa_retained"]) < 0.15
    assert (k["half"]["n_samples"], k["background"]["n_samples"]) == (1, 1)
    X = get_counts(background)
    ratio_v12 = kappa(X[:, :30].toarray(), np.asarray(X.sum(axis=1)).ravel()) / info["k_real"]
    assert ratio_v12 > 0.3


def test_kappa_retained_entry_needs_three_sds_above_the_background():
    """A sample enters when kappa_real - bg_real > max(min_kappa, 3 SD of its random sets' kappa). Samples above
    the v1.2 floor (the q95 of the random sets) but below 3 SD do not enter; n_samples counts the entering
    samples (identity decoy: each enters with ratio 1)."""
    rng = np.random.default_rng(0)
    a2s = [0.0, 0.003, 0.006, 0.009, 0.012, 0.015, 0.02, 0.03]  # program strength per sample
    adata = _dilution_adata([_dilution_sample(rng, 600, a2) for a2 in a2s])
    refs = _reference_kappa(adata, DILUTION_KEYS, DILUTION_SET, None, 20, 0.05, 5)
    gap, entering = [], []
    for s, rows, per_set in refs:
        i = per_set["mp"]
        assert i["kappa_sets"].size >= MIN_FLOOR_SETS
        assert i["bg_real"] == pytest.approx(np.median(i["kappa_sets"]))
        assert i["sd_real"] == pytest.approx(np.std(i["kappa_sets"], ddof=1))
        excess = i["k_real"] - i["bg_real"]
        assert i["enter"] == (excess > max(0.02, EXCESS_ENTRY_SD * i["sd_real"]))
        if i["k_real"] >= max(0.02, np.quantile(i["kappa_sets"], 0.95)) and not i["enter"]:
            gap.append(s)
        entering += [s] if i["enter"] else []
    assert len(gap) >= 2 and 1 <= len(entering) < len(a2s)
    res = kappa_retained(adata, adata, DILUTION_KEYS, DILUTION_SET).iloc[0]
    assert res["n_samples"] == len(entering) and res["kappa_retained"] == pytest.approx(1.0)


def test_random_sets_have_a_stream_per_sample_and_set(cohort, sets):
    """One set's kappa_retained does not change when other sets are added or the order changes."""
    adata, _ = cohort
    dec = make_decoy(adata, KEYS, counts_layer="counts", strata=("doublet_est",), seed=7)
    a = kappa_retained(adata, dec, KEYS, {"cycle": sets["cycle"]}, counts_layer="counts").set_index("mp")
    b = kappa_retained(adata, dec, KEYS, {k: sets[k] for k in ("stress", "CNA_7q", "cycle")},
                       counts_layer="counts").set_index("mp")
    assert a.loc["cycle", "n_samples"] > 0
    assert a.loc["cycle", "kappa_retained"] == b.loc["cycle", "kappa_retained"]
    assert a.loc["cycle", "n_samples"] == b.loc["cycle", "n_samples"]
    # the random sets themselves: the (seed, sample, set) stream
    r1 = _reference_kappa(adata, KEYS, {"cycle": sets["cycle"]}, "counts", 20, 0.05, 5, seed=3)
    r2 = _reference_kappa(adata, KEYS, {"stress": sets["stress"], "cycle": sets["cycle"]}, "counts", 20, 0.05, 5,
                          seed=3)
    s, rows, per = r1[0]
    assert r2[0][0] == s and len(per["cycle"]["sets"]) == len(r2[0][2]["cycle"]["sets"]) == N_NULL_FLOOR
    assert all(np.array_equal(x, y) for x, y in zip(per["cycle"]["sets"], r2[0][2]["cycle"]["sets"]))
    X = get_counts(adata, "counts")[rows]
    mean_expr = np.asarray(X.sum(axis=0)).ravel() / X.sum()
    detect = np.asarray((X > 0).mean(axis=0)).ravel()
    expected = matched_random_sets(per["cycle"]["genes"], mean_expr, detect, N_NULL_FLOOR, _set_rng(3, s, "cycle"),
                                   candidates=detect >= 0.05)
    assert len(per["cycle"]["sets"]) == N_NULL_FLOOR  # every random set has a finite kappa here
    assert all(np.array_equal(x, y) for x, y in zip(per["cycle"]["sets"], expected))


TECH_ARM = ("contamination_est", "doublet_est", "intronic_fraction", "pct_ribo")  # the registered (demo) arm


@pytest.fixture(scope="module")
def default_toy():
    adata, truth = simulate_cohort(seed=0)
    add_qc_obs(adata, counts_layer="counts")
    return adata, truth


def test_kappa_retained_on_the_default_toy(default_toy):
    """Default toy, seed 0, one decoy (seed 1000, kappa_retained seed 1000 as in decoy_recruitment): the technical
    arm keeps most of the technical sets' excess co-expression and little of the biology's; the clone arm keeps
    most of the CNA set's."""
    adata, truth = default_toy
    names = ("abundance", "long_genes", "contamination", "cycle", "hypoxia", "EMT", "IFN", "CNA_7q")
    gs = {k: truth["sets"][k] for k in names}
    tech = make_decoy(adata, KEYS, counts_layer="counts", strata=TECH_ARM, seed=1000)
    kr = kappa_retained(adata, tech, KEYS, gs, counts_layer="counts", seed=1000).set_index("mp")["kappa_retained"]
    assert (kr[["abundance", "long_genes", "contamination"]] >= 0.5).all(), kr
    assert (kr[["cycle", "hypoxia", "EMT", "IFN"]] <= 0.2).all(), kr
    clone = make_decoy(adata, KEYS, counts_layer="counts", strata=("clone_est",), seed=1000)
    cna = kappa_retained(adata, clone, KEYS, {"CNA_7q": gs["CNA_7q"]}, counts_layer="counts", seed=1000)
    assert cna["kappa_retained"].iloc[0] >= 0.5


def test_minority_biology_program_stays_below_the_guard():
    """16-tumour toy, seed 2, EMT kept in its first two carrier tumours plus every tumour without it, technical
    arm (the minority case of scripts/calibration/evaluate_a.py): v1.2's floor let non-carrier samples in and
    EMT reached kappa_retained 0.66 at this seed in the harness (0.71 with this test's single decoy, 7
    non-carrier samples entering; 0.14 at seed 1, where v1.2 passes too); now at most one non-carrier sample
    enters and EMT stays below 0.5."""
    adata, truth = simulate_cohort(n_tumours=16, seed=2)
    add_qc_obs(adata, counts_layer="counts")
    tumour = adata.obs["tumour"].astype(str).to_numpy()
    carriers = sorted(map(str, truth["tumours_emt"]))
    keep = np.isin(tumour, carriers[:2] + [t for t in np.unique(tumour) if t not in carriers])
    a = adata[keep].copy()
    d = make_decoy(a, KEYS, counts_layer="counts", strata=TECH_ARM, seed=1000)  # (the harness decoys all 16)
    gs = {"EMT": truth["sets"]["EMT"]}
    res = kappa_retained(a, d, KEYS, gs, counts_layer="counts", seed=1000).iloc[0]
    assert res["kappa_retained"] < 0.5
    tum_of = dict(zip(a.obs["sample"].astype(str), a.obs["tumour"].astype(str)))
    entering = [s for s, _, per in _reference_kappa(a, KEYS, gs, "counts", 20, 0.05, 5, seed=1000)
                if per.get("EMT", {}).get("enter")]
    assert sum(tum_of[s] not in carriers for s in entering) <= 1
    assert res["n_samples"] == len(entering) > 0
