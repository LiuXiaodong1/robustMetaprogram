"""The toy cohort's planted effects respond as designed (mpaudit.sim)."""
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.special import comb
from scipy.stats import pearsonr, spearmanr

import mpaudit.sim as sim
from mpaudit.coexpr import kappa
from mpaudit.positional import assign_arms, centromeres
from mpaudit.sim import FIB_MARKERS, MYELOID_MARKERS, RIBO, T_MARKERS, _estimate_covariates, simulate_cohort

SMALL = dict(n_tumours=6, n_twin=2, n_region=1, n_mal=100, n_tme_per_type=30, n_genes=800, seed=0)
FLAGS = ("malignant_ambient", "variable_prevalence", "doublet_cycle_coupling", "technical_axes",
         "protocol_capture_bias", "estimated_covariates")
# flags removing the planted nuisance (off together with technical_axes: the null cohort) -> the set each removes
NUISANCE = {"dissociation_stress": "stress", "tme_contamination": "contamination", "cna": "CNA_7q"}


@lru_cache(maxsize=None)
def _cohort(*off: str):
    """The small cohort with the named flags switched off (shared across tests; do not modify)."""
    return simulate_cohort(**SMALL, **{f: False for f in off})


def _dense(adata) -> np.ndarray:
    return adata.layers["counts"].toarray()


def _cpm(X: np.ndarray) -> np.ndarray:
    return X / X.sum(axis=1, keepdims=True) * 1e6


def _idx(adata, genes) -> list[int]:
    pos = {g: i for i, g in enumerate(adata.var_names)}
    return [pos[g] for g in genes]


def _median_within_sample_spearman(x, y, sample) -> float:
    return float(np.median([spearmanr(x[sample == s], y[sample == s])[0] for s in np.unique(sample)]))


def _ari(a, b) -> float:
    ct = pd.crosstab(np.asarray(a), np.asarray(b)).to_numpy()
    pairs, rows, cols = comb(ct, 2).sum(), comb(ct.sum(axis=1), 2).sum(), comb(ct.sum(axis=0), 2).sum()
    expected = rows * cols / comb(ct.sum(), 2)
    return float((pairs - expected) / ((rows + cols) / 2 - expected))


def test_new_columns_and_truth():
    adata, truth = _cohort()
    obs, var, sets = adata.obs, adata.var, truth["sets"]
    for col in ["sample", "tumour", "protocol", "region", "cell_type", "malignant", "clone", "doublet_score",
                "contamination", "study", "cancer_type", "clone_est", "doublet_est", "contamination_est",
                "intronic_fraction"]:
        assert col in obs
    assert (var["end"] == var["start"] + var["gene_length"]).all()
    assert var["gene_length"].between(1_000, 2_500_000).all()
    assert 15_000 < var["gene_length"].median() < 40_000
    assert var.loc[sets["long_genes"], "gene_length"].between(150_000, 800_000).all()
    assert var.loc[sets["abundance"], "gene_length"].between(2_000, 8_000).all()
    assert set(RIBO) <= set(sets["abundance"]) and len(sets["long_genes"]) >= 30
    planted = set().union(*(sets[k] for k in sets if k not in ("abundance", "long_genes")),
                          *(p for ps in truth["private"].values() for p in ps))
    axis_genes = (set(sets["abundance"]) - set(RIBO)) | set(sets["long_genes"])
    assert not axis_genes & planted
    assert not set(truth["malignant_markers"]) & (axis_genes | planted)
    # TME-enriched genes: recorded, generic, never axis or lineage genes, but they may be program genes
    enriched = truth["tme_enriched"]
    assert set(enriched) == {"T", "Myeloid", "Fibroblast"} and all(len(g) == 30 for g in enriched.values())
    tme_enriched = set().union(*enriched.values())
    assert all(g.startswith("G") for g in tme_enriched)
    assert not tme_enriched & (axis_genes | set(truth["malignant_markers"]))
    assert tme_enriched & planted  # seed 0 puts planted program genes among them (module docstring)
    cpm, ct = _cpm(_dense(adata)), obs["cell_type"].to_numpy()
    for t, g in enriched.items():
        idx = _idx(adata, g)
        fold = (cpm[ct == t][:, idx].mean(axis=0) + 1) / (cpm[ct == "Malignant"][:, idx].mean(axis=0) + 1)
        assert np.median(fold) > 3
    assert truth["kinds"]["abundance"] == "technical (abundance / %ribo)"
    assert truth["kinds"]["long_genes"] == "technical (gene length / nuclear fraction)"
    assert set(truth["kinds"]) == set(sets) and truth["expected_tme_shared"] == ["IFN"]

    mal = obs["malignant"].to_numpy(bool)
    assert set(obs.loc[~mal, "clone_est"]) == {"none"} and set(obs.loc[mal, "clone_est"]) == {"e1", "e2"}
    assert obs["doublet_est"].between(0, 1).all() and obs["intronic_fraction"].between(0, 1).all()
    assert obs["contamination_est"].notna().all() and obs["contamination_est"].between(0, 1).all()
    assert spearmanr(obs["intronic_fraction"], truth["cells"]["nuclear"])[0] > 0.9


def test_hypoxia_prevalence_varies_between_samples():
    adata, truth = _cohort()
    mal = adata.obs["malignant"].to_numpy(bool)
    frac = pd.Series(truth["hypoxia_fraction_by_sample"])
    assert frac.between(0.1, 0.5).all() and frac.std() > 0.05
    observed = truth["cells"]["hypoxia"][mal].groupby(adata.obs["sample"][mal]).mean()
    assert pearsonr(observed[frac.index], frac)[0] > 0.8
    assert set(_cohort("variable_prevalence")[1]["hypoxia_fraction_by_sample"].values()) == {0.30}


def test_flags_give_paired_ablations():
    on, t_on = _cohort()
    off, t_off = _cohort(*FLAGS)
    mal = on.obs["malignant"].to_numpy(bool)
    for col in ["sample", "cell_type", "clone"]:
        assert (on.obs[col] == off.obs[col]).all()
    assert np.array_equal(on.obs["contamination"][mal], off.obs["contamination"][mal])
    pd.testing.assert_frame_equal(t_on["cells"][["cycle", "EMT", "IFN", "doublet"]],
                                  t_off["cells"][["cycle", "EMT", "IFN", "doublet"]])
    assert all(t_on["sets"][k] == t_off["sets"][k] for k in t_off["sets"]) and t_on["private"] == t_off["private"]
    # everything new is absent with the flags off
    assert (off.obs.loc[~mal, "contamination"] == 0).all() and "abundance" not in t_off["sets"]
    assert (t_off["capture_efficiency"] == 1).all() and "clone_est" not in off.obs
    # an effect acting on some samples only leaves the others' counts unchanged
    no_capture, _ = _cohort("protocol_capture_bias")
    fresh = (on.obs["protocol"] == "fresh").to_numpy()
    assert (on.layers["counts"][fresh] != no_capture.layers["counts"][fresh]).nnz == 0


def _share(X: np.ndarray, rows: np.ndarray, idx) -> float:
    """Mean per-cell share of counts on genes ``idx`` among ``rows``."""
    return float((X[rows][:, idx].sum(axis=1) / X[rows].sum(axis=1)).mean())


def test_nuisance_flags_give_paired_ablations():
    """Each nuisance flag alone, and all nine flags together, keep the all-on cohort's cells and states."""
    on, t_on = _cohort()
    obs = on.obs
    # cells whose doublet score is the plain singlet draw in every cohort (TME cells included)
    plain = ~t_on["cells"]["doublet"].to_numpy() & ~t_on["cells"]["cycle"].to_numpy()
    assert t_on["cells"]["doublet"].sum() >= 10
    for off_flags in [(f,) for f in NUISANCE] + [FLAGS + tuple(NUISANCE)]:
        off, t_off = _cohort(*off_flags)
        states = ["cycle", "EMT", "IFN", "cytoplasmic", "nuclear"]
        if "variable_prevalence" not in off_flags:  # that flag sets the hypoxic fraction itself
            states.append("hypoxia")
        pd.testing.assert_frame_equal(t_on["cells"][states], t_off["cells"][states])
        for col in ["sample", "cell_type"] + (["clone"] if "cna" not in off_flags else []):
            assert (obs[col] == off.obs[col]).all(), (off_flags, col)
        assert np.array_equal(obs["doublet_score"][plain], off.obs["doublet_score"][plain]), off_flags
        assert t_on["private"] == t_off["private"]
        assert all(t_on["sets"][k] == t_off["sets"][k] for k in t_off["sets"])
        assert not {NUISANCE[f] for f in off_flags if f in NUISANCE} & set(t_off["sets"])


def test_tme_contamination_flag_removes_doublets_and_ambient_only():
    on, t_on = _cohort()
    off, t_off = _cohort("tme_contamination")
    mal = on.obs["malignant"].to_numpy(bool)
    X_on, X_off = _dense(on), _dense(off)
    # TME cells untouched: bit-identical counts and covariates
    assert (X_on[~mal] != X_off[~mal]).sum() == 0
    for col in ["doublet_score", "contamination", "intronic_fraction"]:
        assert np.array_equal(on.obs[col][~mal], off.obs[col][~mal]), col
    # no doublets, and the removed ones score like singlets ...
    assert t_on["cells"]["doublet"].any() and not t_off["cells"]["doublet"].any()
    removed = (t_on["cells"]["doublet"] & ~t_on["cells"]["cycle"]).to_numpy()
    assert on.obs["doublet_score"][removed].mean() > 0.65 and off.obs["doublet_score"][removed].mean() < 0.35
    # ... and rho = 0: the recorded fraction is the (paired) noise alone, never above the planted one
    c_on, c_off = on.obs["contamination"][mal].to_numpy(), off.obs["contamination"][mal].to_numpy()
    assert (c_off <= c_on).all() and c_off.max() < 0.05 and c_on.mean() > 0.03
    singlet = mal & ~t_on["cells"]["doublet"].to_numpy()
    markers = _idx(on, dict.fromkeys(T_MARKERS + MYELOID_MARKERS + FIB_MARKERS))
    assert _share(X_on, singlet, markers) > 1.2 * _share(X_off, singlet, markers)


def test_dissociation_stress_flag_removes_stress_scaling():
    on, t_on = _cohort()
    off, _ = _cohort("dissociation_stress")
    mal = on.obs["malignant"].to_numpy(bool)
    fixed = (on.obs["protocol"] == "fixed").to_numpy()
    X_on, X_off = _dense(on), _dense(off)
    assert (X_on[fixed] != X_off[fixed]).sum() == 0  # no stress under the fixed protocol either way
    stress = _idx(on, t_on["sets"]["stress"])
    for grp in (mal, ~mal):
        fresh_on, fresh_off = _share(X_on, grp & ~fixed, stress), _share(X_off, grp & ~fixed, stress)
        assert fresh_on > 2 * fresh_off
        # unscaled: fresh cells now carry about the fixed cells' stress share (capture bias aside)
        assert 0.67 < fresh_off / _share(X_off, grp & fixed, stress) < 1.5
        assert fresh_on / _share(X_on, grp & fixed, stress) > 2


def test_cna_flag_removes_clones_and_dosage():
    on, t_on = _cohort()
    off, t_off = _cohort("cna")
    mal = on.obs["malignant"].to_numpy(bool)
    X_on, X_off = _dense(on), _dense(off)
    assert set(off.obs.loc[mal, "clone"]) == {"c1"} and t_off["tumours_cna"] == []
    cna_t = on.obs["tumour"].isin(t_on["tumours_cna"]).to_numpy()
    assert cna_t.any() and (X_on[~cna_t] != X_off[~cna_t]).sum() == 0  # other tumours untouched
    # the planted subclone (labels from the all-on cohort) no longer differs on the CNA arm
    c2 = mal & (on.obs["clone"] == "c2").to_numpy()
    c1 = mal & cna_t & ~c2
    arm = _idx(on, t_on["sets"]["CNA_7q"])
    assert _share(X_on, c2, arm) / _share(X_on, c1, arm) > 2
    assert 0.8 < _share(X_off, c2, arm) / _share(X_off, c1, arm) < 1.25


def test_malignant_ambient_leaks_malignant_genes_into_tme():
    on, truth = _cohort()
    off, _ = _cohort("malignant_ambient")
    mal = on.obs["malignant"].to_numpy(bool)
    c_on, c_off = _cpm(_dense(on)), _cpm(_dense(off))
    m_off, t_off = c_off[mal].mean(axis=0), c_off[~mal].mean(axis=0)
    # malignant-enriched genes, leaving out the lineage markers whose TME level the flag also sets
    enriched = (m_off > 2 * t_off) & ~np.isin(on.var_names, truth["malignant_markers"])
    assert enriched.sum() >= 20
    assert c_on[~mal][:, enriched].mean() / c_off[~mal][:, enriched].mean() > 1.05
    tme_cont = on.obs["contamination"][~mal]
    assert 0.02 < tme_cont.mean() < 0.1 and (off.obs["contamination"][~mal] == 0).all()
    # the lineage markers in TME cells come from the soup: their share follows the true fraction
    lineage = _idx(on, truth["malignant_markers"])
    X = _dense(on)[~mal]
    assert spearmanr(X[:, lineage].sum(axis=1) / X.sum(axis=1), tme_cont)[0] > 0.5


def test_pct_ribo_tracks_abundance_axis():
    on, truth = _cohort()
    off, _ = _cohort("technical_axes")
    ribo = _idx(on, RIBO)
    abundance = _idx(on, [g for g in truth["sets"]["abundance"] if g not in RIBO])
    long_genes = _idx(on, truth["sets"]["long_genes"])
    sample = on.obs["sample"].to_numpy()

    def rho(adata):
        X = _dense(adata)
        total = X.sum(axis=1)
        score = np.log1p(X / total[:, None] * 1e4)[:, abundance].mean(axis=1)
        return (_median_within_sample_spearman(X[:, ribo].sum(axis=1) / total, score, sample),
                _median_within_sample_spearman(X[:, long_genes].sum(axis=1) / total,
                                               adata.obs["intronic_fraction"].to_numpy(), sample))

    (ribo_on, long_on), (ribo_off, long_off) = rho(on), rho(off)
    # calibrated magnitudes (sim module docstring): 0.92 and 0.75 here; the first ones gave 0.75 and 0.55
    assert ribo_on > 0.85 and ribo_on > ribo_off + 0.2
    assert long_on > 0.65 and abs(long_off) < 0.15


def _shared_variance(L: np.ndarray, idx) -> float:
    """Sum over the genes ``idx`` of their mean covariance with the set's other genes (columns of
    ``L``): the variance a factor common to the set contributes, summed over its genes."""
    C = np.cov(L[:, idx], rowvar=False)
    return float((C.sum() - np.trace(C)) / (len(idx) - 1))


def test_technical_axes_rival_a_planted_program():
    """Cheap proxy (no NMF) for the calibration target in the sim module docstring: within each
    sample's malignant cells, each axis's shared variance in log2(CPM/10 + 1), the scale the pipeline
    factorises, relative to that of hypoxia (robust in 13 of 16 samples of the default cohort).

    Of the magnitudes tried on the default cohort, those whose ratios here fall below the bounds
    (CYTO_GAMMA 1.75: abundance 0.68; LONG_K 10: long genes 1.28) failed the target and those above
    them met it; the calibrated ones give 0.97 and 1.73. The first magnitudes gave 0.23 and 0.58,
    with neither axis in any NMF program. A guard against the axes slipping back below what NMF can
    see, not a substitute for the NMF check.
    """
    adata, truth = _cohort()
    X = _dense(adata)
    L = np.log2(_cpm(X) / 10 + 1)
    sample, mal = adata.obs["sample"].to_numpy(), adata.obs["malignant"].to_numpy(bool)
    sets = {name: _idx(adata, truth["sets"][name]) for name in ("hypoxia", "abundance", "long_genes")}
    ratios = {"abundance": [], "long_genes": []}
    for s in np.unique(sample):
        rows = mal & (sample == s)
        hypoxia = _shared_variance(L[rows], sets["hypoxia"])
        for name in ratios:
            ratios[name].append(_shared_variance(L[rows], sets[name]) / hypoxia)
    assert np.median(ratios["abundance"]) > 0.75 and np.median(ratios["long_genes"]) > 1.4, ratios


def test_axis_magnitudes_change_no_random_draw(monkeypatch):
    """``CYTO_GAMMA`` and ``LONG_K`` scale planted draws without changing them: with the first
    magnitudes every draw is as before, and the cohort without technical axes (hence the null
    cohort) is bit-identical."""
    on, t_on = _cohort()
    off, _ = _cohort("technical_axes")
    monkeypatch.setattr(sim, "CYTO_GAMMA", 1.0)
    monkeypatch.setattr(sim, "LONG_K", 2.0)
    first, t_first = simulate_cohort(**SMALL)
    first_off, _ = simulate_cohort(**SMALL, technical_axes=False)
    assert (off.layers["counts"] != first_off.layers["counts"]).nnz == 0
    pd.testing.assert_frame_equal(off.obs, first_off.obs)
    assert (on.layers["counts"] != first.layers["counts"]).nnz > 0  # the magnitudes do act
    drawn = [c for c in on.obs if not c.endswith("_est")]  # the estimates are computed from the counts
    pd.testing.assert_frame_equal(on.obs[drawn], first.obs[drawn])
    pd.testing.assert_frame_equal(t_on["cells"], t_first["cells"])
    assert t_on["sets"] == t_first["sets"] and t_on["private"] == t_first["private"]


def test_capture_bias_shifts_means_not_coexpression():
    on, truth = _cohort()
    off, _ = _cohort("protocol_capture_bias")
    X_on, X_off = _dense(on), _dense(off)
    fixed = (on.obs["protocol"] == "fixed").to_numpy() & on.obs["malignant"].to_numpy(bool)
    ok = X_off[fixed].mean(axis=0) >= 0.5
    lfc = np.log(_cpm(X_on)[fixed][:, ok].mean(axis=0) / _cpm(X_off)[fixed][:, ok].mean(axis=0))
    assert pearsonr(lfc, np.log(truth["capture_efficiency"].to_numpy()[ok]))[0] > 0.9
    sample = on.obs["sample"].to_numpy()
    # The kappa difference is estimation noise (100 cells; capture redraws the counts): over seeds 0-5
    # its mean is 0.01 and its SD 0.06 (0.045 before the technical axes were calibrated up).
    for s in np.unique(sample[fixed]):
        rows = fixed & (sample == s)
        for name in ["cycle", "hypoxia"]:
            idx = _idx(on, truth["sets"][name])
            k_on, k_off = (kappa(X[rows][:, idx], X[rows].sum(axis=1)) for X in (X_on, X_off))
            assert k_off > 0.3 and abs(k_on - k_off) < 0.15


def test_estimated_covariates_track_truth():
    adata, truth = _cohort()
    obs = adata.obs
    mal = obs["malignant"].to_numpy(bool)
    aris = [_ari(obs.loc[rows, "clone"], obs.loc[rows, "clone_est"])
            for s in obs["sample"].unique()
            if obs.loc[rows := mal & (obs["sample"] == s).to_numpy(), "tumour"].iloc[0] in truth["tumours_cna"]]
    assert np.mean(aris) > 0.2
    assert spearmanr(obs.loc[mal, "contamination"], obs.loc[mal, "contamination_est"])[0] > 0.5
    assert spearmanr(obs.loc[~mal, "contamination"], obs.loc[~mal, "contamination_est"])[0] > 0.5
    doublet = truth["cells"]["doublet"].to_numpy()
    assert obs.loc[doublet, "doublet_est"].median() > 0.3 > 0.1 > obs.loc[mal & ~doublet, "doublet_est"].median()


def test_estimates_use_only_counts_labels_and_positions():
    adata, truth = _cohort()
    arm = assign_arms(adata.var, centromeres(truth["cytoband"])).to_numpy(dtype=object)
    est = _estimate_covariates(adata.layers["counts"], adata.obs["sample"].to_numpy(),
                               adata.obs["cell_type"].to_numpy(), arm)
    pd.testing.assert_frame_equal(est.set_axis(adata.obs.index), adata.obs[list(est.columns)])


def test_cycling_raises_doublet_scores():
    on, truth = _cohort()
    off, _ = _cohort("doublet_cycle_coupling")
    mal = on.obs["malignant"].to_numpy(bool)
    cycling = truth["cells"]["cycle"].to_numpy()
    singlet = mal & ~truth["cells"]["doublet"].to_numpy()

    def gap(adata, col):
        return adata.obs.loc[singlet & cycling, col].mean() - adata.obs.loc[singlet & ~cycling, col].mean()

    assert gap(on, "doublet_score") > 0.2 and abs(gap(off, "doublet_score")) < 0.05
    totals = _dense(on).sum(axis=1)
    assert np.median(totals[singlet & cycling]) > 1.15 * np.median(totals[singlet & ~cycling])
    # the expression-based estimate inherits the coupling through library size
    r_on, r_off = (spearmanr(a.obs.loc[singlet, "doublet_est"], cycling[singlet])[0] for a in (on, off))
    assert r_on > 0.05 and r_on > r_off + 0.05
