"""A1 / TME co-variation: spillover subtraction, matched random-set null, rho_hat sensitivity, kappa_ratio."""
import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
import pytest

import mpaudit.tme as tme
from mpaudit import Keys
from mpaudit.scoring import control_genes
from mpaudit.sim import simulate_cohort
from mpaudit.tme import PER_SET_COLUMNS, tme_covariation

KEYS = Keys(sample="sample", malignant="malignant", cell_type="cell_type", cancer_type="cancer_type", tumour="tumour")
# 16 fresh samples: no protocol twins, whose capture bias moves every set in both compartments; 1,500 genes so
# that few stress genes are drawn as control genes (they would put part of the stress axis into every score)
COHORT = dict(n_tumours=16, n_twin=0, n_region=0, n_mal=100, n_tme_per_type=50, n_genes=1500, seed=0)
SETS = ("hypoxia", "CNA_7q", "stress", "IFN", "contamination", "cycle")
RHO_COLUMNS = ["partial_r2_rho", "slope_rho", "null_q95_rho", "null_q99_rho"]
NOT_ADJUSTED = "spillover is NOT adjusted"
_NOT_ADJUSTED_WARNINGS: dict[bool, list[str]] = {}


def _not_adjusted_warnings(fn, *args, **kwargs):
    """``fn(*args, **kwargs)`` and the messages of the 'spillover is NOT adjusted' warnings it raised."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = fn(*args, **kwargs)
    return out, [str(w.message) for w in caught if NOT_ADJUSTED in str(w.message)]


@lru_cache(maxsize=None)
def _cohort(malignant_ambient: bool):
    return simulate_cohort(**COHORT, malignant_ambient=malignant_ambient)


@lru_cache(maxsize=None)
def _run(malignant_ambient: bool):
    """A1 tables on the cohort with or without malignant->TME ambient RNA (shared; do not modify)."""
    adata, truth = _cohort(malignant_ambient)
    (per_set, per_sample), _NOT_ADJUSTED_WARNINGS[malignant_ambient] = _not_adjusted_warnings(
        tme_covariation, adata, KEYS, {k: truth["sets"][k] for k in SETS}, counts_layer="counts",
        n_null=200, n_boot=50, kappa_ratio=False)
    return per_set.set_index("mp"), per_sample


def _effective(r2, slope):
    """Partial R^2 counted as 0 for a non-positive slope, as in the null."""
    return np.where(np.asarray(slope) > 0, r2, 0.0)


def _exceeds_q99(row) -> bool:
    return bool(row["slope"] > 0 and row["partial_r2"] > row["null_q99"])


def _exceeds_q99_rho(row) -> bool:
    """The same test in the ``rho_sensitivity`` re-run."""
    return bool(row["slope_rho"] > 0 and row["partial_r2_rho"] > row["null_q99_rho"])


def test_tables_have_the_documented_columns():
    per_set, per_sample = _run(True)
    assert list(per_set.reset_index().columns) == PER_SET_COLUMNS
    assert (per_set["n_samples"] == 16).all() and per_set["spill_adjusted"].all()
    for k in SETS:
        assert {f"{k}::malignant", f"{k}::tme", f"{k}::tme_spill"} <= set(per_sample.columns)
    assert {"rho_hat::T", "rho_hat::Myeloid", "rho_hat::Fibroblast"} <= set(per_sample.columns)
    assert ((per_set["null_q95"] <= per_set["null_q99"]) & per_set["null_p"].between(1 / 201, 1)).all()
    assert (per_set["partial_r2_lo"] <= per_set["partial_r2_hi"]).all()
    assert np.isfinite(per_set[RHO_COLUMNS].to_numpy(float)).all()
    assert (per_set["null_q95_rho"] <= per_set["null_q99_rho"]).all()
    assert not _NOT_ADJUSTED_WARNINGS[True]  # 10 planted lineage markers qualify


def test_spillover_inflates_unadjusted_a1_and_subtraction_removes_it():
    on, _ = _run(True)
    off, _ = _run(False)
    unadj_on = _effective(on["partial_r2_unadjusted"], on["slope_unadjusted"])
    unadj_off = _effective(off["partial_r2_unadjusted"], off["slope_unadjusted"])
    hyp, cna = SETS.index("hypoxia"), SETS.index("CNA_7q")
    # malignant-only programs whose strength varies between samples leak into the TME (seed 0: hypoxia
    # 0.22 vs 0.02, CNA 0.73 vs 0.00 without ambient RNA) ...
    assert unadj_on[hyp] > unadj_off[hyp] + 0.1
    assert unadj_on[cna] > unadj_off[cna] + 0.4
    # ... and are not TME co-variation once the expected leak is subtracted
    assert not _exceeds_q99(on.loc["hypoxia"]) and not _exceeds_q99(on.loc["CNA_7q"])
    assert on.loc["CNA_7q", "partial_r2"] < unadj_on[cna] - 0.4
    # shared by both compartments: dissociation stress (technical) and IFN (TME-shared biology) survive
    for mp in ("stress", "IFN"):
        assert _exceeds_q99(on.loc[mp]) and on.loc[mp, "null_p"] <= 0.01
        assert on.loc[mp, "partial_r2"] > 0.6
    for mp in ("contamination", "cycle"):
        assert not _exceeds_q99(on.loc[mp])
    # at 1.5 x rho_hat the shared programs stay above their null and the leaks stay below it (seed 0: stress
    # 0.95 vs 0.34, IFN 0.84 vs 0.36; CNA and hypoxia over-subtracted to negative slopes)
    for mp in ("stress", "IFN"):
        assert _exceeds_q99_rho(on.loc[mp])
    for mp in ("hypoxia", "CNA_7q", "contamination", "cycle"):
        assert not _exceeds_q99_rho(on.loc[mp])


def test_rho_hat_tracks_the_planted_ambient_fraction():
    adata, _ = _cohort(True)
    _, per_sample = _run(True)
    tme_obs = adata.obs[~adata.obs["malignant"]]
    planted = tme_obs.groupby(["sample", "cell_type"], observed=True)["contamination"].mean().unstack()
    for ct in ("T", "Myeloid", "Fibroblast"):
        est = per_sample.set_index("sample")[f"rho_hat::{ct}"]
        true = planted.loc[est.index, ct]
        assert np.corrcoef(est, true)[0, 1] > 0.8
        # markers carry a little native TME expression (0.1x), so rho_hat errs upward, not downward
        assert 0.9 < float(np.median(est / true)) < 1.6


def test_without_markers_nothing_is_adjusted():
    off, per_sample = _run(False)  # no malignant lineage genes without malignant ambient RNA
    assert not off["spill_adjusted"].any()
    # spillover_adjust was requested, so the fallback to the unadjusted statistic is not silent
    assert len(_NOT_ADJUSTED_WARNINGS[False]) == 1 and "criterion (b)" in _NOT_ADJUSTED_WARNINGS[False][0]
    assert off[RHO_COLUMNS].isna().all().all()
    np.testing.assert_allclose(off["partial_r2"], off["partial_r2_unadjusted"])
    np.testing.assert_allclose(off["slope"], off["slope_unadjusted"])
    assert per_sample["hypoxia::tme_spill"].isna().all()
    assert not any(c.startswith("rho_hat::") for c in per_sample)


def test_unadjusted_run_reproduces_the_reference_columns():
    adata, truth = _cohort(True)
    on, _ = _run(True)
    (plain, _), warned = _not_adjusted_warnings(
        tme_covariation, adata, KEYS, {k: truth["sets"][k] for k in SETS}, counts_layer="counts",
        spillover_adjust=False, n_null=0, n_boot=0, kappa_ratio=False)
    plain = plain.set_index("mp")
    assert not warned  # switched off on purpose: no warning
    assert not plain["spill_adjusted"].any() and plain["null_q99"].isna().all()
    assert plain[RHO_COLUMNS].isna().all().all()
    np.testing.assert_allclose(plain["partial_r2"], on["partial_r2_unadjusted"])
    np.testing.assert_allclose(plain["slope"], on["slope_unadjusted"])


def test_audited_genes_are_not_used_as_markers():
    adata, truth = _cohort(True)
    lineage = {"lineage": truth["malignant_markers"]}
    kw = dict(counts_layer="counts", n_null=0, n_boot=0, kappa_ratio=False)
    # the only malignant-specific genes are the audited set's own: excluded by default, so no adjustment
    with pytest.warns(UserWarning, match=rf"only \d gene.*{NOT_ADJUSTED}.*criterion \(b\)"):
        assert not tme_covariation(adata, KEYS, lineage, **kw)[0]["spill_adjusted"].any()
    adjusted, warned = _not_adjusted_warnings(tme_covariation, adata, KEYS, lineage, exclude_genes=[], **kw)
    assert adjusted[0]["spill_adjusted"].all() and not warned


def test_rho_sensitivity_reruns_the_adjusted_statistic_at_scaled_rho_hat(monkeypatch):
    adata, truth = _cohort(True)
    sets = {k: truth["sets"][k] for k in SETS}
    kw = dict(counts_layer="counts", n_null=40, n_boot=0, kappa_ratio=False)
    ref = tme_covariation(adata, KEYS, sets, rho_sensitivity=1.5, **kw)[0]
    main = [c for c in PER_SET_COLUMNS if c not in RHO_COLUMNS]
    for off in (None, 1, 1.0):
        other = tme_covariation(adata, KEYS, sets, rho_sensitivity=off, **kw)[0]
        assert other[RHO_COLUMNS].isna().all().all()
        pd.testing.assert_frame_equal(other[main], ref[main])  # the re-run draws nothing
    for bad in (0, -1.5, float("nan")):
        with pytest.raises(ValueError, match="rho_sensitivity"):
            tme_covariation(adata, KEYS, sets, rho_sensitivity=bad, **kw)
    # identical to a plain run whose rho_hat is 1.5 times larger, null included (same random sets and controls)
    rho_hat = tme._rho_hat
    monkeypatch.setattr(tme, "_rho_hat", lambda *a: 1.5 * rho_hat(*a))
    scaled = tme_covariation(adata, KEYS, sets, rho_sensitivity=None, **kw)[0]
    for col in ("partial_r2", "slope", "null_q95", "null_q99"):
        np.testing.assert_allclose(ref[f"{col}_rho"], scaled[col], rtol=1e-12, atol=1e-12)
    assert not np.allclose(ref["partial_r2_rho"], ref["partial_r2"])


def test_rho_sensitivity_catches_a_leak_left_by_rho_hat_biased_down(monkeypatch):
    adata, truth = _cohort(True)
    on, _ = _run(True)
    rho_hat = tme._rho_hat
    monkeypatch.setattr(tme, "_rho_hat", lambda *a: 0.7 * rho_hat(*a))  # under-subtraction by 30%
    low = tme_covariation(adata, KEYS, {k: truth["sets"][k] for k in SETS}, counts_layer="counts",
                          n_null=200, n_boot=0, kappa_ratio=False)[0].set_index("mp")
    # seed 0: CNA_7q (malignant only) 0.41 vs null q99 0.27, against 0.15 vs 0.24 at the estimated rho_hat: the
    # residual leak raises the MP's partial R^2 but hardly the null's, a false criterion-(b) flag ...
    assert not _exceeds_q99(on.loc["CNA_7q"]) and _exceeds_q99(low.loc["CNA_7q"])
    assert low.loc["CNA_7q", "partial_r2"] > on.loc["CNA_7q", "partial_r2"] + 0.15
    assert low.loc["CNA_7q", "null_q99"] < on.loc["CNA_7q", "null_q99"] + 0.1
    # ... that the re-run at 1.5 x the biased rho_hat does not repeat (0.11 vs null_q99_rho 0.23), while the
    # shared programs still exceed it (stress 0.96 vs 0.39, IFN 0.85 vs 0.41)
    assert not _exceeds_q99_rho(low.loc["CNA_7q"])
    assert low.loc["CNA_7q", "partial_r2_rho"] < low.loc["CNA_7q", "null_q99_rho"]
    for mp in ("stress", "IFN"):
        assert _exceeds_q99(low.loc[mp]) and _exceeds_q99_rho(low.loc[mp])


def test_null_match_options():
    adata, truth = _cohort(True)
    kw = dict(counts_layer="counts", n_null=40, n_boot=0, kappa_ratio=False)
    sets = {"stress": truth["sets"]["stress"], "tiny": truth["sets"]["stress"][:3]}
    per_set, _ = tme_covariation(adata, KEYS, sets, null_match="expression", **kw)
    assert per_set.loc[0, "null_p"] == pytest.approx(1 / 41)
    assert per_set.loc[1, "n_samples"] == 0 and np.isnan(per_set.loc[1, "partial_r2"])  # < 5 genes
    with pytest.raises(ValueError, match="null_match"):
        tme_covariation(adata, KEYS, sets, null_match="specificity", **kw)


def test_matched_random_sets_and_control_genes():
    rng = np.random.default_rng(0)
    expr, spec = rng.normal(size=600), rng.normal(size=600)
    bins = tme._bins(expr, 10)
    bin_pools = tme._pools(bins)
    match = np.empty(600, dtype=int)
    for b, members in bin_pools.items():
        match[members] = b * 3 + tme._bins(spec[members], 3)
    idx = np.sort(rng.choice(600, 25, replace=False))
    r = tme._matched_random_set(idx, match, tme._pools(match), bins, bin_pools, np.random.default_rng(1))
    assert r.size == idx.size and not set(r) & set(idx)
    assert sorted(match[r]) == sorted(match[idx])  # joint bins of 20 genes are not used up by 25 draws
    ours = tme._control_genes(idx, bins, bin_pools, np.random.default_rng(2), 100)
    np.testing.assert_array_equal(ours, control_genes(idx, expr, np.random.default_rng(2), 10, 100))


def test_kappa_ratio_subsamples_at_random(monkeypatch):
    adata, truth = simulate_cohort(n_tumours=3, n_twin=0, n_region=0, n_mal=60, n_tme_per_type=40, n_genes=500, seed=1)
    sets = {"stress": truth["sets"]["stress"]}
    kw = dict(counts_layer="counts", n_boot=0, max_kappa_cells=25)
    a = tme_covariation(adata, KEYS, sets, n_null=0, **kw)[0]
    b = tme_covariation(adata, KEYS, sets, n_null=5, **kw)[0]
    assert a.loc[0, "kappa_ratio_n"] == 3 and np.isfinite(a.loc[0, "kappa_ratio_median"])
    assert a.loc[0, "kappa_ratio_median"] == b.loc[0, "kappa_ratio_median"]  # own random stream

    seen = []
    monkeypatch.setattr(tme, "kappa_with_null", lambda c, idx, n_null: seen.append(c) or {"kappa": 0.5})
    tme_covariation(adata, KEYS, sets, n_null=0, **kw)
    counts = adata.layers["counts"]
    first = [counts[np.flatnonzero((adata.obs["sample"] == s) & (adata.obs["cell_type"] == ct))[:25]]
             for s in adata.obs["sample"].unique() for ct in adata.obs["cell_type"].unique()]
    assert len(seen) == 12 and all(c.shape[0] == 25 for c in seen)
    assert not any(any((c != f).nnz == 0 for f in first) for c in seen)  # not the first cells of any group
