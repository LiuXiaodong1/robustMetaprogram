"""Axis A4 (mpaudit.protocol): paired Δlog κ, protocol-specific intensity controls and the two nulls."""
from functools import lru_cache

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy import stats

from mpaudit.data import Keys, get_counts, log2_cpm10, malignant_mask
from mpaudit.protocol import (KAPPA_FLOOR, KAPPA_MIN_EFFECT, per_sample_kappa, protocol_shift, protocol_shift_kappa,
                              replicate_bound, tumour_shifts)
from mpaudit.scoring import score_gene_sets
from mpaudit.sim import simulate_cohort

KEYS = Keys(sample="sample", tumour="tumour", protocol="protocol")
# 4 paired tumours (fresh + fixed) and 3 tumours with two fresh pieces (3 replicate pairs);
# protocol_capture_bias is on (the default)
SMALL = dict(n_tumours=7, n_twin=4, n_region=3, n_mal=120, n_tme_per_type=10, n_genes=800, seed=0)


@lru_cache(maxsize=None)
def _cohort():
    """The small cohort and the sets audited on it (shared across tests; do not modify)."""
    adata, truth = simulate_cohort(**SMALL)
    sets = {k: truth["sets"][k] for k in ("cycle", "hypoxia", "stress")}
    # the cycle genes captured worst under the fixed protocol: a program the protocol leaves alone,
    # whose genes nevertheless look lower under it
    cap = truth["capture_efficiency"]
    sets["cycle_low_capture"] = list(cap[truth["sets"]["cycle"]].sort_values().index[:12])
    return adata, sets


@lru_cache(maxsize=None)
def _shift(min_pairs: int):
    adata, sets = _cohort()
    return protocol_shift_kappa(adata, KEYS, sets, "fresh", "fixed", min_pairs=min_pairs, n_null_sets=40)


def _levels(paired: dict[str, float], replicates: dict[str, float], metric: str = "intensity",
            base: float = 0.0, scale=lambda b, d: b + d) -> pd.DataFrame:
    """One MP: per tumour a fresh sample at ``base`` and a fixed (paired) or second fresh (replicate) sample."""
    rows = []
    for t, d in paired.items():
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", metric: base},
                 {"tumour": t, "protocol": "fixed", "mp": "m", metric: scale(base, d)}]
    for t, d in replicates.items():
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", metric: base},
                 {"tumour": t, "protocol": "fresh", "mp": "m", metric: scale(base, d)}]
    return pd.DataFrame(rows)


T95_3 = stats.t.ppf(0.95, 3)  # 2.353: prediction bound factor for 3 replicate tumours


def test_loss_is_one_sided_against_replicates():
    lv = _levels({"A": -3.0, "B": 2.5, "C": 0.1, "D": -0.1}, {"E": 1.0, "F": -1.0, "G": 1.0})
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["null_type"] == "replicates" and res["n_null_pairs"] == 3 and res["n_null_tumours"] == 3
    assert res["n_paired_tumours"] == 4 and res["null_type_intensity"] == "replicates"
    sh = tumour_shifts(lv, "fresh", "fixed")
    assert sh["null_hi_intensity"].to_numpy() == pytest.approx(T95_3)  # t_{0.95,3} * sqrt(mean d^2), d = +-1
    assert sh["null_lo_intensity"].to_numpy() == pytest.approx(-T95_3)
    assert res["frac_loss_beyond_null_intensity"] == 0.25  # only A drops below -2.35
    assert res["frac_gain_beyond_null_intensity"] == 0.25  # only B
    assert res["mean_delta_intensity"] == pytest.approx(-0.125)
    rev = protocol_shift(lv, "fixed", "fresh").iloc[0]
    assert rev["mean_delta_intensity"] == pytest.approx(0.125)
    assert rev["frac_loss_beyond_null_intensity"] == 0.25 and rev["frac_gain_beyond_null_intensity"] == 0.25


def test_kappa_compared_on_floored_log_scale():
    lv = _levels({"A": 0.002, "B": 0.4}, {"C": 0.75, "D": 0.75, "E": 1 / 0.75}, metric="kappa", base=0.4,
                 scale=lambda b, f: b * f)
    sh = tumour_shifts(lv, "fresh", "fixed").set_index("tumour")
    assert sh.loc["A", "delta_log_kappa"] == pytest.approx(np.log(KAPPA_FLOOR / 0.4))  # 0.0008 floored to 0.01
    assert sh.loc["B", "delta_log_kappa"] == pytest.approx(np.log(0.4))
    assert sh["null_hi_kappa"].to_numpy() == pytest.approx(T95_3 * np.log(1 / 0.75))  # replicates on the log scale
    assert sh["loss_kappa"].tolist() == [True, True] and sh["gain_kappa"].tolist() == [False, False]
    # several samples of a protocol: kappa is averaged before the log
    two = pd.concat([lv, pd.DataFrame([{"tumour": "B", "protocol": "fixed", "mp": "m", "kappa": 0.4}])])
    sh2 = tumour_shifts(two, "fresh", "fixed").set_index("tumour")
    assert sh2.loc["B", "delta_log_kappa"] == pytest.approx(np.log(0.28 / 0.4))


def test_kappa_minimum_effect_margin():
    """A κ change inside +-log 1.5 is never a loss or gain, however tight the null."""
    lv = _levels({"A": 1 / 1.3, "B": 1 / 1.6, "C": 1.3, "D": 1.6}, {"E": 1.01, "F": 0.99, "G": 1.01},
                 metric="kappa", base=0.5, scale=lambda b, f: b * f)
    sh = tumour_shifts(lv, "fresh", "fixed").set_index("tumour")
    assert (sh["null_hi_kappa"] < 0.03).all()  # the replicates alone would flag every tumour
    assert sh["loss_kappa"].tolist() == [False, True, False, False]
    assert sh["gain_kappa"].tolist() == [False, False, False, True]
    none = tumour_shifts(lv, "fresh", "fixed", kappa_min_effect=0.0).set_index("tumour")
    assert none["loss_kappa"].tolist() == [True, True, False, False]
    assert KAPPA_MIN_EFFECT == pytest.approx(np.log(1.5))
    # intensity has no margin: a change of 0.3 beyond a tight null counts
    li = _levels({"A": -0.3, "B": 0.3}, {"E": 0.01, "F": -0.01, "G": 0.01})
    assert tumour_shifts(li, "fresh", "fixed")["loss_intensity"].tolist() == [True, False]


def _random_null(tumours: str, deltas: np.ndarray, metric: str = "intensity", base: float = 0.0) -> pd.DataFrame:
    rows = []
    for t in tumours:
        for j, d in enumerate(deltas):
            rows += [{"mp": "m", "tumour": t, "protocol": "fresh", "set": j, metric: base},
                     {"mp": "m", "tumour": t, "protocol": "fixed", "set": j, metric: base + d}]
    return pd.DataFrame(rows)


def test_random_set_null_when_too_few_replicate_pairs():
    lv = _levels({"A": -3.0, "B": -0.5}, {"E": 1.0, "F": -1.0})  # 2 replicate pairs < min_pairs
    none = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert none["null_type"] == "none" and none["n_paired_tumours"] == 0
    assert np.isnan(none["frac_loss_beyond_null_intensity"])
    null = _random_null("AB", np.linspace(-1, 1, 41))
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=null)
    assert (sh["null_type"] == "random_sets").all()
    assert sh["null_lo_intensity"].to_numpy() == pytest.approx(-0.9)  # each tumour's 5th percentile
    assert sh["loss_intensity"].tolist() == [True, False]
    res = protocol_shift(lv, "fresh", "fixed", null_levels=null).iloc[0]
    assert res["frac_loss_beyond_null_intensity"] == 0.5 and res["n_paired_tumours"] == 2
    # too few random sets in a tumour: that tumour is not tested
    short = pd.concat([null[null["tumour"] == "A"], _random_null("B", np.linspace(-1, 1, 10))])
    assert tumour_shifts(lv, "fresh", "fixed", null_levels=short)["loss_intensity"].isna().tolist() == [False, True]
    # with min_pairs replicate pairs the replicates are the null
    lv3 = _levels({"A": -3.0, "B": -0.5}, {"E": 1.0, "F": -1.0, "G": 0.2})
    assert protocol_shift(lv3, "fresh", "fixed", null_levels=null).iloc[0]["null_type"] == "replicates"


def test_replicate_null_counts_distinct_tumours():
    # three pairs from one tumour with three pieces share pieces: one tumour, not three
    rows = [{"tumour": t, "protocol": p, "mp": "m", "intensity": v}
            for t in "ABCD" for p, v in (("fresh", 0.0), ("fixed", -3.0))]
    rows += [{"tumour": "E", "protocol": "fresh", "mp": "m", "intensity": v} for v in (0.0, 1.0, -1.0)]
    lv = pd.DataFrame(rows)
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["n_null_pairs"] == 3 and res["n_null_tumours"] == 1
    assert res["null_type"] == "none" and res["n_paired_tumours"] == 0
    assert protocol_shift(lv, "fresh", "fixed", min_pairs=1).iloc[0]["null_type"] == "replicates"
    # the bound: t with df = distinct tumours, s = sqrt(mean d^2) over every pair (1, -1, -2)
    sh = tumour_shifts(lv, "fresh", "fixed", min_pairs=1)
    assert sh["null_hi_intensity"].to_numpy() == pytest.approx(stats.t.ppf(0.95, 1) * np.sqrt(2.0))
    assert replicate_bound(np.array([1.0, -1.0, np.nan]), 2) == pytest.approx(stats.t.ppf(0.95, 2))
    assert np.isnan(replicate_bound(np.array([np.nan]), 3)) and np.isnan(replicate_bound(np.array([1.0]), 0))


def test_replicate_null_is_calibrated_with_three_tumours():
    """No protocol effect, paired and replicate changes exchangeable: about 5% of tumours are losses.

    The empirical 95th percentile of 3 absolute replicate changes (the previous rule) gives about
    13% per tumour and 11% for >= 2 of 4 tumours in the same simulation.
    """
    rng = np.random.default_rng(1)
    n_mp = 2500
    paired, reps = rng.normal(size=(n_mp, 4)), rng.normal(size=(n_mp, 3))
    mp = np.repeat([f"m{i}" for i in range(n_mp)], 14)
    tumour = np.tile(np.repeat(list("ABCDEFG"), 2), n_mp)
    protocol = np.tile(["fresh", "fixed"] * 4 + ["fresh"] * 6, n_mp)
    value = np.zeros((n_mp, 14))
    value[:, 1:8:2], value[:, 9::2] = paired, reps
    lv = pd.DataFrame({"mp": mp, "tumour": tumour, "protocol": protocol, "intensity": value.ravel()})
    sh = tumour_shifts(lv, "fresh", "fixed")
    assert (sh["null_type"] == "replicates").all()
    per_tumour = sh["loss_intensity"].astype(bool).mean()
    two_of_four = (sh.groupby("mp")["loss_intensity"].sum() >= 2).mean()
    assert 0.04 < per_tumour < 0.06
    assert two_of_four < 0.05
    old_rule = (paired < -np.quantile(np.abs(reps), 0.95, axis=1)[:, None]).mean()
    assert old_rule > 0.11


def test_collapsed_random_null_is_untestable():
    lv = _levels({"A": 0.1, "B": 0.1}, {}, metric="kappa", base=0.5, scale=lambda b, f: b * f)
    null = _random_null("AB", np.zeros(30), metric="kappa", base=0.001)  # every random set changes identically
    res = protocol_shift(lv, "fresh", "fixed", null_levels=null).iloc[0]
    assert res["null_type"] == "random_sets" and res["n_paired_tumours"] == 0
    assert res["mean_delta_log_kappa"] == pytest.approx(np.log(0.1))
    assert np.isnan(res["frac_loss_beyond_null_kappa"])


def _kappa_null(rng, kappa_ref: float, kappa_test: float, sd: float, mp_ref: float, mp_test: float,
                tumours: str = "ABCD", n_sets: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One MP at ``mp_ref -> mp_test`` and random sets at ``kappa_ref -> kappa_test`` (+- ``sd``) per tumour."""
    lv, nl = [], []
    for t in tumours:
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": mp_ref + rng.normal(0, 0.003)},
               {"mp": "m", "tumour": t, "protocol": "fixed", "kappa": mp_test + rng.normal(0, 0.003)}]
        for j in range(n_sets):
            for p, k in (("fresh", kappa_ref), ("fixed", kappa_test)):
                nl.append({"mp": "m", "tumour": t, "protocol": p, "set": j, "kappa": k + rng.normal(0, sd)})
    return pd.DataFrame(lv), pd.DataFrame(nl)


def test_random_set_kappa_null_moved_to_the_mp_level():
    rng = np.random.default_rng(1)
    # background κ rises additively (0.012 -> 0.027) and the MP is unchanged at 0.30: on the log
    # scale of the random sets this was a loss in every tumour; carried to κ = 0.30 it is not
    lv, nl = _kappa_null(rng, 0.012, 0.027, 0.006, 0.30, 0.30)
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl)
    assert (sh["null_type"] == "random_sets").all() and sh["loss_kappa"].notna().all()
    assert not sh["loss_kappa"].any() and not sh["gain_kappa"].any()
    lin = nl.pivot_table(index=["tumour", "set"], columns="protocol", values="kappa")
    q05 = (lin["fixed"] - lin["fresh"]).groupby(level="tumour").quantile(0.05)
    base = sh.set_index("tumour")["kappa_reference"]
    assert sh.set_index("tumour")["null_lo_kappa"].to_numpy() == pytest.approx(np.log((base + q05) / base).to_numpy())
    assert (sh["null_lo_kappa"].abs() < 0.05).all()  # the log-scale null had q05 between +0.05 and +0.3
    # random sets near the floor spread widely on the log scale and hid a twofold loss; not any more
    lv, nl = _kappa_null(rng, 0.012, 0.012, 0.01, 0.60, 0.30)
    res = protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert res["n_paired_tumours"] == 4 and res["frac_loss_beyond_null_kappa"] == 1.0
    # an MP at the floor under the reference cannot lose κ
    lv, nl = _kappa_null(rng, 0.012, 0.012, 0.01, 0.005, 0.001)
    assert protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]["frac_loss_beyond_null_kappa"] == 0


def test_null_type_per_statistic_ignores_nan_kappa_pairs():
    """κ NaN in the replicate tumours: κ falls back to random sets, intensity keeps its replicates."""
    rng = np.random.default_rng(0)
    lv, nl = [], []
    for t in "ABCD":
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": 0.6, "intensity": 0.0},
               {"mp": "m", "tumour": t, "protocol": "fixed", "kappa": 0.05, "intensity": 0.0}]
        for j in range(100):
            for p in ("fresh", "fixed"):
                nl.append({"mp": "m", "tumour": t, "protocol": p, "set": j, "kappa": 0.05 + rng.normal(0, 0.02),
                           "intensity": rng.normal(0, 0.1)})
    for t in "EFG":  # panel too small in the replicate tumours: κ NaN, intensity finite
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": np.nan, "intensity": 0.0},
               {"mp": "m", "tumour": t, "protocol": "fresh", "kappa": np.nan, "intensity": 0.1}]
    lv, nl = pd.DataFrame(lv), pd.DataFrame(nl)
    res = protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert res["null_type"] == "random_sets" and res["null_type_kappa"] == "random_sets"
    assert res["n_null_pairs"] == 0 and res["n_null_tumours"] == 0
    assert res["n_paired_tumours"] == 4 and res["frac_loss_beyond_null_kappa"] == 1.0
    assert res["null_type_intensity"] == "replicates"
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl)
    assert sh["null_hi_intensity"].to_numpy() == pytest.approx(T95_3 * 0.1)
    assert (sh["null_type_intensity"] == "replicates").all() and (sh["null_type_kappa"] == "random_sets").all()


def test_random_sets_drawn_for_mps_without_finite_kappa_replicates():
    adata, sets = _cohort()
    mod = adata.copy()
    X = mod.X.tolil()
    rep_cells = np.flatnonzero(mod.obs["tumour"].isin(["T05", "T06", "T07"]).to_numpy())
    cols = [mod.var_names.get_loc(g) for g in sets["stress"]]
    for c in cols:  # stress genes undetected in the replicate tumours: κ NaN there, intensity finite
        X[rep_cells, c] = 0
    mod.X = X.tocsr()
    per_mp, per_tumour, levels = protocol_shift_kappa(mod, KEYS, {"stress": sets["stress"], "cycle": sets["cycle"]},
                                                      "fresh", "fixed", n_null_sets=40)
    assert levels.loc[levels["tumour"].isin(["T05", "T06", "T07"]) & (levels["mp"] == "stress"), "kappa"].isna().all()
    per_mp = per_mp.set_index("mp")
    stress = per_mp.loc["stress"]
    assert stress["null_type"] == "random_sets" and stress["null_type_intensity"] == "replicates"
    assert stress["n_paired_tumours"] == 4 and stress["frac_loss_beyond_null_kappa"] == 1.0
    assert per_mp.loc["cycle", "null_type"] == "replicates" and per_mp.loc["cycle", "n_null_tumours"] == 3
    assert per_mp.loc["cycle", "frac_loss_beyond_null_kappa"] == 0


def test_kappa_flags_stress_not_cycle_with_replicate_null():
    per_mp = _shift(3)[0].set_index("mp")
    assert (per_mp["null_type"] == "replicates").all() and (per_mp["n_null_pairs"] == 3).all()
    assert (per_mp["n_null_tumours"] == 3).all() and (per_mp["null_type_intensity"] == "replicates").all()
    assert (per_mp["n_paired_tumours"] == 4).all() and (per_mp["iota"] > 0.9).all()
    assert per_mp.loc["stress", "frac_loss_beyond_null_kappa"] == 1.0  # stress is absent under fixed
    assert per_mp.loc["stress", "mean_delta_log_kappa"] < -np.log(2)
    for mp in ("cycle", "hypoxia", "cycle_low_capture"):
        assert per_mp.loc[mp, "frac_loss_beyond_null_kappa"] == 0
        assert abs(per_mp.loc[mp, "mean_delta_log_kappa"]) < np.log(1.5)
    for mp in ("cycle", "hypoxia"):  # unaffected by protocol: no intensity loss with protocol-specific controls
        assert per_mp.loc[mp, "frac_loss_beyond_null_intensity"] == 0


def test_random_set_null_used_below_min_pairs():
    per_mp, per_tumour, _ = _shift(4)  # 3 replicate pairs < 4
    per_mp = per_mp.set_index("mp")
    assert (per_mp["null_type"] == "random_sets").all() and (per_mp["n_paired_tumours"] == 4).all()
    assert (per_tumour["null_lo_kappa"] < per_tumour["null_hi_kappa"]).all()  # per-tumour random-set bounds
    assert per_mp.loc["stress", "frac_loss_beyond_null_kappa"] >= 0.5
    for mp in ("cycle", "hypoxia", "cycle_low_capture"):
        assert per_mp.loc[mp, "frac_loss_beyond_null_kappa"] == 0
    assert per_mp.loc["cycle", "frac_loss_beyond_null_intensity"] == 0
    # hypoxic fractions differ between pieces, which random sets do not model; still below the P1 bar
    assert per_mp.loc["hypoxia", "frac_loss_beyond_null_intensity"] < 0.5
    assert per_mp["frac_loss_beyond_null_prevalence"].isna().all()  # no random-set null for prevalence


def test_protocol_specific_controls_absorb_gene_wise_capture():
    adata, sets = _cohort()
    mal = malignant_mask(adata, KEYS)
    # with controls pooled across protocols, the low-capture cycle genes look lost under fixed ...
    score = score_gene_sets(log2_cpm10(get_counts(adata)), adata.var_names, {"m": sets["cycle_low_capture"]},
                            np.flatnonzero(mal))["m"].to_numpy()
    obs = adata.obs[mal]
    pooled = (pd.DataFrame({"sample": obs["sample"].to_numpy(), "tumour": obs["tumour"].to_numpy(),
                            "protocol": obs["protocol"].to_numpy(), "intensity": score[mal]})
              .groupby(["sample", "tumour", "protocol"], as_index=False)["intensity"].mean().assign(mp="m"))
    assert protocol_shift(pooled, "fresh", "fixed").iloc[0]["frac_loss_beyond_null_intensity"] >= 0.5
    # ... not with protocol-specific controls, and never in kappa
    res = _shift(3)[0].set_index("mp").loc["cycle_low_capture"]
    assert res["frac_loss_beyond_null_intensity"] < 0.5 and res["frac_loss_beyond_null_kappa"] == 0


def test_panel_and_iota():
    adata, sets = _cohort()
    cycle = sets["cycle"]
    fake = [f"NOT_A_GENE{i}" for i in range(20)]
    gs = {"cycle": cycle, "half": cycle[:10] + fake[:10], "mostly_missing": cycle[:4] + fake[:6]}
    base = per_sample_kappa(adata, KEYS, gs, "fresh", "fixed")
    assert (base["kappa"].notna() == ((base["iota"] >= 0.5) & (base["n_genes"] >= 5))).all()
    assert (base.loc[base["mp"] == "half", "iota"] <= 0.5).all()
    assert base.loc[base["mp"] == "mostly_missing", "kappa"].isna().all()
    res = protocol_shift(base, "fresh", "fixed").set_index("mp")
    assert res.loc["mostly_missing", "n_paired_tumours"] == 0 and res.loc["mostly_missing", "iota"] <= 0.4
    assert res.loc["cycle", "n_paired_tumours"] == 4

    # a gene no longer detected under fixed in T01 leaves T01's panel only
    X = adata.X.toarray()
    j = int(np.argmax(X[:, [adata.var_names.get_loc(g) for g in cycle]].sum(axis=0)))
    gene = cycle[j]
    t01_fixed = ((adata.obs["tumour"] == "T01") & (adata.obs["protocol"] == "fixed")).to_numpy()
    X[t01_fixed, adata.var_names.get_loc(gene)] = 0
    mod = adata.copy()
    mod.X = sp.csr_matrix(X)
    b = base[base["mp"] == "cycle"].set_index("sample")["n_genes"]
    m = per_sample_kappa(mod, KEYS, {"cycle": cycle}, "fresh", "fixed").set_index("sample")["n_genes"]
    changed = (b - m)[lambda d: d != 0]
    assert changed.to_dict() == {"T01_fresh_R1": 1, "T01_fixed_R1": 1}


def test_levels_match_and_cell_cap():
    adata, sets = _cohort()
    lv = per_sample_kappa(adata, KEYS, sets, "fresh", "fixed")
    pd.testing.assert_frame_equal(lv, _shift(3)[2])
    assert set(lv["protocol"]) == {"fresh", "fixed"} and (lv["n_cells"] == SMALL["n_mal"]).all()
    capped = per_sample_kappa(adata, KEYS, sets, "fresh", "fixed", max_cells=60)
    assert capped["kappa"].notna().all() and not np.allclose(capped["kappa"], lv["kappa"])
    pd.testing.assert_series_equal(capped["intensity"], lv["intensity"])  # intensity uses every cell


def test_input_errors():
    adata, sets = _cohort()
    with pytest.raises(ValueError, match="protocol"):
        per_sample_kappa(adata, Keys(sample="sample", tumour="tumour"), sets, "fresh", "fixed")
    with pytest.raises(ValueError, match="differ"):
        per_sample_kappa(adata, KEYS, sets, "fresh", "fresh")
    mixed = adata.copy()
    mixed.obs["sample"] = mixed.obs["tumour"].astype(str)  # T01 now holds fresh and fixed cells
    with pytest.raises(ValueError, match="both protocols"):
        per_sample_kappa(mixed, KEYS, sets, "fresh", "fixed")
