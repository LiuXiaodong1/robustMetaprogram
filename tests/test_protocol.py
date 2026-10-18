"""Axis A4 (mpaudit.protocol, P1-GO v1.3): paired Δlog κ, protocol-specific intensity controls, the two nulls, the
self bound, equal κ-cell numbers and depths (no activity gate), and draws that do not depend on the MP list or the
cell order."""
import warnings
from contextlib import nullcontext
from functools import lru_cache

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from scipy import stats

from mpaudit.coexpr import kappa, resample_kappas
from mpaudit.data import Keys, get_counts, log2_cpm10, malignant_mask
from mpaudit.protocol import (KAPPA_ACTIVE, KAPPA_BOOT, KAPPA_FLOOR, KAPPA_MIN_EFFECT, KAPPA_SELF_Z, MAX_KAPPA_CELLS,
                              MIN_CELLS, MIN_DETECT, MIN_GENES, MIN_IOTA, MIN_PAIRS, MIN_RANDOM_SETS, NULL_ALPHA,
                              _boot_rng, _depth_rng, _thin, per_sample_kappa, protocol_shift, protocol_shift_kappa,
                              replicate_bound, tumour_shifts)
from mpaudit.scoring import score_gene_sets
from mpaudit.sim import FIRST_MAGNITUDES, axis_magnitudes, simulate_cohort

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


@lru_cache(maxsize=None)
def _carrier_shift():
    """The small cohort with EMT and IFN, each planted in some tumours only; and which tumours carry them.

    The toy's technical axes are at their first magnitudes (``mpaudit.sim.FIRST_MAGNITUDES``), where
    a program absent from a tumour sits near the background of matched random sets.
    """
    with axis_magnitudes(*FIRST_MAGNITUDES):
        adata, truth = simulate_cohort(**SMALL)
    sets = {k: truth["sets"][k] for k in ("cycle", "hypoxia", "stress", "EMT", "IFN")}
    per_mp, per_tumour, _ = protocol_shift_kappa(adata, KEYS, sets, "fresh", "fixed", n_null_sets=100)
    carriers = {"EMT": set(map(str, truth["tumours_emt"])), "IFN": set(map(str, truth["tumours_ifn"]))}
    return per_mp, per_tumour, carriers


def _levels(paired: dict[str, float], replicates: dict[str, float], metric: str = "intensity",
            base: float = 0.0, scale=lambda b, d: b + d, se: float = 0.0) -> pd.DataFrame:
    """One MP: per tumour a fresh sample at ``base`` and a fixed (paired) or second fresh (replicate) sample.

    For κ every sample gets ``kappa_se`` = ``se`` (0: no self bound in effect)."""
    rows = []
    for t, d in paired.items():
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", metric: base},
                 {"tumour": t, "protocol": "fixed", "mp": "m", metric: scale(base, d)}]
    for t, d in replicates.items():
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", metric: base},
                 {"tumour": t, "protocol": "fresh", "mp": "m", metric: scale(base, d)}]
    lv = pd.DataFrame(rows)
    return lv.assign(kappa_se=se) if metric == "kappa" else lv


T95_3 = stats.t.ppf(0.95, 3)  # 2.353: prediction bound factor for 3 replicate tumours


def _carry(base: float, delta_lin):
    """A linear change of κ carried to an MP's floored reference κ ``base``, as a change in log κ."""
    base = max(base, KAPPA_FLOOR)
    return np.log(np.maximum(base + np.asarray(delta_lin, dtype=float), KAPPA_FLOOR)) - np.log(base)


def test_constants():
    """The rule's fixed constants (stated in report.rule_text; changing one is a new rule version)."""
    assert KAPPA_FLOOR == 0.01 and KAPPA_MIN_EFFECT == pytest.approx(np.log(1.5)) and MIN_PAIRS == 3
    assert MIN_RANDOM_SETS == 20 and NULL_ALPHA == 0.05 and KAPPA_BOOT == 50 and KAPPA_SELF_Z == 1.645
    assert (MAX_KAPPA_CELLS, MIN_CELLS, MIN_DETECT, MIN_IOTA, MIN_GENES) == (2000, 20, 0.05, 0.5, 5)
    assert KAPPA_ACTIVE == 0.05  # informational only


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
    # replicate bound: the lower of the linear bound carried to κ 0.4 and the log-scale bound (t, 3 tumours)
    b_log = T95_3 * np.log(1 / 0.75)
    b_lin = T95_3 * np.sqrt(np.mean(np.array([-0.1, -0.1, 0.4 / 0.75 - 0.4]) ** 2))
    assert sh["null_hi_kappa"].to_numpy() == pytest.approx(max(_carry(0.4, b_lin), b_log))
    assert sh["null_lo_kappa"].to_numpy() == pytest.approx(min(_carry(0.4, -b_lin), -b_log))
    assert _carry(0.4, -b_lin) < -b_log  # here the linear bound decides
    # A (to the floor) is a loss; B (x 0.4) is within the linear bound (v1.2's log-only bound called it a loss)
    assert sh["loss_kappa"].tolist() == [True, False] and sh["gain_kappa"].tolist() == [False, False]
    # several samples of a protocol: kappa is averaged before the log
    two = pd.concat([lv, pd.DataFrame([{"tumour": "B", "protocol": "fixed", "mp": "m", "kappa": 0.4,
                                        "kappa_se": 0.0}])])
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
    # each tumour's floor(0.05 x 42) = 2nd smallest of the 41 changes (not the interpolated 5th percentile, -0.9)
    assert sh["null_lo_intensity"].to_numpy() == pytest.approx(-0.95)
    assert sh["null_hi_intensity"].to_numpy() == pytest.approx(0.95)
    assert sh["loss_intensity"].tolist() == [True, False]
    res = protocol_shift(lv, "fresh", "fixed", null_levels=null).iloc[0]
    assert res["frac_loss_beyond_null_intensity"] == 0.5 and res["n_paired_tumours"] == 2
    # too few random sets in a tumour: that tumour is not tested
    short = pd.concat([null[null["tumour"] == "A"], _random_null("B", np.linspace(-1, 1, 10))])
    assert tumour_shifts(lv, "fresh", "fixed", null_levels=short)["loss_intensity"].isna().tolist() == [False, True]
    # with min_pairs replicate pairs the replicates are the null
    lv3 = _levels({"A": -3.0, "B": -0.5}, {"E": 1.0, "F": -1.0, "G": 0.2})
    assert protocol_shift(lv3, "fresh", "fixed", null_levels=null).iloc[0]["null_type"] == "replicates"


def test_random_set_bound_is_an_order_statistic():
    """The κ random-set bound is the floor(0.05 (n + 1))-th smallest linear change of the tumour's n random sets,
    carried to the MP's κ level: the 5th of 100, the smallest of 20; with 19 the tumour is not tested, nor when
    the k-th smallest is not below the k-th largest."""
    rng = np.random.default_rng(5)
    lv = _levels({t: 0.5 for t in "ABCDE"}, {}, metric="kappa", base=0.6, scale=lambda b, f: b * f)
    rows = []
    for t, n in {"A": 100, "B": 20, "C": 19, "D": 20, "E": 21}.items():
        steps = rng.permutation(n) * 0.002 - 0.05  # the linear changes -0.05, -0.048, ..., shuffled
        if t == "D":
            steps[:] = 0.0  # every set changed alike
        for j, d in enumerate(steps):
            rows += [{"mp": "m", "tumour": t, "protocol": "fresh", "set": j, "kappa": 0.1},
                     {"mp": "m", "tumour": t, "protocol": "fixed", "set": j, "kappa": 0.1 + d}]
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=pd.DataFrame(rows)).set_index("tumour")
    assert sh.loc["A", "null_lo_kappa"] == pytest.approx(_carry(0.6, -0.05 + 4 * 0.002))  # 5th smallest
    assert sh.loc["A", "null_hi_kappa"] == pytest.approx(_carry(0.6, -0.05 + 95 * 0.002))  # 5th largest
    assert sh.loc["B", "null_lo_kappa"] == pytest.approx(_carry(0.6, -0.05))  # smallest of 20
    assert sh.loc["E", "null_lo_kappa"] == pytest.approx(_carry(0.6, -0.05))  # k = floor(1.1) = 1
    assert sh[["null_lo_kappa", "null_hi_kappa"]].loc[["C", "D"]].isna().all().all()
    assert sh["loss_kappa"].isna().tolist() == [False, False, True, True, False]
    assert sh.loc[["A", "B", "E"], "loss_kappa"].all()  # x 0.5 is beyond every bound here


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
    assert np.isnan(res["mean_delta_log_kappa"])  # taken over the tested tumours: none
    assert np.isnan(res["frac_loss_beyond_null_kappa"])


def _kappa_null(rng, kappa_ref: float, kappa_test: float, sd: float, mp_ref: float, mp_test: float,
                tumours: str = "ABCD", n_sets: int = 100) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One MP at ``mp_ref -> mp_test`` (+- 0.003, its ``kappa_se``) and random sets at ``kappa_ref -> kappa_test``
    (+- ``sd``) per tumour."""
    lv, nl = [], []
    for t in tumours:
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": mp_ref + rng.normal(0, 0.003), "kappa_se": 0.003},
               {"mp": "m", "tumour": t, "protocol": "fixed", "kappa": mp_test + rng.normal(0, 0.003), "kappa_se": 0.003}]
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
    fifth = (lin["fixed"] - lin["fresh"]).groupby(level="tumour").apply(lambda x: np.sort(x.to_numpy())[4])
    base = sh.set_index("tumour")["kappa_reference"]
    assert sh.set_index("tumour")["null_lo_kappa"].to_numpy() == pytest.approx(np.log((base + fifth) / base).to_numpy())
    assert (sh["null_lo_kappa"].abs() < 0.05).all()  # the log-scale null had its q05 between +0.05 and +0.3
    # the self bound: -1.645 x sqrt(SE_ref^2 + SE_test^2), carried to the MP's κ
    self_lo = np.log((base - KAPPA_SELF_Z * np.sqrt(2) * 0.003) / base).to_numpy()
    assert sh.set_index("tumour")["self_lo_kappa"].to_numpy() == pytest.approx(self_lo)
    assert sh["kappa_se_reference"].to_numpy() == pytest.approx(0.003)
    # random sets near the floor spread widely on the log scale and hid a twofold loss; not any more
    lv, nl = _kappa_null(rng, 0.012, 0.012, 0.01, 0.60, 0.30)
    res = protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert res["n_paired_tumours"] == 4 and res["frac_loss_beyond_null_kappa"] == 1.0
    # an MP at the floor in both protocols is tested like any other (no activity gate) and is not a loss
    lv, nl = _kappa_null(rng, 0.012, 0.012, 0.01, 0.005, 0.001)
    res = protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert res["n_paired_tumours"] == 4 and res["frac_loss_beyond_null_kappa"] == 0


def test_null_type_per_statistic_ignores_nan_kappa_pairs():
    """κ NaN in the replicate tumours: κ falls back to random sets, intensity keeps its replicates."""
    rng = np.random.default_rng(0)
    lv, nl = [], []
    for t in "ABCD":
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": 0.6, "kappa_se": 0.0, "intensity": 0.0},
               {"mp": "m", "tumour": t, "protocol": "fixed", "kappa": 0.05, "kappa_se": 0.0, "intensity": 0.0}]
        for j in range(100):
            for p in ("fresh", "fixed"):
                nl.append({"mp": "m", "tumour": t, "protocol": p, "set": j, "kappa": 0.05 + rng.normal(0, 0.02),
                           "intensity": rng.normal(0, 0.1)})
    for t in "EFG":  # panel too small in the replicate tumours: κ NaN, intensity finite
        lv += [{"mp": "m", "tumour": t, "protocol": "fresh", "kappa": np.nan, "kappa_se": np.nan, "intensity": 0.0},
               {"mp": "m", "tumour": t, "protocol": "fresh", "kappa": np.nan, "kappa_se": np.nan, "intensity": 0.1}]
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
    rep = levels["tumour"].isin(["T05", "T06", "T07"]) & (levels["mp"] == "stress")
    assert levels.loc[rep, "kappa"].isna().all() and levels.loc[rep, "kappa_se"].isna().all()
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


def _active_levels(paired: dict[str, tuple[float, float]], replicates: dict[str, tuple[float, float]],
                   se: float = 0.0) -> pd.DataFrame:
    """One MP with κ (``kappa_se`` = ``se``) and intensity: paired tumours at ``(reference, test)`` κ; replicates
    with two fresh pieces."""
    rows = []
    for t, (k_ref, k_test) in paired.items():
        rows += [{"sample": f"{t}_fresh", "tumour": t, "protocol": "fresh", "mp": "m", "kappa": k_ref,
                  "intensity": 0.0},
                 {"sample": f"{t}_fixed", "tumour": t, "protocol": "fixed", "mp": "m", "kappa": k_test,
                  "intensity": -3.0}]
    for t, (k1, k2) in replicates.items():
        rows += [{"sample": f"{t}_R1", "tumour": t, "protocol": "fresh", "mp": "m", "kappa": k1, "intensity": 0.0},
                 {"sample": f"{t}_R2", "tumour": t, "protocol": "fresh", "mp": "m", "kappa": k2, "intensity": 0.1}]
    return pd.DataFrame(rows).assign(kappa_se=se)


def _sample_null(lv: pd.DataFrame, kappa_by_sample: dict[str, float], n_sets: int = 40) -> pd.DataFrame:
    """Random sets in every sample of ``lv``: κ about ``kappa_by_sample`` (default 0.01), 0.001 apart."""
    rows = []
    for _, r in lv.iterrows():
        k = kappa_by_sample.get(r["sample"], 0.01)
        rows += [{"mp": "m", "tumour": r["tumour"], "sample": r["sample"], "protocol": r["protocol"], "set": j,
                  "kappa": k + 0.001 * (j - n_sets / 2), "intensity": 0.0} for j in range(n_sets)]
    return pd.DataFrame(rows)


def test_every_paired_tumour_with_defined_change_and_null_is_tested():
    """No activity gate: a tumour where the MP sits below the activity threshold in both pieces (D) is tested;
    kappa_active only fills the informational columns and never changes n_paired_tumours."""
    lv = _active_levels({"A": (0.5, 0.1), "B": (0.5, 0.1), "C": (0.5, 0.1), "D": (0.03, 0.02), "E": (0.02, 0.3)},
                        {"F": (0.5, 0.45), "G": (0.5, 0.55), "H": (0.5, 0.5)})
    sh = tumour_shifts(lv, "fresh", "fixed").set_index("tumour")
    assert sh["kappa_active_reference"].tolist() == [True, True, True, False, False]
    assert sh["kappa_active_test"].tolist() == [True, True, True, False, True]
    assert sh["kappa_active"].tolist() == [True, True, True, False, True]  # D: below 0.05 in both pieces
    assert sh[["kappa_reference_q95_random", "kappa_test_q95_random"]].isna().all().all()
    assert (sh[["kappa_active_threshold_reference", "kappa_active_threshold_test"]] == KAPPA_ACTIVE).all().all()
    # D drops x 2/3 but within the replicate bound carried to its κ of 0.03; E, present only under test, gains
    assert sh["loss_kappa"].tolist() == [True, True, True, False, False]
    assert sh["gain_kappa"].tolist() == [False, False, False, False, True]
    assert sh["loss_intensity"].notna().all()
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["n_paired_tumours"] == 5 and res["frac_loss_beyond_null_kappa"] == 0.6
    assert res["mean_delta_log_kappa"] == pytest.approx((3 * np.log(0.2) + np.log(2 / 3) + np.log(15)) / 5)
    for kappa_active in (0.0, 0.6):
        other = protocol_shift(lv, "fresh", "fixed", kappa_active=kappa_active).iloc[0]
        assert other["n_paired_tumours"] == 5 and other["frac_loss_beyond_null_kappa"] == 0.6
    # with random sets each piece's activity threshold is the larger of the floor and the q95 of the random
    # sets' κ under its protocol (mean over the tumour's samples): informational, B is still tested
    nl = _sample_null(lv, {"B_fresh": 0.55, "B_fixed": 0.15})
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl).set_index("tumour")
    q = 0.001 * (np.quantile(np.arange(40), 0.95) - 20)
    assert sh.loc["B", "kappa_reference_q95_random"] == pytest.approx(0.55 + q)
    assert sh.loc["B", "kappa_active_threshold_reference"] == pytest.approx(0.55 + q)
    assert sh.loc["B", "kappa_test_q95_random"] == pytest.approx(0.15 + q)
    assert sh.loc["B", "kappa_active_threshold_test"] == pytest.approx(0.15 + q)
    assert sh.loc["A", "kappa_active_threshold_reference"] == KAPPA_ACTIVE  # random-set q95 0.029 is below the floor
    assert sh["kappa_active"].tolist() == [True, False, True, False, True]
    assert protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]["n_paired_tumours"] == 5


def test_replicate_null_uses_every_pair_with_finite_kappa():
    """Every same-protocol pair with κ finite in both pieces enters the κ null, whatever the MP's activity; a
    pair with κ NaN in one piece leaves the κ null only (intensity keeps it)."""
    lv = _active_levels({"A": (0.5, 0.1), "B": (0.5, 0.1)},
                        {"E": (0.5, 0.4), "F": (0.5, 0.6), "G": (0.5, 0.02), "H": (0.03, 0.02), "I": (0.5, np.nan)})
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["n_null_pairs"] == 4 and res["n_null_tumours"] == 4 and res["null_type"] == "replicates"
    sh = tumour_shifts(lv, "fresh", "fixed")
    d_log = np.log([0.4 / 0.5, 0.6 / 0.5, 0.02 / 0.5, 0.02 / 0.03])
    d_lin = np.array([-0.1, 0.1, -0.48, -0.01])
    t4 = stats.t.ppf(0.95, 4)
    b_log, b_lin = t4 * np.sqrt(np.mean(d_log**2)), t4 * np.sqrt(np.mean(d_lin**2))
    assert sh["null_lo_kappa"].to_numpy() == pytest.approx(min(_carry(0.5, -b_lin), -b_log))
    assert sh["null_hi_kappa"].to_numpy() == pytest.approx(max(_carry(0.5, b_lin), b_log))
    assert sh["null_hi_intensity"].to_numpy() == pytest.approx(stats.t.ppf(0.95, 5) * 0.1)  # intensity keeps I
    # random sets do not change which pairs enter (there is no activity gate)
    both = _sample_null(lv, {"F_R1": 0.65, "F_R2": 0.65})
    assert protocol_shift(lv, "fresh", "fixed", null_levels=both, min_pairs=1).iloc[0]["n_null_pairs"] == 4


def test_replicate_kappa_bound_is_the_lower_of_linear_and_log():
    """The replicate κ bound is min(C(-t s_lin), -t s_log): the linear bound decides for an MP whose κ is near the
    replicate pieces' level, the log bound for one far above it."""
    lv = _active_levels({"A": (0.5, 0.4), "B": (0.9, 0.8)},
                        {"E": (0.5, 0.6), "F": (0.5, 0.4), "G": (0.5, 0.6)})
    sh = tumour_shifts(lv, "fresh", "fixed").set_index("tumour")
    d_lin, d_log = np.array([0.1, -0.1, 0.1]), np.log([1.2, 0.8, 1.2])
    b_lin, b_log = T95_3 * np.sqrt(np.mean(d_lin**2)), T95_3 * np.sqrt(np.mean(d_log**2))
    assert _carry(0.5, -b_lin) < -b_log < _carry(0.9, -b_lin)
    assert sh.loc["A", "null_lo_kappa"] == pytest.approx(_carry(0.5, -b_lin))
    assert sh.loc["B", "null_lo_kappa"] == pytest.approx(-b_log)


def _synthetic(rng, law, n_mp: int, n_paired: int, n_rep: int = 0, n_sets: int = 0, rep_law=None):
    """κ levels (and random-set levels with ``n_sets``) with no protocol effect: per tumour two pieces, fresh and
    fixed (paired ``P*``) or both fresh (replicate ``R*``); the MP and each of its random sets draw every piece's
    κ independently from ``law`` (``rep_law`` for the MP in replicate tumours)."""
    def draw(lw, shape):
        z = rng.standard_normal(shape)
        return lw[1] + lw[2] * z if lw[0] == "normal" else np.exp(np.log(lw[1]) + lw[2] * z)

    tum = [f"P{i}" for i in range(n_paired)] + [f"R{i}" for i in range(n_rep)]
    T = len(tum)
    k = np.concatenate([draw(law, (n_mp, n_paired, 2)), draw(rep_law or law, (n_mp, n_rep, 2))], axis=1)
    prot = np.array([["fresh", "fixed" if t.startswith("P") else "fresh"] for t in tum])
    lv = pd.DataFrame({"mp": np.repeat(np.arange(n_mp), 2 * T), "tumour": np.tile(np.repeat(tum, 2), n_mp),
                       "protocol": np.tile(prot.ravel(), n_mp), "kappa": k.ravel()})
    if not n_sets:
        return lv, None
    ks = draw(law, (n_mp, T, n_sets, 2))
    nl = pd.DataFrame({"mp": np.repeat(np.arange(n_mp), T * n_sets * 2),
                       "tumour": np.tile(np.repeat(tum, n_sets * 2), n_mp),
                       "protocol": np.tile(np.tile(prot[:, None, :], (1, n_sets, 1)).ravel(), n_mp),
                       "set": np.tile(np.tile(np.repeat(np.arange(n_sets), 2), T), n_mp), "kappa": ks.ravel()})
    return lv, nl


@pytest.mark.parametrize("law", [("lognormal", 0.06, 0.35), ("normal", 0.10, 0.03)])
@pytest.mark.parametrize("path", ["random_sets", "replicates"])
def test_no_gate_calibration(law, path):
    """No protocol effect and the MP drawn like its random sets (or like its replicate pieces): at most about 5% of
    the tested tumours are losses, and few MPs meet criterion (c) on 4 paired tumours.

    v1.2's activity gate tested a tumour only where the MP's κ was above the random sets' q95 in one piece:
    0.24-0.31 of the tested tumours were then losses here on the random-set path (scripts/calibration). The
    synthetic levels have no cells, hence no kappa_se: the self bound is off (self_bound=False), so this tests
    the null bounds alone (the self bound only removes losses)."""
    rng = np.random.default_rng(11)
    lv, nl = _synthetic(rng, law, 300, 4, n_rep=0 if path == "random_sets" else 6,
                        n_sets=100 if path == "random_sets" else 0)
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl, self_bound=False)
    sh = sh[sh["tumour"].str.startswith("P")]
    assert sh["loss_kappa"].notna().all() and (sh["null_type_kappa"] == path).all()  # every tumour tested
    loss = sh["loss_kappa"].astype(bool)
    assert loss.mean() <= 0.065
    per = protocol_shift(lv, "fresh", "fixed", null_levels=nl, self_bound=False)
    crit_c = (per["n_paired_tumours"] >= 4) & (per["frac_loss_beyond_null_kappa"] >= 0.5)
    assert crit_c.mean() <= 0.03


@pytest.mark.parametrize("sdlog", [0.25, 0.35])
def test_replicates_at_the_floor_and_a_heterogeneous_paired_mp(sdlog):
    """Replicate tumours at the floor, paired tumours carrying an MP whose κ varies between pieces (lognormal,
    ``sdlog``), no protocol effect: the linear bound from near-zero replicate changes alone would call many of
    these proportional changes losses; with the lower of the linear and the log bound at most 5% are."""
    rng = np.random.default_rng(int(sdlog * 100))
    lv, _ = _synthetic(rng, ("lognormal", 0.3, sdlog), 500, 4, n_rep=6, rep_law=("normal", 0.02, 0.015))
    sh = tumour_shifts(lv, "fresh", "fixed", self_bound=False)
    sh = sh[sh["tumour"].str.startswith("P")]
    assert sh["loss_kappa"].notna().all()
    assert sh["loss_kappa"].astype(bool).mean() <= 0.05
    # the linear bound alone (carried to each tumour's κ, with the margin)
    r = lv[lv["tumour"].str.startswith("R")]
    d = r.groupby(["mp", "tumour"])["kappa"].diff().dropna()
    b_lin = stats.t.ppf(0.95, 6) * np.sqrt((d**2).groupby(r.loc[d.index, "mp"]).mean())
    base = np.maximum(sh["kappa_reference"].to_numpy(), KAPPA_FLOOR)
    lin_lo = np.log(np.maximum(base - b_lin.reindex(sh["mp"]).to_numpy(), KAPPA_FLOOR)) - np.log(base)
    assert (sh["delta_log_kappa"] < np.minimum(lin_lo, -KAPPA_MIN_EFFECT)).mean() > 0.1


def test_a_patchy_program_widens_the_replicate_null_instead_of_counting():
    """No protocol effect; the program is absent from some pieces (regional patchiness): from one piece of
    replicate tumour R2 and from the test pieces of paired tumours P3 and P4. R2's pair is in the κ null, as
    P3 and P4 are tested, so their drops are within the null."""
    lv = _active_levels({"P1": (0.50, 0.47), "P2": (0.48, 0.52), "P3": (0.50, 0.03), "P4": (0.52, 0.03)},
                        {"R1": (0.50, 0.45), "R2": (0.50, 0.03), "R3": (0.47, 0.52), "R4": (0.50, 0.55)})
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["null_type"] == "replicates" and res["n_null_pairs"] == 4 and res["n_paired_tumours"] == 4
    assert res["frac_loss_beyond_null_kappa"] == 0


def test_patchy_program_on_three_replicate_tumours_is_a_documented_gap():
    """A program absent from 20% of the pieces (present: lognormal(0.4, 0.15); absent: N(0.02, 0.015)), paired and
    replicate tumours exchangeable, no protocol effect, 3 replicate tumours (the registered minimum): per tumour the
    t bound holds about 5%, but every paired tumour shares one bound estimated from a few heavy-tailed changes, so
    losses cluster within the MP and P(c) on 4 paired tumours is about 0.04, above the at most 0.025 of the other
    synthetic laws
    (scripts/calibration ``C2-info``; docs/05_pilot_P1.md, Caveats). This pins the documented size of the gap."""
    rng = np.random.default_rng(7)
    n_mp, tum = 1500, [f"P{i}" for i in range(4)] + [f"R{i}" for i in range(3)]
    present = rng.random((n_mp, len(tum), 2)) >= 0.2
    k = np.where(present, np.exp(np.log(0.4) + 0.15 * rng.standard_normal(present.shape)),
                 0.02 + 0.015 * rng.standard_normal(present.shape))
    prot = np.array([["fresh", "fixed" if t.startswith("P") else "fresh"] for t in tum])
    lv = pd.DataFrame({"mp": np.repeat(np.arange(n_mp), 2 * len(tum)), "tumour": np.tile(np.repeat(tum, 2), n_mp),
                       "protocol": np.tile(prot.ravel(), n_mp), "kappa": k.ravel()})
    sh = tumour_shifts(lv, "fresh", "fixed", self_bound=False)
    assert sh.loc[sh["tumour"].str.startswith("P"), "loss_kappa"].astype(bool).mean() <= 0.065
    per = protocol_shift(lv, "fresh", "fixed", self_bound=False)
    crit_c = ((per["n_paired_tumours"] >= 4) & (per["frac_loss_beyond_null_kappa"] >= 0.5)).mean()
    assert 0.025 < crit_c < 0.07


def test_kappa_tested_in_tumours_without_the_program():
    """EMT and IFN are planted in some tumours; elsewhere their κ is background. There is no activity gate, so
    those tumours are tested like the carriers, and none is a loss (v1.2 left them untested: its gate selected
    tumours on the κ it then tested)."""
    per_mp, pt, carriers = _carrier_shift()
    per_mp = per_mp.set_index("mp")
    for mp in ("EMT", "IFN"):
        t = pt[pt["mp"] == mp].set_index("tumour")
        carrier = t.index.isin(sorted(carriers[mp]))
        assert carrier.any() and (~carrier).any() and not carrier[0]  # T01 carries neither
        assert t["loss_kappa"].notna().all() and not t["loss_kappa"].any()
        # informational activity: every carrier is active; tumours without the program sit at background (one
        # piece can exceed its random sets' q95 by chance, as in about 5% of pieces)
        assert t.loc[carrier, "kappa_active"].all()
        k = t[["kappa_reference", "kappa_test"]]
        assert k[~carrier].max().max() < 0.5 * k[carrier].min().min()
        assert per_mp.loc[mp, "n_paired_tumours"] == 4 and per_mp.loc[mp, "frac_loss_beyond_null_kappa"] == 0
        for side in ("reference", "test"):
            threshold = np.fmax(KAPPA_ACTIVE, t[f"kappa_{side}_q95_random"].to_numpy())
            assert t[f"kappa_active_threshold_{side}"].to_numpy() == pytest.approx(threshold)
        assert np.isfinite(t[["self_lo_kappa", "kappa_se_reference", "kappa_se_test"]].to_numpy()).all()
    # replicate tumours: every pair with a finite κ enters the null, carriers or not
    assert (per_mp["n_null_tumours"] == 3).all() and (per_mp["null_type"] == "replicates").all()
    for mp in ("cycle", "hypoxia", "stress"):
        assert per_mp.loc[mp, "n_paired_tumours"] == 4
    assert per_mp.loc["stress", "frac_loss_beyond_null_kappa"] == 1.0
    assert (per_mp.loc[["cycle", "hypoxia"], "frac_loss_beyond_null_kappa"] == 0).all()


def test_no_random_sets_when_no_mp_passes_the_panel():
    """An MP whose genes are mostly missing (ι < 0.5 in every tumour) gets no random sets; nothing is tested."""
    adata, sets = _cohort()
    absent = {"mostly_absent": sets["cycle"][:10] + [f"NOT_A_GENE{i}" for i in range(30)]}
    per_mp, per_tumour, levels = protocol_shift_kappa(adata, KEYS, absent, "fresh", "fixed", n_null_sets=40)
    assert levels["kappa"].isna().all() and levels["kappa_se"].isna().all()
    assert per_mp.iloc[0]["n_paired_tumours"] == 0 and per_tumour["loss_kappa"].isna().all()
    assert per_tumour["kappa_reference_q95_random"].isna().all()


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
    paired = base["tumour"].isin(["T01", "T02", "T03", "T04"])
    assert (base["kappa_se"].notna() == (base["kappa"].notna() & paired)).all()
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
    assert list(lv.columns) == ["sample", "tumour", "protocol", "mp", "n_cells", "n_kappa_cells", "median_umi",
                                "kappa_thin", "n_genes", "iota", "kappa", "kappa_se", "intensity", "prevalence"]
    assert set(lv["protocol"]) == {"fresh", "fixed"} and (lv["n_cells"] == SMALL["n_mal"]).all()
    assert (lv["n_kappa_cells"] == SMALL["n_mal"]).all()
    # the self bound's SE only where it is used: samples of paired tumours (T01-T04)
    paired = lv["tumour"].isin(["T01", "T02", "T03", "T04"])
    assert (lv.loc[paired, "kappa_se"] > 0).all() and lv.loc[~paired, "kappa_se"].isna().all()
    capped = per_sample_kappa(adata, KEYS, sets, "fresh", "fixed", max_cells=60)
    assert capped["kappa"].notna().all() and not np.allclose(capped["kappa"], lv["kappa"])
    assert (capped["n_kappa_cells"] == 60).all() and (capped["n_cells"] == SMALL["n_mal"]).all()
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


# --------------------------------------------------------------------------- self bound

def _kappa_counts(adata, tumour: str) -> dict[str, sp.csr_matrix]:
    """Per sample of the tumour (equal pieces, so every malignant cell is a κ cell): its malignant cells' counts in
    barcode order, binomially thinned to the tumour's lowest median total count per cell (stream
    mpaudit.protocol._depth_rng), as mpaudit.protocol prepares them."""
    X, obs, mal = get_counts(adata), adata.obs, malignant_mask(adata, KEYS)
    rows = {s: np.flatnonzero((obs["sample"] == s).to_numpy() & mal)
            for s in obs.loc[(obs["tumour"] == tumour).to_numpy() & mal, "sample"].unique()}
    rows = {s: r[np.argsort(adata.obs_names.to_numpy()[r].astype(str), kind="stable")] for s, r in rows.items()}
    med = {s: np.median(np.asarray(X[r].sum(axis=1)).ravel()) for s, r in rows.items()}
    low = min(med.values())
    return {s: _thin(X[r], low / med[s], _depth_rng(0, s)) if med[s] > low else X[r] for s, r in rows.items()}


def _panel_genes(adata, tumour: str, genes: list[str]) -> np.ndarray:
    """Positions of ``genes`` in the tumour's panel (detected in >= 5% of its κ cells, thinned, under each
    protocol), as mpaudit.protocol builds it."""
    counts = _kappa_counts(adata, tumour)
    protocol = adata.obs.groupby("sample", observed=True)["protocol"].first()
    ok = np.ones(adata.n_vars, bool)
    for p in {protocol[s] for s in counts}:
        Y = sp.vstack([Y for s, Y in counts.items() if protocol[s] == p])
        ok &= np.asarray((Y > 0).sum(axis=0)).ravel() / Y.shape[0] >= 0.05  # as protocol._detected
    pos = {g: i for i, g in enumerate(adata.var_names)}
    return np.array([pos[g] for g in dict.fromkeys(genes) if g in pos and ok[pos[g]]], dtype=int)


def test_kappa_se_is_the_bootstrap_sd_from_its_own_stream():
    """kappa_se of a (sample, MP): the SD (ddof 1) of κ over KAPPA_BOOT resamples, with replacement, of the
    sample's κ cells, drawn from default_rng([seed, crc32(sample + '\\x1f' + MP)]); the same values as κ resample
    by resample. Reproducible, and independent of the other MPs in the call."""
    adata, sets = _cohort()
    a = per_sample_kappa(adata, KEYS, {"cycle": sets["cycle"], "stress": sets["stress"]}, "fresh", "fixed")
    b = per_sample_kappa(adata, KEYS, {k: sets[k] for k in ("stress", "hypoxia", "cycle")}, "fresh", "fixed")
    key = ["sample", "mp"]
    a, b = a.set_index(key), b.set_index(key)
    assert a["kappa_se"].notna().sum() == 16  # 8 samples of the 4 paired tumours x 2 MPs (NaN in replicate tumours)
    np.testing.assert_array_equal(a["kappa_se"], b.loc[a.index, "kappa_se"])
    np.testing.assert_array_equal(a["kappa_se"], per_sample_kappa(
        adata, KEYS, {"cycle": sets["cycle"], "stress": sets["stress"]}, "fresh", "fixed").set_index(key)["kappa_se"])
    other = per_sample_kappa(adata, KEYS, {"cycle": sets["cycle"]}, "fresh", "fixed", seed=1).set_index(key)
    se0, se1 = a.xs("cycle", level="mp")["kappa_se"].dropna(), other.xs("cycle", level="mp")["kappa_se"].dropna()
    assert se0.index.equals(se1.index) and not np.allclose(se0.to_numpy(), se1.to_numpy())
    # by hand for one sample (equal pieces: its κ cells are all its malignant cells, in barcode order, thinned to
    # its tumour's lowest median depth; T02_fixed_R1 is the deeper of T02's two samples, so it is thinned)
    sample, tumour = "T02_fixed_R1", "T02"
    assert a.loc[(sample, "cycle"), "kappa_thin"] < 1
    X = _kappa_counts(adata, tumour)[sample]
    Y, s = X[:, _panel_genes(adata, tumour, sets["cycle"])].toarray(), np.asarray(X.sum(axis=1)).ravel()
    rng = _boot_rng(0, sample, "cycle")
    boot = np.array([rng.integers(0, X.shape[0], X.shape[0]) for _ in range(KAPPA_BOOT)])
    one_by_one = np.array([kappa(Y[r], s[r]) for r in boot])
    np.testing.assert_allclose(resample_kappas(Y, s, boot), one_by_one, rtol=1e-10, atol=1e-12)
    assert a.loc[(sample, "cycle"), "kappa_se"] == pytest.approx(np.std(one_by_one, ddof=1), rel=1e-9)
    assert a.loc[(sample, "cycle"), "kappa"] == pytest.approx(kappa(Y, s), rel=1e-12)


def test_kappa_se_halves_with_four_times_the_cells():
    adata, sets = _cohort()
    gs = {k: sets[k] for k in ("cycle", "hypoxia", "stress")}
    full = per_sample_kappa(adata, KEYS, gs, "fresh", "fixed")
    quarter = per_sample_kappa(adata, KEYS, gs, "fresh", "fixed", max_cells=SMALL["n_mal"] // 4)
    ratio = quarter["kappa_se"] / full["kappa_se"]
    assert 1.6 < ratio.median() < 2.5


def _self_case(se_ref, se_test, n_sets: int = 40):
    """Paired tumour A, κ 0.5 -> 0.3 (Δlog κ -0.51, beyond the 1.5-fold margin), with random sets whose changes
    are tiny, and the samples' ``kappa_se`` as given (a list gives several samples per protocol)."""
    se_ref, se_test = np.atleast_1d(se_ref), np.atleast_1d(se_test)
    lv = pd.DataFrame([{"mp": "m", "tumour": "A", "protocol": p, "sample": f"A_{p}{i}", "kappa": k, "kappa_se": se}
                       for p, k, ses in (("fresh", 0.5, se_ref), ("fixed", 0.3, se_test)) for i, se in enumerate(ses)])
    nl = _sample_null(lv, {s: 0.1 for s in lv["sample"]}, n_sets=n_sets).assign(
        kappa=lambda d: 0.1 + 1e-4 * (d["set"] - n_sets / 2) * np.where(d["protocol"] == "fixed", 1.0, 0.5))
    return lv, nl


def test_self_bound_blocks_a_loss_within_the_mps_own_noise():
    """Beyond the random-set bound and the margin, but within -1.645 SE of the MP's own change: not a loss."""
    lv, nl = _self_case(0.1, 0.1)
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert sh["delta_log_kappa"] < min(sh["null_lo_kappa"], -KAPPA_MIN_EFFECT)
    assert sh["self_lo_kappa"] == pytest.approx(_carry(0.5, -KAPPA_SELF_Z * np.sqrt(0.1**2 + 0.1**2)))
    assert sh["delta_log_kappa"] > sh["self_lo_kappa"] and pd.notna(sh["loss_kappa"]) and not bool(sh["loss_kappa"])
    lv, nl = _self_case(0.01, 0.01)
    assert bool(tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]["loss_kappa"])
    # several samples per protocol: SE of the mean = sqrt(sum SE_s^2) / n
    lv, nl = _self_case([0.1, 0.2], [0.05])
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert sh["kappa_se_reference"] == pytest.approx(np.sqrt(0.1**2 + 0.2**2) / 2)
    assert sh["kappa_se_test"] == pytest.approx(0.05)


def test_self_bound_edge_cases():
    # the self bound never makes an untested tumour tested: 19 random sets, no null bound
    lv, nl = _self_case(0.01, 0.01, n_sets=19)
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert np.isfinite(sh["self_lo_kappa"]) and pd.isna(sh["loss_kappa"]) and np.isnan(sh["null_lo_kappa"])
    # an undefined SE in any sample leaves the tumour untested
    lv, nl = _self_case([0.01, np.nan], [0.01])
    sh = tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]
    assert np.isnan(sh["kappa_se_reference"]) and pd.isna(sh["loss_kappa"])
    assert protocol_shift(lv, "fresh", "fixed", null_levels=nl).iloc[0]["n_paired_tumours"] == 0
    # ... but a sample whose κ is undefined is in neither the mean nor its SE
    lv, nl = _self_case([0.01, 0.01], [0.01])
    lv.loc[1, ["kappa", "kappa_se"]] = np.nan
    assert bool(tumour_shifts(lv, "fresh", "fixed", null_levels=nl).iloc[0]["loss_kappa"])
    # a κ test without kappa_se raises; self_bound=False (not the registered rule) drops the self bound
    lv, nl = _self_case(0.1, 0.1)
    with pytest.raises(ValueError, match="kappa_se"):
        tumour_shifts(lv.drop(columns="kappa_se"), "fresh", "fixed", null_levels=nl)
    with pytest.raises(ValueError, match="kappa_se"):
        protocol_shift(lv.drop(columns="kappa_se"), "fresh", "fixed", null_levels=nl)
    off = tumour_shifts(lv.drop(columns="kappa_se"), "fresh", "fixed", null_levels=nl, self_bound=False).iloc[0]
    assert bool(off["loss_kappa"]) and np.isnan(off["self_lo_kappa"])
    # intensity alone needs no kappa_se
    li = _levels({"A": -3.0}, {"E": 1.0, "F": -1.0, "G": 1.0})
    assert bool(tumour_shifts(li, "fresh", "fixed").iloc[0]["loss_intensity"])


# --------------------------------------------------------------------------- equal κ cells per tumour

def _cut(adata, sample: str, n: int):
    """``adata`` without all but the first ``n`` cells of ``sample``."""
    s = (adata.obs["sample"] == sample).to_numpy()
    keep = ~s | (np.cumsum(s) <= n)
    return adata[keep].copy()


def test_kappa_cells_equal_within_a_tumour():
    """Every sample of a tumour is scored on the tumour's smallest sample's number of cells (at most max_cells),
    each drawn from its own stream: the draw does not depend on the order of the samples."""
    adata, sets = _cohort()
    small = _cut(adata, "T01_fixed_R1", 50)
    gs = {"cycle": sets["cycle"]}
    lv = per_sample_kappa(small, KEYS, gs, "fresh", "fixed").set_index("sample")
    assert lv.loc["T01_fixed_R1", "n_cells"] == 50 and lv.loc["T01_fresh_R1", "n_cells"] == SMALL["n_mal"]
    assert (lv.loc[["T01_fresh_R1", "T01_fixed_R1"], "n_kappa_cells"] == 50).all()
    assert (lv.drop(index=["T01_fresh_R1", "T01_fixed_R1"])["n_kappa_cells"] == SMALL["n_mal"]).all()
    capped = per_sample_kappa(small, KEYS, gs, "fresh", "fixed", max_cells=40)
    assert (capped["n_kappa_cells"] == 40).all()
    # the samples in another order: the same κ cells, so the same κ and SE
    order = np.argsort(-small.obs["tumour"].astype(str).str[1:].astype(int).to_numpy(), kind="stable")
    rev = per_sample_kappa(small[order].copy(), KEYS, gs, "fresh", "fixed").set_index("sample")
    pd.testing.assert_frame_equal(rev.loc[lv.index, ["kappa", "kappa_se"]], lv[["kappa", "kappa_se"]])


def test_equal_pieces_use_every_cell():
    """With equal pieces nothing is subsampled: κ is that of all the sample's malignant cells (thinned to the
    tumour's lowest median depth: T01_fresh_R1 is the shallower of its tumour's samples, T06_fresh_R2 the deeper)."""
    adata, sets = _cohort()
    lv = per_sample_kappa(adata, KEYS, {"stress": sets["stress"]}, "fresh", "fixed").set_index("sample")
    assert lv.loc["T01_fresh_R1", "kappa_thin"] == 1 and lv.loc["T06_fresh_R2", "kappa_thin"] < 1
    for sample in ("T01_fresh_R1", "T06_fresh_R2"):
        X = _kappa_counts(adata, sample[:3])[sample]
        g = _panel_genes(adata, sample[:3], sets["stress"])
        assert lv.loc[sample, "kappa"] == pytest.approx(kappa(X[:, g].toarray(), np.asarray(X.sum(axis=1)).ravel()),
                                                        rel=1e-12)


# --------------------------------------------------------------------------- equal depth; draws independent of order

def _thinned_cohort(sample: str, frac: float, seed: int = 9):
    """The small cohort with one sample's counts binomially thinned to ``frac`` of their depth (every cell)."""
    adata, sets = _cohort()
    a = adata.copy()
    X = get_counts(a).tolil()
    rows = np.flatnonzero((a.obs["sample"] == sample).to_numpy())
    X[rows] = _thin(get_counts(a)[rows], frac, np.random.default_rng(seed))
    a.X = sp.csr_matrix(X).astype(np.float32)
    return a, sets


def test_kappa_counts_thinned_to_the_tumours_lowest_median_depth():
    """κ's moment estimate falls as counts per cell fall, so every sample of a tumour is binomially thinned to the
    tumour's lowest median total count per cell before κ; with the fixed piece of T01 at a quarter of its depth, the
    fresh piece is thinned to match and the unchanged programs keep their κ."""
    adata, sets = _cohort()
    gs = {k: sets[k] for k in ("cycle", "hypoxia")}
    base = per_sample_kappa(adata, KEYS, gs, "fresh", "fixed").set_index(["sample", "mp"])
    shallow, _ = _thinned_cohort("T01_fixed_R1", 0.25)
    lv = per_sample_kappa(shallow, KEYS, gs, "fresh", "fixed").set_index(["sample", "mp"])
    fresh, fixed = lv.xs("T01_fresh_R1", level="sample"), lv.xs("T01_fixed_R1", level="sample")
    assert (fixed["kappa_thin"] == 1).all() and (fresh["kappa_thin"] < 0.3).all()
    assert fresh["kappa_thin"].iloc[0] == pytest.approx(fixed["median_umi"].iloc[0] / fresh["median_umi"].iloc[0])
    assert fresh["median_umi"].iloc[0] == base.loc[("T01_fresh_R1", "cycle"), "median_umi"]  # before thinning
    # the other tumours are untouched
    other = ~lv.index.get_level_values("sample").str.startswith("T01")
    pd.testing.assert_series_equal(lv.loc[other, "kappa"], base.loc[other, "kappa"])
    # Δlog κ of T01 stays within the 1.5-fold margin for the unchanged programs (thinning one side alone: about
    # x 0.7-0.9 on the toy, the reviewer's probe)
    d = np.log(fixed["kappa"] / fresh["kappa"])
    assert (d.abs() < KAPPA_MIN_EFFECT).all(), d
    # thinning needs whole-number counts
    scaled = adata.copy()
    scaled.X = (get_counts(adata) * 1.5).astype(np.float32)
    with pytest.raises(ValueError, match="whole numbers"):
        per_sample_kappa(scaled, KEYS, gs, "fresh", "fixed")


def test_normalised_mixed_or_duplicated_input_is_refused():
    """The whole-number check covers every analysed cell, not only thinned ones: with CPM or normalize_total input
    (equal per-cell totals) no sample is thinned, and with only the shallowest piece normalised it is the one left
    unthinned. Repeated barcodes within a sample would make the κ draws depend on cell order."""
    adata, sets = _cohort()
    gs = {k: sets[k] for k in ("cycle", "hypoxia")}
    X = get_counts(adata)
    tot = np.asarray(X.sum(axis=1)).ravel()
    tot[tot == 0] = 1
    for scale in (1e6, 1e4):  # CPM / normalize_total in float64
        a = adata.copy()
        a.X = sp.csr_matrix(sp.diags(scale / tot) @ X)
        with pytest.raises(ValueError, match="whole numbers"):
            per_sample_kappa(a, KEYS, gs, "fresh", "fixed")
    a = adata.copy()  # fixed pieces at 1,000 per cell (the shallowest), fresh pieces raw
    a.X = sp.csr_matrix(sp.diags(np.where(a.obs["protocol"].to_numpy() == "fixed", 1000 / tot, 1.0)) @ X)
    with pytest.raises(ValueError, match="whole numbers"):
        per_sample_kappa(a, KEYS, gs, "fresh", "fixed")
    a = adata.copy()
    rows = np.flatnonzero(malignant_mask(a, KEYS) & (a.obs["sample"] == a.obs["sample"].iloc[0]).to_numpy())
    names = a.obs_names.to_numpy().copy()
    names[rows[1]] = names[rows[0]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # anndata warns about non-unique obs_names
        a.obs_names = names
    with pytest.raises(ValueError, match="repeated obs_names"):
        per_sample_kappa(a, KEYS, gs, "fresh", "fixed")


def test_draws_do_not_depend_on_the_mp_list_or_the_cell_order():
    """An MP's random sets come from its (tumour, MP) stream, and every per-sample draw (κ cells, thinning,
    bootstrap) from the sample's own stream over its cells in barcode order: the MP's bounds and calls do not
    change when MPs are added, dropped or reordered, nor when the cohort's cells are shuffled."""
    adata, sets = _cohort()
    small = _cut(adata, "T02_fixed_R1", 80)  # an unequal tumour, so that κ cells are drawn
    cols = ["kappa_reference", "kappa_test", "null_lo_kappa", "null_hi_kappa", "self_lo_kappa", "loss_kappa",
            "gain_kappa"]

    def run(a, names):
        pt = protocol_shift_kappa(a, KEYS, {k: sets[k] for k in names}, "fresh", "fixed", min_pairs=4,
                                  n_null_sets=40)[1]
        return pt.set_index(["mp", "tumour"])[cols].sort_index()

    ref = run(small, ("cycle", "hypoxia", "stress"))
    assert (ref["null_lo_kappa"] < 0).all()
    shuffled = small[np.random.default_rng(3).permutation(small.n_obs)].copy()
    for other in (run(small, ("stress", "hypoxia", "cycle")), run(small, ("cycle",)),
                  run(small, ("cycle", "cycle_low_capture", "hypoxia", "stress")), run(shuffled, ("cycle", "stress"))):
        common = other.index.intersection(ref.index, sort=False)
        assert len(common) >= 4
        pd.testing.assert_frame_equal(other.loc[common], ref.loc[common])


# --------------------------------------------------------------------------- toy calibration and power

NULL_KW = dict(n_tumours=10, n_twin=0, n_region=10, n_mal=200, n_tme_per_type=10, estimated_covariates=False)


@lru_cache(maxsize=None)
def _toy_null(seed: int, mag: str = "cal", n_mal: int = 200):
    """The toy null of scripts/calibration/evaluate_c.py: 10 tumours with two fresh pieces (R1, R2) each and no
    protocol effect; its 11 gene sets (cycle, hypoxia, stress, EMT, IFN, CNA_7q, 3 random 40-gene sets of
    unplanted genes, the first private programs of T01 and T06); which tumours carry which set."""
    with axis_magnitudes(*FIRST_MAGNITUDES) if mag == "first" else nullcontext():
        adata, truth = simulate_cohort(seed=seed, **{**NULL_KW, "n_mal": n_mal})
    sets = {k: list(truth["sets"][k]) for k in ("cycle", "hypoxia", "stress", "EMT", "IFN", "CNA_7q")}
    planted = {g for s in truth["sets"].values() for g in s} | {g for ps in truth["private"].values() for p in ps
                                                                for g in p}
    free = np.array([g for g in adata.var_names.astype(str) if g.startswith("G") and g not in planted])
    rng = np.random.default_rng(1000 + seed)
    for j in range(3):
        sets[f"random{j}"] = list(rng.choice(free, 40, replace=False))
    tumours = sorted(adata.obs["tumour"].astype(str).unique())
    for t in (tumours[0], tumours[5]):
        sets[f"private_{t}"] = list(truth["private"][t][0])
    car = {"EMT": set(truth["tumours_emt"]), "IFN": set(truth["tumours_ifn"]), "CNA_7q": set(truth["tumours_cna"])}

    def carrier(mp, t):
        if mp in car:
            return t in car[mp]
        return mp.split("_", 1)[1] == t if mp.startswith("private_") else not mp.startswith("random")

    return adata, sets, tumours, carrier


def _pseudo_paired(adata, paired: list[str], cut_test: int | None = None, thin_test: float | None = None,
                   seed: int = 0):
    """R2 of the ``paired`` tumours relabelled as the test protocol (and cut to ``cut_test`` malignant cells, or with
    every count of those pieces binomially thinned to ``thin_test``, as scripts/calibration/evaluate_c.py does)."""
    a = adata.copy()
    hit = a.obs["tumour"].astype(str).isin(paired).to_numpy() & (a.obs["region"].astype(str) == "R2").to_numpy()
    a.obs["protocol"] = np.where(hit, "fixed", a.obs["protocol"].astype(str))
    if thin_test is not None:
        X = get_counts(a)
        rows = np.flatnonzero(hit)
        sub = X[rows].tocoo()
        sub.data = np.random.default_rng([seed, 77]).binomial(sub.data.astype(np.int64), thin_test).astype(np.float64)
        put = sp.csr_matrix((sub.data, (rows[sub.row], sub.col)), shape=X.shape)
        a.X = sp.csr_matrix(sp.diags((~hit).astype(np.float64)) @ X + put).astype(np.float32)
        a.X.eliminate_zeros()
    if cut_test is not None:
        mal = a.obs["malignant"].to_numpy(bool)
        keep = np.ones(a.n_obs, bool)
        for s in pd.unique(a.obs["sample"][hit]):
            rows = np.flatnonzero((a.obs["sample"] == s).to_numpy() & mal)
            keep[rows[cut_test:]] = False
        a = a[keep].copy()
    return a


@lru_cache(maxsize=None)
def _toy_null_run(seed: int, mag: str, regime: str, n_mal: int = 200, cut_test: int | None = None,
                  thin_test: float | None = None):
    adata, sets, tumours, carrier = _toy_null(seed, mag, n_mal)
    paired = tumours if regime == "A" else tumours[:5]
    per_mp, pt, _ = protocol_shift_kappa(_pseudo_paired(adata, paired, cut_test, thin_test, seed), KEYS, sets, "fresh",
                                         "fixed")
    pt = pt.assign(carrier=[carrier(m, t) for m, t in zip(pt["mp"], pt["tumour"])])
    return per_mp, pt


def _absent_loss_rate(pts: list[pd.DataFrame], carrier: bool = False) -> tuple[float, int]:
    """Losses among the tested (set, tumour) pairs whose program is not planted in the tumour (``carrier``: is)."""
    pt = pd.concat(pts)
    pairs = pt[(pt["carrier"] == carrier) & pt["loss_kappa"].notna()]
    return float(pairs["loss_kappa"].astype(bool).mean()), len(pairs)


def _guard(n: int) -> float:
    """Largest absent-pair loss rate these few-seed tests accept: 5% plus 2 binomial SEs of a rate over ``n`` pairs
    (0.09 for 120 pairs). The calibration itself is scripts/calibration/evaluate_c.py's C3 (<= 0.065 over 12 seeds,
    about 720 pairs); a guard at 0.065 over 2 seeds would fail a calibrated rule about a quarter of the time, and this
    one still fails v1.2 (0.13-0.17)."""
    return 0.05 + 2 * np.sqrt(0.05 * 0.95 / n)


def _crit_c(pt: pd.DataFrame) -> pd.Series:
    loss = pt.dropna(subset=["loss_kappa"]).groupby("mp")["loss_kappa"]
    return (loss.count() >= 4) & (loss.apply(lambda x: x.astype(bool).mean()) >= 0.5)


TOY_SEEDS = (0, 1)  # each toy run takes about 20 s; scripts/calibration runs 12 seeds per case


@pytest.mark.parametrize("mag", ["cal", "first"])
def test_toy_null_calibration(mag):
    """Toy null (no protocol effect; scripts/calibration/evaluate_c.py part b, seeds 0-1 here, 0-11 there): about 5%
    of the tested pairs without the program are losses (``_guard``), carriers are never losses, and no random or
    private set meets criterion (c), on all 10 pseudo-paired tumours or on either half of them. v1.2 called
    0.13-0.17 of such pairs losses, and an unplanted random set met (c)."""
    runs = [_toy_null_run(seed, mag, "A") for seed in TOY_SEEDS]
    rate, n = _absent_loss_rate([pt for _, pt in runs])
    assert n >= 100 and rate <= _guard(n)
    for per_mp, pt in runs:
        assert pt["loss_kappa"].notna().mean() > 0.95  # no gate: (almost) every pair is tested
        assert not pt.loc[pt["carrier"], "loss_kappa"].fillna(False).astype(bool).any()
        tumours = sorted(pt["tumour"].unique())
        for part in (tumours, tumours[:5], tumours[5:]):
            crit = _crit_c(pt[pt["tumour"].isin(part)])
            assert not crit[[m for m in crit.index if m.startswith(("random", "private"))]].any()
        assert per_mp.set_index("mp")["null_type"].eq("random_sets").all()


def test_toy_null_calibration_replicate_regime():
    """The same toy null (seed 0) with T01-T05 pseudo-paired and T06-T10 as replicate tumours: every MP with a
    finite κ in the replicate pieces gets the replicate null."""
    per_mp, pt = _toy_null_run(0, "cal", "B")
    rate, n = _absent_loss_rate([pt])
    assert n >= 25 and rate <= _guard(n)
    assert per_mp["null_type"].eq("replicates").all()
    assert not pt.loc[pt["carrier"], "loss_kappa"].fillna(False).astype(bool).any()


def test_unequal_pieces_toy_null_calibration():
    """Test pieces cut to 50 malignant cells against 200-cell reference pieces: κ falls with the number of cells,
    and without equal κ cells per tumour 0.078 of the absent pairs were losses (judge probe, 6 seeds); with them
    at most about 5%."""
    runs = [_toy_null_run(seed, "cal", "A", cut_test=50) for seed in TOY_SEEDS]
    rate, n = _absent_loss_rate([pt for _, pt in runs])
    assert n >= 100 and rate <= _guard(n)
    # carriers: small pieces make their own κ noisier (the harness reports this as C4-info)
    rate, n = _absent_loss_rate([pt for _, pt in runs], carrier=True)
    assert n >= 100 and rate <= 0.05


def test_few_cells_toy_null_calibration():
    """60 malignant cells per piece (v1.2: 0.20 of the absent pairs were losses in the judge probe)."""
    runs = [_toy_null_run(seed, "cal", "A", n_mal=60) for seed in TOY_SEEDS]
    rate, n = _absent_loss_rate([pt for _, pt in runs])
    assert n >= 100 and rate <= _guard(n)
    rate, n = _absent_loss_rate([pt for _, pt in runs], carrier=True)
    assert n >= 100 and rate <= 0.05


def test_unequal_depth_toy_null_calibration():
    """Test pieces at a quarter of their depth (no other protocol effect): with every sample of a tumour thinned to
    its lowest median depth before κ, the pairs without the program stay near 5%, carriers are rarely losses and no
    random or private set meets criterion (c). Without equal depth the review's probe (6 seeds) gave 0.12 of the
    absent pairs, 0.06 of the carriers and a private program meeting (c)."""
    runs = [_toy_null_run(seed, "cal", "A", thin_test=0.25) for seed in TOY_SEEDS]
    rate, n = _absent_loss_rate([pt for _, pt in runs])
    assert n >= 100 and rate <= _guard(n)
    rate, n = _absent_loss_rate([pt for _, pt in runs], carrier=True)
    assert n >= 100 and rate <= 0.05
    for _, pt in runs:
        crit = _crit_c(pt)
        assert not crit[[m for m in crit.index if m.startswith(("random", "private"))]].any()


def test_default_cohort_stress_lost_in_every_paired_tumour():
    """Default toy (seed 0; the fixed protocol has no dissociation stress): stress loses κ in 4 of 4 paired
    tumours and meets criterion (c)."""
    adata, truth = simulate_cohort(seed=0)
    per_mp, pt, _ = protocol_shift_kappa(adata, KEYS, {"stress": truth["sets"]["stress"]}, "fresh", "fixed")
    res = per_mp.iloc[0]
    assert res["n_paired_tumours"] == 4 and res["frac_loss_beyond_null_kappa"] == 1.0


def _inject(adata, genes: list[str], frac: float, seed: int):
    """In every test-protocol piece, permute the counts of a random ``frac`` of ``genes`` across the piece's malignant
    cells, each gene on its own (their correlations go, their counts stay); as in scripts/calibration."""
    import zlib
    pos = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    gidx = np.array([pos[g] for g in genes if g in pos])
    X = adata.X.toarray()
    obs = adata.obs
    for s in sorted(obs.loc[obs["protocol"] == "fixed", "sample"].astype(str).unique()):
        rows = np.flatnonzero((obs["sample"].astype(str) == s).to_numpy() & obs["malignant"].to_numpy(bool))
        rng = np.random.default_rng([seed, zlib.crc32(f"cycle|{s}".encode())])
        for g in rng.choice(gidx, int(round(frac * gidx.size)), replace=False):
            X[rows, g] = X[rng.permutation(rows), g]
    adata.X = sp.csr_matrix(X.astype(np.float32))
    return adata


def test_injected_cycle_loss_is_detected():
    """Toy null, every tumour pseudo-paired, a quarter of the cycle genes permuted in each test piece (a κ drop of
    about 45%): detected in >= 90% of the tumours, and cycle meets criterion (c)."""
    adata, sets, tumours, _ = _toy_null(0)
    inj = _inject(_pseudo_paired(adata, tumours), sets["cycle"], 0.25, 0)
    per_mp, pt, _ = protocol_shift_kappa(inj, KEYS, {"cycle": sets["cycle"]}, "fresh", "fixed")
    assert pt["loss_kappa"].astype(bool).mean() >= 0.9
    assert per_mp.iloc[0]["n_paired_tumours"] == 10 and per_mp.iloc[0]["frac_loss_beyond_null_kappa"] >= 0.5
