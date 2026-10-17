"""Axis A4 (mpaudit.protocol): paired Δlog κ, protocol-specific intensity controls and the two nulls."""
from functools import lru_cache

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit.data import Keys, get_counts, log2_cpm10, malignant_mask
from mpaudit.protocol import KAPPA_FLOOR, per_sample_kappa, protocol_shift, protocol_shift_kappa, tumour_shifts
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


def test_loss_is_one_sided_against_replicates():
    lv = _levels({"A": -3.0, "B": 2.0, "C": 0.1, "D": -0.1}, {"E": 1.0, "F": -1.0, "G": 1.0})
    res = protocol_shift(lv, "fresh", "fixed").iloc[0]
    assert res["null_type"] == "replicates" and res["n_null_pairs"] == 3 and res["n_paired_tumours"] == 4
    assert res["frac_loss_beyond_null_intensity"] == 0.25  # only A drops by more than the null |delta| = 1
    assert res["frac_gain_beyond_null_intensity"] == 0.25  # only B
    assert res["mean_delta_intensity"] == pytest.approx(-0.25)
    rev = protocol_shift(lv, "fixed", "fresh").iloc[0]
    assert rev["mean_delta_intensity"] == pytest.approx(0.25)
    assert rev["frac_loss_beyond_null_intensity"] == 0.25 and rev["frac_gain_beyond_null_intensity"] == 0.25


def test_kappa_compared_on_floored_log_scale():
    lv = _levels({"A": 0.002, "B": 0.5}, {"C": 0.75, "D": 0.75, "E": 1 / 0.75}, metric="kappa", base=0.4,
                 scale=lambda b, f: b * f)
    sh = tumour_shifts(lv, "fresh", "fixed").set_index("tumour")
    assert sh.loc["A", "delta_log_kappa"] == pytest.approx(np.log(KAPPA_FLOOR / 0.4))  # 0.0008 floored to 0.01
    assert sh.loc["B", "delta_log_kappa"] == pytest.approx(np.log(0.5))
    assert sh["null_hi_kappa"].to_numpy() == pytest.approx(np.log(1 / 0.75))
    assert sh["loss_kappa"].tolist() == [True, True] and sh["gain_kappa"].tolist() == [False, False]
    # several samples of a protocol: kappa is averaged before the log
    two = pd.concat([lv, pd.DataFrame([{"tumour": "B", "protocol": "fixed", "mp": "m", "kappa": 0.4}])])
    sh2 = tumour_shifts(two, "fresh", "fixed").set_index("tumour")
    assert sh2.loc["B", "delta_log_kappa"] == pytest.approx(np.log(0.3 / 0.4))


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


def test_collapsed_random_null_is_untestable():
    lv = _levels({"A": 0.1, "B": 0.1}, {}, metric="kappa", base=0.5, scale=lambda b, f: b * f)
    null = _random_null("AB", np.zeros(30), metric="kappa", base=0.001)  # every random set under the floor
    res = protocol_shift(lv, "fresh", "fixed", null_levels=null).iloc[0]
    assert res["null_type"] == "random_sets" and res["n_paired_tumours"] == 0
    assert res["mean_delta_log_kappa"] == pytest.approx(np.log(0.1))
    assert np.isnan(res["frac_loss_beyond_null_kappa"])


def test_kappa_flags_stress_not_cycle_with_replicate_null():
    per_mp = _shift(3)[0].set_index("mp")
    assert (per_mp["null_type"] == "replicates").all() and (per_mp["n_null_pairs"] == 3).all()
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
