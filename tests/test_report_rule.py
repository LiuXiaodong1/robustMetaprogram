"""The frozen P1 GO rule (P1-GO v1.0) on synthetic audit tables."""
import pandas as pd

from mpaudit.report import decide

GENES = [f"G{i}" for i in range(400)]
MPS = {f"MP{i}": GENES[50 * i:50 * (i + 1)] for i in range(6)}


def tables(**over):
    mps = list(MPS)
    t = {
        "decoy_recruitment": pd.DataFrame(
            [{"mp": m, "arm": arm, "n_decoys": 5, "hits": 0, "recruitment": 0.0, "kappa_retained": 0.0}
             for arm in ("technical", "dosage") for m in mps]),
        "tme_covariation": pd.DataFrame({"mp": mps, "partial_r2": 0.1, "slope": 1.0, "null_q95": 0.3,
                                         "null_q99": 0.4, "n_samples": 20}),
        "protocol_shift": pd.DataFrame({"mp": mps, "n_paired_tumours": 6, "frac_loss_beyond_null_kappa": 0.0}),
        "positional": pd.DataFrame({"mp": mps, "best_arm": "1p", "q_excluding_known_clusters": 1.0}),
    }
    for name, fn in over.items():
        fn(t[name])
    return t


def set_rows(col, values, arm=None):
    def fn(df):
        for mp, v in values.items():
            mask = df["mp"] == mp
            if arm is not None:
                mask &= df["arm"] == arm
            df.loc[mask, col] = v
    return fn


def test_null_cohort_pivots():
    d, flags = decide(MPS, tables())
    assert d["verdict"] == "PIVOT" and d["n_counted"] == 0


def test_three_supported_mps_go():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9}),
               protocol_shift=set_rows("frac_loss_beyond_null_kappa", {"MP1": 0.75}))
    dec = t["decoy_recruitment"]
    sel = (dec["mp"] == "MP2") & (dec["arm"] == "technical")
    dec.loc[sel, ["hits", "kappa_retained"]] = [3, 0.8]
    d, flags = decide(MPS, t)
    assert d["verdict"] == "GO" and set(d["counted_mps"]) == {"MP0", "MP1", "MP2"}


def test_exclusions_and_thresholds():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9, "MP1": 0.35}),  # MP1 below q99
               positional=set_rows("q_excluding_known_clusters", {"MP3": 1e-6}))
    dec = t["decoy_recruitment"]
    dec.loc[(dec["mp"] == "MP2") & (dec["arm"] == "technical"), ["hits", "kappa_retained"]] = [1, 0.9]  # 1 hit
    dec.loc[(dec["mp"] == "MP3") & (dec["arm"] == "technical"), ["hits", "kappa_retained"]] = [5, 0.9]  # dosage MP
    dec.loc[(dec["mp"] == "MP4") & (dec["arm"] == "technical"), ["hits", "kappa_retained"]] = [5, 0.2]  # low kappa
    dec.loc[(dec["mp"] == "MP5") & (dec["arm"] == "dosage"), ["hits", "kappa_retained"]] = [4, 0.9]
    rule = {"expected_tme_shared_sets": {"IFN": MPS["MP0"][:20]}}
    d, flags = decide(MPS, t, rule)
    f = flags.set_index("mp")
    assert f.loc["MP0", "expected_tme_shared"] and not f.loc["MP0", "crit_b_tme"]
    assert not f.loc["MP1", "crit_b_tme"]
    assert not f.loc["MP2", "crit_a_tech_decoy"]
    assert f.loc["MP3", "dosage"] and f.loc["MP3", "crit_a_tech_decoy"] and not f.loc["MP3", "counted"]
    assert not f.loc["MP4", "crit_a_tech_decoy"]
    assert f.loc["MP5", "dosage"] and f.loc["MP5", "dosage_reason"] == "dosage decoy arm"
    assert d["n_counted"] == 0 and d["verdict"] == "PIVOT"


def test_curated_artefacts_do_not_count_and_negative_slope_ignored():
    t = tables(tme_covariation=set_rows("partial_r2", {m: 0.9 for m in MPS}))
    t["tme_covariation"].loc[t["tme_covariation"]["mp"] == "MP5", "slope"] = -1.0
    labels = pd.Series({"MP0": 1, "MP1": 1, "MP2": 1})
    d, flags = decide(MPS, t, curated_artefacts=labels)
    assert set(d["counted_mps"]) == {"MP3", "MP4"} and d["verdict"] == "PIVOT"


def test_underpowered_is_incomplete():
    t = tables()
    t["decoy_recruitment"]["n_decoys"] = 3
    t["tme_covariation"]["n_samples"] = 6
    t["protocol_shift"]["n_paired_tumours"] = 2
    d, _ = decide(MPS, t)
    assert d["verdict"] == "INCOMPLETE"
    d, _ = decide(MPS, {})
    assert d["verdict"] == "INCOMPLETE"
