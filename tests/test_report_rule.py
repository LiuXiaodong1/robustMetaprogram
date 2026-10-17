"""The frozen P1 GO rule (P1-GO v1.1) on synthetic audit tables."""
import numpy as np
import pandas as pd
import pytest

from mpaudit.data import write_gmt
from mpaudit.report import decide, rule_text, validate_rule, write_summary

GENES = [f"G{i}" for i in range(400)]
MPS = {f"MP{i}": GENES[50 * i:50 * (i + 1)] for i in range(6)}


def tables(**over):
    mps = list(MPS)
    t = {
        "decoy_recruitment": pd.DataFrame(
            [{"mp": m, "arm": arm, "n_decoys": 5, "hits": 0, "recruitment": 0.0, "kappa_retained": 0.0}
             for arm in ("technical", "dosage") for m in mps]),
        "tme_covariation": pd.DataFrame({"mp": mps, "partial_r2": 0.1, "slope": 1.0, "null_q95": 0.3,
                                         "null_q99": 0.4, "n_samples": 20, "spill_adjusted": True,
                                         "partial_r2_rho": np.nan, "slope_rho": np.nan, "null_q95_rho": np.nan,
                                         "null_q99_rho": np.nan}),
        "protocol_shift": pd.DataFrame({"mp": mps, "n_paired_tumours": 6, "frac_loss_beyond_null_kappa": 0.0}),
        "positional": pd.DataFrame({"mp": mps, "best_arm": "1p", "best_arm_excluding_known_clusters": "1p",
                                    "q_excluding_known_clusters": 1.0}),
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


def decoy(t, mp, arm, hits, kappa, n=None):
    dec = t["decoy_recruitment"]
    sel = (dec["mp"] == mp) & (dec["arm"] == arm)
    dec.loc[sel, ["hits", "kappa_retained"]] = [hits, kappa]
    if n is not None:
        dec.loc[sel, "n_decoys"] = n


def test_null_cohort_pivots():
    d, flags = decide(MPS, tables())
    assert d["verdict"] == "PIVOT" and d["n_counted"] == 0
    assert all(d["powered"].values())


def test_three_supported_mps_go():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9}),
               protocol_shift=set_rows("frac_loss_beyond_null_kappa", {"MP1": 0.75}))
    decoy(t, "MP2", "technical", 3, 0.8)
    d, flags = decide(MPS, t)
    assert d["verdict"] == "GO" and set(d["counted_mps"]) == {"MP0", "MP1", "MP2"}


def test_exclusions_and_thresholds():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9, "MP1": 0.35}),  # MP1 below q99
               positional=set_rows("q_excluding_known_clusters", {"MP3": 1e-6}))
    decoy(t, "MP2", "technical", 1, 0.9)  # 1 hit
    decoy(t, "MP3", "technical", 5, 0.9)  # dosage MP
    decoy(t, "MP4", "technical", 5, 0.2)  # low kappa
    decoy(t, "MP5", "dosage", 4, 0.9)
    rule = {"expected_tme_shared_sets": {"IFN": MPS["MP0"][:20]}}
    d, flags = decide(MPS, t, rule)
    f = flags.set_index("mp")
    assert f.loc["MP0", "expected_tme_shared"] and not f.loc["MP0", "crit_b_tme"]
    assert not f.loc["MP1", "crit_b_tme"]
    assert not f.loc["MP2", "crit_a_tech_decoy"]
    assert f.loc["MP3", "dosage"] and f.loc["MP3", "crit_a_tech_decoy"] and not f.loc["MP3", "counted"]
    assert f.loc["MP3", "dosage_reason"] == "arm enrichment (1p)"
    assert not f.loc["MP4", "crit_a_tech_decoy"]
    assert f.loc["MP5", "dosage"] and f.loc["MP5", "dosage_reason"] == "dosage decoy arm (dosage)"
    assert d["n_counted"] == 0 and d["verdict"] == "PIVOT"


def test_curated_artefacts_do_not_count_and_negative_slope_ignored():
    t = tables(tme_covariation=set_rows("partial_r2", {m: 0.9 for m in MPS}))
    t["tme_covariation"].loc[t["tme_covariation"]["mp"] == "MP5", "slope"] = -1.0
    labels = pd.Series({"MP0": 1, "MP1": 1, "MP2": 1, "MP3": 0, "not_audited": 1})
    d, flags = decide(MPS, t, curated_artefacts=labels)
    assert set(d["counted_mps"]) == {"MP3", "MP4"} and d["verdict"] == "PIVOT"
    assert d["curated_artefact_mps"] == ["MP0", "MP1", "MP2"]


def test_underpowered_is_incomplete():
    t = tables()
    t["decoy_recruitment"]["n_decoys"] = 3
    t["tme_covariation"]["n_samples"] = 6
    t["protocol_shift"]["n_paired_tumours"] = 2
    d, _ = decide(MPS, t)
    assert d["verdict"] == "INCOMPLETE"
    d, _ = decide(MPS, {})
    assert d["verdict"] == "INCOMPLETE"


def test_threshold_boundaries():
    t = tables()
    decoy(t, "MP0", "technical", 2, 0.5)       # exactly at both thresholds: a hit
    decoy(t, "MP1", "technical", 1, 0.9)       # one hit short
    decoy(t, "MP2", "technical", 2, 0.49)      # kappa just short
    decoy(t, "MP3", "technical", 4, 0.9, n=4)  # four decoys: below min_decoys
    f = decide(MPS, t)[1].set_index("mp")
    assert f["crit_a_tech_decoy"].to_dict() == {"MP0": True, "MP1": False, "MP2": False, "MP3": False,
                                                "MP4": False, "MP5": False}

    t = tables(protocol_shift=set_rows("frac_loss_beyond_null_kappa", {"MP0": 0.5, "MP1": 1.0, "MP2": 0.49}))
    t["protocol_shift"]["n_paired_tumours"] = [4, 3, 4, 4, 4, 4]
    d, flags = decide(MPS, t)
    f = flags.set_index("mp")
    assert f.loc["MP0", "crit_c_protocol"]       # 50% of exactly 4 paired tumours
    assert not f.loc["MP1", "crit_c_protocol"]   # 3 paired tumours
    assert not f.loc["MP2", "crit_c_protocol"]
    assert d["powered"]["protocol"]
    t["protocol_shift"]["n_paired_tumours"] = 3
    assert not decide(MPS, t)[0]["powered"]["protocol"]

    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.4, "MP1": 0.41, "MP2": 0.41}))
    t["tme_covariation"]["n_samples"] = [20, 10, 9, 10, 10, 10]
    f = decide(MPS, t)[1].set_index("mp")
    assert not f.loc["MP0", "crit_b_tme"]        # equal to q99 is not above it
    assert f.loc["MP1", "crit_b_tme"]            # exactly min_a1_samples
    assert not f.loc["MP2", "crit_b_tme"]


def test_expected_tme_shared_mp_still_counts_through_a_and_c():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9, "MP1": 0.9, "MP2": 0.9}),
               protocol_shift=set_rows("frac_loss_beyond_null_kappa", {"MP1": 0.75}))
    decoy(t, "MP0", "technical", 3, 0.8)
    rule = {"expected_tme_shared_sets": {"IFN": MPS["MP0"][:10]},  # exactly min_overlap genes
            "expected_tme_shared": ["MP1", "MP2"]}
    d, flags = decide(MPS, t, rule)
    f = flags.set_index("mp")
    assert d["expected_tme_shared_mps"] == ["MP0", "MP1", "MP2"]
    assert not f.loc[["MP0", "MP1", "MP2"], "crit_b_tme"].any()
    assert f.loc["MP0", "crit_a_tech_decoy"] and f.loc["MP1", "crit_c_protocol"]
    assert d["counted_mps"] == ["MP0", "MP1"]    # MP2 had only (b)
    rule["expected_tme_shared_sets"] = {"IFN": MPS["MP0"][:9]}  # one gene short of min_overlap
    assert "MP0" not in decide(MPS, t, rule)[0]["expected_tme_shared_mps"]


def test_unadjusted_a1_cannot_count_or_power():
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9, "MP1": 0.9, "MP2": 0.9}))
    t["tme_covariation"]["spill_adjusted"] = False
    d, flags = decide(MPS, t)
    assert not flags["crit_b_tme"].any() and d["n_counted"] == 0
    assert not d["powered"]["tme"] and d["a1_unadjusted_mps"] == list(MPS)
    t["decoy_recruitment"]["n_decoys"] = 3
    t["protocol_shift"]["n_paired_tumours"] = 2
    assert decide(MPS, t)[0]["verdict"] == "INCOMPLETE"  # with spill_adjusted True this would be GO
    t["tme_covariation"]["spill_adjusted"] = True
    assert decide(MPS, t)[0]["verdict"] == "GO"
    # one MP unadjusted: it alone loses (b); A1 stays powered
    t["tme_covariation"]["spill_adjusted"] = [False, True, True, True, True, True]
    d, flags = decide(MPS, t)
    assert d["counted_mps"] == ["MP1", "MP2"] and d["powered"]["tme"] and d["a1_unadjusted_mps"] == ["MP0"]
    # a table without the column cannot show the adjustment
    t["tme_covariation"] = t["tme_covariation"].drop(columns="spill_adjusted")
    d, _ = decide(MPS, t)
    assert d["n_counted"] == 0 and not d["powered"]["tme"]


def test_a1_without_a_null_is_not_powered():
    t = tables(tme_covariation=set_rows("partial_r2", {m: 0.9 for m in MPS}))
    t["tme_covariation"]["null_q99"] = np.nan  # e.g. tme.n_null: 0
    t["decoy_recruitment"]["n_decoys"] = 3
    t["protocol_shift"]["n_paired_tumours"] = 2
    d, _ = decide(MPS, t)
    assert not d["powered"]["tme"] and d["verdict"] == "INCOMPLETE"


def test_rho_sensitivity_blocks_b():
    t = tables(tme_covariation=set_rows("partial_r2", {m: 0.9 for m in ["MP0", "MP1", "MP2", "MP3", "MP4"]}))
    tme = t["tme_covariation"]
    tme["null_q95_rho"], tme["null_q99_rho"] = 0.3, 0.35
    rho = {"MP0": (0.3, 1.0),     # falls below its null at the larger rho_hat
           "MP1": (0.8, 1.0),     # holds
           "MP2": (0.8, -0.5),    # slope turns negative
           "MP3": (np.nan, np.nan),  # no re-run: the estimated rho_hat decides
           "MP4": (0.33, 1.0)}    # above null_q95_rho, below null_q99_rho
    for mp, (r2, slope) in rho.items():
        tme.loc[tme["mp"] == mp, ["partial_r2_rho", "slope_rho"]] = [r2, slope]
    d, flags = decide(MPS, t)
    f = flags.set_index("mp")
    assert f["crit_b_tme"].to_dict() == {"MP0": False, "MP1": True, "MP2": False, "MP3": True, "MP4": False,
                                         "MP5": False}
    assert d["rho_sensitive_mps"] == ["MP0", "MP2", "MP4"] and d["a1_rho_sensitivity_run"]
    assert f.loc["MP1", "rho_sensitivity"] == "holds" and f.loc["MP3", "rho_sensitivity"] == ""
    d, flags = decide(MPS, t, {"a1_null": "null_q95"})  # the same quantile of the re-run's null
    assert flags.set_index("mp").loc["MP4", "crit_b_tme"]


def test_dosage_arm_needs_no_minimum_number_of_decoys():
    # spec A: a dosage-arm hit is >= 2 hits and kappa_retained >= 0.5, whatever the number of decoys
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9, "MP1": 0.9, "MP2": 0.9}))
    t["decoy_recruitment"]["n_decoys"] = 3
    decoy(t, "MP0", "dosage", 3, 0.9)
    decoy(t, "MP3", "dosage", 2, 0.5)
    decoy(t, "MP4", "technical", 3, 0.9)  # the technical arm still needs >= 5 decoys
    d, flags = decide(MPS, t)
    f = flags.set_index("mp")
    assert d["dosage_mps"] == ["MP0", "MP3"] and not f.loc["MP4", "crit_a_tech_decoy"]
    assert d["counted_mps"] == ["MP1", "MP2"] and d["verdict"] == "PIVOT"
    assert not d["powered"]["technical_decoys"]


def test_positional_dosage_names_the_arm_of_its_q():
    t = tables(positional=set_rows("q_excluding_known_clusters", {"MP0": 0.01}))
    pos = t["positional"]
    pos.loc[pos["mp"] == "MP0", ["best_arm", "best_arm_excluding_known_clusters"]] = ["6p", "1q"]
    f = decide(MPS, t)[1].set_index("mp")
    assert f.loc["MP0", "positional_arm"] == "6p" and f.loc["MP0", "positional_arm_excluding_known"] == "1q"
    assert f.loc["MP0", "dosage"] and f.loc["MP0", "dosage_reason"] == "arm enrichment (1q)"


def test_decoy_arm_outside_the_rule_is_reported():
    t = tables()
    extra = t["decoy_recruitment"].query("arm == 'dosage'").assign(arm="cna")
    t["decoy_recruitment"] = pd.concat([t["decoy_recruitment"], extra], ignore_index=True)
    decoy(t, "MP0", "cna", 5, 0.9)
    with pytest.warns(UserWarning, match="cna"):
        d, flags = decide(MPS, t)
    assert d["ignored_decoy_arms"] == ["cna"] and "MP0" not in d["dosage_mps"]
    d, _ = decide(MPS, t, {"dosage_arms": ["dosage", "cna"]})
    assert d["dosage_mps"] == ["MP0"] and d["ignored_decoy_arms"] == []


def test_validate_rule(tmp_path):
    with pytest.raises(ValueError, match=r"unknown decision key\(s\) \['min_hit', 'typo'\]"):
        validate_rule({"min_hit": 2, "typo": 1})
    with pytest.raises(ValueError, match=r"P1-GO v1\.0.*recruitment by min_hits"):
        validate_rule({"recruitment": 0.3, "partial_r2": 0.3, "protocol_frac_beyond_null": 0.75})
    r = validate_rule({"technical_arms": "technical", "expected_tme_shared": "MP0", "min_hits": "2",
                       "protocol_frac": "0.5", "min_decoys": 5.0, "expected_tme_shared_sets": None})
    assert r["technical_arms"] == ["technical"] and r["expected_tme_shared"] == ["MP0"]
    assert r["min_hits"] == 2 and isinstance(r["min_hits"], int) and r["protocol_frac"] == 0.5
    assert isinstance(r["min_decoys"], int) and r["expected_tme_shared_sets"] == {}
    for bad in [{"min_hits": 2.5}, {"min_hits": True}, {"protocol_frac": "half"}, {"protocol_frac": 1.5},
                {"min_decoys": 0}, {"a1_null": "q99"}, {"technical_arms": [1]}, {"positional_q": None},
                {"technical_arms": ["clone"], "dosage_arms": ["clone"]},
                {"expected_tme_shared_sets": ["G1", "G2"]},
                {"expected_tme_shared_gmt": str(tmp_path / "missing.gmt")}]:
        with pytest.raises(ValueError):
            validate_rule(bad)
    # the review's scenarios: a YAML scalar name excludes that MP, and a scalar arm name works
    t = tables(tme_covariation=set_rows("partial_r2", {m: 0.9 for m in MPS}))
    d, _ = decide(MPS, t, {"expected_tme_shared": "MP0", "technical_arms": "technical"})
    assert d["expected_tme_shared_mps"] == ["MP0"] and "MP0" not in d["counted_mps"]
    gmt = tmp_path / "tme_shared.gmt"
    write_gmt({"IFN": MPS["MP1"][:12]}, gmt)
    d, _ = decide(MPS, t, {"expected_tme_shared_gmt": str(gmt)})
    assert d["expected_tme_shared_mps"] == ["MP1"]


def test_summary_lists_counted_mps_and_notes_truncated_tables(tmp_path):
    mps = {f"MP{i}": [f"G{i}"] for i in range(70)}
    tme = pd.DataFrame({"mp": list(mps), "partial_r2": 0.1, "slope": 1.0, "null_q95": 0.3, "null_q99": 0.4,
                        "n_samples": 20, "spill_adjusted": True})
    tme.loc[tme["mp"].isin(["MP65", "MP66", "MP67"]), "partial_r2"] = 0.9
    tme.loc[tme["mp"] == "MP3", "spill_adjusted"] = False
    d, flags = decide(mps, {"tme_covariation": tme})
    assert d["verdict"] == "GO"
    runs = pd.DataFrame({"kind": "seed", "recall": np.linspace(0, 1, 80), "precision": 1.0,
                         "mean_jaccard_matched": 0.5})
    write_summary(tmp_path / "summary.md", d, flags, {"instability_runs": runs}, {"decision": {}})
    text = (tmp_path / "summary.md").read_text()
    assert "nuisance support (GO needs >= 3): MP65, MP66, MP67." in text
    assert "| MP69 |" in text                   # the per-MP table is not cut
    assert "(first 60 of 80 rows; all rows in instability_runs.csv)" in text
    assert "Spillover NOT adjusted for 1 MP(s)" in text and "through (b): MP3." in text
    assert rule_text({}) in text


def test_rule_text_follows_the_settings():
    text = rule_text({"min_hits": 3, "min_decoys": 6, "dosage_arms": ["clone"], "a1_null": "null_q95"},
                     kappa_min_effect=float(np.log(2)), rho_sensitivity=2.0)
    for part in ["recruits the MP in >= 3 of >= 6 decoys", "dosage decoy arm (clone) with >= 3 hits",
                 "random-gene-set null_q95", "multiplied by 2", "partial_r2_rho > null_q95_rho",
                 "by more than 0.693 in log kappa (a 2-fold loss)", "whatever the number of decoys"]:
        assert part in text
    assert "re-run is off" in rule_text({}, rho_sensitivity=None)
    assert "beyond the null under" in rule_text({}, kappa_min_effect=0.0)
