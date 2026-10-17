"""The frozen P1 GO rule (P1-GO v1.2) on synthetic audit tables, its registered settings and the summary."""
import numpy as np
import pandas as pd
import pytest

from mpaudit.data import write_gmt
from mpaudit.pipeline import merge_config
from mpaudit.report import (REGISTERED_DECISION_KEYS, REGISTERED_SETTINGS, RULE_VERSION, decide, rule_deviations,
                            rule_label, rule_text, validate_rule, write_summary)

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
                     kappa_min_effect=float(np.log(2)), rho_sensitivity=2.0, min_pairs=5, kappa_active=0.1,
                     min_jaccard=0.4)
    for part in ["recruits the MP (a decoy MP at Jaccard >= 0.4 to it) in >= 3 of its 5 decoys (>= 6 needed)",
                 "dosage decoy arm (clone) with >= 3 hits", "random-gene-set null_q95", "multiplied by 2",
                 "partial_r2_rho > null_q95_rho", "by more than 0.693 in log kappa (a 2-fold loss)",
                 "whatever the number of decoys", "replicate pairs of >= 5 distinct tumours", "above max(0.1, "]:
        assert part in text
    # neither the decision-block thresholds nor the settings given here are the registered ones: the text says so
    assert text.startswith(f"Rule {RULE_VERSION} (unregistered variant: decision.min_decoys=6, decision.min_hits=3, "
                           "decision.a1_null='null_q95', min_jaccard=0.4, tme.rho_sensitivity=2, "
                           "protocol.kappa_min_effect=0.693147, protocol.kappa_active=0.1, protocol.min_pairs=5):")
    assert "This is NOT the registered rule" in text
    off = rule_text({}, rho_sensitivity=None)
    assert "re-run is off, so the estimated rho_hat alone decides" in off
    assert off.startswith(f"Rule {RULE_VERSION} (unregistered variant: tme.rho_sensitivity=None):")
    no_margin = rule_text({}, kappa_min_effect=0.0)
    assert "below the null bound (no minimum effect)." in no_margin
    assert "NOT the registered rule P1-GO v1.2; deviations: protocol.kappa_min_effect = 0 (registered: 0.405465 " \
           "(log 1.5))." in no_margin
    # settings stated by the caller and passed as deviations as well are listed once, in the caller's order
    devs = [("decoy.n", 3, "5"), ("tme.n_null", 40, "500")]
    assert rule_text({}, n_decoys=3, deviations=devs).startswith(
        f"Rule {RULE_VERSION} (unregistered variant: decoy.n=3, tme.n_null=40):")


def test_rule_text_states_nulls_gate_jaccard_power_and_informational_items():
    text = rule_text({})
    assert text.startswith(f"Rule {RULE_VERSION}: GO if >= 3 retained") and "NOT the registered" not in text
    for part in ["(a decoy MP at Jaccard >= 0.3 to it)",
                 "which also holds with every rho_hat multiplied by 1.5",
                 "in >= 0.5 of >= 4 paired tumours where the MP is active",
                 "and by more than 0.405 in log kappa (a 1.5-fold loss)",
                 "one-sided 95% t prediction bound from same-protocol replicate pairs of >= 3 distinct tumours, "
                 "otherwise matched random gene sets carried to the MP's kappa level on the linear scale",
                 "above max(0.05, the q95 of the matched random gene sets' kappa in the same samples)",
                 "Active: the MP's kappa under the reference or the test protocol (mean over the tumour's samples "
                 "under it) is above max(0.05, the q95",
                 "a replicate pair enters the null only where the MP is active in at least one piece",
                 "in >= 2 of its 5 decoys (>= 5 needed) with kappa_retained >= 0.5",
                 "for MPs not pre-declared as expected TME-shared (named, or overlapping a pre-declared gene set by "
                 ">= 10 genes)",
                 "dosage MPs (arm enrichment q, known positional gene clusters excluded, < 0.05, or a dosage decoy "
                 "arm (dosage)",
                 "Curated artefacts (labelled, or carried onto re-derived MPs at Jaccard >= 0.3) are not counted",
                 "(a) >= 5 decoys in some technical arm",
                 "(b) for some MP, >= 10 samples with TME cells, a finite null and the spillover adjusted",
                 "(c) for some MP, >= 4 paired tumours where it is active and its kappa change and null bound are "
                 "defined",
                 "Instability, positional flags outside the dosage rule, the within-sample kappa ratio "
                 "(kappa_ratio_median, not spillover-adjusted), A4 intensity and prevalence and the trivial typer "
                 "are informational"]:
        assert part in text, part


def test_decide_labels_deviations_from_the_registered_rule():
    d, _ = decide(MPS, tables())
    assert d["rule_version"] == RULE_VERSION == "P1-GO v1.2"
    assert d["registered_rule"] and d["deviations"] == []
    # arm names and pre-declared lists are filled in by the user, not registered
    d, _ = decide(MPS, tables(), {"dosage_arms": ["dosage", "cna"], "expected_tme_shared": ["MP0"],
                                  "unstable_recovery": 0.8})
    assert d["registered_rule"]
    # a decision threshold is found by decide itself; settings outside the block come from the caller
    d, _ = decide(MPS, tables(), {"a1_null": "null_q95"},
                  deviations=[("decision.a1_null", "null_q95", "'null_q99'"), ("tme.rho_sensitivity", None, "1.5")])
    assert d["rule_version"] == ("P1-GO v1.2 (unregistered variant: decision.a1_null='null_q95', "
                                 "tme.rho_sensitivity=None)")
    assert not d["registered_rule"]
    assert d["deviations"] == [{"setting": "decision.a1_null", "value": "null_q95", "registered": "'null_q99'"},
                               {"setting": "tme.rho_sensitivity", "value": None, "registered": "1.5"}]
    assert rule_label(d["deviations"]) == d["rule_version"]


def test_rule_deviations_cover_every_registered_setting():
    assert rule_deviations(merge_config({})) == []
    assert rule_deviations({"decision": {}}) == []   # settings missing from the config count as registered
    deviating = {"min_jaccard": 0.5, "seed": 1, "decoy.enabled": False, "decoy.n": 10, "decoy.method": "permute",
                 "decoy.n_depth_bins": 3, "decoy.min_kappa": 0.05, "decoy.n_null_floor": 0,
                 "tme.spillover_adjust": False, "tme.rho_sensitivity": None, "tme.null_match": "expression",
                 "tme.n_null": 499, "tme.min_cells": 10, "tme.n_markers": 30, "tme.min_marker_fold": 4.0,
                 "protocol.kappa_min_effect": 0.0, "protocol.kappa_active": 0.1, "protocol.min_pairs": 2,
                 "protocol.n_null_sets": 50, "positional.enabled": False}
    assert set(deviating) == {s for s, _ in REGISTERED_SETTINGS}
    for setting, value in deviating.items():
        section, _, key = setting.rpartition(".")
        devs = rule_deviations(merge_config({section: {key: value}} if section else {key: value}))
        assert [d[:2] for d in devs] == [(setting, value)]
    for key in REGISTERED_DECISION_KEYS:
        value = "null_q95" if key == "a1_null" else 0.01 if key in ("positional_q", "protocol_frac",
                                                                    "min_kappa_retained") else 7
        assert [d[:2] for d in rule_deviations(merge_config({"decision": {key: value}}))] == [(f"decision.{key}",
                                                                                             value)]
    # larger values deviate too: more null sets re-draw the random sets, and more replicate tumours send an MP to
    # the narrower random-set null; the log 1.5 margin as written in the example config (7 decimals) is registered
    assert rule_deviations(merge_config({"protocol": {"kappa_min_effect": 0.4054651}})) == []
    larger = merge_config({"tme": {"n_null": 1000}, "protocol": {"min_pairs": 4, "n_null_sets": 200},
                           "decoy": {"n_null_floor": 50}})
    assert [d[0] for d in rule_deviations(larger)] == ["decoy.n_null_floor", "tme.n_null", "protocol.min_pairs",
                                                       "protocol.n_null_sets"]
    shown = dict((s, (v, r)) for s, v, r in rule_deviations(merge_config({"tme": {"n_null": 40},
                                                                         "protocol": {"kappa_min_effect": 0}})))
    assert shown == {"tme.n_null": (40, "500"), "protocol.kappa_min_effect": (0, "0.405465 (log 1.5)")}


def test_summary_of_an_unregistered_variant_without_labels(tmp_path):
    t = tables(tme_covariation=set_rows("partial_r2", {"MP0": 0.9}))
    t["tme_covariation"]["kappa_ratio_median"] = 0.4
    cfg = {"decision": {}, "tme": {"rho_sensitivity": None}, "protocol": {"min_pairs": 1}}
    devs = rule_deviations(cfg)
    assert [d[0] for d in devs] == ["tme.rho_sensitivity", "protocol.min_pairs"]
    d, flags = decide(MPS, t, cfg["decision"], deviations=devs)
    assert not d["curated_labels_given"]
    nmf_fits = {"reference": {"n_fits": 60, "n_capped_fits": 2, "n_failed_fits": 0},
                "decoys": {"n_runs": 10, "n_fits": 600, "n_capped_fits": 5, "n_failed_fits": 1}, "instability": None}
    write_summary(tmp_path / "summary.md", d, flags, {}, cfg, nmf_fits=nmf_fits)
    text = (tmp_path / "summary.md").read_text()
    assert text.startswith("# P1 audit summary\n\n**This verdict is NOT the registered rule P1-GO v1.2: the run is an "
                           "unregistered variant (unregistered_variant: true). Deviations: tme.rho_sensitivity = None "
                           "(registered: 1.5); protocol.min_pairs = 1 (registered: 3).**")
    assert ("(P1-GO v1.2 (unregistered variant: tme.rho_sensitivity=None, protocol.min_pairs=1))" in text
            and "replicate pairs of >= 1 distinct tumours" in text and "re-run is off" in text)
    assert "Curated artefacts: not labelled (typer.labels not given; no MP excluded)." in text
    assert ("NMF fits stopped at max_iter (of all fits) / by the restart limit: reference 2/60, 0; decoy runs (10) "
            "5/600, 1." in text)
    assert "kappa_ratio_median (within-sample TME/malignant kappa ratio) is not spillover-adjusted" in text
    # registered rule with labels: no warning line, labels listed
    d, flags = decide(MPS, t, curated_artefacts=pd.Series({"MP1": 1}))
    assert d["curated_labels_given"] and d["curated_artefact_mps"] == ["MP1"]
    write_summary(tmp_path / "summary.md", d, flags, {}, {"decision": {}})
    text = (tmp_path / "summary.md").read_text()
    assert "NOT the registered rule" not in text and "Curated artefacts: MP1." in text
    assert f"**Verdict: PIVOT** ({RULE_VERSION}):" in text and "NMF fits" not in text


def test_decoy_table_without_a_column_is_refused():
    """A decoy table also carries the dosage-arm exclusion, so a missing column raises instead of leaving it out."""
    t = tables()
    t["decoy_recruitment"] = t["decoy_recruitment"].drop(columns="kappa_retained")
    with pytest.raises(ValueError, match=r"decoy_recruitment table lacks column\(s\) \['kappa_retained'\]"):
        decide(MPS, t)
    t["decoy_recruitment"] = t["decoy_recruitment"].iloc[:0]
    d, _ = decide(MPS, t)  # an empty table leaves (a) unpowered
    assert not d["powered"]["technical_decoys"]
