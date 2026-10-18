import hashlib
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mpaudit import Keys, cli, pipeline
from mpaudit.coexpr import kappa
from mpaudit.compare import jaccard_matrix
from mpaudit.data import get_counts, malignant_mask, write_gmt
from mpaudit.decoy import PER_DECOY_COLUMNS, make_decoy
from mpaudit.pipeline import (check_arm_values, check_names, config_deviations, curated_labels, merge_config, run_p1,
                              validate_config)
from mpaudit.sim import simulate_cohort

EXAMPLE_CONFIG = Path(__file__).resolve().parents[1] / "configs" / "p1_3ca_example.yaml"

KEYS = dict(sample="sample", malignant="malignant", tumour="tumour", study="study", cancer_type="cancer_type",
            cell_type="cell_type", protocol="protocol", clone="clone", doublet="doublet_score")


@pytest.fixture(scope="module")
def small():
    return simulate_cohort(n_tumours=6, n_twin=2, n_region=1, n_mal=90, n_tme_per_type=25, n_genes=700, seed=3)


def _sample_rows(adata, keys, s):
    mal = malignant_mask(adata, keys)
    return np.flatnonzero(mal & (adata.obs["sample"].astype(str).to_numpy() == s))


def test_permute_decoy_keeps_stratum_marginals_and_tme(small):
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", strata=("clone",), n_depth_bins=1, method="permute", seed=1)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    tme = ~malignant_mask(adata, keys)
    assert abs(X[tme] - D[tme]).sum() == 0
    s = adata.obs["sample"].iloc[0]
    rows = _sample_rows(adata, keys, s)
    for c in adata.obs["clone"].iloc[rows].unique():
        r = rows[adata.obs["clone"].iloc[rows].to_numpy() == c]
        assert np.allclose(np.asarray(X[r].sum(axis=0)), np.asarray(D[r].sum(axis=0)))


def test_residual_decoy_destroys_program_coexpression(small):
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", method="residual", seed=2)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    pos = {g: i for i, g in enumerate(adata.var_names)}
    idx = [pos[g] for g in truth["sets"]["cycle"]]
    rows = _sample_rows(adata, keys, adata.obs["sample"].iloc[0])

    def mean_corr(M):  # depth-robust latent co-expression; raw log counts would keep depth effects
        Y = M[rows]
        return kappa(Y[:, idx].toarray(), np.asarray(Y.sum(axis=1)).ravel())

    assert mean_corr(X) > 0.15
    assert abs(mean_corr(D)) < 0.05
    tot_x, tot_d = X[rows].sum(), D[rows].sum()
    assert abs(tot_d - tot_x) / tot_x < 0.15


def _labels(tmp_path, truth):
    """Curated labels keyed by published (truth) set names, and the GMT of those sets (the fixture
    re-derives a cycle MP; labelling cycle 1 checks that a label is inherited)."""
    labels = {name: int(name != "stress") for name in truth["sets"]}
    pd.DataFrame({"mp": list(labels), "label": list(labels.values())}).to_csv(tmp_path / "labels.csv", index=False)
    write_gmt(truth["sets"], tmp_path / "labelled.gmt")
    return labels


def test_run_p1_end_to_end(small, tmp_path):
    adata, truth = small
    labels = _labels(tmp_path, truth)
    no_names = tmp_path / "no_names.gtf"   # a GTF without gene_name attributes gives no positions
    no_names.write_text('chr1\tsrc\tgene\t11\t500\t.\t+\t.\tgene_id "ENSG1";\n')
    cfg = {
        "keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour",
        "gavish": {"ranks": [3, 4], "nrun": 1, "max_iter": 200, "min_genes_per_cell": 0},
        "instability": {"seeds": [1], "n_boot": 1},
        "decoy": {"n": 1, "arms": {"technical": ["doublet_est", "contamination_est"], "dosage": ["clone_est"]}},
        "tme": {"n_boot": 20, "n_null": 40},
        "positional": {"cytoband": truth["cytoband"], "gtf": str(no_names)},
        "typer": {"labels": str(tmp_path / "labels.csv"), "labelled_gmt": str(tmp_path / "labelled.gmt")},
        "protocol": {"reference": "fresh", "test": "fixed", "n_null_sets": 20},
        "decision": {"expected_tme_shared_sets": {"IFN": truth["sets"]["IFN"]}},
        "unregistered_variant": True,  # few decoys and null sets, for speed
    }
    out_dir = tmp_path / "out"
    with pytest.warns(UserWarning) as rec:
        out = run_p1(adata, cfg, out_dir, log=lambda *a: None)
    messages = [str(w.message) for w in rec]
    assert any("gave no gene positions; falling back to adata.var" in m for m in messages)
    assert any(m.startswith("unregistered variant: the verdict will be labelled") for m in messages)
    assert not any("expected-TME-shared" in m or "no typer.labels" in m for m in messages)
    d = out["decision"]
    variant = "P1-GO v1.3 (unregistered variant: decoy.n=1, tme.n_null=40, protocol.n_null_sets=20)"
    assert d["rule_version"] == variant and not d["registered_rule"] and d["curated_labels_given"]
    assert [x["setting"] for x in d["deviations"]] == ["decoy.n", "tme.n_null", "protocol.n_null_sets"]
    # 1 decoy per arm (< 5), 9 samples with TME cells (< 10), 2 paired tumours (< 4): nothing is powered
    assert d["verdict"] == "INCOMPLETE" and not any(d["powered"].values())
    for f in ["summary.md", "mp_flags.csv", "run_config.json", "typer_features.csv", "decoy_recruitment.csv",
              "protocol_shift.csv", "tme_covariation.csv", "positional.csv", "curated_label_mapping.csv",
              "typer_eval.csv"]:
        assert (out_dir / f).exists()
    assert set(out["tables"]["decoy_recruitment"]["arm"]) == {"technical", "dosage"}
    flags = out["tables"]["mp_flags"].set_index("mp")
    for c in ["spill_adjusted", "partial_r2_rho", "null_q99_rho", "rho_sensitivity", "positional_arm_excluding_known",
              "n_null_tumours"]:
        assert c in flags
    assert out["tables"]["positional"]["n_placed"].gt(0).all()  # fell back to adata.var positions
    # re-derived MPs inherit the label of their best-matching labelled set (Jaccard >= 0.3), by genes only
    names = list(truth["sets"])
    J = jaccard_matrix(out["mp_sets"], truth["sets"])
    expected = [mp for i, mp in enumerate(out["mp_sets"]) if J[i].max() >= 0.3 and labels[names[J[i].argmax()]]]
    assert expected and d["curated_artefact_mps"] == expected
    mapping = out["tables"]["curated_label_mapping"]
    assert (mapping["source"] == "gene overlap").any() and not (mapping["source"] == "name").any()
    run_cfg = json.loads((out_dir / "run_config.json").read_text())
    assert "multiplied by 1.5" in run_cfg["rule"] and "a 1.5-fold loss" in run_cfg["rule"]
    assert run_cfg["rule"].startswith(f"Rule {variant}:")
    assert "replicate pairs from >= 3 distinct tumours" in run_cfg["rule"] and "no activity gate" in run_cfg["rule"]
    assert run_cfg["rule_version"] == variant and run_cfg["registered_rule"] is False
    assert run_cfg["deviations"] == [{"setting": "decoy.n", "value": 1, "registered": "5"},
                                     {"setting": "tme.n_null", "value": 40, "registered": "500"},
                                     {"setting": "protocol.n_null_sets", "value": 20, "registered": "100"}]
    assert run_cfg["config"]["protocol"]["kappa_active"] == 0.05 and run_cfg["config"]["unregistered_variant"] is True
    assert {"kappa_active", "kappa_active_reference", "kappa_active_test", "kappa_active_threshold_reference",
            "kappa_active_threshold_test"} <= set(out["tables"]["protocol_tumours"])
    # NMF fit diagnostics: per decoy run, per perturbation, and totals in run_config.json and summary.md
    runs = out["tables"]["decoy_runs"]
    assert list(runs.columns) == PER_DECOY_COLUMNS and len(runs) == 2
    n_samples = len(out["reference"].top) // sum(cfg["gavish"]["ranks"])  # programs per sample: 3 + 4
    assert (runs["n_fits"] == n_samples * 2).all() and runs["n_capped_fits"].between(0, runs["n_fits"]).all()
    inst = out["tables"]["instability_runs"]
    assert {"n_fits", "n_capped_fits", "n_failed_fits"} <= set(inst.columns) and (inst["n_fits"] > 0).all()
    fits = run_cfg["nmf_fits"]
    assert fits["reference"]["n_fits"] == out["reference"].n_fits == n_samples * 2
    assert fits["decoys"] == {"n_runs": 2, **{c: int(runs[c].sum()) for c in ("n_fits", "n_capped_fits",
                                                                               "n_failed_fits")}}
    assert fits["instability"]["n_runs"] == len(inst) and fits["instability"]["n_fits"] == int(inst["n_fits"].sum())
    summary = (out_dir / "summary.md").read_text()
    assert "INCOMPLETE" in summary and summary.startswith("# P1 audit summary\n\n**This verdict is NOT the registered")
    assert f"decoy runs (2) {fits['decoys']['n_capped_fits']}/{fits['decoys']['n_fits']}, " in summary
    assert "pct_long" in adata.obs  # adata.var has gene_length (the GTF placed no gene)
    # the input files behind the verdict, for a check against the registered ones
    assert set(run_cfg["input_files"]) == {"typer.labels", "typer.labelled_gmt", "positional.gtf"}
    digest = hashlib.sha256((tmp_path / "labels.csv").read_bytes()).hexdigest()
    assert run_cfg["input_files"]["typer.labels"] == {"path": str(tmp_path / "labels.csv"), "sha256": digest}
    assert f"typer.labels labels.csv {digest[:16]}" in summary


def _no_nmf(*args, **kwargs):
    raise AssertionError("NMF ran before the configuration was checked")


BASE = {"keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour",
        "decoy": {"arms": {"technical": ["doublet_est"], "dosage": ["clone_est"]}}}


@pytest.fixture(scope="module")
def registered(small, tmp_path_factory):
    """BASE plus what the registered rule needs with mps: rederive and positional on: curated labels with the
    gene sets of the labelled MPs, and a cytoband (the cohort's adata.var has gene positions)."""
    _, truth = small
    d = tmp_path_factory.mktemp("registered")
    _labels(d, truth)
    return {**BASE, "typer": {"labels": str(d / "labels.csv"), "labelled_gmt": str(d / "labelled.gmt")},
            "positional": {"cytoband": truth["cytoband"]}, "protocol": {"reference": "fresh", "test": "fixed"}}


@pytest.mark.parametrize("change, message", [
    ({"decoy": {"strata": ["clone", "doublet_score"], "arms": {}}}, "decoy.strata is no longer accepted"),
    ({"decoy": {"arms": {"tech": ["doublet_est"], "cna": ["clone_est"]}}}, r"\['cna', 'tech'\] are in neither"),
    ({"decoy": {"arms": {"technical": ["doublet_est", "no_such_column"]}}}, "not columns of adata.obs"),
    ({"decision": {"recruitment": 0.3, "partial_r2": 0.3}}, "before P1-GO v1.0; P1-GO v1.3 replaces"),
    ({"decision": {"technical_arms": ["technical"], "dosage_arms": ["technical"]}}, "both technical_arms"),
    ({"decision": {"expected_tme_shared": ["MP17"]}}, "cannot be pre-declared with mps: rederive"),
    ({"decision": {"expected_tme_shared_gmt": "missing.gmt"}}, "is not a file"),
    ({"typer": {"labels": "LABELS"}}, "cannot match re-derived MPs"),
    ({"typer": {"labelled_gmt": "LABELLED"}}, "needs typer.labels"),
    ({"tme": {"rho_sensitivity": 0}}, "rho_sensitivity must be > 0"),
    ({"protocol": {"kappa_min_effect": "large"}}, "kappa_min_effect must be a number"),
    ({"protocol": {"kappa_active": "high"}}, "kappa_active must be a number"),
    ({"protocol": {"kappa_active": True}}, "kappa_active must be a number"),
    ({"protocol": {"kappa_active": -0.1}}, "kappa_active must be >= 0"),
    ({"protocol": {"kappa_active": float("inf")}}, "kappa_active must be >= 0"),
    ({"protocol": {"min_pairs": 0}}, "min_pairs must be >= 1"),
    ({"decoy": {"n_null_floor": 19}, "unregistered_variant": True}, "n_null_floor must be >= 20"),
    # fewer than MIN_RANDOM_SETS random sets would leave every tumour on the random-set null untested
    ({"protocol": {"n_null_sets": 10}, "unregistered_variant": True},
     r"protocol.n_null_sets must be 0 \(no random-set null\) or >= 20, got 10"),
    ({"protocol": {"n_null_sets": -1}}, "protocol.n_null_sets must be >= 0"),
    ({"protocol": {"min_pairs": 3.0}}, "min_pairs must be a whole number"),
    ({"protocol": {"min_pairs": True}}, "min_pairs must be a whole number"),
    ({"protocol": {"min_pairs": "3"}}, "min_pairs must be a whole number"),
    ({"tme": {"spillover_adjust": "false"}}, "tme.spillover_adjust must be true or false"),
    ({"tme": {"null_match": "genes"}}, "tme.null_match must be one of"),
    ({"tme": {"n_null": 500.0}}, "tme.n_null must be a whole number"),
    ({"decoy": {"method": "shuffle"}}, "decoy.method must be one of"),
    ({"decoy": {"enabled": "no"}}, "decoy.enabled must be true or false"),
    ({"positional": {"enabled": 0}}, "positional.enabled must be true or false"),
    ({"tme": {"enabled": "no"}}, "tme.enabled must be true or false"),
    ({"protocol": {"enabled": "false"}}, "protocol.enabled must be true or false"),
    ({"instability": {"enabled": 1}}, "instability.enabled must be true or false"),
    # misspelt keys would run with the default silently (here as registered, and with criterion c unpowered)
    ({"tme": {"rho_sensitvity": None}}, r"unknown config key\(s\).*'tme.rho_sensitvity' \(valid: enabled, min_cells"),
    ({"protocol": {"referance": "fresh", "test": "fixed"}}, r"unknown config key\(s\).*'protocol.referance'"),
    ({"sed": 3}, r"unknown config key\(s\), probably misspelt \(no part of mpaudit reads them\): 'sed'"),
    ({"tme": None}, "config section tme must be a mapping"),
    ({"keys": {**KEYS, "protocol": "Protocol"}}, r"Not a column of adata.obs: keys.protocol='Protocol'"),
    ({"protocol": {"reference": "fresh", "test": "fresh"}}, "protocol.reference and protocol.test must differ"),
    # files the positional axis reads, checked before any NMF (not when the axis runs, after it)
    ({"positional": {"cytoband": "/no/such/cytoBand.txt.gz"}}, "positional.cytoband '/no/such/cytoBand.txt.gz' is not "
                                                               "a file"),
    ({"positional": {"cytoband": ""}}, "positional.cytoband '' is not a file"),
    ({"positional": {"gtf": "/no/such.gtf"}}, "positional.gtf '/no/such.gtf' is not a file"),
    ({"min_jaccard": 1.5}, "min_jaccard must be <= 1"),
    ({"unregistered_variant": "yes"}, "unregistered_variant must be true or false"),
    ({"decoy": {"arms": {"technical": ["doublet_est", "pct_long"]}}, "TINY_VAR": True},
     r"pct_long needs gene lengths: positional.gtf or adata.var\['gene_length'\]"),
    # deviations from the registered rule, refused without unregistered_variant: true
    ({"tme": {"rho_sensitivity": None}}, r"tme.rho_sensitivity = None \(registered: 1.5\)"),
    ({"tme": {"rho_sensitivity": 0.5}}, r"tme.rho_sensitivity = 0.5 \(registered: 1.5\)"),
    ({"protocol": {"kappa_min_effect": 0.0}}, r"protocol.kappa_min_effect = 0.0 \(registered: 0.405465 \(log 1.5\)\)"),
    ({"protocol": {"min_pairs": 2}}, r"protocol.min_pairs = 2 \(registered: 3\)"),
    ({"protocol": {"min_pairs": 4}}, r"protocol.min_pairs = 4 \(registered: 3\)"),
    ({"seed": 7}, r"seed = 7 \(registered: 0\)"),
    ({"seed": 1.5}, "seed must be a whole number"),
    # BASE re-derives MPs without curated labels and has no cytoband: two more deviations, listed after the rule's
    ({"decision": {"min_mps": 2}}, r"deviates from the registered rule P1-GO v1.3 in 3 setting\(s\): decision.min_mps "
                                   r"= 2 \(registered: 3\); typer.labels = None \(registered: curated labels with "
                                   r"typer.labelled_gmt \(needed with mps: rederive\)\); positional.cytoband = None"),
    ({"typer": {"labels": "RLABELS"}}, r"typer.labelled_gmt = None \(registered: a GMT of the labelled MPs"),
    ({"mps": "missing.gmt"}, "mps must be 'rederive' or a GMT file"),
    ({"mps": "PUBLISHED", "typer": {"labels": "LABELS"}, "unregistered_variant": True},  # BASE has no cytoband
     "typer.labels names match none of the 1 audited MPs"),
])
def test_bad_config_fails_before_any_nmf(small, tmp_path, monkeypatch, change, message):
    adata, truth = small
    pd.DataFrame({"mp": ["MP1 Cell Cycle - G2/M", "MP5 Stress"], "label": [0, 1]}).to_csv(tmp_path / "l.csv",
                                                                                         index=False)
    write_gmt({"MP5 Stress": truth["sets"]["stress"]}, tmp_path / "l.gmt")
    write_gmt({"Stress": truth["sets"]["stress"]}, tmp_path / "published.gmt")
    pd.DataFrame({"mp": ["MP1", "MP2"], "label": [0, 1]}).to_csv(tmp_path / "rl.csv", index=False)
    files = {"LABELS": str(tmp_path / "l.csv"), "LABELLED": str(tmp_path / "l.gmt"),
             "PUBLISHED": str(tmp_path / "published.gmt"), "RLABELS": str(tmp_path / "rl.csv")}

    def sub(x):
        return {k: sub(v) for k, v in x.items()} if isinstance(x, dict) else files.get(x, x) if isinstance(x, str) else x

    change = dict(change)
    if change.pop("TINY_VAR", False):  # a cohort without gene lengths in var (and no pct_long from an earlier run)
        adata = adata.copy()
        adata.var = adata.var.drop(columns="gene_length")
        adata.obs = adata.obs.drop(columns="pct_long", errors="ignore")
    cfg = {k: dict(v) if isinstance(v, dict) else v for k, v in BASE.items()}
    for k, v in sub(change).items():
        cfg[k] = {**cfg[k], **v} if isinstance(cfg.get(k), dict) else v
    monkeypatch.setattr(pipeline, "run_gavish", _no_nmf)
    with pytest.raises(ValueError, match=message):
        run_p1(adata, cfg, tmp_path / "out", log=lambda *a: None)
    assert not (tmp_path / "out").exists()


def test_decoy_arm_checks(small, registered):
    adata, _ = small
    cfg = merge_config({**registered, "decoy": {"arms": {"technical": "doublet_est"}}})
    # no dosage arm: the dosage-arm exclusion cannot run, a deviation from the registered rule
    with pytest.raises(ValueError, match=r"decoy.arms = \['technical'\] \(registered: an arm named in "
                                         r"decision.dosage_arms"):
        _quiet(validate_config, cfg, adata, log=lambda *a: None)
    cfg = _quiet(validate_config, merge_config({**cfg, "unregistered_variant": True}), adata, log=lambda *a: None)
    assert cfg["decoy"]["arms"] == {"technical": ["doublet_est"]}  # a scalar column becomes a list
    cfg = merge_config({**registered, "decoy": {"arms": {"technical": ["pct_ribo"], "cna": ["clone_est"]}},
                        "decision": {"dosage_arms": ["cna"]}})
    _quiet(validate_config, cfg, adata, log=lambda *a: None)  # renamed arm listed in the rule; QC columns come later
    off = {**registered, "decoy": {"enabled": False, "arms": {"tech": ["doublet_est"]}}}
    with pytest.raises(ValueError, match="decoy.enabled = False"):  # decoys off drop the dosage-arm exclusion
        _quiet(validate_config, merge_config(off), adata, log=lambda *a: None)
    cfg = _quiet(validate_config, merge_config({**off, "unregistered_variant": True}), adata, log=lambda *a: None)
    assert [d[0] for d in config_deviations(cfg)] == ["decoy.enabled"]  # arm names unused, not a second deviation


def test_protocol_settings_are_normalised(small, registered):
    adata, _ = small
    cfg = _quiet(validate_config, merge_config(registered), adata, log=lambda *a: None)
    assert cfg["protocol"]["kappa_active"] == 0.05 and cfg["protocol"]["min_pairs"] == 3 and cfg["seed"] == 0
    cfg = _quiet(validate_config, merge_config({**registered, "protocol": {"kappa_active": "0.05",
                                                                            "min_pairs": np.int64(3)}}),
                 adata, log=lambda *a: None)
    assert type(cfg["protocol"]["kappa_active"]) is float and type(cfg["protocol"]["min_pairs"]) is int
    assert config_deviations(cfg) == []
    # kappa_active only fills informational columns (no gate in P1-GO v1.3): another value is not a deviation
    cfg = _quiet(validate_config, merge_config({**registered, "protocol": {"kappa_active": 0.2}}), adata,
                 log=lambda *a: None)
    assert config_deviations(cfg) == [] and cfg["decoy"]["n_null_floor"] == 50


def test_informational_keys_naming_no_column_are_warned_about(small, registered):
    """keys.clone, region, doublet and intronic feed no part of the rule (the typer skips a missing column): one that
    names no obs column is warned about, not refused; any other key is refused
    (test_bad_config_fails_before_any_nmf)."""
    adata, _ = small
    cfg = merge_config({**registered, "keys": {**KEYS, "intronic": "intronic_fraction_missing", "doublet": "dbl"}})
    with pytest.warns(UserWarning, match=r"ignored \(informational keys: the rule does not read them\): "
                                         r"keys.doublet='dbl', keys.intronic='intronic_fraction_missing'"):
        validate_config(cfg, adata, log=lambda *a: None)


def _quiet(fn, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fn(*args, **kwargs)


def test_registered_defaults_and_the_example_config_pass(small, registered, tmp_path):
    adata, _ = small
    cfg = _quiet(validate_config, merge_config(registered), adata, log=lambda *a: None)
    assert config_deviations(cfg) == [] and cfg["unregistered_variant"] is False
    yaml = pytest.importorskip("yaml")
    example = yaml.safe_load(EXAMPLE_CONFIG.read_text())
    assert example["unregistered_variant"] is False and example["min_jaccard"] == 0.3 and example["seed"] == 0
    # as shipped, the reference files (cytoband, GTF) are not in the repository...
    with pytest.raises(ValueError, match="positional.cytoband 'ref/cytoBand.hg38.txt.gz' is not a file"):
        _quiet(validate_config, merge_config(example), log=lambda *a: None)  # no cohort: obs not read
    cyto, gtf = tmp_path / "cytoBand.txt", tmp_path / "genes.gtf"
    small[1]["cytoband"].to_csv(cyto, sep="\t", header=False, index=False)
    gtf.write_text('chr1\tsrc\tgene\t11\t500\t.\t+\t.\tgene_id "G1"; gene_name "G1";\n')
    example["positional"] = {**example["positional"], "cytoband": str(cyto), "gtf": str(gtf)}
    # ...and the curated labels are still to be filled in, which mps: rederive needs
    with pytest.raises(ValueError, match=r"in 1 setting\(s\): typer.labels = None"):
        _quiet(validate_config, merge_config(example), log=lambda *a: None)
    example["typer"] = {**example["typer"], **registered["typer"]}
    cfg = _quiet(validate_config, merge_config(example), log=lambda *a: None)
    assert config_deviations(cfg) == []


@pytest.mark.parametrize("change", [
    {"tme": {"rho_sensitivity": None}}, {"tme": {"rho_sensitivity": 1.0}}, {"tme": {"spillover_adjust": False}},
    {"tme": {"null_match": "expression"}}, {"tme": {"n_null": 100}}, {"tme": {"n_null": 1000}},
    {"protocol": {"kappa_min_effect": 0.2}}, {"protocol": {"min_pairs": 1}}, {"decoy": {"n_null_floor": 20}},
    {"protocol": {"min_pairs": 4}}, {"protocol": {"n_null_sets": 20}}, {"protocol": {"n_null_sets": 200}},
    {"min_jaccard": 0.2}, {"seed": 3}, {"gavish": {"seed": 3}}, {"decoy": {"n": 3}}, {"positional": {"enabled": False}},
    {"decision": {"min_hits": 1}}, {"decision": {"a1_null": "null_q95"}},
    # inputs without which an exclusion does not run although its axis is on
    {"typer": {"labels": None, "labelled_gmt": None}}, {"decoy": {"arms": {"technical": ["doublet_est"]}}},
    {"decision": {"dosage_arms": []}}, {"positional": {"cytoband": None}}])
def test_each_deviation_is_refused_by_default_and_labelled_when_allowed(small, registered, change):
    adata, _ = small
    cfg = {**registered, **{k: {**registered.get(k, {}), **v} if isinstance(v, dict) else v
                            for k, v in change.items()}}
    if change == {"decision": {"dosage_arms": []}}:  # the arm 'dosage' is then named nowhere: rename it technical
        cfg["decoy"] = {"arms": {"technical": ["doublet_est", "clone_est"]}}
    (section, value), = change.items()
    setting = f"{section}.{next(iter(value))}" if isinstance(value, dict) else section
    setting = {"decision.dosage_arms": "decoy.arms"}.get(setting, setting)
    with pytest.raises(ValueError, match=rf"registered rule P1-GO v1.3 in 1 setting\(s\): {setting} = "):
        _quiet(validate_config, merge_config(cfg), adata, log=lambda *a: None)
    log = []
    with pytest.warns(UserWarning, match=rf"labelled 'P1-GO v1.3 \(unregistered variant: {setting}="):
        cfg = validate_config(merge_config({**cfg, "unregistered_variant": True}), adata, log=log.append)
    assert [d[0] for d in config_deviations(cfg)] == [setting]


def test_failed_fits_are_warned_about():
    log = []
    with pytest.warns(UserWarning, match="decoy runs: 2 of 600 NMF fits were stopped by the restart limit"):
        pipeline._warn_failed_fits("decoy runs", {"n_runs": 10, "n_fits": 600, "n_capped_fits": 9,
                                                  "n_failed_fits": 2}, log.append)
    assert log and log[0].startswith("[p1] warning: decoy runs: 2 of 600")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        pipeline._warn_failed_fits("decoy runs", {"n_fits": 600, "n_failed_fits": 0}, log.append)
        pipeline._warn_failed_fits("instability runs", None, log.append)
    runs = pd.DataFrame({"n_fits": [10, 20], "n_capped_fits": [1, 0], "n_failed_fits": [0, 3]})
    assert pipeline._fit_totals(runs) == {"n_runs": 2, "n_fits": 30, "n_capped_fits": 1, "n_failed_fits": 3}
    assert pipeline._fit_totals(None) is None and pipeline._fit_totals(pd.DataFrame({"x": [1]})) is None


def test_cli_refuses_a_deviation_before_reading_the_cohort(tmp_path, capsys):
    missing = str(tmp_path / "no_such_cohort.h5ad")
    path = tmp_path / "cfg.json"
    path.write_text(json.dumps({**BASE, "tme": {"rho_sensitivity": None}}))
    with pytest.raises(SystemExit):
        cli.main(["p1", missing, "--config", str(path), "--out", str(tmp_path / "out")])
    err = capsys.readouterr().err
    assert "tme.rho_sensitivity = None (registered: 1.5)" in err and "unregistered_variant: true" in err
    # allowed as a variant, the check passes and the run stops only at the missing cohort file
    path.write_text(json.dumps({**BASE, "tme": {"rho_sensitivity": None}, "unregistered_variant": True}))
    with pytest.raises((FileNotFoundError, OSError)):
        cli.main(["p1", missing, "--config", str(path), "--out", str(tmp_path / "out")])


def test_tme_and_label_warnings(small, registered, tmp_path):
    adata, truth = small
    with pytest.warns(UserWarning) as rec:
        validate_config(merge_config(registered), adata, log=lambda *a: None)
    assert [str(w.message).split(" (")[0] for w in rec] == ["no expected-TME-shared MPs are pre-declared"]
    # pre-declared gene sets silence it
    declared = {**registered, "decision": {"expected_tme_shared_sets": {"IFN": truth["sets"]["IFN"]}}}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        validate_config(merge_config(declared), adata, log=lambda *a: None)
    # MPs from a GMT without labels: a warning (the published 3CA v2 list already leaves the curated artefacts out);
    # re-derived MPs without labels: a deviation (test_each_deviation_is_refused_by_default_and_labelled_when_allowed)
    write_gmt(truth["sets"], tmp_path / "published.gmt")
    gmt = {**declared, "mps": str(tmp_path / "published.gmt"), "typer": {"labels": None, "labelled_gmt": None}}
    with pytest.warns(UserWarning, match="fine only if the mps GMT already leaves them out"):
        validate_config(merge_config(gmt), adata, log=lambda *a: None)
    no_cells = {**registered, "keys": {k: v for k, v in KEYS.items() if k != "cell_type"}}
    with warnings.catch_warnings():  # without cell types the TME axis cannot run: no TME warning
        warnings.simplefilter("error")
        validate_config(merge_config(no_cells), adata, log=lambda *a: None)


def test_exclusions_that_cannot_run_are_deviations(small, registered, tmp_path, monkeypatch):
    """Positional on without gene positions: refused (before any NMF) or labelled, like positional.enabled false."""
    adata, _ = small
    bare = adata.copy()
    bare.var = bare.var.drop(columns=["chromosome", "start"])
    with pytest.raises(ValueError, match=r"positional.gtf = None \(registered: gene positions"):
        _quiet(validate_config, merge_config(registered), bare, log=lambda *a: None)
    # a GTF that places no var_name is found when it is read, still before any NMF and before the output directory
    no_names = tmp_path / "no_names.gtf"
    no_names.write_text('chr1\tsrc\tgene\t11\t500\t.\t+\t.\tgene_id "ENSG1";\n')
    cfg = {**registered, "positional": {**registered["positional"], "gtf": str(no_names)}}
    monkeypatch.setattr(pipeline, "run_gavish", _no_nmf)
    with pytest.raises(ValueError, match=r"no gene positions for the positional axis: positional.gtf = "):
        _quiet(run_p1, bare, cfg, tmp_path / "out", log=lambda *a: None)
    assert not (tmp_path / "out").exists()
    cfg = {**cfg, "unregistered_variant": True}
    with pytest.raises(AssertionError, match="NMF ran"), pytest.warns(UserWarning) as rec:
        run_p1(bare, cfg, tmp_path / "out", log=lambda *a: None)
    assert any("unregistered variant: the verdict will be labelled 'P1-GO v1.3 (unregistered variant: "
               "positional.gtf=" in str(w.message) for w in rec)


def test_positions_and_cytoband_must_place_genes_on_arms(small, registered, tmp_path, monkeypatch):
    """Inputs that are given but place no gene on a chromosome arm also leave the arm-enrichment exclusion inert:
    refused before any NMF, or labelled as a variant."""
    adata, truth = small
    monkeypatch.setattr(pipeline, "run_gavish", _no_nmf)
    # a GTF whose gene names match no var_name (e.g. Ensembl IDs as var_names), and no positions in adata.var
    bare = adata.copy()
    bare.var = bare.var.drop(columns=["chromosome", "start"])
    other = tmp_path / "other_names.gtf"
    other.write_text('chr7\tsrc\tgene\t11\t500\t.\t+\t.\tgene_id "ENSG1"; gene_name "NOT_A_VAR_NAME";\n')
    cfg = {**registered, "positional": {**registered["positional"], "gtf": str(other)}}
    with pytest.raises(ValueError, match=r"no gene positions for the positional axis: positional.gtf = "), \
            pytest.warns(UserWarning, match="places no var_name"):
        run_p1(bare, cfg, tmp_path / "out", log=lambda *a: None)
    # a cytoband without centromere (acen) bands, e.g. a gap table passed by mistake: no gene gets an arm
    no_acen = {**registered, "positional": {"cytoband": truth["cytoband"].assign(stain="gneg")}}
    log = []
    with pytest.raises(ValueError, match=r"the cytoband places no gene on a chromosome arm .*: positional.cytoband = "
                                         r"'<table>' \(registered: a cytoband whose acen bands place the genes"):
        _quiet(run_p1, adata, no_acen, tmp_path / "out", log=log.append)
    assert "[p1] chromosome arms: 0 of 700 var_names placed" in log and not (tmp_path / "out").exists()
    with pytest.raises(AssertionError, match="NMF ran"), pytest.warns(UserWarning) as rec:
        run_p1(adata, {**no_acen, "unregistered_variant": True}, tmp_path / "out", log=lambda *a: None)
    assert any("(unregistered variant: positional.cytoband='<table>')" in str(w.message) for w in rec)
    # the registered inputs place every gene, still before any NMF
    log = []
    with pytest.raises(AssertionError, match="NMF ran"):
        _quiet(run_p1, adata, registered, tmp_path / "out", log=log.append)
    assert "[p1] chromosome arms: 700 of 700 var_names placed" in log


def test_protocol_axis_that_cannot_pair_a_tumour_is_warned_about(small, registered):
    adata, _ = small
    cases = {"protocol.reference and protocol.test not set": {"reference": None, "test": None},
             "protocol.test not set": {"test": None},
             "protocol.reference 'Fresh' not a value of adata.obs['protocol']": {"reference": "Fresh"}}
    for message, change in cases.items():
        cfg = {**registered, "protocol": {**registered["protocol"], **change}}
        with pytest.warns(UserWarning, match=re.escape(message) + ".*criterion \\(c\\) cannot be powered"):
            validate_config(merge_config(cfg), adata, log=lambda *a: None)
    # names present but never in one tumour
    unpaired = adata.copy()
    tum = unpaired.obs["tumour"].astype(str)
    unpaired.obs["protocol"] = np.where(tum.isin(sorted(tum.unique())[:3]), "fresh", "fixed")
    with pytest.warns(UserWarning, match="no tumour has malignant cells under both protocol.reference 'fresh'"):
        validate_config(merge_config(registered), unpaired, log=lambda *a: None)
    # the protocol axis off, or no protocol key: nothing to warn about
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        declared = {**registered, "decision": {"expected_tme_shared_sets": {"X": ["G1"]}}}
        validate_config(merge_config({**declared, "protocol": {"enabled": False}}), adata, log=lambda *a: None)
        no_key = {**declared, "keys": {k: v for k, v in KEYS.items() if k != "protocol"}}
        validate_config(merge_config(no_key), adata, log=lambda *a: None)


def test_skipped_axes_are_logged_and_noted_in_the_summary(small, registered, tmp_path):
    """An axis switched on that cannot run says so in the log, summary.md and run_config.json (it is not compared
    with the registered rule, but leaves its criterion unpowered)."""
    adata, truth = small
    write_gmt({"cycle": truth["sets"]["cycle"]}, tmp_path / "mps.gmt")
    cfg = {**registered, "mps": str(tmp_path / "mps.gmt"), "keys": {k: v for k, v in KEYS.items() if k != "cell_type"},
           "decoy": {"enabled": False}, "instability": {"enabled": False}, "typer": {"enabled": False},
           "protocol": {"reference": None, "test": None}, "unregistered_variant": True}
    log = []
    _quiet(run_p1, adata, cfg, tmp_path / "out", log=log.append)
    assert "[p1] TME axis skipped (needs keys.cell_type)" in log
    assert any(m.startswith("[p1] protocol axis skipped (needs keys.protocol, protocol.reference and protocol.test; "
                            "got 'protocol', None, None)") for m in log)
    summary = (tmp_path / "out" / "summary.md").read_text()
    assert "The TME axis did not run (needs keys.cell_type), so criterion (b) is not powered." in summary
    assert "The protocol axis did not run" in summary and "Powered criteria: none." in summary
    assert "place 700 of 700 var_names (100%) on chromosome arms" in summary
    run_cfg = json.loads((tmp_path / "out" / "run_config.json").read_text())
    assert any(n.startswith("The protocol axis did not run") for n in run_cfg["notes"])
    assert set(run_cfg["nmf_cache_code"]) == {"real", "decoy"}


def test_decoy_arm_columns_with_nan_are_refused(small):
    adata, _ = small
    adata = adata.copy()
    keys = Keys(**KEYS)
    mal = malignant_mask(adata, keys)
    samples = adata.obs["sample"].astype(str)
    first = samples[mal].iloc[0]
    frac = adata.obs["doublet_est"].to_numpy(dtype=float).copy()
    frac[mal & (samples == first).to_numpy()] = np.nan   # e.g. intronic_fraction missing for one study
    adata.obs["intronic_fraction"] = frac
    cfg = merge_config({**BASE, "decoy": {"arms": {"technical": ["intronic_fraction"], "dosage": ["clone_est"]}}})
    with pytest.raises(ValueError, match=rf"intronic_fraction \({int((mal & (samples == first)).sum())} of "
                                         rf"{int(mal.sum())} malignant cells, in 1 sample\(s\): {first}\)"):
        _quiet(validate_config, cfg, adata, log=lambda *a: None)
    # a constant for every cell of the sample without the measurement is accepted (it conditions nothing there)
    frac[mal & (samples == first).to_numpy()] = 0.0
    adata.obs["intronic_fraction"] = frac
    check_arm_values(adata, keys, {"technical": ["intronic_fraction"]})
    # NaN outside malignant cells (e.g. clone calls only for malignant cells) and in a categorical column
    clone = adata.obs["clone_est"].astype(object).where(mal, None)
    adata.obs["clone_est"] = clone
    check_arm_values(adata, keys, {"dosage": ["clone_est"]})
    clone[np.flatnonzero(mal)[:3]] = None
    adata.obs["clone_est"] = clone
    with pytest.raises(ValueError, match=r"clone_est \(3 of .* explicit label \(e.g. 'unassigned'\)"):
        check_arm_values(adata, keys, {"dosage": ["clone_est"]})
    # text placeholders of a missing value, e.g. after astype(str) or a CSV round trip
    for placeholder in ("nan", " NaN ", "None", ""):
        text = adata.obs["doublet_est"].astype(str).where(mal, "x")
        text.iloc[np.flatnonzero(mal)[:2]] = placeholder
        adata.obs["clone_text"] = text
        with pytest.raises(ValueError, match=r"text placeholder .*: clone_text \(2 of "):
            check_arm_values(adata, keys, {"dosage": ["clone_text"]})
    adata.obs["clone_text"] = adata.obs["doublet_est"].astype(str)
    check_arm_values(adata, keys, {"dosage": ["clone_text"]})


def test_label_and_name_checks_against_audited_mps(tmp_path):
    mps = {"MP1 Cell Cycle": ["A"], "MP2 Stress": ["B"]}
    lab = tmp_path / "labels.csv"
    log = []
    pd.DataFrame({"mp": ["MP9 Ribosome", "MP10 Mito"], "label": [1, 1]}).to_csv(lab, index=False)
    cfg = merge_config({"mps": "published.gmt", "typer": {"labels": str(lab)}})
    with pytest.raises(ValueError, match="match none of the 2 audited MPs"):
        check_names(mps, cfg, log.append)
    pd.DataFrame({"mp": ["MP2 Stress", "MP10 Mito"], "label": [1, 1]}).to_csv(lab, index=False)
    cfg["decision"]["expected_tme_shared"] = ["MP17 IFN"]
    with pytest.warns(UserWarning) as rec:
        check_names(mps, cfg, log.append)
    assert any("MP10 Mito" in str(w.message) for w in rec)
    assert any("expected_tme_shared name(s) match no audited MP: ['MP17 IFN']" in str(w.message) for w in rec)
    # re-derived MPs: labels by name only for this exact run
    pd.DataFrame({"mp": ["MP2"], "label": [1]}).to_csv(lab, index=False)
    with pytest.warns(UserWarning, match="assigned by the run in no fixed order"):
        check_names({"MP1": ["A"], "MP2": ["B"]}, merge_config({"typer": {"labels": str(lab)}}), log.append)


def test_curated_labels_carried_by_gene_overlap():
    genes = [f"G{i}" for i in range(200)]
    ribo, cycle = genes[:40], genes[100:140]
    labelled = {"MP_Ribo": ribo, "MP_Cycle": cycle, "MP_Unlabelled": genes[150:190]}
    labels = pd.Series({"MP_Ribo": 1, "MP_Cycle": 0})
    rederived = {"MP1": ribo[:20] + genes[60:65],   # half of the ribosomal artefact
                 "MP2": ribo[20:] + genes[65:70],   # the other half
                 "MP3": cycle[:30],
                 "MP4": genes[150:180],             # matches a set without a label
                 "MP_Ribo": genes[80:100]}          # shares a name with a labelled set, not its genes
    lab, table = curated_labels(rederived, labels, labelled, rederive=True, min_jaccard=0.3, log=lambda *a: None)
    assert lab["MP1"] == 1 and lab["MP2"] == 1    # both halves: best match per MP, not one-to-one
    assert lab["MP3"] == 0 and np.isnan(lab["MP4"]) and np.isnan(lab["MP_Ribo"])
    assert table.set_index("mp").loc["MP1", "labelled_mp"] == "MP_Ribo"
    # MPs from a GMT: a name in the labels wins, the others go by gene overlap
    lab, table = curated_labels({"MP_Cycle": ribo, "X": ribo}, labels, labelled, rederive=False,
                                log=lambda *a: None)
    assert lab["MP_Cycle"] == 0 and lab["X"] == 1
    assert table.set_index("mp")["source"].to_dict() == {"MP_Cycle": "name", "X": "gene overlap"}
    with pytest.raises(ValueError, match="no set of typer.labelled_gmt has a label"):
        curated_labels(rederived, pd.Series({"other": 1}), labelled, rederive=True, log=lambda *a: None)


def test_positions_from_an_empty_or_sparse_gtf(small, tmp_path):
    adata, _ = small
    pc = dict(merge_config({})["positional"])
    empty = tmp_path / "empty.gtf"
    empty.write_text('chr1\tsrc\tgene\t11\t500\t.\t+\t.\tgene_id "ENSG1";\n')
    with pytest.warns(UserWarning, match="gave no gene positions"):
        assert pipeline._gtf_positions(adata, {**pc, "gtf": str(empty)}, lambda *a: None) is None
    sparse = tmp_path / "sparse.gtf"
    sparse.write_text("".join(f'chr1\tsrc\tgene\t{100 * i + 1}\t{100 * i + 50}\t.\t+\t.\tgene_name "{g}";\n'
                              for i, g in enumerate(adata.var_names[:100])))
    with pytest.warns(UserWarning, match="places only 14% of var_names"):
        table = pipeline._gtf_positions(adata, {**pc, "gtf": str(sparse)}, lambda *a: None)
    assert len(table) == 100
    with pytest.warns(UserWarning, match="gene lengths from the GTF cover only 14%"):
        pipeline._gene_lengths(adata, table, lambda *a: None)
    var = adata.var.copy()
    try:
        adata.var["chromosome"] = adata.var["chromosome"].where(np.arange(adata.n_vars) < 300)
        with pytest.warns(UserWarning, match="place only 43% of var_names"):
            pos = pipeline._var_positions(adata, pc, lambda *a: None)
        assert len(pos) == 300
    finally:
        adata.var = var


def test_residual_decoy_keeps_covariate_driven_coexpression(small):
    # estimated covariates, as an analyst would have them; the planted doublet_score is coupled to
    # cycling by design (see test_decoy_arms.py for that false-recruitment check)
    adata, truth = small
    keys = Keys(**KEYS)
    dec = make_decoy(adata, keys, counts_layer="counts", strata=("clone_est", "doublet_est", "contamination_est"),
                     seed=4)
    X, D = get_counts(adata, "counts"), get_counts(dec, "counts")
    pos = {g: i for i, g in enumerate(adata.var_names)}

    def k(M, name):
        idx = [pos[g] for g in truth["sets"][name]]
        vals = []
        for s in adata.obs["sample"].unique()[:4]:
            Y = M[_sample_rows(adata, keys, s)]
            vals.append(kappa(Y[:, idx].toarray(), np.asarray(Y.sum(axis=1)).ravel()))
        return np.nanmean(vals)

    assert k(D, "contamination") > 0.6 * k(X, "contamination")  # continuous covariate kept
    assert k(D, "cycle") < 0.2 * k(X, "cycle")                  # real program destroyed
