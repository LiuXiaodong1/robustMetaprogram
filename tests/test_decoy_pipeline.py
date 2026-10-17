import json

import numpy as np
import pandas as pd
import pytest

from mpaudit import Keys, pipeline
from mpaudit.coexpr import kappa
from mpaudit.compare import jaccard_matrix
from mpaudit.data import get_counts, malignant_mask, write_gmt
from mpaudit.decoy import make_decoy
from mpaudit.pipeline import check_names, curated_labels, merge_config, run_p1, validate_config
from mpaudit.sim import simulate_cohort

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
    }
    out_dir = tmp_path / "out"
    with pytest.warns(UserWarning, match="gave no gene positions; falling back to adata.var"):
        out = run_p1(adata, cfg, out_dir, log=lambda *a: None)
    d = out["decision"]
    assert d["rule_version"] == "P1-GO v1.1"
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
    assert "INCOMPLETE" in (out_dir / "summary.md").read_text()


def _no_nmf(*args, **kwargs):
    raise AssertionError("NMF ran before the configuration was checked")


BASE = {"keys": KEYS, "counts_layer": "counts", "recurrence_key": "tumour",
        "decoy": {"arms": {"technical": ["doublet_est"], "dosage": ["clone_est"]}}}


@pytest.mark.parametrize("change, message", [
    ({"decoy": {"strata": ["clone", "doublet_score"], "arms": {}}}, "decoy.strata is no longer accepted"),
    ({"decoy": {"arms": {"tech": ["doublet_est"], "cna": ["clone_est"]}}}, r"\['cna', 'tech'\] are in neither"),
    ({"decoy": {"arms": {"technical": ["doublet_est", "no_such_column"]}}}, "not columns of adata.obs"),
    ({"decision": {"recruitment": 0.3, "partial_r2": 0.3}}, "before P1-GO v1.0; P1-GO v1.1 replaces"),
    ({"decision": {"technical_arms": ["technical"], "dosage_arms": ["technical"]}}, "both technical_arms"),
    ({"decision": {"expected_tme_shared": ["MP17"]}}, "cannot be pre-declared with mps: rederive"),
    ({"decision": {"expected_tme_shared_gmt": "missing.gmt"}}, "is not a file"),
    ({"typer": {"labels": "LABELS"}}, "cannot match re-derived MPs"),
    ({"typer": {"labelled_gmt": "LABELLED"}}, "needs typer.labels"),
    ({"tme": {"rho_sensitivity": 0}}, "rho_sensitivity must be > 0"),
    ({"protocol": {"kappa_min_effect": "large"}}, "kappa_min_effect must be a number"),
    ({"mps": "missing.gmt"}, "mps must be 'rederive' or a GMT file"),
    ({"mps": "PUBLISHED", "typer": {"labels": "LABELS"}}, "typer.labels names match none of the 1 audited MPs"),
])
def test_bad_config_fails_before_any_nmf(small, tmp_path, monkeypatch, change, message):
    adata, truth = small
    pd.DataFrame({"mp": ["MP1 Cell Cycle - G2/M", "MP5 Stress"], "label": [0, 1]}).to_csv(tmp_path / "l.csv",
                                                                                         index=False)
    write_gmt({"MP5 Stress": truth["sets"]["stress"]}, tmp_path / "l.gmt")
    write_gmt({"Stress": truth["sets"]["stress"]}, tmp_path / "published.gmt")
    files = {"LABELS": str(tmp_path / "l.csv"), "LABELLED": str(tmp_path / "l.gmt"),
             "PUBLISHED": str(tmp_path / "published.gmt")}

    def sub(x):
        return {k: sub(v) for k, v in x.items()} if isinstance(x, dict) else files.get(x, x) if isinstance(x, str) else x

    cfg = {k: dict(v) if isinstance(v, dict) else v for k, v in BASE.items()}
    for k, v in sub(change).items():
        cfg[k] = {**cfg[k], **v} if isinstance(cfg.get(k), dict) else v
    monkeypatch.setattr(pipeline, "run_gavish", _no_nmf)
    with pytest.raises(ValueError, match=message):
        run_p1(adata, cfg, tmp_path / "out", log=lambda *a: None)
    assert not (tmp_path / "out").exists()


def test_decoy_arm_checks(small):
    adata, _ = small
    cfg = merge_config({**BASE, "decoy": {"arms": {"technical": "doublet_est"}}})
    with pytest.warns(UserWarning, match="no decoy arm is a dosage arm"):
        cfg = validate_config(cfg, adata, log=lambda *a: None)
    assert cfg["decoy"]["arms"] == {"technical": ["doublet_est"]}  # a scalar column becomes a list
    cfg = merge_config({**BASE, "decoy": {"arms": {"technical": ["pct_ribo"], "cna": ["clone_est"]}},
                        "decision": {"dosage_arms": ["cna"]}})
    validate_config(cfg, adata, log=lambda *a: None)  # renamed arm listed in the rule; QC columns come later
    cfg = merge_config({**BASE, "decoy": {"enabled": False, "arms": {"tech": ["doublet_est"]}}})
    validate_config(cfg, adata, log=lambda *a: None)  # decoys off: arm names are not used


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
