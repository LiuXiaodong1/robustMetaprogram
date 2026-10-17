"""Per-cell QC covariates, GTF gene lengths and the typer's gene-length / %ribo / intronic features."""
import gzip
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit import typer
from mpaudit.data import LONG_GENE_BP, Keys, add_qc_obs, qc_metrics
from mpaudit.positional import assign_arms, gene_lengths, gene_positions_from_gtf
from mpaudit.sim import simulate_cohort
from mpaudit.typer import typer_features

QC_GENES = ["MT-CO1", "RPL3", "RPLP0", "RPS6KA1", "RPL13AP5", "GAPDH"]
QC_COUNTS = np.array([[2, 3, 1, 4, 5, 5],   # 20 counts, 6 genes; mito 2, ribo RPL3 + RPLP0 = 4
                      [0, 0, 0, 0, 0, 0],   # empty cell
                      [0, 5, 0, 0, 0, 5]])  # 10 counts, 2 genes; ribo 5
SIM_KEYS = Keys(sample="sample", malignant="malignant", doublet="doublet_score", intronic="intronic_fraction")


def _tiny_adata(genes=QC_GENES) -> ad.AnnData:
    a = ad.AnnData(X=sp.csr_matrix(np.zeros(QC_COUNTS.shape)), var=pd.DataFrame(index=genes))
    a.layers["counts"] = sp.csr_matrix(QC_COUNTS.astype(float))
    return a


def test_qc_metrics_values():
    qc = qc_metrics(sp.csr_matrix(QC_COUNTS), QC_GENES)
    np.testing.assert_allclose(qc["log_umi"], np.log1p([20, 0, 10]))
    np.testing.assert_allclose(qc["log_genes"], np.log1p([6, 0, 2]))
    np.testing.assert_allclose(qc["pct_mito"], [10.0, 0.0, 0.0])
    # RPS6KA1 (a kinase) and RPL13AP5 (a pseudogene) are not ribosomal-protein genes
    np.testing.assert_allclose(qc["pct_ribo"], [20.0, 0.0, 50.0])


def test_add_qc_obs_in_place_and_overwrite():
    a = _tiny_adata()
    a.obs["pct_mito"] = -1.0
    out = add_qc_obs(a, counts_layer="counts")
    assert out is a
    np.testing.assert_allclose(a.obs["pct_mito"], -1.0)  # existing column kept
    np.testing.assert_allclose(a.obs["pct_ribo"], [20.0, 0.0, 50.0])
    np.testing.assert_allclose(a.obs["log_umi"], np.log1p([20, 0, 10]))
    add_qc_obs(a, counts_layer="counts", overwrite=True)
    np.testing.assert_allclose(a.obs["pct_mito"], [10.0, 0.0, 0.0])
    add_qc_obs(a)  # all columns present: adata.X (zeros) is not even read
    np.testing.assert_allclose(a.obs["log_genes"], np.log1p([6, 0, 2]))


def test_add_qc_obs_warns_without_matching_genes():
    a = _tiny_adata([f"ENSG{i:011d}" for i in range(len(QC_GENES))])
    with pytest.warns(UserWarning) as record:
        add_qc_obs(a, counts_layer="counts")
    assert sorted(str(w.message).split("; ")[1].split()[0] for w in record) == ["pct_mito", "pct_ribo"]
    assert (a.obs["pct_mito"] == 0).all() and (a.obs["pct_ribo"] == 0).all()


# genomic spans: GAPDH and the kinase long, RPL3 short, the others unknown (not long)
QC_LENGTHS = pd.Series({"GAPDH": 150_000, "RPS6KA1": LONG_GENE_BP + 1, "RPL3": 5_000, "UNUSED": 900_000})


def test_qc_metrics_pct_long():
    assert typer.LONG_GENE_BP is LONG_GENE_BP  # one constant for the typer and pct_long
    assert "pct_long" not in qc_metrics(sp.csr_matrix(QC_COUNTS), QC_GENES)
    qc = qc_metrics(sp.csr_matrix(QC_COUNTS), QC_GENES, gene_length=QC_LENGTHS)
    np.testing.assert_allclose(qc["pct_long"], [45.0, 0.0, 50.0])  # (4 + 5) / 20; empty cell; 5 / 10


def test_add_qc_obs_pct_long_needs_gene_lengths():
    a = _tiny_adata()
    add_qc_obs(a, counts_layer="counts")
    assert "pct_long" not in a.obs  # no lengths: no column
    a.var["gene_length"] = QC_LENGTHS.reindex(QC_GENES).to_numpy()
    add_qc_obs(a, counts_layer="counts")  # lengths in var are found
    np.testing.assert_allclose(a.obs["pct_long"], [45.0, 0.0, 50.0])
    b = _tiny_adata()
    add_qc_obs(b, counts_layer="counts", gene_length=pd.Series({"RPL3": 200_000}))  # given lengths (e.g. a GTF)
    np.testing.assert_allclose(b.obs["pct_long"], [15.0, 0.0, 50.0])
    c = _tiny_adata()
    with pytest.warns(UserWarning, match="pct_long is 0 for every cell"):
        add_qc_obs(c, counts_layer="counts", gene_length=pd.Series({"ENSG1": 200_000}))  # names do not match


def test_keys_intronic():
    assert Keys().intronic is None
    assert Keys.from_dict({"sample": "s", "intronic": "intronic_fraction"}).intronic == "intronic_fraction"
    with pytest.raises(ValueError):
        Keys.from_dict({"intronic_col": "x"})


GTF = """\
#!genome-build GRCh38
chr1\tHAVANA\tgene\t1001\t2000\t.\t+\t.\tgene_id "G1"; gene_name "AAA";
chr1\tHAVANA\ttranscript\t1001\t1500\t.\t+\t.\tgene_id "G1"; gene_name "AAA"; transcript_id "T1";
chr1\tHAVANA\texon\t1001\t1100\t.\t+\t.\tgene_id "G1"; gene_name "AAA"; transcript_id "T1";
chr7\tHAVANA\tgene\t60000001\t60250000\t.\t-\t.\tgene_id "G2"; gene_type "protein_coding"; gene_name "LONG1";
chrX\tHAVANA\tgene\t500\t799\t.\t+\t.\tgene_id "G3"; gene_name "PAR1";
chrY\tHAVANA\tgene\t900\t1999\t.\t+\t.\tgene_id "G3_PAR_Y"; gene_name "PAR1";
chr2\tHAVANA\tgene\t10\t20\t.\t+\t.\tgene_id "G4";
"""


@pytest.mark.parametrize("gz", [False, True])
def test_gtf_positions_and_lengths(tmp_path, gz):
    path = tmp_path / ("genes.gtf.gz" if gz else "genes.gtf")
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(GTF)
    else:
        path.write_text(GTF)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # gene records present: no fallback, no warning
        pos = gene_positions_from_gtf(path)
    assert list(pos.columns) == ["chromosome", "start", "end", "length"]
    assert list(pos.index) == ["AAA", "LONG1", "PAR1"]  # gene records only; first PAR1 record kept; G4 unnamed
    assert pos.loc["LONG1", "chromosome"] == "7"
    assert (pos.loc["LONG1", "start"], pos.loc["LONG1", "end"], pos.loc["LONG1", "length"]) == (60000001, 60250000, 250000)
    assert pos.loc["AAA", "length"] == 1000 and pos.loc["PAR1", "length"] == 300
    pd.testing.assert_series_equal(gene_lengths(pos), pos["length"].astype(float).rename("gene_length"))
    arms = assign_arms(pos, {"1": 5e3, "7": 5e7})
    assert arms["AAA"] == "1p" and arms["LONG1"] == "7q" and pd.isna(arms["PAR1"])


# UCSC ncbiRefSeq style: transcript and exon lines only, no gene records
UCSC_GTF_2LINE = """\
chr1\tncbiRefSeq.2022-10-28\ttranscript\t11874\t14409\t.\t+\t.\tgene_id "DDX11L1"; transcript_id "NR_046018.2"; gene_name "DDX11L1";
chr1\tncbiRefSeq.2022-10-28\texon\t11874\t12227\t.\t+\t.\tgene_id "DDX11L1"; transcript_id "NR_046018.2"; exon_number "1"; gene_name "DDX11L1";
"""
UCSC_GTF = UCSC_GTF_2LINE + """\
chr7\tncbiRefSeq\ttranscript\t300\t900\t.\t-\t.\tgene_id "TWO"; transcript_id "NM_2.1"; gene_name "TWO";
chr7\tncbiRefSeq\texon\t300\t400\t.\t-\t.\tgene_id "TWO"; transcript_id "NM_2.1"; gene_name "TWO";
chr7\tncbiRefSeq\ttranscript\t100\t500\t.\t-\t.\tgene_id "TWO"; transcript_id "NM_2.2"; gene_name "TWO";
chr7_KI270803v1_alt\tncbiRefSeq\ttranscript\t5\t50\t.\t+\t.\tgene_id "ALT"; transcript_id "NM_3.1"; gene_name "ALT";
chr6\tncbiRefSeq\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ALT"; transcript_id "NM_3.1"; gene_name "ALT";
chr6_GL000250v2_alt\tncbiRefSeq\ttranscript\t1\t99999\t.\t+\t.\tgene_id "ALT"; transcript_id "NM_3.1"; gene_name "ALT";
chr6\tncbiRefSeq\texon\t1900\t2100\t.\t+\t.\tgene_id "ALT"; transcript_id "NM_3.1"; gene_name "ALT";
chrX\tncbiRefSeq\ttranscript\t500\t799\t.\t+\t.\tgene_id "PAR1"; transcript_id "NM_4.1"; gene_name "PAR1";
chrY\tncbiRefSeq\ttranscript\t900\t1999\t.\t+\t.\tgene_id "PAR1"; transcript_id "NM_4.1"; gene_name "PAR1";
chrUn_KI270302v1\tncbiRefSeq\texon\t10\t40\t.\t+\t.\tgene_id "UN"; transcript_id "NM_5.1"; gene_name "UN";
chr2\tncbiRefSeq\tCDS\t70\t90\t.\t+\t0\tgene_id "EXONLY"; transcript_id "NM_6.1"; gene_name "EXONLY";
chr2\tncbiRefSeq\texon\t60\t95\t.\t+\t.\tgene_id "EXONLY"; transcript_id "NM_6.1"; gene_name "EXONLY";
"""


def test_gtf_without_gene_records_aggregates_transcripts(tmp_path):
    path = tmp_path / "ucsc2.gtf"
    path.write_text(UCSC_GTF_2LINE)
    with pytest.warns(UserWarning, match="no 'gene' records"):
        pos = gene_positions_from_gtf(path)
    assert pos.loc["DDX11L1"].tolist() == ["1", 11874, 14409, 2536]

    path = tmp_path / "ucsc.gtf.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(UCSC_GTF)
    with pytest.warns(UserWarning, match="spans of 6 genes"):
        pos = gene_positions_from_gtf(path)
    assert list(pos.columns) == ["chromosome", "start", "end", "length"]
    assert set(pos.index) == {"DDX11L1", "TWO", "ALT", "PAR1", "UN", "EXONLY"}
    assert pos.loc["TWO"].tolist() == ["7", 100, 900, 801]  # min start, max end over both transcripts
    # the primary chromosome wins over alt contigs listed before and after it
    assert pos.loc["ALT"].tolist() == ["6", 1000, 2100, 1101]
    assert pos.loc["PAR1"].tolist() == ["X", 500, 799, 300]  # first chromosome; chrY record ignored
    assert pos.loc["UN"].tolist() == ["Un_KI270302v1", 10, 40, 31]  # only on a contig: kept there
    assert pos.loc["EXONLY", ["start", "end"]].tolist() == [60, 95]  # exon/CDS records only
    pd.testing.assert_series_equal(gene_lengths(pos), pos["length"].astype(float).rename("gene_length"))
    arms = assign_arms(pos, {"7": 5e2, "6": 5e3})
    assert arms["TWO"] == "7p" and arms["ALT"] == "6p" and pd.isna(arms["DDX11L1"])


def test_gtf_without_names_warns(tmp_path):
    path = tmp_path / "ids_only.gtf"
    path.write_text(UCSC_GTF.replace("gene_name", "gene_symbol"))
    with pytest.warns(UserWarning, match="no record carries a gene_name"):
        pos = gene_positions_from_gtf(path)
    assert pos.empty and list(pos.columns) == ["chromosome", "start", "end", "length"]
    with pytest.warns(UserWarning, match="no 'gene' records"):
        assert gene_positions_from_gtf(path, name_attr="gene_symbol").shape[0] == 6


def test_gene_lengths_sources():
    var = pd.DataFrame({"start": [100, 1, 50], "end": [199, 1, 10], "gene_length": [5_000, 7, 0]},
                       index=["a", "b", "c"])
    np.testing.assert_array_equal(gene_lengths(var).to_numpy(), [5_000.0, 7.0, np.nan])  # column wins; 0 -> NaN
    np.testing.assert_array_equal(gene_lengths(var.drop(columns="gene_length")).to_numpy(), [100.0, 1.0, np.nan])
    bm = var.drop(columns="gene_length").rename(columns={"start": "start_position", "end": "end_position"})
    assert gene_lengths(bm, start_col="start_position", end_col="end_position")["a"] == 100.0
    with pytest.raises(ValueError):
        gene_lengths(var[["start"]])


@pytest.fixture(scope="module")
def sim_features():
    adata, truth = simulate_cohort(seed=0)
    feats = typer_features(adata, SIM_KEYS, truth["sets"], counts_layer="counts").set_index("mp")
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(feats.round(3).T.to_string())
    return feats


def test_typer_technical_axes_on_sim(sim_features):
    f = sim_features
    assert f["abs_cor_pct_ribo"].idxmax() == "abundance"
    for col in ("median_gene_length_kb", "frac_long_genes", "abs_cor_intronic"):
        assert f[col].idxmax() == "long_genes", col
    assert f.loc["abundance", "abs_cor_pct_ribo"] > 0.5 and f.loc["long_genes", "abs_cor_intronic"] > 0.4
    assert f.loc["long_genes", "median_gene_length_kb"] > 150 and f.loc["long_genes", "frac_long_genes"] == 1.0
    for bio in ("cycle", "hypoxia"):
        assert f.loc[bio, "median_gene_length_kb"] < 100 and f.loc[bio, "frac_long_genes"] < 0.5
        assert f.loc[bio, "abs_cor_pct_ribo"] < 0.2 and f.loc[bio, "abs_cor_intronic"] < 0.2


def test_typer_inputs_and_qc_consistency():
    adata, truth = simulate_cohort(n_tumours=3, n_twin=0, n_region=0, n_mal=60, n_tme_per_type=10,
                                   n_genes=400, seed=1, estimated_covariates=False)
    sets = {k: truth["sets"][k] for k in ("cycle", "abundance", "long_genes")}
    sets["unplaced"] = ["NOT_A_GENE1", "NOT_A_GENE2"]
    lengths = gene_lengths(adata.var)

    bare = adata.copy()
    bare.var = bare.var.drop(columns=["start", "end", "gene_length"])
    f0 = typer_features(bare, Keys(sample="sample", malignant="malignant"), sets, counts_layer="counts")
    assert not {"median_gene_length_kb", "frac_long_genes", "abs_cor_doublet", "abs_cor_intronic"} & set(f0.columns)

    # explicit lengths on a bare var == lengths taken from adata.var; NaN for a set without known lengths
    f1 = typer_features(bare, SIM_KEYS, sets, counts_layer="counts", gene_length=lengths).set_index("mp")
    f2 = typer_features(adata, SIM_KEYS, sets, counts_layer="counts").set_index("mp")
    pd.testing.assert_frame_equal(f1, f2)
    assert np.isnan(f1.loc["unplaced", "median_gene_length_kb"]) and np.isnan(f1.loc["unplaced", "frac_long_genes"])

    # the typer's internal %ribo is the add_qc_obs column: correlating with it via obs gives the same feature
    add_qc_obs(adata, counts_layer="counts")
    f3 = typer_features(adata, Keys(sample="sample", malignant="malignant", intronic="pct_ribo"), sets,
                        counts_layer="counts").set_index("mp")
    np.testing.assert_allclose(f3["abs_cor_intronic"], f3["abs_cor_pct_ribo"], equal_nan=True)
