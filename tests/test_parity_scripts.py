"""Round trip of the R-versus-port parity scripts (scripts/parity) on a tiny simulated cohort.

R is not needed: the "R output" is mpaudit's own fit written in run_r_snmf.R's file format, so
compare.py must report perfect parity. The last test runs the real R script when Rscript and the
NMF package are installed.
"""
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit.data import Keys, malignant_mask
from mpaudit.gavish import GavishParams, compute_programs, prepare_matrix
from mpaudit.sim import simulate_cohort

PARITY = Path(__file__).resolve().parents[1] / "scripts" / "parity"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(f"parity_{name}", PARITY / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


export_matrices = _load("export_matrices")
parity_compare = _load("compare")
sys.path.insert(0, str(PARITY))  # (the scripts add it too)
import parity_io as pio  # noqa: E402

CONFIG = {
    "keys": {"sample": "sample", "malignant": "malignant", "study": "study"},
    "counts_layer": "counts",
    "recurrence_key": "tumour",
    "gavish": {"ranks": [4, 5, 6], "nrun": 2, "n_universe": 400, "min_cells": 20, "min_genes_per_cell": 0},
}


def _write_r_style_outputs(r_out: Path, manifest: dict, programs) -> None:
    """Bases, fits, input check and versions exactly as run_r_snmf.R writes them (write.csv quoting)."""
    for e in manifest["samples"]:
        d = r_out / e["dir"]
        d.mkdir(parents=True, exist_ok=True)
        spg = programs[e["sample"]]
        for r, W in spg.W.items():
            lines = [",".join(['"gene"'] + [f'"{j + 1}"' for j in range(W.shape[1])])]
            lines += [",".join([f'"{g}"'] + ["%.17g" % v for v in row]) for g, row in zip(spg.universe, np.asarray(W, float))]
            (d / f"W_rank{r}.csv").write_text("\n".join(lines) + "\n")
        fits = ['"dir","rank","nrun","objective","n_iter","elapsed_sec","seed"']
        fits += [f'"{e["dir"]}",{r},{manifest["gavish"]["nrun"]},"NA",NA,NA,{manifest["r_seed"]}' for r in spg.W]
        (d / "fits.csv").write_text("\n".join(fits) + "\n")
        (d / "input_check.csv").write_text('"n_genes","n_cells","sum","max"\n'
                                           f'{e["n_genes"]},{e["n_cells"]},"{e["sum"]:.17g}","{e["max"]:.17g}"\n')
    (r_out / "versions.csv").write_text('"key","value"\n"R","R version 0.0.0 (mpaudit stand-in)"\n'
                                        '"NMF","none"\n"seed","123456"\n')


@pytest.fixture(scope="module")
def parity_run(tmp_path_factory):
    root = tmp_path_factory.mktemp("parity")
    adata, _ = simulate_cohort(n_tumours=6, n_twin=0, n_region=0, n_mal=100, n_tme_per_type=5, n_genes=600, seed=0)
    adata.write_h5ad(root / "cohort.h5ad")
    (root / "config.json").write_text(json.dumps(CONFIG))
    manifest = export_matrices.main([str(root / "cohort.h5ad"), "--config", str(root / "config.json"),
                                     "--out", str(root / "export"), "--n", "5"])
    params = pio.params_from_manifest(manifest)
    programs = compute_programs(adata, Keys.from_dict(manifest["keys"]), params, counts_layer="counts",
                                samples=[e["sample"] for e in manifest["samples"]])
    _write_r_style_outputs(root / "r_out", manifest, programs)
    summary = parity_compare.main(["--export", str(root / "export"), "--r-out", str(root / "r_out"),
                                   "--out", str(root / "report")])
    return {"root": root, "adata": adata, "manifest": manifest, "params": params, "summary": summary}


def test_export_writes_the_exact_nmf_input(parity_run):
    root, adata, manifest, params = (parity_run[k] for k in ("root", "adata", "manifest", "params"))
    assert len(manifest["samples"]) == 5 and manifest["gavish"]["ranks"] == [4, 5, 6]
    assert (root / "export" / "samples.tsv").exists() and (root / "export" / "r_params.tsv").exists()
    kv = pd.read_csv(root / "export" / "r_params.tsv", sep="\t", dtype=str).set_index("key")["value"]
    assert kv["ranks"] == "4,5,6" and kv["nrun"] == "2" and kv["method"] == "snmf/r"
    counts = adata.layers["counts"]
    for e in manifest["samples"]:
        mal = malignant_mask(adata, Keys.from_dict(manifest["keys"]))
        idx = np.flatnonzero((adata.obs["sample"].astype(str) == e["sample"]).to_numpy() & mal)
        A, universe = prepare_matrix(sp.csr_matrix(counts[idx]), np.asarray(adata.var_names), params)
        B = pio.read_matrix(root / "export", e)
        assert np.array_equal(A, B)  # shortest round-trip decimals: identical doubles
        assert pio.read_lines(root / "export" / e["dir"] / "genes.txt") == list(universe)
        cells = pio.read_lines(root / "export" / e["dir"] / "cells.txt")
        assert len(cells) == B.shape[1] == e["n_cells"] and set(cells) <= set(adata.obs_names[idx])
        assert e["unit"] == adata.obs["tumour"].iloc[idx[0]]


def test_self_round_trip_reports_perfect_parity(parity_run):
    s, root = parity_run["summary"], parity_run["root"]
    assert s["verdict"] == "PASS" and not s["issues"]
    m = s["r_vs_port"]
    assert m["median_overlap"] == m["min_overlap"] == 50 and m["median_cosine"] == pytest.approx(1.0)
    assert m["n_robust_a"] == m["n_robust_b"] == m["n_robust_reproduced"] > 0
    assert m["mp_tested"] and m["n_mps_a"] == m["n_mps_b"] == m["n_mps_identical"] >= 1
    assert s["objective"]["median_py_vs_r_pct"] == pytest.approx(0.0) and s["objective"]["n_port_worse"] == 0
    assert m["robust_agreement"] == m["mp_agreement"] == 1.0
    assert s["r_versions"]["R"].startswith("R version")
    mps = pd.read_csv(root / "report" / "mps.csv")
    assert (mps["member_jaccard"] == 1).all() and (mps["n_members_a"] == mps["n_members_b"]).all()
    programs = pd.read_csv(root / "report" / "programs.csv")
    assert len(programs) == 5 * (4 + 5 + 6) and (programs["overlap"] == 50).all()
    report = (root / "report" / "report.md").read_text()
    assert "**Verdict: PASS**" in report and "How to read this report" in report
    # the port's bases are written in the R layout and load back unchanged
    reloaded = pio.read_bases(root / "report" / "python_bases", parity_run["manifest"], [4, 5, 6])
    assert set(reloaded) == {e["sample"] for e in parity_run["manifest"]["samples"]}


def test_scrambled_r_bases_fail(parity_run, tmp_path):
    root, manifest = parity_run["root"], parity_run["manifest"]
    r_bad = tmp_path / "r_bad"
    shutil.copytree(root / "r_out", r_bad)
    rng = np.random.default_rng(1)
    for e in manifest["samples"]:
        for r in (4, 5, 6):
            path = r_bad / e["dir"] / f"W_rank{r}.csv"
            df = pd.read_csv(path)
            df.iloc[:, 1:] = df.iloc[rng.permutation(len(df)), 1:].to_numpy()  # genes no longer match weights
            df.to_csv(path, index=False)
    s = parity_compare.main(["--export", str(root / "export"), "--r-out", str(r_bad), "--out", str(tmp_path / "rep"),
                             "--python-bases", str(root / "report" / "python_bases"), "--no-objective"])
    assert s["verdict"] == "FAIL" and not s["issues"]
    assert s["r_vs_port"]["median_overlap"] < 45


def test_matrix_mismatch_is_invalid(parity_run, tmp_path):
    root, manifest = parity_run["root"], parity_run["manifest"]
    r_bad = tmp_path / "r_bad"
    shutil.copytree(root / "r_out", r_bad)
    check = r_bad / manifest["samples"][0]["dir"] / "input_check.csv"
    c = pd.read_csv(check)
    c["sum"] = c["sum"] * (1 + 1e-6)
    c.to_csv(check, index=False)
    s = parity_compare.main(["--export", str(root / "export"), "--r-out", str(r_bad), "--out", str(tmp_path / "rep"),
                             "--python-bases", str(root / "report" / "python_bases")])
    assert s["verdict"] == "INVALID" and "sum/max" in s["issues"][0]


def test_verdict_rules():
    def summ(median, robust=1.0, mp=1.0, mp_tested=True):
        return {"summary": {"median_overlap": median, "robust_agreement": robust, "mp_agreement": mp,
                            "program_ok": median >= 45, "robust_ok": robust == 1.0, "mp_ok": mp == 1.0,
                            "mp_tested": mp_tested}}

    verdict = parity_compare.verdict
    assert verdict(summ(47), {}, [])[0] == "PASS"
    assert "not tested" in verdict(summ(47, mp=0.0, mp_tested=False), {}, [])[1]
    rep = {"port vs port": summ(46, robust=0.8, mp=0.5)}
    assert verdict(summ(44, robust=0.8, mp=0.5), rep, [])[0] == "WITHIN NOISE"
    assert verdict(summ(44, robust=0.7, mp=0.5), rep, [])[0] == "FAIL"     # robust worse than seed noise
    assert verdict(summ(43, robust=0.9, mp=0.6), rep, [])[0] == "FAIL"     # > 2 genes below the replicate
    assert "no seed replicate" in verdict(summ(44), {}, [])[1]
    assert verdict(summ(50), rep, ["R's matrix sum/max for S differ"])[0] == "INVALID"


def test_kept_cells_names_the_columns_prepare_matrix_keeps():
    rng = np.random.default_rng(0)
    counts = rng.poisson(rng.gamma(2, 2, 30)[None, :], (40, 30)).astype(float)
    counts[:, 29] = 0
    counts[5] = 0
    counts[5, 29] = 50  # cell 5 expresses only a gene outside the universe: all zero after centring
    genes = np.array([f"g{i}" for i in range(30)])
    params = GavishParams(n_universe=20, min_genes_per_cell=0)
    A, universe = prepare_matrix(sp.csr_matrix(counts), genes, params)
    keep = export_matrices.kept_cells(sp.csr_matrix(counts), genes, universe, A, params)
    assert "g29" not in universe and not keep[5] and keep.sum() == A.shape[1] < 40


def test_r_script_matches_the_export_format():
    """run_r_snmf.R cannot run here; check that it uses the file names and call the Python side expects."""
    text = (PARITY / "run_r_snmf.R").read_text()
    for name in (pio.MATRIX, pio.GENES, pio.CELLS, pio.SAMPLES_TSV, pio.R_PARAMS_TSV, pio.FITS, pio.INPUT_CHECK,
                 "W_rank%d.csv", "versions.csv"):
        assert f'"{name}"' in text, name
    assert 'method = "snmf/r"' in text and "nrun = nrun" in text and "seed = seed" in text


def test_scripts_run_from_the_command_line():
    for script in ("export_matrices.py", "compare.py"):
        res = subprocess.run([sys.executable, str(PARITY / script), "--help"], capture_output=True, text=True)
        assert res.returncode == 0 and "usage" in res.stdout, res.stderr


def _r_with_nmf() -> bool:
    if shutil.which("Rscript") is None:
        return False
    return subprocess.run(["Rscript", "-e", "library(NMF); library(Matrix)"], capture_output=True).returncode == 0


@pytest.mark.skipif(not _r_with_nmf(), reason="Rscript with the NMF and Matrix packages is not installed")
def test_real_r_side_reads_the_export(parity_run, tmp_path):
    root = parity_run["root"]
    res = subprocess.run(["Rscript", str(PARITY / "run_r_snmf.R"), str(root / "export"), "--out", str(tmp_path / "r"),
                          "--cores", "1"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    s = parity_compare.main(["--export", str(root / "export"), "--r-out", str(tmp_path / "r"),
                             "--out", str(tmp_path / "rep"), "--python-bases", str(root / "report" / "python_bases")])
    assert not s["issues"] and s["verdict"] in {"PASS", "WITHIN NOISE", "FAIL"}
