"""File formats shared by the R-versus-port parity scripts (``export_matrices.py``, ``run_r_snmf.R``,
``compare.py``).

Export directory (written by ``export_matrices.py``)::

    manifest.json            parameters, input paths, versions, one entry per sample
    samples.tsv, r_params.tsv  the parts of the manifest run_r_snmf.R needs (no JSON parser in base R)
    samples/<dir>/matrix.mtx.gz  NMF input, genes x cells, Matrix Market (shortest round-trip doubles)
    samples/<dir>/genes.txt      row names, one per line
    samples/<dir>/cells.txt      column names, one per line

Basis directory (written by ``run_r_snmf.R`` for R, by ``compare.py`` for the port)::

    <dir>/W_rank<k>.csv      basis W, genes x k: column ``gene``, then ``1`` ... ``k`` (program index)
    <dir>/fits.csv           one row per rank (R: ``objective`` = best-of-nrun sNMF/R objective)
    <dir>/input_check.csv    R only: dimensions, sum and max of the matrix R read

``<dir>`` is the sample's ``dir`` in the manifest (a filesystem-safe name; the manifest maps it back
to the sample name).
"""
from __future__ import annotations

import dataclasses
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from mpaudit.gavish import GavishParams, SamplePrograms

FORMAT = "mpaudit-parity-export"
FORMAT_VERSION = 1
MANIFEST = "manifest.json"
SAMPLES_TSV = "samples.tsv"
R_PARAMS_TSV = "r_params.tsv"
MATRIX = "matrix.mtx.gz"
GENES = "genes.txt"
CELLS = "cells.txt"
FITS = "fits.csv"
INPUT_CHECK = "input_check.csv"
# The port's convergence tolerance for parity runs (config ``gavish.tol: 1.0e-8``; the audit default
# is 1e-7). Looser exports are flagged by export_matrices.py and compare.py.
PARITY_TOL = 1e-8


def basis_file(rank: int) -> str:
    return f"W_rank{rank}.csv"


def load_config(path: str | Path | None) -> dict:
    """A YAML or JSON config as read by ``mpaudit p1`` (empty when no path is given)."""
    if not path:
        return {}
    text = Path(path).read_text()
    if str(path).endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(text) or {}
    return json.loads(text)


def read_manifest(export_dir: str | Path) -> dict:
    manifest = json.loads((Path(export_dir) / MANIFEST).read_text())
    if manifest.get("format") != FORMAT:
        raise ValueError(f"{export_dir}/{MANIFEST} is not a {FORMAT} manifest")
    if manifest.get("format_version") != FORMAT_VERSION:
        raise ValueError(f"Unsupported export format version {manifest.get('format_version')} "
                         f"(this script reads version {FORMAT_VERSION})")
    return manifest


def params_from_manifest(manifest: dict) -> GavishParams:
    """The exported :class:`GavishParams`; fields this mpaudit version does not know are dropped
    with a warning (the export came from different code)."""
    names = {f.name for f in dataclasses.fields(GavishParams)}
    given = dict(manifest["gavish"])
    unknown = sorted(set(given) - names)
    if unknown:
        warnings.warn(f"Ignoring GavishParams fields unknown to this mpaudit version: {unknown}", RuntimeWarning)
    return GavishParams.from_dict({k: v for k, v in given.items() if k in names})


def read_lines(path: str | Path) -> list[str]:
    return Path(path).read_text().split("\n")[:-1]


def write_lines(lines, path: str | Path) -> None:
    lines = [str(x) for x in lines]
    bad = [x for x in lines if "\n" in x or "\r" in x]
    if bad:
        raise ValueError(f"Names with line breaks cannot be exported: {bad[:3]}")
    Path(path).write_text("".join(f"{x}\n" for x in lines))


def read_matrix(export_dir: str | Path, entry: dict) -> np.ndarray:
    """The exported NMF input of one manifest sample entry (genes x cells, float64)."""
    from scipy.io import mmread

    M = mmread(Path(export_dir) / entry["dir"] / MATRIX)
    return np.asarray(M.toarray() if sp.issparse(M) else M, dtype=np.float64)


def write_basis(W: np.ndarray, genes, path: str | Path) -> None:
    """One basis in the R output layout, values at full precision (``%.17g``, as run_r_snmf.R)."""
    df = pd.DataFrame(np.asarray(W, dtype=np.float64), columns=[str(j + 1) for j in range(W.shape[1])])
    df.insert(0, "gene", [str(g) for g in genes])
    df.to_csv(path, index=False, float_format="%.17g")


def read_basis(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """``(W, genes)`` from a basis CSV; columns are put in program-index order."""
    df = pd.read_csv(path, dtype={"gene": str}, keep_default_na=False)  # gene symbols such as "NA" stay
    if "gene" not in df.columns:
        raise ValueError(f"{path}: first column must be 'gene'")
    cols = sorted((c for c in df.columns if c != "gene"), key=int)
    W = df[cols].to_numpy(dtype=np.float64)
    if not np.isfinite(W).all():
        raise ValueError(f"{path}: non-finite basis values")
    return W, df["gene"].to_numpy(dtype=object)


def available_ranks(base_dir: str | Path, manifest: dict) -> list[int]:
    """Ranks with a basis file for every manifest sample under ``base_dir``."""
    common = None
    for entry in manifest["samples"]:
        found = {int(p.stem.removeprefix("W_rank")) for p in (Path(base_dir) / entry["dir"]).glob("W_rank*.csv")}
        common = found if common is None else common & found
    return sorted(common or [])


def read_bases(base_dir: str | Path, manifest: dict, ranks) -> dict[str, SamplePrograms]:
    """Bases of every manifest sample at ``ranks`` as :class:`mpaudit.gavish.SamplePrograms`
    (manifest order). Every rank of a sample must share one gene universe."""
    out = {}
    for entry in manifest["samples"]:
        spg = None
        for r in ranks:
            path = Path(base_dir) / entry["dir"] / basis_file(r)
            if not path.exists():
                raise FileNotFoundError(f"Missing basis {path}")
            W, genes = read_basis(path)
            if W.shape[1] != r:
                raise ValueError(f"{path}: {W.shape[1]} columns for rank {r}")
            if spg is None:
                spg = SamplePrograms(sample=entry["sample"], universe=genes)
            elif not np.array_equal(spg.universe, genes):
                raise ValueError(f"{path}: gene rows differ from the sample's other ranks")
            spg.W[r] = W
        out[entry["sample"]] = spg
    return out


def write_bases(programs: dict[str, SamplePrograms], base_dir: str | Path, manifest: dict) -> None:
    """Write port bases in the R output layout (plus the port's fit diagnostics in ``fits.csv``)."""
    for entry in manifest["samples"]:
        spg = programs[entry["sample"]]
        d = Path(base_dir) / entry["dir"]
        d.mkdir(parents=True, exist_ok=True)
        rows = []
        for r in sorted(spg.W):
            write_basis(spg.W[r], spg.universe, d / basis_file(r))
            rows.append({"sample": entry["sample"], "rank": r, **spg.fits.get(r, {})})
        pd.DataFrame(rows).to_csv(d / FITS, index=False)


def md_table(df: pd.DataFrame, digits: int = 3) -> str:
    """Plain markdown table (no tabulate dependency)."""
    if df.empty:
        return "(none)\n"

    def fmt(v) -> str:
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else f"{v:.{digits}g}"
        return str(v).replace("|", "\\|")

    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * df.shape[1]]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines) + "\n"
