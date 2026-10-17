"""Export the per-sample NMF input matrices of the mpaudit port, for an R-versus-port parity run.

    python scripts/parity/export_matrices.py cohort.h5ad --config configs/p1_3ca_example.yaml \\
        --out parity/export --n 5            # or: --samples S1 S2 S3

For each selected sample this writes exactly the matrix mpaudit factorises,
:func:`mpaudit.gavish.prepare_matrix` applied to the sample's NMF cells (log2(CPM/10+1), top genes
by mean, centred, negatives set to zero, all-zero genes and cells dropped), as a genes x cells
Matrix Market file with ``genes.txt`` and ``cells.txt``; see ``parity_io`` for the layout. Values
are written as shortest round-trip decimals, so Python reads back the identical doubles; R's parser
can differ in the last bit, which ``compare.py`` tolerates (it checks R's reported sum and max).

``--n`` picks the samples with the most NMF cells: these take longest in R but give the most stable
programs. ``manifest.json`` records the keys, counts layer and :class:`mpaudit.gavish.GavishParams`
(from the ``mpaudit p1`` config), the seed R should use and, per sample, the matrix size, sum, max
and SHA-1, the recurrence unit and the study.

Cell selection re-implements the rule of :func:`mpaudit.gavish.compute_programs` from public
functions: malignant cells with >= ``min_genes_per_cell`` detected genes, in samples with >=
``min_cells`` of them (no subsampling). ``compare.py`` checks that the gene universe the port fits
equals the exported one, which catches any drift between the two. ``prepare_matrix`` does not say
which cells it dropped as all zero, so the exporter rebuilds the matrix over all selected cells and
asserts that the kept columns reproduce ``prepare_matrix``'s output before naming them.
"""
from __future__ import annotations

import argparse
import dataclasses
import gzip
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parity_io as pio  # noqa: E402

from mpaudit.data import Keys, get_counts, log2_cpm10, malignant_mask  # noqa: E402
from mpaudit.gavish import GavishParams, prepare_matrix  # noqa: E402

DEFAULT_R_SEED = 123456


def nmf_cells(adata, keys: Keys, params: GavishParams, counts: sp.csr_matrix) -> dict[str, np.ndarray]:
    """Row indices of the cells each sample contributes to NMF, for samples with >= ``min_cells``."""
    mal = malignant_mask(adata, keys)
    if params.min_genes_per_cell > 0:
        detected = np.diff((counts > 0).tocsr().indptr)
        mal = mal & (detected >= params.min_genes_per_cell)
    obs_sample = adata.obs[keys.sample].astype(str).to_numpy()
    cells = {}
    for s in pd.unique(obs_sample[mal]):
        idx = np.flatnonzero(mal & (obs_sample == s))
        if idx.size >= params.min_cells:
            cells[s] = idx
    return cells


def choose_samples(cells: dict[str, np.ndarray], samples: list[str] | None, n: int) -> list[str]:
    """The requested samples (all must be eligible), else the ``n`` with most NMF cells."""
    if samples:
        missing = [s for s in samples if s not in cells]
        if missing:
            raise ValueError(f"Not eligible for NMF (unknown, or fewer than min_cells NMF cells): {missing}")
        return list(dict.fromkeys(samples))
    return sorted(cells, key=lambda s: (-cells[s].size, s))[:n]


def kept_cells(counts: sp.csr_matrix, genes: np.ndarray, universe: np.ndarray, A: np.ndarray,
               params: GavishParams) -> np.ndarray:
    """Mask over the input cells for the columns of ``A = prepare_matrix(counts, genes, params)[0]``.

    Rebuilds the processed matrix over all input cells (centring uses every input cell, as in
    ``prepare_matrix``), keeps the columns with a positive value and asserts they reproduce ``A``.
    """
    col = pd.Index(genes).get_indexer(universe)
    L = log2_cpm10(counts)[:, col].toarray().T
    if params.center:
        L = L - L.mean(axis=1, keepdims=True)
        L[L < 0] = 0.0
    keep = np.ones(L.shape[1], dtype=bool) if A.shape[1] == L.shape[1] else L.any(axis=0)
    if L[:, keep].shape != A.shape or not np.allclose(L[:, keep], A, rtol=0.0, atol=1e-12):
        raise RuntimeError("Cannot map prepare_matrix's columns back to cells: prepare_matrix has changed in a "
                           "way kept_cells does not reproduce; update export_matrices.kept_cells")
    return keep


def safe_dir(i: int, sample: str) -> str:
    return f"{i + 1:02d}_{re.sub(r'[^A-Za-z0-9._-]+', '_', sample)[:60]}"


def export(adata, cfg: dict, out: Path, *, samples: list[str] | None = None, n: int = 5,
           r_seed: int = DEFAULT_R_SEED, sources: dict | None = None, log=print) -> dict:
    """Write the export for ``adata`` under the ``mpaudit p1`` config ``cfg``; returns the manifest.

    Only ``keys``, ``counts_layer``, ``recurrence_key`` and ``gavish`` are read; their defaults
    (:class:`mpaudit.data.Keys`, None, None, :class:`mpaudit.gavish.GavishParams`) are those of
    ``mpaudit p1``.
    """
    keys = Keys.from_dict(cfg.get("keys"))
    params = GavishParams.from_dict(cfg.get("gavish"))
    layer, rk = cfg.get("counts_layer"), cfg.get("recurrence_key")
    if not adata.var_names.is_unique:
        raise ValueError("var_names are not unique; call adata.var_names_make_unique() first")
    counts = get_counts(adata, layer)
    genes = np.asarray(adata.var_names.astype(str))
    obs_names = np.asarray(adata.obs_names.astype(str))
    cells = nmf_cells(adata, keys, params, counts)
    chosen = choose_samples(cells, samples, n)
    obs = adata.obs
    first = {s: int(cells[s][0]) for s in chosen}

    entries = []
    for i, s in enumerate(chosen):
        idx = cells[s]
        A, universe = prepare_matrix(counts[idx], genes, params)
        keep = kept_cells(counts[idx], genes, universe, A, params)
        d = out / "samples" / safe_dir(i, s)
        d.mkdir(parents=True, exist_ok=True)
        with gzip.open(d / pio.MATRIX, "wb", compresslevel=4) as fh:
            mmwrite(fh, sp.coo_matrix(A), comment=f" mpaudit parity export: sample {s}, genes x cells", field="real")
        pio.write_lines(universe, d / pio.GENES)
        pio.write_lines(obs_names[idx][keep], d / pio.CELLS)
        entries.append({
            "sample": s, "dir": f"samples/{d.name}", "n_genes": int(A.shape[0]), "n_cells": int(A.shape[1]),
            "n_cells_selected": int(idx.size), "n_cells_dropped_all_zero": int((~keep).sum()),
            "sum": float(A.sum()), "max": float(A.max()) if A.size else 0.0,
            "sha1": hashlib.sha1(np.ascontiguousarray(A, dtype=np.float64).tobytes()).hexdigest(),
            "unit": str(obs[rk].iloc[first[s]]) if rk else s,
            "study": str(obs[keys.study].iloc[first[s]]) if keys.study else None,
        })
        log(f"[export] {s}: {A.shape[0]} genes x {A.shape[1]} cells -> {d}")

    try:
        from mpaudit.gavish import code_version
        code = code_version()
    except ImportError:  # older mpaudit without a code hash
        code = None
    manifest = {
        "format": pio.FORMAT, "format_version": pio.FORMAT_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": sources or {}, "mpaudit_version": metadata.version("mpaudit"), "code_version": code,
        "keys": dataclasses.asdict(keys), "counts_layer": layer, "recurrence_key": rk,
        "gavish": {k: list(v) if isinstance(v, tuple) else v for k, v in dataclasses.asdict(params).items()},
        "r_seed": int(r_seed),
        "r_call": (f'nmf(x = A, rank = c({", ".join(map(str, params.ranks))}), method = "snmf/r", '
                   f"nrun = {params.nrun}, seed = {int(r_seed)})"),
        "samples": entries,
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / pio.MANIFEST).write_text(json.dumps(manifest, indent=2))
    pd.DataFrame([{k: e[k] for k in ("dir", "n_genes", "n_cells")} for e in entries]).to_csv(
        out / pio.SAMPLES_TSV, sep="\t", index=False)
    pd.DataFrame({"key": ["ranks", "nrun", "seed", "method"],
                  "value": [",".join(map(str, params.ranks)), str(params.nrun), str(int(r_seed)), "snmf/r"]}
                 ).to_csv(out / pio.R_PARAMS_TSV, sep="\t", index=False)
    return manifest


def main(argv=None) -> dict:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad")
    ap.add_argument("--config", help="the YAML or JSON config used for 'mpaudit p1'")
    ap.add_argument("--out", default="parity/export")
    pick = ap.add_mutually_exclusive_group()
    pick.add_argument("--samples", nargs="+", help="samples to export (default: the --n largest)")
    pick.add_argument("--n", type=int, default=5, help="export the N samples with most NMF cells [5]")
    ap.add_argument("--seed", type=int, default=DEFAULT_R_SEED, help=f"seed for R's nmf() [{DEFAULT_R_SEED}]")
    args = ap.parse_args(argv)

    import anndata as ad

    cfg = pio.load_config(args.config)
    sources = {"h5ad": str(Path(args.h5ad).resolve()),
               "config": str(Path(args.config).resolve()) if args.config else None}
    manifest = export(ad.read_h5ad(args.h5ad), cfg, Path(args.out), samples=args.samples, n=args.n,
                      r_seed=args.seed, sources=sources)
    print(f"[export] wrote {len(manifest['samples'])} samples to {args.out}; next: "
          f"Rscript scripts/parity/run_r_snmf.R {args.out} --out <r_out>")
    return manifest


if __name__ == "__main__":
    main()
