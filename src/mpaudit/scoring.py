"""Gene-set scores with expression-matched control genes (Tirosh et al. 2016 style).

Scores are *not* centred per sample, so they can be compared across samples and compartments:
``score_c = mean_{g in set} L_cg - mean_{g in controls} L_cg`` with ``L = log2(CPM/10+1)``.
Control genes are drawn once, from expression bins computed on a reference group of cells
(by default the malignant cells), and the same controls are used for every cell.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp


def control_genes(set_idx: np.ndarray, ref_mean: np.ndarray, rng: np.random.Generator,
                  n_bins: int = 30, n_ctrl: int = 100) -> np.ndarray:
    order = np.argsort(ref_mean, kind="stable")
    bins = np.empty_like(order)
    bins[order] = np.floor(np.arange(order.size) * n_bins / order.size).astype(int)
    in_set = np.zeros(ref_mean.size, bool)
    in_set[set_idx] = True
    ctrl = []
    for g in set_idx:
        pool = np.flatnonzero((bins == bins[g]) & ~in_set)
        if pool.size:
            ctrl.append(rng.choice(pool, size=min(n_ctrl, pool.size), replace=False))
    return np.concatenate(ctrl) if ctrl else np.array([], dtype=int)


def score_gene_sets(
    L: sp.csr_matrix,
    var_names,
    gene_sets: dict[str, list[str]],
    ref_rows: np.ndarray,
    *,
    n_bins: int = 30,
    n_ctrl: int = 100,
    seed: int = 0,
    min_genes: int = 5,
) -> pd.DataFrame:
    """Score every row of ``L`` (log-normalised, cells x genes) for every gene set."""
    rng = np.random.default_rng(seed)
    L = sp.csr_matrix(L)
    pos = {g: i for i, g in enumerate(map(str, var_names))}
    ref_mean = np.asarray(L[ref_rows].mean(axis=0)).ravel()
    scores = {}
    for name, genes in gene_sets.items():
        idx = np.array([pos[g] for g in genes if g in pos], dtype=int)
        if idx.size < min_genes:
            scores[name] = np.full(L.shape[0], np.nan)
            continue
        ctrl = control_genes(idx, ref_mean, rng, n_bins, n_ctrl)
        set_mean = np.asarray(L[:, idx].mean(axis=1)).ravel()
        if ctrl.size:
            w = np.bincount(ctrl, minlength=L.shape[1]).astype(float) / ctrl.size
            ctrl_mean = L @ w
        else:
            ctrl_mean = 0.0
        scores[name] = set_mean - ctrl_mean
    return pd.DataFrame(scores)
