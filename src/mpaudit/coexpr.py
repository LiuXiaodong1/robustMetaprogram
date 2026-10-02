"""Depth-robust co-expression of a gene set (``kappa``).

Latent gene-gene correlations are estimated with the moment estimator behind CS-CORE
(Su et al. 2023), with its iteratively re-weighted step: counts ``y_cg | x_cg ~ Poisson(s_c x_cg)`` with
``s_c`` the cell's total UMI count, ``E[x_g] = mu_g``, so that

    E[(y_cg - s_c mu_g)^2] = s_c^2 sigma_g^2 + s_c mu_g,
    E[(y_cg - s_c mu_g)(y_cg' - s_c mu_g')] = s_c^2 sigma_gg'.

``kappa`` is the mean off-diagonal latent correlation of a gene set, compared with a null band
from random gene sets matched on mean expression and detection rate.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def latent_corr(Y: np.ndarray, s: np.ndarray, n_iter: int = 2) -> np.ndarray:
    """Latent correlation matrix of the columns of ``Y`` (cells x genes counts).

    Starts from the unweighted moment estimates, then re-weights each cell by the inverse of
    its modelled variance ``v_cg = s_c mu_g + s_c^2 sigma_g^2`` (CS-CORE's IRLS step), which
    keeps a few very deep cells from dominating when depth varies widely.
    """
    Y = np.asarray(Y, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)
    mu = Y.sum(axis=0) / max(s.sum(), 1e-12)
    R = Y - s[:, None] * mu[None, :]
    s2 = s**2
    m = s[:, None] * mu[None, :]
    D = R**2 - m
    var = (D * s2[:, None]).sum(axis=0) / float((s2**2).sum())
    cov = (R * s2[:, None]).T @ R / float((s2**2).sum())
    for _ in range(n_iter):
        sig = np.clip(var, 1e-12 * np.maximum(mu, 1e-12) ** 2, None)
        V = m + s2[:, None] * sig[None, :]
        V = np.maximum(V, 1e-12)
        var = ((D * s2[:, None]) / V**2).sum(axis=0) / ((s2**2)[:, None] / V**2).sum(axis=0)
        Rv = R / V
        cov = (Rv * s2[:, None]).T @ Rv / (((s2**2)[:, None] / V).T @ (1.0 / V))
    np.fill_diagonal(cov, var)
    with np.errstate(invalid="ignore", divide="ignore"):
        sd = np.sqrt(np.where(var > 0, var, np.nan))
        C = cov / np.outer(sd, sd)
    return np.clip(C, -1.0, 1.0)


def kappa(Y: np.ndarray, s: np.ndarray) -> float:
    """Mean off-diagonal latent correlation (NaN-aware)."""
    if Y.shape[1] < 2:
        return np.nan
    C = latent_corr(Y, s)
    iu = np.triu_indices(C.shape[0], k=1)
    v = C[iu]
    return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan


def matched_random_sets(
    target: np.ndarray,
    mean_expr: np.ndarray,
    detect: np.ndarray,
    n_sets: int,
    rng: np.random.Generator,
    n_mean_bins: int = 20,
    n_detect_bins: int = 3,
    candidates: np.ndarray | None = None,
) -> list[np.ndarray]:
    """Random gene sets matched gene-by-gene on mean-expression and detection-rate bins."""
    G = mean_expr.size
    cand = np.ones(G, bool) if candidates is None else candidates.copy()
    cand[target] = False
    mb = np.searchsorted(np.quantile(mean_expr, np.linspace(0, 1, n_mean_bins + 1)[1:-1]), mean_expr)
    db = np.searchsorted(np.quantile(detect, np.linspace(0, 1, n_detect_bins + 1)[1:-1]), detect)
    cell = mb * n_detect_bins + db
    pools = {c: np.flatnonzero(cand & (cell == c)) for c in np.unique(cell[target])}
    fallback = {c: np.flatnonzero(cand & (mb == c // n_detect_bins)) for c in pools}
    sets = []
    for _ in range(n_sets):
        picks = []
        for t in target:
            c = cell[t]
            pool = pools[c] if pools[c].size else fallback[c]
            if pool.size == 0:
                pool = np.flatnonzero(cand)
            picks.append(rng.choice(pool))
        sets.append(np.unique(np.array(picks)))
    return sets


def kappa_with_null(
    counts: sp.spmatrix,
    gene_idx: np.ndarray,
    *,
    n_null: int = 50,
    rng: np.random.Generator | None = None,
    min_detect: float = 0.05,
    min_genes: int = 5,
) -> dict:
    """``kappa`` of a gene set in one group of cells, with a matched-random null band.

    Only set genes detected in at least ``min_detect`` of the cells are used.
    """
    rng = rng or np.random.default_rng(0)
    counts = sp.csr_matrix(counts)
    s = np.asarray(counts.sum(axis=1)).ravel()
    keep = s > 0
    counts, s = counts[keep], s[keep]
    n = counts.shape[0]
    out = {"n_cells": int(n), "n_genes": 0, "kappa": np.nan, "null_mean": np.nan, "null_q95": np.nan}
    if n < 3:
        return out
    detect = np.asarray((counts > 0).sum(axis=0)).ravel() / n
    mean_expr = np.asarray(counts.sum(axis=0)).ravel() / max(s.sum(), 1e-12)
    eligible = detect >= min_detect
    target = np.array([g for g in gene_idx if eligible[g]], dtype=int)
    out["n_genes"] = int(target.size)
    if target.size < min_genes:
        return out
    out["kappa"] = kappa(counts[:, target].toarray(), s)
    if n_null > 0:
        nulls = [kappa(counts[:, r].toarray(), s) for r in
                 matched_random_sets(target, mean_expr, detect, n_null, rng, candidates=eligible)]
        nulls = np.array([v for v in nulls if np.isfinite(v)])
        if nulls.size:
            out["null_mean"] = float(nulls.mean())
            out["null_q95"] = float(np.quantile(nulls, 0.95))
    return out
