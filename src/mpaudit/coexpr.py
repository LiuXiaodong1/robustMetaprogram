"""Co-expression of a gene set (``kappa``): the mean latent correlation, corrected for Poisson sampling.

Latent gene-gene correlations are estimated with the moment estimator behind CS-CORE
(Su et al. 2023), with its iteratively re-weighted step: counts ``y_cg | x_cg ~ Poisson(s_c x_cg)`` with
``s_c`` the cell's total UMI count, ``E[x_g] = mu_g``, so that

    E[(y_cg - s_c mu_g)^2] = s_c^2 sigma_g^2 + s_c mu_g,
    E[(y_cg - s_c mu_g)(y_cg' - s_c mu_g')] = s_c^2 sigma_gg'.

``kappa`` is the mean off-diagonal latent correlation of a gene set, compared with a null band
from random gene sets matched on mean expression and detection rate. The estimand does not depend
on depth or on a gene-wise capture efficiency; the estimate does: with fewer counts per cell more
pairs get a non-positive variance estimate (NaN, left out) and the rest shrink, so ``kappa`` falls,
more for weakly co-expressed sets (toy: x 0.75-0.92 for planted programs and x 0.48-0.60 for random
40-gene sets at a quarter of the depth, :mod:`mpaudit.protocol`). Compare ``kappa`` across samples
only at equal depth and on equal numbers of cells, as :mod:`mpaudit.protocol` does.
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


def _gene_moments(Y: np.ndarray, s: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per gene (column of ``Y``), the quantities of :func:`latent_corr` (``n_iter=2``, as in :func:`kappa`)
    that involve that gene alone: ``mu``, the clipped variance ``sig`` behind the weights of the last IRLS
    step, and the final variance ``var``."""
    Y = np.asarray(Y, dtype=np.float64)
    mu = Y.sum(axis=0) / max(s.sum(), 1e-12)
    R = Y - s[:, None] * mu[None, :]
    s2 = s**2
    m = s[:, None] * mu[None, :]
    D = R**2 - m
    var = (D * s2[:, None]).sum(axis=0) / float((s2**2).sum())
    sig = var
    for _ in range(2):
        sig = np.clip(var, 1e-12 * np.maximum(mu, 1e-12) ** 2, None)
        V = m + s2[:, None] * sig[None, :]
        V = np.maximum(V, 1e-12)
        var = ((D * s2[:, None]) / V**2).sum(axis=0) / ((s2**2)[:, None] / V**2).sum(axis=0)
    return mu, sig, var


def resample_kappas(Y: np.ndarray, s: np.ndarray, rows: np.ndarray, max_elements: int = 2_000_000) -> np.ndarray:
    """``kappa(Y[r], s[r])`` for every row ``r`` of ``rows`` (resamples of cells, e.g. a bootstrap), computed in
    batches of resamples.

    The operations of :func:`latent_corr` (``n_iter=2``) and :func:`kappa`, applied to a stack of resamples
    at once (``max_elements`` bounds the size of a stack: resamples x cells x genes), so the values are
    those of :func:`kappa` resample by resample, to rounding (``tests/test_protocol.py``).
    """
    Y = np.asarray(Y, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)
    rows = np.atleast_2d(np.asarray(rows, dtype=int))
    n_res, k = rows.shape[0], Y.shape[1]
    out = np.full(n_res, np.nan)
    if k < 2 or n_res == 0:
        return out
    iu = np.triu_indices(k, k=1)
    diag = np.arange(k)
    step = max(1, max_elements // max(rows.shape[1] * k, 1))
    for a in range(0, n_res, step):
        b = rows[a:a + step]
        Yb, sb = Y[b], s[b]  # (resamples, cells, genes), (resamples, cells)
        mu = Yb.sum(axis=1) / np.maximum(sb.sum(axis=1), 1e-12)[:, None]
        R = Yb - sb[:, :, None] * mu[:, None, :]
        s2 = sb**2
        m = sb[:, :, None] * mu[:, None, :]
        D = R**2 - m
        var = (D * s2[:, :, None]).sum(axis=1) / (s2**2).sum(axis=1)[:, None]
        for _ in range(2):
            sig = np.clip(var, 1e-12 * np.maximum(mu, 1e-12) ** 2, None)
            V = m + s2[:, :, None] * sig[:, None, :]
            V = np.maximum(V, 1e-12)
            var = ((D * s2[:, :, None]) / V**2).sum(axis=1) / ((s2**2)[:, :, None] / V**2).sum(axis=1)
        Rv = R / V
        cov = (np.matmul((Rv * s2[:, :, None]).transpose(0, 2, 1), Rv)
               / np.matmul(((s2**2)[:, :, None] / V).transpose(0, 2, 1), 1.0 / V))
        cov[:, diag, diag] = var
        with np.errstate(invalid="ignore", divide="ignore"):
            sd = np.sqrt(np.where(var > 0, var, np.nan))
            C = np.clip(cov / (sd[:, :, None] * sd[:, None, :]), -1.0, 1.0)
        v = C[:, iu[0], iu[1]]
        fin = np.isfinite(v)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[a:a + len(b)] = np.where(fin.any(axis=1), np.where(fin, v, 0.0).sum(axis=1) / fin.sum(axis=1), np.nan)
    return out


def set_kappas(counts, s: np.ndarray, sets: list[np.ndarray], chunk: int = 1024,
               max_dense: int = 10_000_000) -> np.ndarray:
    """``kappa(counts[:, r], s)`` for every gene set ``r`` in ``sets`` (column indices), computed together.

    In :func:`latent_corr` everything but the covariance of a pair of genes is computed gene by gene,
    and only the last IRLS step's covariance reaches the result. So the per-gene moments
    (:func:`_gene_moments`) are computed once over the union of the sets' genes (``chunk`` genes at a
    time, fewer when cells x ``chunk`` exceeds ``max_dense // 8``, since that step holds about eight
    arrays of that size), and per set only the covariances of its own genes. The union's counts are
    read once into a dense array when cells x union genes is at most ``max_dense`` (80 MB at the
    default), else column by column for each set. Memory therefore stays at about two arrays of
    ``max_dense`` float64 values plus the per-set arrays, whatever the number of cells. Same values as
    :func:`kappa` of each set (the same operations in the same order; equal to rounding at worst, see
    ``tests/test_decoy_arms.py``).
    ``counts``: cells x genes, dense or sparse (CSC slices fastest); ``s``: the cells' size factors.
    """
    s = np.asarray(s, dtype=np.float64)
    out = np.full(len(sets), np.nan)
    sets = [np.asarray(r, dtype=int) for r in sets]
    union = np.unique(np.concatenate(sets)) if sets else np.array([], int)
    if not union.size:
        return out
    mu, sig, var = (np.empty(union.size) for _ in range(3))

    def dense(cols):
        Y = counts[:, cols]
        return np.asarray(Y.toarray() if sp.issparse(Y) else Y, dtype=np.float64)

    s2 = s**2

    def weights(Y, q):  # the last IRLS step's arrays of latent_corr for genes at union positions q
        R = Y - s[:, None] * mu[q][None, :]
        m = s[:, None] * mu[q][None, :]
        V = m + s2[:, None] * sig[q][None, :]
        V = np.maximum(V, 1e-12)
        Rv = R / V
        return Rv * s2[:, None], Rv, (s2**2)[:, None] / V, 1.0 / V

    chunk = max(1, min(chunk, max_dense // (8 * max(counts.shape[0], 1))))
    for a in range(0, union.size, chunk):
        b = slice(a, a + chunk)
        mu[b], sig[b], var[b] = _gene_moments(dense(union[b]), s)
    pos = np.full(int(union[-1]) + 1, -1)
    pos[union] = np.arange(union.size)
    Yu = dense(union) if counts.shape[0] * union.size <= max_dense else None
    triu = {}
    for i, r in enumerate(sets):
        if r.size < 2:
            continue
        p = pos[r]
        Rvs2, Rv, s4V, iV = weights(Yu[:, p] if Yu is not None else dense(r), p)
        cov = Rvs2.T @ Rv / (s4V.T @ iV)
        np.fill_diagonal(cov, var[p])
        with np.errstate(invalid="ignore", divide="ignore"):
            sd = np.sqrt(np.where(var[p] > 0, var[p], np.nan))
            C = np.clip(cov / np.outer(sd, sd), -1.0, 1.0)
        if r.size not in triu:
            triu[r.size] = np.triu_indices(r.size, k=1)
        v = C[triu[r.size]]
        nan = np.isnan(v)
        if not nan.all():  # np.nanmean, written out
            out[i] = float(np.sum(np.where(nan, 0.0, v)) / np.sum(~nan))
    return out


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
