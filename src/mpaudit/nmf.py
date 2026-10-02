"""Sparse NMF with the sNMF/R objective (Kim & Park 2007), the method 3CA runs via R's NMF package.

Minimises ``||A - W H||_F^2 + eta * ||W||_F^2 + beta * sum_c (sum_j H[j, c])^2`` with
``A`` genes x cells, ``W`` genes x k (the gene programs) and ``H`` k x cells.

The objective and defaults (``beta = 0.01``; ``eta = max(A)``) follow R's ``snmf/r``. The
optimiser differs: this uses HALS coordinate descent rather than alternating NNLS, so results
are not bit-identical to the R code. Use ``mpaudit`` results for perturbation contrasts within
the same engine; check parity with the original R code before quoting absolute numbers.
"""
from __future__ import annotations

import numpy as np


def snmf_r(
    A: np.ndarray,
    k: int,
    *,
    beta: float = 0.01,
    eta: float | None = None,
    max_iter: int = 1000,
    tol: float = 1e-5,
    check_every: int = 10,
    rng: np.random.Generator | int | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit one sNMF/R factorisation. Returns ``(W, H, objective)``."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or (A < 0).any():
        raise ValueError("A must be a 2-D non-negative matrix")
    m, n = A.shape
    if eta is None or eta < 0:
        eta = float(A.max()) if A.size else 1.0
    rng = np.random.default_rng(rng)
    scale = np.sqrt(max(float(A.mean()), 1e-12) / k)
    W = rng.random((m, k)) * 2.0 * scale
    H = rng.random((k, n)) * 2.0 * scale
    normA2 = float(np.einsum("ij,ij->", A, A))
    prev = np.inf
    obj = np.inf
    for it in range(max_iter):
        WtA = W.T @ A
        G = W.T @ W + beta
        for j in range(k):
            H[j] = np.maximum(H[j] + (WtA[j] - G[j] @ H) / G[j, j], 0.0)
        AHt = A @ H.T
        S = H @ H.T
        S[np.diag_indices(k)] += eta
        for j in range(k):
            W[:, j] = np.maximum(W[:, j] + (AHt[:, j] - W @ S[:, j]) / S[j, j], 0.0)
        if (it + 1) % check_every == 0 or it == max_iter - 1:
            HHt = S - eta * np.eye(k)
            resid = normA2 - 2.0 * float(np.sum(W * AHt)) + float(np.sum((W.T @ W) * HHt))
            obj = resid + eta * float(np.sum(W * W)) + beta * float(np.sum(H.sum(axis=0) ** 2))
            if np.isfinite(prev) and (prev - obj) <= tol * max(abs(prev), 1e-12):
                break
            prev = obj
    return W, H, obj


def best_of_runs(
    A: np.ndarray,
    k: int,
    *,
    nrun: int,
    seeds: list[np.random.SeedSequence] | None = None,
    **kwargs,
) -> np.ndarray:
    """Run ``nrun`` restarts and return ``W`` of the run with the lowest objective
    (mirrors ``nmf(..., nrun=10)`` keeping the best fit)."""
    seeds = seeds or [np.random.SeedSequence(i) for i in range(nrun)]
    best_W, best_obj = None, np.inf
    for ss in seeds[:nrun]:
        W, _, obj = snmf_r(A, k, rng=np.random.default_rng(ss), **kwargs)
        if obj < best_obj:
            best_W, best_obj = W, obj
    return best_W
