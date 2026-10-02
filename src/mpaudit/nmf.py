"""Sparse NMF with the sNMF/R objective (Kim & Park 2007), the method 3CA runs via R's NMF package.

Mirrors ``NMF::nmf(method = "snmf/r")`` (renozao/NMF, ``R/algorithms-snmf.R``):

* minimises ``||A - W H||_F^2 + eta^2 ||W||_F^2 + beta * sum_c (sum_j H[j, c])^2``
  with ``A`` genes x cells, ``W`` genes x k (programs), ``H`` k x cells; defaults
  ``beta = 0.01`` and ``eta = max(A)`` (R squares ``eta`` in the W step);
* initialises ``W`` uniformly at random with unit-norm columns and starts with the H step;
* restarts from a new random ``W`` whenever a row of ``H`` is all zero (at most 10 restarts);
* reports R's ``.snmf.objective`` (``0.5 * (resid + eta ||W||^2 + beta sum colSums(H)^2)``,
  with ``eta`` unsquared), which is what picks the best of ``nrun`` restarts.

The optimiser differs: HALS coordinate descent instead of exact alternating NNLS, and a relative
objective-change stopping rule instead of R's assignment-stability plus KKT test. Results are
therefore not bit-identical to R. Contrasts between perturbations inside ``mpaudit`` are
like-for-like; check parity with the R code before quoting absolute numbers.
"""
from __future__ import annotations

import warnings

import numpy as np


def _init_W(m: int, k: int, rng: np.random.Generator) -> np.ndarray:
    W = rng.random((m, k))
    return W / np.sqrt((W**2).sum(axis=0, keepdims=True))


def snmf_r(
    A: np.ndarray,
    k: int,
    *,
    beta: float = 0.01,
    eta: float | None = None,
    max_iter: int = 1000,
    tol: float = 1e-5,
    check_every: int = 10,
    max_restarts: int = 10,
    rng: np.random.Generator | int | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit one sNMF/R factorisation. Returns ``(W, H, r_objective)``."""
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or (A < 0).any():
        raise ValueError("A must be a 2-D non-negative matrix")
    m, n = A.shape
    if eta is None or eta < 0:
        eta = float(A.max()) if A.size else 1.0
    eta2 = eta**2
    rng = np.random.default_rng(rng)
    normA2 = float(np.einsum("ij,ij->", A, A))

    W = _init_W(m, k, rng)
    H = np.zeros((k, n))
    prev = np.inf
    restarts = 0
    it = 0
    while it < max_iter:
        it += 1
        WtA = W.T @ A
        G = W.T @ W + beta
        for j in range(k):
            H[j] = np.maximum(H[j] + (WtA[j] - G[j] @ H) / G[j, j], 0.0)
        if (H.sum(axis=1) == 0).any():
            restarts += 1
            if restarts >= max_restarts:
                warnings.warn("snmf_r: too many restarts (zero rows in H); stopping", RuntimeWarning)
                break
            W, H, prev = _init_W(m, k, rng), np.zeros((k, n)), np.inf
            continue
        AHt = A @ H.T
        S = H @ H.T
        S[np.diag_indices(k)] += eta2
        for j in range(k):
            W[:, j] = np.maximum(W[:, j] + (AHt[:, j] - W @ S[:, j]) / S[j, j], 0.0)
        if it % check_every == 0:
            HHt = S - eta2 * np.eye(k)
            resid = normA2 - 2.0 * float(np.sum(W * AHt)) + float(np.sum((W.T @ W) * HHt))
            obj = resid + eta2 * float(np.sum(W * W)) + beta * float(np.sum(H.sum(axis=0) ** 2))
            if np.isfinite(prev) and (prev - obj) <= tol * max(abs(prev), 1e-12):
                break
            prev = obj
    resid = float(np.sum((A - W @ H) ** 2))
    r_obj = 0.5 * (resid + eta * float(np.sum(W * W)) + beta * float(np.sum(H.sum(axis=0) ** 2)))
    return W, H, r_obj


def best_of_runs(
    A: np.ndarray,
    k: int,
    *,
    nrun: int,
    seeds: list[np.random.SeedSequence] | None = None,
    **kwargs,
) -> np.ndarray:
    """Run ``nrun`` restarts and return ``W`` of the run with the lowest R objective
    (mirrors ``nmf(..., nrun=10)`` keeping the best fit)."""
    seeds = seeds or [np.random.SeedSequence(i) for i in range(nrun)]
    best_W, best_obj = None, np.inf
    for ss in seeds[:nrun]:
        W, _, obj = snmf_r(A, k, rng=np.random.default_rng(ss), **kwargs)
        if obj < best_obj:
            best_W, best_obj = W, obj
    return best_W
