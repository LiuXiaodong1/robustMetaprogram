"""Sparse NMF with the sNMF/R objective (Kim & Park 2007), the method 3CA runs via R's NMF package.

Mirrors ``NMF::nmf(method = "snmf/r")`` (renozao/NMF, ``R/algorithms-snmf.R``):

* minimises ``||A - W H||_F^2 + eta^2 ||W||_F^2 + beta * sum_c (sum_j H[j, c])^2``
  with ``A`` genes x cells, ``W`` genes x k (programs), ``H`` k x cells; defaults
  ``beta = 0.01`` and ``eta = max(A)`` (R squares ``eta`` in the W step);
* initialises ``W`` uniformly at random with unit-norm columns and starts with the H step;
* restarts from a new random ``W`` whenever a row of ``H`` is all zero (at most 10 restarts);
  iterations are counted across restarts, as in R. A fit that reaches the restart limit stops
  with a warning, as in R, and is unusable (its last ``H`` has an all-zero row); ``return_info``
  flags it as ``failed``, and :mod:`mpaudit.gavish` counts such runs (``n_failed_fits``);
* reports R's ``.snmf.objective`` (``0.5 * (resid + eta ||W||^2 + beta sum colSums(H)^2)``,
  with ``eta`` unsquared), which is what picks the best of ``nrun`` restarts.

The optimiser differs: HALS coordinate descent instead of exact alternating NNLS. So does the
stopping rule. Both codes test convergence every 5 iterations, but this port stops when the
objective has dropped by less than a relative ``tol`` [1e-5] since the previous test, while R stops
when the row/column argmax assignments of ``W`` and ``H`` have not changed for 10 consecutive tests
and the KKT residual has fallen to 1e-4 of its first value. The iteration cap is ``max_iter``
[5,000]; R's default is 20,000. Fits that stop at the cap are counted (``return_info``, and
``GavishResult.n_capped_fits`` in :mod:`mpaudit.gavish`); :func:`mpaudit.pipeline.run_p1` reports
capped and restart-limit counts for the reference run, each decoy cohort (``decoy_runs.csv``) and
each instability perturbation (``instability_runs.csv``), with totals in ``run_config.json``
(``nmf_fits``) and ``summary.md``. Results are therefore not
bit-identical to R. Contrasts between perturbations inside ``mpaudit`` are like-for-like; check
parity with the R code before quoting absolute numbers.

The tolerance, not the cap, decides where fits stop. On the toy cohort (:mod:`mpaudit.sim`,
1,500 x 200 matrices, ranks 4-9) fits stopped after 135-305 iterations at ``tol=1e-5``, far from
the cap, yet their top-50 genes kept moving: continued on the same trajectory to ``tol=1e-8``,
a program kept a median of 44 of its top 50 genes (48 at ``tol=1e-6``, 49 at ``1e-7``, at about
2x and 3x the runtime). R's rule, applied to the same HALS iterates, stopped after a median of
about 900 iterations. :class:`mpaudit.gavish.GavishParams` therefore defaults to ``tol=1e-7``
(about 3x the runtime of 1e-5). Parity runs against the R code (``scripts/parity``) use
``tol=1e-8``, i.e. ``gavish.tol: 1.0e-8`` in the config (YAML reads ``1e-8`` without a decimal
point as a string; :class:`~mpaudit.gavish.GavishParams` converts it, but other tools may not).
``snmf_r`` itself keeps 1e-5 as its low-level default.
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
    max_iter: int = 5000,
    tol: float = 1e-5,
    check_every: int = 5,
    max_restarts: int = 10,
    rng: np.random.Generator | int | None = None,
    return_info: bool = False,
) -> tuple[np.ndarray, np.ndarray, float] | tuple[np.ndarray, np.ndarray, float, dict]:
    """Fit one sNMF/R factorisation. Returns ``(W, H, r_objective)``.

    With ``return_info`` a fourth element is returned: ``{"n_iter", "converged", "restarts",
    "failed"}``. ``n_iter`` counts iterations across restarts; ``converged`` is true only when the
    tolerance test stopped the fit (not ``max_iter`` or the restart limit); ``failed`` is true when
    the restart limit stopped it (``max_restarts`` all-zero rows of ``H``; the returned factors are
    then not a usable fit).
    """
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
    converged = failed = False
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
                failed = True
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
                converged = True
                break
            prev = obj
    resid = float(np.sum((A - W @ H) ** 2))
    r_obj = 0.5 * (resid + eta * float(np.sum(W * W)) + beta * float(np.sum(H.sum(axis=0) ** 2)))
    if return_info:
        return W, H, r_obj, {"n_iter": it, "converged": converged, "restarts": restarts, "failed": failed}
    return W, H, r_obj


def best_of_runs(
    A: np.ndarray,
    k: int,
    *,
    nrun: int,
    seeds: list[np.random.SeedSequence] | None = None,
    return_info: bool = False,
    **kwargs,
) -> np.ndarray | tuple[np.ndarray, list[dict]]:
    """Run ``nrun`` restarts and return ``W`` of the run with the lowest R objective
    (mirrors ``nmf(..., nrun=10)`` keeping the best fit; the first run wins ties).

    With ``return_info``, returns ``(W, infos)`` where ``infos`` holds every run's
    :func:`snmf_r` info dict, in run order.
    """
    seeds = seeds or [np.random.SeedSequence(i) for i in range(nrun)]
    best_W, best_obj, infos = None, np.inf, []
    for ss in seeds[:nrun]:
        W, _, obj, info = snmf_r(A, k, rng=np.random.default_rng(ss), return_info=True, **kwargs)
        infos.append(info)
        if obj < best_obj:
            best_W, best_obj = W, obj
    return (best_W, infos) if return_info else best_W
