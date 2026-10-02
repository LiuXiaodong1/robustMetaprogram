"""Pilot P1b: nuisance-matched decoys ("how much can nuisance manufacture?").

A decoy keeps, per sample and per gene, the count distribution conditional on declared
nuisance strata (sequencing-depth bins plus any obs columns such as clone, doublet score or an
ambient-contamination estimate) and destroys all other gene-gene co-expression. Running the
unmodified MP pipeline on decoys shows which MPs nuisance alone can reproduce.

Two methods:

* ``"residual"`` (default): stratified conditional residual permutation (SCRP). Per gene a
  negative-binomial model with mean ``s_c * p_{g,stratum}`` (``s_c`` = size factor) gives
  randomised-quantile residuals; residuals are permuted within strata and mapped back through
  each receiving cell's own distribution, so the continuous depth effect is kept.
* ``"permute"``: raw counts permuted within strata (exact stratum marginals, coarser depth).

Strata leak a little real signal: per-stratum gene means carry the stratum's share of every
program, so decoy co-expression grows roughly with (number of strata / cells per sample). Depth
is therefore modelled continuously through ``s_c`` by default (``n_depth_bins=1``), and binned
strata are capped so each has about ``min_stratum_cells`` cells. Categorical strata (clone,
doublet class) are used as given; report their sizes.

This is an audit of pipelines that do not model nuisance. It carries no FDR claim.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .compare import jaccard_matrix
from .data import Keys, get_counts, malignant_mask, sample_map
from .gavish import GavishParams, NMFCache, assemble, compute_programs


def _strata(obs: pd.DataFrame, totals: np.ndarray, strata: tuple[str, ...], n_depth_bins: int, n_cont_bins: int,
            min_stratum_cells: int = 50) -> np.ndarray:
    cap = max(1, len(totals) // max(min_stratum_cells, 1))
    n_depth_bins, n_cont_bins = min(n_depth_bins, cap), min(n_cont_bins, cap)
    labels = []
    if n_depth_bins > 1:
        labels.append(pd.qcut(np.log1p(totals), n_depth_bins, labels=False, duplicates="drop").astype(str))
    for col in strata:
        v = obs[col]
        if pd.api.types.is_numeric_dtype(v) and v.nunique() > 10:
            if n_cont_bins < 2:
                continue
            labels.append(pd.qcut(v.rank(method="first"), n_cont_bins, labels=False, duplicates="drop").astype(str).to_numpy())
        else:
            labels.append(v.astype(str).to_numpy())
    if not labels:
        return np.zeros(len(totals), dtype=int)
    joined = pd.Series(["|".join(t) for t in zip(*labels)])
    return pd.factorize(joined)[0]


def _permute_within(U: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = U.copy()
    for h in np.unique(strata):
        rows = np.flatnonzero(strata == h)
        if rows.size < 2:
            continue
        perm = np.argsort(rng.random((rows.size, U.shape[1])), axis=0)
        out[rows] = np.take_along_axis(U[rows], perm, axis=0)
    return out


def _nb_dispersion(Y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    num = (mu**2).sum(axis=0)
    den = ((Y - mu) ** 2 - mu).sum(axis=0)
    theta = np.where(den > 0, num / np.maximum(den, 1e-12), np.inf)
    return np.clip(theta, 1e-2, np.inf)


def _scrp_block(Y: np.ndarray, s: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    n, g = Y.shape
    mu = np.empty_like(Y)
    for h in np.unique(strata):
        rows = strata == h
        p = Y[rows].sum(axis=0) / max(s[rows].sum(), 1e-12)
        mu[rows] = s[rows, None] * p[None, :]
    theta = _nb_dispersion(Y, mu)
    pois = ~np.isfinite(theta) | (theta > 1e6)
    th = np.where(pois, 1.0, theta)[None, :]
    prob = th / (th + mu)

    def cdf(x):
        out = np.where(pois[None, :], stats.poisson.cdf(x, mu), stats.nbinom.cdf(x, th, prob))
        return np.where(x < 0, 0.0, out)

    lo, hi = cdf(Y - 1), cdf(Y)
    U = lo + rng.random(Y.shape) * (hi - lo)
    U = np.clip(_permute_within(U, strata, rng), 0.0, 1.0 - 1e-12)
    out = np.where(pois[None, :], stats.poisson.ppf(U, mu), stats.nbinom.ppf(U, th, prob))
    out = np.nan_to_num(out, nan=0.0, posinf=0.0)
    return np.maximum(out, 0.0)


def make_decoy(
    adata,
    keys: Keys,
    *,
    counts_layer: str | None = None,
    strata: tuple[str, ...] = (),
    n_depth_bins: int = 1,
    n_cont_bins: int = 3,
    min_stratum_cells: int = 50,
    method: str = "residual",
    seed: int = 0,
    gene_chunk: int = 2000,
):
    """Return a copy of ``adata`` whose malignant-cell counts are replaced by a decoy."""
    if method not in ("residual", "permute"):
        raise ValueError("method must be 'residual' or 'permute'")
    rng = np.random.default_rng(seed)
    X = get_counts(adata, counts_layer)
    mal = malignant_mask(adata, keys)
    samples = adata.obs[keys.sample].astype(str).to_numpy()
    keep_rows = np.ones(X.shape[0])
    add_r, add_c, add_v = [], [], []
    for s in pd.unique(samples[mal]):
        idx = np.flatnonzero(mal & (samples == s))
        if idx.size < 2:
            continue
        Xs = X[idx]
        totals = np.asarray(Xs.sum(axis=1)).ravel()
        st = _strata(adata.obs.iloc[idx], totals, strata, n_depth_bins, n_cont_bins, min_stratum_cells)
        sf = totals / max(np.median(totals), 1e-12)
        sf[sf == 0] = 1e-6
        expressed = np.flatnonzero(np.asarray(Xs.sum(axis=0)).ravel() > 0)
        keep_rows[idx] = 0.0
        for start in range(0, expressed.size, gene_chunk):
            cols = expressed[start:start + gene_chunk]
            Y = Xs[:, cols].toarray()
            block = _permute_within(Y, st, rng) if method == "permute" else _scrp_block(Y, sf, st, rng)
            coo = sp.coo_matrix(block)
            add_r.append(idx[coo.row]); add_c.append(cols[coo.col]); add_v.append(coo.data)
    X_new = (sp.diags(keep_rows) @ X).tocsr()
    if add_r:
        X_new = X_new + sp.csr_matrix((np.concatenate(add_v), (np.concatenate(add_r), np.concatenate(add_c))),
                                      shape=X.shape)
    out = adata.copy()
    X_new = X_new.tocsr()
    if counts_layer:
        out.layers[counts_layer] = X_new
    else:
        out.X = X_new
    return out


def decoy_recruitment(
    adata,
    keys: Keys,
    params: GavishParams,
    real_mps: dict[str, list[str]],
    *,
    n_decoys: int = 3,
    strata: tuple[str, ...] = (),
    n_depth_bins: int = 1,
    method: str = "residual",
    counts_layer: str | None = None,
    cache: NMFCache | None = None,
    recurrence_key: str | None = None,
    min_jaccard: float = 0.3,
    seed: int = 1000,
    n_jobs: int = 1,
    log=print,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, list[str]]]]:
    """Run the unmodified pipeline on ``n_decoys`` decoy cohorts.

    Returns ``(per_mp, per_decoy, decoy_mps)``. ``per_mp.recruitment`` is the fraction of decoy
    cohorts containing an MP with Jaccard >= ``min_jaccard`` to the real MP.
    """
    cache = cache or NMFCache()
    unit = sample_map(adata, keys, recurrence_key) if recurrence_key else None
    study = sample_map(adata, keys, keys.study) if keys.study else None
    best = {m: [] for m in real_mps}
    per_decoy, decoy_mps = [], {}
    for d in range(n_decoys):
        tag = f"decoy-{method}-{'+'.join(strata)}-d{n_depth_bins}-s{seed + d}"
        dec = make_decoy(adata, keys, counts_layer=counts_layer, strata=strata, n_depth_bins=n_depth_bins,
                         method=method, seed=seed + d)
        spg = compute_programs(dec, keys, params, counts_layer=counts_layer, cache=cache, data_tag=tag, n_jobs=n_jobs)
        res = assemble(spg, params, unit_of_sample=unit, study_of_sample=study)
        sets = res.mp_sets()
        decoy_mps[tag] = sets
        per_decoy.append({"decoy": tag, "n_robust": len(res.robust), "n_mps": len(sets)})
        J = jaccard_matrix(real_mps, sets) if sets else np.zeros((len(real_mps), 0))
        for i, m in enumerate(real_mps):
            best[m].append(float(J[i].max()) if J.shape[1] else 0.0)
        log(f"[decoy] {tag}: {len(res.robust)} robust programs, {len(sets)} decoy MPs")
    rows = [{"mp": m, "recruitment": float(np.mean(np.array(v) >= min_jaccard)) if v else np.nan,
             "mean_best_jaccard": float(np.mean(v)) if v else np.nan,
             "max_best_jaccard": float(np.max(v)) if v else np.nan} for m, v in best.items()]
    return pd.DataFrame(rows), pd.DataFrame(per_decoy), decoy_mps
