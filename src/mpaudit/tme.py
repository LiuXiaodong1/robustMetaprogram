"""Pilot P1c / axis A1: does an MP co-vary with the same genes in co-resident TME cells?

Non-malignant cells from the same sample share ambient RNA, dissociation and ischaemia time and
RNA quality with the malignant cells, but not malignant programs. Two statistics per MP:

* ``partial_r2``: across samples, the share of variance in the per-sample malignant MP score
  explained by the same gene set's score in that sample's TME cells, beyond cancer type (and
  optionally study). TME scores are centred per cell type, so TME composition cancels. Only a
  positive association counts (``slope > 0``). It is compared with ``null_q95``, the same
  statistic for random gene sets matched gene-by-gene on expression: sample-level effects that
  move every gene set (compositional shifts when a strong program takes up reads, nuisance genes
  among the control genes) inflate the null too, so they cannot make one MP look TME-driven.
* ``kappa_ratio``: within samples, co-expression (``kappa``) of the MP genes in TME cell types
  divided by that in malignant cells. Near 0: malignant-specific co-expression; near 1 or
  above: the same genes co-vary in TME cells too (technical or niche-wide).

Scores are linear in log expression, so per-sample scores are computed from per-sample mean
expression profiles; the random-set null is therefore cheap even at atlas scale.

Neither statistic alone says "artefact": IFN response and hypoxia are expected to be
TME-shared biology. They feed the provenance typing, together with the protocol axis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .coexpr import kappa_with_null
from .data import Keys, get_counts, log2_cpm10, malignant_mask
from .scoring import control_genes


def _design(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    parts = [np.ones((len(df), 1))]
    for c in cols:
        d = pd.get_dummies(df[c].astype(str), drop_first=True, dtype=float)
        if d.shape[1]:
            parts.append(d.to_numpy())
    return np.hstack(parts)


def _fit(y: np.ndarray, x: np.ndarray, Z: np.ndarray) -> tuple[float, float]:
    """Partial R^2 of ``x`` given ``Z``, and the slope of ``x``."""
    def sse(M):
        beta, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ beta
        return float(r @ r), beta
    red, _ = sse(Z)
    full, beta = sse(np.hstack([Z, x[:, None]]))
    return ((red - full) / red if red > 0 else np.nan), float(beta[-1])


def _group_means(L: sp.csr_matrix, groups: list[np.ndarray]) -> np.ndarray:
    rows, cols, vals = [], [], []
    for i, idx in enumerate(groups):
        rows += [i] * idx.size
        cols += list(idx)
        vals += [1.0 / idx.size] * idx.size
    G = sp.csr_matrix((vals, (rows, cols)), shape=(len(groups), L.shape[0]))
    return np.asarray((G @ L).todense())


def _matched_random_set(idx: np.ndarray, bins: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    excl = np.zeros(bins.size, bool)
    excl[idx] = True
    out = []
    for g in idx:
        pool = np.flatnonzero((bins == bins[g]) & ~excl)
        if pool.size:
            out.append(rng.choice(pool))
    return np.unique(np.array(out, dtype=int))


def tme_covariation(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    min_cells: int = 20,
    covariates: tuple[str, ...] | None = None,
    n_boot: int = 200,
    n_null: int = 200,
    kappa_ratio: bool = True,
    max_kappa_cells: int = 2000,
    n_bins: int = 30,
    n_ctrl: int = 100,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(per_set, per_sample)`` tables for the A1 axis."""
    if keys.cell_type is None:
        raise ValueError("keys.cell_type is required for the TME axis")
    rng = np.random.default_rng(seed)
    counts = get_counts(adata, counts_layer)
    L = log2_cpm10(counts)
    mal = malignant_mask(adata, keys)
    ctype = adata.obs[keys.cell_type].astype(str).to_numpy()
    sample = adata.obs[keys.sample].astype(str).to_numpy()
    tme = ~mal & (adata.obs[keys.cell_type].notna().to_numpy())
    covariates = tuple(c for c in (covariates if covariates is not None else (keys.cancer_type,)) if c)

    # per-sample malignant and per-(sample, cell type) TME mean expression profiles
    samples, mal_groups, meta_rows = [], [], []
    for s in pd.unique(sample):
        r = np.flatnonzero(mal & (sample == s))
        if r.size >= min_cells:
            samples.append(s)
            mal_groups.append(r)
            meta_rows.append({"sample": s, **{c: str(adata.obs[c].iloc[r[0]]) for c in covariates}})
    meta = pd.DataFrame(meta_rows)
    tme_groups, tme_keys = [], []
    for i, s in enumerate(samples):
        t_rows = np.flatnonzero(tme & (sample == s))
        for ct in pd.unique(ctype[t_rows]):
            r = t_rows[ctype[t_rows] == ct]
            if r.size >= min_cells:
                tme_groups.append(r)
                tme_keys.append((i, ct))
    Mm = _group_means(L, mal_groups) if mal_groups else np.zeros((0, L.shape[1]))
    Mt = _group_means(L, tme_groups) if tme_groups else np.zeros((0, L.shape[1]))
    t_sample = np.array([k[0] for k in tme_keys], dtype=int)
    t_type = np.array([k[1] for k in tme_keys], dtype=object)

    ref_mean = np.asarray(L[np.flatnonzero(mal)].mean(axis=0)).ravel()
    order = np.argsort(ref_mean, kind="stable")
    bins = np.empty_like(order)
    bins[order] = np.floor(np.arange(order.size) * n_bins / order.size).astype(int)
    pos = {g: i for i, g in enumerate(map(str, adata.var_names))}
    Z_full = _design(meta, list(covariates)) if len(meta) else np.zeros((0, 1))

    def weights(idx: np.ndarray) -> np.ndarray:
        v = np.zeros(L.shape[1])
        v[idx] += 1.0 / idx.size
        ctrl = control_genes(idx, ref_mean, rng, n_bins, n_ctrl)
        if ctrl.size:
            v -= np.bincount(ctrl, minlength=L.shape[1]) / ctrl.size
        return v

    def xy(v: np.ndarray):
        y = Mm @ v
        x = np.full(len(samples), np.nan)
        if len(tme_keys):
            ts = Mt @ v
            centred = ts.copy()
            for ct in np.unique(t_type):  # centre each cell type across samples
                m = t_type == ct
                centred[m] -= ts[m].mean()
            sums = np.bincount(t_sample, weights=centred, minlength=len(samples))
            cnt = np.bincount(t_sample, minlength=len(samples))
            x = np.where(cnt > 0, sums / np.maximum(cnt, 1), np.nan)
        return y, x

    per_sample = meta.copy()
    out_rows = []
    for name, genes in gene_sets.items():
        idx = np.array(sorted({pos[g] for g in genes if g in pos}), dtype=int)
        rec = {"mp": name, "n_samples": 0, "partial_r2": np.nan, "slope": np.nan, "partial_r2_lo": np.nan,
               "partial_r2_hi": np.nan, "null_q95": np.nan, "spearman": np.nan}
        if idx.size >= 5 and len(samples):
            y, x = xy(weights(idx))
            per_sample[f"{name}::malignant"] = y
            per_sample[f"{name}::tme"] = x
            ok = np.isfinite(x) & np.isfinite(y)
            Z = Z_full[ok]
            if ok.sum() > Z.shape[1] + 2:
                yv, xv = y[ok], x[ok]
                rec["n_samples"] = int(ok.sum())
                rec["partial_r2"], rec["slope"] = _fit(yv, xv, Z)
                rec["spearman"] = float(pd.Series(yv).corr(pd.Series(xv), method="spearman"))
                boots = []
                for _ in range(n_boot):
                    b = rng.integers(0, ok.sum(), ok.sum())
                    Zb = _design(meta.loc[ok].iloc[b], list(covariates))
                    if np.linalg.matrix_rank(Zb) < Zb.shape[1] or b.size <= Zb.shape[1] + 2:
                        continue
                    boots.append(_fit(yv[b], xv[b], Zb)[0])
                boots = np.array([v for v in boots if np.isfinite(v)])
                if boots.size:
                    rec["partial_r2_lo"], rec["partial_r2_hi"] = (float(q) for q in np.quantile(boots, [0.05, 0.95]))
                nulls = []
                for _ in range(n_null):
                    ridx = _matched_random_set(idx, bins, rng)
                    if ridx.size < 5:
                        continue
                    ny, nx = xy(weights(ridx))
                    k = np.isfinite(nx) & np.isfinite(ny)
                    if k.sum() > Z_full[k].shape[1] + 2:
                        r2, slope = _fit(ny[k], nx[k], Z_full[k])
                        nulls.append(r2 if slope > 0 else 0.0)
                if nulls:
                    rec["null_q95"] = float(np.quantile(nulls, 0.95))
        if kappa_ratio:
            ratios = []
            for s in samples:
                m_rows = np.flatnonzero(mal & (sample == s))
                t_rows = np.flatnonzero(tme & (sample == s))
                if m_rows.size < min_cells or t_rows.size < min_cells:
                    continue
                if m_rows.size > max_kappa_cells:
                    m_rows = rng.choice(m_rows, max_kappa_cells, replace=False)
                km = kappa_with_null(counts[m_rows], idx, n_null=0)["kappa"]
                kt, wt = [], []
                for ct in pd.unique(ctype[t_rows]):
                    r = t_rows[ctype[t_rows] == ct]
                    if r.size >= min_cells:
                        if r.size > max_kappa_cells:
                            r = rng.choice(r, max_kappa_cells, replace=False)
                        v = kappa_with_null(counts[r], idx, n_null=0)["kappa"]
                        if np.isfinite(v):
                            kt.append(v)
                            wt.append(r.size)
                if np.isfinite(km) and km > 0 and kt:
                    ratios.append(float(np.average(kt, weights=wt)) / km)
            rec["kappa_ratio_median"] = float(np.median(ratios)) if ratios else np.nan
            rec["kappa_ratio_n"] = len(ratios)
        out_rows.append(rec)
    return pd.DataFrame(out_rows), per_sample
