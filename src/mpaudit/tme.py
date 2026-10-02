"""Pilot P1c / axis A1: does an MP co-vary with the same genes in co-resident TME cells?

Non-malignant cells from the same sample share ambient RNA, dissociation and ischaemia time and
RNA quality with the malignant cells, but not malignant programs. Two statistics per MP:

* ``partial_r2``: across samples, the share of variance in the per-sample malignant MP score
  explained by the same gene set's score in that sample's TME cells, beyond cancer type
  (and optionally study). TME scores are centred per cell type, so TME composition cancels.
* ``kappa_ratio``: within samples, co-expression (``kappa``) of the MP genes in TME cell types
  divided by that in malignant cells. Near 0: malignant-specific co-expression; near 1 or
  above: the same genes co-vary in TME cells too (technical or niche-wide).

Neither statistic alone says "artefact": IFN response and hypoxia are expected to be
TME-shared biology. They feed the provenance typing, together with the protocol axis.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .coexpr import kappa_with_null
from .data import Keys, get_counts, log2_cpm10, malignant_mask
from .scoring import score_gene_sets


def _design(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    parts = [np.ones((len(df), 1))]
    for c in cols:
        d = pd.get_dummies(df[c].astype(str), drop_first=True, dtype=float)
        if d.shape[1]:
            parts.append(d.to_numpy())
    return np.hstack(parts)


def _partial_r2(y: np.ndarray, x: np.ndarray, Z: np.ndarray) -> float:
    def sse(M):
        beta, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ beta
        return float(r @ r)
    red = sse(Z)
    full = sse(np.hstack([Z, x[:, None]]))
    return (red - full) / red if red > 0 else np.nan


def tme_covariation(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    min_cells: int = 20,
    covariates: tuple[str, ...] | None = None,
    n_boot: int = 200,
    kappa_ratio: bool = True,
    max_kappa_cells: int = 2000,
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
    scores = score_gene_sets(L, adata.var_names, gene_sets, np.flatnonzero(mal), seed=seed)
    scores.index = adata.obs_names

    covariates = tuple(c for c in (covariates if covariates is not None else (keys.cancer_type,)) if c)
    meta = adata.obs[[c for c in covariates]].astype(str).copy() if covariates else pd.DataFrame(index=adata.obs_names)
    meta["sample"] = sample

    per_sample_rows = []
    for s in pd.unique(sample):
        m_rows = mal & (sample == s)
        if m_rows.sum() < min_cells:
            continue
        row = {"sample": s, **{c: meta.loc[m_rows, c].iloc[0] for c in covariates}}
        for name in gene_sets:
            row[f"{name}::malignant"] = float(np.nanmean(scores[name].to_numpy()[m_rows]))
        for ct in pd.unique(ctype[tme & (sample == s)]):
            t_rows = tme & (sample == s) & (ctype == ct)
            if t_rows.sum() < min_cells:
                continue
            for name in gene_sets:
                row[f"{name}::tme::{ct}"] = float(np.nanmean(scores[name].to_numpy()[t_rows]))
        per_sample_rows.append(row)
    per_sample = pd.DataFrame(per_sample_rows)

    out_rows = []
    for name, genes in gene_sets.items():
        rec = {"mp": name, "n_samples": 0, "partial_r2": np.nan, "partial_r2_lo": np.nan,
               "partial_r2_hi": np.nan, "spearman": np.nan}
        if not per_sample.empty and f"{name}::malignant" in per_sample:
            tme_cols = [c for c in per_sample.columns if c.startswith(f"{name}::tme::")]
            if tme_cols:
                T = per_sample[tme_cols]
                T = T - T.mean(axis=0)  # centre each cell type across samples
                x = T.mean(axis=1, skipna=True).to_numpy()
                y = per_sample[f"{name}::malignant"].to_numpy()
                ok = np.isfinite(x) & np.isfinite(y)
                df = per_sample.loc[ok]
                Z = _design(df, list(covariates))
                if ok.sum() > Z.shape[1] + 2:
                    yv, xv = y[ok], x[ok]
                    rec["n_samples"] = int(ok.sum())
                    rec["partial_r2"] = _partial_r2(yv, xv, Z)
                    rec["spearman"] = float(pd.Series(yv).corr(pd.Series(xv), method="spearman"))
                    boots = []
                    for _ in range(n_boot):
                        b = rng.integers(0, ok.sum(), ok.sum())
                        Zb = _design(df.iloc[b], list(covariates))
                        if np.linalg.matrix_rank(Zb) < Zb.shape[1] or b.size <= Zb.shape[1] + 2:
                            continue
                        boots.append(_partial_r2(yv[b], xv[b], Zb))
                    boots = np.array([v for v in boots if np.isfinite(v)])
                    if boots.size:
                        rec["partial_r2_lo"], rec["partial_r2_hi"] = (float(q) for q in np.quantile(boots, [0.05, 0.95]))
        if kappa_ratio:
            pos = {g: i for i, g in enumerate(map(str, adata.var_names))}
            idx = np.array([pos[g] for g in genes if g in pos], dtype=int)
            ratios = []
            for s in pd.unique(sample):
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
                        v = kappa_with_null(counts[r[:max_kappa_cells]], idx, n_null=0)["kappa"]
                        if np.isfinite(v):
                            kt.append(v)
                            wt.append(r.size)
                if np.isfinite(km) and km > 0 and kt:
                    ratios.append(float(np.average(kt, weights=wt)) / km)
            rec["kappa_ratio_median"] = float(np.median(ratios)) if ratios else np.nan
            rec["kappa_ratio_n"] = len(ratios)
        out_rows.append(rec)
    return pd.DataFrame(out_rows), per_sample
