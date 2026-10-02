"""Pilot P4: the trivial artefact typer that any provenance method must beat.

Cheap per-MP features (gene-name content, TME specificity, correlation of the MP score with
QC covariates), and their ability to separate labelled artefact MPs (e.g. the 11 MPs removed
by hand in 3CA v2) from retained MPs. Every feature is oriented so that higher = more
artefact-like; AUROC is reported with a Hanley-McNeil standard error, per feature and for a
leave-one-out ridge-logistic combination.

The IEG/heat-shock and TME-marker lists below are short illustrative core lists; replace them
with curated lists (e.g. the van den Brink 2017 dissociation signature) for real analyses.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

from .data import Keys, get_counts, log2_cpm10, malignant_mask
from .scoring import score_gene_sets

IEG_HSP_CORE = ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3",
                "NR4A1", "DUSP1", "ZFP36", "KLF2", "KLF4", "KLF6", "HSPA1A", "HSPA1B", "HSPA8",
                "HSP90AA1", "DNAJB1", "HSPH1", "CCN1"]
TME_MARKERS_CORE = ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GZMA", "CD79A", "MS4A1", "LYZ", "CD68",
                    "CD14", "C1QA", "C1QB", "TYROBP", "FCER1G", "AIF1", "PECAM1", "VWF", "CDH5",
                    "CLDN5", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "HBB", "HBA1", "HBA2"]
_RIBO = re.compile(r"^RP[SL]\d")
_MITO = re.compile(r"^MT-")


def typer_features(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    tme_fold: float = 4.0,
    seed: int = 0,
) -> pd.DataFrame:
    counts = get_counts(adata, counts_layer)
    L = log2_cpm10(counts)
    mal = malignant_mask(adata, keys)
    genes = np.asarray(adata.var_names.astype(str))
    pos = {g: i for i, g in enumerate(genes)}

    tme_specific = set()
    if (~mal).sum() >= 10:
        m_mean = np.asarray(L[mal].mean(axis=0)).ravel()
        t_mean = np.asarray(L[~mal].mean(axis=0)).ravel()
        tme_specific = set(genes[(t_mean - m_mean) >= np.log2(tme_fold)])

    totals = np.asarray(counts.sum(axis=1)).ravel()
    detected = np.asarray((counts > 0).sum(axis=1)).ravel()
    mito_idx = [pos[g] for g in genes if _MITO.match(g)]
    pct_mito = (np.asarray(counts[:, mito_idx].sum(axis=1)).ravel() / np.maximum(totals, 1)) if mito_idx else np.zeros_like(totals)
    qc = {"log_umi": np.log1p(totals), "log_genes": np.log1p(detected), "pct_mito": pct_mito}
    if keys.doublet and keys.doublet in adata.obs:
        qc["doublet"] = pd.to_numeric(adata.obs[keys.doublet], errors="coerce").to_numpy()

    scores = score_gene_sets(L, genes, gene_sets, np.flatnonzero(mal), seed=seed)
    sample = adata.obs[keys.sample].astype(str).to_numpy()

    rows = []
    for name, gs in gene_sets.items():
        n = max(len(gs), 1)
        rec = {
            "mp": name,
            "frac_ribo": sum(bool(_RIBO.match(g)) for g in gs) / n,
            "frac_mito": sum(bool(_MITO.match(g)) for g in gs) / n,
            "frac_ieg_hsp": len(set(gs) & set(IEG_HSP_CORE)) / n,
            "frac_tme_marker": len(set(gs) & set(TME_MARKERS_CORE)) / n,
            "frac_tme_specific": len(set(gs) & tme_specific) / n,
        }
        sc = scores[name].to_numpy()
        for q, v in qc.items():
            cors = []
            for s in pd.unique(sample[mal]):
                r = mal & (sample == s)
                a, b = sc[r], v[r]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() >= 20 and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
                    cors.append(np.corrcoef(a[ok], b[ok])[0, 1])
            rec[f"abs_cor_{q}"] = float(np.median(np.abs(cors))) if cors else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def auroc(score: np.ndarray, label: np.ndarray) -> tuple[float, float]:
    """AUROC (positives = label 1) with Hanley-McNeil standard error."""
    score, label = np.asarray(score, float), np.asarray(label).astype(bool)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n1, n0 = int(label.sum()), int((~label).sum())
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan
    ranks = pd.Series(score).rank().to_numpy()
    auc = (ranks[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    q1, q2 = auc / (2 - auc), 2 * auc**2 / (1 + auc)
    se = np.sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc**2) + (n0 - 1) * (q2 - auc**2)) / (n1 * n0))
    return float(auc), float(se)


def _ridge_logit(X: np.ndarray, y: np.ndarray, lam: float = 1.0, iters: int = 50) -> np.ndarray:
    Xb = np.hstack([np.ones((X.shape[0], 1)), X])
    w = np.zeros(Xb.shape[1])
    P = lam * np.eye(Xb.shape[1])
    P[0, 0] = 0.0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-Xb @ w))
        g = Xb.T @ (p - y) + P @ w
        Hs = Xb.T @ (Xb * (p * (1 - p))[:, None]) + P
        w -= np.linalg.solve(Hs, g)
    return w


def evaluate_typer(features: pd.DataFrame, labels: pd.Series, lam: float = 1.0) -> pd.DataFrame:
    """Per-feature AUROC and a leave-one-out ridge-logistic combination."""
    df = features.set_index("mp").join(labels.rename("label"), how="inner").dropna(subset=["label"])
    feats = [c for c in features.columns if c != "mp"]
    y = df["label"].astype(int).to_numpy()
    rows = []
    for f in feats:
        a, se = auroc(df[f].to_numpy(), y)
        rows.append({"feature": f, "auroc": a, "se": se, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())})
    X = df[feats].to_numpy(dtype=float)
    X = np.where(np.isfinite(X), X, np.nanmean(np.where(np.isfinite(X), X, np.nan), axis=0))
    X = np.nan_to_num(X)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    X = (X - X.mean(axis=0)) / sd
    loo = np.full(len(y), np.nan)
    if y.sum() >= 2 and (1 - y).sum() >= 2:
        for i in range(len(y)):
            tr = np.arange(len(y)) != i
            w = _ridge_logit(X[tr], y[tr], lam)
            loo[i] = w[0] + X[i] @ w[1:]
    a, se = auroc(loo, y)
    rows.append({"feature": "combined_loo_logistic", "auroc": a, "se": se, "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())})
    return pd.DataFrame(rows)
