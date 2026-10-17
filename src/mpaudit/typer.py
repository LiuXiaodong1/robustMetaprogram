"""Pilot P4: the trivial artefact typer that any provenance method must beat.

Cheap per-MP features, and their ability to separate labelled artefact MPs (e.g. the 11 MPs
removed by hand in 3CA v2) from retained MPs. Every feature is oriented so that higher = more
artefact-like; AUROC is reported with a Hanley-McNeil standard error, per feature and for a
leave-one-out ridge-logistic combination. Features (:func:`typer_features`):

* gene content: fractions of the MP's genes that are ribosomal-protein, mitochondrial,
  immediate-early/heat-shock, TME markers, or TME-specific in this cohort (``frac_*``);
* gene length: ``median_gene_length_kb`` and ``frac_long_genes`` (genomic span > 100 kb);
* QC co-variation: in each sample's malignant cells (>= 20), the Pearson correlation of the MP
  score with a per-cell covariate; the feature is the median |r| over samples (``abs_cor_*``).
  Covariates: log UMIs, log detected genes, %mito and %ribo (computed as
  :func:`mpaudit.data.add_qc_obs` does), and, when their ``obs`` columns are given, the doublet
  score (``keys.doublet``) and the intronic fraction (``keys.intronic``).

The length and %ribo/intronic features target two technical axes that the other features miss.
Long genes collect many intronic reads, so in single-nucleus data (and wherever pre-mRNA is
captured) their counts rise with a cell's nuclear/intronic share; abundant short cytoplasmic
transcripts (ribosomal proteins and similar) rise with the cytoplasmic share retained, which
%ribo tracks. Either axis varies between cells of one sample, so per-sample NMF returns it as a
program, and recurrence across samples makes it a metaprogram. User-reported evidence (a breast
cancer snRNA-seq cohort; not reproduced here): two pre-registered de novo MP searches returned
mostly technical axes, among them several "luminal" programs made of long genes (median gene
length 317-437 kb) and a short-gene abundance axis whose score tracked %ribo in 29 of 30 patients
(rho 0.68-0.77), labelled "EMT" by reference matching. Matching to reference programs therefore
does not screen these axes out; the features here are a first screen, not a verdict.

Caveats:

* The length features are oriented for long-gene artefacts. A short-gene abundance program has a
  LOW median length; it is caught by ``abs_cor_pct_ribo`` (and ``frac_ribo``), not by length.
* Genomic span is not neutral: neuronal, synaptic and cell-adhesion genes are long, so some real
  programs are long-gene rich. A length flag asks for a check against the intronic fraction.
* For an MP containing ribosomal-protein genes, ``abs_cor_pct_ribo`` is partly mechanical (the
  score and %ribo share those counts); ``frac_ribo`` flags such MPs directly.
* Correlations are within sample, so a nuisance that acts only between samples is not seen, and
  the intronic fraction must be measured per cell (e.g. STARsolo Velocyto, DropletQC).
* Features whose input is missing are left out of the table: the length features without gene
  lengths, ``abs_cor_doublet`` / ``abs_cor_intronic`` without the ``obs`` column.

The IEG/heat-shock and TME-marker lists below are short illustrative core lists; replace them
with curated lists (e.g. the van den Brink 2017 dissociation signature) for real analyses.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .data import MITO_PATTERN, RIBO_PATTERN, Keys, get_counts, log2_cpm10, malignant_mask, qc_metrics
from .positional import gene_lengths
from .scoring import score_gene_sets

IEG_HSP_CORE = ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3",
                "NR4A1", "DUSP1", "ZFP36", "KLF2", "KLF4", "KLF6", "HSPA1A", "HSPA1B", "HSPA8",
                "HSP90AA1", "DNAJB1", "HSPH1", "CCN1"]
TME_MARKERS_CORE = ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GZMA", "CD79A", "MS4A1", "LYZ", "CD68",
                    "CD14", "C1QA", "C1QB", "TYROBP", "FCER1G", "AIF1", "PECAM1", "VWF", "CDH5",
                    "CLDN5", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "HBB", "HBA1", "HBA2"]
LONG_GENE_BP = 100_000
MIN_CELLS_COR = 20


def typer_features(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    gene_length: pd.Series | None = None,
    tme_fold: float = 4.0,
    seed: int = 0,
) -> pd.DataFrame:
    """One row per gene set with the features described in the module docstring.

    ``gene_length`` is in bp, indexed by gene name (e.g. :func:`mpaudit.positional.gene_lengths` of
    a GTF table); by default it is taken from ``adata.var`` when that has ``gene_length`` or
    ``start``/``end`` columns. Length features use the set's genes with a known length (NaN if none).
    """
    counts = get_counts(adata, counts_layer)
    L = log2_cpm10(counts)
    mal = malignant_mask(adata, keys)
    genes = np.asarray(adata.var_names.astype(str))

    tme_specific = set()
    if (~mal).sum() >= 10:
        m_mean = np.asarray(L[mal].mean(axis=0)).ravel()
        t_mean = np.asarray(L[~mal].mean(axis=0)).ravel()
        tme_specific = set(genes[(t_mean - m_mean) >= np.log2(tme_fold)])

    if gene_length is None:
        try:
            gene_length = gene_lengths(adata.var)
        except ValueError:
            pass
    if gene_length is not None:
        gene_length = pd.Series(pd.to_numeric(gene_length, errors="coerce").to_numpy(dtype=float),
                                index=gene_length.index.astype(str))
        gene_length = gene_length[~gene_length.index.duplicated()]

    qc = {c: v.to_numpy() for c, v in qc_metrics(counts, genes).items()}
    for name, col in (("doublet", keys.doublet), ("intronic", keys.intronic)):
        if col and col in adata.obs:
            qc[name] = pd.to_numeric(adata.obs[col], errors="coerce").to_numpy()

    scores = score_gene_sets(L, genes, gene_sets, np.flatnonzero(mal), seed=seed)
    sample = adata.obs[keys.sample].astype(str).to_numpy()
    groups = [mal & (sample == s) for s in pd.unique(sample[mal])]

    rows = []
    for name, gs in gene_sets.items():
        n = max(len(gs), 1)
        rec = {
            "mp": name,
            "frac_ribo": sum(bool(RIBO_PATTERN.match(g)) for g in gs) / n,
            "frac_mito": sum(bool(MITO_PATTERN.match(g)) for g in gs) / n,
            "frac_ieg_hsp": len(set(gs) & set(IEG_HSP_CORE)) / n,
            "frac_tme_marker": len(set(gs) & set(TME_MARKERS_CORE)) / n,
            "frac_tme_specific": len(set(gs) & tme_specific) / n,
        }
        if gene_length is not None:
            bp = gene_length.reindex(list(dict.fromkeys(gs))).dropna().to_numpy()
            rec["median_gene_length_kb"] = float(np.median(bp)) / 1e3 if bp.size else np.nan
            rec["frac_long_genes"] = float(np.mean(bp > LONG_GENE_BP)) if bp.size else np.nan
        sc = scores[name].to_numpy()
        for q, v in qc.items():
            cors = []
            for r in groups:
                a, b = sc[r], v[r]
                ok = np.isfinite(a) & np.isfinite(b)
                if ok.sum() >= MIN_CELLS_COR and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
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
