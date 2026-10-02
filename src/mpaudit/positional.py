"""Pilot P5 / axis A3 (expression-only part): are an MP's genes clustered on one chromosome arm?

Recurrent arm-level copy-number changes (e.g. chr7 gain in GBM) can produce "programs" whose
genes sit on one arm. A positional flag is a prompt for allele-aware follow-up, not proof of a
dosage artefact: some real programs are positionally clustered gene families. Those families
are reported separately (``KNOWN_CLUSTERS``) and enrichment is also computed without them.

Gene positions come from ``adata.var`` columns or a GTF (:func:`gene_positions_from_gtf`); arms
come from a UCSC ``cytoBand.txt`` (centromere = ``acen`` bands; :func:`load_cytoband`).
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

KNOWN_CLUSTERS = {
    "keratin": r"^KRT\d",
    "HLA_class_II": r"^HLA-D",
    "histone": r"^(H1-\d|H2AC|H2BC|H3C|H4C|HIST\d)",
    "metallothionein": r"^MT1[A-Z]*$|^MT2A$",
    "EDC_S100_SPRR_LCE": r"^(S100A\d|SPRR\d|LCE\d)",
    "protocadherin": r"^PCDH[ABG]",
}


def _norm_chrom(c) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def load_cytoband(path: str | Path) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, names=["chrom", "start", "end", "band", "stain"])
    df["chrom"] = df["chrom"].map(_norm_chrom)
    return df


def centromeres(cyto: pd.DataFrame) -> dict[str, float]:
    acen = cyto[cyto["stain"] == "acen"].groupby("chrom").agg(start=("start", "min"), end=("end", "max"))
    return {c: (r.start + r.end) / 2 for c, r in acen.iterrows()}


def gene_positions_from_gtf(path: str | Path, name_attr: str = "gene_name") -> pd.DataFrame:
    """Gene name, chromosome and start from the ``gene`` records of a GTF (optionally gzipped)."""
    pat = re.compile(rf'{name_attr} "([^"]+)"')
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = pat.search(f[8])
            if m:
                rows.append((m.group(1), _norm_chrom(f[0]), int(f[3])))
    return pd.DataFrame(rows, columns=["gene", "chromosome", "start"]).drop_duplicates("gene").set_index("gene")


def assign_arms(positions: pd.DataFrame, cen: dict[str, float],
                chrom_col: str = "chromosome", start_col: str = "start") -> pd.Series:
    """Arm label (e.g. ``7p``/``7q``) per gene; NaN where position or centromere is unknown."""
    chrom = positions[chrom_col].map(_norm_chrom)
    start = pd.to_numeric(positions[start_col], errors="coerce")
    arms = []
    for c, s in zip(chrom, start):
        if c not in cen or not np.isfinite(s):
            arms.append(np.nan)
        else:
            arms.append(f"{c}{'p' if s < cen[c] else 'q'}")
    return pd.Series(arms, index=positions.index, name="arm")


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def known_cluster_genes(genes) -> set[str]:
    pats = [re.compile(p) for p in KNOWN_CLUSTERS.values()]
    return {g for g in genes if any(p.search(g) for p in pats)}


def arm_enrichment(gene_sets: dict[str, list[str]], gene_arm: pd.Series, universe=None,
                   exclude_known: bool = False) -> pd.DataFrame:
    """Hypergeometric over-representation of each set's genes on each arm (BH within set).

    ``universe`` defaults to all genes with an arm. Returns one row per set: best arm, its
    count, expected count, fold, p and q, plus how many of the set's genes were placed.
    """
    arm = gene_arm.dropna()
    if universe is not None:
        arm = arm[arm.index.isin(list(universe))]
    if exclude_known:
        arm = arm[~arm.index.isin(known_cluster_genes(arm.index))]
    N = arm.size
    arm_sizes = arm.value_counts()
    rows = []
    for name, genes in gene_sets.items():
        placed = arm[arm.index.isin(genes)]
        n = placed.size
        rec = {"mp": name, "n_placed": n, "best_arm": "", "k": 0, "expected": np.nan,
               "fold": np.nan, "p": 1.0, "q": 1.0,
               "known_cluster_genes": ",".join(sorted(known_cluster_genes(genes)))}
        if n >= 3:
            counts = placed.value_counts()
            arms = arm_sizes.index.to_numpy()
            k = counts.reindex(arms).fillna(0).to_numpy()
            K = arm_sizes.to_numpy()
            p = hypergeom.sf(k - 1, N, K, n)
            q = _bh(p)
            i = int(np.argmin(p))
            rec.update(best_arm=str(arms[i]), k=int(k[i]), expected=float(K[i] * n / N),
                       fold=float(k[i] / (K[i] * n / N)), p=float(p[i]), q=float(q[i]))
        rows.append(rec)
    return pd.DataFrame(rows)
