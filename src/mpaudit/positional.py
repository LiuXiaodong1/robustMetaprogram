"""Pilot P5 / axis A3 (expression-only part): are an MP's genes clustered on one chromosome arm?

Recurrent arm-level copy-number changes (e.g. chr7 gain in GBM) can produce "programs" whose
genes sit on one arm. A positional flag is a prompt for allele-aware follow-up, not proof of a
dosage artefact: some real programs are positionally clustered gene families. Those families
are reported separately (``KNOWN_CLUSTERS``) and enrichment is also computed without them.

Gene positions come from ``adata.var`` columns or a GTF (:func:`gene_positions_from_gtf`); arms
come from a UCSC ``cytoBand.txt`` (centromere = ``acen`` bands; :func:`load_cytoband`).

The same sources give gene lengths (:func:`gene_lengths`) for the typer's gene-length features.
Length here is the genomic span of the GTF ``gene`` record (first to last base of any transcript,
introns included), which is what governs how many intronic reads a gene collects in
single-nucleus data; it is not the mature transcript length. GTFs without ``gene`` records (e.g.
UCSC's ``ncbiRefSeq.gtf`` and ``knownGene.gtf``) give the same span from the gene's transcript,
exon and other records, with a warning.
"""
from __future__ import annotations

import gzip
import re
import warnings
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
    """Chromosome, start, end and length per gene from a GTF (optionally gzipped).

    Indexed by ``name_attr``. Coordinates are GTF's (1-based, closed), so ``length = end - start + 1``
    bp. Spans come from the ``gene`` records (GENCODE, Ensembl); a name with several ``gene``
    records (e.g. PAR genes on X and Y) keeps the first.

    For a GTF without any ``gene`` record (e.g. UCSC ``ncbiRefSeq.gtf``/``knownGene.gtf``, which
    have only transcript, exon, CDS and codon lines) the span is aggregated over all records of
    each name instead (min start, max end), on one chromosome: the first one the name appears on,
    preferring primary chromosomes over contigs whose name has an underscore
    (``chr7_KI270803v1_alt``, ``chrUn_...``); records on other chromosomes are ignored. This
    fallback warns (``UserWarning``). So does a GTF in which no record carries ``name_attr``: the
    table is then empty, which leaves every gene without a position or length.
    """
    pat = re.compile(rf'{name_attr} "([^"]+)"')
    genes = []
    spans: dict[str, list] = {}  # name -> [chrom, start, end] from non-gene records, while no gene record is seen
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            is_gene = f[2] == "gene"
            if not is_gene and genes:
                continue
            m = pat.search(f[8])
            if not m:
                continue
            name, chrom, start, end = m.group(1), _norm_chrom(f[0]), int(f[3]), int(f[4])
            if is_gene:
                genes.append((name, chrom, start, end))
                spans = {}
                continue
            cur = spans.get(name)
            if cur is None or ("_" in cur[0] and "_" not in chrom):
                spans[name] = [chrom, start, end]
            elif chrom == cur[0]:
                cur[1], cur[2] = min(cur[1], start), max(cur[2], end)
    rows = genes
    if not genes and spans:
        warnings.warn(f"{path}: no 'gene' records; spans of {len(spans)} genes taken from their transcript/exon "
                      f"records (min start, max end per {name_attr})", UserWarning, stacklevel=2)
        rows = [(name, *span) for name, span in spans.items()]
    elif not genes:
        warnings.warn(f"{path}: no record carries a {name_attr} attribute; no gene positions", UserWarning,
                      stacklevel=2)
    df = pd.DataFrame(rows, columns=["gene", "chromosome", "start", "end"]).drop_duplicates("gene").set_index("gene")
    df["length"] = df["end"] - df["start"] + 1
    return df


def gene_lengths(table: pd.DataFrame, start_col: str = "start", end_col: str = "end") -> pd.Series:
    """Gene length in bp (genomic span), indexed like ``table`` (``adata.var`` or a positions table).

    Taken from a ``gene_length`` column, else ``length`` (as from :func:`gene_positions_from_gtf`),
    else ``end - start + 1`` (1-based closed coordinates, as in GTF and Ensembl BioMart; the 1 bp
    convention difference is immaterial at gene scale). Missing or non-positive values are NaN.
    Raises ``ValueError`` when none of these columns exists.
    """
    for col in ("gene_length", "length"):
        if col in table:
            bp = pd.to_numeric(table[col], errors="coerce")
            break
    else:
        if start_col not in table or end_col not in table:
            raise ValueError(f"no gene_length, length or {start_col}/{end_col} columns to take gene lengths from")
        bp = pd.to_numeric(table[end_col], errors="coerce") - pd.to_numeric(table[start_col], errors="coerce") + 1
    return bp.where(bp > 0).astype(float).rename("gene_length")


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
