"""Cohort data handling: column keys, counts access, normalisation, per-cell QC, gene-set I/O."""
from __future__ import annotations

import dataclasses
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

# HGNC symbols (human). Cytoplasmic ribosomal-protein genes, including RPLP0-2, RPSA and the -like
# paralogs (RPL3L, RPS27L, ...), but not the RPS6K kinases, RPS19BP1, pseudogenes (RPL13AP5,
# RPSAP58) or read-throughs (RPS10-NUDT3), which a plain RPL/RPS prefix would include.
RIBO_PATTERN = re.compile(r"^(RP[LS]\d+[A-OQ-Z]?\d?|RPL36AL|RPLP[0-2]|RPSA)$")
MITO_PATTERN = re.compile(r"^MT-")
QC_COLUMNS = ("log_umi", "log_genes", "pct_mito", "pct_ribo")
LONG_GENE_BP = 100_000  # genomic span above which a gene is "long" (pct_long; the typer's frac_long_genes)


@dataclass
class Keys:
    """Names of the ``adata.obs`` columns the audit uses.

    ``sample`` is the unit NMF is run on (one tumour sample, as in 3CA). ``tumour`` groups
    several samples of the same tumour (protocol twins, regions); when it is unset each sample
    is its own tumour. Optional keys switch on the analyses that need them. ``intronic`` is a
    per-cell intronic (unspliced, or nuclear) read fraction, e.g. from STARsolo Velocyto or
    DropletQC; the typer correlates MP scores with it.
    """

    sample: str = "sample"
    malignant: str = "malignant"
    malignant_label: str | None = None
    tumour: str | None = None
    study: str | None = None
    cancer_type: str | None = None
    cell_type: str | None = None
    protocol: str | None = None
    region: str | None = None
    clone: str | None = None
    doublet: str | None = None
    intronic: str | None = None

    @classmethod
    def from_dict(cls, d: dict | None) -> "Keys":
        d = dict(d or {})
        names = {f.name for f in dataclasses.fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"Unknown key(s): {sorted(unknown)}")
        return cls(**d)

    def tumour_col(self) -> str:
        return self.tumour or self.sample


def get_counts(adata, layer: str | None = None) -> sp.csr_matrix:
    """Raw counts (cells x genes) as CSR float64."""
    X = adata.layers[layer] if layer else adata.X
    X = sp.csr_matrix(X, dtype=np.float64) if not sp.issparse(X) else X.tocsr().astype(np.float64)
    return X


def malignant_mask(adata, keys: Keys) -> np.ndarray:
    """Boolean mask of malignant cells from a bool/0-1 column or a label column."""
    col = adata.obs[keys.malignant]
    if keys.malignant_label is not None:
        return (col.astype(str) == str(keys.malignant_label)).to_numpy()
    if pd.api.types.is_bool_dtype(col) or pd.api.types.is_numeric_dtype(col):
        return col.fillna(0).astype(bool).to_numpy()
    s = col.astype(str).str.lower()
    return s.isin(["true", "yes", "malignant", "1"]).to_numpy()


def log2_cpm10(counts: sp.spmatrix) -> sp.csr_matrix:
    """log2(CPM/10 + 1), inferred to be the 3CA transform (the 3CA MP code log-normalises inside
    the private ``Return_All_cells_Separated``; see :mod:`mpaudit.gavish`). CPM uses each cell's
    total over all genes, which is an assumption."""
    counts = sp.csr_matrix(counts, dtype=np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    scale = sp.diags(1e5 / totals)
    out = scale @ counts
    out.data = np.log2(out.data + 1.0)
    return out.tocsr()


def _long_genes(genes: np.ndarray, gene_length: pd.Series) -> np.ndarray:
    """Positions in ``genes`` of the genes longer than :data:`LONG_GENE_BP` (no length: not long)."""
    bp = pd.to_numeric(pd.Series(gene_length), errors="coerce")
    bp.index = bp.index.astype(str)
    bp = bp[~bp.index.duplicated()]
    return np.flatnonzero(bp.reindex(genes).to_numpy(dtype=float) > LONG_GENE_BP)


def qc_metrics(counts: sp.spmatrix, var_names, gene_length: pd.Series | None = None) -> pd.DataFrame:
    """Per-cell QC covariates from raw counts (cells x genes), one row per cell.

    ``log_umi`` = log(1 + total counts), ``log_genes`` = log(1 + detected genes); ``pct_mito`` and
    ``pct_ribo`` = percentage (0-100) of a cell's counts on mitochondrial (:data:`MITO_PATTERN`) and
    cytoplasmic ribosomal-protein genes (:data:`RIBO_PATTERN`); 0 for cells without counts. With
    ``gene_length`` (bp, indexed by gene name, e.g. :func:`mpaudit.positional.gene_lengths`), also
    ``pct_long``: the percentage of a cell's counts on genes longer than :data:`LONG_GENE_BP`
    (genes without a length count as not long).
    """
    counts = sp.csr_matrix(counts, dtype=np.float64)
    genes = np.asarray(var_names, dtype=str)
    totals = np.asarray(counts.sum(axis=1)).ravel()

    def pct(idx: np.ndarray) -> np.ndarray:
        return 100.0 * np.asarray(counts[:, idx].sum(axis=1)).ravel() / np.maximum(totals, 1.0)

    def matching(pattern: re.Pattern) -> np.ndarray:
        return np.flatnonzero([bool(pattern.match(g)) for g in genes])

    out = pd.DataFrame({"log_umi": np.log1p(totals),
                        "log_genes": np.log1p(np.asarray((counts > 0).sum(axis=1)).ravel()),
                        "pct_mito": pct(matching(MITO_PATTERN)), "pct_ribo": pct(matching(RIBO_PATTERN))})
    if gene_length is not None:
        out["pct_long"] = pct(_long_genes(genes, gene_length))
    return out


def add_qc_obs(adata, counts_layer: str | None = None, overwrite: bool = False,
               gene_length: pd.Series | None = None):
    """Add the :func:`qc_metrics` columns to ``adata.obs`` in place and return ``adata``.

    Columns already in ``obs`` are kept unless ``overwrite`` is set, so repeated calls are harmless
    and an analyst's own QC values win. Warns when no gene matches the mitochondrial or ribosomal
    pattern (e.g. Ensembl IDs as ``var_names``, or mouse symbols): that column is then all zero.

    ``pct_long`` is added only when gene lengths are available: ``gene_length`` (bp, indexed by
    gene name; :func:`mpaudit.pipeline.run_p1` passes the lengths of ``positional.gtf``), else an
    ``adata.var['gene_length']`` column. It can serve as a continuous technical decoy covariate
    for the long-gene (nuclear/intronic) axis where no measured ``intronic_fraction`` exists; it is
    in no decoy arm by default. Caveat (as for ``pct_ribo``; see :mod:`mpaudit.decoy`): it is
    computed from counts that include an audited MP's own long genes, so a decoy conditioned on it
    keeps part of a long-gene MP's sampling noise as co-expression and can inflate that MP's
    ``kappa_retained``.
    """
    if gene_length is None and "gene_length" in adata.var:
        gene_length = pd.Series(adata.var["gene_length"].to_numpy(), index=adata.var_names.astype(str))
    columns = QC_COLUMNS + (("pct_long",) if gene_length is not None else ())
    todo = [c for c in columns if overwrite or c not in adata.obs]
    if not todo:
        return adata
    genes = adata.var_names.astype(str)
    for col, pattern in (("pct_mito", MITO_PATTERN), ("pct_ribo", RIBO_PATTERN)):
        if col in todo and not any(pattern.match(g) for g in genes):
            warnings.warn(f"no gene matches {pattern.pattern!r}; {col} is 0 for every cell", stacklevel=2)
    qc = qc_metrics(get_counts(adata, counts_layer), genes, gene_length=gene_length)
    if "pct_long" in todo and not (qc["pct_long"] > 0).any():
        warnings.warn(f"no counts on genes longer than {LONG_GENE_BP} bp with a known length; pct_long is 0 for "
                      "every cell (do the gene lengths match var_names?)", stacklevel=2)
    for col in todo:
        adata.obs[col] = qc[col].to_numpy()
    return adata


def sample_order(adata, keys: Keys, mask: np.ndarray | None = None) -> list[str]:
    """Samples in order of first appearance (deterministic)."""
    s = adata.obs[keys.sample].astype(str)
    if mask is not None:
        s = s[mask]
    return list(pd.unique(s))


def sample_map(adata, keys: Keys, col: str | None) -> dict[str, str]:
    """Map sample -> value of ``col`` (first occurrence); identity when ``col`` is None."""
    samples = adata.obs[keys.sample].astype(str)
    if col is None:
        return {s: s for s in pd.unique(samples)}
    vals = adata.obs[col].astype(str)
    df = pd.DataFrame({"s": samples.to_numpy(), "v": vals.to_numpy()}).drop_duplicates("s")
    return dict(zip(df["s"], df["v"]))


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for line in Path(path).read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def write_gmt(sets: dict[str, list[str]], path: str | Path, description: str = "mpaudit") -> None:
    lines = ["\t".join([name, description, *genes]) for name, genes in sets.items()]
    Path(path).write_text("\n".join(lines) + "\n")


def load_3ca_study(path: str | Path):
    """Load one study folder downloaded from the 3CA website into AnnData.

    Assumed layout (check against your download): a Matrix Market file of counts
    (``*UMIcounts*.mtx`` preferred, otherwise any ``*.mtx``), ``Genes.txt`` (one gene per line)
    and ``Cells.csv`` (one row per cell; first column the cell name). Orientation of the matrix
    is inferred from the gene and cell counts. Smart-seq studies ship TPM rather than counts;
    those are loaded as-is and flagged in ``uns['3ca_value_type']``. The folders hold gene-level
    counts only: no spliced/unspliced or intronic information, so a per-cell intronic fraction
    needs the raw reads re-processed (e.g. STARsolo Velocyto, velocyto, DropletQC).
    """
    import anndata as ad
    from scipy.io import mmread

    path = Path(path)
    mtx = sorted(path.glob("*UMIcounts*.mtx")) or sorted(path.glob("*.mtx"))
    if not mtx:
        raise FileNotFoundError(f"No .mtx file in {path}")
    genes = [g.strip() for g in (path / "Genes.txt").read_text().splitlines() if g.strip()]
    cells = pd.read_csv(path / "Cells.csv")
    cells.index = cells.iloc[:, 0].astype(str)
    M = sp.csr_matrix(mmread(mtx[0]))
    if M.shape == (len(genes), len(cells)):
        M = M.T.tocsr()
    elif M.shape != (len(cells), len(genes)):
        raise ValueError(f"Matrix shape {M.shape} matches neither genes x cells nor cells x genes")
    adata = ad.AnnData(X=M, obs=cells, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    adata.uns["3ca_value_type"] = "counts" if "UMI" in mtx[0].name else "unknown (possibly TPM)"
    return adata
