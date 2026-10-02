"""Cohort data handling: column keys, counts access, normalisation, gene-set I/O."""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass
class Keys:
    """Names of the ``adata.obs`` columns the audit uses.

    ``sample`` is the unit NMF is run on (one tumour sample, as in 3CA). ``tumour`` groups
    several samples of the same tumour (protocol twins, regions); when it is unset each sample
    is its own tumour. Optional keys switch on the analyses that need them.
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
    """log2(CPM/10 + 1), the 3CA transform. CPM uses each cell's total over all genes."""
    counts = sp.csr_matrix(counts, dtype=np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    scale = sp.diags(1e5 / totals)
    out = scale @ counts
    out.data = np.log2(out.data + 1.0)
    return out.tocsr()


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
    those are loaded as-is and flagged in ``uns['3ca_value_type']``.
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
