"""Comparing MP gene sets between runs."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def jaccard_matrix(a: dict[str, list[str]], b: dict[str, list[str]]) -> np.ndarray:
    sa = [set(v) for v in a.values()]
    sb = [set(v) for v in b.values()]
    J = np.zeros((len(sa), len(sb)))
    for i, x in enumerate(sa):
        for j, y in enumerate(sb):
            u = len(x | y)
            J[i, j] = len(x & y) / u if u else 0.0
    return J


def match_sets(ref: dict[str, list[str]], alt: dict[str, list[str]], min_jaccard: float = 0.3) -> pd.DataFrame:
    """One-to-one (Hungarian) matching of reference sets to alternative sets by Jaccard.

    Columns: ``ref``, ``alt`` (matched set or empty), ``jaccard`` (of the assigned pair),
    ``best_jaccard`` (best over all alternative sets, ignoring one-to-one), ``matched``.
    """
    names_r, names_a = list(ref), list(alt)
    rows = []
    if not names_r:
        return pd.DataFrame(columns=["ref", "alt", "jaccard", "best_jaccard", "matched"])
    J = jaccard_matrix(ref, alt) if names_a else np.zeros((len(names_r), 0))
    assigned = {}
    if names_a:
        ri, ci = linear_sum_assignment(-J)
        assigned = {r: c for r, c in zip(ri, ci)}
    for i, r in enumerate(names_r):
        c = assigned.get(i)
        jac = float(J[i, c]) if c is not None else 0.0
        best = float(J[i].max()) if names_a else 0.0
        rows.append({"ref": r, "alt": names_a[c] if c is not None and jac >= min_jaccard else "",
                     "jaccard": jac, "best_jaccard": best, "matched": jac >= min_jaccard})
    return pd.DataFrame(rows)


def summarize_match(m: pd.DataFrame, n_alt: int) -> dict:
    n_ref = len(m)
    n_matched = int(m["matched"].sum()) if n_ref else 0
    return {
        "n_ref": n_ref,
        "n_alt": n_alt,
        "n_matched": n_matched,
        "recall": n_matched / n_ref if n_ref else np.nan,
        "precision": n_matched / n_alt if n_alt else np.nan,
        "mean_jaccard_matched": float(m.loc[m["matched"], "jaccard"].mean()) if n_matched else np.nan,
    }
