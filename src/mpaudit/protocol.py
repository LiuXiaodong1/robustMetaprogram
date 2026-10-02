"""Axis A4: does an MP change between protocols applied to the same tumour?

For tumours profiled with two protocols (e.g. fresh dissociation vs fixed/Flex vs nuclei),
compares per-tumour MP intensity (mean score) and prevalence (fraction of malignant cells above
a high-score threshold) between protocols. The null is the same contrast between two samples of
the same tumour under the *same* protocol (different regions or pieces), which captures regional
heterogeneity. A protocol effect counts only if it exceeds that null.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

from .data import Keys, get_counts, log2_cpm10, malignant_mask
from .scoring import score_gene_sets


def per_sample_levels(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    high_quantile: float = 0.9,
    min_cells: int = 20,
    seed: int = 0,
) -> pd.DataFrame:
    """Per malignant sample: intensity and prevalence of each MP, with tumour/protocol labels."""
    if keys.protocol is None:
        raise ValueError("keys.protocol is required")
    L = log2_cpm10(get_counts(adata, counts_layer))
    mal = malignant_mask(adata, keys)
    scores = score_gene_sets(L, adata.var_names, gene_sets, np.flatnonzero(mal), seed=seed)
    obs = adata.obs
    tumour = obs[keys.tumour_col()].astype(str).to_numpy()
    sample = obs[keys.sample].astype(str).to_numpy()
    protocol = obs[keys.protocol].astype(str).to_numpy()
    thr = {n: np.nanquantile(scores[n].to_numpy()[mal], high_quantile) for n in gene_sets}
    rows = []
    for s in pd.unique(sample[mal]):
        r = mal & (sample == s)
        if r.sum() < min_cells:
            continue
        for n in gene_sets:
            v = scores[n].to_numpy()[r]
            rows.append({"sample": s, "tumour": tumour[r][0], "protocol": protocol[r][0], "mp": n,
                         "intensity": float(np.nanmean(v)), "prevalence": float(np.nanmean(v > thr[n]))})
    return pd.DataFrame(rows)


def protocol_shift(levels: pd.DataFrame, protocol_a: str, protocol_b: str) -> pd.DataFrame:
    """Paired ``b - a`` differences per MP, against same-protocol replicate differences.

    ``frac_loss_beyond_null`` is the fraction of paired tumours whose drop exceeds the 95th
    percentile of absolute same-protocol differences (one-sided); likewise for gains.
    """
    out = []
    for mp, d in levels.groupby("mp"):
        deltas, nulls = {"intensity": [], "prevalence": []}, {"intensity": [], "prevalence": []}
        for _, t in d.groupby("tumour"):
            a, b = t[t["protocol"] == protocol_a], t[t["protocol"] == protocol_b]
            if len(a) and len(b):
                for m in deltas:
                    deltas[m].append(b[m].mean() - a[m].mean())
            for _, same in t.groupby("protocol"):
                for (_, x), (_, y) in itertools.combinations(same.iterrows(), 2):
                    for m in nulls:
                        nulls[m].append(y[m] - x[m])
        rec = {"mp": mp, "n_paired_tumours": len(deltas["intensity"]), "n_null_pairs": len(nulls["intensity"])}
        for m in deltas:
            dv, nv = np.array(deltas[m]), np.abs(np.array(nulls[m]))
            q = float(np.quantile(nv, 0.95)) if nv.size else np.nan
            ok = dv.size and nv.size
            rec[f"mean_delta_{m}"] = float(dv.mean()) if dv.size else np.nan
            rec[f"null_abs_q95_{m}"] = q
            rec[f"frac_loss_beyond_null_{m}"] = float(np.mean(dv < -q)) if ok else np.nan
            rec[f"frac_gain_beyond_null_{m}"] = float(np.mean(dv > q)) if ok else np.nan
        out.append(rec)
    return pd.DataFrame(out)
