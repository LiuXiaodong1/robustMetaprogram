"""Axis A4: does an MP change between two protocols applied to the same tumour?

For tumours profiled with a reference protocol (e.g. fresh dissociation) and a test protocol (e.g.
fixed/Flex or nuclei), each MP is compared between the two within each tumour, and the change is
set against a null. Changes are ``test - reference``: a loss is a drop under the test protocol.
Loss and gain are tested separately, each one-sided.

**Co-expression, Δlog κ (the primary statistic, criterion (c) of the P1 rule).** Per malignant
sample, ``kappa`` (:func:`mpaudit.coexpr.kappa`, on at most ``max_cells`` random cells) of the MP's
genes, restricted to the tumour's shared panel: genes detected (count > 0) in >= ``min_detect`` of
the tumour's malignant cells under *each* of its protocols. Panel retention ``iota`` is the fraction
of the MP's genes in the panel; below ``min_iota`` (0.5) the MP's kappa in that tumour is NaN. Per
paired tumour, ``Δlog κ = log κ_test - log κ_ref``, with κ of a protocol the mean over the tumour's
samples under it, floored at ``kappa_floor`` (0.01) before the log. Why κ: it is a correlation, so a
gene-wise capture efficiency (probe efficiency in Flex, transcript-length bias in nuclei) multiplies
a gene's latent expression without changing its correlations; the panel keeps genes absent from one
protocol (no probe, not detected) from changing the set being compared. Caveats: κ also moves with
the fraction of cells in the state (a selection shift, e.g. fragile cells lost in one protocol,
changes κ as well as intensity); its estimation noise grows as depth falls; a program absent under
one protocol has κ near 0, so its Δlog κ is set by the floor rather than measured.

**Intensity and prevalence, with protocol-specific control genes (informational).** Scores follow
:mod:`mpaudit.scoring` (set mean minus control mean of log2(CPM/10+1), controls from 30 equal-size
bins of mean expression), except that a set gene's controls are all other genes of its bin (the
limit of drawing many: no control-sampling noise, which would differ between the protocols) and
that bins, controls and the prevalence threshold (the ``high_quantile`` of the scores) are computed
separately within each protocol, on the malignant cells of the *paired* tumours under it. With bins
pooled across protocols, gene-wise capture moves set and control genes differently and can look like
loss (on the toy, the 12 cycle genes captured worst under the fixed protocol lose intensity beyond
the replicate null in 3 of 4 tumours with pooled controls, 1 of 4 with protocol-specific ones).
Each protocol's pool is restricted to the paired tumours because pools of different tumours let
tumour composition masquerade as a protocol shift: with all of each protocol's samples on the toy,
EMT (more prevalent under the fixed protocol) came out as a loss in 3 of 4 tumours. The price of
protocol-specific bins: every gene is matched on its own level within the protocol, so any shift
common to all paired tumours, technical or biological, is absorbed. Intensity and prevalence
therefore measure a tumour's shift *relative to the other paired tumours*, and a cohort-wide loss
shows only in κ. On the toy, dissociation stress (absent under the fixed protocol) loses κ in every
paired tumour while its mean intensity change is about 0; the tumours whose intensity drops beyond
the null are those whose fresh piece was more stressed than the other paired tumours'. Intensity
uses all of the MP's genes in ``var_names``, not the panel.

**Nulls.** ``replicates``: same-protocol pairs of samples of one tumour (pieces from different
regions), used when an MP has >= ``min_pairs`` (3) such pairs. A change counts as beyond the null
when it exceeds the 95th percentile of the replicates' absolute changes, in the tested direction
(with 3 pairs, close to the largest replicate change). This null contains regional heterogeneity,
but most public protocol twins have no replicate pieces.
Otherwise ``random_sets``: per paired tumour, ``n_null_sets`` random gene sets drawn from the
tumour's panel, matched gene by gene to the MP's panel genes on mean-expression and detection bins
in the tumour's reference-protocol cells (:func:`mpaudit.coexpr.matched_random_sets`), and evaluated
like the MP in the same samples (κ on the same cells, intensity with the same protocol-specific
control machinery). A tumour's change counts as a loss when it is below the 5th percentile of that
tumour's random-set changes (a gain: above the 95th; at least ``MIN_RANDOM_SETS`` finite changes
are needed, otherwise the tumour is not tested). This null captures shifts shared by all genes
(capture, depth, estimation noise of κ) but not regional heterogeneity: two pieces that differ in a
program's activity look like a protocol effect, so for a program whose activity varies between
pieces (hypoxia on the toy) intensity changes beyond this null are expected without any protocol
effect. Random sets have κ near 0, often under the floor, so their Δlog κ spreads widely on the log
scale: the null is wide for κ (only large losses pass) and, when nearly all random sets are floored
in both protocols, collapses to 0 (then reported as untestable: NaN). Prevalence has no random-set
null (its threshold would need every random set scored on the whole pool).

**Data.** Oliveira et al. (Nat Genet 2025) is NOT a fresh-vs-fixed twin: its Chromium Flex, Visium HD
and Xenium data come from the same FFPE blocks and no fresh-dissociated malignant cells were found,
so it offers a fixed-dissociated versus in situ comparison, not a fresh-vs-fixed malignant contrast.

Effects on the toy cohort (:mod:`mpaudit.sim`) are true by construction; they show that the code
responds as designed, not how often real MPs change between protocols.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .coexpr import kappa, matched_random_sets
from .data import Keys, get_counts, log2_cpm10, malignant_mask

KAPPA_FLOOR = 0.01
MIN_PAIRS = 3
MIN_RANDOM_SETS = 20  # finite random-set changes needed for a tumour's 5th/95th percentiles
N_BINS = 30  # expression bins for control genes, as in mpaudit.scoring
METRICS = ("kappa", "intensity", "prevalence")


@dataclass
class _Sample:
    name: str
    tumour: str
    protocol: str
    rows: np.ndarray  # the sample's analysed malignant cells (rows of _Cohort.counts)
    kappa_rows: np.ndarray  # at most max_cells of them, for kappa


@dataclass
class _Cohort:
    """Malignant cells of the two protocols, prepared once for the MPs and the random-set null."""

    reference: str
    test: str
    genes: np.ndarray
    counts: sp.csr_matrix  # analysed cells x genes
    L: sp.csr_matrix  # log2(CPM/10 + 1); CPM over all genes
    protocol: np.ndarray  # per analysed cell
    samples: list[_Sample]
    panel: dict[str, np.ndarray]  # tumour -> genes detected in every protocol the tumour has
    paired: list[str]
    pool: dict[str, np.ndarray]  # protocol -> rows defining control bins and prevalence thresholds
    bins: dict[str, np.ndarray]  # protocol -> expression bin of each gene over its pool
    _kappa_data: dict = field(default_factory=dict)  # sample -> (CSC counts of its kappa cells, depth)

    def set_kappa(self, s: _Sample, idx: np.ndarray) -> float:
        """``kappa`` of genes ``idx`` in the sample's kappa cells (depth = total counts over all genes)."""
        if s.name not in self._kappa_data:
            Y = self.counts[s.kappa_rows]
            self._kappa_data[s.name] = (Y.tocsc(), np.asarray(Y.sum(axis=1)).ravel())
        Y, depth = self._kappa_data[s.name]
        return kappa(Y[:, idx].toarray(), depth)

    def score_weights(self, protocol: str, idx: np.ndarray) -> np.ndarray:
        """Gene weights ``v`` such that ``L @ v`` is the set mean minus the control mean under ``protocol``.

        Each set gene's controls are all genes outside the set in its expression bin, weighted so
        that every set gene with a non-empty bin counts equally.
        """
        bins = self.bins[protocol]
        in_set = np.zeros(bins.size, bool)
        in_set[idx] = True
        pool = np.bincount(bins[~in_set], minlength=N_BINS)
        k = np.bincount(bins[idx], minlength=N_BINS) * (pool > 0)
        v = np.where(in_set, 0.0, -(k / np.maximum(pool, 1))[bins] / max(k.sum(), 1))
        v[idx] = 1.0 / idx.size
        return v


def _rngs(seed: int) -> tuple[np.random.Generator, np.random.Generator]:
    """Streams for the cell subsample (κ) and the random gene sets."""
    cells, null = np.random.SeedSequence(seed).spawn(2)
    return np.random.default_rng(cells), np.random.default_rng(null)


def _expression_bins(mean: np.ndarray) -> np.ndarray:
    """``N_BINS`` equal-size bins of genes ranked by ``mean`` (as in :func:`mpaudit.scoring.control_genes`)."""
    bins = np.empty(mean.size, dtype=int)
    bins[np.argsort(mean, kind="stable")] = np.arange(mean.size) * N_BINS // max(mean.size, 1)
    return bins


def _detected(X: sp.csr_matrix) -> np.ndarray:
    return np.asarray((X > 0).sum(axis=0)).ravel() / max(X.shape[0], 1)


def _prepare(adata, keys: Keys, reference: str, test: str, *, counts_layer: str | None, min_cells: int,
             max_cells: int, min_detect: float, rng: np.random.Generator) -> _Cohort:
    if keys.protocol is None:
        raise ValueError("keys.protocol is required")
    if reference == test:
        raise ValueError("reference and test protocols must differ")
    obs = adata.obs
    protocol = obs[keys.protocol].astype(str).to_numpy()
    sample = obs[keys.sample].astype(str).to_numpy()
    tumour = obs[keys.tumour_col()].astype(str).to_numpy()
    counts = get_counts(adata, counts_layer)
    cells = (malignant_mask(adata, keys) & np.isin(protocol, [reference, test])
             & (np.asarray(counts.sum(axis=1)).ravel() > 0))
    n_per_sample = pd.Series(sample[cells]).value_counts()
    cells &= np.isin(sample, n_per_sample.index[n_per_sample >= min_cells])
    rows = np.flatnonzero(cells)
    counts = counts[rows]
    protocol, sample, tumour = protocol[rows], sample[rows], tumour[rows]

    samples = []
    for s in pd.unique(sample):
        r = np.flatnonzero(sample == s)
        if np.unique(protocol[r]).size > 1:
            raise ValueError(f"sample {s!r} has malignant cells under both protocols")
        k = np.sort(rng.choice(r, max_cells, replace=False)) if r.size > max_cells else r
        samples.append(_Sample(s, tumour[r[0]], protocol[r[0]], r, k))

    panel, paired = {}, []
    for t in pd.unique(tumour):
        protocols = {s.protocol for s in samples if s.tumour == t}
        ok = np.ones(counts.shape[1], bool)
        for p in protocols:
            ok &= _detected(counts[np.flatnonzero((tumour == t) & (protocol == p))]) >= min_detect
        panel[t] = ok
        if protocols == {reference, test}:
            paired.append(t)

    L = log2_cpm10(counts)
    pool, bins = {}, {}
    for p in (reference, test):
        in_p = protocol == p
        from_paired = in_p & np.isin(tumour, paired)
        pool[p] = np.flatnonzero(from_paired if from_paired.any() else in_p)
        mean = np.asarray(L[pool[p]].mean(axis=0)).ravel() if pool[p].size else np.zeros(counts.shape[1])
        bins[p] = _expression_bins(mean)
    return _Cohort(reference, test, adata.var_names.astype(str).to_numpy(), counts, L, protocol, samples,
                   panel, paired, pool, bins)


def _levels(cohort: _Cohort, gene_sets: dict[str, list[str]], *, min_iota: float, min_genes: int,
            high_quantile: float) -> pd.DataFrame:
    pos = {g: i for i, g in enumerate(cohort.genes)}
    members = {n: list(dict.fromkeys(map(str, genes))) for n, genes in gene_sets.items()}
    idx = {n: np.array([pos[g] for g in genes if g in pos], dtype=int) for n, genes in members.items()}
    n_cells = cohort.counts.shape[0]
    score = {n: np.full(n_cells, np.nan) for n in gene_sets}
    thr = {}
    for p in (cohort.reference, cohort.test):
        rows = np.flatnonzero(cohort.protocol == p)
        if rows.size == 0:
            continue
        Lp = cohort.L[rows]
        for n in gene_sets:
            if idx[n].size >= min_genes:
                score[n][rows] = Lp @ cohort.score_weights(p, idx[n])
                thr[n, p] = np.quantile(score[n][cohort.pool[p]], high_quantile)

    out = []
    for s in cohort.samples:
        panel = cohort.panel[s.tumour]
        for n in gene_sets:
            target = idx[n][panel[idx[n]]]
            iota = target.size / len(members[n]) if members[n] else np.nan
            v, t = score[n][s.rows], thr.get((n, s.protocol), np.nan)
            out.append({"sample": s.name, "tumour": s.tumour, "protocol": s.protocol, "mp": n,
                        "n_cells": int(s.rows.size), "n_genes": int(target.size), "iota": float(iota),
                        "kappa": cohort.set_kappa(s, target) if iota >= min_iota and target.size >= min_genes
                        else np.nan,
                        "intensity": float(np.mean(v)),
                        "prevalence": float(np.mean(v > t)) if np.isfinite(t) else np.nan})
    return pd.DataFrame(out, columns=["sample", "tumour", "protocol", "mp", "n_cells", "n_genes", "iota",
                                      "kappa", "intensity", "prevalence"])


def _random_set_levels(cohort: _Cohort, gene_sets: dict[str, list[str]], mps: list[str], n_sets: int, *,
                       min_iota: float, min_genes: int, rng: np.random.Generator) -> pd.DataFrame:
    """κ and intensity of matched random gene sets, per sample of each paired tumour."""
    pos = {g: i for i, g in enumerate(cohort.genes)}
    out = []
    for t in cohort.paired:
        ss = [s for s in cohort.samples if s.tumour == t]
        ref_rows = np.concatenate([s.rows for s in ss if s.protocol == cohort.reference])
        Y = cohort.counts[ref_rows]
        mean_expr = np.asarray(Y.sum(axis=0)).ravel() / max(Y.sum(), 1e-12)
        detect = _detected(Y)
        L_s = {s.name: cohort.L[s.rows] for s in ss}
        panel = cohort.panel[t]
        for n in mps:
            genes = list(dict.fromkeys(map(str, gene_sets[n])))
            target = np.array([pos[g] for g in genes if g in pos and panel[pos[g]]], dtype=int)
            if not genes or target.size / len(genes) < min_iota or target.size < min_genes:
                continue
            for j, r in enumerate(matched_random_sets(target, mean_expr, detect, n_sets, rng, candidates=panel)):
                v = {p: cohort.score_weights(p, r) for p in (cohort.reference, cohort.test)}
                for s in ss:
                    out.append({"mp": n, "tumour": t, "sample": s.name, "protocol": s.protocol, "set": j,
                                "kappa": cohort.set_kappa(s, r),
                                "intensity": float(np.mean(L_s[s.name] @ v[s.protocol]))})
    return pd.DataFrame(out, columns=["mp", "tumour", "sample", "protocol", "set", "kappa", "intensity"])


def per_sample_kappa(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    reference: str,
    test: str,
    *,
    counts_layer: str | None = None,
    min_cells: int = 20,
    max_cells: int = 2000,
    min_detect: float = 0.05,
    min_iota: float = 0.5,
    min_genes: int = 5,
    high_quantile: float = 0.9,
    seed: int = 0,
) -> pd.DataFrame:
    """Per malignant sample of the two protocols and MP: κ on the tumour's panel, intensity and prevalence.

    Samples need >= ``min_cells`` malignant cells; see the module docstring for the statistics.

    Columns: ``sample, tumour, protocol, mp, n_cells, n_genes`` (panel genes used for κ), ``iota``,
    ``kappa`` (NaN when ``iota < min_iota`` or fewer than ``min_genes`` panel genes), ``intensity``,
    ``prevalence``. Each sample must have one protocol.
    """
    cohort = _prepare(adata, keys, reference, test, counts_layer=counts_layer, min_cells=min_cells,
                      max_cells=max_cells, min_detect=min_detect, rng=_rngs(seed)[0])
    return _levels(cohort, gene_sets, min_iota=min_iota, min_genes=min_genes, high_quantile=high_quantile)


def _transform(metric: str, x: pd.Series | np.ndarray, kappa_floor: float) -> np.ndarray:
    """κ is compared on the log scale after flooring; intensity and prevalence as they are."""
    x = np.asarray(x, dtype=float)
    return np.log(np.maximum(x, kappa_floor)) if metric == "kappa" else x


def _delta_name(metric: str) -> str:
    return "delta_log_kappa" if metric == "kappa" else f"delta_{metric}"


def _replicate_deltas(levels: pd.DataFrame, metrics: list[str], kappa_floor: float) -> pd.DataFrame:
    """Changes between every two samples of one tumour under the same protocol, per MP."""
    out = []
    for (mp, _, _), g in levels.groupby(["mp", "tumour", "protocol"], sort=False):
        v = {m: _transform(m, g[m], kappa_floor) for m in metrics}
        for i, j in itertools.combinations(range(len(g)), 2):
            out.append({"mp": mp, **{m: v[m][j] - v[m][i] for m in metrics}})
    return pd.DataFrame(out, columns=["mp", *metrics])


def _paired_deltas(df: pd.DataFrame, by: list[str], metrics: list[str], reference: str, test: str,
                   kappa_floor: float) -> pd.DataFrame:
    """``f(mean under test) - f(mean under reference)`` per group ``by`` that has both protocols."""
    agg = df.groupby([*by, "protocol"], sort=False)[metrics].mean().reset_index()
    ref = agg[agg["protocol"] == reference].set_index(by)[metrics]
    tst = agg[agg["protocol"] == test].set_index(by)[metrics]
    both = ref.index.intersection(tst.index, sort=False)
    d = pd.DataFrame({m: _transform(m, tst.loc[both, m], kappa_floor) - _transform(m, ref.loc[both, m], kappa_floor)
                      for m in metrics}, index=both)
    for m in metrics:
        d[f"{m}_reference"], d[f"{m}_test"] = ref.loc[both, m].to_numpy(), tst.loc[both, m].to_numpy()
    return d


def _null_types(levels: pd.DataFrame, reps: pd.DataFrame, min_pairs: int, random_sets: bool) -> pd.DataFrame:
    """Per MP: ``null_type`` and ``n_null_pairs`` (same-protocol replicate pairs)."""
    mps = pd.unique(levels["mp"])
    n = reps["mp"].value_counts().reindex(mps, fill_value=0).astype(int)
    kind = np.where(n >= min_pairs, "replicates", "random_sets" if random_sets else "none")
    return pd.DataFrame({"null_type": kind, "n_null_pairs": n.to_numpy()}, index=pd.Index(mps, name="mp"))


def _shifts(levels: pd.DataFrame, reference: str, test: str, null_levels: pd.DataFrame | None, min_pairs: int,
            kappa_floor: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(per_tumour, null types)``; see :func:`tumour_shifts`."""
    metrics = [m for m in METRICS if m in levels]
    reps = _replicate_deltas(levels, metrics, kappa_floor)
    types = _null_types(levels, reps, min_pairs, null_levels is not None)
    d = _paired_deltas(levels, ["mp", "tumour"], metrics, reference, test, kappa_floor)

    lo = pd.DataFrame(np.nan, index=d.index, columns=metrics)
    hi = lo.copy()
    by_replicates = types.index[types["null_type"] == "replicates"]
    for mp in by_replicates:
        rows = d.index.get_level_values("mp") == mp
        for m in metrics:
            v = np.abs(reps.loc[reps["mp"] == mp, m].to_numpy(float))
            v = v[np.isfinite(v)]
            if v.size >= min_pairs:
                q = float(np.quantile(v, 0.95))
                lo.loc[rows, m], hi.loc[rows, m] = -q, q
    if null_levels is not None and not null_levels.empty:
        nm = [m for m in metrics if m in null_levels]
        nd = _paired_deltas(null_levels[~null_levels["mp"].isin(by_replicates)], ["mp", "tumour", "set"], nm,
                            reference, test, kappa_floor)[nm]
        g = nd.groupby(level=["mp", "tumour"])
        enough = g.count() >= MIN_RANDOM_SETS
        q05, q95 = g.quantile(0.05).where(enough), g.quantile(0.95).where(enough)
        collapsed = q05 >= q95  # e.g. every random set under the kappa floor in both protocols
        q05, q95 = q05.mask(collapsed), q95.mask(collapsed)
        common = d.index.intersection(q05.index, sort=False)
        lo.loc[common, nm] = q05.loc[common, nm].to_numpy()
        hi.loc[common, nm] = q95.loc[common, nm].to_numpy()

    out = d.reset_index().rename(columns={m: _delta_name(m) for m in metrics})
    out.insert(2, "null_type", out["mp"].map(types["null_type"]))
    if "iota" in levels:
        iota = levels.groupby(["mp", "tumour"])["iota"].first()
        out.insert(3, "iota", iota.reindex(d.index).to_numpy())
    for m in metrics:
        delta, lo_m, hi_m = d[m].to_numpy(), lo[m].to_numpy(), hi[m].to_numpy()
        tested = np.isfinite(delta) & np.isfinite(lo_m) & np.isfinite(hi_m)
        out[f"null_lo_{m}"], out[f"null_hi_{m}"] = lo_m, hi_m
        with np.errstate(invalid="ignore"):
            out[f"loss_{m}"] = pd.array(np.where(tested, delta < lo_m, None), dtype="boolean")
            out[f"gain_{m}"] = pd.array(np.where(tested, delta > hi_m, None), dtype="boolean")
    return out, types


def tumour_shifts(
    levels: pd.DataFrame,
    reference: str,
    test: str,
    *,
    null_levels: pd.DataFrame | None = None,
    min_pairs: int = MIN_PAIRS,
    kappa_floor: float = KAPPA_FLOOR,
) -> pd.DataFrame:
    """Per MP and paired tumour: the ``test - reference`` change of each statistic, and whether it is beyond the null.

    ``levels`` has one row per sample and MP (``mp, tumour, protocol`` and any of ``kappa``,
    ``intensity``, ``prevalence``; e.g. from :func:`per_sample_kappa`). ``null_levels`` holds random
    gene sets evaluated in the paired tumours (``mp, tumour, protocol, set`` and statistics); it is
    used for MPs with fewer than ``min_pairs`` same-protocol replicate pairs. Columns: ``mp, tumour,
    null_type``, ``iota`` (when in ``levels``) and per statistic the change (``delta_log_kappa``,
    ``delta_<metric>``), the protocol means (``<metric>_reference``, ``<metric>_test``), the null
    bounds ``null_lo_<metric>``, ``null_hi_<metric>`` (replicates: -/+ the 95th percentile of
    absolute replicate changes; random sets: the tumour's 5th/95th percentiles) and
    ``loss_<metric>``, ``gain_<metric>`` (change below/above the bound; NA when the change or the null
    is undefined).
    """
    return _shifts(levels, reference, test, null_levels, min_pairs, kappa_floor)[0]


def _summarise(shifts: pd.DataFrame, types: pd.DataFrame) -> pd.DataFrame:
    metrics = [m for m in METRICS if f"loss_{m}" in shifts]
    columns = ["mp", "null_type", "n_null_pairs", "n_paired_tumours", "iota"]
    for m in metrics:
        columns += [f"mean_{_delta_name(m)}", f"frac_loss_beyond_null_{m}", f"frac_gain_beyond_null_{m}"]
    out = []
    for mp, row in types.iterrows():
        t = shifts[shifts["mp"] == mp]
        rec = {"mp": mp, "null_type": row["null_type"], "n_null_pairs": int(row["n_null_pairs"]),
               "n_paired_tumours": int(t[f"loss_{metrics[0]}"].notna().sum()) if metrics else 0,
               "iota": float(t["iota"].median()) if "iota" in t and len(t) else np.nan}
        for m in metrics:
            delta = t[_delta_name(m)].to_numpy(float)
            loss, gain = t[f"loss_{m}"].dropna(), t[f"gain_{m}"].dropna()
            rec[f"mean_{_delta_name(m)}"] = float(np.nanmean(delta)) if np.isfinite(delta).any() else np.nan
            rec[f"frac_loss_beyond_null_{m}"] = float(loss.mean()) if len(loss) else np.nan
            rec[f"frac_gain_beyond_null_{m}"] = float(gain.mean()) if len(gain) else np.nan
        out.append(rec)
    return pd.DataFrame(out, columns=columns)


def protocol_shift(
    levels: pd.DataFrame,
    reference: str,
    test: str,
    *,
    null_levels: pd.DataFrame | None = None,
    min_pairs: int = MIN_PAIRS,
    kappa_floor: float = KAPPA_FLOOR,
) -> pd.DataFrame:
    """Per MP: mean change and fraction of paired tumours with a loss (gain) beyond the null.

    Arguments as for :func:`tumour_shifts`. Columns: ``mp, null_type`` (``replicates``,
    ``random_sets``, or ``none`` when neither null is available), ``n_null_pairs`` (same-protocol
    replicate pairs), ``n_paired_tumours`` (paired tumours tested for the first statistic present,
    κ when given: change and null defined), ``iota`` (median over paired tumours), and per statistic
    ``mean_delta_log_kappa`` / ``mean_delta_<metric>``, ``frac_loss_beyond_null_<metric>`` and
    ``frac_gain_beyond_null_<metric>`` (over the paired tumours where that statistic was tested).
    """
    return _summarise(*_shifts(levels, reference, test, null_levels, min_pairs, kappa_floor))


def protocol_shift_kappa(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    reference: str,
    test: str,
    *,
    counts_layer: str | None = None,
    min_pairs: int = MIN_PAIRS,
    n_null_sets: int = 100,
    kappa_floor: float = KAPPA_FLOOR,
    min_cells: int = 20,
    max_cells: int = 2000,
    min_detect: float = 0.05,
    min_iota: float = 0.5,
    min_genes: int = 5,
    high_quantile: float = 0.9,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Axis A4 end to end. Returns ``(per_mp, per_tumour, levels)``.

    ``levels`` is :func:`per_sample_kappa` (same arguments, same values); MPs with fewer than
    ``min_pairs`` same-protocol replicate pairs get the random-set null (``n_null_sets`` sets per
    paired tumour, at least ``MIN_RANDOM_SETS``; 0 switches it off and leaves their null ``none``); ``per_tumour`` is
    :func:`tumour_shifts` and ``per_mp`` :func:`protocol_shift`. Run time is dominated by the
    random-set κ: about ``n_null_sets`` x MPs x samples of paired tumours evaluations of
    :func:`mpaudit.coexpr.kappa`.
    """
    rng_cells, rng_null = _rngs(seed)
    cohort = _prepare(adata, keys, reference, test, counts_layer=counts_layer, min_cells=min_cells,
                      max_cells=max_cells, min_detect=min_detect, rng=rng_cells)
    levels = _levels(cohort, gene_sets, min_iota=min_iota, min_genes=min_genes, high_quantile=high_quantile)
    reps = _replicate_deltas(levels, [m for m in METRICS if m in levels], kappa_floor)
    n_pairs = reps["mp"].value_counts()
    need = [mp for mp in gene_sets if n_pairs.get(mp, 0) < min_pairs]
    null_levels = None
    if need and n_null_sets > 0:
        null_levels = _random_set_levels(cohort, gene_sets, need, n_null_sets, min_iota=min_iota,
                                         min_genes=min_genes, rng=rng_null)
    per_tumour, types = _shifts(levels, reference, test, null_levels, min_pairs, kappa_floor)
    return _summarise(per_tumour, types), per_tumour, levels
