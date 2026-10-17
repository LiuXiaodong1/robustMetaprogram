"""Pilot P1c / axis A1: does an MP co-vary with the same genes in co-resident TME cells?

Non-malignant cells from the same sample share dissociation and ischaemia time, RNA quality and
capture with the malignant cells, but not malignant programs. Two statistics per MP:

* ``partial_r2``: across samples, the share of variance in the per-sample malignant MP score
  explained by the same gene set's TME score, beyond cancer type (and optionally study), after
  the TME score expected from malignant-derived ambient RNA (spillover, below) is subtracted.
  Malignant score: mean over the sample's malignant cells of MP genes minus expression-matched
  control genes on log2(CPM/10+1). TME score: the same weights on pseudo-bulk profiles (log2 of a
  (sample, cell type) group's mean CPM/10, plus 1; the scale on which ambient RNA mixes linearly),
  centred per cell type across samples and averaged per sample, so TME composition cancels.
  ``slope`` is the coefficient of the (adjusted) TME score; ``spearman`` the rank correlation of
  the malignant and the (adjusted) TME score.
* ``kappa_ratio``: within samples, co-expression (``kappa``) of the MP genes in TME cell types
  divided by that in malignant cells. Near 0: malignant-specific co-expression; near 1 or
  above: the same genes co-vary in TME cells too (technical or niche-wide). Groups larger than
  ``max_kappa_cells`` are subsampled at random (own random stream, so the result does not depend
  on ``n_null`` or ``n_boot``). Not adjusted for spillover: varying ambient fractions add
  co-expression of the leaked genes in TME cells, so a malignant program can raise it.

**Spillover.** Most ambient RNA in a tumour suspension comes from malignant cells, so a malignant
program whose strength varies between samples leaks into the same samples' TME cells and makes
their score follow the malignant one with no TME biology at all. With ``spillover_adjust``:

1. Markers: the ``n_markers`` genes with the highest ratio of malignant to TME mean CPM
   (cohort-wide; TME = the TME cell type with the highest mean, so a marker is low in every TME
   type; pseudo-count 1 CPM) among those at least ``min_marker_fold`` times higher, leaving out
   ``exclude_genes`` (default: the genes of all audited sets, so no MP calibrates its own leak).
   Ambient RNA caps this ratio: a TME type holding a fraction ``rho`` of malignant RNA expresses
   every malignant gene at least at ``rho`` times its malignant level, so no gene exceeds about
   ``1 / rho``-fold (toy cohort below, 5% mean ambient: the planted lineage markers reach 11.6-15.8
   fold). A TME type averaging ``rho`` above about ``1 / min_marker_fold``, or one that shares the
   malignant lineage genes (e.g. normal epithelium next to a carcinoma), leaves too few markers.
   With fewer than ``MIN_MARKERS`` markers the adjustment is skipped with a warning
   (``spill_adjusted`` False; ``partial_r2`` and ``slope`` are then the unadjusted statistic, not
   the spillover-adjusted one that criterion (b) of the P1 GO rule is defined on).
2. Per (sample, TME cell type) group, the malignant-derived ambient fraction
   ``rho_hat = median over markers of TME mean / malignant mean`` (CPM, same sample).
3. Expected profile: TME cells hold a fraction ``rho`` of the sample's malignant profile ``m``
   and ``1 - rho`` of their own, so the expected pseudo-bulk expression of gene ``g`` is
   ``log2((1 - rho_hat) n_g + rho_hat m_g + 1)``, with ``n_g`` the cell type's native profile
   (cohort average over samples of ``max(TME mean - rho_hat m_g, 0) / (1 - rho_hat)``). Between
   samples it moves with the leak of malignant programs and with the dilution of the cell type's
   own genes. Weighted like the TME score, centred per cell type and averaged per sample, it gives
   the spill score (``per_sample`` column ``<set>::tme_spill``), which is subtracted from the TME
   score. ``partial_r2_unadjusted`` and ``slope_unadjusted`` use the TME score as observed, on
   the same samples.

Subtraction, not a covariate: the spill score is ``rho`` times a function of the malignant means,
i.e. of the response itself. As a covariate with a free coefficient in the reduced and full models
it absorbs the malignant score wherever a program varies more between samples than ``rho`` does,
and conditioning on it (a common effect of the malignant score and ``rho``) links any
``rho``-dependent TME score to the malignant score. On toy cohorts (16 fresh samples with 80
malignant and 3 x 40 TME cells, 1,000 genes; seeds 0-2) the covariate version removed the genuine
IFN co-variation (0.09-0.20, below null q99 in all three), left stress below q99 once and pushed
the contamination set (TME markers, diluted by malignant ambient) above q99 twice; subtraction kept
IFN (0.63-0.86) and stress (0.80-0.92) above q99, left hypoxia at 0.03-0.04 and never flagged
contamination. On the default toy cohort the CNA set (malignant only, 70% of tumours) drops from
0.60 unadjusted to 0.07. Malignant-only sets (hypoxia, cycle, EMT, CNA) have observed TME scores
that regress on their spill score with slope 0.8-1.3 (medians over seeds), i.e. the prediction is
roughly calibrated. The price: subtraction relies on the absolute size of ``rho_hat``. Native TME
expression of the markers biases it up (toy: 1.0-1.4 times the planted fraction), which
over-subtracts (a pure leak then gets a negative slope, which does not count); a ``rho_hat``
biased down leaves part of the leak in the TME score.

**Sensitivity to rho_hat.** The random-set null does not absorb a leak left in the TME score (see
Null), so a ``rho_hat`` biased down can lift a malignant-only program above q99. With
``rho_sensitivity`` (default 1.5) the adjusted statistic and its null are computed a second time
with every ``rho_hat`` multiplied by that factor (capped at ``RHO_MAX``), on the same samples,
random sets and control genes: ``partial_r2_rho``, ``slope_rho``, ``null_q95_rho`` and
``null_q99_rho`` (NaN when spillover is not adjusted or ``rho_sensitivity`` is None or 1). An
association that exceeds q99 only because ``rho_hat`` is too small falls back at the larger
``rho_hat``; co-variation that both compartments share natively hardly moves. Toy cohort (16 fresh
samples, 1,500 genes, seed 0, 200 null sets), with every ``rho_hat`` multiplied by 0.7 to mimic an
estimate biased down: CNA_7q (malignant only) rises from 0.15 to 0.41 against null q99 0.27, while
at 1.5 times that ``rho_hat`` it is 0.11 against ``null_q99_rho`` 0.23; stress (0.96) and IFN
(0.85) stay above both. The re-run covers a ``rho_hat`` biased down by up to about
``1 / rho_sensitivity`` (with 0.5 times ``rho_hat``, CNA_7q is still 0.37 against 0.27 at 1.5
times that). The price is over-subtraction where ``rho_hat`` is already biased up, which can push
weak genuine co-variation of a program that varies strongly in malignant cells below its null.

**Null.** ``null_q95``, ``null_q99`` and ``null_p`` come from ``n_null`` random gene sets matched
gene by gene on expression (``n_bins`` bins of malignant mean log expression) and, with
``null_match="expression+specificity"``, malignant specificity (log ratio of malignant to TME mean
CPM as in step 1; ``n_spec_bins`` quantile bins within each expression bin), drawn without
replacement and avoiding the set's own genes (falling back to the expression bin when the joint
bin is used up). Each null set gets its own control genes and exactly the same spillover
subtraction. Sample-level effects that move every gene set (compositional shifts when a strong
program takes up reads, protocol-specific capture, nuisance genes among the controls) inflate the
null too, so they cannot make one MP look TME-driven. Matching on specificity gives the null genes
about the same per-gene exposure to spillover as the MP's genes, but it does not make the null
absorb an error of the subtraction: a random set has none of the between-sample variation of a real
program in malignant cells, so a leak left over by a ``rho_hat`` biased down barely correlates its
two scores, whereas for a program that varies strongly between samples it does. In the toy run
above (``rho_hat`` times 0.7), CNA_7q's partial R^2 rose from 0.15 to 0.41 but its null q99 only
from 0.24 to 0.27. Finer specificity bins do not address this (they refine per-gene exposure, not
program-level variation; with 20 bins CNA_7q is still above null q99); hence the
``rho_sensitivity`` re-run.
Partial R^2 with ``slope <= 0`` counts as 0, for null sets and for ``null_p`` =
(1 + #{null >= observed}) / (1 + number of null sets).

**Use in the P1 GO rule** (``P1-GO v1.1`` in :mod:`mpaudit.report`, criterion b): only a positive
association (``slope > 0``) with ``partial_r2 > null_q99`` counts, and no fixed R^2 threshold.
A shared technical or niche-wide factor moves both compartments the same way; a negative
association has no nuisance reading (competition for reads, over-subtracted spillover). A fixed
bar is not interpretable: on the default toy cohort (16 samples, 4 under the fixed protocol) the
random-set null has q95 0.5-0.8 and q99 0.7-0.9, because protocol-specific capture and stress
genes among the control genes move every set in both compartments; with 16 fresh samples and
1,500 genes q95 is 0.1-0.45. q99 rather than q95 because the rule screens dozens of MPs and the
random-set null is approximate. IFN response and hypoxia can be genuinely TME-shared biology, so
MPs on the pre-declared expected-TME-shared list cannot count through this criterion. The rule is
defined on the spillover-adjusted statistic, i.e. on rows with ``spill_adjusted`` True. The
``*_rho`` columns test whether an association depends on the size of ``rho_hat``: an MP above
``null_q99`` at the estimated ``rho_hat`` that fails ``slope_rho > 0`` and
``partial_r2_rho > null_q99_rho`` at the scaled one is consistent with an under-subtracted leak
(criterion (b) in :mod:`mpaudit.report` requires both wherever the re-run exists).

**Caveats.** ``rho_hat`` comes from the same cells as the scores (shared noise). It assumes the
markers are not expressed natively by TME cells and that the ambient profile is the sample's mean
malignant profile; ambient RNA released by dying or lysed cells can differ (e.g. more nuclear or
mitochondrial transcripts), and markers it under-represents bias ``rho_hat`` down. Ambient RNA
from other TME cells, malignant-TME doublets labelled as TME cells and TME-derived ambient RNA in
malignant cells are not modelled here (the decoy arms condition on contamination estimates). The
random-set null lacks the large between-sample variance of a real program, so it calibrates shared
sample-level effects, not every way a strong program can correlate with noise.

Scores are linear in (log) group-mean profiles, which are computed once, as are the expected
spillover profiles (a second set at the scaled ``rho_hat`` for ``rho_sensitivity``); the
random-set null is therefore cheap even at atlas scale.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .coexpr import kappa_with_null
from .data import Keys, get_counts, log2_cpm10, malignant_mask

PER_SET_COLUMNS = ["mp", "n_samples", "partial_r2", "slope", "partial_r2_lo", "partial_r2_hi", "null_q95",
                   "null_q99", "null_p", "partial_r2_unadjusted", "slope_unadjusted", "spill_adjusted",
                   "partial_r2_rho", "slope_rho", "null_q95_rho", "null_q99_rho", "spearman",
                   "kappa_ratio_median", "kappa_ratio_n"]
MIN_MARKERS = 5  # fewer malignant-specific markers: spillover is not adjusted (with a warning)
RHO_MAX = 0.99  # cap on rho_hat (also after rho_sensitivity scaling); < 1 keeps the native profile finite
NULL_MATCH = ("expression+specificity", "expression")
_PSEUDO_CPM10 = 0.1  # 1 CPM on the CPM/10 scale


def _design(df: pd.DataFrame, cols: list[str]) -> np.ndarray:
    parts = [np.ones((len(df), 1))]
    for c in cols:
        d = pd.get_dummies(df[c].astype(str), drop_first=True, dtype=float)
        if d.shape[1]:
            parts.append(d.to_numpy())
    return np.hstack(parts)


def _fit(y: np.ndarray, x: np.ndarray, Z: np.ndarray) -> tuple[float, float]:
    """Partial R^2 of ``x`` given ``Z``, and the slope of ``x``."""
    def sse(M):
        beta, *_ = np.linalg.lstsq(M, y, rcond=None)
        r = y - M @ beta
        return float(r @ r), beta
    red, _ = sse(Z)
    full, beta = sse(np.hstack([Z, x[:, None]]))
    return ((red - full) / red if red > 0 else np.nan), float(beta[-1])


def _group_means(X: sp.csr_matrix, groups: list[np.ndarray]) -> np.ndarray:
    """Mean row of ``X`` over each group of row indices (groups x genes, dense)."""
    if not groups:
        return np.zeros((0, X.shape[1]))
    rows = np.concatenate([np.full(idx.size, i) for i, idx in enumerate(groups)])
    vals = np.concatenate([np.full(idx.size, 1.0 / idx.size) for idx in groups])
    G = sp.csr_matrix((vals, (rows, np.concatenate(groups))), shape=(len(groups), X.shape[0]))
    return np.asarray((G @ X).todense())


def _bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    """Equal-count bins of ``values`` by rank (ties broken by position)."""
    order = np.argsort(values, kind="stable")
    out = np.empty(values.size, dtype=int)
    out[order] = np.floor(np.arange(values.size) * n_bins / values.size).astype(int)
    return out


def _pools(labels: np.ndarray) -> dict[int, np.ndarray]:
    return {int(b): np.flatnonzero(labels == b) for b in np.unique(labels)}


def _control_genes(idx: np.ndarray, bins: np.ndarray, pools: dict[int, np.ndarray], rng: np.random.Generator,
                   n_ctrl: int) -> np.ndarray:
    """:func:`mpaudit.scoring.control_genes` with precomputed bins (same draws for the same ``rng`` state)."""
    in_set = np.zeros(bins.size, bool)
    in_set[idx] = True
    ctrl = []
    for g in idx:
        pool = pools[bins[g]]
        pool = pool[~in_set[pool]]
        if pool.size:
            ctrl.append(rng.choice(pool, size=min(n_ctrl, pool.size), replace=False))
    return np.concatenate(ctrl) if ctrl else np.array([], dtype=int)


def _matched_random_set(idx: np.ndarray, match: np.ndarray, match_pools: dict[int, np.ndarray], bins: np.ndarray,
                        bin_pools: dict[int, np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """One random gene per gene of ``idx`` from its matching bin, without replacement, avoiding ``idx``."""
    taken = np.zeros(bins.size, bool)
    taken[idx] = True
    out = []
    for g in idx:
        pool = match_pools[match[g]]
        pool = pool[~taken[pool]]
        if not pool.size:
            pool = bin_pools[bins[g]]
            pool = pool[~taken[pool]]
        if pool.size:
            pick = int(rng.choice(pool))
            taken[pick] = True
            out.append(pick)
    return np.sort(np.array(out, dtype=int))


def _rho_hat(lin_mal: np.ndarray, lin_tme: np.ndarray, t_sample: np.ndarray, markers: np.ndarray) -> np.ndarray:
    """``rho_hat`` per TME group: median over ``markers`` of TME / malignant mean (module docstring, step 2).

    ``lin_mal`` (samples x genes) and ``lin_tme`` (TME groups x genes) are mean CPM/10 profiles;
    ``t_sample`` indexes the malignant row of each TME group. NaN for a group whose sample expresses
    none of the markers, otherwise clipped to ``[0, RHO_MAX]``.
    """
    m = lin_mal[t_sample[:, None], markers[None, :]]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(m > 0, lin_tme[:, markers] / m, np.nan)
    rho = np.full(len(t_sample), np.nan)
    ok = np.isfinite(ratio).any(axis=1)
    rho[ok] = np.clip(np.nanmedian(ratio[ok], axis=1), 0.0, RHO_MAX)
    return rho


def _expected_profiles(lin_mal: np.ndarray, lin_tme: np.ndarray, t_sample: np.ndarray, t_type: np.ndarray,
                       rho: np.ndarray) -> np.ndarray:
    """Expected log2(CPM/10+1) pseudo-bulk profile of each TME group for a given (finite) ``rho`` (step 3).

    The native profile of a cell type is re-estimated from ``rho``, so the profiles are not linear in
    ``rho``; a rescaled ``rho`` needs one more call, nothing else.
    """
    expected = np.empty_like(lin_tme)
    for ct in np.unique(t_type):  # one cell type at a time keeps memory at one extra groups x genes array
        sel = np.flatnonzero(t_type == ct)
        r = rho[sel][:, None]
        m = lin_mal[t_sample[sel]]
        native = (np.clip(lin_tme[sel] - r * m, 0.0, None) / (1.0 - r)).mean(axis=0)
        expected[sel] = np.log2((1.0 - r) * native + r * m + 1.0)
    return expected


def tme_covariation(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    min_cells: int = 20,
    covariates: tuple[str, ...] | None = None,
    spillover_adjust: bool = True,
    null_match: str = "expression+specificity",
    n_boot: int = 200,
    n_null: int = 500,
    n_markers: int = 50,
    min_marker_fold: float = 8.0,
    exclude_genes: list[str] | None = None,
    rho_sensitivity: float | None = 1.5,
    kappa_ratio: bool = True,
    max_kappa_cells: int = 2000,
    n_bins: int = 30,
    n_spec_bins: int = 5,
    n_ctrl: int = 100,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return ``(per_set, per_sample)`` tables for the A1 axis (see the module docstring).

    ``per_set`` has the columns ``PER_SET_COLUMNS``. ``per_sample`` has one row per sample with at
    least ``min_cells`` malignant cells: the covariates, ``rho_hat::<cell type>`` (when spillover is
    adjusted; unscaled) and per gene set ``<set>::malignant``, ``<set>::tme`` (as observed) and
    ``<set>::tme_spill`` (the spill score; the adjusted TME score is ``tme - tme_spill``). Random
    draws for a set come from streams seeded by ``(seed, position of the set in gene_sets)``; the
    ``rho_sensitivity`` re-run draws nothing, so the other columns do not depend on it.

    Warns when ``spillover_adjust`` is set but fewer than ``MIN_MARKERS`` markers qualify (the
    statistic is then unadjusted). ``rho_sensitivity`` must be positive or None.
    """
    if keys.cell_type is None:
        raise ValueError("keys.cell_type is required for the TME axis")
    if null_match not in NULL_MATCH:
        raise ValueError(f"null_match must be one of {NULL_MATCH}, got {null_match!r}")
    if rho_sensitivity is not None and not (np.isfinite(rho_sensitivity) and rho_sensitivity > 0):
        raise ValueError(f"rho_sensitivity must be a positive number or None, got {rho_sensitivity!r}")
    sensitivity = rho_sensitivity is not None and rho_sensitivity != 1
    counts = get_counts(adata, counts_layer)
    L = log2_cpm10(counts)
    lin = L.copy()
    lin.data = np.exp2(lin.data) - 1.0  # CPM/10
    genes = np.array(list(map(str, adata.var_names)))
    n_genes = genes.size
    mal = malignant_mask(adata, keys)
    ctype = adata.obs[keys.cell_type].astype(str).to_numpy()
    sample = adata.obs[keys.sample].astype(str).to_numpy()
    tme = ~mal & (adata.obs[keys.cell_type].notna().to_numpy())
    covariates = tuple(c for c in (covariates if covariates is not None else (keys.cancer_type,)) if c)

    # per-sample malignant and per-(sample, cell type) TME groups
    samples, mal_groups, meta_rows = [], [], []
    for s in pd.unique(sample):
        r = np.flatnonzero(mal & (sample == s))
        if r.size >= min_cells:
            samples.append(s)
            mal_groups.append(r)
            meta_rows.append({"sample": s, **{c: str(adata.obs[c].iloc[r[0]]) for c in covariates}})
    meta = pd.DataFrame(meta_rows, columns=["sample", *covariates])
    tme_groups, tme_keys = [], []
    for i, s in enumerate(samples):
        t_rows = np.flatnonzero(tme & (sample == s))
        for ct in pd.unique(ctype[t_rows]):
            r = t_rows[ctype[t_rows] == ct]
            if r.size >= min_cells:
                tme_groups.append(r)
                tme_keys.append((i, ct))
    t_sample = np.array([k[0] for k in tme_keys], dtype=int)
    t_type = np.array([k[1] for k in tme_keys], dtype=object)

    # malignant specificity (cohort-wide): malignant vs the highest TME cell type, linear CPM
    lin_mal_all = np.asarray(lin[np.flatnonzero(mal)].mean(axis=0)).ravel() if mal.any() else np.zeros(n_genes)
    tme_types = pd.unique(ctype[tme])
    lin_tme_max = (np.max([np.asarray(lin[np.flatnonzero(tme & (ctype == ct))].mean(axis=0)).ravel()
                           for ct in tme_types], axis=0) if tme_types.size else np.zeros(n_genes))
    specificity = np.log2((lin_mal_all + _PSEUDO_CPM10) / (lin_tme_max + _PSEUDO_CPM10))

    per_sample = meta.copy()
    lin_tme = _group_means(lin, tme_groups)
    adjust, P, P_rho = False, None, None
    if spillover_adjust and tme_keys:
        excluded = set(exclude_genes) if exclude_genes is not None else {g for gs in gene_sets.values() for g in gs}
        cand = np.flatnonzero(~np.isin(genes, list(excluded)) & (specificity >= np.log2(min_marker_fold)))
        markers = cand[np.argsort(-specificity[cand], kind="stable")[:n_markers]]
        adjust = markers.size >= MIN_MARKERS
        if not adjust:
            warnings.warn(
                f"tme_covariation: only {markers.size} gene(s) outside the excluded genes are >= "
                f"{min_marker_fold:g}-fold higher in malignant than in every TME cell type (MIN_MARKERS = "
                f"{MIN_MARKERS}), so spillover is NOT adjusted: partial_r2 and slope are the unadjusted statistic "
                "(spill_adjusted False), and these sets cannot count through criterion (b) of the P1 GO rule",
                stacklevel=2)
    if adjust:
        lin_mal = _group_means(lin, mal_groups)
        rho = _rho_hat(lin_mal, lin_tme, t_sample, markers)
        keep = np.isfinite(rho)  # groups whose sample has no marker expression cannot be adjusted
        t_sample, t_type, rho, lin_tme = t_sample[keep], t_type[keep], rho[keep], lin_tme[keep]
        P = _expected_profiles(lin_mal, lin_tme, t_sample, t_type, rho)
        if sensitivity:
            rho_scaled = np.minimum(rho * rho_sensitivity, RHO_MAX)
            P_rho = _expected_profiles(lin_mal, lin_tme, t_sample, t_type, rho_scaled)
        for ct in pd.unique(t_type):
            col = np.full(len(samples), np.nan)
            col[t_sample[t_type == ct]] = rho[t_type == ct]
            per_sample[f"rho_hat::{ct}"] = col
    Mm = _group_means(L, mal_groups)
    Mt = np.log2(lin_tme + 1.0)  # pseudo-bulk TME profiles

    ref_mean = np.asarray(L[np.flatnonzero(mal)].mean(axis=0)).ravel() if mal.any() else np.zeros(n_genes)
    bins = _bins(ref_mean, n_bins)
    bin_pools = _pools(bins)
    if null_match == "expression+specificity":
        match = np.empty(n_genes, dtype=int)
        for b, members in bin_pools.items():  # specificity quantiles within each expression bin
            match[members] = b * n_spec_bins + _bins(specificity[members], n_spec_bins)
    else:
        match = bins
    match_pools = _pools(match)
    pos = {g: i for i, g in enumerate(genes)}
    Z_full = _design(meta, list(covariates)) if len(meta) else np.zeros((0, 1))

    def centre_average(group_scores: np.ndarray) -> np.ndarray:
        """Centre each TME cell type across samples, then average the cell types of each sample."""
        centred = group_scores.copy()
        for ct in np.unique(t_type):
            m = t_type == ct
            centred[m] -= group_scores[m].mean()
        sums = np.bincount(t_sample, weights=centred, minlength=len(samples))
        cnt = np.bincount(t_sample, minlength=len(samples))
        return np.where(cnt > 0, sums / np.maximum(cnt, 1), np.nan)

    def scores(idx: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, ...]:
        """Per sample (set minus control genes): malignant score, observed TME score, spill score and the
        spill score at the scaled ``rho_hat`` (None without ``rho_sensitivity``)."""
        v = np.zeros(n_genes)
        v[idx] += 1.0 / idx.size
        ctrl = _control_genes(idx, bins, bin_pools, rng, n_ctrl)
        if ctrl.size:
            v -= np.bincount(ctrl, minlength=n_genes) / ctrl.size
        x = centre_average(Mt @ v) if t_sample.size else np.full(len(samples), np.nan)
        spill = centre_average(P @ v) if adjust and t_sample.size else np.zeros(len(samples))
        spill_rho = None
        if P_rho is not None:
            spill_rho = centre_average(P_rho @ v) if t_sample.size else np.zeros(len(samples))
        return Mm @ v, x, spill, spill_rho

    def stat(y: np.ndarray, x: np.ndarray, spill: np.ndarray) -> tuple[float, float, np.ndarray, np.ndarray]:
        """Partial R^2 and slope of the adjusted TME score, the samples used and that score."""
        xa = x - spill
        ok = np.isfinite(xa) & np.isfinite(y)
        if ok.sum() <= Z_full.shape[1] + 2:
            return np.nan, np.nan, ok, xa
        return (*_fit(y[ok], xa[ok], Z_full[ok]), ok, xa)

    out_rows, set_columns = [], {}
    for i, (name, genes_i) in enumerate(gene_sets.items()):
        rng_sets, rng_boot, rng_kappa = map(np.random.default_rng, np.random.SeedSequence([seed, i]).spawn(3))
        idx = np.array(sorted({pos[g] for g in genes_i if g in pos}), dtype=int)
        rec = dict.fromkeys(PER_SET_COLUMNS, np.nan)
        rec.update(mp=name, n_samples=0, spill_adjusted=adjust, kappa_ratio_n=0)
        if idx.size >= 5 and len(samples):
            y, x, spill, spill_rho = scores(idx, rng_sets)
            set_columns.update({f"{name}::malignant": y, f"{name}::tme": x,
                                f"{name}::tme_spill": spill if adjust else np.full(len(samples), np.nan)})
            r2, slope, ok, xa = stat(y, x, spill)
            if np.isfinite(r2):
                yv, xv, Zv = y[ok], xa[ok], Z_full[ok]
                rec.update(n_samples=int(ok.sum()), partial_r2=r2, slope=slope)
                rec["partial_r2_unadjusted"], rec["slope_unadjusted"] = _fit(yv, x[ok], Zv)
                rec["spearman"] = float(pd.Series(yv).corr(pd.Series(xv), method="spearman"))
                if spill_rho is not None:
                    rec["partial_r2_rho"], rec["slope_rho"], *_ = stat(y, x, spill_rho)
                boots = []
                for _ in range(n_boot):
                    b = rng_boot.integers(0, yv.size, yv.size)
                    Zb = _design(meta.loc[ok].iloc[b], list(covariates))
                    if np.linalg.matrix_rank(Zb) < Zb.shape[1] or b.size <= Zb.shape[1] + 2:
                        continue
                    boots.append(_fit(yv[b], xv[b], Zb)[0])
                boots = np.array([v for v in boots if np.isfinite(v)])
                if boots.size:
                    rec["partial_r2_lo"], rec["partial_r2_hi"] = (float(q) for q in np.quantile(boots, [0.05, 0.95]))
                nulls, nulls_rho = [], []  # the scaled-rho_hat null reuses each random set and its control genes
                for _ in range(n_null):
                    ridx = _matched_random_set(idx, match, match_pools, bins, bin_pools, rng_sets)
                    if ridx.size < 5:
                        continue
                    n_y, n_x, n_spill, n_spill_rho = scores(ridx, rng_sets)
                    for spill_i, acc in ((n_spill, nulls), (n_spill_rho, nulls_rho)):
                        if spill_i is None:
                            continue
                        n_r2, n_slope, *_ = stat(n_y, n_x, spill_i)
                        if np.isfinite(n_r2):
                            acc.append(n_r2 if n_slope > 0 else 0.0)
                if nulls:
                    nulls = np.array(nulls)
                    observed = r2 if slope > 0 else 0.0
                    rec["null_q95"], rec["null_q99"] = (float(q) for q in np.quantile(nulls, [0.95, 0.99]))
                    rec["null_p"] = float((1 + np.sum(nulls >= observed)) / (1 + nulls.size))
                if nulls_rho:
                    q95, q99 = np.quantile(nulls_rho, [0.95, 0.99])
                    rec["null_q95_rho"], rec["null_q99_rho"] = float(q95), float(q99)
        if kappa_ratio:
            ratios = []
            for s in samples:
                m_rows = np.flatnonzero(mal & (sample == s))
                t_rows = np.flatnonzero(tme & (sample == s))
                if m_rows.size < min_cells or t_rows.size < min_cells:
                    continue
                if m_rows.size > max_kappa_cells:
                    m_rows = np.sort(rng_kappa.choice(m_rows, max_kappa_cells, replace=False))
                km = kappa_with_null(counts[m_rows], idx, n_null=0)["kappa"]
                kt, wt = [], []
                for ct in pd.unique(ctype[t_rows]):
                    r = t_rows[ctype[t_rows] == ct]
                    if r.size >= min_cells:
                        if r.size > max_kappa_cells:
                            r = np.sort(rng_kappa.choice(r, max_kappa_cells, replace=False))
                        v = kappa_with_null(counts[r], idx, n_null=0)["kappa"]
                        if np.isfinite(v):
                            kt.append(v)
                            wt.append(r.size)
                if np.isfinite(km) and km > 0 and kt:
                    ratios.append(float(np.average(kt, weights=wt)) / km)
            rec["kappa_ratio_median"] = float(np.median(ratios)) if ratios else np.nan
            rec["kappa_ratio_n"] = len(ratios)
        out_rows.append(rec)
    per_sample = pd.concat([per_sample, pd.DataFrame(set_columns, index=per_sample.index)], axis=1)
    return pd.DataFrame(out_rows, columns=PER_SET_COLUMNS), per_sample
