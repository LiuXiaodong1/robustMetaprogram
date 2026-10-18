"""Axis A4: does an MP change between two protocols applied to the same tumour?

For tumours profiled with a reference protocol (e.g. fresh dissociation) and a test protocol (e.g.
fixed/Flex or nuclei), each MP is compared between the two within each tumour, and the change is
set against a null. Changes are ``test - reference``: a loss is a drop under the test protocol.
Loss and gain are tested separately, each one-sided. This module implements criterion (c) of the
frozen rule P1-GO v1.3 (:mod:`mpaudit.report`).

**Co-expression, Δlog κ (the primary statistic, criterion (c) of the P1 rule).** Per malignant
sample with >= ``min_cells`` (``MIN_CELLS``, 20) malignant cells, ``kappa`` (:func:`mpaudit.coexpr.kappa`)
of the MP's genes in the tumour's shared panel, computed on the sample's κ cells:

* *Equal cells.* Every sample of a tumour is scored on the same number of cells, ``n_t =
  min(max_cells, the tumour's smallest sample)`` (``MAX_KAPPA_CELLS``, 2000; reported as
  ``n_kappa_cells``), drawn without replacement from each sample's own stream (:func:`_cells_rng`).
  κ falls with the number of cells it is computed on (toy: gene sets without a program had κ about
  0.19 on 200 cells and 0.15 on 60), and a 50-cell test piece against 200-cell reference pieces made
  0.078 of the absent gene sets' tests losses on the toy null instead of 0.044 (design comparison of
  2026-10-17, 6 seeds; numbers from its probe runs, not re-run here).
* *Equal depth.* The κ cells' counts of every sample of a tumour are binomially thinned to the
  tumour's lowest median total count per cell (each count ``y`` becomes a ``Binomial(y, p)`` draw,
  ``p`` = the tumour's lowest ``median_umi`` / the sample's, reported as ``kappa_thin``; stream
  :func:`_depth_rng`; whole-number counts are required). κ's moment estimate falls as counts per cell
  fall, more for weakly co-expressed sets: on the toy null (seed 0, every sample) thinning the same
  cells to a quarter of their depth gave κ x 0.75-0.92 for cycle, hypoxia and stress and x 0.48-0.60
  for random 40-gene sets (median ratio over the 20 samples; run of 2026-10-18). Without equal
  depth, test pieces at a quarter of the reference depth made 0.12 of the toy null's absent pairs and
  0.06 of its carrier pairs losses, and an unplanted private program met criterion (c) (review probe
  of 2026-10-18 on the code before equal depth, 6 seeds; not re-run here);
  ``scripts/calibration/evaluate_c.py`` now runs this case.
* *Panel.* Genes detected (count > 0) in >= ``min_detect`` (``MIN_DETECT``, 0.05) of the tumour's κ
  cells, after thinning, under *each* of its protocols (on the counts κ is computed on, so both
  protocols select genes alike). Panel retention ``iota`` is the fraction of the MP's genes in the
  panel; below ``min_iota`` (``MIN_IOTA``, 0.5), or with fewer than ``min_genes`` (``MIN_GENES``, 5)
  panel genes, the MP's kappa in that tumour is NaN.

Each sample's cells are taken in barcode order, and every per-sample draw (κ cells, thinning,
bootstrap) comes from a stream seeded by the seed and the sample's name (the bootstrap also by the
MP's), so no draw depends on the other samples or MPs of the call or on the cohort's cell order.
Per paired tumour, ``Δlog κ = log κ_test - log κ_ref``, with κ of a protocol the mean over the
tumour's samples under it, floored at ``kappa_floor`` (0.01) before the log. Why κ: its estimand, a
latent correlation, is unchanged by a gene-wise capture efficiency (probe efficiency in Flex,
transcript-length bias in nuclei), which multiplies a gene's latent expression without changing its
correlations, and by depth; the estimate is not (equal depth, above, makes the two sides
comparable). The panel keeps genes absent from one protocol (no probe, not detected) from changing
the set being compared. Caveats: κ also moves with the fraction of cells in the state (a selection
shift, e.g. fragile cells lost in one protocol, changes κ as well as intensity); a gene-wise capture
loss still lowers the affected genes' counts and so their estimated correlations; a program absent
under the test protocol has κ near 0 there, so the size of its Δlog κ is set by the floor rather
than measured (it is still a loss); a program absent from a tumour under both protocols has κ at
background in both, and that tumour is tested like any other (no activity gate, below).

**Intensity and prevalence, with protocol-specific control genes (informational).** Scores follow
:mod:`mpaudit.scoring` (set mean minus control mean of log2(CPM/10+1), controls from 30 equal-size
bins of mean expression), except that a set gene's controls are all other genes of its bin (the
limit of drawing many: no control-sampling noise, which would differ between the protocols) and
that bins, controls and the prevalence threshold (the ``high_quantile`` of the scores) are computed
separately within each protocol, on the malignant cells of the *paired* tumours under it. With bins
pooled across protocols, gene-wise capture moves set and control genes differently and can look like
loss (on the small test cohort of ``tests/test_protocol.py``, 4 paired and 3 replicate tumours, the
12 cycle genes captured worst under the fixed protocol lose intensity beyond the replicate null in
3 of 4 paired tumours with pooled controls, 0 of 4 with protocol-specific
ones and 0 of 4 in κ). Each protocol's pool is restricted to the paired tumours because pools of
different tumours let tumour composition masquerade as a protocol shift: with all of each protocol's
samples (toy run below), EMT (more prevalent under the fixed protocol) comes out as an intensity loss
in 2 of 4 paired tumours, against 1 of 4 with the paired-tumour pools. The price of protocol-specific
bins: every gene is matched on its own level within the protocol, so any shift common to all paired
tumours, technical or biological, is absorbed. Intensity and prevalence therefore measure a tumour's
shift *relative to the other paired tumours*, and a cohort-wide loss shows only in κ. On the toy run,
dissociation stress (absent under the fixed protocol) loses κ in every paired tumour while its mean
intensity change is -0.03; its intensity drops beyond the null in 2 of the 4 tumours, drops relative
to the other paired tumours. Intensity uses all of the MP's genes in ``var_names``, not the panel,
and all of a sample's malignant cells.

**Which tumours are tested (κ): all of them.** Every paired tumour with a defined Δlog κ, null
bound and self bound (below) is tested; there is no activity gate. P1-GO v1.2 tested a paired tumour
only where the MP's κ under one of the two protocols was above ``max(kappa_active, q95 of the
matched random sets' κ)``, and let a replicate pair into the null on the same condition. That selects
tumours on the noisy κ that is then tested: a tumour whose reference κ is high by chance passes the
gate and tends to drop. With every piece's κ drawn like its random sets, up to 0.31 of the tested
tumours were losses on the random-set path (0.05-0.31 across the exchangeable synthetic laws),
0.13-0.17 of the tested (set, tumour) pairs without the program on the toy null, and an unplanted
random gene set met criterion (c) on the toy (``scripts/calibration/README.md``, v1.2 reference). Without the gate, which tumours count no
longer depends on the noise of the statistic tested. The activity columns are kept for information
(``kappa_active_reference``, ``kappa_active_test``, ``kappa_active``: κ above ``max(kappa_active,
q95)``, ``kappa_active`` 0.05, ``KAPPA_ACTIVE``). The price: tumours without the program count in
the denominator of criterion (c), so an MP present in fewer than half of the paired tumours cannot
meet it; such subset programs can reach nuisance support only through (a) or (b).

**Bounds (κ).** A tested tumour is a loss iff ``Δlog κ < min(null_lo, self_lo, -kappa_min_effect)``
and a gain iff ``Δlog κ > max(null_hi, self_hi, kappa_min_effect)`` (gains are informational). Each
of the null bound, the self bound and the margin is a one-sided bound under its own assumption
(exchangeable replicate pairs or exchangeable random sets; cell-sampling noise; a change too small to
matter), and a loss must beat all of them, so the 5% level holds when any one of the null's and the
self bound's assumptions does.

``replicates`` (null bound): same-protocol pairs of samples of one tumour (pieces from different
regions) with the statistic finite in both pieces, used when such pairs come from >= ``min_pairs``
(``MIN_PAIRS``, 3) distinct tumours; every such pair counts. A bound is a one-sided 95% prediction
bound for a new change, ``t(0.95, df) * s`` (:func:`replicate_bound`): ``s = sqrt(mean d_i^2)`` over
the signed replicate changes ``d_i`` (the order within a pair is arbitrary, so their mean is 0), and
``df`` is the number of distinct tumours behind the pairs, since pairs of a tumour with three or more
pieces share pieces and are not independent (counting tumours is conservative when each tumour has
few pieces; with very unbalanced piece counts the effective df can be slightly below the tumour
count, a negligible effect on ``t(0.95, df)``). Calibration: with no protocol effect and paired and
replicate changes exchangeable and normal, a tumour is a loss by chance in 5% of cases whatever the
number of tumours (Monte Carlo with 3, 5 and 10 replicate tumours: 0.050 per tumour; 0.034 for >= 2
of 4 paired tumours with 3); the empirical 95th percentile of 3 absolute changes, used before v1.1,
gave 0.13 per tumour and 0.11 for >= 2 of 4. For κ there are two such bounds on the same pairs and
df, one on the log changes ``d_i = log max(κ_j, floor) - log max(κ_i, floor)`` and one on the linear
changes ``d_i = κ_j - κ_i``, carried to the MP's κ level like the random-set bound below, and
``null_lo`` is the lower of the two (``null_hi`` the higher). The log bound alone is too narrow for
a tumour whose MP sits near the floor when the replicate pieces carry it well above the floor (an
absolute noise near 0 is a large log change; synthetic ``replicates_carrier_reps`` case of
``scripts/calibration``); the linear bound alone is too narrow for a tumour carrying an MP whose κ
varies between pieces when the replicate pieces sit at the floor (a proportional change is a large
linear change; 0.13-0.20 of such tumours were losses in the design comparison,
``tests/test_protocol.py`` checks <= 0.05 for the lower of the two). This null contains regional
heterogeneity, but most public protocol twins have no replicate pieces, and the t bounds assume
near-normal changes: for a program absent from a share of the pieces (patchy), the changes are a
mixture of small ones and drops to the floor, all paired tumours share one bound estimated from a
few such changes, and losses cluster within the MP. Per tumour the rate stays near 5%, but with 3
replicate tumours the chance that such an MP meets criterion (c) on 4 paired tumours was 0.024-0.042
(absent from 10-30% of pieces; synthetic, ``scripts/calibration`` ``C2-info``), against at most
about 0.02 for normal laws.

``random_sets`` (null bound), otherwise: per paired tumour, ``n_null_sets`` random gene sets drawn
from the tumour's panel, matched gene by gene to the MP's panel genes on mean-expression and
detection bins in the tumour's reference-protocol cells (:func:`mpaudit.coexpr.matched_random_sets`;
from the (tumour, MP)'s own stream, :func:`_sets_rng`, so an MP's sets do not depend on the other
MPs of the call), and evaluated like the MP in the same samples (κ on the same cells, intensity with the same
protocol-specific control machinery). Per tumour, the bounds are the ``k``-th smallest and the
``k``-th largest of the n random sets' changes with a finite value, ``k = floor(0.05 (n + 1))`` (the
5th of 100; the smallest with the minimum ``MIN_RANDOM_SETS`` = 20, :func:`_order_bounds`): for a
change exchangeable with the sets' changes, ``P(change < k-th smallest) = k / (n + 1) <= 0.05``
whatever their distribution (v1.2 used the interpolated 5th percentile, which is not such a
bound). With fewer than 20 finite changes, or a k-th smallest not below the k-th largest (every set
changed alike, e.g. all at κ 0), the tumour is not tested. Intensity is compared with these bounds
directly. κ is not: random sets have κ near 0 while an MP's κ is usually far higher, and an additive
change in background κ is a large change on the log scale near 0 and a negligible one at the MP's
level. Compared on the log scale, the random sets' Δlog κ called unchanged MPs losses when background
κ rose and hid twofold losses when the sets sat near the floor. The order statistics ``L_lo``,
``L_hi`` of the random sets' *linear* changes ``κ_test - κ_ref`` are therefore carried to the MP's
level: ``null_lo = log(max(κ_ref + L_lo, floor)) - log κ_ref`` and likewise ``null_hi``, with
``κ_ref`` the MP's floored reference κ. This null holds what all genes share (ambient co-expression,
estimation noise of near-zero correlations) as an additive change in κ; a shared proportional
change is not in it, which is why depth, whose effect on the estimate scales with the set's κ
level, is equalised before κ rather than left to this null. It holds neither the MP's own sampling
noise (for a coherent
program κ depends on how many sampled cells are in the state; the self bound covers it) nor
regional heterogeneity: two pieces that differ in a program's activity look like a protocol effect,
and only the margin protects against that. On the toy run the random-set bounds lie within +-0.23 for
stress, cycle, hypoxia, EMT, abundance and long genes, and reach 1.04 for CNA_7q (in T01, which
lacks it), 0.46 for IFN and 0.66 for contamination. Prevalence has no random-set null (its
threshold would need every random set scored on the whole pool). The same random sets give the
informational activity threshold, so they are drawn for every MP, also in the tumours that only
have replicate pieces (matched on the tumour's reference-protocol cells, or on all its cells when
it has none).

**Self bound (κ).** Per sample, ``kappa_se`` is the bootstrap SE of κ: the SD (ddof 1) over
``KAPPA_BOOT`` (50) resamples, with replacement, of the sample's κ cells, drawn from a stream seeded
by the seed, the sample and the MP (:func:`_boot_rng`; so an MP's SE does not depend on the other
MPs of the call). The SE of a protocol's mean κ is ``sqrt(sum SE_s^2) / n`` over its n samples with a
finite κ (independent samples), ``sd = sqrt(SE_ref^2 + SE_test^2)`` is that of the change, and
``self_lo = log(max(κ_ref - 1.645 sd, floor)) - log κ_ref`` (``KAPPA_SELF_Z`` = 1.645; ``self_hi``
likewise with ``+``). A tumour where any of those samples has an undefined SE is not tested. On the
toy the bootstrap SE was 0.90-1.11 of the true SD of κ for absent gene sets at 60 and 200 cells
(design comparison); 50 resamples leave about 10% noise in the SE, and the bootstrap ignores
dependence between cells within a piece (doublets, ambient RNA). Only samples of paired tumours get
a ``kappa_se`` (NaN elsewhere: replicate tumours enter the null with their κ alone). On the toy run
the self bound lies between about -0.08 and -0.20 where a program is planted (contamination, κ
about 0.3: -0.39 to -0.54) and between -0.68 and -1.39 where it is absent (CNA_7q in T01 and T03,
IFN in T03).

**Minimum effect (κ).** ``kappa_min_effect = log 1.5`` (``KAPPA_MIN_EFFECT``, the proposal's
equivalence margin for Δlog κ). A null can be narrower than any change that matters (many cells, a
tight replicate set, or the random-set null above), and the margin keeps such changes, shared
proportional shifts below 1.5-fold and regional heterogeneity below 1.5-fold out of criterion (c).
On the toy run, among the 32 κ tests of the eight planted programs other than stress (whose
co-expression the toy protocol does not change by design; EMT's prevalence is planted higher under
the fixed protocol), the null bound alone calls 5 losses and 7 gains, the null and the self bound 3
losses and 3 gains, the null and the margin 1 loss (CNA_7q in T03, a tumour without the CNA: Δlog κ
-0.611 against a null bound of -0.550; its self bound is -0.675), and all three none; stress is a
loss in 4 of 4 paired tumours. Intensity and prevalence are informational and have no margin and no
self bound.

**Calibration (toy and synthetic only, not real data).** ``scripts/calibration/evaluate_c.py``; the
v1.3 reference results are in its README: on exchangeable synthetic laws 0.018-0.054 of the tested
tumours are losses (v1.2: up to 0.31), the chance that an unaffected MP meets (c) is at most 0.025 on
4 paired tumours (0.003 on 10), and on the toy null 0.019-0.063 of the tested pairs without the
program are losses, none with it, and no gene set meets (c); with the test pieces at a half or a
quarter of the reference depth 0.022-0.050 of the pairs without the program and at most 0.017 of the
carriers are losses, and with 50-60-cell pieces at most 0.043 and 0.020. Residual gaps: regional
heterogeneity of a present program on the random-set path (only the margin covers it), a program
absent from part of the pieces on the replicate path (P(c) up to 0.047 with 3 replicate tumours,
above), background co-variation that differs between pieces beyond cell sampling and beyond what
matched random sets share, and depth differences between the paired and the replicate tumours,
which equal depth within a tumour does not remove (on the toy, paired pieces at a quarter of the
replicate pieces' depth stayed within the target above).

**Data.** Oliveira et al. (Nat Genet 2025) is NOT a fresh-vs-fixed twin: its Chromium Flex, Visium HD
and Xenium data come from the same FFPE blocks and no fresh-dissociated malignant cells were found,
so it offers a fixed-dissociated versus in situ comparison, not a fresh-vs-fixed malignant contrast.

Effects on the toy cohort (:mod:`mpaudit.sim`) are true by construction; they show that the code
responds as designed, not how often real MPs change between protocols. "Toy run" figures above come
from ``simulate_cohort()`` with its defaults (seed 0; calibrated magnitudes; 4 paired tumours
T01-T04 and 2 replicate tumours, every sample 200 malignant cells), the planted ``truth["sets"]`` and
``protocol_shift_kappa(adata, keys, sets, "fresh", "fixed")`` with its defaults (``n_null_sets=100``,
``seed=0``, ``min_pairs=3``; every MP gets the random-set null, as only 2 tumours have replicate
pieces), re-run on 2026-10-18 (with equal depth and the random sets drawn per tumour and MP); the
variants without the margin or the self bound change only those settings. The small-cohort figures
use ``SMALL`` of ``tests/test_protocol.py``.
"""
from __future__ import annotations

import itertools
import zlib
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .coexpr import kappa, matched_random_sets, resample_kappas, set_kappas
from .data import Keys, get_counts, log2_cpm10, malignant_mask

KAPPA_FLOOR = 0.01
KAPPA_MIN_EFFECT = float(np.log(1.5))  # smallest |Δlog κ| counted as a loss or gain (the proposal's margin)
KAPPA_ACTIVE = 0.05  # informational: floor of the activity threshold max(KAPPA_ACTIVE, random-set q95)
MIN_PAIRS = 3  # distinct tumours with same-protocol replicate pairs needed for the replicate null
MIN_RANDOM_SETS = 20  # finite random-set changes needed for a tumour's random-set bounds (k = 1 at 20)
NULL_ALPHA = 0.05  # one-sided level of each null bound (replicate t bound, random-set order statistic)
KAPPA_BOOT = 50  # bootstrap resamples of a sample's κ cells for the SE of the MP's own κ (self bound)
KAPPA_SELF_Z = 1.645  # self bound: a κ loss must also exceed KAPPA_SELF_Z x the SE of the MP's own change
MAX_KAPPA_CELLS = 2000  # κ cells per sample: the tumour's smallest sample, at most this many
MIN_CELLS = 20  # malignant cells a sample needs to enter the comparison
MIN_DETECT = 0.05  # panel: genes detected in >= this share of the tumour's malignant cells under each protocol
MIN_IOTA = 0.5  # the MP's κ in a tumour is undefined when less than this share of its genes is in the panel
MIN_GENES = 5  # ... or when fewer than this many of its genes are
N_BINS = 30  # expression bins for control genes, as in mpaudit.scoring
METRICS = ("kappa", "intensity", "prevalence")


def _named_rng(seed: int, *names: str) -> np.random.Generator:
    """Stream fixed by ``seed`` and ``names``: ``default_rng([seed, crc32(names joined by '\\x1f')])``, so a
    draw does not depend on the other samples or MPs of the call or on their order."""
    return np.random.default_rng([int(seed), zlib.crc32("\x1f".join(map(str, names)).encode())])


def _cells_rng(seed: int, sample: str) -> np.random.Generator:
    """Stream of the κ-cell subsample of one sample."""
    return _named_rng(seed, "kappa-cells", sample)


def _depth_rng(seed: int, sample: str) -> np.random.Generator:
    """Stream of the binomial thinning of one sample's κ cells to the tumour's lowest median depth."""
    return _named_rng(seed, "kappa-depth", sample)


def _boot_rng(seed: int, sample: str, mp: str) -> np.random.Generator:
    """Stream of the bootstrap resamples of one (sample, MP): ``default_rng([seed, crc32(sample + '\\x1f' + mp)])``."""
    return _named_rng(seed, sample, mp)


def _sets_rng(seed: int, tumour: str, mp: str) -> np.random.Generator:
    """Stream of the random gene sets of one (tumour, MP), so an MP's random sets do not depend on the other MPs
    of the call, on the tumours or on their order."""
    return _named_rng(seed, "random-sets", tumour, mp)


def _thin(Y: sp.csr_matrix, p: float, rng: np.random.Generator) -> sp.csr_matrix:
    """``Y`` with every count ``y`` replaced by a ``Binomial(y, p)`` draw (counts must be whole numbers >= 0)."""
    Y = Y.copy()
    Y.sum_duplicates()  # canonical (sorted) entries, so the draw does not depend on the input's storage order
    n = np.rint(Y.data)
    if np.any(np.abs(Y.data - n) > 1e-6) or np.any(n < 0):
        raise ValueError("the κ comparison thins counts to equal depth (binomial thinning) and needs counts that are "
                         "whole numbers >= 0 (UMI counts); give the raw counts (counts_layer)")
    Y.data = rng.binomial(n.astype(np.int64), p).astype(np.float64)
    Y.eliminate_zeros()
    return Y


@dataclass
class _Sample:
    name: str
    tumour: str
    protocol: str
    rows: np.ndarray  # the sample's analysed malignant cells (rows of _Cohort.counts), in barcode order
    kappa_rows: np.ndarray  # the tumour's κ-cell number of them (_prepare), for kappa, in barcode order
    median_umi: float = np.nan  # median total count of the κ cells, before thinning
    thin: float = 1.0  # binomial thinning of the κ cells' counts: the tumour's lowest median_umi / median_umi


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
    kappa_data: dict  # sample -> (CSC counts of its κ cells, thinned by its ``thin``; their depth after thinning)

    def _kappa_cells(self, s: _Sample) -> tuple[sp.csc_matrix, np.ndarray]:
        """Counts of the sample's κ cells (CSC), binomially thinned by ``s.thin`` (:func:`_prepare`), and their
        depth (total counts over all genes, after thinning)."""
        return self.kappa_data[s.name]

    def set_kappa(self, s: _Sample, idx: np.ndarray) -> float:
        """``kappa`` of genes ``idx`` in the sample's kappa cells (depth = total counts over all genes)."""
        Y, depth = self._kappa_cells(s)
        return kappa(Y[:, idx].toarray(), depth)

    def set_kappas(self, s: _Sample, sets: list[np.ndarray]) -> np.ndarray:
        """:meth:`set_kappa` of each gene set in ``sets``, from one :func:`mpaudit.coexpr.set_kappas` call."""
        Y, depth = self._kappa_cells(s)
        return set_kappas(Y, depth, sets)

    def set_kappa_se(self, s: _Sample, idx: np.ndarray, seed: int, mp: str) -> float:
        """Bootstrap SE of :meth:`set_kappa`: the SD (ddof 1) of ``kappa`` over ``KAPPA_BOOT`` resamples,
        with replacement, of the sample's κ cells, drawn from :func:`_boot_rng` (``seed``, sample, ``mp``);
        resamples with an undefined ``kappa`` are left out, and the SE is NaN with fewer than 2 left."""
        Y, depth = self._kappa_cells(s)
        rng = _boot_rng(seed, s.name, mp)
        n = Y.shape[0]
        boot = np.array([rng.integers(0, n, n) for _ in range(KAPPA_BOOT)])
        vals = resample_kappas(Y[:, idx].toarray(), depth, boot)  # kappa() of each resample, computed together
        vals = vals[np.isfinite(vals)]
        return float(np.std(vals, ddof=1)) if vals.size >= 2 else np.nan

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


def _expression_bins(mean: np.ndarray) -> np.ndarray:
    """``N_BINS`` equal-size bins of genes ranked by ``mean`` (as in :func:`mpaudit.scoring.control_genes`)."""
    bins = np.empty(mean.size, dtype=int)
    bins[np.argsort(mean, kind="stable")] = np.arange(mean.size) * N_BINS // max(mean.size, 1)
    return bins


def _detected(X: sp.csr_matrix) -> np.ndarray:
    return np.asarray((X > 0).sum(axis=0)).ravel() / max(X.shape[0], 1)


def _prepare(adata, keys: Keys, reference: str, test: str, *, counts_layer: str | None, min_cells: int,
             max_cells: int, min_detect: float, seed: int) -> _Cohort:
    if keys.protocol is None:
        raise ValueError("keys.protocol is required")
    if reference == test:
        raise ValueError("reference and test protocols must differ")
    obs = adata.obs
    protocol = obs[keys.protocol].astype(str).to_numpy()
    sample = obs[keys.sample].astype(str).to_numpy()
    tumour = obs[keys.tumour_col()].astype(str).to_numpy()
    counts = get_counts(adata, counts_layer)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    cells = malignant_mask(adata, keys) & np.isin(protocol, [reference, test]) & (totals > 0)
    n_per_sample = pd.Series(sample[cells]).value_counts()
    cells &= np.isin(sample, n_per_sample.index[n_per_sample >= min_cells])
    rows = np.flatnonzero(cells)
    counts, totals = counts[rows], totals[rows]
    protocol, sample, tumour = protocol[rows], sample[rows], tumour[rows]
    barcode = adata.obs_names.astype(str).to_numpy()[rows]

    samples = []
    for s in pd.unique(sample):
        r = np.flatnonzero(sample == s)
        if np.unique(protocol[r]).size > 1:
            raise ValueError(f"sample {s!r} has malignant cells under both protocols")
        r = r[np.argsort(barcode[r], kind="stable")]  # barcode order: the draws below ignore the cohort's cell order
        samples.append(_Sample(s, tumour[r[0]], protocol[r[0]], r, r))
    # κ falls with the number of cells it is computed on, so every sample of a tumour is scored on the same
    # number n_t = min(max_cells, the tumour's smallest sample), drawn without replacement from the sample's
    # own stream (_cells_rng: independent of the other samples and their order)
    n_t: dict[str, int] = {}
    for smp in samples:
        n_t[smp.tumour] = min(n_t.get(smp.tumour, max_cells), smp.rows.size)
    for smp in samples:
        if smp.rows.size > n_t[smp.tumour]:
            pick = _cells_rng(seed, smp.name).choice(smp.rows.size, n_t[smp.tumour], replace=False)
            smp.kappa_rows = smp.rows[np.sort(pick)]
    # κ's moment estimate also falls with depth (counts per cell), so the κ cells of every sample of a tumour are
    # binomially thinned to the tumour's lowest median depth (_Cohort._kappa_cells, stream _depth_rng)
    low: dict[str, float] = {}
    for smp in samples:
        smp.median_umi = float(np.median(totals[smp.kappa_rows]))
        low[smp.tumour] = min(low.get(smp.tumour, np.inf), smp.median_umi)
    kappa_counts = {}
    for smp in samples:
        smp.thin = min(1.0, low[smp.tumour] / smp.median_umi)
        Y = counts[smp.kappa_rows]
        kappa_counts[smp.name] = _thin(Y, smp.thin, _depth_rng(seed, smp.name)) if smp.thin < 1 else Y

    # the panel, on the counts κ is computed on (thinned κ cells), so that both protocols select genes alike
    panel, paired = {}, []
    for t in pd.unique(tumour):
        protocols = {s.protocol for s in samples if s.tumour == t}
        ok = np.ones(counts.shape[1], bool)
        for p in protocols:
            ok &= _detected(sp.vstack([kappa_counts[s.name] for s in samples
                                       if s.tumour == t and s.protocol == p]).tocsr()) >= min_detect
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
    kappa_data = {k: (Y.tocsc(), np.asarray(Y.sum(axis=1)).ravel()) for k, Y in kappa_counts.items()}
    return _Cohort(reference, test, adata.var_names.astype(str).to_numpy(), counts, L, protocol, samples,
                   panel, paired, pool, bins, kappa_data)


LEVEL_COLUMNS = ["sample", "tumour", "protocol", "mp", "n_cells", "n_kappa_cells", "median_umi", "kappa_thin",
                 "n_genes", "iota", "kappa", "kappa_se", "intensity", "prevalence"]


def _levels(cohort: _Cohort, gene_sets: dict[str, list[str]], *, min_iota: float, min_genes: int,
            high_quantile: float, seed: int) -> pd.DataFrame:
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
    paired = set(cohort.paired)
    for s in cohort.samples:
        panel = cohort.panel[s.tumour]
        for n in gene_sets:
            target = idx[n][panel[idx[n]]]
            iota = target.size / len(members[n]) if members[n] else np.nan
            v, t = score[n][s.rows], thr.get((n, s.protocol), np.nan)
            ok = iota >= min_iota and target.size >= min_genes
            out.append({"sample": s.name, "tumour": s.tumour, "protocol": s.protocol, "mp": n,
                        "n_cells": int(s.rows.size), "n_kappa_cells": int(s.kappa_rows.size),
                        "median_umi": s.median_umi, "kappa_thin": s.thin,
                        "n_genes": int(target.size), "iota": float(iota),
                        "kappa": cohort.set_kappa(s, target) if ok else np.nan,
                        # only paired tumours are tested, so only their samples need the self bound's SE
                        "kappa_se": cohort.set_kappa_se(s, target, seed, n) if ok and s.tumour in paired else np.nan,
                        "intensity": float(np.mean(v)),
                        "prevalence": float(np.mean(v > t)) if np.isfinite(t) else np.nan})
    return pd.DataFrame(out, columns=LEVEL_COLUMNS)


def _null_tumours(cohort: _Cohort) -> list[str]:
    """Paired tumours, then the other tumours with two or more samples under one protocol (replicate pieces)."""
    n = pd.Series([(s.tumour, s.protocol) for s in cohort.samples]).value_counts()
    replicate = {t for (t, _), k in n.items() if k > 1}
    return cohort.paired + [t for t in dict.fromkeys(s.tumour for s in cohort.samples)
                            if t in replicate and t not in cohort.paired]


def _random_set_levels(cohort: _Cohort, gene_sets: dict[str, list[str]], mps: list[str], n_sets: int, *,
                       min_iota: float, min_genes: int, seed: int) -> pd.DataFrame:
    """κ and intensity of matched random gene sets, per sample of each tumour in :func:`_null_tumours`.

    The sets are matched on the tumour's reference-protocol cells (all its cells when it has none) and
    drawn from the (tumour, MP)'s own stream (:func:`_sets_rng`). Each sample's κ of all the tumour's
    random sets comes from one :meth:`_Cohort.set_kappas` call (the same values as set by set).
    """
    pos = {g: i for i, g in enumerate(cohort.genes)}
    out = []
    for t in _null_tumours(cohort):
        ss = [s for s in cohort.samples if s.tumour == t]
        ref = [s for s in ss if s.protocol == cohort.reference] or ss
        Y = cohort.counts[np.concatenate([s.rows for s in ref])]
        mean_expr = np.asarray(Y.sum(axis=0)).ravel() / max(Y.sum(), 1e-12)
        detect = _detected(Y)
        L_s = {s.name: cohort.L[s.rows] for s in ss}
        panel = cohort.panel[t]
        drawn = []  # (MP, set number, genes), in the order of the draws
        for n in mps:
            genes = list(dict.fromkeys(map(str, gene_sets[n])))
            target = np.array([pos[g] for g in genes if g in pos and panel[pos[g]]], dtype=int)
            if not genes or target.size / len(genes) < min_iota or target.size < min_genes:
                continue
            drawn += [(n, j, r) for j, r in enumerate(matched_random_sets(target, mean_expr, detect, n_sets,
                                                                          _sets_rng(seed, t, n), candidates=panel))]
        k = {s.name: cohort.set_kappas(s, [r for _, _, r in drawn]) for s in ss}
        protocols = sorted({s.protocol for s in ss})
        for a in range(0, len(drawn), 256):  # intensities of a block of sets at once: L_s @ (genes x sets weights)
            block = drawn[a:a + 256]
            W = {p: np.column_stack([cohort.score_weights(p, r) for _, _, r in block]) for p in protocols}
            inten = {s.name: np.ascontiguousarray((L_s[s.name] @ W[s.protocol]).T).mean(axis=1) for s in ss}
            for i, (n, j, r) in enumerate(block, start=a):
                for s in ss:
                    out.append({"mp": n, "tumour": t, "sample": s.name, "protocol": s.protocol, "set": j,
                                "kappa": float(k[s.name][i]), "intensity": float(inten[s.name][i - a])})
    return pd.DataFrame(out, columns=["mp", "tumour", "sample", "protocol", "set", "kappa", "intensity"])


def per_sample_kappa(
    adata,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    reference: str,
    test: str,
    *,
    counts_layer: str | None = None,
    min_cells: int = MIN_CELLS,
    max_cells: int = MAX_KAPPA_CELLS,
    min_detect: float = MIN_DETECT,
    min_iota: float = MIN_IOTA,
    min_genes: int = MIN_GENES,
    high_quantile: float = 0.9,
    seed: int = 0,
) -> pd.DataFrame:
    """Per malignant sample of the two protocols and MP: κ on the tumour's panel, its bootstrap SE, intensity
    and prevalence.

    Samples need >= ``min_cells`` malignant cells; see the module docstring for the statistics.

    Columns: ``sample, tumour, protocol, mp, n_cells`` (malignant cells analysed: intensity and
    prevalence use all of them), ``n_kappa_cells`` (the cells κ is computed on: the same number for
    every sample of a tumour, ``min(max_cells, the tumour's smallest sample)``), ``median_umi`` (median
    total count of those cells) and ``kappa_thin`` (the share their counts are binomially thinned to:
    the tumour's lowest ``median_umi`` / the sample's, 1 for its shallowest sample), ``n_genes`` (panel
    genes used for κ), ``iota``, ``kappa`` (NaN when ``iota < min_iota`` or fewer than ``min_genes``
    panel genes), ``kappa_se`` (bootstrap SE of ``kappa`` over ``KAPPA_BOOT`` resamples of the κ cells;
    computed for samples of paired tumours only, the only ones tested, and NaN elsewhere, where
    ``kappa`` is, or with fewer than 2 finite resamples), ``intensity``, ``prevalence``. Each sample
    must have one protocol, and the counts must be whole numbers (UMI counts) for the thinning.
    """
    cohort = _prepare(adata, keys, reference, test, counts_layer=counts_layer, min_cells=min_cells,
                      max_cells=max_cells, min_detect=min_detect, seed=seed)
    return _levels(cohort, gene_sets, min_iota=min_iota, min_genes=min_genes, high_quantile=high_quantile, seed=seed)


def _transform(metric: str, x: pd.Series | np.ndarray, kappa_floor: float) -> np.ndarray:
    """κ is compared on the log scale after flooring; intensity and prevalence as they are."""
    x = np.asarray(x, dtype=float)
    return np.log(np.maximum(x, kappa_floor)) if metric == "kappa" else x


def _delta_name(metric: str) -> str:
    return "delta_log_kappa" if metric == "kappa" else f"delta_{metric}"


def _replicate_deltas(levels: pd.DataFrame, metrics: list[str], kappa_floor: float) -> pd.DataFrame:
    """Changes between every two samples of one tumour under the same protocol, per MP (one row per pair).

    κ's change is on the floored log scale (``kappa``) and, as ``kappa_lin``, on the linear scale
    (``κ_j - κ_i``); NaN unless κ is finite in both samples.
    """
    keys = ["mp", "tumour", "protocol"]
    multi = levels[levels.groupby(keys, sort=False)["mp"].transform("size") > 1]
    lin = ["kappa_lin"] if "kappa" in metrics else []
    out = []
    for (mp, tumour, _), g in multi.groupby(keys, sort=False):
        v = {m: _transform(m, g[m], kappa_floor) for m in metrics}
        k = g["kappa"].to_numpy(float) if lin else None
        for i, j in itertools.combinations(range(len(g)), 2):
            rec = {m: v[m][j] - v[m][i] for m in metrics}
            if lin:
                rec["kappa_lin"] = k[j] - k[i]
            out.append({"mp": mp, "tumour": tumour, **rec})
    return pd.DataFrame(out, columns=["mp", "tumour", *metrics, *lin])


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


def replicate_bound(deltas: np.ndarray, n_tumours: int) -> float:
    """One-sided 95% prediction bound for a new change from signed replicate changes: ``t_{0.95, n_tumours} * s``.

    ``s = sqrt(mean(deltas**2))``: the order within a replicate pair is arbitrary, so the changes
    have mean 0 and ``s`` estimates their standard deviation without subtracting a mean. The degrees
    of freedom are the number of distinct tumours behind the pairs (pairs of one tumour with three or
    more pieces share pieces and are not independent). NaN without finite changes or tumours.
    """
    d = np.asarray(deltas, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0 or n_tumours < 1:
        return np.nan
    return float(stats.t.ppf(1 - NULL_ALPHA, n_tumours) * np.sqrt(np.mean(d**2)))


def _null_types(levels: pd.DataFrame, reps: pd.DataFrame, metrics: list[str], min_pairs: int,
                random_metrics: list[str]) -> pd.DataFrame:
    """Per MP and statistic ``m``: ``null_type_m``, ``n_null_pairs_m``, ``n_null_tumours_m`` and ``bound_m``.

    Only replicate pairs whose change in ``m`` is finite count, so an MP whose κ is NaN in its
    replicate tumours falls back to the random-set null for κ while intensity can keep its replicates.
    ``bound_m`` is the replicate prediction bound (NaN unless ``null_type_m == "replicates"``); for κ it is
    the bound on the log scale, and ``bound_kappa_lin`` the one on the linear scale (same pairs, same
    degrees of freedom).
    """
    mps = pd.Index(pd.unique(levels["mp"]), name="mp")
    out = pd.DataFrame(index=mps)
    need = max(int(min_pairs), 1)
    for m in metrics:
        r = reps.loc[np.isfinite(reps[m].to_numpy(float))]
        g = r.groupby("mp", sort=False)
        n_pairs = g.size().reindex(mps, fill_value=0).astype(int)
        n_tumours = g["tumour"].nunique().reindex(mps, fill_value=0).astype(int)
        enough = (n_tumours >= need).to_numpy()
        out[f"null_type_{m}"] = np.where(enough, "replicates", "random_sets" if m in random_metrics else "none")
        out[f"n_null_pairs_{m}"], out[f"n_null_tumours_{m}"] = n_pairs, n_tumours
        for col, name in ((m, f"bound_{m}"), *((("kappa_lin", "bound_kappa_lin"),) if m == "kappa" else ())):
            bound = {mp: replicate_bound(x.to_numpy(), int(n_tumours[mp])) for mp, x in g[col]}
            out[name] = np.where(enough, pd.Series(bound, dtype=float).reindex(mps).to_numpy(), np.nan)
    return out


def _order_bounds(x: pd.Series, by: list[str], alpha: float = NULL_ALPHA) -> tuple[pd.Series, pd.Series]:
    """Per group: the k-th smallest and the k-th largest finite value of ``x``, ``k = floor(alpha (n + 1))``.

    For a new value exchangeable with the group's n values, ``P(new < k-th smallest) = k / (n + 1) <=
    alpha`` (a distribution-free one-sided prediction bound). Groups with ``k < 1`` are left out.
    """
    s = x[np.isfinite(x.to_numpy(float))]
    g = s.groupby(level=by, sort=False)
    n = g.transform("size").to_numpy()
    k = np.floor(alpha * (n + 1) + 1e-9).astype(int)
    r = g.rank(method="first").to_numpy()
    lo = s[(r == k) & (k >= 1)].groupby(level=by, sort=False).first()
    hi = s[(n - r + 1 == k) & (k >= 1)].groupby(level=by, sort=False).first()
    return lo, hi


def _random_set_bounds(nd: pd.DataFrame, metric: str) -> tuple[pd.Series, pd.Series]:
    """Per (MP, tumour): the ``floor(0.05 (n + 1))``-th smallest and largest of the n random sets' changes
    (κ: linear ``κ_test - κ_ref``; :func:`_order_bounds`).

    NaN with fewer than ``MIN_RANDOM_SETS`` finite changes, or when the k-th smallest is not below the
    k-th largest (e.g. every random set at κ 0 in both protocols).
    """
    x = nd["kappa_test"] - nd["kappa_reference"] if metric == "kappa" else nd[metric]
    enough = x.groupby(level=["mp", "tumour"], sort=False).count() >= MIN_RANDOM_SETS
    lo, hi = _order_bounds(x, ["mp", "tumour"])
    lo, hi = lo.reindex(enough.index).where(enough), hi.reindex(enough.index).where(enough)
    collapsed = lo >= hi
    return lo.mask(collapsed), hi.mask(collapsed)


def _q95(x: pd.Series, by) -> pd.Series:
    """95th percentile of ``x`` per group; NaN with fewer than ``MIN_RANDOM_SETS`` finite values."""
    g = x.groupby(by, sort=False)
    return g.quantile(0.95).where(g.count() >= MIN_RANDOM_SETS)


def _protocol_kappa_q95(null_levels: pd.DataFrame | None, protocol: str) -> pd.Series:
    """Per (MP, tumour): q95 over the random sets of their κ under ``protocol`` (mean over its samples).

    Empty without random sets (no ``null_levels``, an empty frame, e.g. when no MP passes the panel
    in any tumour, or no ``kappa`` column).
    """
    if null_levels is None or null_levels.empty or "kappa" not in null_levels:
        return pd.Series(dtype=float)
    sel = null_levels[null_levels["protocol"] == protocol]
    per_set = sel["kappa"].astype(float).groupby([sel["mp"], sel["tumour"], sel["set"]], sort=False).mean()
    return _q95(per_set, ["mp", "tumour"])


def _protocol_mean_se(levels: pd.DataFrame, index: pd.MultiIndex, protocol: str) -> np.ndarray:
    """Per (MP, tumour) of ``index``: SE of the MP's mean κ under ``protocol``, ``sqrt(sum SE_s^2) / n`` over
    the n samples with a finite κ (those the mean is taken over; samples are independent); NaN when any
    of them has an undefined ``kappa_se``."""
    lv = levels[(levels["protocol"] == protocol) & np.isfinite(levels["kappa"].to_numpy(float))]
    v = lv["kappa_se"].astype(float) ** 2
    g = v.groupby([lv["mp"], lv["tumour"]], sort=False)
    se = np.sqrt(g.sum().where(g.count() == g.size())) / g.size()
    return se.reindex(index).to_numpy(float)


def _shifts(levels: pd.DataFrame, reference: str, test: str, null_levels: pd.DataFrame | None, min_pairs: int,
            kappa_floor: float, kappa_min_effect: float, kappa_active: float,
            self_bound: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(per_tumour, null types)``; see :func:`tumour_shifts`."""
    metrics = [m for m in METRICS if m in levels]
    if "kappa" in metrics and self_bound and "kappa_se" not in levels:
        raise ValueError("levels has kappa but no kappa_se column: the registered κ test needs the bootstrap SE "
                         "of each sample's κ (per_sample_kappa or protocol_shift_kappa give it); "
                         "self_bound=False drops the self bound, which is not the registered rule")
    reps = _replicate_deltas(levels, metrics, kappa_floor)
    random_metrics = [m for m in metrics if null_levels is not None and m in null_levels]
    types = _null_types(levels, reps, metrics, min_pairs, random_metrics)
    d = _paired_deltas(levels, ["mp", "tumour"], metrics, reference, test, kappa_floor)
    row_mp = d.index.get_level_values("mp")
    kind = {m: types[f"null_type_{m}"].reindex(row_mp).to_numpy() for m in metrics}
    base = np.maximum(d["kappa_reference"].to_numpy(float), kappa_floor) if "kappa" in metrics else None

    def carry(delta_lin):  # a linear change of κ carried to the MP's floored reference κ, as a change in log κ
        with np.errstate(invalid="ignore"):
            return np.log(np.maximum(base + delta_lin, kappa_floor)) - np.log(base)

    lo = pd.DataFrame(np.nan, index=d.index, columns=metrics)
    hi = lo.copy()
    for m in metrics:
        rep = kind[m] == "replicates"
        b = np.where(rep, types[f"bound_{m}"].reindex(row_mp).to_numpy(float), np.nan)
        lo[m], hi[m] = -b, b
        if m == "kappa":  # the lower of the linear bound carried to the MP's κ and the log-scale bound
            b_lin = np.where(rep, types["bound_kappa_lin"].reindex(row_mp).to_numpy(float), np.nan)
            lo[m], hi[m] = np.minimum(carry(-b_lin), -b), np.maximum(carry(b_lin), b)
    if random_metrics and not null_levels.empty:
        nd = _paired_deltas(null_levels, ["mp", "tumour", "set"], random_metrics, reference, test, kappa_floor)
        for m in random_metrics:
            q_lo, q_hi = (q.reindex(d.index).to_numpy(float) for q in _random_set_bounds(nd, m))
            if m == "kappa":  # carry the linear change of background κ to the MP's own κ level
                q_lo, q_hi = carry(q_lo), carry(q_hi)
            use = kind[m] == "random_sets"
            lo.loc[use, m], hi.loc[use, m] = q_lo[use], q_hi[use]

    out = d.reset_index().rename(columns={m: _delta_name(m) for m in metrics})
    out.insert(2, "null_type", kind[metrics[0]] if metrics else "none")
    if "iota" in levels:
        iota = levels.groupby(["mp", "tumour"])["iota"].first()
        out.insert(3, "iota", iota.reindex(d.index).to_numpy())
    for m in metrics:
        delta, lo_m, hi_m = d[m].to_numpy(), lo[m].to_numpy(), hi[m].to_numpy()
        tested = np.isfinite(delta) & np.isfinite(lo_m) & np.isfinite(hi_m)
        cut_lo, cut_hi = lo_m, hi_m
        if m == "kappa":
            # the self bound: the MP's own sampling noise, SE of the change of the protocol means
            se = {side: (_protocol_mean_se(levels, d.index, p) if self_bound else np.full(len(d), np.nan))
                  for side, p in (("reference", reference), ("test", test))}
            sd = np.sqrt(se["reference"] ** 2 + se["test"] ** 2)
            self_lo, self_hi = carry(-KAPPA_SELF_Z * sd), carry(KAPPA_SELF_Z * sd)
            out["kappa_se_reference"], out["kappa_se_test"] = se["reference"], se["test"]
            out["self_lo_kappa"], out["self_hi_kappa"] = self_lo, self_hi
            if self_bound:  # an undefined SE leaves the tumour untested
                tested &= np.isfinite(self_lo) & np.isfinite(self_hi)
                cut_lo, cut_hi = np.fmin(cut_lo, self_lo), np.fmax(cut_hi, self_hi)
            # and the minimum effect
            cut_lo, cut_hi = np.minimum(cut_lo, -kappa_min_effect), np.maximum(cut_hi, kappa_min_effect)
            # activity, informational only (no gate): above max(kappa_active, q95 of the random sets' κ)
            either = np.zeros(len(d), bool)
            for side, protocol in (("reference", reference), ("test", test)):
                k = d[f"kappa_{side}"].to_numpy(float)
                q95 = _protocol_kappa_q95(null_levels, protocol).reindex(d.index).to_numpy(float)
                threshold = np.fmax(kappa_active, q95)
                with np.errstate(invalid="ignore"):
                    active = k > threshold
                either |= active
                out[f"kappa_{side}_q95_random"], out[f"kappa_active_threshold_{side}"] = q95, threshold
                out[f"kappa_active_{side}"] = pd.array(np.where(np.isfinite(k), active, None), dtype="boolean")
            out["kappa_active"] = pd.array(np.where(np.isfinite(delta), either, None), dtype="boolean")
        out[f"null_type_{m}"] = kind[m]
        out[f"null_lo_{m}"], out[f"null_hi_{m}"] = lo_m, hi_m
        with np.errstate(invalid="ignore"):
            out[f"loss_{m}"] = pd.array(np.where(tested, delta < cut_lo, None), dtype="boolean")
            out[f"gain_{m}"] = pd.array(np.where(tested, delta > cut_hi, None), dtype="boolean")
    return out, types


def tumour_shifts(
    levels: pd.DataFrame,
    reference: str,
    test: str,
    *,
    null_levels: pd.DataFrame | None = None,
    min_pairs: int = MIN_PAIRS,
    kappa_floor: float = KAPPA_FLOOR,
    kappa_min_effect: float = KAPPA_MIN_EFFECT,
    kappa_active: float = KAPPA_ACTIVE,
    self_bound: bool = True,
) -> pd.DataFrame:
    """Per MP and paired tumour: the ``test - reference`` change of each statistic, and whether it is beyond the null.

    ``levels`` has one row per sample and MP (``mp, tumour, protocol`` and any of ``kappa``,
    ``intensity``, ``prevalence``; with ``kappa`` also ``kappa_se``, and ``sample`` for the
    informational per-piece columns; e.g. from :func:`per_sample_kappa`). ``null_levels`` holds random
    gene sets evaluated in the paired tumours and the tumours with replicate pieces (``mp, tumour,
    sample, protocol, set`` and statistics); for each statistic it is the null of the MPs whose
    same-protocol replicate pairs with a finite change come from fewer than ``min_pairs`` distinct
    tumours. Every paired tumour with a defined change and null bounds is tested (for κ also a defined
    self bound); there is no activity gate (module docstring).

    Columns: ``mp, tumour``, ``null_type`` (of the first statistic present, κ when given), ``iota``
    (when in ``levels``) and per statistic the change (``delta_log_kappa``, ``delta_<metric>``) and the
    protocol means (``<metric>_reference``, ``<metric>_test``). For κ: ``kappa_se_reference``,
    ``kappa_se_test`` (SE of the protocol means: ``sqrt(sum of the samples' kappa_se^2) / n``),
    ``self_lo_kappa``, ``self_hi_kappa`` (the self bounds ``-/+ KAPPA_SELF_Z x sqrt(kappa_se_reference^2
    + kappa_se_test^2)`` carried to the MP's κ level; NaN with ``self_bound=False``), and, for
    information only, per side (``reference``, ``test``) ``kappa_<side>_q95_random`` (NaN without random
    sets), ``kappa_active_threshold_<side>`` (``max(kappa_active, q95)``) and ``kappa_active_<side>``
    (``kappa_<side>`` above it; NA when that κ is NaN), and ``kappa_active`` (either side; NA when the
    change is undefined). Per statistic: ``null_type_<metric>`` (``replicates``, ``random_sets`` or
    ``none``), the null bounds ``null_lo_<metric>``, ``null_hi_<metric>`` (replicates: -/+
    :func:`replicate_bound`, for κ the lower (upper) of the linear bound carried to the MP's κ level and
    the log-scale bound; random sets: the tumour's ``floor(0.05 (n + 1))``-th smallest and largest
    change, for κ the linear changes carried to the MP's κ level) and ``loss_<metric>``,
    ``gain_<metric>`` (NA when not tested). A κ loss needs ``delta_log_kappa < min(null_lo_kappa,
    self_lo_kappa, -kappa_min_effect)``, a gain ``> max(null_hi_kappa, self_hi_kappa,
    kappa_min_effect)``; intensity and prevalence are compared with the null bounds alone.

    ``levels`` with ``kappa`` but without ``kappa_se`` raises ``ValueError`` unless ``self_bound=False``,
    which drops the self bound; that is not the registered rule (for κ levels without cells, e.g.
    synthetic tests). ``kappa_active`` (``KAPPA_ACTIVE``) only sets the informational columns.
    """
    return _shifts(levels, reference, test, null_levels, min_pairs, kappa_floor, kappa_min_effect, kappa_active,
                   self_bound)[0]


def _summarise(shifts: pd.DataFrame, types: pd.DataFrame) -> pd.DataFrame:
    metrics = [m for m in METRICS if f"loss_{m}" in shifts]
    first = metrics[0] if metrics else None
    columns = ["mp", "null_type", "n_null_pairs", "n_null_tumours", "n_paired_tumours", "iota"]
    for m in metrics:
        columns += [f"null_type_{m}", f"mean_{_delta_name(m)}", f"frac_loss_beyond_null_{m}",
                    f"frac_gain_beyond_null_{m}"]
    by_mp = dict(tuple(shifts.groupby("mp", sort=False)))
    out = []
    for mp, row in types.iterrows():
        t = by_mp.get(mp, shifts.iloc[:0])
        rec = {"mp": mp, "null_type": row[f"null_type_{first}"] if first else "none",
               "n_null_pairs": int(row[f"n_null_pairs_{first}"]) if first else 0,
               "n_null_tumours": int(row[f"n_null_tumours_{first}"]) if first else 0,
               "n_paired_tumours": int(t[f"loss_{first}"].notna().sum()) if first else 0,
               "iota": float(t["iota"].median()) if "iota" in t and len(t) else np.nan}
        for m in metrics:
            delta = t[_delta_name(m)].to_numpy(float)
            if m == "kappa":  # over the tested tumours
                delta = delta[t["loss_kappa"].notna().to_numpy(bool)]
            loss, gain = t[f"loss_{m}"].dropna(), t[f"gain_{m}"].dropna()
            rec[f"null_type_{m}"] = row[f"null_type_{m}"]
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
    kappa_min_effect: float = KAPPA_MIN_EFFECT,
    kappa_active: float = KAPPA_ACTIVE,
    self_bound: bool = True,
) -> pd.DataFrame:
    """Per MP: mean change and fraction of tested paired tumours with a loss (gain) beyond the null.

    Arguments as for :func:`tumour_shifts`. Columns: ``mp``; for the first statistic present (κ when
    given): ``null_type`` (``replicates``, ``random_sets``, or ``none`` when neither null is
    available), ``n_null_pairs`` and ``n_null_tumours`` (same-protocol replicate pairs with a finite
    change, and the distinct tumours they come from), ``n_paired_tumours`` (paired tumours tested: change
    and null bounds defined and, for κ, the self bound defined); ``iota`` (median over paired tumours);
    and per statistic ``null_type_<metric>``, ``mean_delta_log_kappa`` (over the tested paired tumours)
    / ``mean_delta_<metric>`` (over the paired tumours with a finite change),
    ``frac_loss_beyond_null_<metric>`` and ``frac_gain_beyond_null_<metric>`` (over the paired tumours
    where that statistic was tested).
    """
    return _summarise(*_shifts(levels, reference, test, null_levels, min_pairs, kappa_floor, kappa_min_effect,
                               kappa_active, self_bound))


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
    kappa_min_effect: float = KAPPA_MIN_EFFECT,
    kappa_active: float = KAPPA_ACTIVE,
    min_cells: int = MIN_CELLS,
    max_cells: int = MAX_KAPPA_CELLS,
    min_detect: float = MIN_DETECT,
    min_iota: float = MIN_IOTA,
    min_genes: int = MIN_GENES,
    high_quantile: float = 0.9,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Axis A4 end to end. Returns ``(per_mp, per_tumour, levels)``.

    ``levels`` is :func:`per_sample_kappa` (same arguments, same values, ``kappa_se`` included). Random
    sets are drawn for every MP (``n_null_sets`` per paired tumour and per tumour with same-protocol
    replicate pieces; at least ``MIN_RANDOM_SETS`` finite changes are needed for a tumour's bounds):
    they are the null of the MPs whose replicate pairs with a finite change come from fewer than
    ``min_pairs`` distinct tumours, and they set the informational activity columns.
    ``n_null_sets=0`` switches them off: such nulls are then ``none``. ``per_tumour`` is
    :func:`tumour_shifts` and ``per_mp`` :func:`protocol_shift`. Run time is dominated by the
    random-set κ (about ``n_null_sets`` x MPs x samples of these tumours evaluations of
    :func:`mpaudit.coexpr.kappa`) and the bootstrap (``KAPPA_BOOT`` x MPs x samples of paired tumours).
    """
    cohort = _prepare(adata, keys, reference, test, counts_layer=counts_layer, min_cells=min_cells,
                      max_cells=max_cells, min_detect=min_detect, seed=seed)
    levels = _levels(cohort, gene_sets, min_iota=min_iota, min_genes=min_genes, high_quantile=high_quantile,
                     seed=seed)
    null_levels = None
    if n_null_sets > 0:
        null_levels = _random_set_levels(cohort, gene_sets, list(gene_sets), n_null_sets, min_iota=min_iota,
                                         min_genes=min_genes, seed=seed)
    per_tumour, types = _shifts(levels, reference, test, null_levels, min_pairs, kappa_floor, kappa_min_effect,
                                kappa_active)
    return _summarise(per_tumour, types), per_tumour, levels
