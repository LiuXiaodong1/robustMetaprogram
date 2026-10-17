"""Toy tumour cohort with known ground truth, for tests and demonstrations.

Planted in malignant cells:

* biology: ``cycle`` and ``hypoxia`` (all tumours; with ``variable_prevalence`` the hypoxic
  fraction varies between samples, U(0.1, 0.5), so the sample-level hypoxia score does too),
  ``EMT`` (half the tumours; higher prevalence under the fixed protocol, i.e. selection against
  fragile cells in fresh samples);
* TME-shared biology: ``IFN`` (active in malignant *and* TME cells of the same samples); the only
  set in ``truth["expected_tme_shared"]``;
* private programs: two per tumour;
* genetic: a recurrent ``CNA`` (3x focal amplification of ~40 genes on one arm in a subclone of
  70% of tumours; a 2x gain at these expression levels was too weak for 3CA-style NMF to pick up);
* nuisance: dissociation ``stress`` (immediate-early/heat-shock genes; sample-level intensity,
  shared by malignant and TME cells, absent under the fixed protocol), ambient RNA from TME cells
  (per-cell fraction Beta(1, 20)), and malignant-TME doublets (4% of malignant cells). With
  ``doublet_cycle_coupling`` cycling cells have 1.3x depth and a doublet score shifted up by 0.3,
  so a real program is correlated with the doublet covariate (cycling cells are larger and called
  doublets more often).

TME cells (T, myeloid, fibroblast) have their own profiles: canonical markers at 40x and 30
TME-enriched genes per type at 10x (``truth["tme_enriched"]``). The TME-enriched genes are drawn
from all generic genes, so they can include planted program genes (never the technical-axis or
lineage genes).

Technical effects:

* ``malignant_ambient``: TME cells take up ambient RNA from their own sample's malignant cells
  (soup = mean malignant profile of the sample; per-cell fraction Beta(1, 25) times a per-sample
  factor U(0.5, 2)). A malignant program whose strength varies between samples (hypoxia, above)
  therefore leaks into the same samples' TME cells and can look TME-shared. Malignant cells also
  express a lineage gene set (``truth["malignant_markers"]``: 10x in malignant, 0.1x in TME cells;
  the analogue of epithelial markers in carcinoma), from which the malignant soup is estimable.
* ``technical_axes`` (all cells): an abundance axis (``abundance``: ~40 short genes plus the
  ribosomal genes, all set to comparable high levels as real RP genes are, together about a third
  of the library; scaled per cell by c ** ``CYTO_GAMMA`` = c ** 2, c ~ lognormal(0, 0.35) the
  planted cytoplasmic factor, so %ribo tracks it) and a gene-length axis (``long_genes``: ~40 genes
  of 150-800 kb scaled by 1 + ``LONG_K`` n = 1 + 20 n, n ~ Beta(2, 5) the nuclear/intronic
  fraction: exonic reads plus intronic reads, which these long genes collect at 21 times their
  exonic rate, in proportion to n; ``obs["intronic_fraction"]`` = n plus N(0, 0.03) noise). Both
  use genes outside the other planted sets and the CNA arm. The magnitudes were calibrated so that
  the standard pipeline can see the axes (below).
* ``protocol_capture_bias``: under the fixed protocol every gene's rate is multiplied by a
  per-gene capture efficiency (lognormal(0, 0.4), the same vector in all fixed samples). It moves
  gene-set intensities between protocols but not co-expression.

``var`` holds ``chromosome``, ``start``, ``gene_length`` (lognormal, median ~25 kb, clipped to
1 kb-2.5 Mb; drawn from the 2-8 kb and 150-800 kb ranges for the two axes) and ``end``.

Calibration of the technical axes. With the first magnitudes (c ** 1 and 1 + 2 n) the standard
pipeline could not see either axis: on ``simulate_cohort(seed=0)`` no per-sample NMF program, robust
or not, had more than 1 abundance or 6 long genes in its top 50, while each planted biology program
was robust in 9-16 of the 16 samples. The planted change was too small for log2(CPM/10 + 1): the
abundance genes make up a third of the library, so normalisation removes about a third of their
log-fold (and spreads it thinly over all other genes), and the long genes, drawn from all
expression levels, are mostly sparse, so counting noise swamps a 1.9-fold (p95/p5) change. The SD
across a sample's malignant cells of the planted log2 rate was 0.33 per abundance gene and 0.34 per
long gene, against 0.6-0.8 for the genes of cycle, hypoxia, EMT and IFN. The magnitudes were then
raised, with every random draw unchanged, against a target fixed before any audit was run: with the
reference pipeline at the demo settings (``GavishParams(nrun=10, min_genes_per_cell=0)``, recurrence
by tumour), each axis must put >= 10 of its genes into one robust program in >= 8 of the 16
samples. ``CYTO_GAMMA`` 2 is the smallest value tried that meets it (1.75: 5 samples); with it,
``LONG_K`` 10 does not (6 samples) and 14-24 do. At 2 and 20 the abundance and long-gene axes
appear in 13 and 11 samples, each forms its own MP, and the planted log2-rate SD is 0.65 and 0.80
per gene, in the range of the biology programs. The cost is competition for the fixed ranks 4-9:
robust appearances went from 16, 13, 9 and 10 samples to 16, 13, 7 and 9 for cycle, hypoxia, EMT and
IFN, from 11 to 8 for contamination and from 8 to 5 for the CNA; stress stayed at 9. Neighbouring
magnitudes move these counts by 2-3 samples either way; at ``LONG_K`` 14, EMT fell to 6 samples and
lost its MP. All counts are toy results at seed 0. Resulting magnitudes (malignant cells, medians
over samples): %ribo p5-p95 10-42% within a sample, SD of log %ribo 0.45 (previously 16-35% and
0.24); within-sample Spearman of the abundance score (mean log expression of its non-ribosomal
genes) with %ribo 0.90, against 0.68-0.77 in the user-reported breast cancer snRNA-seq cohort
(docs/05_pilot_P1.md; the toy axis is cleaner because only c and counting noise move %ribo here);
long-gene rate p95/p5 within a sample 5.7 (previously 1.9), the long genes' share of counts about 6%
(previously 1.6%), and Spearman of that share with ``intronic_fraction`` 0.84. Side effect: with a
third of the library on one factor, every other gene is diluted by the same per-cell factor, so
after total-count normalisation all genes co-vary a little: the latent correlation (``kappa``) of
random 30-gene sets of unplanted genes in malignant cells is about 0.23 (0.09 before; medians over
samples), more in cohorts with fewer genes. Seed-0 numbers derived from the default cohort before
the calibration (e.g. in tests and docstrings of other modules) changed with it. Mechanism tests of
other axes that should not depend on the technical axes (the TME spillover tests and the protocol
activity-gate test) build their cohorts at the first magnitudes (:func:`axis_magnitudes` with
``FIRST_MAGNITUDES``).

Covariates in ``obs``: ``clone``, ``doublet_score`` and ``contamination`` are TRUE values (the
planted clone; a score drawn from the true doublet status; the planted ambient fraction plus
N(0, 0.01), from TME cells for malignant cells and from malignant cells for TME cells).
``intronic_fraction`` is the planted nuclear fraction plus N(0, 0.03), a near-true readout standing
in for a covariate measured from reads (STARsolo Velocyto, velocyto or DropletQC), so recruitment of
the long genes by a decoy arm conditioned on it is true by construction. ``pct_ribo``, added by
:func:`mpaudit.data.add_qc_obs`, is computed from the counts. ``estimated_covariates`` adds what an
analyst would compute from counts, cell-type labels and gene positions only, never from the truth:
``clone_est`` (inferCNV-like two-group split on the arm whose genes co-vary most),
``contamination_est`` (DecontX-like: share of counts on genes specific to the other compartment,
relative to the soup) and ``doublet_est`` (expression-based: contamination estimate and library
size). They are crude on purpose: they err, and the doublet
estimate tracks cell size and hence cycling. Because TME-enriched genes can include planted
program genes, ``contamination_est`` (and ``doublet_est``, built from it) can also track a program
by gene overlap; whether and how much depends on the seed (compare ``truth["tme_enriched"]`` with
``truth["sets"]``). Planted per-cell states are in ``truth["cells"]``.

Outcomes of an audit on this cohort are true by construction. Use it to check that code runs
and that each axis responds as designed, not as evidence about real data. Each optional effect
either has its own random stream or makes its draws whatever its flag, and each block of counts
(sample x compartment) has its own stream, so a flag set to False gives a paired ablation: same
cells, program activities, clones (``cna=False`` merges them into one) and depths, and identical
counts wherever the effect does not act. With ``malignant_ambient``, ``variable_prevalence``,
``doublet_cycle_coupling``, ``technical_axes``, ``protocol_capture_bias`` and
``estimated_covariates`` off (``dissociation_stress``, ``tme_contamination`` and ``cna`` on) the
model is the original one; its random draws differ from versions before the flags were added, and
with few genes the programs are smaller because genes for the technical axes are set aside.
Turning off ``dissociation_stress``, ``tme_contamination``, ``cna`` and ``technical_axes`` instead
gives the null cohort for the P1 rule (see :func:`simulate_cohort`).
"""
from __future__ import annotations

import contextlib

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .positional import assign_arms, centromeres
from .typer import IEG_HSP_CORE

T_MARKERS = ["PTPRC", "CD3D", "CD3E", "CD2", "NKG7", "GZMA"]
MYELOID_MARKERS = ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "TYROBP", "FCER1G", "AIF1"]
FIB_MARKERS = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]
RIBO = [f"RPL{i}" for i in range(3, 13)] + [f"RPS{i}" for i in range(2, 12)]
MITO = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ND4", "MT-CYB", "MT-ATP6"]
AXIS_SIZE = 40  # genes on each technical axis (abundance, long genes)
CYTO_GAMMA = 2.0  # abundance-axis genes scale with cytoplasmic factor ** CYTO_GAMMA (calibrated)
LONG_K = 20.0  # long-axis genes scale with 1 + LONG_K * nuclear fraction (calibrated)
FIRST_MAGNITUDES = (1.0, 2.0)  # (CYTO_GAMMA, LONG_K) before the 2026-10-17 recalibration
N_LINEAGE = 10  # malignant lineage markers
SPECIFIC_FOLD = 8.0  # mean-CPM ratio defining compartment-specific genes for contamination_est


@contextlib.contextmanager
def axis_magnitudes(cyto_gamma: float, long_k: float):
    """Within the block, cohorts are simulated with ``CYTO_GAMMA = cyto_gamma`` and ``LONG_K = long_k``
    (restored on exit). ``axis_magnitudes(*FIRST_MAGNITUDES)`` gives the technical axes at their first
    magnitudes, for ablations and for mechanism tests of other axes that should not depend on them
    (no random draw changes, see :func:`simulate_cohort`)."""
    global CYTO_GAMMA, LONG_K
    old = CYTO_GAMMA, LONG_K
    CYTO_GAMMA, LONG_K = float(cyto_gamma), float(long_k)
    try:
        yield
    finally:
        CYTO_GAMMA, LONG_K = old


def _counts(rate: np.ndarray, depth: np.ndarray, theta: float, rng: np.random.Generator) -> np.ndarray:
    rate = rate / rate.sum(axis=1, keepdims=True)
    mu = depth[:, None] * rate
    lam = rng.gamma(theta, mu / theta)
    return rng.poisson(lam)


def _log_uniform(rng: np.random.Generator, lo: float, hi: float, n: int) -> np.ndarray:
    return np.exp(rng.uniform(np.log(lo), np.log(hi), n))


def _zscore(x: np.ndarray) -> np.ndarray:
    """Standardised ``x``; non-finite entries and constant input give 0."""
    out = np.zeros(x.shape)
    ok = np.isfinite(x)
    if ok.sum() >= 2 and x[ok].std() > 0:
        out[ok] = (x[ok] - x[ok].mean()) / x[ok].std()
    return out


def _otsu_upper(x: np.ndarray) -> np.ndarray:
    """Upper group of the two-group split of ``x`` maximising between-group variance (exact 1-D 2-means)."""
    n = x.size
    upper = np.zeros(n, bool)
    if n < 2:
        return upper
    order = np.argsort(x, kind="stable")
    xs = x[order]
    k = np.arange(1, n)
    cs = np.cumsum(xs)[:-1]
    between = k * (n - k) * (cs / k - (xs.sum() - cs) / (n - k)) ** 2
    between[xs[1:] <= xs[:-1]] = -1.0  # cut only between distinct values
    if between.max() > 0:
        upper[order[k[np.argmax(between)]:]] = True
    return upper


def _clone_split(counts: sp.csr_matrix, arm: np.ndarray, min_mean: float = 0.1, min_arm_genes: int = 5) -> np.ndarray:
    """inferCNV-like two-clone call for one sample's malignant cells.

    log1p CPM of genes with a known arm and mean count >= ``min_mean`` is centred per gene (sample
    mean) and per cell; per-arm means give a cells x arms matrix. The arm whose mean varies most
    across cells relative to what independent genes would give (variance of the arm mean over mean
    gene variance / number of genes; a stand-in for inferCNV's positional smoothing, without which
    arms with few genes win on noise) is split in two (Otsu). Labels: ``e2`` for the higher group,
    ``e1`` otherwise.
    """
    totals = np.maximum(np.asarray(counts.sum(axis=1)).ravel(), 1.0)
    keep = pd.notna(arm) & (np.asarray(counts.mean(axis=0)).ravel() >= min_mean)
    names, codes, sizes = np.unique(arm[keep].astype(str), return_inverse=True, return_counts=True)
    use = np.flatnonzero(sizes >= min_arm_genes)
    if counts.shape[0] < 2 or use.size == 0:
        return np.full(counts.shape[0], "e1", dtype=object)
    L = np.log1p((sp.diags(1e6 / totals) @ counts[:, np.flatnonzero(keep)]).toarray())
    L -= L.mean(axis=0)
    L -= L.mean(axis=1, keepdims=True)
    member = sp.csr_matrix((1.0 / sizes[codes], (np.arange(codes.size), codes)), shape=(codes.size, names.size))
    S = np.asarray(L @ member[:, use])
    independent = (L.var(axis=0) @ member[:, use]).ravel() / sizes[use]
    with np.errstate(divide="ignore", invalid="ignore"):
        excess = np.where(independent > 0, S.var(axis=0) / independent, 0.0)
    return np.where(_otsu_upper(S[:, int(np.argmax(excess))]), "e2", "e1").astype(object)


def _estimate_covariates(counts: sp.spmatrix, sample: np.ndarray, cell_type: np.ndarray,
                         arm: np.ndarray) -> pd.DataFrame:
    """``clone_est``, ``contamination_est`` and ``doublet_est`` from counts, labels and gene arms only.

    * ``contamination_est``: compartment-specific genes are defined cohort-wide from mean CPM
      (pseudo-count 1): TME-specific if the mean in some TME cell type exceeds the malignant mean
      ``SPECIFIC_FOLD``-fold, malignant-specific if the malignant mean exceeds the mean of every TME
      cell type that much. (Pooling TME types would lose markers of one type: ambient RNA alone
      caps their pooled ratio near 8.) Per cell, the share of counts on the other compartment's
      specific genes, divided by that share in the other compartment's cells of the same sample
      (cohort-wide if the sample has none), clipped to [0, 1]. The cell's own low expression of
      those genes biases it up. NaN when no gene qualifies (TME cells without
      ``malignant_ambient``, which plants the malignant markers).
    * ``doublet_est``: per sample and compartment, logistic(2 (z(contamination_est) + 0.5 z(log
      total counts) - 2)).
    * ``clone_est``: :func:`_clone_split` per sample on malignant cells; ``none`` for TME cells.
    """
    counts = sp.csr_matrix(counts, dtype=np.float64)
    n = counts.shape[0]
    totals = np.maximum(np.asarray(counts.sum(axis=1)).ravel(), 1.0)
    mal = cell_type == "Malignant"
    cpm = sp.diags(1e6 / totals) @ counts
    contamination = np.full(n, np.nan)
    if mal.any() and (~mal).any():
        mean_mal = np.asarray(cpm[mal].mean(axis=0)).ravel() + 1
        mean_tme = np.vstack([np.asarray(cpm[cell_type == ct].mean(axis=0)).ravel() + 1
                              for ct in pd.unique(cell_type[~mal])])
        tme_specific = mean_tme.max(axis=0) / mean_mal > SPECIFIC_FOLD
        mal_specific = mean_mal / mean_tme.max(axis=0) > SPECIFIC_FOLD
        for cells, specific in ((mal, tme_specific), (~mal, mal_specific)):
            if not specific.any():
                continue
            share = np.asarray(counts[:, np.flatnonzero(specific)].sum(axis=1)).ravel() / totals
            for s in pd.unique(sample[cells]):
                own, other = cells & (sample == s), ~cells & (sample == s)
                soup = share[other].mean() if other.any() else share[~cells].mean()
                if soup > 0:
                    contamination[own] = np.clip(share[own] / soup, 0.0, 1.0)

    clone = np.full(n, "none", dtype=object)
    doublet = np.zeros(n)
    log_total = np.log1p(totals)
    for s in pd.unique(sample):
        for grp in (mal & (sample == s), ~mal & (sample == s)):
            if grp.any():
                score = _zscore(contamination[grp]) + 0.5 * _zscore(log_total[grp])
                doublet[grp] = 1.0 / (1.0 + np.exp(-2.0 * (score - 2.0)))
        m = mal & (sample == s)
        if m.any():
            clone[m] = _clone_split(counts[np.flatnonzero(m)], arm)
    return pd.DataFrame({"clone_est": clone, "contamination_est": contamination, "doublet_est": doublet})


def simulate_cohort(
    n_tumours: int = 10,
    n_twin: int = 4,
    n_region: int = 2,
    n_mal: int = 200,
    n_tme_per_type: int = 40,
    n_genes: int = 1500,
    program_size: int = 40,
    n_studies: int = 3,
    seed: int = 0,
    *,
    malignant_ambient: bool = True,
    variable_prevalence: bool = True,
    doublet_cycle_coupling: bool = True,
    technical_axes: bool = True,
    protocol_capture_bias: bool = True,
    estimated_covariates: bool = True,
    dissociation_stress: bool = True,
    tme_contamination: bool = True,
    cna: bool = True,
):
    """Return ``(adata, truth)``.

    The flags switch the effects described in the module docstring. ``dissociation_stress``,
    ``tme_contamination`` (TME-derived ambient RNA and malignant-TME doublets in malignant cells) and
    ``cna`` switch off the planted nuisance programs; with these and ``technical_axes`` off the cohort
    holds only biology plus the confounders that can fake nuisance (spillover, doublet-cycle
    coupling, capture bias), which is the null cohort for the P1 rule. ``estimated_covariates=False``
    leaves out ``clone_est``, ``contamination_est`` and ``doublet_est``, and ``technical_axes=False``
    the ``abundance`` and ``long_genes`` sets (``intronic_fraction`` is then uninformative noise).
    ``intronic_fraction`` is always the planted nuclear fraction plus N(0, 0.03): a near-true
    readout standing in for a read-level measurement, not an estimate from the counts. The axis
    magnitudes ``CYTO_GAMMA`` and ``LONG_K`` are module constants calibrated for the default
    arguments (module docstring); they act only when ``technical_axes`` is on and change no random
    draw. ``truth`` holds the planted gene sets (``sets``) and their kinds (``kinds``), a cytoband
    table, the private programs, which tumours carry IFN/EMT/CNA, ``expected_tme_shared``,
    ``hypoxia_fraction_by_sample``, ``tme_ambient_scale_by_sample``, ``malignant_markers``,
    ``tme_enriched`` (TME cell type -> its 30 TME-enriched genes; may overlap planted programs),
    ``capture_efficiency`` (per gene; ones when off) and ``cells``: per-cell planted state (program
    activity, doublet, cytoplasmic factor c as drawn (abundance genes scale with c **
    ``CYTO_GAMMA``), nuclear fraction), indexed like ``adata.obs``.
    """
    import anndata as ad

    rng = np.random.default_rng(seed)
    # Separate streams: each optional effect has its own, and counts come from one stream per block
    # (sample x compartment), because Poisson draws consume a rate-dependent number of variates.
    # Switching an effect off thus keeps every structural draw (cells, programs, clones, depths)
    # and the counts of blocks the effect does not touch: flags give paired ablations.
    *effect_seeds, count_seeds = np.random.SeedSequence(seed).spawn(6)
    rng_len, rng_prev, rng_amb, rng_tech, rng_cap = map(np.random.default_rng, effect_seeds)
    named = list(dict.fromkeys(IEG_HSP_CORE + T_MARKERS + MYELOID_MARKERS + FIB_MARKERS + RIBO + MITO))
    generic = [f"G{i:04d}" for i in range(n_genes - len(named))]
    genes = np.array(named + generic)
    G = genes.size
    gi = {g: i for i, g in enumerate(genes)}

    # positions: 22 chromosomes of 100 Mb, centromere at 40-45 Mb; shuffle genes over the genome
    chrom = np.array([str(1 + (i % 22)) for i in rng.permutation(G)])
    start = rng.integers(1, 100_000_000, G)
    cyto = pd.DataFrame([(c, 40_000_000, 45_000_000, "cen", "acen") for c in map(str, range(1, 23))] +
                        [(c, 0, 40_000_000, "p1", "gneg") for c in map(str, range(1, 23))],
                        columns=["chrom", "start", "end", "band", "stain"])
    arm = np.where(start < 42_500_000, "p", "q")
    cna_arm_genes = np.flatnonzero((chrom == "7") & (arm == "q") & np.char.startswith(genes, "G"))

    base = rng.lognormal(-1.0, 1.5, G)
    base[[gi[g] for g in RIBO]] *= 30
    base[[gi[g] for g in MITO]] *= 20
    pool = rng.permutation([gi[g] for g in generic if gi[g] not in set(cna_arm_genes)])
    n_parts = 5 + 2 * n_tumours
    # genes for the technical axes and lineage markers come off the end of the pool, whatever the
    # flags; programs keep their original genes whenever the pool is large enough (it is at the defaults)
    reserve = min(2 * AXIS_SIZE + N_LINEAGE, pool.size // 4)
    prog_pool = pool[: pool.size - reserve]
    take = iter(np.array_split(prog_pool[: program_size * n_parts], n_parts))
    programs = {name: next(take) for name in ["cycle", "hypoxia", "EMT", "IFN", "stress_extra"]}
    stress_idx = np.concatenate([[gi[g] for g in IEG_HSP_CORE], programs.pop("stress_extra")[:16]])
    private = {f"T{t + 1:02d}": [next(take), next(take)] for t in range(n_tumours)}
    free = np.concatenate([prog_pool[program_size * n_parts:], pool[pool.size - reserve:]])

    tme_types = {"T": T_MARKERS, "Myeloid": ["PTPRC"] + MYELOID_MARKERS, "Fibroblast": FIB_MARKERS}
    tme_profile, tme_enriched = {}, []
    for ct, markers in tme_types.items():
        prof = base * rng.lognormal(0, 0.5, G)
        prof[[gi[g] for g in markers]] *= 40
        enriched = rng.choice(len(generic), 30, replace=False) + len(named)
        prof[enriched] *= 10
        tme_enriched.append(enriched)
        tme_profile[ct] = prof

    # technical-axis and lineage genes: free generic genes, not TME-enriched
    cand = np.setdiff1d(free, np.concatenate(tme_enriched))
    n_lin = min(N_LINEAGE, cand.size // 5)
    n_ax = min(AXIS_SIZE, (cand.size - n_lin) // 2)
    axes = rng_tech.choice(cand, 2 * n_ax, replace=False)
    abundance, long_genes = axes[:n_ax], axes[n_ax:]
    lineage = rng_amb.choice(np.setdiff1d(cand, axes), n_lin, replace=False)
    abund_idx = np.concatenate([abundance, [gi[g] for g in RIBO]])
    axis_spread = rng_tech.lognormal(-0.5**2 / 2, 0.5, abund_idx.size)  # mean 1

    gene_length = np.clip(rng_len.lognormal(np.log(25_000), 1.3, G), 1_000, 2_500_000)
    if technical_axes:
        gene_length[abund_idx] = _log_uniform(rng_len, 2_000, 8_000, abund_idx.size)
        gene_length[long_genes] = _log_uniform(rng_len, 150_000, 800_000, long_genes.size)
        # axis genes at comparable, high levels in every cell type (as the ~80 RP genes are), so that
        # no single gene's noise dominates %ribo or the axis score: abundance genes around the 95th
        # percentile of base expression, ribosomal genes around their previous mean level
        level = np.where(np.arange(abund_idx.size) < abundance.size, np.quantile(base, 0.95),
                         30 * np.exp(-1 + 1.5**2 / 2))
        factor = level * axis_spread / base[abund_idx]
        for prof in [base, *tme_profile.values()]:
            prof[abund_idx] *= factor
    gene_length = np.round(gene_length).astype(np.int64)
    if malignant_ambient:
        for prof in tme_profile.values():
            prof[lineage] *= 0.1
    capture = rng_cap.lognormal(0.0, 0.4, G) if protocol_capture_bias else np.ones(G)

    def technical(r: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Scale the technical axes in ``r`` in place; return cytoplasmic factor, nuclear fraction, its readout."""
        cyto_f = rng_tech.lognormal(0.0, 0.35, r.shape[0])
        nuc = rng_tech.beta(2.0, 5.0, r.shape[0])
        if technical_axes:  # magnitudes act on the draws, never change them (paired ablations, null cohort)
            r[:, abund_idx] *= (cyto_f ** CYTO_GAMMA)[:, None]
            r[:, long_genes] *= (1 + LONG_K * nuc)[:, None]
        return cyto_f, nuc, np.clip(nuc + rng_tech.normal(0, 0.03, r.shape[0]), 0, 1)

    tumours = [f"T{t + 1:02d}" for t in range(n_tumours)]
    t_ifn = set(rng.choice(tumours, int(round(0.6 * n_tumours)), replace=False))
    t_emt = set(rng.choice(tumours, int(round(0.5 * n_tumours)), replace=False))
    t_cna = set(rng.choice(tumours, int(round(0.7 * n_tumours)), replace=False))
    samples = [(t, "fresh", "R1") for t in tumours]
    samples += [(t, "fixed", "R1") for t in tumours[:n_twin]]
    samples += [(t, "fresh", "R2") for t in tumours[n_twin:n_twin + n_region]]

    blocks, obs_rows, cell_rows = [], [], []
    hyp_frac, amb_scale = {}, {}
    for t, protocol, region in samples:
        s_name = f"{t}_{protocol}_{region}"
        stress_level = 0.0 if protocol == "fixed" else rng.uniform(0.5, 2.0)
        if not dissociation_stress:  # draw anyway so the other effects keep their random streams
            stress_level = 0.0
        capture_on = protocol_capture_bias and protocol == "fixed"
        mal_base = base * rng.lognormal(0, 0.3, G)
        if malignant_ambient:
            mal_base[lineage] *= 10
        rate = np.tile(mal_base, (n_mal, 1))
        off = np.zeros(n_mal, bool)

        def activate(idx: np.ndarray, frac: float, fc: float) -> np.ndarray:
            on = rng.random(n_mal) < frac
            rate[np.ix_(np.flatnonzero(on), idx)] *= fc
            return on

        hyp_frac[s_name] = float(rng_prev.uniform(0.1, 0.5)) if variable_prevalence else 0.30
        cycling = activate(programs["cycle"], 0.25, 4.0)
        hypoxic = activate(programs["hypoxia"], hyp_frac[s_name], 3.0)
        emt = activate(programs["EMT"], 0.35 if protocol == "fixed" else 0.20, 3.0) if t in t_emt else off
        ifn = activate(programs["IFN"], 0.30, 3.0) if t in t_ifn else off
        for p in private[t]:
            activate(p, 0.30, 3.0)
        clone = np.where(rng.random(n_mal) < 0.5, "c2", "c1") if t in t_cna else np.full(n_mal, "c1")
        if not cna:
            clone = np.full(n_mal, "c1")
        if t in t_cna and cna:
            rate[np.ix_(clone == "c2", cna_arm_genes)] *= 3.0
        stress_act = stress_level * rng.gamma(2.0, 0.5, n_mal)
        rate[:, stress_idx] *= (1 + 3 * stress_act)[:, None]
        cyto_f, nuc, intronic = technical(rate)
        soup = (rate / rate.sum(axis=1, keepdims=True)).mean(axis=0)

        # doublets with TME cells (drawn whatever the flag, like every draw below)
        planted = rng.random(n_mal) < 0.04
        cts = list(tme_types)
        for r in np.flatnonzero(planted):
            prof = tme_profile[cts[rng.integers(len(cts))]]
            if tme_contamination:
                rate[r] = rate[r] / rate[r].sum() + prof / prof.sum()
        doublet = planted if tme_contamination else np.zeros(n_mal, bool)

        # ambient RNA from TME profiles
        amb = sum(tme_profile.values())
        amb = amb / amb.sum()
        rho = rng.beta(1.0, 20.0, n_mal)
        if not tme_contamination:
            rho = np.zeros(n_mal)
        rate = (1 - rho)[:, None] * rate / rate.sum(axis=1, keepdims=True) + rho[:, None] * amb[None, :]
        if capture_on:
            rate = rate * capture

        depth = rng.lognormal(np.log(4000), 0.5, n_mal)
        if doublet_cycle_coupling:
            depth[cycling] *= 1.3
        blocks.append(_counts(rate, depth, 10.0, np.random.default_rng(count_seeds.spawn(1)[0])))
        for c in range(n_mal):
            # The draw's shapes follow the planted status, not the flag: Beta sampling uses rejection,
            # so how many variates it takes depends on the shapes, and following the flag would shift
            # every later draw when tme_contamination is off. A doublet the flag removed gets
            # 1 - Beta(8, 2), which is Beta(2, 8) like any singlet's score.
            score = rng.beta(8, 2) if planted[c] else rng.beta(2, 8)
            if planted[c] and not doublet[c]:
                score = 1.0 - score
            if doublet_cycle_coupling and cycling[c] and not doublet[c]:
                score = min(score + 0.3, 1.0)
            obs_rows.append({"sample": s_name, "tumour": t, "protocol": protocol, "region": region,
                             "cell_type": "Malignant", "malignant": True, "clone": clone[c],
                             "doublet_score": float(score),
                             "contamination": float(np.clip(rho[c] + rng.normal(0, 0.01), 0, 1)),
                             "intronic_fraction": float(intronic[c])})
        cell_rows.append(pd.DataFrame({"cycle": cycling, "hypoxia": hypoxic, "EMT": emt, "IFN": ifn,
                                       "doublet": doublet, "cytoplasmic": cyto_f, "nuclear": nuc}))

        amb_scale[s_name] = float(rng_amb.uniform(0.5, 2.0)) if malignant_ambient else 0.0
        for ct, prof in tme_profile.items():
            r = np.tile(prof, (n_tme_per_type, 1))
            ifn_t = np.zeros(n_tme_per_type, bool)
            if t in t_ifn:
                ifn_t = rng.random(n_tme_per_type) < 0.3
                r[np.ix_(np.flatnonzero(ifn_t), programs["IFN"])] *= 3.0
            stress_t = stress_level * rng.gamma(2.0, 0.5, n_tme_per_type)
            r[:, stress_idx] *= (1 + 3 * stress_t)[:, None]
            cyto_t, nuc_t, intronic_t = technical(r)
            contamination = np.zeros(n_tme_per_type)
            if malignant_ambient:  # ambient RNA from this sample's malignant cells
                rho_t = np.clip(rng_amb.beta(1.0, 25.0, n_tme_per_type) * amb_scale[s_name], 0, 1)
                r = (1 - rho_t)[:, None] * r / r.sum(axis=1, keepdims=True) + rho_t[:, None] * soup[None, :]
                contamination = np.clip(rho_t + rng_amb.normal(0, 0.01, n_tme_per_type), 0, 1)
            if capture_on:
                r = r * capture
            depth = rng.lognormal(np.log(3000), 0.5, n_tme_per_type)
            blocks.append(_counts(r, depth, 10.0, np.random.default_rng(count_seeds.spawn(1)[0])))
            for i in range(n_tme_per_type):
                obs_rows.append({"sample": s_name, "tumour": t, "protocol": protocol, "region": region,
                                 "cell_type": ct, "malignant": False, "clone": "none",
                                 "doublet_score": float(rng.beta(2, 8)), "contamination": float(contamination[i]),
                                 "intronic_fraction": float(intronic_t[i])})
            cell_rows.append(pd.DataFrame({"cycle": False, "hypoxia": False, "EMT": False, "IFN": ifn_t,
                                           "doublet": False, "cytoplasmic": cyto_t, "nuclear": nuc_t}))

    obs = pd.DataFrame(obs_rows)
    obs["study"] = obs["tumour"].map({t: f"S{(i % n_studies) + 1}" for i, t in enumerate(tumours)})
    obs["cancer_type"] = obs["tumour"].map({t: ("A" if i % 2 == 0 else "B") for i, t in enumerate(tumours)})
    obs.index = [f"cell{i}" for i in range(len(obs))]
    X = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    var = pd.DataFrame({"chromosome": chrom, "start": start, "gene_length": gene_length, "end": start + gene_length},
                       index=pd.Index(genes, name="gene"))
    if estimated_covariates:
        est = _estimate_covariates(X, obs["sample"].to_numpy(), obs["cell_type"].to_numpy(),
                                   assign_arms(var, centromeres(cyto)).to_numpy(dtype=object))
        for col in est:
            obs[col] = est[col].to_numpy()
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()

    truth_sets = {name: list(genes[idx]) for name, idx in programs.items()}
    kinds = {"cycle": "biology", "hypoxia": "biology", "EMT": "biology (protocol-selected)",
             "IFN": "TME-shared biology"}
    if dissociation_stress:
        truth_sets["stress"] = list(genes[stress_idx])
        kinds["stress"] = "technical (dissociation)"
    if cna:
        truth_sets["CNA_7q"] = list(genes[cna_arm_genes])
        kinds["CNA_7q"] = "genetic (dosage)"
    if tme_contamination:
        truth_sets["contamination"] = list(dict.fromkeys(T_MARKERS + MYELOID_MARKERS + FIB_MARKERS))
        kinds["contamination"] = "technical (ambient/doublet)"
    if technical_axes:
        truth_sets["abundance"] = list(genes[abund_idx])
        truth_sets["long_genes"] = list(genes[long_genes])
        kinds["abundance"] = "technical (abundance / %ribo)"
        kinds["long_genes"] = "technical (gene length / nuclear fraction)"
    cells = pd.concat(cell_rows, ignore_index=True)
    cells.index = obs.index
    truth = {"sets": truth_sets, "kinds": kinds, "cytoband": cyto,
             "private": {t: [list(genes[p]) for p in ps] for t, ps in private.items()},
             "tumours_ifn": sorted(t_ifn), "tumours_emt": sorted(t_emt),
             "tumours_cna": sorted(t_cna) if cna else [],
             "expected_tme_shared": ["IFN"],
             "hypoxia_fraction_by_sample": hyp_frac,
             "tme_ambient_scale_by_sample": amb_scale,
             "malignant_markers": list(genes[lineage]) if malignant_ambient else [],
             "tme_enriched": {ct: list(genes[idx]) for ct, idx in zip(tme_types, tme_enriched)},
             "capture_efficiency": pd.Series(capture, index=genes, name="capture_efficiency"),
             "cells": cells}
    return adata, truth
