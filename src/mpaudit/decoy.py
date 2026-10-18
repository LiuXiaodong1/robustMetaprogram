"""Pilot P1b: nuisance-matched decoys ("how much can nuisance manufacture?").

A decoy keeps, per sample and per gene, the count distribution conditional on declared nuisance
(sequencing depth, categorical strata such as CNA clone, and continuous covariates such as a
doublet score or an ambient-contamination estimate) and destroys all other gene-gene
co-expression. Running the MP pipeline unchanged on decoys shows which MPs nuisance alone can
reproduce. (The pipeline is the Python port in :mod:`mpaudit.gavish`, not the 3CA R code.)

Two methods:

* ``"residual"`` (default): stratified conditional residual permutation (SCRP). Per gene a
  negative-binomial model with mean ``s_c * exp(a_{g,stratum} + b_g' z_c)`` (``s_c`` = size
  factor, ``z_c`` = standardised continuous covariates) gives randomised-quantile residuals;
  residuals are permuted within categorical strata and mapped back through each receiving cell's
  own distribution, so depth and covariate effects are kept cell by cell. The mean is a Poisson
  fit (Newton steps, halved per gene when they lower the objective) with a small ridge on the
  slopes ``b_g`` only: the stratum intercepts are free, so the fitted means of each stratum sum to
  its observed counts for every gene and low-count genes are not pulled up.
* ``"permute"``: raw counts permuted within strata (continuous columns are binned).

**Decoy arms.** Nuisance is declared in arms, each a set of obs columns conditioned on together;
each arm runs ``n_decoys`` decoys (decoy ``d`` of every arm uses seed ``seed + d``, so arms differ
only in what they condition on). Per real MP and arm, ``hits`` counts the decoys containing an MP
with Jaccard >= ``min_jaccard`` to it and ``recruitment = hits / n_decoys``; one match among three
decoys is not evidence (the P1 rule in :mod:`mpaudit.report` asks for >= 2 hits among >= 5
decoys). Reporting per arm types the recruitment: an arm of inferred CNA clones recruiting an MP
marks it **dosage** (genetic, not technical nuisance); an arm of technical covariates (doublet
score, ambient-contamination estimate, %ribo, intronic fraction) marks it **nuisance**. Recruitment
is reported per arm, not per covariate. The example config conditions on the four technical
covariates jointly in one arm, because each arm costs ``n_decoys`` (5) more pipeline passes, so a
technical hit is not attributed to a single covariate. To attribute it (e.g. to see whether
``doublet_est``, which tracks cell size and cycling, drives a hit), declare one arm per covariate
and list each in ``decision.technical_arms``. Criterion (a) then counts a hit in any of them, which
adds chances to hit, so fix the arm list in the registered config before running.

**kappa_retained** (P1-GO v1.3). A decoy reproduces an MP to the extent that the declared nuisance
explains the MP's own co-expression, beyond what every gene set of the sample shares. Per MP, arm,
decoy and malignant sample with >= ``min_cells`` cells (of nonzero total): ``G``, the MP's genes
detected in >= 5% of the sample's cells in the real data (at least 5); ``kappa_real``, the mean
latent correlation of ``G`` (:func:`mpaudit.coexpr.kappa`, size factors = total counts); and
``n_null_floor`` [50] random gene sets drawn from the sample's other eligible genes, matched gene by
gene to ``G`` on mean expression and detection (:func:`mpaudit.coexpr.matched_random_sets`), from a
stream of their own per (seed, sample, MP) (:func:`_set_rng`: an MP's sets do not depend on the
other MPs or samples of the call, or on their order). ``bg_real`` and ``s`` are the median and the SD
(ddof 1) of the sets' kappa. The sample enters when >= ``MIN_FLOOR_SETS`` [20] sets have a finite
kappa and ``kappa_real - bg_real > max(min_kappa [0.02], EXCESS_ENTRY_SD [3] x s)``; its ratio is
``(kappa_decoy - bg_decoy) / (kappa_real - bg_real)``, the same genes' and the same sets' kappa on
the decoy's counts of the same cells (``bg_decoy``: the median over the sets with a finite kappa).
``kappa_retained`` is the median ratio over the entering samples, averaged over the arm's decoys;
NaN when no sample enters, and the MP then can neither meet criterion (a) in that arm nor be typed
dosage by it. Every kappa of a sample on one count matrix comes from one
:func:`mpaudit.coexpr.set_kappas` call (the values of :func:`~mpaudit.coexpr.kappa` set by set; on
the default toy, 8 sets, it halved the time of a ``kappa_retained`` call, 22 s to 11 s on one core).
It is not clipped: noise moves single ratios below 0 or above 1, and see below for systematic ratios
above 1.

Why the excess: every gene set of a sample co-varies a little (shared nuisance such as a per-cell
dilution or ambient RNA), and a decoy keeps that background. v1.2's ratio ``kappa_decoy /
kappa_real`` credited the decoy with it: under v1.2 the technical arm kept 0.12-0.19 of the cycle,
hypoxia, EMT and IFN sets' kappa on the default toy (mean of 3 seeds, v1.2 reference of
``scripts/calibration``), and a decoy that keeps only the background of a hand-built sample has a
ratio of about 0.4 and an excess ratio of about 0 (``tests/test_decoy_arms.py``).
Why the 3-SD entry: v1.2's entry (``kappa_real`` above the q95 of 20 random sets) let in samples
where the MP is absent, which carry the nuisance's background but not the program and have high
ratios: with EMT kept in 2 carrier tumours plus the 8 tumours without it (minority cohort of
``scripts/calibration/evaluate_a.py``, seed 2) it reached 0.66 under the technical arm. The 3-SD
entry keeps the denominator's coefficient of variation at most 1/3, and its residual selection bias
is downward, which is conservative for (a); reference results with it are in
``scripts/calibration/README.md``. On the default toy (below) no sample of a tumour without EMT,
IFN or CNA_7q enters for that set. The price: weakly co-expressed MPs lose entering samples
(contamination enters in 7 of the 16 samples of that toy) and can become NaN, i.e. unpowered in the
arm; :func:`decoy_recruitment` reports the entering samples (``kappa_retained_n_samples``) and
``summary.md`` lists the MPs with none in an arm of the rule.

**Caveat: covariates computed from expression.** Doublet scores, contamination estimates, %ribo,
``pct_long``, the intronic fraction and inferred clones are computed from the same counts (the
intronic fraction from the same reads). A decoy conditioned on one keeps whatever co-expression the
covariate predicts, including that of a *real* program correlated with it: cycling cells are larger
and get higher doublet scores, so a doublet arm keeps part of the cycle program, and with enough
cells NMF can find that part in decoys; %ribo tracks translation, ribosome biogenesis and MYC
activity, so an arm on it keeps such programs, and an MP containing ribosomal-protein genes shares
counts with %ribo. Recruitment alone would then call real biology nuisance; ``kappa_retained``
grades how much of the excess co-expression the covariate explains. On the toy cohort
(:func:`mpaudit.sim.simulate_cohort`, seed 0, 16 samples; decoys with seeds 1000 and 1001,
``kappa_retained`` seed 1000 as in :func:`decoy_recruitment`, mean of the two decoys; re-run for
v1.3 on 2026-10-17) an arm on the planted doublet score (cycling cells score 0.3 higher;
within-sample Spearman 0.61 with cycling) keeps 0.42 of the cycle program's excess, an arm on the
expression-based ``doublet_est`` (Spearman 0.26 on this cohort) 0.00; ``clone_est`` keeps 1.00 of
the CNA set's excess. On the null cohort of the same seed (no TME ambient RNA, so
``contamination_est``, from which ``doublet_est`` is built, tracks cycling), ``doublet_est`` has
Spearman 0.51 with cycling and its arm keeps 0.26 of cycle's excess (re-run 2026-10-18). Back on
the default cohort, an arm on ``pct_ribo`` keeps 0.87 of the abundance set's excess, 0.06 of the
long genes' and between -0.07 and 0.04 of cycle, hypoxia, EMT, IFN and CNA; the four-covariate
technical arm of the demo keeps 0.88 of abundance, 0.92 of long genes and 1.39 of contamination,
and between -0.04 and 0.07 of those five biology sets. Two limits remain:

* a real program strongly tied to a covariate's estimator (cell size and the doublet score, a
  transcriptional state that an expression-based clone caller splits on, translation and %ribo)
  keeps most of its co-expression and cannot be told from nuisance here; if the decoys recruit it,
  it counts as nuisance support under criterion (a) (in ``scripts/calibration``, under the technical
  arm with the planted doublet score in place of ``doublet_est``, cycle keeps 0.497, 0.532 and 0.530
  on the default toy at seeds 0, 1 and 2, at the 0.5 guard, and 0.64-0.67 on the null toy, above
  it). Only covariates measured independently of
  expression (genotype demultiplexing, cell hashing, protocol twins) settle such cases;
* a covariate computed from the MP's own genes inflates ``kappa_retained`` above 1: the covariate
  carries those genes' sampling noise as a shared factor, which the decoy turns into
  co-expression, and the excess ratio, with a smaller denominator than v1.2's ratio, inflates it
  further. ``contamination_est`` is the share of counts on TME-specific genes, from which the
  contamination set is drawn, and ``doublet_est`` is built from it: on the toy above an arm on
  ``contamination_est`` keeps 1.97 of the contamination set's excess and one on ``doublet_est`` 1.75.
  ``pct_ribo`` and ``pct_long`` (shares of counts on ribosomal-protein and on long genes,
  :func:`mpaudit.data.qc_metrics`) do the same for MPs made of those genes. Estimate covariates from
  genes outside the audited MP where possible.

Strata leak a little real signal: per-stratum gene means carry the stratum's share of every
program, so decoy co-expression grows roughly with (number of strata / cells per sample). Depth
is therefore modelled continuously through ``s_c`` by default (``n_depth_bins=1``), continuous
columns enter the mean model instead of being binned, and binned strata are capped so each has
about ``min_stratum_cells`` cells. Categorical strata are used as given; report their sizes.

This is an audit of pipelines that do not model nuisance. It carries no FDR claim.
"""
from __future__ import annotations

import hashlib
import inspect
import zlib

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from .coexpr import matched_random_sets, set_kappas
from .compare import jaccard_matrix
from .data import Keys, get_counts, malignant_mask, sample_map
from .gavish import DECOY_TAG, GavishParams, NMFCache, assemble, compute_programs

PER_MP_COLUMNS = ["mp", "arm", "n_decoys", "hits", "recruitment", "mean_best_jaccard", "max_best_jaccard",
                  "kappa_retained", "kappa_retained_n_samples"]
# n_fits: NMF runs behind the decoy cohort's programs (samples x ranks x nrun); n_capped_fits: of which stopped
# at max_iter; n_failed_fits: of which stopped by the restart limit (see mpaudit.gavish)
PER_DECOY_COLUMNS = ["arm", "decoy", "n_robust", "n_mps", "n_fits", "n_capped_fits", "n_failed_fits"]


def _is_continuous(v: pd.Series) -> bool:
    return pd.api.types.is_numeric_dtype(v) and not pd.api.types.is_bool_dtype(v) and v.nunique() > 10


def _covariates(obs: pd.DataFrame, cols: tuple[str, ...]) -> np.ndarray:
    """Standardised continuous covariates (cells x k); constant columns are dropped."""
    out = []
    for col in cols:
        v = pd.to_numeric(obs[col], errors="coerce").to_numpy(dtype=float)
        v = np.where(np.isfinite(v), v, np.nanmedian(v) if np.isfinite(v).any() else 0.0)
        sd = v.std()
        if sd > 0:
            out.append((v - v.mean()) / sd)
    return np.column_stack(out) if out else np.zeros((len(obs), 0))


def _strata(obs: pd.DataFrame, totals: np.ndarray, strata: tuple[str, ...], n_depth_bins: int, n_cont_bins: int,
            min_stratum_cells: int = 50) -> np.ndarray:
    cap = max(1, len(totals) // max(min_stratum_cells, 1))
    n_depth_bins, n_cont_bins = min(n_depth_bins, cap), min(n_cont_bins, cap)
    labels = []
    if n_depth_bins > 1:
        labels.append(pd.qcut(np.log1p(totals), n_depth_bins, labels=False, duplicates="drop").astype(str))
    for col in strata:
        v = obs[col]
        if pd.api.types.is_numeric_dtype(v) and v.nunique() > 10:
            if n_cont_bins < 2:
                continue
            labels.append(pd.qcut(v.rank(method="first"), n_cont_bins, labels=False, duplicates="drop").astype(str).to_numpy())
        else:
            labels.append(v.astype(str).to_numpy())
    if not labels:
        return np.zeros(len(totals), dtype=int)
    joined = pd.Series(["|".join(t) for t in zip(*labels)])
    return pd.factorize(joined)[0]


def _permute_within(U: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = U.copy()
    for h in np.unique(strata):
        rows = np.flatnonzero(strata == h)
        if rows.size < 2:
            continue
        perm = np.argsort(rng.random((rows.size, U.shape[1])), axis=0)
        out[rows] = np.take_along_axis(U[rows], perm, axis=0)
    return out


def _nb_dispersion(Y: np.ndarray, mu: np.ndarray) -> np.ndarray:
    num = (mu**2).sum(axis=0)
    den = ((Y - mu) ** 2 - mu).sum(axis=0)
    theta = np.where(den > 0, num / np.maximum(den, 1e-12), np.inf)
    return np.clip(theta, 1e-2, np.inf)


def _glm_mean(Y: np.ndarray, s: np.ndarray, strata: np.ndarray, Z: np.ndarray, n_iter: int = 30,
              ridge: float = 1e-3, tol: float = 1e-10, max_halvings: int = 20) -> np.ndarray:
    """Per-gene Poisson log-linear fit ``mu = s * exp(a_stratum + Z b)``, vectorised over genes.

    Maximises ``sum(Y log mu - mu) - ridge * n * |b|^2 / 2`` per gene: only the covariate slopes
    are penalised, so at the optimum the fitted means of each stratum sum to its observed counts
    (the intercept score equation). Newton steps (IRLS for the canonical link) are halved, per
    gene, until they do not lower the objective; plain IRLS overshoots on genes driven by a few
    high-leverage cells (e.g. doublets at the far end of a contamination covariate) and can return
    means orders of magnitude above the data. Stops when no gene's objective rises by more than
    ``tol`` (relative), or after ``n_iter`` steps.
    """
    levels, codes = np.unique(strata, return_inverse=True)
    X = np.hstack([np.eye(levels.size)[codes], Z])
    n, p = X.shape
    logs = np.log(s)[:, None]
    pen = np.r_[np.zeros(levels.size), np.full(Z.shape[1], ridge * n)]

    def objective(B, cols=slice(None)):
        eta = np.clip(X @ B + logs, -30, 30)
        return (Y[:, cols] * eta - np.exp(eta)).sum(axis=0) - 0.5 * (pen[:, None] * B**2).sum(axis=0), eta

    B = np.zeros((p, Y.shape[1]))
    for h in range(levels.size):
        rows = codes == h
        B[h] = np.log(np.maximum(Y[rows].sum(axis=0), 0.5) / max(s[rows].sum(), 1e-12))
    obj, eta = objective(B)
    for _ in range(n_iter):
        mu = np.exp(eta)
        grad = X.T @ (Y - mu) - pen[:, None] * B
        H = np.einsum("np,ng,nq->gpq", X, mu, X, optimize=True) + np.diag(pen)
        step = np.linalg.solve(H, grad.T[..., None])[..., 0].T
        new_B, new_obj, new_eta = B.copy(), obj.copy(), eta.copy()
        todo, t = np.arange(Y.shape[1]), 1.0
        for _ in range(max_halvings):
            cand = B[:, todo] + t * step[:, todo]
            o, e = objective(cand, todo)
            ok = o >= obj[todo] - 1e-12 * np.abs(obj[todo])  # rounding near the optimum
            acc = todo[ok]
            new_B[:, acc], new_obj[acc], new_eta[:, acc] = cand[:, ok], o[ok], e[:, ok]
            todo, t = todo[~ok], t / 2
            if not todo.size:
                break
        gain = new_obj - obj
        B, obj, eta = new_B, new_obj, new_eta
        if np.all(gain <= tol * (1.0 + np.abs(obj))):
            break
    return np.exp(eta)


def _scrp_block(Y: np.ndarray, s: np.ndarray, strata: np.ndarray, rng: np.random.Generator,
                Z: np.ndarray | None = None) -> np.ndarray:
    n, g = Y.shape
    if Z is not None and Z.shape[1]:
        mu = _glm_mean(Y, s, strata, Z)
    else:
        mu = np.empty_like(Y)
        for h in np.unique(strata):
            rows = strata == h
            p = Y[rows].sum(axis=0) / max(s[rows].sum(), 1e-12)
            mu[rows] = s[rows, None] * p[None, :]
    theta = _nb_dispersion(Y, mu)
    pois = ~np.isfinite(theta) | (theta > 1e6)
    th = np.where(pois, 1.0, theta)[None, :]
    prob = th / (th + mu)

    def apply(x, poisson_fn, nbinom_fn):  # each gene through its own family only (scipy's are slow)
        out = np.empty(x.shape)
        out[:, pois] = poisson_fn(x[:, pois], mu[:, pois])
        out[:, ~pois] = nbinom_fn(x[:, ~pois], th[:, ~pois], prob[:, ~pois])
        return out

    def cdf(x):
        return np.where(x < 0, 0.0, apply(x, stats.poisson.cdf, stats.nbinom.cdf))

    lo, hi = cdf(Y - 1), cdf(Y)
    U = lo + rng.random(Y.shape) * (hi - lo)
    U = np.clip(_permute_within(U, strata, rng), 0.0, 1.0 - 1e-12)
    out = apply(U, stats.poisson.ppf, stats.nbinom.ppf)
    out = np.nan_to_num(out, nan=0.0, posinf=0.0)
    return np.maximum(out, 0.0)


def make_decoy(
    adata,
    keys: Keys,
    *,
    counts_layer: str | None = None,
    strata: tuple[str, ...] = (),
    n_depth_bins: int = 1,
    n_cont_bins: int = 3,
    min_stratum_cells: int = 50,
    method: str = "residual",
    seed: int = 0,
    gene_chunk: int = 2000,
):
    """Return a copy of ``adata`` whose malignant-cell counts are replaced by a decoy.

    ``strata`` lists obs columns to condition on. With ``method="residual"``, numeric columns with
    more than 10 distinct values enter the mean model as continuous covariates; the others define
    permutation strata. The choice is made per sample, on the sample's malignant cells. Missing
    values are not modelled: a continuous column's NaN take the sample's median, a categorical
    column's NaN form a stratum of their own, and a numeric column left with 10 or fewer distinct
    values in a sample becomes a set of strata there, where a cell alone in its stratum keeps its
    real counts (with ``method="permute"``, NaN fall into their own bin).
    :func:`mpaudit.pipeline.validate_config` therefore refuses decoy arm columns with NaN in
    malignant cells.

    The draws come from one stream (``seed``) taken over the samples and their cells in the
    cohort's order, so the same cells in another order give another decoy, a different draw of the
    same law (the random gene sets of :func:`kappa_retained` and of :mod:`mpaudit.protocol` do not
    depend on order). A rerun reproduces the decoys only from the same cohort file, whose checksum
    ``run_config.json`` does not record: keep it with the run's outputs.
    """
    if method not in ("residual", "permute"):
        raise ValueError("method must be 'residual' or 'permute'")
    rng = np.random.default_rng(seed)
    X = get_counts(adata, counts_layer)
    mal = malignant_mask(adata, keys)
    samples = adata.obs[keys.sample].astype(str).to_numpy()
    keep_rows = np.ones(X.shape[0])
    add_r, add_c, add_v = [], [], []
    for s in pd.unique(samples[mal]):
        idx = np.flatnonzero(mal & (samples == s))
        if idx.size < 2:
            continue
        Xs = X[idx]
        totals = np.asarray(Xs.sum(axis=1)).ravel()
        obs_s = adata.obs.iloc[idx]
        if method == "residual":
            cont = tuple(c for c in strata if _is_continuous(obs_s[c]))
            cat = tuple(c for c in strata if c not in cont)
            st = _strata(obs_s, totals, cat, n_depth_bins, n_cont_bins, min_stratum_cells)
            Z = _covariates(obs_s, cont)
        else:
            st = _strata(obs_s, totals, strata, n_depth_bins, n_cont_bins, min_stratum_cells)
            Z = None
        sf = totals / max(np.median(totals), 1e-12)
        sf[sf == 0] = 1e-6
        expressed = np.flatnonzero(np.asarray(Xs.sum(axis=0)).ravel() > 0)
        keep_rows[idx] = 0.0
        for start in range(0, expressed.size, gene_chunk):
            cols = expressed[start:start + gene_chunk]
            Y = Xs[:, cols].toarray()
            block = _permute_within(Y, st, rng) if method == "permute" else _scrp_block(Y, sf, st, rng, Z)
            coo = sp.coo_matrix(block)
            add_r.append(idx[coo.row]); add_c.append(cols[coo.col]); add_v.append(coo.data)
    X_new = (sp.diags(keep_rows) @ X).tocsr()
    if add_r:
        X_new = X_new + sp.csr_matrix((np.concatenate(add_v), (np.concatenate(add_r), np.concatenate(add_c))),
                                      shape=X.shape)
    out = adata.copy()
    X_new = X_new.tocsr()
    if counts_layer:
        out.layers[counts_layer] = X_new
    else:
        out.X = X_new
    return out


# --------------------------------------------------------------------------- kappa_retained

N_NULL_FLOOR = 50  # matched random gene sets per (sample, set): the background's median and SD (decoy.n_null_floor)
MIN_FLOOR_SETS = 20  # random sets with a finite kappa needed for a sample to enter
EXCESS_ENTRY_SD = 3.0  # a sample enters when kappa_real - bg_real > max(min_kappa, EXCESS_ENTRY_SD x SD of the sets)
MIN_CELLS = 20  # malignant cells (of nonzero total) a sample needs for kappa_retained
MIN_DETECT = 0.05  # G: the set's genes detected in >= this share of the sample's cells (real data) ...
MIN_GENES = 5  # ... at least this many of them


def _set_rng(seed: int, sample: str, mp: str) -> np.random.Generator:
    """Random-set stream of one (sample, set), ``default_rng([seed, crc32(sample + '\\x1f' + set)])``: a set's
    random sets do not depend on the other sets or samples of the call, or on their order."""
    return np.random.default_rng([int(seed), zlib.crc32(f"{sample}\x1f{mp}".encode())])


def _reference_kappa(adata, keys: Keys, gene_sets: dict[str, list[str]], counts_layer: str | None,
                     min_cells: int, min_detect: float, min_genes: int, min_kappa: float = 0.02,
                     n_null_floor: int = N_NULL_FLOOR, seed: int = 0) -> list[tuple[str, np.ndarray, dict]]:
    """Per malignant sample with >= ``min_cells`` cells (of nonzero total): ``(sample, rows, {set: info})``.

    A set's genes ``G`` are those detected in >= ``min_detect`` of the sample's cells in the real
    data; samples where fewer than ``min_genes`` qualify are left out for that set. ``n_null_floor``
    random gene sets are drawn from the sample's other eligible genes, matched gene by gene to ``G``
    on mean expression and detection (:func:`mpaudit.coexpr.matched_random_sets`), from the (sample,
    set)'s own stream (:func:`_set_rng`). ``info``: ``genes`` (``G``), ``sets`` (the random sets with a
    finite kappa), ``kappa_sets`` (their kappa), ``k_real`` (kappa of ``G``), ``bg_real`` and
    ``sd_real`` (median and SD, ddof 1, of ``kappa_sets``; NaN with fewer than ``MIN_FLOOR_SETS``) and
    ``enter`` (``k_real - bg_real > max(min_kappa, EXCESS_ENTRY_SD * sd_real)``). Every kappa of a
    sample comes from one :func:`mpaudit.coexpr.set_kappas` call. The decoy is scored on the same
    cells, genes and random sets.
    """
    if n_null_floor < MIN_FLOOR_SETS:
        raise ValueError(f"n_null_floor must be >= {MIN_FLOOR_SETS} (a sample enters only with that many random "
                         f"sets of finite kappa), got {n_null_floor}")
    X = get_counts(adata, counts_layer)
    mal = malignant_mask(adata, keys)
    samples = adata.obs[keys.sample].astype(str).to_numpy()
    pos = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    idx = {m: np.array(sorted({pos[g] for g in genes if g in pos}), dtype=int) for m, genes in gene_sets.items()}
    refs = []
    for s in pd.unique(samples[mal]):
        rows = np.flatnonzero(mal & (samples == s))
        totals = np.asarray(X[rows].sum(axis=1)).ravel()
        rows, totals = rows[totals > 0], totals[totals > 0]
        if rows.size < min_cells:
            continue
        Xs = X[rows]
        detect = np.asarray((Xs > 0).mean(axis=0)).ravel()
        mean_expr = np.asarray(Xs.sum(axis=0)).ravel() / totals.sum()
        eligible = detect >= min_detect
        drawn = {}
        for m, g in idx.items():
            g = g[eligible[g]]
            if g.size >= min_genes:
                drawn[m] = (g, matched_random_sets(g, mean_expr, detect, n_null_floor, _set_rng(seed, s, m),
                                                   candidates=eligible))
        k = iter(set_kappas(Xs.tocsc(), totals, [x for g, rs in drawn.values() for x in (g, *rs)]))
        per_set = {}
        for m, (g, rs) in drawn.items():
            k_real = next(k)
            b = np.array([next(k) for _ in rs])
            ok = np.isfinite(b)
            info = {"genes": g, "sets": [r for r, o in zip(rs, ok) if o], "kappa_sets": b[ok], "k_real": k_real,
                    "bg_real": np.nan, "sd_real": np.nan, "enter": False}
            if ok.sum() >= MIN_FLOOR_SETS:
                info["bg_real"], info["sd_real"] = float(np.median(b[ok])), float(np.std(b[ok], ddof=1))
                with np.errstate(invalid="ignore"):
                    info["enter"] = bool(k_real - info["bg_real"] > max(min_kappa, EXCESS_ENTRY_SD * info["sd_real"]))
            per_set[m] = info
        refs.append((s, rows, per_set))
    return refs


def _retained_ratios(refs: list[tuple[str, np.ndarray, dict]], decoy_counts: sp.csr_matrix) -> dict[str, list[float]]:
    """Per set, over the entering samples: ``(kappa_decoy - bg_decoy) / (kappa_real - bg_real)``, ``bg`` the
    median kappa of the sample's random sets (those of finite real kappa) on the decoy's or the real counts."""
    out: dict[str, list[float]] = {}
    for _, rows, per_set in refs:
        enter = {m: info for m, info in per_set.items() if info["enter"]}
        if not enter:
            continue
        D = decoy_counts[rows]
        totals = np.asarray(D.sum(axis=1)).ravel()
        k = iter(set_kappas(D.tocsc(), totals, [x for info in enter.values() for x in (info["genes"], *info["sets"])]))
        for m, info in enter.items():
            k_dec = next(k)
            b = np.array([next(k) for _ in info["sets"]])
            b = b[np.isfinite(b)]
            if not b.size:
                continue
            ratio = (k_dec - float(np.median(b))) / (info["k_real"] - info["bg_real"])
            if np.isfinite(ratio):
                out.setdefault(m, []).append(float(ratio))
    return out


def kappa_retained(
    adata,
    decoy,
    keys: Keys,
    gene_sets: dict[str, list[str]],
    *,
    counts_layer: str | None = None,
    min_cells: int = MIN_CELLS,
    min_kappa: float = 0.02,
    min_detect: float = MIN_DETECT,
    min_genes: int = MIN_GENES,
    n_null_floor: int = N_NULL_FLOOR,
    seed: int = 0,
) -> pd.DataFrame:
    """Share of each gene set's co-expression, in excess of matched random gene sets, that one decoy keeps
    (no NMF; see the module docstring).

    ``decoy`` is :func:`make_decoy` output for ``adata``. Columns: ``mp`` (set name),
    ``kappa_retained`` (median over the entering samples of ``(kappa_decoy - bg_decoy) / (kappa_real -
    bg_real)``, ``bg`` the median kappa of ``n_null_floor`` matched random gene sets of the sample on
    the same counts; a sample enters when ``kappa_real - bg_real > max(min_kappa, EXCESS_ENTRY_SD x SD
    of the sets' real kappa)``, with >= ``MIN_FLOOR_SETS`` sets of finite kappa; NaN when no sample
    enters) and ``n_samples`` (samples entering the median). ``seed`` seeds the random sets
    (:func:`_set_rng`); :func:`decoy_recruitment` passes its own ``seed`` (1000 in the pipeline).
    """
    refs = _reference_kappa(adata, keys, gene_sets, counts_layer, min_cells, min_detect, min_genes,
                            min_kappa, n_null_floor, seed)
    ratios = _retained_ratios(refs, get_counts(decoy, counts_layer))
    return pd.DataFrame([{"mp": m, "kappa_retained": float(np.median(ratios[m])) if m in ratios else np.nan,
                          "n_samples": len(ratios.get(m, []))} for m in gene_sets])


# --------------------------------------------------------------------------- recruitment

def decoy_tag(adata, keys: Keys, **make_decoy_kwargs) -> str:
    """NMF cache tag of the decoy ``make_decoy(adata, keys, **make_decoy_kwargs)``.

    Readable part: :data:`mpaudit.gavish.DECOY_TAG`, the method, the conditioning columns, the depth
    bins and the seed. Then a hash of every argument of :func:`make_decoy` (its defaults included, so
    passing another argument, or changing a default, gives another tag), of the sample and malignant
    keys, of the malignant cells' barcodes in order, and of the dtype and the values in malignant
    cells of the sample column and of each conditioning column (:func:`make_decoy` treats a numeric
    column and a categorical one with the same values differently). Only malignant cells enter, as
    in :func:`make_decoy`, so adding or removing a TME cell keeps the tag. The cache key adds the
    fingerprint of the decoy counts themselves (:func:`mpaudit.gavish.compute_programs`), which
    changes with any change to them (also with ``method="permute"``, which keeps every total), and
    the decoy generator's code (:func:`mpaudit.gavish.code_version`).
    """
    bound = inspect.signature(make_decoy).bind(adata, keys, **make_decoy_kwargs)
    bound.apply_defaults()
    args = {k: v for k, v in bound.arguments.items() if k not in ("adata", "keys")}
    mal = malignant_mask(adata, keys)
    h = hashlib.sha1(repr(sorted((k, repr(v)) for k, v in args.items())).encode())
    h.update(repr((keys.sample, keys.malignant, keys.malignant_label)).encode())
    h.update("\x1f".join(np.asarray(adata.obs_names[mal]).astype(str)).encode())
    for col in (keys.sample, *args["strata"]):
        h.update(f"{col}\x1f{adata.obs[col].dtype}".encode())
        h.update(pd.util.hash_pandas_object(adata.obs[col][mal], index=False).to_numpy().tobytes())
    return (f"{DECOY_TAG}{args['method']}-{'+'.join(args['strata'])}-d{args['n_depth_bins']}-s{args['seed']}-"
            f"{h.hexdigest()[:12]}")


def decoy_recruitment(
    adata,
    keys: Keys,
    params: GavishParams,
    real_mps: dict[str, list[str]],
    *,
    arms: dict[str, tuple[str, ...]] | None = None,
    strata: tuple[str, ...] = (),
    n_decoys: int = 5,
    n_depth_bins: int = 1,
    method: str = "residual",
    counts_layer: str | None = None,
    cache: NMFCache | None = None,
    recurrence_key: str | None = None,
    min_jaccard: float = 0.3,
    min_cells: int = MIN_CELLS,
    min_kappa: float = 0.02,
    n_null_floor: int = N_NULL_FLOOR,
    seed: int = 1000,
    n_jobs: int = 1,
    log=print,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict[str, list[str]]]]:
    """Run the MP pipeline on ``n_decoys`` decoy cohorts per arm.

    ``arms`` maps an arm name to the obs columns its decoys condition on (e.g. ``{"technical":
    ("doublet_est", "contamination_est"), "dosage": ("clone_est",)}``); ``strata`` alone defines
    one arm named ``"all"``. An arm's columns are used sorted and de-duplicated, and the NMF cache
    tag of a decoy (:func:`decoy_tag`) is built from the arguments of its :func:`make_decoy` call
    and the columns' values, so arms with the same columns share fits.

    Returns ``(per_mp, per_decoy, decoy_mps)``:

    * ``per_mp``, one row per (MP, arm): ``mp``, ``arm``, ``n_decoys``, ``hits`` (decoys with an MP
      at Jaccard >= ``min_jaccard`` to the real MP), ``recruitment`` (``hits / n_decoys``),
      ``mean_best_jaccard``, ``max_best_jaccard`` (best match per decoy, over decoys),
      ``kappa_retained`` (``min_cells``, ``min_kappa``, ``n_null_floor``: see the module docstring) and
      ``kappa_retained_n_samples`` (samples entering it, decided on the real counts, so the same in
      every arm; with 0, ``kappa_retained`` is NaN and the MP can neither meet criterion (a) nor be
      typed dosage in the arm);
    * ``per_decoy``: ``arm``, ``decoy`` (cache tag), ``n_robust``, ``n_mps``, and the NMF fit
      diagnostics of the decoy cohort (``n_fits``, ``n_capped_fits`` stopped at ``max_iter``,
      ``n_failed_fits`` stopped by the restart limit; :class:`mpaudit.gavish.GavishResult`);
    * ``decoy_mps``: cache tag -> the decoy cohort's MPs.
    """
    if arms is None:
        arms = {"all": tuple(strata)}
    elif strata:
        raise ValueError("pass either arms or strata, not both")
    if not arms:
        raise ValueError("arms is empty")
    if n_decoys < 1:
        raise ValueError("n_decoys must be >= 1")
    cache = cache or NMFCache()
    unit = sample_map(adata, keys, recurrence_key) if recurrence_key else None
    study = sample_map(adata, keys, keys.study) if keys.study else None
    refs = _reference_kappa(adata, keys, real_mps, counts_layer, min_cells, min_detect=MIN_DETECT, min_genes=MIN_GENES,
                            min_kappa=min_kappa, n_null_floor=n_null_floor, seed=seed)
    # samples entering kappa_retained per MP (decided on the real counts, so the same in every arm and decoy)
    n_enter = {m: sum(bool(per_set.get(m, {}).get("enter")) for _, _, per_set in refs) for m in real_mps}
    rows, per_decoy, decoy_mps = [], [], {}
    for arm, arm_cols in arms.items():
        cols = tuple(sorted(set(arm_cols)))
        best = {m: [] for m in real_mps}
        retained = {m: [] for m in real_mps}
        for d in range(n_decoys):
            kwargs = dict(counts_layer=counts_layer, strata=cols, n_depth_bins=n_depth_bins, method=method,
                          seed=seed + d)
            tag = decoy_tag(adata, keys, **kwargs)
            dec = make_decoy(adata, keys, **kwargs)
            for m, ratios in _retained_ratios(refs, get_counts(dec, counts_layer)).items():
                retained[m].append(float(np.median(ratios)))
            spg = compute_programs(dec, keys, params, counts_layer=counts_layer, cache=cache, data_tag=tag,
                                   n_jobs=n_jobs)
            res = assemble(spg, params, unit_of_sample=unit, study_of_sample=study)
            sets = res.mp_sets()
            decoy_mps[tag] = sets
            per_decoy.append({"arm": arm, "decoy": tag, "n_robust": len(res.robust), "n_mps": len(sets),
                              "n_fits": res.n_fits, "n_capped_fits": res.n_capped_fits,
                              "n_failed_fits": res.n_failed_fits})
            J = jaccard_matrix(real_mps, sets) if sets else np.zeros((len(real_mps), 0))
            for i, m in enumerate(real_mps):
                best[m].append(float(J[i].max()) if J.shape[1] else 0.0)
            log(f"[decoy] {arm}: {tag}: {len(res.robust)} robust programs, {len(sets)} decoy MPs; "
                f"{res.n_capped_fits}/{res.n_fits} NMF fits stopped at max_iter, {res.n_failed_fits} by the restart "
                "limit")
        for m in real_mps:
            b = np.array(best[m])
            hits = int((b >= min_jaccard).sum())
            rows.append({"mp": m, "arm": arm, "n_decoys": n_decoys, "hits": hits, "recruitment": hits / n_decoys,
                         "mean_best_jaccard": float(b.mean()), "max_best_jaccard": float(b.max()),
                         "kappa_retained": float(np.mean(retained[m])) if retained[m] else np.nan,
                         "kappa_retained_n_samples": n_enter[m]})
    return pd.DataFrame(rows, columns=PER_MP_COLUMNS), pd.DataFrame(per_decoy, columns=PER_DECOY_COLUMNS), decoy_mps
