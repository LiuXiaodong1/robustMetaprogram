"""Python port of the Gavish et al. 2023 / 3CA metaprogram (MP) pipeline.

Steps (3CA defaults in brackets; every constant is a field of :class:`GavishParams`, which checks
its fields when constructed: unknown option names and values raise ``ValueError``, and numbers
given as strings, as PyYAML reads ``1e-8``, are converted):

1. Per sample, malignant cells with >= ``min_genes_per_cell`` detected genes [1,000]:
   log2(CPM/10 + 1) (inferred: the R code log-normalises upstream of ``Process_mat``, inside the
   private ``Return_All_cells_Separated``; see below); then, as in R ``Process_mat``
   (``MP_analysis_new_batch.R:95-106``): keep the top genes by mean expression [7,000]; centre each
   gene; set negatives to zero; drop genes and cells left all zero (NMF cannot fit them).
2. sNMF/R at each rank [4-9], keeping the best of ``nrun`` restarts [10]; each program is
   summarised by its top genes [50]. All-zero basis columns (dead components) are skipped, so
   they never enter the program filters (counted in ``GavishResult.n_dead_programs``). Runs
   stopped by the iteration cap or by the restart limit are counted (``n_capped_fits``,
   ``n_failed_fits``); :func:`mpaudit.pipeline.run_p1` reports the counts per run (below).
3. Robust programs: overlap >= ``intra_min`` [35] with another program of the same sample;
   max overlap >= ``inter_min`` [10] with a program from another sample; greedy within-sample
   redundancy pruning at ``intra_max`` [10].
4. Greedy MP clustering: founder = program with the most partners at >= 10 shared genes
   (needs > ``min_group_size`` [5] partners); programs join while they share >= 10 genes with
   the evolving 50-gene MP list; MP genes = top 50 by frequency, ties broken by NMF weight.
5. Optional 3CA v2 steps (off by default, i.e. the Gavish 2023 core algorithm): MP filters by
   number of studies and programs (:func:`post_filter`), then the v2 final gene order
   (``v2_gene_order``; :func:`v2_gene_order`).

Checked against tiroshlab/3ca (``ITH_hallmarks/Generating_MPs/*.R``, ``3CA_v2/MP_analysis_new_batch.R``)
and renozao/NMF (``R/algorithms-snmf.R``). Remaining differences:

* NMF optimiser and stopping rule (HALS on the same objective; see :mod:`mpaudit.nmf`). The port
  stops on a relative objective change tested every 5 iterations, capped at ``max_iter`` [5,000];
  R tests assignment stability plus a KKT residual every 5 iterations, capped at 20,000. Fits
  that stop at the cap are counted in ``GavishResult.n_capped_fits``, and fits stopped by the
  restart limit (all-zero rows of ``H``; R also warns and stops) in ``n_failed_fits`` (per sample
  and rank in ``SamplePrograms.fits``). The P1 run reports them for the reference run, each decoy
  cohort (``decoy_runs.csv``) and each instability perturbation (``instability_runs.csv``), with
  totals in ``run_config.json`` (``nmf_fits``) and one line in ``summary.md``, and warns when a fit
  was stopped by the restart limit. On the toy cohort the tolerance binds long before the
  cap, and at ``tol=1e-5`` the port stopped about 5x earlier than R's rule would, with top-50
  lists not yet stable; the default is therefore ``tol=1e-7`` (numbers in :mod:`mpaudit.nmf`).
  Parity runs against R use ``tol=1e-8`` (``gavish.tol: 1.0e-8``; see ``scripts/parity``).
* Sample order, which sets the column order of the program matrices and therefore how ties are
  broken in the robust-program filter, the founder ranking and the dendrogram. 3CA v2 iterates
  over NMF output files: the old batch, then the new batch, each by study folder and then file
  name (``list.files`` order, which follows R's collation locale). The file names carry a suffix
  (``<sample>_rank4_9_nrun10.RDS``), so where one sample name is a prefix of another, R's order
  can differ from a sort of the bare names: in the C locale ``P10_rank...`` sorts before
  ``P1_rank...`` (digits and capitals sort before ``_``), so R lists MGH100, MGH10, P10, P1, P2
  where the bare names sort MGH10, MGH100, P1, P10, P2. ``sample_order`` offers ``"alphabetical"``
  (default; sample name by code point), ``"study"`` (study, then sample name, by code point; not
  R's order in general, for the reason above) and ``"input"``. To reproduce v2 exactly, list the
  samples in R's file order (e.g. from the NMF output folders) and use ``"input"``, passing that
  list as ``samples`` to :func:`compute_programs` or ordering the cohort by it. Ties are common:
  overlaps are integers and most program pairs sit at the same distance (50).
* Dead components: R keeps an all-zero basis column; its "top 50 genes" are then the first 50
  universe genes (R's sort is stable), i.e. the sample's highest-mean genes, the same at every
  rank. The port skips such columns (step 2). Both codes restart when a row of ``H`` is all zero, so dead
  columns survive only when the last W update zeroes a column or the restart limit is hit.
* MP gene order: the Gavish 2023 code lists MP genes in clustering order (frequency among member
  programs); 3CA v2 re-sorts each MP by study-deduplicated frequency. Membership is the same,
  but top-k slices differ; set ``v2_gene_order=True`` for the v2 order.
* Gene ties in R ``table()`` follow the collation locale; the port sorts by code point (same as
  R in the C locale; can differ for symbols with punctuation or mixed case).
* Samples are matched exactly (the R code uses ``grep`` on names, which misfires when one
  sample name is a substring of another).
* Founder ties follow an average-linkage dendrogram reordered by column means, as in the R code,
  but scipy's linkage can break exact distance ties differently from R's ``hclust``.
* The cell filter (>= 1,000 detected genes) and the transform log2(CPM/10 + 1), with CPM over
  each cell's total across all genes, are inferred. The function the MP script calls
  (``Return_All_cells_Separated(..., log_normalize = TRUE)``, ``MP_analysis_new_batch.R:111``)
  does cell QC and log-normalisation and is not public. The filter comes from other 3CA scripts;
  the transform from the commented-out ``MP_analysis_new_batch.R:99`` and from
  ``MPs_distribution/MP_distribution.R:47`` (a scoring script); the CPM denominator is an
  assumption.
* The R code upper-cases gene symbols; this port keeps ``var_names`` as given (only matters when
  comparing against gene lists with different capitalisation).
* The 3CA v2 manual removal of 11 MPs cannot be reproduced; label those MPs instead
  (``typer.labels`` in the P1 config).
"""
from __future__ import annotations

import dataclasses
import hashlib
import inspect
import math
import pickle
import zlib
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .data import Keys, get_counts, log2_cpm10, malignant_mask, sample_map
from .nmf import best_of_runs


SAMPLE_ORDERS = ("alphabetical", "study", "input")
FOUNDER_ORDERS = ("hclust", "input")
GENE_UNIVERSES = ("top_mean", "hvg")
_TRUE, _FALSE = {"true", "yes", "on", "1"}, {"false", "no", "off", "0"}


def _as_int(name: str, v) -> int:
    if isinstance(v, (bool, np.bool_)):
        raise ValueError(f"GavishParams.{name} must be an integer, got {v!r}")
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, str):
        try:
            return int(v.strip())
        except ValueError:
            pass
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"GavishParams.{name} must be an integer, got {v!r}") from None
    if not (math.isfinite(f) and f.is_integer()):
        raise ValueError(f"GavishParams.{name} must be an integer, got {v!r}")
    return int(f)


def _as_float(name: str, v) -> float:
    if isinstance(v, (bool, np.bool_)):
        raise ValueError(f"GavishParams.{name} must be a number, got {v!r}")
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"GavishParams.{name} must be a number, got {v!r}") from None
    if not math.isfinite(f):
        raise ValueError(f"GavishParams.{name} must be finite, got {v!r}")
    return f


def _as_bool(name: str, v) -> bool:
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)) and v in (0, 1):
        return bool(v)
    if isinstance(v, str) and v.strip().lower() in _TRUE | _FALSE:
        return v.strip().lower() in _TRUE
    raise ValueError(f"GavishParams.{name} must be true or false, got {v!r}")


def _as_choice(choices: tuple[str, ...]):
    def conv(name: str, v) -> str:
        if not isinstance(v, str) or v not in choices:
            raise ValueError(f"Unknown {name} {v!r}; use one of {', '.join(map(repr, choices))}")
        return v
    return conv


def _as_ranks(name: str, v) -> tuple[int, ...]:
    if isinstance(v, (int, np.integer, str)) and not isinstance(v, (bool, np.bool_)):
        v = [v]
    try:
        ranks = tuple(_as_int(name, r) for r in v)
    except TypeError:
        raise ValueError(f"GavishParams.{name} must be a list of integers, got {v!r}") from None
    if not ranks or min(ranks) < 1 or len(set(ranks)) != len(ranks):
        raise ValueError(f"GavishParams.{name} must be distinct integers >= 1, got {v!r}")
    return ranks


def _as_optional_float(name: str, v) -> float | None:
    return None if v is None else _as_float(name, v)


# Converter per field, chosen by the field's annotation (or by name for the enumerated options).
_CONVERTERS = {"int": _as_int, "float": _as_float, "bool": _as_bool, "float | None": _as_optional_float,
               "tuple[int, ...]": _as_ranks}
_CHOICES = {"gene_universe": GENE_UNIVERSES, "sample_order": SAMPLE_ORDERS, "founder_order": FOUNDER_ORDERS}
_AT_LEAST_ONE = ("nrun", "n_top", "n_universe", "max_iter")
_NON_NEGATIVE = ("beta", "tol")


@dataclass(frozen=True)
class GavishParams:
    """Settings of the port (3CA defaults unless noted).

    Every field is checked when the object is built (directly, by :meth:`replace` or by
    :meth:`from_dict`): numbers given as strings (PyYAML reads ``tol: 1e-8``, without a decimal
    point, as a string) and integral floats are converted to the field's type, ``ranks`` becomes a
    tuple of distinct ints, and an unknown field name, an unknown ``gene_universe``,
    ``sample_order`` or ``founder_order``, or a value that cannot be converted raises
    ``ValueError`` before any NMF runs.
    """

    ranks: tuple[int, ...] = (4, 5, 6, 7, 8, 9)
    nrun: int = 10
    n_top: int = 50
    gene_universe: str = "top_mean"  # "top_mean" (3CA) or "hvg" (top variance of log expression)
    n_universe: int = 7000
    center: bool = True
    min_cells: int = 50  # minimum malignant cells per sample (an assumption; not verified in 3CA code)
    min_genes_per_cell: int = 1000  # low-complexity cell filter used across 3CA scripts; 0 disables
    beta: float = 0.01
    eta: float | None = None
    max_iter: int = 5000  # R's default is 20,000 with a different stopping rule (see mpaudit.nmf)
    tol: float = 1e-7  # 1e-5 left top-50 lists unconverged (median 44/50 vs tol 1e-8); see mpaudit.nmf
    seed: int = 0
    intra_min: int = 35
    intra_max: int = 10
    inter_min: int = 10
    inter_filter: bool = True
    min_intersect_initial: int = 10
    min_intersect_cluster: int = 10
    min_group_size: int = 5
    # Column order of the programs (R iterates over NMF output files): "alphabetical" (sample
    # name), "study" (study, then sample name) or "input" (order of first appearance in the cohort,
    # or the ``samples`` passed to compute_programs). Only "input", with the samples in R's file
    # order, reproduces 3CA v2 (file names carry a suffix; see the module docstring).
    sample_order: str = "alphabetical"
    founder_order: str = "hclust"       # R ranks founders in dendrogram order; or "input"
    # 3CA v2-style steps (off by default, i.e. Gavish 2023 core algorithm)
    min_programs_per_mp: int = 0
    min_studies_per_mp: int = 1
    drop_two_study_singleton: bool = True  # v2: with 2 studies, drop if one study gives a single NMF
    v2_gene_order: bool = False  # v2: re-sort MP genes by study-deduplicated frequency (needs studies)

    def __post_init__(self):
        for f in dataclasses.fields(self):
            conv = _as_choice(_CHOICES[f.name]) if f.name in _CHOICES else _CONVERTERS[f.type]
            object.__setattr__(self, f.name, conv(f.name, getattr(self, f.name)))
        for name in _AT_LEAST_ONE:
            if getattr(self, name) < 1:
                raise ValueError(f"GavishParams.{name} must be >= 1, got {getattr(self, name)}")
        for name in _NON_NEGATIVE:
            if getattr(self, name) < 0:
                raise ValueError(f"GavishParams.{name} must be >= 0, got {getattr(self, name)}")

    def replace(self, **kw) -> "GavishParams":
        unknown = sorted(set(kw) - {f.name for f in dataclasses.fields(self)})
        if unknown:
            raise ValueError(f"Unknown GavishParams field(s) {unknown}; known: "
                             f"{', '.join(f.name for f in dataclasses.fields(self))}")
        return dataclasses.replace(self, **kw)

    @classmethod
    def from_dict(cls, d: dict | None) -> "GavishParams":
        return cls().replace(**dict(d or {}))


@dataclass
class SamplePrograms:
    """NMF basis matrices for one sample, over that sample's gene universe.

    ``fits`` holds per-rank fit diagnostics over the ``nrun`` restarts: ``n_runs``, ``n_capped``
    (runs stopped by ``max_iter`` before the tolerance test passed), ``n_failed`` (runs stopped by
    the restart limit of :func:`mpaudit.nmf.snmf_r`, i.e. unusable fits) and ``mean_iter``.
    """

    sample: str
    universe: np.ndarray
    W: dict[int, np.ndarray] = field(default_factory=dict)  # rank -> (n_universe, rank) float32
    fits: dict[int, dict] = field(default_factory=dict)     # rank -> {"n_runs", "n_capped", "n_failed", "mean_iter"}
    _index: dict | None = field(default=None, repr=False)

    def gene_index(self) -> dict[str, int]:
        if self._index is None:
            self._index = {g: i for i, g in enumerate(self.universe)}
        return self._index


@dataclass
class MP:
    name: str
    genes: list[str]
    members: list[str]


@dataclass
class GavishResult:
    params: GavishParams
    top: dict[str, list[str]]  # program id -> top genes
    program_sample: dict[str, str]
    robust: list[str]
    mps: list[MP]
    n_dead_programs: int = 0
    n_fits: int = 0          # NMF runs behind the programs (samples x ranks x nrun)
    n_capped_fits: int = 0   # of which stopped by max_iter
    n_failed_fits: int = 0   # of which stopped by the restart limit (all-zero rows of H; unusable)

    def mp_sets(self) -> dict[str, list[str]]:
        return {m.name: m.genes for m in self.mps}

    def programs_table(self) -> pd.DataFrame:
        robust = set(self.robust)
        member_of = {p: m.name for m in self.mps for p in m.members}
        rows = [
            {"program": p, "sample": self.program_sample[p], "robust": p in robust,
             "mp": member_of.get(p, ""), "genes": ",".join(g)}
            for p, g in self.top.items()
        ]
        return pd.DataFrame(rows)


def program_id(sample: str, rank: int, j: int) -> str:
    return f"{sample}|r{rank}|{j + 1}"


# --------------------------------------------------------------------------- per-sample NMF

def _universe_index(L: sp.csr_matrix, params: GavishParams) -> np.ndarray:
    """Columns of ``L`` (cells x genes, log expression) forming the sample's gene universe.

    The top ``n_universe`` genes by mean (3CA) or variance of log expression, in decreasing score
    order, minus the genes whose row ``Process_mat`` would drop as all zero: after centring and
    zeroing negatives that is exactly the genes constant across cells (without centring, the genes
    not detected in any cell). Testing constancy on ``L`` is exact, unlike testing the centred row,
    where the rounding of a constant gene's mean can leave tiny positive values.
    """
    mean = np.asarray(L.mean(axis=0)).ravel()
    if params.gene_universe == "top_mean":
        score = mean
    elif params.gene_universe == "hvg":
        sq = np.asarray(L.multiply(L).mean(axis=0)).ravel()
        score = sq - mean**2
    else:
        raise ValueError(f"Unknown gene_universe {params.gene_universe!r}")
    keep = np.argsort(-score, kind="stable")[: params.n_universe]
    top = L[:, keep]
    hi = top.max(axis=0).toarray().ravel()
    if params.center:
        return keep[hi > top.min(axis=0).toarray().ravel()]
    return keep[hi > 0]


def prepare_matrix(counts: sp.spmatrix, genes: np.ndarray, params: GavishParams) -> tuple[np.ndarray, np.ndarray]:
    """3CA per-sample matrix: genes x cells, log2(CPM/10+1) (inferred; done upstream of R
    ``Process_mat``, see the module docstring), then as in R ``Process_mat``: top genes, centred,
    negatives zeroed, all-zero genes and cells dropped. Returns ``(A, universe)``, ``universe``
    naming the rows of ``A``."""
    L = log2_cpm10(counts)
    keep = _universe_index(L, params)
    A = L[:, keep].toarray().T
    if params.center:
        A = A - A.mean(axis=1, keepdims=True)
        A[A < 0] = 0.0
    A = A[:, A.any(axis=0)]  # cells left all zero (at or below every gene's mean when centred)
    return A, np.asarray(genes)[keep]


def _rank_seeds(seed: int, sample: str, rank: int, nrun: int) -> list[np.random.SeedSequence]:
    base = np.random.SeedSequence([seed, zlib.crc32(sample.encode()), rank])
    return base.spawn(nrun)


def _fit_summary(infos: list[dict], max_iter: int) -> dict:
    """Diagnostics of the ``nrun`` restarts at one rank (see :class:`SamplePrograms`)."""
    return {"n_runs": len(infos),
            "n_capped": sum(not i["converged"] and i["n_iter"] >= max_iter for i in infos),
            "n_failed": sum(bool(i.get("failed", False)) for i in infos),
            "mean_iter": float(np.mean([i["n_iter"] for i in infos])) if infos else np.nan}


def _fit_job(args):
    sample, counts, genes, params, ranks, blas_threads = args
    if blas_threads:
        try:
            from threadpoolctl import threadpool_limits
            with threadpool_limits(blas_threads):
                return _fit_job((sample, counts, genes, params, ranks, None))
        except ImportError:
            pass
    A, universe = prepare_matrix(counts, genes, params)
    out, fits = {}, {}
    for r in ranks:
        W, infos = best_of_runs(A, r, nrun=params.nrun, seeds=_rank_seeds(params.seed, sample, r, params.nrun),
                                beta=params.beta, eta=params.eta, max_iter=params.max_iter, tol=params.tol,
                                return_info=True)
        out[r] = W.astype(np.float32)
        fits[r] = _fit_summary(infos, params.max_iter)
    return sample, universe, out, fits


# Layout of a cached value, the (universe, W, fit summary) tuple NMFCache stores; bump it whenever
# that layout or the meaning of its fields changes without a change to the sources code_version()
# hashes. Part of code_version().
CACHE_FORMAT = 2

# Hash of the code that determines a cached fit; part of every NMFCache key. Computed on first use
# (the package version is defined after this module is imported); tests may monkeypatch it.
CODE_VERSION: str | None = None
# Hash of the code that generates decoy counts, added to the key of fits on decoy data (data tags
# starting with DECOY_TAG, as mpaudit.decoy builds them). Computed on first use, like CODE_VERSION.
DECOY_CODE_VERSION: str | None = None
DECOY_TAG = "decoy-"
# The functions that turn a cohort into decoy counts: make_decoy and everything of mpaudit it calls.
DECOY_SOURCES = ("decoy.make_decoy", "decoy._is_continuous", "decoy._covariates", "decoy._strata",
                 "decoy._permute_within", "decoy._scrp_block", "decoy._glm_mean", "decoy._nb_dispersion",
                 "data.get_counts", "data.malignant_mask")


def code_version(decoy_data: bool = False) -> str:
    """Hash of the code behind a cached fit.

    For fits on real data: the package version, :data:`CACHE_FORMAT` and the sources that shape a
    cached fit, i.e. :mod:`mpaudit.nmf`, :func:`prepare_matrix` (with :func:`_universe_index`),
    :func:`_fit_job`, :func:`_fit_summary` (the cached fit diagnostics), :func:`_rank_seeds` and
    :func:`mpaudit.data.log2_cpm10`. With ``decoy_data``, also the functions that generate the
    decoy counts (:data:`DECOY_SOURCES`: :func:`mpaudit.decoy.make_decoy` and every helper of the
    package it calls; their docstrings included). Editing any of them invalidates an on-disk cache
    instead of silently serving old fits or old diagnostics. Editing anything else in
    :mod:`mpaudit.decoy` (its module docstring, :func:`~mpaudit.decoy.decoy_recruitment`,
    ``kappa_retained``) leaves every key unchanged, and editing the decoy generator leaves the keys
    of fits on real data unchanged.
    """
    global CODE_VERSION, DECOY_CODE_VERSION
    if CODE_VERSION is None:
        from . import __version__, data, nmf
        sources = [__version__, f"cache-format-{CACHE_FORMAT}", inspect.getsource(nmf),
                   inspect.getsource(prepare_matrix), inspect.getsource(_universe_index),
                   inspect.getsource(_fit_job), inspect.getsource(_fit_summary), inspect.getsource(_rank_seeds),
                   inspect.getsource(data.log2_cpm10)]
        CODE_VERSION = hashlib.sha1("\x1e".join(sources).encode()).hexdigest()[:16]
    if not decoy_data:
        return CODE_VERSION
    if DECOY_CODE_VERSION is None:
        from . import data, decoy
        modules = {"data": data, "decoy": decoy}
        sources = [inspect.getsource(getattr(modules[m], f)) for m, f in (s.split(".") for s in DECOY_SOURCES)]
        DECOY_CODE_VERSION = hashlib.sha1("\x1e".join(sources).encode()).hexdigest()[:16]
    return hashlib.sha1(f"{CODE_VERSION}\x1e{DECOY_CODE_VERSION}".encode()).hexdigest()[:16]


class NMFCache:
    """Caches per-(sample, rank) NMF bases and fit diagnostics, in memory or (when ``directory``
    is set) on disk.

    The key covers everything that changes the fit: the code version (:func:`code_version`; for
    data tags starting with :data:`DECOY_TAG` it includes the decoy generator), a fingerprint of
    the input (cell barcodes, gene names, counts layer, total counts), a data tag (real or decoy),
    the gene universe, centring, solver settings and seed. Threshold, clustering
    and leave-study-out perturbations therefore reuse fits; a cache directory cannot serve fits
    from a different dataset or from older code. With a directory, fits are not also kept in
    memory (atlas-scale runs would not fit).
    """

    def __init__(self, directory: str | Path | None = None):
        self.mem: dict[str, tuple[np.ndarray, np.ndarray, dict]] = {}
        self.dir = Path(directory) if directory else None
        if self.dir:
            self.dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(data_tag: str, fingerprint: str, sample: str, rank: int, p: GavishParams) -> str:
        parts = (code_version(decoy_data=data_tag.startswith(DECOY_TAG)), data_tag, fingerprint, sample, rank,
                 p.gene_universe, p.n_universe, p.center, p.nrun, p.seed, p.beta, p.eta, p.max_iter, p.tol)
        return hashlib.sha1(repr(parts).encode()).hexdigest()

    def get(self, key: str) -> tuple[np.ndarray, np.ndarray, dict] | None:
        """``(universe, W, fit_summary)`` or None."""
        if key in self.mem:
            return self.mem[key]
        if self.dir and (self.dir / f"{key}.pkl").exists():
            with open(self.dir / f"{key}.pkl", "rb") as fh:
                return pickle.load(fh)
        return None

    def put(self, key: str, universe: np.ndarray, W: np.ndarray, fit: dict) -> None:
        if not self.dir:
            self.mem[key] = (universe, W, fit)
        else:
            with open(self.dir / f"{key}.pkl", "wb") as fh:
                pickle.dump((universe, W, fit), fh, protocol=pickle.HIGHEST_PROTOCOL)


def _sample_cells(
    adata,
    keys: Keys,
    params: GavishParams,
    counts: sp.csr_matrix,
    *,
    samples: list[str] | None = None,
    cell_fraction: float = 1.0,
    subset_seed: int = 0,
) -> dict[str, np.ndarray]:
    """Rows of the cells each sample contributes to NMF: malignant cells with >=
    ``min_genes_per_cell`` detected genes (optionally subsampled), for the samples left with >=
    ``min_cells`` of them; in ``samples`` order, else in order of first appearance."""
    mal = malignant_mask(adata, keys)
    obs_sample = adata.obs[keys.sample].astype(str).to_numpy()
    if params.min_genes_per_cell > 0:
        detected = np.diff((counts > 0).tocsr().indptr)
        mal = mal & (detected >= params.min_genes_per_cell)
    order = samples if samples is not None else list(pd.unique(obs_sample[mal]))
    cells = {}
    for s in order:
        idx = np.flatnonzero(mal & (obs_sample == s))
        if cell_fraction < 1 and idx.size:
            rng = np.random.default_rng([subset_seed, zlib.crc32(s.encode())])
            idx = np.sort(rng.choice(idx, size=max(1, int(round(cell_fraction * idx.size))), replace=False))
        if idx.size >= params.min_cells:
            cells[s] = idx
    return cells


def check_keys(keys: Keys, params: GavishParams, adata=None, *, recurrence_key: str | None = None) -> None:
    """Raise ``ValueError`` for settings that :func:`assemble` would reject only after every NMF
    fit has run.

    ``sample_order="study"``, ``v2_gene_order`` and ``min_studies_per_mp > 1`` need a study per
    sample (``keys.study``). With ``adata``, that study column (when needed) and
    ``recurrence_key`` must be columns of ``adata.obs``. :func:`run_gavish` and
    :func:`compute_programs` call this before fitting; callers that build the study map
    themselves (``study_of_sample``) should call it too.
    """
    need = [opt for opt, on in (("sample_order='study'", params.sample_order == "study"),
                                ("v2_gene_order=True", params.v2_gene_order),
                                (f"min_studies_per_mp={params.min_studies_per_mp}", params.min_studies_per_mp > 1))
            if on]
    if need and not keys.study:
        raise ValueError(f"GavishParams {', '.join(need)} need(s) a study for each sample: set keys.study")
    if adata is not None:
        cols = {"keys.study": keys.study if need else None, "recurrence_key": recurrence_key}
        missing = [f"{k}={v!r}" for k, v in cols.items() if v and v not in adata.obs.columns]
        if missing:
            raise ValueError(f"Not a column of adata.obs: {', '.join(missing)}")


def eligible_universe(adata, keys: Keys, params: GavishParams, counts_layer: str | None = None) -> list[str]:
    """Genes that can enter an MP, sorted: the union over samples of the gene universes
    :func:`prepare_matrix` would use, with the cell and sample selection of
    :func:`compute_programs`. No NMF is run.

    This is the background for enrichment tests of MP genes (e.g. chromosome arms): a gene outside
    every sample's universe can never be an MP gene, so testing against all ``var_names`` is biased.
    """
    counts = get_counts(adata, counts_layer)
    genes = np.asarray(adata.var_names.astype(str))
    eligible = np.zeros(genes.size, dtype=bool)
    for idx in _sample_cells(adata, keys, params, counts).values():
        eligible[_universe_index(log2_cpm10(counts[idx]), params)] = True
    return sorted(genes[eligible].tolist())


def compute_programs(
    adata,
    keys: Keys,
    params: GavishParams,
    *,
    counts_layer: str | None = None,
    samples: list[str] | None = None,
    cache: NMFCache | None = None,
    data_tag: str = "data",
    cell_fraction: float = 1.0,
    subset_seed: int = 0,
    n_jobs: int = 1,
) -> dict[str, SamplePrograms]:
    """Run (or fetch from cache) per-sample NMF for every sample with enough malignant cells.

    Samples come back in ``samples`` order, else in order of first appearance in ``adata``
    (the order ``sample_order="input"`` uses). Settings the later MP stage needs are checked
    first (:func:`check_keys`), so a configuration error fails before any NMF.
    """
    check_keys(keys, params, adata)
    cache = cache or NMFCache()
    counts = get_counts(adata, counts_layer)
    genes = np.asarray(adata.var_names.astype(str))
    obs_names = np.asarray(adata.obs_names.astype(str))
    gene_hash = hashlib.sha1("\x1f".join(genes).encode()).hexdigest()[:16]
    cells = _sample_cells(adata, keys, params, counts, samples=samples, cell_fraction=cell_fraction,
                          subset_seed=subset_seed)

    result: dict[str, SamplePrograms] = {}
    fingerprints: dict[str, str] = {}
    jobs = []
    for s, idx in cells.items():
        fp = hashlib.sha1(("\x1f".join(obs_names[idx]) + gene_hash + str(counts_layer)
                           + f"{counts[idx].sum():.10g}").encode()).hexdigest()[:20]
        fingerprints[s] = fp
        sp_obj = SamplePrograms(sample=s, universe=np.array([], dtype=object))
        missing = []
        for r in params.ranks:
            hit = cache.get(NMFCache.key(data_tag, fp, s, r, params))
            if hit is None:
                missing.append(r)
            else:
                sp_obj.universe, sp_obj.W[r], sp_obj.fits[r] = hit
        result[s] = sp_obj
        if missing:
            jobs.append((s, counts[idx], genes, params, tuple(missing), 1 if n_jobs > 1 else None))

    if jobs:
        if n_jobs > 1:
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                outputs = list(ex.map(_fit_job, jobs))
        else:
            outputs = [_fit_job(j) for j in jobs]
        for s, universe, Ws, fits in outputs:
            for r, W in Ws.items():
                cache.put(NMFCache.key(data_tag, fingerprints[s], s, r, params), universe, W, fits[r])
                result[s].W[r] = W
                result[s].fits[r] = fits[r]
            result[s].universe = universe
            result[s]._index = None
    return result


# --------------------------------------------------------------------------- program filters

def top_genes(weights: np.ndarray, universe: np.ndarray, n: int) -> list[str]:
    order = np.argsort(-weights, kind="stable")[:n]
    return [str(g) for g in universe[order]]


def order_samples(samples: list[str], params: GavishParams, study_of_sample: dict[str, str] | None = None) -> list[str]:
    """Samples in ``params.sample_order``: ``"alphabetical"``, ``"study"`` (study, then sample
    name, both by code point) or ``"input"`` (as given). Neither sort is R's file order in general
    (see the module docstring); to reproduce 3CA v2, give the samples in R's order with
    ``"input"``."""
    if params.sample_order == "alphabetical":
        return sorted(samples)
    if params.sample_order == "study":
        if study_of_sample is None:
            raise ValueError("sample_order='study' needs a study for each sample (set keys.study)")
        return sorted(samples, key=lambda s: (study_of_sample.get(s, ""), s))
    if params.sample_order == "input":
        return list(samples)
    raise ValueError(f"Unknown sample_order {params.sample_order!r}; use 'alphabetical', 'study' or 'input'")


def collect_top(sample_programs: dict[str, SamplePrograms], params: GavishParams,
                study_of_sample: dict[str, str] | None = None) -> tuple[dict, dict, int]:
    """Top genes per program, in R's column order (sample, then rank, then program index).

    All-zero basis columns (dead components) are skipped: their "top genes" would just be the
    sample's highest-mean genes, identical across ranks and samples.
    """
    top, program_sample, dead = {}, {}, 0
    for s in order_samples(list(sample_programs), params, study_of_sample):
        spg = sample_programs[s]
        for r in params.ranks:
            W = spg.W[r]
            for j in range(W.shape[1]):
                if not np.any(W[:, j] > 0):
                    dead += 1
                    continue
                pid = program_id(s, r, j)
                top[pid] = top_genes(W[:, j], spg.universe, params.n_top)
                program_sample[pid] = s
    return top, program_sample, dead


def _overlap(sets_a: list[set], sets_b: list[set]) -> np.ndarray:
    return np.array([[len(a & b) for b in sets_b] for a in sets_a], dtype=int).reshape(len(sets_a), len(sets_b))


def robust_programs(
    top: dict[str, list[str]],
    program_sample: dict[str, str],
    params: GavishParams,
    unit_of_sample: dict[str, str] | None = None,
) -> list[str]:
    """Port of 3CA ``robust_nmf_programs``. ``unit_of_sample`` defines what counts as
    "another sample" for recurrence (default: each sample is its own unit)."""
    unit_of_sample = unit_of_sample or {}
    by_sample: dict[str, list[str]] = {}
    for pid, s in program_sample.items():
        if pid in top:
            by_sample.setdefault(s, []).append(pid)

    # within-sample reproducibility across ranks
    selected: dict[str, list[str]] = {}
    for s, pids in by_sample.items():
        sets = [set(top[p]) for p in pids]
        ov = _overlap(sets, sets)
        np.fill_diagonal(ov, -1)
        best = ov.max(axis=1) if len(pids) > 1 else np.full(len(pids), -1)
        selected[s] = [p for p, b in zip(pids, best) if b >= params.intra_min]

    all_sel = [p for s in selected for p in selected[s]]
    if not all_sel:
        return []
    pos = {p: i for i, p in enumerate(all_sel)}
    sets = [set(top[p]) for p in all_sel]
    inter = _overlap(sets, sets)
    unit = np.array([unit_of_sample.get(program_sample[p], program_sample[p]) for p in all_sel])

    final: list[str] = []
    for s, pids in selected.items():
        if not pids:
            continue
        cols = np.array([pos[p] for p in pids])
        my_unit = unit_of_sample.get(s, s)
        others = np.flatnonzero(unit != my_unit)
        maxes = inter[np.ix_(others, cols)].max(axis=0) if others.size else np.zeros(len(cols), int)
        order = np.argsort(-maxes, kind="stable")
        ranked = [(pids[i], maxes[i]) for i in order]
        if params.inter_filter:
            ranked = [(p, m) for p, m in ranked if m >= params.inter_min]
        chosen: list[str] = []
        for p, _ in ranked:
            if not chosen or max(inter[pos[c], pos[p]] for c in chosen) <= params.intra_max:
                chosen.append(p)
        final.extend(chosen)
    return final


def dendrogram_order(ov: np.ndarray, n_top: int) -> np.ndarray:
    """Leaf order of R's ``reorder(as.dendrogram(hclust(as.dist(n_top - ov), "average")), colMeans(ov))``.

    At every node the child with the smaller summed weight goes first (ties keep the left child
    first). Exact distance ties may resolve differently from R's ``hclust``.
    """
    from scipy.cluster.hierarchy import linkage, to_tree
    from scipy.spatial.distance import squareform

    n = ov.shape[0]
    if n < 3:
        return np.arange(n)
    D = (n_top - ov).astype(float)
    np.fill_diagonal(D, 0.0)
    D = np.maximum(D, 0.0)
    root = to_tree(linkage(squareform(D, checks=False), method="average"))
    w = ov.mean(axis=0)
    leaves, weight = {}, {}
    stack = [(root, False)]
    while stack:  # iterative post-order (deep trees would exceed Python's recursion limit)
        node, done = stack.pop()
        if node.is_leaf():
            leaves[node.id], weight[node.id] = [node.id], float(w[node.id])
        elif not done:
            stack += [(node, True), (node.get_right(), False), (node.get_left(), False)]
        else:
            a, b = node.get_left().id, node.get_right().id
            if weight[b] < weight[a]:
                a, b = b, a
            leaves[node.id] = leaves.pop(a) + leaves.pop(b)
            weight[node.id] = weight.pop(a) + weight.pop(b)
    return np.array(leaves[root.id])


def cluster_mps(
    top: dict[str, list[str]],
    weight_of,
    params: GavishParams,
) -> list[MP]:
    """Port of the 3CA greedy MP clustering.

    ``top`` maps the robust program ids, in R's column order (sample, rank, index), to their top
    genes; ``weight_of(pid, genes)`` returns those genes' NMF weights in the program's full basis
    (NaN if absent). Candidates are scanned in ``top`` order; founders are ranked in dendrogram
    order when ``params.founder_order == "hclust"``, as in the R code.
    """
    ids = list(top)
    n = len(ids)
    if n == 0:
        return []
    n_top = params.n_top
    sets = [set(top[i]) for i in ids]
    ov = _overlap(sets, sets)
    active = np.ones(n, dtype=bool)
    pool = list(range(n))
    seq = dendrogram_order(ov, n_top) if params.founder_order == "hclust" else np.arange(n)

    def ranked_founders():
        idx = seq[active[seq]]
        if idx.size == 0:
            return idx, np.array([], dtype=int)
        cnt = (ov[np.ix_(idx, idx)] >= params.min_intersect_initial).sum(axis=0) - 1
        o = np.argsort(-cnt, kind="stable")
        return idx[o], cnt[o]

    def ranked_pool(genes_mp: list[str]):
        gs = set(genes_mp)
        inter = np.array([len(gs & sets[p]) for p in pool], dtype=int)
        o = np.argsort(-inter, kind="stable")
        return [pool[i] for i in o], inter[o]

    mps: list[MP] = []
    order, cnt = ranked_founders()
    while cnt.size and cnt[0] > params.min_group_size:
        founder = int(order[0])
        cluster = [founder]
        genes_mp = list(top[ids[founder]])
        history = list(genes_mp)
        pool.remove(founder)
        cand, inter = ranked_pool(genes_mp)
        while cand and inter[0] >= params.min_intersect_cluster:
            p = cand[0]
            cluster.append(p)
            counts = Counter(history + top[ids[p]])
            names = sorted(counts)  # R table(): alphabetical
            names.sort(key=lambda g: -counts[g])  # stable sort by decreasing frequency
            c_border = counts[names[n_top - 1]]
            border = [g for g in names if counts[g] == c_border]
            if len(border) > 1:
                scored = []
                for member in cluster:
                    w = np.asarray(weight_of(ids[member], border), dtype=float)
                    pairs = [(s, g) for g, s in zip(border, w) if not np.isnan(s)]
                    pairs.sort(key=lambda t: -t[0])
                    scored.extend(pairs)
                scored.sort(key=lambda t: -t[0])
                by_score = list(dict.fromkeys(g for _, g in scored))
                genes_tmp = [g for g in names if counts[g] > c_border] + by_score
            else:
                genes_tmp = names[:n_top]
            history += top[ids[p]]
            genes_mp = genes_tmp[:n_top]
            pool.remove(p)
            cand, inter = ranked_pool(genes_mp)
        mps.append(MP(name=f"MP{len(mps) + 1}", genes=genes_mp, members=[ids[i] for i in cluster]))
        active[cluster] = False
        order, cnt = ranked_founders()
    return mps


def post_filter(mps: list[MP], program_sample: dict[str, str], study_of_sample: dict[str, str] | None,
                params: GavishParams) -> list[MP]:
    """3CA v2 filters (``MP_analysis_new_batch.R``): minimum studies per MP; with exactly two
    studies, drop the MP if one study contributes a single NMF; minimum NMFs per MP."""
    if params.min_studies_per_mp > 1 and study_of_sample is None:
        raise ValueError("min_studies_per_mp > 1 needs keys.study")
    out = []
    for m in mps:
        if params.min_studies_per_mp > 1:
            per_study = Counter(study_of_sample.get(program_sample[p], "") for p in m.members)
            if len(per_study) < params.min_studies_per_mp:
                continue
            if params.drop_two_study_singleton and len(per_study) == 2 and 1 in per_study.values():
                continue
        if len(m.members) < params.min_programs_per_mp:
            continue
        out.append(m)
    for i, m in enumerate(out):
        m.name = f"MP{i + 1}"
    return out


def v2_gene_order(mps: list[MP], top: dict[str, list[str]], program_sample: dict[str, str],
                  study_of_sample: dict[str, str] | None) -> list[MP]:
    """3CA v2 final gene order: each MP's genes re-sorted by the number of studies in which at
    least one member program has the gene among its top genes (decreasing; ties alphabetical, as
    R's ``table()`` followed by a stable decreasing sort). Membership is unchanged.

    v2 applies it after its MP filters (here: after :func:`post_filter`). The R code takes the
    study from the program name (dataset folder); here it comes from ``study_of_sample``. The
    ``MP_list_sorted_above_1``/``_above_2`` lists R also builds are not used anywhere in the
    script and are not ported.
    """
    if study_of_sample is None:
        raise ValueError("v2_gene_order needs a study for each sample (set keys.study)")
    # MP_analysis_new_batch.R:603-623, per MP I (N1 = study of each NMF in Cluster_list[[I]];
    # nmf_programs_sig_ccle_KEEP = top genes of the robust programs):
    #   N2 <- unique(N1)
    #   for (J in 1:length(N2)) { ## run across studies, and collect unique gene from the NMF programs in each
    #     curr_NMFs <- Cluster_list[[I]][which(N1 == N2[J])]
    #     R         <- match(curr_NMFs , colnames(nmf_programs_sig_ccle_KEEP))
    #     gene_vector <- c(gene_vector , unique(as.vector(nmf_programs_sig_ccle_KEEP[,R])))
    #   }
    #   sorted_gene_vector <- sort(table(gene_vector) , decreasing = T)
    #   R2                 <- match( names(sorted_gene_vector) , MP_list[[I]]  )
    #   MP_genes_sorted          <- names(sorted_gene_vector)[which(!is.na(R2))]
    out = []
    for m in mps:
        genes_by_study: dict[str, set[str]] = {}
        for p in m.members:
            genes_by_study.setdefault(study_of_sample.get(program_sample[p], ""), set()).update(top[p])
        n_studies = Counter(g for genes in genes_by_study.values() for g in genes)
        names = sorted(n_studies)  # R table(): alphabetical
        names.sort(key=lambda g: -n_studies[g])  # stable sort by decreasing study count
        in_mp = set(m.genes)
        out.append(MP(name=m.name, genes=[g for g in names if g in in_mp], members=list(m.members)))
    return out


def assemble(
    sample_programs: dict[str, SamplePrograms],
    params: GavishParams,
    *,
    unit_of_sample: dict[str, str] | None = None,
    study_of_sample: dict[str, str] | None = None,
    drop_samples: set[str] | frozenset = frozenset(),
) -> GavishResult:
    """Robust-program filter, MP clustering and optional v2 steps on precomputed NMF bases."""
    kept = {s: v for s, v in sample_programs.items() if s not in drop_samples}
    top, program_sample, dead = collect_top(kept, params, study_of_sample)
    robust = robust_programs(top, program_sample, params, unit_of_sample)
    robust_set = set(robust)
    robust_top = {p: g for p, g in top.items() if p in robust_set}  # R column order

    def weight_of(pid: str, genes: list[str]) -> np.ndarray:
        s, r, j = pid.rsplit("|", 2)
        spg = kept[s]
        col = spg.W[int(r[1:])][:, int(j) - 1]
        idx = spg.gene_index()
        return np.array([col[idx[g]] if g in idx else np.nan for g in genes], dtype=float)

    mps = cluster_mps(robust_top, weight_of, params)
    mps = post_filter(mps, program_sample, study_of_sample, params)
    if params.v2_gene_order:
        mps = v2_gene_order(mps, top, program_sample, study_of_sample)
    fits = [spg.fits[r] for spg in kept.values() for r in params.ranks if r in spg.fits]
    return GavishResult(params=params, top=top, program_sample=program_sample, robust=list(robust_top),
                        mps=mps, n_dead_programs=dead, n_fits=sum(f["n_runs"] for f in fits),
                        n_capped_fits=sum(f["n_capped"] for f in fits),
                        n_failed_fits=sum(f.get("n_failed", 0) for f in fits))


def run_gavish(
    adata,
    keys: Keys,
    params: GavishParams | None = None,
    *,
    counts_layer: str | None = None,
    cache: NMFCache | None = None,
    data_tag: str = "data",
    recurrence_key: str | None = None,
    n_jobs: int = 1,
) -> GavishResult:
    """Run the full 3CA-style MP pipeline on an AnnData cohort.

    ``recurrence_key`` is the obs column defining independent units for the cross-sample
    recurrence filter (default: sample, as in 3CA). Set it to the tumour column when a cohort
    contains several samples per tumour, so that a tumour cannot "recur" with itself.

    The keys and settings the MP stage needs are checked (:func:`check_keys`) and the recurrence
    and study maps built before any NMF runs, so a configuration error does not cost the fits.
    """
    params = params or GavishParams()
    check_keys(keys, params, adata, recurrence_key=recurrence_key)
    unit = sample_map(adata, keys, recurrence_key) if recurrence_key else None
    study = sample_map(adata, keys, keys.study) if keys.study else None
    spg = compute_programs(adata, keys, params, counts_layer=counts_layer, cache=cache,
                           data_tag=data_tag, n_jobs=n_jobs)
    return assemble(spg, params, unit_of_sample=unit, study_of_sample=study)
