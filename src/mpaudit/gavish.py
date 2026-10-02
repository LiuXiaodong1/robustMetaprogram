"""Python port of the Gavish et al. 2023 / 3CA metaprogram (MP) pipeline.

Steps (3CA defaults in brackets; every constant is a field of :class:`GavishParams`):

1. Per sample, malignant cells with >= ``min_genes_per_cell`` detected genes [1,000]:
   log2(CPM/10 + 1); keep the top genes by mean expression [7,000]; centre each gene; set
   negatives to zero.
2. sNMF/R at each rank [4-9], keeping the best of ``nrun`` restarts [10]; each program is
   summarised by its top genes [50].
3. Robust programs: overlap >= ``intra_min`` [35] with another program of the same sample;
   max overlap >= ``inter_min`` [10] with a program from another sample; greedy within-sample
   redundancy pruning at ``intra_max`` [10].
4. Greedy MP clustering: founder = program with the most partners at >= 10 shared genes
   (needs > ``min_group_size`` [5] partners); programs join while they share >= 10 genes with
   the evolving 50-gene MP list; MP genes = top 50 by frequency, ties broken by NMF weight.

Checked against tiroshlab/3ca (``ITH_hallmarks/Generating_MPs/*.R``, ``3CA_v2/MP_analysis_new_batch.R``)
and renozao/NMF (``R/algorithms-snmf.R``). Remaining differences:

* NMF optimiser and stopping rule (HALS on the same objective; see :mod:`mpaudit.nmf`).
* Samples are matched exactly (the R code uses ``grep`` on names, which misfires when one
  sample name is a substring of another).
* Founder ties follow an average-linkage dendrogram reordered by column means, as in the R code,
  but scipy's linkage can break exact distance ties differently from R's ``hclust``.
* The cell filter (>= 1,000 detected genes) is inferred from other 3CA scripts; the cell-QC
  function the MP script calls (``Return_All_cells_Separated``) is not public.
* The R code upper-cases gene symbols; this port keeps ``var_names`` as given (only matters when
  comparing against gene lists with different capitalisation).
* The 3CA v2 manual removal of 11 MPs cannot be reproduced; label those MPs instead
  (``typer.labels`` in the P1 config).
"""
from __future__ import annotations

import dataclasses
import hashlib
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


@dataclass(frozen=True)
class GavishParams:
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
    max_iter: int = 1000
    tol: float = 1e-5
    seed: int = 0
    intra_min: int = 35
    intra_max: int = 10
    inter_min: int = 10
    inter_filter: bool = True
    min_intersect_initial: int = 10
    min_intersect_cluster: int = 10
    min_group_size: int = 5
    sample_order: str = "alphabetical"  # R iterates over NMF output files in list.files order; or "input"
    founder_order: str = "hclust"       # R ranks founders in dendrogram order; or "input"
    # 3CA v2-style post filters (off by default, i.e. Gavish 2023 core algorithm)
    min_programs_per_mp: int = 0
    min_studies_per_mp: int = 1
    drop_two_study_singleton: bool = True  # v2: with 2 studies, drop if one study gives a single NMF

    def replace(self, **kw) -> "GavishParams":
        if "ranks" in kw:
            kw["ranks"] = tuple(int(r) for r in kw["ranks"])
        return dataclasses.replace(self, **kw)

    @classmethod
    def from_dict(cls, d: dict | None) -> "GavishParams":
        return cls().replace(**dict(d or {}))


@dataclass
class SamplePrograms:
    """NMF basis matrices for one sample, over that sample's gene universe."""

    sample: str
    universe: np.ndarray
    W: dict[int, np.ndarray] = field(default_factory=dict)  # rank -> (n_universe, rank) float32
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

def prepare_matrix(counts: sp.spmatrix, genes: np.ndarray, params: GavishParams) -> tuple[np.ndarray, np.ndarray]:
    """3CA per-sample matrix: genes x cells, log2(CPM/10+1), top genes, centred, negatives zeroed."""
    L = log2_cpm10(counts)
    mean = np.asarray(L.mean(axis=0)).ravel()
    if params.gene_universe == "top_mean":
        score = mean
    elif params.gene_universe == "hvg":
        sq = np.asarray(L.multiply(L).mean(axis=0)).ravel()
        score = sq - mean**2
    else:
        raise ValueError(f"Unknown gene_universe {params.gene_universe!r}")
    keep = np.argsort(-score, kind="stable")[: params.n_universe]
    keep = keep[mean[keep] > 0]
    A = L[:, keep].toarray().T
    if params.center:
        A = A - A.mean(axis=1, keepdims=True)
        A[A < 0] = 0.0
    return A, np.asarray(genes)[keep]


def _rank_seeds(seed: int, sample: str, rank: int, nrun: int) -> list[np.random.SeedSequence]:
    base = np.random.SeedSequence([seed, zlib.crc32(sample.encode()), rank])
    return base.spawn(nrun)


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
    out = {}
    for r in ranks:
        W = best_of_runs(A, r, nrun=params.nrun, seeds=_rank_seeds(params.seed, sample, r, params.nrun),
                         beta=params.beta, eta=params.eta, max_iter=params.max_iter, tol=params.tol)
        out[r] = W.astype(np.float32)
    return sample, universe, out


class NMFCache:
    """Caches per-(sample, rank) NMF bases, in memory or (when ``directory`` is set) on disk.

    The key covers everything that changes the fit: a fingerprint of the input (cell barcodes,
    gene names, counts layer, total counts), a data tag (real or decoy), the gene universe,
    centring, solver settings and seed. Threshold, clustering and leave-study-out perturbations
    therefore reuse fits, and a cache directory cannot serve fits from a different dataset.
    With a directory, fits are not also kept in memory (atlas-scale runs would not fit).
    """

    def __init__(self, directory: str | Path | None = None):
        self.mem: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self.dir = Path(directory) if directory else None
        if self.dir:
            self.dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def key(data_tag: str, fingerprint: str, sample: str, rank: int, p: GavishParams) -> str:
        parts = (data_tag, fingerprint, sample, rank, p.gene_universe, p.n_universe, p.center,
                 p.nrun, p.seed, p.beta, p.eta, p.max_iter, p.tol)
        return hashlib.sha1(repr(parts).encode()).hexdigest()

    def get(self, key: str):
        if key in self.mem:
            return self.mem[key]
        if self.dir and (self.dir / f"{key}.pkl").exists():
            with open(self.dir / f"{key}.pkl", "rb") as fh:
                return pickle.load(fh)
        return None

    def put(self, key: str, universe: np.ndarray, W: np.ndarray) -> None:
        if not self.dir:
            self.mem[key] = (universe, W)
        else:
            with open(self.dir / f"{key}.pkl", "wb") as fh:
                pickle.dump((universe, W), fh, protocol=pickle.HIGHEST_PROTOCOL)


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
    """Run (or fetch from cache) per-sample NMF for every sample with enough malignant cells."""
    cache = cache or NMFCache()
    counts = get_counts(adata, counts_layer)
    mal = malignant_mask(adata, keys)
    obs_sample = adata.obs[keys.sample].astype(str).to_numpy()
    genes = np.asarray(adata.var_names.astype(str))
    obs_names = np.asarray(adata.obs_names.astype(str))
    gene_hash = hashlib.sha1("\x1f".join(genes).encode()).hexdigest()[:16]
    if params.min_genes_per_cell > 0:
        detected = np.diff((counts > 0).tocsr().indptr)
        mal = mal & (detected >= params.min_genes_per_cell)
    order = samples if samples is not None else list(pd.unique(obs_sample[mal]))

    result: dict[str, SamplePrograms] = {}
    fingerprints: dict[str, str] = {}
    jobs = []
    for s in order:
        idx = np.flatnonzero(mal & (obs_sample == s))
        if cell_fraction < 1 and idx.size:
            rng = np.random.default_rng([subset_seed, zlib.crc32(s.encode())])
            idx = np.sort(rng.choice(idx, size=max(1, int(round(cell_fraction * idx.size))), replace=False))
        if idx.size < params.min_cells:
            continue
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
                sp_obj.universe, sp_obj.W[r] = hit
        result[s] = sp_obj
        if missing:
            jobs.append((s, counts[idx], genes, params, tuple(missing), 1 if n_jobs > 1 else None))

    if jobs:
        if n_jobs > 1:
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                outputs = list(ex.map(_fit_job, jobs))
        else:
            outputs = [_fit_job(j) for j in jobs]
        for s, universe, Ws in outputs:
            for r, W in Ws.items():
                cache.put(NMFCache.key(data_tag, fingerprints[s], s, r, params), universe, W)
                result[s].W[r] = W
            result[s].universe = universe
            result[s]._index = None
    return result


# --------------------------------------------------------------------------- program filters

def top_genes(weights: np.ndarray, universe: np.ndarray, n: int) -> list[str]:
    order = np.argsort(-weights, kind="stable")[:n]
    return [str(g) for g in universe[order]]


def collect_top(sample_programs: dict[str, SamplePrograms], params: GavishParams) -> tuple[dict, dict, int]:
    """Top genes per program, in R's column order (sample, then rank, then program index).

    All-zero basis columns (dead components) are skipped: their "top genes" would just be the
    sample's highest-mean genes, identical across ranks and samples.
    """
    top, program_sample, dead = {}, {}, 0
    samples = sorted(sample_programs) if params.sample_order == "alphabetical" else list(sample_programs)
    for s in samples:
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


def assemble(
    sample_programs: dict[str, SamplePrograms],
    params: GavishParams,
    *,
    unit_of_sample: dict[str, str] | None = None,
    study_of_sample: dict[str, str] | None = None,
    drop_samples: set[str] | frozenset = frozenset(),
) -> GavishResult:
    """Robust-program filter and MP clustering on precomputed NMF bases."""
    kept = {s: v for s, v in sample_programs.items() if s not in drop_samples}
    top, program_sample, dead = collect_top(kept, params)
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
    return GavishResult(params=params, top=top, program_sample=program_sample, robust=list(robust_top),
                        mps=mps, n_dead_programs=dead)


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
    """
    params = params or GavishParams()
    spg = compute_programs(adata, keys, params, counts_layer=counts_layer, cache=cache,
                           data_tag=data_tag, n_jobs=n_jobs)
    unit = sample_map(adata, keys, recurrence_key) if recurrence_key else None
    study = sample_map(adata, keys, keys.study) if keys.study else None
    return assemble(spg, params, unit_of_sample=unit, study_of_sample=study)
