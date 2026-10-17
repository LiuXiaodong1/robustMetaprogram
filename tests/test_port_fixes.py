"""Review fixes to the Gavish port: fit diagnostics, cache code version and cache hits, parameter
validation and early key checks, v2 gene order, sample order, Process_mat parity, eligible
universe, and method variants in the instability axis."""
import inspect

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit import Keys, data, decoy, gavish, nmf
from mpaudit.gavish import (MP, GavishParams, GavishResult, NMFCache, SamplePrograms, assemble, collect_top,
                            compute_programs, eligible_universe, order_samples, prepare_matrix, run_gavish,
                            v2_gene_order)
from mpaudit.instability import Perturbation, default_perturbations, run_instability
from mpaudit.nmf import best_of_runs, snmf_r

G = [f"G{i}" for i in range(200)]
KEYS = Keys(sample="sample", malignant="malignant", study="study")


def planted(seed: int = 0, m: int = 300, n: int = 120, k: int = 3) -> np.ndarray:
    rng = np.random.default_rng(seed)
    Wt = np.zeros((m, k))
    for j in range(k):
        Wt[40 * j:40 * (j + 1), j] = rng.uniform(1, 3, 40)
    return Wt @ rng.gamma(0.5, 1, (k, n)) + 0.05 * rng.random((m, n))


def toy_adata(n_samples: int = 3, n_cells: int = 60, n_genes: int = 200, seed: int = 0):
    """Small cohort: Poisson counts around per-sample baselines plus one shared program."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    blocks, obs = [], []
    for i in range(n_samples):
        rate = np.tile(rng.lognormal(-0.5, 1.2, n_genes), (n_cells, 1))
        rate[np.ix_(rng.random(n_cells) < 0.4, np.arange(20))] *= 4.0
        blocks.append(rng.poisson(rate * rng.lognormal(0, 0.3, (n_cells, 1))))
        obs += [{"sample": f"s{i}", "malignant": True, "study": f"st{i % 2}"}] * n_cells
    obs = pd.DataFrame(obs, index=[f"c{j}" for j in range(n_samples * n_cells)])
    return ad.AnnData(X=sp.csr_matrix(np.vstack(blocks).astype(np.float32)), obs=obs,
                      var=pd.DataFrame(index=G[:n_genes]))


# --------------------------------------------------------------------------- NMF diagnostics

def test_snmf_info_reports_convergence_and_cap():
    A = planted()
    assert len(snmf_r(A, 3, rng=1)) == 3  # default return unchanged
    *_, info = snmf_r(A, 3, rng=1, return_info=True)
    assert set(info) == {"n_iter", "converged", "restarts", "failed"}
    assert info["converged"] and not info["failed"] and info["n_iter"] < 5000 and info["n_iter"] % 5 == 0
    *_, capped = snmf_r(A, 3, rng=1, max_iter=8, return_info=True)
    assert capped == {"n_iter": 8, "converged": False, "restarts": 0, "failed": False}


def test_restart_limit_is_counted_as_failed():
    """Runs stopped by the restart limit are neither converged nor capped; they are counted apart."""
    rng = np.random.default_rng(0)
    with pytest.warns(RuntimeWarning, match="restarts"):
        W, infos = best_of_runs(rng.random((100, 40)), 6, nrun=3, beta=1e4, return_info=True)
    assert all(i["failed"] and not i["converged"] and i["restarts"] == 10 and i["n_iter"] < 5000 for i in infos)
    summary = gavish._fit_summary(infos, 5000)
    assert (summary["n_runs"], summary["n_capped"], summary["n_failed"]) == (3, 0, 3)
    spg = {"s": SamplePrograms(sample="s", universe=np.array([f"g{i}" for i in range(100)]),
                               W={6: W.astype(np.float32)}, fits={6: summary})}
    res = assemble(spg, GavishParams(ranks=(6,)))
    assert (res.n_fits, res.n_capped_fits, res.n_failed_fits) == (3, 0, 3)


def test_best_of_runs_reports_every_run():
    A = planted()
    W, infos = best_of_runs(A, 3, nrun=3, max_iter=50, return_info=True)
    assert len(infos) == 3 and all(i["n_iter"] == 50 or i["converged"] for i in infos)
    assert np.array_equal(W, best_of_runs(A, 3, nrun=3, max_iter=50))


def test_defaults_follow_the_review():
    p = GavishParams()
    assert p.max_iter == 5000 and p.sample_order == "alphabetical" and not p.v2_gene_order
    assert inspect.signature(snmf_r).parameters["check_every"].default == 5


def test_fit_diagnostics_reach_sample_programs_and_result():
    adata = toy_adata()
    p = GavishParams(ranks=(2, 3), nrun=2, max_iter=8, min_genes_per_cell=0, min_cells=10)
    cache = NMFCache()
    spg = compute_programs(adata, KEYS, p, cache=cache)
    assert list(spg) == ["s0", "s1", "s2"]
    for s in spg.values():  # max_iter=8 < second convergence test: every run is capped
        assert s.fits == {r: {"n_runs": 2, "n_capped": 2, "n_failed": 0, "mean_iter": 8.0} for r in (2, 3)}
    res = assemble(spg, p)
    assert (res.n_fits, res.n_capped_fits, res.n_failed_fits) == (12, 12, 0)

    p2 = p.replace(max_iter=5000)
    res2 = assemble(compute_programs(adata, KEYS, p2), p2)
    assert (res2.n_fits, res2.n_capped_fits) == (12, 0)


# --------------------------------------------------------------------------- cache hits

def _spy_fit_job(monkeypatch) -> list:
    """Replace gavish._fit_job by a wrapper recording (sample, ranks) of every fitting job."""
    calls, real = [], gavish._fit_job
    monkeypatch.setattr(gavish, "_fit_job", lambda job: calls.append((job[0], job[4])) or real(job))
    return calls


def test_memory_cache_hits_skip_fitting(monkeypatch):
    adata = toy_adata()
    p = GavishParams(ranks=(2, 3), nrun=2, max_iter=8, min_genes_per_cell=0, min_cells=10)
    calls = _spy_fit_job(monkeypatch)
    cache = NMFCache()
    spg = compute_programs(adata, KEYS, p, cache=cache)
    assert calls == [(s, (2, 3)) for s in ("s0", "s1", "s2")]
    again = compute_programs(adata, KEYS, p, cache=cache)
    assert len(calls) == 3  # every (sample, rank) came from the cache
    for s in spg:  # the cached objects themselves, diagnostics included
        assert again[s].universe is spg[s].universe
        assert all(again[s].W[r] is spg[s].W[r] and again[s].fits[r] is spg[s].fits[r] for r in p.ranks)
    more = compute_programs(adata, KEYS, p.replace(ranks=(2, 3, 4)), cache=cache)
    assert calls[3:] == [(s, (4,)) for s in ("s0", "s1", "s2")]  # only the missing rank is fitted
    assert all(more[s].W[2] is spg[s].W[2] and set(more[s].fits) == {2, 3, 4} for s in spg)


def test_disk_cache_hits_skip_fitting(monkeypatch, tmp_path):
    adata = toy_adata()
    p = GavishParams(ranks=(2, 3), nrun=2, max_iter=8, min_genes_per_cell=0, min_cells=10)
    calls = _spy_fit_job(monkeypatch)
    first = compute_programs(adata, KEYS, p, cache=NMFCache(tmp_path))
    assert len(calls) == 3 and len(list(tmp_path.glob("*.pkl"))) == 6
    cache = NMFCache(tmp_path)  # a fresh object: hits can only come from the files
    back = compute_programs(adata, KEYS, p, cache=cache)
    assert len(calls) == 3 and cache.mem == {}
    for s in first:
        assert np.array_equal(back[s].universe, first[s].universe)
        for r in p.ranks:
            assert back[s].W[r].dtype == np.float32 and np.array_equal(back[s].W[r], first[s].W[r])
            assert back[s].fits[r] == first[s].fits[r] == {"n_runs": 2, "n_capped": 2, "n_failed": 0,
                                                            "mean_iter": 8.0}
    monkeypatch.setattr(gavish, "CODE_VERSION", "edited")  # other code: the files are not served
    compute_programs(adata, KEYS, p, cache=NMFCache(tmp_path))
    assert len(calls) == 6


# --------------------------------------------------------------------------- cache code version

def test_cache_key_includes_code_version(monkeypatch):
    p = GavishParams()
    k0 = NMFCache.key("data", "fp", "s", 4, p)
    assert gavish.CODE_VERSION == gavish.code_version()
    monkeypatch.setattr(gavish, "CODE_VERSION", "edited")
    assert NMFCache.key("data", "fp", "s", 4, p) != k0


def _edit_source(monkeypatch, target):
    """Make inspect.getsource report ``target`` as edited, and drop the cached code versions."""
    real = inspect.getsource
    monkeypatch.setattr(gavish, "CODE_VERSION", None)
    monkeypatch.setattr(gavish, "DECOY_CODE_VERSION", None)
    monkeypatch.setattr(inspect, "getsource", lambda obj: real(obj) + ("\n# edit" if obj is target else ""))


@pytest.mark.parametrize("target", [nmf, data.log2_cpm10, prepare_matrix, gavish._fit_job, gavish._fit_summary])
def test_code_version_covers_fit_sources(monkeypatch, target):
    v0, d0 = gavish.code_version(), gavish.code_version(decoy_data=True)
    _edit_source(monkeypatch, target)
    assert gavish.code_version() != v0 and gavish.code_version(decoy_data=True) != d0  # decoy fits too


REAL_TAG, DECOY_TAG = "data", "decoy-residual-clone_est-d1-s1000"


@pytest.mark.parametrize("target", [decoy.decoy_recruitment, decoy.kappa_retained, decoy._reference_kappa, decoy])
def test_editing_decoy_recruitment_changes_no_cache_key(monkeypatch, target):
    # e.g. a docstring edit in decoy.py outside the decoy generator: no fit, real or decoy, is refitted
    p = GavishParams()
    keys0 = [NMFCache.key(tag, "fp", "s", 4, p) for tag in (REAL_TAG, DECOY_TAG)]
    _edit_source(monkeypatch, target)
    assert [NMFCache.key(tag, "fp", "s", 4, p) for tag in (REAL_TAG, DECOY_TAG)] == keys0


@pytest.mark.parametrize("target", [decoy.make_decoy, decoy._glm_mean, decoy._scrp_block, decoy._strata,
                                    data.malignant_mask])
def test_editing_the_decoy_generator_changes_only_decoy_keys(monkeypatch, target):
    p = GavishParams()
    real0, decoy0 = NMFCache.key(REAL_TAG, "fp", "s", 4, p), NMFCache.key(DECOY_TAG, "fp", "s", 4, p)
    assert real0 != decoy0 and gavish.code_version(decoy_data=True) != gavish.code_version()
    _edit_source(monkeypatch, target)
    assert NMFCache.key(REAL_TAG, "fp", "s", 4, p) == real0   # fits on real data are still served
    assert NMFCache.key(DECOY_TAG, "fp", "s", 4, p) != decoy0


def _package_functions_reached(fn) -> set:
    """Functions of mpaudit that ``fn`` calls, directly or through each other (globals named in its
    code and in the code of functions nested in it)."""
    seen, todo = set(), [fn]
    while todo:
        f = todo.pop()
        codes = [f.__code__]
        names = set()
        while codes:
            c = codes.pop()
            names |= set(c.co_names)
            codes += [k for k in c.co_consts if inspect.iscode(k)]
        for name in names:
            obj = f.__globals__.get(name)
            if inspect.isfunction(obj) and obj.__module__.startswith("mpaudit") and obj not in seen:
                seen.add(obj)
                todo.append(obj)
    return seen


def test_decoy_sources_are_everything_make_decoy_reaches():
    modules = {"data": data, "decoy": decoy}
    hashed = {getattr(modules[m], f) for m, f in (s.split(".") for s in gavish.DECOY_SOURCES)}
    assert hashed == _package_functions_reached(decoy.make_decoy) | {decoy.make_decoy}


def test_code_version_covers_cache_format(monkeypatch):
    v0 = gavish.code_version()
    monkeypatch.setattr(gavish, "CODE_VERSION", None)
    monkeypatch.setattr(gavish, "CACHE_FORMAT", gavish.CACHE_FORMAT + 1)
    assert gavish.code_version() != v0


# --------------------------------------------------------------------------- parameter checks

@pytest.mark.parametrize("field,value", [("sample_order", "studies"), ("founder_order", "dendrogram"),
                                         ("gene_universe", "hvgs")])
def test_params_reject_unknown_options(field, value):
    for build in (lambda: GavishParams(**{field: value}), lambda: GavishParams().replace(**{field: value}),
                  lambda: GavishParams.from_dict({field: value})):
        with pytest.raises(ValueError, match=f"Unknown {field}"):
            build()
    with pytest.raises(ValueError, match="Unknown GavishParams field"):
        GavishParams.from_dict({"sampel_order": "study"})


def test_params_convert_yaml_strings():
    yaml = pytest.importorskip("yaml")
    d = yaml.safe_load("tol: 1e-8\nbeta: 1e-2\neta: null\nmax_iter: 5e3\nnrun: '3'\nranks: [4, '5', 6.0]\n"
                       "center: 'false'\n")
    assert d["tol"] == "1e-8" and d["max_iter"] == "5e3"  # YAML 1.1 floats need a decimal point
    p = GavishParams.from_dict(d)
    assert (p.tol, p.beta, p.eta, p.max_iter, p.nrun, p.ranks, p.center) == (1e-8, 0.01, None, 5000, 3, (4, 5, 6), False)
    assert [type(v) for v in (p.tol, p.beta, p.max_iter, p.nrun, p.center, *p.ranks)] == [float, float, int, int, bool,
                                                                                         int, int, int]
    typed = GavishParams(tol=1e-8, beta=0.01, max_iter=5000, nrun=3, ranks=(4, 5, 6), center=False)
    assert p == typed and NMFCache.key("data", "fp", "s", 4, p) == NMFCache.key("data", "fp", "s", 4, typed)
    assert GavishParams(eta="2.5").eta == 2.5 and type(GavishParams(beta=0).beta) is float
    assert GavishParams(ranks=4).ranks == (4,)
    for bad in ({"max_iter": 2.5}, {"tol": "1e-8x"}, {"nrun": True}, {"ranks": [4, 4]}, {"center": "maybe"},
                {"nrun": 0}, {"tol": -1e-8}):
        with pytest.raises(ValueError):
            GavishParams.from_dict(bad)


def test_string_tol_fits():
    """A YAML-string tol used to crash snmf_r at the second convergence test (it=10)."""
    p = GavishParams.from_dict({"tol": "1e-8", "ranks": [2], "nrun": 1, "max_iter": 20, "min_genes_per_cell": 0,
                                "min_cells": 10})
    spg = compute_programs(toy_adata(n_samples=1), KEYS, p)
    assert spg["s0"].fits[2]["mean_iter"] == 20


_NO_STUDY = Keys(sample="sample", malignant="malignant")


@pytest.mark.parametrize("kw,keys", [({"sample_order": "study"}, _NO_STUDY), ({"v2_gene_order": True}, _NO_STUDY),
                                     ({"min_studies_per_mp": 2}, _NO_STUDY),
                                     ({"v2_gene_order": True}, Keys(sample="sample", malignant="malignant",
                                                                    study="no_such_column"))])
def test_study_options_fail_before_any_nmf(monkeypatch, kw, keys):
    adata = toy_adata()
    calls = _spy_fit_job(monkeypatch)
    p = GavishParams(ranks=(2,), nrun=1, max_iter=5, min_genes_per_cell=0, min_cells=10, **kw)
    with pytest.raises(ValueError, match="study"):
        run_gavish(adata, keys, p)
    with pytest.raises(ValueError, match="study"):
        compute_programs(adata, keys, p)
    assert calls == []
    run_gavish(adata, KEYS, p)  # with the study column the same settings run
    assert len(calls) == 3


def test_missing_recurrence_column_fails_before_any_nmf(monkeypatch):
    calls = _spy_fit_job(monkeypatch)
    with pytest.raises(ValueError, match="recurrence_key"):
        run_gavish(toy_adata(), KEYS, GavishParams(ranks=(2,), nrun=1, min_genes_per_cell=0, min_cells=10),
                   recurrence_key="tumour")
    assert calls == []


# --------------------------------------------------------------------------- v2 gene order

def test_v2_gene_order_counts_studies_not_programs():
    top = {"a1|r4|1": ["g1", "g2", "g3", "g5", "g7"],
           "a2|r4|1": ["g1", "g2", "g4", "g5", "g8"],
           "b1|r4|1": ["g2", "g3", "g9"],
           "c1|r4|1": ["g3", "g9", "g6"]}
    sample = {p: p.split("|")[0] for p in top}
    study = {"a1": "A", "a2": "A", "b1": "B", "c1": "C"}
    mp = MP("MP1", ["g1", "g2", "g3", "g9", "g5"], list(top))
    out = v2_gene_order([mp], top, sample, study)
    # studies per gene: g3 3; g2, g9 2; g1, g5 1 (both only in study A); ties alphabetical.
    # Counting programs instead would give g2, g3, g1, g5, g9.
    assert out[0].genes == ["g3", "g2", "g9", "g1", "g5"]
    assert out[0].name == "MP1" and out[0].members == mp.members
    assert mp.genes == ["g1", "g2", "g3", "g9", "g5"]  # input not modified
    with pytest.raises(ValueError, match="study"):
        v2_gene_order([mp], top, sample, None)


def test_v2_gene_order_in_assemble_needs_studies():
    with pytest.raises(ValueError, match="study"):
        assemble({}, GavishParams(v2_gene_order=True))
    assert assemble({}, GavishParams(v2_gene_order=True), study_of_sample={}).mps == []


# --------------------------------------------------------------------------- sample order

def test_sample_order_options():
    W = np.zeros((10, 2), dtype=np.float32)
    W[:5, 0] = 1.0
    W[5:, 1] = 1.0
    spg = {s: SamplePrograms(sample=s, universe=np.array(G[:10]), W={2: W}) for s in ["b2", "a9", "c1", "a1"]}
    study = {"b2": "S1", "a9": "S2", "c1": "S1", "a1": "S2"}

    def order(**kw):
        top, prog_sample, _ = collect_top(spg, GavishParams(ranks=(2,), **kw), study)
        return list(dict.fromkeys(prog_sample[p] for p in top))

    assert order() == ["a1", "a9", "b2", "c1"]
    assert order(sample_order="input") == ["b2", "a9", "c1", "a1"]
    assert order(sample_order="study") == ["b2", "c1", "a1", "a9"]
    with pytest.raises(ValueError, match="study"):
        collect_top(spg, GavishParams(ranks=(2,), sample_order="study"))
    with pytest.raises(ValueError, match="sample_order"):
        GavishParams(ranks=(2,), sample_order="random")


def test_study_order_is_not_r_file_order_but_input_reproduces_it():
    """The module docstring's example: R lists '<sample>_rank4_9_nrun10.RDS' files, so prefixed
    names sort differently from the bare names; only 'input' with R's order reproduces it."""
    samples = ["P1", "P10", "P2", "MGH10", "MGH100"]
    study = dict.fromkeys(samples, "S1")
    r_order = [f.split("_rank")[0] for f in sorted(f"{s}_rank4_9_nrun10.RDS" for s in samples)]  # C locale
    assert r_order == ["MGH100", "MGH10", "P10", "P1", "P2"]
    assert order_samples(samples, GavishParams(sample_order="study"), study) == ["MGH10", "MGH100", "P1", "P10", "P2"]
    assert order_samples(r_order, GavishParams(sample_order="input"), study) == r_order


# --------------------------------------------------------------------------- Process_mat parity

def test_prepare_matrix_drops_constant_and_zero_genes():
    rng = np.random.default_rng(0)
    n = 8
    X = np.zeros((n, 6))
    X[:, 0] = 10                          # every cell totals 100, so G0 is constant after log-CPM
    X[:, 2:] = rng.multinomial(90, [0.25] * 4, size=n)
    genes = np.array(G[:6])               # G1 is never detected
    A, universe = prepare_matrix(sp.csr_matrix(X), genes, GavishParams())
    assert set(universe) == set(G[2:6]) and A.shape[0] == len(universe)
    assert A.any(axis=1).all() and A.any(axis=0).all()
    _, universe_u = prepare_matrix(sp.csr_matrix(X), genes, GavishParams(center=False))
    assert set(universe_u) == {"G0"} | set(G[2:6])  # without centring only the zero gene goes


def test_prepare_matrix_drops_cells_left_all_zero():
    X = np.zeros((6, 5))
    X[:5, :4] = np.array([[5, 1, 3, 2], [1, 4, 2, 3], [3, 3, 1, 4], [2, 2, 5, 1], [4, 1, 2, 2]])
    X[5, 4] = 50                          # this cell only expresses a gene outside the top 4
    for center in (True, False):
        A, universe = prepare_matrix(sp.csr_matrix(X), np.array(G[:5]), GavishParams(n_universe=4, center=center))
        assert list(universe) != [] and "G4" not in universe
        assert A.shape == (len(universe), 5)


# --------------------------------------------------------------------------- eligible universe

def test_eligible_universe_is_union_of_sample_universes():
    adata = toy_adata()
    p = GavishParams(ranks=(2,), nrun=1, max_iter=5, min_genes_per_cell=0, min_cells=10, n_universe=60)
    spg = compute_programs(adata, KEYS, p)
    expected = sorted(set().union(*(set(s.universe) for s in spg.values())))
    got = eligible_universe(adata, KEYS, p)
    assert got == expected and max(len(s.universe) for s in spg.values()) < len(got) < adata.n_vars
    assert eligible_universe(adata, KEYS, p.replace(min_cells=61)) == []


# --------------------------------------------------------------------------- instability

def test_uncentred_is_a_method_variant():
    perts = default_perturbations(GavishParams())
    assert [p.kind for p in perts if p.name == "uncentred"] == ["method_variant"]
    assert all(p.kind != "method_variant" for p in default_perturbations(GavishParams(), method_variants=False))


def test_method_variants_stay_out_of_recovery_all(monkeypatch):
    from mpaudit import instability

    def fake_assemble(spg, params, **kw):  # the MP survives a seed change but not dropping centring
        mps = [MP("MP1", G[:50], [])] if params.center else []
        return GavishResult(params=params, top={}, program_sample={}, robust=[], mps=mps)

    monkeypatch.setattr(instability, "compute_programs", lambda *a, **kw: {})
    monkeypatch.setattr(instability, "assemble", fake_assemble)
    perts = [Perturbation("seed=1", "seed", {"seed": 1}),
             Perturbation("uncentred", "method_variant", {"center": False})]
    res = run_instability(None, Keys(), GavishParams(), perts, log=lambda *_: None)
    row = res.per_mp.iloc[0]
    assert (row["recovery_all"], row["recovery_seed"], row["recovery_method_variant"]) == (1.0, 1.0, 0.0)
    assert row["median_best_jaccard"] == 1.0 and (res.genes["frequency"] == 1.0).all()
    assert list(res.summary["kind"]) == ["seed", "method_variant"]
    assert {"n_fits", "n_capped_fits", "n_failed_fits"} <= set(res.summary.columns)  # fit diagnostics per run

    res = run_instability(None, Keys(), GavishParams(), perts[:1], log=lambda *_: None)
    assert np.isnan(res.per_mp.iloc[0]["recovery_method_variant"])


def test_instability_checks_every_perturbation_before_fitting(monkeypatch):
    from mpaudit import instability

    monkeypatch.setattr(instability, "compute_programs", lambda *a, **kw: pytest.fail("NMF ran"))
    for pert, match in ((Perturbation("typo", "threshold", {"intra_mn": 30}), "Unknown GavishParams field"),
                        (Perturbation("order", "order", {"sample_order": "study"}), "keys.study")):
        with pytest.raises(ValueError, match=match):
            run_instability(None, Keys(), GavishParams(), [pert], log=lambda *_: None)
