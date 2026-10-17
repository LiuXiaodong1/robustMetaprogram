"""Review fixes to the Gavish port: fit diagnostics, cache code version, v2 gene order, sample
order, Process_mat parity, eligible universe, and method variants in the instability axis."""
import inspect

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mpaudit import Keys, data, decoy, gavish, nmf
from mpaudit.gavish import (MP, GavishParams, GavishResult, NMFCache, SamplePrograms, assemble, collect_top,
                            compute_programs, eligible_universe, prepare_matrix, v2_gene_order)
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
    assert set(info) == {"n_iter", "converged", "restarts"}
    assert info["converged"] and info["n_iter"] < 5000 and info["n_iter"] % 5 == 0  # tested every 5
    *_, capped = snmf_r(A, 3, rng=1, max_iter=8, return_info=True)
    assert capped == {"n_iter": 8, "converged": False, "restarts": 0}


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
        assert s.fits == {r: {"n_runs": 2, "n_capped": 2, "mean_iter": 8.0} for r in (2, 3)}
    res = assemble(spg, p)
    assert (res.n_fits, res.n_capped_fits) == (12, 12)
    again = compute_programs(adata, KEYS, p, cache=cache)  # cache hits carry the diagnostics
    assert all(again[s].fits == spg[s].fits for s in spg)

    p2 = p.replace(max_iter=5000)
    res2 = assemble(compute_programs(adata, KEYS, p2), p2)
    assert (res2.n_fits, res2.n_capped_fits) == (12, 0)


# --------------------------------------------------------------------------- cache code version

def test_cache_key_includes_code_version(monkeypatch):
    p = GavishParams()
    k0 = NMFCache.key("data", "fp", "s", 4, p)
    assert gavish.CODE_VERSION == gavish.code_version()
    monkeypatch.setattr(gavish, "CODE_VERSION", "edited")
    assert NMFCache.key("data", "fp", "s", 4, p) != k0


@pytest.mark.parametrize("target", [nmf, decoy, data.log2_cpm10, prepare_matrix, gavish._fit_job])
def test_code_version_covers_fit_sources(monkeypatch, target):
    v0 = gavish.code_version()
    real = inspect.getsource
    monkeypatch.setattr(gavish, "CODE_VERSION", None)
    monkeypatch.setattr(inspect, "getsource", lambda obj: real(obj) + ("\n# edit" if obj is target else ""))
    assert gavish.code_version() != v0


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
        collect_top(spg, GavishParams(ranks=(2,), sample_order="random"), study)


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

    res = run_instability(None, Keys(), GavishParams(), perts[:1], log=lambda *_: None)
    assert np.isnan(res.per_mp.iloc[0]["recovery_method_variant"])
