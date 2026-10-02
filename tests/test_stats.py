import numpy as np
import pandas as pd

from mpaudit.coexpr import kappa, latent_corr
from mpaudit.compare import match_sets
from mpaudit.positional import arm_enrichment, assign_arms, centromeres, known_cluster_genes
from mpaudit.typer import auroc


def test_latent_corr_recovers_correlation_despite_depth():
    rng = np.random.default_rng(0)
    n = 4000
    z = rng.normal(size=n)
    x1 = np.exp(0.6 * z + 0.4 * rng.normal(size=n)) * 2e-3
    x2 = np.exp(0.6 * z + 0.4 * rng.normal(size=n)) * 1e-3
    x3 = np.exp(0.7 * rng.normal(size=n)) * 2e-3
    s = rng.lognormal(np.log(3000), 0.6, n)  # depth varies a lot
    Y = rng.poisson(s[:, None] * np.c_[x1, x2, x3])
    C = latent_corr(Y, s)
    true12 = np.corrcoef(x1, x2)[0, 1]
    assert abs(C[0, 1] - true12) < 0.1
    assert abs(C[0, 2]) < 0.1
    naive = np.corrcoef(Y[:, 0], Y[:, 2])[0, 1]
    assert naive > 0.2  # depth induces spurious raw correlation that the estimator removes
    assert kappa(Y[:, :2], s) > 0.4


def test_match_sets_hungarian():
    ref = {"a": list("abcdefghij"), "b": list("klmnopqrst")}
    alt = {"x": list("klmnopqrsz"), "y": list("abcdefghiz"), "w": list("uvw")}
    m = match_sets(ref, alt, 0.5).set_index("ref")
    assert m.loc["a", "alt"] == "y" and m.loc["b", "alt"] == "x"
    assert m["matched"].all()


def test_arm_assignment_and_enrichment():
    cyto = pd.DataFrame({"chrom": ["1", "2", "3"], "start": [40, 40, 40], "end": [50, 50, 50],
                         "band": "c", "stain": "acen"})
    cen = centromeres(cyto)
    pos = pd.DataFrame({"chromosome": ["chr1"] * 100 + ["2"] * 100 + ["3"] * 100,
                        "start": list(range(100)) * 3}, index=[f"g{i}" for i in range(300)])
    arms = assign_arms(pos, cen)
    assert arms["g0"] == "1p" and arms["g99"] == "1q"
    sets = {"pos": [f"g{i}" for i in range(260, 290)], "rand": [f"g{i}" for i in range(0, 300, 10)]}
    res = arm_enrichment(sets, arms).set_index("mp")
    assert res.loc["pos", "best_arm"] == "3q" and res.loc["pos", "q"] < 1e-6
    assert res.loc["rand", "q"] > 0.05
    assert known_cluster_genes(["KRT5", "HLA-DRA", "MT2A", "MT-CO1", "GAPDH"]) == {"KRT5", "HLA-DRA", "MT2A"}


def test_auroc():
    a, se = auroc(np.array([0.9, 0.8, 0.1, 0.2]), np.array([1, 1, 0, 0]))
    assert a == 1.0 and se >= 0
    a, _ = auroc(np.array([0.1, 0.2, 0.9, 0.8]), np.array([1, 1, 0, 0]))
    assert a == 0.0


def test_protocol_loss_is_one_sided():
    from mpaudit.protocol import protocol_shift

    rows = []
    for t, d in zip("ABCD", [-3.0, 2.0, 0.1, -0.1]):
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", "intensity": 0.0, "prevalence": 0.0},
                 {"tumour": t, "protocol": "fixed", "mp": "m", "intensity": d, "prevalence": 0.0}]
    for t, d in zip("EF", [1.0, -1.0]):  # same-protocol replicate pairs -> null |delta| = 1
        rows += [{"tumour": t, "protocol": "fresh", "mp": "m", "intensity": 0.0, "prevalence": 0.0},
                 {"tumour": t, "protocol": "fresh", "mp": "m", "intensity": d, "prevalence": 0.0}]
    res = protocol_shift(pd.DataFrame(rows), "fresh", "fixed").iloc[0]
    assert res["frac_loss_beyond_null_intensity"] == 0.25
    assert res["frac_gain_beyond_null_intensity"] == 0.25


def test_nmf_cache_key_depends_on_data():
    from mpaudit import Keys
    from mpaudit.gavish import GavishParams, NMFCache, compute_programs
    from mpaudit.sim import simulate_cohort

    a, _ = simulate_cohort(n_tumours=3, n_twin=0, n_region=0, n_mal=60, n_tme_per_type=5, n_genes=300, seed=1)
    b = a.copy()
    b.layers["counts"] = b.layers["counts"] * 2
    keys = Keys(sample="sample", malignant="malignant")
    p = GavishParams(ranks=(2,), nrun=1, max_iter=20, min_genes_per_cell=0)
    cache = NMFCache()
    compute_programs(a, keys, p, counts_layer="counts", cache=cache)
    n1 = len(cache.mem)
    compute_programs(b, keys, p, counts_layer="counts", cache=cache)
    assert len(cache.mem) == 2 * n1  # different data -> new fits, not stale cache hits
