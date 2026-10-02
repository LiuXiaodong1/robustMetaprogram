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
