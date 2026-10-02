import numpy as np

from mpaudit.gavish import GavishParams, cluster_mps, robust_programs
from mpaudit.nmf import snmf_r

G = [f"G{i}" for i in range(2000)]


def genes(*ranges):
    out = []
    for a, b in ranges:
        out += G[a:b]
    return out


def test_snmf_recovers_planted_programs():
    rng = np.random.default_rng(0)
    m, n, k = 600, 300, 4
    Wt = np.zeros((m, k))
    for j in range(k):
        Wt[50 * j:50 * (j + 1), j] = rng.uniform(1, 3, 50)
    A = Wt @ rng.gamma(0.5, 1, (k, n)) + 0.05 * rng.random((m, n))
    W, H, obj = snmf_r(A, k, rng=1)
    assert np.isfinite(obj) and (W >= 0).all() and (H >= 0).all()
    tops = [set(np.argsort(-W[:, j])[:50]) for j in range(k)]
    for j in range(k):
        assert max(len(set(range(50 * j, 50 * (j + 1))) & t) for t in tops) >= 45


def test_robust_programs_intra_inter_and_redundancy():
    top = {
        "A|r4|1": genes((0, 50)),                  # shared program P
        "A|r5|1": genes((0, 45), (500, 505)),      # same program at another rank (redundant)
        "A|r4|2": genes((100, 150)),               # no within-sample partner -> dropped (intra)
        "A|r5|2": genes((200, 250)),               # within-sample partner but not recurrent -> dropped (inter)
        "A|r6|1": genes((205, 255)),
        "B|r4|1": genes((5, 55)),                  # P in sample B
        "B|r5|1": genes((5, 50), (600, 605)),
    }
    sample = {p: p.split("|")[0] for p in top}
    robust = robust_programs(top, sample, GavishParams())
    assert set(robust) == {"A|r4|1", "B|r4|1"}


def test_recurrence_unit_blocks_self_recurrence():
    top = {"A|r4|1": genes((0, 50)), "A|r5|1": genes((0, 48), (900, 902)),
           "A2|r4|1": genes((0, 50)), "A2|r5|1": genes((1, 50), (950, 951))}
    sample = {p: p.split("|")[0] for p in top}
    assert robust_programs(top, sample, GavishParams())
    assert robust_programs(top, sample, GavishParams(), unit_of_sample={"A": "T1", "A2": "T1"}) == []


def test_cluster_mps_groups_and_orders_by_frequency():
    rng = np.random.default_rng(1)
    core = list(range(0, 40))
    top = {}
    for i in range(8):  # 8 programs sharing 40 core genes plus 10 private genes
        top[f"S{i}|r4|1"] = genes(*[(c, c + 1) for c in core]) + genes((1000 + 10 * i, 1010 + 10 * i))
    for i in range(4):  # a small group: too few partners to found an MP
        top[f"T{i}|r4|1"] = genes((300, 345), (1500 + 5 * i, 1505 + 5 * i))
    weights = {p: {g: rng.random() for g in gs} for p, gs in top.items()}
    mps = cluster_mps(top, lambda p, gs: np.array([weights[p].get(g, np.nan) for g in gs]), GavishParams())
    assert len(mps) == 1
    mp = mps[0]
    assert len(mp.members) == 8 and len(mp.genes) == 50 and len(set(mp.genes)) == 50
    assert set(genes((0, 40))) <= set(mp.genes)
