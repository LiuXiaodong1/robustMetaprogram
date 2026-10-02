import numpy as np
import pytest

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


def test_dead_components_are_skipped():
    from mpaudit.gavish import SamplePrograms, collect_top

    W = np.zeros((100, 4), dtype=np.float32)
    W[:10, 0] = 1.0
    W[10:20, 1] = 1.0
    W[20:30, 2] = 1.0  # column 3 is dead
    spg = {"S": SamplePrograms(sample="S", universe=np.array(G[:100]), W={4: W})}
    top, sample, dead = collect_top(spg, GavishParams(ranks=(4,)))
    assert dead == 1 and len(top) == 3


def test_dendrogram_order_puts_lighter_branch_first():
    from mpaudit.gavish import dendrogram_order

    ov = np.full((4, 4), 0)
    np.fill_diagonal(ov, 50)
    ov[0, 1] = ov[1, 0] = 40   # pair (0,1) is similar
    ov[2, 3] = ov[3, 2] = 30   # pair (2,3) is similar, lower column means
    order = list(dendrogram_order(ov, 50))
    assert sorted(order) == [0, 1, 2, 3]
    assert set(order[:2]) == {2, 3}  # lighter (lower colMeans) cluster first, as R's reorder()


def test_v2_post_filters():
    from mpaudit.gavish import MP, post_filter

    mps = [MP("a", G[:50], ["s1|r4|1", "s2|r4|1", "s3|r4|1"]),          # 2 studies, one singleton -> drop
           MP("b", G[:50], ["s1|r4|1", "s4|r4|1", "s2|r4|1", "s3|r4|1"]),  # 2 studies, 2 + 2 -> keep
           MP("c", G[:50], ["s1|r4|1", "s4|r4|1"])]                       # one study -> drop
    prog_sample = {p: p.split("|")[0] for m in mps for p in m.members}
    study = {"s1": "A", "s4": "A", "s2": "B", "s3": "B"}
    out = post_filter(mps, prog_sample, study, GavishParams(min_studies_per_mp=2))
    assert [m.members for m in out] == [mps[1].members]
    with pytest.raises(ValueError):
        post_filter(mps, prog_sample, None, GavishParams(min_studies_per_mp=2))
