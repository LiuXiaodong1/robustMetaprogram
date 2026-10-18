"""Calibration harness for criterion (c) (protocol Δlog κ losses) of an mpaudit package.

    python scripts/calibration/evaluate_c.py PKG_DIR OUT_JSON [--quick] [--workers 4] [--n-null-sets N]
        [--parts a,b,c,d] [--toy-seeds N] [--variant-seeds N] [--baseline C_JSON]

PKG_DIR is src/mpaudit, a modified copy of it, or a directory containing one. Only the public
protocol API is called: ``protocol_shift_kappa`` (toy cohorts), ``per_sample_kappa`` (realised size of
injected losses) and ``tumour_shifts`` / ``protocol_shift`` (synthetic levels; from P1-GO v1.3 with
``self_bound=False``, since synthetic κ levels have no cells to bootstrap). See README.md for the
metrics and pass criteria. The toy null (part b) also runs variants with unequal depth and small pieces
(``VARIANTS``), and part c an informational patchy program on the replicate path (``PATCHY_LAWS``).
"""
from __future__ import annotations

import argparse
import inspect
import sys
import time
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402  (sets the thread env vars before numpy)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REF, TEST = "fresh", "fixed"
NULL_KW = dict(n_tumours=10, n_twin=0, n_region=10, n_mal=200, n_tme_per_type=10, estimated_covariates=False)
CARRIER_ALL = ("cycle", "hypoxia", "stress")  # planted in every tumour; no protocol effect in the null
MIXED = ("EMT", "IFN", "CNA_7q")  # planted in a subset of tumours
OPTS: dict = {}  # filled by main(); inherited by forked workers


# --------------------------------------------------------------------------- package access

def pkg(pkg_dir: str):
    """(protocol, sim, data, report) modules of the candidate package (imported once per process)."""
    if "mpaudit.protocol" not in sys.modules:
        common.load_mpaudit(pkg_dir)
    import mpaudit.data as data
    import mpaudit.protocol as protocol
    import mpaudit.report as report
    import mpaudit.sim as sim
    common.check_modules(pkg_dir)
    return protocol, sim, data, report


def frozen_rule(report) -> dict:
    rule = report.validate_rule(None) if hasattr(report, "validate_rule") else dict(report.DEFAULT_RULE)
    return {"min_paired_tumours": int(rule["min_paired_tumours"]), "protocol_frac": float(rule["protocol_frac"])}


def crit_c_from_tumours(report, rule: dict, pt_one_mp: pd.DataFrame) -> bool:
    """Criterion (c) for one MP from its per-tumour rows (the frozen rule unless the candidate's report
    module provides ``harness_criterion_c(per_tumour_rows, rule) -> bool``)."""
    hook = getattr(report, "harness_criterion_c", None)
    if hook is not None:
        return bool(hook(pt_one_mp, rule))
    loss = pt_one_mp["loss_kappa"].dropna().astype(bool)
    return bool(len(loss) >= rule["min_paired_tumours"] and loss.mean() >= rule["protocol_frac"])


def crit_c_from_per_mp(report, rule: dict, per_mp: pd.DataFrame) -> pd.Series:
    hook = getattr(report, "harness_criterion_c_per_mp", None)
    if hook is not None:
        return pd.Series(hook(per_mp, rule), index=per_mp.index).astype(bool)
    return ((per_mp["n_paired_tumours"] >= rule["min_paired_tumours"])
            & (per_mp["frac_loss_beyond_null_kappa"] >= rule["protocol_frac"])).fillna(False).astype(bool)


def default_n_null_sets(protocol) -> int:
    return int(inspect.signature(protocol.protocol_shift_kappa).parameters["n_null_sets"].default)


def psk_kwargs(protocol) -> dict:
    return {"n_null_sets": OPTS["n_null_sets"]} if OPTS.get("n_null_sets") else {}


def run_psk(protocol, adata, keys, sets, **kw):
    res = protocol.protocol_shift_kappa(adata, keys, sets, REF, TEST, **psk_kwargs(protocol), **kw)
    return res[0], res[1], res[2]


# --------------------------------------------------------------------------- synthetic (parts a, c)

LAWS = {
    # name: (law of the MP, law of the random sets); "exchangeable" when both are the same
    "floor": ({"kind": "normal", "mu": 0.02, "sd": 0.015},) * 2,
    "weak_lognormal": ({"kind": "lognormal", "med": 0.06, "sdlog": 0.35},) * 2,
    "background": ({"kind": "normal", "mu": 0.10, "sd": 0.03},) * 2,
    "high_background": ({"kind": "normal", "mu": 0.18, "sd": 0.04},) * 2,
    # MP above its random sets (present), sets N(0.10, 0.03)
    "present_weak_1.5sd": ({"kind": "normal", "mu": 0.145, "sd": 0.03}, {"kind": "normal", "mu": 0.10, "sd": 0.03}),
    "present_3sd": ({"kind": "normal", "mu": 0.19, "sd": 0.03}, {"kind": "normal", "mu": 0.10, "sd": 0.03}),
    "present_strong": ({"kind": "normal", "mu": 0.30, "sd": 0.05}, {"kind": "normal", "mu": 0.10, "sd": 0.03}),
}
REP_CARRIER_LAW = {"kind": "lognormal", "med": 0.40, "sdlog": 0.15}
# informational (part c only, replicate paths): a program absent from a share q of the pieces (regional patchiness),
# present at lognormal(median 0.4, sdlog 0.15) and absent at N(0.02, 0.015); paired and replicate tumours exchangeable
PATCHY_LAWS = {f"patchy_q{q}": ({"kind": "patchy", "q": q, "present": REP_CARRIER_LAW,
                                 "absent": {"kind": "normal", "mu": 0.02, "sd": 0.015}},) * 2 for q in (0.1, 0.2, 0.3)}


def _draw(rng, law, rho, shape):
    if law["kind"] == "patchy":  # each piece present or absent on its own
        present = rng.random(shape + (2,)) >= law["q"]
        return np.where(present, _draw(rng, law["present"], rho, shape), _draw(rng, law["absent"], rho, shape))
    z = rng.standard_normal(shape + (2,))
    z[..., 1] = rho * z[..., 0] + np.sqrt(1 - rho**2) * z[..., 1]
    if law["kind"] == "normal":
        return law["mu"] + law["sd"] * z
    return np.exp(np.log(law["med"]) + law["sdlog"] * z)


PATHS = {
    # name: (n_rep, random sets given, replicate tumours carry the MP)
    "random_sets": (0, True, False),  # no replicate tumours: random-set null, gate on the random-set q95
    "replicates_floor_gate": (6, False, False),  # docs' calibrated case: replicate null, gate = kappa_active only
    "replicates_3": (3, False, False),  # the registered minimum of replicate tumours (min_pairs), no random sets
    "replicates_with_sets": (6, True, False),  # replicate null where >= min_pairs tumours have an active pair
    "replicates_carrier_reps": (6, True, True),  # replicate tumours carry the MP (κ ~ 0.4), paired tumours lack it
}


def synth_frames(rng, mp_law, set_law, rho, n_paired, n_rep, n_mp, n_sets, rep_mp_law=None, with_sets=True):
    """levels (one row per MP x sample) and null_levels (MP x sample x set; None without sets) for paired
    tumours P* (fresh, fixed) and replicate tumours R* (fresh, fresh). The MP's κ in replicate tumours
    follows ``rep_mp_law`` (default: its law in the paired tumours, i.e. paired and replicate tumours
    exchangeable). Each piece pair (MP or set) has correlation ``rho`` between its two pieces."""
    tum = [f"P{i}" for i in range(n_paired)] + [f"R{i}" for i in range(n_rep)]
    T = len(tum)
    k = np.concatenate([_draw(rng, mp_law, rho, (n_mp, n_paired)),
                        _draw(rng, rep_mp_law or mp_law, rho, (n_mp, n_rep))], axis=1)
    tum_o = np.array(tum, dtype=object)
    samp_o = np.array([[f"{t}_1", f"{t}_2"] for t in tum], dtype=object)  # (T, 2)
    prot_o = np.empty((T, 2), dtype=object)
    prot_o[:, 0] = REF
    prot_o[:n_paired, 1] = TEST
    prot_o[n_paired:, 1] = REF
    lv = pd.DataFrame({"mp": np.repeat(np.arange(n_mp), T * 2), "tumour": np.tile(np.repeat(tum_o, 2), n_mp),
                       "protocol": np.tile(prot_o.ravel(), n_mp), "sample": np.tile(samp_o.ravel(), n_mp),
                       "kappa": k.ravel()})
    if not with_sets:
        return lv, None
    ks = _draw(rng, set_law, rho, (n_mp, T, n_sets))
    per_t = T * n_sets * 2
    nl = pd.DataFrame({"mp": np.repeat(np.arange(n_mp), per_t),
                       "tumour": np.tile(np.repeat(tum_o, n_sets * 2), n_mp),
                       "sample": np.tile(np.broadcast_to(samp_o[:, None, :], (T, n_sets, 2)).ravel(), n_mp),
                       "protocol": np.tile(np.broadcast_to(prot_o[:, None, :], (T, n_sets, 2)).ravel(), n_mp),
                       "set": np.tile(np.tile(np.repeat(np.arange(n_sets), 2), T), n_mp),
                       "kappa": ks.ravel()})
    return lv, nl


def synthetic_kwargs(protocol) -> dict:
    """``self_bound=False`` for packages whose κ test has a self bound (P1-GO v1.3 on): synthetic κ levels have
    no cells, hence no ``kappa_se``, and the synthetic cases test the null bounds alone."""
    return {"self_bound": False} if "self_bound" in inspect.signature(protocol.tumour_shifts).parameters else {}


def job_synthetic(pkg_dir: str, job: dict) -> dict:
    protocol, _, _, report = pkg(pkg_dir)
    rule = frozen_rule(report)
    if getattr(protocol, "HARNESS_SYNTHETIC", True) is False:
        return {"skipped": "candidate sets protocol.HARNESS_SYNTHETIC = False"}
    skw = synthetic_kwargs(protocol)
    n_sets = OPTS.get("n_null_sets") or default_n_null_sets(protocol)
    mp_law, set_law = {**LAWS, **PATCHY_LAWS}[job["law"]]
    n_rep, with_sets, rep_carrier = PATHS[job["path"]]
    T = job["n_paired"] + n_rep
    per_mp_rows = T * (n_sets + 1) * 2 if with_sets else T * 2
    n_mp = int(min(job["n_mp"], max(40, job["row_budget"] // per_mp_rows)))
    batch = int(max(20, min(n_mp, 2_000_000 // per_mp_rows)))
    rng = np.random.default_rng(job["seed"])
    rep_law = REP_CARRIER_LAW if rep_carrier else None
    per_mp_cnt, crit, kinds = [], [], []
    done = 0
    while done < n_mp:
        b = min(batch, n_mp - done)
        lv, nl = synth_frames(rng, mp_law, set_law, job["rho"], job["n_paired"], n_rep, b, n_sets, rep_law, with_sets)
        sh = protocol.tumour_shifts(lv, REF, TEST, null_levels=nl, **skw)
        per = protocol.protocol_shift(lv, REF, TEST, null_levels=nl, **skw)
        sh = sh[sh["tumour"].astype(str).str.startswith("P")]
        t = sh["loss_kappa"].notna().to_numpy()
        lo = _tb(sh["loss_kappa"]).to_numpy()
        ga = _tb(sh["gain_kappa"]).to_numpy()
        g = pd.DataFrame({"mp": sh["mp"].to_numpy(), "n": 1, "t": t, "l": lo & t, "g": ga & t}).groupby("mp").sum()
        g = g.reindex(per["mp"].to_numpy(), fill_value=0)
        per_mp_cnt.append(g[["n", "t", "l", "g"]].to_numpy())
        crit.append(crit_c_from_per_mp(report, rule, per).to_numpy())
        kinds.append(per["null_type"].astype(str).to_numpy())
        done += b
    c = np.concatenate(per_mp_cnt)
    crit = np.concatenate(crit)
    kinds = np.concatenate(kinds)
    res = {
        "part": job["part"], "law": job["law"], "rho": job["rho"], "path": job["path"],
        "n_paired": job["n_paired"], "n_rep": n_rep, "random_sets_given": with_sets, "n_mp": int(n_mp),
        "n_sets": int(n_sets) if with_sets else 0,
        "null_types": {k: int(v) for k, v in zip(*np.unique(kinds, return_counts=True))},
        "tested_frac": common.rate(c[:, 1], c[:, 0]),
        "loss_uncond": common.rate(c[:, 2], c[:, 0]),
        "loss_given_tested": common.rate(c[:, 2], c[:, 1]),
        "gain_given_tested": common.rate(c[:, 3], c[:, 1]),
        "p_crit_c": common.binom(int(crit.sum()), int(crit.size)),
        "p_powered": common.binom(int((c[:, 1] >= rule["min_paired_tumours"]).sum()), int(c.shape[0])),
    }
    p = res["p_crit_c"]["p"]
    res["p_ge3_of_20_mps_meet_c"] = None if p is None else round(common.p_at_least(3, 20, p), 5)
    if n_rep:  # split by the null each MP actually got
        for kind in ("replicates", "random_sets"):
            sel = kinds == kind
            res[f"loss_given_tested_on_{kind}_null"] = common.rate(c[sel, 2], c[sel, 1])
    return res


def synthetic_jobs(quick: bool) -> list[dict]:
    jobs = []
    budget = 4_000_000 if quick else 12_000_000
    n_mp_a = 400 if quick else 1000
    seed = 100
    # part a: per-tumour rates with 10 paired tumours
    for law in LAWS:
        for rho in (0.0, 0.5):
            for path in PATHS:
                if path == "replicates_carrier_reps" and law not in ("floor", "weak_lognormal", "background"):
                    continue
                seed += 1
                jobs.append(dict(part="a", law=law, rho=rho, path=path, n_paired=10, n_mp=n_mp_a,
                                 row_budget=budget, seed=seed,
                                 cost=1.0 if path in ("replicates_floor_gate", "replicates_3") else 3.0,
                                 name=f"synth_a {law} rho{rho} {path}"))
    # part c: P(criterion c) against the number of paired tumours
    for law, rho in (("weak_lognormal", 0.0), ("background", 0.0), ("background", 0.5), ("present_weak_1.5sd", 0.0)):
        for n_paired in (4, 10, 40):
            for path in ("random_sets", "replicates_floor_gate", "replicates_3", "replicates_with_sets"):
                seed += 1
                n_mp = {4: 2000, 10: 1000, 40: 400}[n_paired] // (3 if quick else 1)
                jobs.append(dict(part="c", law=law, rho=rho, path=path, n_paired=n_paired, n_mp=n_mp,
                                 row_budget=budget * (2 if n_paired == 40 else 1), seed=seed,
                                 cost=(0.3 if path.startswith("replicates_") and path != "replicates_with_sets"
                                       else 1.0) * (1.0 + n_paired / 10),
                                 name=f"synth_c {law} rho{rho} {path} P{n_paired}"))
    # informational: a patchy program on the replicate paths (C2-info)
    for law in PATCHY_LAWS:
        for n_paired in (4, 10):
            for path in ("replicates_3", "replicates_floor_gate"):
                seed += 1
                n_mp = {4: 2000, 10: 1000}[n_paired] // (3 if quick else 1)
                jobs.append(dict(part="c", law=law, rho=0.0, path=path, n_paired=n_paired, n_mp=n_mp,
                                 row_budget=budget, seed=seed, cost=0.3 * (1.0 + n_paired / 10),
                                 name=f"synth_c {law} rho0.0 {path} P{n_paired}"))
    return jobs


# --------------------------------------------------------------------------- toy cohorts (parts b, d)

def simulate(sim, seed: int, mag: str, **kw):
    if mag == "first":
        with sim.axis_magnitudes(*sim.FIRST_MAGNITUDES):
            return sim.simulate_cohort(seed=seed, **kw)
    return sim.simulate_cohort(seed=seed, **kw)


def relabel(adata, paired: list[str]):
    """Copy of ``adata`` whose R2 pieces of the ``paired`` tumours are labelled with the test protocol."""
    a = adata.copy()
    proto = a.obs["protocol"].astype(str).to_numpy().copy()
    hit = a.obs["tumour"].astype(str).isin(paired).to_numpy() & (a.obs["region"].astype(str).to_numpy() == "R2")
    proto[hit] = TEST
    a.obs["protocol"] = proto
    return a


def null_sets(truth, adata, seed: int, owners: list[str]) -> dict:
    """Planted sets used in the toy null, 3 random 40-gene sets of unplanted genes and private programs."""
    sets = {k: list(truth["sets"][k]) for k in (*CARRIER_ALL, *MIXED) if k in truth["sets"]}
    planted = {g for s in truth["sets"].values() for g in s}
    for ps in truth["private"].values():
        for p in ps:
            planted |= set(p)
    free = np.array([g for g in adata.var_names.astype(str) if g.startswith("G") and g not in planted])
    rng = np.random.default_rng(1000 + seed)
    for j in range(3):
        sets[f"random{j}"] = list(rng.choice(free, 40, replace=False))
    for t in owners:
        sets[f"private_{t}"] = list(truth["private"][t][0])
    return sets


def carrier_of(truth, mp: str, tumour: str) -> bool:
    car = {"EMT": truth["tumours_emt"], "IFN": truth["tumours_ifn"], "CNA_7q": truth["tumours_cna"]}
    if mp in car:
        return tumour in car[mp]
    if mp.startswith("private_"):
        return tumour == mp.split("_", 1)[1]
    if mp.startswith("random"):
        return False
    return True


def kind_of(mp: str) -> str:
    if mp.startswith("random"):
        return "random"
    if mp.startswith("private_"):
        return "private"
    if mp in MIXED:
        return "mixed"
    return "carrier_all"


KEEP_COLS = ["mp", "tumour", "null_type", "null_type_kappa", "kappa_reference", "kappa_test", "delta_log_kappa",
             "null_lo_kappa", "kappa_active", "kappa_active_reference", "kappa_active_test", "loss_kappa", "gain_kappa"]


def _slim(pt: pd.DataFrame) -> list[dict]:
    cols = [c for c in KEEP_COLS if c in pt]
    d = pt[cols].copy()
    for c in ("kappa_active", "kappa_active_reference", "kappa_active_test", "loss_kappa", "gain_kappa"):
        if c in d:
            d[c] = d[c].astype(object).where(d[c].notna(), None)
    return d.to_dict("records")


# toy-null variants (part b, calibrated magnitudes): unequal depth or small pieces, with no protocol effect
VARIANTS = ("test_depth_0.5", "test_depth_0.25", "n_mal_60", "test_cut_50")


def vary(adata, variant: str | None, seed: int):
    """The pseudo-paired cohort under a variant: the test pieces' counts binomially thinned to a share of their depth
    (every cell of those pieces, gene by gene), or their malignant cells cut to 50 (``n_mal_60`` is simulated)."""
    import scipy.sparse as sp
    if not variant or variant == "n_mal_60":
        return adata
    test = (adata.obs["protocol"].astype(str) == TEST).to_numpy()
    if variant == "test_cut_50":
        mal = adata.obs["malignant"].to_numpy(bool)
        keep = np.ones(adata.n_obs, bool)
        for s in pd.unique(adata.obs["sample"][test]):
            rows = np.flatnonzero((adata.obs["sample"] == s).to_numpy() & mal)
            keep[rows[50:]] = False
        return adata[keep].copy()
    frac = float(variant.rsplit("_", 1)[1])
    X = adata.X.tocsr().astype(np.float64)
    rows = np.flatnonzero(test)
    sub = X[rows].tocoo()
    sub.data = np.random.default_rng([seed, 77]).binomial(sub.data.astype(np.int64), frac).astype(np.float64)
    keep = sp.diags((~test).astype(np.float64))
    put = sp.csr_matrix((sub.data, (rows[sub.row], sub.col)), shape=X.shape)
    a = adata.copy()
    a.X = sp.csr_matrix(keep @ X + put).astype(np.float32)
    a.X.eliminate_zeros()
    if "counts" in a.layers:
        a.layers["counts"] = a.X.copy()
    return a


def job_toy_null(pkg_dir: str, job: dict) -> dict:
    protocol, sim, data, report = pkg(pkg_dir)
    rule = frozen_rule(report)
    seed, mag, variant = job["seed"], job["mag"], job.get("variant")
    adata, truth = simulate(sim, seed, mag, **{**NULL_KW, **({"n_mal": 60} if variant == "n_mal_60" else {})})
    keys = data.Keys(sample="sample", tumour="tumour", protocol="protocol")
    tumours = sorted(adata.obs["tumour"].astype(str).unique())
    out = {"seed": seed, "mag": f"{mag}+{variant}" if variant else mag, "regimes": {}, "rule": rule}
    # A: every tumour pseudo-paired (R2 relabelled fixed) -> random-set path; B: T01-T05 paired, T06-T10 replicates
    for regime, paired, owners in (("A_random_sets", tumours, [tumours[0], tumours[5]]),
                                   ("B_replicates", tumours[:5], [tumours[0], tumours[5]])):
        sets = null_sets(truth, adata, seed, owners)
        per_mp, pt, _ = run_psk(protocol, vary(relabel(adata, paired), variant, seed), keys, sets)
        pt = pt.copy()
        pt["carrier"] = [carrier_of(truth, m, t) for m, t in zip(pt["mp"], pt["tumour"])]
        crit_pkg = dict(zip(per_mp["mp"], crit_c_from_per_mp(report, rule, per_mp)))
        crit_mine = {m: crit_c_from_tumours(report, rule, g) for m, g in pt.groupby("mp")}
        halves = {}
        if regime == "A_random_sets":
            for h, ts in (("h1", tumours[:5]), ("h2", tumours[5:])):
                halves[h] = {m: crit_c_from_tumours(report, rule, g[g["tumour"].isin(ts)]) for m, g in pt.groupby("mp")}
        out["regimes"][regime] = {
            "rows": [dict(r, carrier=bool(c)) for r, c in zip(_slim(pt), pt["carrier"])],
            "per_mp": per_mp[[c for c in ("mp", "null_type", "n_null_tumours", "n_paired_tumours",
                                          "frac_loss_beyond_null_kappa") if c in per_mp]].to_dict("records"),
            "crit_c": {m: bool(v) for m, v in crit_pkg.items()},
            "crit_c_from_tumours": crit_mine, "crit_c_halves": halves,
            "aggregation_matches": all(crit_pkg.get(m) == v for m, v in crit_mine.items()),
        }
    return out


def inject_loss(adata, truth, program: str, tumours: list[str], frac: float, seed: int):
    """In the test-protocol malignant pieces of ``tumours``, permute the counts of a random ``frac`` of the
    program's genes across the piece's malignant cells (each gene independently), in X and layers['counts']."""
    import scipy.sparse as sp
    genes = np.array(truth["sets"][program])
    pos = {g: i for i, g in enumerate(adata.var_names.astype(str))}
    gidx = np.array([pos[g] for g in genes if g in pos])
    obs = adata.obs
    X = adata.X.toarray() if sp.issparse(adata.X) else np.array(adata.X)
    for s in sorted(obs.loc[obs["tumour"].astype(str).isin(tumours) & (obs["protocol"].astype(str) == TEST),
                            "sample"].astype(str).unique()):
        rows = np.flatnonzero((obs["sample"].astype(str) == s).to_numpy() & obs["malignant"].to_numpy(bool))
        rng = np.random.default_rng([seed, zlib.crc32(f"{program}|{s}".encode())])
        chosen = rng.choice(gidx, int(round(frac * gidx.size)), replace=False)
        for g in chosen:
            X[rows, g] = X[rng.permutation(rows), g]
    Xs = sp.csr_matrix(X.astype(np.float32))
    adata.X = Xs
    if "counts" in adata.layers:
        adata.layers["counts"] = Xs.copy()
    return adata


def job_inject(pkg_dir: str, job: dict) -> dict:
    """Power: no protocol effect except an injected loss of cycle (all tumours) and EMT (carriers) in the
    pseudo-fixed pieces, on the toy-null design (calibrated magnitudes)."""
    protocol, sim, data, report = pkg(pkg_dir)
    rule = frozen_rule(report)
    seed, frac = job["seed"], job["frac"]
    adata, truth = simulate(sim, seed, "cal", **NULL_KW)
    keys = data.Keys(sample="sample", tumour="tumour", protocol="protocol")
    tumours = sorted(adata.obs["tumour"].astype(str).unique())
    sets = {"cycle": truth["sets"]["cycle"], "EMT": truth["sets"]["EMT"]}
    out = {"seed": seed, "frac": frac, "regimes": {}}
    for regime, paired in (("A_random_sets", tumours), ("B_replicates", tumours[:5])):
        base = relabel(adata, paired)
        inj = base.copy()
        inject_loss(inj, truth, "cycle", paired, frac, seed)
        inject_loss(inj, truth, "EMT", [t for t in paired if t in truth["tumours_emt"]], frac, seed)
        lv0 = protocol.per_sample_kappa(base, keys, sets, REF, TEST)
        lv1 = protocol.per_sample_kappa(inj, keys, sets, REF, TEST)
        m = lv0.merge(lv1, on=["sample", "mp"], suffixes=("_0", "_1"))
        m = m[m["protocol_0"] == TEST]
        m["tumour"] = m["tumour_0"]
        m["carrier"] = [carrier_of(truth, a, t) for a, t in zip(m["mp"], m["tumour"])]
        ratio = {mp_: m.loc[(m["mp"] == mp_) & m["carrier"], "kappa_1"].to_numpy()
                 / m.loc[(m["mp"] == mp_) & m["carrier"], "kappa_0"].to_numpy() for mp_ in sets}
        per_mp, pt, _ = run_psk(protocol, inj, keys, sets)
        pt = pt.copy()
        pt["carrier"] = [carrier_of(truth, a, t) for a, t in zip(pt["mp"], pt["tumour"])]
        out["regimes"][regime] = {
            "kappa_ratio_injected": {k: [float(x) for x in v] for k, v in ratio.items()},
            "rows": [dict(r, carrier=bool(c)) for r, c in zip(_slim(pt), pt["carrier"])],
            "crit_c": {a: bool(v) for a, v in zip(per_mp["mp"], crit_c_from_per_mp(report, rule, per_mp))},
            "per_mp": per_mp[[c for c in ("mp", "null_type", "n_paired_tumours", "frac_loss_beyond_null_kappa")
                              if c in per_mp]].to_dict("records"),
        }
    return out


def job_default(pkg_dir: str, job: dict) -> dict:
    """The default toy cohort (real protocol effect: stress absent under the fixed protocol)."""
    protocol, sim, data, report = pkg(pkg_dir)
    rule = frozen_rule(report)
    adata, truth = sim.simulate_cohort(seed=job["seed"])
    keys = data.Keys(sample="sample", tumour="tumour", protocol="protocol")
    per_mp, pt, _ = run_psk(protocol, adata, keys, dict(truth["sets"]))
    pt = pt.copy()
    pt["carrier"] = [carrier_of(truth, m, t) for m, t in zip(pt["mp"], pt["tumour"])]
    return {"seed": job["seed"], "rows": [dict(r, carrier=bool(c)) for r, c in zip(_slim(pt), pt["carrier"])],
            "crit_c": {m: bool(v) for m, v in zip(per_mp["mp"], crit_c_from_per_mp(report, rule, per_mp))},
            "per_mp": per_mp[[c for c in ("mp", "null_type", "n_paired_tumours", "frac_loss_beyond_null_kappa")
                              if c in per_mp]].to_dict("records")}


# --------------------------------------------------------------------------- aggregation

def _rows_df(results: list[dict], key="regimes") -> pd.DataFrame:
    recs = []
    for r in results:
        for regime, v in r[key].items():
            for row in v["rows"]:
                recs.append(dict(row, regime=regime, seed=r["seed"], mag=r.get("mag", "cal")))
    return pd.DataFrame(recs)


def _tb(s: pd.Series) -> pd.Series:
    """True where the value is True; False for False, None and NA."""
    return s.astype(object).eq(True).astype(bool)


def _rates(d: pd.DataFrame, cluster="seed") -> dict:
    t = d["loss_kappa"].notna()
    g = pd.DataFrame({"c": d[cluster].to_numpy(), "n": 1, "t": t.to_numpy(),
                      "l": (_tb(d["loss_kappa"]) & t).to_numpy(),
                      "g": (_tb(d["gain_kappa"]) & t).to_numpy()}).groupby("c").sum()
    return {"n_pairs": int(len(d)), "tested_frac": common.rate(g.t, g.n), "loss_given_tested": common.rate(g.l, g.t),
            "gain_given_tested": common.rate(g.g, g.t), "loss_uncond": common.rate(g.l, g.n)}


def crit_category(m: str, regime: str) -> str:
    """random | private_absent_everywhere (owner is a replicate tumour) | private_owner_paired | mixed | carrier_all"""
    if m.startswith("private_"):
        return "private_absent_everywhere" if (regime.startswith("B") and m.endswith("T06")) else "private_owner_paired"
    return kind_of(m)


CRIT_CATEGORIES = ("random", "private_absent_everywhere", "private_owner_paired", "mixed", "carrier_all")


def summarise_toy_null(results: list[dict]) -> dict:
    d = _rows_df(results)
    d["kind"] = d["mp"].map(kind_of)
    out = {}
    for (mag, regime), x in d.groupby(["mag", "regime"]):
        absent = x[~x["carrier"]]
        res = {"absent_all": _rates(absent), "carrier_all": _rates(x[x["carrier"]])}
        for kind, lab in (("random", "absent_random_sets"), ("private", "absent_private"),
                          ("mixed", "absent_noncarrier_EMT_IFN_CNA")):
            res[lab] = _rates(absent[absent["kind"] == kind])
        if regime == "B_replicates":
            nt = absent["null_type_kappa"] if "null_type_kappa" in absent else absent["null_type"]
            res["absent_on_replicate_null"] = _rates(absent[nt == "replicates"])
            res["absent_fallback_random_sets"] = _rates(absent[nt == "random_sets"])
        # side of the gate among absent tested pairs
        if {"kappa_active_reference", "kappa_active_test"} <= set(absent):
            tt = absent[absent["loss_kappa"].notna()]
            ref_only = tt["kappa_active_reference"].eq(True) & ~tt["kappa_active_test"].eq(True)
            res["absent_ref_side_only"] = _rates(tt[ref_only])
            res["absent_both_sides"] = _rates(tt[tt["kappa_active_reference"].eq(True) & tt["kappa_active_test"].eq(True)])
        out[f"{mag}/{regime}"] = res
    # criterion (c)
    crit = {}
    for r in results:
        for regime, v in r["regimes"].items():
            c = crit.setdefault(f"{r['mag']}/{regime}", {})
            for m, val in v["crit_c"].items():
                c.setdefault(crit_category(m, regime), []).append(bool(val))
                c.setdefault("set:" + (m.rstrip("0123456789") if m.startswith("random") else m), []).append(bool(val))
            for h, hv in v.get("crit_c_halves", {}).items():
                for m, val in hv.items():
                    c.setdefault(crit_category(m, regime) + "@5_paired_tumours(halves)", []).append(bool(val))
    crit_out = {k: {kk: common.binom(sum(vv), len(vv)) for kk, vv in sorted(v.items())} for k, v in crit.items()}
    # reference: expected number of (set, run) meeting (c) if each tested tumour were a loss with probability 0.05
    from scipy import stats
    rule = results[0].get("rule", {"min_paired_tumours": 4, "protocol_frac": 0.5}) if results else {}
    exp = {}
    for (mag, regime, seed, mp_), x in d.groupby(["mag", "regime", "seed", "mp"]):
        n_t = int(x["loss_kappa"].notna().sum())
        p0 = 0.0
        if n_t >= rule["min_paired_tumours"]:
            p0 = float(stats.binom.sf(int(np.ceil(rule["protocol_frac"] * n_t)) - 1, n_t, 0.05))
        exp[f"{mag}/{regime}"] = exp.get(f"{mag}/{regime}", 0.0) + p0
        if regime.startswith("A"):  # halves T01-T05 and T06-T10 (5 paired tumours each)
            ts = sorted(x["tumour"].unique())
            for half in (ts[:5], ts[5:]):
                n_h = int(x.loc[x["tumour"].isin(half), "loss_kappa"].notna().sum())
                if n_h >= rule["min_paired_tumours"]:
                    key = f"{mag}/{regime}@5_paired_tumours(halves)"
                    exp[key] = exp.get(key, 0.0) + float(
                        stats.binom.sf(int(np.ceil(rule["protocol_frac"] * n_h)) - 1, n_h, 0.05))
    mismatch = sum(not v.get("aggregation_matches", True) for r in results for v in r["regimes"].values())
    return {"per_tumour": out, "p_crit_c": crit_out, "crit_c_expected_at_5pct": exp,
            "per_mp_aggregation_mismatches": mismatch,
            "n_runs": len(results)}


def summarise_inject(results: list[dict]) -> dict:
    out = {}
    for frac in sorted({r["frac"] for r in results}):
        rs = [r for r in results if r["frac"] == frac]
        d = _rows_df(rs)
        for regime, x in d.groupby("regime"):
            res = {}
            for mp_ in ("cycle", "EMT"):
                y = x[(x["mp"] == mp_)]
                car, non = y[y["carrier"]], y[~y["carrier"]]
                det = _tb(car["loss_kappa"]).groupby(car["seed"]).sum()
                res[mp_] = {"injected_tumours_detected": common.rate(det.to_numpy(),
                                                                     car.groupby("seed").size().to_numpy()),
                            "carriers": _rates(car), "noncarriers": _rates(non) if len(non) else None}
                ratios = np.concatenate([np.array(r["regimes"][regime]["kappa_ratio_injected"][mp_]) for r in rs])
                res[mp_]["kappa_ratio_injected_piece"] = {
                    "median": float(np.median(ratios)) if ratios.size else None,
                    "q10": float(np.quantile(ratios, .1)) if ratios.size else None,
                    "q90": float(np.quantile(ratios, .9)) if ratios.size else None}
                cc = [r["regimes"][regime]["crit_c"].get(mp_, False) for r in rs]
                res[mp_]["p_crit_c"] = common.binom(sum(cc), len(cc))
            out[f"f{frac}/{regime}"] = res
    return out


def summarise_default(results: list[dict]) -> dict:
    d = pd.DataFrame([dict(row, seed=r["seed"]) for r in results for row in r["rows"]])
    st = d[d["mp"] == "stress"]
    det = _tb(st["loss_kappa"]).groupby(st["seed"]).sum()
    others = d[d["mp"] != "stress"]
    false_c = {r["seed"]: sorted(m for m, v in r["crit_c"].items() if v and m != "stress") for r in results}
    nc = others[~others["carrier"]]
    return {
        "stress_tumours_detected": common.rate(det.to_numpy(), st.groupby("seed").size().to_numpy()),
        "stress_p_crit_c": common.binom(sum(r["crit_c"].get("stress", False) for r in results), len(results)),
        "other_planted_carrier_tests": _rates(others[others["carrier"]]),
        "other_planted_noncarrier_tests": _rates(nc),
        "noncarrier_losses": nc[_tb(nc["loss_kappa"])][["seed", "mp", "tumour"]].to_dict("records"),
        "crit_c_non_stress_sets_by_seed": false_c,
        "n_seeds": len(results),
    }


def summarise_synthetic(results: list[dict]) -> dict:
    a, c = {}, {}
    for r in results:
        if "skipped" in r:
            return {"skipped": r["skipped"]}
        key = f"{r['law']}|rho{r['rho']}|{r['path']}|P{r['n_paired']}"
        (a if r.get("part") == "a" else c)[key] = r
    return {"a_per_tumour": a, "c_p_crit_c": c}


# --------------------------------------------------------------------------- headline / checks

EXCHANGEABLE = ("floor", "weak_lognormal", "background", "high_background")


def headline(out: dict) -> dict:
    h = {}
    sa = out.get("synthetic", {}).get("a_per_tumour", {})
    for path in PATHS:
        rows = {k: v for k, v in sa.items() if v["path"] == path}
        if not rows:
            continue
        h[f"synthetic_a[{path}] loss|tested by law|rho"] = {
            f"{v['law']}|rho{v['rho']}": v["loss_given_tested"]["p"] for v in rows.values()}
        h[f"synthetic_a[{path}] tested fraction by law|rho"] = {
            f"{v['law']}|rho{v['rho']}": v["tested_frac"]["p"] for v in rows.values()}
        ex = [v["loss_given_tested"]["p"] for v in rows.values() if v["law"] in EXCHANGEABLE
              and v["loss_given_tested"]["p"] is not None]
        if ex:
            h[f"synthetic_a[{path}] max loss|tested, exchangeable MP"] = max(ex)
        if path == "replicates_with_sets":
            h[f"synthetic_a[{path}] loss|tested on MPs that got the replicate null"] = {
                f"{v['law']}|rho{v['rho']}": (v.get("loss_given_tested_on_replicates_null") or {}).get("p")
                for v in rows.values()}
    tn = out.get("toy_null", {}).get("per_tumour", {})
    for k, v in tn.items():
        h[f"toy_null[{k}] absent: tested fraction / loss|tested / gain|tested"] = [
            v["absent_all"]["tested_frac"]["p"], v["absent_all"]["loss_given_tested"],
            v["absent_all"]["gain_given_tested"]["p"]]
        h[f"toy_null[{k}] carrier loss|tested"] = v["carrier_all"]["loss_given_tested"]
        if "absent_ref_side_only" in v:
            h[f"toy_null[{k}] absent loss|tested when active on the reference side only"] = \
                v["absent_ref_side_only"]["loss_given_tested"]["p"]
    for k, v in out.get("toy_null", {}).get("p_crit_c", {}).items():
        h[f"toy_null[{k}] crit (c) met, k/n by category"] = {kk: f"{vv['k']}/{vv['n']}" for kk, vv in v.items()
                                                             if not kk.startswith("set:")}
    for k, v in out.get("synthetic", {}).get("c_p_crit_c", {}).items():
        h[f"synthetic_c[{k}] P(crit c)"] = v["p_crit_c"]["p"]
    pw = out.get("power", {})
    if "default" in pw:
        h["power default cohort: stress tumours detected"] = pw["default"]["stress_tumours_detected"]
        h["power default cohort: stress P(crit c)"] = pw["default"]["stress_p_crit_c"]
        h["power default cohort: losses in non-carrier tumours"] = pw["default"]["noncarrier_losses"]
    for regime, v in pw.get("inject", {}).items():
        for mp_ in ("cycle", "EMT"):
            h[f"power inject[{regime}] {mp_}: detected / P(crit c) / kappa ratio"] = [
                v[mp_]["injected_tumours_detected"], v[mp_]["p_crit_c"], v[mp_]["kappa_ratio_injected_piece"]["median"]]
    return h


def _chk(res, name, val, thr, op, **extra):
    if val is None or thr is None:
        res[name] = {"value": None, "pass": None, **extra}
        return
    ok = val <= thr if op == "<=" else val >= thr
    res[name] = {"value": val, "threshold": f"{op} {thr}", "pass": bool(ok), **extra}


def checks(out: dict, baseline: dict | None = None) -> dict:
    """Pass/fail against the targets in README.md (pass None when the metric is missing or informational)."""
    res = {}
    sa = out.get("synthetic", {}).get("a_per_tumour", {})
    for path in ("random_sets", "replicates_floor_gate", "replicates_3", "replicates_with_sets",
                 "replicates_carrier_reps"):
        rows = [v for v in sa.values() if v["path"] == path and v["loss_given_tested"]["p"] is not None]
        if not rows:
            continue
        worst = max(rows, key=lambda v: v["loss_given_tested"]["p"])
        tag = "C1" if path != "replicates_carrier_reps" else "C1-info"
        _chk(res, f"{tag} synthetic max loss|tested over laws and rho [{path}]", worst["loss_given_tested"]["p"],
             0.065, "<=", worst=f"{worst['law']}|rho{worst['rho']}", se=worst["loss_given_tested"]["se"])
    strong = [v for v in sa.values() if v["law"] == "present_strong" and v["path"] != "replicates_carrier_reps"]
    if strong:
        _chk(res, "C1p synthetic present_strong MP tested fraction (min over paths)",
             min(v["tested_frac"]["p"] for v in strong), 0.99, ">=")
    sc = out.get("synthetic", {}).get("c_p_crit_c", {})
    for tag, keep in (("C2", lambda law: law not in PATCHY_LAWS), ("C2-info", lambda law: law in PATCHY_LAWS)):
        rows = [v for v in sc.values() if v["p_crit_c"]["p"] is not None and keep(v["law"])]
        if not rows:
            continue
        worst = max(rows, key=lambda v: v["p_crit_c"]["p"])
        what = ("for an unaffected MP (laws x P in 4,10,40 x paths)" if tag == "C2" else
                "for a patchy program (absent from a share q of pieces) on the replicate paths")
        _chk(res, f"{tag} synthetic max P(crit c) {what}", worst["p_crit_c"]["p"], 0.03, "<=",
             worst=f"{worst['law']}|rho{worst['rho']}|{worst['path']}|P{worst['n_paired']}")
    for k, x in out.get("toy_null", {}).get("per_tumour", {}).items():
        _chk(res, f"C3 toy null absent loss|tested [{k}]", x["absent_all"]["loss_given_tested"]["p"], 0.065, "<=",
             se=x["absent_all"]["loss_given_tested"]["se"], n=x["absent_all"]["loss_given_tested"]["n"])
        small = "n_mal_60" in k or "test_cut_50" in k  # small pieces: informational (carriers' own noise is larger)
        _chk(res, f"{'C4-info' if small else 'C4'} toy null carrier loss|tested [{k}]",
             x["carrier_all"]["loss_given_tested"]["p"], 0.02, "<=", n=x["carrier_all"]["loss_given_tested"]["n"])
    exp = out.get("toy_null", {}).get("crit_c_expected_at_5pct", {})
    for k, x in out.get("toy_null", {}).get("p_crit_c", {}).items():
        for suffix in ("", "@5_paired_tumours(halves)"):
            cats = [c + suffix for c in CRIT_CATEGORIES]
            if not any(c in x for c in cats):
                continue
            tot = sum(vv["k"] for kk, vv in x.items() if kk in cats)
            n = sum(vv["n"] for kk, vv in x.items() if kk in cats)
            e0 = exp.get(k + suffix, 0.0)
            thr = int(np.floor(max(1.0, e0 + 2 * np.sqrt(e0))))
            _chk(res, f"C5 toy null set x run combinations meeting (c) [{k}{suffix}]", tot, thr, "<=", n=n,
                 expected_if_each_tested_tumour_were_a_loss_with_p_0_05=round(e0, 3))
    pw = out.get("power", {})
    if "default" in pw:
        _chk(res, "P1 default cohort stress tumours detected", pw["default"]["stress_tumours_detected"]["p"], 0.95, ">=")
        _chk(res, "P2 default cohort stress P(crit c)", pw["default"]["stress_p_crit_c"]["p"], 1.0, ">=")
    base_inj = (baseline or {}).get("power", {}).get("inject", {})
    for regime, v in pw.get("inject", {}).items():
        for mp_ in ("cycle", "EMT"):
            for metric, lab in (("injected_tumours_detected", "detected per tumour"), ("p_crit_c", "P(crit c)")):
                cand = v[mp_][metric]
                b = base_inj.get(regime, {}).get(mp_, {}).get(metric)
                name = f"P3 injected {mp_} {lab} [{regime}]"
                if b is None or b.get("p") is None or cand.get("p") is None:
                    res[name] = {"value": cand.get("p"), "pass": None, "note": "informational without --baseline"}
                    continue
                se = np.sqrt((cand.get("se") or 0) ** 2 + (b.get("se") or 0) ** 2)
                tol = max(2 * se, 0.05)
                res[name] = {"value": cand["p"], "baseline": b["p"], "threshold": f">= baseline - {tol:.3f}",
                             "pass": bool(cand["p"] >= b["p"] - tol)}
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pkg_dir")
    ap.add_argument("out_json")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-null-sets", type=int, default=None,
                    help="override protocol_shift_kappa's n_null_sets (toy) and the synthetic number of random sets")
    ap.add_argument("--parts", default="a,b,c,d", help="subset of a,b,c,d")
    ap.add_argument("--toy-seeds", type=int, default=None, help="seeds of the toy null (default 12, 4 with --quick)")
    ap.add_argument("--variant-seeds", type=int, default=None,
                    help="seeds of each toy-null variant (unequal depth, small pieces; default 12, 2 with --quick)")
    ap.add_argument("--inject-frac", default="0.2,0.25",
                    help="fractions of program genes permuted in the injected pieces (quick: the last only)")
    ap.add_argument("--baseline", default=None,
                    help="OUT_JSON of an earlier run (e.g. reference/c_v1.2.json), for the power comparisons (P3)")
    args = ap.parse_args()
    t0 = time.time()
    OPTS.update(n_null_sets=args.n_null_sets, quick=args.quick)
    protocol, sim, data, report = pkg(args.pkg_dir)
    parts = set(args.parts.split(","))
    n_seeds = args.toy_seeds or (4 if args.quick else 12)
    jobs = []
    if parts & {"a", "c"}:
        jobs += [dict(j, kind="synthetic") for j in synthetic_jobs(args.quick) if j["part"] in parts]
    if "b" in parts:
        for mag in ("cal", "first"):
            for s in range(n_seeds):
                jobs.append(dict(kind="toy_null", seed=s, mag=mag, cost=40.0, name=f"toy_null {mag} s{s}"))
        for variant in VARIANTS:
            for s in range(args.variant_seeds or (2 if args.quick else 12)):
                jobs.append(dict(kind="toy_null", seed=s, mag="cal", variant=variant,
                                 cost=15.0 if variant in ("n_mal_60", "test_cut_50") else 40.0,
                                 name=f"toy_null cal+{variant} s{s}"))
    if "d" in parts:
        for s in range(4 if not args.quick else 2):
            jobs.append(dict(kind="default", seed=s, cost=10.0, name=f"default s{s}"))
        fracs = [float(x) for x in args.inject_frac.split(",")]
        for frac in (fracs[-1:] if args.quick else fracs):
            for s in range(6 if not args.quick else 2):
                jobs.append(dict(kind="inject", seed=s, frac=frac, cost=10.0, name=f"inject f{frac} s{s}"))
    print(f"[evaluate_c] {len(jobs)} jobs on {args.workers} workers; package {common.resolve_pkg(args.pkg_dir)[1]}",
          flush=True)
    res = common.pool_map("job_dispatch", "__main__", args.pkg_dir, jobs, args.workers)
    by = {}
    errors = []
    for r in res:
        if not r["ok"]:
            errors.append({"job": r["job"].get("name"), "error": r["error"]})
            continue
        by.setdefault(r["job"]["kind"], []).append(dict(r["result"], part=r["job"].get("part")))
    out = {"harness": "evaluate_c", "pkg_dir": str(common.resolve_pkg(args.pkg_dir)[1]),
           "modules": common.check_modules(args.pkg_dir), "quick": args.quick,
           "settings": {"n_null_sets_toy": args.n_null_sets or default_n_null_sets(protocol),
                        "toy_seeds": n_seeds, "variant_seeds": args.variant_seeds or (2 if args.quick else 12),
                        "variants": list(VARIANTS), "inject_frac": args.inject_frac, "frozen_rule": frozen_rule(report),
                        "kappa_min_effect": float(getattr(protocol, "KAPPA_MIN_EFFECT", np.nan)),
                        "kappa_active": float(getattr(protocol, "KAPPA_ACTIVE", np.nan)),
                        "synthetic_kwargs": synthetic_kwargs(protocol)},
           "errors": errors}
    if "synthetic" in by:
        out["synthetic"] = summarise_synthetic(by["synthetic"])
    if "toy_null" in by:
        out["toy_null"] = summarise_toy_null(by["toy_null"])
    out["power"] = {}
    if "default" in by:
        out["power"]["default"] = summarise_default(by["default"])
    if "inject" in by:
        out["power"]["inject"] = summarise_inject(by["inject"])
    out["headline"] = headline(out)
    baseline = None
    if args.baseline:
        import json
        baseline = json.loads(Path(args.baseline).read_text())
    out["checks"] = checks(out, baseline)
    out["runtime_seconds"] = round(time.time() - t0, 1)
    out["job_seconds"] = {r["job"].get("name"): round(r["seconds"], 1) for r in res}
    out["raw"] = {k: v for k, v in by.items() if k in ("toy_null", "inject", "default")}
    common.dump(out, args.out_json)
    print(f"[evaluate_c] wrote {args.out_json} in {out['runtime_seconds']:.0f}s; {len(errors)} errors", flush=True)


def job_dispatch(pkg_dir: str, job: dict):
    return {"synthetic": job_synthetic, "toy_null": job_toy_null, "default": job_default,
            "inject": job_inject}[job["kind"]](pkg_dir, job)


if __name__ == "__main__":
    main()
