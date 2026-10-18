"""Evaluation harness for criterion (a)'s kappa_retained guard of a candidate mpaudit package.

    VENV/bin/python evaluate_a.py PKG_DIR OUT_JSON [--quick] [--workers 4] [--decoys 2]

PKG_DIR is a copy of src/mpaudit (possibly modified) or a directory containing one. Only the public
decoy API is called: ``make_decoy`` and ``kappa_retained`` (no NMF, so decoy recruitment itself is not
evaluated). See README.md for the metrics and pass criteria.
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402  (sets the thread env vars before numpy)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ARMS = {  # decoy_recruitment uses each arm's columns sorted and de-duplicated
    "technical": ("contamination_est", "doublet_est", "intronic_fraction", "pct_ribo"),  # registered (demo) arm
    "pct_ribo": ("pct_ribo",),
    "true_doublet": ("contamination_est", "doublet_score", "intronic_fraction", "pct_ribo"),  # demo --true-doublet
    "dosage": ("clone_est",),
}
TECH_SETS = ("abundance", "long_genes", "contamination")
BIO_SETS = ("cycle", "hypoxia", "EMT", "IFN")
MIXED = ("EMT", "IFN", "CNA_7q")
COHORT_ARMS = {"default": ("technical", "pct_ribo", "true_doublet", "dosage"),
               "null": ("technical", "pct_ribo", "true_doublet"),
               "minority": ("technical", "pct_ribo")}
COHORT_SETS = {"default": (*TECH_SETS, *BIO_SETS, "CNA_7q"), "null": BIO_SETS, "minority": ("cycle", "EMT", "IFN")}
DECOY_SEED = 1000  # decoy_recruitment's default seed: decoy d uses DECOY_SEED + d and the floor uses DECOY_SEED
OPTS: dict = {}


def pkg(pkg_dir: str):
    if "mpaudit.decoy" not in sys.modules:
        common.load_mpaudit(pkg_dir)
    import mpaudit.data as data
    import mpaudit.decoy as decoy
    import mpaudit.report as report
    import mpaudit.sim as sim
    common.check_modules(pkg_dir)
    return decoy, sim, data, report


def simulate(sim, data, cohort: str, seed: int):
    if cohort == "default":
        adata, truth = sim.simulate_cohort(seed=seed)
    elif cohort == "null":
        adata, truth = sim.simulate_cohort(seed=seed, dissociation_stress=False, tme_contamination=False, cna=False,
                                           technical_axes=False)
    elif cohort == "minority":
        adata, truth = sim.simulate_cohort(n_tumours=16, seed=seed)
    else:
        raise ValueError(cohort)
    data.add_qc_obs(adata, counts_layer="counts")
    return adata, truth


def carriers(truth, mp: str) -> set | None:
    return {"EMT": set(truth["tumours_emt"]), "IFN": set(truth["tumours_ifn"]),
            "CNA_7q": set(truth["tumours_cna"])}.get(mp)


def minority_tumours(truth, adata, mp: str) -> tuple[list[str], list[str]]:
    """(kept carriers, non-carriers): the first two carrier tumours and every non-carrier tumour."""
    tum = sorted(adata.obs["tumour"].astype(str).unique())
    car = sorted(carriers(truth, mp))
    return car[:2], [t for t in tum if t not in car]


def kr_call(decoy_mod, adata, dec, keys, sets) -> pd.DataFrame:
    """Public kappa_retained (registered settings = the package defaults; floor seed = decoy_recruitment's)."""
    return decoy_mod.kappa_retained(adata, dec, keys, sets, counts_layer="counts", seed=DECOY_SEED)


def job(pkg_dir: str, j: dict) -> dict:
    decoy_mod, sim, data, _ = pkg(pkg_dir)
    cohort, seed, arm = j["cohort"], j["seed"], j["arm"]
    adata, truth = simulate(sim, data, cohort, seed)
    keys = data.Keys(sample="sample", malignant="malignant", tumour="tumour", study="study")
    names = COHORT_SETS[cohort]
    sets = {k: list(truth["sets"][k]) for k in names if k in truth["sets"]}
    mal = adata.obs["malignant"].to_numpy(bool)
    sample = adata.obs["sample"].astype(str).to_numpy()
    tumour = adata.obs["tumour"].astype(str).to_numpy()
    samples = list(dict.fromkeys(sample[mal]))
    tum_of = {s: tumour[np.flatnonzero(sample == s)[0]] for s in samples}
    cols = tuple(sorted(set(ARMS[arm])))
    out = {"cohort": cohort, "seed": seed, "arm": arm, "columns": list(cols), "decoys": []}
    minority = {}
    if cohort == "minority":
        for mp_ in ("EMT", "IFN"):
            kept, non = minority_tumours(truth, adata, mp_)
            minority[mp_] = {"carriers_kept": kept, "noncarriers": non, "n_tumours": len(kept) + len(non)}
    for d in range(j["n_decoys"]):
        dec = decoy_mod.make_decoy(adata, keys, counts_layer="counts", strata=cols, seed=DECOY_SEED + d)
        whole = kr_call(decoy_mod, adata, dec, keys, sets)
        rec = {"decoy": d, "whole": {r.mp: {"kappa_retained": None if pd.isna(r.kappa_retained) else float(r.kappa_retained),
                                             "n_samples": int(r.n_samples) if "n_samples" in whole else None}
                                     for r in whole.itertuples()}}
        # per-sample decomposition (each malignant sample alone; same full set list, so pool exclusions match)
        per_sample = {}
        if any(m in sets for m in MIXED):
            for s in samples:
                cells = mal & (sample == s)
                a_s, d_s = adata[cells].copy(), dec[cells].copy()
                k = kr_call(decoy_mod, a_s, d_s, keys, sets)
                per_sample[s] = {r.mp: (None if pd.isna(r.kappa_retained) else float(r.kappa_retained),
                                        int(r.n_samples) if "n_samples" in k else int(pd.notna(r.kappa_retained)))
                                 for r in k.itertuples()}
        rec["per_sample"] = per_sample
        # minority-MP cohorts, as coded (the whole call on the restricted cohort)
        rec["minority"] = {}
        for mp_, info in minority.items():
            keep = np.isin(tumour, info["carriers_kept"] + info["noncarriers"])
            k = kr_call(decoy_mod, adata[keep].copy(), dec[keep].copy(), keys, sets)
            row = k[k["mp"] == mp_].iloc[0]
            rec["minority"][mp_] = {"kappa_retained": None if pd.isna(row.kappa_retained) else float(row.kappa_retained),
                                    "n_samples": int(row.n_samples) if "n_samples" in k else None}
        out["decoys"].append(rec)
    out["tumour_of_sample"] = tum_of
    out["carriers"] = {m: sorted(carriers(truth, m)) for m in MIXED if m in sets}
    out["minority_design"] = minority
    return out


# --------------------------------------------------------------------------- aggregation

def _mean_se(x: list[float]) -> dict:
    x = [v for v in x if v is not None and np.isfinite(v)]
    if not x:
        return {"mean": None, "se": None, "min": None, "max": None, "n": 0}
    a = np.array(x, float)
    return {"mean": round(float(a.mean()), 4), "se": round(float(a.std(ddof=1) / np.sqrt(a.size)), 4) if a.size > 1 else None,
            "min": round(float(a.min()), 4), "max": round(float(a.max()), 4), "n": int(a.size)}


def decompose(r: dict, mp_: str) -> dict:
    """Per decoy: ratios of entering carrier / non-carrier samples (per-sample calls)."""
    car = set(r["carriers"].get(mp_, []))
    res = []
    for rec in r["decoys"]:
        c, n = {}, {}
        for s, v in rec["per_sample"].items():
            if mp_ not in v:
                continue
            ratio, entered = v[mp_]
            if entered and ratio is not None:
                (c if r["tumour_of_sample"][s] in car else n)[s] = ratio
        res.append((c, n))
    return res


def summarise(results: list[dict], decoys: int) -> dict:
    """Per (cohort, arm) and set; every per-seed list is in the order of ``seeds`` (ascending)."""
    out = {}
    for (cohort, arm), rs in _group(results, ("cohort", "arm")).items():
        rs = sorted(rs, key=lambda r: r["seed"])
        sets = sorted({m for r in rs for rec in r["decoys"] for m in rec["whole"]})
        res = {}
        for mp_ in sets:
            per_seed, n_samp = [], []
            for r in rs:
                v = [rec["whole"][mp_]["kappa_retained"] for rec in r["decoys"]]
                v = [x for x in v if x is not None]
                per_seed.append(float(np.mean(v)) if v else None)  # decoy_recruitment: mean over decoys
                n_samp.append(r["decoys"][0]["whole"][mp_]["n_samples"])
            entry = {"kappa_retained": _mean_se(per_seed), "seeds": [r["seed"] for r in rs], "per_seed": per_seed,
                     "n_samples_whole": n_samp}
            if mp_ in MIXED and any(r["decoys"][0]["per_sample"] for r in rs):
                n_car, n_non, all_med, car_med, non_ratios = [], [], [], [], []
                for r in rs:
                    dd = decompose(r, mp_)
                    n_car.append(len(dd[0][0]))
                    n_non.append(len(dd[0][1]))
                    am = [float(np.median(list(c.values()) + list(n.values()))) for c, n in dd if c or n]
                    cm = [float(np.median(list(c.values()))) for c, n in dd if c]
                    all_med.append(float(np.mean(am)) if am else None)
                    car_med.append(float(np.mean(cm)) if cm else None)
                    non_ratios += [round(x, 3) for c, n in dd for x in n.values()]
                shift = [a - c for a, c in zip(all_med, car_med) if a is not None and c is not None]
                car_ratios = [x for r in rs for c, n in decompose(r, mp_) for x in c.values()]
                entry["decomposition"] = {
                    "n_carrier_samples_entering": n_car, "n_noncarrier_samples_entering": n_non,
                    "n_noncarrier_samples_total": [sum(1 for s, t in r["tumour_of_sample"].items()
                                                       if t not in set(r["carriers"].get(mp_, [])))
                                                   for r in rs],
                    "kappa_retained_all_entering": _mean_se(all_med),
                    "kappa_retained_carriers_only": _mean_se(car_med),
                    "noncarrier_shift_of_median": _mean_se(shift),
                    "noncarrier_ratios": non_ratios,
                    "carrier_ratio_median": round(float(np.median(car_ratios)), 4) if car_ratios else None,
                }
            res[mp_] = entry
        if cohort == "minority":
            res["_minority"] = summarise_minority(rs)
        out[f"{cohort}/{arm}"] = res
    return out


def summarise_minority(rs: list[dict]) -> dict:
    rs = sorted(rs, key=lambda r: r["seed"])
    out = {}
    for mp_ in ("EMT", "IFN"):
        coded, coded_n, pair_all, pair_car, pair_max = [], [], [], [], []
        for r in rs:
            v = [rec["minority"][mp_]["kappa_retained"] for rec in r["decoys"]]
            v = [x for x in v if x is not None]
            coded.append(float(np.mean(v)) if v else None)
            coded_n.append(r["decoys"][0]["minority"][mp_]["n_samples"])
            # every pair of carrier tumours + all non-carriers, from the per-sample decomposition
            car = sorted(r["carriers"][mp_])
            dd = decompose(r, mp_)
            vals, vals_c = [], []
            for pair in itertools.combinations(car, 2):
                meds, meds_c = [], []
                for c, n in dd:
                    cc = [x for s, x in c.items() if r["tumour_of_sample"][s] in pair]
                    if cc or n:
                        meds.append(float(np.median(cc + list(n.values()))))
                    if cc:
                        meds_c.append(float(np.median(cc)))
                if meds:
                    vals.append(float(np.mean(meds)))
                if meds_c:
                    vals_c.append(float(np.mean(meds_c)))
            pair_all.append(float(np.mean(vals)) if vals else None)
            pair_max.append(float(np.max(vals)) if vals else None)
            pair_car.append(float(np.mean(vals_c)) if vals_c else None)
        out[mp_] = {"design_per_seed": {r["seed"]: r["minority_design"][mp_] for r in rs},
                    "kappa_retained_as_coded": _mean_se(coded), "seeds": [r["seed"] for r in rs],
                    "per_seed_as_coded": coded,
                    "n_samples_as_coded": coded_n,
                    "pairs_mean_all_entering": _mean_se(pair_all), "pairs_max_all_entering": _mean_se(pair_max),
                    "pairs_mean_carriers_only": _mean_se(pair_car),
                    "noncarrier_shift_pairs_mean": _mean_se([a - c for a, c in zip(pair_all, pair_car)
                                                             if a is not None and c is not None])}
    return out


def _group(results, keys):
    g = {}
    for r in results:
        g.setdefault(tuple(r[k] for k in keys), []).append(r)
    return g


def checks(s: dict, thr: float) -> dict:
    res = {}

    def mx(entry):
        return entry["kappa_retained"]["max"], entry["kappa_retained"]["min"]

    # A1 technical sets under the registered technical arm (default cohort), every seed
    for m in TECH_SETS:
        e = s.get("default/technical", {}).get(m)
        if e:
            res[f"A1 {m} kappa_retained >= {thr} [default/technical] (min over seeds)"] = {
                "value": e["kappa_retained"]["min"], "pass": bool(e["kappa_retained"]["min"] is not None
                                                                   and e["kappa_retained"]["min"] >= thr)}
    # A2 biology sets < thr under the technical-type arms, every cohort, every seed
    for key, v in s.items():
        cohort, arm = key.split("/")
        if arm == "dosage":
            continue
        for m in BIO_SETS:
            if m in v:
                hi = v[m]["kappa_retained"]["max"]
                tag = "A2" if arm in ("technical", "pct_ribo") else "A2-info(known coupled-covariate limit)"
                res[f"{tag} {m} kappa_retained < {thr} [{key}] (max over seeds; margin)"] = {
                    "value": hi, "margin": None if hi is None else round(thr - hi, 4),
                    "pass": bool(hi is not None and hi < thr)}
    # A3 CNA under the dosage arm
    e = s.get("default/dosage", {}).get("CNA_7q")
    if e:
        res[f"A3 CNA_7q kappa_retained >= {thr} [default/dosage] (min over seeds)"] = {
            "value": e["kappa_retained"]["min"], "pass": bool(e["kappa_retained"]["min"] is not None
                                                               and e["kappa_retained"]["min"] >= thr)}
    # A4 non-carrier influence
    for key, v in s.items():
        cohort, arm = key.split("/")
        if arm not in ("technical", "pct_ribo"):
            continue
        for m in ("EMT", "IFN", "CNA_7q"):
            dec = v.get(m, {}).get("decomposition")
            if dec and dec["noncarrier_shift_of_median"]["mean"] is not None:
                sh = dec["noncarrier_shift_of_median"]["max"]
                res[f"A4 |non-carrier shift of the median| <= 0.05 [{key} {m}] (max over seeds)"] = {
                    "value": sh, "n_noncarrier_entering": dec["n_noncarrier_samples_entering"],
                    "pass": bool(abs(sh) <= 0.05)}
        if "_minority" in v:
            for m, mm in v["_minority"].items():
                hi = mm["kappa_retained_as_coded"]["max"]
                res[f"A5 minority {m} kappa_retained < {thr} as coded [{key}] (max over seeds; margin)"] = {
                    "value": hi, "margin": None if hi is None else round(thr - hi, 4),
                    "pass": bool(hi is not None and hi < thr)}
                hi = mm["pairs_max_all_entering"]["max"]
                res[f"A5b minority {m} kappa_retained < {thr}, worst pair of carrier tumours [{key}] (max over seeds)"] = {
                    "value": hi, "margin": None if hi is None else round(thr - hi, 4),
                    "pass": bool(hi is not None and hi < thr)}
                sh = mm["noncarrier_shift_pairs_mean"]["max"]
                res[f"A6 minority {m} non-carrier shift (pairs of carriers) <= 0.05 [{key}] (max over seeds)"] = {
                    "value": sh, "pass": bool(sh is not None and sh <= 0.05)}
    return res


def headline(s: dict) -> dict:
    h = {}
    for key, v in s.items():
        for m, e in v.items():
            if m.startswith("_"):
                continue
            h[f"kappa_retained[{key}][{m}]"] = e["kappa_retained"]["mean"]
            if "decomposition" in e:
                d = e["decomposition"]
                h[f"noncarrier_entering[{key}][{m}]"] = f"{sum(d['n_noncarrier_samples_entering'])}/" \
                                                         f"{sum(d['n_noncarrier_samples_total'])}"
                h[f"noncarrier_shift[{key}][{m}]"] = d["noncarrier_shift_of_median"]["mean"]
        if "_minority" in v:
            for m, mm in v["_minority"].items():
                h[f"minority_as_coded[{key}][{m}]"] = mm["kappa_retained_as_coded"]["mean"]
                h[f"minority_carriers_only_pairs[{key}][{m}]"] = mm["pairs_mean_carriers_only"]["mean"]
                h[f"minority_all_pairs_max[{key}][{m}]"] = mm["pairs_max_all_entering"]["max"]
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pkg_dir")
    ap.add_argument("out_json")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--decoys", type=int, default=2)
    ap.add_argument("--seeds", type=int, default=None, help="default 3, 1 with --quick")
    ap.add_argument("--cohorts", default="default,null,minority")
    args = ap.parse_args()
    t0 = time.time()
    decoy_mod, sim, data, report = pkg(args.pkg_dir)
    rule = report.validate_rule(None) if hasattr(report, "validate_rule") else report.DEFAULT_RULE
    thr = float(rule["min_kappa_retained"])
    n_seeds = args.seeds or (1 if args.quick else 3)
    cost = {"default": 1.0, "null": 0.8, "minority": 1.5}
    jobs = [dict(cohort=c, seed=s, arm=a, n_decoys=args.decoys, cost=cost[c], name=f"{c} s{s} {a}")
            for c in args.cohorts.split(",") for a in COHORT_ARMS[c] for s in range(n_seeds)]
    print(f"[evaluate_a] {len(jobs)} jobs on {args.workers} workers; package {common.resolve_pkg(args.pkg_dir)[1]}",
          flush=True)
    res = common.pool_map("job", "__main__", args.pkg_dir, jobs, args.workers)
    ok = [r["result"] for r in res if r["ok"]]
    errors = [{"job": r["job"].get("name"), "error": r["error"]} for r in res if not r["ok"]]
    s = summarise(ok, args.decoys)
    out = {"harness": "evaluate_a", "pkg_dir": str(common.resolve_pkg(args.pkg_dir)[1]),
           "modules": common.check_modules(args.pkg_dir), "quick": args.quick,
           "settings": {"seeds": n_seeds, "decoys_per_arm": args.decoys, "decoy_seeds": [DECOY_SEED + d for d in range(args.decoys)],
                        "floor_seed": DECOY_SEED, "arms": {k: list(v) for k, v in ARMS.items()},
                        "min_kappa_retained": thr},
           "errors": errors, "summary": s, "headline": headline(s), "checks": checks(s, thr),
           "runtime_seconds": round(time.time() - t0, 1),
           "job_seconds": {r["job"].get("name"): round(r["seconds"], 1) for r in res},
           "raw": ok}
    common.dump(out, args.out_json)
    print(f"[evaluate_a] wrote {args.out_json} in {out['runtime_seconds']:.0f}s; {len(errors)} errors", flush=True)


if __name__ == "__main__":
    main()
