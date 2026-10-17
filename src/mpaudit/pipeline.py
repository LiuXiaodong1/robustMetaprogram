"""Pilot P1 runner: audit a set of MPs (re-derived or published) along every available axis and
apply the frozen GO rule (:mod:`mpaudit.report`, ``P1-GO v1.0``)."""
from __future__ import annotations

import copy
import json
import platform
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from .compare import match_sets
from .data import Keys, add_qc_obs, read_gmt, sample_map, write_gmt
from .decoy import decoy_recruitment
from .gavish import GavishParams, NMFCache, eligible_universe, run_gavish
from .instability import default_perturbations, run_instability
from .positional import (arm_enrichment, assign_arms, centromeres, gene_lengths, gene_positions_from_gtf,
                         load_cytoband)
from .protocol import protocol_shift_kappa
from .report import DEFAULT_RULE, RULE_VERSION, decide, rule_text, write_summary
from .tme import tme_covariation
from .typer import evaluate_typer, typer_features

DEFAULT_CONFIG = {
    "keys": {"sample": "sample", "malignant": "malignant"},
    "counts_layer": None,
    "recurrence_key": None,       # obs column defining independent units (e.g. tumour)
    "mps": "rederive",            # or path to a GMT of published MPs to audit
    "gavish": {},                 # GavishParams overrides
    "n_jobs": 1,
    "cache_dir": None,
    "min_jaccard": 0.3,
    "qc": {"add_obs": True},      # add log_umi, log_genes, pct_mito, pct_ribo to obs if missing
    "instability": {"enabled": True, "seeds": [1, 2, 3], "n_boot": 3, "boot_fraction": 0.8,
                    "leave_one_study_out": True},
    # arms: name -> obs columns the decoys condition on; categorical columns become permutation
    # strata, continuous ones enter the mean model. Technical and dosage arms are typed separately.
    "decoy": {"enabled": True, "n": 5, "method": "residual", "arms": {}, "strata": [], "n_depth_bins": 1,
              "min_kappa": 0.02, "n_null_floor": 20},
    "tme": {"enabled": True, "min_cells": 20, "n_boot": 200, "n_null": 500, "covariates": None,
            "spillover_adjust": True, "null_match": "expression+specificity", "n_markers": 50,
            "min_marker_fold": 8.0},
    "positional": {"enabled": True, "cytoband": None, "gtf": None,
                   "chrom_col": "chromosome", "start_col": "start"},
    "typer": {"enabled": True, "labels": None},   # CSV with columns mp,label (1 = curated artefact)
    "protocol": {"enabled": True, "reference": None, "test": None, "min_pairs": 3, "n_null_sets": 100},
    "decision": dict(DEFAULT_RULE),
    "seed": 0,
}


def merge_config(user: dict | None) -> dict:
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    for k, v in (user or {}).items():
        if isinstance(v, dict) and isinstance(cfg.get(k), dict):
            cfg[k].update(v)
        else:
            cfg[k] = v
    return cfg


def _decoy_arms(dc: dict, log) -> dict[str, tuple[str, ...]]:
    arms = {name: tuple(cols) for name, cols in (dc.get("arms") or {}).items() if cols}
    if not arms and dc.get("strata"):
        log("[p1] decoy.strata given without arms: treated as one 'technical' arm; put CNA clone columns "
            "in a separate 'dosage' arm so clone-driven recruitment is typed as dosage")
        arms = {"technical": tuple(dc["strata"])}
    return arms


def _wide_decoys(per_mp: pd.DataFrame) -> pd.DataFrame:
    """One row per MP: recruitment and kappa_retained per arm (for the typer table)."""
    wide = per_mp.pivot(index="mp", columns="arm", values=["recruitment", "kappa_retained"])
    wide.columns = [f"{stat}_{arm}" for stat, arm in wide.columns]
    return wide.reset_index()


def run_p1(adata, config: dict | None = None, out_dir: str | Path = "results", log=print) -> dict:
    cfg = merge_config(config)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    keys = Keys.from_dict(cfg["keys"])
    params = GavishParams.from_dict(cfg["gavish"])
    cache = NMFCache(cfg["cache_dir"])
    layer, rk, n_jobs, seed = cfg["counts_layer"], cfg["recurrence_key"], int(cfg["n_jobs"]), int(cfg["seed"])
    tables: dict[str, pd.DataFrame] = {}
    if cfg["qc"].get("add_obs", True):
        add_qc_obs(adata, counts_layer=layer)

    reference = None
    if cfg["mps"] == "rederive" or cfg["instability"]["enabled"]:
        log("[p1] running the reference 3CA-style pipeline")
        reference = run_gavish(adata, keys, params, counts_layer=layer, cache=cache, recurrence_key=rk, n_jobs=n_jobs)
        write_gmt(reference.mp_sets(), out / "reference_mps.gmt")
        reference.programs_table().to_csv(out / "programs.csv", index=False)
        log(f"[p1] reference: {len(reference.robust)} robust programs, {len(reference.mps)} MPs; "
            f"{reference.n_capped_fits}/{reference.n_fits} NMF fits stopped at max_iter")
    mp_sets = reference.mp_sets() if cfg["mps"] == "rederive" else read_gmt(cfg["mps"])
    if not mp_sets:
        log("[p1] no MPs to audit")
    log(f"[p1] auditing {len(mp_sets)} MPs")

    if cfg["instability"]["enabled"] and reference is not None:
        ic = cfg["instability"]
        study_samples = None
        if ic.get("leave_one_study_out") and keys.study:
            sm = sample_map(adata, keys, keys.study)
            study_samples = {}
            for s, st in sm.items():
                study_samples.setdefault(st, []).append(s)
            if len(study_samples) < 2:
                study_samples = None
        perts = default_perturbations(params, study_samples=study_samples, seeds=tuple(ic["seeds"]),
                                      n_boot=int(ic["n_boot"]), boot_fraction=float(ic["boot_fraction"]))
        inst = run_instability(adata, keys, params, perts, counts_layer=layer, cache=cache, recurrence_key=rk,
                               min_jaccard=cfg["min_jaccard"], n_jobs=n_jobs, reference=reference, log=log)
        tables["instability_runs"] = inst.summary
        tables["gene_stability"] = inst.genes
        if cfg["mps"] == "rederive":
            tables["instability_per_mp"] = inst.per_mp
        else:
            # audited MPs come from a GMT: report the recovery of the matching re-derived MP
            m = match_sets(mp_sets, reference.mp_sets(), cfg["min_jaccard"])
            per = inst.per_mp.set_index("mp")
            m["recovery_all"] = [per.loc[a, "recovery_all"] if a else np.nan for a in m["alt"]]
            tables["instability_rederived_per_mp"] = inst.per_mp
            tables["instability_per_mp"] = m.rename(columns={"ref": "mp", "alt": "rederived_mp"})[
                ["mp", "rederived_mp", "jaccard", "recovery_all"]]

    dc = cfg["decoy"]
    arms = _decoy_arms(dc, log)
    if dc["enabled"] and mp_sets and arms:
        per_mp, per_decoy, _ = decoy_recruitment(
            adata, keys, params, mp_sets, arms=arms, n_decoys=int(dc["n"]), n_depth_bins=int(dc["n_depth_bins"]),
            method=dc["method"], counts_layer=layer, cache=cache, recurrence_key=rk, min_jaccard=cfg["min_jaccard"],
            min_kappa=float(dc["min_kappa"]), n_null_floor=int(dc["n_null_floor"]), seed=1000 + seed,
            n_jobs=n_jobs, log=log)
        tables["decoy_recruitment"] = per_mp
        tables["decoy_runs"] = per_decoy
    elif dc["enabled"] and mp_sets:
        log("[p1] decoys skipped: no decoy arms configured")

    if cfg["tme"]["enabled"] and keys.cell_type and mp_sets:
        tc = cfg["tme"]
        cov = tuple(tc["covariates"]) if tc.get("covariates") is not None else None
        per_set, per_sample = tme_covariation(
            adata, keys, mp_sets, counts_layer=layer, min_cells=int(tc["min_cells"]), covariates=cov,
            spillover_adjust=bool(tc["spillover_adjust"]), null_match=tc["null_match"], n_boot=int(tc["n_boot"]),
            n_null=int(tc["n_null"]), n_markers=int(tc["n_markers"]), min_marker_fold=float(tc["min_marker_fold"]),
            seed=seed)
        tables["tme_covariation"] = per_set
        tables["tme_per_sample"] = per_sample
        log("[p1] TME co-variation done")

    pc = cfg["positional"]
    gtf_positions = gene_positions_from_gtf(pc["gtf"]) if pc.get("gtf") else None
    if pc["enabled"] and mp_sets:
        positions = gtf_positions
        if positions is None and pc["chrom_col"] in adata.var and pc["start_col"] in adata.var:
            positions = adata.var[[pc["chrom_col"], pc["start_col"]]].rename(
                columns={pc["chrom_col"]: "chromosome", pc["start_col"]: "start"})
        cyto = pc.get("cytoband")
        if positions is not None and cyto is not None:
            cyto_df = cyto if isinstance(cyto, pd.DataFrame) else load_cytoband(cyto)
            arms_of_gene = assign_arms(positions, centromeres(cyto_df))
            # background: genes that could enter an MP (union of the per-sample NMF universes)
            universe = eligible_universe(adata, keys, params, counts_layer=layer)
            a = arm_enrichment(mp_sets, arms_of_gene, universe)
            b = arm_enrichment(mp_sets, arms_of_gene, universe, exclude_known=True)
            a["q_excluding_known_clusters"] = b["q"].to_numpy()
            tables["positional"] = a
            log(f"[p1] positional enrichment done (background: {len(universe)} eligible genes)")
        else:
            log("[p1] positional axis skipped (needs gene positions and a cytoband file)")

    if cfg["typer"]["enabled"] and mp_sets:
        lengths = gene_lengths(gtf_positions) if gtf_positions is not None else None
        feats = typer_features(adata, keys, mp_sets, counts_layer=layer, gene_length=lengths, seed=seed)
        if "tme_covariation" in tables:
            feats = feats.merge(tables["tme_covariation"][["mp", "partial_r2"]], on="mp", how="left")
        if "decoy_recruitment" in tables:
            feats = feats.merge(_wide_decoys(tables["decoy_recruitment"]), on="mp", how="left")
        if "positional" in tables:
            pos = tables["positional"][["mp", "q"]].copy()
            pos["neg_log10_positional_q"] = -np.log10(pos.pop("q").clip(lower=1e-300))
            feats = feats.merge(pos, on="mp", how="left")
        tables["typer_features"] = feats
        if cfg["typer"].get("labels"):
            lab = pd.read_csv(cfg["typer"]["labels"]).set_index("mp")["label"]
            tables["typer_eval"] = evaluate_typer(feats, lab)

    prc = cfg["protocol"]
    if prc["enabled"] and keys.protocol and prc.get("reference") and prc.get("test") and mp_sets:
        per_mp, per_tumour, levels = protocol_shift_kappa(
            adata, keys, mp_sets, prc["reference"], prc["test"], counts_layer=layer,
            min_pairs=int(prc["min_pairs"]), n_null_sets=int(prc["n_null_sets"]), seed=seed)
        tables["protocol_shift"] = per_mp
        tables["protocol_tumours"] = per_tumour
        tables["protocol_levels"] = levels
        log("[p1] protocol shift done")

    labels = None
    if cfg["typer"].get("labels"):
        labels = pd.read_csv(cfg["typer"]["labels"]).set_index("mp")["label"]
    decision, flags = decide(mp_sets, tables, cfg["decision"], curated_artefacts=labels)
    tables["mp_flags"] = flags
    for name, df in tables.items():
        df.to_csv(out / f"{name}.csv", index=False)
    write_summary(out / "summary.md", decision, flags, tables, cfg)
    versions = {}
    for pkg in ("mpaudit", "numpy", "scipy", "pandas", "anndata"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = "unknown"

    def plain(v):
        if isinstance(v, pd.DataFrame):
            return "<DataFrame>"
        if isinstance(v, dict):
            return {k: plain(x) for k, x in v.items()}
        return v

    (out / "run_config.json").write_text(json.dumps(
        {"config": plain(cfg), "params": vars(params), "rule_version": RULE_VERSION,
         "rule": rule_text(cfg["decision"]), "versions": versions, "python": platform.python_version(),
         "decision": decision}, indent=2, default=str))
    log(f"[p1] decision: {decision['verdict']} ({decision['n_counted']} counted MPs; {RULE_VERSION}); "
        f"see {out / 'summary.md'}")
    return {"decision": decision, "tables": tables, "reference": reference, "mp_sets": mp_sets}
