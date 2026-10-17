"""Pilot P1 runner: audit a set of MPs (re-derived or published) along every available axis and
apply the frozen GO rule (:mod:`mpaudit.report`, ``P1-GO v1.1``).

:func:`validate_config` checks the whole configuration before any NMF runs: the ``decision`` block
(:func:`mpaudit.report.validate_rule`), the decoy arms against the rule's arm lists, the label and
gene-set files, names that cannot match re-derived MPs, and the GavishParams/keys combination.

Curated-artefact labels (``typer.labels``, a CSV with columns ``mp`` and ``label``, 1 = curated
artefact) and expected-TME-shared names are matched to the audited MPs by name. Re-derived MPs
(``mps: rederive``) are named ``MP1``, ``MP2``, ... by the run, in an order that no config written
in advance can know, so for them: expected-TME-shared biology must be given as gene sets
(``decision.expected_tme_shared_gmt`` / ``expected_tme_shared_sets``; names raise), and labels are
carried over by gene overlap from ``typer.labelled_gmt``, a GMT of the labelled (published) MPs
(:func:`curated_labels`). Labels or names that match no audited MP raise or warn.
"""
from __future__ import annotations

import copy
import json
import platform
import re
import warnings
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from .compare import jaccard_matrix, match_sets
from .data import QC_COLUMNS, Keys, add_qc_obs, read_gmt, sample_map, write_gmt
from .decoy import decoy_recruitment
from .gavish import GavishParams, NMFCache, check_keys, eligible_universe, run_gavish
from .instability import default_perturbations, run_instability
from .positional import (arm_enrichment, assign_arms, centromeres, gene_lengths, gene_positions_from_gtf,
                         load_cytoband)
from .protocol import KAPPA_MIN_EFFECT, protocol_shift_kappa
from .report import DEFAULT_RULE, RHO_SENSITIVITY, RULE_VERSION, decide, rule_text, validate_rule, write_summary
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
    # strata, continuous ones enter the mean model. Every arm must be named in
    # decision.technical_arms or decision.dosage_arms. "strata" (one unnamed arm) is no longer accepted.
    "decoy": {"enabled": True, "n": 5, "method": "residual", "arms": {}, "strata": [], "n_depth_bins": 1,
              "min_kappa": 0.02, "n_null_floor": 20},
    "tme": {"enabled": True, "min_cells": 20, "n_boot": 200, "n_null": 500, "covariates": None,
            "spillover_adjust": True, "null_match": "expression+specificity", "n_markers": 50,
            "min_marker_fold": 8.0, "rho_sensitivity": RHO_SENSITIVITY},
    "positional": {"enabled": True, "cytoband": None, "gtf": None,
                   "chrom_col": "chromosome", "start_col": "start"},
    # labels: CSV with columns mp,label (1 = curated artefact); labelled_gmt: gene sets of the
    # labelled MPs, needed to carry the labels onto re-derived MPs
    "typer": {"enabled": True, "labels": None, "labelled_gmt": None},
    "protocol": {"enabled": True, "reference": None, "test": None, "min_pairs": 3, "n_null_sets": 100,
                 "kappa_min_effect": KAPPA_MIN_EFFECT},
    "decision": dict(DEFAULT_RULE),
    "seed": 0,
}
MIN_POSITION_COVERAGE = 0.5  # warn when gene positions or lengths cover fewer var_names than this
_REDERIVED_NAME = re.compile(r"MP\d+")  # names mpaudit.gavish gives re-derived MPs


def merge_config(user: dict | None) -> dict:
    cfg = copy.deepcopy(DEFAULT_CONFIG)
    for k, v in (user or {}).items():
        if isinstance(v, dict) and isinstance(cfg.get(k), dict):
            cfg[k].update(v)
        else:
            cfg[k] = v
    return cfg


def _warn(msg: str, log) -> None:
    log(f"[p1] warning: {msg}")
    warnings.warn(msg, UserWarning, stacklevel=3)


def _decoy_arms(dc: dict) -> dict[str, tuple[str, ...]]:
    """Arm name -> obs columns. Raises for the legacy ``decoy.strata`` and malformed arms."""
    if dc.get("strata"):
        raise ValueError(
            "decoy.strata is no longer accepted: it made every column, CNA clone included, one technical arm, "
            "so clone-driven recruitment counted as technical nuisance. Use decoy.arms, e.g. "
            "{technical: [doublet_score, contamination, pct_ribo, intronic_fraction], dosage: [clone]}, with the "
            "arm names listed in decision.technical_arms / decision.dosage_arms")
    raw = dc.get("arms") or {}
    if not isinstance(raw, dict):
        raise ValueError(f"decoy.arms must map an arm name to obs columns, got {raw!r}")
    arms = {}
    for name, cols in raw.items():
        if isinstance(cols, str):
            cols = [cols]
        if cols is None:
            cols = []
        if not isinstance(cols, (list, tuple)) or not all(isinstance(c, str) and c for c in cols):
            raise ValueError(f"decoy.arms.{name} must be a list of obs columns, got {cols!r}")
        if cols:
            arms[str(name)] = tuple(cols)
    return arms


def _read_labels(path) -> pd.Series:
    df = pd.read_csv(path)
    if not {"mp", "label"} <= set(df.columns):
        raise ValueError(f"typer.labels {path}: needs columns 'mp' and 'label', has {list(df.columns)}")
    dup = df["mp"][df["mp"].duplicated()].unique().tolist()
    if dup:
        raise ValueError(f"typer.labels {path}: duplicated MP name(s) {dup}")
    lab = pd.to_numeric(df["label"], errors="coerce")
    if lab.isna().any() or not lab.isin([0, 1]).all():
        raise ValueError(f"typer.labels {path}: labels must be 0 or 1")
    return pd.Series(lab.astype(int).to_numpy(), index=df["mp"].astype(str), name="label")


def _optional_number(section: str, key: str, v, *, allow_none: bool, minimum: float, strict: bool):
    if v is None and allow_none:
        return None
    try:
        if isinstance(v, bool):
            raise TypeError
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{section}.{key} must be a number{' or null' if allow_none else ''}, got {v!r}") from None
    if not np.isfinite(x) or x < minimum or (strict and x == minimum):
        raise ValueError(f"{section}.{key} must be {'>' if strict else '>='} {minimum:g}, got {v!r}")
    return x


def validate_config(cfg: dict, adata=None, log=print) -> dict:
    """Check a merged config (:func:`merge_config`) before anything runs; returns it normalised.

    Raises ``ValueError`` for: an invalid ``decision`` block (:func:`mpaudit.report.validate_rule`);
    ``decoy.strata`` (legacy; use ``decoy.arms``); a decoy arm named in neither
    ``decision.technical_arms`` nor ``decision.dosage_arms``; decoy arm columns missing from
    ``adata.obs`` (QC columns that ``qc.add_obs`` will add count as present);
    ``decision.expected_tme_shared`` names with ``mps: rederive``; a missing ``mps`` GMT,
    ``typer.labels`` or ``typer.labelled_gmt`` file, a ``labelled_gmt`` without ``labels``, or,
    with ``mps: rederive`` and no ``labelled_gmt``, labels none of which can name a re-derived MP;
    a bad ``tme.rho_sensitivity`` or ``protocol.kappa_min_effect``; and what
    :func:`mpaudit.gavish.check_keys` rejects. Warns when no configured decoy arm is a dosage arm.
    Name checks against the audited MPs themselves are in :func:`check_names`.
    """
    cfg["decision"] = rule = validate_rule(cfg.get("decision"))
    rederive = cfg["mps"] == "rederive"
    if not rederive and not Path(str(cfg["mps"])).is_file():
        raise ValueError(f"mps must be 'rederive' or a GMT file, got {cfg['mps']!r}")
    keys = Keys.from_dict(cfg["keys"])
    check_keys(keys, GavishParams.from_dict(cfg["gavish"]), adata, recurrence_key=cfg["recurrence_key"])

    dc = cfg["decoy"]
    arms = _decoy_arms(dc)
    if dc.get("enabled", True):
        known = set(rule["technical_arms"]) | set(rule["dosage_arms"])
        stray = sorted(set(arms) - known)
        if stray:
            raise ValueError(f"decoy arm(s) {stray} are in neither decision.technical_arms "
                             f"{rule['technical_arms']} nor decision.dosage_arms {rule['dosage_arms']}: the rule "
                             "would ignore them (a CNA clone arm not listed as dosage would not type dosage)")
        if arms and not set(arms) & set(rule["dosage_arms"]):
            _warn(f"no decoy arm is a dosage arm (decision.dosage_arms {rule['dosage_arms']}): recruitment by CNA "
                  "clones cannot type an MP as dosage, only positional arm enrichment can", log)
        if adata is not None:
            have = set(adata.obs.columns) | (set(QC_COLUMNS) if cfg["qc"].get("add_obs", True) else set())
            missing = sorted({c for cols in arms.values() for c in cols} - have)
            if missing:
                raise ValueError(f"decoy arm column(s) {missing} are not columns of adata.obs")
    dc["arms"] = {k: list(v) for k, v in arms.items()}

    if rederive and rule["expected_tme_shared"]:
        raise ValueError(
            f"decision.expected_tme_shared names {rule['expected_tme_shared']} cannot be pre-declared with mps: "
            "rederive (re-derived MPs are named MP1, MP2, ... by the run, in no fixed order); give the expected "
            "TME-shared biology as gene sets (decision.expected_tme_shared_gmt or expected_tme_shared_sets)")

    tc = cfg["typer"]
    if tc.get("labelled_gmt"):
        if not tc.get("labels"):
            raise ValueError("typer.labelled_gmt needs typer.labels (the labels of its gene sets)")
        if not Path(str(tc["labelled_gmt"])).is_file():
            raise ValueError(f"typer.labelled_gmt {tc['labelled_gmt']!r} is not a file")
    if tc.get("labels"):
        if not Path(str(tc["labels"])).is_file():
            raise ValueError(f"typer.labels {tc['labels']!r} is not a file")
        labels = _read_labels(tc["labels"])
        if tc.get("labelled_gmt") and not set(read_gmt(tc["labelled_gmt"])) & set(labels.index):
            raise ValueError("no set of typer.labelled_gmt has a label in typer.labels (their MP names must agree)")
        if rederive and not tc.get("labelled_gmt") and not any(_REDERIVED_NAME.fullmatch(m) for m in labels.index):
            raise ValueError(
                "typer.labels names cannot match re-derived MPs (mps: rederive names them MP1, MP2, ... in no "
                "fixed order); give typer.labelled_gmt, the gene sets of the labelled MPs, to carry the labels "
                "over by gene overlap")

    cfg["tme"]["rho_sensitivity"] = _optional_number("tme", "rho_sensitivity", cfg["tme"].get("rho_sensitivity"),
                                                     allow_none=True, minimum=0.0, strict=True)
    cfg["protocol"]["kappa_min_effect"] = _optional_number(
        "protocol", "kappa_min_effect", cfg["protocol"].get("kappa_min_effect", KAPPA_MIN_EFFECT),
        allow_none=False, minimum=0.0, strict=False)
    return cfg


def check_names(mp_sets: dict[str, list[str]], cfg: dict, log=print) -> None:
    """Check label and expected-TME-shared names against the audited MPs.

    Raises when ``typer.labels`` is given without ``typer.labelled_gmt`` and none of its names is
    an audited MP. With MPs from a GMT, warns about label or ``expected_tme_shared`` names that
    match no audited MP. With re-derived MPs and labels by name, warns that the names are those of
    this run only.
    """
    names = set(mp_sets)
    tc = cfg["typer"]
    rederive = cfg["mps"] == "rederive"
    if tc.get("labels") and not tc.get("labelled_gmt"):
        labels = _read_labels(tc["labels"])
        hit = sorted(set(labels.index) & names)
        if not hit:
            raise ValueError(f"typer.labels names match none of the {len(names)} audited MPs"
                             + ("; with mps: rederive give typer.labelled_gmt" if rederive else ""))
        miss = sorted(set(labels.index) - names)
        if rederive:
            _warn("typer.labels are matched by name to re-derived MPs, whose names are assigned by the run in no "
                  "fixed order: valid only if written for this exact run (reference_mps.gmt); use typer.labelled_gmt "
                  "to carry labels over by gene overlap", log)
        elif miss:
            _warn(f"{len(miss)} typer.labels name(s) match no audited MP: {miss}", log)
    if not rederive:
        miss = sorted(set(cfg["decision"]["expected_tme_shared"]) - names)
        if miss:
            _warn(f"decision.expected_tme_shared name(s) match no audited MP: {miss}", log)


def curated_labels(mp_sets: dict[str, list[str]], labels: pd.Series, labelled_sets: dict[str, list[str]] | None,
                   *, rederive: bool, min_jaccard: float = 0.3, log=print) -> tuple[pd.Series, pd.DataFrame]:
    """Curated-artefact label per audited MP, and how it was assigned.

    By name for an MP whose name is in ``labels`` (never for re-derived MPs when ``labelled_sets``
    is given). Otherwise, with ``labelled_sets`` (gene sets of the labelled MPs, keyed like
    ``labels``), the MP inherits the label of its best-matching labelled set if their Jaccard index
    is >= ``min_jaccard``. Each MP takes its own best match (not a one-to-one assignment), so a
    curated artefact split into two re-derived MPs labels both. MPs with neither get no label
    (retained). Returns ``(labels indexed by audited MP, table mp/source/labelled_mp/jaccard/label)``.
    """
    use_names = not (rederive and labelled_sets)
    rows = []
    sets = dict(labelled_sets or {})
    unlabelled = sorted(set(sets) - set(labels.index))
    if unlabelled:
        log(f"[p1] {len(unlabelled)} set(s) of typer.labelled_gmt have no label (an MP matching one stays "
            f"unlabelled): {unlabelled}")
    if sets and not set(sets) & set(labels.index):
        raise ValueError("no set of typer.labelled_gmt has a label in typer.labels (names must agree)")
    by_set = [m for m in mp_sets if not (use_names and m in labels.index)]
    J = jaccard_matrix({m: mp_sets[m] for m in by_set}, sets) if sets and by_set else np.zeros((len(by_set), 0))
    best = {m: (list(sets)[int(J[i].argmax())], float(J[i].max())) for i, m in enumerate(by_set) if J.shape[1]}
    for m in mp_sets:
        if use_names and m in labels.index:
            rows.append({"mp": m, "source": "name", "labelled_mp": m, "jaccard": np.nan, "label": labels[m]})
        elif m in best and best[m][1] >= min_jaccard:
            s, j = best[m]
            rows.append({"mp": m, "source": "gene overlap", "labelled_mp": s, "jaccard": j,
                         "label": labels.get(s, np.nan)})
        else:
            rows.append({"mp": m, "source": "", "labelled_mp": best.get(m, ("", np.nan))[0],
                         "jaccard": best.get(m, ("", np.nan))[1], "label": np.nan})
    table = pd.DataFrame(rows, columns=["mp", "source", "labelled_mp", "jaccard", "label"])
    if sets:
        n = int((table["source"] == "gene overlap").sum())
        log(f"[p1] curated labels carried over by gene overlap (Jaccard >= {min_jaccard}) to {n} of "
            f"{len(mp_sets)} MPs; {int(table['label'].eq(1).sum())} labelled curated artefacts")
    return table.set_index("mp")["label"], table


def _coverage(index, var_names) -> float:
    return float(pd.Index(var_names.astype(str)).isin(pd.Index(index).astype(str)).mean()) if len(var_names) else 0.0


def _gtf_positions(adata, pc: dict, log) -> pd.DataFrame | None:
    """Gene positions and lengths from ``positional.gtf``, or None. An empty table counts as missing
    (callers fall back to ``adata.var``); one placing fewer than ``MIN_POSITION_COVERAGE`` of the
    var_names warns."""
    if not pc.get("gtf"):
        return None
    table = gene_positions_from_gtf(pc["gtf"])
    if table.empty:
        _warn(f"positional.gtf {pc['gtf']} gave no gene positions; falling back to adata.var columns", log)
        return None
    cov = _coverage(table.index, adata.var_names)
    if cov < MIN_POSITION_COVERAGE:
        _warn(f"positional.gtf {pc['gtf']} places only {cov:.0%} of var_names (they must match its gene_name "
              "values; Ensembl IDs do not)", log)
    return table


def _var_positions(adata, pc: dict, log) -> pd.DataFrame | None:
    """Positions from the ``adata.var`` columns ``chrom_col``/``start_col``, or None."""
    if pc["chrom_col"] not in adata.var or pc["start_col"] not in adata.var:
        return None
    pos = adata.var[[pc["chrom_col"], pc["start_col"]]].rename(
        columns={pc["chrom_col"]: "chromosome", pc["start_col"]: "start"})
    placed = pos["chromosome"].notna() & pd.to_numeric(pos["start"], errors="coerce").notna()
    cov = float(placed.mean()) if len(pos) else 0.0
    if cov < MIN_POSITION_COVERAGE:
        _warn(f"adata.var {pc['chrom_col']}/{pc['start_col']} place only {cov:.0%} of var_names", log)
    return pos[placed]


def _gene_lengths(adata, gtf_positions: pd.DataFrame | None, log) -> pd.Series | None:
    """Gene lengths for the typer from the GTF table, or None (the typer then takes them from
    ``adata.var``). Warns when the lengths used cover fewer than ``MIN_POSITION_COVERAGE`` of the
    var_names."""
    if gtf_positions is not None:
        lengths, source = gene_lengths(gtf_positions), "the GTF"
    else:
        try:
            lengths, source = gene_lengths(adata.var), "adata.var"
        except ValueError:
            return None
    cov = _coverage(lengths.dropna().index, adata.var_names)
    if cov < MIN_POSITION_COVERAGE:
        _warn(f"gene lengths from {source} cover only {cov:.0%} of var_names", log)
    return lengths if gtf_positions is not None else None


def _wide_decoys(per_mp: pd.DataFrame) -> pd.DataFrame:
    """One row per MP: recruitment and kappa_retained per arm (for the typer table)."""
    wide = per_mp.pivot(index="mp", columns="arm", values=["recruitment", "kappa_retained"])
    wide.columns = [f"{stat}_{arm}" for stat, arm in wide.columns]
    return wide.reset_index()


def run_p1(adata, config: dict | None = None, out_dir: str | Path = "results", log=print) -> dict:
    cfg = validate_config(merge_config(config), adata, log)
    rederive = cfg["mps"] == "rederive"
    mp_sets = None if rederive else read_gmt(cfg["mps"])
    if mp_sets is not None:
        check_names(mp_sets, cfg, log)
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
    if rederive or cfg["instability"]["enabled"]:
        log("[p1] running the reference 3CA-style pipeline")
        reference = run_gavish(adata, keys, params, counts_layer=layer, cache=cache, recurrence_key=rk, n_jobs=n_jobs)
        write_gmt(reference.mp_sets(), out / "reference_mps.gmt")
        reference.programs_table().to_csv(out / "programs.csv", index=False)
        log(f"[p1] reference: {len(reference.robust)} robust programs, {len(reference.mps)} MPs; "
            f"{reference.n_capped_fits}/{reference.n_fits} NMF fits stopped at max_iter, "
            f"{reference.n_failed_fits} by the restart limit")
    if rederive:
        mp_sets = reference.mp_sets()
        check_names(mp_sets, cfg, log)
    if not mp_sets:
        log("[p1] no MPs to audit")
    log(f"[p1] auditing {len(mp_sets)} MPs")

    labels = None
    tc = cfg["typer"]
    if tc.get("labels"):
        labelled_sets = read_gmt(tc["labelled_gmt"]) if tc.get("labelled_gmt") else None
        labels, mapping = curated_labels(mp_sets, _read_labels(tc["labels"]), labelled_sets, rederive=rederive,
                                         min_jaccard=cfg["min_jaccard"], log=log)
        tables["curated_label_mapping"] = mapping

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
        if rederive:
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
    arms = {k: tuple(v) for k, v in dc["arms"].items()}
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
        tc_ = cfg["tme"]
        cov = tuple(tc_["covariates"]) if tc_.get("covariates") is not None else None
        per_set, per_sample = tme_covariation(
            adata, keys, mp_sets, counts_layer=layer, min_cells=int(tc_["min_cells"]), covariates=cov,
            spillover_adjust=bool(tc_["spillover_adjust"]), null_match=tc_["null_match"], n_boot=int(tc_["n_boot"]),
            n_null=int(tc_["n_null"]), n_markers=int(tc_["n_markers"]), min_marker_fold=float(tc_["min_marker_fold"]),
            rho_sensitivity=tc_["rho_sensitivity"], seed=seed)
        tables["tme_covariation"] = per_set
        tables["tme_per_sample"] = per_sample
        unadjusted = per_set.loc[~per_set["spill_adjusted"].astype(bool), "mp"].tolist()
        if unadjusted:
            log(f"[p1] TME co-variation: spillover NOT adjusted for {len(unadjusted)} MP(s); criterion (b) cannot "
                "count them")
        log("[p1] TME co-variation done")

    pc = cfg["positional"]
    need_gtf = mp_sets and (pc["enabled"] or cfg["typer"]["enabled"])
    gtf_positions = _gtf_positions(adata, pc, log) if need_gtf else None
    if pc["enabled"] and mp_sets:
        positions = gtf_positions if gtf_positions is not None else _var_positions(adata, pc, log)
        cyto = pc.get("cytoband")
        if positions is not None and cyto is not None:
            cyto_df = cyto if isinstance(cyto, pd.DataFrame) else load_cytoband(cyto)
            arms_of_gene = assign_arms(positions, centromeres(cyto_df))
            # background: genes that could enter an MP (union of the per-sample NMF universes)
            universe = eligible_universe(adata, keys, params, counts_layer=layer)
            a = arm_enrichment(mp_sets, arms_of_gene, universe)
            b = arm_enrichment(mp_sets, arms_of_gene, universe, exclude_known=True).set_index("mp").reindex(a["mp"])
            # the run without known clusters can have another best arm than the run with them
            a["best_arm_excluding_known_clusters"] = b["best_arm"].to_numpy()
            a["q_excluding_known_clusters"] = b["q"].to_numpy()
            tables["positional"] = a
            log(f"[p1] positional enrichment done (background: {len(universe)} eligible genes)")
        else:
            log("[p1] positional axis skipped (needs gene positions and a cytoband file)")

    if cfg["typer"]["enabled"] and mp_sets:
        lengths = _gene_lengths(adata, gtf_positions, log)
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
        if labels is not None:
            tables["typer_eval"] = evaluate_typer(feats, labels.dropna())

    prc = cfg["protocol"]
    if prc["enabled"] and keys.protocol and prc.get("reference") and prc.get("test") and mp_sets:
        per_mp, per_tumour, levels = protocol_shift_kappa(
            adata, keys, mp_sets, prc["reference"], prc["test"], counts_layer=layer,
            min_pairs=int(prc["min_pairs"]), n_null_sets=int(prc["n_null_sets"]),
            kappa_min_effect=prc["kappa_min_effect"], seed=seed)
        tables["protocol_shift"] = per_mp
        tables["protocol_tumours"] = per_tumour
        tables["protocol_levels"] = levels
        log("[p1] protocol shift done")

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

    rule = rule_text(cfg["decision"], kappa_min_effect=prc["kappa_min_effect"],
                     rho_sensitivity=cfg["tme"]["rho_sensitivity"])
    (out / "run_config.json").write_text(json.dumps(
        {"config": plain(cfg), "params": vars(params), "rule_version": RULE_VERSION, "rule": rule,
         "versions": versions, "python": platform.python_version(), "decision": decision}, indent=2, default=str))
    log(f"[p1] decision: {decision['verdict']} ({decision['n_counted']} counted MPs; {RULE_VERSION}); "
        f"see {out / 'summary.md'}")
    return {"decision": decision, "tables": tables, "reference": reference, "mp_sets": mp_sets}
