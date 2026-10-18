"""Pilot P1 runner: audit a set of MPs (re-derived or published) along every available axis and
apply the frozen GO rule (:mod:`mpaudit.report`, ``P1-GO v1.3``).

:func:`validate_config` checks the whole configuration before any NMF runs: unknown keys, the
``decision`` block (:func:`mpaudit.report.validate_rule`), the keys against the cohort's obs
columns, the decoy arms against the rule's arm lists and the cohort's obs columns (no NaN in
malignant cells), the label, gene-set, cytoband and GTF files, names that cannot match re-derived
MPs, the GavishParams/keys combination, the registered settings against their registered values
and the inputs the dosage and curated-artefact exclusions need (:func:`config_deviations`).
:func:`run_p1` then checks, still before any NMF, that the gene positions and the cytoband place
genes on chromosome arms. A config that deviates is refused unless its top-level
``unregistered_variant`` is true; the verdict is then labelled an unregistered variant. An axis
that is switched on but cannot run (e.g. the protocol axis without protocol names) is logged and
noted in ``summary.md``; the protocol case is also warned about by :func:`validate_config`.

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
import hashlib
import json
import platform
import re
import warnings
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from .compare import jaccard_matrix, match_sets
from .data import QC_COLUMNS, Keys, add_qc_obs, malignant_mask, read_gmt, sample_map, write_gmt
from .decoy import MIN_FLOOR_SETS, N_NULL_FLOOR, decoy_recruitment
from .gavish import GavishParams, NMFCache, check_keys, code_version, eligible_universe, run_gavish
from .instability import default_perturbations, run_instability
from .positional import (arm_enrichment, assign_arms, centromeres, gene_lengths, gene_positions_from_gtf,
                         load_cytoband)
from .protocol import KAPPA_ACTIVE, KAPPA_MIN_EFFECT, MIN_PAIRS, MIN_RANDOM_SETS, protocol_shift_kappa
from .report import (DEFAULT_RULE, MIN_JACCARD, RHO_SENSITIVITY, RULE_VERSION, decide, rule_deviations, rule_label,
                     rule_settings, rule_text, validate_rule, write_summary)
from .tme import NULL_MATCH, tme_covariation
from .typer import evaluate_typer, typer_features

DEFAULT_CONFIG = {
    "keys": {"sample": "sample", "malignant": "malignant"},
    "counts_layer": None,
    "recurrence_key": None,       # obs column defining independent units (e.g. tumour)
    "mps": "rederive",            # or path to a GMT of published MPs to audit
    "gavish": {},                 # GavishParams overrides
    "n_jobs": 1,
    "cache_dir": None,
    "min_jaccard": MIN_JACCARD,   # MP match threshold: decoy hits (criterion a), label carry-over, instability
    # true: run a config that deviates from the registered rule (mpaudit.report.rule_deviations); the
    # verdict is then labelled "P1-GO v1.3 (unregistered variant: ...)". false: such a config is refused.
    "unregistered_variant": False,
    "qc": {"add_obs": True},      # add log_umi, log_genes, pct_mito, pct_ribo to obs if missing
    "instability": {"enabled": True, "seeds": [1, 2, 3], "n_boot": 3, "boot_fraction": 0.8,
                    "leave_one_study_out": True},
    # arms: name -> obs columns the decoys condition on; categorical columns become permutation
    # strata, continuous ones enter the mean model. Every arm must be named in
    # decision.technical_arms or decision.dosage_arms. "strata" (one unnamed arm) is no longer accepted.
    # min_kappa, n_null_floor: kappa_retained's entry, an excess over the median of n_null_floor matched random
    # gene sets above max(min_kappa, 3 SD of theirs) (mpaudit.decoy)
    "decoy": {"enabled": True, "n": 5, "method": "residual", "arms": {}, "strata": [], "n_depth_bins": 1,
              "min_kappa": 0.02, "n_null_floor": N_NULL_FLOOR},
    "tme": {"enabled": True, "min_cells": 20, "n_boot": 200, "n_null": 500, "covariates": None,
            "spillover_adjust": True, "null_match": "expression+specificity", "n_markers": 50,
            "min_marker_fold": 8.0, "rho_sensitivity": RHO_SENSITIVITY},
    "positional": {"enabled": True, "cytoband": None, "gtf": None,
                   "chrom_col": "chromosome", "start_col": "start"},
    # labels: CSV with columns mp,label (1 = curated artefact); labelled_gmt: gene sets of the
    # labelled MPs, needed to carry the labels onto re-derived MPs
    "typer": {"enabled": True, "labels": None, "labelled_gmt": None},
    # kappa_active: informational only (not part of the rule, no gate): the floor of the activity threshold
    # max(kappa_active, random-set q95) behind the kappa_active* columns of protocol_tumours (mpaudit.protocol)
    "protocol": {"enabled": True, "reference": None, "test": None, "min_pairs": MIN_PAIRS, "n_null_sets": 100,
                 "kappa_min_effect": KAPPA_MIN_EFFECT, "kappa_active": KAPPA_ACTIVE},
    "decision": dict(DEFAULT_RULE),
    "seed": 0,                    # registered: decoys, random-gene-set nulls, cell subsamples, depth thinning
                                  # and bootstrap resamples
}
MIN_POSITION_COVERAGE = 0.5  # warn when gene positions or lengths cover fewer var_names than this
DECOY_METHODS = ("residual", "permute")
_REDERIVED_NAME = re.compile(r"MP\d+")  # names mpaudit.gavish gives re-derived MPs
# sections whose keys validate_config checks against DEFAULT_CONFIG; gavish (GavishParams), keys (Keys) and
# decision (validate_rule) check their own keys
_CHECKED_SECTIONS = ("qc", "instability", "decoy", "tme", "positional", "typer", "protocol")
_SWITCHES = (("qc", "add_obs"), ("instability", "enabled"), ("instability", "leave_one_study_out"),
             ("tme", "enabled"), ("typer", "enabled"), ("protocol", "enabled"))  # true/false settings
# keys only the informational typer reads (doublet, intronic) or nothing reads (clone, region): one naming no obs
# column is warned about, not refused (the typer skips a missing column)
_INFORMATIONAL_KEYS = ("clone", "region", "doublet", "intronic")


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


def _name(section: str, key: str) -> str:
    return f"{section}.{key}" if section else key


def _optional_number(section: str, key: str, v, *, allow_none: bool, minimum: float, strict: bool):
    if v is None and allow_none:
        return None
    try:
        if isinstance(v, bool):
            raise TypeError
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{_name(section, key)} must be a number{' or null' if allow_none else ''}, "
                         f"got {v!r}") from None
    if not np.isfinite(x) or x < minimum or (strict and x == minimum):
        raise ValueError(f"{_name(section, key)} must be {'>' if strict else '>='} {minimum:g}, got {v!r}")
    return x


def _whole_number(section: str, key: str, v, *, minimum: int) -> int:
    """``v`` as an int >= ``minimum``; bools, floats (also 3.0) and strings are rejected."""
    if isinstance(v, (bool, np.bool_)) or not isinstance(v, (int, np.integer)):
        raise ValueError(f"{_name(section, key)} must be a whole number, got {v!r}")
    if v < minimum:
        raise ValueError(f"{_name(section, key)} must be >= {minimum}, got {v!r}")
    return int(v)


def _flag(section: str, key: str, v) -> bool:
    """``v`` as a bool; only true/false are accepted (a string such as 'false' would be truthy)."""
    if not isinstance(v, (bool, np.bool_)):
        raise ValueError(f"{_name(section, key)} must be true or false, got {v!r}")
    return bool(v)


def _check_known_keys(cfg: dict) -> None:
    """Raise ``ValueError`` for config keys that no part of mpaudit reads: unknown top-level keys and
    unknown keys in the sections of :data:`_CHECKED_SECTIONS` (all of them listed, with the valid
    keys). A misspelt key would otherwise run with the default silently: ``tme.rho_sensitvity:
    null``, meant as a variant, would run as registered, and ``protocol.referance`` would leave
    criterion (c) unpowered. The ``decision``, ``keys`` and ``gavish`` blocks refuse unknown keys
    themselves."""
    bad = [f"{k!r}" for k in cfg if k not in DEFAULT_CONFIG]
    for sec in _CHECKED_SECTIONS:
        if not isinstance(cfg[sec], dict):
            raise ValueError(f"config section {sec} must be a mapping of settings, got {cfg[sec]!r}")
        bad += [f"'{sec}.{k}' (valid: {', '.join(DEFAULT_CONFIG[sec])})" for k in cfg[sec]
                if k not in DEFAULT_CONFIG[sec]]
    if bad:
        raise ValueError(f"unknown config key(s), probably misspelt (no part of mpaudit reads them): {'; '.join(bad)}. "
                         f"Valid top-level keys: {', '.join(DEFAULT_CONFIG)}")


def _is_file(path) -> bool:
    return isinstance(path, (str, Path)) and bool(str(path)) and Path(path).is_file()


def _check_protocol_axis(cfg: dict, keys: Keys, adata, log) -> None:
    """Refuse equal protocol names; warn when the protocol axis is on and ``keys.protocol`` is set but
    the axis cannot pair a tumour (a protocol name missing, or with ``adata``, a name that is not a
    value of the protocol column, or no tumour with both protocols), so criterion (c) would be
    unpowered without a word."""
    prc = cfg["protocol"]
    names = {k: prc.get(k) for k in ("reference", "test")}
    unset = [k for k, v in names.items() if v is None or not str(v)]
    if not unset and str(names["reference"]) == str(names["test"]):
        raise ValueError(f"protocol.reference and protocol.test must differ, got {names['reference']!r} for both")
    if not (prc["enabled"] and keys.protocol):
        return
    if unset:
        _warn(f"protocol.{' and protocol.'.join(unset)} not set: the protocol axis will not run, so criterion (c) "
              "cannot be powered", log)
        return
    if adata is None:
        return
    prot = adata.obs[keys.protocol].astype(str)
    levels = sorted(pd.unique(prot))
    absent = [f"protocol.{k} {v!r}" for k, v in names.items() if str(v) not in levels]
    if absent:
        _warn(f"{' and '.join(absent)} not a value of adata.obs[{keys.protocol!r}] ({levels[:10]}): no tumour can be "
              "paired, so criterion (c) cannot be powered", log)
        return
    mal = malignant_mask(adata, keys)
    tumour = adata.obs[keys.tumour_col()].astype(str)[mal]
    prot = prot[mal]
    if not set(tumour[prot == str(names["reference"])]) & set(tumour[prot == str(names["test"])]):
        _warn(f"no {keys.tumour_col()} has malignant cells under both protocol.reference {names['reference']!r} and "
              f"protocol.test {names['test']!r}: no tumour can be paired, so criterion (c) cannot be powered", log)


def _has_positions(adata, pc: dict) -> bool:
    """Whether ``adata.var`` has the gene position columns the positional axis falls back to."""
    return pc.get("chrom_col", "chromosome") in adata.var and pc.get("start_col", "start") in adata.var


def config_deviations(cfg: dict, adata=None) -> list[tuple[str, object, str]]:
    """Deviations of a validated config from the registered rule, as ``(setting, value, registered)``:
    :func:`mpaudit.report.rule_deviations`, plus inputs without which an exclusion of the rule does
    not run although its axis is on:

    * with ``mps: rederive``, curated labels with the gene sets of the labelled MPs
      (``typer.labels`` and ``typer.labelled_gmt``): without labels no rediscovered curated artefact
      is excluded, and labels by name cannot be registered, since the run names the MPs;
    * with ``decoy.enabled``, an arm named in ``decision.dosage_arms`` (the dosage-arm exclusion);
    * with ``positional.enabled``, a cytoband and gene positions (``positional.gtf``, or with
      ``adata`` given, the ``chrom_col``/``start_col`` columns of ``adata.var``) for the
      arm-enrichment exclusion. Positions that place no var_name (e.g. a GTF whose gene names do not
      match them) and a cytoband whose centromere bands place none of them on an arm are deviations
      too; :func:`run_p1` finds them from the data, still before any NMF.
    """
    devs = rule_deviations(cfg)
    tc, pc, dc = cfg["typer"], cfg["positional"], cfg["decoy"]
    if cfg["mps"] == "rederive" and not tc.get("labels"):
        devs.append(("typer.labels", None, "curated labels with typer.labelled_gmt (needed with mps: rederive)"))
    elif cfg["mps"] == "rederive" and not tc.get("labelled_gmt"):
        devs.append(("typer.labelled_gmt", None, "a GMT of the labelled MPs (needed with mps: rederive and "
                                                  "typer.labels)"))
    if dc.get("enabled", True):
        arms = sorted(_decoy_arms(dc))
        if not set(arms) & set(cfg["decision"]["dosage_arms"]):
            devs.append(("decoy.arms", arms, "an arm named in decision.dosage_arms (e.g. CNA clone)"))
    if pc.get("enabled", True):
        if pc.get("cytoband") is None:
            devs.append(("positional.cytoband", None, "a cytoband table (arm enrichment types dosage MPs)"))
        if adata is not None and not pc.get("gtf") and not _has_positions(adata, pc):
            devs.append(("positional.gtf", None, "gene positions: a GTF, or chromosome/start columns in adata.var"))
    return devs


def _arm_columns(arms: dict) -> list[str]:
    return sorted({c for cols in arms.values() for c in cols})


MISSING_STRINGS = ("", "nan", "none", "na", "<na>", "null")  # placeholders of a missing value in text columns


def check_arm_values(adata, keys: Keys, arms: dict) -> None:
    """Raise ``ValueError`` for decoy arm columns with missing values in malignant cells: NaN, None
    or infinite values, and in text columns placeholders such as ``'nan'`` or ``''``
    (:data:`MISSING_STRINGS`, after stripping, any case), which a column written with
    ``astype(str)`` or read back from CSV carries. Columns not in ``adata.obs`` are skipped.

    :func:`mpaudit.decoy.make_decoy` does not model missing values: per sample, a continuous
    column's NaN take the sample's median, a categorical column's NaN form a stratum of their own,
    and a numeric column left with <= 10 distinct values becomes a set of strata, so a few measured
    cells can each form a stratum whose residuals are never permuted. Decoys would condition on
    values nobody measured.
    """
    mal = malignant_mask(adata, keys)
    samples = adata.obs[keys.sample].astype(str).to_numpy()[mal]
    bad = []
    for col in _arm_columns(arms):
        if col not in adata.obs:
            continue
        v = adata.obs[col][mal]
        if pd.api.types.is_numeric_dtype(v) and not pd.api.types.is_bool_dtype(v):
            missing = ~np.isfinite(v.to_numpy(dtype=float))
        else:
            text = v.astype(object).map(lambda x: x.strip().lower() if isinstance(x, str) else x)
            missing = (v.isna() | text.isin(MISSING_STRINGS)).to_numpy(dtype=bool)
        if missing.any():
            where = pd.unique(samples[missing])
            bad.append(f"{col} ({int(missing.sum())} of {missing.size} malignant cells, in {where.size} sample(s): "
                       f"{', '.join(map(str, where[:5]))}{', ...' if where.size > 5 else ''})")
    if bad:
        raise ValueError(
            f"decoy arm column(s) with missing values (NaN, or a text placeholder such as 'nan') in malignant "
            f"cells: {'; '.join(bad)}. The decoys do "
            "not model missing values (a continuous column's NaN would take the sample's median, a categorical "
            "column's NaN would form a stratum of their own). Drop the column from decoy.arms, or fill it in "
            "obs: impute it, give unassigned cells of a categorical column an explicit label (e.g. 'unassigned'), "
            "or, to condition only where the covariate was measured (e.g. intronic_fraction for some studies), "
            "set it to one constant for every cell of a sample without it (a column constant within a sample "
            "does not condition that sample's decoys)")


def validate_config(cfg: dict, adata=None, log=print) -> dict:
    """Check a merged config (:func:`merge_config`) before anything runs; returns it normalised.

    Raises ``ValueError`` for: unknown keys at the top level or in a section other than
    ``decision``, ``keys`` and ``gavish`` (:func:`_check_known_keys`; those three refuse unknown
    keys themselves); an invalid ``decision`` block (:func:`mpaudit.report.validate_rule`);
    with ``adata``, a key (``keys.*``) naming no column of ``adata.obs`` (warned about instead for
    the informational keys ``clone``, ``region``, ``doublet`` and ``intronic``); a switch
    (``*.enabled``, ``qc.add_obs``, ``instability.leave_one_study_out``) that is not true/false;
    ``decoy.strata`` (legacy; use ``decoy.arms``); a decoy arm named in neither
    ``decision.technical_arms`` nor ``decision.dosage_arms``; decoy arm columns missing from
    ``adata.obs`` (QC columns that ``qc.add_obs`` will add count as present; ``pct_long`` only
    with gene lengths from ``positional.gtf`` or ``adata.var['gene_length']``) or with NaN in
    malignant cells (:func:`check_arm_values`); ``decision.expected_tme_shared`` names with
    ``mps: rederive``; a missing ``mps`` GMT, ``typer.labels``, ``typer.labelled_gmt`` or
    ``positional.gtf`` file, or, with the positional axis on, a ``positional.cytoband`` that is
    neither null (a deviation) nor a table nor a file; a ``labelled_gmt`` without ``labels``, or,
    with ``mps: rederive`` and no ``labelled_gmt``, labels none of which can name a re-derived MP;
    equal ``protocol.reference`` and ``protocol.test``; values of the wrong type or out of range for
    the settings the rule reads (``tme.rho_sensitivity`` a number > 0 or null,
    ``protocol.kappa_min_effect`` a number >= 0, ``protocol.min_pairs`` a whole number >= 1,
    ``protocol.n_null_sets`` 0 or a whole number >= ``mpaudit.protocol.MIN_RANDOM_SETS`` (20),
    ``decoy.n_null_floor`` a whole number >= ``mpaudit.decoy.MIN_FLOOR_SETS`` (20), ``seed`` a whole
    number >= 0, and the other :data:`mpaudit.report.REGISTERED_SETTINGS`), for the informational
    ``protocol.kappa_active`` not a number >= 0 and for ``unregistered_variant`` not true/false;
    what :func:`mpaudit.gavish.check_keys` rejects; and, unless ``unregistered_variant`` is true,
    any deviation from the registered rule (:func:`config_deviations`; all of them listed, with
    their registered values). With ``unregistered_variant`` true the deviations are warned about
    instead.

    Warns when the TME axis can run (``tme.enabled`` and ``keys.cell_type``) but no
    expected-TME-shared MP or gene set is pre-declared; when MPs come from a GMT without
    ``typer.labels`` (no MP is then treated as a curated artefact); and when the protocol axis is
    on and ``keys.protocol`` is set but no tumour can be paired (:func:`_check_protocol_axis`:
    a protocol name not set, or with ``adata``, not a value of the protocol column or never paired
    in one tumour), so criterion (c) would be unpowered. Name checks against the audited MPs
    themselves are in :func:`check_names`; whether the gene positions and the cytoband place genes
    on chromosome arms is checked by :func:`run_p1`, also before any NMF.
    """
    _check_known_keys(cfg)
    cfg["decision"] = rule = validate_rule(cfg.get("decision"))
    rederive = cfg["mps"] == "rederive"
    if not rederive and not Path(str(cfg["mps"])).is_file():
        raise ValueError(f"mps must be 'rederive' or a GMT file, got {cfg['mps']!r}")
    keys = Keys.from_dict(cfg["keys"])
    check_keys(keys, GavishParams.from_dict(cfg["gavish"]), adata, recurrence_key=cfg["recurrence_key"])
    if adata is not None:
        cols = {k: v for k, v in vars(keys).items() if v and k != "malignant_label"}
        missing = [f"keys.{k}={v!r}" for k, v in cols.items() if v not in adata.obs.columns
                   and k not in _INFORMATIONAL_KEYS]
        if missing:
            raise ValueError(f"Not a column of adata.obs: {', '.join(missing)}")
        ignored = [f"keys.{k}={v!r}" for k, v in cols.items()
                   if v not in adata.obs.columns and k in _INFORMATIONAL_KEYS]
        if ignored:
            _warn(f"not a column of adata.obs, ignored (informational keys: the rule does not read them): "
                  f"{', '.join(ignored)}", log)
    for sec, key in _SWITCHES:
        cfg[sec][key] = _flag(sec, key, cfg[sec].get(key, DEFAULT_CONFIG[sec][key]))

    dc = cfg["decoy"]
    arms = _decoy_arms(dc)
    dc["enabled"] = _flag("decoy", "enabled", dc.get("enabled", True))
    if dc["enabled"]:
        known = set(rule["technical_arms"]) | set(rule["dosage_arms"])
        stray = sorted(set(arms) - known)
        if stray:
            raise ValueError(f"decoy arm(s) {stray} are in neither decision.technical_arms "
                             f"{rule['technical_arms']} nor decision.dosage_arms {rule['dosage_arms']}: the rule "
                             "would ignore them (a CNA clone arm not listed as dosage would not type dosage)")
        if adata is not None:
            have = set(adata.obs.columns)
            if cfg["qc"].get("add_obs", True):
                have |= set(QC_COLUMNS)
                if cfg["positional"].get("gtf") or "gene_length" in adata.var:
                    have.add("pct_long")
            missing = sorted(set(_arm_columns(arms)) - have)
            if missing:
                raise ValueError(f"decoy arm column(s) {missing} are not columns of adata.obs"
                                 + (" (pct_long needs gene lengths: positional.gtf or adata.var['gene_length'])"
                                    if "pct_long" in missing else ""))
            check_arm_values(adata, keys, arms)
    dc["arms"] = {k: list(v) for k, v in arms.items()}
    dc["n"] = _whole_number("decoy", "n", dc.get("n", 5), minimum=1)
    if dc.get("method", "residual") not in DECOY_METHODS:
        raise ValueError(f"decoy.method must be one of {DECOY_METHODS}, got {dc.get('method')!r}")
    dc["n_depth_bins"] = _whole_number("decoy", "n_depth_bins", dc.get("n_depth_bins", 1), minimum=1)
    dc["min_kappa"] = _optional_number("decoy", "min_kappa", dc.get("min_kappa", 0.02), allow_none=False,
                                       minimum=0.0, strict=False)
    dc["n_null_floor"] = _whole_number("decoy", "n_null_floor", dc.get("n_null_floor", N_NULL_FLOOR),
                                       minimum=MIN_FLOOR_SETS)

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
    elif not rederive:  # with mps: rederive, missing labels are a deviation (config_deviations)
        _warn("no typer.labels: no audited MP is treated as a curated artefact, so none is excluded from the count "
              "(fine only if the mps GMT already leaves them out, as the published 3CA v2 list does)", log)

    tm = cfg["tme"]
    tm["rho_sensitivity"] = _optional_number("tme", "rho_sensitivity", tm.get("rho_sensitivity"),
                                             allow_none=True, minimum=0.0, strict=True)
    tm["spillover_adjust"] = _flag("tme", "spillover_adjust", tm.get("spillover_adjust", True))
    if tm.get("null_match", NULL_MATCH[0]) not in NULL_MATCH:
        raise ValueError(f"tme.null_match must be one of {NULL_MATCH}, got {tm.get('null_match')!r}")
    tm["n_null"] = _whole_number("tme", "n_null", tm.get("n_null", 500), minimum=0)
    tm["min_cells"] = _whole_number("tme", "min_cells", tm.get("min_cells", 20), minimum=1)
    tm["n_markers"] = _whole_number("tme", "n_markers", tm.get("n_markers", 50), minimum=1)
    tm["min_marker_fold"] = _optional_number("tme", "min_marker_fold", tm.get("min_marker_fold", 8.0),
                                             allow_none=False, minimum=0.0, strict=True)
    if tm.get("enabled", True) and keys.cell_type and not (
            rule["expected_tme_shared"] or rule["expected_tme_shared_sets"] or rule["expected_tme_shared_gmt"]):
        _warn("no expected-TME-shared MPs are pre-declared (decision.expected_tme_shared / expected_tme_shared_gmt / "
              f"expected_tme_shared_sets): IFN/MHC-II or hypoxia MPs can count through criterion (b); {RULE_VERSION} "
              "expects this list to be filled in and registered", log)
    prc = cfg["protocol"]
    prc["kappa_min_effect"] = _optional_number("protocol", "kappa_min_effect",
                                               prc.get("kappa_min_effect", KAPPA_MIN_EFFECT),
                                               allow_none=False, minimum=0.0, strict=False)
    prc["kappa_active"] = _optional_number("protocol", "kappa_active", prc.get("kappa_active", KAPPA_ACTIVE),
                                           allow_none=False, minimum=0.0, strict=False)
    prc["min_pairs"] = _whole_number("protocol", "min_pairs", prc.get("min_pairs", MIN_PAIRS), minimum=1)
    prc["n_null_sets"] = _whole_number("protocol", "n_null_sets", prc.get("n_null_sets", 100), minimum=0)
    if 0 < prc["n_null_sets"] < MIN_RANDOM_SETS:  # also refused as an unregistered variant: it tests nothing
        raise ValueError(f"protocol.n_null_sets must be 0 (no random-set null) or >= {MIN_RANDOM_SETS}, got "
                         f"{prc['n_null_sets']}: a tumour's random-set bound needs {MIN_RANDOM_SETS} random sets "
                         "with a defined change, so every tumour on the random-set null would be left untested")
    _check_protocol_axis(cfg, keys, adata, log)
    pc = cfg["positional"]
    pc["enabled"] = _flag("positional", "enabled", pc.get("enabled", True))
    # the files the positional axis reads, checked now rather than when the axis runs (after the NMF); a cytoband
    # left null is a deviation (config_deviations), and the GTF is read even with the axis off (gene lengths)
    if pc["enabled"] and pc.get("cytoband") is not None and not isinstance(pc["cytoband"], pd.DataFrame) \
            and not _is_file(pc["cytoband"]):
        raise ValueError(f"positional.cytoband {pc['cytoband']!r} is not a file")
    if pc.get("gtf") and not _is_file(pc["gtf"]):
        raise ValueError(f"positional.gtf {pc['gtf']!r} is not a file")
    mj = _optional_number("", "min_jaccard", cfg.get("min_jaccard", MIN_JACCARD), allow_none=False, minimum=0.0,
                          strict=True)
    if mj > 1:
        raise ValueError(f"min_jaccard must be <= 1, got {cfg['min_jaccard']!r}")
    cfg["min_jaccard"] = mj
    cfg["unregistered_variant"] = _flag("", "unregistered_variant", cfg.get("unregistered_variant", False))
    cfg["seed"] = _whole_number("", "seed", cfg.get("seed", 0), minimum=0)

    devs = config_deviations(cfg, adata)
    if devs and not cfg["unregistered_variant"]:
        raise ValueError(
            f"the config deviates from the registered rule {RULE_VERSION} in {len(devs)} setting(s): "
            + "; ".join(f"{k} = {v!r} (registered: {r})" for k, v, r in devs)
            + ". Restore the registered values, or set unregistered_variant: true to run anyway (the verdict is then "
              f"labelled an unregistered variant and is not {RULE_VERSION})")
    if devs:
        _warn(f"unregistered variant: the verdict will be labelled '{rule_label(devs)}', not {RULE_VERSION}", log)
    elif cfg["unregistered_variant"]:
        log(f"[p1] unregistered_variant is true but no setting deviates: the run is the registered {RULE_VERSION}")
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
    """Gene positions and lengths from ``positional.gtf``, or None. An empty table, or one that places
    no var_name (e.g. Ensembl IDs as var_names), counts as missing (callers fall back to
    ``adata.var``); one placing fewer than ``MIN_POSITION_COVERAGE`` of the var_names warns."""
    if not pc.get("gtf"):
        return None
    table = gene_positions_from_gtf(pc["gtf"])
    if table.empty:
        _warn(f"positional.gtf {pc['gtf']} gave no gene positions; falling back to adata.var columns", log)
        return None
    cov = _coverage(table.index, adata.var_names)
    if cov == 0:
        _warn(f"positional.gtf {pc['gtf']} places no var_name (they must match its gene_name values; Ensembl IDs do "
              "not); falling back to adata.var columns", log)
        return None
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


def _fit_totals(runs: pd.DataFrame | None) -> dict | None:
    """Fit diagnostics summed over runs (rows of ``decoy_runs`` / ``instability_runs``), or None."""
    if runs is None or runs.empty or "n_fits" not in runs:
        return None
    return {"n_runs": int(len(runs)), **{c: int(runs[c].sum()) for c in ("n_fits", "n_capped_fits", "n_failed_fits")}}


def input_file_digests(cfg: dict) -> dict[str, dict[str, str]]:
    """Setting -> ``{"path", "sha256"}`` for each input file the run reads (the audited MPs' GMT,
    the curated labels and their gene sets, the expected-TME-shared GMT, the cytoband and GTF), so
    the files behind a verdict can be checked against the registered ones."""
    pc, tc = cfg["positional"], cfg["typer"]
    paths = {"mps": None if cfg["mps"] == "rederive" else cfg["mps"], "typer.labels": tc.get("labels"),
             "typer.labelled_gmt": tc.get("labelled_gmt"),
             "decision.expected_tme_shared_gmt": cfg["decision"].get("expected_tme_shared_gmt"),
             "positional.cytoband": pc.get("cytoband"), "positional.gtf": pc.get("gtf")}
    out = {}
    for setting, path in paths.items():
        if isinstance(path, (str, Path)) and str(path) and Path(path).is_file():
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for block in iter(lambda: fh.read(1 << 20), b""):
                    h.update(block)
            out[setting] = {"path": str(path), "sha256": h.hexdigest()}
    return out


def _warn_failed_fits(what: str, fits: dict | None, log) -> None:
    if fits and fits["n_failed_fits"]:
        _warn(f"{what}: {fits['n_failed_fits']} of {fits['n_fits']} NMF fits were stopped by the restart limit "
              "(all-zero rows of H; such a fit is unusable)", log)


def _late_deviation(dev: tuple[str, object, str], problem: str, cfg: dict, deviations: list, log) -> None:
    """A deviation found by :func:`run_p1` (still before any NMF) from the data rather than the config:
    raise ``ValueError`` unless ``unregistered_variant`` is true, else append it and warn."""
    if not cfg["unregistered_variant"]:
        raise ValueError(f"{problem}: {dev[0]} = {dev[1]!r} (registered: {dev[2]}); this deviates from the registered "
                         f"rule {RULE_VERSION}. Fix the input, or set unregistered_variant: true to run anyway")
    deviations.append(dev)
    _warn(f"unregistered variant: the verdict will be labelled '{rule_label(deviations)}', not {RULE_VERSION}", log)


def run_p1(adata, config: dict | None = None, out_dir: str | Path = "results", log=print) -> dict:
    """Run the P1 audit and write its outputs to ``out_dir``; returns the decision, the tables, the
    reference run and the audited MP sets.

    Besides one CSV per table: ``summary.md`` (:func:`mpaudit.report.write_summary`) and
    ``run_config.json`` (the normalised config, GavishParams, the rule label and text, the
    deviations from the registered rule, the sha256 of the input files (:func:`input_file_digests`),
    package versions, the decision, ``nmf_fits``: capped and failed NMF fits of the reference run
    and summed over the decoy and instability runs, ``nmf_cache_code``: the code hashes in the NMF
    cache keys (:func:`mpaudit.gavish.code_version`), and ``notes``: axes that did not run and why,
    and how many genes were placed on chromosome arms; the notes also go into ``summary.md``).

    Everything that can be checked without NMF is checked first: the config
    (:func:`validate_config`), then the gene positions and the cytoband. With the positional axis
    on, positions that place no var_name, or a cytoband whose centromere (``acen``) bands place none
    of them on a chromosome arm, leave the arm-enrichment exclusion without input: a deviation,
    refused unless ``unregistered_variant`` is true (fewer than ``MIN_POSITION_COVERAGE`` placed
    warns).
    """
    cfg = validate_config(merge_config(config), adata, log)
    deviations = config_deviations(cfg, adata)
    rederive = cfg["mps"] == "rederive"
    mp_sets = None if rederive else read_gmt(cfg["mps"])
    if mp_sets is not None:
        check_names(mp_sets, cfg, log)
    keys = Keys.from_dict(cfg["keys"])
    params = GavishParams.from_dict(cfg["gavish"])
    cache = NMFCache(cfg["cache_dir"])
    layer, rk, n_jobs, seed = cfg["counts_layer"], cfg["recurrence_key"], int(cfg["n_jobs"]), int(cfg["seed"])
    tables: dict[str, pd.DataFrame] = {}
    notes: list[str] = []  # for summary.md and run_config.json: axes that did not run, genes placed on arms
    pc = cfg["positional"]
    # read once, before the QC columns (pct_long), for the positional axis and the typer
    gtf_positions = _gtf_positions(adata, pc, log) if pc.get("gtf") else None
    positions, arms_of_gene = None, None
    if pc["enabled"]:
        positions = gtf_positions if gtf_positions is not None else _var_positions(adata, pc, log)
        if (positions is None or positions.empty) and not any(d[0] == "positional.gtf" for d in deviations):
            # e.g. a GTF that places no var_name, with no positions in adata.var: no arm-enrichment exclusion
            _late_deviation(("positional.gtf", pc.get("gtf"),
                             "gene positions: a GTF, or chromosome/start columns in adata.var"),
                            "no gene positions for the positional axis", cfg, deviations, log)
        cyto = pc.get("cytoband")
        if positions is not None and not positions.empty and cyto is not None:
            # arms now, before any NMF: a cytoband without acen bands on the genes' chromosomes places no gene
            arms_of_gene = assign_arms(positions, centromeres(cyto if isinstance(cyto, pd.DataFrame)
                                                              else load_cytoband(cyto)))
            n_placed = int(pd.Index(adata.var_names.astype(str)).isin(arms_of_gene.dropna().index.astype(str)).sum())
            share = n_placed / max(adata.n_vars, 1)
            notes.append(f"Gene positions and the cytoband place {n_placed} of {adata.n_vars} var_names ({share:.0%}) "
                         "on chromosome arms (positional axis).")
            log(f"[p1] chromosome arms: {n_placed} of {adata.n_vars} var_names placed")
            if n_placed == 0:
                arms_of_gene = None
                _late_deviation(("positional.cytoband", cyto if isinstance(cyto, (str, Path)) else "<table>",
                                 "a cytoband whose acen bands place the genes on chromosome arms"),
                                "the cytoband places no gene on a chromosome arm (no centromere band on the genes' "
                                "chromosomes)", cfg, deviations, log)
            elif share < MIN_POSITION_COVERAGE:
                _warn(f"gene positions and the cytoband place only {share:.0%} of var_names on chromosome arms", log)
    if cfg["qc"].get("add_obs", True):
        add_qc_obs(adata, counts_layer=layer,
                   gene_length=gene_lengths(gtf_positions) if gtf_positions is not None else None)
    arms = {k: tuple(v) for k, v in cfg["decoy"]["arms"].items()}
    if cfg["decoy"]["enabled"]:
        missing = sorted(set(_arm_columns(arms)) - set(adata.obs.columns))
        if missing:  # e.g. pct_long when positional.gtf placed no gene and adata.var has no gene_length
            raise ValueError(f"decoy arm column(s) {missing} are not columns of adata.obs after the QC columns were "
                             "added")
        check_arm_values(adata, keys, arms)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    reference = None
    if rederive or cfg["instability"]["enabled"]:
        log("[p1] running the reference 3CA-style pipeline")
        reference = run_gavish(adata, keys, params, counts_layer=layer, cache=cache, recurrence_key=rk, n_jobs=n_jobs)
        write_gmt(reference.mp_sets(), out / "reference_mps.gmt")
        reference.programs_table().to_csv(out / "programs.csv", index=False)
        log(f"[p1] reference: {len(reference.robust)} robust programs, {len(reference.mps)} MPs; "
            f"{reference.n_capped_fits}/{reference.n_fits} NMF fits stopped at max_iter, "
            f"{reference.n_failed_fits} by the restart limit")
        _warn_failed_fits("reference run", {"n_fits": reference.n_fits, "n_failed_fits": reference.n_failed_fits},
                          log)
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
        _warn_failed_fits("instability runs", _fit_totals(inst.summary), log)
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
    if dc["enabled"] and mp_sets and arms:
        per_mp, per_decoy, _ = decoy_recruitment(
            adata, keys, params, mp_sets, arms=arms, n_decoys=int(dc["n"]), n_depth_bins=int(dc["n_depth_bins"]),
            method=dc["method"], counts_layer=layer, cache=cache, recurrence_key=rk, min_jaccard=cfg["min_jaccard"],
            min_kappa=float(dc["min_kappa"]), n_null_floor=int(dc["n_null_floor"]), seed=1000 + seed,
            n_jobs=n_jobs, log=log)
        tables["decoy_recruitment"] = per_mp
        tables["decoy_runs"] = per_decoy
        _warn_failed_fits("decoy runs", _fit_totals(per_decoy), log)
    elif dc["enabled"] and mp_sets:
        log("[p1] decoys skipped: no decoy arms configured")
        notes.append("The decoy axis did not run (no decoy arms configured), so criterion (a) is not powered.")

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
    elif cfg["tme"]["enabled"] and mp_sets:
        log("[p1] TME axis skipped (needs keys.cell_type)")
        notes.append("The TME axis did not run (needs keys.cell_type), so criterion (b) is not powered.")

    if pc["enabled"] and mp_sets:
        if arms_of_gene is not None:
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
            why = ("no gene placed on a chromosome arm" if positions is not None and not positions.empty
                   and pc.get("cytoband") is not None else "needs gene positions and a cytoband")
            log(f"[p1] positional axis skipped ({why})")
            notes.append(f"The positional axis did not run ({why}), so no MP can be typed dosage by arm enrichment.")

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
            min_pairs=prc["min_pairs"], n_null_sets=int(prc["n_null_sets"]),
            kappa_min_effect=prc["kappa_min_effect"], kappa_active=prc["kappa_active"], seed=seed)
        tables["protocol_shift"] = per_mp
        tables["protocol_tumours"] = per_tumour
        tables["protocol_levels"] = levels
        log("[p1] protocol shift done")
    elif prc["enabled"] and mp_sets:
        why = (f"needs keys.protocol, protocol.reference and protocol.test; got {keys.protocol!r}, "
               f"{prc.get('reference')!r}, {prc.get('test')!r}")
        log(f"[p1] protocol axis skipped ({why})")
        notes.append(f"The protocol axis did not run ({why}), so criterion (c) is not powered.")

    decision, flags = decide(mp_sets, tables, cfg["decision"], curated_artefacts=labels, deviations=deviations)
    tables["mp_flags"] = flags
    for name, df in tables.items():
        df.to_csv(out / f"{name}.csv", index=False)
    nmf_fits = {"reference": None if reference is None else {
                    "n_fits": reference.n_fits, "n_capped_fits": reference.n_capped_fits,
                    "n_failed_fits": reference.n_failed_fits},
                "decoys": _fit_totals(tables.get("decoy_runs")),
                "instability": _fit_totals(tables.get("instability_runs"))}
    inputs = input_file_digests(cfg)
    write_summary(out / "summary.md", decision, flags, tables, cfg, nmf_fits=nmf_fits, input_files=inputs,
                  notes=notes)
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

    rule = rule_text(cfg["decision"], **rule_settings(cfg), deviations=deviations)
    (out / "run_config.json").write_text(json.dumps(
        {"config": plain(cfg), "params": vars(params), "rule_version": decision["rule_version"],
         "registered_rule": decision["registered_rule"], "deviations": decision["deviations"], "rule": rule,
         "input_files": inputs, "nmf_fits": nmf_fits, "notes": notes, "versions": versions,
         "python": platform.python_version(),
         "nmf_cache_code": {"real": code_version(), "decoy": code_version(decoy_data=True)},
         "decision": decision},
        indent=2, default=str))
    log(f"[p1] decision: {decision['verdict']} ({decision['n_counted']} counted MPs; {decision['rule_version']}); "
        f"see {out / 'summary.md'}")
    return {"decision": decision, "tables": tables, "reference": reference, "mp_sets": mp_sets}
