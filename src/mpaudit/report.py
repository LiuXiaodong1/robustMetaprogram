"""P1 GO/NO-GO decision (frozen rule ``P1-GO v1.3``) and the markdown summary.

The rule is frozen before any real cohort is audited; ``docs/02_TESSERA_proposal.md`` and
``docs/05_pilot_P1.md`` state the same rule. Versions (none was registered or run on real data):

* v1.0 (2026-10-16): the first frozen rule.
* v1.1 (2026-10-17): v1.0 plus three requirements added after an adversarial review: the
  ``rho_hat`` x 1.5 condition in (b), the log 1.5 kappa margin in (c) and the >= 3
  distinct-tumour t prediction bound for the replicate null in (c).
* v1.2 (2026-10-17): v1.1 plus the kappa activity gate in (c) (a paired tumour is tested, and a
  replicate pair enters the null, only where the MP is active in at least one of the two pieces),
  and the registered settings outside the ``decision`` block are enforced
  (:data:`REGISTERED_SETTINGS`, :func:`rule_deviations`), as are the inputs without which the
  dosage and curated-artefact exclusions do not run (:func:`mpaudit.pipeline.config_deviations`).
* v1.3 (2026-10-18): v1.2 with criteria (c) and (a) recalibrated after a design comparison found
  that v1.2's activity gate selected tumours on the noisy kappa it then tested (toy null: 0.13-0.17
  of the tested (set, tumour) pairs without the program called losses, up to 0.31 of the tested
  tumours on exchangeable synthetic laws (0.05-0.31 across them), and an unplanted random set met
  (c)) and that
  ``kappa_retained`` admitted samples where the MP is absent (a minority biology MP reached 0.66).
  (c): no activity gate (every paired tumour with a defined change, null bound and self bound is
  tested; ``protocol.kappa_active`` is informational and no longer registered); the replicate bound
  is the lower of the linear-scale t bound carried to the MP's kappa and the log-scale t bound; the
  random-set bound is the ``floor(0.05 (n + 1))``-th order statistic, its sets drawn per tumour and
  MP; a loss must also beat 1.645 x the bootstrap SE of the MP's own change (self bound); every
  sample of a tumour is scored on the same number of cells, at the same depth (binomial thinning to
  the tumour's lowest median count per cell), with the panel taken on those counts; no draw depends
  on the other MPs or on the cohort's cell order. (a): ``kappa_retained`` is the share of the MP's
  co-expression in excess of 50 matched random gene sets (``decoy.n_null_floor`` 50) that the decoys
  keep, over the samples whose excess is above 3 SD of the sets. Every random draw differs from
  v1.2, so no v1.2 output carries over.

v1.3 is the version to register. Its calibration on the toy and on synthetic laws (toy results only,
not real data) is in ``scripts/calibration/README.md``.

It is built so that it can fail:

* **Retained MPs** are audited MPs not labelled as curated artefacts (label 1 in ``typer.labels``,
  e.g. the 11 MPs removed by hand in 3CA v2). Rediscovering curated artefacts is not news.
  (:func:`mpaudit.pipeline.curated_labels` carries the labels onto re-derived MPs by gene overlap,
  Jaccard >= ``min_jaccard``.) Without ``typer.labels`` no MP is excluded and ``summary.md`` says
  so (``curated_labels_given`` False); with ``mps: rederive`` that is a deviation from the
  registered rule, with MPs from a GMT the config check warns.
* **Pre-declared exclusions.**

  - *Expected TME-shared*: MPs named on ``expected_tme_shared``, or overlapping a gene set of
    ``expected_tme_shared_gmt`` / ``expected_tme_shared_sets`` by >= ``min_overlap`` genes, cannot
    count through criterion (b); they can still count through (a) or (c). IFN response or hypoxia
    are expected to co-vary with TME cells while being real. The list switches off (b) only, a
    deliberate departure from the 2026-10-16 review (which asked to drop listed MPs altogether):
    TME co-variation cannot tell niche-driven biology from leakage, but technical-decoy
    recruitment and a kappa loss without dissociation do not rest on TME co-variation, and an
    expected-TME-shared MP that they flag is the nuisance P1 looks for. Names only work for MPs
    read from a GMT; re-derived MPs (``mps: rederive``) are named by the run, so their exclusions
    are gene sets. With the TME axis on and nothing pre-declared, the config check warns.
  - *Dosage*: an MP whose genes are arm-enriched (``q_excluding_known_clusters`` <
    ``positional_q``), or that a dosage decoy arm (e.g. CNA clone) recruits in >= ``min_hits``
    decoys with ``kappa_retained`` >= ``min_kappa_retained``, is labelled dosage and never counted:
    it is genetic, not technical nuisance. A dosage hit needs no minimum number of decoys (a power
    requirement belongs to negative results; here it would switch off an exclusion).

* **Nuisance support**, any of:

  (a) a technical decoy arm recruits the MP (a decoy MP at Jaccard >= ``min_jaccard``, 0.3) in
      >= ``min_hits`` of its decoys (``decoy.n``, registered at exactly 5; >= ``min_decoys`` are
      needed) and ``kappa_retained`` >= ``min_kappa_retained``: the share of the MP's co-expression in
      excess of matched random gene sets that the decoys keep (:func:`mpaudit.decoy.kappa_retained`:
      per malignant sample of >= 20 cells, on the MP's genes detected in >= 5% of them (at least 5),
      ``(kappa_decoy - bg_decoy) / (kappa_real - bg_real)``, ``bg`` the median
      kappa of ``decoy.n_null_floor`` (50) random gene sets matched on mean expression and detection;
      the median over the samples where ``kappa_real - bg_real > max(decoy.min_kappa, 3 SD of the
      sets' kappa)`` (``decoy.min_kappa`` 0.02), averaged over the arm's decoys; NaN, so no hit, when
      no sample qualifies);
  (b) TME co-variation (:mod:`mpaudit.tme`), all of: the statistic is spillover-adjusted
      (``spill_adjusted`` True; with too few malignant-specific markers ``tme_covariation`` falls
      back to the unadjusted statistic, which cannot count); ``slope > 0``; ``partial_r2`` above
      the ``a1_null`` quantile (q99) of expression- and specificity-matched random gene sets;
      >= ``min_a1_samples`` samples; the MP is not expected TME-shared; and, where the
      ``rho_hat`` sensitivity re-run exists (``partial_r2_rho`` finite; with the registered
      ``tme.rho_sensitivity`` 1.5 it exists wherever the adjusted ``partial_r2`` is finite), the
      association also holds with every ``rho_hat`` scaled up: ``slope_rho > 0`` and ``partial_r2_rho`` above
      the same quantile of its own null (``null_q99_rho``). No fixed R^2 threshold;
  (c) protocol (:mod:`mpaudit.protocol`): >= ``min_paired_tumours`` tested paired tumours and a
      one-sided co-expression (kappa) loss in >= ``protocol_frac`` of them
      (``frac_loss_beyond_null_kappa``). Every paired tumour with a defined Δlog κ, null bound and
      self bound is tested; there is no activity gate (the activity columns are informational).
      κ is that of the MP's panel genes (detected in >= 5% of the tumour's κ cells under each
      protocol; undefined below 5 genes or a panel retention of 0.5), in samples of >= 20 malignant
      cells. Every sample of a tumour is scored on the same number of cells (the tumour's smallest
      sample, at most 2000), their counts binomially thinned to the tumour's lowest median count per
      cell. A tumour's loss is called by ``mpaudit.protocol``: Δlog κ below the null bound, below the
      self bound and below ``-kappa_min_effect`` (log 1.5). The null bound is, with same-protocol
      replicate pairs (κ finite in both pieces) from >= ``min_pairs`` (3) distinct tumours, the
      lower of two one-sided 95% t prediction bounds, on the linear κ changes (carried to the MP's κ
      level) and on the log κ changes; otherwise the k-th smallest linear κ change, ``k =
      floor(0.05 (n + 1))``, of the n >= 20 matched random gene sets of the tumour with a defined
      change (``protocol.n_null_sets`` 100 drawn per tumour and MP), carried to the MP's κ level; a
      tumour whose k-th smallest change is not below its k-th largest is not tested. The self bound
      is -1.645 x the bootstrap SE of the MP's own linear κ change (50 resamples of each sample's κ
      cells, drawn per seed, sample and MP), carried to its κ level; a tumour with an undefined SE
      is not tested.

* **Verdict.** GO if >= ``min_mps`` retained, non-dosage MPs have nuisance support; INCOMPLETE
  if no criterion is adequately powered; otherwise PIVOT. Powered: (a) some technical arm has
  >= ``min_decoys`` decoys; (b) some MP has >= ``min_a1_samples`` samples, a finite ``a1_null``
  and ``spill_adjusted`` True; (c) some MP has >= ``min_paired_tumours`` tested paired tumours (its
  κ change, null bound and self bound defined: ``n_paired_tumours`` of :mod:`mpaudit.protocol`; a
  metadata count of paired tumours is only an upper bound on it).

**Informational only**: instability, positional flags outside the dosage rule, the within-sample
kappa ratio (``kappa_ratio_median``; not spillover-adjusted), the kappa activity columns of
:mod:`mpaudit.protocol`, A4 intensity and prevalence, and the trivial typer.

**Registered settings.** The rule reads the ``decision`` block and settings outside it. The code
compares two groups with their registered values: the thresholds of the ``decision`` block
(``DEFAULT_RULE``; not the arm names or the pre-declared lists, which the user fills in, and not
``unstable_recovery``, which is informational) and :data:`REGISTERED_SETTINGS`, each of which must
equal its registered value (a null size is not an "at least" setting: another size moves the
estimated null quantile or background, so choosing it would give a second chance to pass).
:func:`mpaudit.pipeline.config_deviations` adds the inputs the exclusions need (curated labels with
their gene sets under ``mps: rederive``, a dosage decoy arm, a cytoband and gene positions, which
:func:`mpaudit.pipeline.run_p1` also checks place genes on chromosome arms).
:func:`rule_deviations` lists every difference of the first two groups;
:func:`mpaudit.pipeline.validate_config` refuses a config with any deviation unless its top-level
``unregistered_variant`` is true, and then the verdict is labelled
``P1-GO v1.3 (unregistered variant: ...)`` (:func:`rule_label`) in ``decision``, ``rule_text``,
``run_config.json`` and ``summary.md``. :func:`decide` finds the deviations of the ``decision``
block it is given and takes the others from its ``deviations`` argument; :func:`rule_text` also
compares the settings passed to it. The label does not cover what the code does not compare: the
contents of the decoy arms, the pre-declared lists, the label and gene-set files, the ``gavish``
(NMF) block other than ``gavish.seed``, ``counts_layer``, ``tme.covariates``, ``tme.n_boot``, the
keys, ``positional.chrom_col``/``start_col``, ``recurrence_key``, the protocol names and whether the
TME and protocol axes are switched on (off, a criterion is unpowered, which can turn PIVOT into
INCOMPLETE; :func:`mpaudit.pipeline.run_p1` logs a skipped axis and ``summary.md`` notes it), and
every other config key. These are registered as filled in, in the frozen config file, and a run is
checked against that file through ``run_config.json`` (the normalised config) and the input-file
digests; a registered label alone does not prove that the registered file was used.
Constants of the rule that are fixed in code, not settings (:data:`FIXED_CONSTANTS`: the κ floor
0.01, the self bound's 50 bootstrap resamples and z 1.645, the level and the minimum of 20 sets of
the random-set bound, the cell cap 2000 and the panel of (c) (samples of >= 20 cells, genes detected
in >= 5%, >= 5 genes and a retention >= 0.5), the 3-SD entry and the minimum of 20 random sets of
``kappa_retained`` and its gene set (samples of >= 20 cells, genes detected in >= 5%, >= 5 genes)),
are stated in :func:`rule_text` and pinned by the tests, which also check that the functions
applying them take them as defaults. ``protocol.kappa_active`` only sets informational columns
and is not part of the rule.

:func:`validate_rule` checks a ``decision`` block (unknown or pre-v1.0 keys, types, the GMT path)
and :func:`decide` calls it, as does :func:`mpaudit.pipeline.run_p1` before anything runs.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .data import read_gmt
from .decoy import EXCESS_ENTRY_SD, MIN_FLOOR_SETS, N_NULL_FLOOR
from .decoy import MIN_CELLS as DECOY_MIN_CELLS
from .decoy import MIN_DETECT as DECOY_MIN_DETECT
from .decoy import MIN_GENES as DECOY_MIN_GENES
from .protocol import (KAPPA_BOOT, KAPPA_FLOOR, KAPPA_MIN_EFFECT, KAPPA_SELF_Z, MAX_KAPPA_CELLS, MIN_CELLS, MIN_DETECT,
                       MIN_GENES, MIN_IOTA, MIN_PAIRS, MIN_RANDOM_SETS, NULL_ALPHA)

RULE_VERSION = "P1-GO v1.3"
RHO_SENSITIVITY = 1.5  # registered rho_hat factor of the sensitivity re-run (tme.rho_sensitivity)
MIN_JACCARD = 0.3  # registered MP match threshold (top-level min_jaccard): decoy hits, label carry-over

DEFAULT_RULE = {
    "technical_arms": ["technical"],
    "dosage_arms": ["dosage"],
    "min_decoys": 5,
    "min_hits": 2,
    "min_kappa_retained": 0.5,
    "a1_null": "null_q99",
    "min_a1_samples": 10,
    "min_paired_tumours": 4,
    "protocol_frac": 0.5,
    "positional_q": 0.05,
    "expected_tme_shared": [],
    "expected_tme_shared_gmt": None,
    "expected_tme_shared_sets": {},
    "min_overlap": 10,
    "min_mps": 3,
    "unstable_recovery": 0.5,
}
# Thresholds of the decision block that are registered (compared with DEFAULT_RULE); not the arm
# names or the pre-declared lists, which the user fills in, and not unstable_recovery (informational).
REGISTERED_DECISION_KEYS = ("min_decoys", "min_hits", "min_kappa_retained", "a1_null", "min_a1_samples",
                            "min_paired_tumours", "protocol_frac", "positional_q", "min_overlap", "min_mps")
# Settings outside the decision block that the rule or its axes' rule statistics read: (config path,
# registered value). Any other value deviates. Larger null sizes are not "at least" settings: another
# size moves the estimated null quantile or background (tme.n_null, decoy.n_null_floor and
# protocol.n_null_sets add or drop sets, drawn from a stream per gene set, per (sample, set) or per
# (tumour, MP)), so choosing it would give a second chance to pass.
REGISTERED_SETTINGS = (
    ("min_jaccard", MIN_JACCARD),                # (a) decoy hits; curated labels of re-derived MPs
    ("seed", 0),                                 # decoys, random-gene-set nulls, cell subsamples and thinning,
                                                 # bootstrap resamples
    ("gavish.seed", 0),                          # NMF restarts of the reference run and of every decoy fit
    ("decoy.enabled", True),                     # off: no (a) and no dosage-arm exclusion
    ("decoy.n", 5),                              # (a) is >= min_hits of these: more decoys, cheaper hits
    ("decoy.method", "residual"),
    ("decoy.n_depth_bins", 1),                   # depth strata leak real co-expression into decoys
    ("decoy.min_kappa", 0.02),                   # kappa_retained: minimum excess of an entering sample
    ("decoy.n_null_floor", N_NULL_FLOOR),        # kappa_retained: random sets behind bg and its SD
    ("tme.spillover_adjust", True),              # (b) is defined on the adjusted statistic
    ("tme.rho_sensitivity", RHO_SENSITIVITY),
    ("tme.null_match", "expression+specificity"),
    ("tme.n_null", 500),
    ("tme.min_cells", 20),                       # which (sample, cell type) groups enter A1
    ("tme.n_markers", 50),                       # the markers behind rho_hat
    ("tme.min_marker_fold", 8.0),
    ("protocol.kappa_min_effect", KAPPA_MIN_EFFECT),
    ("protocol.min_pairs", MIN_PAIRS),           # more: replicates fall back to the narrower random-set null
    ("protocol.n_null_sets", 100),
    ("positional.enabled", True),                # off: no arm-enrichment dosage exclusion
)
# Constants of the rule fixed in code (not config settings): (module.NAME, value). rule_text states them and
# tests/test_report_rule.py pins them; changing one is a new rule version.
FIXED_CONSTANTS = (
    ("protocol.KAPPA_FLOOR", KAPPA_FLOOR),          # κ floored before the log
    ("protocol.NULL_ALPHA", NULL_ALPHA),            # level of each one-sided null bound (t bound, order statistic)
    ("protocol.MIN_RANDOM_SETS", MIN_RANDOM_SETS),  # finite random-set changes for a tumour's random-set bound
    ("protocol.KAPPA_BOOT", KAPPA_BOOT),            # bootstrap resamples behind the self bound
    ("protocol.KAPPA_SELF_Z", KAPPA_SELF_Z),        # z of the self bound
    ("protocol.MAX_KAPPA_CELLS", MAX_KAPPA_CELLS),  # κ cells per sample at most (else the tumour's smallest sample)
    ("protocol.MIN_CELLS", MIN_CELLS),              # malignant cells a sample needs for (c)
    ("protocol.MIN_DETECT", MIN_DETECT),            # panel: detected in >= this share of cells under each protocol
    ("protocol.MIN_IOTA", MIN_IOTA),                # panel retention below which the MP's κ is undefined
    ("protocol.MIN_GENES", MIN_GENES),              # panel genes needed for κ
    ("decoy.EXCESS_ENTRY_SD", EXCESS_ENTRY_SD),     # kappa_retained entry: excess > this many SDs of the sets
    ("decoy.MIN_FLOOR_SETS", MIN_FLOOR_SETS),       # kappa_retained entry: random sets of finite kappa needed
    ("decoy.MIN_CELLS", DECOY_MIN_CELLS),           # kappa_retained: cells a sample needs
    ("decoy.MIN_DETECT", DECOY_MIN_DETECT),         # kappa_retained: the MP's genes detected in >= this share ...
    ("decoy.MIN_GENES", DECOY_MIN_GENES),           # ... at least this many of them
)
# keyword of rule_text -> the registered setting it states
_RULE_TEXT_SETTINGS = {"n_decoys": "decoy.n", "min_kappa": "decoy.min_kappa", "n_null_floor": "decoy.n_null_floor",
                       "min_jaccard": "min_jaccard", "rho_sensitivity": "tme.rho_sensitivity",
                       "kappa_min_effect": "protocol.kappa_min_effect", "min_pairs": "protocol.min_pairs",
                       "n_null_sets": "protocol.n_null_sets"}
_SHOWN = {"protocol.kappa_min_effect": f"{KAPPA_MIN_EFFECT:.6g} (log 1.5)"}
A1_NULLS = ("null_q95", "null_q99")
# keys of the rule before P1-GO v1.0 and what replaced them
LEGACY_BEFORE = "P1-GO v1.0"
LEGACY_KEYS = {"recruitment": "min_hits / min_decoys / min_kappa_retained", "partial_r2": "a1_null",
               "protocol_frac_beyond_null": "protocol_frac"}
_LIST_KEYS = ("technical_arms", "dosage_arms", "expected_tme_shared")
_INT_KEYS = {"min_decoys": 1, "min_hits": 1, "min_a1_samples": 1, "min_paired_tumours": 1, "min_overlap": 1,
             "min_mps": 1}  # key -> minimum
_FRACTION_KEYS = {"min_kappa_retained": (0.0, np.inf), "protocol_frac": (0.0, 1.0), "positional_q": (0.0, 1.0),
                  "unstable_recovery": (0.0, 1.0)}  # key -> closed range


def _number(key: str, v, kind: type):
    """``v`` as ``kind`` (numeric strings such as YAML's ``'1e-2'`` are accepted), else ValueError."""
    bad = ValueError(f"decision.{key} must be {'a whole number' if kind is int else 'a number'}, got {v!r}")
    if isinstance(v, (bool, np.bool_)) or v is None:
        raise bad
    if isinstance(v, str):
        try:
            v = float(v.strip())
        except ValueError:
            raise bad from None
    if not isinstance(v, (int, float, np.integer, np.floating)) or not np.isfinite(v):
        raise bad
    if kind is int:
        if float(v) != int(v):
            raise bad
        return int(v)
    return float(v)


def validate_rule(rule: dict | None) -> dict:
    """Check a ``decision`` block and return the full rule (defaults filled in, values normalised).

    Raises ``ValueError`` for keys not in ``DEFAULT_RULE`` (all of them listed; keys of the rule
    before ``P1-GO v1.0`` get a pointer to their replacements), for values of the wrong type or out
    of range, for an ``a1_null`` other than ``null_q95``/``null_q99``, for an arm named both
    technical and dosage, for ``expected_tme_shared_sets`` that is not a mapping of name to genes
    and for an ``expected_tme_shared_gmt`` that is not a file. A single name given for a list key
    (``technical_arms``, ``dosage_arms``, ``expected_tme_shared``; YAML ``key: MP0``) becomes a
    one-element list, and numeric strings become numbers.
    """
    rule = dict(rule or {})
    unknown = sorted(set(rule) - set(DEFAULT_RULE))
    if unknown:
        msg = f"unknown decision key(s) {unknown}; valid keys: {sorted(DEFAULT_RULE)}"
        legacy = [k for k in unknown if k in LEGACY_KEYS]
        if legacy:
            msg += (f". {', '.join(legacy)} belong(s) to the rule before {LEGACY_BEFORE}; {RULE_VERSION} replaces "
                    + "; ".join(f"{k} by {LEGACY_KEYS[k]}" for k in legacy)
                    + " (see docs/05_pilot_P1.md, 'The frozen GO rule')")
        raise ValueError(msg)
    r = {**DEFAULT_RULE, **rule}
    for k in _LIST_KEYS:
        v = r[k]
        if v is None:
            v = []
        elif isinstance(v, str):
            v = [v]
        if not isinstance(v, (list, tuple, set)) or not all(isinstance(x, str) and x for x in v):
            raise ValueError(f"decision.{k} must be a list of names, got {r[k]!r}")
        r[k] = list(dict.fromkeys(v))
    for k, lo in _INT_KEYS.items():
        r[k] = _number(k, r[k], int)
        if r[k] < lo:
            raise ValueError(f"decision.{k} must be >= {lo}, got {r[k]}")
    for k, (lo, hi) in _FRACTION_KEYS.items():
        r[k] = _number(k, r[k], float)
        if not lo <= r[k] <= hi:
            raise ValueError(f"decision.{k} must be in [{lo}, {hi}], got {r[k]}")
    if r["a1_null"] not in A1_NULLS:
        raise ValueError(f"decision.a1_null must be one of {A1_NULLS}, got {r['a1_null']!r}")
    both = sorted(set(r["technical_arms"]) & set(r["dosage_arms"]))
    if both:
        raise ValueError(f"decoy arm(s) {both} are listed as both technical_arms and dosage_arms")
    sets = r["expected_tme_shared_sets"]
    if sets is None:
        sets = {}
    if not isinstance(sets, dict) or not all(isinstance(g, (list, tuple, set)) for g in sets.values()):
        raise ValueError("decision.expected_tme_shared_sets must map a name to a list of genes")
    r["expected_tme_shared_sets"] = {str(k): [str(g) for g in v] for k, v in sets.items()}
    gmt = r["expected_tme_shared_gmt"]
    if gmt is not None and gmt != "":
        if not isinstance(gmt, (str, Path)) or not Path(gmt).is_file():
            raise ValueError(f"decision.expected_tme_shared_gmt {gmt!r} is not a file")
    else:
        r["expected_tme_shared_gmt"] = None
    return r


def _matches(value, registered) -> bool:
    """Whether a config value is the registered one (see :data:`REGISTERED_SETTINGS`)."""
    if isinstance(registered, bool):
        return isinstance(value, (bool, np.bool_)) and bool(value) is registered
    if isinstance(registered, str):
        return value == registered
    if isinstance(value, (bool, np.bool_)) or value is None:
        return False
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    if not np.isfinite(x):
        return False
    return bool(np.isclose(x, registered, rtol=1e-6, atol=1e-9))


def _value_text(v) -> str:
    return f"{v:.6g}" if isinstance(v, (float, np.floating)) else repr(v) if isinstance(v, str) else str(v)


def _shown(setting: str, registered) -> str:
    return _SHOWN.get(setting) or _value_text(registered)


def _decision_deviations(rule: dict) -> list[tuple[str, object, str]]:
    """Registered thresholds of a validated ``decision`` block that differ from ``DEFAULT_RULE``."""
    return [(f"decision.{k}", rule[k], _shown(f"decision.{k}", DEFAULT_RULE[k]))
            for k in REGISTERED_DECISION_KEYS if not _matches(rule[k], DEFAULT_RULE[k])]


def rule_deviations(cfg: dict) -> list[tuple[str, object, str]]:
    """Every setting the rule reads whose value differs from its registered value, as
    ``(setting, value, registered)`` with ``registered`` as text (e.g. ``'1.5'``, ``'500'``); an
    empty list for the registered rule.

    Covers the registered thresholds of the ``decision`` block (:data:`REGISTERED_DECISION_KEYS`,
    checked with :func:`validate_rule`) and :data:`REGISTERED_SETTINGS`. ``cfg`` is a merged config
    (:func:`mpaudit.pipeline.merge_config`; :func:`mpaudit.pipeline.validate_config` normalises its
    values first); a setting missing from it counts as registered.
    """
    out = _decision_deviations(validate_rule(cfg.get("decision")))
    for setting, registered in REGISTERED_SETTINGS:
        node = cfg
        for part in setting.split("."):
            node = node.get(part, registered) if isinstance(node, dict) else registered
        if not _matches(node, registered):
            out.append((setting, node, _shown(setting, registered)))
    return out


def _as_deviations(deviations) -> list[tuple[str, object, str]]:
    """``(setting, value, registered)`` tuples from tuples or from ``decision['deviations']`` dicts."""
    out = []
    for d in deviations or ():
        out.append((d["setting"], d["value"], d["registered"]) if isinstance(d, dict) else tuple(d))
    return out


def _all_deviations(rule: dict, deviations, settings: dict | None = None) -> list[tuple[str, object, str]]:
    """The validated ``rule``'s own deviations, then those of ``deviations``, then those of
    ``settings`` (keyword arguments of :func:`rule_text`), each only if not already listed."""
    out = _decision_deviations(rule)
    registered = dict(REGISTERED_SETTINGS)
    order = {name: i for i, (name, _) in enumerate(REGISTERED_SETTINGS)}
    own = sorted(((_RULE_TEXT_SETTINGS[k], v) for k, v in (settings or {}).items()), key=lambda x: order[x[0]])
    for d in _as_deviations(deviations) + [(s, v, _shown(s, registered[s])) for s, v in own
                                           if not _matches(v, registered[s])]:
        if d[0] not in {x[0] for x in out}:
            out.append(d)
    return out


def rule_label(deviations=()) -> str:
    """``RULE_VERSION``, or ``'P1-GO v1.3 (unregistered variant: setting=value, ...)'``."""
    devs = _as_deviations(deviations)
    if not devs:
        return RULE_VERSION
    return f"{RULE_VERSION} (unregistered variant: " + ", ".join(f"{s}={_value_text(v)}" for s, v, _ in devs) + ")"


def _deviation_list(devs) -> str:
    return "; ".join(f"{s} = {_value_text(v)} (registered: {r})" for s, v, r in devs)


def _expected_tme_shared(mp_sets: dict, rule: dict) -> dict[str, str]:
    """MP -> reason, for MPs on the pre-declared expected-TME-shared list."""
    out = {}
    names = set(rule.get("expected_tme_shared") or [])
    ref = dict(rule.get("expected_tme_shared_sets") or {})
    if rule.get("expected_tme_shared_gmt"):
        ref.update(read_gmt(rule["expected_tme_shared_gmt"]))
    for mp, genes in mp_sets.items():
        if mp in names:
            out[mp] = f"name '{mp}'"
            continue
        for ref_name, ref_genes in ref.items():
            if len(set(genes) & set(ref_genes)) >= rule["min_overlap"]:
                out[mp] = f"overlaps '{ref_name}'"
                break
    return out


def _arm_summary(dec: pd.DataFrame, arms: list[str], rule: dict, mps: list[str],
                 require_power: bool = True) -> pd.DataFrame:
    """Per MP: did any of ``arms`` recruit it under the frozen thresholds?

    A hit needs >= ``min_hits`` hits and ``kappa_retained`` >= ``min_kappa_retained`` and, with
    ``require_power`` (technical arms, criterion a), >= ``min_decoys`` decoys. ``powered`` says
    whether any of the arms ran >= ``min_decoys`` decoys. The reported arm is the first hit, else the
    arm with the most hits; ``kappa_retained_n_samples`` is NaN for a table without that column.
    """
    res = pd.DataFrame({"mp": mps, "hit": False, "best_arm": "", "hits": np.nan, "n_decoys": np.nan,
                        "kappa_retained": np.nan, "kappa_retained_n_samples": np.nan, "powered": False})
    if dec is None or dec.empty:
        return res
    missing = sorted({"mp", "arm", "n_decoys", "hits", "kappa_retained"} - set(dec.columns))
    if missing:  # silently unpowered would also drop the dosage-arm exclusion
        raise ValueError(f"decoy_recruitment table lacks column(s) {missing}")
    d = dec[dec["arm"].isin(arms)]
    if d.empty:
        return res
    res["powered"] = bool((d["n_decoys"] >= rule["min_decoys"]).any())
    row_of = {mp: i for i, mp in enumerate(mps)}
    for mp, rows in d.groupby("mp", sort=False):
        if mp not in row_of:
            continue
        ok = (rows["hits"] >= rule["min_hits"]) & (rows["kappa_retained"] >= rule["min_kappa_retained"])
        if require_power:
            ok &= rows["n_decoys"] >= rule["min_decoys"]
        best = rows.loc[ok].iloc[0] if ok.any() else rows.sort_values("hits", ascending=False).iloc[0]
        res.loc[row_of[mp], ["hit", "best_arm", "hits", "n_decoys", "kappa_retained", "kappa_retained_n_samples"]] = [
            bool(ok.any()), best["arm"], best["hits"], best["n_decoys"], best["kappa_retained"],
            best.get("kappa_retained_n_samples", np.nan)]
    res["hit"] = res["hit"].astype(bool)
    return res


def decide(mp_sets: dict, tables: dict[str, pd.DataFrame], rule: dict | None = None,
           curated_artefacts: pd.Series | None = None, deviations=None):
    """Apply ``P1-GO v1.3``. Returns ``(decision, flags)``.

    ``tables`` holds the axis outputs under the names :func:`mpaudit.pipeline.run_p1` uses
    (``decoy_recruitment``, ``tme_covariation``, ``protocol_shift``, ``positional``,
    ``instability_per_mp``); missing tables leave their criterion unpowered. A decoy table without
    one of the columns ``mp, arm, n_decoys, hits, kappa_retained`` raises ``ValueError`` (it also
    carries the dosage-arm exclusion, which must not vanish silently).
    ``curated_artefacts`` maps an audited MP name to its label (1 = curated artefact); None means
    no labels were given (``decision['curated_labels_given']`` False). ``deviations`` lists the
    deviations from the registered settings outside ``rule`` (:func:`rule_deviations` of the run's
    config); those of ``rule`` itself are found here. With any, ``decision['rule_version']`` reads
    ``P1-GO v1.3 (unregistered variant: ...)``, ``decision['registered_rule']`` is False and
    ``decision['deviations']`` lists them. Criterion (c) is powered when some MP has >=
    ``min_paired_tumours`` tested paired tumours (``n_paired_tumours`` of the ``protocol_shift``
    table: change, null bound and self bound defined; no activity gate).
    """
    rule = validate_rule(rule)
    devs = _all_deviations(rule, deviations)
    mps = list(mp_sets)
    flags = pd.DataFrame({"mp": mps})

    def attach(name, cols):
        nonlocal flags
        t = tables.get(name)
        present = [c for c in cols if t is not None and not t.empty and c in t]
        if present:
            flags = flags.merge(t[["mp", *present]], on="mp", how="left")
        for c in cols:
            if c not in flags:
                flags[c] = np.nan

    attach("tme_covariation", ["n_samples", "partial_r2", "slope", "null_q95", "null_q99", "null_p",
                               "spill_adjusted", "partial_r2_rho", "slope_rho", "null_q95_rho", "null_q99_rho",
                               "partial_r2_unadjusted", "kappa_ratio_median"])
    attach("protocol_shift", ["n_paired_tumours", "null_type", "n_null_tumours", "mean_delta_log_kappa",
                              "frac_loss_beyond_null_kappa", "frac_gain_beyond_null_kappa", "mean_delta_intensity"])
    attach("instability_per_mp", ["recovery_all"])
    attach("positional", ["best_arm", "best_arm_excluding_known_clusters", "q_excluding_known_clusters"])
    flags = flags.rename(columns={"best_arm": "positional_arm",
                                  "best_arm_excluding_known_clusters": "positional_arm_excluding_known"})

    dec = tables.get("decoy_recruitment")
    tech = _arm_summary(dec, rule["technical_arms"], rule, mps)
    dose = _arm_summary(dec, rule["dosage_arms"], rule, mps, require_power=False)
    ignored_arms = []
    if dec is not None and not dec.empty and "arm" in dec:
        ignored_arms = sorted(set(dec["arm"]) - set(rule["technical_arms"]) - set(rule["dosage_arms"]))
        if ignored_arms:
            warnings.warn(f"decoy arm(s) {ignored_arms} are in neither decision.technical_arms nor "
                          "decision.dosage_arms and are ignored by the rule", UserWarning, stacklevel=2)
    flags["tech_decoy_hits"] = tech["hits"].to_numpy()
    flags["tech_decoy_n"] = tech["n_decoys"].to_numpy()
    flags["tech_kappa_retained"] = tech["kappa_retained"].to_numpy()
    flags["tech_kappa_retained_n_samples"] = tech["kappa_retained_n_samples"].to_numpy()
    flags["dosage_decoy_hits"] = dose["hits"].to_numpy()
    flags["dosage_decoy_n"] = dose["n_decoys"].to_numpy()
    flags["dosage_kappa_retained"] = dose["kappa_retained"].to_numpy()
    flags["dosage_kappa_retained_n_samples"] = dose["kappa_retained_n_samples"].to_numpy()
    # MPs whose kappa_retained is undefined (no sample entered) in an arm of the rule: no hit possible there
    undefined = []
    if dec is not None and not dec.empty:
        rule_arms = dec[dec["arm"].isin(rule["technical_arms"] + rule["dosage_arms"]) & dec["mp"].isin(mps)]
        undefined = [f"{m} ({a})" for m, a in rule_arms.loc[rule_arms["kappa_retained"].isna(), ["mp", "arm"]]
                     .itertuples(index=False)]

    expected = _expected_tme_shared(mp_sets, rule)
    flags["expected_tme_shared"] = flags["mp"].map(lambda m: expected.get(m, ""))
    flags["curated_artefact"] = False
    if curated_artefacts is not None:
        lab = pd.to_numeric(pd.Series(curated_artefacts), errors="coerce")
        lab = lab[~lab.index.duplicated()]
        flags["curated_artefact"] = flags["mp"].map(lab).eq(1).to_numpy()

    positional = (flags["q_excluding_known_clusters"] < rule["positional_q"]).fillna(False).to_numpy(dtype=bool)
    dose_hit = dose["hit"].to_numpy(dtype=bool)
    flags["dosage"] = positional | dose_hit
    reasons = []
    for p, h, arm, darm in zip(positional, dose_hit, flags["positional_arm_excluding_known"], dose["best_arm"]):
        r = []
        if p:
            r.append(f"arm enrichment ({arm})" if isinstance(arm, str) and arm else "arm enrichment")
        if h:
            r.append(f"dosage decoy arm ({darm})")
        reasons.append("; ".join(r))
    flags["dosage_reason"] = reasons

    # criterion (b): spillover-adjusted, positive, above the null, robust to a larger rho_hat
    null = flags[rule["a1_null"]]
    null_rho = flags[rule["a1_null"] + "_rho"]
    adjusted = flags["spill_adjusted"].eq(True)
    enough = flags["n_samples"] >= rule["min_a1_samples"]
    above = (adjusted & enough & (flags["slope"] > 0) & (flags["partial_r2"] > null)).fillna(False)
    rho_run = flags["partial_r2_rho"].notna()
    rho_ok = ((flags["slope_rho"] > 0) & (flags["partial_r2_rho"] > null_rho)).fillna(False)
    flags["rho_sensitivity"] = np.where(~rho_run, "", np.where(rho_ok, "holds", "fails"))
    not_expected = flags["expected_tme_shared"] == ""
    flags["crit_b_tme"] = (above & not_expected & (~rho_run | rho_ok)).to_numpy(dtype=bool)
    a1_powered = bool((enough & null.notna() & adjusted).any())

    flags["crit_a_tech_decoy"] = tech["hit"].to_numpy(dtype=bool)
    flags["crit_c_protocol"] = ((flags["n_paired_tumours"] >= rule["min_paired_tumours"])
                                & (flags["frac_loss_beyond_null_kappa"] >= rule["protocol_frac"])
                                ).fillna(False).to_numpy(dtype=bool)
    # powered: some MP with >= min_paired_tumours tested paired tumours (n_paired_tumours counts the tested ones)
    protocol_powered = bool((flags["n_paired_tumours"] >= rule["min_paired_tumours"]).any())
    # informational: the gain mirror of (c), a count of MPs meeting (c) by chance when gains are as likely as losses
    flags["info_gain_mirror_c"] = ((flags["n_paired_tumours"] >= rule["min_paired_tumours"])
                                   & (flags["frac_gain_beyond_null_kappa"] >= rule["protocol_frac"])
                                   ).fillna(False).to_numpy(dtype=bool)

    flags["nuisance_supported"] = flags[["crit_a_tech_decoy", "crit_b_tme", "crit_c_protocol"]].any(axis=1)
    flags["counted"] = flags["nuisance_supported"] & ~flags["curated_artefact"] & ~flags["dosage"]
    flags["info_unstable"] = (flags["recovery_all"] < rule["unstable_recovery"]).fillna(False)

    powered = {"technical_decoys": bool(tech["powered"].any()), "tme": a1_powered, "protocol": protocol_powered}
    n = int(flags["counted"].sum())
    if n >= rule["min_mps"]:
        verdict = "GO"
    elif any(powered.values()):
        verdict = "PIVOT"
    else:
        verdict = "INCOMPLETE"
    decision = {
        "rule_version": rule_label(devs), "registered_rule": not devs,
        "deviations": [{"setting": s, "value": v, "registered": r} for s, v, r in devs],
        "verdict": verdict, "n_counted": n, "min_mps": rule["min_mps"],
        "powered": powered, "counted_mps": flags.loc[flags["counted"], "mp"].tolist(),
        "dosage_mps": flags.loc[flags["dosage"], "mp"].tolist(),
        "expected_tme_shared_mps": [m for m in mps if m in expected],
        "curated_labels_given": curated_artefacts is not None,
        "curated_artefact_mps": flags.loc[flags["curated_artefact"], "mp"].tolist(),
        "a1_unadjusted_mps": flags.loc[flags["spill_adjusted"].eq(False), "mp"].tolist(),
        "a1_rho_sensitivity_run": bool(rho_run.any()),
        "rho_sensitive_mps": flags.loc[above & not_expected & rho_run & ~rho_ok, "mp"].tolist(),
        "ignored_decoy_arms": ignored_arms,
        "kappa_retained_undefined": undefined,
        "gain_mirror_mps": flags.loc[flags["info_gain_mirror_c"], "mp"].tolist(),
    }
    return decision, flags


def _md_table(df: pd.DataFrame, max_rows: int | None = 60, csv: str | None = None) -> str:
    """Markdown table; beyond ``max_rows`` rows (None: no limit) it is cut, with a note naming ``csv``."""
    if df is None or df.empty:
        return "_(empty)_\n"
    n = len(df)
    if max_rows is not None:
        df = df.head(max_rows)

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else f"{v:.3g}"
        return str(v)

    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    if len(df) < n:
        lines.append(f"\n_(first {len(df)} of {n} rows; all rows in {csv or 'the CSV next to this file'})_")
    return "\n".join(lines) + "\n"


def rule_settings(cfg: dict) -> dict:
    """The settings outside the ``decision`` block that :func:`rule_text` states, from a merged
    config (keyword arguments of :func:`rule_text`; a missing setting takes its registered value)."""
    dc, tme, prc = cfg.get("decoy") or {}, cfg.get("tme") or {}, cfg.get("protocol") or {}
    rho = tme.get("rho_sensitivity", RHO_SENSITIVITY)
    return {"n_decoys": int(dc.get("n", 5)),
            "min_kappa": float(dc.get("min_kappa", 0.02)),
            "n_null_floor": int(dc.get("n_null_floor", N_NULL_FLOOR)),
            "kappa_min_effect": float(prc.get("kappa_min_effect", KAPPA_MIN_EFFECT)),
            "rho_sensitivity": None if rho is None else float(rho),
            "min_pairs": int(prc.get("min_pairs", MIN_PAIRS)),
            "n_null_sets": int(prc.get("n_null_sets", 100)),
            "min_jaccard": float(cfg.get("min_jaccard", MIN_JACCARD))}


def rule_text(rule: dict | None = None, *, n_decoys: int = 5, min_kappa: float = 0.02,
              n_null_floor: int = N_NULL_FLOOR, kappa_min_effect: float = KAPPA_MIN_EFFECT,
              rho_sensitivity: float | None = RHO_SENSITIVITY, min_pairs: int = MIN_PAIRS, n_null_sets: int = 100,
              min_jaccard: float = MIN_JACCARD, deviations=None) -> str:
    """The rule as applied, in one paragraph (for ``summary.md`` and ``run_config.json``).

    ``n_decoys``, ``min_kappa``, ``n_null_floor`` (``decoy``), ``kappa_min_effect``, ``min_pairs``,
    ``n_null_sets`` (``protocol``), ``rho_sensitivity`` (``tme``) and ``min_jaccard`` are the settings
    the axes ran with (:func:`rule_settings` reads them from a config); they are stated, not applied,
    here, and any that differs from its registered value (:data:`REGISTERED_SETTINGS`) is a deviation.
    The constants fixed in code (:data:`FIXED_CONSTANTS`) are stated as they are. ``deviations`` as for
    :func:`decide`: with any, the paragraph names the rule an unregistered variant and lists them."""
    r = validate_rule(rule)
    devs = _all_deviations(r, deviations, {"n_decoys": n_decoys, "min_kappa": min_kappa, "n_null_floor": n_null_floor,
                                           "kappa_min_effect": kappa_min_effect, "rho_sensitivity": rho_sensitivity,
                                           "min_pairs": min_pairs, "n_null_sets": n_null_sets,
                                           "min_jaccard": min_jaccard})
    a1 = r["a1_null"]
    if rho_sensitivity is not None and rho_sensitivity != 1:
        rho = (f", which also holds with every rho_hat multiplied by {rho_sensitivity:g} (slope_rho > 0 and "
               f"partial_r2_rho > {a1}_rho, where the re-run exists)")
    else:
        rho = " (the rho_hat sensitivity re-run is off, so the estimated rho_hat alone decides)"
    if kappa_min_effect > 0:
        margin = (f" and by more than {kappa_min_effect:.3g} in log kappa (a {np.exp(kappa_min_effect):.3g}-fold "
                  "loss)")
    else:
        margin = " (no minimum effect)"
    level = 1 - NULL_ALPHA
    text = (
        f"Rule {rule_label(devs)}: GO if >= {r['min_mps']} retained, non-dosage MPs have nuisance support: "
        f"(a) a technical decoy arm ({', '.join(r['technical_arms']) or 'none'}) recruits the MP (a decoy MP at "
        f"Jaccard >= {min_jaccard:g} to it) in >= {r['min_hits']} of its {n_decoys} decoys (>= {r['min_decoys']} "
        f"needed) with kappa_retained >= {r['min_kappa_retained']:g}, kappa_retained being the share of the MP's "
        f"co-expression in excess of matched random gene sets that the decoys keep: per malignant sample of >= "
        f"{DECOY_MIN_CELLS} cells, on the MP's genes detected in >= {DECOY_MIN_DETECT:.0%} of its cells (at least "
        f"{DECOY_MIN_GENES}), (kappa_decoy - bg_decoy) / (kappa_real - bg_real), bg the median kappa of "
        f"{n_null_floor} random gene sets matched on mean expression and detection, its median over the samples "
        f"where kappa_real - bg_real > "
        f"max({min_kappa:g}, {EXCESS_ENTRY_SD:g} SD of the sets' kappa) with >= {MIN_FLOOR_SETS} sets of finite "
        f"kappa, averaged over the arm's decoys (undefined, so no hit, when no sample qualifies); or (b) a "
        f"spillover-adjusted (spill_adjusted) TME association with slope > 0 and partial R^2 above the "
        f"random-gene-set {a1} in >= {r['min_a1_samples']} samples{rho}, for MPs not pre-declared as expected "
        f"TME-shared (named, or overlapping a pre-declared gene set by >= {r['min_overlap']} genes); or (c) a "
        f"one-sided kappa loss under the dissociation-free protocol in >= {r['protocol_frac']:g} of >= "
        f"{r['min_paired_tumours']} tested paired tumours, a tumour's loss being a change in log kappa (test - "
        f"reference; kappa floored at {KAPPA_FLOOR:g}) below the null bound, below the self bound{margin}. Kappa "
        f"is that of the MP's genes in the tumour's panel (detected in >= {MIN_DETECT:.0%} of its malignant cells "
        f"under each protocol; undefined with fewer than {MIN_GENES} panel genes or less than {MIN_IOTA:g} of the "
        f"MP's genes in the panel), in samples of >= {MIN_CELLS} malignant cells. Every paired tumour with a "
        f"defined change, null bound and self bound is tested (no activity gate); every sample of a tumour is "
        f"scored on the same number of cells (the tumour's smallest sample, at most {MAX_KAPPA_CELLS}), their "
        f"counts binomially thinned to the tumour's lowest median count per cell. The null bound is, with "
        f"same-protocol replicate pairs from >= {min_pairs} distinct tumours, the lower of two one-sided "
        f"{level:.0%} t prediction bounds, on the linear kappa changes (carried to the MP's kappa level) and on the "
        f"log kappa changes; otherwise the k-th smallest linear kappa change, k = floor({NULL_ALPHA:g} (n + 1)), of "
        f"the n >= {MIN_RANDOM_SETS} random gene sets with a defined change ({n_null_sets} drawn per tumour and MP, "
        f"matched on mean expression and detection), carried to the MP's kappa level, a tumour whose k-th smallest "
        f"change is not below its k-th largest not being tested. The self bound is -{KAPPA_SELF_Z:g} times the "
        f"bootstrap SE ({KAPPA_BOOT} resamples of each sample's cells, drawn per seed, sample and MP) of the MP's "
        f"own linear kappa change, carried to its kappa level; a tumour with an undefined SE is not tested. Curated "
        f"artefacts (labelled, or carried onto re-derived MPs at Jaccard >= {min_jaccard:g}) are not counted; "
        f"dosage MPs (arm enrichment q, known positional gene clusters excluded, < {r['positional_q']:g}, or a "
        f"dosage decoy arm ({', '.join(r['dosage_arms']) or 'none'}) with >= {r['min_hits']} hits and "
        f"kappa_retained >= {r['min_kappa_retained']:g}, whatever the number of decoys) are labelled, not "
        f"counted. INCOMPLETE if no criterion is powered: (a) >= {r['min_decoys']} decoys in some technical arm; "
        f"(b) for some MP, >= {r['min_a1_samples']} samples with TME cells, a finite null and the spillover "
        f"adjusted; (c) for some MP, >= {r['min_paired_tumours']} tested paired tumours. Otherwise PIVOT. "
        f"Instability, positional flags outside the dosage rule, the within-sample kappa ratio "
        f"(kappa_ratio_median, not spillover-adjusted), the kappa activity columns, A4 intensity and prevalence "
        f"and the trivial typer are informational."
    )
    if devs:
        text += f" This is NOT the registered rule {RULE_VERSION}; deviations: {_deviation_list(devs)}."
    return text


def _nmf_fits_line(nmf_fits: dict) -> str:
    """One line of capped and failed NMF fits per kind of run (``run_config.json`` ``nmf_fits``)."""
    parts = []
    for key, name in (("reference", "reference"), ("decoys", "decoy runs"), ("instability", "instability runs")):
        f = nmf_fits.get(key)
        if f:
            runs = f" ({f['n_runs']})" if "n_runs" in f else ""
            parts.append(f"{name}{runs} {f['n_capped_fits']}/{f['n_fits']}, {f['n_failed_fits']}")
    if not parts:
        return ""
    return "NMF fits stopped at max_iter (of all fits) / by the restart limit: " + "; ".join(parts) + ".\n"


def write_summary(path: Path, decision: dict, flags: pd.DataFrame, tables: dict, cfg: dict,
                  nmf_fits: dict | None = None, input_files: dict | None = None,
                  notes: list[str] | None = None) -> None:
    """Write ``summary.md``: the verdict (opening with a bold warning for an unregistered
    variant), the rule as applied (:func:`rule_text`), exclusions, the NMF fit diagnostics
    (``nmf_fits`` as in ``run_config.json``), the input files' sha256 (``input_files``, setting ->
    ``{"path", "sha256"}``), notes (the run's own ``notes``, e.g. axes that did not run, first) and
    the per-MP and per-run tables."""
    rule = validate_rule(cfg.get("decision"))
    devs = _all_deviations(rule, decision.get("deviations"))
    powered = ", ".join(k for k, v in decision["powered"].items() if v) or "none"
    counted = ", ".join(decision["counted_mps"]) or "none"
    if decision.get("curated_labels_given", True):
        curated = ", ".join(decision["curated_artefact_mps"]) or "none"
    else:
        curated = "not labelled (typer.labels not given; no MP excluded)"
    parts = ["# P1 audit summary\n"]
    if devs:
        parts.append(f"**This verdict is NOT the registered rule {RULE_VERSION}: the run is an unregistered variant "
                     f"(unregistered_variant: true). Deviations: {_deviation_list(devs)}.**\n")
    parts += [
        f"**Verdict: {decision['verdict']}** ({decision['rule_version']}): {decision['n_counted']} retained, "
        f"non-dosage MP(s) with nuisance support (GO needs >= {decision['min_mps']}): {counted}. "
        f"Powered criteria: {powered}.\n",
        rule_text(rule, **rule_settings(cfg), deviations=devs) + "\n",
        f"Dosage MPs (labelled, not counted): {', '.join(decision['dosage_mps']) or 'none'}. "
        f"Expected TME-shared (criterion b disabled): {', '.join(decision['expected_tme_shared_mps']) or 'none'}. "
        f"Curated artefacts: {curated}.\n",
    ]
    if nmf_fits and _nmf_fits_line(nmf_fits):
        parts.append(_nmf_fits_line(nmf_fits))
    if input_files:
        parts.append("Input files (sha256, first 16 hex digits; full digests in run_config.json): " + "; ".join(
            f"{k} {Path(v['path']).name} {v['sha256'][:16]}" for k, v in input_files.items()) + ".\n")
    notes = list(notes or [])
    unadjusted = decision.get("a1_unadjusted_mps") or []
    if unadjusted:
        notes.append(f"Spillover NOT adjusted for {len(unadjusted)} MP(s) (too few malignant-specific markers, or "
                     f"tme.spillover_adjust off): their partial R^2 is the unadjusted statistic and cannot count "
                     f"through (b): {', '.join(unadjusted)}.")
    if decision.get("rho_sensitive_mps"):
        notes.append("Above the null at the estimated rho_hat but not at the scaled rho_hat (consistent with an "
                     f"under-subtracted leak; criterion b not met): {', '.join(decision['rho_sensitive_mps'])}.")
    if decision["powered"].get("tme") and not decision.get("a1_rho_sensitivity_run", True):
        notes.append("The rho_hat sensitivity re-run did not run (tme.rho_sensitivity null or 1), so criterion (b) "
                     "used the estimated rho_hat only.")
    if decision.get("ignored_decoy_arms"):
        notes.append(f"Decoy arm(s) ignored by the rule (in neither technical_arms nor dosage_arms): "
                     f"{', '.join(decision['ignored_decoy_arms'])}.")
    if decision.get("kappa_retained_undefined"):
        notes.append("kappa_retained undefined (no sample's excess co-expression passed the entry rule), so no "
                     "criterion (a) hit or dosage typing in that arm: "
                     f"{', '.join(decision['kappa_retained_undefined'])}.")
    if decision["powered"].get("protocol") and "gain_mirror_mps" in decision:
        mirror = decision["gain_mirror_mps"]
        notes.append(f"Gain mirror of criterion (c) (informational null count: a kappa gain in >= "
                     f"{rule['protocol_frac']:g} of >= {rule['min_paired_tumours']} tested paired tumours; with no "
                     f"protocol effect gains are as likely as losses): {len(mirror)} MP(s)"
                     + (f": {', '.join(mirror)}." if mirror else "."))
    if "kappa_ratio_median" in flags and flags["kappa_ratio_median"].notna().any():
        notes.append("kappa_ratio_median (within-sample TME/malignant kappa ratio) is not spillover-adjusted and is "
                     "not used by the rule.")
    if notes:
        parts += ["**Notes.** " + " ".join(notes) + "\n"]
    a1 = rule["a1_null"]
    cols = ["mp", "tech_decoy_hits", "tech_decoy_n", "tech_kappa_retained", "tech_kappa_retained_n_samples",
            "dosage_decoy_hits", "spill_adjusted",
            "partial_r2", a1, "slope", "partial_r2_rho", f"{a1}_rho", "partial_r2_unadjusted", "kappa_ratio_median",
            "n_paired_tumours", "frac_loss_beyond_null_kappa", "frac_gain_beyond_null_kappa", "positional_arm",
            "positional_arm_excluding_known",
            "q_excluding_known_clusters", "recovery_all", "expected_tme_shared", "dosage", "curated_artefact",
            "crit_a_tech_decoy", "crit_b_tme", "crit_c_protocol", "counted"]
    parts += ["## Per-MP flags\n", _md_table(flags[[c for c in dict.fromkeys(cols) if c in flags]], max_rows=None)]
    if "instability_runs" in tables:
        s = tables["instability_runs"]
        parts += ["## Instability of the reference MPs\n",
                  _md_table(s.groupby("kind")[["recall", "precision", "mean_jaccard_matched"]].mean().reset_index()),
                  "\nPer perturbation:\n", _md_table(s, csv="instability_runs.csv")]
    if "decoy_runs" in tables:
        parts += ["## Decoy cohorts (nuisance only)\n", _md_table(tables["decoy_runs"], csv="decoy_runs.csv")]
    if "typer_eval" in tables:
        parts += ["## Trivial typer vs curated labels\n", _md_table(tables["typer_eval"], csv="typer_eval.csv")]
    parts.append("\nAll tables are written as CSV next to this file.\n")
    Path(path).write_text("\n".join(parts))
