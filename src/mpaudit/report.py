"""P1 GO/NO-GO decision (frozen rule ``P1-GO v1.0``) and the markdown summary.

The rule is frozen before any real cohort is audited; ``docs/02_TESSERA_proposal.md`` and
``docs/05_pilot_P1.md`` state the same rule. It is built so that it can fail:

* **Retained MPs** are audited MPs not labelled as curated artefacts (e.g. the 11 MPs removed by
  hand in 3CA v2). Rediscovering curated artefacts is not news.
* **Pre-declared exclusions.** MPs on the expected-TME-shared list (by name, or by an overlap of
  >= ``min_overlap`` genes with a pre-declared gene set) cannot count through the TME axis: IFN response or hypoxia are expected
  to co-vary with TME cells while being real. MPs typed **dosage** (arm-enriched genes, or
  recruited by a dosage decoy arm such as CNA clone) are labelled and never counted: they are
  genetic, not technical nuisance.
* **Nuisance support**, any of:
  (a) a technical decoy arm recruits the MP in >= ``min_hits`` of >= ``min_decoys`` decoys and
      the decoys keep >= ``min_kappa_retained`` of its co-expression;
  (b) a positive, spillover-adjusted TME association whose partial R^2 exceeds the ``a1_null``
      quantile (q99) of expression- and specificity-matched random gene sets, for an MP not on the
      expected-TME-shared list (no fixed R^2 threshold);
  (c) a one-sided co-expression (kappa) loss under the dissociation-free protocol, beyond the
      null, in >= ``protocol_frac`` of >= ``min_paired_tumours`` paired tumours.
* **Verdict.** GO if >= ``min_mps`` retained, non-dosage MPs have nuisance support; INCOMPLETE
  if no criterion is adequately powered; otherwise PIVOT.

Instability, positional flags outside the dosage rule and the trivial typer are informational.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .data import read_gmt

RULE_VERSION = "P1-GO v1.0"

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


def _arm_summary(dec: pd.DataFrame, arms: list[str], rule: dict, mps: list[str]) -> pd.DataFrame:
    """Per MP: did any of ``arms`` recruit it under the frozen thresholds?"""
    res = pd.DataFrame({"mp": mps, "hit": False, "best_arm": "", "hits": np.nan, "n_decoys": np.nan,
                        "kappa_retained": np.nan, "powered": False})
    if dec is None or dec.empty or "arm" not in dec:
        return res
    d = dec[dec["arm"].isin(arms)]
    if d.empty:
        return res
    res["powered"] = bool((d["n_decoys"] >= rule["min_decoys"]).any())
    for i, mp in enumerate(mps):
        rows = d[d["mp"] == mp]
        if rows.empty:
            continue
        ok = ((rows["n_decoys"] >= rule["min_decoys"]) & (rows["hits"] >= rule["min_hits"])
              & (rows["kappa_retained"] >= rule["min_kappa_retained"]))
        best = rows.loc[ok].iloc[0] if ok.any() else rows.sort_values("hits", ascending=False).iloc[0]
        res.loc[i, ["hit", "best_arm", "hits", "n_decoys", "kappa_retained"]] = [
            bool(ok.any()), best["arm"], best["hits"], best["n_decoys"], best["kappa_retained"]]
    return res


def decide(mp_sets: dict, tables: dict[str, pd.DataFrame], rule: dict | None = None,
           curated_artefacts: pd.Series | None = None):
    """Apply ``P1-GO v1.0``. Returns ``(decision, flags)``."""
    rule = {**DEFAULT_RULE, **(rule or {})}
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

    attach("tme_covariation", ["partial_r2", "slope", "null_q95", "null_q99", "null_p", "n_samples",
                               "partial_r2_unadjusted", "kappa_ratio_median"])
    attach("protocol_shift", ["n_paired_tumours", "null_type", "mean_delta_log_kappa",
                              "frac_loss_beyond_null_kappa", "mean_delta_intensity"])
    attach("instability_per_mp", ["recovery_all"])
    attach("positional", ["best_arm", "q_excluding_known_clusters"])
    flags = flags.rename(columns={"best_arm": "positional_arm"})

    dec = tables.get("decoy_recruitment")
    tech = _arm_summary(dec, rule["technical_arms"], rule, mps)
    dose = _arm_summary(dec, rule["dosage_arms"], rule, mps)
    flags["tech_decoy_hits"] = tech["hits"].to_numpy()
    flags["tech_decoy_n"] = tech["n_decoys"].to_numpy()
    flags["tech_kappa_retained"] = tech["kappa_retained"].to_numpy()
    flags["dosage_decoy_hits"] = dose["hits"].to_numpy()

    expected = _expected_tme_shared(mp_sets, rule)
    flags["expected_tme_shared"] = flags["mp"].map(lambda m: expected.get(m, ""))
    flags["curated_artefact"] = False
    if curated_artefacts is not None:
        flags["curated_artefact"] = flags["mp"].map(curated_artefacts).fillna(0).astype(bool)

    positional = (flags["q_excluding_known_clusters"] < rule["positional_q"]).fillna(False)
    flags["dosage"] = positional.to_numpy() | dose["hit"].to_numpy()
    flags["dosage_reason"] = np.where(positional, "arm enrichment", np.where(dose["hit"], "dosage decoy arm", ""))

    null = flags[rule["a1_null"]] if rule["a1_null"] in flags else pd.Series(np.nan, index=flags.index)
    a1_powered = bool((flags["n_samples"] >= rule["min_a1_samples"]).any())
    flags["crit_a_tech_decoy"] = tech["hit"].to_numpy()
    flags["crit_b_tme"] = ((flags["slope"] > 0) & (flags["partial_r2"] > null)
                           & (flags["n_samples"] >= rule["min_a1_samples"])
                           & (flags["expected_tme_shared"] == "")).fillna(False)
    flags["crit_c_protocol"] = ((flags["n_paired_tumours"] >= rule["min_paired_tumours"])
                                & (flags["frac_loss_beyond_null_kappa"] >= rule["protocol_frac"])).fillna(False)
    protocol_powered = bool((flags["n_paired_tumours"] >= rule["min_paired_tumours"]).any())

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
        "rule_version": RULE_VERSION, "verdict": verdict, "n_counted": n, "min_mps": rule["min_mps"],
        "powered": powered, "counted_mps": flags.loc[flags["counted"], "mp"].tolist(),
        "dosage_mps": flags.loc[flags["dosage"], "mp"].tolist(),
        "expected_tme_shared_mps": [m for m in mps if m in expected],
        "curated_artefact_mps": flags.loc[flags["curated_artefact"], "mp"].tolist(),
    }
    return decision, flags


def _md_table(df: pd.DataFrame, max_rows: int = 60) -> str:
    if df is None or df.empty:
        return "_(empty)_\n"
    df = df.head(max_rows)

    def fmt(v):
        if isinstance(v, (float, np.floating)):
            return "" if np.isnan(v) else f"{v:.3g}"
        return str(v)

    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines) + "\n"


def rule_text(rule: dict) -> str:
    r = {**DEFAULT_RULE, **rule}
    return (
        f"Rule {RULE_VERSION}: GO if >= {r['min_mps']} retained, non-dosage MPs have nuisance support: "
        f"(a) a technical decoy arm ({', '.join(r['technical_arms'])}) recruits the MP in >= {r['min_hits']} of "
        f">= {r['min_decoys']} decoys with kappa_retained >= {r['min_kappa_retained']}; or (b) a positive, "
        f"spillover-adjusted TME association above the random-gene-set {r['a1_null']} (>= {r['min_a1_samples']} "
        f"samples), for MPs not pre-declared as expected TME-shared; or (c) a kappa loss beyond the null under the "
        f"dissociation-free protocol in >= {r['protocol_frac']} of >= {r['min_paired_tumours']} paired tumours. "
        "Curated artefacts are not counted; dosage MPs (arm enrichment q < "
        f"{r['positional_q']} or a dosage decoy arm hit) are labelled, not counted. INCOMPLETE if no criterion is "
        "powered; otherwise PIVOT. Instability, other positional flags and the trivial typer are informational."
    )


def write_summary(path: Path, decision: dict, flags: pd.DataFrame, tables: dict, cfg: dict) -> None:
    rule = cfg.get("decision", {})
    powered = ", ".join(k for k, v in decision["powered"].items() if v) or "none"
    parts = [
        "# P1 audit summary\n",
        f"**Verdict: {decision['verdict']}** ({decision['rule_version']}): {decision['n_counted']} retained, "
        f"non-dosage MP(s) with nuisance support (GO needs >= {decision['min_mps']}). Powered criteria: {powered}.\n",
        rule_text(rule) + "\n",
        f"Dosage MPs (labelled, not counted): {', '.join(decision['dosage_mps']) or 'none'}. "
        f"Expected TME-shared (criterion b disabled): {', '.join(decision['expected_tme_shared_mps']) or 'none'}. "
        f"Curated artefacts: {', '.join(decision['curated_artefact_mps']) or 'none'}.\n",
        "## Per-MP flags\n",
        _md_table(flags[[c for c in ["mp", "tech_decoy_hits", "tech_decoy_n", "tech_kappa_retained",
                                     "dosage_decoy_hits", "partial_r2", "null_q99", "slope",
                                     "partial_r2_unadjusted", "kappa_ratio_median", "n_paired_tumours",
                                     "frac_loss_beyond_null_kappa", "positional_arm",
                                     "q_excluding_known_clusters", "recovery_all", "expected_tme_shared",
                                     "dosage", "curated_artefact", "crit_a_tech_decoy", "crit_b_tme",
                                     "crit_c_protocol", "counted"] if c in flags]]),
    ]
    if "instability_runs" in tables:
        s = tables["instability_runs"]
        parts += ["## Instability of the reference MPs\n",
                  _md_table(s.groupby("kind")[["recall", "precision", "mean_jaccard_matched"]].mean().reset_index()),
                  "\nPer perturbation:\n", _md_table(s)]
    if "decoy_runs" in tables:
        parts += ["## Decoy cohorts (nuisance only)\n", _md_table(tables["decoy_runs"])]
    if "typer_eval" in tables:
        parts += ["## Trivial typer vs curated labels\n", _md_table(tables["typer_eval"])]
    parts.append("\nAll tables are written as CSV next to this file.\n")
    Path(path).write_text("\n".join(parts))
