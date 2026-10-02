"""P1 GO/NO-GO decision and the markdown summary.

Decision rule (proposal, pilot P1): GO if at least ``min_mps`` MPs that curators *retained*
show at least one of: decoy recruitment >= ``recruitment``; a positive TME association with
partial R^2 >= ``partial_r2`` that also exceeds the 95th percentile of expression-matched random
gene sets;
a loss under the dissociation-free protocol beyond the same-protocol null in at least
``protocol_frac_beyond_null`` of paired tumours. MPs labelled as curated artefacts (e.g. the
3CA v2 hand-removed MPs) do not count towards GO: rediscovering them is not news.
Instability and positional enrichment are reported but are not part of the rule.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def decide(mp_sets: dict, tables: dict[str, pd.DataFrame], rule: dict, curated_artefacts: pd.Series | None = None):
    flags = pd.DataFrame({"mp": list(mp_sets)})

    def attach(name, cols):
        nonlocal flags
        if name in tables and not tables[name].empty:
            flags = flags.merge(tables[name][["mp", *cols]], on="mp", how="left")
        else:
            for c in cols:
                flags[c] = np.nan

    attach("decoy_recruitment", ["recruitment"])
    attach("tme_covariation", ["partial_r2", "partial_r2_lo", "slope", "null_q95", "kappa_ratio_median"])
    attach("protocol_shift", ["mean_delta_intensity", "frac_loss_beyond_null_intensity"])
    attach("instability_per_mp", ["recovery_all"])
    attach("positional", ["best_arm", "q", "q_excluding_known_clusters"])
    flags = flags.rename(columns={"q": "positional_q"})

    flags["flag_decoy"] = flags["recruitment"] >= rule["recruitment"]
    flags["flag_tme"] = ((flags["partial_r2"] >= rule["partial_r2"]) & (flags["slope"] > 0)
                         & (flags["partial_r2"] > flags["null_q95"].fillna(-np.inf)))
    flags["flag_protocol_loss"] = flags["frac_loss_beyond_null_intensity"] >= rule["protocol_frac_beyond_null"]
    flags["info_unstable"] = flags["recovery_all"] < rule["unstable_recovery"]
    flags["info_positional"] = flags["q_excluding_known_clusters"] < rule["positional_q"]
    flags["curated_artefact"] = False
    if curated_artefacts is not None:
        flags["curated_artefact"] = flags["mp"].map(curated_artefacts).fillna(0).astype(bool)
    flags["flagged"] = flags[["flag_decoy", "flag_tme", "flag_protocol_loss"]].fillna(False).any(axis=1)
    counted = flags["flagged"] & ~flags["curated_artefact"]
    n = int(counted.sum())
    available = [a for a, t in (("decoy", "decoy_recruitment"), ("tme", "tme_covariation"),
                                ("protocol", "protocol_shift")) if t in tables]
    verdict = "GO" if n >= rule["min_mps"] else ("PIVOT" if available else "INCOMPLETE")
    decision = {"verdict": verdict, "n_flagged": n, "min_mps": rule["min_mps"], "axes_available": available,
                "flagged_mps": flags.loc[counted, "mp"].tolist()}
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


def write_summary(path: Path, decision: dict, flags: pd.DataFrame, tables: dict, cfg: dict) -> None:
    r = cfg["decision"]
    parts = [
        "# P1 audit summary\n",
        f"**Verdict: {decision['verdict']}**: {decision['n_flagged']} retained MP(s) flagged "
        f"(GO needs >= {decision['min_mps']}). Axes available: {', '.join(decision['axes_available']) or 'none'}.\n",
        f"Rule: decoy recruitment >= {r['recruitment']}, or a positive TME association with partial R^2 >= "
        f"{r['partial_r2']} above the random-gene-set null (null_q95), or protocol loss "
        f"beyond the same-protocol null in >= {r['protocol_frac_beyond_null']} of paired tumours. "
        "MPs labelled as curated artefacts are excluded from the count. Instability and positional "
        "flags are informational.\n",
        "## Per-MP flags\n",
        _md_table(flags[[c for c in ["mp", "recruitment", "partial_r2", "null_q95", "slope", "kappa_ratio_median",
                                     "mean_delta_intensity", "frac_loss_beyond_null_intensity", "recovery_all",
                                     "best_arm", "q_excluding_known_clusters", "curated_artefact", "flagged"]
                         if c in flags]]),
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
