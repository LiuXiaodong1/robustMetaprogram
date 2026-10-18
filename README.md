# robustMetaprogram

Planned package for robust, calibrated discovery of tumour metaprograms (MPs) from
single-cell RNA-seq. It is meant to replace the threshold-and-curation pipeline
(per-tumour NMF → top-50 lists → overlap filters → greedy clustering → manual removal)
with MPs that have effect sizes, held-out replication, gene-level error control and a
typed provenance: intrinsic biology, TME-shared, genetic dosage, protocol-modulated or technical.

Status: **design stage, plus the pilot audit (P1)**.

- Design: [docs/](docs/README.md), starting with [docs/02_TESSERA_proposal.md](docs/02_TESSERA_proposal.md).
- Pilot audit code: the `mpaudit` package in this repo; usage, the frozen GO rule and caveats in
  [docs/05_pilot_P1.md](docs/05_pilot_P1.md).

## mpaudit: the P1 pilot audit

P1 is the month-0–3 go/no-go gate of the proposal. It asks whether MPs that the standard
pipeline produces, and curators keep, can be driven by technical nuisance, and decides with a
frozen rule (`P1-GO v1.3`, `mpaudit.report`) that can fail.

| Axis | Module | Question |
|---|---|---|
| Standard pipeline | `mpaudit.gavish` | Python port of the Gavish 2023 / 3CA MP pipeline, every constant exposed |
| Instability | `mpaudit.instability` | Do MPs survive other seeds, rank grids, thresholds, gene universes, leave-one-study-out, cell bootstraps? |
| Decoys | `mpaudit.decoy` | Can nuisance alone (per arm: a technical arm conditioning jointly on doublet score, ambient contamination, %ribo and intronic fraction; a dosage arm on CNA clone, whose recruitment types the MP as dosage) make the pipeline produce the MP, and how much of its co-expression does it explain? |
| TME co-variation (A1) | `mpaudit.tme` | Does the MP co-vary with the same genes in co-resident non-malignant cells, beyond malignant→TME spillover and a matched random-gene null? |
| Positional (A3) | `mpaudit.positional` | Are the MP's genes clustered on one chromosome arm? |
| Protocol (A4) | `mpaudit.protocol` | Is the MP's co-expression, compared on equal cell numbers and equal depth, lost under a dissociation-free protocol, beyond the null, the MP's own sampling noise and 1.5-fold, in at least half of the paired tumours? |
| Trivial typer (P4) | `mpaudit.typer` | How well do cheap features (gene content, gene length, QC and %ribo correlations) separate curated artefact MPs? |

```bash
pip install -e ".[test]"
pytest                                                          # unit and end-to-end tests on toy cohorts (about 7 min)
python scripts/demo_sim.py --out results/demo --n-jobs 4        # toy cohort with known truth
python scripts/demo_sim.py --clean --out results/demo_clean     # null cohort: the rule should not say GO
python scripts/demo_sim.py --true-doublet --out results/demo_true_doublet  # technical arm on the planted doublet score
# v1.3 toy results (seed 0): planted GO (4 technical MPs, 0 of 4 biology), --clean PIVOT, --true-doublet GO
# (cycle kept out only by kappa_retained 0.481); docs/05_pilot_P1.md, Demo
mpaudit p1 cohort.h5ad --config my_config.yaml --out results/p1 --n-jobs 8  # a filled-in copy of the template below
```

`configs/p1_3ca_example.yaml` is a template: as shipped it is refused, because its `ref/` cytoband and GTF files
are not in the repository and `typer.labels` / `typer.labelled_gmt` are empty. Copy it and fill these in first. The
labelled GMT of the 11 MPs removed by hand in 3CA v2 has not been rebuilt yet, so a re-derivation run
(`mps: rederive`) is an unregistered variant for now; an audit of a published MP list (`mps: <GMT>`) does not need
it.

Outputs: `summary.md` (verdict, rule and tables), per-axis CSVs, `reference_mps.gmt`, `run_config.json`.
`mpaudit p1` checks the whole config before any NMF runs: unknown keys, the `decision` block, decoy arms and
their `obs` columns, the keys against the `obs` columns, label, gene-set, cytoband and GTF files, the decision
thresholds and `mpaudit.report.REGISTERED_SETTINGS` against their registered values, and the inputs the exclusions
need (curated labels for re-derived MPs, a dosage decoy arm, a cytoband and gene positions that place genes on
chromosome arms). A deviation is refused unless the config sets `unregistered_variant: true`, and the verdict is
then labelled an unregistered variant. The rest of the config (arm contents, pre-declared lists, label files, NMF
settings other than `gavish.seed`, `counts_layer`, keys and every other key) is not compared by the code: freeze
and register the whole config, with the `mpaudit` commit and the label files, before running on real data, and
compare each run's `run_config.json` and input-file digests with it ([docs/05_pilot_P1.md](docs/05_pilot_P1.md)).
Axes that are switched on but cannot run are noted in `summary.md`. The toy demos have no curated labels, so they
run as an unregistered variant.

Calibration of criteria (c) and (a) on toy cohorts and synthetic κ laws (toy and synthetic only, not real data):
[scripts/calibration/README.md](scripts/calibration/README.md).

R-versus-port parity: [scripts/parity/README.md](scripts/parity/README.md) (export the exact NMF inputs,
run `NMF::nmf(method = "snmf/r")` in R, compare). Its R side has not yet been run.
