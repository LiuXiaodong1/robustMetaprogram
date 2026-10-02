# robustMetaprogram

Planned package for robust, calibrated discovery of tumour metaprograms (MPs) from
single-cell RNA-seq. It is meant to replace the threshold-and-curation pipeline
(per-tumour NMF → top-50 lists → overlap filters → greedy clustering → manual removal)
with MPs that have effect sizes, held-out replication, gene-level error control and a
typed provenance: intrinsic biology, TME-shared, genetic dosage, protocol-modulated or technical.

Status: **design stage, plus the pilot audit (P1)**.

- Design: [docs/](docs/README.md), starting with [docs/02_TESSERA_proposal.md](docs/02_TESSERA_proposal.md).
- Pilot audit code: the `mpaudit` package in this repo; usage in [docs/05_pilot_P1.md](docs/05_pilot_P1.md).

## mpaudit: the P1 pilot audit

P1 is the month-0–3 go/no-go gate of the proposal. It asks whether MPs that the standard
pipeline produces, and curators keep, can be driven by nuisance. It does this with:

| Axis | Module | Question |
|---|---|---|
| Standard pipeline | `mpaudit.gavish` | Python port of the Gavish 2023 / 3CA MP pipeline, every constant exposed |
| Instability | `mpaudit.instability` | Do MPs survive other seeds, rank grids, thresholds, gene universes, centring, leave-one-study-out, cell bootstraps? |
| Decoys | `mpaudit.decoy` | Can nuisance alone (depth, clone, doublet, ambient strata) make the pipeline produce the MP? |
| TME co-variation (A1) | `mpaudit.tme` | Does the MP co-vary with the same genes in co-resident non-malignant cells? |
| Positional (A3) | `mpaudit.positional` | Are the MP's genes clustered on one chromosome arm? |
| Protocol (A4) | `mpaudit.protocol` | Is the MP lost under a dissociation-free protocol beyond the same-protocol null? |
| Trivial typer (P4) | `mpaudit.typer` | How well do cheap features separate curated artefact MPs? |

```bash
pip install -e ".[test]"
pytest                                   # unit + end-to-end tests on a toy cohort
python scripts/demo_sim.py --out results/demo --n-jobs 4      # toy cohort with known truth
mpaudit p1 cohort.h5ad --config configs/p1_3ca_example.yaml --out results/p1 --n-jobs 8
```

Outputs: `summary.md` (verdict and tables), per-axis CSVs, `reference_mps.gmt`, `run_config.json`.
