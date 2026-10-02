# Pilot P1: audit of the standard MP pipeline against real nuisance

P1 is the key go/no-go gate in the proposal ([02_TESSERA_proposal.md](02_TESSERA_proposal.md), section 8).
It asks whether metaprograms (MPs) that the standard pipeline produces, *and curators keep*,
can be explained by nuisance. The code is the `mpaudit` package in this repository.

## What it measures

| Axis | Statistic | Flag used in the GO rule |
|---|---|---|
| **Decoys** (`mpaudit.decoy`) | Fraction of nuisance-only decoy cohorts in which the unmodified pipeline produces an MP with Jaccard ≥ 0.3 to the real MP | recruitment ≥ 0.3 |
| **TME co-variation, A1** (`mpaudit.tme`) | Partial R² of per-sample malignant MP score on the same genes' score in co-resident TME cells, beyond cancer type (bootstrap 90% CI); within-sample co-expression ratio TME/malignant (`kappa_ratio_median`) | partial R² ≥ 0.3 |
| **Protocol, A4** (`mpaudit.protocol`) | Paired change in MP intensity and prevalence between protocols of the same tumour, against same-protocol region pairs | loss beyond the null in ≥ 50% of paired tumours |
| **Instability** (`mpaudit.instability`) | Recovery (Jaccard ≥ 0.3) of each reference MP under one-factor perturbations: seed, rank grid, thresholds, gene universe, centring, leave-one-study-out, 80% cell bootstraps | informational |
| **Positional, A3** (`mpaudit.positional`) | Hypergeometric arm enrichment, BH within MP, also with known positional gene families removed | informational |
| **Trivial typer, P4** (`mpaudit.typer`) | Gene-content and QC-correlation features; AUROC against curated artefact labels | baseline any provenance method must beat |

**Decision rule.** GO if at least 3 retained MPs carry at least one GO flag. MPs listed as
curated artefacts (`typer.labels`, e.g. the 11 MPs removed by hand in 3CA v2) do not count,
because rediscovering them is not news. Otherwise PIVOT (re-centre the project on calibrated
membership, gene tiers and projection; see the proposal), or INCOMPLETE when no GO axis could run.

A GO flag is a reason for follow-up, not a verdict on biology. IFN response and hypoxia are
expected to be TME-shared biology, so they can carry the A1 flag while being real. The full
method (TESSERA) separates those cases with protocol contrasts and equivalence tests.

## Running it on real data

1. Build one AnnData with all cells (malignant **and** TME), raw integer counts in a layer, and
   obs columns for sample, malignant status, cell type, tumour/patient, study, cancer type and,
   where available, protocol, CNA clone and doublet score. `mpaudit.data.load_3ca_study` loads one
   3CA study folder (the file layout is assumed; check it against your download).
2. For the positional axis, download UCSC `cytoBand.txt.gz` (hg38) and a GENCODE GTF, or put
   `chromosome` and `start` columns in `adata.var`.
3. To audit the *published* MPs rather than re-derived ones, set `mps:` to a GMT of the Gavish 2023
   or 3CA v2 MP lists. Instability still uses the re-derived reference.
4. Copy and edit `configs/p1_3ca_example.yaml`, then run:

```bash
mpaudit p1 cohort.h5ad --config my_config.yaml --out results/p1 --n-jobs 8
```

**Runtime.** NMF dominates: about 2–4 s per fit for a 7,000-gene × 1,500-cell sample on one
core, so about 3 minutes per sample for ranks 4–9 × 10 restarts. Fits are cached (`cache_dir`), so
threshold and leave-one-study-out perturbations are free. Seed, rank, universe, centring,
bootstrap and decoy runs each cost one more pass over the cohort. Start with 100–200 tumours.

## Caveats

- **NMF parity.** `mpaudit.nmf` minimises the same sNMF/R objective as R's `NMF::nmf(method = "snmf/r")`
  but with a different optimiser (HALS). Seed-to-seed results therefore differ from the R code.
  Contrasts between perturbations inside `mpaudit` are like-for-like. Before quoting absolute numbers
  for the 3CA pipeline, check agreement with the original R code on a few samples.
- **Gavish port details.** Matching samples by exact name (the R code uses `grep`), tie-breaking order
  and the minimum malignant cells per sample (default 50) differ or are assumptions; see the
  `mpaudit.gavish` docstring.
- **Decoy leakage.** Strata leak a little real signal (per-stratum means carry each program's share).
  Depth is therefore modelled continuously by default, and binned strata are capped at ~50 cells.
- **Expression-only CNA.** A positional flag is a prompt for allele-aware follow-up (Numbat), not proof
  of dosage. Some real programs are positionally clustered gene families; they are listed and excluded
  in `q_excluding_known_clusters`.
- **Gene lists.** The IEG/heat-shock and TME-marker lists in `mpaudit.typer` are short illustrative
  cores. Replace them with curated lists for real analyses.
- **Toy cohort.** `mpaudit.sim` results are true by construction. They show that the code runs and that
  each axis responds as designed; they are not evidence about real tumours.
