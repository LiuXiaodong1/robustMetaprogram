# Pilot P1: audit of the standard MP pipeline against real nuisance

P1 is the key go/no-go gate in the proposal ([02_TESSERA_proposal.md](02_TESSERA_proposal.md), section 8).
It asks whether metaprograms (MPs) that the standard pipeline produces, *and curators keep*,
can be explained by technical nuisance. The code is the `mpaudit` package in this repository.

> **Revision note (2026-10-17).** Revised in response to the 2026-10-16 review: one frozen GO rule
> (`P1-GO v1.1`) that can fail, replacing a rule that would have passed on almost any cohort; decoys
> typed per covariate arm with a graded `kappa_retained`; the TME axis adjusted for malignant→TME
> ambient spillover and compared with a specificity-matched null; the protocol axis based on
> co-expression (κ); gene-length, intronic-fraction and %ribo features; an NMF tolerance fix; an R
> parity harness; and a toy simulator with estimated (not true) covariates and a null cohort. A code
> review of the same date tightened the implementation: criterion (b) counts only adjusted statistics
> that survive ρ̂ × 1.5, the protocol nulls are recalibrated with a log 1.5 minimum effect, dosage decoy
> hits need no minimum number of decoys (as frozen), labels reach re-derived MPs by gene overlap, and
> configs are checked before any NMF runs.
>
> **Rule version.** The review's tightening changed the rule, so its label moved from `P1-GO v1.0`
> (2026-10-16) to `P1-GO v1.1` (2026-10-17). v1.1 adds three requirements: in (b) the association must hold
> with every ρ̂ × 1.5; in (c) a κ loss must also exceed 1.5-fold (log 1.5); and the replicate null of (c) is a
> t prediction bound that needs same-protocol pairs from ≥ 3 distinct tumours. No real cohort has been run
> under either version; v1.1 is the version to register.

## The frozen GO rule (`P1-GO v1.1`)

Implemented in `src/mpaudit/report.py` (the single source of truth; the proposal states the same rule).
Freeze the `decision` block of the config — especially the pre-declared lists — before running on real data.
`mpaudit p1` checks the whole config before any NMF runs and stops with the reason: unknown `decision` keys
(the pre-v1.0 keys `recruitment`, `partial_r2` and `protocol_frac_beyond_null` are named with their replacements),
values of the wrong type (a single name for a list key is accepted), a missing `expected_tme_shared_gmt`, the
legacy `decoy.strata`, decoy arms that neither `decision.technical_arms` nor `decision.dosage_arms` names (it warns
when no arm is a dosage arm), and names that cannot match re-derived MPs (below).

- **Retained MPs** are audited MPs not labelled as curated artefacts (`typer.labels`, a CSV of `mp,label` with 1
  for a curated artefact, e.g. the 11 MPs removed by hand in 3CA v2). Rediscovering curated artefacts is not news.
  Labels are matched by MP name. Re-derived MPs (`mps: rederive`) are named `MP1`, `MP2`, ... by the run, in no
  fixed order, so give `typer.labelled_gmt` (the gene sets of the labelled MPs): each re-derived MP inherits the
  label of its best-matching labelled set at Jaccard ≥ 0.3 (`curated_label_mapping.csv`). Labels that match no
  audited MP stop the run (or warn, when some match).
- **Pre-declared exclusions.**
  - *Expected TME-shared biology* (`decision.expected_tme_shared`, by MP name, or `expected_tme_shared_gmt` /
    `expected_tme_shared_sets`, by an overlap of ≥ 10 genes): such MPs cannot count through criterion (b), but
    still can through (a) or (c). IFN response or hypoxia are expected to co-vary with TME cells while being real.
    Names work only when published MPs are audited (`mps:` a GMT); with `mps: rederive` they stop the run, so give
    gene sets.
  - *Dosage*: an MP whose genes are arm-enriched (q excluding known gene families < 0.05) or that a dosage decoy
    arm (CNA clone) recruits (≥ 2 hits and `kappa_retained` ≥ 0.5, whatever the number of decoys) is labelled
    **dosage** and never counted: it is genetic, not technical nuisance.
- **Nuisance support**, any of:
  - **(a) technical decoys**: a technical decoy arm recruits the MP in ≥ 2 of ≥ 5 decoys, and the decoys keep
    ≥ 0.5 of its co-expression (`kappa_retained`);
  - **(b) TME co-variation**: a positive, spillover-adjusted association whose partial R² exceeds the q99 of
    random gene sets matched on expression and malignant specificity (≥ 10 samples; no fixed R² bar). The
    statistic must be adjusted (`spill_adjusted`; with fewer than 5 malignant-specific markers the adjustment is
    skipped with a warning and the MP cannot count here), and the association must hold again when every ambient
    fraction ρ̂ is multiplied by 1.5 (`tme.rho_sensitivity`): `slope_rho` > 0 and `partial_r2_rho` above that
    run's q99 (`null_q99_rho`). With the re-run switched off, the estimated ρ̂ alone decides;
  - **(c) protocol**: a one-sided co-expression (κ) loss under the dissociation-free protocol in ≥ 50% of ≥ 4
    paired tumours. A tumour counts as a loss when Δlog κ is below the null bound and below −log 1.5
    (`protocol.kappa_min_effect`); the nulls are described in the table below.
- **Verdict.** **GO** if ≥ 3 retained, non-dosage MPs have nuisance support. **INCOMPLETE** if no criterion is
  adequately powered (≥ 5 decoys per technical arm; ≥ 10 samples with TME cells, a finite null and the spillover
  adjusted; ≥ 4 paired tumours). Otherwise **PIVOT** (re-centre the project on calibrated membership, gene tiers
  and projection; see the proposal).
- Informational only: instability, positional flags outside the dosage rule, and the trivial typer.

`summary.md` names the counted MPs, the MPs whose A1 statistic was not adjusted, those that pass at the estimated
ρ̂ but not at 1.5 × ρ̂, and any decoy arm the rule ignores; `mp_flags.csv` holds every per-MP column.

A GO is a reason to build the full method, not a verdict on biology. The TME axis alone cannot tell TME-shared
biology from a technical artefact; the full method (TESSERA) types provenance from several axes with
equivalence tests.

## What each axis measures

| Axis | Module | Statistic |
|---|---|---|
| **Decoys** | `mpaudit.decoy` | Per covariate **arm**, the number of nuisance-only decoy cohorts (≥ 5) in which the unmodified pipeline produces an MP with Jaccard ≥ 0.3 to the real MP, and `kappa_retained`: the median over samples of κ(decoy)/κ(real) on the MP's genes, counting only samples where the MP's real κ beats a matched-random null (q95 of 20 sets) and 0.02. Decoys keep depth, categorical strata (e.g. CNA clone) and continuous covariates (doublet score, contamination, %ribo, intronic fraction) cell by cell, and destroy all other co-expression. A clone arm types recruitment as dosage; a technical arm as nuisance. |
| **TME co-variation, A1** | `mpaudit.tme` | Across samples, partial R² (beyond cancer type) of the malignant MP score on the same genes' pseudo-bulk score in co-resident TME cells, after **subtracting the spillover predicted** from malignant-derived ambient RNA (ρ̂ per sample and cell type from malignant-specific marker genes). Null: 500 random gene sets matched on expression and malignant specificity, put through the same adjustment (q95, q99, empirical p). With fewer than 5 markers (ambient RNA caps a marker's fold at about 1/ρ) the adjustment is skipped with a warning (`spill_adjusted` False). A second pass with every ρ̂ × 1.5 gives `partial_r2_rho`, `slope_rho` and its own null (`null_q99_rho`), because the random-set null does not absorb a leak left by a ρ̂ biased low. Also the within-sample co-expression ratio TME/malignant (`kappa_ratio_median`, not spillover-adjusted). |
| **Protocol, A4** | `mpaudit.protocol` | Per paired tumour, Δlog κ of the MP's genes between protocols (genes detected in both; panel retention ι ≥ 0.5). Null, chosen per statistic: same-protocol replicate pairs from ≥ 3 distinct tumours, as a one-sided 95% t prediction bound (t(0.95, number of tumours) × the root mean square of the replicate changes; 5% false losses per tumour by Monte Carlo); otherwise random matched gene sets in the same paired tumours, whose 5th/95th percentiles of the *linear* κ change are carried to the MP's κ level (on the log scale, an additive shift of near-zero background κ looked like a loss). Neither null holds everything (the random-set null misses the MP's own sampling noise and regional heterogeneity), so a κ loss must also exceed log 1.5. Intensity and prevalence use protocol-specific control genes and are tumour-relative and informational (a cohort-wide shift is absorbed). |
| **Instability** | `mpaudit.instability` | Recovery (Jaccard ≥ 0.3) of each reference MP under one-factor perturbations: seed, rank grid, thresholds, gene universe, leave-one-study-out, 80% cell bootstraps. Method variants (no centring) are reported separately and not mixed into `recovery_all`. |
| **Positional, A3** | `mpaudit.positional` | Hypergeometric arm enrichment against the genes that could enter an MP (union of the per-sample NMF universes), BH within MP, also with known positional gene families removed. |
| **Trivial typer, P4** | `mpaudit.typer` | Gene content (ribosomal, mitochondrial, IEG/heat-shock, TME markers), median gene length and fraction > 100 kb, and the correlation of the MP score with log UMI, log genes, %mito, %ribo, doublet score and intronic fraction; AUROC against curated artefact labels. |

Why the gene-length and abundance features: in a user-reported breast cancer snRNA-seq cohort (HER2 project),
two pre-registered de novo MP searches returned mostly technical axes — several "luminal" programs made of
long genes (median gene length 317–437 kb) and a short-gene abundance axis tracking %ribo in 29 of 30 patients
(ρ 0.68–0.77) that reference matching labelled "EMT".

## Running it on real data

1. Build one AnnData with all cells (malignant **and** TME), raw integer counts in a layer, and obs columns for
   sample, malignant status, cell type, tumour/patient, study, cancer type and, where available, protocol, CNA
   clone, doublet score, an ambient-contamination estimate (DecontX/SoupX) and the intronic fraction
   (STARsolo/velocyto). `mpaudit.data.load_3ca_study` loads one 3CA study folder (layout assumed; check it).
   `pct_ribo`, `pct_mito`, `log_umi` and `log_genes` are added automatically (`qc.add_obs`).
2. For the positional axis and gene lengths, download UCSC `cytoBand.txt.gz` (hg38) and a GENCODE GTF, or put
   `chromosome`, `start` and `gene_length` columns in `adata.var`. The GTF's `gene_name` values must match
   `var_names`: a GTF that places no gene falls back to `adata.var` with a warning, and one placing fewer than half
   of `var_names` warns. GTFs without `gene` records (UCSC `ncbiRefSeq.gtf`) are aggregated from their transcript
   and exon records, with a warning.
3. To audit the *published* MPs rather than re-derived ones, set `mps:` to a GMT of the Gavish 2023 or 3CA v2 MP
   lists. Instability then reports the recovery of the matching re-derived MP. With re-derived MPs, give the
   expected-TME-shared biology as gene sets and the curated labels with `typer.labelled_gmt`.
4. Copy `configs/p1_3ca_example.yaml`, fill in the `decision` block (expected TME-shared list, arms), freeze it,
   then run (a bad config stops before the cohort is read):

```bash
mpaudit p1 cohort.h5ad --config my_config.yaml --out results/p1 --n-jobs 8
```

**Runtime.** NMF dominates. With the default `tol = 1e-7` a fit runs several hundred iterations (about 2.5–3×
the time at `tol = 1e-5`); on a 7,000-gene × 2,000-cell sample one iteration takes about 50 ms on one core, so
ranks 4–9 × 10 restarts take on the order of 15–45 minutes per sample. Fits are cached (`cache_dir`; the key
includes a hash of the fitting code), so threshold and leave-one-study-out perturbations are free; seed, rank,
universe, method-variant, bootstrap and decoy runs (5 per arm) each cost one more pass. Start with 100–200
tumours and use `--n-jobs`.

## Demo on the toy cohort

`python scripts/demo_sim.py --out results/demo` (planted nuisance) and `python scripts/demo_sim.py --clean`
(null cohort: no planted nuisance programs or CNA, but spillover, doublet–cycle coupling and capture bias kept)
run P1 with 3CA defaults on `mpaudit.sim`. Results of the revised axes are being regenerated and will be added
here; the earlier demo numbers were produced by the pre-review code and no longer apply.

## Caveats

- **Checked against the original code.** The port was compared line by line with tiroshlab/3ca
  (`ITH_hallmarks/Generating_MPs/*.R`, `3CA_v2/MP_analysis_new_batch.R`) and R's NMF package
  (`R/algorithms-snmf.R`). Matched: transform, per-sample gene universe and R's dropping of constant genes and
  all-zero cells, centring, the sNMF/R objective (eta² ridge on W), restarts on zero rows of H, best-of-10
  selection by R's objective, the robust-program filter, the clustering rules, candidate order (sample, then
  rank), founder ties in dendrogram order, the v2 post-filters (≥ 2 studies; drop two-study MPs where one study
  gives a single NMF; ≥ 5 NMFs) and, optionally, v2's final gene order (`v2_gene_order`).
- **NMF stopping and parity.** The optimiser (HALS) and stopping rule differ from R's exact alternating NNLS.
  At `tol = 1e-5` the port stopped about 5× earlier than R's rule would and top-50 lists were not stable
  (median 44/50 genes kept against a `tol = 1e-8` fit on the toy); the default is therefore `tol = 1e-7` (49/50).
  Fits that hit `max_iter` (5,000; R's default is 20,000) are counted, as are runs stopped by the restart limit. **R-versus-port parity is untested**: the
  harness in `scripts/parity/` exports the exact NMF inputs, runs `NMF::nmf(method = "snmf/r")` and compares
  top-50 lists, robust programs and MPs, but its R side has not been executed (R could not be installed in the
  build environment). Run it on 3–5 real samples, with `gavish.tol: 1.0e-8` (decimal point, so YAML reads a
  number), before quoting absolute numbers.
  Decoys run through this port, not the R code.
- **Not reproducible from public code.** The cell-QC step the MP script calls (`Return_All_cells_Separated`) is
  not public; the ≥ 1,000 detected-genes filter is taken from other 3CA scripts, and the minimum of 50 malignant
  cells per sample is an assumption. The 11 MPs removed by hand in 3CA v2 cannot be reproduced; label them via
  `typer.labels`. Exact distance ties in the dendrogram may resolve differently from R's `hclust`; sample order is
  alphabetical, by study, or as given (`sample_order`). R iterates its per-sample files in `list.files` order,
  suffix included (MGH100 before MGH10 in the C locale), so give the samples in R's file order with
  `sample_order: input` where tie order matters; v2 also runs the old batch before the new batch. Gene
  symbols are not upper-cased as in R. All-zero NMF components (dead columns) are skipped; R would keep them,
  with the sample's highest-mean genes as their "top genes".
- **Covariates computed from expression.** Doublet scores, contamination estimates and inferred clones come
  from the same counts. A decoy conditioned on one keeps any real program correlated with it (cycling cells get
  higher doublet scores); `kappa_retained` grades how much is explained, but a program tightly tied to a
  covariate's estimator cannot be told from nuisance here. A covariate computed from an MP's own genes can push
  `kappa_retained` above 1. Covariates measured independently of expression (genotype demultiplexing, hashing,
  protocol twins) settle such cases.
- **Spillover adjustment.** It subtracts predicted ambient spillover using ρ̂ from malignant-specific markers;
  if real ambient RNA under-represents those markers, part of the leak remains, and the random-set null does not
  absorb it. Criterion (b) therefore also needs the association at 1.5 × ρ̂; that covers ρ̂ biased low by up to
  about a third (toy: at 0.5 × ρ̂ a malignant-only program still passes), and it can push weak genuine
  co-variation below its null where ρ̂ is already biased high. Without 5 malignant-specific markers (ambient RNA
  above about 1/8 in a TME type, or a TME type sharing the malignant lineage genes, e.g. normal epithelium) no
  adjustment is made and A1 cannot support GO. `kappa_ratio` is not adjusted.
- **Protocol nulls.** The random-set null holds changes shared by all genes as an additive change in κ, not
  shared proportional changes, the MP's own sampling noise or regional heterogeneity, so it is narrow and the
  log 1.5 minimum effect governs criterion (c) on that path (toy: without the margin, 3 losses and 7 gains in 32
  tests of unaffected programs; with it, none). The replicate null needs same-protocol pieces of ≥ 3 tumours,
  which public data rarely have; it assumes paired and replicate changes are exchangeable. Oliveira et al.
  2025 is not a fresh-versus-fixed twin (FFPE Flex versus in situ on the same blocks).
- **Decoy leakage.** Strata leak a little real signal (per-stratum means carry each program's share). Depth is
  modelled continuously by default, continuous covariates enter the mean model, and binned strata are capped at
  ~50 cells.
- **Expression-only CNA.** A positional flag is a prompt for allele-aware follow-up (Numbat), not proof of
  dosage. Some real programs are positionally clustered gene families; they are listed and excluded in
  `q_excluding_known_clusters`.
- **Gene lists and features.** The IEG/heat-shock and TME-marker lists in `mpaudit.typer` are short illustrative
  cores; replace them with curated lists. Gene length is genomic span from the GTF; a short-gene abundance
  program has a *low* median length and is caught by %ribo correlation and ribosomal content instead.
- **Toy cohort.** `mpaudit.sim` results are true by construction. They show that the code runs and how each axis
  responds; they are not evidence about real tumours.
