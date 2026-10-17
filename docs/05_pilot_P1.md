# Pilot P1: audit of the standard MP pipeline against real nuisance

P1 is the key go/no-go gate in the proposal ([02_TESSERA_proposal.md](02_TESSERA_proposal.md), section 8).
It asks whether metaprograms (MPs) that the standard pipeline produces, *and curators keep*,
can be explained by technical nuisance. The code is the `mpaudit` package in this repository.

> **Revision note (2026-10-17).** Revised in response to the 2026-10-16 review: one frozen GO rule that can
> fail, replacing a rule that would have passed on almost any cohort; decoys typed per arm (technical covariates
> versus CNA clone) with a graded `kappa_retained`; the TME axis adjusted for malignant→TME ambient spillover
> and compared with a specificity-matched null; the protocol axis based on co-expression (κ); gene-length,
> intronic-fraction and %ribo features; an NMF tolerance fix; an R parity harness; and a toy simulator with
> estimated (not true) clone, doublet and contamination covariates (the intronic fraction is a noisy readout of
> the planted nuclear fraction, as a measured covariate would be) and a null cohort. A code review of the same
> date tightened the implementation: criterion (b) counts only adjusted statistics that survive ρ̂ × 1.5, the
> protocol nulls are recalibrated with a log 1.5 minimum effect and an activity gate, dosage decoy hits need no
> minimum number of decoys (as frozen), labels reach re-derived MPs by gene overlap, the decision thresholds, the
> registered settings and the inputs the exclusions need are checked before any NMF runs, NMF fit diagnostics are
> reported for every run, and the toy's technical axes were recalibrated so that the standard pipeline can see them.
>
> **Rule versions.** None was registered or run on real data; `P1-GO v1.2` is the version to register.
> - `P1-GO v1.0` (2026-10-16): the first frozen rule.
> - `P1-GO v1.1` (2026-10-17): v1.0 plus three requirements: in (b) the association must also hold with every
>   ρ̂ × 1.5; in (c) a κ loss must also exceed 1.5-fold (log 1.5); and the replicate null of (c) is a t
>   prediction bound that needs same-protocol pairs from ≥ 3 distinct tumours.
> - `P1-GO v1.2` (2026-10-17): v1.1 plus two changes. (i) A κ activity gate in (c): a paired tumour is tested,
>   and a replicate pair enters the null, only where the MP is active in at least one of the two pieces, i.e. its
>   κ there exceeds max(0.05, the q95 of matched random gene sets' κ in the same samples), a simplified form of
>   the proposal's expression threshold κ*_t. Where an MP is not expressed its κ is background, so its "change"
>   is noise or the floor; a fixed floor alone does not adapt, because shared nuisance makes any gene set co-vary
>   a little. The condition is the same for paired tumours and replicate pairs so that their changes stay
>   exchangeable: a first draft that gated paired tumours on the reference piece alone, and replicate pairs on
>   both pieces, called about 15% of tested tumours losses in a null simulation (4–5% now). (ii) The registered
>   settings outside the `decision` block are enforced, each at exactly its registered value (seed included), and
>   so are the inputs without which an exclusion cannot run (curated labels for re-derived MPs, a dosage decoy
>   arm, a cytoband and gene positions): a config that deviates is refused unless `unregistered_variant: true`,
>   and the verdict is then labelled an unregistered variant in every output.

## The frozen GO rule (`P1-GO v1.2`)

Implemented in `src/mpaudit/report.py` (the single source of truth; `rule_text()` prints it into `summary.md`,
and the proposal states the same rule).

- **Retained MPs** are audited MPs not labelled as curated artefacts (`typer.labels`, a CSV of `mp,label` with 1
  for a curated artefact, e.g. the 11 MPs removed by hand in 3CA v2). Rediscovering curated artefacts is not news.
  Labels are matched by MP name. Re-derived MPs (`mps: rederive`) are named `MP1`, `MP2`, ... by the run, in no
  fixed order, so give `typer.labelled_gmt` (the gene sets of the labelled MPs): each re-derived MP inherits the
  label of its best-matching labelled set at Jaccard ≥ `min_jaccard` (0.3; `curated_label_mapping.csv`). Labels
  given by name for re-derived MPs are an unregistered variant. Labels that match no audited MP stop the run (or
  warn, when some match). Without `typer.labels` no MP is excluded and `summary.md` reads "Curated artefacts: not
  labelled": with `mps: rederive` that is an unregistered variant (rediscovered curated artefacts could count
  toward GO), with MPs from a GMT the config check warns.
  - *Gene sets of the 11 MPs removed by hand in 3CA v2.* They are absent from the published v2 MP list, so the
    `typer.labelled_gmt` for a re-derivation run has to be rebuilt from `3CA_v2/MP_analysis_new_batch.R` with the
    MPs exported before `inds_rm` is applied. **This has not been done yet**; until it is, a re-derivation run
    cannot exclude them and is an unregistered variant. An audit of the published v2 list (`mps:` that GMT) does
    not need them, since the list already leaves them out.
- **Pre-declared exclusions.**
  - *Expected TME-shared biology* (MPs named on `decision.expected_tme_shared`, or overlapping by ≥ 10 genes,
    `min_overlap`, a gene set of `expected_tme_shared_gmt`, a GMT file, or of `expected_tme_shared_sets`, inline in
    the config): such MPs cannot count through criterion (b), but still can through (a) or (c). IFN response or hypoxia are expected
    to co-vary with TME cells while being real. This departs on purpose from the 2026-10-16 review, which asked to
    drop listed MPs altogether: TME co-variation cannot tell niche-driven biology from leakage, but technical-decoy
    recruitment and a κ loss without dissociation do not rest on TME co-variation, and an expected TME-shared MP
    that the decoys reproduce, or that loses κ without dissociation (e.g. hypoxia induced by warm ischaemia), is the
    nuisance P1 looks for. Names work only when published MPs are audited (`mps:` a GMT); with `mps: rederive` they
    stop the run, so give gene sets. With the TME axis on and nothing pre-declared, the config check warns.
  - *Dosage*: an MP whose genes are arm-enriched (q excluding known gene families < 0.05) or that a dosage decoy
    arm (CNA clone) recruits (≥ 2 hits and `kappa_retained` ≥ 0.5, whatever the number of decoys) is labelled
    **dosage** and never counted, whatever its other support: it is genetic, not technical nuisance.
- **Nuisance support**, any of:
  - **(a) technical decoys**: a technical decoy arm recruits the MP (a decoy MP at Jaccard ≥ 0.3 to it) in ≥ 2
    of its 5 decoys (`decoy.n`, registered at exactly 5; ≥ 5 are needed), and the decoys keep ≥ 0.5 of its
    co-expression (`kappa_retained`);
  - **(b) TME co-variation**: a positive, spillover-adjusted association whose partial R² exceeds the q99 of
    random gene sets matched on expression and malignant specificity, in ≥ 10 samples with TME cells (no fixed R²
    bar). The statistic must be adjusted (`spill_adjusted`; with fewer than 5 malignant-specific markers the
    adjustment is skipped with a warning and the MP cannot count here), and the association must hold again when
    every ambient fraction ρ̂ is multiplied by 1.5 (`tme.rho_sensitivity`): `slope_rho` > 0 and `partial_r2_rho`
    above that run's q99 (`null_q99_rho`). Switching the re-run off (`rho_sensitivity: null`) or changing the
    factor is a deviation from the registered rule;
  - **(c) protocol**: a one-sided co-expression (κ) loss under the dissociation-free protocol in ≥ 50% of ≥ 4
    paired tumours in which the MP is active. A tumour counts as a loss when Δlog κ is below the null bound and
    below −log 1.5 (`protocol.kappa_min_effect`). The null bound is a one-sided 95% t prediction bound from
    same-protocol replicate pairs of ≥ 3 distinct tumours (`protocol.min_pairs`), otherwise matched random gene
    sets carried to the MP's κ level on the linear scale (table below). *Active*: the MP's κ under the reference
    or under the test protocol (mean over the tumour's samples under it) is above max(0.05, the q95 of the matched
    random gene sets' κ in the same samples) (`protocol.kappa_active`); a replicate pair enters the null only where
    the MP is active in at least one of its pieces.
- **Verdict.** **GO** if ≥ 3 retained, non-dosage MPs have nuisance support. **INCOMPLETE** if no criterion is
  adequately powered: (a) ≥ 5 decoys in some technical arm; (b) for some MP, ≥ 10 samples with TME cells, a finite
  null and the spillover adjusted; (c) for some MP, ≥ 4 paired tumours in which it is active and its Δlog κ and
  null bound are defined (`n_paired_tumours`: ≥ 20 malignant cells per sample under each protocol, panel retention
  ι ≥ 0.5, a usable null and the activity gate; a count of paired tumours from metadata is only an upper bound).
  Otherwise **PIVOT** (re-centre the project on calibrated membership, gene tiers and projection; see the
  proposal).
- **Informational only**: instability (`info_unstable` in `mp_flags.csv` marks `recovery_all` below
  `decision.unstable_recovery`, 0.5), positional flags outside the dosage rule, the within-sample κ ratio
  (`kappa_ratio_median`, not spillover-adjusted), A4 intensity and prevalence, and the trivial typer.

A GO is a reason to build the full method, not a verdict on biology. The TME axis alone cannot tell TME-shared
biology from a technical artefact; the full method (TESSERA) types provenance from several axes with
equivalence tests.

### Registration and the config check

Before running on real data, freeze and register the **whole config file** together with the `mpaudit` commit:
the `decision` block (thresholds and pre-declared lists), `decoy.arms`, the rule settings outside the `decision`
block, and the `typer.labels` / `typer.labelled_gmt` files (and an `expected_tme_shared_gmt`, if used). The
`decision` thresholds (`min_decoys`, `min_hits`, `min_kappa_retained`, `a1_null`, `min_a1_samples`,
`min_paired_tumours`, `protocol_frac`, `positional_q`, `min_overlap`, `min_mps`) must equal the values stated
above. The settings outside the block that the rule, or the statistics it reads, depend on must each equal their
registered value (`mpaudit.report.REGISTERED_SETTINGS`; a larger null size also re-draws the random sets, which
would give a second chance to pass):

| Setting | Registered value | Why it is part of the rule |
|---|---|---|
| `min_jaccard` | 0.3 | decoy hits (a); labels of re-derived MPs |
| `seed` | 0 | decoys, random-gene-set nulls and cell subsamples: another seed is another draw |
| `decoy.enabled`, `positional.enabled` | true | off removes a dosage exclusion, which makes GO easier |
| `decoy.n`, `decoy.method`, `decoy.n_depth_bins` | 5, `residual`, 1 | (a) is ≥ 2 hits among these decoys; depth strata leak real co-expression |
| `decoy.min_kappa`, `decoy.n_null_floor` | 0.02, 20 | which samples enter `kappa_retained` |
| `tme.spillover_adjust`, `tme.rho_sensitivity`, `tme.null_match` | true, 1.5, `expression+specificity` | (b) is defined on the adjusted statistic, its re-run and this null |
| `tme.n_null`, `tme.min_cells`, `tme.n_markers`, `tme.min_marker_fold` | 500, 20, 50, 8 | the null; which groups enter A1; the markers behind ρ̂ |
| `protocol.kappa_min_effect`, `protocol.kappa_active` | log 1.5 (0.4054651), 0.05 | the κ margin and the activity floor of (c) |
| `protocol.min_pairs`, `protocol.n_null_sets` | 3, 100 | the replicate null and the random-set null of (c); a larger `min_pairs` sends MPs to the narrower random-set null |

The exclusions must also be able to run: with `mps: rederive`, `typer.labels` and `typer.labelled_gmt` (labels by
name cannot be registered, since the run names the MPs); with decoys on, an arm named in `decision.dosage_arms`;
with the positional axis on, a cytoband and gene positions (`positional.gtf`, or `chromosome`/`start` columns in
`adata.var`; a GTF that places no gene is caught when it is read, before any NMF). A config without one of these is
a deviation like any other.

The arm names and contents, the pre-declared lists, the label and gene-set files, `unstable_recovery`, the `gavish`
block, `tme.covariates`, `tme.n_boot`, the keys, `recurrence_key`, the protocol names and the `enabled` switches of
the TME and protocol axes are filled in by the user and registered as filled in; the code does not compare them with
anything. Switching the TME or protocol axis off leaves its criterion unpowered, which can turn a PIVOT into an
INCOMPLETE. A run therefore carries the plain `P1-GO v1.2` label whenever the compared settings match, and whether it
used the registered file is checked by comparing its `run_config.json` (the normalised config) and input-file digests
with the frozen file.

`mpaudit p1` checks the whole config before any NMF runs and stops with the reason for the following (from the
command line, everything that does not need the cohort's `obs` columns is checked before the cohort is read):
- **deviations from the registered rule**, all listed with their registered values, unless the top-level
  `unregistered_variant` is true. Then it warns instead, and the verdict is labelled
  `P1-GO v1.2 (unregistered variant: setting=value, ...)` in `decision`, `rule_text`, `run_config.json` and
  `summary.md`, which opens with a bold line saying so;
- unknown `decision` keys (the pre-v1.0 keys `recruitment`, `partial_r2` and `protocol_frac_beyond_null` are named
  with their replacements);
- values of the wrong type or out of range: a single name for a list key is accepted; `decoy.enabled`,
  `tme.spillover_adjust`, `positional.enabled` and `unregistered_variant` must be true or false; `decoy.n`,
  `decoy.n_depth_bins`, `decoy.n_null_floor`, `tme.n_null`, `tme.min_cells`, `tme.n_markers`, `protocol.min_pairs`,
  `protocol.n_null_sets` and `seed` must be whole numbers (3.0 and "3" are refused); `protocol.kappa_active`,
  `protocol.kappa_min_effect` and `decoy.min_kappa` must be finite numbers ≥ 0; `tme.rho_sensitivity` a number > 0
  or null; `tme.min_marker_fold` a number > 0; `min_jaccard` must be in (0, 1]; `decoy.method` and
  `tme.null_match` must be known choices;
- a missing `expected_tme_shared_gmt`, label or GMT file; the legacy `decoy.strata`; decoy arms that neither
  `decision.technical_arms` nor `decision.dosage_arms` names; arm columns missing from `obs` (a `pct_long` arm needs
  gene lengths) or with missing values in malignant cells (NaN, None, or in a text column a placeholder such as
  "nan", "None" or ""; the decoys would condition on values nobody measured; the message says how to proceed); and
  names that cannot match re-derived MPs (above).

It warns when the TME axis can run but no expected-TME-shared MP or gene set is pre-declared, when MPs from a GMT
come without `typer.labels`, and, during the run, when NMF fits were stopped by the restart limit.

`summary.md` names the counted MPs, the MPs whose A1 statistic was not adjusted, those that pass at the estimated
ρ̂ but not at 1.5 × ρ̂, any decoy arm the rule ignores, whether curated labels were given, the capped and failed
NMF fits per kind of run and the sha256 of the input files. `mp_flags.csv` holds every per-MP column;
`run_config.json` holds the normalised config, the rule label and text, the deviations, the input-file digests
(`input_files`) and the fit counts (`nmf_fits`).

## What each axis measures

| Axis | Module | Statistic |
|---|---|---|
| **Decoys** | `mpaudit.decoy` | Per **arm** (a set of `obs` columns conditioned on jointly), the number of nuisance-only decoy cohorts (5, registered) in which the unmodified pipeline produces an MP with Jaccard ≥ 0.3 to the real MP, and `kappa_retained`: the median over samples of κ(decoy)/κ(real) on the MP's genes, counting only samples where the MP's real κ beats a matched-random null (q95 of 20 sets) and 0.02. Decoys keep depth, categorical strata (e.g. CNA clone) and continuous covariates (doublet score, contamination, %ribo, intronic fraction) cell by cell, and destroy all other co-expression. A clone arm types recruitment as dosage; a technical arm as nuisance. The example config conditions on the four technical covariates jointly in one arm, because each arm costs 5 more pipeline passes, so a technical hit is not attributed to a single covariate. To attribute it (e.g. to see whether `doublet_est`, which tracks cell size and cycling, drives a hit), declare one arm per covariate and list each in `decision.technical_arms`. Criterion (a) then counts a hit in any of them, which adds chances to hit, so fix the arm list in the registered config before running. |
| **TME co-variation, A1** | `mpaudit.tme` | Across samples, partial R² (beyond cancer type) of the malignant MP score on the same genes' pseudo-bulk score in co-resident TME cells, after **subtracting the spillover predicted** from malignant-derived ambient RNA (ρ̂ per sample and cell type from malignant-specific marker genes). Null: 500 random gene sets matched on expression and malignant specificity, put through the same adjustment (q95, q99, empirical p). With fewer than 5 markers (ambient RNA caps a marker's fold at about 1/ρ) the adjustment is skipped with a warning (`spill_adjusted` False). A second pass with every ρ̂ × 1.5 gives `partial_r2_rho`, `slope_rho` and its own null (`null_q99_rho`), because the random-set null does not absorb a leak left by a ρ̂ biased low. Also, informational only, the within-sample co-expression ratio TME/malignant (`kappa_ratio_median`, not spillover-adjusted). |
| **Protocol, A4** | `mpaudit.protocol` | Per paired tumour, Δlog κ of the MP's genes between protocols (genes detected in both; panel retention ι ≥ 0.5), tested only where the MP is active under at least one protocol (κ above max(0.05, the q95 of the matched random sets' κ in the same samples)). Null, chosen per statistic: same-protocol replicate pairs from ≥ 3 distinct tumours (for κ, a pair counts only where the MP is active in at least one piece), as a one-sided 95% t prediction bound (t(0.95, number of tumours) × the root mean square of the replicate changes; 5% false losses per tumour by Monte Carlo); otherwise random matched gene sets in the same paired tumours, whose 5th/95th percentiles of the *linear* κ change are carried to the MP's κ level (on the log scale, an additive shift of near-zero background κ looked like a loss). Neither null holds everything (the random-set null misses the MP's own sampling noise and regional heterogeneity), so a κ loss must also exceed log 1.5. Intensity and prevalence use protocol-specific control genes and are tumour-relative and informational (a cohort-wide shift is absorbed). |
| **Instability** | `mpaudit.instability` | Recovery (Jaccard ≥ 0.3) of each reference MP under one-factor perturbations: seed, rank grid, thresholds, gene universe, leave-one-study-out, 80% cell bootstraps. Method variants (no centring) are reported separately and not mixed into `recovery_all`. On small cohorts leave-one-study-out can drop MPs mechanically through the founder rule (> `min_group_size` partners); see the toy demo note. |
| **Positional, A3** | `mpaudit.positional` | Hypergeometric arm enrichment against the genes that could enter an MP (union of the per-sample NMF universes), BH within MP, also with known positional gene families removed. |
| **Trivial typer, P4** | `mpaudit.typer` | Gene content (ribosomal, mitochondrial, IEG/heat-shock, TME markers), median gene length and fraction > 100 kb, and the correlation of the MP score with log UMI, log genes, %mito, %ribo, doublet score and intronic fraction; AUROC against curated artefact labels. |

Why the gene-length and abundance features: in a user-reported breast cancer snRNA-seq cohort (HER2 project),
two pre-registered de novo MP searches returned mostly technical axes — several "luminal" programs made of
long genes (median gene length 317–437 kb) and a short-gene abundance axis tracking %ribo in 29 of 30 patients
(ρ 0.68–0.77) that reference matching labelled "EMT".

## Running it on real data

1. Build one AnnData with all cells (malignant **and** TME), raw integer counts in a layer, and obs columns for
   sample, malignant status, cell type, tumour/patient, study, cancer type and, where available, protocol, CNA
   clone, doublet score, an ambient-contamination estimate (DecontX/SoupX) and the intronic fraction.
   `mpaudit.data.load_3ca_study` loads one 3CA study folder (layout assumed; check it). 3CA study folders hold
   gene-level counts only, so the intronic fraction needs the raw reads re-processed (STARsolo with Velocyto
   features, velocyto or DropletQC), which waits for data access; the P0 census records which studies can
   supply it. Before freezing the config, remove it from `decoy.arms.technical` or set it to one constant for
   every cell of a sample without it (a column constant within a sample conditions nothing there): a listed arm
   column that is missing, or has NaN in malignant cells, stops the run. `pct_ribo`, `pct_mito`, `log_umi` and
   `log_genes` are added automatically (`qc.add_obs`), and so is `pct_long`, the percentage of a cell's counts on
   genes longer than 100 kb, when gene lengths exist (step 2). `pct_long` can stand in for the intronic fraction
   on the long-gene axis, but like `pct_ribo` it is computed from counts that include the audited MP's own genes,
   so it can inflate a long-gene MP's `kappa_retained`. It is in no default arm; decide before registration
   whether to add it.
2. For the positional axis and gene lengths (typer features, `pct_long`), download UCSC `cytoBand.txt.gz` (hg38)
   and a GENCODE GTF, or put `chromosome`, `start` and `gene_length` columns in `adata.var`. The GTF's
   `gene_name` values must match `var_names`: a GTF that places no gene falls back to `adata.var` with a warning,
   and one placing fewer than half of `var_names` warns. GTFs without `gene` records (UCSC `ncbiRefSeq.gtf`) are
   aggregated from their transcript and exon records, with a warning.
3. To audit the *published* MPs rather than re-derived ones, set `mps:` to a GMT of the Gavish 2023 or 3CA v2 MP
   lists. Instability then reports the recovery of the matching re-derived MP. With re-derived MPs, give the
   expected-TME-shared biology as gene sets and the curated labels with `typer.labelled_gmt` (see the rule above
   for the 11 MPs removed by hand in 3CA v2).
4. Copy `configs/p1_3ca_example.yaml`, fill in the `decision` block (expected TME-shared list, arms), the decoy
   arm columns and the curated labels, then freeze and register the whole file with the `mpaudit` commit (above)
   and run (a config that fails the checks that need no cohort stops before the cohort is read):

```bash
mpaudit p1 cohort.h5ad --config my_config.yaml --out results/p1 --n-jobs 8
```

**Runtime.** NMF dominates. With the default `tol = 1e-7` a fit runs several hundred iterations (about 2.5–3×
the time at `tol = 1e-5`); on a 7,000-gene × 2,000-cell sample one iteration takes about 50 ms on one core, so
ranks 4–9 × 10 restarts take on the order of 15–45 minutes per sample. Fits are cached (`cache_dir`; the key
includes a hash of the fitting code, and for decoy fits of the decoy generator, the decoy's arguments and the
values of the columns it conditions on), so threshold and leave-one-study-out perturbations are free; seed, rank,
universe, method-variant, bootstrap and decoy runs (5 per arm) each cost one more pass. The protocol axis draws `n_null_sets` random gene sets per MP in every paired tumour
and every tumour with same-protocol replicate pieces (they also give the activity threshold). Start with 100–200
tumours and use `--n-jobs`.

## Demo on the toy cohort

Three runs of P1 with 3CA defaults (`nrun` 10; no minimum of detected genes, as the toy has 1,500 genes) on
`mpaudit.sim` (seed 0; 10 tumours in 3 studies, 16 samples, 4 tumours profiled under both protocols):

- `python scripts/demo_sim.py --out results/demo`: the planted cohort. Biology (cycle, hypoxia, EMT), TME-shared
  biology (IFN, pre-declared as expected TME-shared by its gene set), a recurrent CNA, private programs, and
  nuisance (dissociation stress, ambient RNA and doublets, an abundance/%ribo axis, a long-gene axis), together
  with confounders that can fake nuisance (malignant→TME spillover, cycling cells with higher doublet scores,
  protocol capture bias).
- `python scripts/demo_sim.py --clean --out results/demo_clean`: the null cohort, with no planted nuisance
  programs, technical axes or CNA, but spillover, doublet–cycle coupling and capture bias kept. The frozen rule
  should not say GO there.
- Optional, `--true-doublet`: the technical arm uses the planted `doublet_score` (cycling cells score 0.3 higher)
  instead of `doublet_est`, to push a real program tightly coupled to a covariate through NMF recruitment and the
  `kappa_retained` ≥ 0.5 guard. With estimated covariates the doublet–cycle coupling is weak: `doublet_est` keeps
  about 0.01 of the cycle program's κ (`tests/test_decoy_arms.py`; 01 §7, G3, item (d)), so a `biology_counted` of 0
  in the default run does not test the `kappa_retained` ≥ 0.5 guard under strong coupling; `--true-doublet` does.
  Arm contents are registered as filled in but not compared by the code, so all three runs carry the same rule
  label; the `--true-doublet` arm is not one a registration would fix (the planted score cannot be observed in real
  data).

The toy has no curated artefacts and no curated labels, and its MPs are re-derived, so every demo run is the
unregistered variant `P1-GO v1.2 (unregistered variant: typer.labels=None)` (the demo sets
`unregistered_variant: true`; no other setting deviates), and `summary.md` opens with a bold line saying so.

The decoys condition on `clone_est`, `doublet_est` and `contamination_est`, estimated from the counts (never the
planted values), on `pct_ribo`, computed from the counts, and on `intronic_fraction`, which the toy draws as the
planted nuclear fraction plus N(0, 0.03): a noisy readout, as a measured covariate would be, so long-gene
recruitment by the technical arm is true by construction.

**Technical axes of the toy.** At their first magnitudes the standard pipeline could not see the abundance/%ribo
and long-gene axes: no per-sample NMF program had more than 1 abundance or 6 long genes in its top 50. They were
recalibrated on 2026-10-17 against a target fixed beforehand and judged on standard-pipeline visibility only,
before any audit, decision rule or demo was run on the new magnitudes: with the reference pipeline at the demo
settings, each axis must put ≥ 10 of its genes into one robust program in ≥ 8 of the 16 samples. The abundance
genes now scale with the square of the planted cytoplasmic factor and the long genes with 1 + 20 × the nuclear
fraction (`mpaudit.sim.CYTO_GAMMA`, `LONG_K`; every random draw is unchanged, and the null cohort is identical).
The axes now appear in 13 and 11 of the 16 samples and each forms its own MP. The cost is competition for ranks
4–9: EMT is robust in 7 samples instead of 9, contamination in 8 instead of 11 and the CNA in 5 instead of 8.
With a third of the library on one per-cell factor, every other gene co-varies a little after normalisation
(random 30-gene sets: κ about 0.23, against 0.09 before). Details and the candidates tried are in the
`mpaudit.sim` docstring.

Outputs besides the usual ones: `truth_comparison.csv` matches each audited MP to a planted set (Jaccard ≥ 0.2)
and lists its flags; `truth_summary.json` reports the verdict and rule label, the numbers of planted biology and
technical MPs and of those counted, biology MPs typed dosage (`biology_dosage`), biology MPs that the technical
decoys recruit in ≥ 2 of their 5 decoys whatever their `kappa_retained` (`biology_tech_recruited`, with a per-MP
list that says whether the `kappa_retained` guard passed and whether the MP is typed dosage: one that is not
counted was stopped by the guard, by a dosage label, or by both), how many counted MPs meet each criterion
(`counted_by_criterion`), the technical arm, the rule label with the deviations, and the code commit with a flag for
uncommitted changes.

RESULTS: to be filled from results/demo*/truth_summary.json after the rerun on the final code.

Leave-one-study-out on the toy is a small-cohort stress test, not a stability estimate. The toy has 10 tumours in
3 studies (16 samples), so leaving out a study drops 3–4 tumours, and EMT, IFN and the CNA are planted in only 5, 6
and 7 tumours. A founder needs more than `min_group_size` (5) partner programs at ≥ 10 shared genes, so MPs shared
by a subset of tumours can fall below the founder rule and leave-one-study-out recall is low by construction. It is
not evidence about the stability of 3CA MPs (more than 1,000 tumours in dozens of studies).

## Caveats

- **Checked against the original code.** The port was compared line by line with tiroshlab/3ca
  (`ITH_hallmarks/Generating_MPs/*.R`, `3CA_v2/MP_analysis_new_batch.R`) and R's NMF package
  (`R/algorithms-snmf.R`). Matched: per-sample gene universe and R's dropping of constant genes and all-zero
  cells, centring, the sNMF/R objective (eta² ridge on W), restarts on zero rows of H, best-of-10 selection by R's
  objective, the robust-program filter, the clustering rules, candidate order (sample, then rank), founder ties in
  dendrogram order, the v2 post-filters (≥ 2 studies; drop two-study MPs where one study gives a single NMF; ≥ 5
  NMFs) and, optionally, v2's final gene order (`v2_gene_order`). The transform is inferred, not matched (below).
- **NMF stopping and parity.** The optimiser (HALS) and stopping rule differ from R's exact alternating NNLS.
  At `tol = 1e-5` the port stopped about 5× earlier than R's rule would and top-50 lists were not stable
  (median 44/50 genes kept against a `tol = 1e-8` fit, on the toy before its technical axes were recalibrated);
  the default is therefore `tol = 1e-7` (49/50). Fits that hit `max_iter` (5,000; R's default is 20,000) are
  counted, as are runs stopped by the restart limit, and reported per run: `decoy_runs.csv`,
  `instability_runs.csv`, `run_config.json` (`nmf_fits`, with the reference run) and one line in `summary.md`; a
  fit stopped by the restart limit also raises a warning. **R-versus-port parity is untested**: the harness in
  `scripts/parity/` exports the exact NMF inputs, runs `NMF::nmf(method = "snmf/r")` and compares top-50 lists,
  robust programs and MPs, but its R side has not been executed (R could not be installed in the build
  environment). Run it on 3–5 real samples, with `gavish.tol: 1.0e-8` (decimal point, so YAML reads a number),
  before quoting absolute numbers. Decoys run through this port, not the R code.
- **Not reproducible from public code.** The function the MP script calls (`Return_All_cells_Separated`,
  `MP_analysis_new_batch.R:111`, with `log_normalize = TRUE`) is not public. It does both cell QC and
  log-normalisation. The transform log2(CPM/10 + 1) is inferred from the commented-out `MP_analysis_new_batch.R:99`
  and from `MPs_distribution/MP_distribution.R:47` (a scoring script); computing CPM over each cell's total across
  all genes is an assumption. The ≥ 1,000 detected-genes filter is taken from other 3CA scripts, and the minimum of
  50 malignant cells per sample is an assumption. The 11 MPs removed by hand in 3CA v2 cannot be reproduced; label
  them via `typer.labels` (their gene sets are not yet available, see the rule above). Exact distance ties in the
  dendrogram may resolve differently from R's `hclust`; sample order is alphabetical, by study, or as given
  (`sample_order`). R iterates its per-sample files in `list.files` order, suffix included (MGH100 before MGH10 in
  the C locale), so give the samples in R's file order with `sample_order: input` where tie order matters; v2 also
  runs the old batch before the new batch. Gene symbols are not upper-cased as in R. All-zero NMF components (dead
  columns) are skipped; R would keep them, with the sample's highest-mean genes as their "top genes".
- **Covariates computed from expression.** Doublet scores, contamination estimates, %ribo, `pct_long`, the
  intronic fraction and inferred clones come from the same counts (the intronic fraction from the same reads). A
  decoy conditioned on one keeps any real program correlated with it: cycling cells get higher doublet scores, and
  %ribo tracks translation, ribosome biogenesis and MYC activity, so a technical arm that includes %ribo keeps such
  programs, and an MP containing ribosomal-protein genes shares counts with %ribo and keeps them mechanically (it
  can exceed `kappa_retained` 1). `kappa_retained` grades how much is explained, but a program tightly tied to a
  covariate's estimator cannot be told from nuisance here: if the decoys recruit it and keep ≥ 0.5, it counts as
  nuisance support under (a). A covariate computed from an MP's own genes can push `kappa_retained` above 1.
  Covariates measured independently of expression (genotype demultiplexing, hashing, protocol twins) settle such
  cases. Preview on code before this revision (to be replaced by the rerun): in the `--clean` toy run the technical
  arm recruited the real cycle program in 5 of 5 decoys, and only its `kappa_retained` of 0.35 (< 0.5) kept it from
  counting.
- **Spillover adjustment.** It subtracts predicted ambient spillover using ρ̂ from malignant-specific markers;
  if real ambient RNA under-represents those markers, part of the leak remains, and the random-set null does not
  absorb it. Criterion (b) therefore also needs the association at 1.5 × ρ̂; that covers ρ̂ biased low by up to
  about a third (toy with the technical axes at their first magnitudes: at 0.5 × ρ̂ a malignant-only program still
  passes), and it can push weak genuine co-variation below its null where ρ̂ is already biased high. Without 5
  malignant-specific markers (ambient RNA above about 1/8 in a TME type, or a TME type sharing the malignant lineage
  genes, e.g. normal epithelium) no adjustment is made and A1 cannot support GO. `kappa_ratio` is not adjusted and
  is informational only.
- **Protocol nulls.** The random-set null holds changes shared by all genes as an additive change in κ, not
  shared proportional changes, the MP's own sampling noise or regional heterogeneity, so it is narrow and the
  log 1.5 minimum effect governs criterion (c) on that path. Toy (`simulate_cohort()` defaults, seed 0; planted
  sets; `protocol_shift_kappa(adata, keys, sets, "fresh", "fixed")` with its defaults, 100 random sets): without
  the margin, 6 losses and 7 gains in the 31 κ tests of the eight planted programs other than stress (32
  tumour-program pairs, one inactive); with it, one loss remains: CNA_7q in a tumour that does not carry the CNA
  but passed the activity gate (1 of its 3 tested tumours, short of the 50% of criterion (c)). The replicate null
  needs same-protocol pieces of ≥ 3 tumours, which public data rarely have; it assumes paired and replicate
  changes are exchangeable, which is why the activity gate selects both on the same condition. The gate removes
  tests at the floor; it does not tell which tumours carry a program. Its 5% level per piece holds only where an
  absent program's genes behave like the matched random sets: with the toy's technical axes at their first
  magnitudes no EMT or IFN test in a paired tumour lacking the program passed it (small test cohort), but on the
  recalibrated toy the abundance axis dilutes every gene by one per-cell factor, the planted sets co-vary beyond
  their matched random sets also where they are absent, and all 4 such tests on the small cohort passed, as did 2
  of 3 on the default cohort (one of them the CNA_7q loss above). Oliveira et al. 2025 is not a fresh-versus-fixed
  twin (FFPE Flex versus in situ on the same blocks).
- **Decoy leakage.** Strata leak a little real signal (per-stratum means carry each program's share). Depth is
  modelled continuously by default, continuous covariates enter the mean model, and binned strata are capped at
  ~50 cells.
- **Expression-only CNA.** A positional flag is a prompt for allele-aware follow-up (Numbat), not proof of
  dosage. Some real programs are positionally clustered gene families; they are listed and excluded in
  `q_excluding_known_clusters`.
- **Gene lists and features.** The IEG/heat-shock and TME-marker lists in `mpaudit.typer` are short illustrative
  cores; replace them with curated lists. Gene length is genomic span from the GTF; a short-gene abundance
  program has a *low* median length and is caught by %ribo correlation and ribosomal content instead. In the
  decoys, gene length enters only through a per-cell covariate: the intronic fraction or `pct_long` (median length
  and fraction > 100 kb are per-MP features and cannot enter a per-cell mean model). Where neither is in the
  technical arm (most processed 3CA matrices lack the intronic fraction, and `pct_long` is not in the default arm),
  long-gene MPs are screened only by the typer's length features, and a decoy "survival" of such an MP is not
  evidence against a nuisance origin.
- **Toy cohort.** `mpaudit.sim` results are true by construction. They show that the code runs and how each axis
  responds; they are not evidence about real tumours. The toy's technical axes were recalibrated on 2026-10-17
  (above), so seed-0 numbers measured before that date differ from the current toy.
