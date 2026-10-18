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
> protocol nulls are recalibrated with a log 1.5 minimum effect and an activity gate (removed again in v1.3,
> below), dosage decoy hits need no minimum number of decoys (as frozen), labels reach re-derived MPs by gene
> overlap, the decision thresholds, the registered settings and the inputs the exclusions need are checked before
> any NMF runs, NMF fit diagnostics are reported for every run, and the toy's technical axes were recalibrated so
> that the standard pipeline can see them.
>
> **Revision note (2026-10-18).** A design comparison of fixes for criteria (c) and (a) found that the v1.2
> activity gate of (c) selected tumours on the noisy κ that it then tested, and that `kappa_retained` let in
> samples where the MP is absent. The rule is now `P1-GO v1.3` (below). A second code review led to these
> changes, none of which changes the rule: `gavish.seed` is registered; the config check refuses unknown keys
> and missing cytoband or GTF files; the run checks before any NMF that the gene positions and the cytoband
> place genes on chromosome arms; skipped axes are noted in `summary.md`; and the NMF cache key now hashes each
> sample's count matrix itself. A third review of the same day, still before registration, changed v1.3's (c):
> κ is now compared at equal depth (a depth difference between the protocols made false losses and a false (c)
> on the toy), the panel is taken on the counts κ is computed on, and no draw depends on the other MPs audited
> or on the cohort's cell order. It also made the cell cap and the panel thresholds fixed constants of the rule,
> refused `protocol.n_null_sets` below 20, reported the samples entering `kappa_retained` and the gain mirror of
> (c), and corrected seed labels and several figures in these documents.
>
> **Rule versions.** None was registered or run on real data; `P1-GO v1.3` is the version to register.
> - `P1-GO v1.0` (2026-10-16): the first frozen rule.
> - `P1-GO v1.1` (2026-10-17): v1.0 plus three requirements: in (b) the association must also hold with every
>   ρ̂ × 1.5; in (c) a κ loss must also exceed 1.5-fold (log 1.5); and the replicate null of (c) is a t
>   prediction bound that needs same-protocol pairs from ≥ 3 distinct tumours.
> - `P1-GO v1.2` (2026-10-17): v1.1 plus two changes. (i) A κ activity gate in (c): a paired tumour was tested,
>   and a replicate pair entered the null, only where the MP's κ in at least one of the two pieces exceeded
>   max(0.05, the q95 of matched random gene sets' κ in the same samples). (ii) The registered settings outside
>   the `decision` block are enforced, each at exactly its registered value (seed included), and so are the
>   inputs without which an exclusion cannot run (curated labels for re-derived MPs, a dosage decoy arm, a
>   cytoband and gene positions): a config that deviates is refused unless `unregistered_variant: true`, and the
>   verdict is then labelled an unregistered variant in every output.
> - `P1-GO v1.3` (2026-10-18): v1.2 with criteria (c) and (a) changed. Why: the gate selected tumours on the κ
>   that it then tested. A tumour whose κ was high by chance passed the gate and then tended to drop. In the
>   calibration harness (`scripts/calibration/`, v1.2 reference), when the MP's κ in each piece was drawn like
>   its random sets' κ, up to 0.31 of tested tumours were called losses on synthetic laws (0.05–0.31 across the
>   exchangeable laws on the random-set path). On the toy null,
>   0.13–0.17 of tested (gene set, tumour) pairs without the program were called losses, and an unplanted random
>   gene set met (c). Also, v1.2's `kappa_retained` entry (κ above the q95 of 20 random sets) let in samples
>   where the MP is absent, and a biology MP kept in 2 of 10 tumours reached 0.66 > 0.5 under the technical
>   arm. The changes:
>   - (c) has no activity gate: every paired tumour with a defined change, null bound and self bound is tested,
>     and the activity columns are informational;
>   - the replicate bound is the lower of the linear-scale t bound, carried to the MP's κ, and the log-scale t
>     bound;
>   - the random-set bound is the floor(0.05 (n + 1))-th order statistic of the sets' changes (the 5th of 100);
>   - a loss must also beat 1.645 × the bootstrap SE of the MP's own change (the self bound);
>   - every sample of a tumour is scored on the same number of cells and at the same depth (its counts binomially
>     thinned to the tumour's lowest median count per cell), and the panel is taken on those counts;
>   - an MP's random sets are drawn from a stream of its own per tumour, and every per-sample draw from the
>     sample's own stream over its cells in barcode order, so no result depends on the other MPs audited or on the
>     cohort's cell order;
>   - the cell cap (2000) and the panel thresholds are fixed constants of the rule (below);
>   - `kappa_retained` is the excess ratio over 50 matched random gene sets, with a 3-SD entry
>     (`decoy.n_null_floor` 50).
>
>   Every random draw differs from v1.2, so no v1.2 output carries over.

## The frozen GO rule (`P1-GO v1.3`)

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
    of its 5 decoys (`decoy.n`, registered at exactly 5; ≥ 5 are needed), and `kappa_retained` ≥ 0.5.
    `kappa_retained` is the share of the MP's co-expression in excess of matched random gene sets that the
    decoys keep. Per malignant sample with ≥ 20 cells, G is the set of the MP's genes detected in ≥ 5% of the
    sample's cells (at least 5 genes). κ_R is κ(G) on the real counts and κ_D the same on the decoy's counts of the same cells.
    bg_R and bg_D are the median κ of 50 random gene sets (`decoy.n_null_floor`) matched gene by gene to G on
    mean expression and detection and scored on the same counts. The sample's ratio is
    (κ_D − bg_D)/(κ_R − bg_R). A sample enters only if ≥ 20 of its sets have a finite κ and
    κ_R − bg_R > max(0.02, 3 SD of the sets' real κ) (`decoy.min_kappa` 0.02 is the minimum excess).
    `kappa_retained` is the median ratio over the entering samples, averaged over the arm's decoys. If no sample
    enters it is undefined, and the MP can then neither meet (a) in that arm nor be typed dosage by it;
  - **(b) TME co-variation**: a positive, spillover-adjusted association whose partial R² exceeds the q99 of
    random gene sets matched on expression and malignant specificity, in ≥ 10 samples with TME cells (no fixed R²
    bar). The statistic must be adjusted (`spill_adjusted`; with fewer than 5 malignant-specific markers the
    adjustment is skipped with a warning and the MP cannot count here), and the association must hold again when
    every ambient fraction ρ̂ is multiplied by 1.5 (`tme.rho_sensitivity`): `slope_rho` > 0 and `partial_r2_rho`
    above that run's q99 (`null_q99_rho`). Switching the re-run off (`rho_sensitivity: null`) or changing the
    factor is a deviation from the registered rule;
  - **(c) protocol**: a one-sided co-expression (κ) loss under the dissociation-free protocol in ≥ 50% of ≥ 4
    *tested* paired tumours. Every paired tumour with a defined change, null bound and self bound is tested;
    there is no activity gate.
    - *Change.* κ is the mean latent correlation of the MP's genes in the tumour's panel (genes detected in
      ≥ 5% of its κ cells, after thinning, under each protocol; panel retention ι ≥ 0.5 and ≥ 5 panel genes),
      in samples with ≥ 20 malignant cells. It is averaged over the tumour's samples under each protocol and
      floored at 0.01. Δlog κ = log κ_test − log κ_ref.
    - *Equal cells and depth.* Every sample of a tumour is scored on the same number of cells: the tumour's
      smallest sample, at most 2000, drawn from each sample's own random stream. Their counts are binomially
      thinned to the tumour's lowest median total count per cell (`kappa_thin` in `protocol_levels.csv`),
      because κ's estimate falls as counts per cell fall.
    - *Loss.* A tested tumour is a loss when Δlog κ is below the null bound, below the self bound and below
      −log 1.5 (`protocol.kappa_min_effect`).
    - *Null bound with replicates.* This applies when same-protocol replicate pairs (κ finite in both pieces)
      come from ≥ 3 distinct tumours (`protocol.min_pairs`). The bound is the lower of two one-sided 95% t
      prediction bounds: one on the linear κ changes, carried to the MP's κ level, and one on the log κ changes.
    - *Null bound otherwise.* The floor(0.05 (n + 1))-th smallest linear κ change of the n ≥ 20 matched random
      gene sets of the tumour with a defined change (100 drawn per tumour and MP, `protocol.n_null_sets`; the
      5th of 100), carried to the MP's κ level. If fewer than 20 sets have a defined change, or the k-th
      smallest change is not below the k-th largest, the tumour is not tested.
    - *Self bound.* −1.645 × the bootstrap SE of the MP's own linear κ change, carried to its κ level. The SE
      of each sample's κ comes from 50 resamples of its κ cells, drawn per seed, sample and MP. A tumour with an
      undefined SE is not tested.
- **Verdict.** **GO** if ≥ 3 retained, non-dosage MPs have nuisance support. **INCOMPLETE** if no criterion is
  adequately powered: (a) ≥ 5 decoys in some technical arm; (b) for some MP, ≥ 10 samples with TME cells, a finite
  null and the spillover adjusted; (c) for some MP, ≥ 4 tested paired tumours, i.e. paired tumours with Δlog κ,
  the null bound and the self bound defined (`n_paired_tumours`: ≥ 20 malignant cells per sample under each
  protocol, panel retention ι ≥ 0.5 and a usable null; a count of paired tumours from metadata is only an upper
  bound). Otherwise **PIVOT** (re-centre the project on calibrated membership, gene tiers and projection; see
  the proposal).
- **Informational only**: instability (`info_unstable` in `mp_flags.csv` marks `recovery_all` below
  `decision.unstable_recovery`, 0.5), positional flags outside the dosage rule, the within-sample κ ratio
  (`kappa_ratio_median`, not spillover-adjusted), the κ activity columns of the protocol axis
  (`kappa_active*`, `protocol.kappa_active`), A4 intensity and prevalence, and the trivial typer.

`report.rule_text()` states the same rule in one paragraph; it is printed into every run's `summary.md` and
`run_config.json`.

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
registered value (`mpaudit.report.REGISTERED_SETTINGS`). A null size is not an "at least" setting either. Another
size moves the estimated null quantile or background: it adds or drops sets, since each gene set (`tme.n_null`),
each sample and gene set (`decoy.n_null_floor`) and each tumour and MP (`protocol.n_null_sets`) has its own random
stream. Choosing the size would therefore give a second chance to pass.

| Setting | Registered value | Why it is part of the rule |
|---|---|---|
| `min_jaccard` | 0.3 | decoy hits (a); labels of re-derived MPs |
| `seed` | 0 | decoys, random-gene-set nulls, bootstrap resamples, cell subsamples and depth thinning: another seed is another draw |
| `gavish.seed` | 0 | NMF restarts of the reference run and of every decoy fit: another seed is another draw |
| `decoy.enabled`, `positional.enabled` | true | off removes a dosage exclusion, which makes GO easier |
| `decoy.n`, `decoy.method`, `decoy.n_depth_bins` | 5, `residual`, 1 | (a) is ≥ 2 hits among these decoys; depth strata leak real co-expression |
| `decoy.min_kappa`, `decoy.n_null_floor` | 0.02, 50 | `kappa_retained`: the minimum excess κ_R − bg_R of an entering sample; the random sets behind the background (median) and its SD |
| `tme.spillover_adjust`, `tme.rho_sensitivity`, `tme.null_match` | true, 1.5, `expression+specificity` | (b) is defined on the adjusted statistic, its re-run and this null |
| `tme.n_null`, `tme.min_cells`, `tme.n_markers`, `tme.min_marker_fold` | 500, 20, 50, 8 | the null; which groups enter A1; the markers behind ρ̂ |
| `protocol.kappa_min_effect` | log 1.5 (0.4054651) | the κ margin of (c) |
| `protocol.min_pairs`, `protocol.n_null_sets` | 3, 100 | the replicate null and the random-set null of (c); a larger `min_pairs` sends MPs to the random-set null |

`protocol.kappa_active` (0.05) is no longer registered: it only sets the informational `kappa_active*` columns.
Constants of the rule that are fixed in code, not set in the config (`mpaudit.report.FIXED_CONSTANTS`; stated in
`rule_text()`, pinned by `tests/test_report_rule.py`; changing one is a new rule version):

| Constant | Value | Role |
|---|---|---|
| `protocol.KAPPA_FLOOR` | 0.01 | κ is floored before the log |
| `protocol.NULL_ALPHA` | 0.05 | level of each one-sided null bound: t(0.95) for replicates, order statistic floor(0.05 (n + 1)) for random sets |
| `protocol.MIN_RANDOM_SETS` | 20 | random sets with a defined change needed for a tumour's random-set bound (then k = 1) |
| `protocol.KAPPA_BOOT`, `protocol.KAPPA_SELF_Z` | 50, 1.645 | bootstrap resamples behind the self bound, and its z |
| `protocol.MAX_KAPPA_CELLS` | 2000 | κ cells per sample at most (otherwise the tumour's smallest sample) |
| `protocol.MIN_CELLS`, `protocol.MIN_DETECT`, `protocol.MIN_GENES`, `protocol.MIN_IOTA` | 20, 0.05, 5, 0.5 | (c)'s samples and panel: malignant cells a sample needs, detection of a panel gene under each protocol, panel genes and panel retention needed for κ |
| `decoy.EXCESS_ENTRY_SD`, `decoy.MIN_FLOOR_SETS` | 3, 20 | `kappa_retained` entry: an excess above 3 SD of the sets' κ, from ≥ 20 sets with a finite κ |
| `decoy.MIN_CELLS`, `decoy.MIN_DETECT`, `decoy.MIN_GENES` | 20, 0.05, 5 | `kappa_retained`'s samples and gene set G: cells a sample needs, detection of a gene of G, genes of G needed |

The functions that apply these constants take them as defaults, and the tests check that too. Also fixed in code:
equal κ cells and equal depth per tumour (above) and the random streams: random sets of (c) per seed, tumour and MP;
bootstrap per seed, sample and MP; κ cells and depth thinning per seed and sample, over the sample's cells in
barcode order; `kappa_retained`'s sets per seed, sample and MP. The decoys themselves come from one stream over the
cohort's samples and cells in order, so the same cells in another order give another draw of the decoys (keep the
cohort file's checksum with the run).

The exclusions must also be able to run: with `mps: rederive`, `typer.labels` and `typer.labelled_gmt` (labels by
name cannot be registered, since the run names the MPs); with decoys on, an arm named in `decision.dosage_arms`;
with the positional axis on, a cytoband and gene positions (`positional.gtf`, or `chromosome`/`start` columns in
`adata.var`). A config without one of these is a deviation like any other. Before any NMF, `mpaudit p1` also checks
that the positions and the cytoband place genes on chromosome arms:
- a GTF that places no `var_name` (for example Ensembl IDs as `var_names`) counts as missing, and the run falls back
  to the `adata.var` columns; with none there, it is the `positional.gtf` deviation;
- if no `var_name` is placed on an arm (for example a cytoband without centromere bands for the genes'
  chromosomes), that is a `positional.cytoband` deviation, refused unless `unregistered_variant` is true;
- fewer than 50% placed only warns.

The number of genes placed goes into the log, `summary.md` and `run_config.json`.

The following are filled in by the user and registered as filled in; the code does not compare them with anything:
the arm names and contents, the pre-declared lists, the label and gene-set files, `unstable_recovery`, the `gavish`
block other than `gavish.seed`, `counts_layer` (which matrix is treated as raw UMI counts), `tme.covariates`,
`tme.n_boot`, the keys, `positional.chrom_col`/`start_col` (used when no GTF is given), `recurrence_key`, the
protocol names, the `enabled` switches of the TME and protocol axes, and every other config key. The contents of
the cytoband, GTF and GMT files are covered by the input-file digests. Switching the TME or protocol axis off
leaves its criterion unpowered, which can turn a PIVOT into an INCOMPLETE; the run logs a skipped axis and notes it
in `summary.md`. A run therefore carries the plain `P1-GO v1.3` label whenever the compared settings match. Whether
it used the registered file is checked by comparing its `run_config.json` (the normalised config) and input-file
digests with the frozen file.

`mpaudit p1` checks the whole config before any NMF runs and stops with the reason for the following (from the
command line, everything that does not need the cohort's `obs` columns is checked before the cohort is read):
- **deviations from the registered rule**, all listed with their registered values, unless the top-level
  `unregistered_variant` is true. Then it warns instead, and the verdict is labelled
  `P1-GO v1.3 (unregistered variant: setting=value, ...)` in `decision`, `rule_text`, `run_config.json` and
  `summary.md`, which opens with a bold line saying so;
- unknown keys: at the top level and in `qc`, `instability`, `decoy`, `tme`, `positional`, `typer` and `protocol`
  (a misspelt key would otherwise run the default silently; the message lists the valid keys), a section that is
  not a mapping, and unknown `decision`, `keys` and `gavish` keys (the pre-v1.0 decision keys `recruitment`,
  `partial_r2` and `protocol_frac_beyond_null` are named with their replacements);
- with the cohort given, a `keys.*` entry that names no `obs` column (for `keys.clone`, `region`, `doublet` and
  `intronic`, which only the informational typer reads, or nothing, a warning instead);
- values of the wrong type or out of range: a single name for a list key is accepted; every `enabled` switch,
  `qc.add_obs`, `instability.leave_one_study_out`, `tme.spillover_adjust` and `unregistered_variant` must be true
  or false; `decoy.n`, `decoy.n_depth_bins`, `decoy.n_null_floor`, `tme.n_null`, `tme.min_cells`, `tme.n_markers`,
  `protocol.min_pairs`, `protocol.n_null_sets` and `seed` must be whole numbers (3.0 and "3" are refused);
  `decoy.n_null_floor` must be ≥ 20 even in an unregistered variant (no sample could enter `kappa_retained` with
  fewer sets), and `protocol.n_null_sets` 0 (no random-set null) or ≥ 20 (with fewer, every tumour on the
  random-set null would be untested); `protocol.kappa_min_effect`, `decoy.min_kappa` and the informational
  `protocol.kappa_active` must be finite numbers ≥ 0; `tme.rho_sensitivity` a number > 0 or null; `tme.min_marker_fold` a number > 0;
  `min_jaccard` must be in (0, 1]; `decoy.method` and `tme.null_match` must be known choices; `protocol.reference`
  and `protocol.test` must differ;
- a missing `expected_tme_shared_gmt`, label, GMT or `positional.gtf` file, and, with the positional axis on, a
  `positional.cytoband` that is neither null (a deviation) nor a table nor an existing file (`''` included); the
  legacy `decoy.strata`; decoy arms that neither `decision.technical_arms` nor `decision.dosage_arms` names; arm
  columns missing from `obs` (a `pct_long` arm needs gene lengths) or with missing values in malignant cells (NaN,
  None, or in a text column a placeholder such as "nan", "None" or ""; the decoys would condition on values nobody
  measured; the message says how to proceed); and names that cannot match re-derived MPs (above);
- gene positions and a cytoband that place no gene on a chromosome arm (above; checked by the run once the cohort
  is read, still before any NMF).

It warns when the TME axis can run but no expected-TME-shared MP or gene set is pre-declared, when MPs from a GMT
come without `typer.labels`, and when the protocol axis is on with `keys.protocol` set but no tumour can be paired
(a protocol name not set, not a value of the protocol column, or no tumour with malignant cells under both). During
the run it warns when fewer than half of the `var_names` are placed on chromosome arms and when NMF fits were stopped
by the restart limit.

`summary.md` names the counted MPs, the MPs whose A1 statistic was not adjusted, those that pass at the estimated
ρ̂ but not at 1.5 × ρ̂, any decoy arm the rule ignores, the MPs whose `kappa_retained` is undefined in an arm of
the rule (no sample entered; `kappa_retained_n_samples` in `decoy_recruitment.csv`), the MPs meeting the gain mirror
of (c) (a κ gain in ≥ 50% of ≥ 4 tested paired tumours; informational, a rough count of what (c) gives by chance
when gains are as likely as losses), whether curated labels were given, the capped and failed NMF fits per kind
of run and the sha256 of the input files. Its notes list every axis that was switched on but did
not run (TME, protocol, positional, decoys), with the criterion left unpowered, and how many genes were placed on
chromosome arms. `mp_flags.csv` holds every per-MP column; `run_config.json` holds the normalised config, the rule
label and text, the deviations, the input-file digests (`input_files`), the fit counts (`nmf_fits`), the same notes
(`notes`) and the code hashes in the NMF cache keys (`nmf_cache_code`).

## What each axis measures

| Axis | Module | Statistic |
|---|---|---|
| **Decoys** | `mpaudit.decoy` | Per **arm** (a set of `obs` columns conditioned on jointly), the number of nuisance-only decoy cohorts (5, registered) in which the unmodified pipeline produces an MP with Jaccard ≥ 0.3 to the real MP, and `kappa_retained`: the median over samples of the excess ratio (κ_decoy − bg_decoy)/(κ_real − bg_real) on the MP's genes, bg being the median κ of 50 random gene sets matched on mean expression and detection in the same sample (drawn per seed, sample and MP), counting only samples where κ_real − bg_real exceeds max(0.02, 3 SD of the sets' κ). It is the share of the MP's own co-expression, above what every gene set of the sample shares, that the nuisance reproduces; v1.2's plain ratio κ(decoy)/κ(real) also credited the decoy with that shared background. Decoys keep depth, categorical strata (e.g. CNA clone) and continuous covariates (doublet score, contamination, %ribo, intronic fraction) cell by cell, and destroy all other co-expression. A clone arm types recruitment as dosage; a technical arm as nuisance. The example config conditions on the four technical covariates jointly in one arm, because each arm costs 5 more pipeline passes, so a technical hit is not attributed to a single covariate. To attribute it (e.g. to see whether `doublet_est`, which tracks cell size and cycling, drives a hit), declare one arm per covariate and list each in `decision.technical_arms`. Criterion (a) then counts a hit in any of them, which adds chances to hit, so fix the arm list in the registered config before running. |
| **TME co-variation, A1** | `mpaudit.tme` | Across samples, partial R² (beyond cancer type) of the malignant MP score on the same genes' pseudo-bulk score in co-resident TME cells, after **subtracting the spillover predicted** from malignant-derived ambient RNA (ρ̂ per sample and cell type from malignant-specific marker genes). Null: 500 random gene sets matched on expression and malignant specificity, put through the same adjustment (q95, q99, empirical p). With fewer than 5 markers (ambient RNA caps a marker's fold at about 1/ρ) the adjustment is skipped with a warning (`spill_adjusted` False). A second pass with every ρ̂ × 1.5 gives `partial_r2_rho`, `slope_rho` and its own null (`null_q99_rho`), because the random-set null does not absorb a leak left by a ρ̂ biased low. Also, informational only, the within-sample co-expression ratio TME/malignant (`kappa_ratio_median`, not spillover-adjusted). |
| **Protocol, A4** | `mpaudit.protocol` | Per paired tumour, Δlog κ of the MP's genes between protocols (genes detected in both; panel retention ι ≥ 0.5), with every sample of the tumour scored on the same number of cells (its smallest sample, at most 2000) and at the same depth (counts binomially thinned to the tumour's lowest median count per cell). Every paired tumour with a defined change and bounds is tested; there is no activity gate, so which tumours count does not depend on the noise of the statistic tested (the activity columns are informational). Three one-sided bounds, and a loss must beat all of them. (1) Null bound: from same-protocol replicate pairs of ≥ 3 distinct tumours, the lower of the t prediction bounds on the linear changes (carried to the MP's κ level) and on the log changes, t(0.95, number of tumours) × the root mean square of the replicate changes; otherwise the floor(0.05 (n + 1))-th smallest *linear* κ change of the tumour's n matched random gene sets (the 5th of 100, a distribution-free 5% bound for a change exchangeable with theirs), carried to the MP's κ level (on the log scale, an additive shift of near-zero background κ looked like a loss). (2) Self bound: 1.645 × the bootstrap SE of the MP's own change (50 resamples of each sample's cells), which covers the MP's cell-sampling noise that the random-set null lacks. (3) Margin: log 1.5, the only protection against regional heterogeneity of the MP on the random-set path. Intensity and prevalence use protocol-specific control genes and are tumour-relative and informational (a cohort-wide shift is absorbed). |
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
   `gene_name` values must match `var_names`: a GTF that places no `var_name` (for example with Ensembl IDs as
   `var_names`) falls back to `adata.var` with a warning, and one placing fewer than half of `var_names` warns.
   GTFs without `gene` records (UCSC `ncbiRefSeq.gtf`) are aggregated from their transcript and exon records,
   with a warning. A missing cytoband or GTF file is refused by the config check, and positions plus a cytoband
   that place no gene on a chromosome arm are a deviation found before any NMF (above).
3. To audit the *published* MPs rather than re-derived ones, set `mps:` to a GMT of the Gavish 2023 or 3CA v2 MP
   lists. Instability then reports the recovery of the matching re-derived MP. With re-derived MPs, give the
   expected-TME-shared biology as gene sets and the curated labels with `typer.labelled_gmt` (see the rule above
   for the 11 MPs removed by hand in 3CA v2).
4. Copy `configs/p1_3ca_example.yaml`. It is a template: as shipped it is refused, because its `ref/` cytoband
   and GTF files are not in the repository and `typer.labels` / `typer.labelled_gmt` are empty. Fill in the
   `decision` block (expected TME-shared list, arms), the decoy arm columns, the reference files and the curated
   labels, then freeze and register the whole file with the `mpaudit` commit (above) and run (a config that fails
   the checks that need no cohort stops before the cohort is read):

```bash
mpaudit p1 cohort.h5ad --config my_config.yaml --out results/p1 --n-jobs 8
```

**Runtime.** NMF dominates. With the default `tol = 1e-7` a fit runs several hundred iterations (about 2.5–3×
the time at `tol = 1e-5`); on a 7,000-gene × 2,000-cell sample one iteration takes about 50 ms on one core, so
ranks 4–9 × 10 restarts take on the order of 15–45 minutes per sample. Fits are cached (`cache_dir`), so threshold
and leave-one-study-out perturbations are free; seed, rank, universe, method-variant, bootstrap and decoy runs (5 per
arm) each cost one more pass. The cache key holds:
- a hash of each sample's count matrix itself (not only its total), with its barcodes, gene names and layer;
- hashes of the fitting code, taken when the package is imported, without docstrings or comments (editing a file
  during a run, or editing a docstring, does not change the keys);
- for decoy fits, a hash of the decoy generator, the decoy's arguments, the malignant barcodes and the values and
  dtypes of the columns it conditions on.

Writes are atomic, an unreadable entry is refitted with a warning, and each sample's fits are stored as soon as that
sample finishes, so an interrupted run resumes from the samples already fitted. The cache format changed in this
revision, so caches written before it are not reused. The protocol axis scores `n_null_sets` random gene sets per MP
in every paired tumour and every tumour with same-protocol replicate pieces, and bootstraps each sample's κ 50 times
per MP for the self bound. `kappa_retained` scores 50 random gene sets per malignant sample and MP on the real
counts and on every decoy; each sample's κ values come from one shared computation (`mpaudit.coexpr.set_kappas`).
On the toy this keeps v1.3 about as fast as v1.2 despite the bootstrap and the 50 sets: a toy-null job of the
calibration script takes 41–44 s against about 92 s for v1.2, and the full scripts take 15.3 minutes (with the
unequal-depth and small-piece cases added on 2026-10-18) and 7.7 minutes on 4 CPUs, against 13.8 and 8.0 for v1.2
(`scripts/calibration/README.md`). Thinning a sample to equal depth costs one binomial draw per count of its κ cells. Start with
100–200 tumours and use `--n-jobs`.

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
  `kappa_retained` ≥ 0.5 guard. On the default cohort the estimated doublet–cycle coupling is weak: `doublet_est`
  has a within-sample Spearman correlation of 0.26 with cycling (planted `doublet_score`: 0.61), a `doublet_est`
  arm keeps 0.00 of the cycle program's excess κ and the registered technical arm 0.01. With the planted score the
  technical arm keeps 0.497 of it at seed 0 (0.504 and 0.491 for decoys 1000 and 1001), just below the 0.5 guard,
  and 0.53 at seeds 1 and 2 (`scripts/calibration/reference/a_v1.3.json`). So a `biology_counted` of 0 in the
  default run does not test the guard under strong coupling; `--true-doublet` puts cycle at the guard. Its
  two-decoy value at seed 0 is 0.003 below it, and the demo averages 5 decoys, so whether cycle counts there
  depends on those decoys and on recruitment. On the null cohort, which has no TME ambient RNA,
  `contamination_est` and `doublet_est` track cycling (Spearman 0.48 and 0.51), and a `doublet_est` arm keeps 0.26
  of cycle's excess κ (see Caveats); the coupling there runs through the TME-specific genes behind
  `contamination_est`, from which `doublet_est` is built. Arm contents are registered as filled in but not compared
  by the code, so all three runs carry the same rule label; the `--true-doublet` arm is not one a registration would
  fix (the planted score cannot be observed in real data).

These are toy results at v1.3, re-run on 2026-10-18: `simulate_cohort(seed=0)` (the null cohort with
`dissociation_stress`, `tme_contamination`, `cna` and `technical_axes` off), Spearman correlations within each
sample's malignant cells (median over samples), decoys `make_decoy(..., seed=1000 + d)` for d = 0, 1 and
`kappa_retained(..., seed=1000)`, averaged over the two decoys, as `decoy_recruitment` does with config seed 0.

The toy has no curated artefacts and no curated labels, and its MPs are re-derived, so every demo run is the
unregistered variant `P1-GO v1.3 (unregistered variant: typer.labels=None)` (the demo sets
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
At seed 0 the axes now appear in 14 and 12 of the 16 samples and each forms its own MP. The cost is competition
for ranks 4–9: EMT is robust in 7 samples instead of 9, contamination in 8 instead of 11 and the CNA in 4 instead
of 8. These counts come from a fresh reference run of the demo at commit edb2736, whose NMF and simulator are
unchanged in v1.3. The calibration run had reported 13, 11 and 5, because the NMF cache of that time served one
sample (T08_fresh_R1) the fits of the first magnitudes: its cache key hashed only the sample's total count, which
was the same at both magnitudes. The key now hashes the count matrix itself. Visibility at other seeds was not
checked with NMF. With a third of the library on one per-cell factor, every other gene co-varies a little after
normalisation. For random 30-gene sets of unplanted genes detected in ≥ 5% of a sample's malignant cells, κ is
about 0.21, against 0.08 before (median over 40 sets per sample, then over the 16 samples; 0.26 against 0.10 for
genes detected in ≥ 20%). Details and the candidates tried are in the `mpaudit.sim` docstring.

Outputs besides the usual ones: `truth_comparison.csv` matches each audited MP to a planted set (Jaccard ≥ 0.2)
and lists its flags; `truth_summary.json` reports the verdict and rule label, the numbers of planted biology and
technical MPs and of those counted, biology MPs typed dosage (`biology_dosage`), biology MPs that the technical
decoys recruit in ≥ 2 of their 5 decoys whatever their `kappa_retained` (`biology_tech_recruited`, with a per-MP
list that says whether the `kappa_retained` guard passed and whether the MP is typed dosage: one that is not
counted was stopped by the guard, by a dosage label, or by both), how many counted MPs meet each criterion
(`counted_by_criterion`), the technical arm, the rule label with the deviations, and the code commit with a flag for
uncommitted changes.

RESULTS: to be filled from results/demo*/truth_summary.json after the rerun on the final code.

*Preview, not v1.3 output:* the planted-cohort demo at commit edb2736 (rule v1.2) gave GO, with the 4 planted
technical MPs counted and 0 of the 4 planted biology MPs. It is superseded by the rerun above.

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
  programs, and an MP containing ribosomal-protein genes shares counts with %ribo and keeps them mechanically (its
  `kappa_retained` can exceed 1). `kappa_retained` grades how much is explained, but a program tightly tied to a
  covariate's estimator cannot be told from nuisance here: if the decoys recruit it and keep ≥ 0.5, it counts as
  nuisance support under (a). In the calibration harness, under the technical arm with the planted doublet score,
  cycle keeps 0.497, 0.532 and 0.530 on the default toy at seeds 0, 1 and 2 (at the 0.5 guard: just below it at
  seed 0, above it at seeds 1 and 2) and 0.64–0.67 on the null toy (above it). Only
  recruitment can then block it, and recruitment was not evaluated there. A covariate computed from an MP's own
  genes pushes `kappa_retained` above 1, and the excess ratio inflates this more than v1.2's ratio did: on the
  default toy (seed 0) a `contamination_est` arm keeps 1.97 of the contamination set's excess, a `doublet_est` arm
  1.75 and the technical arm 1.39. Covariates measured independently of expression (genotype demultiplexing,
  hashing, protocol twins) settle such cases. Preview on code before v1.3 (to be replaced by the rerun): in the
  `--clean` toy run the technical arm recruited the real cycle program in 5 of 5 decoys, and only its
  `kappa_retained` of 0.35 under v1.2's plain ratio (< 0.5) kept it from counting. Under v1.3's excess ratio the
  same arm keeps 0.35 of cycle on that cohort at seed 0 (0.348 and 0.347 for decoys 1000 and 1001), so at the
  demo's seed the margin to the 0.5 guard stays about 0.15; at seeds 1 and 2 it keeps 0.07 and 0.08
  (`scripts/calibration`, two decoys).
- **Entry to `kappa_retained`.** A sample enters only if the MP's excess over its matched random sets exceeds 3 SD
  of theirs, so weakly co-expressed MPs lose samples. On the default toy (seed 0) the contamination set enters in 7
  of 16 samples (its excess is 1.7–3.7 SD; 6–15 samples over seeds 0–2), and an MP with no entering sample has an
  undefined `kappa_retained`: it is unpowered in that arm, which can hide a weak technical MP. The background SD
  comes from matched sets drawn from the sample's whole panel. On the default toy at seed 0, samples that carry
  cycle, hypoxia, EMT, IFN or CNA_7q sit 8.5–19.9 SD above their background and samples of tumours without EMT,
  IFN or CNA_7q at −0.6 to 2.8 SD, so none of the latter enters (over seeds 0–2: 1 of 22 EMT and 0 of 17 IFN
  samples without the program entered); on real data the gap may be narrower. The pair of 50 sets and 3 SD was
  chosen in the design comparison of 2026-10-17 after 20 sets failed the minority-cohort check A6 of
  `scripts/calibration` (shift of the median by samples without the program 0.082 > 0.05).
- **Spillover adjustment.** It subtracts predicted ambient spillover using ρ̂ from malignant-specific markers;
  if real ambient RNA under-represents those markers, part of the leak remains, and the random-set null does not
  absorb it. Criterion (b) therefore also needs the association at 1.5 × ρ̂; that covers ρ̂ biased low by up to
  about a third (toy with the technical axes at their first magnitudes: at 0.5 × ρ̂ a malignant-only program still
  passes), and it can push weak genuine co-variation below its null where ρ̂ is already biased high. Without 5
  malignant-specific markers (ambient RNA above about 1/8 in a TME type, or a TME type sharing the malignant lineage
  genes, e.g. normal epithelium) no adjustment is made and A1 cannot support GO. `kappa_ratio` is not adjusted and
  is informational only.
- **Protocol bounds.** Each bound of (c) is a one-sided 5% bound under its own assumption: exchangeable,
  near-normal replicate changes (replicate null; a t bound), a change exchangeable with the matched random sets'
  changes (random-set null), and cell-sampling noise (self bound). A loss must beat all of them and the log 1.5
  margin, so the 5% level per tumour holds if any one of the assumptions holds. κ is compared on equal cell numbers
  and equal depth, because its estimate falls with both (toy: at a quarter of the depth, ×0.75–0.92 for planted
  programs and ×0.48–0.60 for random 40-gene sets; a fall that scales with the set's κ, so no additive null holds
  it). Without equal depth, toy test pieces at a quarter of the reference depth made 0.12 (random-set path) and
  0.27 (replicate path) of the pairs without the program losses, and 0.06–0.07 of the carriers, and an unplanted
  program met (c) (review probe of 2026-10-18, 6 seeds); with it, 0.048 and 0.048, carriers 0.013 and 0.017, and
  no (c) (12 seeds; `scripts/calibration/README.md`). What none of the bounds covers:
  - *Regional heterogeneity of a present MP on the random-set path.* Two pieces that differ in a program's
    activity look like a protocol effect, and only the log 1.5 margin protects against it. Neither the random sets
    nor the bootstrap carry it. With small pieces, carriers are called losses now and then: 0.015 (60 malignant
    cells per piece) and 0.020 (test pieces cut to 50 cells) of the carrier pairs on the toy null's random-set
    path, against none with 200-cell pieces (12 seeds each).
  - *A patchy program on the replicate path.* The replicate null contains regional heterogeneity, but its t bounds
    assume near-normal changes. A program absent from 10–30% of the pieces gives a mixture of small changes and
    drops to the floor; every paired tumour shares one bound estimated from a few such changes, so losses cluster
    within the MP. Per tumour the rate stays near 5%, but P(c) on 4 paired tumours was 0.031–0.047 with 3 replicate
    tumours (the registered minimum) and 0.024–0.030 with 6 (synthetic, `C2-info` of `scripts/calibration`), against
    at most 0.025 for the other laws. At P(c) = 0.047, P(≥ 3 of 20 independent MPs meet (c)) would be 0.065.
  - *Background co-variation that differs between pieces* beyond cell sampling and beyond what the matched random
    sets share. Examples are piece-specific ambient RNA enriched for an MP's genes, or a protocol that removes
    shared stress co-variation. On the default toy at seed 1 (power part of `scripts/calibration/evaluate_c.py`),
    CNA_7q in T03 and EMT in T04, tumours without the program, fell from κ 0.18–0.24 to 0.05 under the real
    protocol and were called losses; no criterion (c) resulted. On the toy null the random-set path is calibrated
    at about the nominal level with little slack (0.063 of the tested pairs without the program over 12 seeds,
    0.058 over 36; target 0.065; `scripts/calibration/README.md`).
  - *Depth differences between tumours.* Equal depth acts within a tumour. The replicate null comes from other
    tumours at their own depth; on the toy, paired tumours at a quarter of the replicate tumours' depth stayed
    within target (above), but a larger gap is not covered by construction.
  - *The bootstrap is approximate.* In the design comparison of 2026-10-17 its SE was 0.90–1.11 of the true SD of
    κ for toy gene sets without a program, at 60 and 200 cells. 50 resamples leave about 10% noise in the SE, and
    the bootstrap ignores dependence between cells within a piece (doublets, ambient RNA).

  Toy run (`simulate_cohort()` defaults, seed 0; the planted sets; `protocol_shift_kappa(adata, keys, sets,
  "fresh", "fixed")` with its defaults, 100 random sets per tumour and MP; re-run on 2026-10-18 with equal depth).
  All 32 κ tests of the eight planted programs other than stress are tested, 3 of them in a paired tumour that
  lacks the program (CNA_7q in T01 and T03, IFN in T03). Under each set of bounds:
  - the null bound alone: 5 losses and 7 gains;
  - the null bound and the self bound: 3 losses and 3 gains;
  - the null bound and the margin: 1 loss, CNA_7q in T03, a tumour without the CNA (Δlog κ −0.611 against a null
    bound of −0.550);
  - all three, the registered rule: none.

  The self bound is −0.675 there. Stress, which the fixed protocol removes, is a loss in 4 of 4 paired tumours. The
  replicate null needs same-protocol pieces of ≥ 3 tumours, which public data rarely have. Oliveira et al. 2025 is
  not a fresh-versus-fixed twin (FFPE Flex versus in situ on the same blocks).
- **What (c) can and cannot detect.** Without the gate, tumours that lack a program count in the denominator of
  (c). An MP present in fewer than half of the paired tumours therefore cannot meet (c), however strongly it is
  lost where present. Such subset nuisance programs can reach support only through (a) or (b). This is conservative
  for GO, but it under-detects real subset nuisance. In the calibration harness a planted κ drop of about 36–48%
  in EMT, in the paired tumours that carry it (about half), met (c) in 0, 0, 2 and 2 of 6 runs (two drop sizes ×
  two null paths; v1.2: 2, 0, 5 and 3). On the replicate path, EMT was detected in 0.18 of the injected tumours at
  the smaller drop and 0.53 at the larger one (random-set path: 0.47 and 0.83); replicate pairs at the floor widen
  the log-scale bound. Equal cells and equal depth cost power when one piece is small or shallow: the whole
  tumour is then scored at that piece's size and depth, so its test is noisy and rarely a loss.
  And the power clause of INCOMPLETE for (c) is met by any cohort with ≥ 4 paired tumours where an MP has defined
  bounds, so on shallow data PIVOT is more likely than INCOMPLETE.
- **False (c) rates (toy and synthetic only).** From `scripts/calibration/evaluate_c.py` (v1.3 reference results in
  its README):
  - Synthetic κ laws with no protocol effect, MP and random sets (or replicate pieces) independent: 0.018–0.054 of
    tested tumours are losses when the MP is drawn like its random sets (v1.2: 0.05–0.31 on the random-set path),
    and 0.005–0.048 when it is present.
  - For an unaffected MP, P(criterion (c)) is at most 0.025 with 4 paired tumours (worst with 3 replicate tumours;
    at most 0.019 on the other paths), at most 0.003 with 10 and 0 with 40. If the M audited MPs were independent,
    each with P(c) = 0.025, P(≥ 3 of them meet (c) by chance) would be 0.013 for M = 20, 0.038 for M = 30 and
    0.083 for M = 41; at P(c) = 0.003 it is ≤ 0.0003. False (a) and (b) support add to these chances, and real MPs
    are not independent.
  - Toy null (pseudo-paired pieces with no protocol effect; 12 seeds × 11 gene sets at both toy magnitudes; 10
    paired tumours on the random-set path, also split into halves of 5, and 5 paired plus 5 replicate tumours on
    the replicate path): 0.019–0.063 of the tested pairs without the program and none of the pairs with it are
    losses, and no gene set met (c). With the test pieces at a half or a quarter of the reference depth: 0.022–0.050
    and at most 0.017 of the carriers; with 50–60-cell pieces: 0.022–0.043 and at most 0.020 of the carriers; no
    gene set met (c).

  These rates hold under the stated assumptions. They do not cover regional heterogeneity of a present MP on the
  random-set path (above), for which the margin is the only protection (for a weak MP, κ 0.05–0.2, a 1.5-fold
  change is only 0.02–0.07 in κ), nor a patchy program on the replicate path (P(c) up to 0.047 with 3 replicate
  tumours, above).
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
- **Toy cohort and validation.** `mpaudit.sim` results are true by construction. They show that the code runs and
  how each axis responds; they are not evidence about real tumours. The v1.3 rule was calibrated on the toy and on
  idealised synthetic κ laws only (`scripts/calibration/`). The synthetic part cannot see the self bound, the
  equal-cell rule or equal depth, because synthetic κ levels have no cells (the toy part covers them, with toy
  pieces at a half or a quarter of the reference depth and with small pieces), and recruitment (the NMF hits of
  (a)) was not evaluated there. The toy's technical axes were recalibrated on 2026-10-17 (above), so seed-0
  numbers measured before that date differ from the current toy, and every random draw of (a) and (c) changed in
  v1.3, so toy outputs of v1.2 do not carry over.
