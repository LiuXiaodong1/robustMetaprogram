# TESSERA: a proposal for calibrated tumour metaprograms with typed provenance

**Draft project proposal (AI-assisted design synthesis; to be reviewed by the PI). Target: *Nature Methods* Article.**
Date: 2026-10-02. Inputs: the landscape and gap analysis ([01_landscape_and_gap_analysis.md](01_landscape_and_gap_analysis.md)), four method concepts (CREDO, CALIPER, InvarMP, CERT-MP), and eight adversarial referee reports ([03_design_panel.md](03_design_panel.md)).

**Revision note (2026-10-17; item 15 added 2026-10-18).** Revised in response to the review of 2026-10-16 (code, citations and referee outputs of commit 7f4867e). What changed and why:
1. **One P1 GO rule.** The text gave three different P1 gates, and the coded one could pass on almost any cohort. §8 states one frozen rule, now **P1-GO v1.3** (version history in §8), as implemented in `src/mpaudit/report.py`; §0 and §7 summarise it. Dosage and instability no longer count toward GO. Decoys run through the Python port, not the original R code, and this is now said.
2. **Null-relative thresholds.** The fixed partial R² bars of P1 and P6 sat below the random-gene-set null on the toy (null q95 0.35–0.66 in the reviewed demo, up to 0.78 in the current toy). They are replaced by "above the matched-null q99" (§2.5, §8).
3. **A1 spillover.** Malignant→TME ambient RNA is now modelled and subtracted from the pseudo-bulk TME score; the within-sample κ ratio is not adjusted and can no longer make A1 "supported" on its own. The claim that TME cells carry no malignant programs is tempered (§1.3, §2.5).
4. **A4.** The axis is κ-based with protocol-specific controls, and intensity and prevalence are tumour-relative and informational. Oliveira 2025 profiles FFPE blocks and has no fresh-dissociated malignant cells, so it is no longer listed as a fresh-versus-fixed twin (§2.5, §6.4).
5. **Stress.** Baron 2020 was mis-cited. Dissociation-induced expression is attributed to van den Brink 2017, O'Flanagan 2019 and Neuschulz 2023, and the novelty is restated as per-tumour attribution (§7).
6. **Specificity claims.** §1.3 and §6.3 now report what the saved simulations show, including the new `specificity_fix_check.py`. The 72%/100% figures are a specificity failure of the κ > 0 estimand, not a type I error.
7. **New nuisance axes.** Gene length / nuclear (intronic) fraction and transcript abundance / %ribo are added, on user-reported evidence (§1.1, §2.1, §2.5, §6.2).
8. **Pilots, timeline and pre-registration.** The pilot schedule now spans months 0–6 (P4 and P5 possibly to month 9), with summed durations. P1-GO v1.3 and the whole filled-in `mpaudit` config (with the `mpaudit` commit and the curated-label files) are registered on OSF at month 0, together with an audit/confirmation split of studies (§8, §9, §11).
9. **Positioning and references.** Updated for starCAT, TCCA, ISOMERIC, cEBMF, DECIPHER, SIMVI, arXiv 2410.06451, the Gene Program Dictionary and scPFA, each with its verification status (§4, §12).
10. **Implementation status.** Covers `mpaudit` (P1), the NMF tolerance finding and the R-parity harness (§8).
11. **Rule implementation tightened (code review, 2026-10-17).** Criterion (b) counts only spillover-adjusted statistics that still hold with the ambient fraction ρ̂ scaled up by 1.5. The A4 nulls are recalibrated: a t prediction bound from replicate pairs of ≥3 distinct tumours, random-set κ changes carried to the MP's κ level, and a log 1.5 minimum effect. A dosage decoy hit needs no minimum number of decoys, as the frozen rule says. Labels and pre-declared exclusions are carried onto re-derived MPs by gene overlap, and the decision config is checked before any computation (§2.5, §8). Because these change the rule, its label moved from P1-GO v1.0 to P1-GO v1.1.
12. **Activity gate and enforced settings (2026-10-17).** In criterion (c) a paired tumour was tested, and a replicate pair entered the null, only where the MP was active in at least one of the two pieces compared (κ above max(0.05, the q95 of matched random gene sets), a simplified form of the expression threshold κ*_t of §1.2). The gate was removed again in v1.3 (item 15). The registered settings outside the `decision` block (seed included) must equal their registered values, and a run in which an exclusion cannot run (re-derived MPs without curated labels, no dosage decoy arm, no cytoband or gene positions) is a deviation too; a deviation is refused unless the run is declared an unregistered variant, which labels its verdict as one. This made the rule **P1-GO v1.2**.
13. **Toy simulator.** Its abundance/%ribo and long-gene axes were too weak for the standard pipeline to see. They were recalibrated against a visibility target fixed before any audit was run (§8, implementation status). Toy numbers below come from the recalibrated toy unless marked otherwise.
14. **Frozen records.** The design-panel record ([03](03_design_panel.md)) and the raw literature sweeps ([04](04_literature_sweeps.md)) are not revised; an erratum at the top of each points to the corrections (TCAT/starCAT scope, Oliveira 2025).
15. **P1-GO v1.3 (2026-10-18).** A design comparison found that the v1.2 activity gate selected tumours on the noisy κ that it then tested. A tumour whose κ was high by chance passed the gate and then tended to drop. In the calibration harness (`scripts/calibration/`), up to 0.31 of tested tumours were losses on synthetic laws with no protocol effect (0.05–0.31 across the exchangeable laws), 0.13–0.17 of tested pairs without the program on the toy null, and an unplanted random gene set met (c). It also found that v1.2's `kappa_retained` let in samples where the MP is absent: a biology MP kept in 2 of 10 toy tumours reached 0.66 > 0.5. In v1.3:
    - (c) has no activity gate;
    - the replicate bound is the lower of the linear and log t bounds;
    - the random-set bound is a distribution-free order statistic;
    - a loss must also beat 1.645 × the bootstrap SE of the MP's own change;
    - every sample of a tumour is scored on the same number of cells and at the same depth (binomial thinning to the tumour's lowest median count per cell), with the panel taken on those counts;
    - an MP's random sets are drawn per tumour and MP, and no draw depends on the other MPs audited or on the cohort's cell order;
    - `kappa_retained` is the excess over 50 matched random gene sets, with a 3-SD entry.

    Equal depth, the panel on the thinned counts and the order-free draws were added after a third review of the same day, still before registration: test pieces at a quarter of the reference depth had made 0.12 of the toy null's absent pairs losses and a false (c). On the harness every calibration check of (c) passes, including toy pieces at a half or a quarter of the reference depth; the gaps that remain are listed in §8 and in [05_pilot_P1.md](05_pilot_P1.md), Caveats. (a) passes all its checks except the known coupled-covariate limit (§2.6). These are toy and synthetic results only. A second code review added `gavish.seed` to the registered settings, a refusal of unknown config keys and missing reference files, a check before any NMF that gene positions and the cytoband place genes on chromosome arms, and an NMF cache key that hashes each sample's count matrix. Toy numbers of criteria (a) and (c) below were re-run at v1.3 or are marked as earlier previews; the other axes did not change. v1.3 is the version to register; every mention below refers to it (§8).

Numbers from the toy simulator (`mpaudit.sim`) are labelled as toy results. They are true by construction and are not evidence about real tumours.

**Working name: TESSERA** (Typed, Effect-Sized, Specific, Split-validated, Externally Replicated, Audited metaprograms). A tessera is one tile of a mosaic.
- Python package `tessera-mp` (`import tessera`); R bridge `tesseraR`.
- The name has **not** been checked for collisions. Fallback names are `mpcert` and `CREDO`.

---

## 0. Bottom line

- **Can it reach *Nature Methods*? Yes, but only if three conditions hold.**
  1. Pilot P1 returns **GO** under the frozen rule **P1-GO v1.3** (stated in full in §8; implemented in `src/mpaudit/report.py`). GO needs at least 3 *retained* 3CA MPs (not curated artefacts, not typed dosage) with nuisance support from at least one of:
     - a technical decoy arm (≥2 of the arm's 5 decoys recruit the MP at Jaccard ≥0.3, with κ_retained ≥0.5, κ_retained being the share of the MP's co-expression in excess of matched random gene sets that the decoys keep);
     - a positive, spillover-adjusted TME association above the q99 of matched random gene sets in ≥10 samples with TME cells that still holds with the ambient fraction ρ̂ scaled up by 1.5, for MPs not pre-declared as expected TME-shared;
     - a one-sided κ loss under a dissociation-free protocol in ≥50% of ≥4 tested paired tumours (every paired tumour with defined bounds is tested; no activity gate; κ compared on equal cell numbers and equal depth). The loss must beat three bounds: the null bound (95% t prediction bounds from same-protocol replicate pairs of ≥3 distinct tumours, else an order statistic of matched random gene sets); the MP's own sampling noise (1.645 bootstrap SEs); and 1.5-fold.

     Rediscovering the 11 MPs already removed by hand does not count: they are excluded through the registered curated labels (`typer.labels` / `typer.labelled_gmt`), whose gene sets for a re-derivation run still have to be rebuilt from the 3CA v2 code (§8); until then a re-derivation run is an unregistered variant. A dosage-typed MP never counts, whatever its other support (it is genetic, not technical nuisance). Instability is informational and is not nuisance support. If no criterion is adequately powered, the verdict is INCOMPLETE rather than GO or PIVOT.
  2. We secure a matched multi-protocol tumour cohort of **≥20 tumours**, either through a wet-lab partner or, less convincingly, by pooling public "twin" datasets.
  3. The certification layer measurably beats simple baselines: module preservation, count splitting plus BH, and a trivial covariate-based typer.

  My estimates of the odds (judgement, not data):
  - all three conditions met: about 35–45%;
  - computational-only route (condition 2 fails): about 15–20%;
  - a strong fallback paper (*Genome Biology*, *Nature Communications*, *Cell Systems*, or a *Nature Methods* Analysis) in either case: above 80%.
- **One-sentence claim.** TESSERA turns a tumour metaprogram from a threshold-defined top-50 list into a falsifiable estimand with four properties:
  - it has a stated minimum effect;
  - it replicates on held-out cells and held-out cohorts;
  - it is specific beyond every other program;
  - its provenance (malignant-intrinsic biology, biology shared with the microenvironment, genetic dosage, protocol-modulated, or technical) is assigned by **equivalence tests against real negative controls**, never by failing to reject.

  The package wraps any factorisation engine and applies this to the pan-cancer 3CA corpus.
- **The core is a hybrid.**
  - Spine: an engine-agnostic, cross-fitted certification layer, taken from CALIPER and CERT-MP with their statistical flaws repaired.
  - Default engine: a reduced count-native engine taken from CREDO.
  - Protocol axis: InvarMP's invariance idea, applied only where a design bridge exists.
  - Sequencing: CERT-MP's "audit first", which gives a publishable fallback by month 10.
- **What I dropped deliberately**, because referees showed each one is circular, invalid or unfundable at this team size:
  - structureless or parametric decoys as an FDR device;
  - knockoff+ over tumour swaps;
  - "error-controlled existence" against nuisance-only nulls;
  - first-failure typing;
  - universal-inference e-values;
  - expression-only CNA *removal*;
  - the power posterior as a headline.
- **Why it is hard to scoop.** The claim does not depend on any model. If a count-based GBCD, a single-cell port of Grabski's model, or a malignant TCAT appears, it becomes one more engine inside TESSERA. Three assets are hard to copy quickly:
  - within-sample negative controls from co-resident microenvironment (TME) cells;
  - order-invariant equivalence typing;
  - a new matched multi-protocol cohort with pre-registered predictions.
- **Biggest risk.** Most canonical MPs (cycle, IFN, hypoxia, pEMT) may simply be real. The audit would then mostly confirm them, and the paper reads as "a more complex method reproduces known programs with q-values".
  - Mitigation: pilot P1 is a hard go/no-go gate. Its rule (P1-GO v1.3, §8) is frozen and registered before any 3CA expression data are opened, and it can fail.
  - The flagship is built around results that are *quantitative per tumour*, not around overturning MPs: the ex vivo share of each MP in each tumour, selection versus induction, and malignant-intrinsic versus TME-shared.
- **Resources.** About 4 FTE for about 26 months to submission; 2,000–4,000 GPU-hours; a wet-lab cohort costing about USD 300–450k at full scope (an estimate to be quoted), with a minimal version at about USD 150–200k.

---

## 1. The conceptual advance

### 1.1 What is wrong with "recurrence"

In the Gavish/3CA pipeline, "robust" means "recurrent". Its cross-tumour filter (≥10 of 50 shared genes from a 7,000-gene universe) guards against a random-gene null with P ≈ 1e-12 (landscape §1.4). Anything that recurs for a systematic reason passes: ambient RNA, doublets, dissociation stress, quality axes, gene-length / nuclear-fraction and transcript-abundance (%ribo) axes, positional programs. Every published variant then ends in manual curation (3CA v2 `inds_rm`, GeneNMF `dropMetaPrograms`).

The referees showed that simply swapping in a *structureless* null does not fix this. Nuisance-only decoys reject for *any* co-expression, biological or artefactual. At atlas sample sizes every candidate reaches the permutation floor (CREDO referees, fatal flaw 1; CERT-MP referees, flaws 1 and 3).

### 1.2 The redefinition

TESSERA splits the user's complaint ("not robust, sometimes not real") into two questions that need different evidence:

- **Q1, statistical reality.** Is the structure real rather than an artefact of overfitting noise or selection? This is answered by **cross-fitting** on disjoint cells and cohorts, with an **effect-size** null.
- **Q2, provenance.** If it is real, what produces it? This is answered by **real negative controls and design contrasts**:
  - co-resident TME cells in the same sample;
  - the same tumour profiled under different protocols;
  - allele-anchored copy number;
  - empirical regional variation.

  The answer is a provenance profile with explicit "indeterminate" outcomes. Simulated nulls do not answer Q2.

**Formal definition.** Setting:
- cells c are nested in tumours t, tumours in cohorts (study × protocol), cohorts in lineages l;
- UMI counts y_cg;
- every tumour's cells are randomly partitioned into three folds: discovery 𝒟 (40%), membership ℳ (30%) and gene-effect 𝒢 (30%).

A candidate program k is a non-negative loading vector w_k learned *only* on 𝒟 cells, by any engine.

**Definition (TESSERA metaprogram at level α).** A candidate k is an MP if all three of the following hold.

1. **Specific tumour expression with a minimum effect (estimand).** For each tumour, compute on ℳ cells
   κ_tk = Σ_{g≠g'} w_g w_g' R^x_t,gg'(· | 𝒩, a_{−k}) / Σ_{g≠g'} w_g w_g'

   where R^x is the gene–gene correlation of *latent* expression. It is estimated with a measurement-error-corrected moment estimator in the style of CS-CORE (Poisson/NB measurement model with size factors). The estimand does not depend on sequencing depth; the moment estimate does: it falls as counts per gene fall, more for weakly co-expressed sets (toy: ×0.75–0.92 for planted programs and ×0.48–0.60 for random 40-gene sets at a quarter of the depth). κ values from samples of different depth are therefore compared only after thinning to equal depth (P1 does this, §2.5).

   The correlation is computed after removing three things:
   - the declared nuisance 𝒩: depth, the ambient profile, and same-sample TME-profile components;
   - the activities a_{−k} of **every other** candidate program;
   - the tumour's private programs.

   Tumour t **expresses** k if we reject H0: κ_tk ≤ κ*_t, where κ*_t = max(κ_min, q_0.95 of κ over random gene sets in the same tumour). The random sets are matched on mean-expression decile, detection rate, chromosome-arm clustering and weight profile. The default κ_min is 0.05, pre-registered.

   For each tumour we also report a prevalence π_tk (fraction of ℳ cells with activity above the matched-random 99th percentile) and an intensity m_tk.

2. **Recurrence that is replicated, not counted.** Both of the following:
   - a partial-conjunction lower confidence bound r_L(k) ≥ r_min (default 5) on the number of expressing tumours, spanning ≥2 cohorts;
   - **leave-cohort-group-out replication**: when k is rediscovered without cohort group G (matched at a fixed, pre-registered cosine tolerance of 0.7), its κ is above the null band in G's tumours.

   Prevalence per lineage is estimated by empirical-Bayes deconvolution of the (κ̂_tk, SE_tk) pairs. It is not a count of rejections, because rejection counts scale with cell number.

3. **Typed provenance.** Each k carries an evidence vector over five axes (§2.5). Every axis returns *supported*, *excluded* (by equivalence) or *indeterminate*. A program is "certified intrinsic biology" only on **positive** evidence: equivalence on the technical axes, plus invariance on ≥1 bridged protocol axis. Lack of power yields *indeterminate*, never "biology".

### 1.3 Why this is a new idea and not "another NMF"

- **The estimand does not depend on the estimator.** The MP is (w_k, {κ_tk}, gene effects, provenance), and validity does not depend on how w_k was found. CALIPER's referees showed that its location estimand changed with the engine. Here the location is whatever the engine proposes; the *certificate* is engine-free, and engine invariance is reported as an empirical result, not claimed.
- **Membership is defined by effect size and specificity, so it does not inflate with sample size.** The toy checks are in [referee_checks/](referee_checks/README.md); their limits are listed at the end of this item.
  - *The failure, ref1.* [ref1_check.py](referee_checks/ref1_check.py) uses 1,000 cells, NB size 3 and an unrelated program w at fold change 2.5 in 30% of cells. A point-null cross-fold test of κ(m) > 0 declared that a tumour "expresses" m in 72% of tumours when m shares 5 of its 40 genes with w, and in 100% at 10 of 40.
  - *How to read it.* κ(m) > 0 is **true** in that setup, because the shared genes really co-vary through w. The 72% and 100% are therefore a specificity failure of the "κ > 0" estimand, not a type I error. The type I error is the k = 0 row, at 0.05.
  - *A weaker program, ref2 sim2.* The same failure appears and grows with the number of cells. Pooled over 3 seeds, false membership was 10–17% at 4/40 shared genes and 28–89% at 8/40, from 300 to 3,000 cells ([ref2_sim2.out](referee_checks/ref2_sim2.out) and seeds 1 and 3). Standardising by random gene sets did not help (8/40, pooled: 0.28 → 0.29, 0.50 → 0.54 and 0.89 → 0.89).
  - *The proposed fix, tested.* [specificity_fix_check.py](referee_checks/specificity_fix_check.py) tests the estimand of definition 1 (200 replicates per row). Conditioning on the overlapping program kept false membership at ≤ 1% up to 30 of 40 shared genes, with power 1.00. The within-tumour matched-random null band contributed nothing. In non-member tumours its q95 stayed below the floor with total-UMI normalisation and with varied baselines. With the true depth and flat baselines it rose above the floor at 20 and 30 of 40 shared genes (mean q95 0.061 and 0.075), because mean-matched sets are drawn from w's genes in that setting. It still stayed below κ(m) (0.087, 0.199), so false membership was 0.98 and 1.00. Without conditioning, only the κ_min = 0.05 floor protects, and it fails once the overlap pushes κ above 0.05. With the true depth that happens from about 20 of 40 shared genes, where false membership is 98–100%.
  - *Limits of that check.* It used a single tumour, a *perfect* dictionary (true loadings and baseline), a 300-gene universe and a membership fold of all 1,000 cells. Power was checked only at fold change 2.5. Specificity therefore rests on the overlapping program being in the dictionary, and P2 must show that with an estimated dictionary.
- **Real negative controls replace synthetic nulls for provenance.** TME cells in the same droplet suspension share dissociation time, ischaemia and RNA quality with the malignant cells. They do not *express* malignant-specific programs, but they are not clean of them. As the 2026-10-16 review argued, most ambient RNA in tumour suspensions is likely to come from malignant cells (not quantified here). A malignant program whose strength varies between samples therefore leaks into the same samples' TME cells and can look TME-shared.
  - A1 subtracts this predicted spillover from the pseudo-bulk TME score (quantity ii) before using TME cells as controls; the within-sample κ ratio (i) is not yet adjusted (§2.5).
  - Some programs (IFN response, hypoxia, stress) are genuinely expressed in both compartments. TME co-variation alone therefore never types a program as technical in TESSERA. (P1 is cruder: its criterion (b) counts TME co-variation alone as nuisance support, except for MPs pre-declared as expected TME-shared; §8.)

  With these caveats, the control exists in essentially every 10x tumour sample. It turns "is this technical?" into a measured quantity rather than a curator's judgement.
- **Provenance is decided by equivalence tests, with an "indeterminate" outcome.** This reverses the burden of proof that the CALIPER and InvarMP referees found fatal.
- **The conceptual inheritance is acknowledged, and so is the gap.** Cross-fitting and synthetic-control ideas exist (count splitting; ClusterDE; Clipper; sc-SHC; module preservation). What has not been done is to define and type a **cross-tumour recurrent program** with these tools at pan-cancer scale.

---

## 2. Method specification

### 2.1 Inputs and nuisance (precomputed, imported)

Required inputs:
- raw integer UMI counts;
- malignant labels;
- **TME cells retained in the same object**;
- tumour, sample, study, protocol and lineage keys.

Optional inputs:
- empty-droplet matrices (CellBender or SoupX ambient profile; otherwise a DecontX-style proxy, flagged as lower quality);
- scDblFinder scores;
- per-cell intronic (unspliced) fraction from STARsolo or velocyto (it needs raw reads: 3CA matrices carry no spliced/unspliced counts), and gene lengths from a GTF (for the gene-length and nuclear-fraction axis, A5);
- Numbat allele-aware CNA calls where BAMs exist;
- pre-analytic covariates: ischaemia time, digestion time, DV200.

An **EnvSchema** (taken from InvarMP) records, per sample:
- preservation;
- molecule (cell, nucleus or in situ);
- dissociation (warm enzymatic, cold protease, fixed-then-dissociated, none);
- chemistry;
- processing pipeline and version (e.g. Cell Ranger ≥7 counts intronic reads by default);
- site.

Smart-seq2 read-count studies cannot be count-thinned. They are routed to cell-split certification only and get their own dispersion.

### 2.2 Default estimation engine (TESSERA-NMF, a reduced CREDO)

Mean model for malignant cell c in tumour t, sample s, cohort u:

μ_cg = s_c · D_cg · [ b_tg + Σ_k θ_ck W_kg + Σ_r θ^U_cr U_ur,g + Σ_j θ^V_cj V_tj,g + Σ_e π_ce M_se,g ] + s_c ρ_c A_sg

Each term:
- **b_tg ≥ 0** is a tumour baseline that replaces centring. The primary estimand is therefore within-tumour programs, as in Gavish. Between-tumour structure is handled by an optional module (§2.7).
- **W** (K_max = 200, with automatic relevance determination (ARD), rows on the simplex) is the shared dictionary.
- **U** is a cohort tier (4 slots per cohort) that absorbs lab-specific axes.
- **V_t** is a tumour-private tier (3 slots per tumour).
- **M_se,g** are fixed mean profiles of the same sample's annotated TME types. Their usage π_ce ≥ 0 is unconstrained, which covers heterotypic doublets *and* whole TME cells misclassified as malignant. CREDO's referees noted that a Beta(5,5) doublet prior cannot represent these. The L1 penalty on π is weighted by 1 / (calibrated doublet score).
- **A_s, ρ_c** are the ambient profile and the per-cell ambient fraction.
- **D_cg** is a dosage offset applied **only to allele-supported segments**: D_cg = (CN_c,seg / P_t)^{γ_g,l}, with gene-specific, lineage-shrunk γ. It equals 1 elsewhere. Expression-only CNA calls are never used as offsets, only for typing (§2.5).
- **QC covariates are not offsets.** %mito, %ribo and complexity are downstream of biology (cell cycle, OXPHOS, translation), so they are used only as a reported "QC-entangled" axis. (P1's decoy audit is an exception; see §2.5 and §2.6.)

Objective, minimised by MAP with GPU multiplicative updates (Poisson/KL working likelihood) on tumour-stratified minibatches:

L = Σ_t w_t Σ_{c∈𝒟_t} Σ_g KL(y_cg ‖ μ_cg) + ARD(θ, θ^U, θ^V) + λ_x Σ_t Σ_j ‖V_tj − Π_cone(W) V_tj‖^c

Here Π_cone(W) is the NNLS projection onto the cone of W. The exclusion term is written as a penalty on the *in-cone* component of private loadings, applied to **all** tumours (not only members, which fixes a referee point) and implemented as the NNLS-fit share. w_t = min(1, 2000/n_t) stops very large tumours from dictating the dictionary.

- Initialisation: NNDSVDa (deterministic).
- Model complexity (K, λ_x, regularisation strength ζ) is chosen by held-out deviance on **molecular cross-validation**, i.e. binomial thinning of 𝒟 cells. Thinning is used here **only for model selection**, where small fold dependence is tolerable, never for inference.
- An optional lineage-specific loading perturbation W_k,l = W_k ⊙ exp(ξ_k,l), with shrinkage, can be switched on. It answers the "one W across lineages" concern and is reported as a sensitivity analysis.

**What the engine does *not* claim.** It claims no posterior, no power-posterior junk suppression and no identifiability theorem. All inferential claims come from §2.3–2.6. The engine is replaceable. Built-in adapters:
- Gavish (exact tiroshlab R code via rpy2);
- cNMF;
- GeneNMF;
- consensus scHPF;
- GBCD (signed output mapped to its positive and negative parts, documented as a transformation);
- DRVI;
- scVI+Hotspot.

### 2.3 Cross-fitting design (the repair for double dipping and thinning validity)

- **Disjoint cell folds per tumour.** Discovery uses 𝒟, membership uses ℳ, gene effects use 𝒢. Under within-tumour cell exchangeability, ℳ and 𝒢 are independent of the selection made on 𝒟.
- **Dispersion estimates are not needed for validity.** This removes the Hivert et al. objection, the asymmetry between decoy and real thinning, and the false "sign control" corollary. It also works for Smart-seq2 read counts.
- **Membership and genes use different folds,** which removes CALIPER's fold-C double use.
- **Cost:** about 40% of cells per tumour for discovery, and 30% each for testing. Tumours with fewer than 60 cells in a test fold are reported as "insufficient", never as negatives.
- **Primary inference uses one pre-specified split.** Five further random splits feed a **stability index** (the fraction of splits in which k is certified at matched cosine ≥0.7). Derandomised e-value averaging is only a sensitivity analysis, because its hypothesis definition is delicate.
- **Cell exchangeability is an assumption.** Sample-level shared signal (ambient, dissociation) is real, not noise. It replicates across folds, which is exactly why Q2 exists separately from Q1.

### 2.4 Membership, recurrence and gene membership

- **Tumour test.** κ̂_tk with a cell-bootstrap SE; a one-sided test against κ*_t; BH across tumours within each MP. Specificity comes from removing all other programs' ℳ-fold activities (NNLS with W fixed) and private programs before computing R^x.
- **Prevalence.** Empirical-Bayes deconvolution of κ̂_tk across tumours, per lineage, gives the estimated fraction of tumours with κ > κ*. This is a power-adjusted quantity that is comparable between deep 10x and shallow single-nucleus studies.
- **Recurrence bound.** A partial-conjunction lower bound r_L. Simes is applied within cohorts and assumes positive regression dependence (PRDS); that is stated explicitly, and a Bonferroni–Simes sensitivity analysis is reported.
- **External replication.** Studies are grouped into 10 folds, stratified by cancer type, rather than leaving out each of ~150 cohorts. This makes the cost tractable.
- **Gene membership** on 𝒢 cells:
  1. compute a leave-gene-out activity ã^{(−g)} from fixed W;
  2. correct it for split-half reliability (against regression dilution);
  3. fit an NB GLM of y_cg on ã with offset log(s_c μ^{(−k)}_cg);
  4. pool tumour effects by three-level REML random-effects (tumours within studies) with Hartung–Knapp–Sidik–Jonkman (HKSJ) intervals;
  5. test a **minimum-effect null** H0: β_g ≤ log 1.2 per SD of reliable activity;
  6. apply BH within each MP.

  Outputs: core genes (q < 0.05, I² < 0.5), context genes (q < 0.05, I² ≥ 0.5) and cancer-type-specific effects. Gene-list length then follows effect size, not cohort size. We will show that core-list length is flat as a function of the number of cells.

### 2.5 Provenance typing: five evidence axes and equivalence decisions

Each axis returns an estimate with a 90% CI. The pre-registered margins below turn it into supported, excluded or indeterminate.

| Axis | Estimand | "Excluded" (equivalence) | "Supported" |
|---|---|---|---|
| **A1 Compartment sharing** (the within-sample negative control) | Two quantities: (i) κ of k's genes in co-resident T/NK, myeloid, fibroblast and endothelial cells, after removing TME identity programs; (ii) partial R² of the per-sample malignant score on the same gene set's **pseudo-bulk TME score after subtracting predicted malignant→TME ambient spillover** (see below), conditioned on cancer type with study random effects, against random gene sets matched on expression **and** malignant specificity | Upper CI of the TME/malignant κ ratio < 0.3 **and** upper CI of the adjusted partial R² below the matched-null q95 | An adjusted partial R² above the matched-null q99 with a positive slope. The unadjusted κ ratio is reported, but it cannot by itself make A1 supported until it is spillover-adjusted. There is no fixed R² bar (see below) |
| **A2 Contamination/doublet** | Shapley share of κ_tk removed by the M_se and A_s terms, **together with** the cell-level association of k-activity with doublet score and ambient fraction | Share upper CI < 0.2 | Share lower CI > 0.5 **and** a positive doublet/ambient association. Loading similarity to TME cells alone is never enough, which protects mimicry programs such as pEMT, MES-like and MHC-II |
| **A3 Genetic** | (a) dosage share: Shapley share of κ removed by **allele-anchored** D_cg; (b) clone η² of activity with a region-aware block permutation; (c) arm enrichment, with known positionally clustered biological programs as negative controls (EDC 1q21, KRT 12q13/17q21, MHC-II 6p21, histones 6p22, metallothioneins 16q13) | Dosage share upper CI < 0.2 (only where allele data exist) | Dosage share lower CI > 0.5 → "dosage-driven"; η² lower CI > 0.3 → "clone-linked". With expression-only CNA, only a "positional flag" is possible |
| **A4 Protocol** | Paired, same-tumour **Δlog κ** between protocols, on the gene panel both protocols retain (panel retention ι_k ≥ 0.5), in tumours that express k under at least one of the two protocols (κ above κ*_t, §1.2). Selecting tumours by a κ*_t test on the same cells and then testing their Δlog κ carries the selection bias that P1-GO v1.2's activity gate showed (below): tumours selected on a κ high by chance tend to drop. The full method therefore needs cross-fitting (select on one cell fold, test on another) or a test conditional on the selection. P1 tests all paired tumours (below). κ's estimand is a correlation, so gene-wise capture efficiency does not move it; this sidesteps InvarMP's κ-versus-δ identifiability problem. Its estimate still falls when counts fall (lower capture or depth), so P1 compares κ at equal depth (below). Null: same-protocol, multi-piece replicates (regional heterogeneity) where pairs from ≥3 distinct tumours exist, as a one-sided 95% t prediction bound; otherwise random expression-matched gene sets in the same paired tumours, whose κ change is carried to the MP's κ level on the linear scale. Δπ (prevalence) and Δm (intensity) use protocol-specific control genes from each protocol's malignant cells in the paired tumours; they are tumour-relative and informational (see below) | \|mean Δlog κ\| ≤ log 1.5 by two one-sided tests (TOST) | κ loss beyond the null with lower CI of loss > log 2 in fixed or in situ data → ex vivo; a tumour-relative Δπ with stable Δm → candidate selection (informational) |
| **A5 QC entanglement** | Correlation of activity with %mito, complexity, %ribo and per-cell intronic (nuclear) fraction; the program's median gene length and fraction of genes > 100 kb; Shapley share under QC terms | Reported only | Flag only; never used for deletion |

**Malignant→TME spillover (A1).** As the review pointed out, a random-gene null cannot absorb spillover, because random genes are not malignant-specific. The within-sample κ ratio can be inflated in the same way. A1 therefore works in four steps:
1. It scores TME cells as **pseudo-bulk** (log of each TME group's mean CPM), the scale on which ambient RNA mixes linearly.
2. It estimates each sample's ambient fraction **ρ̂** per TME cell type from malignant-specific markers. These are genes ≥ 8-fold higher in malignant cells than in every TME type, excluding the audited genes.
3. It **subtracts** the predicted spillover log2((1 − ρ̂) n_g + ρ̂ m_g + 1), where n_g is the cell type's native profile and m_g the sample's malignant profile.
4. It compares the adjusted partial R² with the **expression- plus specificity-matched** random-set null.

Why subtraction rather than a covariate: on the toy (before its technical axes were recalibrated, §8), entering the spill score as a covariate removed the genuine IFN signal in 3 of 3 seeds.

Why there is no fixed R² bar: the null itself is high and depends on the design. On the toy, null q95 ranged from about 0.2 to 0.8 across cohort designs (details under P6, §8).

Open problems:
- Subtraction under-corrects if real ambient RNA under-represents the markers (ρ̂ biased low), and the random-set null does not absorb the leftover leak, because random sets lack a real program's between-sample variation. P1 therefore re-runs A1 with every ρ̂ multiplied by 1.5 and counts an MP only if it stays above that run's null. This covers ρ̂ biased low by up to about a third. On a 16-sample toy with ρ̂ × 0.7 and the technical axes at their first magnitudes, a malignant-only CNA program rose from 0.15 to 0.41 against q99 0.27, and fell to 0.11 against 0.23 in the re-run; stress and IFN stayed above both. On the recalibrated toy the CNA leak is smaller (unadjusted partial R² 0.07–0.53 over seeds 0–4), and at ρ̂ × 0.7 or × 0.5 the CNA program stays below q99 at every one of these seeds, so this toy no longer exercises the re-run.
- With fewer than 5 malignant-specific markers, the adjustment is skipped with a warning. This happens when ambient RNA exceeds about 1/8 in a TME type, because ambient RNA caps a marker's fold change at about 1/ρ, or when a TME type shares the malignant lineage genes (e.g. normal epithelium). Such MPs cannot count through P1's criterion (b).
- The κ ratio (i) is not yet spillover-adjusted. Until it is, it is used only toward "excluded", where inflation is conservative, and never as sufficient for "supported". In P1 it is informational only.

**Capture and nulls (A4).** The choice of control genes matters:
- With control genes matched on expression pooled across protocols, protocol-specific capture (Flex probe coverage, long genes in nuclei) can masquerade as loss. On the small test cohort of the toy simulator (`tests/test_protocol.py`: 4 paired tumours, 3 replicate tumours), the 12 worst-captured cycle genes lost intensity beyond the replicate null in 3 of 4 paired tumours with pooled controls, 0 of 4 with protocol-specific controls and 0 of 4 by κ.
- Protocol-specific controls match every gene on its own level within a protocol. Any shift common to all paired tumours is therefore absorbed, whether it comes from capture or from biology.
- So intensity and prevalence cannot see a cohort-wide loss; on the toy, stress intensity changed by about 0 (mean −0.03) while its κ fell 3.2- to 5.5-fold in every paired tumour. κ is the only A4 statistic for a cohort-wide loss.

P1 tests every paired tumour; there is no activity gate. P1-GO v1.2 tested a paired tumour only where the program's κ under the reference or under the test protocol exceeded max(κ_min, the q95 of matched random gene sets in the same samples), a simplified form of κ*_t (§1.2), and let a replicate pair into the null on the same condition. The idea was that where a program is not expressed in either piece, its κ is background and its "change" is noise or the floor. But the gate selected tumours on the noisy κ that was then tested: a tumour whose κ was high by chance passed and then tended to drop. In the calibration harness (`scripts/calibration/`, v1.2 reference) the effect was large. Even with every piece's κ drawn exactly like its random sets', up to 0.31 of the tested tumours were losses on synthetic laws (0.05–0.31 across the exchangeable laws on the random-set path). On a toy null with no protocol effect, 0.13–0.17 of the tested (gene set, tumour) pairs without the program were losses, and an unplanted random gene set met criterion (c). P1-GO v1.3 therefore tests every paired tumour with a defined change and bounds, so which tumours count no longer depends on the noise of the statistic tested. On the same harness, 0.018–0.054 of the tested tumours are losses on the synthetic laws, and 0.019–0.063 of the tested pairs without the program on the toy null, with no gene set meeting (c) (toy and synthetic results only). The price is that tumours without the program count in the denominator of (c). A program present in fewer than half of the paired tumours cannot meet (c), however strongly it is lost where present, and can reach nuisance support only through (a) or (b). The activity columns are still reported, for information.

The bounds have their own limits, so a κ loss in P1 must beat three of them, each a one-sided 5% bound under its own assumption. The 5% level per tumour then holds if any one assumption holds:
- **Replicate null.** It needs same-protocol region pairs, which public twins mostly lack. It is the lower of two one-sided t prediction bounds over the replicate changes of ≥3 distinct tumours, one on the linear κ changes (carried to the MP's κ level) and one on the log κ changes. Each gives 5% false losses per tumour by Monte Carlo when paired and replicate changes are exchangeable and normal, against 13% for the empirical 95th percentile of 3 changes used before v1.1. The lower of the two covers both a program near the floor whose replicate pieces carry it well above it and a varying program whose replicate pieces sit at the floor.
- **Random-set null.** It holds additive κ changes shared by all genes, carried to the MP's κ level. The bound is the floor(0.05 (n + 1))-th smallest of the n random sets' changes (the 5th of 100), a distribution-free 5% bound for a change exchangeable with theirs; v1.2's interpolated 5th percentile was not such a bound. For the informational intensity statistic, which is compared with these bounds directly, it does not cover piece-to-piece variation.
- **Self bound.** The random-set null does not hold the MP's own sampling noise. So a loss must also exceed 1.645 × the bootstrap SE of the MP's own change, from 50 resamples of each sample's cells.
- **Equal cells and depth.** κ's estimate falls with the number of cells and with the counts per cell, and that fall scales with the set's κ level, so no additive null holds it. Every sample of a tumour is therefore scored on the same number of cells (its smallest sample, at most 2000), with counts binomially thinned to the tumour's lowest median count per cell, and the panel is taken on those counts. Without equal depth, toy test pieces at a quarter of the reference depth made 0.12 of the absent pairs and 0.06 of the carriers losses and an unplanted program met (c) (review probe, code before the step); with it, the harness's depth variants pass (§8).
- **Margin.** The remaining gap is regional heterogeneity of a program that is present: two pieces that differ in its activity look like a protocol effect. Neither the random sets nor the bootstrap carry it, and only the log 1.5 margin (the A4 equivalence margin) protects. For a weak program (κ 0.05–0.2) that margin is a κ change of only 0.02–0.07. On the replicate path, the null contains regional heterogeneity, but its t bounds assume near-normal changes: for a program absent from 10–30% of the pieces, P(c) on 4 paired tumours reached 0.031–0.047 with 3 replicate tumours (synthetic; §8).

Toy run (`simulate_cohort()` defaults, seed 0; the planted sets; `protocol_shift_kappa` with its defaults, 100 random sets per tumour and MP; re-run on 2026-10-18). All 32 κ tests of the eight planted programs other than stress are tested; the toy protocol is not planted to lower any of them. Losses and gains under each set of bounds:
- the null bound alone: 5 losses and 7 gains;
- the null bound and the self bound: 3 losses and 3 gains;
- the null bound and the margin: 1 loss, CNA_7q in T03, a paired tumour that does not carry the CNA;
- all three, the registered rule: none.

The rates of false (c) support measured on the toy and on synthetic laws are in [05_pilot_P1.md](05_pilot_P1.md) (Caveats).

The selection-versus-induction prediction (§7, result 3) needs a contrast that can see cohort-wide prevalence shifts. That is open.

**Gene-length and abundance axes (A5).** Gene length / nuclear (intronic) fraction and transcript abundance / %ribo are added as nuisance axes because of real-data evidence from the user's breast cancer snRNA-seq (HER2) project. This evidence is user-reported and not reproduced here. Two pre-registered de novo MP searches there produced mostly technical axes:
- several "luminal" programs of long genes (median gene length 317–437 kb);
- an abundance axis of short cytoplasmic genes that tracked %ribo in 29 of 30 patients (rho 0.68–0.77) and that reference matching labelled "EMT".

These axes are flags, not offsets (§2.2). In P1, %ribo and, where raw reads allow, intronic fraction also enter the technical decoy arm as continuous covariates, and the trivial typer gets the matching features (§6.2).
- **%ribo in the decoys is a deliberate exception to §2.2.** %ribo is computed from ribosomal-protein counts and tracks translation, ribosome biogenesis and MYC activity. A technical arm that includes it keeps the co-expression of such programs, and an MP containing RP genes shares counts with %ribo and keeps them mechanically (κ_retained can exceed 1). If the decoys recruit such an MP and keep ≥0.5 of its κ, criterion (a) counts it as nuisance support, so a real translation, ribosome-biogenesis or MYC program can count toward GO. On the toy, %ribo also moves every other gene through a shared per-cell dilution. Under v1.2's plain ratio, a %ribo-only arm was credited with that background and kept 0.10–0.21 of the κ of cycle, hypoxia and CNA. Under v1.3's excess ratio it keeps 0.87 of the abundance set's excess κ and between −0.07 and 0.04 of cycle, hypoxia, EMT, IFN and CNA (seed 0, decoy seeds 1000 and 1001). The planned arm keeps %ribo because the abundance axis was a dominant artefact in the user's data (above); whether to keep it is decided before registration, and a recruited MP that contains RP genes or tracks translation is to be read with this caveat.
- **Intronic fraction.** 3CA study folders hold gene-level counts only, so the per-cell intronic fraction needs the raw reads re-processed (STARsolo Velocyto features, velocyto or DropletQC). Which studies can supply it is counted at P0 and fixed at registration (§8).
- **Gene length in the decoys.** It enters only through a per-cell covariate: the intronic fraction or `pct_long`, the share of a cell's counts on genes > 100 kb (median length and fraction > 100 kb are per-MP features and cannot enter a per-cell mean model). Like %ribo, `pct_long` is computed from counts that include the audited MP's own genes and can inflate its κ_retained; it is available in `mpaudit` but not in the registered arm by default. Where neither covariate is in the arm, long-gene MPs are screened only by the typer's length features, and a decoy "survival" of such an MP is not evidence against a nuisance origin.

Shapley attribution runs over the 2^4 nuisance subsets {ambient, TME components, allele dosage, QC}, with bootstrap over tumours. The attribution is therefore **independent of the order** in which nuisance terms are added.

**Type labels** (each shown with its evidence vector):
- **B1 Certified intrinsic biology.** A1 and A2 excluded; A3 dosage excluded or not testable; A4 invariant on ≥1 bridged axis.
- **B2 TME-shared biology.** A1 supported, A2 excluded, A4 invariant (e.g. IFN response, hypoxia).
- **G Genetic.** Dosage-driven or clone-linked; labelled, never deleted.
- **P Protocol-modulated biology.** Invariant loadings with selection or induction shifts.
- **T Technical.** Contamination/ambient (A2 supported), ex vivo (A4 loss together with A1 supported), or low quality.
- **C Cohort-restricted.** Certified in a single cohort. Reported, not called an artefact. This protects lineage programs from cancer types covered by only one study.
- **U Indeterminate**, with the reason: underpowered; no protocol bridge; no allele data.

Each axis decision controls its own error at α (a TOST or one-sided test). The composite label has **no** formal error rate. Its accuracy is measured empirically on independent simulations and on labelled artefact sets. This is stated plainly.

### 2.6 Diagnostics (reported, not used for selection)

- **Decoy audit.** We generate non-parametric decoys by stratified conditional residual permutation: gene-wise permutation of NB randomised-quantile residuals, which keeps each gene's real marginal misfit. Depth and continuous covariates enter the decoy's mean model, and categorical covariates define permutation strata.
  - Decoys are built **per arm**, so that a recruitment can be typed:
    - a *technical* arm keeps depth, doublet score, ambient/contamination estimate, %ribo and (where raw reads allow) intronic fraction, jointly, so a technical hit is not attributed to one covariate (one arm per covariate would attribute it, at 5 pipeline passes per arm and with more chances to hit, so the arm list is fixed at registration);
    - a *dosage* arm keeps CNA clone.

    Recruitment by the clone arm types an MP as "dosage", not "nuisance".
  - The decoys are run label-blind through pipelines that do not model nuisance: Gavish code, GeneNMF, and the TESSERA engine without certification.
    - In P1 they run through the Python port of the Gavish pipeline (`mpaudit.gavish`), not the original R code.
    - R-versus-port parity is untested (§8). Until it passes on real samples, recruitment rates describe the port.
  - The output is a per-MP **decoy recruitment rate** (hits among the decoys of each arm; P1 registers exactly 5 per arm), together with **κ_retained**, the share of the MP's within-sample co-expression that the decoys keep. In P1 that is the excess over 50 matched random gene sets, per sample (κ_decoy − bg_decoy)/(κ_real − bg_real), over the samples whose excess is above 3 SD of the sets' κ. Every gene set of a sample co-varies a little through shared nuisance, and a decoy keeps that background, so a plain ratio would credit the decoy with it. Together they measure how much nuisance can manufacture.
  - It is valid as an audit because the audited pipelines do not condition on the nuisance model.
  - It carries no FDR claim.
  - Caveat: the covariates (doublet score, contamination, %ribo, intronic fraction, inferred clones) are themselves computed from the same counts or reads. A real program correlated with one of them survives the decoy and can be "recruited". Toy example (`mpaudit.sim`, default cohort, seed 0; decoy seeds 1000 and 1001; P1-GO v1.3's κ_retained): a decoy on the planted doublet score, coupled to cycling (within-sample Spearman 0.61), keeps 0.42 of cycle's excess κ; a decoy on the toy's estimated doublet score keeps 0.00. On this cohort that estimate barely tracks cycling (Spearman 0.26). On the null cohort (`--clean`, no TME ambient RNA) it tracks cycling at 0.51 and keeps 0.26. So the low default value does not show that estimated covariates are safe on real data. Two safeguards limit the damage. (1) Clone recruitment is typed as dosage and never counted, because decoys are built per arm. (2) Criterion (a) needs κ_retained ≥ 0.5 as well as hits, which screens out weaker coupling such as the 0.42 case. With the planted score in the four-covariate technical arm, cycle keeps 0.497, 0.532 and 0.530 on the default toy at seeds 0–2 (at the guard: just below it at seed 0) and 0.64–0.67 on the null toy (above it). Neither safeguard helps a real program tightly tied to a covariate's estimator (cell size and the doublet score; a state an expression-based clone caller splits on; translation and %ribo). Such a program keeps most of its κ, passes (a), and cannot be told from nuisance by the decoy audit. Doublet and contamination share the technical arm, so per-arm reporting does not separate them either. Only covariates measured independently of expression (genotype demultiplexing, cell hashing, protocol twins) settle such cases (see [05_pilot_P1.md](05_pilot_P1.md), Caveats). Preview, on code before v1.3 (to be replaced by the demo rerun): in the null toy cohort (`--clean`), the technical arm recruited the real cycle program in 5 of 5 decoys, and only its κ_retained of 0.35 under v1.2's plain ratio (< 0.5) kept it from counting. Under v1.3's excess ratio the arm keeps 0.35 of cycle there at seed 0 (two decoys; 0.07 and 0.08 at seeds 1 and 2), so at the demo's seed the margin to the 0.5 guard stays about 0.15. How often real biology is flagged on the toy (the false-flag rate) will come from that rerun.
- **Validity–stability plane:** certification strength on one axis, stability index on the other. A "stable artefact" shows up as its own quadrant.
- **Identifiability:** anchor genes and the fraction of near-pure cells per MP. The pooling monotonicity of the sufficiently-scattered condition is stated as a *remark*, not a theorem, together with empirical breakdown maps over overlap and private-cone geometry.

### 2.7 Optional between-tumour module

Tumour baselines b_t are decomposed onto W by KL-NNLS. Claims about tumour-wide or subtype programs are restricted to:
- tumours with ≥2 libraries or regions; and
- programs that replicate across cohorts.

Reason: where tumour = library, tumour-wide biology and library batch **cannot be separated**. This is stated explicitly, and the GBCD and Nomura/Tirosh multilayer framing is credited.

### 2.8 How the constants are set

**User-facing constants:**
- α (0.1);
- κ_min (0.05);
- the quantile of the random-set null band (0.95);
- r_min (5) and the number of cohorts (2);
- the equivalence margins in §2.5;
- the gene minimum effect (log 1.2);
- the match tolerance (cosine 0.7).

All are interpretable, frozen and pre-registered on OSF **before** the atlas run. The atlas run tests them on the held-out confirmation half of the studies, which nothing before pre-registration 2 touches (§8).

**Engine settings** (K_max, J, R, λ_x, n0, fold fractions) are chosen by held-out deviance or fixed defaults.

We publish a **full sensitivity table**: MP count, matched identities and type labels across every constant. This is the analysis the concepts asked of Gavish; we accept the referees' point that thresholds are never fully removed.

---

## 3. How each failure of the standard pipeline is fixed

| Step | Failure (landscape §1.2) | TESSERA fix | Evidence it works (pre-registered) |
|---|---|---|---|
| S0 malignant calling | Contamination and doublet MPs (Macro/T/Endo removed by hand in 3CA v2) | Same-sample TME profiles as fixed components with unconstrained usage (covers misclassified whole cells); axis A2 needs both a loading match and a cell-level doublet/ambient association | 3CA v2 contamination MPs typed T with no manual input; zero mimicry programs (pEMT, MES-like, MHC-II) typed T in planted-into-real tests |
| S1 top 7,000 by mean; centring; zeroing negatives | Ribo/mito LQ MPs; tumour-wide programs invisible; no likelihood | Deviance-selected universe plus TFs; count likelihood with depth offsets; tumour baseline b_t instead of centring; latent-correlation estimand independent of depth | LQ MPs fall into A5/T; membership rate flat against depth and cell number (|slope| ≤ 0.05 per decade above 200 cells) |
| S2 rank grid 4–9 | Merged or split programs | Joint dictionary with ARD; K chosen by molecular cross-validation; any engine can be plugged in | Subset-shared recall ≥0.8 at log2FC ≥1 on the independent simulator |
| S2 best of 10 restarts | Non-identifiability; arbitrary optimum | Deterministic initialisation; stability over 5 independent cell splits × seeds; anchor and near-pure diagnostics | Median matched cosine ≥0.9 across splits for certified MPs |
| S2 per-tumour fits | No borrowing of strength | Joint discovery with tumour balancing; empirical-Bayes prevalence across tumours | Rare (3–10 tumour) program recall reported separately |
| S3 top-50 truncation | No gene significance | Gene effects on a disjoint cell fold, minimum-effect null, three-level random effects, core vs context tiers | Empirical gene FDR ≤ nominal + 0.02; core-list length flat against cohort size |
| S4 35/10/10 | Biased toward dominant axes; random-gene null | Cross-rank stability dropped as a filter; effect-size membership plus replicated recurrence (r_L and held-out cohort groups) | Gavish code's false-MP rate on the independent simulator reported side by side |
| S5 greedy clustering | Order dependence; MP count a by-product of thresholds | MPs are rows of a joint dictionary (or engine programs); existence is a test with replication | MP-count coefficient of variation ≤10% under 50% study subsampling, compared with Gavish |
| S6 manual curation; CNA checked afterwards | Not reproducible | Five-axis equivalence typing with an "indeterminate" outcome; allele-anchored dosage | Typing AUROC on ≥40 labelled artefacts beats a trivial covariate typer by ≥0.05 |
| S7 scores without uncertainty | Not comparable across tumours | Projection with Fisher-information SEs, per-sample detection q-values, a novelty test and panel retention ι_k | Coverage of per-cell SEs on planted data within 0.85–0.95 |

---

## 4. Differentiation: what referees will ask, and our answer

| Competitor | Referee question | Answer |
|---|---|---|
| **GBCD** (Nat Genet 2025) | "GBCD already separates shared from patient-specific programs." | Agreed: we do not claim that split. GBCD has no existence or recurrence test, no tumour membership, no nuisance typing and no projector, and runs per cancer type. It is a TESSERA engine. We certify its output and benchmark GBCD *with* environment covariates. Covariate-moderated EBMF (cEBMF; Denault…Stephens, arXiv 2505.11639, NeurIPS 2025) lets side information moderate the factor priors of an EBMF/NMF fit, so known nuisance covariates could enter a GBCD-like model directly. It does not test cross-tumour recurrence or type programs by provenance; we treat it as a candidate engine. |
| **ITHmapper** (bioRxiv 2026) | "There is already a pan-cancer ITH package with a mapper." | Its states are discrete; sample-as-batch can erase tumour-wide signal; it has no uncertainty or provenance. We wrap scVI+Hotspot as an engine, certify its modules, and compare head-to-head on the matched cohort. |
| **TCAT/starCAT** (Nat Methods 2025) | "A reference plus projector is TCAT." | Agreed, the reference-plus-projection idea is not new. TCAT built a reproducible catalogue of 46 T-cell programs. The abstract reportedly presents starCAT as generalising the framework to other cell types and tissues; this is near-verbatim wording from search results, to be checked on the article before quoting. Our contribution is calibrated recurrence and provenance typing. TCAT's reproducibility is cross-dataset recurrence, which systematic artefacts pass, and its README cautions that the references are not established for 10x Flex. Our projector adds SEs, novelty detection, panel retention and protocol typing. |
| **TCCA** (González-Bermejo et al., bioRxiv 2026) | "A newer pan-cancer MP atlas exists." | TCCA has 36 public studies, about 540 patients plus 183 cell-line samples, and SCEVAN malignant calls. It derives 43 MPs by Gavish-style clustering of per-sample NMF (RcppML, single seed) followed by a hand-curated whitelist. The criteria (mito/ribo-rich, single-study, < 10 programs) are stated only in a code comment. It is a second example of manual MP curation and a standard-pipeline baseline. We use its patient tumours, minus cell lines and minus studies shared with 3CA, for replication; the overlap is counted at P0. |
| **Grabski multi-study NMF / Tetris / APAFA / SUFA** | "Subset sharing is solved." | For mutational signatures and bulk data, with MCMC. If ported to single-cell data, these models become engines; our claims do not depend on a correct posterior. A Tetris run on a 30-tumour subsample is held in reserve. |
| **scDEF, consensus scHPF** | "There are already count-native hierarchical models." | They are engines without existence, membership or provenance inference. Benchmarked. |
| **cNMF / GeneNMF / Gavish code** | "Why not just fix the thresholds?" | They are wrapped as engines; the decoy audit and TESSERA certification are applied to their outputs. |
| **Spectra / DRVI** | "Prior-guided or disentangled factors?" | Neither has a recurrence or provenance inference layer. DRVI is benchmarked as the representative deep model. |
| **CanSig benchmark** | "A benchmark already exists." | It has no planted confounders, null cohorts, calibration or typing metrics. We reuse its 12 datasets. |
| **Clipper / ClusterDE / sc-SHC / recall** (referee-reported) | "Target-versus-null contrast and model-based nulls are known." | Correct, and we cite them. We do *not* use synthetic nulls for inference. Decoys are an audit diagnostic only. The inferential contribution is cross-fitted, effect-sized and specific tumour membership, plus real-negative-control provenance. |
| **FDR via data splitting after clustering** (Wang, Lin & Zhao, arXiv 2410.06451) | "Error control after clustering exists." | It splits the data, tests after clustering on each half and combines the results by a mirror statistic, with asymptotic FDR control in Gaussian settings. This is the closest statistical relative of our cross-fitting. Its multiple-split mirror statistic is a candidate for combining our split replicates, but it does not address cross-tumour recurrence or provenance. |
| **Gene Program Dictionary** (Xu…Ma, bioRxiv 2025.10.29.685322) | "A program catalogue with FDR exists." | It reports 1,534 consensus programs from 466 scRNA-seq datasets with permutation-based FDR (details from search snippets). It does not type tumour MPs by provenance against negative controls. It serves as adjacent work and as a reference to project onto. |
| **Count splitting** (Neufeld et al.) / **thinning for Poisson factor models** (Wang et al., JASA, referee-reported) | "Valid inference after any latent-variable estimator already exists." | We use cell splitting for inference, which needs no dispersion assumptions, and thinning only for model selection. No new thinning theory is claimed. |
| **Module preservation (NetRep) / CS-CORE** (referee-reported) | "Audit mode is module preservation; κ is CS-CORE." | Both are baselines. We adopt the CS-CORE-style estimator and add specificity conditioning, within-tumour null bands, cross-tumour recurrence and provenance. Audit mode must beat NetRep Zsummary at flagging labelled artefacts. |
| **ISOMERIC** (bioRxiv 2025), **TreeAlign, CLCC** (referee-reported), **scPFA** (bioRxiv 2026) | "CN-dependent vs CN-independent expression is done." | ISOMERIC uses an adversarial autoencoder to separate copy-number-dependent from copy-number-independent variation in 189,027 tumour cells from 10 cohorts. It notes that CN dependency confounds the search for shared states. Its goal, like TreeAlign's and CLCC's, is a decomposition within tumours. Ours is a calibrated typing axis for recurrent cross-tumour programs, used only with allele anchoring. ISOMERIC is benchmarked where its code runs. scPFA (lineage-informed factor analysis of heritable programs) overlaps the clone/dosage axis; its abstract was not retrieved, so it must be read before we position against it. |
| **SIMVI, DECIPHER** (Nat Commun 2025) | "Intrinsic versus context-induced states are already disentangled." | Both work within spatial datasets: SIMVI separates intrinsic from spatially induced latent variables, and DECIPHER learns separate molecular and spatial-context embeddings at atlas scale. A1 asks the related question for recurrent programs in dissociated scRNA-seq, using co-resident TME cells. We will use both methods on the Xenium/CosMx rung (§6.4). |
| **inVAE / sysVI** (referee-reported) | "Invariant representation learning exists." | Those learn representations. We produce non-negative gene programs with paired, equivalence-tested protocol invariance on real twins. |
| **Nomura/Tirosh multilayer GBM** (Nat Genet 2025, referee-verified) | "Baseline vs within-tumour programs is known." | Credited. Our between-tumour module is optional, with claims restricted to multi-library tumours. |

---

## 5. Package architecture

**Name.** Working name TESSERA (`tessera-mp`; R: `tesseraR`). Alternatives: `mpcert`, `CREDO`. Run a PyPI, CRAN and Bioconductor collision search before any public repository is created.

**Ecosystem decision: Python/scverse core with an R bridge.** Reasons:
- GPU (PyTorch);
- AnnData, MuData and SpatialData for Xenium and CosMx;
- the scvi-tools ecosystem (DRVI, scVI adapters);
- TileDB-SOMA and CELLxGENE for streaming at atlas scale.

R remains essential, because the Gavish code, GeneNMF, GBCD and most of the user base live there:
- `tesseraR` (reticulate plus zellkonverter) provides `CertifyMetaprograms(seu)`, `AuditGeneSets(seu, gmt)` and `ProjectMP(seu, ref)`, and returns a Seurat assay or DimReduc;
- R engines run through rpy2 inside versioned containers.

**Modules:**
- `tessera.io`: loaders for 3CA, TCCA, CanSig, CELLxGENE and 10x Flex/Xenium; a check for integer UMIs.
- `tessera.schema`: EnvSchema, harmonise(), and a `design_graph()` report of which axes are testable.
- `tessera.pp`: imports of ambient (CellBender/SoupX/DecontX), doublet (scDblFinder) and CNA (Numbat; infercnvpy/SCEVAN/CopyKAT flagged "expression-only"); deviance-based gene selection.
- `tessera.split`: cell folds (𝒟/ℳ/𝒢), split replicates, and binomial thinning for molecular cross-validation.
- `tessera.engine`: `TesseraNMF` (PyTorch, minibatched MU, KL-NNLS) plus adapters (`GavishNMF`, `CNMF`, `GeneNMF`, `SCHPF`, `GBCD`, `DRVI`, `ScviHotspot`, `Custom`).
- `tessera.certify`: the κ estimator, matched random-set null bands, interval tests, empirical-Bayes prevalence, partial conjunction and grouped leave-cohort-out replication.
- `tessera.genes`: leave-gene-out NB GLMs, three-level REML/HKSJ, minimum-effect BH, GMT export with tiers.
- `tessera.type`: axes A1–A5, Shapley attribution, TOST decisions, labels plus evidence cards.
- `tessera.audit`: certifies external catalogues (3CA, Gavish 2023, Kinker, Barkley) on independent cohorts.
- `tessera.decoy`: the stratified conditional residual permutation (SCRP) decoy audit and recruitment rates.
- `tessera.ref` / `tessera.project`: a versioned reference; projection with Fisher SEs, novelty test, panel retention ι_k and per-sample q-values.
- `tessera.sim`: two simulators (§6).
- `tessera.bench`: baseline registry with Apptainer/Docker images and metrics.
- `tessera.pl`, `tessera.report`: MP cards, validity–stability plane, karyograms, forest plots, HTML run report.

**Key API:**
```python
import tessera as ts
ts.pp.setup(adata, layer="counts", tumour="tumour", sample="sample",
            study="study", lineage="cancer_type", malignant="is_malignant")
ts.schema.harmonise(adata, meta)                 # protocol/pre-analytic attributes
ts.pp.add_ambient(adata, source="cellbender"); ts.pp.add_cna(adata, source="numbat")
res = ts.run(adata, engine=ts.engine.TesseraNMF(K_max=200),
             folds=(0.4, 0.3, 0.3), alpha=0.1, kappa_min=0.05, r_min=5)
res.mps; res.membership; res.genes(tier="core"); res.types; res.pl.plane()
ts.audit(adata_new, gene_sets=ts.datasets.gavish2023())   # adoption path
ref = ts.ref.build(res, version="1.0")
ts.project(query, ref, panel=xenium5k_genes)              # SEs, novelty, iota_k
```

**GPU backend.**
- PyTorch with sparse CSR minibatches; dense denominators are batched GEMMs.
- Cost estimate: about 3e13 flops per epoch at 4M cells × 8k genes × K = 200, which is minutes per epoch on an H100. I/O dominates.
- Shards are streamed per tumour from zarr or TileDB-SOMA.
- κ estimation and null bands run as batched GPU matrix products.

**Reproducibility, testing and documentation.**
- Every run writes a config hash, seeds, software versions and data hashes.
- Snakemake pipelines `tessera-atlas` and `tessera-bench` regenerate every figure from a separate repository.
- CI runs unit tests plus **calibration tests on every release**: type-I error of membership and gene tests on small independent-generator nulls, and coverage of projection SEs.
- Distribution on pip and bioconda; Read the Docs with four tutorials (audit my MPs, discover in my cohort, project my Flex/Xenium data, run the benchmark).
- The reference and benchmark are versioned on Zenodo with DOIs.

---

## 6. Benchmark and validation plan

### 6.1 Simulators: independent of the model by design

- **Sim-I (primary, independent generator).**
  - Built on scDesign3 and splatPop fitted per real tumour (CanSig's 12 datasets plus about 300 3CA UMI tumours).
  - Confounders come from a **separate injector** whose functional forms are *not* those of TESSERA:
    - ambient that is nonlinear and depends on cell size;
    - doublets added in proportion to RNA content;
    - whole TME cells misclassified as malignant;
    - focal and arm CNAs with tumour-varying dosage compensation, with CNA calls **re-inferred** from the simulated expression;
    - logistic, cell-level stress activation that varies by cohort and *crosses* cohorts;
    - tumour-specific dispersion;
    - non-NB marginals (lognormal-Poisson, bimodal).
- **Sim-R (planted into real data).**
  - Programs are added to real malignant cells by count addition, keeping all real structure.
  - **Artefacts are planted into dissociation-free backgrounds**: measured warm-dissociation signatures injected into Flex or single-nucleus malignant cells, so the ground truth for "ex vivo" is known.
- **Planted design per replicate:**
  - 6 pan-shared, 12 subset-shared and 6 rare programs;
  - 4 pairs of overlapping programs sharing 10–25% of genes (to test specificity);
  - 3 between-tumour programs;
  - 2–3 private programs per tumour (half CNA-linked);
  - effect sizes 1–20% of UMIs; active-cell fractions 3–50%.
- **Null cohorts** from Sim-I (nuisance only).
- **A fractional, space-filling design of about 1,000 cohorts**, not a full factorial (whose cost the referees showed to be infeasible).
- **An adversarial scenario set designed by a blinded third party** (a collaborating statistician who does not see the code).

### 6.2 Baselines and fairness protocol

Baselines:
- Gavish/3CA code (exact);
- GeneNMF;
- cNMF per tumour plus Gavish-style clustering;
- TCAT-style Harmony+cNMF;
- GBCD per cancer type, plain and with environment covariates (subsampled to 20k cells);
- ITHmapper (scVI+Hotspot);
- consensus scHPF or scDEF;
- LIGER online iNMF;
- DRVI;
- Grabski multi-study NMF with tumour as study, on a subsample.

Simple baselines, which every component must beat:
- (a) published MPs plus NetRep module preservation on independent cohorts;
- (b) per-tumour NMF plus count splitting plus BH;
- (c) a trivial typer: MP-score correlation with %mito, %ribo, intronic fraction and doublet score; median gene length and fraction of genes > 100 kb; TME-marker cosine; arm enrichment.

Fairness:
- every method runs with author defaults **and** a documented grid under an equal compute budget;
- tuning on a separate tuning scenario set and evaluation on held-out scenarios;
- every baseline also gets a "same preprocessing" arm (CellBender plus doublet removal plus TME-cell exclusion), so that gains from the model are separated from gains from preprocessing;
- "TESSERA-wrapped Gavish/cNMF" separates the effect of the certification layer from the effect of the engine.

### 6.3 Metrics scorecard

| Axis | Metric | Pre-registered target (TESSERA) |
|---|---|---|
| Recovery | Hungarian-matched cosine; F1 at cosine ≥0.8; gene AUPRC | Subset-shared F1 ≥0.85 at ≥5% of UMIs |
| Null behaviour | MPs certified per null cohort (Sim-I) | Median 0; P(≥1) ≤ 0.1 |
| False-MP rate | Certified MPs aligned to injected nuisance (cosine >0.6) and typed B1/B2 | ≤5% |
| Specificity | False membership of non-expressing tumours that carry an overlapping program, for the conditional κ estimand with an *estimated* dictionary | ≤5%. For reference, from a toy with a perfect dictionary (single tumour, 300 genes): ≤1% up to 30/40 shared genes with conditioning; without conditioning, 0% until the overlap κ crosses the 0.05 floor, then about 100% (`specificity_fix_check.py`). The ref1 figures (72% at 5/40, 100% at 10/40) reject a null that is false. They measure a different estimand and are not compared with this target |
| Depth invariance | Slope of membership rate against log10(cells) for fixed planted effects | \|slope\| ≤ 0.05 per decade above 200 cells |
| Calibration | Type I of tumour and gene tests under non-NB marginals; empirical gene FDR | Within [0.03, 0.07] at α = 0.05; FDR ≤ nominal + 0.02 |
| Typing | Macro-F1 on Sim-R artefacts; "indeterminate" rate as power falls | Macro-F1 ≥0.8 among determinate calls; false B1 ≤5% |
| Leakage | Private programs certified as MPs | ≤5% |
| Stability | Matched cosine across splits; MP-count CV under 50% study subsampling | ≥0.9; ≤10% |
| Compute | Runtime and memory, 10^5–5×10^6 cells | 3CA-scale fit ≤ 24 GPU-hours |

Hypotheses to measure, not assume:
- the Gavish code certifies ≥3 nuisance-aligned MPs per Sim-I cohort and leaks ≥20% of private programs;
- GBCD's shared programs include ≥1 planted cross-cohort artefact.

### 6.4 Evidence ladder: each rung with a primary endpoint compared head-to-head with Gavish/3CA, GBCD and ITHmapper MPs

1. **Semi-synthetic (Sim-I and Sim-R).** §6.3.
2. **Perturbation (Tahoe-100M; Replogle as backup).**
   - QC-matched pseudobulks, with cytotoxicity and viability modelled.
   - The primary endpoint is *not* the canonical pairs, which any gene set passes. It is that B1/B2 MPs respond specifically to at least one mechanism-matched drug more often than T-typed programs (Fisher test).
   - Plus five pre-registered **non-canonical** drug–MP pairs whose outcome is uncertain.
3. **Cross-protocol, the decisive rung.**
   - New matched cohort (§7) plus public twins:
     - Janesick (1 block);
     - Trinks lung: patient-matched fresh, cryopreserved and FFPE; n = 3 per the referees, and the patient count and Flex use are not re-checked;
     - the ccRCC fresh scRNA-seq / frozen snMultiome / FFPE snFlex preprint (Xirenayi et al., bioRxiv 2025). Whether its specimens are patient-matched is unverified; P0 confirms it before the preprint is used as a twin;
     - the multi-platform COAD/HCC/OV and breast/NSCLC Xenium sets (referee-reported).

     Two further candidates, both unchecked, were found in October 2026: a GEM-X Flex performance preprint (bioRxiv 10.64898/2026.03.06.710062) and a Frontiers 2026 brief report on matched FFPE and frozen samples.
   - **Oliveira et al. (Nat Genet 2025) is not a fresh-versus-fixed twin.** Its Chromium Flex (8 FFPE blocks of tumour and normal tissue), Visium HD (5 patients) and Xenium (4 tumours) data all come from the same FFPE blocks, and we found no fresh-dissociated malignant cells in it. It is used only as a fixed-dissociated versus in situ panel comparison.
   - Most public twins lack same-protocol replicate pieces, so their A4 null is the random-set null, which misses regional heterogeneity; the log 1.5 minimum effect then governs (§2.5).
   - Primary endpoint: AUROC of the "ex vivo/protocol-sensitive" calls for predicting loss in fixed or in situ data from the same tumour, against a multi-region null.
   - Secondary endpoint: concordance of B1-typed MP usage across protocols.
4. **Spatial.**
   - Within-malignant-cell co-expression on Xenium 5K / CosMx 6K, with segmentation-spillover controls.
   - Pre-registered niches (hypoxia near necrosis; IFN near T cells), controlling for clone domains (CalicoST).
5. **Clone-orthogonal.**
   - Allele-aware Numbat on ≥6–8 BAM-accessible datasets (≥60 tumours) plus clonealign scDNA sets.
   - Target: agreement of A3 calls between allele-based and expression-based clones, Cohen's κ ≥ 0.6. The power calculation is pre-specified.
6. **Held-out cohort.**
   - Discover on 3CA; replicate on TCCA's patient tumours after removing its 183 cell-line samples and every study shared with 3CA. The overlap is not stated in TCCA's README or abstract, so it is counted from GEO accessions at P0, together with the number of independent patient studies that remain.
   - Target: ≥80% of B1 MPs re-certified, against ≤40% of T programs.
7. **Clinical.**
   - TCGA associations of B1 MPs and of T-typed MSigDB GAVISH_3CA sets, each against 1,000 expression-matched random signatures (Venet 2011).
   - Hypothesis: T-typed sets are no better than random.

**Labelled artefacts.** At least 40 positives:
- the 11 v2 `inds_rm` MPs, minus the mis-annotation;
- Gavish 2023's removed MPs (low quality, doublet, single study);
- TCAT's technical programs;
- protocol-named MPs in the 2026 catalogue;
- Sim-R artefacts.

AUROC is reported with Hanley–McNeil CIs, against the trivial typer, with typing rules frozen by a team member blind to the labels.

---

## 7. Flagship biological application

**Title: "Which pan-cancer metaprograms are malignant biology? A calibrated provenance atlas."** Everything below is a hypothesis.

**Data:**
- every 3CA v1/v2 tumour with raw UMIs and TME cells, plus CELLxGENE tumour sets. P0 sets the final scope; the realistic target is 1,000–2,500 tumours, ≥8 cancer types and 2–5×10^6 malignant cells. This is lower than the concepts' 10^7, as the referees argued.
- TCCA patient tumours (no cell lines, no studies shared with 3CA), held out for replication.
- The 3CA studies are split at month 0 into an audit half (P1, P6, P7, the parity samples, all method development before the M14 registration and the month-10 audit paper) and a confirmation half. Only the confirmation half is used for the pre-registered TESSERA analyses (§8).

**Results:**

1. **A provenance card for every Gavish-41 and 3CA v2 malignant MP**, certified on studies *not* used to derive each list (for the lists derived from 3CA studies, this means TCCA, CELLxGENE sets, the new cohort and confirmation-half studies outside each list's derivation set). It reports r_L, lineage prevalence, type, Shapley shares, stability and core genes. Expectations:
   - cycle, IFN/MHC-II and hypoxia are B1 or B2;
   - LQ and contamination MPs are T;
   - **the P1 go/no-go result:** under P1-GO v1.3 (§8), at least 3 retained, non-dosage MPs have nuisance support. The atlas then types them with the full evidence vector on the confirmation half.
2. **Per-tumour decomposition of "stress"**, the most contested MP. This is quantitative and non-obvious, but neither the existence of dissociation stress nor its protocol dependence is new:
   - van den Brink et al. (2017) showed dissociation-induced expression in a subpopulation of muscle stem cells.
   - O'Flanagan et al. (2019) compared cold-protease (6 °C) with collagenase (37 °C) dissociation in patient tumours, breast-cancer PDXs and cell lines (155,165 cells), and defined a 512-gene collagenase-associated stress core including FOS and JUN.
   - Neuschulz et al. (2023) showed directly, by 4sU-labelling transcripts made during dissociation of mouse tissues, that dissociation adds a shared core stress signature plus cell-type-specific responses.
   - Fixation before dissociation (FixNCut; Jiménez-Gracia et al. 2024) reduces stress-related artefacts.
   - Baron et al. (2020), by contrast, argue from zebrafish melanoma, spatial transcriptomics and human samples that a stress-like cancer-cell state also exists *in vivo*. They acknowledge the dissociation concern but do not establish dissociation stress.

   The stress MP is therefore plausibly a **mixture** of ex vivo induction and an in vivo state. To our knowledge, nobody has measured its per-tumour ex vivo share. Our novelty is that per-tumour attribution, not the existence of the response. Possible pre-emption: a 2026 preprint comparing dissociated with intact tumour was seen only as a search snippet and has not been identified (P8). Our results are new in three ways:
   - (a) the ex vivo share of the stress MP's usage **in each tumour**, compared with the TME-cell sensor;
   - (b) whether a protocol-invariant HSP/UPR component is B1 (malignant-intrinsic) or B2 (TME-shared);
   - (c) which *other* canonical MPs carry measurable ex vivo induction. Candidate: hypoxia induced by warm ischaemia.
3. **Selection versus induction.** Prediction: pEMT/mesenchymal-like prevalence is *higher* in Flex and single-nucleus data than in fresh dissociation of the same tumour. That is selection against fragile cells, typed P-selection; a gene-list pipeline would read it as absent. Testing it needs a contrast that can see a cohort-wide prevalence shift. A4's prevalence, which is tumour-relative with protocol-specific controls, cannot (§2.5).
4. **Malignant-intrinsic versus TME-shared programs.** For each MP, the fraction of its tumour-level variance that co-varies with the same gene set in co-resident TME cells, after subtracting predicted malignant→TME spillover (§2.5). This separates niche-driven programs (IFN, hypoxia) from cell-intrinsic ones. Validated by spatial and CCLE recurrence.
5. **Genetic provenance where alleles exist.** A calibrated estimate of the fraction of recurrent ITH that is dosage-driven per cancer type (e.g. chr7/chr10 in GBM). This replaces the post hoc 24%, 8% and 39% figures (Gavish, Kinker).
6. **Missed or merged programs.** Subset-shared programs from small tumours and split hybrids (EMT variants), each required to replicate in TCCA.

**New matched cohort (wet-lab partner).**
- 24 fresh resections, 3 cancer types (breast, colorectal, lung).
- Each split into: warm 37 °C enzymatic; cold protease; fixed-then-dissociated Flex; snRNA-seq; Xenium 5K on a serial section.
- In 8 tumours, two same-protocol pieces from distinct regions (the regional null).
- Recorded ischaemia time, digestion time and DV200.
- Optional 4sU labelling during dissociation (following Neuschulz et al., MSB 2023) in a 6-tumour subset.
- Power: with n = 24 paired tumours and a pre-registered family of ≤10 MP × axis hypotheses (Bonferroni α = 0.005, 80% power), the detectable paired effect is d_z ≈ 0.75. With n = 12 it is d_z ≈ 1.1. This is my calculation.
- Optional perturbation: CRISPRi of 3–5 anchor genes of one newly certified B1 program in 2 lines, with a small Perturb-seq readout. Off the critical path.

---

## 8. Go/no-go pilots (months 0–6) and the publishable floor

The pilots were first squeezed into months 0–3. Their durations sum to **18 person-weeks** (P0–P6 and P8; P7's 6 weeks overlap), and four of them need components that do not exist yet. The realistic schedule has three groups:
- **Months 0–3** hold what can run now: P0, P1 (`mpaudit` as built), P6 and P8, which is 9 weeks, plus the R-versus-port parity check.
- **Months 3–6** hold P2 and P3 (5 weeks). They need κ cross-fitting and the Sim-I generator, neither of which is built.
- **P4 and P5** (4 weeks) wait for labels and allele data, and may slip to month 9.

| Pilot | Starts | Duration | Decision criterion |
|---|---|---|---|
| **P0 Data and design census** | Month 0. Metadata only until the OSF registration | 2 weeks | GO if ≥1,000 tumours with raw UMIs **and** TME cells across ≥8 types, and ≥15% of candidate MPs are used by tumours that bridge ≥2 protocol levels. Also count: empty droplets; BAM access (submit DAC applications in week 1); TCCA overlap, from GEO accessions; the cohort × lineage × protocol rank; per study, whether a per-cell intronic fraction can be obtained (raw reads or BAMs, or spliced/unspliced layers; 3CA matrices have none); and whether each half of the study split keeps P1's powered minimum (≥4 paired tumours, ≥10 samples with TME cells). A metadata count of paired tumours is only an upper bound on the per-MP count the rule uses (tumours where the MP's κ change, null bound and self bound are defined). The MP-usage count for the bridge share runs after registration, on the audit half. If the bridge share is <15%, A4 depends on the new cohort, so a wet-lab partner becomes mandatory. |
| **P1 Audit of the standard against real nuisance** (key gate) | Month 0–1, after registration. Runs now with `mpaudit` | 4 weeks, plus NMF compute | Audit of 3CA MPs on the audit half: decoy arms, A1 TME co-variation and A4 protocol κ in public twins; instability, positional flags outside the dosage rule and the other reported statistics listed below are informational. The decision follows **P1-GO v1.3** (below): GO, PIVOT or INCOMPLETE. **PIVOT**: re-centre on calibrated membership, gene tiers and projection, and change venue (Genome Biology / NAR GAB). **INCOMPLETE**: close the data gap that left every criterion unpowered and rerun under the same frozen rule. |
| **P2 Specificity and depth invariance of κ** | Months 3–6. Needs κ cross-fitting (`tessera.certify`) and Sim-I v0, neither built | 3 weeks | Overlap leakage ≤5% at 10–25% shared genes with an *estimated* dictionary (the toy check used a perfect one, §1.3); membership slope ≤0.05 per decade; power ≥0.8 for G2/M and hypoxia in tumours with ≥200 cells. Otherwise revise the null band or margins. |
| **P3 Cross-fitting calibration** | Months 3–6. Needs Sim-I nulls | 2 weeks | Type I in [0.03, 0.07] under non-NB marginals on Sim-I nulls; gene FDR ≤0.07 at nominal 0.05. |
| **P4 Typing beats trivial** | About month 6. Needs ≥40 labelled artefacts and the evidence-vector typer | 1 week | Evidence-vector AUROC exceeds the trivial typer by ≥0.05 on the expanded labelled set. Otherwise typing is demoted from the headline. Only 11 labels exist now (the 3CA v2 `inds_rm` MPs), and `mpaudit` has only the trivial typer. |
| **P5 Positional/CNA audit** | When BAM access arrives (DAC at month 0, then 3–6 months), so months 4–9 | 3 weeks (on BAM subset) | If <5% of MPs are dosage-affected, A3 becomes an annotation rather than a headline. Until then only `mpaudit`'s expression-only positional flag exists, and it cannot establish dosage. |
| **P6 TME sensor** | Months 0–3, on the audit half. Runs now (`mpaudit.tme`) | 2 weeks | Per-sample malignant stress usage against the TME IEG level, spillover-adjusted, with a positive slope. Stress must be **above the random-set null q99**; cycle and hypoxia must **not** be above it. Fixed partial R² bars are not used. On the toy, the null alone reached q95 0.46–0.78 and q99 0.63–0.86 (default cohort, seed 0, 500 sets, the nine planted sets), and q95 0.18–0.41 and q99 0.29–0.57 for the same sets in a fresh-only toy (16 samples of 100 malignant and 3 × 50 TME cells, 1,500 genes; seed 0, 200 sets). The old "≤0.05 for hypoxia" bar failed on the earlier toy demo (0.32). On the current toy (default cohort, seed 0), stress passes (0.98 vs q99 0.85) and cycle (0.13 vs 0.65) and hypoxia (0.11 vs 0.67) stay below. Passing on real data validates A1 specificity; the toy run only checks the code. |
| **P7 Engine prototype and scaling** | Month 3 onward; real data from the audit half only, plus simulations | 6 weeks (overlaps) | Matched cosine for subset-shared programs ≥ GBCD + 0.05; a 3CA-v1-scale fit in ≤24 GPU-hours (timed on audit-half data scaled to 3CA v1 size). Otherwise the default engine is TESSERA-wrapped cNMF, with the headline unchanged. |
| **P8 Scoop check** | Month 0 | 1 week, then every 2 weeks | Landscape §6.3 search terms plus: "cross-fitted gene program", "metaprogram provenance", "protocol-invariant program", "decoy metaprogram", and papers citing GBCD, TCAT, Gavish, O'Flanagan, Neuschulz and Marsh. Identify the 2026 dissociated-versus-intact tumour preprint seen only as a snippet, and read scPFA. NO-GO only if a preprint already types recurrent MPs against within-sample negative controls. A search in October 2026 (about 14 web searches, without bioRxiv full-text or PubMed access) found none; this is weak evidence. |

**P1-GO v1.3** (frozen 2026-10-18). It is implemented in `src/mpaudit/report.py`, which is the single source of truth (its `rule_text()` is printed into every run's `summary.md`); [05_pilot_P1.md](05_pilot_P1.md) states the same rule. Versions, none of them registered or run on real data:
- v1.0 (2026-10-16): the first frozen rule.
- v1.1 (2026-10-17): v1.0 plus three requirements added by the code review before any real data were opened: the ρ̂ × 1.5 condition in (b), the log 1.5 κ margin in (c), and the ≥3 distinct-tumour t prediction bound for the replicate null in (c).
- v1.2 (2026-10-17): v1.1 plus two changes. (i) The κ activity gate in (c): a paired tumour was tested, and a replicate pair entered the null, only where the MP's κ in at least one of the two pieces exceeded max(0.05, the q95 of matched random gene sets in the same samples). (ii) The registered settings outside the `decision` block are enforced, each at exactly its registered value, and so are the inputs without which an exclusion cannot run: a config that deviates is refused unless it sets `unregistered_variant: true`, and the verdict is then labelled an unregistered variant.
- v1.3 (2026-10-18): v1.2 with criteria (c) and (a) changed, because the gate selected tumours on the κ that it then tested and `kappa_retained` let in samples where the MP is absent (revision note, item 15; §2.5). (c) has no activity gate, and a loss must beat the null bound (the lower of the linear and log t bounds, or a random-set order statistic), a self bound from the bootstrap SE of the MP's own change, and the margin. Every sample of a tumour is scored on the same number of cells and at the same depth, with the panel on those counts, and no draw depends on the other MPs or on the cell order. (a)'s `kappa_retained` is the excess ratio over 50 matched random gene sets with a 3-SD entry. Every random draw differs from v1.2, so no v1.2 output carries over. v1.3 is the version to register.

- **Unit.** *Retained* MPs are audited MPs not labelled as curated artefacts (`typer.labels` = 1; e.g. the 11 MPs removed by hand in 3CA v2). Rediscovering curated artefacts is not news. Re-derived MPs inherit the label of their best-matching labelled MP (`typer.labelled_gmt`, Jaccard ≥0.3, `min_jaccard`). Both files are registered. Without labels no MP is excluded and the summary says so: with re-derived MPs (`mps: rederive`) that is a deviation from the registered rule, since rediscovered curated artefacts could then count; with MPs from a GMT the config check warns. Labels given by name for re-derived MPs are a deviation. The gene sets of the 11 removed MPs are absent from the published v2 list, so for a re-derivation run they must be rebuilt from `MP_analysis_new_batch.R`, exporting the MPs before `inds_rm` is applied. This is not yet done; until it is, a re-derivation run is an unregistered variant. (An audit of the published v2 list does not need them.)
- **Pre-declared exclusions.** These are written into the config's `decision` block and registered before any 3CA expression data are opened.
  - *Expected TME-shared biology.* An MP is excluded here if it is named on `expected_tme_shared` or overlaps a pre-declared gene set (`expected_tme_shared_gmt`, a GMT file, or `expected_tme_shared_sets`, inline in the config) by ≥10 genes. This disables criterion (b) only; such an MP can still count through (a) or (c). This is a deliberate departure from the 2026-10-16 review, which asked to drop listed MPs altogether. The list exists because TME co-variation (b) cannot tell niche-driven biology from leakage. Technical-decoy recruitment (a) and a κ loss under a dissociation-free protocol (c) do not rest on TME co-variation. An expected-TME-shared MP that the decoys reproduce, or that loses κ without dissociation (e.g. hypoxia induced by warm ischaemia, §7 result 2c), is exactly the nuisance P1 is meant to find; in the typing scheme it is T, not B2 (§2.5). Names work only when published MPs are audited; re-derived MPs are named by the run, so for them the list must be gene sets.
  - *Dosage.* An MP is typed dosage if its genes are arm-enriched (BH q < 0.05 after excluding known positional gene clusters), or if a dosage decoy arm (CNA clone) recruits it (≥2 hits and κ_retained ≥0.5, whatever the number of decoys). Dosage MPs are labelled and never counted, whatever their other support: they are genetic, not technical nuisance.
- **Nuisance support**, any of:
  - (a) *technical decoys*: in a technical decoy arm, ≥2 of its 5 decoys (`decoy.n`, registered at exactly 5; ≥5 are needed) recruit the MP (a decoy MP at Jaccard ≥0.3), and κ_retained ≥0.5. κ_retained is the share of the MP's co-expression in excess of matched random gene sets that the decoys keep:
    - per malignant sample with ≥20 cells, the ratio is (κ_decoy − bg_decoy)/(κ_real − bg_real), on the MP's genes detected in ≥5% of the sample's cells (at least 5). bg is the median κ of 50 random gene sets (`decoy.n_null_floor`) matched gene by gene on mean expression and detection, scored on the same counts;
    - a sample enters only if ≥20 sets have a finite κ and κ_real − bg_real > max(0.02, 3 SD of the sets' κ) (`decoy.min_kappa` 0.02);
    - κ_retained is the median over the entering samples, averaged over the arm's decoys. If no sample enters it is undefined, and the MP cannot meet (a) in that arm or be typed dosage by it;
  - (b) *TME co-variation (A1)*: a positive slope and a spillover-adjusted partial R² above the q99 of random gene sets matched on expression and malignant specificity, in ≥10 samples with TME cells, for an MP not on the expected-TME-shared list. The statistic must actually be adjusted (not the unadjusted fallback used when markers are too few), and the association must hold again (positive slope, partial R² above that run's q99) with every ambient fraction ρ̂ multiplied by 1.5. There is no fixed R² threshold;
  - (c) *protocol (A4)*: a one-sided κ loss under the dissociation-free protocol in ≥50% of ≥4 tested paired tumours. Every paired tumour with a defined change, null bound and self bound is tested; there is no activity gate (the activity columns are informational).
    - κ is the mean latent correlation of the MP's genes in the tumour's shared panel (genes detected in ≥5% of the tumour's κ cells under each protocol; panel retention ≥0.5, ≥5 genes), in samples of ≥20 malignant cells, averaged over the tumour's samples under each protocol and floored at 0.01. Every sample of a tumour is scored on the same number of cells (its smallest sample, at most 2000), with its counts binomially thinned to the tumour's lowest median count per cell.
    - A tested tumour is a loss when Δlog κ (test − reference) is below the null bound, below the self bound and below −log 1.5.
    - Null bound, with same-protocol replicate pairs from ≥3 distinct tumours: the lower of two one-sided 95% t prediction bounds, on the linear κ changes (carried to the MP's κ level) and on the log κ changes.
    - Otherwise: the k-th smallest linear κ change, k = floor(0.05 (n + 1)), of the n ≥20 matched random gene sets of the tumour with a defined change (100 drawn per tumour and MP; the 5th of 100), carried to the MP's κ level. A tumour whose k-th smallest change is not below its k-th largest is not tested.
    - Self bound: −1.645 × the bootstrap SE (50 resamples of each sample's cells, drawn per seed, sample and MP) of the MP's own linear κ change, carried to its κ level. A tumour with an undefined SE is not tested.
- **Verdict.**
  - GO if ≥3 retained, non-dosage MPs have nuisance support.
  - INCOMPLETE if no criterion is adequately powered. Powered means: for (a), ≥5 decoys in some technical arm; for (b), for some MP, ≥10 samples with TME cells, with a null and the spillover adjustment; for (c), for some MP, ≥4 tested paired tumours (its Δlog κ, null bound and self bound defined).
  - Otherwise PIVOT.
- **Informational only.** Instability (method variants such as "uncentred" are reported separately and not mixed into `recovery_all`), positional flags outside the dosage rule, the within-sample κ ratio (`kappa_ratio_median`, not spillover-adjusted), the κ activity columns of A4, A4 intensity and prevalence, and the trivial typer.
- **Calibration (toy and synthetic only).** With no protocol effect, an unaffected MP met (c) with probability at most 0.025 on 4 paired tumours, 0.003 on 10 and 0 on 40 (synthetic κ laws, `scripts/calibration/`). If 20 audited MPs were independent, the chance that ≥3 of them meet (c) by chance alone would be at most 0.013 (0.083 for 41 MPs); false (a) and (b) support add to it. On the toy null no gene set met (c), also with the test pieces at a half or a quarter of the reference depth or with 50–60-cell pieces. These figures do not cover regional heterogeneity of a present MP on the random-set path, against which only the log 1.5 margin protects, nor a program absent from part of the pieces on the replicate path, whose t bounds assume near-normal changes: absent from 10–30% of the pieces, it met (c) on 4 paired tumours with probability 0.031–0.047 with 3 replicate tumours (synthetic; 0.065 for ≥3 of 20 MPs at 0.047) (§2.5; [05_pilot_P1.md](05_pilot_P1.md), Caveats).
- **Config check.** Before any computation, `mpaudit` compares the `decision` thresholds and `mpaudit.report.REGISTERED_SETTINGS` (listed under pre-registration below) with their registered values, and checks that every exclusion can run (with re-derived MPs, curated labels and their gene sets; with decoys on, an arm in `decision.dosage_arms`; with the positional axis on, a cytoband and gene positions). Once the cohort is read, and still before any NMF, it checks that the gene positions and the cytoband place genes on chromosome arms. None placed is a deviation; fewer than half placed warns. A GTF whose gene names match no `var_name` counts as missing, and the run falls back to positions in `adata.var`. It refuses a deviation, listing each with its registered value, unless the config sets `unregistered_variant: true`; such a run is labelled "P1-GO v1.3 (unregistered variant: …)" in every output. It also refuses:
  - unknown keys, at the top level and in every section (a misspelt key would otherwise run the default silently), and pre-v1.0 decision keys;
  - keys naming no column of the cohort's `obs` (for the keys only the informational typer reads, a warning);
  - values of the wrong type or out of range: true/false switches; whole numbers such as `protocol.min_pairs`, with `decoy.n_null_floor` ≥20 even in a variant and `protocol.n_null_sets` 0 or ≥20; numbers such as `protocol.kappa_min_effect` ≥0;
  - equal reference and test protocols;
  - missing label, GMT, cytoband or GTF files;
  - decoy arms that neither arm list names, and arm columns missing from the cohort or with missing values in malignant cells (NaN, or a text placeholder such as "nan");
  - expected-TME-shared names for re-derived MPs.

  It warns when the TME axis can run but no expected-TME-shared MP or gene set is pre-declared, when MPs from a GMT come without curated labels, and when the protocol axis is on but no tumour can be paired. Axes that are on but cannot run are noted in `summary.md` and `run_config.json`. The code does **not** compare the rest of the registered file: the contents of the decoy arms, the pre-declared lists, the label and gene-set files, the `gavish` block other than `gavish.seed`, `counts_layer`, `tme.covariates`, `tme.n_boot`, the keys, `positional.chrom_col`/`start_col`, `recurrence_key`, the protocol names, whether the TME and protocol axes are on (off, a criterion is unpowered, which can turn PIVOT into INCOMPLETE), and every other config key. These are registered as filled in, and a run is checked against the frozen file through its `run_config.json` (the normalised config) and the digests of its input files (cytoband, GTF and GMT contents included). A registered label alone does not show that the registered file was used.
- **What the audit runs through.** Decoys run through the Python port (`mpaudit.gavish`), not the original R code. R-versus-port parity is untested: the harness is in `scripts/parity`, but its R side has not been executed, because R could not be installed in the build environment. The planned Sim-I backbone cross-check of the decoys does not exist yet, because Sim-I is not built.

**Pre-registration and study split (month 0).** Before any 3CA expression data are opened, we register on OSF:
- P1-GO v1.3;
- the whole filled-in `mpaudit` config file, frozen with the `mpaudit` commit:
  - the `decision` block: the expected-TME-shared list (planned as IFN/MHC-II response and hypoxia and fixed at registration), the arm lists and all thresholds;
  - `decoy.arms`: the technical decoy arm (doublet score (scDblFinder), ambient/contamination estimate (DecontX), %ribo and, only for studies where it can be computed, intronic fraction) and the dosage arm (CNA clone). 3CA distributes gene-level count matrices without spliced/unspliced counts, so the intronic fraction needs re-alignment of raw reads or BAMs (STARsolo Velocyto features, velocyto or DropletQC). The list of studies that have it is fixed at registration from the P0 census; studies without it carry one constant value per sample, which conditions nothing. Whether `pct_long` joins the arm is decided here too;
  - the rule settings outside the `decision` block, each at exactly its registered value (`mpaudit.report.REGISTERED_SETTINGS`): `min_jaccard` 0.3; `seed` 0; `gavish.seed` 0; `decoy.enabled` true, `decoy.n` 5, `decoy.method` residual, `decoy.n_depth_bins` 1, `decoy.min_kappa` 0.02, `decoy.n_null_floor` 50; `tme.spillover_adjust` true, `tme.rho_sensitivity` 1.5, `tme.null_match` expression+specificity, `tme.n_null` 500, `tme.min_cells` 20, `tme.n_markers` 50, `tme.min_marker_fold` 8; `protocol.kappa_min_effect` log 1.5, `protocol.min_pairs` 3, `protocol.n_null_sets` 100; `positional.enabled` true, with the cytoband and GTF. A null size is not an "at least" setting: another size moves the estimated null quantile or background, so choosing it would give a second chance to pass. It adds or drops sets, each gene set (`tme.n_null`), sample and gene set (`decoy.n_null_floor`) or tumour and MP (`protocol.n_null_sets`) having its own random stream. `protocol.kappa_active` is informational and not registered. The constants fixed in code (`mpaudit.report.FIXED_CONSTANTS`: 50 bootstrap resamples and z 1.645 for the self bound, the 5% level and the minimum of 20 sets of the random-set bound, the κ floor 0.01, the cell cap 2000 and (c)'s panel thresholds, and the 3-SD entry, the minimum of 20 sets and the gene-set thresholds of κ_retained) are part of the commit;
  - the other NMF settings (`gavish` block);
- the curated-artefact labels (`typer.labels`: the 11 3CA v2 `inds_rm` MPs, label 1) and their gene sets (`typer.labelled_gmt`, still to be built, see above), with the sha256 of each file (every run writes the digests of its input files into `run_config.json`);
- a **split of the 3CA studies** into an audit half and a held-out confirmation half.

The config check refuses a run that deviates from the registered settings unless it is declared an unregistered variant (`unregistered_variant: true`), which labels its verdict as one.

The split is drawn at random by study, stratified by cancer type and protocol. It is by study rather than by tumour because tumours of one study share lab and protocol. P1, P6, P7, the parity samples, the month-10 audit paper and all method development before the M14 registration use only the audit half, plus simulations (TCCA stays held out). The later TESSERA pre-registration (M14: frozen defaults, margins, type rules, orthogonal endpoints) is kept, but it is tested only on the confirmation half, on TCCA and on the new cohort. Nothing before pre-registration 2 touches the confirmation half, so neither the audit nor the method development reuses the data that the pre-registered analyses test. The cost is power: P1 sees about half the twins and samples, which is why P0 checks each half.

**Decision at month 3** (P0, P1):
- P0 and P1 GO → continue to the method build.
- P1 PIVOT → pivot.
- P1 INCOMPLETE → close the data gap and rerun under the same rule.

**Decision at month 6** (P2, P3):
- P2 and P3 pass → full plan.
- P3 fails after one revision → publish the audit plus the benchmark and stop method development.

P4 and P5 can demote typing or A3 from the headline, but they do not stop the project.

**Implementation status (2026-10-18).** `mpaudit` (this repository; [05_pilot_P1.md](05_pilot_P1.md)) implements P1:
- the Python port of the Gavish/3CA pipeline (`gavish.py`, `nmf.py`), checked line by line against tiroshlab/3ca and renozao/NMF;
- per-arm decoys with the excess κ_retained (`decoy.py`), with %ribo, intronic fraction and, optionally, `pct_long` as continuous covariates (`data.py`);
- A1 with spillover subtraction (`tme.py`);
- κ-based A4 with protocol-specific controls, the replicate and random-set null bounds, the bootstrap self bound and equal κ cells per tumour (`protocol.py`);
- instability and positional arm enrichment;
- a trivial typer with gene-length, %ribo and intronic features (`typer.py`);
- the P1-GO v1.3 decision (`report.py`) and the config check against the registered settings (`pipeline.py`); every run records the rule label, any deviations, the sha256 of its input files, the axes that did not run, the code hashes of the NMF cache and the capped or failed NMF fits of every run (`run_config.json`, `summary.md`);
- the calibration harness of criteria (c) and (a) on the toy and on synthetic κ laws, with v1.2 and v1.3 reference results (`scripts/calibration/`).

The toy simulator (`sim.py`) was recalibrated on 2026-10-17. At their first magnitudes its abundance/%ribo and long-gene axes were invisible to the standard pipeline (no per-sample NMF program had more than 1 abundance or 6 long genes in its top 50). The target was fixed before any audit, decision rule or demo was run, and it concerned standard-pipeline visibility only: at the demo settings, each axis must put ≥10 of its genes into one robust program in ≥8 of the 16 samples. The abundance genes now scale with the square of the planted cytoplasmic factor and the long genes with 1 + 20 × the nuclear fraction, with every random draw unchanged. At seed 0, both axes now form their own MPs (robust in 14 and 12 samples, in a fresh reference run of the demo at commit edb2736, whose NMF and simulator v1.3 keeps), and EMT, contamination and the CNA lose 2, 3 and 4 robust samples to the competition for ranks 4–9. The calibration run had reported 13, 11 and a loss of 3 for the CNA, because the NMF cache of that time served one sample stale fits of the first magnitudes; the cache key now hashes the counts themselves. The calibration and its NMF check used seed 0 only. A no-NMF proxy (`tests/test_sim.py`) puts seed 0 highest of seeds 0–5 for the long-gene axis, so visibility at other seeds is not established. Background κ of random gene sets rises from about 0.08 to 0.21, because a third of the library sits on one per-cell factor (30-gene sets of unplanted genes detected in ≥5% of a sample's malignant cells; median over sets, then samples). Details are in the `mpaudit.sim` docstring.

Findings and state:
- **NMF tolerance.** The port now defaults to tol = 1e-7, because at 1e-5 the top-50 lists were not converged. Toy result (on the toy before its technical axes were recalibrated): we continued the same trajectory to tol 1e-8. Against that reference, a program kept a median of 44 of its top 50 genes at 1e-5, 48 at 1e-6 and 49 at 1e-7. Tol 1e-7 costs about 3× the NMF time of 1e-5. R-parity runs use `gavish.tol: 1.0e-8` (written with a decimal point, so YAML reads a number). The iteration cap, raised to 5,000, never bound. R's own stopping rule, applied to the same iterates, would stop after a median of about 920 iterations, about 5× later than tol 1e-5.
- **R parity.** `scripts/parity` holds three steps:
  1. export the exact NMF input;
  2. `run_r_snmf.R` calls `NMF::nmf(method = "snmf/r")`;
  3. `compare.py` gives PASS / WITHIN NOISE / FAIL / INVALID at the program, robust-program and MP level.

  The Python side is tested; the R side has never been run, and there are no parity numbers yet. Until parity passes on 3–5 real samples (months 0–1), absolute numbers from the port (MP count, MP genes, decoy recruitment) are not 3CA-equivalent.
- **Toy behaviour of each axis** (`mpaudit.sim`, recalibrated toy, default cohort, seed 0 unless stated; true by construction, not evidence about tumours):
  - *A1.* After spillover subtraction, only stress and IFN are above the null q99 (values under P6; IFN 0.89 vs 0.78). The leak is largely removed: hypoxia goes from 0.49 unadjusted to 0.11, EMT from 0.67 to 0.26 and CNA_7q from 0.18 to 0.08 (with a negative slope); without malignant ambient RNA, hypoxia and CNA_7q are at 0.00 unadjusted. IFN can only be excluded by the expected-TME-shared list. Stress and IFN also stay above the null of the ρ̂ × 1.5 re-run (0.98 vs 0.82 and 0.88 vs 0.76); on a 16-sample fresh toy, 0.95 vs 0.33 and 0.82 vs 0.31 (0.95 vs 0.34 and 0.84 vs 0.36 with the technical axes at their first magnitudes, the cohort of the spillover tests in `tests/test_tme.py`), while CNA, hypoxia and cycle have negative slopes there.
  - *A4* (P1-GO v1.3, re-run 2026-10-18 with equal depth). Stress loses κ in 4 of 4 paired tumours (κ 0.67–0.77 fresh vs 0.13–0.22 fixed). None of the other eight planted sets (cycle, hypoxia, IFN, EMT, CNA, contamination, abundance, long genes) loses κ, although all 32 of their paired-tumour tests are tested, 3 of them in a tumour that lacks the program (CNA_7q in T01 and T03, IFN in T03). The null bound alone would call 5 losses and 7 gains; the self bound leaves 3 losses and 3 gains; the margin alone leaves 1 loss (CNA_7q in T03: Δlog κ −0.61 against a null bound of −0.55), which the self bound (−0.68) blocks (§2.5).
  - *Decoys* (P1-GO v1.3's excess κ_retained; decoy seeds 1000 and 1001, the `decoy_recruitment` default with config seed 0; mean of the two decoys). The registered technical arm (estimated doublet score and contamination, %ribo, intronic fraction) keeps 1.39 of the contamination set's excess κ, 0.88 of abundance and 0.92 of long genes, and between −0.04 and 0.07 of cycle, hypoxia, EMT, IFN and CNA (v1.2's plain ratio gave 0.12–0.24 for cycle, hypoxia and CNA, most of it the shared dilution background). Values above 1 arise because the contamination covariate is built partly from the set's own genes (18 of its 36 TME-specific genes are in the set): an arm on `contamination_est` alone keeps 1.97 of that set's excess and one on `doublet_est` alone 1.75. The contamination set enters in only 7 of the 16 samples (its excess is 1.7–3.7 SD of the random sets'). The dosage arm (estimated clone) keeps 1.00 of the CNA set's excess and between −0.01 and 0.01 of the others.
  - *Typer* (seeds 0–3). For the abundance set, \|r\| with %ribo is 0.95–0.96 (next highest 0.19–0.37, the long-gene set). The long-gene set has a median gene length of 351–454 kb and \|r\| with intronic fraction of 0.83–0.86 (others ≤0.09).
- **Not built.** Sim-I, κ cross-fitting and the certification layer, the evidence-vector typer, allele-aware dosage (only an expression-only positional flag exists) and the TESSERA engine.

**Publishable even if the full method fails:**
1. An audit paper on the audit half of the studies: "How many pan-cancer metaprograms survive nuisance-matched decoys, within-sample negative controls and protocol contrasts?" Ready by month 10. Venue: *Nature Methods* Analysis, *Genome Biology* or *Cell Systems*.
2. The Sim-I/Sim-R benchmark as a resource.
3. The matched multi-protocol tumour cohort as a data resource (*Scientific Data* or *Genome Biology*).

---

## 9. Timeline, team, compute and preprint strategy

| Months | Milestone |
|---|---|
| 0 | **OSF pre-registration 1**: P1-GO v1.3, the whole filled-in `mpaudit` config with the commit, the curated-label files, the audit/confirmation study split; DAC applications; P0 metadata census |
| 0–3 | P0, P1 (`mpaudit`, audit half), P6 (audit half), P8; R-versus-port parity on 3–5 real samples of the audit half (`scripts/parity`); IRB approval and wet-lab protocol optimisation |
| 3–6 | κ cross-fitting and Sim-I v0, then P2 and P3; P7 starts; P4 once ≥40 labels and the evidence-vector typer exist; P5 when BAM access arrives (possibly to M9) |
| 3–9 | v0.1: split, certify, genes, audit mode, decoy audit (developed on simulations and the audit half only); Sim-I/Sim-R v1; **audit manuscript (audit half) ready internally at M10** |
| 6–16 | Matched cohort: 24 tumours collected at about 2–3 per month, then sequencing |
| 9–15 | TesseraNMF v1; full baseline harness; benchmark on the fractional design; sensitivity tables |
| 14 | **OSF pre-registration 2**: frozen defaults, margins, type rules, all orthogonal endpoints, tested on the held-out confirmation studies |
| 15–20 | Atlas run; TCCA; allele subset; Tahoe; spatial |
| 18–22 | Scoring of matched-cohort predictions (pre-registered) |
| 22–26 | Reference v1.0, package v1.0, manuscript; **submission about M24–26**; revision to about M32 |

The submission target is unchanged, because the prerequisites of P2–P5 (cross-fitting, Sim-I, labels, allele data) were already scheduled in months 3–9. What moves is the method gate, from month 3 to about month 6, which leaves less slack.

**Team (about 4 FTE):**
- PI (20%);
- statistical methods postdoc (100%): certify, genes, typing, theory remarks;
- computational biologist (100%): data, CNA/allele, simulators, validation;
- research software engineer (75%): engine, GPU, package, R bridge, CI;
- bioinformatics or data engineer (50%): ingestion, STARsolo re-processing, containers;
- statistics co-advisor (25%), who also serves as blinded adversary;
- wet-lab scientist (50% for 12 months) plus a pathology/surgical collaborator.

**Compute:**
- 2–4 H100s; 2,000–4,000 GPU-hours;
- 50–100k CPU core-hours (the Gavish pipeline on real and decoy data, through the Python port plus the R code for parity; GBCD; Numbat);
- 20 TB storage.

**Preprint strategy and scoop risk.**
- Default: one preprint at submission.
- **Tripwires** that trigger posting the audit plus the certification layer within 4–6 weeks (possible from M10):
  - any preprint on decoy or synthetic-null certification of gene programs;
  - a single-cell port of multi-study NMF;
  - GBCD with recurrence testing;
  - a malignant TCAT reference;
  - an ITHmapper journal version with uncertainty.
- Keep the package name registered and the repository private until the preprint.

**Fallback venues:**
- a working certification layer with a weak engine → still *Nature Methods*, because the headline is the certification;
- weak protocol evidence but a strong audit → *Nature Methods* Analysis / *Genome Biology*;
- negative audit → *Genome Biology* software, plus *Nature Cancer* for the matched-cohort biology if the stress and selection results hold.

---

## 10. Risks and mitigations, including every referee fatal flaw

| # | Referee fatal flaw or major risk (source) | Resolution in TESSERA |
|---|---|---|
| 1 | Decoys are non-competitive: the "existence" q-value rejects for any co-expression (CREDO F1; CERT-MP rigor) | Synthetic decoys are removed from inference. Existence = effect-size, specific, cross-fitted membership **plus** held-out-cohort replication. Decoys are kept only as an audit of pipelines that do not model nuisance. |
| 2 | Decoys drawn from the conditioning model are circular or vacuous; thinning asymmetry (CERT-MP F1, F2) | Same as row 1. Inference uses cell splitting, so no NB thinning enters inference. |
| 3 | Knockoff+/target-decoy FDR invalid: tumour swaps flip all W; equal-chance assumption violated (CERT-MP F3; CREDO F1) | Knockoff layer dropped. Within-MP BH, partial-conjunction bounds and replication are used; PRDS assumption stated, with a Bonferroni sensitivity analysis. |
| 4 | Self-validating simulation, i.e. the inverse crime (CREDO F2; InvarMP; CERT-MP) | Primary simulator Sim-I uses an independent generator with an injector whose functional forms differ from the model's. Sim-R plants into real data. Adversarial scenarios are designed blind. The SCRP backbone is reported only as a best case. |
| 5 | Thresholds relocated rather than removed (CREDO F3; CALIPER) | Accepted. A minimal set of interpretable constants, frozen and pre-registered, plus a full sensitivity table. Tier-move rules are removed (no accept/reject moves). |
| 6 | Cohort × lineage × protocol aliasing (CREDO F4; InvarMP F2/F4) | P0 design census as a gate. A4 is tested only where a bridge exists. The "C cohort-restricted" type replaces demotion to artefact. The new matched cohort supplies bridges. |
| 7 | Tumour-wide programs vs library batch; decoys keep tumour means (CREDO rigor F2) | The primary estimand is within-tumour programs. The between-tumour module is optional and its claims are restricted to multi-library tumours plus cross-cohort replication. The limitation is stated explicitly. |
| 8 | Count splitting does not break CNA circularity (CREDO F3; CERT-MP) | Dosage offsets only on allele-anchored segments; expression-only calls give a positional flag only. Positionally clustered biological programs serve as negative controls. Gene- and lineage-specific γ. Validation on ≥60 allele or scDNA tumours. |
| 9 | Point-null membership is non-specific and saturates with n (CALIPER F1; CERT-MP F2) | κ estimates a latent correlation, whose estimand is depth-independent (its moment estimate falls with depth, so κ is compared at equal depth), specific given all other programs, tested against a within-tumour matched null band with a floor. In a toy with a perfect dictionary, conditioning did the work (≤1% false membership up to 30/40 shared genes); the null band added nothing, and the floor alone fails from about 20/40. P2 gates leakage ≤5% with an estimated dictionary, and depth invariance. |
| 10 | Typing by non-rejection; order-dependent ladder (CALIPER F2; InvarMP rigor F1) | TOST equivalence per axis; Shapley order-invariant attribution; "indeterminate" outcome; B1 only on positive evidence. |
| 11 | False or non-novel thinning theorems; Theorem 2 vacuous (CALIPER F3) | No thinning theory is claimed. Thinning is used only for model selection. Neufeld et al. and Hivert et al. are cited. |
| 12 | Engine-indexed estimand contradicts "MP as estimand"; δ degenerate; s_min contradiction (CALIPER F4) | Location = engine proposal; the certificate is engine-free; engine invariance is reported empirically. Fixed biological match tolerance, self-similarity calibrated on independent cell folds. Single-cohort programs form a separate type, not a subtype of "exists". |
| 13 | Fold C used twice (CALIPER rigor F3) | Three disjoint cell folds: discovery, membership, genes. |
| 14 | Transport κ absorbs the dissociation artefact through responsive TME calibrators; κ vs δ not identified for malignant-only genes (InvarMP F1/F3) | No gene-wise transport is used for invariance tests. κ's estimand (a correlation) is invariant to gene-wise capture efficiency; its estimate falls where counts fall, so A4 compares κ at equal depth and on the genes both protocols detect. TME cells serve as *negative controls* (A1), not as calibrators, and their responsiveness is the signal being measured. |
| 15 | e-values invalid or powerless; e-BH arithmetic (InvarMP F3) | No universal-inference e-values. The saved referee runs show why: the e-value as specified stayed valid but did not detect a 1.65-fold shift, and the plug-in variant that has power fired under the null whenever gene-wise capture was misestimated ([referee_checks](referee_checks/README.md)). Paired tumour-level tests on a small, pre-registered hypothesis family (≤10–30 MP × axis pairs) with a stated power. |
| 16 | Pan-cancer invariance claim untestable on public data (InvarMP F2/F4) | The headline is restricted to bridged axes plus the new cohort. Elsewhere the label is "U-no-bridge", reported honestly in the design-graph report. |
| 17 | Thinning null twins understate the real null (regional ITH) (InvarMP rigor) | The protocol null comes from same-protocol multi-region pieces (8 tumours) plus public multi-region data. Where replicates are missing, the null is random matched gene sets in the same tumours, which ignores piece-to-piece variation and the MP's own sampling noise, so a κ loss must also exceed log 1.5; this is reported as the weaker null. |
| 18 | QC covariates are post-treatment and absorb biology (CERT-MP novelty) | QC is never an offset. It is a reported A5 axis with Shapley share. |
| 19 | Doublet prior cannot represent misclassified whole TME cells; mimicry programs mis-typed (CREDO rigor; CERT-MP) | Unconstrained usage of TME components. A2 requires a cell-level doublet or ambient association, not loading similarity. Planted mimicry tests in Sim-R. |
| 20 | Underpowered validation (11-positive AUROC; n ≈ 6 matched tumours) | ≥40 labelled artefacts with CIs and a trivial baseline; matched cohort n = 24 with a power calculation; non-canonical Tahoe endpoints. |
| 21 | Validation is non-comparative | Every rung has a pre-specified endpoint, compared head-to-head with Gavish/3CA, GBCD and ITHmapper MPs. |
| 22 | Scope and team about 2× under-resourced (all) | Scope cut: no SVI or IBP model, no power posterior, no normal-lineage nodes, no γ_g resource, no 17-method harness. About 4 FTE; about 26 months to submission. |
| 23 | Data availability (raw UMIs, empty droplets, BAMs) overestimated | P0 gate; DAC applications at M0; ambient proxy flagged; analyses reported with and without each optional module. |
| 24 | The flagship confirms known biology | P1 gate; flagship built on per-tumour quantitative provenance and non-obvious selection/induction predictions. |
| 25 | Wet-lab cohort fails (accrual, cost) | Minimal cohort of 12 tumours; PDX substitute; public twins; a computational-only claim limited to bridged axes. |
| 26 | Scoop | Engine-agnostic claim; tripwire preprint at M10; fortnightly monitoring. |
| 27 | The P1 gate could pass on almost any cohort: expected biology and CNA MPs counted, one decoy enough, fixed R² below the null, three versions of the rule (review 2026-10-16) | One frozen rule, P1-GO v1.3, stated in full in §8 and registered at M0 together with the whole config, the `mpaudit` commit and the curated-label files. The decision thresholds, the registered settings outside the `decision` block and the inputs the exclusions need are checked before any run; a deviation is refused or, if declared, labelled an unregistered variant. The rest of the frozen file (arm contents, lists, label files, NMF settings other than `gavish.seed`) is checked by comparing the run's `run_config.json` with it. |
| 28 | Malignant→TME ambient spillover makes malignant programs look TME-shared (review 2026-10-16) | Predicted spillover (ρ̂ from malignant-specific markers) is subtracted from pseudo-bulk TME scores, with an expression- plus specificity-matched null. P1 counts A1 only where the adjustment ran and survives ρ̂ × 1.5. The κ ratio is not yet adjusted; this is stated as open, and until it is adjusted it cannot make A1 supported (in P1 it is informational). |
| 29 | The audit and the pre-registered analyses reuse the same 3CA data (review 2026-10-16) | Studies are split at M0 into audit and confirmation halves. P1, P6, P7, the parity samples, the audit paper and all method development before M14 use the audit half only (plus simulations), and pre-registration 2 is tested on the confirmation half, TCCA and the new cohort. |
| 30 | Decoy results describe the Python port, not the R pipeline; parity untested | R-versus-port parity on 3–5 real samples (`scripts/parity`) before P1 numbers are quoted as 3CA-equivalent; NMF tolerance tightened to 1e-7. |

---

## 11. Decisions the user must make

1. **Wet-lab partner and budget.** Choose one:
   - full cohort (24 tumours, about USD 300–450k);
   - minimal cohort (12 tumours, about USD 150–200k);
   - none (odds for *Nature Methods* drop to about 15–20%).

   This is the single biggest lever.
2. **Cancer types for the matched cohort.** I recommend breast, colorectal and NSCLC: they overlap public Flex and Xenium data, and resections are large enough to split.
3. **Python-first (recommended) or R-first.** Python-first means a scverse core with the `tesseraR` bridge. If the lab's user base and expertise are R/Seurat-centred, the engine still needs PyTorch; R-first would mean a torch-for-R port and losing scvi-tools adapters.
4. **Controlled-access data.** Approve EGA and dbGaP applications at M0 (O'Flanagan; BAM-bearing 3CA studies) and accept 3–6 month delays.
5. **Collaborators.**
   - A statistics co-advisor: e.g. a data-thinning or selective-inference group.
   - A blinded adversary for simulator scenarios.
   - Ideally a data-sharing contact for raw 3CA matrices.
   - Decide whether to approach the GBCD or Tirosh groups as collaborators, which lowers scoop risk but shares credit.
6. **Scope of the between-tumour module.** Ship it as an optional, limited-claim module (recommended), or omit it.
7. **Pre-registration on OSF, in two steps.**
   - At M0, before any 3CA expression data are opened, register P1-GO v1.3; the whole filled-in `mpaudit` config with its commit (the `decision` block with the expected-TME-shared list, the decoy arms, the rule settings outside the `decision` block and the NMF settings); the curated-label files (`typer.labels`, `typer.labelled_gmt`); and the audit/confirmation study split.
   - At M14, register the TESSERA defaults, margins, type rules and endpoints, to be tested on the held-out confirmation studies.

   Recommended: it constrains the analysis, but it is the best defence against cherry-picking critiques. The user must approve the expected-TME-shared list before M0.
8. **Preprint policy.** Single preprint at submission with tripwires (recommended), or an early audit preprint.
9. **Name** after the collision check.
10. **Team hiring:** whether the RSE is a hire or shared core staff.

---

## 12. References

**A. Verified in the landscape** ([V], [V+] or [P] as tagged there):
- Gavish et al., Nature 618:598–606 (2023) [V]
- Tyler et al., 3CA v2, Nature Cancer 6:1088–1101 (2025) [V]; code [V+]
- Kinker et al., Nat Genet 52:1208–1218 (2020) [V]
- Barkley et al., Nat Genet 54:1192–1201 (2022) [V]
- GeneNMF (CRAN) [V]
- Liu…Stephens, GBCD, Nat Genet 57:263–273 (2025) [V+]
- Kotliar…Raychaudhuri, TCAT/starCAT, Nat Methods 22:1964–1980 (2025) [V] (re-checked in A2)
- Nofech-Mozes…Abelson, ITHmapper, bioRxiv 2026 [V] (the 2026-10-16 review confirmed that it exists but could not verify its counts; October 2026 searches did not find it)
- Grabski, Trippa, Parmigiani, Genome Biol 26:98 (2025) [V]
- Tetris, AoAS 2023 [V]; APAFA, JASA 2026 [V]; SUFA, JASA 2025 [V]
- Xue et al., BayesPowerNMF, PLoS CB 2026 [V]
- Neufeld et al., data thinning, JMLR 2024 [V]; count splitting, Biostatistics 2024 / arXiv 2307.12985 [V]
- Hivert et al., arXiv 2405.13591 [V]
- ClusterDE, Song, Li et al. [V]
- recall, DenAdel…Crawford, AJHG 2025 [V]
- Fu/Huang/Sidiropoulos, arXiv 1709.00614 and 1803.01257 [V; venues P]
- scDEF [V; peer review unknown]
- consensus scHPF [V]
- cNMF, eLife 2019 [V]
- LIGER iNMF [V]
- DRVI [V]
- Spectra, Nat Biotechnol 2023/24 [V]
- sciRED, Nat Commun 2025 [V]
- CSI-GEP, Cell Genomics 2025 [V]
- KINOMO [V]
- CanSig benchmark, Cancer Res 86(4):873 (2026) [V, one sweep]
- TCCA, bioRxiv 2026 [V] (details re-checked in A2)
- Tahoe-100M [V]
- Numbat, CopyKAT, SCEVAN [V]
- CellBender, SoupX, DecontX [V]
- scDblFinder [V]
- clonealign [V]
- CalicoST [V+; journal unconfirmed]
- Janesick 2023 [V; venue P]
- Oliveira, Nat Genet 2025 [V] (re-checked in A2: FFPE blocks, no fresh-dissociated malignant cells found; not a fresh-versus-fixed twin)
- Splatter, splatPop, scDesign3 [V]
- Venet et al., PLoS CB 2011 [V]
- Hotspot [V]
- MOFA+ [V]
- Replogle 2022 [P]
- Greenwald GBM spatial [V; venue P]

**A2. Checked for the 2026-10-17 revision.** These were checked by web search. Publisher, PMC, bioRxiv and arXiv pages could not be reached, so abstract wording comes from search snippets unless stated. GitHub content was read directly. Tags are as in the landscape: [V] found online, [V+] re-checked from source code/README, [P] partly verified, [U] unverified.
- Baron et al., "The Stress-Like Cancer Cell State Is a Consistent Component of Tumorigenesis", Cell Systems 11(5):536–546.e7 (2020) [V]. One abstract sentence was matched by exact-phrase search. The paper argues that the stress-like state is real *in vivo* (zebrafish melanoma, spatial transcriptomics, human samples); it does not establish dissociation stress.
- O'Flanagan et al., Genome Biol 20:210 (2019) [V]. Cold protease (6 °C) vs collagenase (37 °C); 155,165 cells from patient tumours, breast-cancer PDXs and cell lines; a 512-gene stress core. The patient tumour types and the EGA accession were not verified.
- van den Brink et al., Nat Methods 14:935–936 (2017) [V]. Dissociation-induced expression in muscle stem cells.
- Neuschulz et al., Mol Syst Biol 19:e11147 (2023; online December 2022) [V]. 4sU labelling during dissociation; the data seen were from mouse and zebrafish, and no human tumour data were found.
- Jiménez-Gracia et al., FixNCut, Genome Biol 25:81 (2024) [V]. Tumour content and the size of the stress reduction were not checked.
- Denisenko et al., Genome Biol 2020; Slyper et al., Nat Med 2020; Marsh et al., Nat Neurosci 2022 [P]. Their existence and topic were confirmed by search; the papers were not read.
- Trinks et al., Cell Oncol 47:1221–1231 (2024) [P]. Patient-matched fresh, cryopreserved and FFPE lung cancer. The patient count and Flex use were not checked.
- Xirenayi et al., ccRCC fresh / frozen snMultiome / FFPE snFlex, bioRxiv 10.1101/2025.11.25.690485 (2025) [P]. Patient matching, patient numbers and any stress or malignant-program results were not checked.
- Oliveira et al., Nat Genet 57:1512–1523 (2025) [V]. Flex on 8 FFPE blocks (tumour and normal), Visium HD on 5 patients and Xenium on 4 tumours, all from the same FFPE blocks; the GitHub README lists FFPE as the only tissue preparation. No fresh-dissociated malignant cells were found. The material used for a T-cell immune-profiling arm mentioned on a 10x page was not confirmed, and the review's sub-claim that "the only fresh cells are sorted T cells" is unverified.
- Kotliar et al., TCAT/starCAT, Nat Methods 22:1964–1980 (2025) [V]. The sentence on starCAT generalising to other cell types and tissues is near-verbatim from search results and was not checked on the article. The Flex caveat comes from the GitHub README [V+].
- González-Bermejo et al., TCCA, bioRxiv 10.64898/2026.02.18.706316 (2026) [P]. 36 studies, 34 tumour types, about 537–540 patients plus 183 cell-line samples (sources differ), 43 MPs. The NMF pipeline and the hand-coded MP whitelist were read from GitHub code [V+]. The overlap with 3CA is not stated.
- ISOMERIC, bioRxiv 10.1101/2025.04.21.649842 (2025) [V]. Authors and publication status were not checked.
- Denault, Tayeb, Carbonetto, Willwerscheid, Stephens, cEBMF, arXiv 2505.11639, NeurIPS 2025 [V]. No application to tumour MPs was found; calling it a "GBCD successor" is an interpretation, not the authors'.
- Xia, Cao, Gao, DECIPHER, Nat Commun 16 (2025), doi 10.1038/s41467-025-63140-8 [V].
- Dong, Kluger, Fan, Kluger, SIMVI, Nat Commun 16:2990 (2025) [V].
- Wang, Lin, Zhao, FDR control via data splitting for testing after clustering, arXiv 2410.06451 (2024) [V].
- Xu…Ma, Gene Program Dictionary, bioRxiv 10.1101/2025.10.29.685322 (2025) [P]. Details come from search snippets.
- scPFA, lineage-informed factor analysis of heritable programs, bioRxiv 10.64898/2026.09.21.753290 (2026) [U]. Found by search; the abstract was not retrieved.
- Possible extra protocol-twin data [P]: a GEM-X Flex performance preprint, bioRxiv 10.64898/2026.03.06.710062, and a Frontiers 2026 brief report on matched FFPE and frozen samples (reportedly 12 pairs; tumour types unknown). Neither was checked.
- Tyler et al., 3CA v2: the publication in Nature Cancer (2025) was re-confirmed; the citation in A is unchanged.

**B. Cited by referees as verified in their own sessions, not in the landscape and not re-checked by the PI. Verify before writing:**
- CS-CORE, Nat Commun 14:4846 (2023)
- Clipper, Genome Biol 2021
- sc-SHC, Grabski, Street, Irizarry, Nat Methods 20:1196–1202 (2023)
- Langfelder et al., PLoS CB 2011; NetRep, Cell Systems 2016
- Wang, Xu, Zhao, Wang, data thinning for Poisson factor models, JASA (2025/26)
- Nomura et al. / Tirosh lab, Nat Genet 57:1155–1167 (2025) (landscape listed it as [P])
- CLCC, bioRxiv 10.64898/2026.03.01.708901
- subclonal-dosage preprint, bioRxiv 10.64898/2026.03.15.711906
- TreeAlign (venue unverified)
- Marsh et al., Nat Neurosci 25:306–316 (2022) (existence confirmed in A2; volume and pages not re-checked)
- inVAE, bioRxiv 2024
- sysVI, BMC Genomics 2025
- Weine et al., log1p Poisson NMF, arXiv 2601.05845
- Miller & Dunson, coarsened posteriors, JASA 2019
- Bogomolov & Heller, Stat Sci 2023
- Emery et al., arXiv 1907.01458
- Sarkar & Stephens, Nat Genet 2021
- Hong et al. GPAS NSCLC atlas, bioRxiv 2025.05.06.652317
- multiGroupVI, MLCB 2022
- Xenium contamination paper, Nat Methods 2026 (title only)
- the multi-platform COAD/HCC/OV benchmark (Nat Commun 2025) and the breast/NSCLC Xenium benchmark (Nat Methods 2026)

**C. Unverified, from memory; check before citing:**
- Benjamini & Heller 2008 (partial conjunction)
- Goeman & Solari 2011
- TREAT (McCarthy & Smyth)
- DerSimonian–Laird; HKSJ
- Ren & Barber e-value derandomisation (arXiv 2205.15461 listed by a referee as [V])
- empirical-Bayes deconvolution / ashr
- Shah & Peters 2020
- Elias & Gygi 2007
- CytoTRACE

(The dissociation papers van den Brink 2017, Denisenko 2020 and Slyper 2020, previously listed here, moved to A2.)

The scoop search of 2026-10-02 was proxy-limited, so the absence of a direct precedent is weak evidence. A second search in October 2026 found no published or preprint method that does any of the following:
- calibrates error rates for recurrent tumour metaprograms;
- types them by provenance against within-sample negative controls;
- audits the 3CA/Gavish MPs.

That search used about 14 targeted web searches, without bioRxiv full-text search, PubMed or citing-paper lists, so it is still weak evidence. The nearest work does one of three things:
- separates CN-dependent variation within tumours (ISOMERIC);
- disentangles spatial context (SIMVI, DECIPHER);
- controls FDR for post-clustering tests (Wang et al. 2024).

TCCA repeats the standard pipeline with manual curation. The manual check in P8 is mandatory.