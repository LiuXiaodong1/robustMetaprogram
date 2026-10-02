# TESSERA: a proposal for calibrated tumour metaprograms with typed provenance

**Draft project proposal (AI-assisted design synthesis; to be reviewed by the PI). Target: *Nature Methods* Article.**
Date: 2026-10-02. Inputs: the landscape and gap analysis ([01_landscape_and_gap_analysis.md](01_landscape_and_gap_analysis.md)), four method concepts (CREDO, CALIPER, InvarMP, CERT-MP), and eight adversarial referee reports ([03_design_panel.md](03_design_panel.md)).

**Working name: TESSERA** (Typed, Effect-Sized, Specific, Split-validated, Externally Replicated, Audited metaprograms). A tessera is one tile of a mosaic.
- Python package `tessera-mp` (`import tessera`); R bridge `tesseraR`.
- The name has **not** been checked for collisions. Fallback names are `mpcert` and `CREDO`.

---

## 0. Bottom line

- **Can it reach *Nature Methods*? Yes, but only if three conditions hold.**
  1. Pilot P1 shows that at least 3 MPs that curators *retained* in 3CA are materially driven by nuisance, protocol or dosage, or are unstable. Rediscovering the 11 MPs already removed by hand does not count.
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
  - Mitigation: pilot P1 is a hard go/no-go gate.
  - The flagship is built around results that are *quantitative per tumour*, not around overturning MPs: the ex vivo share of each MP in each tumour, selection versus induction, and malignant-intrinsic versus TME-shared.
- **Resources.** About 4 FTE for about 26 months to submission; 2,000–4,000 GPU-hours; a wet-lab cohort costing about USD 300–450k at full scope (an estimate to be quoted), with a minimal version at about USD 150–200k.

---

## 1. The conceptual advance

### 1.1 What is wrong with "recurrence"

In the Gavish/3CA pipeline, "robust" means "recurrent". Its cross-tumour filter (≥10 of 50 shared genes from a 7,000-gene universe) guards against a random-gene null with P ≈ 1e-12 (landscape §1.4). Anything that recurs for a systematic reason passes: ambient RNA, doublets, dissociation stress, quality axes, positional programs. Every published variant then ends in manual curation (3CA v2 `inds_rm`, GeneNMF `dropMetaPrograms`).

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

   where R^x is the gene–gene correlation of *latent* expression. It is estimated with a measurement-error-corrected moment estimator in the style of CS-CORE (Poisson/NB measurement model with size factors), so it does not depend on sequencing depth.

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
- **Membership is defined by effect size and specificity, so it does not inflate with sample size.** A referee toy simulation ([referee_checks/ref1_check.py](referee_checks/ref1_check.py): 1,000 cells, NB size 3, unrelated program at fold change 2.5 in 30% of cells) showed that a point-null cross-fold test declares a tumour to "express" m in 72% of cases when m shares 5 of 40 genes with an unrelated program, and in 100% at 10 of 40. A second referee found the same direction with weaker effects (about 20–50% at 4–8 of 40; [referee_checks/ref2_sim2.out](referee_checks/ref2_sim2.out)). Conditioning on all other programs, plus an effect-size null calibrated within the tumour, removes that failure.
- **Real negative controls replace synthetic nulls for provenance.** TME cells in the same droplet suspension share ambient RNA, dissociation time, ischaemia and RNA quality with the malignant cells, but not malignant programs. This control exists in essentially every 10x tumour sample. It turns "is this technical?" into a measured quantity rather than a curator's judgement.
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
- **QC covariates are not offsets.** %mito, %ribo and complexity are downstream of biology (cell cycle, OXPHOS, translation), so they are used only as a reported "QC-entangled" axis.

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
| **A1 Compartment sharing** (the within-sample negative control) | Two quantities: (i) κ of k's genes in co-resident T/NK, myeloid, fibroblast and endothelial cells, after removing TME identity programs; (ii) partial R² of per-sample malignant prevalence on the same gene set's TME activity, conditioned on cancer type with study random effects | Upper CI of the TME/malignant κ ratio < 0.3 **and** partial R² upper CI < 0.1 | Lower CI of the ratio > 0.6 **or** partial R² lower CI > 0.3 |
| **A2 Contamination/doublet** | Shapley share of κ_tk removed by the M_se and A_s terms, **together with** the cell-level association of k-activity with doublet score and ambient fraction | Share upper CI < 0.2 | Share lower CI > 0.5 **and** a positive doublet/ambient association. Loading similarity to TME cells alone is never enough, which protects mimicry programs such as pEMT, MES-like and MHC-II |
| **A3 Genetic** | (a) dosage share: Shapley share of κ removed by **allele-anchored** D_cg; (b) clone η² of activity with a region-aware block permutation; (c) arm enrichment, with known positionally clustered biological programs as negative controls (EDC 1q21, KRT 12q13/17q21, MHC-II 6p21, histones 6p22, metallothioneins 16q13) | Dosage share upper CI < 0.2 (only where allele data exist) | Dosage share lower CI > 0.5 → "dosage-driven"; η² lower CI > 0.3 → "clone-linked". With expression-only CNA, only a "positional flag" is possible |
| **A4 Protocol** | Paired, same-tumour Δlog κ, Δπ (selection) and Δm (induction) between protocols. The null variance comes from same-protocol, multi-piece replicates (regional heterogeneity). κ is a correlation, so it is invariant to gene-wise capture efficiency; this sidesteps InvarMP's κ-versus-δ identifiability problem, and panel retention ι_k ≥ 0.5 is required | \|mean Δlog κ\| ≤ log 1.5 **and** \|Δπ\| ≤ 0.1, by two one-sided tests (TOST) | Lower CI of loss > log 2 in fixed or in situ data → ex vivo; a significant Δπ with stable Δm → selection |
| **A5 QC entanglement** | Correlation of activity with %mito, complexity and %ribo; Shapley share under QC terms | Reported only | Flag only; never used for deletion |

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

- **Decoy audit.** We generate non-parametric decoys by stratified conditional residual permutation: gene-wise permutation of NB randomised-quantile residuals within clone × depth × doublet strata, which keeps each gene's real marginal misfit. These decoys are run label-blind through *unmodified* pipelines (Gavish code, GeneNMF, and the TESSERA engine without certification).
  - The output is a per-MP **decoy recruitment rate**, a measure of how much nuisance can manufacture.
  - It is valid as an audit because the audited pipelines do not condition on the nuisance model.
  - It carries no FDR claim.
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

All are interpretable, frozen and pre-registered on OSF **before** the atlas run.

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
| **GBCD** (Nat Genet 2025) | "GBCD already separates shared from patient-specific programs." | Agreed: we do not claim that split. GBCD has no existence or recurrence test, no tumour membership, no nuisance typing and no projector, and runs per cancer type. It is a TESSERA engine. We certify its output and benchmark GBCD *with* environment covariates (covariate-moderated EBMF; Stephens lab arXiv 2505.11639, referee-reported). |
| **ITHmapper** (bioRxiv 2026) | "There is already a pan-cancer ITH package with a mapper." | Its states are discrete; sample-as-batch can erase tumour-wide signal; it has no uncertainty or provenance. We wrap scVI+Hotspot as an engine, certify its modules, and compare head-to-head on the matched cohort. |
| **TCAT/starCAT** (Nat Methods 2025) | "A reference plus projector is TCAT." | The reference is a by-product. TCAT's reproducibility is cross-dataset recurrence, which systematic artefacts pass; its authors warn about Flex. Our projector adds SEs, novelty detection, panel retention and protocol typing. |
| **Grabski multi-study NMF / Tetris / APAFA / SUFA** | "Subset sharing is solved." | For mutational signatures and bulk data, with MCMC. If ported to single-cell data, these models become engines; our claims do not depend on a correct posterior. A Tetris run on a 30-tumour subsample is held in reserve. |
| **scDEF, consensus scHPF** | "There are already count-native hierarchical models." | They are engines without existence, membership or provenance inference. Benchmarked. |
| **cNMF / GeneNMF / Gavish code** | "Why not just fix the thresholds?" | They are wrapped as engines; the decoy audit and TESSERA certification are applied to their outputs. |
| **Spectra / DRVI** | "Prior-guided or disentangled factors?" | Neither has a recurrence or provenance inference layer. DRVI is benchmarked as the representative deep model. |
| **CanSig benchmark** | "A benchmark already exists." | It has no planted confounders, null cohorts, calibration or typing metrics. We reuse its 12 datasets. |
| **Clipper / ClusterDE / sc-SHC / recall** (referee-reported) | "Target-versus-null contrast and model-based nulls are known." | Correct, and we cite them. We do *not* use synthetic nulls for inference. Decoys are an audit diagnostic only. The inferential contribution is cross-fitted, effect-sized and specific tumour membership, plus real-negative-control provenance. |
| **Count splitting** (Neufeld et al.) / **thinning for Poisson factor models** (Wang et al., JASA, referee-reported) | "Valid inference after any latent-variable estimator already exists." | We use cell splitting for inference, which needs no dispersion assumptions, and thinning only for model selection. No new thinning theory is claimed. |
| **Module preservation (NetRep) / CS-CORE** (referee-reported) | "Audit mode is module preservation; κ is CS-CORE." | Both are baselines. We adopt the CS-CORE-style estimator and add specificity conditioning, within-tumour null bands, cross-tumour recurrence and provenance. Audit mode must beat NetRep Zsummary at flagging labelled artefacts. |
| **ISOMERIC, TreeAlign, CLCC** (referee-reported) | "CN-dependent vs CN-independent expression is done." | Their goal is a single-tumour decomposition. Ours is a typing axis for recurrent programs that is used only with allele anchoring. ISOMERIC is benchmarked where its code runs. |
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
- (c) a trivial typer (MP-score correlation with %mito/%ribo, doublet score, TME-marker cosine, arm enrichment).

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
| Specificity | False membership of non-expressing tumours that carry an overlapping program | ≤5% (versus 72% at 5/40 and 100% at 10/40 shared genes in a referee's toy simulation of a point-null cross-fold statistic) |
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
   - New matched cohort (§7) plus public twins: Janesick (1 block), Oliveira (Flex n = 8, Xenium n = 4), Trinks lung (n = 3), the ccRCC fresh/sn/snFlex preprint, the multi-platform COAD/HCC/OV and breast/NSCLC Xenium sets (referee-reported).
   - Primary endpoint: AUROC of the "ex vivo/protocol-sensitive" calls for predicting loss in fixed or in situ data from the same tumour, against a multi-region null.
   - Secondary endpoint: concordance of B1-typed MP usage across protocols.
4. **Spatial.**
   - Within-malignant-cell co-expression on Xenium 5K / CosMx 6K, with segmentation-spillover controls.
   - Pre-registered niches (hypoxia near necrosis; IFN near T cells), controlling for clone domains (CalicoST).
5. **Clone-orthogonal.**
   - Allele-aware Numbat on ≥6–8 BAM-accessible datasets (≥60 tumours) plus clonealign scDNA sets.
   - Target: agreement of A3 calls between allele-based and expression-based clones, Cohen's κ ≥ 0.6. The power calculation is pre-specified.
6. **Held-out cohort.**
   - Discover on 3CA; replicate on TCCA after removing overlapping studies. Independent patient studies are counted at P0.
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
- TCCA, held out for replication.

**Results:**

1. **A provenance card for every Gavish-41 and 3CA v2 malignant MP**, certified on studies *not* used to derive each list. It reports r_L, lineage prevalence, type, Shapley shares, stability and core genes. Expectations:
   - cycle, IFN/MHC-II and hypoxia are B1 or B2;
   - LQ and contamination MPs are T;
   - **at least 2 retained MPs change status.** This is the go/no-go result from P1.
2. **Per-tumour decomposition of "stress"**, the most contested MP. This is quantitative and non-obvious. Prior work (O'Flanagan 2019; Baron 2020; Neuschulz 2023, all referee-reported) has already established *that* dissociation stress exists. Our results are new in three ways:
   - (a) the ex vivo share of the stress MP's usage **in each tumour**, compared with the TME-cell sensor;
   - (b) whether a protocol-invariant HSP/UPR component is B1 (malignant-intrinsic) or B2 (TME-shared);
   - (c) which *other* canonical MPs carry measurable ex vivo induction. Candidate: hypoxia induced by warm ischaemia.
3. **Selection versus induction.** Prediction: pEMT/mesenchymal-like prevalence is *higher* in Flex and single-nucleus data than in fresh dissociation of the same tumour. That is selection against fragile cells, typed P-selection; a gene-list pipeline would read it as absent.
4. **Malignant-intrinsic versus TME-shared programs.** For each MP, the fraction of its tumour-level variance that co-varies with the same gene set in co-resident TME cells. This separates niche-driven programs (IFN, hypoxia) from cell-intrinsic ones. Validated by spatial and CCLE recurrence.
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

## 8. Go/no-go pilots (months 0–3) and the publishable floor

| Pilot | Duration | Decision criterion |
|---|---|---|
| **P0 Data and design census** | 2 weeks | GO if ≥1,000 tumours with raw UMIs **and** TME cells across ≥8 types, and ≥15% of candidate MPs are used by tumours that bridge ≥2 protocol levels. Also count: empty droplets; BAM access (submit DAC applications in week 1); TCCA overlap; the cohort × lineage × protocol rank. If the bridge share is <15%, A4 depends on the new cohort, so a wet-lab partner becomes mandatory. |
| **P1 Audit of the standard against real nuisance** (key gate) | 4 weeks | On retained 3CA MPs: leave-study-out instability of the Gavish code; SCRP decoy recruitment through the unmodified code (with a Sim-I backbone cross-check); A1 TME co-variation; protocol shift in public twins against the multi-region null. **GO** if ≥3 retained MPs show recruitment ≥30%, *or* A1 partial R² ≥0.3, *or* a protocol loss beyond the null. **PIVOT** if none: re-centre on calibrated membership, gene tiers and projection, and change venue (Genome Biology / NAR GAB). |
| **P2 Specificity and depth invariance of κ** | 3 weeks | Overlap leakage ≤5% at 10–25% shared genes; membership slope ≤0.05 per decade; power ≥0.8 for G2/M and hypoxia in tumours with ≥200 cells. Otherwise revise the null band or margins. |
| **P3 Cross-fitting calibration** | 2 weeks | Type I in [0.03, 0.07] under non-NB marginals on Sim-I nulls; gene FDR ≤0.07 at nominal 0.05. |
| **P4 Typing beats trivial** | 1 week | Evidence-vector AUROC exceeds the trivial typer by ≥0.05 on the expanded labelled set. Otherwise typing is demoted from the headline. |
| **P5 Positional/CNA audit** | 3 weeks (on BAM subset) | If <5% of MPs are dosage-affected, A3 becomes an annotation rather than a headline. |
| **P6 TME sensor** | 2 weeks | Per-sample malignant stress usage against the TME IEG level: partial R² ≥0.3 for stress, ≤0.05 for cycle and hypoxia. This validates A1 specificity. |
| **P7 Engine prototype and scaling** | 6 weeks (overlaps) | Matched cosine for subset-shared programs ≥ GBCD + 0.05; a 3CA-v1-scale fit in ≤24 GPU-hours. Otherwise the default engine is TESSERA-wrapped cNMF, with the headline unchanged. |
| **P8 Scoop check** | 1 week, then every 2 weeks | Landscape §6.3 search terms plus: "cross-fitted gene program", "metaprogram provenance", "protocol-invariant program", "decoy metaprogram", and papers citing GBCD, TCAT, Gavish, Neuschulz and Marsh. NO-GO only if a preprint already types recurrent MPs against within-sample negative controls. |

**Decision at month 3.**
- P0, P1 and P3 pass → full plan.
- P1 fails → pivot.
- P3 fails after one revision → publish the audit plus the benchmark and stop method development.

**Publishable even if the full method fails:**
1. An audit paper: "How many pan-cancer metaprograms survive nuisance-matched decoys, within-sample negative controls and protocol contrasts?" Ready by month 10. Venue: *Nature Methods* Analysis, *Genome Biology* or *Cell Systems*.
2. The Sim-I/Sim-R benchmark as a resource.
3. The matched multi-protocol tumour cohort as a data resource (*Scientific Data* or *Genome Biology*).

---

## 9. Timeline, team, compute and preprint strategy

| Months | Milestone |
|---|---|
| 0–3 | Pilots P0–P8; DAC applications; IRB approval and wet-lab protocol optimisation; exact reproduction of the Gavish code |
| 3–9 | v0.1: split, certify, genes, audit mode, decoy audit; Sim-I/Sim-R v1; **audit manuscript ready internally at M10** |
| 6–16 | Matched cohort: 24 tumours collected at about 2–3 per month, then sequencing |
| 9–15 | TesseraNMF v1; full baseline harness; benchmark on the fractional design; sensitivity tables |
| 14 | **OSF pre-registration**: frozen defaults, margins, type rules, all orthogonal endpoints |
| 15–20 | Atlas run; TCCA; allele subset; Tahoe; spatial |
| 18–22 | Scoring of matched-cohort predictions (pre-registered) |
| 22–26 | Reference v1.0, package v1.0, manuscript; **submission about M24–26**; revision to about M32 |

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
- 50–100k CPU core-hours (the Gavish R pipeline on real plus decoy data; GBCD; Numbat);
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
| 9 | Point-null membership is non-specific and saturates with n (CALIPER F1; CERT-MP F2) | κ is a depth-independent latent correlation, specific given all other programs, tested against a within-tumour matched null band with a floor. P2 gates leakage ≤5% and depth invariance. |
| 10 | Typing by non-rejection; order-dependent ladder (CALIPER F2; InvarMP rigor F1) | TOST equivalence per axis; Shapley order-invariant attribution; "indeterminate" outcome; B1 only on positive evidence. |
| 11 | False or non-novel thinning theorems; Theorem 2 vacuous (CALIPER F3) | No thinning theory is claimed. Thinning is used only for model selection. Neufeld et al. and Hivert et al. are cited. |
| 12 | Engine-indexed estimand contradicts "MP as estimand"; δ degenerate; s_min contradiction (CALIPER F4) | Location = engine proposal; the certificate is engine-free; engine invariance is reported empirically. Fixed biological match tolerance, self-similarity calibrated on independent cell folds. Single-cohort programs form a separate type, not a subtype of "exists". |
| 13 | Fold C used twice (CALIPER rigor F3) | Three disjoint cell folds: discovery, membership, genes. |
| 14 | Transport κ absorbs the dissociation artefact through responsive TME calibrators; κ vs δ not identified for malignant-only genes (InvarMP F1/F3) | No gene-wise transport is used for invariance tests. κ (a correlation) is invariant to gene-wise capture efficiency. TME cells serve as *negative controls* (A1), not as calibrators, and their responsiveness is the signal being measured. |
| 15 | e-values invalid or powerless; e-BH arithmetic (InvarMP F3) | No universal-inference e-values. Paired tumour-level tests on a small, pre-registered hypothesis family (≤10–30 MP × axis pairs) with a stated power. |
| 16 | Pan-cancer invariance claim untestable on public data (InvarMP F2/F4) | The headline is restricted to bridged axes plus the new cohort. Elsewhere the label is "U-no-bridge", reported honestly in the design-graph report. |
| 17 | Thinning null twins understate the real null (regional ITH) (InvarMP rigor) | The protocol null comes from same-protocol multi-region pieces (8 tumours) plus public multi-region data. |
| 18 | QC covariates are post-treatment and absorb biology (CERT-MP novelty) | QC is never an offset. It is a reported A5 axis with Shapley share. |
| 19 | Doublet prior cannot represent misclassified whole TME cells; mimicry programs mis-typed (CREDO rigor; CERT-MP) | Unconstrained usage of TME components. A2 requires a cell-level doublet or ambient association, not loading similarity. Planted mimicry tests in Sim-R. |
| 20 | Underpowered validation (11-positive AUROC; n ≈ 6 matched tumours) | ≥40 labelled artefacts with CIs and a trivial baseline; matched cohort n = 24 with a power calculation; non-canonical Tahoe endpoints. |
| 21 | Validation is non-comparative | Every rung has a pre-specified endpoint, compared head-to-head with Gavish/3CA, GBCD and ITHmapper MPs. |
| 22 | Scope and team about 2× under-resourced (all) | Scope cut: no SVI or IBP model, no power posterior, no normal-lineage nodes, no γ_g resource, no 17-method harness. About 4 FTE; about 26 months to submission. |
| 23 | Data availability (raw UMIs, empty droplets, BAMs) overestimated | P0 gate; DAC applications at M0; ambient proxy flagged; analyses reported with and without each optional module. |
| 24 | The flagship confirms known biology | P1 gate; flagship built on per-tumour quantitative provenance and non-obvious selection/induction predictions. |
| 25 | Wet-lab cohort fails (accrual, cost) | Minimal cohort of 12 tumours; PDX substitute; public twins; a computational-only claim limited to bridged axes. |
| 26 | Scoop | Engine-agnostic claim; tripwire preprint at M10; fortnightly monitoring. |

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
7. **Pre-registration on OSF at M14.** Recommended; it constrains the analysis but is the best defence against cherry-picking critiques.
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
- Kotliar…Raychaudhuri, TCAT/starCAT, Nat Methods 22:1964–1980 (2025) [V]
- Nofech-Mozes…Abelson, ITHmapper, bioRxiv 2026 [V]
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
- TCCA, bioRxiv 2026 [V]
- Tahoe-100M [V]
- Numbat, CopyKAT, SCEVAN [V]
- CellBender, SoupX, DecontX [V]
- scDblFinder [V]
- clonealign [V]
- CalicoST [V+; journal unconfirmed]
- Janesick 2023 [V; venue P]
- Oliveira, Nat Genet 2025 [V]
- Splatter, splatPop, scDesign3 [V]
- Venet et al., PLoS CB 2011 [V]
- Hotspot [V]
- MOFA+ [V]
- Replogle 2022 [P]
- Greenwald GBM spatial [V; venue P]

**B. Cited by referees as verified in their own sessions, not in the landscape and not re-checked by the PI. Verify before writing:**
- CS-CORE, Nat Commun 14:4846 (2023)
- Clipper, Genome Biol 2021
- sc-SHC, Grabski, Street, Irizarry, Nat Methods 20:1196–1202 (2023)
- Langfelder et al., PLoS CB 2011; NetRep, Cell Systems 2016
- Wang, Xu, Zhao, Wang, data thinning for Poisson factor models, JASA (2025/26)
- Nomura et al. / Tirosh lab, Nat Genet 57:1155–1167 (2025) (landscape listed it as [P])
- ISOMERIC, bioRxiv 10.1101/2025.04.21.649842
- CLCC, bioRxiv 10.64898/2026.03.01.708901
- subclonal-dosage preprint, bioRxiv 10.64898/2026.03.15.711906
- TreeAlign (venue unverified)
- Neuschulz et al., MSB 19:e11147 (2023)
- Marsh et al., Nat Neurosci 25:306–316 (2022)
- O'Flanagan et al., Genome Biol 2019 (landscape listed it as [U])
- Baron et al., Cell Systems 11:536 (2020)
- FixNCut, Genome Biol 25:81 (2024)
- Trinks et al. 2024 (3 LUAD patients)
- ccRCC protocol benchmark, bioRxiv 2025.11.25.690485
- inVAE, bioRxiv 2024
- sysVI, BMC Genomics 2025
- covariate-moderated EBMF, arXiv 2505.11639
- Weine et al., log1p Poisson NMF, arXiv 2601.05845
- Miller & Dunson, coarsened posteriors, JASA 2019
- Bogomolov & Heller, Stat Sci 2023
- Emery et al., arXiv 1907.01458
- Sarkar & Stephens, Nat Genet 2021
- Gene Program Dictionary, bioRxiv 2025.10.29.685322
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
- dissociation papers van den Brink 2017, Denisenko 2020, Slyper 2020
- CytoTRACE

The scoop search of 2026-10-02 was proxy-limited, so the absence of a direct precedent is weak evidence. The manual check in P8 is mandatory.