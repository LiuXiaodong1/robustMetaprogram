# Design panel: four candidate concepts and their referee reports

> AI-assisted design exercise (2026-10-02). Four independent concepts were generated from different design angles, and each was critiqued by two simulated *Nature Methods* referees: novelty/prior art, and statistical rigor/feasibility. The final proposal ([02_TESSERA_proposal.md](02_TESSERA_proposal.md)) is a hybrid that responds to every fatal flaw listed here. Citations named by referees were **not** independently re-verified.
>
> **Erratum (2026-10-17).** This file is a frozen record of the 2026-10-02 design panel and was not revised after the 2026-10-16 review. Corrected citations are in [01](01_landscape_and_gap_analysis.md) (revision note, §2B, §4.1, §6.1) and [02](02_TESSERA_proposal.md) (§4, §6.4, §12 A2). In particular: (1) the statements here that TCAT/starCAT covers "T cells only" (lines 576 and 3133) are superseded: the TCAT catalogue is T-cell, but its authors present starCAT as generalising to other cell types and tissues; (2) Oliveira et al. 2025 profiles FFPE blocks only (Chromium Flex, Visium HD and Xenium on the same blocks; no fresh-dissociated malignant cells were found), so where the panel counts it among matched-protocol tumours, read it as a fixed-dissociated versus in situ comparison, not a fresh-versus-fixed twin; (3) the figure of 55% at fold change 1.5 in Referee 2's CALIPER report (line 1455) has no saved output behind it and must not be cited ([referee_checks/README.md](referee_checks/README.md)).

## Score summary (1–10)

| Angle | Concept | Referee | Novelty | Soundness | Feasibility | Impact | NM fit |
|---|---|---|---|---|---|---|---|
| generative-model | CREDO | novelty | 5 | 4 | 3 | 6 | 6 |
| generative-model | CREDO | rigor | 6 | 4 | 3 | 6 | 7 |
| inference-layer | CALIPER | novelty | 4 | 3.5 | 4 | 5 | 5 |
| inference-layer | CALIPER | rigor | 6 | 4 | 4 | 6 | 6 |
| invariance | InvarMP | novelty | 4 | 3 | 3 | 6 | 6 |
| invariance | InvarMP | rigor | 7 | 3.5 | 3 | 7 | 7 |
| feasibility-first | CERT-MP | novelty | 4 | 3 | 4 | 5 | 5 |
| feasibility-first | CERT-MP | rigor | 6 | 3 | 4 | 6 | 6 |


---

## 1. CREDO: Calibrated REcurrence Decomposition with Offsets. The name has not been checked for collisions; fallback name is "mpGen".  
*Design angle: `generative-model`*

### Claim

CREDO turns tumour metaprograms from threshold-defined top-50 gene lists into a calibrated generative object. A metaprogram becomes a non-negative program used by a credible subset of tumours. It is estimated jointly from the raw counts of about 3,000 tumours and about 10^7 malignant cells, with copy-number dosage, ambient RNA, doublets, cohort/protocol effects and per-tumour baselines modelled explicitly. Each MP comes with an error-controlled existence call, tumour membership and gene set, plus an automatic typing (biological / clone-linked / protocol-sensitive / nuisance / cohort-restricted) that replaces manual curation.

### Formal definition of an MP

SETTING. Malignant cells n are nested in clones h, clones in tumours t (a tumour can have several samples s on different protocols p), tumours in cohorts c (one study x one protocol), and cohorts in lineages l (cancer types; optionally normal-tissue donor nodes). Data are raw UMI counts y_ng for a gene universe G (about 6,000 genes).

DEFINITION. Under the CREDO generative model M (see model_and_math), a metaprogram is a pair (W_k, T_k):
- W_k is a non-negative gene loading vector on a fold-change scale relative to each gene's pooled baseline, normalised so that sum_g wbar_g W_kg = 1.
- T_k = {t : z_tk = 1} is the set of tumours in which program k makes a NON-TRIVIAL contribution to malignant-cell expression. Non-trivial means a tumour-mean usage m_tk >= delta_min. The default delta_min means k accounts for at least 0.5% of the tumour's expected biological UMIs. The contribution is measured AFTER the model accounts for CNA dosage, ambient RNA, doublets, cohort/protocol gene efficiency, QC covariates and the tumour's own baseline.

A slot k is reported as an MP only if all three conditions hold:
- (a) Recurrence: |T_k| >= 2 with members from >= 2 cohorts. Otherwise it belongs to the tumour-private or cohort-private tier by construction, not by threshold.
- (b) Existence: q_k <= alpha (default 0.1). q_k is a target-decoy q-value for H0_k: "the held-out deviance explained by k is no larger in real tumours than in nuisance-preserving decoy tumours carried through the same fit".
- (c) Effect size above delta_min in at least two members.

Every MP is reported with:
- tumour-membership posterior inclusion probabilities PIP_tk, and a credible tumour set at a stated membership FDR;
- a variable-length gene set at a stated gene FDR, from count-split held-out tests;
- a lineage profile pi_{k,l} with credible intervals;
- an ITH index (within-tumour share of usage variance), which separates within-tumour ITH programs from tumour-wide subtype programs;
- a stability index and an identifiability margin;
- a type label.

WHAT CHANGES. The only user-facing constants are an effect size (delta_min) and error rates (alpha, gene FDR, membership FDR), all interpretable. These replace the overlap thresholds 35/10/10/5 and greedy clustering. Recurrence is necessary but not sufficient: a program must also beat decoys that keep every measured nuisance and the empirical marginal noise. The decoys do not preserve unmeasured nuisance co-expression, such as dissociation stress. That is handled by the protocol-sensitivity typing, which is reported rather than hidden.

### Model and math

(1) TEMPERED, TUMOUR-BALANCED NB LIKELIHOOD
y_ng ~ NB(mean mu_ng, inverse-dispersion phi_{g,p(n)}). phi is per gene and protocol, so Smart-seq2 read counts, sn and Flex get their own dispersions.

Generalised (power) posterior:
log p_zeta(Theta | Y) = zeta * sum_t w_t * sum_{n in t} sum_g log NB(y_ng | mu_ng, phi) + log p(Theta) + const
- w_t = min(1, n0 / n_t), with n0 = 2,000. This stops 20k-cell tumours from dictating the shared dictionary and caps each tumour's effective sample size.
- zeta in (0, 1] is the power-posterior exponent. In the BayesPowerNMF style it stops systematic misfit (overdispersion tails, nonlinear ambient leakage, dropout) from spawning extra components, which would otherwise recur across tumours because the misfit is systematic.

(2) MEAN STRUCTURE: biology on an additive non-negative rate scale, nuisance as fixed or anchored terms
mu_ng = s_n * [ (1 - rho_n - delta_n) * nu_{c(n),g} * D_ng * lambda_ng + rho_n * A_{s(n),g} + delta_n * sum_e pi_ne * M_{s(n),e,g} ]

- s_n: the observed UMI total times exp(a small learned jitter).
- Ambient: A_s is fixed per sample, estimated from empty droplets (CellBender or SoupX). rho_n ~ Beta with mean rho_hat_s (the sample-level estimate) and concentration 50.
- Heterotypic doublets / contamination: delta_n = b_n * omega_n.
  - b_n ~ Bernoulli(p_hat_n), where p_hat_n is the scDblFinder score after isotonic calibration.
  - omega_n ~ Beta(5, 5).
  - M_{s,e} is the normalised mean profile of TME type e (macrophage, T/NK, B/plasma, endothelial, fibroblast) in the SAME sample, falling back to the cohort. pi_n ~ Dirichlet.
  - This is the term that should absorb the Macro / T-cell / Endo MPs that 3CA v2 removed by hand.
- Cohort/protocol gene efficiency: log nu_{c,g} = a_{p(c),g} + e_{c,g}, with e ~ N(0, sigma_e^2) and a normalisation constraint.
  - a_{p,g} (sc vs sn vs Flex vs Smart-seq2) is identified from tumours profiled under >= 2 protocols.
- CNA dosage offset: D_ng = sum_h q_n(h) * (CN_{t,h,seg(g)} / P_t)^{gamma_g}.
  - q_n(h) is the clone posterior from Numbat (allele-aware) run on count-thinning fold A only. Where BAMs are unavailable, an expression-only caller is used on fold A.
  - P_t is tumour ploidy.
  - gamma_g in [0, 1] ~ Beta(2, 2) is shared across ALL tumours: a pan-cancer per-gene dosage sensitivity, which is also a by-product resource.
  - Genetics enters by dosage only. Clone-linked non-dosage biology stays in programs and is TYPED afterwards, never deleted.
- Biological rate:
  lambda_ng = wbar_g * [ beta_{t,g} + sum_k z_tk theta_nk W_kg + sum_r x_tr thetaU_nr U_{c,r,g} + sum_j thetaV_nj V_{t,j,g} ] + sum_q psi_nq Q_qg
  - wbar_g is the gene's pooled mean, so loadings are fold changes. W_kg = 0 means "no change". Highly expressed housekeeping, ribosomal and mito genes no longer dominate (fixes S1).
  - Sparse loadings: W_kg ~ Gamma(0.3, 0.3), normalised as sum_g wbar_g W_kg = 1.
  - PATIENT BASELINE LAYER (no centring): beta_{t,g} = b_{l(t),g} * eps_{t,g}, with eps ~ Gamma(tau, tau) i.i.d. over genes. Because eps is idiosyncratic, the tumour level is factor-analytic. A structured shift shared across many tumours costs less as a shared program (m_tk W_k) than as G independent deviations. Tumour-wide subtype programs therefore surface as MPs instead of being centred away.
  - Usage: theta_nk = m_tk * h_nk, with h_nk ~ Gamma(alpha_tk, alpha_tk) (mean 1).
    - m_tk = z_tk * mtilde_tk, with mtilde_tk ~ Gamma(a_m, a_m / mu_tk) truncated below at delta_min.
    - log mu_ts,k = eta_k0 + eta_{k,l(t)} + eta_{k,p(s)}, where eta_{k,p} is protocol sensitivity.
    - A large alpha_tk means a tumour-wide (homogeneous) program.
  - SUBSET SHARING (hierarchical beta-Bernoulli; a finite approximation to a hierarchical beta process, Thibaux and Jordan 2007 [U, from memory]; cf. James, Lee and Pandey arXiv 2304.05244 [V title/authors only]):
    - z_tk ~ Bernoulli(pi_{k,l(t)})
    - pi_{k,l} ~ Beta(c * pi_k, c * (1 - pi_k))
    - pi_k ~ Beta(alpha0 / K, 1), with K = 200 slots.
    - This yields pan-cancer, lineage-restricted and rare MPs, plus a normal-shared probability when normal donors are lineage nodes.
  - COHORT TIER: U_{c,r} with R = 4 slots per cohort and inclusion x_tr ~ Bernoulli(piU_{c,r}). This absorbs study-specific artefacts such as a lab's dissociation protocol.
  - PRIVATE TIER: V_{t,j} with J = 4 slots per tumour. Usage has an ARD scale tau_tj ~ Gamma(0.01, 0.01), so unused slots switch off.
  - QUALITY TIER: psi_nq = softplus(kappa_q^T c_n).
    - c_n holds observed QC covariates: log UMI, genes/UMI, %mito, %ribo, nuclear-lncRNA fraction.
    - Because usage is confined to the span of the covariates, these terms are identifiable as nuisance.
  - CONE-EXCLUSION PRIOR: log p += -lambda_ce * [ sum_{t,j,k} z_tk cos^2(V_tj, W_k) + sum_{c,r,k} cos^2(U_cr, W_k) ].
    - Private and cohort programs may not re-express shared ones.
    - This is the opposite of SUFA's subspace constraint, which is biologically wrong for amplicon-private programs.

(3) DERIVED QUANTITIES
- ITH_k = E[within-tumour variance of theta_nk] / E[total variance].
- Recurrence R_k = sum_t PIP_tk.
- Clone linkage: CL_k = the share of within-tumour usage variance explained by clone, averaged over T_k.
- Protocol sensitivity: eta_{k,p}.
- Lineage profile: pi_{k,l}.
- Positional enrichment: a chromosome-arm permutation test on W_k that preserves gene expression level.

(4) MEMBERSHIP: LOCAL TEMPERED BAYES FACTORS
- LBF_tk = zeta * w_t * [ELBO_t(z_tk = 1) - ELBO_t(z_tk = 0)], computed with tumour-t usages re-optimised and globals fixed.
- PIP_tk = sigmoid(LBF_tk + logit E[pi_{k,l(t)}]).
- Tempering stops a 10,000-cell tumour from returning astronomically large Bayes factors for trivial structure.
- The delta_min truncation makes membership mean meaningful presence.

(5) EXISTENCE: TARGET-DECOY CALIBRATION INSIDE ONE FIT (the target-decoy idea comes from proteomics, Elias and Gygi 2007 [U])
- Decoys: 15% of tumours are duplicated as decoys using STRATIFIED CONDITIONAL RESIDUAL PERMUTATION (SCRP):
  - fit a nuisance-only NB GLM per tumour (covariates s_n, D_n, rho_hat, doublet score, QC);
  - compute randomised quantile residuals u_ng [Dunn and Smyth, U];
  - permute u_.g independently per gene within strata (clone x depth decile x doublet bin);
  - invert back to counts.
- What decoys keep and destroy:
  - They keep each gene's empirical conditional marginal, including zero inflation and heavy tails (so they carry real misspecification).
  - They keep all covariance induced by measured nuisance and the clone structure.
  - They destroy program co-expression.
  - The model never sees the decoy labels.
- For slot k: Delta_tk is the fractional reduction of held-out (fold B) NB deviance from including k in tumour t.
- Statistic: T_k = mean_{t in real} PIP_tk * Delta_tk - mean_{t in decoy} PIP_tk * Delta_tk.
- p_k comes from permuting real/decoy labels at tumour level, stratified by cohort. BH across slots gives q_k.
- Membership FDR for a credible set at threshold thr: FDRhat_k(thr) = [#decoys with PIP >= thr / n_dec] * n_real / #[real with PIP >= thr].
- Secondary global check: 10 full nuisance-only null cohorts, built entirely from SCRP decoys, must yield about 0 MPs at q <= 0.1.

(6) GENE MEMBERSHIP BY NB COUNT SPLITTING
- Thinning: X^A | X ~ BetaBinomial(X, eps * phi, (1 - eps) * phi), with eps = 0.5. This gives X^A ~ NB(eps * mu, eps * phi) and X^B ~ NB((1 - eps) * mu, (1 - eps) * phi), independent (generalised data thinning).
- phi is estimated CONDITIONALLY on a pilot low-rank fit, not marginally. This answers Hivert et al.'s mixture objection.
- Programs, CNA calls and usages are learned on fold A.
- On fold B, for each MP k, gene g and member tumour, a one-sided score test of beta_gk >= 0 in mu^B_ng = (1 - eps) * muhat^(-k)_ng + beta_gk * s_n * thetahat^A_nk * nu * D.
- Tumour-level tests are combined by PIP-weighted Stouffer, then BH within k. The result is FDR-controlled gene sets of any length plus out-of-fold effect sizes.
- 500 decoy genes (real genes permuted within tumour) must give uniform p-values. If they do not, an empirical-null correction estimated from the decoy genes is applied.

(7) CHOOSING ZETA AND K
- Grid: zeta in {1, 0.5, 0.2, 0.1, 0.05, 0.02}, each fitted on a stratified 500-tumour subsample.
- zeta* is the largest zeta with decoy membership FDRhat <= 0.1 AND held-out NB deviance within 1 SE of the minimum.
- K is not chosen by the user. Slots die through the beta-Bernoulli prior and through tier moves, and existence is then decided by q_k.
- A curve of MP count against zeta is always reported.

(8) IDENTIFIABILITY (stated and tested)
- Scale and permutation are fixed by the simplex normalisation and Hungarian matching across seeds.
- Proposition 1, pooling restores NMF identifiability:
  - the sufficiently scattered condition (SSC) is monotone under row union, because cone(Theta_A^T) is contained in cone([Theta_A; Theta_B]^T);
  - so if the stacked usage of the members of T_k satisfies SSC, W is identifiable even when no single tumour does;
  - per-tumour NMF needs SSC inside every tumour, which is strictly stronger.
- Proposition 2, shared vs private separation (noiseless; to be proven). Assume:
  - (C1) each W_k has >= 1 anchor gene that is exclusive among global programs and carries bounded private loading;
  - (C2) the members of T_k jointly satisfy SSC;
  - (C3) |T_k| >= 2, with private cones intersecting only at 0 and excluding W_k.
  Then any alternative that moves W_k into private slots needs >= 2 duplicate private copies, which raises the l0 complexity. Under the sparsity and cone-exclusion priors the decomposition is therefore unique up to permutation and scale.
- Nuisance is identifiable because its loadings are fixed (A_s, M_{s,e}, CN) or its usage is confined to the covariate span.
- Identifiability is weakest where tumour equals library. Tumour-wide programs are identified only through co-variation across many tumours (eps is i.i.d. over genes). The method requires >= 3 tumours per cohort and uses multi-protocol tumours to anchor nu.
- Reported per MP: number of anchor genes, number of near-pure cells (share of theta_nk >= 0.8), and an identifiability margin (the smallest principal angle between W_k and the cone of all other programs).

### Inference and scalability

ENGINE
Stochastic variational inference in Pyro/PyTorch on GPU, optimising the generalised ELBO:
L = zeta * sum_t w_t * E_q[log p(Y_t | Theta)] - KL(q || p).

VARIATIONAL FAMILY
- Global loadings (W, U, V, eps, b, gamma, a, e, phi, eta): log-normal mean-field.
- z and x: Hard-Concrete relaxed Bernoullis [Louizos et al. 2018, U], with temperature annealed from 1 to 0.1. They are replaced by exact local Bayes factors at the end.
- Cell-local h_n: amortised by an encoder (MLP 6,000 -> 512 -> 256 -> 2 x (K + R + J)) whose input is log1p-normalised counts plus the tumour's slot mask. Outputs are masked by z_t and x_t. The encoder never sees tumour identity, so it transfers to projection.
- Tumour-local parameters (eps_t, V_t, m_t, alpha_t, z_t) live in embedding tables updated with sparse Adam.
- VI is used for ESTIMATION only. Inference claims (membership, existence, genes) are calibrated by decoys and count splitting, not by mean-field variances, which are known to be overconfident.

MINIBATCHING
Tumour-stratified: 64 tumours (sampled in proportion to w_t * n_t) x 256 cells = 16,384 cells per step, with unbiased likelihood scaling.

TIER MOVES (how the IBP-style structure is explored under SVI)
Every 5 epochs, birth, death and reassignment proposals are evaluated. Each is accepted only if the tempered ELBO on held-out cells improves.
- A global slot with R_k < 1.5 becomes private to its dominant tumour.
- Private programs from >= 3 tumours in >= 2 cohorts with cosine > 0.7 are promoted to an empty global slot.
- A global slot with > 90% of its members in one cohort moves to the cohort tier.
- Two global slots with cosine > 0.9 and high co-membership are merged.
- A slot is split when its per-tumour loading residuals are bimodal.

INITIALISATION
- 50 k-means metacells per tumour (about 150k in total).
- GPU KL-NMF with K = 200 initialises W.
- Per-tumour residual NMF initialises V.
- 5 seeds; agreement across seeds is reported as a stability measure.

POLISHING
- With globals fixed, exact multiplicative NB/Poisson updates of every cell's usages, parallel across tumours.
- Then 20 local gradient steps per cell (semi-amortised) to close the amortisation gap.
- Then local Bayes factors (LBF) for z.

COST AT ATLAS SCALE (3,000 tumours, 10^7 cells, G = 6,000, about 220 active components per cell)
- Forward pass: about 2 x 10^10 multiply-accumulates per step, an estimated 10–30 ms per step on one H100 including the encoder. One epoch is about 610 steps.
- Data: about 1.5 x 10^10 non-zeros, roughly 100 GB of CSR shards per tumour on NVMe, streamed from backed AnnData or TileDB-SOMA. I/O, not FLOPs, is the bottleneck.
- Parameters: about 0.2 B floats, dominated by V at 3,000 x 4 x 6,000 = 72 M. This is under 4 GB including optimiser state.
- Estimated wall time: 6–24 H100-hours per full fit, to be verified in pilot P5.
- Multi-GPU: data-parallel over tumour shards, with globals all-reduced and tumour-locals sharded.
- Full calibration (zeta grid, 5 seeds, both thinning folds, leave-one-cohort-out, 10 null-cohort refits on 500-tumour subsamples): about 1,500–4,000 GPU-hours.

FALLBACK
A two-stage version of the same model: learn the dictionary on <= 500 cells per tumour, then polish all cells with W fixed. It is roughly 10x cheaper.

PROJECTION
- Globals stay fixed. Each query tumour gets its own z_t, usages, nuisance terms and J_q = 4 novelty slots.
- Runtime is seconds to minutes per tumour.
- Panel transfer to Xenium/Flex restricts G to the panel and uses the protocol effects a_{p,g} learned from matched data.

### Outputs for users

1. MP TABLE. One row per MP with:
   - id and type label with its evidence vector;
   - existence q-value;
   - recurrence (expected number of tumours, cohorts and lineages);
   - lineage profile pi_{k,l} with 90% credible intervals;
   - ITH index (within-tumour ITH program vs tumour-wide subtype program);
   - protocol sensitivity eta_{k,p};
   - clone linkage CL_k and positional (arm) enrichment;
   - nuisance alignment (maximum cosine to ambient, TME-doublet and QC profiles);
   - decoy recruitment rate;
   - stability index (share of seeds x thinning folds x leave-one-cohort-out refits with a matched program at cosine >= 0.8);
   - identifiability margin and anchor genes.
2. MEMBERSHIP. A tumours x MPs PIP matrix, credible tumour sets at a user-chosen membership FDR (decoy-calibrated), and the posterior co-membership network of MPs (repertoires, co-occurrence and exclusivity).
3. GENES.
   - Fold-change loadings with intervals.
   - FDR-controlled gene sets of variable length, in two tiers: core q < 0.01 and extended q < 0.1.
   - GMT/MSigDB-style export that carries q-values.
4. CELLS. Per-cell usages with posterior SD on a scale comparable across tumours (no centring), stored in obsm; tumour-mean usages; per-tumour private programs and cohort programs, catalogued and labelled as such rather than discarded.
5. NUISANCE DECOMPOSITION.
   - Per-cell ambient fraction, doublet posterior and CNA dosage.
   - Pan-cancer gene dosage sensitivity gamma_g, a stand-alone resource.
   - Protocol gene-efficiency table a_{p,g} (sc / sn / Flex / Smart-seq2).
6. REFERENCE BUNDLE AND PROJECTOR.
   - A versioned pan-cancer reference on Zenodo.
   - project() returns membership PIPs, usages with uncertainty, novelty slots that flag programs the reference lacks, and panel transfer.
7. HTML DIAGNOSTIC REPORT. Zeta curve, decoy recruitment, null-cohort result, decoy-gene p-value calibration, stability heatmaps and identifiability diagnostics.
8. CREDO-Sim. A semi-synthetic tumour simulator with planted confounders, plus a containerised benchmark harness running all baselines.

### Failure modes fixed

S0 MALIGNANT CALLING / CONTAMINATION
Heterotypic doublets and contamination are modelled as a mixture with sample-matched TME profiles (M_{s,e}) whose prior is anchored on scDblFinder scores. Ambient RNA is modelled with a fixed empty-droplet profile. Macrophage, T-cell and endothelial "MPs" (removed by hand in 3CA v2) are absorbed or typed "nuisance-aligned" instead.

S1 TOP-7,000 BY MEAN EXPRESSION, CENTRING AND ZEROING NEGATIVES
- The gene universe is chosen by deviance (Townes 2019 [U]), plus all TFs.
- Loadings are on a fold-change scale (wbar_g), so highly expressed housekeeping, ribosomal and mito genes no longer dominate.
- The NB likelihood replaces log, centre and zero.
- A patient-baseline layer (beta_t, factor-analytic) plus tumour-level usage m_tk makes tumour-wide programs visible and comparable across tumours. The ITH index distinguishes them from within-tumour programs. This addresses failure mode 3 in Section 1.4 (centring blocks joint inference) and gap G5.

S2 FIXED RANK GRID 4–9 AND BEST OF 10 RESTARTS
- The number of programs comes from the beta-Bernoulli prior, tier moves and the power posterior. Zeta is chosen by held-out thinned deviance plus decoy calibration.
- Pooling restores identifiability (Proposition 1), and each MP carries an identifiability diagnostic.
- Restart variance is replaced by a stability index across seeds, folds and cohort hold-outs.

S2 PER-TUMOUR FITS (NO BORROWING)
A single joint fit lets small tumours borrow W. Tumour-balanced weighting w_t stops large tumours from dominating.

S3 TOP-50 TRUNCATION
Full loadings are kept. Gene sets of variable length are FDR-controlled through NB count splitting, with decoy genes checking calibration.

S4 35/10/10 FILTERS
- The filters are replaced by membership PIPs, decoy-calibrated credible tumour sets and target-decoy existence q-values.
- There is no cross-rank stability selection, so the bias toward dominant technical and CNA axes disappears.
- The irrelevant random-gene null (p ~ 1e-12) is replaced by a nuisance-preserving null that carries realistic misspecification.

S5 GREEDY CLUSTERING
There is no post-hoc clustering. MPs are model parameters. Dependence on cohort composition is measured by leave-one-cohort-out refits and by adding or removing 10–50% of studies.

S6 MANUAL CURATION AND POST-HOC CNA
- CNA dosage is an in-model offset, with clone calls on an independent thinning fold, so recurrent arm-level CNAs (e.g., chr7 gain) cannot form positional MPs.
- Clone-linked non-dosage programs are typed, not deleted.
- Pre-registered typing rules replace inds_rm. They are evaluated as a held-out AUROC on the 11 hand-removed 3CA v2 MPs.

S7 SCORES WITHOUT UNCERTAINTY
Per-cell posterior usages on a common scale, and a projector that reports uncertainty and novelty.

NAMED CONCEPTUAL FAILURES
- "Recurrence is not validity": decoys and null cohorts.
- "No likelihood, so no uncertainty": NB plus calibrated inference.
- "Filters are biased": no rank grid or overlap filters.
- "Cohort dependence": stability is reported.
- "Validation only post hoc": CNA, protocol and nuisance are inside discovery, and typing uses them.

### Validation plan

A. SEMI-SYNTHETIC GROUND TRUTH (CREDO-Sim; the primary evidence)
Backbone. About 300 real tumours (the 12 CanSig datasets / 185 patients plus a 3CA subset), turned into program-free backbones by SCRP. This keeps real depth, QC, clone structure and marginal noise.

Planted programs (30 per replicate):
- 6 pan-cancer (prevalence 60–90%);
- 12 lineage-restricted (1–2 of 5 lineages, within-lineage prevalence 30–80%);
- 8 rare (3–10 tumours);
- 4 nested or overlapping pairs (cosine 0.4–0.6);
- both tumour-wide (homogeneous) and heterogeneous versions;
- effect sizes 1–20% of cell UMIs.
Each replicate also gets 2–4 private programs per tumour and 1–2 cohort artefacts per cohort.

Planted confounders, anchored in measured quantities:
- recurrent arm-level CNAs copied from real Numbat profiles (in 30–70% of a lineage's tumours);
- ambient RNA at 0–25%, using real empty-droplet profiles;
- heterotypic doublets at 0–15%, using real TME cells;
- a stress/IEG block in 0–40% of cells, varying by cohort and protocol;
- 3-fold depth variation;
- sc/sn gene-efficiency shifts taken from matched data.

Design: 4 difficulty levels x 10 replicates, plus 20 nuisance-only null cohorts.

Metrics (scorecard from the landscape, Section 4.2):
- Hungarian-matched cosine; F1 at cosine >= 0.8;
- gene AUPRC; empirical versus nominal gene FDR;
- membership AUROC and credible-set FDR calibration (reliability curves);
- false-MP rate (MPs aligned to a planted nuisance at cosine > 0.6);
- private-program leakage;
- number of MPs called on null cohorts;
- error in K;
- runtime and memory.

Pre-registered targets (CREDO):
- median <= 1 MP called on null cohorts at q = 0.1;
- false-MP rate <= 5%;
- F1 >= 0.85 for programs at >= 5% of UMIs;
- membership AUROC >= 0.9, with empirical credible-set FDR within nominal +/- 2 points;
- empirical gene FDR <= nominal + 2 points;
- leakage <= 5%;
- K within +/- 10%.
Hypothesis to be measured, not assumed: the Gavish code calls >= 5 MPs on null cohorts and has >= 20% leakage.

Ablations, each expected to break a specific metric:
- zeta = 1 (more junk MPs);
- no w_t;
- no CNA offset (positional MPs);
- no ambient or doublet terms (contamination MPs);
- no private tier (leakage);
- no cohort tier (cohort artefacts promoted to MPs);
- all-shared, no z (membership lost);
- no decoys (VI PIPs overconfident);
- centred input (tumour-wide programs lost);
- Gaussian-log likelihood (isolates the count likelihood against GBCD).

B. CROSS-PROTOCOL INVARIANCE (the answer to "is it real biology?")
Data: Janesick FFPE breast (Flex / Visium / Xenium) and Oliveira CRC (Visium HD / Flex / Xenium), plus any 3CA studies that profiled the same tumours by sc and sn [P].
Tests:
- MPs typed "biological" should replicate in fixed or nuclear data, with panel-restricted loading correlation and matched usage distributions.
- MPs typed "protocol-sensitive" should collapse.
Metrics: per-MP usage concordance across protocols on the same tumour, and agreement of type labels under leave-one-protocol-out.

C. CAUSAL CHECKS
- Project the reference onto Tahoe-100M and test pre-registered mechanism-to-MP predictions against DMSO controls (CDK4/6 inhibitors lower cell-cycle MPs; MEK inhibitors lower MAPK-response MPs; others fixed before analysis).
- Use PerturbNMF-style specificity metrics.
- Replogle Perturb-seq with non-targeting controls.

D. SPATIAL
- Xenium / Visium HD with CalicoST clone domains [preprint].
- Measure MP spatial autocorrelation conditional on clone domain.
- Pre-registered niches: hypoxia near necrosis in GBM (Greenwald data [P]), interface/EMT programs at the edge, IFN near T cells.

E. GENETIC ORTHOGONALITY AND HERITABILITY
- On DNA- or allele-defined clones (clonealign datasets, Numbat on BAM-available cohorts), MPs typed "clone-linked" should show higher clone-explained variance than plastic MPs.
- KP-Tracer lineage trees (mouse): plastic MPs should have low phylogenetic autocorrelation.

F. GENERALITY AND REPLICATION
- Learn on 3CA and replicate on TCCA (studies shared with 3CA removed): replication rate per MP and recovery of type labels.
- Kinker CCLE (198 lines): cell-intrinsic recurrence. Protocol-sensitive in vivo MPs should be absent.
- TCGA survival associations must beat matched random signatures (Venet).

G. ARTEFACT TYPING
AUROC for flagging the 11 hand-removed 3CA v2 MPs plus Gavish's removed LQ, doublet and single-study MPs. Typing rules are frozen before scoring.

H. STABILITY
Matched cosine and variance of the MP count across seeds, thinning folds, leave-one-cohort-out, and adding or removing 10–50% of studies, compared with the Gavish code under identical perturbations.

BASELINE PANEL (tuned fairly, each in its intended setting)
Gavish/3CA code (not just scalop), GeneNMF, cNMF (per tumour and pooled), TCAT-style Harmony+cNMF, GBCD (per cancer type), consensus iNMF, scINSIGHT, MOFA+/MOFA-FLEX, DRVI, Spectra, Hotspot and ITHmapper, CanSig's top early-integration methods, mosaicMPI, consensus scHPF, scDEF, and Tetris/BMSFA on subsamples.

### Flagship biology

GOAL: a calibrated pan-cancer map of which tumours share which programs, plus an audit of the canonical MPs.

DATA
Joint fit on the 3CA v2 corpus (124 datasets / 2,836 samples [P]) plus other public tumour scRNA-seq with raw UMIs, targeting about 3,000 tumours and 5–10 x 10^6 malignant cells. Normal epithelium donors from public atlases are added as lineage nodes. TCCA is held out for replication.

RESULT 1: AUDIT OF EVERY GAVISH-41 AND 3CA v2 MALIGNANT MP
Each gets a type label, q-value, membership breadth, stability and identifiability. These are hypotheses to test, not results:
- (a) "Stress"-type MPs split into a protocol-sensitive ex vivo IEG/heat-shock component that vanishes in Flex/sn/Xenium, and a protocol-invariant hypoxia/UPR component that is spatially tied to necrosis.
- (b) Some EMT/mesenchymal-like carcinoma MPs are real in some tumours but CAF-doublet or ambient-driven in others. Tumour-level membership resolves this heterogeneity, which a single gene list cannot.
- (c) Some canonical gene lists contain arm-level positional components (e.g., chr7/chr10 in GBM, 1q in breast) that disappear under dosage offsets.
- (d) The pan-cancer dosage-sensitivity map gamma_g identifies dosage-compensated genes.

RESULT 2: PROGRAMS INVISIBLE TO CENTRING
Tumour-wide, lineage-restricted programs (generalising GBCD's PDAC classical/basal finding to all lineages, e.g., squamous vs glandular axes across HNSCC, lung and oesophagus). The ITH index places each MP on a continuum from ITH program to subtype program.

RESULT 3: MP REPERTOIRES
- The posterior co-membership network shows which MPs co-occur in or exclude each other across tumours.
- Its link to lineage, driver genotype where available, and clinical outcome, always tested against random-signature nulls.
- A per-MP normal-shared probability separates cancer-enriched programs from inherited lineage programs.

RESULT 4: PROSPECTIVE TEST (the wet-lab component)
- 8–12 fresh specimens from 2 cancer types, each split four ways: warm-protease scRNA-seq, cold-protease scRNA-seq, fixed-then-dissociated Flex, and snRNA-seq.
- Pre-registered CREDO predictions of which MPs are protocol-sensitive are scored after sequencing.
- This turns typing from a post-hoc label into a falsifiable prediction.

### Software architecture

LANGUAGE AND STACK
Python >= 3.11, PyTorch plus Pyro, built as an scvi-tools-compatible model class (the scvi-tools Pyro module base is assumed [U, check the current API]). AnnData/MuData native, with backed AnnData or TileDB-SOMA for out-of-core shards. A thin R package (credoR, through reticulate) imports usages as a Seurat assay/DimReduc and exports MP gene sets.

MODULES
- credo.pp
  - setup_anndata(adata, layer='counts', tumour_key, sample_key, cohort_key, lineage_key, protocol_key);
  - select_genes (deviance-based, plus TFs and optional known-MP genes);
  - add_ambient (CellBender/SoupX import);
  - add_doublets (scDblFinder or scrublet scores plus TME profiles per sample);
  - add_cna (Numbat / CopyKAT / SCEVAN / infercnvpy import into segment-dosage and clone-posterior tensors);
  - thin(eps=0.5, conditional_dispersion=True);
  - add_decoys(frac=0.15, method='scrp', strata=['clone', 'depth', 'doublet']).
- credo.model: the CREDO class, with tiers configured by n_global=200, n_private=4, n_cohort=4, quality_covariates=[...], zeta='auto', tumour_balance=2000.
- credo.infer: the SVI schedule, Hard-Concrete annealing, tier moves, polishing and local Bayes factors.
- credo.calibrate: select_zeta(), decoy_fdr(), null_cohorts(n=10), gene_tests(fdr=0.05), stability(seeds=5, folds=2, leave_cohort_out=True).
- credo.type: classify(rules='preregistered_v1'), which returns a type label plus its evidence vector.
- credo.ref: export_reference(name, version), which writes a Zenodo bundle with a checksum; project(query, ref, panel=None, novelty_slots=4).
- credo.sim: the CREDO-Sim simulator and benchmark harness (Snakemake, with containerised baselines including R tools).
- credo.pl: plotting.

TYPICAL USE
import credo as cr
cr.pp.setup_anndata(adata, ...)
cr.pp.add_cna(adata, source='numbat')
cr.pp.add_decoys(adata)
m = cr.CREDO(adata)
m.train(accelerator='gpu')
mps = m.metaprograms(alpha=0.1)
P = m.membership(fdr=0.1)
gs = m.gene_sets(fdr=0.05)
cr.type.classify(m)
m.export_reference('credo-pancancer-v1')
cr.ref.project(query, 'credo-pancancer-v1', panel=xenium_genes)

ENGINEERING
- Continuous-integration tests against simulator ground truth.
- Deterministic seeds.
- A figure-reproduction repository.
- The reference is versioned with its parent data hashes.

### Feasibility

TEAM
- 1 statistical-ML postdoc (model, inference, theory with a statistics co-advisor);
- 1 computational biologist (data harmonisation, CNA/ambient/doublet preprocessing, CREDO-Sim, baselines);
- 0.5 research software engineer (scaling, packaging, R bridge);
- a PI;
- a wet-lab or clinical collaborator for the prospective cross-protocol experiment (8–12 specimens; roughly USD 60–100k for sequencing and Flex, an estimate).

COMPUTE
- 4–8 A100/H100 GPUs for about 12 months, about 5,000–10,000 GPU-hours in total including calibration, ablations and GPU baselines.
- CPU and high-memory nodes for the R baselines (GBCD, Gavish code) and for Numbat.
- 20–50 TB of storage including allele-count intermediates.

DATA (all public)
3CA v1/v2, CanSig datasets, TCCA, CELLxGENE tumour sets, Kinker CCLE, Tahoe-100M, Replogle Perturb-seq, Janesick and Oliveira multi-platform data, Greenwald/Barkley spatial data, KP-Tracer.

CAVEATS
- Not all 3CA studies provide raw UMIs; some provide TPM or Smart-seq2 reads. Reads get their own protocol dispersion; TPM-only data are used for projection only.
- Allele-aware CNA needs BAMs, which are often controlled access; an expression-only caller on the thinning fold is the fallback.

TIMELINE (about 21 months to submission)
- M0–3: pilots P1–P7 and data and count audit.
- M3–9: core model v0.1, CREDO-Sim and the synthetic benchmark.
- M9–13: atlas fit, calibration and typing.
- M13–18: cross-protocol, Tahoe, spatial, genetic and TCCA validations, plus the prospective experiment.
- M18–21: preprint, package v1.0 and the Zenodo reference.

### Risks and mitigations

1. SHARED VS PRIVATE NOT IDENTIFIABLE IN PRACTICE
Mitigations:
- the cone-exclusion prior and tier moves;
- Proposition 1 and Proposition 2, with simulation sweeps that vary overlap and private-cone geometry;
- an identifiability margin and anchor genes reported per MP;
- low-margin MPs flagged as "not identifiable" rather than reported as confident.

2. PATIENT BASELINE CONFOUNDED WITH LIBRARY BATCH
Mitigations: require >= 3 tumours per cohort; identify protocol effects from multi-protocol tumours; tumour-wide programs are identified only through cross-tumour co-variation (eps is i.i.d.); sensitivity analysis with and without the baseline layer.

3. SVI AND AMORTISATION GAP, OVERCONFIDENT POSTERIORS
Mitigations: polishing plus local Bayes factors; all claims calibrated by decoys and count splitting, never by VI variances; HMC on a 20-tumour subset to check PIPs and loadings.

4. OVER-TEMPERING KILLS RARE MPs
Mitigations: zeta chosen by decoy FDR plus a one-SE deviance rule, with a curve of MP count against zeta always reported; rare-MP recall measured separately in CREDO-Sim; a rare-program mode (smaller n0, local zeta).

5. DECOY / NULL FIDELITY
SCRP keeps measured nuisance and marginal misspecification but not unmeasured artefact co-expression. This is stated openly and handled by protocol typing and cross-protocol data. Decoy recruitment is reported per MP as a diagnostic of an incomplete nuisance model.

6. VALIDITY OF THINNING (HIVERT ET AL.)
Mitigations: dispersion estimated conditional on the low-rank fit; decoy genes check calibration; an empirical-null correction; simulation showing calibration under mixtures. The extension of thinning theory to factor models is a possible theoretical side contribution.

7. CNA CIRCULARITY
Mitigations: allele-aware Numbat where BAMs exist; CNA called on fold A only; validation against DNA-defined clones. Pilot P3 may demote CNA from removal to typing.

8. OVER-CORRECTING REAL BIOLOGY (e.g., hypoxia correlated with %mito)
Mitigations: quality usage is confined to the covariate span; analysis run with and without nuisance terms; nuisance-explained variance reported per MP.

9. SCALE AND I/O
Mitigations: the two-stage fallback; shards per tumour; data-parallel training.

10. SCOOP
Risks: a count-native GBCD from the Stephens lab (log1p NMF already exists); the Parmigiani group porting multi-study NMF to scRNA-seq; ITHmapper publishing.
Mitigations: speed; differentiators no competitor has (decoy calibration, in-model nuisance and CNA, typing, the audit and the prospective test); preprint at M12. A manual bioRxiv / Google Scholar check is still mandatory (the landscape's Section 6.3).

11. "ANOTHER FACTORISATION" CRITIQUE FROM REVIEWERS
Mitigation: the headline is the definition plus error control plus typing plus the audit, backed by ablations showing each component fixes a named failure.

12. UNVERIFIED CITATIONS IN THIS DESIGN
Hierarchical beta process (Thibaux and Jordan), Hard-Concrete, randomised quantile residuals, deviance gene selection, target-decoy (Elias and Gygi) and the scvi-tools Pyro API are cited from memory and must be checked before writing. James/Lee/Pandey (arXiv 2304.05244), Weine et al. (arXiv 2601.05845) and Moran and Krishnan (arXiv 2601.18128) were confirmed only by search snippet.

### Go/no-go pilots

P1. FALSE-MP AND INSTABILITY AUDIT OF THE STANDARD (about 2 weeks, CPU)
- Run the Gavish/3CA code on 10 SCRP null cohorts built from 100 real tumours, and on planted-program cohorts.
- GO if it calls >= 3 MPs per null cohort OR leaks >= 20% of planted private programs.
- Otherwise pivot the headline to membership calibration, typing and edge MPs.

P2. PROTOTYPE VS GBCD / GAVISH / GeneNMF (about 6 weeks)
- NB plus membership plus private tier plus zeta; no CNA and no decoys yet.
- Run on CanSig (about 185 patients) with planted subset-shared programs at 20–40% prevalence.
- GO if membership AUROC >= 0.9 AND matched cosine >= GBCD + 0.05 for subset-shared programs AND runtime <= 4 GPU-hours.

P3. CNA / POSITIONAL AUDIT (about 3 weeks)
- Arm enrichment and clone-explained usage variance for 3CA robust programs and MPs.
- If < 5% of MPs are affected, keep the dosage offset but present CNA as typing rather than removal.

P4. DECOY SANITY (about 3 weeks)
- On a real cohort with 15% SCRP decoys at zeta = 1, decoys should be recruited by some slots, which demonstrates misspecification. Tuned zeta should bring the decoy membership FDR to <= 0.1 while keeping >= 90% of the zeta = 1 MPs that match known Gavish MPs.
- If decoys are never recruited, drop the power posterior from the headline (not fatal).

P5. SCALING (about 4 weeks)
- Synthetic 3,000 tumours x 3,000 cells x 6,000 genes on one 80 GB GPU in <= 24 h.
- If not, use the two-stage fallback.

P6. DATA AND COUNT AUDIT (about 2 weeks)
- Measure the fraction of 3CA v2 tumours with raw UMIs and BAMs.
- If < 60% have UMIs, rebuild the discovery cohort from CELLxGENE plus TCCA raw data.

P7. ARTEFACT PREDICTABILITY (about 1 week)
- Do simple nuisance features (ambient cosine, doublet score, %mito/ribo, sc vs sn) separate the 11 hand-removed 3CA v2 MPs from retained ones (AUROC >= 0.8)?
- If yes, the typing target is achievable. If no, refine the nuisance model before committing to the claim that typing replaces curation.

### Why Nature Methods

It turns a concept used across the field into a statistical object with a stated estimand and error rates. 3CA MPs sit in MSigDB and drive module scoring in hundreds of papers, yet they rest on overlap thresholds and hand curation. CREDO provides:
- a definition: a program shared by a credible subset of tumours that beats nuisance-preserving decoys;
- an estimator: a count-native, multi-level, subset-sharing NMF at atlas scale;
- calibrated inference: target-decoy existence and membership FDR, count-split gene FDR;
- automatic typing that replaces manual curation.

It meets every point of the Nature Methods bar distilled in the landscape:
- a conceptual advance generic tools cannot match;
- simulations that plant the failure modes;
- more than 10 fairly tuned competitors, including GBCD, TCAT-style, GeneNMF and the Gavish code;
- orthogonal validation (cross-protocol, perturbation, spatial, DNA-clone, lineage, replication cohort, clinical against random nulls);
- an atlas-scale biological audit with a prospective, pre-registered experiment;
- a pip/scverse package, an R bridge, a versioned reference, a projector and an open benchmark.

The framework is general. "Programs shared by subsets of donors, with nuisance modelled and error control" applies directly to multi-donor disease atlases, organoid and perturbation screens, and immune programs across patients. That widens its readership beyond cancer genomics, which is what separates a Nature Methods paper from a Genome Biology or Cancer Research resource.

### Differentiation from closest competitors

| Competitor | Differentiator |
|---|---|
| GBCD (Liu…Stephens, Nat Genet 2025) | GBCD is a Gaussian covariance decomposition of log data with near-binary cell memberships, a user-set Kmax and one cancer type per analysis. It has no tumour-membership variables, no recurrence or existence test, no nuisance or CNA model, no projector, and no demonstrated atlas scale. CREDO differs on each point: an NB likelihood on raw UMIs; explicit tumour-level inclusion z_tk with PIPs and FDR-calibrated credible tumour sets; a hierarchy across lineages (one pan-cancer fit); separate cohort-private and tumour-private dictionaries instead of 'patient-specific factors' that silently absorb nuisance; CNA, ambient, doublet and protocol terms; target-decoy existence q-values; and SVI to 10^7 cells. A count-likelihood ablation and a Gaussian-log ablation of CREDO isolate how much the likelihood contributes compared with the structure. Watch: Weine, Carbonetto, Irizarry and Stephens' Poisson 'log1p NMF' (arXiv 2601.05845 [V title/authors only]) signals that the GBCD group is moving to count models. |
| ITHmapper (Abelson lab, bioRxiv 2026) | ITHmapper uses scVI with sample as batch, then Hotspot modules and discrete states. Treating sample as a batch can erase tumour-wide programs, there is no uncertainty, and CNA is used only for malignant calling. CREDO models continuous, overlapping programs, keeps tumour-wide programs through a baseline layer with an ITH index, attaches calibrated membership and gene FDR, and models CNA, ambient and doublets in-model. Benchmarked head-to-head on the same tumours. |
| TCAT/starCAT (Nature Methods 2025) | TCAT runs cNMF on Harmony-corrected data per dataset, then clusters programs. Singletons are rejected by heuristic plus judgement, the scope is T cells only, and NNLS projection has no uncertainty. In CREDO, singleton and cohort programs are model tiers, recurrence and existence are inferred, and the projector returns PIPs, usage uncertainty, novelty slots and protocol-aware panel transfer. The reference is a deliverable, not the claim. |
| Grabski–Trippa–Parmigiani Bayesian multi-study NMF (Genome Biol 2025) and Tetris/BMSFA | Those models are Poisson or Gaussian with a flat study-by-factor inclusion matrix, use MCMC, and have about 96 features and hundreds of samples. CREDO nests cells within clones within tumours within cohorts within lineages under a hierarchical beta-Bernoulli. It adds an NB likelihood with protocol dispersions, private and cohort tiers with cone exclusion, dosage, ambient and doublet offsets, a tumour-balanced power posterior, amortised SVI with tier moves for 10^7 x 6,000 matrices, local tempered Bayes factors, and target-decoy calibration. Tetris/BMSFA are run on subsamples to document their scaling limits. |
| scDEF (Beerenwinkel lab, bioRxiv) | scDEF's hierarchy runs over factor abstraction (coarse-to-fine Gamma layers) and is batch-aware, but shared versus patient-specific is only implicit and the demonstration covers 9 patients. CREDO's hierarchy follows the sampling design (lineage, cohort, tumour, clone, cell), with explicit inclusion variables, nuisance and CNA terms, error control and atlas scale. |
| LIGER iNMF / UINMF / consensus iNMF / scINSIGHT | In iNMF, V_i is a perturbation of the same K factors, controlled by a hand-set lambda. It uses a Frobenius loss, has no subset sharing, needs a hand-picked K, and the W/V split becomes weak with about 1,000 datasets. scINSIGHT needs condition labels and allows only 'all or one' sharing. CREDO has distinct private dictionaries with cone exclusion, sharing by arbitrary subsets of tumours with posteriors, a count likelihood and learned K. |
| MOFA+ / MuVI / MOFA-FLEX | These models have real-valued (not non-negative) loadings shared by all groups, with group-wise ARD only on factor activity. They are designed for few groups, not 3,000. Patient-private programs need a global slot, and there are no membership FDRs, no nuisance or CNA terms, and no count-native non-negative programs. |
| BayesPowerNMF / bayesNMF / compressive NMF | These are single-study models of mutational signatures fitted by MCMC. CREDO borrows the power-posterior idea but makes it tumour-balanced (w_t), calibrates zeta with in-fit decoys rather than a grid by eye, and embeds it in a multi-level subset-sharing model at single-cell scale. |
| Gavish/3CA code, GeneNMF, scalop, cNMF, consensus scHPF, CSI-GEP | These are per-tumour or pooled factorisations followed by thresholds or clustering. They have no likelihood-based recurrence, no error control and no nuisance model, and GeneNMF forces every program into an MP. CREDO is a single generative model with inferred tiers and calibrated calls. Speed or determinism is not the claim. |
| Moran & Krishnan nonlinear multi-study factor analysis / multi-study sparse VAE (arXiv 2601.18128; ICASSP 2026) | It is nonlinear and not non-negative, separates shared from specific factors across a few disease groups, and gives no uncertainty on sharing. CREDO is linear, non-negative and count-native, with interpretable loadings, subset sharing over thousands of tumours, and calibrated membership. |
| DRVI / contrastiveVI / DeltaNMF | DRVI gives disentangled nonlinear factors with no shared-versus-private structure. contrastiveVI and DeltaNMF need a matched background. In CREDO, the normal cell of origin enters as extra lineage nodes in the same beta-Bernoulli hierarchy, giving a per-MP 'normal-shared' probability as a soft annotation instead of subtraction. Foundation-model priors are deliberately not part of the claim. |
| ClusterDE / recall / countsplit | These calibrate clustering and differential expression. CREDO adapts their ideas to factor programs and recurrence across tumours. SCRP decoys keep measured nuisance and the empirical marginal misspecification but destroy co-expression, and they are fitted jointly as decoy tumours, a new use. NB count splitting uses dispersion conditional on the low-rank fit to address Hivert et al., with decoy genes as an online calibration check. |

### Referee 1: novelty and prior art

**Predicted editorial outcome:** My prediction is that the paper is probably sent to review, roughly 65-70% versus about 30% desk reject, and is then most likely rejected or returned for a major revision that amounts to a resubmission. The editor will send it out for three reasons: 3CA/Gavish MPs are used across the field and sit in MSigDB, there is no calibrated alternative, and the package-plus-benchmark-plus-reference bundle fits the Nature Methods template (TCAT/starCAT, CellBender). Referees will then find the following. (i) The headline 'error-controlled existence' rests on decoys that are not exchangeable with real non-member tumours, so the q-values are close to vacuous at 3,000 tumours. (ii) The calibration evidence comes from a simulator built with the same operator (SCRP) and the same functional forms as the model. (iii) The 'no thresholds' claim is contradicted by the tier-move rules, whose cosine 0.7/0.9, >=3-tumour, >=2-cohort and >90%-one-cohort cut-offs restate the Gavish inter-tumour and single-study filters. (iv) In the 3CA corpus, cohort, lineage and protocol are largely aliased, so protocol typing and lineage-restricted MPs may not be identifiable. (v) Every module has precedent: Grabski multi-study NMF with a Tetris-like prior, BayesPowerNMF, count splitting, ClusterDE, CellBender/scDblFinder, and ISOMERIC for CN-dependent variation. The contribution is integration, so the paper lives or dies on a large, orthogonally validated biological reversal. If the audit mostly re-confirms the canonical MPs (cell cycle, hypoxia, EMT, IFN, stress and others), which is the most likely outcome, I expect a decision of the form 'a more complex method reproduces known programs with q-values' and redirection to Genome Biology or Nature Communications. A path to acceptance exists only if the authors fix the decoy construction, validate on an independent generator and real ground truth, cut scope, and win head-to-head on cross-protocol and DNA-clone criteria.

**Fatal flaws**

- THE EXISTENCE AND MEMBERSHIP 'ERROR CONTROL' DOES NOT CONTROL THE ERROR USERS CARE ABOUT, AND ITS FDR ESTIMATOR IS BIASED DOWNWARD. SCRP decoys keep only the covariance induced by measured nuisance and destroy all other co-expression. Real non-member tumours, by contrast, carry every other program. (1) The null H0_k ('k explains no more held-out deviance in real than in decoy tumours') is rejected by any slot that captures any non-nuisance co-expression, whether biological or unmeasured artefact (dissociation IEG/heat-shock, residual ambient misfit, histone-cluster coupling). With ~2,500 real tumours and 10^7 cells every surviving slot will have q close to 0, so q_k cannot separate biology from artefact. This repeats the proposal's own critique of inter_min=10 one level up: a null 'that no real dataset follows'. (2) The target-decoy assumption (Elias-Gygi: an incorrect target match is as likely as a decoy match) fails. A slot can partially absorb correlated programs in real non-members (hypoxia/glycolysis/stress, or the nested pairs at cosine 0.4-0.6 that CREDO-Sim itself plants), but has nothing to latch onto in a structureless decoy. Decoy PIPs will therefore sit near 0, and FDRhat_k(thr) will underestimate the false-membership rate. (3) Decoys are paired copies of real tumours, so permuting real/decoy labels at tumour level is not exchangeable unless done within pairs. (4) zeta is chosen by decoy FDR, so this flaw spreads to the power posterior and to every PIP.
- SELF-VALIDATING SIMULATION (INVERSE CRIME). CREDO-Sim builds its backbones with SCRP, the same operator that makes the decoys, and the 20 null cohorts are SCRP decoys, which is exactly the null CREDO is calibrated against. Planted programs are additive and non-negative on the model's own rate scale. Ambient (a fixed per-sample profile times a fraction), doublets (a sample-matched TME mixture), CNA (a dosage-power offset), depth and sc/sn efficiency are injected with the model's own functional forms. The pre-registered targets (<=1 MP on null cohorts, calibrated membership and gene FDR, false-MP rate <=5%) are then close to guaranteed by construction, and losses by Gavish, GBCD or ITHmapper on this simulator say nothing about real data. Without an independent generator (scDesign3/splatPop-type, or real-data spike-ins with multiplicative or nonlinear effects, non-additive ambient and misspecified dispersion) and real ground truth, the primary evidence ('A') is not evidence.
- THE CENTRAL CONCEPTUAL CLAIM ('only interpretable constants; replaces 35/10/10 and greedy clustering') IS FALSE AS DESIGNED. The thresholds are moved, not removed. The tier moves promote private programs from '>=3 tumours in >=2 cohorts with cosine > 0.7', which is the Gavish inter-tumour overlap filter in cosine form. 'Global slot with >90% of members in one cohort -> cohort tier' is the 3CA single-study filter. The rules also merge at cosine > 0.9, demote at R_k < 1.5 and split on bimodality. Further unreported settings: K=200, R=4, J=4, n0=2,000, Beta concentration 50, Beta(5,5), Gamma(0.3,0.3), lambda_ce, c and alpha0, the Hard-Concrete temperature schedule, the zeta grid with a 1-SE rule, and K=200 KL-NMF initialisation on metacells. Under non-convex mean-field SVI these settings, not the likelihood, will decide the catalogue. 'Pre-registered typing rules' are curation written down as a decision tree. As stated, the 'MP as a statistical object' headline does not hold.
- THE TARGET CORPUS PROBABLY CANNOT IDENTIFY THE TYPINGS THAT CARRY THE FLAGSHIP. In 3CA, study (cohort), cancer type (lineage) and protocol (sc/sn/Smart-seq2) are largely nested: most studies cover one cancer type with one protocol, single-nucleus is concentrated in particular tissues, and tumours profiled under >=2 protocols are rare. eta_{k,p} (protocol sensitivity), pi_{k,l} (lineage profile) and the cohort tier U_c are then aliased. The >90%-one-cohort move and the >=2-cohort recurrence rule will push genuine lineage-restricted programs from cancer types covered by a single study into the 'cohort artefact' tier by construction. That cuts directly against Result 2 (lineage-restricted tumour-wide programs) and the protocol-sensitivity typing. Pilot P6 checks UMI and BAM availability but not the rank of the cohort x lineage x protocol design, and that check should be a go/no-go.

**Major concerns**

- NOVELTY IS INTEGRATIVE, AND EACH MODULE HAS A NAMED PRECEDENT. Subset-sharing inclusion: Grabski-Trippa-Parmigiani multi-study NMF, which uses a Tetris-like prior for shared ownership; Tetris; and James/Lee/Pandey's hierarchical IBP for within- and across-group sharing. Count-native hierarchical NMF: scHPF and scDEF. Shared vs patient-specific: GBCD and multiGroupVI. Tempering: BayesPowerNMF. Held-out folds: data thinning and NB count splitting. Synthetic nulls: ClusterDE. Ambient and doublets: CellBender, SoupX and scDblFinder. CN-dependent vs independent variation: ISOMERIC. Baseline vs within-tumour programs: GBCD and Nomura et al. (Nat Genet 2025). Catalogue plus projector, including explicitly catalogued technical programs: TCAT/starCAT. The genuinely new pieces (decoy tumours inside the fit, the w_t-balanced tempering, cone exclusion, sample-matched TME-doublet terms) are modest twists. Reviewers will ask which single idea produces a large real-data effect. The ablation plan answers only on CREDO-Sim, which is the self-validating simulator.
- CNA CIRCULARITY IS NOT SOLVED BY COUNT SPLITTING. Thinning makes the noise independent, not the signal. Most 3CA tumours lack accessible BAMs, so expression-only callers (inferCNV/CopyKAT/SCEVAN) will supply D_ng for most of the corpus. These callers smooth along the genome and flag any positionally clustered program as a CNA. Major biological programs are positionally clustered: the epidermal differentiation complex (1q21), type I and II keratins (17q21, 12q13), MHC-II (6p21), replication-dependent histones (6p22), protocadherins and HOX clusters. The dosage offset will therefore attenuate squamous differentiation, antigen presentation and S-phase histone components, in exactly the lineages Result 2 relies on (squamous axes in HNSCC, lung and oesophagus). Result 1c ('1q in breast') is confounded with the method's own mechanism. In addition, a single pan-cancer gamma_g contradicts the 2026 joint scDNA+RNA preprint reporting cancer-type-specific dosage compensation, and misfit from a wrong gamma flows back into the programs.
- THE BASELINE-VERSUS-TUMOUR-WIDE-PROGRAM SPLIT IS DECIDED BY THE PRIOR. For a homogeneous program (alpha_tk large), the likelihood does not change when m_tk W_k is moved into beta_{t,g}, which has about 18M free parameters. Whether a shared shift surfaces as an MP or disappears into baselines depends entirely on Gamma(tau,tau), the sparsity prior and slot availability. beta_t is also entangled with D_ng (clone dosage) and nu_{c,g}, and where tumour equals library (common in 3CA) it is confounded with batch. 'Tumour-wide subtype programs surface as MPs' is therefore an assumption, not a result. It needs simulations where the truth is baseline-only versus a shared program, sensitivity to tau, and evidence on matched multi-sample tumours.
- THE MEMBERSHIP ESTIMAND IS BIASED AGAINST CLASSICAL ITH PROGRAMS, AND THE PIPs ARE NOT PROBABILITIES. Membership requires a tumour-mean usage >= 0.5% of biological UMIs. A program used by 3% of cells at 15% usage (persister-like, cycling cells in slow tumours, an IFN-high minority) falls below that and is declared absent. Membership should be defined on the cell-level usage distribution (prevalence x effect). The tempered local Bayes factor (zeta*w_t*dELBO, with zeta possibly 0.02-0.05 and w_t about 0.2 for a 10k-cell tumour) makes the PIPs arbitrary functions of the tuning. 'Credible tumour sets' then rest entirely on the decoy calibration, which is flawed (fatal flaw 1).
- SHARED/PRIVATE IDENTIFIABILITY IS ASSERTED, NOT SHOWN. Proposition 1 (SSC is monotone under row union) is a near-trivial cone-containment observation in the noiseless case. Proposition 2 is 'to be proven' and relies on exclusive anchor genes and a soft cos^2 penalty. Cone exclusion forces tumour-specific variants of a shared program (an amplicon-flavoured MYC or EMT program) either into the shared W or into an orthogonalised private slot, which distorts both. The 'identifiability margin' is a diagnostic, not a guarantee, under sparse overdispersed counts with millions of nuisance parameters.
- GENE-LEVEL FDR WILL BE UNINFORMATIVE AT ATLAS SCALE AND RESTS ON SHAKY THINNING ASSUMPTIONS. A one-sided point null (beta_gk = 0) at 10^7 cells will be rejected for thousands of genes, so the 'variable-length FDR gene sets' will grow far past usefulness for scoring. An interval or effect-size null is required. NB thinning is exact only with known dispersion. Estimating phi from a pilot low-rank fit on the full data, and then thinning, re-introduces fold dependence. Score tests that plug in fold-A usages and PIP weights as fixed ignore their uncertainty. Decoy genes permuted within tumour lack the correlation structure of real non-member genes (the same flaw as recall knockoffs).
- THE INFERENCE MAY BE AS SEED- AND INITIALISATION-DEPENDENT AS THE METHOD IT CRITIQUES. The fit combines mean-field SVI, an amortised encoder, Hard-Concrete relaxation of about 3,000x200 inclusion variables plus cohort inclusions, accept/reject tier moves on held-out ELBO, and initialisation from a pooled K=200 KL-NMF. The final dictionary may largely reflect that initialisation. Polishing with globals fixed ignores global uncertainty. HMC on 20 tumours cannot validate 3,000-tumour behaviour. A stability index from 5 seeds with Hungarian matching of 200 slots is a weak test of the instability the proposal attributes to Gavish.
- THE ORTHOGONAL VALIDATION IS MOSTLY NON-COMPARATIVE OR NON-DISCRIMINATIVE. 'CDK4/6 lowers cell-cycle MPs' in Tahoe-100M will also be passed by Gavish lists. Cross-protocol rests on two vendor datasets with few patients and probe or panel gene sets (Janesick, Oliveira). The prospective four-way split mainly predicts an established result: warm-collagenase IEG/heat-shock signatures are reduced by cold protease and fixation (O'Flanagan et al. 2019). Each arm needs a pre-specified primary endpoint on which CREDO-defined MPs are compared head-to-head with Gavish/3CA, GBCD and ITHmapper MPs, plus non-obvious predictions (which specific memberships or genes flip).
- THE ARTEFACT-TYPING BENCHMARK IS CIRCULAR AND UNDERPOWERED. Eleven hand-removed 3CA v2 MPs, plus a handful from 2023, give perhaps 15-20 positives, labelled by one lab using the same features the typing rules use (ribo/mito share, TME markers). An AUROC on 15-20 positives has very wide confidence intervals. CREDO's ambient and doublet terms are designed to stop these programs forming at all, so there may be nothing to type. The authors must specify how 3CA MPs map into CREDO (projection? matching?) before the metric is defined.
- THE EXPECTED BIOLOGICAL PAYOFF IS UNCERTAIN. Most canonical MPs (G2/M, G1/S, hypoxia, EMT-like, IFN, MHC-II, UPR, stress, respiration) are strong, reproducible signals. The 'not real biology' problem is concentrated in contamination and low-quality MPs that were already curated out. The prior probability that a far more complex model overturns many canonical MPs is low. Pilot P1 tests Gavish on SCRP nulls, where per-tumour NMF on structureless residuals may well produce few recurrent programs. A GO/NO-GO judged on SCRP is again circular. The decisive pilot should use real nuisance contrasts (matched sc/sn/Flex tumours and cold versus warm dissociation).
- THE SCOPE AND FEASIBILITY ARE NOT CREDIBLE FOR ONE PAPER, 2.5 FTE AND 21 MONTHS. The plan covers a new model, two propositions plus 'thinning for factor models' theory, a new simulator, about 17 baselines (several in R, including Tetris), a 10^7-cell atlas fit with 15% decoys x 2 folds x leave-one-cohort-out x a zeta grid x 5 seeds x 10 null cohorts, seven validation arms, a wet-lab prospective study, a versioned reference with projector and panel transfer, a pan-cancer dosage resource, normal-lineage nodes, and an R bridge. The GPU estimates (6-24 H100-hours per fit; 1,500-4,000 GPU-hours for calibration) are unverified, and per-query projection needs CellBender, scDblFinder and Numbat preprocessing, which limits adoption.
- THE ADDITIVE MIXTURE IS MISSPECIFIED FOR KEY NUISANCES. Ambient is a fixed SoupX/CellBender profile with rho concentration 50 around a noisy sample-level estimate. Dissociation stress is not modelled at all (it is left to typing). The QC tier is confined to the covariate span, so hypoxia correlated with %mito and complexity is either absorbed or leaks. Residual misfit is exactly what produces recurrent junk components, and the power posterior is the only safeguard, tuned by the flawed decoys.

**Prior-art threats**

- GBCD (Liu, Carbonetto, Willwerscheid, Oakes, Macleod, Stephens; Nat Genet 57:263-273, 2025; s41588-024-01997-z) [VERIFIED: landscape V+ and my search]: label-free separation of shared, dataset-specific and patient-specific tumour programs, including PDAC subtype programs missed by per-tumour NMF. It pre-empts the 'shared vs patient-specific' and 'tumour-wide subtype programs' (Result 2) framing.
- Grabski, Trippa, Parmigiani, Bayesian multi-study NMF (Genome Biol 2025) [VERIFIED: Springer link]: a Poisson NMF across studies with a Tetris-like prior 'explicitly modeling the shared ownership of signatures'. This is the z_tk subset-inclusion core of CREDO. Probit follow-up (arXiv 2502.01468) [VERIFIED title via search].
- Tetris (Grabski et al., AoAS 2023) [VERIFIED in landscape] and James, Lee, Pandey, 'Bayesian Analysis of Generalized Hierarchical Indian Buffet Processes for Within and Across Group Sharing of Latent Features' (arXiv 2304.05244) [VERIFIED title via search]: hierarchical beta-Bernoulli and IBP subset sharing across groups is already formalised.
- ISOMERIC (Selega, Maan, Yu, Tan, Campbell; bioRxiv 10.1101/2025.04.21.649842, 2025) [VERIFIED via search snippets; full text blocked]: an adversarial VAE that splits single-tumour-cell expression into copy-number-dependent and -independent components, trained on 189,027 cells from 10 cohorts in 4 cancer types. The landscape MISSED this. It directly undercuts 'no packaged method discovers programs jointly with CNA' and the G2 novelty.
- 'Principles of subclonal gene dosage across human cancer' (bioRxiv 10.64898/2026.03.15.711906, 2026) [VERIFIED title and snippet; authors not retrieved]: joint single-cell WGS plus mRNA in 57 patients; dosage is mostly additive, and compensation is cancer-type-specific. This contradicts a single pan-cancer gamma_g and makes the 'dosage-sensitivity resource' weaker than DNA-anchored estimates. Also Fehrmann et al., 'Gene expression analysis identifies global gene dosage sensitivity in cancer' (Nat Genet 2015) [VERIFIED title via search].
- CLCC, 'Contrastive Alignment of Expression and Copy Number Highlights Dosage-Insensitive Genes in Cancer' (bioRxiv 10.64898/2026.03.01.708901) [VERIFIED title]: dosage-insensitive (escape) genes from aligned scRNA and CNV.
- Nomura, Spitzer, Johnson, Garofano et al., 'The multilayered transcriptional architecture of glioblastoma ecosystems' (Nat Genet 2025; s41588-025-02167-5) [VERIFIED via search; the landscape had it as P]: in 121 GBMs, separates within-tumour malignant states (MPs) from three between-tumour 'baseline programs'. This undercuts the claim that baseline/tumour-wide programs are unformalised (G5, Result 2), although it is not a general method.
- ITHmapper / Pan-Cancer Single-Cell Compendium of ITH (Nofech-Mozes & Abelson; bioRxiv 10.64898/2026.01.06.693992) [VERIFIED]: 1,021 tumours, 14 cancer types, 3.7M cells, scVI with sample as batch plus Hotspot, 351 modules, a mapping tool, and Zenodo data. It takes the 'pan-cancer ITH package plus mapper' positioning.
- TCAT/starCAT (Kotliar...Raychaudhuri; Nat Methods 2025, s41592-025-02793-1) [VERIFIED]: a reproducible cross-dataset catalogue plus projector that already catalogues 5 technical-artefact programs alongside 41 biological ones. This is precedent for typing artefacts inside a catalogue rather than deleting them.
- Weine, Carbonetto, Irizarry, Stephens, log1p (shifted-log-link) Poisson NMF (arXiv 2601.05845, Jan 2026) [VERIFIED title and authors]: the GBCD group already has count-native NMF. A count GBCD would remove CREDO's 'count likelihood vs Gaussian GBCD' differentiator.
- BayesPowerNMF (Xue, Miller, Carter, Huggins; PLoS Comput Biol 2026) [VERIFIED]: a power posterior with a sparsity prior for misspecification-robust NMF with automatic K. CREDO's zeta mechanism is borrowed directly.
- multiGroupVI (Weinberger, Lopez, Huetter, Regev; PMLR v200 / MLCB 2022) [VERIFIED proceedings; a Nature Methods 2024 version was claimed by a search summary and is UNVERIFIED]: shared versus group-specific latent factors in grouped scRNA-seq.
- ClusterDE / synthetic null (Song, Li et al.) [VERIFIED]: synthetic-null calibration with scDesign3. Note that it preserves gene-gene correlation, whereas SCRP destroys it. The 'decoy' idea is an adaptation, not a new principle.
- NB count splitting (Neufeld, Popp, Gao, Battle, Witten; arXiv 2307.12985) and data thinning / generalised data thinning (JMLR 2024; JASA 2025) [VERIFIED]; Hivert et al., limits of thinning (arXiv 2405.13591) [VERIFIED].
- scDEF (Ferreira, Kuipers, Beerenwinkel; bioRxiv 2022 v2) [VERIFIED; still appears to be a preprint]: multi-level Gamma-Poisson with automatic size, batch-aware, shared and patient-specific states in HGSOC.
- SCATrEx, 'Mapping single-cell transcriptomes to copy number evolutionary trees' (bioRxiv 2021.11.04.467244) [title VERIFIED via search; authors and venue UNVERIFIED]: joint model of a CNA tree plus non-CNA expression events. Precedent for clone-linked versus non-genetic typing.
- O'Flanagan et al., cold active protease dissociation (Genome Biol 2019) [VERIFIED via search]: a 512-gene collagenase-induced stress set minimised by cold dissociation. The prospective experiment's main 'prediction' is already established.
- CanSig benchmark (Cancer Res 86(4):873, 2026) [VERIFIED in landscape by one sweep only]: an existing cancer-state discovery benchmark (13 methods, 12 datasets). It reduces the novelty of CREDO-Sim and the harness as a contribution.
- KINOMO (bioRxiv 2022; AACR 2024) [VERIFIED in landscape]: KL-NMF with kernel correction for inter-patient genomic variability; must be acknowledged on the CNA angle.
- 3CA v2 (Nature Cancer 2025) [VERIFIED] and the unattributed 2026 MP catalogue (JerryChen06 repository) [VERIFIED in landscape; authors unknown]: protocol-named MPs ('single-nucleus', 'in vitro') show the field already types protocol dependence by hand, and the baseline to beat is moving.

**Required changes**

- REBUILD THE NULL. Replace structureless SCRP decoys with decoys that keep all programs except the one under test (leave-program-out conditional decoys, e.g., permuting only the residual after removing slot k's contribution), or with real negative-control tumours. Show by simulation that decoy PIP distributions match the PIPs of known non-member real tumours when programs are correlated or nested. Make decoy permutation tests paired. Until then, call the q-values 'rejecting measured-nuisance-plus-noise' and drop the term 'error-controlled existence'.
- MAKE THE EXISTENCE ESTIMAND DISCRIMINATIVE. Define existence as invariance (reproducibility across protocol and cohort strata or across thinning folds with an effect-size floor), or test against a null that includes unmeasured artefact co-expression (real cold-vs-warm and sc-vs-sn contrasts). Report power and specificity for planted artefacts that are NOT in the nuisance model.
- VALIDATE ON AN INDEPENDENT GENERATOR AND REAL GROUND TRUTH. CREDO-Sim must not share generative forms or the null operator with the model. Add scDesign3/splatPop-style and real-spike-in simulations with multiplicative and nonlinear effects, non-additive ambient and misspecified dispersion. Have a blinded third party design adversarial scenarios. Report every metric on both simulators, and present the SCRP-based simulator only as a best case.
- AUDIT ALL HIDDEN DECISIONS. Report sensitivity of the MP catalogue to every tier-move threshold, K/R/J, n0, priors, lambda_ce, initialisation and annealing. Replace accept/reject tier moves with a principled split-merge scheme or show that the catalogue does not change. Retract the 'only interpretable constants' claim.
- ADD A PILOT P0 FOR THE CORPUS DESIGN RANK. Tabulate cohorts x lineages x protocols and multi-protocol tumours. Restrict protocol typing and lineage profiles to the strata where they are identifiable. Report which lineage-restricted MPs depend on a single cohort, and do not demote them to the artefact tier by rule.
- FIX THE CNA HANDLING. Evaluate all CNA claims against DNA- or allele-anchored copy number (scDNA, Numbat with BAMs, clonealign). Use known positionally clustered biological programs (EDC 1q21, KRT 17q21/12q13, MHC-II 6p21, histone 6p22) as negative controls that must NOT be removed. Allow cancer-type-specific gamma. Default to CNA typing rather than an offset unless P3 shows a large effect. Benchmark against ISOMERIC and cite the 2026 subclonal-dosage preprint and CLCC.
- REDEFINE TUMOUR MEMBERSHIP on the cell-level usage distribution (prevalence above an effect size) and report recall for minority-subpopulation programs. Stop calling tempered-LBF quantities 'posterior inclusion probabilities' unless their calibration is shown on real held-out labels.
- USE INTERVAL (EFFECT-SIZE) NULLS FOR GENE MEMBERSHIP. Show that gene sets have a usable size and replicate out of sample (3CA to TCCA). Demonstrate thinning validity when dispersion is estimated from the same data, or use a pre-specified external dispersion estimate.
- MAKE EVERY ORTHOGONAL ARM HEAD-TO-HEAD WITH PRE-SPECIFIED PRIMARY ENDPOINTS: CREDO MPs versus Gavish/3CA, GBCD, ITHmapper and TCAT-style MPs on cross-protocol concordance, DNA-clone orthogonality, spatial niche specificity and TCCA replication. Rewrite the prospective experiment's predictions to be non-obvious (specific membership or gene flips, in vivo stress versus ex vivo stress components), not 'IEG falls with cold protease or fixation'.
- STRENGTHEN THE ARTEFACT-TYPING EVALUATION. Expand the labelled set beyond the 11 3CA v2 MPs (TCAT technical programs, protocol-contrast programs, sn-vs-sc programs). Use features the curators did not use, report confidence intervals, and define exactly how 3CA MPs are mapped into CREDO.
- CUT THE SCOPE. Make the core paper the count-native subset-sharing model, the in-model ambient, doublet and dosage terms, calibrated membership, and one decisive real-data audit with cross-protocol and DNA-clone validation. Move the reference and projector with panel transfer, the gamma_g resource, normal-lineage nodes, the R bridge and the full 17-baseline harness to follow-ups or supplements. Report measured P5 scaling numbers before making atlas-scale claims.
- THEORY: prove Proposition 2 under noise or drop it. Present Proposition 1 as an observation. Give empirical identifiability curves under realistic overdispersion and overlap.
- REDESIGN PILOT P1 AROUND REAL NUISANCE (matched protocol tumours, cold-versus-warm datasets, leave-one-study-out instability of the Gavish code) rather than SCRP nulls, and commit in advance to the venue change that follows if most canonical MPs survive.
- CITATIONS AND NAMING. Cite ISOMERIC, Nomura et al. 2025, multiGroupVI, the subclonal-dosage preprint, CLCC, SCATrEx, O'Flanagan 2019 and Fehrmann 2015. Verify all [U] citations (Thibaux-Jordan, Hard-Concrete, randomised quantile residuals, Townes deviance, Elias-Gygi, the scvi-tools Pyro API). Check the name 'CREDO' for collisions (unchecked).

**Ideas worth keeping**

- Explicit tumour-level inclusion variables (z_tk) with an MP x tumour membership matrix. 'Which tumours carry this program, and how sure are we' is the real biological claim and the true gap left by GBCD, 3CA and ITHmapper. It should be the core of the paper, with a corrected calibration.
- In-model, sample-matched TME-doublet (M_{s,e}) and empty-droplet ambient (A_s) terms, which absorb or type the macrophage, T-cell and endothelial 'MPs' that 3CA v2 removed by hand. This is concrete, falsifiable and easy to show head-to-head.
- Tumour-balanced likelihood weighting (w_t = min(1, n0/n_t)) combined with tempering, so that 20k-cell tumours cannot dictate the dictionary. It is cheap and worth keeping even in the two-stage fallback.
- The ITH index (within-tumour share of usage variance) on an uncentred, common scale. It places every MP on a continuum from ITH program to subtype program, which unifies the Gavish-style and Nomura/GBCD-style baseline views.
- Cross-protocol invariance as the operational definition of 'real biology', including the four-way split (warm, cold, Flex, sn) prospective design, provided the predictions are made non-obvious.
- Per-MP diagnostics reported rather than hidden: decoy recruitment rate as a misspecification signal, nuisance alignment (maximum cosine to ambient, doublet and QC profiles), positional-arm enrichment with an expression-matched permutation, stability across seeds x thinning folds x leave-one-cohort-out, and anchor genes and near-pure cells.
- Count-split held-out deviance for choosing model complexity, with decoy genes as a running check of calibration.
- Separate cohort-private and tumour-private tiers, so that study-specific artefacts and tumour-private (e.g., amplicon) programs are catalogued and labelled rather than silently promoted or discarded.
- A semi-synthetic simulator with confounders anchored in measured quantities (real Numbat segments, real empty-droplet profiles, real TME doublets, real protocol shifts) as a community benchmark, if it is decoupled from the model's own generative forms and null operator.

### Referee 2: statistical rigor and feasibility

**Predicted editorial outcome:** Likely sent to review, then major revision trending to reject. The topic and scope fit Nature Methods: 3CA MPs are used everywhere, and GBCD and TCAT show editors want this kind of work. A model-based definition with error rates plus an atlas audit is the right pitch. But the two headline statistical claims do not test what they name. (1) The 'error-controlled existence call' uses SCRP decoys, which keep only nuisance that CREDO already models. (2) The 'decoy-calibrated membership FDR' assumes non-member real tumours are exchangeable with decoys. A statistical referee will show that the existence null is false for any co-expression, so q-values will sit at the permutation floor. The evidence for 'real biology' is mostly self-consistency: the same SCRP generator builds the simulator backbones, the decoys and the null cohorts. Beyond that it rests on about 6 public matched tumours (Janesick: one FFPE block; Oliveira: 5 CRC patients, both checked this session) and a wet-lab experiment that has no collaborator yet. Fixes that would make it a credible NM paper: rename existence as 'co-expression beyond modelled nuisance'; add between-tumour decoys and negative controls from co-resident non-malignant cells; simulate with nuisance forms the model does not assume; anchor CNA on DNA or allele data; cut the scope to about 500–800 fully preprocessed tumours. As a 21-month, 2.5-FTE plan it is not credible.

**Fatal flaws**

- F1. THE EXISTENCE TEST USES A NON-COMPETITIVE NULL. It repeats the random-gene-null error one level up. SCRP decoys come from a per-tumour NB GLM on exactly the covariates CREDO models explicitly: depth, CNA dosage, ambient, doublet score and QC. Residuals are then permuted gene by gene. So a decoy holds (modelled nuisance) + (gene-wise marginal noise) and no co-expression at all. H0_k ('k explains no more held-out deviance in real tumours than in decoys') is therefore false for every slot that fits any real co-expression: biology, dissociation stress, pre-analytical or ischaemia effects, cohort-crossing technical axes, or residue from neighbouring programs. Decoys get PIP of about 0 for every slot, so T_k reduces to roughly 'slot k has members'. My rough estimate: a 1%-UMI program in 2,000 cells of about 3,000 UMIs gives 10^4–10^5 nats of held-out deviance per member tumour against about 0 in decoys. Every non-degenerate slot then hits the permutation floor 1/(B+1) and passes BH. The q-value cannot separate biology from unmodelled artefact, which is the user's actual complaint. The concept's line 'recurrence is necessary but not sufficient; a program must also beat decoys' is not delivered.
- F2. INTERNAL CONTRADICTION ON TUMOUR-WIDE PROGRAMS (Result 2, gap G5). The per-tumour SCRP GLM must carry a gene-specific intercept for each tumour. Decoys therefore keep each source tumour's mean profile, including the contribution m_tk W_k of any tumour-wide program. For a homogeneous (high-alpha_tk) program, real and decoy tumours give the same Delta_tk. Tumour-wide subtype programs, the advertised novelty, can then neither pass existence nor get a certified membership FDR. If tumour means are stripped from decoys instead, library-level batch effects pass as MPs. Underneath this, when tumour = sample = library (the norm in 3CA), a tumour-wide program and a library effect are not identifiable from the likelihood. beta_{t,g} gives 6,000 free parameters per tumour, so the split is set entirely by the learned Gamma(tau, tau) prior on eps. Structured technical effects that vary across tumours in many cohorts will become tumour-wide 'MPs': cold-ischaemia time, RNA quality, enzyme lot, purity-dependent ambient. Neither the cohort tier nor protocol typing can catch them.
- F3. COUNT SPLITTING DOES NOT BREAK CNA CIRCULARITY. Thinning splits noise, not signal. A positional expression program, such as a co-regulated gene cluster on one arm (keratin, MHC, HOX or histone clusters), is equally present in fold A and fold B. An expression-only caller on fold A (CopyKAT, SCEVAN, infercnv) will call the same spurious gain. The dosage offset D_ng then deletes the program on fold B. The circularity is broken only by allele-aware calls (Numbat BAF), which need BAMs, and those are mostly controlled access. The concept itself expects an expression-only fallback for much of the corpus. D_ng also needs absolute copy number CN_{t,h,seg} and ploidy P_t, which scRNA-based callers do not give: they output relative or discrete states. Whole-genome-doubled tumours (a large minority of solid tumours; exact fraction [U]) break the mapping from state to copy number. Finally, the pan-cancer gamma_g 'resource' regresses expression on an expression-derived quantity. Without DNA copy number it is circular and biased toward 1 for exactly the genes that drive the CNA calls.
- F4. VALIDATION IS LARGELY SELF-CONSISTENCY. SCRP builds three things at once: the CREDO-Sim backbones, the decoys and the null cohorts. Decoys are therefore exchangeable with program-free simulated tumours by construction. The pre-registered targets ('median <= 1 MP on null cohorts at q = 0.1', 'credible-set FDR within nominal +/- 2 points') are guaranteed by design and say nothing about real data. The planted confounders also use exactly CREDO's parametric forms: a fixed empty-droplet profile times a mixing fraction, TME-profile doublet mixtures, and multiplicative arm-level dosage taken from true clone labels. The one confounder CREDO does not model (stress/IEG) is the one it admits it cannot remove. Baselines get none of CREDO's extra inputs (ambient, doublet, CNA, protocol, QC and cohort labels), so any win confounds 'model' with 'preprocessing'. Beyond the simulator, 'real biology' rests on: n of about 6 matched public tumours (Janesick: a single FFPE breast block; Oliveira: n = 5 CRC; both confirmed by search this session); trivial perturbation predictions (CDK4/6 inhibitors lower the cell-cycle MP); and TCCA replication, which tests reproducibility rather than biology because dissociation artefacts replicate across cohorts too.

**Major concerns**

- THINNING VALIDITY UNDER HETEROGENEOUS DISPERSION (my derivation, checked by Monte Carlo). For X ~ NB(mu, size b) thinned with BetaBinomial(eps b', (1 - eps) b'): Cov(X^A, X^B) = eps(1 - eps) mu^2 (b'/b - 1) / (1 + b'). Simulated correlations: 0.07 at mu = 1, b = 2, b' = 4; 0.18 at b = 1, b' = 4; 0.17 at mu = 5, b = 2, b' = 4. CREDO uses one phi per gene and protocol across about 3,000 tumours. Every tumour that is noisier than average (true b_t < phi_hat) therefore has positive fold dependence, which makes tests anti-conservative, and these are the low-quality tumours. PIP-weighted Stouffer amplifies a per-tumour Z bias delta to delta * sqrt(T): delta = 0.1 over 500 member tumours gives a bias of 2.2, enough to flood gene sets. The rank of the pilot low-rank fit sets the sign of the bias. The decoy-gene check is blind to this failure. A gene permuted within a tumour has no low-rank structure, so its conditional dispersion equals its marginal dispersion and is easy to estimate. That is exactly not the Hivert regime that real program genes are in.
- GENE FDR TESTS A NULL THAT IS NEVER TRUE. The null is beta_gk = 0, the counts are 10^6–10^7, and evidence is combined over hundreds of member tumours. Almost every gene correlated with a program at 1–5% fold change will be 'significant'. Core sets for cell cycle or IFN will run to thousands of genes and be useless for scoring. A minimum-effect null is needed (TREAT-like, e.g. fold change >= 1.25 [method citation U]), as is a random-effects combination across tumours rather than fixed-effect Stouffer.
- MEMBERSHIP FDR: EXCHANGEABILITY AND GRANULARITY. (a) The target-decoy estimator assumes non-member real tumours score like decoys. Non-members that carry a related, nested or partially overlapping program (the planted pairs at cosine 0.4–0.6; EMT-I to EMT-VI) are 'near-miss' targets with higher PIPs than any decoy, so FDR is underestimated exactly where it matters. (b) With a 15% decoy fraction, n_real / n_dec is about 6.67. One decoy exceedance gives FDRhat = 6.67 / #members, so a credible set needs >= 67 members to tolerate a single decoy at 0.1 (my calculation; the same at 3,000, 300 or 185 tumours). Rare MPs (3–10 tumours) are certified only when zero decoys exceed, which gives no finite-sample guarantee without a +1 correction and becomes impossible with one. (c) Decoys share lineage and cohort nodes, so they deflate pi_{k,l} by about 15%, roughly a -0.3 logit shift in prior PIP at p = 0.5 (my calculation). (d) LBF_tk is a tempered ELBO difference, not a Bayes factor, and the ELBO gap differs between z = 1 and z = 0. Calling it an 'exact local Bayes factor' is incorrect.
- THE POWER POSTERIOR IS CALIBRATED AGAINST THE WRONG FAILURE MODE. zeta is meant to suppress systematic, cross-tumour misfit, but SCRP destroys co-expression, so decoys carry only gene-wise marginal misfit. Decoys will rarely be recruited at any zeta, and the 'largest zeta with decoy FDR <= 0.1' rule will pick zeta near 1, making tempering a no-op; pilot P4 half-anticipates this. zeta chosen on a 500-tumour subsample does not carry over to 3,000 tumours. Under the coarsened-posterior view (Miller and Dunson, JASA 114:1113–1125, 2019 [V this session]) the needed tempering scales roughly as 1/n. The KL term is not tempered, so tempering acts per cell on cell-local usages. Components carried by tumour-level parameters (m_tk, eps_t) cost O(1) nats per tumour but gain zeta * min(n_t, 2000) * (per-cell gain). Tumour-wide junk components are therefore cheap at any zeta on the grid {1 ... 0.02}. Please derive the effective tempering at the cell, tumour and global levels.
- IDENTIFIABILITY CLAIMS DO NOT COVER THE ACTUAL MODEL. Proposition 1 (SSC is monotone under row union) is essentially trivial. It assumes a common Theta W with no tumour-specific additive terms, whereas CREDO gives each tumour 6,000 free baseline parameters plus 4 x 6,000 private loadings plus cohort slots, any of which can absorb parts of W_k. Proposition 2 is unproven. Its anchor condition (C1), 'an exclusive anchor gene among 200 global programs', is implausible for overlapping programs: G1/S vs G2/M, IFN-I vs IFN-II/MHC-II, the EMT variants. The 'cone-exclusion' prior is pairwise cos^2, not cone exclusion. A private slot equal to an equal mixture of 4 shared programs with disjoint supports has cos^2 = 0.25 to each and is barely penalised. The penalty is also weighted by z_tk, so a non-member tumour can host W_k in a private slot at no cost. That makes non-membership self-reinforcing (the 0.4%-usage tumour gets z = 0, and V_tj absorbs W_k) and contradicts the theory's own condition C3. Measured against about 12,000 private and about 600 cohort vectors in 6,000 dimensions, the 'identifiability margin' (smallest angle to the cone of all other programs) will be near zero or uninformative. Measured against global programs only, it ignores the real competitors.
- COHORT, LINEAGE AND PROTOCOL CONFOUNDING. Most cohorts (one study x one protocol) contain a single cancer type. For lineages with 1–2 cohorts, three things cannot be separated: e_{c,g} (cohort gene efficiency), b_{l,g} (lineage baseline), and lineage-restricted tumour-wide programs. Protocol is nested in lab and lineage: sn mostly in brain and frozen biobank tissue, Smart-seq2 mostly in early Tirosh-lab studies. a_{p,g} and eta_{k,p} are identified only from the few tumours profiled under >= 2 protocols, likely well under 1% of the corpus and in 2–3 cancer types. Protocol typing for MPs outside those lineages will be driven by the prior. The 'normal-shared' probability is also confounded, because normal donors come from different studies and protocols.
- THE THRESHOLDS COME BACK UNDER NEW NAMES. 'The only user-facing constants are delta_min and error rates' is not accurate. The MP set also depends on: K = 200, alpha0, c, n0 = 2000, lambda_ce, the tau prior, Beta(5,5), concentration 50, Gamma(0.3,0.3), J = 4, R = 4, the Hard-Concrete schedule and the zeta grid. Above all it depends on the tier-move rules: R_k < 1.5; >= 3 tumours in >= 2 cohorts at cosine > 0.7; > 90% of members in one cohort; merge at cosine > 0.9; split on bimodality. These are 35/10/10-style overlap rules, now applied inside inference. Reviewers will ask for the same sensitivity analysis CREDO demands of Gavish: MP count and identity against every one of these constants.
- AMORTISATION IS ILL-POSED FOR PRIVATE AND COHORT SLOTS, AND THE HMC CHECK IS NOT CREDIBLE. The encoder never sees tumour identity, yet it must output usages for J = 4 private and R = 4 cohort slots whose meaning differs in every tumour and cohort (V_{t,j}, U_{c,r}). A shared MLP cannot map a cell to 'private slot 2 of tumour t' without seeing V_t. 'HMC on a 20-tumour subset' would mean millions of cell-level latents, discrete z (which HMC cannot sample without marginalising over 2^K configurations), and a multimodal NMF posterior with label switching. It will not mix and cannot serve as a gold standard for PIPs.
- ONE W_k ACROSS ALL LINEAGES. W_k has no lineage-specific loading perturbation. A program whose genes differ by lineage (EMT in glioma vs carcinoma, hypoxia, lineage-flavoured IFN) must either split into several global slots, which the tier moves may later merge at cosine > 0.9, or leak into cohort and private tiers, which cone exclusion forbids for members. MP count and identity will then depend on the tier-move thresholds, reintroducing the instability the method claims to remove.
- DEFINITION AND NUISANCE PRIORS. (a) delta_min applies to tumour-mean usage (0.5% of UMIs). That excludes rare-subpopulation programs: 2% of cells at 10% usage gives a 0.2% mean. The truncation pushes below-threshold usage into private and baseline terms. A dual criterion is needed (fraction of cells with usage > x). (b) The doublet mixing prior omega ~ Beta(5,5) puts only 2.0% of its mass below 0.2 (my calculation). Malignant–lymphocyte doublets usually have a small TME fraction because malignant cells carry more RNA. Whole TME cells misclassified as malignant (delta about 1), the likely source of the 3CA v2 macrophage, T-cell and endothelial MPs, are not represented. (c) The QC tier will absorb real programs correlated with %mito or complexity, such as hypoxia, OXPHOS and cycling.
- TYPING AND 'REAL BIOLOGY' VALIDATION ARE UNDERPOWERED OR MISSPECIFIED. (a) 'Protocol-sensitive' mixes up dissociation artefact with measurement compartment: sn loses cytoplasmic and mitochondrial transcripts, and Flex/Xenium are probe-limited. Real cytoplasmic programs will be typed as artefacts. (b) The prediction 'protocol-sensitive MPs absent in CCLE' is wrong: cell lines are trypsinised at 37 C, and the 2026 catalogue already contains 'Stress (in vitro)' [V+ per landscape]. (c) The AUROC on 11 hand-removed MPs has a 95% CI of about +/- 0.12 at AUROC 0.9 and +/- 0.16 at 0.8 (Hanley–McNeil, my calculation). 7 of the 11 are ribo/mito LQ, which gene content alone detects. CREDO never produces these MPs because nuisance terms absorb them, so it is undefined how its evidence vector is computed for an external gene list. The rule designers have also seen the labels. (d) In the prospective study, 10 paired specimens detect only d_z >= 1.0 at alpha = 0.05 unadjusted, or >= 1.8 after Bonferroni over about 40 MPs (my calculation). Pre-registered IEG and heat-shock predictions are not risky because they are already known [dissociation literature U]. (e) Spatial niches (hypoxia near necrosis) and Tahoe predictions (CDK4/6 inhibitors lower cell cycle) confirm MPs nobody disputes.
- DATA AVAILABILITY IS UNDERESTIMATED. The ambient term A_s needs unfiltered droplet matrices, which most GEO submissions do not deposit. The fallback (ambient estimated from cell profiles, DecontX-style) is circular with the TME profiles M_{s,e}. Allele-aware CNA needs BAMs, mostly on dbGaP or EGA: there would be dozens of data-access-committee (DAC) applications, each typically taking months, so a data audit in M0–3 is not achievable. Smart-seq2 read counts (length bias, amplification) are poorly described by NB even with their own dispersion. The 10^7 malignant-cell target looks inflated: ITHmapper assembled 3.6M cells over more than 1,000 tumours and TCCA 1.09M malignant cells over 36 studies (per landscape [V]). Realistic public malignant cells with raw UMIs are likely 2–5 x 10^6 [estimate]. Controlled-access consent may also restrict release of per-cell outputs in the reference bundle.
- SCALABILITY AND COST. Compute per step is plausible: 16,384 x 6,000 x 220 = 2.2 x 10^10 MAC, and 610 steps per epoch. Missing or optimistic: (a) the number of epochs is never stated; (b) I/O is about 120 GB per epoch (value plus index), read as 16,384 random rows from 64 shards per step, and random row access through backed AnnData/h5py is orders of magnitude slower than sequential, so a loader benchmark is needed; (c) evaluating tier moves (O(K^2) merge candidates, re-optimisation for each proposal) is not costed; (d) leave-one-cohort-out over about 130–160 cohorts at 6–24 GPU-h per fit is 0.8k–3.8k GPU-h by itself, which uses up the whole 1,500–4,000 GPU-h calibration budget before the zeta grid, seeds, null cohorts, add/remove-study perturbations, ablations and GPU baselines.
- TEAM AND TIMELINE ARE ABOUT 2x UNDER-RESOURCED. 1 stat-ML postdoc, 1 computational biologist and 0.5 RSE over 21 months must deliver all of the following: a roughly 10-module Pyro model with novel inference (tier moves, Hard-Concrete, semi-amortisation, LBFs); two theorems plus a thinning extension; harmonisation of about 3,000 tumours with CellBender, scDblFinder and Numbat (BAM access and phasing); a new simulator; about 17 fairly tuned baselines, several in R at scale; six validation axes; a prospective wet-lab study with no collaborator yet and IRB lead time; an R bridge; a versioned reference. A realistic plan is 36–48 months, or the same team with the scope cut to the definition, calibration and typing on 500–800 tumours.
- DOWNSTREAM BIOLOGY CLAIMS ARE CONFOUNDED. The posterior co-membership network ('repertoires, exclusivity') will be driven by lineage, cohort, purity, depth and tumour cell number. Membership power grows with n_t through w_t and zeta, so lineages profiled with few cells (Smart-seq2) will look 'MP-poor'. The network needs conditional analyses and per-tumour minimum-detectable-usage reporting.

**Prior-art threats**

- GBCD: Liu, Carbonetto, Willwerscheid, Oakes, Macleod, Stephens, Nat Genet 57:263-273 (2025). Label-free separation of shared vs patient-specific programs; directly pre-empts the 'shared vs private' framing [V].
- Grabski, Trippa, Parmigiani, Bayesian multi-study NMF, Genome Biol 26:98 (2025). Count-native NMF with a study-by-factor inclusion matrix, learned K and covariates; the statistical core of CREDO's subset sharing [V].
- Tetris: Grabski, De Vito, Trippa, Parmigiani, AoAS 17(3) (2023), with an IBP over factor-sharing patterns. Also APAFA (Bortolato and Canale, JASA 2026) and SUFA (Chandra, Dunson, Xu, JASA 2025): subset-shared, private and identifiable multi-study factor models [V].
- BayesPowerNMF: Xue, Miller, Carter, Huggins, PLoS Comput Biol 22(6):e1014372 (2026). Power posterior for NMF against misspecification [V]. Coarsened posterior: Miller and Dunson, JASA 114:1113-1125 (2019) [V this session]. Together they cover the zeta component and imply zeta should scale with n.
- Hong, Klöckner, Mazzilli, Spira, Beane, Yajima, Campbell, 'An atlas of gene programs in non-small cell lung cancer...' bioRxiv 10.1101/2025.05.06.652317 (2025), GPAS framework. Celda count-based modules across 127 NSCLCs from 12 studies, with linear mixed-effect decomposition of technical vs biological variability. Not in the landscape; partly pre-empts the 'nuisance decomposition and typing across studies' angle [V title/authors/abstract-level via search snippet; method details P].
- ClusterDE synthetic-null FDR (Song, Li et al., bioRxiv 2023) [V]; Neufeld et al. count splitting (Biostatistics 2024) and NB count splitting (arXiv 2307.12985) [V]; Hivert et al., arXiv 2405.13591, limits of thinning under mixtures [V]. These cover the calibration toolkit, so the novelty is the adaptation only. Hivert directly undermines the gene-FDR claim unless tumour-specific dispersion is handled.
- scDEF (Ferreira, Kuipers, Beerenwinkel; bioRxiv): multi-level Gamma-Poisson factorisation with shared and patient states [V]. Consensus scHPF (Levitin...Sims, bioRxiv 2023) [V].
- ITHmapper (Nofech-Mozes...Abelson, bioRxiv 2026): pan-cancer, more than 1,000 tumours, a mapping tool, Tahoe and Visium HD validation [V]. TCAT/starCAT (Kotliar...Raychaudhuri, Nat Methods 22:1964-1980, 2025): catalogue plus projector template [V].
- sciRED (Nat Commun 2025): Poisson-GLM regression of known covariates, then factors [V]. DecontX (Genome Biol 2020): Bayesian contamination mixture with cluster profiles [V per landscape]. CellBender (Nat Methods 2023) [V]. These are precedents for CREDO's nuisance terms.
- KINOMO (bioRxiv 2022; AACR 2024 abstract): KL-NMF with kernel correction for genomic variability between patients. Must be acknowledged on the CNA angle [V exists; details not checked].
- Weine, Carbonetto, Irizarry, Stephens, Poisson 'log1p NMF', arXiv 2601.05845. The GBCD group is moving to count models, a scoop risk for the 'count-native GBCD' differentiator [V title/authors only].
- ZINB-WaVE (Risso et al., Nat Commun 2018) and GLM-PCA (Townes et al., Genome Biol 2019): NB/ZINB count factor models with cell and gene covariates, precedent for 'NB factorisation with nuisance covariates' [U, from memory].
- Knockoff filter (Barber and Candès, Ann Stat 2015) and target-decoy FDR (Elias and Gygi, Nat Methods 2007): methodological precedent for decoy-calibrated FDR. Their exchangeability requirements are exactly what CREDO's decoys violate for near-miss non-members [U, from memory].
- Dissociation-stress literature (van den Brink 2017; O'Flanagan 2019; Denisenko 2020; Slyper 2020). Makes Result 1(a), the split of stress MPs into an ex vivo IEG part and an in vivo hypoxia part, largely expected rather than new [U].

**Required changes**

- Reframe and strengthen existence. (a) Rename the current test 'co-expression beyond modelled nuisance' and report the full distribution of q_k; if most slots sit at the permutation floor, say so. (b) Add competitive nulls that carry unmodelled technical co-expression. The cheapest and most general is a WITHIN-SAMPLE NEGATIVE CONTROL from co-resident non-malignant cells: a program whose per-sample usage co-varies in malignant and TME cells of the same sample (IEG stress, ambient, ischaemia) is flagged as technical or sample-level. This is available in every 10x sample, unlike cross-protocol data. (c) Add a cross-cohort held-out replication test: fit on cohorts A, score on B, with cohort-label permutation.
- Use two levels of decoys. Keep SCRP (within-tumour) for ITH programs. Add between-tumour decoys, e.g. tumour-mean deviations from the lineage x cohort mean permuted gene by gene across tumours within strata, for tumour-wide programs. Limit claims about tumour-wide subtype programs to evidence from tumours with >= 2 libraries or regions, plus cross-cohort replication, and state plainly that tumour-wide biology and library batch cannot be separated when tumour = library.
- Membership FDR. Use the +1 (knockoff+/target-decoy) finite-sample form. Raise the decoy fraction (>= 50%), or repeat decoy draws and report Monte Carlo variability. Report the minimum certifiable credible-set size per MP. Show calibration in simulation specifically for near-miss non-members (nested or overlapping programs at cosine 0.4–0.6). Keep decoys out of the lineage nodes used to estimate pi_{k,l}, or correct the prior.
- Gene FDR. Test against a minimum effect (e.g. H0: fold change <= 1.25) rather than zero. Estimate thinning dispersion per tumour or cohort, not per protocol. Check fold independence empirically by repeated thinning on real data, using the covariance formula Cov = eps(1 - eps) mu^2 (b'/b - 1)/(1 + b'). Replace permuted decoy genes with null genes simulated from the fitted model (W_kg = 0, matched mean, dispersion and co-loading) so the check probes the Hivert regime. Combine tumours by random effects instead of fixed-effect Stouffer, to avoid sqrt(T) bias amplification.
- CNA. Drop the claim that fold-A calling breaks circularity. Apply the dosage offset only to allele-supported (Numbat BAF) segments. For expression-only calls, use CNA for typing, never removal. Model relative states with per-tumour dosage slopes instead of absolute CN / ploidy, and handle WGD. Validate on tumours with matched scDNA or WGS. Present gamma_g as a resource only where DNA copy number is measured.
- Make the simulation independent of the model. Build backbones with at least one generator other than SCRP (e.g. scDesign3 or splatPop fitted to real tumours). Plant confounders in forms CREDO does NOT assume: nonlinear or cell-size-dependent ambient; doublet fractions set by RNA content; misclassified whole TME cells; focal CNAs and tumour-varying dosage compensation; CNA calls re-inferred from the simulated expression rather than true clone labels; stress that crosses cohorts; tumour-specific dispersion. Report recovery at 1–5% effect sizes, not only >= 5%. Give every baseline the same preprocessing (CellBender, doublet removal, CNA regression) as a separate arm.
- Power posterior. Derive the effective tempering at the cell, tumour and global levels. Scale zeta with effective n (coarsened-posterior style) and choose it on the full data by held-out deviance. Calibrate it with decoys that carry systematic shared misfit, e.g. data simulated from the fitted NB model plus an injected cross-tumour misfit component, not with SCRP. Otherwise drop it from the headline.
- Identifiability. State conditions for the full tiered model, including the baseline, private, cohort and nuisance tiers. Prove Proposition 2 or demote it to a conjecture backed by simulation sweeps over anchor availability and private-cone geometry. Implement true cone exclusion (NNLS-projection residual of V_tj onto cone{W}) applied to ALL tumours, not weighted by z_tk. Compute the identifiability margin against the actual competitor set.
- Inference. Infer private and cohort usages by local optimisation, not amortisation, or condition the encoder on the tumour's dictionary (e.g. attention over V_t, U_c). Report the amortisation gap. Replace the HMC check with exact Gibbs or SMC on small synthetic problems (<= 5 tumours, K <= 20).
- Let loadings vary by lineage, e.g. hierarchically W_{k,l} = W_k * exp(xi_{k,l}) with shrinkage, or show that MP identities are insensitive to this.
- Publish a sensitivity table for every constant: K, alpha0, c, n0, lambda_ce, tau prior, delta_min, the Beta and Gamma hyperparameters, J, R, and each tier-move threshold. Report MP count, matched identities and type labels, exactly as the concept demands of the Gavish code.
- Before M3, run a data audit with explicit go/no-go: the fraction of tumours with raw UMIs, with unfiltered droplet matrices, and with obtainable BAMs (including expected DAC timelines); a protocol x lineage x cohort table of multi-protocol tumours, before claiming eta_{k,p}; and realistic malignant-cell counts.
- 'Real biology' validation. (a) Power the prospective split-specimen study (n = 10 detects only d_z >= 1.0 unadjusted), or pool cancer types with stated effect-size targets. (b) Pre-register RISKY predictions (EMT/pEMT, hypoxia/UPR, IFN, cancer-type-specific MPs), not IEG/heat shock. (c) Separate dissociation effects (warm vs cold protease, fixed before dissociation) from compartment effects (sn) in the typing rules. (d) Drop the CCLE-absence prediction. (e) Add perturbation tests of contested MPs in organoids or PDX. (f) Give AUROC with CIs, define how the artefact AUROC is computed for 3CA gene lists, and have the typing rules frozen by someone blind to the 11 labels.
- Scope and resources. Either cut the first paper to the definition, calibration, typing and benchmark on about 500–800 fully preprocessed tumours (CanSig + a 3CA subset + multi-protocol and multi-region datasets) and move the atlas audit to a follow-up, or add about 2 FTE (data engineering, wet-lab coordination) and about 12 months. Budget leave-one-cohort-out and add/remove-study refits explicitly, at stated scale.

**Ideas worth keeping**

- Defining an MP as a pair (W_k, T_k), with tumour-level inclusion z_tk and PIPs. This makes 'which tumours use this program' a first-class, testable output. It is the real definitional advance over top-50 lists and over GBCD.
- Three tiers (global / cohort-private / tumour-private), with the cohort tier explicitly absorbing lab- and study-specific artefacts instead of promoting them to MPs.
- A count likelihood on raw UMIs with fold-change loadings relative to each gene's pooled baseline and no centring, plus an ITH index that places each MP between 'within-tumour ITH program' and 'tumour-wide subtype program'.
- Decoy tumours carried blind through the same fit. Kept as a DIAGNOSTIC of how incomplete the nuisance model is (decoy recruitment per MP), not as the existence test, this is useful and new.
- Count-split, out-of-fold gene effect sizes and variable-length gene sets with q-values, once fixed with tumour-specific dispersion and a minimum-effect null.
- Turning protocol typing into a falsifiable, pre-registered prediction scored on four-arm split specimens (warm vs cold protease, fixed-then-dissociated Flex, sn).
- Always reporting the curve of MP count against tempering, stability across seeds, thinning folds and leave-one-cohort-out, and per-MP identifiability diagnostics (anchor genes, near-pure cells).
- Using the 11 hand-removed 3CA v2 MPs, plus Gavish's removed LQ, doublet and single-study MPs, as a labelled artefact set (with CIs and blind rules).
- The cheap kill-switch pilots P1 (false MPs from the standard), P3 (CNA/positional audit), P6 (data audit) and P7 (artefact predictability). They should run before any model code is written.
- A projector with novelty slots and protocol-aware panel transfer, as a deliverable rather than the headline.


---

## 2. CALIPER (CALibrated Inference for Programs Expressed Recurrently). Working name; check for collisions with existing software before release.  
*Design angle: `inference-layer`*

### Claim

CALIPER is an engine-agnostic statistical layer that turns the output of any gene-program factorisation (Gavish NMF, cNMF, GBCD, DRVI, scVI+Hotspot) into calibrated metaprogram inference. Using molecule-level count thinning with provably sign-controlled overdispersion and cross-fold co-expression tests, it reports for every metaprogram a lower confidence bound on how many tumours and cohorts express it, FDR-controlled tumour and gene membership, and a reproducible type (technical, copy-number-driven, clone-linked, protocol-sensitive or shared non-genetic biology) that replaces manual curation. Applied to the pan-cancer 3CA corpus, it yields an audited, engine-robust metaprogram atlas.

### Formal definition of an MP

PRINCIPLE. The metaprogram (MP) is the ESTIMAND. Any factorisation (Gavish NMF, cNMF, scHPF, GeneNMF, GBCD, DRVI, scVI+Hotspot) is an ESTIMATOR of its location. CALIPER supplies inference about the estimand that stays valid whatever estimator is plugged in.

SETTING
- Tumours t = 1..T, nested in cohorts c(t) (study x chemistry), each with a protocol p(t) (fresh 3' sc, sn, Flex, ...).
- Malignant cells i = 1..n_t; genes j in a universe G; UMI counts X_tij.
- Working observation model: X_tij | mu_tij ~ NB(mean mu_tij, size b_tj). Entries are conditionally independent given the latent means mu: cells are independent, and genes are independent given each cell's state.

NUISANCE LADDER
Nested working models for the conditional mean nu^N; each rung adds one measured nuisance class.
- N0, depth: log nu_tij = log l_ti + alpha_tj.
- N1, exogenous technical: N0 plus ambient RNA (per-cell fraction rho_ti times the sample's empty-droplet profile a_t) and heterotypic doublets (per-cell doublet probability pi_ti times d_t, the mean profile of the same sample's TME cells). nu = l_ti[(1 - rho_ti - pi_ti) exp(alpha_tj) + rho_ti a_tj + pi_ti d_tj].
- N2, QC-entangled: N1 plus log complexity, %mito and %ribo as log-linear covariates.
- N3, cis-dosage: N2 plus gamma_tj * log2(CN_tij/2). CN_tij is per-cell, per-gene copy number from allele-aware Numbat, or from an expression caller run with a 2 Mb window around gene j left out.
- N4, clone: N3 plus per-tumour clone indicators (trans effects of genotype).
Cell cycle and all other biology are never treated as nuisance.

EXCESS CO-EXPRESSION
- Pearson residuals: e^N_tij = (X_tij - nu^N_tij) / sqrt(nu^N_tij + (nu^N_tij)^2 / b_tj).
- Omega^N_t = off-diagonal part of the within-tumour cross-gene covariance of residuals, E_i[e^N_ti e^N_ti^T].
- Omega^N_t = 0 exactly when tumour t has no co-expression beyond nuisance class N. Call this H0^N_t.

TUMOUR EXPRESSION OF A PROGRAM
- For a non-negative unit vector m in S_+^{G-1}: kappa^N_t(m) = m^T Omega^N_t m = sum over j != j' of m_j m_j' Omega^N_t[j, j'].
- Tumour t expresses m beyond N iff kappa^N_t(m) > 0.
- Specific version kappa^N_t(m | M): the same quantity after residualising on the activity scores of a set M of other programs (parent and siblings). This is what decides merge versus split.

LOCATION (engine-indexed, stated explicitly)
- Engine E applied to tumour t at infinite depth has population programs W^E_t = {w_tk}.
- Tumour-coverage function: F_delta(m) = sum_t 1[ max_k cos(w_tk, m) >= 1 - delta ]. It counts TUMOURS, not programs.
- MP locations are the modes of F_delta, i.e. the connected components of its upper level sets {F_delta >= rho}. Varying rho gives a cluster tree, so parent/child MPs (for example EMT-I..VI) are part of the estimand rather than arbitrary choices.
- For pooled engines (GBCD, DRVI, scVI+Hotspot), the locations are the engine's global programs.

DEFINITION
A metaprogram at nuisance level N, recurrence r_min and cohort spread s_min is a tuple (m, S, G, tau):
- m: a location as defined above.
- S^N(m) = {t : kappa^N_t(m | siblings) > 0}: the tumours that express it.
- EXISTENCE: |S^N(m)| >= r_min AND S^N(m) spans >= s_min distinct cohorts. Defaults: r_min = 5, s_min = 2.
- G^N(m) = {j : eta^N_j(m) > eta_0}, where eta^N_j(m) is the across-tumour (random-effects) mean of the nuisance-adjusted within-tumour log fold change of gene j per SD of MP activity. eta_0 = log 1.2 by default; this minimum-effect null stops atlas-scale n from making every gene significant.
- tau(m): the TYPE, set by the first rung at which existence fails.
  - Exists at N0 but not N1: technical (ambient, doublet or contamination).
  - Fails only at N2: QC-entangled. This is a flag, not a deletion.
  - Exists at N2 but not N3: cis-dosage (CNA).
  - Exists at N3 but not N4: clone-linked.
  - Exists at N4: beyond-nuisance, sub-typed as protocol-sensitive (tumour membership depends on protocol), single-cohort (fails s_min), or shared non-genetic biology (pan-cancer vs lineage-restricted, by cancer-type entropy).

WHAT IS GUARANTEED
- A test of H0: |S^N(m)| < r_min (existence), with FDR control across MPs.
- A lower confidence bound r_L(m) with P(r_L > |S^N(m)|) <= alpha.
- FDR-controlled estimates of S(m) and G(m).
- Error control at each rung of the type ladder.
- All inference is conditional on the observed cells and on hypotheses selected from statistically independent molecule folds.

SCOPE. The estimand is within-tumour heterogeneity (Gavish's ITH programs). Programs constant across a whole tumour (between-tumour baselines) are outside the core and covered only by an optional pseudobulk extension.

### Model and math

(1) MULTI-FOLD NB THINNING (molecule-level independence)
- For each tumour and gene: (X^A, X^B, X^C) | X ~ DirichletMultinomial(X; eps_A b~, eps_B b~, eps_C b~), with (eps_A, eps_B, eps_C) = (0.5, 0.2, 0.3).
- Fold roles: A = discovery (engine and aggregation); B = selection (rank, filters, delta); C = final inference only.
- If b~ = b, then X^k ~ NB(eps_k mu, eps_k b) and the folds are mutually independent (Neufeld et al. data thinning, JMLR 2024 [V]; NB count splitting, arXiv 2307.12985 [V]).

THEOREM 1 (exact fold covariance under misspecified dispersion; derived here)
- For k != l: Cov(X^k, X^l | mu) = eps_k eps_l mu^2 (b~ - b) / (b (b~ + 1)).
- Proof: E[X^k X^l | X] = X(X-1) eps_k eps_l b~/(b~+1) and E[X^2] - mu = mu^2 (1+b)/b, so Cov = eps_k eps_l mu^2 [ (1+b) b~ / (b (b~+1)) - 1 ].
- Monte Carlo check in this session (empirical vs formula):
  - (mu, b, b~) = (10, 10, 2): -3.97 vs -4.00
  - (10, 10, inf): 1.499 vs 1.500
  - (20, 5, 50): 10.63 vs 10.59
- Consequences:
  - Poisson or binomial splitting of NB data (b~ = inf) gives covariance +eps_k eps_l mu^2 / b, which is anti-conservative.
  - b~ < b gives negative covariance.

COROLLARY 1 (sign control for any latent structure)
- Claim: for ANY latent structure of mu across a tumour's cells (mixture, continuum or factor model), the population NB size fitted while ignoring that structure satisfies b_marg <= b.
- Proof: m^2 / b_marg = E[mu^2]/b + Var(mu) >= m^2/b, where m = E mu.
- The same bound holds for an NB fitted conditionally on an under-fitted model of rank K < K0 (law of total variance).
- Therefore thinning with marginal or under-fitted dispersion makes every cross-fold covariance <= 0. Over-fitted dispersion does not.
- This answers Hivert et al. (arXiv 2405.13591 [V]) for one-sided statistics.
- Toy check: a gamma-heterogeneous mu gave b_marg = 1.54 against b = 10.

DISPERSION ESTIMATOR ("monotone ascent")
- Pass 1: thin with b~ = the lower 10% profile-likelihood bound of the marginal within-tumour size; select K1.
- Pass 2: re-estimate b~ conditional on the rank-K1 NB fit (still under-fitted with high probability); re-thin; select K2 >= K1.
- Repeat until K is stable, at most 3 passes.

DESIGN RULE: no statistic ever pairs a fold-A quantity with a fold-C quantity from the SAME matrix entry. Only off-diagonal (j != j') or leave-gene-out statistics are used, so dispersion error can enter only through the engine's dependence on single entries.

THEOREM 2 (to prove; first-order sketch)
- Let L_t = max |d m~ / d X^A_tij| be the engine's per-entry influence.
- Then |E CF_t| <= C · L_t · |G| · max_j |Cov(X^A_tj, X^C_tj)|.
- The bias has sign <= 0 when b~ <= b and the loadings are non-decreasing in own-gene counts.
- For averaging or consensus engines L_t = O(1/n_t), so the standardised bias is O(n_t^-1/2).

IN-SESSION TOY CHECKS (my own simulations, not published)
Setup: 300 cells, 400 genes, NB size 3; the "engine" is a PCA that keeps its top 30 genes from fold A.
- Null tumours, one-sided alpha = 0.05, CF statistic from step (5):

  | Dispersion used for thinning | Type I | Reps |
  |---|---|---|
  | Exact b | 0.042 | 500 |
  | Estimated marginal b | 0.032 | 500 |
  | 0.8 x marginal b | 0.030 | 500 |
  | Poisson splitting | 0.993 | 150 |

- A within-fold-C U-statistic on the same m gave type I 0.087 with b~ = b/3 and 0.42 with Poisson splitting. Within-fold co-expression tests are sign-blind and are therefore never used.
- Power with a planted 40-gene program (fold change 1.5 in 30% of cells): 0.39 with exact conditional b vs 0.26 with marginal b. Both reach 1.0 at fold change 2. This motivates the two-pass estimator and pooling across tumours.

(2) NUISANCE GLMs
- Fitted in each fold separately, with offset log(eps_k l_ti), by batched IRLS (glmGamPoi-style).
- Residuals: r^{k,N}_tij = (x^k - nu^k) / sqrt(nu^k + (nu^k)^2 / (eps_k b~)), clipped at ±sqrt(n_t).

(3) ENGINE ADAPTER AND RANK
- An engine maps fold-A counts to loadings W^A_t (G x K) and usages Theta^A_t (n_t x K).
- Universal count-scale refit: for each gene, fit an NB GLM on fold A with log mu^A_ij = log(eps_A l_i) + nuisance + Theta^A_i beta_j, then predict mu^B = (eps_B/eps_A) mu^A.
- Held-out deviance: D_t(K) = 2 sum_ij [ x^B log(x^B/mu^B) - (x^B + eps_B b~) log((x^B + eps_B b~) / (mu^B + eps_B b~)) ].
- K_t = the smallest K within 1 SE (over R thinning replicates) of the minimum deviance.
- Secondary criterion, usable for any engine: cross-fold replicability (fit on A and on B, Hungarian matching, require cos >= 1 - delta).
- Program filter: CF statistic on folds A x B under N0, BH at q = 0.1 within each tumour. This replaces intra_min = 35.

(4) AGGREGATION (the location estimator)
- Nodes are the surviving fold-A programs, as cNMF-style z-scored spectra; distance d = 1 - cos.
- Tumour-distinct core distance: core_rho(v) = distance to the rho-th nearest program from a DIFFERENT tumour (rho = 3).
- Mutual-reachability HDBSCAN-type cluster tree, with one representative per tumour per cluster.
- Consensus location m_C = normalised elementwise median of member loadings over genes with positive median.
- delta = 1 - (5th percentile of the cosine between a tumour program and itself re-estimated on an independent thinning replicate). Two programs count as "the same" if they are as similar as one program is to itself on independent molecules.
- There is no founder order and no 10/50 overlap rule.

(5) TUMOUR-LEVEL TEST: CROSS-FOLD CO-EXPRESSION
- Fold scores: s^k_ti = sum_j m~_j r^k_tij.
- Cell term: h_ti = s^A_ti s^C_ti - sum_j m~_j^2 r^A_tij r^C_tij = sum over j != j' of m~_j m~_j' r^A_tij r^C_tij'.
- CF_t = mean_i h_ti; Z_t = sqrt(n_t) CF_t / sd_i(h_ti).
- Under H0^N_t with exact thinning, E[h] = 0 and Z is asymptotically N(0,1).
- Empirical null for residual nuisance misfit (Efron-style [U]): compute Z for 1,000 pseudo-programs that place m's weight profile on fold-A "null genes" of tumour t (bottom half of max loading and in no program significant on fold B). Genes are matched on expression decile, plus ambient-profile weight from rung N1 upward and chromosome arm from rung N3 upward.
- Standardise: Z*_t = (Z_t - mu0) / sigma0; p^N_t(m) = 1 - Phi(Z*_t).

(6) RECURRENCE BY PARTIAL CONJUNCTION (Benjamini & Heller 2008 [U])
- Sort the tumour p-values p_(1) <= ... <= p_(T).
- p^{r/T}(m) = min over i = r..T of (T - r + 1) p_(i) / (i - r + 1).
- Cohort level: p_c = Simes combination over the tumours in cohort c; then p^{s/C} is computed analogously across cohorts.
- Existence: p_exist = max(p^{r_min/T}, p^{s_min/C}) (intersection-union test).
- Lower bound: r_L(m) = max{ r : p^{r/T} <= alpha }. It is valid because the hypotheses H^{r/T} are nested and tested in a fixed sequence.
- Tumour membership: BH over {p_t} at q_T = 0.1, plus a closed-testing true-discovery bound (Goeman & Solari [U]).
- Across MPs: BH by default; e-BH in strict mode (Wang & Ramdas [U]).
- Split/merge: a child MP is kept only if its specific existence, kappa(child | parent, siblings), is rejected.

(7) GENE MEMBERSHIP
- Leave-gene-out activity: u^{-j}_ti = standardise(s^A_ti - m~_j r^A_tij), residualised on the nuisance design.
- One-step NB score test on fold C:
  - S_tj = sum_i u_ti (x^C - nu^C) / (1 + nu^C / (eps_C b~))
  - I_tj = sum_i u_ti^2 nu^C / (1 + nu^C / (eps_C b~))
  - beta_tj = S/I, with a sandwich SE.
- Random-effects pooling over t in S(m) (DerSimonian-Laird or REML [U]) gives eta_j, SE_j, tau_j^2 and I^2_j.
- Minimum-effect statistic z_j = (eta_j - eta_0) / SE_j (TREAT-style [U]); BH within each MP at q = 0.05.
- Gene labels: core (I^2 < 50%), variable, and anchor (significant partial effect conditional on all other MPs and on no other MP).

(8) DERANDOMISATION AND STABILITY
- R = 5 independent thinnings. MPs are matched to consensus labels by the Hungarian method with cos >= 1 - delta.
- Per hypothesis, replicate results are combined by a Cauchy combination of p-values (default [U]), or in strict mode by averaging e-values e = 0.5 p^{-1/2} (an unmatched replicate contributes e = 0) followed by e-BH (Vovk-Wang calibrators; Ren-Barber derandomisation [U]). The hypothesis being tested is then the delta-neighbourhood hypothesis.
- Stability index: pi_m = P(MP matched and existence rejected) over 50 perturbation draws. Each draw varies the thinning seed, engine seed, K ± 1, a tumour bootstrap, leave-one-cohort-out and the gene-universe variant.
- Stability is reported separately from validity.

(9) DECOY COHORTS (keep the nuisance, remove the co-expression)
- For each tumour, simulate X~ ~ NB(nu^{N4}_tij * xi_tij, b_marg,tj) on the same cells, covariates, copy numbers, clones, ambient profile and doublet partners.
- xi injects misspecification: the ambient profile comes from a disjoint empty-droplet subset, dosage gets a quadratic term, and per-cell dispersion is jittered.
- This approximately keeps each gene's marginal variance but removes co-expression beyond nuisance. ClusterDE/scDesign3 nulls instead keep gene-gene correlation, and with it the programs.
- The identical pipeline is run on the decoys. Empirical FDR = mean number of decoy MPs existing at N4 / number of real MPs existing at N4.

(10) PROTOCOL AND COHORT TYPING
- Protocol: logistic mixed model of tumour membership in S(m) on protocol, with a cancer-type random intercept; permutation p-value obtained by shuffling protocol labels within cancer type.
- Matched-sample invariance is tested where the same tumour was profiled with several protocols.
- Positional enrichment: circular genomic permutation of the MP's gene weights.
- Clone association: eta^2 of fold-C usage on clone, with a within-tumour permutation test.

### Inference and scalability

ALGORITHM (per thinning replicate; tumours run independently except steps 5-6)
1. Per tumour: estimate b~ (pass 1), thin into folds A/B/C, and fit the nuisance GLMs for rungs N0..N4 in each fold.
2. Per tumour: run the engine on fold A over a K grid (default 2..15) and choose K on fold B (held-out deviance, 1-SE rule). If K changed, re-estimate b~ conditional on K and re-thin (at most 3 passes).
3. Per tumour: apply the fold-B cross-fold (CF) filter to the programs.
4. Global gene universe: the union of each tumour's fold-A top 3,000 genes by deviance against N0, kept if detected in >= 10% of tumours, capped at about 8,000 genes.
5. Global: tumour-distinct HDBSCAN on about 50k-150k program nodes, giving candidate MPs and their hierarchy.
6. Per tumour x MP x rung: fold-C CF statistic plus its empirical null. Then, globally: partial conjunction, lower bounds, BH/e-BH, gene score tests with random-effects pooling, merge/split tests and typing.

Around these replicates:
- R = 5 thinning replicates, followed by derandomisation.
- 50 cheap perturbation draws that reuse the engine fits.
- B_null = 3 decoy cohorts.

COMPLEXITY
- Thinning: O(nnz).
- GLMs: O(n_t |G| p^2) per fold and rung, batched on GPU.
- All tests reduce to dense products, (R^A M) ∘ (R^C M) - (R^A ∘ R^C)(M ∘ M), costing O(n_t |G| |M|) per tumour.
- Empirical null: O(n_t |G| x 1,000). This is the dominant inference cost; it runs on GPU.
- Gene tests are score tests, so no per-gene refits are needed.
- Partial conjunction and BH: O(T log T) per MP.

PLANNING ESTIMATES AT 3CA SCALE (to be confirmed in pilot P2)
Assumes about 3,000 tumours, 1-2M malignant cells, |G| ≈ 8,000 genes and |M| ≈ 150 MPs.
- Inference layer: about 3-8 A100-hours per thinning replicate.
- Gavish engine (6 ranks x 10 restarts), 5 replicates plus 3 decoy cohorts: about 2 x 10^4 CPU-hours in R, or about 1-2 x 10^3 GPU-hours with a validated torch reimplementation.
- cNMF and consensus scHPF: similar order of magnitude.
- GBCD: run per cancer type on fold A.

MEMORY AND ORCHESTRATION
- Under 4 GB per tumour. The pooled dense matrix is never built.
- Data stream from zarr-backed AnnData.
- Parallelised over tumours with Ray/Dask locally, or Snakemake/Nextflow with containers on HPC.

SMALL STUDIES
- Works from T >= 10 tumours.
- If T < r_min, only within-tumour inference is reported.
- A single tumour gets program-level tests only.

DATA CONSTRAINTS
- UMI counts are required for thinning.
- Read-count data (e.g. Smart-seq2) are detected automatically and routed to audit mode: PCR duplicates of one molecule would be split across folds and create positive fold dependence.

AUDIT MODE (no thinning)
- Tests external gene sets (e.g. the published 3CA MPs) on an independent cohort.
- A within-data off-diagonal U-statistic is valid here because m does not depend on the data.
- It is fast and immediately useful, so it is the main route to adoption.

### Outputs for users

1. MP TABLE (one row per MP):
   - ID and parent/child position in the hierarchy
   - consensus loading vector
   - type and subtype
   - existence q-value at every rung N0..N4
   - recurrence lower bound r_L, overall, per cancer type and per cohort; number of cohorts
   - protocol odds ratio with q-value; positional-enrichment q-value; clone eta^2
   - stability pi_m
   - identifiability flag (number of anchor genes, fraction of near-pure cells)
   - which engine(s) found it

2. GENE TABLE (per MP x gene):
   - pooled log fold change per SD of MP activity (eta), its SE, and the q-value against the minimum effect
   - heterogeneity I^2
   - core / variable / anchor label
   - gene stability pi_jm
   - chromosomal position
   - a ranked top-N list with FDR flags, which drops into existing module-scoring and MSigDB-style workflows

3. TUMOUR MEMBERSHIP (per MP x tumour):
   - CF Z*, q-value and specific (partial) q-value
   - for each MP's member set, a simultaneous bound on the number of true members

4. CELLS:
   - fold-A usages (valid for inference)
   - full-data activity scores adjusted for nuisance (descriptive)

5. PRIVATE PROGRAMS: within-tumour programs that pass the tumour-level test but belong to no MP. These are candidates for clone-specific biology.

6. DIAGNOSTICS:
   - Validity-stability plane: x = -log10 q at rung N4, y = pi_m. Quadrants are real and stable, real but fragile, nuisance but stable (the artefacts that cross-rank filters favour), and noise.
   - Decoy calibration curves (nominal vs empirical FDR)
   - Karyogram of MP gene weights
   - Protocol forest plots
   - Rank-selection deviance curves
   - Dispersion-sign diagnostics (estimated b~ vs conditional b, with Theorem 1 covariance)

7. AUDIT REPORT for any existing catalogue (3CA, GeneNMF output, a lab's own MPs): which MPs exist beyond each nuisance rung in a new cohort, with recurrence bounds.

8. VERSIONED REFERENCE (Zenodo) PLUS PROJECTOR. For a query tumour it returns a per-MP existence p-value and per-cell scores under the same nuisance model. It also detects novelty: significant residual cross-fold co-expression not explained by reference MPs.

9. RUN REPORT (HTML/Quarto) recording every parameter, seed and software version.

### Failure modes fixed

STEP BY STEP, FOR THE GAVISH/3CA PIPELINE

- S0, malignant calling and contamination. Rung N1 models heterotypic doublets and ambient RNA using the same sample's TME cells and empty droplets. Contamination MPs fail existence at N1 and are typed "technical: contamination", with their loading cosine to TME profiles reported. Target: the macrophage, T-cell and endothelial MPs removed by hand in 3CA v2.

- S1, top 7,000 genes by mean expression. Replaced by a gene universe chosen on fold A by deviance against the nuisance model. Ribosomal/mitochondrial "low-quality" MPs are typed at N1/N2 rather than admitted and then removed by hand. Target: the 7 LQ MPs in inds_rm.

- S1, centring plus zeroing of negatives.
  - Inference runs on counts with an explicit nuisance model, not on centred log data.
  - Users can plug in uncentred engines (cNMF-KL, scHPF, GBCD).
  - Even with the centred Gavish engine, inference stays valid because the tests use count-level folds.

- S2, fixed rank grid 4-9. Per-tumour K is chosen by held-out NB deviance on an independent fold, with the 1-SE rule and monotone-ascent dispersion. For any engine, cross-fold replicability is an alternative.

- S2, best of 10 restarts and non-identifiability.
  - R independent thinnings crossed with engine seeds; evidence is derandomised (Cauchy combination or e-values), so results do not hinge on any one seed.
  - Each program gets an identifiability flag (anchor genes, near-pure cells).
  - Stability pi_m is reported.

- S2, no borrowing of strength across tumours. Partial conjunction accumulates weak per-tumour evidence, and gene effects are pooled by random effects across expressing tumours.

- S3, top-50 truncation. Replaced by variable-length, FDR-controlled gene sets with a minimum effect, plus effect sizes, SEs, heterogeneity and core/variable/anchor labels.

- S4, intra_min = 35. Replaced by a within-tumour cross-fold existence test that does not favour dominant axes. Stability is measured but kept separate from validity.

- S4, inter_min = 10. Replaced by a recurrence test against a structured nuisance-only tumour, not random gene sampling (the random-gene null that 10/50 guards against has P ≈ 1e-12).

- S4, intra_max = 10. No longer needed: one K per tumour and one representative program per tumour per cluster.

- S5, greedy clustering. Replaced by an order-free level-set estimand (a tumour-distinct HDBSCAN cluster tree) with formal merge/split tests on fold C. The number of MPs equals the number of existence hypotheses rejected at FDR q.

- S6, manual curation. Replaced by the nuisance ladder plus protocol and cohort tests. Types are reproducible, error-controlled, and validated against the 11 MPs removed by hand in 3CA v2.

- S6, CNA checked only afterwards. Cis-dosage (N3) and clone (N4) enter inference at discovery time. Positional enrichment and clone eta^2 are reported. Dosage-driven programs are labelled, not deleted.

- S7, scores without uncertainty. The projector returns per-tumour existence p-values and nuisance-adjusted scores, and the reference carries gene-level SEs.

CONCEPTUAL FAILURE MODES (landscape §1.4)
1. Recurrence is not sufficient: existence now requires "beyond nuisance" at a stated rung.
2. No uncertainty: q-values and confidence bounds at MP, tumour and gene level.
3. Centring blocks joint inference: count-level folds plus cross-tumour pooling.
4. Biased selection filters: stability is decoupled from validity (the validity-stability plane).
5. Results depend on the cohort: cohort-level partial conjunction, leave-one-cohort-out stability, random effects.
6. Validation only after the fact: typing is built into inference.

SELECTION BIAS (DOUBLE DIPPING) IS REMOVED THROUGHOUT. Discovery, selection and inference use statistically independent molecule folds, so the same noise cannot both create and confirm a program. Decoys run through the raw Gavish pipeline measure how often the standard procedure manufactures MPs from pure nuisance.

### Validation plan

All targets below are pre-registered hypotheses, not results.

TIER 0: THEORY CHECKS (simulation)
Setting: NB factor, mixture and continuum models with known b.
- Theorem 1 (fold covariance): relative error < 5% across (mu, b, b~).
- Type I of tumour and gene tests at alpha 0.05:
  - exact and conservative b~: target 0.03-0.06;
  - Poisson splitting: expected to be inflated (shown as a negative control).
- FDR across MPs at q = 0.1.
- Accuracy of choosing K.
- Power with marginal vs conditional b~: the two-pass target is >= 85% of oracle power.

TIER 1: TUMOURNULL SEMI-SYNTHETIC BENCHMARK (released)
Base data: decoys fitted to real tumours, namely CanSig's 12 datasets (185 patients) plus about 300 3CA UMI tumours. Each decoy carries the measured depth, QC, Numbat copy number and clones, empty-droplet ambient profile, and real TME doublet partners.
Planted programs:
- pan-shared (80% of tumours);
- subset-shared (10-30% of tumours across >= 2 cohorts);
- single-cohort;
- tumour-private;
- clone-coupled;
- protocol-coupled (a dissociation-like signature only in "fresh" cohorts);
- overlapping pairs (cosine 0.3-0.6) and parent/child nests.
Program parameters: fold change 1.3-3, active-cell fraction 5-50%, program size 20-200 genes.
Grid: T in {20, 100, 500}; cells per tumour in {100, 300, 1000}; median UMI in {1k, 3k, 10k}. Pure-null cohorts are included.
Engines, each run raw and CALIPER-wrapped: Gavish code (tiroshlab/3ca), GeneNMF, cNMF, consensus scHPF, GBCD, DRVI, scVI+Hotspot (ITHmapper-style).
Metrics:
- realised MP-level FDR vs nominal;
- recall (Hungarian-matched cosine >= 0.8);
- gene AUPRC and gene-level FDR;
- tumour-membership FDR and power;
- coverage of r_L (target >= 95%);
- typing macro-F1 (target >= 0.85);
- merge/split accuracy;
- rate of false MPs aligned to injected nuisance;
- leakage of patient-private programs into MPs;
- agreement between engines before vs after wrapping (hypothesis: wrapping makes MPs engine-invariant).
Ablations:
- Poisson splitting;
- marginal-only b~;
- within-fold U-statistic;
- no empirical null;
- ClusterDE-style null that preserves correlation;
- recall knockoff genes;
- no derandomisation;
- Gavish 35/10/10 filters applied on top of CALIPER folds.

TIER 2: END-TO-END DECOYS ON REAL COHORTS
- Run the raw Gavish, GeneNMF and scalop pipelines on N4 decoys of CanSig and 3CA, and count how many MPs they call from pure nuisance.
- CALIPER target: empirical FDR within q + 0.02.

TIER 3: LABELLED ARTEFACTS
- Type the 11 MPs removed by hand from 3CA v2 (plus the LQ, doublet and single-study MPs removed in Gavish 2023) against the retained MPs using the nuisance ladder.
- Target AUROC >= 0.85.
- Report every disagreement on retained MPs.

TIER 4: ORTHOGONAL EVIDENCE
(a) Cross-protocol invariance.
- Data: Janesick FFPE breast (Flex / Xenium / Visium), Oliveira CRC (Flex / Visium HD / Xenium), and 3CA studies with matched sc and sn.
- Metrics: AUROC of the "protocol-sensitive" type for predicting failure to replicate in fixed or in-situ data; agreement of projected usage on panel genes for MPs typed as shared biology.
(b) Genetic orthogonality.
- scDNA-matched clonealign datasets and allele-based CalicoST clones. Expected: eta^2 of usage on DNA-defined clone is high for clone-linked MPs and low for non-genetic ones.
- KP-Tracer lineage: phylogenetic autocorrelation (heritability) for each type.
(c) Causal.
- Tahoe-100M: MPs beyond nuisance should respond to mechanism-matched drugs (CDK4/6 inhibitors → cycling; MEK inhibitors → MAPK) relative to DMSO more often than technical MPs do.
- Kinker CCLE: recurrence in cultured cells (cell-intrinsic).
- Replogle Perturb-seq: knocking down an MP's anchor genes should shift that MP.
(d) Spatial.
- Xenium / Visium HD: autocorrelation of usage and pre-registered niche associations, controlling for CalicoST clone domains.

TIER 5: REPLICATION AND ROBUSTNESS
3CA → TCCA replication, excluding shared studies. Targets:
- >= 80% of MPs beyond nuisance replicate, vs a markedly lower rate for technical MPs;
- genes at q < 0.05 replicate at a higher rate than Gavish top-50 genes ranked 26-50;
- r_L is calibrated in the held-out cohort.
Leave-one-study-out analysis.
Perturbing cohort composition by 10-50% (the v1 → v2 analogue): MP survival and ARI, compared with Gavish.

TIER 6: CLINICAL (secondary)
- TCGA associations compared against expression-matched random signatures (Venet 2011).

### Flagship biology

"AN AUDITED, ENGINE-ROBUST ATLAS OF MALIGNANT METAPROGRAMS"

SCOPE AND DEFINITION
- Re-derive MPs from every UMI tumour in 3CA v2. Target: >= 1,000 tumours across >= 8 cancer types; pilot P1 sets the final numbers.
- Use three engines, each wrapped by CALIPER: the published Gavish NMF code, cNMF, and GBCD run per cancer type.
- An MP is "engine-robust" if it exists at N4 (q < 0.05) under >= 2 engines and has pi_m >= 0.7.
- Replicate in TCCA after removing studies that overlap with 3CA.

Four questions, each with a falsifiable expectation:

Q1. Which canonical MPs survive each rung?
- Applies to Gavish 2023's 41 MPs and the 3CA v2 set, each with a per-cancer-type recurrence lower bound.
- Positive control: the low-quality and contamination MPs removed by hand should fail at N1/N2.
- Open question: how many of the retained MPs turn out to be technical, QC-entangled, protocol-sensitive or dosage-explained.

Q2. How much recurrent ITH is genetic?
- Per cancer type, report the fraction of MPs, and of tumour-level program variance, explained at N3 (cis-dosage) vs N4 (clone).
- Hypothesis: cancers dominated by recurrent arm-level events (e.g. GBM chr7 gain / chr10 loss) carry positional, dosage-typed programs, while ITH in carcinomas is mostly non-genetic.
- This turns the post hoc figures of Gavish (24%) and Kinker (8%/39%) into a calibrated, prospective estimate.

Q3. Decomposing "stress", one of the most used and most suspected MPs.
- The cluster tree plus specific tests split stress into child MPs.
- Test: one child should persist in dissociation-free Flex and Xenium data and sit in hypoxic/necrotic niches (in vivo stress); another should vanish outside fresh dissociation (ex vivo artefact).
- Validation: the matched-protocol breast and CRC blocks, plus spatial data.
- A clean split would settle a long-running ambiguity. It is the most persuasive single demonstration of "real biology vs artefact".

Q4. MPs that fixed rank grids and greedy clustering merged or missed.
- Report new subset-shared MPs and splits of hybrid MPs that have r_L >= 5 tumours in >= 2 cohorts and replicate in TCCA.
- Check each for cell-intrinsic recurrence in Kinker CCLE and for drug response in Tahoe-100M.

OPTIONAL WET-LAB TEST (strengthens the paper but is not required)
- Pick one new MP that is beyond nuisance and cell-intrinsic.
- CRISPRi of 3-5 anchor genes in two cell lines, read out by a small Perturb-seq.
- This tests whether the anchor genes causally move the program.

DELIVERABLE
- "audited-3CA" reference v1.0 on Zenodo. Each MP carries its type, recurrence bounds, gene effects with SEs, and stability.
- MSigDB-compatible exports, flagged by q-value.

### Software architecture

LANGUAGE AND ECOSYSTEM
- Core in Python: numpy/scipy, with torch as an optional GPU backend.
- scverse-native: AnnData/MuData with scanpy conventions. Results are written to .uns/.obsm/.varm and returned as pandas tables.
- R users: caliperR, built on reticulate. It converts between Seurat and AnnData and provides Bioconductor-style accessors.
- R engines (Gavish NMF, GeneNMF, GBCD) run through rpy2 inside versioned containers.
- Distribution: pip and bioconda, Docker/Apptainer images, and Snakemake and Nextflow pipelines for atlas-scale runs.

MODULES
- caliper.pp: setup(), the UMI check, and nuisance covariates:
  - ambient RNA, imported from CellBender/SoupX/DecontX or estimated from raw barcodes;
  - doublets via scDblFinder or scrublet;
  - copy number and clones, imported from Numbat or computed with infercnvpy using leave-window-out.
- caliper.thin:
  - Dirichlet-multinomial multi-fold thinning;
  - dispersion estimators 'marginal-lcb', 'underfit' and 'two-pass';
  - replicate manager;
  - fold-covariance diagnostics based on Theorem 1.
- caliper.nuisance: batched NB IRLS for ladder rungs N0..N4. Users can add rungs.
- caliper.engines: an Engine protocol, fit(counts, meta) -> ProgramSet(loadings, usages, mean=None, scope='tumour'|'pooled'). Built-in engines:
  - GavishNMF: a Python/torch reimplementation regression-tested against the tiroshlab/3ca outputs, plus an R-exact mode;
  - CNMF (wraps cnmf);
  - SCHPF and consensus scHPF;
  - GeneNMF (rpy2);
  - GBCD (rpy2);
  - DRVI (scvi-tools);
  - ScviHotspot;
  - Custom.
- caliper.select: held-out deviance, cross-fold replicability, and delta calibration.
- caliper.aggregate: tumour-distinct HDBSCAN (fast_hdbscan or hdbscan with a custom core distance, kNN via pynndescent), the cluster tree, and consensus loadings.
- caliper.infer:
  - CF tests and empirical nulls;
  - partial conjunction, lower bounds and closed testing;
  - gene score tests and random-effects pooling;
  - BH/e-BH, Cauchy combination and derandomisation.
- caliper.typing: ladder logic, protocol mixed models with permutations, positional enrichment, and clone association.
- caliper.stability: perturbation scheduler and stability indices.
- caliper.decoy: decoy generator, misspecification injection, and the spike-in simulator (TumourNull).
- caliper.audit: tests external gene-set catalogues; ships with the 3CA, Gavish, Kinker and Barkley sets.
- caliper.ref and caliper.project: reference building and versioning, projection with p-values, and novelty detection.
- caliper.pl: validity-stability plane, karyogram, calibration curves, and forest plots.

KEY API
  cp.pp.setup(adata, tumour='sample', cohort='study', protocol='chemistry', tme=adata_tme, cna='numbat')
  res = cp.run(adata, engine=cp.engines.GavishNMF(ranks='auto'), folds=(0.5,0.2,0.3), n_thin=5, ladder=['N0','N1','N2','N3','N4'], q=0.05, r_min=5, s_min=2, decoys=3)
  res.mps; res.genes; res.membership; res.private; res.pl.quadrant()
  cp.audit(adata_new, gene_sets=cp.datasets.gavish2023(), q=0.05)
  ref = cp.ref.build(res, version='1.0'); cp.project(adata_query, ref)
  cp.bench.run(engines=[...], scenario='TumourNull-default')

ENGINEERING STANCE
- No new optimiser is written.
- The novel code is thinning, GLMs, matrix-product tests, multiple-testing procedures and orchestration. All of it can be unit-tested against closed-form nulls.
- At every release, continuous-integration tests check null calibration on small simulated cohorts.

### Feasibility

TEAM (about 3.2 FTE)
- 1 statistician (postdoc): theory (proof of Theorem 2; validity of partial conjunction and closed testing under thinning) and the inference code.
- 1 research software engineer / computational biologist: the package, engine adapters, containers and workflows.
- 1 PhD student or postdoc: the TumourNull benchmark, decoys and validations.
- PI at about 20%.
- Collaborators for matched-protocol and spatial data (public 10x datasets suffice at first) and an optional wet lab.

COMPUTE
- About 2-4 x 10^4 CPU-hours plus 2-4 x 10^3 A100-hours over the project, and 10 TB of storage.
- Firm numbers come from the pilot P2 timings.

DATA
- 3CA v2: availability of raw UMI counts to be audited in pilot P1.
- CanSig's 12 datasets: public.
- TCCA (bioRxiv 2026): access still to be confirmed.
- CELLxGENE tumour sets.
- Janesick and Oliveira matched-protocol blocks: public 10x data.
- Kinker CCLE, Tahoe-100M, Replogle Perturb-seq, KP-Tracer: public.
- clonealign and CalicoST datasets: public; sequencing depth still to be checked.

TIMELINE (about 17 months to submission)

| Months | Work |
|---|---|
| 0-2.5 | Pilots P1-P7 (go/no-go). |
| 3-7 | Package v0.1: thinning, ladder GLMs, CF and gene tests, partial conjunction, Gavish and cNMF adapters, decoys, audit mode. Write up Theorem 1 and Corollary 1. |
| 6-10 | TumourNull benchmark with 7 engines plus ablations; proof of Theorem 2; GBCD, DRVI and Hotspot adapters. |
| 8-9 | Preprint of the method plus an audit of the published 3CA lists on independent cohorts (insurance against being scooped). |
| 9-14 | Flagship 3CA re-derivation, TCCA replication, orthogonal validations (protocol, clone, Tahoe, spatial). |
| 14-17 | Reference v1.0, manuscript, reproducibility repository. |

The optional wet-lab experiment runs in parallel from month 10 and takes 4-6 months.

LOW ENGINEERING RISK
- Every component is standard numerical linear algebra or GLM fitting.
- Existing engines are reused, not rewritten.
- Audit mode gives users something useful from month 4.

### Risks and mitigations

1. POWER LOSS FROM THINNING
Fold C holds only 30% of molecules, which hurts small or shallow tumours.
- Off-diagonal statistics pool many gene pairs, and partial conjunction pools tumours.
- The two-pass dispersion estimate recovers most of the oracle power. In the toy run, marginal b kept about two-thirds of oracle power at fold change 1.5.
- Use an adaptive eps_C for small tumours, and cross-fitting with rotated fold roles plus derandomised combination.
- Pilot P4 sets the minimum-cell guidance.

2. NUISANCE-MODEL MISSPECIFICATION
Residual ambient RNA or nonlinear dosage can miscalibrate the tests.
- Negative-control empirical nulls, matched on expression, ambient weight and chromosome arm.
- Decoys are deliberately simulated from a richer model than the one used for inference.
- Sensitivity analysis across ladder variants.
- Every run reports calibration curves.

3. QC COVARIATES ENTANGLED WITH BIOLOGY
For example, complexity tracks cell cycle, and %mito tracks hypoxia/OXPHOS.
- N1 (exogenous) and N2 (QC) are separate rungs.
- An MP that fails only at N2 is labelled 'QC-entangled', a flag rather than a deletion.
- The flag is resolved with protocol, cell-line and perturbation evidence.

4. CIRCULARITY OF EXPRESSION-INFERRED COPY NUMBER
Copy number inferred from expression can absorb programs whose genes cluster in the genome (keratin, HOX, MHC-II, histone clusters).
- Use allele-aware Numbat where BAM files exist; otherwise compute dosage with a 2 Mb window around each gene left out.
- Validate with DNA- or allele-based clones (clonealign, CalicoST).
- Cis-dosage (N3) and clone/trans (N4) effects are tested separately.

5. GAPS IN THE THINNING THEORY
- Theorem 1 and Corollary 1 are exact: proved here and checked numerically.
- Theorem 2 is so far a first-order argument. If a general proof fails, state it under a monotone-engine assumption, check it per engine on decoys, and keep conservative dispersion as the default.

6. NON-UMI DATA (e.g. Smart-seq2) CANNOT BE THINNED
PCR duplicates break independence between folds.
- Such data are detected automatically and routed to audit mode.
- The flagship analysis is restricted to UMI datasets.

7. RAW COUNTS MAY BE UNAVAILABLE FOR MUCH OF 3CA
- Pilot P1 decides.
- Fallback flagship: CanSig + TCCA + CELLxGENE tumours, plus an audit of the published 3CA MPs on independent cohorts.

8. REVIEWER OBJECTION: "just a pipeline of known tests"
Counter with:
- new theory (fold covariance, sign control, never pairing the same matrix entry across folds);
- a formal level-set estimand;
- the decoy demonstration that standard MPs can arise from pure nuisance;
- engines becoming interchangeable once wrapped;
- the flagship stress decomposition.

9. TYPING DISPUTES (e.g. "stress is real biology too")
- Types are hypothesis outcomes with q-values, never deletions.
- Users see the result at every rung.

10. SCOOP RISK
The GBCD, Witten-lab or Tirosh-lab groups could add inference to their methods.
- Release an early preprint (months 8-9) and audit mode early.
- Run the manual scoop check (landscape §6.3) before committing.
- My limited web searches this session found no MP method built on thinning or partial conjunction. This is weak evidence, because arXiv and bioRxiv were partly blocked.

11. USERS WANT SHORT GENE LISTS
- Ship ranked top-N lists with q-value flags and an MSigDB-compatible export, so existing workflows keep working.

### Go/no-go pilots

P1. DATA AUDIT (2 weeks, 1 person)
- GO if: at least 1,000 tumours across at least 8 cancer types have raw UMI counts (pooling 3CA v2, CanSig, TCCA and CELLxGENE), at least 70% have TME cells, and at least 30% have empty droplets or BAMs.
- Otherwise: shrink the flagship to the available sets and estimate ambient RNA from low-UMI barcodes.

P2. DECOY DEMONSTRATION (3 weeks)
- Build N4 decoys for CanSig's 12 datasets and run the unmodified Gavish code, GeneNMF and scalop on them.
- GO if: the standard pipelines call at least 3 MPs, or at least 15% of their real-data MP count, from pure nuisance in at least 2 cancer types.
- If they call about 0: the false-MP headline weakens, so re-centre the paper on gene-level FDR, recurrence bounds and typing.

P3. THINNING CALIBRATION (2 weeks, simulation only; partly done in this session)
- GO if, under mixture or continuum structure with conservative b~:
  - type I of the CF and gene tests is at most alpha + 0.01;
  - the inflation from Poisson splitting is reproduced;
  - two-pass power is at least 70% of oracle.
- If power falls below 50%: make conditional dispersion the default, with recalibration on decoys.

P4. PER-TUMOUR POWER (2 weeks)
- On 200 real tumours, measure the fraction of tumours in which strong canonical programs (G2/M, G1/S, hypoxia) are detected with eps_C = 0.3.
- GO if: the median is at least 0.6 for tumours with at least 200 cells.
- Otherwise: raise eps_C, use cross-fitting, or impose minimum-cell rules.

P5. ARTEFACT TYPING (2 weeks)
- Run the ladder tests on 3CA v2 programs and compute the AUROC for separating the 11 hand-removed MPs from retained ones.
- AUROC at least 0.85: GO.
- 0.70 to 0.85: add contamination and ambient features.
- Below 0.70: redesign rungs N1/N2.

P6. CNA AUDIT (2 weeks)
- Measure the fraction of 3CA MPs that fail at N3/N4, plus positional enrichment.
- If below 5%: keep dosage typing as an annotation rather than a headline (it still serves Q2).

P7. SCOOP CHECK (1 week)
- Search bioRxiv and Google Scholar manually, following landscape §6.3, plus these terms: 'count splitting NMF', 'data thinning gene program', 'partial conjunction metaprogram', 'replicability gene programs single-cell'.
- Also check papers citing GBCD, TCAT and Gavish.
- NO-GO only if a preprint already combines thinning-based inference with a test of recurrence across tumours.

Total: about 10 weeks with 2 people working in parallel. The outcomes of P2 and P5 decide how the headline is framed.

### Why Nature Methods

1. A conceptual advance stated as a definition.
- "Recurrent" becomes "replicable beyond a stated nuisance class at a stated error rate". This works because the MP is treated as an estimand, separate from whatever estimator produces it.
- The reframing is general. Any multi-sample program discovery inherits it: TCAT-style immune catalogues, organoid screens, developmental atlases, spatial niches.

2. Theory the field needs and that reviewers can check.
- The exact fold covariance under misspecified NB thinning.
- Sign control that holds for any latent structure, which answers a published objection to thinning (Hivert et al.).
- The design rule that no test pairs the same matrix entry across folds.
- Validity of partial-conjunction recurrence bounds on thinned folds.

3. Broad adoption with low friction.
- The layer is engine-agnostic: users keep Gavish NMF, GeneNMF, cNMF or GBCD and add one call.
- Audit mode works today on existing gene-set catalogues.
- This fits the Nature Methods pattern of a widely reusable tool, rather than yet another factorisation in a crowded field.

4. A benchmark that measures failure modes, not just planted signal.
- TumourNull includes measured confounders, pure-nuisance cohorts and calibration metrics.
- This meets the bar set by scIB and BEELINE-type contributions.

5. A flagship audit of a widely used resource (the 3CA MPs and the MSigDB GAVISH sets).
- Validated orthogonally: cross-protocol, DNA-defined clones, perturbation and spatial data.
- Concrete biological payoff: decomposing stress into in-vivo and ex-vivo components, and estimating how much recurrent ITH is genetic.

6. Clear separation from the closest work.
- GBCD, ITHmapper, TCAT and Bayesian multi-study NMF are estimators or catalogues.
- None of them offers calibrated existence, recurrence bounds, tumour- and gene-level FDR, or nuisance typing.
- CALIPER can wrap and audit them all.

REFERENCE STATUS
- [V] items come from the verified landscape.
- The [U] references are standard statistics recalled from memory and not verified this session; check them before writing:
  - Benjamini & Heller 2008 (partial conjunction)
  - Goeman & Solari 2011 (closed testing)
  - Wang & Ramdas 2022 (e-BH)
  - Vovk & Wang 2021 (calibrators)
  - Ren & Barber 2024 (derandomised knockoffs)
  - Efron (empirical null)
  - DerSimonian-Laird
  - McCarthy & Smyth (TREAT)
  - Liu & Xie (Cauchy combination)
  - Campello et al. (HDBSCAN)
- arXiv 2409.03069 (a discussion of misspecified NB fission) was seen only as a search snippet.
- Theorem 1 and Corollary 1 were derived in this session and checked by Monte Carlo. The proofs still need independent review and a comparison against the NB count-splitting paper, which may contain an equivalent expression.

### Differentiation from closest competitors

| Competitor | Differentiator |
|---|---|
| GBCD (Liu, Carbonetto, Willwerscheid, Oakes, Macleod, Stephens; Nature Genetics 2025) | GBCD is an estimator: Gaussian empirical-Bayes matrix factorisation (EBMF) on pooled log data. It gives gene-level lfsr within a program. It does not test whether a program exists, or in how many tumours or cohorts; it does not control error on tumour membership; and it does not model CNA, technical or protocol effects. CALIPER runs GBCD as a plug-in engine on fold A and adds: count-level validity despite GBCD's Gaussian model, partial-conjunction recurrence with lower bounds, FDR-controlled tumour and gene membership, nuisance-ladder typing, and decoy calibration. The pitch is 'CALIPER makes GBCD's (and NMF's) programs testable', not 'better shared-program finding than GBCD'. |
| ITHmapper / Pan-Cancer ITH compendium (Nofech-Mozes ... Abelson; bioRxiv 2026) | ITHmapper defines discrete states from scVI with sample as the batch key, plus Hotspot modules. It relies on multi-seed heuristics, uses CopyKAT only for malignant calling, and has no error control. CALIPER can wrap scVI+Hotspot modules as an engine and test whether they exist within tumours beyond nuisance, without the erasure that sample-as-batch causes. It also returns overlapping continuous programs with gene-level FDR, types states as clone-linked or protocol-sensitive, and can audit ITHmapper's own states. |
| TCAT/starCAT (Kotliar ... Raychaudhuri; Nature Methods 2025) | TCAT is a T-cell program catalogue built with Harmony + cNMF; reproducibility across datasets is judged by heuristics and expert judgement, and projection uses NNLS without uncertainty. CALIPER contributes what such catalogues lack: a formal cross-cohort replicability test (cohort-level partial conjunction), nuisance typing, and a projector with per-tumour p-values. Its reference is a deliverable, not the novelty. It also handles malignant cells with CNA, and could audit TCAT programs. |
| Bayesian multi-study NMF (Grabski, Trippa, Parmigiani; Genome Biology 2025), plus Tetris / APAFA / SUFA | These are model-based posteriors over which studies share which factors. Their calibration depends on the model being correct, they were built for mutational signatures or Gaussian bulk data, and they scale poorly. CALIPER is frequentist and model-agnostic. Its validity under misspecification comes from molecule-level thinning, off-diagonal statistics and empirical nulls. A single-cell port of any of these models becomes one more engine, so the relationship is complementary rather than competing. |
| scDEF, consensus scHPF, cNMF, CSI-GEP, staNMF | These are program-discovery engines or stability heuristics. CSI-GEP replicates across splits of cells within one dataset. CALIPER splits molecules, so every cell is kept, which matters for small tumours. It replicates across tumours and cohorts, and pairs stability with null calibration, because stable technical programs pass stability-only criteria. |
| Count splitting / data thinning (Neufeld, Witten et al.), ClusterDE, recall knockoffs | These are generic tools for differential expression after clustering. CALIPER extends thinning to factor programs and to recurrence across tumours, with new theory: (1) the exact fold covariance under misspecified dispersion; (2) sign control for any latent structure; (3) a rule never to pair the same matrix entry across folds. Together these address the Hivert et al. limitation. In-session toy runs show Poisson splitting gives about 99% false positives with a selection-heavy engine. CALIPER's decoys remove co-expression, whereas ClusterDE nulls preserve gene-gene correlation and so preserve programs. |
| CanSig benchmark (Boeva lab; Cancer Research 2026) | CanSig scores rediscovery of known signatures; it has no planted confounders, no null cohorts and no calibration metrics. CALIPER's TumourNull benchmark plants measured confounders and pure-nuisance cohorts, and scores FDR calibration, coverage of the recurrence bound, and typing accuracy. CanSig's 12 datasets are reused as a harness. |
| Gavish/3CA code, GeneNMF, scalop/pyscalop, mpnmf | These become engines. CALIPER replaces their thresholds, greedy or forced clustering, and manual curation without asking users to give up their preferred factorisation. Adding two lines to an existing workflow is the adoption path. |
| KINOMO; sciRED | KINOMO applies kernel correction for genomic variability inside NMF. sciRED regresses known covariates with a Poisson GLM and yields signed factors. CALIPER treats cis-dosage and clone (trans) effects as separate nested hypotheses with error control, and keeps dosage-driven programs as labelled biology rather than removing them. |

### Referee 1: novelty and prior art

**Predicted editorial outcome:** Borderline: desk reject or send to review, and most likely rejected after review in its current form. Auditing the 3CA / MSigDB GAVISH metaprograms is a strong editorial hook, and an editor could send it out on that alone. The paper rests on two things. One is the decoy result (pilot P2), that the standard pipelines call metaprograms from pure nuisance. The other is the stress split. If P2 comes back weak (the standard pipelines call few or no MPs from nuisance), the headline collapses and I expect a desk reject or a transfer to Nature Communications or Genome Biology. If it goes to review, I expect referees to (a) call it a well-built pipeline of known tests: count splitting (Neufeld/Witten), partial conjunction (Benjamini-Heller), BH/TREAT, DerSimonian-Laird random effects, HDBSCAN, ClusterDE-style nulls, module-preservation-style audit; (b) dispute the theory, since Theorem 1 is probably already in the NB count-splitting literature and Corollary 1 is false in exactly the setting Hivert et al. study; and (c) find the central estimand not specific: 'tumour expresses m' (kappa > 0) is true for almost any tumour with co-expression that overlaps m's genes. Best realistic outcome: major revision that requires the tumour-level estimand to be redefined, typing to use equivalence tests, and decoy calibration to be checked against real negative controls. The paper that clears the bar is narrower: calibrated recurrence and audit inference for published MP catalogues, plus a striking false-MP demonstration on real data.

**Fatal flaws**

- The tumour-level estimand is not specific. S^N(m) = {t : kappa_t(m) > 0}, with m non-negative and a point null. Any tumour carrying ANY residual program w that shares genes with m has kappa_t(m) ≈ lambda[(m·w)^2 − Σ m_j^2 w_j^2] > 0. With n_t in the hundreds to thousands, the CF test then rejects. I checked this in my own simulation (1,000 cells, NB size 3, exact thinning, the CF statistic from step 5; script at [referee_checks/ref1_check.py](referee_checks/ref1_check.py)). The tumour carried an unrelated program w (fold change 2.5 in 30% of cells) and never expressed m. Sharing 5 of 40 genes with w, m was declared 'expressed' in 72% of tumours; sharing 10 of 40, in 100%. At fold change 1.5 the 10/40 rate was 55%. The consensus location (median over all genes with positive median loading) will be broad, which makes this worse. The 'specific' version conditions only on cluster-tree siblings, not on every program in the tumour. So r_L, tumour membership and typing all measure overlap with whatever co-expression the tumour has, not whether the tumour expresses the MP. The proposal itself adds a minimum-effect null for genes (eta_0) because atlas-scale n makes everything significant, but does not apply the same logic to kappa.
- Typing does not have the error control it claims. A type is set by 'the first rung at which existence FAILS', i.e. by NOT rejecting a null. Failing to reject is absence of evidence. Richer rungs cost power because they add parameters and their covariates are collinear with biology (%mito with hypoxia/OXPHOS, complexity with cell cycle, expression-inferred copy number with genomically clustered programs, clone with spatial niche). Typing will therefore systematically push real but weaker programs into technical, QC, CNA or clone types. 'Error control at each rung of the type ladder' controls only the rejections, not the type assignments. As written, the claim that typing 'replaces manual curation' with 'error-controlled' types is false. It needs equivalence or interval tests (an upper confidence bound on kappa^{N+1} relative to kappa^N) plus an explicit 'indeterminate' class.
- The headline theory does not survive. (i) Corollary 1 ('for ANY latent structure, b_marg <= b, so every cross-fold covariance <= 0') assumes one NB size shared by all cells. It fails with state-specific dispersion, which is precisely Hivert et al.'s setting (their v3 explicitly attributes inflated type I under marginal thinning to component-specific overdispersion [V]). Counterexample, checked by Monte Carlo: a 50/50 mixture of (mu=100, b=100) and (mu=1, b=0.5) gives moment-matched b_marg ≈ 1.02 > 0.5. Within the second component the fold covariance is +0.125 (Monte Carlo) against +0.129 (the proposal's own Theorem 1), which is anti-conservative. The claim that it 'answers Hivert' is therefore wrong. (ii) Theorem 1 is a two-line moment calculation. NB count splitting (Neufeld, Popp, Gao, Battle, Witten, arXiv 2307.12985 [V exists]) studies misspecified overdispersion and fold dependence and very probably contains an equivalent expression or sign statement. I could not check the formula because arXiv is blocked. The proposal admits this possibility. (iii) Theorem 2, the result that would make the layer valid for 'any engine', is an unproven first-order sketch. The Efron-style empirical-null standardisation also replaces H0 with 'no more co-expression than expression-matched pseudo-programs'. In practice, then, validity rests on a heuristic empirical null, not on the theorems.
- The central principle contradicts itself. 'The MP is the estimand; any factorisation is an estimator' is stated, but the location is then defined as engine-indexed (W^E_t, the modes of F_delta for engine E). The estimand therefore changes with the estimator. All the layer can certify is the existence of whatever directions each engine proposes, and 'engine-robust' becomes a heuristic cosine intersection across engines. 'Wrapping makes MPs engine-invariant' is offered as a pre-registered hypothesis, but nothing in the method makes it true. Without an engine-free definition of the location, the conceptual-advance claim made in the Nature Methods pitch does not hold.

**Major concerns**

- Novelty is mostly in the integration. An engine-agnostic inference layer built on molecule-level thinning is the stated framing of count splitting (Neufeld et al., Biostatistics 2024: valid inference after latent-variable estimation 'for virtually any latent variable estimation technique' [V]). Thinning-based choice of the number of factors for count factor models already exists with consistency theory (Wang, Xu, Zhao, Wang, JASA 2025, plus CRAN package tcv [V]). A partial-conjunction lower bound on the number of non-null studies is textbook replicability analysis. Recent PC work already handles dependent studies with e-values (Dey et al., arXiv 2511.04130 [V]). Audit mode is module-preservation testing (Langfelder et al., PLoS CB 2011; NetRep, Cell Systems 2016 [V]) on GLM residuals. Excess co-expression beyond depth and measurement error, tested on counts, is CS-CORE (Nature Communications 2023 [V]) and its 2026 framework with simulation-based empirical p-values (Shan, Lin, Zhao, Briefings in Bioinformatics 2026 [V]). Density-based clustering of programs goes back to Kinker 2020, which used DBSCAN. What is left as new: applying these tools to MP recurrence across tumours, the specific cross-fold off-diagonal statistic, and the nuisance ladder. Referee 2 will call this incremental.
- Power. Fold C holds 30% of molecules. The author's own toy runs reach only 0.26–0.39 power at fold change 1.5 with 300 cells. Typical 3CA tumours have 100–500 malignant cells at 1–3k UMI. Simes partial conjunction needs the r-th smallest tumour p-value to be small, so r_L will sit far below true recurrence for anything except strong canonical programs (cycle, hypoxia, IFN). The 'audited atlas' could end up as a short list of already undisputed MPs. Under-powered rungs then feed fatal flaw 2.
- The decoy FDR is circular. Decoys are drawn from the same NB nuisance family that the inference fits, perturbed by misspecifications the author chooses. Real co-expression nuisances outside N0–N4 are absent: cell size and mRNA content, dissociation and immediate-early response (classed as biology), transcriptional bursting, global cell-cycle coupling, spatially structured ambient RNA, and malignant–malignant doublets. Empirical FDR on these decoys measures internal consistency, not real-world FDR. The method is also scored on a benchmark (TumourNull) generated by its own decoy simulator.
- Scope. The estimand is explicitly within-tumour heterogeneity only; between-tumour baselines go to an 'optional pseudobulk extension'. That keeps Gavish's main blind spot, the one GBCD (Nature Genetics 2025 [V]) and the Tirosh lab's multilayered GBM paper (Nature Genetics 57:1155–1167, 2025; now [V]) address. 'Replaces the 3CA pipeline' is therefore overstated.
- Data feasibility. N1 needs empty droplets, which public deposits rarely include. Allele-aware Numbat needs BAMs, which are often access-controlled. Much of 3CA is Smart-seq2 or processed matrices. Where data are missing, N3 falls back to expression-inferred copy number with a 2 Mb leave-out. Expression callers smooth over arm-scale windows, so that remains circular for broad events. TCCA access is unconfirmed. A 30%-with-BAMs-or-empty-droplets GO threshold looks optimistic; pilot P1 may well shrink the flagship to CanSig plus CELLxGENE.
- The ladder's order sets the answer. Nested rungs attribute shared variance to whichever nuisance comes first, and the order (ambient/doublet, then QC, then dosage, then clone) is arbitrary. Clone indicators from expression-based callers partly encode expression programs. A hypoxia program in a spatially segregated clone will be typed 'clone-linked'. TreeAlign (allele-specific dosage effects of subclonal CNAs with scDNA+scRNA [V exists; venue U]) already models clone dosage effects more rigorously when DNA is available.
- Statistical details. (a) The library size l_ti is computed on full data and used as the offset in every fold, which leaks fold-A molecules into fold-C residuals; it should be computed per fold. (b) The design rule 'never pair the same entry across folds' is broken through m~, which depends on X^A_ij' and multiplies r^C_ij'. That route is exactly the one Theorem 2 has to control. (c) Within-cohort tumours share ambient profile, chemistry and processing, so cohort-level Simes needs PRDS and tumour p-values are not independent within cohort. (d) Matching MPs across thinning replicates by cosine and then combining with Cauchy or e-values tests a 'delta-neighbourhood' hypothesis that is never defined formally. (e) BH within each tumour on fold-B programs, then BH across MPs, then BH across genes, gives no family-wise statement for the reported atlas.
- Threshold proliferation in statistical form. eps split, rho = 3, delta from the 5th percentile, r_min = 5, s_min = 2, eta_0 = log 1.2, q_T = 0.1, gene q = 0.05, top-3,000 per tumour, detection in >= 10% of tumours, cap of 8,000 genes, pi_m >= 0.7, K grid 2–15, 1-SE rule. A sceptic will see the 35/10/10 problem moved to new places. A sensitivity analysis over these choices is mandatory.
- The validation targets are weak or already known. The 11 hand-removed 3CA v2 MPs give 11 positives, labelled by one lab, 7 of them trivially ribosomal or mitochondrial. AUROC >= 0.85 has a very wide CI and has to beat trivial classifiers (ribo/mito fraction, cosine to TME profiles). The stress split: dissociation signatures are known (van den Brink 2017 [U]); 3CA v2 already compares sc and sn; the 2026 catalogue already carries 'Stress (in vitro)' [V+ per landscape]; and the matched Flex/Xenium blocks cover 2 cancer types and few patients. Q3 is likely to be under-powered or to confirm what is already known.
- Overscoped for 3.2 FTE in 17 months. The plan covers new theory, 7 engine adapters, a benchmark, a pan-cancer re-derivation, TCCA replication, Tahoe, Perturb-seq, spatial, clone and lineage validations, a reference with a projector, and an optional wet lab. The compute multiplies: 5 rungs × 5 thinnings × 3 decoy cohorts × 50 perturbations × 1,000 pseudo-programs × 3 engines. Users will struggle to adopt something this heavy, which undercuts the 'add two lines' adoption argument.

**Prior-art threats**

- Neufeld, Gao, Popp, Battle, Witten. 'Inference after latent variable estimation for single-cell RNA sequencing data', Biostatistics 2024 (count splitting). Explicitly claims valid inference after virtually any latent-variable estimator, i.e. the 'engine-agnostic thinning layer' framing. [V: search snippet and PMC listing]
- Neufeld, Popp, Gao, Battle, Witten. 'Negative binomial count splitting for single-cell RNA sequencing data', arXiv 2307.12985. Dirichlet-multinomial multi-fold NB splitting; discusses how misspecified overdispersion creates fold dependence. Probably anticipates Theorem 1 or an equivalent sign result. [V that the paper exists; exact formula U because arXiv is blocked]
- Hivert, Agniel, Thiébaut, Hejblum. arXiv 2405.13591 v3 (Aug 2025). Marginal NB thinning inflates type I because overdispersion is component-specific. Corollary 1 does not hold in this setting (counterexample in the fatal flaws). [V]
- Wang, Xu, Zhao, Wang. 'Data thinning for Poisson factor models and its applications', JASA 121(553), 2025/26, plus CRAN package tcv. Thinning-based cross-validation and information criterion for the number of factors, with selection consistency. Undercuts the novelty of held-out-deviance rank selection by thinning. [V]
- Langfelder et al. 'Is my network module preserved and reproducible?', PLoS Comput Biol 2011; Ritchie et al., NetRep, Cell Systems 2016. Testing whether co-expression modules replicate in independent data against permutation nulls. Anticipates CALIPER's audit mode. [V]
- CS-CORE (Su et al., Nature Communications 14:4846, 2023) and Shan, Lin, Zhao, Briefings in Bioinformatics 27(3), 2026. Count-level co-expression beyond depth and measurement error, with valid tests and simulation-based empirical p-values. Anticipates Omega^N/kappa as 'excess co-expression under a measurement model'. [V]
- Benjamini & Heller 2008, partial conjunction [U as cited in the concept]; Bogomolov–Heller replicability line [U]; Dey, Banerjee, Bhuyan, Majumdar, arXiv 2511.04130 (e-Filter for partial conjunction under dependent studies) [V]. The recurrence lower bound is standard replicability methodology.
- ClusterDE / synthetic-control null (Song, Li et al.; PMC10401959) with scDesign3 (Nature Biotechnology 2024). Synthetic nulls passed through the identical pipeline. [V]
- Kinker et al. Nature Genetics 2020. Already used DBSCAN on NMF programs, so density-based program clustering is prior work. [V per landscape]
- TCAT/starCAT (Nature Methods 2025). Cross-dataset reproducible program catalogue plus projector; sets the Nature Methods bar for 'catalogue + reproducibility'. [V]
- GBCD (Liu…Stephens, Nature Genetics 2025). Shared vs patient-specific programs, including between-tumour programs that CALIPER puts out of scope. [V]
- Tirosh lab, 'The multilayered transcriptional architecture of glioblastoma ecosystems', Nature Genetics 57:1155–1167 (2025). Baseline programs plus genetic association of the heterogeneity layers. [V, this session via search]
- TreeAlign (Funnell et al.; 'Allele-specific transcriptional effects of subclonal copy number alterations enable genotype-phenotype mapping in cancer cells'). Clone-specific dosage modelling; identifies genomically independent expression programs. Overlaps the N3/N4 typing. [V exists; venue U]
- ITHmapper (Nofech-Mozes…Abelson, bioRxiv 2026.01.06.693992). Over 1,000 tumours, 14 cancer types, mapping tool; holds the positioning of 'the new pan-cancer ITH package'. [V]
- sciRED (Nature Communications 2025), a GLM nuisance model plus factor–covariate matching; KINOMO, kernel correction for genomic variability. [V per landscape]
- CSI-GEP (Cell Genomics 2025) and consensus scHPF. Replication across splits and consensus K. [V per landscape]
- Gene Program Dictionary (Xu…Ma, bioRxiv 10.1101/2025.10.29.685322). 1,534 consensus programs from 466 datasets via partial-correlation networks; a resource-level competitor for 'replicable programs at atlas scale'. [V]
- CanSig benchmark (Cancer Research 2026). An existing cancer-state discovery benchmark harness; reduces how novel a benchmark contribution can be. [V]
- 'Quantifying stability via count splitting to guide model selection in RNA velocity analyses' (Bioinformatics Advances 2026). Count splitting used for stability and model selection; the generic idea is spreading. [V snippet]

**Required changes**

- Redefine 'tumour t expresses MP m' so that it is specific. Condition kappa on ALL other MPs and on the tumour's private programs, not just cluster-tree siblings. Add a minimum-effect or relative-specificity margin (e.g. kappa_t(m | rest) >= kappa_0, or the share of m's co-expression not explained by other directions). Show in simulation that an unrelated program sharing 10–25% of m's genes does not make m 'expressed' (my check: 72–100% false 'expression' under the current definition).
- Replace non-rejection typing with equivalence or interval inference: an upper confidence bound on how much signal is left after adding rung N+1. Add an 'indeterminate / under-powered' type. Show typing accuracy as power falls with fewer cells and lower UMI.
- Fix the theory. Restrict Corollary 1 to a dispersion shared across cells, or bound component-specific heterogeneity. Reproduce Hivert's mixture-with-component-specific-dispersion design with the one-sided CF statistic and report type I. Prove Theorem 2 under stated engine assumptions, or drop the 'any engine' validity claim. Check Theorem 1 against Neufeld et al. 2307.12985 before claiming it.
- State plainly what the empirical-null standardisation tests. Report how often it absorbs real global structure. Show results with and without it.
- Either define an engine-free location estimand (e.g. a well-posed decomposition of the pooled residual covariance, with identifiability conditions) or withdraw the 'MP as estimand' and 'engine-invariance' claims.
- Validate decoy and FDR calibration against real negative controls, not only self-generated decoys: technical-replicate libraries, species or cell-line mixtures, the same sample run fresh vs Flex, sc vs sn of the same tumour, CCLE lines without TME. Include nuisance sources outside the ladder.
- Add simple baselines and ablate the whole architecture: (a) Gavish MPs + NetRep/module preservation on independent cohorts; (b) NB count splitting + per-tumour NMF + standard BH tests; (c) GBCD with lfsr. Show what each CALIPER component adds over the simplest thinning-based pipeline.
- Report real-data power: r_L against Gavish recurrence counts, the fraction of tumours detected for G2/M, G1/S, hypoxia and IFN, and minimum cells/UMI. If most MPs end up with r_L << the Gavish count, say so and reframe the paper.
- The artefact-typing AUROC must beat trivial classifiers (ribo/mito fraction, cosine to TME and ambient profiles, positional-enrichment score alone) and come with confidence intervals. With 11 positives, add independent labelled sets.
- Either handle between-tumour (baseline) programs or explicitly concede that scope to GBCD and the Tirosh multilayered decomposition, and soften 'replaces the 3CA pipeline'.
- Confirm data before committing: the number of 3CA/TCCA/CELLxGENE tumours with raw UMI counts, with empty droplets, and with BAMs for allele-aware Numbat. Restrict N1/N3 claims to that subset.
- Run a full sensitivity analysis over every tuning constant (eps split, rho, delta percentile, r_min, s_min, eta_0, gene universe, K grid).
- Cut scope to fit 3.2 FTE: inference layer + audit of published catalogues + one decisive real-data demonstration (false MPs from nuisance, or a well-powered protocol-invariance test). Make the clinical tier and wet lab optional, and make the reference a by-product.

**Ideas worth keeping**

- Separating discovery, selection and inference by molecule-level folds for MP recurrence, using one-sided cross-fold off-diagonal (j ≠ j') statistics. This removes the diagonal-entry dependence that inflates within-fold co-expression tests (the author's own toy run: within-fold U-statistic type I 0.087–0.42 vs about 0.03–0.04 cross-fold).
- A partial-conjunction lower confidence bound on the number of tumours and cohorts expressing an MP. It is interpretable and replicability-native, and it is the natural replacement for inter_min = 10 and the 'at least 5 programs' rule, provided the tumour-level test is made specific.
- Audit mode: testing published catalogues (3CA, MSigDB GAVISH, Kinker, Barkley) on independent cohorts against an explicit nuisance model. It is low-risk and immediately useful, and probably the publishable core.
- The validity–stability plane: reporting null-calibrated validity separately from stability, since stable technical programs pass stability-only criteria.
- Labelling dosage- and clone-driven programs as biology with a stated genetic basis, rather than deleting them.
- A decoy demonstration that the unmodified Gavish/GeneNMF/scalop pipelines call MPs from nuisance-only tumours built on real cells, real ambient profiles and real doublet partners. If the effect is large, this is the most persuasive motivating figure.
- Gene membership via minimum-effect (TREAT-style) FDR with random-effects pooling across expressing tumours and core/variable labels from I^2. This is a principled replacement for top-50 truncation.
- The practical warning, with quantification, that Poisson splitting of overdispersed data paired with a selection-heavy engine is badly anti-conservative (about 99% type I in the toy run), together with diagnostics for dispersion sign.

### Referee 2: statistical rigor and feasibility

**Predicted editorial outcome:** Likely sent to review, then major revision with a real risk (~50%) of reject-and-resubmit. The editor would send it out because an audit of the widely used 3CA metaprograms (MPs) and MSigDB GAVISH sets, wrapped around any engine, is timely. Referees would then find that several guarantees do not hold as specified: (i) the tumour-level null kappa>0 has no minimum effect; (ii) types are assigned by non-rejection; (iii) fold C is used twice; (iv) the estimand is defined incoherently. They would also find that the 'real biology' evidence is mostly self-consistency (decoys built from the inference model; artefact labels detected using the very covariates that define them), plus matched-protocol datasets too small to support AUROC claims. The 'new theory' (Theorem 1, Corollary 1) looks elementary and is probably already in the NB count-splitting paper; Hivert et al. report the opposite practical conclusion for marginal thinning under mixtures. To reach Nature Methods it needs: a minimum-effect, power-adjusted estimand; order-invariant typing based on interval or equivalence tests; calibration shown on real negative-control data; and orthogonal validation at n>=15-20 tumours per protocol contrast.

**Fatal flaws**

- TUMOUR MEMBERSHIP / EXISTENCE USE A POINT NULL WITH NO MINIMUM EFFECT, SO S(m), r_L, EXISTENCE AND TYPING TRACK POWER AND NUISANCE-MODEL ERROR, NOT 'EXPRESSES THE MP'.
- H0 is kappa_t(m)=0. With non-negative, dense consensus loadings, kappa_t(m)>0 holds whenever m shares genes with ANY program active in tumour t, or whenever any residual nuisance covariance remains.
- My toy run (exact NB thinning, fixed m, CF statistic): a tumour expressing a different program P (fold change 2 in 30% of cells) was declared to express m when m shared only 4/40 genes with P (cosine 0.1) in 9% / 18% / 20-24% of runs at n = 300 / 1000 / 3000 cells. At 8/40 shared genes the rate was 27% at n=300 and 47-52% at n=1000.
- The 'null-gene' empirical null did not reduce this leakage.
- Gavish MPs may share up to 10/50 genes, so leakage across MPs is the norm, not the edge case.
- Residual nuisance is amplified the same way: with equal-weight m over k genes, Z bias is about sqrt(n)·k·(per-pair residual correlation). A residual correlation of 1.2e-3 (n=300, k=50) or 1.6e-4 (n=1000, k=200) already shifts Z by 1. No ambient, doublet or CNA model is that exact.
- The 'specific' version residualises only on parent and siblings, not on all overlapping MPs.
- Gene membership has a minimum effect (eta_0); tumour membership does not. The two are inconsistent.
- TYPES ARE ASSIGNED BY THE FIRST NON-REJECTION, WHICH IS NOT ERROR-CONTROLLED AND DEPENDS ON RUNG ORDER.
- 'Exists at N0, fails at N1 → technical' treats absence of evidence as evidence. Adding covariates costs power through collinearity, extra parameters and fold-level covariate noise, so underpowered real programs (small or shallow tumours, sn data) get labelled technical, QC-entangled or protocol-sensitive.
- FDR controls false rejections, not false non-rejections. The claim of 'error control at each rung of the type ladder' is therefore false for the type labels themselves.
- Attribution follows a fixed order (N1→N2→N3→N4), like Type-I sequential sums of squares, so a program explained jointly by ambient RNA and dosage gets a type that depends on that order.
- FOLD C IS USED TWICE.
- S^N(m) is decided by the fold-C cross-fold (CF) test (step 5).
- Gene score tests are then run on fold C and pooled by random effects over t ∈ S(m) (step 7).
- Tumours are thus selected for high fold-C co-expression of m's genes, and the same fold-C counts are then used to estimate those genes' effects. This inflates eta_j and invalidates gene-level FDR, which is exactly the double dipping the design claims to remove.
- THE ESTIMAND IS INTERNALLY INCOHERENT.
- (a) The location is ENGINE-INDEXED (modes of F_delta over W^E_t), which contradicts 'the MP is the estimand; engines are estimators'. Different engines (centred-log Frobenius NMF, cNMF-KL, Gaussian EBMF) have different population programs, so wrapping cannot make them invariant; that claim is a hypothesis, not a guarantee.
- (b) delta is calibrated from finite-sample self-replicate cosine. At the stated 'infinite depth' population limit, delta→0, F_delta counts only exact matches, and with any biological variation between tumours no MP exists. The population estimand is degenerate unless delta is a fixed biological tolerance.
- (c) The 'independent thinning replicate' used to calibrate delta is not independent: two thinnings of the same matrix share about half their molecules. Self-cosine is therefore inflated, delta too small, and MPs over-split.
- (d) Existence requires >= s_min = 2 cohorts, yet 'single-cohort (fails s_min)' is listed as a subtype of MPs that exist at N4. This is a direct contradiction.

**Major concerns**

- THE NEW THEORY IS PROBABLY NOT NEW, AND THE SIGN CONTROL IS SECOND-MOMENT ONLY.
- Theorem 1: the NB count-splitting paper (Neufeld et al., arXiv 2307.12985) appears to already characterise fold variance and correlation for an arbitrary b'. The search snippet mentions a theorem with the ratio (b+1)/(b'+1) and positive correlation when overdispersion is too small. The algebra is three lines.
- Corollary 1 is the law of total variance.
- Hivert et al. (arXiv 2405.13591) report that marginal thinning that ignores mixture structure INFLATES type I error and loses power. The proposal cites this paper as being answered but never reconciles with its finding.
- Covariance <= 0 does not imply E[X^C | X^A] is decreasing in the upper tail, which is where selection engines act.
- My toy run (n=300, G=400, H0 true; per-gene marginals independent but not NB; moment-matched b~): mean Z was +0.14 for the on/off latent (about 2.3 SE; type I 0.067) and +0.20 for lognormal-Poisson (about 3.6 SE). With exact NB data I reproduced the authors' numbers: exact b 0.053, MoM 0.037, Poisson splitting 0.98.
- Conservatism appeared only when b~ was deliberately deflated (0.5×MoM gave mean Z −0.3), which costs power.
- Theorem 2's bound carries a factor |G| ≈ 8000 against L_t = O(1/n_t), so it is vacuous at n_t of a few hundred.
- Script and outputs: [referee_checks/ref2_sim.py](referee_checks/ref2_sim.py), `ref2_sim1.out`, `ref2_sim2.out`.
- 'b' IS NOT IDENTIFIED.
- Any NB is a Poisson-gamma mixture, so where cell state ends and NB noise begins is a modelling choice.
- Under H0, what matters is each entry's marginal law given nuisance, and real UMI marginals are not NB.
- The two-pass 'conditional' b~ is estimated from means that were themselves fitted. Overfitting at K up to 15 with 300 cells biases b~ UPWARD, toward the anti-conservative side, which the Corollary does not cover.
- POWER IS LIKELY INADEQUATE FOR THE DISCOVERY CLAIMS (Q4).
- Fold C holds 30% of molecules. The authors' own toy gives per-tumour power of 0.26-0.39 at fold change 1.5.
- Simes partial conjunction multiplies by T−r+1. A subset-shared MP carried by only about 5 tumours out of T=3000 needs p_(5) of roughly 1e-7 once multiplicity across about 150 MPs is included.
- Strict-mode e-BH with the 0.5·p^{-1/2} calibrator needs p ≤ 2.8e-4 per MP even when 100 of 150 are rejected.
- r_L counts rejections, so it scales with cells and UMIs per tumour. Comparisons of recurrence between cancer types and protocols are therefore confounded by depth, unless prevalence is estimated with power adjustment (e.g., empirical-Bayes deconvolution of tumour-level effect estimates).
- THE NUISANCE LADDER REMOVES BIOLOGY AND LEAVES RESIDUAL NUISANCE.
- N1 regresses on the same sample's TME profile. That absorbs genuine malignant programs that share genes with TME cells: partial-EMT vs CAFs (HNSCC), MES-like GBM vs myeloid cells, and IFN/MHC-II vs immune cells.
- The per-cell rho and pi are estimated at fold depth (30-50% of UMIs). That is an errors-in-variables problem, so some contamination co-expression remains, and the high-powered off-diagonal test will flag it as 'beyond nuisance'.
- The empirical null is matched on ambient weight but not on doublet-profile weight.
- N2 covariates (complexity, %mito, %ribo) track cell cycle, OXPHOS/hypoxia and MYC/translation.
- CNA/CLONE RUNGS ARE CIRCULAR AND PARTLY COLLINEAR.
- Allele-aware Numbat needs BAMs. For human tumours these are usually controlled-access (dbGaP/EGA), with per-dataset applications that take months, so N3/N4 will mostly rely on expression-only calls.
- Leaving out a 2 Mb window does not protect N4: clone labels are inferred from all genes, program genes included.
- Programs with clustered genes are at risk of being typed as dosage or clone effects: metallothioneins at 16q13 (the 'Metal' MP), the epidermal differentiation complex at 1q21 (squamous differentiation, on a commonly gained arm), keratin clusters at 17q21 and 12q13, MHC-II at 6p21, and histone clusters.
- When copy number is constant within a clone, cis (N3) and trans (N4) effects are separable only for genes inside altered segments.
- Clone indicators also absorb non-genetic states that differ between clones because of epigenetic memory or spatial niche.
- PROTOCOL TYPING IS CONFOUNDED WITH STUDY AND POWER.
- Protocol is nested within study (lab, tissue handling, dissociation time).
- Permuting protocol labels within cancer type treats tumours from the same study as exchangeable, which is anti-conservative. Permutation must be at study level.
- sn data have fewer UMIs, so fewer tumours are rejected and MPs get typed 'protocol-sensitive' through power alone.
- VALIDATION LARGELY SHOWS SELF-CONSISTENCY, NOT REAL BIOLOGY.
- Tier 1/2 decoys are simulated from the fitted N4 inference model plus misspecification the authors choose. Decoys contain no real programs, so hybrid, split and overlap MPs (the main real-data failure) are never counted. NB(nu^{N4}, b_marg) with per-cell means also double-counts variance.
- Tier 3 is tautological. 7/11 hand-removed MPs are LQ (ribo/mito) and are 'detected' by an N2 rung that uses %mito/%ribo as covariates; 3/11 are contamination MPs detected by an N1 rung that uses TME profiles. With 11 positives, an AUROC of 0.85 has a 95% CI of about 0.70-1.0 (Hanley-McNeil).
- Tier 5: artefacts (ribo/mito, stress) replicate across cohorts too, so replication does not discriminate.
- Tier 4: the matched-protocol data are one Janesick breast block plus a few Oliveira CRC patients, too few for an AUROC. Xenium panels may lack key stress genes (unverified per panel). In Tahoe-100M, QC and translation programs also respond to drugs through viability effects.
- SCOPE UNDERMINES THE 'ENGINE-AGNOSTIC' PITCH.
- The core test is within-tumour, so GBCD's distinctive output, between-tumour subtype programs, is undetectable without the optional pseudobulk extension.
- GBCD works on an n×n cell covariance. Running it per cancer type at 3CA scale (around 1e5 cells) times 5 replicates times 4 cohorts is infeasible without subsampling, which weakens the 'robust under >= 2 engines' criterion.
- RANK SELECTION BIAS.
- K is chosen by held-out deviance on fold B, which holds 20% of molecules (about 600 UMIs per cell at a 3k median). Lower depth favours smaller K, so programs get merged into hybrids, the very failure being targeted.
- The universal refit is log-linear in usages, while NMF is additive. The deviance therefore scores a different model from the engine's.
- An SE taken over thinning replicates of the same cells understates sampling variability, so the 1-SE rule is ill-founded.
- THRESHOLDS ARE RELOCATED, NOT REMOVED.
- About 20 tunables: eps×3, r_min, s_min, eta_0, delta percentile, rho, q×3, K grid, 1-SE, top-3,000 / 10% / 8,000 universe rules, residual clipping, null-gene definition, 1,000 pseudo-programs, R, 50 draws, pi_m >= 0.7, and >= 2 engines.
- That compares with about 6 in Gavish.
- No sensitivity analysis of MP count or identity to these knobs is planned.
- DERANDOMISATION IS UNDER-SPECIFIED.
- Each thinning replicate tests a different m. Cauchy-combining p-values across those replicates tests a union of different hypotheses.
- The combination is valid only asymptotically, in the tail, under strong dependence, while BH here operates at 0.05-0.1.
- The 'delta-neighbourhood hypothesis' is never formally defined.
- GENE EFFECTS ARE DEPTH-DEPENDENT.
- eta_j is a log fold change per SD of a noisy fold-A activity score, so it suffers regression dilution that worsens in shallow or small tumours.
- 'Per SD' also mixes how many cells use the program with how strongly it is expressed.
- The minimum-effect filter therefore removes genes depending on the depth mix of the cohort.
- COMPUTE IS UNDERESTIMATED.
- The method's own default K grid is 2-15 (14 ranks vs Gavish's 6), with up to 3 monotone-ascent passes, 5 replicates and 4 cohorts (real plus 3 decoys). That is about 2,800-8,400 NMF fits per tumour, roughly 2-7× the 2×10^4 CPU-hour estimate for ONE engine. The per-fit timing behind this is my rough estimate and unverified.
- '50 cheap perturbation draws reusing engine fits' contradicts drawing new thinning and engine seeds.
- The benchmark (7 engines × 27 grid cells × planted types × 8 ablations) is a project in its own right.
- DATA AVAILABILITY IS OPTIMISTIC.
- Public filtered matrices rarely include empty droplets.
- Many 3CA tumours lack TME cells (sorted samples) or are Smart-seq2 (no UMIs).
- TCCA (36 studies, 537 patients) almost certainly overlaps heavily with 3CA v2's 124 datasets, so the independent replication set may be small. This should be counted before replication is claimed.
- THE TIMELINE IS UNREALISTIC.
- 17 months with about 3.2 FTE is planned to cover: the package, theory including an unproven Theorem 2, a 7-engine benchmark, a re-derived atlas of >= 1,000 tumours, TCCA replication, four orthogonal validation tiers and an optional wet-lab test.
- A realistic estimate for the full scope is 30+ months, or the scope must be cut.

**Prior-art threats**

- Neufeld et al., 'Negative binomial count splitting for single-cell RNA sequencing data', arXiv 2307.12985: the search snippet describes a theorem giving split-count variance under an arbitrary b', with the ratio (b+1)/(b'+1), and positive fold correlation when overdispersion is too small. Theorem 1 is therefore likely not new. [VERIFIED that the paper exists; theorem content PARTLY VERIFIED from the snippet; full text blocked]
- Hivert et al., 'Practical limitations for real-life application of data fission and data thinning in post-clustering differential analysis', arXiv 2405.13591: marginal thinning that ignores mixture structure inflates type I error and lowers power; the paper proposes conditional thinning. This undercuts Corollary 1's practical conclusion that marginal dispersion is conservative, and must be reconciled. [VERIFIED via search abstract snippet]
- Neufeld, discussion of 'Data fission: splitting a single data point', arXiv 2409.03069: discusses NB fission with an estimated overdispersion. [VERIFIED that it exists; content beyond the snippet UNVERIFIED]
- Langfelder, Luo, Oldham, Horvath, 'Is my network module preserved and reproducible?', PLoS Comput Biol 2011: established statistics for whether a gene module exists in another dataset (Zsummary). This undercuts the novelty of audit mode, which must cite it and be benchmarked against it. [VERIFIED]
- Bogomolov & Heller, 'Replicability Across Multiple Studies', Statistical Science 38(4), 2023, plus arXiv 2511.04130 on partial-conjunction (PC) testing across dependent studies: an established r-out-of-n replicability framework, so the PC recurrence bound is an application rather than an invention. [VERIFIED that both exist]
- GBCD (Liu…Stephens, Nature Genetics 2025): shared vs patient-specific programs with gene-level lfsr. [VERIFIED in landscape]
- ITHmapper / pan-cancer ITH compendium (Abelson lab, bioRxiv 2026): pan-cancer replacement for per-sample NMF, with a mapping tool. [VERIFIED in landscape]
- TCAT/starCAT (Nature Methods 2025): reproducible program catalogue plus projection. [VERIFIED in landscape]
- ClusterDE synthetic nulls and recall knockoffs: the closest prior art to the decoy idea; CALIPER's distinctive step is removing co-expression. [VERIFIED in landscape]
- CanSig benchmark (Cancer Research 2026): the existing tumour-state benchmark harness. [VERIFIED by one sweep only]
- CSI-GEP (Cell Genomics 2025): replication of programs across data splits. [VERIFIED in landscape]
- MetaNeighbor-style cross-dataset replicability (Crow & Gillis) and AdaFilter-type power improvements for multiple PC tests: likely reviewer citations. [UNVERIFIED, from memory]
- Name check: a quick web search found no single-cell package called CALIPER. This is weak evidence, because 'caliper' is a common statistical term (e.g., propensity matching). [PARTLY VERIFIED]

**Required changes**

- Redefine tumour expression with a MINIMUM EFFECT that is specific to ALL other MPs, not just siblings. Examples: kappa_t(m | M_{-m}) >= kappa_0, or a fraction of m-gene residual variance explained >= v0. Report tumour-level effect estimates with CIs, and estimate prevalence by power-adjusted empirical-Bayes deconvolution. Keep r_L as a secondary quantity.
- Replace first-failure typing with order-invariant attribution (all-subsets or Shapley decomposition of kappa across nuisance classes), with bootstrap CIs. Use equivalence or interval tests ('at least X% explained by N1') and an explicit 'indeterminate / underpowered' type. Never type by non-rejection.
- Separate folds for membership and gene inference: decide S(m) on fold B and test genes on fold C, or use rotated 4-fold cross-fitting with derandomisation. Re-derive the power budget accordingly.
- Fix the estimand. Use a fixed, depth-independent biological tolerance delta. Calibrate self-similarity on independent folds within one thinning, not across thinning replicates. State plainly that locations are engine-indexed and make engine invariance a tested hypothesis. Formally define the derandomised (delta-neighbourhood) hypothesis. Remove the contradiction between s_min and the 'single-cohort' subtype.
- Reposition the theory. Cite and compare with the NB count-splitting result on misspecified b' and with Hivert et al., and present Theorem 1 and Corollary 1 as lemmas. Either prove Theorem 2 under explicit assumptions (e.g., monotone-likelihood-ratio thinning and bounded engine influence) with a non-vacuous bound, or drop it in favour of empirical calibration.
- Show calibration beyond NB: lognormal-Poisson, bimodal, ZINB-like and real-data-resampled marginals, at realistic n_t (100-2,000), G (8,000) and selection engines (NMF top-k, cNMF). Report the mean-Z shift and the type I error, not only the type I error (my toy runs show +0.14 to +0.20 shifts with moment-matched b~). Add calibration on real negative controls, e.g., homogeneous cell-line or replicate data with no planted programs, and protocol-matched technical replicates.
- Make the validation non-circular. Evaluate N1/N2 typing on labelled MPs whose defining genes are not themselves the covariates (or with those covariates held out). Report CIs given 11 positives. Use matched-protocol cohorts with >= 15-20 tumours per contrast (fresh vs Flex/fixed vs sn, through collaborators). Permute at study level. Obtain DNA- or allele-defined clones in >= 10 tumours. Pre-register spatial-niche and perturbation predictions that technical programs should fail.
- Add a decoy variant that keeps real programs while removing only the planted artefact, so that false hybrid, split and overlap MPs are counted. Report the FDR estimate's sensitivity to the misspecification family xi.
- Model protocol and cohort effects with study-level random effects and power-adjusted membership (effect sizes with SEs, not counts of rejections).
- Run the data audit before committing. Count: tumours with integer UMI counts; tumours with TME cells; samples with raw or unfiltered barcodes; BAM accessibility (open vs dbGaP/EGA); and TCCA studies independent of 3CA.
- Re-cost compute using the method's own defaults (K 2-15, ascent passes, replicates, decoys). Specify how GBCD is subsampled. Remove the 'cheap perturbation draws' claim, or restrict the draws to precomputed fits.
- Add a sensitivity analysis over all ~20 tunables, showing how the number and identity of MPs move.
- Rescope for about 24 months with the current team: an audit mode (benchmarked against module-preservation statistics), the TumourNull benchmark, two engines (Gavish, cNMF), and ONE flagship question (the stress decomposition with adequately powered matched-protocol data). The remainder becomes follow-up work.

**Ideas worth keeping**

- Separating estimand from estimator, with an engine-agnostic inference layer and adapters. This is a low-friction adoption path and the right framing for Nature Methods.
- Audit mode for fixed external catalogues (3CA, Gavish, Kinker, Barkley) on independent cohorts. With m fixed, the off-diagonal within-data statistic is valid. It is the cheapest and most immediately useful deliverable, provided it beats module-preservation statistics.
- Cross-fold off-diagonal co-expression statistics, the rule never to pair the same matrix entry across folds, and deliberately conservative dispersion. My toy run confirmed calibration under NB (type I 0.053 with exact b, 0.037 with MoM) and the near-total failure of Poisson splitting with selection-heavy engines (0.83-0.98).
- Partial-conjunction replicability at both tumour and cohort level. Keep it, but with effect-size thresholds and power-adjusted prevalence.
- The validity-stability plane: it reports stability separately from null calibration, so 'stable artefact' becomes a visible quadrant.
- Running the unmodified Gavish/GeneNMF/scalop pipelines on nuisance-only decoys (pilot P2) as a cheap, decisive test that determines the headline.
- Random-effects pooling of gene effects with a TREAT-style minimum effect and I^2-based core/variable/anchor labels, in place of top-50 truncation.
- Pre-registered go/no-go pilots and a full run report (parameters, seeds, versions).


---

## 3. InvarMP (Invariance-typed MetaPrograms). Python package `invarmp`, with an R/Seurat wrapper `InvarMPr`. Working name only; check for name collisions before release.  
*Design angle: `invariance`*

### Claim

Because no dissociation, preservation or readout protocol is free of artefacts, InvarMP defines a tumour metaprogram as a gene program in counts whose loadings and usage stay the same across technical environments (cell vs nucleus, fresh vs fixed Flex, dissociated vs in situ, chemistry, site), after one shared measurement transport. It learns these programs jointly with programs owned by a technical environment and with copy-number dosage offsets. Matched multi-platform tumours, sensors read from tumour-microenvironment (TME) cells and DNA/allele clones then give every metaprogram an invariance profile and type with error control, which replaces manual curation at pan-cancer scale.

### Formal definition of an MP

SETTING
- Each malignant cell c has:
  - a tumour t(c) and a study/site s(c);
  - a TECHNICAL environment e(c). This is a tuple of attribute levels: preservation {fresh, cryopreserved, frozen, FFPE}; molecule {whole cell, nucleus, in situ}; dissociation {warm enzymatic 37C, cold protease, fixed-then-dissociated (Flex), none (in situ)}; chemistry/readout {3'v2, 3'v3, 5', Smart-seq2, Flex probes, panel P (Xenium 313/422/5K, CosMx 6K)}; site.
  - a BIOLOGICAL environment b(t): cancer type, primary or metastasis, in vivo or in vitro, treatment, spatial niche;
  - a GENOTYPE: clone q(c) and copy-number vector CN_c. Clones come from DNA via clonealign, or from alleles via Numbat or CalicoST.
- Environment attributes may be missing for a sample. They are then latent variables with a posterior (see model).
- A program is a non-negative loading vector w on the gene simplex, with per-cell usage theta_c >= 0.

DEFINITION
A program w with a tumour-membership set T_w is an INVARIANT METAPROGRAM if all of the following hold.
- (D1) Recurrence. T_w contains at least 2 tumours from at least 2 studies. This is necessary but no longer sufficient.
- (D2) Loading invariance. For every technical environment e in which w is testable, the program's expected contribution to counts is kappa_e ⊙ w.
  - kappa_e is a gene-wise measurement operator: capture, probe or nuclear-retention efficiency.
  - kappa_e is SHARED by all programs and by the non-malignant calibrator cells captured with them.
  - Equivalently, the program-specific deviation is zero: delta_w^e = 0.
- (D3) Usage invariance.
  - In matched designs (the same tumour or block profiled in environments e and e'), the within-tumour usage law is independent of the environment: p(theta | t, e) = p(theta | t, e'). The fraction of active cells (pi) and the active intensity (m) are tested separately.
  - In unmatched designs, tumour-level usage is conditionally independent of the technical anchors given biology: rbar_{w,t} ⟂ A_t | B_t. Here A_t is technical metadata plus TME-derived environment sensors, and B_t is biological covariates.
- (D4) Not reducible to cis dosage. The program survives a per-gene copy-number offset (CN/2)^gamma_g.

Invariance across BIOLOGICAL environments and across clones is NOT required. It is measured and reported.

INVARIANCE PROFILE
For each MP and each technical axis a (cell↔nucleus, fresh↔Flex, dissociated↔in situ, warm↔cold, 3'↔5', Smart-seq2↔10x, site), the profile reports:
- transported loading cosine with credible interval (CI), and an e-value;
- selection shift Δpi and induction shift Δm, each with CI and e-value;
- anchor sensitivity (partial R²_A);
- testability flag, from the design graph and panel information.

It also reports:
- genotype dependence: rho_gen = mean over tumours of Var(E[theta | clone]) / Var(theta); chromosome-arm positional enrichment; cis-dosage share;
- biological context: posterior tumour and cancer-type breadth; presence in vitro (CCLE, Tahoe DMSO); spatial coherence;
- loading heterogeneity across biological strata (MP family vs context variants);
- dispersion across seeds.

TYPES
- T1 Invariant plastic: passes D2 and D3 on every testable axis; low rho_gen.
- T2 Invariant clone-linked: passes D2 and D3; rho_gen is high with a significant clone-permutation e-value. Subtype T2p (positional) if residual arm enrichment is significant.
- T3 Protocol-modulated biology: loadings invariant, but usage shifts on an axis. Subtype "selection" if Δpi ≠ 0 (cells lost or enriched by the protocol). Subtype "induction" if Δm ≠ 0 (intensity amplified, e.g., in vivo stress boosted at 37C).
- T4 Context-restricted biology: invariant across technical environments, but usage is concentrated in one biological stratum (e.g., in vivo only, one lineage, one niche).
- T5 Technical artefact: the component is owned by a technical-environment attribute, OR it fails D2 with transported cosine < 0.5 against the reference environment. Sub-labels come from alignment with known nuisance:
  - ex vivo immediate-early genes, using the 512-gene collagenase stress core set (O'Flanagan 2019);
  - nucleus-specific;
  - ambient;
  - heterotypic doublet/contamination;
  - low quality.
- T6 Untestable/confounded: the design graph gives no bridge across technical environments for the tumours using the program. It is reported, never silently kept or dropped.

The principle in one line: biology is what stays fixed when only the measurement changes. Artefacts belong to environments, and clone programs belong to genotypes.

### Model and math

1) LIKELIHOOD
The model is fitted to malignant cells together with subsampled non-malignant "calibrator" cells from the same samples.

X_cg ~ NB(mu_cg, phi_{g,e(c)}), for g in G_{e(c)}. Genes outside a panel are MISSING, not zero.

mu_cg = l_c * [ kappa_{g,e(c)} * D_cg * ( SUM_{k in S} theta_ck w_kg exp(delta_kg^{e(c)}) + SUM_{j in E} theta_cj v_jg 1[lev(j) in e(c)] + SUM_{m in P_t(c)} theta_cm u_mg + SUM_{r in C} theta_cr h_rg ) + alpha_c b_{g,s(c)} ]

Libraries and terms:
- S: shared library of candidate MPs. Truncation 200, with automatic relevance determination (ARD).
- E: environment-owned library. Each component v_j is attached to ONE attribute level (e.g., dissociation=warm, molecule=nucleus, chemistry=5', site=X). It may be used only by cells whose environment has that level. Default: 4 components per level.
- P_t: tumour-private library, 3 per tumour. This absorbs patient-specific axes; it is not the novelty.
- C: cell-type library. It explains the calibrator cells and supplies heterotypic-doublet profiles for malignant cells, with a theta_cr prior set from scDblFinder scores.
- alpha_c b_{g,s}: the sample's ambient profile from empty droplets. For Xenium/CosMx, this is replaced by a neighbour-spillover profile.
- D_cg = (CN_{c,seg(g)}/2)^{gamma_g}: cis-dosage offset with gamma_g ~ LogNormal(log 0.8, 0.5). The prior centre is a design default and is re-estimated. CN comes from allele-aware Numbat/CalicoST calls, or from DNA via clonealign where available.

2) MEASUREMENT TRANSPORT
- log kappa_{g,e} = SUM_a f_a(x_g) + eps_{g,e}, with eps ~ N(0, s²).
  - x_g are gene features: length, intron count, GC, probe count, nuclear-retention index.
  - f_a are smooth functions per attribute.
- Reference environment e0 = fixed-before-dissociation (Flex-FFPE, whole-transcriptome probes): log kappa_{.,e0} = 0.
- Identification rests on assumption A1: TME identity programs (T cell, fibroblast, endothelial, myeloid) are environment-invariant. kappa is therefore pinned by the calibrator cells, never by the malignant programs under test. A1 is checked on matched data.

3) INVARIANCE STRUCTURE
- Program-specific loading deviation: delta_kg^e ~ Horseshoe(tau_delta), with delta^{e0} = 0. D2 holds iff delta_k^e = 0.
- Hurdle usage: theta_ck = z_{k,t} * h_ck * r_ck
  - h_ck ~ Bern(sigmoid(logit pi_{k,t} + eta^pi_{k,e}))
  - r_ck ~ Gamma(a_k, a_k / (m_{k,t} exp(eta^m_{k,e})))
  - eta ~ Horseshoe
  - eta^pi ≠ 0 means SELECTION; eta^m ≠ 0 means INDUCTION.
- Twin tie: pi_{k,t} and m_{k,t} are shared by every environment in which tumour or block t was profiled.
- Tumour membership: z_{k,t} ~ Bern(p_k), p_k ~ Beta(a0/K, 1) (finite Indian buffet process). This gives credible tumour sets.
- Biological strata get their own sparse delta^{b}, which captures MP context variants. It is reported, not penalised.

4) OBJECTIVE
L = E_q[ zeta * log p(X | ·) ] − KL(q || prior) − lambda_A * SUM_k || Pi_A rbar_k ||²

- zeta in (0,1] is a power posterior (cf. BayesPowerNMF). It suppresses junk components that arise from misspecification.
- rbar_k is the tumour-level mean usage of program k, residualised on biology B_t.
- Pi_A projects onto the anchor column space.
- This is anchor regression (Rothenhäusler, Bühlmann, Meinshausen, Peters, JRSS-B 83:215–246, 2021 [V]) transplanted to programs:
  - lambda_A = 0 gives the plain fit;
  - large lambda_A favours shared programs whose tumour-level usage is robust to shifts along technical anchors.
  - A formal distributional-robustness guarantee for NB-NMF is open and is a theory aim.
- zeta and the truncation are chosen by held-out NB deviance on a count-split fold.

5) PARTIALLY OBSERVED ENVIRONMENTS
- A missing attribute level for sample s is latent: q(e_s) ∝ prior(metadata) × p(F_s | e_s).
- F_s is a TME SENSOR fingerprint computed ONLY from that sample's non-malignant cells:
  - (i) nascent IEG index: unspliced/(spliced+unspliced) reads over about 50 immediate-early genes (FOS, JUN, EGR1, ATF3, ...) in T and stromal cells. HYPOTHESIS: transcription in the minutes of ex vivo dissociation is rich in unspliced reads. This is to be validated against 4sU labelling.
  - (ii) intronic read fraction;
  - (iii) mitochondrial fraction;
  - (iv) cosine similarity to the ambient profile;
  - (v) gene-length bias of the calibrator pseudobulk.
- p(F|e) comes from a classifier trained on samples whose attributes are known.
- Environment components enter the model through E_q[1[lev in e_s]].
- Malignant cells never inform their own environment (no circularity) unless joint mode is chosen.

6) TESTS (per shared program k and axis a)
(a) Screening.
- The IRMv1 penalty is the squared gradient with respect to an environment-specific scaling at the invariant solution. Here that gradient is exactly the Rao score for delta at 0:
  S_kg^a = SUM_{c: a(c)=a1} l_c kappa D theta_ck w_kg (X_cg − mu_cg) / (mu_cg + mu_cg²/phi)
  T_k^a = S' I^{-1} S
- This needs one backward pass and no refit.

(b) Confirmation by thinned split-likelihood-ratio e-values (universal inference [U]).
- Split X = X1 + X2 by data thinning (binomial for Poisson; the Neufeld NB scheme with phi estimated on an independent fold).
- Fit the unrestricted model (delta_k^a free) on X1.
- On X2, with every other component and w_k fixed at their X1 estimates, the restricted model (delta=0) is a CONVEX Poisson problem in usage. This guarantees the global null MLE that validity needs.
- E_k^a = L_X2(unrestricted from X1) / L_X2(restricted MLE). Then E[E] ≤ 1 under H0.
- Validity is exact for Poisson and approximate for NB (Hivert et al. caveat). Calibration is checked with "null twins": one environment thinned into two pseudo-environments.
- eta^pi and eta^m are tested the same way.

(c) Usage invariance in the unmatched atlas.
- Test rbar_{k,t} ⟂ A_t | B_t with study random effects, cross-fitted on halves of the studies.
- This is the invariant-causal-prediction logic (Peters, Bühlmann, Meinshausen, JRSS-B 2016 [V]) with programs in place of predictors.

(d) Genotype.
- rho_gen with a size-preserving clone-permutation null;
- arm enrichment of w_k against expression-matched gene sets;
- cis share = fraction of counts absorbed by D.

(e) Biological context.
- breadth (number of tumours and cancer types with posterior z = 1);
- projection onto Kinker CCLE and Tahoe DMSO cells;
- Moran's I on Xenium/CosMx.

Multiplicity: e-BH at q = 0.1 over all (k, axis) pairs [U]. E-values from different evidence sources for the same pair are averaged, which stays valid under any dependence.

7) IDENTIFIABILITY (theory aim, to be proved)
- Design graph G: a bipartite graph linking technical environments to biological units (blocks, tumours or cancer-type strata).
- Proposition. It holds if three conditions hold:
  - (i) kappa is identified from calibrator cells (A1);
  - (ii) the units using program k span at least 2 environments connected in G that differ in attribute a;
  - (iii) usage is sufficiently scattered within each environment, extending the noiseless results of Fu, Huang and Sidiropoulos.
  Then the split of k into "shared" vs "owned by a" is unique up to permutation and scale.
- If (ii) fails, the program is reported as untestable on a (type T6).
- Defaults pending calibration in simulation:
  - T2: rho_gen ≥ 0.3 and clone e-value ≥ 20;
  - T5: transported cosine < 0.5.

### Inference and scalability

ENGINE
- PyTorch + Pyro stochastic variational inference (SVI) on GPU.
- Amortised non-negative (softplus) encoder q(theta_c | log1p x_c on G_e, e_c, t_c). Tumour-membership and environment-ownership masks are applied to the encoder outputs.
- Variational families:
  - loadings: mean-field LogNormal/Dirichlet;
  - tumour membership z: relaxed Bernoulli (straight-through);
  - horseshoe priors: non-centred.
- Minibatches of 4,096 cells, stratified by tumour × environment. Rare environments (Flex, in situ, cold protease) are up-weighted with importance weights so the ELBO stays unbiased.

DATA HANDLING
- Counts are streamed from Zarr-backed AnnData, MuData or SpatialData (lazy, dask).
- Gene universe: genes detected in ≥1% of malignant cells in ≥10% of tumours (about 10–12k), plus all panel genes. Panels enter as gene masks.
- Calibrator TME cells: at most 500 per cell type per sample.

COST (TARGETS, to be benchmarked)
- About 1.5×10^3 non-zeros per cell, 10^7 cells and about 60 active components per cell (ARD-truncated) give about 10^12 multiply-adds per epoch.
- Target: ≤20 min per epoch and 30–60 epochs on one 80 GB A100/H100, so ≤1 GPU-day for a 3CA v2 + TCCA-scale fit.
- Tumour-private loadings: 3 × about 3,000 tumours × 12k genes ≈ 1.1×10^8 floats ≈ 0.4 GB.

MODEL SELECTION
- Grid over zeta ∈ {0.5, 0.7, 0.85, 1} × truncation, scored by held-out NB deviance on a count-split fold.
- Each grid point is fitted on a 10% cell subsample, then one full fit is run at the chosen setting.
- 5 seeds. Components are consensus-matched across seeds (Hungarian matching on cosine), and the dispersion is reported in each MP profile.

TESTS
- Score screens: a single backward pass.
- Split-LR e-values: one additional fit on thinned folds (2× cost), plus convex restricted refits only for screened candidates (about 200 programs × 8 axes = 1,600, seconds each on GPU).
- Anchor/conditional-independence and clone tests run at tumour or clone level (thousands of rows), so they are negligible.

PROJECTION
- An encoder forward pass, then a few steps of per-cell Poisson-NNLS refinement. Minutes for 10^5–10^6 query cells.
- Panel queries (Xenium 313/5K, CosMx 6K) use masked, kappa-transported loadings, with uncertainty inflated in proportion to panel Fisher information: I_k(P) ∝ SUM_{g in P} (kappa_g w_kg)² / mu_g.

SCALING BEYOND MATCHED DATA
Five mechanisms let the method scale past the matched sets:
- (1) TME calibrators identify kappa_e from thousands of unmatched samples.
- (2) TME sensors act as anchors for invariance tests across the whole atlas.
- (3) Natural experiments: studies mixing sc/sn, 3'/5' or sites; Smart-seq2 vs 10x within 3CA.
- (4) Online addition of new environments (e.g., CosMx 6K) with frozen loadings and re-estimated deviations.
- (5) Thinning-generated null twins, which make the invariance null available everywhere.

### Outputs for users

1. Versioned MP catalogue (Zenodo DOI; semantic versions). For each MP:
- loadings in the reference environment (Flex-equivalent units), with 90% credible intervals;
- local false sign rate (lfsr) per gene, and FDR-controlled gene sets (lfsr < 0.05) in place of top-50 lists;
- kappa-transported loadings for each platform/panel (3' fresh, sn, Flex, Visium HD bins, Xenium 313/422/5K, CosMx 6K), with a panel-informativeness score.

2. Invariance profile card per MP (HTML and JSON). For each technical axis:
- transported loading cosine with CI;
- selection Δpi and induction Δm with CIs;
- e-values and testability.
Also on the card:
- anchor sensitivity R²_A;
- genotype block: rho_gen, positional arm enrichment, cis-dosage share;
- biological block: credible tumour and cancer-type breadth, in vitro presence, spatial Moran's I;
- context-variant loadings;
- seed dispersion;
- TYPE label T1–T6 with confidence.

3. Per tumour: posterior membership P(z_{k,t} = 1), plus tumour-level usage pi_{k,t} and m_{k,t} with CIs.

4. Per cell (written to AnnData):
- usage posterior mean and SD (.obsm['invarmp_usage'], ['invarmp_usage_sd']);
- ENVIRONMENT-STANDARDISED usage: the counterfactual usage in the reference environment, with environment-owned components and eta shifts removed (.obsm['invarmp_usage_std']);
- ex vivo/artefact burden: the fraction of counts explained by environment-owned components;
- cis-dosage contribution;
- ambient and doublet fractions.

5. Per sample:
- posterior over imputed environment attributes;
- TME sensor fingerprint and anchor values;
- QC flags.

6. Artefact library: the environment-owned components, each labelled with its owning attribute (e.g., "warm-dissociation immediate-early", "nucleus-specific", "5'-chemistry"). This reusable nuisance model replaces gene blocklists.

7. Design-graph report: which axes are testable for which MPs and cancer types, and where protocol is confounded with biology.

8. `invarmp.project()`: maps a new dataset or platform into reference units, with per-cell uncertainty and novelty detection. Novelty is a structured-residual e-value against a thinning null.

9. Environment simulator and benchmark harness: planted program types, nuisance-only null cohorts, metrics.

10. Exports: AnnData/MuData/SpatialData, a Seurat object via the R wrapper, and gene sets in GMT format with confidence tiers for MSigDB-style use.

### Failure modes fixed

GAVISH/3CA PIPELINE STEPS

- S0 Malignant calling (contamination and doublet MPs).
  - Fix: calibrator TME cells are modelled together with the malignant cells. Heterotypic-doublet components (from the same sample's cell-type library) and ambient components are explicit.
  - Result: macrophage, T-cell and endothelial "MPs" fall into nuisance and are typed T5-contamination.
  - Target: automatically flag the 3 contamination MPs that 3CA v2 removes by hand (inds_rm).

- S1 Gene universe (top 7,000 by mean) and per-tumour centring with negatives zeroed.
  - Fix: an NB likelihood on all detected genes with size factors, and no centring. Tumour-level usage (pi, m) keeps programs that are active across a whole tumour.
  - kappa transport makes loadings comparable across chemistries, nuclei and panels. Panel genes are treated as missing, not zero.

- S2 Fixed ranks 4–9, best of 10 restarts.
  - Fix: ARD plus a finite Indian buffet process.
  - Truncation and power zeta are chosen by held-out count-split deviance. A 5-seed consensus is run and its dispersion reported per MP.

- S3 Top-50 truncation.
  - Fix: posterior loadings, lfsr and FDR-controlled gene sets, plus panel informativeness per platform.

- S4 Recurrence filters (35/10/10).
  - Fix: recurrence is kept as necessary (D1) but is no longer sufficient. Posterior tumour membership plus invariance e-values replace the thresholds.
  - This addresses the core flaw directly. A warm-dissociation program recurs in every fresh dataset, so it passes 10/50 trivially. It is absent from, or owned away from, Flex, in situ and cold-protease data, so it fails D2/D3 or lands in the environment-owned library.

- S5 Greedy, founder-dependent clustering.
  - Fix: one joint shared library across all tumours, with no clustering or ordering.
  - Cohort dependence is quantified as type stability under leave-one-study-out.

- S6 Manual curation (inds_rm; GeneNMF dropMetaPrograms) and CNA checked only afterwards.
  - Fix: invariance-based typing with e-BH replaces curation. The 11 hand-removed 3CA v2 MPs and the protocol-named MPs in the 2026 catalogue ("Cell cycle single-nucleus", "Stress (in vitro)") serve as labelled tests.
  - CNA enters discovery as a cis-dosage offset. Trans clone-linked states are typed T2 rather than deleted.

- S7 Scoring without uncertainty or comparability.
  - Fix: environment-standardised usage with posterior SD, and projection into Flex/Xenium/CosMx with uncertainty specific to each panel.

CONCEPTUAL FAILURES (landscape section 1.4)
- (1) Recurrence ≠ reality: replaced by invariance across technical interventions.
- (2) No likelihood: NB model, e-values and credible sets.
- (3) Per-tumour centring and factorisation: joint model.
- (4) Biased stability filters: removed. Stability is reported, not used as a selector.
- (5) Cohort dependence: design-graph-aware tests plus stability analysis.
- (6) Post-hoc validation: environment, genotype and matched-platform data are used INSIDE the definition.

NEW FAILURE MODE NAMED AND HANDLED
"Protocol confounding": an artefact coincides with a cancer type that was processed with one protocol (e.g., single-nucleus mostly in brain tumours). InvarMP reports it as T6 untestable instead of calling it biology.

### Validation plan

TIER A: SEMI-SYNTHETIC ENVIRONMENT SIMULATOR (ground truth with failure modes)
- Backgrounds are real malignant cells from matched twin sets, so kappa_e and the noise of each environment are empirical:
  - Trinks et al. 2024: lung, patient-matched fresh / cryopreserved 3'v3 / FFPE Flex [V; venue P];
  - ccRCC bioRxiv 2025.11.25.690485: fresh sc / frozen sn-Multiome / FFPE snFlex [V];
  - Slyper 2020: sc/sn [V];
  - Oliveira 2025: Flex [V].
- 20 planted programs per cohort:
  - 6 invariant shared programs, each used by 10–60% of tumours;
  - 3 protocol-modulated (selection Δpi = ±0.2; induction 1.5–3×);
  - 3 environment-owned (warm immediate-early genes with Gamma per-cell intensity; nucleus-specific; 5'-specific);
  - 2 clone-linked (tied to Numbat clones);
  - 2 positional (planted single-copy arm gains);
  - 2 context-restricted;
  - 2 tumour-private.
- Nuisance: ambient 2–20%, heterotypic doublets 2–8%, depth shifts.
- Grid:
  - matched-twin fraction {0, 1, 5, 20%};
  - protocol–cancer-type confounding, Cramér's V {0, 0.3, 0.6, 0.9};
  - panel size {313, 1k, 5k, whole transcriptome};
  - missing-environment-label rate {0, 30, 70%};
  - 50 replicates per cell of the grid.
- Null cohorts contain nuisance only.
- Targets:
  - type macro-F1 ≥ 0.8 when V ≤ 0.6 and twins ≥ 5%;
  - empirical FDR of non-invariance calls ≤ 0.1 at q = 0.1;
  - ≤ 0.1 T1 calls per null cohort on average;
  - coverage of 90% tumour-membership credible sets between 0.85 and 0.95;
  - ≥ 0.9 correct "untestable" calls when bridging fails.
- Report where every baseline converts artefacts into MPs.

TIER B: HELD-OUT TWIN PREDICTION ("invariance predicts the other platform")
- Fit on the atlas plus twins with the second environment hidden. Predict each MP's Flex or in situ loadings and usage, then compare with the observed data.
- Targets:
  - AUROC ≥ 0.85 for predicting MPs that vanish (transported usage < 25%) in fixed or in situ data;
  - Spearman ≥ 0.6 between predicted and observed Δpi and Δm.
- Data:
  - Janesick 2023 Nature Communications 14:8353, GSE243280: one breast FFPE block; Flex, Visium CytAssist, Xenium 313-plex [V];
  - Oliveira 2025 Nature Genetics: CRC FFPE; Visium HD n = 5, Xenium n = 4 (322 + 100 genes), Flex atlas of 260,506 cells from n = 8 (GSE280311) [V];
  - Nature Communications 2025 (s41467-025-64292-3): COAD, HCC and OV with scRNA-seq on the same samples plus Xenium 5K, CosMx 6K, Visium HD and Stereo-seq on serial sections [V; patient n not retrieved];
  - Nature Methods 2026 (s41592-026-03089-8): Xenium on 19 breast + 22 NSCLC FFPE sections including 5K, with matched snRNA-seq for a subset [V];
  - Yan…Marra, Frontiers in Genetics 2026: 12 matched FFPE–frozen tumour pairs [V];
  - sarcoma PDX sc vs sn (BMC Cancer 2023) [V];
  - paired lung sc/sn (PLoS Genetics 2024) [V];
  - O'Flanagan 2019 Genome Biology: 6C vs 37C dissociation on patient tumours, PDX and lines; 155,165 cells; EGA EGAS00001003753 (controlled access) [V].

TIER C: LABELLED ARTEFACTS
- AUROC ≥ 0.85 for flagging the 11 hand-removed 3CA v2 MPs, Gavish's removed MPs and the protocol-named MPs of the 2026 catalogue.
- ≥ 80% of canonical MPs (cell cycle, hypoxia, IFN/MHC-II, pEMT) should be typed T1, T3 or T4.

TIER D: NEW EXPERIMENT, the decisive ground truth
- 4sU added to the dissociation buffer, following the Molecular Systems Biology 2023 19:e11147 method [V], so transcripts made during dissociation are T>C-labelled.
- n = 8–10 fresh tumours (breast, CRC, ccRCC/lung). Each specimen is split into warm 37C and cold protease arms, with Flex on adjacent FFPE and Xenium 5K on a serial section.
- Pre-registered predictions:
  - T5 ex vivo MPs have a labelled fraction ≥ 3× the T1 median;
  - Spearman(predicted ex vivo burden, measured labelled fraction) ≥ 0.6;
  - the TME nascent immediate-early index correlates with the labelled fraction at r ≥ 0.7.

TIER E: SPATIAL
- On Xenium 5K / CosMx 6K, T1 and T4 MPs show Moran's I above permutation nulls, and T5 MPs are not reconstructable.
- Pre-registered niches (hypoxia near necrosis, IFN near T cells), with a clone-domain control via CalicoST.

TIER F: GENOTYPE
- DNA- or allele-defined clones (clonealign data, CalicoST, FFPE subclonal CNAs from Yan 2026).
- T2-vs-T1 calls should agree between expression-derived and DNA-derived clones with Cohen's kappa ≥ 0.7.
- Calibrate rho_gen.

TIER G: GENERALITY AND CAUSALITY
- Kinker CCLE (198 lines) and Tahoe-100M DMSO: T1 recurs in vitro and T4 in vivo-only programs are absent. CCLE confounds a technical with a biological environment, so it is used only as a check.
- Tahoe drug responses with matched mechanisms vs DMSO: CDK4/6 inhibition → cell-cycle MP down; MEK inhibition → MAPK-linked MP down.
- TCCA replication, excluding studies shared with 3CA: type kappa ≥ 0.7; T1 replication ≥ 80%.
- Leave-one-study-out type stability.
- TCGA associations of environment-standardised scores against random-signature nulls (Venet 2011).

BASELINES
Gavish/3CA code, GeneNMF, cNMF and Harmony-cNMF (TCAT-style), GBCD (on subsets), ITHmapper (scVI + Hotspot), LIGER iNMF, scINSIGHT, MOFA-FLEX (groups = environment), inVAE, contrastiveVI, DRVI, consensus scHPF / scDEF, and Grabski multi-study NMF (environment as study).

ABLATIONS
- no kappa transport;
- no environment library;
- no TME sensors / anchors;
- observed-only environments (no latent imputation);
- no cis-dosage offset;
- score-only tests (no e-values);
- no power posterior.

### Flagship biology

"A protocol-invariant pan-cancer map of malignant metaprograms." Data: 3CA v2 (124 datasets, 2,836 samples), TCCA for replication, all matched twin sets, and the new 4sU tumour dataset. The expected outcomes below are HYPOTHESES to test, not results.

(1) Census.
- For each canonical 3CA MP, report its type (T1–T6) with effect sizes.
- Hypothesis: a minority of 3CA MPs, perhaps 10–25%, are environment-owned, protocol-confounded or clone-linked. Pilot 1 sizes this.

(2) Resolving "stress", the most contested MP.
- Prediction: the 3CA Stress MP and the 2026 catalogue's "Stress (in vitro)" decompose into two parts.
  - (a) A warm-dissociation-owned immediate-early component (FOS/JUN/EGR1/ATF3). It should have a high 4sU-labelled fraction, be absent or strongly reduced in Flex, Xenium and cold-protease data, and be predicted by the TME nascent immediate-early index.
  - (b) An invariant in vivo proteotoxic/hypoxia-adjacent component (HSP/DNAJ). It should be spatially organised near necrosis in Xenium 5K and carry prognostic signal in TCGA after environment standardisation.
- Gene-list filtering currently removes both parts. Separating them is a biological result in its own right.

(3) Protocol-named MPs.
- Test whether "single-nucleus" MPs are explained by kappa transport (nuclear retention, gene-length bias) plus nucleus-owned components, or whether they reveal real nuclear-retention biology.

(4) Plastic vs heritable states.
- Using allele-aware and DNA clones, estimate per MP the share of usage explained by genotype.
- Flag positional, arm-gain-driven programs (T2p).
- Identify T1 programs that are spatially mixed across clones: candidate reversible, therapy-relevant states.

(5) Unmasked programs.
- Test whether removing variance owned by environments and ambient RNA reveals shared programs that were previously merged with stress or EMT/ambient-collagen axes.
- Validate any new MP with Tier B and Tier E.

(6) Practical payoff.
- Environment-standardised MP scores let Flex/FFPE archival cohorts and fresh atlases be analysed in one reference frame, which matters as clinical single-cell work moves to FFPE Flex and 5K in situ panels.

### Software architecture

LANGUAGE AND ECOSYSTEM
- Python ≥ 3.10 on PyTorch + Pyro (GPU; CPU fallback).
- scverse-native: AnnData/MuData; SpatialData for Xenium, CosMx and Visium HD; scanpy-compatible.
- R/Seurat wrapper `InvarMPr` via reticulate, with SingleCellExperiment and Seurat converters.
- Distributed via pip and conda. Docker image provided.
- Snakemake/Nextflow pipeline for atlas runs.

MODULES
- invarmp.io: loaders for 3CA/TCCA/GEO matrices, raw empty-droplet matrices and Xenium/CosMx outputs.
- invarmp.env:
  - EnvSchema: an ontology of attributes and levels;
  - harmonize(): maps metadata to attribute tuples;
  - fingerprint(): computes TME sensors (nascent immediate-early index from velocyto/STARsolo spliced and unspliced layers, intronic fraction, mitochondrial fraction, ambient similarity, length bias);
  - impute(): posterior over missing attributes;
  - design_graph(): the testability diagnostic.
- invarmp.nuisance: ambient profiles (raw or CellBender outputs); spillover for in situ data; ingestion of CNA and clone calls (Numbat, CopyKAT, SCEVAN, CalicoST, clonealign); doublet priors (scDblFinder).
- invarmp.transport: fits kappa_e from calibrator cells; panel masks; panel Fisher information.
- invarmp.model.EnvNMF: the generative model and SVI fitting. Count splitting and data thinning live in invarmp.split.
- invarmp.test:
  - score_screen();
  - evalue_thinned_lr();
  - anchor_ci_test();
  - clone_test();
  - positional_test();
  - ebh().
- invarmp.type: assign() and profile cards (HTML/JSON).
- invarmp.reference: export/import of a versioned reference (Zenodo).
- invarmp.project: query mapping, standardisation and novelty detection.
- invarmp.sim: the environment simulator.
- invarmp.bench: metrics and baseline wrappers.

KEY API
```
im.env.harmonize(mdata, meta)
im.env.fingerprint(mdata, calibrators=['T','Fibro','Endo','Myeloid'])
im.env.impute(mdata)
model = im.EnvNMF(n_shared_max=200, env_per_level=4, private_per_tumor=3, likelihood='nb', power='auto', anchor_lambda='cv')
model.fit(mdata, malignant_key='malignant', tumor_key='tumor', env_key='env', clone_key='clone', cn_key='cn', reference_env='flex_ffpe')
res = im.test.invariance(model, axes='all', null='thinning', fdr=0.1)
cat = im.type.assign(res)
cat.export('invarmp-pan-v1')
q = im.project(query, reference='invarmp-pan-v1', env='xenium_5k', standardize_to='flex_ffpe')
```

OUTPUT LOCATIONS IN ANNDATA
- .varm['invarmp_W']
- .obsm['invarmp_usage'], ['invarmp_usage_sd'], ['invarmp_usage_std'], ['invarmp_exvivo_burden']
- .uns['invarmp']['profiles'] and ['types']

REPRODUCIBILITY
- A separate figure-reproduction repository.
- The benchmark is extensible: new methods plug in through a common wrapper interface.

### Feasibility

TEAM
- 1 senior computational postdoc (model and inference, 100%).
- 1 PhD student or research assistant (data harmonisation, simulator, benchmark, 100%).
- 1 statistician (identifiability proposition and e-value validity under NB thinning, 20–30%).
- 1 wet-lab scientist plus a pathology/surgical collaborator for the 4sU dissociation, Flex and Xenium 5K experiment (about 6 months part-time).
- PI.

COMPUTE
- 2–4 A100/H100 80 GB GPUs.
- About 5,000–10,000 GPU-hours in total, including seeds, grids, thinning refits and baselines. GBCD and Tetris-like baselines run on subsets.
- 3–5 TB storage.

DATA
- Public: 3CA v2; TCCA; Kinker CCLE; Tahoe-100M; GSE243280 (Janesick); GSE280311 (Oliveira Flex) plus the Visium HD/Xenium releases; the deposits of the Nature Communications 2025 and Nature Methods 2026 platform benchmarks; Trinks lung; the ccRCC 2025 preprint data.
- Slyper 2020: repository [P].
- O'Flanagan 2019: EGA, controlled access. Budget 1–3 months for the data access committee.
- New: 8–10 tumours × (warm 4sU, cold 4sU, Flex, Xenium 5K). Rough cost USD 80–150k. This is an estimate, to be quoted.
- Spliced/unspliced layers require re-alignment from FASTQ/BAM where only matrices are public. Plan for STARsolo re-processing of about 20–40 key studies.

TIMELINE
- M0–3: pilots, metadata harmonisation, design-graph census.
- M3–9: EnvNMF v1, simulator, theory.
- M6–12: full atlas fit and tests.
- M8–14: wet lab (IRB, collection, sequencing).
- M12–18: validation tiers B–G and baselines.
- M18–22: reference release, package v1.0, writing.

Realistic total: 20–24 months. A computational-only preprint is possible at about 12–14 months.

### Risks and mitigations

1. Matched data are scarce. Janesick is one block, and Oliveira, Trinks and the ccRCC study have few patients.
- Mitigations:
  - Scale through TME calibrators (kappa) and TME-sensor anchors (unmatched invariance tests).
  - Use natural experiments inside 3CA (sc/sn, Smart-seq2/10x, 3'/5', sites).
  - Generate the new 4sU/Flex/Xenium dataset.
  - Report T6 honestly where bridging fails.

2. Protocol is confounded with biology. For example, single-nucleus data are mostly brain or frozen tissue, and Flex data come from FFPE archives of particular cancer types.
- Mitigations: the design-graph identifiability check, tests stratified by B_t, and the simulator's confounding grid, which shows the breakdown point.

3. Assumption A1 (TME identity programs are invariant) may fail. For example, nuclear-retention bias may differ between cell types.
- Mitigations:
  - Test A1 on twins.
  - Allow cell-type-specific kappa deviations with a shrinkage prior.
  - Use gene-feature-based kappa so that transport generalises.

4. Real biology can depend on environment through selection. Fragile mesenchymal-like cells are lost in fresh dissociation.
- Mitigations:
  - The hurdle split separates selection (Δpi) from induction (Δm). Such programs are typed T3 (biology), not T5.
  - Spatial data adjudicate.

5. In vivo stress overlaps ex vivo stress, and the TME may also respond to in vivo stress (e.g., hypoxia).
- Mitigations:
  - Anchors are chosen to be technical (nascent immediate-early fraction, intronic fraction, ambient similarity), not biological.
  - The nascent-fraction hypothesis is validated against 4sU.
  - Loadings are allowed to split into an environment-owned part and an invariant part.

6. Thinning and e-value validity under NB with estimated overdispersion (Hivert et al.).
- Mitigations:
  - Exact Poisson restricted problems.
  - phi estimated on an independent fold.
  - Calibration on null twins and in simulation.
  - Conservative averaging of e-values.

7. Panel restriction (313- or 422-plex) gives low power.
- Mitigations: a panel-informativeness threshold before testing, with priority on 5K/6K panels and Flex/Visium HD.

8. Circular CNA inference.
- Mitigations: allele-aware calls, a DNA-validated subset, and a cis offset (only trans states are typed).

9. Identifiability of shared vs owned components in noisy data.
- Mitigations: a formal proposition, seed-consensus dispersion and a simulation breakdown map.

10. Reviewer objections: "IRM fails in practice" and "this is batch correction".
- Responses:
  - InvarMP does not optimise IRM. The IRM gradient is used only as a score screen. Decisions rest on explicit Bayesian ownership plus calibrated e-values.
  - Unlike batch correction, environment-linked variation is modelled, tested and returned, not removed.

11. Scoop risk: ITHmapper, a TCAT-group extension to malignant cells, a GBCD extension, or inVAE-like work.
- Mitigations: run the mandatory bioRxiv/Scholar check from landscape section 6.3, with added queries "invariant gene programs", "protocol-invariant", "Flex single-nucleus metaprogram". Lead with the invariance definition and the 4sU ground truth.

12. Not enough tumour specimens for the experiment.
- Mitigation: substitute PDX tumours, splitting fresh vs FFPE vs in situ, and keep the cell-line 4sU arm as a fallback.

### Go/no-go pilots

The pilots take 6–8 weeks in total and use public data only.

P1. Invariance audit of the current standard (2–3 weeks).
- Score the 3CA v2 MP gene lists, and per-tumour NMF programs, in matched twins:
  - Trinks lung (fresh / cryopreserved / FFPE Flex);
  - ccRCC (fresh sc / sn-Multiome / snFlex);
  - Slyper sc vs sn;
  - Oliveira Flex vs Xenium;
  - Janesick Flex vs Xenium.
- Use a crude transport: kappa from pseudobulk of matched TME cell types.
- Statistic: transported usage ratio and loading cosine against a null built by thinning one environment into two pseudo-environments.
- GO if ≥ 10% of testable MPs (and at least 3) shift beyond the 99th percentile of the thinning null while ≥ 50% stay invariant. That shows the axis discriminates.
- NO-GO or REFRAME:
  - if every MP is invariant, artefacts are not the problem; pivot to calibration and genotype typing;
  - if every MP shifts, the transport is inadequate; fix kappa before going further.

P2. TME-sensor anchor (2 weeks).
- Across 3CA samples with TME cells, regress per-sample malignant Stress-MP usage on the immediate-early level of co-processed T cells and fibroblasts, controlling for cancer type and study.
- GO if partial R² ≥ 0.3 for Stress and ≤ 0.05 for Hypoxia and Cell-cycle MPs. This is the specificity check.

P3. Nascent-fraction fingerprint (2 weeks).
- In data with spliced/unspliced layers and known warm vs cold dissociation or sc vs sn (O'Flanagan, if EGA access arrives in time; otherwise sarcoma PDX sc vs sn, or Slyper), test whether the immediate-early unspliced fraction in TME cells separates warm-dissociated samples.
- GO if AUROC ≥ 0.8. If not, drop sensor (i) and rely on the metadata, intronic, mitochondrial and ambient sensors.

P4. Design-graph census (1 week).
- Build the environment × tumour/cancer-type graph for 3CA v2 + TCCA + twins.
- GO if ≥ 30% of candidate MPs (from a quick GeneNMF pass) are used by tumours that bridge at least 2 technical levels on at least one axis.
- Otherwise the paper must rest on anchors plus the new experiment; reassess scope.

P5. Prototype simulation (3 weeks, can overlap).
- EnvNMF v0 on 200 tumours with planted warm-owned and nucleus-owned artefacts, Cramér's V = 0.5, 5% twins.
- GO if artefact AUROC ≥ 0.9 and invariant-program recovery cosine ≥ 0.85, while Gavish code and GeneNMF call ≥ 1 planted artefact an MP.

Decision rule: proceed to full build if P1, P4 and P5 pass, and at least one of P2 or P3 passes.

### Why Nature Methods

1. Conceptual advance.
- The paper turns the field's vague complaint ("MPs are not robust, sometimes not real") into an operational, testable definition: biology is what is invariant under technical interventions.
- That definition has a model, an identifiability condition (design-graph bridging), calibrated error control (thinned split-likelihood-ratio e-values with e-BH) and an interpretable output (invariance profile and type).
- Generic factorisation or batch-correction tools cannot reproduce this. They either remove environment-linked variation or ignore it.

2. Timeliness.
- Clinical and atlas single-cell work is moving fast to FFPE Flex, single-nucleus and 5K-plex in situ panels.
- Every lab now has to decide whether a program found in fresh 3' data exists in its fixed or in situ cohort. InvarMP answers that, and provides the transport and projector.

3. Unusually strong ground truth for a programs paper.
- Matched multi-platform tumours (seven or more public sets).
- A new 4sU-labelled tumour dissociation dataset that marks, for each transcript, whether it was made during dissociation (to our knowledge the first such tumour dataset; not yet checked).
- Spatial, DNA-clone and perturbation (Tahoe) checks.
- The labelled 3CA artefacts.

4. A biological payoff that matters to cancer biologists: dissecting ex vivo from in vivo stress, plastic from clone-linked states, and protocol-named MPs, across a 2,836-sample atlas with TCCA replication.

5. Generality beyond tumours. The same machinery applies to any atlas mixing protocols (Human Cell Atlas-scale sc/sn/Flex/in situ integration).

6. Deliverables at TCAT-level polish, positioned as enablers rather than the claim: a scverse package, a versioned malignant-MP reference with uncertainty, a projector into Flex/Xenium/CosMx, and an open simulator and benchmark.

7. Fit with the landscape. The work avoids the saturated novelty claims (faster Gavish, foundation-model priors, reference+projection, shared vs private). It lands squarely in the verified white space: confounder modelling inside discovery, plus error control.

CITATION STATUS
- Verified this session by web search: Slyper 2020; O'Flanagan 2019; van den Brink 2017; Denisenko 2020; MSB 2023 4sU method; FixNCut 2024; Trinks lung; ccRCC preprint; Yan 2026; Nature Communications 2025 and Nature Methods 2026 platform benchmarks; Janesick and Oliveira details; inVAE; scPFA; ICP; anchor regression.
- [U], from memory: IRM (Arjovsky et al.), REx, universal inference (Wasserman et al.), e-BH (Wang & Ramdas), e-value merging (Vovk & Wang). Check before citing.

### Differentiation from closest competitors

| Competitor | Differentiator |
|---|---|
| GBCD (Liu…Stephens, Nature Genetics 2025) | GBCD splits shared vs patient-specific programs with a Gaussian model on log data, per cancer type. It has no notion of technical environment, so an artefact that recurs across patients (warm-dissociation immediate-early genes, ambient) is by construction a 'shared' program. InvarMP adds four things GBCD lacks: (1) technical environments as interventions, with environment-owned libraries and invariance tests; (2) kappa transport across cell/nucleus/Flex/panels; (3) genotype typing; (4) a count likelihood with credible tumour membership at atlas scale. InvarMP does not pitch shared vs private as its novelty. |
| ITHmapper / Pan-Cancer ITH compendium (Abelson lab, bioRxiv 2026) | ITHmapper treats sample as batch in scVI, which removes environment-linked and whole-tumour variation indiscriminately, and then outputs discrete Hotspot states. InvarMP keeps environment as information. It tests which overlapping continuous programs are invariant, keeps artefact programs as typed outputs, and preserves whole-tumour usage. Benchmarked head-to-head on held-out twin prediction. |
| TCAT/starCAT (Kotliar…Raychaudhuri, Nature Methods 2025) | TCAT defines reproducibility as recurrence across Harmony-corrected datasets. Artefacts shared by all fresh-dissociated 3' datasets pass that test, and the authors warn about transfer to Flex. InvarMP requires invariance across protocols, explicitly learns the Flex/nucleus/panel transport, and projects with uncertainty for each panel. It covers malignant cells rather than T cells. The reference/projector is a deliverable, not the claim. |
| Grabski–Trippa–Parmigiani Bayesian multi-study NMF (Genome Biology 2025) and the Tetris/APAFA line | In these models, studies are exchangeable groups with inclusion indicators, and technical vs biological groupings are not distinguished. InvarMP borrows subset-sharing, but its owners are environment ATTRIBUTES (warm, nucleus, 5'). It adds program-specific deviations, hurdle selection/induction terms, TME calibrators, a design-graph identifiability condition and calibrated invariance e-values, all for scRNA-seq nested as cells in tumours in environments. |
| inVAE (Aliee…Teichmann, Vento-Tormo, Theis; bioRxiv 2024) | inVAE is the closest 'invariant vs spurious' idea in single-cell work. It is a nonlinear VAE latent split that needs labelled covariates and produces representations, not non-negative gene programs. It has no multi-tumour sharing, no tests or types, no matched-platform transport, no clone axis, and no handling of partially observed environments. InvarMP's estimands are program loadings and usage with e-values. |
| contrastiveVI (Nature Methods 2023) / bcNMF (arXiv 2026) | These need one designated background dataset. InvarMP has no clean background: every environment owns its own artefacts, and biology is the intersection ('no clean room' principle). Ownership is learned across many environments that are only partly labelled. |
| LIGER iNMF / scINSIGHT | Their dataset-specific factors cannot tell technical from biological dataset specificity, use no count model, and require labels and hand-chosen K. InvarMP's ownership follows environment attributes, uses a design graph and calibrated tests, and handles missing attributes. |
| Gavish 2023 / 3CA v2 / GeneNMF / scalop-family | These are threshold stacks that define robust as recurrent and end in manual curation. InvarMP replaces recurrence-as-sufficient with invariance and genotype typing under error control. It uses their outputs, including the 11 hand-removed MPs, as labelled benchmarks. |
| CanSig benchmark (Cancer Research 2026) | CanSig benchmarks state discovery against known signatures. InvarMP's benchmark adds planted environment-owned artefacts, protocol–biology confounding, panels, missing environment labels, nuisance-only null cohorts and held-out twin prediction. |
| scDEF, consensus scHPF, sciRED, Spectra | These are count or covariate models with batch terms or priors, but none treats protocol as an intervention, tests invariance, or types programs. sciRED regresses known covariates, which also removes biology that is correlated with protocol. |
| scPFA lineage-informed factor analysis (Staklinski & Siepel, 2026, per ResearchGate listing) | scPFA estimates heritable programs from lineage trees. InvarMP uses genotype only as a typing axis inside multi-tumour, multi-environment discovery. The two are complementary, and scPFA could serve as an external heritability check on lineage-traced data. |

### Referee 1: novelty and prior art

**Predicted editorial outcome:** My prediction is "send to review", then reject or major revision. The odds of eventual acceptance at Nature Methods are low, perhaps 15–25%, unless the work is cut down substantially. Three things favour review: the topic is timely (Flex, single-nucleus and 5K-plex panels in clinical work), the experimental ground truth is unusually strong for a programs paper (warm/cold dissociation with 4sU, Flex and Xenium 5K on the same specimens), and the 3CA re-analysis has an obvious audience. Desk rejection becomes likely in two cases. First, if the submitted atlas results show mostly T6 ("untestable") or mostly type labels decided by anchors, because technical environment is confounded with study and cancer type across 3CA. Second, if the only MPs typed as artefacts are stress, low-quality and contamination MPs that 3CA already removes by hand or that are already known. The paper would then read as a Bayesian NMF that accounts for batch, plus a validation suite. I expect reviewers to object on four points: (i) the invariance framing is already present in single-cell work (inVAE) and in the field's own validation practice, which is cross-protocol and in situ checking; (ii) the measurement transport κ, pinned by microenvironment (TME) "calibrator" cells, can absorb the very dissociation artefact it is meant to detect; (iii) the e-values are not valid as specified and will have negligible power under e-BH with so few matched tumours; (iv) the model is over-built, at roughly ten interacting component libraries fitted by mean-field stochastic variational inference (SVI), so its credible sets will not be calibrated. A cut-down version could be competitive: environment-owned libraries, transport identified from data with fixed and in situ measurements, held-out twin prediction, and the 4sU experiment. The full concept as written is not.

**Fatal flaws**

- THE TRANSPORT CAN ABSORB THE ARTEFACT IT IS MEANT TO DETECT. Assumption A1 says TME identity programs do not change with technical environment, and that a gene-wise operator κ_{g,e}, shared by all programs and all cell types, is pinned by TME calibrator cells. Yet section 5 uses those same TME cells as SENSORS of warm dissociation, through their nascent immediate-early index, which assumes they respond to the protocol. Both cannot hold. Collagenase at 37C induces FOS/JUN/EGR1/HSP genes in every compartment: O'Flanagan 2019 [V, verified this session] defines a 512-gene core stress set across cell types, and MSB 2023 19:e11147 [V] reports a dissociation response with a core shared across tissues plus cell-type-specific programs. Because log κ carries a free gene-wise residual ε_{g,e}, κ_{FOS,warm} will rise to fit the calibrators. The ex vivo program in malignant cells will then look loading-invariant and pass D2, i.e. a false T1 or T3 call. That is exactly the failure the method exists to prevent. The fix requires κ to be a low-dimensional function of gene features only, or κ to be identified from fixed or in situ data, plus a simulation showing that a planted warm-dissociation program shared across compartments is not absorbed into κ.
- THE ATLAS-SCALE CLAIM IS NOT TESTABLE ON EXISTING DATA. The one-sentence claim says invariance typing 'replaces manual curation at pan-cancer scale'. The public matched designs are tiny: Janesick is a single breast block, Oliveira has about 5–8 CRC patients, and Trinks and the ccRCC preprint have few patients each. Across 3CA, technical environment is largely nested in study and cancer type (single-nucleus mostly brain or frozen; Smart-seq2 in particular studies; Flex in particular FFPE archives). Under the proposal's own design-graph proposition, most MP × axis pairs will be T6, except where they are typed through the unmatched route: TME-sensor anchors and the test r̄ ⟂ A | B. That route depends on the A1-contradicting sensors above, and on B_t fully capturing the biology that is correlated with protocol, which it never does in observational atlases. The axes that matter most for the 'is it real?' question (warm vs cold, dissociated vs in situ) are precisely the ones with almost no bridging tumours. Pilot P4 (≥30% of candidate MPs bridging on ≥1 axis) can pass on the sc/sn and 3'/5' axes alone while the stress-relevant axes stay untestable. The headline must then shrink to 'typing where a bridge exists'.
- THE ERROR-CONTROL HEADLINE IS NOT VALID AS SPECIFIED. (a) Universal inference needs the denominator to be the likelihood MAXIMISED OVER THE WHOLE COMPOSITE NULL, or at least at least as large as the likelihood at the true null parameter. In the proposal, the restricted model on X2 fixes w_k, κ, φ, every other component and the latent environments at their X1 estimates. Those estimates are not the true nuisance values. Even when δ = 0, L_X2(restricted) can fall below L_X2(θ_true), so E[E] ≤ 1 is not guaranteed and the test rejects because of nuisance-estimation error. The 'convex Poisson problem in usage' guarantees a global optimum only of the wrong null. (b) Negative-binomial thinning leaves the folds dependent when φ is mis-estimated, and cell-level gamma or hurdle latents are shared between folds (Hivert et al., noted in the landscape). (c) Even if valid, power is the problem. With e-BH at q = 0.1 over m ≈ 1,600 (k, axis) pairs, the k-th rejection needs e ≥ m/(qk): about 320 for 50 rejections and about 80 for 200. A split likelihood-ratio statistic on half the counts from a handful of matched tumours will reach that only for gross shifts, such as immediate-early genes, which are already known. Averaging e-values across evidence sources, done 'for validity', dilutes them further. The likely outcome is either anti-conservative calls or almost no calls.

**Major concerns**

- NOVELTY IS MOSTLY A REFRAMING. Cross-protocol invariance (cold vs warm dissociation, sc vs sn, Flex or in situ vs fresh) is already how the field checks whether a program is real. The landscape itself lists it as validation metric C (section 4.2/4.3). Moving that check from validation into the definition is a reframing, not a new capability. Every model component exists: dataset- or condition-specific factors (LIGER iNMF, scINSIGHT); multi-study count NMF with inclusion indicators and covariates (Grabski 2025); IBP sharing (Tetris); private and shared partitions (APAFA); power posteriors (BayesPowerNMF); thinning (Neufeld); universal inference and e-BH; anchor regression and invariant causal prediction (ICP); ambient and doublet models (CellBender, scDblFinder); copy number from Numbat. Invariance-based disentangling in single-cell data is also not new: inVAE (Aliee…Theis, bioRxiv 2024 [V]) splits invariant (biological) from spurious (technical) latents, and ICP- or IRM-style environment-stable features have appeared in single-cell work (CausalCellInfer, medRxiv 2024 [V, snippet only]; domain-invariant patient classification, bioRxiv 2025.09.22.677881 [V, snippet only]). What is new is the combination: non-negative programs, ownership by environment attribute, a selection/induction split and testability reporting. That is a legitimate contribution, but on its own it is closer to Genome Biology or Nature Communications than to a Nature Methods conceptual advance.
- THE STEPHENS LAB IS CLOSING THE GAP THE PROPOSAL USES AGAINST GBCD. The GBCD differentiator rests on 'count likelihood + covariates'. That lab now has covariate-moderated empirical-Bayes matrix factorisation (Denault, Tayeb, Carbonetto, Willwerscheid, Stephens; NeurIPS 2025 / arXiv 2505.11639 [V]) and Poisson NMF with a shifted-log link (Weine, Carbonetto, Irizarry, Stephens; arXiv 2601.05845 [V]). A count-based GBCD that takes protocol as a covariate is a plausible near-term release. The proposal should benchmark against GBCD with protocol or environment covariates, not only plain GBCD.
- 'INVARIANT UNDER MEASUREMENT CHANGE' DOES NOT MEAN 'BIOLOGY'. Every environment shares the same pre-analytical history: surgical warm ischaemia, cold ischaemia before fixation, biopsy handling and transport delay. An ischaemia- or handling-induced program, such as immediate-early, hypoxia or HSP genes induced during time to fixation, will be loading- and usage-invariant across fresh, Flex, single-nucleus and Xenium samples from the same specimen, and will be typed T1. Choosing Flex-FFPE as the reference e0 bakes this in, and it contradicts the proposal's own 'no clean room' principle. Pre-analytical variables need their own axis, with ischaemia time recorded as an anchor, or the claim must be narrowed to 'invariant to dissociation, molecule and readout'.
- THRESHOLDS RETURN AFTER BEING CRITICISED. A horseshoe never yields δ exactly 0. Type assignment uses transported cosine < 0.5 for T5, ρ_gen ≥ 0.3 with a clone e-value ≥ 20 for T2, and 'testability' flags. The six-way typing is a decision tree with tunable cut-offs, which is the same criticism the proposal levels at 35/10/10. T3 vs T5 is especially fragile. A program present at low usage in vivo but amplified at 37C (T3-induction) and a purely ex vivo program (T5) differ only in whether usage in the fixed or in situ environment is small or zero. Flex's lower sensitivity and the coverage of 313/422/5K panels make those two cases nearly indistinguishable.
- THE MODEL IS OVER-PARAMETERISED AND ITS POSTERIORS WILL NOT BE CALIBRATED. The likelihood combines a shared ARD library (200), environment-owned libraries (4 per attribute level), 3 private components × about 3,000 tumours, a cell-type library, ambient, doublet priors, gene-wise cis-dosage exponents, smooth κ functions with gene residuals, horseshoe δ for every program × environment × gene (200 × ~20 × 12k ≈ 5×10^7 parameters), hurdle usage with horseshoe η, finite-IBP membership, latent environments, biological-stratum δ^b, a power posterior and an anchor penalty. Separating these libraries from one another (shared vs environment-owned vs private vs cell-type vs ambient) is weakly identified in noisy data. The proposed identifiability theorem extends noiseless Fu/Huang/Sidiropoulos results and does not cover the full model. Mean-field SVI with relaxed Bernoullis systematically underestimates posterior uncertainty, so the targeted 0.85–0.95 coverage of tumour-membership credible sets is unlikely without heavy recalibration. Five-seed Hungarian consensus is a heuristic layered on a model that claims to be Bayesian.
- THE ANCHOR PENALTY IS USED TWICE. λ_A pushes tumour-level usage away from the technical anchors while the model is fitted. The same anchors are then used in the conditional-independence test (c), on the same data. The penalty biases the test toward 'invariant'. It also suppresses real biology correlated with protocol in exactly the confounded strata, such as neural programs in brain tumours profiled by single-nucleus, before typing can flag them as T6. Calling the tumour-level test 'invariant causal prediction' oversells it: it is a conditional-independence regression with study random effects. Also, the IRMv1 penalty is the gradient with respect to a scalar dummy scale, whereas the gene-wise score S_kg is the Rao score for a vector deviation. They are related but not 'exactly' the same, so that wording should be fixed.
- THE TME 'NASCENT IEG INDEX' SENSOR IS SPECULATIVE AND TAKES A LOT OF ENGINEERING. Unspliced counts for about 50 immediate-early genes in the non-malignant cells of each sample are sparse and noisy. The 'minutes of dissociation are rich in unspliced reads' hypothesis is unvalidated. Using it requires re-aligning 20–40 studies from FASTQ, which are often controlled-access. If P3 fails, the unmatched scaling route loses its best technical anchor.
- THE FLAGSHIP BIOLOGY IS PARTLY ANSWERED ALREADY. Baron et al., Cell Systems 11:536 (2020) [V] validated a FOS/JUN/HSP stress-like cancer cell state with spatial transcriptomics in zebrafish melanoma and human samples, and linked it to tumour seeding and drug resistance. O'Flanagan 2019 [V] defined the collagenase core and showed that cold protease minimises it. FixNCut (Genome Biology 25:81, 2024 [V]) and Flex remove dissociation artefacts experimentally. MSB 2023 [V] already introduced 4sU labelling during dissociation (zebrafish larvae, mouse cardiomyocytes), with core plus cell-type-specific programs and computational removal of labelled transcripts. The proposed 'stress decomposition' is therefore a quantitative pan-cancer extension, not a discovery. The claim of the 'first tumour 4sU dissociation dataset' is unverified and in any case a straightforward application of an existing method. A reviewer will also ask: if 4sU, cold protease and Flex already resolve the dissociation axis experimentally, why is a ten-component model needed for that axis?
- THE CLONE TYPING DUPLICATES EXISTING POST-HOC ANALYSES. The landscape notes that Kinker found only 8% of continuous programs associated with CNA subclones. ρ_gen together with arm enrichment is close to analyses Gavish and Kinker already ran after the fact. Folding them into the profile card is useful but not new, and the cis-offset prior (gamma ~ LogNormal(log 0.8, 0.5)) is a design choice whose sensitivity must be shown.
- THE PLAN DOES NOT FIT THE TEAM. One postdoc, one student, a 20–30% statistician and part-time wet-lab support are proposed for: a new GPU Bayesian engine, a theory paper (identifiability plus validity of e-values under NB thinning), a semi-synthetic simulator with a 4×4×4×3×50 grid, more than 13 baselines (GBCD and Tetris at scale), STARsolo re-processing, an EGA data-access application, IRB approval, a new 4sU/Flex/Xenium experiment, and a versioned atlas release. All of this in 20–24 months is not credible. The computational-only preprint at 12–14 months would lack Tier D, the only decisive ground truth.

**Prior-art threats**

- inVAE: conditionally invariant representation learning (Aliee, Kapl, Pham, Cakir, Jimba, Cranley, Teichmann, Meyer, Vento-Tormo, Theis; bioRxiv 10.1101/2024.12.06.627196; theislab/inVAE). It separates invariant (biological) from spurious (technical) latents in single-cell data, which pre-empts the 'invariance' framing at the level of representations. VERIFIED this session (bioRxiv and GitHub listings).
- Baron, Tagore, Hunter, …, White, Yanai. 'The Stress-Like Cancer Cell State Is a Consistent Component of Tumorigenesis.' Cell Systems 11(5):536–546 (2020), doi:10.1016/j.cels.2020.08.018. Spatial transcriptomics plus human samples validate the stress state as real in situ, which partly pre-empts flagship result (2). VERIFIED (PubMed listing via search).
- 'A single-cell RNA labeling strategy for measuring stress response upon tissue dissociation.' Molecular Systems Biology 19:e11147 (2023). 4sU added during dissociation; zebrafish larvae and mouse cardiomyocytes; a core plus cell-type-specific dissociation response; computational removal of labelled transcripts. This undercuts A1 (TME cells respond to dissociation) and the novelty of Tier D. VERIFIED (title, venue and key findings via search snippets; full text blocked; authors not retrieved).
- O'Flanagan et al., Genome Biology 2019: cold active protease; 512-gene core collagenase stress set, shared across cell types. This is evidence against A1 on the warm/cold axis. VERIFIED.
- FixNCut: reversible tissue fixation before dissociation (Genome Biology 25:81, 2024). An experimental route that removes dissociation artefacts. VERIFIED.
- Covariate-moderated empirical-Bayes matrix factorisation (Denault, Tayeb, Carbonetto, Willwerscheid, Stephens; NeurIPS 2025 / arXiv 2505.11639) and 'A New Family of Poisson NMF Methods Using the Shifted Log Link' (Weine, Carbonetto, Irizarry, Stephens; arXiv 2601.05845). Together they give the GBCD group the pieces for a count-based GBCD that uses covariates, which erodes the 'count likelihood + environment' differentiator. VERIFIED (titles and authors via search; abstracts not read).
- sysVI (Hrovatin, Moinfar, … Theis; BMC Genomics 2025; in scvi-tools). Integration across 'systems', including single-cell vs single-nucleus protocols, with VampPrior and cycle consistency. A baseline the proposal omits for the sc/sn axis. VERIFIED.
- 'Differences in molecular sampling and data processing explain variation among single-cell and single-nucleus RNA-seq experiments', Genome Research 34:179 (2024). Gene-length bias between sc and sn, and length correction. Precedent for κ as a function of gene features. VERIFIED (title and venue; authors not retrieved). Related: Bisque (gene-specific transformations between technologies) [snippet only] and gimVI (arXiv 1905.02269; technology-specific gene-level shift between scRNA-seq and spatial data) [VERIFIED].
- Invariance-based or invariant-causal-prediction methods already in single-cell work: CausalCellInfer, an ICP-inspired procedure for environment-stable marker genes (medRxiv 2024.10.17.24315646), and 'Domain-Invariant Feature Learning for Patient-Level Phenotype Prediction from Single-Cell Data' using IRM/REx across studies (bioRxiv 2025.09.22.677881). Both were VERIFIED as existing via search snippets; methodological details UNVERIFIED.
- ITHmapper / 'A Pan-Cancer Single-Cell Compendium of Intratumoural Heterogeneity' (Nofech-Mozes, Cook, Ouellette, Hafezi-Bakhtiari, Awadalla, Abelson; bioRxiv 10.64898/2026.01.06.693992). More than 3.6M cells, more than 1,000 tumours, 14 cancer types, with spatial data and drug screens. VERIFIED. It pre-empts the positioning as a pan-cancer package with a mapper.
- TCAT/starCAT (Nature Methods 2025; s41592-025-02793-1). Reference catalogues now extend beyond T cells (MYELOID.GLIOMA.V1, BONEMARROW.CD34POS.HSPC.V1 in the starCAT repo). A malignant catalogue from the same group is plausible. VERIFIED.
- Grabski, Trippa, Parmigiani, Bayesian multi-study NMF (Genome Biology 2025) and the probit follow-up (arXiv 2502.01468). The statistical core of subset-shared count NMF with covariates, so treating 'environment as study' is a trivial port. VERIFIED.
- LIGER iNMF / scINSIGHT / MOFA+ groups / scDEF (still a bioRxiv preprint per search, v2 2022) / KINOMO / MultiRFM (study-shared plus study-specific factors for RNA-seq, arXiv 2506.18478). These cover partial shared vs specific splits. VERIFIED (from the landscape, plus a MultiRFM search).
- Matched multi-platform datasets that weaken the 'unique ground truth' angle (the data exist and others are mining them): 'Benchmarking of spatial transcriptomics platforms across six cancer types' (bioRxiv 2024.05.21.593407) [VERIFIED title]; the ccRCC fresh/sn-Multiome/snFlex comparison (bioRxiv 2025.11.25.690485) [VERIFIED]; Janesick 2023 and Oliveira 2025 [VERIFIED in the landscape].
- The search was limited: arXiv, nature.com, NCBI and the Wiley/EMBO full texts were proxy-blocked this session. No direct 'protocol-invariant metaprogram' paper was found, but absence is weak evidence.

**Required changes**

- Resolve the κ circularity before anything else. Either restrict log κ to a low-dimensional smooth function of gene features with no free gene-wise residual on the warm/cold axis, or identify κ from fixed and in situ references (Flex or Xenium on the same specimen), or from spike-ins or matched bulk RNA, rather than from TME cells that respond to dissociation. Add a simulation in which a warm-dissociation program is planted in BOTH malignant and calibrator cells, and show it is typed T5 rather than absorbed into κ. Test A1 directly on twin data for every axis, and report it.
- Replace or prove the e-value construction. Either maximise the X2 restricted likelihood over all nuisance parameters (κ, φ, other components, usages), or switch to a test whose validity does not depend on plugging in nuisance estimates, for example conditional permutation of environment labels within matched tumours, or thinning-based 'null twins' as the PRIMARY calibration with empirical FDR. Report a power analysis for e-BH at m ≈ 1,600 against realistic numbers of matched tumours, and show the empirical FDR on null cohorts with mis-estimated φ.
- Publish the design-graph census up front, before claiming 'pan-cancer replacement of curation': for each 3CA MP and each axis, testable vs T6. If most stress-relevant axes are untestable, narrow the claim to bridged axes and to the new experiment.
- Add pre-analytical variables, at minimum ischaemia or time-to-fixation where recorded, as an environment axis or anchor. Justify why Flex-FFPE is the reference e0, or give the results for the alternative reference choices.
- Cut the model down to a core whose contribution can be shown in ablations: shared plus environment-owned non-negative libraries, the transport, the hurdle selection/induction split, and the cis offset. Move the latent-environment imputation, the anchor penalty inside the fit, the IBP and the biological-stratum δ^b to optional modules, or drop them. Show identifiability empirically with a breakdown map, and show calibration of credible sets under the actual SVI. If coverage fails, report frequentist calibration via thinning or the bootstrap instead.
- Do not use the anchors both inside the fit and in the test. Use λ_A = 0 for any fit whose output is tested against those anchors, or cross-fit by study.
- Add simple, strong baselines that a reviewer will demand. (a) Score 3CA MP gene sets in matched twins with a crude transport (pilot P1 itself), and show what the full model adds beyond it. (b) Run discovery only on fixed or in situ data (Gavish or GeneNMF on Flex and Xenium 5K alone) as a 'dissociation-free MP' baseline. (c) Regress out the O'Flanagan core, or use 4sU-labelled reads, before standard MP discovery. (d) inVAE, sysVI-based programs, covariate-moderated EBMF or GBCD with environment covariates, Grabski with environment as study, and LIGER with environment-defined datasets.
- Pre-register Tier B (held-out twin prediction) and Tier D (4sU) as the primary endpoints, and report them even if they are negative. Cite MSB 2023 and Baron 2020. Remove the 'first tumour 4sU dataset' claim unless it is checked against the literature.
- Reframe the flagship biology. The stress decomposition must go beyond Baron 2020 and O'Flanagan 2019, for example quantifying per tumour the ex vivo share of each canonical MP, or finding a canonical MP other than stress or low-quality that is not invariant. Otherwise the biological payoff is incremental.
- State, in place of the 0.5 cosine, ρ_gen and e-value cut-offs, decision rules justified by loss or calibration, and show that type labels are stable to these choices (sensitivity grid).
- Correct the technical over-claims: IRMv1 vs the vector Rao score; ICP vs conditional-independence regression; 'convex' holds only for the Poisson restricted problem in usage, not for NB or hurdle models.
- Re-plan the scope and team (at least two methods postdocs and a dedicated statistician), or split the work into a methods paper (core model, simulator, twin prediction, 4sU) and a later resource paper (pan-cancer catalogue and projector).

**Ideas worth keeping**

- Ownership by environment ATTRIBUTE, not by dataset: nuisance components attached to protocol levels (warm dissociation, nucleus, 5', site) and usable only by cells with that level. This gives a reusable, typed artefact library that replaces gene blocklists. It is the most distinctive modelling idea here.
- Splitting protocol effects into selection (Δπ, fraction of active cells) and induction (Δm, active intensity). This interpretable distinction keeps protocol-modulated biology instead of discarding it.
- The design-graph testability analysis, with an explicit 'untestable/confounded' (T6) output. Turning protocol–biology confounding into a reported quantity is honest and new for MP work.
- Held-out twin prediction ('invariance predicts the other platform') as a falsifiable benchmark. It is stronger than recurrence-based reproducibility and could stand on its own as a benchmark contribution.
- The paired experiment, split into warm 37C and cold protease arms with 4sU during dissociation, plus Flex and Xenium 5K on adjacent material. Even if not a first, it is a valuable public ground-truth resource for the whole field.
- Using the 11 MPs hand-removed from 3CA v2, and the protocol-named MPs of the 2026 catalogue, as labelled artefacts for automatic QC (AUROC).
- Panel-aware projection, with Fisher-information-based uncertainty for Xenium and CosMx panels, and treating genes outside a panel as missing rather than zero.
- Thinning-generated 'null twins' (one environment split into two pseudo-environments) as an empirical calibration device for any invariance statistic.
- A cis copy-number dosage offset inside discovery, combined with genotype typing using DNA- or allele-defined clones. This labels clone-linked programs rather than deleting them.

### Referee 2: statistical rigor and feasibility

**Predicted editorial outcome:** SEND TO REVIEW, then most likely MAJOR REVISION BORDERING ON REJECT if the paper looks like this concept. Why it would be sent out: the claim is timely and squarely Nature Methods material. It turns 'MPs are not real' into a falsifiable definition: biology is what stays fixed across technical interventions on the same tumour. The deliverables (transport and projection into Flex/Xenium, a design-graph audit, an artefact library) are what the field needs. Why it would not be accepted as is: (1) The confirmatory machinery does not deliver what it advertises. The 'invariant' label (T1) is assigned when a test fails to reject. The split-LR e-value as specified is either valid but powerless, or powerful but invalid (my toy simulation below). (2) Separating kappa from delta is identified only by priors for exactly the genes that matter: lineage genes expressed only in malignant cells, and immediate-early genes that are induced in TME cells too. (3) The matched data behind every invariance claim are tiny. Trinks has 3 LUAD patients (verified this session), Janesick 1 block, Oliveira Flex n=8 / Xenium n=4. Even a perfectly consistent tumour-level test cannot reach e-BH thresholds at m=1,600 hypotheses with fewer than about 8-14 matched tumours per MP per axis. A revised version could reach likely acceptance if it does four things: (a) runs equivalence-based typing at the tumour level; (b) calibrates nulls on same-protocol multi-region replicates; (c) restricts D2 claims to calibrator-identified genes; (d) builds the new experiment into a pre-registered matched cohort of at least 20-30 tumours.

**Fatal flaws**

- REVERSED BURDEN OF PROOF. T1 ('invariant plastic', the label that certifies 'real biology') means 'passes D2 and D3 on every testable axis', i.e. NO e-BH rejection. e-BH controls false NON-invariance calls; it gives no evidence FOR invariance. Power is set by the number of matched tumours per MP per axis, which is mostly 0-5, so under-powered MPs default to T1. A noise program tested on one 3-patient axis becomes 'invariant biology'. As written, the typing rewards a lack of data. D2/D3 must be recast as equivalence hypotheses with pre-specified margins (e.g. lower CI of transported cosine ≥0.9; |Δπ|≤0.1; |log Δm|≤log 1.5) and a minimum of positive, tumour-replicated evidence. Otherwise the MP is 'insufficient evidence', not T1.
- THE CONFIRMATORY E-VALUE IS EITHER INVALID OR POWERLESS AS SPECIFIED. (a) Universal-inference validity needs the denominator to maximise over the WHOLE null set. Fixing w_k, kappa and all other components at their X1 estimates maximises only over usage, so E[E]≤1 is no longer guaranteed. The test then also rejects nuisance mis-estimation. (b) A Poisson working likelihood on NB data is not the model's likelihood, so 'exact for Poisson' does not apply. The true NB likelihood with identity link is not concave in mu: at x=0, d²ℓ/dμ² = φ/(μ+φ)² > 0. So the 'convex, therefore global null MLE' argument fails exactly where it is needed. (c) Per-cell usage parameters are incidental. The denominator refits them on X2 and the numerator cannot, so log E pays a per-cell penalty. (d) Screening by the Rao score on the full data, X1+X2, and then confirming on X2 makes the selection depend on X2, which invalidates e-BH over the screened set. TOY CHECK (my script ref2_invarmp_evalue.py in the scratchpad: 1 tumour, 3,000+800 cells, 300 genes, K=3, NB φ=3, exact NB thinning with true φ, 6 replicates per row, KL-NMF instead of SVI). Results: H0, log E ≈ -8,859 (Poisson) / -2,652 (NB), so valid but degenerate. A real program-specific 1.65-fold shift on 20 of the program's 40 genes gives log E ≈ -5,618 / -1,776, i.e. 0/6 detections. Only a 2.7-fold shift is detected. A plug-in variant that recovers power (X1 usage on both sides) detects the 1.65-fold shift (log E +910 NB). But under H0 with only an unmodelled shared gene-wise capture error (kappa sd 0.2), the same plug-in gives log E ≈ +1,417 (NB), i.e. 6/6 false rejections. The penalty grows with the number of active components per cell: about 60 in the proposed model vs 3 in the toy. Real data will therefore be worse.
- KAPPA vs DELTA IS NOT IDENTIFIED FOR THE GENES THAT DEFINE MANY MPs. mu contains kappa_{g,e}·exp(δ_kg^e)·w_kg. For any gene with negligible expression in the calibrator cells (T, fibroblast, endothelial, myeloid), eps_{g,e} and δ_kg^e multiply the same term. The efficient Fisher information for δ_kg is then ≈0 apart from the priors (s² vs horseshoe τ_δ), so the D2 verdict for these genes is set by the prior. Such genes include epithelial/lineage genes (KRT, SFTP, SCGB, MUC), neural-lineage genes and melanocytic genes. Replication-dependent histone mRNAs are a concrete example (my reasoning): they are non-polyadenylated, captured differently by oligo-dT 3' chemistry and Flex probes, and expressed in few calibrator cells. A canonical cell-cycle MP could therefore be flagged as non-invariant, which defeats the Tier C target. The opposite failure hits the flagship. Neuschulz et al. (MSB 2023, verified) show the dissociation response has a core shared across tissues plus cell-type-specific programs, i.e. TME cells carry the very IEG response. If calibrators lack their own environment-owned components, TME FOS/JUN/HSPA1A induction is absorbed into kappa_{warm}. The malignant stress artefact is then 'transported away' and looks invariant (false T1). A1 is violated precisely in the genes behind the headline result.
- THE PAN-CANCER CLAIM HAS NO DATA TO STAND ON. 'Replaces manual curation at pan-cancer scale' requires bridging across environments for most MP × cancer-type pairs. In 3CA most cancer types come from one technical environment, mostly fresh 3' 10x. Protocol is chosen for biological or clinical reasons: single-nucleus for frozen biobank and brain tissue, Flex for FFPE archives and trials, Smart-seq2 for early cohorts. So the unmatched test (rbar ⟂ A | B) tests that there is no unmeasured confounding by tissue source, grade, treatment or site. It does not test invariance. No uniformly valid conditional-independence test exists in this setting (Shah & Peters, Ann Stat 2020, from memory [U]). The realistic outcomes are that T6 dominates the catalogue, or that confounded differences are typed T3/T5. Either way the 'manual curation replaced' headline fails.

**Major concerns**

- NULL TWINS DO NOT REPRESENT THE REAL NULL. Thinning one environment into two pseudo-environments yields the SAME cells, so loadings, selection and intensity are identical by construction. That calibrates only sampling noise. Real twins are different tissue pieces or sections, and regional ITH shifts π and m with no technical cause, e.g. hypoxia or MES states near necrosis. D3 (p(θ|t,e)=p(θ|t,e')) is therefore false under the technical null, and Δπ/Δm calls will be anti-conservative. The pilot P1 GO rule ('shift beyond 99th percentile of the thinning null') will flag almost every MP, and that will be misread as 'transport inadequate'. Calibration must come from same-protocol multi-region or multi-piece replicates, plus same-specimen technical replicates.
- POINT NULLS ON CELLS AS REPLICATES. With N=20,000 cells, mean count 0.5 and φ=2, the cell-level SE of a gene-wise log fold change is ≈0.011. A 5% kappa misspecification then gives z≈4.4 per gene (my calculation), and aggregating over 50 program genes is decisive. Large matched sets will reject δ=0 for almost every MP from transport error alone, while small ones have no tumour-level power. The estimand must be defined at the tumour or tumour×region level, with equivalence margins.
- E-BH ARITHMETIC AT THE CORRECT UNIT. m = 200 programs × 8 axes = 1,600. At q=0.1, e-BH needs e≥16,000 for 1 rejection, 1,600 for 10, 160 for 100 and 40 for 400. The largest valid e-value from an exact sign-flip test on n paired tumours is 2^n: 32 for n=5, 1,024 for n=10, 16,384 for n=14. Even with perfect consistency across tumours, the public matched data cannot reach these thresholds. Either pre-register a small hypothesis family (about 10-30 MP-axis pairs) or expand the matched cohort substantially. Averaging e-values from different evidence sources is valid but dilutes: one source with E=100 averaged with two at E=1 gives 34.
- TME SENSORS ARE NEITHER CLEAN NOR EXOGENOUS. (i) Several canonical IEGs and heat-shock genes are intronless (JUN, JUNB, IER2, HSPA1A/B; from memory, high confidence [U]), so their unspliced fraction is undefined. With 3' chemistry, intronic reads come largely from internal priming. (ii) Cell Ranger ≥7.0 counts intronic reads by default (verified), and a 2024 bioRxiv documents the effects of Cell Ranger version (title verified). Intronic fraction is therefore a processing-pipeline attribute that is missing from EnvSchema. (iii) Flex and in situ assays produce no intronic reads, and nuclei are intron-rich for all genes, so missingness and level of the sensor encode the environment itself. (iv) In vivo IEG induction (activated T cells, hypoxic or inflamed niches) moves the anchor, so A is a descendant of biology. Anchor regression's guarantees assume exogenous anchors in a linear SEM; with a non-exogenous anchor the penalty removes real biology. (v) Neuschulz 2023 found cell-type-specific dissociation programs, so the TME burden is at best a noisy proxy for malignant burden.
- ANCHOR PENALTY AND POWER POSTERIOR ARE TUNED BY CRITERIA THAT DEFEAT THEM. lambda_A='cv' is ill-defined: held-out deviance always prefers lambda_A=0, and anchor regression needs a pre-specified shift strength γ. Penalising Π_A rbar_k inside a joint factorisation can rotate programs and mix real biology that correlates with the anchors. ζ chosen by count-split deviance has the same problem: thinned folds carry the same systematic artefacts and misspecification, so held-out likelihood rewards keeping 'junk' components and pushes ζ toward 1, which nullifies the BayesPowerNMF rationale.
- OWNED vs SHIFTED IS PRIOR-DRIVEN, AND THIS IS THE FLAGSHIP CASE. A warm-only program can be represented either as an environment-owned component (T5) or as a shared component with π_{k,t}→0 and η^π_{k,warm}→+∞ (T3), and the priors decide which. The proposed split of 'Stress' into a warm-owned IEG part and an invariant HSP/DNAJ part is not identified within warm data for genes induced both in vivo and ex vivo (FOS, ATF3, HSPA1A/B). It requires measurements in a reference environment from the SAME tumours, which public data barely provide. Also, Marsh et al. (Nat Neurosci 2022, verified) found an ex vivo-like activation signature in postmortem human snRNA-seq. Pre-analytic warm ischaemia therefore affects nucleus and fixed arms too, so invariance across those arms does not exclude it. Cold-ischaemia time and digestion duration are absent from the schema.
- HURDLE SELECTION vs DETECTION. Nuclei carry about 1-3k UMIs and Flex depth differs from 3', so a low-intensity program falls below detectability in more cells. Δπ≠0 ('selection') will then be called on the cell↔nucleus and 3'↔Flex axes from depth alone. Amortised encoders with relaxed-Bernoulli and straight-through gradients are biased in a depth-dependent way. A null twin with a depth-mismatched split (e.g. 80/20) is needed to show that spurious Δπ is controlled.
- THE GENOTYPE AXIS IS CIRCULAR WHERE THE DATA ARE. Allele-aware Numbat and CalicoST need BAMs or SNP pileups, which are controlled-access or absent for most 3CA/TCCA studies. The fallback, expression-derived CNA, makes D4 and T2p tautological: an arm-level expression block is what such callers detect. The clone-permutation null treats cells as exchangeable, but clones are spatially segregated, so clone and region are confounded and the null is anti-conservative. Kinker's 8% association rate for continuous programs [V] suggests T2 will be small. DNA-based clone ground truth (clonealign-type paired scDNA/scRNA) covers few tumours in 1-2 cancer types, so the Cohen's kappa ≥0.7 target would rest on a small n.
- MUCH OF THE FULL MODEL CANNOT RUN ON MOST OF THE ATLAS. The ambient term alpha_c·b_{g,s} needs raw matrices with empty droplets, sensor (i) needs FASTQ/BAM re-alignment, and allele CNA needs BAMs. The plan re-processes only 20-40 key studies of 124 (3CA) plus 36 (TCCA). The remaining majority will run without these terms. That needs an explicit missing-module design, and results must be shown to be stable when those terms are dropped.
- VALIDATION IS LARGELY SELF-CONSISTENCY. Tier A plants artefacts in InvarMP's own functional form (gene-wise multiplicative kappa, additive non-negative environment components, hurdle usage), so recovery only shows the optimiser works when the model is true. The real malignant backgrounds already contain real programs and real protocol artefacts, so 'nuisance-only null cohorts' cannot be built from them without synthetic resampling, and simulator fidelity then governs every FDR claim. Tier C labels are 7 low-quality, 3 contamination and 1 mis-annotation, which are not protocol artefacts. The cell-type library or trivial baselines (cosine to TME markers, %mito/ribo, scDblFinder) would likely flag them, so these baselines must be reported. There is also a risk that real TME-mimicking malignant programs are typed as contamination: MHC-II/IFN, fibroblast-like pEMT, myeloid-like GBM-MES.
- THE VALIDATION TARGETS ARE STATISTICALLY UNDER-POWERED. AUROC 0.85 with 5 positives and 40 negatives has a 95% CI of about [0.63, 1.0]; with 15 positives and 60 negatives, [0.72, 0.98] (Hanley-McNeil). A Spearman of 0.6 on n=10 tumours has a 95% CI of [-0.05, 0.89]; n≈30 is needed for a lower bound above 0.3. In Tier D, correlations pooled over warm and cold arms are inflated by the two-group structure: 4sU incorporation at about 6°C is near zero, so the cold arm sits at the origin for both variables.
- THE 4sU EXPERIMENT IS NOT AS NOVEL OR AS CLEAN AS CLAIMED. Neuschulz et al. (MSB 19:e11147, verified) already add 4sU to the dissociation buffer (30 min at 37°C, methanol fixation, IAA, standard 10x) and report a core response plus cell-type-specific responses. The novelty narrows to human tumours with matched Flex/Xenium, and 'first' is unverified. Labelled fractions of non-IEG genes will sit near the T>C background, so 'T5 ≥3× T1 median' is a ratio of near-zero quantities. Use mixture-model (GRAND-SLAM-type) estimates and differences instead. 4sU itself perturbs stress and nucleolar pathways at high doses (from memory [U]). Splitting fresh specimens into two dissociation arms plus FFPE plus a serial section needs large resections. Protocol optimisation alone likely takes more than '6 months part-time'.
- SCALABILITY: PER-EPOCH CLAIMS ARE FINE, THE TOTAL BUDGET IS NOT. For NB, μ must be dense, because the zero term -φ·log(1+μ/φ) does not factorise. That is about 10^7 cells × 1.2×10^4 genes × 60 components ≈ 7×10^12 MACs per epoch, 8× the stated 10^12, but only seconds of H100 GEMM. Sparse I/O of about 120 GB per epoch from Zarr/dask and slow SVI convergence with IBP and horseshoe priors are the real per-fit limits. The budget fails on the number of fits. The simulator grid is 4×4×4×3 cells × 50 replicates = 9,600 cohorts per background, or 38,400 for 4 backgrounds. At 0.5-2 GPU-h per InvarMP fit that is 1.9×10^4 to 7.7×10^4 GPU-h, before about 14 baselines, 5 seeds and leave-one-study-out runs. The stated budget is 5-10k GPU-h.
- CALIBRATION OF THE BAYESIAN OUTPUTS. Mean-field SVI with amortised encoders and relaxed-Bernoulli tumour membership is typically over-confident. The coverage target for 90% membership sets (0.85-0.95) is unlikely without post-hoc recalibration, e.g. thinning-based or bootstrap-based.
- THE TYPING IS A NEW THRESHOLD STACK. ρ_gen ≥0.3, clone e ≥20 and transported cosine <0.5 are new arbitrary cut-offs. Non-negative 12k-dimensional loadings share highly expressed genes, so unrelated programs often have substantial cosine; 0.5 must be calibrated against a null of unrelated programs. The six-way type label has no error control of its own. e-BH covers only the non-invariance calls.
- BATCH IS CONFOUNDED WITH PATIENT, AND WITHIN-LEVEL VARIATION IS IGNORED. Environment levels are discrete, but block age and DV200, fixation time, digestion duration and ischaemia time vary within a level. One tumour usually means one library, so sample-specific gene-wise effects are absorbed by the tumour-private library or by shared components. Membership z_{k,t} for whole-tumour programs is then confounded with batch (landscape G5 risk).
- COUNTERFACTUAL OVERREACH. Per-cell 'environment-standardised usage' is not identified under selection: cells lost to warm dissociation have no counterparts. Only population-level standardisation is defensible.
- REUSE OF DATA ACROSS STAGES. The same few matched sets (Trinks, ccRCC, Slyper, Oliveira) serve as pilot data, simulator backgrounds, training data for the sensor classifier and kappa fits, and Tier B validation. A locked held-out matched set must be designated before any tuning.
- TEAM AND TIMELINE ARE UNREALISTIC. One postdoc, one PhD student and a 20-30% statistician are asked to deliver: a 6-library Bayesian model with latent environments and an anchor penalty; two open theory problems (validity of NB thinning in factor models, the Hivert caveat; noisy identifiability of the shared/owned split); STARsolo re-processing of 20-40 studies; about 14 fairly tuned baselines; a simulator on the scale of 38k cohorts; a wet-lab resource; and a polished Python/R package with a versioned reference. I estimate 36-48 months, or the addition of a research software engineer, a data engineer and a statistician at ≥50%.

**Prior-art threats**

- Neuschulz A, ... Junker JP. 'A single-cell RNA labeling strategy for measuring stress response upon tissue dissociation', Mol Syst Biol 19:e11147 (2023). 4sU added during dissociation; core and cell-type-specific dissociation programs. Undercuts the 'first 4sU dissociation dataset' framing (only the human-tumour part may be new) and the TME-sensor/A1 assumptions. [VERIFIED this session]
- Marsh SE et al. 'Dissection of artifactual and confounding glial signatures...', Nat Neurosci 25:306-316 (2022). Ex vivo activation signature; transcription/translation inhibitor cocktail as ground truth; the signature is also present in postmortem human snRNA-seq, so nucleus/fixed arms are not artefact-free. [VERIFIED this session]
- Trinks et al. 2024 (bioRxiv 2023.04.25.538273; journal per PMC). Fresh / cryopreserved / FFPE Flex on only 3 LUAD patients (P075, P078, P079). Already reports concordance across protocols and is too small to carry invariance tests. [VERIFIED: patient count]
- ccRCC fresh sc / frozen sn-Multiome / FFPE snFlex benchmark, bioRxiv 2025.11.25.690485. Already reports fine-grained tumour states as concordant across protocols, partly pre-empting the 'invariance audit' findings. [VERIFIED existence; patient n NOT retrieved]
- GBCD (Liu…Stephens, Nat Genet 2025): shared vs patient-specific tumour programs with lfsr. [V per landscape]
- ITHmapper / Pan-Cancer ITH compendium (Abelson lab, bioRxiv 2026): >1,000 tumours, scVI + Hotspot, mapping tool. [V per landscape]
- TCAT/starCAT (Kotliar…Raychaudhuri, Nat Methods 2025): catalogue plus projection; README warns about transfer to Flex. The reference/projector part is derivative. [V per landscape]
- Grabski, Trippa, Parmigiani, Bayesian multi-study NMF (Genome Biol 26:98, 2025), plus Tetris/APAFA/SUFA. The statistical core of subset sharing with inclusion indicators. [V per landscape]
- inVAE (Aliee…Theis, bioRxiv 2024): invariant vs spurious latent split in single-cell data. The closest 'invariance' concept. [V per concept text; not re-checked by me]
- sciRED (Pouyabahar, Andrews, Bader, Nat Commun 2025): factors labelled as covariate/technical vs unexplained, i.e. automated typing of technical factors. [V per landscape]
- Anchor regression (Rothenhäusler et al., JRSS-B 2021) and invariant causal prediction (Peters, Bühlmann, Meinshausen, JRSS-B 2016). Their assumptions (exogenous anchors, environments that do not act on the target) are not met here, so citing them invites scrutiny. [V per concept; not re-checked]
- Dissociation-signature literature: van den Brink 2017 (Nat Methods), O'Flanagan 2019 (Genome Biol), Denisenko 2020 (Genome Biol), Slyper 2020 (Nat Med). Known IEG/HSP artefact signatures and protocol contrasts. [V per concept; the landscape lists them as U]
- 'Impacts of Cell Ranger versions on Chromium gene expression data', bioRxiv 2024.08.10.607413. Processing version as a technical environment. [VERIFIED title existence]
- Hivert et al., arXiv 2405.13591: thinning breaks under mixtures when overdispersion is mis-estimated. Directly relevant to NB thinning in a factor model. [V per landscape]
- Shah & Peters, 'The hardness of conditional independence testing...', Ann Stat 2020. No uniformly valid CI test, relevant to test (c). [UNVERIFIED, from memory]
- Wu YE et al., Neuron 2017: actinomycin D during dissociation to suppress IEG artefacts. [UNVERIFIED, from memory]
- No direct scoop on 'invariant metaprograms' or 'protocol-invariant gene programs' turned up in a 2026 web search this session. The search was shallow (bioRxiv/PMC were proxy-blocked), so this is weak evidence.

**Required changes**

- Recast D2/D3 as equivalence hypotheses with pre-registered margins (e.g. lower 95% bound of transported cosine ≥0.9; |Δπ|≤0.1; |log Δm|≤log 1.5). Assign T1 only on positive, tumour-replicated evidence on ≥1 bridged axis; anything else is 'insufficient evidence' (T6), never T1.
- Move inference to the tumour (or tumour × region) level. Use paired NB-GLM pseudobulk tests or tumour-level permutation/sign-flip tests. Provide a power analysis giving the number of matched tumours needed per MP-axis pair, and either shrink and pre-register the hypothesis family (about 10-30 pairs) or justify m=1,600.
- Rebuild the e-value: screen on X1 only; use the NB likelihood throughout, with a global-optimum certificate or a design where the restricted MLE is reliable (pseudobulk NB-GLM); eliminate per-cell incidental usage (e.g. integrate it out under a fitted usage prior, or aggregate). Show empirical validity under ESTIMATED φ and mixture structure (Hivert scenario) and non-trivial power for 1.5-2-fold deviations with realistic per-cell nuisance (about 60 active components).
- Calibrate the technical null with real biological replicates: same-protocol multi-region or multi-piece pairs and same-specimen technical replicates. Report variance components (regional vs technical) and re-specify the pilot P1 GO rule against this null, not against the thinning null.
- State and enforce identifiability for D2. Test only genes whose kappa is pinned by calibrators, or test deviations relative to other programs that share the gene. Report per MP the fraction of loading mass that is calibrator-identified and a per-gene efficient-information diagnostic. Show how sensitive the types are to s² and τ_δ.
- Give calibrator cells their own environment-owned and cell-type-specific ex vivo components; exclude IEG/HSP genes from kappa estimation; test A1 on twins before using it. Show that a planted TME-plus-malignant dissociation response is NOT transported away (no false T1).
- Make the new experiment the evidentiary core: n ≥20-30 tumours in ≥3 cancer types. Include same-protocol replicate pieces, a Marsh-style transcription-inhibitor arm alongside 4sU, and recorded pre-analytic covariates (ischaemia time, digestion time, DV200). Pre-register predictions and lock a held-out set that is never used to tune the sensor classifier, kappa or the simulator.
- Extend EnvSchema with processing pipeline/version and pre-analytic continuous variables, and re-process uniformly where possible. Report results separately for data with and without the ambient, sensor and allele-CNA modules.
- Give Tier C trivial baselines (marker cosine, %mito/ribo, scDblFinder score) and CIs on every AUROC. Separate artefacts that invariance can catch (protocol-specific) from those it cannot (low quality, doublets, pre-analytic ischaemia). Test explicitly for false 'contamination' typing of TME-mimicking malignant programs.
- Build null cohorts and misspecified artefacts outside InvarMP's own generative form: scDesign3/ClusterDE-style nulls; artefacts that interact with programs, act through selection or are nonlinear. Replace the 38,400-cohort full factorial with a fractional or space-filling design of about 10^3 cohorts and re-budget compute honestly.
- Genotype axis: call T2/T2p only with allele- or DNA-based clones; use a region-aware or block permutation null; report circularity sensitivity (expression-CNA vs allele-CNA calls on the same tumours).
- Choose lambda_A via a pre-specified worst-case shift γ, or replace anchor regression with a plain conditional-independence test plus sensitivity analysis. Choose ζ by a misspecification-calibrated criterion, not held-out deviance on thinned folds.
- Calibrate the posterior outputs (tumour-membership sets, lfsr gene sets) empirically and report coverage; drop the per-cell counterfactual 'standardised usage' or restrict it to population-level standardisation.
- Calibrate the T5 cosine threshold against a null of unrelated programs; give error rates for the type labels themselves, not only for the non-invariance calls.
- Re-scope team and timeline: add a research software engineer and a data engineer, raise the statistician to ≥50%, plan 30-36 months, or split the work into a methods paper on public data and a separate resource paper for the 4sU/Flex/Xenium cohort.
- Run the bioRxiv/Scholar scoop search (landscape section 6.3 plus 'protocol-invariant programs', 'Flex single-nucleus metaprogram', citing papers of Neuschulz 2023 and Marsh 2022) and a name-collision check for 'InvarMP'.

**Ideas worth keeping**

- The core reframing: an MP should be judged by stability under technical INTERVENTIONS on the same tumour, not by recurrence. It is falsifiable and attacks the exact failure that the 35/10/10 filters cannot (systematic artefacts recur by construction).
- Returning artefacts as typed, reusable environment-owned components (an 'artefact library') instead of deleting them by gene blocklists or hand curation.
- Separating protocol-driven selection (Δπ) from induction (Δm). This lets fragile-cell biology be kept as biology rather than discarded as artefact.
- The design-graph testability report (T6): an honest map of where protocol is confounded with cancer type. Valuable on its own as an audit of 3CA/TCCA, even if the inferential layer is simplified.
- Gene-wise transport pinned by non-malignant calibrator cells, with panel-specific Fisher information for projection into Flex, Xenium 5K and CosMx 6K. A practical deliverable that starCAT explicitly lacks, as long as claims are restricted to calibrator-identified genes.
- The ex vivo vs in vivo stress decomposition as the flagship biological question. With a properly designed matched cohort (4sU plus inhibitor arm plus Flex plus Xenium on the same specimens) it would be a high-impact, falsifiable result.
- Using the 3CA v2 hand-removed MPs and the protocol-named MPs of the 2026 catalogue as labelled test cases, together with explicit GO/NO-GO pilots, once the pilots are re-calibrated against multi-region same-protocol nulls.


---

## 4. CERT-MP: Certified Estimation of Recurrent Tumour MetaPrograms (Python package `certmp`, scverse-native, with a thin R/Seurat bridge `certmpR`)  
*Design angle: `feasibility-first`*

### Claim

CERT-MP turns tumour metaprograms from threshold-defined top-50 gene lists into certified statistical objects. It runs a label-blind discovery pipeline on real tumours alongside matched "decoy" tumours. Decoys keep each tumour's copy-number dosage, ambient RNA, doublets, quality structure and depth but contain no programs. Comparing held-out program recurrence in real versus decoy tumours controls the false discovery rate of metaprogram existence, gives calibrated tumour and gene membership, and automatically types every program as non-genetic, dosage-driven or technical. The aim is to show, without manual curation, which pan-cancer metaprograms are real.

### Formal definition of an MP

SETTING
- Tumours t = 1..T come from studies s(t). Malignant cells c are in C_t. X_tcg are UMI counts (read counts flagged separately) over a gene universe G of about 6,000-8,000 genes. G is chosen as the union of per-tumour top genes by binomial deviance, not by mean expression (removes the S1 bias).
- A DECLARED NUISANCE SET N_t contains:
  - library size s_tc;
  - QC covariates q_tc (log %mito, log %ribo, complexity residual);
  - clone label k(c) and clone-level CNA log-ratio profile delta_{t,k,g}, from infercnvpy, SCEVAN or CopyKAT, and Numbat where allele data exist;
  - sample ambient profile rho_tg and per-cell contamination alpha_tc;
  - heterotypic doublet / misclassification rate pi_t, using that sample's real non-malignant cells;
  - a tumour baseline floor b_tg;
  - NB dispersions r_tg.
- A gene program is a non-negative vector w in R_+^G with ||w||_2 = 1. Gene sets are allowed as w = indicator/sqrt(n).

(1) ACTIVITY (tumour membership)
- Each tumour is split by negative-binomial data thinning into independent folds X^A (discovery) and X^B (inference).
- Programs are learned only on X^A.
- On X^B, compute Pearson residuals r_tcg under N_t, then the held-out weighted coherence statistic Z_t(w) (eq. M5 in model_and_math).
- The null is H0_t(w): the genes of w are conditionally independent given N_t in tumour t.
- w is ACTIVE in t if H0_t(w) is rejected at Benjamini-Hochberg level alpha (default 0.1) across tumours.
- Membership set: A(w) = {t : q_t(w) <= alpha}. Tumours with fewer than 50 malignant cells are reported as "insufficient cells", not as negatives.

(2) DECOY TUMOURS
- For every real tumour t, generate a decoy t~ in family D^full (all of N_t, no programs).
- Also generate one in D^tech (N_t with CNA dosage removed, i.e. clone-flat).
- Each decoy inherits s(t) and passes through every step, including engine fitting, without the engine knowing which tumours are decoys.

(3) METAPROGRAM AT LEVEL q
w is a (non-genetic) MP if both hold:
- (a) A(w) contains tumours from at least 2 studies. This is the only definitional constant: "meta" means cross-study. It is user-settable.
- (b) w is selected by knockoff+ at FDR q on the antisymmetric contrast W^full(w) = R(w) - R~^full(w).
  - R(w) = sum_s min(1, |A(w) ∩ s| / 2) is study-capped recurrence over real tumours.
  - R~^full is the identical statistic over D^full decoys.

A program that fails (b) but is selected on W^tech(w) = R(w) - R~^tech(w) is a RECURRENT DOSAGE (genetic) PROGRAM. Its recurrence is explained by CNA but not by technical nuisance. It is labelled, not deleted. A program that fails both is TECHNICAL (if rho^tech = R~^tech/R >= 0.5) or UNSUPPORTED.

(4) ATTRIBUTES OF EACH MP
- GENE MEMBERSHIP G_q(w): genes whose held-out association with the leave-one-gene-out program score, meta-analysed over A(w), has BH q <= q_gene (default 0.05). Each gene carries an effect size beta_g and heterogeneity I^2_g (core vs context genes).
- MODIFIERS:
  - CLONE-LINKED: median eta^2 of clone on the held-out score > 0.3.
  - PROTOCOL-SENSITIVE: CMH enrichment of A(w) in one protocol (sc / sn / Smart-seq2 / Flex) stratified by cancer type, odds ratio > 3 at q < 0.05, or loss in matched dissociation-free data.
  - QC-ENTANGLED: 0.2 <= rho^full < 0.5.
  - POSITIONAL: chromosome-arm enrichment of G_q(w).
- STABILITY: matched cosine across thinning replicates.
- ITH FRACTION: within-tumour vs within-study between-tumour variance.

PLAIN-LANGUAGE DEFINITION. An MP is a non-negative gene program whose within-tumour co-variation replicates on independent data in tumours from at least 2 studies, and does so significantly more often than in matched tumours that contain every declared nuisance but no programs, at a stated FDR. Certification is explicitly relative to N, and N is extensible: any measured artefact signature can be added as a decoy family. This makes "is it real biology?" a calibrated, auditable claim rather than a curator's judgement.

### Model and math

M1. NUISANCE MODEL (per tumour, embarrassingly parallel)
- mu0_tcg = s_tc * [ b_tg * exp(kappa_t * delta_{t,k(c),g} + gamma_tg' q_tc) + alpha_tc * rho_tg ]
- X_tcg ~ NB(mean mu0_tcg, size r_tg)
- r_tg is estimated CONDITIONALLY on a rank-10 Poisson-NMF fit of the same tumour, then empirical-Bayes shrunk on log-mean. This avoids the inflated marginal overdispersion under state mixtures that Hivert et al. (arXiv 2405.13591) [V] show breaks thinning.
- gamma_tg is ridge-penalised, with lambda chosen by GCV.
- kappa_t comes from regressing clone pseudobulk log-ratios on delta. Clone-level, NOT cell-level, CNA is used to limit circularity.
- rho_tg is the normalised pooled profile of all cells in the sample (DecontX-like [V]; no empty droplets needed).
- alpha_tc is estimated from TME-marker genes (genes at least 10x higher in that sample's non-malignant cells).

M2. DECOY GENERATOR
- For each real cell, copy (s_tc, q_tc, k(c), alpha_tc). This preserves the joint law of cell-level nuisance covariates.
- Draw X~_tcg ~ NB(mu0_tcg, r_tg).
- For a fraction pi_t of cells (from scDblFinder [V] rates and the TME-marker tail), add a heterotypic doublet / misclassification component: X~_tc <- Thin(X~_tc, 1 - omega) + Thin(X^TME_{t,c'}, omega), where c' is a random REAL non-malignant cell from the same sample and omega ~ Beta(2,2).
- D^tech sets kappa_t = 0 (clone-flat) and keeps everything else.
- Decoys contain zero within-tumour programs by construction, and all declared nuisance co-variation.

M3. NB DATA THINNING (Neufeld et al. JMLR 2024 [V]; NB count splitting arXiv 2307.12985 [V])
- X^A_tcg | X_tcg ~ BetaBinomial(X_tcg, eps*r_tg, (1-eps)*r_tg); X^B = X - X^A.
- When r is correct: X^A ~ NB(eps*mu, eps*r) and X^B ~ NB((1-eps)*mu, (1-eps)*r), independent.
- Defaults: eps = 0.6, with R = 3 independent thinning replicates (different seeds).

M4. DISCOVERY ENGINE CERT-NMF
Fit on fold A of the LABEL-BLIND AUGMENTED COHORT {real} ∪ D^full ∪ D^tech.
- X^A_tcg ~ Poisson(mu_tcg), likelihood tempered by zeta.
- mu_tcg = eps * s_tc * { Omega_tcg * [ b_tg + sum_{k<=K} H_tck W_kg + sum_{j<=J} U_tcj V_tjg ] + alpha_tc * rho_tg }
- Omega_tcg = exp(kappa_t delta_{t,k(c),g} + gamma_tg' q_tc), fixed from M1. Dosage and QC enter as multiplicative offsets, ambient as an additive fixed term.
- Components:
  - W (K x G, K_max = 200, rows sum to 1): shared dictionary;
  - H >= 0: usages;
  - V_t (J = 3) >= 0: tumour-private programs, which absorb clone- or tumour-specific structure so it does not contaminate W;
  - b_t >= 0: tumour baseline floor, which replaces per-tumour centring and zeroing.
- Objective: minimise zeta * sum KL(X^A | mu) + sum_{t,k} [ (||H_{t.k}||_1 + beta0)/lambda_tk + (n_t + a0 + 1) log lambda_tk ], with the same automatic relevance determination (ARD) on private columns (L1-ARD after Tan & Fevotte, IEEE TPAMI 2013 [U]).
- MM / multiplicative updates:
  - H_tck <- H_tck * [sum_g W_kg Omega_tcg X^A_tcg/mu_tcg] / [sum_g W_kg Omega_tcg + 1/(zeta*eps*s_tc*lambda_tk)]
  - W_kg <- W_kg * [sum_{t,c} s_tc H_tck Omega_tcg X^A_tcg/mu_tcg] / [sum_{t,c} s_tc H_tck Omega_tcg], then renormalise rows and rescale H
  - b_tg <- b_tg * [sum_c s_tc Omega_tcg X^A/mu] / [sum_c s_tc Omega_tcg]
  - lambda_tk = (||H_{t.k}||_1 + beta0)/(n_t + a0 + 1)
- Online, tumour-minibatch version with running sufficient statistics for W (in the style of online iNMF [U]).
- zeta in {0.25, 0.5, 1} and K_max in {150, 250} are chosen by held-out NB deviance on a nested thinning of fold A.
- Under MAP, tempering is equivalent to rescaling the prior. It is used as an out-of-sample-tuned regulariser, NOT claimed as power-posterior novelty (cf. BayesPowerNMF [V]).
- Programs are matched across the 3 replicates (Hungarian on cosine) and kept if reproduced in at least 2/3 with cosine >= 0.8. Stability index = mean matched cosine.
- Initialisation is deterministic (NNDSVDa on pooled normalised data).

M5. HELD-OUT COHERENCE TEST (works with any engine and any w)
- r_tcg = (X^B_tcg - m_tcg) / sqrt(m_tcg + m_tcg^2/((1-eps) r_tg)), with m = (1-eps)*mu0.
- S_tc = sum_g w_g r_tcg
- T_t(w) = sum_c [ S_tc^2 - sum_g w_g^2 r_tcg^2 ] = sum_c sum_{g != g'} w_g w_g' r_tcg r_tcg'
- V_t(w) = 2 * sum_c [ (sum_g w_g^2 r_tcg^2)^2 - sum_g w_g^4 r_tcg^4 ]
- Z_t(w) = T_t / sqrt(V_t)
- PROP 1: under H0_t(w) (genes conditionally independent given N_t, cells independent given nuisance), Z -> N(0,1) by the Lindeberg CLT over cells. Variance misspecification affects only the plug-in V, which is empirical.
- RANDOM-PROGRAM GENOMIC CONTROL: lambda_t = median of Z^2 over 500 random programs matched to w on gene-mean deciles and effective gene number 1/sum(w^4), divided by 0.455. Then Z <- Z / sqrt(max(1, lambda_t)). This removes global co-variation from depth or cell-size mis-modelling.
- ROBUSTNESS TO THINNING ERROR: w is learned on A, and the test uses only CROSS-GENE products within B. Entrywise A-B dependence from a mis-estimated r (Hivert et al.) therefore enters only at second order. This is the key design answer to the known critique.
- One-sided p-values. BH is applied SEPARATELY within each arm (real, D^full, D^tech) so the procedure stays symmetric.
- A "specific membership" variant uses only genes where w dominates: W_kg / sum_k' W_k'g > 0.5.

M6. CERTIFICATION (tumour-level knockoffs)
- R_k = sum_s min(1, |A_k ∩ s|/2) over real tumours; R~^f_k is the same over decoy arm f.
- W^f_k = R_k - R~^f_k + 1e-3 * (sum_{real} Z+_tk - sum_{decoy} Z+_tk), the second term only as a tie-break.
- W is averaged over the 3 thinning replicates, which preserves antisymmetry (derandomisation in the spirit of Ren-Wei-Candes [U]).
- Knockoff+ (Barber & Candes 2015 [U]): tau_q = min{tau > 0 : (1 + #{k: W_k <= -tau}) / max(1, #{k: W_k >= tau}) <= q}; select {k: W_k >= tau_q}. The number of MPs is therefore an FDR-controlled output, not a user K.
- PROP 2 (approximate FDR):
  - Assumption (decoy exchangeability): for programs whose recurrence is generated by nuisance, swapping each real tumour with its decoy leaves the joint law of the label-blind pipeline output invariant.
  - Consequence: null W_k have i.i.d. symmetric signs given |W|, and knockoff+ gives E[FDP] <= q.
- Exchangeability is an assumption, because real tumours also contain biology that can interact with nuisance in the fit. We enforce label-blindness, test exchangeability on nuisance features, and measure FDR inflation under deliberate misspecification. Results are reported as "empirical FDR".

M7. TYPING
- rho^f_k = R~^f_k / R_k.
- Labels:
  - CERTIFIED NON-GENETIC MP: selected on W^full;
  - DOSAGE PROGRAM: not selected on W^full, selected on W^tech;
  - TECHNICAL: not selected on W^tech and rho^tech >= 0.5;
  - UNSUPPORTED: R < 1.
- Modifiers:
  - clone-linked: median eta^2 > 0.3 from a one-way ANOVA of held-out S_tc on clone, BH q < 0.05;
  - protocol-sensitive: CMH test as above;
  - QC-entangled: 0.2 <= rho^full < 0.5;
  - positional: arm-level hypergeometric test of G_q (BH 0.05) plus the fraction of weight on the top arm.

M8. GENE MEMBERSHIP WITH FDR (replaces top-50)
- For MP k, member tumour t and gene g: Z_tkg = sum_c r_tcg S^(-g)_tc / sqrt(sum_c r_tcg^2 (S^(-g)_tc)^2), with S^(-g) = sum_{g' != g} w_g' r_tcg'.
- Effect beta_tkg comes from an NB GLM of X^B_tcg on standardised S^(-g), with offset m_tcg.
- DerSimonian-Laird random-effects meta-analysis over A_k gives Z_kg, beta_kg and I^2_kg. BH across all G genes gives q_kg (Benjamini-Yekutieli as a sensitivity check).
- CORE gene: q <= 0.05, beta > 0 in at least 70% of members, and I^2 < 0.5. CONTEXT gene: q <= 0.05 with I^2 >= 0.5; per-cancer-type effects are reported.
- When cells are plentiful, an optional three-way thinning (0.5/0.25/0.25) removes the mild effect of selecting member tumours.

M9. ITH VS BASELINE (keeps what centring discards)
- beta_t = argmin_{beta >= 0} KL(b_t | sum_k beta_k W_k + e_t) gives between-tumour usage.
- pi^ITH_k = V_within / (V_within + V_between), where V_within = mean_t Var_c(H_tck) and V_between is computed AFTER removing study means. This avoids the confound where one patient is one library.

M10. IDENTIFIABILITY DIAGNOSTIC
- Anchor count a_tk = #{c : H_tck / (baseline share + sum_k' H_tck') > 0.5}.
- MP anchor coverage = fraction of member tumours with a_tk >= 10.
- PROP 3 (noiseless, monotonicity; to be written formally): the sufficiently-scattered condition on cone(rows of the stacked H = [H_1; ...; H_T]) is preserved under adding rows. It holds if any single tumour satisfies it, and can hold when every tumour fails it. Pooling heterogeneous tumours can therefore only enlarge the set of identifiable shared programs. The statement is for the shared block with private components fixed. The sweeps did not find this argument in the literature (Fu/Huang/Sidiropoulos theory [V/P]).

M11. PROJECTION
- For a query sample: nuisance fit, then KL-NNLS of counts on [W_ref, b_query, 3 private slots].
- Per-cell SE from the Poisson Fisher information: I_c = sum_g s_c^2 Omega_cg^2 W_.g W_.g' / mu_cg, and SE_ck = sqrt((I_c^-1)_kk).
- Sample-level MP detection q-values from M5. These are valid without thinning because W_ref was learned on independent data.
- Novelty = M5 coherence of the private slots (fit on query fold A, tested on fold B).
- Panel transfer: restrict W to panel genes, report retention iota_k = ||W_k,panel||^2 / ||W_k||^2, and flag iota_k < 0.3 as "not assessable on this panel".

### Inference and scalability

DESIGN PRINCIPLE: ESTIMATE, THEN CERTIFY
- Estimation (M4) is fast MAP with GPU multiplicative updates. It does not need correct posterior uncertainty.
- All inferential claims come from held-out tests (M5, M8) and decoy contrasts (M6). Both are cheap and calibrated by construction or empirically.
- This avoids the two scale killers of the full-Bayes route: MCMC (Tetris and Grabski-style models do not scale), and variational-inference under-coverage of inclusion posteriors that reviewers would probe.

COST PER STAGE (estimates, to be benchmarked)
- (1) Nuisance fitting: batched IRLS in PyTorch per tumour, about 1 CPU-minute per 1,000 cells x 8,000 genes. About 30 CPU-hours for the roughly 1,500 malignant tumours of 3CA v1 plus the v2 additions.
- (2) Decoy generation and NB thinning: O(nnz), minutes for the full atlas.
- (3) CERT-NMF:
  - Training uses at most 1,000 malignant cells per tumour; all cells are scored by projection.
  - The 3x label-blind augmented cohort is roughly 4-4.5M cells x 8,000 genes, sparse.
  - Per epoch, the numerator costs O(nnz*K) at nonzeros only. The denominators sum_g W_kg Omega_tcg are dense because of QC offsets. For K = 200 this is about 2-3e12 flops per epoch, estimated at 10-20 minutes on one H100.
  - 50-100 online epochs gives about 15-30 GPU-hours per thinning replicate.
  - A full atlas run (3 replicates plus a 6-point hyperparameter grid on a 20% subset) is about 100-150 GPU-hours.
  - H is held in pinned CPU memory, and tumour minibatches are streamed. The dictionary W is 200 x 8,000.
  - Private V_t: 4,500 tumour-units x 3 x 8,000 floats, about 430 MB.
- (4) Membership tests: about 4,500 tumour-units x 200 programs = 0.9M statistics, each O(nnz_t), plus 500 random programs per tumour for genomic control. About 1-2 GPU-hours.
- (5) Gene tests: about 60 certified MPs x about 100 member tumours x 8,000 genes = 5e7 vectorised statistics plus random-effects meta-analysis. About 2 GPU-hours.
- (6) Knockoff+ and typing: seconds.
- (7) Projection: GPU KL-NNLS, estimated at more than 1e5 cells per second. Scaling is shown to at least 3M malignant cells (3CA v2 plus TCCA, uncapped projection) and to 1e7 cells on Tahoe pseudo-cells.

TOTAL COMPUTE BUDGET
- About 1,500-3,000 GPU-hours including the simulation grid and replicate runs.
- About 50k CPU core-hours for baselines. GBCD runs per cancer type with subsampling to 20k cells. The Gavish R pipeline runs on real plus decoy tumours times 3 replicates, about 7,000 per-tumour NMF jobs.

ROBUSTNESS KNOBS (reported, not tuned to taste)
- eps (0.5 / 0.6 / 0.7), K_max (150 / 250), J (2 / 3 / 5);
- CNA caller (infercnvpy / SCEVAN / CopyKAT);
- QC-slope family on or off;
- alpha (0.05 / 0.1);
- q (0.05 / 0.1 / 0.2).

### Outputs for users

PER METAPROGRAM: an "MP card" (HTML/PDF) plus tables in adata.uns['certmp']
- Certification:
  - knockoff+ FDR-level membership (the smallest q at which the MP is selected);
  - W^full and W^tech;
  - decoy reproduction ratios rho^full and rho^tech;
  - type label: certified non-genetic / dosage-genetic / technical / unsupported;
  - modifiers: clone-linked, protocol-sensitive, QC-entangled, positional, each with its continuous score.
- Membership: per-tumour Z and q-values; membership heatmap by cancer type, study and protocol; "specific membership" variant.
- Gene table: every gene's held-out effect size beta, SE, meta-analytic Z, q, I^2, core/context flag, and per-cancer-type effects. Exported as a confidence-tiered GMT (core q < 0.01 / core q < 0.05 / context) for MSigDB-style use.
- Quality indicators:
  - stability index (matched cosine across thinning replicates, seeds and leave-one-study-out);
  - anchor coverage (identifiability);
  - ITH fraction pi^ITH (within-tumour vs within-study between-tumour variance);
  - chromosome-arm positional plot;
  - cross-protocol flag where matched data exist.

PER TUMOUR
- Tumour-private programs with held-out coherence (candidate tumour-specific biology, e.g. clone programs), between-tumour baseline usage beta_t, faithfulness diagnostics of its decoy, and the genomic-control lambda_t.

PER CELL
- Usage H in adata.obsm['certmp_usage'], with Fisher-information SE in adata.obsm['certmp_se'] and a novelty residual.

REFERENCE RESOURCE
- A versioned Zenodo DOI, "certmp-pan v1.0". It contains W, the gene universe, nuisance-model templates, decoy-calibrated genomic-control tables, type labels and gene tables.
- A projector certmp.tl.project(query, ref='certmp-pan-v1', panel=None) returns:
  - per-cell usage with SE;
  - per-sample MP detection q-values (valid without thinning);
  - novelty flags (private-slot coherence);
  - panel retention iota_k, so Flex and Xenium users know which MPs are assessable.

AUDIT MODE
- certmp.tl.certify_external(adata, programs=...) certifies ANY MP list or engine output: Gavish/3CA lists, GeneNMF, cNMF, GBCD, ITHmapper modules or a user's own.
- If the list was derived from the same cells, the tool warns that an independent cohort or a re-run on fold A is needed.

BENCHMARK
- certmp-bench: the semi-synthetic tumour simulator (planted programs plus planted confounders plus null cohorts), metric suite, and containerised baselines. Usable by others to evaluate future MP methods.

### Failure modes fixed

Mapped to the Gavish/3CA pipeline steps (synthesis Section 1.2) and the named systemic failure modes (Section 1.4).

S0 MALIGNANT CALLING ERRORS, DOUBLETS, CONTAMINATION
- Decoys inject heterotypic doublets and misclassified cells, built from that sample's REAL non-malignant cells at the estimated rate, plus the sample's ambient profile.
- The engine also has a per-cell ambient term alpha_tc*rho_tg.
- Contamination "MPs" (the macrophage, T-cell and endothelial MPs hand-removed in 3CA v2) should then be equally recurrent in decoys, with W^tech near 0, and are typed TECHNICAL automatically.

S1 MEAN-EXPRESSION GENE UNIVERSE AND PER-TUMOUR CENTRING / ZEROING
- Deviance-based gene selection replaces top-7,000 by mean.
- A count likelihood with depth offsets replaces log2 centring.
- A tumour baseline floor b_t replaces centring and zeroing. It is decomposed on W (M9), so programs active across a whole tumour (subtype or baseline axes) are reported, not erased.
- Addresses failure mode 3 (per-tumour centring blocks joint inference).

S2 FIXED RANK GRID 4-9 AND BEST-OF-10 RESTARTS
- One joint dictionary with per-tumour ARD (programs switch on and off per tumour).
- zeta and K_max are chosen by held-out deviance on nested thinning.
- Stability across 3 thinning replicates (independent noise realisations, not just seeds) is REPORTED.
- An anchor and identifiability diagnostic is reported per MP.
- Pooling across tumours provably relaxes scatteredness for shared programs (Prop. 3).

S3 TOP-50 TRUNCATION
- Replaced by FDR-controlled held-out gene membership with effect sizes and heterogeneity (M8). List length follows the data.

S4 35/10/10 FILTERS
- intra_min=35 (cross-rank stability) is DROPPED as a filter, which removes its bias toward dominant axes, including technical axes (failure mode 4).
- inter_min=10 guards only against a random-gene null (P approx 1e-12). It is replaced by two tests:
  - the held-out coherence test, against overfitting;
  - the decoy knockoff contrast, against structured nuisance recurrence (failure mode 1: recurrence is necessary, not sufficient).
- intra_max=10 greedy pruning is replaced by a shared dictionary, so redundancy is handled in the fit.

S5 GREEDY CLUSTERING AND THE MP COUNT AS A BY-PRODUCT OF THRESHOLDS
- No clustering step: MPs are rows of a jointly fitted dictionary.
- The NUMBER of MPs is the knockoff+ selection at stated FDR.
- Study-capped recurrence and leave-study-out stability address cohort dependence (failure mode 5).

S6 MANUAL CURATION AND CNA CHECKED ONLY AFTERWARDS
- Typing by decoy family (technical vs dosage vs non-genetic), plus clone-linked, protocol-sensitive and QC-entangled modifiers, replaces inds_rm-style curation and post-hoc CNA checks.
- Nuisance is part of the DEFINITION (failure mode 6: validation is post hoc).
- Pre-registered calibration check: AUROC of the type labels against the 11 MPs hand-removed in 3CA v2 and the MPs Gavish 2023 removed (LQ / doublet / single-study).

S7 MODULE SCORES WITHOUT UNCERTAINTY, NOT COMPARABLE ACROSS TUMOURS
- Projection gives per-cell usage with SE, offset-adjusted (comparable across tumours and protocols), per-sample detection q-values, novelty flags and panel-retention warnings.

FAILURE MODE 2 (NO LIKELIHOOD, NO UNCERTAINTY)
- Replaced by a count likelihood, held-out p-values and q-values at three levels (MP existence, tumour membership, gene membership), and stability indices.

WHAT IS DELIBERATELY NOT FIXED
- Unmodelled nuisance, such as dissociation stress with no reference signature, can still pass certification.
- It is caught instead by the protocol-sensitivity modifier and by cross-protocol validation.
- The nuisance list is extensible, so a measured stress axis can be added as a decoy family.

### Validation plan

All primary success criteria are PRE-REGISTERED on OSF before the orthogonal validations (TCCA, cross-protocol, Tahoe) are run (month ~11), to pre-empt cherry-picking critiques.

A. SEMI-SYNTHETIC BENCHMARK (certmp-bench; evidential level A)

Backbone 1 (decoy-based):
- Decoys of 300 real tumours from CanSig (GBM, breast, LUAD, RMS, cSCC) plus a 3CA subset covering 6+ cancer types; 200-2,000 cells per tumour.

Backbone 2 (independent, guards against circularity):
- splatPop/Splatter [V] parameters estimated per real tumour.
- Confounders come from a SEPARATE injector module with different functional forms: logistic cell-level stress activation, cell-level CNA mosaicism, an alternative mean-zero relationship.

Planted programs (30 per cohort):
- 8 pan-shared (60-100% of tumours);
- 12 subset-shared (5-30% of tumours, lineage-restricted);
- 10 tumour-private (half CNA-linked);
- 4 nested, subclone-restricted;
- 3 between-tumour baseline programs.
Loadings are sparse (30-200 genes) and sampled from real 3CA MP and hallmark weight profiles, with pairwise Jaccard 0.05-0.3 to stress identifiability. Usage is active in 5-40% of cells, at log2FC 0.5-2.

Planted confounders:
- recurrent arm gains (chr7, 1q, 8q) in 40% of tumours, from real clone profiles;
- ambient alpha ~ U(0.02, 0.15);
- real-TME doublets at 2-10%;
- ex vivo stress (immediate-early / heat-shock signature) in 0-60% of cells, varying by study;
- 3x depth shifts between studies;
- an sc/sn platform layer;
- a QC axis.

Arms:
- (i) correct specification;
- (ii) MISSPECIFIED: stress absent from the decoy model; CNA estimated by expression callers while the truth is DNA-level; dispersion mis-estimated by a factor of 0.5-2;
- (iii) independent backbone;
- (iv) 50 NULL cohorts (nuisance only, zero planted programs).

Metrics:
- Hungarian-matched cosine; gene AUPRC;
- recall of subset-shared programs;
- false-MP rate (fraction aligned to injected nuisance);
- empirical FDR vs nominal q in {0.05, 0.1, 0.2};
- tumour-membership F1 and q calibration;
- gene-set empirical FDR and power vs top-50;
- typing macro-F1;
- leakage of private programs.

Targets:
- Null cohorts: P(at least 1 MP called) <= 0.1 at q = 0.1.
- Empirical FDR <= 1.25q under correct specification and <= 2q under misspecification (reported whatever the outcome).
- Subset-shared recall >= 0.8 at log2FC >= 1.
- Typing macro-F1 >= 0.85.
- Private leakage <= 5%.
- Gene-set FDR at nominal 0.05 <= 0.075.

B. BASELINES (9, tuned fairly with authors' defaults plus a documented grid)
- Gavish/3CA code itself;
- GeneNMF;
- cNMF per tumour plus Gavish-style clustering;
- TCAT-style pooled Harmony-cNMF;
- GBCD (per cancer type, subsampled);
- LIGER online iNMF;
- consensus scHPF or scDEF;
- DRVI;
- ITHmapper (scVI + Hotspot);
- CanSig-best early integration (Harmony + Leiden + DE).
Also reported: "CERT-wrapped" versions of Gavish and cNMF, to separate gains from the certification layer and from the engine.

Ablations:
- no decoys (held-out tests only);
- no thinning (in-sample coherence);
- depth-only decoys;
- no QC slopes;
- cell-level vs clone-level CNA;
- no private slots;
- centred instead of baseline;
- no genomic control.
Reserve, run only if reviewers ask: Tetris/BMSFA on a 30-tumour subsample (runtime and calibration), Spectra, MOFA-FLEX.

C. REAL-DATA ROBUSTNESS (3CA v1: 77 studies, 1,163 tumours [V]; 3CA v2 additions: 124 datasets, 2,836 samples [V])

Stability:
- Certified MPs across 3 thinning seeds and leave-one-study-out: median matched cosine >= 0.9 and coefficient of variation of the MP count <= 10%, against the Gavish-pipeline instability measured in the audit.

Cohort composition:
- With 50% random study subsampling (20 draws), at least 85% of certified MPs are re-certified (matched cosine >= 0.8). This is the analogue of the 3CA v1-to-v2 matching problem.

Replication:
- Discover on 3CA, replicate on TCCA (36 studies, about 1.09M malignant cells [V]) after removing overlapping studies.
- At least 80% of certified non-genetic MPs are re-certified on TCCA alone, vs at most 40% of technical or unsupported programs passing.

Artefact calibration:
- AUROC >= 0.85 of the type labels for the 11 MPs hand-removed in 3CA v2 (inds_rm [V+]) plus the LQ, doublet and single-study MPs Gavish 2023 removed, with zero manual input.

D. CROSS-PROTOCOL INVARIANCE (evidential level C; the most direct "real biology" test)
Datasets:
- Janesick et al. 2023 FFPE breast (Chromium Flex + Visium + Xenium) [V];
- Oliveira et al. Nature Genetics 2025 CRC (Flex + Visium HD + Xenium) [V];
- 3CA v2 sc vs sn pairs within cancer types [P].
Predictions:
- Certified non-genetic MPs expected in breast or CRC (from 3CA membership) are detected in dissociation-free Flex at rates comparable to fresh sc, with gene-weight Spearman >= 0.5 on detected genes.
- Protocol-sensitive and immediate-early "stress" programs show significantly lower detection in Flex and sn (CMH q < 0.05).
- Xenium is used only for MPs with panel retention iota >= 0.3.

E. PERTURBATION RESPONSIVENESS (evidential level B; replaces a wet lab)
- Tahoe-100M [V] at pseudocell or pseudobulk level per line x drug x dose vs DMSO.
- 10 drug-to-MP pairs are pre-registered, if the compounds are in the library [U]:
  - CDK4/6 inhibitor -> G1/S and G2/M down;
  - MEK inhibitor -> MAPK/immediate-early down;
  - proteasome inhibitor -> proteotoxic / UPR up;
  - HSP90 inhibitor -> heat-shock up;
  - mTOR inhibitor -> translation / MYC down;
  - DNA-damaging agents -> p53 / senescence up;
  - and others.
- Test: specific response at FDR 0.05 against DMSO and against matched random programs.
- Target: at least 8/10 pairs. Certified MPs also show a higher rate of any specific response than uncertified programs (Fisher test).
- Replogle Perturb-seq [P venue] is a backup.

F. GENETIC ORTHOGONALITY (minimal)
- At least 2 public tumour datasets with accessible reads, run with allele-aware Numbat [V]. Test whether typing (dosage vs non-genetic, clone-linked) agrees when using DNA/allele-informed rather than expression-only clones (Cohen's kappa >= 0.6).
- Positional arm-enrichment audit for all programs.

G. GENERALITY AND CLINICAL (light)
- Kinker CCLE 198 lines [V]: cell-intrinsic recurrence of certified MPs; technical types expected absent.
- TCGA: survival associations of certified vs technical MPs, each compared with 1,000 expression-matched random signatures (Venet 2011 [V]).
- Optional: certify the 3CA v2 non-malignant-cell MPs (8 cell types [P]) to show the framework is cell-type agnostic.

H. COMPUTE
- Runtime and memory vs cells (1e5 to 1e7) and tumours (100 to 3,000), against all baselines that finish.

### Flagship biology

"A CERTIFIED ATLAS OF TRANSCRIPTIONAL INTRATUMOUR HETEROGENEITY"
Re-derive malignant MPs across 3CA v1 plus v2 malignant data with CERT-MP and certify the published catalogues. All items below are PRE-SPECIFIED HYPOTHESES to test, not claimed results.

(1) THE STATUS OF EVERY CANONICAL MP
- For each of the 41 Gavish 2023 MPs [V] and the 3CA v2 malignant MPs, report certification q, type, modifiers, rho^full and rho^tech, and stability. Certification is computed on studies NOT used to derive each list (v2-only studies and TCCA), to avoid double dipping.
- Predictions:
  - cell cycle (G2/M, G1/S), interferon/MHC-II and hypoxia certify as non-genetic (positive controls);
  - LQ/ribosomal and contamination MPs type as technical (calibration against inds_rm);
  - at least one lineage MP is partly dosage-driven (e.g. in GBM, chr7 gain / chr10 loss, or 1q/8q in carcinomas). This depends on pilot P5.

(2) SPLITTING "STRESS"
- Hypothesis: the canonical stress MP mixes two components.
  - A protocol-sensitive immediate-early / heat-shock component (FOS/JUN/EGR1/HSPA1A-type), lower in Flex, sn and fixed data.
  - A protocol-invariant integrated-stress / hypoxia component that is spatially organised in Visium HD/Xenium CRC.
- This connects directly to GBCD's prognostic PDAC stress program: CERT-MP tests whether that program's gene set is protocol-invariant.
- This is the headline biology panel because it speaks to the user's complaint ("sometimes not real biology") with orthogonal, dissociation-free evidence.

(3) PROGRAMS THE THRESHOLD PIPELINE MISSES OR MERGES
- Subset-shared lineage programs supported by borrowing strength across small tumours (fewer than about 200 malignant cells), which per-tumour NMF at ranks 4-9 cannot resolve.
- Hybrid MPs (e.g. EMT-like variants) split into separately certified programs with distinct cancer-type memberships.
- Each new program must replicate on TCCA and pass cross-protocol checks where assessable.

(4) ITH AXES VS SUBTYPE AXES
- pi^ITH for every certified MP shows which programs vary mainly within tumours (plastic ITH) and which mainly between tumours within a study (subtype-like). The PDAC classical/basal and GBM state programs are tested against GBCD's subtype findings.
- Centring-based pipelines cannot produce this.

(5) DOWNSTREAM CONSEQUENCE
- Hypothesis: technical-typed MSigDB GAVISH_3CA gene sets show clinical associations in TCGA no better than matched random signatures, while certified MPs do better.
- This would be a concrete reason for users to migrate to certified gene tiers.

OPTIONAL EXTENSIONS (only if time allows)
- Certification of 3CA v2 MPs in non-malignant cells (generality).
- With a wet-lab collaborator: CRISPRi or drug perturbation of one top regulator of a newly certified program in 2 cell lines. This is NOT on the critical path.

### Software architecture

LANGUAGE AND ECOSYSTEM
- Python >= 3.10, scverse-native (AnnData/MuData in and out, scanpy-compatible), PyTorch backend (CUDA, with a CPU fallback).
- Minimal hard dependencies: anndata, numpy/scipy, pandas, torch, numba, statsmodels.
- Upstream tools (infercnvpy, SCEVAN, CopyKAT or Numbat clone calls; scDblFinder scores; cell-type labels) are accepted as PRECOMPUTED inputs in obs/obsm. This keeps installation light and avoids heavy R-in-Python.
- R/Seurat users get `certmpR` (reticulate plus zellkonverter): CertifyMetaprograms(seu, ...), ProjectMP(seu, ref = 'certmp-pan-v1'), and CertifyExternal(seu, gmt).

MODULES
- certmp.io: loaders and harmonisers for 3CA v1/v2, CanSig, TCCA, Kinker and Tahoe pseudobulk. Stores the UMI vs read-count flag and protocol metadata.
- certmp.pp: fit_nuisance(adata, tumour_key, study_key, protocol_key, clone_key, cna_key, tme=adata_tme, doublet_key, qc_keys); select_genes(method='deviance').
- certmp.decoy: make_decoys(adata, families=('full','tech'), seed); faithfulness(adata, decoys), which runs the nuisance-feature classifier AUC, mean / variance / zero-fraction QQ and CNA-block correlation.
- certmp.thin: thin(adata, eps=0.6, n_rep=3, dispersion='conditional').
- certmp.engine: CERTNMF(K_max=200, J=3, zeta='auto'), plus adapters for 'cnmf', 'gavish' (calls the tiroshlab R code), 'genenmf', 'gbcd' and any callable returning W.
- certmp.tl:
  - fit(adata_aug, engine=...);
  - coherence(adata, programs) for the M5 test with genomic control;
  - certify(adata_aug, q=0.1, alpha=0.1, min_studies=2);
  - type_programs();
  - gene_table(q_gene=0.05, three_way=False);
  - partition_ith();
  - certify_external(adata, programs);
  - project(query, ref, panel=None, novelty=True).
- certmp.pl: mp_card(), decoy_qq(), membership_heatmap(), positional_plot(), stability_plot().
- certmp.bench: simulate(backbone='decoy'|'splatpop', design, confounders, null=False); score(results, truth); a baseline registry with Apptainer/Docker containers per competitor.

WORKFLOWS
- Snakemake pipelines: `certmp-atlas` (3CA/TCCA end-to-end, reproduces every figure) and `certmp-bench`.

ENGINEERING
- Continuous integration with unit tests, including calibration tests: Z ~ N(0,1) on simulated nulls, and knockoff+ FDR on toy problems.
- Documentation on Read the Docs with 4 tutorials: certify my MPs, discover in my cohort, project my Flex/Xenium data, run the benchmark.
- Distributed on pip and bioconda. Versioned reference and benchmark data on Zenodo DOIs.
- A separate figure-reproduction repository, following the Spectra and TCAT deliverable pattern.

### Feasibility

TEAM
- Person A (100%): lead computational biologist / statistician, owns nuisance, decoy, thinning, certification and engine, PyTorch.
- Person B (50-100%): bioinformatician or PhD student, owns data ingestion (3CA, TCCA, CanSig, Kinker, Tahoe, 10x datasets), baseline containers, validations and package/docs.
- A wet-lab collaborator is NOT required. Their work is a bonus panel only.

COMPUTE
- 1-2 A100/H100 GPUs (about 1,500-3,000 GPU-hours total).
- About 50k CPU core-hours (baselines, nuisance fits).
- 3-5 TB storage.

DATA
- All public: 3CA v1/v2 (Weizmann portal), CanSig 12 datasets, TCCA (cnio-bu, Zenodo), Kinker CCLE, Tahoe-100M, 10x Janesick breast and Oliveira CRC, TCGA, and at least 2 tumour datasets with accessible reads for Numbat.

CUT LIST (cost or risk exceeds reviewer value for a small team)
- Full Bayesian IBP or spike-and-slab inclusion with MCMC or VI (Tetris, APAFA, SUFA). Replaced by held-out membership tests plus decoy knockoffs.
- Power posterior as a claimed contribution. It is kept only as an implementation detail.
- Normal cell-of-origin contrast (G6), because DeltaNMF and bcNMF reduce novelty and matched-normal data are sparse.
- Foundation-model priors.
- Wet-lab experiments on the critical path. Replaced by pre-registered Tahoe and Perturb-seq tests.
- KP-Tracer lineage heritability (mouse; separate pipeline).
- Selective-inference theory.
- Volume-regularised NMF. Kept only as a diagnostic plus Proposition 3.
- A web portal.
- Numbat on all of 3CA (BAMs unavailable). Limited to at least 2 datasets.
- More than 10 baselines. MOFA+, Spectra, expiMap, mosaicMPI and Tetris are held in reserve for reviewers.
- Exhaustive spatial work. One CRC Visium HD/Xenium dataset doubles as the cross-protocol and spatial panel.

SEQUENCING (the decoy generator is built FIRST because the audit, the method and the benchmark all reuse it)
- M0-1: Manual scoop check (Section 6.3 list). Ingest CanSig plus a 3CA v1 subset. Reproduce the Gavish code exactly.
- M1-4: Decoy generator v1 and faithfulness diagnostics. Pilots P1-P5. GO/NO-GO at M4.
- M4-8: Thinning, coherence tests, knockoff+ and typing. Label-blind decoy audit of the Gavish pipeline and GeneNMF on the full 3CA v1. Simulator v1. By M8 a complete, standalone AUDIT PAPER exists as the fallback ("How many pan-cancer metaprograms survive nuisance-matched decoys?").
- M6-11: CERT-NMF engine (GPU), replicate consensus, gene FDR, ITH partition. Baselines on the simulation grid.
- M10-14: Pan-cancer run (3CA v1+v2). OSF pre-registration at M11, then TCCA replication, cross-protocol, CCLE, artefact AUROC.
- M12-15: Tahoe tests, Numbat check, TCGA vs random signatures. Reference v1.0 and projector. Package v0.9.
- M14-15: Preprint plus Nature Methods submission.
- M15-24: Review and revision. Reserve analyses (Tetris subsample, MOFA-FLEX, Spectra, extra spatial) are run on request.

REALISTIC TOTAL: about 15 months to submission and up to 24 months including revision. That fits the 18-24-month envelope with about 3 months of slack.

### Risks and mitigations

R1. DECOYS UNFAITHFUL (missing or misspecified nuisance), so FDR is anti-conservative.
- Mitigations:
  - nuisance-feature two-sample classifier (target AUC <= 0.65);
  - QQ checks on gene mean, variance and zero fraction;
  - CNA-block correlation;
  - deliberate misspecification arms in the simulation, with inflation reported honestly;
  - random-program genomic control;
  - "empirical FDR" language;
  - an extensible decoy-family API so users can add measured artefacts.

R2. DECOYS OVER-FAITHFUL (nuisance absorbs biology, so power is lost).
- Examples: QC slopes absorb OXPHOS or translation programs; CNA offsets absorb co-regulated genomic gene clusters (MHC-II 6p21, histone 6p22, HOX, keratin 12q13/17q21); the baseline absorbs subtype axes.
- Mitigations:
  - clone-level rather than cell-level CNA;
  - ridge-shrunk QC slopes;
  - positive controls must certify (cell cycle, IFN, hypoxia);
  - a "QC-entangled" modifier instead of rejection;
  - gene-cluster sensitivity analysis;
  - the subtype axis is reported through pi^ITH.

R3. CNA CIRCULARITY (expression-based calls).
- Mitigations: clone-level offsets; comparison across infercnvpy, SCEVAN and CopyKAT; allele-aware Numbat check on at least 2 datasets (typing kappa >= 0.6); dosage programs are LABELLED, never deleted.

R4. THINNING VALIDITY (Hivert et al. critique).
- Mitigations: conditional dispersion estimation; cross-gene within-fold-B statistic (dependence is second order); genomic control; calibration shown on decoys and simulations with dispersion mis-estimated by 0.5-2x.

R5. THE STANDARD PIPELINE IS ALREADY ROBUST TO NUISANCE (pilot P2 negative).
- Pivot the headline to calibrated gene and tumour membership, typing and uncertainty-aware projection, with an adjusted venue (see fallbacks).

R6. THE ENGINE DOES NOT BEAT GBCD OR cNMF ON RECOVERY.
- The certification layer is engine-agnostic, so ship CERT-wrapped cNMF/GBCD as the default and keep the headline unchanged.

R7. HETEROGENEOUS PUBLIC DATA (Smart-seq2 read counts, processed-only matrices, missing TME cells).
- Mitigations: protocol flag with separate dispersions; read-count studies certified by projection rather than engine fitting; when a sample has no TME cells, use the same study's pooled TME for ambient and doublet profiles; at least 50 malignant cells required for testing.

R8. ARTEFACT LABELS REFLECT ONE LAB'S JUDGEMENT.
- AUROC on inds_rm is one of three artefact checks, alongside cross-protocol loss and decoy reproduction in simulation.

R9. REVIEWER ASKS "WHY NOT FULL BAYES?"
- Present the scaling argument, plus a reserve Tetris run on 30 tumours; held-out tests give frequentist guarantees without posterior coverage assumptions.

R10. SCOOP
- Threats:
  - the Stephens lab adds recurrence nulls to GBCD;
  - Kotliar/Raychaudhuri extend TCAT to malignant cells;
  - Abelson adds certification to ITHmapper;
  - Parmigiani/De Vito port multi-study NMF to scRNA-seq;
  - the Tirosh lab or the unattributed 2026 catalogue repo releases an updated MP resource.
- Mitigation: the novelty is the label-blind tumour-decoy certification plus typing, which is engine-agnostic and so survives model and resource scoops.
- Monitor every 2 weeks with Google Scholar alerts on citations of GBCD, TCAT, Gavish 2023 and Grabski 2025, plus bioRxiv keyword feeds ("metaprogram", "synthetic null program", "decoy", "multi-study NMF single-cell").
- Session scoop check (2026-10-02; 8 WebSearch queries; bioRxiv and arXiv fetch blocked, so weak evidence): no preprint combining nuisance-preserving null tumours or knockoffs with MP recurrence was found. Found only:
  - ZIMultiStudyNMF (https://arxiv.org/abs/2510.07518; diet data);
  - the known ITHmapper preprint (https://www.biorxiv.org/content/10.64898/2026.01.06.693992v1.full);
  - ClusterDE (https://pmc.ncbi.nlm.nih.gov/articles/PMC10401959/);
  - GBCD (https://www.nature.com/articles/s41588-024-01997-z).
- The manual check in Section 6.3 remains mandatory at M0.

PREPRINT TIMING
- Default: one preprint posted together with the Nature Methods submission (about M14-15). Nature Methods accepts preprinted work. Avoid splitting the audit off early, because the audit is the main paper's motivating figure, and an early idea disclosure lets larger labs bolt decoys onto GBCD or TCAT within about 6 months.
- TRIPWIRES that accelerate to an "audit + certification" preprint within 4-6 weeks (the content is ready from M8):
  - any preprint using synthetic or decoy nulls for gene programs or MPs;
  - any scRNA port of multi-study count NMF;
  - GBCD v2 with recurrence testing;
  - a new Tirosh-lab robustness analysis.
- Register the package name and keep the repository private until the preprint.

FALLBACK VENUES
- (a) Everything works: Nature Methods Article.
- (b) Certification works, engine only on par: still Nature Methods (the headline is certification), otherwise Nature Communications or Genome Biology (Method).
- (c) FDR calibration is weak under misspecification but the audit shows substantial nuisance-driven recurrence: Nature Methods Analysis, Genome Biology (Research/Analysis) or Cell Systems, as audit plus benchmark plus simulator.
- (d) Audit negative (MPs robust): calibrated membership, typing and projection go to Genome Biology software or NAR Genomics & Bioinformatics, with a Bioinformatics application note for the package; the certified catalogue could go to Nature Cancer as a resource if the biology (stress split, dosage MPs) holds.

### Go/no-go pilots

All pilots use public data, about 10-12 person-weeks in total, and finish by M4. The first two decide the project.

P0 (1 week). Manual scoop check of the Section 6.3 list (bioRxiv, Scholar, citing papers of GBCD, TCAT, Gavish 2023 and Grabski 2025).
- NO-GO / REDESIGN if a preprint already certifies MP recurrence against nuisance-preserving null tumours.

P1 (3 weeks). DECOY FAITHFULNESS on 200-300 tumours (CanSig 12 datasets plus a 3CA v1 subset).
- GO if all of the following hold:
  - nuisance-feature classifier AUC <= 0.65 in at least 80% of tumours (features: log depth, %mito, %ribo, detected genes, TME-marker score, clone);
  - gene mean, variance and zero-fraction QQ within ±10%;
  - CNA segment correlation r >= 0.8;
  - positive controls (G2/M, G1/S, IFN where active) have coherence Z >= 5 in real tumours and |Z| < 2 in decoys for at least 90% of tumours.
- If faithfulness fails: add the missing nuisance terms (2 more weeks), then re-test.

P2 (2-3 weeks; THE KEY DECISION). LABEL-BLIND DECOY AUDIT OF THE STANDARD. Run the tiroshlab/3ca code unchanged on real plus D^full tumours mixed together.
- GO if either:
  - at least 2 MPs have 30% or more of their member programs coming from decoys; or
  - at least 10% of the robust programs that pass 35/10/10 come from decoys.
- NO-GO / PIVOT if fewer than 3% of programs are decoy-derived and no MP is decoy-majority. The standard is then nuisance-robust, so pivot to calibrated membership, typing and projection (fallback d).
- This pilot is also Figure 1 of the paper and the core of the fallback audit paper.

P3 (1-2 weeks). CALIBRATION OF THE HELD-OUT COHERENCE TEST.
- On decoys, with programs learned on another tumour's fold A: Z vs N(0,1) after genomic control should have KS p > 0.01 and type-I error at 0.05 within [0.03, 0.07].
- In real tumours, the genomic-control median lambda_t should be <= 1.3; otherwise revise the nuisance mean model (depth or size terms).
- Repeat with dispersion deliberately mis-estimated by a factor of 0.5-2. Type-I error must stay <= 0.08; otherwise move to three-way thinning or competitive tests.

P4 (1 week). ARTEFACT PREDICTABILITY.
- On a 3CA v2 run with the code's own outputs, compute rho^tech and W^full for the 11 inds_rm MPs vs retained MPs.
- GO if AUROC >= 0.8. Below 0.65, the nuisance model misses the dominant artefacts, so add families (e.g. a measured stress axis from sc/sn pairs) or demote typing.

P5 (1-2 weeks). CNA / POSITIONAL AUDIT.
- Arm-level hypergeometric enrichment (BH 0.05) of robust programs and MPs, plus clone eta^2 of usage.
- If at least 5% of MPs are affected, dosage typing is a headline axis. Otherwise keep it as a reported modifier and do not headline it (as the synthesis recommends, given Kinker's 8% figure for continuous programs).

P6 (about 4 weeks, M6-10; NOT gating the headline). ENGINE SMOKE TEST.
- CERT-NMF should match or beat GBCD and cNMF+clustering on subset-shared recall in the backbone-1 simulation (F1 at least as high as each) and run 3CA v1 in under 48 GPU-hours.
- If not by M10, the default engine becomes CERT-wrapped cNMF, and the paper headline is unchanged.

DECISION MATRIX AT M4
- P1 and P2 pass: full plan.
- P1 passes, P2 fails: pivot to fallback d.
- P1 fails after one revision: publish the audit and benchmark (fallback c) and stop method development.

### Why Nature Methods

(1) A conceptual advance stated as a definition, which is exactly what the synthesis (Section 5) identifies as the bar.
- An MP becomes a statistical object: held-out within-tumour co-variation, recurrent across studies, in excess of nuisance-matched decoy tumours at a stated FDR.
- The number of MPs, tumour membership and gene membership all become calibrated outputs rather than threshold by-products.
- No published or preprinted method found by the sweeps or this session's checks controls error on program EXISTENCE or RECURRENCE, or models CNA and technical nuisance inside discovery (synthesis Section 6.2: gaps ii and iii are uncovered).

(2) A general principle, not a tumour-specific trick.
- Label-blind tumour-level decoys with knockoff+ selection, plus cross-gene held-out coherence tests, apply to any recurrent-program discovery across samples: TME programs, T-cell catalogues such as TCAT, perturbation programs.
- certify_external makes it a layer on top of existing tools, which Nature Methods editors value for broad uptake. The optional non-malignant-cell panel demonstrates this.

(3) Ground truth that includes the failure modes.
- A semi-synthetic simulator with planted confounders, misspecification arms and nuisance-only null cohorts, plus 9 fairly tuned baselines including the mandatory Gavish code, GeneNMF, TCAT-style cNMF and GBCD.

(4) Several orthogonal validations, achievable without a wet lab.
- Pre-registered cross-protocol invariance on dissociation-free Flex.
- Pre-registered perturbation specificity (Tahoe-100M).
- Independent-cohort replication (3CA to TCCA).
- External artefact labels (the 11 hand-removed 3CA v2 MPs).
- An allele-aware clone check.

(5) Atlas-scale application with biology the community will notice.
- A certified version of the most-used malignant-program catalogue (3CA / MSigDB GAVISH_3CA): which MPs are real, dosage-driven or technical; a stress program split by protocol invariance; programs that were missed or merged; and an ITH-vs-subtype partition.
- Shipped as a versioned reference with an uncertainty-aware projector and a reusable benchmark.

(6) Strategic fit for a small team.
- The costly, scoopable parts (full Bayesian model, normal contrast, foundation-model priors, wet lab) are cut.
- The remaining core is cheap, engine-agnostic and therefore robust to competitors' model or resource releases.
- Sequencing makes the M8 audit a publishable fallback, keeping the downside at Genome Biology or Cell Systems level while preserving the Nature Methods upside.

### Differentiation from closest competitors

| Competitor | Differentiator |
|---|---|
| GBCD (Liu…Stephens, Nature Genetics 2025) [V] | GBCD is a Gaussian model of the cell-by-cell covariance of pooled log data. Kmax is user-set, it runs per cancer type, has not been shown at atlas scale, offers no test of program existence or recurrence, no tumour-membership calls, no CNA or technical modelling and no projector. CERT-MP's novelty is NOT the shared-vs-private split, which GBCD owns. It is the certification layer: label-blind decoy tumours, a knockoff+ FDR for MP existence, held-out membership and gene q-values, and typing by decoy family. GBCD is run INSIDE CERT-MP as one engine to show that certification also improves GBCD output. The CERT-NMF engine adds count likelihood, CNA, QC and ambient offsets, a pan-cancer joint dictionary and an ITH/baseline partition. |
| ITHmapper / Pan-Cancer ITH Compendium (Nofech-Mozes…Abelson, bioRxiv 2026, doi 10.64898/2026.01.06.693992) [V] | ITHmapper uses scVI with sample as batch plus Hotspot modules and gives discrete states. CNA is used only for malignant calling, and there is no FDR or artefact typing. Sample-as-batch can erase programs active across a whole tumour. CERT-MP gives continuous overlapping programs, keeps and decomposes tumour baselines (pi^ITH), and certifies and types every program. ITHmapper modules are benchmarked as a competitor AND certified with certify_external, a direct, non-adversarial head-to-head. |
| TCAT/starCAT (Kotliar…Raychaudhuri, Nature Methods 2025) [V] | TCAT covers T cells only and uses Harmony-corrected cNMF. Singleton programs are excluded by heuristic plus judgement, and NNLS projection has no uncertainty; its authors warn about Flex transfer. CERT-MP targets malignant cells. Existence is decoy-certified rather than judged, and projection carries per-cell SE, sample-level q-values, novelty detection and panel-retention flags. The reference is a deliverable, explicitly NOT the headline. The certification layer is cell-type agnostic and could certify TCAT programs (optional generality panel on 3CA v2 non-malignant MPs). |
| Bayesian multi-study NMF (Grabski, Trippa, Parmigiani, Genome Biology 2025) [V]; Tetris (AoAS 2023) [V]; ZIMultiStudyNMF (Hansen, Zorzetto, Edefonti, De Vito, arXiv 2510.07518, Oct 2025; diet data) [V this session] | These are fully Bayesian MCMC models with inclusion matrices for a few hundred features (mutational channels, diet items). None has been applied to tumour scRNA-seq or uses nuisance-matched nulls. CERT-MP does not depend on a correct posterior: inclusion (tumour membership) comes from held-out tests and existence from decoy contrasts, so it scales to millions of cells. If these groups port their model to scRNA-seq, CERT-MP certifies its output too, so the headline survives a model scoop. |
| ClusterDE synthetic null (Song, Li et al.) [V] and recall knockoff genes (DenAdel…Crawford, AJHG 2025) [V] | Both calibrate clustering or differential expression within ONE dataset. ClusterDE's null preserves gene correlation and so can keep real programs; recall's marginal knockoff genes lack correlation. CERT-MP's null unit is a TUMOUR. It preserves declared nuisance (clone CNA, ambient, real-TME doublets, QC slopes, depth) but no programs, is fitted label-blind next to real tumours to get knockoff flip-sign symmetry, and tests cross-study RECURRENCE rather than DE. Decoy FAMILIES additionally give mechanistic typing (technical vs dosage). |
| Gavish 2023 / 3CA v1-v2 pipeline (Nature 2023; Nature Cancer 2025) [V+] and GeneNMF / scalop / mpnmf / sciNMF [V] | These are threshold stacks with greedy or forced clustering and manual curation. CERT-MP does not compete on speed or determinism, a saturated niche. It wraps these pipelines (label-blind on fold A of real plus decoys) and reports, for each MP, the fraction of its recurrence that decoys reproduce. The audit of the standard is therefore a built-in result, and the 11 MPs hand-removed in 3CA v2 serve as an external calibration set. |
| CanSig benchmark (Boeva lab, Cancer Research 2026) [V] | CanSig scores rediscovery of known signatures and integration metrics. It plants no confounders, runs no null cohorts and studies no seed, rank or thinning stability. CERT-MP reuses its 12 datasets and metrics for continuity, and adds planted CNA, ambient, doublet, QC and stress confounders, nuisance-only null cohorts, FDR calibration, typing accuracy and leakage of tumour-private programs. |
| scDEF [V; peer-review status unknown], consensus scHPF (bioRxiv 2023) [V], LIGER iNMF / scINSIGHT [V] | These are count-native or shared-vs-dataset-specific factorisations. None offers error control on program existence, nuisance or CNA typing, or held-out gene FDR. They are used as baseline engines; one count-native Bayesian model (consensus scHPF or scDEF) plus LIGER iNMF are benchmarked. |
| DRVI / MOFA-FLEX / Spectra / DeltaNMF / bcNMF [V] | These are generic disentangling, prior-guided, contrastive or foundation-model-regularised factor models with no recurrence FDR or tumour nuisance model. DRVI is benchmarked as the representative deep model. Foundation-model priors and normal-contrast angles are deliberately NOT pitched (DeltaNMF and bcNMF exist). |
| BayesPowerNMF (PLoS Comput Biol 2026) [V], bayesNMF (JCGS 2026) [V] | These are robust or rank-learning Bayesian NMF for mutational signatures. CERT-MP borrows the tempering idea only as an out-of-sample-tuned regulariser (equivalent to prior rescaling under MAP) and cites them as foundations. No novelty is claimed there. |

### Referee 1: novelty and prior art

**Predicted editorial outcome:** Desk reject or review, then reject or major revision. My best estimate is about a 35-45% chance of desk rejection. An editor will see three things: a statistical wrapper for a pipeline used by one subfield; a 'general principle' (run the same pipeline on target data and on synthetic-null data, take a contrast score, apply a Barber-Candes threshold) that is already published as Clipper (Genome Biology 2021) and ClusterDE, and a model-based null for any algorithm's output already published as sc-SHC (Nature Methods 2023); and an engine that the authors themselves say is not the contribution. The paper gets sent out only if the P2 audit figure is striking, e.g. a large fraction of a named published 3CA MP is reproduced by decoys through the unchanged Gavish code. If it is reviewed, I expect at least one referee to make the circularity argument below: decoys are drawn from the same parametric nuisance model that is conditioned on, so the FDR guarantee is empty exactly where it matters. That leads to reject, or major revision demanding non-parametric nulls and protocol-orthogonal certification. The realistic landing is Genome Biology, Nature Communications or Cell Systems, or a Nature Methods Analysis built around the audit plus benchmark. A Nature Methods Article becomes plausible only if the null construction is redesigned so that certification is demonstrably sensitive to nuisance the model does not include.

**Fatal flaws**

- THE DECOY NULL IS CIRCULAR. Decoys are NB draws from mu0, the same fitted nuisance model (M1) that the engine uses as fixed offsets (M4) and the coherence test conditions on (M5). Case 1, correct specification: decoys are guaranteed nulls under the test, so decoy Z ~ N(0,1) and R~ is just the BH false-positive floor. The decoy arm adds almost nothing beyond 'held-out residual coherence in >=2 studies'. The P1 GO criterion 'positive controls |Z|<2 in decoys' is satisfied by construction and tells you nothing. Case 2, misspecification, the only case relevant to 'is it real biology': every systematic misfit exists only in real tumours. Examples are within-clone CNA mosaicism, gene-specific dosage compensation, cell-size scaling not captured by library size, non-NB tails, Smart-seq2 zero patterns and unmodelled dissociation stress. Each gives W = R - R~ > 0 and is certified as a 'non-genetic MP'. The decoy contrast is blind in exactly the regime it is meant to police.
- THINNING ASYMMETRY. Decoys have exactly NB(mu0, r) counts, so NB data thinning is exact for them, while for real tumours it is approximate (Hivert et al., arXiv 2405.13591 [V]). Any fold dependence from mis-estimated dispersion therefore inflates coherence in real tumours only, and the knockoff contrast reads it as signal. The concept presents decoys plus genomic control as its answer to the Hivert critique; structurally, decoys cannot answer it.
- PROP 2 (FDR) DOES NOT FOLLOW EVEN UNDER ITS OWN ASSUMPTION. Knockoff+ / target-decoy competition needs a flip-sign property for each null hypothesis: swapping a null feature with its knockoff must flip that W_k and leave the joint law unchanged. Swapping a real tumour with its decoy flips every W_k at once, non-null programs included, and is not a symmetry of the data, because real tumours contain biology and decoys do not. One decoy per tumour also serves all ~200 programs, so null signs are dependent. Dictionary rows mix biology and nuisance, so 'null program' is not even defined. The method offers no valid FDR statement. 'Empirical FDR' rests on simulations whose main backbone (backbone 1) is produced by the same decoy generator, which is circular.
- THE TYPING LOGIC CONTRADICTS THE CONDITIONING SET. N_t includes clone CNA (delta, kappa), and activity A(w) is defined by a test conditional on N_t. CNA-driven co-variation is therefore residualised out of R in real tumours, and 'dosage programs' (fail W^full, pass W^tech) can only arise from CNA misfit. The concept does not say which mu0 is used to score D^tech decoys. If it is the real tumour's kappa*delta, D^tech residuals carry anti-dosage structure; if it is kappa = 0, the arms are scored under different nulls. Either way the dosage/non-genetic typing has no coherent definition. rho = R~/R is also undefined when R ~ 0, which is precisely the case for LQ/ribo/mito programs: those are conditioned away by the %mito/%ribo covariates. Yet these programs make up 7 of the 11 inds_rm MPs used for the headline AUROC.

**Major concerns**

- THE HEADLINE OVERCLAIMS. 'Show, without manual curation, which pan-cancer metaprograms are real' really means 'not explained by six declared nuisance terms under a per-tumour NB GLM'. The user's main complaint, dissociation stress, is explicitly NOT in N. So the biology the user cares about most is exactly what certification cannot adjudicate. The certified label must be renamed and reframed as relative to N.
- CONDITIONING ON POST-TREATMENT COVARIATES. Library size, complexity residual, %mito and %ribo are all downstream of biology: cell cycle drives cell size, UMI count and complexity; OXPHOS, hypoxia, apoptosis and stress move %mito; translation and MYC move %ribo. Gene-specific ridge slopes gamma_tg absorb part of these programs into mu0. Decoys then CONTAIN those programs, so R~ > 0 for real biology, and real programs lose certification or are mislabelled 'QC-entangled'. Ridge/GCV does not solve this, because it is a causal identifiability problem, not a variance problem. The cell-cycle positive control is at real risk.
- THE CNA MODEL IS TOO CRUDE TO SUPPORT DOSAGE TYPING. A single per-tumour kappa_t multiplying clone-level delta ignores gene-specific dosage sensitivity; CLCC (bioRxiv 2026, 10.64898/2026.03.01.708901 [V by sweeps]) is built on that variability. It also ignores within-clone mosaicism. Expression-based callers also smooth co-regulated positional clusters into fake CNAs (MHC-II 6p21, histone, keratin, HOX). Residual dosage then leaks into 'non-genetic' programs, and co-regulated clusters leak into the offset. Numbat on 'at least 2 datasets' is far too thin to establish typing accuracy.
- THE 'SECOND-ORDER' ROBUSTNESS CLAIM FOR M5 IS HAND-WAVING. If thinning leaves entrywise A-B correlation rho, then E[r^B_g r^B_g' | A] scales roughly as rho^2 times the fold-A co-fluctuation. The engine selected w to maximise exactly that co-fluctuation, summed over thousands of cells. The resulting bias can be many null SDs. Formal bounds or simulation at realistic n_cells, with dispersion mis-estimated, are required. A rank-10 NMF fitted on the full data, used to condition dispersions, is itself a double-dipping step.
- GENOMIC CONTROL WITH MEAN-MATCHED RANDOM PROGRAMS IS THE WRONG NULL. It corrects global inflation but not structured, gene-set-specific confounding: positional blocks, highly expressed blocks, size-scaled genes. It also deflates power in tumours dominated by strong real programs.
- POWER AND DEPTH BIAS IN MEMBERSHIP. With eps = 0.6, fold B keeps 40% of the counts. Tumours with 50-200 cells, and shallow studies, will be systematically under-called. R is then confounded with study depth and protocol. This also contradicts flagship claim (3), that borrowing strength rescues small tumours: discovery may borrow, but certification and membership do not.
- TARGET-DECOY MECHANICS. With one decoy per tumour, competition-based FDR is known to have high FDP variance; the proteomics literature developed averaging and multiple-decoy schemes for this. Knockoff+ at q = 0.05 needs at least 20 selections with no negatives before it can select anything. Averaging over 3 thinning replicates does not address decoy-draw variance. The study-capped R saturates at 2 tumours per study and makes no distinction between protocol diversity and study diversity: 'at least 2 studies' can mean 2 fresh-dissociated 10x studies.
- THE ENGINE IS NOT NOVEL AND ITS IDENTIFIABILITY IS UNADDRESSED. CERT-NMF is LIGER/online-iNMF structure (shared W, private V_t) with a Poisson likelihood (scHPF/fastTopics), multiplicative offsets (sciRED/GLM-PCA style) and Tan-Fevotte ARD. With 4,500 tumour units x J = 3, there are 13,500 private columns next to K = 200 shared ones. The real identifiability problem is shared-versus-private ambiguity (SUFA's warning), and the regulariser decides it. Prop 3 (adding rows enlarges cone(H)) is a one-line monotonicity remark about the shared block with private components fixed. It does not address that ambiguity and will not impress a theory referee.
- THE VALIDATION SUITE IS WEAKER THAN ADVERTISED. (a) Cross-protocol: Janesick 2023 is a single FFPE breast block from one patient [V]; Oliveira 2025 CRC is n = 5 [V]. Neither has fresh-dissociated sc from the same tumours, so 'Flex detection comparable to fresh sc' becomes a between-patient comparison, confounded by patient. (b) Tahoe-100M: CDK4/6 inhibitor reducing cell-cycle genes and HSP90 inhibitor inducing heat-shock genes would hold for any hallmark gene set, so they do not validate certification. ITHmapper already uses Tahoe-100M [V per synthesis]. The certified-versus-uncertified Fisher test is confounded by program class: contamination and LQ programs cannot respond in TME-free lines. (c) The inds_rm AUROC rests on 11 positives with a wide CI. A trivial covariate score (mito/ribo fraction, TME-marker overlap) is the obvious baseline and probably matches it; that baseline is not planned. (d) Backbone 1 simulations come from the method's own generator, which is circular.
- THE FLAGSHIP BIOLOGY IS PARTLY PRE-EMPTED. The Tirosh lab already separates baseline programs from within-tumour states in GBM ('The multilayered transcriptional architecture of glioblastoma ecosystems', Nature Genetics 57:1155-1167, 2025 [V title/venue]). That undercuts M9 and the claim that 'centring-based pipelines cannot produce this'. sc-versus-sn depletion of stress states and dissociation-induced IEG signatures are known: 3CA v2 compared sc and sn [P], and the dissociation literature exists [U]. The 'stress split' is a plausible but incremental result.
- FEASIBILITY. 1.5-2 FTE over 15 months is asked to deliver: nuisance model, decoy generator, thinning, GPU engine, certification, gene-level FDR, ITH partition, projector with Fisher SEs and panel transfer, an R bridge, a two-backbone simulator with misspecification arms and 50 null cohorts, 9-10 containerised baselines (GBCD per cancer type, ITHmapper scVI+Hotspot, about 7,000 Gavish NMF jobs), ingestion of 3CA v1+v2, TCCA, Tahoe, Flex/Xenium, Numbat with BAMs, CCLE and TCGA, plus OSF pre-registration. That is at least 3-4 FTE-years. The nuisance model also needs raw UMI counts, same-sample TME cells and clone calls. A significant fraction of 3CA (Smart-seq2 TPM, processed-only matrices) does not have these, and how much is not quantified [U].
- GBCD COMPATIBILITY. GBCD signatures are signed: GBCD does not constrain programs to be non-negative [V, GBCD abstract and press]. M5 and M8 assume w >= 0, so 'GBCD run inside CERT-MP' needs an explicit mapping, e.g. truncation. That changes the GBCD output being 'certified'.

**Prior-art threats**

- Clipper (Ge, Chen, Li et al., Genome Biology 2021) [V]. p-value-free FDR via target-versus-background contrast scores thresholded Barber-Candes style; demonstrated across ChIP-seq, scRNA-seq DE, Hi-C and proteomics. This is the general principle M6 claims, and the concept does not cite it.
- ClusterDE (Song, Li et al.; preprint PMC10401959 [V]; software paper PMC13309996 [V]; journal venue [P]). Synthetic null (scDesign3) runs through the identical pipeline; contrast = target score - null score; Clipper threshold. CERT-MP's M6 changes the unit (gene to program) and the null (one cell type to a nuisance-only tumour). That is a domain transfer, not a new principle. ClusterDE also deliberately preserves gene-gene correlation, i.e. it takes the conservative side of the trade-off that CERT-MP resolves in the anti-conservative direction.
- sc-SHC (Grabski, Street, Irizarry, Nature Methods 20:1196-1202, 2023) [V]. Model-based null fitted to the data; significance of clusters from ANY algorithm; batch-aware. This is direct conceptual precedent for certify_external and for 'calibrated existence' claims, published in the target journal.
- CHOIR (Nature Genetics 2025, s41588-025-02148-8) [V title/venue; details not checked]. Significance-based detection of cell types and states.
- recall (DenAdel...Crawford, AJHG 2025) [V]. Knockoff artificial variables to calibrate over-clustering.
- Count splitting for inference after latent-variable estimation (Neufeld et al., arXiv 2207.00554 / Biostatistics 2024 [V arXiv]) and NB count splitting (arXiv 2307.12985) [V]. M3 and M8 (held-out gene membership for estimated programs) apply these directly.
- Spectra (Kunes...Pe'er, Nature Biotechnology 2023/24) [V]. Held-out co-occurrence ('coherence') of factor genes already used as the program-validity criterion. Precedent for M5's held-out coherence idea. CS-CORE (Nature Communications 2023, s41467-023-40503-7) [V title] is precedent for co-expression testing under a measurement model.
- CSI-GEP (Cell Genomics 2025) [V]. Program and K selection by replication across data splits. Precedent for replicate-stability selection.
- LIGER iNMF / online iNMF [V], scHPF/consensus scHPF [V], Tan & Fevotte ARD [U]. Together they make up CERT-NMF; no engine novelty.
- sciRED (Nature Communications 2025) [V]. Poisson-GLM removal of known covariates, then factors matched to covariates, with 'unexplained' factors flagged. Close precedent for nuisance-GLM-then-type. KINOMO (bioRxiv 2022, AACR 2024) [V] is KL-NMF with correction for genomic variability between patients, precedent on the CNA angle.
- GBCD (Liu...Stephens, Nature Genetics 57:263-273, 2025) [V]. Owns shared versus patient-specific; signed programs.
- Tirosh lab, 'The multilayered transcriptional architecture of glioblastoma ecosystems' (Nature Genetics 57:1155-1167, 2025) [V title/venue via search; method details unverified]. Baseline programs separated from within-tumour states, which undercuts the M9 / flagship (4) novelty.
- ITHmapper (Nofech-Mozes...Abelson, bioRxiv 2026) [V]. Already validates with Tahoe-100M, Visium HD and Kinker lines, and ships a query mapper. It overlaps CERT-MP's validation and projection story.
- starCAT reference catalogue already extended beyond T cells (MYELOID.GLIOMA.V1, BONEMARROW.CD34POS.HSPC.V1) [V, starCAT GitHub/website]. The TCAT group has the pipeline to release a malignant reference.
- CLCC, 'Contrastive Alignment of Expression and Copy Number Highlights Dosage-Insensitive Genes in Cancer' (bioRxiv 2026) [V by sweeps]. Gene-specific dosage sensitivity, which challenges the single-kappa CNA offset.
- Target-decoy competition in proteomics: Elias & Gygi 2007 [U, not verified this session]; competition-FDR variance and multiple-decoy literature, e.g. arXiv 1907.01458 [V]. The conceptual ancestor of tumour 'decoys'; its known variance problems apply here.
- scPFA, lineage-informed phylogenetic factor analysis of heritable programs (Staklinski et al., Sept 2026 preprint) [V via ResearchGate snippet only; venue unknown]. Peripheral, but competes on heritable versus plastic program typing.
- Session scoop check (about 20 WebSearch queries; bioRxiv and nature.com fetches blocked) found no preprint certifying MP recurrence against nuisance-preserving null tumours. That is weak evidence of absence.

**Required changes**

- Redesign the null so it is NOT drawn from the conditioning model. Options: (i) non-parametric within-tumour decoys, e.g. gene-wise conditional permutation of cells within nuisance strata (clone x depth x QC bins), which keep each gene's real marginal misfit; (ii) real negative-control data as decoy families, e.g. matched sn/Flex/fixed samples for stress, TME-free cell lines; (iii) score decoys and real tumours with nuisance models re-estimated identically on each, and show that the decoy arm changes conclusions relative to the 'no decoys' ablation. If the no-decoy ablation performs the same, drop the knockoff layer.
- Validate decoy faithfulness on the statistic that matters: the residual gene-gene covariance spectrum of known-nuisance gene sets (ribo, mito, IEG/HSP, ambient/TME markers, positional blocks) in real versus decoy tumours. Marginal-feature classifier AUC and QQ plots are not enough.
- Replace Prop 2 with a correct statement. Either define per-program null hypotheses with a per-hypothesis flip-sign construction, or drop 'knockoff+ FDR control' and report a calibrated empirical contrast. Use multiple decoy draws per tumour (multiple-competition / averaged TDC) and report FDP variability.
- Specify the conditioning set and mu0 for each arm (real, D^full, D^tech). Re-derive typing so that the recurrence statistic used to type a nuisance does not condition on that same nuisance. Define behaviour when R ~ 0.
- Give a formal analysis, or a dense simulation, of M5 bias under fold dependence (dispersion mis-estimated 0.5-2x) at realistic cell numbers, with engine-selected w rather than fixed w. Add real-data negative controls: programs learned in unrelated cancer types or cell lines, tested on tumours; technical replicates.
- Make the independent simulation backbone (splatPop plus a separate confounder injector) the PRIMARY evaluation. Include arms with unmodelled nuisance, cell-level mosaic CNA and gene-specific dosage sensitivity, where the method is expected to fail, and report failures honestly.
- Add the missing competitors and baselines: Clipper/ClusterDE-style synthetic-null certification (scDesign3 with copula); sc-SHC-style model-based significance adapted to programs; 'held-out coherence + >=2 studies' with no decoys; a trivial covariate-score typer (mito/ribo/TME-marker overlap/arm enrichment) for the inds_rm AUROC.
- Cite and differentiate Clipper, ClusterDE, sc-SHC, CHOIR, recall, count-splitting LVM inference, sciRED, CS-CORE, the Tirosh 2025 GBM multilayer paper and CLCC. Drop 'general principle' language unless something beyond a domain transfer is shown.
- Make protocol-orthogonal evidence a primary certification endpoint rather than a validation: matched sc/sn/Flex designs from the same tumours, e.g. 3CA v2 sc-sn pairs and GBM sc-vs-sn cohorts. Acknowledge that Janesick (1 patient) and Oliveira (n = 5) cannot support invariance claims alone. Pre-register the stress split with a power calculation.
- Reframe the claim as 'certified relative to declared nuisance N' throughout, including the package and MP-card labels. Remove 'show which MPs are real'.
- Cut scope. Ship certification and audit on cNMF/Gavish outputs first and make CERT-NMF optional. Move projector SEs and panel transfer to supplementary. Reduce baselines to the mandatory set (Gavish code, GeneNMF, cNMF, GBCD, ITHmapper, one count-native Bayesian model). Before M1, audit how many 3CA/TCCA tumours have raw UMIs, same-sample TME cells and usable clone calls.
- Test the post-treatment-covariate problem directly. Plant programs correlated with depth, complexity, %mito and %ribo (cell cycle, OXPHOS, translation) and report loss of certification. Offer a sensitivity mode without QC slopes as a co-primary result, not an ablation.

**Ideas worth keeping**

- The label-blind decoy AUDIT of the unchanged Gavish/GeneNMF pipeline (P2), i.e. what fraction of each published MP's member programs comes from nuisance-only tumours. Here decoys are meaningful because the audited pipeline does not condition on nuisance. It is falsifiable, cheap, figure-ready and publishable as a stand-alone analysis.
- Per-tumour held-out (thinned) coherence q-values as an engine-agnostic replacement for threshold-based membership, with 'insufficient cells' reported separately from negatives.
- FDR-controlled held-out gene membership with effect size and I^2 (core versus context genes), replacing top-50 truncation. Exportable as confidence-tiered GMTs.
- Typed outputs that label rather than delete (technical / dosage / non-genetic), plus continuous modifiers. The protocol-sensitivity CMH test stratified by cancer type is simple and useful.
- Study-capped recurrence as a guard against single-cohort domination; leave-one-study-out and 50% study-subsampling stability as reported quantities.
- Using the 11 hand-removed 3CA v2 MPs plus Gavish's removed MPs as an external labelled artefact set, with a trivial baseline and CIs.
- Panel-retention index iota_k to flag MPs not assessable on Flex/Xenium panels.
- Project hygiene: pre-registration before orthogonal validations, GO/NO-GO pilots, and a fallback audit paper at M8. This is sound strategy for a small team.

### Referee 2: statistical rigor and feasibility

**Predicted editorial outcome:** If this were submitted as designed, with results that matched the plan, I expect it would go out for review: the topic is high-visibility, it works at atlas scale, it audits a standard pipeline, and it ships a tool plus a benchmark. From a statistical referee I would then expect major revision trending to reject. The certification layer is the claimed contribution, and as written it has three problems. (1) It goes vacuous for its own flagship engine: decoys are drawn from the same parametric nuisance model that the engine and the test condition on. (2) Its membership test is built on a sharp null that saturates once tumours have thousands of cells. (3) Its FDR proposition does not follow from the stated exchangeability assumption. On top of that, the "real biology" headline rests almost entirely on self-consistency (same decoy model, same nuisance model, replication with the same pipeline) and on well-known positive controls. The only dissociation-free evidence comes from 1 FFPE breast block (Janesick 2023) and 5 CRC samples (Oliveira 2025). Redesigned along the required changes, it is a credible Nature Methods paper. The best version decouples decoys from engine offsets, uses an effect-size null with cell-splitting tests, has independent decoy replicates with demonstrated flip-sign behaviour, and adds a real dissociation-free cohort. Without that redesign, fallback (c) is the realistic outcome: the M8 decoy audit of the Gavish pipeline published as a Nature Methods Analysis, Genome Biology or Cell Systems paper.

**Fatal flaws**

- F1. THE DECOY CONTRAST GOES VACUOUS FOR THE CERT-NMF ENGINE (internal inconsistency).
- Decoys are sampled from mu0 (M1/M2).
- The engine (M4) conditions on exactly the same Omega_tcg (kappa*delta + gamma'q) and alpha_tc*rho_tg as fixed offsets.
- The test (M5) computes residuals under the same N_t.
- So a decoy is explained by the offsets up to NB noise. The only exception is the doublet/misclassification injection, which the engine does not model.
- Consequences:
  - R~^full ≈ R~^tech ≈ 0 for every CNA-, QC- and ambient-driven structure.
  - BH within the decoy arm rejects almost nothing.
  - Every W_k ≈ R_k ≥ 0. With no negative W_k, knockoff+ selects every program with R>0 once at least 1/q of them exist (10 at q = 0.1).
- Certification therefore collapses to 'held-out membership in tumours from ≥2 studies'. That is the recurrence criterion the paper says is necessary but not sufficient.
- The DOSAGE label (fails W^full, passes W^tech) cannot fire, because CNA dosage is offset in both the engine and the test, so D^full decoys cannot reproduce dosage programs.
- Worse, any real-data misfit to the parametric nuisance family appears only in real tumours, because decoys are perfectly specified by construction. Examples: cell-level CNA mosaicism, nonlinear QC effects, mis-estimated ambient profiles, cell-size scaling of gene classes. Such misfit gets typed 'CERTIFIED NON-GENETIC'.
- Decoys only have teeth for engines that ignore nuisance (Gavish, cNMF). The authors must choose:
  - either model nuisance in the engine and drop the decoy claim;
  - or generate decoys from nuisance the engine and test do NOT condition on, re-estimating the nuisance on the decoys with the same noisy estimators.
- F2. THE SHARP CONDITIONAL-INDEPENDENCE NULL MAKES MEMBERSHIP A FUNCTION OF CELL NUMBER.
- H0_t(w) says the program genes are exactly independent given N_t. In real tumours this is false for essentially every 50-gene set.
- Z for an equal-weight 50-gene program with average residual correlation rho (my calculation from M5):
  - rho = 0.01: Z ≈ 2.5 at 50 cells, 11 at 1,000, 25 at 5,000;
  - rho = 0.005: Z ≈ 5.5 at 1,000 cells.
- So A(w), R(w) and the gene lists (M8: list length 'follows the data') mostly measure n_c and unmodelled structure, not program activity.
- Random-program genomic control does not fix this. It divides by sqrt(lambda), which assumes null VARIANCE inflation as in GWAS. Here the violation is a NON-CENTRALITY that grows like sqrt(n_c) and depends on each gene set's background correlation, so the correction is miscalibrated in both directions:
  - weak real programs are over-deflated;
  - gene sets in correlated blocks are under-corrected (ribosomal, mito, histone, IEG).
- Pilot P3 will fail as written. lambda_t ≤ 1.3 at 1,000 cells needs the mean residual correlation among random mean-matched genes to be ≤ ~3e-4, smaller than typical unmodelled cell-size or compositional effects. This last point is my estimate and needs measuring.
- F3. PROP 2 (knockoff+ FDR) DOES NOT FOLLOW FROM THE STATED ASSUMPTION.
- Knockoff+ / target-decoy competition needs null signs to be independent fair coins given |W|: the 'equal chance assumption', independent across hypotheses.
- Here the swap is tumour-level. Exchanging t with t~ flips EVERY program's W at once, so null programs driven by the same nuisance (depth/QC/ambient axes) have perfectly coupled signs.
- Swap invariance of 'the joint law of the pipeline output' also requires real tumours to be distributed like their decoys. That fails whenever real biology is present, which is the whole point.
- Decoys are drawn from a model fitted ON THE REAL TUMOUR: ridge-GCV QC slopes, a b_tg baseline, and dispersion conditional on a rank-10 NMF of the same data. They therefore carry the projection of real biology onto nuisance covariates, so the per-pair symmetry is broken in a direction that depends on the program.
- There is one decoy draw per tumour, so W is high-variance.
- The 'empirical FDR' is then validated mainly on Backbone 1, simulated from the same decoy model. Under arm (i) that is calibrated by construction, so it is circular.
- As stated, the method has no valid or even approximately justified FDR claim for MP existence.

**Major concerns**

- THINNING: THE ASYMMETRY BETWEEN ARMS FEEDS DIRECTLY INTO W.
- Decoys are exact NB draws thinned with their true r. Real UMI data are close to Poisson in measurement (Sarkar & Stephens 2021), with overdispersion from biology.
- Beta-binomial thinning with finite r of a Poisson-like truth gives Cov(X^A, X^B) = -eps(1-eps)mu^2/(r+1). With eps = 0.6 the correlation is:
  - -0.20 at mu = 2, r = 3;
  - -0.38 at mu = 5, r = 3.
- For programs learned largely from a few tumours' fold A (exactly the subset-shared programs pitched as the advantage), this gives a positive selection bias in fold-B cross-gene products. That is roughly Z ≈ 2 to 5 at 1,000 to 5,000 cells for |rho_AB| = 0.2 (my calculation).
- It occurs ONLY in the real arm, so W is biased upward. The same applies to private-slot 'novelty' in M11, where w is learned on the same tumour's fold A.
- 'Second order' is true for programs driven by many tumours, not for subset-shared or private ones. Cell splitting would remove the issue entirely.
- CNA CIRCULARITY IS NOT SOLVED BY USING CLONE-LEVEL OFFSETS.
- Clone labels and delta come from expression-smoothing callers (infercnvpy, CopyKAT, SCEVAN). Clone pseudobulks include the program genes, and Gavish reports that 24% of subclones are associated with a single MP.
- So kappa*delta absorbs state programs whose genes lie in altered segments, which in carcinomas means 20-40% of the genome.
- Positional biological clusters (HLA-II 6p21, histones 6p22 tied to S-phase, keratins 12q13/17q21) are known pseudo-CNA signals.
- D^full vs D^tech uses the same contaminated delta, so biology can be typed 'dosage'.
- The Numbat check on 2 datasets is too small. With about 30 typed programs, kappa = 0.6 has SE ≈ 0.15, so the 95% CI spans about 0.3-0.9.
- NUISANCE COVARIATES ARE BIOLOGY.
- The 'complexity residual' tracks developmental potential and stemness (CytoTRACE [U]), and cycling cells detect more genes.
- %ribo tracks translation/MYC; %mito tracks OXPHOS, hypoxia and stress.
- GCV-tuned gene-specific gamma_tg will absorb these programs, in both the engine offsets and the decoys.
- Expect power loss and mis-typing exactly for stemness, differentiation, OXPHOS and translation MPs. Positive controls may only pass because cell cycle is very strong.
- THE DOUBLET INJECTION PENALISES REAL MIMICRY PROGRAMS.
- Mixing real same-sample TME profiles into decoy malignant cells will reproduce recurrence for malignant programs that share genes with TME cell types: pEMT/CAF-like, GBM MES-like (myeloid-like genes), MHC-II/IFN.
- These get labelled 'technical' or 'QC-entangled'.
- The mixing form X~ <- Thin(X~, 1-omega) + Thin(X_TME, omega) keeps depth fixed, so it models contamination rather than doublets, which add counts.
- pi_t from scDblFinder on data that were already doublet-filtered is unidentified.
- THE FAITHFULNESS DIAGNOSTICS AND PILOTS ARE LARGELY TAUTOLOGICAL.
- P1 classifier on depth, %mito, %ribo, detected genes and clone: decoys COPY s_tc, q_tc and k(c), and the NB fit matches the marginals, so AUC ≈ 0.5 by construction.
- The marginal mean/variance/zero-fraction QQ is matched by any NB fit.
- '|Z| < 2 in decoys' holds by construction, and P3 calibration 'on decoys' is guaranteed.
- None of these tests the joint (cross-gene) nuisance structure, which is what drives false MPs. As written, P1 cannot fail.
- THE P2 GO CRITERION IS MET BY WHAT CURATORS ALREADY REMOVED.
- LQ and contamination MPs reproduced by decoys satisfy '≥2 MPs with ≥30% decoy members' without showing any value beyond manual curation.
- Because decoys inherit biology projected onto QC and CNA covariates, the decoy fraction of a retained MP partly measures leakage, not artefact.
- P2 must target MPs that were RETAINED by curators and include a leakage control.
- THE SHARED-VS-PRIVATE SPLIT IS NOT IDENTIFIED BY THE DATA.
- With a free per-tumour baseline b_tg and J = 3 private slots per tumour (about 4,500 private programs), any subset-shared program can be represented either as one W row or as private V_t rows with the same likelihood.
- The split, and therefore R(w) and certification, is set by the ARD hyperparameters (a0, beta0, J) and K_max. Held-out deviance tends to favour larger K, so K_max is likely to sit at its ceiling.
- Prop 3 is an immediate consequence of the definition: cone(H) ⊆ cone([H; H']). It does not address the real issue, which is SUFA-type identifiability of shared vs private components with offsets and noise.
- TYPING ERROR IS UNCONTROLLED.
- 'Fails W^full but passes W^tech' is a difference between a significant and a non-significant result, with no test of R~^full vs R~^tech.
- The two contrasts are run as two separate knockoff selections with no joint error control.
- The modifier thresholds (eta^2 > 0.3, OR > 3, 0.2 ≤ rho < 0.5) are new hand-set constants. This partly re-imports the threshold stack the paper criticises.
- BASELINE IS CONFOUNDED WITH BATCH.
- pi^ITH and between-tumour usage beta_t treat b_t as patient biology. Within a study, each tumour is usually its own library or run, with its own dissociation time, ambient and capture biases.
- Gene-specific protocol capture differences (nuclear transcripts in sn; 3'/5'; Smart-seq2 full-length) land in b_t and are then decomposed onto W.
- H (MM usage) and beta (KL-NNLS on b_t) are on different scales, so the variance partition is not interpretable without multi-library-per-tumour data.
- THE VALIDATION MOSTLY SHOWS SELF-CONSISTENCY, NOT 'REAL BIOLOGY'.
- 3CA to TCCA replication with the same nuisance model and decoys will reproduce the technical typing tautologically.
- Cross-protocol rests on Janesick 2023 (a single FFPE block [V]) and Oliveira 2025 (n = 5 CRC [V]). About 6 tumours gives no power for per-MP detection-rate or CMH claims. The 3CA sc/sn 'pairs' are mostly different tumours confounded with study.
- Tahoe:
  - the 10 pre-registered drug-MP pairs are canonical pathways and close to certain to pass;
  - cytotoxic drugs raise %mito and produce low-quality/apoptotic cells, so 'technical' programs also respond, which confounds the certified-vs-uncertified Fisher test;
  - in vitro responsiveness shows gene-set coherence, not that the program exists as an in vivo ITH state.
- Kinker CCLE has its own ambient, doublet and low-quality artefacts, so 'technical types absent' is not a clean prediction.
- The 11 inds_rm labels:
  - one is a mis-annotation, not an artefact;
  - the 7 LQ MPs are trivially flagged by their correlation with %mito/%ribo;
  - the AUROC 0.85 target has 95% CI about 0.70-0.99 (Hanley-McNeil SE ≈ 0.074 with 11 positives), and there is no trivial-baseline comparator.
- Stability across thinning replicates of the same cells is close to guaranteed and not evidence of robustness.
- The ≥85% re-certification target under 50% study subsampling is unrealistic for subset-shared programs, because power falls; failing it does not mean instability.
- DATA AVAILABILITY IS UNDER-AUDITED.
- The count model and thinning need raw integer UMIs. 3CA distributes harmonised matrices: CPM/log and log2(TPM/10+1) are used in its methods. The availability of raw counts per study, and the share of Smart-seq2/TPM-only studies, is unverified. Many canonical MP-defining studies are Smart-seq2 [U].
- Read counts cannot be data-thinned. 'Certify by projection' still needs the M5 test, which needs thinning.
- There are no empty droplets, so the ambient profile is a pooled-cell proxy dominated by malignant cells. Sorted samples lack TME cells for the doublet injection.
- BAMs for Numbat are mostly controlled-access (dbGaP/EGA).
- TCCA has 36 studies, 537 patients AND 183 cell lines [V snippet]. Its overlap with 3CA v1/v2 is unquantified, and the independent patient cohort may shrink to a few studies.
- THE COMPUTE ESTIMATE IS WRONG, THOUGH NOT FATAL.
- One dense denominator (cells × genes × K = 4.5e6 × 8e3 × 200) is about 1.4e13 flops. H and W updates together need about 3e13 flops per epoch, not 2-3e12.
- This is still minutes per epoch on an H100 with TF32 matmuls, but the 3× augmentation, 3 replicates, the grid and baselines on augmented cohorts should be re-budgeted from a measured prototype.
- The real bottleneck is human time, not GPU time.
- THE TEAM AND TIMELINE ARE NOT REALISTIC.
The plan asks 1.5-2 FTE to deliver, in 15 months:
- the nuisance model, decoys, thinning, tests, knockoffs, a GPU online NMF, gene meta-analysis, a projector with SE, an R bridge and docs;
- ingestion of about 200 heterogeneous datasets plus TCCA, CanSig, Kinker, Tahoe-100M, two spatial sets and TCGA;
- 9-10 containerised, fairly tuned baselines on a 4-arm simulation grid with 50 null cohorts plus real data;
- Numbat from FASTQ;
- pre-registration.
I estimate 3.5-5 FTE-years, or 24-30 months with 2 FTE. The plan also has 648 combinations of robustness knobs, which is a garden-of-forking-paths risk unless the defaults are frozen and pre-registered before the 3CA run, not at M11.
- GENE MEMBERSHIP (M8) INHERITS F2.
- The sharp null plus BH over 8,000 genes, with about 1e5 cells per MP, will flag thousands of genes with negligible effects.
- DerSimonian-Laird over tumours ignores clustering of tumours within studies, which is anti-conservative.
- Choosing member tumours on fold B and then testing genes on the same fold B is a selection effect that is small per gene but systematic.
- THE STUDY RULE IS CONFOUNDED.
- The '≥2 studies' rule is confounded with protocol: protocol-specific artefacts recur across studies that share a chemistry and pass (a).
- Lineage programs from cancer types covered by a single study can never certify.
- The <50-cell exclusion and the 40% B-fold bias membership toward large 10x tumours, and so bias R toward studies with many large tumours.

**Prior-art threats**

- ClusterDE / 'Synthetic control removes spurious discoveries from double dipping…' (Song, Li et al.; bioRxiv 2023 [V]; a 2026 version titled 'Synthetic Control Enables Reliable Cluster Validation and Marker Discovery in Omics Data' [P on venue]). It passes a synthetic null through the identical pipeline and controls FDR with a target-vs-null contrast. CERT-MP's decoy logic is this idea moved to the tumour level. The difference must be argued on the null unit, recurrence and typing, not on the principle.
- Clipper (Ge…Li, Genome Biology 2021) [V this session]. Contrast scores between two conditions with knockoff+-style thresholding and no p-values. M6 (W = R - R~, knockoff+ threshold) is structurally a Clipper instance and must be cited as such.
- Target-decoy competition and its formal link to knockoffs:
- Emery, Hasam, Noble, Keich, 'Multiple competition-based FDR control' (arXiv 1907.01458) [V];
- 'Competition-Based Control of the False Discovery Proportion' (Biometrics 2023) [V];
- He et al., 'theoretical foundation of target-decoy search' [V title].
This literature has already worked out what CERT-MP's Prop 2 skips: the equal-chance assumption must hold independently per hypothesis, and multiple decoys gain power and stability. Not a scoop, but reviewers will say 'decoys' is target-decoy competition, and that CERT-MP violates its key assumption.
- Molecular cross-validation (Batson, Royer, Webber; bioRxiv 2019) [V]. Binomial splitting of UMIs for scRNA model selection and held-out deviance. Prior art for the held-out tuning of zeta and K.
- NB count splitting (arXiv 2307.12985) [V], data thinning (Neufeld et al., JMLR 2024) [V], and Hivert, Agniel, Thiébaut, Hejblum, 'Running in circles: practical limitations…of data fission and data thinning' (arXiv 2405.13591) [V]. The last shows that component-specific overdispersion under mixtures cannot be estimated well enough to make the folds independent. It applies directly to the conditional-r approach.
- Sarkar & Stephens, 'Separating measurement and expression models…' (Nature Genetics 53:770–777, 2021) [V]. Argues for a Poisson measurement model for UMI data. This undercuts NB thinning as the default and supports Poisson thinning or cell splitting.
- CS-CORE (Su et al., Nature Communications 14:4846, 2023) [V]. Co-expression testing aware of depth and measurement error. It documents inflated false positives in naive residual-correlation approaches. M5 is a co-expression test and must be compared against it.
- Derandomised knockoffs:
- Ren & Barber, e-value derandomisation (arXiv 2205.15461) [V];
- Ren, Wei, Candès, 'Derandomizing knockoffs' [U].
Averaging W over 3 thinning replicates is neither scheme, so the claim 'in the spirit of' needs its own justification.
- GBCD (Nature Genetics 2025) [V], ITHmapper (bioRxiv 2026) [V], TCAT/starCAT (Nature Methods 2025) [V]. These are positioned correctly in the concept. The risk is the scoop the authors already identify: any of these groups adding a ClusterDE-style null to its engine reproduces most of CERT-MP's headline, because the decoy layer is engine-agnostic and easy to copy.
- CytoTRACE (Gulati et al., Science 2020) [U]. Gene counts and complexity track developmental potential. This is evidence that the 'complexity residual' nuisance covariate absorbs real biology.
- 'Resolving sensitivity, specificity and signal contamination in Xenium spatial transcriptomics' (Nature Methods 2026) [V title only]. In situ data carry their own contamination and segmentation-spillover artefacts. Xenium-based 'dissociation-free' validation needs explicit controls.
- Data sizes. Janesick et al. 2023 Nature Communications is a SINGLE FFPE breast block [V]. Oliveira et al. 2025 Nature Genetics 57:1512–1523 has n = 5 CRC samples [V]. TCCA has 36 studies, 537 patients plus 183 cell lines [V snippet]. These figures limit the claimed cross-protocol and replication evidence; they do not pre-empt the method.

**Required changes**

- DECOUPLE DECOYS FROM ENGINE AND TEST CONDITIONING. Pick one:
- (a) The engine uses only depth offsets, and nuisance is handled solely by decoys plus typing. This is the audit design, extended to cNMF, GBCD and others.
- (b) Decoys are generated from nuisance that the engine and test do NOT condition on, AND nuisance is re-estimated on decoys with the same noisy estimators: CNA re-called on decoy counts against real TME references; ambient, QC slopes and dispersion re-fitted.
In either case, show on real 3CA data that R~ is clearly > 0 for known technical programs under the flagship engine. Otherwise the knockoff step is vacuous.
- REPLACE THE SHARP NULL WITH AN EFFECT-SIZE (INTERVAL) NULL.
- Example: H0 is that the mean pairwise residual correlation among program genes, or the fraction of residual variance explained by the program score, is ≤ delta0.
- Calibrate delta0 against random gene sets matched for gene-gene background correlation, or against the same sample's TME cells.
- Report membership as effect sizes with CIs.
- Show that membership rate is flat against malignant-cell count.
- Drop sqrt(lambda) genomic control, or justify it against a non-centrality model.
- USE CELL SPLITTING INSTEAD OF COUNT THINNING FOR M5 AND M8.
- Learning w on half the cells and testing coherence on the other half is exact under within-tumour cell independence, needs no dispersion model and works for Smart-seq2 read counts.
- Keep thinning only where the same cells must be scored. There, use Poisson thinning for UMIs and apply the identical estimation procedure in both arms.
- Quantify fold dependence with Cov = -eps(1-eps)mu^2/(r+1) for subset-shared and private programs.
- MAKE THE FDR CLAIM DEFENSIBLE.
- Generate m ≥ 5 independent decoy replicates per tumour and use multiple-competition or mirror statistics.
- Show empirically that null W signs are close to independent across programs that share nuisance (sign-correlation diagnostics).
- Or move to e-value / conformal formulations whose assumptions can be checked.
- Validate FDR on Backbone 2 and on 'plant into REAL tumours' designs, not on the decoy generator itself.
- State plainly which assumptions are untestable.
- ADD A BIOLOGY-LEAKAGE CONTROL.
- Estimate QC slopes gamma and ambient from the same sample's non-malignant cells, which share technical effects but not malignant programs, or from genes outside candidate programs.
- Quantify decoy reproduction of programs PLANTED into real tumours: the target is R~/R ≈ 0.
- Show that complexity, %ribo and %mito covariates do not absorb stemness, differentiation, OXPHOS or translation MPs, or else drop them.
- BREAK CNA CIRCULARITY.
- Estimate delta while leaving out candidate program genes, or use allele-informed calls.
- Expand the allele-aware validation to at least 6-8 datasets and 60 or more tumours, with a pre-specified power calculation for kappa.
- Run a direct test of R~^full vs R~^tech for dosage typing instead of contrasting a significant and a non-significant result.
- REPLACE THE TAUTOLOGICAL FAITHFULNESS CHECKS WITH JOINT-STRUCTURE TESTS.
The current checks (marginal QQ, a classifier on copied features) cannot fail. Candidates:
- residual co-variation among technical-only gene classes, real vs decoy;
- decoys of TME populations vs real TME technical axes;
- barnyard or cell-line-mixture data with known doublets and ambient;
- technical-replicate or low-heterogeneity cell-line data as near-null tumours.
Make P1 and P3 falsifiable.
- REDEFINE PILOT P2'S GO CRITERION. Require decoy reproduction of MPs that curators RETAINED in published catalogues, net of planted-program leakage. Re-discovering inds_rm adds nothing beyond manual curation.
- OBTAIN REAL DISSOCIATION-FREE EVIDENCE AT ADEQUATE n.
- Matched fresh sc / sn / Flex on at least 10-15 tumours. This needs a wet-lab or core collaborator, which contradicts 'no wet lab' and should be budgeted.
- And/or within-malignant-cell in situ co-expression on public Xenium 5K / CosMx whole-transcriptome tumour datasets, with segmentation-spillover controls.
- Pre-specify per-MP power.
- AUDIT DATA AT M0 BEFORE COMMITTING:
- raw-UMI availability per 3CA study;
- Smart-seq2 / TPM-only fraction;
- samples lacking TME cells;
- TCCA overlap with 3CA and its cell-line content;
- controlled-access status of BAMs.
Re-plan the replication cohort if independent TCCA patients are few.
- STRENGTHEN THE ARTEFACT CALIBRATION.
- Report the AUROC with CIs and against a trivial baseline (MP-score correlation with %mito/%ribo/doublet/TME-marker scores). CERT-MP must beat that baseline.
- Expand the labelled set: Gavish 2023 removals, Kinker's complexity programs, GeneNMF dropMetaPrograms examples, and protocol-named MPs in newer catalogues.
- Remove the mis-annotation case from the positives.
- FIX THE TAHOE ANALYSIS.
- Use QC-matched pseudobulks and exclude, or model, cytotoxicity and apoptosis.
- Pre-register non-canonical pairs where the outcome is uncertain.
- Make the primary endpoint 'certified beat technical in specific response', not '8/10 canonical pairs'.
- ADDRESS SHARED-VS-PRIVATE IDENTIFIABILITY DIRECTLY.
- Show that certification and R(w) are stable across a0, beta0, J and K_max.
- Either prove a real identifiability condition for the shared block with per-tumour baselines and private slots, or demote Prop 3 to a remark.
- FIX GENE MEMBERSHIP.
- Use a minimum-effect null.
- Use three-level meta-analysis (cells within tumours within studies) with REML/HKSJ rather than DerSimonian-Laird.
- Report the core gene list's length against cohort size.
- FREEZE AND PRE-REGISTER DEFAULTS before any 3CA-wide run, not at M11: eps or the split fraction, K_max, J, CNA caller, alpha, q and the modifier thresholds. Report the knob grid only as a sensitivity analysis.
- RE-BUDGET COMPUTE AND PEOPLE.
- Measure a prototype epoch (expect about 3e13 flops per epoch before the 3× augmentation).
- Add about 1 FTE, or cut the CERT-NMF engine and projector to a v2. A certification layer on existing engines is the lower-risk headline.
- Plan 24-30 months to submission.

**Ideas worth keeping**

- The tumour as the null unit. Nuisance-matched decoy tumours, run LABEL-BLIND through the UNCHANGED standard pipeline (the P2 audit), are the strongest, cheapest and most publishable component. Doublet and ambient injection built from the same sample's real TME cells is especially concrete.
- certify_external as an engine-agnostic audit layer over Gavish, GeneNMF, cNMF, GBCD and ITHmapper outputs. It survives model scoops and gives broad uptake, provided it is fixed per F1-F3.
- A held-out coherence statistic that uses only cross-gene products, applicable to any fixed gene set. It is truly valid for query projection when W_ref is learned on independent data: per-sample detection without thinning. This is a useful, defensible deliverable once it moves to an interval null.
- Study-capped recurrence R = sum_s min(1, |A ∩ s|/2), a simple, sensible guard against one cohort dominating.
- Replacing top-50 truncation with effect-size plus heterogeneity gene tiers (core vs context, I^2) exported as confidence-tiered GMTs. This would change how MSigDB GAVISH_3CA sets are used, given a minimum-effect null.
- Decoy FAMILIES (full vs clone-flat vs user-added artefact families) for mechanistic typing. The idea is good once decoys carry nuisance the engine does not condition on, and once typing is a direct contrast.
- A benchmark with planted confounders, explicit MISSPECIFICATION arms, nuisance-only null cohorts and an independent backbone. Reported honestly, this is a real contribution beyond CanSig.
- Panel-retention iota_k for Flex/Xenium transfer, a small but practical addition.
- The sequencing discipline: build the decoy generator first, keep the M8 audit paper as a fallback, explicit go/no-go decisions, tripwires that accelerate the preprint, and pre-registration of the orthogonal validations.
- Separating within-tumour from between-tumour (baseline) variation instead of centring it away. Keep it, but estimate it only where multi-library or multi-region tumours can separate patient from batch.
