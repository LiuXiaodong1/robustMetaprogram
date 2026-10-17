# Landscape and Gap Analysis: building a robust, model-based tumour metaprogram (MP) package aimed at Nature Methods

> **Status:** AI-assisted literature synthesis (2026-10-02; revised 2026-10-17, see the revision note below), produced from five parallel literature sweeps.
> Verification tags: [V] found online, [V+] re-checked from source code/README, [P] partly verified, [U] unverified.
> bioRxiv/PubMed/arXiv were partly unreachable during the search, so treat absence-of-prior-art claims as weak and re-check before citing.


As of 2026-10-02. This synthesises 5 independent literature sweeps, plus a re-check pass in which I could reach only GitHub.

> **Revision note (2026-10-17).** Revised in response to the external review of 2026-10-16 (mainly §3 of that review, "Proposal documents", plus the confounders raised in its §2.2–2.4), using a follow-up citation check (WebSearch snippets plus GitHub; publisher, PMC, bioRxiv and arXiv pages were still unreachable, so abstract wording is near-verbatim at best). What changed:
> - **TCAT/starCAT** no longer described as "T cells only": its authors present starCAT as generalising to other cell types and tissues, and further references sit on Zenodo (§0, §2B, §5, §6.1, G7). The Flex caveat stays.
> - **Oliveira 2025** corrected: the same FFPE blocks (CRC and normal tissue) profiled by Chromium Flex, Visium HD and Xenium; no fresh-dissociated malignant cells, so it is a fixed-dissociated vs in situ comparison, **not** a fresh-vs-fixed twin. Whether Janesick 2023 has a fresh-dissociated arm is flagged as unverified (§4.1, §4.3, G3).
> - **Dissociation references** upgraded from [U] to [V] with one-line content; **Baron 2020** added with the correct reading (it argues the stress-like cancer-cell state is real); Neuschulz 2023, FixNCut 2024, Trinks 2024 and the ccRCC fresh/frozen/snFlex preprint added with hedged details (§4.1, §6.1, G3, Appendix B).
> - **TCCA** described from its code: Gavish-style clustering of single-seed per-sample NMF, then a hand-coded whitelist of 43 MPs; a second example of irreproducible manual curation; cell-line samples must be excluded for tumour replication (§0, §1.2, §2A, §4, §6.1, G3).
> - **Scoop and competitor update** (§2, §3, §6): ISOMERIC, cEBMF, DECIPHER, SIMVI, arXiv 2410.06451, Gene Program Dictionary, scPFA, PORCELAN, Spectra's artefact motivation, and the protocol-comparison papers that pre-empt parts of the stress result. A search-engine-only scoop check in October 2026 found no calibrated or provenance-typed MP method and no audit of 3CA (weak evidence).
> - **Confounders** (§0, §1.2, §4, G2, G3): malignant→TME ambient spillover; expression-derived covariates as a caveat for decoys; gene-length / nuclear-fraction and abundance / %ribo axes, with the user-reported HER2 snRNA-seq evidence labelled as such. Numbers from the `mpaudit` toy simulator are labelled as toy results.
> - The P1 GO / PIVOT / INCOMPLETE rule is not restated here; see the pointer under "Pilot analyses" (single source of truth: `src/mpaudit/report.py`, rule `P1-GO v1.0`).

---

## 0. Executive summary

- **The standard is a stack of thresholds, not a model.** The Gavish 2023 / 3CA pipeline is still unchanged in 3CA v2 (Nature Cancer 2025). I read its code and confirmed this. It runs:
  - per-tumour NMF on log2, gene-centred, negative-zeroed expression;
  - ranks 4–9, 10 runs each, keeping the best fit;
  - top-50 gene truncation;
  - fixed filters: 35/50 overlap across ranks, 10/50 across tumours, at most 10 shared genes for redundancy;
  - greedy clustering into MPs, with a 10-gene join rule and at least 5 programs per MP;
  - manual curation at the end. The 3CA v2 code removes 11 MPs by hand: `inds_rm <- c(5,7,11,17,25,30,36,44,58,68,74)`, annotated as 7 low-quality, macrophage, T-cell, endothelial and one mis-annotation. The TCCA atlas (bioRxiv 2026) ends the same way: a hard-coded whitelist of 43 MPs, with its criteria only in a code comment [V+ code].
- **"Robust" means "recurrent", and that is the core flaw.** Anything that recurs systematically passes by design. That includes dissociation stress, ambient RNA (including malignant transcripts leaking into TME cells), doublets, quality and complexity axes, gene-length / nuclear-fraction and abundance / %ribo axes, and positional programs driven by recurrent arm-level copy-number alterations (CNAs). The cross-tumour filter (at least 10 of 50 shared genes from a 7,000-gene universe) guards against a random-gene null that no real dataset follows (my calculation, Section 1.4).
- **The closest principled competitor is GBCD** (Liu…Stephens, Nature Genetics 2025). It is a label-free model that separates shared from patient-specific programs. It works on Gaussian log data, has a user-set Kmax, is fitted per cancer type, and has not been shown at atlas scale. It offers:
  - no calibrated test of whether a program recurs across tumours;
  - no modelling of CNAs or artefacts;
  - no reference or projection tool.
- **The closest Nature Methods template is TCAT/starCAT** (Kotliar…Raychaudhuri, 2025): a reproducible catalogue built across datasets, plus a projection tool. The catalogue is of T cells, but the authors present starCAT as generalising the framework to other cell types and tissues (near-verbatim abstract wording, not yet read on the article), and further references are hosted on Zenodo (which ones was not checked). The authors caution that their references are not established for 10x Flex.
- **The closest new large-scale competitor is ITHmapper** (Abelson lab, bioRxiv Jan 2026). It uses scVI with sample as batch, Hotspot modules, discrete states and a mapping tool across more than 1,000 tumours (counts from one sweep; the review of 2026-10-16 could not verify them, and a 2026-10 search did not find the preprint).
- **Statistical machinery that exists but has never been applied to tumour MPs:**
  - multi-study factor models with sharing by subsets of studies (Tetris, APAFA, SUFA);
  - count-native multi-study NMF with learned rank and covariates (Grabski et al. Genome Biology 2025, for mutational signatures);
  - power posteriors for NMF that is robust to model misspecification (BayesPowerNMF);
  - data thinning and count splitting for held-out rank choice and gene inference, and FDR control by data splitting after clustering (arXiv 2410.06451);
  - nuisance-preserving synthetic nulls (ClusterDE) and knockoff calibration (recall);
  - NMF identifiability theory.
- **The defensible novelty is not "another factorisation".** It is making the MP a statistical object: a definition, an estimator, uncertainty and error control, with built-in deconfounding of genetic (CNA/clone) and technical (artefact or protocol) signal. It should come with a benchmark that plants confounders and with an uncertainty-aware reference and projector.
- **Caveat on the search.** The sweeps' web search budget ran out, and nature.com, bioRxiv, PubMed and arXiv were proxy-blocked. Absence of a scoop is therefore weak evidence. A search-engine-only re-check in October 2026 (about 14 targeted searches) also found no calibrated, error-controlled or provenance-typed MP method and no audit of the 3CA/Gavish MPs, but it could not run bioRxiv full-text, PubMed or citing-paper searches. A manual bioRxiv / Google Scholar check is still required before committing (Section 6.3).

**Verification legend** (used throughout):
- **[V]** At least one sweep (or, for items added on 2026-10-17, the follow-up citation check) found the title, venue and key facts online (live search result, publisher snippet or GitHub). Publisher full text was generally not reachable.
- **[V+]** Also re-checked by me in this pass from live GitHub content.
- **[P]** Partly verified: the work exists, but the venue, authors or a specific detail comes from memory or from a single sweep's snippet.
- **[U]** Not verified online. Included only as context and clearly flagged.

---

## 1. Current standard pipelines and their precise failure modes

### 1.1 The canonical pipeline (Gavish et al., Nature 618:598–606, 2023 [V]; 3CA v2, Tyler et al., Nature Cancer 6:1088–1101, 2025 [V])

| Step | What the code does (tiroshlab/3ca; 3CA_v2/MP_analysis_new_batch.R [V+]) |
|---|---|
| S0. Malignant calling | Classify malignant cells from CNA inference before NMF (3CA practice). The per-tumour minimum cell count was not verified. |
| S1. Gene universe and transform | Top 7,000 genes by mean log2 expression. Each gene is centred by its row mean within the sample, and negative values are set to 0. |
| S2. Per-tumour NMF | `nmf(rank=4:9, method='snmf/r', nrun=10)` keeps the best of 10 restarts. This gives 39 programs per tumour. |
| S3. Truncation | Each program is reduced to its top 50 genes, and the loadings are discarded. |
| S4. "Robust program" filter | `intra_min=35`: shares at least 35/50 genes with a program at another rank in the same tumour. `intra_max=10`: a greedy redundancy pruning. `inter_min=10`: shares at least 10/50 genes with a program in another tumour. |
| S5. Clustering into MPs | Greedy. The founder is the program with the most partners at 10 or more shared genes; programs join while they share at least 10 genes with the evolving MP gene list; clusters need at least 5 programs. MP genes are the top 50 by frequency, with ties broken by NMF score. |
| S6. Curation | Gavish 2023: MPs that look low-quality (ribosomal or mitochondrial), single-study or doublet-like are removed after the fact. v2: single-study MPs and MPs with fewer than 5 programs are dropped, then **11 MPs are removed by index** (LQ×7, Macro, T-cells, Endo, mis-annotation) [V+]. New MPs were matched to old ones by Jaccard above 0.3 with all other overlaps below 0.2 [V]. |
| S7. Downstream use | MP gene lists are used for module scoring and as MSigDB-style gene sets. No uncertainty is attached. |

### 1.2 Map from step to failure

| Step and parameter | Failure mechanism | Observable consequence | Evidence |
|---|---|---|---|
| S0 malignant calling | Doublets and contaminating TME cells are misclassified as malignant | Contamination MPs (macrophage, T-cell, endothelial) | 3CA v2 removed exactly these by hand [V+]. Gavish removed doublet-like MPs [V]. |
| S1 top-7,000 by **mean** expression | Universe is biased toward highly expressed housekeeping, ribosomal, mitochondrial and immediate-early genes | "Low-quality" ribosomal/mito MPs. Programs of lowly expressed transcription factors are missed. | 7 of the 11 hand-removed v2 MPs were "LQ" [V+] |
| S1–S4, no axis for **gene length / nuclear fraction** or **transcript abundance** | Mechanism (hypothesis): in single-nucleus data, long genes (more intronic, nuclear-retained reads) co-vary with each nucleus's intronic fraction, and abundant short cytoplasmic transcripts co-vary with %ribo. Both axes are systematic, so they recur across tumours and pass the recurrence filters. | Technical programs that reference matching names as biology | **User-reported, unpublished, not reproduced here:** in the user's HER2 breast-cancer snRNA-seq project, two pre-registered de novo MP searches gave mostly technical axes: several "luminal" programs of long genes (median gene length 317–437 kb), and an abundance axis of short cytoplasmic genes tracking %ribo in 29/30 patients (rho 0.68–0.77) that reference matching labelled "EMT". Published quantifications of these axes were not searched for in this revision. |
| S1 **centring per tumour plus zeroing negatives** | Only within-tumour relative up-regulation survives. A program active in every cell of a tumour has zero centred variance. Not a likelihood; depth and complexity are not modelled. | Programs that differ between tumours (subtype or baseline states) are invisible, and programs cannot be compared across tumours. A complexity axis becomes a "program". | GBCD found subtype programs that tumour-by-tumour NMF missed [V]. Kinker removed complexity-correlated programs by hand [V]. Sweep 1 reports that the Tirosh lab's 2025 GBM paper models "baseline programs" separately [P]. |
| S2 **fixed rank grid 4–9** | Tumours with many states are under-factorised (programs merged). Simple tumours are over-factorised (split or noise programs). | Hybrid MPs, and dependence on rank | BioBombe: different k capture complementary biology, so no single k is optimal [V] |
| S2 **best of 10 restarts** | NMF is not identifiable when cells are not "sufficiently scattered" (overlapping, correlated programs, few pure cells). Many near-equivalent optima exist, and best-of-n picks one arbitrarily. Restart variance is never reported. | Genes near rank 50 change between runs and between software versions | Kinker README warns NMF output changes with R and package versions [V]. Identifiability theory of Fu, Huang and Sidiropoulos [V/P]. |
| S2 per-tumour fits | No borrowing of strength. Small tumours give noisy programs. | Rare MPs are missed. Noise passes when partners exist by chance. | Structural |
| S3 **top-50 truncation** | Loadings are thrown away. List length is fixed whatever the program's sparsity. No significance per gene. | Gene membership is all-or-nothing and unstable. Large programs are truncated and small ones padded with noise. | Structural. GeneNMF v0.6 moved to cosine similarity on weights for this reason [V]. |
| S4 `intra_min=35` (across ranks) | Favours dominant axes that survive any rank change: cell cycle, stress, quality, clone identity. Penalises fine programs that appear only at higher K. | A systematic selection bias toward strong axes, including technical ones and CNA-driven ones | Structural argument. No sweep found a published quantification, so this is a pilot target. |
| S4 `inter_min=10` (across tumours) | Calibrated against nothing. Under random gene sampling, P(overlap ≥10) ≈ 1e-12 (Section 1.4), so the filter effectively passes **any shared structure**, artefactual or genetic | Recurrent artefacts and positional programs from recurrent arm-level CNAs (e.g., chr7 gain in GBM) can pass | The positional-program mechanism is a hypothesis; no sweep found a positional-enrichment test applied to MPs |
| S4 `intra_max=10` | Greedy and order-dependent pruning | Which redundant program survives is arbitrary | Code [V+] |
| S5 **greedy clustering** with an evolving gene list | Depends on order and founder. Chaining. No merge/split test. The number of MPs is a by-product of thresholds. | MPs change with cohort composition; 3CA v1→v2 needed a matching rule | v2 matching criterion (Jaccard >0.3, others <0.2) [V] |
| S6 **manual curation** | Expert judgement replaces inference and cannot be reproduced | Every major pipeline ends in hand curation | 3CA v2 `inds_rm` [V+]. TCCA (bioRxiv 2026) keeps a hard-coded whitelist of 43 MPs, with its removal criteria stated only in a code comment [V+ code]. Kinker hand-made metaclusters [V]. Barkley chose final modules by hand [V]. GeneNMF ships `dropMetaPrograms()` [V]. |
| S6 CNA checked only afterwards | Genetic programs are neither removed nor labelled | Clone-driven states are mixed with plastic states | Gavish: only 24% of genomic subclones were associated with a single MP [V]. Kinker: 39% of discrete clusters but only 8% of continuous NMF programs were associated with CNA subclones [V]. |
| S7 scoring | Module scores carry no uncertainty and are not comparable across tumours (because of centring) | Clinical and spatial associations are noisy | Structural |

### 1.3 Variant pipelines and where each one fails

- **Kinker 2020** (Nature Genetics 52:1208–1218 [V]). Cell lines; k = 6–9; the same 35/10/10 filter. Programs are hierarchically clustered (distance = 50 − overlap) alongside a DBSCAN analysis. Complexity-correlated and cell-cycle programs are removed by hand, and metaclusters are defined by hand (genes in at least 25% of member programs). The result is 12 recurrent heterogeneous programs (RHPs).
- **Barkley 2022** (Nature Genetics 54:1192–1201 [V]; sweep 4 marked it unverified, but sweeps 1 and 5 verified it). Correction to the task framing: this is **not cNMF**. It uses nsNMF with `seed='ica'` and `nrun=1`, so it is deterministic but untested for robustness. Ranks run 2–25 with a "rank minus number of modules" rule. A module is kept if its Jaccard similarity is above 0.05 with at least 3 others, applied iteratively. An infomap gene co-occurrence graph is clustered, and the 16 final modules (including "Stress" and "Metal") are chosen by hand.
- **GeneNMF** (Carmona lab; CRAN; no methods paper, cited through Yerly et al. bioRxiv 2024 [V]):
  - RcppML, one fixed seed (123), 2,000 HVGs with an optional blocklist;
  - cosine similarity on weights, Ward hclust, then `cutree` at a **user-set nMP**, which forces every program into some MP;
  - genes kept with `min.confidence` 0.5 and `weight.explained` 0.5;
  - the authors claim per-sample NMF removes the need for batch correction.
- **Lightweight re-implementations**:
  - scalop (Tirosh lab, R) and pyscalop (Python port, May 2026) [V];
  - mpnmf (Python, Apr 2026, deterministic NNDSVDa initialisation) [V];
  - sciNMF (Brief Bioinform 2025) [V];
  - SSpMosaic (Cell Genomics 2025; differential co-expression programs linked by network propagation and used as integration anchors) [V].

  All are threshold-driven, with no uncertainty and no confounder model. The "faster or deterministic Gavish" niche is saturated.

### 1.4 Systemic diagnosis: the conceptual failure modes a new method must name

1. **Recurrence is necessary but not sufficient.** Under a random-gene null with a 7,000-gene universe, expected overlap between two top-50 lists is 0.36 genes, and P(overlap ≥10) ≈ 1.1e-12. With 2,000 HVGs the figures are 1.25 genes and 1.8e-7 (my hypergeometric calculation, not from a paper). The filter therefore rejects only random noise, which is not the relevant null. Structured nulls (shared quality axes, highly expressed gene blocks, recurrent CNAs, stress) pass. **The null that matters is a "nuisance-only" tumour, and no pipeline simulates one.**
2. **No likelihood means no uncertainty.** Thresholds stand in for inference. No MP has a p-value or FDR for existence, gene membership, tumour membership or number of MPs.
3. **Per-tumour centring and factorisation blocks joint inference.** It cannot borrow strength, cannot see programs active across a whole tumour, and can only reject patient-private programs (genotype or clone) indirectly, through thresholds.
4. **The selection filters are biased, not just noisy.** Cross-rank stability enriches for dominant axes, technical ones included.
5. **Results depend on the cohort.** The MP set is a function of which tumours are included and of greedy ordering.
6. **Validation is post hoc.** CNA association, spatial analysis, cell lines and clinical associations are all checked after MPs are defined, never used to define them.

---

## 2. Competitor table

Notes: "Count model" means a Poisson, negative-binomial (NB) or multinomial likelihood on raw counts. "Confounders" means explicit modelling, not post-hoc filtering. Unverified details are flagged inline.

### 2A. Direct multi-tumour malignant-program / MP methods

| Method (status) | What it solves | What it does NOT solve | Multi-sample? | Count model? | Confounder handling | Validation |
|---|---|---|---|---|---|---|
| **Gavish 2023 / 3CA** (Nature 2023) [V+ code] | Pan-cancer recurrent ITH programs at scale (77 studies, 1,163 tumours, 41 MPs [V]); simple and interpretable | Seed, rank and threshold dependence; uncertainty; programs that differ between tumours; artefacts; CNAs | Per-tumour, then greedy clustering | No (log2, centred, negatives zeroed, Frobenius snmf/r) | Manual and post hoc (LQ, doublets, single-study); CNA only post hoc | Recurrence; comparison with normal-epithelium MPs and Kinker RHPs; CNA-subclone association post hoc |
| **3CA v2** (Tyler et al., Nature Cancer 6:1088–1101, 2025) [V; code V+] | Expanded cohort (124 datasets, 2,836 samples; 149 MPs across 8 cell types [P: count from preprint and press]); some sc vs sn comparisons | Same method; 11 MPs removed by hand | Same | No | Single-study filter plus manual index removal | Matching to v1 MPs; sc vs sn for some MPs |
| **Kinker 2020** (Nature Genetics) [V] | Cell-intrinsic programs without TME or dissociation | Hand-made metaclusters; thresholds | Per cell line, then hierarchical clustering | No | Complexity and cell-cycle programs removed by hand | Senescence-like subpopulations, drug sensitivity, comparison with tumours |
| **Barkley 2022** (Nature Genetics) [V] | Graph/community aggregation instead of greedy; links states to the TME | Robustness (single run); manual selection; CNA | Per tumour, then infomap on gene co-occurrence | No | None | Visium spatial; mouse models (IFN module depends on location and drops with lymphocyte depletion) |
| **GeneNMF** (CRAN; no methods paper) [V] | Fast, Seurat-native, weight-based (cosine) MP clustering with diagnostics (silhouette, sample coverage) | Data-driven nMP; seed variance; CNA, doublets, stress (blocklist only); forced assignment | Per sample, then Ward hclust | No | Optional gene blocklist; manual drop | Demos; CanSig found its signatures **less rediscoverable** than top early-integration methods [V] |
| **scalop / pyscalop / mpnmf / sciNMF / SSpMosaic** [V] | Re-implementations or variants; determinism (mpnmf NNDSVDa) | All of the inferential gaps | Per sample, then clustering | No | None | Minimal; none systematically benchmarked |
| **CanSig method** (Boeva lab; bioRxiv 2022) [V] | Integrates patients, clusters, derives DE signatures, places them next to CNV profiles | Continuous or overlapping programs; CNV handled only for malignant calling and side-by-side display | Early integration | No | CNV used for malignant calling only | See the CanSig benchmark (Section 4) |
| **ITHmapper / Pan-Cancer ITH compendium** (Nofech-Mozes…Abelson, bioRxiv 2026) [V; my re-check: the manuscript repo has a placeholder for the preprint link, and the mapping tool supports 12 cancer types; cell, tumour and cancer-type counts P, not verifiable per the 2026-10-16 review] | 3.6M cells, more than 1,000 tumours, 14 cancer types; scVI (sample as batch, 10 seeds) plus Hotspot modules; multi-resolution clustering; query mapping | Discrete states rather than overlapping programs; sample-as-batch can erase programs active across a whole tumour; no uncertainty; CNA used only for malignant calling (CopyKAT) | Integration-based | scVI uses NB; modules do not | Sample batch key only | Survival, Visium HD, metastasis enrichment, Kinker lines, Tahoe-100M |
| **GBCD** (Liu, Carbonetto, Willwerscheid, Oakes, Macleod, Stephens; Nature Genetics 57:263–273, 2025) [V+ README citation] | Label-free separation of patient-specific, dataset-specific and shared programs; gene signatures with posterior uncertainty (lfsr); beats tumour-by-tumour NMF on shared subtype programs | Kmax user-set (authors suggest 10–40 [P]); cell-by-cell covariance limits scale; one cancer type per analysis; no tumour-membership credible sets or recurrence test; no CNA or artefact modelling; no projection tool | Pooled, joint | No (Gaussian on log data) | Implicit (nuisance falls into "patient-specific" components) | HNSCC plus 3 PDAC studies; simulations; stress program prognostic of PDAC survival independent of stage and subtype |
| **TCCA MPs** (González-Bermejo…Al-Shahrour, CNIO; bioRxiv 2026, DOI 10.64898/2026.02.18.706316) [V; code V+ (cnio-bu/tcca)] | 43 MPs over 36 public studies, with a portal and drug-vulnerability annotation | Same threshold stack, a single NMF seed, and a hand-coded MP whitelist (see §4.1) | Per sample, then clustering with parameters identical to the Gavish defaults | No (RcppML NMF on log2(CP100K+1), row-centred, negatives zeroed; ranks 4–9; seed 123, no restarts) | SCEVAN for malignant calling only; artefact criteria (mito/ribo-rich, single-study, <10 programs) stated in a code comment and applied by hand | Not checked (full preprint unreachable) |
| **ISOMERIC** (bioRxiv 2025, DOI 10.1101/2025.04.21.649842; authors and publication status not checked) [V] | Splits tumour-cell expression variation into copy-number-dependent and CN-independent latent parts (expression autoencoder with a CN decoder and an adversarially trained network); 189,027 tumour cells from 10 public cohorts (pancreatic, breast, ovarian, lung); argues that CN dependency confounds the search for shared states | Not MP discovery across tumours; no recurrence test or error control; no technical-artefact typing | Pooled across cohorts | Not checked | CN, inferred (e.g. inferCNV) or measured (e.g. G&T-seq) | Some CN-dependent dimensions track mutational signatures and transcriptional subtypes (ovarian, pancreatic); both parts linked to TME composition |
| **KINOMO** (bioRxiv 2022; AACR 2024 abstract) [V] | KL-divergence NMF with kernel correction for genomic variability between patients | Not a full MP framework; low visibility; unpublished | Pooled | KL (Poisson-like) | Variability between patients (kernel) | Not checked in detail |
| **UVG** (Yoon & Nam, Brief Bioinform 2024) [V] | Gene selection that reduces malignant clusters specific to one sample | Not program discovery | n/a | n/a | Heuristic gene selection | n/a |

### 2B. General program-discovery engines (baselines and components)

| Method (status) | What it solves | What it does NOT solve | Multi-sample? | Count model? | Confounders | Validation |
|---|---|---|---|---|---|---|
| **cNMF** (Kotliar, eLife 2019) [V] | Consensus over seeds within one dataset; identity vs activity programs | Aggregation across tumours; K picked by eye; confounders | No | No (Frobenius default) | None | Simulations with planted programs |
| **TCAT / starCAT** (Kotliar…Raychaudhuri, Nature Methods 22:1964–1980, 2025; DOI 10.1038/s41592-025-02793-1) [V] | Reproducible consensus programs across 7 datasets (1.7M T cells, ~700 individuals, 38 tissues; 46 consensus programs, 52 in the reference); fixed reference plus NNLS projection. The authors present starCAT as generalising this to other cell types and tissues [near-verbatim abstract wording; check on the article before quoting], and the README lists several references besides the T-cell one, hosted on Zenodo (which ones was not checked) | The catalogue itself is T-cell only, and we know of no malignant-cell reference (the Zenodo references were not checked); singleton programs rejected by heuristic plus judgement (49 of 52 excluded); not a count model; no per-cell uncertainty; authors warn about transfer to 10x Flex | cNMF per dataset on Harmony-corrected non-negative data, then clustering | No | Batch correction (adapted Harmony) | Cross-dataset reproducibility (9 programs in all 7 datasets, mean R=0.81); CITE-seq and activation [P, from memory in sweep 5] |
| **scHPF** (MSB 2019) [V] and **consensus scHPF** (Levitin…Sims, bioRxiv 2023) [V] | Count-native Bayesian Poisson factorisation; consensus over seeds and K gives automatic K | Patient, batch and CNA structure; still a preprint | Pooled | **Yes** (Poisson) | None | 52 glioma slice cultures from 19 patients |
| **scDEF** (Ferreira, Kuipers, Beerenwinkel; bioRxiv 2022, v2) [V; authors from sweep 2; peer-review status **unknown**, the README shows no journal citation] | Multi-level Gamma-Poisson (hierarchical Bayesian NMF); automatic size; batch-aware; optional gene-set guidance | Atlas scale unproven; shared vs patient-specific only implicit; no CNA | Batch-aware | **Yes** | Batch | 9 HGSOC patients (patient-specific and shared states) |
| **LIGER iNMF / UINMF / consensus iNMF** (Cell 2019; Nature Communications 2022; rliger 2024) [V] | Explicit split into shared W and dataset-specific V_i; scalable (online) | Count likelihood; K and λ chosen by hand; V_i not sparse; weak identifiability with about 1,000 "datasets"; no subset sharing | Yes (joint) | No | Dataset-specific term | Integration benchmarks |
| **scINSIGHT** (Genome Biology 2022) [V] | Gene modules common to all conditions or specific to one | Needs condition labels; no subset sharing | Yes | No | Condition term | Simulations and condition data |
| **MOFA+** (Genome Biology 2020) [V by sweep 2; sweep 3 did not verify] / **MuVI** (AISTATS 2023) [V] / **MOFA-FLEX** (bioRxiv 2025) [V] | ARD and structured sparsity; group-wise factor activity; priors from gene sets (MuVI, MOFA-FLEX) | Loadings are not non-negative; MOFA+ groups share loadings, so patient-private programs cannot be represented; no recurrence test | Groups | MOFA+ offers Poisson; MuVI is Gaussian | Groups or views | Multi-omics demos |
| **Spectra** (Kunes…Pe'er, Nature Biotechnology; online Sept 2023, print 2024 42:1084–1095) [V] | Programs guided by priors; cell-type-specific factors; η score for dependence on the prior. Motivated partly by artefacts: the authors describe unsupervised factorisation as "prone to technical artifacts" [snippet wording] | Needs priors and cell-type labels; patient and CNA structure; likelihood not verified | No | Not verified | Can absorb known nuisance programs through priors | Beat Slalom, scHPF, NMF and expiMap; immunotherapy and TIL datasets |
| **cEBMF** (Denault, Tayeb, Carbonetto, Willwerscheid, Stephens; arXiv 2505.11639, NeurIPS 2025) [V] | Covariate-moderated empirical Bayes matrix factorisation: side information (covariates, images, text, graphs) shapes adaptive priors on the factors through any probabilistic model or neural network; NMF and semi-NMF are special cases; extends the EBMF/flashier family that GBCD builds on; PyTorch code | Not a tumour-MP method: no cross-tumour recurrence test, no provenance typing, no application to tumour MPs found | Via side information | Not checked | Known covariates moderate the priors (could take nuisance covariates directly) | Simulations, spatial transcriptomics, collaborative filtering |
| **Gene Program Dictionary** (Xu…Ma, bioRxiv 2025, DOI 10.1101/2025.10.29.685322) [P: search snippets only] | 1,534 consensus programs from 466 scRNA-seq datasets, with permutation-based FDR | Not tumour-MP provenance typing; nuisance handling not checked | Many datasets | Not checked | Not checked | Not checked |
| **expiMap** (Nature Cell Biology 2023) [V] | Program-tied VAE plus reference mapping | Few de novo programs; database-bound | CVAE covariates | scArches family (default not verified) | Covariates | Atlas query mapping |
| **DRVI** (Moinfar & Theis, bioRxiv 2024 and later; in scvi-tools) [V; no journal version found] | Disentangled, interpretable factors **and** integration in one model; inactive dimensions give a soft automatic K | Nonlinear programs (no single sparse loading); no shared-vs-patient structure; no CNA | Batch covariates | **Yes** (NB) | Batch | HLCA decomposition; integration benchmarks |
| **LDVAE** (Bioinformatics 2020) [V] | Scalable NB linear factor model | Signed, dense loadings; rotation | Batch covariate | **Yes** (NB) | Batch | — |
| **Hotspot** (Cell Systems 2021) [V] | Gene modules by local autocorrelation with no K | Hard, disjoint modules; depends on the graph (pooled graphs are dominated by patient) | Inherits from its latent space | Null-model Z-scores | None | vs WGCNA and GRNBoost |
| **mosaicMPI** (NAR 2024) [V] | Multi-resolution cNMF across cohorts and platforms; network "themes"; annotation propagation | Heuristic thresholds; no uncertainty; no CNA or artefacts; not benchmarked against Gavish | Yes (network) | No | None | GBM cross-platform |
| **MetaTiME** (Nature Communications 2023) [V] | ICA meta-components across 79 TME datasets, plus annotation tool | TME only; signed components; fixed correlation thresholds | Per dataset, then clustering | No | None | Annotation of new data |
| **CSI-GEP** (Cell Genomics 2025) [V] | GPU, atlas-scale NMF; K and programs chosen by replication across data splits | No aggregation across tumours; splits by cell; stable technical programs still pass | No | No | None | Recovery of simulated ground-truth programs |
| **sciRED** (Nature Communications 2025) [V] | Poisson-GLM regression of known covariates, then varimax factors; matches factors to covariates; flags "unexplained" factors | Signed factors; regressing out sample also removes real biology; no uncertainty | Covariates | Poisson GLM for nuisance | **Known covariates** | Interpretability metrics |
| **EcoTyper** (Cell 2021) [V] | Automatic rank (cophenetic coefficient) and automatic false-positive state filter (AFI) | Bulk-deconvolution setting; discrete states | Yes | No | AFI filter | Clinical, bulk, sc and spatial |
| **DeltaNMF** (Karpurapu, Gersbach, Singh; bioRxiv 2026; RECOMB 2026 listing [P]) [V] | Baseline (control) programs vs case-specific programs; foundation-model gene-similarity Laplacian prior; about 20× faster than cNMF on GPU | No multi-tumour hierarchy, CNA handling or uncertainty | Case/control | No | None | CAD Perturb-seq; 2.4× coherence on PPI |
| **contrastiveVI** (Nature Methods 2023) [V] / **bcNMF** (arXiv 2026) [V] | Variation enriched in target vs background (NB VAE; contrastive NMF) | Needs a matched background; not applied to tumour vs normal cell of origin | Target/background | contrastiveVI NB; bcNMF Frobenius | Background subtraction | Case-control datasets |
| **fastTopics plus GoM DE** (Stat Comput 2026; Genome Biology 2023) [V]; **Topyfic** (bioRxiv 2024) [V] | Count-native topics; per-topic log fold-changes with uncertainty (principled "which genes define a program"); reproducible topics over many seeds | K chosen by user; no patient terms | No | **Yes** | None | — |
| **GSFA** (Nature Methods 2023) [V] | Spike-and-slab sparse factors with calibrated PIPs | Gaussian on transformed data; MCMC scale; perturbation designs only | Covariates | No | Covariates | CRISPR screens |
| **Tripso** (Lotfollahi lab; bioRxiv 2026) [V; authors not retrieved] | Cell embeddings per program from a self-supervised transformer | Compute; not malignant ITH; initialisation from predefined programs [P] | — | No | — | TRM program checked in spatial data |
| **scGPT gene programs; FM sparse autoencoders** (Nature Methods 2024; bioRxiv 2025 / ICLR 2026; arXiv 2603.02952) [V] | Gene-embedding priors; SAE features include technical axes | Programs averaged over contexts; zero-shot embeddings underperform simple baselines (Kedzierska, Genome Biology 2025 [V]) | Implicit | No | SAE "technical" features (hint only) | Enrichment |

### 2C. Peripheral (sample-level or multicellular; useful downstream)

- DIALOGUE (Nature Biotechnology 2022) [V]
- scITD (Nature Biotechnology 2024/2025) [V]
- MOFAcell (eLife 2023) [V]
- Tensor-cell2cell (Nature Communications 2022) [V]
- LEMUR (Nature Genetics 2025) [V]: a covariate-aware embedding template

None of these resolves programs inside malignant cells. They are useful for validating that MPs couple to TME programs.

**Not to be cited as a competitor without checking.**
- **Glioma myeloid programs (Nature 2025).** Sweep 1 mentions it in passing as using cNMF per cohort and then hierarchical clustering into 14 consensus programs. No citation was retrieved [P].
- **Tirosh lab GBM paper, Nature Genetics 57:1155–1167 (2025), "baseline programs".** Only sweep 1 reports it [P].

### 2D. Adjacent disentangling and lineage methods (added 2026-10-17)

None of these discovers or types recurrent tumour MPs, but each answers a neighbouring question that reviewers will raise.
- **SIMVI** (Dong, Kluger H, Fan, Kluger Y; Nature Communications 16:2990, 2025; DOI 10.1038/s41467-025-58089-7) [V]. Annotation-free deep generative model (built on scvi-tools) that separates cell-intrinsic from spatially induced latent variables. Data include MERFISH cortex, Slide-seq v2, Slide-tags tonsil, spatial multiome melanoma and a CosMx melanoma cohort.
- **DECIPHER** (Xia, Cao, Gao; Nature Communications 2025; DOI 10.1038/s41467-025-63140-8) [V]. Separate embeddings for a cell's molecular identity and its spatial context (neighbouring cells), trained with cross-scale contrastive learning; scales to atlases of millions of cells.
  - Both are spatial counterparts of separating malignant-intrinsic from TME-shared variation. That split is asked here of recurrent programs in dissociated scRNA-seq; SIMVI and DECIPHER are candidate tools for the Xenium/CosMx rung.
- **PORCELAN** (Schlüter & Uhler, Nature Communications 2025) [P: from search snippets] finds lineage-informative genes and subtrees. It is a possible competitor or tool for clone/heritability typing (G2).
- **scPFA**, "Lineage-informed factor analysis reveals heritable programs of single-cell gene expression" (bioRxiv, Sept 2026; DOI 10.64898/2026.09.21.753290) [P: title only; **abstract not retrieved**]. Going by the title, it separates heritable from transient programs, which would overlap clone/dosage typing (G2). Read it before claiming novelty on that axis.

---

## 3. Statistical tools available but not yet applied to tumour MP discovery

None of the sweeps found any application of the tools below to tumour scRNA-seq MPs.

| Tool (status) | What it provides | MP failure it addresses | Adaptation needed / caveat |
|---|---|---|---|
| **BMSFA** (De Vito et al., AoAS 2021) [V] | Shared loadings Φ plus study-specific Λ_s; number of factors learned by multiplicative-gamma shrinkage | Patient-private programs absorbed by Λ_s instead of becoming MPs | Gaussian, not non-negative; a factor is shared by all or private to one; Gibbs sampling scales poorly |
| **Tetris** (Grabski et al., AoAS 2023) [V] | An **IBP prior over which studies own which factors**; factors shared by any subset; credible balls for the sharing pattern | Directly matches "an MP recurs in a subset of tumours"; replaces greedy clustering with a posterior over tumour membership | Gaussian; very slow (Stat Med 2026 tutorial: "takes far more time to run") |
| **Bayesian multi-study NMF** (Grabski, Trippa, Parmigiani, Genome Biology 26:98, 2025) [V]; probit follow-up (arXiv 2502.01468) [P] | **Count-native (Poisson) NMF** across studies with a binary factor-inclusion matrix; learned number of signatures; "recovery-discovery" (known plus new); covariates on exposures | Blueprint for multi-tumour NMF: Gavish MPs as soft priors, new MPs discovered, nuisance as covariates | Built for 96 channels; needs SVI or minibatching, NB likelihood, and nesting of cells within tumours within cohorts |
| **APAFA** (Bortolato & Canale, JASA 2026) [V] | Shrinkage priors that allow factors private to sub-study partitions, private to one study, or shared, with identifiability conditions | Natural nesting: subclone-private, tumour-private, cohort-private (likely technical), shared across cohorts (MP candidate) | Gaussian, Gibbs; not single-cell |
| **SUFA** (Chandra, Dunson, Xu, JASA 2025) [V] | Identifiability of the shared vs specific split (specific loadings in the shared span) | Without an identifiability condition, "shared" MPs depend on the prior and initialisation | The subspace assumption is doubtful for amplicon-private programs; the method needs its own identifiability condition |
| **PFA** (Roy et al., AoAS 2021) [V] | Study-specific perturbation of shared factors, plus a divergence statistic | Platform and cohort distortions modelled as perturbations, not new MPs; flags outlier datasets | Dense p×p perturbations; Gaussian |
| **Multi-study FA benchmark** (Stat Med 2026; BLAST arXiv 2502.14600; CAVI arXiv 2305.13188) [V; arXiv author lists unconfirmed] | Simulation designs and scalable comparators | Answers the reviewer question "why not just use MSFA or Tetris?" | Bulk scale and Gaussian |
| **Multi-study sparse VAE** (Moran & Krishnan, ICASSP 2026) [V] | Amortised, sparse decoder; shared vs group-specific factors | Scalable engine for many tumours | Not non-negative; bulk demo; no uncertainty on sharing |
| **BayesPowerNMF** (Xue, Miller, Carter, Huggins, PLoS Comput Biol 2026) [V] | **Power (coarsened) posterior** that suppresses spurious components caused by misspecification; automatic K | scRNA-seq misfit (overdispersion, ambient RNA, doublets) creates "junk" programs that **recur because the misfit is systematic** | ζ needs calibration (e.g., held-out likelihood on count-split folds); MCMC cost |
| **bayesNMF** (Landy, Basava, Parmigiani, JCGS 2026) [V]; **compressive NMF** (Zito & Miller, arXiv 2404.10974) [V] | Learned rank with posterior uncertainty; informed priors | Replaces "run NMF at ranks 4–9 and intersect" | Single study; low-dimensional demos |
| **Data thinning** (Neufeld et al., JMLR 2024) [V]; **generalised data thinning** (arXiv 2303.12931 [V]; JASA [P]) | Independent folds X = X1 + X2 for convolution-closed families; held-out evaluation for low-rank models | Rank chosen by held-out NB deviance; validating programs on independent noise | Needs known or estimated nuisance parameters (overdispersion) |
| **Count splitting** (Neufeld et al., Biostatistics 2024; NB version arXiv 2307.12985) [V] | Estimate latent structure on one fold, test genes on the other, avoiding double dipping | **FDR-controlled MP gene membership** instead of "top 50"; held-out K and number of MPs | Theory covers clusters, pseudotime and PCs, not NMF or aggregation across tumours |
| **Limits of thinning** (Hivert et al., arXiv 2405.13591) [V] | Shows fold dependence under mixtures when overdispersion is mis-estimated | A reviewer will raise this; calibration must be demonstrated | Extension to factor models is open (a possible theoretical contribution) |
| **ClusterDE synthetic null** (Song, Li et al., bioRxiv 2023; software paper indexed in PubMed) [V] | Simulated null through the identical pipeline (scDesign3), then target-vs-null contrasts with FDR | **Empirical FDR for each MP** using "nuisance-only" null cohorts (depth, QC, CNA blocks, cell cycle kept; other co-expression removed) | Validity depends on how faithful the null is (the main risk for any FDR claim) |
| **recall knockoffs** (DenAdel…Crawford, AJHG 2025) [V] | Artificial null genes give a knockoff+ threshold | Knockoff genes added before NMF give a null for gene loadings; detects K too large | Knockoffs drawn from marginals lack correlation, so validity is heuristic |
| **Selective inference** (Gao, Bien, Witten JASA [V]; Chen & Witten JMLR 2023 [P authors]; pipelines arXiv 2603.18413 [V]) | Exact p-values after clustering | Testing MP merges and splits | NMF with restarts and top-k truncation gives intractable selection events; best treated as a theoretical companion |
| **FDR via data splitting for testing after clustering** (Wang, Lin, Zhao; arXiv 2410.06451, 2024) [V] | Split the data in halves, run the usual test-after-clustering on each half, and combine the results into a mirror statistic per feature that is symmetric about zero under the null; multiple splits improve stability and power; scRNA-seq applications included | FDR-controlled MP gene membership, and a split-and-combine template for tumour-level cross-fitting | Asymptotic FDR control proven only in Gaussian settings; no journal version found; not yet extended to factorisation or aggregation across tumours |
| **NMF identifiability** (Fu, Huang, Sidiropoulos; arXiv 1709.00614, 1803.01257) [V; venues P] | Sufficiently scattered and anchor conditions; volume-regularised NMF | Explains seed instability mechanistically. Diagnostic per program; a possible new theory that **pooling heterogeneous tumours restores identifiability of shared programs** | Noiseless theory; the condition is hard to check on noisy data |
| **staNMF** (Wu…Yu, PNAS 2016) [V]; stability selection (Meinshausen & Bühlmann JRSS-B 2010) [P] | Choosing K by stability; subsampling-based error control | Stability across combined perturbations (seeds, cell bootstrap, thinning, cohort hold-out) | Stability is not validity: technical programs are very stable, so stability must be paired with null calibration |
| **GoM DE** (Genome Biology 2023) [V] | Per-program log fold-changes with uncertainty under grade of membership | Principled gene membership | Single dataset |
| **Random-signature nulls** (Venet et al., PLoS CB 2011) [V] | Most random signatures associate with outcome | Clinical validation of MPs must beat matched random signatures | Bulk-era |
| **HDP / IBP nonparametric factor models** (Knowles & Ghahramani, arXiv 1011.6293 [V; AoAS venue P]; HDP across scRNA datasets arXiv 2212.02505 [V; authors P]) | Global program pool reused by groups in different proportions | "HDP-Poisson factorisation over tumours" is the cleanest generative definition of an MP | Over-splits under misspecification (pair with power posteriors); scalable implementation needed |

**Highest-leverage combination** (each piece verified to exist; the combination was not found anywhere):
- count-native multi-study NMF with subset-sharing inclusion variables (Grabski + Tetris);
- tumour-private and cohort-private slots (APAFA);
- power posterior (BayesPowerNMF);
- ζ and K chosen by count-split held-out likelihood;
- per-MP empirical FDR from nuisance-preserving synthetic nulls (ClusterDE-style);
- an identifiability diagnostic for each program.

---

## 4. Ground-truth and benchmark resources, and metrics

### 4.1 Resources

| Resource (status) | What it provides | Use in evaluation | Caveat |
|---|---|---|---|
| **CanSig benchmark** (Boeva lab, Cancer Research 86(4):873, 2026, DOI 10.1158/0008-5472.CAN-25-0940) [V by sweep 1 only; I could not find the repo under BoevaLab in a re-check] | 13 methods; 12 datasets; 5 cancer types; 185 patients; 174k malignant cells. Metrics: bio-conservation, batch removal, signature correlation, rediscovery, cross-dataset reproducibility, clinical relevance. Harmony, BBKNN and fastMNN scored best; scalop and GeneNMF were less rediscoverable; late integration improved with more samples. | Established harness: reuse its datasets and metrics and beat both early- and late-integration methods | Ground truth limited to known signatures; seed and rank stability not studied; Gavish code not benchmarked directly |
| **3CA v1/v2 corpus plus 11 hand-removed MPs** [V+] | Largest curated pan-cancer scRNA-seq set; MP lists | Flagship re-derivation; **the 11 removed MPs are a labelled artefact set** for automatic QC; cohort-composition stability tests | Labels reflect one lab's judgement |
| **TCCA** (González-Bermejo…Al-Shahrour, CNIO; bioRxiv 2026, DOI 10.64898/2026.02.18.706316; code cnio-bu/tcca) [V; code V+; 3CA overlap U] | 36 public scRNA-seq studies, 34 tumour types, about 1.84M cells (1,089,024 malignant, 750,218 non-malignant) in 853 samples: 670 from about 540 patients (537 in another source) plus **183 cell-line samples**. SCEVAN malignant calls; scANVI integration; drug annotation. **43 MPs**: per-sample RcppML NMF (log2(CP100K+1), row-centred, negatives zeroed, ranks 4–9, a single seed, no restarts; genes kept if log2 CP100K > 3.5 in more than 2% of cells, not the top 7,000), clustered with parameters identical to the Gavish defaults, then a hard-coded whitelist of 31 "cancer_states" and 12 "lineage_specific" MP IDs (out of at least 57), renumbered MP1–MP43. The removal criteria (mito/ribo-rich, single-study, fewer than 10 programs) appear only in a code comment and are not coded as a test. | Independent replication cohort: learn on 3CA, replicate on TCCA's **patient tumours only**, excluding studies shared with 3CA and all cell-line samples. Its 43 MPs are a standard-pipeline baseline, and its whitelist is a **second example of irreproducible manual curation** (after 3CA v2 `inds_rm`) | Overlap with 3CA is not stated in the README or in the abstract text retrieved; count it from GEO accessions at P0. A single NMF seed makes its programs seed-dependent |
| **Kinker CCLE scRNA** (198 lines) [V] | Cell-intrinsic programs with no TME or dissociation | Test of recurrence in cell lines | Cultured cells lack in vivo programs |
| **Tahoe-100M** (bioRxiv 2025) [V] | About 100M cells, cancer lines × drugs × doses | Causal test: do MPs respond to drugs with matched mechanisms (CDK4/6 → cycle, MEK → MAPK, etc.) against DMSO nulls? | In vitro; dominated by drug response |
| **Replogle 2022 Perturb-seq** (K562, RPE1) [P: venue from memory; datasets confirmed] and **scPerturb** (Nature Methods 2024) [V] | Causal interventions with non-targeting nulls | Perturbation-specificity metric | Atypical lines; mostly essential genes |
| **PerturbNMF** (Engreitz lab, Zenodo 2025–26) [V] | Nine program-quality metrics with calibrated perturbation nulls | Template for the biological-validity panel | Not peer reviewed; not tumour-specific |
| **Matched multi-platform tumours (fixed-dissociated vs in situ).** Janesick et al. 2023 (FFPE breast: Chromium Flex, Visium, Xenium; Nature Communications [P venue]) [V; not re-checked on 2026-10-17]. Oliveira et al., Nature Genetics 57:1512–1523, 2025 (DOI 10.1038/s41588-025-02193-3; GEO GSE280318) [V]: the **same FFPE blocks** profiled by Chromium Flex (260,506 cells dissociated from serial FFPE sections of 8 CRC and normal-tissue blocks), Visium HD (5 patients, plus adjacent normal tissue from 3) and Xenium (4 tumours). | Fixed-dissociated (Flex) scRNA plus in situ data on the same blocks | **Fixed-dissociated vs in situ invariance**: real MPs replicate; ambient and doublet MPs should not replicate in situ (which has its own segmentation spillover) | **Neither is a verified fresh-vs-fixed twin.** Oliveira has no fresh-dissociated malignant cells (all tumour tissue is FFPE; the material used for its immune-profiling arm was not confirmed), so it cannot test loss of ex vivo dissociation programs. Whether Janesick has a fresh-dissociated arm from the same tumour was not verified; do not use it as a fresh-vs-fixed twin until checked. Two cancer types, few patients, probe panels, and no same-protocol replicate pairs for a protocol null |
| **Candidate fresh-vs-fixed tumour twins** (added 2026-10-17). Trinks et al., Cell Oncology 47:1221–1231, 2024 [P]: patient-matched fresh, cryopreserved and FFPE lung cancer; FFPE blocks showed higher cell diversity than fresh tissue, and a senior-author interview says stress signatures were higher in fresh epithelial cells. Xirenayi et al., bioRxiv 2025, DOI 10.1101/2025.11.25.690485 [P]: ccRCC fresh scRNA-seq vs flash-frozen snMultiome (RNA) vs FFPE snFlex, with "broadly consistent" cell-type profiles. Possible extras, from search snippets only: a 10x GEM-X Flex performance preprint (bioRxiv, DOI 10.64898/2026.03.06.710062; FFPE Flex vs frozen-biopsy snRNA-seq in a clinical-trial set) and a Frontiers 2026 brief report on matched FFPE and frozen samples (reportedly 12 pairs; tumour types unknown) [P]. | Fresh vs fixed (or frozen) profiles from the same tumours | Loss of ex vivo programs (stress) under fixation, if the specimens are patient-matched | Not checked: Trinks's patient count and Flex use; whether the ccRCC modalities come from the same patients, how many, and any stress or malignant-program result. Confirm at P0 before counting any of these as twins |
| **Greenwald et al. GBM spatial** (title and authors V; Cell 2024 [P]) [V] | Visium plus CODEX; MP spatial coherence | Spatial coherence metric; expected niches (hypoxia/MES near necrosis) | Spots mix cell types |
| **Barkley Visium** (GSE203612) [V] | Spatial co-localisation of states | Spatial validation | — |
| **KP-Tracer lineage tumours** (Yang, Jones et al.; bioRxiv 2021 [V]; Cell 2022 [P]) plus **Cassiopeia** (Genome Biology 2020) [V] | Lineage measured independently of expression | Heritability of an MP (phylogenetic autocorrelation): heritable vs plastic | Mouse only |
| **clonealign** (Genome Biology 2019) [V]; **CalicoST** (Ma…Raphael; bioRxiv 2024 per README [V+]; a journal version, believed to be Nature Methods by sweep 4, is **not confirmed**) | Clone labels from DNA or alleles, not circular with expression; spatial clone domains | Clone-independence test; checks whether an "MP domain" is just a clone domain | Needs scDNA or adequate SNP depth |
| **CNA callers: Numbat** (Nature Biotechnology 2022), **CopyKAT** (Nature Biotechnology 2021), **SCEVAN** (Nature Communications 2023), infercna [V] | Clone labels and segment dosage per cell | Nuisance covariates; planted CNA confounders in simulation | Expression-based calls can be circular; Numbat's allele awareness helps |
| **Ambient RNA: CellBender** (Nature Methods 2023), **SoupX** (GigaScience 2020), **DecontX** (Genome Biology 2020) [V] | Ambient profile per sample from empty droplets | Fixed nuisance factor; injection at known fractions, in **both** directions (TME→malignant and malignant→TME) | Correction is incomplete. Ambient RNA comes from the whole suspension, so in malignant-rich samples it carries malignant transcripts into TME cells: a malignant-only program whose strength varies between samples then looks "TME-shared" (review of 2026-10-16, §2.2; a design argument, not a published measurement). TME co-variation tests need a per-sample estimate of malignant-derived ambient and nulls matched on malignant specificity |
| **Doublets: scDblFinder** (F1000Research 2021) plus the Xi & Li benchmark (Cell Systems) [V] | Doublet scores; artificial doublets | Inject malignant–TME doublets at known rates | Homotypic doublets are hard to detect |
| **Dissociation-artefact literature** [V, re-checked 2026-10-16/17 from search snippets; details marked "not checked" are open]: **van den Brink 2017** (Nature Methods 14:935–936; DOI 10.1038/nmeth.4437): a muscle stem-cell subpopulation strongly affected by the dissociation protocol, so some published studies may need reinterpreting. **O'Flanagan 2019** (Genome Biology 20:210): cold active protease (6 °C) vs collagenase (37 °C) on patient tumours, breast-cancer PDXs and cell lines (155,165 cells); a 512-gene collagenase-associated heat-shock/stress core including FOS and JUN (patient tumour types not checked). **Denisenko 2020** (Genome Biology): dissociation and storage biases in sc and sn workflows, in mouse kidney (title-level check). **Slyper 2020** (Nature Medicine): sc and sn toolbox for fresh and frozen human tumours (title-level check). **Marsh 2022** (Nature Neuroscience): enzymatic dissociation induces an ex vivo signature, mainly in microglia. **Neuschulz 2023** (Mol Syst Biol 19:e11147): 4sU labelling during dissociation shows a shared core stress signature plus cell-type-specific responses in mouse and zebrafish tissues; no human tumour data found. **FixNCut** (Jiménez-Gracia et al., Genome Biology 25:81, 2024): reversible DSP fixation before dissociation; fewer stress artefacts reported, with cell composition preserved (tumour content and size of the effect not checked). **Counterpoint, Baron et al. 2020** (Cell Systems 11:536–546.e7) [V]: they acknowledge the dissociation concern but argue, from zebrafish melanoma, spatial transcriptomics and human samples, that a stress-like cancer-cell state (fos, jun, heat-shock genes) is real in vivo and present from early tumorigenesis; their reanalysis of TNBC, oligodendroglioma and PDAC data suggests it is conserved. It does **not** establish that dissociation stress exists. | Immediate-early and heat-shock signatures; protocol contrasts; direct labelling of ex vivo transcription | Stress signature for injection; protocol contrasts as natural ground truth; 4sU design for a labelled subset | The stress MP is plausibly a mixture of ex vivo and in vivo stress, and its per-tumour ex vivo share has not been measured. O'Flanagan already shows in patient tumours that the response exists and depends on protocol, so what remains novel is per-tumour attribution, not existence |
| **Simulators**: Splatter (Genome Biology 2017), splatPop (Genome Biology 2021), SymSim (Nature Communications 2019), scDesign3 (Nature Biotechnology 2024) [V]; SERGIO [P venue]; scMultiSim (bioRxiv 2023 per README [V+]; journal **unconfirmed**) | Planted programs and usage | Backbone of a semi-synthetic simulator | **None jointly models segment CNAs plus subclones, ambient RNA in both directions (TME→malignant and malignant→TME), malignant–TME doublets, cell-varying stress, gene-length / abundance axes and platform effects.** (The `mpaudit.sim` toy used for the P1 pilot now plants simplified versions of these, e.g. a single CNA arm and no TME–TME doublets. It exists to test the audit code and is not a benchmark simulator) |
| **Metric templates**: BEELINE (Nature Methods 2020), scIB (Nature Methods 2022), BioBombe (Genome Biology 2020), momix (Nature Communications 2021) [V]; Saelens 2018 module-detection benchmark [P] | Recognised metric vocabularies | AUPRC plus stability (BEELINE); conservation vs removal (scIB) | Not tumour-MP-specific |
| **Benchmark guidelines**: Weber 2019, Sonrel 2023 (Genome Biology) [V] | Neutral, extensible benchmark design | Reviewer checklist | — |

### 4.2 Metrics: a proposed scorecard

| Axis | Metric | Ground truth |
|---|---|---|
| **Recovery** | Hungarian-matched cosine between estimated and planted loadings; program-gene AUPRC; Spearman of per-cell usage; F1 of recovered programs at a cosine threshold | Semi-synthetic |
| **False-MP rate** | Fraction of reported MPs aligned (cosine above τ) to injected nuisance (CNA block, ambient profile in either direction, doublet mixture, stress, depth or complexity, gene-length / nuclear-fraction axis, abundance / %ribo axis) | Semi-synthetic |
| **Empirical FDR / null behaviour** | Number of MPs called on **null cohorts** that carry only nuisance and no planted programs; calibration of the claimed FDR | Synthetic null |
| **Patient-private leakage** | Fraction of planted tumour-private (e.g., CNA-driven) programs reported as MPs | Semi-synthetic |
| **Calibration** | Coverage of credible intervals or sets for gene membership, tumour membership and number of MPs; nominal vs empirical FDR of gene lists | Semi-synthetic |
| **Stability** | Matched cosine or Jaccard and the variance of the number of MPs across seeds, rank grid, cell bootstrap, count-thinning folds, tumour subsampling and leave-one-study-out; ARI of program-to-MP assignment | Real (3CA, CanSig) |
| **Cohort-composition stability** | MP survival when 10–50% of studies are added or removed (analogue of the v1→v2 matching) | Real |
| **Held-out replication** | MPs learned on half the studies are recovered in the other half; 3CA → TCCA (patient tumours only, minus studies shared with 3CA) | Real |
| **Artefact classification** | AUROC for flagging the 11 hand-removed 3CA v2 MPs (and Gavish's removed MPs) | Real, labelled |
| **Cross-protocol invariance** | Concordance of gene weights (restricted to panel genes) and usage between fresh sc, Flex, sn and Xenium on the same tumour; stress and ambient MPs should fall away | Matched platforms |
| **Genetic orthogonality** | Variance in MP usage explained by DNA- or allele-defined clone; chromosome-arm positional enrichment of MP genes | clonealign/Numbat/CalicoST |
| **Spatial coherence** | Moran's I or Hotspot autocorrelation of usage; pre-registered niche associations (hypoxia near necrosis, EMT at the edge, IFN near T cells), controlling for clone domain | Xenium, Visium HD |
| **Perturbation responsiveness** | Fraction of MPs specifically responsive to at least one perturbation against non-targeting or DMSO nulls; coherence of responding regulators | Tahoe-100M, Perturb-seq |
| **Cell-intrinsic recurrence** | MP usage variance and recurrence in CCLE lines | Kinker |
| **Clinical** | TCGA survival or response association, beating matched random signatures | TCGA |
| **Compute** | Runtime and memory at 10^6–10^7 cells and about 3,000 tumours | — |

### 4.3 Recommended evaluation ladder (in order of evidential weight)

- **A. Semi-synthetic spike-ins on real malignant cells.** Plant programs that are shared, shared by a subset, tumour-private or nested. Add measured confounders: Numbat segments, empty-droplet ambient profiles (including malignant→TME leakage), doublets with real TME cells, a stress signature at varying fractions, depth shifts, gene-length / intronic-fraction and abundance / %ribo axes, and a platform layer. Plant at least one real program that correlates with an estimated nuisance covariate (e.g. cycling with doublet score) to measure false nuisance calls. Also run null cohorts. This is where Gavish-style pipelines should visibly fail on false MPs and leakage.
- **B. Causal checks.** Tahoe-100M and Perturb-seq, with calibrated nulls.
- **C. Cross-protocol invariance** on matched platforms. This is the most direct answer to "is it real biology?". Public data so far support fixed-dissociated vs in situ comparisons (Janesick, Oliveira); a fresh-vs-fixed contrast needs the candidate twins of §4.1 to be confirmed as patient-matched, or new matched data.
- **D. Spatial coherence**, with a clone-domain control.
- **E. Genetic orthogonality and heritability**: DNA-defined clones and KP-Tracer lineage.
- **F. Generality**: CCLE, held-out cohorts (TCCA patient tumours, without its cell lines or studies shared with 3CA), and clinical associations against random-signature nulls.

**Baseline panel to beat:**
- Gavish/3CA code itself (not just scalop)
- GeneNMF
- cNMF, per tumour and pooled
- TCAT-style multi-dataset cNMF with Harmony
- GBCD
- LIGER iNMF / consensus iNMF and scINSIGHT
- MOFA+/MOFA-FLEX
- DRVI
- Spectra or expiMap
- Hotspot and ITHmapper (scVI plus Hotspot)
- CanSig's top early-integration methods (Harmony, BBKNN, fastMNN, each followed by clustering and DE)
- mosaicMPI
- consensus scHPF and scDEF
- TCCA's 43 MPs (standard pipeline plus hand whitelist), and ISOMERIC on the CN-dependence axis where its code runs
- Tetris or BMSFA on a subsample, to show their scaling limits

---

## 5. Publication bar: what comparable papers delivered

| Paper (venue, status) | Conceptual claim | Benchmarks | Orthogonal validation | Deliverables |
|---|---|---|---|---|
| **TCAT/starCAT** (Nature Methods 2025) [V] | Programs reproducible across datasets as a fixed catalogue, plus projection; built on T cells, with starCAT presented as generalising to other cell types and tissues [near-verbatim abstract wording] | Reproducibility across 7 datasets; annotation accuracy | CITE-seq, activation, immunotherapy [P] | pip/R packages, CLI, Zenodo versioned reference, custom-reference builder |
| **Spectra** (Nature Biotechnology 2023/24) [V] | Program discovery guided by priors and specific to cell types | vs Slalom, scHPF, NMF, expiMap; several real datasets (Bassez breast immunotherapy etc.) | Immunotherapy response and TIL factors | pip/bioconda, AnnData-native, full figure-reproduction repo |
| **GSFA** (Nature Methods 2023) [V] | Calibrated sparse factors (PIPs) guided by perturbations | Simulations plus CRISPR screens | Perturbation effects | R package |
| **GBCD** (Nature Genetics 2025) [V+] | Shared vs patient-specific tumour programs without integration | Simulations; vs NMF; HNSCC, PDAC | Prognostic stress program (survival independent of stage and subtype) | R package, workflowr site |
| **contrastiveVI** (Nature Methods 2023) [V] | Salient variation specific to the target | Case-control datasets | — | scvi-tools |
| **CellBender** (Nature Methods 2023) [V] | Generative removal of ambient RNA | Simulated and real data | — | Widely adopted tool |
| **scIB / BEELINE** (Nature Methods 2022 / 2020) [V] | Neutral benchmarks with metric frameworks | Many methods × datasets | — | Extensible benchmark code |
| **scITD, DIALOGUE** (Nature Biotechnology) [V] | Multicellular programs across donors | Cohort stratification | Phenotype replication | R packages |
| **Calibration points below Nature Methods**: MetaTiME → Nature Communications; Hotspot → Cell Systems; CSI-GEP → Cell Genomics; sciRED → Nature Communications; mosaicMPI → NAR; CanSig → Cancer Research; GeneNMF → no methods paper [V] | Resource plus tool, or a single engineering advance, lands below Nature Methods | — | — | — |

**The bar, distilled:**
1. **A conceptual advance a generic tool cannot match**, stated as a definition (here: "an MP is a program shared by a credible subset of tumours, distinguishable from nuisance-only nulls at a stated FDR").
2. **Simulations with ground truth that include the failure modes** (confounders and nulls), not only planted signal.
3. **5–10 strong competitors**, tuned fairly. GBCD, TCAT-style cNMF, GeneNMF and the Gavish code are mandatory.
4. **Several orthogonal validations**: perturbation, cross-protocol, spatial, genetic orthogonality, clinical against random nulls.
5. **Application at atlas scale with a new biological result**:
   - re-derive the 3CA MPs;
   - show which are artefactual, clone-linked or unstable;
   - show which new or split MPs were missed;
   - ideally one wet-lab or perturbation follow-up.
6. **Software and resource**:
   - a pip/conda and/or R package compatible with scverse or Seurat;
   - a reproducibility repo;
   - a versioned MP reference on Zenodo;
   - a projection function with uncertainty and novelty detection;
   - GPU scaling.

---

## 6. Scoop risk: the closest recent work (2024–2026)

### 6.1 Risk table

| Work | Status / date | What it pre-empts | Risk | How to differentiate |
|---|---|---|---|---|
| **GBCD** (Nature Genetics 2025) [V+] | Published | "Shared vs patient-specific tumour programs without per-tumour thresholds"; gene-level uncertainty | **High on that claim alone** | Count likelihood; subset-sharing with tumour-membership credible sets; calibrated recurrence/FDR; CNA and artefact nuisance; pan-cancer scale; projection |
| **ITHmapper / Pan-Cancer ITH compendium** (bioRxiv Jan 2026) [V] | Preprint | "New pan-cancer ITH package replacing per-sample NMF, with a mapping tool", at more than 1,000 tumours (count [P], not verified) | **High on positioning** | Overlapping continuous programs vs discrete states; no erasure by sample-as-batch; uncertainty; CNA and artefact typing; benchmark against it directly |
| **TCAT/starCAT** (Nature Methods 2025) [V] | Published. The catalogue is of T cells, but the authors present starCAT as generalising to other cell types and tissues [near-verbatim abstract wording], and the README lists further references on Zenodo (which ones not checked; sweep 2 named glioma myeloid and bone marrow [P]) | "Reproducible program catalogue plus projection"; the framework is already pitched beyond T cells, so a malignant-cell reference from the same group is plausible | **Medium-high** for the reference/projection part | Make the reference a deliverable, not the novelty; add per-cell uncertainty, novelty detection, panel transfer |
| **Grabski–Trippa–Parmigiani multi-study NMF** (Genome Biology 2025) [V] plus Tetris/BMSFA line | Published (mutational signatures) | The statistical core (count NMF across studies with inclusion variables and covariates) | **Medium**: if this group ports it to scRNA-seq, the model novelty disappears | Move fast; contribute scRNA-specific pieces (cells within tumours within cohorts, NB, SVI, CNA offsets, power posterior, thinning-based calibration) |
| **scDEF** (bioRxiv; peer-review status unknown) [V] | Preprint | Count, hierarchical, automatic size, shared and patient states in tumours | **Medium** (a pan-cancer extension would be close) | Explicit subset-sharing, nuisance modelling, error control, scale |
| **3CA v2** (Tyler et al., Nature Cancer 6:1088–1101, 2025) [V; code V+] and **unattributed 2026 MP catalogue** (JerryChen06 repo, created 2026-08-03, updated 2026-09-30; my re-check: 70+ MPs, 16 cancer-type groups, EMT-I to EMT-VI, includes "Cell cycle single-nucleus" and "Stress (in vitro)"; no authors) [V+] | Published / unknown | A **resource** scoop (an expanded official catalogue), not a methods scoop | Low (methods), high (baseline moves) | Benchmark against both the 2023 and the latest catalogue; use protocol-named MPs as motivating examples |
| **CanSig benchmark** (Cancer Research 2026) [V] | Published | "First cancer-state discovery benchmark" | **Medium** for a benchmark-only contribution | Plant confounders and nulls; target MPs specifically; add stability, calibration and cross-protocol metrics |
| **DeltaNMF** (bioRxiv 2026) [V] | Preprint | NMF with foundation-model priors; baseline-vs-new program decomposition | Medium for an "FM-informed NMF" or "normal-baseline" angle | Do not pitch FM priors as the novelty; if normal-contrast is used, do it within the multi-tumour hierarchy |
| **bcNMF** (arXiv 2026), **contrastiveVI** (2023) [V] | Preprint / published | Contrastive machinery | Low–medium | Apply to malignant vs matched cell-of-origin inside the MP model |
| **BayesPowerNMF** (PLoS CB 2026), **bayesNMF** (JCGS 2026) [V] | Published | Robust NMF with automatic K (signatures) | Low (different domain) | Cite as foundations |
| **DRVI** (bioRxiv; in scvi-tools), **MOFA-FLEX** (bioRxiv 2025), **Tripso** (bioRxiv 2026) [V] | Preprints | Reviewer requests: "why not DRVI / MOFA-FLEX / FM?" | Low–medium | Benchmark them |
| **Consensus scHPF** (2023), **CSI-GEP** (2025), **KINOMO** (2022/AACR 2024) [V] | Mixed | Automatic K by consensus; replication across splits; partial CNV/patient correction | Low–medium | Cite; KINOMO must be acknowledged on the CNA angle |
| **pyscalop, mpnmf, sciNMF, GeneNMF, CellProgramMapper** [V] | Software | The "faster or deterministic Gavish" and "projection" niches | Saturated | Do not compete there |
| **TCCA** (González-Bermejo…Al-Shahrour, bioRxiv 2026) [V; code V+] | Preprint | Resource built with the standard method: Gavish-style clustering of single-seed per-sample NMF, then rule-guided but hand-coded artefact removal (43 MPs kept by whitelist) | Low (methods); it strengthens the case that curation is manual and irreproducible | Use its patient tumours (minus 3CA-shared studies and cell lines) for replication and its MPs as a baseline; count the 3CA overlap at P0 |
| **ISOMERIC** (bioRxiv 2025) [V] | Preprint | Separating copy-number-dependent from CN-independent tumour-cell variation; the argument that CN dependency confounds shared-state discovery | **Medium** for the CNA/dosage axis (G2) | Type *recurrent cross-tumour* programs against allele-anchored clones with a calibrated test; benchmark ISOMERIC where its code runs |
| **cEBMF** (Stephens lab; arXiv 2505.11639, NeurIPS 2025) [V] | Published (conference) | Covariate-moderated EBMF; known nuisance covariates could enter a GBCD-like factor model directly. Not a tumour-MP method, and no application to tumour MPs was found | Low–medium (a candidate engine rather than a competitor; calling it a "GBCD successor", as the review did, is an interpretation no source states) | Use as a candidate engine inside the model; the novelty remains the recurrence test, error control and provenance typing |
| **SIMVI, DECIPHER** (Nature Communications 2025) [V] | Published | Disentangling intrinsic from spatially induced variation, within spatial datasets | Low (spatial, not recurrent MPs) | Cite as the spatial counterpart of the malignant-intrinsic vs TME-shared split; use on the in situ rung |
| **FDR via data splitting after clustering** (Wang, Lin, Zhao; arXiv 2410.06451) [V] | Preprint | Asymptotic FDR control for features tested after clustering, by split-and-mirror statistics | Low–medium for the "error-controlled gene membership" claim | Cite; adapt to tumour-level cross-fitting of MP membership; extend beyond Gaussian settings and to factorisation |
| **Gene Program Dictionary** (bioRxiv 2025) [P]; **scPFA** (bioRxiv 2026) [P, abstract not retrieved]; **PORCELAN** (Nature Communications 2025) [P] | Preprints / published | Atlas-scale consensus programs with permutation FDR (Dictionary); heritable vs transient programs (scPFA); lineage-informative genes (PORCELAN) | Low–medium; scPFA unknown until read | Read scPFA before claiming novelty on clone/heritability typing |
| **Spectra** (Nature Biotechnology 2023/24) [V] | Published | Motivated partly by technical artefacts in unsupervised factorisation; absorbs known nuisance through priors | Low–medium for "artefact-aware program discovery" | Typing without curated priors, against within-sample negative controls |
| **Protocol and stress studies**: O'Flanagan 2019, Neuschulz 2023, FixNCut 2024, Trinks 2024, the ccRCC fresh/frozen/snFlex preprint 2025 (see §4.1) | Published / preprint | That dissociation stress exists, depends on protocol and can be reduced by cold protease or fixation; protocol comparisons on patient tumours | **Medium** for any stress or protocol-cohort result framed as a discovery | Claim only per-tumour attribution of the stress MP's ex vivo share, which none of them reports measuring as far as checked (full texts were not read) |

### 6.2 White space that remains, given the evidence gathered

No sweep found a 2024–2026 work, published or preprint, that does all of the following:
- (i) defines MPs probabilistically with sharing by a subset of tumours and credible tumour membership;
- (ii) controls error for MP existence and gene membership;
- (iii) models CNA/clone and technical/protocol nuisance inside discovery;
- (iv) ships a malignant-MP reference with a projector that carries uncertainty.

GBCD covers part of (i). ITHmapper and TCAT cover part of (iv). Nothing covers (ii) or (iii).

**Update, October 2026.** About 14 targeted web searches found no published or preprint method that (a) does calibrated or error-controlled MP discovery across tumours, (b) types recurrent tumour programs by nuisance or provenance against negative controls or decoys, or (c) audits or reanalyses the 3CA/Gavish MPs. The nearest work covers parts of the problem: separating CN-dependent variation within tumours (ISOMERIC), disentangling spatial context (SIMVI, DECIPHER), letting covariates moderate a factor model (cEBMF), controlling FDR for tests after clustering (Wang et al., arXiv 2410.06451), and standard-pipeline atlases with rule-guided hand curation (TCCA). This is **weak evidence**: only a search engine was available, without bioRxiv/medRxiv full-text search, PubMed or citing-paper lists, and abstract wording came from snippets.

### 6.3 Mandatory manual scoop check before committing

The sweeps could not reach bioRxiv, Google Scholar or PubMed, and the October 2026 re-check used a search engine only (see §6.2), so this check is still open. Search, 2025–2026:
- "metaprogram robust", "consensus meta-program", "probabilistic metaprogram"
- "pan-cancer NMF uncertainty", "copy number aware gene programs", "clone-aware NMF"
- "multi-study NMF single-cell", "Tetris single-cell"
- citing papers of GBCD, TCAT and Gavish 2023
- any new Stephens, Parmigiani, Raychaudhuri/Kotliar, Boeva, Tirosh or Abelson preprints
- items flagged in Appendix B as found but not read: the scPFA abstract, the full TCCA preprint (3CA overlap), the starCAT abstract sentence, and an unidentified 2026 dissociated-vs-intact tumour preprint

---

## 7. The most promising, defensible novelty gaps

They are ranked by defensibility and impact. G1 to G4 form the recommended scientific core. G5 and G6 are optional biological axes. G7 and G8 are deliverables that strengthen the paper but are risky as stand-alone novelty.

### G1. The MP as a statistical object: count-native, hierarchical, subset-shared program model with credible tumour membership and a data-driven number of MPs

- **Why it matters biologically.** MPs recur in *subsets* of tumours: lineage-restricted, subtype-restricted, or pan-cancer. Saying which tumours express an MP, with confidence, is the biological claim, and it drives downstream associations (subtype, response, survival). It also separates shared programs from patient-private ones (genotype, clone) by construction rather than by threshold.
- **Why existing tools fail.**
  - Gavish, GeneNMF and Barkley rely on thresholds and greedy or forced clustering.
  - GBCD is Gaussian, has a user-set Kmax, gives no credible sets for which tumours share a program, and has not been shown at atlas scale.
  - LIGER/iNMF splits shared from dataset-specific but has no subset sharing, no count likelihood and a hand-picked K.
  - MOFA+ shares loadings across groups.
  - Tetris and BMSFA are Gaussian and do not scale.
  - Grabski multi-study NMF is a mutational-signature model.
  - scDEF is small-scale, with an implicit shared/private split.
- **Evidence that it is unsolved.**
  - No sweep found subset-sharing (IBP or spike-and-slab per program × tumour) count NMF applied to tumour scRNA-seq.
  - Sweep 3 explicitly found no application of Tetris, BMSFA, APAFA, SUFA or multi-study NMF to tumour scRNA-seq.
  - The number of MPs is user-set in every competitor (nMP in GeneNMF, Kmax in GBCD) or falls out of thresholds.
- **Risk.**
  - (a) Engineering: SVI or amortised inference over about 3,000 tumours and about 10^7 cells.
  - (b) Identifiability of the shared/private split (SUFA's warning). An explicit condition is needed (sparsity, anchors, or a scattered condition across tumours).
  - (c) Reviewers will see "another factorisation" unless it is paired with G3 and G4.
  - (d) Overlap with GBCD has to be addressed head-on.
- **De-risking pilot.** Fit a prototype on CanSig datasets plus a 3CA subset (about 200 tumours). Show recovery of planted subset-shared programs where GBCD and Gavish merge or miss them, and coverage of the tumour-membership credible sets.

### G2. Genetic deconfounding inside discovery: CNA dosage and clone as model terms, a positional-enrichment test, and typing of each MP as clone-linked vs non-genetic

- **Why it matters biologically.**
  - Separating heritable subclonal identity from plastic, non-genetic states is central to ITH biology and to therapy resistance (plastic states are reversible and targetable differently).
  - Recurrent arm-level CNAs (e.g., chr7 gain in GBM) can produce positional "programs" that recur across tumours.
  - Dosage-driven programs are real biology, so the goal is to **label** them, not delete them.
- **Why existing tools fail.**
  - CNA is used only to call malignant cells (CanSig, ITHmapper, TCCA) or checked post hoc (Gavish: 24% of subclones linked to a single MP; Kinker: association rates of 39% for clusters vs 8% for programs).
  - sciRED-style GLM nuisance regression exists but has never used CNA.
  - KINOMO and UVG are partial heuristics.
  - CLCC (bioRxiv 2026) aligns expression with CNA but is not program discovery.
  - ISOMERIC (bioRxiv 2025) comes closest: an adversarial autoencoder that splits tumour-cell variation into CN-dependent and CN-independent parts. It does not test cross-tumour recurrence or type MPs. scPFA (bioRxiv 2026, heritable vs transient programs) has not been read.
- **Evidence that it is unsolved.** All 5 sweeps independently list it as a gap. No packaged method was found that discovers recurrent cross-tumour programs jointly with CNA or clone inference and types them (the October 2026 search agrees; weak evidence, §6.2).
- **Risk.**
  - (a) **Circularity**: CNA is inferred from the same expression. Mitigate with allele-aware Numbat, held-out gene sets for CNA inference, and DNA-orthogonal validation (clonealign, CalicoST).
  - (a′) **Expression-derived covariates in decoys or regressions.** A null that holds inferred clones fixed keeps CNA-driven programs by construction (to the extent the inferred clones capture the CNA). Recruitment by a clone-conditioned decoy therefore means "dosage", not "technical nuisance", and recruitment must be reported per covariate. Toy result (`mpaudit.sim`, seed 0): a decoy stratified on estimated clones kept 0.84–0.88 of a planted CNA program's co-expression, against at most 0.04 under a technical (estimated doublet/contamination) decoy; two decoys each.
  - (b) **The effect may be small for continuous programs.** Kinker's 8% figure suggests CNA confounding mostly affects discrete states. **Quantify first.**
  - (c) Real dosage-driven programs must not be over-corrected.
- **De-risking pilot.** On 3CA, compute for all robust programs and MPs: chromosome-arm positional enrichment of top genes, and usage variance explained by Numbat clone. If fewer than about 5% of MPs are affected, reposition G2 as "MP typing" (genetic vs plastic annotation) rather than "removing artefactual MPs".

### G3. Technical-artefact and protocol-aware discovery: explicit nuisance factors plus automatic MP typing, replacing manual curation

- **Why it matters biologically.**
  - "Stress", "Metal", ribosomal/mito and contamination MPs recur in every cancer type.
  - Ex vivo dissociation stress overlaps real in vivo stress (hypoxia, UPR, heat shock near necrosis), so gene-list filters remove both. Baron et al. 2020 argue that a stress-like cancer-cell state is real in vivo, while O'Flanagan 2019 and Neuschulz 2023 show that dissociation induces such a response (§4.1).
  - **Gene-length / nuclear-fraction and abundance / %ribo axes** are not covered by any MP pipeline or QC filter reviewed here. User-reported evidence (HER2 breast-cancer snRNA-seq, unpublished; §1.2): two pre-registered de novo MP searches gave mostly technical axes, namely "luminal" programs of long genes (median 317–437 kb) and a %ribo-tracking abundance axis (29/30 patients, rho 0.68–0.77) that reference matching called "EMT".
  - **Malignant→TME ambient spillover.** In malignant-rich suspensions, a malignant-only program whose strength varies between samples can leak into the same samples' TME cells and look "TME-shared" (the 2026-10-16 review's argument; no published measurement found). Toy result (`mpaudit.sim`, seed 0, 16 samples): adding malignant→TME ambient raised the unadjusted TME partial R² of a CNA-driven malignant-only set from 0.00 to 0.60, and of hypoxia from 0.07 to 0.23. A per-sample spillover adjustment brought them to 0.07 and 0.00. Even unadjusted, both stayed below the toy's matched-random null q99 (about 0.70, from the adjusted run). The hypoxia rise depends on the seed and toy configuration; the real-data size is unknown.
  - Protocol-named MPs persist in current catalogues ("single-nucleus" and "in vitro" MPs in the 2026 catalogue [V+]; sweep 5 reported "Glioma single-nucleus", while my re-check saw "Cell cycle single-nucleus", so the sweeps disagree on the exact name).
- **Why existing tools fail.**
  - Every pipeline ends in manual curation: 3CA v2 `inds_rm` (11 MPs), TCCA's hard-coded whitelist of 43 MPs (criteria only in a code comment), GeneNMF `dropMetaPrograms`, Kinker's and Barkley's manual choices.
  - Spectra is motivated partly by artefacts in unsupervised factorisation, but can absorb nuisance programs only with curated priors and cell-type labels.
  - cEBMF lets known covariates moderate factor priors, but neither tests recurrence nor types programs.
  - None of the pipelines above types a program by gene length, intronic fraction or correlation with %ribo; QC usually stops at UMIs, genes, %mito and doublet scores.
  - sciRED needs known covariates and gives signed factors.
  - SAE technical features are only a hint.
- **Evidence that it is unsolved.** The labelled artefact set exists (the 11 hand-removed v2 MPs plus Gavish's removed LQ, doublet and single-study MPs), and no method has been shown to flag these automatically.
- **Risk.**
  - (a) In vivo vs ex vivo stress needs matched protocol data (fresh vs Flex vs sn vs Xenium), which is scarce. The two vendor datasets are fixed-dissociated vs in situ comparisons: Oliveira 2025 has no fresh-dissociated malignant cells, and whether Janesick 2023 has a fresh-dissociated arm is unverified. The candidate fresh-vs-fixed twins (Trinks 2024 lung; the ccRCC fresh/frozen/snFlex preprint) still need patient matching confirmed, and public twins mostly lack same-protocol replicate pairs for a protocol null.
  - (b) Ambient correction is imperfect, and spillover runs in both directions (malignant→TME as well as TME→malignant).
  - (c) Labels are one lab's judgement.
  - (d) **Expression-derived covariates.** Doublet scores, contamination estimates and inferred clones are computed from expression, so a decoy or regression that holds them fixed also keeps any real program correlated with them (the 2026-10-16 review's example: cycling cells tend to get higher doublet scores). Toy result (`mpaudit.sim`, seed 0): with cycling coupled to the planted doublet score, a decoy conditioned on that score kept about 0.37 of the cycle program's co-expression (per-sample 0.23–0.55), against about 0.01 for the toy's estimated doublet score. That estimate barely tracks cycling in the toy (within-sample Spearman 0.27, against 0.61 for the planted score), so the low value does not show that estimated covariates are safe on real data.

  Mitigation:
  - Model per-sample ambient profiles, doublet scores, depth/complexity, %ribo, intronic fraction, cell cycle and protocol as covariates or fixed factors (Spectra-style anchors, sciRED-style GLM).
  - Type programs with gene-length (median length, fraction of genes > 100 kb), %ribo-correlation and intronic-fraction features. Test TME co-variation after a per-sample malignant-spillover adjustment, against nulls matched on expression and malignant specificity. Report decoy recruitment per covariate arm (technical vs clone/dosage), and plant real programs correlated with nuisance covariates in simulation. (The P1 pilot code, `mpaudit`, implements these on the toy; on real data the length features need a GTF.)
  - Validate typing with leave-one-protocol-out invariance.
  - Report the AUROC for the 3CA v2 hand-removed MPs.

### G4. Error control and calibrated robustness: empirical FDR for MP existence, FDR-controlled gene membership, held-out rank choice, and stability as a reported quantity

- **Why it matters biologically.**
  - MP gene lists are reused everywhere (MSigDB GAVISH_3CA, scoring, clinical correlations). Users need to know which genes and which MPs are trustworthy.
  - Clinical claims need nulls (random signatures associate with outcome; Venet 2011).
- **Why existing tools fail.**
  - No null model for recurrence. The `inter_min=10` threshold guards against a random null with p ≈ 1e-12, which is irrelevant to the structured nulls that actually occur.
  - Top-50 truncation has no significance attached.
  - Stability is assessed at best over restarts (cNMF) or data splits (CSI-GEP), never over combined perturbations and never paired with null calibration.
  - GBCD's lfsr covers genes within a program, not program existence or recurrence.
- **Evidence that it is unsolved.** Sweep 3 found no application of count splitting, data thinning, synthetic nulls (ClusterDE), knockoffs (recall) or selective inference to NMF programs or MP recurrence. Wang, Lin & Zhao (arXiv 2410.06451) control FDR for features tested after clustering, but not for factorisation or cross-tumour recurrence; the October 2026 search found no MP method with error control (weak evidence, §6.2).
- **Risk.**
  - (a) Hivert et al. show thinning loses validity under mixtures with mis-estimated overdispersion. Overdispersion must be estimated conditionally on a preliminary low-rank fit, with calibration shown by simulation. Extending this result to factor models is itself a possible theoretical contribution.
  - (b) How faithful the synthetic null is drives every FDR claim. The null must preserve nuisance covariance (CNA blocks, depth, QC) without preserving real programs.
  - (c) Thinning loses power for small tumours.
  - (d) Stability alone is not publishable (technical programs are stable); it must be combined with null calibration.
- **Built-in identifiability diagnostic.**
  - Report near-pure cells and anchor genes for each program.
  - Optionally use volume or anchor regularisation (Fu/Huang/Sidiropoulos).
  - Argue, with theory and simulation, that pooling heterogeneous tumours restores identifiability of shared programs. No sweep found this argument in the literature.

### G5. Joint within-tumour and between-tumour ("baseline") decomposition without per-tumour centring

- **Why it matters biologically.**
  - Subtype-like states (e.g., PDAC classical vs basal) vary both within and between tumours.
  - Per-tumour centring and zeroing negatives makes programs active across a whole tumour invisible.
  - Treating sample as a batch (ITHmapper's scVI) erases them.
- **Why existing tools fail.**
  - Gavish centres away anything active across a whole tumour.
  - GBCD handles it in part but is per cancer type, not count-based, and not atlas scale.
  - The field appears to be moving toward multi-layer decomposition (sweep 1 reports the Tirosh lab's 2025 GBM "baseline programs" [P]), but no general method formalises it.
- **Evidence that it is unsolved.** GBCD found subtype programs that tumour-by-tumour NMF missed [V]. No atlas-scale count model jointly estimates within-tumour programs, patient baselines, study/chemistry effects and CNA effects.
- **Risk.** **Patient baseline is confounded with batch** whenever one patient equals one library or study. Identifiability needs multi-sample-per-study designs or cross-protocol anchors, and reviewers will push hard on this. This is best built as a layer of G1 rather than a stand-alone claim.

### G6. Contrast with the normal cell of origin: cancer-specific ITH programs vs lineage programs reused from normal tissue

- **Why it matters biologically.** Gavish showed that many carcinoma MPs mirror programs in normal epithelium [V]. Knowing which programs are cancer-specific (candidate targets, drivers of progression) and which are inherited physiological variability changes how MPs are interpreted.
- **Why existing tools fail.**
  - No MP pipeline uses a normal baseline.
  - The machinery exists elsewhere: contrastiveVI (Nature Methods 2023), bcNMF (arXiv 2026), DeltaNMF (baseline vs new programs), and the Human Cell Census NMF reference (bioRxiv 2025) [V].
- **Evidence that it is unsolved.** The factorisation and mp-methods sweeps both list it, and no application was found.
- **Risk.**
  - (a) Matched normal cell-of-origin data are not always available.
  - (b) Adjacent "normal" tissue may be field-cancerised.
  - (c) Contrastive objectives can suppress real programs shared with normal tissue that still matter in tumours. A soft contrast is preferable to subtraction.
  - (d) DeltaNMF and bcNMF reduce novelty if this axis is pitched alone. It works best as an annotation layer ("normal-shared" vs "tumour-enriched" score per MP).

### G7. A malignant-MP reference with uncertainty and a projection tool with novelty detection and panel transfer

- **Why it matters biologically and practically.**
  - Comparable MP usage across studies and platforms (including Flex and Xenium panels) replaces ad hoc module scores.
  - Novelty detection flags programs in a query cohort that the reference lacks.
- **Why existing tools fail.**
  - starCAT's published catalogue is of T cells (the authors present the framework as generalising to other cell types and tissues; we know of no malignant-cell reference, though the Zenodo references were not checked), it is NNLS without uncertainty, and its authors warn about transfer to Flex.
  - ITHmapper maps to discrete states.
  - CellProgramMapper is software without uncertainty.
  - 3CA offers fixed gene lists only.
- **Evidence that it is unsolved.** No sweep found a versioned malignant-MP reference with per-cell uncertainty, QC flags and confidence per gene.
- **Risk.** **High "derivative" risk on its own** (TCAT pattern; the TCAT group or Tirosh lab could release one). Treat it as the deliverable that makes G1–G4 usable, not as the headline.

### G8. An open MP benchmark with a tumour-confounder simulator and an evidence ladder

- **Why it matters.** Without ground truth, "MPs are not robust or real" cannot be measured. A benchmark turns the user's complaint into numbers.
- **Why existing tools fail.**
  - CanSig (Cancer Research 2026) benchmarks state discovery against known signatures. It neither plants confounders and nulls nor studies seed or rank stability.
  - No published simulator jointly models segment CNAs plus subclones, ambient RNA in both directions, malignant–TME doublets, cell-varying stress, gene-length / abundance axes and platform effects.
  - PerturbNMF metrics are not tumour-specific.
- **Evidence that it is unsolved.** All 5 sweeps flag it. Sweep 4 found no dedicated benchmark of gene-program methods on tumour scRNA-seq; CanSig is the nearest, and that sweep's search was limited.
- **Risk.**
  - A benchmark alone fits Genome Biology or a Nature Methods "Analysis" article better than a methods paper.
  - Semi-synthetic realism will be scrutinised.

  Mitigation: anchor confounders in measured quantities from the same tumours (Numbat segments, empty-droplet profiles, real TME doublets), and pair with real-data evidence (cross-protocol, perturbation).

### Recommended composite claim and what to avoid

- **Headline (G1 + G2 + G3 + G4).** A count-native hierarchical model of shared, subset-shared, tumour-private and cohort-private programs, with:
  - CNA/clone and technical/protocol nuisance terms;
  - power-posterior robustness;
  - K, ζ and the number of MPs chosen by held-out (count-split) likelihood;
  - empirical FDR per MP from nuisance-preserving null cohorts;
  - FDR-controlled gene membership;
  - credible tumour membership;
  - a per-MP type label (shared biology / clone-linked / protocol-sensitive / artefact) plus a stability index.
- **Flagship biology.** Re-derive 3CA. Show which canonical MPs are unstable, artefactual, clone-linked or protocol-dependent, and which shared programs were merged or missed. Validate with cross-protocol, perturbation, spatial and clone-orthogonal evidence. Ideally include one experimental follow-up.
- **Deliverables (G7 + G8, with G5 and G6 as layers).** A versioned malignant-MP reference plus projector, and an open benchmark and simulator.
- **Do NOT pitch any of these as the main novelty:**
  - another per-tumour NMF variant, or a faster or deterministic Gavish (saturated);
  - foundation-model-informed NMF (DeltaNMF exists; FM embeddings underperform baselines zero-shot);
  - reference plus projection (TCAT pattern);
  - "shared vs patient-specific" (GBCD).

### Pilot analyses to run before committing (cheap, decisive)

1. **Instability audit of the standard.** Run the Gavish code on 3CA subsets varying seed, rank grid, the 35/10/10/5 thresholds, gene universe (top-7,000 vs HVG) and leave-one-study-out. Quantify how MP identity, gene lists and MP count move. If the top MPs prove highly stable, refocus the pitch on calibration, typing and edge MPs.
2. **CNA and positional audit.** Arm-level enrichment of robust programs and MPs; variance in usage explained by clone (Numbat), per tumour.
3. **Artefact predictability.** Can simple nuisance covariates (ambient profile similarity, doublet score, %mito/ribo, correlation with %ribo and intronic fraction, median gene length and fraction of genes > 100 kb, sc vs sn) separate the 11 hand-removed 3CA v2 MPs from retained ones? TCCA's whitelist (43 kept of at least 57) gives a second, independent set of curation decisions to test against.
4. **Effect of centring.** Compare centred and uncentred per-tumour NMF against GBCD on CanSig datasets to size G5.

Pilots 1–3 are implemented, in reduced form, as the P1 audit (`mpaudit`; see `docs/05_pilot_P1.md`). It runs a Python port of the 3CA code, not the R code itself; R-versus-port parity on real samples has not been tested yet (R was unavailable; the comparison scripts are in `scripts/parity/`). Its GO / PIVOT / INCOMPLETE decision follows the frozen rule **P1-GO v1.0**. The single source of truth is `src/mpaudit/report.py` (module docstring and `rule_text()`), and `docs/02` and `docs/05` state the same rule. It is deliberately not restated here, so that no fourth version can drift.

---

## Appendix A. Where the sweeps disagreed, and how this was resolved

| Topic | Disagreement | Resolution |
|---|---|---|
| 3CA v2 code | Sweep 4: the code building v2 MPs is not in the repo (MPs loaded precomputed). Sweeps 1 and 5: it is in `MP_analysis_new_batch.R`. | **Resolved in favour of sweeps 1 and 5.** I fetched the script, which contains ranks 4–9, nrun 10, top 7,000, centring and zeroing, 35/10/10/5, the single-study filter and `inds_rm` with annotations [V+]. |
| Barkley 2022 | Sweep 4 marked it unverified. | Verified by sweeps 1 and 5 (Nature Genetics 54:1192–1201; DOI 10.1038/s41588-022-01141-9). Not cNMF (sweep 1, from code). |
| Gavish MP count | Sweep 1: 41 [V]. Sweep 5: "recalled 41, not re-checked". | Use 41, per sweep 1. |
| GBCD authors and K | Sweep 3 listed 3 authors. Sweeps 1/2 differ on Kmax ("user-set" vs "EB-adaptive with a user maximum"). | The README lists 6 authors [V+]. The two K descriptions agree: empirical-Bayes shrinkage under a user-set maximum. |
| GBCD comparison | Sweep 1: beats *tumour-by-tumour* NMF. Sweeps 2/3: *pooled* NMF finds mostly patient-specific programs. | Probably both comparisons are in the paper; full text not accessible. Cite both cautiously. |
| Spectra year | 2023 vs 2024. | Online Sept 2023; print 2024, Nature Biotechnology 42:1084–1095. |
| MOFA+ | Sweep 2 verified; sweep 3 did not. | Treat as verified (Genome Biology 21:111, 2020). |
| CalicoST venue | Sweep 4 believed Nature Methods. | README cites bioRxiv 2024 [V+]. Journal status unconfirmed; cite as a preprint. |
| scMultiSim venue | Sweep 4 believed Nature Methods. | README cites bioRxiv 2023 [V+]. Journal status unconfirmed. |
| UCE venue | "Nature" snippet vs s41588 URL. | Likely Nature Genetics [P]; peripheral. |
| 2026 MP catalogue MP name | Sweep 5: "Glioma single-nucleus". | My re-check of `programs_per_cancer_type.html` showed "Cell cycle single-nucleus" and "Stress (in vitro)" [V+]. The page may be truncated, so treat the exact name as uncertain. The point stands: protocol-named MPs persist. |
| ITHmapper coverage | Sweep 1: 14 cancer types, more than 1,000 tumours. | The mapping tool README lists 12 supported cancer types. The manuscript repo has a placeholder for the preprint link; data on Zenodo 10.5281/zenodo.16852297 [V+]. |
| Existence of an MP benchmark | Sweep 4: none found (search unavailable). Sweep 1: CanSig benchmark (Cancer Research 2026). | CanSig exists [V by sweep 1] but is a state-discovery benchmark on known signatures, without planted confounders, nulls or seed/rank stability. The MP-specific gap stands, narrowed. |

## Appendix B. Items not verified, or only partly verified (do not cite without checking)

- **Moved out of this list on 2026-10-17** (now [V] in the main text, with content): the dissociation papers van den Brink 2017, O'Flanagan 2019, Denisenko 2020 and Slyper 2020 (the last two at title level only). Also verified and added: Marsh 2022, Baron 2020, Neuschulz 2023, FixNCut 2024, TCCA (code read on GitHub), ISOMERIC, cEBMF, SIMVI, DECIPHER, and arXiv 2410.06451. Open details on some of them are listed below.
- **[U]**
  - scPFA (bioRxiv 2026): title only, abstract not retrieved.
  - An unidentified 2026 bioRxiv preprint comparing dissociated with intact tumour (HSF1, RELA, NFKB1 and JUN reportedly activated by dissociation), seen only as a fragmentary search snippet; title not found. If real, it pre-empts part of the stress-decomposition result. Find it before P0.
- **[P]**
  - Tirosh lab 2025 GBM "baseline programs" paper (Nature Genetics 57:1155–1167; single sweep).
  - Glioma myeloid cNMF paper (Nature 2025; no citation retrieved).
  - TCAT CITE-seq and activation validation details.
  - starCAT: the exact abstract sentence on generalising to other cell types and tissues (near-verbatim; read the article before quoting), and which references exist on Zenodo.
  - Greenwald GBM venue (Cell 2024).
  - KP-Tracer venue (Cell 2022).
  - Replogle 2022 venue.
  - Janesick 2023 venue, and whether it has a fresh-dissociated arm from the same tumour.
  - Oliveira 2025: the material used for the immune-profiling arm. The review's sub-claim that "the only fresh cells are sorted T cells" is unverified; what is verified is that there are no fresh-dissociated malignant cells.
  - O'Flanagan 2019 patient tumour types; FixNCut tumour content and size of the stress reduction.
  - Trinks 2024 patient count and Flex use (the stress finding comes from an interview, not the abstract).
  - ccRCC fresh/frozen/snFlex preprint (Xirenayi et al. 2025): patient matching, patient numbers, any stress or malignant-program result.
  - 10x GEM-X Flex performance preprint (2026) and Frontiers 2026 matched FFPE/frozen brief report: snippets only.
  - TCCA: overlap with 3CA studies (count from GEO accessions at P0); patient count 537 vs 540 by source.
  - ISOMERIC authors and publication status; Gene Program Dictionary and PORCELAN details; Spectra's "prone to technical artifacts" wording (snippet).
  - ITHmapper cell, tumour and cancer-type counts (the 2026-10-16 review could not verify them, and a 2026-10 search did not find the preprint).
  - Saelens 2018.
  - Meinshausen & Bühlmann 2010.
  - Fu/Huang/Sidiropoulos journal venues.
  - Generalised data thinning JASA venue.
  - Chen & Witten and Yun & Barber author lists.
  - BLAST and CAVI arXiv authors.
  - Probit multi-study NMF authors.
  - DeltaNMF RECOMB acceptance.
  - Tripso authors and initialisation.
  - GeneNMF-citing BCC paper's J Invest Dermatol 2025 publication.
  - 3CA v2 "149 MPs" (spans 8 cell types).
  - scDEF peer-review status.
  - Spectra and expiMap likelihoods.
  - GBCD practical maximum scale.
- **Hypergeometric figures** in Sections 1.2 and 1.4 are my own calculation (top-50 lists; 7,000- and 2,000-gene universes), not a published result.
- **User-reported evidence** (§1.2, G3): the HER2 breast-cancer snRNA-seq results (long-gene "luminal" programs, a %ribo-tracking "EMT" abundance axis) are the user's own unpublished, pre-registered analyses. They are not reproduced here and are not a citation.
- **Toy results** (G2, G3): the numbers from `mpaudit.sim` come from a small synthetic cohort built to test the audit code. They show that a mechanism can occur, not how large it is in real data.
- **Search limitation.** This synthesis pass had no WebSearch budget left (200/200 used). Publisher sites, bioRxiv, arXiv, PubMed, Crossref and OpenAlex were proxy-blocked. Only GitHub content could be re-checked. Coverage of 2025–2026 preprints is therefore incomplete (Section 6.3). In the 2026-10-16/17 re-check, WebSearch worked but publisher, PMC, bioRxiv, arXiv and Europe PMC pages still did not resolve. Abstract wording therefore comes from search snippets and is not guaranteed verbatim, except where an exact-phrase search matched. GitHub content (TCCA code, starCAT and Oliveira READMEs) came through a summarising fetcher, but the code lines quoted from it are verbatim.