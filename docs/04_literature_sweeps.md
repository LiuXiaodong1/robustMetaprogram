# Raw literature sweeps

> Raw per-item output from five parallel literature sweeps (2026-10-02), the input to [01_landscape_and_gap_analysis.md](01_landscape_and_gap_analysis.md). `verified_online` means the sweep found the item online during the search; it does not guarantee every detail. Check before citing.


## Direct competitors: multi-tumour metaprogram methods (24 items)

### Gavish et al. 2023 meta-programs (MPs) / 3CA ITH_hallmarks pipeline

*Gavish A, Tyler M, Greenwald AC, ... Tirosh I. Nature 618:598-606 (2023). Code: github.com/tiroshlab/3ca (ITH_hallmarks)* — ✅ verified online — <https://www.nature.com/articles/s41586-023-06130-4>

- **Core idea:** 77 studies, 1,163 tumours, 24 tumour types gave 41 malignant MPs and 11 'hallmarks' of transcriptional ITH. The steps below come from the repo code. (1) Per sample, take the top 7,000 genes by mean expression, centre each gene and set negative values to 0. Run NMF (R 'NMF' package, method snmf/r) at ranks K=4..9, nrun=10, which gives 39 programs per tumour. (2) Each program is reduced to its top 50 genes. (3) Robust-program filter (robust_nmf_programs): intra_min=35 means a program must share at least 35/50 genes (70%) with a program from another rank in the same tumour. inter_min=10 means it must share at least 10/50 genes (20%) with a program in another tumour. intra_max=10 removes redundancy: programs are ranked by cross-tumour similarity and greedily kept, and any program sharing more than 10 genes with one already kept is dropped. (4) Programs are clustered greedily and iteratively into MPs by gene overlap. The founder is the program with the most partners at 10 or more shared genes (Min_intersect_initial=10). New programs are added while they share at least 10 genes with the evolving MP gene list (Min_intersect_cluster=10), and a cluster needs at least 5 programs (Min_group_size=5). The MP is the top 50 genes by frequency among member programs, with ties broken by NMF score. (5) MPs are then removed after the fact if they look low-quality (enriched in ribosomal-protein or mitochondrial genes), come from a single study, or look like doublets (resembling T cells or macrophages).
- **Strengths:** Very large curated scale. Simple and interpretable. Code is public. MPs are widely reused (e.g., MSigDB GAVISH_3CA gene sets). Malignant MPs were compared with normal-epithelium MPs and with Kinker cell-line RHPs, showing that carcinoma MPs often mirror programs in non-malignant epithelium.
- **Limitations:** Rank: a fixed grid (4-9). 'Robustness' only means recurring across neighbouring ranks; no rank is chosen per tumour. Seed: nrun=10 in the NMF package keeps the best fit, not a consensus, and seed variance is never reported. Thresholds: 35/10/10/5, the top-50 list length and the 7,000-gene universe are all arbitrary, and binarised gene lists throw away the loadings. Clustering is greedy and order-dependent, with no probabilistic model. Technical artefacts are removed by manual post-hoc curation. CNV is only checked afterwards (only 24% of genomic subclones were associated with a single MP), and CNV-driven programs are not removed. Recurrent CNAs (e.g., chr7 gain in GBM) could create recurrent 'positional' programs, and no positional-enrichment test is run. Centring plus zeroing negatives keeps only relative upregulation within each tumour, so programs present tumour-wide are lost. MP membership has no uncertainty or statistical significance attached.
- **Relevance:** This is the de facto standard the user wants to replace, so it must be the primary baseline. Every threshold and the manual-curation step are concrete targets to improve on. The robustness, CNV and validation analyses in the paper are all post hoc.

### Curated Cancer Cell Atlas v2 (3CA) updated MPs

*Tyler M, ... Tirosh I. Nature Cancer 6:1088-1101 (2025); bioRxiv 2024.10.11.617836. Code: tiroshlab/3ca/3CA_v2* — ✅ verified online — <https://www.nature.com/articles/s43018-025-00957-8>

- **Core idea:** 124 datasets, more than 40 cancer types and 2,836 samples. The preprint and press summaries report 149 MPs across 8 cell types. The Gavish pipeline was re-run on the expanded cohort with the same settings according to MP_analysis_new_batch.R: ranks 4-9, nrun 10, top 7,000 genes, centring and zeroing negatives, top-50 genes, overlap of 10 or more genes. New filters: MPs drawn from a single study are dropped, and so are MPs with fewer than 5 programs. Then 11 MPs are removed by hand. The code line is inds_rm <- c(5,7,11,17,25,30,36,44,58,68,74), annotated 'LQ, LQ, Macro, T-cells, Endo, LQ, LQ, LQ, miss-annotation, LQ, LQ'. New MPs were matched to old ones when Jaccard was above 0.3 and every other overlap was below 0.2. Single-cell and single-nucleus versions were compared for some MPs (Cilia, Cell cycle).
- **Strengths:** Largest curated, harmonised pan-cancer scRNA-seq resource, free to download. Includes non-malignant cell types. The web portal lets users compare their own gene sets with the MPs.
- **Limitations:** The method is unchanged from 2023. Seven MPs were hand-flagged as low quality, three as contaminating non-malignant cell types (macrophage, T cell, endothelial) and one as mis-annotated. This shows the pipeline routinely produces artefactual MPs that need expert curation. Study, chemistry (sc vs sn) and dissociation effects are not modelled explicitly. The atlas offers gene lists and scores, not a projection model with uncertainty.
- **Relevance:** This is the latest and largest instance of the standard, and the natural benchmark resource (data plus MP lists). The 11 hand-removed MPs make a ready-made 'artefact' test set for any new method that claims automatic QC.

### Kinker et al. 2020 recurrent heterogeneous programs (RHPs) in cancer cell lines

*Kinker GS, Greenwald AC, Tal R, Orlova Z, Cuoco MS, McFarland JM, ... Tirosh I. Nature Genetics 52:1208-1218 (2020). Code: github.com/gabrielakinker/CCLE_heterogeneity* — ✅ verified online — <https://www.nature.com/articles/s41588-020-00726-6>

- **Core idea:** Multiplexed scRNA-seq of 198 cell lines from 22 cancer types. Per cell line: NMF at k=6-9, top-50 genes, then the same robust_nmf_programs filter (35/10/10), giving 1,445 robust programs. These were hierarchically clustered (average linkage, distance = 50 - number of overlapping genes), using only programs with at least 10/50 overlap in another line. A DBSCAN discrete-cluster analysis ran in parallel. In module2_rhp.R, programs correlated with cell complexity and cell-cycle clusters are removed by hand, and metaclusters are defined by hand (genes present in at least 25% of member programs). The result is 12 RHPs.
- **Strengths:** Clean system with no TME or dissociation confounders. Functional follow-up included senescence-like subpopulations and their drug sensitivities. Programs were compared with tumour programs.
- **Limitations:** Fixed rank grid, single-run seeds and threshold heuristics, with metaclusters curated by hand. The genetic analysis is post hoc. Cell lines have no microenvironment, dissociation or patient effects.
- **Relevance:** This is the direct ancestor of the Gavish procedure. It gives an empirical estimate of CNV confounding: 39% of discrete clusters but only 8% of continuous NMF programs were significantly associated with CNA-defined subclones. A new method could turn this into a built-in test rather than a post-hoc check.

### Barkley et al. 2022 recurrent cancer cell states (consensus modules)

*Barkley D, Moncada R, Pour M, Liberman DA, ... Yanai I. Nature Genetics 54:1192-1201 (2022). Code: github.com/yanailab/PanCancer* — ✅ verified online — <https://www.nature.com/articles/s41588-022-01141-9>

- **Core idea:** About 20,000 malignant cells from 62 untreated primary tumours across 15 cancer types gave 16 consensus modules (stress, IFN, EMT, metal response, basal, ciliated, etc.). From the repo code: per-tumour NMF with the NMF package, method nsNMF, seed='ica' (deterministic), nrun=1, ranks 2-25. Each gene is assigned to the factor where its normalised score ranks highest, and factors with fewer than 5 genes are dropped (gmin=5). The rank chosen minimises (rank - number of modules), with ties going to the higher rank. Across tumours, modules are filtered iteratively (Jaccard above 0.05 with at least 3 other modules, 5 iterations). A gene co-occurrence graph is then clustered with infomap (100 trials). The final modules are chosen and renamed by hand in Finding.R. Correction to the task framing: this is not cNMF.
- **Strengths:** Validated with Visium spatial transcriptomics and mouse models: the interferon module depends on tumour location and drops when lymphocytes are depleted. Links states to TME interactions.
- **Limitations:** Single NMF run per rank with ICA initialisation, so it is deterministic but not tested for robustness. The rank heuristic is ad hoc. Final module selection is manual. CNV and technical confounders are not modelled. The cohort is small (62 tumours).
- **Relevance:** An alternative to Gavish-style clustering that uses a graph of gene co-occurrence plus community detection. It shows the field has no agreed cross-tumour aggregation step.

### GeneNMF (multi-sample NMF meta-programs, R/Seurat)

*Andreatta M, Carmona SJ (Univ. Geneva). CRAN package GeneNMF (GitHub DESCRIPTION v0.9.6). No standalone methods paper found; the authors ask users to cite Yerly L, Andreatta M, ... Carmona SJ, Kuonen F, 'Wounding triggers invasive progression in human basal cell carcinoma', bioRxiv 2024 (reported as published in J Invest Dermatol 2025, which I did not check against the journal).* — ✅ verified online — <https://github.com/carmonalab/GeneNMF>

- **Core idea:** multiNMF runs NMF per sample with RcppML over a k range (default k=5:6; the cancer demo uses 4:9). It uses 2,000 HVGs with an optional hvg.blocklist (e.g., mitochondrial or ribosomal genes) and min.exp/max.exp filters, does not centre or scale by default (since v0.6), and makes one run per k with a fixed seed=123. getMetaPrograms computes cosine (or Jaccard) similarity between program gene-weight vectors, clusters with hclust ward.D2, and cuts the tree into a user-specified nMP (default 10). MP genes must appear in at least min.confidence (0.5) of member programs and fall within weight.explained (0.5), with specificity.weight=5 and at most 200 genes. Diagnostics are sample coverage, silhouette and mean similarity. dropMetaPrograms() removes MPs by hand; the demo drops an MP with negative silhouette. The authors say per-sample NMF 'bypasses the need to perform batch effect correction'.
- **Strengths:** Fast (RcppML), easy Seurat integration, on CRAN, uses gene weights through cosine similarity rather than binarised top-50 lists, and gives useful MP-level diagnostics.
- **Limitations:** The user must pick the number of MPs (nMP), and there is no data-driven choice. A single fixed seed makes results reproducible but hides seed variance. Ward hclust plus cutree forces every program into some MP. There is no CNV, doublet or stress handling beyond an optional gene blocklist, and MPs are dropped by hand. The CanSig benchmark (Cancer Res 2026) found GeneNMF signatures less rediscoverable across datasets than those from the best early-integration methods.
- **Relevance:** The most direct packaged competitor and the most-used user-facing implementation of the Gavish idea. A new package must clearly beat it on stability, rank and nMP selection, and artefact rejection.

### CanSig: benchmark for malignant shared-state discovery (and the original CanSig method)

*Boeva lab (ETH Zurich). Benchmark paper: 'CanSig Benchmarks Methods for Reproducible Cancer Cell State Discovery from Single-Cell Transcriptomic Data', Cancer Research 86(4):873 (Feb 2026), DOI 10.1158/0008-5472.CAN-25-0940. Earlier method preprints: bioRxiv 10.1101/2022.04.14.488324 (v1-v3).* — ✅ verified online — <https://doi.org/10.1158/0008-5472.CAN-25-0940>

- **Core idea:** The benchmark scores 13 methods on 12 scRNA-seq datasets from 5 cancer types (GBM, breast, LUAD, rhabdomyosarcoma, cSCC), covering 185 patients and 174k malignant cells. Eleven are 'early integration' methods (BBKNN, CCA, Combat, Dhaka, Harmony, ldvae, fastMNN, MNC, rPCA, Scanorama, scVI). Two are 'late integration' methods that build per-patient programs and then aggregate them: scalop (the Tirosh-lab R toolkit) and GeneNMF. Metrics combine batch-correction and bio-conservation scores with a signature-correlation metric, and cover rediscovery of known signatures, cross-dataset reproducibility and clinical relevance. Harmony, BBKNN and fastMNN scored best. Signatures from CCA, scalop and GeneNMF were less rediscoverable. Late-integration methods improved as more samples were added. The original CanSig pipeline calls CNVs (infercnvpy) to pick malignant cells, integrates patients, clusters cells, derives DE signatures and sets them side by side with CNV profiles.
- **Strengths:** Standardised, multi-cancer and open source. It brings clinical relevance and cross-dataset reproducibility into evaluation.
- **Limitations:** Ground truth is limited to known signatures in a few cancer types. Rank and seed stability of NMF pipelines are not the focus. CNV is handled only for malignant calling and side-by-side comparison, not as a confounder that gets removed. The Gavish/3CA code itself was not benchmarked directly (scalop is related but different).
- **Relevance:** The only cancer-specific benchmark harness, and it already shows that MP-style late integration underperforms on reproducibility. A Nature Methods submission should reuse or extend its datasets and metrics and beat both groups of methods.

### GBCD (generalized binary covariance decomposition)

*Liu Y, Carbonetto P, Willwerscheid J, Oakes SA, Macleod KF, Stephens M. Nature Genetics 57:263-273 (2025); bioRxiv 2023.08.15.553436. R package stephenslab/gbcd* — ✅ verified online — <https://www.nature.com/articles/s41588-024-01997-z>

- **Core idea:** Fits the cell-by-cell covariance of pooled log expression across all tumours as LL^T. It uses empirical Bayes matrix factorisation (flashier/ebnm) with 'generalized binary' priors on membership L and keeps programs near-orthogonal. This stops shared programs from being absorbed into patient-specific ones. Gene signatures are log-fold changes with posterior uncertainty (lfsr). Without using patient or study labels, it splits variation into patient-specific, dataset-specific and shared programs. In HNSCC and three PDAC studies it recovered subtype programs missed by tumour-by-tumour NMF, plus a stress-related program prognostic of PDAC survival independent of stage and subtype.
- **Strengths:** Principled probabilistic model, gene-level uncertainty, no reliance on per-sample NMF thresholds. Shows a direct win over tumour-by-tumour NMF.
- **Limitations:** Kmax is still user-set (the authors suggest 10-40; they used 24 for HNSCC and 32 for PDAC). The cell-by-cell covariance limits scale, and the method has not been shown at around 1,000-tumour atlas scale. It uses a Gaussian model on log counts. CNV-driven and technical variation is not modelled explicitly; it simply ends up in 'patient-specific' components. Each analysis covers a single cancer type.
- **Relevance:** The strongest statistical competitor: a joint, label-free model that tackles the intra- versus inter-tumour confound and gives gene-level uncertainty, both of which the MP pipeline lacks. A new method should be compared against it and ideally combine its strengths with atlas-scale, CNV-aware and QC-aware modelling.

### TCAT / starCAT (multi-dataset consensus GEP catalog plus projection)

*Kotliar D, Curtis M, Agnew R, ... Raychaudhuri S. Nature Methods 22:1964-1980 (2025). github.com/immunogenomics/starCAT* — ✅ verified online — <https://www.nature.com/articles/s41592-025-02793-1>

- **Core idea:** cNMF was run separately on 7 T-cell datasets (1.7M cells, about 700 individuals, 38 tissues). Each dataset was batch-corrected with an adapted Harmony that outputs nonnegative gene-level data. GEPs were then clustered across datasets and each cluster averaged into a consensus GEP: 46 cGEPs, with 52 GEPs in the TCAT.V1 reference. Nine cGEPs appeared in all 7 datasets (mean R=0.81), and 68-96.8% of GEPs had a match in another dataset. Of 52 singleton GEPs, 49 were excluded as likely dataset artefacts and 3 kept with a biological justification. starCAT projects the fixed catalog onto new data to estimate per-cell usages and derived scores.
- **Strengths:** Rigorous reproducibility analysis, a reusable reference with versioning on Zenodo, and Python and R tooling.
- **Limitations:** T cells only; no malignant catalog. Artefacts are still rejected by a reproducibility heuristic plus manual judgement. Uses cNMF's K selection per dataset. Reliability on targeted panels (e.g., 10x Flex) is untested.
- **Relevance:** The key Nature Methods precedent and template: cross-dataset reproducibility as the main criterion, a fixed reference catalog and a projection tool. A malignant-cell analogue with uncertainty and automatic artefact and CNV handling would be a natural positioning, but it must show more than just TCAT applied to tumours.

### cNMF (consensus NMF)

*Kotliar D, Veres A, Nagy MA, ... Sabeti PC. eLife 8:e43803 (2019). github.com/dylkot/cNMF* — ✅ verified online — <https://elifesciences.org/articles/43803>

- **Core idea:** Runs NMF many times on one dataset, pools the components, removes outliers by local-density filtering, clusters them and takes cluster medians as consensus GEPs. Per-cell usages are then refit. K is picked from a stability (silhouette) versus reconstruction-error plot. It outperformed other factorisations in simulations at recovering identity and activity programs.
- **Strengths:** Well validated, widely used, and gives interpretable spectra and usages.
- **Limitations:** Single-dataset method. K is still picked by eye from a plot. Nothing aggregates programs across samples. Batch, CNV or technical confounders are not modelled.
- **Relevance:** The usual per-sample engine for the 'robust' variant of MP pipelines (it fixes seed instability within a sample). Examples include the glioma myeloid programs in Nature 2025 (cNMF per cohort, then hierarchical clustering into 14 consensus programs) and TCAT. A new package could adopt consensus over seeds as one layer of its uncertainty model.

### Pan-Cancer Single-Cell Compendium of ITH + ITHmapper

*Nofech-Mozes I, Cook A, Ouellette TW, Hafezi-Bakhtiari S, Awadalla P, Abelson S. bioRxiv 2026 (10.64898/2026.01.06.693992). Code: github.com/abelson-lab/ITHmapper and ITHmapper-manuscript* — ✅ verified online — <https://www.biorxiv.org/content/10.64898/2026.01.06.693992v1>

- **Core idea:** Covers 3.6M cells from more than 1,000 primary and metastatic tumours across 14 cancer types. From the manuscript scripts: per cancer type, scVI is trained on malignant cells with batch_key='sample', 10,000 HVGs and 10 seeds. Hotspot autocorrelation gene modules are computed on X_scVI for each seed (n_neighbors=30, FDR<0.05, minimum module size 20). Cells are clustered at multiple resolutions and filtered by silhouette. Meta-clusters are merged using correlated deconvolved fractions, and a pan-cancer network of states is built with GSEA, averaging about 18.5 states per cancer type. ITHmapper (Python) maps query cells to reference states using multi-seed module scoring, a neighbour graph and a silhouette cutoff of 0.2. Validation covers survival associations, Visium HD spatial data, metastasis enrichment, cell lines (Kinker) and Tahoe-100M drug perturbations.
- **Strengths:** Scale, multi-seed design, spatial and drug-response validation, and a query-mapping tool.
- **Limitations:** Preprint. States are discrete clusters rather than continuous, overlapping programs. Treating sample as the scVI batch can erase programs that are present tumour-wide. CopyKAT is used for malignant calling, not to remove CNV-driven modules. Seed and resolution choices are still heuristic.
- **Relevance:** A very recent large-scale competitor (2026) that replaces per-sample NMF with integration plus Hotspot and ships a mapping tool. It overlaps directly with a 'new package' pitch and must be cited and benchmarked.

### MetaTiME (meta-components across many tumour datasets)

*Zhang Y et al. Nature Communications 14:2634 (2023). github.com/yi-zhang/MetaTiME. The full author list and senior author were not checked in this session.* — ✅ verified online — <https://www.nature.com/articles/s41467-023-38333-8>

- **Core idea:** Runs ICA on the log expression of each of 79 tumour scRNA-seq datasets (1.7M TME cells). It keeps components that are reproducible across cohorts (maximum Pearson r of at least 0.3 with an IC from another dataset), clusters them with Louvain and averages each cluster into 86 annotated meta-components (MeCs). A tool scores and annotates new data.
- **Strengths:** Large scale, reusable reference and projection tool, and correlation on weights rather than binarised gene sets.
- **Limitations:** Covers non-malignant TME only. Uses ICA with signed components and fixed correlation thresholds. No explicit treatment of seed or rank uncertainty or of technical confounders.
- **Relevance:** A program-meta-clustering analogue that uses continuous correlation of component weights and graph clustering instead of top-50 Jaccard, applied across many datasets. A useful baseline and source of design ideas.

### mosaicMPI (multi-resolution program network across cohorts)

*Verhey TB, Seo H, Gillmor A, Thoppey-Manoharan V, Schriemer D, Morrissy S. Nucleic Acids Research 52(12):e53 (2024). github.com/MorrissyLab/mosaicMPI* — ✅ verified online — <https://academic.oup.com/nar/article/52/12/e53/7685157>

- **Core idea:** Runs cNMF on each dataset across many k values (multi-resolution), without needing shared features or HVG intersection. Programs from all resolutions and datasets are linked in a correlation network, community detection groups them into 'themes', and clinical and molecular annotations are statistically propagated onto the network. It integrates single-cell and bulk data across platforms (GBM demo).
- **Strengths:** Flexible, cross-modal and incremental, with annotation propagation.
- **Limitations:** Correlation and community thresholds are heuristic, and there is no probabilistic uncertainty. CNV and technical artefacts are not modelled. It has not been benchmarked against the Gavish pipeline on malignant ITH.
- **Relevance:** Handles rank by keeping all resolutions and letting the network merge them, which is a direct conceptual alternative to the 'robust across K' filter. It supports incremental addition of cohorts.

### EcoTyper (cell states and ecotypes across tumours)

*Luca BA, Steen CB, ... Gentles AJ, Newman AM. Cell (2021), 'Atlas of clinically distinct cell states and ecosystems across human solid tumors'. github.com/digitalcytometry/EcoTyper* — ✅ verified online — <https://github.com/digitalcytometry/EcoTyper>

- **Core idea:** Purifies cell-type-specific profiles with CIBERSORTx. For each cell type it runs NMF over ranks 2-20 and picks the rank from the cophenetic coefficient, with a heuristic for cases where the coefficient curve has several peaks. Quality control drops states with fewer than 10 markers and states whose Adjusted False Positive Index (AFI) is 1 or higher. It defines multicellular ecotypes and recovers states in bulk, scRNA-seq and spatial data.
- **Strengths:** Automated QC, clinical associations, and transfer across modalities.
- **Limitations:** Mostly a bulk-deconvolution setting. States are discrete. It predates the 2023-2026 window.
- **Relevance:** An example of built-in automatic false-positive filtering of factorisation-derived states (AFI) and of automated rank selection. Both are missing from MP pipelines and could be adapted.

### Spectra (supervised, prior-guided gene programs)

*Kunes RZ, Walle T, ... Pe'er D. Nature Biotechnology (published 21 Sep 2023)* — ✅ verified online — <https://www.nature.com/articles/s41587-023-01940-3>

- **Core idea:** Factorisation that combines user-supplied gene sets (as a gene-gene knowledge graph) and cell-type labels with discovery of new programs. It learns both global and cell-type-specific factors, and a parameter (eta) measures how far each factor moved from its prior. It outperformed NMF and scHPF at finding perturbation-related programs in the correct cell type.
- **Strengths:** Interpretable, handles cell-type-specific programs, and shows strong results in tumour-immune settings.
- **Limitations:** Needs gene-set priors and cell-type labels. It does not model recurrence across many tumours or patient effects. Results depend on the input gene sets.
- **Relevance:** Offers a way to explicitly 'absorb' known nuisance programs (dissociation stress, cell cycle, ribosomal/mitochondrial, contamination signatures) into labelled factors so they do not contaminate MPs, and to anchor interpretation.

### Joint shared-plus-specific NMF: LIGER/iNMF (and consensus iNMF), UINMF, scINSIGHT

*Welch JD et al., Cell 2019 (LIGER/iNMF); Kriebel AR & Welch JD, Nature Communications 13:780 (2022) (UINMF); Qian K, Fu S, Li H, Li WV, Genome Biology 2022 (scINSIGHT)* — ✅ verified online — <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02649-3>

- **Core idea:** iNMF factorises all datasets jointly into shared metagenes plus dataset-specific metagenes. rliger also offers a consensus iNMF (runCINMF) over random starts. scINSIGHT jointly learns common and condition-specific gene modules across multiple samples, so biological differences between conditions are not removed as technical effects.
- **Strengths:** Mature software, scalable implementations (RcppPlanc), and an explicit shared/specific split.
- **Limitations:** K and regularisation (lambda) must be tuned. Designed for integrating cell types or conditions rather than malignant ITH with hundreds of patients. Sample-specific factors can soak up shared signal. No CNV or QC modelling.
- **Relevance:** These are joint rather than per-sample models: a 'shared versus patient-specific program' decomposition is exactly what MP pipelines approximate with overlap thresholds. Natural baselines and building blocks.

### scDEF (deep exponential family, hierarchical shared and patient-specific states)

*Beerenwinkel lab (ETH Zurich), cbg-ethz/scDEF. bioRxiv 10.1101/2022.10.15.512383 (v2, 'Identifying hierarchical cell states and gene signatures with deep exponential families for single-cell transcriptomics'). I did not check the author list or whether a peer-reviewed version exists.* — ✅ verified online — <https://github.com/cbg-ethz/scDEF>

- **Core idea:** Multi-level Bayesian (Gamma) matrix factorisation that learns a hierarchy of cell states and gene signatures. It also performs batch integration and can be guided by known gene sets. In 9 high-grade serous ovarian cancer patients it found both patient-specific and shared states.
- **Strengths:** Count-based Bayesian model, hierarchical interpretability, and maintained Python package (updated 2025).
- **Limitations:** Small cohorts demonstrated. Number of factors and layers still to be chosen or pruned. CNV is not modelled. Peer-reviewed status unclear.
- **Relevance:** A hierarchical model (shared and specific programs, coarse and fine) is a principled replacement for flat MP clustering and could represent 'MP families'.

### Consensus scHPF (multi-patient hierarchical Poisson factorization)

*Levitin HM, Zhao W, Bruce JN, Canoll P, Sims PA. bioRxiv 10.1101/2023.12.05.570193 (Dec 2023)* — ✅ verified online — <https://www.biorxiv.org/content/10.1101/2023.12.05.570193v1>

- **Core idea:** Runs scHPF (count-based Bayesian factorisation with no normalisation) from many random initialisations, keeps the most robust factors and determines the number of factors automatically. Applied to 52 drug-treated glioma slice-culture samples from 19 patients across 3 cell types.
- **Strengths:** Count likelihood and automatic K.
- **Limitations:** Preprint. Focused on drug response rather than pan-cancer ITH. No explicit handling of patient or CNV confounders.
- **Relevance:** Combines automatic K selection with consensus over seeds on pooled multi-patient data, which addresses two of the user's complaints (seed and rank dependence).

### CSI-GEP (GPU, consistent, atlas-scale GEP recovery)

*Liu X, Chapple RH, Bennett D, Wright WC, Sanjali A, Culp E, Zhang Y, Pan M, Geeleher P. Cell Genomics (online Jan 2025; DOI 10.1016/j.xgen.2024.100739)* — ✅ verified online — <https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00368-9>

- **Core idea:** GPU-based unsupervised NMF framework for atlas-scale scRNA-seq. It targets three weaknesses of existing GEP methods: poor scalability, inconsistency across repeated runs, and no established way to choose key parameters such as K. It reports better recovery of ground-truth GEPs than competing methods.
- **Strengths:** Scales to atlas size and evaluates consistency explicitly.
- **Limitations:** Works on a single atlas or matrix, with no cross-tumour aggregation and no CNV or technical modelling. I did not read the details of how it selects K.
- **Relevance:** Shows that consistency across runs and principled K selection at scale are recognised, publishable problems. GPU acceleration is needed if a new MP method runs consensus or bootstrap analysis over thousands of samples.

### sciRED (confounder-aware interpretable factor decomposition)

*Pouyabahar D, Andrews T, Bader GD. Nature Communications 16 (2025)* — ✅ verified online — <https://www.nature.com/articles/s41467-025-57157-2>

- **Core idea:** Regresses out known confounders (e.g., sample, library size) before factor analysis, rotates factors (varimax) for interpretability, matches factors to covariates with an ensemble of classifiers, and flags 'unexplained' factors as candidate new biology.
- **Strengths:** Interpretability and systematic matching of factors to covariates.
- **Limitations:** Linear (PCA-based) factors. Confounders must be known. It does not address recurrence across tumours.
- **Relevance:** A template for an explicit step that labels each program as technical, covariate-driven (e.g., CNV, study, chemistry, dissociation) or unexplained. MP pipelines currently do this by hand.

### DeltaNMF (baseline versus case-specific programs)

*Karpurapu A, Gersbach C, Singh R. bioRxiv 10.64898/2026.01.22.701049 (Jan 2026); listed in RECOMB 2026 materials. github.com/rohitsinghlab/deltanmf* — ✅ verified online — <https://www.biorxiv.org/content/10.64898/2026.01.22.701049v1>

- **Core idea:** Two-stage neural NMF. It first learns baseline programs from control cells, then learns case-specific programs, separating changed usage of existing programs from genuinely new programs. GPU implementation is about 20x faster than cNMF with comparable accuracy. Applied to CAD Perturb-seq data.
- **Strengths:** Clean conceptual decomposition and speed.
- **Limitations:** Designed for case-control and perturbation data, not multi-tumour recurrence. Preprint.
- **Relevance:** Gavish found that many carcinoma MPs mirror programs in normal epithelium. A DeltaNMF-like normal-baseline stage could separate cancer-specific programs from re-used normal-lineage programs, which no MP pipeline does formally.

### Human Cell Census NMF reference embedding

*Liu Y-T, Triche TJ Jr, DeBruine ZJ. bioRxiv 10.64898/2025.12.30.697085 (Dec 2025)* — ✅ verified online — <https://www.biorxiv.org/content/10.64898/2025.12.30.697085v1>

- **Core idea:** Fixed NMF reference fitted to 28.5M healthy cells and about 60k genes from the Human Cell Census. New datasets can be projected onto its additive programs without retraining.
- **Strengths:** Atlas scale and reusable projection.
- **Limitations:** Healthy tissue only. A single fixed factorisation without uncertainty. Preprint.
- **Relevance:** Could serve as a normal-program coordinate system: malignant programs could be expressed as 'normal program plus deviation', and lineage or contamination programs flagged automatically.

### DIALOGUE (multicellular programs), brief

*Jerby-Arnon L, Regev A. Nature Biotechnology (2022). github.com/livnatje/DIALOGUE* — ✅ verified online — <https://www.nature.com/articles/s41587-022-01288-0>

- **Core idea:** Penalised matrix decomposition (multi-CCA-like) on sample- or niche-level representations of several cell types, combined with multilevel modelling. It finds multicellular programs that co-vary across samples or niches and links them to phenotypes (Alzheimer's disease, ulcerative colitis, immunotherapy resistance).
- **Strengths:** Phenotype association and cross-cohort replication.
- **Limitations:** Captures variation between samples or niches across cell types, not programs within malignant cells. Number of programs and covariates are user-set.
- **Relevance:** Peripheral. It could add a layer linking MP activity to coordinated TME programs, and its multilevel modelling of cell-quality covariates is a useful idea.

### scITD (tensor decomposition of multicellular inter-individual variation), brief

*Mitchel J et al. Nature Biotechnology (2024)* — ✅ verified online — <https://www.nature.com/articles/s41587-024-02411-z>

- **Core idea:** Builds a donor x gene x cell-type pseudobulk tensor and applies Tucker decomposition to find coordinated axes of inter-individual variation (lupus, COVID-19).
- **Strengths:** Stratifies cohorts and is interpretable.
- **Limitations:** Works at donor (pseudobulk) level, so it does not resolve intratumour programs.
- **Relevance:** Peripheral. It models inter-patient structure that MP pipelines remove by centring. Useful if the new method also outputs patient-level 'baseline' programs.

### MOFAcell and Tensor-cell2cell (sample-level multicellular factors), brief

*Ramirez Flores RO et al. (Saez-Rodriguez lab), eLife (2023) 'Multicellular factor analysis of single-cell data for a tissue-centric understanding of disease'; Armingol E et al., Nature Communications (2022) 'Context-aware deconvolution of cell-cell communication with Tensor-cell2cell'* — ✅ verified online — <https://elifesciences.org/articles/93161>

- **Core idea:** MOFAcell runs MOFA on pseudobulk profiles per cell type, treated as views, to get multicellular factors across samples. Tensor-cell2cell decomposes a context x ligand-receptor x sender x receiver tensor to find communication modules.
- **Strengths:** Mature tooling (MOFAcellulaR, LIANA) and flexible designs across conditions.
- **Limitations:** Sample-level or communication-level, not intratumour malignant programs.
- **Relevance:** Peripheral. Possible downstream modules for interpreting MPs in their tissue context.

**Gaps identified by this sweep**

- Every major pipeline ends with manual curation, and nothing does it automatically. Gavish 2023 removed MPs flagged as low quality, single-study or doublet-like. 3CA v2 code removes 11 MPs by hand ('LQ' x7, Macro, T-cells, Endo, miss-annotation). Kinker removed complexity-correlated and cell-cycle programs by hand and defined metaclusters by hand. Barkley picked the final modules by name. GeneNMF ships dropMetaPrograms(). No method automatically flags MPs as low-quality (ribosomal/mitochondrial, low complexity), contamination or doublet, ambient RNA, dissociation stress, or chemistry-specific (sc vs sn).
- Seed robustness is never quantified. Gavish keeps the best of 10 restarts (no consensus), Barkley uses a single ICA-initialised run, and GeneNMF uses one fixed seed. None reports how much MPs or their genes vary across seeds. cNMF-style consensus exists within one sample but has not been carried through to MP-level uncertainty.
- Rank and number of MPs are not chosen from the data. K is fixed (4-9, 6-9) or picked by heuristic (Barkley's rank-minus-modules rule). nMP is user-set in GeneNMF, and Kmax in GBCD. The final number of MPs comes from arbitrary overlap thresholds. No method gives a calibrated, data-adaptive number of recurrent programs with confidence.
- Hard, arbitrary thresholds sit on binarised top-50 gene lists: 35/50 within tumour, 10/50 across tumours, 10 genes for clustering, 5 programs minimum, a 7,000-gene universe. Gene loadings are thrown away, gene-list length is fixed, and gene membership has no statistical significance (no FDR). Sensitivity to these parameters has not been systematically reported.
- Cross-tumour aggregation is greedy and order-dependent (Gavish), forces every program into a cluster (Ward plus cutree in GeneNMF), or uses ad hoc graph thresholds (Barkley, MetaTiME, mosaicMPI). There is no generative or probabilistic model of programs being shared across tumours, cancer types and studies, and no hierarchy of MP families versus cancer-type-specific variants.
- Genetic (CNV) confounding is only examined after the fact. Kinker found 39% of discrete clusters and 8% of continuous programs associated with CNA subclones, and Gavish found only 24% of subclones associated with a single MP. No pipeline removes or labels clone-driven programs, residualises expression on inferred clone or CNV state, or tests MPs for positional (chromosome-arm) gene enrichment, which matters for recurrent arm-level CNAs.
- Variation within and between tumours is confounded. Per-sample centring plus zeroing negatives (Gavish), or treating sample as the scVI batch (ITHmapper), removes programs active tumour-wide in some patients. GBCD tackles this but needs Kmax, does not scale to atlas size and has been run on one cancer type at a time. The Tirosh lab's 2025 GBM paper separately describes 'baseline programs' between tumours. No unified atlas-scale model jointly estimates within-tumour programs, patient-level baselines, study and chemistry effects, and CNV effects.
- Study, platform and protocol effects are handled only by heuristics (e.g., 'at least 2 studies'). Chemistry (10x 3' vs 5', Smart-seq2, single-nucleus), dissociation protocol and sample handling are not modelled as covariates, even though the 3CA v2 code compares sc and sn versions of the same MPs.
- Validation is weak and there is no ground truth. The only cancer-specific benchmark (CanSig, Cancer Research 2026) found late-integration signatures (GeneNMF, scalop) less rediscoverable than the best early-integration methods. There is no realistic simulator that combines patient baselines, subclonal CNVs, doublets and ambient RNA, dissociation stress and study effects with known recurrent programs. There is also no standard held-out-cohort replication metric, and spatial, perturbation or lineage evidence is rarely built into method evaluation.
- There is no uncertainty-aware reference catalog of malignant programs with a projection model. TCAT/starCAT covers T cells, ITHmapper maps to discrete states per cancer type, and 3CA provides fixed gene lists scored with ad hoc module scores. No versioned malignant-program reference exists with per-cell usages that are comparable across samples, QC flags and per-gene confidence.
- Normal versus malignant programs are not separated formally. Gavish showed that many carcinoma MPs resemble normal-epithelium programs, but no MP method uses a normal-tissue baseline (DeltaNMF-like, or a Human Cell Census NMF reference) to separate cancer-specific programs from reused lineage and physiological programs.
- Scale and compute: per-sample NMF across ranks and restarts over about 2,800+ samples, especially with bootstrap or consensus layers, needs GPU or streaming implementations (CSI-GEP, DeltaNMF, RcppML show this is possible). No MP tool provides this with stability analysis.
- Continuous versus discrete: tools mix continuous programs (NMF) with discrete states (ITHmapper, CanSig clustering), and none models overlapping program usage while also reporting state-level structure and programs that switch together.
- Clinical and biological grounding is post hoc (TCGA correlations, survival). Few methods treat external replication, spatial co-localisation or functional relevance (e.g., dependency or perturbation data) as part of program selection or ranking.

**Notes:** Access constraints in this session:
- nature.com, biorxiv.org, PMC/PubMed/EuropePMC, CRAN, rdrr.io, researchgate, semanticscholar, weizmann.ac.il, zenodo and aacrjournals were blocked by the egress proxy.
- Exact procedures were therefore checked from the authors' public code on GitHub where possible:
  - tiroshlab/3ca: Generate_Meta_Programs.R, robust_nmf_programs.R, 3CA_v2/MP_analysis_new_batch.R
  - gabrielakinker/CCLE_heterogeneity (module2_rhp.R)
  - yanailab/PanCancer (part2.R, Finding.R, seurat_functions_public.R)
  - carmonalab/GeneNMF (R/main.R, DESCRIPTION, cancer demo Rmd)
  - abelson-lab/ITHmapper-manuscript (scVI and Hotspot scripts)
- Everything else comes from search-result snippets. The session's web-search budget (200 calls) ran out near the end, so a final search for papers that explicitly test MP seed or initialisation instability was not run.

Not verified:
- The minimum malignant-cell count per tumour used by 3CA.
- The exact count of malignant MPs in 3CA v2 (149 is the total across 8 cell types).
- Full author lists for MetaTiME, scDEF and scITD.
- Whether scDEF has a peer-reviewed version.
- The JID 2025 publication of the GeneNMF-citing BCC paper.

Corrections to the task framing:
- Barkley 2022 did not use cNMF. It used nsNMF with seed='ica', nrun=1, ranks 2-25, and a rank-minus-modules rule, then an infomap gene co-occurrence graph and manual module selection.
- I found no peer-reviewed paper (2023-2026) whose main purpose is to critique the Gavish MP procedure.
- The strongest published evidence against it is indirect:
  - CanSig (Cancer Research 2026): MP-style late-integration signatures (GeneNMF, scalop) were less rediscoverable across datasets.
  - GBCD (Nature Genetics 2025): tumour-by-tumour NMF misses shared subtype programs.
  - The Tirosh lab's own 3CA v2 code needed 11 MPs removed by hand.
- Several small independent reimplementations exist, a sign of demand for a proper package: pyscalop (Python port of scalop, May 2026) and mpnmf (Python, Apr 2026, deterministic NNDSVDa initialisation, ranks 7-12, thresholds 0.7/0.2/0.3). Neither is published.

Positioning advice for a Nature Methods-level package:
1. The precedent is TCAT/starCAT (Nature Methods 2025): cross-dataset reproducibility as the core criterion, a fixed versioned catalog and a projection tool. A malignant-cell analogue alone may look incremental unless it solves problems TCAT does not:
   - automatic artefact and CNV labelling
   - joint modelling of within-tumour, patient-baseline and study effects
   - calibrated uncertainty for MP existence and gene membership
   - data-driven number of MPs
2. The closest direct competitors to benchmark against:
   - the Gavish/3CA code itself (not only scalop)
   - GeneNMF
   - multi-dataset cNMF with Harmony-nonnegative correction (TCAT-style)
   - GBCD
   - CanSig's top early-integration methods (Harmony, BBKNN, fastMNN plus clustering and DE)
   - ITHmapper (scVI plus Hotspot)
   - mosaicMPI
   - LIGER/iNMF or scINSIGHT
3. Evaluate on:
   - the CanSig datasets and metrics (an established harness)
   - the 3CA v2 corpus, using its 11 hand-removed MPs as a labelled artefact set
   - a new realistic simulator with known recurrent programs plus patient baselines, subclonal and recurrent arm-level CNVs, doublets and ambient RNA, dissociation-stress modules, chemistry and study effects
   - external validation: held-out cohorts, spatial co-localisation, Kinker and Tahoe-100M cell lines, perturbation data
4. Design ideas from the literature:
   - consensus or bootstrap over seeds, ranks and samples (cNMF, consensus scHPF, mosaicMPI multi-resolution)
   - a probabilistic model of sharing and hierarchy (GBCD, scDEF)
   - nuisance factors anchored by priors (Spectra) and covariate matching (sciRED)
   - automatic false-positive state filtering (EcoTyper AFI)
   - a normal-baseline contrast (DeltaNMF; Human Cell Census NMF reference)
   - clone- and positional-aware residualisation using inferred CNV (infercnv, CopyKAT, Numbat)
   - GPU scaling (CSI-GEP, RcppML)
5. The Tirosh lab's 2025 GBM paper (Nature Genetics 57:1155-1167) separates malignant states from three between-tumour 'baseline programs'. It shows the field is already moving toward multilayer decomposition, which a new method should formalise.


## Gene-program / latent-factor methods (25 items)

### Current standard: Tirosh-lab / 3CA metaprogram pipeline (Gavish 2023; Tyler 2025), plus GeneNMF and SSpMosaic implementations

*Gavish A, Tyler M, ..., Tirosh I. Nature 2023 ('Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours'); Tyler M, Gavish A, ..., Tirosh I. Nature Cancer 6:1088-1101 (2025) ('The Curated Cancer Cell Atlas provides a comprehensive characterization of tumors at single-cell resolution'); GeneNMF R package (Andreatta M, Carmona SJ; CRAN, cited through Yerly et al. bioRxiv 2024); SSpMosaic (compbioNJU), Cell Genomics 2025* — ✅ verified online — <https://www.nature.com/articles/s41586-023-06130-4>

- **Core idea:** Run NMF separately for each tumor over a range of ranks (K=4..9 in Gavish). Keep programs that recur across ranks within a tumor (top-50 gene overlap) and that overlap at least 20% with a program from another tumor. Then cluster programs across tumors into metaprograms (MPs): 41 MPs from 1,163 tumours in 2023, updated with 124 datasets and 2,836 samples in 2025. GeneNMF v0.6 replaced top-gene Jaccard with cosine similarity on full gene-weight vectors, averages weights into consensus MP weights, and adds a 'specificity.weight' to make loadings sparser. SSpMosaic builds differential co-expression programs inside each dataset and links them with network propagation and hierarchical clustering into 'metaprograms', which it uses as integration and annotation anchors. Likelihood: none (Frobenius NMF on centered or relative expression with negatives clipped). K: fixed range plus heuristics. Multi-sample: per-tumor analysis, then clustering. Interpretability: binary gene lists. Scale: embarrassingly parallel. Tumor MPs: yes, the defining use case.
- **Strengths:** Simple, transparent, and widely adopted. It has a huge curated resource (3CA) and is biologically validated in many follow-up studies.
- **Limitations:** No generative or statistical model and no uncertainty. The result rests on many hard thresholds (K range, top-50 genes, overlap cutoffs, minimum cluster size, greedy clustering order) and typically one seed per K. MPs and their gene lists change with cohort composition. Patient-specific and CNV-driven programs are removed only through the recurrence heuristic. The method cannot borrow strength across tumors when a tumor has few cells. Gene membership is all-or-nothing.
- **Relevance:** This is the baseline a new package must beat. A strong paper would show that the new method recovers the real Gavish/3CA MPs and is more stable across seeds, ranks, and cohort subsampling. It should also flag or split MPs that are artifacts or patient/CNV-specific.

### cNMF (consensus NMF)

*Kotliar D, Veres A, Nagy MA, Tabrizi S, Hodis E, Melton DA, Sabeti PC. eLife 8:e43803 (2019)* — ✅ verified online — <https://elifesciences.org/articles/43803>

- **Core idea:** Run NMF many times (typically 100+ replicates) on variance-scaled highly variable gene (HVG) counts. Pool all components, remove outliers by local density, cluster with k-means, and take cluster medians as consensus gene expression programs (GEPs). Then refit cell usages. Likelihood: Frobenius by default (beta-divergence options), not a count model. K: chosen by eye from a stability (silhouette) versus reconstruction-error plot. Multi-sample: not modeled. Batch effects produce duplicated, dataset-specific GEPs unless the data are corrected first. Interpretability: non-negative identity and activity programs with gene spectra scores. Scale: CPU-heavy because of replicates; GPU ports exist. Tumor MPs: commonly applied per tumor in cancer studies.
- **Strengths:** Robust to random initialization, well benchmarked, widely used, and separates identity programs from activity programs.
- **Limitations:** Single-dataset model. K is picked by visual inspection. The outlier-density threshold is arbitrary. There is no per-gene or per-program uncertainty and no handling of patient or CNV structure. Consensus removes seed dependence but not rank or threshold dependence at the cross-tumor MP step.
- **Relevance:** Introduced the consensus-over-seeds idea that most later tools reuse. A new MP method should include it as a baseline and should solve the cross-tumor consensus step that cNMF leaves to ad hoc clustering.

### starCAT / TCAT (multi-dataset consensus GEP catalog plus projection)

*Kotliar D, et al. (Raychaudhuri lab). Nature Methods 2025 (s41592-025-02793-1); bioRxiv 2024 (10.1101/2024.05.03.592310)* — ✅ verified online — <https://www.nature.com/articles/s41592-025-02793-1>

- **Core idea:** Run cNMF on 7 datasets (1.7M T cells, 38 tissues, about 700 individuals). The data are first passed through a Harmony variant adapted to output batch-corrected, non-negative gene-level values, because standard correction is incompatible with NMF. GEPs are then clustered across datasets into 46 reproducible consensus GEPs. starCAT scores new datasets against the fixed catalog using non-negative least squares on the fixed GEPs. Likelihood: Frobenius NMF. K: per-dataset cNMF selection. Multi-sample: batch correction first, then cross-dataset clustering. Interpretability: high, with annotated program catalog. Scale: millions of cells. Tumor MPs: no (T cells; there are also glioma myeloid and bone marrow references).
- **Strengths:** Shows the evaluation bar Nature Methods accepted: reproducibility across datasets, annotation accuracy, and reusable references and software.
- **Limitations:** Cross-dataset consensus still relies on correlation or clustering thresholds. Harmony-corrected values are not a count model. The catalog is fixed, so new programs require rerunning. Shared and dataset-specific programs are not separated explicitly. It does not address malignant ITH or CNV confounding.
- **Relevance:** The closest Nature Methods precedent for the paper the user has in mind: a reproducible multi-dataset program catalog plus a projection tool, benchmarked on cross-dataset reproducibility. A principled, model-based malignant-MP catalog with projection (an 'MP-CAT') would be a natural and clearly differentiated successor.

### scHPF and consensus scHPF (hierarchical Poisson factorization)

*Levitin HM, et al. (Sims lab). Molecular Systems Biology 15:e8557 (2019); consensus scHPF: Levitin HM, Zhao W, Bruce JN, Canoll P, Sims PA. bioRxiv 2023 (10.1101/2023.12.05.570193)* — ✅ verified online — <https://link.springer.com/article/10.15252/msb.20188557>

- **Core idea:** A Bayesian Poisson factorization with gamma priors and cell- and gene-specific 'capacity' hyperpriors that model variable sparsity. It runs directly on raw UMI counts without normalization, using coordinate-ascent variational inference. Consensus scHPF fits many models across random initializations and K values, clusters factors by gene scores into recurrent modules, and initializes a final model from module medians. The number of modules becomes K, so K is chosen automatically. It was applied to drug-treated glioma slices from 19 patients and 52 samples. Likelihood: Poisson (count-native). K: heuristic (original) or consensus-automated. Multi-sample: pooled, with no explicit patient terms. Interpretability: sparse non-negative gene scores. Scale: moderate (CPU, mean-field VI). Tumor: yes (glioma).
- **Strengths:** Statistically matched to UMI data and interpretable. Handles both discrete and continuous programs, and held-out likelihood can be used for model comparison.
- **Limitations:** Patient, batch, and CNV structure is not modeled. Consensus K still depends on clustering heuristics. Mean-field VI underestimates uncertainty. The consensus version is still an unpublished preprint.
- **Relevance:** Shows that a count likelihood plus consensus-based automatic K is feasible for multi-patient tumor data. It is the nearest NMF-family precedent to a count-native MP method.

### fastTopics / Poisson NMF with GoM DE, LDA topic models (CountClust), and Topyfic (reproducible LDA)

*Carbonetto P, Sarkar A, Wang Z, Stephens M. Statistics and Computing 36:131 (2026) (arXiv 2105.13440); Carbonetto P, Luo K, Sarkar A, Hung A, Tayeb K, Pott S, Stephens M. Genome Biology 24:236 (2023) (GoM DE); Dey KK, Hsiao CJ, Stephens M. PLoS Genetics 13:e1006599 (2017) (CountClust); Topyfic (Mortazavi lab), bioRxiv 2024 (10.1101/2024.02.26.582178)* — ✅ verified online — <https://link.springer.com/article/10.1007/s11222-026-10866-0>

- **Core idea:** fastTopics uses the equivalence between the multinomial topic model and Poisson NMF to fit topic models quickly with NMF algorithms, which often gives much better fits. GoM DE then estimates per-topic differential expression (log fold-changes with uncertainty) while allowing grades of membership. CountClust brought grade-of-membership and LDA models to RNA-seq. Topyfic runs LDA over many seeds, clusters topics with Leiden (using Harmony across groups), and drops low-participation topics to get reproducible topics. Likelihood: Poisson/multinomial. K: user-chosen. Multi-sample: none (Topyfic: Harmony at the topic-aggregation step). Interpretability: high, with statistically grounded distinctive genes. Scale: good (fastTopics), moderate (LDA). Tumor MPs: not targeted.
- **Strengths:** Count-native, fast, and statistically rigorous about which genes define a program.
- **Limitations:** K must be chosen by the user. There are no batch or patient terms. The LDA posterior is multimodal, and multi-seed aggregation is heuristic. Topics mix lineage and state signals when patient effects dominate.
- **Relevance:** GoM DE is a principled replacement for 'top-50 genes' when deciding MP gene membership. Topyfic shows that multi-seed aggregation reduces topic irreproducibility, which is the same issue as MP seed dependence.

### LIGER iNMF / UINMF / consensus iNMF (and scINSIGHT)

*Welch JD, Kozareva V, Ferreira A, Vanderburg C, Martin C, Macosko EZ. Cell 2019 (LIGER); Kriebel AR, Welch JD. Nature Communications 13:780 (2022) (UINMF); rliger runCINMF consensus iNMF (released Oct 2024); Qian K, Fu S, Li H, Li WV. Genome Biology 23 (2022) (scINSIGHT)* — ✅ verified online — <https://www.nature.com/articles/s41467-022-28431-4>

- **Core idea:** Integrative NMF models each dataset as X_i ≈ (W + V_i)H_i, where W holds metagenes shared across datasets and V_i holds dataset-specific metagenes, controlled by a penalty λ. Online iNMF scales to large data, and UINMF adds features not shared across datasets. Consensus iNMF (rliger 2.x, 2024) runs iNMF from several random starts (nRandomStarts, default 10), clusters the pooled W/V factors, and takes the consensus of the largest clusters, following cNMF. scINSIGHT models common and condition-specific gene modules. Likelihood: Frobenius on scaled, non-centered data. K and λ: user-set. Multi-sample: explicit shared versus dataset-specific structure. Interpretability: metagenes. Scale: millions of cells (online). Tumor MPs: not targeted.
- **Strengths:** Mature software and very scalable. The shared/specific decomposition is explicit, and a consensus variant now exists.
- **Limitations:** No count likelihood, and K and λ are chosen by hand. V_i is not sparse or interpretable by design and can absorb biology and batch effects alike. Built for cell-type integration rather than within-tumor heterogeneity. With about 1,000 tumors as 'datasets', the number of parameters grows fast and identifiability between W and V_i is weak.
- **Relevance:** Its shared-plus-specific decomposition is exactly the structural prior MPs need: shared MPs versus patient-specific programs. A count-likelihood, sparse, hierarchical version of iNMF with automatic K is an obvious and still unbuilt design.

### Spectra (supervised, prior-informed gene program discovery)

*Kunes RZ, Walle T, ..., Pe'er D. Nature Biotechnology (online Sept 2023), s41587-023-01940-3* — ✅ verified online — <https://www.nature.com/articles/s41587-023-01940-3>

- **Core idea:** Spectra factorizes expression into global and cell-type-specific factors. Input gene sets are encoded as a gene-gene knowledge graph and act as a graph penalty (stochastic-block-model-style prior with eta/kappa/rho parameters). Gene sets can be edited by the data, and new factors can be learned. Each factor gets an η score measuring how much it depends on the prior. The default is one factor per gene set plus one per cell type, and Cytopus provides 231 mostly immune gene sets. Likelihood: works on normalized expression; to my understanding not a full count model (not re-verified in this session). K: set by the gene sets plus user additions. Multi-sample: no explicit batch or patient model. Interpretability: very high. Scale: hundreds of thousands of cells, faster than scHPF and Slalom. Tumor: tumor-immune contexts (CD8 T cell exhaustion and tumor reactivity, macrophages); not malignant MPs.
- **Strengths:** High interpretability. Combines known and new programs and quantifies how much each program relies on prior knowledge.
- **Limitations:** Needs cell-type labels and curated priors, which barely exist for malignant ITH states. Priors can pull factors toward known biology and hide new malignant programs. Patient and CNV effects are not modeled.
- **Relevance:** Shows that soft knowledge priors improve coherence and interpretability (benchmarks show unconstrained NMF and scHPF factors match annotated gene sets poorly). For MPs, a soft prior from pathway graphs or foundation-model gene embeddings, with a per-program 'novelty' score like η, is attractive.

### expiMap and pathway-informed VAEs (VEGA, pmVAE, OntoVAE)

*Lotfollahi M, Rybakov S, Hrovatin K, Hediyeh-zadeh S, Talavera-López C, Misharin AV, Theis FJ. Nature Cell Biology 25:337-350 (2023); VEGA: Seninge L, et al. Nature Communications 2021; pmVAE: Gut G, et al. ICML CompBio 2021 / bioRxiv 2021* — ✅ verified online — <https://www.nature.com/articles/s41556-022-01072-x>

- **Core idea:** expiMap is a conditional VAE built on scArches with a masked linear decoder. Each latent node corresponds to one database gene program. Group lasso switches off irrelevant programs, an L1 penalty lets programs gain genes, and a few unconstrained nodes learn new programs. It supports reference mapping of query data. VEGA, pmVAE, and OntoVAE likewise wire the decoder to pathway or ontology memberships. Likelihood: count-type reconstruction available in the scArches family (exact default not re-verified). K: number of input gene sets plus a few new nodes. Multi-sample: conditions or batches as CVAE covariates. Interpretability: tied to annotations. Scale: GPU, atlas-scale. Tumor MPs: no.
- **Strengths:** Supports reference mapping, models covariates, and links programs to known pathways.
- **Limitations:** Interpretability is limited by the gene-set databases, and only a few new programs can be learned. Latent nodes for overlapping gene sets are correlated. Not designed for ITH programs or tumor-specific confounders.
- **Relevance:** Provides a template for combining priors, new programs, and query mapping across cohorts. However, malignant MPs are largely data-driven, so priors must stay soft.

### MuVI (structured-sparsity Bayesian factor model) and MOFA+ (group factor analysis with ARD)

*Qoku A, Buettner F. AISTATS 2023 (PMLR 206); Argelaguet R, Arnol D, Bredikhin D, et al. Genome Biology 21:111 (2020) (MOFA+)* — ✅ verified online — <https://proceedings.mlr.press/v206/qoku23a/qoku23a.pdf>

- **Core idea:** MuVI is a multi-view Bayesian factor model with a modified regularized horseshoe prior. Noisy prior gene sets encode structured sparsity on loadings, the model can correct false-positive and false-negative annotations, and dense unconstrained factors are added. MOFA+ is Bayesian group factor analysis with automatic relevance determination (ARD) across views and sample groups, fitted by GPU stochastic VI. Likelihood: Gaussian for MuVI; Gaussian, Bernoulli, or Poisson for MOFA+. K: user-set upper bound with sparsity or ARD shrinkage. Multi-sample: MOFA+ groups (factor activity varies by group, loadings are shared). Interpretability: sparse loadings tied to priors. Scale: moderate (MuVI) to large (MOFA+). Tumor MPs: no.
- **Strengths:** Principled sparsity, robust to noisy priors, and gives automatic relevance or factor pruning.
- **Limitations:** Gaussian likelihood on normalized data. MOFA+ groups share loadings, so patient-specific programs cannot be represented. MuVI targets multi-omics views rather than many patients. Neither models CNV.
- **Relevance:** Horseshoe or ARD priors are the right tools for automatic K and sparse, calibrated gene membership. Combined with count likelihoods and patient-level hierarchies, they would remove the K-range and top-N thresholds.

### Hotspot (local-autocorrelation gene modules)

*DeTomaso D, Yosef N. Cell Systems 12:446-456 (2021)* — ✅ verified online — <https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00114-9>

- **Core idea:** Computes each gene's local autocorrelation on a cell kNN graph built from any latent space, using depth-adjusted null models, to find informative genes. Pairwise local correlations are then clustered hierarchically into gene modules and scored per cell. Likelihood: null-model-based Z-scores, not a generative model. K: emerges from clustering under minimum-gene and FDR thresholds. Multi-sample: inherited from the upstream latent space (for example scVI). Interpretability: disjoint gene modules. Scale: O(G²) gene-pair computation limits the number of genes. Tumor MPs: not designed for them; usable per tumor.
- **Strengths:** No K to choose, sensitive for sparse genes, and modality-agnostic.
- **Limitations:** Modules are hard and disjoint, so a gene cannot belong to two programs. Depends on the kNN graph, which for pooled malignant cells is dominated by patient identity. No mixed membership and no explicit sharing across tumors.
- **Relevance:** An alternative per-tumor program detector that avoids rank selection. Its graph-based autocorrelation statistic could serve as an independent check of whether MP genes really co-vary within each tumor.

### scETM embedded topic model and successors scE2TM and scConcept

*Zhao Y, Cai H, Zhang Z, Tang J, Li Y. Nature Communications 12:5261 (2021) (scETM); scE2TM: Nature Communications 2026 (s41467-026-76825-5; arXiv 2507.08355); scConcept: Chen H, Li Y. bioRxiv 2026 (10.64898/2026.04.21.719959)* — ✅ verified online — <https://www.nature.com/articles/s41467-021-25534-2>

- **Core idea:** scETM pairs a VAE encoder with a linear tri-factorized decoder: softmax(θ·α·ρᵀ + batch intercepts). Topics and genes share an embedding space, the likelihood is multinomial, and batches enter as linear intercepts. It scales beyond 10⁶ cells, supports zero-shot transfer, and can take pathway-informed gene embeddings. scE2TM adds embedding-clustering regularization so each topic is the center of its own gene cluster, which fixes 'interpretation collapse' (redundant topics). It also brings in single-cell foundation model knowledge through a cross-view encoder. scConcept uses an LLM to turn neural-topic gene programs into named 'concepts' (label, description, gene set) and maps them back onto cells. K: user-set. Interpretability: high. Tumor MPs: no.
- **Strengths:** Scalable, transferable, and the gene embedding space makes priors easy to inject.
- **Limitations:** Batch is modeled only as an additive intercept, with no patient-level program structure. Topic redundancy needs explicit regularization. LLM-generated names risk hallucination. None target malignant ITH.
- **Relevance:** A count-likelihood factor model whose genes live in an embedding space is the natural place to plug in foundation-model or text gene embeddings as soft priors. The anti-redundancy regularizer speaks directly to non-redundant MPs.

### LDVAE (linearly decoded scVI)

*Svensson V, Gayoso A, Yosef N, Pachter L. Bioinformatics 36(11):3418-3421 (2020)* — ✅ verified online — <https://academic.oup.com/bioinformatics/article/36/11/3418/5807606>

- **Core idea:** scVI with a linear decoder, so decoder weights form a gene-loading matrix. Uses a negative binomial likelihood, an amortized nonlinear encoder, batch as a covariate, and scales to more than 1M cells. Likelihood: NB. K: user-set. Multi-sample: batch covariate. Interpretability: loadings, but signed and dense. Scale: excellent. Tumor MPs: no.
- **Strengths:** Count-native, fast, and integrated into scvi-tools.
- **Limitations:** Loadings are neither non-negative nor sparse, factors are rotation-ambiguous and correlated, and the model has less capacity for complex batch effects. No automatic K and no shared versus specific structure.
- **Relevance:** Shows that count-likelihood linear factor models scale with amortized variational inference. An MP model needs this scalability plus non-negativity, sparsity, and a tumor hierarchy.

### DRVI (disentangled representation variational inference)

*Moinfar AA, Theis FJ. bioRxiv 2024, latest version 2026 (10.1101/2024.11.06.622266); ICLR 2025 LMRL workshop; available in scvi-tools. No journal version found as of this search.* — ✅ verified online — <https://www.biorxiv.org/content/10.1101/2024.11.06.622266v3>

- **Core idea:** An NB VAE with scVI-compatible batch covariates whose decoder is a sum of per-latent-dimension subnetworks pooled with log-sum-exp. This makes each latent dimension a disentangled, possibly nonlinear 'gene program'. Built-in analyses map dimensions back to genes. n_latent sets an upper bound (default 32), and inactive dimensions are flagged as 'vanished', which acts as a soft automatic K. Integration quality matches entangled scVI-type methods, and it decomposed the Human Lung Cell Atlas into interpretable factors. Likelihood: NB. Multi-sample: batch covariates. Scale: atlas-scale on GPU. Tumor MPs: not specifically targeted.
- **Strengths:** Integration and interpretability together, count-native, and software is available.
- **Limitations:** Programs are nonlinear, so gene contributions come from decoder interrogation rather than a single sparse loading vector. Identifiability is not guaranteed in general. There is no explicit patient-specific versus shared program structure and no CNV handling.
- **Relevance:** The strongest deep-learning competitor: it integrates data and yields interpretable programs in one model. Any MP package must benchmark against DRVI on reproducibility and biological validity.

### GSFA (guided sparse factor analysis)

*Zhou Y, Luo K, Liang L, Chen M, He X. Nature Methods 20:1693-1703 (2023)* — ✅ verified online — <https://www.nature.com/articles/s41592-023-02017-4>

- **Core idea:** Bayesian sparse factor analysis with spike-and-slab priors on gene loadings. Factors are regressed on perturbation indicators, which 'guides' them. Fitted by Gibbs sampling, it outputs posterior inclusion probabilities (PIPs) for gene-factor and perturbation-factor effects and propagates them to gene-level perturbation effects. Likelihood: Gaussian on transformed or residualized data. K: user-set with sparsity shrinkage. Multi-sample: covariates can be regressed out. Interpretability: calibrated PIPs. Scale: Gibbs sampling limits size (tens of thousands of cells). Tumor MPs: no (CRISPR screens).
- **Strengths:** Full Bayesian uncertainty, interpretable sparse modules, and a high-impact methods precedent.
- **Limitations:** Gaussian on transformed data, slow MCMC, built for perturbation designs, and no multi-tumor hierarchy.
- **Relevance:** Shows Nature Methods values calibrated gene membership (PIPs) over ad hoc top-N lists. 'Guiding' factors by sample-level covariates (CNV, subtype, treatment) is a useful design idea for MPs.

### scDEF (single-cell deep exponential families)

*Ferreira PF, Kuipers J, Beerenwinkel N. bioRxiv 2022, revised v2 ('Identifying hierarchical cell states and gene signatures with deep exponential families for single-cell transcriptomics'); no journal version found* — ✅ verified online — <https://www.biorxiv.org/content/10.1101/2022.10.15.512383v2>

- **Core idea:** A multi-level Gamma-Poisson deep exponential family, essentially hierarchical Bayesian NMF. It enforces non-negativity and sparsity, links fine and coarse factors in a hierarchy, and selects model size automatically through relevance-type pruning. It also supports batch integration, and the informed variant (iscDEF) accepts known gene lists. On 9 high-grade serous ovarian cancer patients it found both patient-specific and shared cell states. Likelihood: Poisson. K: automatic. Multi-sample: batch-aware. Interpretability: multi-resolution signatures. Scale: GPU VI; atlas-scale behavior unclear. Tumor: yes (ovarian cancer).
- **Strengths:** Principled, sparse, hierarchical, and selects model size automatically.
- **Limitations:** Still a preprint. Scaling to thousands of tumors is unproven. Patient-specific versus shared programs emerge implicitly rather than as an explicit decomposition. CNV is not modeled.
- **Relevance:** The existing model closest to the ideal MP engine: count likelihood, automatic K, hierarchy, and shared plus patient-specific states in tumors. It shows the idea is plausible but not yet delivered at pan-cancer scale.

### Bayesian NMF with evidence-based or consensus rank: CoGAPS 3 and ccfindR

*CoGAPS 3: Sherman TD, Gao T, Fertig EJ. BMC Bioinformatics (2020); protocol: Johnson JAI, ..., Fertig EJ, Stein-O'Brien GL. Nature Protocols 18:3690-3731 (2023); ccfindR: Woo J, Aliferis C, Wang J et al., Bioconductor package* — ✅ verified online — <https://www.nature.com/articles/s41596-023-00892-x>

- **Core idea:** CoGAPS is MCMC Bayesian NMF with an atomic sparsity prior. Version 3 adds asynchronous parallel updates, sparse data structures, and distributed runs on random subsets of genes or cells whose patterns are matched afterwards (a consensus across data subsets). It handles about 1000 times more cells than before. ccfindR is variational Bayesian NMF on counts with gamma priors that selects the rank by maximizing the marginal likelihood (evidence), aimed at finding cancer clones and clusters. Likelihood: Gaussian with uncertainty (CoGAPS), Poisson (ccfindR). K: user-set or tested (CoGAPS), evidence-selected (ccfindR). Multi-sample: none. Tumor: ccfindR named for cancer clones.
- **Strengths:** Bayesian uncertainty and principled rank choice; CoGAPS has strong documentation and protocols.
- **Limitations:** CoGAPS is slow compared with modern optimizers. ccfindR focuses on clustering and has small-scale tooling. Neither models patient or CNV structure.
- **Relevance:** Evidence-based K selection (ccfindR) and consensus over data subsets (CoGAPS) are direct, principled answers to the rank and seed dependence the user complains about.

### GBCD / EBMF (generalized binary covariance decomposition, empirical Bayes matrix factorization via flashier)

*Liu Y, Carbonetto P, Willwerscheid J, ..., Stephens M. Nature Genetics 57:263-273 (2025); flashier R package (Willwerscheid J, Carbonetto P, Stephens M; CRAN)* — ✅ verified online — <https://www.nature.com/articles/s41588-024-01997-z>

- **Core idea:** Factorizes the cell-by-cell covariance as XXᵀ ≈ LLᵀ, using a 'generalized binary' empirical Bayes prior on memberships L (near 0/1 but allowing continuous values) fitted with flashier. It then estimates per-program gene signatures with uncertainty (log fold-changes and lfsr). Applied jointly to malignant cells from many tumors without any integration, it separates patient-specific, dataset-specific, and shared GEPs. Tested on HNSCC and three PDAC studies, it refined PDAC subtypes and found a stress-related program prognostic of survival independent of stage and subtype. The authors show that NMF on pooled malignant cells mostly finds patient-specific programs and misses shared ones. Likelihood: Gaussian on log-transformed counts. K: adaptive through EB shrinkage with a user-supplied maximum. Multi-sample: implicit, because patient effects become their own factors. Interpretability: high. Scale: covariance-based, so atlas scale (millions of cells) is unclear (not re-verified). Tumor MPs: yes, the explicit target.
- **Strengths:** Statistically principled, needs no integration step, and separates shared from patient-specific programs. Clinically validated (prognostic program) and published in a high-impact journal.
- **Limitations:** Gaussian on transformed data rather than counts. Scaling is likely limited for very large pan-cancer cohorts. It does not use prior knowledge, does not model CNV or technical nuisance explicitly, and has no projection or reference tool for new cohorts.
- **Relevance:** The most direct principled competitor to Gavish-style MPs for multi-tumor malignant cells. A new package must benchmark against GBCD and go beyond it on count likelihood, scale, CNV and artifact modeling, priors, and calibrated recurrence.

### sciRED (single-cell interpretable residual decomposition)

*Pouyabahar D, Andrews T, Bader GD. Nature Communications 16 (2025), s41467-025-57157-2* — ✅ verified online — <https://www.nature.com/articles/s41467-025-57157-2>

- **Core idea:** First regresses user-specified technical covariates (library size, sample, protocol) out of counts with a Poisson GLM, then factorizes the Pearson residuals with PCA plus varimax rotation. Factor-covariate association metrics (FCAT) and interpretability scores (separability, effect size, homogeneity, bimodality) label factors as covariate-linked or 'unexplained' biology. Likelihood: Poisson GLM for nuisance, Gaussian PCA for factors. K: user-set. Multi-sample: covariates regressed out. Interpretability: designed for it. Scale: good. Tumor MPs: no.
- **Strengths:** Explicit handling of technical confounders plus quantitative interpretability diagnostics.
- **Limitations:** Factors are signed, not non-negative. Regressing out sample also removes real inter-patient biology. K is user-chosen and factors carry no uncertainty.
- **Relevance:** Shows how to remove nuisance covariates inside a count GLM before factorization; the same slot could take inferred CNV dosage, cell quality, or dissociation stress. Its factor-interpretability metrics could inform an 'MP quality' or artifact score.

### MetaTiME (pan-cancer tumor-microenvironment meta-components via ICA)

*Zhang Y, Xiang G, et al. (Dana-Farber / Harvard). Nature Communications 2023, s41467-023-38333-8* — ✅ verified online — <https://www.nature.com/articles/s41467-023-38333-8>

- **Core idea:** Runs ICA on about 1.7M tumor-microenvironment cells from many scRNA-seq datasets across cancer types. Stable components are consolidated across datasets into annotated meta-components (cell types, states, signaling activities), and new data can be projected into this space for annotation. Likelihood: none (ICA on normalized data). K: per-dataset component number, then clustering. Multi-sample: per-dataset decomposition, then cross-dataset consolidation. Interpretability: annotated meta-components. Scale: million-cell. Tumor: microenvironment, not malignant cells.
- **Strengths:** Large scale, a reusable projection resource, and pan-cancer coverage.
- **Limitations:** Same 'decompose per dataset, then cluster' paradigm as MPs, with the same threshold dependence. Components are signed and not count-based.
- **Relevance:** Shows the meta-component approach and a projection resource at million-cell scale. A new MP method could unify malignant MPs and microenvironment meta-components in one model.

### Contrastive factorization: contrastiveVI and background-contrastive NMF (bcNMF)

*Weinberger E, Lin C, Lee SI. Nature Methods 20:1336-1345 (2023) (contrastiveVI); Li Y, Yang AY, Li Y. arXiv 2602.22387 (2026) (bcNMF; li-lab-mcgill)* — ✅ verified online — <https://www.nature.com/articles/s41592-023-01955-3>

- **Core idea:** contrastiveVI is an NB VAE with background latents shared by target and background cells, plus 'salient' latents active only in target cells. bcNMF jointly factorizes a target dataset and a matched background with shared non-negative bases under a contrastive objective that suppresses structure expressed in the background. It uses GPU minibatch multiplicative updates and was applied to disease brain scRNA-seq, leukemia treatment, and TP53-dependent drug responses. Likelihood: NB (contrastiveVI), Frobenius NMF (bcNMF). K: user-set. Multi-sample: via the target/background design. Tumor MPs: not specifically.
- **Strengths:** Principled separation of shared from enriched variation; bcNMF keeps NMF interpretability.
- **Limitations:** Needs a suitable matched background, such as normal cells of the same lineage, which is not always available. contrastiveVI's salient latents are not linear programs. bcNMF is an early preprint.
- **Relevance:** Gavish et al. found that many carcinoma MPs also vary in normal epithelium. A contrastive tumor-versus-matched-normal-lineage design could separate cancer-specific ITH programs from inherited lineage variability. No MP pipeline offers this axis today.

### Gene programs from single-cell foundation-model gene embeddings (scGPT, Geneformer, scFoundation, UCE, GenePT) and critiques

*Cui H, Wang C, Maan H, Pang K, Luo F, Duan N, Wang B. Nature Methods 21:1470-1480 (2024) (scGPT); Theodoris CV, et al. Nature 618 (2023) (Geneformer); Hao M, et al. Nature Methods 21:1481-1491 (2024) (scFoundation); Rosen Y, Roohani Y, Agarwal A, Samotorčan L, Tabula Sapiens Consortium, Quake SR, Leskovec J. bioRxiv 2023, journal version 2026 (nature.com s41588-026-02718-4; venue appears to be Nature Genetics) (UCE); Chen Y, Zou J. Nature Biomedical Engineering 2025 (GenePT); critiques: Kedzierska KZ, Crawford L, Amini AP, Lu AX. Genome Biology 2025; Kendiukhov I. arXiv 2602.17532 (2026); CancerFoundation bioRxiv 2024 (10.1101/2024.11.01.621087)* — ✅ verified online — <https://www.nature.com/articles/s41592-024-02201-0>

- **Core idea:** scGPT builds a similarity graph on its learned gene embeddings and clusters it with Leiden into 'gene programs' (for example HLA and CD gene programs) with cell-type-specific activation. Geneformer, scFoundation, UCE (which uses ESM2 protein embeddings), and GenePT (GPT text embeddings of NCBI/UniProt gene summaries) also provide gene embeddings usable as gene-gene similarity priors. CancerFoundation is a small foundation model trained only on malignant cells (about 1M cells, more than 40 tumor types). Likelihood: none for program discovery. K: set by Leiden resolution. Multi-sample: implicit through pretraining. Interpretability: moderate. Tumor MPs: not targeted.
- **Strengths:** Cheap, broad prior knowledge, including for genes poorly covered by curated databases.
- **Limitations:** Programs from embeddings are averaged over contexts and do not reflect co-variation within a given tumor. Zero-shot foundation-model embeddings underperform HVG, Harmony, and scVI baselines (Kedzierska 2025). Attention largely reflects co-expression and adds no value over gene-level baselines for regulatory inference (Kendiukhov 2026). Most models are pretrained on non-malignant cells.
- **Relevance:** Best used as an optional soft prior (gene-gene graph or embedding regularizer) or for annotation, not as the source of MPs. Reviewers will demand ablations showing any foundation-model component beats PPI, co-expression, or random-graph priors.

### Sparse autoencoders (SAEs) on single-cell foundation models

*Pedrocchi F, Barkmann F, Joudaki A, Boeva V. bioRxiv 2025 (10.1101/2025.10.22.681631), ICLR 2026; 'Sparse autoencoders reveal interpretable cell-type programs in single-cell foundation model representations', Journal of Biomedical Informatics 2026 (authors not retrieved); arXiv 2603.02952 (Geneformer/scGPT SAE atlas, 2026)* — ✅ verified online — <https://openreview.net/forum?id=5AuN6RD072>

- **Core idea:** Train sparse autoencoders on intermediate token or residual-stream activations of scGPT, scFoundation, and Geneformer. The resulting sparse features capture biological programs (pathways, cell-type programs) and also technical signals such as batch. Steering SAE features can reduce technical effects while keeping biological signal, and later layers carry more biologically coherent features. Likelihood: reconstruction of activations. K: dictionary size with L1/top-k sparsity. Multi-sample: technical features found implicitly. Tumor MPs: no.
- **Strengths:** Novel mechanistic-interpretability angle and can explicitly separate technical from biological features.
- **Limitations:** Features describe the model rather than the dataset. Running foundation-model forward passes over millions of cells is costly. Features can split or absorb one another. Validation is mostly enrichment-based, and malignant ITH has not been tested.
- **Relevance:** A possible source of learned 'technical artifact detectors' (batch, stress, quality) that could serve as nuisance covariates for an MP model. Not yet a program discovery method in its own right.

### Tripso (program-centric self-supervised transformer)

*Lotfollahi lab (Wellcome Sanger Institute). bioRxiv March 2026 (10.64898/2026.03.24.713961), 'Self-supervised learning for a gene program-centric view of cell states'; full author list not retrieved* — ✅ verified online — <https://tripso-docs.cog.sanger.ac.uk/index.html>

- **Core idea:** A self-supervised transformer that represents each cell with several gene-program-specific embeddings instead of one entangled embedding. It learns contextualized gene, program, and cell embeddings, quantifies each gene's contribution to a program and each program's influence on cell identity, and supports discovery of new programs. Examples include JAK-STAT activity in pediatric hematopoiesis and a new tissue-resident memory T cell program in atopic dermatitis validated spatially. Likelihood: self-supervised transformer objective. K: number of programs; to my understanding initialized from defined gene programs (not re-verified). Tumor MPs: no.
- **Strengths:** Context-dependent gene weighting with program-level cell representations.
- **Limitations:** Heavy compute, unpublished preprint, and not applied to malignant ITH. Programs depend on how the initial programs are defined.
- **Relevance:** Defines the 'AI program discovery' competitor a Nature Methods editor will think of. An MP package should explain why a statistically principled, count-native, multi-tumor model is better suited to reproducible MP definition.

### DeltaNMF (two-stage neural NMF with foundation-model gene-similarity regularization)

*Karpurapu A, Gersbach C, Singh R. bioRxiv January 2026 (10.64898/2026.01.22.701049); listed in RECOMB 2026 materials (not independently confirmed)* — ✅ verified online — <https://www.biorxiv.org/content/10.64898/2026.01.22.701049v1>

- **Core idea:** Stage 1 learns baseline programs from control cells. Stage 2 learns case-specific programs, separating changed usage of existing programs from genuinely new programs. A graph-Laplacian regularizer built from foundation-model-derived gene similarities improves program coherence 2.4-fold on protein-interaction networks. The GPU implementation is about 20 times faster than cNMF with similar accuracy. Likelihood: NMF objective. K: user-set per stage. Multi-sample: case/control. Tumor MPs: no.
- **Strengths:** Simple, effective prior injection, fast, and interpretable.
- **Limitations:** Designed for case-control or perturbation settings, with no multi-tumor hierarchy, CNV handling, or calibrated recurrence. Early preprint.
- **Relevance:** Direct precedent for foundation-model embeddings as soft graph priors inside NMF, plus GPU speed. An MP package must go clearly beyond it (tumor hierarchy, nuisance modeling, uncertainty) to claim novelty on the prior-informed axis.

### Evidence on tumor-specific confounding: UVG, patient-unique malignant identity, CNV-expression contrastive alignment (CLCC), KINOMO

*Yoon SH, Nam JW. Briefings in Bioinformatics 25(1):bbad460 (2024) (UVG); 'A unique malignant cell type per patient tumor encoded in each cancer cell transcriptome', iScience 2026 (bioRxiv 10.1101/2024.10.30.620989); CLCC: 'Contrastive Alignment of Expression and Copy Number Highlights Dosage-Insensitive Genes in Cancer', bioRxiv 2026 (10.64898/2026.03.01.708901); KINOMO: bioRxiv 2022 (10.1101/2022.05.02.490362), AACR 2024 abstract* — ✅ verified online — <https://academic.oup.com/bib/article/25/1/bbad460/7469351>

- **Core idea:** UVG normalizes gene variances to pick 'universally variable genes', which reduces sample-specific clusters of malignant cells. The iScience study, covering 160 tumors and 9 cancer types, shows that patient-specific malignant identity is spread across the whole transcriptome, is not caused by technical bias, and cannot be reduced to tumor-specific gene sets. CLCC aligns expression with inferred CNV profiles contrastively to classify CNV-concordant versus discordant cells and genes. KINOMO is a KL-divergence NMF with kernel-based correction aimed at inter-patient genomic variability.
- **Strengths:** Recent, quantitative evidence that patient and CNV effects are pervasive and that heuristic gene selection only partly helps.
- **Limitations:** None of these is a complete MP discovery framework with CNV modeled as an explicit nuisance (KINOMO comes closest, as an unpublished preprint).
- **Relevance:** Strong motivation for modeling patient identity and CNV dosage explicitly inside the factor model, rather than avoiding them through per-tumor analysis, which also gives up borrowing strength across tumors.

**Gaps identified by this sweep**

- No method combines a count likelihood (Poisson/NB) with an explicit shared-vs-patient-specific program hierarchy (iNMF-like W + V_i), automatic choice of K (ARD, horseshoe, or empirical Bayes), and calibrated uncertainty at pan-cancer scale (thousands of tumors, millions of malignant cells). GBCD is principled but Gaussian and covariance-based. scDEF is count-based, hierarchical, and picks K automatically, but it is a small-scale preprint. iNMF has the shared/specific structure but no count likelihood and a hand-picked K.
- No program discovery method uses inferred copy-number profiles (inferCNV, CopyKAT, Numbat) or subclone structure as an offset or covariate, so CNV dosage effects cannot be separated from real transcriptional programs. sciRED-style GLM nuisance regression exists but has never been applied to CNV, and UVG and KINOMO are partial heuristics.
- No method models technical artifacts (dissociation or stress response, ambient RNA, mitochondrial and ribosomal content, doublets, low-quality cells, sequencing depth or chemistry) as dedicated nuisance factors with priors. As a result, artifact MPs cannot be flagged in a principled way; SAE-detected technical features are only a hint.
- There is no statistical definition or test of program 'recurrence across tumors'. Current MPs rely on fixed overlap cutoffs (top-50 genes, 20% overlap, minimum cluster size) with no null model that accounts for gene expression level, HVG selection, or gene-set size. Nobody reports per-MP or per-gene confidence such as posterior inclusion probabilities or bootstrap stability over tumors, cells, seeds, and ranks.
- Robustness is not evaluated as a first-class output. No tool reports, for each MP, its stability under resampling tumors and cohorts, changing seeds and K, or splitting by technology or platform. Reproducibility in held-out cohorts is shown informally at best (starCAT is the closest precedent, but for T cells).
- Prior knowledge is used either as hard architecture (expiMap, VEGA) or a strong penalty (Spectra), or not at all. Soft, ablatable priors from pathway graphs or foundation-model/text gene embeddings exist only in narrow settings (DeltaNMF for case-control, scE2TM for topic models), and none come with a per-program 'novelty versus prior dependence' score in multi-tumor MP discovery.
- No contrastive design separates cancer-specific ITH programs from variability inherited from the normal cell of origin, for example malignant cells versus matched normal epithelium from the same lineage, even though Gavish et al. showed that most carcinoma MPs are also present in normal epithelium. contrastiveVI and bcNMF provide the machinery but have not been applied this way.
- Multi-resolution and context-dependent MPs (generic versus lineage-specific variants such as EMT-like or stress subtypes) are not modeled jointly. Hierarchical models (scDEF, scGEM nested HDP) have not been scaled to pan-cancer malignant data.
- No ground-truth benchmark exists for metaprograms. A useful suite would include: semi-synthetic tumors with planted shared and patient-specific programs plus CNV and batch confounding; held-out cohorts; Perturb-seq or drug-perturbation recovery; spatial co-localization; lineage or clone tracing; and clinical association. The Engreitz-lab gene_network_evaluation framework and Spectra or GBCD simulations are partial starts.
- There is no count-model projection or reference-mapping tool for malignant MPs comparable to starCAT for T cells, with uncertainty on per-cell MP usage and detection of new programs in a query cohort.
- Malignant MPs and microenvironment states are not modeled jointly (Barkley 2022 and MetaTiME analyze them separately), so links between cancer-cell programs and immune or stromal states are inferred after the fact.
- Scalability: count-native, hierarchical Bayesian factorization with GPU stochastic VI over about 10M malignant cells from about 3,000 tumors (3CA scale) has not been demonstrated. Existing consensus methods multiply cost by the number of replicates.
- LLM-based program naming (GeneAgent, scConcept) is not built into MP pipelines with verification, so MP annotation remains manual and inconsistent across studies.

**Notes:** LANDSCAPE SUMMARY (as of 2026-10-02):
(1) The tumor-MP standard (Gavish 2023 Nature; 3CA update by Tyler et al., Nature Cancer 2025; GeneNMF v0.6 cosine-based consensus; SSpMosaic, Cell Genomics 2025) is still a threshold-based heuristic built on per-tumor NMF.
(2) The most principled direct competitor is GBCD (Liu, Carbonetto, Willwerscheid, ..., Stephens; Nature Genetics 2025). It was built specifically to separate shared from patient-specific programs across tumors without integration, and the authors show that NMF on pooled malignant cells misses shared programs. Any new package must benchmark against GBCD, plus cNMF, GeneNMF, consensus scHPF, consensus iNMF, scDEF and DRVI.
(3) The best Nature Methods precedent for a 'reproducible program catalog + projection' paper is starCAT/TCAT (Kotliar et al., Nature Methods 2025), but it covers T cells only.
(4) Single-cell foundation models give gene embeddings and SAE features, but zero-shot they underperform simple baselines (Kedzierska 2025, Genome Biology). Attention mostly reflects co-expression (Kendiukhov 2026). They should be optional soft priors with ablations, not the core. DeltaNMF (2026 preprint) already uses foundation-model gene-similarity Laplacian priors in NMF, so 'FM-informed NMF' alone is not novel.

LIKELY WHITE SPACE FOR A NATURE METHODS-LEVEL PACKAGE: a count-native (Poisson/NB) hierarchical Bayesian non-negative factor model combining these components:
- shared MP loadings plus sparse tumor-specific deviations (iNMF/GBCD-like);
- automatic K via ARD, horseshoe, or beta-process priors;
- explicit nuisance terms for inferred CNV dosage (offset from inferCNV/Numbat), cell quality and depth, and an artifact library (dissociation stress, ambient RNA);
- an optional soft gene-graph prior (pathway, PPI or foundation-model embeddings) with a per-program prior-dependence or novelty score;
- posterior inclusion probabilities for MP gene membership and a statistical recurrence test across tumors;
- an optional contrastive term against matched normal cell-of-origin;
- GPU stochastic VI to reach 3CA scale (about 2,800 samples);
- a starCAT-like projection tool with per-cell uncertainty.

Pair the package with a benchmark suite (semi-synthetic planted programs with CNV and batch confounding, cohort-split reproducibility, seed/rank stability, perturbation recovery, spatial validation, prognostic replication). Re-deriving 3CA MPs and flagging the unstable or artifactual ones would be the flagship biological result.

VERIFICATION CAVEATS:
- WebFetch was blocked by the egress proxy for biorxiv.org, pubmed and github.io, so most details come from search-result snippets and publisher landing-page summaries.
- Details I could not fully re-check: Spectra's exact loss/likelihood; expiMap's default reconstruction loss; GBCD's maximum practical scale; whether Tripso initializes from predefined programs; the full scHPF, Spectra, MetaTiME and Tripso author lists (abbreviated with 'et al.' on purpose); the authors of the 2026 J Biomed Inform SAE paper.
- UCE's journal venue: the URL prefix s41588 indicates Nature Genetics, while a snippet said 'Nature'.
- DeltaNMF appears in RECOMB 2026 listings, but I did not confirm this independently.
- scDEF, consensus scHPF, KINOMO, Tripso, scConcept, bcNMF, CLCC and DRVI were not found as journal publications; they appear to be preprints, workshop papers or software only.
- All items marked verified_online=true were found online in this session.


## Statistical foundations (23 items)

### BMSFA: Bayesian multi-study factor analysis (shared + study-specific loadings)

*De Vito R, Bellio R, Trippa L, Parmigiani G. Annals of Applied Statistics 15(4), 2021. DOI 10.1214/21-AOAS1456* — ✅ verified online — <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-15/issue-4/Bayesian-multistudy-factor-analysis-for-high-throughput-biological-data/10.1214/21-AOAS1456.full>

- **Core idea:** Each study s is modelled as x = Phi f + Lambda_s l_s + e. Phi holds loadings shared by every study and Lambda_s holds loadings unique to study s. Multiplicative-gamma shrinkage priors (a multi-study extension of the sparse Bayesian infinite factor model) learn how many shared and how many specific factors there are. Fitted with a Gibbs sampler; the paper's application was 7 breast-cancer expression studies.
- **Strengths:** The first principled Bayesian treatment of separating shared from study-specific factors. The number of factors is learned rather than fixed, and the shared factors replicated across studies.
- **Limitations:** Gaussian likelihood with real-valued loadings: not NMF and not count-based. A factor is either shared by all studies or private to one, with nothing in between. Gibbs sampling scales poorly to about 20k genes and to hundreds of studies or tumours. Rotational non-identifiability needs post-processing (varimax/Procrustes).
- **Relevance:** Treat each tumour as a 'study'. MPs become columns of Phi. Patient-specific programs, for example CNV-driven blocks, are absorbed by Lambda_s, so they never get counted as recurrent. This addresses the user's complaint that patient-specific programs leak into MPs. The all-or-one sharing assumption does not fit MPs, which recur in subsets of tumours; that motivates Tetris (next item). In my searches I found no application of BMSFA to tumour scRNA-seq MPs.

### Tetris: Bayesian combinatorial multi-study factor analysis (IBP over sharing patterns)

*Grabski IN, De Vito R, Trippa L, Parmigiani G. Annals of Applied Statistics 17(3), 2023. DOI 10.1214/22-AOAS1715 (arXiv 2007.12616)* — ✅ verified online — <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-3/Bayesian-combinatorial-MultiStudy-factor-analysis/10.1214/22-AOAS1715.short>

- **Core idea:** Latent factors can be shared by any subset of studies. A binary factor-by-study 'ownership' matrix gets an Indian buffet process prior, which infers the number of factors and their sharing patterns together. Uncertainty in the sharing pattern is summarised with credible balls.
- **Strengths:** Matches the definition of an MP as 'a program recurring in a subset of tumours'. It gives a posterior over which studies own each factor, plus an uncertainty summary.
- **Limitations:** Gaussian rather than count or non-negative data. Computationally heavy: the 2026 Stat Med tutorial's repository says 'Tetris takes far more time to run'. Shown only on modest numbers of studies and features. Rotation and label-switching issues remain.
- **Relevance:** Conceptually the closest statistical template for replacing greedy overlap clustering of NMF programs. One would fit a single joint model in which an MP is a factor plus a posterior set of tumours expressing it, with credible sets for membership. A Poisson-NMF version that scales by stochastic variational inference to hundreds or thousands of tumours, with genomic and QC covariates, would be new. I found no application to tumour scRNA MPs.

### Bayesian multi-study NMF for mutational signatures (discovery-only and recovery-discovery), plus probit multi-study NMF follow-up

*Grabski IN, Trippa L, Parmigiani G. Genome Biology 26:98, 2025 (DOI 10.1186/s13059-025-03563-0). Follow-up: 'Bayesian Probit Multi-Study Non-negative Matrix Factorization for Mutational Signatures', arXiv 2502.01468, 2025 (authors not confirmed in this session).* — ✅ verified online — <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03563-0>

- **Core idea:** Poisson NMF fitted jointly across studies (cancer types). A binary factor-inclusion matrix learns which signatures are present in which study and the total number of signatures. Two modes: a fully unsupervised 'discovery-only' model, and a 'recovery-discovery' model that estimates known signatures (informative priors) and novel ones together. Both extend to covariate effects on exposures. The probit follow-up changes how inclusion is modelled.
- **Strengths:** The only count-native, non-negative, multi-study model with subset sharing and a learned number of factors that I found peer-reviewed. Recovery-discovery is a principled way to combine prior knowledge with new discovery.
- **Limitations:** Built for about 96 mutation channels and hundreds of samples, not about 20k genes by millions of cells. MCMC would not scale to tumour scRNA without big engineering (SVI, minibatching, gene filtering). The Poisson likelihood ignores overdispersion.
- **Relevance:** A direct blueprint for 'multi-tumour NMF'. In recovery-discovery, the 41 Gavish MPs could serve as priors while novel MPs are discovered. Covariates could absorb nuisance: inferCNV profiles, technology, dissociation protocol, QC metrics. Moving this framework from mutational signatures to scRNA-seq (cells nested in tumours nested in cohorts, NB likelihood, scalable inference) is a clear unbuilt gap.

### SUFA: subspace factor analysis (identifiability of shared vs study-specific covariance)

*Chandra NK, Dunson DB, Xu J. Journal of the American Statistical Association, 2025 (arXiv 2305.04113). R package SUFA (HMC-within-Gibbs, C++/OpenMP).* — ✅ verified online — <https://github.com/noirritchandra/SUFA>

- **Core idea:** Study-specific loadings are constrained to lie in the column space of the shared loading matrix. That constraint makes the split between shared and study-specific covariance components identifiable in high dimensions.
- **Strengths:** Tackles the identifiability problem directly: plain MSFA models can trade variance between the shared and specific parts. Comes with theory and a fast sampler.
- **Limitations:** The subspace assumption is biologically questionable for tumours. Truly private programs, for example from amplicons, need not lie in the shared span. Gaussian, not non-negative.
- **Relevance:** Supplies identifiability theory and a cautionary point. Any 'shared vs tumour-specific' model must say what makes the split identifiable, or 'shared' MPs will depend on the prior and initialisation, which is exactly the robustness complaint. A new method should state an explicit identifiability condition, for example anchor genes, sparsity, or a non-negative sufficiently-scattered condition across tumours.

### PFA: perturbed factor analysis

*Roy A, Lavine I, Herring AH, Dunson DB. Annals of Applied Statistics 15(3), 2021. DOI 10.1214/20-AOAS1435* — ✅ verified online — <https://projecteuclid.org/journals/annals-of-applied-statistics/volume-15/issue-3/Perturbed-factor-analysis--Accounting-for-group-differences-in-exposure/10.1214/20-AOAS1435.full>

- **Core idea:** All factors are common, but each study's data is first multiplied by a study-specific perturbation matrix near the identity. A divergence statistic computed from the perturbation matrices quantifies how far each study, or pair of studies, departs from the common signal.
- **Strengths:** Produces one interpretable set of factors that generalise across studies, together with a quantitative per-study distortion measure.
- **Limitations:** Gaussian. Perturbations are dense p-by-p matrices, so it scales poorly to many genes. It does not represent genuinely private programs.
- **Relevance:** Useful for platform and cohort effects in 3CA-type compendia (10x vs Smart-seq2, fresh vs frozen, dissociation). Model these as gene-level perturbations of a shared MP loading rather than letting them produce new 'MPs'. The divergence statistic could flag outlier datasets or tumours before MP calling. Not applied to scRNA-seq.

### APAFA: adaptive partition factor analysis

*Bortolato E, Canale A. Journal of the American Statistical Association, 2026 (online; DOI 10.1080/01621459.2026.2621517); arXiv 2410.18939 (v3, Jan 2026). Code: github.com/elenabortolato/APAFA* — ✅ verified online — <https://www.tandfonline.com/doi/full/10.1080/01621459.2026.2621517>

- **Core idea:** New shrinkage priors for multi-study factor analysis cover the full range: no study-specific factors at all, factors private to small subgroups nested within studies, and factors shared across studies. The goal is explicitly to separate shared latent traits from 'artifactual or population-specific' variation. Gives identifiability conditions and a Gibbs sampler.
- **Strengths:** The most recent statistics paper aimed squarely at telling artifactual or subgroup factors apart from shared ones. Allows sub-study partitions, which matter because subclones inside a tumour are nested subgroups. Provides identifiability conditions.
- **Limitations:** Gaussian and Gibbs-based. Applications are bird co-occurrence and ovarian-cancer gene expression, not single-cell. Scalability to scRNA-seq is undemonstrated.
- **Relevance:** Its partition idea maps naturally onto cells nested in subclones, nested in tumours, nested in cohorts. It suggests a prior in which a program can be subclone-private, tumour-private, cohort-private (probably technical) or cross-cohort (an MP candidate). That gives a statistical definition of 'real' recurrent biology versus cohort artifact. Not done for scRNA-seq.

### Benchmark and tutorial of Bayesian multi-study FA methods (PFA, MOM-SS, SUFA, BMSFA, CAVI-BMSFA, BLAST, Tetris)

*Liang M, et al. 'A Tutorial on Bayesian Multi-Study Factor Analysis With Applications in Nutrition and Genomics', Statistics in Medicine 45(10-12), 2026 (DOI 10.1002/sim.70531). R package bmfaToolkits. Related scalable methods: 'Spectral decomposition-assisted multi-study factor analysis' (BLAST, arXiv 2502.14600) and 'Fast Variational Inference for Bayesian Factor Analysis in Single and Multi-Study Settings' (arXiv 2305.13188). Author lists for the two arXiv papers were not confirmed in this session.* — ✅ verified online — <https://onlinelibrary.wiley.com/doi/10.1002/sim.70531>

- **Core idea:** A head-to-head simulation comparison of seven Bayesian integrative factor models on computing cost, loading accuracy and recovery of the number of factors, with a full R workflow and code.
- **Strengths:** A ready-made baseline suite and simulation designs for the 'shared vs specific' part of a new method. It also shows which methods recover the number of factors and which are impractically slow (Tetris).
- **Limitations:** Gaussian models on bulk-scale data (nutrition and genomics). Nothing on count data or single-cell scale. Only the tutorial's existence and method list were confirmed; detailed results were not read in this session.
- **Relevance:** For a Nature Methods paper, reviewers will ask why you didn't just use MSFA or Tetris. This toolkit lets you show where they fail on scRNA-seq (scale, counts, non-negativity), which justifies a new count-based, scalable multi-tumour model. BLAST and CAVI are the scalable comparators to beat.

### Nonlinear multi-study factor analysis (multi-study sparse VAE)

*Moran GE, Krishnan A. ICASSP 2026; arXiv 2601.18128 (Jan 2026)* — ✅ verified online — <https://arxiv.org/abs/2601.18128>

- **Core idea:** A sparse variational autoencoder whose decoder lets each gene depend on few latent factors. Factors may be shared across disease groups or specific to one. Sparsity implicitly penalises the number of factors and helps separate shared from group-specific ones. Applied to platelet gene expression across diseases.
- **Strengths:** Amortised inference scales to large numbers of samples or cells. The sparse decoder gives identifiability-style guarantees, building on earlier sparse-VAE identifiability work by Moran and colleagues (recalled, not re-verified here). Nonlinear.
- **Limitations:** A short conference paper on bulk data. Not non-negative or NMF-style. Interpreting nonlinear factors as gene lists is less direct, and there is no FDR or uncertainty on sharing.
- **Relevance:** A plausible engine for a scalable 'multi-tumour' model: the encoder handles millions of cells and the sparse decoder yields gene-program loadings. Its group structure could be tumour, cohort or technology. Not yet applied to tumour scRNA-seq.

### GBCD: generalized binary covariance decomposition (multi-tumour scRNA: patient-specific, dataset-specific, shared GEPs)

*Liu Y, Carbonetto P, Stephens M. Nature Genetics 57:263-273, 2025 (bioRxiv 10.1101/2023.08.15.553436). R package stephenslab/gbcd* — ✅ verified online — <https://www.nature.com/articles/s41588-024-01997-z>

- **Core idea:** Pools malignant cells from many tumours and decomposes the cell-by-cell covariance under a 'generalized binary' prior on cell memberships, so each GEP is roughly on or off in groups of cells. Components come out patient-specific, dataset-specific or shared. The authors show that pooled NMF mostly finds patient-specific GEPs and misses shared ones. Applied to HNSCC and PDAC, where it found a stress-related GEP prognostic of survival.
- **Strengths:** The most direct published competitor that statistically separates shared from patient-specific programs in tumour scRNA-seq without per-tumour NMF followed by heuristic clustering. Published in a high-impact journal with software.
- **Limitations:** These points are my reading of the method; the paper's own limitations section could not be opened (domain blocked). It works on transformed (log-normalised) data rather than counts. The generalized-binary assumption may fit continuous gradients poorly. Component number and labels come from greedy empirical-Bayes fitting. There is no formal error rate or FDR for 'shared' calls and no uncertainty on which tumours share a GEP. Scaling to about 1,000 tumours and over a million cells is undemonstrated as far as I could verify.
- **Relevance:** A must-benchmark baseline. A new package has to beat the Gavish/3CA pipeline, GeneNMF and GBCD on reproducibility, recovery of planted shared programs, and rejection of planted patient-specific or artifact programs. GBCD also shows Nature Genetics or Nature Methods-level interest in this problem.

### Shared + dataset/condition-specific NMF for single cells: iNMF/LIGER, UINMF, scINSIGHT

*Welch JD et al., Cell 177:1873-1887, 2019 (LIGER/iNMF). Kriebel AR, Welch JD, 2022 (UINMF, per the LIGER repo). Qian K, Li WV, Genome Biology 2022 (scINSIGHT, DOI 10.1186/s13059-022-02649-3)* — ✅ verified online — <https://github.com/welch-lab/liger>

- **Core idea:** iNMF factorises each dataset as (W + V_i) H_i: W is shared and V_i are dataset-specific gene loadings, with a penalty lambda controlling how much they may differ. scINSIGHT adds sparse non-negative gene modules that are common to all conditions or specific to one condition.
- **Strengths:** Count-adjacent non-negative models already used at single-cell scale, with mature software.
- **Limitations:** In iNMF, V_i is a deviation of the same k factors, not a separate private program. k and lambda are tuned heuristically. Nothing models sharing by a subset of tumours, and there is no uncertainty or FDR. The goal is integration (alignment), not discovering recurrent programs. scINSIGHT needs known condition labels and assumes modules are common to all or specific to one condition.
- **Relevance:** The natural 'existing-tool' baselines for one-model-across-tumours MP discovery. They show that single-cell-scale shared or specific NMF is technically feasible. A new package could generalise them with a subset-sharing (IBP or Tetris-like) structure and count likelihoods.

### Multi-group latent variable models: MOFA+ multi-group ARD, multiGroupVI

*Argelaguet R et al. 'MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data', Genome Biology 21:111, 2020 (from memory; not re-verified in this session). Weinberger E et al. 'Disentangling shared and group-specific variations in single-cell transcriptomics data with multiGroupVI', bioRxiv 10.1101/2022.12.13.520349, 2022 (title and URL confirmed; author list from memory).* — ⚠️ not verified — <https://www.biorxiv.org/content/10.1101/2022.12.13.520349v1>

- **Core idea:** MOFA+ places automatic relevance determination (ARD) priors on factor activity per sample group, so a factor can be active in only some groups, and fits by stochastic variational inference. multiGroupVI is a VAE with one shared latent space plus group-specific latent spaces.
- **Strengths:** Scalable (SVI or amortised inference) and widely used in single-cell work.
- **Limitations:** MOFA+ uses Gaussian or Poisson likelihoods with real-valued, non-NMF loadings, and group-wise ARD gives factor variance per group, not calibrated presence calls. multiGroupVI's latent spaces are not interpretable gene programs, and it handles few groups, not about 1,000 tumours. Neither gives FDR for 'shared'.
- **Relevance:** Shows that group-wise ARD (factor-by-tumour relevance) is practical at scale. A metaprogram model could use ARD or spike-and-slab per (program, tumour), which is essentially a scalable Tetris. MOFA+ details here are not verified online in this session.

### NMF identifiability and uniqueness theory (sufficiently scattered condition, minimum/maximum-volume NMF)

*Fu X, Huang K, Sidiropoulos ND. 'On Identifiability of Nonnegative Matrix Factorization', IEEE Signal Processing Letters 2018 (arXiv 1709.00614). Fu X, Huang K, Sidiropoulos ND, Ma W-K. 'Nonnegative Matrix Factorization for Signal and Data Analytics: Identifiability, Algorithms, and Applications', IEEE Signal Processing Magazine 2019 (arXiv 1803.01257). The arXiv entries were found; venues and the full author list of the second paper are from memory. Also 'Identification of NMF by choosing maximum-volume basis vectors', arXiv 2603.24227 (2026; authors not confirmed).* — ✅ verified online — <https://arxiv.org/pdf/1803.01257>

- **Core idea:** NMF is unique up to permutation and scaling only under geometric conditions. The separability or anchor condition is the strict version. The sufficiently scattered condition (SSC) is milder: rows of one factor must spread well over the non-negative orthant. Volume-regularised (min-vol or max-vol) criteria restore identifiability under SSC-type conditions.
- **Strengths:** Gives a precise mathematical account of why NMF programs move with seed and rank, and a principled cure.
- **Limitations:** Theory is for noiseless or low-noise models. Hard to check SSC on real noisy scRNA-seq. I found no single-cell-specific identifiability analysis.
- **Relevance:** Explains the seed and rank instability. Tumour cells often express overlapping, correlated programs and few 'pure' cells exist, so SSC fails and NMF has a continuum of equivalent solutions. Top-50 gene lists then change between runs. A new package could (a) use volume-regularised or anchor-gene NMF to make programs identifiable, (b) report a per-tumour, per-program identifiability diagnostic such as near-pure cells and anchor genes, and (c) argue that pooling across tumours restores SSC (tumours differ in state mix). As far as I found, (c) is a new theoretical angle.

### staNMF: stability-driven NMF model selection (and cNMF's stability/error K selection)

*Wu S, Joseph A, Hammonds AS, Celniker SE, Yu B, Frise E. PNAS 113(16), 2016. Related: Kotliar D et al., eLife 2019 (cNMF; repo confirms the k-selection plot and local-density outlier filter)* — ✅ verified online — <https://www.pnas.org/doi/10.1073/pnas.1521171113>

- **Core idea:** Run NMF many times at each K and pick the K whose dictionaries are most reproducible across random initialisations (an Amari-type instability score). cNMF combines a stability measure with reconstruction error and filters outlier components before forming consensus.
- **Strengths:** Simple, model-agnostic, and part of the PCS / veridical data science tradition from the Yu lab. Produced biologically validated patterns (21 principal patterns in Drosophila).
- **Limitations:** Stability is not validity. Technical programs (ambient RNA, mitochondrial or ribosomal, dissociation stress, cell cycle) are highly stable. Restart-only perturbation ignores sampling variability. There is no error-rate guarantee. Consensus clustering can show spurious stability, a known pitfall (e.g. Senbabaoglu et al., Sci Rep 2014; not re-verified here).
- **Relevance:** Gavish-style pipelines use a fixed rank range (around 4 to 9) plus overlap thresholds. Swapping in per-tumour, stability-chosen K is an easy improvement but not a publishable advance on its own. Stability should be measured across several perturbation types (restarts, cell bootstrap, count thinning, gene subsets, rank) and combined with a null calibration (see the synthetic-null items).

### Data thinning / generalized data thinning (convolution-closed splitting)

*Neufeld A, Dharamshi A, Gao LL, Witten D. 'Data thinning for convolution-closed distributions', JMLR 25(57):1-35, 2024. Dharamshi A, Neufeld A, Motwani K, Gao LL, Witten D, Bien J. 'Generalized data thinning using sufficient statistics', arXiv 2303.12931 (published in JASA, from memory)* — ✅ verified online — <http://www.jmlr.org/papers/volume25/23-0446/23-0446.pdf>

- **Core idea:** Split each observation X into independent parts X1 + X2 that follow the same family with a scaled parameter (Poisson, NB with known overdispersion, Gaussian with known variance, gamma, binomial, ...). Fit on X1 and validate or infer on X2. Works when there is one parameter per observation, as in low-rank matrix estimation, where ordinary sample splitting fails.
- **Strengths:** Valid out-of-sample evaluation for unsupervised matrix models; the papers name choosing the rank of a low-rank approximation as a use case.
- **Limitations:** Needs the nuisance parameters (NB overdispersion, Gaussian variance) to be known or well estimated, otherwise the folds are dependent. Splitting loses signal, which hurts for lowly expressed genes and small tumours.
- **Relevance:** Gives a principled per-tumour rank selection: fit NMF at each K on X1 and choose K by held-out NB deviance on X2. It also enables program validation: a program found on X1 must explain covariance in X2. Repeated multi-fold thinning gives stability across independent noise realisations rather than just restarts. I found no application to tumour MP discovery.

### Count splitting for scRNA-seq latent variable models (Poisson and negative binomial)

*Neufeld A, Gao LL, Popp J, Battle A, Witten D. 'Inference after latent variable estimation for single-cell RNA sequencing data', Biostatistics 25(1):270-287, 2024. Neufeld A, Popp J, Gao LL, Battle A, Witten D. 'Negative binomial count splitting for single-cell RNA sequencing data', arXiv 2307.12985 (2023). R package countsplit (CRAN v4).* — ✅ verified online — <https://github.com/anna-neufeld/countsplit>

- **Core idea:** Thin the cell-by-gene count matrix into independent folds. Estimate latent structure (clusters, pseudotime, PCs or factors) on one fold and test gene associations on the other, which avoids double dipping. The authors also propose fold-wise cross-validation to choose the number of latent dimensions or clusters.
- **Strengths:** Designed for scRNA-seq with a software package. A 2026 Bioinformatics Advances paper ('Quantifying stability via count splitting to guide model selection in RNA velocity analyses') shows count splitting already used for stability-guided model selection in scRNA-seq (authors not confirmed).
- **Limitations:** The NB version needs per-gene overdispersion. If overdispersion is misestimated, the folds are not independent and false positives return. The theory was built for clusters, pseudotime and PCs, not NMF programs or cross-tumour aggregation.
- **Relevance:** Concrete adaptation: (1) count-split every tumour; (2) discover programs and MPs only on the training folds; (3) on the test folds, test each gene's association with the MP score (cell-level usage) to get calibrated p-values and FDR-controlled MP gene sets, replacing the arbitrary 'top 50 genes'; (4) choose per-tumour K and the global MP number by held-out likelihood. I did not find this applied to NMF metaprograms, so it is a publishable methodological contribution.

### Limits of data fission / data thinning for structured (mixture) data

*Hivert B, Agniel D, Thiebaut R, Hejblum BP. 'Running in circles: practical limitations for real-life application of data fission and data thinning in post-clustering differential analysis', arXiv 2405.13591, 2024 (v3)* — ✅ verified online — <https://arxiv.org/abs/2405.13591>

- **Core idea:** Data fission and thinning assume a non-mixture parametric model. For clustered data, the marginal scale or overdispersion estimate is biased relative to the within-component value, which creates dependence between folds and inflates type I error. The correct conditional thinning needs the unknown clustering, which is circular.
- **Strengths:** An important negative result that reviewers will know about.
- **Limitations:** Focuses on post-clustering differential expression with discrete clusters. Continuous-factor (NMF) settings are not analysed.
- **Relevance:** Malignant cells form mixtures of states, so a naive marginal NB fit per gene would mis-specify overdispersion. A credible new method must (a) estimate overdispersion conditionally on a preliminary low-rank fit or use robust estimators, (b) show calibration by simulation and synthetic nulls, or (c) use Poisson thinning on overdispersion-robust statistics. Extending this analysis to factor models is an open question.

### Selective inference after clustering (hierarchical, k-means, unknown variance, pipelines)

*Gao LL, Bien J, Witten D. 'Selective inference for hierarchical clustering', JASA (accepted 2022). Chen YT, Witten D. 'Selective inference for k-means clustering', JMLR 2023 (PMC10805457; authors from memory). Yun Y-J, Barber RF. 'Selective inference for clustering with unknown variance', Electronic Journal of Statistics 17:1923-1946, 2023 (authors from memory). Miyata Y, Shiraishi T, Nishino S, Takeuchi I. 'Statistical Testing Framework for Clustering Pipelines by Selective Inference', arXiv 2603.18413, 2026. Also 'Selective Inference for Deep Clustering in Latent Spaces', arXiv 2609.28756, 2026 (authors not confirmed).* — ✅ verified online — <https://arxiv.org/pdf/2012.02936>

- **Core idea:** Test whether two estimated clusters really differ in mean by conditioning on the event that the clustering algorithm produced them. This gives exact selective p-values under Gaussian models, and the 2026 work extends it to whole pipelines (outlier removal, feature selection, clustering) and to deep-latent clustering.
- **Strengths:** Rigorous fix for double dipping. Pipeline-level selective inference is an emerging (2026) direction.
- **Limitations:** Needs Gaussian noise with known or estimable covariance and tractable (polyhedral or piecewise) selection events. NMF with random restarts, top-gene truncation and overlap clustering gives non-tractable selection events. I found nothing for NMF or gene programs.
- **Relevance:** The cross-tumour clustering of programs into MPs is the step where 'is MP A really distinct from MP B?' (should they merge?) and 'is this cluster more than noise?' arise. Selective inference could test merges if programs were represented as Gaussian summaries, such as program-gene loadings with estimated covariance from bootstrap or thinning. In practice, count splitting or synthetic nulls are more feasible for MPs; selective inference is more of a theoretical companion.

### Synthetic null / in-silico negative control for double dipping (ClusterDE)

*Song D, Li JJ and colleagues. 'Synthetic control removes spurious discoveries from double dipping in single-cell and spatial transcriptomics data analyses', bioRxiv 10.1101/2023.07.21.550107 (2023). Software paper 'ClusterDE: A Statistical Software Package for Removing Double-Dipping Bias in Post-Clustering Differential Expression Analysis' (indexed in PubMed 41060830 / PMC13309996; journal not confirmed).* — ✅ verified online — <https://www.biorxiv.org/content/10.1101/2023.07.21.550107v2.full.pdf>

- **Core idea:** Simulate a synthetic null dataset (via scDesign3) that keeps gene means, variances and gene-gene correlations but contains a single cell type. Run the same clustering and DE pipeline on target and null, and use target-vs-null contrasts (a knockoff-style statistic) to control FDR of the discovered DE genes.
- **Strengths:** Pipeline-agnostic: the whole black-box pipeline is calibrated empirically. It controlled FDR where count splitting and the truncated-normal test did not, per the authors' benchmark.
- **Limitations:** Validity depends on how faithful the null simulator is. Preserving gene-gene correlation can keep real program structure in the 'null', which makes it conservative. Designed for clusters and DE, not programs.
- **Relevance:** The most directly adaptable FDR idea. For each tumour, generate nulls that keep nuisance structure (library size, QC, CNV-driven gene-block covariance from inferCNV, maybe cell cycle) but remove all other co-expression. Run the full MP pipeline on real and null cohorts. MPs whose recurrence or coherence beats null-derived MPs get an empirical FDR. This makes 'MPs not reflecting real biology' a calibrated quantity rather than a judgement call. I found no such application.

### Knockoff-style artificial variables to calibrate clustering (recall)

*DenAdel A, Ramseier M, Navia A, Shalek A, Raghavan S, Winter P, Amini A, Crawford L. 'Artificial variables help to avoid over-clustering in single-cell RNA sequencing', American Journal of Human Genetics, 2025 (bioRxiv title: 'A knockoff calibration method to avoid over-clustering in single-cell RNA-sequencing', 10.1101/2024.03.08.584180)* — ✅ verified online — <https://github.com/lcrawlab/recall>

- **Core idea:** Augment the data with synthetic 'fake genes' that carry no cell-type signal, then run preprocessing and clustering on the augmented matrix. A knockoff+-style data-dependent threshold on DE statistics (real vs artificial genes) decides whether clusters are real. Clusters without significant DE get merged and re-clustered.
- **Strengths:** Works with any clustering algorithm, scales to large data, and gives an FDR-flavoured threshold without exact distributional theory.
- **Limitations:** Artificial genes are drawn from marginals, so they lack real genes' correlation structure. Knockoff validity in this setting is heuristic, not model-X exact. Built for clusters, not programs.
- **Relevance:** Easy to transfer. Add knockoff genes to each tumour before NMF. A real program's top-gene list should exclude knockoffs, and knockoff loadings give a null distribution for gene-program weights. That supports a knockoff+ threshold for MP gene membership (FDR-controlled MP gene lists) and for deciding when K is too large: extra programs made of knockoff-like genes. Not applied to NMF programs.

### CSI-GEP: split-replication consensus NMF for GEPs at atlas scale

*Liu X et al. (St. Jude). 'CSI-GEP: A GPU-based unsupervised machine learning approach for recovering gene expression programs in atlas-scale single-cell RNA-seq data', Cell Genomics, Jan 2025 (DOI 10.1016/j.xgen.2024.100739)* — ✅ verified online — <https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00368-9>

- **Core idea:** Repeated GPU NMF on independent splits of the data across a sweep of K (starting at 2), about 100 replicates per split per K. Programs are aggregated by distance, and the K or programs that replicate across independent splits are kept, separating signal programs from noise programs.
- **Strengths:** The closest published 'replicate across independent data halves' criterion for GEP discovery. Scales to atlases on GPUs and recovers ground-truth GEPs in simulations, per the authors.
- **Limitations:** Splits by cells rather than by count thinning, so cell-level programs in small tumours lose power. Built for single datasets or atlases, not cross-tumour recurrence. No error-rate guarantee, and the replication criterion still passes stable technical programs.
- **Relevance:** A baseline and a design precedent. A new MP package should go further: split by count thinning (keeps every cell) and replicate across tumours and cohorts (held-out cohorts), with null calibration. Should be cited to show the difference.

### BayesPowerNMF: power posteriors for misspecification-robust NMF with automatic number of components

*Xue C, Miller JW, Carter SL, Huggins JH. 'Robust discovery of mutational signatures using power posteriors', PLoS Computational Biology 22(6):e1014372, 2026 (bioRxiv 10.1101/2024.10.23.619958). Code: TARPS-group/BayesPowerNMF* — ✅ verified online — <https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1014372>

- **Core idea:** Bayesian NMF with sparsity priors that switch off unneeded components. The likelihood is raised to a power zeta < 1 (a power or coarsened posterior), so small departures from the Poisson model do not create spurious extra components. zeta is chosen over a grid. In simulations it recovered more true signatures with lower error than SigProfilerExtractor and SignatureAnalyzer.
- **Strengths:** Tackles the specific failure mode the user describes: misspecification producing spurious programs. Gives principled automatic K and uncertainty quantification.
- **Limitations:** Shown on mutational catalogues (96 channels). MCMC cost in scRNA-seq dimensions is unclear. Choosing zeta needs calibration.
- **Relevance:** scRNA-seq NMF is heavily misspecified (overdispersion, ambient RNA, doublets, dropout). Standard Bayesian or ARD NMF then adds 'junk' programs to absorb misfit, which then recur across tumours because the misfit is systematic. A power-posterior multi-tumour NMF with zeta chosen by count-split held-out likelihood would be a new and well-justified combination. Not applied to scRNA-seq.

### Bayesian rank learning in NMF: compressive hyperpriors and Bayesian factor inclusion (bayesNMF)

*Zito A, Miller JW. 'Compressive Bayesian non-negative matrix factorization for mutational signatures analysis', arXiv 2404.10974, 2024. Landy JM, Basava N, Parmigiani G. 'bayesNMF: Fast Bayesian Poisson NMF with Automatically Learned Rank Applied to Mutational Signatures', Journal of Computational and Graphical Statistics, 2026 (arXiv 2502.18674)* — ✅ verified online — <https://www.tandfonline.com/doi/abs/10.1080/10618600.2026.2657487>

- **Core idea:** Compressive hyperpriors shrink unneeded factors toward zero while barely shrinking needed ones, with data-dependent hyperprior strength, and allow biologically informed priors on signatures. bayesNMF uses MH-within-Gibbs with Normal-likelihood proposals to avoid Poisson augmentation, plus Bayesian factor inclusion (BFI/SBFI, a BIC-based sparsity prior) to learn the rank while keeping posterior uncertainty.
- **Strengths:** Principled automatic rank with posterior uncertainty, replacing heuristic rank sweeps. Informed priors allow 'recovery-discovery' with a known catalogue.
- **Limitations:** Tested on mutational signatures (low-dimensional). Scaling to scRNA-seq would need variational or minibatch inference. Single-study models.
- **Relevance:** Replace the 'run NMF at ranks 4-9 and intersect' step with a posterior over per-tumour rank and programs. The informed-prior idea lets known MPs (Gavish 41) act as soft priors while novel MPs are learned. Not applied to scRNA-seq MPs.

### Bayesian nonparametric program models: IBP sparse factor models, HDP topic models, hierarchical Poisson factorization

*Knowles D, Ghahramani Z. 'Nonparametric Bayesian sparse factor models with application to gene expression modeling', Annals of Applied Statistics 2011 (arXiv 1011.6293; venue from memory). Gerber GK et al., GeneProgram (HDP; MIT CSAIL tech report 2007 found online). Liu J, Wade S, Bochkina N (author list from memory). 'Shared differential clustering across single-cell RNA sequencing datasets with the hierarchical Dirichlet process', arXiv 2212.02505. Levitin HM et al. 'De novo gene signature identification from single-cell RNA-seq with hierarchical Poisson factorization' (scHPF), Molecular Systems Biology 2019.* — ✅ verified online — <https://arxiv.org/pdf/1011.6293>

- **Core idea:** An IBP prior on a sparse loading matrix lets the number of factors grow with the data. An HDP gives a global pool of 'topics' (programs) that each group (tumour) reuses in different proportions, so the global program count is inferred. scHPF is count-native Poisson factorization with gamma hierarchies on cells and genes, though with fixed K.
- **Strengths:** Their generative story matches MPs: a global program pool, a tumour-specific subset, and cells mixing programs. Count likelihoods are available.
- **Limitations:** Mixing is slow and label switching is common. DP and HDP posteriors on the number of components are known to over-split under misspecification (results by Miller and Harrison; not re-verified here). None gives calibrated error rates. No modern scalable implementation for over a million malignant cells across about 1,000 tumours was found.
- **Relevance:** An 'HDP-Poisson factorization over tumours' is the cleanest probabilistic definition of a metaprogram, with tumour-private topics representing patient-specific programs. It would need power posteriors or count-split calibration (previous items) to avoid over-splitting into artifact programs. I found no tumour scRNA MP implementation.

**Gaps identified by this sweep**

- No count-native (Poisson or NB), non-negative, multi-tumour factor model yet combines (i) programs shared by any subset of tumours (Tetris/IBP-style ownership), (ii) explicit tumour-private and cohort-private programs (APAFA-style partitions), and (iii) scalable inference (SVI or amortised) for hundreds to about 1,000 tumours and over a million cells. Tetris, BMSFA and multi-study NMF exist only at bulk or mutational-signature scale, and GBCD/LIGER/MOFA+ lack subset-sharing with uncertainty.
- There is no error-rate or FDR definition for 'a discovered metaprogram'. Current pipelines (Gavish/3CA, GeneNMF) and GBCD give no null distribution for recurrence across tumours. Nuisance-preserving synthetic nulls (ClusterDE-style) or knockoff genes (recall-style) have not been applied to NMF programs or MP recurrence.
- MP gene lists are arbitrary (top 50 genes). Nobody has produced FDR-controlled, count-split-validated gene membership for programs. Count splitting or data thinning has not been used for per-tumour NMF rank selection or MP gene inference in tumour scRNA-seq, as far as I found.
- Data thinning validity under mixtures and misestimated overdispersion (Hivert et al. 2024) has been studied for clusters, not for continuous-factor NMF. A calibrated thinning scheme for factor models (conditional overdispersion from a preliminary low-rank fit, or robust statistics) is open.
- Selective inference for NMF-derived programs and for program-clustering steps (merge or split of MPs) does not exist. Selection events from NMF with restarts and top-k truncation are not tractable for current theory, even the 2026 pipeline selective-inference work.
- Identifiability: nobody has measured whether per-tumour NMF solutions on malignant scRNA-seq satisfy SSC or anchor conditions. A per-program identifiability diagnostic is missing, as is theory showing that pooling heterogeneous tumours restores identifiability of shared programs.
- Misspecification-induced spurious programs: power or coarsened posteriors (BayesPowerNMF) and compressive priors exist for mutational signatures but have not been carried over to scRNA-seq program discovery, where systematic misfit (ambient RNA, doublets, dissociation stress) recurs across tumours and creates recurrent artifact MPs.
- Known nuisance structure is not modelled as fixed or known factors in MP discovery: inferCNV-derived CNV gene blocks, QC metrics, cell cycle, technology and protocol. Multi-study NMF with covariates (Grabski 2025) offers a template but has not been applied to scRNA-seq.
- Stability is measured only across random restarts or rank, never across combined perturbations (cell bootstrap, count thinning, gene subsets, cohort hold-out, technology). No PCS-style stability-plus-null-calibration framework for MPs exists.
- There is no standard benchmark for MP methods: no realistic multi-tumour simulator with planted shared, subset-shared, tumour-private (CNV-driven) and technical programs, and no external replication protocol (held-out cohorts, cross-platform, spatial or perturbation validation). A Nature Methods paper would need to build one.
- No posterior uncertainty is reported on which tumours express each MP, or on MP hierarchy and nesting (e.g. EMT-I/II/III variants). Tree-structured or hierarchical sharing priors are untested for MPs.

**Notes:** VERIFICATION NOTES: The session's WebSearch budget (200 calls, apparently shared) ran out partway through. WebFetch could not reach arxiv.org, pmc/pubmed, nature.com, biorxiv, cell.com, elifesciences, genomebiology or github.io (proxy-blocked); github.com worked. 'verified_online=true' therefore means the paper's title, venue and key description appeared in live search results or a GitHub README in this session. Author lists marked 'from memory' or 'not confirmed' should be re-checked. MOFA+ (Argelaguet et al. 2020) was not re-verified, so that item is marked false. Two other works are cited only from memory and were not verified: the consensus-clustering pitfalls paper (Senbabaoglu et al., Sci Rep 2014) and stability selection (Meinshausen and Buhlmann, JRSSB 2010). GeneNMF's README asks users to cite Yerly L et al., 'Wounding triggers invasive progression in human basal cell carcinoma', bioRxiv 2024 (10.1101/2024.05.31.596823). In my searches I found no paper applying BMSFA, Tetris, APAFA, SUFA, multi-study NMF, data thinning, knockoffs or selective inference to tumour scRNA-seq metaprograms. The search for multi-study FA on tumour scRNA-seq returned GBCD and multiGroupVI, not MSFA-family methods.

STRATEGIC TAKE FOR A NATURE METHODS-LEVEL PACKAGE: The novelty is unlikely to be 'another factorisation'. It is more likely to be making 'metaprogram' a statistical object with a stated definition, an estimator, uncertainty and error control. A coherent design that the literature supports:
(1) Model. A hierarchical NB or Poisson NMF: cells nested in tumours nested in cohorts. A global program dictionary has per-(program, tumour) inclusion variables (IBP or spike-and-slab, as in Tetris or Grabski multi-study NMF) plus tumour-private and cohort-private slots (APAFA idea). Known nuisance factors enter as fixed covariates or offsets: inferCNV gene-block profiles, QC/ambient/doublet scores, cell cycle, technology. Fitted by SVI or amortised inference (as in MOFA+ or multi-study sparse VAEs) to reach 3CA scale.
(2) Robustness to misspecification. A power or coarsened posterior (BayesPowerNMF), with zeta and K chosen by held-out likelihood on count-split folds (Neufeld et al.).
(3) Identifiability. Volume or anchor regularisation (Fu/Huang/Sidiropoulos), plus a reported identifiability diagnostic per program.
(4) Error control. (a) FDR-controlled MP gene membership by testing gene-MP association on the held-out thinned fold. (b) An empirical FDR for each MP from nuisance-preserving synthetic null cohorts (ClusterDE/scDesign3 style) or knockoff genes (recall style), run through the identical pipeline. (c) Posterior credible sets for which tumours express each MP.
(5) Stability as a reported quantity across combined perturbations, not just seeds.
(6) Validation. Planted-program simulations including CNV-driven private programs and recurrent technical programs. Cross-cohort and cross-platform replication on 3CA. Prognostic and spatial follow-up.
Baselines to beat: the Gavish/3CA pipeline, GeneNMF, cNMF, GBCD (Nat Genet 2025, the main competitor), LIGER/iNMF, scINSIGHT, MOFA+ multi-group, CSI-GEP, and ideally Tetris or BMSFA on a subsample to show their scaling limits.

RISKS: Hivert et al. show that thinning can lose validity under mixtures with misestimated dispersion, so calibration must be shown empirically. Synthetic-null fidelity is the crux of any FDR claim. A Bayesian IBP model at 3CA scale is a real engineering challenge. GBCD already covers much of the 'shared vs patient-specific' story, so the main selling points should be error control, uncertainty, robustness to misspecification and modelling of nuisance and CNV.


## Confounders and ground-truth resources (23 items)

### Gavish et al. 3CA metaprogram pipeline (the baseline to beat)

*Gavish A, Tyler M, Greenwald AC, ... Tirosh I. 2023, Nature ('Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours'). The README of the code repo tiroshlab/3ca titles it 'The transcriptional hallmarks of intra-tumor heterogeneity across a thousand tumors'. Volume and DOI were not checked online.* — ✅ verified online — <https://github.com/tiroshlab/3ca/tree/main/ITH_hallmarks/Generating_MPs>

- **Core idea:** I read the code itself (robust_nmf_programs.R, Generate_Meta_Programs.R). Per-sample NMF: nmf(rank = 4:9, method = 'snmf/r', nrun = 10). Each program is reduced to its top 50 genes. A program is kept if it shares at least 35 of those 50 genes with a program at another rank in the same sample (intra_min = 35). Programs sharing more than 10 genes with an already-kept program in the same sample are dropped as redundant (intra_max = 10). A kept program must share at least 10 genes with a program from another sample (inter_min = 10). Clustering is greedy: the program with the most overlaps of at least 10 genes seeds a cluster, and programs keep joining while they share at least 10 genes with the cluster's gene set (Min_intersect_initial and Min_intersect_cluster = 10). Clusters need at least 5 members (Min_group_size = 5). The 50 MP genes are the cluster's genes ranked by how often they appear, with ties broken by NMF score. A 3CA_v2 folder holds an updated MP set, but the code that builds it is not in the repo; the new MPs are loaded precomputed from RDS files.
- **Strengths:** Simple and transparent. Validated by recurrence across about 1,000 tumors and many cancer types. It is the community reference and gene lists are available.
- **Limitations:** Every step is a hard threshold on overlaps between top-50 gene sets (35, 10, 10, 10, 5). This throws away gene weights and makes the result depend on rank, seed and threshold. The greedy cluster growth depends on order. There is no uncertainty on MP membership or on how many MPs there are. Nothing models confounders (CNA dosage, ambient RNA, doublets, dissociation stress, library size). 'Robustness' means recurrence, which an artifact that recurs across samples (e.g., stress) also satisfies.
- **Relevance:** This is the exact pipeline the new package has to beat. Its thresholds are what a benchmark should vary, to show the new method is stable where this one is not.

### Kinker et al. pan-cancer cell-line heterogeneity (recurrent heterogeneous programs)

*Kinker GS, et al. 2020, Nature Genetics 52:1208-1218, 'Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity', DOI 10.1038/s41588-020-00726-6* — ✅ verified online — <https://github.com/gabrielakinker/CCLE_heterogeneity>

- **Core idea:** scRNA-seq of about 200 cancer cell lines. NMF finds heterogeneity patterns that recur across cell lines and also appear in human tumors (recurrent heterogeneous programs, RHPs), and the paper compares the in vitro and in vivo patterns. The repo README warns that NMF gives stochastic outputs that can change with R and package versions.
- **Strengths:** Cell lines have no tissue microenvironment and no tissue dissociation, so a program that recurs in them is likely cell-intrinsic. The README's warning is first-hand evidence of the reproducibility problem the user describes.
- **Limitations:** Cultured cells lack hypoxia gradients, immune signals and stromal contact. Clonal and CNA structure within each line is still a confounder.
- **Relevance:** Gives an orthogonal 'cell-intrinsic' check: score each candidate MP on CCLE scRNA. The version-dependence warning motivates a seed-stable, deterministic design.

### Barkley et al. pan-cancer recurrent cell states (Yanai lab)

*Barkley D, et al. 2022, Nature Genetics ('Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment'). I found the code repo yanailab/PanCancer (its data links are under Dalia Barkley's account), but the README does not print the citation. Title and venue are from memory.* — ⚠️ not verified — <https://github.com/yanailab/PanCancer>

- **Core idea:** Finding.R converts per-sample NMF results into gene modules. A module is kept if its Jaccard similarity is above 0.05 with at least 3 other modules (sum(x > 0.05) >= 3), and this filter is applied iteratively. A gene co-occurrence graph is then clustered with cluster_infomap(nb.trials = 100), giving 16 programs, including Cycle, Stress, Interferon, Hypoxia, Oxphos and Metal. The repo also contains Visium.R for spatial analysis.
- **Strengths:** Uses graph community detection instead of greedy overlap. Programs were placed in space with Visium.
- **Limitations:** Still built on gene-set overlaps with arbitrary thresholds (Jaccard 0.05, at least 3 neighbors). 'Stress' and 'Metal' (metallothionein) programs are typical candidates for dissociation or ex vivo artifacts and need orthogonal checks.
- **Relevance:** A second baseline with a different clustering approach. Comparing Gavish-style and Barkley-style MPs on the same data is a cheap way to show the field's MPs depend on the method.

### GeneNMF (carmonalab), current R implementation of MP discovery

*Yerly L, et al. 2024, bioRxiv 10.1101/2024.05.31.596823 ('Wounding triggers invasive progression in human basal cell carcinoma'). This is the citation the package asks for. Version 0.6 is on GitHub.* — ✅ verified online — <https://github.com/carmonalab/GeneNMF>

- **Core idea:** multiNMF runs NMF across samples and several values of k. getMetaPrograms clusters the programs into MPs. From version 0.6 it uses cosine similarity of full gene-weight vectors instead of gene-set overlap. New options: a specificity.weight re-weighting, a weight.explained cumulative cutoff for which genes go into an MP, and min.confidence (the fraction of programs that contain a gene). It reports MP metrics such as silhouette and sample coverage.
- **Strengths:** Fast (uses RcppML) and integrates with Seurat. Moving from gene sets to weights already fixes part of the Gavish design.
- **Limitations:** Still heuristic clustering with a user-set number of MPs (nMP). No confounder model and no CNA awareness. No published ground-truth benchmark. Robustness comes from consensus across runs, without calibrated uncertainty.
- **Relevance:** The closest packaged competitor. A Nature Methods paper must benchmark against it as well as against the original Gavish code.

### cNMF, plus Harmony-corrected cNMF and starCAT/TCAT reference programs

*Kotliar D, et al. 2019, eLife (cNMF; elifesciences.org/articles/43803). Batch-corrected cNMF preprint: bioRxiv 10.1101/2024.05.03.592310 (2024). starCAT/TCAT: Nature Methods 2025, per the starCAT README.* — ✅ verified online — <https://github.com/immunogenomics/starCAT>

- **Core idea:** cNMF runs NMF many times, filters outlier components and clusters the rest into consensus programs. Rank is chosen by the trade-off between stability and error. The newer work corrects batch effects with Harmony before cNMF, then builds fixed reference programs (TCAT) and projects new datasets onto them (starCAT), so programs are reproducible across datasets. The starCAT README warns that whole-transcriptome references may not transfer reliably to probe-based assays such as 10x Flex.
- **Strengths:** The consensus-over-seeds idea is mature. The 'fixed reference plus projection' approach directly addresses reproducibility across cohorts. Kotliar 2019 benchmarked on simulated data with known identity and activity programs (the repo has an analyze_simulated_example_data notebook).
- **Limitations:** Designed for one dataset or a set of integrated datasets, not for tumor-by-tumor heterogeneity. Harmony correction can erase patient-specific biology. No CNA model. The platform-transfer caveat is stated by the authors themselves.
- **Relevance:** (i) A strong baseline (pooled cNMF versus per-tumor NMF). (ii) A design pattern: reference MPs plus projection, with measured reproducibility across datasets and platforms. (iii) Evidence that platform (3' vs Flex) is a real confounder.

### Spectra, supervised gene-program discovery

*Kunes RZ, Walle T, ... Pe'er D. Preprint bioRxiv 10.1101/2022.12.20.521311. Published in Nature Biotechnology 2024 according to my memory; the README cites only the preprint.* — ✅ verified online — <https://github.com/dpeerlab/spectra>

- **Core idea:** Factor analysis seeded with prior gene sets through a gene-graph loss weighted by lambda. It fits cell-type-specific factors alongside global ones and labels factors by overlap with the input gene sets (Szymkiewicz-Simpson coefficient).
- **Strengths:** Interpretable. Lets prior knowledge anchor factors, which reduces seed dependence. Separates factors specific to a cell type from shared ones.
- **Limitations:** Depends on how good and how biased the prior gene sets are. Not designed for per-patient malignant heterogeneity or CNA confounding.
- **Relevance:** A baseline and a design idea: putting known confounder signatures (stress, IFN, cell cycle, ambient profile) in as explicit 'nuisance factors' is a natural extension for MPs.

### MetaTiME and DIALOGUE, cross-dataset and multicellular program frameworks

*MetaTiME: Zhang Y, et al. 2023, Nature Communications (s41467-023-38333-8). DIALOGUE: Jerby-Arnon L & Regev A. 2022, Nature Biotechnology (s41587-022-01288-0).* — ✅ verified online — <https://github.com/yi-zhang/MetaTiME>

- **Core idea:** MetaTiME runs ICA on millions of tumor-microenvironment cells from hundreds of scRNA-seq datasets to learn reusable 'meta-components' that can be applied to new data. DIALOGUE finds multicellular programs, genes co-regulated across cell types, and can use spatial data to support them (repo: livnatje/DIALOGUE).
- **Strengths:** MetaTiME is pretrained and can be transferred. DIALOGUE uses coordination between cell types as evidence that a program is biological.
- **Limitations:** MetaTiME targets microenvironment cells, not malignant-cell MPs. DIALOGUE needs many samples and assumes the coupling is linear.
- **Relevance:** A real malignant MP (e.g., IFN response, hypoxia, EMT) should co-vary with matching states in microenvironment cells across patients. That cross-cell-type coupling is orthogonal to the malignant cells' own expression and can be part of a biological-validity score.

### CNA inference and clone tools: Numbat, CopyKAT, SCEVAN, infercna

*Numbat: Gao T, ... Kharchenko PV. 2022, Nature Biotechnology (s41587-022-01468-y). CopyKAT: Gao R, ... Navin NE. 2021, Nature Biotechnology (10.1038/s41587-020-00795-2). SCEVAN: De Falco A, et al. 2023, Nature Communications (10.1038/s41467-023-36790-9). infercna: Tirosh lab R package (GitHub).* — ✅ verified online — <https://github.com/kharchenkolab/numbat>

- **Core idea:** Infer copy-number profiles and subclones for each cell from scRNA. Numbat is haplotype- and allele-aware and also handles spatial data. CopyKAT separates aneuploid from diploid cells and builds subclones from 220 kb bins. SCEVAN classifies malignant cells and characterizes clones variationally (reported F1 0.90 across 106 samples).
- **Strengths:** Gives per-cell genetic covariates (clone label, dosage of each segment) from the same data used to find programs, with no extra assay needed.
- **Limitations:** Expression-based CNA calls get circular if used naively: a program concentrated in one chromosomal region can look like a CNA, and the reverse. Resolution is limited, and copy-neutral events are missed without allelic information.
- **Relevance:** This is the main tool for confounder (1). A CNA-aware MP method should (a) condition program discovery on clone and segment dosage, or regress them out, and (b) report each MP's association with clones, so genetic (subclone) programs are separated from non-genetic (plastic) ones. I could not find any packaged method that discovers programs jointly with CNA inference; this is a gap.

### clonealign and CalicoST, orthogonal genetic ground truth for clones

*clonealign: Campbell KR, et al. 2019, Genome Biology (10.1186/s13059-019-1645-z). CalicoST: Ma C, Balaban M, Liu J, Chen S, Ding L, Raphael BJ. 2024, bioRxiv per the README ('Inferring allele-specific copy number aberrations and tumor phylogeography from spatially resolved transcriptomics'). I believe it was later published in Nature Methods, but did not confirm that.* — ✅ verified online — <https://github.com/raphael-group/CalicoST>

- **Core idea:** clonealign assigns scRNA cells to clones defined independently by scDNA-seq. CalicoST infers allele-specific integer copy numbers and clone labels for each spatial spot, and reconstructs the phylogeography of clones.
- **Strengths:** Clone labels come from DNA (clonealign) or allelic data (CalicoST), so they are not circular with expression programs. CalicoST makes it possible to check whether a spatial MP domain is simply a clone domain.
- **Limitations:** Needs matched scDNA (clonealign) or enough spatial depth and SNP coverage (CalicoST). Few public tumors have scRNA, scDNA and spatial data together.
- **Relevance:** The strongest test of 'non-genetic program versus subclone artifact': an MP that stays after conditioning on DNA-defined clones, and is spread across clones in space, is a credible non-genetic state.

### Ambient RNA removal: CellBender, SoupX, DecontX

*CellBender: Fleming SJ, ... Babadi M. 2023, Nature Methods (10.1038/s41592-023-01943-7). SoupX: Young MD & Behjati S. 2020, GigaScience 9(12):giaa151. DecontX: Yang S, et al. 2020, Genome Biology (10.1186/s13059-020-1950-6).* — ✅ verified online — <https://github.com/broadinstitute/CellBender>

- **Core idea:** Estimate the ambient (soup) expression profile from empty droplets and remove contaminating counts from each cell. CellBender uses a generative model that also covers barcode swapping. SoupX reports contamination typically ranging from 2-50%.
- **Strengths:** Mature tools, and the ambient profile is measurable for each sample.
- **Limitations:** Correction is incomplete and can remove real signal. Ambient fractions differ by sample and are partly confounded with patient.
- **Relevance:** In tumors, ambient RNA from microenvironment cells (immunoglobulins, collagens, hemoglobin, myeloid genes) can show up as low-level 'programs' inside malignant cells that look like EMT, immune or stromal programs. A new method could use each sample's measured ambient profile as a fixed nuisance factor. Semi-synthetic tests can inject ambient profiles at known fractions.

### Doublet detection: scDblFinder and the Xi & Li benchmark

*scDblFinder: Germain PL, et al. 2021, F1000Research (10.12688/f1000research.73600.2). Benchmark: Xi NM & Li JJ, Cell Systems (10.1016/j.cels.2020.11.008; addendum arXiv:2101.08860).* — ✅ verified online — <https://github.com/plger/scDblFinder>

- **Core idea:** Detect heterotypic doublets by simulating artificial doublets and classifying cells. In the independent Xi & Li benchmark, scDblFinder had the highest mean AUPRC and AUROC across 16 datasets.
- **Strengths:** Validated, with an independent benchmark.
- **Limitations:** Homotypic doublets (malignant with malignant) are hard to detect. Malignant-microenvironment doublets may be partly removed or flagged.
- **Relevance:** A malignant-fibroblast or malignant-macrophage doublet produces a mixed 'program' (e.g., partial-EMT-like or immune-like) and can pass the recurrence filter. The benchmark should inject doublets at known rates and score false MPs.

### Dissociation, ex vivo stress and protocol artifacts (classic evidence)

*Not checked online in this session (from memory, high confidence that these exist): van den Brink SC, et al. 2017, Nature Methods ('Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations'). O'Flanagan CH, et al. 2019, Genome Biology (cold-active protease reduces collagenase-associated stress responses in tumors). Denisenko E, et al. 2020, Genome Biology (biases from dissociation and storage in scRNA and snRNA workflows). Slyper M, et al. 2020, Nature Medicine (scRNA/snRNA toolbox for fresh and frozen human tumors).* — ⚠️ not verified

- **Core idea:** Enzymatic dissociation at 37°C triggers immediate-early and heat-shock genes (FOS, JUN, EGR1, HSPA1A/B, etc.) within minutes. This signature differs by cell type and protocol and can be reduced by cold proteases, fixation before dissociation, or single-nucleus profiling. Barkley's verified code lists 'Stress' and 'Metal' programs among its 16 states, the typical places where such artifacts land.
- **Strengths:** Explains mechanistically why 'stress' MPs recur in every cancer type. Gives protocol contrasts (cold versus warm, cells versus nuclei, fresh versus fixed) that can act as natural ground truth.
- **Limitations:** Real in vivo stress biology (hypoxia, UPR, heat shock near necrosis) overlaps the artifact signature, so a gene-list filter alone cannot separate them.
- **Relevance:** The key test is whether an MP persists in fixed (10x Flex), single-nucleus or in situ (Xenium, MERFISH, CosMx) data from the same tumor. Artifact programs should drop sharply there, and real stress programs should stay spatially organized (e.g., near necrosis).

### Matched multi-platform tumor datasets (Xenium, Visium and Chromium Flex on the same FFPE breast cancer; Visium HD, Flex and Xenium on colorectal cancer)

*Janesick A, Shelansky R, Gottscho A, Wagner F, Rouault M, et al. 2023 (published 20 Oct 2023; Nature Communications per the companion repo), 'High resolution mapping of the breast cancer tumor microenvironment using integrated single cell, spatial and in situ analysis of FFPE tissue'. Oliveira MF, Romero JP, Chung M, et al. 2025, Nature Genetics, 'High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer' (10.1038/s41588-025-02193-3).* — ✅ verified online — <https://github.com/10XGenomics/HumanColonCancer_VisiumHD>

- **Core idea:** The same tumor blocks were profiled with dissociation-free fixed scRNA (Chromium Flex), spot or high-definition spatial data (Visium; Visium HD at 2 µm bins with nucleus segmentation) and imaging in situ (Xenium), with registered images.
- **Strengths:** Public, matched across platforms, and partly free of dissociation. They let you test whether an MP found in one protocol reappears in another, both in the cells and in space.
- **Limitations:** Vendor-generated, with few patients and two cancer types. Probe panels (Flex, Xenium) cover a subset of genes, which limits gene-level comparison of programs (see the starCAT caveat).
- **Relevance:** Concrete data for a cross-protocol invariance benchmark. A method that recovers the same malignant MPs (by gene-weight and usage concordance) across Flex, Visium HD and Xenium, while artifact MPs drop out, is persuasive evidence of real biology.

### Greenwald et al. integrative spatial analysis of glioblastoma (spatial MP coherence)

*Greenwald AC, Galili Darnell N, Hoefflin R, et al. ('Integrative spatial analysis reveals a multi-layered organization of glioblastoma'). Title and authors verified via the tiroshlab/Spatial_Glioma repo; venue Cell 2024 from memory.* — ✅ verified online — <https://github.com/tiroshlab/Spatial_Glioma>

- **Core idea:** Visium and CODEX analysis of glioblastoma. Builds metaprograms per spot (extending Kinker 2020 and Gavish 2023), then measures spatial coherence, associations between programs, interactions, and CNAs in space, with protein imaging as support.
- **Strengths:** Shows MPs are spatially organized (e.g., hypoxia and mesenchymal-like near necrosis), and introduces spatial coherence as a quantitative property of a program.
- **Limitations:** Visium spots mix cells, so malignant programs are confounded by cell-type composition. One tumor type.
- **Relevance:** Gives a metric (spatial coherence or autocorrelation of program usage) and expected spatial organizations that a new method's MPs should reproduce, ideally on single-cell-resolution platforms.

### Genome-scale Perturb-seq (Replogle 2022) and transcriptome-wide effect calibration (TRADE)

*Replogle JM, et al. 2022, Cell ('Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq'). Title and venue from memory; the existence of the Replogle 2022 K562 and RPE1 essential-gene datasets was confirmed via the GEARS repo. TRADE: bioRxiv 10.1101/2024.07.03.601903 (2024; first author Nadig A per my memory, verified only via the TRADEtools repo).* — ✅ verified online — <https://github.com/snap-stanford/GEARS>

- **Core idea:** Pooled CRISPRi Perturb-seq across thousands of genes in K562 (genome-wide and essential) and RPE1 (essential). TRADE estimates how DE effects of each perturbation are distributed across the transcriptome, separating broad 'transcriptome-wide impact' from specific effects.
- **Strengths:** Causal interventions: a program whose usage moves consistently when known regulators are knocked down (e.g., MYC targets, mTOR or UPR, cell-cycle checkpoints, IFN pathway genes) is a real co-regulated unit. Non-targeting guides give an empirical null.
- **Limitations:** Few cell lines (K562 and RPE1 are atypical for solid tumors). Mostly essential genes. Programs that exist only in vivo (hypoxia gradients, invasion, immune-driven states) are absent.
- **Relevance:** Gives ground truth for 'is this a coherent, causally regulated program'. Run each MP method on control cells or all cells, then test how many programs respond specifically to perturbations (calibrated with non-targeting guides), and whether the responding regulators make biological sense.

### Tahoe-100M drug perturbation atlas across cancer cell lines

*Tahoe Therapeutics with Arc Institute (Vevo/Tahoe team). 2025, bioRxiv 10.1101/2025.02.20.639398, 'Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent Gene Function and Cellular Modeling'.* — ✅ verified online — <https://github.com/ArcInstitute/arc-virtual-cell-atlas/tree/main/tahoe-100M>

- **Core idea:** 100,648,790 cells in 1,344 samples: cancer cell lines (Cellosaurus IDs) treated with many drugs at several concentrations, profiled by multiplexed scRNA-seq. From memory, about 50 cell lines and about 1,100 drug treatments; the README I read gives neither count.
- **Strengths:** The largest perturbation resource with many cancer cell lines. Drugs with known mechanisms (CDK4/6, MEK, mTOR, HDAC, proteasome, IFN-related) should induce or suppress specific MPs (cell cycle, MAPK, stress, UPR, EMT-like, IFN) across many genetic backgrounds.
- **Limitations:** In vitro only. Effects are dominated by drug response, not native heterogeneity. Pooled multiplexing can add ambient and doublet artifacts of its own.
- **Relevance:** Test whether MPs from untreated or DMSO cells (i) recur across cell lines and (ii) respond as expected to mechanism-matched drugs. This is a cancer-specific, causal check beyond Replogle.

### scPerturb harmonized perturbation compendium

*Peidli S, Green TD, et al. 2024, Nature Methods ('scPerturb: harmonized single-cell perturbation data'; README title 'A resource and a python/R tool for single-cell perturbation data'; s41592-023-02144-y).* — ✅ verified online — <https://github.com/sanderlab/scPerturb>

- **Core idea:** Harmonized collection of single-cell perturbation datasets (RNA, ATAC, CITE-seq) with E-distance and E-test statistics for perturbation effects.
- **Strengths:** One uniform format for pulling cancer-relevant perturbation datasets (e.g., CRISPR in melanoma or cancer lines and compound screens), with a statistical framework for effect size.
- **Limitations:** Datasets differ in quality and power. Mostly cell lines.
- **Relevance:** Makes it practical to run the perturbation-response metric across many datasets instead of only Replogle.

### PerturbNMF (Engreitz lab), a modern evaluation suite for gene programs

*Engreitz Lab, 2026, software (MIT license; Zenodo concept DOI 10.5281/zenodo.21878267). A fork of Kotliar's cNMF. No paper is cited in the README. A companion repo analyzes TeloHAEC genome-wide Perturb-seq with cNMF (K = 50 and 150).* — ✅ verified online — <https://github.com/EngreitzLab/PerturbNMF>

- **Core idea:** Three-stage pipeline: cNMF inference on CPU or GPU, then evaluation with nine metrics, then interpretation with LLM-based annotation. The metrics are categorical association, perturbation sensitivity (U-test, CRT or matched-cell DE, calibrated against non-targeting guide nulls), motif enrichment, GWAS/OpenTargets trait enrichment, GO and MSigDB enrichment, explained variance, reconstruction error, and stability.
- **Strengths:** The most complete program-evaluation metric set I found, and it uses calibrated perturbation nulls.
- **Limitations:** Not tumor-specific. No CNA, ambient or dissociation confounder metrics. Only cNMF-based programs. Not peer-reviewed.
- **Relevance:** A template for the benchmark's 'biological validity' panel. A new MP package could adopt or extend these metrics and add tumor-specific ones: clone-independence, cross-protocol invariance, spatial coherence, and recurrence in cell lines.

### Lineage-traced tumors (KP-Tracer) and Cassiopeia, ground truth for heritable versus plastic programs

*KP-Tracer: Yang D, Jones MG, et al. ('Lineage recording reveals the phylodynamics, plasticity and paths of tumor evolution'). bioRxiv 2021 per the repo; published in Cell 2022 from memory. Cassiopeia: Jones MG, et al. 2020, Genome Biology (10.1186/s13059-020-02000-8).* — ✅ verified online — <https://github.com/mattjones315/KPTracer-release>

- **Core idea:** CRISPR-based evolving lineage recorders in KP mouse lung tumors, with scRNA-seq readout. Phylogenies are reconstructed with Cassiopeia, which allows quantifying how plastic or heritable each transcriptional state is.
- **Strengths:** Lineage is measured independently of expression and genetics, so heritability of a program can be computed (e.g., phylogenetic autocorrelation with Hotspot).
- **Limitations:** Mouse model only. Lineage recording is not possible in human tumors (approximations: mtDNA variants, CNA clones).
- **Relevance:** Lets the method separate an MP that is a heritable lineage or subclone state from one that is a transient state. That is a biological annotation Gavish-style MPs cannot give, and a strong validation for a 'program type' classifier.

### Hotspot, testing a program's coherence on spatial or lineage graphs

*DeTomaso D & Yosef N. 2021, Cell Systems (S2405471221001149; preprint bioRxiv 10.1101/2020.02.06.937805).* — ✅ verified online — <https://github.com/YosefLab/Hotspot>

- **Core idea:** Local autocorrelation of genes and modules on a cell-similarity graph that can be transcriptional, spatial or a lineage tree. Clusters locally correlated genes into modules and scores them per cell.
- **Strengths:** One statistic for 'is this program coherent in space' and 'is it heritable along lineage', with p-values.
- **Limitations:** Needs a good graph. Module detection itself still depends on parameters.
- **Relevance:** Use it as an evaluator, not a competitor. For each MP, compute spatial autocorrelation on Xenium or Visium HD and phylogenetic autocorrelation on lineage data. Real MPs should score high on whichever axis fits the biology.

### SCENIC+, multiome regulatory support for programs

*Bravo González-Blas C, De Winter S, et al. 2023, Nature Methods (10.1038/s41592-023-01938-4).* — ✅ verified online — <https://github.com/aertslab/scenicplus>

- **Core idea:** Infers enhancer-driven regulons (eRegulons) from scRNA plus scATAC or 10x Multiome by linking TF motifs, accessible enhancers and target genes.
- **Strengths:** Chromatin evidence is independent of RNA artifacts such as ambient RNA and dissociation immediate-early genes.
- **Limitations:** Multiome tumor data are rarer. The inferred regulons are themselves predictions, not ground truth.
- **Relevance:** An orthogonal 'regulatory coherence' check. A real MP's genes should be enriched for shared TF-enhancer regulation, with matching accessibility changes in the same cells. Dissociation IEG artifacts should lack that chromatin support because they arise within minutes.

### Simulators that can encode gene programs: Splatter/splatPop, SymSim, SERGIO, scDesign3, scMultiSim

*Splatter: Zappia L, Phipson B, Oshlack A. 2017, Genome Biology (10.1186/s13059-017-1305-0). splatPop: Azodi CB, et al. 2021, Genome Biology (10.1186/s13059-021-02546-1). SymSim: Zhang X, Xu C, Yosef N. 2019, Nature Communications (10.1038/s41467-019-10500-w). SERGIO: Dibaeinia P & Sinha S (Cell Systems 2020 from memory; the README cites no paper). scDesign3: Song D, et al. 2024, Nature Biotechnology 42:247-252 (10.1038/s41587-023-01772-1). scMultiSim: Li H, Zhang Z, Squires M, Chen X, Zhang X (bioRxiv 2023 per README; I believe since published in Nature Methods, not confirmed).* — ✅ verified online — <https://github.com/SONGDONGYUAN1994/scDesign3>

- **Core idea:** Splatter is a parametric NB simulator; Kotliar 2019 extended it to add activity programs, and splatPop adds a population or eQTL structure. SymSim models extrinsic variability factors and co-expression modules with technical noise and batch. SERGIO and scMultiSim are driven by GRNs with technical noise, and scMultiSim adds cell-cell interactions, spatial layout and batch. scDesign3 fits real reference data with transparent, editable parameters.
- **Strengths:** Exact ground truth for program genes and per-cell usage. scDesign3 and splatPop can mimic real data and patient structure.
- **Limitations:** None models tumor-specific confounders together: segment-level CNA dosage with subclones, ambient RNA from the microenvironment, malignant-microenvironment doublets, dissociation IEG bursts that vary by cell, and differences between platforms. Purely synthetic benchmarks are weak evidence for a Nature Methods claim about 'real biology'.
- **Relevance:** Best used for semi-synthetic spike-ins: start from real malignant cells (or scDesign3 fits of them), add known programs at controlled usage, and add confounders measured from the same tumors (inferred CNA segments, empty-droplet ambient profiles, doublets with real microenvironment cells, a stress signature). Then measure recovery and the false-program rate.

### Benchmark precedents and their metrics: BEELINE, BioBombe, momix, scIB (plus Saelens module-detection benchmark)

*BEELINE: Pratapa A, et al. 2020, Nature Methods (10.1038/s41592-019-0690-6). BioBombe: Way GP, et al. 2020, Genome Biology (10.1186/s13059-020-02021-3). momix: Cantini L, et al. 2021, Nature Communications (10.1038/s41467-020-20430-7). scIB: Luecken MD, et al. 2022, Nature Methods (10.1038/s41592-021-01336-8). Saelens W, Cannoodt R, Saeys Y. 2018, Nature Communications (module-detection evaluation): repo saeyslab/moduledetection-evaluation exists, citation from memory.* — ✅ verified online — <https://github.com/Murali-group/Beeline>

- **Core idea:** BEELINE scores against synthetic (BoolODE) and curated ground-truth networks using AUPRC, AUROC, early precision ratio, stability (Jaccard of top-k) and runtime. BioBombe compares PCA, ICA, NMF, DAE and VAE at 28 latent dimensions on TCGA, GTEx and TARGET, and finds that different k values capture complementary biology, so a single 'optimal' k is misleading. momix scores factor methods by clustering on simulated data, association with survival and clinical variables in TCGA, and pathway enrichment of metagenes. scIB defines biological-conservation metrics (cell-type ASW, NMI/ARI, isolated labels, cell-cycle conservation, trajectory conservation) versus batch-removal metrics (iLISI, kBET, graph connectivity, PCR). Saelens 2018 used known regulatory modules as ground truth for module detection.
- **Strengths:** Established metric vocabularies and benchmark designs that Nature Methods reviewers recognize.
- **Limitations:** None targets malignant MPs in tumor scRNA. Despite trying, I did not find a dedicated 2022-2026 benchmark of gene-program methods on tumor scRNA (but literature search was severely limited in this session; see notes).
- **Relevance:** (i) BioBombe directly supports moving from a fixed rank to a multi-resolution or rank-free model. (ii) scIB's split between conservation and removal maps onto 'remove patient, CNA and artifact signal while keeping biology'. (iii) BEELINE-style AUPRC plus stability is the template for scoring recovery of known programs.

**Gaps identified by this sweep**

- I found no dedicated community benchmark for gene-program or metaprogram discovery in tumor scRNA-seq (no shared ground truth, datasets or metrics), unlike GRN inference (BEELINE) or integration (scIB). Caveat: web and literature search were unavailable in this session, so a 2025-2026 preprint could exist and must be checked before claiming novelty.
- No simulator jointly models tumor-specific confounders: segment-level CNA dosage with subclones, ambient RNA from microenvironment cells, malignant-microenvironment doublets, dissociation immediate-early bursts that vary by cell, library-size and dropout differences, and platform differences (3' vs 5' vs Flex vs snRNA vs probe panels). A semi-synthetic 'tumor confounder simulator' would be a contribution in its own right.
- No verified method discovers programs jointly with CNA or clone inference, or conditions program discovery on clone and segment dosage. Current practice is post hoc (inferCNV/Numbat, then manual inspection) or nothing.
- MP 'robustness' has no statistical definition. Gavish thresholds (top-50 overlaps of 35/10/10, cluster size of at least 5) and Barkley thresholds (Jaccard > 0.05 with at least 3 neighbors) are ad hoc. There is no uncertainty quantification for gene membership, MP number or per-cell usage, and no explicit stability over seeds, ranks or subsamples (the Kinker README itself warns that NMF output changes with R and package versions).
- Patient effect versus biology is unresolved. Integration (e.g., Harmony before cNMF) can erase patient-specific but real biology. Recurrence filters keep artifacts that recur (stress, ambient). A hierarchical model of shared versus patient-specific programs, with explicit nuisance factors, is missing.
- Causal ground truth comes mostly from cell lines (K562 and RPE1 in Replogle; cancer lines in Tahoe-100M). In vivo-specific programs (hypoxia gradients, invasion front or partial EMT, immune-induced IFN, stromal-contact states) have little perturbational ground truth. In vivo or organoid Perturb-seq with spatial readout in tumors is still scarce.
- Matched same-tumor, multi-protocol datasets (fresh 3' vs 5' vs Flex-fixed vs snRNA vs Xenium/MERFISH/CosMx vs Visium HD) exist only as a few vendor datasets (Janesick breast 2023, Oliveira CRC 2025). A purpose-built, multi-cancer 'protocol-invariance' reference dataset is missing and would make the strongest ground truth for separating artifact from biology.
- Lineage ground truth for heritable versus plastic programs exists only in mouse models (KP-Tracer). Human tumors need proxies (mtDNA variants, DNA-defined clones), and no MP method reports program heritability.
- Probe-based assays (Flex, Xenium or CosMx panels) cover only part of the transcriptome. Methods to transfer or validate whole-transcriptome MPs onto targeted panels are immature (the starCAT authors flag this).
- Separating ex vivo dissociation stress from real in vivo stress (hypoxia, UPR, heat shock near necrosis) has no validated computational solution. Gene-list filters remove both.
- Annotating and naming MPs is still manual (enrichment plus expert curation). LLM-based annotation (as in PerturbNMF) is not validated against ground truth.
- Benchmark metrics are not agreed. Recovery (AUPRC of program genes, correlation of usage), stability (matched cosine or Jaccard across seeds, ranks and bootstraps), false-program rate against injected confounders, perturbation responsiveness, spatial coherence, cross-protocol reproducibility, clone-independence and recurrence in cell lines have never been combined into one scoring framework for MPs.

**Notes:** IMPORTANT LIMITATIONS OF THIS SEARCH: Web search was unavailable from the start of this subagent: the session's WebSearch budget (200 of 200) was already spent. The egress proxy also blocked nature.com, biorxiv.org, pubmed and ncbi, europepmc, crossref, openalex, arxiv, elifesciences, bioconductor, CRAN, huggingface and readthedocs. Everything above was verified by reading GitHub repositories (README files and source code via github.com and raw.githubusercontent.com) and PyPI. So verified_online=true means I confirmed the work or tool exists and checked the citation details printed in its repo. Where a journal or year is marked 'from memory', it was not confirmed. Coverage of 2025-2026 preprints is therefore incomplete. Before writing a novelty claim, re-run with literature search to look for: (a) any 2024-2026 benchmark of gene-program discovery; (b) CNV-aware or clone-aware program methods; (c) newer genome-scale Perturb-seq atlases (e.g., HepG2 and Jurkat, X-Atlas/Orion, primary T cells); (d) newer spatial pan-cancer MP validation work (Xenium, CosMx, MERFISH); (e) whether scMultiSim and CalicoST have since been published in Nature Methods (I believe so, not confirmed).

Verified facts worth quoting: the Gavish code uses nmf(rank = 4:9, method = 'snmf/r', nrun = 10), top-50 genes, intra_min 35, intra_max 10, inter_min 10, Min_intersect 10, Min_group_size 5 and greedy clustering. Barkley uses Jaccard > 0.05 with at least 3 neighbors, iterated, then infomap clustering, giving 16 programs including Stress and Metal. GeneNMF v0.6 moved to cosine similarity of full weight vectors. starCAT warns that whole-transcriptome references may not transfer to 10x Flex. The Kinker README warns that NMF output varies by R and package version.

WHY TUMOR PROGRAMS ARE CONFOUNDED (synthesis): (1) CNA dosage. Genes in amplified or deleted segments co-vary within subclones, so NMF finds 'programs' that are really chromosomal blocks or clone identity, and these pass the within-tumor recurrence filter. (2) Patient or clone identity. Per-tumor NMF absorbs subclone differences, and the cross-tumor recurrence filter only partly removes them. (3) Dissociation and ex vivo stress (immediate-early and heat-shock genes, metallothioneins): technically driven, yet they recur across tumors and so pass recurrence filters. (4) Ambient RNA and doublets: microenvironment-derived signal creates EMT-, immune- or stroma-like 'malignant' programs. (5) Cell cycle: real but dominant; it can take over a rank and split other programs. (6) Library size, complexity and dropout: a low-quality or complexity axis often appears as a mito, ribosomal or 'housekeeping' program. (7) Platform: 3' vs 5' vs Flex probes vs Smart-seq2 vs snRNA change gene detectability and remove dissociation signals in fixed or nuclear protocols. Because the Gavish recipe defines 'robust' as 'recurs', artifacts that recur systematically are kept by design. This is the conceptual weakness a new method should target.

MOST CONVINCING GROUND-TRUTH STRATEGY (recommended evidence ladder for a Nature Methods paper; no single source is enough):
(A) Semi-synthetic spike-ins on real tumor data. Inject known programs at controlled usage plus measured confounders: CNA segments from Numbat or CopyKAT, empty-droplet ambient profiles from the same samples, doublets with real microenvironment cells, a dissociation signature at varying fractions, library-size shifts. Report program-gene AUPRC, usage correlation, stability over seeds, ranks and subsamples, and a FALSE-MP rate (MPs aligned to the injected confounders). This is where Gavish-style NMF should visibly fail.
(B) Causal perturbation ground truth in cancer cells: Tahoe-100M (mechanism-matched drug response across many cancer lines), Replogle K562 and RPE1, scPerturb datasets. Report the fraction of MPs specifically responsive to at least one perturbation, calibrated against non-targeting or DMSO nulls as in PerturbNMF, and the biological coherence of the responding regulators.
(C) Cross-protocol invariance on the same tumors. Fresh 3' vs Flex (fixed before dissociation) vs snRNA vs Xenium or Visium HD (Janesick breast, Oliveira CRC; ideally a new purpose-built dataset). Real MPs replicate across protocols; dissociation and ambient MPs vanish in fixed and in situ data. This is the most direct answer to 'is it real biology'.
(D) Spatial coherence in situ: spatial autocorrelation (Hotspot or Moran's I) and pre-registered niche associations (hypoxia near necrosis, EMT at the invasive edge, IFN near T cells), plus checking with CalicoST that the spatial MP domains are not just clone domains.
(E) Genetic orthogonality and heritability: MP usage should be independent of DNA-defined clones (clonealign or Numbat), or the method should flag the MP as clone-linked. Heritability should be measured on lineage-traced tumors (KP-Tracer).
(F) Cell-intrinsic recurrence in cancer cell lines (Kinker CCLE scRNA) and replication on held-out cohorts and cancer types (3CA), and optionally orthogonal regulatory support (SCENIC+ on multiome).
Strategies A, B and C carry the most weight. A gives exact truth for the confounder problem the user cares about, B gives causal truth, and C gives artifact-versus-biology truth on real tumors. D and E add tumor-specific biological credibility, and F adds generality. A combined scorecard (recovery, stability, false-MP rate, perturbation responsiveness, cross-protocol reproducibility, spatial coherence, clone-independence) released as an open benchmark would itself address the field's biggest gap and strengthen the Nature Methods case.

PACKAGE DESIGN IMPLICATIONS (for the orchestrator): use a hierarchical, weight-based (not top-50 gene-set) model of shared versus patient-specific programs with explicit nuisance factors (CNA segment dosage or clone, per-sample ambient profile, doublet score, stress, library size and complexity, cell cycle as a known factor). Use rank-free or multi-resolution inference (supported by BioBombe), report uncertainty for membership and MP number via bootstrap or Bayesian methods, classify each MP's type (genetic or clone-linked, plastic, protocol-sensitive or artifact), and use reference projection for new datasets (starCAT pattern).


## Publication bar and scoop risk (22 items)

### Gavish et al. 2023 'Hallmarks of transcriptional ITH' (the current standard), plus the Tirosh lab's unpublished 3CA_v2 metaprogram update

*Gavish A, Tyler M, Greenwald AC, ... Suvà ML, Puram SV, Tirosh I. Nature 618:598–606 (2023). doi:10.1038/s41586-023-06130-4. 3CA_v2 code is in github.com/tiroshlab/3ca (3CA_v2/MP_analysis_new_batch.R, mp_new.R); I found no paper attached to it.* — ✅ verified online — <https://github.com/tiroshlab/3ca>

- **Core idea:** 77 studies, 1,163 tumour samples, 24 tumour types. NMF is run on each sample's malignant cells, programs are defined by their top-50 genes, programs that recur across ranks and across samples are kept, and greedy overlap clustering turns them into consensus MPs (I recall 41 MPs; I did not re-check that number). The 3CA_v2 script runs the same greedy algorithm with 64 new datasets added and uses explicit hard thresholds: top 7,000 genes by mean expression; centred data with negatives set to 0; ranks 4–9 with 10 runs each; within-sample overlap of at least 35 genes; a founder program needs at least 10 overlapping partners; a program joins an MP at 10 or more shared genes; an MP needs at least 5 programs; single-study clusters are removed.
- **Strengths:** Huge, well-curated resource (3CA). The MPs are widely reused; an igordot/sns fgsea script lists '3CA' as an MSigDB category. The procedure is simple and biologically interpretable.
- **Limitations:** Every step uses hard cutoffs (top-50 truncation, gene-overlap counts, minimum group size). There is no uncertainty, no significance test and no null model. Greedy clustering depends on order and founder choice. CNV, clone, protocol and dissociation effects are not modelled explicitly; a 'single-nucleus' MP in a 2026 MP list I found hints that protocol-driven MPs exist.
- **Relevance:** This is the baseline the user must beat. Expect the Tirosh lab to publish an expanded MP catalogue from 3CA_v2, so benchmarks should cover both the 2023 and v2 MPs. The 3CA_v2 thresholds are exactly the parameters a robustness or sensitivity analysis should sweep.

### cNMF (consensus NMF) and its 2024–2026 updates

*Kotliar D et al. eLife 2019 (8:e43803). Recent repo updates v1.5–1.7 (github.com/dylkot/cNMF).* — ✅ verified online — <https://github.com/dylkot/cNMF>

- **Core idea:** Many NMF replicates are clustered into consensus programs, and K is chosen by a stability vs error trade-off. Recent releases add a Harmony-based Preprocess class for batch-correcting before cNMF (described in the 2024 TCAT preprint), building starCAT references in the consensus step, sparse-matrix and numpy-2 support, and h5ad input. A GPU 'torch-cNMF' fork with minibatching now exists (used by PerturbNMF).
- **Strengths:** De facto standard (about 470 GitHub stars), pip install, scanpy and AnnData compatible, actively maintained.
- **Limitations:** Works per dataset. Consensus covers random seeds, not datasets or patients. K selection is still heuristic, and there is no notion of programs recurring across tumours.
- **Relevance:** Reviewers will expect cNMF as the per-tumour engine baseline, or as a competitor run on pooled data. Its consensus-over-replicates idea is what the user's robustness framework must generalise: consensus over seeds × ranks × tumours × gene-set choices.

### starCAT / TCAT: a reproducible reference of gene expression programs plus a projection tool

*Kotliar D, Curtis M, Agnew R, Weinand K, Nathan A, Baglaenko Y, Slowikowski K, Zhao Y, Sabeti PC, Rao DA, Raychaudhuri S. 'Reproducible single-cell annotation of programs underlying T cell subsets, activation states and functions.' Nature Methods 22(9), 2025 (s41592-025-02793-1). Preprint: bioRxiv 10.1101/2024.05.03.592310.* — ✅ verified online — <https://github.com/immunogenomics/starCAT>

- **Core idea:** Batch-corrected cNMF is run on several T-cell datasets, and programs that recur across datasets are merged into a fixed reference catalogue (TCAT.V1 on Zenodo). starCAT then projects any new dataset onto the fixed programs using NNLS, so program usage is comparable across studies. From memory, not re-checked: about 46 GEPs from roughly 1.7M T cells across many tissues and diseases, validated with CITE-seq surface proteins and in-vitro activation, and applied to immunotherapy response.
- **Strengths:** Ships Python (pip starcatpy) and R tutorials, a CLI, Zenodo references, and custom-reference building. Third-party R ports already exist (CellProgramMapper), which shows uptake.
- **Limitations:** T cells only. The catalogue is fixed, so new biology needs a rebuild. Projection reports no per-cell uncertainty, and cross-dataset recurrence is decided heuristically.
- **Relevance:** The closest blueprint for a 2025 Nature Methods 'program' paper: discovery across many datasets → reproducible reference → projection tool → biological validation. A malignant-cell analogue ('MP reference + projection') is an obvious next step, which raises the novelty bar: 'cNMF across tumours plus a reference' alone may look derivative.

### Spectra: supervised (prior-knowledge-guided) discovery of interpretable gene programs

*Kunes RZ, Walle T, Land M, Nawy T, Pe'er D. Nature Biotechnology 42:1084–1095 (2024). doi:10.1038/s41587-023-01940-3 (preprint bioRxiv 10.1101/2022.12.20.521311).* — ✅ verified online — <https://github.com/dpeerlab/spectra>

- **Core idea:** Factor analysis that starts from cell-type-specific and global gene sets (a gene–gene knowledge graph) and lets factors depart from those priors to find new programs, with cell-type-specific factors.
- **Strengths:** pip scSpectra and bioconda. Fully AnnData-native: results go to adata.uns/obsm. Full figure-reproduction repo. Extended by insitro's PertSpectra for Perturb-seq.
- **Limitations:** Needs curated gene sets and cell-type labels, and the results depend on the prior. It is not built for many patients or tumours, or for separating patient-specific from shared signal.
- **Relevance:** Shows the benchmark style Nature Biotech accepted. Competitors were Slalom, scHPF, NMF and expiMap (SpectraReproducibility repo). Real datasets were Bassez breast-cancer immunotherapy (main), Kartha PBMC, Caushi TIL, Salcher TIL and Zhang. Biological payoff came from immunotherapy-response and TIL factors. Known MPs could serve as soft priors in a 'semi-supervised MP' design.

### GBCD: generalized binary covariance decomposition for multi-tumour scRNA-seq (closest published competitor)

*Liu Y, Carbonetto P, Willwerscheid J, Oakes SA, Macleod KF, Stephens M. 'Dissecting tumor transcriptional heterogeneity from single-cell RNA-seq data by generalized binary covariance decomposition.' Nature Genetics 57:263–273 (2025). doi:10.1038/s41588-024-01997-z.* — ✅ verified online — <https://github.com/stephenslab/gbcd>

- **Core idea:** Malignant cells from many tumours are factorised jointly, without batch correction, using empirical-Bayes matrix factorisation (flashier/ebnm) with near-binary memberships. The covariance is decomposed into patient-specific and shared gene expression programs. The workflow repo covers HNSCC, PDAC and simulations.
- **Strengths:** Principled statistics with no batch correction needed. R package plus a reproducible workflowr site.
- **Limitations:** R only. The README gives no numbers on scale or speed, and I did not verify how it handles >1,000 tumours. Interpretation relies on a near-binary membership assumption. I did not find a projection or reference tool.
- **Relevance:** Directly tackles the user's complaint (patient-specific vs shared programs) and is already in a top journal. A new MP method must benchmark against GBCD and state clearly what it adds, e.g. scale, uncertainty, CNV handling, cross-study reproducibility or a reference tool.

### expiMap: biologically informed deep learning to query gene programs in atlases

*Lotfollahi M, Rybakov S, Hrovatin K, Hediyeh-Zadeh S, Talavera-López C, Misharin AV, Theis FJ. Nature Cell Biology (2023). doi:10.1038/s41556-022-01072-x.* — ✅ verified online — <https://github.com/theislab/scarches>

- **Core idea:** A VAE whose latent dimensions are tied to known gene programs through a masked linear decoder, with extra de novo programs. It supports scArches reference mapping, so query data can be read in terms of program activities.
- **Strengths:** Inside the scArches/scverse ecosystem and reused by other groups (e.g., lipiMap from the La Manno lab).
- **Limitations:** Depends on curated gene sets. Interpreting de novo programs is harder. Results depend on the deep-learning setup (seeds, architecture).
- **Relevance:** A competitor to benchmark against (it was benchmarked in Spectra), and a template for 'reference model plus query mapping' of programs within scverse.

### Hotspot: informative gene modules from local autocorrelation

*DeTomaso D, Yosef N. Cell Systems 12(5):446–456.e9 (2021). doi:10.1016/j.cels.2021.04.005.* — ✅ verified online — <https://github.com/YosefLab/Hotspot>

- **Core idea:** Genes are tested for non-random variation over a cell-similarity graph (kNN, lineage or spatial), and informative genes are grouped into modules by local pairwise autocorrelation. Benchmarked against WGCNA and GRNBoost.
- **Strengths:** Works across modalities (lineage, spatial, transcriptome). Python package that is widely reused (Stereopy, TCCA).
- **Limitations:** Works on one dataset. Modules are hard partitions. No mechanism for recurrence across tumours.
- **Relevance:** A non-NMF module-discovery baseline, and a possible gene-filtering step to drop uninformative genes before factorisation. It was a Cell Systems methods paper, which is a realistic fallback venue.

### DRVI: unsupervised deep disentangled representations of single-cell omics

*Moinfar AA, Theis FJ. bioRxiv 2024 (doi:10.1101/2024.11.06.622266). Model merged into scvi-tools ≥1.5.0 as scvi.external.DRVI. I found no journal version.* — ✅ verified online — <https://github.com/theislab/drvi>

- **Core idea:** A structured latent space with additive decoders that splits cellular heterogeneity into interpretable factors (gene programs), with a built-in disentanglement benchmark (DiscreteDisentanglementBenchmark, SMI/SPN metrics).
- **Strengths:** First-class scverse integration (scvi-tools), with documentation and metrics.
- **Limitations:** Still a preprint. Factors are model-dependent and may vary with seed. Recurrence across patients is not modelled explicitly.
- **Relevance:** A competitor reviewers will likely request, and a source of disentanglement metrics. A master's thesis repo already applies DRVI to tumour epithelium together with metaprogram signature collections.

### scITD: tensor decomposition of donor × cell type × gene variation

*Mitchel J, Gordon MG, Perez RK, Biederstedt E, Bueno R, Ye CJ, Kharchenko PV. Nature Biotechnology 43(7):1192–1201 (2025; online 2024). doi:10.1038/s41587-024-02411-z.* — ✅ verified online — <https://doi.org/10.1038/s41587-024-02411-z>

- **Core idea:** Tucker-type tensor decomposition of pseudobulk data (donors × cell types × genes) finds multicellular programs that co-vary across individuals and stratify patient cohorts.
- **Strengths:** Nature Biotech precedent for a 'patient-level program' method with cohort stratification as validation.
- **Limitations:** Works on pseudobulk, so variation within a tumour is lost. Needs many donors.
- **Relevance:** A different view: variation between patients is treated as signal rather than nuisance. It could link malignant MPs to TME programs across patients, a strong biological validation angle.

### DIALOGUE: multicellular programs

*Jerby-Arnon L, Regev A. Nature Biotechnology 40:1467–1477 (2022). doi:10.1038/s41587-022-01288-0.* — ✅ verified online — <https://github.com/livnatje/DIALOGUE>

- **Core idea:** Penalised multi-view CCA plus multilevel modelling finds programs that co-vary across cell types within the same samples or niches.
- **Strengths:** Spatial and single-cell inputs. Ported into scverse (pertpy).
- **Limitations:** Needs multi-cell-type samples. The R implementation is heavy; the Python port in scverse/pertpy is labelled work in progress.
- **Relevance:** A tool for showing that robust MPs couple to TME states, which addresses 'does it reflect real biology?'.

### MetaTiME: pan-cancer meta-components of the tumour immune microenvironment

*Zhang Y et al. (X. Shirley Liu lab). Nature Communications 2023, PMID 37149682 (bio.tools record).* — ✅ verified online — <https://github.com/yi-zhang/MetaTiME>

- **Core idea:** ICA on many tumour scRNA-seq datasets learns reusable 'meta-components' (MeCs). New TME data is annotated by scoring those MeCs.
- **Strengths:** Pip package with a pretrained component reference and tutorials.
- **Limitations:** Covers immune and TME cells, not malignant cells. ICA components are signed and need to be filtered by hand for interpretability.
- **Relevance:** The closest precedent for 'pan-cancer meta-components plus annotation tool'. It landed in Nature Communications, not Nature Methods, which is a useful calibration of how much novelty a resource-plus-tool earns.

### LEMUR: multi-condition latent embedding multivariate regression

*Ahlmann-Eltze C, Huber W. Nature Genetics 2025. doi:10.1038/s41588-024-01996-0.* — ✅ verified online — <https://github.com/const-ae/lemur>

- **Core idea:** Models how a cell's latent embedding changes with sample-level covariates (conditions), giving cluster-free differential expression.
- **Strengths:** Clean statistics with R and Python (pylemur) implementations.
- **Limitations:** Not a gene-program method. Covariates must be known.
- **Relevance:** A statistical template for putting tumour, study and protocol covariates into the model rather than fitting each tumour separately. Relevant if the user builds a hierarchical model with shared vs tumour-specific loadings.

### MOFA-FLEX (successor to MuVI): factor models with prior knowledge and groups

*Qoku A, Rohbeck M, Walter FC, Kats I, Stegle O, Buettner F. bioRxiv 2025, doi:10.1101/2025.11.03.686250. Predecessor MuVI: Qoku & Buettner, arXiv 2022 / AISTATS.* — ✅ verified online — <https://github.com/bioFAM/mofaflex>

- **Core idea:** A modular Bayesian factor framework with structured sparsity, gene-set priors and groups or views (samples, conditions).
- **Strengths:** From the MOFA lineage (bioFAM), with strong scverse ties.
- **Limitations:** Preprint. General purpose, not tailored to many tumours, CNV or recurrence.
- **Relevance:** A possible backbone or competitor: groups could encode tumours, giving shared vs group-specific factors. Reviewers may ask why a general group factor model is not enough.

### Foundation-model routes to programs: scGPT gene embeddings, and a 2026 sparse-autoencoder atlas of Geneformer and scGPT

*Cui H et al. 'scGPT: toward building a foundation model for single-cell multi-omics using generative AI.' Nature Methods 2024 (s41592-024-02201-0). Kendiukhov (2026) 'Sparse autoencoders reveal organized biological knowledge but minimal regulatory logic in single-cell foundation models,' arXiv:2603.02952.* — ✅ verified online — <https://arxiv.org/abs/2603.02952>

- **Core idea:** scGPT clusters gene embeddings into gene programs and modules. The SAE study trains TopK sparse autoencoders on Geneformer V2-316M (82,525 features) and scGPT (24,527 features) residual streams, finding pathway-organised features but little causal regulatory logic.
- **Strengths:** Nature Methods precedent for foundation-model papers that include program extraction. SAE atlases are public.
- **Limitations:** Foundation-model programs inherit biases in the training data, are hard to tie to individual tumours, and the 2026 SAE work finds minimal regulatory logic. Little is known about how reproducible foundation-model features are across tumours.
- **Relevance:** Reviewers in 2026 may ask whether foundation-model embeddings or SAE features give more robust MPs than NMF. Including a foundation-model baseline, or a foundation-model-informed variant, pre-empts this.

### PerturbNMF / cNMF_benchmarking (Engreitz lab, IGVF): perturbation-grounded evaluation of gene programs

*Engreitz lab, Zenodo concept DOI 10.5281/zenodo.21878267 (2025–2026). The earlier cNMF_benchmarking repo is now superseded. Related IGVF tf_perturb_seq WG2 gene-program working group (2026).* — ✅ verified online — <https://github.com/EngreitzLab/PerturbNMF>

- **Core idea:** End-to-end pipeline: CPU or GPU cNMF, then nine program-quality metrics (perturbation sensitivity with calibration against non-targeting-guide nulls, GO, motif and trait enrichment, categorical association, stability), then K-selection and LLM-assisted annotation. The IGVF working group classifies programs as shared_all, shared_partial or dataset-specific by meta-consensus clustering.
- **Strengths:** Open, GPU-scalable, and offers a ready metric suite to adopt or extend.
- **Limitations:** Not published as a peer-reviewed paper (Zenodo only). Aimed at Perturb-seq and differentiation systems, not tumours.
- **Relevance:** Shows how to give 'real biology' an objective, falsifiable definition: programs should respond to the right perturbations, with calibrated statistics. A Nature Methods MP paper could use Perturb-seq in cancer cell lines or organoids (and Kinker-style cell-line heterogeneity) as ground truth.

### Existing re-implementations of the Gavish MP pipeline: GeneNMF, sciNMF, mpnmf, pyscalop/scalop

*GeneNMF (Carmona lab, R; cited via Yerly L et al. bioRxiv 2024, doi:10.1101/2024.05.31.596823). sciNMF (Zhang Z et al., Brief Bioinform 2024/2025, doi:10.1093/bib/bbaf028). mpnmf (K.-W. Kim, Python, April 2026). pyscalop (J. Laffy, Python port of the Tirosh lab's scalop, May 2026).* — ✅ verified online — <https://github.com/carmonalab/GeneNMF>

- **Core idea:** All of them do per-sample NMF over multiple ranks, then robust-program filtering (overlap or IQR thresholds), then overlap or cosine clustering into MPs. GeneNMF uses cosine similarity of full weight vectors with cumulative-weight gene selection (v0.6 'more robust' MP calculation). mpnmf replaces random restarts with deterministic NNDSVDa initialisation and offers HEG or HVG gene-selection modes.
- **Strengths:** Easy to run in Seurat/R and Python, and widely adopted.
- **Limitations:** None models uncertainty, patient-specific vs shared signal, CNV or technical confounders. All remain threshold-driven, and none comes with a systematic benchmark.
- **Relevance:** This space is crowded: another 'faster or deterministic Gavish' package will be judged incremental. These tools are good, easy baselines to include (GeneNMF has about 210 stars and is widely used), and their disagreement with each other can motivate the paper (the same data giving different MPs).

### CellProgramMapper: R projection of single-cell data onto reference GEPs

*Liu Z. R package, 2026 (software citation only).* — ✅ verified online — <https://github.com/Zaoqu-Liu/CellProgramMapper>

- **Core idea:** Projects cells onto reference programs (e.g., TCAT.V1) by NNLS (C++ coordinate-descent or active-set solvers), with Seurat input and custom references built from several cNMF runs.
- **Strengths:** Seurat-native.
- **Limitations:** Software only, with no paper and no uncertainty quantification. I found no pre-built malignant MP reference.
- **Relevance:** Shows the 'project onto reference programs' niche is filling up. A malignant-MP reference and projector would need added value such as uncertainty, novelty detection or out-of-reference programs.

### TCCA (Therapeutic Cancer Cell Atlas): a 2026 pan-cancer MP resource built on the standard method

*González-Bermejo M, Serrano-Ron L, García-Martín S, et al. bioRxiv 2026, doi:10.64898/2026.02.18.706316 (CNIO). Also an ECCB 2026 poster (A-T.29).* — ✅ verified online — <https://github.com/cnio-bu/tcca>

- **Core idea:** 1,089,024 malignant and 750,218 non-malignant cells from 36 studies (the poster says about 1.8M transcriptomes, 537 patients, 183 cell lines, 34 cancer types). NMF MPs follow Gavish 2023, with UCell scoring, Beyondcell drug sensitivity, SCEVAN CNV calls and 12 TME archetypes.
- **Strengths:** Large, harmonised data with drug-response annotation and Snakemake reproducibility.
- **Limitations:** Uses the standard threshold pipeline and does not claim to improve robustness.
- **Relevance:** Not a methods scoop, but a competing resource and a ready test bed. Re-deriving MPs on TCCA or 3CA with the new method and comparing would be a natural benchmark.

### Unattributed 2026 pan-cancer MP catalogue (JerryChen06/Metaprogram-Supplemental-Information)

*Public GitHub repo created 2026-08-03 and updated 2026-09-30. No README, authors, lab or preprint identified.* — ✅ verified online — <https://github.com/JerryChen06/Metaprogram-Supplemental-Information>

- **Core idea:** HTML supplements: MP gene lists ranked by NMF frequency, MP heatmaps (including a carcinoma subset), programs per cancer type across 16 categories (about 60 or more MPs; 2–22 per type), and 'interesting hybrids'. MP names include cell-cycle variants, stress, EMT variants, hypoxia, MYC, IFN/MHC-II, respiration, PDAC-classical, skin pigmentation and 'Glioma single-nucleus'.
- **Strengths:** Concrete current snapshot of what an expanded MP catalogue looks like.
- **Limitations:** Provenance unknown. It may be a student project or the supplement of a forthcoming paper.
- **Relevance:** A sign that a larger, updated pan-cancer MP manuscript may be close. The 'Glioma single-nucleus' MP suggests protocol-driven artefacts still survive the standard pipeline, which is a concrete motivating example for the user.

### Prior pan-cancer program discovery: Kinker 2020 (cell lines) and Barkley 2022 (tumours)

*Kinker GS, Greenwald AC, Tal R, ... Tirosh I. 'Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity.' Nature Genetics 52:1208–1218 (2020), doi:10.1038/s41588-020-00726-6. Barkley D, Moncada R, ... Yanai I. 'Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment.' Nature Genetics 54:1192–1201 (2022), doi:10.1038/s41588-022-01141-9.* — ✅ verified online — <https://doi.org/10.1038/s41588-022-01141-9>

- **Core idea:** Kinker: recurrent heterogeneous programs (RHPs) across cancer cell lines, compared with tumours. Barkley: NMF-derived cancer cell states recurring across tumour types, with spatial co-localisation to TME states (data in GSE203612).
- **Strengths:** Established biology that reviewers know well.
- **Limitations:** Same threshold- and overlap-based logic. No formal robustness or uncertainty analysis.
- **Relevance:** Cell lines give a controlled system for perturbation validation. Spatial data (Barkley) gives orthogonal validation of 'real' states. Both are expected citations and comparison points.

### Benchmarking standards reviewers apply: Weber et al. 2019 and Sonrel et al. 2023

*Weber LM et al. 'Essential guidelines for computational method benchmarking.' Genome Biology 2019, doi:10.1186/s13059-019-1738-8. Sonrel A et al. 'Meta-analysis of (single-cell method) benchmarks reveals the need for extensibility and interoperability.' Genome Biology 24 (2023), doi:10.1186/s13059-023-02962-5.* — ✅ verified online — <https://doi.org/10.1186/s13059-019-1738-8>

- **Core idea:** Benchmarks should define scope; use simulated and real data with ground truth; choose competitors neutrally with fair parameter tuning; use several metrics, including runtime and scalability; and be reproducible and extensible (e.g., Open Problems-style).
- **Strengths:** Widely cited and used by editors and reviewers as a checklist.
- **Limitations:** These are general guidelines, not specific to gene programs. No neutral benchmark for program or MP discovery exists yet.
- **Relevance:** An open, extensible MP-discovery benchmark (simulations with known shared, patient-specific, CNV-driven and artefact programs, plus perturbation and spatial ground truth) would itself be a strong Nature Methods contribution alongside the method.

### Reproducibility of gene signatures: random-signature nulls and stability selection

*Venet D, Dumont JE, Detours V. 'Most random gene expression signatures are significantly associated with breast cancer outcome.' PLoS Comput Biol 7:e1002240 (2011). Ein-Dor L, Zuk O, Domany E. PNAS 103:5923–5928 (2006). Meinshausen N, Bühlmann P. 'Stability selection.' JRSS-B 72:417–473 (2010). Citations seen in online reference lists and bioconda/SigCheck metadata.* — ✅ verified online — <https://doi.org/10.1371/journal.pcbi.1002240>

- **Core idea:** Signatures confounded with proliferation or other dominant axes look 'significant' by chance. Robust gene lists need far more samples than are usually used. Stability selection controls false selections through subsampling.
- **Strengths:** Simple, well-understood statistics that transfer directly to MP gene-membership calls.
- **Limitations:** From the bulk and biomarker era, not single-cell.
- **Relevance:** Gives a principled basis for the user's robustness claim: MP membership with subsampling-based error control (genes and tumours), random and matched-null programs, and checks for confounding by cycling, CNV or quality axes. These are classic reviewer objections to signature papers.

**Gaps identified by this sweep**

- No current MP method gives calibrated uncertainty: no FDR for a gene belonging to an MP, no significance for a program recurring across tumours, no confidence intervals on MP usage. Gavish, 3CA_v2, GeneNMF, sciNMF and mpnmf all rely on hard overlap or top-N cutoffs (e.g., top-50 genes, at least 10 or 35 shared genes, at least 5 programs).
- No method explicitly removes or labels copy-number-, clone- or patient-genotype-driven programs during MP discovery. GBCD separates patient-specific from shared programs statistically but does not use inferred CNV or clone structure (inferCNV, numbat, SCEVAN) as a covariate.
- No neutral, extensible benchmark for metaprogram discovery exists: simulations with known shared, cancer-type-specific, patient-specific, CNV-driven, protocol and dissociation-artefact programs, plus held-out-study reproducibility, perturbation ground truth (as in PerturbNMF) and spatial co-localisation. That benchmark could be a co-contribution.
- Technical-artefact programs (dissociation and stress response, ribosomal or mitochondrial, ambient RNA, doublets, single-nucleus vs single-cell protocol) are handled ad hoc. A 2026 MP catalogue still contains a 'Glioma single-nucleus' MP.
- No published malignant-cell MP reference with a projection tool, uncertainty and novelty detection (a 'starCAT for malignant MPs'). Projection tools (starCAT, CellProgramMapper) exist only for T-cell or immune references.
- Rank (K) selection, gene-universe choice (top 7,000 by mean vs HVG), centring and top-50 truncation strongly shape MPs, but no study I found systematically measures MP sensitivity to these choices or reports a reproducibility index.
- Hierarchy (pan-cancer → lineage → cancer type → tumour) and 'hybrid' MPs are not modelled formally; current outputs are flat lists produced by greedy clustering.
- Cross-cohort generalisation is rarely tested, e.g. learn MPs on 3CA and replicate them on TCCA or other independent cohorts, or show that held-out studies recover the same MPs.
- Foundation-model-derived programs (scGPT embeddings, SAE features) have not been rigorously compared with NMF MPs for robustness across tumours.
- Search limitation: the WebSearch budget was used up before this agent started, and PubMed, Europe PMC, Nature, bioRxiv, arXiv, Crossref and OpenAlex were blocked by egress policy. The scoop search therefore covered only GitHub-indexed content (code search over READMEs, bib files and digest repos). Unindexed 2025–2026 bioRxiv preprints on robust or consensus MP discovery may exist and must be checked by hand before committing.

**Notes:** VERIFICATION METHOD: WebSearch was exhausted (200/200) when this agent started, and publisher and literature APIs were egress-blocked. All verified_online=true items were confirmed this session from live GitHub content: official repos, CITATION files, READMEs, or bib entries in third-party repos that give journal, volume and DOI. Details marked 'from memory' (e.g., TCAT ~46 GEPs / ~1.7M cells, Gavish 41 MPs) were not re-verified and should be checked before citing. DRVI and MOFA-FLEX are confirmed only as preprints; I found no journal versions. I could not confirm where scDEF (cbg-ethz, hierarchical deep exponential families with batch handling) was published; it is a possibly relevant competitor to check manually.

WHAT IT TAKES (pattern across 2022–2026 Nature Methods, Nature Biotechnology, Nature Genetics and Cell Systems program or factor papers):
(1) A clear conceptual advance that a generic tool cannot match: TCAT = a cross-dataset reproducible reference plus projection; Spectra = knowledge-guided, cell-type-specific factors; GBCD = shared vs patient-specific programs without batch correction; scITD = donor-level multicellular programs. A faster or deterministic Gavish re-implementation would not qualify; GeneNMF, sciNMF, mpnmf and pyscalop already fill that niche.
(2) Benchmarks: simulations with known ground truth plus several real datasets, with 4–8 competitors. For an MP method, expect GBCD, cNMF (pooled or per-tumour), Gavish/3CA, GeneNMF, Spectra or expiMap, DRVI or scVI-family, MOFA-FLEX/MuVI, Hotspot, and possibly a foundation-model baseline. Metrics should cover recovery of planted programs, stability over seeds, ranks, subsampled cells and tumours, held-out-study replication, and runtime and memory at 10^6–10^7 cells.
(3) Orthogonal biological validation: perturbation (Perturb-seq or CRISPR in cancer lines or organoids, scored with PerturbNMF-style calibrated tests), spatial co-localisation (Barkley-style), protein (CITE-seq or IF), clinical association (immunotherapy response, survival) with matched-random-signature nulls (Venet). Ideally include one new experiment showing that a robust MP the old pipeline missed, or a spurious MP it removes, behaves as predicted.
(4) Software and resource: pip or conda package, AnnData/scverse- or Seurat-compatible (Spectra writes to adata.uns/obsm; DRVI was merged into scvi-tools), documentation and tutorials, a reproducibility repo (as SpectraReproducibility), a downloadable reference catalogue on Zenodo (as TCAT.V1), and a projection function.
(5) Scale: re-derive MPs on 3CA (about 1,163+ tumours, with 3CA_v2 adding 64 datasets) and replicate on an independent atlas such as TCCA (about 1.1M malignant cells, 36 studies).

COMPETITIVE LANDSCAPE / SCOOP RISK (as of 2026-10-02, GitHub-visible only):
(a) GBCD (Nature Genetics 2025) already holds the 'shared vs patient-specific tumour programs' claim, so differentiation is essential: uncertainty-calibrated MPs, CNV- and artefact-aware modelling, hierarchy, scale, and a reference plus projection with novelty detection.
(b) The Tirosh lab has unpublished 3CA_v2 MP code using the same greedy procedure, so an expanded official MP catalogue is likely coming; benchmark against it.
(c) An unattributed 2026 repo (JerryChen06) posts a ~60-MP pan-cancer catalogue, suggesting another manuscript is close.
(d) TCCA (bioRxiv 2026) is a large resource built on the standard method, not a methods competitor.
(e) Many lightweight MP packages (GeneNMF, sciNMF, mpnmf, pyscalop, CellProgramMapper) appeared in 2024–2026, so incremental tooling is saturated.
I found no 2025–2026 GitHub-visible preprint claiming a statistically principled, uncertainty-aware, confounder-aware pan-cancer MP method. Absence in GitHub search is weak evidence, so do a manual bioRxiv and Google Scholar check, e.g. 'metaprogram robust', 'consensus meta-program', 'pan-cancer NMF uncertainty', 'GBCD'.

SUGGESTED FRAMING: 'MP discovery as a statistical inference problem.' A hierarchical, probabilistic consensus across seeds, ranks, tumours and studies, with error-controlled gene membership (stability selection), explicit covariates for CNV, clone, protocol and quality, MP families organised by cancer type, held-out-study replication, a public benchmark suite, a versioned MP reference with a projector (uncertainty plus novelty detection), and perturbation and spatial validation. That combination is not covered by any single existing tool found here.
