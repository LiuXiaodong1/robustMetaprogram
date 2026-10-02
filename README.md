# robustMetaprogram

Planned package for robust, calibrated discovery of tumour metaprograms (MPs) from
single-cell RNA-seq. It is meant to replace the threshold-and-curation pipeline
(per-tumour NMF → top-50 lists → overlap filters → greedy clustering → manual removal)
with MPs that have effect sizes, held-out replication, gene-level error control and a
typed provenance: intrinsic biology, TME-shared, genetic dosage, protocol-modulated or technical.

Status: **design stage**. See [docs/](docs/README.md), starting with
[docs/02_TESSERA_proposal.md](docs/02_TESSERA_proposal.md).
