# R versus port parity check

The mpaudit port of the Gavish/3CA pipeline (`mpaudit.gavish`, `mpaudit.nmf`) minimises the same
sNMF/R objective as R's `NMF::nmf(method = "snmf/r")`, but with a different optimiser (HALS
instead of alternating NNLS) and a different stopping rule. Its programs are therefore not
bit-identical to R's. These scripts measure how far apart they are on real samples, which has to
be done before quoting absolute numbers (number of MPs, MP genes) as 3CA-equivalent (review of
2026-10-16, section 1 and next step 4).

R could not be installed where these scripts were written, so the R side has never been executed
by the author: `run_r_snmf.R` follows the NMF 0.30.1 sources (renozao/NMF) but expect to fix small
things on first use. The Python side is tested end to end (`tests/test_parity_scripts.py`), with
"R output" produced by mpaudit itself in the R file format.

## Files

| file | what it does |
|---|---|
| `export_matrices.py` | writes each sample's NMF input exactly as mpaudit builds it (`prepare_matrix`) as Matrix Market + `genes.txt` + `cells.txt`, and `manifest.json` with the parameters |
| `run_r_snmf.R` | runs the 3CA call `nmf(x, rank = 4:9, method = "snmf/r", nrun = 10, seed = <seed>)` on each exported matrix; writes every rank's basis W, the best-fit objective, `sessionInfo()` and package versions |
| `compare.py` | fits (or loads) the port on the same samples and parameters and compares programs, robust programs and MPs; writes `report.md` and CSVs with a verdict |
| `parity_io.py` | the file formats shared by the three |

## Running it

Requirements: the mpaudit environment (Python) and R >= 3.6 with `install.packages(c("NMF", "Matrix"))`
(if NMF's installation asks for Bioconductor's Biobase: `BiocManager::install("Biobase")`).

```bash
# 1. Export the 5 samples with most NMF cells (or --samples S1 S2 S3), using the 'mpaudit p1' config.
python scripts/parity/export_matrices.py cohort.h5ad --config configs/p1_3ca_example.yaml \
    --out parity/export --n 5

# 2. R. Smoke test first (one sample, one rank, one run) to time it and catch setup problems:
Rscript scripts/parity/run_r_snmf.R parity/export --out parity/r_smoke --samples 01_<name> --ranks 4 --nrun 1
#    then the real run (ranks, nrun and seed come from the export):
Rscript scripts/parity/run_r_snmf.R parity/export --out parity/r_out --cores 8
#    optional noise floor for R: a second run at another seed
Rscript scripts/parity/run_r_snmf.R parity/export --out parity/r_out_seed2 --seed 2 --cores 8

# 3. Compare (fits the port on the same samples; --cache-dir reuses an 'mpaudit p1' NMF cache).
python scripts/parity/compare.py --export parity/export --r-out parity/r_out --out parity/report \
    --h5ad cohort.h5ad --python-replicate-seed 1 --r-replicate parity/r_out_seed2 --n-jobs 4
```

With a multithreaded BLAS (OpenBLAS, MKL), set `OPENBLAS_NUM_THREADS=1` (or `MKL_NUM_THREADS=1`)
when using `--cores`, to avoid oversubscribing the machine; `versions.csv` records the BLAS, which
can change the last digits of every fit.

`compare.py` writes the port's bases to `<out>/python_bases`; pass `--python-bases parity/report/python_bases`
to re-run the comparison without refitting. `run_r_snmf.R` skips samples that already have a
`fits.csv`, so an interrupted run can simply be restarted (`--force` refits).

Sample choice: the largest samples give the most stable programs and are what matters for an
atlas, but they are also the slowest; 3-5 samples from different studies are enough for the
program level. Use the same config as the audit (the export records it; `min_genes_per_cell`,
`n_universe`, ranks and nrun all change the matrix or the fits).

## Expected runtime

- Export: about 10 s and 100 MB (gzipped) per 7,000 x 3,000 sample here (7.9 million non-zero
  entries written in 7.4 s), plus reading the .h5ad.
- Port: one sNMF/R run on a 7,000 x 2,000 matrix took 6-8 s on one core here (115-135 iterations at
  ~50 ms each, `tol=1e-5`), so ranks 4-9 x nrun 10 is roughly 5-15 min per sample per core;
  `--n-jobs` parallelises over samples. A replicate doubles it. Lowering `tol` multiplies it.
- R: not measured (R unavailable here). R's snmf/r solves exact NNLS problems at each iteration
  and its stopping rule ran ~900 iterations when applied to the port's iterates on the toy cohort
  (see `mpaudit/nmf.py`), so expect it to be several times slower than the port per run: plan for
  hours per sample on one core, divided by `--cores` (the nrun restarts run in parallel). The smoke
  test's `elapsed_sec` in `fits.csv`, times 6 ranks x 10 runs, gives your own estimate.
- Compare: the port's fits plus seconds; each profile objective (H re-fitted for one basis) took
  0.08 s at 7,000 x 3,000, rank 9, plus 1.6 s per sample to read the matrix (`--no-objective` skips them).

## Reading the report

`report.md` starts with the verdict and the input checks, then one table per level. The CSVs
hold the details: `programs.csv` (every paired program: overlap, cosine), `ranks.csv` (per sample and
rank, with objectives), `robust.csv`, `mps.csv`, `summary.json` (all numbers and thresholds).

Levels, all computed on R's bases and on the port's bases for the same samples and parameters:

1. **Programs**: per sample and rank, R's k programs and the port's k programs are paired
   one-to-one by the Hungarian assignment maximising the summed top-50 overlap. The statistic is
   the matched overlap (0-50). The **profile objective** is the sNMF/R objective R reports, with H
   re-fitted to the exported matrix for each basis; lower is a better fit.
2. **Robust programs**: the 3CA robust-program filter (`mpaudit.gavish.assemble`) on each side;
   per sample, robust programs are paired, and one counts as reproduced when its partner shares
   >= `intra_min` (35) genes.
3. **MPs**: MP clustering on each side, MPs paired by Jaccard (`mpaudit.compare.match_sets`,
   >= 0.3); per pair, shared genes, member counts and member agreement (R's member programs
   mapped to their level-1 partners).

Verdict:

| verdict | rule | meaning |
|---|---|---|
| PASS | median matched program overlap >= 45/50; every robust program reproduced, same count per sample; same number of MPs, each paired with >= 45/50 shared genes and the same number of members | absolute numbers from the port can be quoted as 3CA-equivalent for this setting |
| WITHIN NOISE | not PASS, but R vs port agrees at least as well as every seed replicate (port vs port, R vs R) at every level: median overlap within 2 genes, robust-program agreement (reproduced / larger count) and MP agreement (identical / larger count) at least the replicate's | the solvers differ no more than reruns of either; quote numbers with that caveat |
| FAIL | otherwise | see below |
| INVALID | R did not read the exported matrix, or a side's genes differ from the export | nothing else in the report is meaningful; re-export with the current code |

On a FAIL, look at the profile objectives. If the port's basis fits worse at most sample-ranks
(`py_vs_r_pct` > 0), however small the difference, the port stops short of R: lower `tol` (e.g.
1e-7) or raise `max_iter` in the config, re-export and compare again. If the signs are mixed, the
two solvers reach different but equally good optima (NMF is not identifiable); robust-program and
MP agreement then say whether MP-level conclusions survive the solver change, and absolute counts
must not be quoted as 3CA-equivalent.

A FAIL is only interpretable against a seed replicate: always pass `--python-replicate-seed`
(cheap) and, if possible, `--r-replicate`. Toy illustration (5 simulated samples, 400-gene
universe, ranks 4-6, nrun 2; standing in for R, the port itself at `tol=1e-8` and another seed):
the port against itself at another seed gave a median overlap of 49 (10th percentile 43), 19 of 22
robust programs reproduced and 2 of 5 MPs identical; the stand-in gave 48, 17 of 22 and 0 of 5, a
FAIL, and its bases fitted better at 12 of 15 sample-ranks, by only 0.002-0.03% (worse at the other 3). Early stopping
costs little objective but moves top genes, robust programs and MPs.

Caveats:

- **MP stage settings.** With 3-5 samples the 3CA rule (a founder needs > 5 partners, i.e.
  >= 7 samples) forms no MP. By default `compare.py` lowers `min_group_size` so that an MP must recur
  in all but one exported sample (at least 3) and switches off the 3CA v2 post-filters and gene order;
  the report states the values used, both sides get the same ones, and `--mp-as-configured`
  restores the exported values. MP parity on 5 samples tests the clustering on these inputs, not
  the atlas MPs.
- **Sample-set dependence.** The robust filter needs recurrence in another sample, so robust
  programs and MPs depend on which samples were exported; that is the same for both sides.
- **Overlaps are not independent** across one sample's ranks (the same program reappears at
  several ranks); the median is a summary, not a test.
- **What the export fixes.** Cell selection, transform, gene universe and centring are taken from
  mpaudit, so the comparison isolates the NMF solver and stopping rule. It says nothing about the
  parts of 3CA that are not public (`Return_All_cells_Separated` cell QC); see `mpaudit/gavish.py`.
- **Float parsing.** R's decimal parser may differ from Python's in the last bit; the input check
  allows a relative difference of 1e-9 in the matrix sum and max.
- **Cell selection is re-implemented** in `export_matrices.py` from public functions (malignant
  cells with >= `min_genes_per_cell` detected genes, samples with >= `min_cells`). If
  `mpaudit.gavish.compute_programs` changes its selection, the input check reports the port's
  genes as differing from the export (INVALID) instead of comparing different matrices.
