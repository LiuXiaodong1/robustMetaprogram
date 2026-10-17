#!/usr/bin/env Rscript
# Run R's NMF::nmf(method = "snmf/r") on the matrices written by export_matrices.py: the R side of
# the mpaudit R-versus-port parity check (see README.md in this folder).
#
#   Rscript scripts/parity/run_r_snmf.R <export_dir> [--out DIR] [--seed N] [--cores N]
#       [--samples DIR1,DIR2] [--ranks 4,5,6,7,8,9] [--nrun 10] [--force]
#
# The call is the one given in the Gavish et al. 3CA README (tiroshlab/3ca, ITH_hallmarks/README.md):
#   NMFs_per_sample = nmf(x = Expression_matrix, rank = 4:9, method="snmf/r", nrun = 10)
# with a numeric seed added so the run can be repeated (NMF then draws one L'Ecuyer RNG stream per
# run; the result does not depend on the parallel backend or the number of cores). With several
# ranks, nmf() hands over to nmfEstimateRank(), which fits each rank in turn with the same seed and
# also computes rank-survey measures (consensus, silhouettes; not used here). Ranks, nrun and seed
# default to the export's r_params.tsv (the mpaudit config); --ranks/--nrun are for smoke tests only.
#
# Output, per sample (<out>/<sample dir>/, the sample dir as in the export):
#   W_rank<k>.csv    basis W of the best of nrun runs (genes x k): column "gene", then 1..k
#   fits.csv         per rank: nrun, objective (residuals() of the best fit: NMF's .snmf.objective,
#                    0.5 * (||A - WH||^2 + eta ||W||^2 + beta sum colSums(H)^2), eta = max(A)),
#                    iterations of the best run, elapsed seconds for all nrun runs
#   input_check.csv  dimensions, sum and max of the matrix as R read it
#   warnings.txt     warnings raised by nmf(), if any
# and once per run: sessionInfo.txt and versions.csv (R, NMF and Matrix versions, BLAS, settings).
# A sample whose fits.csv exists is skipped unless --force, so an interrupted run can be resumed.
# Requires the NMF and Matrix packages. Written against the NMF 0.30.1 sources (renozao/NMF) and never
# executed by its author (R was unavailable); expect small fixes on first use.

suppressPackageStartupMessages({
  library(NMF)
  library(Matrix)
})

usage <- "Usage: Rscript run_r_snmf.R <export_dir> [--out DIR] [--seed N] [--cores N] [--samples D1,D2] [--ranks 4,5] [--nrun N] [--force]"

parse_args <- function(argv) {
  opts <- list(out = NULL, seed = NULL, cores = NULL, samples = NULL, ranks = NULL, nrun = NULL, force = FALSE)
  pos <- character(0)
  i <- 1
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (a == "--force") {
      opts$force <- TRUE
    } else if (startsWith(a, "--")) {
      key <- substring(a, 3)
      if (!key %in% names(opts) || i == length(argv)) stop(usage, call. = FALSE)
      opts[[key]] <- argv[[i + 1]]
      i <- i + 1
    } else {
      pos <- c(pos, a)
    }
    i <- i + 1
  }
  if (length(pos) != 1) stop(usage, call. = FALSE)
  opts$export <- pos[[1]]
  opts
}

split_ints <- function(x) as.integer(strsplit(x, ",", fixed = TRUE)[[1]])

opts <- parse_args(commandArgs(trailingOnly = TRUE))
export_dir <- opts$export
kv <- read.delim(file.path(export_dir, "r_params.tsv"), colClasses = "character", quote = "")
kv <- setNames(kv$value, kv$key)
ranks <- split_ints(if (is.null(opts$ranks)) kv[["ranks"]] else opts$ranks)
nrun <- as.integer(if (is.null(opts$nrun)) kv[["nrun"]] else opts$nrun)
seed <- as.integer(if (is.null(opts$seed)) kv[["seed"]] else opts$seed)
out_dir <- if (is.null(opts$out)) file.path(export_dir, "..", "r_out") else opts$out
samples <- read.delim(file.path(export_dir, "samples.tsv"), colClasses = "character", quote = "")
if (!is.null(opts$samples)) {
  wanted <- strsplit(opts$samples, ",", fixed = TRUE)[[1]]
  unknown <- setdiff(wanted, c(samples$dir, basename(samples$dir)))
  if (length(unknown)) stop("Unknown sample dirs: ", paste(unknown, collapse = ", "), call. = FALSE)
  samples <- samples[samples$dir %in% wanted | basename(samples$dir) %in% wanted, , drop = FALSE]
}
# --cores 1 runs the nrun restarts sequentially; N > 1 uses N doParallel workers; unset keeps
# nmf()'s default backend (as in the 3CA call).
pbackend <- if (is.null(opts$cores)) NULL else if (as.integer(opts$cores) <= 1) "seq" else as.integer(opts$cores)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt17 <- function(x) sprintf("%.17g", x)

write_basis <- function(W, path) {
  vals <- matrix(fmt17(W), nrow = nrow(W))
  out <- data.frame(gene = rownames(W), vals, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(out) <- c("gene", seq_len(ncol(W)))
  write.csv(out, path, row.names = FALSE, quote = 1)
}

# Called with explicit arguments (not do.call) so that NMF's warnings, which deparse the argument
# name, stay short.
run_nmf <- function(A) {
  if (is.null(pbackend)) {
    nmf(x = A, rank = ranks, method = "snmf/r", nrun = nrun, seed = seed)
  } else {
    nmf(x = A, rank = ranks, method = "snmf/r", nrun = nrun, seed = seed, .pbackend = pbackend)
  }
}

# A diagnostic number, or NA when the accessor fails or returns nothing (keeps fits.csv rectangular).
safely <- function(expr) {
  v <- tryCatch(expr, error = function(e) NULL)
  if (length(v) == 1) as.numeric(v) else NA_real_
}

for (i in seq_len(nrow(samples))) {
  sdir <- samples$dir[[i]]
  odir <- file.path(out_dir, sdir)
  if (file.exists(file.path(odir, "fits.csv")) && !opts$force) {
    message(sprintf("[r] %s: fits.csv exists, skipping (use --force to refit)", sdir))
    next
  }
  dir.create(odir, recursive = TRUE, showWarnings = FALSE)
  A <- as.matrix(readMM(gzfile(file.path(export_dir, sdir, "matrix.mtx.gz"))))
  rownames(A) <- readLines(file.path(export_dir, sdir, "genes.txt"))
  colnames(A) <- readLines(file.path(export_dir, sdir, "cells.txt"))
  write.csv(data.frame(n_genes = nrow(A), n_cells = ncol(A), sum = fmt17(sum(A)), max = fmt17(max(A))),
            file.path(odir, "input_check.csv"), row.names = FALSE)
  message(sprintf("[r] %s: %d genes x %d cells; ranks %s, nrun %d, seed %d", sdir, nrow(A), ncol(A),
                  paste(ranks, collapse = ","), nrun, seed))
  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  res <- withCallingHandlers(run_nmf(A), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  elapsed_total <- proc.time()[["elapsed"]] - t0
  if (length(warns)) writeLines(warns, file.path(odir, "warnings.txt"))
  fits <- if (inherits(res, "NMF.rank")) res$fit else setNames(list(res), as.character(ranks))
  rows <- list()
  for (r in ranks) {
    f <- fits[[as.character(r)]]
    if (!isNMFfit(f)) stop(sprintf("%s: rank %d failed: %s", sdir, r, paste(format(f), collapse = " ")), call. = FALSE)
    W <- basis(f)
    if (is.null(rownames(W))) rownames(W) <- rownames(A)
    write_basis(W, file.path(odir, sprintf("W_rank%d.csv", r)))
    rows[[length(rows) + 1]] <- data.frame(
      dir = sdir, rank = r, nrun = nrun(f), objective = fmt17(residuals(f)),
      n_iter = safely(niter(minfit(f))), elapsed_sec = safely(runtime.all(f)[["elapsed"]]), seed = seed)
  }
  write.csv(do.call(rbind, rows), file.path(odir, "fits.csv"), row.names = FALSE)
  message(sprintf("[r] %s: done in %.0f s", sdir, elapsed_total))
}

si <- sessionInfo()
writeLines(capture.output(print(si)), file.path(out_dir, "sessionInfo.txt"))
versions <- data.frame(
  key = c("R", "NMF", "Matrix", "BLAS", "LAPACK", "platform", "seed", "ranks", "nrun", "pbackend", "date"),
  value = c(R.version.string, as.character(packageVersion("NMF")), as.character(packageVersion("Matrix")),
            if (is.null(si$BLAS)) NA else si$BLAS, if (is.null(si$LAPACK)) NA else si$LAPACK,
            R.version$platform, seed, paste(ranks, collapse = ","), nrun,
            if (is.null(pbackend)) "default" else as.character(pbackend), format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  stringsAsFactors = FALSE)
write.csv(versions, file.path(out_dir, "versions.csv"), row.names = FALSE)
message(sprintf("[r] wrote %s", out_dir))
