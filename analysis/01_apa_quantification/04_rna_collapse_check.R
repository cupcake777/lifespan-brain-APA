#!/usr/bin/env Rscript
# Within-batch collapse check for an expression matrix. Usage: Rscript rna_collapse_check.R matrix.tsv meta.tsv [n_features]

suppressPackageStartupMessages(library(data.table))

min_unique_frac <- 0.5
a <- commandArgs(trailingOnly = TRUE)
n_features <- if (length(a) > 2L) as.integer(a[[3L]]) else 2000L

M <- fread(a[[1L]], nrows = n_features, check.names = FALSE)
meta <- fread(a[[2L]], check.names = FALSE)
X <- as.matrix(M[, -1L, with = FALSE])
storage.mode(X) <- "double"
b <- meta$batch[match(colnames(X), meta$sample)]
stopifnot(!anyNA(b), all(is.finite(X)))

res <- rbindlist(lapply(split(seq_along(b), b), function(j) data.table(
  n_samples = length(j),
  median_unique_frac = median(apply(round(X[, j, drop = FALSE], 5), 1,
                                    function(v) length(unique(v)) / length(v))),
  median_within_sd = median(apply(X[, j, drop = FALSE], 1, sd))
)), idcol = "batch")
res[, `:=`(
  matrix = a[[1L]],
  collapse_fail = median_unique_frac < min_unique_frac | !is.finite(median_within_sd) | median_within_sd < 1e-6
)]

print(res)
if (any(res$collapse_fail)) {
  cat("FAIL: within-batch collapse detected\n")
  quit(status = 1L)
}
cat("PASS\n")
