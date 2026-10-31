#!/usr/bin/env Rscript
# Bridge-pair agreement of batch-correction methods, with a donor-paired bootstrap of ComBat_REF vs limma.

suppressPackageStartupMessages(library(data.table))

# --- paths ---
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.local.sh)")
batch_dir <- file.path(Sys.getenv("APA_ROOT"), "data/01_quant/batch")
checkpoint <- readRDS(file.path(batch_dir, "batch_correction_checkpoint.rds"))
method_files <- c(
  Raw = "method_matrix_Raw.rds",
  ComBat_RAW = "method_matrix_ComBat_RAW.rds",
  ComBat_REF = "method_matrix_ComBat_REF.rds",
  ComBat_REF_NP = "method_matrix_ComBat_REF_NP.rds",
  ComBat_REF_meanonly = "method_matrix_ComBat_REF_meanonly.rds",
  limma = "method_matrix_limma.rds",
  `RUV-III` = "method_matrix_RUV-III.rds",
  RUVs = "method_matrix_RUVs.rds"
)
bootstrap_iterations <- 10000L
bootstrap_seed <- 20261012L

# --- helpers ---
ccc <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  2 * cov(x, y) / (var(x) + var(y) + (mean(x) - mean(y))^2)
}

measure_pairs <- function(mat, pairs, method) {
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    x <- mat[, pairs$original_sample[i]]
    y <- mat[, pairs$batch9_replicate[i]]
    data.table(
      pair_id = i,
      original_sample = pairs$original_sample[i],
      batch9_replicate = pairs$batch9_replicate[i],
      method = method,
      pearson = cor(x, y),
      ccc = ccc(x, y),
      mae = mean(abs(x - y)),
      rmse = sqrt(mean((x - y)^2))
    )
  }))
}

# --- load ---
pairs <- as.data.table(checkpoint$bridge_df_for_eval)
common_features <- checkpoint$common_eval_genes
stopifnot(nrow(pairs) > 1L, !anyDuplicated(pairs$original_sample),
          !anyDuplicated(pairs$batch9_replicate), length(common_features) > 1L)

# --- analysis 1: bridge agreement ---
pair_metrics <- rbindlist(lapply(names(method_files), function(method) {
  mat <- readRDS(file.path(batch_dir, method_files[[method]]))
  stopifnot(all(common_features %in% rownames(mat)), all(unlist(pairs) %in% colnames(mat)))
  measure_pairs(mat[common_features, , drop = FALSE], pairs, method)
}))
method_summary <- pair_metrics[, .(
  n_pairs = .N,
  mean_pearson = mean(pearson),
  mean_ccc = mean(ccc),
  mean_mae = mean(mae),
  mean_rmse = mean(rmse)
), by = method]

combat <- copy(pair_metrics[method == "ComBat_REF"])
limma_result <- copy(pair_metrics[method == "limma"])
setnames(combat, c("pearson", "ccc", "mae", "rmse"), paste0(c("pearson", "ccc", "mae", "rmse"), "_ComBat_REF"))
setnames(limma_result, c("pearson", "ccc", "mae", "rmse"), paste0(c("pearson", "ccc", "mae", "rmse"), "_limma"))

comparison <- merge(
  combat[, setdiff(names(combat), "method"), with = FALSE],
  limma_result[, setdiff(names(limma_result), "method"), with = FALSE],
  by = c("pair_id", "original_sample", "batch9_replicate")
)
comparison[, `:=`(
  delta_pearson = pearson_ComBat_REF - pearson_limma,
  delta_ccc = ccc_ComBat_REF - ccc_limma,
  delta_mae = mae_limma - mae_ComBat_REF,
  delta_rmse = rmse_limma - rmse_ComBat_REF
)]
delta_columns <- c("delta_pearson", "delta_ccc", "delta_mae", "delta_rmse")
stopifnot(nrow(comparison) == nrow(pairs),
          all(vapply(comparison[, ..delta_columns], function(x) all(is.finite(x)), logical(1))))

# --- analysis 2: donor-paired bootstrap ---
set.seed(bootstrap_seed)
bootstrap_index <- replicate(bootstrap_iterations, sample.int(nrow(comparison), nrow(comparison), replace = TRUE))
bootstrap_values <- vapply(seq_len(bootstrap_iterations),
                           function(i) colMeans(comparison[bootstrap_index[, i], ..delta_columns]),
                           numeric(length(delta_columns)))

bootstrap_summary <- rbindlist(lapply(seq_along(delta_columns), function(i) {
  values <- bootstrap_values[i, ]
  data.table(
    estimand = delta_columns[i],
    direction_favoring_ComBat_REF = "> 0",
    observed_difference = mean(comparison[[delta_columns[i]]]),
    observed_median_difference = median(comparison[[delta_columns[i]]]),
    ci95_lower = quantile(values, 0.025, names = FALSE),
    ci95_upper = quantile(values, 0.975, names = FALSE),
    bootstrap_probability_favoring_ComBat_REF = mean(values > 0),
    n_pairs_favoring_ComBat_REF = sum(comparison[[delta_columns[i]]] > 0),
    n_pairs = nrow(comparison),
    bootstrap_iterations = bootstrap_iterations
  )
}))

# --- save ---
fwrite(pair_metrics, file.path(batch_dir, "bridge_agreement_pair_metrics.tsv"), sep = "\t")
fwrite(method_summary, file.path(batch_dir, "bridge_agreement_method_summary.tsv"), sep = "\t")
fwrite(comparison, file.path(batch_dir, "bridge_agreement_paired_comparison.tsv"), sep = "\t")
fwrite(bootstrap_summary, file.path(batch_dir, "bridge_agreement_paired_bootstrap.tsv"), sep = "\t")

print(method_summary)
print(bootstrap_summary)
