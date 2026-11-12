#!/usr/bin/env Rscript
# Bridge-pair agreement of batch-correction methods, with a donor-paired bootstrap of ComBat_REF vs limma.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.local.sh)")

suppressPackageStartupMessages({
  library(data.table)
})

# --- config ---

data_dir <- file.path(Sys.getenv("APA_ROOT"), "data/01_quant")
checkpoint_file <- file.path(data_dir, "batch_correction_checkpoint.rds")
matrix_files <- c(
  Raw = file.path(data_dir, "method_matrix_Raw.rds"),
  ComBat_RAW = file.path(data_dir, "method_matrix_ComBat_RAW.rds"),
  ComBat_REF = file.path(data_dir, "method_matrix_ComBat_REF.rds"),
  ComBat_REF_NP = file.path(data_dir, "method_matrix_ComBat_REF_NP.rds"),
  ComBat_REF_meanonly = file.path(data_dir, "method_matrix_ComBat_REF_meanonly.rds"),
  limma = file.path(data_dir, "method_matrix_limma.rds"),
  `RUV-III` = file.path(data_dir, "method_matrix_RUV-III.rds"),
  RUVs = file.path(data_dir, "method_matrix_RUVs.rds")
)
bootstrap_iterations <- 10000L
bootstrap_seed <- 20261012L

# --- loading data ---

checkpoint <- readRDS(checkpoint_file)
pairs <- as.data.table(checkpoint$bridge_df_for_eval)
common_features <- checkpoint$common_eval_genes

stopifnot(nrow(pairs) > 1L, !anyDuplicated(pairs$original_sample),
          !anyDuplicated(pairs$batch9_replicate), length(common_features) > 1L)

# --- analysis 1: bridge agreement ---

pair_metrics <- rbindlist(lapply(names(matrix_files), function(method) {
  mat <- readRDS(matrix_files[[method]])
  stopifnot(all(common_features %in% rownames(mat)),
            all(unlist(pairs) %in% colnames(mat)))
  mat <- mat[common_features, , drop = FALSE]
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    x <- mat[, pairs$original_sample[i]]
    y <- mat[, pairs$batch9_replicate[i]]
    keep <- is.finite(x) & is.finite(y)
    xk <- x[keep]
    yk <- y[keep]
    data.table(
      pair_id = i,
      original_sample = pairs$original_sample[i],
      batch9_replicate = pairs$batch9_replicate[i],
      method = method,
      pearson = cor(x, y),
      ccc = 2 * cov(xk, yk) / (var(xk) + var(yk) + (mean(xk) - mean(yk))^2),
      mae = mean(abs(x - y)),
      rmse = sqrt(mean((x - y)^2))
    )
  }))
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
setnames(combat, c("pearson", "ccc", "mae", "rmse"),
         paste0(c("pearson", "ccc", "mae", "rmse"), "_ComBat_REF"))
setnames(limma_result, c("pearson", "ccc", "mae", "rmse"),
         paste0(c("pearson", "ccc", "mae", "rmse"), "_limma"))

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
stopifnot(
  nrow(comparison) == nrow(pairs),
  all(vapply(comparison[, ..delta_columns], function(x) all(is.finite(x)), logical(1)))
)

# --- analysis 2: donor-paired bootstrap ---

set.seed(bootstrap_seed)
bootstrap_index <- replicate(
  bootstrap_iterations,
  sample.int(nrow(comparison), nrow(comparison), replace = TRUE)
)
bootstrap_values <- vapply(
  seq_len(bootstrap_iterations),
  function(i) colMeans(comparison[bootstrap_index[, i], ..delta_columns]),
  numeric(length(delta_columns))
)

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

# --- save outputs ---

fwrite(pair_metrics, file.path(data_dir, "bridge_agreement_pair_metrics.tsv"), sep = "\t")
fwrite(method_summary, file.path(data_dir, "bridge_agreement_method_summary.tsv"), sep = "\t")
fwrite(comparison, file.path(data_dir, "bridge_agreement_paired_comparison.tsv"), sep = "\t")
fwrite(bootstrap_summary, file.path(data_dir, "bridge_agreement_paired_bootstrap.tsv"), sep = "\t")

print(method_summary)
print(bootstrap_summary)
