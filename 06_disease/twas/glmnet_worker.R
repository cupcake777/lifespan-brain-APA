# Per-phenotype glmnet elastic-net training with CV performance, one model per feature.
suppressPackageStartupMessages(library(glmnet))

# --- input ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: glmnet_worker.R <x.tsv> <metadata.tsv> <phenotype.tsv> <weights.tsv> <detail.tsv> <n_folds> <seed>")
}

x_path <- args[[1]]
metadata_path <- args[[2]]
phenotype_path <- args[[3]]
weights_out <- args[[4]]
detail_out <- args[[5]]
n_folds <- as.integer(args[[6]])
seed_base <- as.integer(args[[7]])

min_abs_weight <- 1e-12
pmax_mode <- "control=list(pmax=200)"

x_df <- read.table(x_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.table(metadata_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- read.table(phenotype_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

if (!("sample_id" %in% colnames(x_df))) {
  stop("x.tsv must include sample_id")
}
sample_ids <- as.character(x_df$sample_id)
x_mat <- as.matrix(x_df[, setdiff(colnames(x_df), "sample_id"), drop = FALSE])
storage.mode(x_mat) <- "double"

weight_rows <- list()
detail_rows <- list()

# --- fit per phenotype ---
for (i in seq_len(nrow(phenotype))) {
  model_id <- as.character(phenotype$model_id[i])
  feature_index <- as.integer(phenotype$feature_index[i])
  y <- as.numeric(phenotype[i, sample_ids])
  finite_y <- is.finite(y)
  n_training <- sum(finite_y)
  status <- "not_started"
  n_nonzero <- 0L
  pred_r <- NA_real_
  pred_r2 <- NA_real_
  pred_p <- NA_real_
  lambda_min <- NA_real_
  lambda_min_idx <- NA_integer_
  n_complete_cv <- NA_integer_
  na_fit_preval <- NA_integer_
  warning_text <- ""
  elapsed <- NA_real_

  if (ncol(x_mat) < 2) {
    status <- "too_few_snps"
  } else if (n_training < max(n_folds + 1L, 10L)) {
    status <- "too_few_samples"
  } else if (stats::sd(y[finite_y]) == 0 || length(unique(y[finite_y])) < 2) {
    status <- "constant_phenotype"
  } else {
    y_used <- y[finite_y]
    x_used <- x_mat[finite_y, , drop = FALSE]
    set.seed(seed_base + feature_index)
    foldid <- sample(rep(seq_len(n_folds), length.out = length(y_used)))
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      withCallingHandlers(
        cv.glmnet(
          x = x_used,
          y = y_used,
          nfolds = n_folds,
          foldid = foldid,
          keep = TRUE,
          alpha = 0.5,
          nlambda = 50,
          control = list(pmax = 200)
        ),
        warning = function(w) {
          warning_text <<- paste(c(warning_text, conditionMessage(w)), collapse = " | ")
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) e
    )
    elapsed <- proc.time()[["elapsed"]] - started

    if (inherits(fit, "error")) {
      status <- paste0("fit_failed:", class(fit)[[1]])
      warning_text <- conditionMessage(fit)
    } else {
      lambda_min <- fit$lambda.min
      lambda_min_idx <- as.integer(fit$index["min", ])
      beta <- as.numeric(fit$glmnet.fit$beta[, lambda_min_idx])
      names(beta) <- colnames(x_mat)
      cv_pred <- fit$fit.preval[, lambda_min_idx]
      na_fit_preval <- sum(!is.finite(cv_pred))
      complete <- complete.cases(y_used, cv_pred)
      n_complete_cv <- sum(complete)
      if (n_complete_cv >= 3 && stats::sd(cv_pred[complete]) > 0 && stats::sd(y_used[complete]) > 0) {
        performance <- suppressWarnings(stats::cor.test(y_used[complete], cv_pred[complete]))
        pred_r <- as.numeric(performance$estimate)
        pred_r2 <- pred_r ^ 2
        pred_p <- as.numeric(performance$p.value)
      }

      nonzero <- which(abs(beta) > min_abs_weight)
      n_nonzero <- length(nonzero)
      if (n_nonzero == 0L) {
        status <- "trained_zero"
      } else if (is.finite(pred_r) && is.finite(pred_p) && pred_r > 0.1 && pred_p < 0.05) {
        status <- "written"
        for (idx in nonzero) {
          raw_col <- names(beta)[[idx]]
          meta_idx <- match(raw_col, metadata$raw_col)
          weight_rows[[length(weight_rows) + 1L]] <- data.frame(
            model_id = model_id,
            raw_col = raw_col,
            rsid = metadata$rsid[meta_idx],
            gene = model_id,
            weight = beta[[idx]],
            ref_allele = metadata$ref_allele[meta_idx],
            eff_allele = metadata$eff_allele[meta_idx],
            stringsAsFactors = FALSE
          )
        }
      } else {
        status <- "trained_nonzero_failed_perf"
      }
    }
  }

  detail_rows[[length(detail_rows) + 1L]] <- data.frame(
    model_id = model_id,
    feature_index = feature_index,
    status = status,
    n_training_samples = n_training,
    n_snps_available = ncol(x_mat),
    n_nonzero_weights = n_nonzero,
    pred_perf_r = pred_r,
    pred_perf_R2 = pred_r2,
    pred_perf_pval = pred_p,
    lambda_min = lambda_min,
    lambda_min_idx = lambda_min_idx,
    cv_complete_observations = n_complete_cv,
    fit_preval_na = na_fit_preval,
    glmnet_seconds = elapsed,
    pmax_mode = pmax_mode,
    glmnet_warnings = warning_text,
    stringsAsFactors = FALSE
  )
}

# --- write outputs ---
weights_df <- if (length(weight_rows)) {
  do.call(rbind, weight_rows)
} else {
  data.frame(model_id = character(), raw_col = character(), rsid = character(), gene = character(), weight = numeric(), ref_allele = character(), eff_allele = character())
}
detail_df <- if (length(detail_rows)) {
  do.call(rbind, detail_rows)
} else {
  data.frame()
}

write.table(weights_df, weights_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(detail_df, detail_out, sep = "\t", quote = FALSE, row.names = FALSE)
