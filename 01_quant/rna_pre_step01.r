#!/usr/bin/env Rscript
# RNA-seq preprocessing: alias fix, sex check, filtering, TMM/VST, outliers, ComBat, INT, QTL covariates.
# Usage: Rscript rna_pre_step01.r [gene|transcript]  (default gene)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(WGCNA)
  library(tibble)
  library(dplyr)
  library(sva)
  library(limma)
  library(DESeq2)
  library(PCAForQTL)
})

# --- settings ---
root <- Sys.getenv("APA_ROOT")
feature <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(feature)) feature <- "gene"
count_file <- file.path(root, "raw", if (feature == "gene") "all.912.expr.count" else "all.912.isoform.count")
tpm_file <- file.path(root, "raw", if (feature == "gene") "all.912.expr.tpm" else "all.912.isoform.tpm")
meta_file <- file.path(root, "raw", "all.912.meta")
gtf_file <- file.path(root, "ref", "gencode.v40.annotation.gtf.gz")
out_dir <- file.path(root, "data", "01_quant", feature)
combat_file <- file.path(root, "data", "01_quant", if (feature == "gene") "expr_after_combat.tsv" else "isoform_after_combat.tsv")
geno_file <- ""  # optional genotype PCs; not set
sex_col <- "sex"
batch_col <- "batch"
stage_col <- "LifeStage"
age_col <- "age"
ref_batch <- "batch9"
adjust_col <- ""
tpm_threshold <- 0.1
count_threshold <- 6
sample_frac <- 0.20
outlier_threshold <- 3
set.seed(2026)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- helpers ---
canonical_sample_ids <- function(ids) {
  duplicate_20171125 <- which(ids == "XYA20171125")
  if (length(duplicate_20171125) > 2) stop("More than two XYA20171125 columns were found")
  if (length(duplicate_20171125) == 2) ids[duplicate_20171125[2]] <- "20171125C"
  aliases <- c(
    XYA20180516B = "XYA20180516",
    XYA20210813B = "XYA20210813",
    XYA20230626B = "XYA20230626",
    XYA20231026B = "XYA20231026",
    XYA20240412A = "XYA20240412"
  )
  matched <- match(ids, names(aliases))
  ids[!is.na(matched)] <- aliases[matched[!is.na(matched)]]
  ids
}

save_data <- function(data, filepath, id_col = "ID") {
  df <- as.data.frame(data)
  if (is.character(rownames(df))) {
    df <- tibble::rownames_to_column(df, id_col)
  }
  fwrite(df, filepath, sep = "\t", row.names = FALSE, quote = FALSE)
}

# ComBat (mean-only, LifeStage protected), optional covariate removal, INT
process_path <- function(norm_matrix, meta_data, method_name) {
  clean_matrix <- norm_matrix[, rownames(meta_data), drop = FALSE]
  save_data(clean_matrix, file.path(out_dir, paste0("norm_clean_", method_name, ".tsv")))

  combat_matrix <- as.matrix(clean_matrix)
  if (batch_col %in% colnames(meta_data)) {
    batch <- factor(meta_data[[batch_col]])
    mod <- model.matrix(reformulate(stage_col), data = meta_data)
    ref_batch_use <- if (ref_batch %in% levels(batch)) ref_batch else NULL
    # Mean-only ComBat keeps within-batch residual variation when batch and LifeStage are aliased.
    combat_matrix <- ComBat(dat = combat_matrix, batch = batch, mod = mod, ref.batch = ref_batch_use,
                            par.prior = TRUE, mean.only = TRUE, prior.plots = FALSE)
  }
  save_data(combat_matrix, file.path(out_dir, paste0("combat_", method_name, ".tsv")))
  if (method_name == "tmm") save_data(combat_matrix, combat_file, feature_id_col)

  adjust_cols <- trimws(unlist(strsplit(adjust_col, ",")))
  adjust_cols <- adjust_cols[adjust_cols %in% colnames(meta_data)]
  final_matrix <- if (length(adjust_cols) > 0) {
    covariates <- model.matrix(as.formula(paste("~ 0 +", paste(adjust_cols, collapse = " + "))), data = meta_data)
    removeBatchEffect(combat_matrix, covariates = covariates)
  } else {
    combat_matrix
  }
  # inverse-normal transform, per feature
  int_matrix <- t(apply(as.matrix(final_matrix), 1, function(x) qnorm(rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1))))
  save_data(int_matrix, file.path(out_dir, paste0(method_name, ".final.int.tsv")))
  int_matrix
}

# --- load ---
count_header <- colnames(fread(count_file, nrows = 0, check.names = FALSE))
tpm_header <- colnames(fread(tpm_file, nrows = 0, check.names = FALSE))
feature_id_col <- count_header[[1]]
id_col <- if (feature == "gene") "gene_id" else "transcript_id"
if (feature_id_col != id_col || tpm_header[[1]] != id_col) {
  stop("Count/TPM ID columns do not match feature type ", feature)
}
count_header[-1] <- canonical_sample_ids(count_header[-1])
tpm_header[-1] <- canonical_sample_ids(tpm_header[-1])
if (anyDuplicated(count_header[-1]) || anyDuplicated(tpm_header[-1])) {
  stop("Count/TPM matrix contains duplicated sample IDs after alias normalization")
}
if (!identical(count_header, tpm_header)) stop("Count and TPM matrices must have identical sample order")
count <- fread(count_file, data.table = FALSE, check.names = FALSE)
tpm <- fread(tpm_file, data.table = FALSE, check.names = FALSE)
meta <- fread(meta_file, data.table = FALSE, check.names = FALSE)
colnames(count) <- count_header
colnames(tpm) <- tpm_header
if (!identical(count[[1]], tpm[[1]])) stop("Count and TPM matrices must have identical feature order")
if (anyDuplicated(meta[[1]])) stop("Metadata contains duplicated sample IDs")
rownames(count) <- count[[1]]; count[[1]] <- NULL
rownames(tpm) <- tpm[[1]]; tpm[[1]] <- NULL
rownames(meta) <- meta[[1]]; meta[[1]] <- NULL
if (!stage_col %in% colnames(meta)) {
  if (!age_col %in% colnames(meta)) stop("Metadata must contain LifeStage or numeric age")
  meta[[stage_col]] <- ifelse(meta[[age_col]] < 0, "Prenatal", "Postnatal")
}
meta[[stage_col]] <- factor(meta[[stage_col]], levels = c("Prenatal", "Postnatal"))
if (anyNA(meta[[stage_col]])) stop("LifeStage contains values outside Prenatal/Postnatal")
meta <- meta[meta$RIN>=5,]
common_samples <- intersect(rownames(meta), colnames(count))
meta <- meta[common_samples, , drop = FALSE]
count <- count[, common_samples]
tpm <- tpm[, common_samples]

# --- sex check ---
meta$sex_imputed <- if (sex_col %in% colnames(meta)) {
  reported <- tolower(trimws(as.character(meta[[sex_col]])))
  reported[reported %in% c("1", "female", "f")] <- "female"
  reported[reported %in% c("0", "male", "m")] <- "male"
  reported[reported %in% c("unknown", "", "na", "n/a")] <- NA_character_
  reported
} else {
  rep(NA_character_, nrow(meta))
}

gtf <- as.data.frame(rtracklayer::import(gtf_file))
annotation_ids <- sub("\\.[0-9]+$", "", as.character(gtf[[id_col]]))
feature_ids <- sub("\\.[0-9]+$", "", rownames(count))
xist_ids <- rownames(count)[feature_ids %in% unique(annotation_ids[gtf$gene_name == "XIST"])]
uty_ids <- rownames(count)[feature_ids %in% unique(annotation_ids[gtf$gene_name == "UTY"])]
if (!(length(xist_ids) && length(uty_ids))) stop("XIST or UTY features were not found in the count matrix")

xist_expr <- colSums(as.matrix(count[xist_ids, , drop = FALSE]))
uty_expr <- colSums(as.matrix(count[uty_ids, , drop = FALSE]))
km_xist <- kmeans(log1p(xist_expr), centers = 2, nstart = 25)
km_uty <- kmeans(log1p(uty_expr), centers = 2, nstart = 25)
xist_threshold <- expm1(mean(km_xist$centers))
uty_threshold <- expm1(mean(km_uty$centers))
cat(sprintf("XIST k-means threshold: %.2f; UTY k-means threshold: %.2f\n", xist_threshold, uty_threshold))

xist_pos <- xist_expr > xist_threshold
uty_pos <- uty_expr > uty_threshold
imputed_sex <- rep(NA_character_, length(xist_pos))
names(imputed_sex) <- names(xist_pos)
imputed_sex[xist_pos & !uty_pos] <- "female"
imputed_sex[uty_pos] <- "male"
is_ambiguous <- !xist_pos & !uty_pos
is_unknown <- is.na(meta$sex_imputed)
imputed_sex <- imputed_sex[rownames(meta)]
meta$sex_imputed[is_unknown] <- imputed_sex[is_unknown]

sex_mismatch <- rownames(meta)[!is_unknown & !is_ambiguous &
                                !is.na(imputed_sex) &
                                meta$sex_imputed != imputed_sex]
ambiguous_samples <- rownames(meta)[is_ambiguous]
samples_to_filter <- unique(c(sex_mismatch, ambiguous_samples))

sex_df <- data.frame(
  Sample = rownames(meta), XIST = xist_expr, UTY = uty_expr,
  Sex = ifelse(is.na(imputed_sex), "Ambiguous", imputed_sex),
  Ambiguous = rownames(meta) %in% ambiguous_samples,
  stringsAsFactors = FALSE
)
save_data(sex_df, file.path(out_dir, "sex_check_plot_data.tsv"))
cat("Total samples:", nrow(meta), "\n")
cat("Sex mismatches:", length(sex_mismatch), "\n")
cat("Ambiguous samples:", length(ambiguous_samples), "\n")

if (length(samples_to_filter) > 0) {
  keep <- setdiff(rownames(meta), samples_to_filter)
  meta <- meta[keep, , drop = FALSE]
  count <- count[, keep, drop = FALSE]
  tpm <- tpm[, keep, drop = FALSE]
}

meta$sex <- meta$sex_imputed
meta$sex_imputed <- NULL
save_data(meta, file.path(out_dir, "meta_after_sexcheck.tsv"))

# --- gene filtering ---
autosomal_features <- unique(annotation_ids[
  as.character(gtf$seqnames) %in% paste0("chr", 1:22) & !is.na(annotation_ids)
])
threshold <- sample_frac * ncol(count)
tpm_filter <- rowSums(tpm > tpm_threshold, na.rm = TRUE) >= threshold
count_filter <- rowSums(count >= count_threshold, na.rm = TRUE) >= threshold
keep_features <- tpm_filter & count_filter & (feature_ids %in% autosomal_features)
count_filtered <- count[keep_features, ]
if (!nrow(count_filtered)) stop("No features retained after expression/autosome filtering")
cat("Retained", nrow(count_filtered), feature, "features\n")
save_data(count_filtered, file.path(out_dir, "count_filtered.tsv"))

# --- normalization and outlier detection ---
count_tmm <- cpm(calcNormFactors(DGEList(count_filtered)), log = TRUE)
save_data(count_tmm, file.path(out_dir, "norm_tmm_before_outlier_removal.tsv"))
count_vst <- DESeq2::varianceStabilizingTransformation(as.matrix(round(count_filtered)), blind = TRUE)
save_data(count_vst, file.path(out_dir, "norm_vst_before_outlier_removal.tsv"))

outlier_list <- list()
for (stage in levels(meta[[stage_col]])) {
  stage_samples <- rownames(meta)[meta[[stage_col]] == stage]
  if (length(stage_samples) < 10) next
  for (norm_name in c("tmm", "vst")) {
    norm_matrix <- if (norm_name == "tmm") count_tmm else count_vst
    stage_matrix <- norm_matrix[, stage_samples, drop = FALSE]
    top_genes <- order(apply(stage_matrix, 1, var), decreasing = TRUE)[seq_len(min(5000, nrow(stage_matrix)))]
    z_k <- scale(fundamentalNetworkConcepts(adjacency(t(as.matrix(stage_matrix[top_genes, , drop = FALSE])), type = "signed", corFnc = "bicor"))$Connectivity)
    outlier_list[[paste(norm_name, stage, sep = "_")]] <- stage_samples[is.finite(z_k[, 1]) & z_k[, 1] < -outlier_threshold]
  }
}
detected_outliers <- unique(unlist(outlier_list))
save_data(data.frame(Sample = detected_outliers), file.path(out_dir, "outlier_samples.txt"))

meta_clean <- if (length(detected_outliers) > 0) meta[!rownames(meta) %in% detected_outliers, ] else meta
save_data(meta_clean, file.path(out_dir, "metadata_final.tsv"))

# --- batch adjustment and INT ---
tmm_int <- process_path(count_tmm, meta_clean, "tmm")
vst_int <- process_path(count_vst, meta_clean, "vst")

# --- PCA covariates for QTL (gene only) ---
if (feature == "gene") {
  expr_for_qtl <- tmm_int
  pca_result <- prcomp(t(expr_for_qtl), center = TRUE, scale. = TRUE)
  sdev_sq <- (pca_result$sdev)^2
  PVEs <- sdev_sq / sum(sdev_sq)
  cumPVEs <- cumsum(PVEs)

  K_elbow <- tryCatch({
    k_e <- PCAForQTL::runElbow(prcompResult = pca_result)
    if (is.numeric(k_e) && length(k_e) == 1 && !is.na(k_e)) as.integer(k_e) else NA_integer_
  }, error = function(e) {
    NA_integer_
  })

  K_BE <- tryCatch({
    set.seed(2025)
    be_res <- PCAForQTL::runBE(t(expr_for_qtl), B = 20)
    if (!is.null(be_res$numOfPCsChosen) && is.numeric(be_res$numOfPCsChosen)) as.integer(be_res$numOfPCsChosen) else NA_integer_
  }, error = function(e) {
    NA_integer_
  })

  k_cum <- NA_integer_
  if (is.na(K_elbow) && is.na(K_BE) && any(!is.na(cumPVEs))) {
    k_idx <- which(cumPVEs >= 0.95)[1]
    if (!is.na(k_idx) && length(k_idx) == 1) k_cum <- as.integer(k_idx)
  }

  available_pcs <- ncol(pca_result$x)
  capK <- min(30L, available_pcs)
  chosen_values <- na.omit(c(K_elbow, K_BE))
  if (length(chosen_values) == 2) {
    K_candidate <- max(chosen_values)
  } else if (length(chosen_values) == 1) {
    K_candidate <- chosen_values[1]
  } else if (!is.na(k_cum)) {
    K_candidate <- k_cum
  } else {
    K_candidate <- min(10L, available_pcs)
  }
  K_final <- as.integer(max(1L, min(K_candidate, capK)))
  cat("K_elbow=", K_elbow, " K_BE=", K_BE, " k_cum(>=0.95)=", k_cum, " => using K_final=", K_final, "\n")

  save_data(pca_result$x, file.path(out_dir, "pca_results_for_qtl.tsv"))
  diag_df <- data.frame(
    method = c("K_elbow", "K_BE", "k_cum_0.95"),
    K_value = c(
      if (is.na(K_elbow)) NA_integer_ else as.integer(K_elbow),
      if (is.na(K_BE)) NA_integer_ else as.integer(K_BE),
      if (is.na(k_cum)) NA_integer_ else as.integer(k_cum)
    )
  )
  write.table(diag_df, file = file.path(out_dir, "pca_for_qtl_k_values.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  known_covariates <- meta_clean[colnames(expr_for_qtl), , drop = FALSE]
  cols_to_remove <- c(batch_col, trimws(unlist(strsplit(adjust_col, ","))))
  cols_to_remove <- cols_to_remove[cols_to_remove %in% colnames(known_covariates)]
  if (length(cols_to_remove) > 0) {
    known_covariates <- known_covariates[, !colnames(known_covariates) %in% cols_to_remove, drop = FALSE]
  }
  known_covariates <- known_covariates %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(across(where(is.factor), as.numeric))
  known_covariates <- known_covariates[, sapply(known_covariates, is.numeric), drop = FALSE]

  if (geno_file != "" && file.exists(geno_file)) {
    geno <- fread(geno_file, data.table = FALSE)
    rownames(geno) <- geno[, 1]
    geno_cols <- if (ncol(geno) >= 3) 3:min(7, ncol(geno)) else integer(0)
    if (length(geno_cols) >= 1) {
      geno_sub <- geno[, geno_cols, drop = FALSE]
      colnames(geno_sub) <- paste0("genoPC", seq_len(ncol(geno_sub)))
      common_geno_samples <- intersect(rownames(geno_sub), rownames(known_covariates))
      if (length(common_geno_samples) > 0) {
        known_covariates <- known_covariates[common_geno_samples, , drop = FALSE]
        known_covariates <- cbind(known_covariates, geno_sub[common_geno_samples, , drop = FALSE])
      }
    }
  }

  K_use <- min(K_final, ncol(pca_result$x))
  PCsTop <- as.data.frame(pca_result$x)[, seq_len(K_use), drop = FALSE]
  common_samples <- intersect(rownames(known_covariates), rownames(PCsTop))
  known_covariates <- known_covariates[common_samples, , drop = FALSE]
  PCsTop <- PCsTop[common_samples, , drop = FALSE]
  if (ncol(known_covariates) > 0) {
    knownCovariatesFiltered <- tryCatch(
      PCAForQTL::filterKnownCovariates(known_covariates, PCsTop, unadjustedR2_cutoff = 0.9),
      error = function(e) {
        warning("filterKnownCovariates failed: ", conditionMessage(e))
        known_covariates
      }
    )
  } else {
    knownCovariatesFiltered <- known_covariates
  }
  stopifnot(identical(rownames(knownCovariatesFiltered), rownames(PCsTop)))
  covariatesToUse <- cbind(knownCovariatesFiltered, PCsTop)
  save_data(t(covariatesToUse), file.path(out_dir, "covariates_for_qtl.txt"))
  final_expr_no_version <- expr_for_qtl[, common_samples, drop = FALSE]
  rownames(final_expr_no_version) <- sub("\\..*$", "", rownames(final_expr_no_version))
  save_data(final_expr_no_version, file.path(out_dir, "expression_for_qtl.tsv"))
}
