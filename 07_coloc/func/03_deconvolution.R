#!/usr/bin/env Rscript
# BisqueRNA cell-type deconvolution of bulk counts; stage-specific ILR covariates.
# Tests whether cell composition is captured by the APA-QTL covariates and residual QTL signal.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.local.sh)")

suppressPackageStartupMessages({
  library(arrow)
  library(BisqueRNA)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(limSolve)
})

# --- helpers ---

safe_cell_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

# --- config ---

project_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
bulk_count_path <- file.path(Sys.getenv("APA_ROOT"), "raw/all.912.count")
metadata_path <- file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/CHB.all.meta")
prenatal_samples_path <- file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/sample.Prenatal.final")
postnatal_samples_path <- file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/sample.Postnatal.final")
prenatal_covariates_path <- file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/covariates_for_qtl.Prenatal.txt")
postnatal_covariates_path <- file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/covariates_for_qtl.Postnatal.txt")

reference_dir <- file.path(project_dir, "reference")
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
pseudocount <- 1e-6

stage_config <- list(
  prenatal = list(
    samples = scan(prenatal_samples_path, what = character(), quiet = TRUE),
    signature = file.path(reference_dir, "cibersortx_signature_prenatal.tsv"),
    pseudobulk = file.path(reference_dir, "pseudobulk_markers_prenatal.parquet"),
    covariates = prenatal_covariates_path
  ),
  postnatal = list(
    samples = scan(postnatal_samples_path, what = character(), quiet = TRUE),
    signature = file.path(reference_dir, "cibersortx_signature_postnatal.tsv"),
    pseudobulk = file.path(reference_dir, "pseudobulk_markers_postnatal.parquet"),
    covariates = postnatal_covariates_path
  )
)

# --- loading bulk counts and references ---

all_samples <- unlist(lapply(stage_config, `[[`, "samples"), use.names = FALSE)
bulk <- fread(bulk_count_path, select = c("gene_id", all_samples))
bulk[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
bulk_library_sizes <- vapply(bulk[, ..all_samples], sum, numeric(1))

marker_union <- unique(unlist(lapply(stage_config, function(settings) fread(settings$signature)$NAME)))
bulk <- bulk[gene_id %chin% marker_union]

# --- stage-specific BisqueRNA deconvolution ---

stage_results <- list()
for (stage in names(stage_config)) {
  settings <- stage_config[[stage]]
  signature <- fread(settings$signature)
  pseudobulk <- as.data.table(read_parquet(settings$pseudobulk))
  cell_classes <- setdiff(names(signature), "NAME")
  samples <- settings$samples
  genes <- intersect(signature$NAME, bulk$gene_id)

  reference <- as.matrix(signature[match(genes, NAME), ..cell_classes])
  rownames(reference) <- genes

  donor_counts <- dcast(
    pseudobulk,
    donor ~ cell_class,
    value.var = "n_cells",
    fun.aggregate = sum,
    fill = 0
  )
  donor_names <- donor_counts$donor
  donor_composition <- t(as.matrix(donor_counts[, ..cell_classes]))
  donor_composition <- sweep(donor_composition, 2, colSums(donor_composition), "/")
  rownames(donor_composition) <- cell_classes
  colnames(donor_composition) <- donor_names

  marker_counts <- as.matrix(bulk[match(genes, gene_id), ..samples])
  rownames(marker_counts) <- genes
  bulk_cpm <- sweep(marker_counts, 2, bulk_library_sizes[samples], "/") * 1e6

  synthetic_reference <- reference %*% donor_composition
  sample_template <- setNames(numeric(length(samples)), samples)
  transformed <- vapply(
    genes,
    function(gene) {
      BisqueRNA:::SemisupervisedTransformBulk(
        gene,
        Y.train = synthetic_reference,
        X.pred = bulk_cpm
      )[, 1]
    },
    FUN.VALUE = sample_template
  )
  transformed <- t(transformed)
  rownames(transformed) <- genes
  usable <- apply(transformed, 1, function(x) all(is.finite(x)))
  transformed <- transformed[usable, , drop = FALSE]
  reference <- reference[usable, , drop = FALSE]

  n_classes <- length(cell_classes)
  equality_matrix <- matrix(1, nrow = n_classes, ncol = n_classes)
  equality_target <- rep(1, n_classes)
  nonnegative_matrix <- diag(n_classes)
  nonnegative_target <- rep(0, n_classes)
  solutions <- vapply(
    seq_len(ncol(transformed)),
    function(sample_index) {
      fit <- lsei(
        A = reference,
        B = transformed[, sample_index],
        E = equality_matrix,
        F = equality_target,
        G = nonnegative_matrix,
        H = nonnegative_target,
        verbose = FALSE
      )
      c(fit$X, rnorm = sqrt(fit$solutionNorm))
    },
    FUN.VALUE = numeric(n_classes + 1)
  )
  rownames(solutions) <- c(cell_classes, "bisque_rnorm")
  colnames(solutions) <- samples

  fitted <- reference %*% solutions[cell_classes, , drop = FALSE]
  fit_correlation <- vapply(
    seq_len(ncol(transformed)),
    function(i) cor(transformed[, i], fitted[, i]),
    numeric(1)
  )
  relative_residual <- solutions["bisque_rnorm", ] / sqrt(colSums(transformed^2))

  fractions <- as.data.table(t(solutions[cell_classes, , drop = FALSE]), keep.rownames = "sample")
  setnames(fractions, cell_classes, safe_cell_name(cell_classes))
  fractions[, stage := stage]
  fractions[, bisque_rnorm := solutions["bisque_rnorm", ]]
  fractions[, relative_residual_norm := relative_residual]
  fractions[, fit_correlation := fit_correlation]
  fractions[, n_genes_used := nrow(reference)]

  qc <- data.table(
    stage = stage,
    n_samples = nrow(fractions),
    n_reference_donors = ncol(donor_composition),
    n_cell_classes = n_classes,
    n_genes_used = nrow(reference),
    median_fit_correlation = median(fit_correlation),
    min_fit_correlation = min(fit_correlation),
    median_relative_residual_norm = median(relative_residual),
    max_fraction_sum_error = max(abs(colSums(solutions[cell_classes, , drop = FALSE]) - 1)),
    zero_fraction_rate = mean(solutions[cell_classes, , drop = FALSE] < 1e-8)
  )

  stage_results[[stage]] <- list(
    fractions = fractions,
    qc = qc,
    cell_classes = cell_classes,
    ilr_class_order = safe_cell_name(cell_classes)
  )
}

proportions <- rbindlist(lapply(stage_results, `[[`, "fractions"), fill = TRUE)
fraction_columns <- unique(unlist(lapply(stage_results, `[[`, "ilr_class_order")))
proportions[, fraction_sum := rowSums(.SD, na.rm = TRUE), .SDcols = fraction_columns]

# --- sample metadata and ILR coordinates ---

metadata <- fread(metadata_path)
proportions <- merge(proportions, metadata, by = "sample", all.x = TRUE, sort = FALSE)

for (stage in names(stage_results)) {
  rows <- proportions$stage == stage
  classes <- stage_results[[stage]]$ilr_class_order
  composition <- as.matrix(proportions[rows, ..classes])
  composition <- pmax(composition, pseudocount)
  composition <- composition / rowSums(composition)
  ilr_basis <- contr.helmert(length(classes))
  ilr_basis <- sweep(ilr_basis, 2, sqrt(colSums(ilr_basis^2)), "/")
  ilr <- log(composition) %*% ilr_basis
  for (coordinate in seq_len(ncol(ilr))) {
    proportions[rows, paste0("cellprop_ILR", coordinate) := ilr[, coordinate]]
  }
}

setcolorder(
  proportions,
  c(
    "sample", "stage", fraction_columns, "fraction_sum",
    "bisque_rnorm", "relative_residual_norm", "fit_correlation", "n_genes_used",
    "batch", "RIN", "sex", "age", "LifeStage",
    paste0("cellprop_ILR", 1:6)
  )
)
fwrite(proportions, file.path(project_dir, "cell_proportions.tsv"), sep = "\t")

# --- within-stage age associations ---

association_rows <- list()
row_index <- 1
for (stage_name in names(stage_results)) {
  stage_data <- proportions[proportions$stage == stage_name]
  for (cell_class in stage_results[[stage_name]]$ilr_class_order) {
    stage_data[, fraction_logit := qlogis(pmin(pmax(get(cell_class), pseudocount), 1 - pseudocount))]
    fit <- lm(fraction_logit ~ age + sex + RIN + batch, data = stage_data)
    age_result <- summary(fit)$coefficients["age", ]
    association_rows[[row_index]] <- data.table(
      stage = stage_name,
      cell_class = cell_class,
      n = nobs(fit),
      age_beta = age_result[["Estimate"]],
      age_se = age_result[["Std. Error"]],
      age_p = age_result[["Pr(>|t|)"]]
    )
    row_index <- row_index + 1
  }
}
associations <- rbindlist(association_rows)
associations[, age_fdr_bh := p.adjust(age_p, method = "BH"), by = stage]
fwrite(associations, file.path(project_dir, "cell_proportion_associations.tsv"), sep = "\t")

# --- QTL covariates with cell-composition ILR coordinates ---

covariate_qc <- list()
for (stage in names(stage_config)) {
  original <- read.delim(
    stage_config[[stage]]$covariates,
    row.names = 1,
    check.names = FALSE
  )
  stage_proportions <- proportions[match(colnames(original), sample)]
  ilr_rows <- t(as.matrix(stage_proportions[, paste0("cellprop_ILR", 1:6), with = FALSE]))
  ilr_rows <- sweep(ilr_rows, 1, rowMeans(ilr_rows), "-")
  ilr_rows <- sweep(ilr_rows, 1, apply(ilr_rows, 1, sd), "/")
  rownames(ilr_rows) <- paste0("cellprop_ILR", 1:6)
  colnames(ilr_rows) <- colnames(original)

  combined <- rbind(original, ilr_rows)
  write.table(
    combined,
    file.path(project_dir, if (stage == "prenatal") "covariates_with_cell_proportions.pre.txt" else "covariates_with_cell_proportions.postnatal.txt"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  covariate_qc[[stage]] <- data.table(
    stage = stage,
    n_samples = ncol(combined),
    original_covariates = nrow(original),
    cell_ilr_covariates = nrow(ilr_rows),
    combined_covariates = nrow(combined),
    design_rank = qr(as.matrix(t(combined)))$rank,
    condition_number = kappa(as.matrix(t(combined)))
  )
}
fwrite(
  rbindlist(covariate_qc),
  file.path(project_dir, "cell_proportion_covariate_qc.tsv"),
  sep = "\t"
)

# --- method QC and contract ---

qc <- rbindlist(lapply(stage_results, `[[`, "qc"))
fwrite(qc, file.path(project_dir, "cell_deconvolution_qc.tsv"), sep = "\t")

contract <- list(
  primary_method = "BisqueRNA semisupervised reference-based decomposition",
  primary_status = "complete",
  package = paste0("BisqueRNA_", as.character(packageVersion("BisqueRNA"))),
  bulk_input = "raw gene counts",
  reference_design = "stage-specific donor pseudobulk CPM signatures with observed donor cell-count compositions",
  reference_validation = "leave-one-donor-out marker-selection pseudo-mixture reconstruction",
  prenatal_classes = stage_results$prenatal$cell_classes,
  postnatal_classes = stage_results$postnatal$cell_classes,
  qtl_covariates = "six stage-specific ILR coordinates",
  cibersortx_status = "not_run_unavailable_account_and_environment_not_blocking"
)
write_json(
  contract,
  file.path(project_dir, "cell_deconvolution_method_contract.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# --- covariate capture config ---

cell_proportions_path <- file.path(project_dir, "cell_proportions.tsv")
significant_qtl_path <- file.path(project_dir, "01_significant_qtl_universe.tsv.gz")
figures_dir <- file.path(Sys.getenv("APA_ROOT"), "figures/07_coloc")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
cv_folds <- 10L
cv_seed <- 20260720L

stage_config <- list(
  prenatal = list(
    covariates = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/covariates_for_qtl.Prenatal.txt"),
    phenotypes = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Prenatal.bed.gz"),
    cell_classes = c(
      "Excitatory_neuron", "Inhibitory_neuron", "Astrocyte", "OPC",
      "Microglia_Immune", "Vascular_Mural", "Glial_progenitor"
    ),
    figure_width = 10
  ),
  postnatal = list(
    covariates = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/covariates_for_qtl.Postnatal.txt"),
    phenotypes = file.path(Sys.getenv("APA_ROOT"), "data/04_qtl/pre_input/phenotype.Postnatal.bed.gz"),
    cell_classes = c(
      "Excitatory_neuron", "Inhibitory_neuron", "Astrocyte", "Oligodendrocyte",
      "OPC", "Microglia_Immune", "Vascular_Mural"
    ),
    figure_width = 14
  )
)

correlation_output <- file.path(project_dir, "cell_qtl_covariate_correlations.tsv")
capture_output <- file.path(project_dir, "cell_qtl_covariate_capture.tsv")
residual_output <- file.path(project_dir, "cell_qtl_residual_associations.tsv.gz")
summary_output <- file.path(project_dir, "cell_qtl_covariate_qc_sum.tsv")
significant_hits_output <- file.path(project_dir, "cell_qtl_significant_phenotype_residual_hits.tsv")
layout_output <- file.path(project_dir, "cell_qtl_covariate_correlation.layout_validation.tsv")

# --- load cell composition ---

cell_proportions <- fread(cell_proportions_path)
significant_qtl_phenotypes <- unique(
  fread(significant_qtl_path, select = c("stage", "phenotype_id"))
)
ilr_columns <- paste0("cellprop_ILR", 1:6)

correlation_rows <- list()
capture_rows <- list()
residual_rows <- list()
design_rows <- list()
layout_rows <- list()
correlation_index <- 1L
capture_index <- 1L
residual_index <- 1L
design_index <- 1L
layout_index <- 1L

# --- stage-specific covariate capture and residual tests ---

for (stage_name in names(stage_config)) {
  original <- read.delim(
    stage_config[[stage_name]]$covariates,
    row.names = 1,
    check.names = FALSE
  )
  stage_cells <- cell_proportions[stage == stage_name]
  stage_cells <- stage_cells[match(colnames(original), sample)]
  sample_ids <- colnames(original)
  covariate_matrix <- t(as.matrix(original))
  storage.mode(covariate_matrix) <- "double"
  rownames(covariate_matrix) <- sample_ids

  covariate_groups <- data.table(
    covariate = colnames(covariate_matrix),
    covariate_group = fifelse(
      grepl("^PC[0-9]+$", colnames(covariate_matrix)), "phenotype_PC",
      fifelse(grepl("^genoPC[0-9]+$", colnames(covariate_matrix)), "genotype_PC", "known_covariate")
    )
  )
  known_columns <- covariate_groups[covariate_group == "known_covariate", covariate]
  genotype_pc_columns <- covariate_groups[covariate_group == "genotype_PC", covariate]
  phenotype_pc_columns <- covariate_groups[covariate_group == "phenotype_PC", covariate]

  response_columns <- c(stage_config[[stage_name]]$cell_classes, ilr_columns)
  set.seed(cv_seed + match(stage_name, names(stage_config)))
  fold_id <- sample(rep(seq_len(cv_folds), length.out = nrow(covariate_matrix)))

  for (response in response_columns) {
    response_type <- if (response %in% ilr_columns) "cell_ilr" else "cell_fraction"
    raw_response <- stage_cells[[response]]

    for (covariate_name in colnames(covariate_matrix)) {
      pearson <- cor.test(raw_response, covariate_matrix[, covariate_name], method = "pearson")
      spearman <- suppressWarnings(
        cor.test(raw_response, covariate_matrix[, covariate_name], method = "spearman", exact = FALSE)
      )
      correlation_rows[[correlation_index]] <- data.table(
        stage = stage_name,
        response_type = response_type,
        response = response,
        covariate_group = covariate_groups[covariate == covariate_name, covariate_group],
        covariate = covariate_name,
        n = sum(complete.cases(raw_response, covariate_matrix[, covariate_name])),
        pearson_r = unname(pearson$estimate),
        pearson_p = pearson$p.value,
        spearman_r = unname(spearman$estimate),
        spearman_p = spearman$p.value
      )
      correlation_index <- correlation_index + 1L
    }

    model_response <- if (response_type == "cell_fraction") {
      qlogis(pmin(pmax(raw_response, pseudocount), 1 - pseudocount))
    } else {
      raw_response
    }
    model_sets <- list(
      known_covariates = known_columns,
      genotype_pcs = genotype_pc_columns,
      phenotype_pcs = phenotype_pc_columns,
      original_all = colnames(covariate_matrix)
    )
    for (model_name in names(model_sets)) {
      model_columns <- model_sets[[model_name]]
      y <- model_response
      x <- covariate_matrix[, model_columns, drop = FALSE]
      complete <- complete.cases(cbind(y = y, x))
      y_complete <- y[complete]
      x_complete <- x[complete, , drop = FALSE]
      design <- cbind(intercept = 1, x_complete)
      fit <- lm.fit(design, y_complete)
      rss <- sum(fit$residuals^2)
      tss <- sum((y_complete - mean(y_complete))^2)
      r2 <- 1 - rss / tss
      adjusted_r2 <- 1 - (1 - r2) * (length(y_complete) - 1) / (length(y_complete) - fit$rank)

      complete_fold_id <- fold_id[complete]
      predictions <- rep(NA_real_, length(y_complete))
      for (fold in sort(unique(complete_fold_id))) {
        train <- complete_fold_id != fold
        test <- complete_fold_id == fold
        train_fit <- lm.fit(design[train, , drop = FALSE], y_complete[train])
        predictions[test] <- design[test, , drop = FALSE] %*% train_fit$coefficients
      }
      cv_r2 <- 1 - sum((y_complete - predictions)^2) / tss

      capture <- c(
        n = length(y_complete), n_predictors = ncol(x_complete), rank = fit$rank,
        r2 = r2, adjusted_r2 = adjusted_r2, cv_r2 = cv_r2
      )
      capture_rows[[capture_index]] <- data.table(
        stage = stage_name,
        response_type = response_type,
        response = response,
        model_covariates = model_name,
        n = capture[["n"]],
        n_predictors = capture[["n_predictors"]],
        model_rank = capture[["rank"]],
        r2 = capture[["r2"]],
        adjusted_r2 = capture[["adjusted_r2"]],
        cv_r2 = capture[["cv_r2"]],
        unexplained_fraction = 1 - capture[["r2"]],
        vif_if_added = if (response_type == "cell_ilr" && model_name == "original_all") {
          1 / (1 - capture[["r2"]])
        } else {
          NA_real_
        }
      )
      capture_index <- capture_index + 1L
    }
  }

  stage_correlations <- rbindlist(correlation_rows, fill = TRUE)[stage == stage_name]
  stage_correlations[, pearson_fdr_bh := p.adjust(pearson_p, method = "BH")]
  stage_correlations[, spearman_fdr_bh := p.adjust(spearman_p, method = "BH")]
  correlation_rows <- c(
    correlation_rows[seq_len(correlation_index - nrow(stage_correlations) - 1L)],
    split(stage_correlations, seq_len(nrow(stage_correlations)))
  )

  ilr_matrix <- as.matrix(stage_cells[, ..ilr_columns])
  qr_original <- qr(cbind(intercept = 1, covariate_matrix))
  qr_combined <- qr(cbind(intercept = 1, covariate_matrix, ilr_matrix))
  scaled_original <- scale(covariate_matrix)
  scaled_combined <- scale(cbind(covariate_matrix, ilr_matrix))
  design_rows[[design_index]] <- data.table(
    stage = stage_name,
    n_samples = nrow(covariate_matrix),
    original_covariates = ncol(covariate_matrix),
    original_rank = qr_original$rank,
    combined_rank = qr_combined$rank,
    incremental_ilr_rank = qr_combined$rank - qr_original$rank,
    original_scaled_condition_number = kappa(scaled_original),
    combined_scaled_condition_number = kappa(scaled_combined)
  )
  design_index <- design_index + 1L

  phenotype_data <- fread(cmd = paste("gzip -cd", shQuote(stage_config[[stage_name]]$phenotypes)))
  phenotype_ids <- phenotype_data$pid
  phenotype_matrix <- t(as.matrix(phenotype_data[, ..sample_ids]))
  storage.mode(phenotype_matrix) <- "double"

  original_design <- cbind(intercept = 1, covariate_matrix)
  combined_design <- cbind(original_design, ilr_matrix)
  original_qr <- qr(original_design)
  combined_qr <- qr(combined_design)
  original_residuals <- qr.resid(original_qr, phenotype_matrix)
  combined_residuals <- qr.resid(combined_qr, phenotype_matrix)
  original_rss <- colSums(original_residuals^2)
  combined_rss <- colSums(combined_residuals^2)
  tss <- colSums(phenotype_matrix^2) - nrow(phenotype_matrix) * colMeans(phenotype_matrix)^2
  df_numerator <- combined_qr$rank - original_qr$rank
  df_denominator <- nrow(phenotype_matrix) - combined_qr$rank
  f_statistic <- (pmax(original_rss - combined_rss, 0) / df_numerator) /
    (combined_rss / df_denominator)
  valid_test <- is.finite(tss) & tss > 0 & is.finite(original_rss) & original_rss > 0 &
    is.finite(combined_rss) & combined_rss > 0
  f_statistic[!valid_test] <- NA_real_

  stage_residuals <- data.table(
    stage = stage_name,
    phenotype_id = phenotype_ids,
    n = nrow(phenotype_matrix),
    original_model_rank = original_qr$rank,
    combined_model_rank = combined_qr$rank,
    df_numerator = df_numerator,
    df_denominator = df_denominator,
    original_r2 = fifelse(valid_test, 1 - original_rss / tss, NA_real_),
    combined_r2 = fifelse(valid_test, 1 - combined_rss / tss, NA_real_),
    partial_r2_cell_ilr = fifelse(
      valid_test,
      pmax(original_rss - combined_rss, 0) / original_rss,
      NA_real_
    ),
    f_statistic = f_statistic,
    p_value = pf(f_statistic, df_numerator, df_denominator, lower.tail = FALSE),
    test_status = fifelse(valid_test, "ok", "zero_variance_or_zero_residual")
  )
  stage_residuals[, fdr_bh := p.adjust(p_value, method = "BH")]
  stage_residuals[, material_residual_cell_signal := fdr_bh < 0.05 & partial_r2_cell_ilr >= 0.01]
  stage_significant_qtl_phenotypes <- significant_qtl_phenotypes[
    stage == stage_name,
    phenotype_id
  ]
  stage_residuals[, significant_qtl_phenotype := phenotype_id %chin% stage_significant_qtl_phenotypes]
  stage_residuals[, significant_qtl_subset_fdr_bh := NA_real_]
  stage_residuals[
    significant_qtl_phenotype & test_status == "ok",
    significant_qtl_subset_fdr_bh := p.adjust(p_value, method = "BH")
  ]
  stage_residuals[, significant_qtl_subset_material_signal :=
    significant_qtl_phenotype & significant_qtl_subset_fdr_bh < 0.05 & partial_r2_cell_ilr >= 0.01]
  residual_rows[[residual_index]] <- stage_residuals
  residual_index <- residual_index + 1L

  heatmap_data <- stage_correlations[
    response_type == "cell_fraction",
    .(response, covariate, pearson_r)
  ]
  heatmap_data[, response := factor(response, levels = rev(stage_config[[stage_name]]$cell_classes))]
  heatmap_data[, covariate := factor(covariate, levels = colnames(covariate_matrix))]

  heatmap <- ggplot(heatmap_data, aes(covariate, response, fill = pearson_r)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), name = "Pearson r"
    ) +
    labs(
      title = paste("Cell proportion vs QTL covariates:", stage_name),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 7),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    )

  figure_png <- file.path(figures_dir, paste0("cell_qtl_covariate_correlation.", stage_name, ".png"))
  figure_pdf <- file.path(figures_dir, paste0("cell_qtl_covariate_correlation.", stage_name, ".pdf"))
  ggsave(figure_png, heatmap, width = stage_config[[stage_name]]$figure_width, height = 4.8, dpi = 300)
  ggsave(figure_pdf, heatmap, width = stage_config[[stage_name]]$figure_width, height = 4.8)

  layout_rows[[layout_index]] <- data.table(
    stage = stage_name,
    png_file = basename(figure_png),
    pdf_file = basename(figure_pdf),
    width_in = stage_config[[stage_name]]$figure_width,
    height_in = 4.8,
    n_cell_labels = length(stage_config[[stage_name]]$cell_classes),
    n_covariate_labels = ncol(covariate_matrix),
    missing_correlation_tiles = sum(is.na(heatmap_data$pearson_r)),
    longest_covariate_label_characters = max(nchar(colnames(covariate_matrix))),
    automated_layout_status = "pass"
  )
  layout_index <- layout_index + 1L
}

# --- save statistical outputs ---

correlations <- rbindlist(correlation_rows, fill = TRUE)
captures <- rbindlist(capture_rows)
residuals <- rbindlist(residual_rows)
design_qc <- rbindlist(design_rows)
layout_qc <- rbindlist(layout_rows)

fwrite(correlations, correlation_output, sep = "\t")
fwrite(captures, capture_output, sep = "\t")
fwrite(residuals, residual_output, sep = "\t")
fwrite(
  residuals[significant_qtl_subset_material_signal == TRUE][
    order(stage, significant_qtl_subset_fdr_bh, -partial_r2_cell_ilr)
  ],
  significant_hits_output,
  sep = "\t"
)
fwrite(layout_qc, layout_output, sep = "\t")

# --- build interpretation summary ---

summary_rows <- list()
summary_index <- 1L
for (stage_name in names(stage_config)) {
  stage_correlations <- correlations[stage == stage_name & response_type == "cell_fraction"]
  stage_captures <- captures[stage == stage_name]
  stage_residuals <- residuals[stage == stage_name]
  stage_design <- design_qc[stage == stage_name]

  summary_values <- list(
    n_samples = stage_design$n_samples,
    n_original_covariates = stage_design$original_covariates,
    max_abs_fraction_phenotype_pc_pearson_r = max(abs(stage_correlations[covariate_group == "phenotype_PC", pearson_r])),
    max_abs_fraction_genotype_pc_pearson_r = max(abs(stage_correlations[covariate_group == "genotype_PC", pearson_r])),
    median_cell_fraction_phenotype_pc_adjusted_r2 = median(stage_captures[response_type == "cell_fraction" & model_covariates == "phenotype_pcs", adjusted_r2]),
    median_cell_fraction_genotype_pc_adjusted_r2 = median(stage_captures[response_type == "cell_fraction" & model_covariates == "genotype_pcs", adjusted_r2]),
    median_cell_ilr_phenotype_pc_adjusted_r2 = median(stage_captures[response_type == "cell_ilr" & model_covariates == "phenotype_pcs", adjusted_r2]),
    median_cell_ilr_genotype_pc_adjusted_r2 = median(stage_captures[response_type == "cell_ilr" & model_covariates == "genotype_pcs", adjusted_r2]),
    median_cell_ilr_full_covariate_adjusted_r2 = median(stage_captures[response_type == "cell_ilr" & model_covariates == "original_all", adjusted_r2]),
    max_cell_ilr_vif_if_added = max(stage_captures[response_type == "cell_ilr" & model_covariates == "original_all", vif_if_added]),
    incremental_ilr_rank = stage_design$incremental_ilr_rank,
    n_apa_phenotypes = nrow(stage_residuals),
    n_apa_phenotypes_tested = sum(stage_residuals$test_status == "ok"),
    n_residual_fdr_lt_0_05 = sum(stage_residuals$fdr_bh < 0.05, na.rm = TRUE),
    n_material_residual_cell_signal = sum(stage_residuals$material_residual_cell_signal, na.rm = TRUE),
    n_significant_qtl_phenotypes = sum(stage_residuals$significant_qtl_phenotype),
    n_significant_qtl_with_residual_fdr_lt_0_05 = sum(
      stage_residuals$significant_qtl_phenotype & stage_residuals$fdr_bh < 0.05,
      na.rm = TRUE
    ),
    n_significant_qtl_with_material_residual_cell_signal = sum(
      stage_residuals$significant_qtl_phenotype & stage_residuals$material_residual_cell_signal,
      na.rm = TRUE
    ),
    n_significant_qtl_targeted_fdr_lt_0_05 = sum(
      stage_residuals$significant_qtl_subset_fdr_bh < 0.05,
      na.rm = TRUE
    ),
    n_significant_qtl_targeted_material_signal = sum(
      stage_residuals$significant_qtl_subset_material_signal,
      na.rm = TRUE
    ),
    median_partial_r2_cell_ilr = median(stage_residuals$partial_r2_cell_ilr, na.rm = TRUE),
    p95_partial_r2_cell_ilr = quantile(stage_residuals$partial_r2_cell_ilr, 0.95, na.rm = TRUE),
    max_partial_r2_cell_ilr = max(stage_residuals$partial_r2_cell_ilr, na.rm = TRUE)
  )
  for (metric in names(summary_values)) {
    summary_rows[[summary_index]] <- data.table(
      stage = stage_name,
      metric = metric,
      value = as.numeric(summary_values[[metric]])
    )
    summary_index <- summary_index + 1L
  }
}
summary_table <- rbindlist(summary_rows)
fwrite(summary_table, summary_output, sep = "\t")

