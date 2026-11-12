#!/usr/bin/env Rscript
# Batch correction of APA loci: Raw, ComBat variants, limma, RUV-III and RUVs.
# Bridge pairs (batch9 replicates of donors measured in other batches) score each method.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.local.sh)")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tibble)
  library(ggrepel)
  library(limma)
  library(sva)
  library(dplyr)
  library(ruv)
  library(patchwork)
  library(variancePartition)
  library(RUVSeq)
  library(tidyr)
  library(cluster)
  library(purrr)
  library(softImpute)
  library(flashier)
})

# --- helpers ---
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

named_vector <- function(x) {
  if (is.null(x)) return(character())
  unlist(x, use.names = TRUE)
}

fig_size <- function(name) {
  fig_cfg$sizes[[name]] %||% list(width = 8, height = 5)
}

save_plot_png <- function(filename, plot_obj, ...) {
  args <- list(...)
  args$dpi <- args$dpi %||% fig_cfg$dpi
  formats <- fig_cfg$formats %||% c("png")
  base <- sub("\\.png$", "", filename)
  written <- character()
  for (fmt in formats) {
    out <- paste0(base, ".", fmt)
    device <- switch(
      fmt,
      png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
      pdf = grDevices::cairo_pdf,
      # NULL, not the "svg" string: ggsave's svg device needs svglite and would abort
      # the run. Unavailable formats are skipped and recorded in layout_validation.tsv.
      svg = if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else NULL,
      NULL
    )
    if (is.null(device)) next
    do.call(ggsave, c(list(filename = out, plot = plot_obj, device = device), args))
    written <- c(written, out)
  }
  plot_validation_records[[length(plot_validation_records) + 1]] <<- data.frame(
    artifact = filename,
    written_files = paste(written, collapse = ","),
    width = args$width %||% NA_real_,
    height = args$height %||% NA_real_,
    dpi = args$dpi %||% NA_real_,
    formats = paste(formats, collapse = ","),
    status = ifelse(length(written) > 0, "saved", "not_saved"),
    stringsAsFactors = FALSE
  )
  invisible(written)
}

write_tsv <- function(df, filename) {
  data.table::fwrite(df, filename, sep = "\t", quote = FALSE, na = "NA")
}

save_rds <- function(object, filename, compress = "gzip") {
  out <- file.path(out_dir, filename)
  saveRDS(object, out, compress = compress)
  out
}

theme_batch <- function(base_size = fig_cfg$base_size) {
  theme_minimal(base_size = base_size, base_family = fig_cfg$font_family) +
    theme(
      plot.title = element_text(face = "bold", size = fig_cfg$title_size,
                                color = fig_cfg$text,
                                hjust = fig_cfg$title_hjust %||% 0.5),
      plot.subtitle = element_text(size = fig_cfg$subtitle_size,
                                   color = fig_cfg$muted_text,
                                   hjust = fig_cfg$title_hjust %||% 0.5),
      axis.title = element_text(size = fig_cfg$label_size, color = fig_cfg$text),
      axis.text = element_text(size = fig_cfg$tick_size, color = fig_cfg$text),
      strip.text = element_text(face = "bold", size = fig_cfg$strip_size,
                                color = fig_cfg$text),
      strip.background = element_rect(fill = "#F7F7F7", color = NA),
      panel.grid.major = element_line(color = fig_cfg$grid, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = fig_cfg$panel_border, fill = NA,
                                  linewidth = 0.35),
      legend.position = fig_cfg$layout$legend_position %||% "bottom",
      legend.title = element_text(face = "bold", size = fig_cfg$tick_size),
      legend.text = element_text(size = fig_cfg$tick_size),
      plot.margin = margin(7, 8, 7, 8)
    )
}

rotate_x <- function(angle = fig_cfg$layout$axis_text_angle %||% 35) {
  theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1))
}

scale_method_fill <- function(...) {
  scale_fill_manual(values = method_palette, drop = FALSE, ...)
}

scale_batch_fill <- function(batch_levels, ...) {
  batch_colors <- named_vector(fig_cfg$batches)
  missing_batches <- setdiff(batch_levels, names(batch_colors))
  if (length(missing_batches) > 0) {
    fallback <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(missing_batches))
    names(fallback) <- missing_batches
    batch_colors <- c(batch_colors, fallback)
  }
  scale_fill_manual(values = batch_colors[batch_levels], drop = FALSE, ...)
}

as_numeric_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  mat
}

finite_feature_ids <- function(mat, min_var = 1e-8) {
  mat <- as_numeric_matrix(mat)
  finite_rows <- rowSums(!is.finite(mat)) == 0
  row_vars <- rep(NA_real_, nrow(mat))
  if (any(finite_rows)) {
    row_vars[finite_rows] <- apply(mat[finite_rows, , drop = FALSE], 1, var)
  }
  rownames(mat)[finite_rows & !is.na(row_vars) & row_vars > min_var]
}

filter_mat_for_eval <- function(mat, method_name, feature_ids = NULL,
                                sample_ids = NULL, min_var = 1e-8) {
  mat <- as_numeric_matrix(mat)
  if (!is.null(feature_ids)) {
    feature_ids <- intersect(feature_ids, rownames(mat))
    mat <- mat[feature_ids, , drop = FALSE]
  }
  if (!is.null(sample_ids)) {
    sample_ids <- intersect(sample_ids, colnames(mat))
    mat <- mat[, sample_ids, drop = FALSE]
  }
  keep <- rowSums(!is.finite(mat)) == 0
  if (any(keep)) {
    vars <- apply(mat[keep, , drop = FALSE], 1, var)
    keep[keep] <- !is.na(vars) & vars > min_var
  }
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) {
    stop(sprintf("%s has no finite variable loci available for evaluation.", method_name))
  }
  mat
}

find_elbow <- function(data) {
  df <- data.frame(k = data$k, var = data$batch)
  if (max(df$var, na.rm = TRUE) == min(df$var, na.rm = TRUE)) { return(df$k[1]) }
  df_norm <- data.frame(
    k_norm   = (df$k   - min(df$k))   / (max(df$k)   - min(df$k)),
    var_norm = (df$var - min(df$var)) / (max(df$var) - min(df$var))
  )
  p1 <- as.numeric(df_norm[1, ]); p_last <- as.numeric(df_norm[nrow(df_norm), ])
  distances <- sapply(1:nrow(df_norm), function(i) {
    p_i <- as.numeric(df_norm[i, ])
    v1 <- p_last - p1; v2 <- p_i - p1
    abs(v1[1] * v2[2] - v1[2] * v2[1]) / sqrt(sum(v1^2))
  })
  return(df$k[which.max(distances)])
}

mean_ranks <- function(metrics) {
  rowMeans(ranking_df[, paste0(metrics, "_rank"), drop = FALSE])
}

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
apa_file <- file.path(apa_root, "raw/all.apa.count")
meta_file <- file.path(apa_root, "data/04_qtl/CHB.all.meta")
out_dir <- file.path(apa_root, "data/01_quant")
fig_dir <- file.path(apa_root, "figures/01_quant")
selected_output_file <- file.path(apa_root, "data/04_qtl/apa_after_combat.tsv")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
plot_validation_records <- list()

config_path <- Sys.getenv("APA_CONFIG", file.path(apa_root, "config.yml"))
defaults <- list(
  dpi = 320,
  base_size = 10,
  font_family = "Arial",
  title_size = 13,
  title_hjust = 0.5,
  subtitle_size = 9,
  label_size = 9,
  tick_size = 8,
  strip_size = 9,
  stroke = 0.28,
  point_alpha = 0.82,
  grid = "#E6E6E6",
  panel_border = "#BDBDBD",
  text = "#222222",
  muted_text = "#5A5A5A",
  highlight = "#D55E00",
  selected = "#009E73",
  reference_line = "#CC79A7",
  bridge_arrow = "#D55E00",
  gradient_low = "#F7FBFF",
  gradient_mid = "#F7F7F7",
  gradient_high = "#D55E00",
  diverging_low = "#2166AC",
  diverging_mid = "#F7F7F7",
  diverging_high = "#B2182B",
  formats = c("png", "pdf", "svg"),
  stages = list(Prenatal = "#E69F00", Postnatal = "#0072B2"),
  pair_type = list(cross_batch = "#0072B2", within_batch = "#E69F00"),
  asw = list(ASW_batch = "#D55E00", ASW_LifeStage = "#0072B2"),
  methods = list(
    Raw = "#4D4D4D",
    ComBat_RAW = "#E69F00",
    ComBat_REF = "#009E73",
    ComBat_REF_NP = "#56B4E9",
    ComBat_REF_meanonly = "#CC79A7",
    limma = "#0072B2",
    `RUV-III` = "#F0E442",
    RUVs = "#D55E00"
  ),
  batches = list(),
  sizes = list(
    outlier = list(width = 7.0, height = 5.2),
    ruv_k = list(width = 8.0, height = 5.8),
    pca = list(width = 14.8, height = 10.2),
    residual = list(width = 11.5, height = 7.2),
    pc_heatmap = list(width = 12.8, height = 10.2),
    age_concordance = list(width = 8.8, height = 5.4),
    age_gam = list(width = 8.8, height = 5.4),
    bridge = list(width = 10.5, height = 5.6),
    variance_partition = list(width = 9.2, height = 5.8),
    asw = list(width = 9.2, height = 5.4),
    method_choice = list(width = 13.2, height = 7.8)
  ),
  layout = list(
    pca_ncol = 4,
    diagnostic_ncol = 4,
    pc_heatmap_ncol = 4,
    legend_position = "bottom",
    axis_text_angle = 35,
    compact_axis_text_angle = 45
  )
)
cfg <- list()
if (requireNamespace("yaml", quietly = TRUE) && file.exists(config_path)) {
  yaml_data <- yaml::read_yaml(config_path)
  plot_settings <- yaml_data$plot_settings %||% list()
  theme_settings <- yaml_data$theme %||% list()
  colors <- yaml_data$colors %||% list()
  cfg <- utils::modifyList(
    list(
      dpi = plot_settings$dpi,
      font_family = plot_settings$font_family,
      formats = plot_settings$formats,
      text = colors$base$ink,
      muted_text = colors$base$muted,
      grid = colors$base$grid,
      panel_border = colors$base$line,
      stages = colors$stages,
      title_hjust = theme_settings$title_hjust,
      layout = list(legend_position = theme_settings$legend_position)
    ),
    yaml_data$batch_correction_fig %||% list(),
    keep.null = TRUE
  )
}
fig_cfg <- utils::modifyList(defaults, cfg, keep.null = TRUE)
fig_cfg$font_family <- "sans"

method_palette <- named_vector(fig_cfg$methods)
stage_col <- named_vector(fig_cfg$stages)
pair_palette <- named_vector(fig_cfg$pair_type)
asw_palette <- named_vector(fig_cfg$asw)

# --- load ---
# APA count file: Gene | fit_value | Predicted_Proximal_APA | Loci | samples...
# Gene is "transcript|gene|chr|strand"; phenotype_id = transcript|gene|Loci|strand.
dt <- fread(apa_file, check.names = FALSE)
if (grepl("\\|", dt$Gene[1])) {
  dt[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed = TRUE)]
} else {
  stop("Gene column does not contain '|' — unexpected format.")
}
dt[, phenotype_id := paste(transcript, gene, Loci, strand, sep = "|")]

prefix_cols <- c(names(dt)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
loci <- dt[, .(Gene, fit_value, Predicted_Proximal_APA, Loci, chr, strand, phenotype_id)]
setDF(loci); rownames(loci) <- loci$phenotype_id

apa <- dt[, setdiff(names(dt), prefix_cols), with = FALSE]
setDF(apa); rownames(apa) <- dt$phenotype_id

meta <- fread(meta_file, data.table = FALSE)
rownames(meta) <- meta[, 1]
meta <- meta[, -1, drop = FALSE]

common_samples <- intersect(rownames(meta), colnames(apa))
meta <- meta[common_samples, , drop = FALSE]
apa <- apa[, common_samples]

# --- filter and impute ---
# QC thresholds follow 01_QTL/preprocess.r and 03_interaction/preprocess.r:
# sample dropped if NA rate > 80%; locus dropped if NA rate > 50% or SD < 0.001; autosomal only.
keep_samples <- colMeans(is.na(apa)) <= 0.8
apa  <- apa[, keep_samples, drop = FALSE]
meta <- meta[keep_samples, , drop = FALSE]

row_na_rate <- rowMeans(is.na(apa))
row_sd      <- apply(apa, 1, sd, na.rm = TRUE)
autosomal   <- !loci[rownames(apa), "chr"] %in% c("chrX", "chrY", "chrM")
keep_rows   <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & !is.na(row_sd) & autosomal
apa <- apa[keep_rows, ]

# Imputation (softImpute, then flashier), cached as apa_imputed.rds.
# The cache is keyed on the dimnames of the post-filter matrix, so a change in
# samples or loci invalidates it.
impute_cache <- file.path(out_dir, "apa_imputed.rds")
cached <- NULL
if (file.exists(impute_cache)) {
  cached <- tryCatch(readRDS(impute_cache), error = function(e) NULL)
  if (!is.matrix(cached) ||
      !identical(dimnames(cached), dimnames(as.matrix(apa)))) {
    message("Imputation cache does not match current matrix; recomputing.")
    cached <- NULL
  }
}

if (!is.null(cached)) {
  apa_imputed <- cached
} else {
  lambda <- if (nrow(apa) < 1500) 5 else 30
  X_mis_C <- as(as.matrix(apa), "Incomplete")
  fit_soft <- softImpute::softImpute(X_mis_C, rank = 50, lambda = lambda, type = "svd")
  pheno_soft <- softImpute::complete(as.matrix(apa), fit_soft)

  min_sd   <- min(apply(pheno_soft, 1, sd))
  pca_init <- irlba::irlba(pheno_soft, nv = 10)
  fl_pca <- flash_init(as.matrix(apa), S = min_sd, var_type = 1) |>
    flash_factors_init(list(d = pca_init$d, u = pca_init$u, v = pca_init$v),
                       ebnm_fn = ebnm::ebnm_point_laplace) |>
    flash_backfit(maxiter = 300)

  apa_imputed <- as.matrix(apa)
  na_idx <- is.na(apa_imputed)
  apa_imputed[na_idx] <- fitted(fl_pca)[na_idx]
  saveRDS(apa_imputed, impute_cache)
}

if (any(is.na(apa_imputed))) {
  message("Post-imputation NAs remain; filling with row means.")
  row_means <- rowMeans(apa_imputed, na.rm = TRUE)
  for (i in seq_len(nrow(apa_imputed)))
    apa_imputed[i, is.na(apa_imputed[i, ])] <- row_means[i]
}

# --- outliers ---
# Mahalanobis distance on the first 5 PCs within each LifeStage (as in 01_QTL/03_interaction).
n_pcs <- 5
p_cutoff <- 0.001
vars         <- apply(apa_imputed, 1, var)
expr_for_pca <- apa_imputed[vars > 1e-10, ]
pca_obj      <- prcomp(t(expr_for_pca), center = TRUE, scale. = TRUE)
pcs          <- pca_obj$x[, 1:n_pcs]

pca_data_list <- list()
for (stage in unique(meta$LifeStage)) {
  stage_samples <- rownames(meta)[meta$LifeStage == stage]
  stage_samples <- intersect(stage_samples, rownames(pcs))
  if (length(stage_samples) <= n_pcs) {
    warning(sprintf("Stage %s: too few samples for Mahalanobis (n=%d, df=%d). Skipping.",
                    stage, length(stage_samples), n_pcs))
    next
  }
  gp     <- pcs[stage_samples, ]
  m_dist <- mahalanobis(gp, colMeans(gp), cov(gp))
  p_vals <- 1 - pchisq(m_dist, df = n_pcs)
  pca_data_list[[stage]] <- data.frame(
    Sample    = stage_samples,
    PC1       = gp[, 1], PC2 = gp[, 2],
    M_Dist    = m_dist,
    P_Value   = p_vals,
    LifeStage = stage,
    IsOutlier = p_vals < p_cutoff
  )
}
pca_data <- do.call(rbind, pca_data_list)

outlier_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
  geom_point(aes(size = M_Dist, alpha = IsOutlier), shape = 16) +
  scale_color_manual(values = stage_col, drop = FALSE, name = "LifeStage") +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
  geom_text_repel(data = subset(pca_data, IsOutlier),
                  aes(label = Sample), color = fig_cfg$reference_line, size = 3,
                  fontface = "bold") +
  labs(title = "Outlier Mahalanobis Distance in PC Space",
       subtitle = sprintf("Labels: p < %s (chi-sq df=%d)", p_cutoff, n_pcs),
       x = "PC1", y = "PC2", size = "Mahalanobis Dist") +
  theme_batch()

size_outlier <- fig_size("outlier")
save_plot_png(file.path(fig_dir, "outlier_detection.png"), outlier_plot,
              width = size_outlier$width, height = size_outlier$height)

outlier_ids  <- pca_data$Sample[pca_data$IsOutlier]
samples_keep <- setdiff(colnames(apa_imputed), outlier_ids)
apa_imputed <- apa_imputed[, samples_keep]
meta_clean  <- meta[samples_keep, , drop = FALSE]

# Log2-logit transform (epsilon = 0.001, matching 01_QTL/03_interaction)
epsilon     <- 0.001
dat_bounded <- pmin(pmax(apa_imputed, epsilon), 1 - epsilon)
apa_scaled  <- log2(dat_bounded / (1 - dat_bounded))

# --- bridge pairs ---
batch9_samples <- rownames(meta_clean)[meta_clean$batch == "batch9"]
batch9_original_names <- gsub("_rep$", "", batch9_samples)

bridge_info <- list()
for (b in unique(meta_clean$batch)[unique(meta_clean$batch) != "batch9"]) {
  batch_samples  <- rownames(meta_clean)[meta_clean$batch == b]
  bridge_samples <- intersect(batch9_original_names, batch_samples)
  if (length(bridge_samples) > 0) {
    bridge_info[[b]] <- list(
      batch          = b,
      bridge_samples = bridge_samples,
      # match() keeps each replicate paired with its own donor
      batch9_reps    = batch9_samples[match(bridge_samples, batch9_original_names)]
    )
  }
}

# One row per donor pair; the pairing, control selection and bridge metrics read this.
bridge_pairs <- do.call(rbind, lapply(bridge_info, function(b) {
  data.frame(batch = b$batch,
             original_sample = b$bridge_samples,
             batch9_replicate = b$batch9_reps,
             stringsAsFactors = FALSE)
}))
rownames(bridge_pairs) <- NULL
stopifnot(nrow(bridge_pairs) > 0,
          all(bridge_pairs$original_sample %in% colnames(apa_scaled)),
          all(bridge_pairs$batch9_replicate %in% colnames(apa_scaled)))

# --- RUV negative controls ---
# A control locus must (a) not track the factor of interest and (b) respond to batch.
# (b) is estimated on bridge pairs, where biology cancels within a donor. (a) is
# estimated on batch9 only, where batch is constant. Batch is nearly aliased with
# LifeStage in this cohort, so the two conditions need different sample sets.

# (b) batch responsiveness per locus
d <- as.matrix(apa_scaled[, bridge_pairs$original_sample, drop = FALSE]) -
     as.matrix(apa_scaled[, bridge_pairs$batch9_replicate, drop = FALSE])
batches <- unique(bridge_pairs$batch)
delta_b <- vapply(batches, function(b) {
  rowMeans(d[, bridge_pairs$batch == b, drop = FALSE], na.rm = TRUE)
}, numeric(nrow(d)))
ss <- rep(0, nrow(d)); dfn <- rep(0, nrow(d))
for (b in batches) {
  dd <- d[, bridge_pairs$batch == b, drop = FALSE]
  if (ncol(dd) < 2) next
  ss  <- ss + rowSums((dd - rowMeans(dd, na.rm = TRUE))^2, na.rm = TRUE)
  dfn <- dfn + pmax(rowSums(is.finite(dd)) - 1, 0)
}
noise_var <- ss / pmax(dfn, 1)
signal    <- apply(delta_b, 1, var, na.rm = TRUE)
n_bar     <- nrow(bridge_pairs) / length(batches)
br <- data.frame(
  feature_id   = rownames(d),
  batch_signal = signal,
  batch_noise  = noise_var,
  batch_F      = signal / pmax(noise_var / n_bar, .Machine$double.eps),
  stringsAsFactors = FALSE, row.names = NULL
)

# (a) biological inertness on batch9 (age spline, sex and RIN adjusted)
ids <- rownames(meta_clean)[meta_clean$batch == "batch9"]
md  <- droplevels(meta_clean[ids, , drop = FALSE])
if (length(unique(md$LifeStage)) < 2) {
  stop("batch9 does not span both LifeStages; cannot estimate biological inertness.")
}
y  <- as.matrix(apa_scaled[, ids, drop = FALSE])
la <- log10(pmax(md$age + 1, 1e-6))
bs <- splines::ns(la, df = 3)
colnames(bs) <- paste0("age_ns", seq_len(ncol(bs)))
mdf <- cbind(md, as.data.frame(bs))
d_full <- model.matrix(as.formula(paste("~", paste(c(colnames(bs), "sex", "RIN"),
                                                   collapse = " + "))), data = mdf)
d_red  <- model.matrix(~ sex + RIN, data = mdf)
f_full <- limma::lmFit(y, d_full)
f_red  <- limma::lmFit(y, d_red)
ss_full <- f_full$sigma^2 * f_full$df.residual
ss_red  <- f_red$sigma^2  * f_red$df.residual
bio <- data.frame(
  feature_id  = rownames(y),
  bio_partial_r2 = pmax((ss_red - ss_full) / pmax(ss_red, .Machine$double.eps), 0),
  stringsAsFactors = FALSE, row.names = NULL
)

# control selection: batch-responsive (top half by batch F) and biologically inert
tab <- merge(br, bio, by = "feature_id", sort = FALSE)
q_batch <- quantile(tab$batch_F, 0.50, na.rm = TRUE)
q_bio   <- quantile(tab$bio_partial_r2, 0.40, na.rm = TRUE)
tab$batch_responsive <- is.finite(tab$batch_F) & tab$batch_F >= q_batch
tab$bio_inert        <- is.finite(tab$bio_partial_r2) & tab$bio_partial_r2 <= q_bio
tab$is_control       <- tab$batch_responsive & tab$bio_inert

# relax the inertness quantile rather than return a set too small to define k factors
relaxed <- 0.40
while (sum(tab$is_control) < 300 && relaxed < 0.9) {
  relaxed <- relaxed + 0.1
  q_bio <- quantile(tab$bio_partial_r2, relaxed, na.rm = TRUE)
  tab$bio_inert  <- is.finite(tab$bio_partial_r2) & tab$bio_partial_r2 <= q_bio
  tab$is_control <- tab$batch_responsive & tab$bio_inert
  message(sprintf("RUV controls: relaxed biological-inertness quantile to %.2f (n=%d).",
                  relaxed, sum(tab$is_control)))
}
control_table <- tab

control_table$mean_pdui <- rowMeans(dat_bounded[control_table$feature_id, , drop = FALSE],
                                    na.rm = TRUE)
stable_apa_loci <- control_table$feature_id[control_table$is_control]
if (length(stable_apa_loci) < 50) {
  stop("Too few RUV negative controls; refusing to run RUV with an ill-conditioned control set.")
}

write_tsv(control_table, file.path(out_dir, "ruv_control_selection.tsv"))
fwrite(data.frame(apa_locus_id = stable_apa_loci),
       file.path(out_dir, "stable_apa_loci.txt"),
       col.names = FALSE)

# control-set check: if the controls work, leading PCs track batch and not age or LifeStage
y <- as.matrix(apa_scaled[stable_apa_loci, , drop = FALSE])
v <- apply(y, 1, var, na.rm = TRUE)
y <- y[is.finite(v) & v > 1e-10, , drop = FALSE]
pcs <- prcomp(t(y), center = TRUE, scale. = TRUE)$x[, seq_len(5), drop = FALSE]
md  <- meta_clean[rownames(pcs), , drop = FALSE]
bo  <- model.matrix(~ batch + 0, data = md)
control_pc_check <- purrr::map_dfr(seq_len(5), function(i) {
  data.frame(
    pc = paste0("PC", i),
    max_abs_r_batch = max(abs(cor(pcs[, i], bo))),
    abs_r_age       = abs(cor(pcs[, i], log10(pmax(md$age + 1, 1e-6)),
                              use = "pairwise.complete.obs")),
    abs_r_lifestage = abs(cor(pcs[, i], as.numeric(factor(md$LifeStage)))),
    stringsAsFactors = FALSE
  )
})
write_tsv(control_pc_check, file.path(out_dir, "ruv_control_pc_check.tsv"))

# --- correction methods ---
count_raw  <- apa_scaled
meta_final <- meta_clean
meta_final$batch     <- as.factor(meta_final$batch)
meta_final$LifeStage <- as.factor(meta_final$LifeStage)
meta_final$sex       <- as.factor(meta_final$sex)

# ComBat model protects LifeStage. Age, sex and RIN enter the downstream design
# matrices instead of the ComBat adjustment model.
mod <- model.matrix(~ LifeStage, data = meta_final)
design_batch <- model.matrix(~ batch, data = meta_final)
full_design <- cbind(mod[, -1, drop = FALSE], design_batch[, -1, drop = FALSE])
if (qr(full_design)$rank != ncol(full_design)) {
  message("ComBat: LifeStage confounded with batch; using intercept-only model.")
  mod <- model.matrix(~ 1, data = meta_final)
}

combat_specs <- list(
  ComBat_RAW = list(par.prior = TRUE),
  ComBat_REF = list(ref.batch = "batch9", par.prior = TRUE),
  ComBat_REF_NP = list(ref.batch = "batch9", par.prior = FALSE),
  ComBat_REF_meanonly = list(ref.batch = "batch9", par.prior = TRUE, mean.only = TRUE)
)
combat_worker_request <- suppressWarnings(as.integer(
  Sys.getenv("APA_COMBAT_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
))
if (!is.finite(combat_worker_request) || combat_worker_request < 1L) combat_worker_request <- 1L
combat_workers <- max(1L, min(length(combat_specs), combat_worker_request))
combat_methods <- parallel::mclapply(
  names(combat_specs),
  function(method_label) {
    spec <- combat_specs[[method_label]]
    args <- list(
      dat = as.matrix(count_raw),
      batch = meta_final$batch,
      mod = mod,
      par.prior = spec$par.prior %||% TRUE,
      prior.plots = FALSE,
      mean.only = spec$mean.only %||% FALSE
    )
    if (!is.null(spec$ref.batch)) args$ref.batch <- spec$ref.batch
    res <- tryCatch(
      do.call(ComBat, args),
      error = function(e) {
        warning(sprintf("%s failed: %s", method_label, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    as.matrix(res)
  },
  mc.cores = combat_workers,
  mc.preschedule = FALSE
)
names(combat_methods) <- names(combat_specs)
combat_methods <- combat_methods[!vapply(combat_methods, is.null, logical(1))]

eval_worker_request <- suppressWarnings(as.integer(
  Sys.getenv("APA_EVAL_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
))
if (!is.finite(eval_worker_request) || eval_worker_request < 1L) eval_worker_request <- 1L
eval_workers <- max(1L, min(4L, eval_worker_request))

count_limma <- removeBatchEffect(
  count_raw,
  batch  = meta_final$batch,
  design = mod
)

# --- RUV-III ---
Y   <- t(count_raw)
ctl <- rownames(count_raw) %in% stable_apa_loci
k_range <- 0:10
varpart_formula_ruviii <- ~ (1|batch) + (1|LifeStage) + RIN + (1|sex)
M <- ruv::replicate.matrix(gsub("_rep$", "", rownames(meta_final)))
ruviii_results <- parallel::mclapply(k_range, function(k_val) {
  if (k_val == 0) {
    corrected_counts <- as.matrix(count_raw)
  } else {
    fit_ruviii <- RUVIII(Y = Y, M = M, ctl = ctl, k = k_val)
    if (is.null(fit_ruviii)) return(NULL)
    corrected_counts <- t(fit_ruviii)
  }

  set.seed(2025)
  gene_subset         <- sample(rownames(corrected_counts), min(2000, nrow(corrected_counts)))
  common_samples_loop <- intersect(colnames(corrected_counts), rownames(meta_final))
  varPart <- fitExtractVarPartModel(
    corrected_counts[gene_subset, common_samples_loop],
    varpart_formula_ruviii,
    meta_final[common_samples_loop, ]
  )
  result <- data.frame(
    method = "RUV-III", k = k_val, t(colMeans(varPart, na.rm = TRUE))
  )
  result
}, mc.cores = eval_workers, mc.preschedule = FALSE)
varpart_results_ruviii <- dplyr::bind_rows(ruviii_results)

# --- RUVs ---
# One row per donor pair. A per-batch row would tell RUVs that different donors are
# replicates of each other, which pushes between-donor biology into the residual.
scIdx <- as.matrix(data.frame(
  original = match(bridge_pairs$original_sample,  colnames(count_raw)),
  replicate = match(bridge_pairs$batch9_replicate, colnames(count_raw))
))
stopifnot(!anyNA(scIdx), ncol(scIdx) == 2, nrow(scIdx) == nrow(bridge_pairs))
ctl_character <- rownames(count_raw)[rownames(count_raw) %in% stable_apa_loci]
ruvs_results <- parallel::mclapply(k_range, function(k_val) {
  if (k_val == 0) {
    count_corr      <- count_raw
    metadata_for_vp <- meta_final
  } else {
    set_corr        <- RUVs(as.matrix(count_raw), cIdx = ctl_character,
                            k = k_val, scIdx = scIdx, isLog = TRUE)
    count_corr      <- set_corr$normalizedCounts
    metadata_for_vp <- cbind(meta_final, set_corr$W)
  }

  random_effects     <- "(1|batch) + (1|LifeStage) + (1|sex)"
  fixed_effects_base <- c("RIN")
  if (k_val > 0) {
    w_factors         <- paste0("W_", 1:k_val)
    fixed_effects_all <- c(fixed_effects_base, w_factors)
  } else {
    fixed_effects_all <- fixed_effects_base
  }
  varpart_formula_final <- as.formula(
    paste("~", random_effects, "+", paste(fixed_effects_all, collapse = " + "))
  )

  set.seed(2025)
  gene_subset_vp <- sample(rownames(count_corr), min(2000, nrow(count_corr)))
  varPart <- fitExtractVarPartModel(count_corr[gene_subset_vp, ], varpart_formula_final, metadata_for_vp)
  result <- data.frame(
    method = "RUVs", k = k_val, t(colMeans(varPart, na.rm = TRUE))
  )
  result
}, mc.cores = eval_workers, mc.preschedule = FALSE)
varpart_results_ruvs <- dplyr::bind_rows(ruvs_results)

varpart_results <- dplyr::bind_rows(varpart_results_ruviii, varpart_results_ruvs)

k_opt_ruviii <- find_elbow(varpart_results %>% filter(method == "RUV-III"))
k_opt_ruvs   <- find_elbow(varpart_results %>% filter(method == "RUVs"))

varpart_long <- varpart_results %>%
  pivot_longer(cols = all_of(c("batch", "LifeStage", "RIN", "sex")),
               names_to = "variable", values_to = "variance")

p_diagnostic <- ggplot(varpart_long, aes(x = k, y = variance, color = method, group = method)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_vline(xintercept = k_opt_ruviii, linetype = "dashed", color = method_palette[["RUV-III"]]) +
  geom_vline(xintercept = k_opt_ruvs,   linetype = "dashed", color = method_palette[["RUVs"]]) +
  facet_wrap(~variable, scales = "free_y") +
  scale_color_manual(values = method_palette, drop = FALSE, name = "Method") +
  labs(
    title = "Variance Partitioning Diagnostics for K-Selection",
    x     = "Number of Unwanted Factors (k)",
    y     = "Mean Proportion of Variance Explained"
  ) +
  scale_x_continuous(breaks = k_range) +
  theme_batch()

size_ruv <- fig_size("ruv_k")
save_plot_png(file.path(fig_dir, "ruv_k_selection.png"), p_diagnostic,
              width = size_ruv$width, height = size_ruv$height)

# --- apply optimal k ---
count_ruv_III <- as.data.frame(t(RUVIII(Y = Y, M = M, ctl = ctl, k = k_opt_ruviii)))
set_corr_final <- RUVs(as.matrix(count_raw), cIdx = ctl_character,
                       k = k_opt_ruvs, scIdx = scIdx, isLog = TRUE)
count_ruv_s <- as.data.frame(set_corr_final$normalizedCounts)

methods <- c(
  list(Raw = count_raw),
  combat_methods,
  list(
    limma     = count_limma,
    `RUV-III` = count_ruv_III,
    `RUVs`    = count_ruv_s
  )
)
method_levels <- c("Raw", setdiff(names(methods), "Raw"))

method_nonfinite_summary <- purrr::imap_dfr(methods, function(mat, method_name) {
  mat <- as_numeric_matrix(mat)
  nonfinite <- !is.finite(mat)
  data.frame(
    method = method_name,
    n_features = nrow(mat),
    n_samples = ncol(mat),
    n_values = length(mat),
    n_nonfinite_values = sum(nonfinite),
    n_na_values = sum(is.na(mat)),
    n_pos_inf_values = sum(mat == Inf, na.rm = TRUE),
    n_neg_inf_values = sum(mat == -Inf, na.rm = TRUE),
    n_features_with_nonfinite = sum(rowSums(nonfinite) > 0),
    n_samples_with_nonfinite = sum(colSums(nonfinite) > 0),
    n_finite_variable_features = length(finite_feature_ids(mat)),
    stringsAsFactors = FALSE
  )
})
write_tsv(method_nonfinite_summary, file.path(out_dir, "method_nonfinite_summary.tsv"))

nonfinite_file <- file.path(out_dir, "method_nonfinite_examples.tsv")
nonfinite_examples <- purrr::imap_dfr(methods, function(mat, method_name) {
  mat <- as_numeric_matrix(mat)
  bad_rows <- rownames(mat)[rowSums(!is.finite(mat)) > 0]
  if (length(bad_rows) == 0) return(NULL)
  bad_rows <- head(bad_rows, 200)
  purrr::map_dfr(bad_rows, function(feature_id) {
    bad_samples <- colnames(mat)[!is.finite(mat[feature_id, ])]
    data.frame(
      method = method_name,
      feature = feature_id,
      n_bad_samples = length(bad_samples),
      bad_samples = paste(head(bad_samples, 20), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
})
if (nrow(nonfinite_examples) > 0) {
  write_tsv(nonfinite_examples, nonfinite_file)
} else if (file.exists(nonfinite_file)) {
  file.remove(nonfinite_file)
}

if (any(method_nonfinite_summary$n_nonfinite_values > 0)) {
  warning(
    paste(
      "Non-finite values detected in corrected matrices; evaluation will use finite loci.",
      paste(
        method_nonfinite_summary$method,
        method_nonfinite_summary$n_nonfinite_values,
        sep = "=",
        collapse = "; "
      )
    )
  )
}

common_eval_genes <- Reduce(intersect, lapply(methods, finite_feature_ids))
common_eval_samples <- Reduce(intersect, lapply(methods, colnames))
if (length(common_eval_genes) < 2 || length(common_eval_samples) < 2) {
  stop("Too few common finite loci or samples for evaluation.")
}

bridge_df_for_eval <- bridge_pairs[, c("original_sample", "batch9_replicate")]

method_matrix_manifest <- purrr::imap_dfr(methods, function(mat, method_name) {
  out <- save_rds(mat, paste0("method_matrix_", gsub("[^A-Za-z0-9_.-]+", "_", method_name), ".rds"))
  data.frame(
    method = method_name,
    rds_file = out,
    n_features = nrow(mat),
    n_samples = ncol(mat),
    stringsAsFactors = FALSE
  )
})
write_tsv(method_matrix_manifest, file.path(out_dir, "method_matrix_rds_manifest.tsv"))

save_rds(
  list(
    meta_final = meta_final,
    bridge_info = bridge_info,
    bridge_df_for_eval = bridge_df_for_eval,
    stable_apa_loci = stable_apa_loci,
    k_opt_ruviii = k_opt_ruviii,
    k_opt_ruvs = k_opt_ruvs,
    varpart_results = varpart_results,
    method_nonfinite_summary = method_nonfinite_summary,
    common_eval_genes = common_eval_genes,
    common_eval_samples = common_eval_samples,
    method_matrix_manifest = method_matrix_manifest,
    fig_cfg = fig_cfg
  ),
  "batch_correction_checkpoint.rds"
)

# --- evaluation: PCA ---
# Unified PCA over all methods (same axes), with bridge-pair arrows per panel.
all_data <- do.call(cbind, lapply(names(methods), function(method_name) {
  mat <- filter_mat_for_eval(
    methods[[method_name]],
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_eval_samples
  )
  colnames(mat) <- paste0(colnames(mat), "___", method_name)
  mat
}))

pca_data <- t(all_data)
pca_data <- pca_data[, apply(pca_data, 2, var) > 1e-8, drop = FALSE]
pca    <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample_method") %>%
  tidyr::separate(sample_method, into = c("sample", "method"), sep = "___")
var_exp <- summary(pca)$importance[2, 1:2] * 100

plot_df <- meta_final %>%
  rownames_to_column("sample") %>%
  right_join(pca_df, by = "sample") %>%
  mutate(method = factor(method, levels = names(methods)))

bridge_conn <- bridge_df_for_eval %>%
  tidyr::crossing(method = unique(plot_df$method)) %>%
  left_join(plot_df %>% select(sample, method, PC1_o = PC1, PC2_o = PC2),
            by = c("original_sample" = "sample", "method" = "method")) %>%
  left_join(plot_df %>% select(sample, method, PC1_r = PC1, PC2_r = PC2),
            by = c("batch9_replicate" = "sample", "method" = "method")) %>%
  filter(!is.na(PC1_o) & !is.na(PC1_r)) %>%
  mutate(dist = sqrt((PC1_o - PC1_r)^2 + (PC2_o - PC2_r)^2))

write_tsv(plot_df, file.path(out_dir, "pca_comparison_scores.tsv"))
write_tsv(bridge_conn, file.path(out_dir, "pca_bridge_segments.tsv"))

method_dist <- bridge_conn %>%
  group_by(method) %>%
  summarise(mean_dist = mean(dist, na.rm = TRUE))

batch_levels <- sort(unique(plot_df$batch))
group_shapes <- c(Prenatal = 21, Postnatal = 22)

plots_no_legend <- list()
for (method_name in names(methods)) {
  method_data   <- plot_df     %>% filter(method == method_name)
  method_bridge <- bridge_conn %>% filter(method == method_name)
  mean_dist     <- method_dist %>% filter(method == method_name) %>% pull(mean_dist)

  p <- ggplot(method_data, aes(x = PC1, y = PC2)) +
    geom_segment(data = method_bridge,
                 aes(x = PC1_o, y = PC2_o, xend = PC1_r, yend = PC2_r),
                 arrow = arrow(length = unit(0.1, "cm")),
                 color = fig_cfg$bridge_arrow, alpha = 0.45, linewidth = 0.45) +
    geom_point(aes(fill = batch, shape = LifeStage), size = 2.45, alpha = fig_cfg$point_alpha,
               color = "black", stroke = fig_cfg$stroke) +
    scale_batch_fill(batch_levels, name = "Batch") +
    scale_shape_manual(values = group_shapes, name = "LifeStage") +
    labs(
      title    = method_name,
      subtitle = sprintf("Bridge distance: %.1f", mean_dist),
      x = sprintf("PC1 (%.1f%%)", var_exp[1]),
      y = sprintf("PC2 (%.1f%%)", var_exp[2])
    ) +
    theme_batch(base_size = max(fig_cfg$base_size - 1, 10)) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = fig_cfg$bridge_arrow)
    )

  plots_no_legend[[method_name]] <- p
}

legend_plot <- ggplot() +
  geom_point(data = plot_df %>% distinct(batch, LifeStage),
             aes(x = 0, y = 0, fill = batch, shape = LifeStage),
             size = 4, stroke = 0.5, color = "black") +
  scale_batch_fill(batch_levels, name = "Batch") +
  scale_shape_manual(values = group_shapes, name = "LifeStage") +
  theme_void() +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, size = 4), ncol = 5),
    shape = guide_legend(override.aes = list(fill = "grey70", size = 4))
  ) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",
    legend.margin   = margin(5, 5, 5, 5),
    legend.title    = element_text(face = "bold")
  )

legend <- cowplot::get_legend(legend_plot)

combined_plot <- patchwork::wrap_plots(plots_no_legend, ncol = fig_cfg$layout$pca_ncol %||% 4) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Comparison of Batch Correction Methods (APA)",
    theme = theme(plot.title = element_text(face = "bold", size = fig_cfg$title_size + 2,
                                            hjust = fig_cfg$title_hjust %||% 0.5,
                                            family = fig_cfg$font_family))
  )
final_plot <- cowplot::plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(0.9, 0.1))

size_pca <- fig_size("pca")
save_plot_png(file.path(fig_dir, "pca_comparison.png"), final_plot,
              width = size_pca$width, height = size_pca$height)

# --- diagnostics: residual batch offsets ---
residual_batch_results <- parallel::mclapply(names(methods), function(method_name) {
  mat <- methods[[method_name]]
  common_samples <- intersect(colnames(mat), rownames(meta_final))
  meta_aligned <- meta_final[common_samples, , drop = FALSE]
  mat <- filter_mat_for_eval(
    mat,
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_samples
  )
  design <- model.matrix(~ 0 + batch + LifeStage + age + sex + RIN, data = meta_aligned)
  batch_cols <- grep("^batch", colnames(design), value = TRUE)
  if (length(batch_cols) == 0) return(data.frame() %>% mutate(method = method_name, .before = 1))
  fit <- limma::lmFit(mat, design = design)
  coefs <- fit$coefficients[, batch_cols, drop = FALSE]
  batch_offsets <- sweep(coefs, 1, rowMeans(coefs, na.rm = TRUE), "-")
  data.frame(
    batch_level = sub("^batch", "", batch_cols),
    median_abs_coef = apply(abs(batch_offsets), 2, median, na.rm = TRUE),
    p95_abs_coef = apply(abs(batch_offsets), 2, function(x) quantile(x, 0.95, na.rm = TRUE)),
    p99_abs_coef = apply(abs(batch_offsets), 2, function(x) quantile(x, 0.99, na.rm = TRUE)),
    n_features_coef_gt_1 = colSums(abs(batch_offsets) > 1, na.rm = TRUE),
    row.names = NULL
  ) %>% mutate(method = method_name, .before = 1)
}, mc.cores = eval_workers, mc.preschedule = FALSE)
residual_batch_df <- dplyr::bind_rows(residual_batch_results)
residual_batch_df <- residual_batch_df %>%
  mutate(batch_level = factor(batch_level, levels = levels(meta_final$batch)))
write_tsv(residual_batch_df, file.path(out_dir, "residual_batch_coefficients.tsv"))

summary_residual <- residual_batch_df %>%
  group_by(method) %>%
  summarise(
    ResidualCoef_median = median(median_abs_coef, na.rm = TRUE),
    ResidualCoef_p95 = median(p95_abs_coef, na.rm = TRUE),
    ResidualCoef_max_p95 = max(p95_abs_coef, na.rm = TRUE),
    ResidualCoef_gt1_total = sum(n_features_coef_gt_1, na.rm = TRUE),
    .groups = "drop"
  )

p_residual <- ggplot(residual_batch_df, aes(x = batch_level, y = median_abs_coef, fill = method)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.25) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = fig_cfg$reference_line) +
  facet_wrap(~ method, ncol = fig_cfg$layout$diagnostic_ncol %||% 4, scales = "free_x") +
  scale_method_fill(name = "Method") +
  labs(title = "Residual Batch Coefficients", x = "Batch", y = "Median |batch offset|") +
  theme_batch() +
  rotate_x(fig_cfg$layout$compact_axis_text_angle %||% 45) +
  theme(legend.position = "none")
size_residual <- fig_size("residual")
save_plot_png(file.path(fig_dir, "residual_batch_coefficients.png"), p_residual,
              width = size_residual$width, height = size_residual$height)

# --- diagnostics: PC correlation with batch and covariates ---
pc_batch_results <- parallel::mclapply(names(methods), function(method_name) {
  mat <- methods[[method_name]]
  common_samples <- intersect(colnames(mat), rownames(meta_final))
  meta_aligned <- meta_final[common_samples, , drop = FALSE]
  mat <- filter_mat_for_eval(
    mat,
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_samples
  )
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  n_pcs <- min(10, ncol(pca$x))
  pcs <- pca$x[, seq_len(n_pcs), drop = FALSE]
  colnames(pcs) <- paste0("PC", seq_len(n_pcs))
  batch_onehot <- model.matrix(~ batch + 0, data = meta_aligned)

  batch_rows <- tidyr::expand_grid(
    pc = colnames(pcs),
    variable = colnames(batch_onehot)
  ) %>%
    mutate(pearson_r = purrr::map2_dbl(pc, variable, ~ cor(pcs[, .x], batch_onehot[, .y])))

  covar_rows <- lapply(c("age", "RIN"), function(var) {
    data.frame(
      pc = colnames(pcs),
      variable = var,
      pearson_r = apply(pcs, 2, cor, y = meta_aligned[[var]], use = "pairwise.complete.obs"),
      row.names = NULL
    )
  }) %>% bind_rows()

  sex_onehot <- model.matrix(~ sex + 0, data = meta_aligned)
  sex_rows <- tidyr::expand_grid(
    pc = colnames(pcs),
    variable = colnames(sex_onehot)
  ) %>%
    mutate(pearson_r = purrr::map2_dbl(pc, variable, ~ cor(pcs[, .x], sex_onehot[, .y])))

  bind_rows(batch_rows, covar_rows, sex_rows) %>% mutate(method = method_name, .before = 1)
}, mc.cores = eval_workers, mc.preschedule = FALSE)
pc_batch_df <- dplyr::bind_rows(pc_batch_results)
write_tsv(pc_batch_df, file.path(out_dir, "pc_batch_covariate_correlations.tsv"))

summary_pc_batch <- pc_batch_df %>%
  filter(grepl("^batch", variable), pc %in% paste0("PC", 1:3)) %>%
  group_by(method) %>%
  summarise(PC_batch_max_abs_r = max(abs(pearson_r), na.rm = TRUE), .groups = "drop")

p_pc_batch <- ggplot(pc_batch_df %>% filter(pc %in% paste0("PC", 1:10)),
                     aes(x = variable, y = pc, fill = pearson_r)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = fig_cfg$diverging_low, mid = fig_cfg$diverging_mid,
                       high = fig_cfg$diverging_high,
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  facet_wrap(~ method, ncol = fig_cfg$layout$pc_heatmap_ncol %||% 4) +
  labs(title = "PC Correlation With Batch and Covariates", x = NULL, y = NULL) +
  theme_batch(base_size = max(fig_cfg$base_size - 2, 9)) +
  rotate_x(60)
size_pc <- fig_size("pc_heatmap")
save_plot_png(file.path(fig_dir, "pc_batch_covariate_correlations.png"), p_pc_batch,
              width = size_pc$width, height = size_pc$height)

# --- diagnostics: postnatal age-matched concordance ---
age_concordance_results <- parallel::mclapply(names(methods), function(method_name) {
  mat <- methods[[method_name]]
  common_samples <- intersect(colnames(mat), rownames(meta_final))
  meta_aligned <- meta_final[common_samples, , drop = FALSE]
  post_meta <- meta_aligned[meta_aligned$LifeStage == "Postnatal" & !is.na(meta_aligned$age), , drop = FALSE]
  post_meta <- post_meta[order(post_meta$age), , drop = FALSE]

  pair_df <- NULL
  if (nrow(post_meta) >= 2) {
    pair_df <- utils::combn(rownames(post_meta), 2, simplify = FALSE) %>%
      purrr::map_dfr(function(pair) {
        m1 <- post_meta[pair[1], , drop = FALSE]
        m2 <- post_meta[pair[2], , drop = FALSE]
        if (abs(m1$age - m2$age) > 5) return(NULL)
        data.frame(
          sample1 = pair[1],
          sample2 = pair[2],
          pair_type = ifelse(m1$batch == m2$batch, "within_batch", "cross_batch"),
          stringsAsFactors = FALSE
        )
      })
  }

  if (is.null(pair_df) || nrow(pair_df) == 0) {
    res <- data.frame(pair_type = c("cross_batch", "within_batch"), n_pairs = 0,
                      mean_cor = NA_real_, sd_cor = NA_real_)
  } else {
    mat <- filter_mat_for_eval(
      mat,
      method_name,
      feature_ids = common_eval_genes,
      sample_ids = common_samples
    )
    pair_df$correlation <- purrr::map2_dbl(pair_df$sample1, pair_df$sample2, ~ cor(mat[, .x], mat[, .y]))
    res <- pair_df %>%
      group_by(pair_type) %>%
      summarise(n_pairs = n(), mean_cor = mean(correlation, na.rm = TRUE),
                sd_cor = sd(correlation, na.rm = TRUE), .groups = "drop")
  }
  res %>% mutate(method = method_name, .before = 1)
}, mc.cores = eval_workers, mc.preschedule = FALSE)
age_concordance_df <- dplyr::bind_rows(age_concordance_results)
write_tsv(age_concordance_df, file.path(out_dir, "postnatal_age_matched_concordance.tsv"))

summary_age_concordance <- age_concordance_df %>%
  select(method, pair_type, mean_cor) %>%
  pivot_wider(names_from = pair_type, values_from = mean_cor)
if (!"within_batch" %in% names(summary_age_concordance)) summary_age_concordance$within_batch <- NA_real_
if (!"cross_batch" %in% names(summary_age_concordance)) summary_age_concordance$cross_batch <- NA_real_
summary_age_concordance <- summary_age_concordance %>%
  mutate(AgeMatched_gap = within_batch - cross_batch)

p_age_concordance <- ggplot(age_concordance_df, aes(x = method, y = mean_cor, fill = pair_type)) +
  geom_col(position = "dodge", color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pair_palette, drop = FALSE, name = "Pair type") +
  labs(title = "Postnatal Age-Matched Concordance", x = NULL, y = "Mean Pearson correlation", fill = "Pair type") +
  theme_batch() +
  rotate_x()
size_age <- fig_size("age_concordance")
save_plot_png(file.path(fig_dir, "postnatal_age_matched_concordance.png"), p_age_concordance,
              width = size_age$width, height = size_age$height)

# --- diagnostics: age-GAM preservation ---
# Per locus, a REML GAM on age (k ladder, k-index check) against a reduced model
# without age. Supported grid points are the ages with enough local samples.
age_gam_results <- parallel::mclapply(names(methods), function(method_name) {
  mat <- methods[[method_name]]
  common_samples <- intersect(colnames(mat), rownames(meta_final))
  meta_aligned <- meta_final[common_samples, , drop = FALSE]
  meta_aligned$age_model <- log10(pmax(as.numeric(meta_aligned$age) + 1, 1e-6))
  mat <- filter_mat_for_eval(
    mat,
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_samples
  )
  set.seed(2025)
  feat_subset <- sample(rownames(mat), min(5000, nrow(mat)))

  age_range <- diff(range(meta_aligned$age_model, na.rm = TRUE))
  grid_model <- seq(min(meta_aligned$age_model, na.rm = TRUE),
                    max(meta_aligned$age_model, na.rm = TRUE),
                    length.out = 128)
  half_w_outer <- 0.08 * age_range
  half_w_inner <- 0.04 * age_range
  local_n <- vapply(grid_model, function(g) sum(abs(meta_aligned$age_model - g) <= half_w_outer), integer(1))
  local_n_inner <- vapply(grid_model, function(g) sum(abs(meta_aligned$age_model - g) <= half_w_inner), integer(1))
  supported <- local_n >= 5L & local_n_inner >= 3L
  support_df <- data.frame(
    grid_index = seq_along(grid_model) - 1L,
    age_model = grid_model,
    grid_age = 10^grid_model - 1,
    local_n = local_n,
    local_n_inner = local_n_inner,
    is_supported = supported
  )

  model_df <- as.data.frame(meta_aligned[, c("age_model", "sex", "RIN", "batch")])
  grid_df <- data.frame(
    age_model = grid_model,
    sex = factor(levels(meta_aligned$sex)[which.max(tabulate(meta_aligned$sex))],
                 levels = levels(meta_aligned$sex)),
    RIN = median(meta_aligned$RIN, na.rm = TRUE),
    batch = factor(levels(meta_aligned$batch)[which.max(tabulate(meta_aligned$batch))],
                   levels = levels(meta_aligned$batch))
  )
  reduced_formula <- y ~ sex + RIN + batch
  k_ladder <- c(20L, 30L, 40L, 60L)
  k_index_p_thr <- 0.05
  edf_margin <- 1.0

  res_df <- purrr::map_dfr(feat_subset, function(f) {
    model_df$y <- as.numeric(mat[f, rownames(model_df)])
    if (anyNA(model_df$y)) model_df$y[is.na(model_df$y)] <- mean(model_df$y, na.rm = TRUE)
    fit_red <- tryCatch(mgcv::gam(reduced_formula, data = model_df, method = "REML"), error = function(e) NULL)
    if (is.null(fit_red)) {
      return(data.frame(
        feature = f, final_k = NA_integer_, smooth_edf = NA_real_,
        k_index = NA_real_, k_index_pvalue = NA_real_, age_pvalue = NA_real_,
        partial_deviance = NA_real_, delta_logit = NA_real_,
        amplitude_logit = NA_real_, row.names = NULL
      ))
    }

    best <- NULL
    for (k in k_ladder) {
      full_formula <- as.formula(sprintf("y ~ sex + RIN + batch + s(age_model, k = %d, bs = 'tp')", k))
      fit_full <- tryCatch(mgcv::gam(full_formula, data = model_df, method = "REML"), error = function(e) NULL)
      if (is.null(fit_full)) next
      s_tab <- summary(fit_full)$s.table
      edf <- suppressWarnings(as.numeric(s_tab[1, "edf"]))
      kc <- tryCatch(
        mgcv::k.check(fit_full, subsample = 5000, n.rep = 400),
        error = function(e) matrix(NA_real_, 1, 4,
          dimnames = list(NULL, c("k'", "edf", "k-index", "p-value")))
      )
      k_idx <- suppressWarnings(as.numeric(kc[1, "k-index"]))
      k_p <- suppressWarnings(as.numeric(kc[1, "p-value"]))
      dev_res <- tryCatch(anova(fit_red, fit_full, test = "F"), error = function(e) NULL)
      age_p <- if (!is.null(dev_res)) suppressWarnings(as.numeric(dev_res$`Pr(>F)`[2])) else
        suppressWarnings(as.numeric(s_tab[1, "p-value"]))
      if (!is.finite(age_p)) age_p <- suppressWarnings(as.numeric(s_tab[1, "p-value"]))
      best <- list(fit = fit_full, final_k = k, smooth_edf = edf,
                   k_index = k_idx, k_index_pvalue = k_p, age_pvalue = age_p)
      k_ok <- is.finite(k_p) && k_p > k_index_p_thr &&
        !(is.finite(edf) && edf >= k - edf_margin)
      if (k_ok) break
    }
    if (is.null(best)) {
      return(data.frame(
        feature = f, final_k = NA_integer_, smooth_edf = NA_real_,
        k_index = NA_real_, k_index_pvalue = NA_real_, age_pvalue = NA_real_,
        partial_deviance = NA_real_, delta_logit = NA_real_,
        amplitude_logit = NA_real_, row.names = NULL
      ))
    }
    pred <- as.numeric(predict(best$fit, newdata = grid_df))
    sup_idx <- which(supported)
    if (length(sup_idx) >= 2) {
      pred_sup <- pred[sup_idx]
      delta_logit <- pred_sup[length(pred_sup)] - pred_sup[1]
      amplitude_logit <- max(pred_sup, na.rm = TRUE) - min(pred_sup, na.rm = TRUE)
    } else {
      delta_logit <- NA_real_
      amplitude_logit <- NA_real_
    }
    dev_red <- deviance(fit_red)
    dev_full <- deviance(best$fit)
    data.frame(
      feature = f,
      final_k = best$final_k,
      smooth_edf = best$smooth_edf,
      k_index = best$k_index,
      k_index_pvalue = best$k_index_pvalue,
      age_pvalue = best$age_pvalue,
      partial_deviance = (dev_red - dev_full) / dev_red,
      delta_logit = delta_logit,
      amplitude_logit = amplitude_logit,
      row.names = NULL
    )
  })

  write_tsv(res_df, file.path(out_dir, paste0("age_gam_preservation_features_", method_name, ".tsv")))
  write_tsv(support_df %>% mutate(method = method_name, .before = 1),
            file.path(out_dir, paste0("age_support_profile_", method_name, ".tsv")))
  data.frame(
    n_tested = nrow(res_df),
    n_sig_p001 = sum(res_df$age_pvalue < 0.001, na.rm = TRUE),
    n_sig_p01 = sum(res_df$age_pvalue < 0.01, na.rm = TRUE),
    median_partial_deviance = median(res_df$partial_deviance, na.rm = TRUE),
    median_abs_delta_logit = median(abs(res_df$delta_logit), na.rm = TRUE),
    median_amplitude_logit = median(res_df$amplitude_logit, na.rm = TRUE),
    p01_fraction = mean(res_df$age_pvalue < 0.01, na.rm = TRUE),
    k_index_fail_fraction = mean(res_df$k_index_pvalue <= 0.05, na.rm = TRUE),
    median_smooth_edf = median(res_df$smooth_edf, na.rm = TRUE),
    supported_grid_points = sum(support_df$is_supported),
    total_grid_points = nrow(support_df)
  ) %>% mutate(method = method_name, .before = 1)
}, mc.cores = eval_workers, mc.preschedule = FALSE)
summary_age_gam <- dplyr::bind_rows(age_gam_results)
write_tsv(summary_age_gam, file.path(out_dir, "age_gam_preservation_summary.tsv"))

p_age_gam <- ggplot(summary_age_gam, aes(x = method, y = median_partial_deviance, fill = n_sig_p001)) +
  geom_col(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = fig_cfg$gradient_low, high = fig_cfg$gradient_high, name = "p<0.001 loci") +
  labs(title = "Age-GAM Signal Preservation", x = NULL, y = "Median partial deviance") +
  theme_batch() +
  rotate_x()
size_gam <- fig_size("age_gam")
save_plot_png(file.path(fig_dir, "age_gam_preservation_summary.png"), p_age_gam,
              width = size_gam$width, height = size_gam$height)

# --- bridge-pair agreement ---
bridge_results_list <- parallel::mclapply(names(methods), function(method_name) {
  method_results <- list()
  expr_mat <- filter_mat_for_eval(
    methods[[method_name]],
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_eval_samples
  )
  for (i in 1:nrow(bridge_df_for_eval)) {
    sample1_id <- bridge_df_for_eval$original_sample[i]
    sample2_id <- bridge_df_for_eval$batch9_replicate[i]
    vec1 <- expr_mat[, sample1_id]
    vec2 <- expr_mat[, sample2_id]
    method_results[[i]] <- data.frame(
      method      = method_name,
      correlation = cor(vec1, vec2, method = "pearson"),
      ccc         = 2 * cov(vec1, vec2) /
        (var(vec1) + var(vec2) + (mean(vec1) - mean(vec2))^2),
      mae         = mean(abs(vec1 - vec2)),
      mse         = mean((vec1 - vec2)^2),
      rmse        = sqrt(mean((vec1 - vec2)^2))
    )
  }
  dplyr::bind_rows(method_results)
}, mc.cores = eval_workers, mc.preschedule = FALSE)

bridge_results_df <- dplyr::bind_rows(bridge_results_list)
bridge_results_df$method <- factor(bridge_results_df$method,
                                   levels = method_levels)
write_tsv(bridge_results_df, file.path(out_dir, "bridge_pair_metrics.tsv"))

p_corr <- ggplot(bridge_results_df, aes(x = method, y = correlation, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.82, linewidth = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = fig_cfg$text) +
  scale_method_fill(name = "Method") +
  labs(title = "Bridge Pair Correlation", x = NULL, y = "Pearson Correlation") +
  theme_batch() +
  rotate_x() +
  theme(legend.position = "none")

p_mse <- ggplot(bridge_results_df, aes(x = method, y = mse, fill = method)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.82, linewidth = 0.35) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.1, color = fig_cfg$text) +
  scale_method_fill(name = "Method") +
  scale_y_log10() +
  labs(title = "Bridge Pair MSE", x = NULL, y = "MSE (log10 scale)") +
  theme_batch() +
  rotate_x() +
  theme(legend.position = "none")

p_bridge_eval <- (p_corr | p_mse) +
  plot_annotation(title = "Evaluation of Bridge Pair Concordance",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          size = fig_cfg$title_size + 2,
                                                          hjust = fig_cfg$title_hjust %||% 0.5,
                                                          family = fig_cfg$font_family)))

size_bridge <- fig_size("bridge")
save_plot_png(file.path(fig_dir, "bridge_pair_evaluation.png"), p_bridge_eval,
              width = size_bridge$width, height = size_bridge$height)

# --- variance partition ---
vp_formula <- ~ (1|LifeStage) + (1|batch) + (1|sex) + RIN

vp_method_results <- parallel::mclapply(names(methods), function(method_name) {
  expr_mat       <- filter_mat_for_eval(
    methods[[method_name]],
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_eval_samples
  )
  common_samples <- intersect(colnames(expr_mat), rownames(meta_final))
  set.seed(2025)
  gene_subset <- sample(rownames(expr_mat), min(5000, nrow(expr_mat)))

  varPart <- fitExtractVarPartModel(
    expr_mat[gene_subset, common_samples],
    vp_formula,
    meta_final[common_samples, ]
  )

  vp_mean <- colMeans(varPart)
  vp_results <- data.frame(
    method             = method_name,
    variable           = names(vp_mean),
    variance_explained = vp_mean
  )

  varpart_df <- as.data.frame(varPart) %>%
    rownames_to_column("feature") %>%
    mutate(method = method_name, .before = 1)
  write_tsv(varpart_df, file.path(out_dir, paste0("variance_partition_features_", method_name, ".tsv")))

  if ("batch" %in% colnames(varPart)) {
    batch_var <- varPart[, "batch"]
    vp_tail <- data.frame(
      method = method_name,
      Var_batch_median = median(batch_var, na.rm = TRUE),
      Var_batch_p75 = quantile(batch_var, 0.75, na.rm = TRUE),
      Var_batch_p95 = quantile(batch_var, 0.95, na.rm = TRUE),
      n_features_var_batch_gt_0.1 = sum(batch_var > 0.1, na.rm = TRUE),
      n_features_var_batch_gt_0.3 = sum(batch_var > 0.3, na.rm = TRUE)
    )
  } else {
    vp_tail <- data.frame(
      method = method_name,
      Var_batch_median = NA_real_,
      Var_batch_p75 = NA_real_,
      Var_batch_p95 = NA_real_,
      n_features_var_batch_gt_0.1 = NA_integer_,
      n_features_var_batch_gt_0.3 = NA_integer_
    )
  }
  list(results = vp_results, tail = vp_tail)
}, mc.cores = eval_workers, mc.preschedule = FALSE)

vp_results_list <- lapply(vp_method_results, `[[`, "results")
vp_tail_list <- lapply(vp_method_results, `[[`, "tail")

vp_results_df <- do.call(rbind, vp_results_list)
vp_tail_df <- do.call(rbind, vp_tail_list)
write_tsv(vp_tail_df, file.path(out_dir, "variance_partition_tail_summary.tsv"))

vp_results_df$method <- factor(vp_results_df$method,
                               levels = method_levels)

p_vp_eval <- ggplot(vp_results_df, aes(x = variable, y = variance_explained, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.25, width = 0.8) +
  scale_method_fill(name = "Method") +
  labs(
    title = "Variance Explained by Key Factors After Correction",
    x     = "Variable",
    y     = "Mean Proportion of Variance Explained"
  ) +
  theme_batch() +
  rotate_x()

size_vp <- fig_size("variance_partition")
save_plot_png(file.path(fig_dir, "variance_partition_evaluation.png"), p_vp_eval,
              width = size_vp$width, height = size_vp$height)

# --- ASW ---
asw_results_list <- parallel::mclapply(names(methods), function(method_name) {
  expr_mat <- filter_mat_for_eval(
    methods[[method_name]],
    method_name,
    feature_ids = common_eval_genes,
    sample_ids = common_eval_samples
  )
  set.seed(2025)
  variances         <- apply(expr_mat, 1, var)
  expr_mat_filtered <- expr_mat[variances > 1e-8, , drop = FALSE]
  pca_res           <- prcomp(t(expr_mat_filtered), center = TRUE, scale. = TRUE)
  n_pcs             <- min(50, ncol(pca_res$x))
  pca_scores        <- pca_res$x[, 1:n_pcs]
  dist_matrix       <- dist(pca_scores)
  meta_aligned      <- meta_final[rownames(pca_scores), ]

  batch_numeric <- as.numeric(meta_aligned$batch)
  sil_batch     <- silhouette(batch_numeric, dist_matrix)
  asw_batch     <- mean(sil_batch[, "sil_width"])

  lifestage_numeric <- as.numeric(meta_aligned$LifeStage)
  sil_lifestage     <- silhouette(lifestage_numeric, dist_matrix)
  asw_lifestage     <- mean(sil_lifestage[, "sil_width"])

  data.frame(
    method        = method_name,
    ASW_batch     = asw_batch,
    ASW_LifeStage = asw_lifestage
  )
}, mc.cores = eval_workers, mc.preschedule = FALSE)

asw_results_df <- do.call(rbind, asw_results_list)
asw_long_df <- asw_results_df %>%
  pivot_longer(cols = c("ASW_batch", "ASW_LifeStage"),
               names_to = "metric", values_to = "score") %>%
  mutate(
    metric = factor(metric, levels = c("ASW_batch", "ASW_LifeStage")),
    Goal   = ifelse(metric == "ASW_batch", "Batch Effect", "LifeStage Boundary")
  )
asw_long_df$method <- factor(asw_long_df$method,
                             levels = method_levels)

p_asw_eval <- ggplot(asw_long_df, aes(x = method, y = score, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", color = "white", linewidth = 0.25) +
  facet_wrap(~ Goal, scales = "free_x") +
  scale_fill_manual(
    values = asw_palette,
    labels = c("Batch Effect", "LifeStage Boundary")
  ) +
  labs(
    title = "Average Silhouette Width",
    x     = "Correction Method",
    y     = "Average Silhouette Score",
    fill  = "Metric"
  ) +
  theme_batch() +
  rotate_x()
size_asw <- fig_size("asw")
save_plot_png(file.path(fig_dir, "asw_evaluation.png"), p_asw_eval,
              width = size_asw$width, height = size_asw$height)

# --- method selection ---
# Four equally weighted domains; age-derived metrics are diagnostic only.
summary_bridge <- bridge_results_df %>%
  group_by(method) %>%
  summarise(Correlation = mean(correlation, na.rm = TRUE),
            CCC         = mean(ccc,         na.rm = TRUE),
            MAE         = mean(mae,         na.rm = TRUE),
            MSE         = mean(mse,         na.rm = TRUE),
            RMSE        = mean(rmse,        na.rm = TRUE),
            .groups = "drop")

summary_vp <- vp_results_df %>%
  filter(.data$variable %in% c("batch", "LifeStage")) %>%
  mutate(metric = dplyr::recode(.data$variable,
                                batch = "Var_batch",
                                LifeStage = "Var_LifeStage")) %>%
  select(.data$method, .data$metric, .data$variance_explained) %>%
  pivot_wider(names_from = "metric", values_from = "variance_explained")

if (!"Var_batch" %in% names(summary_vp)) summary_vp$Var_batch <- NA_real_
if (!"Var_LifeStage" %in% names(summary_vp)) summary_vp$Var_LifeStage <- NA_real_
summary_vp <- summary_vp %>% select(.data$method, .data$Var_batch, .data$Var_LifeStage)
summary_vp <- summary_vp %>% left_join(vp_tail_df %>% select(method, Var_batch_p95), by = "method")

summary_asw <- asw_results_df

summary_age_concordance_rank <- summary_age_concordance %>%
  select(method, AgeMatched_gap)

summary_age_gam_rank <- summary_age_gam %>%
  select(
    method,
    AgeGAM_sig_p001 = n_sig_p001,
    AgeGAM_median_partial_deviance = median_partial_deviance,
    AgeGAM_median_abs_delta_logit = median_abs_delta_logit,
    AgeGAM_k_index_fail_fraction = k_index_fail_fraction
  )

summary_all <- list(
  summary_bridge,
  summary_vp,
  summary_asw,
  summary_residual,
  summary_pc_batch,
  summary_age_concordance_rank,
  summary_age_gam_rank
) %>%
  reduce(left_join, by = "method")

ranking_metrics <- c(
  "Correlation", "CCC", "MAE", "RMSE",
  "ResidualCoef_median", "ResidualCoef_p95", "ResidualCoef_max_p95",
  "ResidualCoef_gt1_total", "Var_batch", "Var_batch_p95", "ASW_batch",
  "Var_LifeStage", "ASW_LifeStage"
)
higher_is_better <- c("Correlation", "CCC", "Var_LifeStage", "ASW_LifeStage")
ranking_df <- summary_all %>% dplyr::filter(.data$method != "Raw")

for (metric in ranking_metrics) {
  values <- ranking_df[[metric]]
  ranking_df[[paste0(metric, "_rank")]] <- rank(
    if (metric %in% higher_is_better) -values else values,
    ties.method = "average"
  )
}

ranking_df$Bridge_concordance_score <- mean_ranks(c("Correlation", "CCC"))
ranking_df$Bridge_error_score <- mean_ranks(c("MAE", "RMSE"))
ranking_df$Residual_batch_score <- mean_ranks(c(
  "ResidualCoef_median", "ResidualCoef_p95", "ResidualCoef_max_p95",
  "ResidualCoef_gt1_total"
))
ranking_df$Batch_variance_score <- mean_ranks(c("Var_batch", "Var_batch_p95"))
ranking_df$Batch_removal_score <- rowMeans(cbind(
  ranking_df$Residual_batch_score,
  ranking_df$Batch_variance_score,
  ranking_df$ASW_batch_rank
))
ranking_df$Biological_preservation_score <- mean_ranks(c("Var_LifeStage", "ASW_LifeStage"))
domain_score_cols <- c(
  "Bridge_concordance_score", "Bridge_error_score",
  "Batch_removal_score", "Biological_preservation_score"
)
ranking_df$Composite_score <- rowMeans(ranking_df[, domain_score_cols, drop = FALSE])

candidate_ranking <- ranking_df %>%
  dplyr::arrange(.data$Composite_score, .data$method)
selected_method <- candidate_ranking$method[[1]]
write_tsv(
  data.table::as.data.table(methods[[selected_method]], keep.rownames = "site_id"),
  selected_output_file
)
ranking_df$selection <- ifelse(ranking_df$method == selected_method, "top_ranked", "")
ranking_df$selection_rationale <- ifelse(
  ranking_df$method == selected_method,
  paste(
    "Selected as the lowest equally weighted four-domain composite.",
    "Age-derived metrics were excluded from selection."
  ),
  ""
)
rank_cols <- paste0(ranking_metrics, "_rank")
ranking_df <- ranking_df %>%
  dplyr::select(
    .data$method,
    dplyr::all_of(rank_cols),
    .data$Residual_batch_score,
    .data$Batch_variance_score,
    dplyr::all_of(domain_score_cols),
    .data$Composite_score,
    .data$selection,
    .data$selection_rationale
  )
write_tsv(summary_all, file.path(out_dir, "method_summary_metrics.tsv"))
write_tsv(ranking_df, file.path(out_dir, "method_domain_composite.tsv"))

choice_scores <- ranking_df %>%
  dplyr::select(
    .data$method,
    .data$Bridge_concordance_score,
    .data$Bridge_error_score,
    .data$Batch_removal_score,
    .data$Biological_preservation_score,
    .data$Composite_score,
    .data$selection
  ) %>%
  dplyr::arrange(.data$Composite_score)
write_tsv(choice_scores, file.path(out_dir, "method_choice_scores.tsv"))

method_labels <- choice_scores %>%
  dplyr::mutate(
    method_label = .data$method,
    method_label = factor(.data$method_label, levels = rev(.data$method_label)),
    method_role = dplyr::case_when(
      .data$selection == "top_ranked" ~ "Selected",
      TRUE ~ "Candidate"
    )
  )

choice_long <- method_labels %>%
  dplyr::select(
    .data$method,
    .data$method_label,
    dplyr::all_of(domain_score_cols)
  ) %>%
  tidyr::pivot_longer(
    cols = -c(.data$method, .data$method_label),
    names_to = "criterion",
    values_to = "domain_score"
  ) %>%
  dplyr::mutate(
    criterion = factor(
      .data$criterion,
      levels = c(
        "Bridge_concordance_score",
        "Bridge_error_score",
        "Batch_removal_score",
        "Biological_preservation_score"
      ),
      labels = c(
        "Bridge concordance",
        "Bridge error",
        "Batch removal",
        "Biological preservation"
      )
    )
  ) %>%
  dplyr::mutate(domain_contribution = .data$domain_score / length(domain_score_cols))

p_method_choice <- ggplot(choice_long, aes(x = domain_contribution, y = method_label, fill = criterion)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25,
           position = position_stack(reverse = TRUE)) +
  geom_col(
    data = method_labels,
    aes(x = .data$Composite_score, y = .data$method_label, color = .data$method_role),
    inherit.aes = FALSE,
    fill = NA,
    width = 0.72,
    linewidth = 1.0
  ) +
  geom_text(
    data = method_labels,
    aes(x = .data$Composite_score + 0.2, y = .data$method_label,
        label = sprintf("%.2f", .data$Composite_score)),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5.0,
    color = fig_cfg$text
  ) +
  scale_fill_manual(
    values = c(
      "Bridge concordance" = "#009E73",
      "Bridge error" = "#56B4E9",
      "Batch removal" = "#D55E00",
      "Biological preservation" = "#CC79A7"
    ),
    name = "Domain"
  ) +
  scale_color_manual(
    values = c(Selected = "#000000", Candidate = "transparent"),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, max(method_labels$Composite_score, na.rm = TRUE) + 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Batch-Correction Method Selection by Domain Composite",
    x = "Composite domain-rank score (lower is better)",
    y = NULL
  ) +
  theme_batch(base_size = max(fig_cfg$base_size + 5, 15)) +
  theme(
    plot.title = element_text(face = "bold", size = fig_cfg$title_size + 5,
                              hjust = fig_cfg$title_hjust %||% 0.5,
                              color = fig_cfg$text),
    axis.text.y = element_text(size = fig_cfg$tick_size + 6, color = fig_cfg$text),
    axis.text.x = element_text(size = fig_cfg$tick_size + 5, color = fig_cfg$text),
    axis.title.x = element_text(size = fig_cfg$label_size + 6, color = fig_cfg$text,
                                margin = margin(t = 8)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = fig_cfg$tick_size + 5),
    legend.text = element_text(size = fig_cfg$tick_size + 5),
    legend.key.size = unit(18, "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

size_choice <- fig_size("method_choice")
save_plot_png(file.path(fig_dir, "method_choice_evidence.png"), p_method_choice,
              width = size_choice$width, height = size_choice$height)

write_tsv(
  dplyr::bind_rows(plot_validation_records),
  file.path(out_dir, "layout_validation.tsv")
)
