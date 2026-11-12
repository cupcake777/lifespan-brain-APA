#!/usr/bin/env Rscript
# Batch correction of expression for eQTL mapping.
# Sex QC, outlier removal, seven correction methods scored on bridge pairs and variance partition.
# Writes QTL-ready phenotype, covariate and sample files per life stage.

for (.apa_v in c("APA_HPC_HOME")) if (!nzchar(Sys.getenv(.apa_v))) stop("Set ", .apa_v, " (source config/paths.sh)")

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(sva)
  library(dplyr)
  library(ruv)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(variancePartition)
  library(PCAForQTL)
  library(RUVSeq)
  library(splines)
})

# --- helpers ---

log_message <- function(message) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), message))
}

write_tsv <- function(df, filename, row_names = FALSE) {
  data.table::fwrite(df, filename, sep = "\t", quote = FALSE, na = "NA", row.names = row_names)
}

load_matrix <- function(filepath) {
  df <- data.table::fread(filepath, data.table = FALSE, check.names = FALSE)
  rownames(df) <- df[[1]]
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

# k with the largest distance from the line joining the first and last points
select_k <- function(df_k) {
  if (max(df_k$batch, na.rm = TRUE) == min(df_k$batch, na.rm = TRUE)) {
    return(df_k$k[1])
  }
  df_norm <- data.frame(
    k = (df_k$k - min(df_k$k)) / (max(df_k$k) - min(df_k$k)),
    batch = (df_k$batch - min(df_k$batch, na.rm = TRUE)) /
      (max(df_k$batch, na.rm = TRUE) - min(df_k$batch, na.rm = TRUE))
  )
  p1 <- as.numeric(df_norm[1, ])
  p2 <- as.numeric(df_norm[nrow(df_norm), ])
  d <- sapply(seq_len(nrow(df_norm)), function(i) {
    p <- as.numeric(df_norm[i, ])
    v1 <- p2 - p1
    v2 <- p - p1
    abs(v1[1] * v2[2] - v1[2] * v2[1]) / sqrt(sum(v1^2))
  })
  df_k$k[which.max(d)]
}

# --- config ---

root <- Sys.getenv("APA_ROOT")
raw_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project", "raw_data")

count_file <- file.path(root, "raw", "all.912.expr.count")
tpm_file <- file.path(root, "raw", "all.912.expr.tpm")
meta_file <- file.path(raw_dir, "all.827.meta")
gtf_file <- file.path(root, "ref", "gencode.v40.annotation.gtf.gz")
geno_prefix <- file.path(raw_dir, "genotype")

out_dir <- file.path(root, "data", "04_qtl")
fig_dir <- file.path(root, "figures", "04_qtl")
bc_dir <- file.path(out_dir, "eqtl", "batch_correct")
if (!dir.exists(bc_dir)) dir.create(bc_dir, recursive = TRUE)

rin_min <- 5
tpm_min <- 0.1
count_min <- 6
expressed_fraction <- 0.20

primary_batch <- "batch9"
control_batch_response_min_q <- 0.50
control_bio_inert_max_q <- 0.40
control_min_n <- 300
control_age_df <- 3
ruv_k_range <- 0:10
test_methods <- c(
  "ComBat_RAW",
  "ComBat_REF",
  "ComBat_REF_NP",
  "ComBat_REF_meanonly",
  "limma",
  "RUV-III",
  "RUVs"
)
method_family <- c(
  Raw = "raw",
  ComBat_RAW = "ComBat",
  ComBat_REF = "ComBat",
  ComBat_REF_NP = "ComBat",
  ComBat_REF_meanonly = "ComBat",
  limma = "limma",
  `RUV-III` = "RUV",
  RUVs = "RUV"
)

vp_n_genes <- 5000
pca_n_genes <- 5000
geno_pca_n <- 5
geno_pca_maf <- 0.05
geno_pca_indep_window <- 50
geno_pca_indep_step <- 5
geno_pca_indep_r2 <- 0.2

method_palette <- c(
  Raw = "#4D4D4D",
  ComBat_RAW = "#E69F00",
  ComBat_REF = "#009E73",
  ComBat_REF_NP = "#56B4E9",
  ComBat_REF_meanonly = "#CC79A7",
  limma = "#0072B2",
  `RUV-III` = "#F0E442",
  RUVs = "#D55E00"
)

stage_palette <- c(Prenatal = "#E69F00", Postnatal = "#0072B2")
stage_shapes <- c(Prenatal = 21, Postnatal = 22)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- loading expression and annotation ---

log_message("Loading count, TPM, metadata, and GTF annotation.")

count <- load_matrix(count_file)
tpm <- load_matrix(tpm_file)

meta <- data.table::fread(meta_file, data.table = FALSE, check.names = FALSE)
rownames(meta) <- meta[[1]]
meta <- meta[, -1, drop = FALSE]

gtf <- data.table::fread(cmd = sprintf("zcat %s | grep -v '^#'", shQuote(gtf_file)),
                         sep = "\t", header = FALSE, data.table = FALSE, quote = "")
colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end",
                        "score", "strand", "frame", "attribute")
gtf <- gtf[gtf$feature == "gene", , drop = FALSE]

annotation <- data.frame(
  chr = gtf$chr,
  start = as.integer(gtf$start),
  end = as.integer(gtf$end),
  strand = gtf$strand,
  gene_id = sub('.*gene_id "([^"]+)".*', "\\1", gtf$attribute),
  gene_name = ifelse(grepl('gene_name "', gtf$attribute),
                     sub('.*gene_name "([^"]+)".*', "\\1", gtf$attribute),
                     NA_character_),
  gene_type = ifelse(grepl('gene_type "', gtf$attribute),
                     sub('.*gene_type "([^"]+)".*', "\\1", gtf$attribute),
                     NA_character_),
  stringsAsFactors = FALSE
)
annotation$gene_id_base <- sub("\\.[0-9]+$", "", annotation$gene_id)
annotation <- annotation[!duplicated(annotation$gene_id_base), ]

# --- sample and gene alignment ---

meta <- meta[meta$RIN >= rin_min, , drop = FALSE]

common_samples <- Reduce(intersect, list(rownames(meta), colnames(count), colnames(tpm)))
meta <- meta[common_samples, , drop = FALSE]
count <- count[, common_samples, drop = FALSE]
tpm <- tpm[, common_samples, drop = FALSE]

common_genes <- intersect(rownames(count), rownames(tpm))
count <- count[common_genes, , drop = FALSE]
tpm <- tpm[common_genes, , drop = FALSE]

log_message(sprintf("After RIN/sample/gene alignment: %d genes x %d samples.",
                    nrow(count), ncol(count)))

# --- sex validation ---

set.seed(2026)

xist_id <- grep("^ENSG00000229807", rownames(count), value = TRUE)[1]
uty_id <- grep("^ENSG00000183878", rownames(count), value = TRUE)[1]

xist_expr <- setNames(as.numeric(count[xist_id, ]), colnames(count))
uty_expr <- setNames(as.numeric(count[uty_id, ]), colnames(count))

xist_threshold <- expm1(mean(kmeans(log1p(xist_expr), centers = 2)$centers))
uty_threshold <- expm1(mean(kmeans(log1p(uty_expr), centers = 2)$centers))

predicted <- rep("Ambiguous", nrow(meta))
predicted[xist_expr[rownames(meta)] > xist_threshold &
            uty_expr[rownames(meta)] <= uty_threshold] <- "Female"
predicted[xist_expr[rownames(meta)] <= xist_threshold &
            uty_expr[rownames(meta)] > uty_threshold] <- "Male"

reported <- as.character(meta$sex)
reported[grepl("^[Ff]", reported)] <- "Female"
reported[grepl("^[Mm]", reported)] <- "Male"

unknown <- is.na(reported) | !reported %in% c("Female", "Male")
ambiguous <- predicted == "Ambiguous"
mismatch <- !unknown & reported != predicted & !ambiguous
status <- ifelse(ambiguous, "ambiguous", ifelse(mismatch, "mismatch", "kept"))

sex_df <- data.frame(
  sample = rownames(meta),
  XIST = xist_expr[rownames(meta)],
  UTY = uty_expr[rownames(meta)],
  reported_sex = reported,
  predicted_sex = predicted,
  status = factor(status, levels = c("kept", "mismatch", "ambiguous")),
  stringsAsFactors = FALSE
)

write_tsv(sex_df, file.path(out_dir, "sex_validation.tsv"))

df_plot <- sex_df[order(sex_df$status), ]
p_sex <- ggplot(df_plot, aes(x = log1p(XIST), y = log1p(UTY))) +
  geom_vline(xintercept = log1p(xist_threshold), linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = log1p(uty_threshold), linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_point(aes(color = predicted_sex, shape = status, size = status, alpha = status),
             stroke = 0.6) +
  geom_text_repel(
    data = subset(df_plot, status != "kept"),
    aes(label = sample, color = predicted_sex),
    size = 2.8,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.12,
    box.padding = 0.6,
    point.padding = 0.4,
    force = 6,
    force_pull = 0.5,
    max.overlaps = Inf,
    segment.color = "grey40",
    segment.size = 0.35,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Predicted sex",
    values = c(Female = "#D55E00", Male = "#0072B2", Ambiguous = "#7F7F7F")
  ) +
  scale_shape_manual(
    name = "QC status",
    values = c(kept = 16, mismatch = 17, ambiguous = 8)
  ) +
  scale_size_manual(
    name = "QC status",
    values = c(kept = 1.6, mismatch = 3.2, ambiguous = 3.5)
  ) +
  scale_alpha_manual(
    name = "QC status",
    values = c(kept = 0.45, mismatch = 0.95, ambiguous = 0.95)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 9.5),
    axis.title = element_text(color = "black", size = 10.5, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.spacing.y = grid::unit(0.2, "cm"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0)
  ) +
  labs(
    title = "Sample sex concordance validation",
    x = expression(bold(log(italic(XIST)~count + 1))),
    y = expression(bold(log(italic(UTY)~count + 1)))
  )

ggsave(file.path(fig_dir, "sex_validation.png"), p_sex,
       width = 8.5, height = 5.8, dpi = 300)

drop_samples <- sex_df$sample[sex_df$status != "kept"]
meta <- meta[setdiff(rownames(meta), drop_samples), , drop = FALSE]
count <- count[, rownames(meta), drop = FALSE]
tpm <- tpm[, rownames(meta), drop = FALSE]

log_message(sprintf("Sex QC removed %d samples.", length(drop_samples)))

# --- expression filtering and TMM ---

threshold_n <- ceiling(expressed_fraction * ncol(count))
gene_id_base <- sub("\\.[0-9]+$", "", rownames(count))
autosomal_ids <- annotation$gene_id_base[annotation$chr %in% paste0("chr", 1:22)]

keep_gene <- rowSums(count) > 0 &
  rowSums(tpm > tpm_min) >= threshold_n &
  rowSums(count >= count_min) >= threshold_n &
  gene_id_base %in% autosomal_ids

count <- count[keep_gene, , drop = FALSE]
tpm <- tpm[keep_gene, , drop = FALSE]

dge <- edgeR::DGEList(count)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
logcpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)

log_message(sprintf("After expression filtering: %d genes x %d samples.",
                    nrow(logcpm), ncol(logcpm)))

# --- sample outlier detection ---

gene_mad <- apply(as.matrix(logcpm), 1, mad, na.rm = TRUE)
keep_genes <- names(sort(gene_mad, decreasing = TRUE))[1:min(8000, length(gene_mad))]
logcpm_sub <- logcpm[keep_genes, , drop = FALSE]

outlier_list <- list()

for (stage in unique(meta$LifeStage)) {
  stage_samples <- intersect(rownames(meta)[meta$LifeStage == stage], colnames(logcpm_sub))
  if (length(stage_samples) < 3) next

  sample_cor <- cor(as.matrix(logcpm_sub[, stage_samples, drop = FALSE]),
                    use = "pairwise.complete.obs")
  adj <- 0.5 * (1 + sample_cor)
  conn <- rowSums(adj) - 1
  z_conn <- as.numeric(scale(conn))

  outlier_list[[stage]] <- data.frame(
    sample = stage_samples,
    LifeStage = stage,
    connectivity_z = z_conn,
    is_outlier = z_conn < -3,
    stringsAsFactors = FALSE
  )
}

outliers <- bind_rows(outlier_list)
write_tsv(outliers, file.path(out_dir, "sample_outlier_connectivity.tsv"))

outliers <- outliers[order(outliers$LifeStage), ]
outliers$sample_index <- seq_len(nrow(outliers))

p_outlier <- ggplot(outliers, aes(x = sample_index, y = connectivity_z, color = LifeStage)) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "grey40", linewidth = 0.45) +
  geom_point(size = 3.6, alpha = 0.8) +
  geom_text_repel(
    data = subset(outliers, is_outlier),
    aes(label = sample),
    size = 2.8,
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.12,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 50,
    force = 4,
    show.legend = FALSE,
    min.segment.length = 0,
    segment.color = "grey40",
    segment.size = 0.4,
    segment.alpha = 0.8
  ) +
  scale_color_manual(values = stage_palette) +
  scale_y_continuous(breaks = seq(-8, 3, 2)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0)
  ) +
  labs(
    title = "Sample outlier detection",
    x = "Sample index",
    y = "Connectivity Z-score"
  )

ggsave(file.path(fig_dir, "sample_outlier_connectivity.png"), p_outlier,
       width = 8, height = 4.8, dpi = 300)

drop_outliers <- outliers$sample[outliers$is_outlier]
keep_samples <- setdiff(colnames(logcpm), drop_outliers)

meta <- meta[keep_samples, , drop = FALSE]
count <- count[, keep_samples, drop = FALSE]
logcpm <- logcpm[, keep_samples, drop = FALSE]

meta$batch <- factor(meta$batch)
meta$LifeStage <- factor(meta$LifeStage)
meta$sex <- factor(meta$sex)
meta$RIN <- as.numeric(meta$RIN)
meta$age <- as.numeric(meta$age)

log_message(sprintf("Outlier QC removed %d samples; retained %d genes x %d samples.",
                    length(drop_outliers), nrow(logcpm), ncol(logcpm)))

# --- bridge pairs ---

batch9_samples <- rownames(meta)[meta$batch == primary_batch]
batch9_base <- sub("[-_.](rep|replicate)$", "", batch9_samples, ignore.case = TRUE)

pair_list <- list()

for (batch_id in setdiff(unique(as.character(meta$batch)), primary_batch)) {
  if (is.na(batch_id)) next

  batch_samples <- rownames(meta)[meta$batch == batch_id]
  batch_base <- sub("[-_.](rep|replicate)$", "", batch_samples, ignore.case = TRUE)
  common_base <- intersect(batch_base, batch9_base)

  if (length(common_base) == 0) next

  pair_list[[batch_id]] <- data.frame(
    batch = batch_id,
    original_sample = batch_samples[match(common_base, batch_base)],
    batch9_replicate = batch9_samples[match(common_base, batch9_base)],
    stringsAsFactors = FALSE
  )
}

bridge_pairs <- bind_rows(pair_list)
write_tsv(bridge_pairs, file.path(out_dir, "bridge_pairs.tsv"))

log_message(sprintf("Bridge pairs: %d.", nrow(bridge_pairs)))

# --- RUV negative-control genes ---

stable_genes <- character()

if (nrow(bridge_pairs) > 0) {
  mat <- as.matrix(logcpm)
  d <- as.matrix(mat[, bridge_pairs$original_sample, drop = FALSE]) -
    as.matrix(mat[, bridge_pairs$batch9_replicate, drop = FALSE])
  batches <- unique(bridge_pairs$batch)

  delta_b <- vapply(batches, function(b) {
    rowMeans(d[, bridge_pairs$batch == b, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(d)))

  ss <- rep(0, nrow(d))
  dfn <- rep(0, nrow(d))

  for (b in batches) {
    dd <- d[, bridge_pairs$batch == b, drop = FALSE]
    if (ncol(dd) < 2) next
    ss <- ss + rowSums((dd - rowMeans(dd, na.rm = TRUE))^2, na.rm = TRUE)
    dfn <- dfn + pmax(rowSums(is.finite(dd)) - 1, 0)
  }

  noise_var <- ss / pmax(dfn, 1)
  signal <- apply(delta_b, 1, var, na.rm = TRUE)
  n_bar <- nrow(bridge_pairs) / length(batches)

  df_batch_res <- data.frame(
    gene_id = rownames(d),
    batch_signal = signal,
    batch_noise = noise_var,
    batch_F = signal / pmax(noise_var / n_bar, .Machine$double.eps),
    stringsAsFactors = FALSE
  )

  batch9_ids <- rownames(meta)[meta$batch == primary_batch]
  md_batch9 <- droplevels(meta[batch9_ids, , drop = FALSE])

  y_batch9 <- as.matrix(mat[, batch9_ids, drop = FALSE])
  age_basis <- splines::ns(log10(pmax(md_batch9$age + 1, 1e-6)), df = control_age_df)
  colnames(age_basis) <- paste0("age_ns", seq_len(ncol(age_basis)))
  md_age <- cbind(md_batch9, as.data.frame(age_basis))

  d_full <- model.matrix(
    as.formula(paste("~", paste(c(colnames(age_basis), "sex", "RIN"), collapse = " + "))),
    data = md_age
  )
  d_red <- model.matrix(~ sex + RIN, data = md_age)

  fit_full <- limma::lmFit(y_batch9, d_full)
  fit_red <- limma::lmFit(y_batch9, d_red)
  ss_full <- fit_full$sigma^2 * fit_full$df.residual
  ss_red <- fit_red$sigma^2 * fit_red$df.residual

  df_bio_inert <- data.frame(
    gene_id = rownames(y_batch9),
    bio_partial_r2 = pmax((ss_red - ss_full) / pmax(ss_red, .Machine$double.eps), 0),
    stringsAsFactors = FALSE
  )

  control_tab <- merge(df_batch_res, df_bio_inert, by = "gene_id", sort = FALSE)
  q_batch <- quantile(control_tab$batch_F, control_batch_response_min_q, na.rm = TRUE)
  control_tab$batch_responsive <- is.finite(control_tab$batch_F) &
    control_tab$batch_F >= q_batch

  relaxed <- control_bio_inert_max_q

  repeat {
    q_bio <- quantile(control_tab$bio_partial_r2, relaxed, na.rm = TRUE)
    control_tab$bio_inert <- is.finite(control_tab$bio_partial_r2) &
      control_tab$bio_partial_r2 <= q_bio
    control_tab$is_control <- control_tab$batch_responsive & control_tab$bio_inert

    if (sum(control_tab$is_control) >= control_min_n || relaxed >= 0.9) break

    relaxed <- relaxed + 0.1
    log_message(sprintf("Control genes: relaxing biological-inertness quantile to %.2f.",
                        relaxed))
  }

  control_tab$mean_logcpm <- rowMeans(mat[control_tab$gene_id, , drop = FALSE],
                                      na.rm = TRUE)
  stable_genes <- control_tab$gene_id[control_tab$is_control]
  write_tsv(data.frame(gene_id = stable_genes),
            file.path(out_dir, "ruv_control_genes.tsv"))
  write_tsv(control_tab, file.path(out_dir, "ruv_control_selection.tsv"))
}

log_message(sprintf("RUV negative-control genes: %d.", length(stable_genes)))

# --- seven correction tests ---

combat_mod <- model.matrix(~ LifeStage, data = meta)

combat_raw <- sva::ComBat(
  dat = as.matrix(logcpm),
  batch = meta$batch,
  mod = combat_mod,
  par.prior = TRUE
)

combat_ref <- sva::ComBat(
  dat = as.matrix(logcpm),
  batch = meta$batch,
  mod = combat_mod,
  ref.batch = primary_batch,
  par.prior = TRUE
)

combat_ref_np <- sva::ComBat(
  dat = as.matrix(logcpm),
  batch = meta$batch,
  mod = combat_mod,
  ref.batch = primary_batch,
  par.prior = FALSE
)

combat_ref_meanonly <- sva::ComBat(
  dat = as.matrix(logcpm),
  batch = meta$batch,
  mod = combat_mod,
  ref.batch = primary_batch,
  par.prior = TRUE,
  mean.only = TRUE
)

count_limma <- limma::removeBatchEffect(
  logcpm,
  batch = meta$batch,
  design = combat_mod
)

log_message(sprintf("Testing correction methods: %s", paste(test_methods, collapse = ", ")))

Y <- t(logcpm)
ctl <- rownames(logcpm) %in% stable_genes
M <- ruv::replicate.matrix(gsub("_rep$", "", rownames(meta)))

sc_idx <- as.matrix(data.frame(
  original = match(bridge_pairs$original_sample, colnames(logcpm)),
  replicate = match(bridge_pairs$batch9_replicate, colnames(logcpm))
))
ctl_character <- rownames(logcpm)[rownames(logcpm) %in% stable_genes]

ruv_k_rows <- list()

for (k_val in ruv_k_range) {
  if (k_val == 0) {
    count_corr <- as.matrix(logcpm)
  } else {
    count_corr <- t(RUVIII(Y = Y, M = M, ctl = ctl, k = k_val))
  }

  set.seed(2026)
  gene_subset <- sample(rownames(count_corr), min(2000, nrow(count_corr)))
  varpart <- variancePartition::fitExtractVarPartModel(
    count_corr[gene_subset, rownames(meta), drop = FALSE],
    ~ (1 | LifeStage) + (1 | batch) + (1 | sex) + RIN,
    meta
  )

  ruv_k_rows[[length(ruv_k_rows) + 1]] <- data.frame(
    method = "RUV-III",
    k = k_val,
    t(colMeans(varpart)),
    stringsAsFactors = FALSE
  )
}

for (k_val in ruv_k_range) {
  if (k_val == 0) {
    count_corr <- as.matrix(logcpm)
    meta_ruv <- meta
  } else {
    set_corr <- RUVSeq::RUVs(as.matrix(logcpm), cIdx = ctl_character,
                             k = k_val, scIdx = sc_idx, isLog = TRUE)
    count_corr <- as.matrix(set_corr$normalizedCounts)
    meta_ruv <- cbind(meta, set_corr$W)
  }

  fixed_effects <- "RIN"
  if (k_val > 0) {
    fixed_effects <- paste(c("RIN", paste0("W_", seq_len(k_val))), collapse = " + ")
  }

  set.seed(2026)
  gene_subset <- sample(rownames(count_corr), min(2000, nrow(count_corr)))
  varpart <- variancePartition::fitExtractVarPartModel(
    count_corr[gene_subset, rownames(meta), drop = FALSE],
    as.formula(paste("~ (1 | LifeStage) + (1 | batch) + (1 | sex) +", fixed_effects)),
    meta_ruv
  )

  ruv_k_rows[[length(ruv_k_rows) + 1]] <- data.frame(
    method = "RUVs",
    k = k_val,
    t(colMeans(varpart)),
    stringsAsFactors = FALSE
  )
}

ruv_k_table <- bind_rows(ruv_k_rows)

k_opt_ruviii <- select_k(ruv_k_table[ruv_k_table$method == "RUV-III", ])
k_opt_ruvs <- select_k(ruv_k_table[ruv_k_table$method == "RUVs", ])

log_message(sprintf("RUV k selected: RUV-III=%d, RUVs=%d.", k_opt_ruviii, k_opt_ruvs))

if (k_opt_ruviii == 0) {
  count_ruv_iii <- as.matrix(logcpm)
} else {
  count_ruv_iii <- t(RUVIII(Y = Y, M = M, ctl = ctl, k = k_opt_ruviii))
}

if (k_opt_ruvs == 0) {
  count_ruvs <- as.matrix(logcpm)
} else {
  set_corr_final <- RUVSeq::RUVs(as.matrix(logcpm), cIdx = ctl_character,
                                 k = k_opt_ruvs, scIdx = sc_idx, isLog = TRUE)
  count_ruvs <- as.matrix(set_corr_final$normalizedCounts)
}

methods <- list(
  Raw = logcpm,
  ComBat_RAW = combat_raw,
  ComBat_REF = combat_ref,
  ComBat_REF_NP = combat_ref_np,
  ComBat_REF_meanonly = combat_ref_meanonly,
  limma = count_limma,
  `RUV-III` = count_ruv_iii,
  RUVs = count_ruvs
)

# --- evaluation genes and samples ---

raw_var <- apply(logcpm, 1, var)
eval_genes <- rownames(logcpm)[raw_var > 1e-8]
eval_samples <- rownames(meta)

log_message(sprintf("Evaluation set: %d genes x %d samples.",
                    length(eval_genes), length(eval_samples)))

# --- bridge-pair metrics ---

bridge_summary <- data.frame(method = character(), correlation = numeric(), mse = numeric())

if (nrow(bridge_pairs) > 0) {
  bridge_rows <- list()

  for (method_name in names(methods)) {
    mat <- as.matrix(methods[[method_name]][eval_genes, eval_samples, drop = FALSE])

    for (i in seq_len(nrow(bridge_pairs))) {
      original_sample <- bridge_pairs$original_sample[i]
      batch9_replicate <- bridge_pairs$batch9_replicate[i]

      bridge_rows[[length(bridge_rows) + 1]] <- data.frame(
        method = method_name,
        correlation = cor(mat[, original_sample], mat[, batch9_replicate]),
        mse = mean((mat[, original_sample] - mat[, batch9_replicate])^2),
        stringsAsFactors = FALSE
      )
    }
  }

  bridge_summary <- bind_rows(bridge_rows)
}

# --- variance partition ---

set.seed(2026)
vp_genes <- sample(eval_genes, min(vp_n_genes, length(eval_genes)))
vp_formula <- ~ (1 | LifeStage) + (1 | batch) + (1 | sex) + RIN
vp_rows <- list()

for (method_name in names(methods)) {
  log_message(sprintf("Variance partition for %s.", method_name))
  mat <- as.matrix(methods[[method_name]][vp_genes, eval_samples, drop = FALSE])

  varpart <- variancePartition::fitExtractVarPartModel(
    mat,
    vp_formula,
    meta[eval_samples, , drop = FALSE]
  )

  vp_rows[[method_name]] <- data.frame(
    method = method_name,
    variable = names(colMeans(varpart)),
    variance = colMeans(varpart),
    stringsAsFactors = FALSE
  )
}

vp_summary <- bind_rows(vp_rows)

# --- pick method from seven tests ---

vp_rank <- vp_summary |>
  filter(variable %in% c("batch", "LifeStage")) |>
  mutate(variable = paste0("var_", variable)) |>
  pivot_wider(names_from = variable, values_from = variance)

if (nrow(bridge_summary) > 0) {
  bridge_rank <- bridge_summary |>
    group_by(method) |>
    summarise(
      bridge_correlation = mean(correlation, na.rm = TRUE),
      bridge_mse = mean(mse, na.rm = TRUE),
      .groups = "drop"
    )

  method_score <- left_join(bridge_rank, vp_rank, by = "method") |>
    mutate(
      rank_bridge_cor = dense_rank(dplyr::desc(bridge_correlation)),
      rank_bridge_mse = dense_rank(bridge_mse),
      rank_batch = dense_rank(var_batch),
      rank_stage = dense_rank(dplyr::desc(var_LifeStage)),
      total_rank = rowMeans(
        cbind(rank_bridge_cor, rank_bridge_mse, rank_batch, rank_stage),
        na.rm = TRUE
      )
    ) |>
    arrange(method == "Raw", total_rank, match(method, test_methods))
} else {
  method_score <- vp_rank |>
    mutate(
      bridge_correlation = NA_real_,
      bridge_mse = NA_real_,
      rank_bridge_cor = NA_integer_,
      rank_bridge_mse = NA_integer_,
      rank_batch = dense_rank(var_batch),
      rank_stage = dense_rank(dplyr::desc(var_LifeStage)),
      total_rank = rowMeans(cbind(rank_batch, rank_stage), na.rm = TRUE)
    ) |>
    arrange(method == "Raw", total_rank, match(method, test_methods))
}

method_score$selection <- ""
selected_method <- method_score$method[method_score$method %in% test_methods][1]
method_score$selection[method_score$method == selected_method] <- "selected"
method_score$family <- unname(method_family[method_score$method])
method_score$k <- NA_integer_
method_score$k[method_score$method == "RUV-III"] <- k_opt_ruviii
method_score$k[method_score$method == "RUVs"] <- k_opt_ruvs
method_score <- method_score |>
  select(method, family, k, everything())

write_tsv(method_score, file.path(out_dir, "method_selection.tsv"))

# --- method plots ---

set.seed(2026)
pca_genes <- sample(eval_genes, min(pca_n_genes, length(eval_genes)))

all_mat <- do.call(cbind, lapply(names(methods), function(method_name) {
  mat <- as.matrix(methods[[method_name]][pca_genes, eval_samples, drop = FALSE])
  colnames(mat) <- paste0(colnames(mat), "___", method_name)
  mat
}))

method_pca <- prcomp(t(all_mat), center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(method_pca$x[, 1:2, drop = FALSE]) |>
  rownames_to_column("sample_method") |>
  tidyr::separate(sample_method, into = c("sample", "method"), sep = "___") |>
  left_join(meta |> rownames_to_column("sample"), by = "sample")

var_exp <- summary(method_pca)$importance[2, 1:2] * 100

p_method_pca <- ggplot(pca_df, aes(PC1, PC2, fill = batch, shape = LifeStage)) +
  geom_point(size = 2.2, alpha = 0.8, color = "black", stroke = 0.25) +
  facet_wrap(~ method, ncol = 2) +
  scale_shape_manual(values = stage_shapes) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom") +
  labs(
    title = "Batch-correction method PCA",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  )

ggsave(file.path(fig_dir, "method_pca.png"), p_method_pca,
       width = 8.5, height = 7, dpi = 300)

p_vp <- vp_summary |>
  filter(variable %in% c("batch", "LifeStage")) |>
  mutate(method = factor(method, levels = c("Raw", test_methods))) |>
  ggplot(aes(method, variance, fill = variable)) +
  geom_col(position = "dodge", linewidth = 0.25, color = "white") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(title = "Variance partition", x = NULL, y = "Mean variance fraction")

if (nrow(bridge_summary) > 0) {
  p_bridge <- bridge_summary |>
    mutate(method = factor(method, levels = c("Raw", test_methods))) |>
    ggplot(aes(method, correlation, fill = method)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 1.1, alpha = 0.45) +
    scale_fill_manual(values = method_palette[names(method_palette) %in% names(methods)]) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(title = "Bridge pair concordance", x = NULL, y = "Pearson correlation")

  ggsave(file.path(fig_dir, "method_summary.png"), p_bridge | p_vp,
         width = 10, height = 4.8, dpi = 300)
} else {
  ggsave(file.path(fig_dir, "method_summary.png"), p_vp,
         width = 6.5, height = 4.8, dpi = 300)
}

# --- selected matrix ---

selected_expr <- methods[[selected_method]]

write_tsv(
  data.frame(gene_id = rownames(selected_expr), selected_expr, check.names = FALSE),
  file.path(bc_dir, paste0("expression.", selected_method, ".tsv"))
)

# --- qtl-ready files ---

gene_id_base <- sub("\\.[0-9]+$", "", rownames(selected_expr))
coords <- annotation[match(gene_id_base, annotation$gene_id_base), ]
keep_coord <- !is.na(coords$chr)

selected_expr <- selected_expr[keep_coord, , drop = FALSE]
coords <- coords[keep_coord, , drop = FALSE]
rownames(coords) <- rownames(selected_expr)

genotype_samples <- data.table::fread(paste0(geno_prefix, ".psam"), data.table = FALSE)
genotype_samples <- as.character(genotype_samples[["#IID"]])

for (stage in unique(meta$LifeStage)) {
  stage_name <- as.character(stage)
  stage_id <- tolower(stage_name)
  log_message(sprintf("Writing QTL inputs for %s.", stage_name))

  stage_samples <- rownames(meta)[meta$LifeStage == stage]
  stage_meta <- meta[stage_samples, , drop = FALSE]
  stage_expr <- selected_expr[, stage_samples, drop = FALSE]

  stage_expr_int <- t(apply(as.matrix(stage_expr), 1, function(x) {
    qnorm(rank(x, ties.method = "average") / (length(x) + 1))
  }))
  rownames(stage_expr_int) <- rownames(stage_expr)
  colnames(stage_expr_int) <- colnames(stage_expr)

  vars <- apply(stage_expr_int, 1, var)
  expr_for_pca <- stage_expr_int[vars > 1e-10, , drop = FALSE]
  expr_pca <- prcomp(t(expr_for_pca), center = TRUE, scale. = TRUE)

  k_be <- as.integer(PCAForQTL::runBE(t(expr_for_pca), B = 20)$numOfPCsChosen)
  k_final <- min(max(k_be, 10), 50, ncol(expr_pca$x))

  known_covs <- data.frame(
    RIN = as.numeric(stage_meta$RIN),
    sex = as.numeric(stage_meta$sex),
    age = as.numeric(stage_meta$age),
    row.names = rownames(stage_meta)
  )

  geno_samples <- intersect(rownames(known_covs), genotype_samples)
  known_covs <- known_covs[geno_samples, , drop = FALSE]
  stage_expr_int <- stage_expr_int[, geno_samples, drop = FALSE]
  expr_pca$x <- expr_pca$x[geno_samples, , drop = FALSE]

  geno_out_prefix <- file.path(bc_dir, paste0("genotype.", stage_id))
  keep_file <- paste0(geno_out_prefix, ".keep")
  write_tsv(data.frame("#IID" = geno_samples, stringsAsFactors = FALSE, check.names = FALSE),
            keep_file)

  subset_cmd <- c("--pfile", geno_prefix, "--keep", keep_file, "--make-pgen", "--out", geno_out_prefix)
  if (system2("plink2", subset_cmd) != 0) {
    stop(sprintf("plink2 subset failed for %s", geno_out_prefix))
  }

  prune_prefix <- paste0(geno_out_prefix, ".pca_pruned")
  prune_cmd <- c(
    "--pfile", geno_out_prefix,
    "--maf", sprintf("%.4f", geno_pca_maf),
    "--indep-pairwise", geno_pca_indep_window, geno_pca_indep_step,
    sprintf("%.3f", geno_pca_indep_r2),
    "--out", prune_prefix
  )
  if (system2("plink2", prune_cmd) != 0) {
    stop(sprintf("plink2 LD pruning failed for %s", geno_out_prefix))
  }

  prune_in <- paste0(prune_prefix, ".prune.in")
  pca_n <- min(geno_pca_n, length(geno_samples) - 1)
  pca_prefix <- paste0(geno_out_prefix, ".pca")
  pca_cmd <- c("--pfile", geno_out_prefix, "--extract", prune_in, "--pca", pca_n, "--out", pca_prefix)
  if (system2("plink2", pca_cmd) != 0) {
    stop(sprintf("plink2 PCA failed for %s", geno_out_prefix))
  }

  eigenvec <- data.table::fread(paste0(pca_prefix, ".eigenvec"), data.table = FALSE)
  rownames(eigenvec) <- as.character(eigenvec[["#IID"]])
  pc_cols <- grep("^PC[0-9]+$", colnames(eigenvec), value = TRUE)
  geno_pcs <- eigenvec[, pc_cols, drop = FALSE]
  colnames(geno_pcs) <- paste0("genoPC", seq_len(ncol(geno_pcs)))

  known_covs <- cbind(known_covs[rownames(geno_pcs), , drop = FALSE], geno_pcs)
  stage_expr_int <- stage_expr_int[, rownames(geno_pcs), drop = FALSE]
  expr_pca$x <- expr_pca$x[rownames(geno_pcs), , drop = FALSE]

  pcs_top <- as.data.frame(expr_pca$x)[rownames(known_covs), seq_len(k_final), drop = FALSE]
  colnames(pcs_top) <- paste0("exprPC", seq_len(ncol(pcs_top)))

  covs_filtered <- PCAForQTL::filterKnownCovariates(
    known_covs,
    pcs_top,
    unadjustedR2_cutoff = 0.9
  )
  covs_to_use <- cbind(covs_filtered, pcs_top)

  cov_file <- file.path(bc_dir, sprintf("covariates_for_qtl.%s.%s.txt",
                                         selected_method, stage_name))
  write_tsv(as.data.frame(t(covs_to_use)), cov_file, row_names = TRUE)

  final_samples <- rownames(covs_to_use)
  stage_expr_int <- stage_expr_int[, final_samples, drop = FALSE]
  stage_coords <- coords[rownames(stage_expr_int), , drop = FALSE]

  bed_df <- data.frame(
    `#Chr` = stage_coords$chr,
    start = stage_coords$start,
    end = stage_coords$end,
    pid = rownames(stage_expr_int),
    stage_expr_int,
    check.names = FALSE
  )

  chrom_order <- paste0("chr", c(1:22, "X", "Y", "M"))
  bed_df$sort_chr <- factor(bed_df$`#Chr`, levels = chrom_order, ordered = TRUE)
  bed_df <- bed_df[order(bed_df$sort_chr, bed_df$start, bed_df$end), ]
  bed_df$sort_chr <- NULL

  bed_file <- file.path(bc_dir, sprintf("phenotype.%s.%s.bed", selected_method, stage_name))
  write_tsv(bed_df, bed_file)

  if (system2("bgzip", c("-f", bed_file)) != 0) {
    stop(sprintf("bgzip failed for %s", bed_file))
  }
  if (system2("tabix", c("-p", "bed", "-f", paste0(bed_file, ".gz"))) != 0) {
    stop(sprintf("tabix failed for %s.gz", bed_file))
  }

  writeLines(final_samples, file.path(bc_dir, sprintf("sample.%s.%s.final",
                                                       selected_method, stage_name)))
}

log_message(sprintf("Finished. Selected correction method: %s", selected_method))
