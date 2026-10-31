#!/usr/bin/env Rscript
# Preprocess stage-specific APA phenotypes for QTL mapping (module 04).
# Writes per-stage phenotype BED, covariates, sample lists and genotype PCs to data/04_qtl/pre_input/.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

suppressMessages({
library(data.table)
library(dplyr)
library(flashier)
library(PCAForQTL)
library(sva)
library(softImpute)
library(irlba)
library(ebnm)
library(ggplot2)
library(ggrepel)
})

# --- helpers ---

run_or_stop <- function(cmd) {
  code <- system(cmd)
  if (code != 0) stop(sprintf("Command failed with exit code %d: %s", code, cmd))
}

# --- settings ---

set.seed(2026)

apa_root <- Sys.getenv("APA_ROOT")
apa_file <- file.path(apa_root, "raw", "all.apa.count")
meta_file <- file.path(apa_root, "data", "04_qtl", "CHB.all.meta")
geno_prefix <- file.path(apa_root, "data", "00_geno", "genotype")
out_dir <- file.path(apa_root, "data", "04_qtl", "pre_input")
fig_dir <- file.path(apa_root, "figures", "04_qtl")
script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
stage_col <- c(Stage1="#56B4E9", Stage2="#FFB482", Stage3="#8DE5A1")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- phenotype table ---

data <- fread(apa_file, check.names=FALSE)

if (grepl("\\|", data$Gene[1])) {
  data[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed=TRUE)]
} else {
  warning("Gene column does not contain '|'.")
}
data[, phenotype_id := paste(transcript, gene, Loci, strand, sep="|")]
saving_names <- data$phenotype_id

info_cols_indices <- 1:4
info_data <- cbind(data[, ..info_cols_indices], data[, .(phenotype_id)])
fwrite(info_data, file.path(out_dir, "proximal.apa.info"), sep='\t', quote=FALSE)
cols_to_exclude <- c(names(data)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
expr_data <- data[, !cols_to_exclude, with=FALSE]
expr_data <- as.data.frame(expr_data)
rownames(expr_data) <- saving_names

# --- metadata and sample / feature QC ---

meta <- fread(meta_file, data.table=FALSE)
rownames(meta) <- meta[,1]
meta <- meta[,-1]

common_samples <- intersect(colnames(expr_data), rownames(meta))
expr_data <- expr_data[, common_samples]
meta <- meta[common_samples, , drop=FALSE]

sample_na_rate <- colMeans(is.na(expr_data))
keep_samples <- sample_na_rate <= 0.8
expr_data_sample_qc <- expr_data[, keep_samples]
meta_clean <- meta[keep_samples, , drop=FALSE]

row_na_rate <- rowMeans(is.na(expr_data_sample_qc))
row_sd <- apply(expr_data_sample_qc, 1, sd, na.rm = TRUE)
keep_rows <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & (!is.na(row_sd))
expr_data_qc <- expr_data_sample_qc[keep_rows, ]

# --- imputation (softImpute start, flash backfit) ---

message("Imputation running...")
lambda <- if (nrow(expr_data_qc) < 1500) 5 else 30
rank <- 50
X_mis_C <- as(as.matrix(expr_data_qc), "Incomplete")
fit1 <- softImpute::softImpute(X_mis_C, rank = rank, lambda = lambda, type = "svd")
pheno_soft <- softImpute::complete(as.matrix(expr_data_qc), fit1)

min_sd <- min(apply(pheno_soft, 1, sd))
pca_res <- irlba::irlba(pheno_soft, nv = 10)
pca_res_list <- list(d = pca_res$d, u = pca_res$u, v = pca_res$v)

fl_pca <- flash_init(as.matrix(expr_data_qc), S = min_sd, var_type = 1) |>
  flash_factors_init(pca_res_list, ebnm_fn = ebnm::ebnm_point_laplace) |>
  flash_backfit(maxiter = 300)

expr_data_imputed <- ifelse(is.na(as.matrix(expr_data_qc)), fitted(fl_pca), as.matrix(expr_data_qc))

if (any(is.na(expr_data_imputed))) {
  message("Post-imputation NAs found. Filling with row means.")
  row_means <- rowMeans(expr_data_imputed, na.rm = TRUE)
  for (i in 1:nrow(expr_data_imputed)) {
    expr_data_imputed[i, is.na(expr_data_imputed[i, ])] <- row_means[i]
  }
}

# --- outlier detection (Mahalanobis distance in PC space, per stage) ---

n_pcs <- 5
p_value_cutoff <- 0.001
message(sprintf("Running Outlier Detection with Mahalanobis distance (PCs=1:%d)...", n_pcs))

vars <- apply(expr_data_imputed, 1, var)
expr_for_pca <- expr_data_imputed[vars > 1e-10, ]

pca_obj <- prcomp(t(expr_for_pca), center = TRUE, scale. = TRUE)
pcs <- pca_obj$x[, 1:n_pcs]

pca_data_list <- list()
for (stage in unique(meta_clean$LifeStage)) {
  stage_samples <- rownames(meta_clean[meta_clean$LifeStage == stage, ])
  stage_samples <- intersect(stage_samples, rownames(pcs))

  if (length(stage_samples) <= n_pcs) {
    warning(sprintf("Sample size for stage %s is too small for df=%d. Skipping.", stage, n_pcs))
    next
  }

  group_pcs <- pcs[stage_samples, ]
  group_center <- colMeans(group_pcs)
  group_cov <- cov(group_pcs)

  m_dist <- mahalanobis(group_pcs, center = group_center, cov = group_cov)
  p_vals <- 1 - pchisq(m_dist, df = n_pcs)

  df <- data.frame(
    sample = stage_samples,
    PC1 = group_pcs[,1],
    PC2 = group_pcs[,2],
    M_Dist = m_dist,
    P_Value = p_vals,
    LifeStage = stage,
    IsOutlier = p_vals < p_value_cutoff
  )
  pca_data_list[[stage]] <- df
}

pca_data <- do.call(rbind, pca_data_list)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
  geom_point(aes(size = M_Dist, alpha = IsOutlier), shape = 16) +
  scale_color_manual(values = stage_col) +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
  geom_text_repel(data = subset(pca_data, IsOutlier),
                  aes(label = sample), color = "red", size = 3, fontface = "bold") +
  labs(title = "Outlier Detection: Mahalanobis Distance in PC Space",
       subtitle = sprintf("Red labels: p < %s (Chi-square dist, df=%d)", p_value_cutoff, n_pcs),
       x = "PC1", y = "PC2", size = "Mahalanobis Dist") +
  theme_bw() + theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(file.path(fig_dir, "outlier_detect_mahalanobis.png"), p, width = 8, height = 6)

outliers <- pca_data$sample[pca_data$IsOutlier]
if (length(outliers) > 0) {
  message(sprintf("Removing %d outlier samples based on Mahalanobis distance.", length(outliers)))
  samples_keep <- setdiff(colnames(expr_data_imputed), outliers)
  expr_data_imputed <- expr_data_imputed[, samples_keep]
  meta_clean <- meta_clean[samples_keep, ]
}

# --- ComBat batch correction ---

epsilon <- 0.001
dat_bounded <- pmin(pmax(expr_data_imputed, epsilon), 1 - epsilon)
expr_data_scaled <- log2(dat_bounded / (1 - dat_bounded))

if ("batch" %in% colnames(meta_clean)) {
  batch_vec <- as.factor(meta_clean$batch)
  mod <- model.matrix(~ LifeStage, data = meta_clean)
  # reference batch
  ref_batch <- "batch9"
  # This cohort has LifeStage-batch confounding and repeated-sequencing samples.
  # Estimate ComBat on the full post-QC cohort before genotype overlap, downsampling,
  # and stage-specific grouping so batch correction does not depend on the sampled subset.
  message("Running ComBat with ", ref_batch, " as reference batch...")
  expr_clean <- ComBat(dat = as.matrix(expr_data_scaled), batch = batch_vec, mod = mod,
                       ref.batch = if(ref_batch %in% levels(batch_vec)) ref_batch else NULL)
} else {
  expr_clean <- as.matrix(expr_data_scaled)
}

# --- genotype alignment (after ComBat, before rank normalization) ---

message("Performing global genotype alignment...")
psam_file <- paste0(geno_prefix, ".psam")
psam_data <- fread(psam_file, data.table = FALSE)
geno_samples <- psam_data[[1]]  # First column is sample ID

common_samples_global <- intersect(colnames(expr_clean), geno_samples)
message(sprintf("Global alignment: %d expression samples -> %d samples with genotype",
              ncol(expr_clean), length(common_samples_global)))

if (length(common_samples_global) < 10) {
  stop("Insufficient overlapping samples between expression and genotype data!")
}

expr_clean <- expr_clean[, common_samples_global, drop=FALSE]
expr_data_imputed <- expr_data_imputed[, common_samples_global, drop=FALSE]
meta_clean <- meta_clean[common_samples_global, , drop=FALSE]

fwrite(as.data.frame(expr_data_imputed), file.path(out_dir, "after.impute.all.raw"), sep='\t', quote=FALSE, row.names=TRUE)
fwrite(as.data.frame(expr_clean), file.path(out_dir, "after.combat.txt"), sep='\t', quote=FALSE, row.names=TRUE)

stages <- unique(meta_clean$LifeStage)
message(sprintf("After genotype alignment, processing %d stages: %s",
              length(stages), paste(stages, collapse=", ")))

# --- per-stage phenotype, genotype PCs, covariates ---

for (stage in stages) {
  message(sprintf("[%s] processing", stage))

  stage_samples <- rownames(meta_clean)[meta_clean$LifeStage == stage]
  if (length(stage_samples) < 5) {
    message(sprintf("[%s] Skipping: Only %d samples available (< 5 minimum required)", stage, length(stage_samples)))
    next
  }

  sub_expr_clean <- expr_clean[, stage_samples, drop=FALSE]
  sub_meta <- meta_clean[stage_samples, , drop=FALSE]

  sub_expr_norm <- apply(as.matrix(sub_expr_clean), 1, function(x) {
    n_val <- sum(!is.na(x))
    qnorm(rank(x, na.last = "keep", ties.method = "average") / (n_val+1))
  }) %>% t()
  sub_expr_norm <- sub_expr_norm[complete.cases(sub_expr_norm), ]
  out_filename <- sprintf("pre.phenotype.combat.int.%s.txt", stage)
  fwrite(as.data.frame(sub_expr_norm), file.path(out_dir, out_filename), sep='\t', quote=FALSE, row.names=TRUE)

  # stage sample list, then genotype extraction, LD pruning and PCA (plink2)
  stage_sample_file <- file.path(out_dir, sprintf("samples.%s.txt", tolower(stage)))
  writeLines(stage_samples, stage_sample_file)

  pgen_prefix <- file.path(out_dir, paste0("genotype.", tolower(stage)))
  run_or_stop(sprintf("plink2 --pfile %s --keep %s --make-pgen --out %s",
                      geno_prefix, stage_sample_file, pgen_prefix))

  # LD-prune before PCA so PCs are not dominated by correlated variants in long haplotypes or high-LD regions.
  prune_prefix <- file.path(out_dir, paste0("genotype_prune.", tolower(stage)))
  run_or_stop(sprintf("plink2 --pfile %s --indep-pairwise 50 5 0.2 --out %s",
                      pgen_prefix, prune_prefix))
  prune_in_file <- paste0(prune_prefix, ".prune.in")

  stage_pca_prefix <- file.path(out_dir, paste0("genotype_pca.", tolower(stage)))
  run_or_stop(sprintf("plink2 --pfile %s --extract %s --pca 20 --output-chr chrM --out %s",
                      pgen_prefix, prune_in_file, stage_pca_prefix))
  stage_eigenvec <- paste0(stage_pca_prefix, ".eigenvec")

  # covariates: PCAForQTL auto-K on the stage phenotype PCs, filtered known covariates, genotype PCs
  vars_norm <- apply(sub_expr_norm, 1, var)
  sub_expr_for_pca <- sub_expr_norm[vars_norm > 1e-10, ]

  if (ncol(sub_expr_for_pca) >= 2 && nrow(sub_expr_for_pca) >= 2) {
    sub_pca_final <- prcomp(t(sub_expr_for_pca), center = TRUE, scale. = TRUE)
    K_BE <- as.integer(PCAForQTL::runBE(t(sub_expr_for_pca), B = 20)$numOfPCsChosen)
    K_final <- min(max(K_BE, 10), 50)

    remove_col <- c("batch", "sample_id", "LifeStage")
    known_covs <- sub_meta %>% mutate(across(where(is.character), as.factor)) %>%
      mutate(across(where(is.factor), as.numeric)) %>% select(where(is.numeric)) %>% select(-any_of(remove_col))

    # stage-specific genotype PCs
    geno_df <- fread(stage_eigenvec, data.table = FALSE)
    if (!"#IID" %in% colnames(geno_df)) {
      stop(sprintf("[%s] Unexpected eigenvec format: missing #IID column in %s", stage, stage_eigenvec))
    }
    rownames(geno_df) <- geno_df[["#IID"]]
    geno_sub <- geno_df[, 2:min(6, ncol(geno_df)), drop = FALSE]
    colnames(geno_sub) <- paste0("genoPC", seq_len(ncol(geno_sub)))
    common_cov_samples <- intersect(rownames(geno_sub), rownames(known_covs))
    if (length(common_cov_samples) > 0) {
      known_covs <- cbind(known_covs[common_cov_samples,], geno_sub[common_cov_samples,])
    }

    PCsTop <- as.data.frame(sub_pca_final$x)[rownames(known_covs), 1:K_final, drop = FALSE]
    covsFiltered <- PCAForQTL::filterKnownCovariates(known_covs, PCsTop, unadjustedR2_cutoff = 0.9)
    covsToUse <- cbind(covsFiltered, PCsTop)

    out_filename <- sprintf("covariates_for_qtl.%s.txt", stage)
    fwrite(as.data.frame(t(covsToUse)), file.path(out_dir, out_filename), sep = "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

    # phenotype BED (gene, loci and strand parsed from the phenotype id)
    sub_expr_final <- as.data.frame(sub_expr_norm)
    sub_expr_final$pid <- rownames(sub_expr_norm)
    setDT(sub_expr_final)
    sub_expr_final[, c("transcript", "gene", "Loci", "strand") := tstrsplit(pid, "|", fixed=TRUE)]
    sub_expr_final[, c("chr", "region") := tstrsplit(Loci, ":", fixed=TRUE)]
    sub_expr_final[, c("start", "end") := tstrsplit(region, "-", fixed=TRUE)]
    sub_expr_final[, `:=`(start = as.numeric(start), end = as.numeric(end), `#Chr` = chr)]
    final_samples <- rownames(covsToUse)
    format_cols <- c("#Chr", "start", "end", "pid", final_samples)

    format_data <- sub_expr_final[, ..format_cols]
    row_count_before <- nrow(format_data)
    format_data <- unique(format_data, by = "pid")
    row_count_after <- nrow(format_data)

    if (row_count_before > row_count_after) {
      message(sprintf("[%s] Warning: Removed %d duplicated phenotype rows.",
                      stage, row_count_before - row_count_after))
    }
    setorder(format_data, `#Chr`, start, end)

    bed_file <- file.path(out_dir, sprintf("phenotype.%s.bed", stage))
    fwrite(format_data, bed_file, sep='\t', quote=FALSE)

    # CHR prefix must match the genotype files
    geno_chr_sample <- fread(paste0(pgen_prefix, ".pvar"), nrows = 1, header = FALSE)$V1
    has_chr_prefix <- grepl("^chr", geno_chr_sample)

    format_data_check <- fread(bed_file)
    pheno_chr_sample <- format_data_check$`#Chr`[1]
    pheno_has_prefix <- grepl("^chr", pheno_chr_sample)

    if (has_chr_prefix != pheno_has_prefix) {
      if (has_chr_prefix) {
        message(sprintf("[%s] Adding 'chr' prefix to phenotype CHR to match genotype...", stage))
        format_data_check[, `#Chr` := paste0("chr", `#Chr`)]
      } else {
        message(sprintf("[%s] Removing 'chr' prefix from phenotype CHR to match genotype...", stage))
        format_data_check[, `#Chr` := sub("^chr", "", `#Chr`)]
      }
      fwrite(format_data_check, bed_file, sep='\t', quote=FALSE)
    }

    system(paste("bgzip -f", bed_file))
    system(paste("tabix -p bed -f", paste0(bed_file, ".gz")))

    sample_file <- file.path(out_dir, sprintf("sample.%s.final", stage))
    writeLines(final_samples, sample_file)
  } else {
    warning(sprintf("[%s] Insufficient samples or phenotypes for covariate analysis. Skipping.", stage))
  }
}

# --- group files and GeneRange.bed ---

format_cmd <- sprintf("python3 %s %s %s", file.path(script_dir, "format.py"), out_dir, out_dir)
format_exit <- system(format_cmd)
if (format_exit != 0) {
  warning("format.py failed with exit code ", format_exit)
} else {
  message("format.py completed successfully")
}
