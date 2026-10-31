# Preprocess APA counts for G x Stage: imputation, outlier removal, ComBat, INT,
# covariates and phenotype BED. Writes data/05_gxs/ and figures/05_gxs/.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(sva)
  library(softImpute)
  library(irlba)
  library(flashier)
  library(ebnm)
  library(PCAForQTL)
  library(ggplot2)
  library(ggrepel)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
set.seed(2026)
raw_dir <- file.path(APA_ROOT, "raw")
out_dir <- file.path(APA_ROOT, "data", "05_gxs")
fig_dir <- file.path(APA_ROOT, "figures", "05_gxs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
apa_file <- file.path(raw_dir, "all.apa.count")
meta_file <- file.path(APA_ROOT, "data", "04_qtl", "CHB.all.meta")
geno_prefix <- file.path(raw_dir, "genotype")
ref_batch <- "batch9"
stage_col <- c(Prenatal = "#56B4E9", Postnatal = "#FFB482")

# --- read APA counts and metadata ---
data <- fread(apa_file, check.names = FALSE)
data[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed = TRUE)]
data[, phenotype_id := paste(transcript, gene, Loci, strand, sep = "|")]
saving_names <- data$phenotype_id
info_data <- cbind(data[, 1:4, with = FALSE], data[, .(phenotype_id)])
fwrite(info_data, file.path(out_dir, "proximal.apa.info"), sep = "\t", quote = FALSE)

cols_to_exclude <- c(names(data)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
expr_data <- as.data.frame(data[, !cols_to_exclude, with = FALSE])
rownames(expr_data) <- saving_names

meta <- fread(meta_file, data.table = FALSE)
rownames(meta) <- meta[, 1]
meta <- meta[, -1]
common_samples <- intersect(colnames(expr_data), rownames(meta))
expr_data <- expr_data[, common_samples]
meta <- meta[common_samples, , drop = FALSE]

# --- sample and feature QC ---
sample_na_rate <- colMeans(is.na(expr_data))
keep_samples <- sample_na_rate <= 0.8
expr_data_sample_qc <- expr_data[, keep_samples]
meta_clean <- meta[keep_samples, , drop = FALSE]

row_na_rate <- rowMeans(is.na(expr_data_sample_qc))
row_sd <- apply(expr_data_sample_qc, 1, sd, na.rm = TRUE)
keep_rows <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & (!is.na(row_sd))
expr_data_qc <- expr_data_sample_qc[keep_rows, ]

# --- imputation (softImpute, then flash backfit) ---
lambda <- if (nrow(expr_data_qc) < 1500) 5 else 30
rank <- 50
maxiter_flash <- 300
X_mis_C <- as(as.matrix(expr_data_qc), "Incomplete")
fit1 <- softImpute::softImpute(X_mis_C, rank = rank, lambda = lambda, type = "svd")
pheno_soft <- softImpute::complete(as.matrix(expr_data_qc), fit1)

min_sd <- min(apply(pheno_soft, 1, sd))
pca_res <- irlba::irlba(pheno_soft, nv = 10)
pca_res_list <- list(d = pca_res$d, u = pca_res$u, v = pca_res$v)
fl_pca <- flash_init(as.matrix(expr_data_qc), S = min_sd, var_type = 1) |>
  flash_factors_init(pca_res_list, ebnm_fn = ebnm::ebnm_point_laplace) |>
  flash_backfit(maxiter = maxiter_flash)
imputed_vals <- fitted(fl_pca)
original_vals <- as.matrix(expr_data_qc)
expr_data_imputed <- original_vals
na_indices <- is.na(original_vals)
expr_data_imputed[na_indices] <- imputed_vals[na_indices]

if (any(is.na(expr_data_imputed))) {
  message("Warning: Post-imputation NAs found. Filling with row means.")
  row_means <- rowMeans(expr_data_imputed, na.rm = TRUE)
  for (i in 1:nrow(expr_data_imputed)) {
    expr_data_imputed[i, is.na(expr_data_imputed[i, ])] <- row_means[i]
  }
}

# --- outlier detection (Mahalanobis distance in PC space, per life stage) ---
n_pcs <- 5
p_value_cutoff <- 0.001
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
  pca_data_list[[stage]] <- data.frame(
    sample = stage_samples,
    PC1 = group_pcs[, 1],
    PC2 = group_pcs[, 2],
    M_Dist = m_dist,
    P_Value = p_vals,
    LifeStage = stage,
    IsOutlier = p_vals < p_value_cutoff
  )
}
pca_data <- do.call(rbind, pca_data_list)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = LifeStage)) +
  geom_point(aes(size = M_Dist, alpha = IsOutlier), shape = 16) +
  scale_color_manual(values = stage_col) +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
  geom_text_repel(data = subset(pca_data, IsOutlier), aes(label = sample),
                  color = "red", size = 3, fontface = "bold") +
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

# --- ComBat (batch, with intercept-only fallback when confounded) ---
fwrite(as.data.frame(expr_data_imputed), file.path(out_dir, "before.combat.txt"), sep = "\t", quote = FALSE, row.names = TRUE)
epsilon <- 0.001
expr_data_imputed <- as.matrix(expr_data_imputed)
dat_bounded <- pmin(pmax(expr_data_imputed, epsilon), 1 - epsilon)
expr_data_scaled <- log2(dat_bounded / (1 - dat_bounded))

if ("batch" %in% colnames(meta_clean)) {
  batch_vec <- as.factor(meta_clean$batch)
  life_stage_vec <- as.factor(meta_clean$LifeStage)
  if (length(unique(batch_vec)) < 2 || length(unique(life_stage_vec)) < 2) {
    message("Only one batch or LifeStage level after filtering. Skipping ComBat.")
    expr_clean <- as.matrix(expr_data_scaled)
  } else {
    mod <- model.matrix(~ LifeStage, data = meta_clean)
    design_batch <- model.matrix(~ batch_vec)
    full_design <- cbind(mod[, -1, drop = FALSE], design_batch[, -1, drop = FALSE])
    if (qr(full_design)$rank < ncol(full_design)) {
      message("ComBat: covariate(s) confounded with batch; falling back to intercept-only model.")
      mod <- model.matrix(~ 1, data = meta_clean)
    }
    expr_clean <- ComBat(dat = as.matrix(expr_data_scaled), batch = batch_vec, mod = mod,
                         ref.batch = if (ref_batch %in% levels(batch_vec)) ref_batch else NULL)
  }
} else {
  expr_clean <- as.matrix(expr_data_scaled)
}

# --- align to genotype samples and compute genotype PCs ---
psam_data <- fread(paste0(geno_prefix, ".psam"), data.table = FALSE)
geno_samples <- psam_data[[1]]
survivors <- colnames(expr_clean)
targets <- intersect(geno_samples, survivors)
if (length(targets) == 0) {
  stop("No overlapping samples between genotype and expression/meta after outlier removal.")
}
dropped_from_geno <- setdiff(geno_samples, targets)
dropped_from_expr <- setdiff(survivors, targets)
if (length(dropped_from_geno) > 0) {
  message(sprintf("Dropping %d genotype samples not in expression/meta: %s", length(dropped_from_geno),
                  paste(head(dropped_from_geno, 10), collapse = ", ")))
}
if (length(dropped_from_expr) > 0) {
  message(sprintf("Dropping %d expression/meta samples not in genotype: %s", length(dropped_from_expr),
                  paste(head(dropped_from_expr, 10), collapse = ", ")))
}
expr_clean <- expr_clean[, targets, drop = FALSE]
meta_clean <- meta_clean[targets, , drop = FALSE]
fwrite(as.data.frame(expr_clean), file.path(out_dir, "after.combat.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# LD pruning and PCA on aligned samples
pca_sample_file <- file.path(out_dir, "pca_samples.txt")
writeLines(targets, pca_sample_file)
prune_prefix <- file.path(out_dir, "genotype_pruned")
stopifnot(system(sprintf("plink2 --pfile %s --keep %s --indep-pairwise 50 5 0.2 --out %s",
                         geno_prefix, pca_sample_file, prune_prefix)) == 0)
stopifnot(system(sprintf("plink2 --pfile %s --keep %s --extract %s --pca 20 --output-chr chrM --out %s",
                         geno_prefix, pca_sample_file, paste0(prune_prefix, ".prune.in"),
                         file.path(out_dir, "genotype_pca"))) == 0)

# --- inverse-normal transform ---
apply_int <- function(expr_matrix) {
  apply(as.matrix(expr_matrix), 1, function(x) {
    n_val <- sum(!is.na(x))
    qnorm(rank(x, na.last = "keep", ties.method = "average") / (n_val + 1))
  }) %>% t()
}
expr_data_norm <- apply_int(expr_clean)
expr_data_norm <- expr_data_norm[complete.cases(expr_data_norm), ]
expr_data_norm <- expr_data_norm[, targets, drop = FALSE]
fwrite(as.data.frame(expr_data_norm), file.path(out_dir, "pre.phenotype.combat.int.txt"),
       sep = "\t", quote = FALSE, row.names = TRUE)

# --- expression PCs (runBE chooses K) ---
vars_norm <- apply(expr_data_norm, 1, var)
expr_for_final_pca <- expr_data_norm[vars_norm > 1e-10, ]
pca_final <- prcomp(t(expr_for_final_pca), center = TRUE, scale. = TRUE)
K_BE <- as.integer(PCAForQTL::runBE(t(expr_for_final_pca), B = 20)$numOfPCsChosen)
K_final <- min(max(K_BE, 10), 50)

# --- covariates: clinical covariates + genotype PCs + expression PCs ---
remove_col <- c("batch", "sample_id", "LifeStage")
known_covs <- meta_clean %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.factor), as.numeric)) %>%
  select(where(is.numeric)) %>%
  select(-any_of(remove_col))
geno_df <- fread(file.path(out_dir, "genotype_pca.eigenvec"), data.table = FALSE)
rownames(geno_df) <- geno_df[["#IID"]]
geno_sub <- geno_df[, 2:min(6, ncol(geno_df)), drop = FALSE]
colnames(geno_sub) <- paste0("genoPC", seq_len(ncol(geno_sub)))
geno_sub <- geno_sub[targets, , drop = FALSE]
stopifnot(identical(rownames(known_covs), targets))
known_covs <- cbind(known_covs, geno_sub)

PCsTop <- as.data.frame(pca_final$x)[rownames(known_covs), 1:K_final, drop = FALSE]
PCsTop <- PCsTop[targets, , drop = FALSE]
covsToUse <- cbind(known_covs, PCsTop)
stopifnot(identical(rownames(covsToUse), targets))
fwrite(as.data.frame(t(covsToUse)), file.path(out_dir, "covariates_for_qtl.txt"),
       sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# --- phenotype BED (sorted, CHR prefix matched to genotype) ---
expr_final <- as.data.frame(expr_data_norm)
expr_final$pid <- rownames(expr_data_norm)
setDT(expr_final)
expr_final[, c("transcript", "gene", "Loci", "strand") := tstrsplit(pid, "|", fixed = TRUE)]
expr_final[, c("chr", "region") := tstrsplit(Loci, ":", fixed = TRUE)]
expr_final[, c("start", "end") := tstrsplit(region, "-", fixed = TRUE)]
expr_final[, `:=`(start = as.numeric(start), end = as.numeric(end), `#Chr` = chr)]
stopifnot(identical(rownames(pca_final$x), targets))
format_data <- expr_final[, c("#Chr", "start", "end", "pid", targets), with = FALSE]
setorder(format_data, `#Chr`, start, end)

bed_file <- file.path(out_dir, "phenotype.bed")
fwrite(format_data, bed_file, sep = "\t", quote = FALSE)

geno_chr_sample <- fread(paste0(geno_prefix, ".pvar"), nrows = 1, header = FALSE)$V1
has_chr_prefix <- grepl("^chr", geno_chr_sample)
format_data_check <- fread(bed_file)
pheno_has_prefix <- grepl("^chr", format_data_check$`#Chr`[1])
if (has_chr_prefix != pheno_has_prefix) {
  if (has_chr_prefix) {
    format_data_check[, `#Chr` := paste0("chr", `#Chr`)]
  } else {
    format_data_check[, `#Chr` := sub("^chr", "", `#Chr`)]
  }
  fwrite(format_data_check, bed_file, sep = "\t", quote = FALSE)
}
system(paste("bgzip -f", bed_file))
system(paste("tabix -p bed -f", paste0(bed_file, ".gz")))

# --- gene-range phenotype groups (phenotype.group, GeneRange BED) ---
gene_range <- copy(format_data_check)
bed_cols <- names(gene_range)
chr_col <- bed_cols[1]
start_col <- bed_cols[2]
end_col <- bed_cols[3]
pid_col <- bed_cols[4]

gene_range[, gid_raw := tstrsplit(get(pid_col), "|", fixed = TRUE, keep = 2)]
gene_range[is.na(gid_raw), gid_raw := get(pid_col)]
locus_coords <- gene_range[, .(new_start = min(get(start_col)), new_end = max(get(end_col))),
                           by = c(chr_col, "gid_raw")]
gene_range <- merge(gene_range, locus_coords, by = c(chr_col, "gid_raw"), all.x = TRUE)
set(gene_range, j = start_col, value = gene_range$new_start)
set(gene_range, j = end_col, value = gene_range$new_end)
gene_range[, locus_id := paste(get(chr_col), as.integer(get(start_col)), sep = ":")]
locus_counts <- gene_range[, .(n_loci = uniqueN(locus_id)), by = gid_raw]
gene_range <- merge(gene_range, locus_counts, by = "gid_raw", all.x = TRUE)
gene_range[, gid := gid_raw]
gene_range[n_loci > 1, gid := paste(gid_raw, locus_id, sep = "|")]
setorderv(gene_range, c(chr_col, start_col, end_col, "gid", pid_col))

phenotype_group <- gene_range[, .(phenotype_id = get(pid_col), gid = gid)]
fwrite(phenotype_group, file.path(out_dir, "phenotype.group"), sep = "\t", col.names = FALSE)

gene_range_file <- file.path(out_dir, "phenotype.GeneRange.bed")
fwrite(gene_range[, ..bed_cols], gene_range_file, sep = "\t", quote = FALSE)
system(paste("bgzip -f", gene_range_file))
system(paste("tabix -p bed -f", paste0(gene_range_file, ".gz")))

writeLines(targets, file.path(out_dir, "sample.final"))
