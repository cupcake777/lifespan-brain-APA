#!/usr/bin/env Rscript
for (v in c("APA_ROOT", "APA_HPC_HOME")) if (!nzchar(Sys.getenv(v))) stop("Set ", v, " (source config/paths.sh)")
# Stage-specific APA-QTL downsampling: expression QC, imputation, ComBat, genotype alignment.
# Per seed: sex-stratified postnatal sampling, LD-pruned genotype PCA, covariates, phenotype BED.

suppressPackageStartupMessages({
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

normalize_sex <- function(x) {
  y <- trimws(as.character(x))
  y_low <- tolower(y)
  fifelse(y_low %in% c("f", "female", "woman", "2"), "Female",
    fifelse(y_low %in% c("m", "male", "man", "1"), "Male", y)
  )
}

largest_remainder_counts <- function(stratum_counts, target_n) {
  raw <- target_n * stratum_counts / sum(stratum_counts)
  base <- floor(raw)
  remainder <- target_n - sum(base)
  if (remainder > 0) {
    ord <- order(raw - base, stratum_counts, decreasing = TRUE)
    base[ord[seq_len(remainder)]] <- base[ord[seq_len(remainder)]] + 1L
  }
  base
}

# --- inputs ---

set.seed(2026)

qtl_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
apa_file <- file.path(Sys.getenv("APA_ROOT"), "raw", "all.apa.count")
meta_file <- file.path(qtl_dir, "CHB.all.meta")
geno_prefix <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/00_raw_data/genotype")
downsample_base_dir <- file.path(qtl_dir, "downsample")
shared_out_dir <- file.path(downsample_base_dir, "shared")
audit_path <- file.path(qtl_dir, "downsample_audit_seed_composition.tsv")
seed_env <- Sys.getenv("DOWNSAMPLE_SEEDS", "")
downsample_seeds <- if (nzchar(seed_env)) {
  as.integer(strsplit(seed_env, "[,[:space:]]+")[[1]])
} else {
  2026:2045
}
downsample_seeds <- downsample_seeds[!is.na(downsample_seeds)]
if (length(downsample_seeds) == 0) stop("No valid downsample seeds were provided.")
stage_col <- c(Stage1="#56B4E9", Stage2="#FFB482", Stage3="#8DE5A1")

dir.create(shared_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(Sys.getenv("APA_ROOT"), "figures", "04_qtl"), recursive = TRUE, showWarnings = FALSE)

data <- fread(apa_file, check.names=FALSE)

if (grepl("\\|", data$Gene[1])) {
  data[, c("transcript", "gene", "chr", "strand") := tstrsplit(Gene, "|", fixed=TRUE)]
}
data[, phenotype_id := paste(transcript, gene, Loci, strand, sep="|")]
saving_names <- data$phenotype_id

# --- proximal APA info ---

info_cols_indices <- 1:4
info_data <- cbind(data[, ..info_cols_indices], data[, .(phenotype_id)])
fwrite(info_data, file.path(shared_out_dir, "proximal.apa.info"), sep='\t', quote=FALSE)
cols_to_exclude <- c(names(data)[1:4], "transcript", "gene", "chr", "strand", "phenotype_id")
expr_data <- data[, !cols_to_exclude, with=FALSE]
expr_data <- as.data.frame(expr_data)
rownames(expr_data) <- saving_names

meta <- fread(meta_file, data.table=FALSE)
rownames(meta) <- meta[,1]
meta <- meta[,-1]

common_samples <- intersect(colnames(expr_data), rownames(meta))
expr_data <- expr_data[, common_samples]
meta <- meta[common_samples, , drop=FALSE]

# --- QC ---

sample_na_rate <- colMeans(is.na(expr_data))
keep_samples <- sample_na_rate <= 0.8
expr_data_sample_qc <- expr_data[, keep_samples]
meta_clean <- meta[keep_samples, , drop=FALSE]

row_na_rate <- rowMeans(is.na(expr_data_sample_qc))
row_sd <- apply(expr_data_sample_qc, 1, sd, na.rm = TRUE)
keep_rows <- (row_na_rate <= 0.5) & (row_sd >= 0.001) & (!is.na(row_sd))
expr_data_qc <- expr_data_sample_qc[keep_rows, ]

# --- imputation ---

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
  row_means <- rowMeans(expr_data_imputed, na.rm = TRUE)
  for (i in 1:nrow(expr_data_imputed)) {
    expr_data_imputed[i, is.na(expr_data_imputed[i, ])] <- row_means[i]
  }
}

# --- outlier samples (Mahalanobis distance in top PCs) ---

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

  if (length(stage_samples) <= n_pcs) next

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
ggsave(file.path(Sys.getenv("APA_ROOT"), "figures", "04_qtl", "outlier_detect_mahalanobis.png"), p, width = 8, height = 6)

outliers <- pca_data$sample[pca_data$IsOutlier]
if (length(outliers) > 0) {
  samples_keep <- setdiff(colnames(expr_data_imputed), outliers)
  expr_data_imputed <- expr_data_imputed[, samples_keep]
  meta_clean <- meta_clean[samples_keep, ]
}

# --- scaling and ComBat ---

epsilon <- 0.001
dat_bounded <- pmin(pmax(expr_data_imputed, epsilon), 1 - epsilon)
expr_data_scaled <- log2(dat_bounded / (1 - dat_bounded))

if ("batch" %in% colnames(meta_clean)) {
  batch_vec <- as.factor(meta_clean$batch)
  mod <- model.matrix(~ LifeStage, data = meta_clean)
  ref_batch <- "batch9"
  # Estimate ComBat on the full post-QC cohort, before genotype overlap and downsampling,
  # so batch correction does not depend on the sampled subset.
  expr_clean <- ComBat(dat = as.matrix(expr_data_scaled), batch = batch_vec, mod = mod,
                       ref.batch = if(ref_batch %in% levels(batch_vec)) ref_batch else NULL)
} else {
  expr_clean <- as.matrix(expr_data_scaled)
}

# --- genotype alignment ---

psam_file <- paste0(geno_prefix, ".psam")
psam_data <- fread(psam_file, data.table = FALSE)
geno_samples <- psam_data[[1]]

common_samples_global <- intersect(colnames(expr_clean), geno_samples)
if (length(common_samples_global) < 10) {
  stop("Insufficient overlapping samples between expression and genotype data!")
}

expr_clean <- expr_clean[, common_samples_global, drop=FALSE]
expr_data_imputed <- expr_data_imputed[, common_samples_global, drop=FALSE]
meta_clean <- meta_clean[common_samples_global, , drop=FALSE]

if (any(is.na(meta_clean$LifeStage) | meta_clean$LifeStage == "")) {
  missing_stage_samples <- rownames(meta_clean)[is.na(meta_clean$LifeStage) | meta_clean$LifeStage == ""]
  stop(
    "Downsample requires non-missing LifeStage for all samples. Missing LifeStage samples: ",
    paste(missing_stage_samples, collapse = ", ")
  )
}
stage_sample_list <- split(rownames(meta_clean), meta_clean$LifeStage)
stage_counts <- lengths(stage_sample_list)
if (length(stage_counts) < 2) {
  stop("Downsample requires at least two LifeStage groups after genotype alignment.")
}
downsample_n <- min(stage_counts)

postnatal_stage_name <- names(stage_counts)[tolower(names(stage_counts)) == "postnatal"]
if (length(postnatal_stage_name) != 1) {
  stop("Expected exactly one Postnatal LifeStage after genotype alignment.")
}
full_postnatal_meta <- meta_clean[meta_clean$LifeStage == postnatal_stage_name, , drop = FALSE]
full_postnatal_age <- suppressWarnings(as.numeric(full_postnatal_meta$age))
if (all(!is.finite(full_postnatal_age))) {
  stop("Postnatal age audit requires numeric age values.")
}
full_postnatal_sex <- normalize_sex(full_postnatal_meta$sex)
full_postnatal_sex_counts <- table(full_postnatal_sex)
if (!all(c("Female", "Male") %in% names(full_postnatal_sex_counts))) {
  stop("Full postnatal samples must include both Female and Male strata.")
}

# --- per-seed downsampling and QTL inputs ---

audit_v2_rows <- list()

for (downsample_seed in downsample_seeds) {
  out_dir <- file.path(downsample_base_dir, paste0("seed_", downsample_seed), "pre_input")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(info_data, file.path(out_dir, "proximal.apa.info"), sep='\t', quote=FALSE)

  set.seed(downsample_seed)
  downsampled_samples <- unlist(lapply(stage_sample_list, function(samples) {
    stage_name <- unique(meta_clean[samples, "LifeStage", drop = TRUE])
    if (length(stage_name) != 1) {
      stop("Each downsample group must contain exactly one LifeStage.")
    }
    if (tolower(stage_name) == "postnatal") {
      # sex-stratified draw of postnatal samples
      post_meta <- copy(as.data.table(meta_clean[samples, , drop = FALSE], keep.rownames = "sample"))
      post_meta[, sex_stratum := normalize_sex(sex)]
      if (any(is.na(post_meta$sex_stratum) | post_meta$sex_stratum == "")) {
        stop("Postnatal sex-stratified downsample found missing sex labels.")
      }
      sex_counts <- post_meta[, .N, by = sex_stratum]
      if (!all(c("Female", "Male") %in% sex_counts$sex_stratum)) {
        stop("Postnatal sex-stratified downsample requires both Female and Male strata.")
      }
      sex_counts <- sex_counts[match(c("Female", "Male"), sex_stratum)]
      draw_counts <- largest_remainder_counts(sex_counts$N, downsample_n)
      names(draw_counts) <- sex_counts$sex_stratum

      unlist(lapply(names(draw_counts), function(one_sex) {
        pool <- post_meta[sex_stratum == one_sex, sample]
        sample(pool, draw_counts[[one_sex]])
      }), use.names = FALSE)
    } else if (length(samples) > downsample_n) {
      sample(samples, downsample_n)
    } else {
      samples
    }
  }), use.names = FALSE)
  downsampled_samples <- intersect(common_samples_global, downsampled_samples)

  expr_clean_seed <- expr_clean[, downsampled_samples, drop=FALSE]
  meta_seed <- meta_clean[downsampled_samples, , drop=FALSE]
  downsampled_counts <- table(meta_seed$LifeStage)
  post_seed_meta <- meta_seed[meta_seed$LifeStage == postnatal_stage_name, , drop = FALSE]
  post_seed_sex <- normalize_sex(post_seed_meta$sex)
  post_seed_age <- suppressWarnings(as.numeric(post_seed_meta$age))
  age_ks_p <- suppressWarnings(tryCatch(
    ks.test(post_seed_age[is.finite(post_seed_age)], full_postnatal_age[is.finite(full_postnatal_age)])$p.value,
    error = function(e) NA_real_
  ))
  audit_v2_seed <- data.table(
    seed = downsample_seed,
    n_female = unname(sum(post_seed_sex == "Female", na.rm = TRUE)),
    n_male = unname(sum(post_seed_sex == "Male", na.rm = TRUE)),
    age_mean = mean(post_seed_age, na.rm = TRUE),
    age_median = median(post_seed_age, na.rm = TRUE),
    age_sd = sd(post_seed_age, na.rm = TRUE),
    age_ks_p_vs_full = age_ks_p,
    K_PC_selected = NA_integer_
  )

  downsample_manifest <- data.frame(sample = rownames(meta_seed), meta_seed, check.names = FALSE)
  fwrite(
    downsample_manifest,
    file.path(out_dir, "downsampled.samples.txt"),
    sep='\t', quote=FALSE
  )
  downsample_summary <- data.frame(
    LifeStage = names(stage_counts),
    before = as.integer(stage_counts),
    after = as.integer(downsampled_counts[names(stage_counts)]),
    target = downsample_n,
    seed = downsample_seed
  )
  fwrite(
    downsample_summary,
    file.path(out_dir, "downsample_summary.txt"),
    sep='\t', quote=FALSE
  )

  stages <- unique(meta_seed$LifeStage)

  for (stage in stages) {
    # skip stages with fewer than 5 samples
    stage_samples <- rownames(meta_seed)[meta_seed$LifeStage == stage]
    if (length(stage_samples) < 5) next

    sub_expr_clean <- expr_clean_seed[, stage_samples, drop=FALSE]
    sub_meta <- meta_seed[stage_samples, , drop=FALSE]

    sub_expr_norm <- apply(as.matrix(sub_expr_clean), 1, function(x) {
      n_val <- sum(!is.na(x))
      qnorm(rank(x, na.last = "keep", ties.method = "average") / (n_val+1))
    }) %>% t()
    sub_expr_norm <- sub_expr_norm[complete.cases(sub_expr_norm), ]

    # genotype of this stage, LD-pruned before PCA
    stage_sample_file <- file.path(out_dir, sprintf("samples.%s.txt", tolower(stage)))
    writeLines(stage_samples, stage_sample_file)

    pgen_prefix <- file.path(out_dir, paste0("genotype.", tolower(stage)))
    plink_cmd <- sprintf("plink2 --pfile %s --keep %s --make-pgen --out %s",
                         geno_prefix, stage_sample_file, pgen_prefix)
    if (system(plink_cmd) != 0) stop(sprintf("[%s] plink2 genotype extraction failed", stage))

    prune_prefix <- file.path(out_dir, paste0("genotype_prune.", tolower(stage)))
    prune_cmd <- sprintf("plink2 --pfile %s --indep-pairwise 50 5 0.2 --out %s",
                         pgen_prefix, prune_prefix)
    if (system(prune_cmd) != 0) stop(sprintf("[%s] plink2 LD pruning failed", stage))
    prune_in_file <- paste0(prune_prefix, ".prune.in")

    stage_pca_prefix <- file.path(out_dir, paste0("genotype_pca.", tolower(stage)))
    pca_cmd <- sprintf("plink2 --pfile %s --extract %s --pca 20 --output-chr chrM --out %s",
                       pgen_prefix, prune_in_file, stage_pca_prefix)
    if (system(pca_cmd) != 0) stop(sprintf("[%s] plink2 PCA failed", stage))
    stage_eigenvec <- paste0(stage_pca_prefix, ".eigenvec")

    # covariates and phenotype BED
    vars_norm <- apply(sub_expr_norm, 1, var)
    sub_expr_for_pca <- sub_expr_norm[vars_norm > 1e-10, ]

    if (ncol(sub_expr_for_pca) >= 2 && nrow(sub_expr_for_pca) >= 2) {
      sub_pca_final <- prcomp(t(sub_expr_for_pca), center = TRUE, scale. = TRUE)
      K_BE <- as.integer(PCAForQTL::runBE(t(sub_expr_for_pca), B = 20)$numOfPCsChosen)
      K_final <- min(max(K_BE, 10), 50)
      if (tolower(stage) == "postnatal") {
        audit_v2_seed[, K_PC_selected := K_final]
      }

      remove_col <- c("batch", "sample_id", "LifeStage")
      known_covs <- sub_meta %>% mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), as.numeric)) %>% select(where(is.numeric)) %>% select(-any_of(remove_col))

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
      format_data <- unique(format_data, by = "pid")
      setorder(format_data, `#Chr`, start, end)

      bed_file <- file.path(out_dir, sprintf("phenotype.%s.bed", stage))
      fwrite(format_data, bed_file, sep='\t', quote=FALSE)

      # match chr prefix of the phenotype BED to the genotype
      geno_chr_sample <- fread(paste0(pgen_prefix, ".pvar"), nrows = 1, header = FALSE)$V1
      has_chr_prefix <- grepl("^chr", geno_chr_sample)

      format_data_check <- fread(bed_file)
      pheno_chr_sample <- format_data_check$`#Chr`[1]
      pheno_has_prefix <- grepl("^chr", pheno_chr_sample)

      if (has_chr_prefix != pheno_has_prefix) {
        if (has_chr_prefix) {
          format_data_check[, `#Chr` := paste0("chr", `#Chr`)]
        } else {
          format_data_check[, `#Chr` := sub("^chr", "", `#Chr`)]
        }
        fwrite(format_data_check, bed_file, sep='\t', quote=FALSE)
      }

      system(paste("bgzip -f", bed_file))
      system(paste("tabix -p bed -f", paste0(bed_file, ".gz")))

      sample_file <- file.path(out_dir, sprintf("sample.%s.final", stage))
      writeLines(final_samples, sample_file)
    }
  }

  audit_v2_rows[[as.character(downsample_seed)]] <- audit_v2_seed
  format_cmd <- sprintf("python3 %s %s %s", file.path(getwd(), "format.py"), out_dir, out_dir)
  format_exit <- system(format_cmd)
  if (format_exit != 0) warning("format.py failed with exit code ", format_exit)
}

# --- audit ---

audit_v2 <- rbindlist(audit_v2_rows, fill = TRUE)
fwrite(audit_v2, audit_path, sep = "\t", quote = FALSE)
