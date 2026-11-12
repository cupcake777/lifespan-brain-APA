#!/usr/bin/env Rscript
# Preprocess expression data for stage-specific eQTL and apaQTL mapping.
# Sample QC, sex validation, gene filtering, TMM-CPM, outlier screen, ComBat, inverse normal transform.

suppressMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(WGCNA)
  library(sva)
  library(dplyr)
  library(PCAForQTL)
})

# --- helpers ---
load_data <- function(filepath) {
  df <- fread(filepath, data.table = FALSE)
  ids <- as.character(df[, 1])
  if (anyDuplicated(ids)) {
    stop(sprintf("Error: %d duplicate row ids in %s (first: %s)",
                 sum(duplicated(ids)), filepath, ids[anyDuplicated(ids)]))
  }
  rownames(df) <- ids
  df[, -1, drop = FALSE]
}

log_message <- function(message, level = "INFO") {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, message))
}

run_or_stop <- function(cmd, what) {
  message(sprintf("Running: %s", cmd))
  code <- system(cmd)
  if (code != 0) stop(sprintf("%s failed with exit code %d", what, code))
  invisible(code)
}

# --- parameters ---
set.seed(2026)  # kmeans (sex thresholds) and runBE are stochastic
WGCNA::allowWGCNAThreads(nThreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))

RIN_MIN               <- 5      # minimum RIN to keep a sample
EXPR_MIN_FRAC         <- 0.20   # fraction of samples a gene must be expressed in
TPM_MIN               <- 0.1    # TPM threshold for "expressed"
COUNT_MIN             <- 6      # raw count threshold for "expressed"
OUTLIER_Z             <- 3      # |Z| of network connectivity flagging an outlier
MIN_SAMPLES_PER_STAGE <- 5      # stages below this are skipped entirely
K_PC_MIN              <- 1      # lower bound on expression PCs
K_PC_MAX              <- 50     # upper bound on expression PCs
RUNBE_B               <- 20     # permutations for PCAForQTL::runBE
KNOWN_COV_R2_CUTOFF   <- 0.9    # drop known covariates this well explained by PCs
COV_MIN_VARIANCE      <- 1e-10  # drop constant covariates (rank deficiency)
COMBAT_REF_BATCH      <- "batch9"  # the only batch spanning both life stages
GENO_PCA_MAF          <- 0.05
GENO_PCA_INDEP_WINDOW <- 50
GENO_PCA_INDEP_STEP   <- 5
GENO_PCA_INDEP_R2     <- 0.2
GENO_N_PCS            <- 5

# --- paths ---
root <- Sys.getenv("APA_ROOT")
count_file <- file.path(root, "raw", "all.912.expr.count")
tpm_file   <- file.path(root, "raw", "all.912.expr.tpm")
gtf_file   <- file.path(root, "ref", "gencode.v40.annotation.gtf.gz")
geno_file  <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/raw_data/genotype")
meta_file  <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/raw_data/all.827.meta")
out_dir <- file.path(root, "data", "04_qtl")
fig_dir <- file.path(root, "figures", "04_qtl")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
pre_dir <- file.path(out_dir, "eqtl", "pre_input")
if (!dir.exists(pre_dir)) dir.create(pre_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --- load data ---
log_message("Loading input data files")
count <- load_data(count_file)
tpm <- load_data(tpm_file)

# Align TPM to counts by gene id; the expression filter is derived from counts.
if (!identical(rownames(count), rownames(tpm))) {
  shared_genes <- intersect(rownames(count), rownames(tpm))
  log_message(sprintf("Re-aligning TPM to counts by gene id: %d -> %d genes",
                      nrow(count), length(shared_genes)), "WARNING")
  count <- count[shared_genes, , drop = FALSE]
  tpm <- tpm[shared_genes, , drop = FALSE]
}

meta <- load_data(meta_file)

# GENCODE GTF, gene rows only; one entry per unversioned gene id.
gtf <- fread(gtf_file, sep = "\t", header = FALSE, data.table = FALSE, skip = "#", quote = "")
colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end",
                        "score", "strand", "frame", "attribute")
gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
annotation_df <- data.frame(
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
annotation_df$gene_id_base <- sub("\\.[0-9]+$", "", annotation_df$gene_id)
annotation_df <- annotation_df[!duplicated(annotation_df$gene_id_base), ]

log_message(sprintf("Loaded: %d genes x %d samples", nrow(count), ncol(count)))

# --- sample QC ---
required_cols <- c("RIN", "sex", "LifeStage", "batch")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Missing required columns in metadata: %s", paste(missing_cols, collapse = ", ")))
}

# which() keeps samples with NA RIN out.
meta$RIN <- as.numeric(meta$RIN)
meta <- meta[which(meta$RIN >= RIN_MIN), , drop = FALSE]

common_samples <- intersect(rownames(meta), colnames(count))
meta <- meta[common_samples, , drop = FALSE]
count <- count[, common_samples, drop = FALSE]
tpm <- tpm[, common_samples, drop = FALSE]

# --- sex validation (XIST / UTY) ---
meta$sex_imputed <- meta$sex
xist_id <- grep("^ENSG00000229807", rownames(count), value = TRUE)
uty_id <- grep("^ENSG00000183878", rownames(count), value = TRUE)
common_samples_meta <- intersect(rownames(meta), colnames(count))

if (length(xist_id) == 0 || length(uty_id) == 0) {
  warning("XIST or UTY gene not found in count matrix. Skipping sex validation.")
} else if (length(common_samples_meta) == 0) {
  warning("No common samples between meta and expression data. Skipping sex validation.")
} else {
  xist_expr <- setNames(as.numeric(count[xist_id[1], ]), colnames(count))
  uty_expr <- setNames(as.numeric(count[uty_id[1], ]), colnames(count))

  meta <- meta[common_samples_meta, , drop = FALSE]
  xist_expr <- xist_expr[common_samples_meta]
  uty_expr <- uty_expr[common_samples_meta]

  is_unknown <- is.na(meta$sex_imputed) | meta$sex_imputed %in% c("unknown", "")

  # k-means midpoint as the expression threshold
  km_xist <- kmeans(log1p(xist_expr), centers = 2, nstart = 25)
  km_uty <- kmeans(log1p(uty_expr), centers = 2, nstart = 25)
  xist_threshold <- expm1(mean(km_xist$centers))
  uty_threshold <- expm1(mean(km_uty$centers))

  log_message(sprintf("XIST threshold: %.2f (midpoint)", xist_threshold))
  log_message(sprintf("UTY threshold: %.2f (midpoint)", uty_threshold))

  # One precedence rule for imputation and mismatch: UTY positive is Male,
  # otherwise XIST positive is Female, both low is Ambiguous.
  xist_pos <- xist_expr > xist_threshold
  uty_pos  <- uty_expr > uty_threshold

  imputed_sex <- rep(NA_character_, length(xist_pos))
  names(imputed_sex) <- names(xist_pos)
  imputed_sex[xist_pos & !uty_pos] <- "Female"
  imputed_sex[uty_pos] <- "Male"
  is_ambiguous <- !xist_pos & !uty_pos

  # Impute unknown sex; undetermined stays NA and is removed below as ambiguous.
  meta$sex_imputed[is_unknown] <- imputed_sex[is_unknown]

  sex_mismatch <- rownames(meta)[!is_unknown & !is_ambiguous &
                                   !is.na(imputed_sex) & meta$sex != imputed_sex]
  ambiguous <- rownames(meta)[is_ambiguous]
  samples_to_filter <- unique(c(sex_mismatch, ambiguous))

  # plot
  sex_colors <- c(Male = "#3498db", Female = "#e74c3c", Unknown = "#7f8c8d")
  status_shapes <- c("Match" = 16, "Mismatch" = 4, "Unknown-Imputed" = 17, "Ambiguous" = 8)

  sex_df <- data.frame(
    Sample = rownames(meta),
    XIST = xist_expr,
    UTY = uty_expr,
    ReportedSex = ifelse(is_unknown, "Unknown", meta$sex),
    ImputedSex = ifelse(is.na(imputed_sex), "Ambiguous", imputed_sex),
    stringsAsFactors = FALSE
  )

  sex_df$Status <- "Match"
  sex_df$Status[sex_df$Sample %in% sex_mismatch] <- "Mismatch"
  sex_df$Status[sex_df$Sample %in% ambiguous] <- "Ambiguous"
  sex_df$Status[is_unknown] <- "Unknown-Imputed"

  p_scatter <- ggplot(sex_df, aes(XIST, UTY, color = ReportedSex, shape = Status)) +
    geom_point(alpha = 0.7, size = 3) +
    scale_color_manual(values = sex_colors, name = "Reported Sex") +
    scale_shape_manual(values = status_shapes, name = "Status") +
    geom_vline(xintercept = xist_threshold, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = uty_threshold, linetype = "dashed", alpha = 0.5) +
    labs(title = "Sex Validation: XIST vs UTY Expression", x = "XIST Expression", y = "UTY Expression") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "bottom")

  ggsave(file.path(fig_dir, "sex_validation_scatter.png"), p_scatter, width = 10, height = 8, dpi = 300)

  # sum(..., na.rm = TRUE): one NA sex must not turn the summary into NA
  n_total <- nrow(sex_df)
  counts <- c(n_total,
              sum(meta$sex == "Male", na.rm = TRUE),
              sum(meta$sex == "Female", na.rm = TRUE),
              sum(is_unknown),
              length(sex_mismatch),
              length(ambiguous),
              length(samples_to_filter))
  summary_df <- data.frame(
    Category = c("Total samples", "Reported Male", "Reported Female", "Unknown sex imputed",
                 "Sex mismatch (will filter)", "Ambiguous (will filter)", "Total to filter"),
    Count = counts,
    Percentage = round(100 * counts / n_total, 1)
  )
  summary_df$Percentage[1] <- 100
  write.csv(summary_df, file.path(out_dir, "sex_validation_summary.csv"), row.names = FALSE)

  if (length(samples_to_filter) > 0) {
    filtered_samples <- sex_df[sex_df$Sample %in% samples_to_filter,
                               c("Sample", "ReportedSex", "ImputedSex", "XIST", "UTY", "Status")]
    write.csv(filtered_samples, file.path(out_dir, "sex_filtered_samples.csv"), row.names = FALSE)
  }

  log_message(sprintf("Sex validation: %d samples processed", n_total))
  log_message(sprintf("  - Sex mismatch: %d, Ambiguous: %d, Total to filter: %d",
                      length(sex_mismatch), length(ambiguous), length(samples_to_filter)))

  meta <- meta[!(rownames(meta) %in% samples_to_filter), , drop = FALSE]

  # an undetermined sex would drop rows silently in model.matrix
  if (any(is.na(meta$sex_imputed))) {
    stop(sprintf("Error: %d retained samples still have an undetermined sex",
                 sum(is.na(meta$sex_imputed))))
  }
}

common_samples <- rownames(meta)
count <- count[, common_samples, drop = FALSE]
tpm <- tpm[, common_samples, drop = FALSE]

meta$sex <- meta$sex_imputed
meta$sex_imputed <- NULL

# --- gene filtering ---
log_message("Performing gene filtering")

keep <- rowSums(count) > 0
count <- count[keep, , drop = FALSE]
tpm <- tpm[keep, , drop = FALSE]

threshold <- EXPR_MIN_FRAC * ncol(count)
tpm_filter <- rowSums(tpm > TPM_MIN) >= threshold
count_filter <- rowSums(count >= COUNT_MIN) >= threshold

# whitelist chr1-22 so unplaced scaffolds cannot pass as autosomal
autosomes <- c(as.character(1:22), paste0("chr", 1:22))
autosomal_ids <- annotation_df$gene_id_base[annotation_df$chr %in% autosomes]
count_gene_ids_base <- sub("\\.[0-9]+$", "", rownames(count))

keep <- tpm_filter & count_filter & (count_gene_ids_base %in% autosomal_ids)
log_message(sprintf("After autosomal filtering: %d genes", sum(keep)))
count <- count[keep, , drop = FALSE]
log_message(sprintf("After gene filtering: %d genes x %d samples", nrow(count), ncol(count)))
if (nrow(count) == 0) stop("Error: gene filtering removed every gene")

# --- normalization ---
log_message("Performing TMM normalization")

dge <- DGEList(count)
dge <- calcNormFactors(dge, method = "TMM")
# normalized.lib.sizes = TRUE applies lib.size * norm.factors
count_norm <- as.data.frame(cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))

# --- outlier detection (network connectivity within life stage) ---
log_message("Performing outlier detection using network connectivity")

outlier_list <- list()
for (g in unique(na.omit(meta$LifeStage))) {
  samples <- rownames(meta)[which(meta$LifeStage == g)]
  if (length(samples) < 3) {
    log_message(sprintf("Skipping group %s: insufficient samples (%d)", g, length(samples)), "WARNING")
    next
  }

  adj <- adjacency(as.matrix(count_norm[, samples]), type = "signed", corFnc = "bicor")
  k <- fundamentalNetworkConcepts(adj)$Connectivity

  if (all(is.na(k))) {
    log_message(sprintf("Skipping group %s: all NA connectivity values", g), "WARNING")
    next
  }

  Z.K <- scale(k)[, 1]
  outlier_list[[g]] <- data.frame(
    Sample = samples,
    Z.K = Z.K,
    IsOutlier = abs(Z.K) > OUTLIER_Z,
    Group = g,
    stringsAsFactors = FALSE
  )
}

if (length(outlier_list) == 0) {
  warning("No valid groups for outlier detection")
  outliers <- data.frame()
} else {
  outliers <- do.call(rbind, outlier_list) %>%
    arrange(Group, Z.K) %>%
    mutate(PlotIndex = 1:n())

  p_outlier <- ggplot(outliers, aes(x = PlotIndex, y = Z.K, color = Group)) +
    geom_point(size = 2) +
    geom_hline(yintercept = c(-OUTLIER_Z, OUTLIER_Z), linetype = "dashed") +
    geom_text_repel(data = subset(outliers, IsOutlier), aes(label = Sample),
                    vjust = 0, hjust = 1.2, color = "red", size = 3) +
    labs(title = "Sample Outlier Detection",
         y = "Network Connectivity (Z-score)",
         x = "Sample Index") +
    theme_classic() +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "outlier_detection.png"), p_outlier, width = 10, height = 5)
  fwrite(outliers, file.path(out_dir, "outlier_detection.tsv"), sep = "\t")
}

# Exclude flagged samples; samples in unscreened groups are kept.
samples_drop <- if (nrow(outliers) > 0) outliers$Sample[outliers$IsOutlier] else character()
samples_keep <- setdiff(colnames(count_norm), samples_drop)

if (nrow(outliers) > 0 && nrow(outliers) < ncol(count_norm)) {
  log_message(sprintf("Outlier detection covered %d of %d samples; the remainder were not screened",
                      nrow(outliers), ncol(count_norm)), "WARNING")
}

count_clean <- count_norm[, samples_keep, drop = FALSE]
meta_clean <- meta[samples_keep, , drop = FALSE]
log_message(sprintf("Removed %d outlier samples", length(samples_drop)))
log_message(sprintf("After outlier removal: %d genes x %d samples", nrow(count_clean), ncol(count_clean)))

# --- batch correction (ComBat) ---
log_message("Starting batch correction and normalization pipeline")

meta_clean$batch <- as.factor(meta_clean$batch)
meta_clean$LifeStage <- as.factor(meta_clean$LifeStage)
meta_clean$sex <- as.factor(meta_clean$sex)

batch_vec <- meta_clean$batch

mod <- model.matrix(~ LifeStage, data = meta_clean)
if (nrow(mod) != ncol(count_clean)) {
  stop(sprintf("Error: ComBat design has %d rows but the expression matrix has %d samples (NA in LifeStage?)",
               nrow(mod), ncol(count_clean)))
}

ref_batch <- COMBAT_REF_BATCH
if (!ref_batch %in% levels(batch_vec)) {
  log_message(sprintf("Reference batch '%s' not present; running ComBat without a reference", ref_batch), "WARNING")
  ref_batch <- NULL
}
expr_clean <- suppressWarnings(
  ComBat(dat = as.matrix(count_clean), batch = batch_vec, mod = mod, ref.batch = ref_batch)
)
fwrite(as.data.frame(expr_clean), file.path(pre_dir, "after.combat.txt"), sep = '\t', quote = FALSE, row.names = TRUE)

# --- genotype alignment ---
if (!all(file.exists(paste0(geno_file, c(".pgen", ".pvar", ".psam"))))) {
  stop(sprintf("Error: no genotype pfile at %s (.pgen/.pvar/.psam)", geno_file))
}

geno_psam <- fread(paste0(geno_file, ".psam"), data.table = FALSE)
iid_col <- intersect(c("#IID", "IID"), colnames(geno_psam))[1]
geno_samples <- as.character(geno_psam[[iid_col]])

common_samples_global <- intersect(colnames(expr_clean), geno_samples)
log_message(sprintf("Global alignment: %d expression samples -> %d samples with genotype",
                    ncol(expr_clean), length(common_samples_global)))
if (length(common_samples_global) < MIN_SAMPLES_PER_STAGE) {
  stop(sprintf(paste("Error: only %d samples have both expression and genotype.",
                     "Check that sample ids match between %s and %s.psam"),
               length(common_samples_global), count_file, geno_file))
}

expr_clean <- expr_clean[, common_samples_global, drop = FALSE]
meta_clean <- meta_clean[common_samples_global, , drop = FALSE]

stages <- unique(na.omit(meta_clean$LifeStage))
log_message(sprintf("After filtering, processing %d stages: %s",
                    length(stages), paste(stages, collapse = ", ")))

# --- per-stage QTL inputs ---
# Chromosome naming is checked against the per-stage pvar that tensorQTL reads,
# since plink2 --make-pgen re-codes chr1 to 1.
stage_summary <- list()

for (stage in stages) {
  stage_samples <- rownames(meta_clean)[which(meta_clean$LifeStage == stage)]
  if (length(stage_samples) < MIN_SAMPLES_PER_STAGE) {
    log_message(sprintf("[%s] Skipping: Only %d samples available (< %d minimum required)",
                        stage, length(stage_samples), MIN_SAMPLES_PER_STAGE))
    next
  }

  sub_expr_clean <- expr_clean[, stage_samples, drop = FALSE]
  sub_meta <- meta_clean[stage_samples, , drop = FALSE]

  # rank-based inverse normal transform per gene
  sub_expr_norm <- t(apply(as.matrix(sub_expr_clean), 1, function(x) {
    if (length(unique(x)) <= 1) {
      # setNames keeps sample names when apply() receives an unnamed result
      return(setNames(rep(0, length(x)), names(x)))
    }
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    n <- length(x)
    qnorm(ranks / (n + 1))
  }))
  sub_expr_norm <- sub_expr_norm[complete.cases(sub_expr_norm), , drop = FALSE]
  out_filename <- sprintf("pre.phenotype.combat.int.%s.txt", stage)
  fwrite(as.data.frame(sub_expr_norm), file.path(pre_dir, out_filename), sep = '\t', quote = FALSE, row.names = TRUE)

  # expression PCs and runBE
  vars_norm <- apply(sub_expr_norm, 1, var)
  sub_expr_for_pca <- sub_expr_norm[vars_norm > 1e-10, , drop = FALSE]

  if (ncol(sub_expr_for_pca) < 2 || nrow(sub_expr_for_pca) < 2) {
    warning(sprintf("[%s] Insufficient samples or phenotypes for covariate analysis. Skipping.", stage))
    next
  }

  sub_pca_final <- prcomp(t(sub_expr_for_pca), center = TRUE, scale. = TRUE)
  K_BE <- suppressWarnings(
    as.integer(PCAForQTL::runBE(t(sub_expr_for_pca), B = RUNBE_B)$numOfPCsChosen)
  )
  if (is.na(K_BE)) {
    log_message(sprintf("[%s] runBE returned NA; falling back to %d PCs", stage, K_PC_MAX), "WARNING")
    K_BE <- K_PC_MAX
  }
  K_final <- min(max(K_BE, K_PC_MIN), K_PC_MAX, ncol(sub_pca_final$x))
  log_message(sprintf("[%s] runBE chose %d expression PCs; using %d (min %d, max %d)",
                      stage, K_BE, K_final, K_PC_MIN, K_PC_MAX))

  # numeric known covariates; sex is the only categorical after batch and LifeStage are removed
  known_covs <- as.data.frame(lapply(sub_meta, function(x) {
    as.numeric(if (is.character(x)) as.factor(x) else x)
  }))
  rownames(known_covs) <- rownames(sub_meta)
  known_covs <- known_covs[, !colnames(known_covs) %in% c("batch", "sample_id", "LifeStage"), drop = FALSE]

  # stage genotype and genotype PCs (plink2)
  pgen_prefix <- file.path(pre_dir, paste0("genotype.", tolower(stage)))
  keep_file <- paste0(pgen_prefix, ".keep")
  fwrite(data.frame("#IID" = rownames(known_covs), stringsAsFactors = FALSE, check.names = FALSE),
         keep_file, sep = "\t", quote = FALSE)

  pca_n <- min(GENO_N_PCS, nrow(known_covs) - 1)

  run_or_stop(
    sprintf("plink2 --pfile %s --keep %s --make-pgen --out %s",
            shQuote(geno_file), shQuote(keep_file), shQuote(pgen_prefix)),
    sprintf("plink2 --make-pgen for %s", pgen_prefix)
  )

  prune_prefix <- paste0(pgen_prefix, ".pca_pruned")
  run_or_stop(
    sprintf("plink2 --pfile %s --maf %.4f --indep-pairwise %d %d %.3f --out %s",
            shQuote(pgen_prefix), GENO_PCA_MAF,
            GENO_PCA_INDEP_WINDOW, GENO_PCA_INDEP_STEP, GENO_PCA_INDEP_R2,
            shQuote(prune_prefix)),
    sprintf("plink2 LD pruning for %s", prune_prefix)
  )

  prune_in_file <- paste0(prune_prefix, ".prune.in")
  if (!file.exists(prune_in_file) || file.info(prune_in_file)$size == 0) {
    stop(sprintf("plink2 LD pruning produced no variants for PCA: %s", prune_in_file))
  }

  pca_prefix <- paste0(pgen_prefix, ".pca")
  run_or_stop(
    sprintf("plink2 --pfile %s --extract %s --pca %d --out %s",
            shQuote(pgen_prefix), shQuote(prune_in_file), pca_n, shQuote(pca_prefix)),
    sprintf("plink2 PCA for %s", pca_prefix)
  )

  eigenvec <- fread(paste0(pca_prefix, ".eigenvec"), data.table = FALSE)
  iid_col <- intersect(c("#IID", "IID"), colnames(eigenvec))[1]
  rownames(eigenvec) <- as.character(eigenvec[[iid_col]])
  geno_pcs <- eigenvec[, grep("^PC[0-9]+$", colnames(eigenvec), value = TRUE), drop = FALSE]
  colnames(geno_pcs) <- paste0("genoPC", seq_len(ncol(geno_pcs)))

  # chromosome coding of the stage pvar (the file tensorQTL is given)
  con <- file(paste0(pgen_prefix, ".pvar"), open = "rt")
  stage_geno_chroms <- character()
  n_left <- 20000
  repeat {
    lines <- readLines(con, n = 5000, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!startsWith(lines, "##") & !startsWith(lines, "#CHROM")]
    if (length(lines)) {
      stage_geno_chroms <- unique(c(stage_geno_chroms, sub("\t.*$", "", lines)))
    }
    n_left <- n_left - 5000
    if (n_left <= 0) break
  }
  close(con)
  stage_geno_has_chr_prefix <- length(stage_geno_chroms) > 0 &&
    all(startsWith(stage_geno_chroms, "chr"))

  common_cov_samples <- intersect(rownames(known_covs), rownames(geno_pcs))
  if (length(common_cov_samples) < MIN_SAMPLES_PER_STAGE) {
    stop(sprintf("[%s] Error: only %d samples remain after joining genotype PCs",
                 stage, length(common_cov_samples)))
  }
  known_covs <- cbind(
    known_covs[common_cov_samples, , drop = FALSE],
    geno_pcs[common_cov_samples, , drop = FALSE]
  )
  sub_expr_norm <- sub_expr_norm[, common_cov_samples, drop = FALSE]
  sub_pca_final$x <- sub_pca_final$x[common_cov_samples, , drop = FALSE]

  # drop constant covariates (rank deficiency in tensorQTL)
  cov_var <- apply(known_covs, 2, var, na.rm = TRUE)
  const_covs <- names(which(!is.finite(cov_var) | cov_var <= COV_MIN_VARIANCE))
  if (length(const_covs) > 0) {
    log_message(sprintf("[%s] Dropping %d constant covariate(s): %s",
                        stage, length(const_covs), paste(const_covs, collapse = ", ")), "WARNING")
    known_covs <- known_covs[, setdiff(colnames(known_covs), const_covs), drop = FALSE]
  }
  if (anyNA(known_covs)) {
    stop(sprintf("[%s] Error: %d missing values in the known covariates",
                 stage, sum(is.na(known_covs))))
  }

  PCsTop <- as.data.frame(sub_pca_final$x)[rownames(known_covs), 1:K_final, drop = FALSE]
  covsFiltered <- PCAForQTL::filterKnownCovariates(known_covs, PCsTop,
                                                   unadjustedR2_cutoff = KNOWN_COV_R2_CUTOFF)
  covsToUse <- cbind(covsFiltered, PCsTop)
  log_message(sprintf("[%s] Covariates: %d known kept (of %d) + %d expression PCs = %d total",
                      stage, ncol(covsFiltered), ncol(known_covs), ncol(PCsTop), ncol(covsToUse)))

  out_filename <- sprintf("covariates_for_qtl.%s.txt", stage)
  fwrite(as.data.frame(t(covsToUse)), file.path(pre_dir, out_filename),
         sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

  # phenotype BED: 1 bp strand-aware TSS window (GTEx / tensorQTL convention)
  gene_ids <- rownames(sub_expr_norm)
  gene_ids_base <- sub("\\.[0-9]+$", "", gene_ids)
  gene_coords <- annotation_df[match(gene_ids_base, annotation_df$gene_id_base),
                               c("chr", "start", "end", "strand", "gene_id", "gene_id_base")]

  missing_coords <- is.na(gene_coords$chr)
  if (sum(missing_coords) > 0) {
    warning(sprintf("[%s] %d genes without coordinates found, removing them",
                    stage, sum(missing_coords)))
    gene_ids <- gene_ids[!missing_coords]
    gene_ids_base <- gene_ids_base[!missing_coords]
    gene_coords <- gene_coords[!missing_coords, ]
    sub_expr_norm <- sub_expr_norm[gene_ids, , drop = FALSE]
  }

  final_samples <- rownames(covsToUse)

  # unknown strand is treated as "+" but reported
  unknown_strand <- !gene_coords$strand %in% c("+", "-")
  if (any(unknown_strand)) {
    warning(sprintf("%d genes with unusable strand; TSS taken as the start coordinate",
                    sum(unknown_strand)))
  }
  is_minus <- !is.na(gene_coords$strand) & gene_coords$strand == "-"
  tss <- ifelse(is_minus, gene_coords$end, gene_coords$start)
  bed_pos <- data.frame(start = as.integer(tss) - 1L, end = as.integer(tss))

  format_data <- data.frame(
    `#Chr` = gene_coords$chr,
    start = bed_pos$start,
    end = bed_pos$end,
    pid = gene_ids,
    sub_expr_norm[, final_samples, drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # match chromosome naming to the stage genotype before sorting
  pheno_has_chr_prefix <- all(startsWith(as.character(format_data$`#Chr`), "chr"))
  target_has_chr_prefix <- if (length(stage_geno_chroms) > 0) stage_geno_has_chr_prefix else pheno_has_chr_prefix
  if (target_has_chr_prefix && !pheno_has_chr_prefix) {
    log_message(sprintf("[%s] Adding 'chr' prefix to phenotype chromosomes to match the stage genotype", stage), "WARNING")
    format_data$`#Chr` <- paste0("chr", format_data$`#Chr`)
  } else if (!target_has_chr_prefix && pheno_has_chr_prefix) {
    log_message(sprintf("[%s] Stripping 'chr' prefix from phenotype chromosomes to match the stage genotype", stage), "WARNING")
    format_data$`#Chr` <- sub("^chr", "", format_data$`#Chr`)
  }

  # sort by chromosome then position (tensorQTL requires sorted BED)
  chrom_core <- c(1:22, "X", "Y", "M")
  chrom_order <- if (target_has_chr_prefix) paste0("chr", chrom_core) else as.character(chrom_core)
  format_data$sort_chr <- factor(format_data$`#Chr`, levels = chrom_order, ordered = TRUE)
  if (anyNA(format_data$sort_chr)) {
    bad <- unique(as.character(format_data$`#Chr`)[is.na(format_data$sort_chr)])
    stop(sprintf("[%s] Error: %d genes on unrecognised chromosomes (%s); sorting would be wrong",
                 stage, sum(is.na(format_data$sort_chr)),
                 paste(head(bad, 10), collapse = ", ")))
  }
  format_data <- format_data[order(format_data$sort_chr, format_data$start, format_data$end), ]
  format_data$sort_chr <- NULL
  rownames(format_data) <- NULL

  bed_file <- file.path(pre_dir, sprintf("phenotype.%s.bed", stage))
  # header kept: tensorQTL expects the first column named '#Chr'
  fwrite(format_data, bed_file, sep = '\t', quote = FALSE, col.names = TRUE)
  run_or_stop(paste("bgzip -f", shQuote(bed_file)), sprintf("bgzip for %s", bed_file))
  run_or_stop(paste("tabix -p bed -f", shQuote(paste0(bed_file, ".gz"))),
              sprintf("tabix for %s.gz", bed_file))

  sample_file <- file.path(pre_dir, sprintf("sample.%s.final", stage))
  writeLines(final_samples, sample_file)

  stage_summary[[as.character(stage)]] <- data.frame(
    stage = as.character(stage),
    n_samples = length(final_samples),
    n_genes = nrow(format_data),
    n_known_covs = ncol(covsFiltered),
    n_expr_pcs = ncol(PCsTop),
    runBE_K = K_BE,
    coord_mode = "tss",
    stringsAsFactors = FALSE
  )
}

if (length(stage_summary) == 0) {
  stop("Error: no stage produced QTL inputs")
}

summary_tbl <- do.call(rbind, stage_summary)
fwrite(summary_tbl, file.path(out_dir, "preprocess_summary.tsv"), sep = "\t")
print(summary_tbl)
