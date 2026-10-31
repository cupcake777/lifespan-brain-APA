# Expression preprocessing for eQTL: RIN filter, sex validation, gene filter, TMM,
# outlier screen, ComBat, inverse-normal transform, genotype PCs, covariates and BED.
# Output: data/04_qtl/<PRE_PREFIX>* (default prefix eqtl.); figures/04_qtl/.

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

# --- settings ---
RNG_SEED <- 2026
set.seed(RNG_SEED)
WGCNA::allowWGCNAThreads(nThreads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")))

RIN_MIN <- 5
EXPR_MIN_FRAC <- 0.20
TPM_MIN <- 0.1
COUNT_MIN <- 6
OUTLIER_Z <- 3
MIN_SAMPLES_PER_STAGE <- 5
K_PC_MIN <- 1          # 10 reproduces the 2026-05-06 covariate files
K_PC_MAX <- 50
RUNBE_B <- 20
KNOWN_COV_R2_CUTOFF <- 0.9
COV_MIN_VARIANCE <- 1e-10
COMBAT_REF_BATCH <- "batch9"  # only batch spanning both life stages
PHENO_COORD_MODE <- "tss"     # tss: 1 bp TSS window (GTEx); gene_body: 2026-05-06 files
GENO_N_PCS <- 5
GENO_PCA_MAF <- 0.05
GENO_PCA_INDEP_WINDOW <- 50
GENO_PCA_INDEP_STEP <- 5
GENO_PCA_INDEP_R2 <- 0.2
stopifnot(PHENO_COORD_MODE %in% c("tss", "gene_body"))

# --- paths ---
APA_ROOT <- Sys.getenv("APA_ROOT")
APA_HPC_HOME <- Sys.getenv("APA_HPC_HOME")
stopifnot(nzchar(APA_ROOT), nzchar(APA_HPC_HOME))
RAW <- file.path(APA_HPC_HOME, "project", "raw_data")
REF <- file.path(APA_HPC_HOME, "project", ".ref", "gene")
count_file <- file.path(RAW, "all.912.count")
tpm_file <- file.path(RAW, "all.912.tpm")
meta_file <- file.path(RAW, "all.827.meta")
geno_file <- file.path(RAW, "genotype")
gtf_file <- file.path(REF, "gencode.v40.annotation.gtf.gz")

out_dir <- file.path(APA_ROOT, "data", "04_qtl")
fig_dir <- file.path(APA_ROOT, "figures", "04_qtl")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
QP <- Sys.getenv("PRE_PREFIX", unset = "eqtl.")

# --- helpers used more than once ---
load_data <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Error: File not found: %s", filepath))
  }
  df <- fread(filepath, data.table = FALSE)
  if (ncol(df) < 2) {
    stop(sprintf("Error: File %s has insufficient columns", filepath))
  }
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
  if (code != 0) {
    stop(sprintf("%s failed with exit code %d", what, code))
  }
  invisible(code)
}

# --- load counts, TPM, metadata, annotation ---
log_message("Loading input data files")
count <- load_data(count_file)
if (file.exists(tpm_file)) {
  tpm <- load_data(tpm_file)
} else {
  log_message("TPM file not found, creating from count data", "WARNING")
  tpm <- count
}

# Align TPM to counts by gene id, not position: the expression filter is
# derived from counts and applied to both matrices.
if (!identical(rownames(count), rownames(tpm))) {
  shared_genes <- intersect(rownames(count), rownames(tpm))
  if (length(shared_genes) == 0) {
    stop("Error: count and TPM matrices share no gene ids")
  }
  log_message(sprintf("Re-aligning TPM to counts by gene id: %d -> %d genes",
                      nrow(count), length(shared_genes)), "WARNING")
  count <- count[shared_genes, , drop = FALSE]
  tpm <- tpm[shared_genes, , drop = FALSE]
}

meta <- load_data(meta_file)

# GENCODE gene records, one row per version-free gene id (first occurrence kept).
stopifnot(file.exists(gtf_file))
gtf <- fread(gtf_file, sep = "\t", header = FALSE, data.table = FALSE, skip = "#", quote = "")
if (ncol(gtf) < 9) {
  stop(sprintf("Error: GTF annotation has fewer than 9 columns: %s", gtf_file))
}
colnames(gtf)[1:9] <- c("chr", "source", "feature", "start", "end",
                        "score", "strand", "frame", "attribute")
gtf <- gtf[gtf$feature == "gene", , drop = FALSE]
annotation_df <- data.frame(
  chr = gtf$chr,
  start = as.integer(gtf$start),
  end = as.integer(gtf$end),
  strand = gtf$strand,
  gene_id = sub('.*gene_id "([^"]+)".*', "\\1", gtf$attribute),
  stringsAsFactors = FALSE
)
annotation_df$gene_id_base <- sub("\\.[0-9]+$", "", annotation_df$gene_id)
annotation_df <- annotation_df[!duplicated(annotation_df$gene_id_base), ]
log_message(sprintf("Loaded annotation: %d genes", nrow(annotation_df)))
log_message(sprintf("Loaded: %d genes x %d samples", nrow(count), ncol(count)))

# --- sample QC: required columns, RIN filter, alignment ---
required_cols <- c("RIN", "sex", "LifeStage", "batch")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Error: Missing required columns in metadata: %s", paste(missing_cols, collapse = ", ")))
}

meta$RIN <- as.numeric(meta$RIN)
n_rin_missing <- sum(is.na(meta$RIN))
if (n_rin_missing > 0) {
  log_message(sprintf("%d samples have missing RIN and will be dropped", n_rin_missing), "WARNING")
}
meta <- meta[which(meta$RIN >= RIN_MIN), , drop = FALSE]
log_message(sprintf("Filtered by RIN >= %g: %d samples remaining", RIN_MIN, nrow(meta)))

common_samples <- intersect(rownames(meta), colnames(count))
if (length(common_samples) == 0) {
  stop("Error: No common samples between metadata and count matrix")
}
meta <- meta[common_samples, , drop = FALSE]
count <- count[, common_samples, drop = FALSE]
tpm <- tpm[, common_samples, drop = FALSE]
log_message(sprintf("After sample alignment: %d genes x %d samples", nrow(count), ncol(count)))

# --- sex validation: XIST/UTY k-means thresholds, one precedence rule ---
meta$sex_imputed <- meta$sex
xist_id <- grep("^ENSG00000229807", rownames(count), value = TRUE)
uty_id <- grep("^ENSG00000183878", rownames(count), value = TRUE)

if (length(xist_id) == 0 || length(uty_id) == 0) {
  warning("XIST or UTY gene not found in count matrix. Skipping sex validation.")
} else {
  xist_id <- xist_id[1]
  uty_id <- uty_id[1]
  xist_expr <- setNames(as.numeric(count[xist_id, ]), colnames(count))
  uty_expr <- setNames(as.numeric(count[uty_id, ]), colnames(count))

  common_samples_meta <- intersect(rownames(meta), names(xist_expr))
  if (length(common_samples_meta) == 0) {
    warning("No common samples between meta and expression data. Skipping sex validation.")
  } else {
    meta <- meta[common_samples_meta, , drop = FALSE]
    xist_expr <- xist_expr[common_samples_meta]
    uty_expr <- uty_expr[common_samples_meta]

    is_unknown <- is.na(meta$sex_imputed) | meta$sex_imputed %in% c("unknown", "")

    # k-means midpoint between the two clusters is the expression threshold
    km_xist <- kmeans(log1p(xist_expr), centers = 2, nstart = 25)
    km_uty <- kmeans(log1p(uty_expr), centers = 2, nstart = 25)
    xist_threshold <- expm1(mean(km_xist$centers))
    uty_threshold <- expm1(mean(km_uty$centers))
    log_message(sprintf("XIST threshold: %.2f (midpoint)", xist_threshold))
    log_message(sprintf("UTY threshold: %.2f (midpoint)", uty_threshold))

    # UTY positive -> Male; else XIST positive -> Female; both low -> ambiguous
    xist_pos <- xist_expr > xist_threshold
    uty_pos <- uty_expr > uty_threshold
    imputed_sex <- rep(NA_character_, length(xist_pos))
    names(imputed_sex) <- names(xist_pos)
    imputed_sex[xist_pos & !uty_pos] <- "Female"
    imputed_sex[uty_pos] <- "Male"
    is_ambiguous <- !xist_pos & !uty_pos

    meta$sex_imputed[is_unknown] <- imputed_sex[is_unknown]

    sex_mismatch <- rownames(meta)[!is_unknown & !is_ambiguous &
                                     !is.na(imputed_sex) & meta$sex != imputed_sex]
    ambiguous <- rownames(meta)[is_ambiguous]
    samples_to_filter <- unique(c(sex_mismatch, ambiguous))

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
      scale_color_manual(values = c(Male = "#3498db", Female = "#e74c3c", Unknown = "#7f8c8d"),
                         name = "Reported Sex") +
      scale_shape_manual(values = c("Match" = 16, "Mismatch" = 4, "Unknown-Imputed" = 17, "Ambiguous" = 8),
                         name = "Status") +
      geom_vline(xintercept = xist_threshold, linetype = "dashed", alpha = 0.5) +
      geom_hline(yintercept = uty_threshold, linetype = "dashed", alpha = 0.5) +
      labs(title = "Sex Validation: XIST vs UTY Expression", x = "XIST Expression", y = "UTY Expression") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
    ggsave(file.path(fig_dir, paste0(QP, "sex_validation_scatter.png")), p_scatter,
           width = 10, height = 8, dpi = 300)

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
    write.csv(summary_df, file.path(out_dir, paste0(QP, "sex_validation_summary.csv")), row.names = FALSE)

    if (length(samples_to_filter) > 0) {
      filtered_samples <- sex_df[sex_df$Sample %in% samples_to_filter,
                                 c("Sample", "ReportedSex", "ImputedSex", "XIST", "UTY", "Status")]
      write.csv(filtered_samples, file.path(out_dir, paste0(QP, "sex_filtered_samples.csv")), row.names = FALSE)
    }

    log_message(sprintf("Sex validation: %d samples processed", n_total))
    log_message(sprintf("  - Sex mismatch: %d, Ambiguous: %d, Total to filter: %d",
                        length(sex_mismatch), length(ambiguous), length(samples_to_filter)))

    if (length(samples_to_filter) > 0) {
      log_message(sprintf("Removing %d samples with sex validation issues", length(samples_to_filter)))
      meta <- meta[!(rownames(meta) %in% samples_to_filter), , drop = FALSE]
    }

    if (any(is.na(meta$sex_imputed))) {
      stop(sprintf("Error: %d retained samples still have an undetermined sex",
                   sum(is.na(meta$sex_imputed))))
    }
  }
}

common_samples <- rownames(meta)
count <- count[, common_samples, drop = FALSE]
tpm <- tpm[, common_samples, drop = FALSE]
meta$sex <- meta$sex_imputed
meta$sex_imputed <- NULL

# --- gene filtering: zero counts, expression thresholds, autosomes ---
log_message("Performing gene filtering")
keep <- rowSums(count) > 0
count <- count[keep, , drop = FALSE]
tpm <- tpm[keep, , drop = FALSE]
log_message(sprintf("After removing zero-count genes: %d genes", nrow(count)))

threshold <- EXPR_MIN_FRAC * ncol(count)
tpm_filter <- rowSums(tpm > TPM_MIN) >= threshold
count_filter <- rowSums(count >= COUNT_MIN) >= threshold
log_message(sprintf("TPM filter: %d genes, Count filter: %d genes", sum(tpm_filter), sum(count_filter)))

# whitelist chr1-22 rather than blacklisting X/Y/M, so unplaced scaffolds cannot pass
autosomes <- c(as.character(1:22), paste0("chr", 1:22))
autosomal_ids <- annotation_df$gene_id_base[annotation_df$chr %in% autosomes]
count_gene_ids_base <- sub("\\.[0-9]+$", "", rownames(count))
keep <- tpm_filter & count_filter & (count_gene_ids_base %in% autosomal_ids)
log_message(sprintf("After autosomal filtering: %d genes", sum(keep)))
count <- count[keep, , drop = FALSE]
log_message(sprintf("After gene filtering: %d genes x %d samples", nrow(count), ncol(count)))
if (nrow(count) == 0) stop("Error: gene filtering removed every gene")

# --- TMM normalisation (log CPM) ---
log_message("Performing TMM normalization")
dge <- DGEList(count)
dge <- calcNormFactors(dge, method = "TMM")
count_norm <- as.data.frame(cpm(dge, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1))

# --- outlier screen: signed bicor network connectivity, Z per life stage ---
log_message("Performing outlier detection using network connectivity")
outlier_list <- list()
life_stages <- unique(na.omit(meta$LifeStage))
for (g in life_stages) {
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
  outliers <- do.call(rbind, outlier_list)
  outliers <- outliers %>%
    arrange(Group, Z.K) %>%
    mutate(PlotIndex = 1:n())

  p_outlier <- ggplot(outliers, aes(x = PlotIndex, y = Z.K, color = Group)) +
    geom_point(size = 2) +
    geom_hline(yintercept = c(-OUTLIER_Z, OUTLIER_Z), linetype = "dashed") +
    geom_text_repel(data = subset(outliers, IsOutlier), aes(label = Sample),
                    vjust = 0, hjust = 1.2, color = "red", size = 3) +
    labs(title = "Sample Outlier Detection", y = "Network Connectivity (Z-score)", x = "Sample Index") +
    theme_classic() +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0(QP, "outlier_detection.png")), p_outlier, width = 10, height = 5)
  fwrite(outliers, file.path(out_dir, paste0(QP, "outlier_detection.tsv")), sep = "\t")
}

# remove by exclusion: a sample is only lost if it was flagged, never because its
# life-stage group was skipped above
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

# --- ComBat (batch9 reference) ---
log_message("Starting batch correction and normalization pipeline")
meta_final <- meta_clean
meta_final$batch <- as.factor(meta_final$batch)
meta_final$LifeStage <- as.factor(meta_final$LifeStage)
meta_final$sex <- as.factor(meta_final$sex)
meta_clean <- meta_final

# batch1-4 are Postnatal-only and batch5-8 Prenatal-only in this cohort; batch9 is the
# only batch that separates stage from batch, so the cross-stage adjustment rests on it
batch_vec <- as.factor(meta_clean$batch)
message("Batch x LifeStage structure:")
print(table(batch_vec, meta_clean$LifeStage))
pure_batches <- names(which(rowSums(table(batch_vec, meta_clean$LifeStage) > 0) == 1))
if (length(pure_batches) > 0) {
  message(sprintf("  %d of %d batches contain a single life stage: %s",
                  length(pure_batches), nlevels(batch_vec), paste(pure_batches, collapse = ", ")))
}

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
message(sprintf("Running ComBat [reference batch: %s]...",
                if (is.null(ref_batch)) "none" else ref_batch))
expr_clean <- suppressWarnings(
  ComBat(dat = as.matrix(count_clean), batch = batch_vec, mod = mod, ref.batch = ref_batch)
)
fwrite(as.data.frame(expr_clean), file.path(out_dir, paste0(QP, "after.combat.txt")),
       sep = "\t", quote = FALSE, row.names = TRUE)

# --- genotype alignment: keep samples that have a genotype ---
if (!all(file.exists(paste0(geno_file, c(".pgen", ".pvar", ".psam"))))) {
  stop(sprintf("Error: no genotype pfile at %s (.pgen/.pvar/.psam)", geno_file))
}
message("Performing global genotype alignment...")
geno_psam_file <- paste0(geno_file, ".psam")
geno_psam <- fread(geno_psam_file, data.table = FALSE)
iid_col <- intersect(c("#IID", "IID"), colnames(geno_psam))[1]
if (is.na(iid_col)) {
  iid_col <- colnames(geno_psam)[1]
}
geno_samples <- as.character(geno_psam[[iid_col]])

common_samples_global <- intersect(colnames(expr_clean), geno_samples)
message(sprintf("Global alignment: %d expression samples -> %d samples with genotype",
                ncol(expr_clean), length(common_samples_global)))
if (length(common_samples_global) < MIN_SAMPLES_PER_STAGE) {
  stop(sprintf(paste("Error: only %d samples have both expression and genotype.",
                     "Check that sample ids match between %s and %s.psam"),
               length(common_samples_global), count_file, geno_file))
}
expr_clean <- expr_clean[, common_samples_global, drop = FALSE]
meta_clean <- meta_clean[common_samples_global, , drop = FALSE]
stages <- unique(na.omit(meta_clean$LifeStage))
message(sprintf("After filtering, processing %d stages: %s", length(stages), paste(stages, collapse = ", ")))

# --- per life stage: INT phenotype, genotype PCs, covariates, BED ---
# Chromosome naming is checked against the stage pvar that plink2 writes, because
# tensorQTL reads that file; plink2 --make-pgen re-codes chr1 to 1.
stage_summary <- list()

for (stage in stages) {
  message(sprintf("\n========== Processing %s ==========", stage))

  stage_samples <- rownames(meta_clean)[which(meta_clean$LifeStage == stage)]
  if (length(stage_samples) < MIN_SAMPLES_PER_STAGE) {
    message(sprintf("[%s] Skipping: Only %d samples available (< %d minimum required)",
                    stage, length(stage_samples), MIN_SAMPLES_PER_STAGE))
    next
  }

  sub_expr_clean <- expr_clean[, stage_samples, drop = FALSE]
  sub_meta <- meta_clean[stage_samples, , drop = FALSE]

  # inverse-normal transform per gene; constant genes map to 0 with names kept
  sub_expr_norm <- apply(as.matrix(sub_expr_clean), 1, function(x) {
    if (length(unique(x)) <= 1) {
      return(setNames(rep(0, length(x)), names(x)))
    }
    ranks <- rank(x, na.last = "keep", ties.method = "average")
    n <- length(x)
    qnorm(ranks / (n + 1))
  })
  sub_expr_norm <- t(sub_expr_norm)
  sub_expr_norm <- sub_expr_norm[complete.cases(sub_expr_norm), , drop = FALSE]
  if (is.null(colnames(sub_expr_norm))) {
    stop(sprintf("[%s] Error: inverse normal transform lost the sample names", stage))
  }
  out_filename <- sprintf("%spre.phenotype.combat.int.%s.txt", QP, stage)
  fwrite(as.data.frame(sub_expr_norm), file.path(out_dir, out_filename), sep = "\t", quote = FALSE, row.names = TRUE)
  message(sprintf("[%s] Saved: %s", stage, out_filename))

  # expression PCs: runBE chooses K
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
  message(sprintf("[%s] runBE chose %d expression PCs; using %d (min %d, max %d)",
                  stage, K_BE, K_final, K_PC_MIN, K_PC_MAX))

  # known covariates: every column coerced to numeric, batch/sample_id/LifeStage removed
  # (a multi-level categorical would become ordinal here; sex is binary, so this holds)
  remove_col <- c("batch", "sample_id", "LifeStage")
  sub_meta_numeric <- as.data.frame(lapply(sub_meta, function(x) {
    if (is.factor(x)) {
      as.numeric(x)
    } else if (is.character(x)) {
      as.numeric(as.factor(x))
    } else {
      as.numeric(x)
    }
  }))
  rownames(sub_meta_numeric) <- rownames(sub_meta)
  known_covs <- sub_meta_numeric[, !colnames(sub_meta_numeric) %in% remove_col, drop = FALSE]
  rownames(known_covs) <- rownames(sub_meta)

  # stage genotype: plink2 keep, LD-pruned PCA -> genoPC1..n
  pgen_prefix <- file.path(out_dir, paste0(QP, "genotype.", tolower(stage)))
  keep_file <- paste0(pgen_prefix, ".keep")
  keep_df <- data.frame("#IID" = rownames(known_covs), stringsAsFactors = FALSE, check.names = FALSE)
  fwrite(keep_df, keep_file, sep = "\t", quote = FALSE)

  pca_n <- min(GENO_N_PCS, nrow(known_covs) - 1)
  geno_sub <- NULL
  if (pca_n >= 1) {
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

    eigenvec_file <- paste0(pca_prefix, ".eigenvec")
    if (!file.exists(eigenvec_file)) {
      stop(sprintf("plink2 did not create eigenvec file: %s", eigenvec_file))
    }
    eigenvec <- fread(eigenvec_file, data.table = FALSE)
    iid_col <- intersect(c("#IID", "IID"), colnames(eigenvec))[1]
    if (is.na(iid_col)) {
      iid_col <- colnames(eigenvec)[min(2, ncol(eigenvec))]
    }
    rownames(eigenvec) <- as.character(eigenvec[[iid_col]])
    pc_cols <- grep("^PC[0-9]+$", colnames(eigenvec), value = TRUE)
    if (length(pc_cols) == 0) {
      pc_cols <- tail(colnames(eigenvec), pca_n)
    }
    geno_sub <- eigenvec[, pc_cols, drop = FALSE]
    colnames(geno_sub) <- paste0("genoPC", seq_len(ncol(geno_sub)))
  }

  # chromosome coding of the stage pvar, which tensorQTL is pointed at
  pvar_file <- paste0(pgen_prefix, ".pvar")
  stage_geno_chroms <- character()
  if (file.exists(pvar_file)) {
    pvar_lines <- readLines(pvar_file, n = 20000, warn = FALSE)
    pvar_lines <- pvar_lines[!startsWith(pvar_lines, "##") & !startsWith(pvar_lines, "#CHROM")]
    stage_geno_chroms <- unique(sub("\t.*$", "", pvar_lines))
  }
  stage_geno_has_chr_prefix <- length(stage_geno_chroms) > 0 && all(startsWith(stage_geno_chroms, "chr"))
  if (length(stage_geno_chroms) > 0) {
    message(sprintf("[%s] Stage genotype chromosome naming: %s (e.g. %s)", stage,
                    if (stage_geno_has_chr_prefix) "chr-prefixed" else "unprefixed",
                    paste(head(stage_geno_chroms, 3), collapse = ", ")))
  }

  if (!is.null(geno_sub)) {
    common_cov_samples <- intersect(rownames(known_covs), rownames(geno_sub))
    if (length(common_cov_samples) < MIN_SAMPLES_PER_STAGE) {
      stop(sprintf("[%s] Error: only %d samples remain after joining genotype PCs",
                   stage, length(common_cov_samples)))
    }
    known_covs <- cbind(
      known_covs[common_cov_samples, , drop = FALSE],
      geno_sub[common_cov_samples, , drop = FALSE]
    )
    sub_expr_norm <- sub_expr_norm[, common_cov_samples, drop = FALSE]
    sub_pca_final$x <- sub_pca_final$x[common_cov_samples, , drop = FALSE]
  }

  # constant covariates make the tensorQTL design rank deficient
  cov_var <- apply(known_covs, 2, var, na.rm = TRUE)
  const_covs <- names(which(!is.finite(cov_var) | cov_var <= COV_MIN_VARIANCE))
  if (length(const_covs) > 0) {
    log_message(sprintf("[%s] Dropping %d constant covariate(s): %s",
                        stage, length(const_covs), paste(const_covs, collapse = ", ")), "WARNING")
    known_covs <- known_covs[, setdiff(colnames(known_covs), const_covs), drop = FALSE]
  }
  if (anyNA(known_covs)) {
    stop(sprintf("[%s] Error: %d missing values in the known covariates", stage, sum(is.na(known_covs))))
  }

  PCsTop <- as.data.frame(sub_pca_final$x)[rownames(known_covs), 1:K_final, drop = FALSE]
  covsFiltered <- PCAForQTL::filterKnownCovariates(known_covs, PCsTop,
                                                   unadjustedR2_cutoff = KNOWN_COV_R2_CUTOFF)
  covsToUse <- cbind(covsFiltered, PCsTop)
  message(sprintf("[%s] Covariates: %d known kept (of %d) + %d expression PCs = %d total",
                  stage, ncol(covsFiltered), ncol(known_covs), ncol(PCsTop), ncol(covsToUse)))

  cov_out <- sprintf("%scovariates_for_qtl.%s.txt", QP, stage)
  fwrite(as.data.frame(t(covsToUse)), file.path(out_dir, cov_out), sep = "\t",
         row.names = TRUE, col.names = TRUE, quote = FALSE)
  message(sprintf("[%s] Saved: %s", stage, cov_out))

  # phenotype BED: annotation coordinates, chromosome naming matched to the stage genotype
  gene_ids <- rownames(sub_expr_norm)
  gene_ids_base <- sub("\\.[0-9]+$", "", gene_ids)
  gene_coords <- annotation_df[match(gene_ids_base, annotation_df$gene_id_base),
                               c("chr", "start", "end", "strand", "gene_id", "gene_id_base")]

  missing_coords <- is.na(gene_coords$chr)
  if (sum(missing_coords) > 0) {
    warning(sprintf("[%s] %d genes without coordinates found, removing them", stage, sum(missing_coords)))
    gene_ids <- gene_ids[!missing_coords]
    gene_coords <- gene_coords[!missing_coords, ]
    sub_expr_norm <- sub_expr_norm[gene_ids, , drop = FALSE]
  }

  final_samples <- rownames(covsToUse)
  if (!all(final_samples %in% colnames(sub_expr_norm))) {
    stop("Critical Error: Mismatch between covariate samples and expression samples!")
  }

  # BED intervals: tss = 1 bp window at the strand-aware TSS (GTEx);
  # gene_body = geneStart-1 to geneEnd. GTF coordinates are 1-based inclusive.
  if (PHENO_COORD_MODE == "gene_body") {
    bed_pos <- data.frame(start = gene_coords$start - 1L, end = gene_coords$end)
  } else {
    unknown_strand <- !gene_coords$strand %in% c("+", "-")
    if (any(unknown_strand)) {
      warning(sprintf("%d genes with unusable strand; TSS taken as the start coordinate",
                      sum(unknown_strand)))
    }
    is_minus <- !is.na(gene_coords$strand) & gene_coords$strand == "-"
    tss <- ifelse(is_minus, gene_coords$end, gene_coords$start)
    bed_pos <- data.frame(start = as.integer(tss) - 1L, end = as.integer(tss))
  }

  format_data <- data.frame(
    `#Chr` = gene_coords$chr,
    start = bed_pos$start,
    end = bed_pos$end,
    pid = gene_ids,
    sub_expr_norm[, final_samples, drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # match phenotype chromosome naming to the stage genotype before sorting; otherwise
  # tensorQTL fails its sorted-positions check
  pheno_has_chr_prefix <- all(startsWith(as.character(format_data$`#Chr`), "chr"))
  target_has_chr_prefix <- if (length(stage_geno_chroms) > 0) stage_geno_has_chr_prefix else pheno_has_chr_prefix
  if (target_has_chr_prefix && !pheno_has_chr_prefix) {
    log_message(sprintf("[%s] Adding 'chr' prefix to phenotype chromosomes to match the stage genotype", stage), "WARNING")
    format_data$`#Chr` <- paste0("chr", format_data$`#Chr`)
  } else if (!target_has_chr_prefix && pheno_has_chr_prefix) {
    log_message(sprintf("[%s] Stripping 'chr' prefix from phenotype chromosomes to match the stage genotype", stage), "WARNING")
    format_data$`#Chr` <- sub("^chr", "", format_data$`#Chr`)
  }

  # natural chromosome order, then position (tensorQTL requires sorted positions)
  chrom_core <- c(1:22, "X", "Y", "M")
  chrom_order <- if (target_has_chr_prefix) paste0("chr", chrom_core) else as.character(chrom_core)
  format_data$sort_chr <- factor(format_data$`#Chr`, levels = chrom_order, ordered = TRUE)
  if (anyNA(format_data$sort_chr)) {
    bad <- unique(as.character(format_data$`#Chr`)[is.na(format_data$sort_chr)])
    stop(sprintf("[%s] Error: %d genes on unrecognised chromosomes (%s); sorting would be wrong",
                 stage, sum(is.na(format_data$sort_chr)), paste(head(bad, 10), collapse = ", ")))
  }
  format_data <- format_data[order(format_data$sort_chr, format_data$start, format_data$end), ]
  format_data$sort_chr <- NULL
  rownames(format_data) <- NULL

  # header kept: tensorQTL expects the first column to be '#Chr'
  bed_file <- file.path(out_dir, sprintf("%sphenotype.%s.bed", QP, stage))
  fwrite(format_data, bed_file, sep = "\t", quote = FALSE, col.names = TRUE)
  run_or_stop(paste("bgzip -f", shQuote(bed_file)), sprintf("bgzip for %s", bed_file))
  run_or_stop(paste("tabix -p bed -f", shQuote(paste0(bed_file, ".gz"))),
              sprintf("tabix for %s.gz", bed_file))
  message(sprintf("[%s] Created: %s(.gz + .tbi)", stage, bed_file))

  sample_file <- file.path(out_dir, sprintf("%ssample.%s.final", QP, stage))
  writeLines(final_samples, sample_file)
  message(sprintf("[%s] Created: %s", stage, sample_file))

  stage_summary[[as.character(stage)]] <- data.frame(
    stage = as.character(stage),
    n_samples = length(final_samples),
    n_genes = nrow(format_data),
    n_known_covs = ncol(covsFiltered),
    n_expr_pcs = ncol(PCsTop),
    runBE_K = K_BE,
    coord_mode = PHENO_COORD_MODE,
    stringsAsFactors = FALSE
  )
  message(sprintf("[%s] Processing complete.", stage))
}

if (length(stage_summary) == 0) {
  stop("Error: no stage produced QTL inputs")
}
summary_tbl <- do.call(rbind, stage_summary)
fwrite(summary_tbl, file.path(out_dir, paste0(QP, "preprocess_summary.tsv")), sep = "\t")
message("\n========== Summary ==========")
print(summary_tbl)
message(sprintf("Seed: %d | Outputs: %s", RNG_SEED, out_dir))
