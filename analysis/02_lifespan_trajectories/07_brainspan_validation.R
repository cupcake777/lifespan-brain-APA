# BrainSpan and BrainGVEX fixed-template external validation of the lifespan trajectories.
suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(splines)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
RAW <- file.path(root, "raw")
REF <- file.path(root, "ref")

SITES_FILE <- file.path(DATA, "02_sites_clustered.tsv")
SUPPORT_FILE <- file.path(DATA, "01_support.tsv")
TRAJ_LOGIT_FILE <- file.path(DATA, "traj_logit.csv")
DISCOVERY_MATRIX_FILE <- file.path(DATA, "apa_after_combat_donor_mean.tsv")
DISCOVERY_COV_FILE <- file.path(DATA, "00_cov.tsv")
DISCOVERY_PROXIMAL_FILE <- file.path(root, "data/04_qtl/proximal.apa.info")

BRAINSPAN_PROXIMAL_FILE <- file.path(RAW, "brainspan.proximal.apa.info")
BRAINSPAN_META_FILE <- file.path(RAW, "all.apa.brainspan.meta")
BRAINSPAN_LOGIT_FILE <- file.path(RAW, "all.apa.brainspan.count")
BRAINSPAN_RAW_FILE <- file.path(RAW, "all.apa.brainspan.count")
BRAINSPAN_LIFTOVER_BIN <- Sys.getenv(
  "LIFTOVER_BIN", paste0(Sys.getenv("APA_HPC_HOME"), "/bin/liftOver")
)
BRAINSPAN_CHAIN_FILE <- Sys.getenv(
  "BRAINSPAN_CHAIN_FILE",
  file.path(REF, "liftover/hg19ToHg38.over.chain.gz")
)
GENCODE_VERSION <- "GENCODE v40"

ANNO_DIR <- Sys.getenv("ANNO_DIR", file.path(REF, "gencode"))
HG38_TRANSCRIPT_BED <- file.path(ANNO_DIR, "gencode.v40.ref.bed")
HG19_TRANSCRIPT_BED <- file.path(ANNO_DIR, "gencode.v40_hg19.ref.bed")

BRAINGVEX_MATRIX_FILE <- file.path(RAW, "after.combat.txt")
BRAINGVEX_SAMPLE_FILE <- file.path(RAW, "sample.Adult.final")
BRAINGVEX_META_FILE <- file.path(RAW, "eur.adult.meta")
BRAINGVEX_MATRIX_SCALE <- "ComBat-corrected log2-PDUI-logit"

OUT_MAPPING <- file.path(DATA, "08_external_mapping_qc.tsv")
OUT_IDENTIFIABILITY <- file.path(DATA, "08_template_identifiability.tsv")
OUT_BRAINSPAN <- file.path(DATA, "08_brainspan_event_replication.tsv")
OUT_BRAINGVEX <- file.path(DATA, "08_braingvex_event_replication.tsv")
OUT_CROSSFIT <- file.path(DATA, "08_discovery_crossfit.tsv")
OUT_SENSITIVITY <- file.path(DATA, "08_external_replication_sensitivity.tsv")
OUT_OVERLAY <- file.path(DATA, "08_external_trajectory_overlay.tsv")
OUT_LIFTOVER_QC <- file.path(DATA, "08_brainspan_liftover_qc.tsv")

# --- settings ---

REFERENCE_PAS_TOLERANCE <- 100
PAS_SENSITIVITY <- c(50, 100, 250, 500, Inf)
PDUI_EPSILON <- 0.001
MIN_SEGMENT_SAMPLES <- 5L
MIN_TEMPLATE_RANGE <- 0.05
MIN_TEMPLATE_SD <- 0.02
MIN_TEMPLATE_RESIDUAL_SS <- 0.02
MAX_TEMPLATE_VIF <- 10
MIN_BATCH_N <- 20L
MIN_BATCH_AGE_SPAN <- 10
MIN_META_BATCHES <- 2L
ADULT_MIN_AGE <- 18
ADULT_MAX_AGE <- 90
BRAINSPAN_SMOOTH_DF <- 2L
BRAINGVEX_SMOOTH_DF <- 3L
MIN_SMOOTH_BOOT_SUCCESS <- 0.90
MAX_SHAPE_CI_WIDTH <- 1.50
SMOOTH_BOOT_ITER <- as.integer(Sys.getenv("APA_EXTERNAL_SMOOTH_BOOT", "200"))
CROSSFIT_FOLDS <- 5L
CROSSFIT_BASIS_K <- 20L
CROSSFIT_Q <- 0.05
CROSSFIT_WORKERS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
CROSSFIT_MAX_SITES <- as.integer(Sys.getenv("APA_EXTERNAL_CROSSFIT_MAX_SITES", "0"))
RUN_DISCOVERY_CROSSFIT <- Sys.getenv("APA_EXTERNAL_RUN_CROSSFIT", "1") == "1"
SEED <- 2026L

set.seed(SEED)

# --- helper functions ---

strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))

parse_field <- function(ids, pos) {
  vapply(strsplit(as.character(ids), "|", fixed = TRUE), function(x) {
    if (length(x) >= pos) x[pos] else ""
  }, character(1))
}

parse_locus <- function(loci) {
  loci <- as.character(loci)
  chr <- sub(":.*$", "", loci)
  coords <- sub("^[^:]+:", "", loci)
  start <- as.integer(sub("-.*$", "", coords))
  end <- as.integer(sub("^.*-", "", coords))
  data.table(
    chr = chr,
    start = pmin(start, end),
    end = pmax(start, end)
  )
}

run_liftover_bed <- function(coords, label) {
  if (!file.exists(BRAINSPAN_LIFTOVER_BIN)) {
    stop("Missing liftOver binary: ", BRAINSPAN_LIFTOVER_BIN)
  }
  if (!file.exists(BRAINSPAN_CHAIN_FILE)) {
    stop("Missing hg19-to-hg38 chain: ", BRAINSPAN_CHAIN_FILE)
  }
  prefix <- tempfile(paste0("brainspan_", label, "_"), tmpdir = tempdir())
  input_file <- paste0(prefix, ".bed")
  mapped_file <- paste0(prefix, ".mapped.bed")
  unmapped_file <- paste0(prefix, ".unmapped.bed")
  log_file <- paste0(prefix, ".log")
  on.exit(unlink(c(input_file, mapped_file, unmapped_file, log_file)), add = TRUE)

  valid <- coords[
    !is.na(chr) & nzchar(chr) & is.finite(start) & is.finite(end) &
      start >= 1L & end >= start
  ]
  bed <- valid[, .(
    chr, bed_start = start - 1L, bed_end = end,
    event_index
  )]
  fwrite(bed, input_file, sep = "\t", col.names = FALSE, quote = FALSE)
  status <- system2(
    BRAINSPAN_LIFTOVER_BIN,
    c(input_file, BRAINSPAN_CHAIN_FILE, mapped_file, unmapped_file),
    stdout = log_file, stderr = log_file
  )
  if (!identical(as.integer(status), 0L)) {
    detail <- if (file.exists(log_file)) {
      paste(tail(readLines(log_file, warn = FALSE), 20L), collapse = " | ")
    } else ""
    stop("liftOver failed: ", detail)
  }

  mapped <- if (file.exists(mapped_file) && file.info(mapped_file)$size > 0) {
    fread(mapped_file, header = FALSE, col.names = c(
      "lift_chr", "lift_start0", "lift_end", "event_index"
    ), showProgress = FALSE)
  } else data.table(
    lift_chr = character(), lift_start0 = integer(), lift_end = integer(),
    event_index = integer()
  )
  if (nrow(mapped)) {
    mapped[, `:=`(
      lift_start = as.integer(lift_start0) + 1L,
      lift_end = as.integer(lift_end),
      lift_width = as.integer(lift_end) - as.integer(lift_start0)
    )]
    setorder(mapped, event_index, -lift_width, lift_chr, lift_start)
    mapped <- mapped[, .SD[1L], by = event_index]
  }
  list(
    mapped = mapped,
    valid_input_n = nrow(valid),
    mapped_n = uniqueN(mapped$event_index),
    unmapped_n = nrow(valid) - uniqueN(mapped$event_index)
  )
}

read_transcript_bed <- function(path, prefix) {
  dt <- fread(
    path, header = FALSE,
    select = c(1L, 2L, 3L, 4L, 6L, 10L, 11L, 12L),
    col.names = c("chr", "start", "end", "transcript_id", "strand",
                  "block_n", "block_sizes", "block_starts")
  )
  setnames(dt, setdiff(names(dt), "transcript_id"),
           paste0(prefix, "_", setdiff(names(dt), "transcript_id")))
  dt[!duplicated(transcript_id)]
}

read_transcript_interval_bed <- function(path, prefix) {
  dt <- fread(
    path, header = FALSE,
    select = c(1L, 2L, 3L, 4L, 6L, 10L, 11L, 12L),
    col.names = c("bed_chr", "bed_start", "bed_end", "transcript_id",
                  "bed_strand", "block_n", "block_sizes", "block_starts")
  )
  setnames(dt, names(dt), paste0(prefix, "_", names(dt)))
  dt[!duplicated(get(paste0(prefix, "_transcript_id")))]
}

transcript_offset <- function(coord, transcript_start, sizes_text,
                              starts_text, strand) {
  sizes <- scan(text = sizes_text, sep = ",", quiet = TRUE)
  starts <- scan(text = starts_text, sep = ",", quiet = TRUE)
  if (!length(sizes) || length(sizes) != length(starts)) return(NA_real_)
  genomic_starts <- transcript_start + starts
  block <- which(
    coord - 1L >= genomic_starts & coord - 1L < genomic_starts + sizes
  )
  if (length(block) != 1L) return(NA_real_)
  offset_plus <- sum(sizes[seq_len(block - 1L)]) +
    (coord - 1L - genomic_starts[block])
  if (strand == "+") offset_plus else sum(sizes) - 1L - offset_plus
}

add_transcript_interval <- function(dt, prefix) {
  event_start <- paste0(prefix, "_event_start")
  event_end <- paste0(prefix, "_event_end")
  bed_start <- paste0(prefix, "_bed_start")
  block_sizes <- paste0(prefix, "_block_sizes")
  block_starts <- paste0(prefix, "_block_starts")
  event_strand <- paste0(prefix, "_event_strand")
  offset_start <- paste0(prefix, "_offset_start")
  offset_end <- paste0(prefix, "_offset_end")
  interval_low <- paste0(prefix, "_interval_low")
  interval_high <- paste0(prefix, "_interval_high")
  dt[, (offset_start) := mapply(
    transcript_offset, get(event_start), get(bed_start), get(block_sizes),
    get(block_starts), get(event_strand)
  )]
  dt[, (offset_end) := mapply(
    transcript_offset, get(event_end), get(bed_start), get(block_sizes),
    get(block_starts), get(event_strand)
  )]
  dt[, (interval_low) := pmin(get(offset_start), get(offset_end))]
  dt[, (interval_high) := pmax(get(offset_start), get(offset_end))]
  dt
}

map_bed12_coordinate <- function(coord, source_start, source_sizes,
                                 source_starts, target_start, target_sizes,
                                 target_starts) {
  source_size <- scan(text = source_sizes, sep = ",", quiet = TRUE)
  source_rel <- scan(text = source_starts, sep = ",", quiet = TRUE)
  target_size <- scan(text = target_sizes, sep = ",", quiet = TRUE)
  target_rel <- scan(text = target_starts, sep = ",", quiet = TRUE)
  if (length(source_size) != length(target_size)) return(NA_real_)
  pos0 <- coord - 1L
  block_start <- source_start + source_rel
  block <- which(pos0 >= block_start & pos0 < block_start + source_size)
  if (length(block) != 1L) return(NA_real_)
  fraction <- (pos0 - block_start[block]) / source_size[block]
  target_start + target_rel[block] +
    round(fraction * target_size[block]) + 1L
}

full_rank <- function(x) {
  x[, qr(x)$pivot[seq_len(qr(x)$rank)], drop = FALSE]
}

braingvex_nuisance <- function(dt, include_diagnosis = TRUE) {
  columns <- list(intercept = rep(1, nrow(dt)))
  if (uniqueN(dt$Sex) > 1L) columns$sex <- as.integer(factor(dt$Sex)) - 1L
  if (sd(dt$RIN) > 0) columns$RIN_z <- as.numeric(scale(dt$RIN))
  if (sd(dt$PMI) > 0) columns$PMI_z <- as.numeric(scale(dt$PMI))
  if (include_diagnosis && uniqueN(dt$Diagnosis) > 1L) {
    diagnosis <- model.matrix(~ factor(Diagnosis), data = dt)[, -1L, drop = FALSE]
    for (j in seq_len(ncol(diagnosis))) {
      columns[[paste0("diagnosis_", j)]] <- diagnosis[, j]
    }
  }
  full_rank(do.call(cbind, columns))
}

interpolation_weights <- function(grid_age, target_age) {
  out <- matrix(0, nrow = length(grid_age), ncol = length(target_age))
  for (j in seq_along(target_age)) {
    age <- target_age[j]
    upper <- which(grid_age >= age)[1]
    lower <- tail(which(grid_age <= age), 1)
    if (!length(lower) || is.na(lower)) lower <- 1L
    if (!length(upper) || is.na(upper)) upper <- length(grid_age)
    if (lower == upper) {
      out[lower, j] <- 1
    } else {
      weight_upper <- (age - grid_age[lower]) /
        (grid_age[upper] - grid_age[lower])
      out[lower, j] <- 1 - weight_upper
      out[upper, j] <- weight_upper
    }
  }
  out
}

fit_template_rows <- function(y, template, nuisance) {
  nuisance <- full_rank(nuisance)
  nuisance_inv <- solve(crossprod(nuisance))
  residualize <- function(x) {
    x - (x %*% nuisance %*% nuisance_inv) %*% t(nuisance)
  }
  y_res <- residualize(y)
  t_res <- residualize(template)
  template_ss <- rowSums(t_res^2)
  beta <- rowSums(t_res * y_res) / template_ss
  residual <- y_res - beta * t_res
  h_nuisance <- rowSums((nuisance %*% nuisance_inv) * nuisance)
  leverage <- sweep(t_res^2, 1L, template_ss, "/") +
    matrix(h_nuisance, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
  hc3_term <- t_res * residual / pmax(1 - leverage, 1e-8)
  se <- sqrt(rowSums(hc3_term^2) / template_ss^2)
  df <- ncol(y) - ncol(nuisance) - 1L
  pvalue <- 2 * pt(abs(beta / se), df = df, lower.tail = FALSE)
  y_ss <- rowSums(y_res^2)
  partial_r <- rowSums(t_res * y_res) / sqrt(template_ss * y_ss)
  partial_r2 <- pmax(0, pmin(1, partial_r^2))
  centered_ss <- rowSums((template - rowMeans(template))^2)
  data.table(
    beta = beta,
    se_hc3 = se,
    ci_low = beta - qt(0.975, df) * se,
    ci_high = beta + qt(0.975, df) * se,
    p_hc3 = pvalue,
    partial_r = partial_r,
    partial_r2 = partial_r2,
    residual_df = df,
    template_residual_ss = template_ss,
    template_vif = centered_ss / template_ss,
    residual_sd = sqrt(rowSums(residual^2) / df),
    max_leverage = apply(leverage, 1L, max)
  )
}

fit_smooth_rows <- function(y, age, sex, prediction_age, df) {
  basis <- ns(age, df = df)
  design <- full_rank(cbind(intercept = 1, sex = sex, basis))
  prediction_basis <- predict(basis, prediction_age)
  prediction <- cbind(intercept = 1, sex = mean(sex), prediction_basis)
  prediction <- prediction[, colnames(design), drop = FALSE]
  coefficients <- y %*% design %*% solve(crossprod(design))
  coefficients %*% t(prediction)
}

row_shape_metrics <- function(template, external) {
  template_z <- t(scale(t(template)))
  external_z <- t(scale(t(external)))
  pearson <- rowMeans(template_z * external_z)
  spearman <- vapply(seq_len(nrow(template)), function(i) {
    suppressWarnings(cor(template[i, ], external[i, ], method = "spearman"))
  }, numeric(1))
  affine_rmse <- vapply(seq_len(nrow(template)), function(i) {
    fit <- lm(external[i, ] ~ template[i, ])
    sqrt(mean(residuals(fit)^2)) / sd(external[i, ])
  }, numeric(1))
  data.table(shape_pearson = pearson, shape_spearman = spearman,
             normalized_rmse = affine_rmse)
}

bootstrap_shape_metrics <- function(y, age, sex, prediction_age,
                                    template_grid, df, iterations) {
  n_event <- nrow(y)
  pearson <- matrix(NA_real_, n_event, iterations)
  spearman <- matrix(NA_real_, n_event, iterations)
  rmse <- matrix(NA_real_, n_event, iterations)
  sex_groups <- split(seq_along(sex), sex)
  for (b in seq_len(iterations)) {
    index <- unlist(lapply(sex_groups, sample, replace = TRUE), use.names = FALSE)
    smooth <- tryCatch(
      fit_smooth_rows(y[, index, drop = FALSE], age[index], sex[index],
                      prediction_age, df),
      error = function(e) NULL
    )
    if (!is.null(smooth)) {
      metric <- row_shape_metrics(template_grid, smooth)
      pearson[, b] <- metric$shape_pearson
      spearman[, b] <- metric$shape_spearman
      rmse[, b] <- metric$normalized_rmse
    }
  }
  data.table(
    shape_pearson_ci_low = apply(pearson, 1L, quantile, 0.025, na.rm = TRUE),
    shape_pearson_ci_high = apply(pearson, 1L, quantile, 0.975, na.rm = TRUE),
    shape_spearman_ci_low = apply(spearman, 1L, quantile, 0.025, na.rm = TRUE),
    shape_spearman_ci_high = apply(spearman, 1L, quantile, 0.975, na.rm = TRUE),
    normalized_rmse_ci_low = apply(rmse, 1L, quantile, 0.025, na.rm = TRUE),
    normalized_rmse_ci_high = apply(rmse, 1L, quantile, 0.975, na.rm = TRUE),
    smooth_boot_success = rowSums(is.finite(pearson))
  )
}

assign_mismatched_template <- function(dt, variables) {
  result <- integer(nrow(dt))
  for (cluster in sort(unique(dt$cluster_id))) {
    index <- which(dt$cluster_id == cluster)
    x <- scale(as.matrix(dt[index, ..variables]))
    x[!is.finite(x)] <- 0
    informative <- apply(x, 2L, sd) > 0
    x <- x[, informative, drop = FALSE]
    score <- if (ncol(x) <= 1L) {
      if (ncol(x) == 1L) x[, 1L] else seq_along(index)
    } else {
      prcomp(x, center = FALSE, scale. = FALSE, rank. = 1L)$x[, 1L]
    }
    n_groups <- min(100L, max(1L, floor(length(index) / 30L)))
    groups <- ceiling(rank(score, ties.method = "first") /
                      length(index) * n_groups)
    for (group in unique(groups)) {
      local <- which(groups == group)
      if (length(local) == 1L) {
        candidate <- setdiff(seq_along(index), local)
        distance <- rowSums((x[candidate, , drop = FALSE] -
          matrix(x[local, ], nrow = length(candidate), ncol = ncol(x),
                 byrow = TRUE))^2)
        result[index[local]] <- index[candidate[which.min(distance)]]
      } else {
        distance <- as.matrix(dist(x[local, , drop = FALSE]))
        diag(distance) <- Inf
        result[index[local]] <- index[local[max.col(-distance)]]
      }
    }
  }
  result
}

row_min_finite <- function(x) {
  apply(x, 1L, function(y) if (any(is.finite(y))) min(y, na.rm = TRUE) else NA_real_)
}

row_max_finite <- function(x) {
  apply(x, 1L, function(y) if (any(is.finite(y))) max(y, na.rm = TRUE) else NA_real_)
}

fit_gvex_smooth_rows <- function(y, meta, prediction_age, df) {
  factors <- model.matrix(
    ~ factor(batch) + factor(Sex) + factor(Diagnosis), data = meta
  )
  basis <- ns(meta$AgeDeath, df = df)
  design <- full_rank(cbind(factors, RIN_z = as.numeric(scale(meta$RIN)),
                            PMI_z = as.numeric(scale(meta$PMI)), basis))
  prediction <- matrix(0, nrow = length(prediction_age), ncol = ncol(design),
                       dimnames = list(NULL, colnames(design)))
  prediction[, "(Intercept)"] <- 1
  basis_columns <- intersect(colnames(design), colnames(basis))
  prediction[, basis_columns] <- predict(basis, prediction_age)[, basis_columns,
                                                                  drop = FALSE]
  coefficients <- y %*% design %*% solve(crossprod(design))
  list(prediction = coefficients %*% t(prediction),
       design = design, prediction_design = prediction)
}

bootstrap_gvex_shape <- function(y, meta, design, prediction_design,
                                 template_grid, iterations) {
  pearson <- matrix(NA_real_, nrow(y), iterations)
  spearman <- matrix(NA_real_, nrow(y), iterations)
  rmse <- matrix(NA_real_, nrow(y), iterations)
  batch_groups <- split(seq_len(nrow(meta)), meta$batch)
  for (b in seq_len(iterations)) {
    index <- unlist(lapply(batch_groups, sample, replace = TRUE),
                    use.names = FALSE)
    design_b <- full_rank(design[index, , drop = FALSE])
    prediction_b <- prediction_design[, colnames(design_b), drop = FALSE]
    smooth <- tryCatch({
      coefficient <- y[, index, drop = FALSE] %*% design_b %*%
        solve(crossprod(design_b))
      coefficient %*% t(prediction_b)
    }, error = function(e) NULL)
    if (!is.null(smooth)) {
      metric <- row_shape_metrics(template_grid, smooth)
      pearson[, b] <- metric$shape_pearson
      spearman[, b] <- metric$shape_spearman
      rmse[, b] <- metric$normalized_rmse
    }
  }
  data.table(
    shape_pearson_ci_low = apply(pearson, 1L, quantile, 0.025, na.rm = TRUE),
    shape_pearson_ci_high = apply(pearson, 1L, quantile, 0.975, na.rm = TRUE),
    shape_spearman_ci_low = apply(spearman, 1L, quantile, 0.025, na.rm = TRUE),
    shape_spearman_ci_high = apply(spearman, 1L, quantile, 0.975, na.rm = TRUE),
    normalized_rmse_ci_low = apply(rmse, 1L, quantile, 0.025, na.rm = TRUE),
    normalized_rmse_ci_high = apply(rmse, 1L, quantile, 0.975, na.rm = TRUE),
    smooth_boot_success = rowSums(is.finite(pearson))
  )
}

meta_analyze <- function(beta, se, estimable) {
  weight <- ifelse(estimable, 1 / se^2, 0)
  weight[!is.finite(weight)] <- 0
  n_batch <- rowSums(weight > 0)
  fixed <- rowSums(weight * beta, na.rm = TRUE) / rowSums(weight)
  fixed_se <- sqrt(1 / rowSums(weight))
  q <- rowSums(weight * (beta - fixed)^2, na.rm = TRUE)
  q_df <- pmax(n_batch - 1L, 0L)
  c_value <- rowSums(weight) - rowSums(weight^2) / rowSums(weight)
  tau2 <- pmax(0, (q - q_df) / c_value)
  random_weight <- ifelse(estimable, 1 / (se^2 + tau2), 0)
  random_weight[!is.finite(random_weight)] <- 0
  random <- rowSums(random_weight * beta, na.rm = TRUE) /
    rowSums(random_weight)
  random_se <- sqrt(1 / rowSums(random_weight))
  data.table(
    n_estimable_batches = n_batch,
    beta_meta = fixed,
    se_meta = fixed_se,
    ci_low = fixed - 1.96 * fixed_se,
    ci_high = fixed + 1.96 * fixed_se,
    p_meta = 2 * pnorm(abs(fixed / fixed_se), lower.tail = FALSE),
    cochran_q = q,
    cochran_q_p = pchisq(q, q_df, lower.tail = FALSE),
    i2 = fifelse(q > 0, pmax(0, (q - q_df) / q), 0),
    tau2 = tau2,
    beta_random = random,
    se_random = random_se,
    random_ci_low = random - 1.96 * random_se,
    random_ci_high = random + 1.96 * random_se
  )
}

# --- frozen discovery events and trajectories ---

sites <- fread(SITES_FILE)
dynamic <- sites[directional_pass == TRUE & !is.na(cluster_id), .(
  site_id, transcript_id, tx_base = strip_ver(transcript_id), gene_id,
  delta_logit, amplitude_logit, cluster_id,
  cluster_label = paste0("C", cluster_id), delta_logit_threshold
)]
support <- fread(SUPPORT_FILE)
supported <- copy(support[is_supported == TRUE])
supported[, segment_id := cumsum(c(TRUE, diff(grid_index) != 1L))]
trajectory_cols <- sprintf("grid_%03d", supported$grid_index)
trajectory <- fread(
  TRAJ_LOGIT_FILE, select = c("site_id", trajectory_cols), showProgress = FALSE
)
trajectory <- trajectory[match(dynamic$site_id, site_id)]
trajectory_mat <- as.matrix(trajectory[, ..trajectory_cols])
storage.mode(trajectory_mat) <- "numeric"
rownames(trajectory_mat) <- trajectory$site_id
dynamic <- dynamic[match(rownames(trajectory_mat), site_id)]

# --- BrainSpan gene + version-stripped ENST remapping ---

brainspan_meta <- fread(BRAINSPAN_META_FILE)
brainspan_meta[, age_birth_aligned := as.numeric(age)]
brainspan_meta[, segment_id := NA_integer_]
for (segment in unique(supported$segment_id)) {
  age_range <- range(supported[segment_id == segment, grid_age])
  brainspan_meta[
    age_birth_aligned >= age_range[1] & age_birth_aligned <= age_range[2],
    segment_id := segment
  ]
}
brainspan_meta <- brainspan_meta[!is.na(segment_id)]
brainspan_meta[, sex_binary := as.integer(sex == "Male")]

brainspan_raw <- fread(BRAINSPAN_RAW_FILE, showProgress = FALSE)
raw_annotation <- cbind(
  brainspan_raw[, .(
    event_index = .I, Gene, Predicted_Proximal_APA, Loci
  )],
  parse_locus(brainspan_raw$Loci)
)
raw_annotation[, `:=`(
  brainspan_transcript = parse_field(Gene, 1L),
  tx_base = strip_ver(parse_field(Gene, 1L)),
  gene_id = parse_field(Gene, 2L),
  brainspan_strand = parse_field(Gene, 4L),
  brainspan_site_id = paste(
    parse_field(Gene, 1L), parse_field(Gene, 2L), Loci,
    parse_field(Gene, 4L), sep = "|"
  )
)]
event_lift <- run_liftover_bed(
  raw_annotation[, .(event_index, chr, start, end)], "events"
)
raw_annotation <- merge(raw_annotation, event_lift$mapped,
                        by = "event_index", all.x = TRUE)
pas_lift <- run_liftover_bed(
  raw_annotation[, .(
    event_index, chr, start = as.integer(Predicted_Proximal_APA),
    end = as.integer(Predicted_Proximal_APA)
  )], "pas"
)
raw_annotation <- merge(
  raw_annotation,
  pas_lift$mapped[, .(event_index, lifted_pas_hg38 = lift_start)],
  by = "event_index", all.x = TRUE
)

brainspan_rows <- raw_annotation[, .(
  brainspan_site_id, brainspan_transcript, tx_base, gene_id,
  brainspan_strand, raw_chr = chr, raw_start = start, raw_end = end,
  lift_chr, lift_start, lift_end, lifted_pas_hg38
)]
mapping <- merge(dynamic, brainspan_rows, by = c("tx_base", "gene_id"),
                 allow.cartesian = TRUE)
mapping[, matching_key := paste(gene_id, tx_base, sep = "|")]
mapping[, exact_transcript_version := transcript_id == brainspan_transcript]
mapping[, `:=`(
  discovery_locus = parse_field(site_id, 3L),
  discovery_strand = parse_field(site_id, 4L),
  discovery_chr = sub(":.*$", "", parse_field(site_id, 3L))
)]
discovery_interval <- parse_locus(mapping$discovery_locus)
mapping[, `:=`(
  discovery_start = discovery_interval$start,
  discovery_end = discovery_interval$end,
  overlap_bp = pmax(
    0L, pmin(discovery_interval$end, lift_end) -
      pmax(discovery_interval$start, lift_start) + 1L
  ),
  midpoint_distance = abs(
    (discovery_interval$start + discovery_interval$end) / 2 -
      (lift_start + lift_end) / 2
  ),
  strand_pass = discovery_strand == brainspan_strand,
  chromosome_pass = discovery_chr == lift_chr
)]
mapping[, coordinate_match :=
  strand_pass & chromosome_pass & is.finite(overlap_bp) & overlap_bp > 0
]
mapping <- mapping[coordinate_match == TRUE]
setorder(mapping, site_id, -overlap_bp, midpoint_distance,
         -exact_transcript_version, brainspan_site_id)
mapping <- mapping[, .SD[1L], by = site_id]

# PAS distance remains a diagnostic sensitivity and is not the mapping gate.
discovery_proximal <- fread(
  DISCOVERY_PROXIMAL_FILE, select = c("Predicted_Proximal_APA", "phenotype_id")
)
setnames(discovery_proximal, c("discovery_pas_hg38", "site_id"))
mapping <- merge(mapping, discovery_proximal, by = "site_id", all.x = TRUE)
mapping[, `:=`(
  mapped_brainspan_pas_hg38 = lifted_pas_hg38,
  bridge_pass = coordinate_match,
  pas_distance_bp = abs(lifted_pas_hg38 - discovery_pas_hg38),
  mapped_qc = coordinate_match
)]

liftover_qc <- data.table(
  source_matrix = BRAINSPAN_RAW_FILE,
  source_coordinate_build = "hg19",
  target_coordinate_build = "GRCh38",
  gencode_version = GENCODE_VERSION,
  n_raw_events = nrow(raw_annotation),
  n_coordinate_valid_events = event_lift$valid_input_n,
  n_lifted_events = event_lift$mapped_n,
  n_liftover_lost_events = event_lift$unmapped_n,
  n_dynamic_sites = nrow(dynamic),
  n_strictly_matched_dynamic_sites = nrow(mapping),
  n_matched_raw_events = uniqueN(mapping$brainspan_site_id),
  n_lifted_unmatched_raw_events = event_lift$mapped_n -
    uniqueN(mapping$brainspan_site_id),
  match_definition = paste(
    "gene_id + version-stripped ENST + same strand + same chromosome +",
    "positive hg38 interval overlap; maximum overlap then midpoint distance"
  ),
  liftOver_binary = BRAINSPAN_LIFTOVER_BIN,
  liftOver_chain = BRAINSPAN_CHAIN_FILE
)

# --- BrainSpan fixed-template models ---

raw_sample_columns <- setdiff(
  names(brainspan_raw),
  c("Gene", "fit_value", "Predicted_Proximal_APA", "Loci")
)
brainspan_logit <- brainspan_raw[, c("Gene", "Loci", raw_sample_columns),
                                 with = FALSE]
brainspan_logit[, brainspan_site_id := paste(
  parse_field(Gene, 1L), parse_field(Gene, 2L), Loci,
  parse_field(Gene, 4L), sep = "|"
)]
brainspan_logit[, c("Gene", "Loci") := NULL]
setcolorder(brainspan_logit, c("brainspan_site_id", raw_sample_columns))
for (sample_name in raw_sample_columns) {
  pdui <- as.numeric(brainspan_logit[[sample_name]])
  pdui <- pmin(pmax(pdui, PDUI_EPSILON), 1 - PDUI_EPSILON)
  set(brainspan_logit, j = sample_name, value = log2(pdui / (1 - pdui)))
}

samples <- intersect(brainspan_meta$sample, names(brainspan_logit))
brainspan_meta <- brainspan_meta[match(samples, sample)]
brainspan_keep <- which(mapping$mapped_qc)
mapping_primary <- mapping[brainspan_keep]
trajectory_index <- match(mapping_primary$site_id, rownames(trajectory_mat))
brainspan_template <- matrix(NA_real_, nrow(mapping_primary), length(samples))
for (segment in unique(brainspan_meta$segment_id)) {
  sample_index <- which(brainspan_meta$segment_id == segment)
  grid_index <- which(supported$segment_id == segment)
  weight <- interpolation_weights(
    supported$grid_age[grid_index],
    brainspan_meta$age_birth_aligned[sample_index]
  )
  brainspan_template[, sample_index] <-
    trajectory_mat[trajectory_index, grid_index, drop = FALSE] %*% weight
}
logit_index <- match(mapping_primary$brainspan_site_id,
                     brainspan_logit$brainspan_site_id)
brainspan_y <- as.matrix(brainspan_logit[logit_index, ..samples])
storage.mode(brainspan_y) <- "numeric"
finite_count <- rowSums(is.finite(brainspan_y))
finite_prenatal <- rowSums(is.finite(
  brainspan_y[, brainspan_meta$segment_id == 1L, drop = FALSE]
))
finite_adult <- rowSums(is.finite(
  brainspan_y[, brainspan_meta$segment_id == 2L, drop = FALSE]
))
brainspan_y[!is.finite(brainspan_y)] <- rowMeans(brainspan_y, na.rm = TRUE)[
  row(brainspan_y)[!is.finite(brainspan_y)]
]

global_nuisance <- model.matrix(~ sex_binary, data = brainspan_meta)
global_fit <- fit_template_rows(brainspan_y, brainspan_template,
                                global_nuisance)
segment_fits <- list()
for (segment in sort(unique(brainspan_meta$segment_id))) {
  index <- which(brainspan_meta$segment_id == segment)
  nuisance <- model.matrix(~ sex_binary, data = brainspan_meta[index])
  segment_fits[[as.character(segment)]] <- fit_template_rows(
    brainspan_y[, index, drop = FALSE],
    brainspan_template[, index, drop = FALSE], nuisance
  )
}
prenatal_fit <- segment_fits[["1"]]
adult_fit <- segment_fits[["2"]]
heterogeneity_z <- (prenatal_fit$beta - adult_fit$beta) /
  sqrt(prenatal_fit$se_hc3^2 + adult_fit$se_hc3^2)

brainspan_events <- cbind(
  mapping_primary[, .(site_id, brainspan_site_id, transcript_id, gene_id,
                       cluster_id, cluster_label, delta_logit,
                       amplitude_logit, pas_distance_bp, mapped_qc)],
  global_fit
)
setnames(brainspan_events,
         c("beta", "se_hc3", "ci_low", "ci_high", "p_hc3", "partial_r",
           "partial_r2", "residual_df", "template_residual_ss", "template_vif",
           "residual_sd", "max_leverage"),
         paste0(c("beta", "se_hc3", "ci_low", "ci_high", "p_hc3", "partial_r",
                  "partial_r2", "residual_df", "template_residual_ss", "template_vif",
                  "residual_sd", "max_leverage"), "_global"))
brainspan_events[, `:=`(
  beta_prenatal = prenatal_fit$beta,
  se_prenatal = prenatal_fit$se_hc3,
  prenatal_ci_low = prenatal_fit$ci_low,
  prenatal_ci_high = prenatal_fit$ci_high,
  prenatal_p = prenatal_fit$p_hc3,
  beta_adult = adult_fit$beta,
  se_adult = adult_fit$se_hc3,
  adult_ci_low = adult_fit$ci_low,
  adult_ci_high = adult_fit$ci_high,
  adult_p = adult_fit$p_hc3,
  segment_heterogeneity_z = heterogeneity_z,
  segment_heterogeneity_p = 2 * pnorm(abs(heterogeneity_z), lower.tail = FALSE),
  n_finite = finite_count,
  n_finite_prenatal = finite_prenatal,
  n_finite_adult = finite_adult,
  template_range = apply(brainspan_template, 1L, function(x) diff(range(x))),
  template_sd = apply(brainspan_template, 1L, sd)
)]
brainspan_events[, estimable :=
  n_finite_prenatal >= MIN_SEGMENT_SAMPLES &
  n_finite_adult >= MIN_SEGMENT_SAMPLES &
  template_range >= MIN_TEMPLATE_RANGE & template_sd >= MIN_TEMPLATE_SD &
  template_residual_ss_global >= MIN_TEMPLATE_RESIDUAL_SS &
  template_vif_global <= MAX_TEMPLATE_VIF]

# --- BrainSpan smooth shape uncertainty ---

smooth_point <- matrix(NA_real_, nrow(brainspan_events), nrow(supported))
shape_rows <- list()
for (segment in sort(unique(brainspan_meta$segment_id))) {
  sample_index <- which(brainspan_meta$segment_id == segment)
  grid_index <- which(supported$segment_id == segment)
  external_smooth <- fit_smooth_rows(
    brainspan_y[, sample_index, drop = FALSE],
    brainspan_meta$age_birth_aligned[sample_index],
    brainspan_meta$sex_binary[sample_index],
    supported$grid_age[grid_index], BRAINSPAN_SMOOTH_DF
  )
  template_grid <- trajectory_mat[trajectory_index, grid_index, drop = FALSE]
  smooth_point[, grid_index] <- external_smooth
  metric <- row_shape_metrics(template_grid, external_smooth)
  uncertainty <- bootstrap_shape_metrics(
    brainspan_y[, sample_index, drop = FALSE],
    brainspan_meta$age_birth_aligned[sample_index],
    brainspan_meta$sex_binary[sample_index],
    supported$grid_age[grid_index], template_grid,
    BRAINSPAN_SMOOTH_DF, SMOOTH_BOOT_ITER
  )
  shape_rows[[as.character(segment)]] <- cbind(metric, uncertainty)
}
brainspan_events[, `:=`(
  shape_pearson_prenatal = shape_rows[["1"]]$shape_pearson,
  shape_pearson_prenatal_ci_low = shape_rows[["1"]]$shape_pearson_ci_low,
  shape_pearson_prenatal_ci_high = shape_rows[["1"]]$shape_pearson_ci_high,
  shape_spearman_prenatal = shape_rows[["1"]]$shape_spearman,
  shape_spearman_prenatal_ci_low = shape_rows[["1"]]$shape_spearman_ci_low,
  shape_spearman_prenatal_ci_high = shape_rows[["1"]]$shape_spearman_ci_high,
  nrmse_prenatal = shape_rows[["1"]]$normalized_rmse,
  nrmse_prenatal_ci_low = shape_rows[["1"]]$normalized_rmse_ci_low,
  nrmse_prenatal_ci_high = shape_rows[["1"]]$normalized_rmse_ci_high,
  shape_pearson_adult = shape_rows[["2"]]$shape_pearson,
  shape_pearson_adult_ci_low = shape_rows[["2"]]$shape_pearson_ci_low,
  shape_pearson_adult_ci_high = shape_rows[["2"]]$shape_pearson_ci_high,
  shape_spearman_adult = shape_rows[["2"]]$shape_spearman,
  shape_spearman_adult_ci_low = shape_rows[["2"]]$shape_spearman_ci_low,
  shape_spearman_adult_ci_high = shape_rows[["2"]]$shape_spearman_ci_high,
  nrmse_adult = shape_rows[["2"]]$normalized_rmse,
  nrmse_adult_ci_low = shape_rows[["2"]]$normalized_rmse_ci_low,
  nrmse_adult_ci_high = shape_rows[["2"]]$normalized_rmse_ci_high,
  smooth_boot_success = pmin(shape_rows[["1"]]$smooth_boot_success,
                             shape_rows[["2"]]$smooth_boot_success)
)]
brainspan_events[, `:=`(
  shape_pearson = rowMeans(cbind(shape_pearson_prenatal,
                                 shape_pearson_adult), na.rm = TRUE),
  shape_spearman = rowMeans(cbind(shape_spearman_prenatal,
                                  shape_spearman_adult), na.rm = TRUE),
  normalized_rmse = rowMeans(cbind(nrmse_prenatal, nrmse_adult), na.rm = TRUE),
  endpoint_concordance_prenatal =
    sign(smooth_point[, max(which(supported$segment_id == 1L))] -
         smooth_point[, min(which(supported$segment_id == 1L))]) ==
    sign(trajectory_mat[trajectory_index, max(which(supported$segment_id == 1L))] -
         trajectory_mat[trajectory_index, min(which(supported$segment_id == 1L))]),
  endpoint_concordance_adult =
    sign(smooth_point[, max(which(supported$segment_id == 2L))] -
         smooth_point[, min(which(supported$segment_id == 2L))]) ==
    sign(trajectory_mat[trajectory_index, max(which(supported$segment_id == 2L))] -
         trajectory_mat[trajectory_index, min(which(supported$segment_id == 2L))]),
  uninformative_smooth = smooth_boot_success <
    MIN_SMOOTH_BOOT_SUCCESS * SMOOTH_BOOT_ITER |
    (shape_pearson_prenatal_ci_high - shape_pearson_prenatal_ci_low) >
      MAX_SHAPE_CI_WIDTH |
    (shape_pearson_adult_ci_high - shape_pearson_adult_ci_low) >
      MAX_SHAPE_CI_WIDTH
)]
brainspan_events[, event_q_hc3 := p.adjust(p_hc3_global, method = "BH")]

# --- BrainSpan phenotype-specific mismatched-template null ---

brainspan_match_dt <- brainspan_events[, .(
  cluster_id,
  outcome_mean = rowMeans(brainspan_y),
  outcome_variance = apply(brainspan_y, 1L, var),
  missing_fraction = 1 - n_finite / ncol(brainspan_y),
  coverage = n_finite,
  pas_distance_bp
)]
brainspan_wrong <- assign_mismatched_template(
  brainspan_match_dt,
  c("outcome_mean", "outcome_variance", "missing_fraction", "coverage",
    "pas_distance_bp")
)
wrong_fit <- fit_template_rows(brainspan_y, brainspan_template[brainspan_wrong, ],
                               global_nuisance)
brainspan_events[, `:=`(
  mismatched_site_id = site_id[brainspan_wrong],
  same_similarity = partial_r_global,
  mismatched_similarity = wrong_fit$partial_r,
  same_minus_mismatched = partial_r_global - wrong_fit$partial_r
)]
brainspan_events[, prenatal_classification := fcase(
  !estimable | !is.finite(prenatal_ci_low) | !is.finite(prenatal_ci_high),
    "not_estimable",
  prenatal_ci_high < 0, "reversal",
  prenatal_ci_low > 0, "preservation",
  default = "non_replication"
)]
brainspan_events[, adult_classification := fcase(
  !estimable | !is.finite(adult_ci_low) | !is.finite(adult_ci_high),
    "not_estimable",
  adult_ci_high < 0, "reversal",
  adult_ci_low > 0, "preservation",
  default = "non_replication"
)]
brainspan_events[, classification := paste0(
  "prenatal:", prenatal_classification,
  "|adult:", adult_classification
)]

# --- BrainGVEX outcome-blind identifiability and event models ---

adult_index <- which(supported$grid_age >= ADULT_MIN_AGE &
                     supported$grid_age <= ADULT_MAX_AGE)
braingvex_samples <- scan(BRAINGVEX_SAMPLE_FILE, what = "character", quiet = TRUE)
braingvex_meta <- fread(BRAINGVEX_META_FILE)[match(braingvex_samples, BID)]
braingvex_header <- names(fread(BRAINGVEX_MATRIX_FILE, nrows = 0))
braingvex_select <- c(1L, match(braingvex_samples, braingvex_header))
braingvex <- fread(BRAINGVEX_MATRIX_FILE, select = braingvex_select,
                   showProgress = FALSE)
setnames(braingvex, 1L, "site_id")
braingvex <- merge(dynamic, braingvex, by = "site_id")
gvex_traj_index <- match(braingvex$site_id, rownames(trajectory_mat))
braingvex_y <- as.matrix(braingvex[, ..braingvex_samples])
storage.mode(braingvex_y) <- "numeric"
stopifnot(any(braingvex_y < 0, na.rm = TRUE),
          any(braingvex_y > 1, na.rm = TRUE))

gvex_beta <- gvex_se <- gvex_similarity <- gvex_estimable <- matrix(
  NA_real_, nrow(braingvex), uniqueN(braingvex_meta$batch),
  dimnames = list(NULL, sort(unique(braingvex_meta$batch)))
)
identifiability_rows <- list()
batch_templates <- list()
for (batch in colnames(gvex_beta)) {
  sample_index <- which(braingvex_meta$batch == batch)
  weight <- interpolation_weights(
    supported$grid_age[adult_index],
    braingvex_meta$AgeDeath[sample_index]
  )
  template <- trajectory_mat[gvex_traj_index, adult_index, drop = FALSE] %*% weight
  batch_templates[[batch]] <- template
  nuisance <- braingvex_nuisance(braingvex_meta[sample_index])
  fit <- fit_template_rows(braingvex_y[, sample_index, drop = FALSE],
                           template, nuisance)
  template_range <- apply(template, 1L, function(x) diff(range(x)))
  template_sd <- apply(template, 1L, sd)
  age_span <- diff(range(braingvex_meta$AgeDeath[sample_index]))
  estimable <- length(sample_index) >= MIN_BATCH_N &
    age_span >= MIN_BATCH_AGE_SPAN &
    template_range >= MIN_TEMPLATE_RANGE & template_sd >= MIN_TEMPLATE_SD &
    fit$template_residual_ss >= MIN_TEMPLATE_RESIDUAL_SS &
    fit$template_vif <= MAX_TEMPLATE_VIF
  gvex_beta[, batch] <- fit$beta
  gvex_se[, batch] <- fit$se_hc3
  gvex_similarity[, batch] <- fit$partial_r
  gvex_estimable[, batch] <- estimable
  identifiability_rows[[batch]] <- data.table(
    cohort = "BrainGVEX", site_id = braingvex$site_id,
    cluster_id = braingvex$cluster_id, batch = batch,
    n_samples = length(sample_index), age_span = age_span,
    template_range = template_range, template_sd = template_sd,
    template_residual_ss = fit$template_residual_ss,
    template_vif = fit$template_vif, estimable = estimable,
    reason = fcase(
      rep(length(sample_index) < MIN_BATCH_N, nrow(braingvex)), "batch_n",
      rep(age_span < MIN_BATCH_AGE_SPAN, nrow(braingvex)), "age_span",
      template_range < MIN_TEMPLATE_RANGE, "template_range",
      template_sd < MIN_TEMPLATE_SD, "template_sd",
      fit$template_residual_ss < MIN_TEMPLATE_RESIDUAL_SS, "residual_ss",
      fit$template_vif > MAX_TEMPLATE_VIF, "collinearity",
      default = "pass"
    )
  )
}
gvex_meta <- meta_analyze(gvex_beta, gvex_se, gvex_estimable == 1)
braingvex_events <- cbind(
  braingvex[, .(site_id, transcript_id, gene_id, cluster_id, cluster_label,
                 delta_logit, amplitude_logit)], gvex_meta
)
braingvex_events[, `:=`(
  mapped_qc = TRUE,
  estimable = n_estimable_batches >= MIN_META_BATCHES,
  beta = beta_meta
)]

# --- BrainGVEX secondary smooth shape uncertainty ---

gvex_shape_age <- seq(ADULT_MIN_AGE, ADULT_MAX_AGE, length.out = 60L)
gvex_smooth_fit <- fit_gvex_smooth_rows(
  braingvex_y, braingvex_meta, gvex_shape_age, BRAINGVEX_SMOOTH_DF
)
gvex_template_grid <- trajectory_mat[gvex_traj_index, adult_index, drop = FALSE] %*%
  interpolation_weights(supported$grid_age[adult_index], gvex_shape_age)
gvex_shape <- row_shape_metrics(gvex_template_grid, gvex_smooth_fit$prediction)
gvex_shape_uncertainty <- bootstrap_gvex_shape(
  braingvex_y, braingvex_meta, gvex_smooth_fit$design,
  gvex_smooth_fit$prediction_design, gvex_template_grid, SMOOTH_BOOT_ITER
)
braingvex_events[, `:=`(
  shape_pearson = gvex_shape$shape_pearson,
  shape_spearman = gvex_shape$shape_spearman,
  normalized_rmse = gvex_shape$normalized_rmse,
  shape_pearson_ci_low = gvex_shape_uncertainty$shape_pearson_ci_low,
  shape_pearson_ci_high = gvex_shape_uncertainty$shape_pearson_ci_high,
  shape_spearman_ci_low = gvex_shape_uncertainty$shape_spearman_ci_low,
  shape_spearman_ci_high = gvex_shape_uncertainty$shape_spearman_ci_high,
  normalized_rmse_ci_low = gvex_shape_uncertainty$normalized_rmse_ci_low,
  normalized_rmse_ci_high = gvex_shape_uncertainty$normalized_rmse_ci_high,
  smooth_boot_success = gvex_shape_uncertainty$smooth_boot_success,
  uninformative_smooth = gvex_shape_uncertainty$smooth_boot_success <
    MIN_SMOOTH_BOOT_SUCCESS * SMOOTH_BOOT_ITER |
    (gvex_shape_uncertainty$shape_pearson_ci_high -
      gvex_shape_uncertainty$shape_pearson_ci_low) > MAX_SHAPE_CI_WIDTH,
  endpoint_concordance =
    sign(gvex_smooth_fit$prediction[, ncol(gvex_smooth_fit$prediction)] -
         gvex_smooth_fit$prediction[, 1L]) ==
    sign(gvex_template_grid[, ncol(gvex_template_grid)] -
         gvex_template_grid[, 1L]),
  event_q_meta = p.adjust(p_meta, method = "BH")
)]

# --- BrainGVEX controls and leave-one-batch-out sensitivity ---

control_beta <- control_se <- matrix(NA_real_, nrow(braingvex), ncol(gvex_beta))
colnames(control_beta) <- colnames(control_se) <- colnames(gvex_beta)
for (batch in colnames(gvex_beta)) {
  sample_index <- which(braingvex_meta$batch == batch &
                        braingvex_meta$Diagnosis == "Control")
  template <- trajectory_mat[gvex_traj_index, adult_index, drop = FALSE] %*%
    interpolation_weights(supported$grid_age[adult_index],
                          braingvex_meta$AgeDeath[sample_index])
  nuisance <- braingvex_nuisance(
    braingvex_meta[sample_index], include_diagnosis = FALSE
  )
  fit <- fit_template_rows(braingvex_y[, sample_index, drop = FALSE],
                           template, nuisance)
  control_beta[, batch] <- fit$beta
  control_se[, batch] <- fit$se_hc3
}
control_meta <- meta_analyze(control_beta, control_se,
                             is.finite(control_beta) & is.finite(control_se))
lobo <- sapply(seq_len(ncol(gvex_beta)), function(drop_batch) {
  meta_analyze(gvex_beta[, -drop_batch, drop = FALSE],
               gvex_se[, -drop_batch, drop = FALSE],
               gvex_estimable[, -drop_batch, drop = FALSE] == 1)$beta_meta
})
braingvex_events[, `:=`(
  controls_beta_meta = control_meta$beta_meta,
  controls_ci_low = control_meta$ci_low,
  controls_ci_high = control_meta$ci_high,
  lobo_beta_min = row_min_finite(lobo),
  lobo_beta_max = row_max_finite(lobo)
)]

# --- BrainGVEX phenotype-specific mismatched-template null ---

gvex_match_dt <- data.table(
  cluster_id = braingvex$cluster_id,
  outcome_mean = rowMeans(braingvex_y),
  outcome_variance = apply(braingvex_y, 1L, var),
  missing_fraction = rowMeans(!is.finite(braingvex_y)),
  coverage = rowSums(is.finite(braingvex_y)),
  effect_size = abs(braingvex$delta_logit)
)
gvex_wrong <- assign_mismatched_template(
  gvex_match_dt,
  c("outcome_mean", "outcome_variance", "missing_fraction", "coverage",
    "effect_size")
)
wrong_beta <- wrong_se <- matrix(NA_real_, nrow(braingvex), ncol(gvex_beta))
colnames(wrong_beta) <- colnames(wrong_se) <- colnames(gvex_beta)
for (batch in colnames(gvex_beta)) {
  sample_index <- which(braingvex_meta$batch == batch)
  nuisance <- braingvex_nuisance(braingvex_meta[sample_index])
  wrong_template <- batch_templates[[batch]][gvex_wrong, , drop = FALSE]
  fit <- fit_template_rows(braingvex_y[, sample_index, drop = FALSE],
                           wrong_template, nuisance)
  wrong_beta[, batch] <- fit$partial_r
  wrong_se[, batch] <- fit$se_hc3
}
same_similarity <- rowSums(gvex_similarity * gvex_estimable, na.rm = TRUE) /
  rowSums(gvex_estimable, na.rm = TRUE)
wrong_similarity <- rowSums(wrong_beta * gvex_estimable, na.rm = TRUE) /
  rowSums(gvex_estimable, na.rm = TRUE)
braingvex_events[, `:=`(
  mismatched_site_id = site_id[gvex_wrong],
  same_similarity = same_similarity,
  mismatched_similarity = wrong_similarity,
  same_minus_mismatched = same_similarity - wrong_similarity,
  classification = fcase(
    !estimable, "not_estimable_brainGVEX_adult",
    ci_high < 0, "reversal",
    ci_low > 0 & same_similarity > wrong_similarity, "phenotype_preserved",
    ci_low > 0, "cluster_shape_only",
    default = "non_replication"
  )
)]

# --- mapping and identifiability outputs ---

mapping_output <- rbindlist(list(
  mapping[, .(
    cohort = "BrainSpan", site_id, external_site_id = brainspan_site_id,
    transcript_id, gene_id, cluster_id, exact_transcript_version,
    strand_pass, chromosome_pass,
    coordinate_bridge_pass = bridge_pass, pas_distance_bp, mapped_qc,
    matching_key,
    transcript_overlap_bp = overlap_bp,
    transcript_midpoint_distance = midpoint_distance,
    mapping_reason = "gene_txbase_max_transcript_overlap"
  )],
  braingvex[, .(
    cohort = "BrainGVEX", site_id, external_site_id = site_id,
    transcript_id, gene_id, cluster_id,
    exact_transcript_version = TRUE, strand_pass = TRUE,
    chromosome_pass = TRUE, coordinate_bridge_pass = TRUE,
    pas_distance_bp = 0, mapped_qc = TRUE, mapping_reason = "exact_site_id"
  )]
), fill = TRUE)
identifiability <- rbindlist(c(
  list(brainspan_events[, .(
    cohort = "BrainSpan", site_id, cluster_id, batch = NA_character_,
    n_samples = ncol(brainspan_y),
    age_span = diff(range(brainspan_meta$age_birth_aligned)),
    template_range, template_sd,
    template_residual_ss = template_residual_ss_global,
    template_vif = template_vif_global, estimable,
    reason = fifelse(estimable, "pass", "template_or_design_gate")
  )]),
  identifiability_rows
), fill = TRUE)

# --- artifact sensitivities ---

brainspan_raw <- fread(BRAINSPAN_RAW_FILE, showProgress = FALSE)
brainspan_raw[, brainspan_site_id := paste(
  parse_field(Gene, 1L), parse_field(Gene, 2L), Loci,
  parse_field(Gene, 4L), sep = "|"
)]
sensitivity_rows <- list()
sens_index <- 0L
for (tolerance in PAS_SENSITIVITY) {
  feature <- mapping[
    bridge_pass & (is.infinite(tolerance) | pas_distance_bp <= tolerance)
  ]
  feature_trajectory_index <- match(feature$site_id, rownames(trajectory_mat))
  feature_template <- matrix(NA_real_, nrow(feature), length(samples))
  for (segment in unique(brainspan_meta$segment_id)) {
    sample_index <- which(brainspan_meta$segment_id == segment)
    grid_index <- which(supported$segment_id == segment)
    feature_template[, sample_index] <- trajectory_mat[
      feature_trajectory_index, grid_index, drop = FALSE
    ] %*% interpolation_weights(
      supported$grid_age[grid_index],
      brainspan_meta$age_birth_aligned[sample_index]
    )
  }
  imputed_y <- as.matrix(brainspan_logit[
    match(feature$brainspan_site_id, brainspan_logit$brainspan_site_id), ..samples
  ])
  raw_y <- as.matrix(brainspan_raw[
    match(feature$brainspan_site_id, brainspan_raw$brainspan_site_id), ..samples
  ])
  storage.mode(imputed_y) <- storage.mode(raw_y) <- "numeric"
  raw_y <- pmin(pmax(raw_y, PDUI_EPSILON), 1 - PDUI_EPSILON)
  raw_y <- log2(raw_y / (1 - raw_y))
  matrix_values <- list(imputed_logit = imputed_y, raw_pdui_logit = raw_y)
  for (matrix_id in c("imputed_logit", "raw_pdui_logit")) {
    matrix_y <- matrix_values[[matrix_id]]
    fit <- fit_template_rows(matrix_y, feature_template, global_nuisance)
    for (cluster in sort(unique(feature$cluster_id))) {
      index <- which(feature$cluster_id == cluster & is.finite(fit$beta))
      gene <- data.table(gene_id = feature$gene_id[index],
                         beta = fit$beta[index])[, .(beta = median(beta)), by = gene_id]
      sens_index <- sens_index + 1L
      sensitivity_rows[[sens_index]] <- data.table(
        cohort = "BrainSpan", analysis = "mapping_preprocessing",
        matrix_id = matrix_id, pas_tolerance_bp = tolerance,
        cluster_id = cluster, n_events = length(index), n_genes = nrow(gene),
        estimate = median(gene$beta), detail = "gene_equal_median_beta"
      )
    }
  }
}
for (cluster in sort(unique(brainspan_events$cluster_id))) {
  event <- brainspan_events[cluster_id == cluster & estimable == TRUE]
  top_cutoff <- quantile(abs(event$delta_logit), 0.95)
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainSpan", analysis = "drop_top_discovery_effect",
    matrix_id = "imputed_logit", pas_tolerance_bp = REFERENCE_PAS_TOLERANCE,
    cluster_id = cluster, n_events = sum(abs(event$delta_logit) < top_cutoff),
    n_genes = uniqueN(event[abs(delta_logit) < top_cutoff, gene_id]),
    estimate = median(event[abs(delta_logit) < top_cutoff, beta_global]),
    detail = "drop_top_5pct_abs_discovery_delta"
  )
  primary_index <- which(mapping_primary$cluster_id == cluster)
  nonsaturated <- rowMeans(abs(brainspan_y[primary_index, , drop = FALSE]) >= 9.9) < 0.5
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainSpan", analysis = "exclude_saturated_events",
    matrix_id = "imputed_logit", pas_tolerance_bp = REFERENCE_PAS_TOLERANCE,
    cluster_id = cluster, n_events = sum(nonsaturated),
    n_genes = uniqueN(mapping_primary$gene_id[primary_index][nonsaturated]),
    estimate = median(brainspan_events$beta_global[primary_index][nonsaturated]),
    detail = "exclude_events_ge50pct_abs_logit_9.9"
  )
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainSpan", analysis = "event_equal_weighting",
    matrix_id = "imputed_logit", pas_tolerance_bp = REFERENCE_PAS_TOLERANCE,
    cluster_id = cluster, n_events = nrow(event), n_genes = uniqueN(event$gene_id),
    estimate = median(event$beta_global), detail = "event_equal_median_beta"
  )
}
for (cluster in sort(unique(braingvex_events$cluster_id))) {
  event <- braingvex_events[cluster_id == cluster & estimable == TRUE]
  top_cutoff <- quantile(abs(event$delta_logit), 0.95)
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", analysis = "controls_only", matrix_id = NA_character_,
    pas_tolerance_bp = NA_real_, cluster_id = cluster,
    n_events = nrow(event), n_genes = uniqueN(event$gene_id),
    estimate = median(event$controls_beta_meta, na.rm = TRUE),
    detail = "controls_fixed_effect_meta"
  )
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", analysis = "leave_one_batch_out",
    matrix_id = NA_character_, pas_tolerance_bp = NA_real_,
    cluster_id = cluster, n_events = nrow(event), n_genes = uniqueN(event$gene_id),
    estimate = median(event$lobo_beta_min, na.rm = TRUE),
    detail = sprintf("median_min_to_max=%g:%g",
      median(event$lobo_beta_min, na.rm = TRUE),
      median(event$lobo_beta_max, na.rm = TRUE))
  )
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", analysis = "drop_top_discovery_effect",
    matrix_id = NA_character_, pas_tolerance_bp = NA_real_,
    cluster_id = cluster, n_events = sum(abs(event$delta_logit) < top_cutoff),
    n_genes = uniqueN(event[abs(delta_logit) < top_cutoff, gene_id]),
    estimate = median(event[abs(delta_logit) < top_cutoff, beta], na.rm = TRUE),
    detail = "drop_top_5pct_abs_discovery_delta"
  )
  gvex_index <- which(braingvex_events$cluster_id == cluster &
                      braingvex_events$estimable == TRUE)
  nonsaturated <- rowMeans(abs(braingvex_y[gvex_index, , drop = FALSE]) >= 9.9) < 0.5
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", analysis = "exclude_saturated_events",
    matrix_id = "imputed_pdui_logit", pas_tolerance_bp = NA_real_,
    cluster_id = cluster, n_events = sum(nonsaturated),
    n_genes = uniqueN(braingvex_events$gene_id[gvex_index][nonsaturated]),
    estimate = median(braingvex_events$beta[gvex_index][nonsaturated]),
    detail = "exclude_events_ge50pct_abs_logit_9.9"
  )
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", analysis = "event_equal_weighting",
    matrix_id = "imputed_pdui_logit", pas_tolerance_bp = NA_real_,
    cluster_id = cluster, n_events = nrow(event), n_genes = uniqueN(event$gene_id),
    estimate = median(event$beta), detail = "event_equal_median_beta"
  )
}
sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)

# --- selection-aware discovery cross-fitting ---

run_discovery_crossfit <- function() {
  cov <- fread(DISCOVERY_COV_FILE)
  cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]
  apa <- if (CROSSFIT_MAX_SITES > 0L) {
    fread(DISCOVERY_MATRIX_FILE, nrows = CROSSFIT_MAX_SITES,
          showProgress = FALSE)
  } else {
    fread(DISCOVERY_MATRIX_FILE, showProgress = FALSE)
  }
  setnames(apa, 1L, "site_id")
  sample_ids <- intersect(cov$sample, names(apa))
  cov <- cov[match(sample_ids, sample)]
  y_all <- as.matrix(apa[, ..sample_ids])
  rownames(y_all) <- apa$site_id
  individuals <- unique(cov[, .(individual_id, age, sex_binary)])
  setorder(individuals, sex_binary, age)
  individuals[, fold := rep(seq_len(CROSSFIT_FOLDS), length.out = .N),
              by = sex_binary]
  cov[, fold := individuals$fold[match(individual_id, individuals$individual_id)]]
  delta_threshold <- unique(na.omit(dynamic$delta_logit_threshold))[1]
  grid_df <- data.frame(
    age_model = supported$age_model,
    sex_binary = weighted.mean(cov$sex_binary, cov$profile_weight),
    RIN_z = 0
  )
  fold_rows <- list()
  for (fold in seq_len(CROSSFIT_FOLDS)) {
    train <- which(cov$fold != fold)
    test <- which(cov$fold == fold)
    model_df <- as.data.frame(cov[train, .(
      age_model, sex_binary, RIN_z, profile_weight
    )])
    fit_one <- function(i) {
      dat <- model_df
      dat$y <- y_all[i, train]
      fit <- tryCatch(gam(
        y ~ sex_binary + RIN_z + s(age_model, k = CROSSFIT_BASIS_K, bs = "tp"),
        data = dat, weights = profile_weight, method = "REML"
      ), error = function(e) NULL)
      if (is.null(fit)) return(c(p = NA, rep(NA, nrow(supported))))
      c(p = summary(fit)$s.table[1, "p-value"],
        predict(fit, newdata = grid_df))
    }
    fit_list <- if (CROSSFIT_WORKERS > 1L) {
      parallel::mclapply(seq_len(nrow(y_all)), fit_one,
                         mc.cores = CROSSFIT_WORKERS)
    } else {
      lapply(seq_len(nrow(y_all)), fit_one)
    }
    fit_mat <- do.call(rbind, fit_list)
    qvalue <- p.adjust(fit_mat[, 1L], method = "BH")
    pred <- fit_mat[, -1L, drop = FALSE]
    delta <- pred[, ncol(pred)] - pred[, 1L]
    selected <- which(qvalue < CROSSFIT_Q & abs(delta) >= delta_threshold)
    pred_z <- t(scale(t(pred[selected, , drop = FALSE])))
    pca <- prcomp(pred_z, rank. = min(12L, ncol(pred_z)))
    cluster_raw <- kmeans(pca$x, centers = 2L, nstart = 8L,
                          iter.max = 500L, algorithm = "Lloyd")$cluster
    rank_label <- rank(-tapply(delta[selected], cluster_raw, median))
    cluster <- as.integer(rank_label[as.character(cluster_raw)])
    selected_ids <- rownames(y_all)[selected]
    held_template <- matrix(NA_real_, length(selected), length(test))
    for (segment in unique(supported$segment_id)) {
      test_segment <- which(cov$age[test] >= min(supported[segment_id == segment,
        grid_age]) & cov$age[test] <= max(supported[segment_id == segment, grid_age]))
      if (length(test_segment)) {
        grid_index <- which(supported$segment_id == segment)
        held_template[, test_segment] <- pred[selected, grid_index, drop = FALSE] %*%
          interpolation_weights(supported$grid_age[grid_index],
                                cov$age[test[test_segment]])
      }
    }
    valid_test <- which(colSums(is.finite(held_template)) == nrow(held_template))
    held_y <- y_all[selected, test[valid_test], drop = FALSE]
    held_template <- held_template[, valid_test, drop = FALSE]
    nuisance <- model.matrix(~ sex_binary + RIN_z,
                             data = cov[test[valid_test]])
    held_fit <- fit_template_rows(held_y, held_template, nuisance)
    full_match <- dynamic[match(selected_ids, site_id)]
    for (cluster_id in 1:2) {
      index <- which(cluster == cluster_id)
      gene <- data.table(
        gene_id = parse_field(selected_ids[index], 2L),
        beta = held_fit$beta[index], rmse = held_fit$residual_sd[index]
      )[, .(beta = median(beta), rmse = median(rmse)), by = gene_id]
      fold_rows[[length(fold_rows) + 1L]] <- data.table(
        fold = fold, cluster_id = cluster_id, n_train_individuals =
          uniqueN(cov$individual_id[train]), n_test_individuals =
          uniqueN(cov$individual_id[test]), n_selected_events = length(index),
        n_selected_genes = nrow(gene), heldout_beta = median(gene$beta),
        heldout_rmse = median(gene$rmse),
        heldout_direction_concordance = mean(gene$beta > 0),
        cluster_stability = mean(full_match$cluster_id[index] == cluster_id,
                                 na.rm = TRUE)
      )
    }
  }
  rbindlist(fold_rows)
}

discovery_crossfit <- if (RUN_DISCOVERY_CROSSFIT) {
  run_discovery_crossfit()
} else {
  data.table(
    fold = NA_integer_, cluster_id = NA_integer_,
    n_train_individuals = NA_integer_, n_test_individuals = NA_integer_,
    n_selected_events = NA_integer_, n_selected_genes = NA_integer_,
    heldout_beta = NA_real_, heldout_rmse = NA_real_,
    heldout_direction_concordance = NA_real_, cluster_stability = NA_real_,
    status = "not_run_APA_EXTERNAL_RUN_CROSSFIT_0"
  )
}

# --- representative trajectory overlays ---

representative_brainspan <- brainspan_events[estimable == TRUE][
  order(cluster_id, -same_minus_mismatched, p_hc3_global), .SD[1],
  by = cluster_id
]
overlay_rows <- list()
for (i in seq_len(nrow(representative_brainspan))) {
  event_index <- match(representative_brainspan$site_id[i],
                       brainspan_events$site_id)
  discovery_z <- as.numeric(scale(trajectory_mat[trajectory_index[event_index], ]))
  external_point <- smooth_point[event_index, ]
  external_scale <- sd(external_point)
  external_center <- mean(external_point)
  external_z <- (external_point - external_center) / external_scale
  external_ci_low <- external_ci_high <- rep(NA_real_, nrow(supported))
  for (segment in sort(unique(brainspan_meta$segment_id))) {
    sample_index <- which(brainspan_meta$segment_id == segment)
    grid_index <- which(supported$segment_id == segment)
    prediction_boot <- matrix(NA_real_, SMOOTH_BOOT_ITER, length(grid_index))
    sex_groups <- split(sample_index, brainspan_meta$sex_binary[sample_index])
    for (b in seq_len(SMOOTH_BOOT_ITER)) {
      boot_index <- unlist(lapply(sex_groups, sample, replace = TRUE),
                           use.names = FALSE)
      prediction_boot[b, ] <- tryCatch(
        fit_smooth_rows(
          brainspan_y[event_index, boot_index, drop = FALSE],
          brainspan_meta$age_birth_aligned[boot_index],
          brainspan_meta$sex_binary[boot_index],
          supported$grid_age[grid_index], BRAINSPAN_SMOOTH_DF
        ), error = function(e) rep(NA_real_, length(grid_index))
      )
    }
    external_ci_low[grid_index] <-
      apply(prediction_boot, 2L, quantile, 0.025, na.rm = TRUE)
    external_ci_high[grid_index] <-
      apply(prediction_boot, 2L, quantile, 0.975, na.rm = TRUE)
  }
  overlay_rows[[length(overlay_rows) + 1L]] <- data.table(
    cohort = "BrainSpan", cluster_id = representative_brainspan$cluster_id[i],
    cluster_label = representative_brainspan$cluster_label[i],
    site_id = representative_brainspan$site_id[i],
    segment_id = supported$segment_id, age = supported$grid_age,
    discovery_z = discovery_z, external_z = external_z,
    external_ci_low = (external_ci_low - external_center) / external_scale,
    external_ci_high = (external_ci_high - external_center) / external_scale
  )
}
# The manuscript trajectory overlay is prespecified to BrainSpan; BrainGVEX
# remains available in its dedicated event table but is not plotted here.
representative_gvex <- braingvex_events[FALSE]
for (i in seq_len(nrow(representative_gvex))) {
  event_index <- match(representative_gvex$site_id[i], braingvex_events$site_id)
  y <- as.numeric(braingvex_y[event_index, ])
  plot_data <- copy(braingvex_meta)
  plot_data[, y := y]
  fit <- gam(y ~ factor(batch) + factor(Sex) + scale(RIN) + scale(PMI) +
               factor(Diagnosis) + s(AgeDeath, k = 4), data = plot_data,
             method = "REML")
  age_grid <- sort(unique(seq(ADULT_MIN_AGE, ADULT_MAX_AGE, length.out = 60)))
  newdata <- data.frame(
    batch = factor(names(sort(table(plot_data$batch), decreasing = TRUE))[1],
                   levels = levels(factor(plot_data$batch))),
    Sex = factor(names(sort(table(plot_data$Sex), decreasing = TRUE))[1],
                 levels = levels(factor(plot_data$Sex))),
    RIN = mean(plot_data$RIN), PMI = mean(plot_data$PMI),
    Diagnosis = factor("Control", levels = levels(factor(plot_data$Diagnosis))),
    AgeDeath = age_grid
  )
  prediction <- predict(fit, newdata = newdata, se.fit = TRUE)
  template <- trajectory_mat[gvex_traj_index[event_index], adult_index] %*%
    interpolation_weights(supported$grid_age[adult_index], age_grid)
  external_center <- mean(prediction$fit)
  external_scale <- sd(prediction$fit)
  overlay_rows[[length(overlay_rows) + 1L]] <- data.table(
    cohort = "BrainGVEX", cluster_id = representative_gvex$cluster_id[i],
    cluster_label = representative_gvex$cluster_label[i],
    site_id = representative_gvex$site_id[i], segment_id = 2L, age = age_grid,
    discovery_z = as.numeric(scale(as.numeric(template))),
    external_z = (prediction$fit - external_center) / external_scale,
    external_ci_low = (prediction$fit - 1.96 * prediction$se.fit -
      external_center) / external_scale,
    external_ci_high = (prediction$fit + 1.96 * prediction$se.fit -
      external_center) / external_scale
  )
}
overlay <- rbindlist(overlay_rows)

# --- save outputs ---

brainspan_events[, smooth_boot_iter := SMOOTH_BOOT_ITER]
braingvex_events[, smooth_boot_iter := SMOOTH_BOOT_ITER]
discovery_crossfit[, `:=`(
  crossfit_requested = RUN_DISCOVERY_CROSSFIT,
  crossfit_max_sites = CROSSFIT_MAX_SITES,
  crossfit_scope = fifelse(
    !RUN_DISCOVERY_CROSSFIT, "not_run",
    fifelse(CROSSFIT_MAX_SITES > 0L, "smoke_subset", "full_discovery_matrix")
  )
)]
sensitivity[, smooth_boot_iter := SMOOTH_BOOT_ITER]
overlay[, smooth_boot_iter := SMOOTH_BOOT_ITER]

add_input_provenance <- function(dt) {
  dt[, `:=`(
    validation_version = "brainspan_v2_raw_hg19_liftover",
    brainspan_matrix_file = BRAINSPAN_RAW_FILE,
    brainspan_matrix_scale = "raw PDUI transformed to log2(PDUI/(1-PDUI))",
    brainspan_source_build = "hg19",
    brainspan_target_build = "GRCh38",
    brainspan_gencode_version = GENCODE_VERSION,
    brainspan_liftover_chain = BRAINSPAN_CHAIN_FILE,
    braingvex_matrix_file = BRAINGVEX_MATRIX_FILE,
    braingvex_matrix_scale = BRAINGVEX_MATRIX_SCALE
  )]
  dt
}

mapping_output <- add_input_provenance(mapping_output)
identifiability <- add_input_provenance(identifiability)
brainspan_events <- add_input_provenance(brainspan_events)
braingvex_events <- add_input_provenance(braingvex_events)
discovery_crossfit <- add_input_provenance(discovery_crossfit)
sensitivity <- add_input_provenance(sensitivity)
overlay <- add_input_provenance(overlay)

setorder(mapping_output, cohort, cluster_id, site_id)
setorder(identifiability, cohort, cluster_id, site_id, batch)
setorder(brainspan_events, cluster_id, gene_id, site_id)
setorder(braingvex_events, cluster_id, gene_id, site_id)
setorder(sensitivity, cohort, cluster_id, analysis)
setorder(overlay, cohort, cluster_id, segment_id, age)

fwrite(mapping_output, OUT_MAPPING, sep = "\t")
fwrite(identifiability, OUT_IDENTIFIABILITY, sep = "\t")
fwrite(brainspan_events, OUT_BRAINSPAN, sep = "\t")
fwrite(braingvex_events, OUT_BRAINGVEX, sep = "\t")
fwrite(discovery_crossfit, OUT_CROSSFIT, sep = "\t")
fwrite(sensitivity, OUT_SENSITIVITY, sep = "\t")
fwrite(overlay, OUT_OVERLAY, sep = "\t")
fwrite(liftover_qc, OUT_LIFTOVER_QC, sep = "\t")

message("done.")
