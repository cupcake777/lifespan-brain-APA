# Cluster supported-grid trajectories with seeded Lloyd k-means and robustness checks.
suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
  library(bluster)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")

# --- settings ---

MIN_K <- 2L
MAX_K <- 8L
MIN_CLUSTER_N <- as.integer(Sys.getenv("APA_CLUSTER_MIN_N", "30"))
N_PC <- 12L
KMEANS_STARTS <- 8L
KMEANS_ITER <- 500L
KMEANS_ALGORITHM <- "Lloyd"
BOOT_ITER <- as.integer(Sys.getenv("APA_CLUSTER_BOOT_ITER", "100"))
BOOT_FRAC <- 0.8
SIL_SAMPLE_PRIMARY <- 5000L
SIL_SAMPLE_SENSITIVITY <- 2000L
PRIMARY_SEED <- 2026L
SENSITIVITY_SEEDS <- c(2026L, 2027L, 2028L)
ROBUSTNESS_MIN <- 0.90

# --- helper functions ---

fit_kmeans <- function(x, k) {
  kmeans(x, centers = k, nstart = KMEANS_STARTS, iter.max = KMEANS_ITER,
         algorithm = KMEANS_ALGORITHM)
}

gene_states <- function(assignments) {
  assignments[, .(gene_state = paste(sort(unique(cluster_id)), collapse = "+")), by = gene_id]
}

classify_warning <- function(message) {
  if (grepl("did not converge|empty cluster|Quick-TRANSfer", message, ignore.case = TRUE)) {
    "result_affecting"
  } else {
    "harmless_numerical"
  }
}

run_cluster_analysis <- function(site_dt, mat, label, seed, boot_iter, sil_sample) {
  captured_warnings <- character()
  result <- withCallingHandlers({
    set.seed(seed)
    pca <- prcomp(mat, center = TRUE, scale. = FALSE, rank. = min(N_PC, ncol(mat)))
    pc <- pca$x[, seq_len(min(N_PC, ncol(pca$x))), drop = FALSE]
    k_range <- MIN_K:min(MAX_K, nrow(pc) - 1L)
    fits <- setNames(lapply(k_range, function(k) fit_kmeans(pc, k)), as.character(k_range))
    labels_by_k <- lapply(fits, `[[`, "cluster")

    sil_idx <- sample(nrow(pc), min(sil_sample, nrow(pc)))
    sil_scores <- sapply(k_range, function(k) {
      lab <- labels_by_k[[as.character(k)]][sil_idx]
      mean(silhouette(lab, dist(pc[sil_idx, , drop = FALSE]))[, "sil_width"])
    })

    if (boot_iter > 0L) {
      boot_agree <- sapply(k_range, function(k) {
        base <- labels_by_k[[as.character(k)]]
        median(replicate(boot_iter, {
          boot_n <- min(nrow(pc), max(round(nrow(pc) * BOOT_FRAC), k + 1L))
          idx <- sample(nrow(pc), boot_n)
          km <- fit_kmeans(pc[idx, , drop = FALSE], k)
          distance <- sapply(seq_len(k), function(j) {
            rowSums((pc[idx, , drop = FALSE] - matrix(
              km$centers[j, ], nrow = length(idx), ncol = ncol(pc), byrow = TRUE
            ))^2)
          })
          query_lab <- max.col(-distance)
          base_lab <- base[idx]
          tab <- table(query_lab, base_lab)
          mapped <- apply(tab, 1L, which.max)[as.character(query_lab)]
          mean(mapped == base_lab, na.rm = TRUE)
        }))
      })
    } else {
      boot_agree <- rep(NA_real_, length(k_range))
    }

    kdiag <- data.table(
      analysis = label,
      seed = seed,
      k = k_range,
      sample_silhouette = sil_scores,
      sample_bootstrap_median_agreement = boot_agree,
      min_cluster_size = sapply(k_range, function(k) min(tabulate(labels_by_k[[as.character(k)]])))
    )
    kdiag[, size_gate_pass := min_cluster_size >= MIN_CLUSTER_N]
    selected_k <- kdiag[size_gate_pass == TRUE][order(-sample_silhouette, k)]$k[1]
    if (length(selected_k) == 0L || is.na(selected_k)) {
      selected_k <- kdiag[order(-sample_silhouette, k)]$k[1]
    }
    kdiag[, auto_select := k == selected_k]

    labels0 <- labels_by_k[[as.character(selected_k)]]
    median_delta <- tapply(site_dt$delta_logit, labels0, median, na.rm = TRUE)
    label_rank <- rank(-median_delta, ties.method = "first")
    labels <- as.integer(label_rank[as.character(labels0)])
    silhouette_site <- as.numeric(bluster::approxSilhouette(pc, labels0)$width)

    agreement <- rep(NA_real_, nrow(pc))
    if (boot_iter > 0L) {
      agreement <- numeric(nrow(pc))
      for (b in seq_len(boot_iter)) {
        boot_n <- min(nrow(pc), max(round(nrow(pc) * BOOT_FRAC), selected_k + 1L))
        idx <- sample(nrow(pc), boot_n)
        km <- fit_kmeans(pc[idx, , drop = FALSE], selected_k)
        distance <- sapply(seq_len(selected_k), function(j) {
          rowSums((pc - matrix(
            km$centers[j, ], nrow = nrow(pc), ncol = ncol(pc), byrow = TRUE
          ))^2)
        })
        boot_lab <- max.col(-distance)
        tab <- table(boot_lab[idx], labels0[idx])
        mapped <- apply(tab, 1L, which.max)[as.character(boot_lab)]
        agreement <- agreement + (mapped == labels0)
      }
      agreement <- agreement / boot_iter
    }

    assignments <- site_dt[, .(site_id, gene_id, delta_logit, amplitude_logit)]
    assignments[, `:=`(
      cluster_id = labels,
      cluster_silhouette = silhouette_site,
      cluster_consensus = agreement
    )]

    summary <- assignments[, .(
      n_sites = .N,
      n_genes = uniqueN(gene_id),
      median_delta_logit = median(delta_logit),
      median_amplitude = median(amplitude_logit),
      median_silhouette = median(cluster_silhouette),
      median_consensus = if (all(is.na(cluster_consensus))) NA_real_ else
        median(cluster_consensus, na.rm = TRUE),
      direction = fifelse(median(delta_logit) > 0, "lengthening",
                          fifelse(median(delta_logit) < 0, "shortening", "transient")),
      locked_k = selected_k
    ), by = cluster_id][order(cluster_id)]

    pca_var <- data.table(analysis = label, pc = seq_along(pca$sdev), variance = pca$sdev^2)
    pca_var[, `:=`(
      variance_explained = variance / sum(variance),
      cumulative_variance_explained = cumsum(variance / sum(variance))
    )]

    list(assignments = assignments, summary = summary, kdiag = kdiag, pca_var = pca_var)
  }, warning = function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  result$warnings <- unique(captured_warnings)
  result
}

# --- inputs ---

sites <- fread(file.path(DATA, "02_sites.tsv"))
traj <- fread(file.path(DATA, "traj_z.csv"), showProgress = FALSE)
support <- fread(file.path(DATA, "01_support.tsv"))
mask_flags <- fread(file.path(DATA, "02_mask_site_flags.tsv"))
selected_masks <- fread(file.path(DATA, "01_support_selected_masks.tsv"))
selection <- fread(file.path(DATA, "01_support_selection_probability.tsv"))

all_traj_cols <- grep("^grid_", names(traj), value = TRUE)
grid_index <- as.integer(sub("grid_", "", all_traj_cols)) + 1L
baseline_id <- selected_masks[is_baseline == TRUE, candidate_id][1]
baseline_cols <- all_traj_cols[support$is_supported[grid_index]]

# --- primary baseline clustering ---

primary_sites <- sites[directional_pass == TRUE]
primary <- merge(primary_sites, traj[, c("site_id", baseline_cols), with = FALSE], by = "site_id")
primary_result <- run_cluster_analysis(
  primary, as.matrix(primary[, ..baseline_cols]), "directional",
  PRIMARY_SEED, BOOT_ITER, SIL_SAMPLE_PRIMARY
)

allshape_sites <- sites[dynamic_any_pass == TRUE]
allshape <- merge(allshape_sites, traj[, c("site_id", baseline_cols), with = FALSE], by = "site_id")
allshape_result <- run_cluster_analysis(
  allshape, as.matrix(allshape[, ..baseline_cols]), "all_shape",
  PRIMARY_SEED, BOOT_ITER, SIL_SAMPLE_PRIMARY
)

# --- primary outputs ---

sites[, `:=`(
  cluster_id = NA_integer_,
  cluster_silhouette = NA_real_,
  cluster_consensus = NA_real_,
  allshape_cluster_id = NA_integer_
)]
primary_index <- match(primary_result$assignments$site_id, sites$site_id)
sites[primary_index, `:=`(
  cluster_id = primary_result$assignments$cluster_id,
  cluster_silhouette = primary_result$assignments$cluster_silhouette,
  cluster_consensus = primary_result$assignments$cluster_consensus
)]
allshape_index <- match(allshape_result$assignments$site_id, sites$site_id)
sites[allshape_index, allshape_cluster_id := allshape_result$assignments$cluster_id]

fwrite(sites, file.path(DATA, "02_sites_clustered.tsv"), sep = "\t")
fwrite(primary_result$summary, file.path(DATA, "03_clusters.tsv"), sep = "\t")
fwrite(allshape_result$summary, file.path(DATA, "03_allshape_clusters.tsv"), sep = "\t")
fwrite(primary_result$kdiag, file.path(DATA, "03_k_diagnostics.tsv"), sep = "\t")
fwrite(allshape_result$kdiag, file.path(DATA, "03_allshape_k_diagnostics.tsv"), sep = "\t")
fwrite(rbind(primary_result$pca_var, allshape_result$pca_var),
       file.path(DATA, "03_pca_variance.tsv"), sep = "\t")
fwrite(primary_result$assignments, file.path(DATA, "03_consensus_annotation.tsv"), sep = "\t")

# --- neighbouring-mask and seed robustness ---

baseline_assignments <- primary_result$assignments[, .(site_id, gene_id, cluster_id)]
baseline_genes <- gene_states(baseline_assignments)
baseline_site_ids <- baseline_assignments$site_id
warning_rows <- list()
if (length(primary_result$warnings)) {
  warning_rows[[length(warning_rows) + 1L]] <- data.table(
    analysis = "primary_directional", candidate_id = baseline_id,
    seed = PRIMARY_SEED, message = primary_result$warnings
  )
}
if (length(allshape_result$warnings)) {
  warning_rows[[length(warning_rows) + 1L]] <- data.table(
    analysis = "primary_all_shape", candidate_id = baseline_id,
    seed = PRIMARY_SEED, message = allshape_result$warnings
  )
}

robustness_rows <- list()
row_index <- 0L
for (candidate in selected_masks$candidate_id) {
  supported <- selection[candidate_id == candidate][order(grid_index)]$full_supported
  supported_cols <- all_traj_cols[supported[grid_index]]
  candidate_flags <- mask_flags[candidate_id == candidate & directional_pass == TRUE]
  candidate_data <- merge(
    candidate_flags[, .(site_id, delta_logit, amplitude_logit)],
    sites[, .(site_id, gene_id)],
    by = "site_id"
  )
  candidate_data <- merge(candidate_data, traj[, c("site_id", supported_cols), with = FALSE],
                          by = "site_id")

  for (seed in SENSITIVITY_SEEDS) {
    if (candidate == baseline_id && seed == PRIMARY_SEED) {
      result <- primary_result
    } else {
      result <- run_cluster_analysis(
        candidate_data, as.matrix(candidate_data[, ..supported_cols]),
        paste(candidate, seed, sep = "_"), seed, 0L, SIL_SAMPLE_SENSITIVITY
      )
    }

    common <- merge(
      baseline_assignments[, .(site_id, baseline_cluster = cluster_id)],
      result$assignments[, .(site_id, candidate_cluster = cluster_id)],
      by = "site_id"
    )
    candidate_genes <- gene_states(result$assignments)
    gene_common <- merge(
      baseline_genes[, .(gene_id, baseline_state = gene_state)],
      candidate_genes[, .(gene_id, candidate_state = gene_state)],
      by = "gene_id"
    )

    row_index <- row_index + 1L
    robustness_rows[[row_index]] <- data.table(
      candidate_id = candidate,
      seed = seed,
      n_directional = nrow(result$assignments),
      dynamic_site_jaccard = {
        union_n <- length(union(baseline_site_ids, result$assignments$site_id))
        if (union_n == 0L) 1 else
          length(intersect(baseline_site_ids, result$assignments$site_id)) / union_n
      },
      selected_k = result$summary$locked_k[1],
      cluster_ari = {
        tab <- table(common$baseline_cluster, common$candidate_cluster)
        choose2 <- function(x) x * (x - 1) / 2
        sum_cells <- sum(choose2(tab))
        sum_rows <- sum(choose2(rowSums(tab)))
        sum_cols <- sum(choose2(colSums(tab)))
        expected <- sum_rows * sum_cols / choose2(sum(tab))
        denominator <- (sum_rows + sum_cols) / 2 - expected
        if (denominator == 0) 1 else (sum_cells - expected) / denominator
      },
      c1_proportion = mean(result$assignments$cluster_id == 1L),
      c2_proportion = mean(result$assignments$cluster_id == 2L),
      gene_state_agreement = mean(gene_common$baseline_state == gene_common$candidate_state),
      n_common_sites = nrow(common),
      n_common_genes = nrow(gene_common),
      n_warnings = length(result$warnings),
      n_result_affecting_warnings = sum(vapply(
        result$warnings, classify_warning, character(1)) == "result_affecting")
    )

    if (length(result$warnings)) {
      warning_rows[[length(warning_rows) + 1L]] <- data.table(
        analysis = "mask_seed_sensitivity", candidate_id = candidate,
        seed = seed, message = result$warnings
      )
    }
  }
}

robustness <- rbindlist(robustness_rows)
primary_k <- primary_result$summary$locked_k[1]
robustness[, gate_pass :=
  selected_k == primary_k & dynamic_site_jaccard >= ROBUSTNESS_MIN &
  cluster_ari >= ROBUSTNESS_MIN & n_result_affecting_warnings == 0L]
fwrite(robustness, file.path(DATA, "03_mask_seed_robustness.tsv"), sep = "\t")

gate <- data.table(
  primary_k = primary_k,
  selected_k_consistent = all(robustness$selected_k == primary_k),
  median_dynamic_site_jaccard = median(robustness$dynamic_site_jaccard),
  median_cluster_ari = median(robustness$cluster_ari),
  median_gene_state_agreement = median(robustness$gene_state_agreement),
  n_combinations = nrow(robustness),
  n_gate_pass = sum(robustness$gate_pass),
  robustness_gate_pass =
    all(robustness$selected_k == primary_k) &
    median(robustness$dynamic_site_jaccard) >= ROBUSTNESS_MIN &
    median(robustness$cluster_ari) >= ROBUSTNESS_MIN &
    sum(robustness$n_result_affecting_warnings) == 0L
)
fwrite(gate, file.path(DATA, "03_robustness_gate.tsv"), sep = "\t")

warning_table <- if (length(warning_rows)) {
  rbindlist(warning_rows, fill = TRUE)
} else {
  data.table(analysis = character(), candidate_id = character(), seed = integer(),
             message = character())
}
warning_table[, warning_class := vapply(message, classify_warning, character(1))]
fwrite(warning_table, file.path(DATA, "03_cluster_warnings.tsv"), sep = "\t")
