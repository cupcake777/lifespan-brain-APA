# Individual-weighted functional analysis of the converged trajectory clusters.
suppressPackageStartupMessages({
  library(data.table)
  library(splines)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
UPSTREAM <- file.path(root, "data/01_quant")
REF <- file.path(root, "ref")

# --- settings ---

MIN_SET_SIZE <- 10L
MAX_SET_SIZE <- 500L
CANDIDATE_LIBRARY_FDR <- 0.10
REPORT_ADJUSTED_FDR <- 0.05
ENRICHMENT_ONLY <- identical(Sys.getenv("APA_ENRICHMENT_ONLY", "0"), "1")

# --- helpers ---

strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))

safe_spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 5L || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

row_weighted_cor <- function(a, b, weights) {
  weights <- weights / sum(weights)
  a <- a - as.numeric(a %*% weights)
  b <- b - as.numeric(b %*% weights)
  denom <- sqrt(rowSums(a^2 * rep(weights, each = nrow(a))) *
                  rowSums(b^2 * rep(weights, each = nrow(b))))
  out <- rowSums(a * b * rep(weights, each = nrow(a))) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

fit_rank_grid <- function(age_vec, y_vec, weight_vec, grid_vec, max_deg = 3L) {
  ok <- is.finite(age_vec) & is.finite(y_vec) & is.finite(weight_vec)
  if (sum(ok) < 5L) return(rep(NA_real_, length(grid_vec)))
  y_rank <- rank(y_vec[ok], ties.method = "average")
  weights <- weight_vec[ok]
  degree <- min(max_deg, sum(ok) - 2L)
  x <- outer(age_vec[ok], 0:degree, `^`)
  coef <- solve(crossprod(x, x * weights), crossprod(x, y_rank * weights))
  pred <- Reduce("+", lapply(0:degree, function(i) coef[i + 1L] * grid_vec^i))
  as.numeric(scale(pred))
}

load_rows <- function(path, wanted, strip_versions = FALSE) {
  wanted <- unique(as.character(wanted))
  id_file <- tempfile()
  out_file <- tempfile(fileext = ".tsv")
  writeLines(wanted, id_file)
  awk_program <- if (strip_versions) {
    "FNR==NR{ids[$1]=1;next}FNR==1{print;next}{key=$1;sub(/\\.[0-9]+$/ ,\"\",key);if(key in ids)print}"
  } else {
    "FNR==NR{ids[$1]=1;next}FNR==1||($1 in ids)"
  }
  system(sprintf("awk %s %s %s > %s", shQuote(awk_program),
                 shQuote(id_file), shQuote(path), shQuote(out_file)))
  out <- fread(out_file, showProgress = FALSE)
  file.remove(id_file, out_file)
  out
}

donor_mean_rows <- function(dt, id_col, profile_cov, donor_order) {
  profile_samples <- intersect(profile_cov$sample, names(dt))
  sample_cov <- profile_cov[match(profile_samples, sample)]
  mat <- as.matrix(dt[, ..profile_samples])
  storage.mode(mat) <- "numeric"
  stopifnot(all(is.finite(mat)))

  donor_counts <- sample_cov[, .N, by = individual_id]
  donor_by_feature <- rowsum(t(mat), sample_cov$individual_id, reorder = FALSE)
  donor_by_feature <- donor_by_feature / donor_counts$N
  donor_mat <- t(donor_by_feature)
  colnames(donor_mat) <- donor_counts$individual_id
  donors <- intersect(donor_order, donor_counts$individual_id)

  out <- cbind(data.table(feature_id = dt[[id_col]]),
               as.data.table(donor_mat[, donors, drop = FALSE]))
  setnames(out, "feature_id", id_col)
  out
}

residualize_rows <- function(mat, sample_cov) {
  x <- model.matrix(~ ns(age_model, df = 5) + sex_binary + RIN_z, data = sample_cov)
  weights <- sample_cov$profile_weight
  beta <- mat %*% (x * weights) %*% solve(crossprod(x, x * weights))
  mat - beta %*% t(x)
}

read_gmt <- function(path) {
  lines <- readLines(path)
  sets <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][-(1:2)])
  setNames(sets, vapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][1],
                        character(1)))
}

run_enrichment <- function(fg_sets, background, gene_sets, library_name, source_group, level) {
  bg <- intersect(unique(background), unique(unlist(gene_sets)))
  eligible <- gene_sets[vapply(gene_sets, function(x) {
    n <- length(intersect(x, bg))
    n >= MIN_SET_SIZE && n <= MAX_SET_SIZE
  }, logical(1))]

  out <- rbindlist(lapply(names(fg_sets), function(group_name) {
    fg <- intersect(unique(fg_sets[[group_name]]), bg)
    ref <- setdiff(bg, fg)
    rbindlist(lapply(names(eligible), function(term_name) {
      term_genes <- intersect(eligible[[term_name]], bg)
      a <- length(intersect(fg, term_genes))
      b <- length(fg) - a
      c <- length(intersect(ref, term_genes))
      d <- length(ref) - c
      test <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      data.table(
        library_name = library_name,
        source_group = source_group,
        level = level,
        gene_set = group_name,
        term = term_name,
        pvalue = test$p.value,
        fisher_odds_ratio = unname(test$estimate),
        n_fg = a,
        n_bg = c,
        n_fg_total = length(fg),
        n_bg_total = length(bg)
      )
    }))
  }))
  out[, library_fdr := p.adjust(pvalue, "BH"), by = gene_set]
  out
}

# --- clusters and gene mapping ---
# APA_ENRICHMENT_ONLY=1 reuses the saved gene table and runs enrichment only.

if (!ENRICHMENT_ONLY) {
  sites <- fread(file.path(DATA, "02_sites_clustered.tsv"))
  support <- fread(file.path(DATA, "01_support.tsv"))
  traj <- fread(file.path(DATA, "traj_logit.csv"))
  cov <- fread(file.path(DATA, "00_cov.tsv"))
  profile_cov <- fread(file.path(DATA, "00_cov_profiles.tsv"))
  cov[, age_model := log10(pmax(as.numeric(age) + 1, 1e-6))]

  clustered <- sites[directional_pass == TRUE & !is.na(cluster_id)]
  clustered[, transcript_id := sub("\\|.*$", "", site_id)]
  clustered[, tx_base := strip_ver(transcript_id)]

  tx2g <- fread(file.path(REF, "gencode.v40.tx2gene.tsv"))
  tx2g[, `:=`(
    tx_base = strip_ver(tx_name),
    gene_base = strip_ver(gene_id),
    gene_length = end - start + 1
  )]
  tx_map <- unique(tx2g[, .(tx_base, transcript_id = tx_name, gene_symbol = gene_name,
                             gene_id, gene_base)])
  gene_span <- tx2g[, .(gene_length = max(end) - min(start) + 1),
                    by = .(gene_symbol = gene_name, gene_id_base = gene_base)]

  clustered <- merge(clustered, tx_map, by = "tx_base", all.x = TRUE, suffixes = c("", "_map"))
  clustered[is.na(gene_symbol) | gene_symbol == "", gene_symbol := gene_id_map]
  clustered[, gene_base := strip_ver(fifelse(
    !is.na(gene_id_map) & gene_id_map != "", gene_id_map, gene_id
  ))]

  # --- two-level gene architecture ---

  gene_agg <- clustered[gene_symbol != "", .(
    gene_id_base = gene_base[1],
    n_sites = uniqueN(site_id),
    n_c1 = sum(cluster_id == 1L),
    n_c2 = sum(cluster_id == 2L)
  ), by = gene_symbol]
  gene_agg[, gene_state := fcase(
    n_c1 > 0 & n_c2 == 0, "C1-spec",
    n_c1 > n_c2 & n_c2 > 0, "C1-main",
    n_c1 == n_c2 & n_c1 > 0, "Balanced",
    n_c2 > n_c1 & n_c1 > 0, "C2-main",
    n_c2 > 0 & n_c1 == 0, "C2-spec"
  )]
  gene_agg[, membership := fifelse(gene_state %chin% c("C1-spec", "C2-spec"), "Specific", "Shared")]
  fwrite(gene_agg, file.path(DATA, "05_gene_aggregation.tsv"), sep = "\t")

  # --- expression rows and mean expression ---

  expr_rows <- load_rows(file.path(UPSTREAM, "expr_after_combat.tsv"),
                         unique(gene_agg$gene_id_base), strip_versions = TRUE)
  setnames(expr_rows, 1L, "gene_id")
  expr_rows[, gene_base := strip_ver(gene_id)]
  expr_rows <- donor_mean_rows(expr_rows, "gene_id", profile_cov, cov$sample)
  expr_rows[, gene_base := strip_ver(gene_id)]
  expr_samples <- intersect(names(expr_rows), cov$sample)
  expr_weights <- cov[match(expr_samples, sample), profile_weight]
  expr_matrix <- as.matrix(expr_rows[, ..expr_samples])
  storage.mode(expr_matrix) <- "numeric"
  expr_rows[, mean_expression := as.numeric(expr_matrix %*% expr_weights / sum(expr_weights))]
  expr_mean <- unique(expr_rows[, .(gene_base, mean_expression)], by = "gene_base")

  # --- supported-grid age-trend coupling ---

  traj_cols <- grep("^grid_", names(traj), value = TRUE)
  grid_index <- as.integer(sub("grid_", "", traj_cols)) + 1L
  supported_cols <- traj_cols[support$is_supported[grid_index]]
  supported_grid <- support$age_model[support$is_supported]

  iso_rows <- load_rows(file.path(UPSTREAM, "isoform_after_combat.tsv"),
                        unique(clustered$transcript_id))
  setnames(iso_rows, 1L, "transcript_id")
  iso_rows <- donor_mean_rows(iso_rows, "transcript_id", profile_cov, cov$sample)
  iso_samples <- intersect(names(iso_rows), cov$sample)

  iso_grid <- setNames(lapply(seq_len(nrow(iso_rows)), function(i) {
    fit_rank_grid(
      cov[match(iso_samples, sample), age_model],
      as.numeric(iso_rows[i, ..iso_samples]),
      cov[match(iso_samples, sample), profile_weight],
      supported_grid
    )
  }), iso_rows$transcript_id)

  expr_grid <- setNames(lapply(seq_len(nrow(expr_rows)), function(i) {
    fit_rank_grid(
      cov[match(expr_samples, sample), age_model],
      as.numeric(expr_rows[i, ..expr_samples]),
      cov[match(expr_samples, sample), profile_weight],
      supported_grid
    )
  }), expr_rows$gene_base)

  apa_supported <- traj[, c("site_id", supported_cols), with = FALSE]
  apa_grid <- setNames(lapply(seq_len(nrow(apa_supported)), function(i) {
    as.numeric(apa_supported[i, ..supported_cols])
  }), apa_supported$site_id)

  trend <- rbindlist(lapply(seq_len(nrow(clustered)), function(i) {
    row <- clustered[i]
    apa_t <- apa_grid[[row$site_id]]
    tx_t <- iso_grid[[row$transcript_id]]
    gene_t <- expr_grid[[row$gene_base]]
    data.table(
      site_id = row$site_id,
      cluster_id = row$cluster_id,
      transcript_id = row$transcript_id,
      gene_symbol = row$gene_symbol,
      gene_base = row$gene_base,
      apa_tx_spearman = if (is.null(tx_t)) NA_real_ else safe_spearman(apa_t, tx_t),
      tx_gene_spearman = if (is.null(tx_t) || is.null(gene_t)) NA_real_ else safe_spearman(tx_t, gene_t),
      apa_gene_spearman = if (is.null(gene_t)) NA_real_ else safe_spearman(apa_t, gene_t)
    )
  }))
  fwrite(trend, file.path(DATA, "06_cross_layer_trend_outputs.tsv"), sep = "\t")

  # --- sample-level covariate-residual coupling ---

  common_samples <- Reduce(intersect, list(
    cov$sample,
    names(fread(file.path(UPSTREAM, "apa_combat_donor_mean.tsv"), nrows = 0)),
    names(iso_rows),
    names(expr_rows)
  ))
  sample_cov <- cov[match(common_samples, sample)]

  apa_rows <- load_rows(file.path(UPSTREAM, "apa_combat_donor_mean.tsv"), clustered$site_id)
  setnames(apa_rows, 1L, "site_id")

  apa_order <- match(clustered$site_id, apa_rows$site_id)
  iso_order <- match(clustered$transcript_id, iso_rows$transcript_id)
  expr_order <- match(clustered$gene_base, expr_rows$gene_base)
  complete <- is.finite(apa_order) & is.finite(iso_order) & is.finite(expr_order)

  apa_mat <- as.matrix(apa_rows[apa_order[complete], ..common_samples])
  iso_mat <- as.matrix(iso_rows[iso_order[complete], ..common_samples])
  expr_mat <- as.matrix(expr_rows[expr_order[complete], ..common_samples])
  storage.mode(apa_mat) <- storage.mode(iso_mat) <- storage.mode(expr_mat) <- "numeric"

  apa_resid <- residualize_rows(apa_mat, sample_cov)
  iso_resid <- residualize_rows(iso_mat, sample_cov)
  expr_resid <- residualize_rows(expr_mat, sample_cov)

  residual <- clustered[complete, .(site_id, cluster_id, transcript_id, gene_symbol, gene_base)]
  residual[, `:=`(
    apa_tx_residual_r = row_weighted_cor(apa_resid, iso_resid, sample_cov$profile_weight),
    tx_gene_residual_r = row_weighted_cor(iso_resid, expr_resid, sample_cov$profile_weight),
    apa_gene_residual_r = row_weighted_cor(apa_resid, expr_resid, sample_cov$profile_weight),
    n_profiles = length(common_samples),
    n_individuals = uniqueN(sample_cov$individual_id)
  )]
  fwrite(residual, file.path(DATA, "06_cross_layer_residual.tsv"), sep = "\t")

  sync_summary <- merge(
    trend[, .(
      n_sites = .N,
      med_apa_tx_trend = median(apa_tx_spearman, na.rm = TRUE),
      med_tx_gene_trend = median(tx_gene_spearman, na.rm = TRUE),
      med_apa_gene_trend = median(apa_gene_spearman, na.rm = TRUE)
    ), by = cluster_id],
    residual[, .(
      n_residual_sites = .N,
      med_apa_tx_residual = median(apa_tx_residual_r, na.rm = TRUE),
      med_tx_gene_residual = median(tx_gene_residual_r, na.rm = TRUE),
      med_apa_gene_residual = median(apa_gene_residual_r, na.rm = TRUE)
    ), by = cluster_id],
    by = "cluster_id"
  )
  fwrite(sync_summary, file.path(DATA, "04_synchrony_summary.tsv"), sep = "\t")

  neutral_sensitivity <- rbindlist(lapply(c(0.05, 0.10, 0.20), function(cutoff) {
    rbindlist(lapply(c("apa_tx_spearman", "tx_gene_spearman", "apa_gene_spearman"), function(metric) {
      trend[, .(
        n = sum(is.finite(get(metric))),
        fraction_negative = mean(get(metric) < -cutoff, na.rm = TRUE),
        fraction_neutral = mean(abs(get(metric)) <= cutoff, na.rm = TRUE),
        fraction_positive = mean(get(metric) > cutoff, na.rm = TRUE)
      ), by = cluster_id][, `:=`(metric = metric, neutral_cutoff = cutoff)]
    }))
  }))
  fwrite(neutral_sensitivity, file.path(DATA, "06_cross_layer_neutral_sens.tsv"), sep = "\t")

  # --- expression direction by gene architecture ---

  expr_age <- expr_grid[gene_agg$gene_id_base]
  expr_rho <- vapply(expr_age, function(x) safe_spearman(supported_grid, x), numeric(1))
  expr_apa <- gene_agg[, .(gene_symbol, gene_id_base, gene_state, membership)]
  expr_apa[, expr_age_rho := expr_rho]
  expr_apa[, expr_direction := fifelse(expr_age_rho > 0, "up", fifelse(expr_age_rho < 0, "down", "flat"))]
  fwrite(expr_apa, file.path(DATA, "08_expr_apa_consistency.tsv"), sep = "\t")

  # --- covariate-adjusted LOEUF ---

  constraint <- fread(
    cmd = sprintf("zcat %s", shQuote(file.path(REF, "gnomad.v4.1.1.constraint_metrics.tsv.bgz"))),
    select = c("gene", "canonical", "mane_select", "lof.oe_ci.upper")
  )
  setnames(constraint, c("gene", "canonical", "mane_select", "lof.oe_ci.upper"),
           c("gene_symbol", "canonical", "mane_select", "loeuf"))
  constraint[, loeuf := as.numeric(loeuf)]
  constraint <- constraint[is.finite(loeuf)]
  constraint[, priority := fifelse(as.logical(mane_select), 1L, fifelse(as.logical(canonical), 2L, 3L))]
  constraint <- constraint[, .SD[priority == min(priority)], by = gene_symbol]
  constraint <- constraint[, .(loeuf = median(unique(loeuf))), by = gene_symbol]

  loeuf_dt <- merge(gene_agg, constraint, by = "gene_symbol")
  loeuf_dt <- merge(loeuf_dt, gene_span, by = c("gene_symbol", "gene_id_base"), all.x = TRUE)
  loeuf_dt <- merge(loeuf_dt, expr_mean, by.x = "gene_id_base", by.y = "gene_base", all.x = TRUE)
  loeuf_dt <- loeuf_dt[is.finite(loeuf) & is.finite(gene_length) & is.finite(mean_expression)]
  loeuf_dt[, `:=`(
    gene_state = factor(gene_state, levels = c("C1-spec", "C1-main", "Balanced", "C2-main", "C2-spec")),
    membership = factor(membership, levels = c("Specific", "Shared"))
  )]
  fwrite(loeuf_dt, file.path(DATA, "ed4c_loeuf.tsv"), sep = "\t")

  membership_model <- lm(
    log1p(loeuf) ~ membership + scale(log1p(n_sites)) + scale(log1p(gene_length)) + scale(mean_expression),
    data = loeuf_dt
  )
  state_model <- lm(
    log1p(loeuf) ~ gene_state + scale(log1p(n_sites)) + scale(log1p(gene_length)) + scale(mean_expression),
    data = loeuf_dt
  )
  membership_coef <- as.data.table(coef(summary(membership_model)), keep.rownames = "term")
  membership_coef[, model := "membership"]
  state_coef <- as.data.table(coef(summary(state_model)), keep.rownames = "term")
  state_coef[, model := "five_state"]
  model_table <- rbindlist(list(membership_coef, state_coef), fill = TRUE)
  setnames(model_table, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"),
           c("estimate", "std_error", "t_value", "pvalue"))
  model_table[, fdr := p.adjust(pvalue, "BH"), by = model]
  fwrite(model_table, file.path(DATA, "11_loeuf_adj_models_outputs.tsv"), sep = "\t")
} else {
  gene_agg <- fread(file.path(DATA, "05_gene_aggregation.tsv"))
  cov <- fread(file.path(DATA, "00_cov.tsv"))
  profile_cov <- fread(file.path(DATA, "00_cov_profiles.tsv"))

  tx2g <- fread(file.path(REF, "gencode.v40.tx2gene.tsv"))
  tx2g[, `:=`(
    gene_base = strip_ver(gene_id),
    gene_length = end - start + 1
  )]
  gene_span <- tx2g[, .(
    gene_length = max(end) - min(start) + 1
  ), by = .(gene_symbol = gene_name, gene_id_base = gene_base)]

  expr_rows <- load_rows(file.path(UPSTREAM, "expr_after_combat.tsv"),
                         unique(gene_agg$gene_id_base), strip_versions = TRUE)
  setnames(expr_rows, 1L, "gene_id")
  expr_rows[, gene_base := strip_ver(gene_id)]
  expr_rows <- donor_mean_rows(expr_rows, "gene_id", profile_cov, cov$sample)
  expr_rows[, gene_base := strip_ver(gene_id)]
  expr_samples <- intersect(names(expr_rows), cov$sample)
  expr_weights <- cov[match(expr_samples, sample), profile_weight]
  expr_matrix <- as.matrix(expr_rows[, ..expr_samples])
  storage.mode(expr_matrix) <- "numeric"
  expr_rows[, mean_expression := as.numeric(expr_matrix %*% expr_weights / sum(expr_weights))]
  expr_mean <- unique(expr_rows[, .(gene_base, mean_expression)], by = "gene_base")
}

# --- unadjusted and covariate-adjusted enrichment ---

gwas_lines <- readLines(file.path(REF, "GWAS_Catalog_2025.txt"))
gwas_sets <- setNames(
  lapply(gwas_lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]][-(1:2)]),
  vapply(gwas_lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]][1], character(1))
)
annotation_libraries <- list(
  "GO BP" = read_gmt(file.path(REF, "c5.go.bp.v2024.1.Hs.symbols.gmt")),
  "GO MF" = read_gmt(file.path(REF, "c5.go.mf.v2024.1.Hs.symbols.gmt")),
  "KEGG" = read_gmt(file.path(REF, "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")),
  "GWAS Catalog" = gwas_sets,
  "HPO" = read_gmt(file.path(REF, "Human_Phenotype_Ontology.gmt"))
)
source_groups <- c(
  "GO BP" = "GO/KEGG", "GO MF" = "GO/KEGG", "KEGG" = "GO/KEGG",
  "GWAS Catalog" = "GWAS/HPO", "HPO" = "GWAS/HPO"
)

membership_sets <- split(gene_agg$gene_symbol, gene_agg$membership)
state_sets <- split(gene_agg$gene_symbol, gene_agg$gene_state)
background <- gene_agg$gene_symbol

enrichment <- rbindlist(lapply(names(annotation_libraries), function(library_name) {
  gene_sets <- annotation_libraries[[library_name]]
  source_group <- unname(source_groups[library_name])
  rbindlist(list(
    run_enrichment(membership_sets, background, gene_sets, library_name, source_group, "Membership"),
    run_enrichment(state_sets, background, gene_sets, library_name, source_group, "Five-state")
  ))
}))

covariates <- merge(gene_agg, gene_span, by = c("gene_symbol", "gene_id_base"), all.x = TRUE)
covariates <- merge(covariates, expr_mean, by.x = "gene_id_base", by.y = "gene_base", all.x = TRUE)
covariates <- covariates[is.finite(gene_length) & is.finite(mean_expression)]

candidate_terms <- unique(enrichment[library_fdr < CANDIDATE_LIBRARY_FDR,
  .(library_name, source_group, level, gene_set, term)])
adjusted_rows <- lapply(seq_len(nrow(candidate_terms)), function(i) {
  hit <- candidate_terms[i]
  gene_sets <- annotation_libraries[[hit$library_name]]
  dt <- copy(covariates)
  dt[, in_group := as.integer(if (hit$level == "Membership") {
    membership == hit$gene_set
  } else {
    gene_state == hit$gene_set
  })]
  dt[, in_term := as.integer(gene_symbol %chin% gene_sets[[hit$term]])]
  fit <- try(glm(
    in_group ~ in_term + scale(log1p(n_sites)) + scale(log1p(gene_length)) + scale(mean_expression),
    family = binomial(), data = dt
  ), silent = TRUE)
  if (inherits(fit, "try-error") || !"in_term" %in% names(coef(fit))) return(NULL)
  coef_table <- coef(summary(fit))
  if (!"in_term" %in% rownames(coef_table)) return(NULL)
  data.table(
    library_name = hit$library_name,
    source_group = hit$source_group,
    level = hit$level,
    gene_set = hit$gene_set,
    term = hit$term,
    n_genes = nrow(dt),
    adjusted_odds_ratio = exp(coef(fit)["in_term"]),
    adjusted_pvalue = coef_table["in_term", "Pr(>|z|)"]
  )
})
adjusted <- rbindlist(adjusted_rows, fill = TRUE)
if (nrow(adjusted) > 0L) {
  adjusted[, adjusted_fdr := p.adjust(adjusted_pvalue, "BH"), by = .(library_name, level, gene_set)]
}
result_keys <- c("library_name", "source_group", "level", "gene_set", "term")
adjusted <- merge(
  enrichment[, c(result_keys, "fisher_odds_ratio", "library_fdr"), with = FALSE],
  adjusted, by = result_keys
)
adjusted[, report_pass := adjusted_fdr < REPORT_ADJUSTED_FDR]
enrichment <- merge(
  enrichment,
  adjusted[, c(result_keys, "adjusted_odds_ratio", "adjusted_pvalue", "adjusted_fdr", "report_pass"),
           with = FALSE],
  by = result_keys, all.x = TRUE
)
fwrite(enrichment, file.path(DATA, "08_enrich_two_levels.tsv"), sep = "\t")
fwrite(adjusted, file.path(DATA, "08_enrich_adjusted_outputs.tsv"), sep = "\t")
