# Gene architecture, directional coherence against shuffled labels, and APA
# trajectories by gene state.
suppressPackageStartupMessages({
  library(patchwork)
})
source("fig_common.R")

# --- config ---

GENE_AGG_FILE <- file.path(OUT_DIR, "05_gene_aggregation.tsv")
EXPR_APA_FILE <- file.path(OUT_DIR, "08_expr_apa_consistency.tsv")
SITES_FILE <- file.path(OUT_DIR, "02_sites_clustered.tsv")
SUPPORT_FILE <- file.path(OUT_DIR, "01_support.tsv")
TRAJ_Z_FILE <- file.path(INTERMEDIATE_DIR, "traj_z.csv")
TX2GENE_FILE <- file.path(APA_ROOT, "ref", "gencode.v40.tx2gene.tsv")

cluster_levels <- c("C1", "C2")
concordance_levels <- c("Specific", "Shared")
gene_class_levels <- c("C1-spec", "C1-main", "Balanced",
                       "C2-main", "C2-spec")
expr_trend_levels <- c("Expression up", "Expression down")

SPECIFICITY_PERM_ITER <- 10000L
SPECIFICITY_PERM_SEED <- 20260914L

membership_colors <- c("C1-specific" = "#56B4E9",
                       "Shared" = "#9A9A9A",
                       "C2-specific" = "#E69F00")
gene_class_colors <- c("C1-spec" = "#56B4E9",
                       "C1-main" = "#0072B2",
                       "Balanced" = "#7A7A7A",
                       "C2-main" = "#D55E00",
                       "C2-spec" = "#E69F00")

gene_class_full_names <- c(
  "C1-spec" = "C1-specific",
  "C1-main" = "Shared, C1-dominant",
  "Balanced" = "Shared, balanced",
  "C2-main" = "Shared, C2-dominant",
  "C2-spec" = "C2-specific"
)

out_png <- file.path(PANEL_DIR, "07_gene_direction.png")
out_pdf <- file.path(PANEL_DIR, "07_gene_direction.pdf")
out_cluster <- file.path(OUT_DIR, "07_gene_direction_cluster_composition.tsv")
out_gene <- file.path(OUT_DIR, "07_gene_direction_gene_composition.tsv")
out_trajectory <- file.path(OUT_DIR, "07_gene_direction_traj.tsv")
out_stats <- file.path(OUT_DIR, "07_gene_direction_stats.tsv")
out_specificity <- file.path(OUT_DIR, "07_gene_direction_specificity_null.tsv")
out_specificity_perm <- file.path(OUT_DIR, "07_gene_direction_specificity_perms.tsv")

# --- loading data ---

gene_agg <- fread(GENE_AGG_FILE)
expr_apa <- fread(EXPR_APA_FILE, select = c("gene_symbol", "expr_age_rho"))
expr_apa <- unique(expr_apa[is.finite(expr_age_rho)], by = "gene_symbol")

gene_agg[, gene_class := fcase(
  n_c1 > 0 & n_c2 == 0, "C1-spec",
  n_c1 == 0 & n_c2 > 0, "C2-spec",
  n_c1 > n_c2, "C1-main",
  n_c1 == n_c2 & n_c1 > 0 & n_c2 > 0, "Balanced",
  n_c1 < n_c2, "C2-main"
)]
gene_agg[, gene_class := factor(gene_class, levels = gene_class_levels)]

# --- panel A data: cluster-centric composition, count-first ---

c1_specific_n <- gene_agg[n_c1 > 0 & n_c2 == 0, .N]
c2_specific_n <- gene_agg[n_c2 > 0 & n_c1 == 0, .N]
shared_n <- gene_agg[n_c1 > 0 & n_c2 > 0, .N]

cluster_comp <- data.table(
  cluster = factor(rep(cluster_levels, each = 2), levels = cluster_levels),
  concordance = factor(rep(concordance_levels, times = 2),
                       levels = concordance_levels),
  n_genes = c(c1_specific_n, shared_n, c2_specific_n, shared_n)
)
cluster_comp[, total_cluster_genes := sum(n_genes), by = cluster]
cluster_comp[, pct := n_genes / total_cluster_genes]
cluster_comp[, label := comma(n_genes)]
cluster_comp[, membership_state := factor(
  fcase(
    concordance == "Shared", "Shared",
    cluster == "C1", "C1-specific",
    cluster == "C2", "C2-specific"
  ),
  levels = c("C1-specific", "Shared", "C2-specific")
)]

# --- panel B data: 3'UTR trajectory by gene state ---

gene_comp <- gene_agg[, .(
  n_genes = .N,
  n_sites = sum(n_sites),
  n_c1_sites = sum(n_c1),
  n_c2_sites = sum(n_c2)
), by = gene_class]
complete_gene <- data.table(gene_class = factor(gene_class_levels,
                                                levels = gene_class_levels))
gene_comp <- merge(complete_gene, gene_comp, by = "gene_class", all.x = TRUE)
gene_comp[is.na(n_genes), `:=`(
  n_genes = 0L,
  n_sites = 0L,
  n_c1_sites = 0L,
  n_c2_sites = 0L
)]
gene_comp[, pct_all_genes := n_genes / sum(n_genes)]

# --- panel D data: specificity against fixed-margin label shuffles ---

gene_index <- rep(seq_len(nrow(gene_agg)), gene_agg$n_sites)
n_events <- length(gene_index)
n_c1_events <- sum(gene_agg$n_c1)

set.seed(SPECIFICITY_PERM_SEED)
specificity_perm <- data.table(
  iteration = seq_len(SPECIFICITY_PERM_ITER),
  c1_spec = integer(SPECIFICITY_PERM_ITER),
  c2_spec = integer(SPECIFICITY_PERM_ITER)
)
for (i in seq_len(SPECIFICITY_PERM_ITER)) {
  shuffled_c1 <- tabulate(
    gene_index[sample.int(n_events, n_c1_events)],
    nbins = nrow(gene_agg)
  )
  specificity_perm[i, `:=`(
    c1_spec = sum(shuffled_c1 == gene_agg$n_sites),
    c2_spec = sum(shuffled_c1 == 0L)
  )]
}
specificity_perm[, any_spec := c1_spec + c2_spec]

null_means <- specificity_perm[, .(
  any_spec = mean(any_spec),
  c1_spec = mean(c1_spec),
  c2_spec = mean(c2_spec)
)]
specificity_perm[, c2_c1_relative :=
  (c2_spec / null_means$c2_spec) / (c1_spec / null_means$c1_spec)]

observed_metrics <- c(
  any_spec = c1_specific_n + c2_specific_n,
  c1_spec = c1_specific_n,
  c2_spec = c2_specific_n,
  c2_c1_relative = (c2_specific_n / null_means$c2_spec) /
    (c1_specific_n / null_means$c1_spec)
)
null_metrics <- list(
  any_spec = specificity_perm$any_spec / null_means$any_spec,
  c1_spec = specificity_perm$c1_spec / null_means$c1_spec,
  c2_spec = specificity_perm$c2_spec / null_means$c2_spec,
  c2_c1_relative = specificity_perm$c2_c1_relative
)

specificity_null <- rbindlist(lapply(names(null_metrics), function(metric_id) {
  null_ratio <- null_metrics[[metric_id]]
  observed_ratio <- if (metric_id == "c2_c1_relative") {
    observed_metrics[[metric_id]]
  } else {
    observed_metrics[[metric_id]] / null_means[[metric_id]]
  }
  data.table(
    metric_id = metric_id,
    observed_count = if (metric_id == "c2_c1_relative") NA_real_ else
      observed_metrics[[metric_id]],
    null_mean_count = if (metric_id == "c2_c1_relative") NA_real_ else
      null_means[[metric_id]],
    observed_to_null = observed_ratio,
    null_ratio_q025 = quantile(null_ratio, 0.025),
    null_ratio_q975 = quantile(null_ratio, 0.975),
    empirical_p_greater = (1 + sum(null_ratio >= observed_ratio)) /
      (SPECIFICITY_PERM_ITER + 1)
  )
}))

specificity_labels <- c(
  any_spec = "Any-spec",
  c1_spec = "C1-spec",
  c2_spec = "C2-spec",
  c2_c1_relative = "C2:C1 specificity"
)
specificity_colors <- c(
  any_spec = "#333333",
  c1_spec = "#56B4E9",
  c2_spec = "#E69F00",
  c2_c1_relative = "#8E6BBE"
)
specificity_null[, metric_label := factor(
  specificity_labels[metric_id],
  levels = rev(unname(specificity_labels))
)]
specificity_null[, ratio_label := sprintf("%.2fx", observed_to_null)]
specificity_null[, is_primary := metric_id %chin%
                   c("any_spec", "c2_c1_relative")]
specificity_null[, plot_label := fifelse(
  is_primary,
  sprintf("%s\nPemp < 1e-4", ratio_label),
  ratio_label
)]

expr_df <- merge(
  gene_agg[, .(gene_symbol, n_c1, n_c2, gene_class)],
  expr_apa,
  by = "gene_symbol"
)
expr_df[, expr_trend := fifelse(expr_age_rho > 0,
                                "Expression up", "Expression down")]
expr_df[, expr_trend := factor(expr_trend, levels = expr_trend_levels)]
expr_df[, utr_trajectory := factor(as.character(gene_class),
                                   levels = gene_class_levels)]
expr_df[, apa_side := fcase(
  gene_class %chin% c("C1-spec", "C1-main"), "C1 3'UTR lengthening",
  gene_class %chin% c("C2-main", "C2-spec"), "C2 3'UTR shortening",
  gene_class == "Balanced", "Balanced"
)]
expr_df[, coupling := fcase(
  apa_side == "C1 3'UTR lengthening" & expr_trend == "Expression up", "Coupled",
  apa_side == "C2 3'UTR shortening" & expr_trend == "Expression down", "Coupled",
  apa_side == "Balanced", "Balanced",
  default = "Decoupled"
)]

expr_summary <- expr_df[, .(
  n_expr_genes = .N,
  pct_expression_up = mean(expr_trend == "Expression up"),
  median_expr_rho = median(expr_age_rho),
  q25_expr_rho = quantile(expr_age_rho, 0.25),
  q75_expr_rho = quantile(expr_age_rho, 0.75)
), by = gene_class]
gene_comp <- merge(gene_comp, expr_summary, by = "gene_class", all.x = TRUE)
gene_comp[is.na(n_expr_genes), `:=`(
  n_expr_genes = 0L,
  pct_expression_up = NA_real_,
  median_expr_rho = NA_real_,
  q25_expr_rho = NA_real_,
  q75_expr_rho = NA_real_
)]

expr_table <- table(expr_df$gene_class, expr_df$expr_trend)
expr_chisq <- chisq.test(expr_table)

directional_df <- copy(expr_df[apa_side %in% c("C1 3'UTR lengthening",
                                               "C2 3'UTR shortening")])
directional_df[, coupling_binary := factor(coupling, levels = c("Coupled", "Decoupled"))]
directional_table <- table(directional_df$apa_side, directional_df$coupling_binary)
directional_fisher <- fisher.test(directional_table)
directional_odds <- if (length(directional_fisher$estimate) > 0) {
  unname(directional_fisher$estimate)
} else {
  NA_real_
}

support <- fread(SUPPORT_FILE)
sites_for_traj <- fread(SITES_FILE, select = c("site_id", "directional_pass",
                                               "cluster_id"))
sites_for_traj[, transcript_id := tstrsplit(site_id, "|", fixed = TRUE,
                                             keep = 1L)]
sites_for_traj[, tx_base := sub("\\.[0-9]+$", "", transcript_id)]

tx2g <- fread(TX2GENE_FILE, select = c("tx_name", "gene_name", "gene_id"))
tx2g[, `:=`(
  tx_base = sub("\\.[0-9]+$", "", tx_name),
  gene_id_base = sub("\\.[0-9]+$", "", gene_id)
)]
tx_map <- unique(tx2g[, .(tx_base, gene_symbol = gene_name, gene_id_base)])

traj_z <- fread(TRAJ_Z_FILE, showProgress = FALSE)
traj_cols <- grep("^grid_", names(traj_z), value = TRUE)
grid_idx <- as.integer(sub("grid_", "", traj_cols))
age_lookup <- support[match(grid_idx, grid_index), .(grid = traj_cols,
                                                     age_model,
                                                     grid_age,
                                                     is_supported)]

site_gene_map <- merge(
  sites_for_traj[directional_pass == TRUE & !is.na(cluster_id),
                 .(site_id, transcript_id, tx_base, cluster_id)],
  tx_map,
  by = "tx_base",
  all.x = TRUE
)
site_gene_map <- unique(site_gene_map,
                        by = c("site_id", "cluster_id", "tx_base",
                               "gene_symbol", "gene_id_base"))
site_gene_class <- merge(
  site_gene_map,
  gene_agg[, .(gene_symbol, gene_id_base, gene_class)],
  by = c("gene_symbol", "gene_id_base")
)
site_gene_class[, gene_class := factor(gene_class, levels = gene_class_levels)]

stopifnot(setequal(site_gene_class$gene_symbol, gene_agg$gene_symbol))

traj_dt <- merge(site_gene_class, traj_z, by = "site_id")
traj_long <- melt(
  traj_dt,
  id.vars = c("site_id", "gene_id_base", "cluster_id", "gene_class"),
  measure.vars = traj_cols,
  variable.name = "grid",
  value.name = "z"
)
traj_long <- merge(traj_long, age_lookup, by = "grid")

gene_trajectory <- traj_long[, .(
  z = median(z, na.rm = TRUE)
), by = .(gene_id_base, gene_class, age_model, grid_age, is_supported)]

trajectory_summary <- gene_trajectory[, .(
  n_genes = uniqueN(gene_id_base),
  median_z = median(z, na.rm = TRUE),
  q25_z = quantile(z, 0.25, na.rm = TRUE),
  q75_z = quantile(z, 0.75, na.rm = TRUE)
), by = .(gene_class, age_model, grid_age, is_supported)]

set.seed(SEED)
trace_genes <- unique(site_gene_class[, .(gene_id_base, gene_class)])[
  ,
  .SD[sample(.N, min(.N, 90L))],
  by = gene_class
]
trace_long <- merge(
  trace_genes,
  gene_trajectory[, .(gene_id_base, gene_class, age_model, z)],
  by = c("gene_id_base", "gene_class")
)

# --- panel C data: unique-gene states ---

gene_comp[, state_label := sprintf("%s (%.1f%%)", comma(n_genes),
                                    100 * pct_all_genes)]
gene_comp[, state_full_name := gene_class_full_names[as.character(gene_class)]]
gene_comp[, membership_level := fifelse(
  gene_class %chin% c("C1-spec", "C2-spec"), "Specific", "Shared"
)]

gene_value_max <- max(gene_comp$n_genes)
label_char_width <- 0.020 * gene_value_max
gene_comp[, estimated_label_right :=
            n_genes + nchar(state_label) * label_char_width]
all_label_right <- max(gene_comp$estimated_label_right)
gene_axis_max <- ceiling(max(1.5 * gene_value_max,
                             1.12 * all_label_right) / 500) * 500
shared_text_y <- 0.84 * gene_axis_max
shared_text_x <- 3

# --- figure 07A: compact cluster composition ---

p_cluster <- ggplot(cluster_comp, aes(cluster, n_genes, fill = membership_state,
                                      group = concordance)) +
  geom_col(width = 0.58, color = "white", linewidth = 0.35) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 3.25, lineheight = 0.92, color = "white",
            fontface = "bold") +
  geom_text(data = unique(cluster_comp[, .(cluster, total_cluster_genes)]),
            aes(x = cluster, y = total_cluster_genes,
                label = paste0("total ", comma(total_cluster_genes))),
            inherit.aes = FALSE, vjust = -0.45, size = 3.1,
            color = "#222222", fontface = "bold") +
  scale_fill_manual(values = membership_colors, name = NULL) +
  scale_y_continuous(labels = comma,
                     limits = c(0, max(cluster_comp$total_cluster_genes) * 1.14),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Genes",
       title = "A  Cluster membership") +
  theme_pub(base = 11) +
  theme(legend.position = "right",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 11, face = "bold"),
        plot.title = element_text(hjust = 0, face = "bold"),
        plot.margin = margin(4, 8, 4, 4))

# --- figure 07B: fitted 3'UTR trajectories ---

end_labels <- trajectory_summary[
  is_supported == TRUE,
  .SD[which.max(age_model)],
  by = gene_class
]
end_labels <- merge(
  end_labels,
  gene_comp[, .(gene_class, pct_expression_up)],
  by = "gene_class"
)
end_labels[, short_label := as.character(gene_class)]
end_labels[, label := sprintf("%s\nn=%s", short_label, comma(n_genes))]

trajectory_y_min <- floor((min(trajectory_summary$q25_z, na.rm = TRUE) - 0.12) * 10) / 10
trajectory_y_max <- ceiling((max(trajectory_summary$q75_z, na.rm = TRUE) + 0.12) * 10) / 10
trace_plot <- trace_long[z >= trajectory_y_min & z <= trajectory_y_max]

label_positions <- data.table(
  gene_class = factor(gene_class_levels, levels = gene_class_levels),
  label_y = c(-0.80, -0.52, 0.05, 0.54, 0.92)
)
end_labels <- merge(end_labels, label_positions, by = "gene_class")
end_labels[, label_x := max(age_model) + 0.06]
end_labels[, segment_x := max(age_model) + 0.015]
trajectory_x_max <- max(age_lookup$age_model) + 0.55

p_trajectory <- ggplot() +
  geom_hline(yintercept = 0, color = "grey78", linewidth = 0.28) +
  geom_line(data = trace_plot,
            aes(age_model, z, group = gene_id_base, color = gene_class),
            alpha = 0.018, linewidth = 0.14, show.legend = FALSE) +
  geom_ribbon(data = trajectory_summary,
              aes(age_model, ymin = q25_z, ymax = q75_z, fill = gene_class),
              alpha = 0.16, color = NA, show.legend = FALSE) +
  geom_line(data = trajectory_summary,
            aes(age_model, median_z, color = gene_class),
            linewidth = 1.08, show.legend = FALSE) +
  geom_point(data = end_labels,
             aes(age_model, median_z, color = gene_class),
             size = 1.7, show.legend = FALSE) +
  geom_segment(data = end_labels,
               aes(x = age_model, xend = segment_x,
                   y = median_z, yend = label_y, color = gene_class),
               linewidth = 0.35, alpha = 0.8, show.legend = FALSE) +
  geom_text(data = end_labels,
            aes(label_x, label_y, label = label, color = gene_class),
            hjust = 0, size = 3.15, lineheight = 0.86,
            fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = gene_class_colors, name = "3'UTR state") +
  scale_fill_manual(values = gene_class_colors, name = "3'UTR state") +
  scale_x_continuous(breaks = age_breaks, labels = age_labels,
                     limits = c(min(age_lookup$age_model), trajectory_x_max),
                     expand = expansion(mult = c(0.01, 0))) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(ylim = c(trajectory_y_min, trajectory_y_max), clip = "off") +
  labs(x = "Age (years)", y = "Fitted APA Z-score",
       title = "B  APA trajectories by gene state") +
  theme_pub(base = 12) +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.key.size = unit(0.34, "cm"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
        axis.text = element_text(size = 10.5),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(hjust = 0, face = "bold", size = 13.5),
        plot.margin = margin(8, 34, 8, 8))

# --- figure 07C: unique-gene state composition ---

p_gene <- ggplot(gene_comp, aes(gene_class, n_genes, fill = gene_class)) +
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = "grey94", color = NA) +
  geom_col(width = 0.62, color = "white", linewidth = 0.3) +
  geom_text(aes(label = state_label), hjust = -0.06,
            size = 2.8, fontface = "bold", color = "#222222") +
  annotate("text", x = shared_text_x, y = shared_text_y, label = "Shared",
           size = 2.7, fontface = "bold", color = "#555555") +
  scale_fill_manual(values = gene_class_colors, guide = "none") +
  scale_x_discrete(limits = rev(gene_class_levels)) +
  scale_y_continuous(labels = comma,
                     breaks = seq(0, 6000, by = 2000),
                     limits = c(0, gene_axis_max),
                     expand = expansion(mult = c(0, 0))) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Unique genes", title = "C  Gene-state hierarchy") +
  theme_pub(base = 10.5) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 9.2, face = "bold"),
        axis.text.x = element_text(size = 8.8),
        axis.title.x = element_text(size = 10.5, face = "bold"),
        plot.title = element_text(hjust = 0, face = "bold", size = 13),
        plot.margin = margin(4, 18, 4, 4))

# --- figure 07D: directional coherence against shuffled labels ---

specificity_text_char_width <- 0.018
specificity_null[, max_line_chars := vapply(
  strsplit(plot_label, "\n", fixed = TRUE),
  function(x) max(nchar(x)), integer(1L)
)]
specificity_null[, estimated_label_right :=
                   observed_to_null + max_line_chars * specificity_text_char_width]
ratio_x_max <- ceiling((max(specificity_null$estimated_label_right) + 0.10) * 20) / 20

p_specificity <- ggplot(specificity_null, aes(observed_to_null, metric_label)) +
  geom_vline(xintercept = 1, color = "grey55", linewidth = 0.45,
             linetype = "dashed") +
  geom_errorbar(aes(xmin = null_ratio_q025, xmax = null_ratio_q975),
                width = 0, orientation = "y",
                linewidth = 2.2, color = "grey80") +
  geom_point(aes(color = metric_id), size = 3.2) +
  geom_text(aes(label = plot_label), hjust = -0.18,
            size = 2.65, lineheight = 0.88,
            fontface = "bold", color = "#222222") +
  scale_color_manual(values = specificity_colors, guide = "none") +
  scale_x_continuous(limits = c(0.9, ratio_x_max),
                     breaks = c(1, 1.2, 1.4, 1.6),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Observed / shuffled null", y = NULL,
       title = "D  Gene-level directionality") +
  theme_pub(base = 10.5) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
        axis.text.y = element_text(size = 8.8, face = "bold"),
        axis.text.x = element_text(size = 8.8),
        axis.title.x = element_text(size = 9.8, face = "bold"),
        plot.title = element_text(hjust = 0, face = "bold", size = 13),
        plot.margin = margin(10, 28, 18, 8))

# --- assemble figure 07 ---

design <- "
ABBB
CBBB
DBBB
"
p07 <- p_cluster + p_trajectory + p_gene + p_specificity +
  plot_layout(design = design, widths = c(1.20, 1, 1, 1),
              heights = c(1, 1, 1.05)) &
  theme(plot.title = element_text(hjust = 0),
        plot.margin = margin(8, 12, 10, 10))

ggsave(out_png, p07, width = 11.4, height = 8.6, dpi = 300, bg = "white")
ggsave(out_pdf, p07, width = 11.4, height = 8.6, bg = "white")

# --- save companion summaries ---

stats_summary <- rbindlist(list(
  data.table(
    test_id = "gene_level_directional_coherence",
    scope = "All dynamic APA genes",
    method = sprintf("%s fixed-margin site-label permutations",
                     comma(SPECIFICITY_PERM_ITER)),
    contrast = "Observed C1-specific plus C2-specific genes vs shuffled labels",
    statistic = specificity_null[metric_id == "any_spec", observed_to_null],
    odds_ratio = NA_real_,
    pvalue = specificity_null[metric_id == "any_spec", empirical_p_greater],
    interpretation = sprintf(
      "Observed %s specific genes; shuffled mean %.1f; observed/null %.2f",
      comma(specificity_null[metric_id == "any_spec", observed_count]),
      specificity_null[metric_id == "any_spec", null_mean_count],
      specificity_null[metric_id == "any_spec", observed_to_null]
    )
  ),
  data.table(
    test_id = "c2_relative_specificity",
    scope = "C1-specific and C2-specific gene states",
    method = sprintf("%s fixed-margin site-label permutations",
                     comma(SPECIFICITY_PERM_ITER)),
    contrast = "C2-specific observed/null enrichment relative to C1-specific observed/null enrichment",
    statistic = specificity_null[metric_id == "c2_c1_relative", observed_to_null],
    odds_ratio = NA_real_,
    pvalue = specificity_null[metric_id == "c2_c1_relative", empirical_p_greater],
    interpretation = sprintf(
      "C1-specific observed/null %.2f; C2-specific observed/null %.2f; C2:C1 %.2f",
      specificity_null[metric_id == "c1_spec", observed_to_null],
      specificity_null[metric_id == "c2_spec", observed_to_null],
      specificity_null[metric_id == "c2_c1_relative", observed_to_null]
    )
  ),
  data.table(
    test_id = "expression_direction_by_utr_trajectory",
    scope = "All dynamic APA gene classes with expression-age rho",
    method = "Pearson chi-square test",
    contrast = "Expression up vs down distribution across five 3'UTR trajectory states",
    statistic = unname(expr_chisq$statistic),
    odds_ratio = NA_real_,
    pvalue = expr_chisq$p.value,
    interpretation = sprintf("Expression-up range %.1f%% to %.1f%% across 3'UTR states",
                             100 * min(gene_comp$pct_expression_up, na.rm = TRUE),
                             100 * max(gene_comp$pct_expression_up, na.rm = TRUE))
  ),
  data.table(
    test_id = "directional_coupling_by_c1_c2_side",
    scope = "C1-side and C2-side genes with expression-age rho",
    method = "Fisher exact test",
    contrast = "Coupled vs decoupled expression direction in C1-side vs C2-side genes",
    statistic = NA_real_,
    odds_ratio = directional_odds,
    pvalue = directional_fisher$p.value,
    interpretation = sprintf("Coupled genes: C1-side %.1f%%; C2-side %.1f%%",
                             100 * directional_df[apa_side == "C1 3'UTR lengthening",
                                                   mean(coupling == "Coupled")],
                             100 * directional_df[apa_side == "C2 3'UTR shortening",
                                                   mean(coupling == "Coupled")])
  )
))
stats_summary[, fdr := p.adjust(pvalue, "BH")]

fwrite(cluster_comp[, .(cluster, concordance, membership_state, n_genes,
                        total_cluster_genes, pct)],
       out_cluster, sep = "\t")
fwrite(gene_comp[, .(gene_class, n_genes, pct_all_genes, n_sites,
                     n_c1_sites, n_c2_sites, n_expr_genes,
                     pct_expression_up, median_expr_rho,
                     q25_expr_rho, q75_expr_rho, state_full_name,
                     membership_level)],
       out_gene, sep = "\t")
fwrite(trajectory_summary[, .(gene_class, age_model, grid_age, is_supported,
                              n_genes, median_z, q25_z, q75_z)],
       out_trajectory, sep = "\t")
fwrite(stats_summary, out_stats, sep = "\t")
fwrite(specificity_null, out_specificity, sep = "\t")
fwrite(specificity_perm, out_specificity_perm, sep = "\t")
