# Cross-layer summary panels: RNA and protein amplitude overview, decoupling ranking, APA coefficients and trajectories.
# Writes each panel as PDF and PNG, with layout checks.

library(data.table)
setDTthreads(1L)

# --- config ---

if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
mod <- file.path(Sys.getenv("APA_ROOT"), "data/03_omics")
figure_dir <- file.path(Sys.getenv("APA_ROOT"), "figures/03_omics")
primary_df <- fread(file.path(mod, "02_df_calibration.tsv"))$fixed_df[1]
fixed_model <- sprintf("fixed_ns_df%d", primary_df)
dpi <- 300L
font_size <- 9
layer_colors <- c(APA = "#009E73", RNA = "#0072B2", protein = "#D55E00")
model_colors <- c("Age" = "#777777", "Age + RNA" = "#0072B2")
background_color <- "#BFC4C9"
scatter_size <- c(4.3, 3.8)
ranking_size <- c(4.5, 4.7)
forest_size <- c(5.8, 4.7)
trajectory_size <- c(4.1, 3.3)

# --- helpers ---

panel_legend <- function(...) {
  key <- legend(...)
  horizontal <- grconvertX(
    c(key$rect$left, key$rect$left + key$rect$w), "user", "ndc"
  )
  vertical <- grconvertY(
    c(key$rect$top - key$rect$h, key$rect$top), "user", "ndc"
  )
  all(c(horizontal, vertical) >= 0 & c(horizontal, vertical) <= 1)
}

exported_panels <- character()
save_panel <- function(stem, size, draw) {
  for (extension in c("pdf", "png")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    if (extension == "pdf") {
      pdf(path, width = size[1], height = size[2], family = "Helvetica",
          pointsize = font_size, useDingbats = FALSE)
    } else {
      png(path, width = size[1], height = size[2], units = "in", res = dpi,
          type = "cairo", pointsize = font_size, bg = "white")
    }
    tryCatch({
      par(family = "sans", las = 1, bty = "l", tcl = -0.22,
          mgp = c(2.4, 0.65, 0), cex.axis = 0.9, cex.lab = 0.95)
      checks <- draw()
      cat(sprintf("LAYOUT %s.%s: %s\n", stem, extension,
                  paste(names(checks), ifelse(checks, "PASS", "FAIL"),
                        sep = "=", collapse = "; ")))
      stopifnot(all(checks))
    }, finally = dev.off())
  }
  exported_panels <<- c(exported_panels, stem)
}

# --- loading data ---

rna <- fread(file.path(mod, "02_rna_trajectories.tsv"))[
  analysis_df == primary_df,
  .(gencode_gene_base, gene_symbol, protein_accession,
    rna_amplitude = standardized_amplitude)
]
protein <- fread(file.path(mod, "02_protein_traj.tsv"))[
  analysis_df == primary_df,
  .(gencode_gene_base, protein_accession,
    protein_amplitude = fitted_amplitude / residual_sd)
]
overview <- merge(rna, protein, by = c("gencode_gene_base", "protein_accession"))
ranked <- fread(file.path(mod, "03_decoupling_ranked_pairs.tsv"))
setorder(ranked, decoupling_rank)
ranked[, pair_label := sprintf("P%02d %s", decoupling_rank, gene_symbol)]
ranking_size[2] <- max(ranking_size[2], 2.2 + 0.19 * nrow(ranked))
forest_size[2] <- max(forest_size[2], 2.2 + 0.19 * nrow(ranked))
pairs <- fread(file.path(mod, "03_site_protein_pairs.tsv"))[
  analysis_df == primary_df,
  .(site_id, protein_accession, se_age_adjusted, se_rna_adjusted)
]
coefficients <- merge(ranked, pairs, by = c("site_id", "protein_accession"),
                      all.x = TRUE, sort = FALSE)
setorder(coefficients, decoupling_rank)
sample_qc <- fread(file.path(mod, "01_sample_qc.tsv"))
n_samples <- as.integer(sample_qc[metric == "n_samples", value])
curves <- fread(file.path(mod, "04_gam_curve_data.tsv"))
support <- fread(file.path(mod, "01_support_profile.tsv"))[is_supported == TRUE]
setorder(support, grid_index)
support[, support_segment := cumsum(c(TRUE, diff(grid_index) > 1))]
curves[, support_segment := support$support_segment[match(grid_index, support$grid_index)]]
trajectory_groups <- unique(curves[, .(gene_symbol, layer, feature_id, decoupling_ranks)])
trajectory_groups[, first_rank := vapply(strsplit(decoupling_ranks, ";", fixed = TRUE),
                                         function(x) min(as.integer(x)), integer(1))]
trajectory_groups[, layer_order := match(layer, c("APA", "RNA", "protein"))]
setorder(trajectory_groups, first_rank, layer_order)

expected_trajectory_groups <- nrow(ranked) +
  2L * uniqueN(ranked[, .(gencode_gene_base, protein_accession)])
stopifnot(nrow(overview) > 0L,
          !anyDuplicated(overview$gencode_gene_base),
          all(is.finite(overview$rna_amplitude) & overview$rna_amplitude > 0),
          all(is.finite(overview$protein_amplitude) & overview$protein_amplitude > 0),
          nrow(ranked) > 0L, nrow(coefficients) == nrow(ranked),
          identical(coefficients$site_id, ranked$site_id),
          all(is.finite(coefficients$se_age_adjusted)),
          all(is.finite(coefficients$se_rna_adjusted)),
          all(ranked$q_age_adjusted < 0.05 & ranked$q_rna_adjusted < 0.05),
          !anyNA(curves$support_segment),
          nrow(trajectory_groups) == expected_trajectory_groups)
cat(sprintf("INPUTS: %d genes; %d ranked pairs; %d trajectory panels; %d donors\n",
            nrow(overview), nrow(ranked), nrow(trajectory_groups), n_samples))

# --- RNA and protein amplitude overview ---

overview[, beta_hit := gencode_gene_base %in% ranked$gencode_gene_base]
save_panel("rna_protein_amplitudes", scatter_size, function() {
  par(mar = c(4.1, 4.3, 1.2, 0.8))
  x_range <- range(log10(overview$rna_amplitude))
  y_range <- range(log10(overview$protein_amplitude))
  plot(log10(overview$rna_amplitude), log10(overview$protein_amplitude),
       pch = 16, cex = 0.37, col = adjustcolor(background_color, alpha.f = 0.4),
       xlim = x_range + c(-0.04, 0.05) * diff(x_range),
       ylim = y_range + c(-0.04, 0.26) * diff(y_range),
       xlab = "RNA standardized amplitude", ylab = "Protein standardized amplitude",
       xaxt = "n", yaxt = "n")
  ticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
  x_ticks <- ticks[log10(ticks) >= par("usr")[1] & log10(ticks) <= par("usr")[2]]
  y_ticks <- ticks[log10(ticks) >= par("usr")[3] & log10(ticks) <= par("usr")[4]]
  axis(1, at = log10(x_ticks), labels = x_ticks)
  axis(2, at = log10(y_ticks), labels = y_ticks)
  abline(a = 0, b = 1, col = "#B7B7B7", lty = 2, lwd = 0.8)
  points(log10(overview[beta_hit == TRUE, rna_amplitude]),
         log10(overview[beta_hit == TRUE, protein_amplitude]),
         pch = 21, bg = "#0072B2", col = "white", lwd = 0.5, cex = 0.95)
  key_ok <- panel_legend("topleft", inset = 0.01,
                         legend = c(sprintf("Measured genes (%s)",
                                            format(nrow(overview), big.mark = ",")),
                                    sprintf("Beta-hit genes (%d)",
                                            uniqueN(ranked$gencode_gene_base))),
                         col = c(background_color, "#0072B2"), pch = 16,
                         bty = "n", cex = 0.85, pt.cex = 0.85)
  c(legend_bounds = key_ok,
    x_label_width = strwidth("RNA standardized amplitude", "inches", cex = 0.95) < par("pin")[1],
    y_label_width = strwidth("Protein standardized amplitude", "inches", cex = 0.95) < par("pin")[2])
})

# --- continuous decoupling ranking ---

row_positions <- rev(seq_len(nrow(ranked)))
save_panel("decoupling_ranking", ranking_size, function() {
  par(mar = c(4.4, 8.5, 1.3, 0.8))
  ratio_ticks <- c(0.2, 0.5, 1, 2, 5, 10)
  plot(NA, xlim = range(log2(ratio_ticks)), ylim = c(0.4, nrow(ranked) + 1.5),
       xlab = "RNA / protein amplitude ratio\n(both in residual SD units)", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(h = row_positions, col = "#EFEFEF", lwd = 0.6)
  abline(v = 0, col = "#999999", lty = 2, lwd = 0.8)
  axis(1, at = log2(ratio_ticks), labels = ratio_ticks)
  axis(2, at = row_positions, labels = ranked$pair_label, las = 2, tick = FALSE,
       cex.axis = 0.9)
  segments(0, row_positions, log2(ranked$decoupling_ratio), row_positions,
           col = "#B9CBD8", lwd = 1.2)
  points(log2(ranked$decoupling_ratio), row_positions,
         pch = ifelse(ranked$high_collinearity, 1, 16), col = "#0072B2", cex = 0.85)
  key_ok <- panel_legend("topleft", inset = 0.01, legend = c("VIF < 5", "VIF >= 5"),
                         pch = c(16, 1), col = "#0072B2", bty = "n", horiz = TRUE,
                         cex = 0.85)
  c(legend_bounds = key_ok,
    row_label_margin = max(strwidth(ranked$pair_label, "inches", cex = 0.9)) + 0.08 < par("mai")[2],
    row_spacing = par("pin")[2] / diff(par("usr")[3:4]) > strheight("Mg", "inches", cex = 0.9) * 1.5)
})

# --- APA coefficients before and after RNA adjustment ---

critical_m0 <- qt(0.975, df = n_samples - primary_df - 2L)
critical_m1 <- qt(0.975, df = n_samples - primary_df - 3L)
coefficients[, `:=`(
  lower_m0 = beta_age_adjusted - critical_m0 * se_age_adjusted,
  upper_m0 = beta_age_adjusted + critical_m0 * se_age_adjusted,
  lower_m1 = beta_rna_adjusted - critical_m1 * se_rna_adjusted,
  upper_m1 = beta_rna_adjusted + critical_m1 * se_rna_adjusted
)]
save_panel("apa_protein_coefficients", forest_size, function() {
  par(mar = c(4.4, 8.5, 1.3, 8.7))
  x_range <- range(coefficients[, c(lower_m0, upper_m0, lower_m1, upper_m1)])
  plot(NA, xlim = x_range + c(-0.05, 0.05) * diff(x_range),
       ylim = c(0.4, nrow(ranked) + 1.5),
       xlab = "APA coefficient\n(log2 protein per log2 APA odds)", ylab = "", yaxt = "n")
  abline(h = row_positions, col = "#EFEFEF", lwd = 0.6)
  abline(v = 0, col = "#999999", lty = 2, lwd = 0.8)
  axis(2, at = row_positions, labels = ranked$pair_label, las = 2, tick = FALSE,
       cex.axis = 0.9)
  segments(coefficients$lower_m0, row_positions + 0.13,
           coefficients$upper_m0, row_positions + 0.13, col = model_colors[1], lwd = 1.1)
  segments(coefficients$lower_m1, row_positions - 0.13,
           coefficients$upper_m1, row_positions - 0.13, col = model_colors[2], lwd = 1.1)
  points(coefficients$beta_age_adjusted, row_positions + 0.13,
         col = model_colors[1], pch = 16, cex = 0.7)
  points(coefficients$beta_rna_adjusted, row_positions - 0.13,
         col = model_colors[2], pch = 16, cex = 0.7)
  q_x <- par("usr")[2] + diff(par("usr")[1:2]) * 0.07
  q_labels <- sprintf("%.3f / %.3f", coefficients$q_age_adjusted,
                      coefficients$q_rna_adjusted)
  text(q_x, row_positions, q_labels, adj = c(0, 0.5), cex = 0.82, xpd = NA)
  text(q_x, nrow(ranked) + 1.05, "Global BH q", adj = c(0, 0.5), cex = 0.85, xpd = NA)
  text(q_x, nrow(ranked) + 0.55, "Age / Age + RNA", adj = c(0, 0.5), cex = 0.82, xpd = NA)
  key_ok <- panel_legend("topleft", inset = 0.01, legend = names(model_colors),
                         col = model_colors, pch = 16, lty = 1, bty = "n", horiz = TRUE,
                         cex = 0.85, seg.len = 1.1)
  q_right <- grconvertX(q_x, "user", "inches") +
    max(strwidth(c(q_labels, "Age / Age + RNA"), "inches", cex = 0.82))
  c(legend_bounds = key_ok,
    q_column_bounds = q_right < par("din")[1] - 0.05,
    row_label_margin = max(strwidth(ranked$pair_label, "inches", cex = 0.9)) + 0.08 < par("mai")[2],
    row_spacing = par("pin")[2] / diff(par("usr")[3:4]) > strheight("Mg", "inches", cex = 0.9) * 1.5)
})

# --- individual gene, layer and site trajectories ---

axis_labels <- c(APA = "Adjusted log2 odds", RNA = "Adjusted log2 CPM",
                 protein = "Adjusted log2 abundance")
for (index in seq_len(nrow(trajectory_groups))) {
  item <- trajectory_groups[index]
  dat <- curves[gene_symbol == item$gene_symbol & layer == item$layer &
                  feature_id == item$feature_id]
  site_suffix <- if (item$layer == "APA") sprintf("_P%02d", item$first_rank) else ""
  stem <- paste0("trajectory_", item$gene_symbol, "_", tolower(item$layer), site_suffix)
  title <- paste(item$gene_symbol, item$layer, sep = " | ")
  if (item$layer == "APA") title <- paste(title, sprintf("P%02d", item$first_rank), sep = " | ")
  color <- layer_colors[[item$layer]]
  y_label <- axis_labels[[item$layer]]
  save_panel(stem, trajectory_size, function() {
    fixed <- dat[model == fixed_model]
    adaptive <- dat[model == "adaptive_gam"]
    y_range <- range(fixed$lower_95, fixed$upper_95, adaptive$fitted)
    y_span <- diff(y_range)
    age_ticks <- c(-0.5, 0, 1, 12, 40, 100)
    par(mar = c(4.1, 4.3, 1.8, 0.8))
    plot(NA, xlim = range(support$age_model),
         ylim = y_range + c(-0.06, 0.27) * y_span,
         xlab = "Age (years; log10(age + 1))", ylab = y_label,
         xaxt = "n")
    axis(1, at = log10(age_ticks + 1), labels = age_ticks, cex.axis = 0.9)
    mtext(title, side = 3, line = 0.55, adj = 0, cex = 1, font = 2)
    for (segment_id in unique(fixed$support_segment)) {
      band <- fixed[support_segment == segment_id][order(grid_index)]
      comparison <- adaptive[support_segment == segment_id][order(grid_index)]
      polygon(c(band$age_model, rev(band$age_model)),
              c(band$lower_95, rev(band$upper_95)),
              col = adjustcolor(color, alpha.f = 0.18), border = NA)
      lines(band$age_model, band$fitted, col = color, lwd = 1.5)
      lines(comparison$age_model, comparison$fitted,
            col = "#666666", lwd = 1.1, lty = 2)
    }
    key_ok <- panel_legend("topleft", inset = 0.015,
                           legend = c("Fixed spline (95% CI)", "Adaptive GAM"),
                           col = c(color, "#666666"), lty = c(1, 2),
                           lwd = c(1.5, 1.1), horiz = TRUE, bty = "n", cex = 0.82,
                           seg.len = 1.5, x.intersp = 0.6)
    label_width <- strwidth(as.character(age_ticks), units = "inches", cex = 0.9)
    tick_positions <- grconvertX(log10(age_ticks + 1), "user", "inches")
    c(legend_bounds = key_ok,
      age_tick_spacing = all(diff(tick_positions) >
                              (head(label_width, -1) + tail(label_width, -1)) / 2),
      axis_label_width = strwidth("Age (years; log10(age + 1))", "inches",
                                 cex = 0.95) < par("pin")[1],
      title_width = strwidth(title, "inches") < par("pin")[1],
      supported_segments = all(vapply(split(fixed$grid_index, fixed$support_segment),
                                       function(x) all(diff(sort(x)) == 1), logical(1))))
  })
}

# --- export summary ---

output_paths <- file.path(figure_dir, paste0(rep(exported_panels, each = 2),
                                           rep(c(".pdf", ".png"), times = length(exported_panels))))
stopifnot(length(exported_panels) == 3L + expected_trajectory_groups,
          !anyDuplicated(output_paths),
          all(file.info(output_paths)$size > 0))
cat(sprintf("COMPLETE: %d standalone panels; %d PDF/PNG files; all layout checks passed\n",
            length(exported_panels), length(output_paths)))
