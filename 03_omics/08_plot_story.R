#!/usr/bin/env Rscript
# Story figures: trajectory triplets, nested delta R2 with VIF, INPP4A motif track,
# and protein robustness with repeated donor CV.

suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
omics_dir <- file.path(apa_root, "data", "03_omics")
figure_dir <- file.path(apa_root, "figures", "03_omics")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
primary_df <- fread(file.path(omics_dir, "02_df_calibration.tsv"))$fixed_df[1]
fixed_model <- sprintf("fixed_ns_df%d", primary_df)
dpi <- 300L
font_size <- 9
layer_colors <- c(APA = "#009E73", RNA = "#0072B2", protein = "#D55E00")
candidate_color <- "#0072B2"
triplet_genes <- data.table(
  gene_symbol = c("PIK3CA", "KPNA4", "INPP4A", "ARIH1", "RHOT1", "NAPRT"),
  site_id = c(
    "ENST00000675467.1|PIK3CA|chr3:179234094-179235371|+",
    "ENST00000334256.9|KPNA4|chr3:160495007-160502202|-",
    "ENST00000463367.1|INPP4A|chr2:98520687-98521766|+",
    "ENST00000563310.1|ARIH1|chr15:72583208-72586559|+",
    "ENST00000354266.7|RHOT1|chr17:32224616-32225666|+",
    "ENST00000340490.7|NAPRT|chr8:143574787-143575093|-"
  )
)
triplet_size <- c(7.6, 5.0)
delta_r2_size <- c(7.0, 4.9)
motif_size <- c(7.4, 4.2)
sensitivity_files <- c(
  exclude_PTB752 = "04_robustness_exclude_PTB752.tsv.gz",
  imputed_observed_covariates = "04_robustness_imputed_observed_covariates.tsv.gz",
  all_unambiguous_proteins = "04_robustness_all_unambiguous_proteins.tsv.gz",
  prenatal_continuous_age = "04_robustness_prenatal_continuous_age.tsv.gz",
  postnatal_continuous_age = "04_robustness_postnatal_continuous_age.tsv.gz",
  post_combat_batch = "04_robustness_post_combat_batch.tsv.gz",
  observed_only = "04_robust_obs_only.tsv.gz"
)

# --- helpers ---
save_panel <- function(stem, size, draw, point_size = font_size) {
  for (extension in c("pdf", "png")) {
    path <- file.path(figure_dir, paste0(stem, ".", extension))
    if (extension == "pdf") {
      pdf(path, width = size[1], height = size[2], family = "Helvetica",
          pointsize = point_size, useDingbats = FALSE)
    } else {
      png(path, width = size[1], height = size[2], units = "in", res = dpi,
          type = "cairo", pointsize = point_size, bg = "white")
    }
    par(family = "sans", las = 1, bty = "l", tcl = -0.22,
        mgp = c(2.4, 0.65, 0), cex.axis = 0.9, cex.lab = 0.95)
    draw()
    dev.off()
  }
}

format_q <- function(q) {
  if (q < 1e-3) sprintf("%.1e", q) else sprintf("%.3f", q)
}

# --- load data ---
ranked <- fread(file.path(omics_dir, "03_decoupling_ranked_pairs.tsv"))
setorder(ranked, decoupling_rank)
curves <- fread(file.path(omics_dir, "04_gam_curve_data.tsv"))
naprt_control_site <- triplet_genes[gene_symbol == "NAPRT", site_id]
control_pairs <- fread(
  file.path(omics_dir, "03_site_protein_pairs.tsv"),
  select = c("site_id", "gene_symbol", "gencode_gene_base", "protein_accession",
             "analysis_df", "q_rna_adjusted", "q_rna_independent")
)[analysis_df == primary_df & site_id == naprt_control_site]
control_curves <- fread(file.path(omics_dir, "07_naprt_ctl_curves.tsv"))
curves <- rbindlist(list(curves, control_curves), use.names = TRUE, fill = TRUE)
support <- fread(file.path(omics_dir, "01_support_profile.tsv"))[is_supported == TRUE]
setorder(support, grid_index)
support[, support_segment := cumsum(c(TRUE, diff(grid_index) > 1))]
curves[, support_segment := support$support_segment[match(grid_index, support$grid_index)]]
model_comparison <- fread(file.path(omics_dir, "07_model_comp_delta_r2.tsv"))
motif_regions <- fread(file.path(omics_dir, "07_motif_track_regions.tsv"))
motif_sites <- fread(file.path(omics_dir, "07_motif_track_sites.tsv"))
mirna_rank <- fread(file.path(omics_dir, "07_inpp4a_mirna_brain_rank.tsv"))

stopifnot(nrow(control_pairs) == 1L, !anyNA(curves$support_segment))

# --- figure 1: apa / rna / protein trajectory triplets ---
# centered original-unit curves on a compressed age axis: the unsupported
# 0.16-11.6 year interval is squeezed to a narrow band and bridged with dashed
# lines plus a translucent interpolated band
gap_start <- support[support_segment == 1L, max(age_model)]
gap_end <- support[support_segment == 2L, min(age_model)]
gap_compress <- 0.12
compress_age <- function(x) {
  ifelse(x <= gap_start, x,
         ifelse(x >= gap_end,
                gap_start + (gap_end - gap_start) * gap_compress + (x - gap_end),
                gap_start + (x - gap_start) * gap_compress))
}
triplet_ticks <- c(-0.5, 0, 12, 40, 100)

save_panel("trajectory_triplets", triplet_size, function() {
  layout(matrix(c(1L, 1L, 1L, 2L, 3L, 4L, 5L, 6L, 7L), nrow = 3L,
                byrow = TRUE),
         heights = c(0.11, 1, 1))
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = names(layer_colors), col = layer_colors,
         lty = 1, lwd = 2, horiz = TRUE, bty = "n", cex = 1,
         seg.len = 1.6, x.intersp = 0.7)
  for (index in seq_len(nrow(triplet_genes))) {
    gene_symbol <- triplet_genes$gene_symbol[index]
    selected_site_id <- triplet_genes$site_id[index]
    show_xlab <- index > 3L

    pair_row <- ranked[site_id == selected_site_id]
    if (nrow(pair_row) == 0L) {
      pair_row <- control_pairs[site_id == selected_site_id]
    }
    features <- c(APA = pair_row$site_id,
                  RNA = pair_row$gencode_gene_base,
                  protein = pair_row$protein_accession)
    layer_data <- lapply(names(features), function(layer_name) {
      feature_value <- features[[layer_name]]
      gene_sym <- gene_symbol
      dat <- curves[gene_symbol == gene_sym & layer == layer_name &
                      feature_id == feature_value & model == fixed_model]
      setorder(dat, grid_index)
      center <- mean(dat$fitted)
      dat[, `:=`(y = fitted - center,
                 y_low = lower_95 - center,
                 y_high = upper_95 - center,
                 x = compress_age(age_model))]
      dat
    })
    names(layer_data) <- names(features)
    y_values <- unlist(lapply(layer_data, function(d) d[, c(y, y_low, y_high)]))
    y_range <- range(y_values)
    y_span <- diff(y_range)

    par(mar = c(3.0, 4.0, 1.7, 0.6))
    plot(NA, xlim = range(compress_age(support$age_model)),
         ylim = y_range + c(-0.08, 0.36) * y_span,
         xlab = "", ylab = "Deviation from mean", xaxt = "n")
    axis(1, at = compress_age(log10(triplet_ticks + 1)), labels = triplet_ticks,
         cex.axis = 0.92)
    if (show_xlab) {
      mtext("Age (years)", side = 1, line = 1.9, cex = 0.8)
    }

    band_left <- compress_age(gap_start)
    band_right <- compress_age(gap_end)
    rect(band_left, par("usr")[3], band_right, par("usr")[4],
         col = adjustcolor("#888888", alpha.f = 0.10), border = NA)

    for (layer in names(layer_data)) {
      dat <- layer_data[[layer]]
      color <- layer_colors[[layer]]
      for (segment_id in unique(dat$support_segment)) {
        band <- dat[support_segment == segment_id][order(grid_index)]
        polygon(c(band$x, rev(band$x)), c(band$y_low, rev(band$y_high)),
                col = adjustcolor(color, alpha.f = 0.16), border = NA)
      }
      seg1 <- dat[support_segment == 1L][.N]
      seg2 <- dat[support_segment == 2L][1L]
      if (!is.na(seg1$grid_index) && !is.na(seg2$grid_index)) {
        polygon(c(seg1$x, seg2$x, seg2$x, seg1$x),
                c(seg1$y_low, seg2$y_low, seg2$y_high, seg1$y_high),
                col = adjustcolor(color, alpha.f = 0.10), border = NA)
      }
      for (segment_id in unique(dat$support_segment)) {
        band <- dat[support_segment == segment_id][order(grid_index)]
        lines(band$x, band$y, col = color, lwd = 1.6)
      }
      if (!is.na(seg1$grid_index) && !is.na(seg2$grid_index)) {
        lines(c(seg1$x, seg2$x), c(seg1$y, seg2$y), col = color, lwd = 1.3,
              lty = 3)
      }
    }

    mtext(gene_symbol, side = 3, line = 0.4, adj = 0, cex = 1.1, font = 2)
    stats_lines <- c(sprintf("APA q = %s", format_q(pair_row$q_rna_adjusted)),
                     sprintf("RNA q = %s", format_q(pair_row$q_rna_independent)))
    text(par("usr")[2], par("usr")[4] - 0.03 * y_span, stats_lines[1],
         adj = c(1, 1), cex = 0.95)
    text(par("usr")[2], par("usr")[4] - 0.16 * y_span, stats_lines[2],
         adj = c(1, 1), cex = 0.95)
  }
}, point_size = 12)

# --- figure 2: nested-model delta R2 with partial F-test and VIF ---
gene_comparison <- merge(
  model_comparison[site_id %chin% ranked$site_id],
  ranked[, .(site_id, high_collinearity, high_vif_source,
             q_rna_adjusted, q_rna_independent)],
  by = "site_id"
)
naprt_comparison <- merge(
  model_comparison[site_id == naprt_control_site],
  control_pairs[, .(site_id, q_rna_adjusted, q_rna_independent)],
  by = "site_id"
)
gene_comparison <- gene_comparison[
  order(decoupling_rank),
  .SD[1L],
  by = gene_symbol
]
stopifnot(nrow(naprt_comparison) == 1L)
gene_comparison <- gene_comparison[gene_symbol != "NAPRT"]
# VIF-clean genes take the visual lead; high-VIF genes go to the right as caveat
setorder(gene_comparison, high_collinearity, -delta_r2)
n_genes <- nrow(gene_comparison)
n_vif_high <- sum(gene_comparison$high_collinearity)
high_vif_color <- "#BFBFBF"
delta_r2_size[1] <- max(delta_r2_size[1], 2.4 + 0.38 * n_genes)

save_panel("model_comparison_delta_r2", delta_r2_size, function() {
  par(oma = c(0.4, 0.3, 0.2, 0.2))
  layout(matrix(c(1L, 2L, 3L, 3L), nrow = 2L),
         widths = c(3.6, 1.1), heights = c(3.1, 1.3))

  bar_space <- rep(0.3, n_genes)
  if (n_vif_high > 0L) bar_space[n_genes - n_vif_high + 1L] <- 0.9
  bar_left <- cumsum(bar_space) + seq_len(n_genes) - 1L
  bar_positions <- bar_left + 0.5
  plot_xlim <- c(0, bar_left[n_genes] + 1.3)
  bar_fill <- ifelse(gene_comparison$high_collinearity,
                     high_vif_color, candidate_color)

  # main panel: candidate genes only
  par(mar = c(0.4, 4.6, 1.6, 0.3))
  plot(NA, xlim = plot_xlim,
       ylim = c(0, max(gene_comparison$delta_r2) * 1.22),
       ylab = expression(Delta * R^2 * " (APA | RNA, age)"),
       xlab = "", xaxt = "n", cex.lab = 1.05)
  rect(bar_left, 0, bar_left + 1, gene_comparison$delta_r2,
       col = bar_fill, border = NA)
  text(bar_positions, gene_comparison$delta_r2,
       sprintf("%.3f", gene_comparison$delta_r2),
       pos = 3, cex = 0.85, offset = 0.25)
  legend("topleft", inset = 0.015, legend = c("VIF clean", "VIF high"),
         fill = c(candidate_color, high_vif_color), bty = "n", cex = 0.95,
         horiz = TRUE, x.intersp = 0.6)

  # annotation track
  par(mar = c(4.6, 4.6, 0.3, 0.3))
  plot(NA, xlim = plot_xlim, ylim = c(0.5, 2.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(2, at = c(2, 1), labels = c("F-test", "VIF"), las = 1, tick = FALSE,
       cex.axis = 0.95)
  text(bar_positions, rep(2, n_genes),
       vapply(gene_comparison$lrt_p, function(p) {
         if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
       }, character(1)), cex = 1.05)
  points(bar_positions, rep(1, n_genes), pch = 22, bg = bar_fill,
         col = bar_fill, cex = 1.05)
  text(bar_positions, 0.34, gene_comparison$gene_symbol, srt = 45,
       adj = c(1, 1), cex = 0.92, xpd = NA)

  # side panel: NAPRT coefficients in their original predictor units
  par(mar = c(4.6, 3.4, 1.6, 0.3))
  beta_bars <- c(abs(naprt_comparison$beta_rna_full),
                 abs(naprt_comparison$beta_apa_full))
  beta_q <- c(naprt_comparison$q_rna_independent,
              naprt_comparison$q_rna_adjusted)
  plot(NA, xlim = c(0.2, 3.1), ylim = c(0, 1),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(2, at = c(0, 0.5, 1), cex.axis = 0.9)
  mtext(expression("|" * beta * "|"), side = 2, line = 2.0,
        cex = 0.85, las = 0)
  mtext("NAPRT", side = 3, line = 0.4, cex = 0.9, font = 2, adj = 1)
  rect(c(0.6, 1.9), 0, c(1.4, 2.7), beta_bars,
       col = c(layer_colors[["RNA"]], layer_colors[["APA"]]), border = NA)
  text(c(1, 2.3), beta_bars,
       sprintf("%.2f\nq = %s", beta_bars, vapply(beta_q, format_q, character(1))),
       pos = 3, cex = 0.8, offset = 0.4)
  axis(1, at = c(1, 2.3), labels = c("RNA", "APA"), tick = FALSE,
       cex.axis = 0.95)
}, point_size = 12.5)

# --- figure 3: 3'UTR motif track mechanism ---
# display the same INPP4A site across reranked families
display_region <- motif_regions[gene_symbol == "INPP4A"]
display_sites <- motif_sites[gene_symbol == "INPP4A"]
brain_mirnas <- mirna_rank[selected_for_label == TRUE][
  order(brain_rank), .(label, display_label, brain_rank)
]
rug_sites <- display_sites[class != "miRNA_7mer_A1"]
labeled_sites <- merge(
  display_sites[class != "miRNA_7mer_A1"], brain_mirnas, by = "label"
)
labeled_sites[, tier := match(class, c("miRNA_8mer", "miRNA_7mer_m8"))]
setorder(labeled_sites, label, tier, pos_tx)
labeled_sites <- labeled_sites[, .SD[1L], by = label]
setorder(labeled_sites, brain_rank, pos_tx)
stopifnot(nrow(labeled_sites) == nrow(brain_mirnas))

save_panel("motif_track_mechanism", motif_size, function() {
  len_kb <- display_region$locus_length_tx / 1000
  x_pad <- 0.06 * len_kb
  xlim <- c(-x_pad, len_kb + x_pad)
  lengthening <- display_region$apa_direction == "lengthening"
  site_colors <- c(miRNA_8mer = "#C0392B", miRNA_7mer_m8 = "#D9822B")

  par(mar = c(3.6, 1.2, 0.6, 1.2))
  plot(NA, xlim = xlim, ylim = c(-2.0, 3.4), xaxt = "n", yaxt = "n",
       xlab = sprintf("%s tandem 3'UTR (kb, transcript orientation)",
                      display_region$gene_symbol), ylab = "", cex.lab = 1)
  kb_ticks <- pretty(c(0, len_kb), n = 6)
  axis(1, at = kb_ticks, labels = kb_ticks, cex.axis = 0.95)

  # gained / lost region shading
  region_fill <- if (lengthening) "#009E73" else "#D55E00"
  rect(0, -0.32, len_kb, 0.32,
       col = adjustcolor(region_fill, alpha.f = 0.10), border = NA)
  segments(-x_pad, 0, len_kb + x_pad, 0, col = "#333333", lwd = 1.6)
  region_label <- if (lengthening) {
    sprintf("gained region (n = %d matches)", nrow(rug_sites))
  } else {
    sprintf("lost region (n = %d matches)", nrow(rug_sites))
  }
  text(len_kb / 2, -0.72, region_label, cex = 0.78, col = region_fill)

  # proximal / distal locus ends; label as PAS only when PolyA_DB anchored
  pas_x <- c(0, len_kb)
  end_sources <- c(display_region$proximal_pa_source,
                   display_region$distal_pa_source)
  polya_anchored <- end_sources == "PolyA_DB"
  points(pas_x, c(0, 0), pch = ifelse(polya_anchored, 17, 3),
         cex = 1.6, col = "#333333")
  pas_labels <- ifelse(
    polya_anchored,
    c(paste("proximal PAS", display_region$proximal_pas_signal, sep = "\n"),
      paste("distal PAS", display_region$distal_pas_signal, sep = "\n")),
    c("proximal boundary", "distal boundary")
  )
  text(pas_x[1], -1.18, pas_labels[1], cex = 0.9, adj = 0)
  text(pas_x[2], -1.18, pas_labels[2], cex = 0.9, adj = 1)
  pas_note <- ifelse(polya_anchored, "PolyA_DB", "locus boundary")
  text(pas_x[1], -1.82, pas_note[1], cex = 0.78, col = "#777777", adj = 0)
  text(pas_x[2], -1.82, pas_note[2], cex = 0.78, col = "#777777", adj = 1)

  # lollipops and labels for atlas-selected seed families
  n_lab <- nrow(labeled_sites)
  x_kb <- labeled_sites$pos_tx / 1000
  lab_ok <- rep(TRUE, n_lab)
  kb_per_inch <- diff(xlim) / par("pin")[1]
  label_w <- strwidth(labeled_sites$display_label, "inches", cex = 0.85) * kb_per_inch
  label_x <- pmin(pmax(x_kb, label_w / 2 + 0.05), len_kb - label_w / 2 - 0.05)
  label_levels <- integer(n_lab)
  level_right <- numeric(0)
  for (i in which(lab_ok)) {
    left_edge <- label_x[i] - label_w[i] / 2 - 0.1
    free <- which(level_right < left_edge)
    level <- if (length(free) == 0L) length(level_right) + 1L else free[1L]
    if (level > length(level_right)) level_right <- c(level_right, -Inf)
    label_levels[i] <- level
    level_right[level] <- label_x[i] + label_w[i] / 2 + 0.1
  }
  heights <- ifelse(lab_ok, 0.85 + label_levels * 0.55, 0.55)
  segments(x_kb, 0.32, x_kb, heights, col = "#999999", lwd = 0.8)
  points(x_kb, heights, pch = 19, col = site_colors[labeled_sites$class],
         cex = ifelse(lab_ok, 1.1, 0.85))
  text(label_x[lab_ok], heights[lab_ok] + 0.22, labeled_sites$display_label[lab_ok],
       cex = 0.85, col = site_colors[labeled_sites$class[lab_ok]])

  legend("topright", inset = 0.015,
         legend = c("miRNA 8mer", "miRNA 7mer-m8"),
         col = site_colors, pch = 19, bty = "n", cex = 0.9,
         pt.cex = 1)
}, point_size = 12)

# --- figure 4 data: coefficient sensitivity and internal prediction ---
sensitivity <- sapply(names(sensitivity_files), function(name) {
  x <- fread(file.path(omics_dir, sensitivity_files[[name]]),
             select = c("site_id", "q_age_adjusted", "q_rna_adjusted", "beta_rna_adjusted"))
  x <- x[match(ranked$site_id, site_id)]
  stopifnot(identical(x$site_id, ranked$site_id))
  ifelse(sign(x$beta_rna_adjusted) != sign(ranked$beta_rna_adjusted), 2L,
         as.integer(x$q_age_adjusted < .05 & x$q_rna_adjusted < .05))
})
bootstrap <- fread(file.path(omics_dir, "04_donor_boot_ranked.tsv"))[
  match(ranked$site_id, site_id)]
cv <- fread(file.path(omics_dir, "04_cv_ranked_pairs.tsv"))[
  match(ranked$site_id, site_id)]
stopifnot(!anyNA(sensitivity), identical(bootstrap$site_id, ranked$site_id),
          identical(cv$site_id, ranked$site_id))
pair_ids <- sprintf("P%02d", ranked$decoupling_rank)
row_y <- rev(seq_len(nrow(ranked)))
row_top <- nrow(ranked) + 1.8
robustness_height <- max(5.1, 2.3 + 0.19 * nrow(ranked))
sensitivity_colors <- c("#E6E6E6", "#0072B2", "#D55E00")

# --- figure 4: full-family refitting and candidate-conditional prediction ---
# CV figure intervals show repeat-split spread, not confidence intervals.
# Pair IDs follow the current decoupling rank; S1-S7 follow sensitivity_names.
save_panel("protein_robustness_cv", c(9.4, robustness_height), function() {
  par(mfrow = c(1, 3), mar = c(3.8, 3.0, 2.5, 0.7), mgp = c(2.2, .6, 0))
  plot(NA, xlim = c(.5, 7.5), ylim = c(.5, row_top), xlab = "Sensitivity analysis",
       ylab = "", xaxt = "n", yaxt = "n", bty = "n")
  for (j in seq_len(ncol(sensitivity))) {
    rect(j-.42, row_y-.37, j+.42, row_y+.37,
         col = sensitivity_colors[sensitivity[, j]+1L], border = NA)
  }
  axis(1, at = seq_len(7), labels = paste0("S", seq_len(7)), tick = FALSE)
  axis(2, at = row_y, labels = pair_ids, tick = FALSE)
  mtext("A  Coefficient sensitivity", side = 3, line = .7, adj = 0, font = 2)
  legend("top", inset = .005, legend = c("Both q < 0.05", "Below gate"),
         fill = sensitivity_colors[c(2, 1)], bty = "n", cex = .78, horiz = TRUE, x.intersp = .5)

  plot(NA, xlim = c(-.035, 1.04), ylim = c(.5, row_top), xlab = "Bootstrap selection frequency",
       ylab = "", yaxt = "n", xaxs = "i", bty = "l")
  segments(0, row_y, 1, row_y, col = "#EEEEEE")
  points(bootstrap$beta_selection_frequency, row_y+.13, pch = 16, cex = .85, col = "#0072B2")
  points(bootstrap$hard_cut_frequency, row_y-.13, pch = 17, cex = .85, col = "#D55E00")
  mtext("B  Conditional donor bootstrap", side = 3, line = .7, adj = 0, font = 2)
  legend("top", inset = .005, legend = c("Beta gate", "Auxiliary cut"),
         col = c("#0072B2", "#D55E00"), pch = c(16, 17), bty = "n", cex = .78,
         horiz = TRUE, x.intersp = .5)

  cv_range <- range(c(cv$split_low95, cv$split_high95, cv$mean_delta_r2_training_selected))*100
  plot(NA, xlim = cv_range+c(-.3, .3), ylim = c(.5, row_top),
       xlab = expression(Delta*R^2*" (percentage points)"), ylab = "", yaxt = "n")
  abline(v = 0, lty = 3, col = "#777777")
  segments(cv$split_low95*100, row_y+.13, cv$split_high95*100, row_y+.13,
           col = "#0072B2", lwd = 1.1)
  points(cv$mean_delta_r2*100, row_y+.13, pch = 16, cex = .85, col = "#0072B2")
  points(cv$mean_delta_r2_training_selected*100, row_y-.13,
         pch = 17, cex = .85, col = "#D55E00")
  mtext("C  Repeated donor CV", side = 3, line = .7, adj = 0, font = 2)
  legend("top", inset = .005, legend = c("Frozen pair", "Training selected"),
         col = c("#0072B2", "#D55E00"), pch = c(16, 17), bty = "n", cex = .75,
         horiz = TRUE, x.intersp = .45)
}, point_size = 10)
