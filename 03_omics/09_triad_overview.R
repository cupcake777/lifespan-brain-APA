#!/usr/bin/env Rscript
# ATP2B1 triad overview: aligned donor tracks for APA, RNA and protein, with raw observations.

suppressPackageStartupMessages(library(data.table))

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
omics_dir <- file.path(apa_root, "data", "03_omics")
figure_dir <- file.path(apa_root, "figures", "03_omics")
apa_path <- file.path(omics_dir, "apa_matched.rds")
rna_path <- file.path(omics_dir, "rna_matched.rds")
protein_path <- file.path(omics_dir, "protein_matched.rds")
manifest_path <- file.path(omics_dir, "matched_sample_manifest.tsv")
support_path <- file.path(omics_dir, "01_support_profile.tsv")
summary_path <- file.path(omics_dir, "09_atp2b1_triad_sum_results2.tsv")
validation_path <- file.path(omics_dir, "09_global_triad_layout_validation_results3.tsv")
figure_path <- file.path(figure_dir, "09_global_triad_overview")

apa_event <- "ENST00000551310.1|ATP2B1|chr12:89655583-89656107|-"
rna_gene <- "ENSG00000070961"
protein_accession <- "P20020"
display_limit <- 2.5
expected_donor_count <- 371L
expected_sparse_donor_count <- 4L

# --- helpers ---
age_model <- function(x) log10(x + 1)

draw_strip <- function(values, colors, label, show_age_axis = FALSE, ages = NULL) {
  n <- length(values)
  plot.new()
  plot.window(xlim = c(0.5, n + 0.5), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  rect(seq_len(n) - 0.5, 0, seq_len(n) + 0.5, 1,
       col = colors[values], border = NA)
  box(col = "#666666", lwd = 0.6)
  mtext(label, side = 2, line = 1.2, las = 1, cex = 0.76)
  if (show_age_axis) {
    tick_ages <- c(-0.5, 0, 50, 75, 100)
    tick_pos <- vapply(tick_ages, function(x) which.min(abs(ages - x)), integer(1))
    axis(1, at = tick_pos,
         labels = c("Prenatal", "Birth", "50", "75", "100"),
         cex.axis = 0.65, tck = -0.03, mgp = c(1.2, 0.25, 0))
    mtext("Donors ordered by age (years)", side = 1, line = 1.25, cex = 0.72)
  }
}

draw_scatter <- function(x, y, color, ylab, gap, show_x_axis = FALSE) {
  x_range <- range(x)
  y_pad <- diff(range(y)) * 0.06
  plot.new()
  plot.window(xlim = x_range, ylim = range(y) + c(-y_pad, y_pad), xaxs = "i")
  rect(age_model(gap[1]), par("usr")[3], age_model(gap[2]), par("usr")[4],
       col = "#F2F2F2", border = NA)
  points(x, y, pch = 16, cex = 0.42, col = adjustcolor(color, alpha.f = 0.40))
  fitted_segments <- list()
  for (keep in list(x <= age_model(gap[1]), x >= age_model(gap[2]))) {
    if (sum(keep) >= 10L) {
      fit_data <- data.frame(x = x[keep], y = y[keep])
      fit <- loess(y ~ x, data = fit_data, span = 0.65, degree = 1)
      x_new <- seq(min(x[keep]), max(x[keep]), length.out = 150)
      y_new <- predict(fit, newdata = data.frame(x = x_new))
      lines(x_new, y_new, col = color, lwd = 1.4)
      fitted_segments[[length(fitted_segments) + 1L]] <- list(x = x_new, y = y_new)
    }
  }
  if (length(fitted_segments) == 2L) {
    lines(
      c(tail(fitted_segments[[1]]$x, 1), head(fitted_segments[[2]]$x, 1)),
      c(tail(fitted_segments[[1]]$y, 1), head(fitted_segments[[2]]$y, 1)),
      col = color, lwd = 1.3, lty = 2
    )
  }
  axis(2, las = 1, cex.axis = 0.64, tck = -0.03, mgp = c(1.2, 0.28, 0))
  mtext(ylab, side = 2, line = 1.65, cex = 0.7)
  if (show_x_axis) {
    tick_ages <- c(-0.5, 0, 1, 10, 25, 50, 75, 100)
    axis(1, at = age_model(tick_ages), labels = tick_ages,
         cex.axis = 0.64, tck = -0.03, mgp = c(1.2, 0.25, 0))
    mtext("Age (years; log-scaled)", side = 1, line = 1.25, cex = 0.72)
  }
  box(col = "#666666", lwd = 0.6)
}

draw_figure <- function() {
  layout(matrix(seq_len(8), nrow = 4, byrow = TRUE),
         widths = c(3.9, 2.6), heights = c(0.95, 1, 1, 1))
  par(family = "sans", mar = c(0.4, 4.2, 1.2, 0.5), oma = c(1.0, 0.5, 1.0, 0.4))

  plot.new()
  plot.window(xlim = c(0.5, length(donor_ids) + 0.5), ylim = c(0, 2),
              xaxs = "i", yaxs = "i")
  rect(seq_along(donor_ids) - 0.5, 1.05, seq_along(donor_ids) + 0.5, 1.85,
       col = stage_colors[stage_index], border = NA)
  rect(seq_along(donor_ids) - 0.5, 0.10, seq_along(donor_ids) + 0.5, 0.90,
       col = age_palette[age_index], border = NA)
  sparse_columns <- which(ages >= unsupported[1] & ages <= unsupported[2])
  if (length(sparse_columns)) {
    rect(min(sparse_columns) - 0.5, 0.06, max(sparse_columns) + 0.5, 0.94,
         border = "#000000", lwd = 1.2)
  }
  axis(2, at = c(0.5, 1.45), labels = c("Age", "Stage"), las = 1,
       tick = FALSE, cex.axis = 0.72, mgp = c(1.2, 0.4, 0))
  box(col = "#666666", lwd = 0.6)
  mtext("ATP2B1", side = 3, line = 0.2, adj = 0, cex = 0.9)

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  legend("topleft", legend = names(stage_colors), fill = stage_colors,
         border = NA, bty = "n", cex = 0.72, horiz = TRUE, xpd = NA)
  rect(0.08, 0.56, 0.13, 0.64, border = "#000000", lwd = 1.1)
  text(0.16, 0.60, "Sparse 0.2-11 y", adj = 0, cex = 0.66)
  bar_x <- seq(0.08, 0.92, length.out = length(z_palette) + 1)
  rect(head(bar_x, -1), 0.30, tail(bar_x, -1), 0.46,
       col = z_palette, border = NA)
  axis(1, at = c(0.08, 0.50, 0.92), labels = c(-display_limit, 0, display_limit),
       pos = 0.30, tck = -0.15, cex.axis = 0.64, mgp = c(1, 0.2, 0))
  text(0.5, 0.12, "Row z-score", cex = 0.7)

  par(mar = c(0.2, 4.2, 0.2, 0.5))
  draw_strip(z_index$APA, z_palette, "APA")
  par(mar = c(0.2, 3.9, 0.2, 0.7))
  draw_scatter(age_x, raw_values$APA, layer_colors["APA"], "Logit PDUI",
               unsupported)

  par(mar = c(0.2, 4.2, 0.2, 0.5))
  draw_strip(z_index$RNA, z_palette, "RNA")
  par(mar = c(0.2, 3.9, 0.2, 0.7))
  draw_scatter(age_x, raw_values$RNA, layer_colors["RNA"], "TMM logCPM",
               unsupported)

  par(mar = c(2.0, 4.2, 0.2, 0.5))
  draw_strip(z_index$Protein, z_palette, "Protein", TRUE, ages)
  par(mar = c(2.0, 3.9, 0.2, 0.7))
  draw_scatter(age_x, raw_values$Protein, layer_colors["Protein"],
               "Adjusted abundance", unsupported, TRUE)
}

# --- loading and aligning donors ---
apa <- readRDS(apa_path)
rna <- readRDS(rna_path)
protein <- readRDS(protein_path)
manifest <- fread(manifest_path)
support <- fread(support_path)

source_order_identical <- identical(colnames(apa), colnames(rna)) &&
  identical(colnames(apa), colnames(protein))

shared_ids <- Reduce(intersect, list(colnames(apa), colnames(rna), colnames(protein)))
manifest <- manifest[match(shared_ids, sample_id)]
stopifnot(!anyNA(manifest$sample_id), !anyDuplicated(manifest$sample_id))

donor_order <- order(manifest$age_years, manifest$sample_id)
manifest <- manifest[donor_order]
donor_ids <- manifest$sample_id

apa_index <- match(apa_event, rownames(apa))
rna_index <- match(rna_gene, sub("\\.[0-9]+$", "", rownames(rna)))
protein_index <- match(protein_accession, rownames(protein))
stopifnot(!is.na(apa_index), !is.na(rna_index), !is.na(protein_index))

raw_values <- list(
  APA = as.numeric(apa[apa_index, donor_ids]),
  RNA = as.numeric(rna[rna_index, donor_ids]),
  Protein = as.numeric(protein[protein_index, donor_ids])
)
z_values <- lapply(raw_values, function(x) as.numeric((x - mean(x)) / sd(x)))

unsupported <- support[is_supported == FALSE, range(grid_age)]
ages <- manifest$age_years
age_x <- age_model(ages)
sparse_donor_count <- sum(ages >= unsupported[1] & ages <= unsupported[2])

# --- summary table ---
summary <- data.table(
  gene_symbol = "ATP2B1",
  apa_event = apa_event,
  rna_gene = rna_gene,
  protein_accession = protein_accession,
  n_source_donors = ncol(apa),
  n_shared_donors = length(shared_ids),
  n_excluded = ncol(apa) - length(shared_ids),
  exclusion_reason = "none; source matrices have identical donor columns",
  n_prenatal = sum(manifest$life_stage == "Prenatal"),
  n_postnatal = sum(manifest$life_stage == "Postnatal"),
  unsupported_age_lower = unsupported[1],
  unsupported_age_upper = unsupported[2],
  n_sparse_donors_shown = sparse_donor_count,
  sparse_visualization = "observed points shown; dashed bridge excludes sparse donors from fitting",
  spearman_apa_rna = cor(raw_values$APA, raw_values$RNA, method = "spearman"),
  spearman_apa_protein = cor(raw_values$APA, raw_values$Protein, method = "spearman"),
  spearman_rna_protein = cor(raw_values$RNA, raw_values$Protein, method = "spearman")
)
fwrite(summary, summary_path, sep = "\t")

# --- figure 09: aligned ATP2B1 donor tracks and raw observations ---
stage_colors <- c(Prenatal = "#56B4E9", Postnatal = "#E69F00")
stage_index <- match(manifest$life_stage, names(stage_colors))
age_palette <- hcl.colors(100, "Viridis")
age_index <- pmax(1L, pmin(100L, 1L + floor(
  99 * (age_x - min(age_x)) / diff(range(age_x))
)))
z_palette <- colorRampPalette(c("#0072B2", "#F7F7F7", "#D55E00"))(101)
z_breaks <- seq(-display_limit, display_limit, length.out = 102)
z_index <- lapply(z_values, function(x) {
  findInterval(pmax(-display_limit, pmin(display_limit, x)), z_breaks,
               all.inside = TRUE)
})
layer_colors <- c(APA = "#009E73", RNA = "#0072B2", Protein = "#D55E00")

png(paste0(figure_path, ".png"), width = 7.2, height = 5.8,
    units = "in", res = 300, type = "cairo")
draw_figure()
dev.off()

pdf(paste0(figure_path, ".pdf"), width = 7.2, height = 5.8,
    useDingbats = FALSE)
draw_figure()
dev.off()

# --- layout validation ---
validation <- data.table(
  check = c(
    "source_column_order_identical", sprintf("shared_donor_count_%d", expected_donor_count),
    "manifest_mapping_complete", "stable_age_sample_order",
    "apa_event_mapped", "rna_gene_mapped", "protein_accession_mapped",
    "focal_values_complete", sprintf("sparse_donors_%d_shown", expected_sparse_donor_count),
    "unsupported_interval_marked", "sparse_bridge_defined",
    "png_written", "pdf_written"
  ),
  value = c(
    source_order_identical, length(shared_ids) == expected_donor_count,
    !anyNA(manifest$sample_id),
    identical(order(manifest$age_years, manifest$sample_id), seq_len(nrow(manifest))),
    !is.na(apa_index), !is.na(rna_index), !is.na(protein_index),
    all(vapply(raw_values, function(x) all(is.finite(x)), logical(1))),
    sparse_donor_count == expected_sparse_donor_count,
    length(which(ages >= unsupported[1] & ages <= unsupported[2])) > 0L,
    sum(age_x <= age_model(unsupported[1])) >= 10L &&
      sum(age_x >= age_model(unsupported[2])) >= 10L,
    file.exists(paste0(figure_path, ".png")),
    file.exists(paste0(figure_path, ".pdf"))
  )
)
validation[, status := fifelse(value, "PASS", "FAIL")]
fwrite(validation, validation_path, sep = "\t")
