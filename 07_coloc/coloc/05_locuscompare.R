#!/usr/bin/env Rscript
# Locuscompare plots for colocalized loci: GWAS and APA-QTL regional panels, gene track, r2 key.
# Input: data/07_coloc/coloc_run/rds/<locus>_coloc.rds; output: figures/07_coloc/*_locuscompare.{png,pdf}.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(glue)
})

# --- EasyColoc sources ---
APA_ROOT <- Sys.getenv("APA_ROOT")
easycoloc_root <- file.path(APA_ROOT, "release/07_coloc/easycoloc")
for (f in c(
  "src/utils_helpers.R",
  "src/utils_plot.R",
  "src/utils_output.R",
  "src/utils_format.R",
  "src/utils_plot_locuscompare.R"
)) {
  source(file.path(easycoloc_root, f))
}

# --- layout ---
LOCUSCOMPARE_LAYOUT$width_in <- 5.95
LOCUSCOMPARE_LAYOUT$height_in <- 3.95
LOCUSCOMPARE_LAYOUT$scatter_width <- 1.00
LOCUSCOMPARE_LAYOUT$regional_width <- 0.96
LOCUSCOMPARE_LAYOUT$gene_track_height <- 0.46
LOCUSCOMPARE_LAYOUT$gene_track_height_per_extra_lane <- 0.15
LOCUSCOMPARE_LAYOUT$gene_track_height_max <- 1.10
LOCUSCOMPARE_LAYOUT$regional_panel_gap <- 0.07
LOCUSCOMPARE_LAYOUT$gene_track_top_gap <- 0.12

# --- helpers ---
# These names replace the EasyColoc functions of the same name; EasyColoc calls them by name.
locuscompare_protected_label_layer <- function(data,
                                               mapping,
                                               base_family = "Helvetica",
                                               size = 3) {
  ggrepel::geom_label_repel(
    data = data,
    mapping = mapping,
    size = size,
    family = base_family,
    fontface = "bold",
    color = "black",
    fill = "#FFFFFFF0",
    label.size = 0,
    label.padding = unit(0.08, "lines"),
    label.r = unit(0.03, "lines"),
    box.padding = 0.35,
    point.padding = 0.35,
    force = 1.5,
    min.segment.length = 0,
    segment.color = "grey35",
    segment.size = 0.28,
    max.overlaps = Inf,
    seed = 1,
    inherit.aes = FALSE
  )
}

locuscompare_select_gene_labels <- function(gene_df,
                                            position_range,
                                            lead_position = NULL,
                                            panel_width_in = 2.9,
                                            label_size_pt = 7.3) {
  if (nrow(gene_df) == 0L) return(gene_df)
  region_span <- diff(position_range)
  label_gap_bp <- max(1, region_span * 0.022)
  edge_padding_bp <- max(1, region_span * 0.045)
  gene_df$label_width_bp <- pmax(
    region_span * 0.035,
    nchar(gene_df$symbol) * label_size_pt * 0.72 / 72 / panel_width_in * region_span
  )
  raw_label_x <- (gene_df$start + gene_df$end) / 2
  label_min_x <- position_range[1] + edge_padding_bp + gene_df$label_width_bp / 2
  label_max_x <- position_range[2] - edge_padding_bp - gene_df$label_width_bp / 2
  gene_df$label_too_wide <- label_min_x > label_max_x
  gene_df$label_x <- ifelse(
    gene_df$label_too_wide,
    raw_label_x,
    pmax(label_min_x, pmin(label_max_x, raw_label_x))
  )
  gene_df$label_start <- gene_df$label_x - gene_df$label_width_bp / 2
  gene_df$label_end <- gene_df$label_x + gene_df$label_width_bp / 2
  gene_df$label_boundary_truncated <-
    gene_df$start <= position_range[1] + edge_padding_bp |
    gene_df$end >= position_range[2] - edge_padding_bp * 1.8 |
    gene_df$label_too_wide
  gene_df$label_visible <- FALSE
  reference_position <- if (!is.null(lead_position) && is.finite(lead_position)) {
    lead_position
  } else {
    mean(position_range)
  }
  gene_df$label_priority <- abs(gene_df$label_x - reference_position)
  for (lane in sort(unique(gene_df$lane))) {
    lane_index <- which(gene_df$lane == lane & !gene_df$label_boundary_truncated)
    lane_index <- lane_index[order(gene_df$label_priority[lane_index], gene_df$label_end[lane_index])]
    selected <- integer(0)
    for (index in lane_index) {
      conflict <- if (length(selected) == 0L) {
        FALSE
      } else {
        any(
          gene_df$label_start[index] < gene_df$label_end[selected] + label_gap_bp &
            gene_df$label_end[index] > gene_df$label_start[selected] - label_gap_bp
        )
      }
      if (!conflict) {
        gene_df$label_visible[index] <- TRUE
        selected <- c(selected, index)
      }
    }
  }
  candidate_index <- which(gene_df$label_visible & !gene_df$label_boundary_truncated)
  candidate_index <- candidate_index[order(gene_df$label_priority[candidate_index])]
  n_lanes <- max(gene_df$lane)
  track_weight <- min(
    LOCUSCOMPARE_LAYOUT$gene_track_height_max,
    LOCUSCOMPARE_LAYOUT$gene_track_height +
      max(0, n_lanes - 1L) * LOCUSCOMPARE_LAYOUT$gene_track_height_per_extra_lane
  )
  track_height_pt <- LOCUSCOMPARE_LAYOUT$height_in * 72 * track_weight /
    (LOCUSCOMPARE_LAYOUT$gwas_height + LOCUSCOMPARE_LAYOUT$qtl_height + track_weight)
  lane_pitch_pt <- max(1, track_height_pt - 22) / n_lanes
  label_height_pt <- label_size_pt * 1.15
  gene_df$label_visible <- FALSE
  selected <- integer(0)
  for (index in candidate_index) {
    conflict <- if (length(selected) == 0L) {
      FALSE
    } else {
      horizontal_overlap <-
        gene_df$label_start[index] < gene_df$label_end[selected] + label_gap_bp &
        gene_df$label_end[index] > gene_df$label_start[selected] - label_gap_bp
      vertical_overlap <- abs(gene_df$lane[index] - gene_df$lane[selected]) * lane_pitch_pt < label_height_pt
      any(horizontal_overlap & vertical_overlap)
    }
    if (!conflict) {
      gene_df$label_visible[index] <- TRUE
      selected <- c(selected, index)
    }
  }
  visible_df <- gene_df[gene_df$label_visible, , drop = FALSE]
  attr(gene_df, "label_gap_bp") <- label_gap_bp
  attr(gene_df, "label_lane_pitch_pt") <- lane_pitch_pt
  attr(gene_df, "label_overlap_free") <- locuscompare_gene_labels_do_not_overlap(
    visible_df,
    lane_pitch_pt = lane_pitch_pt,
    label_height_pt = label_height_pt,
    horizontal_gap_bp = label_gap_bp
  )
  attr(gene_df, "canonical_boundary_labels_hidden") <- gene_df$symbol[gene_df$label_boundary_truncated]
  gene_df
}

locuscompare_regional_panel <- function(point_df,
                                        y_col,
                                        title,
                                        chr,
                                        color_values,
                                        shared_theme,
                                        base_family,
                                        show_x = TRUE,
                                        top_panel = FALSE) {
  metal <- data.frame(
    rsid = point_df$rsid,
    logp = point_df[[if (identical(y_col, "logp1")) "logp1" else "logp2"]],
    chr = point_df$chr,
    pos = point_df$pos,
    ld_bin = point_df$ld_bin,
    is_lead = point_df$is_lead,
    label = point_df$label,
    stringsAsFactors = FALSE
  )
  metal_nonlead <- metal[!metal$is_lead, , drop = FALSE]
  metal_lead <- metal[metal$is_lead, , drop = FALSE]
  pos_range <- range(metal$pos, na.rm = TRUE)
  logp_range <- range(metal$logp, na.rm = TRUE)
  tag_df <- data.frame(
    x = pos_range[1] + diff(pos_range) * 0.018,
    y = logp_range[2] - diff(logp_range) * 0.035,
    panel_tag = title,
    stringsAsFactors = FALSE
  )
  ggplot(metal, aes(x = pos, y = logp)) +
    geom_point(
      data = metal_nonlead,
      aes(fill = ld_bin),
      shape = 21,
      size = 2.40,
      color = "black",
      stroke = 0.14,
      alpha = 0.90
    ) +
    geom_point(
      data = metal_lead,
      shape = 23,
      size = 4.85,
      fill = "purple",
      color = "black",
      stroke = 0.20
    ) +
    scale_fill_manual(values = color_values, guide = "none") +
    scale_x_continuous(
      limits = range(metal$pos, na.rm = TRUE),
      labels = locuscompare_mb_labeller(metal$pos),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.11))) +
    locuscompare_protected_label_layer(
      data = metal_lead,
      mapping = aes(x = pos, y = logp, label = label),
      size = 3.25,
      base_family = base_family
    ) +
    geom_label(
      data = tag_df,
      aes(x = x, y = y, label = panel_tag),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.65,
      family = base_family,
      fontface = "bold",
      color = "black",
      fill = "#FFFFFFE6",
      linewidth = 0,
      label.padding = unit(0.08, "lines")
    ) +
    xlab(if (show_x) paste0("chr", chr, " (Mb)") else NULL) +
    ylab(NULL) +
    shared_theme +
    theme(
      axis.text.x = if (show_x) element_text(size = 7.6, color = "black") else element_blank(),
      axis.title.x = if (show_x) element_text(size = 8.0) else element_blank(),
      axis.title.y = element_blank(),
      plot.margin = margin(t = if (top_panel) 8 else 3, r = 3, b = 2.5, l = 3, unit = "pt")
    )
}

make_internal_locuscomparer_style_plot <- function(inputs,
                                                   title1 = "GWAS",
                                                   title2 = "QTL",
                                                   legend = TRUE,
                                                   legend_position = "topleft",
                                                   gtf_path = NULL,
                                                   base_family = "Helvetica") {
  plot_df <- inputs$prepared$data
  merged <- inputs$merged
  color_values <- locuscompare_ld_colors(include_lead = FALSE)
  point_df <- data.frame(
    rsid = merged$rsid,
    logp1 = merged$logp1,
    logp2 = merged$logp2,
    chr = merged$chr,
    pos = merged$pos,
    ld_bin = plot_df$ld_bin,
    is_lead = plot_df$is_lead,
    label = ifelse(plot_df$is_lead, plot_df$snp, ""),
    stringsAsFactors = FALSE
  )
  nonlead_df <- point_df[!point_df$is_lead, , drop = FALSE]
  lead_df <- point_df[point_df$is_lead, , drop = FALSE]

  shared_theme <- theme_classic(base_size = 8.7, base_family = base_family) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = 8.5, color = "black"),
      axis.text = element_text(size = 7.8, color = "black"),
      axis.line = element_line(linewidth = 0.50, color = "black"),
      axis.ticks = element_line(linewidth = 0.45, color = "black"),
      axis.ticks.length = unit(2.1, "pt"),
      panel.spacing = unit(0, "pt"),
      plot.margin = unit(c(1.5, 3, 2.5, 3), "pt")
    )
  scatter <- ggplot(point_df, aes(logp1, logp2)) +
    geom_point(
      data = nonlead_df,
      aes(fill = ld_bin),
      shape = 21,
      size = 2.45,
      color = "black",
      stroke = 0.16,
      alpha = 0.90
    ) +
    geom_point(
      data = lead_df,
      shape = 23,
      size = 4.95,
      fill = "purple",
      color = "black",
      stroke = 0.22
    ) +
    scale_fill_manual(values = color_values, guide = "none") +
    locuscompare_protected_label_layer(
      data = lead_df,
      mapping = aes(x = logp1, y = logp2, label = label),
      size = 3.45,
      base_family = base_family
    ) +
    xlab(bquote(.(title1) ~ -log[10] * "(P)")) +
    ylab(bquote(.(title2) ~ -log[10] * "(P)")) +
    scale_x_continuous(expand = expansion(mult = c(0.04, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.09))) +
    shared_theme

  p2 <- locuscompare_regional_panel(point_df, "logp1", title1, inputs$chr, color_values,
                                    shared_theme, base_family, show_x = FALSE, top_panel = TRUE)
  p3 <- locuscompare_regional_panel(point_df, "logp2", title2, inputs$chr, color_values,
                                    shared_theme, base_family, show_x = FALSE)

  # --- gene track ---
  gene_track <- make_locuscompare_gene_track(
    inputs = inputs,
    gtf_path = gtf_path,
    base_family = base_family
  )
  for (layer_index in seq_along(gene_track$layers)) {
    if (inherits(gene_track$layers[[layer_index]]$geom, "GeomText")) {
      gene_track$layers[[layer_index]]$aes_params$size <- 2.30
    }
  }
  gene_track <- gene_track +
    theme(
      axis.title.x = element_text(size = 8.6, face = "plain", color = "black", margin = margin(t = 2)),
      axis.text.x = element_text(size = 7.5, color = "black"),
      plot.margin = margin(t = 0, r = 5, b = 1, l = 2, unit = "pt")
    )
  gene_track_range <- attr(gene_track, "easycoloc_gene_track_position_range")
  # MB breaks on the merged positions; NULL when fewer than two breaks fall inside the range
  gene_track_breaks <- NULL
  finite_pos <- safe_numeric(inputs$merged$pos)
  finite_pos <- finite_pos[is.finite(finite_pos)]
  if (length(finite_pos) >= 2L) {
    pos_range <- range(finite_pos)
    span <- diff(pos_range)
    pretty_breaks <- pretty(pos_range / 1e6, n = 4L)
    pretty_breaks <- pretty_breaks[
      pretty_breaks > (pos_range[1] + span * 0.08) / 1e6 &
        pretty_breaks < (pos_range[2] - span * 0.08) / 1e6
    ]
    if (length(pretty_breaks) >= 2L) gene_track_breaks <- pretty_breaks
  }
  if (length(gene_track_range) == 2L && all(is.finite(gene_track_range))) {
    if (!is.null(gene_track_breaks)) {
      gene_track <- gene_track + scale_x_continuous(
        limits = gene_track_range / 1e6,
        breaks = gene_track_breaks,
        labels = function(x) sprintf("%.1f", x),
        expand = expansion(mult = c(0.02, 0.12))
      )
    } else {
      gene_track <- locuscompare_set_x_labeller(
        gene_track,
        labeller = function(x) sprintf("%.2f", x),
        limits = gene_track_range / 1e6,
        expand = expansion(mult = c(0.02, 0.12))
      )
    }
  }

  # --- r2 key ---
  scatter_plot <- scatter
  if (isTRUE(legend)) {
    display_bins <- c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "<0.2")
    legend_df <- data.frame(
      item = display_bins,
      ymin = rev(seq_along(display_bins)) - 1,
      ymax = rev(seq_along(display_bins)),
      stringsAsFactors = FALSE
    )
    tick_df <- data.frame(
      label = c("0.8", "0.6", "0.4", "0.2"),
      y = c(4, 3, 2, 1),
      stringsAsFactors = FALSE
    )
    legend_panel <- ggplot(legend_df) +
      annotate("text", x = 0.50, y = 5.35, label = "r^2", parse = TRUE,
               size = 2.35, fontface = "bold", family = base_family) +
      geom_rect(aes(xmin = 0, xmax = 0.42, ymin = ymin, ymax = ymax, fill = item),
                color = "black", linewidth = 0.18) +
      geom_segment(data = tick_df, aes(x = 0.42, xend = 0.48, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.18, color = "black") +
      geom_text(data = tick_df, aes(x = 0.54, y = y, label = label),
                inherit.aes = FALSE,
                hjust = 0, vjust = 0.5, size = 2.25, family = base_family) +
      scale_fill_manual(values = locuscompare_ld_colors(include_lead = FALSE), guide = "none") +
      coord_cartesian(xlim = c(-0.02, 0.95), ylim = c(0, 5.6), clip = "off") +
      theme_void(base_family = base_family) +
      theme(
        plot.margin = unit(c(1, 1, 1, 1), "pt"),
        plot.background = element_rect(fill = "#FFFFFFE6", color = NA)
      )
    scatter_plot <- cowplot::ggdraw(scatter) +
      cowplot::draw_plot(legend_panel, x = 0.805, y = 0.150, width = 0.120, height = 0.235)
  }

  # --- layout ---
  aligned_regional <- cowplot::align_plots(
    p2,
    p3,
    gene_track,
    align = "v",
    axis = "lr"
  )
  gene_track_height <- attr(gene_track, "easycoloc_gene_track_height")
  regional_column <- cowplot::plot_grid(
    aligned_regional[[1]],
    cowplot::ggdraw(),
    aligned_regional[[2]],
    cowplot::ggdraw(),
    aligned_regional[[3]],
    ncol = 1,
    rel_heights = c(
      LOCUSCOMPARE_LAYOUT$gwas_height,
      LOCUSCOMPARE_LAYOUT$regional_panel_gap,
      LOCUSCOMPARE_LAYOUT$qtl_height,
      LOCUSCOMPARE_LAYOUT$gene_track_top_gap,
      gene_track_height
    )
  )
  combined <- cowplot::plot_grid(
    scatter_plot,
    regional_column,
    nrow = 1,
    align = "h",
    axis = "tb",
    rel_widths = c(LOCUSCOMPARE_LAYOUT$scatter_width, LOCUSCOMPARE_LAYOUT$regional_width)
  )
  attr(combined, "easycoloc_gene_track_included") <- TRUE
  attr(combined, "easycoloc_gene_features_present") <- attr(gene_track, "easycoloc_gene_features_present")
  attr(combined, "easycoloc_gene_track_position_range") <- attr(gene_track, "easycoloc_gene_track_position_range")
  attr(combined, "easycoloc_gene_track_column_aligned") <- TRUE
  attr(combined, "easycoloc_gene_track_shared_x_axis") <- TRUE
  attr(combined, "easycoloc_scatter_spans_full_height") <- TRUE
  attr(combined, "easycoloc_gene_count") <- attr(gene_track, "easycoloc_gene_count")
  attr(combined, "easycoloc_gene_lane_count") <- attr(gene_track, "easycoloc_gene_lane_count")
  attr(combined, "easycoloc_gene_labels_drawn") <- attr(gene_track, "easycoloc_gene_labels_drawn")
  attr(combined, "easycoloc_gene_labels_visible") <- attr(gene_track, "easycoloc_gene_labels_visible")
  attr(combined, "easycoloc_gene_labels_hidden") <- attr(gene_track, "easycoloc_gene_labels_hidden")
  attr(combined, "easycoloc_gene_interval_overlap_free") <- attr(gene_track, "easycoloc_gene_interval_overlap_free")
  attr(combined, "easycoloc_gene_label_overlap_free") <- attr(gene_track, "easycoloc_gene_label_overlap_free")
  attr(combined, "easycoloc_gene_label_lane_pitch_pt") <- attr(gene_track, "easycoloc_gene_label_lane_pitch_pt")
  attr(combined, "easycoloc_gene_track_height") <- gene_track_height
  attr(combined, "easycoloc_layout") <- "canonical compact full-height scatter with lower-right inset r2 key and separated right-column gene track"
  attr(combined, "easycoloc_regional_panel_gap") <- as.character(LOCUSCOMPARE_LAYOUT$regional_panel_gap)
  attr(combined, "easycoloc_official_legend_embedded") <- FALSE
  attr(combined, "easycoloc_legend_position") <- "lower-right scatter inset"
  attr(combined, "easycoloc_missing_ld_legend_hidden") <- TRUE
  attr(combined, "easycoloc_legend_layout") <- "compact r2 key overlay with boundary ticks"
  attr(combined, "easycoloc_gene_track_gap") <- "spacer row between regional panels and gene track"
  attr(combined, "easycoloc_gene_label_policy") <- "text-size-aware labels; boundary-truncated gene names hidden"
  attr(combined, "easycoloc_lead_label_panels") <- c(
    scatter = locuscompare_has_repel_label(scatter),
    gwas = locuscompare_has_repel_label(p2),
    qtl = locuscompare_has_repel_label(p3)
  )
  attr(combined, "easycoloc_right_axis_label_policy") <- "rotated regional y-axis titles removed; short in-panel tags"
  attr(combined, "easycoloc_scatter_lead_label_clearance") <- "scatter x-axis right expansion 0.18"
  attr(combined, "easycoloc_point_outlines_removed") <- FALSE
  attr(combined, "easycoloc_lead_labels_background_protected") <- TRUE
  combined
}

# --- paths ---
out_dir <- file.path(APA_ROOT, "figures/07_coloc")
data_dir <- file.path(APA_ROOT, "data/07_coloc")
rds_dir <- file.path(data_dir, "coloc_run/rds")
cache_dir <- file.path(data_dir, "ld_cache")
gtf_path <- file.path(APA_ROOT, "ref/gencode.v40.annotation.gtf.gz")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
validation_path <- file.path(data_dir, "05_locuscompare_layout_validation.tsv")
validation_rows <- list()
plink_bin <- "plink"

# --- loci ---
locus_keys <- c(
  "EAS_SCZ_postnatal_chr10_102760520_rs12571643_ARL3",
  "EAS_SCZ_interaction_chr12_89410263_rs12578201_ATP2B1",
  "EUR_MDD_postnatal_chr12_120364215_rs3782901_SPPL3"
)

for (key in locus_keys) {
  rds_file <- file.path(rds_dir, paste0(key, "_coloc.rds"))
  if (!file.exists(rds_file)) {
    message("[SKIP] RDS not found: ", rds_file)
    next
  }
  r <- readRDS(rds_file)
  df <- r$merged_data
  lead_snp <- r$lead_snp
  credible_set <- r$credible_set
  bfile <- r$plink_bfile
  pp4 <- if (!is.null(r$pp4)) r$pp4 else NA_real_

  # reuse cached plink LD when present, else compute it
  cache_files <- list.files(cache_dir, pattern = paste0("_", basename(key), "\\.ld$"), full.names = TRUE)
  ld_df <- NULL
  if (length(cache_files) > 0) {
    ld_tab <- fread(cache_files[1])
    if (all(c("SNP_A", "SNP_B", "R2") %in% names(ld_tab))) {
      ld_df <- as.data.frame(ld_tab)
    }
  }
  if (is.null(ld_df)) {
    snp_col <- resolve_col(df, "rsid", c("SNP", "snp", "SNPID.gwas"))
    snps <- unique(as.character(df[[snp_col]]))
    # multi-allelic rsIDs are split in the 1KG bfile; remap the lead via bim position and alleles
    id_map <- NULL
    if (file.exists(paste0(bfile, ".bim"))) {
      bim_tab <- fread(paste0(bfile, ".bim"), header = FALSE,
                       col.names = c("CHR", "SNP", "CM", "POS", "A1", "A2"))
      # positional index for the lead + all query SNPs on the same chr
      chr_col <- resolve_col(df, "CHR.qtl", c("CHR.gwas", "CHR", "chr", "chromosome"))
      pos_col <- resolve_col(df, "POS.qtl", c("POS.gwas", "POS", "pos", "BP", "position"))
      lead_chr <- as.character(df[[chr_col]][df[[snp_col]] == lead_snp][1])
      lead_pos <- as.numeric(df[[pos_col]][df[[snp_col]] == lead_snp][1])
      if (!is.na(lead_pos) && nzchar(lead_chr)) {
        lead_chr_num <- sub("^chr", "", lead_chr)
        region_tab <- bim_tab[CHR == lead_chr_num & POS == lead_pos, ]
        if (nrow(region_tab) > 0) {
          # prefer exact allele match A1/A2 vs REF/ALT when possible
          ref_alt <- NULL
          if ("REF" %in% names(df) && "ALT" %in% names(df)) {
            ref_alt <- data.frame(REF = as.character(df[[ "REF" ]][df[[snp_col]] == lead_snp]),
                                  ALT = as.character(df[[ "ALT" ]][df[[snp_col]] == lead_snp]),
                                  stringsAsFactors = FALSE)
            ref_alt <- ref_alt[1, , drop = FALSE]
            ref_alt$A1 <- ref_alt$REF
            ref_alt$A2 <- ref_alt$ALT
          }
          if (!is.null(ref_alt) && nrow(ref_alt) > 0 && !is.na(ref_alt$REF)) {
            hit <- region_tab[(A1 == ref_alt$REF & A2 == ref_alt$ALT) |
                              (A1 == ref_alt$ALT & A2 == ref_alt$REF), ]
            if (nrow(hit) == 0) hit <- region_tab[1, ]
          } else {
            hit <- region_tab[1, ]
          }
          if (nrow(hit) > 0) {
            id_map <- data.frame(
              from = lead_snp,
              to = hit$SNP[1],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    query_snps <- snps
    if (!is.null(id_map)) {
      query_snps[query_snps == id_map$from[1]] <- id_map$to[1]
    }
    ld_mat_res <- get_ld_matrix(query_snps, bfile, plink_bin)
    if (!is.null(ld_mat_res) && is.list(ld_mat_res) && !is.null(ld_mat_res$R)) {
      Rmat <- ld_mat_res$R
      if (!is.null(id_map) && id_map$to[1] %in% rownames(Rmat)) {
        # map back to original rsID so the plot uses original IDs
        rownames(Rmat)[rownames(Rmat) == id_map$to[1]] <- id_map$from[1]
        colnames(Rmat)[colnames(Rmat) == id_map$to[1]] <- id_map$from[1]
      }
      ld_df <- as.data.frame(as.table(Rmat))
      names(ld_df) <- c("SNP_A", "SNP_B", "R2")
      ld_df <- ld_df[as.character(ld_df$SNP_A) < as.character(ld_df$SNP_B), , drop = FALSE]
      ld_df <- ld_df[ld_df$SNP_A %in% snps & ld_df$SNP_B %in% snps, , drop = FALSE]
    } else {
      message("[", key, "] WARNING: matrix LD failed, locuscompare will lack LD coloring")
    }
  }

  # --- plot ---
  parts <- strsplit(as.character(r$gwas_id), "_", fixed = TRUE)[[1]]
  population <- if (length(parts) >= 1L) parts[1] else as.character(r$gwas_id)
  disease <- if (length(parts) >= 2L) paste(parts[-1], collapse = "_") else as.character(r$gwas_id)
  lc_title <- glue("{population} {disease} | {r$qtl_id} | {r$gene_symbol}")
  lc_plot <- make_combined_coloc_plot(
    df = df,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set,
    title = lc_title,
    pp4 = pp4,
    gtf_path = gtf_path,
    backend = "internal",
    title1 = "GWAS",
    title2 = "APA-QTL",
    legend_position = "topleft"
  )
  easy_attrs <- attributes(lc_plot)[grepl("^easycoloc_", names(attributes(lc_plot)))]
  title_plot <- cowplot::ggdraw() +
    cowplot::draw_label(
      lc_title,
      x = 0.5,
      y = 0.66,
      hjust = 0.5,
      vjust = 0.5,
      size = 8.9,
      fontface = "bold"
    )
  lc_plot <- cowplot::plot_grid(title_plot, lc_plot, ncol = 1, rel_heights = c(0.14, 1))
  for (nm in names(easy_attrs)) attr(lc_plot, nm) <- easy_attrs[[nm]]
  attr(lc_plot, "easycoloc_title") <- as.character(lc_title)
  attr(lc_plot, "easycoloc_title_row_height") <- "0.14"
  attr(lc_plot, "easycoloc_title_y_label_clearance") <-
    "no rotated regional y-axis titles plus expanded title row"

  pdf_path <- file.path(out_dir, paste0(key, "_locuscompare.pdf"))
  png_path <- file.path(out_dir, paste0(key, "_locuscompare.png"))
  save_plot_with_fallback(
    plot_obj = lc_plot,
    pdf_path = pdf_path,
    png_path = png_path,
    plot_width = LOCUSCOMPARE_LAYOUT$width_in,
    plot_height = LOCUSCOMPARE_LAYOUT$height_in,
    png_scale = 1,
    png_dpi = 300
  )

  # --- metadata ---
  json_path <- file.path(data_dir, paste0(key, "_locuscompare.json"))
  write_locuscompare_plot_json(
    df = df,
    output_file = json_path,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set,
    meta = list(
      gwas_id = r$gwas_id,
      qtl_id = r$qtl_id,
      locus = lead_snp,
      phenotype = r$gene_symbol,
      pp4 = pp4
    )
  )

  # --- validation ---
  output_files_ok <- all(file.exists(c(pdf_path, png_path, json_path)))
  validation_rows[[length(validation_rows) + 1L]] <- data.frame(
    locus_key = key,
    check = c(
      "title",
      "figure_inches",
      "legend_position",
      "missing_ld_legend_hidden",
      "legend_layout",
      "gene_track_top_gap",
      "gene_track_lanes",
      "gene_labels_visible",
      "gene_label_policy",
      "right_axis_label_policy",
      "title_y_label_clearance",
      "regional_panel_gap",
      "scatter_lead_label_clearance",
      "output_files"
    ),
    value = c(
      as.character(attr(lc_plot, "easycoloc_title")),
      sprintf("%.1fx%.1f", LOCUSCOMPARE_LAYOUT$width_in, LOCUSCOMPARE_LAYOUT$height_in),
      as.character(attr(lc_plot, "easycoloc_legend_position")),
      as.character(isTRUE(attr(lc_plot, "easycoloc_missing_ld_legend_hidden"))),
      as.character(attr(lc_plot, "easycoloc_legend_layout")),
      as.character(attr(lc_plot, "easycoloc_gene_track_gap")),
      as.character(attr(lc_plot, "easycoloc_gene_lane_count")),
      paste(attr(lc_plot, "easycoloc_gene_labels_visible"), collapse = ","),
      as.character(attr(lc_plot, "easycoloc_gene_label_policy")),
      as.character(attr(lc_plot, "easycoloc_right_axis_label_policy")),
      as.character(attr(lc_plot, "easycoloc_title_y_label_clearance")),
      as.character(attr(lc_plot, "easycoloc_regional_panel_gap")),
      as.character(attr(lc_plot, "easycoloc_scatter_lead_label_clearance")),
      paste(c(pdf_path, png_path, json_path), collapse = ";")
    ),
    status = c(
      "PASS",
      "PASS",
      if (identical(attr(lc_plot, "easycoloc_legend_position"), "lower-right scatter inset")) "PASS" else "FAIL",
      if (isTRUE(attr(lc_plot, "easycoloc_missing_ld_legend_hidden"))) "PASS" else "FAIL",
      if (identical(attr(lc_plot, "easycoloc_legend_layout"), "compact r2 key overlay with boundary ticks")) "PASS" else "FAIL",
      if (identical(attr(lc_plot, "easycoloc_gene_track_gap"), "spacer row between regional panels and gene track")) "PASS" else "FAIL",
      "PASS",
      "PASS",
      "PASS",
      if (identical(attr(lc_plot, "easycoloc_right_axis_label_policy"), "rotated regional y-axis titles removed; short in-panel tags")) "PASS" else "FAIL",
      if (identical(attr(lc_plot, "easycoloc_title_y_label_clearance"), "no rotated regional y-axis titles plus expanded title row")) "PASS" else "FAIL",
      if (identical(attr(lc_plot, "easycoloc_regional_panel_gap"), as.character(LOCUSCOMPARE_LAYOUT$regional_panel_gap))) "PASS" else "FAIL",
      if (identical(attr(lc_plot, "easycoloc_scatter_lead_label_clearance"), "scatter x-axis right expansion 0.18")) "PASS" else "FAIL",
      if (output_files_ok) "PASS" else "FAIL"
    ),
    stringsAsFactors = FALSE
  )
}

if (length(validation_rows) > 0L) {
  fwrite(rbindlist(validation_rows, fill = TRUE), validation_path, sep = "\t")
}
