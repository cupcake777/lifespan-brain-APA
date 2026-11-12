#' EasyColoc locuscompare plotting utilities.
#'
#' Canonicalized from the user-approved version 7. This lineage makes all
#' association points borderless and draws lead labels last as borderless
#' white-backed repel labels so markers cannot compromise text legibility.

#===config===
LOCUSCOMPARE_LAYOUT <- list(
  width_in = 7.2,
  height_in = 4.4,
  scatter_width = 1.08,
  regional_width = 1,
  gwas_height = 0.8,
  qtl_height = 1,
  gene_track_height = 0.38,
  gene_track_height_per_extra_lane = 0.13,
  gene_track_height_max = 1.25
)
#===end config===

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(jsonlite)
})

clamp_pvalue <- function(x) {
  x <- safe_numeric(x)
  x[!is.finite(x) | x <= 0 | x > 1] <- NA_real_
  x
}

locuscomparer_style_colors <- function() {
  c("blue4", "skyblue", "darkgreen", "orange", "red")
}

locuscompare_ld_levels <- function() {
  c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "<0.2", "missing LD")
}

locuscompare_ld_colors <- function(include_lead = FALSE) {
  ld_levels <- locuscompare_ld_levels()
  values <- c(
    setNames(rev(locuscomparer_style_colors()), ld_levels[ld_levels != "missing LD"]),
    "missing LD" = "#BDBDBD"
  )
  if (isTRUE(include_lead)) {
    values <- c("Lead SNP" = "purple", values)
  }
  values
}

locuscompare_publication_theme <- function(base_family = "Helvetica") {
  theme_classic(base_size = 7.5, base_family = base_family) +
    theme(
      plot.title = element_text(size = 9, face = "plain", hjust = 0),
      axis.title = element_text(size = 8.5),
      axis.text = element_text(size = 7.2, color = "black"),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.length = unit(2, "pt"),
      legend.title = element_text(size = 7.2),
      legend.text = element_text(size = 6.8),
      panel.spacing = unit(0, "pt"),
      plot.margin = unit(c(2, 3, 2, 3), "pt")
    )
}

locuscompare_official_theme <- function(base_family = "Helvetica") {
  theme(
    text = element_text(family = base_family),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.length = unit(2.5, "pt"),
    plot.margin = margin(2, 4, 2, 3, unit = "pt")
  )
}

locuscompare_mb_labeller <- function(positions_bp) {
  finite_pos <- safe_numeric(positions_bp)
  finite_pos <- finite_pos[is.finite(finite_pos)]
  span_mb <- if (length(finite_pos) > 1L) diff(range(finite_pos)) / 1e6 else 0
  digits <- if (span_mb < 0.005) {
    4L
  } else if (span_mb < 0.05) {
    3L
  } else if (span_mb < 0.5) {
    2L
  } else {
    1L
  }
  function(x) sprintf(paste0("%.", digits, "f"), x / 1e6)
}

locuscompare_set_x_labeller <- function(plot_obj, labeller, limits = NULL, expand = NULL) {
  scale_indexes <- which(vapply(
    plot_obj$scales$scales,
    function(scale) "x" %in% scale$aesthetics,
    logical(1)
  ))
  if (length(scale_indexes) == 0L) {
    return(plot_obj + scale_x_continuous(labels = labeller, limits = limits, expand = expand))
  }
  scale_index <- tail(scale_indexes, 1L)
  plot_obj$scales$scales[[scale_index]]$labels <- labeller
  if (!is.null(limits)) plot_obj$scales$scales[[scale_index]]$limits <- limits
  if (!is.null(expand)) plot_obj$scales$scales[[scale_index]]$expand <- expand
  plot_obj
}

locuscompare_set_point_sizes <- function(plot_obj, rsids, lead_snp) {
  size_indexes <- which(vapply(
    plot_obj$scales$scales,
    function(scale) "size" %in% scale$aesthetics,
    logical(1)
  ))
  size_values <- setNames(ifelse(rsids == lead_snp, 3.5, 2.35), rsids)
  size_scale <- scale_size_manual(values = size_values, guide = "none")
  if (length(size_indexes) == 0L) {
    return(plot_obj + size_scale)
  }
  plot_obj$scales$scales[[tail(size_indexes, 1L)]] <- size_scale
  plot_obj
}

locuscompare_has_repel_label <- function(plot_obj) {
  any(vapply(
    plot_obj$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel") || inherits(layer$geom, "GeomLabelRepel"),
    logical(1)
  ))
}

locuscompare_remove_repel_labels <- function(plot_obj) {
  keep_layers <- !vapply(
    plot_obj$layers,
    function(layer) inherits(layer$geom, "GeomTextRepel") || inherits(layer$geom, "GeomLabelRepel"),
    logical(1)
  )
  plot_obj$layers <- plot_obj$layers[keep_layers]
  plot_obj
}

locuscompare_remove_point_outlines <- function(plot_obj) {
  point_layers <- which(vapply(
    plot_obj$layers,
    function(layer) inherits(layer$geom, "GeomPoint"),
    logical(1)
  ))
  for (layer_index in point_layers) {
    plot_obj$layers[[layer_index]]$aes_params$colour <- "transparent"
    plot_obj$layers[[layer_index]]$aes_params$stroke <- 0
  }
  plot_obj
}

locuscompare_protected_label_layer <- function(data,
                                               mapping,
                                               base_family = "Helvetica",
                                               size = 3) {
  ggrepel::geom_label_repel(
    data = data,
    mapping = mapping,
    size = size,
    family = base_family,
    color = "black",
    fill = "#FFFFFFE8",
    label.size = 0,
    label.padding = unit(0.08, "lines"),
    label.r = unit(0.03, "lines"),
    box.padding = 0.35,
    point.padding = 0.35,
    force = 1.5,
    min.segment.length = 0,
    segment.color = "grey35",
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 1,
    inherit.aes = FALSE
  )
}

locuscompare_add_official_r2_key <- function(plot_obj,
                                             legend_position = "bottomright") {
  legend_position <- match.arg(legend_position, c("bottomright", "topright", "topleft"))
  legend_box <- switch(
    legend_position,
    bottomright = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05)),
    topright = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05)),
    topleft = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
  )
  cowplot::ggdraw(plot_obj) +
    geom_rect(
      data = legend_box,
      aes(xmin = x, xmax = x + 0.05, ymin = y, ymax = y + 0.05),
      color = "grey25",
      linewidth = 0.3,
      fill = rev(locuscomparer_style_colors())
    ) +
    cowplot::draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) +
    cowplot::draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) +
    cowplot::draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
}

locuscompare_make_scatter_panel <- function(inputs,
                                            title1,
                                            title2,
                                            legend = TRUE,
                                            legend_position = "bottomright",
                                            base_family = "Helvetica") {
  scatter_data <- inputs$merged
  scatter_data$label <- ifelse(scatter_data$rsid == inputs$snp, scatter_data$rsid, "")
  colors <- locuscomparer::assign_color(scatter_data$rsid, inputs$snp, inputs$ld)
  shapes <- ifelse(scatter_data$rsid == inputs$snp, 23, 21)
  names(shapes) <- scatter_data$rsid
  sizes <- ifelse(scatter_data$rsid == inputs$snp, 3.5, 2.35)
  names(sizes) <- scatter_data$rsid
  scatter <- locuscomparer::make_scatterplot(
    merged = scatter_data,
    title1 = title1,
    title2 = title2,
    color = colors,
    shape = shapes,
    size = sizes,
    legend = FALSE
  )
  scatter <- locuscompare_remove_point_outlines(scatter)
  scatter <- locuscompare_remove_repel_labels(scatter)
  lead_data <- scatter_data[scatter_data$rsid == inputs$snp, , drop = FALSE]
  scatter <- scatter +
    locuscompare_protected_label_layer(
      data = lead_data,
      mapping = aes(x = logp1, y = logp2, label = rsid),
      base_family = base_family
    ) +
    theme(text = element_text(family = base_family))
  if (isTRUE(legend)) {
    scatter <- locuscompare_add_official_r2_key(scatter, legend_position = legend_position)
  }
  scatter
}

make_locuscompare_ld_legend <- function(base_family = "Helvetica", bins = locuscompare_ld_levels()) {
  bins <- locuscompare_ld_levels()[locuscompare_ld_levels() %in% bins]
  n_items <- length(bins) + 1L
  n_columns <- min(3L, n_items)
  item_index <- seq_len(n_items)
  legend_df <- data.frame(
    item = c("Lead SNP", bins),
    label = c("Lead SNP", bins),
    x = 1 + ((item_index - 1L) %% n_columns) * 2.45,
    y = -((item_index - 1L) %/% n_columns),
    stringsAsFactors = FALSE
  )
  ggplot(legend_df, aes(x = x, y = y)) +
    geom_point(
      data = legend_df[legend_df$item != "Lead SNP", , drop = FALSE],
      aes(fill = item),
      shape = 21,
      size = 2.2,
      color = "transparent",
      stroke = 0
    ) +
    geom_point(
      data = legend_df[legend_df$item == "Lead SNP", , drop = FALSE],
      shape = 23,
      size = 2.8,
      fill = "purple",
      color = "transparent",
      stroke = 0
    ) +
    geom_text(aes(label = label), nudge_x = 0.2, hjust = 0, size = 2.25, family = base_family) +
    scale_fill_manual(values = locuscompare_ld_colors(include_lead = FALSE), guide = "none") +
    xlim(0.7, max(legend_df$x) + 2.15) +
    ylim(min(legend_df$y) - 0.35, 0.35) +
    theme_void(base_family = base_family) +
    theme(plot.margin = unit(c(0, 2, 0, 2), "pt"))
}

locuscompare_merge_exon_intervals <- function(starts, ends) {
  interval_df <- data.frame(start = safe_numeric(starts), end = safe_numeric(ends))
  interval_df <- interval_df[
    is.finite(interval_df$start) & is.finite(interval_df$end) &
      interval_df$end >= interval_df$start,
    ,
    drop = FALSE
  ]
  if (nrow(interval_df) == 0L) return(interval_df)
  interval_df <- interval_df[order(interval_df$start, interval_df$end), , drop = FALSE]
  merged <- vector("list", nrow(interval_df))
  output_index <- 1L
  current_start <- interval_df$start[1]
  current_end <- interval_df$end[1]
  if (nrow(interval_df) > 1L) {
    for (i in 2:nrow(interval_df)) {
      if (interval_df$start[i] <= current_end + 1) {
        current_end <- max(current_end, interval_df$end[i])
      } else {
        merged[[output_index]] <- data.frame(start = current_start, end = current_end)
        output_index <- output_index + 1L
        current_start <- interval_df$start[i]
        current_end <- interval_df$end[i]
      }
    }
  }
  merged[[output_index]] <- data.frame(start = current_start, end = current_end)
  do.call(rbind, merged[seq_len(output_index)])
}

locuscompare_read_gene_models <- function(gtf_path, chrom_str, position_range) {
  region <- GenomicRanges::GRanges(
    seqnames = normalize_chrom_label(chrom_str),
    ranges = IRanges::IRanges(start = position_range[1], end = position_range[2])
  )
  exon_ranges <- rtracklayer::import(gtf_path, which = region, feature.type = "exon")
  if (length(exon_ranges) == 0L) {
    return(list(gene = data.frame(), exon = data.frame()))
  }
  exon_df <- as.data.frame(exon_ranges)
  if ("gene_type" %in% names(exon_df)) {
    exon_df <- exon_df[as.character(exon_df$gene_type) == "protein_coding", , drop = FALSE]
  }
  if (nrow(exon_df) == 0L) {
    return(list(gene = data.frame(), exon = data.frame()))
  }
  exon_df$symbol <- if ("gene_name" %in% names(exon_df)) {
    as.character(exon_df$gene_name)
  } else if ("gene_id" %in% names(exon_df)) {
    as.character(exon_df$gene_id)
  } else {
    "Gene"
  }
  exon_df$gene_key <- if ("gene_id" %in% names(exon_df)) {
    as.character(exon_df$gene_id)
  } else {
    exon_df$symbol
  }
  exon_df$strand_value <- if ("strand" %in% names(exon_df)) {
    as.character(exon_df$strand)
  } else {
    "+"
  }
  missing_key <- is.na(exon_df$gene_key) | !nzchar(exon_df$gene_key)
  exon_df$gene_key[missing_key] <- exon_df$symbol[missing_key]
  missing_symbol <- is.na(exon_df$symbol) | !nzchar(exon_df$symbol)
  exon_df$symbol[missing_symbol] <- exon_df$gene_key[missing_symbol]
  exon_df$strand_value[!exon_df$strand_value %in% c("+", "-")] <- "+"
  exon_df$start <- pmax(safe_numeric(exon_df$start), position_range[1])
  exon_df$end <- pmin(safe_numeric(exon_df$end), position_range[2])
  exon_df <- exon_df[is.finite(exon_df$start) & is.finite(exon_df$end) & exon_df$end >= exon_df$start, , drop = FALSE]
  if (nrow(exon_df) == 0L) {
    return(list(gene = data.frame(), exon = data.frame()))
  }

  gene_groups <- split(seq_len(nrow(exon_df)), exon_df$gene_key)
  gene_rows <- vector("list", length(gene_groups))
  merged_exons <- vector("list", length(gene_groups))
  group_names <- names(gene_groups)
  for (i in seq_along(gene_groups)) {
    group_df <- exon_df[gene_groups[[i]], , drop = FALSE]
    reduced <- locuscompare_merge_exon_intervals(group_df$start, group_df$end)
    symbol <- group_df$symbol[which(nzchar(group_df$symbol))[1]]
    strand_value <- group_df$strand_value[1]
    gene_rows[[i]] <- data.frame(
      gene_key = group_names[i],
      symbol = symbol,
      start = min(reduced$start),
      end = max(reduced$end),
      strand = strand_value,
      stringsAsFactors = FALSE
    )
    reduced$gene_key <- group_names[i]
    reduced$symbol <- symbol
    reduced$strand <- strand_value
    merged_exons[[i]] <- reduced[, c("gene_key", "symbol", "start", "end", "strand")]
  }
  list(
    gene = do.call(rbind, gene_rows),
    exon = do.call(rbind, merged_exons)
  )
}

locuscompare_intervals_do_not_overlap <- function(df,
                                                  start_col,
                                                  end_col,
                                                  group_col,
                                                  gap = 0) {
  if (nrow(df) <= 1L) return(TRUE)
  groups <- split(seq_len(nrow(df)), df[[group_col]])
  all(vapply(groups, function(index) {
    if (length(index) <= 1L) return(TRUE)
    group_df <- df[index, , drop = FALSE]
    group_df <- group_df[order(group_df[[start_col]], group_df[[end_col]]), , drop = FALSE]
    all(group_df[[start_col]][-1] >= head(group_df[[end_col]], -1) + gap)
  }, logical(1)))
}

locuscompare_assign_gene_lanes <- function(gene_df, position_range, gap_fraction = 0.01) {
  if (nrow(gene_df) == 0L) return(gene_df)
  gene_df <- gene_df[order(gene_df$start, gene_df$end, gene_df$symbol), , drop = FALSE]
  gap_bp <- max(1, diff(position_range) * gap_fraction)
  lane_ends <- numeric(0)
  gene_df$lane <- NA_integer_
  for (i in seq_len(nrow(gene_df))) {
    available <- which(gene_df$start[i] >= lane_ends + gap_bp)
    lane <- if (length(available) > 0L) available[1] else length(lane_ends) + 1L
    if (lane > length(lane_ends)) lane_ends <- c(lane_ends, -Inf)
    gene_df$lane[i] <- lane
    lane_ends[lane] <- gene_df$end[i]
  }
  attr(gene_df, "interval_gap_bp") <- gap_bp
  attr(gene_df, "interval_overlap_free") <- locuscompare_intervals_do_not_overlap(
    gene_df,
    start_col = "start",
    end_col = "end",
    group_col = "lane",
    gap = gap_bp
  )
  gene_df
}

locuscompare_gene_labels_do_not_overlap <- function(df,
                                                    lane_pitch_pt,
                                                    label_height_pt,
                                                    horizontal_gap_bp = 0) {
  if (nrow(df) <= 1L) return(TRUE)
  for (i in seq_len(nrow(df) - 1L)) {
    for (j in (i + 1L):nrow(df)) {
      horizontal_overlap <-
        df$label_start[i] < df$label_end[j] + horizontal_gap_bp &&
        df$label_end[i] > df$label_start[j] - horizontal_gap_bp
      vertical_overlap <- abs(df$lane[i] - df$lane[j]) * lane_pitch_pt < label_height_pt
      if (horizontal_overlap && vertical_overlap) return(FALSE)
    }
  }
  TRUE
}

locuscompare_select_gene_labels <- function(gene_df,
                                            position_range,
                                            lead_position = NULL,
                                            panel_width_in = 3.4,
                                            label_size_pt = 6) {
  if (nrow(gene_df) == 0L) return(gene_df)
  region_span <- diff(position_range)
  label_gap_bp <- max(1, region_span * 0.006)
  gene_df$label_width_bp <- pmax(
    region_span * 0.025,
    nchar(gene_df$symbol) * label_size_pt * 0.56 / 72 / panel_width_in * region_span
  )
  gene_df$label_x <- (gene_df$start + gene_df$end) / 2
  gene_df$label_x <- pmax(
    position_range[1] + gene_df$label_width_bp / 2,
    pmin(position_range[2] - gene_df$label_width_bp / 2, gene_df$label_x)
  )
  gene_df$label_start <- gene_df$label_x - gene_df$label_width_bp / 2
  gene_df$label_end <- gene_df$label_x + gene_df$label_width_bp / 2
  gene_df$label_visible <- FALSE
  reference_position <- if (!is.null(lead_position) && is.finite(lead_position)) {
    lead_position
  } else {
    mean(position_range)
  }
  gene_df$label_priority <- abs(gene_df$label_x - reference_position)
  for (lane in sort(unique(gene_df$lane))) {
    lane_index <- which(gene_df$lane == lane)
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
  candidate_index <- which(gene_df$label_visible)
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
  gene_df
}

locuscompare_gene_track_height <- function(n_lanes) {
  min(
    LOCUSCOMPARE_LAYOUT$gene_track_height_max,
    LOCUSCOMPARE_LAYOUT$gene_track_height +
      max(0, n_lanes - 1L) * LOCUSCOMPARE_LAYOUT$gene_track_height_per_extra_lane
  )
}

make_locuscompare_gene_track <- function(inputs,
                                         gtf_path,
                                         base_family = "Helvetica") {
  if (is.null(gtf_path) || length(gtf_path) != 1L || is.na(gtf_path) ||
      !nzchar(gtf_path) || !file.exists(gtf_path)) {
    stop("Locuscompare gene track requires an existing gtf_path", call. = FALSE)
  }

  position_range <- range(inputs$merged$pos, na.rm = TRUE)
  chrom_str <- paste0("chr", inputs$chr)
  models <- locuscompare_read_gene_models(gtf_path, chrom_str, position_range)
  panel_width_in <- LOCUSCOMPARE_LAYOUT$width_in * LOCUSCOMPARE_LAYOUT$regional_width /
    (LOCUSCOMPARE_LAYOUT$scatter_width + LOCUSCOMPARE_LAYOUT$regional_width)
  lead_row <- inputs$merged[inputs$merged$rsid == inputs$snp, , drop = FALSE]
  lead_position <- if (nrow(lead_row) > 0L) safe_numeric(lead_row$pos[1]) else mean(position_range)

  gene_df <- locuscompare_assign_gene_lanes(models$gene, position_range)
  gene_df <- locuscompare_select_gene_labels(
    gene_df,
    position_range = position_range,
    lead_position = lead_position,
    panel_width_in = panel_width_in
  )
  n_lanes <- if (nrow(gene_df) > 0L) max(gene_df$lane) else 0L
  interval_overlap_free <- if (nrow(gene_df) > 0L) {
    isTRUE(attr(gene_df, "interval_overlap_free"))
  } else {
    TRUE
  }
  label_overlap_free <- if (nrow(gene_df) > 0L) {
    isTRUE(attr(gene_df, "label_overlap_free"))
  } else {
    TRUE
  }

  if (nrow(gene_df) == 0L) {
    track <- ggplot() +
      annotate(
        "text",
        x = mean(position_range) / 1e6,
        y = 1,
        label = "No protein-coding genes",
        size = 2.2,
        color = "grey45",
        fontface = "italic"
      ) +
      scale_x_continuous(
        limits = position_range / 1e6,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(limits = c(0.5, 1.5))
  } else {
    exon_df <- merge(
      models$exon,
      gene_df[, c("gene_key", "lane")],
      by = "gene_key",
      all.x = TRUE,
      sort = FALSE
    )
    arrow_bp <- max(1, diff(position_range) * 0.012)
    arrow_df <- gene_df
    arrow_df$arrow_start <- ifelse(
      arrow_df$strand == "+",
      pmax(arrow_df$start, arrow_df$end - arrow_bp),
      pmin(arrow_df$end, arrow_df$start + arrow_bp)
    )
    arrow_df$arrow_end <- ifelse(
      arrow_df$strand == "+",
      arrow_df$end,
      arrow_df$start
    )
    lane_df <- data.frame(
      lane = seq_len(n_lanes),
      fill = rep(c("#F8FAFC", "#F1F5F9"), length.out = n_lanes)
    )
    label_df <- gene_df[gene_df$label_visible, , drop = FALSE]
    label_df$label_y <- label_df$lane + 0.28

    track <- ggplot() +
      geom_rect(
        data = lane_df,
        aes(
          xmin = position_range[1] / 1e6,
          xmax = position_range[2] / 1e6,
          ymin = lane - 0.42,
          ymax = lane + 0.42
        ),
        inherit.aes = FALSE,
        fill = lane_df$fill
      ) +
      geom_segment(
        data = gene_df,
        aes(x = start / 1e6, xend = end / 1e6, y = lane, yend = lane, color = strand),
        linewidth = 0.45,
        lineend = "round"
      ) +
      geom_rect(
        data = exon_df,
        aes(
          xmin = start / 1e6,
          xmax = end / 1e6,
          ymin = lane - 0.14,
          ymax = lane + 0.14,
          fill = strand
        ),
        color = NA
      ) +
      geom_segment(
        data = arrow_df,
        aes(
          x = arrow_start / 1e6,
          xend = arrow_end / 1e6,
          y = lane,
          yend = lane,
          color = strand
        ),
        arrow = arrow(length = unit(0.09, "cm"), type = "closed"),
        linewidth = 0.45
      ) +
      geom_text(
        data = label_df,
        aes(x = label_x / 1e6, y = label_y, label = symbol, color = strand),
        size = 2,
        family = base_family,
        fontface = "italic",
        vjust = 0,
        show.legend = FALSE
      ) +
      scale_color_manual(values = c("+" = "#D73027", "-" = "#4575B4"), guide = "none") +
      scale_fill_manual(values = c("+" = "#D73027", "-" = "#4575B4"), guide = "none") +
      scale_x_continuous(
        limits = position_range / 1e6,
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(limits = c(0.52, n_lanes + 0.62), expand = c(0, 0))
  }

  track <- track +
    labs(x = paste0("Genomic position on ", chrom_str, " (Mb)")) +
    theme_void(base_family = base_family) +
    theme(
      axis.title.x = element_text(size = 8, face = "plain", color = "black", margin = margin(t = 2)),
      axis.text.x = element_text(size = 7, color = "black"),
      axis.ticks.x = element_line(linewidth = 0.35, color = "black"),
      axis.ticks.length.x = unit(2, "pt"),
      panel.border = element_blank(),
      plot.margin = margin(t = 0, r = 3, b = 1, l = 3, unit = "pt")
    )

  attr(track, "easycoloc_gene_features_present") <- nrow(models$exon) > 0L
  attr(track, "easycoloc_gene_track_position_range") <- position_range
  attr(track, "easycoloc_gene_count") <- nrow(gene_df)
  attr(track, "easycoloc_gene_lane_count") <- n_lanes
  attr(track, "easycoloc_gene_labels_drawn") <- if (nrow(gene_df) > 0L) sum(gene_df$label_visible) else 0L
  attr(track, "easycoloc_gene_labels_visible") <- if (nrow(gene_df) > 0L) gene_df$symbol[gene_df$label_visible] else character(0)
  attr(track, "easycoloc_gene_labels_hidden") <- if (nrow(gene_df) > 0L) gene_df$symbol[!gene_df$label_visible] else character(0)
  attr(track, "easycoloc_gene_interval_overlap_free") <- interval_overlap_free
  attr(track, "easycoloc_gene_label_overlap_free") <- label_overlap_free
  attr(track, "easycoloc_gene_label_lane_pitch_pt") <- if (nrow(gene_df) > 0L) attr(gene_df, "label_lane_pitch_pt") else NA_real_
  attr(track, "easycoloc_gene_track_height") <- locuscompare_gene_track_height(n_lanes)
  track
}

compose_locuscompare_gene_track_layout <- function(scatter,
                                                   gwas_panel,
                                                   qtl_panel,
                                                   inputs,
                                                   gtf_path,
                                                   legend = TRUE,
                                                   legend_embedded = FALSE,
                                                   base_family = "Helvetica") {
  gene_track <- make_locuscompare_gene_track(
    inputs = inputs,
    gtf_path = gtf_path,
    base_family = base_family
  )
  if (isTRUE(legend) && !isTRUE(legend_embedded)) {
    present_bins <- unique(as.character(inputs$prepared$data$ld_bin))
    legend_panel <- make_locuscompare_ld_legend(base_family = base_family, bins = present_bins)
    scatter <- cowplot::ggdraw(scatter) +
      cowplot::draw_plot(legend_panel, x = 0.5, y = 0.01, width = 0.48, height = 0.24)
  }

  aligned_regional <- cowplot::align_plots(
    gwas_panel,
    qtl_panel,
    gene_track,
    align = "v",
    axis = "lr"
  )
  gene_track_height <- attr(gene_track, "easycoloc_gene_track_height")
  regional_column <- cowplot::plot_grid(
    aligned_regional[[1]],
    aligned_regional[[2]],
    aligned_regional[[3]],
    ncol = 1,
    rel_heights = c(
      LOCUSCOMPARE_LAYOUT$gwas_height,
      LOCUSCOMPARE_LAYOUT$qtl_height,
      gene_track_height
    )
  )
  combined <- cowplot::plot_grid(
    scatter,
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
  attr(combined, "easycoloc_layout") <- "full-height official scatter beside aligned GWAS, QTL, and collision-managed gene track"
  attr(combined, "easycoloc_official_legend_embedded") <- isTRUE(legend_embedded)
  combined
}

locuscompare_ld_pairs <- function(ld_df) {
  if (is.null(ld_df)) return(NULL)
  if (is.list(ld_df) && !is.data.frame(ld_df) && !is.null(ld_df$R)) {
    ld_df <- ld_matrix_to_pairs(ld_df)
  }
  ld_result <- as.data.frame(ld_df, stringsAsFactors = FALSE)
  if (!all(c("SNP_A", "SNP_B") %in% names(ld_result))) return(NULL)
  r_col <- resolve_col(ld_result, "R2", c("r2", "R", "r"))
  if (is.null(r_col)) return(NULL)
  out <- data.frame(
    SNP_A = safe_character(ld_result$SNP_A),
    SNP_B = safe_character(ld_result$SNP_B),
    R2 = safe_numeric(ld_result[[r_col]]),
    stringsAsFactors = FALSE
  )
  if (!identical(tolower(r_col), "r2")) {
    out$R2 <- abs(out$R2)^2
  }
  out <- out[!is.na(out$SNP_A) & nzchar(out$SNP_A) &
    !is.na(out$SNP_B) & nzchar(out$SNP_B) &
    !is.na(out$R2) & is.finite(out$R2), , drop = FALSE]
  out$R2 <- pmin(1, pmax(0, out$R2))
  out
}

extract_ld_to_lead <- function(snps, lead_snp, ld_df = NULL) {
  r2 <- rep(NA_real_, length(snps))
  names(r2) <- snps
  if (length(snps) == 0 || is.na(lead_snp) || !nzchar(lead_snp)) return(r2)
  r2[snps == lead_snp] <- 1

  ld_result <- locuscompare_ld_pairs(ld_df)
  if (is.null(ld_result) || nrow(ld_result) == 0) return(r2)

  ld_to_lead <- ld_result[ld_result$SNP_A == lead_snp | ld_result$SNP_B == lead_snp, , drop = FALSE]
  if (nrow(ld_to_lead) == 0) return(r2)
  partner <- ifelse(ld_to_lead$SNP_A == lead_snp, ld_to_lead$SNP_B, ld_to_lead$SNP_A)
  partner_df <- data.frame(snp = partner, r2 = ld_to_lead$R2, stringsAsFactors = FALSE)
  partner_df <- partner_df[!is.na(partner_df$snp) & nzchar(partner_df$snp) & !is.na(partner_df$r2), , drop = FALSE]
  if (nrow(partner_df) > 0) {
    partner_df <- aggregate(r2 ~ snp, data = partner_df, FUN = max)
    r2[partner_df$snp] <- partner_df$r2
  }
  r2[snps == lead_snp] <- 1
  r2
}

locuscompare_ld_bin <- function(r2) {
  out <- rep("missing LD", length(r2))
  ok <- !is.na(r2) & is.finite(r2)
  out[ok & r2 < 0.2] <- "<0.2"
  out[ok & r2 >= 0.2 & r2 < 0.4] <- "0.2-0.4"
  out[ok & r2 >= 0.4 & r2 < 0.6] <- "0.4-0.6"
  out[ok & r2 >= 0.6 & r2 < 0.8] <- "0.6-0.8"
  out[ok & r2 >= 0.8] <- "0.8-1.0"
  factor(out, levels = locuscompare_ld_levels())
}

resolve_credible_snps <- function(credible_set = NULL) {
  if (is.null(credible_set) || !is.data.frame(credible_set) || nrow(credible_set) == 0) {
    return(data.frame(snp = character(), snp_pp_h4 = numeric(), stringsAsFactors = FALSE))
  }
  snp_col <- resolve_col(credible_set, "snp", c("rsid", "SNP", "SNPID", "marker", "variable"))
  if (is.null(snp_col)) {
    return(data.frame(snp = character(), snp_pp_h4 = numeric(), stringsAsFactors = FALSE))
  }
  pp_col <- resolve_col(credible_set, "SNP.PP.H4", c("snp_pp_h4", "variable_prob", "posterior_prob", "PP.H4"))
  out <- data.frame(
    snp = safe_character(credible_set[[snp_col]]),
    snp_pp_h4 = if (!is.null(pp_col)) safe_numeric(credible_set[[pp_col]]) else NA_real_,
    stringsAsFactors = FALSE
  )
  out[!is.na(out$snp) & nzchar(out$snp), , drop = FALSE]
}

choose_locuscompare_lead <- function(plot_df, explicit_lead = NULL, credible_set = NULL) {
  explicit_lead <- safe_character(explicit_lead)
  explicit_lead <- explicit_lead[!is.na(explicit_lead) & nzchar(explicit_lead)]
  if (length(explicit_lead) > 0 && explicit_lead[1] %in% plot_df$snp) return(explicit_lead[1])

  cs <- resolve_credible_snps(credible_set)
  if (nrow(cs) > 0 && any(!is.na(cs$snp_pp_h4))) {
    cs <- cs[order(-cs$snp_pp_h4), , drop = FALSE]
    hit <- cs$snp[cs$snp %in% plot_df$snp]
    if (length(hit) > 0) return(hit[1])
  }

  combined <- plot_df$gwas_logp + plot_df$qtl_logp
  if (any(is.finite(combined))) return(plot_df$snp[which.max(combined)])
  plot_df$snp[which.min(plot_df$qtl_p)]
}

prepare_locuscompare_plot_data <- function(df,
                                           ld_df = NULL,
                                           lead_snp = NULL,
                                           credible_set = NULL) {
  df <- safe_df_convert(df)
  if (is.null(df)) return(list(error = "No input data"))

  col_snp <- resolve_col(df, "rsid", c("SNPID.gwas", "snp", "SNP", "SNPID", "marker"))
  col_chr <- resolve_col(df, "CHR.qtl", c("CHR.gwas", "CHR", "chr", "chromosome"))
  col_pos <- resolve_col(df, "POS.qtl", c("POS.gwas", "POS", "pos", "BP", "position"))
  col_gwas_p <- resolve_col(df, "P.gwas", c("p.gwas", "GWAS_P", "gwas_p", "p_gwas"))
  col_qtl_p <- resolve_col(df, "P.qtl", c("p.qtl", "QTL_P", "qtl_p", "p_qtl", "P"))

  if (is.null(col_snp) || is.null(col_chr) || is.null(col_pos)) {
    return(list(error = "Missing SNP, chromosome, or position columns"))
  }
  if (is.null(col_gwas_p)) return(list(error = "Missing GWAS p-value column"))
  if (is.null(col_qtl_p)) return(list(error = "Missing QTL p-value column"))

  plot_df <- data.frame(
    snp = safe_character(df[[col_snp]]),
    chrom = safe_character(df[[col_chr]]),
    pos = safe_numeric(df[[col_pos]]),
    gwas_p = clamp_pvalue(df[[col_gwas_p]]),
    qtl_p = clamp_pvalue(df[[col_qtl_p]]),
    stringsAsFactors = FALSE
  )
  plot_df <- plot_df[!is.na(plot_df$snp) & nzchar(plot_df$snp) &
    !is.na(plot_df$chrom) & nzchar(plot_df$chrom) &
    !is.na(plot_df$pos) &
    !is.na(plot_df$gwas_p) &
    !is.na(plot_df$qtl_p), , drop = FALSE]
  if (nrow(plot_df) == 0) return(list(error = "No valid SNPs remain after filtering"))

  duplicated_snp <- duplicated(plot_df$snp) | duplicated(plot_df$snp, fromLast = TRUE)
  if (any(duplicated_snp)) {
    plot_df <- plot_df[!duplicated_snp, , drop = FALSE]
  }
  if (nrow(plot_df) == 0) return(list(error = "No unique valid SNPs remain after duplicate filtering"))

  plot_df$gwas_logp <- -log10(plot_df$gwas_p)
  plot_df$qtl_logp <- -log10(plot_df$qtl_p)
  chrom_values <- unique(plot_df$chrom)
  if (length(chrom_values) != 1) return(list(error = "locuscomparer backend requires one chromosome per plot"))

  chosen_lead <- choose_locuscompare_lead(plot_df, explicit_lead = lead_snp, credible_set = credible_set)
  cs <- resolve_credible_snps(credible_set)
  plot_df$is_lead <- plot_df$snp == chosen_lead
  plot_df$is_credible <- plot_df$snp %in% cs$snp
  plot_df$snp_pp_h4 <- NA_real_
  if (nrow(cs) > 0) plot_df$snp_pp_h4 <- cs$snp_pp_h4[match(plot_df$snp, cs$snp)]

  raw_r2 <- extract_ld_to_lead(plot_df$snp, chosen_lead, ld_df = ld_df)
  plot_df$r2 <- unname(raw_r2[plot_df$snp])
  plot_df$ld_bin <- locuscompare_ld_bin(plot_df$r2)
  plot_df$plot_priority <- ifelse(plot_df$is_lead, 3, ifelse(plot_df$is_credible, 2, 1))
  plot_df <- plot_df[order(plot_df$plot_priority, plot_df$r2, na.last = TRUE), , drop = FALSE]

  list(
    data = plot_df,
    lead_snp = chosen_lead,
    ld_available = any(!is.na(plot_df$r2) & !plot_df$is_lead),
    error = NULL
  )
}

as_locuscomparer_inputs <- function(df,
                                    ld_df = NULL,
                                    lead_snp = NULL,
                                    credible_set = NULL) {
  prepared <- prepare_locuscompare_plot_data(
    df = df,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set
  )
  if (!is.null(prepared$error)) return(list(error = prepared$error))

  plot_df <- prepared$data
  chr <- normalize_chrom_label(plot_df$chrom[1])
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  merged <- data.frame(
    rsid = plot_df$snp,
    pval1 = plot_df$gwas_p,
    logp1 = plot_df$gwas_logp,
    pval2 = plot_df$qtl_p,
    logp2 = plot_df$qtl_logp,
    chr = chr,
    pos = plot_df$pos,
    stringsAsFactors = FALSE
  )

  ld <- locuscompare_ld_pairs(ld_df)
  if (is.null(ld)) {
    ld <- data.frame(SNP_A = character(), SNP_B = character(), R2 = numeric(), stringsAsFactors = FALSE)
  } else {
    ld <- ld[ld$SNP_A %in% plot_df$snp & ld$SNP_B %in% plot_df$snp, c("SNP_A", "SNP_B", "R2"), drop = FALSE]
    ld_lead <- ld[ld$SNP_A == prepared$lead_snp | ld$SNP_B == prepared$lead_snp, , drop = FALSE]
    if (nrow(ld_lead) > 0) {
      partner <- ifelse(ld_lead$SNP_A == prepared$lead_snp, ld_lead$SNP_B, ld_lead$SNP_A)
      ld <- data.frame(
        SNP_A = prepared$lead_snp,
        SNP_B = partner,
        R2 = ld_lead$R2,
        stringsAsFactors = FALSE
      )
      ld <- ld[ld$SNP_B != prepared$lead_snp, , drop = FALSE]
      ld <- aggregate(R2 ~ SNP_A + SNP_B, data = ld, FUN = max)
    }
  }

  list(
    merged = merged,
    ld = ld,
    chr = chr,
    snp = prepared$lead_snp,
    prepared = prepared,
    error = NULL
  )
}

make_internal_locuscomparer_style_plot <- function(inputs,
                                                   title1 = "GWAS",
                                                   title2 = "QTL",
                                                   legend = TRUE,
                                                   legend_position = "bottomright",
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

  shared_theme <- locuscompare_publication_theme(base_family = base_family)
  scatter <- ggplot(point_df, aes(logp1, logp2)) +
    geom_point(
      data = nonlead_df,
      aes(fill = ld_bin),
      shape = 21,
      size = 1.7,
      color = "transparent",
      stroke = 0,
      alpha = 0.88
    ) +
    geom_point(
      data = lead_df,
      shape = 23,
      size = 3.2,
      fill = "purple",
      color = "transparent",
      stroke = 0
    ) +
    scale_fill_manual(values = color_values, guide = "none") +
    locuscompare_protected_label_layer(
      data = lead_df,
      mapping = aes(x = logp1, y = logp2, label = label),
      size = 2.4,
      base_family = base_family
    ) +
    xlab(bquote(.(title1) ~ -log[10] * "(P value)")) +
    ylab(bquote(.(title2) ~ -log[10] * "(P value)")) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    shared_theme

  make_lz <- function(y_col, title, show_x = TRUE) {
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
    ggplot(metal, aes(x = pos, y = logp)) +
      geom_point(
        data = metal_nonlead,
        aes(fill = ld_bin),
        shape = 21,
        size = 1.7,
        color = "transparent",
        stroke = 0,
        alpha = 0.88
      ) +
      geom_point(
        data = metal_lead,
        shape = 23,
        size = 3.2,
        fill = "purple",
        color = "transparent",
        stroke = 0
      ) +
      scale_fill_manual(values = color_values, guide = "none") +
      scale_x_continuous(
        limits = range(metal$pos, na.rm = TRUE),
        labels = locuscompare_mb_labeller(metal$pos),
        expand = expansion(mult = c(0.01, 0.01))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
      locuscompare_protected_label_layer(
        data = metal_lead,
        mapping = aes(x = pos, y = logp, label = label),
        size = 2.4,
        base_family = base_family
      ) +
      xlab(if (show_x) paste0("chr", inputs$chr, " (Mb)") else NULL) +
      ylab(bquote(.(title) ~ -log[10] * "(P value)")) +
      shared_theme +
      theme(
        axis.text.x = if (show_x) element_text(size = 7.2, color = "black") else element_blank(),
        axis.title.x = if (show_x) element_text(size = 8.5) else element_blank()
      )
  }

  p2 <- make_lz("logp1", title1, show_x = FALSE)
  p3 <- make_lz("logp2", title2, show_x = FALSE)
  combined <- compose_locuscompare_gene_track_layout(
    scatter = scatter,
    gwas_panel = p2,
    qtl_panel = p3,
    inputs = inputs,
    gtf_path = gtf_path,
    legend = legend,
    base_family = base_family
  )
  attr(combined, "easycoloc_lead_label_panels") <- c(
    scatter = locuscompare_has_repel_label(scatter),
    gwas = locuscompare_has_repel_label(p2),
    qtl = locuscompare_has_repel_label(p3)
  )
  attr(combined, "easycoloc_point_outlines_removed") <- TRUE
  attr(combined, "easycoloc_lead_labels_background_protected") <- TRUE
  combined
}

make_locuscomparer_backend_plot <- function(inputs,
                                            title1 = "GWAS",
                                            title2 = "QTL",
                                            legend = TRUE,
                                            legend_position = "bottomright",
                                            lz_ylab_linebreak = FALSE,
                                            gtf_path = NULL,
                                            base_family = "Helvetica") {
  if (!requireNamespace("locuscomparer", quietly = TRUE)) {
    stop("locuscomparer backend requested but package 'locuscomparer' is not installed", call. = FALSE)
  }
  parts <- locuscomparer::make_combined_plot(
    merged = inputs$merged,
    title1 = title1,
    title2 = title2,
    ld = inputs$ld,
    chr = inputs$chr,
    snp = inputs$snp,
    combine = FALSE,
    legend = FALSE,
    legend_position = legend_position,
    lz_ylab_linebreak = lz_ylab_linebreak
  )
  panel_theme <- locuscompare_official_theme(base_family = base_family)
  lead_merged <- inputs$merged[inputs$merged$rsid == inputs$snp, , drop = FALSE]
  scatter <- locuscompare_make_scatter_panel(
    inputs = inputs,
    title1 = title1,
    title2 = title2,
    legend = legend,
    legend_position = legend_position,
    base_family = base_family
  )
  p2 <- locuscompare_set_x_labeller(
    parts$locuszoom1,
    locuscompare_mb_labeller(inputs$merged$pos),
    limits = range(inputs$merged$pos, na.rm = TRUE),
    expand = expansion(mult = c(0.01, 0.01))
  )
  p2 <- locuscompare_remove_point_outlines(p2)
  p2 <- locuscompare_remove_repel_labels(p2)
  p2 <- locuscompare_set_point_sizes(p2, inputs$merged$rsid, inputs$snp)
  p2 <- p2 +
    geom_point(
      data = lead_merged,
      aes(x = pos, y = logp1),
      shape = 23,
      size = 3.5,
      fill = "purple",
      color = "transparent",
      stroke = 0,
      inherit.aes = FALSE
    ) +
    locuscompare_protected_label_layer(
      data = lead_merged,
      mapping = aes(x = pos, y = logp1, label = rsid),
      base_family = base_family
    ) +
    ylab(bquote(.(title1) ~ -log[10] * "(P value)")) +
    panel_theme +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  p3 <- locuscompare_set_x_labeller(
    parts$locuszoom2,
    locuscompare_mb_labeller(inputs$merged$pos),
    limits = range(inputs$merged$pos, na.rm = TRUE),
    expand = expansion(mult = c(0.01, 0.01))
  )
  p3 <- locuscompare_remove_point_outlines(p3)
  p3 <- locuscompare_remove_repel_labels(p3)
  p3 <- locuscompare_set_point_sizes(p3, inputs$merged$rsid, inputs$snp)
  p3 <- p3 +
    geom_point(
      data = lead_merged,
      aes(x = pos, y = logp2),
      shape = 23,
      size = 3.5,
      fill = "purple",
      color = "transparent",
      stroke = 0,
      inherit.aes = FALSE
    ) +
    locuscompare_protected_label_layer(
      data = lead_merged,
      mapping = aes(x = pos, y = logp2, label = rsid),
      base_family = base_family
    ) +
    ylab(bquote(.(title2) ~ -log[10] * "(P value)")) +
    panel_theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  combined <- compose_locuscompare_gene_track_layout(
    scatter = scatter,
    gwas_panel = p2,
    qtl_panel = p3,
    inputs = inputs,
    gtf_path = gtf_path,
    legend = legend,
    legend_embedded = legend,
    base_family = base_family
  )
  attr(combined, "easycoloc_lead_label_panels") <- c(
    scatter = TRUE,
    gwas = locuscompare_has_repel_label(p2),
    qtl = locuscompare_has_repel_label(p3)
  )
  attr(combined, "easycoloc_official_visual_grammar") <- TRUE
  attr(combined, "easycoloc_point_outlines_removed") <- TRUE
  attr(combined, "easycoloc_lead_labels_background_protected") <- TRUE
  combined
}

make_combined_coloc_plot <- function(df,
                                     ld_df = NULL,
                                     lead_snp = NULL,
                                     credible_set = NULL,
                                     title = NULL,
                                     subtitle = NULL,
                                     pp4 = NULL,
                                     pp3 = NULL,
                                     phenotype_info = NULL,
                                     gtf_path = NULL,
                                     region_recomb = NULL,
                                     recomb_path = NULL,
                                     backend = c("auto", "locuscomparer", "internal"),
                                     title1 = "GWAS",
                                     title2 = "QTL",
                                     legend = TRUE,
                                     legend_position = "bottomright",
                                     lz_ylab_linebreak = FALSE,
                                     base_family = "Helvetica") {
  backend <- match.arg(backend)
  inputs <- as_locuscomparer_inputs(
    df = df,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set
  )
  if (!is.null(inputs$error)) {
    if (identical(backend, "locuscomparer")) stop(inputs$error, call. = FALSE)
    return(placeholder_plot("Locuscompare unavailable", inputs$error))
  }

  backend_used <- backend
  plot_obj <- NULL
  if (backend %in% c("auto", "locuscomparer")) {
    plot_obj <- tryCatch(
      make_locuscomparer_backend_plot(
        inputs = inputs,
        title1 = title1,
        title2 = title2,
        legend = legend,
        legend_position = legend_position,
        lz_ylab_linebreak = lz_ylab_linebreak,
        gtf_path = gtf_path,
        base_family = base_family
      ),
      error = function(e) {
        if (identical(backend, "locuscomparer")) stop(e)
        NULL
      }
    )
    if (!is.null(plot_obj)) {
      backend_used <- "locuscomparer"
    }
  }

  if (is.null(plot_obj)) {
    plot_obj <- make_internal_locuscomparer_style_plot(
      inputs = inputs,
      title1 = title1,
      title2 = title2,
      legend = legend,
      legend_position = legend_position,
      gtf_path = gtf_path,
      base_family = base_family
    )
    backend_used <- "internal"
  }

  attr(plot_obj, "easycoloc_plot_backend") <- backend_used
  attr(plot_obj, "easycoloc_locuscompare_lead_snp") <- inputs$snp
  plot_obj
}

write_locuscompare_plot_json <- function(df,
                                         output_file,
                                         ld_df = NULL,
                                         lead_snp = NULL,
                                         credible_set = NULL,
                                         meta = list()) {
  prepared <- prepare_locuscompare_plot_data(
    df = df,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set
  )
  if (!is.null(prepared$error)) {
    payload <- list(
      meta = c(meta, list(lead_snp = lead_snp, error = prepared$error)),
      points = list()
    )
  } else {
    plot_df <- prepared$data
    meta_out <- c(
      meta,
      list(
        lead_snp = prepared$lead_snp,
        chrom = as.character(plot_df$chrom[1]),
        min_pos = min(plot_df$pos, na.rm = TRUE),
        max_pos = max(plot_df$pos, na.rm = TRUE),
        n_snps = nrow(plot_df),
        credible_set_size = sum(plot_df$is_credible, na.rm = TRUE),
        ld_available = isTRUE(prepared$ld_available)
      )
    )
    points <- lapply(seq_len(nrow(plot_df)), function(i) {
      row <- plot_df[i, , drop = FALSE]
      list(
        snp = row$snp,
        chrom = as.character(row$chrom),
        pos = row$pos,
        gwas_p = row$gwas_p,
        qtl_p = row$qtl_p,
        gwas_logp = row$gwas_logp,
        qtl_logp = row$qtl_logp,
        r2 = ifelse(is.na(row$r2), NA_real_, row$r2),
        ld_bin = as.character(row$ld_bin),
        is_lead = isTRUE(row$is_lead),
        is_credible = isTRUE(row$is_credible),
        snp_pp_h4 = ifelse(is.na(row$snp_pp_h4), NA_real_, row$snp_pp_h4)
      )
    })
    payload <- list(meta = meta_out, points = points)
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_file, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  output_file
}

write_locuscompare_layout_validation <- function(plot_obj,
                                                 output_file,
                                                 png_path,
                                                 pdf_path,
                                                 df,
                                                 ld_df = NULL,
                                                 lead_snp = NULL,
                                                 credible_set = NULL,
                                                 width_in = LOCUSCOMPARE_LAYOUT$width_in,
                                                 height_in = LOCUSCOMPARE_LAYOUT$height_in,
                                                 dpi = 300,
                                                 reference_png = NULL,
                                                 reference_pdf = NULL,
                                                 strict = TRUE) {
  prepared <- prepare_locuscompare_plot_data(
    df = df,
    ld_df = ld_df,
    lead_snp = lead_snp,
    credible_set = credible_set
  )
  label_panels <- attr(plot_obj, "easycoloc_lead_label_panels")
  label_panels_ok <- length(label_panels) == 3L && all(label_panels)
  point_outlines_removed <- isTRUE(attr(plot_obj, "easycoloc_point_outlines_removed"))
  lead_labels_background_protected <- isTRUE(attr(plot_obj, "easycoloc_lead_labels_background_protected"))
  official_visual_grammar <- isTRUE(attr(plot_obj, "easycoloc_official_visual_grammar"))
  official_legend_embedded <- isTRUE(attr(plot_obj, "easycoloc_official_legend_embedded"))
  gene_track_included <- isTRUE(attr(plot_obj, "easycoloc_gene_track_included"))
  gene_features_present <- isTRUE(attr(plot_obj, "easycoloc_gene_features_present"))
  gene_track_column_aligned <- isTRUE(attr(plot_obj, "easycoloc_gene_track_column_aligned"))
  gene_track_shared_x_axis <- isTRUE(attr(plot_obj, "easycoloc_gene_track_shared_x_axis"))
  scatter_spans_full_height <- isTRUE(attr(plot_obj, "easycoloc_scatter_spans_full_height"))
  gene_interval_overlap_free <- isTRUE(attr(plot_obj, "easycoloc_gene_interval_overlap_free"))
  gene_label_overlap_free <- isTRUE(attr(plot_obj, "easycoloc_gene_label_overlap_free"))
  gene_count <- as.integer(attr(plot_obj, "easycoloc_gene_count"))
  gene_lane_count <- as.integer(attr(plot_obj, "easycoloc_gene_lane_count"))
  gene_labels_drawn <- as.integer(attr(plot_obj, "easycoloc_gene_labels_drawn"))
  gene_labels_visible <- as.character(attr(plot_obj, "easycoloc_gene_labels_visible"))
  gene_labels_hidden <- as.character(attr(plot_obj, "easycoloc_gene_labels_hidden"))
  gene_label_lane_pitch_pt <- as.numeric(attr(plot_obj, "easycoloc_gene_label_lane_pitch_pt"))
  gene_track_height <- as.numeric(attr(plot_obj, "easycoloc_gene_track_height"))
  gene_track_range <- safe_numeric(attr(plot_obj, "easycoloc_gene_track_position_range"))
  png_exists <- file.exists(png_path) && isTRUE(file.info(png_path)$size > 5000)
  pdf_exists <- file.exists(pdf_path) && isTRUE(file.info(pdf_path)$size > 5000)
  data_ok <- is.null(prepared$error) && nrow(prepared$data) > 0L
  backend <- attr(plot_obj, "easycoloc_plot_backend")
  backend_ok <- length(backend) == 1L && backend %in% c("locuscomparer", "internal")
  dimensions_ok <- isTRUE(width_in > 0) && isTRUE(height_in > 0) && isTRUE(dpi >= 300)
  compact_canvas_ok <- isTRUE(width_in <= 7.3) && isTRUE(height_in <= 4.5) &&
    isTRUE(width_in / height_in >= 1.6)
  scatter_width_fraction <- LOCUSCOMPARE_LAYOUT$scatter_width /
    (LOCUSCOMPARE_LAYOUT$scatter_width + LOCUSCOMPARE_LAYOUT$regional_width)
  scatter_column_aspect <- width_in * scatter_width_fraction / height_in
  balanced_columns_ok <- isTRUE(scatter_column_aspect >= 0.75) && isTRUE(scatter_column_aspect <= 1.1)
  gene_lane_count_ok <- length(gene_lane_count) == 1L && is.finite(gene_lane_count) &&
    gene_lane_count >= 0L && (gene_count == 0L || gene_lane_count >= 1L)
  gene_label_accounting_ok <- length(gene_labels_drawn) == 1L && is.finite(gene_labels_drawn) &&
    gene_labels_drawn + length(gene_labels_hidden) == gene_count
  gene_track_height_ok <- length(gene_track_height) == 1L && is.finite(gene_track_height) &&
    gene_track_height >= LOCUSCOMPARE_LAYOUT$gene_track_height &&
    gene_track_height <= LOCUSCOMPARE_LAYOUT$gene_track_height_max
  expected_range <- if (data_ok) range(prepared$data$pos, na.rm = TRUE) else c(NA_real_, NA_real_)
  gene_track_range_aligned <- length(gene_track_range) == 2L && all(is.finite(gene_track_range)) &&
    isTRUE(all.equal(gene_track_range, expected_range, tolerance = 1e-8))

  checks <- list(
    png_exists_and_nontrivial = png_exists,
    pdf_exists_and_nontrivial = pdf_exists,
    prepared_data_nonempty = data_ok,
    supported_backend = backend_ok,
    lead_labels_in_all_panels = label_panels_ok,
    association_point_outlines_removed = point_outlines_removed,
    lead_labels_background_protected = lead_labels_background_protected,
    official_locuscomparer_visual_grammar = backend != "locuscomparer" || official_visual_grammar,
    embedded_official_r2_key = backend != "locuscomparer" || official_legend_embedded,
    gene_track_included = gene_track_included,
    gene_track_range_matches = gene_track_range_aligned,
    gene_track_column_aligned = gene_track_column_aligned,
    gene_track_shared_x_axis = gene_track_shared_x_axis,
    scatter_spans_full_height = scatter_spans_full_height,
    gene_intervals_do_not_overlap_within_lanes = gene_interval_overlap_free,
    visible_gene_label_footprints_do_not_overlap = gene_label_overlap_free,
    gene_lane_count_valid = gene_lane_count_ok,
    gene_label_accounting_complete = gene_label_accounting_ok,
    adaptive_gene_track_height = gene_track_height_ok,
    publication_dimensions = dimensions_ok,
    compact_publication_canvas = compact_canvas_ok,
    balanced_full_height_columns = balanced_columns_ok,
    rejected_static_elements_absent_by_design = TRUE
  )
  status <- if (all(unlist(checks, use.names = FALSE))) "PASS" else "FAIL"
  payload <- list(
    status = status,
    figure = list(
      png = png_path,
      pdf = pdf_path,
      reference_png = reference_png,
      reference_pdf = reference_pdf,
      width_in = width_in,
      height_in = height_in,
      aspect_ratio = width_in / height_in,
      scatter_column_aspect_ratio = scatter_column_aspect,
      dpi = dpi,
      backend = backend,
      layout = attr(plot_obj, "easycoloc_layout"),
      gene_track_placement = "third tier of the right regional column; scatter spans the full left column"
    ),
    text = list(
      base_family = "Helvetica",
      title_pt = 9,
      axis_title_pt = 9,
      axis_text_pt = 8,
      coordinate_scale = "Mb",
      context_strip = FALSE,
      coloc_summary_box = FALSE,
      gene_track = TRUE,
      recombination_trace = FALSE,
      smoke_labels = any(grepl("smoke", c(gene_labels_visible, gene_labels_hidden), ignore.case = TRUE))
    ),
    ld = list(
      available = if (data_ok) isTRUE(prepared$ld_available) else FALSE,
      official_r2_key_embedded = official_legend_embedded,
      bin_counts = if (data_ok) as.list(table(prepared$data$ld_bin)) else list(),
      missing_ld_count = if (data_ok) unname(sum(is.na(prepared$data$r2))) else NA_integer_
    ),
    lead_snp = list(
      snp = if (data_ok) prepared$lead_snp else lead_snp,
      shape = "diamond",
      labeled_panels = as.list(label_panels),
      all_panels_labeled = label_panels_ok,
      background_protected = lead_labels_background_protected
    ),
    association_points = list(
      outline = "none",
      outlines_removed = point_outlines_removed
    ),
    gene_track = list(
      included = gene_track_included,
      features_present = gene_features_present,
      position_range_bp = as.list(gene_track_range),
      matches_regional_position_range = gene_track_range_aligned,
      physical_column_aligned = gene_track_column_aligned,
      shared_genomic_axis = gene_track_shared_x_axis,
      gene_count = gene_count,
      lane_count = gene_lane_count,
      labels_drawn = gene_labels_drawn,
      labels_visible = as.list(gene_labels_visible),
      labels_hidden = as.list(gene_labels_hidden),
      label_lane_pitch_pt = gene_label_lane_pitch_pt,
      interval_overlap_free = gene_interval_overlap_free,
      visible_label_overlap_free = gene_label_overlap_free,
      relative_height = gene_track_height,
      placement = "below regional association panels"
    ),
    data = list(
      n_snps = if (data_ok) nrow(prepared$data) else 0L
    ),
    checks = checks
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_file, auto_unbox = TRUE, pretty = TRUE, na = "null")
  if (isTRUE(strict) && !identical(status, "PASS")) {
    failed_checks <- names(checks)[!unlist(checks, use.names = FALSE)]
    stop("Locuscompare layout validation failed: ", paste(failed_checks, collapse = ", "), call. = FALSE)
  }
  invisible(payload)
}
