#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
# Figure 3 panel C: prenatal vs postnatal SuSiE credible-set overlap per gene.
# Gene labels sit in fixed left/right lanes outside the coordinate range,
# on the exact chromosome row, joined to overlap centres by dashed leaders.
suppressPackageStartupMessages({
    library(arrow)
    library(data.table)
    library(ggplot2)
    library(scales)
    library(rtracklayer)
})

# --- config ---

qtl_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
fig_dir <- file.path(Sys.getenv("APA_ROOT"), "figures", "04_qtl")
layout_validation_path <- file.path(qtl_dir, "figure3_panelC_layout_validation.tsv")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- list(
    prenatal_susie = file.path(qtl_dir, "susie", "prenatal.SuSiE_summary.parquet"),
    ds_susie_glob  = file.path(qtl_dir, "downsample", "seed_*", "qtl", "Postnatal", "susie", "*postnatal.SuSiE_summary.parquet"),
    gtf            = file.path(Sys.getenv("APA_ROOT"), "ref", "gencode.v40.annotation.gtf.gz")
)

panel_width <- 12.0
panel_height <- 7.2
dpi <- 320
base_size <- 10.5
min_block_width_mb <- 0.55
top_overlap_label_n <- 12L
max_labels_per_chr <- 1L
row_spacing <- 1.35
left_label_x <- -15
right_label_x <- 265

# --- helper functions ---

log_msg <- function(message) {
    cat(sprintf("[%s] %s\n", Sys.time(), message))
}

add_pair_keys <- function(dt) {
    dt <- data.table::as.data.table(data.table::copy(dt))

    # gene = second '|' field of the phenotype id, else the whole id
    phenotype_id <- as.character(dt$phenotype_id)
    parts <- data.table::tstrsplit(phenotype_id, "|", fixed = TRUE, fill = NA_character_)
    raw_genes <- parts[[2]]
    bad <- is.na(raw_genes) | raw_genes == ""
    raw_genes[bad] <- phenotype_id[bad]

    cleaned_genes <- sub("_.*$", "", sub("\\.[0-9]+$", "", raw_genes))
    is_ensg <- startsWith(cleaned_genes, "ENSG")
    translated_names <- id2name[cleaned_genes]
    dt[, gene_key := fifelse(is_ensg & !is.na(translated_names), translated_names, cleaned_genes)]

    # variant key = chr:pos:ref:alt, build suffix _b<N> dropped
    clean_ids <- as.character(dt$variant_id)
    clean_ids <- sub("_b[0-9]+$", "", clean_ids)
    m <- regexec("^(chr[^:_]+)[_:]([0-9]+)[_:]([^_:]+)[_:]([^_:]+)$", clean_ids)
    mm <- regmatches(clean_ids, m)
    ok <- lengths(mm) == 5
    if (any(ok)) {
        clean_ids[ok] <- vapply(mm[ok], \(z) paste(z[2], z[3], z[4], z[5], sep = ":"), character(1))
    }
    dt[, variant_key := clean_ids]

    dt[, pair_key := paste(gene_key, variant_key, sep = "||")]
    dt[]
}

cs_intervals <- function(cs_dt, stage_label) {
    parsed <- tstrsplit(cs_dt$variant_key, ":", fixed = TRUE)
    cs_dt <- cbind(cs_dt, chr = parsed[[1]], pos = as.numeric(parsed[[2]]))
    cs_dt <- cs_dt[is.finite(pos)]
    cs_dt[, .(
        start = min(pos), end = max(pos),
        center = mean(pos), max_pip = max(pip, na.rm = TRUE),
        n_variants = .N, stage = stage_label
    ), by = .(gene_key, cs_id, chr)]
}

# --- loading data ---

log_msg("Loading gene annotation")
gtf_data <- rtracklayer::import(inputs$gtf)
gene_idx <- gtf_data$type == "gene"
gtf_dt <- data.table::as.data.table(mcols(gtf_data)[gene_idx, c("gene_id", "gene_name")])
gtf_dt[, gene_id_clean := sub("\\.[0-9]+$", "", as.character(gene_id))]
id2name <- setNames(as.character(gtf_dt$gene_name), gtf_dt$gene_id_clean)

log_msg("Loading prenatal SuSiE")
prenatal_susie <- add_pair_keys(as.data.table(read_parquet(inputs$prenatal_susie)))

log_msg("Loading postnatal downsample SuSiE")
ds_susie_paths <- Sys.glob(inputs$ds_susie_glob)
postnatal_ds_susie <- rbindlist(lapply(ds_susie_paths, function(path) {
    dt <- add_pair_keys(as.data.table(read_parquet(path)))
    dt[, .(pip = max(pip, na.rm = TRUE)), by = .(gene_key, variant_key)]
}), fill = TRUE)
postnatal_ds_susie <- postnatal_ds_susie[, .(pip = max(pip, na.rm = TRUE)), by = .(gene_key, variant_key)]

# --- figure 3C data: credible set overlap ---

pre_cs_input <- prenatal_susie[, .(pip = max(pip, na.rm = TRUE)),
    by = .(gene_key, variant_key, cs_id)]
post_cs_input <- copy(postnatal_ds_susie)
post_cs_input[, cs_id := gene_key]
pre_intervals <- cs_intervals(pre_cs_input, "Prenatal")
post_intervals <- cs_intervals(post_cs_input, "Postnatal")
pre_intervals[, chr_num := suppressWarnings(as.integer(sub("^chr", "", chr)))]
post_intervals[, chr_num := suppressWarnings(as.integer(sub("^chr", "", chr)))]
pre_intervals <- pre_intervals[chr_num %between% c(1L, 22L)]
post_intervals <- post_intervals[chr_num %between% c(1L, 22L)]

pre_intervals <- pre_intervals[, .(
    start = min(start),
    end = max(end),
    center = median(center, na.rm = TRUE),
    max_pip = max(max_pip, na.rm = TRUE),
    n_cs = uniqueN(cs_id)
), by = .(gene_key, chr, chr_num, stage)]
post_intervals <- post_intervals[, .(
    start = min(start),
    end = max(end),
    center = median(center, na.rm = TRUE),
    max_pip = max(max_pip, na.rm = TRUE),
    n_cs = uniqueN(cs_id)
), by = .(gene_key, chr, chr_num, stage)]

cs_plot_dt <- rbindlist(list(
    pre_intervals[, .(gene_key, chr, chr_num, start, end, y = 0.48, evidence = "Prenatal CS")],
    post_intervals[, .(gene_key, chr, chr_num, start, end, y = 0.48, evidence = "Postnatal DS CS")]
), fill = TRUE)

shared_dt <- merge(
    pre_intervals,
    post_intervals,
    by = c("gene_key", "chr", "chr_num"),
    suffixes = c("_pre", "_post"),
    allow.cartesian = TRUE
)
shared_dt <- shared_dt[start_pre <= end_post & start_post <= end_pre]
if (nrow(shared_dt) > 0) {
    shared_dt[, `:=`(
        start = pmax(start_pre, start_post),
        end = pmin(end_pre, end_post),
        y = 0.48,
        evidence = "Shared overlap",
        shared_score = pmin(max_pip_pre, max_pip_post)
    )]
    cs_plot_dt <- rbindlist(list(
        cs_plot_dt,
        shared_dt[, .(gene_key, chr, chr_num, start, end, y, evidence)]
    ), fill = TRUE)
}

cs_plot_dt[, chr := factor(chr, levels = paste0("chr", 1:22))]
chr_counts <- unique(rbindlist(list(pre_intervals[, .(chr, chr_num)], post_intervals[, .(chr, chr_num)])))
shared_genes <- unique(shared_dt$gene_key)

chr_backbone_dt <- cs_plot_dt[, .(chr_end = max(end / 1e6, na.rm = TRUE)), by = .(chr, chr_num)]
chr_backbone_dt[, `:=`(
    chr = factor(chr, levels = paste0("chr", 1:22)),
    xstart = 0,
    xend = pmax(chr_end, 1),
    chr_y = (23 - chr_num) * row_spacing
)]
chr_band_dt <- copy(chr_backbone_dt)
chr_band_dt[, `:=`(
    xmin = -Inf,
    xmax = Inf,
    ymin = chr_y - 0.38,
    ymax = chr_y + 0.38,
    band_fill = fifelse(chr_num %% 2L == 0L, "#F7F8FA", "#FFFFFF")
)]

block_dt <- rbindlist(list(
    pre_intervals[, .(gene_key, chr, chr_num,
                      xmin = start / 1e6, xmax = end / 1e6,
                      ymin = (23 - chr_num) * row_spacing - 0.24,
                      ymax = (23 - chr_num) * row_spacing + 0.24,
                      stage = "Prenatal CS")],
    post_intervals[, .(gene_key, chr, chr_num,
                       xmin = start / 1e6, xmax = end / 1e6,
                       ymin = (23 - chr_num) * row_spacing - 0.24,
                       ymax = (23 - chr_num) * row_spacing + 0.24,
                       stage = "Postnatal DS CS")]
), fill = TRUE)
block_dt[, stage := factor(stage, levels = c("Prenatal CS", "Postnatal DS CS"))]
block_dt[, center_mb := (xmin + xmax) / 2]
block_dt[, width_mb := pmax(xmax - xmin, min_block_width_mb)]
block_dt[, `:=`(
    xmin_plot = center_mb - width_mb / 2,
    xmax_plot = center_mb + width_mb / 2
)]

overlap_label_dt <- shared_dt[is.finite(shared_score) & !grepl("^ENSG", gene_key)][order(chr_num, -shared_score, gene_key)]
overlap_label_dt <- overlap_label_dt[, head(.SD, max_labels_per_chr), by = chr_num]
overlap_label_dt <- overlap_label_dt[order(-shared_score, chr_num, gene_key)][seq_len(min(.N, top_overlap_label_n))]
overlap_label_dt[, `:=`(
    x = (start + end) / 2e6,
    y_target = (23 - chr_num) * row_spacing,
    y = (23 - chr_num) * row_spacing,
    x_label = fifelse((start + end) / 2e6 < 125, left_label_x, right_label_x),
    hjust_label = fifelse((start + end) / 2e6 < 125, 0, 1),
    label = gene_key
)]

log_msg(sprintf("Figure 3C plot rows: intervals=%s, shared=%s, chromosomes=%s",
    comma(nrow(cs_plot_dt)), comma(nrow(shared_dt)), comma(nrow(chr_counts))))

# --- figure 3C plot ---

p3c <- ggplot() +
    geom_rect(data = chr_band_dt,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#F8FAFC",
              inherit.aes = FALSE, color = NA, alpha = 0.75) +
    geom_segment(data = chr_backbone_dt,
                 aes(x = xstart, xend = xend, y = chr_y, yend = chr_y),
                 linewidth = 4.8, lineend = "round", color = "#ECEFF3") +
    geom_rect(data = block_dt[stage == "Postnatal DS CS"],
              aes(xmin = xmin_plot, xmax = xmax_plot, ymin = ymin, ymax = ymax, fill = stage),
              inherit.aes = FALSE, color = NA, alpha = 0.38) +
    geom_rect(data = block_dt[stage == "Prenatal CS"],
              aes(xmin = xmin_plot, xmax = xmax_plot, ymin = ymin, ymax = ymax, fill = stage),
              inherit.aes = FALSE, color = NA, alpha = 0.68) +
    geom_segment(data = overlap_label_dt,
        aes(x = x_label, xend = x, y = y, yend = y_target),
        inherit.aes = FALSE, color = "grey35", linewidth = 0.20,
        linetype = "dashed") +
    geom_point(data = overlap_label_dt,
        aes(x = x, y = y_target),
        inherit.aes = FALSE, color = "grey25", size = 0.55, alpha = 0.75) +
    geom_label(data = overlap_label_dt,
        aes(x = x_label, y = y, label = label, hjust = hjust_label),
        inherit.aes = FALSE, size = 2.65, color = "grey12",
        fontface = "bold", vjust = 0.5, fill = alpha("white", 0.86),
        linewidth = 0.12, label.r = unit(0.06, "lines"),
        label.padding = unit(0.06, "lines")) +
    scale_fill_manual(values = c(
        `Prenatal CS` = "#E69F00",
        `Postnatal DS CS` = "#0072B2"
    )) +
    scale_x_continuous(labels = label_number(accuracy = 1),
                       limits = c(left_label_x - 1, right_label_x + 1),
                       breaks = seq(0, 250, by = 50),
                       expand = c(0, 0)) +
    scale_y_continuous(
        limits = c(0.35, 22.95 * row_spacing),
        breaks = chr_backbone_dt$chr_y,
        labels = as.character(chr_backbone_dt$chr),
        expand = c(0, 0)
    ) +
    labs(x = "Genomic position (Mb)", y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(size = 7, color = "grey35"),
        axis.text.y = element_text(size = 7.2, face = "bold", color = "grey35"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "grey70", linewidth = 0.25),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.width = unit(1.10, "cm"),
        legend.text = element_text(size = 8),
        plot.margin = margin(0, 8, 3, 8)
    )

# --- save figure 3C ---

log_msg("Saving standalone Figure 3C")
ggsave(file.path(fig_dir, "figure3_panelC.pdf"), p3c,
       width = panel_width, height = panel_height, device = cairo_pdf, bg = "white")
ggsave(file.path(fig_dir, "figure3_panelC.png"), p3c,
       width = panel_width, height = panel_height, dpi = dpi, bg = "white")

layout_summary <- data.table(
    figure = "figure3_panelC",
    width_in = panel_width,
    height_in = panel_height,
    prenatal_genes = uniqueN(pre_intervals$gene_key),
    postnatal_ds_genes = uniqueN(post_intervals$gene_key),
    shared_overlap_genes = length(shared_genes),
    chromosome_rows = nrow(chr_counts),
    direct_gene_labels = nrow(overlap_label_dt),
    high_overlap_labels = nrow(overlap_label_dt),
    figure4_example_labels = 0L,
    figure4_examples = "",
    min_block_width_mb = min_block_width_mb,
    row_spacing = row_spacing,
    left_label_x = left_label_x,
    right_label_x = right_label_x,
    notes = "Standalone blue/yellow CS overlay with sparse high-overlap labels placed in fixed side lanes on the exact chromosome row, connected to overlap centers by dashed leaders."
)
fwrite(layout_summary, layout_validation_path, sep = "\t")
log_msg("Finished standalone Figure 3C")
