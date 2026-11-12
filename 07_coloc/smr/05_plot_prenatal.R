#!/usr/bin/env Rscript

# Plot all prenatal SMR candidates that pass FDR and GWAS-threshold filters.
# Tables go to data/07_coloc, figures to figures/07_coloc.
APA_ROOT <- Sys.getenv("APA_ROOT")

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(svglite)
})

# --- config ---
input_file <- file.path(APA_ROOT, "data/07_coloc/smr_prenatal_candidates.tsv")
table_dir <- file.path(APA_ROOT, "data/07_coloc")
fig_dir <- file.path(APA_ROOT, "figures/07_coloc")
output_stem <- "prenatal_smr_candidates"
max_candidates <- NA_integer_
figure_width_mm <- 183
figure_height_mm <- 140

# --- loading prenatal candidates ---

candidates <- fread(input_file)
if (!is.na(max_candidates)) {
    candidates <- head(candidates, max_candidates)
}

candidates[, analysis_id := paste(gsub("_", "-", gwas_id), Gene, sep = " | ")]
candidates[, q_smr := as.numeric(p_SMR_BH)]
candidates[, logq := -log10(q_smr)]
candidates[, q_label := sprintf("q=%.3f", q_smr)]
candidates[, candidate_status := factor(
    candidate_status,
    levels = c("FDR + HEIDI pass", "FDR + HEIDI fail", "FDR + HEIDI unavailable")
)]
candidates[, candidate_status := droplevels(candidate_status)]
candidates[, ci_low := as.numeric(b_SMR) - 1.96 * as.numeric(se_SMR)]
candidates[, ci_high := as.numeric(b_SMR) + 1.96 * as.numeric(se_SMR)]
candidates[, analysis_id := factor(
    analysis_id,
    levels = rev(analysis_id[order(logq)])
)]

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- figure source data ---

source_data <- candidates[, .(
    prenatal_candidate_rank,
    analysis_id = as.character(analysis_id),
    qtl_type,
    gwas_id,
    Gene,
    probeID,
    p_SMR = p_SMR_num,
    p_SMR_BH = q_smr,
    p_SMR_bonferroni,
    n_smr_tests,
    gwas_threshold,
    p_GWAS = p_GWAS_num,
    gwas_threshold_pass,
    b_SMR,
    se_SMR,
    ci_low,
    ci_high,
    p_HEIDI = p_HEIDI_num,
    candidate_status
)]
fwrite(source_data, file.path(table_dir, paste0(output_stem, "_source.tsv")), sep = "\t")

# --- figure theme ---

status_colors <- c(
    "FDR + HEIDI pass" = "#0072B2",
    "FDR + HEIDI fail" = "#C44E52",
    "FDR + HEIDI unavailable" = "#777777"
)
status_shapes <- c(
    "FDR + HEIDI pass" = 16,
    "FDR + HEIDI fail" = 4,
    "FDR + HEIDI unavailable" = 1
)

base_theme <- theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
        axis.line = element_line(linewidth = 0.35, color = "black"),
        axis.ticks = element_line(linewidth = 0.35, color = "black"),
        axis.title = element_text(size = 7.8),
        axis.text = element_text(size = 7.0, color = "black"),
        plot.title = element_text(size = 8.2, face = "bold", hjust = 0),
        plot.subtitle = element_blank(),
        legend.title = element_text(size = 7.2),
        legend.text = element_text(size = 6.8),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        plot.margin = margin(4, 5, 4, 5)
    )

# --- figure a: within-analysis FDR ---

fdr_threshold <- -log10(0.05)

p_fdr <- ggplot(candidates, aes(y = analysis_id)) +
    geom_vline(xintercept = fdr_threshold, color = "grey45", linewidth = 0.4, linetype = "22") +
    geom_segment(
        aes(x = fdr_threshold, xend = logq, yend = analysis_id),
        color = "grey75", linewidth = 0.45
    ) +
    geom_point(
        aes(x = logq, color = candidate_status, shape = candidate_status),
        size = 2.4, stroke = 0.8
    ) +
    geom_text(
        aes(x = logq, label = q_label),
        hjust = -0.35, size = 2.15, color = "grey20"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.03, 0.30))) +
    scale_color_manual(values = status_colors, name = NULL, drop = TRUE) +
    scale_shape_manual(values = status_shapes, name = NULL, drop = TRUE) +
    labs(
        title = "a  Within-analysis FDR",
        x = expression(-log[10](q[SMR])),
        y = NULL
    ) +
    base_theme +
    theme(
        axis.text.y = element_text(size = 7.0),
        panel.grid.major.x = element_blank()
    )

# --- figure b: smr effects ---

p_effect <- ggplot(candidates, aes(y = analysis_id)) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4, linetype = "22") +
    geom_segment(
        aes(x = ci_low, xend = ci_high, yend = analysis_id, color = candidate_status),
        linewidth = 0.65
    ) +
    geom_point(
        aes(x = b_SMR, color = candidate_status, shape = candidate_status),
        size = 2.4, stroke = 0.8
    ) +
    scale_color_manual(values = status_colors, name = NULL, drop = FALSE) +
    scale_shape_manual(values = status_shapes, name = NULL, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) +
    labs(
        title = "b  SMR effect",
        x = expression(beta[SMR]~"(95% CI)"),
        y = NULL
    ) +
    guides(color = "none", shape = "none") +
    base_theme +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major.x = element_blank()
    )

# --- assemble and export figure ---

figure <- p_fdr + p_effect +
    plot_layout(widths = c(1.45, 1)) +
    plot_annotation(
        title = "Prenatal FDR-significant SMR results",
        theme = theme(
            plot.title = element_text(size = 9.5, face = "bold", hjust = 0)
        )
    ) &
    theme(legend.position = "bottom")

figure_width_in <- figure_width_mm / 25.4
figure_height_in <- figure_height_mm / 25.4

ggsave(
    file.path(fig_dir, paste0(output_stem, ".pdf")),
    figure, width = figure_width_in, height = figure_height_in,
    device = cairo_pdf
)
ggsave(
    file.path(fig_dir, paste0(output_stem, ".svg")),
    figure, width = figure_width_in, height = figure_height_in,
    device = svglite
)
ggsave(
    file.path(fig_dir, paste0(output_stem, ".png")),
    figure, width = figure_width_in, height = figure_height_in,
    dpi = 400, bg = "white"
)

# --- layout validation ---

n_labels <- nrow(candidates)
max_label_characters <- max(nchar(as.character(candidates$analysis_id)))
n_legend_items <- uniqueN(candidates$candidate_status)
n_legend_blocks <- 1L
overall_status <- if (
    n_labels <= 15 && max_label_characters <= 32 &&
        n_legend_items <= 3 && n_legend_blocks == 1
) "PASS" else "REVIEW"

layout_validation <- data.table(
    metric = c(
        "visible_labels", "max_label_characters", "legend_items",
        "legend_blocks", "unique_analysis_ids", "overall"
    ),
    value = c(
        n_labels, max_label_characters, n_legend_items,
        n_legend_blocks, uniqueN(candidates$analysis_id), overall_status
    ),
    status = c(
        ifelse(n_labels <= 15, "PASS", "REVIEW"),
        ifelse(max_label_characters <= 32, "PASS", "REVIEW"),
        ifelse(n_legend_items <= 3, "PASS", "REVIEW"),
        ifelse(n_legend_blocks == 1, "PASS", "REVIEW"),
        ifelse(uniqueN(candidates$analysis_id) == n_labels, "PASS", "REVIEW"),
        overall_status
    )
)
fwrite(
    layout_validation,
    file.path(table_dir, paste0(output_stem, "_layout_validation.tsv")),
    sep = "\t"
)

cat(file.path(fig_dir, paste0(output_stem, ".png")), "\n")
