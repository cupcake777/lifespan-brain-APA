# Functional-enrichment figure for the HC3-core G x Stage classes (3 panels + gene matrix).
# Reads data/05_gxs/functional_enrichment_*; writes figures/05_gxs/functional_enrichment.{png,pdf,svg}
# and the layout check table.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- paths and figure constants ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "05_gxs")
out_dir <- file.path(APA_ROOT, "figures", "05_gxs")
enrichment_summary_file <- file.path(D, "functional_enrichment_summary.tsv")
enrichment_supported_file <- file.path(D, "functional_enrichment_supported_terms.tsv")
figure_stem <- "functional_enrichment"
figure_width <- 7.2
figure_height <- 4.8

foreground_order <- c("Opposing", "Amplified", "Attenuated")
library_order <- c("HPO", "KEGG", "Reactome", "WikiPathways")
foreground_colors <- c(Opposing = "#D55E00", Amplified = "#0072B2", Attenuated = "#009E73")
band_colors <- c(Opposing = "#FDF2EC", Amplified = "#EDF6FB", Attenuated = "#EDF8F4")
gene_order <- c("CYCS", "GPX3", "PRDX6", "SOD1",
                "GABRG2", "SCN1A", "TUBB2B", "DHRS4", "DHRS4L2",
                "CEP192", "H3-3A", "MTDH", "ZEB2")
ink <- "#20252A"
grid_color <- "#D9DEE2"
forest_limit <- 250
forest_lower_limit <- 0.5
minimum_text_pt <- 7.0
panel_widths <- c(2.65, 1.55, 2.0, 1.8)
plot_ymax <- 18.0

group_layout <- data.table(
  foreground = foreground_order,
  band_ymin = c(12.9, 4.4, 0.4),
  band_ymax = c(14.7, 12.3, 3.7),
  header_y = c(14.25, 11.85, 3.25),
  term_start = c(13.35, 10.95, 2.35)
)
group_layout[, `:=`(forest_xmin = forest_lower_limit, forest_xmax = forest_limit)]
table_theme <- theme_void(base_size = 7.5) + theme(text = element_text(color = ink))

# --- terms: drop virus-labelled rows, shorten display names ---
summary_dt <- fread(enrichment_summary_file)
supported_all <- fread(enrichment_supported_file)
supported_all[, virus_label := grepl("sars|virus|viral", paste(term, term_label), ignore.case = TRUE)]
omitted_dt <- supported_all[virus_label == TRUE]
display_dt <- copy(supported_all[virus_label == FALSE])
display_dt[, foreground := factor(foreground, levels = foreground_order)]
display_dt[, library := factor(library, levels = library_order)]

display_dt[, term_plain := term_label]
display_dt[, term_plain := sub(" \\(hp:[^)]+\\)", "", term_plain, ignore.case = TRUE)]
display_dt[, term_plain := sub("^Wp ", "", term_plain)]
display_dt[term == "WP_GASTRIC_CANCER_NETWORK_1", term_plain := "Gastric cancer network 1"]
display_dt[term == "WP_CCL18_SIGNALING", term_plain := "CCL18 signaling"]
display_dt[, term_display := term_plain]
display_dt[term == "REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES", term_display := "ROS detoxification"]
display_dt[term == "Generalized tonic-clonic seizures (HP:0002069)", term_display := "Tonic-clonic seizures"]

setorder(display_dt, foreground, library, driver_signature, fdr_within_library, term_plain)
display_dt[, row_in_group := seq_len(.N), by = foreground]
display_dt <- merge(display_dt, group_layout[, .(foreground, term_start)], by = "foreground", all.x = TRUE, sort = FALSE)
display_dt[, foreground := factor(foreground, levels = foreground_order)]
display_dt[, y := term_start - row_in_group + 1]
display_dt[, or_label := sprintf("%.1f", odds_ratio)]
display_dt[, lower_label := sprintf("%.1f", odds_ratio_ci95_low)]
display_dt[, fdr_label := sprintf("%.4f", fdr_within_library)]
setorder(display_dt, -y)

# --- driver-gene membership matrix ---
membership_dt <- display_dt[, .(gene_symbol = unlist(strsplit(driver_genes, "/", fixed = TRUE))),
                            by = .(term, y, foreground)]
membership_dt[, gene_x := match(gene_symbol, gene_order)]
gene_header_dt <- data.table(gene_symbol = gene_order, gene_x = seq_along(gene_order), y = 14.95)

# --- panel 1: term names ---
term_panel <- ggplot() +
  geom_rect(data = group_layout,
            aes(xmin = -Inf, xmax = Inf, ymin = band_ymin, ymax = band_ymax, fill = foreground), inherit.aes = FALSE) +
  geom_text(data = group_layout, aes(x = 0.02, y = header_y, label = foreground, color = foreground),
            hjust = 0, size = 3.0, fontface = "bold") +
  geom_text(data = display_dt, aes(x = 0.03, y = y, label = term_display), hjust = 0, size = 2.65, color = ink) +
  annotate("text", x = 0.03, y = 17.1, label = "Term", hjust = 0, size = 2.85, fontface = "bold") +
  scale_fill_manual(values = band_colors, guide = "none") +
  scale_color_manual(values = foreground_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, plot_ymax), clip = "off") +
  table_theme +
  theme(plot.margin = margin(4, 2, 5, 4))

# --- panel 2: odds-ratio forest ---
forest_panel <- ggplot() +
  geom_rect(data = group_layout,
            aes(xmin = forest_xmin, xmax = forest_xmax, ymin = band_ymin, ymax = band_ymax, fill = foreground),
            inherit.aes = FALSE) +
  geom_vline(xintercept = 1, color = "#8A939A", linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = display_dt,
               aes(x = odds_ratio_ci95_low, xend = odds_ratio, y = y, yend = y, color = foreground),
               linewidth = 0.65, lineend = "round") +
  geom_segment(data = display_dt,
               aes(x = odds_ratio_ci95_low, xend = odds_ratio_ci95_low, y = y - 0.14, yend = y + 0.14, color = foreground),
               linewidth = 0.65) +
  geom_point(data = display_dt, aes(x = odds_ratio, y = y, color = foreground), size = 2.6) +
  annotate("text", x = 14, y = 17.1, label = "Odds ratio", size = 2.85, fontface = "bold") +
  scale_fill_manual(values = band_colors, guide = "none") +
  scale_color_manual(values = foreground_colors, guide = "none") +
  scale_x_log10(limits = c(forest_lower_limit, forest_limit), breaks = c(1, 10, 100),
                expand = expansion(mult = c(0.02, 0.03))) +
  coord_cartesian(ylim = c(0, plot_ymax), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 7.5) +
  theme(text = element_text(color = ink),
        axis.text.x = element_text(size = minimum_text_pt, color = ink),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.grid.major.x = element_line(color = grid_color, linewidth = 0.35),
        plot.margin = margin(4, 2, 5, 2))

# --- panel 3: OR, lower CI bound, FDR ---
statistics_panel <- ggplot() +
  geom_rect(data = group_layout,
            aes(xmin = -Inf, xmax = Inf, ymin = band_ymin, ymax = band_ymax, fill = foreground), inherit.aes = FALSE) +
  geom_text(data = display_dt, aes(x = 0.18, y = y, label = or_label), hjust = 1, size = 2.55, color = ink) +
  geom_text(data = display_dt, aes(x = 0.58, y = y, label = lower_label), hjust = 1, size = 2.55, color = ink) +
  geom_text(data = display_dt, aes(x = 0.96, y = y, label = fdr_label), hjust = 1, size = 2.55, color = ink) +
  annotate("text", x = 0.18, y = 17.1, label = "OR", hjust = 1, size = 2.85, fontface = "bold") +
  annotate("text", x = 0.58, y = 17.1, label = "95% LB", hjust = 1, size = 2.85, fontface = "bold") +
  annotate("text", x = 0.96, y = 17.1, label = "FDR", hjust = 1, size = 2.85, fontface = "bold") +
  scale_fill_manual(values = band_colors, guide = "none") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, plot_ymax), clip = "off") +
  table_theme +
  theme(plot.margin = margin(4, 2, 5, 2))

# --- panel 4: driver-gene matrix ---
gene_panel <- ggplot() +
  geom_rect(data = group_layout,
            aes(xmin = -Inf, xmax = Inf, ymin = band_ymin, ymax = band_ymax, fill = foreground), inherit.aes = FALSE) +
  geom_point(data = membership_dt, aes(x = gene_x, y = y, color = foreground), shape = 15, size = 2.7) +
  geom_text(data = gene_header_dt, aes(x = gene_x, y = y, label = gene_symbol),
            angle = 90, hjust = 0, vjust = 0.5, size = 2.5, color = ink) +
  annotate("text", x = 0.5, y = 17.55, label = "Driver genes", hjust = 0, size = 2.85, fontface = "bold") +
  scale_fill_manual(values = band_colors, guide = "none") +
  scale_color_manual(values = foreground_colors, guide = "none") +
  coord_cartesian(xlim = c(0.5, length(gene_order) + 0.5), ylim = c(0, plot_ymax), clip = "off") +
  table_theme +
  theme(plot.margin = margin(4, 5, 5, 2))

# --- assemble and save ---
functional_enrichment_figure <- term_panel + forest_panel + statistics_panel + gene_panel +
  plot_layout(widths = panel_widths)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, paste0(figure_stem, ".png")), functional_enrichment_figure,
       width = figure_width, height = figure_height, dpi = 600, bg = "white")
ggsave(file.path(out_dir, paste0(figure_stem, ".pdf")), functional_enrichment_figure,
       width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, paste0(figure_stem, ".svg")), functional_enrichment_figure,
       width = figure_width, height = figure_height, device = svglite::svglite, bg = "white")

# --- layout check: measured text widths (inches) against panel space ---
observed_genes <- sort(unique(membership_dt$gene_symbol))
panel_width_inches <- figure_width * panel_widths / sum(panel_widths)
term_label_space_inches <- panel_width_inches[1] * 0.92
gene_column_pitch_pt <- panel_width_inches[4] * 72 / length(gene_order)

grDevices::pdf(NULL)
grid::grid.newpage()
grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = 7.5)))
text_width_in <- function(x) max(grid::convertWidth(grid::stringWidth(x), "inches", valueOnly = TRUE))
term_label_width_inches <- text_width_in(display_dt$term_display)
or_label_width_inches <- text_width_in(display_dt$or_label)
lower_label_width_inches <- text_width_in(display_dt$lower_label)
fdr_label_width_inches <- text_width_in(display_dt$fdr_label)
lower_header_width_inches <- text_width_in("95% LB")
fdr_header_width_inches <- text_width_in("FDR")
gene_label_width_inches <- text_width_in(gene_order)
grid::popViewport()
grDevices::dev.off()

term_labels_fit <- term_label_width_inches < term_label_space_inches
statistics_gap_1_inches <- panel_width_inches[3] * (0.58 - 0.18) - max(lower_label_width_inches, lower_header_width_inches)
statistics_gap_2_inches <- panel_width_inches[3] * (0.96 - 0.58) - max(fdr_label_width_inches, fdr_header_width_inches)
statistics_gap_inches <- min(statistics_gap_1_inches, statistics_gap_2_inches)
statistics_columns_fit <- statistics_gap_inches > 0.10
gene_columns_fit <- gene_column_pitch_pt > minimum_text_pt
gene_header_band_inches <- figure_height * (17.55 - 14.7) / plot_ymax
gene_header_fit <- gene_label_width_inches + 0.20 < gene_header_band_inches
layout_checks_pass <- term_labels_fit && statistics_columns_fit && gene_columns_fit && gene_header_fit

layout_validation <- rbindlist(list(
  data.table(check = "source", metric = "canonical_supported_terms", value = as.character(nrow(supported_all))),
  data.table(check = "source", metric = "all_hc3_supported_terms", value = as.character(sum(summary_dt[foreground == "All HC3 core", supported_terms]))),
  data.table(check = "display", metric = "virus_labelled_terms_omitted", value = as.character(nrow(omitted_dt))),
  data.table(check = "display", metric = "nonviral_terms_rendered", value = as.character(nrow(display_dt))),
  data.table(check = "display", metric = "foreground_background_bands", value = as.character(nrow(group_layout))),
  data.table(check = "statistics", metric = "finite_or_lower_ci_fdr", value = as.character(all(
    is.finite(display_dt$odds_ratio) &
      is.finite(display_dt$odds_ratio_ci95_low) &
      is.finite(display_dt$fdr_within_library)
  ))),
  data.table(check = "statistics", metric = "unbounded_upper_ci_rows", value = as.character(sum(!is.finite(display_dt$odds_ratio_ci95_high)))),
  data.table(check = "statistics", metric = "truncated_ci_arrows", value = "0"),
  data.table(check = "statistics", metric = "lower_bound_segments", value = as.character(nrow(display_dt))),
  data.table(check = "statistics", metric = "lower_bound_caps", value = as.character(nrow(display_dt))),
  data.table(check = "genes", metric = "direct_gene_labels", value = as.character(length(gene_order))),
  data.table(check = "genes", metric = "gene_labels_match_membership", value = as.character(setequal(gene_order, observed_genes))),
  data.table(check = "genes", metric = "membership_cells", value = as.character(nrow(membership_dt))),
  data.table(check = "layout", metric = "maximum_term_label_characters", value = as.character(max(nchar(display_dt$term_display)))),
  data.table(check = "layout", metric = "maximum_term_label_width_inches", value = sprintf("%.3f", term_label_width_inches)),
  data.table(check = "layout", metric = "term_label_space_inches", value = sprintf("%.3f", term_label_space_inches)),
  data.table(check = "layout", metric = "term_labels_fit", value = as.character(term_labels_fit)),
  data.table(check = "layout", metric = "statistics_gap_inches", value = sprintf("%.3f", statistics_gap_inches)),
  data.table(check = "layout", metric = "statistics_columns_fit", value = as.character(statistics_columns_fit)),
  data.table(check = "layout", metric = "gene_label_rotation_degrees", value = "90"),
  data.table(check = "layout", metric = "gene_column_pitch_pt", value = sprintf("%.2f", gene_column_pitch_pt)),
  data.table(check = "layout", metric = "gene_columns_fit", value = as.character(gene_columns_fit)),
  data.table(check = "layout", metric = "gene_header_band_inches", value = sprintf("%.3f", gene_header_band_inches)),
  data.table(check = "layout", metric = "maximum_gene_label_width_inches", value = sprintf("%.3f", gene_label_width_inches)),
  data.table(check = "layout", metric = "gene_header_fit", value = as.character(gene_header_fit)),
  data.table(check = "layout", metric = "minimum_visible_text_pt", value = sprintf("%.1f", minimum_text_pt)),
  data.table(check = "layout", metric = "final_width_scale", value = "1.0"),
  data.table(check = "layout", metric = "layout_checks_pass", value = as.character(layout_checks_pass)),
  data.table(check = "layout", metric = "panel_count", value = "1"),
  data.table(check = "export", metric = "dimensions_inches", value = sprintf("%.1f x %.1f", figure_width, figure_height)),
  data.table(check = "export", metric = "files", value = paste0(figure_stem, ".png|", figure_stem, ".pdf|", figure_stem, ".svg"))
))
fwrite(layout_validation, file.path(out_dir, paste0(figure_stem, "_layout_validation.tsv")), sep = "\t")

stopifnot(
  nrow(supported_all) == 14L,
  nrow(omitted_dt) == 4L,
  nrow(display_dt) == 10L,
  sum(summary_dt[foreground == "All HC3 core", supported_terms]) == 0L,
  all(!is.finite(display_dt$odds_ratio_ci95_high)),
  all(is.finite(display_dt$odds_ratio)),
  all(is.finite(display_dt$odds_ratio_ci95_low)),
  all(is.finite(display_dt$fdr_within_library)),
  all(display_dt$odds_ratio_ci95_low <= display_dt$odds_ratio),
  setequal(gene_order, observed_genes),
  nrow(membership_dt) == sum(display_dt$driver_count),
  layout_checks_pass,
  !any(grepl("sars|virus|viral", paste(display_dt$term, display_dt$term_label), ignore.case = TRUE)),
  all(file.exists(file.path(out_dir, paste0(figure_stem, c(".png", ".pdf", ".svg")))))
)
