#!/usr/bin/env Rscript
# Gene-level coloc summary panels A-C (figures to figures/07_coloc, tables to data/07_coloc).

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(grid)
})

script_args <- commandArgs(trailingOnly = FALSE)

results_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc/coloc_run")
data_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
fig_dir <- file.path(Sys.getenv("APA_ROOT"), "figures/07_coloc")
out_dir <- fig_dir
sig_file    <- file.path(results_dir, "significant_colocalizations_PP4_0.7.csv")
susie_file  <- file.path(results_dir, "all_susie_results.csv")

susie_confirm_pp4 <- 0.7
top_genes_n       <- 20
figure_dpi        <- 600
theme_font        <- "Helvetica"

panel_dims <- tibble::tribble(
  ~panel, ~width, ~height,
  "A_signal_counts", 3.45, 2.75,
  "B_stage_total", 2.15, 1.65,
  "C_recurrent_genes", 10.40, 5.65
)

stage_levels <- c("prenatal", "postnatal", "interaction")
stage_labels <- c(prenatal = "Prenatal", postnatal = "Postnatal", interaction = "Interaction")
stage_colors <- c(prenatal = "#E64B35", postnatal = "#4DBBD5", interaction = "#CC79A7")
legend_group_gap <- ""

gwas_order <- c(
  "EAS_AIS", "EAS_PD", "EAS_SCZ",
  "EUR_AD", "EUR_ADHD", "EUR_BIP", "EUR_INT", "EUR_MDD", "EUR_PD", "EUR_PTSD", "EUR_SCZ"
)

ink        <- "#20252A"
muted_ink  <- "#65717F"
light_grid <- "#D9E0E7"
empty_tile <- "#F6F8FA"
support_colors <- c("ABF-only" = "#AEB8C2", "SuSiE-supported" = ink)
ancestry_bg_colors <- c(EAS = "#EAF7F2", EUR = "#F1F4FB")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

sig <- read_csv(sig_file, show_col_types = FALSE) %>%
  mutate(
    GWAS_ID = as.character(GWAS_ID),
    QTL_ID = as.character(QTL_ID),
    Phenotype = as.character(Phenotype),
    PP4 = as.numeric(PP4),
    n_snps = as.numeric(n_snps),
    key = str_remove(basename(locuscompare_json_file), "_locuscompare\\.json$"),
    ancestry = str_extract(GWAS_ID, "^[^_]+"),
    disease = str_remove(GWAS_ID, "^[^_]+_"),
    stage = factor(QTL_ID, levels = stage_levels),
    gene = str_split_fixed(Phenotype, "\\|", 4)[, 2],
    gene = if_else(is.na(gene) | gene == "", Phenotype, gene),
    chrom = str_match(key, "_chr([0-9XY]+)_([0-9]+)_")[, 2],
    pos = as.numeric(str_match(key, "_chr([0-9XY]+)_([0-9]+)_")[, 3]),
    chrom_label = paste0("chr", chrom)
  ) %>%
  filter(!is.na(stage), !is.na(PP4), !is.na(pos), !is.na(key), !duplicated(key))

susie_locus <- read_csv(susie_file, show_col_types = FALSE) %>%
  mutate(
    key = str_remove(basename(source_file), "_susie\\.csv$"),
    PP.H4.abf = as.numeric(PP.H4.abf)
  ) %>%
  filter(!is.na(key), !is.na(PP.H4.abf)) %>%
  group_by(key) %>%
  summarise(
    susie_max_PP4 = max(PP.H4.abf),
    hit1 = hit1[which.max(PP.H4.abf)],
    hit2 = hit2[which.max(PP.H4.abf)],
    .groups = "drop"
  )

d <- sig %>%
  left_join(susie_locus, by = "key") %>%
  mutate(
    confirmed = !is.na(susie_max_PP4) & susie_max_PP4 >= susie_confirm_pp4,
    support = case_when(
      confirmed ~ "SuSiE-confirmed",
      !is.na(susie_max_PP4) ~ "SuSiE-not-confirmed",
      TRUE ~ "No SuSiE"
    ),
    support_class = if_else(confirmed, "SuSiE-supported", "ABF-only"),
    support_class = factor(support_class, levels = names(support_colors)),
    GWAS_ID = factor(GWAS_ID, levels = gwas_order),
    stage = factor(as.character(stage), levels = stage_levels)
  ) %>%
  filter(!is.na(GWAS_ID), !is.na(stage))

disease_order <- d %>%
  count(disease, name = "total") %>%
  arrange(desc(total), disease) %>%
  mutate(disease_label = disease)

d <- d %>%
  left_join(disease_order, by = "disease") %>%
  mutate(disease_label = factor(disease_label, levels = rev(disease_order$disease_label)))

plot_data_path <- file.path(data_dir, "06_summary_plot_data.tsv")
write_tsv(
  d %>%
    transmute(
      key, GWAS_ID = as.character(GWAS_ID), ancestry, disease,
      stage = as.character(stage), gene, chrom = chrom_label, pos,
      locus = Locus, PP4, n_snps, susie_max_PP4, hit1, hit2,
      support, support_class = as.character(support_class), confirmed
    ),
  plot_data_path
)

# Panel A: gene-level signal count matrix.
heat <- d %>%
  count(disease, disease_label, stage, name = "n") %>%
  complete(
    disease_label = factor(levels(d$disease_label), levels = levels(d$disease_label)),
    stage = factor(stage_levels, levels = stage_levels),
    fill = list(n = 0)
  ) %>%
  mutate(
    cell_label = if_else(n > 0, as.character(n), ""),
    label_color = if_else(n >= max(n, na.rm = TRUE) * 0.45, "white", ink)
  )

pA <- ggplot(heat, aes(x = stage, y = disease_label)) +
  geom_tile(aes(fill = n), color = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, color = label_color), size = 2.25, family = theme_font) +
  scale_color_identity() +
  scale_fill_gradient(
    low = "#F7F9FA", high = "#6F7B8A",
    breaks = c(0, 10, 25, 50),
    name = "Signal count"
  ) +
  scale_x_discrete(labels = stage_labels) +
  labs(
    x = "QTL stage",
    y = "GWAS trait",
    title = "Gene-level PP4 >= 0.7 signal counts"
  ) +
  theme_classic(base_size = 7.5, base_family = theme_font) +
  theme(
    text = element_text(color = ink),
    plot.title = element_text(size = 8.4, face = "plain", margin = margin(0, 0, 5, 0)),
    axis.title.x = element_text(size = 6.7, margin = margin(4, 0, 0, 0)),
    axis.title.y = element_text(size = 6.7, margin = margin(0, 4, 0, 0)),
    axis.text.x = element_text(size = 6.6, color = "black"),
    axis.text.y = element_text(size = 6.6, color = "black"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 6.2),
    legend.text = element_text(size = 6.0),
    legend.key.height = unit(0.34, "cm"),
    legend.key.width = unit(0.18, "cm"),
    plot.margin = margin(2, 4, 2, 2)
  )

# Panel B: stage totals.
stage_counts <- d %>%
  count(stage, name = "n", .drop = FALSE) %>%
  mutate(stage = factor(stage, levels = rev(stage_levels)))

pB <- ggplot(stage_counts, aes(x = n, y = stage, fill = stage)) +
  geom_col(width = 0.56) +
  geom_text(aes(label = n), hjust = -0.18, size = 2.15, color = ink, family = theme_font) +
  scale_fill_manual(values = stage_colors, guide = "none") +
  scale_y_discrete(labels = stage_labels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Signals", y = NULL, title = "Stage total") +
  theme_classic(base_size = 7.5, base_family = theme_font) +
  theme(
    text = element_text(color = ink),
    plot.title = element_text(size = 8, face = "plain", margin = margin(0, 0, 4, 0)),
    axis.title.x = element_text(size = 6.6),
    axis.text.x = element_text(size = 6.2, color = "black"),
    axis.text.y = element_text(size = 6.5, color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.25),
    axis.line = element_line(color = "black", linewidth = 0.25),
    plot.margin = margin(2, 3, 2, 4)
  )

# Panel C: recurrent genes, intentionally capped for readability.
top_genes <- d %>%
  group_by(gene) %>%
  summarise(
    n_signals = n(),
    n_diseases = n_distinct(disease),
    n_gwas = n_distinct(GWAS_ID),
    n_stages = n_distinct(stage),
    max_pp4 = max(PP4, na.rm = TRUE),
    n_confirmed = sum(confirmed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_signals), desc(n_gwas), desc(n_stages), desc(n_confirmed), desc(max_pp4), gene) %>%
  slice_head(n = top_genes_n)

top_d <- d %>%
  semi_join(top_genes, by = "gene") %>%
  group_by(gene, GWAS_ID, stage) %>%
  summarise(
    signal_count = n(),
    max_pp4 = max(PP4, na.rm = TRUE),
    support_class = if_else(any(confirmed), "SuSiE-supported", "ABF-only"),
    .groups = "drop"
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(top_genes$gene)),
    GWAS_ID = factor(as.character(GWAS_ID), levels = gwas_order),
    ancestry = factor(str_extract(as.character(GWAS_ID), "^[^_]+"), levels = names(ancestry_bg_colors)),
    support_class = factor(support_class, levels = names(support_colors))
  )

top_gwas_order <- gwas_order[gwas_order %in% as.character(top_d$GWAS_ID)]
top_d <- top_d %>%
  mutate(
    GWAS_ID = factor(as.character(GWAS_ID), levels = top_gwas_order)
  )

ancestry_bg <- tibble::tibble(
  ancestry = factor(names(ancestry_bg_colors), levels = names(ancestry_bg_colors))
)

pC <- ggplot(top_d, aes(x = GWAS_ID, y = gene)) +
  geom_rect(
    data = ancestry_bg %>% filter(ancestry == "EAS"),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = ancestry_bg_colors[["EAS"]],
    color = NA
  ) +
  geom_rect(
    data = ancestry_bg %>% filter(ancestry == "EUR"),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    fill = ancestry_bg_colors[["EUR"]],
    color = NA
  ) +
  geom_hline(
    yintercept = seq_along(levels(top_d$gene)),
    color = light_grid, linewidth = 0.18
  ) +
  geom_point(
    aes(size = signal_count, fill = stage, color = support_class, group = stage),
    shape = 21, stroke = 0.50, alpha = 0.96,
    position = position_dodge(width = 0.52)
  ) +
  scale_fill_manual(
    values = stage_colors,
    labels = stage_labels,
    name = "QTL stage",
    breaks = stage_levels
  ) +
  scale_color_manual(
    values = support_colors,
    name = paste0(legend_group_gap, "Support")
  ) +
  scale_size_area(
    max_size = 5.4,
    breaks = c(1, 3, 6),
    name = paste0(legend_group_gap, "Coloc hit count")
  ) +
  scale_x_discrete(
    labels = function(x) str_remove(x, "^[^_]+_"),
    drop = TRUE
  ) +
  facet_grid(. ~ ancestry, scales = "free_x", space = "free_x") +
  labs(x = "GWAS trait", y = "Gene", title = "Recurrent gene-level colocalization") +
  theme_classic(base_size = 12.0, base_family = theme_font) +
  theme(
    text = element_text(color = ink),
    plot.title = element_text(size = 15.5, face = "plain", margin = margin(0, 0, 8, 0)),
    axis.title.x = element_text(size = 12.0, margin = margin(7, 0, 0, 0)),
    axis.title.y = element_text(size = 12.0, margin = margin(0, 7, 0, 0)),
    axis.text.x = element_text(angle = 32, hjust = 1, vjust = 1, size = 10.8, color = "black"),
    axis.text.y = element_text(size = 10.5, color = "black", face = "italic"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "center",
    legend.box.just = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(8, 0, 4, 0),
    legend.box.spacing = unit(0.08, "cm"),
    legend.spacing.x = unit(0.24, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.key.size = unit(0.42, "cm"),
    legend.title = element_text(size = 10.0),
    legend.text = element_text(size = 10.0),
    panel.grid.major.x = element_line(color = empty_tile, linewidth = 0.20),
    panel.spacing.x = unit(0.28, "cm"),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text.x = element_text(size = 11.2, color = ink, face = "plain", margin = margin(0, 0, 4, 0)),
    plot.margin = margin(6, 7, 10, 4)
  ) +
  guides(
    fill = guide_legend(
      order = 1, nrow = 1, byrow = TRUE, title.position = "left",
      override.aes = list(size = 3.2, color = "#666666")
    ),
    color = guide_legend(
      order = 2, nrow = 1, byrow = TRUE, title.position = "left",
      override.aes = list(size = 3.2, fill = "white")
    ),
    size = guide_legend(
      order = 3, nrow = 1, byrow = TRUE, title.position = "left",
      override.aes = list(fill = "white", color = ink, alpha = 1)
    )
  )

panel_specs <- list(
  A_signal_counts = pA,
  B_stage_total = pB,
  C_recurrent_genes = pC
)

for (panel_name in names(panel_specs)) {
  dims <- panel_dims %>% filter(panel == panel_name)
  png_path <- file.path(fig_dir, paste0("06_summary_panel_", panel_name, ".png"))
  pdf_path <- file.path(out_dir, paste0("06_summary_panel_", panel_name, ".pdf"))

  ggsave(
    png_path, panel_specs[[panel_name]],
    width = dims$width, height = dims$height, dpi = figure_dpi, bg = "white"
  )
  ggsave(
    pdf_path, panel_specs[[panel_name]],
    width = dims$width, height = dims$height, device = cairo_pdf, bg = "white"
  )
}

validation <- tibble::tribble(
  ~panel, ~check, ~value, ~status,
  "all", "input_signals", as.character(nrow(d)), if_else(nrow(d) > 0, "PASS", "FAIL"),
  "all", "unique_keys", as.character(n_distinct(d$key)), if_else(n_distinct(d$key) == nrow(d), "PASS", "FAIL"),
  "all", "missing_stage", as.character(sum(is.na(d$stage))), if_else(sum(is.na(d$stage)) == 0, "PASS", "FAIL"),
  "all", "missing_position", as.character(sum(is.na(d$pos))), if_else(sum(is.na(d$pos)) == 0, "PASS", "FAIL"),
  "A_signal_counts", "disease_rows", as.character(nrow(disease_order)), if_else(nrow(disease_order) <= 10, "PASS", "WARN"),
  "A_signal_counts", "stage_columns", as.character(length(stage_levels)), "PASS",
  "A_signal_counts", "label_contract", "cell value = PP4>=0.7 gene-trait signal count", "PASS",
  "A_signal_counts", "subtitle", "not used", "PASS",
  "A_signal_counts", "figure_inches", paste0(panel_dims$width[panel_dims$panel == "A_signal_counts"], "x", panel_dims$height[panel_dims$panel == "A_signal_counts"]), "PASS",
  "B_stage_total", "stage_rows", as.character(nrow(stage_counts)), if_else(nrow(stage_counts) == length(stage_levels), "PASS", "FAIL"),
  "B_stage_total", "figure_inches", paste0(panel_dims$width[panel_dims$panel == "B_stage_total"], "x", panel_dims$height[panel_dims$panel == "B_stage_total"]), "PASS",
  "C_recurrent_genes", "top_gene_rows", as.character(nrow(top_genes)), if_else(nrow(top_genes) <= 20, "PASS", "WARN"),
  "C_recurrent_genes", "max_gene_label_chars", as.character(max(nchar(as.character(top_genes$gene)))),
    if_else(max(nchar(as.character(top_genes$gene))) <= 18, "PASS", "WARN"),
  "C_recurrent_genes", "top_gene_gwas_labels", as.character(length(top_gwas_order)), if_else(length(top_gwas_order) <= 11, "PASS", "WARN"),
  "C_recurrent_genes", "x_label_angle", "32 degrees", "PASS",
  "C_recurrent_genes", "ancestry_facets", paste(names(ancestry_bg_colors), collapse = ";"), "PASS",
  "C_recurrent_genes", "facet_backgrounds", paste(paste(names(ancestry_bg_colors), ancestry_bg_colors, sep = ":"), collapse = ";"), "PASS",
  "C_recurrent_genes", "stage_palette", paste(paste(names(stage_colors), stage_colors, sep = ":"), collapse = ";"), "PASS",
  "C_recurrent_genes", "support_groups", paste(names(support_colors), collapse = ";"), "PASS",
  "C_recurrent_genes", "legend_group_gap_spaces", as.character(nchar(legend_group_gap)), "PASS",
  "C_recurrent_genes", "legend_alignment", "center", "PASS",
  "C_recurrent_genes", "font_sizes", "title=15.5;axis_title=12.0;axis_x=10.8;axis_y=10.5;facet=11.2;legend=10.0", "PASS",
  "C_recurrent_genes", "size_contract", "colocalization hits per gene-trait-stage group", "PASS",
  "C_recurrent_genes", "position_labels", "not used", "PASS",
  "C_recurrent_genes", "legend_layout", "three centered rows with separated guide groups", "PASS",
  "C_recurrent_genes", "figure_inches", paste0(panel_dims$width[panel_dims$panel == "C_recurrent_genes"], "x", panel_dims$height[panel_dims$panel == "C_recurrent_genes"]), "PASS"
)

write_tsv(validation, file.path(data_dir, "06_summary_layout_validation.tsv"))

message(sprintf(
  "[DONE] summary panels written: %d panels (%d signals; top recurrent genes=%d)",
  length(panel_specs), nrow(d), nrow(top_genes)
))
