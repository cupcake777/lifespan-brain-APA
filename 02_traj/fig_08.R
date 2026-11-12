# Directional-membership functional enrichment (GO/KEGG, GWAS/HPO) as one panel.
suppressPackageStartupMessages({
  library(patchwork)
})
source("fig_common.R")

# --- helpers ---

read_gmt <- function(path) {
  lines <- readLines(path)
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    list(name = fields[1], genes = fields[-(1:2)])
  })
  setNames(lapply(sets, `[[`, "genes"),
           vapply(sets, `[[`, character(1), "name"))
}

run_enrichment <- function(fg_sets, background, gene_sets, library_name) {
  bg <- intersect(unique(background), unique(unlist(gene_sets)))
  eligible_sets <- gene_sets[
    vapply(gene_sets, function(x) {
      n <- length(intersect(x, bg))
      n >= MIN_SET_SIZE && n <= MAX_SET_SIZE
    }, logical(1))
  ]

  rows <- lapply(names(fg_sets), function(group_name) {
    fg <- intersect(unique(fg_sets[[group_name]]), bg)
    ref <- setdiff(bg, fg)

    rbindlist(lapply(names(eligible_sets), function(term_name) {
      term_genes <- intersect(eligible_sets[[term_name]], bg)
      a <- length(intersect(fg, term_genes))
      b <- length(fg) - a
      c <- length(intersect(ref, term_genes))
      d <- length(ref) - c
      test <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                          alternative = "greater")

      data.table(
        library_name = library_name,
        gene_set = group_name,
        term = term_name,
        pvalue = test$p.value,
        fisher_odds_ratio = unname(test$estimate),
        n_fg = a,
        n_bg = c,
        n_fg_total = length(fg),
        n_bg_total = length(bg)
      )
    }))
  })

  out <- rbindlist(rows)
  out[, library_fdr := p.adjust(pvalue, "BH"), by = gene_set]
  out
}

# --- config ---

GENE_AGG_FILE <- file.path(OUT_DIR, "05_gene_aggregation.tsv")
ADJUSTED_FILE <- file.path(OUT_DIR, "08_enrich_adjusted.tsv")
GO_BP_FILE <- file.path(APA_ROOT, "ref", "c5.go.bp.v2024.1.Hs.symbols.gmt")
GO_MF_FILE <- file.path(APA_ROOT, "ref", "c5.go.mf.v2024.1.Hs.symbols.gmt")
KEGG_FILE <- file.path(APA_ROOT, "ref", "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")
GWAS_FILE <- file.path(APA_ROOT, "ref", "GWAS_Catalog_2025.txt")
HPO_FILE <- file.path(APA_ROOT, "ref", "Human_Phenotype_Ontology.gmt")

MIN_SET_SIZE <- 10L
MAX_SET_SIZE <- 500L
TOP_PER_GROUP <- 6L
CANDIDATE_LIBRARY_FDR <- 0.10
REPORT_ADJUSTED_FDR <- 0.05

OUT_RESULTS <- file.path(OUT_DIR, "08_enrich_two_levels.tsv")
OUT_PLOT_DATA <- file.path(OUT_DIR, "08_enrich_plot_data.tsv")
OUT_LAYOUT <- file.path(OUT_DIR, "08_enrich_layout_validation.tsv")
OUT_PNG <- file.path(PANEL_DIR, "08_enrichment.png")
OUT_PDF <- file.path(PANEL_DIR, "08_enrichment.pdf")

level_order <- c("Directional membership", "Five-state")
source_order <- c("GO/KEGG", "GWAS/HPO")
state_order <- c("C1-spec", "C1-main", "Balanced", "C2-main", "C2-spec")
directional_order <- c("C1-spec", "Shared", "C2-spec")

group_colors <- c(
  "C1-spec" = "#56B4E9",
  "C1-main" = "#0072B2",
  "Balanced" = "#7A7A7A",
  "C2-main" = "#D55E00",
  "C2-spec" = "#E69F00",
  "Shared" = "#A65628"
)

group_shapes <- c(
  "C1-spec" = 21,
  "C1-main" = 22,
  "Balanced" = 23,
  "C2-main" = 22,
  "C2-spec" = 21,
  "Shared" = 24
)

group_linetypes <- c(
  "C1-spec" = "solid",
  "C1-main" = "longdash",
  "Balanced" = "dotted",
  "C2-main" = "longdash",
  "C2-spec" = "solid",
  "Shared" = "dotdash"
)

source_band_colors <- c(
  "GO/KEGG" = "#EAF7F2",
  "GWAS/HPO" = "#FBEAF3"
)

# --- loading data and defining gene sets ---

gene_agg <- fread(GENE_AGG_FILE)
gene_agg[, gene_state := fcase(
  n_c1 > 0 & n_c2 == 0, "C1-spec",
  n_c1 > n_c2 & n_c2 > 0, "C1-main",
  n_c1 == n_c2 & n_c1 > 0, "Balanced",
  n_c2 > n_c1 & n_c1 > 0, "C2-main",
  n_c2 > 0 & n_c1 == 0, "C2-spec"
)]
gene_agg[, directional_membership := fcase(
  gene_state == "C1-spec", "C1-spec",
  gene_state == "C2-spec", "C2-spec",
  default = "Shared"
)]

background_genes <- unique(gene_agg$gene_symbol)
gene_sets_l1 <- split(gene_agg$gene_symbol, gene_agg$directional_membership)
gene_sets_l2 <- split(gene_agg$gene_symbol, gene_agg$gene_state)
gene_sets_l1 <- gene_sets_l1[directional_order]
gene_sets_l2 <- gene_sets_l2[state_order]

# --- loading annotation sets ---

gwas_lines <- readLines(GWAS_FILE)
gwas_sets <- setNames(
  lapply(gwas_lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][-(1:2)]),
  vapply(gwas_lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][1],
         character(1))
)
annotation_libraries <- list(
  "GO BP" = read_gmt(GO_BP_FILE),
  "GO MF" = read_gmt(GO_MF_FILE),
  "KEGG" = read_gmt(KEGG_FILE),
  "GWAS Catalog" = gwas_sets,
  "HPO" = read_gmt(HPO_FILE)
)
source_groups <- c(
  "GO BP" = "GO/KEGG", "GO MF" = "GO/KEGG", "KEGG" = "GO/KEGG",
  "GWAS Catalog" = "GWAS/HPO", "HPO" = "GWAS/HPO"
)

# --- enrichment at both hierarchy levels ---

enrichment_results <- rbindlist(lapply(names(annotation_libraries), function(library_name) {
  gene_sets <- annotation_libraries[[library_name]]
  source_group <- unname(source_groups[library_name])
  l1 <- run_enrichment(gene_sets_l1, background_genes, gene_sets, library_name)
  l1[, `:=`(source_group = source_group, level = "Directional membership")]
  l2 <- run_enrichment(gene_sets_l2, background_genes, gene_sets, library_name)
  l2[, `:=`(source_group = source_group, level = "Five-state")]
  rbindlist(list(l1, l2))
}))
enrichment_results[, `:=`(
  source_group = factor(source_group, levels = source_order),
  level = factor(level, levels = level_order),
  gene_set = factor(gene_set,
                    levels = unique(c(directional_order, state_order)))
)]

adjusted_original <- fread(ADJUSTED_FILE)
adjusted_directional <- rbindlist(list(
  adjusted_original[level == "Five-state" & gene_set == "C1-spec"],
  adjusted_original[level == "Membership" & gene_set == "Shared"],
  adjusted_original[level == "Five-state" & gene_set == "C2-spec"]
))
adjusted_directional[, level := "Directional membership"]
adjusted_directional[, adjusted_fdr := p.adjust(adjusted_pvalue, "BH"),
                     by = .(library_name, gene_set)]
adjusted <- rbindlist(list(
  adjusted_directional,
  adjusted_original[level == "Five-state"]
), use.names = TRUE)
result_keys <- c("library_name", "source_group", "level", "gene_set", "term")
enrichment_results <- merge(
  enrichment_results,
  adjusted[, c(result_keys, "adjusted_odds_ratio", "adjusted_pvalue",
               "adjusted_fdr"), with = FALSE],
  by = result_keys,
  all.x = TRUE
)
enrichment_results[, report_pass := adjusted_fdr < REPORT_ADJUSTED_FDR]
fwrite(enrichment_results, OUT_RESULTS, sep = "\t")

# --- significant terms for Figure 08 ---

significant_terms <- enrichment_results[
  library_fdr < CANDIDATE_LIBRARY_FDR & report_pass == TRUE
]
setorder(significant_terms, source_group, level, gene_set,
         adjusted_fdr, -adjusted_odds_ratio)
plot_data <- significant_terms[, head(.SD, TOP_PER_GROUP),
                               by = .(source_group, level, gene_set)]
plot_data[, source := library_name]
plot_data[, term_label := {
  label <- gsub("^(GOBP_|GOMF_|GOCC_|KEGG_|HP_[0-9]+_)", "", term)
  label <- tools::toTitleCase(tolower(gsub("_", " ", label)))
  label <- gsub("^g Protein Coupled Receptor Activity$", "GPCR activity", label)
  label <- gsub("Oxidative Phosphorylation", "Oxidative phosphorylation", label)
  label <- gsub("Aerobic Respiration", "Aerobic respiration", label)
  label <- gsub("Biosynthetic Process", "biosynthesis", label)
  label <- gsub("Nucleoside Triphosphate", "NTP", label)
  label <- gsub("Atp", "ATP", label)
  label <- gsub("Rna", "RNA", label)
  label <- gsub("Dna", "DNA", label)
  label <- gsub("\\(Hp:", "(HP:", label)
  label <- gsub("Alzheimers Disease", "Alzheimer's disease", label)
  fifelse(nchar(label) <= 48L, label, paste0(substr(label, 1L, 47L), "..."))
}]
plot_data[, term_label := fifelse(
  duplicated(term_label) | duplicated(term_label, fromLast = TRUE),
  paste0(term_label, " [", source, "]"),
  term_label
), by = .(source_group, level)]
plot_data[, neg_log10_fdr := -log10(adjusted_fdr)]
plot_data[, odds_label := sprintf("adjusted OR %.1f", adjusted_odds_ratio)]
plot_data <- plot_data[level == "Directional membership"]

fwrite(plot_data, OUT_PLOT_DATA, sep = "\t")

# --- figure 08A: directional-membership lollipop ---

dt <- copy(plot_data[level == "Directional membership"])
dt[, `:=`(
  gene_set = factor(as.character(gene_set), levels = directional_order),
  source_group = factor(as.character(source_group), levels = source_order)
)]
legend_dt <- data.table(
  x = NA_real_, y = NA_real_,
  gene_set = factor(directional_order, levels = directional_order)
)

out <- copy(dt)
setorder(out, source_group, gene_set, adjusted_fdr, -adjusted_odds_ratio)
out[, y := NA_real_]
y_cursor <- nrow(out) + length(source_order) - 1
for (block in source_order) {
  idx <- which(as.character(out$source_group) == block)
  if (length(idx) == 0L) next
  out$y[idx] <- seq(y_cursor, by = -1, length.out = length(idx))
  y_cursor <- y_cursor - length(idx) - 1
}
dt <- out

bands <- dt[, .(
  ymin = min(y) - 0.5,
  ymax = max(y) + 1.50,
  label_y = max(y) + 1.22,
  n_terms = .N
), by = source_group]
x_max <- ceiling(max(dt$neg_log10_fdr, na.rm = TRUE) + 2.0)

p_membership <- ggplot() +
  geom_rect(data = bands,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax,
                fill = source_group),
            alpha = 0.68, inherit.aes = FALSE) +
  geom_segment(data = dt,
               aes(x = 0, xend = neg_log10_fdr, y = y, yend = y,
                   color = gene_set, linetype = gene_set),
               linewidth = 0.72, lineend = "round") +
  geom_point(data = dt,
             aes(x = neg_log10_fdr, y = y, color = gene_set,
                 shape = gene_set),
             fill = "white", size = 3.6,
             stroke = 0.75) +
  geom_point(data = legend_dt,
             aes(x = x, y = y, color = gene_set, shape = gene_set),
             fill = "white", size = 3.6,
             stroke = 0.75, na.rm = TRUE, inherit.aes = FALSE) +
  geom_text(data = dt,
            aes(x = neg_log10_fdr + 0.10, y = y, label = odds_label),
            hjust = 0, size = 3.6, color = "#111827", fontface = "bold") +
  geom_text(data = bands,
            aes(x = 0.985 * x_max, y = label_y, label = source_group),
            hjust = 1, size = 4.0, color = "#374151",
            fontface = "bold", inherit.aes = FALSE) +
  scale_y_continuous(breaks = dt$y, labels = dt$term_label,
                     expand = expansion(mult = c(0.025, 0.07))) +
  scale_x_continuous(limits = c(0, x_max),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_color_manual(values = group_colors, breaks = directional_order,
                     limits = directional_order, drop = FALSE, name = NULL,
                     guide = guide_legend(
                       override.aes = list(size = 3.6, linewidth = 1.0)
                     )) +
  scale_shape_manual(values = group_shapes, breaks = directional_order,
                     limits = directional_order, drop = FALSE, name = NULL) +
  scale_linetype_manual(values = group_linetypes, breaks = directional_order,
                        limits = directional_order, drop = FALSE, name = NULL) +
  scale_fill_manual(values = source_band_colors, guide = "none") +
  labs(x = expression(-log[10](adjusted~FDR)), y = NULL,
       title = "A  Directional membership") +
  theme_pub(base = 12.0) +
  theme(
    panel.grid.major.x = element_line(color = "#E5E7EB", linewidth = 0.30),
    panel.grid.major.y = element_line(color = "#F2F2F2", linewidth = 0.22),
    axis.text.y = element_text(size = 10.8, color = "black", face = "bold"),
    axis.text.x = element_text(size = 10.6, color = "black", face = "bold"),
    axis.title.x = element_text(size = 11.6, face = "bold"),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left",
    legend.key.size = unit(0.45, "cm"),
    legend.text = element_text(size = 9.8, face = "bold"),
    plot.title = element_text(size = 13.2, face = "bold", hjust = 0),
    plot.margin = margin(5, 15, 5, 5)
  )

p08 <- p_membership

ggsave(OUT_PNG, p08, width = 11.6, height = 6.5, dpi = 300, bg = "white")
ggsave(OUT_PDF, p08, width = 11.6, height = 6.5, bg = "white")

# --- layout validation ---

expected_groups <- data.table(
  level = "Directional membership",
  gene_set = directional_order
)
layout_validation <- plot_data[, .(
  n_terms = .N,
  max_label_chars = max(nchar(term_label)),
  max_neg_log10_fdr = max(neg_log10_fdr),
  max_or_label_chars = max(nchar(odds_label))
), by = .(level = as.character(level), gene_set = as.character(gene_set))]
layout_validation <- merge(expected_groups, layout_validation,
                           by = c("level", "gene_set"), all.x = TRUE)
layout_validation[is.na(n_terms), `:=`(
  n_terms = 0L,
  max_label_chars = 0L,
  max_neg_log10_fdr = NA_real_,
  max_or_label_chars = 0L
)]
layout_validation[, `:=`(
  figure_width_in = 11.6,
  figure_height_in = 6.5,
  group_directly_labeled = FALSE,
  empty_group_displayed = FALSE,
  source_band_labeled_once = TRUE,
  point_label_content = "odds ratio only",
  redundant_group_encoding = "color+shape+linetype",
  row_overlap_flag = anyDuplicated(plot_data[, .(source_group, term_label)]) > 0L,
  right_label_clearance = ceiling(max(plot_data$neg_log10_fdr, na.rm = TRUE) + 2.0) -
    max(plot_data$neg_log10_fdr, na.rm = TRUE),
  validation_note = "Single directional-membership panel; all three groups use one compact legend; source bands are labeled once; row labels are unique within source bands and OR labels retain positive x-axis clearance."
)]
fwrite(layout_validation, OUT_LAYOUT, sep = "\t")
