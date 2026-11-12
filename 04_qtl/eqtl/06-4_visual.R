# Candidate landscape (mediation FDR vs coloc PP4) and case locus figures.
# Writes figures/04_qtl/06-1_landscape and 06-2_case_locus (png and pdf).

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(jsonlite)
})

# --- config ---

APA_ROOT <- Sys.getenv("APA_ROOT")
res_dir <- file.path(APA_ROOT, "data", "04_qtl")
fig_dir <- file.path(APA_ROOT, "figures", "04_qtl")
gwas_coloc_dir <- file.path(APA_ROOT, "data", "07_coloc")
gtf_path <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")

star_genes <- c("NPIPA5", "MUC20", "GATD3", "RPS26", "CKMT1A")

fig_w15 <- 5.35
fig_w2 <- 7.09
png_dpi <- 600
base_font <- 10
font_sm <- 9

col_bg <- "#D9D9D9"
col_focus <- "#1F4E5F"
col_gwas <- "#D55E00"
col_lead <- "#A020F0"
col_gene <- "#4A7BA6"

# LD palette: classic locuscompare/locuszoom colours (as 04_coloc plots/v2)
ld_cols <- c(
  "<0.2" = "#1F2E8C",
  "0.2-0.4" = "#A6CEE3",
  "0.4-0.6" = "#33A02C",
  "0.6-0.8" = "#FF9F1C",
  "0.8-1.0" = "#E31A1C",
  "missing LD" = "#D3D3D3"
)

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ---

as_flag <- function(x) {
  y <- as.logical(x)
  fifelse(is.na(y), FALSE, y)
}

theme_06 <- function() {
  theme_classic(base_size = base_font) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = font_sm),
      axis.title = element_text(colour = "black", size = base_font),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      legend.title = element_text(size = font_sm),
      legend.text = element_text(size = font_sm),
      legend.key.size = grid::unit(0.32, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base_font + 1, face = "bold"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

save_plot <- function(plot, stem, width, height) {
  png_file <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
  unlink(c(png_file, pdf_file))
  ggsave(png_file, plot, width = width, height = height,
         dpi = png_dpi, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height,
         device = grDevices::cairo_pdf, bg = "white")
}

regional_panel <- function(y_col, tag) {
  ggplot(locus_points, aes(pos_mb, .data[[y_col]])) +
    geom_point(aes(colour = ld_bin), size = 2.0, stroke = 0) +
    geom_point(data = lead_pts,
               shape = 23, fill = col_lead, colour = "black",
               size = 3.6, stroke = 0.5) +
    geom_text(data = lead_pts, aes(label = snp),
              fontface = "bold", size = 3.2, hjust = 1.12, vjust = -0.7) +
    annotate("text", x = -Inf, y = Inf, label = tag,
             hjust = -0.08, vjust = 1.5, size = 3.4) +
    scale_colour_manual(values = ld_cols, drop = FALSE, name = expression(r^2)) +
    scale_x_continuous(limits = pos_limits, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = NULL) +
    theme_06() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

# --- loading data ---

priority <- fread(file.path(res_dir, "06_candidate.tsv"))
gwas_hits <- fread(file.path(gwas_coloc_dir, "coloc_run", "all_colocalization.csv"))[as.numeric(PP4) >= 0.7]

display_pool <- priority[
  gene_type == "protein_coding" &
    priority_class %in% c("APA_expression_disease", "APA_expression")
]
display_pool[, core_rank_05 := fcase(
  core_status_05 == "Strict core", 1L,
  core_status_05 == "Consensus core", 2L,
  core_status_05 == "Any-seed core", 3L,
  default = 4L
)]
display_pool[, abs_effect := abs(representative_beta_alt_PDUI)]
setorder(display_pool, -delta_utr_preferred, core_rank_05,
         min_fdr_ab_fwd, -max_abs_prop_med, -abs_effect,
         -disease_coloc, -gwas_best_pp4, na.last = TRUE)
display_pool[, display_rank := .I]
top15 <- display_pool[!duplicated(gene_symbol)][1:min(15L, .N)]
top15[, display_rank := seq_len(.N)]

priority[, key := paste(stage, phenotype_id, gene_symbol, sep = "\t")]
top15[, key := paste(stage, phenotype_id, gene_symbol, sep = "\t")]

priority[, plot_pph4 := fifelse(is.na(max_PP_H4), 0, pmax(0, pmin(1, max_PP_H4)))]
priority[, plot_med := fifelse(is.na(min_fdr_ab_fwd), 0, -log10(pmax(min_fdr_ab_fwd, 1e-300)))]
# only candidates tested by mediation (untested are not drawn)
landscape_dt <- priority[!is.na(min_fdr_ab_fwd)]
landscape_dt[, tier := fcase(
  as_flag(disease_coloc), "Disease-colocalized",
  as_flag(apa_to_expr), "Expression-mediated",
  default = "Tested background"
)]
landscape_dt[, tier := factor(tier, levels = c("Tested background",
                                               "Expression-mediated",
                                               "Disease-colocalized"))]
tier_cols <- c("Tested background" = col_bg,
               "Expression-mediated" = col_focus,
               "Disease-colocalized" = col_gwas)

star_labels <- priority[key %in% top15[gene_symbol %in% star_genes, key]]
star_labels <- merge(
  star_labels,
  top15[, .(key, display_rank, label = gene_symbol)],
  by = "key",
  all.x = TRUE
)
star_label_offsets <- data.table(
  label = star_genes,
  label_nudge_x = c(-0.010, 0.006, -0.015, 0.010, 0.012),
  label_nudge_y = c(0.75, 0.35, 0.08, -0.20, -0.48)
)
star_labels <- merge(star_labels, star_label_offsets, by = "label", all.x = TRUE)
star_labels[, label_order := match(label, star_genes)]
setorder(star_labels, label_order)

# --- figure 06-1: candidate landscape ---

fdr_line <- -log10(0.05)

landscape <- ggplot(landscape_dt, aes(plot_pph4, plot_med)) +
  geom_vline(xintercept = 0.8, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_hline(yintercept = fdr_line, linetype = "dashed",
             colour = "grey70", linewidth = 0.3) +
  geom_point(data = landscape_dt[tier == "Tested background"],
             aes(colour = tier), alpha = 0.45, size = 0.85, stroke = 0) +
  geom_point(data = landscape_dt[tier == "Expression-mediated"],
             aes(colour = tier), alpha = 0.88, size = 1.55, stroke = 0) +
  geom_point(data = landscape_dt[tier == "Disease-colocalized"],
             aes(colour = tier), alpha = 0.95, size = 2.1, stroke = 0) +
  geom_text_repel(
    data = star_labels,
    aes(label = label),
    colour = "black",
    size = 3.1,
    min.segment.length = 0,
    segment.size = 0.22,
    segment.colour = "grey55",
    box.padding = 0.8,
    point.padding = 0.16,
    nudge_x = star_labels$label_nudge_x,
    nudge_y = star_labels$label_nudge_y,
    force = 3,
    force_pull = 0.4,
    max.overlaps = Inf,
    seed = 6
  ) +
  scale_colour_manual(values = tier_cols, name = NULL) +
  scale_x_continuous("Coloc PP.H4 (apaQTL vs eQTL)",
                     limits = c(0.75, 1.02),
                     breaks = seq(0.75, 1.00, by = 0.05),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous("Mediation -log10(FDR)",
                     expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Candidate landscape") +
  theme_06() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(colour = NA, fill = alpha("white", 0.75))) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1)))

save_plot(landscape, "06-1_landscape", fig_w15, 4.1)

# --- case candidate selection ---

top15[, disease_flag := as.integer(as_flag(disease_coloc))]
top15[, strict_flag := as.integer(core_status_05 == "Strict core")]
setorder(top15, -disease_flag, -strict_flag, display_rank,
         min_fdr_ab_fwd, -max_PP_H4, -abs_effect)
case <- top15[1]

case_hits <- gwas_hits[
  GWAS_ID == case$gwas_best_trait &
    QTL_ID == case$stage &
    Phenotype == case$phenotype_id &
    locuscompare_json_file != ""
]
if (nrow(case_hits) == 0L) {
  stop("No locuscompare row for selected case candidate: ", case$gene_symbol)
}
setorder(case_hits, -PP4)
case_hit <- case_hits[1]
json_file <- file.path(gwas_coloc_dir, case_hit$locuscompare_json_file)

locus_json <- jsonlite::fromJSON(json_file)
locus_points <- as.data.table(locus_json$points)
locus_points[, ld_bin := factor(ld_bin, levels = names(ld_cols))]
locus_points[, pos_mb := pos / 1e6]
locus_chr <- locus_points[!is.na(chrom), chrom[1]]
pos_limits <- range(locus_points$pos_mb, na.rm = TRUE)

# --- figure 06-2: case locus ---

lead_pts <- locus_points[is_lead == TRUE]

gwas_panel <- regional_panel("gwas_logp", "GWAS")
qtl_panel  <- regional_panel("qtl_logp", "APA-QTL")

# --- gene track (gencode v40, protein_coding) ---
win_bp <- range(locus_points$pos, na.rm = TRUE)
exon_dt <- fread(
  cmd = sprintf(
    "zcat %s | awk -v c='%s' -v s=%d -v e=%d 'BEGIN{FS=\"\t\";OFS=\"\t\"} $1==c && $3==\"exon\" && $4<=e && $5>=s {if(match($9,/gene_name \"([^\"]+)\"/,a) && match($9,/gene_type \"([^\"]+)\"/,b)) print a[1],$4,$5,$7,b[1]}'",
    shQuote(gtf_path), locus_chr, win_bp[1], win_bp[2]
  ),
  col.names = c("gene_symbol", "start", "end", "strand", "gene_type")
)
if (nrow(exon_dt[gene_type == "protein_coding"]) > 0) {
  exon_dt <- exon_dt[gene_type == "protein_coding"]
}
gene_dt <- exon_dt[, .(start = min(start), end = max(end)),
                   by = gene_symbol][order(start)]
# interval packing: non-overlapping genes (30 kb gap) share a lane
lane_end <- numeric(0)
gene_dt[, lane := 0L]
for (i in seq_len(nrow(gene_dt))) {
  placed <- FALSE
  for (ln in seq_along(lane_end)) {
    if (gene_dt$start[i] - 30000 > lane_end[ln]) {
      gene_dt$lane[i] <- ln
      lane_end[ln] <- gene_dt$end[i]
      placed <- TRUE
      break
    }
  }
  if (!placed) {
    lane_end <- c(lane_end, gene_dt$end[i])
    gene_dt$lane[i] <- length(lane_end)
  }
}
exon_dt <- merge(exon_dt, gene_dt[, .(gene_symbol, lane)], by = "gene_symbol")

# dense loci: greedy label placement by estimated text width; named genes before ENSG*
char_kb <- 20
gene_dt[, `:=`(mid = (start + end) / 2,
               half_w = nchar(gene_symbol) * char_kb * 1000 / 2,
               is_ensg = grepl("^ENSG", gene_symbol),
               is_case = gene_symbol == case$gene_symbol,
               show_label = FALSE)]
setorder(gene_dt, is_ensg, mid)
last_right <- -Inf
for (i in which(!gene_dt$is_case)) {
  left <- gene_dt$mid[i] - gene_dt$half_w[i]
  if (left > last_right + 40000) {
    gene_dt$show_label[i] <- TRUE
    last_right <- max(last_right, gene_dt$mid[i] + gene_dt$half_w[i])
  }
}
case_i <- which(gene_dt$is_case)
case_left <- gene_dt$mid[case_i] - gene_dt$half_w[case_i] - 20000
case_right <- gene_dt$mid[case_i] + gene_dt$half_w[case_i] + 20000
gene_dt$show_label[!gene_dt$is_case &
                     gene_dt$mid - gene_dt$half_w < case_right &
                     gene_dt$mid + gene_dt$half_w > case_left] <- FALSE
gene_dt$show_label[case_i] <- TRUE

gene_panel <- ggplot(gene_dt) +
  geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = lane, yend = lane),
               colour = col_gene, linewidth = 0.4) +
  geom_rect(data = exon_dt,
            aes(xmin = start / 1e6, xmax = end / 1e6,
                ymin = lane - 0.16, ymax = lane + 0.16),
            colour = col_gene, fill = col_gene, linewidth = 0.3) +
  geom_text(data = gene_dt[show_label == TRUE],
            aes(x = mid / 1e6, y = lane + 0.45, label = gene_symbol),
            fontface = "italic", size = 2.9, colour = col_gene) +
  scale_x_continuous(sprintf("Genomic position on %s (Mb)", locus_chr),
                     limits = pos_limits, expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0.3, max(gene_dt$lane) + 0.9),
                     expand = expansion(mult = c(0, 0))) +
  labs(y = NULL) +
  theme_06() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

# --- left: GWAS vs APA-QTL scatter ---
compare_plot <- ggplot(locus_points, aes(gwas_logp, qtl_logp)) +
  geom_point(aes(colour = ld_bin), size = 2.0, stroke = 0) +
  geom_point(data = lead_pts,
             shape = 23, fill = col_lead, colour = "black",
             size = 3.8, stroke = 0.5) +
  geom_text(data = lead_pts, aes(label = snp),
            fontface = "bold", size = 3.2, hjust = 1.12, vjust = -0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.35,
           size = 3.0, colour = "grey20",
           label = sprintf("APA-GWAS PP4 = %.3f\nAPA-eQTL PP.H4 = %.3f",
                           as.numeric(case_hit$PP4), as.numeric(case$max_PP_H4))) +
  scale_colour_manual(values = ld_cols, drop = FALSE, name = expression(r^2)) +
  labs(x = "GWAS −log10(P)", y = "APA-QTL −log10(P)") +
  theme_06() +
  theme(legend.position = c(0.97, 0.03),
        legend.justification = c(1, 0),
        legend.background = element_rect(colour = NA, fill = alpha("white", 0.75)),
        legend.title = element_text(size = font_sm, face = "italic"))

right_col <- (gwas_panel / qtl_panel / gene_panel) +
  plot_layout(heights = c(1, 1, 0.5))

case_fig <- (compare_plot | right_col) +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = sprintf("%s | %s | %s",
                    gsub("_", " ", case$gwas_best_trait), case$stage, case$gene_symbol),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                            size = base_font + 2))
  )

save_plot(case_fig, "06-2_case_locus", fig_w2, 5.2)
