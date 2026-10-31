# Coloc summary per life stage: PP4 definition, class composition, LocusZoom loci selection,
# and figures 1-3. Reads coloc_*_apa_eqtl.tsv; writes data/04_qtl and figures/04_qtl.
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
OUT_DIR <- file.path(APA_ROOT, "figures", "04_qtl")
GTF_FILE <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
STAGES <- c("prenatal", "postnatal")
PP_CALL <- 0.8
PP_REPORT <- 0.7
LZ_MAX_LOCI <- 12
FIG_W15 <- 5.35
FIG_W2 <- 7.09
FIG_HMAX <- 9.45
BASE_FONT_SIZE <- 8
FONT_SM <- 7
TITLE_FONT_SIZE <- 9
PNG_DPI <- 600
CLASS_LEVELS <- c("H4 shared causal variant", "H3 distinct causal variants", "H2 APA signal only",
                  "H1 eQTL signal only", "H0 no signal", "undetermined")
DETERMINED <- CLASS_LEVELS[1:5]
CLASS_COL <- c("H4 shared causal variant" = "#4C72B0", "H3 distinct causal variants" = "#C44E52",
               "H2 APA signal only" = "#DD8452", "H1 eQTL signal only" = "#8C8C8C",
               "H0 no signal" = "#CCCCCC", "undetermined" = "#EDEDED")
DIR_COL <- c(concordant = "#009E73", discordant = "#CC79A7")
pt_mm <- function(pt) pt / .pt

# --- helpers ---
theme_sci <- function(base_size = BASE_FONT_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black", size = FONT_SM),
      axis.title = element_text(colour = "black", size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      legend.title = element_text(size = FONT_SM),
      legend.text = element_text(size = FONT_SM),
      legend.key.size = grid::unit(0.30, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      strip.text = element_text(face = "bold", size = base_size, margin = margin(b = 1.5)),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold", size = TITLE_FONT_SIZE, hjust = 0),
      plot.margin = margin(2, 3, 2, 2)
    )
}

save_plot <- function(plot, filename, width, height) {
  if (height > FIG_HMAX) warning(sprintf("%s: height %.2f in exceeds the %.2f in page limit", filename, height, FIG_HMAX))
  ggsave(file.path(OUT_DIR, filename), plot, width = width, height = height, dpi = PNG_DPI, bg = "white")
}

# --- coloc results and loci per stage ---
coloc_all <- rbindlist(lapply(STAGES, function(stage) {
  fread(file.path(D, sprintf("coloc_%s_apa_eqtl.tsv", stage)))
}), use.names = TRUE, fill = TRUE)
coloc_all[, stage := factor(stage, levels = STAGES)]
coloc_all[, coloc_class := factor(coloc_class, levels = CLASS_LEVELS)]
fwrite(coloc_all, file.path(D, "coloc_all_stages.tsv"), sep = "\t", quote = FALSE)

locus_best <- coloc_all[order(stage, apa_pheno, -PP_H4)][, .SD[1], by = .(stage, apa_pheno)]

# --- PP4 definition: share of significant APA loci (eQTL-tested genes) with PP4 >= 0.7 ---
gtf <- fread(GTF_FILE, sep = "\t", header = FALSE, select = c(3, 9), quote = "")
gtf <- gtf[V3 == "gene"]
gtf_fields <- tstrsplit(gtf$V9, "\"", fixed = TRUE)
gene_map <- unique(data.table(ensg0 = sub("\\..*", "", gtf_fields[[2]]), gene_name = gtf_fields[[6]]), by = "gene_name")
gene_to_ensg <- setNames(gene_map$ensg0, gene_map$gene_name)

pp4_summary <- rbindlist(lapply(STAGES, function(stage_value) {
  sig_apa <- unique(fread(file.path(D, sprintf("%s.sig_QTL.txt.gz", stage_value)), select = "phenotype_id"))
  sig_apa[, gene_raw := tstrsplit(phenotype_id, "|", fixed = TRUE, keep = 2L)]
  sig_apa[, ensg0 := fifelse(grepl("^ENSG", gene_raw), sub("\\..*", "", gene_raw), unname(gene_to_ensg[gene_raw]))]
  eqtl_tested <- unique(fread(file.path(D, sprintf("nominal_eqtl_%s.tsv.gz", stage_value)), select = "phenotype_id"))
  eqtl_tested[, ensg0 := sub("\\..*", "", phenotype_id)]
  denominator <- sig_apa[ensg0 %in% eqtl_tested$ensg0, .(apa_pheno = phenotype_id, ensg0)]
  best <- locus_best[stage == stage_value, .(apa_pheno, PP_H4)]
  audit <- merge(denominator, best, by = "apa_pheno", all.x = TRUE)
  numerator <- sum(audit$PP_H4 >= PP_REPORT, na.rm = TRUE)
  data.table(stage = stage_value, pp4_threshold = PP_REPORT, numerator = numerator,
             denominator = nrow(audit), fraction = numerator / nrow(audit),
             n_with_coloc_result = sum(!is.na(audit$PP_H4)), n_without_coloc_result = sum(is.na(audit$PP_H4)))
}))
fwrite(pp4_summary, file.path(D, "coloc_pp4_definition_summary.tsv"), sep = "\t", quote = FALSE)

# --- class composition ---
class_counts <- locus_best[, .(n = .N), by = .(stage, coloc_class)]
stage_tot <- locus_best[, .(n_tested = .N, n_determinable = sum(coloc_class %in% DETERMINED)), by = stage]
class_summary <- merge(class_counts, stage_tot, by = "stage")
class_summary[, frac_tested := n / n_tested]
class_summary[, frac := fifelse(coloc_class %in% DETERMINED, n / n_determinable, NA_real_)]
setorder(class_summary, stage, coloc_class)
fwrite(class_summary, file.path(D, "coloc_class_summary.tsv"), sep = "\t", quote = FALSE)

# --- LocusZoom selection: shared H4, independent H3, discordant, and stage-switch loci ---
pick_loci <- function(lb) {
  take <- function(dt, n, tag) {
    if (nrow(dt) == 0) return(NULL)
    dt <- dt[!duplicated(ensg0)]
    out <- head(dt, n)
    out[, lz_tag := tag]
    out[]
  }
  sel <- list()
  for (st in STAGES) {
    d <- lb[stage == st]
    sel[[length(sel) + 1]] <- take(d[PP_H4 >= PP_CALL][order(-PP_H4)], 2, "shared-H4")
    sel[[length(sel) + 1]] <- take(d[PP_H3 >= PP_CALL][order(-PP_H3)], 2, "independent-H3")
    sel[[length(sel) + 1]] <- take(d[PP_H4 >= PP_CALL & direction == "discordant"][order(-PP_H4)], 1, "discordant-direction")
  }
  gw <- lb[, .(PP_H4 = max(PP_H4), PP_H3 = max(PP_H3)), by = .(stage, ensg0)]
  wide <- dcast(gw, ensg0 ~ stage, value.var = c("PP_H4", "PP_H3"), fill = 0)
  need <- c("PP_H4_prenatal", "PP_H4_postnatal", "PP_H3_prenatal", "PP_H3_postnatal")
  if (all(need %in% names(wide))) {
    sw <- wide[(PP_H4_prenatal >= PP_CALL & PP_H3_postnatal >= PP_CALL) | (PP_H4_postnatal >= PP_CALL & PP_H3_prenatal >= PP_CALL)]
    if (nrow(sw) > 0) {
      g <- sw[order(-pmax(PP_H4_prenatal, PP_H4_postnatal))][1]$ensg0
      d <- lb[ensg0 == g][order(stage, -pmax(PP_H4, PP_H3))]
      d <- d[!duplicated(stage)]
      d[, lz_tag := "stage-switch"]
      sel[[length(sel) + 1]] <- d
    }
  }
  out <- rbindlist(Filter(Negate(is.null), sel), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0) return(out)
  out <- out[!duplicated(paste(stage, apa_pheno, lz_tag))]
  head(out, LZ_MAX_LOCI)
}
selected <- pick_loci(locus_best)
fwrite(selected, file.path(D, "coloc_locuszoom_selection.tsv"), sep = "\t", quote = FALSE)

# --- figure 1: class composition (fractions of determinable loci) and determinable rate ---
plot_dt <- class_summary[coloc_class %in% DETERMINED & n > 0]
plot_dt[, coloc_class := factor(as.character(coloc_class), levels = rev(intersect(CLASS_LEVELS, unique(coloc_class))))]
plot_dt[, stage_lab := sprintf("%s\n(n = %s)", stage, format(n_determinable, big.mark = ","))]
plot_dt[, stage_lab := factor(stage_lab, levels = unique(stage_lab[order(stage)]))]
fig1a <- ggplot(plot_dt, aes(stage_lab, frac, fill = coloc_class)) +
  geom_col(width = 0.58, colour = "white", linewidth = 0.25) +
  geom_text(data = plot_dt[frac >= 0.08],
            aes(label = sprintf("%s\n(%.1f%%)", format(n, big.mark = ","), 100 * frac)),
            position = position_stack(vjust = 0.5), size = pt_mm(FONT_SM), colour = "black", lineheight = 0.9) +
  scale_fill_manual(values = CLASS_COL, breaks = CLASS_LEVELS, name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.02))) +
  labs(x = NULL, y = "Fraction of determinable APA-QTL loci") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_sci() + theme(legend.position = "right")

det_dt <- copy(stage_tot)
det_dt[, frac := n_determinable / n_tested]
det_dt[, stage := factor(stage, levels = STAGES)]
fig1b <- ggplot(det_dt, aes(stage, frac)) +
  geom_col(width = 0.5, fill = "grey55", colour = "grey25", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%s / %s", format(n_determinable, big.mark = ","), format(n_tested, big.mark = ","))),
            vjust = -0.4, size = pt_mm(FONT_SM)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = expansion(c(0, 0.06))) +
  labs(x = NULL, y = "Loci reaching a call") +
  theme_sci()
fig1 <- (fig1a | fig1b) + plot_layout(widths = c(1, 0.42)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = TITLE_FONT_SIZE))
save_plot(fig1, "fig1_coloc_class_composition.png", width = FIG_W2, height = 2.9)

# --- figure 2: PP(H3) vs PP(H4) ---
fig2 <- ggplot(locus_best, aes(PP_H4, PP_H3, colour = coloc_class)) +
  geom_point(alpha = 0.55, size = 0.6, stroke = 0) +
  geom_hline(yintercept = PP_CALL, linetype = 2, colour = "grey40", linewidth = 0.25) +
  geom_vline(xintercept = PP_CALL, linetype = 2, colour = "grey40", linewidth = 0.25) +
  scale_colour_manual(values = CLASS_COL, breaks = CLASS_LEVELS, name = NULL) +
  facet_wrap(~ stage) +
  labs(x = "PP(H4) shared causal variant", y = "PP(H3) distinct causal variants") +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 2)) +
  theme_sci() + theme(legend.position = "bottom")
save_plot(fig2, "fig2_coloc_pp_landscape.png", width = FIG_W15, height = 3.1)

# --- figure 3: eQTL vs APA effect direction at PP(H4) classes ---
quad <- locus_best[coloc_class == CLASS_LEVELS[1] & is.finite(h4_slope_eqtl) & is.finite(h4_slope_apa)]
if (nrow(quad) > 0) {
  lim <- max(abs(c(quad$h4_slope_eqtl, quad$h4_slope_apa))) * 1.10
  ann <- quad[, .(n = .N), by = .(stage, direction)]
  ann[, `:=`(x = ifelse(direction == "concordant", 0.70 * lim, -0.70 * lim), y = 0.93 * lim)]
  lab_dt <- quad[order(stage, -abs(h4_slope_apa))][, head(.SD, 8), by = stage]
  fig3 <- ggplot(quad, aes(h4_slope_eqtl, h4_slope_apa, colour = direction)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.25) +
    geom_point(alpha = 0.7, size = 0.8, stroke = 0) +
    geom_text_repel(data = lab_dt, aes(label = gene_name), size = pt_mm(FONT_SM), max.overlaps = 20,
                    min.segment.length = 0.2, segment.size = 0.2, box.padding = 0.18,
                    segment.colour = "grey60", show.legend = FALSE) +
    geom_text(data = ann, aes(x, y, label = sprintf("%s\nn = %s", direction, n)),
              size = pt_mm(FONT_SM), lineheight = 0.9, show.legend = FALSE) +
    scale_colour_manual(values = DIR_COL, name = NULL) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    facet_wrap(~ stage) +
    labs(x = "eQTL effect on expression (slope per alt allele)",
         y = "APA-QTL effect on 3'UTR usage (slope per alt allele)") +
    guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
    theme_sci() + theme(legend.position = "bottom")
  save_plot(fig3, "fig3_effect_direction_quadrant.png", width = FIG_W15, height = 3.2)
}
