#!/usr/bin/env Rscript
# Fig 3b: PDUI effects (slope +/- 1.96 SE) of apaQTL variants that disrupt a canonical or variant poly(A) hexamer,
# with a directional prediction. Hexamer annotation: pas.tsv.gz; slopes: {prenatal,postnatal}.sig_QTL.txt.gz.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

APA_ROOT <- Sys.getenv("APA_ROOT")
if (!nzchar(APA_ROOT)) stop("Set APA_ROOT (source config/paths.sh)")
D        <- file.path(APA_ROOT, "data", "04_qtl")
FIG_OUT  <- file.path(APA_ROOT, "figures", "04_qtl")

source(file.path(APA_ROOT, "release", "figures", "theme.R"))

# --- helpers ---
# one panel per predicted direction
mk <- function(d, ttl, ANN_X, ANN_Y, ANN_H, ANN_V) {
  d <- copy(d)
  setorder(d, slope)
  d[, row_id := factor(row_id, levels = row_id)]
  ggplot(d, aes(slope, row_id, colour = stage, shape = concordant)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 0.25) +
    geom_point(size = 1.4, stroke = 0.4, fill = "white") +
    scale_colour_manual(values = STAGE_COLORS[c("Prenatal", "Postnatal")], name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c("Discordant", "Concordant"), name = NULL) +
    scale_x_continuous(limits = xr, breaks = c(-1, 0, 1)) +
    scale_y_discrete(labels = function(x) trimws(x)) +
    labs(x = NULL, y = NULL) +
    annotate("text", x = ANN_X, y = ANN_Y, label = ttl, hjust = ANN_H, vjust = ANN_V, size = 5.5 / .pt, colour = INK, lineheight = 0.95) +
    guides(shape = guide_legend(override.aes = list(colour = INK))) +
    theme(axis.text.y = element_text(face = "italic", size = 5.5))
}

# --- hexamer-disrupting candidates with a predicted PDUI direction ---
pas <- fread(cmd = paste("zcat", file.path(D, "pas.tsv.gz")))
pas[, stage := ifelse(grepl("^prenatal", candidate_id), "Prenatal", "Postnatal")]
pas[, phen := sub("^[a-z]+::(.*)::chr[^:]+:[0-9]+:[ACGT]+:[ACGT]+$", "\\1", candidate_id)]
pas[, variant_id := sub("^.*::(chr[^:]+:[0-9]+:[ACGT]+:[ACGT]+)$", "\\1", candidate_id)]
h <- pas[predicted_pdui_direction %in% c("lengthening", "shortening") & pas_core_change != "not_eligible"]
h <- unique(h[, .(stage, phen, variant_id, pas_label, pas_core_change, predicted_pdui_direction)])

# --- slopes of the candidates in the stage-specific apaQTL maps ---
sq <- rbindlist(lapply(c("Prenatal", "Postnatal"), function(st) {
  f <- if (st == "Prenatal") "prenatal.sig_QTL.txt.gz" else "postnatal.sig_QTL.txt.gz"
  x <- fread(cmd = paste("zcat", file.path(D, f)), select = c("phenotype_id", "variant_id", "slope", "slope_se", "pval_nominal"))
  x[, stage := st]; x
}))
m <- merge(h, sq, by.x = c("stage", "phen", "variant_id"), by.y = c("stage", "phenotype_id", "variant_id"), all.x = TRUE)
stopifnot(all(!is.na(m$slope)))
m[, gene := tstrsplit(phen, "|", fixed = TRUE)[[2]]]
m[, concordant := (predicted_pdui_direction == "lengthening" & slope > 0) | (predicted_pdui_direction == "shortening" & slope < 0)]
m[, lo := slope - 1.96 * slope_se]; m[, hi := slope + 1.96 * slope_se]
setorder(m, predicted_pdui_direction, slope)
m[, row_id := paste0(gene, strrep(" ", seq_len(.N))), by = predicted_pdui_direction]
xr <- range(c(m$lo, m$hi, 0)) + c(-0.1, 0.1)
conc_txt <- sprintf("%d of %d (%.1f%%) concordant", sum(m$concordant), nrow(m), 100 * mean(m$concordant))

# --- panels ---
dl <- m[predicted_pdui_direction == "lengthening"]
ds <- m[predicted_pdui_direction == "shortening"]
pl <- mk(dl, "Predicted lengthening", xr[1] + 0.05, Inf, 0, 1.8)
ps <- mk(ds, paste0("Predicted shortening\n", conc_txt), xr[2] - 0.05, -Inf, 1, -0.4) +
  labs(x = "PDUI effect per alternative allele (95% CI)")
pl <- pl + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
p <- (pl / ps) + plot_layout(heights = c(nrow(dl), nrow(ds)), guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(1, "mm"))
p <- p + plot_annotation(caption = NULL)

# --- save (88 mm = single-column panel width, 80 mm height) ---
dir.create(FIG_OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FIG_OUT, "fig3b_hexamer_pdui.png"), p, width = 88, height = 80, units = "mm", dpi = 600, bg = "white")
