# Fig 5 d-e: regional tracks (GWAS / apaQTL / eQTL) with LD colours and gene model; stage and gene set per call.
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/locus.R"))
# --- helpers ---
locus_panel <- function(gene, gwas_lab, stage = "prenatal", half = 150e3, pp4 = NA, xlim = NULL, min_gap = 30e3) {
  TG <- gene
  D  <- fread(file.path(DATA, sprintf("figures/fig5/fig5_%s_tracks.tsv", gene)))
  GN <- fread(file.path(DATA, sprintf("figures/fig5/fig5_%s_genes.tsv", gene)))
  lp <- D$lead_pos[1]; lsnp <- D$lead_snp[1]; x0 <- lp - half; x1 <- lp + half
  if (!is.null(xlim)) { x0 <- xlim[1]; x1 <- xlim[2] }
  D <- D[pos >= x0 & pos <= x1]
  # intervals without any cohort variant (union of tracks): grey band + dashed bridge, never blank
  up <- sort(unique(c(x0, D$pos, x1))); gi <- which(diff(up) > min_gap)
  GAP <- data.table(gs = up[gi], ge = up[gi + 1])
  D[, ldb := cut(r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = LDB, right = FALSE)]
  TR <- c("GWAS", "apaQTL", "eQTL")
  tl <- c(GWAS = sprintf("%s GWAS", sub("^(EUR|EAS) ", "", gwas_lab)), apaQTL = "apaQTL", eQTL = "eQTL")
  D[, track := factor(track, TR, tl[TR])]
  D[, isref := pos == lp]
  setorder(D, r2)
  top <- D[, .SD[which.max(logp)], by = track]
  top[, lab := fifelse(grepl("^rs", snp), snp, sub("^(chr[^:]+:[0-9]+).*", "\\1", snp))]
  ymax <- D[, .(ym = max(logp) * 1.22 + 0.3), by = track]
  # gene model: protein-coding genes + target, packed in rows
  g <- GN[feat == "gene" & ((type == "protein_coding" & !grepl("^ENSG", gene)) | gene == TG) & end >= x0 & start <= x1]
  g <- unique(g, by = "gene")[order(start)]
  g[, `:=`(s = pmax(start, x0), e = pmin(end, x1))]
  g[, row := 0L]; ends <- numeric(0)
  for (k in seq_len(nrow(g))) {
    lab_end <- g$s[k] + max(g$e[k] - g$s[k], nchar(g$gene[k]) * 4.6e3) + 6e3
    r <- which(ends < g$s[k])[1]
    if (is.na(r)) { ends <- c(ends, lab_end); r <- length(ends) } else ends[r] <- lab_end
    g$row[k] <- r
  }
  ex <- merge(GN[feat == "exon"], g[, .(gene, row)], by = "gene")[end >= x0 & start <= x1]
  g[, tgt := gene == TG]
  base <- ggplot(D, aes(pos / 1e6, logp)) +
    geom_blank(data = ymax, aes(x = lp / 1e6, y = ym), inherit.aes = FALSE) +
    geom_rect(data = GAP, aes(xmin = gs / 1e6, xmax = ge / 1e6, ymin = -Inf, ymax = Inf), fill = "#EEEEEE", inherit.aes = FALSE) +
    geom_segment(data = GAP, aes(x = gs / 1e6, xend = ge / 1e6, y = 0, yend = 0), linetype = "22", linewidth = 0.3, colour = NEUTRAL, inherit.aes = FALSE) +
    geom_vline(xintercept = lp / 1e6, linewidth = 0.25, colour = LEAD_COL, linetype = "dotted") +
    geom_point(data = D[isref == FALSE], aes(fill = ldb), shape = 21, colour = "white", stroke = 0.12, size = 1.05) +
    geom_point(data = D[isref == TRUE], shape = 23, fill = LEAD_COL, colour = "white", stroke = 0.3, size = 2.2) +
    geom_text(data = top, aes(label = lab), size = 5.2 / .pt, vjust = -0.9, colour = INK, fontface = "italic") +
    facet_grid(track ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = LDC, drop = FALSE, name = expression(italic(r)^2), guide = guide_legend(override.aes = list(size = 1.8), nrow = 1)) +
    scale_x_continuous(limits = c(x0, x1) / 1e6, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = scales::breaks_pretty(3), limits = c(0, NA)) +
    labs(x = NULL, y = expression(-log[10] * italic(P))) +
    theme(strip.text.y.left = element_text(size = 6, angle = 90, face = "bold", margin = margin(0, 0.6, 0, 0.6, "mm")),
          strip.placement = "outside", strip.background.y = element_rect(fill = "#EFEFEF", colour = NA),
          panel.spacing.y = unit(1.5, "mm"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "top", legend.justification = c(1, 1), legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.5),
          legend.title = element_text(size = 6), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -1.5, 0, "mm"),
          plot.title = element_text(size = 7, face = "bold.italic", hjust = 0, margin = margin(0, 0, 0.5, 0, "mm")),
          plot.subtitle = element_text(size = 5.8, colour = MUTED, hjust = 0, margin = margin(0, 0, 0.5, 0, "mm")),
          plot.margin = margin(1, 5, 1, 1.5, "mm"))
  nr <- max(g$row)
  gm <- ggplot() +
    geom_segment(data = g, aes(x = s / 1e6, xend = e / 1e6, y = -row, yend = -row, colour = tgt), linewidth = 0.3) +
    geom_rect(data = ex, aes(xmin = pmax(start, x0) / 1e6, xmax = pmin(end, x1) / 1e6, ymin = -row - 0.18, ymax = -row + 0.18,
                             fill = gene == TG), colour = NA) +
    geom_text(data = g, aes(x = s / 1e6, y = -row - 0.48, label = gene, colour = tgt, fontface = ifelse(tgt, "bold.italic", "italic")), size = 5 / .pt, hjust = 0) +
    geom_vline(xintercept = lp / 1e6, linewidth = 0.25, colour = LEAD_COL, linetype = "dotted") +
    scale_colour_manual(values = c(`FALSE` = "#7A7A7A", `TRUE` = APA_COL), guide = "none") +
    scale_fill_manual(values = c(`FALSE` = "#7A7A7A", `TRUE` = APA_COL), guide = "none") +
    scale_x_continuous(limits = c(x0, x1) / 1e6, expand = c(0, 0), breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(limits = c(-nr - 0.8, -0.6), expand = c(0, 0)) +
    labs(x = sprintf("%s position (Mb)", sub("^chr", "Chr ", GN$chr[1])), y = NULL) +
    coord_cartesian(xlim = c(x0, x1) / 1e6, clip = "off") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), plot.margin = margin(0, 5, 1, 1.5, "mm"))
  if (nrow(GAP)) base <- base + geom_text(data = data.table(track = factor(tl[["GWAS"]], tl[TR]), x = (GAP$gs + GAP$ge) / 2e6), aes(x = x, y = Inf, label = "no\nvariants"), lineheight = 0.85,
                                          vjust = 1.6, size = 5 / .pt, colour = MUTED, inherit.aes = FALSE)
  subt <- sprintf("%s \u00D7 %s apaQTL, PP4 = %.2f; \u25C6 apaQTL lead %s", gwas_lab, stage, pp4, lsnp)
  list(top = base + labs(title = gene, subtitle = subt), gm = gm, nr = nr)
}

# --- setup ---
DATA <- file.path(ROOT, "data")
LDB <- c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8")
LDC <- setNames(unname(LOCUS_LD_COLS[c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")]), LDB)
LEAD_COL <- "#6A4C93"
