# Fig 4 composite: G x stage classes (a), exemplar PDUI boxes (b), adult-map detection (c), mechanism schematic (d).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(readxl); library(ggrepel); library(png)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
# --- helpers ---
negu <- function(x) sub("-", "\u2212", x, fixed = TRUE)
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1))
mk_b <- function(g, cl) {
  d <- copy(B4[gene == g]); lv <- d[order(dosage), unique(genotype)]
  d[, gt := factor(genotype, lv)][, xpos := as.integer(gt) + fifelse(stage == "Prenatal", -0.19, 0.19)]
  ol <- d[!is.na(outliers) & nzchar(outliers), .(v = as.numeric(strsplit(outliers, ",")[[1]])), by = .(stage, gt, xpos)]
  rng <- range(c(d$whisker_lo, d$whisker_hi)); pad <- diff(rng); bb <- gx[gene_id == g]
  ggplot(d) +
    geom_segment(aes(x = xpos, xend = xpos, y = whisker_lo, yend = whisker_hi, colour = stage), linewidth = 0.3) +
    geom_rect(aes(xmin = xpos - 0.16, xmax = xpos + 0.16, ymin = q1, ymax = q3, fill = stage, colour = stage), linewidth = 0.3, alpha = 0.75) +
    geom_segment(aes(x = xpos - 0.16, xend = xpos + 0.16, y = median, yend = median), colour = "white", linewidth = 0.5) +
    geom_point(data = ol, aes(xpos, v, colour = stage), size = 0.3, stroke = 0) +
    geom_line(aes(xpos, mean, colour = stage, group = stage), linewidth = 0.4) +
    geom_point(aes(xpos, mean, colour = stage), shape = 21, fill = "white", size = 0.9, stroke = 0.4) +
    geom_text(aes(x = as.integer(gt), y = rng[1] - pad * fifelse(stage == "Prenatal", 0.08, 0.17), label = n, colour = stage), size = 5.2 / .pt) +
    scale_colour_manual(values = STAGE_TEXT, guide = "none") +
    scale_fill_manual(values = SC, guide = "none") +
    scale_x_continuous(breaks = 1:3, labels = lv, expand = expansion(add = 0.4)) +
    scale_y_continuous(labels = function(x) sub("-", "\u2212", x)) +
    coord_cartesian(ylim = c(rng[1] - pad * 0.22, rng[2] + pad * 0.04), clip = "off") +
    labs(x = NULL, y = NULL, title = bquote(bolditalic(.(g)) ~ bold(.(cl))),
         subtitle = bquote(beta[pre] == .(sub("-", "\u2212", sprintf("%.2f", bb$beta_pre))) * "," ~ beta[post] == .(sub("-", "\u2212", sprintf("%.2f", bb$beta_post))))) +
    theme(plot.title = element_text(size = 6.5, colour = INK, hjust = 0, margin = margin(0, 0, 0.5, 0)),
          plot.subtitle = element_text(size = 6, colour = MUTED, hjust = 0, margin = margin(0, 0, 1, 0)),
          axis.text.x = element_text(size = 6), plot.margin = margin(1, 1.2, 1, 1.2, "mm"))
}

# --- setup ---
D4 <- file.path(ROOT, "data/04_qtl")
DF <- file.path(ROOT, "data/figures/fig4")
theme_set(theme_nature(base_size = 7) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), axis.ticks.length = unit(0.8, "mm"),
  plot.margin = margin(1, 1.5, 1, 1.5, "mm")))
T6 <- 6 / .pt
CLS_ORDER <- c("Reversed", "Amplified", "Attenuated")
CC <- CLASS_COLORS[CLS_ORDER]
TINT <- CLASS_TINT[CLS_ORDER]
SC <- STAGE_COLORS[c("Prenatal", "Postnatal")]

# --- a: prenatal vs postnatal effect per gene, classes by sign and magnitude ---
gx <- fread(file.path(D4, "gxstage_strict_classes_local_checks.tsv"))
gx[, cls := fifelse(sign(beta_pre) != sign(beta_post), "Reversed", fifelse(abs(beta_post) > abs(beta_pre), "Amplified", "Attenuated"))]
stopifnot(nrow(gx) == 193, all(table(gx$cls)[CLS_ORDER] == c(122, 40, 31)))
gx[, cls := factor(cls, CLS_ORDER)]
L <- 1.9  # axis half-width; CI bounds reach +1.82 / -1.79
lab_genes <- c("PRDX6", "MIR124-1HG", "ZEB2", "ATP2B1", "TESPA1", "SLCO1A2")
hl <- gx[gene_id %in% lab_genes]; bg <- gx[!gene_id %in% lab_genes]
badge <- data.table(cls = factor(CLS_ORDER, CLS_ORDER), x = c(-L + 0.06, 0.06, L - 0.06), y = c(L - 0.06, L - 0.06, 0.62),
                    hj = c(0, 0, 1), vj = c(1, 1, 0), lab = sprintf("%s  %d", CLS_ORDER, as.integer(table(gx$cls)[CLS_ORDER])))
reg <- rbind(
  data.table(cls = "Reversed",   grp = "q2", x = c(-L, 0, 0, -L), y = c(0, 0, L, L)),
  data.table(cls = "Reversed",   grp = "q4", x = c(0, L, L, 0),   y = c(-L, -L, 0, 0)),
  data.table(cls = "Amplified",  grp = "a1", x = c(0, L, 0),      y = c(0, L, L)),
  data.table(cls = "Amplified",  grp = "a3", x = c(0, -L, 0),     y = c(0, -L, -L)),
  data.table(cls = "Attenuated", grp = "t1", x = c(0, L, L),      y = c(0, 0, L)),
  data.table(cls = "Attenuated", grp = "t3", x = c(0, -L, -L),    y = c(0, 0, -L)))
reg[, cls := factor(cls, CLS_ORDER)]
pa <- ggplot(gx, aes(beta_pre, beta_post)) +
  geom_polygon(data = reg, aes(x, y, group = grp, fill = cls), colour = NA, alpha = 0.7, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = MUTED) + geom_vline(xintercept = 0, linewidth = 0.25, colour = MUTED) +
  geom_abline(slope = c(1, -1), intercept = 0, linewidth = 0.2, linetype = "22", colour = NEUTRAL) +
  geom_segment(data = bg, aes(x = hc3_pre_ci95_low, xend = hc3_pre_ci95_high, yend = beta_post, colour = cls), linewidth = 0.15, alpha = 0.1) +
  geom_segment(data = bg, aes(y = hc3_post_ci95_low, yend = hc3_post_ci95_high, xend = beta_pre, colour = cls), linewidth = 0.15, alpha = 0.1) +
  geom_point(data = bg, aes(colour = cls), size = 0.7, shape = 16, alpha = 0.6) +
  geom_segment(data = hl, aes(x = hc3_pre_ci95_low, xend = hc3_pre_ci95_high, yend = beta_post), linewidth = 0.35, colour = INK) +
  geom_segment(data = hl, aes(y = hc3_post_ci95_low, yend = hc3_post_ci95_high, xend = beta_pre), linewidth = 0.35, colour = INK) +
  geom_point(data = hl, aes(colour = cls), shape = 21, fill = "white", size = 1.6, stroke = 0.6) +
  geom_text_repel(data = hl, aes(label = gene_id), size = T6, fontface = "italic", family = FONT, colour = INK, min.segment.length = 0,
                  segment.size = 0.2, box.padding = 0.3, point.padding = 0.25, seed = 7, max.overlaps = 50) +
  geom_text(data = badge, aes(x, y, label = lab, colour = cls, hjust = hj, vjust = vj), size = T6, fontface = "bold", family = FONT, inherit.aes = FALSE) +
  scale_colour_manual(values = CC, guide = "none") + scale_fill_manual(values = TINT, guide = "none") +
  coord_equal(xlim = c(-L, L), ylim = c(-L, L), expand = FALSE) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), labels = function(x) negu(as.character(x))) + scale_y_continuous(breaks = seq(-1, 1, 0.5), labels = function(x) negu(as.character(x))) +
  labs(x = expression("Prenatal effect " * beta[pre]), y = expression("Postnatal effect " * beta[post]))

# --- b: exemplar genes, stage-overlaid boxplots of adjusted PDUI (summary box statistics) ---
B4 <- fread(file.path(DF, "fig4b_boxstats.tsv")); B4[, stage := factor(stage, c("Prenatal", "Postnatal"))]
stopifnot(all(fread(file.path(DF, "fig4b_signcheck.tsv"))$same_sign))
pb1 <- mk_b("PRDX6", "Reversed"); pb2 <- mk_b("MIR124-1HG", "Amplified"); pb3 <- mk_b("ZEB2", "Attenuated")
leg <- ggplot(data.table(x = c(1, 2.1), s = factor(c("Prenatal", "Postnatal"), c("Prenatal", "Postnatal")))) +
  geom_tile(aes(x, 1, fill = s), width = 0.16, height = 0.55, alpha = 0.75) + geom_text(aes(x + 0.14, 1, label = s), hjust = 0, size = T6) +
  scale_fill_manual(values = SC, guide = "none") +
  annotate("text", x = 3.2, y = 1, hjust = 0, size = T6, colour = MUTED, label = "box, IQR; circle, mean; numbers, n per genotype") +
  coord_cartesian(xlim = c(0.9, 6.2), ylim = c(0.5, 1.5)) + theme_void()
ylab_b <- wrap_elements(full = grid::textGrob("Adjusted PDUI", rot = 90, gp = grid::gpar(fontsize = 7, fontfamily = FONT, col = INK)))
pb <- (ylab_b | ((pb1 | pb2 | pb3) / leg + plot_layout(heights = c(1, 0.07)))) + plot_layout(widths = c(0.035, 1))

# --- d (mechanism schematic from fig4_panel_c.py) ---
cimg <- readPNG(file.path(ROOT, "figures/fig4_c.png"))
pc <- wrap_elements(full = grid::rasterGrob(cimg, interpolate = TRUE)) + labs(tag = "d") +
  theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1))

# --- c: adult maps, stacked detection bars (bulk PFC and aged snRNA-seq) ---
bgv <- as.data.table(read_excel(file.path(D4, "braingvex_adult_bulk_comp.xlsx"), sheet = "GxStage_class"))
bgv <- bgv[cls %in% CLS_ORDER, .(cls, n = n_genes, k = round(lead_pair_detected * n_genes))]
ex <- fread(file.path(D4, "R4_intx_nan_exact_pairs.tsv"))[matched == TRUE | matched == "TRUE"]
ex <- merge(ex, gx[, .(local_id = gene_id, cls)], by = "local_id")
sn <- ex[, .(n = .N, k = sum(replicated == TRUE | replicated == "TRUE")), by = cls]
E <- rbind(cbind(bgv, map = "Adult bulk PFC"), cbind(sn, map = "Aged snRNA-seq"))
E[, cls := factor(as.character(cls), CLS_ORDER)]; E[, map := factor(map, c("Adult bulk PFC", "Aged snRNA-seq"))]
S <- rbind(E[, .(map, cls, part = "det", v = 100 * k / n)], E[, .(map, cls, part = "not", v = 100 * (n - k) / n)])
S[, part := factor(part, c("not", "det"))]
S[, fillkey := ifelse(part == "det", as.character(cls), paste0(cls, "_t"))]
FILLS <- c(setNames(scales::alpha(CC, 0.85), names(CC)), setNames(TINT, paste0(names(TINT), "_t")))
pexact <- fread(file.path(DF, "r4_fig4e_exact_match_models.tsv"))$P[1]
# adult-bulk bracket: workbook row 'Reversed vs other classes, lead pair detected' (adjusted logistic, postnatal |z| and MAF)
tst <- as.data.table(read_excel(file.path(D4, "braingvex_adult_bulk_comp.xlsx"), sheet = "tests", col_names = FALSE))
txt <- as.character(tst[[2]])[grepl("Reversed vs other classes, lead pair detected", as.character(tst[[1]]), fixed = TRUE)]
stopifnot(length(txt) == 1)
p_bulk <- as.numeric(sub(".*adjusted \\(postnatal \\|z\\|, MAF\\) OR [0-9.]+, P ([0-9.eE+-]+).*", "\\1", txt))
stopifnot(is.finite(p_bulk), p_bulk > 0, p_bulk < 1)
# exponent is a quoted plotmath string so that the minus is U+2212 (a bare U+2212 would not parse)
pe_e <- floor(log10(pexact)); pe_m <- round(pexact / 10^pe_e, 1)
if (pe_m >= 10) { pe_m <- 1; pe_e <- pe_e + 1 }
pann <- data.table(map = factor(levels(E$map), levels(E$map)),
                   lab = c(sprintf("italic(P) == %.3f", p_bulk), sprintf("GEE~italic(P) == %.1f %%*%% 10^\"%s\"", pe_m, negu(as.character(pe_e)))))
BR <- rbind(pann[, .(map, x = 1, xend = 2.5, y = 109, yend = 109)], pann[, .(map, x = 2, xend = 3, y = 105, yend = 105)],
            pann[, .(map, x = 2.5, xend = 2.5, y = 105, yend = 109)])
E[, `:=`(inside = k / n > 0.2, pct = 100 * k / n)]
pd <- ggplot(S, aes(cls, v)) +
  geom_col(aes(fill = fillkey, group = part), width = 0.72, colour = NA) +
  geom_text(data = E[inside == TRUE], aes(cls, pct, label = sprintf("%d/%d", k, n), colour = "white"), vjust = 1.5, size = T6) +
  scale_colour_identity() +
  geom_text(data = E[inside == FALSE], aes(cls, pct, label = sprintf("%d/%d", k, n)), vjust = -0.5, size = T6, colour = INK) +
  geom_segment(data = BR, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.25, colour = INK, inherit.aes = FALSE) +
  geom_text(data = pann, aes(x = 1.75, y = 111, label = lab), parse = TRUE, vjust = 0, size = T6, inherit.aes = FALSE) +
  facet_wrap(~map, nrow = 1) +
  scale_fill_manual(values = FILLS, guide = "none") +
  scale_y_continuous(limits = c(0, 122), breaks = c(0, 50, 100), expand = c(0, 0)) +
  scale_x_discrete(labels = c(Reversed = "Rev.", Amplified = "Amp.", Attenuated = "Att.")) +
  labs(x = NULL, y = "Detected in adult map (%)") +
  theme(strip.text = element_text(size = 6.5, hjust = 0.5, margin = margin(b = 1)), panel.spacing = unit(3, "mm"),
        axis.ticks.x = element_blank())

# --- assemble ---
top <- (tg(pa, "a") | (wrap_elements(full = pb) + labs(tag = "b") + theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1)))) + plot_layout(widths = c(36, 64))
bot <- (tg(pd, "c") | pc) + plot_layout(widths = c(84, 96))
fig <- (top / bot) + plot_layout(heights = c(64, 54))
ggsave(file.path(FIG_DIR, "fig4.png"), fig, width = 180, height = 120, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig4.pdf"), fig, width = 180, height = 120, units = "mm", device = cairo_pdf)
