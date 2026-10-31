# Fig 3 composite: apaQTL/eQTL gene-set UpSet (a), PAS enrichment (b), hexamer effects (c), CEP41 locus (c image), PP4, mediation, ACSF3.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(ggrastr); library(ragg); library(ggrepel) })
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
D <- file.path(ROOT, "data/figures/fig3")
D4 <- file.path(ROOT, "data/04_qtl")
theme_set(theme_nature(base_size = 7) + theme(
  text = element_text(family = FONT, colour = INK),
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 6.5, colour = INK),
  axis.line = element_line(linewidth = 0.35, colour = INK), axis.ticks = element_line(linewidth = 0.35, colour = INK),
  legend.text = element_text(size = 6), plot.margin = margin(1, 1.5, 1, 1, "mm")))
T6 <- 6 / .pt; T55 <- 5.5 / .pt; MINUS <- "\u2212"
neg <- function(x) sub("-", MINUS, x, fixed = TRUE)
tg <- function(p, tag) p + labs(tag = tag) + theme(plot.tag = element_text(size = 8, face = "bold", family = FONT, colour = INK))
PRE <- STAGE_COLORS[["Prenatal"]]; POST <- STAGE_COLORS[["Postnatal"]]; POSTTXT <- STAGE_TEXT[["Postnatal"]]
SC <- c(Prenatal = PRE, Postnatal = POST)
APA_C <- LAYER_COLORS[["APA"]]; EXP_C <- "#2E5B7E"

# --- gene sets (a) ---
apa <- fread(file.path(D, "fig3_apa_genes.tsv"))
eg <- fread(file.path(D4, "fig3_egenes.tsv"))
SETS <- list("Prenatal eQTL" = eg[stage == "prenatal", unique(ensg0)], "Prenatal apaQTL" = apa[stage == "prenatal", unique(ensg0)],
             "Postnatal apaQTL" = apa[stage == "postnatal", unique(ensg0)], "Postnatal eQTL" = eg[stage == "postnatal", unique(ensg0)])
stopifnot(length(intersect(SETS[[1]], SETS[[2]])) == 321, length(intersect(SETS[[3]], SETS[[4]])) == 1686)
ORD <- c("Prenatal apaQTL", "Postnatal apaQTL", "Prenatal eQTL", "Postnatal eQTL")
LAB <- c("Pre apaQTL", "Post apaQTL", "Pre eQTL", "Post eQTL")
U4 <- unique(unlist(SETS[ORD])); M4 <- sapply(SETS[ORD], function(x) U4 %in% x)
ix <- data.table(kk = apply(M4, 1, function(r) paste0(as.integer(r), collapse = "")))[, .N, by = kk][order(-N)]
ix[, col := seq_len(.N)]
ix[, apa := substr(kk, 1, 1) == "1" | substr(kk, 2, 2) == "1"][, eq := substr(kk, 3, 3) == "1" | substr(kk, 4, 4) == "1"]
ix[, grp := factor(fifelse(apa & !eq, "apaQTL only", fifelse(apa & eq, "apaQTL + eQTL", "eQTL only")), c("apaQTL only", "apaQTL + eQTL", "eQTL only"))]
GC <- c("apaQTL only" = APA_C, "apaQTL + eQTL" = "#6F8F86", "eQTL only" = EXP_C)
dots <- ix[, .(row = 1:4, on = as.integer(strsplit(kk, "")[[1]]) == 1L), by = .(col, grp)]
segs <- dots[on == TRUE, .(r0 = min(row), r1 = max(row)), by = .(col, grp)][r1 > r0]
rowlab <- sprintf("%s (%s)", LAB, trimws(format(lengths(SETS[ORD]), big.mark = ",")))
pa_bar <- ggplot(ix, aes(col, N, fill = grp)) + geom_col(width = 0.72) +
  geom_text(aes(label = format(N, big.mark = ",")), angle = 90, hjust = -0.15, size = 4.8 / .pt, colour = INK) +
  scale_fill_manual(values = GC, name = NULL) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1k", "10k"), expand = expansion(mult = c(0, 0.28))) +
  scale_x_continuous(limits = c(0.4, nrow(ix) + 0.6), expand = c(0, 0)) +
  labs(x = NULL, y = "Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
        legend.position = "inside", legend.position.inside = c(1, 1), legend.justification = c(1, 1), legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5.2), legend.background = element_blank(), plot.margin = margin(1, 1, 0, 1, "mm"))
pa_mat <- ggplot() +
  geom_tile(data = data.table(row = 1:4), aes(x = (nrow(ix) + 1) / 2, y = row), width = nrow(ix) + 0.2, height = 0.92,
            fill = rep(c("#F4F4F4", "white"), 2), inherit.aes = FALSE) +
  geom_point(data = dots[on == FALSE], aes(col, row), colour = "#D6D6D6", size = 1.3) +
  geom_segment(data = segs, aes(x = col, xend = col, y = r0, yend = r1, colour = grp), linewidth = 0.55) +
  geom_point(data = dots[on == TRUE], aes(col, row, colour = grp), size = 1.6) +
  scale_colour_manual(values = GC, guide = "none") +
  scale_y_reverse(breaks = 1:4, labels = rowlab, expand = c(0, 0.5)) +
  scale_x_continuous(limits = c(0.4, nrow(ix) + 0.6), expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 5.6, colour = c(APA_C, APA_C, EXP_C, EXP_C)), plot.margin = margin(0, 1, 1, 1, "mm"))

# --- b: PAS-distance enrichment of fine-mapping signal ---
pfp <- fread(file.path(D, "fig3_pip_profile.tsv"))[order(left)][, xi := seq_len(.N)]
brk <- c(pfp$left, tail(pfp$right, 1)); brk_lab <- neg(ifelse(abs(brk) >= 1000, paste0(brk / 1000, "k"), as.character(brk)))
nt2x <- function(v) approx(brk, seq(0.5, nrow(pfp) + 0.5, 1), xout = v)$y
el <- data.table(lab = c("UGUA", "AAUAAA", "GU-rich"), s = c(-100, -35, 10), e = c(-40, -10, 30))
el[, `:=`(xs = nt2x(s), xe = nt2x(e))]
ytop <- max(pfp$hi) * 1.10
pb1 <- ggplot(pfp, aes(xi)) +
  geom_rect(data = el, aes(xmin = xs, xmax = xe, ymin = -Inf, ymax = Inf), fill = "#EDEDED", inherit.aes = FALSE) +
  geom_segment(data = el, aes(x = xs, xend = xe, y = ytop, yend = ytop), linewidth = 0.5, colour = MUTED, inherit.aes = FALSE) +
  geom_text(data = el, aes(x = (xs + xe) / 2, y = ytop, label = lab), vjust = -0.5, size = 5.4 / .pt, colour = MUTED, inherit.aes = FALSE) +
  geom_hline(yintercept = 1, linetype = "22", linewidth = 0.3, colour = NEUTRAL) +
  geom_vline(xintercept = nt2x(0), linewidth = 0.3, colour = INK) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = APA_C, alpha = 0.28) +
  geom_line(aes(y = enrich), colour = APA_C, linewidth = 0.9) + geom_point(aes(y = enrich), colour = APA_C, size = 1.6) +
  annotate("text", x = nrow(pfp) + 0.4, y = max(pfp$hi) * 0.93, label = sprintf("%s credible sets", format(pfp$n_cs[1], big.mark = ",")), hjust = 1, size = T55, colour = MUTED) +
  scale_x_continuous(breaks = seq(0.5, nrow(pfp) + 0.5, 1), labels = brk_lab, expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(0, max(pfp$hi) * 1.24), breaks = seq(0, 8, 2)) +
  labs(x = "Distance to cleavage site (nt)", y = "Fine-mapping enrichment") +
  theme(axis.text.x = element_text(size = 5.4))

# --- c: hexamer-changing pairs, predicted vs observed direction ---
hx <- fread(file.path(D, "fig3_hexamer_screen.tsv"))[pred != 0]
hx[, stage := factor(tools::toTitleCase(stage), c("Prenatal", "Postnatal"))]
hx[, grp := factor(fifelse(pred > 0, "Longer", "Shorter"), c("Longer", "Shorter"))]
hx[, conc := sign(slope) == sign(pred)]
stopifnot(nrow(hx) == 28, sum(hx$conc) == 21)
gn <- hx[, .(n = .N, k = sum(conc)), by = grp]
set.seed(4); hx[, xj := as.integer(grp) + runif(.N, -0.22, 0.22)]
pb2 <- ggplot(hx, aes(xj, slope)) +
  annotate("rect", xmin = 0.55, xmax = 1.45, ymin = 0, ymax = Inf, fill = APA_C, alpha = 0.08) +
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = 0, fill = APA_C, alpha = 0.08) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = INK) +
  geom_point(data = hx[conc == TRUE], aes(fill = stage), shape = 21, colour = "white", stroke = 0.3, size = 2.6) +
  geom_point(data = hx[conc == FALSE], aes(colour = stage), shape = 21, fill = "white", stroke = 0.6, size = 2.3) +
  geom_text(data = gn, aes(x = as.integer(grp), y = c(Longer = 1.38, Shorter = -1.38)[as.character(grp)], label = sprintf("%d/%d", k, n)), size = T6, inherit.aes = FALSE) +
  annotate("text", x = 2.45, y = 1.62, label = sprintf("%d/%d, P = %.3f", sum(hx$conc), nrow(hx), binom.test(sum(hx$conc), nrow(hx))$p.value), hjust = 1, size = 5.2 / .pt, colour = MUTED) +
  geom_text_repel(data = hx[gene == "CEP41" & stage == "Prenatal"], aes(xj, slope, label = "CEP41"), fontface = "bold.italic", size = T55,
                  nudge_x = 0.30, nudge_y = 0.35, segment.size = 0.3, segment.colour = INK, min.segment.length = 0, seed = 2, inherit.aes = FALSE) +
  scale_fill_manual(values = SC, name = NULL) + scale_colour_manual(values = SC, guide = "none") +
  scale_x_continuous(breaks = 1:2, labels = c("Longer", "Shorter"), limits = c(0.5, 2.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.5, 1.7), breaks = c(-1, 0, 1), labels = neg) +
  labs(x = "Predicted 3\u2032UTR", y = "Observed PDUI effect") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size = 5.6, lineheight = 0.9),
        legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0), legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 5.4), legend.background = element_blank())

# --- CEP41 locus image (from fig3_panel_c.py) ---
img_c <- png::readPNG(file.path(FIG_DIR, "figure_3d.png"))
pc1 <- ggplot() + annotation_raster(img_c, 0, 1, 0, 1) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void()

# --- d: PP4 distribution by stage ---
pl <- fread(file.path(D, "fig3d_coloc_pp_by_gene.tsv")); br <- seq(0, 1, 0.05)
hd <- pl[, .(n = .N), by = .(stage, bin = cut(PP_H4, br, include.lowest = TRUE, labels = FALSE))]
hd <- hd[CJ(stage = c("Prenatal", "Postnatal"), bin = 1:20), on = .(stage, bin)][is.na(n), n := 0]
hd[, frac := n / sum(n) * 100, by = stage][, mid := (bin - 0.5) * 0.05][, yv := fifelse(stage == "Prenatal", frac, -frac)][, fillk := fifelse(mid > 0.8, "coloc", stage)]
sh <- pl[, .(k = sum(PP_H4 > 0.8), n = .N), by = stage]; stopifnot(sh[stage == "Prenatal", k] == 164, sh[stage == "Postnatal", k] == 608)
ymx <- max(abs(hd$yv)) * 1.12
pd <- ggplot(hd, aes(mid, yv, fill = fillk)) + geom_col(width = 0.05, colour = "white", linewidth = 0.12) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_vline(xintercept = 0.8, linetype = "22", linewidth = 0.3) +
  scale_fill_manual(values = c(Prenatal = PRE, Postnatal = POST, coloc = EXP_C), guide = "none") +
  annotate("text", x = 0.78, y = ymx * 0.80, label = sprintf("%.1f%%\n(%s/%s)", 100 * sh[stage == "Prenatal", k] / sh[stage == "Prenatal", n], format(sh[stage == "Prenatal", k], big.mark = ","), format(sh[stage == "Prenatal", n], big.mark = ",")), hjust = 1, size = T55, colour = EXP_C, lineheight = 0.9) +
  annotate("text", x = 0.78, y = -ymx * 0.80, label = sprintf("%.1f%%\n(%s/%s)", 100 * sh[stage == "Postnatal", k] / sh[stage == "Postnatal", n], format(sh[stage == "Postnatal", k], big.mark = ","), format(sh[stage == "Postnatal", n], big.mark = ",")), hjust = 1, size = T55, colour = EXP_C, lineheight = 0.9) +
  annotate("text", x = 0.02, y = ymx * 0.95, label = "Prenatal", hjust = 0, size = T6, colour = PRE) +
  annotate("text", x = 0.02, y = -ymx * 0.95, label = "Postnatal", hjust = 0, size = T6, colour = POSTTXT) +
  scale_x_continuous(breaks = c(0, 0.5, 0.8, 1), labels = c("0", "0.5", "0.8", "1"), expand = expansion(add = 0.01)) +
  scale_y_continuous(labels = function(x) abs(x), limits = c(-ymx, ymx)) + labs(x = "Colocalization with eQTL (PP4)", y = "apaQTL genes (%)")

# --- e: mediation in both directions ---
md <- fread(file.path(D4, "04_mediation_results_module06.tsv")); stopifnot(nrow(md) == 772)
md[, cls := fcase(mediation_call %in% c("APA -> EXP", "Partial mediation"), "APA-mediated", mediation_call == "EXP -> APA", "Expression-mediated", default = "None")]
md[, `:=`(fwd = ab_fwd / c_total, rev = (a - d_prime) / a)]
nm <- md[, .N, by = cls]; stopifnot(nm[cls == "APA-mediated", N] == 162)
pe <- ggplot(md, aes(fwd, rev)) + geom_hline(yintercept = 0, linewidth = 0.25, colour = NEUTRAL) + geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL) +
  rasterise(geom_point(data = md[cls == "None"], colour = NEUTRAL_LIGHT, size = 0.6, stroke = 0), dpi = 600) +
  rasterise(geom_point(data = md[cls == "APA-mediated"], colour = APA_C, shape = 16, size = 1.4, alpha = 0.85), dpi = 600) +
  rasterise(geom_point(data = md[cls == "Expression-mediated"], colour = EXP_C, shape = 17, size = 1.7, alpha = 0.9), dpi = 600) +
  geom_point(data = md[gene_name == "ACSF3" & stage == "prenatal"], shape = 21, fill = NA, colour = INK, size = 2.2, stroke = 0.5) +
  annotate("text", x = md[gene_name == "ACSF3" & stage == "prenatal", fwd] + 0.06, y = md[gene_name == "ACSF3" & stage == "prenatal", rev] - 0.06, label = "ACSF3", fontface = "italic", size = T6, hjust = 0) +
  annotate("text", x = 0.98, y = -0.22, label = sprintf("APA-mediated %d", nm[cls == "APA-mediated", N]), hjust = 1, size = T6, colour = APA_C) +
  annotate("text", x = -0.28, y = 1.12, label = sprintf("Expression-mediated %d", nm[cls == "Expression-mediated", N]), hjust = 0, size = T6, colour = EXP_C) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = neg(c("-1", "-0.5", "0", "0.5", "1"))) + scale_y_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = neg(c("-1", "-0.5", "0", "0.5", "1"))) +
  coord_cartesian(xlim = c(-1.2, 1), ylim = c(-1.2, 1), clip = "off") +
  labs(x = "SNP \u2192 expression effect via APA", y = "SNP \u2192 APA effect via expression")

# --- f: ACSF3 mediation path (both stages) and split of the total SNP -> expression effect ---
mf <- md[gene_name == "ACSF3"][, stage := factor(tools::toTitleCase(stage), c("Prenatal", "Postnatal"))][order(stage)]
stopifnot(nrow(mf) == 2, all(abs(mf$ab_fwd + mf$c_prime - mf$c_total) < 1e-6))
f2 <- function(v) sprintf("%.2f", v)
cl <- function(v) sprintf("<span style='color:%s'>%s</span> / <span style='color:%s'>%s</span>", STAGE_TEXT[["Prenatal"]], f2(v[1]), STAGE_TEXT[["Postnatal"]], f2(v[2]))
nodes <- data.table(x = c(0.12, 0.50, 0.88), y = c(0.30, 0.82, 0.30), lab = c("SNP", "PDUI", "Expression"), col = c(INK, APA_C, EXP_C))
arr <- data.table(x = c(0.18, 0.60, 0.30), y = c(0.40, 0.72, 0.24), xe = c(0.40, 0.82, 0.72), ye = c(0.72, 0.40, 0.24))
alab <- data.table(x = c(0.27, 0.73, 0.50), y = c(0.60, 0.60, 0.12), hj = c(1, 0, 0.5),
                   lab = c(paste0("*a* = ", cl(mf$a)), paste0("*b* = ", cl(mf$b)), paste0("*c*\u2032 = ", cl(mf$c_prime))))
pf1 <- ggplot() +
  geom_segment(data = arr, aes(x, y, xend = xe, yend = ye), arrow = arrow(length = unit(1.2, "mm"), type = "closed"), linewidth = 0.65, colour = c(APA_C, APA_C, NEUTRAL)) +
  geom_label(data = nodes, aes(x, y, label = lab), size = 6 / .pt, lineheight = 0.85, linewidth = 0.3, label.r = unit(1, "mm"), colour = nodes$col, fill = "white") +
  ggtext::geom_richtext(data = alab, aes(x, y, label = lab, hjust = hj), size = 5.2 / .pt, fill = NA, label.colour = NA, label.padding = unit(0, "mm")) +
  coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(0, 0.95), expand = FALSE, clip = "off") + theme_void() +
  labs(title = "ACSF3") + theme(plot.title = element_text(size = 6.5, face = "bold.italic", colour = INK, hjust = 0, margin = margin(0, 0, 0.5, 0, "mm")))
dec <- melt(mf[, .(stage, Indirect = ab_fwd, Direct = c_prime, pm = prop_med)], id.vars = c("stage", "pm"), variable.name = "part", value.name = "v")
dec[, part := factor(part, c("Direct", "Indirect"))]
pf2 <- ggplot(dec, aes(v, stage, fill = part)) + geom_col(width = 0.62, colour = "white", linewidth = 0.2) +
  geom_text(data = dec[part == "Indirect"], aes(x = v / 2, label = sprintf("%.0f%%", 100 * pm)), colour = "white", size = 5.6 / .pt, fontface = "bold") +
  scale_fill_manual(values = c(Indirect = APA_C, Direct = "#C9CED4"), labels = c(Indirect = "via APA (a\u00D7b)", Direct = "direct (c\u2032)"), name = NULL, breaks = c("Indirect", "Direct")) +
  scale_y_discrete(limits = rev) + scale_x_continuous(expand = c(0, 0), limits = c(0, 0.85), breaks = c(0, 0.4, 0.8)) +
  labs(x = "SNP \u2192 expression total effect (c)", y = NULL) +
  theme(axis.text.y = element_text(colour = unname(STAGE_TEXT[c("Postnatal", "Prenatal")]), size = 6), axis.ticks.y = element_blank(),
        legend.position = "top", legend.key.size = unit(2, "mm"), legend.text = element_text(size = 5.4), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -2, 0, "mm"))

# --- e (PDUI and expression boxes) and CEP41 box helpers ---
mix <- function(col, f) { m <- grDevices::col2rgb(col) / 255; grDevices::rgb(t(m * (1 - f) + f)) }
SC_LIGHT <- c(Prenatal = mix(PRE, 0.45), Postnatal = mix(POST, 0.45))
CGT <- c("C/C", "C/T", "T/T")
bx6 <- fread(file.path(D, "fig3_locus_boxstats.tsv"))
EXF <- file.path(D4, "fig3_expr_boxstats.tsv")
box2 <- function(d, levs, ylab, fillv, note, show_n, show_x) {
  d <- copy(d); d[, stage := factor(stage, c("Prenatal", "Postnatal"))][, gt := factor(genotype, levs)]
  d[, xpos := as.integer(gt) + fifelse(stage == "Prenatal", -0.19, 0.19)]
  ol <- d[n >= 5 & nzchar(as.character(outliers)), .(v = as.numeric(strsplit(as.character(outliers), ",")[[1]])), by = .(stage, gt, xpos)]
  rng <- range(c(d$whisker_lo, d$whisker_hi, ol$v), na.rm = TRUE); pad <- diff(rng)
  bands <- data.table(x = seq_along(levs), f = rep(c("#F3F3F3", "white"), length.out = length(levs)))
  p <- ggplot(d) +
    geom_rect(data = bands, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf, fill = I(f)), inherit.aes = FALSE) +
    geom_segment(aes(x = xpos, xend = xpos, y = whisker_lo, yend = whisker_hi, colour = stage), linewidth = 0.35) +
    geom_rect(aes(xmin = xpos - 0.17, xmax = xpos + 0.17, ymin = q1, ymax = q3, fill = I(fillv[as.character(stage)]), colour = stage), linewidth = 0.35) +
    geom_segment(aes(x = xpos - 0.17, xend = xpos + 0.17, y = median, yend = median), colour = "white", linewidth = 0.6) +
    geom_point(data = ol, aes(xpos, v, colour = stage), size = 0.35, stroke = 0) +
    geom_line(aes(xpos, mean, colour = stage, group = stage), linewidth = 0.5) +
    geom_point(aes(xpos, mean, colour = stage), shape = 21, fill = "white", size = 1.1, stroke = 0.45) +
    annotate("text", x = 0.55, y = rng[2] + pad * 0.16, label = note, hjust = 0, vjust = 1, size = T55, parse = TRUE) +
    scale_colour_manual(values = c(Prenatal = PRE, Postnatal = POSTTXT), guide = "none") +
    scale_x_continuous(breaks = seq_along(levs), labels = if (show_x) levs else NULL, expand = c(0, 0), limits = c(0.5, length(levs) + 0.5)) +
    scale_y_continuous(labels = neg) + labs(x = NULL, y = ylab) +
    coord_cartesian(ylim = c(rng[1] - pad * (if (show_n) 0.30 else 0.04), rng[2] + pad * 0.18), clip = "off") +
    theme(axis.ticks.x = element_blank())
  if (show_n) p <- p + geom_text(aes(x = as.integer(gt), y = rng[1] - pad * fifelse(stage == "Prenatal", 0.11, 0.23), label = n, colour = stage), size = 5.2 / .pt)
  p
}
ePD <- bx6[gene == "CEP41"]; eEX <- fread(EXF)[gene == "CEP41"]; eEX[, genotype := CGT[dosage_alt + 1]]
bnote <- function(b) sprintf("list(beta[Pre]==%.2f, ~~beta[Post]==%.2f)", b[1], b[2])
mb <- fread(file.path(D4, "05_mash_posterior_apa_handoff.tsv"))[phenotype_id == "ENST00000223208.10|CEP41|chr7:130393771-130399039|-" & variant_id == "chr7:130393789:C:T"]
stopifnot(nrow(mb) == 1)
pbx1 <- box2(ePD, CGT, "Adjusted PDUI", SC, bnote(c(mb$b_prenatal, mb$b_postnatal)), FALSE, FALSE)
pbx2 <- box2(eEX, CGT, "Adjusted expression", SC_LIGHT, bnote(eEX[, .(b = slope[1]), by = stage][order(-(stage == "Prenatal")), b]), TRUE, TRUE)

# --- assemble ---
r1 <- (wrap_elements(full = tg(pa_bar, "a") / pa_mat + plot_layout(heights = c(1.25, 1))) | tg(pb1, "b") | tg(pb2, "c")) + plot_layout(widths = c(62, 68, 50))
r2 <- (tg(pc1, "d") | wrap_elements(full = (tg(pbx1, "e") / pbx2) + plot_layout(heights = c(1, 1.15)))) + plot_layout(widths = c(124, 56))
r3 <- (tg(pd, "f") | tg(pe, "g") | wrap_elements(full = tg(pf1, "h") / pf2 + plot_layout(heights = c(1.5, 1)))) + plot_layout(widths = c(50, 62, 68))
fig <- wrap_elements(full = r1) / wrap_elements(full = r2) / wrap_elements(full = r3) + plot_layout(heights = c(54, 66, 52))
ggsave(file.path(FIG_DIR, "fig3.png"), fig, width = 180, height = 172, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig3.pdf"), fig, width = 180, height = 172, units = "mm", device = cairo_pdf)
