# ED Fig 1: sample flow (a), donor age density with supported windows (b), PolyASite 2.0 distance to proximal sites (c).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(readxl); library(ragg); library(grid)
})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
DF <- file.path(ROOT, "data/figures/ed1"); D2 <- file.path(ROOT, "data/02_traj")
th <- theme_nature(base_size = 6) + theme(
  axis.text   = element_text(size = 6, colour = INK),
  axis.title  = element_text(size = 7, colour = INK),
  axis.line   = element_line(linewidth = 0.25, colour = INK),
  axis.ticks  = element_line(linewidth = 0.25, colour = INK),
  axis.ticks.length = unit(0.8, "mm"),
  legend.text  = element_text(size = 5.5, colour = INK),
  legend.title = element_text(size = 5.5, colour = INK),
  strip.text   = element_text(size = 6, colour = INK),
  plot.margin  = margin(1, 1.5, 1, 1.5, "mm"))
theme_set(th)
tg <- function(p, t) p + labs(tag = t) + theme(
  plot.tag = element_text(size = 8, face = "bold", family = FONT),
  plot.tag.position = c(0, 1))
fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

# --- flow counts ---
flow <- fread(file.path(DF, "ed1a_flow.tsv"))
FN <- setNames(flow$Total, flow$step)
n_pfc  <- FN[["PFC samples (raw metadata)"]]
n_rin  <- FN[["RIN >= 5 samples"]]
n_apa  <- FN[["APA profiles after QC"]]
n_geno <- FN[["Genotyped (QTL)"]]
n_expr <- FN[["Expression samples"]]
n_prot <- FN[["Matched proteome donors"]]
n_rin_fail <- n_pfc - n_rin
n_qc_fail  <- n_rin - n_apa

# --- a: CONSORT-style branching flow diagram ---
gl <- gList()
bw <- 0.54; xcen <- 0.38
xfilt <- 0.80
top <- 0.985
main <- list(
  list(lab = sprintf("CHB cohort\nPFC bulk RNA-seq\nn = %s", fmt_n(n_pfc)), fill = "#F0F0F0", nl = 3),
  list(lab = sprintf("RIN \u2265 5\nn = %s", fmt_n(n_rin)),                  fill = "#F0F0F0", nl = 2),
  list(lab = sprintf("APA profiles\n(DaPars2 QC)\nn = %s", fmt_n(n_apa)),   fill = "#E8F4E8", nl = 3))
filters <- c(NA, sprintf("RIN < 5 (n = %s)", fmt_n(n_rin_fail)), sprintf("QC fail (n = %s)", fmt_n(n_qc_fail)))
bh_for  <- function(nl) if (nl >= 3) 0.130 else 0.095
gap_main <- 0.150
yc_main <- numeric(3); ycur <- top
for (i in 1:3) {
  bh <- bh_for(main[[i]]$nl)
  yc_main[i] <- ycur - bh / 2
  ycur <- ycur - bh - gap_main
}
for (i in 1:3) {
  bh <- bh_for(main[[i]]$nl)
  gl <- gList(gl,
    roundrectGrob(x = xcen, y = yc_main[i], width = bw, height = bh,
                  r = unit(1.2, "mm"),
                  gp = gpar(fill = main[[i]]$fill, col = INK, lwd = 0.5)))
  gl <- gList(gl,
    textGrob(main[[i]]$lab, x = xcen, y = yc_main[i],
             gp = gpar(fontsize = 5, fontfamily = FONT, col = INK, lineheight = 0.85)))
  if (i > 1) {
    y_prev_bh <- bh_for(main[[i-1]]$nl)
    y_from <- yc_main[i-1] - y_prev_bh / 2
    y_to   <- yc_main[i] + bh / 2
    gl <- gList(gl,
      segmentsGrob(x0 = xcen, y0 = y_from, x1 = xcen, y1 = y_to + 0.004,
                   arrow = arrow(length = unit(1, "mm"), type = "closed"),
                   gp = gpar(col = INK, lwd = 0.5, fill = INK)))
    if (!is.na(filters[i])) {
      y_mid <- (y_from + y_to) / 2
      gl <- gList(gl,
        segmentsGrob(x0 = xcen + bw/2 + 0.01, y0 = y_mid,
                     x1 = xfilt - 0.10, y1 = y_mid,
                     gp = gpar(col = NEUTRAL, lwd = 0.4, lty = "dotted")))
      gl <- gList(gl,
        textGrob(filters[i], x = xfilt, y = y_mid, hjust = 0.5,
                 gp = gpar(fontsize = 4.5, fontfamily = FONT, col = MUTED, fontface = "italic")))
    }
  }
}
branches <- list(
  list(lab = sprintf("Genotyped\nn = %s", fmt_n(n_geno)),  x = 0.13),
  list(lab = sprintf("Expression\nn = %s", fmt_n(n_expr)), x = 0.38),
  list(lab = sprintf("Proteome\nn = %s", fmt_n(n_prot)),   x = 0.63))
bh3 <- bh_for(3)
y_branch_top <- yc_main[3] - bh3 / 2
bh2 <- 0.095; bw2 <- 0.22
gap_branch <- 0.165
y_branch <- y_branch_top - gap_branch - bh2 / 2
for (j in seq_along(branches)) {
  bx <- branches[[j]]$x
  gl <- gList(gl,
    segmentsGrob(x0 = xcen, y0 = y_branch_top,
                 x1 = bx,   y1 = y_branch + bh2/2 + 0.004,
                 arrow = arrow(length = unit(0.8, "mm"), type = "closed"),
                 gp = gpar(col = INK, lwd = 0.4, fill = INK)))
  gl <- gList(gl,
    roundrectGrob(x = bx, y = y_branch, width = bw2, height = bh2,
                  r = unit(1.2, "mm"),
                  gp = gpar(fill = "#E8F4E8", col = INK, lwd = 0.5)))
  gl <- gList(gl,
    textGrob(branches[[j]]$lab, x = bx, y = y_branch,
             gp = gpar(fontsize = 5, fontfamily = FONT, col = INK, lineheight = 0.9)))
}
pa <- wrap_elements(full = gTree(children = gl))

# --- b: donor age density (log10(age + 1)) with KDE-supported windows ---
raw <- fread(file.path(ROOT, "raw/all.912.meta"))
apa_samples <- names(fread(file.path(ROOT, "data/04_qtl/apa_after_combat.tsv"), nrows = 0))[-1]
donor_ids   <- names(fread(file.path(ROOT, "data/01_quant/apa_combat_donor_mean.tsv"), nrows = 0))[-1]
stopifnot(length(apa_samples) == n_apa, length(donor_ids) == 750, all(donor_ids %in% apa_samples),
          sum(raw$RIN < 5) == n_rin_fail, sum(raw$RIN >= 5) == n_rin)
raw <- raw[sample %in% donor_ids]
stopifnot(nrow(raw) == length(donor_ids))
N_DON <- nrow(raw)

sup <- fread(file.path(D2, "01_support.tsv"))
sup <- sup[order(grid_index)]
sup[, run := cumsum(c(1, diff(as.integer(is_supported)) != 0))]
segs <- sup[is_supported == TRUE, .(lo = min(grid_age), hi = max(grid_age)), by = run][order(lo)]
A_LO <- segs$lo[1]; A_HI <- segs$hi[1]; B_LO <- segs$lo[2]; B_HI <- segs$hi[2]

la <- function(a) log10(a + 1)
raw[, x := la(as.numeric(age))]
dd <- density(raw$x, adjust = 0.6, from = la(-0.55), to = la(110), n = 600)
D <- data.table(x = dd$x, y = dd$y)
YMAX <- max(D$y) * 1.15
gap_lo <- la(A_HI); gap_hi <- la(B_LO)
win_pre_lo <- la(A_LO); win_pre_hi <- la(A_HI)
win_post_lo <- la(B_LO); win_post_hi <- la(B_HI)

pb <- ggplot() +
  annotate("rect", xmin = win_pre_lo,  xmax = win_pre_hi,
           ymin = 0, ymax = YMAX, fill = STAGE_COLORS[["Prenatal"]], alpha = 0.15) +
  annotate("rect", xmin = win_post_lo, xmax = win_post_hi,
           ymin = 0, ymax = YMAX, fill = STAGE_COLORS[["Postnatal"]], alpha = 0.15) +
  annotate("rect", xmin = gap_lo, xmax = gap_hi,
           ymin = 0, ymax = YMAX, fill = NEUTRAL_LIGHT, alpha = 0.45) +
  annotate("text", x = (win_pre_lo + win_pre_hi) / 2, y = YMAX * 0.97,
           label = "Prenatal", size = 4 / .pt, colour = STAGE_TEXT[["Prenatal"]],
           family = FONT, fontface = "italic", hjust = 0.5, vjust = 1) +
  annotate("text", x = (win_post_lo + win_post_hi) / 2, y = YMAX * 0.97,
           label = "Postnatal", size = 4 / .pt, colour = STAGE_TEXT[["Postnatal"]],
           family = FONT, fontface = "italic", hjust = 0.5, vjust = 1) +
  annotate("text", x = (gap_lo + gap_hi) / 2, y = YMAX * 0.55,
           label = "Low\ndensity", size = 4 / .pt, colour = MUTED,
           family = FONT, fontface = "italic", hjust = 0.5, vjust = 0.5, lineheight = 0.85) +
  geom_line(data = D, aes(x, y), colour = INK, linewidth = 0.4) +
  geom_area(data = D, aes(x, y), fill = NEUTRAL, alpha = 0.12) +
  geom_rug(data = raw, aes(x), colour = INK, alpha = 0.25,
           linewidth = 0.2, length = unit(0.8, "mm")) +
  geom_vline(xintercept = c(win_pre_lo, win_pre_hi), colour = STAGE_TEXT[["Prenatal"]],
             linewidth = 0.25, linetype = "21") +
  geom_vline(xintercept = c(win_post_lo, win_post_hi), colour = STAGE_TEXT[["Postnatal"]],
             linewidth = 0.25, linetype = "21") +
  scale_x_continuous(
    breaks = la(c(0, 1, 10, 100)),
    labels = c("Birth", "1 y", "10 y", "100 y")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Age (log scale)", y = sprintf("Donor density (n = %d)", N_DON)) +
  theme(plot.margin = margin(1, 2, 2, 1.5, "mm"))

# --- c: PolyASite 2.0 distance, DaPars2 proximal vs random 3' UTR site ---
pe <- as.data.table(read_excel(file.path(ROOT, "ref/polyasite2/polyasite_validation.xlsx"), sheet = "per_event"))
pe[, d_prox := as.numeric(d_prox)]; pe[, d_rand := as.numeric(d_rand)]
pe[, is_apaQTL_lead := as.logical(is_apaQTL_lead)]
ev <- pe[!is.na(d_prox) & !is.na(d_rand)]
L2 <- ev[is_apaQTL_lead == TRUE]
kA <- sum(ev$d_prox < ev$d_rand); kL <- sum(L2$d_prox < L2$d_rand)
d <- rbind(
  data.table(set = "All events",         type = "DaPars2 proximal",    dist = ev$d_prox),
  data.table(set = "All events",         type = "Random 3\u2032 UTR",    dist = ev$d_rand),
  data.table(set = "apaQTL lead events", type = "DaPars2 proximal",    dist = L2$d_prox),
  data.table(set = "apaQTL lead events", type = "Random 3\u2032 UTR",    dist = L2$d_rand))
d[, set := factor(set, c("All events", "apaQTL lead events"))]
d[, logdist := log10(dist + 1)]
ann <- data.table(
  set = factor(c("All events", "apaQTL lead events"), levels = levels(d$set)),
  lab = c(sprintf("All events\nn = %s, %.1f%% closer", fmt_n(nrow(ev)), 100 * kA / nrow(ev)),
          sprintf("apaQTL lead\nn = %s, %.1f%% closer", fmt_n(nrow(L2)), 100 * kL / nrow(L2))))
pc <- ggplot(d, aes(logdist, fill = type, colour = type)) +
  geom_density(alpha = 0.18, linewidth = 0.4, adjust = 1.0) +
  facet_wrap(~set, nrow = 1, labeller = labeller(set = setNames(ann$lab, ann$set))) +
  scale_fill_manual(values = c("DaPars2 proximal" = LAYER_COLORS[["APA"]],
                               "Random 3\u2032 UTR" = NEUTRAL_LIGHT), name = NULL) +
  scale_colour_manual(values = c("DaPars2 proximal" = LAYER_COLORS[["APA"]],
                                 "Random 3\u2032 UTR" = NEUTRAL), name = NULL) +
  scale_x_continuous(
    breaks = log10(c(1, 11, 101, 1001)),
    labels = c("0", "10", "100", "1,000")) +
  labs(x = "Distance to PolyASite 2.0 cluster (nt)", y = "Density") +
  theme(plot.margin = margin(3, 2, 2, 1.5, "mm"),
        strip.text = element_text(size = 6, colour = INK, lineheight = 0.9,
                                  margin = margin(t = 1.5, b = 2, unit = "pt")),
        legend.position = "bottom", legend.justification = "center",
        legend.direction = "horizontal",
        legend.key.size = unit(2, "mm"), legend.key.height = unit(2, "mm"),
        legend.key.width = unit(3, "mm"),
        legend.spacing.x = unit(1, "mm"), legend.margin = margin(0, 0, 0, 0))

fig <- (tg(pa, "a") | tg(pb, "b") | tg(pc, "c")) +
  plot_layout(widths = c(46, 56, 78))
ggsave(file.path(FIG_DIR, "ed1.png"), fig,
       width = 180, height = 75, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed1.pdf"), fig,
       width = 180, height = 75, units = "mm", device = cairo_pdf)
