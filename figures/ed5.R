# ED Fig 5: RNA-protein correlation before and after covariate removal (a), pi1 and neurite OR (b), protein beyond RNA vs dPDUI (c).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(splines); library(ragg)
})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
D2 <- file.path(ROOT, "data/02_traj"); D3 <- file.path(ROOT, "data/03_omics"); DF <- file.path(ROOT, "data/figures/ed5")

TAG_SIZE <- 12; AXIS_TITLE <- 8; AXIS_TEXT <- 7
theme_set(theme_nature(base_size = AXIS_TITLE, base_family = FONT) + theme(
  axis.text  = element_text(size = AXIS_TEXT, colour = INK),
  axis.title = element_text(size = AXIS_TITLE, colour = INK),
  axis.line  = element_line(linewidth = 0.3, colour = INK),
  axis.ticks = element_line(linewidth = 0.3, colour = INK),
  plot.margin = margin(2, 2, 2, 2, "mm")))
MINUS <- function(x) sub("-", "\u2212", x, fixed = TRUE)
tg <- function(p, t) p + labs(tag = t) + theme(
  plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FONT),
  plot.tag.position = c(0.005, 0.97))
LEG <- theme(legend.position = "top", legend.key.size = unit(2, "mm"),
             legend.margin = margin(0, 0, 0, 0), legend.text = element_text(size = 6),
             legend.box.spacing = unit(0.3, "mm"), legend.background = element_blank(),
             legend.key = element_blank())

# --- a: RNA-protein correlation per gene-protein pair, before and after age spline removal ---
man  <- fread(file.path(D3, "matched_sample_manifest_raw.tsv"))
rna  <- readRDS(file.path(D3, "rna_matched_raw.rds"))
prot <- readRDS(file.path(D3, "protein_matched_raw.rds"))
ids  <- man$sample_id
rna <- rna[, ids, drop = FALSE]; prot <- prot[, ids, drop = FALSE]
age <- man$age_years
Bs  <- ns(log10(age + 1), df = 2)
Xall <- cbind(1, as.matrix(Bs))
pp <- unique(fread(file.path(D3, "02_protein_traj_results.tsv"), select = c("gencode_gene_base", "protein_accession")))
rbase <- sub("[.].*$", "", rownames(rna))
pairs <- pp[gencode_gene_base %in% rbase & protein_accession %in% rownames(prot)]
pairs <- pairs[!duplicated(protein_accession)]
out <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  gi <- which(rbase == pairs$gencode_gene_base[i])[1]
  pi <- which(rownames(prot) == pairs$protein_accession[i])
  x <- as.numeric(rna[gi, ]); y <- as.numeric(prot[pi, ])
  ok <- is.finite(x) & is.finite(y)
  Xo <- Xall[ok, , drop = FALSE]; Ho <- solve(crossprod(Xo))
  rx <- x[ok] - Xo %*% (Ho %*% crossprod(Xo, x[ok]))
  ry <- y[ok] - Xo %*% (Ho %*% crossprod(Xo, y[ok]))
  out[[i]] <- c(cor(x[ok], y[ok]), cor(as.numeric(rx), as.numeric(ry)))
}
R <- do.call(rbind, out)
da <- rbind(data.table(state = "Before", r = R[, 1]),
            data.table(state = "After",  r = R[, 2]))
da <- da[is.finite(r)]
da[, state := factor(state, c("Before", "After"))]
meds <- da[, .(med = median(r)), by = state]
mb <- meds[state == "Before", med]; ma <- meds[state == "After", med]
thr <- 0.3
da[, seg := ifelse(r >= thr, "Coupled", "Uncoupled")]
fr <- da[, .(n = .N), by = .(state, seg)]
fr[, frac := n / sum(n), by = state]
fw <- dcast(fr, state ~ seg, value.var = "frac")
rects <- rbind(
  fw[, .(state, seg = "Uncoupled", frac = Uncoupled * 100, xmin = 0, xmax = Uncoupled * 100)],
  fw[, .(state, seg = "Coupled",   frac = Coupled * 100,   xmin = Uncoupled * 100, xmax = 100)]
)
rects[, lab_col := ifelse(seg == "Uncoupled", "white", INK)]

pa_top <- ggplot(da, aes(x = r, colour = state, fill = state)) +
  geom_vline(xintercept = thr, linewidth = 0.3, colour = NEUTRAL, linetype = "dotted") +
  geom_density(alpha = 0.18, linewidth = 0.5, colour = INK) +
  geom_vline(xintercept = mb, linewidth = 0.35, colour = MUTED, linetype = "22") +
  geom_vline(xintercept = ma, linewidth = 0.35, colour = C1_COLOR, linetype = "solid") +
  annotate("text", x = mb + 0.03, y = 4.65, label = sprintf("%.2f", mb),
           size = 6.5 / .pt, family = FONT, colour = MUTED, hjust = 0) +
  annotate("text", x = ma + 0.03, y = 5.25, label = sprintf("%.2f", ma),
           size = 6.5 / .pt, family = FONT, colour = C1_COLOR, fontface = "bold", hjust = 0) +
  coord_cartesian(xlim = c(-0.85, 1.05), ylim = c(0, 5.6), expand = FALSE) +
  scale_fill_manual(values = c(Before = NEUTRAL_LIGHT, After = C1_COLOR), guide = "none") +
  scale_colour_manual(values = c(Before = MUTED, After = C1_COLOR), guide = "none") +
  scale_x_continuous(breaks = seq(-0.75, 0.75, 0.25),
                     labels = function(b) MINUS(sprintf("%.2f", b))) +
  scale_y_continuous(breaks = seq(0, 4, 2)) +
  labs(x = "RNA\u2013protein correlation (r)", y = "Density")
pa_top <- pa_top + theme(axis.text.x = element_text(size = AXIS_TEXT),
                         plot.margin = margin(1, 1.5, 0.5, 2, "mm"))
pa_top <- tg(pa_top, "a")

pa_bot <- ggplot(rects, aes(xmin = xmin, xmax = xmax,
                            ymin = as.numeric(state) - 0.35,
                            ymax = as.numeric(state) + 0.35)) +
  geom_rect(aes(fill = seg), colour = "white", linewidth = 0.3) +
  geom_text(data = rects[frac >= 10],
            aes(x = (xmin + xmax) / 2, y = state, label = sprintf("%.0f%%", frac), colour = lab_col),
            size = 7 / .pt, family = FONT, fontface = "bold") +
  scale_colour_identity(guide = "none") +
  geom_text(data = rects[frac < 10],
            aes(x = xmax, y = state, label = sprintf("%.1f%%", frac)),
            size = 7 / .pt, family = FONT, colour = INK, fontface = "bold",
            hjust = -0.12) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0.55, 2.45), expand = FALSE, clip = "off") +
  scale_fill_manual(values = c(Uncoupled = MUTED, Coupled = SHARED_COLOR),
                    name = NULL,
                    labels = c("Coupled", "Uncoupled")) +
  scale_x_continuous(breaks = c(0, 50, 100), labels = c("0%", "50%", "100%"),
                     name = "Proportion") +
  scale_y_discrete(expand = expansion(add = 0.2)) +
  LEG + theme(legend.position = "top", legend.direction = "horizontal",
              legend.key.size = unit(1.8, "mm"), legend.text = element_text(size = 6),
              legend.margin = margin(t = 0, b = 1, unit = "mm"),
              legend.box.spacing = unit(0.2, "mm"),
              legend.key = element_blank(), legend.background = element_blank(),
              plot.margin = margin(0.5, 4, 1.5, 2, "mm"))
pa <- (pa_top / pa_bot) + plot_layout(heights = c(3.2, 1))

# --- b: pi1 of apaQTL effects (forest, left) and neurite-vs-soma OR (right) ---
s <- fread(file.path(D3, "S1_apa_prot_cis_vs_trans_371.tsv"))
mk <- function(cl, met) {
  r <- s[metric == met]; stopifnot(nrow(r) == 1)
  data.table(cluster = cl, obs = r$observed, null_mean = r$null_mean,
             null_q05 = r$null_q05, null_q95 = r$null_q95, emp = r$emp_p)
}
db <- rbind(mk("C1", "pi1_C1"), mk("C2", "pi1_C2"), mk("All", "pi1"))
db[, cluster := factor(cluster, levels = c("C2", "C1", "All"))]
db[, fold := obs / null_mean]
db[, ynum := as.numeric(cluster)]
xmax_b <- max(db$obs) + 0.04

pb_left <- ggplot(db, aes(x = obs, y = ynum)) +
  geom_segment(aes(x = null_q05, xend = null_q95, y = ynum, yend = ynum),
               colour = NEUTRAL, linewidth = 0.4) +
  geom_segment(aes(x = null_q05, xend = null_q05, y = ynum - 0.12, yend = ynum + 0.12),
               colour = NEUTRAL, linewidth = 0.4) +
  geom_segment(aes(x = null_q95, xend = null_q95, y = ynum - 0.12, yend = ynum + 0.12),
               colour = NEUTRAL, linewidth = 0.4) +
  geom_segment(aes(x = null_mean, xend = null_mean, y = ynum - 0.22, yend = ynum + 0.22),
               colour = MUTED, linewidth = 0.5) +
  geom_point(aes(fill = cluster), shape = 23, size = 2.3, stroke = 0.4, colour = INK) +
  geom_text(aes(label = sprintf("%.1f\u00d7", fold)), hjust = -0.5, nudge_x = 0.003,
            size = 6.5 / .pt, family = FONT, colour = INK) +
  scale_fill_manual(values = c(All = INK, C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_x_continuous(limits = c(0.075, xmax_b), breaks = seq(0.08, 0.22, 0.04),
                     name = expression(bold(pi[1]))) +
  scale_y_continuous(limits = c(0.05, 3.6), breaks = 1:3, labels = c("C2", "C1", "All"),
                     expand = expansion(mult = 0)) +
  labs(x = expression(bold(pi[1]))) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = AXIS_TEXT),
        plot.margin = margin(1.5, 1, 1.5, 2, "mm"))
pb_left <- tg(pb_left, "b")

dd <- fread(file.path(DF, "ed5d_adjusted_or.tsv"))
dd[, cls_all := factor(cls, levels = c("C2-specific", "C1-specific", "Shared"),
                       labels = c("C2", "C1", "All"))]
dd[, ynum := as.numeric(cls_all)]
dd[, set := factor(set, c("Ciolli", "Tushev"))]
dd[, sig := P < 0.05]
dd[, y_plot := ynum + ifelse(set == "Ciolli", -0.14, 0.14)]

pb_right <- ggplot(dd, aes(x = OR, y = y_plot)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = NEUTRAL, linetype = "dashed") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = cls_all), width = 0, linewidth = 0.3,
                na.rm = TRUE) +
  geom_point(aes(fill = cls_all, shape = set, alpha = sig), size = 1.9, stroke = 0.4,
             colour = INK) +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR, All = NEUTRAL), guide = "none") +
  scale_fill_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR, All = NEUTRAL), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_shape_manual(values = c(Ciolli = 21, Tushev = 24), name = NULL) +
  scale_x_log10(limits = c(0.4, 2.6), breaks = c(0.5, 1, 2), labels = c("0.5", "1", "2"),
                name = "Neurite OR") +
  scale_y_continuous(limits = c(0.05, 3.6), breaks = 1:3, labels = c("C2", "C1", "All"),
                     expand = expansion(mult = 0)) +
  LEG + theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              plot.margin = margin(1.5, 2, 0.5, 3.5, "mm"))
pb_right <- tg(pb_right, "c") + theme(plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FONT, hjust = 0))

pb_ind <- ggplot() +
  annotate("text", x = 0.46, y = 1, hjust = 0, vjust = 0.5,
           label = "\u2190 Soma", size = 6 / .pt, family = FONT, colour = MUTED) +
  annotate("text", x = 2.5, y = 1, hjust = 1, vjust = 0.5,
           label = "Neurite \u2192", size = 6 / .pt, family = FONT, colour = MUTED) +
  scale_x_log10(limits = c(0.4, 2.6)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 2, 0.5, 1, "mm"))
pb <- (pb_left | pb_right) / pb_ind + plot_layout(heights = c(1, 0.055))

# --- c: protein beyond RNA vs dPDUI, adjusted for mean RNA, mean protein and |dR| ---
d6 <- fread(file.path(D3, "S6_prot_beyond_rna_vs_dPDUI.tsv"))
rm_rna <- rowMeans(rna); names(rm_rna) <- sub("[.].*$", "", names(rm_rna))
rm_prot <- rowMeans(prot)
d6[, mean_rna := rm_rna[match(ens, names(rm_rna))]]
d6[, mean_prot := rm_prot[match(acc, names(rm_prot))]]
d6[, abs_dR := abs(dR)]
d6 <- d6[is.finite(mean_dA) & is.finite(beyond) & is.finite(mean_rna) & is.finite(mean_prot) & is.finite(abs_dR)]
fit <- lm(beyond ~ scale(mean_dA) + mean_rna + mean_prot + abs_dR, data = d6)
cf <- summary(fit)$coefficients
beta_sd <- cf[2, 1]; p_sd <- cf[2, 4]
covs <- lm(beyond ~ mean_rna + mean_prot + abs_dR, data = d6)$residuals
xres <- lm(mean_dA ~ mean_rna + mean_prot + abs_dR, data = d6)$residuals
dp <- data.table(sym = d6$sym, x = xres, y = covs)
dp[, decile := cut(x, breaks = quantile(x, seq(0, 1, 0.1)), include.lowest = TRUE)]
bins <- dp[, .(x_m = mean(x), y_m = mean(y), y_se = sd(y)/sqrt(.N)), by = decile]
setorder(bins, x_m)
CX <- range(bins$x_m) + c(-0.35, 0.35)
# y window from the decile summaries (mean +/- s.e.), not from individual genes
CY <- range(c(bins$y_m - bins$y_se, bins$y_m + bins$y_se)) * 1.25

pc_main <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "dotted") +
  geom_errorbar(data = bins, aes(x = x_m, ymin = y_m - y_se, ymax = y_m + y_se),
                width = 0, linewidth = 0.5, colour = C2_COLOR) +
  geom_point(data = bins, aes(x_m, y_m), size = 2.2, colour = C2_COLOR) +
  geom_smooth(data = bins, aes(x_m, y_m), method = "lm", formula = y ~ x,
              se = TRUE, colour = INK, linewidth = 0.55,
              fill = INK, alpha = 0.10) +
  annotate("text", x = CX[2], y = CY[2], vjust = 1,
           label = MINUS(sprintf("\u03b2 = %.2f, P = %.1f \u00d7 10\u207b\u2074 (n = %s)",
                           beta_sd, p_sd * 1e4, format(nrow(dp), big.mark = ","))),
           size = 6.5 / .pt, family = FONT, colour = INK, hjust = 1, lineheight = 0.9) +
  scale_x_continuous(breaks = seq(-2, 2, 1), labels = function(b) MINUS(sprintf("%.0f", b)),
                     name = "\u0394PDUI (residual)") +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2), labels = function(b) MINUS(sprintf("%.1f", b)),
                     name = "Protein beyond RNA (residual)") +
  coord_cartesian(xlim = CX, ylim = CY) +
  theme(plot.margin = margin(4, 1.5, 1.5, 2, "mm"),
        axis.title.y = element_text(size = AXIS_TITLE, margin = margin(r = 1.2, unit = "mm")),
        axis.title.x = element_text(size = AXIS_TITLE, margin = margin(t = 1.2, unit = "mm")),
        axis.text = element_text(size = AXIS_TEXT))
pc <- pc_main + labs(tag = "d") +
  theme(plot.tag = element_text(size = TAG_SIZE, face = "bold", family = FONT),
        plot.tag.position = c(0.005, 0.97))

# --- assemble ---
fig <- (pa | pb) / pc + plot_layout(widths = c(1, 1), heights = c(1, 0.8))
ggsave(file.path(FIG_DIR, "ed5.png"), fig, width = 180, height = 132,
       units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed5.pdf"), fig, width = 180, height = 132,
       units = "mm", device = cairo_pdf)
