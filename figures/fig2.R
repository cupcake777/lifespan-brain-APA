# Fig 2: APA-expression coupling (a), QQ enrichment (b), motif and RBP enrichment (c), gene model (d), layer trajectories (e).
# Reads data/02_traj, data/03_omics, data/figures/fig2 and figures/fig2_d.png; writes figures/fig2.png and fig2.pdf.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(MASS)
  library(ggrastr)
})

source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
D2 <- file.path(ROOT, "data/02_traj"); D3 <- file.path(ROOT, "data/03_omics"); DF <- file.path(ROOT, "data/figures/fig2")
tg <- function(p, tag) p + labs(tag = tag) + theme(plot.tag = element_text(size = 8, face = "bold", family = FONT, colour = INK))
MINUS <- "\u2212"
fmt2 <- function(x) sub("-", MINUS, sprintf("%.2f", x), fixed = TRUE)
fmt_minus <- function(x) sub("-", MINUS, trimws(format(x)), fixed = TRUE)

# --- helpers ---
kde_gemini <- function(d, grp_label, n = 120) {
  span_x <- diff(lim_x); span_y <- diff(lim_y)
  ext_x <- c(lim_x[1] - 0.15 * span_x, lim_x[2] + 0.15 * span_x)
  ext_y <- c(lim_y[1] - 0.15 * span_y, lim_y[2] + 0.15 * span_y)

  ks <- MASS::kde2d(d$expr_dev, d$pdui_dev, n = n, lims = c(ext_x, ext_y))
  z_max <- max(ks$z)

  fill_levs <- z_max * c(0.35, 0.60, 0.82)
  cl_fill <- grDevices::contourLines(ks$x, ks$y, ks$z, levels = fill_levs)
  df_fill <- rbindlist(lapply(seq_along(cl_fill), function(i) {
    l <- cl_fill[[i]]
    data.table(x = l$x, y = l$y, level = l$level, id = paste0(grp_label, "_f_", i), grp = grp_label)
  }))

  line_levs <- z_max * c(0.18, 0.38, 0.62, 0.85)
  cl_line <- grDevices::contourLines(ks$x, ks$y, ks$z, levels = line_levs)
  df_line <- rbindlist(lapply(seq_along(cl_line), function(i) {
    l <- cl_line[[i]]
    data.table(x = l$x, y = l$y, level = l$level, id = paste0(grp_label, "_l_", i), grp = grp_label)
  }))

  list(fill = df_fill, line = df_line)
}

pv <- function(m) s1p[metric == m, observed]; nv <- function(m) s1p[metric == m, null_mean]

make_layer_plot <- function(l_key, is_top = FALSE, is_bottom = FALSE) {
  sub_raw <- tl[layer == l_key]
  sub_cv <- curves_empirical[layer == l_key]
  col_main <- COLORS[l_key]
  col_tint <- TINTS[l_key]
  lbl <- LAYER_NAMES[l_key]
  gap_bridge <- sub_cv[seg == "gap"]

  p <- ggplot() +
    annotate("rect", xmin = 0.52, xmax = 0.62, ymin = -Inf, ymax = Inf,
             fill = "#D9D9D9", alpha = 0.35) +
    geom_point(data = sub_raw, aes(x, val), colour = NEUTRAL, size = 0.38, alpha = 0.35, shape = 16) +
    geom_ribbon(data = sub_cv[seg != "gap"], aes(x, ymin = lower_95, ymax = upper_95, group = seg),
                fill = col_tint, alpha = 0.65) +
    geom_line(data = sub_cv[seg != "gap"], aes(x, fitted, group = seg),
              colour = col_main, linewidth = 0.70) +
    geom_line(data = gap_bridge, aes(x, fitted),
              colour = col_main, linewidth = 0.60, linetype = "22") +
    scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0)) +
    scale_y_continuous(labels = fmt_minus) +
    labs(x = NULL, y = lbl) +
    theme_nature(base_size = 5.8) +
    theme(
      plot.margin = margin(if (is_top) 1 else 0.5, 2, if (is_bottom) 5 else 0.5, 1, "mm"),
      axis.title.y = element_text(size = 5.5, family = FONT, colour = col_main, face = "bold"),
      axis.text.y = element_text(size = 5.0, family = FONT, colour = INK)
    )

  if (!is_bottom) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   axis.line.x = element_blank())
  } else {
    ticks_x <- c(0.00, 0.52, 0.62, 1.00)
    ticks_lab <- c("GW14", "Birth", "11 y", "100 y")
    p <- p +
      scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = ticks_x, labels = ticks_lab) +
      theme(axis.text.x = element_text(size = 5.2, family = FONT, colour = INK),
            axis.title.x = element_text(size = 5.4, family = FONT, colour = INK),
            axis.line.x = element_line(linewidth = 0.35, colour = INK)) +
      labs(x = "Age (developmental stage)")
  }

  if (is_top) {
    p <- p + labs(tag = "e") + theme(plot.tag = element_text(size = 8, face = "bold", family = FONT, colour = INK))
  }
  p
}

# --- a: APA vs mRNA coupling, C1 and C2 density contours ---
fa <- fread(file.path(D2, "fig2a_event_coupling.tsv"), select = c("cl", "apa_early", "apa_gap", "ex_early", "ex_gap"))
fa[, `:=`(pdui_dev = apa_early + apa_gap, expr_dev = ex_early + ex_gap, grp = fifelse(cl == 1, "C1", "C2"))]
rho <- fa[, .(n = .N, rho = cor(pdui_dev, expr_dev, method = "spearman")), by = grp][order(grp)]

lim_x <- quantile(fa$expr_dev, c(0.001, 0.999))
lim_y <- quantile(fa$pdui_dev, c(0.001, 0.999))
fa_v <- copy(fa[expr_dev > lim_x[1] & expr_dev < lim_x[2] & pdui_dev > lim_y[1] & pdui_dev < lim_y[2]])

res_c1 <- kde_gemini(fa_v[grp == "C1"], "C1")
res_c2 <- kde_gemini(fa_v[grp == "C2"], "C2")

fills <- rbind(res_c1$fill, res_c2$fill)
lines <- rbind(res_c1$line, res_c2$line)
fills[, lev_rank := as.numeric(factor(level)), by = grp]

pa <- ggplot() +
  ggrastr::rasterise(geom_point(data = fa_v, aes(expr_dev, pdui_dev, colour = grp),
                                size = 0.35, alpha = 0.18, shape = 16), dpi = 600) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "22") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "22") +
  geom_polygon(data = fills, aes(x, y, group = id, fill = grp, alpha = lev_rank), colour = NA) +
  geom_path(data = lines, aes(x, y, group = id, colour = grp), linewidth = 0.35, alpha = 0.9) +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_fill_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_alpha_continuous(range = c(0.14, 0.28), guide = "none") +
  annotate("text", x = lim_x[1] + 0.06 * diff(lim_x), y = lim_y[2] - 0.05 * diff(lim_y), hjust = 0, vjust = 1, size = 5.6 / .pt,
           colour = C1_COLOR, fontface = "bold", family = FONT, label = paste0("C1 \u03C1 = ", fmt2(rho$rho[1]))) +
  annotate("text", x = lim_x[1] + 0.06 * diff(lim_x), y = lim_y[2] - 0.16 * diff(lim_y), hjust = 0, vjust = 1, size = 5.6 / .pt,
           colour = C2_COLOR, fontface = "bold", family = FONT, label = paste0("C2 \u03C1 = ", fmt2(rho$rho[2]))) +
  scale_x_continuous(labels = fmt_minus) +
  scale_y_continuous(labels = fmt_minus) +
  coord_cartesian(xlim = lim_x, ylim = lim_y, expand = FALSE) +
  labs(x = "\u0394 mRNA (log2 CPM)", y = "\u0394 PDUI (logit)") +
  theme(plot.margin = margin(1, 2, 1, 1, "mm"))

# --- b: cis vs trans proteome QQ-plot (merged with pi1) ---
qq <- fread(file.path(DF, "fig2d_qq_quantiles.tsv"))
s1p <- fread(file.path(D3, "S1_apa_prot_cis_vs_trans_371.tsv"))
N85 <- as.integer(pv("n_q05"))
xmax <- max(qq$null_med) * 1.03; ymax <- max(qq$null_q95, qq$obs_nlogp) * 1.03

pb <- ggplot(qq, aes(null_med, obs_nlogp)) +
  geom_ribbon(aes(ymin = null_q05, ymax = null_q95), fill = NEUTRAL_LIGHT, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "22") +
  geom_point(data = qq[rank > N85], colour = NEUTRAL, size = 0.25, alpha = 0.5, shape = 16) +
  geom_point(data = qq[rank <= N85], colour = INK, size = 0.45, shape = 16) +
  annotate("text", x = 0.06, y = ymax * 0.99, hjust = 0, vjust = 1, size = 5.2 / .pt, colour = INK, family = FONT,
           label = sprintf("\u03C01 = %.3f (null %.3f)", pv("pi1"), nv("pi1"))) +
  annotate("text", x = 0.06, y = ymax * 0.86, hjust = 0, vjust = 1, size = 5.2 / .pt, colour = INK, family = FONT,
           label = sprintf("%d pairs FDR<0.05 (null %.1f)", N85, nv("n_q05"))) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
  labs(x = "Null \u2212log10 P", y = "Obs. \u2212log10 P") +
  theme(plot.margin = margin(1, 1, 1, 2, "mm"))

# --- c: bubble matrix, motif/RBP x contrast ---
aU <- fread(file.path(D3, "S2_aUTR_vs_cUTR_elem.tsv"))
el <- fread(file.path(D3, "S2_C1C2_vs_nondynamic_elem.tsv")); setnames(el, c("motif", "group", "est", "lo", "hi", "p"))
cl2 <- fread(file.path(D3, "S2_C1C2_vs_nondynamic_clip_adj.tsv")); setnames(cl2, c("rbp", "group", "est", "lo", "hi", "p"))
LAB <- c(CPE_UUUUAU = "CPE", PUM_UGUAHAUA = "PUM", ARE_AUUUA = "ARE", mir7_brain = "miR-7",
         mir7_brain_ctrl = "miR-7 ctrl", PUM1 = "PUM1", PUM2 = "PUM2", FMR1 = "FMR1", AGO2_brain = "AGO2")
mc <- rbind(
  aU[, .(feat = LAB[motif], contrast = "aUTR vs cUTR", rr = rr_AUadj, p = p_AUadj)],
  el[, .(feat = LAB[motif], contrast = fifelse(group == "C1", "C1 vs non-dyn.", "C2 vs non-dyn."), rr = est, p = p)],
  cl2[, .(feat = LAB[rbp], contrast = fifelse(group == "C1", "C1 vs non-dyn.", "C2 vs non-dyn."), rr = est, p = p)])
mc <- mc[!is.na(feat)]
mc[, feat := factor(feat, rev(c("CPE", "PUM", "ARE", "miR-7", "miR-7 ctrl", "PUM1", "PUM2", "FMR1", "AGO2")))]
mc[, contrast := factor(contrast, c("aUTR vs cUTR", "C1 vs non-dyn.", "C2 vs non-dyn."))]
mc[, block := fifelse(feat %in% c("CPE", "PUM", "ARE", "miR-7", "miR-7 ctrl"), "Sequence\nmotifs", "CLIP RBP\n(brain)")]
mc[, block := factor(block, c("Sequence\nmotifs", "CLIP RBP\n(brain)"))]
mc[, eff := pmin(abs(log2(rr)), 1.2)][, nlp := pmin(-log10(p), 60)]
mc[, star := fifelse(p < 0.001, "***", fifelse(p < 0.01, "**", fifelse(p < 0.05, "*", "")))]

pc <- ggplot(mc, aes(contrast, feat)) +
  geom_point(aes(size = eff, fill = nlp), shape = 21, colour = INK, stroke = 0.28) +
  geom_text(aes(label = star), size = 2.4, family = FONT, colour = INK, hjust = -0.45, vjust = 0.5) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#D9E2E8", high = "#1D4E75", limits = c(0, 60), breaks = c(0, 20, 40, 60),
                      name = expression(-log[10](italic(p)))) +
  scale_size_continuous(range = c(1.6, 4.6), breaks = c(0.3, 0.6, 0.9, 1.2), limits = c(0, 1.2),
                        labels = c("0.3", "0.6", "0.9", "1.2"), name = expression("|"*log[2]*"(RR)|")) +
  scale_x_discrete(labels = c("aUTR vs cUTR" = "aUTR vs cUTR",
                              "C1 vs non-dyn." = "C1 vs non-dyn.",
                              "C2 vs non-dyn." = "C2 vs non-dyn.")) +
  guides(fill = guide_colorbar(order = 1, barwidth = unit(2.0, "mm"), barheight = unit(13, "mm"),
                               title.theme = element_text(size = 5.2, family = FONT),
                               label.theme = element_text(size = 4.8, family = FONT)),
         size = guide_legend(order = 2, override.aes = list(fill = "#1D4E75", colour = INK, stroke = 0.28),
                             title.theme = element_text(size = 5.2, family = FONT),
                             label.theme = element_text(size = 4.8, family = FONT))) +
  labs(x = NULL, y = NULL) +
  theme(strip.background = element_rect(fill = "#EFEEEA", colour = NA),
        strip.text.y = element_text(size = 5.6, face = "bold", family = FONT, colour = INK, lineheight = 0.9, margin = margin(r = 2.5, l = 2.5)),
        axis.text.x = element_text(size = 5.8, family = FONT, colour = INK),
        axis.text.y = element_text(size = 5.8, family = FONT, colour = INK),
        legend.title = element_text(size = 5.2, family = FONT),
        legend.text = element_text(size = 4.8, family = FONT),
        panel.spacing.y = unit(1.8, "mm"),
        plot.margin = margin(1, 1, 1, 1, "mm"))

# --- d: gene model ---
img_d <- png::readPNG(file.path(FIG_DIR, "fig2_d.png"))
pd <- ggplot() +
  annotation_raster(img_d, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
  theme_void() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(tag = "d") +
  theme(plot.tag = element_text(size = 8, face = "bold", family = FONT, colour = INK),
        plot.margin = margin(0, 1, 0, 1, "mm"))

# --- e: layer trajectories, empirical LOESS with 95% CI ---
apa <- readRDS(file.path(D3, "apa_matched_raw.rds"))
rna <- readRDS(file.path(D3, "rna_matched_raw.rds"))
prt <- readRDS(file.path(D3, "protein_matched_raw.rds"))
man <- fread(file.path(D3, "matched_sample_manifest_raw.tsv"))

rep_site <- "ENST00000675467.1|PIK3CA|chr3:179234094-179235371|+"
rid_rna <- grep("^ENSG00000121879", rownames(rna), value = TRUE)

td <- data.table(
  sample = colnames(apa),
  age_years = man$age_years,
  stage = man$life_stage,
  apa = as.numeric(apa[rep_site, ]),
  rna = as.numeric(rna[rid_rna, ]),
  prot = as.numeric(prt["P42336", ])
)
ay <- td$age_years
x_age <- numeric(length(ay))
is_pre <- ay <= 0
x_age[is_pre] <- 0.00 + 0.52 * (ay[is_pre] - (-0.5005)) / (0 - (-0.5005))
is_post <- ay >= 11.6
x_age[is_post] <- 0.62 + 0.38 * (log10(ay[is_post]) - log10(11.6)) / (log10(102) - log10(11.6))
is_gap <- ay > 0 & ay < 11.6
x_age[is_gap] <- 0.52 + 0.10 * (ay[is_gap] - 0) / (11.6 - 0)
td[, x := x_age]
rm(apa, rna, prt); invisible(gc())

tl <- melt(td, id.vars = c("sample", "age_years", "stage", "x"),
           measure.vars = c("apa", "rna", "prot"), variable.name = "layer", value.name = "val")
tl <- tl[is.finite(val)]

# loess fit per layer; prenatal and postnatal (>= 11.6 y) segments, dashed bridge across the gap
curves_list <- lapply(c("apa", "rna", "prot"), function(l) {
  sub_dt <- tl[layer == l]
  n_grid <- 80
  d_pre <- sub_dt[stage == "Prenatal"]
  fit_pre <- loess(val ~ x, data = d_pre, span = 0.65, degree = 2)
  grid_pre <- data.table(x = seq(min(d_pre$x), max(d_pre$x), length.out = n_grid), seg = "pre")
  pr_pre <- predict(fit_pre, newdata = grid_pre, se = TRUE)
  grid_pre[, `:=`(fitted = as.numeric(pr_pre$fit), se = as.numeric(pr_pre$se.fit))]

  d_post <- sub_dt[stage == "Postnatal" & age_years >= 11.6]
  fit_post <- loess(val ~ x, data = d_post, span = 0.65, degree = 2)
  grid_post <- data.table(x = seq(min(d_post$x), max(d_post$x), length.out = n_grid), seg = "post")
  pr_post <- predict(fit_post, newdata = grid_post, se = TRUE)
  grid_post[, `:=`(fitted = as.numeric(pr_post$fit), se = as.numeric(pr_post$se.fit))]

  grid_gap <- data.table(x = seq(max(d_pre$x), min(d_post$x), length.out = 30), seg = "gap")
  y_start <- grid_pre[.N, fitted]; y_end <- grid_post[1, fitted]
  grid_gap[, fitted := y_start + (x - min(x)) / (max(x) - min(x)) * (y_end - y_start)]
  grid_gap[, se := mean(c(grid_pre[.N, se], grid_post[1, se]))]

  d <- rbind(grid_pre, grid_gap, grid_post)
  d[, `:=`(lower_95 = fitted - 1.96 * se, upper_95 = fitted + 1.96 * se)]
  d[, layer := l]
  d
})
curves_empirical <- rbindlist(curves_list)

COLORS <- c(
  apa  = "#009E73",  # APA Teal-green
  rna  = "#2E5B7E",  # mRNA Slate-blue
  prot = "#D17A4A"   # Protein Warm-orange
)
TINTS <- c(
  apa  = "#A3D9C9",
  rna  = "#C9D4DE",
  prot = "#F2D2BE"
)
LAYER_NAMES <- c(
  apa  = "APA (PDUI logit)",
  rna  = "mRNA (log2 CPM)",
  prot = "Protein (abundance)"
)

pe1 <- make_layer_plot("apa", is_top = TRUE)
pe2 <- make_layer_plot("rna")
pe3 <- make_layer_plot("prot", is_bottom = TRUE)

pe_stack <- (pe1 / pe2 / pe3) + plot_layout(heights = c(1, 1, 1.15))

# --- assemble ---
row1 <- (tg(pa, "a") | tg(pb, "b")) + plot_layout(widths = c(1.38, 1.00))
row2 <- tg(pc, "c")
row3 <- (pd | pe_stack) + plot_layout(widths = c(1.25, 1.00))

fig <- (row1 / row2 / row3) + plot_layout(heights = c(0.78, 0.85, 1.45))

FIG_W <- 168; FIG_H <- 176
ggsave(file.path(FIG_DIR, "fig2.png"), fig, width = FIG_W, height = FIG_H, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig2.pdf"), fig, width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
