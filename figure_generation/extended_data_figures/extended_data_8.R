# ED Fig 8: genotype-by-stage interaction genes: detection route (a), gating scatter (b), expression shift by class (c), ATP2B1 saturation (d), postnatal resolvability (e).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel); library(ragg); library(ggbeeswarm)})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
D4 <- file.path(ROOT, "data/04_qtl"); D5 <- file.path(ROOT, "data/05_gxs")
CLS <- c("Reversed", "Amplified", "Attenuated"); CC <- CLASS_COLORS[CLS]
SC_P <- STAGE_COLORS[["Prenatal"]]; SC_Q <- STAGE_COLORS[["Postnatal"]]
theme_set(theme_nature(base_size = 6.5) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 6.5, colour = INK),
  legend.text = element_text(size = 6, colour = INK), legend.key.size = unit(2.2, "mm"),
  plot.margin = margin(1, 1.5, 1, 1, "mm"),
  plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1),
  strip.text = element_text(size = 6.5, face = "bold")))
T6 <- 6 / .pt
tg <- function(p, t) p + labs(tag = t)
lab_m <- function(x) sub("^-", "\u2212", format(x, trim = TRUE, drop0trailing = TRUE))
fmt_p <- function(p) { if (p >= 0.01) return(sprintf("P = %s", formatC(signif(p, 2), format = "fg", digits = 2, flag = "#")))
  e <- floor(log10(p)); m <- round(p / 10^e, 1); if (m >= 10) { m <- 1; e <- e + 1 }
  sup <- c("0"="\u2070","1"="\u00b9","2"="\u00b2","3"="\u00b3","4"="\u2074","5"="\u2075","6"="\u2076","7"="\u2077","8"="\u2078","9"="\u2079","-"="\u207b")
  sprintf("P = %.1f \u00d7 10%s", m, paste(sup[strsplit(as.character(e), "")[[1]]], collapse = "")) }

# --- data ---
core <- fread(file.path(D5, "intx_hc3_core_groups.tsv"))
core[, cls := fifelse(sign(beta_pre) != sign(beta_post), "Reversed",
               fifelse(abs(beta_post) > abs(beta_pre), "Amplified", "Attenuated"))]
stopifnot(nrow(core) == 193)
core[, route := factor(c(`Shared group` = "Both maps", `Postnatal detected group` = "Postnatal map",
                         `Prenatal detected group` = "Prenatal map", `Interaction-only group` = "Neither")[stage_group_status],
                       c("Neither", "Prenatal map", "Postnatal map", "Both maps"))]
stopifnot(!anyNA(core$route))
s4 <- fread(file.path(D4, "S4_gating_eqtl_fixed_RNA.tsv"))[!is.na(dexpr_v3)]
d  <- merge(s4[cls %in% CLS], core[, .(gene_id, route)], by = "gene_id")
stopifnot(nrow(d) == 191)
d[, cls := factor(cls, CLS)]

# --- a: detection route of the 193 genes, by class ---
ra <- core[, .N, by = .(route, cls)][, cls := factor(cls, rev(CLS))]
tot <- core[, .N, by = route]
pa <- ggplot(ra, aes(N, route, fill = cls)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(data = ra[N >= 9], aes(label = N), position = position_stack(vjust = 0.5),
            size = T6, colour = "white", fontface = "bold") +
  geom_text(data = tot, aes(N, route, label = N), inherit.aes = FALSE, hjust = -0.25, size = T6, colour = INK) +
  scale_fill_manual(values = CC, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)), breaks = c(0, 25, 50, 75)) +
  labs(x = "G\u00d7Stage genes", y = "Detected in") +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
pa <- tg(pa, "a")

# --- b: gating scatter, change in |beta| vs change in expression ---
rt <- cor.test(d$dabs, d$dexpr_v3, method = "spearman", exact = FALSE)
lab <- d[gene_id %in% c("TESPA1", "SLCO1A2")]
leg <- data.table(x = 6.3, y = c(-0.92, -1.06, -1.20), k = CLS)
pb <- ggplot(d, aes(dexpr_v3, dabs)) +
  geom_hline(yintercept = 0, colour = GRID_COL, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = GRID_COL, linewidth = 0.3) +
  geom_point(aes(fill = cls), shape = 21, colour = "white", stroke = 0.25, size = 1.5) +
  geom_text_repel(data = lab, aes(label = gene_id), fontface = "italic", size = T6, colour = INK,
                  nudge_y = 0.22, min.segment.length = 0, segment.size = 0.25, segment.colour = MUTED, seed = 1) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = T6, colour = INK,
           label = sprintf("\u03c1 = %.2f, %s", unname(rt$estimate), fmt_p(rt$p.value))) +
  geom_point(data = leg, aes(x, y, fill = k), shape = 21, colour = "white", stroke = 0.25, size = 1.6, inherit.aes = FALSE) +
  geom_text(data = leg, aes(x + 0.35, y, label = k), hjust = 0, size = T6, colour = INK, inherit.aes = FALSE) +
  scale_fill_manual(values = CC, guide = "none") + scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(breaks = seq(-2, 8, 2), labels = lab_m) + scale_y_continuous(labels = lab_m) +
  labs(x = "\u0394log\u2082CPM (post \u2212 pre)", y = "\u0394|\u03b2| (post \u2212 pre)")
pb <- tg(pb, "b")

# --- c: expression shift by class vs shared apaQTL genes ---
s4c <- copy(s4)
s4c[, grp := factor(fifelse(group == "Shared apaQTL", "Shared", group), rev(c("Shared", CLS)))]
stopifnot(!anyNA(s4c$grp))
wam <- wilcox.test(s4c[grp == "Amplified", dexpr_v3], s4c[grp == "Attenuated", dexpr_v3], exact = FALSE)
s4c <- s4c[grp != "Shared"][, grp := droplevels(grp)]
med <- s4c[, .(m = median(dexpr_v3)), by = grp]
nn <- s4c[, .N, by = grp]
GC <- c(CC, Shared = NEUTRAL)
ylab_c <- setNames(sprintf("%s (%s)", nn$grp, formatC(nn$N, big.mark = ",", format = "d")), as.character(nn$grp))
pc <- ggplot(s4c, aes(dexpr_v3, grp)) +
  geom_vline(xintercept = 0, colour = GRID_COL, linewidth = 0.3) +
  geom_quasirandom(aes(colour = grp), groupOnX = FALSE, width = 0.32, size = 0.75, alpha = 0.85) +
  geom_segment(data = med, aes(x = m, xend = m, y = as.numeric(grp) - 0.38, yend = as.numeric(grp) + 0.38),
               linewidth = 0.7, colour = INK, inherit.aes = FALSE) +
  annotate("segment", x = 7.6, xend = 7.6, y = 1, yend = 2, linewidth = 0.25, colour = INK) +
  annotate("text", x = 7.85, y = 1.5, label = fmt_p(wam$p.value), hjust = 0, size = T6, colour = INK) +
  scale_colour_manual(values = GC, guide = "none") +
  scale_y_discrete(labels = ylab_c) +
  scale_x_continuous(breaks = seq(-4, 8, 4), labels = lab_m, expand = expansion(mult = c(0.03, 0.2))) +
  labs(x = "\u0394log\u2082CPM (post \u2212 pre)", y = NULL) +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
pc <- tg(pc, "c")

# --- d: ATP2B1 saturation by genotype and stage ---
diag <- fread(file.path(D4, "atp2b1_stage_diagnostics.tsv"))
gt <- c("C/C", "C/A", "A/A")
dd <- rbindlist(lapply(seq_len(nrow(diag)), function(i) data.table(
  stage = diag$stage[i], g = gt, n = unlist(diag[i, .(geno0_n, geno1_n, geno2_n)]),
  f = unlist(diag[i, .(geno0_frac_saturated, geno1_frac_saturated, geno2_frac_saturated)]))))
dd[, stage := factor(stage, c("prenatal", "postnatal"), c("Prenatal", "Postnatal"))]
dd[, g := factor(g, gt)]
stopifnot(dd[stage == "Postnatal", sum(n)] == 508, dd[stage == "Prenatal", sum(n)] == 189)
dl <- rbind(dd[, .(stage, g, n, part = "Complete", pct = 100 * f)], dd[, .(stage, g, n, part = "Intermediate", pct = 100 * (1 - f))])
dl[, part := factor(part, c("Intermediate", "Complete"))]
dl[, fillc := ifelse(part == "Complete", ifelse(stage == "Prenatal", SC_P, SC_Q),
                     ifelse(stage == "Prenatal", "#F6E6C6", "#DCEBF5"))]
pd_lab <- data.table(stage = factor(c("Prenatal", "Postnatal"), levels(dd$stage)),
                     lab = c(fmt_p(diag[stage == "prenatal", saturated_p]), fmt_p(diag[stage == "postnatal", saturated_p])))
pdd <- ggplot(dl, aes(g, pct)) +
  geom_col(aes(fill = I(fillc), group = part), width = 0.7, colour = "white", linewidth = 0.25) +
  geom_text(data = dd, aes(g, 100 * f / 2, label = sprintf("%.0f", 100 * f)), size = T6, colour = "white", fontface = "bold") +
  geom_text(data = dd, aes(g, -9, label = paste0("n = ", n)), size = T6, colour = MUTED) +
  geom_text(data = pd_lab, aes(x = 2, y = 108, label = lab), size = T6, colour = INK) +
  facet_wrap(~stage, nrow = 1) +
  scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = 0)) +
  coord_cartesian(ylim = c(-14, 113), clip = "off") +
  labs(x = expression(italic(ATP2B1)~genotype), y = "Complete distal usage (%)") +
  theme(panel.spacing = unit(2, "mm"))
pdd <- pdd + theme(strip.text = element_text(size = 6.5, face = "bold", colour = INK))
pdd <- tg(pdd, "d")

# --- e: postnatal individuals with intermediate PDUI, per class ---
bd <- fread(file.path(D5, "boundary_refit_hc3_core.tsv")); bd <- merge(bd, core[, .(gene_id, cls)], by = "gene_id")
stopifnot(nrow(bd) == 193)
bd[, cls := factor(cls, rev(CLS))]
bd[, sat := post_n_kept < 10]
gs <- bd[, .(k = sum(sat), n = .N), by = cls]
m  <- matrix(c(bd[cls == "Amplified", sum(sat)], bd[cls == "Amplified", sum(!sat)],
               bd[cls != "Amplified", sum(sat)], bd[cls != "Amplified", sum(!sat)]), 2, byrow = TRUE)
ft <- fisher.test(m)
pe <- ggplot(bd, aes(post_n_kept, cls)) +
  annotate("rect", xmin = -Inf, xmax = 9.5, ymin = -Inf, ymax = Inf, fill = "#F3F3F3") +
  geom_vline(xintercept = 9.5, linetype = 2, linewidth = 0.3, colour = MUTED) +
  geom_quasirandom(aes(colour = cls), groupOnX = FALSE, width = 0.36, size = 0.85, alpha = 0.9) +
  annotate("text", x = 0, y = 3.62, label = "Saturated", hjust = 0, size = T6, colour = MUTED, fontface = "italic") +
  annotate("text", x = 40, y = 3.62, label = fmt_p(ft$p.value), hjust = 0.5, size = T6, colour = INK) +
  scale_colour_manual(values = CC, guide = "none") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 3, 10, 30, 100, 300, 508), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_discrete(expand = expansion(add = c(0.55, 0.85)),
                   labels = setNames(sprintf("%s\n%d/%d saturated", gs$cls, gs$k, gs$n), as.character(gs$cls))) +
  labs(x = "Postnatal individuals with 0 < PDUI < 1", y = NULL) +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_text(lineheight = 0.95, hjust = 1))
pe <- tg(pe, "e")

# --- assemble (180 x 108 mm) ---
row1 <- pa + pb + pc + plot_layout(widths = c(36, 80, 52))
row2 <- pdd + pe + plot_layout(widths = c(70, 100))
fig  <- row1 / row2 + plot_layout(heights = c(62, 40))
ggsave(file.path(FIG_DIR, "ed8.png"), fig, width = 180, height = 108, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed8.pdf"), fig, width = 180, height = 108, units = "mm", device = cairo_pdf)
