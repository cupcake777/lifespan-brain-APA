# ED Fig 7: apaQTL and eQTL cross-stage effects (a, b), prenatal strength by replication (c), validation rates (d), lead distance to APA boundary (e).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)
                                library(scattermore); library(ragg); library(readxl)})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
D4 <- file.path(ROOT, "data/04_qtl"); DF <- file.path(ROOT, "data/figures/ed7")
theme_set(theme_nature(base_size = 7) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
  axis.ticks.length = unit(0.8, "mm"), legend.text = element_text(size = 6, colour = INK),
  legend.title = element_text(size = 6, colour = INK), legend.key.size = unit(2, "mm"),
  legend.margin = margin(0, 0, 0, 0), plot.margin = margin(1, 2, 1, 1, "mm"),
  plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1)))
T55 <- 5.5 / .pt
TSTAT <- 7.5 / .pt          # stat boxes in a and b
SC  <- c(Prenatal = STAGE_COLORS[["Prenatal"]], Postnatal = STAGE_COLORS[["Postnatal"]])
SCT <- c(Prenatal = STAGE_TEXT[["Prenatal"]], Postnatal = STAGE_TEXT[["Postnatal"]])
tg <- function(p, t) p + labs(tag = t)
src <- list()

# shared builder for the two cross-stage scatters (a, b): same binning, colour ramp and stat box
cross_stage_scatter <- function(x, y, ylab_x, ylab_y, box_lines, grid_n = 46L) {
  lim <- max(abs(c(x, y))) * 1.05
  gx <- cut(x, breaks = seq(-lim, lim, length.out = grid_n + 1), include.lowest = TRUE)
  gy <- cut(y, breaks = seq(-lim, lim, length.out = grid_n + 1), include.lowest = TRUE)
  dens <- as.data.table(as.data.frame(table(gx, gy)))
  dens[, `:=`(xmid = (as.numeric(gx) - 0.5) * (2 * lim / grid_n) + -lim,
              ymid = (as.numeric(gy) - 0.5) * (2 * lim / grid_n) + -lim)]
  dens[, cnt := as.numeric(Freq)]
  dens <- dens[cnt > 0]
  dens[, w := sqrt(cnt / max(cnt))]         # sqrt: single-pair cells visible, dense cells not saturated
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = NEUTRAL_LIGHT) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = NEUTRAL_LIGHT) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.25, colour = MUTED) +
    geom_tile(data = dens, aes(x = xmid, y = ymid, fill = w),
              width = 2 * lim / grid_n * 1.02, height = 2 * lim / grid_n * 1.02) +
    scale_fill_gradientn(colours = c("#F7F9FB", "#C3D6E8", "#7FA9CE", "#2E6699", "#123B63"),
                         limits = c(0, 1), guide = "none") +
    geom_scattermore(pointsize = 0.22, alpha = 0.10, stroke = 0, colour = INK,
                     pixels = c(600, 600)) +
    annotate("text", x = -lim * 0.97, y = lim * 0.97, hjust = 0, vjust = 1,
             size = TSTAT, colour = INK, lineheight = 1.15, label = box_lines) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("\u22121", "0", "1")) +
    scale_y_continuous(breaks = c(-1, 0, 1), labels = c("\u22121", "0", "1")) +
    labs(x = ylab_x, y = ylab_y)
}

# --- a: apaQTL cross-stage effects (pairs significant in both stages) ---
ga <- fread(file.path(DF, "ed7a_variant_gene_pairs.tsv"))
n_a <- nrow(ga)
same_a  <- mean(sign(ga$slope_pre) == sign(ga$slope_post))
r2_a <- mean(pmax(abs(ga$slope_pre) / abs(ga$slope_post),
                  abs(ga$slope_post) / abs(ga$slope_pre)) <= 2)   # within 2-fold: ratio of slopes (no s.e. in file)
rho_a <- cor(ga$slope_pre, ga$slope_post, method = "spearman")
stopifnot(all(ga$pval_nominal_pre < 0.05, ga$pval_nominal_post < 0.05))
stopifnot(n_a == 49957, abs(same_a - 0.978181) < 1e-4)
src$a <- "data/figures/ed7/ed7a_variant_gene_pairs.tsv (49,957 pairs significant in both stages; statistics computed here)"
ms <- fread(file.path(D4, "05_mash_summary.tsv"))
msh <- ms[metric == "pairwise_sharing_by_magnitude"]
SH_A <- msh[qtl == "apaQTL", value]; stopifnot(length(SH_A) == 1)
sh_e <- msh[qtl == "eQTL", value]
pa <- cross_stage_scatter(
  ga$slope_pre, ga$slope_post,
  "Prenatal apaQTL effect (slope)", "Postnatal apaQTL effect (slope)",
  sprintf("apaQTL  n = %s\nSame sign: %.1f%%\nWithin 2-fold: %.1f%%\n\u03c1 = %.2f\nmashr shared: %.1f%%",
          format(n_a, big.mark = ","), 100 * same_a, 100 * r2_a, rho_a, 100 * SH_A))
pa <- tg(pa, "a")

# --- b: eQTL cross-stage effects, same construction as a ---
eb <- fread(file.path(D4, "05_mash_posterior_eqtl.tsv"))
eb <- eb[abs(b_prenatal) / s_prenatal > 1.96 & abs(b_postnatal) / s_postnatal > 1.96]
n_b <- nrow(eb)
same_b <- mean(sign(eb$b_prenatal) == sign(eb$b_postnatal))
r2_b   <- mean(pmax(abs(eb$b_prenatal) / abs(eb$b_postnatal),
                    abs(eb$b_postnatal) / abs(eb$b_prenatal)) <= 2)
rho_b  <- cor(eb$b_prenatal, eb$b_postnatal, method = "spearman")
stopifnot(n_b == 14396)
src$b <- sprintf("data/04_qtl/05_mash_posterior_eqtl.tsv (%s pairs significant in both stages, |b|/s > 1.96; statistics computed here); data/04_qtl/05_mash_summary.tsv (mashr pairwise sharing by magnitude, strong sets)", format(n_b, big.mark = ","))
stopifnot(length(sh_e) == 1, round(100 * SH_A, 1) == 85.0, round(100 * sh_e, 1) == 63.3)
pb <- cross_stage_scatter(
  eb$b_prenatal, eb$b_postnatal,
  "Prenatal eQTL effect (\u03b2)", "Postnatal eQTL effect (\u03b2)",
  sprintf("eQTL  n = %s\nSame sign: %.1f%%\nWithin 2-fold: %.1f%%\n\u03c1 = %.2f\nmashr shared: %.1f%%",
          format(n_b, big.mark = ","), 100 * same_b, 100 * r2_b, rho_b, 100 * sh_e))
pb <- tg(pb, "b")

# --- c: prenatal |z| by postnatal replication frequency ---
rc <- fread(file.path(D4, "ed7_per_gene_recovery.tsv"))
rc[, grp := fifelse(n_recovered_seeds == 0, "0",
                    fifelse(n_recovered_seeds == 20, "20", "1-19"))]
rc[, grp := factor(grp, c("0", "1-19", "20"))]
nn <- rc[, .N, by = grp]
stopifnot(nn$N[match(c("0", "1-19", "20"), nn$grp)] == c(144, 267, 258), nrow(rc) == 669)
rho_c <- cor(rc$n_recovered_seeds, rc$abs_z_prenatal, method = "spearman")
src$c <- "data/04_qtl/ed7_per_gene_recovery.tsv (669 prenatal apaQTL genes; group counts 144 / 267 / 258 verified)"
pc <- ggplot(rc, aes(grp, abs_z_prenatal)) +
  geom_boxplot(outlier.shape = NA, width = 0.45, linewidth = 0.25, colour = INK, fill = NA) +
  geom_point(position = position_jitter(width = 0.18, height = 0, seed = 20261026), size = 0.4, alpha = 0.25, stroke = 0,
             colour = adjustcolor(STAGE_COLORS[["Prenatal"]], alpha.f = 0.55)) +
  geom_text(data = nn, aes(grp, Inf, label = paste0("n = ", N)), vjust = 1.5,
            size = T55, colour = MUTED, inherit.aes = FALSE) +
  annotate("text", x = 0.55, y = max(rc$abs_z_prenatal) * 0.97, hjust = 0, vjust = 1,
           size = T55, colour = INK, label = sprintf("\u03c1 = %.2f", rho_c)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "Postnatal replication frequency", y = "Prenatal |z|")
pc <- tg(pc, "c")

# --- d: validation rate of prenatal-specific and shared genes, 100% stacked ---
bx <- read_excel(file.path(D4, "braingvex_adult_bulk_comp.xlsx"), sheet = "prenatal_specific")
nx <- read_excel(file.path(D4, "nan2026_comparison.xlsx"), sheet = "any_celltype")
b_sh <- bx[bx$group == "shared", ]; b_pr <- bx[bx$group %in% c("absent_postnatal", "never_recovered"), ]
bgv <- data.table(map = "BrainGVEX", grp = c("Shared", "Prenatal-specific"),
                  num = c(round(b_sh$gene_detected_BrainGVEX * b_sh$n_genes),
                          sum(round(b_pr$gene_detected_BrainGVEX * b_pr$n_genes))),
                  den = c(sum(b_sh$n_genes), sum(b_pr$n_genes)))
n_sh <- nx[nx$group == "prenatal: shared", ]
n_ab <- nx[nx$group == "prenatal: absent from postnatal map (106)", ]
n_nv <- nx[nx$group == "prenatal: never recovered (144)", ]
agd <- data.table(map = "Aged snRNA", grp = c("Shared", "Prenatal-specific"),
                  num = c(round(n_sh$replicated_any * n_sh$genes_matched),
                          round(n_ab$replicated_any * n_ab$genes_matched) + round(n_nv$replicated_any * n_nv$genes_matched)),
                  den = c(n_sh$genes_matched, n_ab$genes_matched + n_nv$genes_matched))
stopifnot(bgv$num[1] == 381, bgv$den[1] == 518, bgv$num[2] == 43, bgv$den[2] == 151,
          agd$num[1] == 119, agd$den[1] == 325, agd$num[2] == 13, agd$den[2] == 83)
src$d <- "data/04_qtl/braingvex_adult_bulk_comp.xlsx (prenatal_specific); data/04_qtl/nan2026_comparison.xlsx (any_celltype)"
dd <- rbind(bgv, agd)
ci <- rbindlist(lapply(seq_len(nrow(dd)), function(i) {
  b <- binom.test(dd$num[i], dd$den[i])$conf.int
  data.table(lo = b[1], hi = b[2]) }))
dd <- cbind(dd, ci)
dd[, frac := num / den]
dd[, grp := factor(grp, c("Shared", "Prenatal-specific"))]
dd[, map := factor(map, c("BrainGVEX", "Aged snRNA"),
                   labels = c("Adult bulk PFC\n(BrainGVEX)", "Aged snRNA-seq\n(3' aQTL atlas)"))]
dd[, status := factor(ifelse(num > 0, "Detected", "Undetected"), c("Detected", "Undetected"))]
dd[, cnt_det := num]
dd[, cnt_und := den - num]
stk <- rbind(dd[, .(map, grp, status = "Detected", cnt = frac, num, den)],
             dd[, .(map, grp, status = "Undetected", cnt = 1 - frac, num, den)])
stk[, status := factor(status, c("Detected", "Undetected"))]
stk[, mcode := as.integer(map)]
stk[, gcode := as.integer(grp)]
stk[, key := paste0("m", mcode, "g", gcode)]
dd2 <- copy(dd)
dd2[, mcode := as.integer(map)]
dd2[, gcode := as.integer(grp)]
dd2[, key := paste0("m", mcode, "g", gcode)]
ordv <- c(paste0("m", 1, "g", 1), paste0("m", 1, "g", 2),
          paste0("m", 2, "g", 1), paste0("m", 2, "g", 2))
stk[, key := factor(key, levels = rev(ordv))]
dd2[, key := factor(key, levels = rev(ordv))]
pd <- ggplot(stk, aes(x = cnt, y = key, fill = status)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  # labels sit right of the 100% end (in axis units) so no segment boundary clips them
  geom_text(data = dd2, aes(x = 1.02, y = key,
                            label = sprintf("%.1f%% (%d/%d)", 100 * frac, num, den)),
            hjust = 0, size = 5.2 / .pt, colour = INK, inherit.aes = FALSE) +
  scale_fill_manual(values = c(Detected = SCT[["Prenatal"]], Undetected = NEUTRAL_LIGHT), name = NULL) +
  # axis to 1.31: the longest label needs the panel to reach ~1.25
  scale_x_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0, 1.31)) +
  scale_y_discrete(labels = c("m1g1" = "BrainGVEX\nShared", "m1g2" = "BrainGVEX\nPrenatal-specific",
                              "m2g1" = "Aged snRNA\nShared", "m2g2" = "Aged snRNA\nPrenatal-specific")) +
  labs(x = "Validation rate", y = NULL) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.4, "mm"), legend.text = element_text(size = 5.6),
        legend.margin = margin(-1, 0, 0, 0),
        axis.text.y = element_text(size = 5.6, lineheight = 0.95),
        plot.margin = margin(1, 2, 1, 1, "mm"))
pd <- tg(pd, "d")

# --- e: lead-variant distance to APA-interval end, by stage ---
de <- fread(file.path(DF, "ed7e_lead_distance_by_stage.tsv"))
de[, stage := factor(stage, c("Prenatal", "Postnatal"))]
ecd <- de[order(dist), .(dist = pmax(dist, 1), ecdf = seq_len(.N) / .N), by = stage]
src$e <- "data/figures/ed7/ed7e_lead_distance_by_stage.tsv"
pe <- ggplot(ecd, aes(dist, ecdf, colour = stage)) +
  geom_vline(xintercept = 5000, linetype = 2, linewidth = 0.25, colour = MUTED) +
  annotate("text", x = 5000, y = 0.06, label = "5 kb", hjust = 0, vjust = 0,
           size = 5 / .pt, colour = MUTED, lineheight = 0.9) +
  geom_step(linewidth = 0.35) +
  scale_colour_manual(values = SC, name = NULL) +
  scale_x_log10(breaks = c(1, 10, 100, 1e3, 1e4, 1e5, 1e6),
                labels = c("1", "10", "100", "1 kb", "10 kb", "100 kb", "1 Mb"),
                limits = c(1, max(ecd$dist) * 1.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Distance to nearest APA boundary (bp)", y = "Cumulative fraction") +
  theme(legend.position = c(0.03, 0.97), legend.justification = c(0, 1),
        legend.key.height = unit(2, "mm"),
        legend.background = element_rect(fill = "white", colour = NA))
pe <- tg(pe, "e")

# --- assemble ---
row1 <- (pa | pb | pc) + plot_layout(widths = c(1, 1, 1.05))
row2 <- (pd | pe) + plot_layout(widths = c(0.92, 1.08))
fig <- (row1 / plot_spacer() / row2) + plot_layout(heights = c(1, 0.075, 0.92))
ggsave(file.path(FIG_DIR, "ed7.png"), fig, width = 180, height = 126, units = "mm",
       dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed7.pdf"), fig, width = 180, height = 126, units = "mm",
       device = cairo_pdf)
writeLines(c("panel\tsource", paste(names(src), unlist(src), sep = "\t")),
           file.path(DF, "ed7_sources.tsv"))
