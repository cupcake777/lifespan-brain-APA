# Fig 6: ATP2B1, MTRFR, OGFOD2 and ARL3 loci, boxplots and PP4 strips (repaired coloc run).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
source(file.path(ROOT, "release/figures/locus.R"))
# --- helpers ---
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1))
df_of <- function(p, logp_col, lead) data.frame(pos = p$pos, logp = p[[logp_col]], ld_bin = p$ld_bin, is_lead = p$snp == lead)
mk_track <- function(df, ylab, pp4 = NULL, thr = FALSE) list(df = df, ylab = ylab, pp4 = pp4, threshold = thr)
# regional points for one comparison; lead_rs must be in r2 >= 0.99 with the apaQTL lead of the exported table
reg <- function(g, gw, st, lead_rs) {
  b <- CO[gene == g & GWAS_ID == gw & QTL_ID == st][which.max(PP4)]; stopifnot(nrow(b) == 1)
  d <- RT[gene == g & gwas_id == gw & qtl_stage == st]
  stopifnot(nrow(d) > 0, d$locus[1] == b$Locus, lead_rs %in% d$snp, d[snp == lead_rs, r2_to_apaqtl_lead] >= 0.99)
  p <- d[, .(snp, pos, gwas_logp = -log10(gwas_p), qtl_logp = -log10(qtl_p), r2 = r2_to_apaqtl_lead)]
  ld <- factor(as.character(cut(p$r2, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"), right = FALSE)), LOCUS_LD_LEVELS)
  ld[is.na(p$r2)] <- factor("missing LD", LOCUS_LD_LEVELS)
  p[, ld_bin := ld][, is_lead := snp == lead_rs]
  list(points = as.data.frame(p), pp4 = b$PP4,
       meta = list(chrom = sub(":.*", "", tstrsplit(b$Phenotype, "|", fixed = TRUE)[[3]]), min_pos = min(p$pos), max_pos = max(p$pos), phenotype = b$Phenotype, pp4 = b$PP4))
}
boxes <- function(gene_, layer_, ylab, ann_pre, ann_post, cap = NULL, ann_y = 0.27, top = 0.30) {
  b <- BX[gene == gene_ & layer == layer_]
  b[, stage := factor(stage, c("Prenatal", "Postnatal"))][, x := dosage + 1 + fifelse(stage == "Prenatal", -0.19, 0.19)]
  gl <- unique(b[, .(dosage, genotype)])[order(dosage)]
  lo <- min(b$whisker_lo); hi <- max(b$whisker_hi); rng <- hi - lo
  ggplot(b) +
    { if (!is.null(cap)) annotate("rect", xmin = -Inf, xmax = Inf, ymin = cap - rng * 0.012, ymax = cap + rng * 0.012, fill = "#E2E2E2") } +
    geom_segment(aes(x = x, xend = x, y = whisker_lo, yend = whisker_hi, colour = stage), linewidth = 0.3) +
    geom_rect(aes(xmin = x - 0.15, xmax = x + 0.15, ymin = q1, ymax = q3, fill = stage, colour = stage), linewidth = 0.3, alpha = 0.85) +
    geom_segment(aes(x = x - 0.15, xend = x + 0.15, y = median, yend = median), colour = "white", linewidth = 0.5) +
    geom_line(aes(x = x, y = mean, colour = stage, group = stage), linewidth = 0.45) +
    geom_point(aes(x = x, y = mean, colour = stage), shape = 21, fill = "white", size = 1.3, stroke = 0.45) +
    geom_text(aes(x = x, y = lo - rng * 0.10, label = n, colour = stage), size = 5 / .pt, show.legend = FALSE) +
    annotate("text", x = 0.55, y = hi + rng * ann_y, label = ann_pre, parse = TRUE, hjust = 0, vjust = 0.5, size = 5.2 / .pt, colour = unname(STX["Prenatal"])) +
    annotate("text", x = 0.55, y = hi + rng * (ann_y - 0.125), label = ann_post, parse = TRUE, hjust = 0, vjust = 0.5, size = 5.2 / .pt, colour = unname(STX["Postnatal"])) +
    scale_fill_manual(values = SC, name = NULL) + scale_colour_manual(values = STX, name = NULL) +
    scale_x_continuous(breaks = 1:3, labels = gl$genotype, expand = c(0, 0), limits = c(0.5, 3.5)) +
    scale_y_continuous(breaks = scales::breaks_pretty(4)) +
    coord_cartesian(ylim = c(lo - rng * 0.16, hi + rng * top), clip = "off") +
    labs(x = NULL, y = ylab) + theme(legend.position = "none", plot.title = element_text(size = 7, face = "bold.italic", hjust = 0))
}
lc2 <- function(p, title, xlab, ylab, corner, rs) {
  p <- as.data.table(p)[order(!is.na(r2), r2)]
  ggplot(p[is_lead == FALSE], aes(qtl_logp, gwas_logp)) +
    geom_point(aes(colour = ld_bin), size = 0.7, stroke = 0) +
    geom_point(data = p[is_lead == TRUE], shape = 21, fill = "white", colour = INK, size = 1.8, stroke = 0.55) +
    geom_label(data = p[is_lead == TRUE], label = rs, hjust = 0.85, vjust = -0.75, size = 4.8 / .pt, colour = INK, fill = "white", linewidth = 0, label.padding = unit(0.3, "mm")) +
    annotate("text", x = Inf, y = -Inf, label = corner, hjust = 1.05, vjust = -0.35, size = 5.2 / .pt, lineheight = 0.9) +
    scale_colour_manual(values = LOCUS_LD_COLS, limits = LOCUS_LD_LEVELS, guide = "none") +
    scale_x_continuous(breaks = scales::breaks_pretty(3), expand = expansion(mult = c(0.03, 0.06))) +
    scale_y_continuous(breaks = scales::breaks_pretty(3), expand = expansion(mult = c(0.03, 0.26))) +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(size = 7, face = "bold.italic", hjust = 0, margin = margin(0, 0, 0.4, 0, "mm")),
          axis.title = element_text(size = 6), aspect.ratio = 1)
}
SC  <- c(Prenatal = unname(STAGE_COLORS["Prenatal"]), Postnatal = unname(STAGE_COLORS["Postnatal"]))
EXP_COL <- "#2E5B7E"; APA_COL <- "#009E73"
STX <- c(Prenatal = unname(STAGE_TEXT["Prenatal"]), Postnatal = unname(STAGE_TEXT["Postnatal"]))
DATA <- file.path(ROOT, "data")
theme_set(theme_nature(base_size = 7) + theme(axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), plot.margin = margin(1, 1.5, 1, 1.5, "mm")))

# --- repaired coloc run and regional exports ---
CO <- fread(file.path(DATA, "07_coloc/coloc_run/all_colocalization.csv"), select = c("GWAS_ID", "QTL_ID", "Locus", "Phenotype", "PP4"))
CO[, gene := tstrsplit(Phenotype, "|", fixed = TRUE)[[2]]]
RT <- fread(file.path(DATA, "figures/fig6/fig6_reg_tables_all.tsv"))

# --- b: ATP2B1 x EUR ADHD (interaction coloc) ---
ADHD <- reg("ATP2B1", "EUR_ADHD", "interaction", "rs2279574")
lead_adhd <- "rs2279574"
b_tr <- list(
  gwas = mk_track(df_of(ADHD$points, "gwas_logp", lead_adhd), "EUR ADHD\n-log10(P)", sprintf("PP4 = %.2f", ADHD$pp4), TRUE),
  qtl  = mk_track(df_of(ADHD$points, "qtl_logp", lead_adhd), "G\u00d7Stage\napaQTL\n-log10(P)"))
b_tgt <- parse_target(ADHD$meta); b_tgt$chrom <- ADHD$meta$chrom
b_genes <- load_window_genes(ADHD$meta$chrom, ADHD$meta$min_pos, ADHD$meta$max_pos)
b_win <- c(ADHD$meta$min_pos, ADHD$meta$max_pos)
b_lead <- ADHD$points$pos[ADHD$points$snp == lead_adhd]
b_xlim <- b_win / 1e6; b_leadx <- b_lead / 1e6
b_gl <- gene_lanes(b_genes, b_tgt, b_lead, b_win, pw_mm = 98)
b_plots <- lapply(seq_along(b_tr), function(i) track_plot(b_tr[[i]], b_xlim, b_leadx, base_size = 6.5))
b_plots[[1]] <- tg(b_plots[[1]], "b")
b_plots[[1]] <- b_plots[[1]] + annotate("text", x = b_leadx - diff(b_xlim) * 0.015, y = max(b_tr[[1]]$df$logp), label = lead_adhd, hjust = 1, vjust = 0.3, size = 4.8 / .pt, colour = INK)
b_pg <- gene_plot(b_gl, b_tgt, b_xlim, base_size = 6.5)
b_hr <- c(rep(1, length(b_plots)), max(0.6, 0.17 * max(b_gl$lane)))
pb_locus <- wrap_plots(c(b_plots, list(b_pg)), ncol = 1, heights = b_hr)

# --- d/e: 12q24.31 ---
MT <- reg("MTRFR", "EAS_SCZ", "postnatal", "rs28362667")
OG <- reg("OGFOD2", "EUR_INT", "prenatal", "rs10773921")

# --- a: ATP2B1 genotype x stage boxes ---
BX <- fread(file.path(DATA, "figures/fig6/fig6_boxstats.tsv"))
cap_atp <- max(BX[gene == "ATP2B1" & layer == "PDUI", q3])   # PDUI ceiling = largest box q3 (9.964)
se <- fread(file.path(DATA, "07_coloc/fig6a_atp2b1_stage_effects.tsv"))
pa <- boxes("ATP2B1", "PDUI", "PDUI (log2 logit)", sprintf("beta[pre]==%.2f", se$beta_pre), sprintf("beta[post]==%.2f", se$beta_post),
            cap = cap_atp, ann_y = 0.32, top = 0.38) +
  labs(title = "ATP2B1", x = "rs2279574 genotype") + theme(axis.title.x = element_text(size = 6))

# --- c: PP4 strip at the ATP2B1 x EUR ADHD signal (G x Stage apaQTL vs cohort and 13 GTEx brain eQTLs) ---
EL <- rbind(fread(file.path(DATA, "figures/fig6/cohort_eqtl_gwas_coloc_long.tsv")), fread(file.path(DATA, "figures/fig6/gtex_eqtl_gwas_coloc_long.tsv")))
ea <- EL[signal_id == "EUR_ADHD__interaction__ATP2B1" & status == "ok"]
stopifnot(nrow(ea) == 15)
GRP <- c("G\u00d7Stage apaQTL", "Cohort eQTL", "GTEx brain eQTL")
cs <- rbind(data.table(grp = GRP[1], pp4 = ADHD$pp4, tissue = "apaQTL"),
            ea[source_family == "cohort", .(grp = GRP[2], pp4, tissue = source)],
            ea[source_family == "GTEx", .(grp = GRP[3], pp4, tissue = source)])
cs[, grp := factor(grp, rev(GRP))]
setorder(cs, grp, pp4)
cs[, yo := as.integer(grp) + (((seq_len(.N) - 1) %% 3) - 1) * 0.16 * (grp == GRP[3]), by = grp]
cs[, kind := fifelse(grp == GRP[1], "apa", "eqtl")]
lab_c <- cs[tissue == "Brain_Cerebellum"]
pc <- ggplot(cs, aes(pp4, yo)) +
  annotate("rect", xmin = 0.7, xmax = 1.02, ymin = -Inf, ymax = Inf, fill = "#F4F4F4") +
  geom_vline(xintercept = 0.7, linetype = "22", linewidth = 0.25, colour = MUTED) +
  geom_point(data = cs[kind == "eqtl"], shape = 21, fill = EXP_COL, colour = "white", stroke = 0.15, size = 1.35, alpha = 0.85) +
  geom_point(data = cs[kind == "apa"], shape = 21, fill = APA_COL, colour = "white", stroke = 0.2, size = 2.4) +
  geom_text(data = lab_c, aes(x = 0.86, label = "cerebellum"), hjust = 0.5, vjust = 1.9, size = 4.8 / .pt, colour = MUTED) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:3, labels = rev(GRP), limits = c(0.5, 3.6), expand = c(0, 0)) +
  labs(x = "Colocalization PP4", y = NULL, title = "EUR ADHD") +
  theme(plot.title = element_text(size = 6, hjust = 0, margin = margin(0, 0, 0.5, 0, "mm")), axis.text.y = element_text(size = 5.6))

# --- d-g: regional scatters and ARL3 boxes ---
pd <- lc2(MT$points, "MTRFR", expression("Postnatal apaQTL " * -log[10] * italic(P)), expression("EAS SCZ " * -log[10] * italic(P)), sprintf("PP4 = %.2f", MT$pp4), "rs28362667")
pe <- lc2(OG$points, "OGFOD2", expression("Prenatal apaQTL " * -log[10] * italic(P)), expression("EUR INT " * -log[10] * italic(P)), sprintf("PP4 = %.2f", OG$pp4), "rs10773921")
AJ <- reg("ARL3", "EAS_SCZ", "prenatal", "rs8354")
eur_max <- max(CO[gene == "ARL3" & GWAS_ID == "EUR_SCZ", PP4])
eur_bound <- ceiling(eur_max * 100) / 100
pfl <- lc2(AJ$points, "ARL3", expression("Prenatal apaQTL " * -log[10] * italic(P)), expression("EAS SCZ " * -log[10] * italic(P)),
           sprintf("PP4 = %.2f\nEUR SCZ \u2264 %.2f", AJ$pp4, eur_bound), "rs8354")

RG <- fread(file.path(DATA, "07_coloc/Fig6c_ARL3_eqtl_reg.tsv"))
eq8 <- RG[rsid == "rs8354"]; stopifnot(nrow(eq8) == 2)
b_eq_pre <- eq8[stage == "prenatal", beta]; b_eq_post <- eq8[stage == "postnatal", beta]
# apaQTL betas: rs8354 = chr10:102676884:C:T
E7 <- fread(file.path(DATA, "figures/ed7/ed7a_variant_gene_pairs.tsv"))[gene == "ARL3" & variant_id == "chr10:102676884:C:T"]
stopifnot(nrow(E7) == 1)
ph1 <- boxes("ARL3", "PDUI", "PDUI (log2 logit)", sprintf("beta[apaQTL]==%.2f", E7$slope_pre), sprintf("beta[apaQTL]==%.2f", E7$slope_post), ann_y = 0.27, top = 0.30) +
  labs(x = "rs8354 genotype") + theme(axis.title.x = element_text(size = 6))
ph2 <- boxes("ARL3", "Expression", "Expression (log2CPM)", sprintf("beta[eQTL]==%.2f", b_eq_pre), sprintf("beta[eQTL]==%.2f", b_eq_post)) +
  labs(x = "rs8354 genotype") + theme(axis.title.x = element_text(size = 6))

# --- shared keys ---
key_ld <- ggplot(data.table(l = factor(LOCUS_LD_LEVELS, levels = LOCUS_LD_LEVELS), x = seq_along(LOCUS_LD_LEVELS)), aes(x, 1)) +
  geom_point(aes(colour = l), size = 1.8, stroke = 0) +
  geom_text(aes(label = l), hjust = -0.35, size = 5.0 / .pt, colour = INK) +
  annotate("point", x = 7.15, y = 1, shape = 21, fill = "white", colour = INK, size = 2.0, stroke = 0.6) +
  annotate("text", x = 7.35, y = 1, label = "lead", hjust = 0, size = 5.0 / .pt) +
  scale_colour_manual(values = LOCUS_LD_COLS, limits = LOCUS_LD_LEVELS, guide = "none") +
  coord_cartesian(xlim = c(0.4, 8.1), clip = "off") + theme_void()
key_st <- ggplot(data.table(s = factor(c("Prenatal", "Postnatal"), c("Prenatal", "Postnatal")), x = c(1, 2.4)), aes(x, 1)) +
  geom_point(aes(fill = s, colour = s), shape = 22, size = 2.2) + geom_text(aes(label = s, colour = s), hjust = -0.35, size = 5.4 / .pt) +
  scale_fill_manual(values = SC, guide = "none") + scale_colour_manual(values = STX, guide = "none") +
  coord_cartesian(xlim = c(0.7, 3.6), clip = "off") + theme_void()

# --- assemble ---
ph1s <- ph1 + labs(title = "ARL3 rs8354")
row1 <- (tg(pa, "a") | pb_locus | tg(pc, "c")) + plot_layout(widths = c(44, 98, 38))
row2 <- (tg(pd, "d") | tg(pe, "e") | tg(pfl, "f") | tg(ph1s, "g") | ph2) + plot_layout(widths = c(36, 36, 36, 36, 36))
keys <- (wrap_elements(full = key_ld) | wrap_elements(full = key_st)) + plot_layout(widths = c(120, 60))
fig <- wrap_elements(full = row1) / wrap_elements(full = row2) / keys + plot_layout(heights = c(58, 48, 6))
ggsave(file.path(FIG_DIR, "fig6.png"), fig, width = 180, height = 116, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig6.pdf"), fig, width = 180, height = 116, units = "mm", device = cairo_pdf)
