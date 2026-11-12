# ED Fig 9: S-LDSC fold enrichment (a), MAGMA background robustness (b), AD trajectory gene sets (c).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg); library(ggrepel)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
DD <- file.path(ROOT, "data/06_disease"); DF <- file.path(ROOT, "data/figures/ed9")
theme_set(theme_nature(base_size = 6.5) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 6.5, colour = INK),
  legend.text = element_text(size = 6, colour = INK), legend.title = element_text(size = 6, colour = INK),
  legend.key.size = unit(2.2, "mm"), plot.margin = margin(1, 1.5, 1, 1, "mm"),
  plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1)))
T6 <- 6 / .pt
tg <- function(p, t) p + labs(tag = t)

# --- a: S-LDSC fold enrichment, GWAS with h2 z >= 7 ---
pr  <- fread(file.path(DD, "sldsc_enrich_processed.tsv"))
fl  <- fread(file.path(DF, "sldsc_enrich_sum_flagged.tsv"))
ad2 <- fread(file.path(DD, "sldsc_ad2_summary.tsv"))
ANN <- c(prenatal_apa_qtl = "Prenatal", postnatal_apa_qtl = "Postnatal")
pr[, gw := paste(Population, Trait)]; fl[, gw := paste(Population, Trait)]
h2z <- unique(fl[, .(gw, h2_z)])
d1 <- pr[Annotation %in% names(ANN), .(gw, pop = Population, trait = Trait, Annotation, enr = Enrichment, p = Enrichment_p)]
d1[, h2_z := h2z$h2_z[match(gw, h2z$gw)]]
d1 <- d1[h2_z >= 7 & !(gw == "EUR AD")]          # Kunkle AD has h2 z < 7; Bellenguez row added below
ad2f <- ad2[apoe_excluded == "F" & annotation %in% names(ANN)]
ad2_z <- ad2[gwas_id == "EUR_AD2", h2 / h2_se][1]; stopifnot(ad2_z >= 7)
d2 <- data.table(gw = "EUR AD2", pop = "EUR", trait = "AD", Annotation = ad2f$annotation, enr = ad2f$enrichment,
                 p = ad2f$enrichment_p, h2_z = ad2_z)
d <- rbind(d1, d2)
d[, fdr := p.adjust(p, "BH"), by = gw]            # BH across the two stages within each GWAS
stopifnot(all(d[, .N, by = gw]$N == 2))
d[, stage := factor(ANN[Annotation], c("Prenatal", "Postnatal"))]
ORD <- c("SCZ", "BIP", "MDD", "ADHD", "ASD", "PTSD", "INT", "PD", "AD")
d[, row := fifelse(pop == "EAS", sprintf("%s (EAS)", trait), trait)]
d[, row := factor(row, rev(unique(d[order(pop == "EAS", match(trait, ORD)), row])))]
d[, l2 := log2(enr)]
LIM <- 1.1; stopifnot(max(abs(d$l2)) <= LIM)
pa <- ggplot(d, aes(stage, row, fill = l2)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = paste0(sprintf("%.2f", enr), ifelse(fdr < 0.05, "*", ""))), size = T6, colour = INK) +
  scale_fill_gradient2(low = DIVERGING_RAMP[1], mid = DIVERGING_RAMP[2], high = DIVERGING_RAMP[3], midpoint = 0,
                       limits = c(-LIM, LIM), breaks = log2(c(0.5, 1, 2)), labels = c("0.5", "1", "2"),
                       name = "Fold enrichment",
                       guide = guide_colourbar(title.position = "left", title.vjust = 1,
                                               barwidth = unit(16, "mm"), barheight = unit(1.6, "mm"))) +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        legend.margin = margin(t = 0), legend.box.spacing = unit(0.5, "mm"), axis.text.x.top = element_text(size = 6.5, face = "bold",
        colour = c(STAGE_TEXT[["Prenatal"]], STAGE_TEXT[["Postnatal"]])))
pa <- tg(pa, "a")

# --- b: MAGMA gene-set P, genome vs APA-tested background ---
g0 <- fread(file.path(DD, "magma_gsa_summary.tsv"))[set_name %in% c("prenatal", "postnatal")]
g1 <- fread(file.path(DD, "magma_condition_hide_sum.tsv"))[set_name %in% c("prenatal", "postnatal")]
mb <- merge(g0[, .(gwas_id, set_name, p_all = p)], g1[, .(gwas_id, set_name, p_tst = p)], by = c("gwas_id", "set_name"))
mb <- mb[!gwas_id %in% c("EUR_BIP", "EAS_BIP")]
mb[, key := fifelse(gwas_id == "EUR_SCZ", paste("SCZ", set_name), "other")]
ml <- melt(mb, id.vars = c("gwas_id", "set_name", "key"), measure.vars = c("p_all", "p_tst"), variable.name = "bg", value.name = "p")
ml[, x := fifelse(bg == "p_all", 1, 2)]; ml[, lp := -log10(p)]
ml[, col := fifelse(set_name == "prenatal", STAGE_COLORS[["Prenatal"]], STAGE_COLORS[["Postnatal"]])]
hl <- ml[key != "other"]
lab <- hl[x == 1][, txt := fifelse(set_name == "prenatal", "SCZ prenatal", "SCZ postnatal")]
pb <- ggplot(ml, aes(x, lp, group = interaction(gwas_id, set_name))) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = MUTED) +
  geom_line(data = ml[key == "other"], aes(colour = I(col)), linewidth = 0.3, alpha = 0.45) +
  geom_point(data = ml[key == "other"], aes(colour = I(col)), size = 0.8, alpha = 0.6) +
  geom_line(data = hl, aes(colour = I(col)), linewidth = 0.9) +
  geom_point(data = hl, aes(fill = I(col)), shape = 21, colour = "white", stroke = 0.3, size = 2.1) +
  geom_text(data = lab, aes(x = 0.93, y = lp, label = txt, colour = I(fifelse(set_name == "prenatal", STAGE_TEXT[["Prenatal"]], STAGE_TEXT[["Postnatal"]]))),
            hjust = 1, size = T6, fontface = "bold") +
  annotate("text", x = 1.5, y = -log10(0.05), label = "P = 0.05", hjust = 0.5, vjust = -0.4, size = T6, colour = MUTED) +
  scale_x_continuous(breaks = 1:2, labels = c("Genome", "APA-tested"), limits = c(0.12, 2.12), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Background genes", y = expression("Gene-set"~-log[10]~italic(P))) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
pb <- tg(pb, "b")

# --- c: AD trajectory gene sets, both backgrounds, with and without APOE ---
tr <- fread(file.path(DF, "magma_traj_sum_EUR_AD2.tsv"))
tr[, set := factor(c(traj_C1_specific = "C1", traj_C2_specific = "C2", traj_Shared = "Bidirectional")[set_name], c("C1", "C2", "Bidirectional"))]
tr[, rowl := paste0(set, fifelse(apoe_excluded == "T", " \u2212APOE", ""))]
tr[, rowl := factor(rowl, rev(c("Bidirectional", "Bidirectional \u2212APOE", "C2", "C2 \u2212APOE", "C1", "C1 \u2212APOE")))]
tr[, bg := factor(fifelse(background == "default", "Genome", "APA-tested"), c("Genome", "APA-tested"))]
tr[, lp := -log10(p)]
stopifnot(nrow(tr) == 12, !anyNA(tr$rowl))
wc <- dcast(tr, rowl + set ~ bg, value.var = "lp")
ccol <- setNames(CLUSTER_COLORS[c("C1", "C2", "Shared")], c("C1", "C2", "Bidirectional"))
pc <- ggplot(tr, aes(lp, rowl)) +
  geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = MUTED) +
  geom_segment(data = wc, aes(x = `Genome`, xend = `APA-tested`, y = rowl, yend = rowl, colour = set),
               linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(aes(colour = set, fill = ifelse(bg == "Genome", ccol[as.character(set)], "white")),
             shape = 21, size = 2, stroke = 0.6) +
  geom_text(data = tr[bg == "Genome" & fdr_within_gwas < 0.05], aes(label = sprintf("FDR %.3f", fdr_within_gwas)),
            nudge_y = 0.42, size = T6, colour = INK) +
  annotate("point", x = c(1.30, 1.30), y = c(2.0, 1.25), shape = 21, size = 2, stroke = 0.6,
           colour = INK, fill = c(INK, "white")) +
  annotate("text", x = 1.38, y = c(2.0, 1.25), label = c("Genome", "APA-tested"), hjust = 0, size = T6) +
  scale_colour_manual(values = ccol, guide = "none") + scale_fill_identity() +
  scale_x_continuous(limits = c(0, 2.2), breaks = 0:2, expand = c(0, 0.04)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.8))) +
  labs(x = expression("Gene-set"~-log[10]~italic(P)*", AD GWAS"), y = NULL) +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())
pc <- tg(pc, "c")

fig <- pa + pb + pc + plot_layout(widths = c(40, 76, 52))
ggsave(file.path(FIG_DIR, "ed9.png"), fig, width = 180, height = 66, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed9.pdf"), fig, width = 180, height = 66, units = "mm", device = cairo_pdf)
