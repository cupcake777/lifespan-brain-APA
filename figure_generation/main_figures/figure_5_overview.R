# Fig 5 a-c: S-LDSC coefficient z-scatter (a), MAGMA gene-set lollipops (b), apaQTL-GWAS coloc Manhattan (c).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(ggh4x)})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
DATA <- file.path(ROOT, "data")

theme_set(theme_nature(base_size = 7) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3), axis.ticks.length = unit(0.8, "mm"),
  plot.margin = margin(1, 1.5, 1, 1.5, "mm")))
T5 <- 5.5 / .pt; T6 <- 6 / .pt
SC  <- c(prenatal = unname(STAGE_COLORS["Prenatal"]), postnatal = unname(STAGE_COLORS["Postnatal"]))
STX <- c(prenatal = unname(STAGE_TEXT["Prenatal"]), postnatal = unname(STAGE_TEXT["Postnatal"]))
APA_COL <- "#009E73"; EXP_COL <- "#2E5B7E"
TRAIT_ORDER <- c("SCZ","BIP","MDD","ADHD","ASD","PTSD","ANX","INS","INT","AD","PD","EP","ALS","AIS")
ZT <- qnorm(0.95)

# --- a: joint S-LDSC coefficient z, apaQTL vs eQTL ---
tau <- fread(file.path(DATA, "04_qtl/03_sldsc_tau_star_handoff.tsv"))
a <- dcast(tau, population + trait + stage ~ annotation, value.var = c("coef_z", "coef_fdr"))
setnames(a, c("coef_z_apaQTL","coef_z_eQTL","coef_fdr_apaQTL","coef_fdr_eQTL"), c("za","ze","fa","fe"))
stopifnot(nrow(a) == 38)
a[, short := sub("^(EUR|EAS)_", "", trait)]
a[, short_d := sub("^BIP$", "BD", short)]
a[, sig := fifelse(za > ZT & ze > ZT, "both", fifelse(za > ZT, "apa", fifelse(ze > ZT, "eqtl", "none")))]
a[, lab := fifelse(sig != "none", sprintf("%s %s", population, short_d), NA_character_)]
lim <- c(min(-1.6, floor(min(a$za, a$ze) * 2) / 2), ceiling(max(a$za, a$ze) * 2) / 2 + 0.2)
# heritability filter (h2 z >= 7); lim is computed from all 38 rows so the axes do not change
h2z <- unique(tau[, .(population, trait, stage, h2z = h2 / h2_se)])
a <- merge(a, h2z, by = c("population", "trait", "stage"), all.x = TRUE)
stopifnot(!anyNA(a$h2z))
a <- a[h2z >= 7]
sigtr <- a[sig != "none", unique(trait)]
a[, issig := sig != "none"]
pa <- ggplot(a, aes(ze, za)) +
  annotate("rect", xmin = lim[1], xmax = ZT, ymin = ZT, ymax = lim[2], fill = APA_COL, alpha = 0.12) +
  annotate("rect", xmin = ZT, xmax = lim[2], ymin = ZT, ymax = lim[2], fill = "#8C8C8C", alpha = 0.10) +
  annotate("rect", xmin = ZT, xmax = lim[2], ymin = lim[1], ymax = ZT, fill = EXP_COL, alpha = 0.10) +
  annotate("text", x = lim[1] + 0.08, y = lim[2] - 0.08, label = "apaQTL only", hjust = 0, vjust = 1, size = T6, colour = APA_COL, fontface = "bold") +
  annotate("text", x = lim[2] - 0.06, y = ZT + 0.06, label = "Both", hjust = 1, vjust = 0, size = T6, colour = "#5A5A5A", fontface = "bold") +
  annotate("text", x = lim[2] - 0.06, y = lim[1] + 0.08, label = "eQTL only", hjust = 1, vjust = 0, size = T6, colour = EXP_COL, fontface = "bold") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = NEUTRAL, linetype = "dotted") +
  geom_hline(yintercept = ZT, linewidth = 0.25, linetype = "dashed", colour = MUTED) +
  geom_vline(xintercept = ZT, linewidth = 0.25, linetype = "dashed", colour = MUTED) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = GRID_COL) + geom_vline(xintercept = 0, linewidth = 0.2, colour = GRID_COL) +
  geom_point(data = a[issig == FALSE], aes(fill = stage, shape = population), colour = "white", stroke = 0.1, size = 1.3, alpha = 0.5) +
  geom_point(data = a[issig == TRUE], aes(fill = stage, shape = population), colour = "white", stroke = 0.3, size = 2.2) +
  geom_text_repel(aes(label = lab, colour = stage), size = T5, min.segment.length = 0.1, segment.size = 0.2, segment.colour = MUTED,
                  box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf, seed = 11, na.rm = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = SC, labels = c(prenatal = "Prenatal", postnatal = "Postnatal"), name = NULL, breaks = c("prenatal", "postnatal")) +
  scale_colour_manual(values = STX, guide = "none") +
  scale_shape_manual(values = c(EUR = 21, EAS = 24), name = NULL) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(x = "eQTL coefficient z (joint S-LDSC)", y = "apaQTL coefficient z (joint S-LDSC)") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2, colour = NA), order = 1), shape = guide_legend(override.aes = list(fill = NEUTRAL, colour = NA, size = 1.8), order = 2)) +
  theme(legend.position = "top", legend.justification = c(0, 1), legend.box = "horizontal", legend.box.just = "left",
        legend.spacing.x = unit(0.5, "mm"), legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.5), legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, -1, 0, "mm"))

# --- b: MAGMA gene-set P by GWAS and stage ---
m <- fread(file.path(DATA, "06_disease/magma_gsa_summary.tsv"))[set_name %in% c("prenatal", "postnatal")]
upd <- fread(file.path(DATA, "figures/fig5/S_MAGMA.tsv"))[gwas_id %in% c("EUR_BIP", "EAS_BIP") & set_name %in% c("prenatal", "postnatal") & model == "unconditioned",
  .(gwas_id, set_name, p_new = p, fdr_new = fdr_within_gwas)]
stopifnot(nrow(upd) == 4)
m <- merge(m, upd, by = c("gwas_id", "set_name"), all.x = TRUE)
m[!is.na(p_new), `:=`(p = p_new, fdr_within_gwas = fdr_new)][, c("p_new", "fdr_new") := NULL]
m[, short := factor(sub("^(EUR|EAS)_", "", gwas_id), TRAIT_ORDER)]
stopifnot(!anyNA(m$short))
m[, stage := factor(fifelse(set_name == "prenatal", "Prenatal", "Postnatal"), c("Prenatal", "Postnatal"))]
m[, population := factor(population, c("EUR", "EAS"))]
m[, lp := -log10(p)]
m[, short := factor(as.character(short), rev(TRAIT_ORDER))]
m[, short_d := sub("^BIP$", "BD", as.character(short))]
m[, fdrsig := fdr_within_gwas < 0.05]
m[, stg := fifelse(stage == "Prenatal", "prenatal", "postnatal")]
GAP <- 0.42; XM <- 3.1
brk <- c(0, 1, 2, 3)
m[, sgn := fifelse(stg == "prenatal", -1, 1)]
m[, x0 := sgn * GAP][, x1 := sgn * (GAP + lp)]
rowband <- unique(m[, .(population, short)])[order(population, -as.integer(short))][, k := seq_len(.N), by = population][, bf := fifelse(k %% 2 == 0, "#F3F3F3", "white")]
lab_b <- m[, .(anysig = any(fdrsig), short_d = short_d[1]), by = .(population, short)]
pb <- ggplot(m, aes(y = short)) +
  geom_tile(data = rowband, aes(x = 0, y = short, fill = I(bf)), width = 2 * (GAP + XM), height = 1, inherit.aes = FALSE) +
  geom_vline(xintercept = c(-1, 1) * (GAP + -log10(0.05)), linewidth = 0.25, linetype = "dashed", colour = MUTED) +
  geom_vline(xintercept = c(-1, 1) * GAP, linewidth = 0.3, colour = "#BDBDBD") +
  geom_segment(aes(x = x0, xend = x1, yend = short, colour = stg), linewidth = 0.5) +
  geom_point(aes(x = x1, fill = stg, alpha = p < 0.05), shape = 21, colour = "white", stroke = 0.25, size = 2.1) +
  geom_text(data = m[fdrsig == TRUE], aes(x = x1 + sgn * 0.16, label = "*"), size = 8 / .pt, vjust = 0.75, colour = INK) +
  geom_text(data = lab_b, aes(x = 0, y = short, label = short_d, fontface = ifelse(anysig, "bold", "plain")), size = T6, colour = INK, inherit.aes = FALSE) +
  facet_grid(population ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_text(data = data.table(population = factor("EUR", c("EUR", "EAS")), x = c(-1, 1) * (GAP + XM / 2), lab = c("Prenatal", "Postnatal"), col = unname(STX[c("prenatal", "postnatal")])),
            aes(x = x, y = Inf, label = lab, colour = I(col)), vjust = -0.6, size = 6.5 / .pt, fontface = "bold", inherit.aes = FALSE) +
  scale_colour_manual(values = SC, guide = "none") + scale_fill_manual(values = SC, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.45, `TRUE` = 1), guide = "none") +
  scale_x_continuous(limits = c(-(GAP + XM), GAP + XM), breaks = c(-(GAP + rev(brk)), GAP + brk), labels = c(rev(brk), brk), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = expression("MAGMA gene-set " * -log[10] * italic(P)), y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        panel.spacing.y = unit(1.2, "mm"), plot.margin = margin(4, 1.5, 1, 1.5, "mm"),
        strip.text.y.left = element_text(size = 6.5, face = "bold", colour = INK, angle = 90, margin = margin(0, 0.8, 0, 0.6, "mm")),
        strip.background.y = element_rect(fill = "#E6E6E6", colour = NA), strip.placement = "outside")

# --- c: apaQTL-GWAS coloc Manhattan, one point per GWAS x stage x APA event ---
co <- fread(file.path(DATA, "07_coloc/coloc_run/all_colocalization.csv"), select = c("GWAS_ID", "QTL_ID", "Locus", "Phenotype", "PP4"))
co[, gene := tstrsplit(Phenotype, "|", fixed = TRUE)[[2]]]
co[, coord := tstrsplit(Phenotype, "|", fixed = TRUE)[[3]]]
co[, chr := sub("^chr", "", sub(":.*", "", coord))]
co[, pos := as.numeric(sub("-.*", "", sub("^[^:]*:", "", coord)))]
co <- co[chr %in% as.character(1:22)]
# interaction-term tests only for the 193 genes with a significant genotype-by-stage interaction (Fig 4 definition)
G193 <- fread(file.path(DATA, "05_gxs/intx_hc3_core_groups.tsv"))$gene_id
stopifnot(length(G193) == 193)
co <- co[QTL_ID != "interaction" | gene %in% G193]
cc <- co[, .(PP4 = max(PP4)), by = .(GWAS_ID, QTL_ID, Phenotype, gene, chr, pos)]
cls <- fread(file.path(DATA, "figures/fig5/eqtl_gwas_signal_sum.tsv"))
sc <- unique(cls[, .(gene, GWAS = gwas_id, stage = apa_stage,
                     cl = fifelse(is.na(class) | class == "", "eQTL not testable", fifelse(class == "Shared", "eQTL also colocalizes", "No eQTL colocalization")))])
stopifnot(!anyDuplicated(sc[, .(gene, GWAS, stage)]))
cc <- merge(cc, sc, by.x = c("gene", "GWAS_ID", "QTL_ID"), by.y = c("gene", "GWAS", "stage"), all.x = TRUE)
cc[, hit := PP4 >= 0.7]
# GRCh38 chromosome lengths for a true genome axis
GL <- c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,
        114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468)
off <- data.table(chr = as.character(1:22), len = GL, start = c(0, cumsum(GL)[-22]))
off[, mid := start + len / 2]
cc <- merge(cc, off[, .(chr, start)], by = "chr"); cc[, gx := start + pos]
CL_COL <- c("No eQTL colocalization" = APA_COL, "eQTL also colocalizes" = EXP_COL, "eQTL not testable" = NEUTRAL)
STG_SHP <- c(prenatal = 21, postnatal = 22, interaction = 23)
key <- c("SCAF1", "KLC1", "ATP2B1", "MTRFR", "OGFOD2", "ARL3")
hits <- cc[hit == TRUE]
top_by_gene <- hits[order(-PP4)][, .SD[1], by = gene]
lab <- top_by_gene[gene %in% key]; lab[, keyg := TRUE]
pc <- ggplot() +
  geom_rect(data = off[as.integer(chr) %% 2 == 0], aes(xmin = start, xmax = start + len, ymin = -Inf, ymax = Inf), fill = "#F5F5F5", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.7, linewidth = 0.25, linetype = "dashed", colour = MUTED) +
  geom_point(data = cc[hit == FALSE], aes(gx, PP4), colour = "#C4C4C4", size = 0.35, alpha = 0.45, shape = 16) +
  geom_point(data = hits, aes(gx, PP4, fill = cl, shape = QTL_ID), size = 1.7, colour = "white", stroke = 0.25) +
  geom_text_repel(data = lab, aes(gx, PP4, label = gene, fontface = ifelse(keyg, "bold.italic", "italic")), size = T5, colour = INK,
                  nudge_y = 1.10 - lab$PP4, direction = "x", hjust = 0.5, min.segment.length = 0, segment.size = 0.25, segment.colour = MUTED,
                  box.padding = 0.25, force = 2, max.overlaps = Inf, seed = 7, ylim = c(1.06, 1.14)) +
  scale_colour_manual(values = c(`0` = "#BDBDBD", `1` = "#9A9A9A"), guide = "none") +
  scale_fill_manual(values = CL_COL, name = NULL, breaks = c("eQTL also colocalizes", "No eQTL colocalization", "eQTL not testable"),
                    labels = c("No eQTL colocalization" = "apaQTL-only", "eQTL also colocalizes" = "Shared", "eQTL not testable" = "NA")) +
  scale_shape_manual(values = STG_SHP, breaks = c("prenatal", "postnatal", "interaction"), labels = c(prenatal = "Prenatal", postnatal = "Postnatal", interaction = "G\u00D7Stage"), name = NULL) +
  scale_x_continuous(breaks = off$mid, labels = c(1:16, "", 18, "", 20, "", 22), expand = c(0.005, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1.16), expand = c(0, 0.01)) +
  labs(x = "Chromosome", y = "apaQTL\u2013GWAS PP4") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), order = 1), shape = guide_legend(override.aes = list(fill = NEUTRAL, size = 1.8), order = 2)) +
  theme(legend.position = "inside", legend.position.inside = c(0.995, 0.36), legend.justification = c(1, 0.5), legend.box = "horizontal", legend.box.just = "top",
        legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
        legend.key.size = unit(2.6, "mm"),
        legend.text = element_text(size = 5.8), legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(0, "mm"), legend.margin = margin(0, 5, 0, 0, "mm"),
        legend.box.margin = margin(1, 1, 1, 1, "mm"), axis.ticks.x = element_blank())
