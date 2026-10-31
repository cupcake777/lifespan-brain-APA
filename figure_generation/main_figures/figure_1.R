# Fig 1 composite: cohort (a), lifespan landscape (b), trajectories (c), gene-set overlap (d), enrichment (e), LOEUF (f).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrastr); library(ragg); library(cowplot)
})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
D2 <- file.path(ROOT, "data/02_traj"); D3 <- file.path(ROOT, "data/03_omics"); D4 <- file.path(ROOT, "data/04_qtl")
BS <- 7

# --- Fig 1 palette (overrides INK and MUTED from theme.R for this figure) ---
C1_COLOR     <- "#4A7BA7"
C2_COLOR     <- "#D17A4A"
SHARED_COLOR <- "#7A9E9F"
OTHER_COLOR  <- "#C8C8C8"
INK          <- "#2A2A2A"
MUTED        <- "#6B7280"
LINE_COL     <- "#E0E0E0"
CLS  <- c("C1-specific" = C1_COLOR, "Shared" = SHARED_COLOR, "C2-specific" = C2_COLOR)
TINT <- c("C1-specific" = scales::alpha(C1_COLOR, 0.22),
          "Shared"      = scales::alpha(SHARED_COLOR, 0.28),
          "C2-specific" = scales::alpha(C2_COLOR, 0.22))
SEXC <- c(Male = "#3E434A", Female = "#A8ADB5")
OMC  <- c(Transcriptome = C1_COLOR, Genotype = SHARED_COLOR, Proteome = C2_COLOR)

theme_set(theme_nature(base_size = BS) + theme(
  text             = element_text(family = FONT, colour = INK),
  axis.text        = element_text(size = 6, colour = INK),
  axis.title       = element_text(size = 6.5, colour = INK),
  axis.line        = element_line(linewidth = 0.35, colour = INK),
  axis.ticks       = element_line(linewidth = 0.35, colour = INK),
  panel.background = element_rect(fill = "#FFFFFF", colour = NA),
  plot.background  = element_rect(fill = "#FFFFFF", colour = NA),
  plot.margin      = margin(1, 1, 1, 1, "mm")
))

T6 <- 6 / .pt; T55 <- 5.5 / .pt
MINUS <- "\u2212"
fmt_minus <- function(x) sub("-", MINUS, trimws(format(x)), fixed = TRUE)

# --- a: cohort ages, sex and omics availability ---
meta <- fread(file.path(ROOT, "raw/all.912.meta"))
# genotyped QTL individuals (189 + 508 = 697)
geno <- unique(trimws(unlist(lapply(file.path(D4, c("apa.sample.Prenatal.final", "apa.sample.Postnatal.final")), readLines))))
stopifnot(length(geno) == 697L)
prot <- fread(file.path(D3, "matched_sample_manifest_raw.tsv"), select = "rna_sample_id")[[1]]
meta[, `:=`(sex_l = factor(ifelse(sex == 0, "Male", "Female"), c("Male", "Female")),
            stage = ifelse(age < 0, "Prenatal", "Postnatal"),
            gw = 40 + age * 365 / 7)]
meta[, `:=`(has_rna = TRUE, has_dna = as.character(sample) %in% geno, has_prot = as.character(sample) %in% prot)]

mk_stage <- function(d, xvar, step, lims, brks, xlab, show_y) {
  d <- copy(d); d[, xv := get(xvar)]
  d[, bin := floor((xv - lims[1]) / step)]; d[, xc := lims[1] + (bin + 0.5) * step]
  bs <- d[, .N, by = .(xc, sex_l)]; bt <- d[, .N, by = xc]
  pb <- ggplot() +
    geom_col(data = bs, aes(xc, N, fill = sex_l), width = step * 0.86, colour = "white", linewidth = 0.2, alpha = 0.92) +
    geom_text(data = bt, aes(xc, N, label = N), vjust = -0.3, size = 5 / .pt, colour = INK) +
    scale_fill_manual(values = SEXC, name = NULL,
                      labels = c(Male = sprintf("Male (%d)", sum(meta$sex_l == "Male")),
                                 Female = sprintf("Female (%d)", sum(meta$sex_l == "Female")))) +
    scale_x_continuous(limits = lims, breaks = brks, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
    labs(x = NULL, y = if (show_y) "Samples" else NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
          legend.position = "none", plot.margin = margin(1, 1, 0, 1, "mm"))
  r <- rbind(d[has_rna == TRUE, .(xv, om = "Transcriptome")],
             d[has_dna == TRUE, .(xv, om = "Genotype")],
             d[has_prot == TRUE, .(xv, om = "Proteome")])
  r[, om := factor(om, rev(names(OMC)))]
  pr <- ggplot(r, aes(x = xv, y = om, colour = om)) +
    geom_point(shape = "|", size = 2.1, alpha = 0.85) +
    scale_colour_manual(values = OMC, guide = "none") +
    scale_x_continuous(limits = lims, breaks = brks, expand = c(0, 0)) +
    labs(x = xlab, y = NULL) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = if (show_y) element_text(size = 5.5, colour = INK) else element_blank(),
          axis.text.x = element_text(size = 5.5),
          plot.margin = margin(0, 1, 1, 1, "mm"))
  list(pb, pr)
}

pre  <- mk_stage(meta[stage == "Prenatal"],  "gw",  2,  c(12, 40), seq(14, 38, 8), "Gestational week", TRUE)
pre[[1]]  <- pre[[1]]  + annotate("text", x = 12.6, y = Inf, vjust = 1.1, hjust = 0, size = T6,
                                  fontface = "bold", colour = INK,
                                  label = sprintf("Prenatal, n = %d", sum(meta$stage == "Prenatal")))
post <- mk_stage(meta[stage == "Postnatal"], "age", 10, c(0, 110), seq(0, 100, 25), "Age (years)", FALSE)
post[[1]] <- post[[1]] + annotate("text", x = 109, y = Inf, vjust = 1.1, hjust = 1, size = T6,
                                  fontface = "bold", colour = INK,
                                  label = sprintf("Postnatal, n = %d", sum(meta$stage == "Postnatal")))
sexleg <- get_legend(
  ggplot(meta, aes(age, fill = sex_l)) + geom_bar() +
    scale_fill_manual(values = SEXC, name = NULL,
                      labels = c(Male = sprintf("Male (%d)", sum(meta$sex_l == "Male")),
                                 Female = sprintf("Female (%d)", sum(meta$sex_l == "Female")))) +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 5.5))
)
pre[[1]] <- pre[[1]] + labs(tag = "a") + theme(plot.tag = element_text(size = 8, face = "bold"))
pa <- wrap_plots(pre[[1]], post[[1]], pre[[2]], post[[2]], ncol = 2, widths = c(1, 1.45), heights = c(2.2, 1))
pa <- wrap_plots(pa, wrap_elements(full = sexleg), ncol = 1, heights = c(1, 0.06))

# --- b: lifespan landscape, -log10 q by signed change in PDUI ---
s <- fread(file.path(D2, "02_sites_clustered.tsv"), select = c("age_qvalue", "delta_logit", "cluster_id"))
s <- s[is.finite(age_qvalue) & is.finite(delta_logit)]
# q = 0 (underflow) is floored at the smallest non-zero q; no cap on the y-scale
Q_MIN_POS <- min(s$age_qvalue[s$age_qvalue > 0])
s[, y := -log10(pmax(age_qvalue, Q_MIN_POS))]
YC <- max(s$y)   # label anchor = top of the data
s[, grp := factor(fifelse(cluster_id %in% 1, "C1", fifelse(cluster_id %in% 2, "C2", "Other")), c("Other", "C1", "C2"))]
nn <- s[, .N, by = grp]; lab <- function(g) format(nn[grp == g, N], big.mark = ",")
pb <- ggplot(mapping = aes(delta_logit, y)) +
  rasterise(geom_point(data = s[grp == "Other"], colour = OTHER_COLOR, size = 0.2, alpha = 0.5, shape = 16), dpi = 600) +
  rasterise(geom_point(data = s[grp != "Other"], aes(colour = grp), size = 0.25, alpha = 0.7, shape = 16), dpi = 600) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = "dashed", colour = LINE_COL) +
  annotate("text", x = 10.4, y = YC, label = paste0("C1 distal\n", lab("C1")), colour = C1_COLOR,
           size = T6, fontface = "bold", hjust = 1, vjust = 1, lineheight = 0.9) +
  annotate("text", x = -10.4, y = YC, label = paste0("C2 proximal\n", lab("C2")), colour = C2_COLOR,
           size = T6, fontface = "bold", hjust = 0, vjust = 1, lineheight = 0.9) +
  annotate("text", x = 0, y = YC, label = paste0("Non-directional\n", lab("Other")), colour = MUTED,
           size = T6, vjust = 1, lineheight = 0.9) +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5, base = 2), limits = c(-11, 11),
                     breaks = c(-8, -2, 0, 2, 8), labels = fmt_minus, expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, YC * 1.03), breaks = seq(0, 300, 100), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Lifespan change, \u0394logit(PDUI)", y = expression(-log[10]~italic(q)[age]))

# --- c: lifespan trajectories of directional sites (line panel + heatmap) ---
sites <- fread(file.path(D2, "02_sites_clustered.tsv"), select = c("cluster_id", "directional_pass"))
sup   <- fread(file.path(D2, "01_support.tsv"))
tz    <- fread(file.path(D2, "traj_z.csv"))
lg    <- fread(file.path(D2, "traj_logit.csv"))

dir_idx <- which(sites$directional_pass)
gcols   <- sprintf("grid_%03d", 0:127)
M_z     <- as.matrix(tz[dir_idx, ..gcols])          # all 128 grids, incl. interpolation span
M_lg    <- as.matrix(lg[dir_idx, ..gcols])
age_all <- sup$grid_age
supp    <- sup$is_supported
cl      <- sites$cluster_id[dir_idx]

# piecewise-linear age -> x warp: prenatal [-0.5, 0] gets PRE_FRAC of the width
PRE_FRAC <- 0.26
x_of_age <- function(a) {
  a <- pmin(pmax(a, min(age_all)), max(age_all))
  ifelse(a <= 0,
         PRE_FRAC * (a - min(age_all)) / (0 - min(age_all)),
         PRE_FRAC + (1 - PRE_FRAC) * a / max(age_all))
}
xa    <- x_of_age(age_all)
GAPR  <- 0.012
n_pre <- sum(cl == 1); n_post <- sum(cl == 2)
h1 <- n_pre / nrow(M_z); h2 <- n_post / nrow(M_z)
yb <- list(c(1 - h1 * (1 - GAPR), 1), c(0, h2 * (1 - GAPR)))

# row order: peak-change timing, then total magnitude
dM     <- t(apply(M_lg, 1, diff))
peak_g <- apply(abs(dM), 1, which.max)
i_pre  <- which(age_all <= 0)
fr <- (M_lg[, max(i_pre)] - M_lg[, 1]) / (M_lg[, ncol(M_lg)] - M_lg[, 1] + 1e-6)
ord <- lapply(1:2, function(k) {
  ii <- which(cl == k)
  ii[order(peak_g[ii], fr[ii], decreasing = (k == 2))]
})

# quantile-clipped continuous bipolar palette (data-driven endpoints)
Q_LO <- as.numeric(quantile(M_z, 0.02, na.rm = TRUE))
Q_HI <- as.numeric(quantile(M_z, 0.98, na.rm = TRUE))
HEAT_CONTINUOUS <- grDevices::colorRampPalette(c(
  "#2E5B7E", "#4A7BA7", "#7FA3C4", "#B8CEDF", "#E2E8EA",
  "#EDEBE7",
  "#F0E0D5", "#E5BEA4", "#DA9C77", "#D17A4A", "#A9501F"
))(255)
toras_cont <- function(A) {
  A <- pmax(pmin(A, Q_HI), Q_LO)
  idx <- round((A - Q_LO) / (Q_HI - Q_LO) * 254) + 1
  matrix(HEAT_CONTINUOUS[idx], nrow = nrow(A))
}

# solid = directly supported grid ages; dashed = interpolated span not used in modelling
dsp <- seq(1, ncol(M_z), by = 3)
supp_dsp <- supp[dsp]
pts_solid <- dsp[supp_dsp]
pts_uns   <- dsp[!supp_dsp]
if (length(pts_uns) && min(pts_uns) > 1) pts_uns <- c(pts_solid[max(which(pts_solid < min(pts_uns)))], pts_uns)
L <- rbind(
  data.table(cl = rep(ifelse(cl == 1, "C1", "C2"), times = length(pts_solid)),
             x = rep(xa[pts_solid], each = nrow(M_z)), interp = FALSE,
             z = as.vector(M_z[, pts_solid, drop = FALSE])),
  data.table(cl = rep(ifelse(cl == 1, "C1", "C2"), times = length(pts_uns)),
             x = rep(xa[pts_uns], each = nrow(M_z)), interp = TRUE,
             z = as.vector(M_z[, pts_uns, drop = FALSE])))
S <- L[, .(med = median(z), q25 = unname(quantile(z, .25)), q75 = unname(quantile(z, .75))),
       by = .(cl, interp, x)]
setorder(S, cl, x)
S[, seg := paste(cl, interp, sep = "_")]
# no-donor interval: last supported grid age before the gap to the first supported grid age after it
GAP_LO <- max(age_all[supp & age_all > 0 & age_all < 5])
GAP_HI <- min(age_all[supp & age_all > 5])
UNS_XC <- c(x_of_age(GAP_LO), x_of_age(GAP_HI))

pc1 <- ggplot(S, aes(x, med)) +
  annotate("rect", xmin = UNS_XC[1], xmax = UNS_XC[2], ymin = -Inf, ymax = Inf, fill = "#F1EFEA") +
  annotate("text", x = mean(UNS_XC), y = 1.08, size = 4.8 / .pt, colour = MUTED,
           fontface = "italic", label = "no direct sampling", lineheight = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = LINE_COL, linetype = "dashed") +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = cl, group = seg), alpha = 0.13, colour = NA) +
  geom_line(data = S[interp == FALSE], aes(colour = cl, group = seg), linewidth = 0.5) +
  geom_line(data = S[interp == TRUE],  aes(colour = cl, group = seg), linewidth = 0.5, linetype = "11") +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_fill_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), guide = "none") +
  scale_x_continuous(limits = c(-0.075, 1), expand = c(0, 0), breaks = NULL) +
  scale_y_continuous(limits = c(-1.25, 1.25), breaks = c(-1, 0, 1), labels = fmt_minus, expand = c(0, 0)) +
  labs(x = NULL, y = "Median z (IQR)") +
  theme(axis.line.x = element_blank(), plot.margin = margin(1, 1, 0, 1, "mm"))

pc2 <- ggplot()
pc2 <- pc2 + annotation_raster(toras_cont(M_z[ord[[1]], ]), xmin = min(xa), xmax = max(xa),
                               ymin = yb[[1]][1], ymax = yb[[1]][2], interpolate = FALSE)
pc2 <- pc2 + annotation_raster(toras_cont(M_z[ord[[2]], ]), xmin = min(xa), xmax = max(xa),
                               ymin = yb[[2]][1], ymax = yb[[2]][2], interpolate = FALSE)

tk <- data.table(age = c((14 - 40) * 7 / 365, 0, 11.4, 30, 60, 100),   # GW14 on gw = 40 + age * 365 / 7
                 lab = c("GW14", "Birth", "11 y", "30 y", "60 y", "100 y"))
tk[, key := lab %in% c("Birth", "11 y")]
tk[, x := x_of_age(age)]
tk[, hj := fifelse(lab == "Birth", 0.5, fifelse(lab == "GW14", 0, 0.5))]

pc2 <- pc2 +
  annotate("rect", xmin = UNS_XC[1], xmax = UNS_XC[2], ymin = 0, ymax = 1, fill = scales::alpha("#F1EFEA", 0.55)) +
  annotate("segment", x = UNS_XC[1], xend = UNS_XC[1], y = 0, yend = 1, linetype = "dotted", colour = "#C9C4B8", linewidth = 0.3) +
  annotate("segment", x = UNS_XC[2], xend = UNS_XC[2], y = 0, yend = 1, linetype = "dotted", colour = "#C9C4B8", linewidth = 0.3) +
  annotate("text", x = -0.068, y = 0.5, label = "APA events", angle = 90, size = T6, colour = INK, vjust = 0) +
  annotate("rect", xmin = -0.062, xmax = -0.024, ymin = c(yb[[1]][1], yb[[2]][1]), ymax = c(yb[[1]][2], yb[[2]][2]),
           fill = c(C1_COLOR, C2_COLOR)) +
  annotate("text", x = -0.043, y = c(mean(yb[[1]]), mean(yb[[2]])), angle = 90, colour = "white", fontface = "bold",
           size = 6.2 / .pt, label = c(paste0("C1 \u00b7 ", format(n_pre, big.mark = ",")),
                                       paste0("C2 \u00b7 ", format(n_post, big.mark = ",")))) +
  annotate("segment", x = tk$x, xend = tk$x, y = 0, yend = -0.013, linewidth = 0.35, colour = INK) +
  annotate("text", x = tk$x, y = -0.02, label = tk$lab, vjust = 1, hjust = tk$hj,
           size = ifelse(tk$key, 6.2, 5.8) / .pt,
           colour = INK, fontface = ifelse(tk$key, "bold", "plain")) +
  coord_cartesian(xlim = c(-0.075, 1), ylim = c(-0.058, 1), expand = FALSE, clip = "off") +
  geom_tile(data = data.table(x = 0, y = 0, z = c(Q_LO, Q_HI)), aes(x, y, fill = z), alpha = 0) +
  scale_fill_gradientn(colours = HEAT_CONTINUOUS, limits = c(Q_LO, Q_HI),
                       breaks = round(c(Q_LO, 0, Q_HI), 1),
                       labels = c(sprintf("%.1f", Q_LO), "0", sprintf("%.1f", Q_HI)),
                       name = "PDUI (z-score)") +
  theme_void(base_size = BS) +
  theme(legend.position = "bottom", legend.key.width = unit(8, "mm"), legend.key.height = unit(1.5, "mm"),
        legend.margin = margin(3, 0, 0, 0, "mm"), legend.title = element_text(size = 5.8, vjust = 0.8),
        legend.text = element_text(size = 5.4), plot.margin = margin(0, 1, 1, 1, "mm"))

# --- d: gene-set overlap (C1-specific, shared, C2-specific) ---
g  <- fread(file.path(D2, "05_gene_aggregation.tsv"))
g[, cls := factor(fcase(gene_state == "C1-spec", "C1-specific", gene_state == "C2-spec", "C2-specific", membership == "Shared", "Shared"), names(CLS))]
g <- g[!is.na(cls)]
ng <- g[, .N, by = cls]

circ <- function(x0, r = 1, n = 360) { t <- seq(0, 2 * pi, length.out = n); data.table(x = x0 + r * cos(t), y = r * sin(t)) }
CA <- circ(-0.55); CB <- circ(0.55)
lens <- rbind(CA[(x - 0.55)^2 + y^2 <= 1][order(atan2(y, x - (-0.55)))],
              CB[(x + 0.55)^2 + y^2 <= 1][order(atan2(y, x - 0.55))])
lens <- lens[order(atan2(y, x))]

ngv <- setNames(ng$N, as.character(ng$cls))
vt <- data.table(x = c(-1.0, 0, 1.0), cls = c("C1-specific", "Shared", "C2-specific"),
                 lab = c("C1-spec.", "Shared", "C2-spec."))
vt[, n := ngv[cls]]

pdv <- ggplot() +
  geom_polygon(data = CA, aes(x, y), fill = scales::alpha(C1_COLOR, 0.18), colour = NA) +
  geom_polygon(data = CB, aes(x, y), fill = scales::alpha(C2_COLOR, 0.18), colour = NA) +
  geom_polygon(data = lens, aes(x, y), fill = scales::alpha(SHARED_COLOR, 0.26), colour = NA) +
  geom_path(data = CA, aes(x, y), colour = C1_COLOR, linewidth = 0.5) +
  geom_path(data = CB, aes(x, y), colour = C2_COLOR, linewidth = 0.5) +
  geom_text(data = vt, aes(x, 0.2, label = lab, colour = cls), size = T6, fontface = "bold") +
  geom_text(data = vt, aes(x, -0.12, label = format(n, big.mark = ",")), size = 6.8 / .pt, fontface = "bold", colour = INK) +
  scale_colour_manual(values = CLS, guide = "none") +
  coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.22, 1.16), clip = "off") +
  theme_void(base_size = BS) +
  theme(plot.margin = margin(0, 0, 0, 0, "mm"))

# --- e: enrichment heatmap (Fisher OR; asterisk = covariate-adjusted FDR < 0.05) ---
ea <- fread(file.path(D2, "08_enrich_adjusted.tsv"))
eg <- ea[adjusted_fdr < 0.05 & source_group == "GO/KEGG" & level == "Five-state"]
c2named <- c("KEGG_PROTEASOME", "GOBP_OXIDATIVE_PHOSPHORYLATION")
de <- rbind(eg[gene_set == "C1-spec"], eg[gene_set == "C2-spec" & term %in% c2named],
            eg[gene_set == "C2-spec" & !(term %in% c2named)][order(adjusted_fdr, -adjusted_odds_ratio)][1:2])
clean <- tolower(gsub("_", " ", sub("^(GOBP|GOMF|KEGG)_", "", de$term)))
clean <- paste0(toupper(substring(clean, 1, 1)), substring(clean, 2))
clean <- sub("G protein coupled receptor activity", "GPCR activity", clean)
clean <- sub("Structural constituent of ribosome", "Ribosome constituent", clean)
de[, lab := clean]
de[, cls := fifelse(gene_set == "C1-spec", "C1-specific", "C2-specific")]

two <- fread(file.path(D2, "08_enrich_two_levels_1016.tsv"))
setnames(two, "odds_ratio", "fisher_odds_ratio")
TERMS <- data.table(term = c(de$term, "Alzheimers Disease"), lab = c(de$lab, "AD risk genes"),
                    blk = c(ifelse(de$cls == "C1-specific", "C1 terms", "C2 terms"), "Disease"))
mx <- two[term %in% TERMS$term & ((level == "Five-state" & gene_set %in% c("C1-spec", "C2-spec")) | (level == "Directional membership" & gene_set == "Shared"))]
mx[, cls := factor(c("C1-spec" = "C1-specific", "C2-spec" = "C2-specific", "Shared" = "Shared")[gene_set], names(CLS))]
mx <- merge(mx, TERMS, by = "term")
adr <- ea[source_group == "GWAS/HPO" & level == "Membership" & gene_set == "Shared" & term == "Alzheimers Disease"]
adjs <- rbind(ea[adjusted_fdr < 0.05, .(term, gs = gene_set)])
mx[, star := paste(term, gene_set) %in% paste(adjs$term, adjs$gs) | (term == "Alzheimers Disease" & gene_set == "Shared" & adr$adjusted_fdr < 0.05)]
mx[, show_or := fisher_odds_ratio]

labord <- c("Ribosome", "Cytoplasmic translation", "Ribosome constituent", "GPCR activity",
            "Oxidative phosphorylation", "Aerobic respiration", "Proteasome", "Ligase activity", "AD risk genes")
mx[, lab := factor(as.character(lab), rev(labord))]
mx[, l2 := pmax(pmin(log2(show_or), 3), -3)]

CLS_TBL <- rbindlist(lapply(c("C1spec", "C2spec"), function(gn)
  fread(file.path(D2, sprintf("T1_term_clusters_%s_fdr10.tsv", gn)))))
CLU_LAB <- c(C1spec_01 = "Translation", C1spec_02 = "GPCR", C2spec_01 = "Energy", C2spec_08 = "Ligase", C2spec_14 = "Proteasome")
mx[, cat := {cid <- CLS_TBL$cluster_id[match(term, CLS_TBL$term)]; ifelse(is.na(cid), "Disease", CLU_LAB[cid])}]
BAND_CLS <- c(Translation = "C1-specific", GPCR = "C1-specific", Energy = "C2-specific", Proteasome = "C2-specific", Ligase = "C2-specific", Disease = "Shared")
bands <- unique(mx[, .(lab, cat, cls = BAND_CLS[cat])])[, y := as.integer(lab)][, .(y0 = min(y) - 0.5, y1 = max(y) + 0.5), by = .(cat, cls)]

mx[, fill_val := ifelse(star, l2, NA_real_)]
pe_ <- ggplot(mx, aes(cls, lab)) +
  geom_tile(fill = "#FAFAFA", colour = "#FFFFFF", linewidth = 0.9) +
  geom_tile(aes(fill = fill_val), colour = "#FFFFFF", linewidth = 0.9) +
  geom_text(aes(label = sprintf("%.1f%s", show_or, ifelse(star, "*", "")),
                fontface = ifelse(star, "bold", "plain"),
                colour = ifelse(star, INK, "#A0A0A0")), size = 6.3 / .pt) +
  geom_hline(yintercept = c(1.5, 5.5), colour = "#E8E8E8", linewidth = 0.7) +
  annotate("rect", xmin = 0.38, xmax = 0.46, ymin = bands$y0 + 0.08, ymax = bands$y1 - 0.08,
           fill = scales::alpha(CLS[bands$cls], 0.25)) +
  annotate("text", x = 0.32, y = (bands$y0 + bands$y1) / 2, label = sub(" ", "\n", bands$cat),
           hjust = 1, size = 6 / .pt, colour = CLS[bands$cls], fontface = "bold", lineheight = 0.85) +
  scale_fill_gradientn(colours = HEAT_CONTINUOUS, limits = c(-3, 3), breaks = c(-3, 0, 3),
                       labels = c("1/8", "1", "8"), name = "Fisher OR (sig.)  ", na.value = "#FFFFFF") +
  scale_colour_identity() +
  scale_x_discrete(position = "top", labels = c("C1-specific" = "C1-spec.", "Shared" = "Shared", "C2-specific" = "C2-spec.")) +
  scale_y_discrete(expand = c(0, 0), position = "right") +
  coord_cartesian(xlim = c(0.5, 3.5), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(face = "bold", size = 7, colour = CLS[c("C1-specific", "Shared", "C2-specific")]),
        axis.text.y = element_text(size = 6, colour = INK),
        legend.position = "bottom", legend.key.width = unit(4.5, "mm"), legend.key.height = unit(1.5, "mm"),
        legend.title = element_text(size = 5.8), legend.text = element_text(size = 5.4), plot.margin = margin(1, 1, 1, 14, "mm"))

# --- f: LOEUF by gene set ---
ld <- unique(fread(file.path(D2, "ed4c_loeuf.tsv"))[, .(gene_id_base, loeuf)])
gl <- merge(g, ld, by = "gene_id_base")[!is.na(loeuf)]
# current LOEUF model (refit; supersedes the root-level copy)
mt <- fread(file.path(D2, "11_loeuf_adj_models_outputs.tsv")); b <- mt[term == "membershipShared"]
pef <- function(x) { e <- floor(log10(x)); m <- round(x / 10^e, 1); sprintf("%.1f %%*%% 10^\"%s\"", m, sub("-", MINUS, sprintf("%d", e), fixed = TRUE)) }

pf1 <- ggplot(gl, aes(cls, loeuf)) +
  geom_violin(aes(fill = cls), colour = NA, scale = "width", width = 0.85, trim = TRUE) +
  geom_boxplot(aes(colour = cls), width = 0.2, outlier.shape = NA, coef = 0, fill = "#FFFFFF", linewidth = 0.4) +
  annotate("segment", x = c(1, 3, 1), xend = c(1, 3, 3), y = c(2.45, 2.45, 2.7), yend = c(2.7, 2.7, 2.7), linewidth = 0.25, colour = INK) +
  annotate("segment", x = 2, xend = 2, y = 2.7, yend = 2.95, linewidth = 0.25, colour = INK) +
  annotate("text", x = 2, y = 4.35, vjust = 0, size = 5.6 / .pt, colour = INK, label = "Shared vs specific") +
  annotate("text", x = 2, y = 3.08, vjust = 0, size = 5.6 / .pt, colour = INK,
           label = paste0("beta == \"", sub("-", MINUS, sprintf("%.3f", b$estimate), fixed = TRUE), "\"*','~FDR == ", pef(b$fdr)), parse = TRUE) +
  scale_fill_manual(values = TINT, guide = "none") +
  scale_colour_manual(values = CLS, guide = "none") +
  scale_y_log10(breaks = c(0.2, 0.5, 1, 2)) +
  scale_x_discrete(labels = c("C1-specific" = "C1-spec.", "Shared" = "Shared", "C2-specific" = "C2-spec.")) +
  coord_cartesian(ylim = c(0.15, 5.7), clip = "on") +
  labs(x = NULL, y = "LOEUF") +
  theme(axis.text.x = element_text(size = 5.8, face = "bold"))

# --- assemble ---
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0.005, 0.995))
row1 <- (wrap_elements(full = pa) | tg(pb, "b")) + plot_layout(widths = c(1.5, 1))
pc   <- (tg(pc1, "c") / pc2) + plot_layout(heights = c(0.78, 3.6))
rcol <- (wrap_elements(full = tg(pdv, "d")) / wrap_elements(full = tg(pe_, "e")) / wrap_elements(full = tg(pf1, "f"))) + plot_layout(heights = c(0.92, 1.48, 0.95))
row2 <- (wrap_elements(full = pc) | wrap_elements(full = rcol)) + plot_layout(widths = c(1.12, 1))
fig  <- (row1 / row2) + plot_layout(heights = c(0.6, 1.62))
ggsave(file.path(FIG_DIR, "fig1.png"), fig, width = 180, height = 200, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig1.pdf"), fig, width = 180, height = 200, units = "mm", device = cairo_pdf)
