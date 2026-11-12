# ED Fig 6: apaQTL replication pi1 across four cohorts (a), cell-type proportions (b1), partial R2 of cell composition (b2).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(readxl); library(ragg)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
D4 <- file.path(ROOT, "data/04_qtl"); DF <- file.path(ROOT, "data/figures/ed6")
theme_set(theme_nature(base_size = 6.5, base_family = FONT) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.25, colour = INK), axis.ticks = element_line(linewidth = 0.25, colour = INK),
  plot.margin = margin(1, 2, 1, 1, "mm")))
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold", family = FONT),
                                               plot.tag.position = c(0, 1))
LEG <- theme(legend.position = "top", legend.key.size = unit(2, "mm"), legend.margin = margin(0, 0, 0, 0),
             legend.text = element_text(size = 5.5), legend.box.spacing = unit(0.5, "mm"))
SC <- STAGE_COLORS[c("Prenatal", "Postnatal")]
SX <- file.path(D4, "Supplementary_Tables_draft.xlsx")

# --- a: replication heatmap, discovery x target (pi1 of apaQTL effects) ---
COH <- c("prenatal", "postnatal", "BrainGvex adult", "GTEx Brain Cortex")
SHORT <- c(prenatal = "Prenatal", postnatal = "Postnatal",
           `BrainGvex adult` = "BrainGVEX", `GTEx Brain Cortex` = "GTEx")
rep <- as.data.table(read_excel(SX, sheet = "S_replication_pi1"))
hm <- copy(rep)
stopifnot(nrow(hm) == 12)
# unname() keeps the supplied level order (factor() on a named vector keeps the source strings)
hm[, dl := factor(unname(SHORT[discovery]), rev(unname(SHORT)))]
hm[, tl := factor(unname(SHORT[target]), unname(SHORT))]
stopifnot(nlevels(hm$dl) == 4, nlevels(hm$tl) == 4)
hm[, vlab := sprintf("%.3f", pi1)]
hm[, nlab := prettyNum(matched_pairs, big.mark = ",")]
hm[, n_missing := is.na(pi1)]
# cell text is white on dark fills; luminance from the same ramp that scale_fill_gradientn uses
PI1_RAMP <- c("#F5F8FB", "#C6D9EC", "#6BA3D6", "#1F5C99")
hm[, fill_lum := NA_real_]
hm[!is.na(pi1), fill_lum := { m <- col2rgb(scales::gradient_n_pal(PI1_RAMP)(scales::rescale(pi1, from = c(0.5, 1)))) / 255; 0.299 * m[1, ] + 0.587 * m[2, ] + 0.114 * m[3, ] }]
hm[, dark := !is.na(fill_lum) & fill_lum < 0.5]
hm[, vcol := ifelse(dark, "white", INK)]
hm[, ncol := ifelse(dark, "white", MUTED)]
# diagonal (self) cells are not in the sheet; drawn as labelled grey tiles
diag_cells <- data.table(discovery = COH, target = COH)
diag_cells[, dl := factor(unname(SHORT[discovery]), rev(unname(SHORT)))]
diag_cells[, tl := factor(unname(SHORT[target]), unname(SHORT))]
pa <- ggplot() +
  geom_tile(data = diag_cells, aes(tl, dl),
            width = 0.96, height = 0.96, fill = NEUTRAL_LIGHT, colour = "white", linewidth = 0.7) +
  geom_text(data = diag_cells, aes(tl, dl, label = "self"),
            size = 5 / .pt, family = FONT, colour = MUTED) +
  geom_tile(data = hm, aes(tl, dl, fill = pi1),
            width = 0.96, height = 0.96, colour = "white", linewidth = 0.7) +
  # cells with no estimable pi1 stay blank but framed: "no value", not "no cell"
  geom_tile(data = hm[n_missing == TRUE], aes(tl, dl),
            width = 0.96, height = 0.96, fill = "grey93",
            colour = "white", linewidth = 0.7) +
  geom_text(data = hm[!is.na(pi1)], aes(tl, dl, label = vlab, colour = vcol),
            size = 5.6 / .pt, family = FONT, nudge_y = 0.08) +
  geom_text(data = hm[!is.na(pi1)], aes(tl, dl, label = nlab, colour = ncol),
            size = 4.3 / .pt, family = FONT, nudge_y = -0.20) +
  geom_text(data = hm[n_missing == TRUE], aes(tl, dl, label = "n.a."),
            size = 5 / .pt, family = FONT, colour = MUTED, nudge_y = 0.06) +
  scale_colour_identity() +
  scale_fill_gradientn(colours = PI1_RAMP,
                       limits = c(0.5, 1), name = expression(pi[1]),
                       na.value = "grey93") +
  scale_x_discrete(labels = SHORT, expand = expansion(add = 0.02)) +
  scale_y_discrete(labels = SHORT, expand = expansion(add = 0.02)) +
  labs(x = "Target cohort", y = "Discovery cohort") +
  theme(
    axis.text = element_text(size = 5.6),
    axis.title.x = element_text(size = 6, margin = margin(b = 1.2)),
    axis.title.y = element_text(size = 6, margin = margin(l = 1.2)),
    legend.position = "right",
    legend.key.height = unit(4.5, "mm"),
    legend.key.width = unit(1.6, "mm"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    legend.margin = margin(l = 0.8, t = 1, b = 1)
  )
pa <- tg(pa, "a")

# --- b1: cell-type proportions per stage ---
prop <- fread(file.path(DF, "ED6_cell_proportions.tsv"))
prop[, life_stage := factor(life_stage, c("Prenatal", "Postnatal"))]
# Oligodendrocyte has no prenatal value and Glial_progenitor no postnatal value: absent, not zero.
# Each absent slot gets an explicit "n.m." marker.
cell_levels <- c("Excitatory_neuron", "Inhibitory_neuron", "Astrocyte", "OPC",
                 "Oligodendrocyte", "Glial_progenitor", "Microglia_Immune",
                 "Vascular_Mural")
have <- prop[, .(have = TRUE), by = .(cell_class, life_stage)]
grid <- CJ(cell_class = cell_levels, life_stage = c("Prenatal", "Postnatal"))
grid <- merge(grid, have, by = c("cell_class", "life_stage"), all.x = TRUE)
grid[is.na(have), have := FALSE]
nm_rows <- grid[have == FALSE]
# position_dodge() centres a lone marker, so the n.m. layer carries both stage rows per class;
# the measured twin carries y = NA and is not drawn.
nm_draw <- CJ(cell_class = cell_levels, life_stage = c("Prenatal", "Postnatal"), unique = TRUE)
nm_draw <- merge(nm_draw, nm_rows[, .(cell_class, life_stage)],
                 by = c("cell_class", "life_stage"), all.x = TRUE)
nm_draw[is.na(life_stage), life_stage := "Postnatal"]
nm_draw[, life_stage := factor(life_stage, c("Prenatal", "Postnatal"))]
nm_key <- paste(nm_rows$cell_class, nm_rows$life_stage, sep = "|")
nm_draw[, is_nm := paste(cell_class, life_stage, sep = "|") %in% nm_key]
nm_draw[, y := ifelse(is_nm, 0.04, NA_real_)]
nm_draw[, cell_class := factor(cell_class, levels = cell_levels)]
alph <- sort(unique(as.character(prop$cell_class)))
xpos_of <- function(cc, st) match(as.character(cc), rev(alph)) + ifelse(as.character(st) == "Postnatal", 0.19, -0.19)
prop[, xpos := xpos_of(cell_class, life_stage)]
nm_draw <- nm_draw[is_nm == TRUE][, xpos := xpos_of(cell_class, life_stage)]
b1 <- ggplot() +
  geom_boxplot(data = prop, aes(xpos, proportion, fill = life_stage, group = interaction(cell_class, life_stage)),
               linewidth = 0.25, width = 0.34, outlier.size = 0.3, outlier.stroke = 0.1) +
  geom_point(data = nm_draw, aes(xpos, y = y, shape = life_stage),
             inherit.aes = FALSE, size = 1.2, stroke = 0.3, colour = MUTED, fill = "white") +
  geom_text(data = nm_draw, aes(xpos, y = y + 0.025, label = "n.m."),
            inherit.aes = FALSE, size = 4.6 / .pt, colour = MUTED, family = FONT, hjust = 0) +
  scale_shape_manual(values = c(Prenatal = 21, Postnatal = 24), guide = "none") +
  scale_fill_manual(values = SC, name = NULL) +
  scale_x_continuous(breaks = seq_along(alph), labels = gsub("_", " ", rev(alph)), expand = expansion(add = 0.5)) +
  coord_flip() +
  labs(x = NULL, y = "Estimated cell-type proportion") +
  theme(axis.text.y = element_text(size = 5.5)) + LEG
b1 <- tg(b1, "b")

# --- b2: partial R2 of cell composition, ECDF by stage ---
cc <- as.data.table(read_excel(SX, sheet = "S_cell_composition"))
cc[, stage := factor(stage, c("prenatal", "postnatal"), c("Prenatal", "Postnatal"))]
cs <- cc[, .(n = .N), by = stage]
lab_st <- setNames(sprintf("%s (n = %s)", cs$stage, prettyNum(cs$n, big.mark = ",")), as.character(cs$stage))
# x limit: the ECDF tail beyond ~0.14 carries almost no mass and flattens both curves
XMAX_B2 <- 0.15
b2 <- ggplot(cc, aes(partial_r2_cell_ilr, colour = stage)) +
  stat_ecdf(linewidth = 0.45) +
  scale_colour_manual(values = SC, name = NULL, labels = lab_st[c("Prenatal", "Postnatal")]) +
  scale_x_continuous(limits = c(0, XMAX_B2), breaks = c(0, 0.05, 0.10, 0.15),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = expression("Partial " * R^2 * " of cell composition"), y = "Cumulative fraction of phenotypes") +
  LEG
b2 <- tg(b2, "")

# --- assemble ---
fig <- pa | b1 | b2
fig <- fig + plot_layout(widths = c(0.34, 0.38, 0.28))
ggsave(file.path(FIG_DIR, "ed6.png"), fig, width = 180, height = 62, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed6.pdf"), fig, width = 180, height = 62, units = "mm", device = cairo_pdf)
