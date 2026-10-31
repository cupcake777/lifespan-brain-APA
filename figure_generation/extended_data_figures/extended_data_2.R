# ED Fig 2: batch x stage structure (a), method rank map (b), PCA before and after ComBat_REF (c).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg); library(scales)})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
DD <- file.path(ROOT, "data/figures/ed2")
theme_set(theme_nature(base_size = 7) + theme(
  axis.text = element_text(size = 7, colour = INK),
  axis.title = element_text(size = 8, colour = INK),
  axis.line   = element_line(linewidth = 0.25, colour = INK),
  axis.ticks  = element_line(linewidth = 0.25, colour = INK),
  axis.ticks.length = unit(0.8, "mm"),
  legend.text  = element_text(size = 6.5, colour = INK),
  legend.title = element_text(size = 6.5, colour = INK),
  strip.text   = element_text(size = 7, colour = INK),
  plot.margin  = margin(1, 1.5, 1, 1.5, "mm")))
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 9, face = "bold", family = FONT),
                                               plot.tag.position = c(0, 1))
BATCH_LV <- paste0("batch", 1:9)
BATCH_SHAPE <- setNames(c(1, 0, 2, 5, 6, 21, 22, 23, 24), BATCH_LV)

# --- a: batch x stage structure ---
bx <- fread(file.path(DD, "ed2a_batch_stage.tsv"))
L <- melt(bx, id.vars = "batch", variable.name = "stage", value.name = "n")
L <- L[n > 0]
L[, stage := factor(stage, c("Prenatal", "Postnatal"))]
L[, batch := factor(batch, rev(BATCH_LV))]
L[, shp := BATCH_SHAPE[as.character(batch)]]
SYM_X <- -10   # glyph column, close to the bars
pa <- ggplot(L, aes(n, batch)) +
  geom_col(aes(fill = stage), width = 0.5, position = position_stack(), colour = NA) +
  geom_text(aes(label = n, group = stage), position = position_stack(vjust = 0.5), size = 6.5 / .pt,
            colour = INK, family = FONT) +
  geom_point(data = unique(L[, .(batch, shp)]), aes(x = SYM_X, y = batch, shape = shp),
             size = 1.3, colour = INK, fill = INK, stroke = 0.25, inherit.aes = FALSE) +
  annotate("text", x = SYM_X, y = Inf, label = "Symbol", vjust = -1.1, hjust = 0.5,
           size = 6.5 / .pt, fontface = "bold", colour = INK, family = FONT) +
  scale_shape_identity() +
  scale_fill_manual(values = STAGE_COLORS[c("Prenatal", "Postnatal")], name = "Life stage") +
  scale_y_discrete(labels = function(x) sub("batch9", "batch 9 (bridge)", sub("batch", "batch ", x))) +
  scale_x_continuous(breaks = c(0, 100, 200), limits = c(-18, 200), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "APA profiles", y = NULL) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 7, hjust = 1), legend.position = "bottom",
        legend.key.size = unit(2, "mm"), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(1, 1.5, 1, 1.5, "mm"))

# --- b: method x metric rank map (greyscale) ---
ms <- fread(file.path(DD, "ed2b_method_scores.tsv"))
meth_lv <- c("ComBat_REF", "ComBat_REF_NP", "limma", "RUVs", "RUV-III", "ComBat_RAW", "ComBat_REF_meanonly")
wt <- unique(ms[, .(method, weighted_total)])
short <- c("ComBat_REF" = "ComBat\nREF", "ComBat_REF_NP" = "ComBat\nREF-NP", "limma" = "limma", "RUVs" = "RUVs",
           "RUV-III" = "RUV-III", "ComBat_RAW" = "ComBat\nRAW", "ComBat_REF_meanonly" = "ComBat\nmean-only")
mlab <- wt[, .(method, lab = sprintf("%s\n%.2f", short[method], weighted_total))]
ms <- ms[method %in% meth_lv]
ms[, method := factor(method, meth_lv)]
mlabs <- setNames(mlab$lab, mlab$method)[meth_lv]
glv <- c("Residual batch", "Age signal", "Biological variance", "Replicate agreement")
mtxt <- c(ASW_batch = "ASW, batch", PC_batch_max_abs_r = "max |r| PC vs batch", Var_batch = "batch variance",
          Var_batch_p95 = "batch variance, p95", ResidualCoef_gt1_total = "residual coef. > 1, total",
          ResidualCoef_max_p95 = "residual coef., max p95", ResidualCoef_median = "residual coef., median",
          ResidualCoef_p95 = "residual coef., p95", AgeGAM_sig_p001 = "age GAM, p < 0.01 count",
          AgeGAM_median_partial_deviance = "age GAM, deviance", AgeGAM_median_abs_delta_logit = "age GAM, |dlogit|",
          AgeGAM_k_index_fail_fraction = "age GAM, k fail", AgeMatched_gap = "age-matched gap",
          ASW_LifeStage = "ASW, life stage", Var_LifeStage = "life-stage variance",
          CCC = "CCC, replicates", Correlation = "correlation, replicates", MAE = "MAE, replicates",
          RMSE = "RMSE, replicates", MSE = "MSE, replicates")
ms[, metric_lab := mtxt[metric]]
ms[, group := factor(group, glv)]
ms[, metric_lab := factor(metric_lab, rev(unique(mtxt[unique(ms[order(group, metric)]$metric)])))]
RANK_RAMP <- c(INK, NEUTRAL, "#FFFFFF")
pb <- ggplot(ms, aes(method, metric_lab, fill = rank)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = format(rank, nsmall = 0, drop0trailing = TRUE),
                colour = I(fifelse(rank <= 2.5, "#FFFFFF", INK))), size = 6.5 / .pt, family = FONT) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(colours = RANK_RAMP, values = rescale(c(1, 4, 7)), limits = c(1, 7), breaks = c(1, 7),
                       labels = c("1 (best)", "7 (worst)"), name = "Rank",
                       guide = guide_colourbar(direction = "horizontal", barwidth = unit(36, "mm"),
                                               barheight = unit(1.6, "mm"), title.position = "top", title.hjust = 0.5)) +
  scale_x_discrete(position = "top", labels = mlabs) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x.top = element_text(size = 6.5, lineheight = 0.95, colour = INK),
        axis.text.y = element_text(size = 6.5, colour = INK),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6, face = "bold", colour = INK),
        panel.spacing = unit(1.2, "mm"),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.5, "mm"),
        legend.box.margin = margin(t = 1.5, b = 0.5, unit = "mm"),
        plot.margin = margin(1, 1.5, 2, 1.5, "mm"))

# --- c: PCA before and after ComBat_REF ---
pc <- fread(file.path(DD, "ed2c_pca.tsv"))
pc[, condition := factor(condition, c("before", "after_ComBat_REF"), c("Before correction", "After ComBat_REF"))]
pc[, stage := factor(life_stage, c("Prenatal", "Postnatal"))]
pc[, batch := factor(batch, BATCH_LV)]
pc[, shp := BATCH_SHAPE[as.character(batch)]]
pcp <- ggplot(pc, aes(PC1, PC2, colour = stage, fill = stage, shape = batch)) +
  geom_point(size = 0.65, alpha = 0.7, stroke = 0.2) +
  facet_wrap(~condition, nrow = 1, scales = "free") +
  scale_x_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE))) +
  scale_y_continuous(labels = function(b) sub("^-", "\u2212", format(b, trim = TRUE))) +
  scale_shape_manual(values = BATCH_SHAPE, guide = "none") +
  scale_colour_manual(values = STAGE_COLORS[c("Prenatal", "Postnatal")], guide = "none") +
  scale_fill_manual(values = STAGE_COLORS[c("Prenatal", "Postnatal")], guide = "none") +
  labs(x = "PC1", y = "PC2") +
  theme(strip.text = element_text(size = 7.5, colour = INK))

# --- assemble ---
top <- (tg(pa, "a") | tg(pb, "b")) + plot_layout(widths = c(44, 136))
fig <- (top / tg(pcp, "c")) + plot_layout(heights = c(100, 62))
ggsave(file.path(FIG_DIR, "ed2.png"), fig, width = 180, height = 156, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed2.pdf"), fig, width = 180, height = 156, units = "mm", device = cairo_pdf)
