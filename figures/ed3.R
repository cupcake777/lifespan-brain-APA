# ED Fig 3: threshold sensitivity (a), bootstrap and seed ECDFs (b), BrainSpan trajectories (c).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg); library(scales)})
source(file.path(Sys.getenv("APA_ROOT"), "release/figures/theme.R"))
DD <- file.path(ROOT, "data/figures/ed3")
th <- theme_nature(base_size = 7) + theme(
  axis.text = element_text(size = 7, colour = INK),
  axis.title = element_text(size = 8, colour = INK),
  axis.line   = element_line(linewidth = 0.25, colour = INK),
  axis.ticks  = element_line(linewidth = 0.25, colour = INK),
  axis.ticks.length = unit(0.8, "mm"),
  legend.text  = element_text(size = 6.5, colour = INK),
  legend.title = element_text(size = 6.5, colour = INK),
  strip.text   = element_text(size = 7, colour = INK, face = "bold"),
  plot.margin  = margin(1, 1.5, 1, 1.5, "mm"))
theme_set(th)
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 9, face = "bold", family = FONT),
                                               plot.tag.position = c(0, 1))
PRIMARY_CAND <- "I04N3_O08N5"

# --- a: directional events and Jaccard across effect-size quantiles ---
m <- fread(file.path(DD, "ed3a_masks.tsv"))
m[, qf := factor(sprintf("%.2f", quantile), c("0.90", "0.95", "0.99"))]
m[, prim := candidate_id == PRIMARY_CAND & seed == 2026]
LA <- melt(m, id.vars = c("qf", "prim"), measure.vars = c("n_directional", "directional_jaccard_vs_primary"),
           variable.name = "metric", value.name = "v", na.rm = TRUE)
LA[, metric := factor(metric, c("n_directional", "directional_jaccard_vs_primary"),
                      c("Directional events", "Jaccard vs primary"))]
MA <- LA[, .(v = mean(v)), by = .(qf, metric)]
PA <- LA[prim == TRUE]
PA[, txt := fifelse(metric == "Directional events", format(round(v), big.mark = ","), sprintf("%.2f", v))]
# label placement: peaks above, valleys and steep endpoints below
PA[, vj := fifelse(metric == "Directional events",
                   fifelse(qf == "0.90", -1.25, fifelse(qf == "0.95", -1.5, 1.9)),
                   fifelse(qf == "0.95", -1.5, 1.9))]
pa <- ggplot() +
  annotate("rect", xmin = 1.62, xmax = 2.38, ymin = -Inf, ymax = Inf, fill = NEUTRAL_LIGHT, alpha = 0.28) +
  geom_line(data = MA, aes(qf, v, group = metric), linewidth = 0.55, colour = INK) +
  geom_point(data = LA, aes(qf, v), colour = NEUTRAL, size = 0.9, alpha = 0.55, shape = 16) +
  geom_point(data = PA, aes(qf, v), colour = INK, size = 2.1, shape = 16) +
  geom_text(data = PA, aes(qf, v, label = txt, vjust = vj), size = 6.2 / .pt, colour = INK,
            family = FONT) +
  facet_wrap(~metric, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.26, 0.22)),
                     labels = function(x) format(x, big.mark = ",", trim = TRUE)) +
  labs(x = "Effect-size quantile", y = NULL) +
  theme(axis.text.x = element_text(size = 6.5), panel.spacing = unit(3, "mm"))

# --- b: bootstrap (direction, assignment) and seed ECDFs ---
rb <- fread(file.path(DD, "09_sample_boot_site_robust.tsv"))
rb[, cls := fifelse(original_cluster == 1, "C1", "C2")]
LB <- melt(rb, id.vars = "cls", measure.vars = c("direction_reproducibility", "assignment_reproducibility"),
           variable.name = "set", value.name = "r")
LB[, set := factor(set, c("direction_reproducibility", "assignment_reproducibility"), c("Direction", "Assignment"))]
s <- fread(file.path(DD, "ed3b_seed.tsv"))
LS <- rbindlist(list(
  data.table(set = "Seeds", r = s$cluster_ari),
  data.table(set = "Seeds", r = s$gene_state_agreement)))
LB[, set := factor(set, c("Direction", "Assignment", "Seeds"))]
LS[, set := factor(set, c("Direction", "Assignment", "Seeds"))]
marks <- data.table(set = factor(c("Direction", "Assignment"), c("Direction", "Assignment", "Seeds")), x = 0.8)
ANN_DIR <- mean(rb$direction_reproducibility >= 0.8)
ANN_ASG <- rb[, .(f = mean(assignment_reproducibility >= 0.8)), by = cls]$f
note <- data.table(
  set = factor(c("Direction", "Assignment", "Seeds"), c("Direction", "Assignment", "Seeds")),
  x = c(0.47, 0.47, min(LS$r)), y = c(0.95, 0.95, 0.95),
  txt = c(sprintf("%.1f%% \u2265 0.8", 100 * ANN_DIR),
          sprintf("%.1f\u2013%.1f%% \u2265 0.8", 100 * min(ANN_ASG), 100 * max(ANN_ASG)),
          sprintf("min %.5f", min(LS$r))))
# bootstrap facets on 0.6/0.8/1.0; seeds facet on 0.996-1.0
pb <- ggplot() +
  geom_vline(data = marks, aes(xintercept = x), linetype = "22", colour = NEUTRAL, linewidth = 0.3) +
  geom_step(data = LB, aes(r, colour = cls), stat = "ecdf", direction = "hv", linewidth = 0.5) +
  geom_step(data = LS, aes(r), stat = "ecdf", direction = "hv", linewidth = 0.5, colour = INK) +
  geom_text(data = note, aes(x, y, label = txt), hjust = 0, size = 6 / .pt, colour = MUTED, family = FONT) +
  facet_wrap(~set, scales = "free_x") +
  scale_colour_manual(values = c(C1 = C1_COLOR, C2 = C2_COLOR), name = NULL, breaks = c("C1", "C2")) +
  scale_y_continuous(limits = c(0, 1.04), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_x_continuous(breaks = function(lims) if (lims[1] > 0.9) c(0.996, 0.998, 1.0) else c(0.6, 0.8, 1.0), expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "Reproducibility (r)", y = "Cumulative fraction") +
  theme(axis.text.x = element_text(size = 6.5), panel.spacing = unit(3, "mm"),
        axis.title.x = element_text(size = 7), legend.position = "none")

# --- c: BrainSpan trajectories vs discovery, C1 and C2 ---
ev <- fread(file.path(DD, "ed3c_events.tsv"))
tp <- fread(file.path(DD, "ed3c_traj.tsv"))
ev[, cls := factor(ifelse(discovery_cluster_id == 1, "C1", "C2"), c("C1", "C2"))]
tp[, cls := factor(ifelse(discovery_cluster_id == 1, "C1", "C2"), c("C1", "C2"))]
PRE <- tp$segment_id == 1
tp[, x := ifelse(PRE, 0.38 * (40 + age * 52.14 - 16) / 16, 0.44 + 0.56 * (age - 13) / 27)]
sm <- ev[, .(n = .N, nq = sum(external_age_q < 0.05, na.rm = TRUE),
             nconc = sum(discovery_cluster_direction_concordant == TRUE, na.rm = TRUE)), by = cls]
sm[, txt := sprintf("Concordant: %d/%d (%.0f%%)  [q < 0.05: %d/%d]", nconc, n, 100*nconc/n, nq, n)]
br <- rbindlist(lapply(levels(tp$cls), function(cl) {
  p <- tp[cls == cl & PRE][which.max(age)]; q <- tp[cls == cl & !PRE][which.min(age)]
  rbind(data.table(cls = cl, series = "ext", x0 = p$x, x1 = q$x, y0 = p$external_z, y1 = q$external_z),
        data.table(cls = cl, series = "disc", x0 = p$x, x1 = q$x, y0 = p$discovery_z, y1 = q$discovery_z))
}))
cols <- c(C1 = C1_COLOR, C2 = C2_COLOR)
# GW24 dropped (collides at panel width); last/first ticks edge-aligned per position
x_ticks <- c(0, 0.38, 0.44, 0.586, 0.793, 1.0)
x_labels <- c("GW16", "GW32", "13 y", "20 y", "30 y", "40 y")
x_hj <- c(0, 1, 0, 0.62, 0.5, 1)   # 20 y nudged right to clear '13 y'
pcc <- ggplot(tp, aes(x, colour = cls, fill = cls)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = GRID_COL, linewidth = 0.3) +
  annotate("rect", xmin = 0.38, xmax = 0.44, ymin = -Inf, ymax = Inf, fill = GRID_COL, alpha = 0.5) +
  annotate("segment", x = 0.38, xend = 0.38, y = -Inf, yend = Inf, linetype = "dotted", colour = NEUTRAL_LIGHT, linewidth = 0.3) +
  annotate("segment", x = 0.44, xend = 0.44, y = -Inf, yend = Inf, linetype = "dotted", colour = NEUTRAL_LIGHT, linewidth = 0.3) +
  geom_ribbon(aes(ymin = external_ci_low, ymax = external_ci_high, group = interaction(cls, segment_id)),
              alpha = 0.16, colour = NA) +
  geom_line(aes(y = external_z, group = interaction(cls, segment_id), linetype = "BrainSpan (95% CI)"), linewidth = 0.5) +
  geom_line(aes(y = discovery_z, group = interaction(cls, segment_id), linetype = "Discovery"), linewidth = 0.35, colour = INK) +
  geom_segment(data = br[series == "ext"], aes(x = x0, xend = x1, y = y0, yend = y1, colour = cls),
               linewidth = 0.4, linetype = "22", inherit.aes = FALSE) +
  geom_segment(data = br[series == "disc"], aes(x = x0, xend = x1, y = y0, yend = y1),
               linewidth = 0.35, linetype = "22", colour = INK, inherit.aes = FALSE) +
  geom_text(data = sm, aes(x = 0.02, y = Inf, label = txt), inherit.aes = FALSE, hjust = 0, vjust = 1.1,
            size = 6.2 / .pt, colour = INK, family = FONT, lineheight = 1.12) +
  facet_wrap(~cls, nrow = 1) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_linetype_manual(values = c("BrainSpan (95% CI)" = "solid", "Discovery" = "22"), name = NULL) +
  geom_text(data = data.table(x = x_ticks, lab = x_labels, hj = x_hj),
            aes(x, y = -Inf, label = lab, hjust = hj), vjust = 1.6,
            size = 6.5 / .pt, colour = INK, family = FONT, inherit.aes = FALSE) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(-0.01, 1.01), breaks = x_ticks, labels = NULL) +
  scale_y_continuous(limits = c(-1.3, 1.78), breaks = c(-1, 0, 1),
                     labels = function(x) sub("^-", "\u2212", as.character(x))) +
  labs(x = "Age (GW / years)", y = "Trajectory (z)") +
  theme(
    legend.position = c(0.85, 0.88),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.spacing = unit(5, "mm"),
    axis.text.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 4.5, unit = "mm"))
  )

# --- assemble ---
top    <- (tg(pa, "a") | tg(pb, "b")) + plot_layout(widths = c(45, 55))
bottom <- tg(pcc, "c")
fig    <- top / bottom + plot_layout(heights = c(38, 56))
ggsave(file.path(FIG_DIR, "ed3.png"), fig, width = 180, height = 97, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed3.pdf"), fig, width = 180, height = 97, units = "mm", device = cairo_pdf)
