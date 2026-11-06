# ED Fig 4: adjusted enrichment of gene classes (a), background and expression-matched synaptic overlap (b), LOEUF by five-level class (c).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(ragg)})
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/shared/theme.R"))
DE <- file.path(ROOT, "data/figures/ed4"); D2 <- file.path(ROOT, "data/02_traj")
theme_set(theme_nature(base_size = 6.5, base_family = FONT) + theme(
  axis.text = element_text(size = 6, colour = INK), axis.title = element_text(size = 7, colour = INK),
  axis.line = element_line(linewidth = 0.25, colour = INK), axis.ticks = element_line(linewidth = 0.25, colour = INK),
  plot.margin = margin(1, 2, 1, 1, "mm")))
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold", family = FONT),
                                               plot.tag.position = c(0, 1))
CLS_COL <- c("C1-spec" = C1_COLOR, "C2-spec" = C2_COLOR, "Shared" = SHARED_COLOR)

# --- a: adjusted enrichment, one row per significant term ---
ea <- fread(file.path(DE, "ed4a_adjusted_sig.tsv"))
ea <- ea[gene_set %in% c("C1-spec", "C2-spec", "Shared")]
stopifnot(nrow(ea) == 39, all(ea$adjusted_fdr < 0.05))
lab_one <- function(term, src) {
  if (grepl("HP:", term)) return(sub("\\s*\\(HP:[0-9]+\\)", "", term))
  if (src == "GWAS/HPO") return("GWAS Alzheimer's disease")
  body <- tolower(gsub("_", " ", sub("^(GOBP|GOMF|KEGG)_", "", term)))
  body <- paste0(toupper(substring(body, 1, 1)), substring(body, 2))
  body <- gsub("\\b[Aa]tp\\b", "ATP", body, perl = TRUE)
  body <- sub("G protein coupled", "G protein-coupled", body)
  body <- sub("Iron sulfur", "iron-sulfur", body)
  body <- sub("Metal sulfur", "metal-sulfur", body)
  body <- sub("Disease$", "disease", body)
  body
}
rows <- ea[, .(gene_set, term, adjusted_odds_ratio, adjusted_fdr, source_group)]
stopifnot(nrow(rows) == 39, uniqueN(rows$term) == 39)
rows[, lab := mapply(lab_one, term, source_group)]
rows[, src := fcase(grepl("^GOBP_|^GOMF_", term), "GO",
                    grepl("^KEGG_", term), "KEGG",
                    default = "GWAS / HPO")]
rows[, src := factor(src, c("GO", "KEGG", "GWAS / HPO"))]
rows[, cls := factor(gene_set, c("C1-spec", "C2-spec", "Shared"))]
rows[, mlogfdr := -log10(adjusted_fdr)]
rows <- rows[order(src, cls, -adjusted_odds_ratio)]
rows[, ynum := seq_len(.N)]
band <- rows[, .(ymin = min(ynum) - 0.5, ymax = max(ynum) + 0.5), by = .(src, cls)]
brk_lab <- setNames(rows$lab, as.character(rows$ynum))
# keep only rows inside each facet's own band; a neighbouring facet's row sits outside it
brk_fn <- function(lims) rows[ynum >= min(lims) + 0.6 - 1e-6 & ynum <= max(lims) - 0.6 + 1e-6, ynum]

pa <- ggplot() +
  geom_rect(data = band, aes(xmin = 0.92, xmax = 17.7, ymin = ymin, ymax = ymax, fill = cls),
            alpha = 0.08, colour = NA) +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = NEUTRAL, linetype = "22") +
  geom_segment(data = rows, aes(x = 1, xend = adjusted_odds_ratio, y = ynum, yend = ynum, colour = cls),
               linewidth = 0.25) +
  geom_point(data = rows, aes(adjusted_odds_ratio, ynum, colour = cls, size = mlogfdr)) +
  scale_colour_manual(values = CLS_COL, guide = "none") +
  scale_fill_manual(values = CLS_COL, guide = "none") +
  scale_size_continuous(range = c(0.8, 2.6), limits = c(0, 8), name = expression(-log[10]~adj.~FDR), breaks = c(2, 4, 6, 8)) +
  scale_x_log10(limits = c(0.9, 18), breaks = c(1, 2, 4, 8, 16), labels = c("1", "2", "4", "8", "16")) +
  scale_y_continuous(trans = "reverse", breaks = brk_fn,
                     labels = function(x) brk_lab[as.character(x)], expand = expansion(add = 0.6)) +
  facet_grid(src ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Adjusted odds ratio (log)", y = NULL) +
  theme(axis.text.y = element_text(size = 5.8, colour = INK),
        strip.text.y.right = element_text(angle = 0, size = 5.8, margin = margin(l = 1.2, r = 1.2)),
        panel.spacing.y = unit(1.4, "mm"),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_text(size = 5.5), legend.text = element_text(size = 5.5),
        legend.key.size = unit(2.2, "mm"),
        legend.margin = margin(t = -0.5, r = 0, b = 0, l = 0, unit = "mm"),
        legend.background = element_blank(), legend.key = element_blank())

# --- b: brain-expressed/enriched backgrounds and SynGO overlap, raw and expression-matched ---
bg <- fread(file.path(DE, "ed4b_background.tsv"))
bgd <- bg[, .(grp = background, lab = cls, OR, lo, hi, p = fdr, variant = "Raw")]
sm <- fread(file.path(D2, "T1_synaptic_models.tsv"))
sy <- sm[gene_class == "Shared" & gene_set == "SynGO_all" & model == "unadjusted"]
stopifnot(nrow(sy) == 1)
mt <- fread(file.path(D2, "T1_matched_synaptic.tsv"))
a_ <- sum(mt$syn_shared == TRUE & mt$control_shared == FALSE)
b_ <- sum(mt$syn_shared == FALSE & mt$control_shared == TRUE)
pm <- binom.test(a_, a_ + b_, 0.5)
or_m <- a_ / b_; ci_m <- pm$conf.int / (1 - pm$conf.int)
synd <- rbind(
  data.table(grp = "SynGO", lab = "Shared", OR = sy$odds_ratio, lo = sy$ci_low, hi = sy$ci_high,
             p = sy$p, variant = "Raw"),
  data.table(grp = "SynGO", lab = "Shared", OR = or_m, lo = ci_m[1], hi = ci_m[2],
             p = pm$p.value, variant = "Expr.-matched"))
db <- rbind(bgd, synd, fill = TRUE)
db[, grp := factor(grp, c("Brain-expressed", "Brain-enriched", "SynGO"))]
db[, cls := factor(lab, c("C1-spec", "C2-spec", "Shared"))]
db[, yv := as.numeric(match(lab, c("C1-spec", "C2-spec", "Shared")))]
db[grp == "SynGO", yv := yv + fifelse(variant == "Raw", 0.16, -0.16)]
db[, sig := p < 0.05]
db[, colr := fifelse(grp == "SynGO", SHARED_COLOR, CLS_COL[as.character(lab)])]

pb <- ggplot(db, aes(OR, yv, colour = colr)) +
  geom_vline(xintercept = 1, linewidth = 0.2, colour = NEUTRAL, linetype = "22") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.3) +
  geom_point(aes(shape = variant, fill = variant), size = 1.8, stroke = 0.35) +
  scale_colour_identity() +
  scale_fill_manual(values = c("Raw" = SHARED_COLOR, "Expr.-matched" = "white"), guide = "none") +
  scale_shape_manual(values = c("Raw" = 16, "Expr.-matched" = 21), name = NULL,
                     breaks = c("Raw", "Expr.-matched")) +
  facet_wrap(~grp, ncol = 3) +
  scale_x_log10(limits = c(0.2, 32), breaks = c(0.3, 1, 3, 10, 30), labels = c("0.3", "1", "3", "10", "30")) +
  scale_y_continuous(breaks = 1:3, labels = c("C1-spec", "C2-spec", "Bidirectional"), expand = expansion(add = 0.7)) +
  labs(x = "Adjusted odds ratio (log)", y = NULL) +
  theme(strip.text = element_text(size = 5.8, margin = margin(b = 1)),
        strip.clip = "off",   # \"Brain-expressed\" overflows its facet; the gap is empty there
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.2, "mm"), legend.text = element_text(size = 5.5),
        legend.margin = margin(t = 0, r = 0, b = 2, l = 0, unit = "mm"),
        legend.background = element_blank(), legend.key = element_blank(),
        plot.margin = margin(t = 1, r = 2, b = 2, l = 1, unit = "mm"),
        panel.spacing.x = unit(3, "mm"))

# --- c: LOEUF by five-level gene class ---
lo <- fread(file.path(D2, "ed4c_loeuf.tsv"))
lv <- c("C1-spec", "C1-main", "Balanced", "C2-main", "C2-spec")
lo[, gs := factor(gene_state, lv)]
ns <- lo[, .N, by = gs]
cols <- c("C1-spec" = C1_COLOR, "C1-main" = scales::alpha(C1_COLOR, 0.5), "Balanced" = NEUTRAL,
          "C2-main" = scales::alpha(C2_COLOR, 0.5), "C2-spec" = C2_COLOR)
ylabs <- c("C1-\nspec", "C1-\nmain", "Balanced", "C2-\nmain", "C2-\nspec")
xlab5 <- sprintf("%s\n(n=%s)", ylabs, format(ns$N[match(lv, ns$gs)], big.mark = ","))
# the three sub-groups are subsets of Shared: one green tint for all three shadows
shd <- data.table(xmin = c(1.5, 2.5, 3.5), xmax = c(2.5, 3.5, 4.5),
                  ymin = 0.04, ymax = 2.1,
                  fill = rep(scales::alpha(SHARED_COLOR, 0.10), 3))
pc <- ggplot(lo, aes(gs, loeuf, colour = gs)) +
  geom_rect(data = shd, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE, colour = NA) +
  geom_point(position = position_jitter(width = 0.25, height = 0, seed = 20261026), size = 0.12, alpha = 0.10, stroke = 0) +
  geom_boxplot(fill = NA, width = 0.5, outlier.shape = NA, linewidth = 0.3) +
  scale_y_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2), labels = c("0.05", "0.1", "0.25", "0.5", "1", "2")) +
  coord_cartesian(ylim = c(0.05, 2.0)) +
  scale_fill_identity() +
  scale_colour_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = xlab5) +
  labs(x = "Five-level gene class", y = "LOEUF (gnomAD v4.1.1)") +
  theme(axis.text.x = element_text(size = 5.5, lineheight = 0.95),
        plot.margin = margin(t = 3, r = 2, b = 1, l = 1, unit = "mm"))

# --- assemble ---
right <- tg(pb, "b") / tg(pc, "c") + plot_layout(heights = c(0.72, 1))
fig <- tg(pa, "a") | right
fig <- fig + plot_layout(widths = c(1, 1.12))
ggsave(file.path(FIG_DIR, "ed4.png"), fig, width = 180, height = 112, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "ed4.pdf"), fig, width = 180, height = 112, units = "mm", device = cairo_pdf)
