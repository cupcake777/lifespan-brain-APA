# S-LDSC summary: tau* per annotation, FDR within population x stage x trait, and four figures.
# Reads data/04_qtl/sldsc/ (logs, results) and annot_sd.tsv; writes tables to data/04_qtl and figures to figures/04_qtl.
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# --- paths ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
sldsc_dir <- file.path(D, "sldsc")
fig_dir <- file.path(APA_ROOT, "figures", "04_qtl")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- helpers ---
save_panel <- function(p, stem, w = 5.5, h = 4.8) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 320, bg = "white")
}

theme_sldsc <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(colour = "black", size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1)
    )
}

trait_dict <- data.table(
  code = c("SCZ", "BIP", "MDD", "ADHD", "ASD", "ANX", "PTSD", "INT", "INS", "AD", "PD", "ALS", "EP", "AIS"),
  label = c("Schizophrenia", "Bipolar disorder", "Major depressive disorder", "ADHD", "Autism spectrum disorder",
            "Anxiety disorder", "PTSD", "Intelligence", "Insomnia", "Alzheimer's disease", "Parkinson's disease",
            "Amyotrophic lateral sclerosis", "Epilepsy", "Ischemic stroke")
)

# --- heritability per trait (LDSC log: Total Observed scale h2) ---
h2_dt <- rbindlist(lapply(list.files(file.path(sldsc_dir, "logs"), pattern = "^h2\\.maxcpp\\..*\\.log$", full.names = TRUE), function(f) {
  p <- strsplit(sub("\\.log$", "", basename(f)), ".", fixed = TRUE)[[1]]
  hit <- grep("^Total Observed scale h2:", readLines(f, warn = FALSE), value = TRUE)
  if (!length(hit)) return(NULL)
  m <- regmatches(hit[1], regexec("Total Observed scale h2:\\s*([0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", hit[1]))[[1]]
  if (length(m) < 3L) return(NULL)
  data.table(population = p[3], stage = p[4], trait = p[5], h2 = as.numeric(m[2]), h2_se = as.numeric(m[3]))
}))
h2_dt <- unique(h2_dt, by = c("population", "stage", "trait"))
h2_dt[, well_powered := (h2 / h2_se) >= 7.0]

# --- regression results: 38 trait x stage files = 76 annotation tests ---
result_files <- list.files(file.path(sldsc_dir, "results", "maxcpp"), pattern = "\\.results$", recursive = TRUE, full.names = TRUE)
stopifnot(length(result_files) == 38L)
raw <- rbindlist(lapply(result_files, function(f) {
  data.table(model = "maxcpp", population = basename(dirname(dirname(f))), stage = basename(dirname(f)),
             trait = sub("\\.results$", "", basename(f)), fread(f))
}))
legacy <- fread(file.path(D, "sldsc_qtl_results.tsv"))[model != "maxcpp"]
combined <- rbindlist(list(legacy, raw), use.names = TRUE)
fwrite(combined, file.path(D, "sldsc_qtl_results.tsv"), sep = "\t")
fwrite(combined[grepl("^(apaQTL|eQTL)L2_[0-9]+$", Category)], file.path(D, "sldsc_qtl_annotations.tsv"), sep = "\t")

setnames(raw, "Coefficient_z-score", "coef_z", skip_absent = TRUE)
qtl <- raw[grepl("^(apaQTL|eQTL)L2_[0-9]+$", Category) & model == "maxcpp"]
qtl[, annotation := fifelse(grepl("^apaQTL", Category), "apaQTL", "eQTL")]
qtl[, coef_p_1s := pnorm(coef_z, lower.tail = FALSE)]
qtl[, coef_fdr := p.adjust(coef_p_1s, method = "BH"), by = .(population, stage, trait)]
qtl[, n_tests_family := .N, by = .(population, stage, trait)]
stopifnot(nrow(qtl) == 76L, all(qtl$n_tests_family == 2L), !anyDuplicated(qtl[, .(population, stage, trait, annotation)]))

sd_dt <- unique(fread(file.path(D, "annot_sd.tsv"))[annot_set == "maxCPP", .(population, stage, annotation, sd_annot = sd, M_common = n_snps)],
                by = c("population", "stage", "annotation"))
qtl <- merge(qtl, h2_dt, by = c("population", "stage", "trait"), all.x = TRUE)
qtl <- merge(qtl, sd_dt, by = c("population", "stage", "annotation"), all.x = TRUE)

# --- standardised tau* and 95% CI ---
qtl[, tau_star := fifelse(h2 > 0, Coefficient * M_common * sd_annot / h2, NA_real_)]
qtl[, tau_star_se := fifelse(h2 > 0, Coefficient_std_error * M_common * sd_annot / h2, NA_real_)]
qtl[, `:=`(tau_lo = tau_star - 1.96 * tau_star_se, tau_hi = tau_star + 1.96 * tau_star_se)]
fwrite(qtl[order(population, stage, annotation, trait)], file.path(D, "03_sldsc_tau_star.tsv"), sep = "\t")

qtl[, code := sub("^(EAS|EUR)_", "", trait)]
qtl <- merge(qtl, trait_dict, by = "code", all.x = TRUE)
qtl[is.na(label), label := code]
qtl[, stage := factor(stage, levels = c("prenatal", "postnatal"), labels = c("Prenatal", "Postnatal"))]
qtl[, label := factor(code, levels = rev(trait_dict$code))]
dodge_forest <- position_dodge(width = 0.6)

# --- figure 1: tau* forest ---
p_forest <- ggplot(qtl[(well_powered)], aes(x = tau_star, y = label, colour = annotation, group = annotation)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_linerange(aes(xmin = tau_lo, xmax = tau_hi), position = dodge_forest, linewidth = 0.6) +
  geom_point(aes(shape = coef_fdr < 0.05), position = dodge_forest, size = 2.2) +
  facet_grid(population ~ stage, scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c(apaQTL = "#0F6FA8", eQTL = "#C46A28"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 21), labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(x = expression(paste("Standardised coefficient ", tau, "*")), y = NULL) +
  theme_sldsc() + theme(legend.position = "top", strip.background = element_rect(fill = "grey95", colour = NA))
save_panel(p_forest, "03-1_forest", w = 6.5, h = 5.2)

# --- figure 2: apaQTL vs eQTL z-score ---
scat_dt <- unique(qtl[(well_powered), .(trait, code, label, stage, population, annotation, coef_z, coef_p_1s, coef_fdr)],
                  by = c("trait", "stage", "population", "annotation"))
scat <- dcast(scat_dt, trait + code + label + stage + population ~ annotation, value.var = c("coef_z", "coef_p_1s", "coef_fdr"))
setnames(scat, c("coef_z_apaQTL", "coef_z_eQTL", "coef_p_1s_apaQTL", "coef_p_1s_eQTL", "coef_fdr_apaQTL", "coef_fdr_eQTL"),
         c("z_apa", "z_eqtl", "p_apa", "p_eqtl", "fdr_apa", "fdr_eqtl"))
scat[, call := fifelse(fdr_apa < 0.05 & fdr_eqtl >= 0.05, "apaQTL only",
                fifelse(fdr_apa >= 0.05 & fdr_eqtl < 0.05, "eQTL only",
                fifelse(fdr_apa < 0.05 & fdr_eqtl < 0.05, "Both", "Neither")))]
scat[, call := factor(call, levels = c("apaQTL only", "Both", "eQTL only", "Neither"))]
scat[, plot_label := sprintf("%s-%s", population, code)]
lab <- scat[call != "Neither"]
lab[, nudge_x := fifelse(z_eqtl < 1.645, -0.35, 0.35)]
lab[, nudge_y := fifelse(z_apa < 1.645, -0.35, 0.35)]
lim <- max(abs(c(scat$z_apa, scat$z_eqtl)), na.rm = TRUE) * 1.3
p_scat <- ggplot(scat, aes(x = z_eqtl, y = z_apa)) +
  geom_hline(yintercept = 1.645, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 1.645, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, colour = "grey80") +
  geom_point(aes(fill = call, shape = stage), size = 3.0, colour = "grey20", stroke = 0.4) +
  geom_text_repel(data = lab, aes(label = plot_label), size = 2.8, colour = "grey15",
                  nudge_x = lab$nudge_x, nudge_y = lab$nudge_y, direction = "both",
                  box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf, seed = 123) +
  scale_fill_manual(values = c("apaQTL only" = "#0F6FA8", "Both" = "#4F7942", "eQTL only" = "#C46A28", "Neither" = "grey85"), name = "Category") +
  scale_shape_manual(values = c(Prenatal = 21, Postnatal = 24), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(fill = "grey60", size = 3))) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "eQTL coefficient z-score", y = "apaQTL coefficient z-score") +
  theme_sldsc() + theme(legend.position = "right")
save_panel(p_scat, "03-2_zscore", w = 6.6, h = 5.2)

# --- figure 3: GWAS h2 z-score QC ---
qc_h2 <- unique(qtl[, .(code, label, population, stage, h2_z = h2 / h2_se, well_powered)], by = c("code", "population", "stage"))
p_h2_faceted <- ggplot(qc_h2, aes(x = h2_z, y = reorder(label, h2_z), fill = well_powered)) +
  geom_vline(xintercept = 7.0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_col(width = 0.65, colour = "grey30", linewidth = 0.3) +
  facet_grid(population ~ stage, scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#3E6E8E", `FALSE` = "grey82"), labels = c(`TRUE` = "Retained", `FALSE` = "Underpowered"), name = NULL) +
  labs(x = expression(paste("GWAS Total ", h^2, " z-score")), y = NULL) +
  theme_sldsc() + theme(legend.position = "top", strip.background = element_rect(fill = "grey95", colour = NA))
save_panel(p_h2_faceted, "03-3_h2", w = 7.0, h = 6.0)

# --- figure 4: apaQTL tau* heatmap ---
heat_dt <- qtl[annotation == "apaQTL", .(population, stage, label, tau_star, coef_p_1s, coef_fdr, well_powered)]
heat_dt[, usable := (well_powered) & !is.na(tau_star)]
heat_dt[, sig := fifelse(!usable, "", paste0(fifelse(coef_fdr < 0.05, "*", ""),
                                             fifelse(coef_p_1s < 0.05 & coef_fdr >= 0.05, "\u2020", "")))]
heat_dt[, cell_txt := fifelse(is.na(tau_star), "\u2014", sprintf("%.1f%s", tau_star, sig))]
heat_dt[usable == FALSE, tau_star := NA_real_]
heat_dt[, population := factor(population, levels = c("EUR", "EAS"))]
p_heat <- ggplot(heat_dt, aes(x = stage, y = label)) +
  geom_tile(aes(fill = tau_star), colour = "white", linewidth = 0.6) +
  geom_text(aes(label = cell_txt, colour = usable), size = 2.8) +
  facet_grid(. ~ population, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-3, 3),
                       oob = scales::squish, na.value = "grey88",
                       name = expression(paste("Standardised coefficient ", tau, "*"))) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey55"), guide = "none") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = "* FDR < 0.05; \u2020 nominal P < 0.05; grey: h\u00b2 z < 7 or h\u00b2 \u2264 0") +
  theme_sldsc() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        plot.caption = element_text(size = 8, colour = "grey30"))
save_panel(p_heat, "03-4_heatmap", w = 7.2, h = 5.4)
