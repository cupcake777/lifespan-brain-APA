#!/usr/bin/env Rscript
# Cohort overview: donor-weighted age histograms and counts for trajectory and matched-protein cohorts.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
setDTthreads(1L)

# --- paths and settings ---
root <- Sys.getenv("APA_ROOT")
trajectory_path <- file.path(root, "data/02_traj/00_cov_profiles.tsv")
protein_path <- file.path(root, "data/03_omics/matched_sample_manifest.tsv")
out_tsv <- file.path(root, "data/01_quant/cohort_overview.tsv")
out_fig <- file.path(root, "figures/01_quant/cohort_overview")
age_ticks <- c(-0.5, 0, 1, 5, 12, 40, 100)
figure_size <- c(7.2, 4.6)
stage_colors <- c(Prenatal = "#0072B2", Postnatal = "#D55E00")
dir.create(dirname(out_fig), recursive = TRUE, showWarnings = FALSE)

# --- load the two analysis cohorts ---
trajectory <- fread(trajectory_path)
protein <- fread(protein_path)
age_data <- rbindlist(list(
  trajectory[, .(cohort = "Trajectory: 815 profiles / 750 individuals", sample_id = sample,
                 individual_id, age, weight = profile_weight)],
  protein[, .(cohort = "Matched proteomics: 366 donors", sample_id,
              individual_id = canonical_id, age = age_years, weight = 1)]))
age_data[, `:=`(age_model = log10(age + 1),
                life_stage = factor(ifelse(age < 0, "Prenatal", "Postnatal"), levels = names(stage_colors)))]
age_data[, cohort := factor(cohort, levels = unique(cohort))]

summary_dt <- age_data[, .(n_profiles = .N, n_individuals = uniqueN(individual_id),
                           total_donor_weight = sum(weight), n_prenatal = sum(age < 0), n_postnatal = sum(age >= 0),
                           min_age = min(age), max_age = max(age),
                           n_profiles_age_1_to_10 = sum(age >= 1 & age <= 10),
                           donor_weight_age_1_to_10 = sum(weight[age >= 1 & age <= 10])), by = cohort]
stopifnot(identical(summary_dt$n_profiles, c(815L, 366L)),
          all(abs(summary_dt$total_donor_weight - c(750, 366)) < 1e-8))
fwrite(summary_dt, paste0(out_tsv), sep = "\t")

# --- age histograms on a shared log10(age + 1) grid ---
binwidth <- diff(range(age_data$age_model)) / 50
p_age <- ggplot(age_data, aes(x = age_model, weight = weight, fill = life_stage)) +
  geom_histogram(binwidth = binwidth, boundary = 0, color = "white", linewidth = 0.15) +
  facet_wrap(~ cohort, ncol = 1) +
  scale_x_continuous(breaks = log10(age_ticks + 1), labels = age_ticks,
                     expand = expansion(mult = c(0.015, 0.015))) +
  scale_fill_manual(values = stage_colors, drop = FALSE) +
  labs(x = "Age (years; log10[age + 1] spacing)", y = "Donor-weighted count", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
        strip.background = element_blank(), strip.text = element_text(hjust = 0),
        panel.spacing = grid::unit(0.65, "lines"), plot.margin = margin(5, 10, 5, 5))
for (extension in c("png", "pdf")) {
  ggsave(paste0(out_fig, ".", extension), p_age, width = figure_size[1], height = figure_size[2], dpi = 300,
         device = if (extension == "pdf") grDevices::cairo_pdf else "png", bg = "white")
}
print(summary_dt)
