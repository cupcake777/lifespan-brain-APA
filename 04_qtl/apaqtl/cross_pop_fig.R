#!/usr/bin/env Rscript
# Cross-population apaQTL tables and figures 1-4 (EAS postnatal vs EUR adult).
# Inputs from data/04_qtl; tables to data/04_qtl; figures to figures/04_qtl.
if (!nzchar(Sys.getenv("APA_ROOT")) || !nzchar(Sys.getenv("APA_HPC_HOME"))) stop("Set APA_ROOT and APA_HPC_HOME (source config/paths.sh)")

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(arrow)
  library(patchwork)
  library(scales)
  library(yaml)
  library(ggforce)
  library(qvalue)
})

# --- helpers ---

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

theme_crosspop <- function(base_size = base_font_size) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 1.2, color = ink_color),
      plot.subtitle = element_blank(),
      axis.line = element_line(linewidth = 0.38, color = ink_color),
      axis.ticks = element_line(linewidth = 0.32, color = ink_color),
      axis.title = element_text(face = "plain", color = ink_color),
      axis.text = element_text(color = ink_color),
      legend.title = element_text(face = "bold", color = ink_color),
      legend.text = element_text(color = ink_color),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(7, 8, 7, 8)
    )
}

save_figure <- function(plot_obj, basename, width = presentation_width, height = presentation_height) {
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))
  ggsave(pdf_path, plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot_obj, width = width, height = height, dpi = plot_dpi, device = ragg::agg_png, bg = "white")
  list(pdf = pdf_path, png = png_path)
}

read_sample_n <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  length(readLines(path, warn = FALSE))
}

count_unique <- function(x) uniqueN(x[!is.na(x)])

parse_phenotype_fields <- function(dt) {
  dt <- copy(dt)
  parts <- tstrsplit(dt$phenotype_id, "|", fixed = TRUE, fill = NA_character_)
  dt[, transcript_id := parts[[1]]]
  dt[, raw_gene_label := parts[[2]]]
  dt[, apa_locus := parts[[3]]]
  dt[, strand := parts[[4]]]
  dt[, feature_key := paste(transcript_id, apa_locus, strand, sep = "|")]
  dt[]
}

parse_variant_fields <- function(dt) {
  dt <- copy(dt)
  dt[, variant_key := as.character(variant_id)]
  parts <- tstrsplit(dt$variant_key, ":", fixed = TRUE, fill = NA_character_)
  dt[, variant_chr := parts[[1]]]
  dt[, variant_pos := suppressWarnings(as.integer(parts[[2]]))]
  dt[, variant_ref := parts[[3]]]
  dt[, variant_alt := parts[[4]]]
  dt[]
}

load_apa_map <- function(path, cohort) {
  dt <- fread(path)
  dt <- parse_phenotype_fields(dt)
  parts_gene <- tstrsplit(dt$Gene, "|", fixed = TRUE, fill = NA_character_)
  dt[, gene_id_from_gene_col := parts_gene[[2]]]
  dt[, cohort := cohort]
  unique(dt[, .(
    cohort,
    phenotype_id_map = phenotype_id,
    feature_key,
    transcript_id,
    gene_id_map = gene_id_from_gene_col,
    apa_locus,
    strand,
    Predicted_Proximal_APA
  )])
}

add_apa_map <- function(dt, apa_map) {
  dt <- merge(
    dt,
    unique(apa_map[, .(feature_key, gene_id_map, phenotype_id_map, Predicted_Proximal_APA)]),
    by = "feature_key",
    all.x = TRUE
  )
  dt[, gene_id := fifelse(!is.na(gene_id_map), gene_id_map, raw_gene_label)]
  dt[]
}

load_qtl_table <- function(path, cohort, type, apa_map) {
  dt <- fread(cmd = sprintf("gzip -cd %s", shQuote(path)))
  dt[, cohort := cohort]
  dt[, qtl_table_type := type]
  dt <- parse_phenotype_fields(dt)
  dt <- parse_variant_fields(dt)
  dt <- add_apa_map(dt, apa_map)
  dt[]
}

load_susie_table <- function(path, cohort, apa_map) {
  dt <- as.data.table(read_parquet(path))
  dt[, cohort := cohort]
  dt <- parse_phenotype_fields(dt)
  dt <- parse_variant_fields(dt)
  dt <- add_apa_map(dt, apa_map)
  dt[]
}

compute_directional_pair_pi1 <- function(discovery_dt, discovery_cohort, counterpart_cohort, counterpart_input) {
  key_path <- tempfile(pattern = "crosspop_pair_keys_", fileext = ".tsv")
  key_dt <- unique(discovery_dt[, .(phenotype_id, variant_id)])
  fwrite(key_dt, key_path, sep = "\t", col.names = FALSE)
  on.exit(unlink(key_path), add = TRUE)
  awk_script <- paste0(
    "awk 'BEGIN{FS=OFS=\"\\t\"} ",
    "FNR==NR{k[$1 FS $2]=1; next} ",
    "FNR==1{for(i=1;i<=NF;i++){",
    "if($i==\"phenotype_id\") pid=i; ",
    "if($i==\"variant_id\") vid=i; ",
    "if($i==\"pval_nominal\") pval=i; ",
    "if($i==\"slope\") slope=i}; ",
    "if(!pid || !vid || !pval || !slope){exit 2}; ",
    "print \"phenotype_id\",\"variant_id\",\"pval_nominal\",\"slope\"; next} ",
    "{key=$pid FS $vid; if(key in k){print $pid,$vid,$pval,$slope}}' ",
    shQuote(key_path),
    " -"
  )
  cmd <- paste("gzip -cd", shQuote(counterpart_input), "|", awk_script)
  matched <- fread(cmd = cmd)
  if (nrow(matched) == 0) {
    stop("No matched QTL pairs found for directional pi1 in ", counterpart_input)
  }
  matched <- parse_phenotype_fields(matched)
  matched <- parse_variant_fields(matched)
  matched[, pval_nominal := as.numeric(pval_nominal)]
  matched[, slope := as.numeric(slope)]
  pvals <- as.numeric(matched$pval_nominal)
  pvals <- pvals[is.finite(pvals) & pvals >= 0 & pvals <= 1]
  if (length(pvals) == 0) stop("No valid pval_nominal values found for directional pi1")
  pi0 <- qvalue::pi0est(pvals, lambda = seq(0.1, 0.9, 0.1), method = "smoother")$pi0
  pi0 <- max(0, min(1, pi0))
  pi <- list(pi0 = pi0, pi1 = 1 - pi0, n = length(pvals))
  discovery_effect <- unique(discovery_dt[, .(
    feature_key,
    variant_key,
    discovery_slope = slope,
    discovery_pval_nominal = pval_nominal
  )])
  matched_effect <- merge(
    matched[, .(feature_key, variant_key, counterpart_slope = slope, counterpart_pval_nominal = pval_nominal)],
    discovery_effect,
    by = c("feature_key", "variant_key"),
    all.x = TRUE
  )
  matched_effect <- matched_effect[is.finite(discovery_slope) & is.finite(counterpart_slope)]
  same_direction <- if (nrow(matched_effect) > 0) {
    mean(matched_effect$discovery_slope * matched_effect$counterpart_slope > 0, na.rm = TRUE)
  } else {
    NA_real_
  }
  data.table(
    cohort = discovery_cohort,
    counterpart_cohort = counterpart_cohort,
    pi1 = pi$pi1,
    pi0 = pi$pi0,
    discovery_qtl_pairs_for_pi1 = uniqueN(paste(discovery_dt$feature_key, discovery_dt$variant_key, sep = "||")),
    matched_qtl_pairs_for_pi1 = pi$n,
    matched_fraction_for_pi1 = pi$n / uniqueN(paste(discovery_dt$feature_key, discovery_dt$variant_key, sep = "||")),
    counterpart_p_lt_0_05_fraction = mean(matched$pval_nominal < 0.05, na.rm = TRUE),
    same_direction_fraction_for_pi1 = same_direction,
    pi1_input_file = counterpart_input,
    pi1_input_column = "pval_nominal",
    pi1_match_key = "phenotype_id + variant_id",
    pi1_estimand = paste0(
      cohort_plot_labels[discovery_cohort],
      " significant QTL pairs -> ",
      cohort_plot_labels[counterpart_cohort],
      " matched nominal QTL-pair p-values"
    )
  )
}

append_manifest <- function(figure, paths, source_rows, filter, note) {
  manifest <<- rbind(
    manifest,
    data.table(
      figure = figure,
      pdf = paths$pdf,
      png = paths$png,
      source_rows = as.integer(source_rows),
      filter = filter,
      note = note
    )
  )
}

summarize_numeric <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    list(min = NA_real_, median = NA_real_, max = NA_real_, n_unique = NA_integer_)
  } else {
    list(
      min = min(x),
      median = as.numeric(median(x)),
      max = max(x),
      n_unique = uniqueN(x)
    )
  }
}

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
qtl_dir <- file.path(apa_root, "data", "04_qtl")
table_dir <- qtl_dir
figure_dir <- file.path(apa_root, "figures", "04_qtl")
global_config_path <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/config.yml")

inputs <- list(
  eas_sig = file.path(qtl_dir, "postnatal.sig_QTL.txt.gz"),
  eur_sig = file.path(qtl_dir, "eur_adult.sig_pairs.txt.gz"),
  eas_perm = file.path(qtl_dir, "postnatal.cis_qtl.txt.gz"),
  eur_perm = file.path(qtl_dir, "BrainGvex", "permutation", "eur_adult.cis_qtl.txt.gz"),
  eas_nominal = file.path(qtl_dir, "nominal", "postnatal.cis_qtl_pairs.txt.gz"),
  eur_nominal = file.path(qtl_dir, "BrainGvex", "nominal", "eur_adult.cis_qtl_pairs.txt.gz"),
  eas_susie = file.path(qtl_dir, "susie", "postnatal.SuSiE_summary.parquet"),
  eur_susie = file.path(qtl_dir, "BrainGvex", "susie", "eur_adult.SuSiE_summary.parquet"),
  eas_apa = file.path(qtl_dir, "pre_input", "proximal.apa.info"),
  eur_apa = file.path(qtl_dir, "BrainGvex", "pre_input", "proximal.apa.info"),
  eas_samples = file.path(qtl_dir, "apa.sample.Postnatal.final"),
  eur_samples = file.path(apa_root, "data", "00_geno", "sample.Adult.final")
)

cohort_levels <- c("EAS postnatal", "EUR adult")
cohort_plot_labels <- c("EAS postnatal" = "EAS", "EUR adult" = "EUR")

cfg <- if (file.exists(global_config_path)) yaml::read_yaml(global_config_path) else list()
plot_cfg <- cfg$plot_settings %||% list()
fig_cfg <- plot_cfg$fig_size %||% list()
color_cfg <- cfg$colors %||% list()
palette_cfg <- color_cfg$palette %||% list()
base_cfg <- color_cfg$base %||% list()

plot_dpi <- plot_cfg$dpi %||% 320
font_family <- plot_cfg$font_family %||% "Arial"
base_font_size <- plot_cfg$base_font_size %||% 12
presentation_width <- fig_cfg$presentation_width %||% 6.0
presentation_height <- fig_cfg$presentation_height %||% 4.5

cool_color <- palette_cfg$cool %||% "#0072B2"
warm_color <- palette_cfg$warm %||% "#E69F00"
teal_color <- palette_cfg$teal %||% "#009E73"
gold_color <- palette_cfg$gold %||% "#F0E442"
shared_color <- teal_color
opposite_color <- palette_cfg$pink %||% "#CC79A7"
absent_color <- palette_cfg$absent %||% "#BDBDBD"
faint_color <- base_cfg$faint %||% "#F7F7F7"
line_color <- base_cfg$line %||% "#9E9E9E"
grid_color <- base_cfg$grid %||% "#E5E5E5"
ink_color <- base_cfg$ink %||% "#222222"

signal_colors <- c(
  "Significant QTL pairs" = teal_color,
  "Significant genes" = cool_color,
  "Significant variants" = warm_color
)
overlap_colors <- c(
  "EAS" = absent_color,
  "Shared" = shared_color,
  "EUR" = gold_color
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# --- manifest ---
manifest <- data.table(
  figure = character(),
  pdf = character(),
  png = character(),
  source_rows = integer(),
  filter = character(),
  note = character()
)

# --- load inputs ---
eas_apa <- load_apa_map(inputs$eas_apa, "EAS postnatal")
eur_apa <- load_apa_map(inputs$eur_apa, "EUR adult")
eas_sig <- load_qtl_table(inputs$eas_sig, "EAS postnatal", "significant_pairs", eas_apa)
eur_sig <- load_qtl_table(inputs$eur_sig, "EUR adult", "significant_pairs", eur_apa)
eas_perm <- load_qtl_table(inputs$eas_perm, "EAS postnatal", "permutation", eas_apa)
eur_perm <- load_qtl_table(inputs$eur_perm, "EUR adult", "permutation", eur_apa)
eas_susie <- load_susie_table(inputs$eas_susie, "EAS postnatal", eas_apa)
eur_susie <- load_susie_table(inputs$eur_susie, "EUR adult", eur_apa)

sample_summary <- data.table(
  cohort = cohort_levels,
  sample_size = c(read_sample_n(inputs$eas_samples), read_sample_n(inputs$eur_samples))
)

sig_all <- rbindlist(list(eas_sig, eur_sig), fill = TRUE)
perm_all <- rbindlist(list(eas_perm, eur_perm), fill = TRUE)
susie_all <- rbindlist(list(eas_susie, eur_susie), fill = TRUE)

# --- summaries and tables ---
discovery_summary <- sig_all[, .(
  significant_pairs = .N,
  significant_apa_features = count_unique(feature_key),
  significant_genes = count_unique(gene_id),
  significant_variants = count_unique(variant_key)
), by = cohort]

susie_summary <- susie_all[, .(
  susie_cs_rows = .N,
  susie_credible_sets = uniqueN(paste(feature_key, cs_id, sep = "||")),
  fine_mapped_variants = count_unique(variant_key),
  variants_pip_ge_0_5 = count_unique(variant_key[pip >= 0.5]),
  variants_pip_ge_0_8 = count_unique(variant_key[pip >= 0.8])
), by = cohort]

summary_table <- Reduce(
  function(x, y) merge(x, y, by = "cohort", all = TRUE),
  list(discovery_summary, susie_summary, sample_summary)
)
pi1_table <- rbindlist(list(
  compute_directional_pair_pi1(eas_sig, "EAS postnatal", "EUR adult", inputs$eur_nominal),
  compute_directional_pair_pi1(eur_sig, "EUR adult", "EAS postnatal", inputs$eas_nominal)
), fill = TRUE)
summary_table <- merge(summary_table, pi1_table, by = "cohort", all.x = TRUE)
summary_table[, cohort := factor(cohort, levels = cohort_levels)]
setorder(summary_table, cohort)

feature_counts <- sig_all[, .(sig_pair_count = .N), by = .(cohort, feature_key)]
feature_wide <- dcast(feature_counts, feature_key ~ cohort, value.var = "sig_pair_count", fill = 0)
setnames(feature_wide, make.names(names(feature_wide)))
feature_wide[, feature_overlap_class := fifelse(
  EAS.postnatal > 0 & EUR.adult > 0,
  "Shared",
  fifelse(
    EAS.postnatal > 0,
    "EAS",
    "EUR"
  )
)]

pair_counts <- sig_all[, .(sig_pair_count = .N), by = .(cohort, feature_key, variant_key)]
pair_wide <- dcast(pair_counts, feature_key + variant_key ~ cohort, value.var = "sig_pair_count", fill = 0)
setnames(pair_wide, make.names(names(pair_wide)))
pair_wide[, pair_overlap_class := fifelse(
  EAS.postnatal > 0 & EUR.adult > 0,
  "Shared",
  fifelse(EAS.postnatal > 0, "EAS", "EUR")
)]

pair_overlap <- pair_wide[, .(pair_count = .N), by = pair_overlap_class]
feature_overlap <- feature_wide[, .(feature_count = .N), by = feature_overlap_class]
overlap_table <- merge(
  data.table(overlap_class = names(overlap_colors)),
  pair_overlap,
  by.x = "overlap_class",
  by.y = "pair_overlap_class",
  all.x = TRUE
)
overlap_table <- merge(
  overlap_table,
  feature_overlap,
  by.x = "overlap_class",
  by.y = "feature_overlap_class",
  all.x = TRUE
)
overlap_table[is.na(pair_count), pair_count := 0L]
overlap_table[is.na(feature_count), feature_count := 0L]
overlap_table[, pair_fraction_of_all_significant_pairs := pair_count / sum(pair_count)]
overlap_table[, feature_fraction_of_all_significant_features := feature_count / sum(feature_count)]
overlap_table[, overlap_class := factor(overlap_class, levels = names(overlap_colors))]
setorder(overlap_table, overlap_class)

eas_effect <- eas_sig[, .(
  feature_key, variant_key, transcript_id, gene_id, apa_locus, strand,
  slope_eas = slope,
  slope_se_eas = slope_se,
  pval_nominal_eas = pval_nominal,
  af_eas = af,
  phenotype_id_eas = phenotype_id,
  variant_id = variant_id
)]
eur_effect <- eur_sig[, .(
  feature_key, variant_key,
  slope_eur = slope,
  slope_se_eur = slope_se,
  pval_nominal_eur = pval_nominal,
  af_eur = af,
  phenotype_id_eur = phenotype_id
)]

effect_concordance <- merge(eas_effect, eur_effect, by = c("feature_key", "variant_key"), all = FALSE)
if (nrow(effect_concordance) > 0) {
  effect_concordance[, direction := fifelse(
    slope_eas * slope_eur >= 0,
    "Same direction",
    "Opposite direction"
  )]
  effect_concordance[, abs_delta_af := abs(af_eas - af_eur)]
  effect_concordance[, effect_delta := slope_eas - slope_eur]
}

susie_all[, cs_key := paste(cohort, feature_key, cs_id, sep = "||")]
susie_cs_sizes <- susie_all[, .(
  cs_size = .N,
  max_pip = max(pip, na.rm = TRUE),
  mean_pip = mean(pip, na.rm = TRUE),
  variants_pip_ge_0_5 = sum(pip >= 0.5, na.rm = TRUE),
  variants_pip_ge_0_8 = sum(pip >= 0.8, na.rm = TRUE)
), by = .(cohort, feature_key, transcript_id, gene_id, apa_locus, strand, cs_id, cs_key)]

susie_pip <- susie_all[, .(
  cohort,
  feature_key,
  variant_key,
  cs_id,
  pip,
  af,
  transcript_id,
  gene_id,
  apa_locus,
  strand
)]

susie_feature_cs_counts <- susie_all[, .(
  cs_count_per_feature = uniqueN(cs_id)
), by = .(cohort, feature_key, gene_id)]
susie_feature_cs_counts[, cohort_label := cohort_plot_labels[cohort]]

susie_pip_wide <- dcast(
  susie_pip,
  feature_key + variant_key ~ cohort,
  value.var = "pip",
  fun.aggregate = max,
  fill = NA_real_
)
setnames(susie_pip_wide, make.names(names(susie_pip_wide)))
susie_pip_compare <- susie_pip_wide[
  is.finite(EAS.postnatal) & is.finite(EUR.adult),
  .(
    feature_key,
    variant_key,
    pip_eas = EAS.postnatal,
    pip_eur = EUR.adult
  )
]

susie_feature_cs_wide <- dcast(
  susie_feature_cs_counts,
  feature_key + gene_id ~ cohort,
  value.var = "cs_count_per_feature",
  fill = 0
)
setnames(susie_feature_cs_wide, make.names(names(susie_feature_cs_wide)))
susie_feature_cs_wide[, cs_count_eas := EAS.postnatal]
susie_feature_cs_wide[, cs_count_eur := EUR.adult]
susie_feature_cs_wide[, cs_bin_eas := fifelse(cs_count_eas >= 3, "3+", as.character(cs_count_eas))]
susie_feature_cs_wide[, cs_bin_eur := fifelse(cs_count_eur >= 3, "3+", as.character(cs_count_eur))]
cs_bin_levels <- c("0", "1", "2", "3+")
susie_feature_cs_wide[, cs_bin_eas := factor(cs_bin_eas, levels = cs_bin_levels)]
susie_feature_cs_wide[, cs_bin_eur := factor(cs_bin_eur, levels = cs_bin_levels)]
susie_feature_cs_matrix <- susie_feature_cs_wide[, .(
  feature_count = .N
), by = .(cs_bin_eas, cs_bin_eur)]
susie_feature_cs_matrix <- merge(
  CJ(cs_bin_eas = factor(cs_bin_levels, levels = cs_bin_levels), cs_bin_eur = factor(cs_bin_levels, levels = cs_bin_levels), unique = TRUE),
  susie_feature_cs_matrix,
  by = c("cs_bin_eas", "cs_bin_eur"),
  all.x = TRUE
)
susie_feature_cs_matrix[is.na(feature_count), feature_count := 0L]

pip_bins <- c(0.1, 0.5, 0.8, 0.95, Inf)
pip_labels <- c("0.1-0.5", "0.5-0.8", "0.8-0.95", ">0.95")
susie_pip_bins <- susie_all[pip > 0.1]
susie_pip_bins[, pip_bin := cut(pip, breaks = pip_bins, labels = pip_labels, right = FALSE, include.lowest = FALSE)]
pip_bin_summary <- susie_pip_bins[!is.na(pip_bin), .(
  variant_row_count = .N,
  unique_variant_count = count_unique(variant_key)
), by = .(cohort, pip_bin)]
pip_bin_grid <- CJ(cohort = cohort_levels, pip_bin = pip_labels, unique = TRUE)
pip_bin_summary <- merge(pip_bin_grid, pip_bin_summary, by = c("cohort", "pip_bin"), all.x = TRUE)
pip_bin_summary[is.na(variant_row_count), variant_row_count := 0L]
pip_bin_summary[is.na(unique_variant_count), unique_variant_count := 0L]
pip_bin_summary[, cohort_label := cohort_plot_labels[cohort]]
pip_bin_summary[, pip_bin := factor(pip_bin, levels = pip_labels)]

pip_threshold_summary <- rbindlist(lapply(c(0.1, 0.5, 0.8), function(thresh) {
  susie_all[pip > thresh, .(
    pip_threshold = paste0(">", thresh),
    variant_row_count = .N,
    unique_variant_count = count_unique(variant_key)
  ), by = cohort]
}), use.names = TRUE)
pip_threshold_summary[, cohort_label := cohort_plot_labels[cohort]]

cs_summary_table <- rbindlist(list(
  susie_all[, .(
    metric = "n_credible_sets",
    value = uniqueN(cs_key)
  ), by = cohort],
  susie_all[, .(
    metric = "n_features_with_cs",
    value = count_unique(feature_key)
  ), by = cohort],
  susie_all[, .(
    metric = "n_genes_with_cs",
    value = count_unique(gene_id)
  ), by = cohort],
  susie_feature_cs_counts[, .(
    metric = "cs_count_per_feature_min",
    value = as.numeric(min(cs_count_per_feature))
  ), by = cohort],
  susie_feature_cs_counts[, .(
    metric = "cs_count_per_feature_median",
    value = as.numeric(median(cs_count_per_feature))
  ), by = cohort],
  susie_feature_cs_counts[, .(
    metric = "cs_count_per_feature_max",
    value = as.numeric(max(cs_count_per_feature))
  ), by = cohort],
  pip_bin_summary[, .(
    metric = paste0("pip_bin_", as.character(pip_bin), "_rows"),
    value = as.numeric(variant_row_count)
  ), by = cohort],
  pip_threshold_summary[, .(
    metric = paste0("pip_", pip_threshold, "_rows"),
    value = as.numeric(variant_row_count)
  ), by = cohort],
  susie_cs_sizes[, .(
    metric = "true_cs_size_min",
    value = as.numeric(min(cs_size))
  ), by = cohort],
  susie_cs_sizes[, .(
    metric = "true_cs_size_median",
    value = as.numeric(median(cs_size))
  ), by = cohort],
  susie_cs_sizes[, .(
    metric = "true_cs_size_max",
    value = as.numeric(max(cs_size))
  ), by = cohort]
), fill = TRUE)

cs_count_stats <- summarize_numeric(susie_feature_cs_counts$cs_count_per_feature)
cs_size_stats <- summarize_numeric(susie_cs_sizes$cs_size)
pip_stats <- summarize_numeric(susie_all$pip)
susie_variable_audit <- rbindlist(list(
  data.table(
    variable_name = "cs_count_per_feature",
    definition_used = "number of credible sets assigned to one APA feature",
    source_file = paste(c(inputs$eas_susie, inputs$eur_susie), collapse = "; "),
    source_column = "derived from unique cs_id per cohort + feature_key",
    unit = "credible sets per APA feature",
    is_integer = TRUE,
    min = cs_count_stats$min,
    median = cs_count_stats$median,
    max = cs_count_stats$max,
    n_unique = cs_count_stats$n_unique,
    notes = "Fig4B APA feature level; not CS size."
  ),
  data.table(
    variable_name = "cs_size",
    definition_used = "number of variants in one credible set",
    source_file = paste(c(inputs$eas_susie, inputs$eur_susie), collapse = "; "),
    source_column = "derived from row count per cohort + feature_key + cs_id",
    unit = "variants per credible set",
    is_integer = TRUE,
    min = cs_size_stats$min,
    median = cs_size_stats$median,
    max = cs_size_stats$max,
    n_unique = cs_size_stats$n_unique,
    notes = "Fig4C credible-set level; true CS membership rows are available."
  ),
  data.table(
    variable_name = "pip",
    definition_used = "posterior inclusion probability for one variant in one credible set row",
    source_file = paste(c(inputs$eas_susie, inputs$eur_susie), collapse = "; "),
    source_column = "pip",
    unit = "probability",
    is_integer = FALSE,
    min = pip_stats$min,
    median = pip_stats$median,
    max = pip_stats$max,
    n_unique = pip_stats$n_unique,
    notes = "Fig4A feature-variant row level."
  )
), use.names = TRUE)

summary_table[, pi1_direction := pi1_estimand]
summary_table[, pi1_value := pi1]
summary_table[, pi1_input_columns := pi1_input_column]

fwrite(summary_table, file.path(table_dir, "01_discovery_sample_size_pi1.tsv"), sep = "\t")
fwrite(overlap_table, file.path(table_dir, "02_overlap.tsv"), sep = "\t")
fwrite(effect_concordance, file.path(table_dir, "03_effect.tsv"), sep = "\t")
fwrite(cs_summary_table, file.path(table_dir, "04_susie_cs_summary.tsv"), sep = "\t")
fwrite(susie_variable_audit, file.path(table_dir, "04_susie_variable_audit.tsv"), sep = "\t")
fwrite(susie_cs_sizes, file.path(table_dir, "04_susie_cs_sizes.tsv"), sep = "\t")
fwrite(susie_feature_cs_counts, file.path(table_dir, "04_susie_feature_cs_counts.tsv"), sep = "\t")
fwrite(susie_pip_compare, file.path(table_dir, "04_susie_pip_compare.tsv"), sep = "\t")

# --- figure 1: discovery landscape and directional pi1 ---

fig1_count_dt <- copy(summary_table)[, .(
  cohort,
  cohort_label = cohort_plot_labels[as.character(cohort)],
  sample_size,
  `Significant QTL pairs` = significant_pairs,
  `Significant genes` = significant_genes,
  `Significant variants` = significant_variants
)]
fig1_count_dt <- melt(
  fig1_count_dt,
  id.vars = c("cohort", "cohort_label", "sample_size"),
  variable.name = "signal_dimension",
  value.name = "signal_count"
)
fig1_count_dt[, signal_dimension := factor(signal_dimension, levels = names(signal_colors))]
fig1_count_dt[, cohort_x := match(cohort_label, c("EAS", "EUR"))]
fig1_count_dt[, x_plot := cohort_x + c(
  "Significant QTL pairs" = -0.055,
  "Significant genes" = 0,
  "Significant variants" = 0.055
)[as.character(signal_dimension)]]
fig1_count_dt[, label := comma(signal_count)]
fig1_count_dt[, label_y := signal_count * fifelse(signal_dimension == "Significant genes", 1.35, 1.16)]

fig1_pi1_dt <- copy(summary_table)[, .(
  cohort_label = cohort_plot_labels[as.character(cohort)],
  pi1_value = pi1
)]
fig1_pi1_dt[, direction_label := fifelse(cohort_label == "EAS", "EAS -> EUR", "EUR -> EAS")]
fig1_arrow_dt <- data.table(
  direction_label = c("EAS -> EUR", "EUR -> EAS"),
  x = c(1, 2),
  xend = c(2, 1),
  y = c(5.0e5, 4.2e5),
  yend = c(5.0e5, 4.2e5),
  label_x = c(1.40, 1.60),
  label_y = c(8.5e5, 5.8e5)
)
fig1_arrow_dt <- merge(fig1_arrow_dt, fig1_pi1_dt[, .(direction_label, pi1_value)], by = "direction_label", all.x = TRUE)

fig1 <- ggplot(fig1_count_dt, aes(x = x_plot, y = signal_count, color = signal_dimension, shape = signal_dimension)) +
  geom_curve(
    data = fig1_arrow_dt[direction_label == "EAS -> EUR"],
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    curvature = 0.22,
    arrow = arrow(length = unit(0.09, "inches"), type = "closed"),
    linetype = "dashed",
    linewidth = 0.48,
    color = shared_color,
    alpha = 0.75
  ) +
  geom_curve(
    data = fig1_arrow_dt[direction_label == "EUR -> EAS"],
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    curvature = -0.22,
    arrow = arrow(length = unit(0.09, "inches"), type = "closed"),
    linetype = "dashed",
    linewidth = 0.48,
    color = shared_color,
    alpha = 0.75
  ) +
  geom_text(
    data = fig1_arrow_dt,
    aes(x = label_x, y = label_y, label = paste0(direction_label, "  π1 = ", sprintf("%.3f", pi1_value))),
    inherit.aes = FALSE,
    color = ink_color,
    size = 3.0,
    fontface = "bold"
  ) +
  geom_point(size = 3.7, stroke = 0.7) +
  geom_text(
    aes(y = label_y, label = label),
    size = 2.85,
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(values = signal_colors, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      paste0("EAS\nn = ", comma(summary_table[cohort == "EAS postnatal", sample_size])),
      paste0("EUR\nn = ", comma(summary_table[cohort == "EUR adult", sample_size]))
    ),
    limits = c(0.68, 2.32)
  ) +
  scale_y_log10(
    breaks = c(1e3, 1e4, 1e5, 1e6),
    labels = label_number(scale_cut = cut_short_scale()),
    limits = c(1e3, 1.15e6),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = NULL, y = "Significant signals (log scale)", title = "Discovery landscape and pair-level sharing") +
  theme_crosspop(10.5) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10.5, lineheight = 0.9),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
fig1_paths <- save_figure(fig1, "01_sample_size_discovery_pi1", width = 7.2, height = 4.5)
append_manifest(
  "01_sample_size_discovery_pi1",
  fig1_paths,
  nrow(fig1_count_dt) + nrow(fig1_pi1_dt),
  "significant signal counts and directional cross-population nominal pval_nominal pi1",
  "Pair-level pi1 was computed from other-population nominal p-values matched to discovery significant QTL pairs."
)

# --- figure 2: observed pair-level overlap ---

fig2_dt <- copy(overlap_table)
fig2_dt[, overlap_class := factor(overlap_class, levels = names(overlap_colors))]
venn_circle_dt <- data.table(
  group = c("EAS", "EUR"),
  x0 = c(0, 1.02),
  y0 = c(0, 0),
  r = c(0.78, 0.78),
  fill = c(absent_color, gold_color)
)
venn_label_dt <- data.table(
  overlap_class = c("EAS", "Shared", "EUR"),
  x = c(-0.38, 0.51, 1.40),
  y = c(0.06, -0.10, 0.06)
)
venn_label_dt <- merge(
  venn_label_dt,
  fig2_dt[, .(
    overlap_class = as.character(overlap_class),
    pair_count,
    pair_fraction_of_all_significant_pairs
  )],
  by = "overlap_class",
  all.x = TRUE
)
venn_label_dt[, label := paste0(
  overlap_class,
  "\n",
  comma(pair_count),
  "\n",
  percent(pair_fraction_of_all_significant_pairs, accuracy = 0.1)
)]

venn_panel <- ggplot() +
  geom_circle(
    data = venn_circle_dt,
    aes(x0 = x0, y0 = y0, r = r, fill = group),
    color = ink_color,
    linewidth = 0.4,
    alpha = 0.30
  ) +
  geom_text(
    data = venn_label_dt,
    aes(x = x, y = y, label = label),
    size = 3.15,
    fontface = "bold",
    lineheight = 0.90
  ) +
  scale_fill_manual(values = c("EAS" = absent_color, "EUR" = gold_color), guide = "none") +
  coord_equal(xlim = c(-0.90, 1.92), ylim = c(-0.87, 0.87), clip = "off") +
  labs(title = "Exact significant pairs") +
  theme_void(base_family = font_family) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5, color = ink_color),
    plot.margin = margin(7, 7, 7, 7)
  )

percent_panel <- ggplot(fig2_dt, aes(x = overlap_class, y = pair_fraction_of_all_significant_pairs, fill = overlap_class)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = percent(pair_fraction_of_all_significant_pairs, accuracy = 0.1)),
    vjust = -0.35,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = overlap_colors, guide = "none") +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 0.56),
    breaks = seq(0, 0.5, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Share of union", title = "Composition of the union") +
  theme_crosspop(10.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = grid_color, linewidth = 0.25)
  )

fig2 <- (venn_panel | percent_panel) +
  plot_layout(widths = c(1.05, 0.95)) +
  plot_annotation(
    title = "Observed pair-level overlap",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5, color = ink_color))
  )
fig2_paths <- save_figure(fig2, "02_overlap", width = 7.6, height = 3.9)
append_manifest(
  "02_overlap",
  fig2_paths,
  nrow(overlap_table),
  "exact significant feature_key+variant_key pair categories; feature counts retained in table",
  "Figure shows pair-level observed overlap; feature-level counts are retained in 02_overlap.tsv."
)

# --- figure 3: effect concordance ---

if (nrow(effect_concordance) >= 3) {
  fig3_rho <- cor(effect_concordance$slope_eas, effect_concordance$slope_eur, method = "spearman", use = "complete.obs")
  same_rate <- mean(effect_concordance$direction == "Same direction", na.rm = TRUE)
  opposite_n <- sum(effect_concordance$direction == "Opposite direction", na.rm = TRUE)
  lims <- range(c(effect_concordance$slope_eas, effect_concordance$slope_eur), na.rm = TRUE)
  pad <- diff(lims) * 0.08
  if (!is.finite(pad) || pad == 0) pad <- 0.1
  lims <- c(lims[1] - pad, lims[2] + pad)
  fig3_note <- paste0(
    "Spearman rho = ", sprintf("%.2f", fig3_rho),
    "\n", percent(same_rate, accuracy = 0.1), " same direction",
    "\n", comma(nrow(effect_concordance)), " shared pairs"
  )
  fig3 <- ggplot() +
    geom_bin_2d(
      data = effect_concordance[direction == "Same direction"],
      aes(x = slope_eas, y = slope_eur, fill = after_stat(count)),
      bins = 58,
      alpha = 0.96
    ) +
    geom_point(
      data = effect_concordance[direction == "Opposite direction"],
      aes(x = slope_eas, y = slope_eur, color = "Opposite direction"),
      size = 1.15,
      alpha = 0.62
    ) +
    geom_hline(yintercept = 0, color = line_color, linewidth = 0.32) +
    geom_vline(xintercept = 0, color = line_color, linewidth = 0.32) +
    geom_abline(intercept = 0, slope = 1, color = ink_color, linetype = "dashed", linewidth = 0.48) +
    annotate(
      "label",
      x = lims[1] + diff(lims) * 0.035,
      y = lims[2] - diff(lims) * 0.035,
      hjust = 0,
      vjust = 1,
      label = fig3_note,
      size = 3.15,
      fontface = "bold",
      label.size = 0,
      fill = "white",
      color = ink_color,
      lineheight = 1.05
    ) +
    coord_equal(xlim = lims, ylim = lims) +
    scale_fill_gradient(
      low = faint_color,
      high = teal_color,
      trans = "sqrt",
      name = "Same-direction\npairs per bin"
    ) +
    scale_color_manual(values = c("Opposite direction" = opposite_color), name = NULL) +
    guides(
      fill = guide_colorbar(order = 1, title.position = "top", barwidth = unit(1.45, "inches"), barheight = unit(0.10, "inches")),
      color = guide_legend(order = 2, override.aes = list(size = 2.4, alpha = 1))
    ) +
    labs(
      x = "EAS slope",
      y = "EUR slope",
      title = "Shared apaQTL effects are strongly concordant"
    ) +
    theme_crosspop(10.5) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
} else {
  fig3 <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Fewer than 3 exact shared significant apaQTL pairs", size = 4) +
    theme_void(base_family = font_family)
}
fig3_paths <- save_figure(fig3, "03_effect", width = 6.4, height = 5.4)
append_manifest(
  "03_effect",
  fig3_paths,
  nrow(effect_concordance),
  "exact shared significant feature_key+variant_key pairs",
  "Effect concordance only."
)

# --- figure 4: fine-mapping agreement and resolution ---

susie_pip_compare[, high_pip_class := fifelse(
  pip_eas >= 0.8 & pip_eur >= 0.8,
  "Both cohorts PIP >= 0.8",
  fifelse(
    pip_eas >= 0.8 | pip_eur >= 0.8,
    "One cohort PIP >= 0.8",
    "Below 0.8 in both"
  )
)]
susie_pip_compare[, high_pip_class := factor(
  high_pip_class,
  levels = c("Below 0.8 in both", "One cohort PIP >= 0.8", "Both cohorts PIP >= 0.8")
)]
pip_high_either <- susie_pip_compare[pmax(pip_eas, pip_eur) >= 0.8, .N]
pip_high_both <- susie_pip_compare[pip_eas >= 0.8 & pip_eur >= 0.8, .N]

panel_a <- ggplot() +
  geom_point(
    data = susie_pip_compare[high_pip_class == "Below 0.8 in both"],
    aes(x = pip_eas, y = pip_eur),
    color = absent_color,
    size = 0.72,
    alpha = 0.34
  ) +
  geom_point(
    data = susie_pip_compare[high_pip_class != "Below 0.8 in both"],
    aes(x = pip_eas, y = pip_eur, fill = high_pip_class),
    shape = 21,
    color = ink_color,
    stroke = 0.20,
    size = 2.15,
    alpha = 0.90
  ) +
  geom_abline(intercept = 0, slope = 1, color = ink_color, linetype = "dashed", linewidth = 0.42) +
  geom_vline(xintercept = 0.8, color = line_color, linetype = "dotted", linewidth = 0.4) +
  geom_hline(yintercept = 0.8, color = line_color, linetype = "dotted", linewidth = 0.4) +
  annotate(
    "label",
    x = 0.04,
    y = 0.97,
    hjust = 0,
    vjust = 1,
    label = paste0("Both high: ", pip_high_both, " / ", pip_high_either, " exact rows high in either"),
    size = 2.7,
    fontface = "bold",
    label.size = 0,
    fill = "white",
    color = ink_color
  ) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = number_format(accuracy = 0.1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(labels = number_format(accuracy = 0.1), breaks = seq(0, 1, by = 0.25)) +
  scale_fill_manual(
    values = c(
      "One cohort PIP >= 0.8" = opposite_color,
      "Both cohorts PIP >= 0.8" = teal_color
    ),
    name = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(size = 2.8))) +
  labs(x = "PIP in EAS", y = "PIP in EUR", title = "Feature-variant PIP") +
  theme_crosspop(9.5) +
  theme(
    plot.title = element_text(size = 10.5, hjust = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 7.6),
    plot.margin = margin(8, 8, 6, 8)
  )

panel_b <- ggplot(susie_feature_cs_matrix, aes(x = cs_bin_eas, y = cs_bin_eur, fill = feature_count)) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(aes(label = ifelse(feature_count > 0, comma(feature_count), "")), size = 2.8, fontface = "bold") +
  scale_fill_gradient(low = faint_color, high = teal_color, trans = "sqrt", labels = comma, name = "Features") +
  labs(x = "CS count in EAS", y = "CS count in EUR", title = "CS count per APA feature") +
  coord_fixed() +
  theme_crosspop(9.2) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 10.5, hjust = 0),
    plot.margin = margin(8, 6, 4, 8)
  )

susie_cs_sizes[, cohort_label := cohort_plot_labels[cohort]]
cs_median_dt <- susie_cs_sizes[, .(median_cs_size = median(cs_size)), by = cohort_label]
panel_c <- ggplot(susie_cs_sizes, aes(x = cs_size, linetype = cohort_label)) +
  stat_ecdf(linewidth = 0.82, color = ink_color) +
  scale_x_continuous(
    labels = comma_format(accuracy = 1),
    trans = "log10",
    breaks = c(1, 3, 10, 30, 100, 300, 1000)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0.01, 0.02))) +
  scale_linetype_manual(values = c("EAS" = "solid", "EUR" = "longdash"), name = NULL) +
  guides(linetype = guide_legend(nrow = 1, byrow = TRUE)) +
  annotate(
    "text",
    x = 1.1,
    y = 0.97,
    hjust = 0,
    vjust = 1,
    label = paste0(
      "Median: EAS ", cs_median_dt[cohort_label == "EAS", median_cs_size],
      "; EUR ", cs_median_dt[cohort_label == "EUR", median_cs_size]
    ),
    size = 2.6,
    color = ink_color
  ) +
  labs(x = "Variants per credible set", y = "Cumulative fraction", title = "Credible-set resolution") +
  theme_crosspop(9.2) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(0.46, "inches"),
    plot.title = element_text(size = 10.5, hjust = 0),
    panel.grid.major.y = element_line(color = grid_color, linewidth = 0.24),
    plot.margin = margin(4, 6, 5, 8)
  )

fig4 <- (panel_a | (panel_b / panel_c + plot_layout(heights = c(0.95, 1.05)))) +
  plot_layout(widths = c(1.22, 1.0)) +
  plot_annotation(
    title = "Fine-mapping agreement and resolution",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13.5, color = ink_color),
      plot.tag = element_text(face = "bold", size = 11, color = ink_color)
    )
  ) &
  theme(
    plot.tag.position = c(0.01, 0.99)
  )
fig4_paths <- save_figure(fig4, "04_susie_cs", width = 8.6, height = 5.8)
append_manifest(
  "04_susie_cs",
  fig4_paths,
  nrow(susie_pip_compare),
  "SuSiE rows with phenotype_id, variant_id, cs_id, and pip",
  "Panels show feature-variant PIP, APA feature CS counts, and true variants per credible set."
)

# --- layout validation ---

layout_validation <- data.table(
  figure = c(
    "01_sample_size_discovery_pi1",
    "02_overlap",
    "03_effect",
    "04_susie_cs"
  ),
  width_in = c(7.2, 7.6, 6.4, 8.6),
  height_in = c(4.5, 3.9, 5.4, 5.8),
  dpi = plot_dpi,
  png = c(fig1_paths$png, fig2_paths$png, fig3_paths$png, fig4_paths$png),
  pdf = c(fig1_paths$pdf, fig2_paths$pdf, fig3_paths$pdf, fig4_paths$pdf),
  required_content = c(
    "log-scale signal counts; cohort n; directional pi1 labels",
    "Venn counts; union percentages; EAS/Shared/EUR labels",
    "effect axes; rho; same-direction fraction; density and opposite-direction legend",
    "panel tags; PIP threshold; CS-count matrix; CS-size ECDF"
  ),
  visual_overlap = "pending rendered inspection",
  visual_clipping = "pending rendered inspection",
  final_size_readability = "pending rendered inspection"
)
layout_validation[, expected_width_px := round(width_in * dpi)]
layout_validation[, expected_height_px := round(height_in * dpi)]
layout_validation[, png_bytes := file.info(png)$size]
layout_validation[, pdf_bytes := file.info(pdf)$size]
layout_validation[, export_status := fifelse(png_bytes > 0 & pdf_bytes > 0, "pass", "fail")]
fwrite(layout_validation, file.path(table_dir, "crosspop_layout_validation.tsv"), sep = "\t")

fwrite(manifest, file.path(table_dir, "crosspop_figure_manifest.tsv"), sep = "\t")
