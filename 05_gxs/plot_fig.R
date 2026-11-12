# Stage-specific interaction figures (panels A-D, combined main figure) and summary tables.
# Reads data/05_gxs/sig.interaction.result and the stage, phenotype and genotype inputs; writes figures/05_gxs and data/05_gxs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

# --- helpers ---
msg <- function(...) cat(sprintf(...), "\n")

pair_key <- function(phenotype_id, variant_id) {
  paste(phenotype_id, variant_id, sep = "||")
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

save_plot <- function(plot_obj, filename, out_dir, width, height) {
  stem <- basename(filename)
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 320, bg = "white")
  msg("saved: %s", png_file)
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  ggsave(pdf_file, plot = plot_obj, width = width, height = height, device = cairo_pdf, bg = "white")
  msg("saved: %s", pdf_file)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_file <- file.path(out_dir, paste0(stem, ".svg"))
    ggsave(svg_file, plot = plot_obj, width = width, height = height, device = "svg", bg = "white")
    msg("saved: %s", svg_file)
  } else {
    msg("skip svg: svglite package is not installed")
  }
}

theme_main_minimal <- function(base_size = 7.6) {
  theme_classic(base_size = base_size, base_family = "") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, color = base_colors[["ink"]]),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(face = "bold", color = base_colors[["ink"]]),
      axis.text = element_text(color = base_colors[["ink"]]),
      axis.line = element_line(color = base_colors[["ink"]], linewidth = 0.35),
      axis.ticks = element_line(color = base_colors[["ink"]], linewidth = 0.3),
      legend.title = element_text(face = "bold", size = base_size * 0.85),
      legend.text = element_text(size = base_size * 0.85),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = base_colors[["ink"]])
    )
}

add_panel_tag <- function(plot_obj, tag) {
  plot_obj +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(face = "bold", size = 11.5, color = base_colors[["ink"]]),
      plot.tag.position = c(0, 1)
    )
}

format_variant_label <- function(variant_id) {
  parts <- strsplit(variant_id, ":", fixed = TRUE)[[1]]
  if (length(parts) >= 4) {
    return(sprintf("%s:%s %s>%s", parts[1], parts[2], parts[3], parts[4]))
  }
  variant_id
}

read_signif_pairs <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }
  fread(cmd = sprintf("gzip -dc %s", shQuote(file)), select = c("phenotype_id", "variant_id"))
}

read_lifestage <- function(file) {
  if (!file.exists(file)) {
    return(NULL)
  }
  stage_dt <- fread(file)
  if (!all(c("sample", "LifeStage") %in% names(stage_dt))) {
    stop("LifeStage file must contain sample and LifeStage columns")
  }
  stage_dt[, stage := fifelse(as.character(LifeStage) == "0", "Prenatal", "Postnatal")]
  stage_dt[, stage := factor(stage, levels = c("Prenatal", "Postnatal"))]
  stage_dt[, .(sample, stage)]
}

read_phenotype_subset <- function(file, phenotype_ids) {
  if (!file.exists(file)) {
    return(NULL)
  }
  phenotype_ids <- unique(phenotype_ids)
  pheno_dt <- fread(file)
  id_col <- if ("pid" %in% names(pheno_dt)) "pid" else names(pheno_dt)[4]
  pheno_dt <- pheno_dt[get(id_col) %in% phenotype_ids]
  if (nrow(pheno_dt) == 0) {
    return(NULL)
  }
  sample_cols <- setdiff(names(pheno_dt), names(pheno_dt)[1:4])
  long_dt <- melt(
    pheno_dt,
    id.vars = id_col,
    measure.vars = sample_cols,
    variable.name = "sample",
    value.name = "phenotype"
  )
  setnames(long_dt, id_col, "phenotype_id")
  long_dt[, phenotype := as.numeric(phenotype)]
  long_dt
}

load_genotype_subset <- function(pfile_prefix, variant_ids) {
  variant_ids <- unique(variant_ids)
  if (length(variant_ids) == 0 || Sys.which("plink2") == "") {
    return(NULL)
  }
  expected_files <- paste0(pfile_prefix, c(".pgen", ".pvar", ".psam"))
  if (any(!file.exists(expected_files))) {
    return(NULL)
  }
  tmp_snplist <- tempfile("example_snps_")
  tmp_out <- tempfile("example_geno_")
  writeLines(variant_ids, tmp_snplist)
  system2(
    "plink2",
    c("--pfile", pfile_prefix, "--extract", tmp_snplist, "--export", "A", "--out", tmp_out),
    stdout = TRUE,
    stderr = TRUE
  )
  raw_file <- paste0(tmp_out, ".raw")
  on.exit(unlink(c(tmp_snplist, paste0(tmp_out, c(".raw", ".log")))), add = TRUE)
  if (!file.exists(raw_file)) {
    msg("skip examples: PLINK2 genotype export failed")
    return(NULL)
  }
  geno_dt <- fread(raw_file)
  if (!"IID" %in% names(geno_dt) || ncol(geno_dt) < 7) {
    return(NULL)
  }
  value_cols <- names(geno_dt)[7:ncol(geno_dt)]
  long_dt <- melt(
    geno_dt,
    id.vars = "IID",
    measure.vars = value_cols,
    variable.name = "variant_id_raw",
    value.name = "genotype_numeric"
  )
  long_dt[, sample := as.character(IID)]
  long_dt[, counted_allele := sub("^.*_([^_]+)$", "\\1", as.character(variant_id_raw))]
  long_dt[, variant_id := sub("_[^_]+$", "", as.character(variant_id_raw))]
  long_dt[, genotype_numeric := as.numeric(genotype_numeric)]
  long_dt[, .(sample, variant_id, genotype_numeric, counted_allele)]
}

format_genotype <- function(variant_id, genotype_numeric, counted_allele) {
  parts <- strsplit(variant_id, ":", fixed = TRUE)[[1]]
  ref <- if (length(parts) >= 3) parts[3] else "Ref"
  alt <- if (length(parts) >= 4) parts[4] else "Alt"
  counted <- unique(na.omit(counted_allele))
  counted <- if (length(counted) > 0) counted[1] else alt
  uncounted <- if (identical(counted, ref)) alt else ref
  labels <- c(
    paste0(uncounted, "/", uncounted),
    paste0(uncounted, "/", counted),
    paste0(counted, "/", counted)
  )
  factor(
    cut(genotype_numeric, breaks = c(-Inf, 0.5, 1.5, Inf), labels = labels),
    levels = labels,
    ordered = TRUE
  )
}

format_genotype_dosage <- function(variant_id, genotype_numeric, counted_allele) {
  parts <- strsplit(variant_id, ":", fixed = TRUE)[[1]]
  ref <- if (length(parts) >= 3) parts[3] else "Ref"
  alt <- if (length(parts) >= 4) parts[4] else "Alt"
  counted <- unique(na.omit(counted_allele))
  counted <- if (length(counted) > 0) counted[1] else alt
  uncounted <- if (identical(counted, ref)) alt else ref
  dosage_labels <- c(
    paste0(uncounted, "/", uncounted),
    paste0(uncounted, "/", counted),
    paste0(counted, "/", counted)
  )
  ref_alt_levels <- c(
    paste0(ref, "/", ref),
    paste0(ref, "/", alt),
    paste0(alt, "/", alt)
  )
  dosage_class <- as.character(cut(genotype_numeric, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("G0", "G1", "G2")))
  if (identical(counted, ref)) {
    genotype <- dplyr::case_when(
      dosage_class == "G0" ~ paste0(alt, "/", alt),
      dosage_class == "G1" ~ paste0(ref, "/", alt),
      dosage_class == "G2" ~ paste0(ref, "/", ref),
      TRUE ~ NA_character_
    )
  } else if (identical(counted, alt)) {
    genotype <- dplyr::case_when(
      dosage_class == "G0" ~ paste0(ref, "/", ref),
      dosage_class == "G1" ~ paste0(ref, "/", alt),
      dosage_class == "G2" ~ paste0(alt, "/", alt),
      TRUE ~ NA_character_
    )
  } else {
    genotype <- as.character(cut(genotype_numeric, breaks = c(-Inf, 0.5, 1.5, Inf), labels = dosage_labels))
  }
  factor(genotype, levels = ref_alt_levels, ordered = TRUE)
}

build_example_candidate_pool <- function(inter_dt, max_candidates) {
  inter_dt <- as.data.table(inter_dt)
  per_pattern <- ceiling(max_candidates / max(1L, length(pattern_order)))
  candidate_dt <- rbindlist(lapply(pattern_order, function(pattern_name) {
    inter_dt[pattern == pattern_name][order(pval_gi, -abs_interaction)][seq_len(min(per_pattern, .N))]
  }), use.names = TRUE, fill = TRUE)
  if (nrow(candidate_dt) < max_candidates) {
    fill_dt <- inter_dt[!interaction_key %in% candidate_dt$interaction_key][order(pval_gi, -abs_interaction)]
    candidate_dt <- rbindlist(
      list(candidate_dt, fill_dt[seq_len(min(max_candidates - nrow(candidate_dt), nrow(fill_dt)))]),
      use.names = TRUE,
      fill = TRUE
    )
  }
  candidate_dt[seq_len(min(max_candidates, .N))]
}

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
D <- file.path(APA_ROOT, "data", "05_gxs")
fig_dir <- file.path(APA_ROOT, "figures", "05_gxs")
example_dir <- file.path(fig_dir, "examples")
interaction_file <- file.path(D, "sig.interaction.result")
phenotype_group_file <- "pre_input/phenotype.group"
prenatal_signif_file <- "../04_coloc/input/prenatal.signif_pairs.txt.gz"
postnatal_signif_file <- "../04_coloc/input/postnatal.signif_pairs.txt.gz"
phenotype_file <- file.path(D, "phenotype.bed.gz")
lifestage_file <- file.path(D, "LifeStage.txt")
genotype_pfile <- file.path(APA_ROOT, "raw", "genotype")

pattern_order <- c("Flip", "Amplifying", "Attenuating")
pattern_colors <- c(Flip = "#B44E45", Amplifying = "#4C78A8", Attenuating = "#6A8F6B")
status_colors <- c(
  `Interaction only` = "#4D4D4D",
  Both = "#BDBDBD",
  `Prenatal only` = "#D99000",
  `Postnatal only` = "#2F78B7"
)
stage_colors <- c(Prenatal = "#D99000", Postnatal = "#2F78B7")
base_colors <- c(ink = "#222222", muted = "#666666", grid = "#E6E6E6", line = "#A6A6A6")
main_base <- 7.6

dir.create(D, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_bw(base_size = 11, base_family = "") +
  theme(
    plot.title = element_text(face = "bold", size = 11, color = base_colors[["ink"]]),
    plot.subtitle = element_blank(),
    plot.caption = element_text(size = 11 * 0.68, color = base_colors[["muted"]], hjust = 0),
    panel.grid.major = element_line(color = base_colors[["grid"]], linewidth = 0.25),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = base_colors[["ink"]], linewidth = 0.35),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11 * 0.77),
    legend.text = element_text(size = 11 * 0.77),
    axis.title = element_text(face = "bold", color = base_colors[["ink"]]),
    axis.text = element_text(color = base_colors[["ink"]]),
    axis.ticks = element_line(color = base_colors[["ink"]], linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = base_colors[["ink"]])
  )

# --- interaction results ---
inter_dt <- fread(interaction_file)
inter_dt[, gene_symbol := vapply(
  strsplit(as.character(phenotype_id), "\\|"),
  function(x) if (length(x) >= 2) x[2] else x[1],
  character(1)
)]
if (file.exists(phenotype_group_file)) {
  group_dt <- fread(
    phenotype_group_file,
    header = FALSE,
    sep = "\t",
    col.names = c("phenotype_id", "group_gene")
  )
  inter_dt <- merge(inter_dt, group_dt, by = "phenotype_id", all.x = TRUE)
  inter_dt[, gene_symbol := ifelse(is.na(group_gene) | group_gene == "", gene_symbol, group_gene)]
  inter_dt[, group_gene := NULL]
}
inter_dt[, pattern := ifelse(
  b_g * (b_g + b_gi) < 0,
  "Flip",
  ifelse(abs(b_g + b_gi) > abs(b_g), "Amplifying", "Attenuating")
)]
inter_dt[, slope_pre := b_g]
inter_dt[, slope_post := b_g + b_gi]
inter_dt[, abs_interaction := abs(b_gi)]
inter_dt[, interaction_key := pair_key(phenotype_id, variant_id)]

# --- separate-mapping status (prenatal / postnatal significant pairs) ---
prenatal_signif <- read_signif_pairs(prenatal_signif_file)
postnatal_signif <- read_signif_pairs(postnatal_signif_file)
pre_keys <- if (!is.null(prenatal_signif)) pair_key(prenatal_signif$phenotype_id, prenatal_signif$variant_id) else character()
post_keys <- if (!is.null(postnatal_signif)) pair_key(postnatal_signif$phenotype_id, postnatal_signif$variant_id) else character()
annotated_dt <- copy(inter_dt)
annotated_dt[, prenatal_sig := interaction_key %in% pre_keys]
annotated_dt[, postnatal_sig := interaction_key %in% post_keys]
annotated_dt[, separate_status := fifelse(
  prenatal_sig & postnatal_sig, "Both",
  fifelse(prenatal_sig, "Prenatal only",
          fifelse(postnatal_sig, "Postnatal only", "Neither"))
)]
annotated_dt[, separate_status := factor(separate_status, levels = c("Both", "Prenatal only", "Postnatal only", "Neither"))]

# --- panel A: stage-specific effects ---
dt <- copy(as.data.table(inter_dt))
dt[, pattern := factor(pattern, levels = c("Flip", "Amplifying", "Attenuating"))]
label_dt <- dt[order(pval_gi, -abs_interaction)]
label_dt <- label_dt[!duplicated(gene_symbol)][seq_len(min(6L, .N))]
max_abs <- max(abs(c(dt$slope_pre, dt$slope_post)), na.rm = TRUE)
lim <- c(-max_abs, max_abs) * 1.06

p <- ggplot(dt, aes(x = slope_pre, y = slope_post, color = pattern)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = base_colors[["line"]], linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = base_colors[["line"]], linewidth = 0.35) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#5F5F5F", linewidth = 0.4) +
  geom_point(size = 1.95, alpha = 0.66, stroke = 0) +
  geom_text_repel(
    data = label_dt,
    aes(label = gene_symbol),
    size = 2.25,
    color = base_colors[["ink"]],
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.12,
    max.overlaps = Inf,
    seed = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(values = pattern_colors, drop = FALSE) +
  coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
  labs(
    title = "Stage-specific effects",
    x = expression(beta[G]~"prenatal"),
    y = expression(beta[G]~"postnatal"),
    color = "Pattern"
  ) +
  theme_main_minimal(base_size = main_base) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = main_base + 1.1, color = base_colors[["ink"]]),
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.12),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.82), color = NA),
    legend.key.width = unit(0.7, "lines"),
    plot.margin = margin(4, 6, 4, 4)
  )
panel_a <- add_panel_tag(p, "A")

# --- panel B: interaction-only discovery tiles ---
status_dt <- annotated_dt[, .N, by = separate_status]
status_dt[, status_label := fifelse(separate_status == "Neither", "Interaction only", as.character(separate_status))]
status_dt[, status_label := factor(status_label, levels = c("Interaction only", "Postnatal only", "Prenatal only", "Both"))]
setorder(status_dt, status_label)
tile_dt <- status_dt[, .(tile_id = seq_len(N)), by = status_label]
tile_dt[, global_id := seq_len(.N)]
tile_dt[, x := ((global_id - 1L) %% 40L) + 1L]
tile_dt[, y := ceiling(global_id / 40L)]
tile_dt[, status_label := factor(status_label, levels = c("Interaction only", "Postnatal only", "Prenatal only", "Both"))]
main_n <- status_dt[status_label == "Interaction only", N]
total_n <- sum(status_dt$N)
main_fraction <- main_n / total_n

p <- ggplot(tile_dt, aes(x = x, y = -y, fill = status_label)) +
  geom_tile(width = 0.86, height = 0.86, color = "white", linewidth = 0.08) +
  scale_fill_manual(values = status_colors, drop = FALSE, name = NULL) +
  coord_equal(clip = "off") +
  labs(
    title = "Interaction-only discovery",
    subtitle = sprintf("%s/%s (%s) interaction-only", comma(main_n), comma(total_n), percent(main_fraction, accuracy = 0.1)),
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = main_base, base_family = "") +
  theme(
    plot.title = element_text(face = "bold", size = main_base + 1, color = base_colors[["ink"]], hjust = 0.5),
    plot.subtitle = element_text(face = "bold", size = 7.6, color = base_colors[["ink"]], hjust = 0.5, margin = margin(t = 1, b = 3)),
    legend.position = "bottom",
    legend.text = element_text(size = 6.3),
    legend.key.size = unit(0.3, "lines"),
    plot.margin = margin(4, 3, 2, 3)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
panel_b <- add_panel_tag(p, "B")

# --- panel C: pattern composition and effect size ---
dt <- copy(as.data.table(inter_dt))
dt[, pattern := factor(pattern, levels = c("Flip", "Amplifying", "Attenuating"))]
pattern_dt <- dt[, .(
  N = .N,
  fraction = .N / nrow(dt),
  median_abs_beta_gi = median(abs_interaction, na.rm = TRUE),
  q25_abs_beta_gi = quantile(abs_interaction, 0.25, na.rm = TRUE),
  q75_abs_beta_gi = quantile(abs_interaction, 0.75, na.rm = TRUE)
), by = pattern][order(pattern)]
pattern_dt[, xmax := cumsum(fraction)]
pattern_dt[, xmin := shift(xmax, fill = 0)]
pattern_dt[, xmid := (xmin + xmax) / 2]
pattern_dt[, label := sprintf("%s\n%s", pattern, percent(fraction, accuracy = 0.1))]

composition <- ggplot(pattern_dt, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = pattern)) +
  geom_rect(color = "white", linewidth = 0.35) +
  geom_text(aes(x = xmid, y = 0.5, label = label), size = 2.05, fontface = "bold", color = "white", lineheight = 0.88) +
  scale_fill_manual(values = pattern_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Dynamic pattern composition", x = NULL, y = NULL) +
  theme_void(base_size = 8, base_family = "") +
  theme(
    plot.title = element_text(face = "bold", size = 9, color = base_colors[["ink"]], hjust = 0),
    plot.margin = margin(4, 3, 0, 3)
  )

effect <- ggplot(pattern_dt, aes(x = median_abs_beta_gi, y = pattern, color = pattern)) +
  geom_linerange(aes(xmin = q25_abs_beta_gi, xmax = q75_abs_beta_gi), linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.0) +
  geom_text(aes(label = sprintf("n=%s", N)), x = Inf, hjust = 1.05, size = 1.95, color = base_colors[["muted"]]) +
  scale_color_manual(values = pattern_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = expression(median~"|"*beta[GI]*"|"~"(IQR)"), y = NULL) +
  theme_main_minimal(base_size = 7.2) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(0, 3, 4, 3)
  )
panel_c <- add_panel_tag(composition / effect + plot_layout(heights = c(0.5, 1)), "C")

# --- panel D: representative loci (genotype x stage, main examples) ---
cand_main <- build_example_candidate_pool(inter_dt, 80L)
stage_dt <- read_lifestage(lifestage_file)
pheno_dt <- read_phenotype_subset(phenotype_file, cand_main$phenotype_id)
geno_dt <- load_genotype_subset(genotype_pfile, cand_main$variant_id)
main_dt <- NULL
if (is.null(stage_dt) || is.null(pheno_dt) || is.null(geno_dt)) {
  msg("skip main examples: unable to load stage, phenotype, or genotype subset")
} else {
  merged_dt <- merge(cand_main, pheno_dt, by = "phenotype_id", allow.cartesian = TRUE)
  merged_dt <- merge(merged_dt, geno_dt, by = c("sample", "variant_id"))
  merged_dt <- merge(merged_dt, stage_dt, by = "sample")
  merged_dt <- merged_dt[is.finite(phenotype) & is.finite(genotype_numeric)]
  if (nrow(merged_dt) > 0) {
    merged_dt[, genotype := format_genotype(variant_id[1], genotype_numeric, counted_allele), by = variant_id]
    merged_dt[, genotype_dosage := format_genotype_dosage(variant_id[1], genotype_numeric, counted_allele), by = variant_id]
    merged_dt[, genotype_class := factor(
      as.character(cut(genotype_numeric, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("G0", "G1", "G2"))),
      levels = c("G0", "G1", "G2"),
      ordered = TRUE
    )]
    merged_dt <- merged_dt[!is.na(genotype)]
    count_dt <- merged_dt[, .N, by = .(interaction_key, stage, genotype)]
    valid_keys <- count_dt[, .(
      n_cells = .N,
      min_n = min(N)
    ), by = interaction_key][n_cells == 6 & min_n >= 3L, interaction_key]
    selected_keys <- cand_main[0]
    if (length(valid_keys) > 0) {
      available <- copy(cand_main[interaction_key %in% valid_keys])
      available <- available[order(pval_gi, -abs_interaction)]
      available[, display_key := paste(pattern, gene_symbol, variant_id, sep = "|")]
      available <- available[!duplicated(display_key)]
      picked <- rbindlist(lapply(pattern_order, function(pattern_name) {
        one_pattern <- available[pattern == pattern_name][seq_len(min(2L, .N))]
        one_pattern[, pattern_rank := seq_len(.N)]
        one_pattern
      }), use.names = TRUE, fill = TRUE)
      picked[, pattern_order_i := match(pattern, pattern_order)]
      picked <- picked[order(pattern_rank, pattern_order_i)]
      selected_keys <- picked[!duplicated(interaction_key)]
    }
    if (nrow(selected_keys) > 0) {
      selected_keys[, example_order := seq_len(.N)]
      out_dt <- merged_dt[interaction_key %in% selected_keys$interaction_key]
      out_dt <- merge(out_dt, selected_keys[, .(interaction_key, example_order)], by = "interaction_key")
      out_dt[, variant_label := vapply(variant_id, format_variant_label, character(1))]
      out_dt[, example_label := sprintf("%s | %s\n%s", gene_symbol, pattern, variant_label)]
      out_dt[, genotype_label := sprintf("%s: %s", genotype_class, genotype_dosage)]
      main_dt <- out_dt
    }
  }
}

panel_d_plot <- NULL
panel_d_data <- NULL
if (!is.null(main_dt)) {
  summary_dt <- main_dt[, .(
    mean_phenotype = mean(phenotype, na.rm = TRUE),
    se_phenotype = sd(phenotype, na.rm = TRUE) / sqrt(.N),
    n = .N
  ), by = .(example_order, example_label, interaction_key, pattern, stage, genotype_class, genotype_label, genotype_dosage)]
  summary_dt[, example_label := factor(example_label, levels = unique(example_label[order(example_order)]))]

  p <- ggplot(summary_dt, aes(x = genotype_class, y = mean_phenotype, color = stage, group = stage)) +
    geom_errorbar(aes(ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype),
                  width = 0.1, linewidth = 0.4, position = position_dodge(width = 0.18)) +
    geom_line(linewidth = 0.65, position = position_dodge(width = 0.18)) +
    geom_point(shape = 21, fill = "white", stroke = 0.6, size = 2.0, position = position_dodge(width = 0.18)) +
    facet_wrap(~ example_label, scales = "free_y", ncol = 3L) +
    scale_color_manual(values = stage_colors, name = "Life stage") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) +
    labs(title = "Representative loci", x = "Genotype", y = "Normalized PDUI") +
    theme_main_minimal(base_size = main_base) +
    theme(
      legend.position = "top",
      strip.text = element_text(size = main_base - 0.4, face = "bold", color = base_colors[["ink"]]),
      axis.text.x = element_text(size = main_base - 0.9),
      plot.margin = margin(1, 4, 3, 4)
    )
  panel_d_plot <- add_panel_tag(p, "D")
  panel_d_data <- unique(main_dt[, .(example_order, phenotype_id, variant_id, gene_symbol, pattern, pval_gi, b_g, b_gi, slope_pre, slope_post)])[order(example_order)]
}

# --- combined main figure ---
save_plot(panel_a, "panel_a_stage_specific_effects", fig_dir, 4.9, 4.6)
save_plot(panel_b, "panel_b_interaction_only_discovery", fig_dir, 4.4, 2.0)
save_plot(panel_c, "panel_c_dynamic_pattern_summary", fig_dir, 4.4, 2.3)
if (!is.null(panel_d_plot)) save_plot(panel_d_plot, "panel_d_representative_loci", fig_dir, 8.2, 3.0)

if (!is.null(panel_d_plot)) {
  right_stack <- panel_b / panel_c + plot_layout(heights = c(0.92, 1.08))
  top_row <- panel_a | right_stack + plot_layout(widths = c(1.05, 0.95))
  combined <- top_row / panel_d_plot +
    plot_layout(heights = c(1.12, 0.88), guides = "collect") &
    theme(legend.position = "bottom")
  combined <- combined +
    plot_annotation(theme = theme(plot.tag = element_blank()))
  save_plot(combined, "main_fig", fig_dir, width = 8.6, height = 9.4)
}

# --- extra examples: one plot per locus and a summary facet ---
cand_ex <- build_example_candidate_pool(inter_dt, 120L)
stage_dt <- read_lifestage(lifestage_file)
pheno_dt <- read_phenotype_subset(phenotype_file, cand_ex$phenotype_id)
geno_dt <- load_genotype_subset(genotype_pfile, cand_ex$variant_id)
ex_dt <- NULL
if (is.null(stage_dt) || is.null(pheno_dt) || is.null(geno_dt)) {
  msg("skip examples: unable to load stage, phenotype, or genotype subset")
} else {
  merged_dt <- merge(cand_ex, pheno_dt, by = "phenotype_id", allow.cartesian = TRUE)
  merged_dt <- merge(merged_dt, geno_dt, by = c("sample", "variant_id"))
  merged_dt <- merge(merged_dt, stage_dt, by = "sample")
  merged_dt <- merged_dt[is.finite(phenotype) & is.finite(genotype_numeric)]
  if (nrow(merged_dt) > 0) {
    merged_dt[, genotype := format_genotype(variant_id[1], genotype_numeric, counted_allele), by = variant_id]
    merged_dt[, genotype_dosage := format_genotype_dosage(variant_id[1], genotype_numeric, counted_allele), by = variant_id]
    merged_dt[, genotype_refalt_code := {
      parts <- strsplit(variant_id[1], ":", fixed = TRUE)[[1]]
      ref <- if (length(parts) >= 3) parts[3] else "Ref"
      alt <- if (length(parts) >= 4) parts[4] else "Alt"
      counted <- unique(na.omit(counted_allele))
      counted <- if (length(counted) > 0) counted[1] else alt
      if (identical(counted, ref)) 2 - genotype_numeric else genotype_numeric
    }, by = variant_id]
    merged_dt[, genotype_class := factor(
      as.character(cut(genotype_numeric, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("G0", "G1", "G2"))),
      levels = c("G0", "G1", "G2"),
      ordered = TRUE
    )]
    merged_dt <- merged_dt[!is.na(genotype)]
    count_dt <- merged_dt[, .N, by = .(interaction_key, stage, genotype)]
    valid_keys <- count_dt[, .(
      n_cells = .N,
      min_n = min(N)
    ), by = interaction_key][n_cells == 6 & min_n >= 3L, interaction_key]
    if (length(valid_keys) == 0) {
      msg("skip examples: no candidate has all genotypes in both stages")
    } else {
      available <- copy(cand_ex[interaction_key %in% valid_keys])
      available[, pattern_rank := match(pattern, pattern_order)]
      available[is.na(pattern_rank), pattern_rank := length(pattern_order) + 1L]
      available <- available[order(pattern_rank, pval_gi, -abs_interaction)]
      available <- available[!duplicated(paste(gene_symbol, variant_id, pattern, sep = "||"))]
      picked <- rbindlist(lapply(pattern_order, function(pattern_name) {
        available[pattern == pattern_name][seq_len(min(2L, .N))]
      }), use.names = TRUE, fill = TRUE)
      picked <- picked[!duplicated(interaction_key)]
      if (nrow(picked) < 6L) {
        fill_dt <- available[!interaction_key %in% picked$interaction_key]
        picked <- rbindlist(
          list(picked, fill_dt[seq_len(min(6L - nrow(picked), nrow(fill_dt)))]),
          use.names = TRUE,
          fill = TRUE
        )
      }
      selected_keys <- picked[seq_len(min(6L, .N))]
      selected_keys[, example_order := seq_len(.N)]
      out_dt <- merged_dt[interaction_key %in% selected_keys$interaction_key]
      out_dt <- merge(out_dt, selected_keys[, .(interaction_key, example_order)], by = "interaction_key")
      out_dt[, counted_allele_plot := unique(na.omit(counted_allele))[1], by = interaction_key]
      out_dt[, variant_label := vapply(variant_id, format_variant_label, character(1))]
      out_dt[, example_label := sprintf("%s | %s\n%s; counted %s", gene_symbol, pattern, variant_label, counted_allele_plot)]
      ex_dt <- out_dt
    }
  }
}

if (is.null(ex_dt) || nrow(ex_dt) == 0) {
  msg("skip extra examples: no valid genotype-by-stage examples")
} else {
  dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
  index_dt <- unique(ex_dt[, .(
    example_order, phenotype_id, variant_id, gene_symbol, pattern, pval_gi, b_g, b_gi, slope_pre, slope_post
  )])
  fwrite(index_dt, file.path(D, "example_plot_index.tsv"), sep = "\t")

  limit <- min(12L, nrow(index_dt))
  for (i in seq_len(limit)) {
    key <- index_dt$example_order[i]
    one_dt <- ex_dt[example_order == key]
    stem <- sprintf("%02d_%s_%s_%s", i, safe_filename(index_dt$pattern[i]), safe_filename(index_dt$gene_symbol[i]), safe_filename(index_dt$variant_id[i]))

    summary_dt <- one_dt[, .(
      mean_phenotype = mean(phenotype, na.rm = TRUE),
      se_phenotype = sd(phenotype, na.rm = TRUE) / sqrt(.N),
      n = .N,
      genotype_numeric = median(genotype_numeric, na.rm = TRUE),
      genotype_refalt_code = median(genotype_refalt_code, na.rm = TRUE)
    ), by = .(gene_symbol, variant_id, pattern, pval_gi, stage, genotype, genotype_dosage)]
    title_dt <- unique(one_dt[, .(gene_symbol, variant_id, pattern, pval_gi)])
    parts <- strsplit(title_dt$variant_id[1], ":", fixed = TRUE)[[1]]
    ref <- if (length(parts) >= 3) parts[3] else "Ref"
    alt <- if (length(parts) >= 4) parts[4] else "Alt"
    axis_labels <- c(paste0(ref, "/", ref), paste0(ref, "/", alt), paste0(alt, "/", alt))
    title_text <- sprintf(
      "%s | %s\n%s; Pgi=%s",
      title_dt$gene_symbol[1],
      title_dt$pattern[1],
      format_variant_label(title_dt$variant_id[1]),
      formatC(title_dt$pval_gi[1], format = "e", digits = 2)
    )

    single_plot <- ggplot() +
      geom_jitter(
        data = one_dt,
        aes(x = genotype_refalt_code, y = phenotype),
        color = base_colors[["muted"]],
        alpha = 0.16,
        size = 0.85,
        width = 0.07,
        height = 0,
        show.legend = FALSE
      ) +
      geom_errorbar(
        data = summary_dt,
        aes(x = genotype_refalt_code, ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype, color = stage, group = stage),
        width = 0.13,
        linewidth = 0.65,
        position = position_dodge(width = 0.35)
      ) +
      geom_line(
        data = summary_dt,
        aes(x = genotype_refalt_code, y = mean_phenotype, color = stage, group = stage),
        linewidth = 0.95,
        position = position_dodge(width = 0.35)
      ) +
      geom_point(
        data = summary_dt,
        aes(x = genotype_refalt_code, y = mean_phenotype, color = stage, size = n),
        shape = 21,
        fill = "white",
        stroke = 0.9,
        position = position_dodge(width = 0.35)
      ) +
      scale_color_manual(values = stage_colors, name = "Life stage") +
      scale_size_continuous(range = c(2.6, 5.2), name = "n") +
      scale_x_continuous(breaks = c(0, 1, 2), labels = axis_labels, limits = c(-0.45, 2.45), expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
      labs(
        title = title_text,
        x = "Genotype",
        y = "Normalized PDUI"
      ) +
      theme_classic(base_size = 11, base_family = "") +
      guides(size = "none") +
      theme(
        axis.line = element_line(color = base_colors[["ink"]], linewidth = 0.45),
        axis.ticks = element_line(color = base_colors[["ink"]], linewidth = 0.35),
        axis.ticks.length = unit(2.2, "pt"),
        axis.title = element_text(face = "bold", color = base_colors[["ink"]]),
        axis.text = element_text(color = base_colors[["ink"]]),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = base_colors[["grid"]], linewidth = 0.22),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 11, face = "bold", color = base_colors[["ink"]]),
        plot.margin = margin(8, 10, 8, 8)
      )
    save_plot(single_plot, stem, example_dir, width = 6, height = 5)
  }

  summary_ncol <- 3L
  summary_nrow <- ceiling(nrow(index_dt) / summary_ncol)
  summary_width <- 3.8 * summary_ncol
  summary_height <- 3.2 * summary_nrow + 1.1
  summary_dt <- ex_dt[, .(
    mean_phenotype = mean(phenotype, na.rm = TRUE),
    se_phenotype = sd(phenotype, na.rm = TRUE) / sqrt(.N),
    n = .N,
    genotype_numeric = median(genotype_numeric, na.rm = TRUE),
    genotype_refalt_code = median(genotype_refalt_code, na.rm = TRUE)
  ), by = .(example_order, gene_symbol, variant_id, pattern, stage, genotype_dosage)]
  summary_dt[, label := sprintf("%s | %s\n%s", gene_symbol, pattern, vapply(variant_id, format_variant_label, character(1)))]
  summary_dt[, label := factor(label, levels = unique(label[order(example_order)]))]
  summary_plot <- ggplot(summary_dt, aes(x = genotype_refalt_code, y = mean_phenotype, color = stage, group = stage)) +
    geom_errorbar(aes(ymin = mean_phenotype - se_phenotype, ymax = mean_phenotype + se_phenotype),
                  width = 0.15, linewidth = 0.45, position = position_dodge(width = 0.35)) +
    geom_line(linewidth = 0.75, position = position_dodge(width = 0.35)) +
    geom_point(aes(size = n), shape = 21, fill = "white", stroke = 0.65, position = position_dodge(width = 0.35)) +
    facet_wrap(~ label, scales = "free_y", ncol = summary_ncol) +
    scale_color_manual(values = stage_colors, name = "Life stage") +
    scale_size_continuous(range = c(1.8, 4.2), name = "n") +
    scale_x_continuous(breaks = c(0, 1, 2), labels = c("Ref/Ref", "Ref/Alt", "Alt/Alt"), limits = c(-0.45, 2.45), expand = c(0, 0)) +
    labs(
      title = "Representative Genotype-by-Stage Examples",
      x = "Genotype",
      y = "Normalized PDUI"
    ) +
    plot_theme +
    guides(size = "none") +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
      strip.text = element_text(size = 8.3, face = "bold", color = base_colors[["ink"]]),
      plot.title = element_text(size = 12, face = "bold", color = base_colors[["ink"]])
    )
  save_plot(summary_plot, "example_summary_facets", example_dir, width = summary_width, height = summary_height)
}

# --- summary tables ---
overview_dt <- data.table(
  metric = c(
    "significant_interaction_pairs",
    "unique_genes",
    "representative_examples"
  ),
  value = c(
    nrow(inter_dt),
    uniqueN(inter_dt$gene_symbol),
    ifelse(is.null(panel_d_data), NA_integer_, nrow(panel_d_data))
  )
)
fwrite(overview_dt, file.path(D, "figure_overview.tsv"), sep = "\t")

pattern_dt <- inter_dt[, .N, by = pattern][order(-N)]
pattern_dt[, fraction := N / sum(N)]
pattern_dt[, median_abs_beta_gi := inter_dt[.SD, on = "pattern", median(abs_interaction, na.rm = TRUE), by = .EACHI]$V1]
fwrite(pattern_dt, file.path(D, "pattern_counts.tsv"), sep = "\t")

status_dt <- annotated_dt[, .N, by = separate_status][order(separate_status)]
status_dt[, fraction := N / sum(N)]
fwrite(status_dt, file.path(D, "separate_status_counts.tsv"), sep = "\t")

recovery_dt <- annotated_dt[, .N, by = .(pattern, recovery = ifelse(postnatal_sig, "Recovered in postnatal", "Missed without prenatal"))]
recovery_dt[, fraction := N / sum(N)]
fwrite(recovery_dt, file.path(D, "prenatal_contribution_counts.tsv"), sep = "\t")

if (!is.null(panel_d_data)) {
  fwrite(panel_d_data, file.path(D, "panel_2a_representative_examples.tsv"), sep = "\t")
}
