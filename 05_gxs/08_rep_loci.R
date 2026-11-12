# Representative-loci panel: six G x Stage loci (2 per pattern), genotype-stratified PDUI
# means by life stage. Reads data/05_gxs and raw/genotype; writes figures/05_gxs.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
D <- file.path(APA_ROOT, "data", "05_gxs")
interaction_file <- file.path(D, "intx_hc3_core_groups.tsv")
phenotype_file <- file.path(D, "phenotype.bed.gz")
lifestage_file <- file.path(D, "LifeStage.txt")
genotype_pfile <- file.path(APA_ROOT, "raw", "genotype")
out_dir <- file.path(APA_ROOT, "figures", "05_gxs")
figure_stem <- "panel_representative_loci"

pval_g_cutoff <- 0.05
example_candidates_per_pattern <- 27L
examples_per_pattern <- 2L
min_n_per_cell <- 20L
figure_width <- 8.8
figure_height <- 8.4
dpi <- 420

ink <- "#1E252B"
muted_ink <- "#56616A"
grid_color <- "#DCE3E8"
band_fill <- "#F7F9FA"
stage_cols <- c(Prenatal = "#E69F00", Postnatal = "#0072B2")
base_size <- 11
pt_size <- 2.7
line_width <- 0.95
err_width <- 0.16
dodge_width <- 0.16

# --- interaction leads ---
inter_dt <- fread(interaction_file)
pattern_order <- c("Opposing", "Amplified", "Attenuated")
inter_dt[, gene_symbol := gene_id]
inter_dt[is.na(gene_symbol) | gene_symbol == "",
         gene_symbol := vapply(strsplit(lead_phenotype_id, "|", fixed = TRUE), function(x) x[2], character(1))]
inter_dt[, interaction_key := paste(lead_phenotype_id, lead_variant_id, sep = "||")]
inter_dt[, pattern := factor(point_estimate_class, levels = pattern_order)]
inter_dt[, abs_interaction := abs(beta_gxstage)]
inter_dt <- inter_dt[is.finite(hc3_pval_gxstage) & hc3_pval_gxstage < pval_g_cutoff]

pool_dt <- rbindlist(lapply(pattern_order, function(pat) {
  dt <- inter_dt[pattern == pat]
  dt[order(hc3_pval_gxstage, -abs_interaction)][seq_len(min(example_candidates_per_pattern, .N))]
}))

# --- phenotype values for the candidate pool (awk filter on the BED) ---
tmp_ids <- tempfile()
writeLines(unique(pool_dt$lead_phenotype_id), tmp_ids)
bed <- fread(cmd = sprintf(
  "gzip -dc %s | awk 'BEGIN{while((getline<\"%s\")>0) keep[$1]=1} NR==1 || ($4 in keep)'",
  shQuote(phenotype_file), tmp_ids))
unlink(tmp_ids)
pheno_dt <- melt(bed, id.vars = "pid", measure.vars = setdiff(names(bed), names(bed)[1:4]),
                 variable.name = "sample", value.name = "phenotype")
pheno_dt[, phenotype := as.numeric(phenotype)]
setnames(pheno_dt, "pid", "phenotype_id")
pheno_dt <- pheno_dt[, .(sample = as.character(sample), phenotype_id, phenotype)]

merged <- merge(pool_dt, pheno_dt, by.x = "lead_phenotype_id", by.y = "phenotype_id", allow.cartesian = TRUE)

# --- genotype dosages via plink2 --export A ---
tmp_snp <- tempfile()
tmp_out <- tempfile()
writeLines(unique(pool_dt$lead_variant_id), tmp_snp)
system2("plink2", c("--pfile", genotype_pfile, "--extract", tmp_snp, "--export", "A", "--out", tmp_out),
        stdout = FALSE, stderr = FALSE)
raw <- fread(paste0(tmp_out, ".raw"))
unlink(c(tmp_snp, paste0(tmp_out, c(".raw", ".log", ".nosex"))))
geno_dt <- melt(raw, id.vars = "IID", measure.vars = names(raw)[7:ncol(raw)],
                variable.name = "raw_id", value.name = "dosage")
geno_dt[, `:=`(
  sample = as.character(IID),
  counted_allele = sub("^.*_([^_]+)$", "\\1", as.character(raw_id)),
  variant_id = sub("_[^_]+$", "", as.character(raw_id)),
  dosage = as.numeric(dosage)
)]
geno_dt <- geno_dt[, .(sample, variant_id, dosage, counted_allele)]

merged <- merge(merged, geno_dt, by.x = c("sample", "lead_variant_id"), by.y = c("sample", "variant_id"))
stage_dt <- fread(lifestage_file)[, .(sample = as.character(sample),
  stage = factor(ifelse(LifeStage == 0, "Prenatal", "Postnatal"), levels = c("Prenatal", "Postnatal")))]
merged <- merge(merged, stage_dt, by = "sample")
merged <- merged[is.finite(phenotype) & is.finite(dosage)]

# --- genotype codes and labels (dosage counted on the REF/ALT allele of the variant ID) ---
merged[, g_code := {
  parts <- strsplit(lead_variant_id[1], ":", fixed = TRUE)[[1]]
  ref <- if (length(parts) >= 3) parts[3] else "Ref"
  counted <- unique(na.omit(counted_allele))[1]
  if (identical(counted, ref)) 2 - dosage else dosage
}, by = lead_variant_id]
merged[, g_label := {
  parts <- strsplit(lead_variant_id[1], ":", fixed = TRUE)[[1]]
  ref <- if (length(parts) >= 3) parts[3] else "Ref"
  alt <- if (length(parts) >= 4) parts[4] else "Alt"
  lv <- c(paste0(ref, "/", ref), paste0(ref, "/", alt), paste0(alt, "/", alt))
  factor(as.character(cut(g_code, breaks = c(-Inf, 0.5, 1.5, Inf), labels = lv)),
         levels = lv, ordered = TRUE)
}, by = lead_variant_id]
merged <- merged[!is.na(g_label)]

# --- six cells per locus (3 genotypes x 2 stages), each n >= min_n_per_cell ---
valid <- merged[, .N, by = .(interaction_key, stage, g_code)][
  , .(n_cells = .N, min_n = min(N)), by = interaction_key][
  n_cells == 6 & min_n >= min_n_per_cell, interaction_key]
available <- pool_dt[interaction_key %in% valid][order(hc3_pval_gxstage, -abs_interaction)][
  !duplicated(paste(pattern, gene_symbol, lead_variant_id, sep = "|"))]
selected <- rbindlist(lapply(pattern_order, function(pat) {
  available[pattern == pat][seq_len(min(examples_per_pattern, .N))]
}))[, example_order := seq_len(.N)]

example_dt <- merge(merged[interaction_key %in% selected$interaction_key],
                    selected[, .(interaction_key, example_order, pattern)], by = c("interaction_key", "pattern"))
example_dt[, strip_label := factor(sprintf("%s | %s", gene_symbol, pattern),
                                   levels = unique(sprintf("%s | %s", gene_symbol, pattern)[order(example_order)]))]

summary_d <- example_dt[, .(mean_p = mean(phenotype), se_p = sd(phenotype) / sqrt(.N), n = .N),
                        by = .(example_order, strip_label, interaction_key, stage, g_code, g_label)]
summary_d[, strip_label := factor(strip_label, levels = unique(strip_label[order(example_order)]))]

pval_d <- unique(example_dt[, .(example_order, strip_label, hc3_pval_gxstage)])
pval_d[, strip_label := factor(strip_label, levels = levels(summary_d$strip_label))]
pval_d[, pval_label := sprintf("P[GxS] = %s", ifelse(is.na(hc3_pval_gxstage), "NA",
  ifelse(hc3_pval_gxstage < 1e-4, formatC(hc3_pval_gxstage, format = "e", digits = 2),
         sprintf("%.3g", hc3_pval_gxstage))))]

n_d <- dcast(summary_d[, .(n = sum(n)), by = .(example_order, strip_label, g_code, stage)],
             example_order + strip_label + g_code ~ stage, value.var = "n")
n_d[, strip_label := factor(strip_label, levels = levels(summary_d$strip_label))]
n_d[, n_label := sprintf("Pre n=%s\nPost n=%s", Prenatal, Postnatal)]

axis_labels_dt <- unique(summary_d[, .(example_order, strip_label, g_code, g_label)])
axis_labels_dt[, strip_label := factor(strip_label, levels = levels(summary_d$strip_label))]

y_min <- min(summary_d$mean_p - summary_d$se_p, na.rm = TRUE)
y_max <- max(summary_d$mean_p + summary_d$se_p, na.rm = TRUE)
y_rng <- y_max - y_min
n_d[, label_y := y_min - 0.095 * y_rng]
axis_labels_dt[, label_y := y_min - 0.235 * y_rng]
band_d <- unique(summary_d[, .(example_order, strip_label)])
band_d[, `:=`(band_ymin = y_min - 0.285 * y_rng, band_ymax = y_min - 0.015 * y_rng)]

# --- plot ---
p <- ggplot(summary_d, aes(x = g_code, y = mean_p, color = stage, group = stage)) +
  geom_rect(data = band_d, aes(xmin = -Inf, xmax = Inf, ymin = band_ymin, ymax = band_ymax),
            inherit.aes = FALSE, fill = band_fill, color = NA) +
  geom_hline(yintercept = 0, color = grid_color, linewidth = 0.35, linetype = "dashed") +
  geom_errorbar(aes(ymin = mean_p - se_p, ymax = mean_p + se_p), width = err_width, linewidth = 0.55,
                position = position_dodge(dodge_width)) +
  geom_line(aes(linetype = stage), linewidth = line_width, position = position_dodge(dodge_width)) +
  geom_point(aes(shape = stage), size = pt_size, stroke = 0.9, fill = "white",
             position = position_dodge(dodge_width)) +
  geom_text(data = pval_d, aes(x = Inf, y = Inf, label = pval_label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.35, size = 3.15, fontface = "bold", color = ink) +
  geom_text(data = n_d, aes(x = g_code, y = label_y, label = n_label), inherit.aes = FALSE,
            color = muted_ink, size = 2.85, lineheight = 0.82, vjust = 0.5) +
  geom_text(data = axis_labels_dt, aes(x = g_code, y = label_y, label = g_label), inherit.aes = FALSE,
            color = "black", fontface = "bold", size = 3.25, angle = 12, hjust = 1, vjust = 0.5) +
  facet_wrap(~strip_label, nrow = 3, ncol = 2) +
  scale_color_manual(values = stage_cols, name = "Life stage") +
  scale_shape_manual(values = c(Prenatal = 21, Postnatal = 24), name = "Life stage") +
  scale_linetype_manual(values = c(Prenatal = "solid", Postnatal = "dashed"), name = "Life stage") +
  scale_x_continuous(breaks = 0:2, labels = NULL, limits = c(-0.35, 2.35), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y_min - 0.315 * y_rng, y_max + 0.16 * y_rng), expand = c(0, 0)) +
  labs(title = "Representative loci", x = "Genotype", y = "Normalized PDUI") +
  theme_classic(base_size = base_size) +
  theme(
    text = element_text(color = ink),
    plot.title = element_text(face = "bold", size = base_size + 5, hjust = 0.5, margin = margin(b = 8)),
    axis.title = element_text(face = "bold", size = base_size + 1),
    axis.text = element_text(size = base_size - 1, color = "black"),
    axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.y = element_line(color = ink, linewidth = 0.45),
    axis.ticks.y = element_line(color = ink, linewidth = 0.45),
    axis.title.x = element_text(margin = margin(t = 18)),
    axis.title.y = element_text(margin = margin(r = 6)),
    strip.text = element_text(face = "bold", size = base_size + 1, color = ink, margin = margin(t = 4, b = 4)),
    strip.background = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = base_size + 1),
    legend.text = element_text(size = base_size),
    legend.key.size = unit(0.9, "lines"), legend.margin = margin(t = 4),
    plot.margin = margin(8, 8, 8, 8)
  )

# --- save and layout check ---
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(out_dir, paste0(figure_stem, ".png"))
pdf_file <- file.path(out_dir, paste0(figure_stem, ".pdf"))
svg_file <- file.path(out_dir, paste0(figure_stem, ".svg"))
ggsave(png_file, p, width = figure_width, height = figure_height, dpi = dpi, bg = "white")
ggsave(pdf_file, p, width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(svg_file, p, width = figure_width, height = figure_height, device = svglite::svglite, bg = "white")
}

panel_order <- summary_d[, .(n_genotypes = uniqueN(g_code), min_n = min(n), max_n = max(n)), by = strip_label]
stage_n <- summary_d[, .(n_stage_labels = uniqueN(stage)), by = .(strip_label, g_code)]
errorbar_summary <- summary_d[, .(min_se = min(se_p, na.rm = TRUE), median_se = median(se_p, na.rm = TRUE),
                                  max_se = max(se_p, na.rm = TRUE), zero_or_na_se = sum(!is.finite(se_p) | se_p <= 0))]
selected_summary <- selected[, paste(sprintf("%s:%s", gene_symbol, pattern), collapse = "|")]
layout_checks_pass <- nrow(panel_order) == examples_per_pattern * length(pattern_order) &&
  all(panel_order$n_genotypes == 3L) && all(stage_n$n_stage_labels == 2L) &&
  nrow(n_d) == examples_per_pattern * length(pattern_order) * 3L && errorbar_summary$zero_or_na_se == 0L

layout_validation <- rbindlist(list(
  data.table(check = "source", metric = "interaction_file", value = interaction_file),
  data.table(check = "source", metric = "pvalue_label_source", value = "hc3_pval_gxstage"),
  data.table(check = "selection", metric = "selected_loci", value = selected_summary),
  data.table(check = "selection", metric = "examples_per_pattern", value = as.character(examples_per_pattern)),
  data.table(check = "selection", metric = "min_genotype_stage_n_threshold", value = as.character(min_n_per_cell)),
  data.table(check = "panel", metric = "facet_cells", value = as.character(nrow(panel_order))),
  data.table(check = "panel", metric = "all_facets_have_three_ordered_allele_labels", value = as.character(all(panel_order$n_genotypes == 3L))),
  data.table(check = "panel", metric = "all_genotype_cells_have_two_stage_labels", value = as.character(all(stage_n$n_stage_labels == 2L))),
  data.table(check = "panel", metric = "genotype_stage_n_labels", value = as.character(nrow(n_d))),
  data.table(check = "panel", metric = "cell_n_min_max", value = sprintf("%s/%s", min(panel_order$min_n), max(panel_order$max_n))),
  data.table(check = "layout", metric = "facet_grid", value = "2 columns x 3 rows"),
  data.table(check = "statistics", metric = "errorbar_se_min_median_max", value = sprintf("%.4g/%.4g/%.4g", errorbar_summary$min_se, errorbar_summary$median_se, errorbar_summary$max_se)),
  data.table(check = "statistics", metric = "zero_or_na_se", value = as.character(errorbar_summary$zero_or_na_se)),
  data.table(check = "layout", metric = "count_label_band", value = "inside each facet bottom band"),
  data.table(check = "layout", metric = "minimum_visible_text_pt", value = sprintf("%.1f", min(2.85, 3.15, 3.25) * 2.845)),
  data.table(check = "layout", metric = "layout_checks_pass", value = as.character(layout_checks_pass)),
  data.table(check = "export", metric = "dimensions_inches", value = sprintf("%.1f x %.1f", figure_width, figure_height)),
  data.table(check = "export", metric = "files", value = paste(basename(c(png_file, pdf_file, svg_file[file.exists(svg_file)])), collapse = "|"))
))
fwrite(layout_validation, file.path(out_dir, paste0(figure_stem, "_layout_validation.tsv")), sep = "\t")

stopifnot(
  nrow(inter_dt) == 193L,
  nrow(selected) == examples_per_pattern * length(pattern_order),
  layout_checks_pass,
  all(file.exists(c(png_file, pdf_file))),
  file.exists(svg_file) || !requireNamespace("svglite", quietly = TRUE)
)
