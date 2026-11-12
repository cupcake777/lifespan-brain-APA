#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
  library(jsonlite)
})

source("src/utils_helpers.R")
source("src/utils_plot.R")
source("src/utils_plot_locuscompare.R")
source("src/utils_output.R")
source("examples/synthetic_plot_data_v2.R")

cfg_global <- yaml::read_yaml("config/global.yml")
plot_inputs <- resolve_synthetic_plot_inputs(cfg_global)

make_dense_corner_points <- function() {
  data.frame(
    x = c(seq(0.02, 0.24, length.out = 28), seq(0.68, 0.86, length.out = 4)),
    y = c(seq(0.66, 0.98, length.out = 28), seq(0.12, 0.22, length.out = 4)),
    is_lead = c(TRUE, rep(FALSE, 31)),
    is_credible = c(rep(TRUE, 4), rep(FALSE, 28)),
    stringsAsFactors = FALSE
  )
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

outline_probe <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_point(shape = 21, fill = "red", color = "black", stroke = 1)
outline_probe <- locuscompare_remove_point_outlines(outline_probe)
outline_layer <- outline_probe$layers[[1]]
assert_true(identical(outline_layer$aes_params$colour, "transparent"), "point-outline remover should make point outlines transparent")
assert_true(identical(outline_layer$aes_params$stroke, 0), "point-outline remover should set stroke to zero")
lowercase_r2_probe <- locuscompare_ld_pairs(data.frame(SNP_A = "lead", SNP_B = "proxy", r2 = 0.5))
correlation_probe <- locuscompare_ld_pairs(data.frame(SNP_A = "lead", SNP_B = "proxy", R = 0.5))
assert_true(isTRUE(all.equal(lowercase_r2_probe$R2, 0.5)), "lowercase r2 input must not be squared again")
assert_true(isTRUE(all.equal(correlation_probe$R2, 0.25)), "correlation R input must be converted to R2")
plot_pair_probe <- data.frame(SNP_A = "lead", SNP_B = "proxy", R = 0.5)
plot_pair_normalized <- normalize_plot_ld_pairs(plot_pair_probe)
assert_true(
  nrow(plot_pair_normalized) == 1L && identical(plot_pair_normalized$SNP_A, "lead"),
  "LD_plot pair-table normalization must preserve lead-oriented R pairs"
)
label_probe <- ggplot(data.frame(x = 1, y = 1, label = "lead"), aes(x, y)) +
  locuscompare_protected_label_layer(
    data = data.frame(x = 1, y = 1, label = "lead"),
    mapping = aes(x = x, y = y, label = label)
  )
assert_true(inherits(label_probe$layers[[1]]$geom, "GeomLabelRepel"), "protected labels should use a repel label layer")

stress_genes <- data.frame(
  gene_key = paste0("stress_", seq_len(9)),
  symbol = paste0("DENSEGENE", seq_len(9)),
  start = c(100000, 150000, 200000, 250000, 300000, 350000, 400000, 205000, 255000),
  end = c(105000, 155000, 205000, 255000, 305000, 355000, 405000, 245000, 295000),
  strand = rep(c("+", "-"), length.out = 9),
  stringsAsFactors = FALSE
)
stress_lanes <- locuscompare_assign_gene_lanes(stress_genes, c(0, 1e6))
assert_true(isTRUE(attr(stress_lanes, "interval_overlap_free")), "dense gene structures should not overlap within lanes")
assert_true(max(stress_lanes$lane) >= 2L, "overlapping stress genes should be assigned to separate lanes")
stress_labels <- locuscompare_select_gene_labels(
  stress_lanes,
  position_range = c(0, 1e6),
  lead_position = 250000,
  panel_width_in = 1.2
)
assert_true(isTRUE(attr(stress_labels, "label_overlap_free")), "visible dense gene labels should not overlap within lanes")
assert_true(any(!stress_labels$label_visible), "dense label stress test should exercise audited label suppression")

vertical_stress <- data.frame(
  gene_key = paste0("vertical_", seq_len(20)),
  symbol = paste0("STACKED", seq_len(20)),
  start = rep(220000, 20),
  end = rep(280000, 20),
  strand = rep(c("+", "-"), 10),
  stringsAsFactors = FALSE
)
vertical_lanes <- locuscompare_assign_gene_lanes(vertical_stress, c(0, 1e6))
vertical_labels <- locuscompare_select_gene_labels(
  vertical_lanes,
  position_range = c(0, 1e6),
  lead_position = 250000,
  panel_width_in = 3.4
)
assert_true(isTRUE(attr(vertical_labels, "label_overlap_free")), "stacked gene labels should pass two-dimensional collision checks")
assert_true(any(!vertical_labels$label_visible), "deeply stacked genes should suppress labels that exceed vertical print capacity")

assert_true(!is.null(plot_inputs$gtf_path), "smoke plotting test requires a resolvable gene annotation fixture")
assert_true(!is.null(plot_inputs$recomb_path), "smoke plotting test requires a resolvable CHB recombination fixture")

example_data <- make_synthetic_locus()
assert_true(nrow(example_data$merged_data) >= 400L, "official-style demo requires a dense synthetic locus")
dense_corner_points <- make_dense_corner_points()

feature_region <- parse_feature_region(example_data$phenotype_info)
assert_true(!is.null(feature_region), "parse_feature_region() returned NULL")
assert_true(feature_region$chrom == "chr1", "feature_region chrom parse failed")
assert_true(feature_region$start == 42482663, "feature_region start parse failed")
assert_true(feature_region$end == 42484158, "feature_region end parse failed")

window_info <- resolve_plot_window(
  leadSNP_DF = example_data$merged_data,
  pos_col = "POS.qtl",
  chr_col = "CHR.qtl",
  lead_snp_val = example_data$lead_snp,
  plot_window_bp = 200000,
  phenotype_info = example_data$phenotype_info
)
assert_true(!is.null(window_info), "resolve_plot_window() returned NULL")
assert_true(window_info$chrom_label == "chr1", "resolve_plot_window() chrom label failed")
assert_true(isTRUE(window_info$expanded), "resolve_plot_window() should expand to cover phenotype transcript")
assert_true(window_info$min_pos == 42482663, "resolve_plot_window() expanded min_pos is incorrect")
assert_true(window_info$max_pos == 43766000, "resolve_plot_window() expanded max_pos is incorrect")

legend_layout <- choose_internal_legend_position(
  points_df = dense_corner_points,
  xlim = c(0, 1),
  ylim = c(0, 1),
  legend_nrow = 6,
  legend_labels = c("Lead SNP", "95% CS", "0.8-1.0", "0.6-0.8", "0.4-0.6", "< 0.2"),
  plot_width = 6.6,
  plot_height = 4.9
)
assert_true(
  !(legend_layout$position[1] < 0.5 && legend_layout$position[2] > 0.5),
  "dynamic legend selector should avoid dense top-left corner"
)

external_legend_right <- choose_external_legend_layout(
  legend_nrow = 6,
  plot_width = 6.6,
  plot_height = 4.9
)
assert_true(
  identical(external_legend_right$position, "right"),
  "standard locus plot should place legend on the right outside the panel"
)

external_legend_bottom <- choose_external_legend_layout(
  legend_nrow = 7,
  plot_width = 5.2,
  plot_height = 4.4
)
assert_true(
  identical(external_legend_bottom$position, "bottom"),
  "narrow plot should place legend at the bottom outside the panel"
)

assign("lead_SNP", example_data$lead_snp, envir = .GlobalEnv)
assign("geneSymbol", example_data$phenotype_info, envir = .GlobalEnv)
assign("plink_bfile", NULL, envir = .GlobalEnv)
assign("trait", "SMOKE_GWAS", envir = .GlobalEnv)

plot_obj <- plot_qtl_association(
  qtl_all_chrom = "CHR.qtl",
  qtl_all_pvalue = "P.qtl",
  leadSNP_DF = example_data$merged_data,
  ld_df = NULL,
  gtf_path = plot_inputs$gtf_path,
  region_recomb = NULL,
  recomb_path = plot_inputs$recomb_path,
  show_lead_line = FALSE,
  qtl_type = example_data$qtl_type,
  phenotype_info = example_data$phenotype_info,
  significance_threshold = 5e-8,
  significance_label = NULL,
  title_phenotype_field = "gene",
  plot_width = 6.6,
  plot_height = 4.9,
  credible_set = example_data$credible_set,
  plot_window_bp = 200000
)

tmp_pdf <- tempfile(fileext = ".pdf")
tmp_png <- tempfile(fileext = ".png")
save_res <- save_plot_with_fallback(
  plot_obj = plot_obj,
  pdf_path = tmp_pdf,
  png_path = tmp_png,
  plot_width = 6.6,
  plot_height = 4.9,
  png_scale = 1,
  png_dpi = 200
)

output_path <- save_res$path
assert_true(file.exists(output_path), "smoke plot output file was not created")
assert_true(file.info(output_path)$size > 5000, "smoke plot output file is suspiciously small")

locuscompare_obj <- make_combined_coloc_plot(
  df = example_data$merged_data,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set,
  title = "Colocalization evidence",
  pp4 = 0.92,
  pp3 = 0.04,
  phenotype_info = example_data$phenotype_info,
  gtf_path = plot_inputs$gtf_path,
  backend = "auto",
  title1 = "GWAS",
  title2 = "QTL"
)
assert_true(
  identical(attr(locuscompare_obj, "easycoloc_plot_backend"), "locuscomparer"),
  "auto backend should use locuscomparer when installed and conversion succeeds"
)
direct_label_panels <- attr(locuscompare_obj, "easycoloc_lead_label_panels")
assert_true(
  length(direct_label_panels) == 3L && all(direct_label_panels),
  "direct locuscomparer backend should retain lead-SNP repel labels in all three panels"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_gene_track_included")),
  "direct locuscomparer backend should include the aligned gene track"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_gene_track_column_aligned")),
  "direct locuscomparer backend should physically align the gene track with the regional column"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_gene_track_shared_x_axis")),
  "direct locuscomparer backend should use the gene track as the shared genomic x axis"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_scatter_spans_full_height")),
  "direct locuscomparer backend should not leave an empty lower-left grid cell"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_gene_interval_overlap_free")) &&
    isTRUE(attr(locuscompare_obj, "easycoloc_gene_label_overlap_free")),
  "direct locuscomparer gene track should pass interval and label collision checks"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_point_outlines_removed")),
  "direct locuscomparer backend should remove association-point outlines"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_lead_labels_background_protected")),
  "direct locuscomparer backend should protect lead labels from point occlusion"
)
assert_true(
  isTRUE(attr(locuscompare_obj, "easycoloc_official_visual_grammar")),
  "direct backend should retain the official locuscomparer visual grammar"
)

internal_locuscompare_obj <- make_combined_coloc_plot(
  df = example_data$merged_data,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set,
  gtf_path = plot_inputs$gtf_path,
  backend = "internal",
  title1 = "GWAS",
  title2 = "QTL"
)
assert_true(
  identical(attr(internal_locuscompare_obj, "easycoloc_plot_backend"), "internal"),
  "internal backend should be explicitly selectable"
)
assert_true(
  isTRUE(attr(internal_locuscompare_obj, "easycoloc_gene_track_included")),
  "internal backend should include the aligned gene track"
)

tmp_lc_pdf <- tempfile(fileext = ".pdf")
tmp_lc_png <- tempfile(fileext = ".png")
lc_save <- save_plot_with_fallback(
  plot_obj = locuscompare_obj,
  pdf_path = tmp_lc_pdf,
  png_path = tmp_lc_png,
  plot_width = LOCUSCOMPARE_LAYOUT$width_in,
  plot_height = LOCUSCOMPARE_LAYOUT$height_in,
  png_scale = 1,
  png_dpi = 200
)
assert_true(file.exists(lc_save$path), "locuscompare plot output file was not created")
assert_true(file.info(lc_save$path)$size > 5000, "locuscompare plot output file is suspiciously small")

tmp_internal_lc_pdf <- tempfile(fileext = ".pdf")
tmp_internal_lc_png <- tempfile(fileext = ".png")
internal_lc_save <- save_plot_with_fallback(
  plot_obj = internal_locuscompare_obj,
  pdf_path = tmp_internal_lc_pdf,
  png_path = tmp_internal_lc_png,
  plot_width = LOCUSCOMPARE_LAYOUT$width_in,
  plot_height = LOCUSCOMPARE_LAYOUT$height_in,
  png_scale = 1,
  png_dpi = 200
)
assert_true(file.exists(internal_lc_save$path), "internal locuscompare output file was not created")
assert_true(file.info(internal_lc_save$path)$size > 5000, "internal locuscompare output file is suspiciously small")

tmp_lc_json <- tempfile(fileext = ".json")
write_locuscompare_plot_json(
  df = example_data$merged_data,
  output_file = tmp_lc_json,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set,
  meta = list(gwas_id = "SMOKE_GWAS", qtl_id = example_data$qtl_type, pp4 = 0.92)
)
lc_json <- jsonlite::read_json(tmp_lc_json, simplifyVector = FALSE)
assert_true(is.list(lc_json$meta), "locuscompare JSON missing meta")
assert_true(length(lc_json$points) > 0, "locuscompare JSON missing points")
assert_true(identical(lc_json$meta$lead_snp, example_data$lead_snp), "locuscompare JSON lead SNP mismatch")

dir.create("plots", recursive = TRUE, showWarnings = FALSE)
dir.create("plot_data", recursive = TRUE, showWarnings = FALSE)
reference_lc_pdf <- file.path("plots", "reference_locuscomparer_example.pdf")
reference_lc_png <- file.path("plots", "reference_locuscomparer_example.png")
demo_lc_pdf <- file.path("plots", "demo_locuscompare.pdf")
demo_lc_png <- file.path("plots", "demo_locuscompare.png")
demo_lc_json <- file.path("plot_data", "demo_locuscompare.json")
reference_lc_save <- save_plot_with_fallback(
  plot_obj = locuscompare_obj,
  pdf_path = reference_lc_pdf,
  png_path = reference_lc_png,
  plot_width = LOCUSCOMPARE_LAYOUT$width_in,
  plot_height = LOCUSCOMPARE_LAYOUT$height_in,
  png_scale = 1,
  png_dpi = 300
)
demo_lc_save <- save_plot_with_fallback(
  plot_obj = locuscompare_obj,
  pdf_path = demo_lc_pdf,
  png_path = demo_lc_png,
  plot_width = LOCUSCOMPARE_LAYOUT$width_in,
  plot_height = LOCUSCOMPARE_LAYOUT$height_in,
  png_scale = 1,
  png_dpi = 300
)
write_locuscompare_plot_json(
  df = example_data$merged_data,
  output_file = demo_lc_json,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set,
  meta = list(gwas_id = "DEMO_GWAS", qtl_id = example_data$qtl_type, pp4 = 0.92)
)
demo_prepared <- prepare_locuscompare_plot_data(
  df = example_data$merged_data,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set
)
demo_layout_validation <- write_locuscompare_layout_validation(
  plot_obj = locuscompare_obj,
  output_file = file.path("plot_data", "demo_locuscompare_layout_validation.json"),
  png_path = demo_lc_png,
  pdf_path = demo_lc_pdf,
  df = example_data$merged_data,
  ld_df = example_data$ld_df,
  lead_snp = example_data$lead_snp,
  credible_set = example_data$credible_set,
  width_in = LOCUSCOMPARE_LAYOUT$width_in,
  height_in = LOCUSCOMPARE_LAYOUT$height_in,
  dpi = 300,
  reference_png = reference_lc_png,
  reference_pdf = reference_lc_pdf
)
assert_true(identical(demo_layout_validation$status, "PASS"), "demo locuscompare layout validation failed")
assert_true(isTRUE(demo_layout_validation$gene_track$included), "demo validation should record the gene track")
assert_true(isTRUE(demo_layout_validation$gene_track$matches_regional_position_range), "gene track range should match the regional panels")
assert_true(isTRUE(demo_layout_validation$gene_track$physical_column_aligned), "gene track should align with the regional-panel x coordinates")
assert_true(isTRUE(demo_layout_validation$gene_track$shared_genomic_axis), "gene track should carry the shared regional genomic axis")
assert_true(isTRUE(demo_layout_validation$gene_track$interval_overlap_free), "gene structures should be collision-free within lanes")
assert_true(isTRUE(demo_layout_validation$gene_track$visible_label_overlap_free), "visible gene labels should be collision-free within lanes")
assert_true(isTRUE(demo_layout_validation$association_points$outlines_removed), "association points should be borderless")
assert_true(isTRUE(demo_layout_validation$lead_snp$background_protected), "lead labels should have protected backgrounds")
assert_true(isTRUE(demo_layout_validation$checks$compact_publication_canvas), "demo should use the compact publication canvas")
assert_true(isTRUE(demo_layout_validation$checks$official_locuscomparer_visual_grammar), "demo should retain official locuscomparer styling")
assert_true(isTRUE(demo_layout_validation$checks$embedded_official_r2_key), "demo should retain the embedded official r2 key")
assert_true(isTRUE(demo_layout_validation$checks$scatter_spans_full_height), "demo should not retain a blank lower-left grid cell")
assert_true(isTRUE(demo_layout_validation$checks$balanced_full_height_columns), "demo should retain balanced full-height columns")
reported_gene_labels <- c(
  unlist(demo_layout_validation$gene_track$labels_visible, use.names = FALSE),
  unlist(demo_layout_validation$gene_track$labels_hidden, use.names = FALSE)
)
assert_true(
  identical(
    demo_layout_validation$text$smoke_labels,
    any(grepl("smoke", reported_gene_labels, ignore.case = TRUE))
  ),
  "demo validation should derive smoke-label metadata from reported gene labels"
)
assert_true(file.exists(demo_lc_png), "demo locuscompare PNG was not created")
assert_true(file.exists(demo_lc_pdf), "demo locuscompare PDF was not created")
assert_true(file.exists(reference_lc_png), "reference locuscomparer PNG was not created")
assert_true(file.exists(reference_lc_pdf), "reference locuscomparer PDF was not created")
assert_true(file.exists(demo_lc_json), "demo locuscompare JSON was not created")
assert_true(file.info(demo_lc_png)$size > 5000, "demo locuscompare PNG is suspiciously small")
assert_true(file.info(demo_lc_pdf)$size > 5000, "demo locuscompare PDF is suspiciously small")
assert_true(file.info(reference_lc_png)$size > 5000, "reference locuscomparer PNG is suspiciously small")
assert_true(file.info(reference_lc_pdf)$size > 5000, "reference locuscomparer PDF is suspiciously small")
assert_true(file.info(demo_lc_json)$size > 500, "demo locuscompare JSON is suspiciously small")
assert_true(isTRUE(demo_prepared$ld_available), "demo locuscompare LD should be available")
assert_true(sum(as.character(demo_prepared$data$ld_bin) == "missing LD") == 0, "demo locuscompare should not be dominated by missing LD")

cat("[SMOKE] plotting smoke test passed\n")
cat("[SMOKE] output:", output_path, "\n")
cat("[SMOKE] locuscompare output:", lc_save$path, "\n")
cat("[SMOKE] internal locuscompare output:", internal_lc_save$path, "\n")
cat("[SMOKE] reference locuscomparer PNG:", reference_lc_png, "\n")
cat("[SMOKE] reference locuscomparer PDF:", reference_lc_pdf, "\n")
cat("[SMOKE] demo locuscompare PNG:", demo_lc_png, "\n")
cat("[SMOKE] demo locuscompare PDF:", demo_lc_pdf, "\n")
cat("[SMOKE] demo locuscompare JSON:", demo_lc_json, "\n")
