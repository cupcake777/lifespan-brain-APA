#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

utils_files <- list.files("src", pattern = "^utils_.*\\.R$", full.names = TRUE)
invisible(lapply(utils_files, source))

cfg_global <- yaml::read_yaml("config/global.yml")
source("examples/synthetic_plot_data.R")

plot_inputs <- resolve_synthetic_plot_inputs(cfg_global)
example_data <- make_synthetic_locus()

assign("lead_SNP", example_data$lead_snp, envir = .GlobalEnv)
assign("geneSymbol", example_data$phenotype_info, envir = .GlobalEnv)
assign("plink_bfile", NULL, envir = .GlobalEnv)
assign("trait", "SMOKE_GWAS", envir = .GlobalEnv)

plot_obj <- plot_qtl_association(
  qtl_all_chrom = "CHR.qtl",
  qtl_all_pvalue = "P.qtl",
  leadSNP_DF = example_data$merged_data,
  ld_df = example_data$ld_df,
  gtf_path = plot_inputs$gtf_path,
  region_recomb = NULL,
  recomb_path = plot_inputs$recomb_path,
  show_lead_line = FALSE,
  qtl_type = example_data$qtl_type,
  phenotype_info = example_data$phenotype_info,
  credible_set = example_data$credible_set,
  plot_window_bp = 200000
)

output_dir <- "examples/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_res <- save_plot_with_fallback(
  plot_obj = plot_obj,
  pdf_path = file.path(output_dir, "synthetic_locus_demo.pdf"),
  png_path = file.path(output_dir, "synthetic_locus_demo.png"),
  plot_width = 6.6,
  plot_height = 4.9,
  png_scale = 1,
  png_dpi = 240
)

message("[DEMO] Saved minimal synthetic panel to: ", save_res$path)
