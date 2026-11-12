#!/usr/bin/env Rscript
# eQTL colocalization at APA-QTL candidate phenotypes (optional; run_analysis.sh gates it).
# Writes data/07_coloc/eqtl_coloc.tsv.
for (.apa_v in c("APA_ROOT", "APA_HPC_HOME")) if (!nzchar(Sys.getenv(.apa_v))) stop("Set ", .apa_v, " (source config/paths.sh)")

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# --- helpers ---

summarize_empty <- function(candidates, reason) {
  data.table(
    candidate_id = candidates$candidate_id,
    apa_phenotype_id = candidates$phenotype_id,
    egene = NA_character_,
    stage = candidates$stage,
    n_snps = 0L,
    PP.H0 = NA_real_,
    PP.H1 = NA_real_,
    PP.H2 = NA_real_,
    PP.H3 = NA_real_,
    PP.H4 = NA_real_,
    top_apa_snp = NA_character_,
    top_eqtl_snp = NA_character_,
    eQTL_coloc = FALSE,
    coloc_status = reason
  )
}

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
apa_hpc <- Sys.getenv("APA_HPC_HOME")
root <- file.path(apa_root, "data/07_coloc")
candidates_path <- file.path(root, "candidates_prioritized.tsv")
output_path <- file.path(root, "eqtl_coloc.tsv")
pp4_threshold <- 0.75
min_shared_snps <- 100

nominal_paths <- list(
  prenatal = list(
    apa = paste0(apa_hpc, "/project/work/apa/01_QTL/nominal/prenatal.cis_qtl_pairs.txt.gz"),
    eqtl = paste0(apa_hpc, "/project/work/effect/EXP/stage/data/nominal/prenatal.cis_qtl_pairs.txt.gz")
  ),
  postnatal = list(
    apa = paste0(apa_hpc, "/project/work/apa/01_QTL/nominal/postnatal.cis_qtl_pairs.txt.gz"),
    eqtl = paste0(apa_hpc, "/project/work/effect/EXP/stage/data/nominal/postnatal.cis_qtl_pairs.txt.gz")
  )
)

# --- coloc ---

candidates <- fread(candidates_path)
rows <- list()
for (stage_name in names(nominal_paths)) {
  stage_candidates <- candidates[stage == stage_name]
  if (nrow(stage_candidates) == 0) next
  apa <- fread(nominal_paths[[stage_name]]$apa, sep = "\t", select = c("phenotype_id", "variant_id", "af", "ma_samples", "slope", "slope_se"))
  apa <- apa[phenotype_id %chin% unique(stage_candidates$phenotype_id)]
  if (nrow(apa) == 0) {
    rows[[length(rows) + 1]] <- summarize_empty(stage_candidates, "no_apa_nominal_rows")
    next
  }
  eqtl_all <- fread(nominal_paths[[stage_name]]$eqtl, sep = "\t", select = c("phenotype_id", "variant_id", "af", "ma_samples", "slope", "slope_se"))
  for (i in seq_len(nrow(stage_candidates))) {
    cand <- stage_candidates[i]
    apa_i <- apa[phenotype_id == cand$phenotype_id]
    if (nrow(apa_i) == 0) {
      rows[[length(rows) + 1]] <- summarize_empty(cand, "no_apa_nominal_rows")
      next
    }
    shared_eqtl <- eqtl_all[variant_id %chin% apa_i$variant_id]
    if (nrow(shared_eqtl) == 0) {
      rows[[length(rows) + 1]] <- summarize_empty(cand, "no_eqtl_overlap")
      next
    }
    for (egene in unique(shared_eqtl$phenotype_id)) {
      eqtl_i <- shared_eqtl[phenotype_id == egene]
      common <- intersect(apa_i$variant_id, eqtl_i$variant_id)
      if (length(common) < min_shared_snps) next
      apa_m <- apa_i[variant_id %chin% common][match(common, variant_id)]
      eqtl_m <- eqtl_i[variant_id %chin% common][match(common, variant_id)]
      coloc_res <- coloc.abf(
        dataset1 = list(beta = apa_m$slope, varbeta = apa_m$slope_se^2, snp = common, type = "quant", N = max(apa_m$ma_samples, na.rm = TRUE), MAF = apa_m$af),
        dataset2 = list(beta = eqtl_m$slope, varbeta = eqtl_m$slope_se^2, snp = common, type = "quant", N = max(eqtl_m$ma_samples, na.rm = TRUE), MAF = eqtl_m$af),
        p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
      )
      summary <- coloc_res$summary
      rows[[length(rows) + 1]] <- data.table(
        candidate_id = cand$candidate_id,
        apa_phenotype_id = cand$phenotype_id,
        egene = egene,
        stage = stage_name,
        n_snps = length(common),
        PP.H0 = summary[["PP.H0.abf"]],
        PP.H1 = summary[["PP.H1.abf"]],
        PP.H2 = summary[["PP.H2.abf"]],
        PP.H3 = summary[["PP.H3.abf"]],
        PP.H4 = summary[["PP.H4.abf"]],
        top_apa_snp = apa_m$variant_id[which.max(abs(apa_m$slope / apa_m$slope_se))],
        top_eqtl_snp = eqtl_m$variant_id[which.max(abs(eqtl_m$slope / eqtl_m$slope_se))],
        eQTL_coloc = summary[["PP.H4.abf"]] >= pp4_threshold,
        coloc_status = "ok"
      )
    }
  }
}

# --- write ---

out <- if (length(rows)) rbindlist(rows, fill = TRUE) else summarize_empty(candidates[0], "no_candidates")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
fwrite(out, output_path, sep = "\t")
