#!/usr/bin/env Rscript
# Age trajectories of protein-minus-mRNA against PDUI for C1 and C2 gene groups.
# Reads 03_omics and 02_traj tables; writes T3_traj_protein_mrna.tsv to data/02_traj/.
stopifnot(nzchar(Sys.getenv(c("APA_ROOT", "APA_HPC_HOME"))))
suppressPackageStartupMessages({
  library(data.table)
  library(splines)
  library(readxl)
})

# --- helpers ---
fit_curve <- function(y, design, grid) {
  coefficients <- lm.fit(design, y)$coefficients
  grid %*% coefficients
}

# --- config ---
omics_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "03_omics")
traj_dir <- file.path(Sys.getenv("APA_ROOT"), "data", "02_traj")
out <- file.path(traj_dir, "T3_traj_protein_mrna.tsv")
ref <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment")
permutations <- 10000L
set.seed(20261017L)

# --- loading data ---
apa <- readRDS(file.path(omics_dir, "apa_matched.rds"))
rna <- readRDS(file.path(omics_dir, "rna_matched.rds"))
protein <- readRDS(file.path(omics_dir, "protein_matched.rds"))
samples <- fread(file.path(omics_dir, "matched_sample_manifest.tsv"))
stopifnot(identical(colnames(apa), samples$sample_id),
          identical(colnames(rna), samples$sample_id),
          identical(colnames(protein), samples$sample_id))
u <- fread(file.path(omics_dir, "01_analysis_universe.tsv"),
           select = c("site_id", "gencode_gene_base", "gene_symbol",
                      "protein_accession", "rna_available", "protein_available",
                      "gene_symbol_unambiguous", "protein_mapping_unambiguous"))
u <- u[rna_available & protein_available & gene_symbol_unambiguous &
         protein_mapping_unambiguous]
genes <- fread(file.path(traj_dir, "05_gene_aggregation.tsv"))[
  gene_state %in% c("C1-spec", "C2-spec")]
mapping <- unique(u[, .(gencode_gene_base, gene_symbol, protein_accession)])
mapping <- merge(mapping, genes[, .(gene_id_base, gene_state)],
                 by.x = "gencode_gene_base", by.y = "gene_id_base")
stopifnot(!anyDuplicated(mapping$gencode_gene_base))

mito <- read_excel(file.path(ref, "Human.MitoCarta3.0.xls"),
                   sheet = "C MitoPathways")
oxphos <- trimws(strsplit(mito$Genes[which(mito$MitoPathway == "OXPHOS subunits")[1L]],
                          ",", fixed = TRUE)[[1L]])
proteasome <- fread(file.path(ref, "HGNC_Proteasome.tsv"))$`Approved symbol`
mapping[, group := fifelse(gene_state == "C1-spec", "C1",
                   fifelse(gene_symbol %in% oxphos, "C2_OXPHOS",
                   fifelse(gene_symbol %in% proteasome, "C2_proteasome", "C2_other")))]

# --- gene-level pdui and matched protein/rna ---
events <- u[gencode_gene_base %in% mapping$gencode_gene_base]
event_index <- match(events$site_id, rownames(apa))
stopifnot(!anyNA(event_index))
apa_gene <- rowsum(plogis(apa[event_index, , drop = FALSE]),
                   group = events$gencode_gene_base, reorder = TRUE)
apa_gene <- apa_gene / as.integer(table(events$gencode_gene_base)[rownames(apa_gene)])
rna_index <- match(mapping$gencode_gene_base, sub("\\.[0-9]+$", "", rownames(rna)))
protein_index <- match(mapping$protein_accession, rownames(protein))
stopifnot(!anyNA(rna_index), !anyNA(protein_index),
          !anyNA(match(mapping$gencode_gene_base, rownames(apa_gene))))
apa_y <- t(apa_gene[mapping$gencode_gene_base, , drop = FALSE])
rna_y <- t(rna[rna_index, , drop = FALSE])
protein_y <- t(protein[protein_index, , drop = FALSE])

# --- trajectory correlations ---
stage_index <- list(all = seq_len(nrow(samples)),
                    prenatal = which(samples$life_stage == "Prenatal"),
                    postnatal = which(samples$life_stage == "Postnatal"))
results <- rbindlist(lapply(names(stage_index), function(stage) {
  i <- stage_index[[stage]]
  age <- log10(samples$age_years[i] + 1)
  basis <- ns(age, df = 2L)
  grid_age <- sort(unique(quantile(age, seq(0, 1, length.out = 101L))))
  grid_basis <- predict(basis, grid_age)
  sex <- as.integer(samples$sex[i] == "M")
  rin <- samples$rna_rin[i] - mean(samples$rna_rin[i])
  design <- cbind(1, sex, rin, basis)
  grid <- cbind(1, 0, 0, grid_basis)
  protein_design <- cbind(1, basis)
  protein_grid <- cbind(1, grid_basis)
  pd <- fit_curve(apa_y[i, , drop = FALSE], design, grid)
  mr <- fit_curve(rna_y[i, , drop = FALSE], design, grid)
  pr <- fit_curve(protein_y[i, , drop = FALSE], protein_design, protein_grid)
  ratio <- pr - mr
  correlation <- vapply(seq_len(ncol(pd)), function(j) {
    if (sd(pd[, j]) == 0 || sd(ratio[, j]) == 0) return(NA_real_)
    cor(pd[, j], ratio[, j])
  }, numeric(1))
  cbind(mapping[, .(gencode_gene_base, gene_symbol, gene_state, group)],
        data.table(stage = stage, n_samples = length(i),
                   n_grid = length(grid_age), correlation = correlation))
}))
results <- results[is.finite(correlation)]

# --- empirical comparisons ---
summary <- results[, .(n_genes = .N, n_samples = first(n_samples),
                       n_grid = first(n_grid),
                       median_correlation = median(correlation),
                       q25_correlation = quantile(correlation, 0.25),
                       q75_correlation = quantile(correlation, 0.75)),
                   by = .(stage, group)]
for (j in seq_len(nrow(summary))) {
  row <- summary[j]
  for (null_group in c("C2_other", "C1")) {
    control <- results[stage == row$stage & group == null_group]
    p <- if (row$group == "C1" || row$group == null_group) NA_real_ else {
      null_medians <- replicate(permutations,
        median(sample(control$correlation, row$n_genes, replace = FALSE)))
      (1 + sum(null_medians <= row$median_correlation)) / (permutations + 1)
    }
    set(summary, j, paste0("empirical_p_vs_", null_group), p)
  }
}
results[, `:=`(n_genes = NA_integer_, median_correlation = NA_real_,
               q25_correlation = NA_real_, q75_correlation = NA_real_,
               empirical_p_vs_C2_other = NA_real_, empirical_p_vs_C1 = NA_real_,
               record_type = "gene")]
summary[, `:=`(gencode_gene_base = NA_character_, gene_symbol = NA_character_,
                gene_state = NA_character_, correlation = NA_real_,
                record_type = "summary")]
fwrite(rbindlist(list(results, summary), fill = TRUE), out, sep = "\t", na = "NA")
