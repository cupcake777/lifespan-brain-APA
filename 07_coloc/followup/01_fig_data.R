#!/usr/bin/env Rscript
# Figure source tables for PRDX6 and HID1: sample-level data, coloc and SMR evidence, locuscompare JSONs.
# Inputs from data/<module>/ under APA_ROOT; outputs to data/07_coloc/.
if (!nzchar(Sys.getenv("APA_ROOT")) || !nzchar(Sys.getenv("APA_HPC_HOME"))) stop("Set APA_ROOT and APA_HPC_HOME (source config/paths.sh)")

# --- load packages ---
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(coloc)
})

# --- helpers ---
read_matrix_row <- function(path, key, compressed = FALSE) {
  reader <- if (compressed) "gzip -cd" else "cat"
  cmd <- sprintf(
    "{ %s %s | head -n 1; %s %s | grep -F -m1 -- %s; }",
    reader, shQuote(path), reader, shQuote(path), shQuote(key)
  )
  x <- fread(cmd = cmd, check.names = FALSE)
  stopifnot(nrow(x) == 1L)
  x
}

long_row <- function(x, value_name, first_columns = 1L) {
  id_columns <- names(x)[seq_len(first_columns)]
  melt(x, id.vars = id_columns, variable.name = "sample", value.name = value_name,
       variable.factor = FALSE)[,
    c("sample", value_name), with = FALSE
  ]
}

# --- config ---
root <- Sys.getenv("APA_ROOT")
out_dir <- file.path(root, "data/07_coloc")
bam_dir <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/bam")  # cluster-only BAMs

events <- data.table(
  gene = c("PRDX6", "HID1"),
  event = c(
    "ENST00000470017.1|PRDX6|chr1:173487735-173488015|+",
    "ENST00000528902.1|HID1|chr17:74961597-74961989|-"
  ),
  variant = c("chr1:173460856:G:A", "chr17:74961620:T:C"),
  tpm_gene_id = c("ENSG00000117592.9", "ENSG00000167861.16"),
  chrom = c("chr1", "chr17"),
  window_start = c(173485735L, 74959597L),
  window_end = c(173490015L, 74963989L)
)

bam_samples <- sub("\\.bam$", "", list.files(bam_dir, pattern = "\\.bam$"))
genes <- c("SCAF1", "NRBP2", "ATP2B1", "MTRFR", "OGFOD2", "ARL3")

# --- R2a: genotypes and exact model phenotypes ---
plink_out <- file.path(out_dir, "r2_genotypes")
system2("plink2", c(
  "--pfile", file.path(root, "data/04_qtl/genotype"),
  "--extract", file.path(root, "data/04_qtl/r2_variants.txt"),
  "--export", "A", "--out", plink_out
))
dosage <- fread(paste0(plink_out, ".raw"))
setnames(dosage, "#FID", "FID", skip_absent = TRUE)
dosage_columns <- setdiff(names(dosage), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
stopifnot(length(dosage_columns) == 2L)

meta <- fread(file.path(root, "data/04_qtl/CHB.all.meta"))
tpm <- rbindlist(lapply(seq_len(nrow(events)), function(i) {
  x <- read_matrix_row(file.path(root, "raw/all.912.expr.tpm"), events$tpm_gene_id[i])
  y <- long_row(x, "tpm")
  y[, gene := events$gene[i]]
  y
}))
tpm <- tpm[, .SD[1L], by = .(sample, gene)]

corrected <- rbindlist(lapply(seq_len(nrow(events)), function(i) {
  x <- read_matrix_row(file.path(root, "data/05_gxs/after.combat.txt"), events$event[i])
  y <- long_row(x, "corrected_logit_pdui")
  y[, gene := events$gene[i]]
  y
}))

prdx6_model <- long_row(
  read_matrix_row(file.path(root, "data/05_gxs/pre.phenotype.combat.int.txt"), events[gene == "PRDX6", event]),
  "model_phenotype"
)
prdx6_model[, `:=`(gene = "PRDX6", model_source = "GxStage full-cohort inverse-normal phenotype")]

hid1_model <- rbindlist(lapply(c("Prenatal", "Postnatal"), function(stage) {
  path <- file.path(root, sprintf("data/04_qtl/phenotype.%s.GeneRange.bed.gz", stage))
  x <- read_matrix_row(path, events[gene == "HID1", event], compressed = TRUE)
  y <- long_row(x, "model_phenotype", first_columns = 4L)
  y[, `:=`(gene = "HID1", model_source = sprintf("stage-specific apaQTL %s inverse-normal phenotype", tolower(stage)))]
  y
}))

model_values <- rbindlist(list(prdx6_model, hid1_model), use.names = TRUE)
all_r2 <- merge(model_values, corrected, by = c("sample", "gene"), all.x = TRUE)
all_r2 <- merge(all_r2, tpm, by = c("sample", "gene"), all.x = TRUE)
all_r2 <- merge(all_r2, meta[, .(sample, stage = LifeStage, age, batch)], by = "sample", all.x = TRUE)

for (i in seq_len(nrow(events))) {
  variant_column <- dosage_columns[startsWith(dosage_columns, events$variant[i])]
  stopifnot(length(variant_column) == 1L)
  counted_allele <- sub(paste0("^", events$variant[i], "_"), "", variant_column)
  alt_allele <- tstrsplit(events$variant[i], ":", fixed = TRUE)[[4L]]
  d <- dosage[, .(genotype_sample = IID, exported_dosage = get(variant_column))]
  d[, genotype_dosage := if (counted_allele == alt_allele) exported_dosage else 2 - exported_dosage]
  d[, `:=`(exported_counted_allele = counted_allele, dosage_allele = alt_allele)]
  d[, sample := genotype_sample]
  d[, sample := sub("_rep$", "", sample)]
  x <- merge(all_r2[gene == events$gene[i]], d, by = "sample", all.x = TRUE)
  x[, `:=`(event = events$event[i], variant = events$variant[i])]
  setcolorder(x, c("gene", "event", "variant", "sample", "genotype_sample", "genotype_dosage",
                   "dosage_allele", "exported_dosage", "exported_counted_allele",
                   "stage", "age", "batch", "model_phenotype", "corrected_logit_pdui", "tpm", "model_source"))
  fwrite(x, file.path(out_dir, sprintf("R2a_%s_sample_data.tsv", events$gene[i])), sep = "\t", na = "NA")
}

# --- R2b: select representatives among supplied BAMs ---
r2 <- rbindlist(lapply(events$gene, function(gene_name) {
  fread(file.path(out_dir, sprintf("R2a_%s_sample_data.tsv", gene_name)), na.strings = "NA")
}))
r2[, genotype_group := pmax(0L, pmin(2L, round(genotype_dosage)))]
r2[, stage_median_tpm := median(tpm, na.rm = TRUE), by = .(gene, stage)]
r2[, phenotype_group_mean := mean(model_phenotype, na.rm = TRUE), by = .(gene, stage, genotype_group)]
candidates <- r2[sample %chin% bam_samples & tpm >= stage_median_tpm]
candidates[, phenotype_distance := abs(model_phenotype - phenotype_group_mean)]
setorder(candidates, gene, stage, genotype_group, phenotype_distance, sample)
selected <- candidates[, .SD[1L], by = .(gene, stage, genotype_group)]
selected <- merge(selected, events[, .(gene, chrom, window_start, window_end)], by = "gene")
selected[, bam := file.path(bam_dir, paste0(sample, ".bam"))]
selected[, original_selection_change := "not assessable: wig_sample_picks_draft.tsv absent"]
fwrite(selected[, .(gene, event, variant, stage, genotype_group, sample, genotype_dosage,
                    model_phenotype, phenotype_group_mean, age, batch, tpm, stage_median_tpm,
                    chrom, window_start, window_end, bam, original_selection_change)],
       file.path(out_dir, "R2b_selected_bam_samples.tsv"), sep = "\t", na = "NA", quote = FALSE)

# --- R3: coloc and SMR evidence ---
summary_data <- fread(file.path(root, "data/07_coloc/06_summary_plot_data.tsv"))
coloc_rows <- summary_data[gene %chin% genes]
fwrite(coloc_rows, file.path(out_dir, "R3_coloc_susie_summary.tsv"), sep = "\t", na = "NA")

smr <- fread(file.path(root, "data/07_coloc/smr_all_results.tsv"),
             select = c("probeID", "Gene", "topSNP", "p_SMR", "p_HEIDI", "nsnp_HEIDI",
                        "qtl_type", "gwas_id", "population", "p_SMR_BH", "p_SMR_bonferroni"))
smr_rows <- smr[Gene %chin% genes]
fwrite(smr_rows, file.path(out_dir, "R3_smr_heidi_results.tsv"), sep = "\t", na = "NA")

signal_keys <- coloc_rows[, .(key, GWAS_ID, stage, gene, locus)]
signal_keys[, probeID := vapply(key, function(k) {
  x <- readRDS(file.path(out_dir, paste0(k, "_coloc.rds")))
  unique(x$merged_data$phenotype_id)[1L]
}, character(1))]
signal_smr <- merge(
  signal_keys,
  smr,
  by.x = c("probeID", "GWAS_ID", "stage"),
  by.y = c("probeID", "gwas_id", "qtl_type"),
  all.x = TRUE
)
setnames(signal_smr, "stage", "qtl_type")
fwrite(signal_smr, file.path(out_dir, "R3_target_signal_smr_heidi.tsv"), sep = "\t", na = "NA")

scaf1_locuscompare <- file.path(root, "data/figures/fig5/EUR_SCZ_pre_chr19_rs7251_SCAF1_locuscompare.json")
scaf1_json <- fromJSON(scaf1_locuscompare, simplifyDataFrame = TRUE)
scaf1_check <- as.data.table(scaf1_json$points)[snp %chin% c("rs7047", "rs7251")]
setnames(scaf1_check, "snp_pp_h4", "json_snp_pp_h4")
scaf1_rds <- readRDS(file.path(out_dir, "EUR_SCZ_prenatal_chr19_49659652_rs7251_SCAF1_coloc.rds"))
scaf1_data <- scaf1_rds$merged_data
scaf1_abf <- coloc.abf(
  list(beta = scaf1_data$BETA.gwas, varbeta = scaf1_data$SE.gwas^2,
       type = "cc", s = 0.4086, N = scaf1_data$N, MAF = scaf1_data$EAF,
       snp = scaf1_data$SNPID),
  list(beta = scaf1_data$BETA.qtl, varbeta = scaf1_data$SE.qtl^2,
       type = "quant", N = 189, MAF = rep(0.1, nrow(scaf1_data)),
       snp = scaf1_data$SNPID)
)
exact_pp <- as.data.table(scaf1_abf$results)[snp %chin% c("rs7047", "rs7251"), .(snp, recomputed_snp_pp_h4 = SNP.PP.H4)]
scaf1_check <- merge(scaf1_check, exact_pp, by = "snp", all.x = TRUE)
scaf1_check[, `:=`(
  coloc_PP4 = scaf1_abf$summary[["PP.H4.abf"]],
  actual_locus_index = scaf1_json$meta$locus,
  qtl_maf_contract = 0.1,
  qtl_n = 189L,
  gwas_case_fraction = 0.4086
)]
fwrite(scaf1_check, file.path(out_dir, "R3_SCAF1_index_check.tsv"), sep = "\t", na = "NA")

nr <- fromJSON(file.path(out_dir, "EUR_SCZ_prenatal_chr8_143788531_rs11136325_NRBP2_locuscompare.json"), simplifyDataFrame = FALSE)
nr$points[[length(nr$points) + 1L]] <- list(
  snp = "rs11136325", chrom = "chr8", pos = 143788531L,
  gwas_p = NULL, qtl_p = NULL, gwas_logp = NULL, qtl_logp = NULL,
  r2 = NULL, ld_bin = "missing LD", is_lead = FALSE,
  is_credible = FALSE, snp_pp_h4 = NULL,
  availability = "index SNP absent from merged QTL data; association fields unavailable"
)
nr$meta$requested_index_snp <- "rs11136325"
nr$meta$requested_index_present <- TRUE
nr$meta$requested_index_has_qtl <- FALSE
write_json(nr, file.path(out_dir, "Fig5e_NRBP2_locuscompare.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
file.copy(scaf1_locuscompare, file.path(out_dir, "Fig5d_SCAF1_locuscompare.json"), overwrite = TRUE)

eur_ad <- summary_data[GWAS_ID == "EUR_AD"]
eur_ad[, `:=`(
  configured_dataset = "EUR_AD",
  configured_source = "AD_hg19_eur.txt.gz",
  study_mapping = "Kunkle 2019 per requested EUR_AD/EUR_AD2 naming convention"
)]
fwrite(eur_ad, file.path(out_dir, "R3_EUR_AD_provenance.tsv"), sep = "\t", na = "NA")
