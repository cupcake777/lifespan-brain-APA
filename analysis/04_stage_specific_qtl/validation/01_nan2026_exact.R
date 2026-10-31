#!/usr/bin/env Rscript
# Exact Nan 2026 replication: interaction leads (R4.1) and stage-specific apaQTL genes (R4.2).
# Nominal Nan 2026 rows are matched on transcript, position and allele pair by streaming awk.
# Writes R4_* tables to data/04_qtl/.

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# --- config ---
APA_ROOT <- Sys.getenv("APA_ROOT")
if (!nzchar(APA_ROOT)) stop("Set APA_ROOT (source config/paths.sh)")
out_dir <- file.path(APA_ROOT, "data", "04_qtl")
nominal_dir <- file.path(APA_ROOT, "ref", "nan2026")
workbook <- file.path(APA_ROOT, "ref", "nan2026", "41588_2026_2758_MOESM3_ESM.xlsx")
interaction_path <- file.path(APA_ROOT, "data", "05_gxs", "interaction_hc3_core_groups.tsv")
prenatal_path <- file.path(out_dir, "prenatal.cis_qtl.txt.gz")
postnatal_path <- file.path(out_dir, "postnatal.cis_qtl.txt.gz")
gtf_path <- file.path(APA_ROOT, "ref", "gencode", "gencode.v40.annotation.gtf.gz")
cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc")

# --- helpers ---
split_variant <- function(x) {
  parts <- tstrsplit(x, ":", fixed = TRUE)
  data.table(chromosome = parts[[1L]], position = as.integer(parts[[2L]]),
             reference = parts[[3L]], alternate = parts[[4L]])
}

allele_key <- function(transcript, position, a1, a2) {
  paste(transcript, position, pmin(a1, a2), pmax(a1, a2), sep = ":")
}

read_stage_leads <- function(path, prefix) {
  x <- fread(path, select = c("phenotype_id", "variant_id", "group_id", "qval", "af", "pval_nominal", "slope"))
  setorder(x, group_id, pval_nominal)
  x <- x[, .SD[1L], by = group_id]
  setnames(x, setdiff(names(x), "group_id"), paste0(prefix, "_", setdiff(names(x), "group_id")))
  x
}

model_table <- function(full, reduced, model_id) {
  lrt <- anova(reduced, full, test = "Chisq")
  coefficients <- as.data.table(summary(full)$coefficients, keep.rownames = "term")
  setnames(coefficients, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
           c("estimate", "std_error", "statistic", "p_value"))
  overall <- data.table(
    term = "overall_class", estimate = NA_real_, std_error = NA_real_,
    statistic = lrt$Deviance[2L], p_value = lrt$`Pr(>Chi)`[2L], test = "likelihood-ratio"
  )
  coefficients[, test := "Wald"]
  rbind(overall, coefficients, fill = TRUE)[, `:=`(
    model_id = model_id, n_observations = nobs(full),
    n_positive = sum(model.response(model.frame(full))),
    formula = paste(deparse(formula(full)), collapse = " ")
  )]
}

# --- local interaction pairs ---
interaction <- fread(interaction_path)
interaction[, transcript_id := tstrsplit(lead_phenotype_id, "|", fixed = TRUE)[[1L]]]
interaction_variant <- split_variant(interaction$lead_variant_id)
interaction <- cbind(interaction, interaction_variant)
interaction[, match_key := allele_key(transcript_id, position, reference, alternate)]
interaction_local <- interaction[, .(
  analysis = "interaction", local_id = group_id, transcript_id, position,
  reference, alternate, match_key, local_maf = lead_af,
  local_slope = beta_pre, point_estimate_class
)]

# --- local separate-stage pairs ---
stage <- merge(read_stage_leads(prenatal_path, "prenatal"),
               read_stage_leads(postnatal_path, "postnatal"), by = "group_id", all = TRUE)
stage[, qtl_class := "Neither"]
stage[prenatal_qval <= 0.05 & (is.na(postnatal_qval) | postnatal_qval > 0.05), qtl_class := "Prenatal-specific"]
stage[(is.na(prenatal_qval) | prenatal_qval > 0.05) & postnatal_qval <= 0.05, qtl_class := "Postnatal-specific"]
stage[prenatal_qval <= 0.05 & postnatal_qval <= 0.05, qtl_class := "Shared"]
stage <- stage[qtl_class != "Neither"]
stage[, selected_stage := fifelse(qtl_class == "Prenatal-specific", "prenatal", "postnatal")]
stage[, phenotype_id := fifelse(selected_stage == "prenatal", prenatal_phenotype_id, postnatal_phenotype_id)]
stage[, variant_id := fifelse(selected_stage == "prenatal", prenatal_variant_id, postnatal_variant_id)]
stage[, local_maf := fifelse(selected_stage == "prenatal", prenatal_af, postnatal_af)]
stage[, local_slope := fifelse(selected_stage == "prenatal", prenatal_slope, postnatal_slope)]
stage[, transcript_id := tstrsplit(phenotype_id, "|", fixed = TRUE)[[1L]]]
stage_variant <- split_variant(stage$variant_id)
stage <- cbind(stage, stage_variant)
stage[, match_key := allele_key(transcript_id, position, reference, alternate)]
stage_local <- stage[, .(
  analysis = "stage_qtl", local_id = group_id, transcript_id, position,
  reference, alternate, match_key, local_maf, local_slope, qtl_class
)]

# --- stream exact Nan nominal rows ---
all_local <- unique(rbind(
  interaction_local[, .(match_key)], stage_local[, .(match_key)]
))
key_path <- file.path(out_dir, "R4_nan2026_exact_match_keys.tsv")
fwrite(all_local, key_path, col.names = FALSE, sep = "\t")

matched_list <- lapply(cell_types, function(cell_type) {
  source_path <- file.path(nominal_dir, paste0(cell_type, ".3aQTL.sumstats.txt.gz"))
  stopifnot(file.exists(source_path), file.info(source_path)$size > 0L)
  output_path <- file.path(out_dir, paste0(".R4_", cell_type, "_matches.tsv"))
  awk <- paste0(
    "BEGIN{FS=OFS=\"\\t\"} ",
    "NR==FNR{key[$1]=1;next} ",
    "FNR==1{print \"celltype\",$0,\"match_key\";next} ",
    "{split($6,p,\"|\"); a1=$4; a2=$5; if(a1>a2){z=a1;a1=a2;a2=z}; ",
    "k=p[1] \":\" $3 \":\" a1 \":\" a2; if(k in key) print \"", cell_type, "\",$0,k}"
  )
  command <- sprintf("gzip -cd %s | awk -f <(printf '%%s\\n' %s) %s - > %s",
                     shQuote(source_path), shQuote(awk), shQuote(key_path), shQuote(output_path))
  status <- system2("bash", c("-c", shQuote(command)))
  if (status != 0L) stop("Nominal extraction failed for ", cell_type)
  fread(output_path)
})
nominal <- rbindlist(matched_list, fill = TRUE)
nominal[, match_key := allele_key(tstrsplit(phenotype_id, "|", fixed = TRUE)[[1L]],
                                  position, Reference_allele, Effect_allele)]
fwrite(nominal, file.path(out_dir, "R4_nan2026_nominal_matches.tsv.gz"), sep = "\t", compress = "gzip")

# --- R4.1 interaction replication model ---
interaction_long <- CJ(local_row = seq_len(nrow(interaction_local)), celltype = cell_types)
interaction_long <- cbind(interaction_local[interaction_long$local_row], interaction_long[, .(celltype)])
interaction_long <- merge(interaction_long, nominal[, .(
  celltype, match_key, nan_maf = maf, nan_p = pval_nominal, nan_slope = slope
)], by = c("celltype", "match_key"), all.x = TRUE, sort = FALSE)
interaction_long[, matched := is.finite(nan_p)]
interaction_long[, replicated := nan_p < 0.05]
interaction_long[, same_direction := sign(local_slope) == sign(nan_slope)]
fwrite(interaction_long, file.path(out_dir, "R4_interaction_nan_exact_pairs.tsv"), sep = "\t")

interaction_model_data <- interaction_long[matched == TRUE]
interaction_model_data[, point_estimate_class := factor(
  point_estimate_class, levels = c("Opposing", "Amplified", "Attenuated")
)]
interaction_full <- glm(replicated ~ point_estimate_class + nan_maf + factor(celltype),
                        data = interaction_model_data, family = binomial())
interaction_reduced <- update(interaction_full, . ~ . - point_estimate_class)
interaction_model <- model_table(interaction_full, interaction_reduced,
                                 "R4.1 exact pair x cell type")
interaction_rates <- interaction_model_data[, .(
  n_tested = .N, n_replicated = sum(replicated), replication_rate = mean(replicated),
  same_direction_rate_among_replicated = mean(same_direction[replicated])
), by = point_estimate_class]
fwrite(interaction_rates, file.path(out_dir, "R4_interaction_nan_replication_rates.tsv"), sep = "\t")

# --- R4.2 separate-stage gene comparison ---
nan <- as.data.table(read_excel(workbook, sheet = "Table 11", skip = 1L))
setnames(nan, "Phenotype_id (transcript ID|gene ID|chromosome|strand)", "phenotype_id")
nan <- nan[!is.na(SNP) & !is.na(gene) & !is.na(celltype)]
nan_genes <- unique(nan$gene)
gtf <- fread(cmd = sprintf("gzip -cd %s | awk 'BEGIN{FS=OFS=\"\\t\"} $3==\"gene\"{print $9}'", shQuote(gtf_path)),
             sep = "\t", header = FALSE, quote = "", fill = TRUE)
gtf[, gene_id := sub(".*gene_id \\\"([^\\\"]+)\\\".*", "\\1", V1)]
gtf[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
gtf[, gene_symbol := sub(".*gene_name \\\"([^\\\"]+)\\\".*", "\\1", V1)]
gene_map <- unique(gtf[, .(gene_id, gene_symbol)], by = "gene_id")
stage_gene <- stage[, .(group_id, qtl_class, maf = local_maf)]
stage_gene[, gene_symbol := fifelse(startsWith(group_id, "ENSG"),
                                    gene_map$gene_symbol[match(group_id, gene_map$gene_id)], group_id)]
stage_gene[, replicated := gene_symbol %chin% nan_genes]
stage_gene[, qtl_class := factor(qtl_class, levels = c("Postnatal-specific", "Prenatal-specific", "Shared"))]
stage_full <- glm(replicated ~ qtl_class + maf, data = stage_gene, family = binomial())
stage_reduced <- update(stage_full, . ~ . - qtl_class)
stage_model <- model_table(stage_full, stage_reduced, "R4.2 gene any Nan cell type")
stage_rates <- stage_gene[, .(
  n_genes = .N, n_replicated = sum(replicated), replication_rate = mean(replicated)
), by = qtl_class]
fwrite(stage_rates, file.path(out_dir, "R4_stage_nan_replication_rates.tsv"), sep = "\t")

models <- rbind(interaction_model, stage_model, fill = TRUE)
setcolorder(models, c("model_id", "test", "term", "estimate", "std_error", "statistic",
                      "p_value", "n_observations", "n_positive", "formula"))
fwrite(models, file.path(out_dir, "R4_nan2026_exact_model_outputs.tsv"), sep = "\t")

cat(sprintf(
  "R4 complete: interaction matches=%d/%d; R4.1 LRT P=%.8g; R4.2 LRT P=%.8g\n",
  interaction_model_data[, uniqueN(paste(local_id, celltype))], nrow(interaction_long),
  interaction_model[test == "likelihood-ratio", p_value],
  stage_model[test == "likelihood-ratio", p_value]
))
