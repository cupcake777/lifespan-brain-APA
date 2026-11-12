#!/usr/bin/env Rscript
# CHB eQTL coloc for representative signal-stage pairs missing from the matched-stage table.
# Writes representative_chb_crossstage.tsv to data/07_coloc.
if (!all(nzchar(Sys.getenv(c("APA_ROOT", "APA_HPC_HOME"))))) stop("Source config/paths.sh first")

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(yaml)
})

# --- config ---

root <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa")
context <- file.path(root, "14_context")
coloc_dir <- paste0(Sys.getenv("APA_ROOT"), "/data/07_coloc")
genes <- c("SCAF1", "NRBP2", "ATP2B1", "MTRFR", "OGFOD2", "ARL3")
stages <- c("prenatal", "postnatal")
gwas_dir <- paste0(Sys.getenv("APA_HPC_HOME"), "/share_group_folder/.ref/GWAS/harmony")
output <- file.path(coloc_dir, "representative_chb_crossstage.tsv")

# --- helpers ---

strip_version <- function(x) sub("\\.[0-9]+$", "", x)
safe_maf <- function(x) {
  x <- pmin(as.numeric(x), 1 - as.numeric(x))
  x[!is.finite(x)] <- 0.05
  pmax(pmin(x, 1 - 1e-6), 1e-6)
}

# --- loading data ---

signals <- fread(file.path(coloc_dir, "apa_signal_eqtl_classification.tsv"))[
  gene %chin% genes, .(signal_id, gene, GWAS, stage, apa_event, index_chrom, index_position)
]
existing <- fread(file.path(coloc_dir, "chb_stage_eqtl_coloc_context.tsv"))
targets <- CJ(signal_id = signals$signal_id, chb_stage = stages)[
  !existing, on = .(signal_id, chb_stage)
]
targets <- merge(targets, signals, by = "signal_id")
stopifnot(nrow(targets) == 10L)

tx <- fread(file.path(root, "02_pattern/data/process/gencode.v40.tx2gene.tsv"),
            select = c("tx_name", "gene_name", "gene_id"))
tx[, tx_base := strip_version(tx_name)]
tx[, gene_id_base := strip_version(gene_id)]
targets[, tx_base := strip_version(sub("\\|.*", "", apa_event))]
targets <- merge(targets, unique(tx[, .(tx_base, gene = gene_name, gene_id_base)]),
                 by = c("tx_base", "gene"), all.x = TRUE)
stopifnot(!anyNA(targets$gene_id_base), !anyDuplicated(targets[, .(signal_id, chb_stage)]))

config <- yaml.load_file(file.path(coloc_dir, "gwas.yml"))$datasets
config <- setNames(config, vapply(config, function(x) x$id, character(1)))

# --- complete GWAS windows ---

gwas_windows <- list()
for (gwas_id in unique(targets$GWAS)) {
  gwas_file <- list.files(gwas_dir, pattern = sprintf("^%s_b(?:19|38)to38_harmonized\\.tsv$", gwas_id), full.names = TRUE)
  stopifnot(length(gwas_file) == 1L)
  regions <- unique(targets[GWAS == gwas_id, .(index_chrom, start = pmax(1L, index_position - 500000L), end = index_position + 500000L)])
  region_file <- tempfile(fileext = ".tsv")
  extracted <- tempfile(fileext = ".tsv")
  fwrite(regions, region_file, sep = "\t", col.names = FALSE)
  awk <- paste0('NR==FNR{n++;rchr[n]=$1;rstart[n]=$2;rend[n]=$3;next}',
                'FNR==1{for(j=1;j<=NF;j++){if($j=="CHR")ci=j;if($j=="POS")pi=j};print;next}',
                '{for(j=1;j<=n;j++){if($ci==rchr[j]&&$pi>=rstart[j]&&$pi<=rend[j]){print;break}}}')
  command <- sprintf("awk -F '\\t' '%s' %s %s > %s", awk, shQuote(region_file), shQuote(gwas_file), shQuote(extracted))
  stopifnot(system(command) == 0L)
  gwas_windows[[gwas_id]] <- fread(extracted, showProgress = FALSE)
  unlink(c(region_file, extracted))
}

# --- CHB opposite-stage coloc ---

gene_file <- tempfile(fileext = ".txt")
writeLines(unique(targets$gene_id_base), gene_file)
results <- list()
for (chb_stage in stages) {
  nominal <- file.path(context, ".staging_eqtl_rerun_20261015/nominal", paste0(chb_stage, ".cis_qtl_pairs.txt.gz"))
  sample_file <- file.path(context, ".staging_eqtl_rerun_20261015/pre_input", paste0("sample.", if (chb_stage == "prenatal") "Prenatal" else "Postnatal", ".final"))
  extracted <- tempfile(fileext = ".tsv")
  awk <- sprintf('BEGIN{while((getline line < "%s")>0) ids[line]=1}NR==1{print;next}{key=$1;sub(/\\.[0-9]+$/,"",key);if(key in ids)print}', gene_file)
  command <- sprintf("gzip -cd %s | awk -F '\\t' '%s' > %s", shQuote(nominal), awk, shQuote(extracted))
  stopifnot(system(command) == 0L)
  eqtl <- fread(extracted, showProgress = FALSE)
  unlink(extracted)
  eqtl[, `:=`(gene_id_base = strip_version(phenotype_id), variant_key = as.character(variant_id), maf = safe_maf(af), n = length(readLines(sample_file)))]
  setorder(eqtl, gene_id_base, variant_key, pval_nominal)
  eqtl <- eqtl[!duplicated(eqtl[, .(gene_id_base, variant_key)])]
  for (i in which(targets$chb_stage == chb_stage)) {
    signal <- targets[i]
    window <- gwas_windows[[signal$GWAS]][CHR == signal$index_chrom &
      POS >= signal$index_position - 500000L & POS <= signal$index_position + 500000L]
    cfg <- config[[signal$GWAS]]
    gwas <- data.table(variant_key = as.character(window$variant_id), beta = as.numeric(window$BETA),
                       se = as.numeric(window$SE), maf = safe_maf(window$EAF), n = as.numeric(window$N),
                       gwas_type = as.character(cfg$type), case_fraction = if (!is.null(cfg$prop)) as.numeric(cfg$prop) else NA_real_)
    gwas <- gwas[!duplicated(variant_key)]
    gene_eqtl <- eqtl[gene_id_base == signal$gene_id_base]
    fit <- list(pp4 = NA_real_, n_snps = 0L, status = "gene_not_found")
    if (nrow(gene_eqtl)) {
      common <- intersect(gwas$variant_key, gene_eqtl$variant_key)
      fit <- list(pp4 = NA_real_, n_snps = length(common), status = "fewer_than_10_shared_variants")
      if (length(common) >= 10L) {
        g <- gwas[match(common, variant_key)]
        e <- gene_eqtl[match(common, variant_key)]
        ok <- is.finite(g$beta) & is.finite(g$se) & g$se > 0 &
          is.finite(e$slope) & is.finite(e$slope_se) & e$slope_se > 0 &
          is.finite(g$maf) & is.finite(e$maf)
        g <- g[ok]
        e <- e[ok]
        fit <- list(pp4 = NA_real_, n_snps = nrow(g), status = "fewer_than_10_complete_variants")
        if (nrow(g) >= 10L) {
          d1 <- list(snp = g$variant_key, beta = g$beta, varbeta = g$se^2,
                     MAF = g$maf, N = g$n, type = g$gwas_type[1L])
          if (identical(d1$type, "cc")) d1$s <- g$case_fraction[1L]
          d2 <- list(snp = e$variant_key, beta = e$slope, varbeta = e$slope_se^2,
                     MAF = e$maf, N = e$n[1L], type = "quant")
          res <- try(suppressWarnings(coloc.abf(d1, d2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)), silent = TRUE)
          fit <- if (inherits(res, "try-error")) list(pp4 = NA_real_, n_snps = nrow(g), status = "coloc_error") else list(pp4 = as.numeric(res$summary["PP.H4.abf"]), n_snps = nrow(g), status = "ok")
        }
      }
    }
    results[[length(results) + 1L]] <- data.table(signal_id = signal$signal_id, chb_stage = chb_stage,
                                                   pp4 = fit$pp4, n_snps = fit$n_snps, status = fit$status)
  }
}
unlink(gene_file)
result <- rbindlist(results)
stopifnot(nrow(result) == nrow(targets), all(is.na(result$pp4) | result$pp4 >= 0 & result$pp4 <= 1))
fwrite(result, output, sep = "\t")
message(sprintf("Wrote %d cross-stage CHB rows to %s", nrow(result), output))
