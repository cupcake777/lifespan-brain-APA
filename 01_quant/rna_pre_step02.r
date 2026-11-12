#!/usr/bin/env Rscript
# Align genotype, expression and covariate samples; write PGEN, tensorQTL BED and covariate table.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--expr_file", type = "character", help = "expression matrix (rows=genes, cols=samples)"),
  make_option("--anno", type = "character", help = "gene annotation bed (gencode)"),
  make_option("--geno", type = "character", help = "genotype VCF"),
  make_option("--cov", type = "character", default = NULL, help = "covariates (rows=covs, cols=samples), optional"),
  make_option("--outdir", type = "character", default = file.path(Sys.getenv("APA_ROOT"), "data/01_quant"))
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- samples shared by all inputs ---
geno_samples <- fread(cmd = sprintf("bcftools query -l %s", opt$geno), header = FALSE)$V1
expr_samples <- colnames(fread(opt$expr_file, nrows = 0))[-1]
shared_list <- list(geno_samples, expr_samples)
if (!is.null(opt$cov)) {
  cov_samples <- names(fread(opt$cov, nrows = 0, header = TRUE))[-1]
  shared_list <- list(geno_samples, expr_samples, cov_samples)
}
final_shared <- Reduce(intersect, shared_list)
final_samples <- geno_samples[geno_samples %in% final_shared]
if (length(final_samples) == 0) stop("No shared samples found across input files. Check sample IDs.")
cat(sprintf("Found %d shared samples. Using VCF sample order as reference.\n", length(final_samples)))
sample_list_path <- file.path(opt$outdir, "final_shared_samples.txt")
writeLines(final_samples, sample_list_path)

# --- genotype: PGEN for the shared samples ---
plink_prefix <- file.path(opt$outdir, "genotype.aligned")
system(sprintf("plink2 --vcf %s --keep %s --make-pgen --out %s --output-chr chrM",
               opt$geno, sample_list_path, plink_prefix))
if (!file.exists(paste0(plink_prefix, ".psam"))) stop("PLINK2 failed to create .psam file.")
plink_samples <- fread(paste0(plink_prefix, ".psam"), data.table = FALSE)$IID
if (!all(plink_samples == final_samples)) {
  warning("Final PLINK sample order does not match VCF order. This is unusual.")
  final_samples <- plink_samples
}

# --- phenotype: tensorQTL BED ---
expr_data <- fread(opt$expr_file, data.table = FALSE)
rownames(expr_data) <- expr_data[, 1]
expr_data <- expr_data[, -1]
if (!all(final_samples %in% colnames(expr_data))) {
  stop("Not all final samples were found in expression file headers.")
}
expr_aligned <- expr_data[, final_samples]

geneinfo <- fread(opt$anno, header = FALSE, stringsAsFactors = FALSE, select = 1:5)
colnames(geneinfo) <- c("chr", "start", "end", "strand", "gene_id")
geneinfo <- geneinfo[geneinfo$chr %in% paste0("chr", 1:22), ]
geneinfo$gene_id <- sub("\\..*$", "", geneinfo$gene_id)

common_genes <- intersect(rownames(expr_aligned), geneinfo$gene_id)
if (length(common_genes) == 0) stop("No common genes found between expression matrix and annotation file.")
expr_final <- expr_aligned[common_genes, ]
geneinfo_final <- geneinfo[match(common_genes, geneinfo$gene_id), ]
TSS <- ifelse(geneinfo_final$strand == "+", geneinfo_final$start, geneinfo_final$end)

bed_dt <- as.data.table(expr_final)
bed_dt[, `:=`(`#chr` = geneinfo_final$chr, start = TSS, end = TSS + 1, gene_id = geneinfo_final$gene_id)]
setcolorder(bed_dt, c("#chr", "start", "end", "gene_id"))
setkeyv(bed_dt, c("#chr", "start"))

phenotype_file <- file.path(opt$outdir, "phenotype.bed")
write.table(bed_dt, phenotype_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
system(sprintf("bgzip -f %s", phenotype_file))
system(sprintf("tabix -p bed %s.gz", phenotype_file))

# --- covariates (covs x samples) ---
if (!is.null(opt$cov)) {
  cov_data <- fread(opt$cov, data.table = FALSE)
  cov_id_name <- colnames(cov_data)[1]
  rownames(cov_data) <- cov_data[, 1]
  cov_data <- cov_data[, -1]
  if (!all(final_samples %in% colnames(cov_data))) stop("Not all final samples were found in covariates file.")
  cov_aligned_dt <- as.data.table(cov_data[, final_samples], keep.rownames = TRUE)
  setnames(cov_aligned_dt, "rn", cov_id_name)
  fwrite(cov_aligned_dt, file.path(opt$outdir, "covariates.aligned.txt"), sep = "\t", quote = FALSE)
}

cat("All tasks completed successfully.\n")
