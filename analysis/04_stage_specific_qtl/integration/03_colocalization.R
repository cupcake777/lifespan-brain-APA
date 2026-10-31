# Colocalization (coloc.abf) of each significant APA locus with the eQTL of its gene, per life stage.
# Reads nominal_*_<stage>.tsv.gz from 02-1_extract.py; writes coloc_* to data/04_qtl.
suppressMessages({
  library(data.table)
  library(coloc)
  library(parallel)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
GTF_FILE <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")
STAGES <- c("prenatal", "postnatal")
MIN_SHARED_SNPS <- 50
PP_CALL <- 0.8
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
P1 <- 1e-4
P2 <- 1e-4
P12 <- 1e-5
CLASS_LEVELS <- c("H4 shared causal variant", "H3 distinct causal variants", "H2 APA signal only",
                  "H1 eQTL signal only", "H0 no signal", "undetermined")
setDTthreads(N_CORES)
NOMINAL_COLS <- c("phenotype_id", "variant_id", "af", "ma_count", "pval_nominal", "slope", "slope_se")

# --- helpers ---
estimate_n <- function(dt) {
  s <- dt[!is.na(ma_count) & af > 0 & af < 1]
  if (nrow(s) == 0) return(NA_integer_)
  s <- s[sample(.N, min(.N, 20000L))]
  maf <- pmin(s$af, 1 - s$af)
  as.integer(round(median(s$ma_count / (2 * maf))))
}

# --- gene annotation (ENSG <-> symbol) ---
gene_meta <- fread(cmd = paste0(
  "zcat ", GTF_FILE, " | awk -F'\\t' '$3==\"gene\" {",
  "match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10);",
  "match($9,/gene_name \"[^\"]+\"/); s=substr($9,RSTART+11,RLENGTH-12);",
  "print $1\"\\t\"$4\"\\t\"$5\"\\t\"$7\"\\t\"g\"\\t\"s}'"),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("chr", "gene_start", "gene_end", "strand", "ensg", "sym"))
gene_meta[, ensg0 := sub("\\..*", "", ensg)]
gene_meta <- unique(gene_meta, by = "ensg0")
ensg0_to_sym <- setNames(gene_meta$sym, gene_meta$ensg0)
sym_map <- unique(gene_meta[!is.na(sym) & sym != ""], by = "sym")
sym_to_ensg0 <- setNames(sym_map$ensg0, sym_map$sym)

apa_gene_to_ensg0 <- function(gene_raw) {
  out <- ifelse(grepl("^ENSG", gene_raw), sub("\\..*", "", gene_raw), NA_character_)
  need <- is.na(out)
  out[need] <- unname(sym_to_ensg0[gene_raw[need]])
  out
}

parse_apa_pheno <- function(dt) {
  parts <- tstrsplit(dt$phenotype_id, "|", fixed = TRUE, fill = NA_character_)
  dt[, `:=`(transcript_id = parts[[1]], gene_raw = parts[[2]],
            apa_interval = parts[[3]], apa_strand = parts[[4]])]
  dt[, ensg0 := apa_gene_to_ensg0(gene_raw)]
  dt[, gene_name := unname(ensg0_to_sym[ensg0])]
  dt[is.na(gene_name) | gene_name == "", gene_name := gene_raw]
  dt[]
}

# --- coloc per gene: all eQTL-phenotype x APA-phenotype pairs sharing >= MIN_SHARED_SNPS variants ---
coloc_one_gene <- function(g, dt_e, dt_a, n_e, n_a) {
  sub_e <- dt_e[.(g), nomatch = 0L]
  sub_a <- dt_a[.(g), nomatch = 0L]
  if (nrow(sub_e) == 0L || nrow(sub_a) == 0L) return(NULL)
  out <- list()
  for (e_tr in unique(sub_e$phenotype_id)) {
    d_e0 <- sub_e[phenotype_id == e_tr][!duplicated(variant_id)]
    for (a_tr in unique(sub_a$phenotype_id)) {
      d_a0 <- sub_a[phenotype_id == a_tr][!duplicated(variant_id)]
      shared <- intersect(d_e0$variant_id, d_a0$variant_id)
      if (length(shared) < MIN_SHARED_SNPS) next
      d_e <- d_e0[match(shared, variant_id)]
      d_a <- d_a0[match(shared, variant_id)]
      ds_e <- list(beta = d_e$slope, varbeta = d_e$slope_se^2, snp = shared, type = "quant", sdY = 1, N = n_e)
      ds_a <- list(beta = d_a$slope, varbeta = d_a$slope_se^2, snp = shared, type = "quant", sdY = 1, N = n_a)
      res <- NULL
      capture.output(res <- tryCatch(coloc.abf(dataset1 = ds_e, dataset2 = ds_a, p1 = P1, p2 = P2, p12 = P12),
                                     error = function(e) NULL))
      if (is.null(res) || is.null(res$summary)) next
      rr <- as.data.table(res$results)
      best <- rr[which.max(SNP.PP.H4)]
      i_h4 <- match(best$snp, shared)
      i_le <- which.min(d_e$pval_nominal)
      i_la <- which.min(d_a$pval_nominal)
      out[[length(out) + 1L]] <- data.table(
        ensg0 = g, eqtl_pheno = e_tr, apa_pheno = a_tr,
        nsnps = as.integer(res$summary[["nsnps"]]),
        PP_H0 = res$summary[["PP.H0.abf"]], PP_H1 = res$summary[["PP.H1.abf"]],
        PP_H2 = res$summary[["PP.H2.abf"]], PP_H3 = res$summary[["PP.H3.abf"]],
        PP_H4 = res$summary[["PP.H4.abf"]],
        h4_snp = best$snp, h4_snp_pp = best$SNP.PP.H4,
        h4_slope_eqtl = d_e$slope[i_h4], h4_slope_apa = d_a$slope[i_h4],
        h4_p_eqtl = d_e$pval_nominal[i_h4], h4_p_apa = d_a$pval_nominal[i_h4],
        eqtl_lead_snp = d_e$variant_id[i_le], eqtl_lead_p = d_e$pval_nominal[i_le],
        eqtl_lead_slope = d_e$slope[i_le],
        apa_lead_snp = d_a$variant_id[i_la], apa_lead_p = d_a$pval_nominal[i_la],
        apa_lead_slope = d_a$slope[i_la]
      )
    }
  }
  if (length(out) == 0L) return(NULL)
  rbindlist(out)
}

# --- per stage ---
run_stage_coloc <- function(stage) {
  sig_apa <- unique(fread(file.path(D, paste0(stage, ".sig_QTL.txt.gz")), select = "phenotype_id"))
  sig_apa <- parse_apa_pheno(sig_apa)
  sig_apa <- sig_apa[!is.na(ensg0)]
  sig_eqtl <- unique(fread(file.path(D, paste0("eqtl.", stage, ".sig_QTL.txt.gz")), select = "phenotype_id"))
  sig_eqtl[, ensg0 := sub("\\..*", "", phenotype_id)]

  dt_e <- fread(file.path(D, paste0("nominal_eqtl_", stage, ".tsv.gz")), select = NOMINAL_COLS)
  dt_e[, ensg0 := sub("\\..*", "", phenotype_id)]
  dt_e <- dt_e[is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  setkey(dt_e, ensg0)
  dt_a <- fread(file.path(D, paste0("nominal_apa_", stage, ".tsv.gz")), select = NOMINAL_COLS)
  dt_a <- merge(dt_a, sig_apa[, .(phenotype_id, ensg0)], by = "phenotype_id", all.x = FALSE)
  dt_a <- dt_a[is.finite(slope) & is.finite(slope_se) & slope_se > 0]
  setkey(dt_a, ensg0)

  n_e <- estimate_n(dt_e)
  n_a <- estimate_n(dt_a)
  genes <- intersect(unique(dt_a$ensg0), unique(dt_e$ensg0))
  res <- rbindlist(mclapply(genes, coloc_one_gene, dt_e = dt_e, dt_a = dt_a, n_e = n_e, n_a = n_a,
                            mc.cores = N_CORES, mc.preschedule = TRUE), use.names = TRUE, fill = TRUE)
  if (nrow(res) == 0) stop("no coloc test succeeded for stage ", stage)

  res <- merge(res, unique(sig_apa[, .(apa_pheno = phenotype_id, gene_name, transcript_id, apa_interval, apa_strand)]),
               by = "apa_pheno", all.x = TRUE)
  res[, stage := stage]
  res[, is_egene := ensg0 %in% sig_eqtl$ensg0]
  res[, coloc_class := fcase(
    PP_H4 >= PP_CALL, CLASS_LEVELS[1],
    PP_H3 >= PP_CALL, CLASS_LEVELS[2],
    PP_H2 >= PP_CALL, CLASS_LEVELS[3],
    PP_H1 >= PP_CALL, CLASS_LEVELS[4],
    PP_H0 >= PP_CALL, CLASS_LEVELS[5],
    default = CLASS_LEVELS[6])]
  res[, direction := fifelse(h4_slope_eqtl * h4_slope_apa > 0, "concordant", "discordant")]
  setcolorder(res, c("stage", "ensg0", "gene_name", "eqtl_pheno", "apa_pheno", "transcript_id", "nsnps",
                     "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "coloc_class", "direction"))
  setorder(res, -PP_H4)
  fwrite(res, file.path(D, sprintf("coloc_%s_apa_eqtl.tsv", stage)), sep = "\t", quote = FALSE)
  res
}

stage_runs <- lapply(STAGES, run_stage_coloc)
names(stage_runs) <- STAGES
coloc_all <- rbindlist(stage_runs, use.names = TRUE, fill = TRUE)
coloc_all[, stage := factor(stage, levels = STAGES)]
coloc_all[, coloc_class := factor(coloc_class, levels = CLASS_LEVELS)]
fwrite(coloc_all, file.path(D, "coloc_all_stages.tsv"), sep = "\t", quote = FALSE)
