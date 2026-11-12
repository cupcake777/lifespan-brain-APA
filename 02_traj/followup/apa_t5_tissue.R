#!/usr/bin/env Rscript
# Tissue specificity of APA-class genes (HPA elevated tissues, GTEx tau) and MAGMA gene-level models.
# Reads 02_traj tables and HPA from ref/; writes T5 tables to data/02_traj/.
suppressPackageStartupMessages(library(data.table))

# --- config ---
out <- file.path(Sys.getenv("APA_ROOT"), "data", "02_traj")
hpa_zip <- file.path(Sys.getenv("APA_ROOT"), "ref", "hpa", "proteinatlas.tsv.zip")
gtex_gct <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment/GTEx_v8_gene_median_tpm.gct.gz")
magma_dir <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/10_magma")

# --- loading data ---
features <- fread(file.path(out, "T1_gene_annotation_features.tsv"))
reference <- fread(file.path(out, "T2B_nondynamic_reference.tsv"))
hpa <- fread(cmd = sprintf("unzip -p %s proteinatlas.tsv", shQuote(hpa_zip)),
             select = c("Ensembl", "RNA tissue specificity", "RNA tissue distribution",
                        "RNA tissue specific nTPM"))
setnames(hpa, c("Ensembl", "RNA tissue specificity", "RNA tissue distribution",
                "RNA tissue specific nTPM"),
         c("gene_id_base", "specificity", "distribution", "specific_ntpm"))
hpa <- unique(hpa, by = "gene_id_base")
gtex <- fread(cmd = sprintf("gzip -cd %s", shQuote(gtex_gct)), skip = 2L)
setnames(gtex, 1:2, c("gene_id", "gene_symbol_gtex"))
gtex[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
gtex <- unique(gtex, by = "gene_id_base")
tissues <- setdiff(names(gtex), c("gene_id", "gene_symbol_gtex", "gene_id_base"))
brain_tissues <- grep("^Brain - ", tissues, value = TRUE)
stopifnot(length(tissues) == 54L, length(brain_tissues) == 13L,
          !anyDuplicated(features$gene_id_base), !anyDuplicated(reference$gene_id_base))

# --- hpa elevated tissues ---
elevated <- hpa[specificity %in% c("Tissue enriched", "Group enriched", "Tissue enhanced") &
                  !is.na(specific_ntpm) & nzchar(specific_ntpm)]
listed <- strsplit(elevated$specific_ntpm, ";", fixed = TRUE)
elevated_tissue <- data.table(gene_id_base = rep(elevated$gene_id_base, lengths(listed)),
                             tissue = trimws(sub(":.*$", "", unlist(listed))))
elevated_tissue <- unique(elevated_tissue[nzchar(tissue)])
stopifnot(all(elevated_tissue$gene_id_base %in% hpa$gene_id_base))

direct <- merge(features[, .(gene_id_base, gene_state, membership, n_sites,
                             gene_length, mean_expression)],
                hpa[, .(gene_id_base, distribution)], by = "gene_id_base")
direct[, gene_class := fifelse(gene_state == "C1-spec", "C1spec",
                              fifelse(gene_state == "C2-spec", "C2spec", "Shared"))]
control <- merge(reference[has_expression == TRUE & is.finite(gene_length),
                           .(gene_id_base, gene_length, mean_expression)],
                 hpa[, .(gene_id_base, distribution)], by = "gene_id_base")
control[, `:=`(gene_class = "NonDirectional", n_sites = NA_integer_)]
universe <- rbindlist(list(direct, control), fill = TRUE)
stopifnot(!anyDuplicated(universe$gene_id_base), all(is.finite(universe$gene_length)),
          all(is.finite(universe$mean_expression)))

comparisons <- list(
  list(background = "A_C1_vs_C2", gene_class = "C1spec", reference_class = "C2spec"),
  list(background = "A_C1_vs_C2", gene_class = "C2spec", reference_class = "C1spec"),
  list(background = "B_nondirectional", gene_class = "C1spec", reference_class = "NonDirectional"),
  list(background = "B_nondirectional", gene_class = "Shared", reference_class = "NonDirectional"),
  list(background = "B_nondirectional", gene_class = "C2spec", reference_class = "NonDirectional")
)

tissue_results <- rbindlist(lapply(comparisons, function(cmp) {
  comparison <- universe[gene_class %in% c(cmp$gene_class, cmp$reference_class)]
  comparison[, in_class := as.integer(gene_class == cmp$gene_class)]
  rbindlist(lapply(sort(unique(elevated_tissue$tissue)), function(tissue) {
    this_tissue <- tissue
    comparison[, in_tissue := as.integer(gene_id_base %in%
                                         elevated_tissue[tissue == this_tissue, gene_id_base])]
    tab <- table(factor(comparison$in_class, levels = 0:1),
                 factor(comparison$in_tissue, levels = 0:1))
    f <- fisher.test(tab)
    formula <- if (cmp$background == "A_C1_vs_C2")
      in_class ~ in_tissue + scale(log1p(n_sites)) +
        scale(log1p(gene_length)) + scale(mean_expression)
    else in_class ~ in_tissue + scale(log1p(gene_length)) + scale(mean_expression)
    fit <- suppressWarnings(glm(formula, data = comparison, family = binomial()))
    est <- coef(summary(fit))["in_tissue", ]
    data.table(background = cmp$background, gene_class = cmp$gene_class,
               reference_class = cmp$reference_class, tissue = tissue,
               is_brain = tissue %in% c("brain", "cerebral cortex"),
               n_class = sum(comparison$in_class), n_reference = sum(1L - comparison$in_class),
               n_elevated_class = tab[2L, 2L], n_elevated_reference = tab[1L, 2L],
               n_detected_all_class = sum(comparison$in_class == 1L &
                                          comparison$distribution == "Detected in all"),
               n_detected_all_reference = sum(comparison$in_class == 0L &
                                              comparison$distribution == "Detected in all"),
               fisher_or = unname(f$estimate), fisher_ci_low = f$conf.int[1L],
               fisher_ci_high = f$conf.int[2L], fisher_p = f$p.value,
               adjusted_or = exp(est["Estimate"]),
               adjusted_ci_low = exp(est["Estimate"] - 1.96 * est["Std. Error"]),
               adjusted_ci_high = exp(est["Estimate"] + 1.96 * est["Std. Error"]),
               adjusted_p = est["Pr(>|z|)"])
  }))
}))
tissue_results[, `:=`(fisher_fdr = p.adjust(fisher_p, "BH"),
                      adjusted_fdr = p.adjust(adjusted_p, "BH")),
               by = .(background, gene_class)]
fwrite(tissue_results, file.path(out, "T5_tissue_specificity.tsv"), sep = "\t")

# --- gtex tissue specificity ---
matrix_tpm <- as.matrix(gtex[, ..tissues])
storage.mode(matrix_tpm) <- "double"
max_tpm <- apply(matrix_tpm, 1L, max)
tau <- rowSums(1 - sweep(matrix_tpm, 1L, pmax(max_tpm, .Machine$double.eps), "/")) /
  (length(tissues) - 1L)
tau[max_tpm <= 0] <- NA_real_
brain_max <- apply(matrix_tpm[, brain_tissues, drop = FALSE], 1L, max)
other_max <- apply(matrix_tpm[, setdiff(tissues, brain_tissues), drop = FALSE], 1L, max)
gtex[, `:=`(tau = tau, peak_tissue = tissues[max.col(matrix_tpm, ties.method = "first")],
           brain_peak = brain_max > other_max & max_tpm > 0)]
gtex_gene <- merge(universe[, .(gene_id_base, gene_class)],
                   gtex[, .(gene_id_base, tau, peak_tissue, brain_peak)],
                   by = "gene_id_base")
gtex_gene <- gtex_gene[is.finite(tau)]

tau_comparisons <- c(comparisons, list(
  list(background = "C1_vs_Shared", gene_class = "C1spec", reference_class = "Shared"),
  list(background = "Shared_vs_C2", gene_class = "Shared", reference_class = "C2spec")))
tau_results <- rbindlist(lapply(tau_comparisons, function(cmp) {
  target <- gtex_gene[gene_class == cmp$gene_class]
  baseline <- gtex_gene[gene_class == cmp$reference_class]
  f <- fisher.test(matrix(c(sum(target$brain_peak), sum(!target$brain_peak),
                            sum(baseline$brain_peak), sum(!baseline$brain_peak)), nrow = 2L,
                          byrow = TRUE))
  data.table(background = cmp$background, gene_class = cmp$gene_class,
             reference_class = cmp$reference_class,
             n_class = nrow(target), n_reference = nrow(baseline),
             tau_median_class = median(target$tau),
             tau_q1_class = unname(quantile(target$tau, 0.25)),
             tau_q3_class = unname(quantile(target$tau, 0.75)),
             tau_median_reference = median(baseline$tau),
             tau_q1_reference = unname(quantile(baseline$tau, 0.25)),
             tau_q3_reference = unname(quantile(baseline$tau, 0.75)),
             tau_wilcox_p = wilcox.test(target$tau, baseline$tau, exact = FALSE)$p.value,
             n_brain_peak_class = sum(target$brain_peak),
             n_brain_peak_reference = sum(baseline$brain_peak),
             brain_peak_fraction_class = mean(target$brain_peak),
             brain_peak_fraction_reference = mean(baseline$brain_peak),
             brain_peak_fisher_or = unname(f$estimate),
             brain_peak_ci_low = f$conf.int[1L],
             brain_peak_ci_high = f$conf.int[2L], brain_peak_fisher_p = f$p.value,
             brain_tissues = paste(brain_tissues, collapse = ";"))
}))
fwrite(tau_results, file.path(out, "T5_gtex_tau_by_class.tsv"), sep = "\t")

# --- magma gene z models ---
stopifnot(!anyNA(features$polyasite_pas_count),
          all(features$polyasite_pas_count >= 0L))
gene_loc <- fread(file.path(magma_dir, "ref/NCBI38.gene.loc"), header = FALSE,
                  col.names = c("gene_id", "chr", "start", "stop", "strand", "gene_symbol"))
gene_loc[, gene_id := as.character(gene_id)]
gwas_ids <- unique(fread(file.path(out, "magma_traj_summary.tsv"),
                         select = "gwas_id")[[1L]])
models <- list(
  base = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length)",
  sites = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites))",
  pas = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites)) + scale(log1p(polyasite_pas_count))",
  ends = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites)) + scale(log1p(polyasite_pas_count)) + scale(log1p(n_ends))",
  expression = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites)) + scale(log1p(polyasite_pas_count)) + scale(log1p(n_ends)) + scale(mean_expression)"
)
magma <- rbindlist(lapply(gwas_ids, function(gwas_id) {
  raw <- fread(file.path(magma_dir, "enrichment", gwas_id, "baseline.genes.out"))
  setnames(raw, "GENE", "gene_id")
  raw[, gene_id := as.character(gene_id)]
  raw <- merge(raw, unique(gene_loc[, .(gene_id, gene_symbol)], by = "gene_id"),
               by = "gene_id")
  raw <- merge(raw, features, by = "gene_symbol")
  raw <- unique(raw[NSNPS > 0 & NPARAM > 0 & gene_length > 0 & is.finite(ZSTAT)],
                by = "gene_symbol")
  rbindlist(lapply(c("C1spec", "Shared", "C2spec"), function(class_name) {
    model_data <- copy(raw)
    model_data[, in_class := as.integer(if (class_name == "Shared")
      membership == "Shared" else gene_state == sub("spec$", "-spec", class_name))]
    rbindlist(lapply(names(models), function(model_name) {
      fit <- lm(as.formula(models[[model_name]]), data = model_data)
      est <- coef(summary(fit))["in_class", ]
      data.table(gwas_id = gwas_id, gene_class = class_name,
                 model = model_name, n = nobs(fit),
                 n_class = sum(model.frame(fit)$in_class),
                 beta = est["Estimate"], se = est["Std. Error"],
                 ci_low = est["Estimate"] - 1.96 * est["Std. Error"],
                 ci_high = est["Estimate"] + 1.96 * est["Std. Error"],
                 p = est["Pr(>|t|)"])
    }))
  }))
}))
magma[, fdr := p.adjust(p, "BH"), by = .(gene_class, model)]
fwrite(magma, file.path(out, "T5_magma_gene_z_models.tsv"), sep = "\t")