#!/usr/bin/env Rscript
# APA function follow-up: c1 versus c2 tests, direction index, synaptic and MAGMA z models.
# Inputs from data/01_quant, data/02_traj, data/04_qtl and ref; tables go to data/06_disease.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(readxl)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
if (!nzchar(apa_root)) stop("Set APA_ROOT (source config/paths.local.sh)")
out <- file.path(apa_root, "data", "06_disease")
traj_dir <- file.path(apa_root, "data", "02_traj")
quant_dir <- file.path(apa_root, "data", "01_quant")
qtl_dir <- file.path(apa_root, "data", "04_qtl")
ref <- file.path(apa_root, "ref")

# --- helpers ---

descendants <- function(root_id) {
  ids <- root_id
  repeat {
    updated <- unique(c(ids, syngo_ont$id[syngo_ont$parent_id %in% ids]))
    if (length(updated) == length(ids)) return(ids)
    ids <- updated
  }
}

# --- loading data ---

genes <- fread(file.path(traj_dir, "05_gene_aggregation.tsv"))
tx2g <- fread(file.path(ref, "gencode.v40.tx2gene.tsv"),
              select = c("gene_id", "gene_name", "start", "end"))
tx2g[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
span <- tx2g[, .(gene_length = max(end) - min(start) + 1),
             by = .(gene_symbol = gene_name, gene_id_base)]
genes <- merge(genes, span, by = c("gene_symbol", "gene_id_base"))
cov <- fread(file.path(traj_dir, "00_cov.tsv"))
expr <- fread(file.path(quant_dir, "expr_after_combat.tsv"), showProgress = FALSE)
setnames(expr, 1L, "gene_id")
expr[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
samples <- intersect(names(expr), cov$sample)
weights <- cov[match(samples, sample), profile_weight]
expr[, mean_expression := as.numeric(as.matrix(.SD) %*% weights / sum(weights)),
     .SDcols = samples]
genes <- merge(genes, unique(expr[, .(gene_id_base, mean_expression)], by = "gene_id_base"),
               by = "gene_id_base")
genes <- unique(genes[is.finite(gene_length) & is.finite(mean_expression),
                      .(gene_symbol, gene_id_base, n_sites, n_c1, n_c2,
                        gene_state, membership, gene_length, mean_expression)],
                by = "gene_symbol")
stopifnot(nrow(genes) == 11322L, all(genes$n_c1 + genes$n_c2 == genes$n_sites))

gmt <- do.call(c, lapply(file.path(ref, c(
  "c5.go.bp.v2024.1.Hs.symbols.gmt",
  "c5.go.mf.v2024.1.Hs.symbols.gmt",
  "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")), function(path) {
  rows <- strsplit(readLines(path), "\t", fixed = TRUE)
  setNames(lapply(rows, function(x) unique(x[-c(1L, 2L)])),
           vapply(rows, `[`, character(1), 1L))
}))

mito <- read_excel(file.path(ref, "Human.MitoCarta3.0.xls"),
                   sheet = "A Human MitoCarta3.0")
syngo_dir <- tempfile("syngo_")
dir.create(syngo_dir)
unzip(file.path(ref, "syngo1.3_complete_data.zip"), exdir = syngo_dir)
syngo_genes <- read_excel(file.path(syngo_dir, "genes.xlsx"))
syngo_ann <- read_excel(file.path(syngo_dir, "annotations.xlsx"))
syngo_ont <- read_excel(file.path(syngo_dir, "ontologies.xlsx"))
rpl <- fread(file.path(ref, "HGNC_L_ribosomal_proteins.txt"))
rps <- fread(file.path(traj_dir, "HGNC_S_ribosomal_prots.tsv"))
proteasome <- fread(file.path(ref, "HGNC_Proteasome.tsv"))
stopifnot(all(rps$`Group ID` == 728L))
broad <- list(
  MitoCarta3 = unique(mito$Symbol),
  SynGO_all = unique(c(syngo_genes$hgnc_symbol, syngo_ann$hgnc_symbol)),
  SynGO_presynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098793")]),
  SynGO_postsynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098794")]),
  HGNC_RPL_RPS = unique(c(rpl$`Approved symbol`, rps$`Approved symbol`)),
  HGNC_Proteasome = unique(proteasome$`Approved symbol`)
)
sets <- c(gmt, broad)
sets <- lapply(sets, function(x) unique(x[!is.na(x) & nzchar(x)]))

# --- t2a: c1 versus c2 ---
enrichment <- fread(file.path(traj_dir, "08_enrich_two_levels.tsv"))
selected <- unique(enrichment[source_group == "GO/KEGG" &
                                fdr_by_gene_set < 0.05 &
                                gene_set %in% c("C1-spec", "C2-spec"), term])
test_sets <- c(gmt[intersect(names(gmt), selected)], broad)
direct <- genes[gene_state %in% c("C1-spec", "C2-spec")]
direct[, in_c1 := as.integer(gene_state == "C1-spec")]
tests <- list()
for (set_name in names(test_sets)) {
  members <- test_sets[[set_name]]
  in_set <- direct$gene_symbol %in% members
  tab <- table(factor(direct$in_c1, levels = 0:1), factor(in_set, levels = c(FALSE, TRUE)))
  fisher <- fisher.test(tab)
  fit <- glm(in_c1 ~ in_set + scale(log1p(n_sites)) +
               scale(log1p(gene_length)) + scale(mean_expression),
             data = direct, family = binomial())
  est <- coef(summary(fit))["in_setTRUE", ]
  tests[[set_name]] <- data.table(gene_set = set_name, source = if (set_name %in% names(broad)) "broad" else "GO/KEGG",
             n_c1 = sum(direct$in_c1), n_c2 = sum(!direct$in_c1),
             n_set_c1 = sum(direct$in_c1 == 1L & in_set),
             n_set_c2 = sum(direct$in_c1 == 0L & in_set),
             fisher_or = unname(fisher$estimate), fisher_ci_low = fisher$conf.int[1L],
             fisher_ci_high = fisher$conf.int[2L], fisher_p = fisher$p.value,
             adjusted_or = exp(est["Estimate"]),
             adjusted_ci_low = exp(est["Estimate"] - 1.96 * est["Std. Error"]),
             adjusted_ci_high = exp(est["Estimate"] + 1.96 * est["Std. Error"]),
             adjusted_p = est["Pr(>|z|)"])
}
t2a <- rbindlist(tests)
t2a[, `:=`(fisher_fdr = p.adjust(fisher_p, "BH"),
           adjusted_fdr = p.adjust(adjusted_p, "BH"))]
fwrite(t2a, file.path(out, "T2A_c1_vs_c2.tsv"), sep = "\t")

# --- t2c: residualized direction index ---
direct[, direction_index := (n_c1 - n_c2) / (n_c1 + n_c2)]
direct[, direction_residual := resid(lm(direction_index ~
  scale(log1p(n_sites)) + scale(log1p(gene_length)) + scale(mean_expression), data = direct))]
ranks <- setNames(direct$direction_residual, direct$gene_symbol)
ranks <- sort(ranks, decreasing = TRUE)
fg <- fgseaMultilevel(pathways = sets, stats = ranks, minSize = 10L,
                      maxSize = 1000L, eps = 0)
fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ",", character(1))]
fg[, source := fifelse(pathway %in% names(broad), "broad", "GO/KEGG")]
setcolorder(fg, c("pathway", "source", "pval", "padj", "NES", "size", "leadingEdge"))
fwrite(fg, file.path(out, "T2C_direction_index_fgsea.tsv"), sep = "\t")

# --- t1: gencode transcript features ---
gtf <- file.path(ref, "gencode.v40.annotation.gtf.gz")
gtf_cmd <- sprintf("gzip -cd %s | awk -F '\t' '$3==\"transcript\" || $3==\"exon\" || $3==\"CDS\" || $3==\"UTR\" || $3==\"stop_codon\"'",
                   shQuote(gtf))
annotation <- fread(cmd = gtf_cmd, header = FALSE, sep = "\t", showProgress = FALSE)
setnames(annotation, c("seqname", "source", "feature", "start", "end", "score",
                       "strand", "frame", "attributes"))
annotation[, `:=`(
  gene_id_base = sub("\\.[0-9]+$", "", sub('.*gene_id "([^"]+)".*', '\\1', attributes)),
  transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', attributes)
)]
annotation <- annotation[gene_id_base %in% genes$gene_id_base]
tx <- annotation[feature == "transcript", .(gene_id_base, transcript_id, seqname,
                                             strand, start, end, attributes)]
tx[, `:=`(tx_base = sub("\\.[0-9]+$", "", transcript_id),
          tx_end = fifelse(strand == "+", end, start),
          canonical = grepl('tag "Ensembl_canonical"', attributes, fixed = TRUE))]
cds <- annotation[feature == "CDS", .(cds_length = sum(end - start + 1L),
                                       cds_min = min(start), cds_max = max(end)),
                  by = transcript_id]
utr <- merge(annotation[feature == "UTR", .(transcript_id, start, end, strand)],
             cds, by = "transcript_id")
utr <- utr[(strand == "+" & start > cds_max) |
             (strand == "-" & end < cds_min),
           .(utr3_length = sum(end - start + 1L),
             utr3_start = min(start), utr3_end = max(end)), by = transcript_id]
tx <- merge(tx, cds, by = "transcript_id", all.x = TRUE)
tx <- merge(tx, utr, by = "transcript_id", all.x = TRUE)
ends <- tx[, .(n_ends = uniqueN(tx_end), n_transcripts = .N,
               canonical_utr3_length = if (any(canonical & !is.na(utr3_length)))
                 as.numeric(max(utr3_length[canonical], na.rm = TRUE)) else NA_real_,
               longest_utr3_length = if (all(is.na(utr3_length))) NA_real_
                 else as.numeric(max(utr3_length, na.rm = TRUE)),
               canonical_cds_length = if (any(canonical & !is.na(cds_length)))
                 as.numeric(max(cds_length[canonical], na.rm = TRUE)) else NA_real_),
           by = gene_id_base]
features <- merge(genes, ends, by = "gene_id_base", all.x = TRUE)
features[, polyasite_pas_count := NA_integer_]
setcolorder(features, c("gene_id_base", "gene_symbol", "gene_state", "membership",
                        "n_sites", "n_c1", "n_c2", "polyasite_pas_count", "n_ends",
                        "n_transcripts", "canonical_utr3_length", "longest_utr3_length",
                        "gene_length", "canonical_cds_length", "mean_expression"))
fwrite(features, file.path(out, "T1_gene_annotation_features.tsv"), sep = "\t")

# --- t1: synaptic models ---
syn_sets <- broad[c("SynGO_all", "SynGO_presynapse", "SynGO_postsynapse")]
models <- list(
  unadjusted = "in_class ~ in_set",
  n_sites = "in_class ~ in_set + scale(log1p(n_sites))",
  n_ends = "in_class ~ in_set + scale(log1p(n_ends))",
  gene_span = "in_class ~ in_set + scale(log1p(gene_length))",
  utr3_length = "in_class ~ in_set + scale(log1p(canonical_utr3_length))",
  expression = "in_class ~ in_set + scale(mean_expression)"
)
class_mask <- list(C1spec = features$gene_state == "C1-spec",
                   Shared = features$membership == "Shared",
                   C2spec = features$gene_state == "C2-spec")
t1_models <- rbindlist(lapply(names(class_mask), function(class_name) {
  rbindlist(lapply(names(syn_sets), function(set_name) {
    rbindlist(lapply(names(models), function(model_name) {
      model_data <- copy(features)
      model_data[, `:=`(in_class = as.integer(class_mask[[class_name]]),
                       in_set = as.integer(gene_symbol %in% syn_sets[[set_name]]))]
      fit <- glm(as.formula(models[[model_name]]), data = model_data, family = binomial())
      est <- coef(summary(fit))["in_set", ]
      data.table(gene_class = class_name, gene_set = set_name, model = model_name,
                 n = nobs(fit), n_class = sum(model.frame(fit)$in_class),
                 n_set = sum(model.frame(fit)$in_set),
                 odds_ratio = exp(est["Estimate"]),
                 ci_low = exp(est["Estimate"] - 1.96 * est["Std. Error"]),
                 ci_high = exp(est["Estimate"] + 1.96 * est["Std. Error"]),
                 p = est["Pr(>|z|)"])
    }))
  }))
}))
t1_models[, fdr := p.adjust(p, "BH"), by = .(gene_class, model)]
fwrite(t1_models, file.path(out, "T1_synaptic_models.tsv"), sep = "\t")

# --- t1: shared isoform pairs ---
stops <- annotation[feature == "stop_codon", .(stop_start = min(start),
                                                stop_end = max(end)), by = transcript_id]
terminal <- annotation[feature == "exon", .(terminal_start = if (strand[1L] == "+")
  start[which.max(end)] else min(start),
  terminal_end = if (strand[1L] == "+") max(end) else end[which.min(start)]),
  by = transcript_id]
tx <- merge(tx, stops, by = "transcript_id", all.x = TRUE)
tx <- merge(tx, terminal, by = "transcript_id", all.x = TRUE)
sites <- fread(file.path(traj_dir, "03_consensus_annotation.tsv"))
sites[, `:=`(transcript_id = sub("\\|.*", "", site_id),
             site_coord = sub("^.*\\|(chr[^|]+)\\|[+-]$", "\\1", site_id))]
sites[, `:=`(site_start = as.integer(sub("^.*:([0-9]+)-.*$", "\\1", site_coord)),
             site_end = as.integer(sub("^.*-([0-9]+)$", "\\1", site_coord)))]
sites <- merge(sites, tx[, .(transcript_id, gene_id_base, utr3_start, utr3_end,
                             stop_start, stop_end, terminal_start, terminal_end)],
               by = "transcript_id")
sites <- merge(sites, genes[, .(gene_id_base, membership)], by = "gene_id_base")
sites <- sites[membership == "Shared" & cluster_id %in% c(1L, 2L)]
setnames(sites, c("site_id", "transcript_id", "site_start", "site_end", "stop_start",
                  "stop_end", "terminal_start", "terminal_end"),
         paste0(c("site_id", "transcript_id", "site_start", "site_end", "stop_start",
                  "stop_end", "terminal_start", "terminal_end"), "_c1"))
one <- sites[cluster_id == 1L]
setnames(sites, names(sites), sub("_c1$", "_c2", names(sites)))
two <- sites[cluster_id == 2L]
pairs <- merge(one, two, by = "gene_id_base", suffixes = c("_c1", "_c2"),
               allow.cartesian = TRUE)
pairs <- pairs[site_end_c1 < site_start_c2 | site_end_c2 < site_start_c1]
pairs[, `:=`(both_coding = !is.na(stop_start_c1) & !is.na(stop_start_c2),
             same_stop_codon = fifelse(is.na(stop_start_c1) | is.na(stop_start_c2),
                                       NA, stop_start_c1 == stop_start_c2 &
                                         stop_end_c1 == stop_end_c2),
             different_terminal_exon = fifelse(is.na(terminal_start_c1) |
                                                is.na(terminal_start_c2), NA,
                                                terminal_start_c1 != terminal_start_c2 |
                                                  terminal_end_c1 != terminal_end_c2))]
pair_output <- pairs[, .(gene_id_base, gene_symbol = gene_id_c1,
                         transcript_id_c1, transcript_id_c2, site_id_c1, site_id_c2,
                         site_start_c1, site_end_c1, site_start_c2, site_end_c2,
                         utr3_start_c1, utr3_end_c1, utr3_start_c2, utr3_end_c2,
                         both_coding, same_stop_codon, different_terminal_exon)]
fwrite(pair_output, file.path(out, "T1_shared_isoform_pairs.tsv"), sep = "\t")

# --- t1: transcript expression matching ---
iso <- fread(file.path(quant_dir, "isoform_after_combat.tsv"),
             showProgress = FALSE)
setnames(iso, 1L, "transcript_id")
iso[, tx_base := sub("\\.[0-9]+$", "", transcript_id)]
iso <- merge(iso, unique(tx[!is.na(utr3_length), .(tx_base, gene_id_base)],
                         by = "tx_base"), by = "tx_base")
prenatal <- intersect(names(iso), cov[age < 0, sample])
postnatal <- intersect(names(iso), cov[age >= 0, sample])
stopifnot(length(prenatal) > 0L, length(postnatal) > 0L)
iso[, `:=`(pre_cpm = rowMeans(2^as.matrix(.SD)),
           post_cpm = NA_real_), .SDcols = prenatal]
iso[, post_cpm := rowMeans(2^as.matrix(.SD)), .SDcols = postnatal]
tx_expr <- iso[, .(pre_cpm = sum(pre_cpm), post_cpm = sum(post_cpm)),
               by = gene_id_base]
tx_expr[, min_stage_cpm := pmin(pre_cpm, post_cpm)]
match_data <- merge(features[, .(gene_id_base, gene_symbol, membership)], tx_expr,
                    by = "gene_id_base")
match_data[, synaptic := gene_symbol %in% syn_sets$SynGO_all]
match_data[, decile := pmin(10L, as.integer(ceiling(10 * frank(min_stage_cpm,
                                                               ties.method = "average") / .N)))]
matched <- rbindlist(lapply(1:10, function(bin) {
  syn <- match_data[decile == bin & synaptic == TRUE][order(min_stage_cpm)]
  control <- match_data[decile == bin & synaptic == FALSE][order(min_stage_cpm)]
  if (!nrow(syn) || !nrow(control)) return(NULL)
  chosen <- integer(nrow(syn))
  available <- seq_len(nrow(control))
  for (i in seq_len(nrow(syn))) {
    if (!length(available)) break
    pick <- available[which.min(abs(log1p(control$min_stage_cpm[available]) -
                                    log1p(syn$min_stage_cpm[i])))]
    chosen[i] <- pick
    available <- available[available != pick]
  }
  syn <- syn[chosen > 0L]
  control <- control[chosen[chosen > 0L]]
  data.table(decile = bin, syn_gene = syn$gene_symbol,
             control_gene = control$gene_symbol,
             syn_shared = syn$membership == "Shared",
             control_shared = control$membership == "Shared",
             syn_pre_cpm = syn$pre_cpm, control_pre_cpm = control$pre_cpm,
             syn_post_cpm = syn$post_cpm, control_post_cpm = control$post_cpm,
             syn_min_cpm = syn$min_stage_cpm,
             control_min_cpm = control$min_stage_cpm)
}))
fwrite(matched, file.path(out, "T1_matched_synaptic.tsv"), sep = "\t")

# --- t2b: dapars2-quantifiable nondirectional genes ---
quant_ids <- fread(cmd = sprintf("cut -f1 %s", shQuote(file.path(qtl_dir,
  "apa_after_combat.tsv"))), header = FALSE)[[1L]][-1L]
quant <- data.table(tx_base = unique(sub("\\.[0-9]+$", "", sub("\\|.*", "", quant_ids))))
full_tx <- fread(file.path(ref, "gencode.v40.tx2gene.tsv"),
                 select = c("tx_name", "gene_id", "gene_name", "start", "end"))
full_tx[, `:=`(tx_base = sub("\\.[0-9]+$", "", tx_name),
               gene_id_base = sub("\\.[0-9]+$", "", gene_id))]
quant <- merge(quant, unique(full_tx[, .(tx_base, gene_id_base, gene_symbol = gene_name)],
                             by = "tx_base"), by = "tx_base")
quant <- unique(quant[, .(gene_symbol, gene_id_base)], by = "gene_symbol")
non_directional <- quant[!gene_symbol %in% fread(file.path(traj_dir,
  "05_gene_aggregation.tsv"), select = "gene_symbol")[[1L]]]
full_span <- full_tx[, .(gene_length = max(end) - min(start) + 1),
                     by = .(gene_symbol = gene_name, gene_id_base)]
non_directional <- merge(non_directional, full_span,
                         by = c("gene_symbol", "gene_id_base"), all.x = TRUE)
non_directional <- merge(non_directional,
                         unique(expr[, .(gene_id_base, mean_expression)],
                                by = "gene_id_base"), by = "gene_id_base", all.x = TRUE)
non_directional[, has_expression := is.finite(mean_expression)]
fwrite(non_directional, file.path(out, "T2B_nondynamic_reference.tsv"), sep = "\t")

reference <- non_directional[has_expression == TRUE & is.finite(gene_length)]
t2b <- rbindlist(lapply(names(class_mask), function(class_name) {
  target <- features[class_mask[[class_name]],
                     .(gene_symbol, gene_length, mean_expression)]
  comparison <- rbindlist(list(target[, in_class := 1L],
                               reference[, .(gene_symbol, gene_length,
                                             mean_expression, in_class = 0L)]))
  rbindlist(lapply(names(test_sets), function(set_name) {
    comparison[, in_set := as.integer(gene_symbol %in% test_sets[[set_name]])]
    tab <- table(factor(comparison$in_class, levels = 0:1),
                 factor(comparison$in_set, levels = 0:1))
    f <- fisher.test(tab)
    m <- glm(in_class ~ in_set + scale(log1p(gene_length)) +
               scale(mean_expression), data = comparison, family = binomial())
    e <- coef(summary(m))["in_set", ]
    data.table(gene_class = class_name, gene_set = set_name,
               n_class = nrow(target), n_reference = nrow(reference),
               n_set_class = sum(comparison$in_class == 1L & comparison$in_set == 1L),
               n_set_reference = sum(comparison$in_class == 0L & comparison$in_set == 1L),
               fisher_or = unname(f$estimate), fisher_p = f$p.value,
               adjusted_or = exp(e["Estimate"]),
               adjusted_ci_low = exp(e["Estimate"] - 1.96 * e["Std. Error"]),
               adjusted_ci_high = exp(e["Estimate"] + 1.96 * e["Std. Error"]),
               adjusted_p = e["Pr(>|z|)"])
  }))
}))
t2b[, `:=`(fisher_fdr = p.adjust(fisher_p, "BH"),
            adjusted_fdr = p.adjust(adjusted_p, "BH")), by = gene_class]
fwrite(t2b, file.path(out, "T2B_nondynamic_tests.tsv"), sep = "\t")

# --- t5: magma gene z models ---
gene_loc <- fread(file.path(ref, "NCBI38.gene.loc"), header = FALSE,
                  col.names = c("gene_id", "chr", "start", "stop", "strand", "gene_symbol"))
gene_loc[, gene_id := as.character(gene_id)]
gwas_ids <- unique(fread(file.path(traj_dir, "magma_traj_summary.tsv"),
                         select = "gwas_id")[[1L]])
z_formulas <- list(
  base = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length)",
  sites = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites))",
  ends_without_pas = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites)) + scale(log1p(n_ends))",
  expression_without_pas = "ZSTAT ~ in_class + log(NSNPS) + log(NPARAM) + log(gene_length) + scale(log1p(n_sites)) + scale(log1p(n_ends)) + scale(mean_expression)"
)
t5 <- rbindlist(lapply(gwas_ids, function(gwas_id) {
  path <- file.path(out, "enrichment", gwas_id, "baseline.genes.out")
  raw <- fread(path)
  setnames(raw, "GENE", "gene_id")
  raw[, gene_id := as.character(gene_id)]
  raw <- merge(raw, unique(gene_loc[, .(gene_id, gene_symbol)], by = "gene_id"),
               by = "gene_id")
  raw <- merge(raw, features, by = "gene_symbol")
  raw <- unique(raw[NSNPS > 0 & NPARAM > 0 & gene_length > 0 & is.finite(ZSTAT)],
                by = "gene_symbol")
  rbindlist(lapply(names(class_mask), function(class_name) {
    model_data <- copy(raw)
    model_data[, in_class := as.integer(if (class_name == "Shared")
      membership == "Shared" else gene_state == sub("spec$", "-spec", class_name))]
    rbindlist(lapply(names(z_formulas), function(model_name) {
      fit <- lm(as.formula(z_formulas[[model_name]]), data = model_data)
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
t5[, fdr := p.adjust(p, "BH"), by = .(gene_class, model)]
fwrite(t5, file.path(out, "T5_magma_gene_z_models.tsv"), sep = "\t")
