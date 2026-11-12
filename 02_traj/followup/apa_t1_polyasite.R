#!/usr/bin/env Rscript
# PolyA-site counts per gene and SynGO enrichment models for the C1, Shared and C2 gene classes.
# Reads and updates T1_gene_annotation_features.tsv in data/02_traj/.
stopifnot(nzchar(Sys.getenv(c("APA_ROOT", "APA_HPC_HOME"))))
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

# --- helpers ---
count_pas <- function(intervals, column) {
  x <- unique(intervals[, .(seqname, start, end, strand, gene_id_base)])
  setkey(x, seqname, start, end)
  hits <- foverlaps(pas, x, nomatch = 0L)[strand == i.strand,
                                           unique(.SD), .SDcols = c("gene_id_base", "pas_id")]
  setnames(hits[, .N, by = gene_id_base], "N", column)
}
descendants <- function(root_id) {
  ids <- root_id
  repeat {
    updated <- unique(c(ids, syngo_ont$id[syngo_ont$parent_id %in% ids]))
    if (length(updated) == length(ids)) return(ids)
    ids <- updated
  }
}

# --- config ---
out <- file.path(Sys.getenv("APA_ROOT"), "data", "02_traj")
gtf <- file.path(Sys.getenv("APA_ROOT"), "ref", "gencode.v40.annotation.gtf.gz")
pas_bed <- file.path(Sys.getenv("APA_ROOT"), "ref", "polyasite2", "atlas.clusters.2.0.GRCh38.96.bed.gz")
syngo_zip <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment/syngo1.3_complete_data.zip")

# --- loading data ---
features <- fread(file.path(out, "T1_gene_annotation_features.tsv"))
pas <- fread(cmd = paste("gzip -cd", shQuote(pas_bed)), header = FALSE,
             select = c(1L, 4L, 6L), col.names = c("seqname", "pas_id", "strand"))
pas[, `:=`(seqname = paste0("chr", seqname),
           start = as.integer(sub("^.*:([0-9]+):[+-]$", "\\1", pas_id)))]
pas[, end := start]
stopifnot(!anyNA(pas$start), uniqueN(pas$pas_id) == nrow(pas))

gtf_cmd <- sprintf("gzip -cd %s | awk -F '\t' '$3==\"exon\" || $3==\"CDS\" || $3==\"UTR\"'",
                   shQuote(gtf))
annotation <- fread(cmd = gtf_cmd, header = FALSE, sep = "\t", showProgress = FALSE)
setnames(annotation, c("seqname", "source", "feature", "start", "end", "score",
                       "strand", "frame", "attributes"))
annotation[, `:=`(
  gene_id_base = sub("\\.[0-9]+$", "", sub('.*gene_id "([^"]+)".*', '\\1', attributes)),
  transcript_id = sub('.*transcript_id "([^"]+)".*', '\\1', attributes)
)]
annotation <- annotation[gene_id_base %in% features$gene_id_base]
cds <- annotation[feature == "CDS", .(cds_min = min(start), cds_max = max(end)),
                  by = transcript_id]
utr3 <- merge(annotation[feature == "UTR",
                         .(seqname, start, end, strand, gene_id_base, transcript_id)],
              cds, by = "transcript_id")
utr3 <- utr3[(strand == "+" & start > cds_max) |
               (strand == "-" & end < cds_min)]
exons <- annotation[feature == "exon",
                    .(seqname, start, end, strand, gene_id_base, transcript_id)]
terminal <- exons[, .SD[if (strand[1L] == "+") which.max(end) else which.min(start)],
                  by = transcript_id]

# --- t1: polyasite counts ---
utr_counts <- count_pas(utr3, "polyasite_utr3_count")
terminal_counts <- count_pas(terminal, "polyasite_terminal_exon_count")
union_intervals <- unique(rbindlist(list(utr3[, .(seqname, start, end, strand, gene_id_base)],
                                         terminal[, .(seqname, start, end, strand, gene_id_base)])))
pas_counts <- count_pas(union_intervals, "polyasite_pas_count")
features[, intersect(names(features), c("polyasite_pas_count", "polyasite_utr3_count",
                                       "polyasite_terminal_exon_count")) := NULL]
features <- merge(features, pas_counts, by = "gene_id_base", all.x = TRUE)
features <- merge(features, utr_counts, by = "gene_id_base", all.x = TRUE)
features <- merge(features, terminal_counts, by = "gene_id_base", all.x = TRUE)
for (column in c("polyasite_pas_count", "polyasite_utr3_count",
                 "polyasite_terminal_exon_count")) set(features, which(is.na(features[[column]])),
                                                     column, 0L)
fwrite(features, file.path(out, "T1_gene_annotation_features.tsv"), sep = "\t")

# --- t1: synaptic models ---
syngo_dir <- tempfile("syngo_")
dir.create(syngo_dir)
unzip(syngo_zip, exdir = syngo_dir)
syngo_genes <- read_excel(file.path(syngo_dir, "genes.xlsx"))
syngo_ann <- read_excel(file.path(syngo_dir, "annotations.xlsx"))
syngo_ont <- read_excel(file.path(syngo_dir, "ontologies.xlsx"))
syn_sets <- list(
  SynGO_all = unique(c(syngo_genes$hgnc_symbol, syngo_ann$hgnc_symbol)),
  SynGO_presynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098793")]),
  SynGO_postsynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098794")])
)
models <- list(
  unadjusted = in_class ~ in_set,
  n_sites = in_class ~ in_set + scale(log1p(n_sites)),
  polyasite_pas = in_class ~ in_set + scale(log1p(polyasite_pas_count)),
  n_ends = in_class ~ in_set + scale(log1p(n_ends)),
  gene_span = in_class ~ in_set + scale(log1p(gene_length)),
  utr3_length = in_class ~ in_set + scale(log1p(canonical_utr3_length)),
  expression = in_class ~ in_set + scale(mean_expression),
  full = in_class ~ in_set + scale(log1p(n_sites)) +
    scale(log1p(polyasite_pas_count)) + scale(log1p(n_ends)) +
    scale(log1p(gene_length)) + scale(log1p(canonical_utr3_length)) +
    scale(mean_expression)
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
      fit <- glm(models[[model_name]], data = model_data, family = binomial())
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
