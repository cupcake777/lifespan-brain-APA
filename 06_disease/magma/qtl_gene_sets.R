#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

# Build MAGMA row-based APA-QTL gene-set definitions from significant APA-QTL
# genes. These sets are used with genome-wide MAGMA gene results:
#   magma --gene-results <baseline.genes.raw> --set-annot <sets> --out <prefix>

# --- load packages ---

library(data.table)
library(org.Hs.eg.db)

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
out_dir <- file.path(apa_root, "data", "06_disease")
gene_loc_path <- file.path(apa_root, "ref", "magma", "NCBI38.gene.loc")

qtl_tasks <- list(
  prenatal = file.path(apa_root, "data", "04_qtl", "prenatal.sig_QTL.txt.gz"),
  postnatal = file.path(apa_root, "data", "04_qtl", "postnatal.sig_QTL.txt.gz"),
  interaction = file.path(apa_root, "data", "05_gxs", "sig.interaction.result")
)

out_sets <- file.path(out_dir, "apa_qtl_gene_sets.sets")
out_summary <- file.path(out_dir, "apa_qtl_gene_sets.sum.tsv")
out_gene_table <- file.path(out_dir, "apa_qtl_gene_sets.genes.tsv")
out_overlap <- file.path(out_dir, "apa_qtl_gene_sets.overlap.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- load gene ID to Entrez map ---

symbol_map <- as.data.table(stack(as.list(org.Hs.egSYMBOL2EG)))
setnames(symbol_map, c("gene_id", "gene_key"))

ensembl_map <- as.data.table(stack(as.list(org.Hs.egENSEMBL2EG)))
setnames(ensembl_map, c("gene_id", "gene_key"))

db_map <- unique(rbindlist(list(symbol_map, ensembl_map), use.names = TRUE))
db_map[, gene_id := as.character(gene_id)]
setkey(db_map, gene_key)

# --- load valid MAGMA genes ---

loc_dict <- fread(gene_loc_path, col.names = c("gene_id", "chr", "start", "stop", "strand", "symbol"))
loc_dict[, gene_id := as.character(gene_id)]
loc_dict <- loc_dict[!is.na(start) & !is.na(stop) & stop > start]
loc_dict <- loc_dict[, .(gene_id, symbol)]
setkey(loc_dict, gene_id)

# --- build APA-QTL gene sets ---

set_rows <- list()
summary_rows <- list()
gene_rows <- list()

for (set_name in names(qtl_tasks)) {
  message("Processing ", set_name)

  qtl_data <- fread(qtl_tasks[[set_name]])
  qtl_data[, gene_key := tstrsplit(phenotype_id, "|", fixed = TRUE)[[2]]]

  mapped <- merge(qtl_data[, .(gene_key)], db_map, by = "gene_key", all.x = TRUE)
  mapped <- merge(mapped[!is.na(gene_id)], loc_dict, by = "gene_id", all.x = FALSE)

  genes <- sort(unique(mapped$gene_id))

  set_rows[[set_name]] <- data.table(
    set_name = set_name,
    genes = paste(genes, collapse = " ")
  )

  summary_rows[[set_name]] <- data.table(
    set_name = set_name,
    n_input_qtl_rows = nrow(qtl_data),
    n_mapped_rows = nrow(mapped),
    n_unique_entrez_genes = length(genes)
  )

  gene_rows[[set_name]] <- unique(mapped[, .(set_name, gene_id, symbol)])
}

sets <- rbindlist(set_rows, use.names = TRUE)
summary <- rbindlist(summary_rows, use.names = TRUE)
gene_table <- rbindlist(gene_rows, use.names = TRUE)

fwrite(sets, out_sets, sep = " ", col.names = FALSE, quote = FALSE)
fwrite(summary, out_summary, sep = "\t")
fwrite(gene_table, out_gene_table, sep = "\t")

# --- summarize gene-set overlap ---

wide <- dcast(unique(gene_table[, .(gene_id, set_name)]), gene_id ~ set_name, fun.aggregate = length)
for (col in setdiff(names(wide), "gene_id")) {
  wide[[col]] <- wide[[col]] > 0
}

overlap <- CJ(set_a = names(qtl_tasks), set_b = names(qtl_tasks))
overlap[, n_overlap := mapply(
  function(a, b) sum(wide[[a]] & wide[[b]]),
  set_a,
  set_b
)]

fwrite(overlap, out_overlap, sep = "\t")

message("Done.")
