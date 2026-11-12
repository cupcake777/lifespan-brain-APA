# Export accepted interaction pairs and the MAGMA group table.
# Reads nominal interaction hits and HC3 core groups; writes data/05_gxs/ and data/07_coloc/intx.signif_pairs.txt.gz.
# --- load packages ---
suppressMessages(library(data.table))

# --- config ---
data_dir <- file.path(Sys.getenv("APA_ROOT"), "data/05_gxs")
core_file <- file.path(data_dir, "intx_hc3_core_groups.tsv")
hits_file <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/03_interaction/results/01_tensorqtl_nominal/sig.interaction.result")
out_hits <- file.path(data_dir, "sig.interaction.result")
out_pairs <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc/intx.signif_pairs.txt.gz")
out_groups <- file.path(data_dir, "accepted_intx_core_groups.tsv")

# --- load accepted group set and nominal hits ---
core <- fread(core_file)
hits <- fread(hits_file)
stopifnot(nrow(core) == 193L, uniqueN(core$group_id) == nrow(core))

# --- retain every nominal variant-phenotype hit within accepted groups ---
hits[, group_id := tstrsplit(phenotype_id, "|", fixed = TRUE)[[2L]]]
accepted <- hits[group_id %in% core$group_id]
group_hits <- accepted[, .(n_nominal_hits = .N), by = group_id]
magma_groups <- merge(
  core[, .(group_id, phenotype_id = lead_phenotype_id)],
  group_hits,
  by = "group_id",
  all.x = TRUE
)
magma_groups[is.na(n_nominal_hits), n_nominal_hits := 0L]
stopifnot(!anyDuplicated(accepted[, .(phenotype_id, variant_id)]))

# --- format interaction pairs for MAGMA and sLDSC ---
variant_parts <- tstrsplit(accepted$variant_id, ":", fixed = TRUE)
stopifnot(length(variant_parts) == 4L, all(vapply(variant_parts, function(x) !anyNA(x), logical(1))))
stopifnot(all(is.finite(accepted$pval_gi)), all(is.finite(accepted$b_gi)),
          all(is.finite(accepted$b_gi_se)), all(accepted$b_gi_se > 0))

signif_pairs <- data.table(
  chr = variant_parts[[1L]],
  pos = as.integer(variant_parts[[2L]]),
  ref = variant_parts[[3L]],
  alt = variant_parts[[4L]],
  phenotype_id = accepted$phenotype_id,
  variant_id = accepted$variant_id,
  start_distance = accepted$start_distance,
  end_distance = accepted$end_distance,
  af = accepted$af,
  ma_samples = accepted$ma_samples,
  ma_count = accepted$ma_count,
  pval_nominal = accepted$pval_gi,
  slope = accepted$b_gi,
  slope_se = accepted$b_gi_se
)
stopifnot(all(is.finite(signif_pairs$pos)))

# --- write outputs ---
hit_columns <- setdiff(names(hits), "group_id")
fwrite(accepted[, ..hit_columns], out_hits, sep = "\t")
fwrite(signif_pairs, out_pairs, sep = "\t", compress = "gzip")
fwrite(magma_groups, out_groups, sep = "\t")

cat(sprintf(
  "accepted groups=%d; groups with nominal hits=%d; groups without nominal hits=%d; retained interaction pairs=%d; unique variants=%d\n",
  nrow(core), uniqueN(accepted$group_id), sum(magma_groups$n_nominal_hits == 0L),
  nrow(accepted), uniqueN(accepted$variant_id)
))
