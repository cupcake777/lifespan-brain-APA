#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glue)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(as.character(x[[1]]))) y else x
}
script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE), error = function(e) NA_character_)
repo_root <- normalizePath(file.path(dirname(script_path %||% "tools"), ".."), mustWork = FALSE)
if (!file.exists(file.path(repo_root, "src", "utils_config.R"))) {
  repo_root <- normalizePath(getwd(), mustWork = FALSE)
}

global_cfg_path <- normalizePath(arg_value("--global", "config/global.yml"), mustWork = FALSE)
gwas_cfg_path <- normalizePath(arg_value("--gwas", "config/gwas.yml"), mustWork = FALSE)
qtl_cfg_path <- normalizePath(arg_value("--qtl", "config/qtl.yml"), mustWork = FALSE)
output_dir <- normalizePath(arg_value("--output-dir", "results/harmony_qc"), mustWork = FALSE)
output_html <- normalizePath(arg_value("--output", file.path(output_dir, "harmony_QC_report.html")), mustWork = FALSE)
sample_n <- as.integer(arg_value("--sample-n", "200000"))
dbsnp_sample_n <- as.integer(arg_value("--dbsnp-sample-n", "5000"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(repo_root, "src", "utils_config.R"))
source(file.path(repo_root, "src", "utils_format.R"))

cfg_global <- yaml::read_yaml(global_cfg_path)
cfg_gwas <- yaml::read_yaml(gwas_cfg_path)
cfg_qtl <- if (file.exists(qtl_cfg_path)) yaml::read_yaml(qtl_cfg_path) else list()
base_dir <- easycoloc_resolve_project_root(cfg_global, global_cfg_path)
cfg_global <- easycoloc_resolve_global_paths(cfg_global, base_dir)
cfg_gwas <- easycoloc_resolve_gwas_paths(cfg_gwas, base_dir)
target_build <- cfg_qtl$qtl_info$build %||% "hg38"
target_build <- tolower(as.character(target_build))
target_build_num <- easycoloc_normalize_build(target_build)
dbsnp_vcf <- if (identical(target_build_num, "19")) {
  cfg_global$dbsnp_hg19 %||% ""
} else {
  cfg_global$dbsnp_hg38 %||% ""
}

harmony_dir <- cfg_global$harmonize_dir
if (is.null(harmony_dir) || !dir.exists(harmony_dir)) {
  stop(glue("harmonize_dir not found: {harmony_dir}"))
}

dataset_by_id <- list()
for (ds in cfg_gwas$datasets) dataset_by_id[[ds$id]] <- ds

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

status_for <- function(value, pass, warn = NULL, higher_is_better = TRUE) {
  if (is.na(value)) return("fail")
  if (higher_is_better) {
    if (value >= pass) return("pass")
    if (!is.null(warn) && value >= warn) return("warn")
    return("fail")
  }
  if (value <= pass) return("pass")
  if (!is.null(warn) && value <= warn) return("warn")
  "fail"
}

metric_long <- function(summary_dt) {
  metrics <- c(
    "raw_schema_complete", "canonical_schema_complete", "n_complete_rate", "invalid_eaf_rate", "invalid_p_rate",
    "invalid_se_rate", "liftover_mapped_rate", "fasta_ref_match_rate", "dbsnp_position_match_rate",
    "rsid_available_rate", "variant_id_available_rate", "palindromic_rate"
  )
  metric_cols <- intersect(metrics, names(summary_dt))
  metric_dt <- copy(summary_dt[, c("dataset", metric_cols), with = FALSE])
  for (col in metric_cols) metric_dt[, (col) := as.numeric(get(col))]
  long <- melt(metric_dt, id.vars = "dataset", measure.vars = metric_cols,
               variable.name = "metric", value.name = "value")
  long[, status := fifelse(metric %in% c("invalid_eaf_rate", "invalid_p_rate", "invalid_se_rate"),
    fifelse(value <= 0, "pass", fifelse(value <= 0.001, "warn", "fail")),
    fifelse(metric %in% c("raw_schema_complete", "canonical_schema_complete"),
      fifelse(value >= 1, "pass", "fail"),
      fifelse(value >= 0.99, "pass", fifelse(value >= 0.95, "warn", "fail"))
    )
  )]
  long
}

save_plot_png <- function(plot, filename, width = 9, height = 5) {
  path <- file.path(fig_dir, filename)
  ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = 140)
  path
}

embed_png <- function(path) {
  encoded <- jsonlite::base64_enc(path)
  glue('<img src="data:image/png;base64,{encoded}" style="max-width:100%; height:auto;" alt="{html_escape(basename(path))}">')
}

check_dbsnp_position_concordance <- function(dt, dbsnp_vcf, build, dataset) {
  if (is.null(dbsnp_vcf) || !nzchar(dbsnp_vcf) || !file.exists(dbsnp_vcf) || dbsnp_sample_n <= 0) {
    return(list(match_rate = NA_real_, checked = 0L, mismatches = data.table()))
  }
  if (!all(c("CHR", "POS", "SNPID") %in% names(dt))) {
    return(list(match_rate = NA_real_, checked = 0L, mismatches = data.table()))
  }
  rs_ok <- !is.na(dt$SNPID) & nzchar(as.character(dt$SNPID)) &
    grepl("^rs[0-9]+$", as.character(dt$SNPID), ignore.case = TRUE) &
    !is.na(dt$CHR) & !is.na(dt$POS)
  candidates <- unique(dt[rs_ok, .(CHR = as.character(CHR), POS = as.integer(POS), SNPID = as.character(SNPID))])
  if (nrow(candidates) == 0) {
    return(list(match_rate = NA_real_, checked = 0L, mismatches = data.table()))
  }
  candidates <- candidates[seq_len(min(nrow(candidates), dbsnp_sample_n))]
  ref <- easycoloc_query_vcf_regions(
    dt = candidates,
    vcf_file = dbsnp_vcf,
    fields = c("%CHROM", "%POS", "%ID"),
    col_names = c("CHR", "POS", "rsID"),
    build = build,
    verbose = FALSE,
    label = "dbSNP position QC"
  )
  if (is.null(ref) || nrow(ref) == 0) {
    mismatches <- head(candidates[, .(
      dataset = dataset,
      SNPID,
      query_CHR = gsub("^chr", "", as.character(CHR), ignore.case = TRUE),
      query_POS = as.integer(POS),
      dbsnp_positions = NA_character_
    )], 25L)
    return(list(match_rate = 0, checked = nrow(candidates), mismatches = mismatches))
  }
  ref <- ref[grepl("^rs[0-9]+$", as.character(rsID), ignore.case = TRUE)]
  ref <- unique(ref[, .(SNPID = as.character(rsID), dbsnp_CHR = as.character(CHR), dbsnp_POS = as.integer(POS))])
  merged <- merge(candidates, ref, by = "SNPID", all.x = TRUE, allow.cartesian = TRUE)
  merged[, query_CHR := gsub("^chr", "", as.character(CHR), ignore.case = TRUE)]
  merged[, dbsnp_CHR_norm := gsub("^chr", "", as.character(dbsnp_CHR), ignore.case = TRUE)]
  merged[, position_match := !is.na(dbsnp_POS) & query_CHR == dbsnp_CHR_norm & as.integer(POS) == as.integer(dbsnp_POS)]
  by_snpid <- merged[, .(
    position_match = any(position_match, na.rm = TRUE),
    query_CHR = first(query_CHR),
    query_POS = first(as.integer(POS)),
    dbsnp_positions = if (all(is.na(dbsnp_POS))) NA_character_ else paste(unique(paste0(dbsnp_CHR_norm[!is.na(dbsnp_POS)], ":", dbsnp_POS[!is.na(dbsnp_POS)])), collapse = ";")
  ), by = SNPID]
  mismatches <- head(by_snpid[position_match == FALSE][,
    .(dataset = dataset, SNPID, query_CHR, query_POS, dbsnp_positions)
  ], 25L)
  list(
    match_rate = if (nrow(by_snpid) > 0) sum(by_snpid$position_match) / nrow(by_snpid) else NA_real_,
    checked = nrow(by_snpid),
    mismatches = mismatches
  )
}

qc_one_file <- function(file) {
  dataset <- sub("_b[0-9]+to[0-9]+_harmonized\\.tsv(\\.gz)?$", "", basename(file))
  ds <- dataset_by_id[[dataset]]
  header <- names(fread(file, nrows = 0, showProgress = FALSE))
  required <- easycoloc_harmonized_gwas_output_cols()
  raw_schema_complete <- all(required %in% header)
  raw_missing_required <- paste(setdiff(required, header), collapse = ",")
  dt <- fread(file, nrows = sample_n, showProgress = FALSE)
  dt <- easycoloc_standardize_harmonized_gwas(dt, sample_size_n = ds$sample_size_n %||% NA_real_)
  n <- nrow(dt)
  canonical_schema_complete <- all(required %in% names(dt))
  canonical_missing_required <- paste(setdiff(required, names(dt)), collapse = ",")

  invalid_eaf <- if ("EAF" %in% names(dt)) sum(is.na(dt$EAF) | dt$EAF < 0 | dt$EAF > 1) else n
  invalid_p <- if ("P" %in% names(dt)) sum(is.na(dt$P) | dt$P < 0 | dt$P > 1) else n
  invalid_se <- if ("SE" %in% names(dt)) sum(is.na(dt$SE) | dt$SE <= 0) else n
  missing_n <- if ("N" %in% names(dt)) sum(is.na(dt$N)) else n
  has_rsid <- if ("SNPID" %in% names(dt)) sum(!is.na(dt$SNPID) & nzchar(as.character(dt$SNPID)) & grepl("^rs[0-9]+$", as.character(dt$SNPID), ignore.case = TRUE)) else 0L
  has_varid <- if ("variant_id" %in% names(dt)) sum(!is.na(dt$variant_id) & nzchar(as.character(dt$variant_id))) else 0L
  pal <- if (all(c("EA", "NEA") %in% names(dt))) sum(easycoloc_is_palindromic(dt$EA, dt$NEA), na.rm = TRUE) else NA_integer_
  liftover_mapped <- if ("LIFTOVER_MAPPED" %in% names(dt)) sum(dt$LIFTOVER_MAPPED %in% TRUE, na.rm = TRUE) else NA_integer_
  ref_match <- if ("REF_MATCH" %in% names(dt)) sum(dt$REF_MATCH %in% TRUE, na.rm = TRUE) else NA_integer_
  ref_checkable <- if ("REF_MATCH" %in% names(dt)) sum(!is.na(dt$REF_MATCH)) else NA_integer_
  dbsnp_qc <- check_dbsnp_position_concordance(dt, dbsnp_vcf = dbsnp_vcf, build = target_build_num, dataset = dataset)

  summary <- data.table(
    dataset = dataset,
    file = basename(file),
    rows_sampled = n,
    original_schema = paste(header, collapse = ", "),
    raw_schema_complete = as.integer(raw_schema_complete),
    raw_missing_required = raw_missing_required,
    canonical_schema_complete = as.integer(canonical_schema_complete),
    canonical_missing_required = canonical_missing_required,
    n_missing = missing_n,
    n_complete_rate = if (n > 0) 1 - missing_n / n else NA_real_,
    invalid_eaf_rate = if (n > 0) invalid_eaf / n else NA_real_,
    invalid_p_rate = if (n > 0) invalid_p / n else NA_real_,
    invalid_se_rate = if (n > 0) invalid_se / n else NA_real_,
    rsid_available_rate = if (n > 0) has_rsid / n else NA_real_,
    variant_id_available_rate = if (n > 0) has_varid / n else NA_real_,
    palindromic_rate = if (n > 0) pal / n else NA_real_,
    liftover_mapped_rate = if (!is.na(liftover_mapped) && n > 0) liftover_mapped / n else NA_real_,
    fasta_ref_match_rate = if (!is.na(ref_match) && ref_checkable > 0) ref_match / ref_checkable else NA_real_,
    dbsnp_position_checked = dbsnp_qc$checked,
    dbsnp_position_match_rate = dbsnp_qc$match_rate
  )

  flags <- dt[, .(
    dataset = dataset,
    SNPID, variant_id, CHR, POS, EA, NEA, EAF, BETA, SE, P, N,
    invalid_EAF = is.na(EAF) | EAF < 0 | EAF > 1,
    invalid_P = is.na(P) | P < 0 | P > 1,
    invalid_SE = is.na(SE) | SE <= 0,
    missing_N = is.na(N),
    palindromic = easycoloc_is_palindromic(EA, NEA)
  )]
  list(summary = summary, flags = flags, dbsnp_mismatches = dbsnp_qc$mismatches)
}

files <- sort(unique(c(
  Sys.glob(file.path(harmony_dir, "*_harmonized.tsv.gz")),
  Sys.glob(file.path(harmony_dir, "*_harmonized.tsv"))
)))
if (length(files) == 0) stop(glue("No harmonized TSV/TSV.GZ files found in {harmony_dir}"))

qc <- lapply(files, qc_one_file)
summary_dt <- rbindlist(lapply(qc, `[[`, "summary"), fill = TRUE)
flags_dt <- rbindlist(lapply(qc, `[[`, "flags"), fill = TRUE)
dbsnp_mismatch_dt <- rbindlist(lapply(qc, `[[`, "dbsnp_mismatches"), fill = TRUE)
summary_dt[, overall_status := "pass"]
summary_dt[raw_schema_complete < 1 | canonical_schema_complete < 1 | n_complete_rate < 1 | invalid_eaf_rate > 0 | invalid_p_rate > 0 | invalid_se_rate > 0, overall_status := "fail"]
summary_dt[overall_status != "fail" & (
  (!is.na(liftover_mapped_rate) & liftover_mapped_rate < 0.99) |
    (!is.na(fasta_ref_match_rate) & fasta_ref_match_rate < 0.99) |
    (!is.na(dbsnp_position_match_rate) & dbsnp_position_match_rate < 0.99) |
    rsid_available_rate < 0.95 |
    variant_id_available_rate < 1
), overall_status := "warn"]

summary_file <- file.path(output_dir, "harmony_QC_summary.tsv")
flags_file <- file.path(output_dir, "harmony_QC_per_variant_flags.tsv.gz")
dbsnp_mismatch_file <- file.path(output_dir, "harmony_QC_dbsnp_position_mismatches.tsv")
fwrite(summary_dt, summary_file, sep = "\t", na = "NA", quote = FALSE)
fwrite(flags_dt, flags_file, sep = "\t", na = "NA", quote = FALSE)
fwrite(dbsnp_mismatch_dt, dbsnp_mismatch_file, sep = "\t", na = "NA", quote = FALSE)

long <- metric_long(summary_dt)
p_heat <- ggplot(long, aes(metric, dataset, fill = status)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(pass = "#2E7D32", warn = "#F9A825", fail = "#C62828"), drop = FALSE) +
  labs(x = NULL, y = NULL, fill = "QC", title = "Harmony QC Status Matrix") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
heat_png <- save_plot_png(p_heat, "01_qc_heatmap.png", width = 11, height = max(4, 0.32 * nrow(summary_dt) + 2))

retention_dt <- summary_dt[, .(
  dataset,
  retained = as.numeric(rows_sampled),
  missing_N = as.numeric(n_missing),
  invalid_EAF = as.numeric(round(invalid_eaf_rate * rows_sampled)),
  invalid_P = as.numeric(round(invalid_p_rate * rows_sampled)),
  invalid_SE = as.numeric(round(invalid_se_rate * rows_sampled))
)]
ret_long <- melt(retention_dt, id.vars = "dataset", variable.name = "category", value.name = "count")
p_ret <- ggplot(ret_long, aes(dataset, count, fill = category)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Sampled variants / QC counts", fill = NULL, title = "Variant Retention and Basic QC Counts") +
  theme_minimal(base_size = 11)
ret_png <- save_plot_png(p_ret, "02_variant_retention.png", width = 9, height = max(4, 0.28 * nrow(summary_dt) + 2))

eaf_plot_dt <- flags_dt[!is.na(EAF)]
eaf_plot_dt[, dataset := factor(dataset, levels = summary_dt$dataset)]
p_eaf <- ggplot(eaf_plot_dt, aes(EAF, fill = palindromic)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  facet_wrap(~dataset, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#E15759")) +
  labs(x = "Effect allele frequency", y = "Variants", fill = "Palindromic", title = "EAF Distribution by Dataset") +
  theme_minimal(base_size = 10)
eaf_png <- save_plot_png(p_eaf, "03_eaf_distribution.png", width = 11, height = max(5, ceiling(nrow(summary_dt) / 3) * 2.1))

schema_dt <- summary_dt[, .(dataset, raw_schema_complete, raw_missing_required, canonical_schema_complete, canonical_missing_required, original_schema)]
table_rows <- apply(summary_dt[, .(
  dataset, overall_status, rows_sampled, raw_schema_complete, canonical_schema_complete, n_complete_rate,
  invalid_eaf_rate, invalid_p_rate, invalid_se_rate, rsid_available_rate,
  variant_id_available_rate, liftover_mapped_rate, fasta_ref_match_rate,
  dbsnp_position_checked, dbsnp_position_match_rate, raw_missing_required
)], 1, function(row) {
  status <- row[["overall_status"]]
  paste0(
    "<tr class='", html_escape(status), "'>",
    paste(sprintf("<td>%s</td>", html_escape(row)), collapse = ""),
    "</tr>"
  )
})
table_header <- paste(sprintf("<th>%s</th>", names(summary_dt[, .(
  dataset, overall_status, rows_sampled, raw_schema_complete, canonical_schema_complete, n_complete_rate,
  invalid_eaf_rate, invalid_p_rate, invalid_se_rate, rsid_available_rate,
  variant_id_available_rate, liftover_mapped_rate, fasta_ref_match_rate,
  dbsnp_position_checked, dbsnp_position_match_rate, raw_missing_required
)])), collapse = "")

schema_rows <- apply(schema_dt, 1, function(row) {
  paste0("<tr>", paste(sprintf("<td>%s</td>", html_escape(row)), collapse = ""), "</tr>")
})

generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
html <- glue('
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>EasyColoc Harmony QC Report</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 28px; color: #202124; }}
h1, h2 {{ margin-bottom: 0.35rem; }}
.meta {{ color: #5f6368; margin-top: 0; }}
.summary {{ display: grid; grid-template-columns: repeat(4, minmax(120px, 1fr)); gap: 12px; margin: 18px 0; }}
.card {{ border: 1px solid #dadce0; border-radius: 6px; padding: 12px; }}
.value {{ font-size: 24px; font-weight: 650; }}
table {{ border-collapse: collapse; width: 100%; font-size: 12px; margin: 12px 0 24px; }}
th, td {{ border: 1px solid #e0e0e0; padding: 6px 8px; text-align: left; vertical-align: top; }}
th {{ background: #f8f9fa; position: sticky; top: 0; }}
tr.pass td:first-child, .pass-badge {{ border-left: 5px solid #2E7D32; }}
tr.warn td:first-child, .warn-badge {{ border-left: 5px solid #F9A825; }}
tr.fail td:first-child, .fail-badge {{ border-left: 5px solid #C62828; }}
.note {{ background: #f8f9fa; border-left: 4px solid #5f6368; padding: 10px 12px; margin: 12px 0; }}
.figure {{ margin: 22px 0; border: 1px solid #e0e0e0; border-radius: 6px; padding: 12px; overflow-x: auto; }}
a {{ color: #174ea6; }}
code {{ background: #f1f3f4; padding: 1px 4px; border-radius: 3px; }}
</style>
</head>
<body>
<h1>EasyColoc Harmony QC Report</h1>
<p class="meta">Generated: {generated_at}<br>Harmony directory: <code>{html_escape(harmony_dir)}</code><br>Target build: <code>{html_escape(target_build)}</code><br>Rows sampled per dataset: <code>{sample_n}</code><br>dbSNP position check sample per dataset: <code>{dbsnp_sample_n}</code></p>

<div class="summary">
  <div class="card pass-badge"><div class="value">{sum(summary_dt$overall_status == "pass")}</div><div>Pass datasets</div></div>
  <div class="card warn-badge"><div class="value">{sum(summary_dt$overall_status == "warn")}</div><div>Warn datasets</div></div>
  <div class="card fail-badge"><div class="value">{sum(summary_dt$overall_status == "fail")}</div><div>Fail datasets</div></div>
  <div class="card"><div class="value">{nrow(summary_dt)}</div><div>Total datasets</div></div>
</div>

<div class="note">
This report checks whether harmonized GWAS caches are suitable as shared upstream inputs for coloc, SMR, and sLDSC. 
It evaluates the on-disk schema, EasyColoc canonical readability, N/EAF/P/SE validity, SNPID and chr:pos:ref:alt availability, optional liftover/reference-match markers when present, and sampled SNPID-to-position concordance against the configured target-build dbSNP reference.
</div>

<h2>QC Status Matrix</h2>
<div class="figure">{embed_png(heat_png)}</div>

<h2>Variant Retention and Basic QC</h2>
<div class="figure">{embed_png(ret_png)}</div>

<h2>EAF Distribution</h2>
<div class="figure">{embed_png(eaf_png)}</div>

<h2>Dataset Summary</h2>
<p>Download: <a href="{html_escape(basename(summary_file))}">harmony_QC_summary.tsv</a>, <a href="{html_escape(basename(flags_file))}">harmony_QC_per_variant_flags.tsv.gz</a>, and <a href="{html_escape(basename(dbsnp_mismatch_file))}">harmony_QC_dbsnp_position_mismatches.tsv</a></p>
<table><thead><tr>{table_header}</tr></thead><tbody>{paste(table_rows, collapse = "\n")}</tbody></table>

<h2>Original Schema Audit</h2>
<table><thead><tr><th>dataset</th><th>raw_schema_complete</th><th>raw_missing_required</th><th>canonical_schema_complete</th><th>canonical_missing_required</th><th>original_schema</th></tr></thead><tbody>{paste(schema_rows, collapse = "\n")}</tbody></table>

</body>
</html>
')

writeLines(html, output_html)
message(glue("[Harmony QC] Wrote HTML report: {output_html}"))
message(glue("[Harmony QC] Wrote summary: {summary_file}"))
message(glue("[Harmony QC] Wrote per-variant flags: {flags_file}"))
message(glue("[Harmony QC] Wrote dbSNP position mismatches: {dbsnp_mismatch_file}"))
