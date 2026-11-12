#!/usr/bin/env Rscript
# Screen exact APA-QTL pairs with a strong effect in one stage and a near-zero effect in the other.
# Reads the significant-pair, permutation and nominal QTL tables in data/04_qtl; writes the screen, shortlist and raw PDUI summary.

# --- load packages ---

suppressPackageStartupMessages(library(data.table))

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
qtl_dir <- file.path(apa_root, "data", "04_qtl")
pre_dir <- file.path(qtl_dir, "pre_input")
prenatal_sig_path <- file.path(qtl_dir, "prenatal.sig_QTL.txt.gz")
postnatal_sig_path <- file.path(qtl_dir, "postnatal.sig_QTL.txt.gz")
prenatal_nominal_path <- file.path(qtl_dir, "nominal", "prenatal.cis_qtl_pairs.txt.gz")
postnatal_nominal_path <- file.path(qtl_dir, "nominal", "postnatal.cis_qtl_pairs.txt.gz")
prenatal_perm_path <- file.path(qtl_dir, "prenatal.cis_qtl.txt.gz")
postnatal_perm_path <- file.path(qtl_dir, "postnatal.cis_qtl.txt.gz")
pas_info_path <- file.path(pre_dir, "proximal.apa.info")
output_path <- file.path(qtl_dir, "figure4_stage_specific_pair_screen.tsv")
shortlist_path <- file.path(qtl_dir, "figure4_stage_specific_pair_shortlist.tsv")
raw_pdui_summary_path <- file.path(qtl_dir, "figure4_stage_specific_raw_pdui_summary.tsv")
raw_pdui_root <- file.path(apa_root, "raw", "quant", "apa", "gencode", "all")
prenatal_pfile <- file.path(pre_dir, "genotype.prenatal")
postnatal_pfile <- file.path(pre_dir, "genotype.postnatal")
prenatal_samples_path <- file.path(pre_dir, "samples.prenatal.txt")
postnatal_samples_path <- file.path(pre_dir, "samples.postnatal.txt")

target_abs_beta_min <- 0.50
other_abs_beta_max <- 0.10
other_p_min <- 0.50
phenotype_q_max <- 0.05
effect_difference_fdr_max <- 0.05
figure_extension_min_bp <- 300L
shortlist_events_per_stage <- 30L

# --- helpers ---

make_phenotype_key <- function(x) {
    vapply(strsplit(as.character(x), "|", fixed = TRUE), function(parts) {
        paste(parts[c(1L, 3L, 4L)], collapse = "|")
    }, character(1))
}

add_pair_keys <- function(dt) {
    dt[, phenotype_key := make_phenotype_key(phenotype_id)]
    x <- sub("_b[0-9]+$", "", as.character(dt$variant_id))
    m <- regexec("^(chr[^:_]+)[_:]([0-9]+)[_:]([^_:]+)[_:]([^_:]+)$", x)
    mm <- regmatches(x, m)
    ok <- lengths(mm) == 5L
    x[ok] <- vapply(mm[ok], function(z) paste(z[2], z[3], z[4], z[5], sep = ":"), character(1))
    dt[, variant_key := x]
    dt[, allele_phenotype_key := paste(phenotype_key, variant_key, sep = "||")]
    dt
}

stream_nominal_matches <- function(path, key_dt) {
    key_path <- tempfile(pattern = "figure4_pair_keys_", fileext = ".tsv")
    transcript_path <- tempfile(pattern = "figure4_transcripts_", fileext = ".txt")
    fwrite(unique(key_dt[, .(allele_phenotype_key)]), key_path, sep = "\t", col.names = FALSE)
    transcript_ids <- unique(tstrsplit(key_dt$phenotype_key, "|", fixed = TRUE)[[1]])
    writeLines(paste0(transcript_ids, "|"), transcript_path)
    on.exit(unlink(c(key_path, transcript_path)), add = TRUE)

    header <- names(fread(cmd = paste("gzip -cd", shQuote(path), "| head -n 2"), nrows = 0L))
    required <- c("phenotype_id", "variant_id", "pval_nominal", "slope", "slope_se",
                  "af", "ma_samples", "ma_count")
    column_index <- match(required, header)

    awk_script <- paste0(
        "awk 'BEGIN{FS=OFS=\"\\t\"} ",
        "FNR==NR{k[$1]=1; next} ",
        "{key=phenotype_key($", column_index[1], ") \"||\" variant_key($", column_index[2],
        "); if(key in k){print key,$", column_index[3], ",$", column_index[4], ",$", column_index[5],
        ",$", column_index[6], ",$", column_index[7], ",$", column_index[8], "}} ",
        "function phenotype_key(pid,parts){split(pid,parts,\"\\\\|\"); ",
        "return parts[1] \"|\" parts[3] \"|\" parts[4]} ",
        "function variant_key(v,x,parts){x=v; sub(/_b[0-9]+$/,\"\",x); ",
        "gsub(/_/,\":\",x); split(x,parts,\":\"); ",
        "return parts[1] \":\" parts[2] \":\" parts[3] \":\" parts[4]}' ",
        shQuote(key_path), " -"
    )
    data_command <- paste(
        "gzip -cd", shQuote(path), "| grep -F -f", shQuote(transcript_path), "|", awk_script
    )
    matched <- fread(
        cmd = data_command,
        header = FALSE,
        sep = "\t",
        col.names = c("allele_phenotype_key", "other_p", "other_beta", "other_se",
                      "other_af", "other_ma_samples", "other_ma_count")
    )
    setorder(matched, other_p)
    unique(matched, by = "allele_phenotype_key")
}

raw_pdui_cache <- new.env(parent = emptyenv())

value_at_dosage <- function(dosage, value, target_dosage) {
    hit <- value[dosage == target_dosage]
    if (length(hit) == 0L) NA_real_ else hit[1L]
}

summarize_raw_pdui <- function(candidate, stage, samples, pfile_prefix) {
    chromosome <- sub(":.*$", "", candidate$region)
    raw_path <- file.path(
        raw_pdui_root,
        paste0("Dapars2_out_", chromosome),
        paste0("Dapars2_result_temp.", chromosome, ".txt")
    )
    raw <- fread(raw_path)
    raw_gene_key <- paste(
        candidate$transcript, candidate$gene, chromosome, candidate$strand,
        sep = "|"
    )
    if (exists(raw_gene_key, envir = raw_pdui_cache, inherits = FALSE)) {
        pdui <- get(raw_gene_key, envir = raw_pdui_cache, inherits = FALSE)
    } else {
        raw <- raw[Gene == raw_gene_key]
        if (nrow(raw) != 1L) stop("Expected one raw PDUI row for ", raw_gene_key)
        sample_columns <- names(raw)[-(1:4)]
        sample_ids <- sub("_PDUI$", "", basename(sample_columns))
        pdui <- data.table(
            sample_id = sample_ids,
            raw_pdui = as.numeric(unlist(raw[1, ..sample_columns]))
        )
        assign(raw_gene_key, pdui, envir = raw_pdui_cache)
    }

    fields <- strsplit(candidate$variant_id, ":", fixed = TRUE)[[1]]
    if (length(fields) != 4L) stop("Unexpected variant ID: ", candidate$variant_id)
    tmp_bed <- tempfile(fileext = ".bed")
    tmp_out <- tempfile(pattern = "figure4_dosage_")
    writeLines(sprintf("%s\t%s\t%s", fields[1], fields[2], fields[2]), tmp_bed)
    on.exit(unlink(c(tmp_bed, paste0(tmp_out, c(".vcf", ".log")))), add = TRUE)
    status <- system2(
        "plink2",
        c("--pfile", pfile_prefix, "--extract", "bed1", tmp_bed,
          "--export", "vcf", "--out", tmp_out, "--threads", "2"),
        stdout = FALSE,
        stderr = FALSE
    )
    vcf_path <- paste0(tmp_out, ".vcf")
    dosage <- data.table()
    if (status == 0L && file.exists(vcf_path)) {
        vcf <- fread(vcf_path, skip = "#CHROM", colClasses = "character")
        if (nrow(vcf) == 1L) {
            sample_columns <- setdiff(
                names(vcf),
                c("#CHROM", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
            )
            dosage <- vapply(sample_columns, function(sample_column) {
                genotype <- strsplit(vcf[[sample_column]][1], ":", fixed = TRUE)[[1]][1]
                if (is.na(genotype) || genotype %in% c("", ".", "./.", ".|.")) return(NA_integer_)
                alleles <- as.integer(strsplit(genotype, "[|/]")[[1]])
                if (anyNA(alleles)) return(NA_integer_)
                sum(alleles)
            }, integer(1))
            if (vcf$REF[1] == fields[4] && vcf$ALT[1] == fields[3]) dosage <- 2L - dosage
            if (!(
                (vcf$REF[1] == fields[3] && vcf$ALT[1] == fields[4]) ||
                (vcf$REF[1] == fields[4] && vcf$ALT[1] == fields[3])
            )) stop("VCF alleles do not match ", candidate$variant_id)
            dosage <- data.table(sample_id = sample_columns, dosage = dosage)
        }
    }

    matched <- merge(
        data.table(sample_id = samples),
        merge(pdui, dosage, by = "sample_id"),
        by = "sample_id"
    )[is.finite(raw_pdui) & !is.na(dosage)]
    if (nrow(matched) == 0L) return(data.table())
    group_summary <- matched[, .(
        n = .N,
        median_raw_pdui = median(raw_pdui),
        mean_raw_pdui = mean(raw_pdui),
        sd_raw_pdui = sd(raw_pdui)
    ), by = dosage]
    model <- lm(raw_pdui ~ dosage, data = matched)
    model_coefficients <- coef(base::summary(model))
    dosage_p <- if ("dosage" %in% rownames(model_coefficients)) {
        model_coefficients["dosage", "Pr(>|t|)"]
    } else {
        NA_real_
    }
    group_summary[, `:=`(
        target_stage = candidate$target_stage,
        gene = candidate$gene,
        transcript = candidate$transcript,
        phenotype_id = candidate$phenotype_id,
        variant_id = candidate$variant_id,
        measured_stage = stage,
        raw_pdui_slope = unname(coef(model)["dosage"]),
        raw_pdui_p = dosage_p
    )]
    group_summary
}

# --- load or build the exact cross-stage pair screen ---

if (file.exists(output_path)) {
    screen <- fread(output_path)
} else {

    # --- load significant pairs and phenotype-level q values ---

    prenatal_sig <- add_pair_keys(fread(prenatal_sig_path))
    postnatal_sig <- add_pair_keys(fread(postnatal_sig_path))

    prenatal_sig <- unique(prenatal_sig[order(pval_nominal)], by = "allele_phenotype_key")
    postnatal_sig <- unique(postnatal_sig[order(pval_nominal)], by = "allele_phenotype_key")

    prenatal_q <- fread(prenatal_perm_path, select = c("phenotype_id", "qval"))
    postnatal_q <- fread(postnatal_perm_path, select = c("phenotype_id", "qval"))
    prenatal_q[, phenotype_key := make_phenotype_key(phenotype_id)]
    postnatal_q[, phenotype_key := make_phenotype_key(phenotype_id)]
    stage_q <- merge(
        unique(prenatal_q[, .(phenotype_key, prenatal_q = qval)], by = "phenotype_key"),
        unique(postnatal_q[, .(phenotype_key, postnatal_q = qval)], by = "phenotype_key"),
        by = "phenotype_key"
    )

    # --- match prenatal-significant pairs in postnatal nominal QTL ---

    message("Scanning postnatal nominal QTL for prenatal-significant pairs")
    prenatal_in_postnatal <- stream_nominal_matches(postnatal_nominal_path, prenatal_sig)
    prenatal_screen <- merge(prenatal_sig, prenatal_in_postnatal, by = "allele_phenotype_key")
    prenatal_screen <- merge(prenatal_screen, stage_q, by = "phenotype_key")
    prenatal_screen[, `:=`(
        target_stage = "prenatal",
        target_p = pval_nominal,
        target_beta = slope,
        target_se = slope_se,
        target_af = af,
        target_ma_samples = ma_samples,
        target_ma_count = ma_count,
        target_q = prenatal_q,
        other_q = postnatal_q,
        prenatal_beta = slope,
        prenatal_se = slope_se,
        postnatal_beta = other_beta,
        postnatal_se = other_se
    )]

    # --- match postnatal-significant pairs in prenatal nominal QTL ---

    message("Scanning prenatal nominal QTL for postnatal-significant pairs")
    postnatal_in_prenatal <- stream_nominal_matches(prenatal_nominal_path, postnatal_sig)
    postnatal_screen <- merge(postnatal_sig, postnatal_in_prenatal, by = "allele_phenotype_key")
    postnatal_screen <- merge(postnatal_screen, stage_q, by = "phenotype_key")
    postnatal_screen[, `:=`(
        target_stage = "postnatal",
        target_p = pval_nominal,
        target_beta = slope,
        target_se = slope_se,
        target_af = af,
        target_ma_samples = ma_samples,
        target_ma_count = ma_count,
        target_q = postnatal_q,
        other_q = prenatal_q,
        prenatal_beta = other_beta,
        prenatal_se = other_se,
        postnatal_beta = slope,
        postnatal_se = slope_se
    )]

    # --- cross-stage effect contrast ---

    screen <- rbindlist(list(prenatal_screen, postnatal_screen), fill = TRUE)
    screen[, delta_beta_pre_minus_post := prenatal_beta - postnatal_beta]
    screen[, delta_se := sqrt(prenatal_se^2 + postnatal_se^2)]
    screen[, z_effect_difference := delta_beta_pre_minus_post / delta_se]
    screen[, p_effect_difference := 2 * pnorm(-abs(z_effect_difference))]
    screen[, p_effect_difference_fdr := p.adjust(p_effect_difference, method = "BH"), by = target_stage]
    screen[, `:=`(
        target_abs_beta = abs(target_beta),
        other_abs_beta = abs(other_beta),
        effect_contrast = abs(target_beta) - abs(other_beta),
        target_phenotype_significant = target_q < phenotype_q_max,
        other_phenotype_nonsignificant = other_q >= phenotype_q_max,
        other_near_zero = abs(other_beta) <= other_abs_beta_max & other_p >= other_p_min,
        strong_target = abs(target_beta) >= target_abs_beta_min,
        effect_difference_significant = p_effect_difference_fdr < effect_difference_fdr_max
    )]
    screen[, evidence_class := fcase(
        strong_target & other_near_zero & target_phenotype_significant &
            other_phenotype_nonsignificant & effect_difference_significant,
            "strict_stage_specific",
        strong_target & other_near_zero & effect_difference_significant,
            "exact_pair_specific",
        effect_difference_significant,
            "effect_contrast_only",
        default = "no_cross_stage_evidence"
    )]

    # --- annotate APA event length ---

    parts <- tstrsplit(screen$phenotype_id, "|", fixed = TRUE)
    screen[, `:=`(
        transcript = parts[[1]],
        gene = parts[[2]],
        region = parts[[3]],
        strand = parts[[4]]
    )]
    region_text <- sub("^chr[^:]+:", "", screen$region)
    region_parts <- tstrsplit(region_text, "-", fixed = TRUE)
    screen[, `:=`(
        region_start = as.integer(region_parts[[1]]),
        region_end = as.integer(region_parts[[2]])
    )]
    pas_info <- fread(pas_info_path, select = c("phenotype_id", "Predicted_Proximal_APA"))
    setnames(pas_info, "Predicted_Proximal_APA", "proximal_pas")
    screen <- merge(screen, unique(pas_info, by = "phenotype_id"), by = "phenotype_id", all.x = TRUE)
    screen[, `:=`(
        common_utr_bp = fifelse(strand == "+", proximal_pas - region_start, region_end - proximal_pas),
        long_extension_bp = fifelse(strand == "+", region_end - proximal_pas, proximal_pas - region_start)
    )]
    screen[, figure_extension_sufficient := long_extension_bp >= figure_extension_min_bp]

    # --- rank and save evidence table ---

    screen[, evidence_rank := fcase(
        evidence_class == "strict_stage_specific", 1L,
        evidence_class == "exact_pair_specific", 2L,
        evidence_class == "effect_contrast_only", 3L,
        default = 4L
    )]
    setorder(screen, target_stage, evidence_rank, -figure_extension_sufficient,
             other_abs_beta, -target_abs_beta, p_effect_difference_fdr)

    keep_columns <- c(
        "target_stage", "evidence_class", "gene", "transcript", "phenotype_id", "variant_id",
        "target_p", "target_beta", "target_se", "target_q", "target_af", "target_ma_samples", "target_ma_count",
        "other_p", "other_beta", "other_se", "other_q", "other_af", "other_ma_samples", "other_ma_count",
        "prenatal_beta", "prenatal_se", "postnatal_beta", "postnatal_se",
        "delta_beta_pre_minus_post", "z_effect_difference", "p_effect_difference", "p_effect_difference_fdr",
        "strong_target", "other_near_zero", "target_phenotype_significant",
        "other_phenotype_nonsignificant", "effect_difference_significant",
        "region", "strand", "proximal_pas", "common_utr_bp", "long_extension_bp",
        "figure_extension_sufficient", "allele_phenotype_key"
    )
    fwrite(screen[, ..keep_columns], output_path, sep = "\t")

}

# --- shortlist events and verify genotype trends with raw DaPars2 PDUI ---

strict_events <- screen[
    evidence_class == "strict_stage_specific" & figure_extension_sufficient == TRUE
]
strict_events[, min_cross_stage_maf := pmin(
    target_af, 1 - target_af, other_af, 1 - other_af
)]
strict_events[, qtl_visual_score :=
    abs(target_beta) - abs(other_beta) +
    1.5 * min_cross_stage_maf +
    0.1 * pmin(log10(long_extension_bp / figure_extension_min_bp), 1)
]
event_best <- strict_events[
    order(-qtl_visual_score, p_effect_difference_fdr),
    .SD[1L],
    by = .(target_stage, phenotype_id)
]
shortlist <- event_best[
    order(target_stage, -qtl_visual_score),
    head(.SD, shortlist_events_per_stage),
    by = target_stage
]

prenatal_samples <- fread(prenatal_samples_path, header = FALSE)[[1]]
postnatal_samples <- fread(postnatal_samples_path, header = FALSE)[[1]]
raw_pdui_summary <- rbindlist(lapply(seq_len(nrow(shortlist)), function(i) {
    candidate <- shortlist[i]
    rbindlist(list(
        summarize_raw_pdui(candidate, "prenatal", prenatal_samples, prenatal_pfile),
        summarize_raw_pdui(candidate, "postnatal", postnatal_samples, postnatal_pfile)
    ))
}))
setcolorder(raw_pdui_summary, c(
    "target_stage", "gene", "transcript", "phenotype_id", "variant_id",
    "measured_stage", "dosage", "n", "median_raw_pdui", "mean_raw_pdui",
    "sd_raw_pdui", "raw_pdui_slope", "raw_pdui_p"
))
fwrite(raw_pdui_summary, raw_pdui_summary_path, sep = "\t")

raw_metrics <- raw_pdui_summary[, .(
    raw_pdui_slope = first(raw_pdui_slope),
    raw_pdui_p = first(raw_pdui_p),
    raw_pdui_median_dosage0 = value_at_dosage(dosage, median_raw_pdui, 0L),
    raw_pdui_median_dosage2 = value_at_dosage(dosage, median_raw_pdui, 2L),
    raw_pdui_median_delta_2_minus_0 =
        value_at_dosage(dosage, median_raw_pdui, 2L) -
        value_at_dosage(dosage, median_raw_pdui, 0L),
    raw_pdui_min_genotype_n = min(n),
    raw_pdui_genotype_groups = .N
), by = .(target_stage, phenotype_id, variant_id, measured_stage)]
raw_metrics_wide <- dcast(
    raw_metrics,
    target_stage + phenotype_id + variant_id ~ measured_stage,
    value.var = c(
        "raw_pdui_slope", "raw_pdui_p", "raw_pdui_median_dosage0",
        "raw_pdui_median_dosage2", "raw_pdui_median_delta_2_minus_0",
        "raw_pdui_min_genotype_n", "raw_pdui_genotype_groups"
    )
)
shortlist <- merge(
    shortlist,
    raw_metrics_wide,
    by = c("target_stage", "phenotype_id", "variant_id"),
    all.x = TRUE
)
shortlist_keep <- c(
    "target_stage", "gene", "transcript", "phenotype_id", "variant_id",
    "target_beta", "target_p", "target_q", "other_beta", "other_p", "other_q",
    "p_effect_difference_fdr", "target_af", "other_af", "min_cross_stage_maf",
    "target_ma_samples", "other_ma_samples", "long_extension_bp", "qtl_visual_score",
    "raw_pdui_slope_prenatal", "raw_pdui_p_prenatal",
    "raw_pdui_median_delta_2_minus_0_prenatal", "raw_pdui_min_genotype_n_prenatal",
    "raw_pdui_slope_postnatal", "raw_pdui_p_postnatal",
    "raw_pdui_median_delta_2_minus_0_postnatal", "raw_pdui_min_genotype_n_postnatal"
)
setorder(shortlist, target_stage, -qtl_visual_score)
fwrite(shortlist[, ..shortlist_keep], shortlist_path, sep = "\t")

message("Saved: ", output_path)
message("Saved: ", shortlist_path)
message("Saved: ", raw_pdui_summary_path)
print(screen[evidence_class == "strict_stage_specific",
    .N, by = .(target_stage, figure_extension_sufficient)])
print(screen[evidence_class == "strict_stage_specific",
    .(target_stage, gene, variant_id, target_beta, target_q, other_beta, other_p,
      other_q, p_effect_difference_fdr, long_extension_bp)][1:min(.N, 20L)])