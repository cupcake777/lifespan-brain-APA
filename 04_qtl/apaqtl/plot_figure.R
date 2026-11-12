# Figures 1-4 of the stage-specific apaQTL analysis.
# Reads QTL tables from data/04_qtl and writes figures to figures/04_qtl.
suppressMessages({
    library(arrow)
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(qvalue)
    library(scales)
    library(ComplexHeatmap)
    library(rtracklayer)
    library(GenomicRanges)
    library(parallel)
    library(ggrepel)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
data_dir <- file.path(apa_root, "data", "04_qtl")
fig_dir <- file.path(apa_root, "figures", "04_qtl")
hpc_qtl <- file.path(Sys.getenv("APA_HPC_HOME"), "project/work/apa/01_QTL")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1L)
inputs <- list(
    prenatal_sig      = file.path(data_dir, "prenatal.sig_QTL.txt.gz"),
    postnatal_sig     = file.path(data_dir, "postnatal.sig_QTL.txt.gz"),
    prenatal_nominal  = file.path(data_dir, "prenatal.cis_qtl.txt.gz"),
    postnatal_nominal = file.path(data_dir, "postnatal.cis_qtl.txt.gz"),
    pas_info          = file.path(data_dir, "pre_input", "proximal.apa.info"),
    brain_sig         = file.path(data_dir, "eur_adult.sig_pairs.txt.gz"),
    brain_nominal     = file.path(hpc_qtl, "BrainGvex", "nominal", "eur_adult.cis_qtl_pairs.txt.gz"),
    gtex_sig          = file.path(Sys.getenv("APA_GROUP_SHARE"), ".ref/GTEx_v11/GTEx_Analysis_v11_apaQTL/Brain_Frontal_Cortex_BA9.v11.apaQTLs.signif_pairs.parquet"),
    gtex_nominal      = file.path(apa_root, "ref/gtex/Brain_Frontal_Cortex_BA9.v11.cis_apaqtl.allpairs.txt.gz"),
    prenatal_susie    = file.path(hpc_qtl, "susie", "prenatal.SuSiE_summary.parquet"),
    postnatal_susie   = file.path(hpc_qtl, "susie", "postnatal.SuSiE_summary.parquet"),
    gtf               = file.path(apa_root, "ref/gencode.v40.annotation.gtf.gz"),
    prenatal_bw_dir   = file.path(Sys.getenv("APA_GROUP_SHARE"), "data/CHFB/bigwig"),
    postnatal_bam_root = file.path(Sys.getenv("APA_GROUP_HOME"), "guomh/RNA/bam"),
    prenatal_pfile    = file.path(data_dir, "pre_input", "genotype.prenatal"),
    postnatal_pfile   = file.path(data_dir, "pre_input", "genotype.postnatal"),
    prenatal_bed      = file.path(data_dir, "pre_input", "phenotype.Prenatal.bed.gz"),
    postnatal_bed     = file.path(data_dir, "pre_input", "phenotype.Postnatal.bed.gz"),
    prenatal_samples  = file.path(data_dir, "pre_input", "samples.prenatal.txt"),
    postnatal_samples = file.path(data_dir, "pre_input", "samples.postnatal.txt"),
    ds_glob           = file.path(data_dir, "seed_*_post.cis_qtl.txt.gz"),
    ds_susie_glob     = file.path(hpc_qtl, "downsample", "seed_*", "qtl", "Postnatal", "susie", "*postnatal.SuSiE_summary.parquet")
)
figure2_pi1_cache <- file.path(data_dir, "figure2_pair_pi1.tsv")


# --- preprocess annotation ---
gtf_data <- rtracklayer::import(inputs$gtf)
gene_idx <- gtf_data$type == "gene"
gtf_dt <- data.table::as.data.table(mcols(gtf_data)[gene_idx, c("gene_id", "gene_name")])
gtf_dt[, gene_id_clean := sub("\\.[0-9]+$", "", as.character(gene_id))]
id2name <- setNames(as.character(gtf_dt$gene_name), gtf_dt$gene_id_clean)

stage_levels  <- c("prenatal", "postnatal")
set_levels    <- c("prenatal", "postnatal", "BrainGvex", "GTEx Brain Cortex")
ds_q_threshold <- 0.05
dpi <- 320
base_size <- 12

stage_colors <- c(prenatal = "#E69F00", postnatal = "#0072B2")
ink      <- "#222222"
grid_col <- "#E5E5E5"
line_col <- "#9E9E9E"

fig_dims <- list(
    figure1 = c(12.8, 7.2),
    figure2 = c(12.0, 5.8),
    figure3 = c(10.8, 8.4),
    figure4 = c(13.2, 8.0)
)
figure3c_min_block_width_mb <- 0.55
figure3c_top_overlap_label_n <- 12L
figure3c_max_labels_per_chr <- 1L
figure3c_row_spacing <- 1.35
figure3c_left_label_x <- -15
figure3c_right_label_x <- 265
figure4_tracks_per_genotype <- 1L


# --- helpers ---

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

make_phenotype_key <- function(phenotype_id) {
    phenotype_id <- as.character(phenotype_id)
    vapply(strsplit(phenotype_id, "|", fixed = TRUE), function(parts) {
        if (length(parts) >= 4L && all(!is.na(parts[c(1L, 3L, 4L)]))) {
            paste(parts[1L], parts[3L], parts[4L], sep = "|")
        } else {
            paste(parts, collapse = "|")
        }
    }, character(1))
}


add_pair_keys <- function(dt, gene_col = NULL) {
    dt <- data.table::as.data.table(data.table::copy(dt))
    target_col <- c(gene_col, "gene_id", "group_id", "phenotype_id")[c(gene_col, "gene_id", "group_id", "phenotype_id") %in% names(dt)][1]
    if (is.na(target_col)) stop("Required gene or phenotype column not found.")
    raw_genes <- as.character(dt[[target_col]])
    if (target_col == "phenotype_id") {
        parts <- data.table::tstrsplit(raw_genes, "|", fixed = TRUE, fill = NA_character_)
        out <- parts[[2]]
        bad <- is.na(out) | out == ""
        out[bad] <- raw_genes[bad]
        raw_genes <- out
    }
    cleaned_genes <- sub("_.*$", "", sub("\\.[0-9]+$", "", raw_genes))
    is_ensg <- startsWith(cleaned_genes, "ENSG")
    translated_names <- id2name[cleaned_genes]
    dt[, gene_key := fifelse(is_ensg & !is.na(translated_names), translated_names, cleaned_genes)]
    if (!"variant_id" %in% names(dt)) stop("column 'variant_id' not found.")
    vk <- sub("_b[0-9]+$", "", as.character(dt$variant_id))
    vm <- regmatches(vk, regexec("^(chr[^:_]+)[_:]([0-9]+)[_:]([^_:]+)[_:]([^_:]+)$", vk))
    ok <- lengths(vm) == 5
    if (any(ok)) vk[ok] <- vapply(vm[ok], \(z) paste(z[2], z[3], z[4], z[5], sep = ":"), character(1))
    dt[, variant_key := vk]
    dt[, pair_key    := paste(gene_key, variant_key, sep = "||")]
    if ("phenotype_id" %in% names(dt)) {
        dt[, phenotype_key := make_phenotype_key(phenotype_id)]
        dt[, allele_phenotype_key := paste(phenotype_key, variant_key, sep = "||")]
    }
    return(dt[])
}


format_p <- function(p) {
    if (!is.finite(p)) return("NA")
    if (p == 0) return("<2.2e-308")
    if (p < 1e-3) return(sprintf("%.1e", p))
    sprintf("%.3f", p)
}

save_plot <- function(plot, name, width, height) {
    ggplot2::ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
    ggplot2::ggsave(file.path(fig_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = dpi)
}

theme_apa <- function() {
    ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
            axis.text        = ggplot2::element_text(color = "black"),
            axis.title       = ggplot2::element_text(color = ink),
            legend.position  = "bottom",
            legend.title     = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = grid_col, linewidth = 0.25),
            panel.grid.minor = ggplot2::element_blank(),
            plot.tag         = ggplot2::element_text(face = "bold", size = base_size + 1),
            plot.margin      = ggplot2::margin(8, 8, 8, 8)
        )
}
ggplot2::theme_set(theme_apa())

cs_intervals <- function(cs_dt, stage_label) {
    parsed <- tstrsplit(cs_dt$variant_key, ":", fixed = TRUE)
    cs_dt <- cbind(cs_dt, chr = parsed[[1]], pos = as.numeric(parsed[[2]]))
    cs_dt <- cs_dt[is.finite(pos)]
    cs_dt[, .(
        start = min(pos), end = max(pos),
        center = mean(pos), max_pip = max(pip, na.rm = TRUE),
        n_variants = .N, stage = stage_label
    ), by = .(gene_key, cs_id, chr)]
}

stream_nominal_allele_phenotype_matches <- function(path, key_dt) {
    key_path <- tempfile(pattern = "figure4_allele_phenotype_keys_", fileext = ".tsv")
    fwrite(unique(key_dt[, .(allele_phenotype_key)]), key_path, sep = "\t", col.names = FALSE)
    on.exit(unlink(key_path), add = TRUE)

    awk_script <- paste0(
        "awk 'BEGIN{FS=OFS=\"\\t\"} ",
        "FNR==NR{k[$1]=1; next} ",
        "FNR==1{for(i=1;i<=NF;i++){",
        "if($i==\"phenotype_id\") pid=i; ",
        "if($i==\"variant_id\") vid=i; ",
        "if($i==\"pval_nominal\") pval=i; ",
        "if($i==\"slope\") slope=i}; ",
        "if(!pid || !vid || !pval || !slope){exit 2}; ",
        "print \"allele_phenotype_key\",\"other_stage_pval_nominal\",\"other_stage_slope\"; next} ",
        "{key=phenotype_key($pid) \"||\" variant_key($vid); if(key in k){print key,$pval,$slope}} ",
        "function phenotype_key(pid, parts){",
        "split(pid,parts,\"\\\\|\"); ",
        "if(length(parts)>=4){return parts[1] \"|\" parts[3] \"|\" parts[4]} else {return pid}} ",
        "function variant_key(v, x, parts){",
        "x=v; sub(/_b[0-9]+$/, \"\", x); gsub(/_/, \":\", x); split(x,parts,\":\"); ",
        "if(length(parts)>=4){return parts[1] \":\" parts[2] \":\" parts[3] \":\" parts[4]} else {return x}}' ",
        shQuote(key_path), " -"
    )

    matched <- fread(cmd = paste("gzip -cd", shQuote(path), "|", awk_script))
    if (nrow(matched) == 0L) {
        return(data.table(
            allele_phenotype_key = character(),
            other_stage_pval_nominal = numeric(),
            other_stage_slope = numeric()
        ))
    }
    matched[, `:=`(
        other_stage_pval_nominal = as.numeric(other_stage_pval_nominal),
        other_stage_slope = as.numeric(other_stage_slope)
    )]
    matched <- matched[is.finite(other_stage_pval_nominal)]
    setorder(matched, other_stage_pval_nominal)
    unique(matched, by = "allele_phenotype_key")
}


enrich_locus <- function(spec_dt) {
    if (nrow(spec_dt) == 0) return(spec_dt)
    m <- merge(spec_dt, pas_info[, .(phenotype_id, pas_pos)], by = "phenotype_id", all.x = TRUE)
    parts <- tstrsplit(m$phenotype_id, "|", fixed = TRUE)
    loc   <- tstrsplit(parts[[3]], ":", fixed = TRUE)
    m[, `:=`(gene = parts[[2]], chr = loc[[1]], region = parts[[3]])]
    m
}

make_ex <- function(row, window_kb = 12L, stage_specific = NA_character_) {
    list(phenotype_id = row$phenotype_id,
         phenotype_key = row$phenotype_key,
         gene = row$gene, chr = row$chr, pas_pos = row$pas_pos,
         region = row$region, lead_variant = row$variant_id,
         window_kb = window_kb,
         effect = row$effect,
         stage_specific = stage_specific)
}

pull_pdui_values <- function(bed, pid, samples) {
    row_idx <- which(bed$pid == pid)
    if (length(row_idx) == 0L) row_idx <- which(bed$phenotype_key == make_phenotype_key(pid))
    if (length(row_idx) == 0) return(data.table(sample_id = samples, pdui = NA_real_))
    vals <- bed[row_idx[1], ..samples]
    data.table(sample_id = samples, pdui = as.numeric(unlist(vals)))
}

score_pdui_trend <- function(pdui_dt, dosage) {
    x <- merge(pdui_dt, data.table(sample_id = names(dosage), dosage = as.integer(dosage)),
               by = "sample_id")
    x <- x[is.finite(pdui) & !is.na(dosage)]
    med <- x[, .(median_pdui = median(pdui, na.rm = TRUE), n = .N), by = dosage]
    if (nrow(med) < 3L || any(med$n < 3L)) return(data.table(trend_score = -Inf, trend_direction = "weak", trend_range = NA_real_))
    setorder(med, dosage)
    diffs <- diff(med$median_pdui)
    monotone <- all(diffs > 0) || all(diffs < 0)
    trend_range <- abs(med[dosage == 2, median_pdui] - med[dosage == 0, median_pdui])
    data.table(
        trend_score = ifelse(monotone, trend_range, trend_range * 0.25),
        trend_direction = ifelse(all(diffs > 0), "increasing", ifelse(all(diffs < 0), "decreasing", "weak")),
        trend_range = trend_range
    )
}


pick_star_example <- function(candidates, specific_stage, max_candidates = 24L) {
    candidates <- candidates[is.finite(pas_pos) & !is.na(chr) & !is.na(variant_id)]
    if (nrow(candidates) == 0L) stop("No candidate rows available for ", specific_stage)
    candidates <- candidates[order(-abs_effect)][seq_len(min(nrow(candidates), max_candidates))]
    scored <- rbindlist(lapply(seq_len(nrow(candidates)), function(i) {
        row <- candidates[i]
        row_abs_effect <- row$abs_effect[1L]
        vi <- tstrsplit(row$variant_id, ":", fixed = TRUE)
        v_chr <- vi[[1]]
        v_pos <- as.integer(vi[[2]])
        v_ref <- vi[[3]]
        v_alt <- vi[[4]]
        pre_dos <- extract_dosage(inputs$prenatal_pfile, v_chr, v_pos, v_ref, v_alt)
        post_dos <- extract_dosage(inputs$postnatal_pfile, v_chr, v_pos, v_ref, v_alt)
        pre_trend <- score_pdui_trend(
            pull_pdui_values(prenatal_bed_full, row$phenotype_id, prenatal_sample_list),
            pre_dos
        )
        post_trend <- score_pdui_trend(
            pull_pdui_values(postnatal_bed_full, row$phenotype_id, postnatal_sample_list),
            post_dos
        )
        trend <- if (specific_stage == "prenatal") pre_trend else post_trend
        nonspecific_trend <- if (specific_stage == "prenatal") post_trend else pre_trend
        bw_avail <- sum(intersect(names(pre_dos), names(bw_files_by_sample)) %in% prenatal_sample_list)
        bam_avail <- sum(intersect(names(post_dos), names(bam_files_by_sample)) %in% postnatal_sample_list)
        data.table(
            row,
            trend_score = trend$trend_score,
            trend_direction = trend$trend_direction,
            trend_range = trend$trend_range,
            nonspecific_trend_score = nonspecific_trend$trend_score,
            pre_samples_with_dosage = bw_avail,
            post_samples_with_dosage = bam_avail,
            star_score = row_abs_effect + 3 * pmax(trend$trend_score, 0) +
                0.03 * min(bw_avail, 12L) + 0.03 * min(bam_avail, 12L)
        )
    }), fill = TRUE)
    scored <- scored[is.finite(star_score)]
    fwrite(scored[order(-star_score)], file.path(data_dir, sprintf("figure4_%s_star_candidates.tsv", specific_stage)), sep = "\t")
    scored[order(-star_score)][1L]
}

extract_dosage <- function(pfile_prefix, chr, pos, ref, alt) {
    tmp_bed <- tempfile(fileext = ".bed1")
    tmp_out <- tempfile()
    writeLines(sprintf("%s\t%d\t%d", chr, pos, pos), tmp_bed)
    system(sprintf("plink2 --pfile %s --extract bed1 %s --export vcf --out %s --threads 2",
        pfile_prefix, shQuote(tmp_bed), tmp_out), ignore.stdout = TRUE, ignore.stderr = TRUE)
    on.exit(unlink(c(tmp_bed, paste0(tmp_out, c(".vcf", ".log")))), add = TRUE)
    vcf_path <- paste0(tmp_out, ".vcf")
    if (!file.exists(vcf_path)) return(setNames(integer(0), character(0)))
    vcf <- fread(vcf_path, skip = "#CHROM", colClasses = "character")
    sample_cols <- setdiff(names(vcf), c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"))
    dosage <- vapply(sample_cols, function(s) {
        gt <- vcf[[s]][1]
        if (is.na(gt) || gt %in% c("", "./.", ".")) return(NA_integer_)
        sum(as.integer(strsplit(gt, "[|/]")[[1]]))
    }, integer(1))
    setNames(dosage, sample_cols)
}

pick_by_geno <- function(pool, dosage, max_per = figure4_tracks_per_genotype) {
    ds <- dosage[intersect(names(dosage), pool)]
    ds <- ds[!is.na(ds)]
    unlist(lapply(list(`ref/ref` = names(ds)[ds == 0],
                        `ref/alt` = names(ds)[ds == 1],
                        `alt/alt` = names(ds)[ds == 2]),
        function(s) if (length(s) <= max_per) s else sample(s, max_per)), use.names = FALSE)
}

geno_label <- function(dosage) {
    fifelse(is.na(dosage), "ref/alt",
    fifelse(dosage == 0, "ref/ref",
    fifelse(dosage == 2, "alt/alt", "ref/alt")))
}

geno_label_allele <- function(genotype, ref, alt) {
    fcase(
        genotype == "ref/ref", paste(ref, ref, sep = "/"),
        genotype == "ref/alt", paste(ref, alt, sep = "/"),
        genotype == "alt/alt", paste(alt, alt, sep = "/"),
        default = genotype
    )
}


build_cov <- function(file_paths, sample_ids, source, chr, s_pos, e_pos, bin_sz) {
    rbindlist(lapply(seq_along(file_paths), function(i) {
        if (source == "bw") {
            cov <- rtracklayer::import.bw(file_paths[i], which = GRanges(chr, IRanges(s_pos, e_pos)))
            raw <- if (length(cov) == 0) {
                data.table(start = seq.int(s_pos, e_pos), end = seq.int(s_pos, e_pos), score = 0)
            } else {
                data.table(start = start(cov), end = end(cov), score = as.numeric(cov$score))
            }
        } else {
            dt <- tryCatch(
                fread(cmd = sprintf("samtools depth -a -r %s:%d-%d %s", chr, s_pos, e_pos, shQuote(file_paths[i])),
                      header = FALSE, col.names = c("chr", "pos", "depth")),
                error = function(e) data.table(chr = character(), pos = integer(), depth = integer())
            )
            raw <- data.table(start = integer(), end = integer(), score = numeric())
            if (nrow(dt) > 0) {
                dt[, bin := floor((pos - s_pos) / bin_sz) * bin_sz + s_pos]
                raw <- dt[, .(start = bin, end = bin + bin_sz - 1L, score = mean(depth)), by = bin][, bin := NULL]
            }
        }
        if (nrow(raw) == 0) return(NULL)
        bins <- data.table(start = seq(s_pos, e_pos, by = bin_sz))
        bins[, `:=`(end = start + bin_sz - 1L, bin_mid = start + bin_sz / 2)]
        raw[, bin := floor((start - s_pos) / bin_sz) * bin_sz + s_pos]
        agg <- raw[, .(score = sum(score)), by = .(start = bin)]
        g <- merge(bins, agg, by = "start", all.x = TRUE)
        g[is.na(score), score := 0]
        g[, .(bin_mid, score, sample_id = sample_ids[i])]
    }), fill = TRUE)
}


attach_geno <- function(cov_dt, dos) {
    if (is.null(cov_dt) || nrow(cov_dt) == 0) return(NULL)
    gmap <- data.table(sample_id = names(dos), dosage = as.integer(dos))
    gmap[, genotype := geno_label(dosage)]
    merge(cov_dt, gmap[, .(sample_id, genotype)], by = "sample_id", all.x = TRUE)
}


order_samples <- function(cov_g, v_ref, v_alt) {
    if (is.null(cov_g)) return(NULL)
    cov_g[, genotype := factor(genotype, levels = c("ref/ref", "ref/alt", "alt/alt"))]
    cov_g[, genotype_label := factor(geno_label_allele(as.character(genotype), v_ref, v_alt),
        levels = geno_label_allele(c("ref/ref", "ref/alt", "alt/alt"), v_ref, v_alt))]
    cov_g[, sample_max := max(score, na.rm = TRUE), by = sample_id]
    sample_order <- cov_g[, .(genotype = genotype[1L], sample_max = sample_max[1L]), by = sample_id][
        order(genotype, -sample_max), sample_id]
    display_levels <- sample_order
    cov_g[, sample_label := as.character(genotype_label)]
    label_levels <- cov_g[match(display_levels, sample_id), sample_label]
    cov_g[, sample_label := factor(sample_label, levels = label_levels)]
    cov_g
}


build_track <- function(cov_g, stage_label, ref_value, pas_pos, v_pos, s_pos, e_pos, geno_colors) {
    if (is.null(cov_g) || nrow(cov_g) == 0) return(ggplot() + theme_void() + labs(title = stage_label))
    cov_g[, score_plot := pmin(score / ref_value, 1.35)]
    ggplot(cov_g, aes(bin_mid, score_plot)) +
        geom_area(aes(fill = genotype_label), color = NA, alpha = 0.82) +
        geom_line(color = "grey20", linewidth = 0.16) +
        geom_vline(xintercept = pas_pos, color = "#D55E00", linetype = "dashed", linewidth = 0.45) +
        geom_vline(xintercept = v_pos, color = "#CC79A7", linetype = "dotted", linewidth = 0.45) +
        facet_wrap(~sample_label, ncol = 1) +
        scale_fill_manual(values = geno_colors, drop = FALSE) +
        scale_y_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
        scale_x_continuous(expand = expansion(mult = c(0, 0)), labels = comma) +
        coord_cartesian(xlim = c(s_pos, e_pos), clip = "on") +
        labs(x = NULL, y = stage_label, fill = "Genotype") +
        theme(
            strip.text = element_text(size = 7.2, color = "grey25", hjust = 0.5, face = "bold"),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_line(color = grid_col, linewidth = 0.18),
            panel.grid.minor = element_blank(),
            panel.spacing.y = unit(0.34, "lines"),
            plot.margin = margin(2, 6, 3, 4),
            legend.position = "none"
        )
}


build_example_figure <- function(ex, bw_files, bam_files,
                                  prenatal_pfile, postnatal_pfile,
                                  prenatal_samples, postnatal_samples,
                                  prenatal_bed, postnatal_bed, pas_info) {

    vi      <- tstrsplit(ex$lead_variant, ":", fixed = TRUE)
    v_chr   <- vi[[1]]; v_pos <- as.integer(vi[[2]])
    v_ref   <- vi[[3]]; v_alt <- vi[[4]]

    loc_parts <- tstrsplit(sub("^chr[^:]+:", "", ex$region), "-", fixed = TRUE)
    gene_start <- suppressWarnings(as.integer(loc_parts[[1]]))
    gene_end   <- suppressWarnings(as.integer(loc_parts[[2]]))
    locus_min <- min(ex$pas_pos, v_pos, gene_start, gene_end, na.rm = TRUE)
    locus_max <- max(ex$pas_pos, v_pos, gene_start, gene_end, na.rm = TRUE)
    pad <- max(2500L, ceiling((locus_max - locus_min) * 0.20), ex$window_kb * 1000L)
    s_pos  <- max(1L, as.integer(locus_min - pad))
    e_pos  <- as.integer(locus_max + pad)
    bin_sz <- 50L

    bw_dos  <- extract_dosage(prenatal_pfile,  v_chr, v_pos, v_ref, v_alt)
    bam_dos <- extract_dosage(postnatal_pfile, v_chr, v_pos, v_ref, v_alt)

    bw_picks  <- pick_by_geno(union(intersect(prenatal_samples, names(bw_files)), names(bw_files)),  bw_dos)
    bam_picks <- pick_by_geno(union(intersect(postnatal_samples, names(bam_files)), names(bam_files)), bam_dos)
    bw_picks  <- bw_picks[bw_picks %in% names(bw_files)]
    bam_picks <- bam_picks[bam_picks %in% names(bam_files)]

    bw_cov  <- if (length(bw_picks)  > 0) build_cov(bw_files[bw_picks],   bw_picks,  "bw", ex$chr, s_pos, e_pos, bin_sz) else NULL
    bam_cov <- if (length(bam_picks) > 0) build_cov(bam_files[bam_picks], bam_picks, "bam", ex$chr, s_pos, e_pos, bin_sz) else NULL
    bw_cov_g  <- attach_geno(bw_cov,  bw_dos)
    bam_cov_g <- attach_geno(bam_cov, bam_dos)
    bw_cov_g  <- order_samples(bw_cov_g, v_ref, v_alt)
    bam_cov_g <- order_samples(bam_cov_g, v_ref, v_alt)

    geno_labels <- geno_label_allele(c("ref/ref", "ref/alt", "alt/alt"), v_ref, v_alt)
    geno_colors <- setNames(c(stage_colors[["prenatal"]], "#BDBDBD", "#1A1A1A"), geno_labels)

    stage_specific <- if (grepl("prenatal", ex$stage_specific %||% "", fixed = TRUE)) "prenatal" else
        if (grepl("postnatal", ex$stage_specific %||% "", fixed = TRUE)) "postnatal" else
        ifelse(abs(ex$effect) >= 0, "prenatal", "postnatal")
    locus_title <- sprintf("%s-specific: %s", tools::toTitleCase(stage_specific), ex$gene)
    locus_subtitle <- sprintf("%s; lead %s; effect %.2f", ex$region, ex$lead_variant, ex$effect)

    pre_ref <- if (!is.null(bw_cov_g) && nrow(bw_cov_g) > 0) quantile(bw_cov_g$score, 0.995, na.rm = TRUE) else 1
    post_ref <- if (!is.null(bam_cov_g) && nrow(bam_cov_g) > 0) quantile(bam_cov_g$score, 0.995, na.rm = TRUE) else 1
    if (!is.finite(pre_ref) || pre_ref <= 0) pre_ref <- 1
    if (!is.finite(post_ref) || post_ref <= 0) post_ref <- 1


    p_pre <- build_track(bw_cov_g, "prenatal", pre_ref, ex$pas_pos, v_pos, s_pos, e_pos, geno_colors) +
        labs(title = locus_title, subtitle = locus_subtitle) +
        theme(
            plot.title = element_text(size = 12, face = "bold", color = "grey15"),
            plot.subtitle = element_text(size = 8.2, color = "grey35")
        )
    p_post <- build_track(bam_cov_g, "postnatal", post_ref, ex$pas_pos, v_pos, s_pos, e_pos, geno_colors)

    gtf <- rtracklayer::import(inputs$gtf, format = "gtf")
    mc   <- as.data.frame(mcols(gtf))
    hits <- gtf[seqnames(gtf) == ex$chr & start(gtf) <= e_pos & end(gtf) >= s_pos &
                mc$gene_name == ex$gene &
                mc$type %in% c("exon", "CDS", "UTR", "transcript", "gene",
                               "three_prime_UTR", "five_prime_UTR")]
    mc2  <- as.data.frame(mcols(hits))
    gene_dt <- data.table(chr = as.character(seqnames(hits)), start = start(hits), end = end(hits),
               strand = as.character(strand(hits)), type = mc2$type,
               transcript_id = mc2$transcript_id, gene_name = mc2$gene_name)
    if (nrow(gene_dt) > 0) {
        tx_rank <- gene_dt[type == "exon", .(n_exon = .N, total_bp = sum(end - start + 1)), by = transcript_id]
        gene_dt <- gene_dt[transcript_id %in% tx_rank[order(-n_exon, -total_bp), transcript_id][1L] | type == "gene"]
    }
    p_gene <- if (nrow(gene_dt) > 0) {
        feat  <- gene_dt[type %in% c("exon", "CDS", "three_prime_UTR", "five_prime_UTR")]
        back  <- feat[type %in% c("exon", "CDS"), .(xstart = min(start), xend = max(end)), by = transcript_id]
        ggplot(feat) +
            geom_segment(data = back, aes(x = xstart, xend = xend, y = 1, yend = 1),
                         linewidth = 0.24, color = line_col) +
            geom_rect(aes(xmin = start, xmax = end, ymin = 0.88, ymax = 1.12, fill = type), color = NA) +
            scale_fill_manual(values = c(exon = "#9E9E9E", CDS = "#222222",
                three_prime_UTR = "#E69F00", five_prime_UTR = "#56B4E9"), name = "Feature", drop = FALSE) +
            scale_y_continuous(limits = c(0.80, 1.20), expand = c(0, 0)) +
            geom_vline(xintercept = ex$pas_pos, color = stage_colors[["prenatal"]],
                       linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = v_pos, color = "#CC79A7", linetype = "dotted", linewidth = 0.5) +
            scale_x_continuous(expand = expansion(mult = c(0, 0)), labels = comma) +
            coord_cartesian(xlim = c(s_pos, e_pos), clip = "on") +
            labs(x = sprintf("Genomic position on %s (bp)", ex$chr), y = NULL) +
            theme(
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                plot.margin = margin(0, 6, 1, 4),
                legend.position = "none",
                legend.key.size = unit(0.25, "lines"),
                legend.text = element_text(size = 7),
                legend.title = element_text(size = 7)
            )
    } else ggplot() + theme_void()

    pid <- ex$phenotype_id %||% pas_info[grepl(ex$region, phenotype_id, fixed = TRUE), phenotype_id[1]]
    if (is.na(pid)) stop("Cannot find PAS phenotype row for region ", ex$region)

    pdui_pre  <- merge(pull_pdui_values(prenatal_bed,  pid, prenatal_samples),
                       data.table(sample_id = names(bw_dos),  dosage = as.integer(bw_dos)),
                       by = "sample_id", all.x = TRUE)
    pdui_post <- merge(pull_pdui_values(postnatal_bed, pid, postnatal_samples),
                       data.table(sample_id = names(bam_dos), dosage = as.integer(bam_dos)),
                       by = "sample_id", all.x = TRUE)
    pdui_long <- rbindlist(list(
        pdui_pre[,  .(sample_id, pdui, stage = "prenatal",  dosage)],
        pdui_post[, .(sample_id, pdui, stage = "postnatal", dosage)]
    ), fill = TRUE)
    pdui_long <- pdui_long[!is.na(dosage)]
    pdui_long[, genotype := factor(geno_label(dosage), levels = c("ref/ref", "ref/alt", "alt/alt"))]
    pdui_long[, genotype_label := factor(geno_label_allele(as.character(genotype), v_ref, v_alt),
        levels = geno_labels)]
    pdui_long[, stage    := factor(stage, levels = stage_levels)]
    pdui_long[, stage_focus := fifelse(as.character(stage) == stage_specific, "specific stage", "other stage")]

    pdui_medians <- pdui_long[, .(
        median_pdui = median(pdui, na.rm = TRUE),
        n = .N,
        stage_focus = stage_focus[1L]
    ), by = .(genotype_label, stage)]
    
    y_lim <- quantile(pdui_long$pdui, c(0.01, 0.99), na.rm = TRUE)
    y_pad <- diff(y_lim) * 0.10
    if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.2

    p_box <- ggplot(pdui_long, aes(genotype_label, pdui, fill = stage, alpha = stage_focus)) +
        geom_violin(trim = TRUE, scale = "width", linewidth = 0.25,
                    position = position_dodge(width = 0.72)) +
        geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.28,
                     position = position_dodge(width = 0.72)) +
        geom_crossbar(data = pdui_medians,
                      aes(x = genotype_label, y = median_pdui, ymin = median_pdui, ymax = median_pdui,
                          color = stage, alpha = stage_focus),
                      width = 0.28, linewidth = 0.42, middle.linewidth = 0.42,
                      position = position_dodge(width = 0.72), show.legend = FALSE) +
        scale_fill_manual(values = stage_colors, labels = c(prenatal = "pre", postnatal = "post")) +
        scale_color_manual(values = stage_colors, labels = c(prenatal = "pre", postnatal = "post")) +
        scale_alpha_manual(values = c(`specific stage` = 0.95, `other stage` = 0.28), guide = "none") +
        coord_cartesian(ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)) +
        labs(x = sprintf("Genotype at %s", ex$lead_variant),
             y = "APA phenotype score", fill = "Stage", color = "Stage") +
        theme(
            legend.position = "bottom",
            plot.margin = margin(2, 6, 6, 4),
            axis.text.x = element_text(size = 8.5),
            axis.title.x = element_text(size = 9)
        )

    (p_pre / p_post / p_gene / p_box) +
        plot_layout(heights = c(1.75, 1.75, 0.34, 1.35), guides = "collect") &
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 8)
        )
}


# --- load data ---
prenatal_sig  <- add_pair_keys(fread(inputs$prenatal_sig))[,  stage := "prenatal"]
postnatal_sig <- add_pair_keys(fread(inputs$postnatal_sig))[, stage := "postnatal"]
sig_pairs <- unique(rbindlist(list(prenatal_sig, postnatal_sig), fill = TRUE), by = c("stage", "pair_key"))
sig_pairs[, `:=`(stage = factor(stage, levels = stage_levels), variant_pos = suppressWarnings(as.integer(data.table::tstrsplit(variant_key, ":", fixed = TRUE)[[2]])))]
pas_info <- unique(fread(inputs$pas_info, select = c("phenotype_id", "Predicted_Proximal_APA")), by = "phenotype_id")
setnames(pas_info, "Predicted_Proximal_APA", "pas_pos")
pas_info[, pas_pos := as.integer(pas_pos)]

brain_sig <- add_pair_keys(fread(inputs$brain_sig))
gtex_sig <- add_pair_keys(as.data.table(read_parquet(inputs$gtex_sig)), gene_col = "group_id")
prenatal_susie  <- add_pair_keys(as.data.table(read_parquet(inputs$prenatal_susie)))
postnatal_susie <- add_pair_keys(as.data.table(read_parquet(inputs$postnatal_susie)))
ds_susie_paths <- Sys.glob(inputs$ds_susie_glob)
if (length(ds_susie_paths) == 0L) stop("No postnatal DS SuSiE summary files found: ", inputs$ds_susie_glob)
postnatal_ds_susie <- rbindlist(lapply(ds_susie_paths, function(path) {
    dt <- add_pair_keys(as.data.table(read_parquet(path)))
    dt[, seed := as.integer(sub(".*seed_([0-9]+).*", "\\1", path))]
    dt[, .(pip = max(pip, na.rm = TRUE)), by = .(gene_key, variant_key)]
}), fill = TRUE)
postnatal_ds_susie <- postnatal_ds_susie[, .(pip = max(pip, na.rm = TRUE)), by = .(gene_key, variant_key)]

# --- figure 1 data process ---
counts_long <- melt(
    sig_pairs[, .(genes = uniqueN(gene_key), variants = uniqueN(variant_key), pairs = uniqueN(pair_key)), by = stage],
    id.vars = "stage", variable.name = "metric", value.name = "count"
)
counts_long[, metric := factor(metric, levels = c("genes", "variants", "pairs"),
                                labels = c("apaGenes", "variants", "pairs"))]

signal_distance <- merge(sig_pairs[, .(stage, phenotype_id, variant_pos)], pas_info, by = "phenotype_id", all.x = TRUE)
signal_distance <- signal_distance[is.finite(variant_pos) & is.finite(pas_pos)]
signal_distance[, distance_kb := (variant_pos - pas_pos) / 1000]
signal_distance <- signal_distance[abs(distance_kb) <= 1000]

shared_beta <- dcast(sig_pairs[, .(stage, pair_key, slope)], pair_key ~ stage, value.var = "slope")
shared_beta <- shared_beta[is.finite(prenatal) & is.finite(postnatal)]
shared_beta[, `:=`(abs_prenatal = abs(prenatal), abs_postnatal = abs(postnatal))]
shared_beta[, direction := fifelse(sign(prenatal) == sign(postnatal), "Same direction", "Discordant")]

rho <- cor(shared_beta$prenatal, shared_beta$postnatal, method = "spearman", use = "complete.obs")
sign_concordance <- mean(shared_beta$direction == "Same direction")

beta_long <- melt(shared_beta[, .(pair_key, prenatal = abs_prenatal, postnatal = abs_postnatal)],
    id.vars = "pair_key", variable.name = "stage", value.name = "abs_beta")
beta_long[, stage := factor(stage, levels = stage_levels)]

# --- figure 1A: significant qtl counts ---

p1a <- ggplot(counts_long, aes(x = count, y = metric)) +
    geom_line(aes(group = metric), color = "#E5E7EB", linewidth = 1.5) +
    geom_point(aes(color = stage), size = 4.0) +
    geom_text(aes(label = comma(count), color = stage), vjust = -0.8, size = 4.0, fontface = "bold", show.legend = FALSE) +
    scale_color_manual(values = stage_colors) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                  expand = expansion(mult = c(0.1, 0.1))) + 
    labs(x = "Significant count", y = NULL) +
    theme_apa() +
    theme(
        panel.grid.major.x = element_line(color = "#F3F4F6", linewidth = 0.4), 
        panel.grid.major.y = element_blank(),
        legend.position = "top"
    )

# --- figure 1B: distance to proximal PAS ---

p1b <- ggplot(signal_distance, aes(distance_kb, color = stage, fill = stage)) +
    geom_vline(xintercept = 0, color = line_col, linewidth = 0.35, linetype = "dashed") +
    geom_density(alpha = 0.22, linewidth = 0.55, bounds = c(-1000, 1000)) +
    scale_color_manual(values = stage_colors) +
    scale_fill_manual(values = stage_colors) +
    scale_x_continuous(labels = comma) +
    coord_cartesian(xlim = c(-1000, 1000)) +
    labs(x = "Distance to proximal PAS (kb)", y = "Density")

# --- figure 1C: absolute beta by stage ---

beta_wide <- dcast(beta_long, pair_key ~ stage, value.var = "abs_beta")
beta_wide[, delta_beta := postnatal - prenatal]
median_delta_full <- median(beta_wide$delta_beta, na.rm = TRUE)
wilcox_full <- wilcox.test(beta_wide$postnatal, beta_wide$prenatal, paired = TRUE)


# --- Figure 1c: Paired Difference Density Plot (Full Postnatal vs Prenatal) ---
p1c <- ggplot(beta_wide, aes(x = delta_beta)) +
    geom_density(fill = stage_colors[["postnatal"]], alpha = 0.25, 
                 color = stage_colors[["postnatal"]], linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = median_delta_full, color = "#D55E00", linewidth = 0.6) +
    xlim(quantile(beta_wide$delta_beta, 0.005, na.rm = TRUE), 
         quantile(beta_wide$delta_beta, 0.995, na.rm = TRUE)) +
    annotate("text", 
             x = quantile(beta_wide$delta_beta, 0.01, na.rm = TRUE), 
             y = max(density(beta_wide$delta_beta, na.rm = TRUE)$y) * 0.9, 
             hjust = 0,
             label = sprintf("Paired Wilcoxon n = %s\np = %s\nMedian = %.4f", 
                             comma(nrow(beta_wide)), 
                             format_p(wilcox_full$p.value),
                             median_delta_full), size = 3.5, fontface = "bold") +
    labs(x = "Beta (|postnatal| - |prenatal|)", y = "Density") +
    theme_classic(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank())

# --- fig 1D: beta by stage (after ds) ---

discordant_genes <- shared_beta[direction == "Discordant"]
discordant_genes[, gene := sub("\\|\\|.*", "", pair_key)]
discordant_labels <- discordant_genes[, .SD[order(-(abs(prenatal) + abs(postnatal)))][1L], by = gene]

m <- max(abs(c(shared_beta$prenatal, shared_beta$postnatal)), na.rm = TRUE)

p1d <- ggplot() +
    geom_hline(yintercept = 0, color = grid_col, linewidth = 0.3) +
    geom_vline(xintercept = 0, color = grid_col, linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = line_col, linewidth = 0.35) +
    geom_point(data = shared_beta[direction == "Discordant"], 
               aes(prenatal, postnatal, color = direction), size = 1.8, alpha = 0.65) +
    geom_point(data = shared_beta[direction == "Same direction"], 
               aes(prenatal, postnatal, color = direction), size = 1.8, alpha = 0.55) +
    geom_text_repel(
        data = discordant_labels,
        aes(prenatal, postnatal, label = gene), 
        size = 3.0,
        fontface = "bold",
        color = "#CC79A7",
        max.overlaps = Inf,
        box.padding = 0.5,
        point.padding = 0.2,
        segment.size = 0.3,
        segment.color = "#CC79A7", 
        segment.alpha = 0.7,
        segment.linetype = "dashed", 
        min.segment.length = 0,
        show.legend = FALSE) +
    scale_color_manual(values = c(`Same direction` = "#666666", Discordant = "#CC79A7")) +
    coord_equal(xlim = c(-m, m), ylim = c(-m, m), clip = "off") +
    annotate("text", x = -m * 0.95, y = m * 0.93, hjust = 0,
        label = sprintf("rho=%.2f\nsame=%.1f%%", rho, sign_concordance * 100), size = 3.5, fontface = "bold") +
    labs(x = "prenatal beta", y = "postnatal beta") +
    theme_classic(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank(),
        legend.position = "none",
        plot.margin = margin(8, 18, 8, 8)
    )

fig1 <- (p1a + p1b) / (p1c + p1d) +
    plot_layout(heights = c(0.78, 1.35), widths = c(1, 1))

save_plot(fig1, "figure1", fig_dims$figure1[1], fig_dims$figure1[2])


# --- figure 2 data: cross-dataset overlap ---
sets_dt <- rbindlist(list(
    unique(sig_pairs[stage == "prenatal",  .(pair_key, dataset = "prenatal")]),
    unique(sig_pairs[stage == "postnatal", .(pair_key, dataset = "postnatal")]),
    unique(brain_sig[, .(pair_key, dataset = "BrainGvex")]),
    unique(gtex_sig[,  .(pair_key, dataset = "GTEx Brain Cortex")])
), fill = TRUE)

sets_wide <- dcast(sets_dt[, present := TRUE], pair_key ~ dataset, value.var = "present", fill = FALSE)
for (set_name in set_levels) {
    if (!set_name %in% names(sets_wide)) sets_wide[, (set_name) := FALSE]
}
sets_wide[, (set_levels) := lapply(.SD, as.logical), .SDcols = set_levels]
comb_m <- make_comb_mat(as.matrix(sets_wide[, ..set_levels]))

p2a_summary <- data.table(
    combination = comb_name(comb_m),
    degree = comb_degree(comb_m),
    pairs = comb_size(comb_m)
)
p2a_summary <- p2a_summary[order(-pairs, -degree)]
fwrite(p2a_summary, file.path(data_dir, "figure2a_intersection_summary.tsv"), sep = "\t")

p2a <- UpSet(
    comb_m,
    set_order  = set_levels,
    comb_order = order(comb_degree(comb_m), comb_size(comb_m), decreasing = c(FALSE, TRUE)),
    top_annotation = upset_top_annotation(
        comb_m, 
        ylim = c(0, max(comb_size(comb_m)) * 1.1),
        bar_width = 0.6, 
        gp = gpar(fill = "#4D4D4D"), 
        axis_param = list(
            side = "right", 
            gp = gpar(fontsize = 7)
        ),
        annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    right_annotation = upset_right_annotation(
        comb_m, 
        gp = gpar(fill = "#0072B2"),
        axis_param = list(
            gp = gpar(fontsize = 7)
        ),
        annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    row_names_side = "left")


# --- figure 2B: cross-dataset pi1 heatmap ---
discovery_sets <- list(
    prenatal          = unique(sig_pairs[stage == "prenatal",  .(pair_key)]),
    postnatal         = unique(sig_pairs[stage == "postnatal", .(pair_key)]),
    `BrainGvex` = unique(brain_sig[, .(pair_key)]),
    `GTEx Brain Cortex` = unique(gtex_sig[, .(pair_key)])
)
target_inputs <- list(
    prenatal            = list(path = inputs$prenatal_nominal),
    postnatal           = list(path = inputs$postnatal_nominal),
    `BrainGvex`         = list(path = inputs$brain_nominal),
    `GTEx Brain Cortex` = list(path = inputs$gtex_nominal)
)

pi1_pairs <- NULL
if (file.exists(figure2_pi1_cache)) {
    cached <- fread(figure2_pi1_cache)
    required_cols <- c("discovery", "target", "pi1", "pi0", "discovery_pairs", "matched_pairs", "target_input")
    if (all(required_cols %in% names(cached))) {
        cached <- cached[, ..required_cols]
        cached[, `:=`(discovery = as.character(discovery), target = as.character(target))]
        valid_labels <- all(unique(c(cached$discovery, cached$target)) %in% set_levels)
        expected <- CJ(discovery = set_levels, target = set_levels)[discovery != target]
        has_all_pairs <- nrow(merge(expected, cached, by = c("discovery", "target"))) == nrow(expected)
        target_paths <- vapply(target_inputs, \(x) x$path, character(1))
        valid_paths <- all(cached$target_input %in% target_paths)
        if (valid_labels && has_all_pairs && valid_paths) pi1_pairs <- cached[]
    }
}

if (is.null(pi1_pairs)) {
    n_cores <- min(4L, detectCores() - 1L)
    res_list <- mclapply(names(target_inputs), function(target_label) {
        tryCatch({
            key_dt <- rbindlist(lapply(setdiff(names(discovery_sets), target_label), function(disc_label) {
                discovery_sets[[disc_label]][, .(discovery = disc_label, pair_key)]
            }))
            if (nrow(key_dt) == 0) return(NULL)

            target <- target_inputs[[target_label]]
            matches <- data.table(pair_key = character(), pval_nominal = numeric(), slope = numeric())
            is_gtex <- grepl("GTEx_Analysis|v11", target$path)
            if (is_gtex) {
                raw_variants <- unique(tstrsplit(key_dt$pair_key, "||", fixed = TRUE)[[2]])
                gtex_variants <- paste0(gsub(":", "_", raw_variants), "_b38")
                vid_path <- tempfile(pattern = "gtex_vids_", fileext = ".tsv")
                fwrite(data.table(v = gtex_variants), vid_path, sep = "\t", col.names = FALSE)
        awk_gtex <- paste0(
            "awk 'BEGIN{FS=OFS=\"\\t\"} ",
            "FNR==NR{k[$1]=1; next} ",
            "FNR==1{for(i=1;i<=NF;i++){",
            "if($i==\"phenotype_id\") pid=i; ",
            "if($i==\"variant_id\") vid=i; ",
            "if($i==\"pval_nominal\") pval=i; ",
            "if($i==\"slope\") slope=i}; ",
            "if(!pid || !vid || !pval || !slope){exit 2}; ",
            "print \"phenotype_id\",\"variant_id\",\"pval_nominal\",\"slope\"; next} ",
            "$vid in k { ",
            "g=$pid; sub(/_.*/, \"\", g); sub(/\\.[0-9]+$/, \"\", g); ",
            "print g, $vid, $pval, $slope}' ",
            shQuote(vid_path), " -"
        )
                raw_dt <- fread(cmd = paste("gzip -cd", shQuote(target$path), "|", awk_gtex))
                unlink(vid_path)
                matched <- raw_dt[0]
                if (nrow(raw_dt) > 0) {
                    raw_dt[, gene_key := id2name[as.character(phenotype_id)]]
                    raw_dt[is.na(gene_key), gene_key := as.character(phenotype_id)]
                    raw_dt[, variant_key := gsub("_b[0-9]+$", "", as.character(variant_id))]
                    raw_dt[, variant_key := gsub("_", ":", variant_key)]
                    raw_dt[, pair_key := paste(gene_key, variant_key, sep = "||")]
                    matched <- raw_dt[pair_key %in% key_dt$pair_key]
                }
            } else {
                key_path <- tempfile(pattern = "pi1_keys_", fileext = ".tsv")
                fwrite(unique(key_dt[, .(pair_key)]), key_path, sep = "\t", col.names = FALSE)
        awk_script <- paste0(
            "awk 'BEGIN{FS=OFS=\"\\t\"} ",
            "FNR==NR{k[$1]=1; next} ",
            "FNR==1{for(i=1;i<=NF;i++){",
            "if($i==\"phenotype_id\") pid=i; ",
            "if($i==\"variant_id\") vid=i; ",
            "if($i==\"pval_nominal\") pval=i; ",
            "if($i==\"slope\") slope=i}; ",
            "if(!pid || !vid || !pval || !slope){exit 2}; ",
            "print \"pair_key\",\"pval_nominal\",\"slope\"; next} ",
            "{key=gene_key($pid) \"||\" variant_key($vid); if(key in k){print key,$pval,$slope}} ",
            "function gene_key(pid, parts, g){",
            "if(index(pid,\"|\")>0){split(pid,parts,\"\\\\|\"); g=parts[2]} else {split(pid,parts,\"_\"); g=parts[1]}; ",
            "sub(/\\.[0-9]+$/, \"\", g); sub(/_.*/, \"\", g); return g} ",
            "function variant_key(v, x, parts){",
            "x=v; sub(/_b[0-9]+$/, \"\", x); gsub(/_/, \":\", x); split(x,parts,\":\"); ",
            "if(length(parts)>=4){return parts[1] \":\" parts[2] \":\" parts[3] \":\" parts[4]} else {return x}}' ",
            shQuote(key_path), " -"
        )
                matched <- fread(cmd = paste("gzip -cd", shQuote(target$path), "|", awk_script))
                unlink(key_path)
            }
            if (nrow(matched) > 0) {
                matched[, pval_nominal := as.numeric(pval_nominal)]
                matched <- matched[is.finite(pval_nominal)]
                setorder(matched, pval_nominal)
                matches <- unique(matched, by = "pair_key")
            }

            merged <- merge(key_dt, matches, by = "pair_key", allow.cartesian = TRUE)
            res <- rbindlist(lapply(setdiff(names(discovery_sets), target_label), function(disc_label) {
                vals <- merged[discovery == disc_label & is.finite(pval_nominal)]

                if (nrow(vals) < 5) {
                    return(data.table(
                        discovery = disc_label, target = target_label,
                        pi1 = NA_real_, pi0 = NA_real_,
                        discovery_pairs = nrow(discovery_sets[[disc_label]]),
                        matched_pairs = nrow(vals),
                        target_input = target$path
                    ))
                }

                pi0_res <- qvalue::pi0est(vals$pval_nominal, lambda = seq(0.1, 0.9, 0.1), method = "smoother")
                data.table(
                    discovery = disc_label, target = target_label,
                    pi1 = 1 - pi0_res$pi0, pi0 = pi0_res$pi0,
                    discovery_pairs = nrow(discovery_sets[[disc_label]]),
                    matched_pairs = nrow(vals),
                    target_input = target$path
                )
            }))
            return(res)
        }, error = function(e) {
            cat(sprintf("[ERROR] %s failed: %s\n", target_label, e$message))
            return(NULL)
        })
    }, mc.cores = n_cores, mc.preschedule = FALSE)

    pi1_pairs <- rbindlist(res_list[!sapply(res_list, is.null)], fill = TRUE)
    fwrite(pi1_pairs, figure2_pi1_cache, sep = "\t")
}

setDTthreads(0)
pi1_heat <- merge(
    CJ(discovery = set_levels, target = set_levels)[discovery != target],
    pi1_pairs, by = c("discovery", "target"), all.x = TRUE
) 
pi1_heat[, `:=`(
    discovery = factor(discovery, levels = rev(set_levels)),
    target    = factor(target,    levels = set_levels)
)]

p2b <- ggplot(pi1_heat, aes(target, discovery, fill = pi1)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = ifelse(is.finite(pi1),
        sprintf("%.2f\nn=%s", pi1, comma(matched_pairs)), "")),
        size = 2.35, lineheight = 0.9) +
    scale_fill_viridis_c(limits = c(0, 1), option = "C", na.value = "#F2F2F2") +
    labs(x = "Target nominal dataset", y = "Discovery dataset", fill = "pi1") +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.position = "bottom",
        legend.key.width = unit(1.4, "cm"),
        plot.margin = margin(12, 18, 22, 12)
    )

fig2_upset <- wrap_elements(
    grid::grid.grabExpr({
        grid::pushViewport(grid::viewport(x = 0.54, y = 0.54, width = 0.90, height = 0.78))
        draw(p2a, newpage = FALSE)
        grid::popViewport()
    })
)
fig2_heatmap <- wrap_elements(
    grid::grid.grabExpr({
        grid::pushViewport(grid::viewport(x = 0.51, y = 0.48, width = 0.92, height = 0.84))
        grid::grid.draw(ggplotGrob(p2b))
        grid::popViewport()
    })
)
fig2 <- (fig2_upset | fig2_heatmap) +
    plot_layout(widths = c(1.05, 0.95))
save_plot(fig2, "figure2", fig_dims$figure2[1], fig_dims$figure2[2])




# --- figure 3 data: downsampling, concordance, credible sets ---

ds_paths   <- Sys.glob(inputs$ds_glob)
ds_sig <- rbindlist(lapply(ds_paths, function(perm_path) {
    threshold_dt <- fread(perm_path, select = c("phenotype_id", "qval", "pval_nominal_threshold"))
    threshold_dt <- threshold_dt[is.finite(qval) & qval < ds_q_threshold & is.finite(pval_nominal_threshold),
        .(phenotype_id, pval_nominal_threshold)]
    if (nrow(threshold_dt) == 0) return(data.table(phenotype_id = character(), variant_id = character(),
        slope = numeric(), pval_nominal = numeric(), seed = integer()))
    seed <- as.integer(sub(".*seed_([0-9]+).*", "\\1", perm_path))
    nominal_paths <- Sys.glob(file.path(dirname(dirname(perm_path)), "nominal", "*.parquet"))
    nominal_dt <- as.data.table(open_dataset(nominal_paths, format = "parquet") |>
        dplyr::filter(phenotype_id %in% threshold_dt$phenotype_id) |>
        dplyr::select(phenotype_id, variant_id, pval_nominal, slope) |>
        dplyr::collect())
    nominal_dt <- merge(nominal_dt, threshold_dt, by = "phenotype_id")
    nominal_dt <- nominal_dt[is.finite(pval_nominal) & pval_nominal <= pval_nominal_threshold]
    nominal_dt[, seed := seed]
    add_pair_keys(nominal_dt)
}), fill = TRUE)

full_post_keys <- unique(sig_pairs[stage == "postnatal", .(pair_key)])
ds_metrics <- ds_sig[, .(
    ds_pairs = uniqueN(pair_key),
    tp = uniqueN(pair_key[pair_key %in% full_post_keys$pair_key])
), by = seed]
ds_metrics[, `:=`(fp = ds_pairs - tp, fn = nrow(full_post_keys) - tp)]
ds_metrics[, `:=`(precision = tp / (tp + fp), recall = tp / (tp + fn))]
ds_metrics[, F1 := 2 * precision * recall / (precision + recall)]

ds_long <- melt(ds_metrics, id.vars = "seed", measure.vars = c("precision", "recall", "F1"),
    variable.name = "metric", value.name = "score")
ds_long[, metric := factor(metric, levels = c("precision", "recall", "F1"))]
ds_summary <- ds_long[, .(median = median(score, na.rm = TRUE)), by = metric]

full_post_beta  <- unique(sig_pairs[stage == "postnatal", .(pair_key, full_abs_beta = abs(slope))], by = "pair_key")
ds_beta         <- merge(ds_sig[, .(seed, pair_key, ds_abs_beta = abs(slope))], full_post_beta, by = "pair_key")
ds_beta_summary <- ds_beta[, .(ds_median_abs_beta = median(ds_abs_beta, na.rm = TRUE),
                                full_abs_beta = full_abs_beta[1L]), by = pair_key]


# --- figure 3A: downsampling precision/recall/F1 ---

f1_grid <- as.data.table(expand.grid(
  precision = seq(0.85, 1, length.out = 100),
  recall    = seq(0.10, 0.45, length.out = 100)
))
f1_grid[, F1 := 2 * precision * recall / (precision + recall)]
p3a_note <- data.table(
    precision = 0.852,
    recall = 0.435,
    label = sprintf("Median: precision %.1f%%\nrecall %.1f%%, F1 %.2f",
        ds_summary[metric == "precision", median] * 100,
        ds_summary[metric == "recall", median] * 100,
        ds_summary[metric == "F1", median])
)

p3a <- ggplot() +
    geom_contour_filled(data = f1_grid, aes(x = precision, y = recall, z = F1), alpha = 0.12) +
    geom_point(data = ds_metrics, aes(x = precision, y = recall), 
               color = "grey15", size = 2.5, alpha = 0.85, stroke = 0.4, shape = 21, fill = stage_colors[["postnatal"]]) +
    geom_vline(xintercept = ds_summary[metric == "precision", median], linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = ds_summary[metric == "recall", median], linetype = "dashed", color = "grey40", linewidth = 0.4) +
    annotate("point", x = ds_summary[metric == "precision", median], y = ds_summary[metric == "recall", median], 
             color = "#D55E00", size = 3.5, shape = 18) +
    geom_label(data = p3a_note, aes(x = precision, y = recall, label = label),
               inherit.aes = FALSE, size = 2.6, fontface = "bold",
               hjust = 0, vjust = 1, linewidth = 0, fill = alpha("white", 0.72)) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0.85, 1.00), expand = c(0.01, 0.01)) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.10, 0.45), expand = c(0.01, 0.01)) +
    scale_fill_viridis_d(option = "E", name = "F1 Background", guide = "none") +
    labs(x = "Precision", y = "Recall (Power)") +
    theme_classic(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank()
    )

# --- figure 3B: beta compare ---
p3b_summary <- merge(
    ds_beta_summary[, .(pair_key, ds_median_abs_beta)],
    shared_beta[, .(pair_key, abs_prenatal)],
    by = "pair_key"
)
p3b_summary <- p3b_summary[is.finite(abs_prenatal) & is.finite(ds_median_abs_beta)]
wilcox_p3b <- wilcox.test(p3b_summary$abs_prenatal, p3b_summary$ds_median_abs_beta, paired = TRUE)

p3b_summary[, delta_beta := ds_median_abs_beta - abs_prenatal]
p3b_summary <- p3b_summary[is.finite(delta_beta)]
median_delta <- median(p3b_summary$delta_beta, na.rm = TRUE)
p3b_xlim <- quantile(p3b_summary$delta_beta, c(0.005, 0.995), na.rm = TRUE)
p3b <- ggplot(p3b_summary, aes(x = delta_beta)) +
    geom_density(fill = "#009E73", alpha = 0.25, color = "#009E73", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = median_delta, color = "#D55E00", linewidth = 0.6) +
    coord_cartesian(xlim = p3b_xlim) +
    annotate("text", 
             x = quantile(p3b_summary$delta_beta, 0.01, na.rm = TRUE), 
             y = max(density(p3b_summary$delta_beta, na.rm = TRUE)$y) * 0.9, 
             hjust = 0,
             label = sprintf("Paired Wilcoxon n = %s\np = %s\nMedian = %.4f", 
                             comma(nrow(p3b_summary)), 
                             format_p(wilcox_p3b$p.value),
                             median_delta),
             size = 3.2, fontface = "bold") +
    labs(x = "Beta (|postnatal DS| - |prenatal|)", y = "Density") +
    theme_classic(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank()
    )

# --- figure 3C: chromosome-faceted CS interval overlap ---


pre_cs_input <- prenatal_susie[, .(pip = max(pip, na.rm = TRUE)),
    by = .(gene_key, variant_key, cs_id)]
post_cs_input <- copy(postnatal_ds_susie)
post_cs_input[, cs_id := gene_key]
pre_intervals <- cs_intervals(pre_cs_input, "Prenatal")
post_intervals <- cs_intervals(post_cs_input, "Postnatal")
pre_intervals[, chr_num := suppressWarnings(as.integer(sub("^chr", "", chr)))]
post_intervals[, chr_num := suppressWarnings(as.integer(sub("^chr", "", chr)))]
pre_intervals <- pre_intervals[chr_num %between% c(1L, 22L)]
post_intervals <- post_intervals[chr_num %between% c(1L, 22L)]

pre_intervals <- pre_intervals[, .(
    start = min(start),
    end = max(end),
    center = median(center, na.rm = TRUE),
    max_pip = max(max_pip, na.rm = TRUE),
    n_cs = uniqueN(cs_id)
), by = .(gene_key, chr, chr_num, stage)]
post_intervals <- post_intervals[, .(
    start = min(start),
    end = max(end),
    center = median(center, na.rm = TRUE),
    max_pip = max(max_pip, na.rm = TRUE),
    n_cs = uniqueN(cs_id)
), by = .(gene_key, chr, chr_num, stage)]

cs_plot_dt <- rbindlist(list(
    pre_intervals[, .(gene_key, chr, chr_num, start, end, y = 0.48, evidence = "Prenatal CS")],
    post_intervals[, .(gene_key, chr, chr_num, start, end, y = 0.48, evidence = "Postnatal DS CS")]
), fill = TRUE)

shared_dt <- merge(
    pre_intervals,
    post_intervals,
    by = c("gene_key", "chr", "chr_num"),
    suffixes = c("_pre", "_post"),
    allow.cartesian = TRUE
)
shared_dt <- shared_dt[start_pre <= end_post & start_post <= end_pre]
if (nrow(shared_dt) > 0) {
    shared_dt[, `:=`(
        start = pmax(start_pre, start_post),
        end = pmin(end_pre, end_post),
        y = 0.48,
        evidence = "Shared overlap",
        shared_score = pmin(max_pip_pre, max_pip_post)
    )]
    cs_plot_dt <- rbindlist(list(
        cs_plot_dt,
        shared_dt[, .(gene_key, chr, chr_num, start, end, y, evidence)]
    ), fill = TRUE)
}

cs_plot_dt[, evidence := factor(evidence,
    levels = c("Prenatal CS", "Postnatal DS CS", "Shared overlap"))]
cs_plot_dt[, chr := factor(chr, levels = paste0("chr", 1:22))]

chr_backbone_dt <- cs_plot_dt[, .(chr_end = max(end / 1e6, na.rm = TRUE)), by = .(chr, chr_num)]
chr_backbone_dt[, `:=`(
    chr = factor(chr, levels = paste0("chr", 1:22)),
    xstart = 0,
    xend = pmax(chr_end, 1),
    chr_y = (23 - chr_num) * figure3c_row_spacing
)]
chr_band_dt <- copy(chr_backbone_dt)
chr_band_dt[, `:=`(
    xmin = -Inf,
    xmax = Inf,
    ymin = chr_y - 0.42,
    ymax = chr_y + 0.42
)]

block_dt <- rbindlist(list(
    pre_intervals[, .(gene_key, chr, chr_num,
                      xmin = start / 1e6, xmax = end / 1e6,
                      ymin = (23 - chr_num) * figure3c_row_spacing - 0.24,
                      ymax = (23 - chr_num) * figure3c_row_spacing + 0.24,
                      stage = "Prenatal CS")],
    post_intervals[, .(gene_key, chr, chr_num,
                       xmin = start / 1e6, xmax = end / 1e6,
                       ymin = (23 - chr_num) * figure3c_row_spacing - 0.24,
                       ymax = (23 - chr_num) * figure3c_row_spacing + 0.24,
                       stage = "Postnatal DS CS")]
), fill = TRUE)
block_dt[, stage := factor(stage, levels = c("Prenatal CS", "Postnatal DS CS"))]
block_dt[, center_mb := (xmin + xmax) / 2]
block_dt[, width_mb := pmax(xmax - xmin, figure3c_min_block_width_mb)]
block_dt[, `:=`(
    xmin_plot = center_mb - width_mb / 2,
    xmax_plot = center_mb + width_mb / 2
)]

overlap_label_dt <- shared_dt[is.finite(shared_score) & !grepl("^ENSG", gene_key)][order(chr_num, -shared_score, gene_key)]
overlap_label_dt <- overlap_label_dt[, head(.SD, figure3c_max_labels_per_chr), by = chr_num]
overlap_label_dt <- overlap_label_dt[order(-shared_score, chr_num, gene_key)][seq_len(min(.N, figure3c_top_overlap_label_n))]
overlap_label_dt[, `:=`(
    x = (start + end) / 2e6,
    y_target = (23 - chr_num) * figure3c_row_spacing,
    y = (23 - chr_num) * figure3c_row_spacing,
    label_side = fifelse((start + end) / 2e6 < 125, "left", "right"),
    x_label = fifelse((start + end) / 2e6 < 125, figure3c_left_label_x, figure3c_right_label_x),
    hjust_label = fifelse((start + end) / 2e6 < 125, 0, 1),
    label = gene_key,
    label_type = "High overlap"
)]

p3c <- ggplot() +
    geom_rect(data = chr_band_dt,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#F8FAFC",
              inherit.aes = FALSE, color = NA, alpha = 0.75) +
    geom_segment(data = chr_backbone_dt,
                 aes(x = xstart, xend = xend, y = chr_y, yend = chr_y),
                 linewidth = 4.8, lineend = "round", color = "#ECEFF3") +
    geom_rect(data = block_dt[stage == "Postnatal DS CS"],
              aes(xmin = xmin_plot, xmax = xmax_plot, ymin = ymin, ymax = ymax, fill = stage),
              inherit.aes = FALSE, color = NA, alpha = 0.38) +
    geom_rect(data = block_dt[stage == "Prenatal CS"],
              aes(xmin = xmin_plot, xmax = xmax_plot, ymin = ymin, ymax = ymax, fill = stage),
              inherit.aes = FALSE, color = NA, alpha = 0.68) +
    geom_segment(data = overlap_label_dt,
        aes(x = x_label, xend = x, y = y, yend = y_target),
        inherit.aes = FALSE, color = "grey35", linewidth = 0.20,
        linetype = "dashed") +
    geom_point(data = overlap_label_dt,
        aes(x = x, y = y_target),
        inherit.aes = FALSE, color = "grey25", size = 0.55, alpha = 0.75) +
    geom_label(data = overlap_label_dt,
        aes(x = x_label, y = y, label = label, hjust = hjust_label),
        inherit.aes = FALSE, size = 2.65, color = "grey12",
        fontface = "bold", vjust = 0.5, fill = alpha("white", 0.86),
        linewidth = 0.12, label.r = unit(0.06, "lines"),
        label.padding = unit(0.06, "lines")) +
    scale_fill_manual(values = c(
        `Prenatal CS` = "#E69F00",
        `Postnatal DS CS` = "#0072B2"
    )) +
    scale_x_continuous(labels = label_number(accuracy = 1),
                       limits = c(figure3c_left_label_x - 1, figure3c_right_label_x + 1),
                       breaks = seq(0, 250, by = 50),
                       expand = c(0, 0)) +
    scale_y_continuous(
        limits = c(0.35, 22.95 * figure3c_row_spacing),
        breaks = chr_backbone_dt$chr_y,
        labels = as.character(chr_backbone_dt$chr),
        expand = c(0, 0)
    ) +
    labs(x = "Genomic position (Mb)", y = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(size = 7, color = "grey35"),
        axis.text.y = element_text(size = 7.2, face = "bold", color = "grey35"),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(color = "grey70", linewidth = 0.25),
        legend.position = "top",
        legend.title = element_blank(),
        legend.key.width = unit(1.10, "cm"),
        legend.text = element_text(size = 8),
        plot.margin = margin(0, 8, 3, 8)
    )

fig3_top <- (p3a | p3b) + plot_layout(widths = c(0.82, 1.72))
fig3 <- fig3_top / p3c + plot_layout(heights = c(1.42, 1.24))
save_plot(fig3, "figure3", fig_dims$figure3[1], fig_dims$figure3[2])


# --- figure 4: stage-specific example loci ---
# Stage-specific = same allele-phenotype tested in both stages, significant in exactly one stage.
sig_pairs_phenotype <- unique(rbindlist(list(prenatal_sig, postnatal_sig), fill = TRUE),
                              by = c("stage", "allele_phenotype_key"))

pre_sig  <- unique(sig_pairs_phenotype[stage == "prenatal",  .(allele_phenotype_key)])
post_sig <- unique(sig_pairs_phenotype[stage == "postnatal", .(allele_phenotype_key)])

pre_in_post_nominal  <- stream_nominal_allele_phenotype_matches(inputs$postnatal_nominal, pre_sig)
post_in_pre_nominal  <- stream_nominal_allele_phenotype_matches(inputs$prenatal_nominal, post_sig)

pre_only <- setdiff(pre_in_post_nominal$allele_phenotype_key, post_sig$allele_phenotype_key)
post_only <- setdiff(post_in_pre_nominal$allele_phenotype_key, pre_sig$allele_phenotype_key)

pre_dt  <- unique(sig_pairs_phenotype[stage == "prenatal"  & allele_phenotype_key %in% pre_only,
                  .(phenotype_id, phenotype_key, variant_id, variant_key, pair_key, allele_phenotype_key,
                    significant_stage_pval_nominal = pval_nominal,
                    effect = slope, abs_effect = abs(slope))], by = "allele_phenotype_key")
post_dt <- unique(sig_pairs_phenotype[stage == "postnatal" & allele_phenotype_key %in% post_only,
                  .(phenotype_id, phenotype_key, variant_id, variant_key, pair_key, allele_phenotype_key,
                    significant_stage_pval_nominal = pval_nominal,
                    effect = slope, abs_effect = abs(slope))], by = "allele_phenotype_key")

pre_dt <- merge(pre_dt, pre_in_post_nominal, by = "allele_phenotype_key", all.x = TRUE)
post_dt <- merge(post_dt, post_in_pre_nominal, by = "allele_phenotype_key", all.x = TRUE)
pre_dt[, `:=`(significant_stage = "prenatal", other_stage = "postnatal", tested_in_other_stage = TRUE)]
post_dt[, `:=`(significant_stage = "postnatal", other_stage = "prenatal", tested_in_other_stage = TRUE)]


loci_specific <- list(
    prenatal_specific  = enrich_locus(pre_dt[order(-abs_effect)]),
    postnatal_specific = enrich_locus(post_dt[order(-abs_effect)])
)
if (nrow(loci_specific$prenatal_specific) == 0 || nrow(loci_specific$postnatal_specific) == 0)
    stop("No stage-specific loci found — check sig_pairs data")

bam_files <- unlist(lapply(Sys.glob(file.path(inputs$postnatal_bam_root, "chb*_bam")),
    function(d) list.files(d, pattern = "Aligned\\.sortedByCoord\\.out\\.bam$", full.names = TRUE)))
bam_files <- bam_files[file.exists(paste0(bam_files, ".bai"))]
bam_id <- sub("[A-D]$", "", sub("Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(bam_files)))
bam_files_by_sample <- setNames(bam_files, bam_id)
bw_files <- list.files(inputs$prenatal_bw_dir, pattern = "\\.bw$", full.names = TRUE)
bw_files_by_sample <- setNames(bw_files, sub("\\.bw$", "", basename(bw_files)))

prenatal_sample_list  <- fread(inputs$prenatal_samples,  header = FALSE)$V1
postnatal_sample_list <- fread(inputs$postnatal_samples, header = FALSE)$V1

prenatal_bed_full  <- fread(inputs$prenatal_bed,
    colClasses = c(`#Chr` = "integer", start = "integer", end = "integer", pid = "character"))
setnames(prenatal_bed_full,  c("#Chr", "start", "end", "pid"), c("chr", "start", "end", "pid"))
postnatal_bed_full <- fread(inputs$postnatal_bed,
    colClasses = c(`#Chr` = "integer", start = "integer", end = "integer", pid = "character"))
setnames(postnatal_bed_full, c("#Chr", "start", "end", "pid"), c("chr", "start", "end", "pid"))
prenatal_bed_full[, phenotype_key := make_phenotype_key(pid)]
postnatal_bed_full[, phenotype_key := make_phenotype_key(pid)]

star_pre <- pick_star_example(loci_specific$prenatal_specific, "prenatal")
star_post <- pick_star_example(loci_specific$postnatal_specific, "postnatal")
stage_examples <- list(
    prenatal_specific  = make_ex(star_pre, stage_specific = "prenatal"),
    postnatal_specific = make_ex(star_post, stage_specific = "postnatal")
)

fig4a <- build_example_figure(stage_examples$prenatal_specific,
    bw_files_by_sample, bam_files_by_sample,
    inputs$prenatal_pfile, inputs$postnatal_pfile,
    prenatal_sample_list, postnatal_sample_list,
    prenatal_bed_full, postnatal_bed_full, pas_info)
fig4b <- build_example_figure(stage_examples$postnatal_specific,
    bw_files_by_sample, bam_files_by_sample,
    inputs$prenatal_pfile, inputs$postnatal_pfile,
    prenatal_sample_list, postnatal_sample_list,
    prenatal_bed_full, postnatal_bed_full, pas_info)

fig4 <- (fig4a | fig4b) +
    plot_layout(widths = c(1, 1), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.margin = margin(6, 12, 6, 6))
save_plot(fig4, "figure4", fig_dims$figure4[1], fig_dims$figure4[2])
