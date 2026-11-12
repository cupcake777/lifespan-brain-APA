#!/usr/bin/env Rscript
# APA-QTL x GWAS SMR locus figures from SMR --plot raw output.
# Reads data/07_coloc/smr_gene_significant.tsv; writes figures/07_coloc and data/07_coloc tables.
if (!all(nzchar(Sys.getenv(c("APA_ROOT", "APA_HPC_HOME"))))) stop("Source config/paths.sh first")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- helpers ---

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") return(default)
  as.integer(value)
}

env_logical <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (value == "") return(default)
  value %in% c("1", "true", "t", "yes", "y")
}

read_fixed_section <- function(lines, start, n, col_names) {
  if (n <= 0) {
    out <- data.table()
    for (nm in col_names) out[, (nm) := vector()]
    return(out)
  }
  dt <- fread(text = paste(lines[start:(start + n - 1)], collapse = "\n"), sep = " ", header = FALSE)
  setnames(dt, col_names)
  dt
}

theme_apa_smr <- function(base_size = 8.8) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9.5, face = "bold", hjust = 0, margin = ggplot2::margin(0, 0, 2, 0)),
      plot.tag = ggplot2::element_text(size = 10, face = "bold"),
      axis.title = ggplot2::element_text(size = 8.4),
      axis.text = ggplot2::element_text(size = 8.0, color = "black"),
      legend.title = ggplot2::element_text(size = 8.0),
      legend.text = ggplot2::element_text(size = 7.6),
      legend.position = "right",
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.2),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 5, 4, 4)
    )
}

add_ld_bins <- function(dt) {
  dt[, ld_r2 := ld_to_top^2]
  dt[, ld_bin := cut(
    ld_r2,
    breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8"),
    right = FALSE
  )]
  dt[!is.finite(ld_r2), ld_bin := NA_character_]
  dt[, ld_bin := factor(ld_bin, levels = c(">=0.8", "0.6-0.8", "0.4-0.6", "0.2-0.4", "<0.2"))]
  dt
}

# --- config ---

MAX_LOCI <- NA_integer_
PLOT_WINDOW_KB <- 500
PEQTL_SMR <- 1e-3
PEQTL_HEIDI <- 1.5e-3
HEIDI_MIN_M <- 3
HEIDI_MAX_M <- 20
CIS_WIND_KB <- 2000
MAKE_PNG <- TRUE
OVERWRITE_RAW <- FALSE
MAX_LOCI <- env_int("MAX_LOCI", MAX_LOCI)
PLOT_WINDOW_KB <- env_int("PLOT_WINDOW_KB", PLOT_WINDOW_KB)
MAKE_PNG <- env_logical("MAKE_PNG", MAKE_PNG)
OVERWRITE_RAW <- env_logical("OVERWRITE_RAW", OVERWRITE_RAW)

SMR_EXE <- Sys.getenv("SMR", unset = paste0(Sys.getenv("APA_HPC_HOME"), "/bin/tools/smr-1.4.0-linux-x86_64/smr"))
PLINK2_EXE <- Sys.getenv("PLINK2", unset = "")
if (!nzchar(PLINK2_EXE)) PLINK2_EXE <- unname(Sys.which("plink2"))
if (!nzchar(PLINK2_EXE)) PLINK2_EXE <- paste0(Sys.getenv("APA_HPC_HOME"), "/.local/share/mamba/bin/plink2")
LD_KEEP_DIR <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/04_coloc/refs/plink_keep_runtime")
GENCODE_GTF <- paste0(Sys.getenv("APA_ROOT"), "/ref/gencode.v40.annotation.gtf.gz")

DATA_DIR <- paste0(Sys.getenv("APA_ROOT"), "/data/07_coloc")
FIG_DIR <- paste0(Sys.getenv("APA_ROOT"), "/figures/07_coloc")
LOG_DIR <- paste0(Sys.getenv("APA_ROOT"), "/logs")
RAW_DIR <- file.path(DATA_DIR, "raw")
LD_REF <- file.path(DATA_DIR, "refs", "plink_hg38_uniqueid_runtime", "1kgp_hg38_chrpos_ref_alt")
LD_REF_BY_CHR_DIR <- file.path(DATA_DIR, "refs", "plink_hg38_uniqueid_by_chr_runtime")
for (d in c(DATA_DIR, FIG_DIR, LOG_DIR, RAW_DIR, LD_REF_BY_CHR_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# --- significant hits ---

hits <- fread(file.path(DATA_DIR, "smr_gene_significant.tsv"))
if (!nrow(hits)) {
  message("No gene-level significant hits found; no SMR-style figures generated.")
  quit(save = "no")
}
stage_levels <- c("prenatal", "postnatal", "interaction")
hits[, p_SMR_BH_num := as.numeric(p_SMR_BH)]
hits[, gene_min_p_SMR_num := as.numeric(gene_min_p_SMR)]
hits[, stage_order := match(qtl_type, stage_levels)]
setorder(hits, stage_order, p_SMR_BH_num, gene_min_p_SMR_num, Gene, probeID)
if (!is.na(MAX_LOCI)) {
  hits <- hits[seq_len(min(MAX_LOCI, nrow(hits)))]
}
hits[, c("p_SMR_BH_num", "gene_min_p_SMR_num", "stage_order") := NULL]
if (!nrow(hits)) stop("No gene-level significant hits found.")

# --- SMR gene list ---

gene_list <- file.path(DATA_DIR, "gencode.v40.smr_gene_list.txt")
if (!(file.exists(gene_list) && file.info(gene_list)$size > 0)) {
  con <- gzfile(GENCODE_GTF, open = "rt")
  if (file.exists(gene_list)) file.remove(gene_list)
  header_written <- FALSE
  repeat {
    lines <- readLines(con, n = 100000, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[!grepl("^#", lines)]
    if (!length(lines)) next
    dt <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = FALSE)
    dt <- dt[V3 == "gene"]
    if (!nrow(dt)) next
    dt[, chrom := sub("^chr", "", V1)]
    dt[, gene := sub('.*gene_name "([^"]+)".*', "\\1", V9)]
    dt <- dt[gene != V9 & gene != ""]
    fwrite(
      dt[, .(chrom, V4, V5, gene, V7)],
      gene_list,
      sep = "\t",
      col.names = FALSE,
      append = header_written
    )
    header_written <- TRUE
  }
  close(con)
  stopifnot(file.info(gene_list)$size > 0)
}
gene_list <- normalizePath(gene_list, mustWork = TRUE)

# --- per-hit SMR plot data and figures ---

manifest <- vector("list", nrow(hits))
for (i in seq_len(nrow(hits))) {
  hit <- hits[i]
  qtl_type <- hit[["qtl_type"]]
  gwas_id <- hit[["gwas_id"]]
  ancestry <- sub("_.*$", "", gwas_id)
  probe_id <- hit[["probeID"]]
  keep <- file.path(LD_KEEP_DIR, paste0(ancestry, ".keep"))
  gwas_ma <- file.path(DATA_DIR, paste0(gwas_id, ".ma"))
  besd <- file.path(DATA_DIR, "besd", qtl_type)
  base <- sprintf("%s_%s_%s", qtl_type, gwas_id, hit[["Gene"]])
  base <- gsub("[^A-Za-z0-9._-]+", "_", base)
  base <- gsub("_+", "_", base)
  base <- gsub("^_|_$", "", base)
  expected <- file.path(RAW_DIR, "plot", paste0(base, ".", probe_id, ".txt"))

  # chromosome-specific LD reference
  ld_chr <- sub("^chr", "", as.character(hit[["ProbeChr"]]))
  ld_chr <- sub("\\.0$", "", ld_chr)
  if (!nzchar(ld_chr) || is.na(ld_chr)) stop("Missing chromosome for LD reference selection.")
  ld_prefix <- file.path(LD_REF_BY_CHR_DIR, paste0("1kgp_hg38_chrpos_ref_alt_chr", ld_chr))
  ld_files <- paste0(ld_prefix, c(".bed", ".bim", ".fam"))
  if (!all(file.exists(ld_files) & file.info(ld_files)$size > 0)) {
    log_file <- file.path(LOG_DIR, paste0(basename(ld_prefix), ".plink2.log"))
    status <- system2(
      PLINK2_EXE,
      args = shQuote(c(
        "--bfile", LD_REF,
        "--chr", ld_chr,
        "--make-bed",
        "--out", ld_prefix
      )),
      stdout = log_file,
      stderr = log_file
    )
    if (!identical(status, 0L)) {
      stop(sprintf("plink2 chromosome LD extraction failed for chr%s; see %s", ld_chr, log_file))
    }
  }
  ld_ref <- normalizePath(ld_prefix, mustWork = FALSE)

  # SMR --plot raw data, reused if present
  if (!(file.exists(expected) && file.info(expected)$size > 0 && !OVERWRITE_RAW)) {
    args <- c(
      "--bfile", ld_ref,
      "--keep", keep,
      "--gwas-summary", gwas_ma,
      "--beqtl-summary", besd,
      "--out", file.path(RAW_DIR, base),
      "--plot",
      "--probe", probe_id,
      "--probe-wind", as.character(PLOT_WINDOW_KB),
      "--heidi-mtd", "1",
      "--heidi-min-m", as.character(HEIDI_MIN_M),
      "--heidi-max-m", as.character(HEIDI_MAX_M),
      "--cis-wind", as.character(CIS_WIND_KB),
      "--peqtl-smr", as.character(PEQTL_SMR),
      "--peqtl-heidi", as.character(PEQTL_HEIDI),
      "--gene-list", gene_list,
      "--thread-num", "2"
    )
    if (ancestry == "EUR") {
      args <- c(args, "--diff-freq-prop", Sys.getenv("EUR_DIFF_FREQ_PROP", unset = "1.0"))
    }
    log_file <- file.path(LOG_DIR, paste0(base, ".plot.log"))
    status <- system2(SMR_EXE, args = shQuote(args), stdout = log_file, stderr = log_file)
    if (!identical(status, 0L)) stop(sprintf("smr --plot failed for %s; see %s", base, log_file))
  }
  raw_file <- normalizePath(expected, mustWork = TRUE)

  # parse SMR plot file
  lines <- readLines(raw_file, warn = FALSE)
  pos <- 1
  probe_header <- strsplit(lines[pos], "\\s+")[[1]]
  if (probe_header[1] != "$probe") stop("Unexpected plot file: missing $probe")
  n_probe <- as.integer(probe_header[2])
  target_probe <- probe_header[3]
  pos <- pos + 1
  raw_probe <- read_fixed_section(lines, pos, n_probe, c("probeID", "chr", "bp", "gene", "gene_start", "gene_end", "strand", "p_smr", "p_heidi"))
  pos <- pos + n_probe

  snp_header <- strsplit(lines[pos], "\\s+")[[1]]
  if (snp_header[1] != "$SNP") stop("Unexpected plot file: missing $SNP")
  n_snp <- as.integer(snp_header[2])
  pos <- pos + 1
  raw_snp <- read_fixed_section(lines, pos, n_snp, c("SNP", "chr", "bp", "A1", "A2"))
  pos <- pos + n_snp

  gwas_header <- strsplit(lines[pos], "\\s+")[[1]]
  if (gwas_header[1] != "$GWAS") stop("Unexpected plot file: missing $GWAS")
  n_gwas <- as.integer(gwas_header[2])
  pos <- pos + 1
  raw_gwas <- read_fixed_section(lines, pos, n_gwas, c("SNP", "b_gwas", "se_gwas"))
  pos <- pos + n_gwas

  apa_header <- strsplit(lines[pos], "\\s+")[[1]]
  # SMR's raw plot format calls the molecular-QTL block "$eQTL" regardless of
  # whether the QTL trait is expression, APA, methylation, or another feature.
  if (apa_header[1] != "$eQTL") stop("Unexpected plot file: missing SMR molecular-QTL section")
  n_apa_probe <- as.integer(apa_header[2])
  pos <- pos + 1
  apa_list <- vector("list", n_apa_probe)
  for (k in seq_len(n_apa_probe)) {
    sub_header <- strsplit(lines[pos], "\\s+")[[1]]
    apa_probe <- sub_header[1]
    n_apa <- as.integer(sub_header[2])
    pos <- pos + 1
    dt <- read_fixed_section(lines, pos, n_apa, c("SNP", "b_apa", "se_apa", "ld_to_top"))
    dt[, probeID := apa_probe]
    apa_list[[k]] <- dt
    pos <- pos + n_apa
  }
  raw_apa <- rbindlist(apa_list, fill = TRUE)

  raw_gene <- data.table()
  if (pos <= length(lines)) {
    gene_header <- strsplit(lines[pos], "\\s+")[[1]]
    if (length(gene_header) >= 2 && gene_header[1] == "$Gene") {
      n_gene <- as.integer(gene_header[2])
      pos <- pos + 1
      raw_gene <- read_fixed_section(lines, pos, n_gene, c("chr", "start", "end", "gene", "strand"))
    }
  }
  for (cc in c("chr", "bp", "gene_start", "gene_end", "p_smr", "p_heidi")) raw_probe[, (cc) := as.numeric(get(cc))]
  for (cc in c("chr", "bp")) raw_snp[, (cc) := as.numeric(get(cc))]
  for (cc in c("b_gwas", "se_gwas")) raw_gwas[, (cc) := as.numeric(get(cc))]
  for (cc in c("b_apa", "se_apa", "ld_to_top")) raw_apa[, (cc) := as.numeric(get(cc))]
  if (nrow(raw_gene)) for (cc in c("chr", "start", "end")) raw_gene[, (cc) := as.numeric(get(cc))]

  # plot tables
  snp <- copy(raw_snp)
  gwas <- merge(copy(raw_gwas), snp, by = "SNP", all.x = TRUE)
  gwas <- gwas[is.finite(b_gwas) & is.finite(se_gwas) & se_gwas > 0 & is.finite(bp)]
  gwas[, p_gwas := 2 * stats::pnorm(-abs(b_gwas / se_gwas))]
  gwas[, logp := pmin(-log10(pmax(p_gwas, .Machine$double.xmin)), 50)]

  apa <- merge(copy(raw_apa[probeID == target_probe]), snp, by = "SNP", all.x = TRUE)
  apa <- apa[is.finite(b_apa) & is.finite(se_apa) & se_apa > 0 & is.finite(bp)]
  apa[, p_apa := 2 * stats::pnorm(-abs(b_apa / se_apa))]
  apa[, logp := pmin(-log10(pmax(p_apa, .Machine$double.xmin)), 50)]

  top_snp <- hit[["topSNP"]] %||% raw_probe[probeID == target_probe][["probeID"]][1]
  if (!top_snp %in% gwas$SNP || !top_snp %in% apa$SNP) {
    common <- merge(apa[, .(SNP, z_apa = abs(b_apa / se_apa))], gwas[, .(SNP)], by = "SNP")
    if (nrow(common)) top_snp <- common[which.max(z_apa), SNP]
  }
  gwas[, point_role := fifelse(SNP == top_snp, "Top SNP", "Other SNPs")]
  apa[, point_role := fifelse(SNP == top_snp, "Top SNP", "Other SNPs")]

  gwas <- merge(gwas, apa[, .(SNP, ld_to_top)], by = "SNP", all.x = TRUE)
  gwas <- add_ld_bins(gwas)
  apa <- add_ld_bins(apa)

  effect <- merge(
    apa[, .(SNP, b_apa, se_apa)],
    gwas[, .(SNP, b_gwas, se_gwas)],
    by = "SNP"
  )
  effect <- effect[is.finite(b_apa) & is.finite(se_apa) & is.finite(b_gwas) & is.finite(se_gwas)]
  effect[, point_role := fifelse(SNP == top_snp, "Top SNP", "Shared SNPs")]
  effect[, z_apa := abs(b_apa / se_apa)]

  # gene lanes
  genes <- copy(raw_gene)
  if (!nrow(genes)) {
    genes[, lane := integer()]
  } else {
    genes <- copy(genes[order(start, end)])
    lane_end <- numeric()
    genes[, lane := 1L]
    for (g in seq_len(nrow(genes))) {
      placed <- FALSE
      for (l in seq_along(lane_end)) {
        if (genes$start[g] > lane_end[l] + 20000) {
          genes$lane[g] <- l
          lane_end[l] <- genes$end[g]
          placed <- TRUE
          break
        }
      }
      if (!placed) {
        lane_end <- c(lane_end, genes$end[g])
        genes$lane[g] <- length(lane_end)
      }
    }
  }
  if (nrow(genes)) genes[, label_gene := gene == hit[["Gene"]]] else genes[, label_gene := logical()]
  genes[, x_mid := (start + end) / 2]

  # figure
  trait <- sub("^[A-Z]+_", "", gwas_id)
  chr <- unique(na.omit(gwas$chr))[1]
  x_limits <- range(c(gwas$bp, apa$bp, genes$start, genes$end), na.rm = TRUE)
  x_breaks <- pretty(x_limits / 1e6, n = 4)
  top_color <- "#D55E00"
  other_color <- "#7F8D99"
  probe_color <- "#0072B2"
  ld_colors <- c(
    ">=0.8" = "#D55E00",
    "0.6-0.8" = "#E69F00",
    "0.4-0.6" = "#009E73",
    "0.2-0.4" = "#56B4E9",
    "<0.2" = "#7F8D99"
  )
  top_bp <- gwas[point_role == "Top SNP", bp][1] %||% apa[point_role == "Top SNP", bp][1]
  probe_bp <- raw_probe[probeID == target_probe, bp][1]
  ld_legend_seed <- data.table(
    bp = rep(x_limits[1], length(ld_colors)),
    logp = 0,
    ld_bin = factor(names(ld_colors), levels = names(ld_colors))
  )
  hit_label <- sprintf(
    "%s | %s %s | SMR q=%s, HEIDI P=%s",
    hit[["Gene"]],
    ancestry,
    trait,
    formatC(as.numeric(hit[["p_SMR_BH"]]), format = "e", digits = 2),
    formatC(as.numeric(hit[["p_HEIDI_num"]]), format = "e", digits = 2)
  )

  p_gwas <- ggplot2::ggplot(gwas, ggplot2::aes(bp / 1e6, logp)) +
    ggplot2::geom_point(data = ld_legend_seed, ggplot2::aes(color = ld_bin), size = 1.75, alpha = 0, show.legend = TRUE) +
    ggplot2::geom_point(data = gwas[point_role != "Top SNP" & is.na(ld_bin)], color = "#C7CDD2", size = 1.75, alpha = 0.70) +
    ggplot2::geom_point(data = gwas[point_role != "Top SNP" & !is.na(ld_bin)], ggplot2::aes(color = ld_bin), size = 1.75, alpha = 0.82) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey55", linewidth = 0.28) +
    ggplot2::geom_vline(xintercept = probe_bp / 1e6, color = probe_color, linewidth = 0.32) +
    ggplot2::geom_vline(xintercept = top_bp / 1e6, color = top_color, linewidth = 0.32, linetype = "22") +
    ggplot2::geom_point(data = gwas[point_role == "Top SNP"], color = top_color, fill = "white", shape = 21, size = 2.6, stroke = 0.65) +
    ggplot2::scale_color_manual(values = ld_colors, drop = FALSE, name = expression(LD~r^2), guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.4, alpha = 1))) +
    ggplot2::scale_x_continuous(limits = x_limits / 1e6, breaks = x_breaks) +
    ggplot2::labs(title = "GWAS regional association", x = NULL, y = "GWAS -log10(P)") +
    theme_apa_smr() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none")

  p_apa <- ggplot2::ggplot(apa, ggplot2::aes(bp / 1e6, logp)) +
    ggplot2::geom_point(data = ld_legend_seed, ggplot2::aes(color = ld_bin), size = 1.75, alpha = 0, show.legend = TRUE) +
    ggplot2::geom_point(data = apa[point_role != "Top SNP" & is.na(ld_bin)], color = "#C7CDD2", size = 1.75, alpha = 0.70) +
    ggplot2::geom_point(data = apa[point_role != "Top SNP" & !is.na(ld_bin)], ggplot2::aes(color = ld_bin), size = 1.75, alpha = 0.82) +
    ggplot2::geom_vline(xintercept = probe_bp / 1e6, color = probe_color, linewidth = 0.32) +
    ggplot2::geom_vline(xintercept = top_bp / 1e6, color = top_color, linewidth = 0.32, linetype = "22") +
    ggplot2::geom_point(data = apa[point_role == "Top SNP"], color = top_color, fill = "white", shape = 21, size = 2.6, stroke = 0.65) +
    ggplot2::scale_color_manual(values = ld_colors, drop = FALSE, name = expression(LD~r^2), guide = ggplot2::guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 2.4, alpha = 1))) +
    ggplot2::scale_x_continuous(limits = x_limits / 1e6, breaks = x_breaks) +
    ggplot2::labs(title = "APA-QTL regional association", x = NULL, y = "APA-QTL -log10(P)") +
    theme_apa_smr() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

  if (nrow(genes)) {
    genes[, x0 := fifelse(strand == "-", end, start)]
    genes[, x1 := fifelse(strand == "-", start, end)]
    p_gene <- ggplot2::ggplot(genes, ggplot2::aes(y = lane)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = x0 / 1e6, xend = x1 / 1e6, yend = lane, alpha = label_gene),
        linewidth = 0.32,
        arrow = ggplot2::arrow(length = grid::unit(0.045, "inches"), ends = "last", type = "closed")
      ) +
      ggplot2::geom_text(
        data = genes[label_gene == TRUE],
        ggplot2::aes(x = x_mid / 1e6, label = gene),
        size = 2.15,
        vjust = -0.55,
        check_overlap = TRUE
      ) +
      ggplot2::scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.32), guide = "none") +
      ggplot2::scale_x_continuous(limits = x_limits / 1e6, breaks = x_breaks) +
      ggplot2::scale_y_reverse() +
      ggplot2::labs(title = "Gene context", x = sprintf("Chromosome %s position (Mb)", chr), y = "Genes") +
      theme_apa_smr() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), legend.position = "none")
  } else {
    p_gene <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(x = sprintf("Chromosome %s position (Mb)", chr))
  }

  top <- effect[point_role == "Top SNP"]
  slope <- if (nrow(top) && is.finite(top$b_apa[1]) && top$b_apa[1] != 0) top$b_gwas[1] / top$b_apa[1] else NA_real_
  effect_bg <- effect[point_role != "Top SNP"]
  effect_top <- effect[point_role == "Top SNP"]
  p_effect <- ggplot2::ggplot(effect, ggplot2::aes(b_apa, b_gwas)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey86", linewidth = 0.28) +
    ggplot2::geom_vline(xintercept = 0, color = "grey86", linewidth = 0.28) +
    ggplot2::geom_point(data = effect_bg, color = other_color, size = 1.1, alpha = 0.42) +
    ggplot2::geom_errorbar(data = effect_top, ggplot2::aes(ymin = b_gwas - se_gwas, ymax = b_gwas + se_gwas), width = 0, color = top_color, linewidth = 0.55) +
    ggplot2::geom_errorbar(data = effect_top, ggplot2::aes(xmin = b_apa - se_apa, xmax = b_apa + se_apa), orientation = "y", width = 0, color = top_color, linewidth = 0.55) +
    ggplot2::geom_point(data = effect_top, color = top_color, fill = "white", shape = 21, size = 3.3, stroke = 0.8) +
    ggplot2::labs(title = "EFFECT", x = "APA-QTL effect", y = "GWAS effect") +
    theme_apa_smr() +
    ggplot2::theme(legend.position = "none")
  if (is.finite(slope)) {
    p_effect <- p_effect + ggplot2::geom_abline(slope = slope, intercept = 0, color = top_color, linewidth = 0.55, linetype = "22")
  }

  max_gene_lanes <- if (nrow(genes)) max(genes$lane, na.rm = TRUE) else 0
  gene_track_weight <- min(0.95, max(0.45, 0.28 + 0.045 * max_gene_lanes))
  figure_height <- max(4.6, 4.25 + 0.14 * max(max_gene_lanes - 6, 0))
  right <- p_gwas / p_apa / p_gene + patchwork::plot_layout(heights = c(1, 1, gene_track_weight), guides = "collect")
  fig <- p_effect | right
  fig <- fig + patchwork::plot_layout(widths = c(1, 1.8)) &
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.spacing.x = grid::unit(0.16, "inches")
    )
  fig <- fig + patchwork::plot_annotation(title = hit_label, theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 11, hjust = 0)))
  combined_pdf <- file.path(FIG_DIR, paste0(base, "_apa_smr_main.pdf"))
  combined_png <- file.path(FIG_DIR, paste0(base, "_apa_smr_main.png"))
  ggplot2::ggsave(combined_pdf, fig, width = 7.4, height = figure_height, device = grDevices::cairo_pdf)
  if (isTRUE(MAKE_PNG)) {
    ggplot2::ggsave(combined_png, fig, width = 7.4, height = figure_height, dpi = 320, bg = "white")
  }

  manifest[[i]] <- data.table(
    qtl_type = hit[["qtl_type"]],
    gwas_id = hit[["gwas_id"]],
    gene = hit[["Gene"]],
    probeID = hit[["probeID"]],
    topSNP = top_snp,
    raw_plot_file = raw_file,
    custom_pdf = combined_pdf,
    custom_png = if (isTRUE(MAKE_PNG)) combined_png else NA_character_,
    panel_title = hit_label,
    n_gwas_snps = nrow(gwas),
    n_apa_snps = nrow(apa),
    n_shared_snps = nrow(effect),
    max_gene_lanes = if (nrow(genes)) max(genes$lane, na.rm = TRUE) else 0,
    n_labeled_genes = if (nrow(genes)) sum(genes$label_gene, na.rm = TRUE) else 0,
    layout_note = "effect panel left; regional panels right; LD bins use squared LD correlation"
  )
}

manifest <- rbindlist(manifest, fill = TRUE)
manifest[, rank := seq_len(.N)]
fig_tbl <- manifest[, .(
  rank, qtl_type, gwas_id, gene, probeID, topSNP, custom_pdf, custom_png,
  panel_title, n_gwas_snps, n_apa_snps, n_shared_snps, max_gene_lanes,
  n_labeled_genes, layout_note
)]
fwrite(fig_tbl, file.path(DATA_DIR, "smr_locus_plot_manifest.tsv"), sep = "\t")
layout_validation <- manifest[, .(
  figure = basename(custom_png),
  visible_title_characters = nchar(panel_title),
  labeled_genes = n_labeled_genes,
  gene_lanes = max_gene_lanes,
  shared_snps = n_shared_snps,
  status = fifelse(
    nchar(panel_title) <= 90 & n_labeled_genes <= 3 & max_gene_lanes <= 14 & n_shared_snps > 0,
    "PASS",
    "REVIEW"
  )
)]
fwrite(layout_validation, file.path(DATA_DIR, "smr_locus_plot_layout_validation.tsv"), sep = "\t")
fwrite(layout_validation, file.path(DATA_DIR, "smr_style_layout_validation.tsv"), sep = "\t")
