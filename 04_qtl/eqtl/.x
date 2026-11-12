#!/usr/bin/env Rscript
for (.apa_v in c("APA_HPC_HOME")) if (!nzchar(Sys.getenv(.apa_v))) stop("Set ", .apa_v, " (source config/paths.sh)")

# 06-5: 批量导出所有 disease-colocalized 候选的 locuscompare 风格 locus 图
# 池 = 06_candidate.tsv 中 disease_coloc == TRUE 且在 04_coloc
# significant_colocalizations_PP4_0.7.csv 里有 locuscompare JSON 的候选
# 输出: fig/06_locus/<rank>_<gene>_<stage>_<trait>.{png,pdf}

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

# --- config ---

project_dir <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/11_eqtl")
res_dir <- file.path(project_dir, "results")
run_mode <- Sys.getenv("LOCUS_MODE", "followup")
out_dir <- Sys.getenv("LOCUS_OUT_DIR", file.path(project_dir, "fig", "06_locus"))
gwas_coloc_dir <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/04_coloc/results")
gtf_path <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/gene/gencode.v40.annotation.gtf.gz")
setwd(project_dir)

fig_w2 <- 7.09
fig_h <- 5.2
png_dpi <- 600
base_font <- 10
font_sm <- 9

col_lead <- "#A020F0"
col_gene <- "#4A7BA6"

ld_cols <- c(
  "<0.2" = "#1F2E8C",
  "0.2-0.4" = "#A6CEE3",
  "0.4-0.6" = "#33A02C",
  "0.6-0.8" = "#FF9F1C",
  "0.8-1.0" = "#E31A1C",
  "missing LD" = "#D3D3D3"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- helper functions ---

as_flag <- function(x) {
  y <- as.logical(x)
  fifelse(is.na(y), FALSE, y)
}

theme_06 <- function() {
  theme_classic(base_size = base_font) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = font_sm),
      axis.title = element_text(colour = "black", size = base_font),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      legend.title = element_text(size = font_sm),
      legend.text = element_text(size = font_sm),
      legend.key.size = grid::unit(0.32, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.title = element_text(size = base_font + 1, face = "bold"),
      plot.margin = margin(4, 5, 4, 4)
    )
}

# --- loading data: 候选 + locuscompare 匹配 ---

gwas_hits <- fread(file.path(gwas_coloc_dir, "significant_colocalizations_PP4_0.7.csv"))

hits <- gwas_hits[locuscompare_json_file != ""]
hits[, json_abs := file.path(gwas_coloc_dir, locuscompare_json_file)]
hits <- hits[file.exists(json_abs)]
if (run_mode == "interaction") {
  case_dt <- hits[QTL_ID == "interaction"]
  case_dt[, gene_symbol := vapply(strsplit(Phenotype, "|", fixed = TRUE), `[`, character(1), 2L)]
  case_dt[, `:=`(stage = QTL_ID, hit_trait = GWAS_ID, hit_pp4 = as.numeric(PP4),
                 max_PP_H4 = NA_real_, candidate_rank = .I)]
  cat("interaction signals with locuscompare JSON:", nrow(case_dt), "\n")
} else {
  priority <- fread(file.path(res_dir, "06_candidate.tsv"))
  cand <- priority[as_flag(disease_coloc) & gwas_best_trait != ""]
  cand[, max_PP_H4 := as.numeric(max_PP_H4)]
  setorder(hits, -PP4)
  hits_best <- hits[, .SD[1], by = .(QTL_ID, Phenotype)]
  case_dt <- merge(
    cand,
    hits_best[, .(stage = QTL_ID, phenotype_id = Phenotype,
                  hit_trait = GWAS_ID, hit_pp4 = PP4, json_abs)],
    by = c("stage", "phenotype_id")
  )
  setorder(case_dt, candidate_rank)
  keep <- data.table(
    gene_symbol = c("ARL3", "CKMT1A", "GAS5", "RAB29"),
    stage = rep("postnatal", 4L),
    hit_trait = c("EAS_SCZ", "EUR_SCZ", "EUR_MDD", "EAS_PD")
  )
  case_dt <- merge(case_dt, keep, by = c("gene_symbol", "stage", "hit_trait"))
  cat("disease-colocalized candidates:", nrow(cand),
      "; with locuscompare JSON:", nrow(hits_best),
      "; exported after keep-filter:", nrow(case_dt), "\n")
}

# --- 读全部 JSON，收集窗口 ---

locus_list <- vector("list", nrow(case_dt))
for (i in seq_len(nrow(case_dt))) {
  pts <- as.data.table(fromJSON(case_dt$json_abs[i])$points)
  pts[, ld_bin := factor(ld_bin, levels = names(ld_cols))]
  pts[, pos_mb := pos / 1e6]
  locus_list[[i]] <- pts
}
win_dt <- rbindlist(lapply(seq_along(locus_list), function(i) {
  pts <- locus_list[[i]]
  data.table(win_id = i,
             chrom = pts[!is.na(chrom), chrom[1]],
             start = min(pts$pos, na.rm = TRUE),
             end = max(pts$pos, na.rm = TRUE))
}))

# --- 单次 GTF 扫描提取全部窗口的 exon ---

window_prints <- paste(
  sprintf("if ($1==\"%s\" && $4<=%d && $5>=%d) print %d,a[1],$4,$5,$7,b[1];",
          win_dt$chrom, win_dt$end, win_dt$start, win_dt$win_id),
  collapse = " "
)
awk_prog <- sprintf(
  "BEGIN{FS=\"\\t\";OFS=\"\\t\"} $3==\"exon\" && match($9,/gene_name \"([^\"]+)\"/,a) && match($9,/gene_type \"([^\"]+)\"/,b){ %s }",
  window_prints
)
exon_all <- fread(
  cmd = sprintf("zcat %s | awk '%s'", shQuote(gtf_path), awk_prog),
  col.names = c("win_id", "gene_symbol", "start", "end", "strand", "gene_type")
)
cat("exon rows across windows:", nrow(exon_all), "\n")

# --- 单候选绘图函数（与 06-4 同风格） ---

plot_locus <- function(pts, gene_symbol, stage, trait, pp4, pp_h4, exon_dt) {
  lead_pts <- pts[is_lead == TRUE]
  locus_chr <- pts[!is.na(chrom), chrom[1]]
  pos_limits <- range(pts$pos_mb, na.rm = TRUE)

  regional_panel <- function(y_col, tag) {
    ggplot(pts, aes(pos_mb, .data[[y_col]])) +
      geom_point(aes(colour = ld_bin), size = 2.0, stroke = 0) +
      geom_point(data = lead_pts,
                 shape = 23, fill = col_lead, colour = "black",
                 size = 3.6, stroke = 0.5) +
      geom_text(data = lead_pts, aes(label = snp),
                fontface = "bold", size = 3.2, hjust = 1.12, vjust = -0.7) +
      annotate("text", x = -Inf, y = Inf, label = tag,
               hjust = -0.08, vjust = 1.5, size = 3.4) +
      scale_colour_manual(values = ld_cols, drop = FALSE, name = expression(r^2)) +
      scale_x_continuous(limits = pos_limits, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(x = NULL, y = NULL) +
      theme_06() +
      theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }

  gwas_panel <- regional_panel("gwas_logp", "GWAS")
  qtl_panel <- regional_panel("qtl_logp", "APA-QTL")

  if (nrow(exon_dt[gene_type == "protein_coding"]) > 0) {
    exon_dt <- exon_dt[gene_type == "protein_coding"]
  }
  gene_dt <- exon_dt[, .(start = min(start), end = max(end)),
                     by = gene_symbol][order(start)]
  lane_end <- numeric(0)
  gene_dt[, lane := 0L]
  for (i in seq_len(nrow(gene_dt))) {
    placed <- FALSE
    for (ln in seq_along(lane_end)) {
      if (gene_dt$start[i] - 30000 > lane_end[ln]) {
        gene_dt$lane[i] <- ln
        lane_end[ln] <- gene_dt$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      lane_end <- c(lane_end, gene_dt$end[i])
      gene_dt$lane[i] <- length(lane_end)
    }
  }
  exon_dt <- merge(exon_dt, gene_dt[, .(gene_symbol, lane)], by = "gene_symbol")

  char_kb <- 20
  case_gene <- gene_symbol
  gene_dt[, `:=`(mid = (start + end) / 2,
                 half_w = nchar(gene_symbol) * char_kb * 1000 / 2,
                 is_ensg = grepl("^ENSG", gene_symbol),
                 is_case = gene_symbol == case_gene,
                 show_label = FALSE)]
  setorder(gene_dt, is_ensg, mid)
  last_right <- -Inf
  for (i in which(!gene_dt$is_case)) {
    left <- gene_dt$mid[i] - gene_dt$half_w[i]
    if (left > last_right + 40000) {
      gene_dt$show_label[i] <- TRUE
      last_right <- max(last_right, gene_dt$mid[i] + gene_dt$half_w[i])
    }
  }
  case_i <- which(gene_dt$is_case)
  case_left <- gene_dt$mid[case_i] - gene_dt$half_w[case_i] - 20000
  case_right <- gene_dt$mid[case_i] + gene_dt$half_w[case_i] + 20000
  gene_dt$show_label[!gene_dt$is_case &
                       gene_dt$mid - gene_dt$half_w < case_right &
                       gene_dt$mid + gene_dt$half_w > case_left] <- FALSE
  gene_dt$show_label[case_i] <- TRUE

  gene_panel <- ggplot(gene_dt) +
    geom_segment(aes(x = start / 1e6, xend = end / 1e6, y = lane, yend = lane),
                 colour = col_gene, linewidth = 0.4) +
    geom_rect(data = exon_dt,
              aes(xmin = start / 1e6, xmax = end / 1e6,
                  ymin = lane - 0.16, ymax = lane + 0.16),
              colour = col_gene, fill = col_gene, linewidth = 0.3) +
    geom_text(data = gene_dt[show_label == TRUE],
              aes(x = mid / 1e6, y = lane + 0.45, label = gene_symbol),
              fontface = "italic", size = 2.9, colour = col_gene) +
    scale_x_continuous(sprintf("Genomic position on %s (Mb)", locus_chr),
                       limits = pos_limits, expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0.3, max(gene_dt$lane) + 0.9),
                       expand = expansion(mult = c(0, 0))) +
    labs(y = NULL) +
    theme_06() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  compare_plot <- ggplot(pts, aes(gwas_logp, qtl_logp)) +
    geom_point(aes(colour = ld_bin), size = 2.0, stroke = 0) +
    geom_point(data = lead_pts,
               shape = 23, fill = col_lead, colour = "black",
               size = 3.8, stroke = 0.5) +
    geom_text(data = lead_pts, aes(label = snp),
              fontface = "bold", size = 3.2, hjust = 1.12, vjust = -0.7) +
    {if (!is.na(pp_h4)) annotate(
      "text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.35,
      size = 3.0, colour = "grey20",
      label = sprintf("APA-GWAS PP4 = %.3f\nAPA-eQTL PP.H4 = %.3f", pp4, pp_h4)
    )} +
    scale_colour_manual(values = ld_cols, drop = FALSE, name = expression(r^2)) +
    labs(x = "GWAS −log10(P)", y = "APA-QTL −log10(P)") +
    theme_06() +
    theme(legend.position = c(0.97, 0.03),
          legend.justification = c(1, 0),
          legend.background = element_rect(colour = NA, fill = alpha("white", 0.75)),
          legend.title = element_text(size = font_sm, face = "italic"))

  right_col <- (gwas_panel / qtl_panel / gene_panel) +
    plot_layout(heights = c(1, 1, 0.5))

  (compare_plot | right_col) +
    plot_layout(widths = c(1, 1.15)) +
    plot_annotation(
      title = sprintf("%s | %s | %s",
                      gsub("_", " ", trait), stage, gene_symbol),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                              size = base_font + 2))
    )
}

# --- 批量导出 ---

n_ok <- 0L
n_fail <- 0L
for (i in seq_len(nrow(case_dt))) {
  this <- case_dt[i]
  stem <- if (run_mode == "interaction")
    sub("\\.json$", "", basename(this$locuscompare_json_file))
  else sprintf("%04d_%s_%s_%s",
               this$candidate_rank,
               gsub("[^A-Za-z0-9._-]", "_", this$gene_symbol),
               this$stage, this$hit_trait)
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  ok <- tryCatch({
    p <- plot_locus(
      pts = locus_list[[i]],
      gene_symbol = this$gene_symbol,
      stage = this$stage,
      trait = this$hit_trait,
      pp4 = as.numeric(this$hit_pp4),
      pp_h4 = this$max_PP_H4,
      exon_dt = exon_all[win_id == i]
    )
    ggsave(png_file, p, width = fig_w2, height = fig_h, dpi = png_dpi, bg = "white")
    ggsave(pdf_file, p, width = fig_w2, height = fig_h,
           device = grDevices::cairo_pdf, bg = "white")
    TRUE
  }, error = function(e) {
    cat("FAIL", stem, ":", conditionMessage(e), "\n")
    FALSE
  })
  if (ok) n_ok <- n_ok + 1L else n_fail <- n_fail + 1L
}

cat("locus figures written:", n_ok, "; failed:", n_fail, "\n")
cat("output dir:", out_dir, "\n")
