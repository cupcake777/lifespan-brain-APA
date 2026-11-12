#!/usr/bin/env Rscript
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")
# LocusZoom panels (eQTL, APA-QTL, gene track) for the selected coloc loci.
# LD from plink2 against the APA lead variant; gallery of up to 3 loci per stage.

# --- load packages ---
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

# --- config ---
D <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
OUT_DIR <- file.path(Sys.getenv("APA_ROOT"), "figures", "04_qtl")
GTF_FILE <- file.path(Sys.getenv("APA_ROOT"), "ref", "gencode.v40.annotation.gtf.gz")
GENO_PFX <- file.path(D, "pre_input", "genotype.%s")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STAGES <- c("prenatal", "postnatal")

LZ_FLANK_BP <- 200000
FIG_W15 <- 5.35
FIG_W2 <- 7.09
FIG_HMAX <- 9.45
BASE_FONT_SIZE <- 8
FONT_SM <- 7
TITLE_FONT_SIZE <- 9
PNG_DPI <- 600

LD_LEVELS <- c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0.0-0.2", "NA")
LD_COL <- c("0.8-1.0" = "#D62728", "0.6-0.8" = "#FF9E4A", "0.4-0.6" = "#5CB85C",
            "0.2-0.4" = "#7FC7E8", "0.0-0.2" = "#26456E", "NA" = "#BDBDBD")
NOMINAL_COLS <- c("phenotype_id", "variant_id", "af", "ma_count",
                  "pval_nominal", "slope", "slope_se")

# --- helpers ---

theme_sci <- function(base_size = BASE_FONT_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black", size = FONT_SM),
      axis.title = element_text(colour = "black", size = base_size),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      legend.title = element_text(size = FONT_SM),
      legend.text = element_text(size = FONT_SM),
      legend.key.size = grid::unit(0.30, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      strip.text = element_text(face = "bold", size = base_size, margin = margin(b = 1.5)),
      strip.background = element_blank(),
      plot.title = element_text(face = "bold", size = TITLE_FONT_SIZE, hjust = 0),
      plot.margin = margin(2, 3, 2, 2)
    )
}

save_plot <- function(plot, filename, width, height) {
  if (height > FIG_HMAX) {
    warning(sprintf("%s: height %.2f in exceeds the %.2f in page limit",
                    filename, height, FIG_HMAX))
  }
  ggsave(file.path(OUT_DIR, filename), plot, width = width, height = height, dpi = PNG_DPI, bg = "white")
}

split_variant <- function(v) {
  p <- strsplit(v, ":", fixed = TRUE)[[1]]
  list(chr = p[1], pos = as.integer(p[2]))
}

add_pos <- function(d, chr, from_bp, to_bp) {
  d <- copy(d)
  parts <- tstrsplit(d$variant_id, ":", fixed = TRUE)
  d[, `:=`(vchr = parts[[1]], pos = as.integer(parts[[2]]))]
  d[vchr == chr & pos >= from_bp & pos <= to_bp]
}

attach_ld <- function(d, ld, lead) {
  if (is.null(ld)) d[, r2 := NA_real_] else d[ld, r2 := i.r2, on = "variant_id"]
  d[variant_id == lead, r2 := 1]
  d[, ld_bin := as.character(cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                                 labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6",
                                            "0.6-0.8", "0.8-1.0"),
                                 right = FALSE))]
  d[is.na(ld_bin), ld_bin := "NA"]
  d[, ld_bin := factor(ld_bin, levels = LD_LEVELS)]
  d[]
}

assoc_panel <- function(dat, lead, title, lead_pos, other_pos, xlim, legend_title) {
  dat <- dat[order(ld_bin, decreasing = TRUE)]
  p <- ggplot(dat, aes(pos / 1e6, -log10(pmax(pval_nominal, 1e-300)))) +
    geom_vline(xintercept = lead_pos / 1e6, linetype = 3,
               colour = "#7B3294", linewidth = 0.35)
  if (!is.na(other_pos) && other_pos != lead_pos) {
    p <- p + geom_vline(xintercept = other_pos / 1e6, linetype = 3,
                        colour = "#0072B2", linewidth = 0.35)
  }
  p +
    geom_point(aes(fill = ld_bin), shape = 21, size = 1.0, stroke = 0.1,
               colour = "grey30", alpha = 0.9) +
    geom_point(data = dat[variant_id == lead], shape = 23, size = 1.9,
               fill = "#7B3294", colour = "black", stroke = 0.3,
               inherit.aes = TRUE, show.legend = FALSE) +
    scale_fill_manual(values = LD_COL, name = legend_title) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 1.5,
                                                   colour = "grey30", alpha = 1))) +
    scale_x_continuous(limits = xlim / 1e6) +
    scale_y_continuous(expand = expansion(c(0.02, 0.12))) +
    labs(x = NULL, y = expression(-log[10](italic(p))), title = title) +
    theme_sci() +
    theme(plot.title = element_text(size = BASE_FONT_SIZE, face = "plain"))
}

# --- inputs ---

selected <- fread(file.path(D, "coloc_locuszoom_selection.tsv"))
coloc_all <- fread(file.path(D, "coloc_all_stages.tsv"))
apa_map <- unique(coloc_all[, .(phenotype_id = apa_pheno, ensg0)])

exon_dt <- fread(cmd = paste0(
  "zcat ", GTF_FILE, " | awk -F'\\t' '$3==\"exon\" {",
  "match($9,/gene_name \"[^\"]+\"/); s=substr($9,RSTART+11,RLENGTH-12);",
  "print $1\"\\t\"$4\"\\t\"$5\"\\t\"$7\"\\t\"s}'"),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("chr", "start", "end", "strand", "gene_name"))

# --- render loci ---

lz_out <- list()
for (st in STAGES) {
  sel_st <- selected[stage == st]
  if (nrow(sel_st) == 0) next
  dt_e <- fread(file.path(D, sprintf("nominal_eqtl_%s.tsv.gz", st)), select = NOMINAL_COLS)
  dt_e[, ensg0 := sub("\\..*", "", phenotype_id)]
  dt_e <- dt_e[is.finite(slope) & is.finite(slope_se) & slope_se > 0]

  dt_a <- fread(file.path(D, sprintf("nominal_apa_%s.tsv.gz", st)), select = NOMINAL_COLS)
  dt_a <- merge(dt_a, apa_map, by = "phenotype_id", all.x = FALSE)
  dt_a <- dt_a[is.finite(slope) & is.finite(slope_se) & slope_se > 0]

  nom <- list(eqtl = dt_e, apa = dt_a)

  for (i in seq_len(nrow(sel_st))) {
    s <- as.list(sel_st[i])
    r <- tryCatch({
      stage <- as.character(s$stage)
      lead <- s$apa_lead_snp
      lv <- split_variant(lead)
      ev <- split_variant(s$eqtl_lead_snp)

      from_bp <- max(1L, min(lv$pos, ev$pos) - LZ_FLANK_BP)
      to_bp <- max(lv$pos, ev$pos) + LZ_FLANK_BP
      xlim <- c(from_bp, to_bp)

      d_e <- add_pos(nom$eqtl[phenotype_id == s$eqtl_pheno], lv$chr, from_bp, to_bp)
      d_a <- add_pos(nom$apa[phenotype_id == s$apa_pheno], lv$chr, from_bp, to_bp)
      if (nrow(d_e) < 20 || nrow(d_a) < 20) {
        NULL
      } else {
        # LD of the window with the APA lead variant
        pfx <- tempfile(pattern = "ld")
        args <- c("--pfile", sprintf(GENO_PFX, stage),
                  "--chr", sub("^chr", "", lv$chr),
                  "--from-bp", format(from_bp, scientific = FALSE),
                  "--to-bp", format(to_bp, scientific = FALSE),
                  "--ld-snp", lead,
                  "--r2-unphased", "cols=chrom,pos,id",
                  "--ld-window-kb", "5000", "--ld-window-r2", "0",
                  "--threads", "2", "--memory", "8000",
                  "--out", pfx)
        status <- suppressWarnings(system2("plink2", args, stdout = FALSE, stderr = FALSE))
        f <- paste0(pfx, ".vcor")
        ok <- identical(status, 0L) && file.exists(f)
        ld <- if (ok) fread(f) else NULL
        unlink(list.files(dirname(pfx), paste0(basename(pfx), "*"), full.names = TRUE))
        if (is.null(ld) || nrow(ld) == 0) {
          ld <- NULL
        } else {
          setnames(ld, c("ID_B", "UNPHASED_R2"), c("variant_id", "r2"), skip_absent = TRUE)
          ld <- unique(ld[, .(variant_id, r2)], by = "variant_id")
        }

        d_e <- attach_ld(d_e, ld, lead)
        d_a <- attach_ld(d_a, ld, lead)

        legend_title <- bquote(r^2 ~ "with APA lead")
        p_e <- assoc_panel(d_e, lead, sprintf("eQTL   %s   (%s)", s$gene_name, s$eqtl_pheno),
                           lv$pos, ev$pos, xlim, legend_title)
        p_a <- assoc_panel(d_a, lead, sprintf("APA-QTL   %s   %s", s$gene_name, s$transcript_id),
                           lv$pos, ev$pos, xlim, legend_title) +
          theme(legend.position = "none")

        # gene track: genes with exons in the window, packed into rows; locus gene highlighted
        reg <- exon_dt[chr == lv$chr & end >= xlim[1] & start <= xlim[2]]
        if (nrow(reg) == 0) {
          p_g <- ggplot() + scale_x_continuous(limits = xlim / 1e6) +
            labs(x = sprintf("%s position (Mb)", lv$chr), y = NULL) +
            theme_sci() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.line.y = element_blank())
        } else {
          span <- reg[, .(gs = min(start), ge = max(end)), by = gene_name]
          setorder(span, gs)
          span[, row := 0L]
          row_end <- numeric(0)
          pad <- diff(xlim) * 0.10
          for (k in seq_len(nrow(span))) {
            slot <- which(row_end < span$gs[k] - pad)
            row_k <- if (length(slot)) slot[1] else length(row_end) + 1L
            span$row[k] <- row_k
            row_end[row_k] <- span$ge[k]
          }
          span[, is_hl := gene_name %in% s$gene_name]
          span[, mid := (gs + ge) / 2]
          label_span <- span[is_hl == TRUE]
          if (nrow(label_span) == 0) label_span <- span[order(abs(mid - mean(xlim)))][1]
          reg <- merge(reg, span[, .(gene_name, row, is_hl)], by = "gene_name")

          p_g <- ggplot() +
            geom_segment(data = span, aes(x = pmax(gs, xlim[1]) / 1e6,
                                          xend = pmin(ge, xlim[2]) / 1e6,
                                          y = row, yend = row, colour = is_hl),
                         linewidth = 0.35) +
            geom_rect(data = reg, aes(xmin = pmax(start, xlim[1]) / 1e6,
                                      xmax = pmin(end, xlim[2]) / 1e6,
                                      ymin = row - 0.20, ymax = row + 0.20,
                                      fill = is_hl)) +
            geom_text_repel(data = label_span,
                            aes(x = (pmax(gs, xlim[1]) + pmin(ge, xlim[2])) / 2 / 1e6,
                                y = row + 0.52, label = gene_name, colour = is_hl,
                                fontface = ifelse(is_hl, "bold.italic", "italic")),
                            size = FONT_SM / .pt, direction = "y", seed = 1,
                            box.padding = 0.15, point.padding = 0.03,
                            min.segment.length = 0, max.overlaps = Inf,
                            segment.size = 0.2, show.legend = FALSE) +
            scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey35"),
                                guide = "none") +
            scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey35"),
                              guide = "none") +
            scale_x_continuous(limits = xlim / 1e6) +
            scale_y_reverse(expand = expansion(c(0.38, 0.40))) +
            labs(x = sprintf("%s position (Mb)", lv$chr), y = NULL) +
            theme_sci() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                  axis.line.y = element_blank(),
                  plot.margin = margin(2, 8, 2, 2)) +
            coord_cartesian(clip = "off")
        }

        fig <- (p_e / p_a / p_g) +
          plot_layout(heights = c(1, 1, 0.72), guides = "collect") +
          plot_annotation(
            title = sprintf("%s | %s | %s", stage, s$gene_name, s$lz_tag),
            theme = theme(plot.title = element_text(face = "bold", size = TITLE_FONT_SIZE))) &
          theme(legend.position = "right")

        fname <- sprintf("locuszoom_%s_%s_%s.png", stage,
                         gsub("[^A-Za-z0-9]", "-", s$gene_name), s$lz_tag)
        save_plot(fig, fname, width = FIG_W15, height = 4.4)
        list(file = fname, apa_panel = p_a, tag = s$lz_tag, sel = s)
      }
    }, error = function(e) {
      message("locuszoom failed for ", s$gene_name, ": ", conditionMessage(e))
      NULL
    })
    if (!is.null(r)) lz_out[[length(lz_out) + 1]] <- r
  }
  rm(nom)
  invisible(gc(verbose = FALSE))
}

# --- figure 4: LocusZoom gallery ---

if (length(lz_out) > 0) {
  lz_dt <- data.table(i = seq_along(lz_out),
                      stage = vapply(lz_out, function(r) as.character(r$sel$stage), ""),
                      tag = vapply(lz_out, function(r) r$tag, ""))
  lz_dt[, tag_rank := match(tag, c("shared-H4", "independent-H3",
                                   "discordant-direction", "stage-switch"))]
  pick <- lz_dt[order(tag_rank, stage)][, head(.SD, 3), by = stage][order(tag_rank, stage)]$i
  gallery <- patchwork::wrap_plots(lapply(lz_out[pick], function(r) {
    r$apa_panel +
      labs(title = sprintf("%s | %s | %s", r$sel$stage, r$sel$gene_name, r$tag)) +
      theme(plot.title = element_text(size = FONT_SM, face = "bold"),
            legend.position = "none")
  }), ncol = 2)
  save_plot(gallery, "fig4_locuszoom_gallery.png", width = FIG_W2, height = 6.0)
}
