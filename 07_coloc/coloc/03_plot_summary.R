# Colocalization summary: PP4 overview and locuscompare plots
# Selects qualifying loci, computes PLINK LD around each lead SNP, draws the figures

suppressMessages({
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(readr)
  library(stringr)
  library(yaml)
})

# --- paths ---
root <- Sys.getenv("APA_ROOT")
data_dir <- file.path(root, "data", "07_coloc")
fig_dir <- file.path(root, "figures", "07_coloc")
coloc_file <- file.path(data_dir, "coloc_run", "all_colocalization.csv")
rds_dir <- file.path(data_dir, "coloc_run", "rds")
ld_dir <- file.path(data_dir, "plink_ld")

global_cfg <- yaml::read_yaml(file.path(data_dir, "global.yml"))
plink_bin <- Sys.which("plink")
plink_bfile <- path.expand(str_replace(global_cfg$plink_hg38, fixed("$HOME"), Sys.getenv("HOME")))
gene_anno <- path.expand(str_replace(global_cfg$gene_anno, fixed("$HOME"), Sys.getenv("HOME")))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ld_dir, recursive = TRUE, showWarnings = FALSE)

# --- config ---
stage_levels   <- c("prenatal", "postnatal", "interaction")
stage_labels   <- c(prenatal = "Prenatal", postnatal = "Postnatal", interaction = "Interaction")
stage_colors   <- c(prenatal = "#E64B35", postnatal = "#4DBBD5", interaction = "#00A087")

plink_threads  <- 4
figure_dpi     <- 320
pp4_cutoff     <- 0.80   # PIP > 0.8

ld_bin_colors <- c(
  "0.8–1.0" = "#D73027",
  "0.6–0.8" = "#FC8D59",
  "0.4–0.6" = "#FEE090",
  "0.2–0.4" = "#91BFDB",
  "0.0–0.2" = "#4575B4",
  "Missing"  = "#C8C8C8"
)
ld_bin_levels <- c("0.8–1.0", "0.6–0.8", "0.4–0.6", "0.2–0.4", "0.0–0.2", "Missing")

theme_locuscompare <- theme_classic(base_size = 10) +
  theme(
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9, face = "bold", colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 8.5, hjust = 0, colour = "black"),
    plot.subtitle = element_text(size = 7, colour = "grey40"),
    plot.margin = margin(4, 4, 4, 4)
  )

# --- helpers ---
build_locuscompare_data <- function(meta, rds_dir, plink_keep, plink_bfile, plink_bin, ld_dir, gene_models) {
  # RDS file
  gwas_re  <- str_replace_all(as.character(meta$GWAS_ID), "([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1")
  qtl_re   <- str_replace_all(as.character(meta$QTL_ID),  "([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1")
  locus_re <- str_replace_all(meta$Locus,                  "([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1")
  gene_re  <- str_replace_all(meta$gene,                   "([.^$|()\\[\\]{}*+?\\\\])", "\\\\\\1")
  rds_pattern <- paste0("^", gwas_re, "_", qtl_re, "_.*_", locus_re, "_", gene_re, "_coloc\\.rds$")
  rds_hits <- list.files(rds_dir, pattern = rds_pattern, full.names = TRUE)
  if (length(rds_hits) == 0) return(NULL)

  obj <- readRDS(rds_hits[[1]])
  md <- as_tibble(obj$merged_data) %>%
    filter(!is.na(snp), !is.na(P.gwas), !is.na(P.qtl), P.gwas > 0, P.qtl > 0) %>%
    distinct(snp, .keep_all = TRUE)
  if (nrow(md) == 0) return(NULL)

  lead_snp <- obj$lead_snp
  if (is.null(lead_snp) || length(lead_snp) == 0 || !(lead_snp %in% md$snp)) {
    lead_snp <- md %>% arrange(pmin(P.gwas, P.qtl)) %>% slice(1) %>% pull(snp)
  }

  ancestry <- str_extract(as.character(meta$GWAS_ID), "^[^_]+")
  keep_file <- plink_keep[[ancestry]]

  pidx <- meta$plot_index
  safe_stem <- paste(
    sprintf("%02d", pidx),
    as.character(meta$GWAS_ID), as.character(meta$QTL_ID),
    meta$gene, meta$Locus, sep = "_"
  ) %>% str_replace_all("[^A-Za-z0-9._-]+", "_")

  # PLINK LD
  ld_for_plot <- data.frame(SNP_A = lead_snp, SNP_B = md$snp, R2 = NA_real_)
  if (plink_bin != "" && file.exists(paste0(plink_bfile, ".bed")) &&
      !is.null(keep_file) && !is.na(keep_file) && file.exists(keep_file)) {
    ld_prefix <- file.path(ld_dir, safe_stem)
    snp_file  <- paste0(ld_prefix, ".snps")
    write_lines(unique(md$snp), snp_file)
    plink_args <- c(
      "--bfile", plink_bfile, "--allow-extra-chr", "--keep", keep_file,
      "--extract", snp_file, "--r2", "--ld-snp", lead_snp,
      "--ld-window-kb", "1000", "--ld-window", "999999", "--ld-window-r2", "0",
      "--threads", as.character(plink_threads), "--out", ld_prefix
    )
    plink_status <- system2(plink_bin, args = plink_args)
    ld_file <- paste0(ld_prefix, ".ld")
    if (identical(as.integer(plink_status), 0L) && file.exists(ld_file) && file.info(ld_file)$size > 0) {
      ld_raw <- read_table(ld_file, show_col_types = FALSE, progress = FALSE,
                           col_types = cols(.default = col_character()))
      r2_col    <- intersect(names(ld_raw), c("R2", "UNPHASED_R2"))
      snp_a_col <- intersect(names(ld_raw), c("SNP_A", "SNP1", "SNP"))
      snp_b_col <- intersect(names(ld_raw), c("SNP_B", "SNP2"))
      if (length(r2_col) > 0 && length(snp_a_col) > 0 && length(snp_b_col) > 0) {
        ld_for_plot <- ld_raw %>%
          transmute(
            SNP_A = lead_snp,
            SNP_B = if_else(.data[[snp_a_col[[1]]]] == lead_snp, .data[[snp_b_col[[1]]]], .data[[snp_a_col[[1]]]]),
            R2 = as.numeric(.data[[r2_col[[1]]]])
          ) %>%
          filter(!is.na(SNP_B), !is.na(R2)) %>%
          bind_rows(data.frame(SNP_A = lead_snp, SNP_B = lead_snp, R2 = 1)) %>%
          distinct(SNP_A, SNP_B, .keep_all = TRUE)
      }
    }
  }

  locus_data <- md %>%
    transmute(
      snp = snp,
      gwas_logp = -log10(P.gwas),
      qtl_logp  = -log10(P.qtl),
      chr = paste0("chr", str_remove(as.character(CHR.gwas), "^chr")),
      pos = as.numeric(POS.gwas)
    ) %>%
    left_join(ld_for_plot %>% select(snp = SNP_B, R2), by = "snp") %>%
    mutate(
      R2 = if_else(snp == lead_snp, 1, R2),
      ld_bin = case_when(
        snp == lead_snp   ~ NA_character_,
        is.na(R2)         ~ "Missing",
        R2 >= 0.8         ~ "0.8–1.0",
        R2 >= 0.6         ~ "0.6–0.8",
        R2 >= 0.4         ~ "0.4–0.6",
        R2 >= 0.2         ~ "0.2–0.4",
        TRUE              ~ "0.0–0.2"
      ),
      ld_bin = factor(ld_bin, levels = ld_bin_levels)
    ) %>%
    filter(!is.na(gwas_logp), !is.na(qtl_logp))

  # Gene track
  locus_chr  <- locus_data$chr[[1]]
  locus_start <- min(locus_data$pos, na.rm = TRUE)
  locus_end   <- max(locus_data$pos, na.rm = TRUE)

  track_genes <- gene_models %>%
    filter(chr == locus_chr, start <= locus_end, end >= locus_start) %>%
    mutate(
      plot_start = pmax(start, locus_start) / 1e6,
      plot_end   = pmin(end, locus_end) / 1e6,
      midpoint   = (plot_start + plot_end) / 2,
      is_target  = gene_name == meta$gene,
      is_named   = !str_detect(gene_name, "^ENSG"),
      is_pc      = gene_type == "protein_coding",
      span       = plot_end - plot_start
    ) %>%
    arrange(desc(is_target), desc(is_pc), desc(is_named), desc(span), start, end, gene_name) %>%
    slice_head(n = 18) %>%
    arrange(start, end, gene_name)

  if (nrow(track_genes) > 0) {
    row_ends <- numeric()
    track_rows <- integer(nrow(track_genes))
    for (j in seq_len(nrow(track_genes))) {
      available <- which(row_ends < track_genes$start[[j]])
      row_id <- if (length(available) == 0) length(row_ends) + 1 else available[[1]]
      row_ends[[row_id]] <- track_genes$end[[j]]
      track_rows[[j]] <- row_id
    }
    track_genes <- track_genes %>% mutate(track_row = track_rows)
  }

  label_genes <- track_genes %>%
    arrange(desc(is_target), desc(is_pc), desc(is_named), desc(span), gene_name) %>%
    slice_head(n = 6)

  list(
    data = locus_data,
    lead_snp = lead_snp,
    track_genes = track_genes,
    label_genes = label_genes,
    locus_chr = locus_chr,
    locus_start = locus_start,
    locus_end = locus_end,
    pp4 = meta$PP4,
    stage = as.character(meta$QTL_ID),
    gene = meta$gene,
    gwas_id = as.character(meta$GWAS_ID),
    locus = meta$Locus
  )
}

plot_locuscompare <- function(lc) {
  d <- lc$data
  lead_label <- d %>% filter(snp == lc$lead_snp)
  span_x <- range(d$gwas_logp, na.rm = TRUE)
  span_y <- range(d$qtl_logp, na.rm = TRUE)

  # scatter panel
  p_scatter <- ggplot(d, aes(x = gwas_logp, y = qtl_logp)) +
    geom_point(aes(fill = ld_bin), shape = 21, colour = "grey30", stroke = 0.2, size = 1.8, alpha = 0.85) +
    geom_point(data = lead_label, shape = 23, fill = stage_colors[lc$stage], colour = "black", size = 2.8, stroke = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_text_repel(
      data = lead_label, aes(label = snp),
      size = 2.5, fontface = "bold", colour = "black",
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40",
      box.padding = 0.4, point.padding = 0.2, nudge_x = 0.3, nudge_y = 0.3,
      max.overlaps = Inf, force = 2, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = ld_bin_colors, breaks = ld_bin_levels, limits = ld_bin_levels,
      name = expression(LD~r^2), drop = FALSE,
      guide = guide_legend(override.aes = list(shape = 21, size = 2.5, stroke = 0.15))
    ) +
    coord_cartesian(xlim = span_x + c(-0.1, 0.5) * diff(span_x),
                    ylim = span_y + c(-0.1, 0.5) * diff(span_y)) +
    labs(
      title = sprintf("%s · %s · %s", lc$gwas_id, lc$gene, lc$locus),
      subtitle = sprintf("PP4 = %.3f · %s", lc$pp4, stage_labels[lc$stage]),
      x = expression(-log[10](GWAS~italic(P))),
      y = expression(-log[10](APA-QTL~italic(P)))
    ) +
    theme_locuscompare +
    theme(legend.position = "right")

  # gene track panel
  locus_span_mb <- c(lc$locus_start, lc$locus_end) / 1e6
  if (nrow(lc$track_genes) == 0) {
    p_gene <- ggplot(data.frame(x = locus_span_mb, y = 1)) +
      geom_blank(aes(x = x, y = y)) +
      scale_x_continuous(limits = locus_span_mb, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_continuous(limits = c(0.5, 1.5)) +
      labs(x = paste0(lc$locus_chr, " (Mb)"), y = "Genes") +
      theme_locuscompare +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    p_gene <- ggplot(lc$track_genes) +
      geom_segment(aes(x = plot_start, xend = plot_end, y = track_row, yend = track_row),
                   linewidth = 0.6, colour = "grey25") +
      geom_segment(data = lc$track_genes %>% filter(strand == "+"),
                   aes(x = pmax(plot_start, plot_end - 0.01), xend = plot_end, y = track_row, yend = track_row),
                   arrow = arrow(length = unit(0.04, "inches")), linewidth = 0.3, colour = "grey25") +
      geom_segment(data = lc$track_genes %>% filter(strand == "-"),
                   aes(x = pmin(plot_end, plot_start + 0.01), xend = plot_start, y = track_row, yend = track_row),
                   arrow = arrow(length = unit(0.04, "inches")), linewidth = 0.3, colour = "grey25") +
      geom_text_repel(
        data = lc$label_genes,
        aes(x = midpoint, y = track_row, label = gene_name),
        size = 2.2, fontface = "italic", colour = "grey15",
        direction = "both", min.segment.length = 0, segment.size = 0.15, segment.colour = "grey50",
        box.padding = 0.4, point.padding = 0.15, nudge_y = 0.35, max.overlaps = Inf,
        force = 3, force_pull = 0.5, seed = 7, show.legend = FALSE
      ) +
      scale_x_continuous(limits = locus_span_mb, expand = expansion(mult = c(0.01, 0.01))) +
      scale_y_reverse(limits = c(max(lc$track_genes$track_row) + 0.6, -0.2),
                      expand = expansion(mult = c(0.02, 0.05))) +
      labs(x = paste0(lc$locus_chr, " (Mb)"), y = "Genes") +
      theme_locuscompare +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
  }

  p_scatter / p_gene + plot_layout(heights = c(3, 1.2), guides = "collect") &
    theme(plot.margin = margin(3, 3, 3, 3))
}

# --- load coloc results ---
coloc <- read_csv(coloc_file, show_col_types = FALSE) %>%
  mutate(
    GWAS_ID   = as.character(GWAS_ID),
    QTL_ID    = as.character(QTL_ID),
    Locus     = as.character(Locus),
    Phenotype = as.character(Phenotype),
    PP4       = as.numeric(PP4),
    n_snps    = as.numeric(n_snps)
  )

pheno_parts <- str_split_fixed(coloc$Phenotype, "\\|", 4)
coloc <- coloc %>%
  mutate(
    ancestry = str_extract(GWAS_ID, "^[^_]+"),
    disease  = str_remove(GWAS_ID, "^[^_]+_"),
    gene     = pheno_parts[, 2],
    gene     = if_else(is.na(gene) | gene == "", Phenotype, gene),
    stage    = factor(QTL_ID, levels = stage_levels)
  ) %>%
  filter(!is.na(stage), !is.na(PP4), !is.na(n_snps))

# --- qualifying loci ---
# Rule: PP4 > 0.8  OR  (PP4 >= 0.75 AND stage ∈ {prenatal, interaction})
sig_all <- coloc %>%
  arrange(GWAS_ID, QTL_ID, gene, Locus, desc(PP4), desc(n_snps)) %>%
  distinct(GWAS_ID, QTL_ID, gene, Locus, .keep_all = TRUE)

qualifying <- sig_all %>%
  filter(PP4 > pp4_cutoff | (PP4 >= 0.75 & QTL_ID %in% c("prenatal", "interaction")))


write_csv(qualifying, file.path(data_dir, "qualifying_loci.csv"))

# --- summary overview ---
# Use all PP4 >= 0.75 for the overview (broader context)
overview <- sig_all %>% filter(PP4 >= 0.75)

# Order genes by max PP4
gene_order <- overview %>%
  group_by(gene) %>%
  summarise(max_PP4 = max(PP4), .groups = "drop") %>%
  arrange(desc(max_PP4)) %>%
  pull(gene)

overview <- overview %>%
  mutate(gene = factor(gene, levels = rev(gene_order)))

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text.x = element_text(size = 9, colour = "black", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8, colour = "black", margin = margin(r = 3)),
    axis.title = element_text(size = 10, face = "bold", colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.3, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "#f1f5f9", colour = NA),
    strip.text = element_text(face = "bold", size = 9, colour = "black"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0, colour = "black"),
    plot.margin = margin(10, 16, 10, 16)
  )

n_genes <- n_distinct(overview$gene)
overview_height <- max(5, 0.22 * n_genes + 2.0)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

p_overview <- ggplot(overview, aes(x = stage, y = gene)) +
  geom_point(aes(size = PP4, fill = stage), shape = 21, colour = "grey30", stroke = 0.3, alpha = 0.88) +
  geom_text(data = overview %>% filter(PP4 > 0.95),
            aes(label = sprintf("%.2f", PP4)), size = 2.5, fontface = "bold", colour = "black",
            nudge_y = 0.25, show.legend = FALSE) +
  scale_fill_manual(values = stage_colors, name = "Stage", drop = FALSE) +
  scale_size_continuous(range = c(1.5, 5.5), name = "PP4", breaks = c(0.8, 0.9, 1.0), limits = c(0.75, 1.0)) +
  scale_x_discrete(labels = stage_labels, expand = c(0.08, 0.08)) +
  scale_y_discrete(expand = c(0.02, 0.02)) +
  labs(
    title = "Colocalization Landscape: PP4 ≥ 0.75",
    subtitle = sprintf("%d unique locus-level results across %d genes", nrow(overview), n_genes),
    x = NULL, y = NULL
  ) +
  theme_pub +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "summary_overview.pdf"), p_overview,
       width = 8.5, height = overview_height, units = "in", device = pdf_device)
ggsave(file.path(fig_dir, "summary_overview.png"), p_overview,
       width = 8.5, height = overview_height, units = "in", dpi = figure_dpi)



# --- gene annotation ---
gene_models <- read_tsv(
  pipe(paste("zcat", shQuote(gene_anno), "| awk -F '\\t' '$3==\"gene\" {print}'")),
  col_names = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
  col_types = cols(
    chr = col_character(), source = col_character(), feature = col_character(),
    start = col_double(), end = col_double(), score = col_character(),
    strand = col_character(), frame = col_character(), attribute = col_character()
  ),
  progress = FALSE
) %>%
  mutate(
    gene_name = str_match(attribute, "gene_name \\\"([^\\\"]+)\\\"")[, 2],
    gene_type = str_match(attribute, "gene_type \\\"([^\\\"]+)\\\"")[, 2]
  ) %>%
  select(chr, start, end, strand, gene_name, gene_type) %>%
  filter(!is.na(gene_name))

# --- plink keep files ---
raw_keeps <- global_cfg$plink_keep
plink_keep <- list()
for (anc in names(raw_keeps)) {
  keep_path <- path.expand(str_replace(raw_keeps[[anc]], fixed("$HOME"), Sys.getenv("HOME")))
  if (!file.exists(keep_path)) {
    plink_keep[[anc]] <- NA_character_
  } else {
    keep_tbl <- read_table(keep_path, col_names = FALSE, show_col_types = FALSE, progress = FALSE,
                           col_types = cols(.default = col_character()))
    if (ncol(keep_tbl) >= 2) {
      plink_keep[[anc]] <- keep_path
    } else {
      out_keep <- file.path(ld_dir, paste0("keep_", anc, "_plink12.keep"))
      data.frame(FID = "0", IID = keep_tbl[[1]]) %>% write_tsv(out_keep, col_names = FALSE)
      plink_keep[[anc]] <- out_keep
    }
  }
}

# --- locuscompare plots ---
qualifying$plot_index <- seq_len(nrow(qualifying))

pdf_device(file = file.path(fig_dir, "locuscompare_all_qualifying.pdf"), width = 8, height = 6.5)
for (i in seq_len(nrow(qualifying))) {
  meta <- qualifying[i, ]
  lc <- build_locuscompare_data(meta, rds_dir, plink_keep, plink_bfile, plink_bin, ld_dir, gene_models)
  if (is.null(lc)) next
  p <- plot_locuscompare(lc)
  print(p)
}
dev.off()

# --- locuscompare highlights ---
highlight_dir <- file.path(fig_dir, "locuscompare_highlights")
dir.create(highlight_dir, recursive = TRUE, showWarnings = FALSE)

highlights <- qualifying %>% filter(QTL_ID %in% c("prenatal", "interaction") | PP4 > 0.95)
highlights$plot_index <- seq_len(nrow(highlights))

for (i in seq_len(nrow(highlights))) {
  meta <- highlights[i, ]
  lc <- build_locuscompare_data(meta, rds_dir, plink_keep, plink_bfile, plink_bin, ld_dir, gene_models)
  if (is.null(lc)) next

  safe_name <- paste(as.character(meta$GWAS_ID), meta$QTL_ID, meta$gene, meta$Locus, sep = "_") %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_")

  p <- plot_locuscompare(lc)
  ggsave(file.path(highlight_dir, paste0(safe_name, ".pdf")), p,
         width = 8, height = 6.5, units = "in", device = pdf_device)
  ggsave(file.path(highlight_dir, paste0(safe_name, ".png")), p,
         width = 8, height = 6.5, units = "in", dpi = figure_dpi)
}
