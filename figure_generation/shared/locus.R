# Regional locus plots (LD-coloured tracks + collision-free gene track); used by fig5_de.R and fig6.R.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(jsonlite); library(data.table) })

LOCUS_LD_LEVELS <- c("0.8-1.0", "0.6-0.8", "0.4-0.6", "0.2-0.4", "<0.2", "missing LD")
LOCUS_LD_COLS  <- c("0.8-1.0" = "red", "0.6-0.8" = "orange", "0.4-0.6" = "darkgreen",
                    "0.2-0.4" = "skyblue", "<0.2" = "blue4", "missing LD" = "#BDBDBD")
LOCUS_GENE_TABLE <- file.path(ROOT, "ref/gencode_v40_genes.tsv")
LOCUS_P5E8 <- -log10(5e-8)
LAB_MM_PER_CHAR <- 1.5     # italic 5.5 pt Arial, conservative glyph advance
LAB_PAD_MM      <- 2.0
LAB_GAP_MM      <- 1.5
N_LABEL_NEIGHBOURS <- 2

# --- readers ---
read_locus_json <- function(path) {
  j <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  pts <- as.data.frame(j$points, stringsAsFactors = FALSE)
  pts$ld_bin <- factor(pts$ld_bin, levels = LOCUS_LD_LEVELS)
  pts$is_lead <- as.logical(pts$is_lead)
  list(meta = j$meta, points = pts)
}

# target symbol and coordinates from the phenotype field (ENST|SYMBOL|chr:s-e|strand)
parse_target <- function(meta) {
  f <- strsplit(meta$phenotype, "\\|")[[1]]
  coords <- regmatches(f[3], regexec("chr[0-9XYM]+:([0-9]+)-([0-9]+)", f[3]))[[1]]
  list(symbol = f[2], start = as.numeric(coords[2]), end = as.numeric(coords[3]))
}

# genes in window from the cached GENCODE v40 table (protein-coding and lncRNA)
load_window_genes <- function(chrom, from, to) {
  g <- read.delim(LOCUS_GENE_TABLE, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("chrom", "start", "end", "strand", "gene", "gene_type"))
  keep <- g$chrom == chrom & g$end >= from & g$start <= to &
          g$gene_type %in% c("protein_coding", "lncRNA")
  g[keep, ]
}

# --- gene track ---
# label placed beside its own bar (right, else left); bars + labels packed into lanes
# with a greedy interval rule so no two extents share a lane
gene_lanes <- function(genes, target, lead_pos, window, pw_mm) {
  x0 <- window[1] / 1e6; x1 <- window[2] / 1e6; xr <- x1 - x0
  g <- genes[genes$gene_type == "protein_coding" | genes$gene == target$symbol, , drop = FALSE]
  g$xs <- g$start / 1e6; g$xe <- g$end / 1e6
  g$mid <- (g$xs + g$xe) / 2
  g$focal <- g$gene == target$symbol
  oth <- g[!g$focal, , drop = FALSE]
  oth <- oth[order(abs(oth$mid - lead_pos / 1e6)), , drop = FALSE]
  lab_set <- c(g$gene[g$focal], head(oth$gene, N_LABEL_NEIGHBOURS))
  g$lab <- g$gene %in% lab_set
  gap <- LAB_GAP_MM / pw_mm * xr
  g$wmb <- ifelse(g$lab, (nchar(g$gene) * LAB_MM_PER_CHAR + LAB_PAD_MM) / pw_mm * xr, 0)
  g$side <- "none"; g$lx <- NA_real_; g$hj <- NA_real_
  g$iv_s <- g$xs; g$iv_e <- g$xe
  for (i in which(g$lab)) {
    if (g$xe[i] + gap + g$wmb[i] <= x1) {
      g$side[i] <- "right"; g$lx[i] <- g$xe[i] + gap; g$hj[i] <- 0
      g$iv_e[i] <- g$xe[i] + gap + g$wmb[i]
    } else if (g$xs[i] - gap - g$wmb[i] >= x0) {
      g$side[i] <- "left"; g$lx[i] <- g$xs[i] - gap; g$hj[i] <- 1
      g$iv_s[i] <- g$xs[i] - gap - g$wmb[i]
    } else stop("no room for label ", g$gene[i], " in window")
  }
  pad <- 0.003 * xr
  ord <- order(g$iv_s)
  lane_end <- numeric(0); g$lane <- NA_integer_
  for (i in ord) {
    cand <- which(lane_end < g$iv_s[i] - pad)
    if (length(cand) == 0) { lane_end <- c(lane_end, g$iv_e[i]); g$lane[i] <- length(lane_end) }
    else { j <- cand[1]; lane_end[j] <- g$iv_e[i]; g$lane[i] <- j }
  }
  for (ln in unique(g$lane)) {
    s <- g[g$lane == ln, ]; s <- s[order(s$iv_s), ]
    if (nrow(s) > 1 && any(s$iv_s[-1] < s$iv_e[-nrow(s)])) stop("gene-track collision in lane ", ln)
  }
  g
}

gene_plot <- function(gl, target, x_lim, base_size = 6.5) {
  nl <- max(gl$lane)
  gl$xa <- ifelse(gl$strand == "+", gl$xs, gl$xe)
  gl$xb <- ifelse(gl$strand == "+", gl$xe, gl$xs)
  lab <- gl[gl$lab, , drop = FALSE]
  ggplot() +
    geom_segment(data = gl, aes(x = xa, xend = xb, y = lane, yend = lane, colour = focal, linewidth = focal),
                 arrow = arrow(type = "closed", length = unit(0.5, "mm"), ends = "last")) +
    geom_text(data = lab[!lab$focal, , drop = FALSE], aes(x = lx, y = lane, label = gene, hjust = hj),
              size = 5.5 / .pt, fontface = "italic", family = FONT, colour = MUTED) +
    geom_text(data = lab[lab$focal, , drop = FALSE], aes(x = lx, y = lane, label = gene, hjust = hj),
              size = 5.5 / .pt, fontface = "bold.italic", family = FONT, colour = INK) +
    scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = "#7F7F7F"), guide = "none") +
    scale_linewidth_manual(values = c(`TRUE` = 0.5, `FALSE` = 0.3), guide = "none") +
    scale_y_reverse(limits = c(nl + 0.6, 0.4), expand = c(0, 0)) +
    scale_x_continuous(limits = x_lim, breaks = pretty(x_lim, n = 4),
                       labels = function(b) sprintf("%.2f", b), expand = expansion(mult = c(0.01, 0.01))) +
    coord_cartesian(xlim = x_lim, clip = "off") +
    labs(x = sprintf("Chromosome %s position (Mb)", sub("chr", "", target$chrom)), y = NULL) +
    theme_nature(base_size = base_size) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          plot.margin = margin(0, 2, 1, 1, "mm"))
}

# --- LD-coloured track ---
# tr: list(df, ylab, pp4, threshold); df columns pos, logp, ld_bin, is_lead
track_plot <- function(tr, x_lim, lead_x, base_size = 6.5) {
  df <- tr$df; df$x <- df$pos / 1e6
  df <- df[order(df$ld_bin == "missing LD", df$ld_bin == "<0.2", decreasing = TRUE), ]
  ymax <- max(c(df$logp, if (isTRUE(tr$threshold)) LOCUS_P5E8), na.rm = TRUE) * 1.08
  ymax <- max(ymax, 1)
  lead <- df[df$is_lead, , drop = FALSE]
  seed <- data.frame(x = mean(x_lim), y = NA_real_, ld_bin = factor(LOCUS_LD_LEVELS, levels = LOCUS_LD_LEVELS))
  p <- ggplot() +
    geom_point(data = seed, aes(x, y, colour = ld_bin), alpha = 0, na.rm = TRUE) +
    geom_point(data = df[!df$is_lead, ], aes(x, logp, colour = ld_bin), size = 1.1, alpha = 0.82, stroke = 0) +
    geom_vline(xintercept = lead_x, colour = "#7F7F7F", linetype = "22", linewidth = 0.3)
  if (isTRUE(tr$threshold))
    p <- p + geom_hline(yintercept = LOCUS_P5E8, linetype = "dashed", colour = "#7F7F7F", linewidth = 0.25)
  if (nrow(lead) > 0)
    p <- p + geom_point(data = lead, aes(x, logp), shape = 21, fill = "white", colour = INK, size = 2.1, stroke = 0.6)
  if (!is.null(tr$pp4))
    p <- p + annotate("label", x = x_lim[2], y = ymax, label = tr$pp4, hjust = 1, vjust = 1.0,
                      size = 5.5 / .pt, colour = INK, family = FONT, fill = "white",
                      linewidth = 0, label.padding = unit(0.4, "mm"))
  p + scale_colour_manual(values = LOCUS_LD_COLS, drop = FALSE, name = "r\u00b2", limits = LOCUS_LD_LEVELS,
                          guide = "none") +
    scale_x_continuous(limits = x_lim, breaks = pretty(x_lim, n = 4), labels = function(b) sprintf("%.2f", b),
                       expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = tr$ylab) +
    theme_nature(base_size = base_size) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "none", plot.margin = margin(1, 2, 0, 1, "mm"))
}
