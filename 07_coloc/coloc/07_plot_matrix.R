# Gene x GWAS colocalization dot matrix for the PP4 >= 0.7 loci
# Rows are genes ordered by chromosome and position; columns are GWAS traits, EAS | EUR

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# --- paths ---
root <- Sys.getenv("APA_ROOT")
data_dir <- file.path(root, "data", "07_coloc")
fig_dir <- file.path(root, "figures", "07_coloc")
matrix_file <- file.path(data_dir, "matrix_data.tsv")  # not in the path map; no producer in repo

# --- config ---
stage_levels <- c("prenatal", "postnatal", "interaction")
stage_cols   <- c(prenatal = "#E64B35", postnatal = "#4DBBD5", interaction = "#00A087")
ink          <- "#20252A"

gwas_order <- c("EAS_AIS", "EAS_PD", "EAS_SCZ",
                "EUR_AD", "EUR_ADHD", "EUR_BIP", "EUR_INT",
                "EUR_MDD", "EUR_PD", "EUR_PTSD", "EUR_SCZ")
ancestry_group <- ifelse(startsWith(gwas_order, "EAS"), "EAS", "EUR")

figure_width  <- 7.2
figure_height <- 5.6
dpi           <- 600
base_size     <- 7.5
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- load matrix ---
d <- read_tsv(matrix_file, show_col_types = FALSE) %>%
  mutate(
    PP4   = as.numeric(PP4),
    stage = factor(stage, levels = stage_levels),
    GWAS  = factor(GWAS, levels = gwas_order)
  ) %>%
  filter(!is.na(PP4), !is.na(stage))

# gene order: by chromosome then position, key = gene|chr|pos
gene_chr <- d %>% distinct(gene, chrom, pos)
gene_key_order <- gene_chr %>%
  mutate(
    chr_num = as.integer(gsub("^chr", "", chrom)),
    key     = paste(gene, chr_num, pos, sep = "|")
  ) %>%
  arrange(chr_num, pos) %>%
  pull(key)
d <- d %>%
  mutate(
    chr_num = as.integer(gsub("^chr", "", chrom)),
    gene_key = paste(gene, chr_num, pos, sep = "|"),
    gene_key = factor(gene_key, levels = gene_key_order)
  )

# Row bands: odd/even chromosome backgrounds (per unique gene key)
gene_meta <- gene_chr %>%
  mutate(chr_num = as.integer(gsub("^chr", "", chrom))) %>%
  arrange(chr_num, pos) %>%
  mutate(
    key     = paste(gene, chr_num, pos, sep = "|"),
    band    = ifelse(chr_num %% 2 == 1, "odd", "even")
  )
row_band <- data.frame(
  gene_key = factor(gene_meta$key, levels = gene_key_order),
  band     = gene_meta$band
)

# Column bands: EAS | EUR ancestry groups
grp_boundaries <- data.frame(
  xmin = c(0.5, sum(ancestry_group == "EAS") + 0.5),
  xmax = c(sum(ancestry_group == "EAS") + 0.5, length(gwas_order) + 0.5),
  grp  = c("EAS", "EUR")
)

# jitter dots within a gene x GWAS cell by stage (max 3 stages per cell)
stage_offset <- c(prenatal = -0.22, postnatal = 0, interaction = 0.22)
d <- d %>% mutate(
  gwas_x = as.integer(GWAS) + stage_offset[as.character(stage)]
)

# --- plot ---
p <- ggplot() +
  # column ancestry bands (grey tones, drawn first = background)
  geom_rect(
    data = grp_boundaries,
    aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = length(gene_key_order) + 0.5,
        fill = grp),
    inherit.aes = FALSE
  ) +
  # row chromosome bands (white / very light grey stripes)
  geom_rect(
    data = row_band,
    aes(xmin = 0.5, xmax = length(gwas_order) + 0.5,
        ymin = as.integer(gene_key) - 0.5, ymax = as.integer(gene_key) + 0.5,
        fill = band),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(EAS = "#EEF1F5", EUR = "#F7F5F0",
               odd = "white", even = "#F3F5F8"),
    guide = "none"
  ) +
  # dots
  geom_point(
    data = d,
    aes(x = gwas_x, y = gene_key, color = stage, size = PP4),
    alpha = 0.9
  ) +
  scale_color_manual(
    values = stage_cols,
    labels = c(prenatal = "Prenatal", postnatal = "Postnatal", interaction = "Interaction"),
    name = NULL
  ) +
  scale_size_continuous(
    range = c(1.1, 3.4),
    limits = c(0.7, 1.0),
    breaks = c(0.7, 0.8, 0.9, 1.0),
    name = expression(PP[H4])
  ) +
  scale_x_continuous(
    breaks = seq_along(gwas_order),
    labels = gwas_order,
    expand = c(0.01, 0.01),
    limits = c(0.5, length(gwas_order) + 0.5)
  ) +
  scale_y_discrete(
    limits = gene_key_order,
    labels = setNames(sub("\\|.*$", "", gene_key_order), gene_key_order)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = base_size, base_family = "Helvetica") +
  theme(
    text = element_text(color = ink),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.2, color = "black"),
    axis.text.y = element_text(size = 5.6, color = "black", face = "italic"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(0.25, "cm"),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2),
    legend.key.size = unit(0.30, "cm"),
    plot.margin = margin(4, 6, 2, 4)
  ) +
  guides(
    color = guide_legend(order = 1, title = NULL,
                         override.aes = list(size = 2.0)),
    size  = guide_legend(order = 2)
  )

# --- save ---
png_path <- file.path(fig_dir, "07_summary_matrix.png")
pdf_path <- file.path(fig_dir, "07_summary_matrix.pdf")
ggsave(png_path, p, width = figure_width, height = figure_height, dpi = dpi, bg = "white")
ggsave(pdf_path, p, width = figure_width, height = figure_height, device = cairo_pdf, bg = "white")
