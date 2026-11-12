#!/usr/bin/env Rscript
# Figure for EUR APA-QTL enrichment from precomputed GO tables.
if (!nzchar(Sys.getenv("APA_ROOT"))) stop("Set APA_ROOT (source config/paths.sh)")

# --- load packages ---

suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# --- config ---

apa_root <- Sys.getenv("APA_ROOT")
magma_dir <- file.path(apa_root, "data", "06_disease")
figure_dir <- file.path(apa_root, "figures", "06_disease")

# --- load results ---

go_shared <- fread(file.path(magma_dir, "enrichment_GO_shared.tsv"))
go_postnatal <- fread(file.path(magma_dir, "enrichment_GO_postnatal_only.tsv"))

# --- theme ---

theme_enrich <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text.y = element_text(color = "black", size = base_size),
      axis.text.x = element_text(color = "black", size = base_size - 1),
      axis.title = element_text(color = "black", size = base_size + 1),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 2),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size)
    )
}

# --- panel A: gene counts ---

count_data <- data.table(
  group = factor(c("Prenatal", "Postnatal", "Shared", "Prenatal-only", "Postnatal-only"),
                 levels = c("Postnatal-only", "Prenatal-only", "Shared", "Postnatal", "Prenatal")),
  count = c(133, 443, 113, 20, 330),
  category = c("Total", "Total", "Overlap", "Specific", "Specific")
)

p_counts <- ggplot(count_data, aes(x = count, y = group, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = count), hjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c("Total" = "grey65", "Overlap" = "#009E73", "Specific" = "#0072B2"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene count", y = NULL, title = "A  EUR APA-QTL genes (FDR < 0.05)") +
  theme_enrich() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0, size = 12)
  )

# --- panel B: shared GO ---

plot_go <- function(df, title, n = 8) {
  df <- head(df[order(p.adjust)], n)
  df[, GeneRatio_num := as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))]
  df[, Term := substr(Description, 1, 55)]
  
  ggplot(df, aes(x = GeneRatio_num, y = reorder(Term, -log10(p.adjust)))) +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_viridis_c(option = "plasma", direction = -1, name = "Adj. P") +
    scale_size_continuous(range = c(3, 7), name = "Genes") +
    labs(x = "Gene ratio", y = NULL, title = title) +
    theme_enrich()
}

p_shared_go <- plot_go(go_shared, "B  Shared (Prenatal ∩ Postnatal)")

# --- panel C: postnatal-only GO ---

p_post_go <- plot_go(go_postnatal, "C  Postnatal-only")

# --- assemble figure ---

figure <- p_counts / p_shared_go / p_post_go +
  plot_layout(heights = c(1, 1.2, 0.6)) +
  plot_annotation(
    subtitle = "EUR GWAS, genome-wide FDR < 0.05",
    tag_levels = NULL
  )

ggsave(file.path(figure_dir, "enrichment_figure.pdf"), figure,
       width = 8, height = 9, device = cairo_pdf)
ggsave(file.path(figure_dir, "enrichment_figure.png"), figure,
       width = 8, height = 9, dpi = 300)
