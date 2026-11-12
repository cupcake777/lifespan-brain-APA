#!/usr/bin/env Rscript
# Tier-1 (SuSiE-confirmed) colocalization summary and curated locuscompare panels.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(stringr)
  library(forcats)
})

# --- config ---
stage_levels   <- c("prenatal", "postnatal", "interaction")
stage_labels   <- c(prenatal = "Prenatal", postnatal = "Postnatal", interaction = "Interaction")
stage_colors   <- c(prenatal = "#E64B35", postnatal = "#4DBBD5", interaction = "#00A087")

data_dir <- file.path(Sys.getenv("APA_ROOT"), "data/07_coloc")
out_dir <- file.path(Sys.getenv("APA_ROOT"), "figures/07_coloc")
coloc_file     <- file.path(data_dir, "coloc_run", "all_colocalization.csv")
susie_file     <- file.path(data_dir, "all_susie_results.csv")
json_dir       <- data_dir

susie_confirm_pp4 <- 0.7   # same threshold as ABF significance
figure_dpi        <- 320

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# --- load data ---
abf <- read_csv(coloc_file, show_col_types = FALSE) %>%
  mutate(
    GWAS_ID = as.character(GWAS_ID),
    QTL_ID  = as.character(QTL_ID),
    Locus   = as.character(Locus),
    PP4     = as.numeric(PP4),
    n_snps  = as.numeric(n_snps),
    key     = sub("_coloc\\.png$", "", sub("^plots/", "", basename(plot_file)))
  ) %>%
  filter(!is.na(key), !duplicated(key), PP4 >= 0.7)

pheno_parts <- str_split_fixed(abf$Phenotype, "\\|", 4)
abf <- abf %>%
  mutate(
    ancestry = str_extract(GWAS_ID, "^[^_]+"),
    disease  = str_remove(GWAS_ID, "^[^_]+_"),
    gene     = if_else(is.na(pheno_parts[, 2]) | pheno_parts[, 2] == "", Phenotype, pheno_parts[, 2]),
    stage    = factor(QTL_ID, levels = stage_levels)
  ) %>%
  filter(!is.na(stage), !is.na(PP4))

susie <- read_csv(susie_file, show_col_types = FALSE) %>%
  mutate(
    key = sub("_susie\\.csv$", "", basename(source_file)),
    PP.H4.abf = as.numeric(PP.H4.abf)
  ) %>%
  filter(!is.na(key), !is.na(PP.H4.abf))

# Per-locus SuSiE summary: best (max PP.H4.abf) CS pair + number of pairs
susie_locus <- susie %>%
  group_by(key) %>%
  summarise(
    susie_max_PP4 = max(PP.H4.abf),
    susie_n_pairs = n(),
    susie_nsnps   = nsnps[which.max(PP.H4.abf)],
    hit1          = hit1[which.max(PP.H4.abf)],
    hit2          = hit2[which.max(PP.H4.abf)],
    idx1          = idx1[which.max(PP.H4.abf)],
    idx2          = idx2[which.max(PP.H4.abf)],
    .groups = "drop"
  )

# --- tier assignment ---
tiers <- abf %>%
  left_join(susie_locus, by = "key") %>%
  mutate(
    tier = case_when(
      !is.na(susie_max_PP4) & susie_max_PP4 >= susie_confirm_pp4 ~ "Tier1",
      !is.na(susie_max_PP4) & susie_max_PP4 <  susie_confirm_pp4 ~ "Tier2a",
      TRUE ~ "Tier2b"
    ),
    cs_agree = if_else(!is.na(hit1) & !is.na(hit2) & hit1 == hit2, TRUE, FALSE)
  )

tier_counts <- tiers %>% count(tier)
message(sprintf("Tier assignment: %s", paste(sprintf("%s=%d", tier_counts$tier, tier_counts$n), collapse = ", ")))

write_tsv(
  tiers %>%
    select(key, GWAS_ID, QTL_ID, stage, gene, Locus, PP4, n_snps,
           susie_max_PP4, susie_n_pairs, susie_nsnps, hit1, hit2, idx1, idx2,
           cs_agree, tier),
  file.path(data_dir, "04_tier_assignment.tsv")
)

# --- summary figure: tier-1 landscape (no funnel) ---
tier1 <- tiers %>% filter(tier == "Tier1")

# Y order: gene by max SuSiE PP4 (desc); label multi-locus genes with rsid
gene_order <- tier1 %>%
  group_by(gene) %>%
  summarise(max_PP4 = max(susie_max_PP4), .groups = "drop") %>%
  arrange(desc(max_PP4)) %>%
  pull(gene)

gene_n <- tier1 %>% count(gene, name = "n_loci")
tier1 <- tier1 %>%
  left_join(gene_n, by = "gene") %>%
  mutate(
    ylabel = if_else(n_loci > 1, paste0(gene, " (", Locus, ")"), gene),
    ylabel = fct_reorder(ylabel, susie_max_PP4),
    ylabel = fct_rev(ylabel)
  )

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black", margin = margin(r = 3)),
    axis.title = element_text(size = 10, face = "bold", colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.title = element_text(face = "bold", size = 12, hjust = 0, colour = "black"),
    plot.margin = margin(10, 16, 10, 16)
  )

n_t1 <- nrow(tier1)
n_agree <- sum(tier1$cs_agree)
p_summary <- ggplot(tier1, aes(x = stage, y = ylabel)) +
  geom_point(aes(fill = stage, size = susie_max_PP4,
                 shape = cs_agree),
             position = position_dodge(width = 0.35),
             colour = "grey30", stroke = 0.35, alpha = 0.9) +
  geom_text_repel(
    data = tier1 %>% filter(!cs_agree),
    aes(label = Locus), size = 2.4, colour = "black",
    min.segment.length = 0, segment.size = 0.2, segment.colour = "grey50",
    max.overlaps = Inf, nudge_y = 0.35, box.padding = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = stage_colors, name = "Stage", drop = FALSE,
                    labels = stage_labels) +
  scale_shape_manual(
    values = c("TRUE" = 21, "FALSE" = 24),
    name = "CS pair",
    labels = c("TRUE" = "Shared CS pair", "FALSE" = "Different CS pairs")
  ) +
  scale_size_continuous(range = c(2, 5.5), name = "SuSiE PP4",
                        breaks = c(0.7, 0.8, 0.9, 1.0), limits = c(0.7, 1.0)) +
  scale_x_discrete(labels = stage_labels, expand = c(0.15, 0.15)) +
  labs(
    title = sprintf("SuSiE-confirmed colocalizations (n = %d)", n_t1),
    x = NULL, y = NULL
  ) +
  theme_pub

fig_height <- max(5, 0.20 * n_distinct(tier1$ylabel) + 2.0)
pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave(file.path(out_dir, "04_summary_tier1.pdf"), p_summary,
       width = 8, height = fig_height, units = "in", device = pdf_device)
ggsave(file.path(out_dir, "04_summary_tier1.png"), p_summary,
       width = 8, height = fig_height, units = "in", dpi = figure_dpi)

# --- curated locuscompare examples ---
# Interesting loci: cleanest single-CS hits, prenatal examples, multi-signal
# (hit1 != hit2) cases, and one Tier-2a contrast (ABF significant, SuSiE not).
curated <- c(
  "EAS_SCZ_interaction_chr12_89410263_rs12578201_ATP2B1",     # PP4=1, cleanest
  "EUR_ADHD_prenatal_chr3_49879277_rs7613360_GNAI2",          # prenatal, PP4=1
  "EUR_SCZ_postnatal_chr16_29983601_rs3814883_TLCD3B",        # multi-signal (hit1!=hit2)
  "EUR_SCZ_postnatal_chr16_29983601_rs3814883_PPP4C",         # multi-signal (hit1!=hit2)
  "EUR_MDD_postnatal_chr12_120364215_rs3782901_SPPL3",        # MDD, PP4=1
  "EAS_SCZ_postnatal_chr10_102760520_rs12571643_ARL3",        # EAS, PP4=0.998
  "EAS_AIS_postnatal_chr12_110886135_rs7311323_TCHP"          # Tier2a contrast
)

ld_bin_colors <- c(
  "0.8-1.0" = "#D73027",
  "0.6-0.8" = "#FC8D59",
  "0.4-0.6" = "#FEE090",
  "0.2-0.4" = "#91BFDB",
  "0.0-0.2" = "#4575B4",
  "Missing" = "#C8C8C8"
)
ld_bin_levels <- names(ld_bin_colors)

theme_lc <- theme_classic(base_size = 10) +
  theme(
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.3, colour = "black"),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_text(size = 9, face = "bold", colour = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 9, hjust = 0, colour = "black"),
    plot.subtitle = element_text(size = 7.5, colour = "grey25"),
    legend.position = "right",
    plot.margin = margin(5, 8, 5, 5)
  )

# --- curated locuscompare examples ---
for (key in curated) {
  json_path <- file.path(json_dir, paste0(key, "_locuscompare.json"))
  if (!file.exists(json_path)) {
    message(sprintf("  -> JSON missing for %s, skipping", key))
    next
  }
  d <- jsonlite::fromJSON(json_path)
  pts <- as_tibble(d$points) %>%
    mutate(
      ld_bin = case_when(
        is.na(r2) ~ "Missing",
        r2 >= 0.8 ~ "0.8-1.0",
        r2 >= 0.6 ~ "0.6-0.8",
        r2 >= 0.4 ~ "0.4-0.6",
        r2 >= 0.2 ~ "0.2-0.4",
        TRUE      ~ "0.0-0.2"
      ),
      ld_bin = factor(ld_bin, levels = ld_bin_levels)
    )
  lead <- pts %>% filter(is_lead)
  if (nrow(lead) == 0) lead <- pts %>% arrange(desc(gwas_logp + qtl_logp)) %>% slice(1)

  meta <- d$meta
  stage_i <- match(meta$qtl_id, stage_levels)
  stage_col <- if (!is.na(stage_i)) stage_colors[stage_i] else "#333333"

  pp4_text <- sprintf("PP4 = %.3f", as.numeric(meta$pp4))
  if (!is.null(meta$credible_set_size) && meta$credible_set_size > 0) {
    pp4_text <- sprintf("%s | CS = %d SNPs", pp4_text, meta$credible_set_size)
  }

  p <- ggplot(pts, aes(x = gwas_logp, y = qtl_logp)) +
    geom_point(aes(fill = ld_bin), shape = 21, colour = "grey50",
               stroke = 0.2, size = 1.9, alpha = 0.85) +
    # credible set variants: bold outline
    geom_point(data = pts %>% filter(is_credible),
               shape = 21, colour = "black", stroke = 0.7,
               size = 2.3, show.legend = FALSE) +
    geom_point(data = lead, shape = 23, fill = stage_col, colour = "black",
               size = 3.0, stroke = 0.4, show.legend = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey60", linewidth = 0.3) +
    geom_text_repel(
      data = lead, aes(label = snp),
      size = 2.6, fontface = "bold", colour = "black",
      min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40",
      nudge_x = 0.3, nudge_y = 0.3, box.padding = 0.4, max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = ld_bin_colors, breaks = ld_bin_levels, limits = ld_bin_levels,
      name = expression(LD~r^2), drop = FALSE,
      guide = guide_legend(override.aes = list(shape = 21, size = 2.5, stroke = 0.15))
    ) +
    labs(
      title = sprintf("%s | %s | %s", meta$gwas_id,
                      stage_labels[meta$qtl_id],
                      str_split_fixed(meta$phenotype, "\\|", 4)[, 2]),
      subtitle = pp4_text,
      x = expression(-log[10](GWAS~italic(P))),
      y = expression(-log[10](APA-QTL~italic(P)))
    ) +
    theme_lc

  safe_name <- str_replace_all(key, "[^A-Za-z0-9._-]+", "_")
  ggsave(file.path(out_dir, paste0("04_locuscompare_", safe_name, ".pdf")), p,
         width = 6.5, height = 5.2, units = "in", device = pdf_device)
  ggsave(file.path(out_dir, paste0("04_locuscompare_", safe_name, ".png")), p,
         width = 6.5, height = 5.2, units = "in", dpi = figure_dpi)
}
