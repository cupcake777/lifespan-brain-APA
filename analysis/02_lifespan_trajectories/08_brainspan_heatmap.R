# Exploratory BrainSpan heatmap of raw PDUI for genes in one trajectory class.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# --- paths ---

root <- Sys.getenv("APA_ROOT")
DATA <- file.path(root, "data/02_traj")
FIG <- file.path(root, "figures/02_traj")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
MAPPING_FILE <- file.path(DATA, "08_external_mapping_qc.tsv")
SITES_FILE <- file.path(DATA, "02_sites_clustered.tsv")
BRAINSPAN_RAW_FILE <- file.path(root, "raw/all.apa.brainspan.count")
BRAINSPAN_META_FILE <- file.path(root, "raw/all.apa.brainspan.meta")

# --- settings ---

MAX_GENE_MISSING <- 0.10
HEATMAP_BINS_PER_CLUSTER <- 60L

# --- matched events and raw BrainSpan PDUI ---

mapping <- fread(MAPPING_FILE, select = c(
  "cohort", "site_id", "external_site_id", "gene_id", "cluster_id", "pas_distance_bp", "mapped_qc"
))
sites <- fread(SITES_FILE, select = c("site_id", "delta_logit"))
mapping <- merge(mapping[cohort == "BrainSpan" & mapped_qc == TRUE], sites, by = "site_id")
mapping <- mapping[!is.na(gene_id) & gene_id != "" & gene_id != "NA"]

gene_membership <- unique(mapping[, .(gene_id, cluster_id)])[, .(
  n_clusters = uniqueN(cluster_id)
), by = gene_id]
pure_genes <- gene_membership[n_clusters == 1L, gene_id]
mapping <- mapping[gene_id %chin% pure_genes]

gene_discovery <- unique(mapping[, .(gene_id, cluster_id, site_id, delta_logit)])[, .(
  discovery_delta = median(delta_logit),
  n_discovery_sites = uniqueN(site_id)
), by = .(gene_id, cluster_id)]

mapping <- unique(mapping[, .(external_site_id, gene_id, cluster_id)],
                  by = c("external_site_id", "gene_id", "cluster_id"))

brainspan_raw <- fread(BRAINSPAN_RAW_FILE, showProgress = FALSE)
gene_fields <- tstrsplit(brainspan_raw$Gene, "|", fixed = TRUE)
brainspan_raw[, external_site_id := paste(
  gene_fields[[1]], gene_fields[[2]], Loci, gene_fields[[4]], sep = "|"
)]
brainspan_raw <- brainspan_raw[external_site_id %chin% mapping$external_site_id]

brainspan_meta <- fread(BRAINSPAN_META_FILE)
setorder(brainspan_meta, age, sample)
samples <- intersect(brainspan_meta$sample, names(brainspan_raw))
brainspan_meta <- brainspan_meta[match(samples, sample)]

# --- gene-level raw PDUI ---

raw_long <- melt(brainspan_raw, id.vars = "external_site_id", measure.vars = samples,
                 variable.name = "sample", value.name = "pdui")
raw_long <- merge(raw_long, mapping, by = "external_site_id")

gene_event_counts <- mapping[, .(n_brainspan_events = uniqueN(external_site_id)),
                             by = .(gene_id, cluster_id)]

gene_sample <- raw_long[, .(
  pdui = if (all(is.na(pdui))) NA_real_ else median(pdui, na.rm = TRUE)
), by = .(gene_id, cluster_id, sample)]
gene_sample <- merge(gene_sample, brainspan_meta[, .(sample, age, LifeStage)], by = "sample")

# --- gene-level direction summary ---

gene_stats <- gene_sample[, {
  finite <- is.finite(pdui) & is.finite(age)
  prenatal <- pdui[LifeStage == "Prenatal" & is.finite(pdui)]
  postnatal <- pdui[LifeStage == "Postnatal" & is.finite(pdui)]
  list(
    missing_fraction = mean(!is.finite(pdui)),
    rho_age = if (sum(finite) >= 5L && sd(pdui[finite]) > 0) {
      suppressWarnings(cor(pdui[finite], age[finite], method = "spearman"))
    } else {
      NA_real_
    },
    stage_delta = if (length(prenatal) > 0L && length(postnatal) > 0L) {
      median(postnatal) - median(prenatal)
    } else {
      NA_real_
    }
  )
}, by = .(gene_id, cluster_id)]

gene_stats <- merge(gene_stats, gene_event_counts, by = c("gene_id", "cluster_id"))
gene_stats <- merge(gene_stats, gene_discovery, by = c("gene_id", "cluster_id"))
gene_stats[, `:=`(
  rho_direction_concordant = fifelse(cluster_id == 1L, rho_age > 0, rho_age < 0),
  stage_direction_concordant = fifelse(cluster_id == 1L, stage_delta > 0, stage_delta < 0),
  heatmap_pass = missing_fraction <= MAX_GENE_MISSING & is.finite(rho_age)
)]
setorder(gene_stats, cluster_id, gene_id)
fwrite(gene_stats, file.path(DATA, "11_brainspan_gene_heatmap_stats.tsv"), sep = "\t")

cluster_summary <- gene_stats[heatmap_pass == TRUE, .(
  n_genes = .N,
  n_brainspan_events = sum(n_brainspan_events),
  median_rho_age = median(rho_age),
  fraction_rho_direction_concordant = mean(rho_direction_concordant),
  median_postnatal_minus_prenatal_pdui = median(stage_delta),
  fraction_stage_direction_concordant = mean(stage_direction_concordant),
  median_missing_fraction = median(missing_fraction),
  n_samples = length(samples)
), by = cluster_id]
fwrite(cluster_summary, file.path(DATA, "11_brainspan_gene_heatmap_summary.tsv"), sep = "\t")

# --- heatmap rows and age ordering ---

heatmap_genes <- gene_stats[heatmap_pass == TRUE]
heatmap_genes[, order_score := fifelse(cluster_id == 1L, -discovery_delta, discovery_delta)]
setorder(heatmap_genes, cluster_id, order_score, gene_id)
heatmap_genes[, gene_rank := seq_len(.N), by = cluster_id]
heatmap_genes[, heat_bin := pmin(
  HEATMAP_BINS_PER_CLUSTER,
  pmax(1L, ceiling(gene_rank / .N * HEATMAP_BINS_PER_CLUSTER))
), by = cluster_id]

heat_long <- merge(gene_sample, heatmap_genes[, .(gene_id, cluster_id, heat_bin)],
                   by = c("gene_id", "cluster_id"))
heat_long[, z := {
  finite <- is.finite(pdui)
  center <- mean(pdui[finite])
  spread <- sd(pdui[finite])
  fifelse(finite, (pdui - center) / spread, NA_real_)
}, by = .(gene_id, cluster_id)]

sample_axis <- brainspan_meta[, .(sample, sample_index = seq_len(.N), age, LifeStage)]
heat_long <- merge(heat_long, sample_axis, by = c("sample", "age", "LifeStage"))

heat_plot <- heat_long[, .(
  z = median(z, na.rm = TRUE),
  n_genes = uniqueN(gene_id)
), by = .(cluster_id, heat_bin, sample_index, age, LifeStage)]
heat_plot[!is.finite(z), z := NA_real_]
heat_plot[, cluster_label := factor(paste0("C", cluster_id), levels = c("C1", "C2"))]

# --- exploratory BrainSpan heatmap ---

x_breaks <- unique(round(seq(1, nrow(sample_axis), length.out = 8L)))
x_labels <- format(round(sample_axis$age[x_breaks], 2), trim = TRUE)
birth_boundary <- sum(sample_axis$age < 0) + 0.5

p_heatmap <- ggplot(heat_plot, aes(sample_index, heat_bin, fill = z)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = birth_boundary, color = "#4B5563",
             linewidth = 0.35, linetype = "dashed") +
  facet_grid(cluster_label ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B"),
    values = rescale(c(-2.5, -1, 0, 1, 2.5)),
    limits = c(-2.5, 2.5),
    oob = squish,
    na.value = "grey75",
    name = "Raw PDUI\nrow Z-score"
  ) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Age (years; samples ordered youngest to oldest)",
    y = NULL,
    title = "BrainSpan matched-gene APA patterns"
  ) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold"),
    axis.title.x = element_text(color = "black", face = "bold"),
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", color = "black"),
    strip.placement = "outside",
    panel.spacing = unit(0.18, "lines"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  )

ggsave(file.path(FIG, "11_brainspan_gene_heatmap.png"), p_heatmap, width = 6.4, height = 5.3,
       dpi = 300, bg = "white")
ggsave(file.path(FIG, "11_brainspan_gene_heatmap.pdf"), p_heatmap, width = 6.4, height = 5.3,
       device = cairo_pdf, bg = "white")

# --- layout validation ---

layout_validation <- heat_plot[, .(
  figure = "11_brainspan_gene_heatmap",
  artifact_class = "exploratory plot",
  n_genes = uniqueN(heat_long[cluster_id == .BY$cluster_id]$gene_id),
  n_heatmap_bins = uniqueN(heat_bin),
  n_samples = uniqueN(sample_index),
  age_min = min(age),
  age_max = max(age),
  missing_tile_fraction = mean(!is.finite(z)),
  colorbar = "bottom horizontal",
  validation_note = paste(
    "Raw PDUI is aggregated within gene and cluster, row-standardized,",
    "and compressed into discovery-effect-ordered gene bins; grey marks missing values."
  )
), by = cluster_id]
fwrite(layout_validation, file.path(FIG, "11_brainspan_gene_heatmap_layout_validation.tsv"), sep = "\t")
