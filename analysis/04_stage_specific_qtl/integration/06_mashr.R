# mashr on apaQTL and eQTL pair effects per life stage: posterior classes, effect consistency,
# gene-level APA-eQTL coupling, and downsampling replication. Reads data/04_qtl/mash/;
# writes 05_* tables to data/04_qtl and figures to figures/04_qtl.
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(mashr)
})

# --- paths and settings ---
APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
mash_dir <- file.path(D, "mash")
figure_dir <- file.path(APA_ROOT, "figures", "04_qtl")
gtf_file <- file.path(APA_ROOT, "ref", "gencode.v40.annotation.gtf.gz")

qtl_labels <- c(apa = "apaQTL", eqtl = "eQTL")
random_target <- 100000
lfsr_cut <- 0.05
share_factor <- 0.5
seed <- 1
fig_w15 <- 5.35
fig_w2 <- 7.09
fig_hmax <- 9.45
png_dpi <- 600
base_font <- 8
font_sm <- 7
font_tag <- 9
class_levels <- c("Prenatal-specific", "Postnatal-specific", "Shared", "Not significant")
class_cols <- c("Prenatal-specific" = "#B2182B", "Postnatal-specific" = "#2166AC",
                "Shared" = "#4D4D4D", "Not significant" = "grey85")
effect_class_short <- c("Prenatal-specific" = "Prenatal-only", "Postnatal-specific" = "Postnatal-only",
                        "Shared" = "Shared", "Not significant" = "n.s.")
gene_state_levels <- c("Shared", "Prenatal-only", "Postnatal-only")
set.seed(seed)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- helpers ---
pt_mm <- function(pt) pt / .pt

theme_mash <- function() {
  theme_classic(base_size = base_font) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = font_sm),
      axis.title = element_text(colour = "black", size = base_font),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(colour = "black", size = base_font, face = "bold", margin = margin(b = 1.5)),
      legend.title = element_text(size = font_sm),
      legend.text = element_text(size = font_sm),
      legend.key = element_blank(),
      legend.key.size = grid::unit(0.30, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      plot.title = element_text(size = base_font, face = "bold"),
      plot.margin = margin(2, 3, 2, 2)
    )
}

save_plot <- function(plot, stem, width, height) {
  if (height > fig_hmax) warning(sprintf("%s: height %.2f in exceeds %.2f in", stem, height, fig_hmax))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  unlink(c(png_file, pdf_file))
  ggsave(png_file, plot, width = width, height = height, dpi = png_dpi, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
}

read_stage_pairs <- function(qtl, stage) {
  dt <- fread(file.path(mash_dir, sprintf("%s.%s.pairs.tsv", qtl, stage)), sep = "\t", showProgress = FALSE)
  setnames(dt, c("slope", "slope_se"), c(paste0("b_", stage), paste0("s_", stage)))
  dt[, pval_nominal := NULL]
  dt
}

as_mash_matrices <- function(dt) {
  list(Bhat = as.matrix(dt[, .(prenatal = b_prenatal, postnatal = b_postnatal)]),
       Shat = as.matrix(dt[, .(prenatal = s_prenatal, postnatal = s_postnatal)]))
}

gene_state_from_pairs <- function(dt, qtl_name, gene_map) {
  qtl_label <- qtl_labels[[qtl_name]]
  x <- copy(dt[qtl == qtl_label])
  if (qtl_name == "apa") {
    x[, gene_name := vapply(strsplit(phenotype_id, "\\|"), function(z) {
      if (length(z) >= 2L && nzchar(z[2])) z[2] else z[1]
    }, character(1))]
  } else {
    x[, ensg0 := sub("\\..*", "", phenotype_id)]
    x <- merge(x, gene_map, by = "ensg0", all.x = TRUE)
    x[is.na(gene_name) | gene_name == "", gene_name := ensg0]
  }
  gene_dt <- x[, .(
    n_pairs = .N,
    n_detected = sum(lfsr_prenatal < lfsr_cut | lfsr_postnatal < lfsr_cut, na.rm = TRUE),
    prenatal_hit = any(lfsr_prenatal < lfsr_cut, na.rm = TRUE),
    postnatal_hit = any(lfsr_postnatal < lfsr_cut, na.rm = TRUE)
  ), by = gene_name]
  gene_dt[, gene_state := fcase(prenatal_hit & postnatal_hit, "Shared", prenatal_hit, "Prenatal-only",
                                postnatal_hit, "Postnatal-only", default = "Not detected")]
  gene_dt[, qtl := qtl_label]
  gene_dt[]
}

read_perm_summary <- function(path, seed = NA_integer_) {
  dt <- fread(path, select = c("phenotype_id", "variant_id", "group_id", "group_size", "qval"), showProgress = FALSE)
  dt[, seed := seed]
  dt
}

# --- mashr input: prenatal and postnatal effects per pair (both stages present, finite SEs) ---
pairs_by_qtl <- lapply(names(qtl_labels), function(qtl) {
  pre <- read_stage_pairs(qtl, "prenatal")
  post <- read_stage_pairs(qtl, "postnatal")
  setnames(pre, "set", "set_pre")
  setnames(post, "set", "set_post")
  both <- merge(pre, post, by = c("phenotype_id", "variant_id"))
  both[, set := fifelse(set_pre == "strong" | set_post == "strong", "strong", "random")]
  both[, c("set_pre", "set_post") := NULL]
  both[is.finite(b_prenatal) & is.finite(b_postnatal) & is.finite(s_prenatal) & is.finite(s_postnatal) &
         s_prenatal > 0 & s_postnatal > 0]
})
names(pairs_by_qtl) <- names(qtl_labels)

# --- mash fit: data-driven covariance from random pairs, posterior for strong pairs ---
results <- lapply(names(qtl_labels), function(qtl_name) {
  pairs <- pairs_by_qtl[[qtl_name]]
  strong <- pairs[set == "strong"]
  random <- pairs[set == "random"]
  if (nrow(random) > random_target) random <- random[sample.int(nrow(random), random_target)]
  m_random <- as_mash_matrices(random)
  m_strong <- as_mash_matrices(strong)
  d_random0 <- mash_set_data(m_random$Bhat, m_random$Shat)
  vhat <- estimate_null_correlation_simple(d_random0)
  d_random <- mash_set_data(m_random$Bhat, m_random$Shat, V = vhat)
  d_strong <- mash_set_data(m_strong$Bhat, m_strong$Shat, V = vhat)
  u_pca <- tryCatch(cov_pca(d_strong, npc = 2), error = function(e) NULL)
  ulist <- c(u_pca, cov_canonical(d_random))
  not_matrix <- !vapply(ulist, is.matrix, logical(1))
  if (any(not_matrix)) {
    stop(sprintf("%s: non-matrix covariance component: %s", qtl_name, paste(names(ulist)[not_matrix], collapse = ", ")))
  }
  fit <- mash(d_random, Ulist = ulist, outputlevel = 1)
  post <- mash(d_strong, g = get_fitted_g(fit), fixg = TRUE)
  list(qtl = qtl_name, strong = strong, fit = fit, post = post, Vhat = vhat)
})
names(results) <- names(qtl_labels)

# --- posterior classes per strong pair ---
posterior_dt <- rbindlist(lapply(results, function(res) {
  lfsr <- get_lfsr(res$post)
  pm <- get_pm(res$post)
  psd <- get_psd(res$post)
  dt <- copy(res$strong)
  dt[, `:=`(lfsr_prenatal = lfsr[, "prenatal"], lfsr_postnatal = lfsr[, "postnatal"],
            pm_prenatal = pm[, "prenatal"], pm_postnatal = pm[, "postnatal"],
            psd_prenatal = psd[, "prenatal"], psd_postnatal = psd[, "postnatal"])]
  dt[, ratio := pmin(abs(pm_prenatal), abs(pm_postnatal)) / pmax(abs(pm_prenatal), abs(pm_postnatal))]
  dt[, shared_by_magnitude := sign(pm_prenatal) == sign(pm_postnatal) & is.finite(ratio) & ratio >= share_factor]
  sig_pre <- dt$lfsr_prenatal < lfsr_cut
  sig_post <- dt$lfsr_postnatal < lfsr_cut
  dt[, mash_class := fcase(
    sig_pre & sig_post, "Shared",
    sig_pre & !sig_post & !shared_by_magnitude, "Prenatal-specific",
    !sig_pre & sig_post & !shared_by_magnitude, "Postnatal-specific",
    sig_pre & !sig_post & shared_by_magnitude, "Shared",
    !sig_pre & sig_post & shared_by_magnitude, "Shared",
    default = "Not significant")]
  dt[, mash_class := factor(mash_class, levels = class_levels)]
  dt[, qtl := qtl_labels[[res$qtl]]]
  dt[]
}), use.names = TRUE)
posterior_dt[, qtl := factor(qtl, levels = unname(qtl_labels))]
for (qtl_name in names(qtl_labels)) {
  label <- qtl_labels[[qtl_name]]
  fwrite(posterior_dt[qtl == label], file.path(D, sprintf("05_mash_posterior_%s.tsv", qtl_name)), sep = "\t", quote = FALSE)
}

# --- mash summary: class counts, pairwise sharing, null correlation, mixture weights ---
pi_dt <- rbindlist(lapply(names(qtl_labels), function(qtl_name) {
  pi <- get_estimated_pi(results[[qtl_name]]$fit)
  data.table(qtl = qtl_labels[[qtl_name]], component = names(pi), weight = as.numeric(pi))
}))
sharing_dt <- rbindlist(lapply(names(qtl_labels), function(qtl_name) {
  sharing <- get_pairwise_sharing(results[[qtl_name]]$post, factor = share_factor)
  data.table(qtl = qtl_labels[[qtl_name]], sharing = sharing["prenatal", "postnatal"],
             null_cor = results[[qtl_name]]$Vhat["prenatal", "postnatal"])
}))
class_dt <- posterior_dt[, .N, by = .(qtl, mash_class)]
class_dt[, frac := N / sum(N), by = qtl]
summary_dt <- rbindlist(list(
  class_dt[, .(qtl, metric = paste0("class_", mash_class), value = as.numeric(N), fraction = frac)],
  sharing_dt[, .(qtl, metric = "pairwise_sharing_by_magnitude", value = sharing, fraction = NA_real_)],
  sharing_dt[, .(qtl, metric = "null_correlation", value = null_cor, fraction = NA_real_)],
  pi_dt[, .(qtl, metric = paste0("pi_", component), value = weight, fraction = NA_real_)]
), use.names = TRUE)
fwrite(summary_dt, file.path(D, "05_mash_summary.tsv"), sep = "\t", quote = FALSE)

# --- effect consistency among detected pairs ---
class_detected <- posterior_dt[mash_class != "Not significant", .N, by = .(qtl, mash_class)]
class_detected[, frac := N / sum(N), by = qtl]
class_detected[, class_label := fcase(mash_class == "Shared", "Shared",
                                      mash_class == "Postnatal-specific", "Postnatal-only",
                                      mash_class == "Prenatal-specific", "Prenatal-only")]
setorder(class_detected, qtl, -frac)
direction_dt <- posterior_dt[mash_class != "Not significant", .(
  n_detected = .N,
  same_direction = mean(sign(pm_prenatal) == sign(pm_postnatal)),
  pairwise_sharing = sharing_dt[qtl == .BY$qtl, sharing]
), by = qtl]
direction_dt[, badge := sprintf("Same sign %.1f%%\nWithin 2-fold %.1f%%", 100 * same_direction, 100 * pairwise_sharing)]
effect_note <- class_detected[, .(label = paste(sprintf("%s %.1f%%", class_label, 100 * frac), collapse = "\n")), by = qtl]
lim <- posterior_dt[, max(abs(c(pm_prenatal, pm_postnatal)), na.rm = TRUE)] * 1.05
effect_note[, `:=`(x = -0.965 * lim, y = 0.965 * lim)]
badge_dt <- copy(direction_dt)
badge_dt[, `:=`(x = 0.92 * lim, y = -0.92 * lim)]
fwrite(class_detected, file.path(D, "05_effect_consistency_class_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(direction_dt, file.path(D, "05_effect_consistency_direction_summary.tsv"), sep = "\t", quote = FALSE)

effect_fig <- ggplot(posterior_dt[order(mash_class == "Not significant", decreasing = TRUE)],
                     aes(x = pm_prenatal, y = pm_postnatal, colour = mash_class)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.35, colour = "grey70") +
  geom_abline(slope = c(0.5, 2), intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey70") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey88") +
  geom_point(size = 0.45, alpha = 0.72, stroke = 0) +
  geom_text(data = effect_note, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1,
            size = pt_mm(font_sm), lineheight = 0.9) +
  geom_label(data = badge_dt, aes(x = x, y = y, label = badge), inherit.aes = FALSE, hjust = 1, vjust = 0,
             size = pt_mm(font_sm), linewidth = 0.2, label.padding = grid::unit(1.5, "pt"),
             fill = "white", colour = "black", lineheight = 0.9) +
  scale_colour_manual(values = class_cols, labels = effect_class_short, name = NULL, drop = FALSE) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  facet_wrap(~ qtl, nrow = 1) +
  labs(x = "Posterior mean, prenatal", y = "Posterior mean, postnatal") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.8, alpha = 1))) +
  theme_mash() + theme(legend.position = "bottom")
save_plot(effect_fig, "mash_effect_scatter", width = fig_w15, height = 2.8)

# --- gene-level APA-eQTL coupling ---
gene_map <- {
  attrs <- fread(cmd = sprintf("zcat %s | awk -F'\\t' '$3==\"gene\" {print $9}'", shQuote(gtf_file)),
                 header = FALSE, sep = "\t", quote = "", showProgress = FALSE)
  attrs[, ensg0 := sub("\\..*", "", sub('.*gene_id "([^"]+)".*', "\\1", V1))]
  attrs[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", V1)]
  unique(attrs[, .(ensg0, gene_name)])
}
apa_gene_state <- gene_state_from_pairs(posterior_dt, "apa", gene_map)
eqtl_gene_state <- gene_state_from_pairs(posterior_dt, "eqtl", gene_map)
fwrite(apa_gene_state, file.path(D, "05_gene_state_apa.tsv"), sep = "\t", quote = FALSE)
fwrite(eqtl_gene_state, file.path(D, "05_gene_state_eqtl.tsv"), sep = "\t", quote = FALSE)

coupling_dt <- merge(
  apa_gene_state[gene_state != "Not detected", .(gene_name, apa_state = gene_state, apa_pairs = n_pairs, apa_detected_pairs = n_detected)],
  eqtl_gene_state[gene_state != "Not detected", .(gene_name, eqtl_state = gene_state, eqtl_pairs = n_pairs, eqtl_detected_pairs = n_detected)],
  by = "gene_name")
coupling_dt[, `:=`(apa_state = factor(apa_state, levels = gene_state_levels),
                   eqtl_state = factor(eqtl_state, levels = gene_state_levels))]
state_grid <- CJ(eqtl_state = gene_state_levels, apa_state = gene_state_levels, unique = TRUE)
state_grid[, `:=`(eqtl_state = factor(eqtl_state, levels = gene_state_levels),
                  apa_state = factor(apa_state, levels = gene_state_levels))]
coupling_matrix <- coupling_dt[, .N, by = .(eqtl_state, apa_state)]
coupling_matrix <- merge(state_grid, coupling_matrix, by = c("eqtl_state", "apa_state"), all.x = TRUE)
coupling_matrix[is.na(N), N := 0L]
coupling_total <- coupling_matrix[, sum(N)]
coupling_matrix[, global_frac := N / coupling_total]
coupling_matrix[, label := sprintf("%s; %.1f%%", formatC(N, format = "d", big.mark = ","), 100 * global_frac)]
gene_state_summary <- rbindlist(list(
  apa_gene_state[gene_state != "Not detected", .N, by = .(qtl, gene_state)],
  eqtl_gene_state[gene_state != "Not detected", .N, by = .(qtl, gene_state)]), use.names = TRUE)
gene_state_summary[, frac := N / sum(N), by = qtl]
fwrite(coupling_dt, file.path(D, "05_coupling_gene_table.tsv"), sep = "\t", quote = FALSE)
fwrite(coupling_matrix, file.path(D, "05_coupling_gene_matrix.tsv"), sep = "\t", quote = FALSE)
fwrite(gene_state_summary, file.path(D, "05_coupling_gene_state_summary.tsv"), sep = "\t", quote = FALSE)

coupling_fig <- ggplot(coupling_matrix, aes(x = apa_state, y = eqtl_state, fill = N)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = gsub("; ", "\n", label, fixed = TRUE), colour = N > max(N) / 2),
            size = pt_mm(font_sm), lineheight = 0.85, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                      labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE), name = "Count") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE, limits = rev(gene_state_levels)) +
  coord_equal() +
  labs(x = "APA gene state", y = "eQTL gene state") +
  theme_mash() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), axis.ticks = element_blank(),
        axis.line = element_blank(), legend.position = "right")
save_plot(coupling_fig, "mash_coupling", width = fig_w15, height = 3.4)

# --- downsampling replication (postnatal subsamples, n = 189) vs full postnatal and prenatal ---
downsample_files <- sort(list.files(D, pattern = "^seed_[0-9]+_post\\.cis_qtl\\.txt\\.gz$", full.names = TRUE))
full_post <- read_perm_summary(file.path(D, "postnatal.cis_qtl.txt.gz"))
setnames(full_post, c("phenotype_id", "variant_id", "qval"), c("full_phenotype_id", "full_variant_id", "full_qval"))
downsample_dt <- rbindlist(lapply(downsample_files, function(path) {
  read_perm_summary(path, as.integer(sub("^seed_([0-9]+)_.*$", "\\1", basename(path))))
}), use.names = TRUE)
prenatal <- read_perm_summary(file.path(D, "prenatal.cis_qtl.txt.gz"))

downsample_full <- merge(downsample_dt, full_post[, .(group_id, full_phenotype_id, full_variant_id, full_qval)],
                         by = "group_id", all.x = TRUE)
downsample_full[, `:=`(significant = qval < 0.05, full_significant = full_qval < 0.05,
                       same_phenotype = phenotype_id == full_phenotype_id,
                       same_variant = variant_id == full_variant_id)]
downsample_seed <- downsample_full[, .(
  n_groups = .N,
  n_significant = sum(significant, na.rm = TRUE),
  frac_significant = mean(significant, na.rm = TRUE),
  phenotype_concordance_all = mean(same_phenotype, na.rm = TRUE),
  phenotype_concordance_multifeature = mean(same_phenotype[group_size > 1], na.rm = TRUE),
  variant_concordance_all = mean(same_variant, na.rm = TRUE),
  n_full_significant_recovered = sum(significant & full_significant, na.rm = TRUE),
  frac_full_significant_recovered = sum(significant & full_significant, na.rm = TRUE) / sum(full_significant, na.rm = TRUE)
), by = seed]

winner_counts <- downsample_dt[, .N, by = .(group_id, group_size, phenotype_id)]
winner_stability <- winner_counts[order(group_id, -N, phenotype_id),
  .(modal_phenotype_id = first(phenotype_id), modal_seed_count = first(N), n_seeds_tested = sum(N),
    n_distinct_winners = .N), by = .(group_id, group_size)]
winner_stability[, `:=`(modal_fraction = modal_seed_count / n_seeds_tested,
                        frac_seeds_tested = n_seeds_tested / length(downsample_files))]
winner_stability <- merge(winner_stability,
                          unique(full_post[, .(group_id, full_phenotype_id, full_variant_id, full_qval)]),
                          by = "group_id", all.x = TRUE)
winner_stability[, modal_matches_full := modal_phenotype_id == full_phenotype_id]
fwrite(downsample_seed, file.path(D, "05_downsample_seed_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(winner_stability, file.path(D, "05_downsample_group_stability.tsv"), sep = "\t", quote = FALSE)

prenatal_sig_groups <- unique(prenatal[qval < 0.05, .(group_id, prenatal_phenotype_id = phenotype_id,
                                                      prenatal_variant_id = variant_id, prenatal_qval = qval)])
prenatal_sig_n <- nrow(prenatal_sig_groups)

# --- controlled robust core: prenatal groups replicated across downsampled postnatal seeds ---
post_sig_by_seed <- unique(downsample_dt[qval < 0.05, .(seed, group_id, postnatal_phenotype_id = phenotype_id,
                                                        postnatal_variant_id = variant_id, postnatal_qval = qval)])
replicated <- merge(post_sig_by_seed, prenatal_sig_groups, by = "group_id", all.x = FALSE, all.y = FALSE)
replication_seed <- replicated[, .(n_replicated = uniqueN(group_id), replication_rate = uniqueN(group_id) / prenatal_sig_n), by = seed]
all_seeds <- data.table(seed = sort(unique(downsample_dt$seed)))
replication_seed <- merge(all_seeds, replication_seed, by = "seed", all.x = TRUE)
replication_seed[is.na(n_replicated), `:=`(n_replicated = 0L, replication_rate = 0)]
core_groups <- replicated[, .(
  n_replicating_seeds = uniqueN(seed),
  prenatal_phenotype_id = first(prenatal_phenotype_id),
  prenatal_variant_id = first(prenatal_variant_id),
  prenatal_qval = first(prenatal_qval)), by = group_id]
core_groups[, `:=`(replicated_any = n_replicating_seeds >= 1,
                   replicated_majority = n_replicating_seeds >= ceiling(length(all_seeds$seed) / 2),
                   replicated_all = n_replicating_seeds == length(all_seeds$seed))]
core_summary <- data.table(
  metric = c("prenatal_significant_groups", "median_replicated_groups_per_seed", "median_replication_rate_per_seed",
             "core_replicated_any_seed", "core_replicated_at_least_10_of_20", "strict_core_replicated_20_of_20"),
  value = c(prenatal_sig_n, median(replication_seed$n_replicated), median(replication_seed$replication_rate),
            core_groups[replicated_any == TRUE, .N], core_groups[replicated_majority == TRUE, .N],
            core_groups[replicated_all == TRUE, .N]))
fwrite(replication_seed, file.path(D, "05_controlled_core_seed_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(core_groups, file.path(D, "05_controlled_core_groups.tsv"), sep = "\t", quote = FALSE)
fwrite(core_summary, file.path(D, "05_controlled_core_summary.tsv"), sep = "\t", quote = FALSE)

# --- downsample figure ---
downsample_gene <- merge(prenatal_sig_groups, core_groups[, .(group_id, n_replicating_seeds)], by = "group_id", all.x = TRUE)
downsample_gene[is.na(n_replicating_seeds), n_replicating_seeds := 0L]
downsample_gene[, gene_symbol := {
  idx <- match(group_id, gene_map$ensg0)
  mapped <- gene_map$gene_name[idx]
  use_map <- grepl("^ENSG", group_id) & !is.na(mapped) & nzchar(mapped)
  out <- group_id
  out[use_map] <- mapped[use_map]
  out
}]
downsample_gene[, prenatal_strength := -log10(pmax(prenatal_qval, .Machine$double.xmin))]
downsample_gene[, replication_class := fifelse(n_replicating_seeds == 0, "Prenatal-specific", "Replicated")]
replication_cdf <- data.table(k = 0:length(all_seeds$seed))
replication_cdf[, n_genes := vapply(k, function(this_k) downsample_gene[n_replicating_seeds >= this_k, .N], integer(1))]
replication_cdf[, frac_genes := n_genes / prenatal_sig_n]
label_genes <- downsample_gene[replication_class == "Prenatal-specific" & !grepl("^ENSG", gene_symbol)][order(prenatal_qval)]
label_genes <- label_genes[seq_len(min(nrow(label_genes), 5))]
fwrite(replication_cdf, file.path(D, "05_downsample_replication_cdf.tsv"), sep = "\t", quote = FALSE)
fwrite(downsample_gene[replication_class == "Prenatal-specific"], file.path(D, "05_downsample_prenatal_specific.tsv"), sep = "\t", quote = FALSE)

p_repl_cdf <- ggplot(replication_cdf, aes(x = k, y = n_genes)) +
  geom_line(linewidth = 0.45, colour = "#0F6FA8") +
  geom_point(size = 1.2, colour = "#0F6FA8") +
  geom_vline(xintercept = c(10, 20), linetype = "dashed", colour = "grey65", linewidth = 0.3) +
  scale_x_continuous(breaks = seq(0, length(all_seeds$seed), by = 5)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06)), labels = scales::comma) +
  labs(x = "Replicating postnatal n=189 seeds (>= k)", y = "Prenatal QTL genes", title = "Replication frequency") +
  theme_mash()
p_repl_scatter <- ggplot(downsample_gene, aes(x = prenatal_strength, y = n_replicating_seeds, colour = replication_class)) +
  geom_jitter(width = 0, height = 0.2, size = 0.9, alpha = 0.75) +
  geom_label_repel(data = label_genes, aes(x = prenatal_strength, y = n_replicating_seeds, label = gene_symbol),
                   inherit.aes = FALSE, nudge_x = 3, nudge_y = 2, direction = "both", box.padding = 0.5,
                   point.padding = 0.5, min.segment.length = 0, segment.color = "#B2182B", segment.size = 0.2,
                   force = 4, size = pt_mm(font_sm), linewidth = 0.15, label.padding = grid::unit(1.0, "pt"),
                   fill = "white", colour = "#B2182B", seed = 1, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c("Prenatal-specific" = "#B2182B", "Replicated" = "grey55"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  scale_y_continuous(breaks = seq(0, length(all_seeds$seed), by = 5),
                     limits = c(-0.5, length(all_seeds$seed) + 0.5)) +
  labs(x = "Prenatal QTL strength (-log10 q)", y = "Replicating seeds", title = "Prenatal-specific candidates") +
  theme_mash() + theme(legend.position = "bottom")
core_fig <- (p_repl_cdf | p_repl_scatter) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = font_tag))
save_plot(core_fig, "mash_downsample", width = fig_w2, height = 3.0)
