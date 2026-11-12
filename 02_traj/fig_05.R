# Supported-grid clustering diagnostics: k selection, bootstrap consensus,
# cluster silhouette and effect direction.
suppressPackageStartupMessages({
  library(patchwork)
})
source("fig_common.R")

# --- load data ---

sites <- fread(file.path(OUT_DIR, "02_sites_clustered.tsv"))
clustered <- sites[directional_pass == TRUE & !is.na(cluster_id)]
k_diag <- fread(file.path(OUT_DIR, "03_k_diagnostics.tsv"))
k_diag[, selected := auto_select == TRUE]
boot_label_file <- file.path(INTERMEDIATE_DIR, "03_consensus_annotation.tsv")
boot_iter <- 100L
boot_iter_label <- comma(boot_iter)

# --- panel A: k selection ---

p_k <- ggplot(k_diag, aes(k, sample_silhouette)) +
  geom_line(color = "#333333", linewidth = 0.6) +
  geom_point(aes(fill = selected), shape = 21, size = 2.6, color = "black", stroke = 0.35) +
  geom_text(data = k_diag[selected == TRUE],
            aes(x = k + 0.18, y = sample_silhouette + 0.025,
                label = paste0("selected k=", k)),
            hjust = 0, size = 3.0, color = "#D55E00", fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "white"), guide = "none") +
  scale_x_continuous(limits = c(1.8, 8.2), breaks = k_diag$k) +
  labs(x = "Candidate k", y = "Sample silhouette", title = "K selection") +
  theme_pub(base = 11)

# --- panel B: bootstrap consensus ---

n_total <- nrow(clustered)
n_cons1 <- sum(clustered$cluster_consensus == 1.0)
n_cons_lt1 <- sum(clustered$cluster_consensus < 1.0)
consensus_summary <- data.table(
  status = factor(c("Consensus = 1", "Variable"), levels = c("Consensus = 1", "Variable")),
  n = c(n_cons1, n_cons_lt1)
)
consensus_summary[, pct := 100 * n / n_total]
consensus_summary[, label_y := fifelse(status == "Consensus = 1", 88, 6)]
consensus_summary[, label := sprintf("%s\n%.2f%%", scales::comma(n), pct)]

p_cons <- ggplot(consensus_summary, aes(status, pct, fill = status)) +
  geom_col(width = 0.58, color = NA) +
  geom_text(aes(y = label_y, label = label), size = 3.05, color = "black",
            fontface = "bold", lineheight = 0.95) +
  scale_fill_manual(values = c("Consensus = 1" = "#0072B2", "Variable" = "#D55E00"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100),
                     labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Sites", title = sprintf("Bootstrap consensus\nn=%s, %s iter.",
                                              scales::comma(n_total),
                                              boot_iter_label)) +
  theme_pub(base = 11)

# --- panel C: silhouette by cluster ---

p_sil <- ggplot(clustered[is.finite(cluster_silhouette)],
                aes(paste0("C", cluster_id), cluster_silhouette,
                    fill = paste0("C", cluster_id))) +
  geom_violin(alpha = 0.75, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = C_LABEL_COLORS, guide = "none") +
  labs(x = NULL, y = "Silhouette width", title = "Cluster separation") +
  theme_pub(base = 11)

# --- panel D: delta logit by cluster ---

p_delta <- ggplot(clustered[is.finite(delta_logit)],
                  aes(paste0("C", cluster_id), delta_logit,
                      fill = paste0("C", cluster_id))) +
  geom_violin(alpha = 0.75, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.10, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.3) +
  scale_fill_manual(values = C_LABEL_COLORS, guide = "none") +
  labs(x = NULL, y = expression(Delta * "logit(PDUI)"), title = "Effect direction") +
  theme_pub(base = 11)

# --- assemble ---

p05 <- (p_k | p_cons | p_sil | p_delta) +
  plot_layout(widths = c(1.1, 1.0, 1.1, 1.1))
save_png(p05, "05_clu_consensus.png", 13.2, 4.2)
save_pdf(p05, "05_clu_consensus.pdf", 13.2, 4.2)

layout_validation <- rbindlist(list(
  data.table(
    figure = "05_clu_consensus",
    panel = "A_k_selection",
    n_rows = nrow(k_diag),
    n_labels = sum(k_diag$selected, na.rm = TRUE),
    max_label_chars = max(nchar(paste0("selected k=", k_diag[selected == TRUE]$k))),
    validation_note = sprintf(
      "K selection uses sample silhouette across k=%s-%s; selected k=%s.",
      min(k_diag$k), max(k_diag$k), paste(k_diag[selected == TRUE]$k, collapse = ",")
    )
  ),
  data.table(
    figure = "05_clu_consensus",
    panel = "B_bootstrap_consensus",
    n_rows = nrow(consensus_summary),
    n_labels = nrow(consensus_summary),
    max_label_chars = max(nchar(consensus_summary$label)),
    validation_note = sprintf(
      "Consensus summary is computed from %s bootstrap iterations inferred from %s.",
      boot_iter_label, basename(boot_label_file)
    )
  ),
  data.table(
    figure = "05_clu_consensus",
    panel = "C_cluster_separation",
    n_rows = nrow(clustered[is.finite(cluster_silhouette)]),
    n_labels = 2L,
    max_label_chars = 2L,
    validation_note = "Cluster separation panel shows silhouette distributions for C1 and C2."
  ),
  data.table(
    figure = "05_clu_consensus",
    panel = "D_effect_direction",
    n_rows = nrow(clustered[is.finite(delta_logit)]),
    n_labels = 2L,
    max_label_chars = 2L,
    validation_note = "Effect direction panel shows delta-logit distributions for C1 and C2 with a zero reference line."
  )
), fill = TRUE)
fwrite(layout_validation,
       file.path(OUT_DIR, "05_clu_consensus_layout_validation.tsv"),
       sep = "\t")
