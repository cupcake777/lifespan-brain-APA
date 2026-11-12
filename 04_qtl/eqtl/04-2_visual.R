# Mediation overview figure (colocalization rate, expression attenuation, apaQTL retention).
# Reads data/04_qtl mediation and coloc summaries; writes figures/04_qtl/04_mediation_overview.{pdf,png}.
suppressMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
D <- file.path(APA_ROOT, "data", "04_qtl")
fig_dir <- file.path(APA_ROOT, "figures", "04_qtl")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

stage_key <- c(prenatal = "Prenatal", postnatal = "Postnatal")
call_cols <- c("APA -> EXP" = "#1F78B4", "EXP -> APA" = "#E31A1C", "Partial mediation" = "#33A02C", "No mediation" = "#E0E0E0")
stage_cols <- c("Prenatal" = "#E41A1C", "Postnatal" = "#2B83BA")

theme_pub <- function(base_size = 10.5) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = base_size - 1),
      axis.title = element_text(colour = "black", size = base_size, face = "bold"),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.key.size = grid::unit(0.35, "cm"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0, margin = margin(b = 6))
    )
}

# --- data ---
med <- fread(file.path(D, "04_mediation_results.tsv"))
med[, stage_lab := factor(stage_key[stage], levels = c("Prenatal", "Postnatal"))]
coloc <- fread(file.path(D, "coloc_pp4_definition_sum.tsv"))
coloc[, stage_lab := factor(stage_key[stage], levels = c("Prenatal", "Postnatal"))]
med[, call_clean := fcase(
  mediation_call %in% c("APA -> expression", "APA -> EXP"), "APA -> EXP",
  mediation_call %in% c("Expression -> APA", "EXP -> APA"), "EXP -> APA",
  mediation_call %in% c("Independent / partial", "Partial mediation") & fdr_ab_fwd < 0.05, "Partial mediation",
  default = "No mediation")]
med[, call_clean := factor(call_clean, levels = names(call_cols))]

# --- a: colocalization rate ---
coloc[, pct := 100 * fraction]
coloc[, count_label := sprintf("%s/%s (%.1f%%)", format(numerator, big.mark = ",", trim = TRUE),
                               format(denominator, big.mark = ",", trim = TRUE), pct)]
pA <- ggplot(coloc, aes(x = pct, y = stage_lab, fill = stage_lab)) +
  geom_col(width = 0.52, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = count_label), hjust = 1.08, size = 3.1, colour = "white", fontface = "bold") +
  scale_fill_manual(values = stage_cols, guide = "none") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 32), breaks = seq(0, 30, 10), labels = function(x) paste0(x, "%")) +
  labs(x = "Colocalized apaQTL-eQTL pairs", y = NULL, title = "a   Colocalization rate") +
  theme_pub() + theme(axis.text.y = element_text(face = "bold"))

# --- pairs with mediating evidence ---
med_effect <- med[call_clean %in% c("APA -> EXP", "Partial mediation") &
                    is.finite(c_total) & is.finite(c_prime) & is.finite(a) & is.finite(d_prime)]

# --- b: eQTL effect attenuation (total |c| vs conditioned |c'|) ---
max_val_b <- max(c(abs(med_effect$c_total), abs(med_effect$c_prime)), na.rm = TRUE) * 1.05
slope_b <- med_effect[, .(slope = sum(abs(c_total) * abs(c_prime)) / sum(abs(c_total)^2)), by = stage_lab]
slope_b[, `:=`(x = max_val_b * 0.03, y = max_val_b * fifelse(stage_lab == "Prenatal", 0.97, 0.88),
               intercept = 0, label = sprintf("beta==%.2f", slope))]
pB <- ggplot(med_effect, aes(x = abs(c_total), y = abs(c_prime), colour = stage_lab)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.65, size = 1.8, stroke = 0) +
  geom_abline(data = slope_b, aes(intercept = intercept, slope = slope, colour = stage_lab), linewidth = 0.8, show.legend = FALSE) +
  geom_text(data = slope_b, aes(x = x, y = y, label = label, colour = stage_lab), hjust = 0, vjust = 1,
            size = 2.8, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = stage_cols, name = "Stage") +
  coord_fixed(ratio = 1, xlim = c(0, max_val_b), ylim = c(0, max_val_b)) +
  labs(x = "Total eQTL effect (|c|)", y = "eQTL effect conditioned on APA (|c'|)", title = "b   Expression attenuation") +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank())

# --- c: apaQTL effect retention (total |a| vs conditioned |d'|) ---
max_val_c <- max(c(abs(med_effect$a), abs(med_effect$d_prime)), na.rm = TRUE) * 1.05
slope_c <- med_effect[, .(slope = sum(abs(a) * abs(d_prime)) / sum(abs(a)^2)), by = stage_lab]
slope_c[, `:=`(x = max_val_c * 0.03, y = max_val_c * fifelse(stage_lab == "Prenatal", 0.97, 0.88),
               intercept = 0, label = sprintf("beta==%.2f", slope))]
pC <- ggplot(med_effect, aes(x = abs(a), y = abs(d_prime), colour = stage_lab)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.65, size = 1.8, stroke = 0) +
  geom_abline(data = slope_c, aes(intercept = intercept, slope = slope, colour = stage_lab), linewidth = 0.8, show.legend = FALSE) +
  geom_text(data = slope_c, aes(x = x, y = y, label = label, colour = stage_lab), hjust = 0, vjust = 1,
            size = 2.8, fontface = "bold", parse = TRUE, show.legend = FALSE) +
  scale_colour_manual(values = stage_cols, name = "Stage") +
  coord_fixed(ratio = 1, xlim = c(0, max_val_c), ylim = c(0, max_val_c)) +
  labs(x = "Total apaQTL effect (|a|)", y = "apaQTL effect conditioned on expression (|d'|)",
       title = "c   Independent APA retention") +
  theme_pub() + theme(legend.position = "bottom", legend.title = element_blank())

layout <- (pA | pB | pC) + plot_layout(widths = c(0.9, 1.0, 1.0), guides = "collect") +
  plot_annotation(caption = "Colocalization: PP4 >= 0.7; mediation subset: PP4 >= 0.8.",
                  theme = theme(plot.caption = element_text(size = 9, hjust = 0))) &
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(fig_dir, "04_mediation_overview.pdf"), layout, width = 11.2, height = 4.45, device = grDevices::cairo_pdf)
ggsave(file.path(fig_dir, "04_mediation_overview.png"), layout, width = 11.2, height = 4.45, dpi = 300, bg = "white")
