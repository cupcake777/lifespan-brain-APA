# Shared theme, helpers, colors and loader for module 02 figures.
# Paths come from APA_ROOT (source config/paths.local.sh first).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# --- helpers ---

theme_pub <- function(base = 11) {
  theme_classic(base_size = base) +
    theme(
      axis.text = element_text(color = "black", face = "bold"),
      axis.title = element_text(color = "black", face = "bold"),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.3),
      plot.title = element_text(size = base + 1, face = "bold"),
      plot.subtitle = element_blank(),
      plot.caption = element_text(size = base - 1.5, color = "#4D4D4D",
                                  hjust = 0, face = "bold"),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(size = base, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = base - 1, face = "bold"),
      legend.title = element_text(size = base - 1, face = "bold")
    )
}

save_png <- function(plot_obj, filename, width, height) {
  out_file <- file.path(PANEL_DIR, filename)
  ggsave(out_file, plot_obj, width = width, height = height, dpi = 300, bg = "white")
  out_file
}

save_pdf <- function(plot_obj, filename, width, height) {
  out_file <- file.path(PANEL_DIR, filename)
  ggsave(out_file, plot_obj, width = width, height = height,
         device = cairo_pdf, bg = "white")
  out_file
}

cluster_label <- function(x) {
  ifelse(as.integer(x) == 1L, "C1 lengthening", "C2 shortening")
}

safe_spearman <- function(a, b) {
  mask <- is.finite(a) & is.finite(b)
  if (sum(mask) < 5L || sd(a[mask]) == 0 || sd(b[mask]) == 0) return(NA_real_)
  suppressWarnings(cor(a[mask], b[mask], method = "spearman"))
}

load_sites <- function() {
  sites <- fread(file.path(OUT_DIR, "02_sites.tsv"))
  if (!"dynamic_pass" %in% names(sites)) {
    null_file <- file.path(INTERMEDIATE_DIR, "03_null.tsv")
    if (file.exists(null_file)) {
      null <- fread(null_file)
      delta_thr <- null[
        seed == 2026L & quantile == 0.95,
        delta_threshold
      ][1]
    } else {
      delta_thr <- 0.3948
      message("warning: 03_null.tsv not found, using fallback threshold ", delta_thr)
    }
    sites[, dynamic_pass := association_pass == TRUE &
                              is.finite(delta_logit) &
                              abs(delta_logit) >= delta_thr]
    sites[, delta_logit_threshold := delta_thr]
  }
  sites
}

# --- paths ---

APA_ROOT <- Sys.getenv("APA_ROOT")
if (!nzchar(APA_ROOT)) stop("APA_ROOT is not set; source config/paths.local.sh first")
OUT_DIR <- file.path(APA_ROOT, "data", "02_traj")
INTERMEDIATE_DIR <- OUT_DIR
PANEL_DIR <- file.path(APA_ROOT, "figures", "02_traj")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

# --- constants ---

FDR_THR <- 0.05
SEED <- 7L

C_LABEL_COLORS <- c("C1" = "#0072B2", "C2" = "#D55E00")
TRAJ_COLORS <- c("C1 lengthening" = "#0072B2", "C2 shortening" = "#D55E00")

age_breaks <- log10(c(0, 5, 10, 20, 50, 100) + 1)
age_labels <- c("birth", "5", "10", "20", "50", "100")

set.seed(SEED)
