# Nature Genetics theme, palette and font setup shared by all figure scripts.
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales); library(grid)
})

ROOT    <- Sys.getenv("APA_ROOT")
FIG_DIR <- file.path(ROOT, "figures")

# --- font ---
# Arial lives in APA_FONT_DIR, outside the system font path; register it explicitly.
.fdir <- Sys.getenv("APA_FONT_DIR")
if (dir.exists(.fdir)) {
  .pf <- function(f) file.path(.fdir, f)
  if (all(file.exists(c(.pf("arial.ttf"), .pf("arialbd.ttf"), .pf("ariali.ttf"), .pf("arialbi.ttf"))))) {
    try(systemfonts::register_font("Arial", plain = .pf("arial.ttf"), bold = .pf("arialbd.ttf"),
                                   italic = .pf("ariali.ttf"), bolditalic = .pf("arialbi.ttf")), silent = TRUE)
  }
}
.fams <- tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character())
has_font <- function(fam) {
  m <- try(systemfonts::match_fonts(fam), silent = TRUE)
  fam %in% .fams || (!inherits(m, "try-error") && nrow(m) > 0)
}
FONT <- if (has_font("Arial")) "Arial" else if (has_font("Helvetica")) "Helvetica" else
  if ("Nimbus Sans" %in% .fams) "Nimbus Sans" else if ("Liberation Sans" %in% .fams) "Liberation Sans" else "sans"

# --- palette ---
INK           <- "#222222"
MUTED         <- "#5A5A5A"
NEUTRAL       <- "#7F7F7F"
NEUTRAL_LIGHT <- "#D9D9D9"
GRID_COL      <- "#EBEBEB"

C1_COLOR     <- "#4A7BA7"
C2_COLOR     <- "#D17A4A"
SHARED_COLOR <- "#7A9E9F"
CLUSTER_COLORS <- c(C1 = C1_COLOR, C2 = C2_COLOR, Shared = SHARED_COLOR, Other = NEUTRAL_LIGHT)

STAGE_COLORS <- c(
  Prenatal = "#E3A73F", Postnatal = "#6BAED6",
  "Prenatal-specific" = "#E3A73F", "Postnatal-specific" = "#6BAED6", Shared = "#7F7F7F"
)
STAGE_TEXT <- c(Prenatal = "#B07D22", Postnatal = "#3C84B8")

CLASS_COLORS <- c(
  Reversed = "#D17A4A", Flip = "#D17A4A", Amplified = "#4A7BA7", Amplifying = "#4A7BA7",
  Attenuated = "#7A9E9F", Attenuating = "#7A9E9F"
)
CLASS_TINT <- c(Reversed = "#F6E4D9", Amplified = "#DCE6EF", Attenuated = "#E3ECEC")

LAYER_COLORS <- c(APA = "#009E73", RNA = "#0072B2", Protein = "#D55E00", protein = "#D55E00")
DIVERGING_RAMP <- c("#2166AC", "#F7F7F7", "#B2182B")

# --- theme ---
theme_nature <- function(base_size = 6.5, base_family = FONT) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text               = element_text(colour = INK, size = base_size, family = base_family),
      axis.title         = element_text(colour = INK, size = base_size, face = "plain"),
      axis.title.x       = element_text(margin = margin(t = 2.5, unit = "pt")),
      axis.title.y       = element_text(margin = margin(r = 2.5, unit = "pt"), angle = 90),
      axis.text          = element_text(colour = INK, size = base_size - 1),
      axis.text.x        = element_text(margin = margin(t = 1.5, unit = "pt")),
      axis.text.y        = element_text(margin = margin(r = 1.5, unit = "pt")),
      axis.line          = element_line(linewidth = 0.25, colour = INK, lineend = "square"),
      axis.ticks         = element_line(linewidth = 0.25, colour = INK),
      axis.ticks.length  = unit(0.8, "mm"),
      legend.title       = element_text(size = base_size - 0.5, colour = INK),
      legend.text        = element_text(size = base_size - 1, colour = INK),
      legend.key.size    = unit(2.2, "mm"),
      legend.background  = element_blank(),
      legend.box.background = element_blank(),
      legend.key         = element_blank(),
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "mm"),
      strip.background   = element_blank(),
      strip.text         = element_text(size = base_size, colour = INK, face = "plain",
                                        margin = margin(b = 2, unit = "pt")),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      plot.margin        = margin(2, 2.5, 2, 2, "mm"),
      panel.background   = element_blank(),
      plot.background    = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank()
    )
}
theme_set(theme_nature())

update_geom_defaults("point",    list(stroke = 0.2, size = 0.8))
update_geom_defaults("line",     list(linewidth = 0.4))
update_geom_defaults("segment",  list(linewidth = 0.25, colour = INK))
update_geom_defaults("errorbar", list(linewidth = 0.25, colour = INK))
update_geom_defaults("text",     list(family = FONT, size = 5.5 / .pt, colour = INK))
