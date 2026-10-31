# Fig 5 composite: a-c from fig5_abc.R; d (SCAF1, prenatal) and e (KLC1, postnatal) from fig5_de.R.
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/main_figures/figure_5_overview.R"))
source(file.path(Sys.getenv("APA_CODE_ROOT"), "figure_generation/main_figures/figure_5_loci.R"))
tg <- function(p, t) p + labs(tag = t) + theme(plot.tag = element_text(size = 8, face = "bold"), plot.tag.position = c(0, 1))

# max colocalization PP4 for one gene x GWAS x stage (repaired run; co from fig5_abc.R)
pp4_of <- function(g, gwas, qtl) {
  s <- co[GWAS_ID == gwas & QTL_ID == qtl & gene == g]
  stopifnot(nrow(s) > 0)
  max(s$PP4)
}
pp_scaf <- pp4_of("SCAF1", "EUR_SCZ", "prenatal")
pp_klc  <- pp4_of("KLC1", "EUR_SCZ", "postnatal")
dd <- locus_panel("SCAF1", "EUR SCZ", stage = "prenatal", pp4 = pp_scaf)
ee <- locus_panel("KLC1", "EUR SCZ", stage = "postnatal", pp4 = pp_klc)
loc <- function(L, t) (tg(L$top, t) / L$gm) + plot_layout(heights = c(3, 0.16 * L$nr + 0.25))
top <- (tg(pa, "a") | tg(pb, "b")) + plot_layout(widths = c(80, 100))
bot <- (loc(dd, "d") | loc(ee, "e")) + plot_layout(widths = c(1, 1))
fig <- wrap_elements(full = top) / wrap_elements(full = tg(pc, "c")) / wrap_elements(full = bot) + plot_layout(heights = c(68, 50, 84))
ggsave(file.path(FIG_DIR, "fig5.png"), fig, width = 180, height = 202, units = "mm", dpi = 600, bg = "white", device = ragg::agg_png)
ggsave(file.path(FIG_DIR, "fig5.pdf"), fig, width = 180, height = 202, units = "mm", device = cairo_pdf)
