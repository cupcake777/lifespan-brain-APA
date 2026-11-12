# Interaction maxCPP S-LDSC: collect regressions, FDR (within trait and global), tau* per trait.
# Reads data/04_qtl/interaction_sldsc/; writes interaction_maxcpp_tau_star.tsv to data/04_qtl/.
suppressMessages(library(data.table))

# --- config ---
D <- file.path(Sys.getenv("APA_ROOT"), "data", "04_qtl")
sldsc_dir <- file.path(D, "interaction_sldsc")

traits <- c("EAS_AD", "EAS_AIS", "EAS_ALS", "EAS_BIP", "EAS_EP", "EAS_PD", "EAS_SCZ",
            "EUR_AD", "EUR_ADHD", "EUR_ANX", "EUR_ASD", "EUR_BIP", "EUR_EP", "EUR_INS",
            "EUR_INT", "EUR_MDD", "EUR_PD", "EUR_PTSD", "EUR_SCZ")

# --- h2 from regression logs ---
h2_dt <- rbindlist(lapply(traits, function(trait) {
  pop <- sub("_.*$", "", trait)
  f <- file.path(sldsc_dir, "logs", sprintf("h2.interaction_maxcpp.%s.%s.log", pop, trait))
  txt <- readLines(f, warn = FALSE)
  hit <- grep("^Total Observed scale h2:", txt, value = TRUE)
  if (!length(hit)) return(NULL)
  m <- regmatches(hit[1], regexec("Total Observed scale h2:\\s*([0-9.eE+-]+)\\s*\\(([0-9.eE+-]+)\\)", hit[1]))[[1]]
  if (length(m) < 3L) return(NULL)
  data.table(population = pop, trait = trait, h2 = as.numeric(m[2]), h2_se = as.numeric(m[3]))
}))
h2_dt[, well_powered := (h2 / h2_se) >= 7.0]

# --- regression results ---
result_files <- list.files(file.path(sldsc_dir, "results", "interaction_maxcpp"),
                           pattern = "\\.results$", recursive = TRUE, full.names = TRUE)
stopifnot(length(result_files) == 19L)
raw <- rbindlist(lapply(result_files, function(f) {
  data.table(model = "interaction_maxcpp", population = basename(dirname(f)),
             stage = "interaction", trait = sub("\\.results$", "", basename(f)), fread(f))
}))
setnames(raw, "Coefficient_z-score", "coef_z", skip_absent = TRUE)
qtl <- raw[grepl("^L2_[0-9]+$", Category)]
qtl[, annotation := "interactionQTL"]
stopifnot(nrow(qtl) == 19L, !anyDuplicated(qtl[, .(population, trait)]),
          uniqueN(qtl$Category) == 1L, qtl$Category[1] == "L2_1")

qtl[, coef_p_1s := pnorm(coef_z, lower.tail = FALSE)]
qtl[, coef_fdr := p.adjust(coef_p_1s, method = "BH")]
qtl[, coef_fdr_global := p.adjust(coef_p_1s, method = "BH")]
qtl[, coef_fdr_trait := p.adjust(coef_p_1s, method = "BH"), by = .(population, trait)]
qtl[, n_tests_family := .N]

sd_dt <- fread(file.path(D, "annot_sd_interaction.tsv"))[
  , .(population, stage, annotation, sd_annot = sd, M_common = n_snps)]

qtl <- merge(qtl, h2_dt, by = c("population", "trait"), all.x = TRUE)
qtl <- merge(qtl, sd_dt, by = c("population", "stage", "annotation"), all.x = TRUE)
stopifnot(!anyNA(qtl[, .(h2, h2_se, M_common, sd_annot)]))

qtl[, tau_star := fifelse(h2 > 0, Coefficient * M_common * sd_annot / h2, NA_real_)]
qtl[, tau_star_se := fifelse(h2 > 0, Coefficient_std_error * M_common * sd_annot / h2, NA_real_)]
qtl[, `:=`(tau_lo = tau_star - 1.96 * tau_star_se, tau_hi = tau_star + 1.96 * tau_star_se)]

keep_cols <- c("population", "stage", "annotation", "trait", "model", "Category",
               "Prop._SNPs", "Prop._h2", "Prop._h2_std_error", "Enrichment",
               "Enrichment_std_error", "Enrichment_p", "Coefficient", "Coefficient_std_error",
               "coef_z", "coef_p_1s", "coef_fdr", "coef_fdr_global", "coef_fdr_trait",
               "n_tests_family", "h2", "h2_se", "well_powered", "sd_annot", "M_common",
               "tau_star", "tau_star_se", "tau_lo", "tau_hi")
fwrite(qtl[order(population, trait), ..keep_cols],
       file.path(D, "interaction_maxcpp_tau_star.tsv"), sep = "\t")
