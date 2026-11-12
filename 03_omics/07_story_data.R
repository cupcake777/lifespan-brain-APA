#!/usr/bin/env Rscript
# Story-figure inputs: nested model comparison (delta R2) with NAPRT control,
# NAPRT fixed-spline curves, INPP4A motif track and brain-ranked miRNA families.

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(splines)
})
setDTthreads(1L)

# --- config ---
apa_root <- Sys.getenv("APA_ROOT")
omics_dir <- file.path(apa_root, "data", "03_omics")
ref_dir <- file.path(apa_root, "ref")
motif_seq_dir <- file.path(apa_root, "raw", "motif_seq")
primary_df <- fread(file.path(omics_dir, "02_df_calibration.tsv"))$fixed_df[1]

ranked_path <- file.path(omics_dir, "03_decoupling_ranked_pairs.tsv")
delta_r2_output_path <- file.path(omics_dir, "07_model_comp_delta_r2.tsv")
motif_regions_output_path <- file.path(omics_dir, "07_motif_track_regions.tsv")
motif_sites_output_path <- file.path(omics_dir, "07_motif_track_sites.tsv")
mirna_rank_output_path <- file.path(omics_dir, "07_inpp4a_mirna_brain_rank.tsv")
mirna_atlas_path <- file.path(omics_dir, "miRNATissueAtlas2_human_miRNA_rpmm.parquet")
naprt_control_curve_output_path <- file.path(omics_dir, "07_naprt_ctl_curves.tsv")
n_mirna_labels <- 3L

# motif track candidate regions: primary ranked site per core gene
motif_regions <- data.table(
  gene_symbol = c("PIK3CA", "KPNA4", "INPP4A", "ARIH1", "RHOT1"),
  site_id = c(
    "ENST00000675467.1|PIK3CA|chr3:179234094-179235371|+",
    "ENST00000334256.9|KPNA4|chr3:160495007-160502202|-",
    "ENST00000463367.1|INPP4A|chr2:98520687-98521766|+",
    "ENST00000563310.1|ARIH1|chr15:72583208-72586559|+",
    "ENST00000354266.7|RHOT1|chr17:32224616-32225666|+"
  ),
  chrom = c("chr3", "chr3", "chr2", "chr15", "chr17"),
  locus_start = c(179234094L, 160495007L, 98520687L, 72583208L, 32224616L),
  locus_end = c(179235371L, 160502202L, 98521766L, 72586559L, 32225666L),
  strand = c("+", "-", "+", "+", "+"),
  fasta_file = c("PIK3CA_P01", "KPNA4_P02", "INPP4A_P03", "ARIH1_P04", "RHOT1_P05")
)
naprt_control_site <- "ENST00000340490.7|NAPRT|chr8:143574787-143575093|-"
flank <- 500L
pas_hexamers <- c("AATAAA", "ATTAAA", "AGTAAA", "TATAAA", "CATAAA", "GATAAA",
                  "AATATA", "AATACA", "AATAGA", "AAAAAG", "ACTAAA", "AACAAA")

# --- helpers ---
revcomp <- function(x) {
  chartr("ACGTUacgtu", "TGCAAtgcaa", paste(rev(strsplit(x, "")[[1]]), collapse = ""))
}
find_all <- function(pattern, subject) {
  hit <- gregexpr(pattern, subject, perl = TRUE)[[1]]
  if (hit[1] == -1) integer(0) else as.integer(hit)
}

# --- load data ---
manifest <- fread(file.path(omics_dir, "matched_sample_manifest.tsv"))
ranked <- fread(ranked_path)
setorder(ranked, decoupling_rank)
pair_source <- fread(file.path(omics_dir, "03_site_protein_pairs.tsv"), select = c(
  "site_id", "gene_symbol", "gencode_gene_base", "protein_accession",
  "analysis_df", "analysis_role"
))
control_pair <- pair_source[analysis_df == primary_df & site_id == naprt_control_site]
motif_regions[, pair_rank := ranked$decoupling_rank[match(site_id, ranked$site_id)]]
stopifnot(!anyNA(motif_regions$pair_rank), nrow(control_pair) == 1L)
universe <- fread(file.path(omics_dir, "01_analysis_universe.tsv"))
apa <- readRDS(file.path(omics_dir, "apa_matched.rds"))
rna <- readRDS(file.path(omics_dir, "rna_matched.rds"))
protein <- readRDS(file.path(omics_dir, "protein_matched.rds"))
stopifnot(
  identical(colnames(apa), manifest$sample_id),
  identical(colnames(rna), manifest$sample_id),
  identical(colnames(protein), manifest$sample_id)
)
rna_gene_base <- sub("\\.[0-9]+$", "", as.character(rownames(rna)))
age_model <- log10(pmax(manifest$age_years + 1, 1e-6))
age_basis <- ns(age_model, df = primary_df)

# --- nested model comparison: protein ~ RNA + age vs protein ~ APA + RNA + age ---
n_samples <- nrow(manifest)
model_comparison <- rbindlist(lapply(seq_len(nrow(ranked)), function(index) {
  apa_value <- as.numeric(apa[ranked$site_id[index], ])
  rna_value <- as.numeric(rna[match(ranked$gencode_gene_base[index], rna_gene_base), ])
  protein_value <- as.numeric(protein[ranked$protein_accession[index], ])
  reduced <- lm.fit(cbind("(Intercept)" = 1, RNA = rna_value, age_basis), protein_value)
  full <- lm.fit(
    cbind("(Intercept)" = 1, APA = apa_value, RNA = rna_value, age_basis),
    protein_value
  )
  sst <- sum((protein_value - mean(protein_value))^2)
  sse_reduced <- sum(reduced$residuals^2)
  sse_full <- sum(full$residuals^2)
  df_full <- n_samples - full$rank
  f_stat <- (sse_reduced - sse_full) / (sse_full / df_full)
  data.table(
    gene_symbol = ranked$gene_symbol[index],
    site_id = ranked$site_id[index],
    protein_accession = ranked$protein_accession[index],
    decoupling_rank = ranked$decoupling_rank[index],
    r2_rna_age = 1 - sse_reduced / sst,
    r2_apa_rna_age = 1 - sse_full / sst,
    delta_r2 = (sse_reduced - sse_full) / sst,
    lrt_f = f_stat,
    lrt_df = df_full,
    lrt_p = pf(f_stat, 1, df_full, lower.tail = FALSE),
    beta_apa_full = full$coefficients[2L],
    beta_rna_full = full$coefficients[3L]
  )
}))
model_comparison[, lrt_q_bh := p.adjust(lrt_p, method = "BH")]

control_apa <- as.numeric(apa[naprt_control_site, ])
control_rna <- as.numeric(rna[match(control_pair$gencode_gene_base, rna_gene_base), ])
control_protein <- as.numeric(protein[control_pair$protein_accession, ])
control_reduced <- lm.fit(cbind("(Intercept)" = 1, RNA = control_rna, age_basis),
                          control_protein)
control_full <- lm.fit(cbind("(Intercept)" = 1, APA = control_apa,
                             RNA = control_rna, age_basis), control_protein)
control_sst <- sum((control_protein - mean(control_protein))^2)
control_sse_reduced <- sum(control_reduced$residuals^2)
control_sse_full <- sum(control_full$residuals^2)
control_df_full <- n_samples - control_full$rank
control_f <- (control_sse_reduced - control_sse_full) /
  (control_sse_full / control_df_full)
control_model_comparison <- data.table(
  gene_symbol = control_pair$gene_symbol,
  site_id = naprt_control_site,
  protein_accession = control_pair$protein_accession,
  decoupling_rank = NA_integer_,
  r2_rna_age = 1 - control_sse_reduced / control_sst,
  r2_apa_rna_age = 1 - control_sse_full / control_sst,
  delta_r2 = (control_sse_reduced - control_sse_full) / control_sst,
  lrt_f = control_f,
  lrt_df = control_df_full,
  lrt_p = pf(control_f, 1, control_df_full, lower.tail = FALSE),
  lrt_q_bh = NA_real_,
  beta_apa_full = control_full$coefficients[2L],
  beta_rna_full = control_full$coefficients[3L]
)
model_comparison <- rbindlist(list(model_comparison, control_model_comparison),
                              use.names = TRUE)
fwrite(model_comparison, delta_r2_output_path, sep = "\t", quote = FALSE)

# --- NAPRT control: fixed-spline trajectories ---
support <- fread(file.path(omics_dir, "01_support_profile.tsv"))[is_supported == TRUE]
setorder(support, grid_index)
curve_model_dat <- data.table(
  age_model = age_model,
  sex_binary = as.integer(manifest$sex == "M"),
  rin_z = as.numeric(scale(manifest$rna_rin))
)
curve_grid <- data.table(
  grid_index = support$grid_index,
  age_model = support$age_model,
  age_years = support$grid_age,
  sex_binary = 0,
  rin_z = 0
)
control_features <- list(
  APA = naprt_control_site,
  RNA = control_pair$gencode_gene_base,
  protein = control_pair$protein_accession
)
control_curves <- rbindlist(lapply(names(control_features), function(layer) {
  matrix <- switch(layer, APA = apa, RNA = rna, protein = protein)
  feature_id <- control_features[[layer]]
  feature_index <- if (layer == "RNA") {
    match(feature_id, rna_gene_base)
  } else {
    match(feature_id, rownames(matrix))
  }
  stopifnot(!is.na(feature_index))
  dat <- copy(curve_model_dat)
  dat[, y := as.numeric(matrix[feature_index, ])]
  formula <- if (layer == "protein") {
    y ~ ns(age_model, df = primary_df)
  } else {
    y ~ sex_binary + rin_z + ns(age_model, df = primary_df)
  }
  fit <- lm(formula, data = dat)
  prediction <- predict(fit, newdata = curve_grid, se.fit = TRUE)
  fitted <- as.numeric(prediction$fit)
  se <- as.numeric(prediction$se.fit)
  data.table(
    gene_symbol = "NAPRT",
    layer = layer,
    feature_id = feature_id,
    decoupling_ranks = "control",
    curve_grid[, .(grid_index, age_model, age_years)],
    model = sprintf("fixed_ns_df%d", primary_df),
    fitted = fitted,
    se = se,
    lower_95 = fitted - qnorm(0.975) * se,
    upper_95 = fitted + qnorm(0.975) * se,
    fitted_z = as.numeric(scale(fitted))
  )
}))
fwrite(control_curves, naprt_control_curve_output_path, sep = "\t", na = "NA")

# --- motif track: miRNA seed families from miRBase mature.fa ---
mature_lines <- readLines(file.path(ref_dir, "mature.fa"))
header_index <- grep("^>", mature_lines)
headers <- mature_lines[header_index]
sequences <- mature_lines[header_index + 1L]
is_human <- startsWith(headers, ">hsa-")
mirna <- data.table(
  name = sub("^>(\\S+).*", "\\1", headers[is_human]),
  sequence = toupper(sequences[is_human])
)
mirna <- mirna[nchar(sequence) >= 18L]
mirna[, `:=`(seed8 = substr(sequence, 2L, 8L),
             is_3p = grepl("-3p$", name))]
setorder(mirna, seed8, -is_3p, name)
seed_families <- mirna[, .(
  representative = sub("^hsa-", "", name[1L]),
  n_members = .N,
  seed7 = substr(seed8[1L], 1L, 6L)
), by = "seed8"]

# --- motif track: PolyA_DB cleavage sites ---
polya_db <- fread(
  cmd = sprintf("zcat %s", shQuote(file.path(ref_dir, "PolyA_DB_hg38.txt.gz")))
)

# --- motif track: scan candidate tandem-UTR regions ---
protein_trajectories <- fread(file.path(omics_dir, "02_protein_traj.tsv"))[
  analysis_df == primary_df
]

region_summaries <- list()
site_rows <- list()
for (index in seq_len(nrow(motif_regions))) {
  region <- motif_regions[index]
  fasta_lines <- readLines(file.path(motif_seq_dir,
                                     paste0(region$fasta_file, ".fa")))
  genomic_seq <- toupper(paste(fasta_lines[-1L], collapse = ""))
  ucsc_start0 <- as.integer(sub(".*:(\\d+)-(\\d+)$", "\\1", fasta_lines[1L]))
  region_length <- nchar(genomic_seq)
  ucsc_end0 <- ucsc_start0 + region_length
  stopifnot(region_length ==
              region$locus_end - region$locus_start + 1L + 2L * flank)

  # transcript-oriented full region; locus occupies flank+1 .. flank+locus_length
  full_tx <- if (region$strand == "+") genomic_seq else revcomp(genomic_seq)
  locus_length <- region$locus_end - region$locus_start + 1L
  locus_first <- flank + 1L
  locus_last <- flank + locus_length
  transcript_seq <- substr(full_tx, locus_first, locus_last)
  # genomic 1-based coordinate of a full-region position
  to_genomic <- if (region$strand == "+") {
    function(pos) ucsc_start0 + pos
  } else {
    function(pos) ucsc_end0 - pos + 1L
  }

  # miRNA seed matches: best tier per seed family
  hits <- rbindlist(lapply(seq_len(nrow(seed_families)), function(f) {
    rc8 <- revcomp(seed_families$seed8[f])
    rc7 <- revcomp(seed_families$seed7[f])
    pos8 <- find_all(rc8, transcript_seq)
    pos7a1 <- setdiff(find_all(paste0(rc7, "A"), transcript_seq), pos8 + 1L)
    out <- list()
    if (length(pos8) > 0L) {
      next_base <- substr(paste0(transcript_seq, "N"), pos8 + 7L, pos8 + 7L)
      if (any(next_base == "A")) out[["mer8"]] <- data.table(
        class = "miRNA_8mer", label = seed_families$representative[f],
        pos_tx = pos8[which(next_base == "A")[1L]])
      if (any(next_base != "A")) out[["m8"]] <- data.table(
        class = "miRNA_7mer_m8", label = seed_families$representative[f],
        pos_tx = pos8[which(next_base != "A")[1L]])
    }
    if (length(pos7a1) > 0L) out[["a1"]] <- data.table(
      class = "miRNA_7mer_A1", label = seed_families$representative[f],
      pos_tx = pos7a1[1L])
    if (length(out) == 0L) return(NULL)
    rbindlist(out)
  }))
  pre_hits <- find_all("TGTA[ACT]ATA", transcript_seq)
  if (length(pre_hits) > 0L) {
    hits <- rbind(hits,
                  data.table(class = "PUM_PRE", label = "PUM1/2", pos_tx = pre_hits))
  }
  if (nrow(hits) > 0L) {
    hits[, `:=`(gene_symbol = region$gene_symbol,
                genomic_pos = to_genomic(flank + pos_tx))]
    site_rows[[index]] <- hits
  }

  # PAS annotation: nearest PolyA_DB cleavage site to each boundary
  gene_pa <- polya_db[Gene_Symbol == region$gene_symbol & chr == region$chrom]
  pas_rows <- rbindlist(lapply(c("proximal", "distal"), function(which_end) {
    boundary_genomic <- if (which_end == "proximal") {
      to_genomic(locus_first)
    } else {
      to_genomic(locus_last)
    }
    if (nrow(gene_pa) > 0L) {
      distances <- abs(gene_pa$pA_site - boundary_genomic)
      nearest <- which.min(distances)
      if (distances[nearest] <= 2000L) {
        return(data.table(
          end = which_end, pa_site = gene_pa$pA_site[nearest],
          pas_signal = gene_pa$PAS_Signal[nearest],
          source = "PolyA_DB", distance_bp = distances[nearest]
        ))
      }
    }
    data.table(end = which_end, pa_site = boundary_genomic,
               pas_signal = NA_character_, source = "locus_boundary",
               distance_bp = NA_integer_)
  }))
  # PAS hexamer sanity scan: 40 nt upstream of each PAS, transcript orientation
  pas_rows[, hexamer_upstream := vapply(end, function(which_end) {
    pas_full <- if (which_end == "proximal") locus_first else locus_last
    window_seq <- substr(full_tx, max(1L, pas_full - 40L), max(pas_full - 1L, 1L))
    found <- pas_hexamers[vapply(pas_hexamers, function(h) {
      length(find_all(h, window_seq)) > 0L
    }, logical(1))]
    if (length(found) == 0L) "" else paste(found, collapse = ";")
  }, character(1))]

  pair_row <- ranked[site_id == region$site_id]
  stopifnot(nrow(pair_row) == 1L)
  discovery_delta <- universe[site_id == pair_row$site_id, discovery_delta_logit]
  stopifnot(length(discovery_delta) == 1L)
  region_summaries[[index]] <- data.table(
    region[, .(gene_symbol, pair_rank, chrom, locus_start, locus_end, strand)],
    site_id = pair_row$site_id,
    protein_accession = pair_row$protein_accession,
    locus_length_tx = locus_length,
    discovery_delta_logit = discovery_delta,
    apa_direction = if (discovery_delta > 0) "lengthening" else "shortening",
    protein_endpoint_delta = protein_trajectories[
      protein_accession == pair_row$protein_accession, endpoint_delta
    ],
    n_mirna_8mer = if (nrow(hits)) sum(hits$class == "miRNA_8mer") else 0L,
    n_mirna_7mer_m8 = if (nrow(hits)) sum(hits$class == "miRNA_7mer_m8") else 0L,
    n_mirna_7mer_A1 = if (nrow(hits)) sum(hits$class == "miRNA_7mer_A1") else 0L,
    n_pum_pre = if (nrow(hits)) sum(hits$class == "PUM_PRE") else 0L,
    proximal_pa = pas_rows[end == "proximal", pa_site],
    proximal_pa_source = pas_rows[end == "proximal", source],
    proximal_pas_signal = pas_rows[end == "proximal", pas_signal],
    proximal_hexamer_upstream = pas_rows[end == "proximal", hexamer_upstream],
    distal_pa = pas_rows[end == "distal", pa_site],
    distal_pa_source = pas_rows[end == "distal", source],
    distal_pas_signal = pas_rows[end == "distal", pas_signal],
    distal_hexamer_upstream = pas_rows[end == "distal", hexamer_upstream]
  )
}

region_summary <- rbindlist(region_summaries)
motif_sites <- rbindlist(site_rows, fill = TRUE)

# --- INPP4A miRNA families ranked by brain expression ---
atlas <- as.data.table(read_parquet(mirna_atlas_path))
brain_atlas <- atlas[Biotype == "tissue" & Organ_system == "brain" & Metric == "mean"]
atlas_mirnas <- intersect(mirna$name, names(brain_atlas))
brain_long <- melt(
  brain_atlas[, c("Tissue", atlas_mirnas), with = FALSE],
  id.vars = "Tissue", variable.name = "name", value.name = "rpmm",
  variable.factor = FALSE
)
mirna_brain <- brain_long[, .(
  mature_brain_median_rpmm = median(rpmm, na.rm = TRUE),
  n_brain_tissues = sum(is.finite(rpmm))
), by = name]
family_brain <- merge(
  mirna[, .(name, seed8)], mirna_brain, by = "name", all.x = TRUE
)[, .(
  family_brain_rpmm = sum(mature_brain_median_rpmm, na.rm = TRUE),
  n_atlas_members = sum(is.finite(mature_brain_median_rpmm)),
  n_brain_tissues = max(n_brain_tissues, na.rm = TRUE),
  display_label = sub("^hsa-", "", name[which.max(mature_brain_median_rpmm)]),
  family_members = paste(sub("^hsa-", "", name), collapse = ";")
), by = seed8]
inpp4a_rank <- merge(
  motif_sites[
    gene_symbol == "INPP4A" & class %chin% c("miRNA_8mer", "miRNA_7mer_m8"),
    .(label, site_class = class, pos_tx, genomic_pos)
  ],
  seed_families[, .(label = representative, seed8, n_mirbase_members = n_members)],
  by = "label"
)
inpp4a_rank <- merge(inpp4a_rank, family_brain, by = "seed8")
setorder(inpp4a_rank, -family_brain_rpmm, label)
inpp4a_rank[, brain_rank := seq_len(.N)]
inpp4a_rank[, selected_for_label := brain_rank <= n_mirna_labels]
setcolorder(inpp4a_rank, c(
  "brain_rank", "label", "display_label", "seed8", "site_class", "pos_tx", "genomic_pos",
  "family_brain_rpmm", "n_mirbase_members", "n_atlas_members",
  "n_brain_tissues", "family_members", "selected_for_label"
))

fwrite(region_summary, motif_regions_output_path, sep = "\t", quote = FALSE)
fwrite(motif_sites, motif_sites_output_path, sep = "\t", quote = FALSE)
fwrite(inpp4a_rank, mirna_rank_output_path, sep = "\t", quote = FALSE)
