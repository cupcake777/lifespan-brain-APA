resolve_existing_file <- function(paths) {
  paths <- Filter(function(path) !is.null(path) && nzchar(path), paths)
  for (path in paths) {
    if (file.exists(path)) return(path)
  }
  NULL
}

resolve_recomb_prefix <- function(prefixes, chrom = "1") {
  prefixes <- Filter(function(path) !is.null(path) && nzchar(path), prefixes)
  for (prefix in prefixes) {
    txt_path <- paste0(prefix, "_recombination_map_hapmap_format_hg38_chr_", chrom, ".txt")
    bed_path <- paste0(prefix, "_recombination_map_hg38_chr_", chrom, ".bed")
    if (file.exists(txt_path) || file.exists(bed_path)) return(prefix)
  }
  NULL
}

resolve_synthetic_plot_inputs <- function(cfg_global, fixture_dir = "examples") {
  list(
    gtf_path = resolve_existing_file(c(
      file.path(fixture_dir, "smoke_hg38_chr1.gtf"),
      cfg_global$gene_anno
    )),
    recomb_path = resolve_recomb_prefix(c(
      file.path(fixture_dir, "CHB"),
      cfg_global$recom
    ))
  )
}

make_synthetic_locus <- function() {
  set.seed(20260327)

  n_bg <- 80
  lead_snp <- "rs10890255"
  lead_pos <- 43566000
  bg_pos <- sort(sample(seq(43120000, 43580000, by = 500), n_bg))
  signal_pos <- c(43556000, 43561000, lead_pos, 43570000, 43576000, 43582000)
  signal_rs <- c("rs2819340", "rs3791138", lead_snp, "rs3791137", "rs2004899", "rs2819341")
  all_pos <- c(bg_pos, signal_pos)
  all_rsid <- c(paste0("rsBG", seq_len(n_bg)), signal_rs)

  distance_kb <- abs(all_pos - lead_pos) / 1000
  r2_to_lead <- pmin(1, pmax(0.02, exp(-distance_kb / 92) + rnorm(length(all_pos), sd = 0.035)))
  r2_to_lead[all_rsid == lead_snp] <- 1
  r2_to_lead[all_rsid %in% c("rs2819340", "rs3791138", "rs3791137", "rs2004899")] <- c(0.82, 0.91, 0.78, 0.62)

  shared_signal <- 1.25 + 8.5 * r2_to_lead + rnorm(length(all_pos), sd = 0.35)
  gwas_logp <- pmax(0.8, shared_signal + rnorm(length(all_pos), sd = 0.55))
  qtl_logp <- pmax(0.8, shared_signal * 0.93 + rnorm(length(all_pos), sd = 0.45))
  gwas_logp[all_rsid == lead_snp] <- 9.2
  qtl_logp[all_rsid == lead_snp] <- 8.9

  df <- data.frame(
    rsid = all_rsid,
    CHR.qtl = "chr1",
    POS.qtl = all_pos,
    P.qtl = 10^(-qtl_logp),
    P.gwas = 10^(-gwas_logp),
    stringsAsFactors = FALSE
  )

  credible_set <- data.frame(
    snp = c(lead_snp, "rs3791138", "rs2819340", "rs3791137", "rs2004899"),
    SNP.PP.H4 = c(0.36, 0.27, 0.18, 0.11, 0.06),
    stringsAsFactors = FALSE
  )

  ld_pairs <- do.call(
    rbind,
    lapply(seq_along(all_rsid), function(i) {
      data.frame(
        SNP_A = lead_snp,
        SNP_B = all_rsid[i],
        R = sqrt(r2_to_lead[i]),
        stringsAsFactors = FALSE
      )
    })
  )
  ld_pairs <- ld_pairs[ld_pairs$SNP_A != ld_pairs$SNP_B, , drop = FALSE]
  ld_pairs <- rbind(
    ld_pairs,
    data.frame(SNP_A = ld_pairs$SNP_B, SNP_B = ld_pairs$SNP_A, R = ld_pairs$R, stringsAsFactors = FALSE)
  )

  list(
    merged_data = df,
    ld_df = ld_pairs,
    credible_set = credible_set,
    lead_snp = lead_snp,
    phenotype_info = "ENST00000654683.1|CCDC30|chr1:42482663-42484158|+",
    qtl_type = "postnatal"
  )
}
