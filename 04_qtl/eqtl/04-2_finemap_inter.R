# SuSiE (summary-statistics RSS) fine-mapping of the G x Stage z-scores per chromosome, with
# PIP aggregation. Usage: Rscript 04-2_finemap_inter.R <chr> [max_events|Inf] [out_subdir] [restrict_ids]
suppressMessages({
  library(data.table)
  library(susieR)
  library(parallel)
})

# --- paths and arguments ---
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
max_events <- if (length(args) >= 2 && args[2] != "Inf") as.integer(args[2]) else Inf
out_sub <- if (length(args) >= 3) args[3] else "finemap"
restrict_file <- if (length(args) >= 4) args[4] else ""

APA_ROOT <- Sys.getenv("APA_ROOT")
stopifnot(nzchar(APA_ROOT))
base_dir <- file.path(APA_ROOT, "data", "04_qtl", "interaction_sldsc")
z_file <- file.path(base_dir, "z", paste0("chr", chr, ".tsv.gz"))
gt_prefix <- file.path(base_dir, "genotype", paste0("chr", chr))
out_dir <- file.path(base_dir, out_sub)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# --- inputs ---
zt <- fread(z_file)
bim <- fread(paste0(gt_prefix, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
fam <- fread(paste0(gt_prefix, ".fam"), header = FALSE)
n <- nrow(fam)
bpr <- ceiling(n / 4)
bed_raw <- readBin(paste0(gt_prefix, ".bed"), "raw", file.info(paste0(gt_prefix, ".bed"))$size)
stopifnot(length(bed_raw) == 3 + nrow(bim) * bpr,
          identical(as.integer(bed_raw[1:3]), c(108L, 27L, 1L)))  # PLINK 1 bed size and magic

ev_ids <- unique(zt$phenotype_id)
if (nzchar(restrict_file)) ev_ids <- intersect(fread(restrict_file, header = FALSE)$V1, ev_ids)
if (length(ev_ids) > max_events) ev_ids <- ev_ids[seq_len(max_events)]
zt <- zt[phenotype_id %in% ev_ids]

# --- genotype decoding (PLINK 1 bed: 2 bits per sample; dosage = count of bim A2) ---
lut <- t(vapply(0:255, function(b) c(0, NA, 1, 2)[c(b %/% 1, b %/% 4, b %/% 16, b %/% 64) %% 4 + 1], numeric(4)))
lut_flat <- as.numeric(t(lut))
decode_block <- function(lo, hi) {
  bytes <- bed_raw[(4 + (lo - 1) * bpr):(3 + hi * bpr)]
  bi <- as.integer(bytes)
  g <- lut_flat[rep(bi, each = 4L) * 4L + rep(1:4, times = length(bi))]
  m <- matrix(g, nrow = hi - lo + 1, ncol = bpr * 4, byrow = TRUE)
  t(m[, seq_len(n), drop = FALSE])
}

# --- allele orientation: pair z is relative to ALT; the dosage counts bim A2 ---
zu <- unique(zt[, .(variant_id, af)])
sp <- tstrsplit(zu$variant_id, ":", fixed = TRUE)
zu[, `:=`(ref = sp[[3]], alt = sp[[4]])]
midx_u <- match(zu$variant_id, bim$id)
chk <- which(!is.na(midx_u))
chk <- chk[seq(1, length(chk), length.out = min(5000, length(chk)))]
frac_a1_alt <- mean(bim$a1[midx_u[chk]] == zu$alt[chk])
frac_a2_alt <- mean(bim$a2[midx_u[chk]] == zu$alt[chk])
set.seed(1)
chk20 <- sample(which(!is.na(midx_u)), min(20, sum(!is.na(midx_u))))
G20 <- decode_block(min(midx_u[chk20]), max(midx_u[chk20]))[, midx_u[chk20] - min(midx_u[chk20]) + 1, drop = FALSE]
a2_freq <- colMeans(G20, na.rm = TRUE) / 2
err_if_a2_alt <- mean(abs(a2_freq - zu$af[chk20]))
err_if_a1_alt <- mean(abs(1 - a2_freq - zu$af[chk20]))
flip <- frac_a1_alt > frac_a2_alt
stopifnot(abs(err_if_a2_alt - err_if_a1_alt) > 0.05)
if ((err_if_a1_alt < err_if_a2_alt) != flip) stop("string-based and af-based orientation verdicts disagree")

# --- per-event SuSiE (L = 10) ---
event_rows <- split(seq_len(nrow(zt)), factor(zt$phenotype_id, levels = ev_ids))
results <- mclapply(event_rows, function(rows) {
  pid <- zt$phenotype_id[rows[1]]
  t0 <- proc.time()[3]
  ev <- zt[rows]
  midx <- match(ev$variant_id, bim$id)
  ok <- !is.na(midx)
  ev <- ev[ok]; midx <- midx[ok]
  n_tested <- length(midx)
  summ <- data.table(phenotype_id = pid, n_variants = n_tested, max_pip = NA_real_, top_variant = NA_character_,
                     n_cs_95 = NA_integer_, converged = NA, niter = NA_integer_, seconds = NA_real_, status = "ok")
  if (n_tested < 50) {
    summ[, `:=`(status = "skipped_lt50", seconds = proc.time()[3] - t0)]
    return(list(summary = summ))
  }
  lo <- min(midx); hi <- max(midx)
  G <- decode_block(lo, hi)[, midx - lo + 1, drop = FALSE]
  if (flip) G <- 2 - G
  cm <- colMeans(G, na.rm = TRUE)
  nobs <- colSums(!is.na(G))
  csd <- sqrt(colSums((sweep(G, 2, cm, `-`))^2, na.rm = TRUE) / (nobs - 1))
  keep <- is.finite(csd) & csd > 0
  G <- G[, keep, drop = FALSE]; cm <- cm[keep]; csd <- csd[keep]
  z_ev <- ev$z[keep]; midx_k <- midx[keep]
  Gs <- sweep(sweep(G, 2, cm, `-`), 2, csd, `/`)
  Gs[is.na(Gs)] <- 0
  R <- crossprod(Gs) / (nrow(Gs) - 1)
  fit <- tryCatch(susie_rss(z = z_ev, R = R, n = n, L = 10), error = function(e) e)
  if (inherits(fit, "error")) {
    fit <- tryCatch(susie_rss(z = z_ev, R = R + 1e-6 * diag(ncol(R)), n = n, L = 10), error = function(e) e)
  }
  secs <- proc.time()[3] - t0
  if (inherits(fit, "error")) {
    summ[, `:=`(status = paste0("failed:", conditionMessage(fit)), seconds = secs)]
    return(list(summary = summ))
  }
  w <- which.max(fit$pip)
  summ[, `:=`(max_pip = fit$pip[w], top_variant = ev$variant_id[keep][w], n_cs_95 = length(fit$sets$cs),
              converged = fit$converged, niter = fit$niter, seconds = secs)]
  list(idx = midx_k, pip = fit$pip, summary = summ)
}, mc.cores = cores, mc.preschedule = FALSE)

# --- max PIP per variant across events, and event table ---
maxpip <- numeric(nrow(bim))
tested_union <- integer(0)
for (res in results) {
  if (!is.null(res$pip)) {
    maxpip[res$idx] <- pmax(maxpip[res$idx], res$pip)
    tested_union <- c(tested_union, res$idx)
  }
}
tested_union <- sort(unique(tested_union))
events_dt <- rbindlist(lapply(results, `[[`, "summary"))
fwrite(data.table(variant_id = bim$id[tested_union], max_pip = maxpip[tested_union]),
       file.path(out_dir, paste0("chr", chr, ".maxpip.tsv.gz")), sep = "\t")
fwrite(events_dt, file.path(out_dir, paste0("chr", chr, ".events.tsv.gz")), sep = "\t")
