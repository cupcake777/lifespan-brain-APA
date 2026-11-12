# Section 2.2 supplement on the matched APA/RNA/protein cohort: (a) APA-protein coupling beyond RNA and age,
# (c) within-age APA-RNA coupling, (d) trajectory-level decoupling; cis versus trans-gene nulls (100 sets).
suppressPackageStartupMessages({ library(data.table); library(splines) })
ROOT <- Sys.getenv("APA_ROOT")
D3 <- file.path(ROOT, "data/03_omics")
OD <- D3; NN <- 100L

apa <- readRDS(file.path(D3, "apa_matched_raw.rds"))
rna <- readRDS(file.path(D3, "rna_matched_raw.rds"))
pro <- readRDS(file.path(D3, "protein_matched_raw.rds"))
man <- fread(file.path(D3, "matched_sample_manifest_raw.tsv"))
stopifnot(identical(colnames(apa), man$sample_id), identical(colnames(rna), man$sample_id), identical(colnames(pro), man$sample_id))

t3  <- fread(file.path(D3, "03_apa_incremental_prot_tests_results.tsv"))
tri <- fread(file.path(D3, "09_global_triad_correlations.tsv"))
u <- merge(t3, tri[, .(site_id, protein_accession, cluster_id, apa_dynamic)], by = c("site_id", "protein_accession"))
u <- u[!is.na(p_rna_adjusted)]
rna_base <- sub("\\.\\d+$", "", rownames(rna))
u[, `:=`(rna_idx = match(gencode_gene_base, rna_base), apa_idx = match(site_id, rownames(apa)), pro_idx = match(protein_accession, rownames(pro)))]
stopifnot(!anyNA(u$rna_idx), !anyNA(u$apa_idx), !anyNA(u$pro_idx))

n <- ncol(apa); age_model <- log10(man$age_years + 1); B <- cbind(1, ns(age_model, df = 2)); dfres <- n - 5
pi0 <- function(p, lam = seq(0.5, 0.9, 0.05)) min(1, median(sapply(lam, function(l) mean(p > l) / (1 - l))))

apa_t <- function(a_idx, p_idx, r_idx, Amat = apa) {
  out <- numeric(length(a_idx))
  for (i in seq_along(a_idx)) {
    X <- cbind(B, rna[r_idx[i], ]); qx <- qr(X)
    ya <- qr.resid(qx, Amat[a_idx[i], ]); yp <- qr.resid(qx, pro[p_idx[i], ])
    b <- sum(ya * yp) / sum(ya^2); res <- yp - b * ya
    out[i] <- b / sqrt(sum(res^2) / dfres / sum(ya^2))
  }
  out
}
draw_trans <- function(gene, pool_gene) {
  j <- sample(length(pool_gene), length(gene), TRUE); bad <- pool_gene[j] == gene
  while (any(bad)) { j[bad] <- sample(length(pool_gene), sum(bad), TRUE); bad <- pool_gene[j] == gene }
  j
}

# (a) APA -> protein beyond RNA and age
obs_t <- apa_t(u$apa_idx, u$pro_idx, u$rna_idx)
summ_a <- function(t, cl, g) {
  p <- 2 * pt(-abs(t), dfres); q <- p.adjust(p, "BH")
  d <- data.table(t = t, g = g, cl = cl)[, .(mt = mean(t), c2 = mean(cl == 2)), by = g]
  data.table(frac_p05 = mean(p < 0.05), pi1 = 1 - pi0(p), n_q05 = sum(q < 0.05), frac_neg = mean(t < 0), mean_t = mean(t),
             pi1_C1 = 1 - pi0(p[cl == 1]), pi1_C2 = 1 - pi0(p[cl == 2]), neg_C1 = mean(t[cl == 1] < 0), neg_C2 = mean(t[cl == 2] < 0),
             g_frac_neg = mean(d$mt < 0), g_mean_t = mean(d$mt),
             g_frac_neg_C2major = mean(d$mt[d$c2 > 0.5] < 0), g_frac_neg_C1major = mean(d$mt[d$c2 < 0.5] < 0))
}
gp <- unique(u[, .(gencode_gene_base, pro_idx, rna_idx)])
set.seed(7)
NULLT <- vector("list", NN)
null_a <- rbindlist(lapply(seq_len(NN), function(r) {
  j <- draw_trans(u$gencode_gene_base, gp$gencode_gene_base)
  tt <- apa_t(u$apa_idx, gp$pro_idx[j], gp$rna_idx[j]); NULLT[[r]] <<- tt
  summ_a(tt, u$cluster_id, u$gencode_gene_base)
}))
saveRDS(do.call(cbind, NULLT), file.path(OD, "s22a_null_t.rds"))
obs_a <- summ_a(obs_t, u$cluster_id, u$gencode_gene_base)
lower_tail <- c("mean_t", "g_mean_t")
emp <- function(obs, nul) sapply(names(obs), function(k) {
  x <- nul[[k]]; o <- obs[[k]]
  if (k %in% lower_tail) (1 + sum(x <= o)) / (1 + length(x)) else (1 + sum(x >= o)) / (1 + length(x)) })
tab_a <- data.table(metric = names(obs_a), observed = unlist(obs_a), null_mean = sapply(null_a, mean),
                    null_q95 = sapply(null_a, quantile, 0.95), null_q05 = sapply(null_a, quantile, 0.05), emp_p = emp(obs_a, null_a))
fwrite(tab_a, file.path(OD, "s22a_apa_prot_cis_vs_trans.tsv"), sep = "\t")
fwrite(null_a, file.path(OD, "s22a_null_reps.tsv"), sep = "\t")
fwrite(u[, .(site_id, protein_accession, gene_symbol, gencode_gene_base, cluster_id, t_apa_M1 = obs_t)], file.path(OD, "s22a_pair_t.tsv"), sep = "\t")

# (c) within-age APA-RNA coupling
s <- unique(u[, .(site_id, gencode_gene_base, cluster_id, apa_idx, rna_idx)])
qB <- qr(B); RA <- t(qr.resid(qB, t(apa[s$apa_idx, ]))); RR <- t(qr.resid(qB, t(rna)))
rowcor <- function(X, Y) { X <- X - rowMeans(X); Y <- Y - rowMeans(Y); rowSums(X * Y) / sqrt(rowSums(X^2) * rowSums(Y^2)) }
rp <- function(r) 2 * pt(-abs(r * sqrt((n - 5) / (1 - r^2))), n - 5)
summ_c <- function(r, cl) { p <- rp(r); q <- p.adjust(p, "BH")
  data.table(pi1 = 1 - pi0(p), frac_q05 = mean(q < 0.05), frac_pos_q05 = mean(q < 0.05 & r > 0), frac_neg_q05 = mean(q < 0.05 & r < 0),
             median_r = median(r), frac_absr_lt0.1 = mean(abs(r) < 0.1), pi1_C1 = 1 - pi0(p[cl == 1]), pi1_C2 = 1 - pi0(p[cl == 2])) }
r_obs <- rowcor(RA, RR[s$rna_idx, ])
ug <- unique(s[, .(gencode_gene_base, rna_idx)]); set.seed(11)
null_c <- rbindlist(lapply(seq_len(NN), function(i) {
  j <- draw_trans(s$gencode_gene_base, ug$gencode_gene_base); summ_c(rowcor(RA, RR[ug$rna_idx[j], ]), s$cluster_id) }))
obs_c <- summ_c(r_obs, s$cluster_id)
fwrite(null_c, file.path(OD, "s22c_null_reps.tsv"), sep = "\t")
tab_c <- data.table(metric = names(obs_c), observed = unlist(obs_c), null_mean = sapply(null_c, mean),
                    null_q95 = sapply(null_c, quantile, 0.95), null_q05 = sapply(null_c, quantile, 0.05),
                    emp_p_upper = sapply(names(obs_c), function(k) (1 + sum(null_c[[k]] >= obs_c[[k]])) / (1 + NN)),
                    emp_p_lower = sapply(names(obs_c), function(k) (1 + sum(null_c[[k]] <= obs_c[[k]])) / (1 + NN)))
fwrite(tab_c, file.path(OD, "s22c_apa_rna_cis_vs_trans.tsv"), sep = "\t")
fwrite(cbind(s[, .(site_id, gencode_gene_base, cluster_id)], r = r_obs, p = rp(r_obs), q = p.adjust(rp(r_obs), "BH")),
       file.path(OD, "s22c_apa_rna_site_coupling.tsv"), sep = "\t")

# (d) trajectory-level decoupling (C1 versus C2 mRNA trajectories)
at <- fread(file.path(D3, "02_apa_trajectories_results.tsv"))[analysis_role == "primary" & apa_dynamic == TRUE]
atu <- unique(at[, .(site_id, gencode_gene_base, cluster_id, apa_delta = endpoint_delta)], by = "site_id")
genes <- unique(atu$gencode_gene_base); gi <- match(genes, rna_base); ok <- !is.na(gi); genes <- genes[ok]; gi <- gi[ok]
Xa <- cbind(1, ns(age_model, df = 2)); qa <- qr(Xa); Y <- t(rna[gi, ])
fitv <- qr.fitted(qa, Y); rs <- qr.resid(qa, Y); rsd <- sqrt(colSums(rs^2) / (n - 4))
F <- ((colSums((Y - rep(colMeans(Y), each = n))^2) - colSums(rs^2)) / 3) / rsd^2; pF <- pf(F, 3, n - 4, lower.tail = FALSE)
lo <- which.min(age_model); hi <- which.max(age_model)
rtg <- data.table(gencode_gene_base = genes, rna_delta = fitv[hi, ] - fitv[lo, ], rna_std = (apply(fitv, 2, max) - apply(fitv, 2, min)) / rsd, rna_q = p.adjust(pF, "BH"))
m <- merge(atu, rtg, by = "gencode_gene_base")
g2 <- m[, .(cl = round(mean(cluster_id)), apa_delta = median(apa_delta), rna_delta = rna_delta[1], rna_std = rna_std[1], rna_q = rna_q[1]), by = gencode_gene_base]
tab_d <- g2[, .(n_genes = .N, frac_rna_up = mean(rna_delta > 0), median_rna_std_amp = median(rna_std),
                frac_rna_age_q05 = mean(rna_q < 0.05), rho_apa_rna_delta = cor(apa_delta, rna_delta, method = "spearman")), by = cl][order(cl)]
tab_d[, `:=`(p_rna_up_C1_vs_C2 = prop.test(g2[, table(cl, rna_delta > 0)])$p.value,
             p_rna_amp_C1_vs_C2 = wilcox.test(rna_std ~ cl, data = g2)$p.value)]
fwrite(tab_d, file.path(OD, "s22d_traj_decoupling.tsv"), sep = "\t")
fwrite(g2, file.path(OD, "s22d_gene_traj.tsv"), sep = "\t")
