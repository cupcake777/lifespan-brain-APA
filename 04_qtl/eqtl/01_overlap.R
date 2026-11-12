# APA-QTL and eQTL overlap at gene and QTL-pair level, prenatal and postnatal.
# Fig 1: gene Venn; fig 2: pair UpSet; fig 3: lead SNP distance; fig 4: positional enrichment.

suppressMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomicFeatures)
  library(ggplot2)
  library(patchwork)
  library(txdbmaker)
  library(ggVennDiagram)
})

# --- config ---
APA_ROOT <- Sys.getenv("APA_ROOT")
DATA_DIR <- file.path(APA_ROOT, "data", "04_qtl")
REF_DIR <- file.path(APA_ROOT, "ref")
EQTL_DIR <- file.path(Sys.getenv("APA_HPC_HOME"), "project", "work", "eqtl")
FIG_DIR <- file.path(APA_ROOT, "figures", "04_qtl")
dir.create(FIG_DIR, showWarnings = FALSE)
QTL_FDR <- 0.05
FLANK <- 50
QTL_COL <- c("apaQTL" = "#D55E00", "eQTL" = "#0072B2")
VENN4_COL <- c("pre apaQTL" = "#E69F00", "pre eQTL" = "#56B4E9", "post apaQTL" = "#D55E00", "post eQTL" = "#009E73")
MAF_BREAKS <- c(-Inf, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, Inf)
MAF_LABELS <- c("<0.01", "0.01-0.05", "0.05-0.10", "0.10-0.20", "0.20-0.30", "0.30-0.40", ">=0.40")

# --- helpers ---
pt_to_mm <- function(pt) pt / .pt
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

cmh_feature_test <- function(is_case, hit, idx_by_maf) {
  counts <- rbindlist(lapply(idx_by_maf, function(idx) {
    is_case_i <- is_case[idx]
    hit_i <- hit[idx]
    data.table(a = sum(is_case_i & hit_i, na.rm = TRUE),
               b = sum(is_case_i & !hit_i, na.rm = TRUE),
               c0 = sum(!is_case_i & hit_i, na.rm = TRUE),
               d = sum(!is_case_i & !hit_i, na.rm = TRUE))
  }))
  mh_array <- array(as.vector(t(as.matrix(counts[, .(a, b, c0, d)]))),
                    dim = c(2, 2, nrow(counts)))
  mh <- mantelhaen.test(mh_array, correct = FALSE)
  mh$p.value
}

theme_sci <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(colour = "black", size = 7),
      axis.title = element_text(colour = "black", size = 8),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = grid::unit(1.6, "pt"),
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.key.size = grid::unit(0.30, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(2, "pt"),
      strip.text = element_text(face = "bold", size = 8,
                                margin = margin(b = 1.5)),
      plot.title = element_text(face = "bold", size = 9, hjust = 0),
      plot.caption = element_text(size = 7, hjust = 0),
      plot.margin = margin(2, 3, 2, 2)
    )
}

parse_variants <- function(dt) {
  dt[, c("chr_v", "pos_v") := tstrsplit(variant_id, ":", keep = 1:2)]
  dt[, `:=`(
    chr_v = paste0("chr", sub("^chr", "", chr_v)),
    pos_v = as.integer(pos_v)
  )]
  dt[]
}

parse_apa <- function(dt) {
  parts <- tstrsplit(dt$phenotype_id, "|", fixed = TRUE, fill = NA_character_)
  dt[, gene_raw := parts[[2]]]
  dt[, gene_group := if ("group_id" %in% names(dt)) group_id else gene_raw]
  dt[, ensg0 := fifelse(grepl("^ENSG", gene_group),
                        sub("\\..*", "", gene_group),
                        NA_character_)]
  dt[is.na(ensg0) & grepl("^ENSG", gene_raw), ensg0 := sub("\\..*", "", gene_raw)]
  dt[is.na(ensg0), ensg0 := sym_to_ensg[gene_group]]
  dt[is.na(ensg0), ensg0 := sym_to_ensg[gene_raw]]
  dt[]
}

parse_eqtl <- function(dt) {
  dt[, ensg0 := sub("\\..*", "", phenotype_id)]
  dt[]
}

join_gene_meta <- function(dt) {
  merge(dt,
        gene_meta[, .(ensg0, gene_chr = chr, strand, tss,
                      gene_start, gene_end)],
        by = "ensg0", all.x = TRUE)
}

add_gene_body_x <- function(dt) {
  dt <- copy(dt[!is.na(gene_start) & !is.na(gene_end) & !is.na(pos_v)])
  dt[, body_len := pmax(gene_end - gene_start, 1)]
  dt[strand == "+", gene_body_x := fcase(
    pos_v < gene_start, (pos_v - gene_start) / 1000 / FLANK,
    pos_v > gene_end,   1 + (pos_v - gene_end) / 1000 / FLANK,
    default             = (pos_v - gene_start) / body_len
  )]
  dt[strand == "-", gene_body_x := fcase(
    pos_v > gene_end,   -(pos_v - gene_end) / 1000 / FLANK,
    pos_v < gene_start, 1 + (gene_start - pos_v) / 1000 / FLANK,
    default             = (gene_end - pos_v) / body_len
  )]
  dt[gene_body_x >= -1 & gene_body_x <= 2]
}

# --- reference annotation ---
txdb <- txdbmaker::makeTxDbFromGFF(file.path(REF_DIR, "gencode.v40.annotation.gtf.gz"))
gene_gr <- genes(txdb)
id_map <- unique(fread(
  cmd = paste0("zcat ", file.path(REF_DIR, "gencode.v40.annotation.gtf.gz"), " | awk -F'\t' '$3==\"gene\" {print $9}'"),
  header = FALSE, sep = "\t", quote = "")[, .(ensg0 = sub("\\..*", "", sub('.*gene_id "([^"]+)".*', "\\1", V1)), sym   = sub('.*gene_name "([^"]+)".*', "\\1", V1))])
id_map <- id_map[!is.na(sym) & !is.na(ensg0) & !duplicated(ensg0) & !duplicated(sym)]
sym_to_ensg <- setNames(id_map$ensg0, id_map$sym)

gene_meta <- unique(data.table(
  ensg0 = sub("\\..*", "", names(gene_gr)),
  chr = as.character(seqnames(gene_gr)),
  strand = as.character(strand(gene_gr)),
  gene_start = as.integer(start(gene_gr)),
  gene_end = as.integer(end(gene_gr)),
  tss = as.integer(ifelse(as.character(strand(gene_gr)) == "+", start(gene_gr), end(gene_gr)))
))
gene_meta <- merge(gene_meta, id_map, by = "ensg0", all.x = TRUE)

prom <- trim(GRanges(seqnames(gene_gr),
                     IRanges(ifelse(as.character(strand(gene_gr)) == "+",
                                    start(gene_gr) - 2000, end(gene_gr) - 200),
                             ifelse(as.character(strand(gene_gr)) == "+",
                                    start(gene_gr) + 200, end(gene_gr) + 2000))))
utr5 <- reduce(unlist(fiveUTRsByTranscript(txdb)))
cds <- reduce(unlist(cdsBy(txdb, by = "tx")))
utr3 <- reduce(unlist(threeUTRsByTranscript(txdb)))

ccre <- fread(file.path(REF_DIR, "encode/GRCh38-cCREs.bed"), header = FALSE,
              col.names = c("chr", "start", "end", "accession", "group", "label"))
ccre[, ccre_class := fcase(
  grepl("^PLS", label), "Promoter",
  grepl("^pELS|^dELS", label), "Enhancer",
  grepl("^DNase-H3K4me3", label), "Open chromatin",
  grepl("^CTCF-only", label), "CTCF-only",
  default = "Other cCRE"
)]
ccre[, ctcf_bound := grepl("CTCF-bound", label, fixed = TRUE)]
ccre_gr <- GRanges(ccre$chr, IRanges(ccre$start + 1L, ccre$end))
mcols(ccre_gr)$ccre_class <- ccre$ccre_class
mcols(ccre_gr)$ctcf_bound <- ccre$ctcf_bound

# accumulators across stages
gene_sets <- list()
body_dt_all <- data.table()
context_summary_all <- data.table()

# --- prenatal ---
apa_pre <- fread(file.path(DATA_DIR, "prenatal.sig_QTL.txt.gz"))
eq_pre <- fread(file.path(EQTL_DIR, "prenatal.sig_QTL.txt.gz"))
apa_pre_perm <- fread(file.path(DATA_DIR, "prenatal.cis_qtl.txt.gz"))[qval <= QTL_FDR]
eq_pre_perm <- fread(file.path(DATA_DIR, "eqtl", "permutation", "prenatal.cis_qtl.txt.gz"))[qval <= QTL_FDR]

apa_pre <- join_gene_meta(parse_variants(parse_apa(apa_pre)))
eq_pre <- join_gene_meta(parse_variants(parse_eqtl(eq_pre)))
apa_pre_perm <- join_gene_meta(parse_variants(parse_apa(apa_pre_perm)))
eq_pre_perm <- join_gene_meta(parse_variants(parse_eqtl(eq_pre_perm)))

apa_pre <- apa_pre[!is.na(ensg0)]
eq_pre <- eq_pre[!is.na(ensg0)]
apa_pre_perm <- apa_pre_perm[!is.na(ensg0)]
eq_pre_perm <- eq_pre_perm[!is.na(ensg0)]

gene_sets[["pre apaQTL"]] <- unique(apa_pre$ensg0)
gene_sets[["pre eQTL"]] <- unique(eq_pre$ensg0)

pair_sets <- list()
pair_sets[["pre apaQTL"]] <- unique(paste(apa_pre$ensg0, apa_pre$variant_id, sep = "::"))
pair_sets[["pre eQTL"]] <- unique(paste(eq_pre$ensg0, eq_pre$variant_id, sep = "::"))

# --- prenatal body density ---
body_pre <- rbind(
  add_gene_body_x(apa_pre_perm)[, .(qtl = "apaQTL", gene_body_x)],
  add_gene_body_x(eq_pre_perm)[, .(qtl = "eQTL", gene_body_x)]
)
body_pre[, stage := "prenatal"]
body_dt_all <- rbind(body_dt_all, body_pre)

# --- prenatal GARFIELD enrichment ---
chmm_pre <- fread(file.path(REF_DIR, "roadmap/E082_15_coreMarks_hg38lift_mnemonics.bed.gz"),
                  header = FALSE,
                  col.names = c("chr", "start", "end", "state"))
chmm_pre[, state := fcase(
  grepl("^1_TssA", state), "Active TSS",
  grepl("^2_TssAFlnk|^3_TxFlnk|^11_BivFlnk", state), "Promoter-flanking",
  grepl("^10_TssBiv", state), "Bivalent promoter",
  grepl("^6_EnhG|^7_Enh|^12_EnhBiv", state), "Enhancer",
  grepl("^4_Tx|^5_TxWk", state), "Transcribed",
  grepl("^13_ReprPC|^14_ReprPCWk", state), "Repressed",
  grepl("^15_Quies", state), "Quiescent",
  grepl("^8_ZNF/Rpts|^9_Het", state), "Heterochromatin/repeats",
  default = "Other"
)]
chmm_pre_gr <- GRanges(paste0("chr", sub("^chr", "", chmm_pre$chr)),
                       IRanges(chmm_pre$start + 1L, chmm_pre$end))
mcols(chmm_pre_gr)$state <- chmm_pre$state

# load genotype background, then assign context (cCRE + ChromHMM + gene model)
bg_pre <- fread(file.path(DATA_DIR, "pre_input", "genotype.prenatal.pvar"),
                skip = "#CHROM", select = c(1, 2, 3, 7))
setnames(bg_pre, c("chr_v", "pos_v", "variant_id", "info"))
bg_pre[, `:=`(
  chr_v = paste0("chr", sub("^chr", "", chr_v)),
  pos_v = as.integer(pos_v)
)]
bg_pre[, maf := as.numeric(sub(".*MAF=([^;]+).*", "\\1", info))]
bg_pre[, maf_bin := as.character(cut(maf, breaks = MAF_BREAKS, labels = MAF_LABELS,
                                    include.lowest = TRUE, right = FALSE))]
bg_pre[is.na(maf_bin), maf_bin := "missing"]
bg_pre <- unique(bg_pre[, .(variant_id, chr_v, pos_v, maf, maf_bin)])
idx_by_maf_pre <- split(seq_len(nrow(bg_pre)), bg_pre$maf_bin)

bg_pre[, region := fcase(
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), prom), "Promoter/TSS",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), utr5),  "5' UTR",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), cds),   "CDS",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), utr3),  "3' UTR",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), gene_gr), "Intron",
  default = "Intergenic"
)]
bg_pre[, state := "Other"]
ov_pre <- findOverlaps(GRanges(bg_pre$chr_v, IRanges(bg_pre$pos_v, bg_pre$pos_v)),
                       chmm_pre_gr)
if (length(ov_pre)) {
  ov_dt <- data.table(qh = queryHits(ov_pre),
                      state = chmm_pre_gr$state[subjectHits(ov_pre)])
  priority <- c("Active TSS" = 9, "Promoter-flanking" = 8, Enhancer = 7,
                Transcribed = 6, "Bivalent promoter" = 5, Repressed = 4,
                Quiescent = 3, "Heterochromatin/repeats" = 2, Other = 1)
  ov_dt[, w := priority[state]]
  setorder(ov_dt, qh, -w)
  ov_dt <- ov_dt[!duplicated(qh)]
  bg_pre$state[ov_dt$qh] <- ov_dt$state
}
bg_pre[, `:=`(ccre_class = "No cCRE", ctcf_bound = "Not CTCF-bound")]
ccre_ov_pre <- findOverlaps(GRanges(bg_pre$chr_v, IRanges(bg_pre$pos_v, bg_pre$pos_v)),
                            ccre_gr)
if (length(ccre_ov_pre)) {
  ccre_dt <- data.table(
    qh = queryHits(ccre_ov_pre),
    ccre_class = ccre_gr$ccre_class[subjectHits(ccre_ov_pre)],
    ctcf_bound = ccre_gr$ctcf_bound[subjectHits(ccre_ov_pre)]
  )
  ccre_priority <- c(Promoter = 5, Enhancer = 4, "Open chromatin" = 3,
                     "CTCF-only" = 2, "Other cCRE" = 1, "No cCRE" = 0)
  ccre_dt[, w := ccre_priority[ccre_class]]
  setorder(ccre_dt, qh, -w)
  ccre_best <- ccre_dt[!duplicated(qh)]
  bg_pre$ccre_class[ccre_best$qh] <- ccre_best$ccre_class
  ccre_bound <- ccre_dt[, .(ctcf_bound = any(ctcf_bound)), by = qh]
  bg_pre$ctcf_bound[ccre_bound$qh] <- ifelse(ccre_bound$ctcf_bound,
                                             "CTCF-bound", "Not CTCF-bound")
}

# per-feature MH odds ratios
gene_levels <- c("Promoter/TSS", "5' UTR", "CDS", "3' UTR", "Intron", "Intergenic")
hmm_levels  <- c("Active TSS", "Promoter-flanking", "Enhancer", "Transcribed",
                 "Bivalent promoter", "Repressed", "Quiescent",
                 "Heterochromatin/repeats")
ccre_levels <- c("Enhancer", "Promoter", "Open chromatin", "CTCF binding site")
feature_grid_pre <- rbind(
  data.table(dim = "cCRE",       feature = ccre_levels),
  data.table(dim = "ChromHMM",   feature = hmm_levels),
  data.table(dim = "Gene model", feature = gene_levels)
)
feature_grid_pre[, feature_label := paste(dim, feature, sep = ": ")]
feature_grid_pre[, feature_y := rev(seq_len(.N))]

out_pre <- data.table()
qtl_cases <- list(
  apaQTL = unique(apa_pre_perm$variant_id),
  eQTL    = unique(eq_pre_perm$variant_id)
)
for (qi in c("apaQTL", "eQTL")) {
  bg_pre[, is_case := variant_id %in% qtl_cases[[qi]]]
  case_total <- sum(bg_pre$is_case)
  background_total <- nrow(bg_pre) - case_total
  for (i in seq_len(nrow(feature_grid_pre))) {
    di <- feature_grid_pre$dim[i]
    fi <- feature_grid_pre$feature[i]
    hit <- switch(di,
      "Gene model" = bg_pre$region == fi,
      "ChromHMM"   = bg_pre$state == fi,
      "cCRE"       = if (fi == "CTCF binding site")
                       bg_pre$ctcf_bound == "CTCF-bound"
                     else
                       bg_pre$ccre_class == fi
    )
    a  <- sum(bg_pre$is_case & hit, na.rm = TRUE)
    b  <- case_total - a
    c0 <- sum(!bg_pre$is_case & hit, na.rm = TRUE)
    d  <- background_total - c0
    a2 <- a + 0.5
    b2 <- b + 0.5
    c2 <- c0 + 0.5
    d2 <- d + 0.5
    log_or <- log((a2 * d2) / (b2 * c2))
    se <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
    base_row <- data.table(
      qtl = qi, dim = di, feature = fi,
      feature_label = feature_grid_pre$feature_label[i],
      feature_y = feature_grid_pre$feature_y[i],
      n = a, total = case_total,
      pct = a / case_total * 100,
      log2_or = log_or / log(2),
      ci_low  = (log_or - 1.96 * se) / log(2),
      ci_high = (log_or + 1.96 * se) / log(2),
      p_value = cmh_feature_test(bg_pre$is_case, hit, idx_by_maf_pre)
    )
    out_pre <- rbind(out_pre, base_row)
  }
}
out_pre[, dim := factor(dim, levels = c("cCRE", "ChromHMM", "Gene model"))]
out_pre[, qtl := factor(qtl, levels = c("apaQTL", "eQTL"))]
out_pre[, stage := "prenatal"]
context_summary_all <- rbind(context_summary_all, out_pre)

# --- postnatal ---
apa_post <- fread(file.path(DATA_DIR, "postnatal.sig_QTL.txt.gz"))
eq_post <- fread(file.path(EQTL_DIR, "postnatal.sig_QTL.txt.gz"))
apa_post_perm <- fread(file.path(DATA_DIR, "postnatal.cis_qtl.txt.gz"))[qval <= QTL_FDR]
eq_post_perm <- fread(file.path(DATA_DIR, "eqtl", "permutation", "postnatal.cis_qtl.txt.gz"))[qval <= QTL_FDR]

apa_post <- join_gene_meta(parse_variants(parse_apa(apa_post)))
eq_post <- join_gene_meta(parse_variants(parse_eqtl(eq_post)))
apa_post_perm <- join_gene_meta(parse_variants(parse_apa(apa_post_perm)))
eq_post_perm <- join_gene_meta(parse_variants(parse_eqtl(eq_post_perm)))

apa_post <- apa_post[!is.na(ensg0)]
eq_post <- eq_post[!is.na(ensg0)]
apa_post_perm <- apa_post_perm[!is.na(ensg0)]
eq_post_perm <- eq_post_perm[!is.na(ensg0)]

gene_sets[["post apaQTL"]] <- unique(apa_post$ensg0)
gene_sets[["post eQTL"]] <- unique(eq_post$ensg0)

pair_sets[["post apaQTL"]] <- unique(paste(apa_post$ensg0, apa_post$variant_id, sep = "::"))
pair_sets[["post eQTL"]] <- unique(paste(eq_post$ensg0, eq_post$variant_id, sep = "::"))

# --- postnatal body density ---
body_post <- rbind(
  add_gene_body_x(apa_post_perm)[, .(qtl = "apaQTL", gene_body_x)],
  add_gene_body_x(eq_post_perm)[, .(qtl = "eQTL", gene_body_x)]
)
body_post[, stage := "postnatal"]
body_dt_all <- rbind(body_dt_all, body_post)

# --- postnatal GARFIELD enrichment ---
chmm_post <- fread(file.path(REF_DIR, "roadmap/E073_15_coreMarks_hg38lift_mnemonics.bed.gz"),
                   header = FALSE,
                   col.names = c("chr", "start", "end", "state"))
chmm_post[, state := fcase(
  grepl("^1_TssA", state), "Active TSS",
  grepl("^2_TssAFlnk|^3_TxFlnk|^11_BivFlnk", state), "Promoter-flanking",
  grepl("^10_TssBiv", state), "Bivalent promoter",
  grepl("^6_EnhG|^7_Enh|^12_EnhBiv", state), "Enhancer",
  grepl("^4_Tx|^5_TxWk", state), "Transcribed",
  grepl("^13_ReprPC|^14_ReprPCWk", state), "Repressed",
  grepl("^15_Quies", state), "Quiescent",
  grepl("^8_ZNF/Rpts|^9_Het", state), "Heterochromatin/repeats",
  default = "Other"
)]
chmm_post_gr <- GRanges(paste0("chr", sub("^chr", "", chmm_post$chr)),
                        IRanges(chmm_post$start + 1L, chmm_post$end))
mcols(chmm_post_gr)$state <- chmm_post$state

bg_post <- fread(file.path(DATA_DIR, "pre_input", "genotype.postnatal.pvar"),
                 skip = "#CHROM", select = c(1, 2, 3, 7))
setnames(bg_post, c("chr_v", "pos_v", "variant_id", "info"))
bg_post[, `:=`(
  chr_v = paste0("chr", sub("^chr", "", chr_v)),
  pos_v = as.integer(pos_v)
)]
bg_post[, maf := as.numeric(sub(".*MAF=([^;]+).*", "\\1", info))]
bg_post[, maf_bin := as.character(cut(maf, breaks = MAF_BREAKS, labels = MAF_LABELS,
                                     include.lowest = TRUE, right = FALSE))]
bg_post[is.na(maf_bin), maf_bin := "missing"]
bg_post <- unique(bg_post[, .(variant_id, chr_v, pos_v, maf, maf_bin)])
idx_by_maf_post <- split(seq_len(nrow(bg_post)), bg_post$maf_bin)

bg_post[, region := fcase(
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), prom), "Promoter/TSS",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), utr5),  "5' UTR",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), cds),   "CDS",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), utr3),  "3' UTR",
  overlapsAny(GRanges(chr_v, IRanges(pos_v, pos_v)), gene_gr), "Intron",
  default = "Intergenic"
)]
bg_post[, state := "Other"]
ov_post <- findOverlaps(GRanges(bg_post$chr_v, IRanges(bg_post$pos_v, bg_post$pos_v)),
                        chmm_post_gr)
if (length(ov_post)) {
  ov_dt <- data.table(qh = queryHits(ov_post),
                      state = chmm_post_gr$state[subjectHits(ov_post)])
  priority <- c("Active TSS" = 9, "Promoter-flanking" = 8, Enhancer = 7,
                Transcribed = 6, "Bivalent promoter" = 5, Repressed = 4,
                Quiescent = 3, "Heterochromatin/repeats" = 2, Other = 1)
  ov_dt[, w := priority[state]]
  setorder(ov_dt, qh, -w)
  ov_dt <- ov_dt[!duplicated(qh)]
  bg_post$state[ov_dt$qh] <- ov_dt$state
}
bg_post[, `:=`(ccre_class = "No cCRE", ctcf_bound = "Not CTCF-bound")]
ccre_ov_post <- findOverlaps(GRanges(bg_post$chr_v, IRanges(bg_post$pos_v, bg_post$pos_v)),
                             ccre_gr)
if (length(ccre_ov_post)) {
  ccre_dt <- data.table(
    qh = queryHits(ccre_ov_post),
    ccre_class = ccre_gr$ccre_class[subjectHits(ccre_ov_post)],
    ctcf_bound = ccre_gr$ctcf_bound[subjectHits(ccre_ov_post)]
  )
  ccre_priority <- c(Promoter = 5, Enhancer = 4, "Open chromatin" = 3,
                     "CTCF-only" = 2, "Other cCRE" = 1, "No cCRE" = 0)
  ccre_dt[, w := ccre_priority[ccre_class]]
  setorder(ccre_dt, qh, -w)
  ccre_best <- ccre_dt[!duplicated(qh)]
  bg_post$ccre_class[ccre_best$qh] <- ccre_best$ccre_class
  ccre_bound <- ccre_dt[, .(ctcf_bound = any(ctcf_bound)), by = qh]
  bg_post$ctcf_bound[ccre_bound$qh] <- ifelse(ccre_bound$ctcf_bound,
                                              "CTCF-bound", "Not CTCF-bound")
}

feature_grid_post <- rbind(
  data.table(dim = "cCRE",       feature = ccre_levels),
  data.table(dim = "ChromHMM",   feature = hmm_levels),
  data.table(dim = "Gene model", feature = gene_levels)
)
feature_grid_post[, feature_label := paste(dim, feature, sep = ": ")]
feature_grid_post[, feature_y := rev(seq_len(.N))]

out_post <- data.table()
qtl_cases <- list(
  apaQTL = unique(apa_post_perm$variant_id),
  eQTL    = unique(eq_post_perm$variant_id)
)
for (qi in c("apaQTL", "eQTL")) {
  bg_post[, is_case := variant_id %in% qtl_cases[[qi]]]
  case_total <- sum(bg_post$is_case)
  background_total <- nrow(bg_post) - case_total
  for (i in seq_len(nrow(feature_grid_post))) {
    di <- feature_grid_post$dim[i]
    fi <- feature_grid_post$feature[i]
    hit <- switch(di,
      "Gene model" = bg_post$region == fi,
      "ChromHMM"   = bg_post$state == fi,
      "cCRE"       = if (fi == "CTCF binding site")
                       bg_post$ctcf_bound == "CTCF-bound"
                     else
                       bg_post$ccre_class == fi
    )
    a  <- sum(bg_post$is_case & hit, na.rm = TRUE)
    b  <- case_total - a
    c0 <- sum(!bg_post$is_case & hit, na.rm = TRUE)
    d  <- background_total - c0
    a2 <- a + 0.5
    b2 <- b + 0.5
    c2 <- c0 + 0.5
    d2 <- d + 0.5
    log_or <- log((a2 * d2) / (b2 * c2))
    se <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
    base_row <- data.table(
      qtl = qi, dim = di, feature = fi,
      feature_label = feature_grid_post$feature_label[i],
      feature_y = feature_grid_post$feature_y[i],
      n = a, total = case_total,
      pct = a / case_total * 100,
      log2_or = log_or / log(2),
      ci_low  = (log_or - 1.96 * se) / log(2),
      ci_high = (log_or + 1.96 * se) / log(2),
      p_value = cmh_feature_test(bg_post$is_case, hit, idx_by_maf_post)
    )
    out_post <- rbind(out_post, base_row)
  }
}
out_post[, dim := factor(dim, levels = c("cCRE", "ChromHMM", "Gene model"))]
out_post[, qtl := factor(qtl, levels = c("apaQTL", "eQTL"))]
out_post[, stage := "postnatal"]
context_summary_all <- rbind(context_summary_all, out_post)

# --- fig1: gene-level venn ---

gene_venn_sets <- gene_sets[names(VENN4_COL)]

p_gene_venn <- ggVennDiagram(
  gene_venn_sets,
  set_color = unname(VENN4_COL),
  set_size = 3.8,
  label = "count",
  label_alpha = 0,
  label_size = 3.2,
  label_color = "black",
  edge_size = 0.35
) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.10)) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(title = "Sig-QTL genes across stage and modality") +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
    plot.margin = margin(4, 6, 4, 6)
  )

ggsave(file.path(FIG_DIR, "fig1_gene_overlap_combined.png"), p_gene_venn,
       width = 5.6, height = 5.0, dpi = 600, bg = "white")

# --- fig2: QTL pair-level 4-set Upset ---
pair_sets <- pair_sets[names(VENN4_COL)]
set_names <- names(pair_sets)
all_pairs <- unique(unlist(pair_sets))
member_pairs <- data.table(pair = all_pairs)
for (nm in set_names) member_pairs[, (nm) := pair %in% pair_sets[[nm]]]
member_pairs[, pattern := apply(.SD, 1, function(x) paste(set_names[as.logical(x)], collapse = "&")),
             .SDcols = set_names]
counts_pairs <- member_pairs[, .N, by = pattern][pattern != ""]
setorder(counts_pairs, -N, pattern)
counts_pairs[, pattern_id := seq_len(.N)]
counts_pairs[, pattern_f := factor(pattern_id, levels = pattern_id)]

matrix_dt_pairs <- rbindlist(lapply(seq_len(nrow(counts_pairs)), function(i) {
  present_sets <- strsplit(counts_pairs$pattern[i], "&", fixed = TRUE)[[1]]
  data.table(pattern_id = counts_pairs$pattern_id[i],
             pattern_f = counts_pairs$pattern_f[i],
             set = set_names,
             present = set_names %in% present_sets)
}))
set_dt_pairs <- data.table(
  set = set_names,
  y = seq_along(set_names),
  set_label = sprintf("%s (n=%s)", set_names, comma(lengths(pair_sets)))
)
matrix_dt_pairs <- merge(matrix_dt_pairs, set_dt_pairs, by = "set")
seg_dt_pairs <- matrix_dt_pairs[present == TRUE, .(ymin = min(y), ymax = max(y)),
                                by = .(pattern_id, pattern_f)]

p_bar_pairs <- ggplot(counts_pairs, aes(pattern_f, N)) +
  geom_col(fill = "#4D4D4D", width = 0.72) +
  geom_text(aes(label = comma(N)), vjust = -0.25, size = pt_to_mm(7)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(title = "Sig-QTL pairs across stage and modality",
       x = NULL, y = "Intersection QTL pairs") +
  theme_sci() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())

p_matrix_pairs <- ggplot(matrix_dt_pairs, aes(pattern_f, y)) +
  geom_segment(data = seg_dt_pairs,
               aes(x = pattern_f, xend = pattern_f, y = ymin, yend = ymax),
               inherit.aes = FALSE, colour = "grey35", linewidth = 0.7) +
  geom_point(aes(fill = present), shape = 21, size = 1.8, colour = "grey35", stroke = 0.25) +
  scale_fill_manual(values = c("TRUE" = "#1F1F1F", "FALSE" = "#D9D9D9"), guide = "none") +
  scale_y_continuous(breaks = set_dt_pairs$y, labels = set_dt_pairs$set_label,
                     expand = expansion(add = 0.45)) +
  labs(x = "Set intersection", y = NULL) +
  theme_sci() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())

fig2 <- p_bar_pairs / p_matrix_pairs + plot_layout(heights = c(2.5, 1.2))
ggsave(file.path(FIG_DIR, "fig2_qtl_pair_overlap_upset.png"), fig2,
       width = 7.09, height = 4.1, dpi = 600, bg = "white")

# --- fig3: lead SNP body density, faceted by stage ---
fig3 <- ggplot(body_dt_all, aes(gene_body_x, colour = qtl)) +
  geom_density(linewidth = 0.95, adjust = 1.05) +
  geom_vline(xintercept = c(0, 1), linetype = "dotted", colour = "grey35",
             linewidth = 0.55) +
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0,
           colour = "grey45", linewidth = 0.6) +
  annotate("text", x = 0.5, y = 0, label = "gene body", vjust = -0.8,
           size = pt_to_mm(7), colour = "grey25") +
  scale_colour_manual(values = QTL_COL) +
  scale_x_continuous(
    breaks = c(-1, 0, 1, 2),
    labels = c(paste0("-", FLANK, " kb"), "5'-end", "3'-end", paste0("+", FLANK, " kb"))
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Lead SNP positional distribution",
       x = NULL, y = "Variant density") +
  facet_wrap(~stage, nrow = 1) +
  theme_sci() +
  theme(legend.position = "top",
        panel.spacing.x = grid::unit(0.65, "cm"),
        panel.border = element_blank())
ggsave(file.path(FIG_DIR, "fig3_lead_snp_distance.png"), fig3,
       width = 7.09, height = 2.9, dpi = 600, bg = "white")

# --- fig4: GARFIELD-style enrichment, facet_grid(dim ~ stage) ---
stats_dt <- copy(context_summary_all)
stats_dt[, qtl_chr := as.character(qtl)]
stats_dt[, y_plot := feature_y + fifelse(qtl_chr == "apaQTL", 0.16, -0.16)]
stats_dt[, p_adj := p.adjust(p_value, method = "BH"), by = .(stage, qtl_chr)]
stats_dt[, sig_label := fifelse(p_adj < 0.001, "***",
                         fifelse(p_adj < 0.01, "**",
                         fifelse(p_adj < 0.05, "*", "")))]
axis_dt <- unique(stats_dt[, .(feature_y, feature)])
setorder(axis_dt, feature_y)
axis_labels <- setNames(axis_dt$feature, axis_dt$feature_y)
xlim_left <- range(c(stats_dt$ci_low, stats_dt$ci_high), na.rm = TRUE)
if (any(!is.finite(xlim_left))) xlim_left <- c(-1, 1)
xlim_left <- c(floor(min(-0.5, xlim_left[1])), ceiling(max(0.5, xlim_left[2])))
stripe_dt <- unique(stats_dt[, .(stage, dim, feature_y)])
stripe_dt <- stripe_dt[feature_y %% 2 == 0]

p_left <- ggplot(stats_dt, aes(y = y_plot)) +
  geom_rect(data = stripe_dt,
            aes(ymin = feature_y - 0.5, ymax = feature_y + 0.5,
                xmin = -Inf, xmax = Inf),
            inherit.aes = FALSE, fill = "grey93", colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.58) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = y_plot, colour = qtl_chr),
               linewidth = 0.8, alpha = 0.92) +
  geom_point(aes(x = log2_or, colour = qtl_chr), size = 1.4) +
  scale_colour_manual(values = QTL_COL, breaks = c("apaQTL", "eQTL"),
                      guide = guide_legend(override.aes = list(linetype = 0, linewidth = 0,
                                                               shape = 16, size = 1.6))) +
  scale_y_continuous(breaks = axis_dt$feature_y, labels = axis_labels,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(limits = xlim_left) +
  labs(x = "log2 odds ratio", y = NULL) +
  facet_grid(dim ~ stage, scales = "free_y", space = "free_y") +
  theme_sci() +
  theme(legend.position = "top",
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
        legend.key.height = grid::unit(0.30, "cm"),
        legend.key.width  = grid::unit(0.30, "cm"),
        legend.text = element_text(size = 6),
        strip.placement = "outside",
        panel.spacing.y = grid::unit(0.10, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35))

pct_dt <- stats_dt[, .(stage, dim, feature_y, group = qtl_chr, pct, sig_label,
                       total, n,
                       y_plot = feature_y + fifelse(qtl_chr == "apaQTL", 0.16, -0.16))]
pct_dt[, group := factor(group, levels = c("apaQTL", "eQTL"))]
pct_ci <- t(vapply(seq_len(nrow(pct_dt)), function(i) {
  bt <- binom.test(pct_dt$n[i], pct_dt$total[i])
  c(low = bt$conf.int[1], high = bt$conf.int[2])
}, numeric(2)))
pct_dt[, `:=`(pct_low = pct_ci[, "low"], pct_high = pct_ci[, "high"])]
pct_dt[, `:=`(prop = pct / 100, prop_low = pct_low, prop_high = pct_high)]
max_pct <- max(pct_dt$pct_high, na.rm = TRUE)
if (!is.finite(max_pct)) max_pct <- 1
pct_limit <- max(0.60, ceiling(max_pct * 10) / 10)
pct_sig_dt <- pct_dt[sig_label != ""]
if (nrow(pct_sig_dt)) {
  pct_sig_dt[, sig_x := prop_high + pct_limit * 0.035]
  pct_limit <- max(pct_limit, ceiling(max(pct_sig_dt$sig_x, na.rm = TRUE) * 10) / 10)
}
pct_breaks <- c(0, 0.1, 0.5)

p_right <- ggplot() +
  geom_rect(data = pct_dt,
            aes(xmin = 0, xmax = prop, ymin = y_plot - 0.12, ymax = y_plot + 0.12,
                fill = group),
            colour = "grey30", linewidth = 0.2) +
  geom_segment(data = pct_dt,
               aes(x = prop_low, xend = prop_high, y = y_plot, yend = y_plot,
                   colour = group),
               linewidth = 0.55) +
  geom_text(data = pct_sig_dt,
            aes(x = sig_x, y = y_plot, label = sig_label),
            inherit.aes = FALSE, size = pt_to_mm(6.8), colour = "black", hjust = 0) +
  scale_fill_manual(values = QTL_COL, guide = "none") +
  scale_colour_manual(values = QTL_COL, guide = "none") +
  scale_y_continuous(breaks = axis_dt$feature_y, labels = NULL,
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(limits = c(0, pct_limit),
                     breaks = pct_breaks,
                     labels = function(x)
                       ifelse(x == 0, "0", sub("^0\\.", "0.", format(x, trim = TRUE))),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = "Proportion\nof variants", y = NULL) +
  facet_grid(dim ~ stage, scales = "free_y", space = "free_y") +
  theme_sci() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text.y = element_blank(),
        panel.spacing.x = grid::unit(0.14, "cm"),
        panel.spacing.y = grid::unit(0.10, "cm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35))

fig4 <- p_left + p_right + plot_layout(widths = c(3.05, 1.05), guides = "collect") &
  theme(legend.position = "top")
ggsave(file.path(FIG_DIR, "fig4_positional_enrichment.png"), fig4,
       width = 7.09, height = 5.6, dpi = 600, bg = "white")
