#!/usr/bin/env Rscript
# GO/KEGG term clusters, GO slim and BRITE mapping, and broad gene-set enrichment
# Module 02 inputs and outputs live in data/02_traj; cluster-only references stay under APA_HPC_HOME
suppressPackageStartupMessages({
  library(data.table)
  library(GO.db)
  library(AnnotationDbi)
  library(jsonlite)
  library(readxl)
})

# --- config ---

out <- file.path(Sys.getenv("APA_ROOT"), "data/02_traj")
fig_out <- file.path(Sys.getenv("APA_ROOT"), "figures/02_traj")
ref <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/.ref/enrichment")
hpc_pattern <- paste0(Sys.getenv("APA_HPC_HOME"), "/project/work/apa/02_pattern")
gene_file <- file.path(out, "05_gene_aggregation.tsv")
enrichment_file <- file.path(out, "08_enrich_two_levels_1016.tsv")
adjusted_file <- file.path(out, "08_enrich_adjusted_outputs.tsv")
cov_file <- file.path(out, "00_cov.tsv")
expr_file <- file.path(Sys.getenv("APA_ROOT"), "data/01_quant/expr_after_combat.tsv")
tx2g_file <- file.path(hpc_pattern, "data/process/gencode.v40.tx2gene.tsv")
gmt_files <- file.path(hpc_pattern, "data/ref", c(
  "c5.go.bp.v2024.1.Hs.symbols.gmt",
  "c5.go.mf.v2024.1.Hs.symbols.gmt",
  "c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"
))
classes <- c(C1spec = "C1-spec", C2spec = "C2-spec", Shared = "Shared")
stop_words <- c("process", "regulation", "positive", "negative", "cellular", "activity",
                "of", "to", "and", "in", "by", "protein", "metabolic", "biosynthetic")

# --- helpers ---

clean_name <- function(x) {
  x <- toupper(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_+|_+$", "", x)
}

descendants <- function(root) {
  ids <- root
  repeat {
    updated <- unique(c(ids, syngo_ont$id[syngo_ont$parent_id %in% ids]))
    if (length(updated) == length(ids)) return(ids)
    ids <- updated
  }
}

# --- loading data ---

genes <- fread(gene_file)
enr <- fread(enrichment_file)[source_group == "GO/KEGG"]
adj <- fread(adjusted_file)[source_group == "GO/KEGG"]
gene_sets <- c(lapply(gmt_files, function(gmt) {
  lines <- strsplit(readLines(gmt), "\t", fixed = TRUE)
  setNames(lapply(lines, function(x) unique(x[-c(1, 2)])),
           vapply(lines, `[`, character(1), 1L))
}), recursive = FALSE)
gene_sets <- do.call(c, gene_sets)
tx2g <- fread(tx2g_file, select = c("gene_id", "gene_name", "start", "end"))
tx2g[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
gene_span <- tx2g[, .(gene_length = max(end) - min(start) + 1),
                  by = .(gene_symbol = gene_name, gene_id_base)]
expr_ids <- fread(cmd = sprintf("cut -f1 %s", shQuote(expr_file)), header = FALSE)[[1L]]
expr_ids <- sub("\\.[0-9]+$", "", expr_ids[-1L])
covariate_genes <- merge(genes, gene_span, by = c("gene_symbol", "gene_id_base"))[
  is.finite(gene_length) & gene_id_base %in% expr_ids]

cov <- fread(cov_file)
expr <- fread(expr_file, showProgress = FALSE)
setnames(expr, 1L, "gene_id")
expr[, gene_base := sub("\\.[0-9]+$", "", gene_id)]
expr <- expr[gene_base %in% covariate_genes$gene_id_base]
expr_samples <- intersect(names(expr), cov$sample)
expr_weights <- cov[match(expr_samples, sample), profile_weight]
expr_matrix <- as.matrix(expr[, ..expr_samples])
storage.mode(expr_matrix) <- "numeric"
expr[, mean_expression := as.numeric(expr_matrix %*% expr_weights / sum(expr_weights))]
expr_mean <- unique(expr[, .(gene_base, mean_expression)], by = "gene_base")
covariate_genes <- merge(covariate_genes, expr_mean, by.x = "gene_id_base",
                         by.y = "gene_base", all.x = TRUE)
covariate_genes <- covariate_genes[is.finite(mean_expression)]
stopifnot(uniqueN(covariate_genes$gene_symbol) == 11322L)

# --- term clusters ---

adj_directional <- rbindlist(list(
  adj[level == "Five-state" & gene_set == "C1-spec"],
  adj[level == "Membership" & gene_set == "Shared"],
  adj[level == "Five-state" & gene_set == "C2-spec"]
))
adj_directional[, level := "Directional membership"]
adj_directional[, adjusted_fdr := p.adjust(adjusted_pvalue, "BH"), by = source_group]
adj_all <- rbindlist(list(adj[level == "Five-state"], adj_directional))
adj_all <- unique(adj_all, by = c("level", "gene_set", "term"))
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)
cluster_summary <- list()
for (class_name in names(classes)) {
  group <- classes[[class_name]]
  selected_level <- if (class_name == "Shared") "Directional membership" else "Five-state"
  fg <- covariate_genes[if (class_name == "Shared") membership == "Shared" else gene_state == group,
                        gene_symbol]
  for (cutoff in c(0.05, 0.10)) {
    suffix <- if (cutoff == 0.05) "fdr05" else "fdr10"
    dt <- merge(enr[level == selected_level & gene_set == group & fdr_by_gene_set < cutoff],
                adj_all[level == selected_level & gene_set == group,
                        .(term, adjusted_odds_ratio, adjusted_fdr)], by = "term", all.x = TRUE)
    stopifnot(nrow(dt) > 0L, !anyDuplicated(dt$term))

    # gene-set Jaccard distance on foreground genes, average linkage
    setorder(dt, term)
    overlaps <- lapply(dt$term, function(term) intersect(gene_sets[[term]], fg))
    dt[, n_fg_covariate_bg := lengths(overlaps)]
    n <- nrow(dt)
    if (n > 1L) {
      distances <- outer(seq_len(n), seq_len(n), Vectorize(function(i, j) {
        1 - length(intersect(overlaps[[i]], overlaps[[j]])) /
          length(union(overlaps[[i]], overlaps[[j]]))
      }))
      dimnames(distances) <- list(dt$term, dt$term)
      tree <- hclust(as.dist(distances), method = "average")
      memberships <- lapply(c(0.6, 0.7, 0.8), function(h) cutree(tree, h = h))
    } else {
      tree <- NULL
      memberships <- rep(list(1L), 3L)
    }
    for (i in seq_along(memberships)) {
      values <- memberships[[i]]
      values <- match(values, unique(values))
      dt[, (c("cluster_id_d06", "cluster_id", "cluster_id_d08")[i]) :=
           sprintf("%s_%02d", class_name, values)]
    }
    dt[, `:=`(cluster_label_top_term = term[which.max(n_fg)],
              cluster_label_words = {
                words <- unlist(strsplit(tolower(gsub("^(GOBP|GOMF|KEGG)_", "", term)), "_"))
                words <- words[nchar(words) > 2 & !words %in% stop_words]
                if (!length(words)) "" else {
                  counts <- sort(table(words), decreasing = TRUE)
                  paste(names(counts)[seq_len(min(3L, length(counts)))], collapse = ", ")
                }
              }), by = cluster_id]
    columns <- c("term", "cluster_id", "cluster_id_d06", "cluster_id_d08",
                 "cluster_label_top_term", "cluster_label_words", "odds_ratio",
                 "fdr_by_gene_set", "adjusted_odds_ratio", "adjusted_fdr", "n_fg",
                 "n_fg_covariate_bg")
    fwrite(dt[, ..columns], file.path(out, sprintf("T1_term_clusters_%s_%s.tsv", class_name, suffix)), sep = "\t")
    if (suffix == "fdr05") {
      png(file.path(fig_out, sprintf("T1_treeplot_%s.png", class_name)),
          width = 2400, height = max(1100, 38 * n), res = 180)
      par(mar = c(4, 4, 3, 29))
      if (is.null(tree)) {
        plot.new()
        text(0.5, 0.5, dt$term)
      } else {
        labels <- gsub("_", " ", sub("^(GOBP|GOMF|KEGG)_", "", dt$term))
        tree$labels <- ifelse(nchar(labels) > 38,
                              paste0(substr(labels, 1L, 37L), "..."), labels)
        plot(as.dendrogram(tree), horiz = TRUE, cex = if (n > 45) 0.45 else 0.72,
             main = sprintf("%s: %d terms, %d clusters at d = 0.7", class_name,
                            n, uniqueN(dt$cluster_id)), xlab = "Jaccard distance")
        abline(v = 0.7, col = "#c45b42", lty = 2)
      }
      dev.off()
    }
    if (cutoff == 0.05) cluster_summary[[class_name]] <- dt
  }
}

# --- GO slim mapping ---

slim_lines <- readLines(file.path(out, "goslim_generic.obo"))
slim_ids <- unique(sub("^id: ", "", slim_lines[grepl("^id: GO:", slim_lines)]))
go_info <- AnnotationDbi::select(GO.db, keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
                                 columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
setDT(go_info)
go_info[, msigdb_term := paste0(ifelse(ONTOLOGY == "BP", "GOBP_", "GOMF_"), clean_name(TERM))]
go_info <- go_info[ONTOLOGY %in% c("BP", "MF")]
selected <- unique(rbindlist(cluster_summary)[, .(term)])
selected <- unique(c(selected$term,
                     enr[source_group == "GO/KEGG" & fdr_by_gene_set < 0.10, term]))
go_terms <- selected[grepl("^GO(BP|MF)_", selected)]
go_map <- merge(data.table(term = go_terms), go_info[, .(term = msigdb_term, go_id = GOID,
                                                          go_name = TERM, ontology = ONTOLOGY)],
                by = "term", all.x = TRUE, allow.cartesian = TRUE)
go_rows <- rbindlist(lapply(seq_len(nrow(go_map)), function(i) {
  row <- go_map[i]
  if (is.na(row$go_id)) return(data.table(term = row$term, go_id = NA_character_,
                                        slim_id = NA_character_, slim_name = NA_character_,
                                        status = "unmapped_go_name"))
  ancestors <- if (row$ontology == "BP") GOBPANCESTOR[[row$go_id]] else GOMFANCESTOR[[row$go_id]]
  matches <- intersect(c(row$go_id, ancestors), slim_ids)
  if (!length(matches)) return(data.table(term = row$term, go_id = row$go_id,
                                         slim_id = NA_character_, slim_name = NA_character_,
                                         status = "no_generic_slim_ancestor"))
  data.table(term = row$term, go_id = row$go_id, slim_id = matches,
             slim_name = go_info$TERM[match(matches, go_info$GOID)], status = "mapped")
}))
setorder(go_rows, term, slim_id)
fwrite(unique(go_rows), file.path(out, "T2_term_to_goslim.tsv"), sep = "\t")

# --- KEGG BRITE mapping ---

brite <- jsonlite::fromJSON(file.path(out, "kegg_br08901.json"), simplifyVector = FALSE)
brite_rows <- rbindlist(lapply(brite$children, function(level1) {
  rbindlist(lapply(level1$children, function(level2) {
    rbindlist(lapply(level2$children, function(pathway) {
      title <- sub("^[0-9]+ +", "", pathway$name)
      data.table(brite_level1 = level1$name, brite_level2 = level2$name,
                 pathway_id = sub("^([0-9]+).*", "\\1", pathway$name),
                 pathway_name = title, match_key = clean_name(title))
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
kegg_terms <- selected[grepl("^KEGG_", selected)]
kegg_map <- merge(data.table(term = kegg_terms,
                             match_key = sub("S_DISEASE$", "_DISEASE",
                                             sub("^KEGG_", "", kegg_terms))),
                  brite_rows[, match_key := sub("S_DISEASE$", "_DISEASE", match_key)],
                  by = "match_key", all.x = TRUE, allow.cartesian = TRUE)
kegg_map[, status := fifelse(is.na(pathway_id), "unmapped_kegg_name", "mapped")]
setcolorder(kegg_map, c("term", "brite_level1", "brite_level2", "pathway_id",
                        "pathway_name", "status", "match_key"))
setorder(kegg_map, term, brite_level1, brite_level2)
fwrite(kegg_map, file.path(out, "T2_kegg_to_brite.tsv"), sep = "\t")

# --- broad-set enrichment ---

mito <- suppressWarnings(readxl::read_excel(file.path(ref, "Human.MitoCarta3.0.xls"),
                                        sheet = "A Human MitoCarta3.0"))
syngo_dir <- tempfile("syngo_")
dir.create(syngo_dir)
unzip(file.path(ref, "syngo1.3_complete_data.zip"), exdir = syngo_dir)
syngo_genes <- readxl::read_excel(file.path(syngo_dir, "genes.xlsx"))
syngo_ann <- readxl::read_excel(file.path(syngo_dir, "annotations.xlsx"))
syngo_ont <- readxl::read_excel(file.path(syngo_dir, "ontologies.xlsx"))
rpl <- fread(file.path(ref, "HGNC_L_ribosomal_proteins.txt"))
rps <- fread(file.path(out, "HGNC_S_ribosomal_prots.tsv"))
proteasome <- fread(file.path(ref, "HGNC_Proteasome.tsv"))
stopifnot(nrow(mito) == 1136L, nrow(rpl) == 54L,
          nrow(rps) == 35L, nrow(proteasome) == 45L,
          all(rpl$`Group ID` == 729L), all(rps$`Group ID` == 728L),
          all(proteasome$`Group ID` == 690L))
broad_sets <- list(
  MitoCarta3 = unique(mito$Symbol),
  SynGO_all = unique(c(syngo_genes$hgnc_symbol, syngo_ann$hgnc_symbol)),
  SynGO_presynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098793")]),
  SynGO_postsynapse = unique(syngo_ann$hgnc_symbol[
    syngo_ann$go_id %in% descendants("GO:0098794")]),
  HGNC_RPL_RPS = unique(c(rpl$`Approved symbol`, rps$`Approved symbol`)),
  HGNC_Proteasome = unique(proteasome$`Approved symbol`)
)
broad_sets <- lapply(broad_sets, function(x) unique(x[!is.na(x) & x != ""]))
stopifnot(all(lengths(broad_sets) >= 10L))

t3 <- rbindlist(lapply(names(classes), function(class_name) {
  group <- classes[[class_name]]
  group_genes <- if (class_name == "Shared") {
    covariate_genes$membership == "Shared"
  } else {
    covariate_genes$gene_state == group
  }
  rbindlist(lapply(names(broad_sets), function(set_name) {
    set_genes <- broad_sets[[set_name]]
    term_genes <- covariate_genes$gene_symbol %in% set_genes
    a <- sum(group_genes & term_genes)
    b <- sum(group_genes & !term_genes)
    c <- sum(!group_genes & term_genes)
    d <- sum(!group_genes & !term_genes)
    contingency <- matrix(c(a, b, c, d), nrow = 2)
    fisher <- fisher.test(contingency, alternative = "greater")
    fisher_ci <- fisher.test(contingency)$conf.int
    fit_data <- copy(covariate_genes)
    fit_data[, `:=`(in_group = as.integer(group_genes), in_term = as.integer(term_genes))]
    fit <- glm(in_group ~ in_term + scale(log1p(n_sites)) +
                 scale(log1p(gene_length)) + scale(mean_expression),
               family = binomial(), data = fit_data)
    estimate <- coef(summary(fit))["in_term", ]
    data.table(
      gene_class = class_name, gene_set = set_name,
      source_version = if (grepl("^SynGO", set_name)) "SynGO 1.3" else
        if (set_name == "MitoCarta3") "MitoCarta 3.0" else "HGNC group download 2026-10-16",
      n_set_source = length(set_genes), n_set_background = a + c,
      n_fg = a, n_fg_total = a + b, n_background = a + b + c + d,
      fisher_odds_ratio = unname(fisher$estimate),
      fisher_ci95_low = fisher_ci[1L], fisher_ci95_high = fisher_ci[2L],
      fisher_p = fisher$p.value,
      adjusted_odds_ratio = exp(estimate["Estimate"]),
      adjusted_ci95_low = exp(estimate["Estimate"] - 1.96 * estimate["Std. Error"]),
      adjusted_ci95_high = exp(estimate["Estimate"] + 1.96 * estimate["Std. Error"]),
      adjusted_p = estimate["Pr(>|z|)"]
    )
  }))
}))
t3[, `:=`(fisher_fdr = p.adjust(fisher_p, "BH"),
           adjusted_fdr = p.adjust(adjusted_p, "BH"))]
stopifnot(all(is.finite(t3$adjusted_p)))
fwrite(t3, file.path(out, "T3_broad_sets_enrich.tsv"), sep = "\t")

# --- results and provenance ---

report <- c(
  "# Fig 1e enrichment hierarchy",
  "",
  "GO/KEGG terms use Fisher `fdr_by_gene_set < 0.05` for the main clusters. Clustering uses the 11,322 genes with complete covariates, foreground-intersected gene-set Jaccard distance, average linkage, and a 0.7 tree cut. The TSVs also report cuts at 0.6 and 0.8; separate FDR < 0.10 tables are provided. `n_fg` is the original Fisher count; `n_fg_covariate_bg` is the clustering overlap count.",
  "",
  "GO slim mapping uses GO.db ancestors and the downloaded generic slim; terms without a generic slim ancestor remain explicitly unmapped. KEGG mapping uses the downloaded BRITE br08901 hierarchy. The 2024 MSigDB possessive disease names are normalized to current KEGG titles.",
  ""
)
for (class_name in names(cluster_summary)) {
  dt <- cluster_summary[[class_name]]
  groups <- unique(dt[, .(cluster_id, cluster_label_top_term, cluster_label_words)])
  groups[, n_terms := dt[, .N, by = cluster_id]$N[match(cluster_id,
                       dt[, .N, by = cluster_id]$cluster_id)]]
  setorder(groups, cluster_id)
  report <- c(report, sprintf("## %s: %d terms, %d clusters", class_name,
                               nrow(dt), nrow(groups)), "")
  report <- c(report, switch(class_name,
    C1spec = "The C1-specific terms separate translation/ribosome overlap from a GPCR activity singleton.",
    C2spec = "The C2-specific clusters include respiration/energy production, nucleotide-related metabolism, and separate proteostasis terms.",
    Shared = "The Shared clusters span cytoskeleton/junction organization, transport and hormone response, with several neural-development terms."
  ), "")
  for (i in seq_len(nrow(groups))) {
    row <- groups[i]
    title <- gsub("_", " ", sub("^(GOBP|GOMF|KEGG)_", "", row$cluster_label_top_term))
    sentence <- if (row$n_terms == 1L) {
      sprintf("This singleton is the %s term.", tolower(title))
    } else {
      sprintf("These %d overlapping terms center on %s, with frequent title words: %s.",
              row$n_terms, tolower(title), row$cluster_label_words)
    }
    report <- c(report, sprintf("- `%s` (%d %s): %s", row$cluster_id,
                                row$n_terms, ifelse(row$n_terms == 1L, "term", "terms"), sentence))
  }
  report <- c(report, "")
}
report <- c(report,
  "The clusters are descriptive overlap groups; singleton terms and GO/KEGG pathway titles are not independent mechanistic confirmations.",
  "", "## Broad gene sets", "",
  "T3 uses the 11,322-gene complete-covariate universe for one-sided Fisher tests and logistic models with `scale(log1p(n_sites))`, `scale(log1p(gene_length))`, and `scale(mean_expression)`. Fisher P values are one-sided; logistic P values and Wald 95% CIs are two-sided. BH FDR is calculated separately for all 18 Fisher tests and all 18 adjusted tests.",
  "",
  "SynGO input is the supplied version 1.3 archive (the original request named 1.2). The all-synapse set is the union of symbols in `genes.xlsx` and `annotations.xlsx`; pre/post sets use annotations for the specified GO roots and all descendants in `ontologies.xlsx`. The supplied HGNC S file duplicated the L group, so RPS uses a fresh HGNC group 728 snapshot in this handoff.",
  "")
for (class_name in names(classes)) {
  hits <- t3[gene_class == class_name & fisher_fdr < 0.05]
  adjusted_up <- t3[gene_class == class_name & adjusted_fdr < 0.05 & adjusted_odds_ratio > 1]
  adjusted_down <- t3[gene_class == class_name & adjusted_fdr < 0.05 & adjusted_odds_ratio < 1]
  report <- c(report, sprintf("- %s: Fisher enrichment in %s; adjusted enrichment in %s; adjusted depletion in %s (FDR < 0.05).",
                              class_name,
                              if (nrow(hits)) paste(hits$gene_set, collapse = ", ") else "none",
                              if (nrow(adjusted_up)) paste(adjusted_up$gene_set, collapse = ", ") else "none",
                              if (nrow(adjusted_down)) paste(adjusted_down$gene_set, collapse = ", ") else "none"))
}
report <- c(report, "",
  "Shared SynGO_all changes from Fisher enrichment to adjusted depletion; these are different estimands and should not be presented as concordant evidence.",
  "")
writeLines(report, file.path(out, "results_handoff2.md"))

# --- manifest ---

output_files <- c(
  unlist(lapply(names(classes), function(x) c(
    sprintf("T1_term_clusters_%s_fdr05.tsv", x),
    sprintf("T1_term_clusters_%s_fdr10.tsv", x),
    sprintf("T1_treeplot_%s.png", x)
  ))),
  "T2_term_to_goslim.tsv", "T2_kegg_to_brite.tsv",
  "T3_broad_sets_enrich.tsv", "results_handoff2.md",
  "goslim_generic.obo", "kegg_br08901.json",
  "HGNC_S_ribosomal_prots.tsv", "enrich_hier.R"
)
source_paths <- fifelse(grepl("^T1_", output_files),
                       paste(c(gene_file, enrichment_file, adjusted_file, gmt_files),
                             collapse = ";"),
                       fifelse(grepl("^T2_term", output_files),
                               paste("GO.db", file.path(out, "goslim_generic.obo"),
                                     enrichment_file, sep = ";"),
                               fifelse(grepl("^T2_kegg", output_files),
                                       paste(file.path(out, "kegg_br08901.json"),
                                             enrichment_file, sep = ";"), "")))
source_paths[output_files == "goslim_generic.obo"] <-
  "https://current.geneontology.org/ontology/subsets/goslim_generic.obo"
source_paths[output_files == "kegg_br08901.json"] <-
  "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json"
source_paths[output_files == "HGNC_S_ribosomal_prots.tsv"] <-
  "https://www.genenames.org/cgi-bin/genegroup/download?id=728&type=branch"
source_paths[output_files == "T3_broad_sets_enrich.tsv"] <-
  paste(c(gene_file, cov_file, expr_file,
          file.path(ref, c("Human.MitoCarta3.0.xls", "syngo1.3_complete_data.zip",
                           "HGNC_L_ribosomal_proteins.txt", "HGNC_Proteasome.tsv")),
          file.path(out, "HGNC_S_ribosomal_prots.tsv")), collapse = ";")
source_paths[output_files == "results_handoff2.md"] <- paste(gene_file, enrichment_file, sep = ";")
source_paths[output_files == "enrich_hier.R"] <-
  "codex_prompt_enrichment_hierarchy_1016.md"
manifest <- data.table(file = output_files, source_path = source_paths,
                       script = file.path(Sys.getenv("APA_ROOT"), "release/02_traj/followup/enrich_hier.R"),
                       job_id = Sys.getenv("SLURM_JOB_ID"),
                       compute_host = Sys.info()[["nodename"]])
fwrite(manifest, file.path(out, "README_handoff2.tsv"), sep = "\t")
