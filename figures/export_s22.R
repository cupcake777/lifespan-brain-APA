# Export the Fig 2d null t-statistic matrix (RDS, 20334 pairs x 100 null sets) to gzipped TSV.
ROOT <- Sys.getenv("APA_ROOT")
x <- readRDS(file.path(ROOT, "data/03_omics/s22a_null_t.rds"))
stopifnot(is.matrix(x), is.numeric(x))
colnames(x) <- paste0("null_set_", seq_len(ncol(x)))
con <- gzfile(file.path(ROOT, "data/03_omics/s22a_null_t.tsv.gz"), "w")
writeLines(paste(c("pair_row", colnames(x)), collapse = "\t"), con)
for (i in seq_len(nrow(x))) writeLines(paste(c(i, formatC(x[i, ], digits = 15, format = "g")), collapse = "\t"), con)
close(con)
