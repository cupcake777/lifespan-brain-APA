#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("src/utils_plot.R")

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

tmp_dir <- tempfile("easycoloc_ld_timeout_")
dir.create(tmp_dir, recursive = TRUE)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

bfile <- file.path(tmp_dir, "toy")
data.table::fwrite(
  data.table(
    CHR = c(1, 1),
    SNP = c("rs1", "rs2"),
    CM = c(0, 0),
    POS = c(100, 200),
    A1 = c("A", "C"),
    A2 = c("G", "T")
  ),
  paste0(bfile, ".bim"),
  sep = "\t",
  col.names = FALSE
)
invisible(file.create(paste0(bfile, ".bed")))

plink_stub <- file.path(tmp_dir, "plink_stub.sh")
writeLines(c(
  "#!/usr/bin/env bash",
  "sleep 2",
  "exit 0"
), plink_stub)
Sys.chmod(plink_stub, mode = "0755")

old_timeout <- getOption("easycoloc.ld_plink_timeout")
on.exit(options(easycoloc.ld_plink_timeout = old_timeout), add = TRUE)
options(easycoloc.ld_plink_timeout = 1L)

result <- ld_extract(c("rs1", "rs2"), bfile, plink_stub)
assert_true(is.null(result), "ld_extract should return NULL when PLINK times out")

cat("[SMOKE] LD timeout smoke test passed\n")
