#!/usr/bin/env Rscript


library(jsonlite)
library(glue)

HASH_DIR <- Sys.getenv("EASYCOLOC_HASH_JSON_DIR", unset = "~/work/coloc/snp_ref")
if (!dir.exists(HASH_DIR)) stop("Directory not found!")

files <- list.files(HASH_DIR, pattern = "chr_.*_snp151_hash_table\\.json$", full.names = TRUE)

cat(glue("Found {length(files)} JSON files to convert.\n\n"))

for (f in files) {
  f_out <- sub("\\.json$", ".rds", f)
  
  fname <- basename(f)
  
  if (file.exists(f_out)) {
    cat(glue("[SKIP] {fname} -> RDS already exists\n"))
    next
  }
  
  cat(glue("[CONVERT] Reading {fname}... This may take a minute...\n"))
  
  t1 <- Sys.time()
  data <- tryCatch({
    fromJSON(txt = f)
  }, error = function(e) {
    cat(glue("  ERROR reading {fname}: {e$message}\n"))
    return(NULL)
  })
  
  if (is.null(data)) next
  
  t2 <- Sys.time()
  read_time <- round(as.numeric(difftime(t2, t1, units="secs")), 1)
  
  cat(glue("  Loaded in {read_time}s. Saving to RDS...\n"))
  
  saveRDS(data, file = f_out, compress = TRUE) 
  
  t3 <- Sys.time()
  write_time <- round(as.numeric(difftime(t3, t2, units="secs")), 1)
  
  cat(glue("  Saved in {write_time}s. Total: {read_time + write_time}s\n"))
  
  rm(data)
  gc()
}

cat("\nAll done! You can now use the updated utils_hash.R for fast loading.\n")
