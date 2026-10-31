#!/usr/bin/env python3
# Convert a harmonized GWAS table to sLDSC sumstats (SNP, A1, A2, N, Z, FRQ).
# EA is the effect allele for BETA; SNP = SNPID, or variant_id when SNPID is empty.
# Usage: 05_format_gwas.py <harmonized.tsv> <out.sumstats.gz>
import csv
import gzip
import sys

in_tsv, out_gz = sys.argv[1], sys.argv[2]

n_written = n_skipped = 0
with open(in_tsv, "r") as fin, gzip.open(out_gz, "wt") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["SNP", "A1", "A2", "N", "Z", "FRQ"])
    for row in reader:
        # rows with a missing column or non-numeric BETA/SE are skipped
        try:
            snp = row["SNPID"] if row["SNPID"] else row["variant_id"]
            beta = float(row["BETA"])
            se = float(row["SE"])
            if se <= 0:
                n_skipped += 1
                continue
            writer.writerow([snp, row["EA"], row["NEA"], row["N"], f"{beta / se:.6f}", row["EAF"]])
            n_written += 1
        except (KeyError, ValueError):
            n_skipped += 1

print(f"Written {n_written} variants to {out_gz}")
print(f"Skipped {n_skipped} variants due to missing/invalid data")
