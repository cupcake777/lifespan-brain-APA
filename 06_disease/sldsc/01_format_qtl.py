#!/usr/bin/env python3
# Convert APA-QTL signif pairs to sLDSC sumstats (SNP, A1, A2, N, Z, FRQ).
# Z = slope / slope_se, keeping the input slope direction.
# Usage: 01_format_qtl.py <signif_pairs.txt.gz> <out.sumstats.gz>
import csv
import gzip
import sys

in_gz, out_gz = sys.argv[1], sys.argv[2]

n_written = n_skipped = 0
with gzip.open(in_gz, "rt") as fin, gzip.open(out_gz, "wt") as fout:
    reader = csv.DictReader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["SNP", "A1", "A2", "N", "Z", "FRQ"])
    for row in reader:
        # rows with a missing column or non-numeric slope are skipped
        try:
            slope = float(row["slope"])
            slope_se = float(row["slope_se"])
            if slope_se <= 0:
                n_skipped += 1
                continue
            z = slope / slope_se
            writer.writerow([row["variant_id"], row["ref"], row["alt"], row["ma_samples"], f"{z:.6f}", row["af"]])
            n_written += 1
        except (KeyError, ValueError):
            n_skipped += 1

print(f"Written {n_written} variants to {out_gz}")
print(f"Skipped {n_skipped} variants due to missing/invalid data")
