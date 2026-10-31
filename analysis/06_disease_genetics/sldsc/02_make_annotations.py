#!/usr/bin/env python3
# Build +/-window bp BED intervals around QTL variants, merged per chromosome.
# Usage: 02_annot_bed.py <signif_pairs.txt.gz> <out.bed> [window_bp]
# Coordinates come from chr/pos columns, or from variant_id (chr:pos:ref:alt).
import csv
import gzip
import sys
from collections import defaultdict

in_gz, out_bed = sys.argv[1], sys.argv[2]
window = int(sys.argv[3]) if len(sys.argv) > 3 else 500

# --- unique positions per chromosome ---
positions = defaultdict(set)
with gzip.open(in_gz, "rt") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("chr") and row.get("pos"):
            chromosome, position = row["chr"], row["pos"]
        else:
            chromosome, position = row["variant_id"].split(":", 2)[:2]
        if not chromosome.startswith("chr"):
            chromosome = "chr" + chromosome
        positions[chromosome].add(int(position))

# --- merge windows and write BED, chromosomes in numeric order ---
unique_positions = sum(len(values) for values in positions.values())
interval_count = 0
with open(out_bed, "w") as out:
    for chromosome in sorted(positions, key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c[3:])):
        merged = []
        for position in sorted(positions[chromosome]):
            start = max(0, position - 1 - window)
            end = position + window
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        interval_count += len(merged)
        for start, end in merged:
            out.write(f"{chromosome}\t{start}\t{end}\n")

print(f"Created {out_bed}: {unique_positions} unique positions, {interval_count} merged intervals, +/-{window} bp")
