#!/usr/bin/env python3
# Write one chromosome's LDSC annotation (.annot.gz): 1 if a .bim SNP lies in a merged BED window.
# Usage: 03_annot.py <windows.bed> <chr.bim> <annot_name> <out.annot.gz>
import gzip
import sys
from bisect import bisect_right
from collections import defaultdict

bed_path, bim_path, annot_name, out_path = sys.argv[1:5]

# --- merged windows per chromosome: sorted starts and matching ends ---
raw = defaultdict(list)
with open(bed_path) as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3:
            raw[fields[0]].append((int(fields[1]), int(fields[2])))

lookup = {}
for chromosome, intervals in raw.items():
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    lookup[chromosome] = ([start for start, _ in merged], [end for _, end in merged])

# --- annotate each BIM SNP (inclusive window overlap) ---
total = annotated = 0
with open(bim_path) as bim, gzip.open(out_path, "wt") as output:
    output.write(f"CHR\tBP\tSNP\tCM\t{annot_name}\n")
    for line in bim:
        fields = line.split()
        if len(fields) < 6:
            continue
        chromosome, snp, cm, bp = fields[:4]
        position = int(bp)
        value = 0
        if "chr" + chromosome in lookup:
            starts, ends = lookup["chr" + chromosome]
            index = bisect_right(starts, position) - 1
            value = int(index >= 0 and position <= ends[index])
        annotated += value
        total += 1
        output.write(f"{chromosome}\t{position}\t{snp}\t{cm}\t{value}\n")

print(f"Created {out_path}: {total} SNPs, {annotated} annotated")
