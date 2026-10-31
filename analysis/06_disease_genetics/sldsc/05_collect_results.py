#!/usr/bin/env python3
# Collect sLDSC .results files into one table; FDR_p = BH over all rows.
# Usage: 07_collect.py <regression_dir> <out.tsv>
import csv
import sys
from pathlib import Path

input_dir, output_file = sys.argv[1], sys.argv[2]
COLS = ["Coefficient", "Coefficient_std_error", "Coefficient_z-score",
        "Enrichment", "Enrichment_std_error", "Enrichment_p"]

results_files = list(Path(input_dir).glob("*.results"))
if not results_files:
    print(f"No .results files found in {input_dir}")
    sys.exit(0)
print(f"Found {len(results_files)} results files")

# --- one row per .results file: the last row is the custom annotation ---
rows = []
for results_file in sorted(results_files):
    # file name: EUR_<TRAIT>_<annotation>.results
    parts = results_file.stem.split("_")
    if parts[0] not in ("EUR", "EAS"):
        print(f"Warning: Cannot parse filename: {results_file.name}")
        continue
    data = {}
    with open(results_file, "r") as fin:
        for row in csv.DictReader(fin, delimiter="\t"):
            data = {c: row.get(c, "") for c in COLS}
    if not data:
        continue
    rows.append({"Trait": parts[1], "Population": parts[0], "Annotation": "_".join(parts[2:]), **data})

# --- BH over all valid Enrichment_p values ---
pvals = []
for row in rows:
    try:
        pvals.append(float(row["Enrichment_p"]))
    except (TypeError, ValueError):
        pvals.append(None)

valid = [(i, p) for i, p in enumerate(pvals) if p is not None]
valid_sorted = sorted(valid, key=lambda x: x[1])
fdr_values = [None] * len(rows)
m = len(valid_sorted)
prev_q = 1.0
for rank, (idx, pval) in reversed(list(enumerate(valid_sorted, start=1))):
    qval = min(prev_q, pval * m / rank)
    qval = min(qval, 1.0)
    fdr_values[idx] = qval
    prev_q = qval

for row, qval in zip(rows, fdr_values):
    row["FDR_p"] = "" if qval is None else f"{qval:.6g}"

# --- write summary table ---
header = ["Trait", "Population", "Annotation"] + COLS + ["FDR_p"]
with open(output_file, "w", newline="") as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(header)
    for row in rows:
        writer.writerow([row[c] for c in header])
print(f"Summary table written to: {output_file}")
