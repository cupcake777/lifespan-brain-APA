#!/usr/bin/env python3
# Write the SMR qfile of grouped APA-QTL nominal pairs for one QTL type.
# af is treated as the ALT allele frequency: A1 = ALT, A2 = REF, Freq = A1 frequency.
import argparse
import gzip
import math
from pathlib import Path

# --- helpers ---
def parse_finite_number(value):
    if value == "":
        return None
    try:
        num = float(value)
    except ValueError:
        return None
    return value if math.isfinite(num) else None


# --- args ---
p = argparse.ArgumentParser()
p.add_argument("--tabix", required=True, help="tabix .gz file (all pairs)")
p.add_argument("--out-dir", required=True, help="output directory for grouped SMR qfile input")
p.add_argument("--qtl-type", required=True, choices=["prenatal", "postnatal", "interaction"], help="QTL type label")
args = p.parse_args()

out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)
qfile_path = out_dir / f"{args.qtl_type}.qfile.txt"

# --- write qfile ---
n_total = n_written = n_bad_freq = n_bad_effect = 0
probes = set()
with gzip.open(args.tabix, "rt") as fh, open(qfile_path, "w") as qfile:
    qfile.write("SNP Chr BP A1 A2 Freq Probe Probe_Chr Probe_bp Gene Orientation b se p\n")
    header = fh.readline()
    if not header:
        raise SystemExit(f"ERROR: input file is empty: {args.tabix}")

    for line in fh:
        n_total += 1
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 14:
            continue

        chr_, pos, ref, alt = cols[0], cols[1], cols[2], cols[3]
        phenotype_id = cols[4]
        variant_id = cols[5]
        try:
            freq = float(cols[8])
            freq_ok = 0.0 <= freq <= 1.0
        except ValueError:
            freq_ok = False
        pval = parse_finite_number(cols[11])
        slope = parse_finite_number(cols[12])
        slope_se = parse_finite_number(cols[13])

        if not freq_ok:
            n_bad_freq += 1
            continue
        if pval is None or slope is None or slope_se is None:
            n_bad_effect += 1
            continue

        parts = phenotype_id.split("|")
        gene = parts[1] if len(parts) >= 2 else phenotype_id
        locus = parts[2] if len(parts) >= 3 else "chrUn:0-1"
        strand = parts[3] if len(parts) >= 4 and parts[3] in {"+", "-"} else "NA"
        probe_chr, coords = locus.split(":")
        probe_bp = coords.split("-", 1)[0]
        snp = variant_id.replace("chr", "", 1)
        chrom_num = chr_.replace("chr", "")
        probe_chr_num = probe_chr.replace("chr", "")

        qfile.write(
            f"{snp} {chrom_num} {pos} {alt} {ref} {cols[8]} "
            f"{phenotype_id} {probe_chr_num} {probe_bp} {gene} {strand} "
            f"{slope} {slope_se} {pval}\n"
        )
        probes.add(phenotype_id)
        n_written += 1

print(f"Read {n_total} rows, kept {n_written} QTL rows across {len(probes)} probes -> {qfile_path}")
print(f"Skipped: bad freq {n_bad_freq}, missing effect {n_bad_effect}")
