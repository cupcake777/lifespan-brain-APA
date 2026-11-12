#!/usr/bin/env python3
# Convert EasyColoc harmonized GWAS files (*_harmonized.tsv) to SMR .ma input.
import argparse
import gzip
import sys
from pathlib import Path

SNPID_COLS = ["variant_id", "SNPID", "SNP", "MarkerID", "rsID"]
A1_COLS = ["EA", "A1", "ALT", "effect_allele"]
A2_COLS = ["NEA", "A2", "REF", "other_allele"]
FREQ_COLS = ["EAF", "AF", "Freq", "freq", "effect_allele_frequency"]
BETA_COLS = ["BETA", "b", "beta", "Effect"]
SE_COLS = ["SE", "se", "StdErr"]
P_COLS = ["P", "p", "p_value", "Pvalue"]
N_COLS = ["N", "n", "N_analyzed", "TotalN"]
CHR_COLS = ["CHR", "chrom", "chromosome"]
POS_COLS = ["POS", "BP", "bp", "position"]


# --- helpers ---
def first_col(header, candidates):
    for c in candidates:
        if c in header:
            return c
    return None


def get(row, col, c, default="NA"):
    return row[col[c]] if c and c in col and col[c] < len(row) else default


# --- args ---
p = argparse.ArgumentParser()
p.add_argument("--harmony-dir", required=True, help="EasyColoc harmony directory")
p.add_argument("--out-dir", required=True, help="output directory for .ma files")
p.add_argument(
    "--include-prefix",
    action="append",
    default=[],
    help="only convert datasets with this prefix; can be repeated",
)
args = p.parse_args()

harmony_dir = Path(args.harmony_dir)
out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

files = sorted(harmony_dir.glob("*_harmonized.tsv"))
if not files:
    sys.exit(f"No *_harmonized.tsv files found in {harmony_dir}")

# --- convert each GWAS ---
for src in files:
    gwas_id = src.stem.removesuffix("_harmonized").split("_b", 1)[0]
    if args.include_prefix and not any(gwas_id.startswith(prefix) for prefix in args.include_prefix):
        continue
    out_path = out_dir / f"{gwas_id}.ma"
    if out_path.exists():
        print(f"  skip {gwas_id} (exists)")
        continue

    opener = gzip.open if str(src).endswith(".gz") else open
    with opener(src, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}

        snp_col = first_col(col, SNPID_COLS)
        a1_col = first_col(col, A1_COLS)
        a2_col = first_col(col, A2_COLS)
        freq_col = first_col(col, FREQ_COLS)
        beta_col = first_col(col, BETA_COLS)
        se_col = first_col(col, SE_COLS)
        p_col = first_col(col, P_COLS)
        n_col = first_col(col, N_COLS)
        chr_col = first_col(col, CHR_COLS)
        pos_col = first_col(col, POS_COLS)

        missing = [
            name
            for name, c in [
                ("SNP", snp_col),
                ("A1", a1_col),
                ("A2", a2_col),
                ("freq", freq_col),
                ("BETA", beta_col),
                ("SE", se_col),
                ("P", p_col),
            ]
            if c is None
        ]
        if missing:
            print(f"  SKIP {src.name}: missing columns {missing}")
            continue

        # frequency must refer to the A1 allele
        effect_freq = freq_col in {"EAF", "effect_allele_frequency"} and a1_col in {"EA", "effect_allele"}
        freq_is_a1 = effect_freq or freq_col in {"Freq", "freq"} or (freq_col == "AF" and a1_col == "ALT")
        if not freq_is_a1:
            print(f"  SKIP {src.name}: cannot prove {freq_col} is frequency of A1={a1_col}")
            continue

        n_written = n_missing_value = n_bad_freq = n_unmapped = 0
        with open(out_path, "w") as out:
            out.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tn\n")
            for line in fh:
                row = line.rstrip("\n").split("\t")
                snp = get(row, col, snp_col)
                a1 = get(row, col, a1_col)
                a2 = get(row, col, a2_col)
                freq = get(row, col, freq_col)
                b = get(row, col, beta_col)
                se = get(row, col, se_col)
                p = get(row, col, p_col)
                n = get(row, col, n_col, "NA")
                chrom = get(row, col, chr_col, "")
                pos = get(row, col, pos_col, "")

                if "NA" in (snp, a1, a2, freq, b, se, p) or snp in ("", "."):
                    n_missing_value += 1
                    continue
                try:
                    freq_ok = 0.0 <= float(freq) <= 1.0
                except ValueError:
                    freq_ok = False
                if not freq_ok:
                    n_bad_freq += 1
                    continue
                # normalize to chr:pos (keep chr:pos:ref:alt IDs as they are)
                snp = snp.removeprefix("chr")
                if snp.count(":") != 3:
                    chrom = chrom.removeprefix("chr")
                    if chrom and pos:
                        snp = f"{chrom}:{pos}"
                    elif snp.count(":") != 1:
                        n_unmapped += 1
                        continue
                out.write(f"{snp}\t{a1}\t{a2}\t{freq}\t{b}\t{se}\t{p}\t{n}\n")
                n_written += 1

    print(
        f"  {gwas_id}: {n_written:,} SNPs -> {out_path} "
        f"(missing/NA: {n_missing_value:,}; bad freq: {n_bad_freq:,}; unmapped: {n_unmapped:,})"
    )
