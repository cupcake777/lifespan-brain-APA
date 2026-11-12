# Phenotype group file, GeneRange BED and tabix index per stage.
# Called by preprocess.R as: python3 format.py <input_dir> <output_dir>

import glob
import os
import re
import subprocess
import sys

import pandas as pd

# --- inputs ---

APA_ROOT = os.environ["APA_ROOT"]
default_dir = os.path.join(APA_ROOT, "data", "04_qtl", "pre_input")
input_dir = sys.argv[1] if len(sys.argv) > 1 else default_dir
output_dir = sys.argv[2] if len(sys.argv) > 2 else input_dir

# --- stages ---

# stages are detected from the phenotype files present
stage_pattern = re.compile(r'phenotype\.(\w+)\.bed\.gz')
stages = {}
for f in glob.glob(os.path.join(input_dir, "phenotype.*.bed.gz")):
    match = stage_pattern.search(os.path.basename(f))
    if match:
        stages[match.group(1).lower()] = match.group(1)

if not stages:
    stages = {'pre': 'Prenatal', 'post': 'Postnatal'}  # fallback default

# --- per-stage tables ---

for stage_name in stages.values():
    input_file = os.path.join(input_dir, f"phenotype.{stage_name}.bed.gz")
    if not os.path.exists(input_file):
        print(f"Skip: {input_file}")
        continue

    print(f"\n>>> Processing: {stage_name} <<<")
    df = pd.read_csv(input_file, sep='\t', compression='gzip')

    cols = df.columns.tolist()
    chr_col, start_col, end_col, pid_col = cols[0], cols[1], cols[2], cols[3]
    # gene id = second field of the phenotype id, else the whole id
    df['gid_raw'] = df[pid_col].apply(
        lambda pid: pid.split('|')[1]
        if isinstance(pid, str) and len(pid.split('|')) >= 2 else pid)
    locus_coords = df.groupby([chr_col, 'gid_raw']).agg({start_col: 'min', end_col: 'max'}).reset_index()
    locus_coords.columns = [chr_col, 'gid_raw', 'new_start', 'new_end']

    df = df.merge(locus_coords, on=[chr_col, 'gid_raw'], how='left')
    df[start_col] = df['new_start']
    df[end_col] = df['new_end']

    df['_locus_id'] = df[chr_col].astype(str) + ":" + df[start_col].astype(int).astype(str)

    locus_counts = df.groupby('gid_raw')['_locus_id'].nunique()
    multi_locus_genes = locus_counts[locus_counts > 1].index
    df['gid'] = df['gid_raw']
    mask = df['gid_raw'].isin(multi_locus_genes)
    df.loc[mask, 'gid'] = df.loc[mask, 'gid_raw'] + "|" + df.loc[mask, '_locus_id']
    df = df.sort_values(by=[chr_col, start_col, end_col, 'gid', pid_col])

    output_group = os.path.join(output_dir, f"phenotype.{stage_name}.group")
    df[[pid_col, 'gid']].to_csv(output_group, sep='\t', index=False, header=False)

    output_bed = os.path.join(output_dir, f"phenotype.{stage_name}.GeneRange.bed")
    bed_out = df[cols].copy()
    if not bed_out.columns[0].startswith('#'):
        bed_out.rename(columns={bed_out.columns[0]: '#' + bed_out.columns[0]}, inplace=True)

    bed_out.to_csv(output_bed, sep='\t', index=False, header=True)
    subprocess.run(f"bgzip -f {output_bed}", shell=True, check=True)
    subprocess.run(f"tabix -f -p bed {output_bed}.gz", shell=True, check=True)

    print(f"  - Done. Rows: {len(df)}, Multi-locus genes fixed: {len(multi_locus_genes)}")

print("\nAll tasks completed.")
