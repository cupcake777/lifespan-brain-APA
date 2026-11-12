# Mean and SD of the apaQTL and eQTL annotations over common SNPs (MAF 0.05-0.5), per population and stage.
# Usage: 03-2_calc_annot.py <pop> <stage> <sldsc_ref> <ldscore_dir> <out.tsv>
import os
import sys

import numpy as np
import pandas as pd

pop, stage, sldsc_ref, ldscore_dir, out_file = sys.argv[1:6]
chromosomes = [int(value) for value in os.environ.get("CHROMOSOMES", " ".join(str(i) for i in range(1, 23))).split()]
records = []
for chrom in chromosomes:
    frq = pd.read_csv(f"{sldsc_ref}/{pop}_hg38/plink_files/1000G.{pop}.hg38.{chrom}.frq", sep=r"\s+", usecols=["SNP", "MAF"])
    annot = pd.read_csv(f"{ldscore_dir}/{pop}/{stage}_maxCPP.{chrom}.annot.gz", sep="\t")
    merged = pd.merge(annot, frq, on="SNP")
    common = merged[(merged["MAF"] > 0.05) & (merged["MAF"] <= 0.5)]
    for col in ["apaQTL", "eQTL"]:
        records.append({"annot": col, "vals": common[col].values})

res = []
for a in ["apaQTL", "eQTL"]:
    all_vals = np.concatenate([r["vals"] for r in records if r["annot"] == a])
    res.append({"population": pop, "stage": stage, "annot_set": "maxCPP", "annotation": a,
                "n_snps": len(all_vals), "mean": np.mean(all_vals), "sd": np.std(all_vals, ddof=0)})
pd.DataFrame(res).to_csv(out_file, sep="\t", index=False)
