# S-LDSC inputs for the G x Stage maxCPP set: merged max-PIP table, annotation files, annotation SD.
# Usage: merge <out.tsv> | annot <bim> <pip.tsv> <out.annot.gz> | calcsd <pop> <sldsc_ref> <ldscore_dir> <out.tsv>
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(os.environ["APA_ROOT"]) / "data" / "04_qtl" / "interaction_sldsc"
mode = sys.argv[1]

if mode == "merge":
    dfs = [pd.read_csv(BASE / "finemap" / f"chr{chrom}.maxpip.tsv.gz", sep="\t") for chrom in range(1, 23)]
    df = pd.concat(dfs, ignore_index=True).dropna()
    df = df.groupby("variant_id", as_index=False)["max_pip"].max()
    splits = df["variant_id"].str.split(":", expand=True)
    df["chr"] = splits[0].str.replace("chr", "")
    df["pos"] = splits[1]
    df["a1"] = splits[2].str.upper()
    df["a2"] = splits[3].str.upper()
    df = df[df["chr"].str.isdigit()].copy()
    df["chr_num"] = df["chr"].astype(int)
    df["pos_num"] = df["pos"].astype(int)
    df = df.sort_values(["chr_num", "pos_num"])
    df[["chr", "pos", "a1", "a2", "max_pip"]].to_csv(sys.argv[2], sep="\t", index=False, header=False, float_format="%.6g")

elif mode == "annot":
    bim_file, pip_tsv, out_gz = sys.argv[2:5]
    bim = pd.read_csv(bim_file, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    pip = pd.read_csv(pip_tsv, sep="\t", header=None, names=["CHR", "BP", "A1", "A2", "interactionQTL"])
    k1 = pip["CHR"].astype(str) + ":" + pip["BP"].astype(str) + ":" + pip["A1"] + ":" + pip["A2"]
    k2 = pip["CHR"].astype(str) + ":" + pip["BP"].astype(str) + ":" + pip["A2"] + ":" + pip["A1"]
    map1 = dict(zip(k1, pip["interactionQTL"]))
    map2 = dict(zip(k2, pip["interactionQTL"]))
    bim_keys = bim["CHR"].astype(str) + ":" + bim["BP"].astype(str) + ":" + bim["A1"].str.upper() + ":" + bim["A2"].str.upper()
    bim["interactionQTL"] = bim_keys.map(map1).fillna(bim_keys.map(map2)).fillna(0.0)
    bim[["CHR", "BP", "SNP", "CM", "interactionQTL"]].to_csv(out_gz, sep="\t", index=False, compression="gzip", float_format="%.6g")

elif mode == "calcsd":
    pop, sldsc_ref, ldscore_dir, out_file = sys.argv[2:6]
    vals, n_common = [], 0
    for chrom in range(1, 23):
        frq = pd.read_csv(f"{sldsc_ref}/{pop}_hg38/plink_files/1000G.{pop}.hg38.{chrom}.frq", sep=r"\s+", usecols=["SNP", "MAF"])
        annot = pd.read_csv(f"{ldscore_dir}/{pop}/interaction_maxCPP.{chrom}.annot.gz", sep="\t")
        common = pd.merge(annot, frq, on="SNP")
        common = common[(common["MAF"] > 0.05) & (common["MAF"] <= 0.5)]
        vals.append(common["interactionQTL"].values)
        n_common += len(common)
    all_vals = np.concatenate(vals)
    pd.DataFrame([{"population": pop, "stage": "interaction", "annot_set": "interaction_maxCPP",
                   "annotation": "interactionQTL", "n_snps": n_common,
                   "mean": np.mean(all_vals), "sd": np.std(all_vals, ddof=0)}]).to_csv(out_file, sep="\t", index=False)
