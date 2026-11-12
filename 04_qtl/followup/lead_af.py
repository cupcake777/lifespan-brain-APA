#!/usr/bin/env python3
# ALT-allele frequency of apaQTL lead variants in 1000 Genomes EUR and EAS (GRCh38 PLINK reference), vs CHB cohort.
# Reads per-chromosome .bim and .frq files under ref/sLDSC; writes lead_af_eas_eur.tsv and lead_af_eas_eur_summary.tsv.
import os
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(os.environ["APA_ROOT"])
REF = Path(os.environ["APA_GROUP_SHARE"]) / ".ref" / "sLDSC"
POPS = {"EUR": REF / "EUR_hg38/plink_files/1000G.EUR.hg38", "EAS": REF / "EAS_hg38/plink_files/1000G.EAS.hg38"}
D = ROOT / "data" / "04_qtl"

# --- lead variants ---
lv = pd.read_csv(ROOT / "data" / "07_coloc" / "lead_variants_for_af.tsv", sep="\t")
parts = lv["variant_id"].str.split(":", expand=True)
lv["chr"] = parts[0].str.replace("chr", "", regex=False)
lv["pos"] = parts[1].astype(int)
lv["ref"], lv["alt"] = parts[2], parts[3]

# --- 1000G allele frequency, matched on position and allele pair ---
for pop, prefix in POPS.items():
    af = pd.Series(np.nan, index=lv.index)
    for chrom, d in lv.groupby("chr"):
        bim = pd.read_csv(Path(f"{prefix}.{chrom}.bim"), sep=r"\s+", header=None, names=["chr", "SNP", "cm", "bp", "b1", "b2"])
        frq = pd.read_csv(Path(f"{prefix}.{chrom}.frq"), sep=r"\s+")[["SNP", "A1", "A2", "MAF"]]
        m = bim[bim.bp.isin(set(d.pos))].merge(frq, on="SNP", how="inner")
        for i, row in d.iterrows():
            hit = m[(m.bp == row.pos) & (((m.A1 == row.alt) & (m.A2 == row.ref)) | ((m.A1 == row.ref) & (m.A2 == row.alt)))]
            if len(hit):
                h = hit.iloc[0]
                af[i] = h.MAF if h.A1 == row.alt else 1 - h.MAF
    lv[f"af_1KG_{pop}_alt"] = af
    print(f"[{pop}] matched {af.notna().sum()} / {len(lv)} lead rows", flush=True)

# --- MAF, deltas and flags ---
def maf(x): return np.minimum(x, 1 - x)
lv["maf_CHB"] = maf(lv.af_CHB_alt)
lv["maf_EUR"] = maf(lv.af_1KG_EUR_alt)
lv["delta_af_CHB_minus_EUR"] = lv.af_CHB_alt - lv.af_1KG_EUR_alt
lv["EAS_enriched"] = (lv.maf_CHB >= 0.05) & ((lv.maf_EUR < 0.01) | lv.af_1KG_EUR_alt.isna())
lv["large_af_difference"] = lv.delta_af_CHB_minus_EUR.abs() >= 0.2
lv.to_csv(D / "lead_af_eas_eur.tsv", sep="\t", index=False)

# --- summary per lead set ---
rows = []
for s, d in lv.groupby("set"):
    rows.append({"set": s, "n_leads": len(d), "in_1KG_EUR": int(d.af_1KG_EUR_alt.notna().sum()),
                 "EAS_enriched(MAF_CHB>=5%, MAF_EUR<1% or absent)": int(d.EAS_enriched.sum()),
                 "abs_dAF>=0.2": int(d.large_af_difference.sum()),
                 "r_CHB_vs_1KG_EAS": round(d[["af_CHB_alt", "af_1KG_EAS_alt"]].corr().iloc[0, 1], 3)})
pd.DataFrame(rows).to_csv(D / "lead_af_eas_eur_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
