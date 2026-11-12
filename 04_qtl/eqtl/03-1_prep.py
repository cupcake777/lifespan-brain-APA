#!/usr/bin/env python3
# S-LDSC inputs: maxCPP sites from SuSiE summaries and per-chromosome apaQTL/eQTL annotations.
# Usage: 03-1_prep.py extract <susie.parquet> <out.tsv> | annot <bim> <apa.tsv> <eqtl.tsv> <out.annot.gz>
import sys

import pandas as pd
import pyarrow.parquet as pq

# --- helpers ---
def add_keys(d, val_col):
    k1 = d["CHR"].astype(str) + ":" + d["BP"].astype(str) + ":" + d["A1"] + ":" + d["A2"]
    k2 = d["CHR"].astype(str) + ":" + d["BP"].astype(str) + ":" + d["A2"] + ":" + d["A1"]
    return dict(zip(k1, d[val_col])), dict(zip(k2, d[val_col]))

# --- mode ---
mode = sys.argv[1]

# --- extract ---
if mode == "extract":
    susie_parquet, out_tsv = sys.argv[2], sys.argv[3]
    tab = pq.read_table(susie_parquet, columns=["variant_id", "pip"])
    df = pd.DataFrame({"variant_id": tab["variant_id"].to_pylist(), "pip": tab["pip"].to_pylist()})
    df = df.dropna().groupby("variant_id", as_index=False)["pip"].max()
    splits = df["variant_id"].str.split(":", expand=True)
    df["chr"] = splits[0].str.replace("chr", "")
    df["pos"] = splits[1]
    df["a1"] = splits[2].str.upper()
    df["a2"] = splits[3].str.upper()
    df = df[df["chr"].str.isdigit()].copy()
    df["chr_num"] = df["chr"].astype(int)
    df["pos_num"] = df["pos"].astype(int)
    df = df.sort_values(["chr_num", "pos_num"])
    df[["chr", "pos", "a1", "a2", "pip"]].to_csv(out_tsv, sep="\t", index=False, header=False, float_format="%.6g")

# --- annot ---
elif mode == "annot":
    bim_file, apa_tsv, eqtl_tsv, out_gz = sys.argv[2:6]
    bim = pd.read_csv(bim_file, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    apa = pd.read_csv(apa_tsv, sep="\t", header=None, names=["CHR", "BP", "A1", "A2", "apaQTL"])
    eqtl = pd.read_csv(eqtl_tsv, sep="\t", header=None, names=["CHR", "BP", "A1", "A2", "eQTL"])

    apa_map1, apa_map2 = add_keys(apa, "apaQTL")
    eqtl_map1, eqtl_map2 = add_keys(eqtl, "eQTL")

    bim_keys = bim["CHR"].astype(str) + ":" + bim["BP"].astype(str) + ":" + bim["A1"].str.upper() + ":" + bim["A2"].str.upper()

    bim["apaQTL"] = bim_keys.map(apa_map1).fillna(bim_keys.map(apa_map2)).fillna(0.0)
    bim["eQTL"] = bim_keys.map(eqtl_map1).fillna(bim_keys.map(eqtl_map2)).fillna(0.0)

    bim[["CHR", "BP", "SNP", "CM", "apaQTL", "eQTL"]].to_csv(out_gz, sep="\t", index=False, compression="gzip", float_format="%.6g")
