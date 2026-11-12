#!/usr/bin/env python3
import os

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT_DIR = Path(os.environ["APA_ROOT"])
RAW = ROOT_DIR / "data" / "07_coloc" / "intx.signif_pairs.txt.gz"
OUT = ROOT_DIR / "data" / "07_coloc" / "interaction.susie_qtl.standardized.tsv.gz"


def parse_phenotype_fields(df: pd.DataFrame) -> pd.DataFrame:
    phenotype_parts = df["phenotype_id"].str.split("|", expand=True, regex=False)
    if phenotype_parts.shape[1] != 4:
        raise ValueError("phenotype_id does not have four pipe-delimited fields")
    phenotype_parts.columns = ["transcript_id", "gene_name", "phenotype_region", "phenotype_strand"]
    region_parts = phenotype_parts["phenotype_region"].str.extract(
        r"^(?P<phenotype_chr>[^:]+):(?P<phenotype_start>\d+)-(?P<phenotype_end>\d+)$"
    )
    if region_parts.isna().any().any():
        raise ValueError("failed to parse phenotype interval")
    region_parts["phenotype_start"] = region_parts["phenotype_start"].astype(int)
    region_parts["phenotype_end"] = region_parts["phenotype_end"].astype(int)
    return pd.concat([df.reset_index(drop=True), phenotype_parts.drop(columns="phenotype_region"), region_parts], axis=1)


def parse_variant_fields(df: pd.DataFrame) -> pd.DataFrame:
    variant_parts = df["variant_id"].str.extract(
        r"^(?P<variant_chr>[^:]+):(?P<variant_pos>\d+):(?P<ref_allele>[^:]+):(?P<alt_allele>[^:]+)$"
    )
    if variant_parts.isna().any().any():
        raise ValueError("failed to parse variant_id")
    variant_parts["variant_pos"] = variant_parts["variant_pos"].astype(int)
    variant_parts["variant_start0"] = variant_parts["variant_pos"] - 1
    variant_parts["variant_end"] = variant_parts["variant_pos"]
    return pd.concat([df.reset_index(drop=True), variant_parts], axis=1)


def main() -> None:
    df = pd.read_csv(RAW, sep="\t")
    df["stage"] = "interaction_GxStage"
    df = parse_variant_fields(parse_phenotype_fields(df))
    df["pip"] = pd.NA
    df["cs_id"] = "interaction_sig"
    if "slope" in df.columns and "b_g" not in df.columns:
        df["b_g"] = pd.to_numeric(df["slope"], errors="coerce")
    if "pval_nominal" in df.columns and "pval_gi" not in df.columns:
        df["pval_gi"] = pd.to_numeric(df["pval_nominal"], errors="coerce")
    if "slope_se" in df.columns and "b_g_se" not in df.columns:
        df["b_g_se"] = pd.to_numeric(df["slope_se"], errors="coerce")
    for col in ["b_g", "b_gi", "pval_gi", "b_g_se", "b_gi_se"]:
        if col not in df.columns:
            df[col] = pd.NA
    out_cols = [
        "stage",
        "phenotype_id",
        "transcript_id",
        "gene_name",
        "phenotype_chr",
        "phenotype_start",
        "phenotype_end",
        "phenotype_strand",
        "variant_id",
        "variant_chr",
        "variant_start0",
        "variant_pos",
        "variant_end",
        "ref_allele",
        "alt_allele",
        "pip",
        "af",
        "cs_id",
        "slope",
        "b_g",
        "b_gi",
        "pval_gi",
        "b_g_se",
        "b_gi_se",
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df[[col for col in out_cols if col in df.columns]].to_csv(OUT, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()
