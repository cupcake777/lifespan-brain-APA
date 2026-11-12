#!/usr/bin/env python
# Standardize SuSiE QTL tables (prenatal, postnatal, interaction) for colocalization and function steps.
# Inputs and outputs live in data/07_coloc/; reference table names are unchanged.
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

# --- paths and settings ---
ROOT = Path(os.environ["APA_ROOT"])
COLOC_DIR = ROOT / "data" / "07_coloc"
HIGH_PIP_THRESHOLD = 0.8
HIGH_PIP_TAG = f"{HIGH_PIP_THRESHOLD:.1f}"
STAGE_FILES = {
    "prenatal": "prenatal.SuSiE_summary.parquet",
    "postnatal": "postnatal.SuSiE_summary.parquet",
}
INTERACTION_SIG_QTL = COLOC_DIR / "intx.signif_pairs.txt.gz"
STAGE_ORDER = [*STAGE_FILES, "interaction_GxStage"]


# --- helpers ---
def parse_phenotype_fields(df: pd.DataFrame) -> pd.DataFrame:
    phenotype_parts = df["phenotype_id"].str.split("|", expand=True, regex=False)
    if phenotype_parts.shape[1] != 4:
        raise ValueError("phenotype_id does not have four pipe-delimited fields")
    phenotype_parts.columns = [
        "transcript_id",
        "gene_name",
        "phenotype_region",
        "phenotype_strand",
    ]
    region_parts = phenotype_parts["phenotype_region"].str.extract(
        r"^(?P<phenotype_chr>[^:]+):(?P<phenotype_start>\d+)-(?P<phenotype_end>\d+)$"
    )
    if region_parts.isna().any().any():
        raise ValueError("failed to parse phenotype genomic interval")
    region_parts["phenotype_start"] = region_parts["phenotype_start"].astype(int)
    region_parts["phenotype_end"] = region_parts["phenotype_end"].astype(int)
    phenotype_parts = pd.concat(
        [phenotype_parts.drop(columns="phenotype_region"), region_parts], axis=1
    )
    return pd.concat([df.reset_index(drop=True), phenotype_parts], axis=1)


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


def write_table(df: pd.DataFrame, path: Path) -> None:
    if path.suffix == ".gz":
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    else:
        df.to_csv(path, sep="\t", index=False)


def make_top_variant_table(df: pd.DataFrame) -> pd.DataFrame:
    ordered = df.sort_values(
        ["stage", "phenotype_id", "pip", "variant_chr", "variant_pos"],
        ascending=[True, True, False, True, True],
        kind="stable",
    )
    return ordered.groupby(["stage", "phenotype_id"], as_index=False, observed=True).head(1)


def make_phenotype_summary(df: pd.DataFrame) -> pd.DataFrame:
    summary = (
        df.groupby(
            [
                "stage",
                "phenotype_id",
                "transcript_id",
                "gene_name",
                "phenotype_chr",
                "phenotype_start",
                "phenotype_end",
                "phenotype_strand",
            ],
            as_index=False,
            observed=True,
        )
        .agg(
            n_variants=("variant_id", "size"),
            n_cs=("cs_id", "nunique"),
            max_pip=("pip", "max"),
            n_high_pip=("pip", lambda s: int((s >= HIGH_PIP_THRESHOLD).sum())),
        )
        .sort_values(
            ["stage", "max_pip", "gene_name", "phenotype_id"],
            ascending=[True, False, True, True],
            kind="stable",
        )
    )
    summary["phenotype_start0"] = summary["phenotype_start"] - 1
    return summary[
        [
            "stage",
            "phenotype_id",
            "transcript_id",
            "gene_name",
            "phenotype_chr",
            "phenotype_start0",
            "phenotype_start",
            "phenotype_end",
            "phenotype_strand",
            "n_variants",
            "n_cs",
            "max_pip",
            "n_high_pip",
        ]
    ]


def make_variant_bed(df: pd.DataFrame) -> pd.DataFrame:
    return df[
        [
            "variant_chr",
            "variant_start0",
            "variant_end",
            "variant_id",
            "phenotype_id",
            "gene_name",
            "stage",
            "pip",
            "cs_id",
        ]
    ].sort_values(
        ["variant_chr", "variant_start0", "variant_end", "phenotype_id"],
        ascending=[True, True, True, True],
        kind="stable",
    )


def make_phenotype_bed(df: pd.DataFrame) -> pd.DataFrame:
    return df[
        [
            "phenotype_chr",
            "phenotype_start0",
            "phenotype_end",
            "phenotype_id",
            "gene_name",
            "stage",
            "max_pip",
            "n_variants",
            "n_high_pip",
        ]
    ].sort_values(
        ["phenotype_chr", "phenotype_start0", "phenotype_end", "phenotype_id"],
        ascending=[True, True, True, True],
        kind="stable",
    )


# --- load stages ---
COLOC_DIR.mkdir(exist_ok=True)
stage_tables = {}
for stage, filename in STAGE_FILES.items():
    df = pd.read_parquet(COLOC_DIR / filename, engine="fastparquet").copy()
    df["stage"] = stage
    df = parse_phenotype_fields(df)
    df = parse_variant_fields(df)
    df["pip"] = df["pip"].astype(float)
    df["af"] = df["af"].astype(float)
    df["cs_id"] = df["cs_id"].astype(str)
    df["stage"] = pd.Categorical(df["stage"], categories=STAGE_ORDER, ordered=True)

    signif_path = COLOC_DIR / filename.replace("SuSiE_summary.parquet", "signif_pairs.txt.gz")
    if signif_path.exists():
        signif = pd.read_csv(signif_path, sep="\t", usecols=["phenotype_id", "variant_id", "slope"])
        df = df.merge(signif, on=["phenotype_id", "variant_id"], how="left")
    else:
        df["slope"] = float("nan")

    stage_tables[stage] = df[
        [
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
        ]
    ].sort_values(
        ["stage", "pip", "phenotype_id", "variant_chr", "variant_pos"],
        ascending=[True, False, True, True, True],
        kind="stable",
    )

if not INTERACTION_SIG_QTL.exists():
    interaction_df = pd.DataFrame(columns=stage_tables["prenatal"].columns)
else:
    df = pd.read_csv(INTERACTION_SIG_QTL, sep="\t").copy()
    df["stage"] = "interaction_GxStage"
    df = parse_phenotype_fields(df)
    df = parse_variant_fields(df)
    df["pip"] = pd.NA
    df["cs_id"] = "interaction_sig"
    if "slope" not in df.columns and "b_g" in df.columns:
        df["slope"] = df["b_g"]
    if "pval_nominal" in df.columns and "pval_gi" not in df.columns:
        df["pval_gi"] = pd.to_numeric(df["pval_nominal"], errors="coerce")
    if "slope_se" in df.columns and "b_g_se" not in df.columns:
        df["b_g_se"] = pd.to_numeric(df["slope_se"], errors="coerce")
    for col in ["af", "slope"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    interaction_df = df[
        [
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
        ]
    ].sort_values(
        ["stage", "phenotype_id", "variant_chr", "variant_pos"],
        ascending=[True, True, True, True],
        kind="stable",
    )
if not interaction_df.empty:
    stage_tables["interaction_GxStage"] = interaction_df
all_df = pd.concat(stage_tables.values(), ignore_index=True)

# --- per-stage tables ---
for stage, df in stage_tables.items():
    standardized_path = COLOC_DIR / f"{stage}.susie_qtl.standardized.tsv.gz"
    write_table(df, standardized_path)

    top_df = make_top_variant_table(df)
    write_table(top_df, COLOC_DIR / f"{stage}.top_variant_per_phenotype.tsv.gz")

    high_df = df.loc[df["pip"].fillna(1.0) >= HIGH_PIP_THRESHOLD].copy()
    write_table(high_df, COLOC_DIR / f"{stage}.high_pip.tsv.gz")
    write_table(high_df, COLOC_DIR / f"{stage}.high_pip_{HIGH_PIP_TAG}.tsv.gz")
    high_bed = make_variant_bed(high_df)
    write_table(high_bed, COLOC_DIR / f"{stage}.high_pip.variants.bed")
    write_table(high_bed, COLOC_DIR / f"{stage}.high_pip_{HIGH_PIP_TAG}.variants.bed")

    phenotype_summary = make_phenotype_summary(df)
    write_table(phenotype_summary, COLOC_DIR / f"{stage}.phenotype_summary.tsv.gz")
    write_table(make_phenotype_bed(phenotype_summary), COLOC_DIR / f"{stage}.phenotypes.bed")

# --- all-stage tables ---
write_table(all_df, COLOC_DIR / "all_stages.susie_qtl.standardized.tsv.gz")
write_table(make_top_variant_table(all_df), COLOC_DIR / "all_stages.top_variant_per_phenotype.tsv.gz")

all_high_df = all_df.loc[all_df["pip"].fillna(1.0) >= HIGH_PIP_THRESHOLD].copy()
write_table(all_high_df, COLOC_DIR / "all_stages.high_pip.tsv.gz")
write_table(all_high_df, COLOC_DIR / f"all_stages.high_pip_{HIGH_PIP_TAG}.tsv.gz")
all_high_bed = make_variant_bed(all_high_df)
write_table(all_high_bed, COLOC_DIR / "all_stages.high_pip.variants.bed")
write_table(all_high_bed, COLOC_DIR / f"all_stages.high_pip_{HIGH_PIP_TAG}.variants.bed")

all_phenotype_summary = make_phenotype_summary(all_df)
write_table(all_phenotype_summary, COLOC_DIR / "all_stages.phenotype_summary.tsv.gz")
write_table(make_phenotype_bed(all_phenotype_summary), COLOC_DIR / "all_stages.phenotypes.bed")

# --- manifests ---
rows = []
for stage, filename in STAGE_FILES.items():
    df = stage_tables[stage]
    rows.append(
        {
            "stage": stage,
            "summary_parquet": f"data/07_coloc/{filename}",
            "susie_pickle": f"data/07_coloc/{stage}.SuSiE.pickle",
            "tensorqtl_log": f"raw/{stage}.tensorQTL.cis_susie.log",
            "n_rows": int(df.shape[0]),
            "n_phenotypes": int(df["phenotype_id"].nunique()),
            "n_genes": int(df["gene_name"].nunique()),
            "n_high_pip": int((df["pip"].fillna(1.0) >= HIGH_PIP_THRESHOLD).sum()),
            "max_pip": float(pd.to_numeric(df["pip"], errors="coerce").max()) if df["pip"].notna().any() else None,
        }
    )
source_manifest = pd.DataFrame(rows)
source_manifest["stage"] = pd.Categorical(source_manifest["stage"], categories=STAGE_ORDER, ordered=True)
source_manifest = source_manifest.sort_values("stage", kind="stable")
write_table(source_manifest, COLOC_DIR / "qtl_source_manifest.tsv")

tested_universe = (
    all_df[
        [
            "phenotype_id",
            "transcript_id",
            "gene_name",
            "phenotype_chr",
            "phenotype_start",
            "phenotype_end",
            "phenotype_strand",
            "stage",
        ]
    ]
    .drop_duplicates()
    .copy()
)
tested_universe["source"] = "qtl_and_interaction_inputs"
tested_universe = tested_universe.sort_values(
    ["gene_name", "phenotype_id", "stage"],
    ascending=[True, True, True],
    kind="stable",
)
write_table(tested_universe, COLOC_DIR / "tested_apa_gene_universe.tsv")

summary_payload = {
    "high_pip_threshold": HIGH_PIP_THRESHOLD,
    "tested_apa_gene_universe": {
        "path": "data/07_coloc/tested_apa_gene_universe.tsv",
        "n_rows": int(tested_universe.shape[0]),
        "n_phenotypes": int(tested_universe["phenotype_id"].nunique()),
        "n_genes": int(tested_universe["gene_name"].nunique()),
        "source": "pipeline_derived_from_qtl_and_interaction_inputs",
    },
    "stages": {
        stage: {
            "n_rows": int(df.shape[0]),
            "n_phenotypes": int(df["phenotype_id"].nunique()),
            "n_genes": int(df["gene_name"].nunique()),
            "n_high_pip": int((df["pip"].fillna(1.0) >= HIGH_PIP_THRESHOLD).sum()),
            "n_top_variants": int(df["phenotype_id"].nunique()),
        }
        for stage, df in stage_tables.items()
    },
}
with (COLOC_DIR / "qtl_prepare_summary.json").open("w", encoding="utf-8") as handle:
    json.dump(summary_payload, handle, indent=2, ensure_ascii=True)
