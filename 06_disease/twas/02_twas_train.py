#!/usr/bin/env python3
# Train glmnet TWAS models on PDUI 3UTR features of one stage, as a full stage or as one chunk.
# Writes weights, covariance, training detail, the stage DB and summary tables.

from __future__ import annotations

import argparse
import math
import os
import shutil
import sqlite3
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---
APA = Path(os.environ["APA_ROOT"])
INPUT_DIR = APA / "data" / "06_disease" / "twas" / "inputs"
OUT_DIR = APA / "data" / "06_disease" / "twas" / "models"
GENOTYPE_ROOT = APA / "data" / "04_qtl" / "pre_input"
GLMNET_WORKER = Path(__file__).resolve().parent / "glmnet_worker.R"
STAGES = ("prenatal", "postnatal")

DEFAULT_CHUNK_SIZE = 200
DEFAULT_THREADS = 1
DEFAULT_N_FOLDS = 5
RANDOM_SEED = 20260520
MIN_SNPS_FOR_MODEL = 2


# --- helpers ---
@dataclass
class StagePaths:
    annotation: Path
    phenotype: Path
    pfile: Path
    stage_dir: Path
    chunk_dir: Path
    tmp_dir: Path


def median_or_nan(values: list[float]) -> float:
    finite = [float(value) for value in values if np.isfinite(value)]
    return float(np.median(finite)) if finite else np.nan


def stage_paths(stage: str, chunk_id: int | None = None) -> StagePaths:
    stage_dir = OUT_DIR / stage
    chunk_name = f"chunk_{chunk_id:04d}" if chunk_id is not None else "combine"
    chunk_dir = stage_dir / "chunks" / chunk_name
    return StagePaths(
        annotation=INPUT_DIR / f"{stage}.annotation.tsv",
        phenotype=INPUT_DIR / f"{stage}.residualized_phenotype.tsv",
        pfile=GENOTYPE_ROOT / f"genotype.{stage}",
        stage_dir=stage_dir,
        chunk_dir=chunk_dir,
        tmp_dir=chunk_dir / "tmp_plink2_raw",
    )


def load_stage_inputs(stage: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = stage_paths(stage)
    annotation = pd.read_table(paths.annotation, dtype={"model_id": str})
    phenotype = pd.read_table(paths.phenotype, dtype={"model_id": str})
    if annotation.empty or phenotype.empty:
        raise ValueError(f"{stage}: full inputs are empty")
    if set(annotation["model_id"]) != set(phenotype["model_id"]):
        raise ValueError(f"{stage}: annotation and residualized phenotype model_id sets differ")
    return annotation, phenotype


def chunk_bounds(n_items: int, chunk_id: int, chunk_size: int) -> tuple[int, int]:
    if chunk_id < 1:
        raise ValueError("--chunk-id is 1-based")
    start = (chunk_id - 1) * chunk_size
    end = min(chunk_id * chunk_size, n_items)
    if start >= n_items:
        raise ValueError(f"chunk_id {chunk_id} starts after {n_items} features")
    return start, end


def p_adjust_bh(pvalues: np.ndarray) -> np.ndarray:
    adjusted = np.full(len(pvalues), np.nan)
    finite = np.isfinite(pvalues)
    if finite.sum() == 0:
        return adjusted
    order = np.argsort(pvalues[finite])
    finite_values = pvalues[finite][order]
    ranks = np.arange(1, len(finite_values) + 1)
    bh = finite_values * len(finite_values) / ranks
    bh = np.minimum.accumulate(bh[::-1])[::-1]
    bh = np.minimum(bh, 1.0)
    finite_indices = np.where(finite)[0][order]
    adjusted[finite_indices] = bh
    return adjusted


def train_selected_features(
    stage: str,
    selected: pd.DataFrame,
    phenotype: pd.DataFrame,
    paths: StagePaths,
    sample_ids: list[str],
    threads: int,
    n_folds: int,
    chunk_id: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, object]]:
    weights_all: list[pd.DataFrame] = []
    extra_all: list[dict[str, object]] = []
    cov_rows_all: list[dict[str, object]] = []
    detail_all: list[pd.DataFrame] = []
    plink_seconds: list[float] = []
    r_batch_seconds: list[float] = []
    unique_windows = 0

    phenotype_reset = phenotype.reset_index().copy()
    cached_key: tuple[str, int, int] | None = None
    cached_dosage: pd.DataFrame | None = None

    for key, group in selected.groupby(["chr", "training_left_1mb", "training_right_1mb"], sort=False):
        chrom, left, right = str(key[0]), int(key[1]), int(key[2])
        window = f"{chrom}:{left}-{right}"

        # --- export one cis window per group change ---
        if key != cached_key:
            started = time.monotonic()
            numeric_chr = chrom.removeprefix("chr")
            prefix = paths.tmp_dir / f"{stage}.{chrom}_{left}_{right}"
            raw_path = Path(f"{prefix}.raw")
            cmd = [
                "plink2",
                "--pfile",
                str(paths.pfile),
                "--chr",
                numeric_chr,
                "--from-bp",
                str(left),
                "--to-bp",
                str(right),
                "--snps-only",
                "just-acgt",
                "--export",
                "A",
                "--out",
                str(prefix),
                "--threads",
                str(threads),
            ]
            result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
            log_path = Path(f"{prefix}.plink2.log")
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(result.stdout)
            if result.returncode != 0:
                if "No variants remaining after main filters" not in result.stdout:
                    raise RuntimeError(f"Command failed with exit {result.returncode}: {' '.join(cmd)}\n{result.stdout}")
                cached_dosage = pd.DataFrame({"sample_id": sample_ids})
            else:
                cached_dosage = pd.read_table(raw_path, sep=r"\s+", dtype={"IID": str})
                cached_dosage = cached_dosage.rename(columns={"IID": "sample_id"})
            plink_seconds.append(time.monotonic() - started)
            cached_key = key
            unique_windows += 1

        # --- usable genotype matrix for this feature set ---
        dosage = cached_dosage
        variant_cols = list(dosage.columns[6:])
        sample_order_note = "exact"
        if not variant_cols:
            x = pd.DataFrame(index=sample_ids)
            metadata: list[dict[str, str]] = []
            dosage_info = {
                "sample_order": sample_order_note,
                "missing_genotype_values_imputed": 0,
                "missing_genotype_handling": "no usable genotype variants in cis window",
            }
        else:
            dosage = dosage.set_index("sample_id")
            missing_samples = [sample_id for sample_id in sample_ids if sample_id not in dosage.index]
            if missing_samples:
                raise ValueError(f"missing genotype samples: {missing_samples[:5]}")
            if list(dosage.index[: len(sample_ids)]) != sample_ids:
                sample_order_note = "reordered_to_phenotype"
            x = dosage.loc[sample_ids, variant_cols].apply(pd.to_numeric, errors="coerce")

            keep_cols: list[str] = []
            metadata = []
            for col in variant_cols:
                variant_id, counted_allele = col.rsplit("_", 1)
                counted_allele = counted_allele.upper()
                var_chrom, pos, ref, alt = variant_id.split(":")
                pos = int(pos)
                ref = ref.upper()
                alt = alt.upper()
                values = x[col]
                if len(ref) != 1 or len(alt) != 1 or len({ref, alt}) != 2:
                    continue
                if counted_allele not in {ref, alt}:
                    continue
                if values.notna().sum() < 3:
                    continue
                if np.nanvar(values.to_numpy(dtype=float)) <= 0:
                    continue
                keep_cols.append(col)
                metadata.append(
                    {
                        "raw_col": col,
                        "rsid": variant_id,
                        "chr": var_chrom,
                        "pos": str(pos),
                        "ref_allele": alt if counted_allele == ref else ref,
                        "eff_allele": counted_allele,
                    }
                )
            x = x[keep_cols].copy()
            missing_before = int(x.isna().sum().sum()) if not x.empty else 0
            for col in x.columns:
                if x[col].isna().any():
                    x[col] = x[col].fillna(x[col].median())
            dosage_info = {
                "sample_order": sample_order_note,
                "missing_genotype_values_imputed": missing_before,
                "missing_genotype_handling": "median imputation per SNP after deterministic phenotype sample ordering",
            }

        pheno_block = group[["model_id", "feature_index"]].merge(
            phenotype_reset, on="model_id", how="left", validate="one_to_one"
        )
        if pheno_block[sample_ids].isna().all(axis=1).any():
            raise ValueError(f"{stage} {window}: missing phenotype row after merge")

        if x.shape[1] < MIN_SNPS_FOR_MODEL:
            detail = group[["model_id", "feature_index", "original_phenotype_id"]].copy()
            detail["status"] = "too_few_snps"
            detail["n_training_samples"] = pheno_block[sample_ids].notna().sum(axis=1).to_numpy()
            detail["n_snps_available"] = x.shape[1]
            detail["n_nonzero_weights"] = 0
            detail["pred_perf_r"] = np.nan
            detail["pred_perf_R2"] = np.nan
            detail["pred_perf_pval"] = np.nan
            detail["lambda_min"] = np.nan
            detail["lambda_min_idx"] = np.nan
            detail["cv_complete_observations"] = np.nan
            detail["fit_preval_na"] = np.nan
            detail["glmnet_seconds"] = np.nan
            detail["pmax_mode"] = "control=list(pmax=200)"
            detail["glmnet_warnings"] = ""
        else:
            # --- glmnet fits for this window in one R batch ---
            worker_dir = paths.tmp_dir / f"glmnet_{chrom}_{left}_{right}"
            worker_dir.mkdir(parents=True, exist_ok=True)
            x_path = worker_dir / "glmnet_x.tsv"
            metadata_path = worker_dir / "glmnet_metadata.tsv"
            phenotype_path = worker_dir / "glmnet_phenotype.tsv"
            weights_path = worker_dir / "glmnet_weights.tsv"
            detail_path = worker_dir / "glmnet_detail.tsv"
            x_out = x.copy()
            x_out.insert(0, "sample_id", sample_ids)
            x_out.to_csv(x_path, sep="\t", index=False, lineterminator="\n")
            pd.DataFrame(metadata).to_csv(metadata_path, sep="\t", index=False, lineterminator="\n")
            pheno_block[["model_id", "feature_index", *sample_ids]].to_csv(
                phenotype_path, sep="\t", index=False, lineterminator="\n"
            )

            started = time.monotonic()
            cmd = [
                "Rscript",
                str(GLMNET_WORKER),
                str(x_path),
                str(metadata_path),
                str(phenotype_path),
                str(weights_path),
                str(detail_path),
                str(n_folds),
                str(RANDOM_SEED),
            ]
            result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
            (worker_dir / "glmnet_worker.log").write_text(result.stdout)
            if result.returncode != 0:
                raise RuntimeError(f"Command failed with exit {result.returncode}: {' '.join(cmd)}\n{result.stdout}")
            r_batch_seconds.append(time.monotonic() - started)
            weights = pd.read_table(weights_path, dtype={"model_id": str, "raw_col": str, "rsid": str, "gene": str})
            detail = pd.read_table(detail_path, dtype={"model_id": str})

            if not weights.empty:
                weights_all.append(weights[["model_id", "raw_col", "rsid", "gene", "weight", "ref_allele", "eff_allele"]])
                raw_to_rsid = {item["raw_col"]: item["rsid"] for item in metadata}
                for model_id, model_weights in weights.groupby("model_id", sort=False):
                    raw_cols = model_weights["raw_col"].tolist()
                    x_used = x[raw_cols]
                    cov = np.cov(x_used.to_numpy(dtype=float), rowvar=False)
                    if x_used.shape[1] == 1:
                        cov = np.array([[float(cov)]])
                    cols = list(x_used.columns)
                    for i, col_i in enumerate(cols):
                        for j in range(i, len(cols)):
                            col_j = cols[j]
                            cov_rows_all.append(
                                {
                                    "GENE": str(model_id),
                                    "RSID1": raw_to_rsid[col_i],
                                    "RSID2": raw_to_rsid[col_j],
                                    "VALUE": round(float(cov[i, j]), 3),
                                }
                            )
            detail = detail.merge(group[["model_id", "original_phenotype_id"]], on="model_id", how="left", validate="one_to_one")

        detail["stage"] = stage
        if chunk_id is not None:
            detail["chunk_id"] = chunk_id
        detail["window"] = window
        detail["sample_order"] = dosage_info["sample_order"]
        detail["missing_genotype_values_imputed"] = dosage_info["missing_genotype_values_imputed"]
        detail["missing_genotype_handling"] = dosage_info["missing_genotype_handling"]
        detail_all.append(detail)

    weights_df = pd.concat(weights_all, ignore_index=True) if weights_all else pd.DataFrame(columns=["model_id", "raw_col", "rsid", "gene", "weight", "ref_allele", "eff_allele"])
    detail_df = pd.concat(detail_all, ignore_index=True) if detail_all else pd.DataFrame()
    anno_by_model = selected.set_index("model_id")
    for model_id, model_weights in weights_df.groupby("model_id", sort=False):
        drow = detail_df.loc[detail_df["model_id"] == model_id].iloc[0]
        extra_all.append(
            {
                "gene": model_id,
                "genename": anno_by_model.loc[model_id, "gene_id_or_name"],
                "n.snps.in.model": int(len(model_weights)),
                "pred.perf.R2": drow["pred_perf_R2"],
                "pred.perf.pval": drow["pred_perf_pval"],
            }
        )
    extra_df = pd.DataFrame(extra_all, columns=["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval"])
    if not extra_df.empty:
        extra_df["pred.perf.qval"] = p_adjust_bh(extra_df["pred.perf.pval"].to_numpy(dtype=float))
    else:
        extra_df["pred.perf.qval"] = pd.Series(dtype=float)
    cov_df = pd.DataFrame(cov_rows_all, columns=["GENE", "RSID1", "RSID2", "VALUE"])
    summary = {
        "features_attempted": int(len(selected)),
        "features_failed": int(detail_df["status"].astype(str).str.startswith("fit_failed").sum()) if not detail_df.empty else 0,
        "features_with_zero_usable_variants": int((detail_df["n_snps_available"].fillna(0).astype(float) < MIN_SNPS_FOR_MODEL).sum()) if not detail_df.empty else 0,
        "features_with_nonzero_weights": int((detail_df["n_nonzero_weights"].fillna(0).astype(float) > 0).sum()) if not detail_df.empty else 0,
        "features_passing_perf": int((detail_df["status"] == "written").sum()) if not detail_df.empty else 0,
        "features_written": int(extra_df["gene"].nunique()) if not extra_df.empty else 0,
        "total_nonzero_weights": int(len(weights_df)),
        "covariance_rows": int(len(cov_df)),
        "unique_export_windows": int(unique_windows),
        "median_plink_export_seconds": median_or_nan(plink_seconds),
        "median_glmnet_seconds_per_feature": median_or_nan(detail_df["glmnet_seconds"].tolist()) if "glmnet_seconds" in detail_df else np.nan,
        "median_r_batch_seconds_per_window": median_or_nan(r_batch_seconds),
    }
    return weights_df, extra_df, cov_df, detail_df, summary


def write_stage_db(
    stage: str,
    stage_dir: Path,
    weights_path: Path,
    extra_path: Path,
    cov_path: Path,
    detail_path: Path,
    sample_count: int,
    extra_rows: list[dict[str, object]],
    features_attempted: int,
    weight_rows: int,
    covariance_rows: int,
    unique_export_windows: int,
) -> dict[str, object]:
    extra = pd.DataFrame(extra_rows)
    if not extra.empty:
        extra["pred.perf.qval"] = p_adjust_bh(extra["pred.perf.pval"].to_numpy(dtype=float))
    else:
        extra = pd.DataFrame(columns=["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"])

    output_db = stage_dir / f"{stage}.full_twas.db"
    if output_db.exists():
        output_db.unlink()

    with sqlite3.connect(output_db) as con:
        wrote_weights = False
        if weight_rows > 0:
            for weights_chunk in pd.read_table(weights_path, chunksize=200_000):
                weights_chunk[["rsid", "gene", "weight", "ref_allele", "eff_allele"]].to_sql(
                    "weights",
                    con,
                    index=False,
                    if_exists="append" if wrote_weights else "replace",
                )
                wrote_weights = True
        if not wrote_weights:
            pd.DataFrame(columns=["rsid", "gene", "weight", "ref_allele", "eff_allele"]).to_sql("weights", con, index=False)
        extra[["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]].to_sql("extra", con, index=False)
        pd.DataFrame([{"n.samples": sample_count}]).to_sql("sample_info", con, index=False)

    extra.to_csv(stage_dir / f"{stage}.extra.tsv", sep="\t", index=False, lineterminator="\n")
    summary = {
        "stage": stage,
        "features_attempted": features_attempted,
        "models_in_db": int(extra["gene"].nunique()) if not extra.empty else 0,
        "weight_rows": weight_rows,
        "covariance_rows": covariance_rows,
        "unique_export_windows": unique_export_windows,
        "sample_count": sample_count,
        "db_path": str(output_db.relative_to(APA)),
        "cov_path": str(cov_path.relative_to(APA)),
        "detail_path": str(detail_path.relative_to(APA)),
    }
    pd.DataFrame([summary]).to_csv(stage_dir / "full_stage_summary.tsv", sep="\t", index=False, lineterminator="\n")
    print(pd.DataFrame([summary]).to_string(index=False))
    return summary


# --- cli ---
parser = argparse.ArgumentParser()
sub = parser.add_subparsers(dest="command", required=True)

train_stage_parser = sub.add_parser("train-stage", help="train one full stage and write stage-level DB/cov")
train_stage_parser.add_argument("--stage", choices=STAGES, required=True)
train_stage_parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
train_stage_parser.add_argument("--n-folds", type=int, default=DEFAULT_N_FOLDS)
train_stage_parser.add_argument("--clean-tmp", action="store_true")
train_stage_parser.add_argument("--out-dir", type=Path, default=OUT_DIR)

manifest = sub.add_parser("manifest", help="fallback only: write chunk manifest for one stage")
manifest.add_argument("--stage", choices=STAGES, required=True)
manifest.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)

train = sub.add_parser("train-chunk", help="fallback only: train one stage/chunk")
train.add_argument("--stage", choices=STAGES, required=True)
train.add_argument("--chunk-id", type=int, required=True)
train.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
train.add_argument("--threads", type=int, default=DEFAULT_THREADS)
train.add_argument("--n-folds", type=int, default=DEFAULT_N_FOLDS)
train.add_argument("--clean-tmp", action="store_true")
train.add_argument("--out-dir", type=Path, default=OUT_DIR, help="Override model output root for blocking validation runs.")

combine = sub.add_parser("combine", help="fallback only: combine completed chunks into final DB/cov")
combine.add_argument("--stage", choices=STAGES, required=True)
combine.add_argument("--out-dir", type=Path, default=OUT_DIR, help="Override model output root for blocking validation runs.")

subset = sub.add_parser("train-subset", help="train a bounded validation subset only")
subset.add_argument("--stage", choices=STAGES, required=True)
subset.add_argument("--n-features", type=int, required=True)
subset.add_argument("--feature-start", type=int, default=1)
subset.add_argument("--threads", type=int, default=DEFAULT_THREADS)
subset.add_argument("--n-folds", type=int, default=DEFAULT_N_FOLDS)
subset.add_argument("--out-dir", type=Path, required=True)
subset.add_argument("--keep-tmp", action="store_true")

args = parser.parse_args()
stage = args.stage

# --- manifest ---
if args.command == "manifest":
    annotation, _ = load_stage_inputs(stage)
    rows = []
    n_chunks = math.ceil(len(annotation) / args.chunk_size)
    for chunk_id in range(1, n_chunks + 1):
        start, end = chunk_bounds(len(annotation), chunk_id, args.chunk_size)
        rows.append(
            {
                "stage": stage,
                "chunk_id": chunk_id,
                "chunk_size": args.chunk_size,
                "feature_start_1based": start + 1,
                "feature_end_1based": end,
                "n_features": end - start,
            }
        )
    out = OUT_DIR / stage / "chunk_manifest.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, lineterminator="\n")
    print(pd.DataFrame(rows).to_string(index=False))

# --- train-stage ---
elif args.command == "train-stage":
    OUT_DIR = args.out_dir.resolve()
    annotation, phenotype = load_stage_inputs(stage)
    phenotype = phenotype.set_index("model_id")
    paths = stage_paths(stage)
    paths.stage_dir.mkdir(parents=True, exist_ok=True)
    paths.tmp_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = [col for col in phenotype.columns if col != "model_id"]
    output_db = paths.stage_dir / f"{stage}.full_twas.db"
    output_cov = paths.stage_dir / f"{stage}.full_twas.cov"
    output_weights = paths.stage_dir / f"{stage}.weights.tsv"
    output_detail = paths.stage_dir / f"{stage}.training_detail.tsv"
    for output in (output_db, output_cov, output_weights, output_detail):
        if output.exists():
            output.unlink()

    weights, extra, cov, detail, run_summary = train_selected_features(
        stage, annotation, phenotype, paths, sample_ids, args.threads, args.n_folds
    )
    weights[["rsid", "gene", "weight", "ref_allele", "eff_allele"]].to_csv(output_weights, sep="\t", index=False, lineterminator="\n")
    cov[["GENE", "RSID1", "RSID2", "VALUE"]].to_csv(output_cov, sep="\t", index=False, lineterminator="\n")
    detail.to_csv(output_detail, sep="\t", index=False, lineterminator="\n")

    summary = write_stage_db(
        stage,
        paths.stage_dir,
        output_weights,
        paths.stage_dir / f"{stage}.extra.tsv",
        output_cov,
        output_detail,
        len(sample_ids),
        extra.to_dict("records"),
        features_attempted=len(annotation),
        weight_rows=len(weights),
        covariance_rows=len(cov),
        unique_export_windows=run_summary["unique_export_windows"],
    )
    summary.update(run_summary)
    pd.DataFrame([summary]).to_csv(paths.stage_dir / "full_stage_summary.tsv", sep="\t", index=False, lineterminator="\n")
    if args.clean_tmp and paths.tmp_dir.exists():
        shutil.rmtree(paths.tmp_dir)

# --- train-chunk ---
elif args.command == "train-chunk":
    OUT_DIR = args.out_dir.resolve()
    chunk_id = args.chunk_id
    annotation, phenotype = load_stage_inputs(stage)
    start, end = chunk_bounds(len(annotation), chunk_id, args.chunk_size)
    selected = annotation.iloc[start:end].copy()
    phenotype = phenotype.set_index("model_id")
    paths = stage_paths(stage, chunk_id)
    paths.chunk_dir.mkdir(parents=True, exist_ok=True)
    paths.tmp_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = [col for col in phenotype.columns if col != "model_id"]
    weights_df, extra_df, cov_df, detail_df, run_summary = train_selected_features(
        stage, selected, phenotype, paths, sample_ids, args.threads, args.n_folds, chunk_id=chunk_id
    )

    weights_df[["rsid", "gene", "weight", "ref_allele", "eff_allele"]].to_csv(
        paths.chunk_dir / "weights.tsv", sep="\t", index=False, lineterminator="\n"
    )
    extra_df.to_csv(paths.chunk_dir / "extra.tsv", sep="\t", index=False, lineterminator="\n")
    cov_df[["GENE", "RSID1", "RSID2", "VALUE"]].to_csv(paths.chunk_dir / "covariance.tsv", sep="\t", index=False, lineterminator="\n")
    detail_df.to_csv(paths.chunk_dir / "training_detail.tsv", sep="\t", index=False, lineterminator="\n")

    summary = {
        "stage": stage,
        "chunk_id": chunk_id,
        "chunk_dir": str(paths.chunk_dir.relative_to(APA)),
        **run_summary,
    }
    pd.DataFrame([summary]).to_csv(paths.chunk_dir / "chunk_summary.tsv", sep="\t", index=False, lineterminator="\n")
    if args.clean_tmp and paths.tmp_dir.exists():
        shutil.rmtree(paths.tmp_dir)
    print(pd.DataFrame([summary]).to_string(index=False))

# --- combine ---
elif args.command == "combine":
    OUT_DIR = args.out_dir.resolve()
    paths = stage_paths(stage)
    stage_dir = paths.stage_dir
    chunk_dirs = sorted((stage_dir / "chunks").glob("chunk_*"))
    if not chunk_dirs:
        raise ValueError(f"{stage}: no chunk outputs under {stage_dir / 'chunks'}")

    chunk_summary_parts = [pd.read_table(path / "chunk_summary.tsv") for path in chunk_dirs if (path / "chunk_summary.tsv").exists()]
    expected_chunks = pd.read_table(stage_dir / "chunk_manifest.tsv") if (stage_dir / "chunk_manifest.tsv").exists() else pd.DataFrame()
    if not expected_chunks.empty and len(chunk_summary_parts) != len(expected_chunks):
        raise ValueError(f"{stage}: chunk summaries present for {len(chunk_summary_parts)} of {len(expected_chunks)} expected chunks")

    weights_parts = [pd.read_table(path / "weights.tsv") for path in chunk_dirs if (path / "weights.tsv").exists()]
    extra_parts = [pd.read_table(path / "extra.tsv") for path in chunk_dirs if (path / "extra.tsv").exists()]
    cov_parts = [pd.read_table(path / "covariance.tsv") for path in chunk_dirs if (path / "covariance.tsv").exists()]
    detail_parts = [pd.read_table(path / "training_detail.tsv") for path in chunk_dirs if (path / "training_detail.tsv").exists()]
    if not weights_parts or not extra_parts or not cov_parts:
        raise ValueError(f"{stage}: missing weights/extra/covariance chunk outputs")

    weights = pd.concat(weights_parts, ignore_index=True)
    extra = pd.concat(extra_parts, ignore_index=True)
    cov = pd.concat(cov_parts, ignore_index=True)
    detail = pd.concat(detail_parts, ignore_index=True) if detail_parts else pd.DataFrame()

    if not extra.empty:
        extra["pred.perf.qval"] = p_adjust_bh(extra["pred.perf.pval"].to_numpy(dtype=float))
    sample_count = len(pd.read_table(INPUT_DIR / f"{stage}.residualized_phenotype.tsv", nrows=0).columns) - 1

    output_db = stage_dir / f"{stage}.full_twas.db"
    output_cov = stage_dir / f"{stage}.full_twas.cov"
    output_detail = stage_dir / f"{stage}.training_detail.tsv"
    for output in (output_db, output_cov):
        if output.exists():
            output.unlink()

    with sqlite3.connect(output_db) as con:
        if weights.empty:
            pd.DataFrame(columns=["rsid", "gene", "weight", "ref_allele", "eff_allele"]).to_sql("weights", con, index=False)
        else:
            weights[["rsid", "gene", "weight", "ref_allele", "eff_allele"]].to_sql("weights", con, index=False)
        if extra.empty:
            pd.DataFrame(columns=["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]).to_sql("extra", con, index=False)
        else:
            extra[["gene", "genename", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval"]].to_sql("extra", con, index=False)
        pd.DataFrame([{"n.samples": sample_count}]).to_sql("sample_info", con, index=False)
    cov[["GENE", "RSID1", "RSID2", "VALUE"]].to_csv(output_cov, sep="\t", index=False, lineterminator="\n")
    detail.to_csv(output_detail, sep="\t", index=False, lineterminator="\n")

    chunk_summary = pd.concat(chunk_summary_parts, ignore_index=True) if chunk_summary_parts else pd.DataFrame()
    summary = {
        "stage": stage,
        "chunks_combined": len(chunk_dirs),
        "models_in_db": int(extra["gene"].nunique()) if not extra.empty else 0,
        "weight_rows": len(weights),
        "covariance_rows": len(cov),
        "sample_count": sample_count,
        "features_attempted": int(chunk_summary["features_attempted"].sum()) if "features_attempted" in chunk_summary else np.nan,
        "features_failed": int(chunk_summary["features_failed"].sum()) if "features_failed" in chunk_summary else np.nan,
        "features_with_zero_usable_variants": int(chunk_summary["features_with_zero_usable_variants"].sum()) if "features_with_zero_usable_variants" in chunk_summary else np.nan,
        "features_with_nonzero_weights": int(chunk_summary["features_with_nonzero_weights"].sum()) if "features_with_nonzero_weights" in chunk_summary else np.nan,
        "features_passing_perf": int(chunk_summary["features_passing_perf"].sum()) if "features_passing_perf" in chunk_summary else np.nan,
        "features_written": int(chunk_summary["features_written"].sum()) if "features_written" in chunk_summary else np.nan,
        "median_plink_export_seconds": median_or_nan(chunk_summary["median_plink_export_seconds"].tolist()) if "median_plink_export_seconds" in chunk_summary else np.nan,
        "median_glmnet_seconds_per_feature": median_or_nan(chunk_summary["median_glmnet_seconds_per_feature"].tolist()) if "median_glmnet_seconds_per_feature" in chunk_summary else np.nan,
        "db_path": str(output_db.relative_to(APA)),
        "cov_path": str(output_cov.relative_to(APA)),
        "detail_path": str(output_detail.relative_to(APA)),
    }
    pd.DataFrame([summary]).to_csv(stage_dir / "combine_summary.tsv", sep="\t", index=False, lineterminator="\n")
    print(pd.DataFrame([summary]).to_string(index=False))

# --- train-subset ---
elif args.command == "train-subset":
    OUT_DIR = args.out_dir.resolve()
    annotation, phenotype = load_stage_inputs(stage)
    start = args.feature_start - 1
    end = min(start + args.n_features, len(annotation))
    if start < 0 or start >= len(annotation):
        raise ValueError(f"{stage}: feature_start {args.feature_start} outside 1..{len(annotation)}")
    selected = annotation.iloc[start:end].copy()
    phenotype = phenotype.set_index("model_id")

    paths = stage_paths(stage, chunk_id=1)
    paths.chunk_dir.mkdir(parents=True, exist_ok=True)
    paths.tmp_dir.mkdir(parents=True, exist_ok=True)
    sample_ids = [col for col in phenotype.columns if col != "model_id"]
    weights, extra, cov, detail, run_summary = train_selected_features(
        stage, selected, phenotype, paths, sample_ids, args.threads, args.n_folds, chunk_id=1
    )
    stage_dir = OUT_DIR / stage
    stage_dir.mkdir(parents=True, exist_ok=True)
    weights_path = paths.chunk_dir / "weights.tsv"
    extra_path = paths.chunk_dir / "extra.tsv"
    cov_path = paths.chunk_dir / "covariance.tsv"
    detail_path = paths.chunk_dir / "training_detail.tsv"
    stage_cov_path = stage_dir / f"{stage}.full_twas.cov"
    stage_detail_path = stage_dir / f"{stage}.training_detail.tsv"
    weights[["rsid", "gene", "weight", "ref_allele", "eff_allele"]].to_csv(weights_path, sep="\t", index=False, lineterminator="\n")
    extra.to_csv(extra_path, sep="\t", index=False, lineterminator="\n")
    cov[["GENE", "RSID1", "RSID2", "VALUE"]].to_csv(cov_path, sep="\t", index=False, lineterminator="\n")
    detail.to_csv(detail_path, sep="\t", index=False, lineterminator="\n")
    cov[["GENE", "RSID1", "RSID2", "VALUE"]].to_csv(stage_cov_path, sep="\t", index=False, lineterminator="\n")
    detail.to_csv(stage_detail_path, sep="\t", index=False, lineterminator="\n")
    summary = {"stage": stage, "feature_start": args.feature_start, "n_features": len(selected), **run_summary}
    pd.DataFrame([summary]).to_csv(paths.chunk_dir / "chunk_summary.tsv", sep="\t", index=False, lineterminator="\n")
    write_stage_db(
        stage,
        stage_dir,
        weights_path,
        extra_path,
        stage_cov_path,
        stage_detail_path,
        len(sample_ids),
        extra.to_dict("records"),
        features_attempted=len(selected),
        weight_rows=len(weights),
        covariance_rows=len(cov),
        unique_export_windows=run_summary["unique_export_windows"],
    )
    pd.DataFrame([summary]).to_csv(stage_dir / "validation_summary.tsv", sep="\t", index=False, lineterminator="\n")
    if not args.keep_tmp and paths.tmp_dir.exists():
        shutil.rmtree(paths.tmp_dir)
    print(pd.DataFrame([summary]).to_string(index=False))
