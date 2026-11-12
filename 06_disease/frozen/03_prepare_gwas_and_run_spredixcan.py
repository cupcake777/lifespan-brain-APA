#!/usr/bin/env python3
# Align harmonized GWAS alleles to full-TWAS model variants and run local SPrediXcan.
# GWAS N is kept per variant; h2/N calibration is passed only when both options are given.
import argparse
import csv
import gzip
import os
import re
import sqlite3
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
TWAS_DIR = ROOT / "data" / "06_disease" / "twas"
MODEL_DIR = TWAS_DIR / "models"
GWAS_DIR = Path(os.environ["APA_GROUP_HOME"]) / "liyc/share_group_folder/.ref/GWAS/harmony"
DEFAULT_GWAS_ID = "EAS_SCZ"
METAXCAN = Path(os.environ["APA_HPC_HOME"]) / "bin/MetaXcan-master/software/SPrediXcan.py"
SPREDIXCAN_DIR = TWAS_DIR / "spredixcan" / "by_run"
QC_DIR = TWAS_DIR / "qc"
LOG_DIR = ROOT / "logs" / "spredixcan"
OUT_DIR = SPREDIXCAN_DIR
STAGES = ("prenatal", "postnatal")

# --- helpers ---


def stage_qc_dir(stage):
    return QC_DIR / "by_run" / stage


def stage_aligned_dir(stage):
    return QC_DIR / "gwas_aligned" / "by_run" / stage


def stage_log_dir(stage):
    return LOG_DIR / stage


# --- arguments ---

parser = argparse.ArgumentParser(description="Align GWAS to full-TWAS models and run local SPrediXcan.")
parser.add_argument("--stage", choices=STAGES, required=True)
parser.add_argument("--gwas-id", default=DEFAULT_GWAS_ID, help="GWAS identifier such as EAS_SCZ or EAS_AD")
parser.add_argument("--gwas-source", type=Path, default=None, help="Override harmonized GWAS source path")
parser.add_argument("--gwas-dir", type=Path, default=GWAS_DIR, help="Directory containing harmonized GWAS TSVs")
parser.add_argument("--model-dir", type=Path, default=MODEL_DIR)
parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
parser.add_argument("--format-gwas", action="store_true")
parser.add_argument("--run-spredixcan", action="store_true")
parser.add_argument("--calibration-gwas-N", type=int, default=None)
parser.add_argument("--calibration-gwas-h2", type=float, default=None)
parser.add_argument("--additional-output", action="store_true")
parser.add_argument(
    "--stream-covariance",
    action="store_true",
    help="Pass MetaXcan --stream_covariance. Disabled by default because local MetaXcan currently fails on this path with pandas to_numeric(errors='ignore').",
)
args = parser.parse_args()

stage = args.stage
gwas_id = args.gwas_id
MODEL_DIR = args.model_dir.resolve()
out_dir = args.out_dir.resolve()
if not re.fullmatch(r"[A-Za-z0-9_]+", gwas_id):
    raise ValueError(f"Unsafe --gwas-id value: {gwas_id}")
if args.gwas_source is None:
    matches = sorted(args.gwas_dir.glob(f"{gwas_id}_*_harmonized.tsv"))
    if len(matches) != 1:
        raise ValueError(f"Expected one harmonized GWAS file for {gwas_id} in {args.gwas_dir}, found {[str(p) for p in matches]}")
    gwas_source = matches[0]
else:
    gwas_source = args.gwas_source
if not args.format_gwas and not args.run_spredixcan:
    raise ValueError("choose --format-gwas, --run-spredixcan, or both")

# --- format gwas ---

if args.format_gwas:
    db_path = MODEL_DIR / stage / f"{stage}.full_twas.db"
    if not db_path.exists():
        raise FileNotFoundError(db_path)
    with sqlite3.connect(db_path) as con:
        weights = pd.read_sql_query("SELECT DISTINCT rsid, ref_allele, eff_allele FROM weights", con)
    targets = defaultdict(list)
    for _, row in weights.iterrows():
        chrom, pos, ref, alt = str(row["rsid"]).split(":")
        ref_allele = str(row["ref_allele"]).upper()
        eff_allele = str(row["eff_allele"]).upper()
        targets[(chrom.removeprefix("chr"), pos)].append(
            {
                "snp": str(row["rsid"]),
                "ref_allele": ref_allele,
                "eff_allele": eff_allele,
                "allele_pair": frozenset([ref_allele, eff_allele]),
            }
        )

    aligned_dir = stage_aligned_dir(stage)
    audit_dir = stage_qc_dir(stage)
    aligned_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    formatted_path = aligned_dir / f"{gwas_id}.model_snps.allele_aligned.tsv"
    audit_path = audit_dir / f"{gwas_id}.gwas_alignment_audit.tsv"

    rows_written = 0
    rows_seen_at_model_positions = 0
    exact = 0
    flipped = 0
    allele_mismatch = 0
    n_values = set()
    matched_snps = set()

    with gwas_source.open() as src, formatted_path.open("w", newline="") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        required = {"CHR", "POS", "EA", "NEA", "BETA", "SE", "P", "N"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{gwas_source}: missing required columns {sorted(missing)}")
        fieldnames = ["SNP", "EA", "NEA", "BETA", "SE", "P", "N"]
        writer = csv.DictWriter(dst, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            key = (row["CHR"].removeprefix("chr"), row["POS"])
            candidates = targets.get(key)
            if not candidates:
                continue
            rows_seen_at_model_positions += 1
            ea = row["EA"].upper()
            nea = row["NEA"].upper()
            for candidate in candidates:
                if frozenset([ea, nea]) != candidate["allele_pair"]:
                    allele_mismatch += 1
                    continue
                beta = float(row["BETA"])
                if ea == candidate["eff_allele"] and nea == candidate["ref_allele"]:
                    aligned_beta = beta
                    exact += 1
                elif ea == candidate["ref_allele"] and nea == candidate["eff_allele"]:
                    aligned_beta = -beta
                    flipped += 1
                else:
                    allele_mismatch += 1
                    continue
                writer.writerow(
                    {
                        "SNP": candidate["snp"],
                        "EA": candidate["eff_allele"],
                        "NEA": candidate["ref_allele"],
                        "BETA": aligned_beta,
                        "SE": row["SE"],
                        "P": row["P"],
                        "N": row["N"],
                    }
                )
                rows_written += 1
                n_values.add(row["N"])
                matched_snps.add(candidate["snp"])

    audit = {
        "gwas_id": gwas_id,
        "stage": stage,
        "gwas_source": str(gwas_source),
        "model_variant_count": sum(len(v) for v in targets.values()),
        "matched_model_variant_count": len(matched_snps),
        "gwas_rows_seen_at_model_positions": rows_seen_at_model_positions,
        "formatted_gwas_rows": rows_written,
        "exact_effect_allele_matches": exact,
        "flipped_effect_allele_matches": flipped,
        "allele_mismatch_events": allele_mismatch,
        "unique_N_values": len(n_values),
        "min_N": min(map(float, n_values)) if n_values else None,
        "max_N": max(map(float, n_values)) if n_values else None,
        "formatted_gwas_path": str(formatted_path.relative_to(ROOT)),
    }
    pd.DataFrame([audit]).to_csv(audit_path, sep="\t", index=False, lineterminator="\n")
    print(pd.DataFrame([audit]).to_string(index=False))

# --- run spredixcan ---

if args.run_spredixcan:
    db_path = MODEL_DIR / stage / f"{stage}.full_twas.db"
    cov_path = MODEL_DIR / stage / f"{stage}.full_twas.cov"
    gwas_path = stage_aligned_dir(stage) / f"{gwas_id}.model_snps.allele_aligned.tsv"
    stage_out_dir = out_dir / stage
    stage_out_dir.mkdir(parents=True, exist_ok=True)
    stage_log_dir(stage).mkdir(parents=True, exist_ok=True)
    output_path = stage_out_dir / f"{gwas_id}.csv"
    compressed_output_path = stage_out_dir / f"{gwas_id}.csv.gz"
    log_path = stage_log_dir(stage) / f"{gwas_id}.spredixcan.log"

    cmd = [
        sys.executable,
        str(METAXCAN),
        "--model_db_path",
        str(db_path),
        "--gwas_file",
        str(gwas_path),
        "--snp_column",
        "SNP",
        "--effect_allele_column",
        "EA",
        "--non_effect_allele_column",
        "NEA",
        "--beta_column",
        "BETA",
        "--se_column",
        "SE",
        "--pvalue_column",
        "P",
        "--separator",
        "\t",
        "--keep_non_rsid",
        "--covariance",
        str(cov_path),
        "--output_file",
        str(output_path),
        "--overwrite",
        "--throw",
        "--verbosity",
        "8",
    ]
    if args.additional_output:
        cmd.append("--additional_output")
    if args.stream_covariance:
        cmd.append("--stream_covariance")
    calibration_gwas_n = args.calibration_gwas_N
    calibration_gwas_h2 = args.calibration_gwas_h2
    if calibration_gwas_n is not None or calibration_gwas_h2 is not None:
        if calibration_gwas_n is None or calibration_gwas_h2 is None:
            raise ValueError("calibration requires both --calibration-gwas-N and --calibration-gwas-h2")
        cmd.extend(["--gwas_N", str(calibration_gwas_n), "--gwas_h2", str(calibration_gwas_h2)])

    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit {result.returncode}: {' '.join(cmd)}\n{result.stdout}")
    rows = pd.read_csv(output_path).shape[0]
    with output_path.open("rb") as src, gzip.open(compressed_output_path, "wb") as dst:
        dst.writelines(src)
    output_path.unlink()
    summary = {
        "gwas_id": gwas_id,
        "stage": stage,
        "spredixcan_output_rows": rows,
        "spredixcan_output": str(compressed_output_path.relative_to(ROOT)),
        "spredixcan_log": str(log_path.relative_to(ROOT)),
        "calibration_used": int(calibration_gwas_n is not None and calibration_gwas_h2 is not None),
        "stream_covariance": int(args.stream_covariance),
    }
    summary_path = stage_qc_dir(stage) / f"{gwas_id}.spredixcan_run_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(summary_path, sep="\t", index=False, lineterminator="\n")
    print(pd.DataFrame([summary]).to_string(index=False))
