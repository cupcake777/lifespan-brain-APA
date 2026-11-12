#!/usr/bin/env python3
# Prepare transcript-level 3UTR PDUI TWAS inputs (feature IDs kept, shell-safe model_id).
# Residualizes INT_PDUI on stage covariates; writes annotation, model_id map, residual and audit TSVs.
import csv
import gzip
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd

# --- config ---
APA = Path(os.environ["APA_ROOT"])
RAW_DIR = APA / "raw"
OUT_DIR = APA / "data" / "06_disease" / "twas" / "inputs"
CHUNK_SIZE = 1000

STAGES = {
    "prenatal": {
        "phenotype": RAW_DIR / "phenotype.Prenatal.GeneRange.bed.gz",
        "covariates": RAW_DIR / "covariates_for_qtl.Prenatal.txt",
        "psam": RAW_DIR / "genotype.prenatal.psam",
    },
    "postnatal": {
        "phenotype": RAW_DIR / "phenotype.Postnatal.GeneRange.bed.gz",
        "covariates": RAW_DIR / "covariates_for_qtl.Postnatal.txt",
        "psam": RAW_DIR / "genotype.postnatal.psam",
    },
}

ANNOTATION_COLUMNS = [
    "stage",
    "model_id",
    "feature_index",
    "original_phenotype_id",
    "transcript_id",
    "gene_id_or_name",
    "3utr_region",
    "strand",
    "chr",
    "group_left",
    "group_right",
    "training_left_1mb",
    "training_right_1mb",
    "transcript_left",
    "transcript_right",
    "group_interval_differs_from_transcript_interval",
    "n_complete_samples_for_residualization",
]

OUT_DIR.mkdir(parents=True, exist_ok=True)
audit_rows = []

for stage, cfg in STAGES.items():
    prefix = OUT_DIR / stage
    annotation_path = Path(f"{prefix}.annotation.tsv")
    map_path = Path(f"{prefix}.model_id_map.tsv")
    residual_path = Path(f"{prefix}.residualized_phenotype.tsv")

    # --- sample intersection and covariate matrix ---
    with gzip.open(cfg["phenotype"], "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    phenotype_samples = header[4:]
    cov = pd.read_csv(cfg["covariates"], sep="\t", index_col=0)
    cov.index.name = "covariate"
    psam_samples = pd.read_csv(cfg["psam"], sep="\t")["#IID"].astype(str).tolist()
    psam_set = set(psam_samples)
    cov_sample_set = set(cov.columns.astype(str))
    common_samples = [sample for sample in phenotype_samples if sample in cov_sample_set and sample in psam_set]
    if not common_samples:
        raise ValueError(f"{stage}: no common samples across phenotype/covariate/psam")

    cov_t = cov.loc[:, common_samples].T.apply(pd.to_numeric, errors="coerce")
    if cov_t.isna().any().any():
        raise ValueError(f"{stage}: covariate matrix has {int(cov_t.isna().sum().sum())} missing/non-numeric values")
    design_cov = np.column_stack([np.ones(cov_t.shape[0]), cov_t.to_numpy(dtype=float)])

    annotation_rows = []
    map_rows = []
    feature_index = 0
    incomplete_features = 0

    with residual_path.open("w", newline="") as residual_handle:
        residual_writer = csv.DictWriter(
            residual_handle,
            delimiter="\t",
            fieldnames=["model_id", *common_samples],
            lineterminator="\n",
        )
        residual_writer.writeheader()

        for chunk in pd.read_csv(cfg["phenotype"], sep="\t", compression="gzip", chunksize=CHUNK_SIZE):
            # --- residualize the chunk: complete rows in one lstsq, partial rows one by one ---
            values = chunk[common_samples].apply(pd.to_numeric, errors="coerce")
            sample_ids = list(values.columns)
            y_matrix = values.to_numpy(dtype=float)
            complete_rows = ~np.isnan(y_matrix).any(axis=1)
            residualized_chunk = [None] * y_matrix.shape[0]

            if complete_rows.any():
                y_complete = y_matrix[complete_rows, :].T
                beta, *_ = np.linalg.lstsq(design_cov, y_complete, rcond=None)
                resid = (y_complete - design_cov @ beta).T
                for output_index, row_index in enumerate(np.where(complete_rows)[0]):
                    residualized_chunk[row_index] = (sample_ids, resid[output_index, :], len(sample_ids))

            for row_index in np.where(~complete_rows)[0]:
                row = values.iloc[row_index]
                complete = row.notna() & cov_t.notna().all(axis=1)
                complete_samples = row.index[complete].tolist()
                if int(complete.sum()) < (cov_t.shape[1] + 2):
                    residualized_chunk[row_index] = (complete_samples, np.array([], dtype=float), int(complete.sum()))
                    continue
                y = row.loc[complete].to_numpy(dtype=float)
                x = cov_t.loc[complete].to_numpy(dtype=float)
                design = np.column_stack([np.ones(x.shape[0]), x])
                beta, *_ = np.linalg.lstsq(design, y, rcond=None)
                residualized_chunk[row_index] = (complete_samples, y - design @ beta, int(complete.sum()))

            for local_index, (_, row) in enumerate(chunk.iterrows()):
                complete_samples, resid, n_complete = residualized_chunk[local_index]
                feature_index += 1
                phenotype_id = str(row["pid"])
                transcript_id, gene_id_or_name, utr_region, strand = phenotype_id.split("|")
                parsed = {
                    "transcript_id": transcript_id,
                    "gene_id_or_name": gene_id_or_name,
                    "3utr_region": utr_region,
                    "strand": strand,
                }
                digest = hashlib.sha1(phenotype_id.encode("utf-8")).hexdigest()[:10]
                model_id = f"{stage}_F{feature_index:06d}_{digest}"
                transcript_start = int(utr_region.split(":", 1)[1].split("-", 1)[0])
                transcript_end = int(utr_region.rsplit("-", 1)[1])
                group_start = int(row["start"])
                group_end = int(row["end"])

                if n_complete < (cov_t.shape[1] + 2):
                    incomplete_features += 1
                    continue

                annotation_rows.append(
                    {
                        "stage": stage,
                        "model_id": model_id,
                        "feature_index": feature_index,
                        "original_phenotype_id": phenotype_id,
                        **parsed,
                        "chr": f"chr{row['#Chr']}",
                        "group_left": group_start,
                        "group_right": group_end,
                        "training_left_1mb": max(group_start - 1_000_000, 1),
                        "training_right_1mb": group_end + 1_000_000,
                        "transcript_left": transcript_start,
                        "transcript_right": transcript_end,
                        "group_interval_differs_from_transcript_interval": (
                            group_start != transcript_start or group_end != transcript_end
                        ),
                        "n_complete_samples_for_residualization": n_complete,
                    }
                )
                map_rows.append(
                    {
                        "stage": stage,
                        "model_id": model_id,
                        "original_phenotype_id": phenotype_id,
                        "transcript_id": parsed["transcript_id"],
                        "gene_id_or_name": parsed["gene_id_or_name"],
                        "3utr_region": parsed["3utr_region"],
                        "strand": parsed["strand"],
                    }
                )
                residual_row = {"model_id": model_id}
                residual_row.update(dict(zip(complete_samples, resid)))
                residual_writer.writerow(residual_row)

    annotation = pd.DataFrame(annotation_rows, columns=ANNOTATION_COLUMNS)
    if not annotation["chr"].astype(str).str.fullmatch(r"chr([1-9]|1[0-9]|2[0-2])").all():
        raise ValueError(f"{stage}: non-autosomal chr values encountered")
    pd.DataFrame(map_rows).to_csv(map_path, sep="\t", index=False, lineterminator="\n")
    annotation.to_csv(annotation_path, sep="\t", index=False, lineterminator="\n")

    audit_rows.append(
        {
            "stage": stage,
            "phenotype_samples": len(phenotype_samples),
            "covariate_samples": len(cov.columns),
            "psam_samples": len(psam_samples),
            "common_samples": len(common_samples),
            "covariates": len(cov.index),
            "features_seen": feature_index,
            "features_written": len(annotation),
            "features_skipped_low_complete_n": incomplete_features,
            "features_using_group_interval_not_transcript_interval": int(
                annotation["group_interval_differs_from_transcript_interval"].sum()
            ),
            "annotation_path": str(annotation_path.relative_to(APA)),
            "residualized_phenotype_path": str(residual_path.relative_to(APA)),
            "model_id_map_path": str(map_path.relative_to(APA)),
        }
    )

# --- audit table: replace rows for the stages just run ---
audit_path = OUT_DIR / "input_audit.tsv"
audit_df = pd.DataFrame(audit_rows)
if audit_path.exists():
    previous = pd.read_table(audit_path)
    previous = previous.loc[~previous["stage"].isin(audit_df["stage"])]
    audit_df = pd.concat([previous, audit_df], ignore_index=True)
audit_df = audit_df.sort_values("stage")
audit_df.to_csv(audit_path, sep="\t", index=False, lineterminator="\n")
print(audit_df.to_string(index=False))
