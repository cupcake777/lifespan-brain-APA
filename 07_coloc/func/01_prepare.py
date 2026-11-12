#!/usr/bin/env python3
# Record input, resource, and software provenance.
# Assign LD signals, build the significant-QTL universe, and prepare cell references.

from __future__ import annotations

import hashlib
import importlib.metadata
import os
import platform
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pysam
import yaml
from scipy import sparse
from scipy.optimize import nnls


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"
RESULTS = ROOT / "data" / "07_coloc"
LD_DIR = RESULTS / "ld"
REFERENCE_DIR = RESULTS / "reference"
RESULTS.mkdir(parents=True, exist_ok=True)

ANALYSIS_VERSION = "framework_2026-07-25_v5"
CHECKSUM_MAX_BYTES = 1024**3
CHUNK_ROWS = 5000
PLINK2 = f"{os.environ['APA_HPC_HOME']}/.local/share/mamba/bin/plink2"
RESOURCE_OUTPUT = RESULTS / "provenance.tsv"
SOFTWARE_OUTPUT = RESULTS / "software.tsv"
OUTPUT = RESULTS / "signals.tsv.gz"
SUMMARY = RESULTS / "signals_summary.tsv"
OUT_UNIVERSE = RESULTS / "01_significant_qtl_universe.tsv.gz"
OUT_BLOCKED = RESULTS / "05_excluded_or_blocked.tsv.gz"

PRENATAL_LINEAGE_MAP = {
    "ExNeu": "Excitatory neuron",
    "IN": "Inhibitory neuron",
    "AST": "Astrocyte",
    "OPC": "OPC",
    "MG": "Microglia/Immune",
    "VASC": "Vascular/Mural",
    "GLIALPROG": "Glial progenitor",
}


# --- helpers ---

def _expand_env(obj):
    # expand ${VAR} in config values
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


def command_version(command: list[str]) -> str:
    try:
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        return f"not_found_or_failed:{type(exc).__name__}"
    output = (completed.stdout or completed.stderr).strip().splitlines()
    return output[0] if output else "unknown"


def decode(values: np.ndarray) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def categorical_values(group: h5py.Group) -> tuple[list[str], np.ndarray]:
    return decode(group["categories"][:]), np.asarray(group["codes"][:], dtype=np.int32)


def h5ad_index(handle: h5py.File, group_name: str) -> list[str]:
    group = handle[group_name]
    index_key = group.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    return decode(group[index_key][:])


def make_groups(
    donor_codes: np.ndarray,
    donor_names: list[str],
    class_codes: np.ndarray,
    classes: list[str],
    keep: np.ndarray,
    source: str,
) -> tuple[pd.DataFrame, np.ndarray]:
    pair = donor_codes.astype(np.int64) * len(classes) + class_codes.astype(np.int64)
    kept_pairs, counts = np.unique(pair[keep], return_counts=True)
    pair_to_group = {int(value): index for index, value in enumerate(kept_pairs)}
    row_to_group = np.full(donor_codes.shape[0], -1, dtype=np.int32)
    row_to_group[keep] = np.fromiter(
        (pair_to_group[int(value)] for value in pair[keep]), dtype=np.int32
    )
    records = []
    for value, n_cells in zip(kept_pairs, counts):
        donor_code = int(value // len(classes))
        class_code = int(value % len(classes))
        records.append(
            {
                "source": source,
                "donor": donor_names[donor_code],
                "cell_class": classes[class_code],
                "n_cells": int(n_cells),
            }
        )
    return pd.DataFrame(records), row_to_group


def aggregate_csr(
    h5ad_path: Path,
    x_key: str,
    row_to_group: np.ndarray,
    gene_indices: np.ndarray,
    n_groups: int,
) -> np.ndarray:
    sums = np.zeros((n_groups, len(gene_indices)), dtype=np.float64)
    with h5py.File(h5ad_path, "r") as handle:
        x = handle[x_key]
        data = x["data"]
        indices = x["indices"]
        indptr = x["indptr"]
        n_rows, n_genes = map(int, x.attrs["shape"])
        for row_start in range(0, n_rows, CHUNK_ROWS):
            row_end = min(n_rows, row_start + CHUNK_ROWS)
            local_groups = row_to_group[row_start:row_end]
            if np.all(local_groups < 0):
                continue
            ptr = np.asarray(indptr[row_start : row_end + 1], dtype=np.int64)
            data_start, data_end = int(ptr[0]), int(ptr[-1])
            chunk = sparse.csr_matrix(
                (
                    np.asarray(data[data_start:data_end], dtype=np.float64),
                    np.asarray(indices[data_start:data_end], dtype=np.int32),
                    ptr - data_start,
                ),
                shape=(row_end - row_start, n_genes),
            )[:, gene_indices]
            retained = local_groups >= 0
            retained_groups = local_groups[retained]
            unique_groups, compact = np.unique(retained_groups, return_inverse=True)
            design = sparse.csr_matrix(
                (np.ones(len(compact)), (compact, np.arange(len(compact)))),
                shape=(len(unique_groups), len(compact)),
            )
            sums[unique_groups] += (design @ chunk[retained]).toarray()
    return sums


def counts_to_cpm(counts: np.ndarray) -> np.ndarray:
    return counts / np.maximum(counts.sum(axis=1, keepdims=True), 1.0) * 1_000_000.0


def select_markers(
    cpm: np.ndarray,
    group_meta: pd.DataFrame,
    genes: list[str],
    classes: list[str],
    params: dict,
) -> pd.DataFrame:
    class_means = np.vstack(
        [cpm[group_meta["cell_class"].eq(cell_class)].mean(axis=0) for cell_class in classes]
    )
    marker_rows = []
    for class_index, cell_class in enumerate(classes):
        class_mask = group_meta["cell_class"].eq(cell_class).to_numpy()
        detection = (cpm[class_mask] > 0).mean(axis=0)
        other_max = np.delete(class_means, class_index, axis=0).max(axis=0)
        log2_fc = np.log2((class_means[class_index] + 1.0) / (other_max + 1.0))
        eligible = (
            (class_means[class_index] >= params["marker_min_cpm"])
            & (detection >= params["marker_min_detection"])
            & (log2_fc >= params["marker_min_log2_fc"])
        )
        candidates = np.flatnonzero(eligible)
        candidates = candidates[np.argsort(log2_fc[candidates])[::-1]][
            : params["marker_genes_per_class"]
        ]
        for rank, gene_index in enumerate(candidates, start=1):
            marker_rows.append(
                {
                    "cell_class": cell_class,
                    "rank": rank,
                    "gene_id": genes[gene_index],
                    "mean_cpm": class_means[class_index, gene_index],
                    "max_other_cpm": other_max[gene_index],
                    "log2_fc": log2_fc[gene_index],
                    "donor_detection": detection[gene_index],
                }
            )
    return pd.DataFrame(marker_rows)


def build_signature(
    cpm: np.ndarray,
    group_meta: pd.DataFrame,
    genes: list[str],
    classes: list[str],
    markers: pd.DataFrame,
) -> pd.DataFrame:
    marker_genes = markers.sort_values(["cell_class", "rank"])["gene_id"].drop_duplicates()
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    marker_indices = np.array([gene_lookup[gene] for gene in marker_genes], dtype=int)
    signature = {"NAME": marker_genes.tolist()}
    for cell_class in classes:
        mask = group_meta["cell_class"].eq(cell_class).to_numpy()
        signature[cell_class] = cpm[mask][:, marker_indices].mean(axis=0)
    return pd.DataFrame(signature)


def lodo_benchmark(
    cpm: np.ndarray,
    group_meta: pd.DataFrame,
    genes: list[str],
    classes: list[str],
    params: dict,
    stage: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(params["benchmark_seed"] + (0 if stage == "prenatal" else 1))
    donors = sorted(group_meta["donor"].unique())
    estimates = []
    metrics = []
    for donor in donors:
        test_mask = group_meta["donor"].eq(donor).to_numpy()
        train_mask = ~test_mask
        test_meta = group_meta.loc[test_mask].reset_index(drop=True)
        if test_meta["cell_class"].nunique() < 5:
            continue
        train_markers = select_markers(
            cpm[train_mask], group_meta.loc[train_mask].reset_index(drop=True), genes, classes, params
        )
        marker_genes = train_markers["gene_id"].drop_duplicates().tolist()
        gene_lookup = {gene: index for index, gene in enumerate(genes)}
        marker_indices = np.array([gene_lookup[gene] for gene in marker_genes], dtype=int)
        design = np.vstack(
            [
                cpm[train_mask & group_meta["cell_class"].eq(cell_class).to_numpy()][
                    :, marker_indices
                ].mean(axis=0)
                for cell_class in classes
            ]
        ).T
        test_profiles = {
            row.cell_class: cpm[index, marker_indices]
            for index, row in zip(np.flatnonzero(test_mask), test_meta.itertuples())
        }
        available = [cell_class for cell_class in classes if cell_class in test_profiles]
        for replicate in range(1, params["synthetic_mixtures_per_donor"] + 1):
            weights = rng.dirichlet(np.ones(len(available)))
            truth = dict.fromkeys(classes, 0.0)
            truth.update(dict(zip(available, weights)))
            mixture = sum(weight * test_profiles[cell_class] for weight, cell_class in zip(weights, available))
            coefficients, _ = nnls(design, mixture)
            coefficients = coefficients / coefficients.sum()
            mix_id = f"{stage}_{donor}_mix{replicate}"
            for class_index, cell_class in enumerate(classes):
                estimates.append(
                    {
                        "stage": stage,
                        "donor": donor,
                        "mixture": mix_id,
                        "cell_class": cell_class,
                        "truth": truth[cell_class],
                        "estimate": coefficients[class_index],
                    }
                )
            truth_vector = np.array([truth[cell_class] for cell_class in classes])
            metrics.append(
                {
                    "stage": stage,
                    "donor": donor,
                    "mixture": mix_id,
                    "n_available_classes": len(available),
                    "n_markers": len(marker_indices),
                    "mae": np.mean(np.abs(coefficients - truth_vector)),
                    "rmse": np.sqrt(np.mean((coefficients - truth_vector) ** 2)),
                    "pearson": np.corrcoef(coefficients, truth_vector)[0, 1],
                }
            )
    return pd.DataFrame(estimates), pd.DataFrame(metrics)


def write_cibersortx_mixture(
    bulk: pd.DataFrame,
    samples: list[str],
    marker_genes: list[str],
    output: Path,
) -> None:
    mixture = bulk[["gene_id", *samples]].copy()
    mixture["gene_id"] = mixture["gene_id"].str.replace(r"\.\d+$", "", regex=True)
    mixture = mixture.groupby("gene_id", sort=False, as_index=False).mean(numeric_only=True)
    mixture = mixture[mixture["gene_id"].isin(marker_genes)]
    mixture = mixture.rename(columns={"gene_id": "NAME"})
    mixture.to_csv(output, sep="\t", index=False)


# --- load config ---

with CONFIG_PATH.open(encoding="utf-8") as handle:
    config = _expand_env(yaml.safe_load(handle))
inputs = config["inputs"]


# --- resource and input manifest ---

rows = []
for category, mapping in [("reference_resource", config["reference_resources"]), ("input", inputs)]:
    for name, specification in mapping.items():
        if not isinstance(specification, dict):
            continue
        path_fields = (["path"] if "path" in specification else []) + [
            field
            for field in ("h5ad", "metadata", "covariates", "liftover_chain")
            if str(specification.get(field, "")).strip()
        ]
        if not path_fields:
            continue
        for path_field in path_fields:
            path_text = str(specification[path_field]).strip()
            row_name = name if path_field == "path" else f"{name}.{path_field}"
            if not path_text or path_text.startswith("python_module:"):
                rows.append({
                    "category": category,
                    "name": row_name,
                    "configured_path": path_text,
                    "resolved_path": path_text,
                    "release": specification.get("release", "not_recorded"),
                    "object_type": "logical_or_absent",
                    "size_bytes": pd.NA,
                    "mtime": pd.NA,
                    "sha256": pd.NA,
                    "sha256_status": "not_applicable",
                })
                continue
            path = Path(path_text)
            resolved = path.resolve()
            if not resolved.exists():
                rows.append({
                    "category": category,
                    "name": row_name,
                    "configured_path": path_text,
                    "resolved_path": str(resolved),
                    "release": specification.get("release", "not_recorded"),
                    "object_type": "missing",
                    "size_bytes": pd.NA,
                    "mtime": pd.NA,
                    "sha256": pd.NA,
                    "sha256_status": "missing",
                })
                continue
            stat = resolved.stat()
            if resolved.is_file() and stat.st_size <= CHECKSUM_MAX_BYTES:
                digest = hashlib.sha256()
                with resolved.open("rb") as handle:
                    for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                        digest.update(block)
                checksum = digest.hexdigest()
                checksum_status = "computed"
            elif resolved.is_file():
                checksum = pd.NA
                checksum_status = "skipped_above_1GiB_recorded_size_mtime"
            else:
                checksum = pd.NA
                checksum_status = "directory_recorded_path_mtime"
            rows.append({
                "category": category,
                "name": row_name,
                "configured_path": path_text,
                "resolved_path": str(resolved),
                "release": specification.get("release", "not_recorded"),
                "object_type": "file" if resolved.is_file() else "directory",
                "size_bytes": stat.st_size if resolved.is_file() else pd.NA,
                "mtime": pd.Timestamp(stat.st_mtime, unit="s", tz="Asia/Shanghai"),
                "sha256": checksum,
                "sha256_status": checksum_status,
            })

resource_table = pd.DataFrame(rows)
resource_table["analysis_version"] = ANALYSIS_VERSION
resource_table["config_path"] = str(CONFIG_PATH)
resource_table.to_csv(RESOURCE_OUTPUT, sep="\t", index=False)


# --- software versions ---

software_rows = [
    {"software": "Python", "version": platform.python_version(), "source": platform.python_implementation()},
    {"software": "PLINK2", "version": command_version(["plink2", "--version"]), "source": "executable"},
    {"software": "RNAfold", "version": command_version(["RNAfold", "--version"]), "source": "executable"},
    {"software": "R", "version": command_version(["Rscript", "-e", "cat(R.version.string)"]), "source": "Rscript"},
]
for package in ["numpy", "pandas", "pysam", "PyYAML", "scipy", "statsmodels"]:
    software_rows.append({
        "software": package,
        "version": importlib.metadata.version(package),
        "source": "python_package",
    })
for package in ["data.table", "ggplot2", "clusterProfiler", "org.Hs.eg.db", "coloc"]:
    version = command_version(["Rscript", "-e", f"cat(as.character(packageVersion('{package}')))"])
    software_rows.append({"software": package, "version": version, "source": "R_package"})

software_table = pd.DataFrame(software_rows)
software_table["analysis_version"] = ANALYSIS_VERSION
software_table.to_csv(SOFTWARE_OUTPUT, sep="\t", index=False)


# --- LD signal assignment ---

# --- loading inputs ---

frames = []
for stage in ["prenatal", "postnatal"]:
    association_path = Path(inputs[f"{stage}_significant_apaqtl"]["path"])
    fine_path = Path(inputs[f"{stage}_finemapped_apaqtl"]["path"])
    association = pd.read_csv(
        association_path,
        sep="\t",
        usecols=["phenotype_id", "variant_id", "pos", "pval_nominal"],
    )
    fine = pd.read_parquet(fine_path, columns=["phenotype_id", "variant_id", "pip", "cs_id"])
    fine = fine.drop_duplicates(["phenotype_id", "variant_id"])
    association = association.merge(fine, on=["phenotype_id", "variant_id"], how="left")
    association.insert(0, "stage", stage)
    frames.append(association)
associations = pd.concat(frames, ignore_index=True)


# --- stage-matched dosage export and assignment ---

signal_config = config["parameters"]["signal_assignment"]
pfile_prefix = config["reference_resources"]["ld_reference"]["pfile_prefix"]
window_bp = int(signal_config["ld_window_kb"]) * 1000
r2_threshold = float(signal_config["ld_r2"])

assigned_frames = []
summary_rows = []
for stage in ["prenatal", "postnatal"]:
    stage_associations = associations.loc[associations["stage"].eq(stage)].copy()
    variants = sorted(stage_associations["variant_id"].astype(str).unique())
    sample_path = Path(inputs[f"{stage}_qtl_samples"]["path"])

    # export dosages for this stage's variants and samples
    stage_dir = LD_DIR / stage
    stage_dir.mkdir(parents=True, exist_ok=True)
    extract_path = stage_dir / "variants.txt"
    keep_path = stage_dir / "samples.keep"
    prefix = stage_dir / "dosage"
    extract_path.write_text("\n".join(variants) + "\n", encoding="utf-8")
    sample_ids = [line.strip() for line in sample_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    keep_path.write_text("#IID\n" + "\n".join(sample_ids) + "\n", encoding="utf-8")
    subprocess.run(
        [
            PLINK2,
            "--pfile", pfile_prefix,
            "--extract", str(extract_path),
            "--keep", str(keep_path),
            "--export", "Av",
            "--out", str(prefix),
        ],
        check=True,
    )
    dosage_table = pd.read_csv(prefix.with_suffix(".traw"), sep="\t")
    sample_columns = dosage_table.columns[6:]
    dosage = {
        str(row.SNP): row.loc[sample_columns].to_numpy(dtype=np.float32)
        for _, row in dosage_table.iterrows()
    }
    sample_n = len(sample_columns)
    missing = sorted(set(variants) - set(dosage))
    if missing:
        raise ValueError(f"{stage}: {len(missing)} significant variants are absent from the LD reference")

    # SuSiE credible-set leads first, then LD-clumped leads for the rest
    rows = []
    for phenotype_id, group in stage_associations.groupby("phenotype_id", sort=False):
        group = group.sort_values(["pval_nominal", "variant_id"], kind="stable").copy()
        leads: list[tuple[str, int, str]] = []

        credible = group.loc[group["cs_id"].notna() & group["cs_id"].astype(str).str.strip().ne("")].copy()
        for cs_id, members in credible.groupby("cs_id", sort=False):
            member = members.sort_values(["pip", "pval_nominal"], ascending=[False, True], kind="stable").iloc[0]
            signal_id = f"{member.stage}::{phenotype_id}::susie_cs::{cs_id}"
            lead_variant = str(member.variant_id)
            leads.append((lead_variant, int(member.pos), signal_id))
            for item in members.itertuples(index=False):
                rows.append({
                    "stage": item.stage,
                    "phenotype_id": phenotype_id,
                    "variant_id": item.variant_id,
                    "signal_id": signal_id,
                    "signal_assignment_method": "susie_credible_set",
                    "signal_lead_variant": lead_variant,
                    "ld_r2_to_lead": np.nan,
                })

        unassigned = group.loc[~group.index.isin(credible.index)]
        for item in unassigned.itertuples(index=False):
            variant = str(item.variant_id)
            position = int(item.pos)
            candidates: list[tuple[float, str, str]] = []
            for lead_variant, lead_position, signal_id in leads:
                if abs(position - lead_position) > window_bp:
                    continue
                r2 = np.nan
                left, right = dosage[variant], dosage[lead_variant]
                valid = np.isfinite(left) & np.isfinite(right)
                if valid.sum() >= 3:
                    x = left[valid] - left[valid].mean()
                    y = right[valid] - right[valid].mean()
                    denominator = np.sqrt(np.dot(x, x) * np.dot(y, y))
                    if denominator != 0:
                        r2 = float((np.dot(x, y) / denominator) ** 2)
                if np.isfinite(r2) and r2 >= r2_threshold:
                    candidates.append((r2, lead_variant, signal_id))
            if candidates:
                r2, lead_variant, signal_id = max(candidates, key=lambda value: value[0])
                method = "ld_to_susie_signal" if "::susie_cs::" in signal_id else "ld_clump"
            else:
                lead_variant = variant
                signal_id = f"{item.stage}::{phenotype_id}::ldclump::{variant}"
                r2 = 1.0
                method = "ld_clump_lead"
                leads.append((lead_variant, position, signal_id))
            rows.append({
                "stage": item.stage,
                "phenotype_id": phenotype_id,
                "variant_id": variant,
                "signal_id": signal_id,
                "signal_assignment_method": method,
                "signal_lead_variant": lead_variant,
                "ld_r2_to_lead": r2,
            })
    assigned = pd.DataFrame(rows)

    assigned_frames.append(assigned)
    summary_rows.append({
        "stage": stage,
        "association_n": len(stage_associations),
        "variant_n": len(variants),
        "sample_n": sample_n,
        "signal_n": assigned["signal_id"].nunique(),
        "susie_signal_n": assigned.loc[assigned["signal_assignment_method"].eq("susie_credible_set"), "signal_id"].nunique(),
        "ld_clump_signal_n": assigned.loc[assigned["signal_assignment_method"].eq("ld_clump_lead"), "signal_id"].nunique(),
    })


# --- write signal assignments ---

assignments = pd.concat(assigned_frames, ignore_index=True)
if assignments.duplicated(["stage", "phenotype_id", "variant_id"]).any():
    raise ValueError("Signal assignments are not unique by stage, phenotype, and variant")

assignments["ld_reference"] = pfile_prefix
assignments["ld_measure"] = signal_config["ld_measure"]
assignments["ld_window_kb"] = int(signal_config["ld_window_kb"])
assignments["ld_r2_threshold"] = r2_threshold
assignments["analysis_version"] = ANALYSIS_VERSION
assignments.to_csv(OUTPUT, sep="\t", index=False, compression="gzip")
pd.DataFrame(summary_rows).to_csv(SUMMARY, sep="\t", index=False)


# --- candidate universe ---

# --- loading significant association universe ---

stage_inputs = {
    "prenatal": Path(inputs["prenatal_significant_apaqtl"]["path"]),
    "postnatal": Path(inputs["postnatal_significant_apaqtl"]["path"]),
}
frames = []
for stage, path in stage_inputs.items():
    frame = pd.read_csv(path, sep="\t")
    frame.insert(0, "source_row_id", [f"{stage}:{i + 1}" for i in range(len(frame))])
    frame.insert(1, "stage", stage)
    frame["source_file"] = str(path)
    frames.append(frame)
universe = pd.concat(frames, ignore_index=True, sort=False)
universe["candidate_id"] = universe["stage"] + "::" + universe["phenotype_id"].astype(str) + "::" + universe["variant_id"].astype(str)

# phenotype_id and variant_id parsing
parts = universe["phenotype_id"].astype(str).str.extract(
    r"^(?P<transcript_id>[^|]+)\|(?P<gene_name>[^|]+)\|"
    r"(?P<phenotype_chr>[^:|]+):(?P<phenotype_start>\d+)-(?P<phenotype_end>\d+)\|"
    r"(?P<phenotype_strand>[+-])$"
)
universe = pd.concat([universe.reset_index(drop=True), parts], axis=1)
universe["phenotype_parse_ok"] = parts.notna().all(axis=1)
universe["phenotype_start"] = pd.to_numeric(universe["phenotype_start"], errors="coerce").astype("Int64")
universe["phenotype_end"] = pd.to_numeric(universe["phenotype_end"], errors="coerce").astype("Int64")

parts = universe["variant_id"].astype(str).str.extract(
    r"^(?P<variant_chr>[^:]+):(?P<variant_pos>\d+):(?P<variant_ref>[^:]+):(?P<variant_alt>[^:]+)$"
)
universe = pd.concat([universe.reset_index(drop=True), parts], axis=1)
universe["variant_parse_ok"] = parts.notna().all(axis=1)
universe["variant_pos"] = pd.to_numeric(universe["variant_pos"], errors="coerce").astype("Int64")
universe["variant_start0"] = universe["variant_pos"] - 1
universe["variant_end"] = universe["variant_start0"] + universe["variant_ref"].str.len().astype("Int64")
universe["variant_fields_match_input"] = (
    universe["variant_parse_ok"]
    & universe["variant_chr"].eq(universe["chr"].astype(str))
    & universe["variant_pos"].eq(pd.to_numeric(universe["pos"], errors="coerce").astype("Int64"))
    & universe["variant_ref"].str.upper().eq(universe["ref"].astype(str).str.upper())
    & universe["variant_alt"].str.upper().eq(universe["alt"].astype(str).str.upper())
)


# --- reference allele and direction ---

fasta_path = Path(config["reference_resources"]["genome_fasta"]["path"])
allele_variants = universe.loc[universe["variant_parse_ok"], ["variant_id", "variant_chr", "variant_pos", "variant_ref"]].drop_duplicates()
fasta = pysam.FastaFile(str(fasta_path))
references = set(fasta.references)
allele_rows = []
for row in allele_variants.itertuples(index=False):
    chrom = str(row.variant_chr)
    if chrom not in references and chrom.removeprefix("chr") in references:
        chrom = chrom.removeprefix("chr")
    if chrom not in references and f"chr{chrom}" in references:
        chrom = f"chr{chrom}"
    observed = None
    status = "chromosome_not_in_reference"
    if chrom in references:
        start0 = int(row.variant_pos) - 1
        observed = fasta.fetch(chrom, start0, start0 + len(str(row.variant_ref))).upper()
        status = "ok" if observed == str(row.variant_ref).upper() else "reference_allele_mismatch"
    allele_rows.append({"variant_id": row.variant_id, "reference_genome_allele": observed, "reference_allele_status": status})
fasta.close()
reference = pd.DataFrame(allele_rows)
universe = universe.merge(reference, on="variant_id", how="left")
universe["allele_harmonized"] = universe["variant_fields_match_input"] & universe["reference_allele_status"].eq("ok")
universe["beta_alt_PDUI"] = pd.to_numeric(universe["slope"], errors="coerce")
universe["direction_blocked"] = ~universe["allele_harmonized"] | universe["beta_alt_PDUI"].isna() | universe["beta_alt_PDUI"].eq(0)
universe["apa_direction"] = "direction_blocked"
universe.loc[~universe["direction_blocked"] & universe["beta_alt_PDUI"].gt(0), "apa_direction"] = "lengthening"
universe.loc[~universe["direction_blocked"] & universe["beta_alt_PDUI"].lt(0), "apa_direction"] = "shortening"


# --- fine-mapping evidence ---

fine_frames = []
for stage, key in [("prenatal", "prenatal_finemapped_apaqtl"), ("postnatal", "postnatal_finemapped_apaqtl")]:
    fine = pd.read_parquet(inputs[key]["path"])
    fine.insert(0, "stage", stage)
    fine_frames.append(fine)
fine = pd.concat(fine_frames, ignore_index=True).drop_duplicates(["stage", "phenotype_id", "variant_id"])
universe = universe.merge(
    fine[["stage", "phenotype_id", "variant_id", "pip", "cs_id"]],
    on=["stage", "phenotype_id", "variant_id"],
    how="left",
)
thresholds = config["candidate_sets"]["evidence_thresholds"]
universe["credible_set_member"] = universe["cs_id"].notna() & universe["cs_id"].astype(str).str.strip().ne("")
numeric_pip = pd.to_numeric(universe["pip"], errors="coerce")
finemap_tier = pd.Series("not_finemapped", index=universe["pip"].index, dtype="object")
finemap_tier.loc[numeric_pip.gt(thresholds["pip_exploratory"]["value"])] = "exploratory"
finemap_tier.loc[numeric_pip.ge(thresholds["pip_high_confidence"]["value"])] = "high_confidence"
universe["finemap_tier"] = finemap_tier


# --- interaction and disease evidence layers ---

interaction = pd.read_csv(inputs["interaction_apaqtl"]["path"], sep="\t")
interaction_cols = ["phenotype_id", "variant_id", "b_g", "b_gi", "pval_gi", "pval_adj_bh"]
interaction = interaction[interaction_cols].drop_duplicates(["phenotype_id", "variant_id"])
universe = universe.merge(interaction, on=["phenotype_id", "variant_id"], how="left")
universe["interaction_significant"] = pd.to_numeric(universe["pval_adj_bh"], errors="coerce").le(0.05).fillna(False)

interaction_groups = pd.read_csv(inputs["interaction_group_results"]["path"], sep="\t")
interaction_group_columns = [
    "group_id", "gene_id", "lead_phenotype_id", "lead_variant_id", "previous_115", "new_83",
]
interaction_groups = interaction_groups[interaction_group_columns].copy()
if interaction_groups["group_id"].duplicated().any():
    raise ValueError("Current HC3-core interaction table is not unique by group_id")
if interaction_groups.duplicated(["lead_phenotype_id", "lead_variant_id"]).any():
    raise ValueError("Current HC3-core interaction table is not unique by representative pair")
interaction_groups["interaction_representative_phenotype_id"] = interaction_groups["lead_phenotype_id"]
interaction_groups["interaction_representative_variant_id"] = interaction_groups["lead_variant_id"]
interaction_groups = interaction_groups.rename(columns={
    "group_id": "interaction_group_id",
    "gene_id": "interaction_gene_id",
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "previous_115": "interaction_previous_115",
    "new_83": "interaction_new_83",
})
universe = universe.merge(
    interaction_groups,
    on=["phenotype_id", "variant_id"],
    how="left",
    validate="many_to_one",
)

gxstage = pd.read_csv(inputs["gxstage_classes"]["path"], sep="\t")
gxstage = gxstage[[
    "lead_phenotype_id", "lead_variant_id", "final_pattern", "classification_status",
    "pattern_tau", "pattern_delta", "pattern_method",
]].rename(columns={
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "final_pattern": "gxstage_class",
    "classification_status": "gxstage_classification_status",
})
if gxstage.duplicated(["phenotype_id", "variant_id"]).any():
    raise ValueError("Pinned GxStage contract is not unique by representative pair")
valid_gxstage_classes = {"Opposing", "Amplified", "Attenuated", "Uncertain"}
if not gxstage["gxstage_classification_status"].eq("formal").all():
    raise ValueError("Pinned GxStage contract contains non-formal classifications")
if not gxstage["gxstage_class"].isin(valid_gxstage_classes).all():
    raise ValueError("Pinned GxStage contract contains unsupported posterior classes")
gxstage_group_check = gxstage.merge(
    interaction_groups[["phenotype_id", "variant_id", "interaction_group_id"]],
    on=["phenotype_id", "variant_id"],
    how="left",
    validate="one_to_one",
)
if gxstage_group_check["interaction_group_id"].isna().any():
    raise ValueError("Pinned GxStage contract contains representative pairs absent from the current HC3 core")

posterior_source = pd.read_csv(inputs["interaction_posterior_results"]["path"], sep="\t")
posterior_source = posterior_source[[
    "group_id", "lead_phenotype_id", "lead_variant_id", "posterior_pattern",
    "hc3_pattern_tau", "hc3_pattern_delta",
]].rename(columns={
    "group_id": "posterior_interaction_group_id",
    "lead_phenotype_id": "phenotype_id",
    "lead_variant_id": "variant_id",
    "posterior_pattern": "source_posterior_pattern",
    "hc3_pattern_tau": "source_pattern_tau",
    "hc3_pattern_delta": "source_pattern_delta",
})
if posterior_source["posterior_interaction_group_id"].duplicated().any():
    raise ValueError("Interaction posterior source is not unique by group_id")
if posterior_source.duplicated(["phenotype_id", "variant_id"]).any():
    raise ValueError("Interaction posterior source is not unique by representative pair")
gxstage_source_check = gxstage_group_check.merge(
    posterior_source,
    on=["phenotype_id", "variant_id"],
    how="left",
    validate="one_to_one",
    indicator=True,
)
if not gxstage_source_check["_merge"].eq("both").all():
    raise ValueError("Pinned GxStage contract contains representative pairs absent from the posterior source")
if not gxstage_source_check["interaction_group_id"].astype(str).eq(
    gxstage_source_check["posterior_interaction_group_id"].astype(str)
).all():
    raise ValueError("Pinned GxStage contract maps to inconsistent interaction group IDs")
if not gxstage_source_check["gxstage_class"].astype(str).eq(
    gxstage_source_check["source_posterior_pattern"].astype(str)
).all():
    raise ValueError("Pinned GxStage classes differ from the 03_interaction posterior source")
for contract_column, source_column in [
    ("pattern_tau", "source_pattern_tau"),
    ("pattern_delta", "source_pattern_delta"),
]:
    contract_values = pd.to_numeric(gxstage_source_check[contract_column], errors="coerce")
    source_values = pd.to_numeric(gxstage_source_check[source_column], errors="coerce")
    if contract_values.isna().any() or source_values.isna().any() or not np.isclose(
        contract_values.to_numpy(), source_values.to_numpy(), rtol=0, atol=1e-12
    ).all():
        raise ValueError(
            f"Pinned GxStage {contract_column} differs from the 03_interaction posterior source"
        )
universe = universe.merge(gxstage, on=["phenotype_id", "variant_id"], how="left", validate="many_to_one")
current_group_without_pattern = universe["interaction_group_id"].notna() & universe["gxstage_class"].isna()
universe.loc[current_group_without_pattern, "gxstage_classification_status"] = "posterior_class_not_available"
universe["gxstage_class_source"] = np.select(
    [universe["gxstage_class"].notna(), current_group_without_pattern],
    [
        "03_interaction_previous115_posterior_pattern_representative_pair",
        "03_interaction_current_hc3_core_without_formal_posterior_pattern",
    ],
    default="not_interaction_group_representative_pair",
)

coloc = pd.read_csv(inputs["gwas_coloc_results"]["path"])
coloc_summary = (
    coloc.rename(columns={"QTL_ID": "stage", "Phenotype": "phenotype_id", "PP4": "gwas_coloc_max_pp4"})
    .groupby(["stage", "phenotype_id"], as_index=False)
    .agg(gwas_coloc_max_pp4=("gwas_coloc_max_pp4", "max"), gwas_coloc_row_count=("GWAS_ID", "size"))
)
universe = universe.merge(coloc_summary, on=["stage", "phenotype_id"], how="left")
universe["has_gwas_coloc_result"] = universe["gwas_coloc_row_count"].fillna(0).gt(0)


# --- independent signal assignment ---

assignment_path = Path(inputs["independent_signal_assignments"]["path"])
signal_assignments = pd.read_csv(assignment_path, sep="\t")
assignment_columns = [
    "stage", "phenotype_id", "variant_id", "signal_id", "signal_assignment_method",
    "signal_lead_variant", "ld_r2_to_lead", "ld_reference", "ld_measure",
    "ld_window_kb", "ld_r2_threshold",
]
signal_assignments = signal_assignments[assignment_columns]
if signal_assignments.duplicated(["stage", "phenotype_id", "variant_id"]).any():
    raise ValueError("Independent-signal assignments contain duplicate association keys")
universe = universe.merge(signal_assignments, on=["stage", "phenotype_id", "variant_id"], how="left", validate="one_to_one")
if universe["signal_id"].isna().any():
    raise ValueError(f"{universe['signal_id'].isna().sum()} significant associations lack signal assignment")
universe["independent_signal_assigned"] = True
universe["signal_representative"] = universe["variant_id"].eq(universe["signal_lead_variant"])
if universe.groupby("signal_id")["signal_representative"].sum().ne(1).any():
    raise ValueError("Every independent signal must have exactly one representative variant")

universe["analysis_version"] = ANALYSIS_VERSION
universe["config_path"] = str(CONFIG_PATH)
universe["reference_build"] = "GRCh38"
universe["genome_fasta_release"] = config["reference_resources"]["genome_fasta"]["release"]
universe["qtl_significant_definition"] = "permutation_fdr_phenotype_plus_nominal_pair_threshold"
universe["qtl_universe"] = "all_significant_pairs"


# --- write candidate universe ---

universe.sort_values(["stage", "phenotype_id", "variant_pos"], kind="stable").to_csv(
    OUT_UNIVERSE, sep="\t", index=False, compression="gzip"
)

blocked_checks = [
    (~universe["phenotype_parse_ok"], "phenotype_parse", "invalid_phenotype_id", "blocked"),
    (~universe["variant_parse_ok"], "variant_parse", "invalid_variant_id", "blocked"),
    (~universe["variant_fields_match_input"], "variant_fields", "variant_id_disagrees_with_input_columns", "blocked"),
    (~universe["reference_allele_status"].eq("ok"), "reference_allele", "reference_allele_not_verified", "blocked_for_allele_analysis"),
    (universe["direction_blocked"], "direction", "alt_effect_direction_not_interpretable", "kept_with_blocked_direction"),
    (~universe["independent_signal_assigned"], "signal_assignment", "missing_required_signal_assignment", "blocked"),
]
keep = ["candidate_id", "source_row_id", "stage", "phenotype_id", "variant_id"]
blocked_parts = []
for mask, step, reason, severity in blocked_checks:
    part = universe.loc[mask, keep].copy()
    part["step"] = step
    part["reason"] = reason
    part["severity"] = severity
    blocked_parts.append(part)
blocked = pd.concat(blocked_parts, ignore_index=True)
blocked["analysis_version"] = ANALYSIS_VERSION
blocked["config_path"] = str(CONFIG_PATH)
blocked.to_csv(OUT_BLOCKED, sep="\t", index=False, compression="gzip")


# --- cell-reference preparation ---

REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

# --- loading bulk gene namespace ---

resources = config["reference_resources"]
params = config["parameters"]["cell_deconvolution"]

prenatal_samples = Path(inputs["prenatal_qtl_samples"]["path"]).read_text().splitlines()
postnatal_samples = Path(inputs["postnatal_qtl_samples"]["path"]).read_text().splitlines()
all_samples = prenatal_samples + postnatal_samples
bulk_path = Path(inputs["bulk_rna_tpm"]["path"])
bulk = pd.read_csv(bulk_path, sep="\t", usecols=["gene_id", *all_samples])
bulk_genes = set(bulk["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True))


# --- prenatal donor-level reference ---

vel_h5ad = Path(resources["velmeshev_scrna"]["h5ad"])
vel_meta_path = Path(resources["velmeshev_scrna"]["metadata"])
vel_meta = pd.read_csv(
    vel_meta_path,
    sep="\t",
    usecols=["Cell_ID", "Age", "Individual", "Lineage"],
    dtype=str,
)
with h5py.File(vel_h5ad, "r") as handle:
    vel_obs = h5ad_index(handle, "obs")
    vel_genes = decode(handle["raw"]["var"]["features"][:])
if vel_obs != vel_meta["Cell_ID"].tolist():
    raise SystemExit("Velmeshev metadata row order does not match h5ad obs order")

prenatal_classes = params["prenatal_classes"]
gestational = vel_meta["Age"].str.contains(r"(?i)(?:^GW|^ga|\bGW\b)", regex=True, na=False)
mapped_class = vel_meta["Lineage"].map(PRENATAL_LINEAGE_MAP)
vel_donor_codes, vel_donors = pd.factorize(vel_meta["Individual"], sort=True)
vel_class_codes = mapped_class.map({value: index for index, value in enumerate(prenatal_classes)})
vel_keep = gestational.to_numpy() & vel_class_codes.notna().to_numpy()
prenatal_meta, prenatal_row_to_group = make_groups(
    vel_donor_codes,
    vel_donors.astype(str).tolist(),
    vel_class_codes.fillna(-1).to_numpy(dtype=np.int32),
    prenatal_classes,
    vel_keep,
    "Velmeshev_prenatal",
)
prenatal_shared = [gene for gene in vel_genes if gene in bulk_genes]
vel_gene_lookup = {gene: index for index, gene in enumerate(vel_genes)}
prenatal_gene_indices = np.array([vel_gene_lookup[gene] for gene in prenatal_shared])
prenatal_counts = aggregate_csr(
    vel_h5ad,
    "raw/X",
    prenatal_row_to_group,
    prenatal_gene_indices,
    len(prenatal_meta),
)
prenatal_cpm = counts_to_cpm(prenatal_counts)
prenatal_markers = select_markers(
    prenatal_cpm, prenatal_meta, prenatal_shared, prenatal_classes, params
)
prenatal_markers.insert(0, "stage", "prenatal")
prenatal_signature = build_signature(
    prenatal_cpm, prenatal_meta, prenatal_shared, prenatal_classes, prenatal_markers
)


# --- postnatal donor-level reference ---

full_h5ad = Path(resources["fullard_aging_scrna"]["h5ad"])
postnatal_classes = params["postnatal_classes"]
with h5py.File(full_h5ad, "r") as handle:
    full_genes = h5ad_index(handle, "var")
    donor_names, donor_codes = categorical_values(handle["obs"]["donor_id"])
    disease_names, disease_codes = categorical_values(handle["obs"]["disease"])
    stage_names, stage_codes = categorical_values(handle["obs"]["development_stage"])
    class_names, class_codes = categorical_values(handle["obs"]["class"])
    subclass_names, subclass_codes = categorical_values(handle["obs"]["subclass"])

normal_code = disease_names.index(params["postnatal_disease_filter"])
unknown_stage_code = stage_names.index("unknown")
full_class_lookup = {
    "EN": "Excitatory neuron",
    "IN": "Inhibitory neuron",
    "Astro": "Astrocyte",
    "Oligo": "Oligodendrocyte",
    "OPC": "OPC",
    "Immune": "Microglia/Immune",
    "Endo": "Vascular/Mural",
    "Mural": "Vascular/Mural",
}
target_lookup = {value: index for index, value in enumerate(postnatal_classes)}
full_target_codes = np.full(len(class_names), -1, dtype=np.int32)
for source_class, cell_class in full_class_lookup.items():
    full_target_codes[class_names.index(source_class)] = target_lookup[cell_class]
full_cell_codes = full_target_codes[class_codes]
if "PVM" in subclass_names:
    full_cell_codes[subclass_codes == subclass_names.index("PVM")] = target_lookup["Vascular/Mural"]
full_keep = (
    (disease_codes == normal_code)
    & (stage_codes != unknown_stage_code)
    & (full_cell_codes >= 0)
)
postnatal_meta, postnatal_row_to_group = make_groups(
    donor_codes,
    donor_names,
    full_cell_codes,
    postnatal_classes,
    full_keep,
    "Fullard_normal_postnatal",
)
postnatal_shared = [gene for gene in full_genes if gene in bulk_genes]
full_gene_lookup = {gene: index for index, gene in enumerate(full_genes)}
postnatal_gene_indices = np.array([full_gene_lookup[gene] for gene in postnatal_shared])
postnatal_counts = aggregate_csr(
    full_h5ad,
    "X",
    postnatal_row_to_group,
    postnatal_gene_indices,
    len(postnatal_meta),
)
postnatal_cpm = counts_to_cpm(postnatal_counts)
postnatal_markers = select_markers(
    postnatal_cpm, postnatal_meta, postnatal_shared, postnatal_classes, params
)
postnatal_markers.insert(0, "stage", "postnatal")
postnatal_signature = build_signature(
    postnatal_cpm, postnatal_meta, postnatal_shared, postnatal_classes, postnatal_markers
)


# --- leave-one-donor-out signature recovery ---

prenatal_estimates, prenatal_benchmark = lodo_benchmark(
    prenatal_cpm, prenatal_meta, prenatal_shared, prenatal_classes, params, "prenatal"
)
postnatal_estimates, postnatal_benchmark = lodo_benchmark(
    postnatal_cpm, postnatal_meta, postnatal_shared, postnatal_classes, params, "postnatal"
)


# --- CIBERSORTx inputs and reference outputs ---

prenatal_signature.to_csv(
    REFERENCE_DIR / "cibersortx_signature_prenatal.tsv", sep="\t", index=False
)
postnatal_signature.to_csv(
    REFERENCE_DIR / "cibersortx_signature_postnatal.tsv", sep="\t", index=False
)
write_cibersortx_mixture(
    bulk,
    prenatal_samples,
    prenatal_signature["NAME"].tolist(),
    REFERENCE_DIR / "cibersortx_mixture_prenatal.tsv",
)
write_cibersortx_mixture(
    bulk,
    postnatal_samples,
    postnatal_signature["NAME"].tolist(),
    REFERENCE_DIR / "cibersortx_mixture_postnatal.tsv",
)

marker_lookup_prenatal = {gene: index for index, gene in enumerate(prenatal_shared)}
prenatal_marker_indices = [marker_lookup_prenatal[gene] for gene in prenatal_signature["NAME"]]
prenatal_pseudobulk = pd.concat(
    [prenatal_meta, pd.DataFrame(prenatal_counts[:, prenatal_marker_indices], columns=prenatal_signature["NAME"])],
    axis=1,
)
marker_lookup_postnatal = {gene: index for index, gene in enumerate(postnatal_shared)}
postnatal_marker_indices = [marker_lookup_postnatal[gene] for gene in postnatal_signature["NAME"]]
postnatal_pseudobulk = pd.concat(
    [postnatal_meta, pd.DataFrame(postnatal_counts[:, postnatal_marker_indices], columns=postnatal_signature["NAME"])],
    axis=1,
)
prenatal_pseudobulk.to_parquet(REFERENCE_DIR / "pseudobulk_markers_prenatal.parquet", index=False)
postnatal_pseudobulk.to_parquet(REFERENCE_DIR / "pseudobulk_markers_postnatal.parquet", index=False)
reference_tables = {
    "cell_reference_markers.tsv": pd.concat([prenatal_markers, postnatal_markers]),
    "cell_reference_lodo_estimates.tsv": pd.concat([prenatal_estimates, postnatal_estimates]),
    "cell_reference_lodo_benchmark.tsv": pd.concat([prenatal_benchmark, postnatal_benchmark]),
}
for name, table in reference_tables.items():
    table["analysis_version"] = ANALYSIS_VERSION
    table["config_path"] = str(CONFIG_PATH)
    table.to_csv(RESULTS / name, sep="\t", index=False)

reference_qc = pd.DataFrame(
    [
        {
            "stage": "prenatal",
            "source": "Velmeshev gestational-age nuclei",
            "source_h5ad": str(vel_h5ad),
            "source_h5ad_bytes": vel_h5ad.stat().st_size,
            "source_metadata": str(vel_meta_path),
            "n_donors": prenatal_meta["donor"].nunique(),
            "n_profiles": len(prenatal_meta),
            "n_cells": prenatal_meta["n_cells"].sum(),
            "n_classes": len(prenatal_classes),
            "n_shared_genes": len(prenatal_shared),
            "n_markers": prenatal_signature.shape[0],
        },
        {
            "stage": "postnatal",
            "source": "Fullard disease-normal known-age DLPFC nuclei",
            "source_h5ad": str(full_h5ad),
            "source_h5ad_bytes": full_h5ad.stat().st_size,
            "source_metadata": pd.NA,
            "n_donors": postnatal_meta["donor"].nunique(),
            "n_profiles": len(postnatal_meta),
            "n_cells": postnatal_meta["n_cells"].sum(),
            "n_classes": len(postnatal_classes),
            "n_shared_genes": len(postnatal_shared),
            "n_markers": postnatal_signature.shape[0],
        },
    ]
)
reference_qc["bulk_tpm_path"] = str(bulk_path)
reference_qc["analysis_version"] = ANALYSIS_VERSION
reference_qc["config_path"] = str(CONFIG_PATH)
reference_qc.to_csv(RESULTS / "cell_reference_qc.tsv", sep="\t", index=False)
