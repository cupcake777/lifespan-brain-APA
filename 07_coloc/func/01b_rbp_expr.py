#!/usr/bin/env python3
# PFC cohort RBP expression reference for Level 2 RBP support (study's own 912-sample cohort).
# Input: ${APA_HPC_HOME}/project/raw_data/all.912.tpm (cluster). Output: data/07_coloc/pfc_rbp_expression.tsv.gz
from __future__ import annotations

import gzip
import math
import os
from pathlib import Path

import yaml

# --- paths and settings ---
ROOT = Path(os.environ["APA_ROOT"])
RAW = Path(f"{os.environ['APA_HPC_HOME']}/project/raw_data")
APA_RAW = ROOT / "data" / "04_qtl"
TPM_PATH = RAW / "all.912.tpm"
META_912 = RAW / "all.912.meta"
META_CHB = APA_RAW / "CHB.all.meta"
SYMBOL_MAP = ROOT / "ref" / "gencode_v40_gene_symbol.tsv"
OUTPUT = ROOT / "data" / "07_coloc" / "pfc_rbp_expression.tsv.gz"
CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"
DETECT_THRESHOLDS = (0.1, 1.0)


# --- helpers ---
def mean_tpm(entry: dict[str, list[float]], label: str) -> str:
    count, total, _, _ = entry[label]
    return "" if count == 0 else f"{total / count:.4f}"


def pct(entry: dict[str, list[float]], label: str, which: int) -> str:
    count = entry[label][0]
    detected = entry[label][2 + which]
    return "" if count == 0 else f"{detected / count:.4f}"


# --- sample stage labels ---
# CHB.all.meta is canonical; samples missing there fall back to all.912.meta age
stage_map: dict[str, str] = {}
if META_CHB.exists():
    with META_CHB.open(encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        if "sample" not in index:
            raise ValueError(f"{META_CHB} lacks a 'sample' column")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            label = fields[index["LifeStage"]] if "LifeStage" in index else ""
            if label in ("Prenatal", "Postnatal"):
                stage_map[fields[index["sample"]]] = label

age: dict[str, float] = {}
if META_912.exists():
    with META_912.open(encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        if "sample" in index and "age" in index:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                try:
                    age[fields[index["sample"]]] = float(fields[index["age"]])
                except (ValueError, IndexError):
                    pass
for sample, years in age.items():
    if sample not in stage_map:
        stage_map[sample] = "Prenatal" if years < 0 else "Postnatal"

# --- gene symbols and CISBP RBP list ---
symbols: dict[str, str] = {}
with SYMBOL_MAP.open(encoding="utf-8", errors="ignore") as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 2 and fields[0]:
            symbols[fields[0]] = fields[1]

config = yaml.safe_load(CONFIG_PATH.read_text(encoding="utf-8"))
meme_path = Path(os.path.expandvars(config["reference_resources"]["cisbp_rna_motifs"]["path"]))
rbp_symbols: set[str] = set()
with meme_path.open(encoding="utf-8", errors="ignore") as handle:
    for line in handle:
        if line.startswith("MOTIF "):
            parts = line.split()
            if len(parts) < 3:
                continue
            name = parts[2]
            if name.startswith("("):
                name = name[1:name.index(")")] if ")" in name else name
            rbp_symbols.add(name)

# --- TPM accumulation ---
# column indices per stage in the TPM matrix (fields[0] is gene_id,
# sample columns start at fields[1], matching header[1:])
with TPM_PATH.open(encoding="utf-8", errors="ignore") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    samples = header[1:]
    stage_indices: dict[str, list[int]] = {"Prenatal": [], "Postnatal": []}
    seen: set[str] = set()
    for column, sample in enumerate(samples):
        if sample in seen:
            continue  # duplicated sample column: count once
        seen.add(sample)
        label = stage_map.get(sample)
        if label in stage_indices:
            stage_indices[label].append(column + 1)

    # gene_id -> {stage: [n, sum, d01, d1]}
    acc: dict[str, dict[str, list[float]]] = {}
    duplicates = 0
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < len(header) or not fields[0]:
            continue  # ragged/empty row: skip rather than crash the whole run
        gene_id = fields[0]
        if gene_id in acc:
            duplicates += 1
            continue
        entry: dict[str, list[float]] = {}
        for label, indices in stage_indices.items():
            total = 0.0
            count = 0
            detected = [0, 0]
            for column in indices:
                raw = fields[column]
                try:
                    value = float(raw)
                except ValueError:
                    continue
                if math.isnan(value):
                    continue
                total += value
                count += 1
                if value > DETECT_THRESHOLDS[0]:
                    detected[0] += 1
                if value > DETECT_THRESHOLDS[1]:
                    detected[1] += 1
            entry[label] = [float(count), total, float(detected[0]), float(detected[1])]
        acc[gene_id] = entry

n_prenatal = len(stage_indices["Prenatal"])
n_postnatal = len(stage_indices["Postnatal"])
print(f"TPM samples: {len(samples)} total, {n_prenatal} Prenatal, {n_postnatal} Postnatal")
print(f"TPM genes: {len(acc)} (duplicate gene_id rows skipped: {duplicates})")

# --- output ---
with gzip.open(OUTPUT, "wt", encoding="utf-8") as out:
    out.write(
        "gene_id\tgene_symbol\tn_prenatal\tprenatal_mean_tpm\tprenatal_pct_detected_0.1\t"
        "prenatal_pct_detected_1\tn_postnatal\tpostnatal_mean_tpm\tpostnatal_pct_detected_0.1\t"
        "postnatal_pct_detected_1\toverall_mean_tpm\n"
    )
    rows = []
    for gene_id, entry in acc.items():
        symbol = symbols.get(gene_id, "")
        total = sum(entry[label][1] for label in ("Prenatal", "Postnatal"))
        count = sum(entry[label][0] for label in ("Prenatal", "Postnatal"))
        overall = "" if count == 0 else f"{total / count:.4f}"
        rows.append(
            "\t".join(
                [
                    gene_id,
                    symbol,
                    f"{int(entry['Prenatal'][0])}",
                    mean_tpm(entry, "Prenatal"),
                    pct(entry, "Prenatal", 0),
                    pct(entry, "Prenatal", 1),
                    f"{int(entry['Postnatal'][0])}",
                    mean_tpm(entry, "Postnatal"),
                    pct(entry, "Postnatal", 0),
                    pct(entry, "Postnatal", 1),
                    overall,
                ]
            )
        )
    out.write("\n".join(rows) + "\n")

print(f"Written {OUTPUT} ({len(acc)} genes)")

# --- QC: CISBP RBP genes that Level 2 evidence will rely on ---
lookup = {}
for gene_id, entry in acc.items():
    symbol = symbols.get(gene_id)
    if symbol and symbol in rbp_symbols:
        lookup.setdefault(symbol, entry)
print(f"CISBP RBP gene symbols checked: {len(rbp_symbols)}; found in cohort TPM: {len(lookup)}")
missing = sorted(rbp_symbols - set(lookup))
print(f"RBP symbols NOT detected in cohort TPM (unsupported): {missing}")
for symbol in sorted(lookup):
    entry = lookup[symbol]
    print(
        f"  {symbol:<12} prenatal_mean={mean_tpm(entry, 'Prenatal'):>10} "
        f"postnatal_mean={mean_tpm(entry, 'Postnatal'):>10}"
    )
