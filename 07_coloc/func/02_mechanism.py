#!/usr/bin/env python3
# Classify the SNP-to-APA mechanism of each independent QTL signal.
# Direct PAS core changes, PAS-local RBP motifs, RNA structure, APARENT, phastCons and POSTAR3 overlap.

from __future__ import annotations

from bisect import bisect_right
from collections import Counter
import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import pyBigWig
import pysam
import yaml


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data" / "07_coloc"
CONFIG_PATH = Path(__file__).resolve().parent / "config.yml"
ASSOCIATIONS = DATA / "01_significant_qtl_universe.tsv.gz"
OUTPUT = DATA / "02_signal_mechanisms.tsv.gz"
PAS_DETAIL_OUTPUT = DATA / "pas_detail.tsv.gz"
RBP_DETAIL_OUTPUT = DATA / "rbp_detail.tsv.gz"
RBPDB_HUMAN_PROTEINS = ROOT / "ref" / "RBPDB_v1.3.1_proteins_human_2012-11-21.csv"
ANALYSIS_VERSION = "framework_2026-09-12_v2"
RBP_DELTA_THRESHOLD = 0.20
POSTAR_RBP_LIMIT = 20
POSTAR_CHUNK_SIZE = 250_000
APARENT_MODEL_NAME = "aparent_large_lessdropout_all_libs_no_sampleweights.h5"
APARENT_WINDOW_HALF = 102
PIPELINE_RUN_ID = os.environ.get("APA_FUNC_RUN_ID", "manual")
GXSTAGE_FIELDS = [
    "interaction_group_id",
    "interaction_gene_id",
    "interaction_representative_phenotype_id",
    "interaction_representative_variant_id",
    "interaction_previous_115",
    "interaction_new_83",
    "gxstage_class",
    "gxstage_classification_status",
    "gxstage_class_source",
    "pattern_tau",
    "pattern_delta",
    "pattern_method",
]
GXSTAGE_LINK_FIELDS = [
    "interaction_pair_signal_assignment_method",
    "interaction_pair_signal_lead_variant",
    "interaction_pair_ld_r2_to_signal_lead",
    "interaction_pair_cs_id",
    "gxstage_mechanism_variant_link",
    "gxstage_mechanism_variant_link_r2",
]


# --- helpers ---

# expand ${VAR} in config values (paths set via config/paths.sh)
def _expand_env(obj):
    if isinstance(obj, dict):
        return {k: _expand_env(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_env(v) for v in obj]
    return os.path.expandvars(obj) if isinstance(obj, str) else obj


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def transcript_allele(allele: str, strand: str) -> str:
    allele = str(allele).upper()
    return allele if strand == "+" else revcomp(allele)


def transcript_relative_position(variant_pos: int, pas_pos: int, strand: str) -> int:
    return variant_pos - pas_pos if strand == "+" else pas_pos - variant_pos


def fetch_pas_window(fasta: pysam.FastaFile, chrom: str, pas_pos: int, strand: str) -> str:
    reference_chrom = chrom
    if reference_chrom not in fasta.references and reference_chrom.removeprefix("chr") in fasta.references:
        reference_chrom = reference_chrom.removeprefix("chr")
    if strand == "+":
        sequence = fasta.fetch(reference_chrom, pas_pos - 101, pas_pos + 50).upper()
    else:
        sequence = revcomp(fasta.fetch(reference_chrom, pas_pos - 51, pas_pos + 100).upper())
    return sequence


def motif_hits(sequence: str, motifs: set[str]) -> Counter:
    hits = Counter()
    for motif in motifs:
        for index in range(len(sequence) - len(motif) + 1):
            if sequence[index:index + len(motif)] == motif:
                hits[motif] += 1
    return hits


# CISBP-RNA MEME file -> list of PWMs with max/min log-odds scores
def parse_cisbp_meme(path: Path) -> list[dict]:
    motifs = []
    motif_id = None
    rbp_name = None
    width = None
    rows = []

    def add_motif() -> None:
        if motif_id is None or rbp_name is None or width is None or len(rows) != width:
            return
        max_score = sum(math.log2(max(max(row), 1e-6) / 0.25) for row in rows)
        min_score = sum(math.log2(max(min(row), 1e-6) / 0.25) for row in rows)
        motifs.append({
            "motif_id": motif_id,
            "rbp_name": rbp_name,
            "width": width,
            "pwm": rows.copy(),
            "max_score": max_score,
            "min_score": min_score,
        })

    with path.open(encoding="utf-8", errors="ignore") as handle:
        in_matrix = False
        for raw_line in handle:
            line = raw_line.strip()
            if line.startswith("MOTIF "):
                add_motif()
                parts = line.split()
                motif_id = parts[1]
                rbp_name = parts[2]
                width = None
                rows = []
                in_matrix = False
            elif line.startswith("letter-probability matrix"):
                match = re.search(r"w=\s*(\d+)", line)
                width = int(match.group(1)) if match else None
                rows = []
                in_matrix = width is not None
            elif in_matrix and len(line.split()) == 4:
                rows.append([float(value) for value in line.split()])
                if width is not None and len(rows) == width:
                    add_motif()
                    motif_id = None
                    rbp_name = None
                    width = None
                    rows = []
                    in_matrix = False
    add_motif()
    return motifs


def normalized_pwm_score(sequence: str, motif: dict) -> float:
    if len(sequence) != motif["width"]:
        return 0.0
    index = {"A": 0, "C": 1, "G": 2, "T": 3}
    score = 0.0
    for base, probabilities in zip(sequence, motif["pwm"], strict=True):
        if base not in index:
            return 0.0
        score += math.log2(max(probabilities[index[base]], 1e-6) / 0.25)
    denominator = motif["max_score"] - motif["min_score"]
    return (score - motif["min_score"]) / denominator if denominator > 0 else 0.0


def oriented_variant_sequences(row: pd.Series, fasta: pysam.FastaFile, flank: int) -> tuple[str, str, int, str]:
    chrom = str(row["variant_chr"])
    if chrom not in fasta.references and chrom.removeprefix("chr") in fasta.references:
        chrom = chrom.removeprefix("chr")
    start1 = max(1, int(row["variant_pos"]) - flank)
    end1 = int(row["variant_pos"]) + flank
    reference = fasta.fetch(chrom, start1 - 1, end1).upper()
    if row["phenotype_strand"] == "-":
        reference = revcomp(reference)
        offset = end1 - int(row["variant_pos"])
    else:
        offset = int(row["variant_pos"]) - start1
    ref = transcript_allele(row["variant_ref"], row["phenotype_strand"])
    alt = transcript_allele(row["variant_alt"], row["phenotype_strand"])
    if reference[offset:offset + len(ref)] != ref:
        return reference, reference, offset, "reference_conflict"
    alternate = reference[:offset] + alt + reference[offset + len(ref):]
    return reference, alternate, offset, "ok"


# CISBP rbp_name -> human gene symbol; ortholog motifs such as '(HNRNPR)_(Gallus_gallus)_(RBD_0.97)' give the name in parentheses
def rbp_gene_symbol(rbp_name: object) -> str:
    name = str(rbp_name)
    if name.startswith("(") and ")" in name:
        return name[1:name.index(")")]
    return name


def expression_info(symbol: str, expression: dict | None) -> dict | None:
    if expression is None:
        return None
    return expression.get(symbol) or expression.get(symbol.upper())


# Level 2 support: RBP gene expressed in the cohort at the signal stage (mean TPM >= min_mean_tpm).
# No expression table -> legacy behaviour (sequence-only evidence).
def rbp_expression_supported(
    symbol: str, expression: dict | None, min_mean_tpm: float, stage: str
) -> bool:
    if expression is None:
        return True
    info = expression_info(symbol, expression)
    if info is None:
        return False
    mean_tpm = pd.to_numeric(info.get(stage.capitalize()), errors="coerce")
    return bool(pd.notna(mean_tpm) and mean_tpm >= min_mean_tpm)


def stage_of(row: pd.Series) -> str:
    stage = row.get("stage")
    if isinstance(stage, str) and stage:
        return stage
    candidate_id = str(row.get("candidate_id", ""))
    if candidate_id.startswith("prenatal::"):
        return "prenatal"
    return "postnatal"


# pfc_rbp_expression table -> expression lookup keyed by current symbols and RBPDB aliases
def load_rbp_expression(path: Path | None) -> dict[str, dict[str, float]] | None:
    if path is None or not Path(path).exists():
        return None
    table = pd.read_csv(path, sep="\t")
    lookup: dict[str, dict[str, float]] = {}
    by_gene_base: dict[str, dict[str, float]] = {}

    def add_key(key: object, info: dict[str, float]) -> None:
        if pd.isna(key):
            return
        text = str(key).strip()
        if not text or text == "\\N":
            return
        lookup.setdefault(text, info)
        lookup.setdefault(text.upper(), info)

    for record in table.to_dict("records"):
        symbol = record["gene_symbol"]
        if pd.isna(symbol) or not str(symbol).strip():
            continue
        gene_id = str(record["gene_id"])
        info = {
            "expression_gene_id": gene_id,
            "expression_gene_symbol": str(symbol),
            "Prenatal": pd.to_numeric(record["prenatal_mean_tpm"], errors="coerce"),
            "Postnatal": pd.to_numeric(record["postnatal_mean_tpm"], errors="coerce"),
        }
        add_key(symbol, info)
        by_gene_base[gene_id.split(".", 1)[0]] = info

    if RBPDB_HUMAN_PROTEINS.exists():
        with RBPDB_HUMAN_PROTEINS.open(newline="", encoding="utf-8", errors="ignore") as handle:
            for fields in csv.reader(handle):
                if len(fields) < 10:
                    continue
                info = by_gene_base.get(fields[1].split(".", 1)[0])
                if info is None:
                    continue
                add_key(fields[4], info)
                for alias in fields[9].split(";"):
                    add_key(alias, info)
    if "TUT1" in lookup:
        add_key("STAR-PAP", lookup["TUT1"])
    return lookup


def score_local_rbp(
    row: pd.Series,
    fasta: pysam.FastaFile,
    motifs: list[dict],
    expression: dict | None = None,
    min_mean_tpm: float = 1.0,
    stage: str = "postnatal",
) -> dict:
    max_width = max(motif["width"] for motif in motifs)
    reference, alternate, offset, status = oriented_variant_sequences(row, fasta, max_width - 1)
    best_gain = (None, 0.0)
    best_loss = (None, 0.0)
    for motif in motifs:
        # Level 2: only RBP genes with own-cohort, stage-matched expression count as evidence
        if not rbp_expression_supported(rbp_gene_symbol(motif["rbp_name"]), expression, min_mean_tpm, stage):
            continue
        width = motif["width"]
        start_min = max(0, offset - width + 1)
        start_max = min(offset, len(reference) - width)
        ref_best = 0.0
        alt_best = 0.0
        for start in range(start_min, start_max + 1):
            ref_best = max(ref_best, normalized_pwm_score(reference[start:start + width], motif))
            alt_best = max(alt_best, normalized_pwm_score(alternate[start:start + width], motif))
        delta = alt_best - ref_best
        if delta > best_gain[1]:
            best_gain = (motif["rbp_name"], delta)
        if delta < best_loss[1]:
            best_loss = (motif["rbp_name"], delta)
    if best_gain[0] is None and best_loss[0] is None:
        max_delta = 0.0
    else:
        max_delta = max(abs(best_gain[1]), abs(best_loss[1]))
    return {
        "candidate_id": row["candidate_id"],
        "v2_top_rbp_gain": best_gain[0],
        "v2_top_rbp_gain_delta": best_gain[1],
        "v2_top_rbp_loss": best_loss[0],
        "v2_top_rbp_loss_delta": best_loss[1],
        "v2_rbp_variant_delta_support": max_delta >= RBP_DELTA_THRESHOLD,
        "v2_rbp_scoring_status": status,
    }


def score_local_structure(table: pd.DataFrame, fasta: pysam.FastaFile, threshold: float) -> pd.DataFrame:
    records = []
    sequences = []
    for row in table.to_dict("records"):
        reference, alternate, _, status = oriented_variant_sequences(pd.Series(row), fasta, 100)
        records.append((row["candidate_id"], status))
        sequences.extend([(f">{row['candidate_id']}::ref", reference), (f">{row['candidate_id']}::alt", alternate)])
    if not sequences:
        return pd.DataFrame(columns=["candidate_id", "local_structure_delta_mfe", "local_structure_change", "local_structure_status"])
    payload = "\n".join(item for record in sequences for item in record) + "\n"
    process = subprocess.run(["RNAfold", "--noPS"], input=payload, text=True, capture_output=True, check=True)
    energies = {}
    current = None
    for line in process.stdout.splitlines():
        if line.startswith(">"):
            current = line[1:]
        elif current and "(" in line:
            match = re.search(r"\(\s*(-?[0-9.]+)\)\s*$", line)
            if match:
                energies[current] = float(match.group(1))
    rows = []
    for candidate_id, mutation_status in records:
        ref_mfe = energies.get(f"{candidate_id}::ref")
        alt_mfe = energies.get(f"{candidate_id}::alt")
        delta = alt_mfe - ref_mfe if ref_mfe is not None and alt_mfe is not None else np.nan
        rows.append({
            "candidate_id": candidate_id,
            "local_structure_delta_mfe": delta,
            "local_structure_change": pd.notna(delta) and abs(delta) >= threshold,
            "local_structure_status": "ok" if pd.notna(delta) and mutation_status == "ok" else mutation_status,
        })
    return pd.DataFrame(rows)


def direct_pas_change(
    row: pd.Series,
    pas_pos: int,
    pas_name: str,
    fasta: pysam.FastaFile,
    motifs: set[str],
) -> dict:
    relative = transcript_relative_position(int(row["variant_pos"]), int(pas_pos), row["phenotype_strand"])
    result = {
        "pas_name": pas_name,
        "pas_pos": int(pas_pos),
        "variant_relative_to_pas_nt": relative,
        "pas_core_eligible": -50 <= relative <= -1,
        "pas_core_change": "not_eligible",
        "ref_pas_hits": "",
        "alt_pas_hits": "",
        "predicted_pdui_direction": "not_available",
    }
    if not result["pas_core_eligible"]:
        return result

    sequence = fetch_pas_window(fasta, row["phenotype_chr"], int(pas_pos), row["phenotype_strand"])
    index = relative + 100
    ref = transcript_allele(row["variant_ref"], row["phenotype_strand"])
    alt = transcript_allele(row["variant_alt"], row["phenotype_strand"])
    if sequence[index:index + len(ref)] != ref:
        result["pas_core_change"] = "reference_conflict"
        return result
    alt_sequence = sequence[:index] + alt + sequence[index + len(ref):]
    ref_hits = motif_hits(sequence[50:100], motifs)
    alt_hits = motif_hits(alt_sequence[50:100], motifs)
    result["ref_pas_hits"] = ";".join(sorted(ref_hits.elements()))
    result["alt_pas_hits"] = ";".join(sorted(alt_hits.elements()))
    if sum(alt_hits.values()) > sum(ref_hits.values()):
        change = "Gain"
    elif sum(alt_hits.values()) < sum(ref_hits.values()):
        change = "Loss"
    elif alt_hits != ref_hits:
        change = "Modulation"
    else:
        change = "No core change"
    result["pas_core_change"] = change
    if change in {"Gain", "Loss"}:
        proximal_gain = pas_name == "proximal" and change == "Gain"
        distal_loss = pas_name == "distal" and change == "Loss"
        result["predicted_pdui_direction"] = "shortening" if proximal_gain or distal_loss else "lengthening"
    return result


def collapse_pas_results(results: list[dict]) -> dict:
    direct = [item for item in results if item["pas_core_change"] in {"Gain", "Loss"}]
    directions = {item["predicted_pdui_direction"] for item in direct}
    direct_sites = {item["pas_name"] for item in direct}
    return {
        "proximal_pas_core_change": results[0]["pas_core_change"],
        "distal_pas_core_change": results[1]["pas_core_change"],
        "proximal_variant_distance_nt": results[0]["variant_relative_to_pas_nt"],
        "distal_variant_distance_nt": results[1]["variant_relative_to_pas_nt"],
        "direct_pas_evidence": bool(direct),
        "direct_pas_sites": ";".join(sorted(direct_sites)),
        "direct_pas_anchor_high_confidence": "proximal" in direct_sites,
        "direct_pas_predicted_pdui_direction": next(iter(directions)) if len(directions) == 1 else "conflicting" if directions else "not_available",
        "direct_pas_reference_conflict": any(item["pas_core_change"] == "reference_conflict" for item in results),
    }


def point_conservation(bw: pyBigWig.pyBigWig, chrom: str, pos: int) -> float | pd.NA:
    try:
        value = bw.values(chrom, pos - 1, pos)[0]
    except RuntimeError:
        return pd.NA
    return float(value) if value is not None else pd.NA


def annotate_conservation(table: pd.DataFrame, path: Path) -> pd.DataFrame:
    bw = pyBigWig.open(str(path))
    rows = []
    try:
        for row in table.itertuples(index=False):
            rows.append({
                "candidate_id": row.candidate_id,
                "variant_phastcons": point_conservation(bw, row.variant_chr, int(row.variant_pos)),
                "proximal_pas_phastcons": (
                    point_conservation(bw, row.phenotype_chr, int(row.proximal_pas_pos))
                    if pd.notna(row.proximal_pas_pos)
                    else pd.NA
                ),
                "distal_pas_phastcons": (
                    point_conservation(bw, row.phenotype_chr, int(row.distal_pas_pos))
                    if pd.notna(row.distal_pas_pos)
                    else pd.NA
                ),
            })
    finally:
        bw.close()
    return pd.DataFrame(rows)


def build_interval_index(table: pd.DataFrame, start_column: str, end_column: str) -> dict[str, tuple[list[int], list[tuple[int, int, str]]]]:
    index = {}
    valid = table.loc[table[start_column].notna() & table[end_column].notna()]
    for chrom, group in valid.groupby("phenotype_chr", sort=False):
        intervals = sorted(
            (int(row[start_column]), int(row[end_column]), row["candidate_id"])
            for row in group.to_dict("records")
        )
        index[str(chrom)] = ([item[0] for item in intervals], intervals)
    return index


def stream_postar_overlaps(table: pd.DataFrame, path: Path) -> pd.DataFrame:
    variant_index = build_interval_index(table.assign(variant_start=table["variant_pos"], variant_end=table["variant_pos"]), "variant_start", "variant_end")
    candidate_chroms = set(table["phenotype_chr"].astype(str))
    counts = {candidate_id: 0 for candidate_id in table["candidate_id"]}
    rbps = {candidate_id: set() for candidate_id in table["candidate_id"]}
    technologies = {candidate_id: set() for candidate_id in table["candidate_id"]}
    for chunk in pd.read_csv(path, usecols=[0, 4, 20, 21], chunksize=POSTAR_CHUNK_SIZE, encoding="utf-8-sig"):
        chunk = chunk.dropna(subset=["l1"])
        coordinates = chunk["l1"].astype(str).str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
        chunk = pd.concat([chunk, coordinates], axis=1).dropna(subset=["chrom", "start", "end"])
        chunk = chunk.loc[chunk["chrom"].isin(candidate_chroms)]
        for row in chunk.itertuples(index=False):
            rbp = str(row.RBP).strip() if pd.notna(row.RBP) and str(row.RBP).strip() else str(row.RBP_Query).strip()
            technology = str(row.CLIP_Technology).strip() if pd.notna(row.CLIP_Technology) else "unknown_technology"
            starts, intervals = variant_index.get(str(row.chrom), ([], []))
            hits = [
                candidate_id
                for _, interval_end, candidate_id in intervals[:bisect_right(starts, int(row.end))]
                if interval_end >= int(row.start)
            ]
            for candidate_id in hits:
                counts[candidate_id] += 1
                rbps[candidate_id].add(rbp or "unknown_rbp")
                technologies[candidate_id].add(technology or "unknown_technology")
    return pd.DataFrame({
        "candidate_id": list(counts),
        "postar_variant_overlap_count": list(counts.values()),
        "postar_variant_rbps": [";".join(sorted(rbps[candidate_id])[:POSTAR_RBP_LIMIT]) for candidate_id in counts],
        "postar_variant_technologies": [";".join(sorted(technologies[candidate_id])[:POSTAR_RBP_LIMIT]) for candidate_id in counts],
    })


def score_aparent_variants(table: pd.DataFrame, fasta: pysam.FastaFile, model_directory: Path) -> pd.DataFrame:
    sys.path.insert(0, str(model_directory))
    import keras
    from aparent.predictor import get_aparent_encoder

    model = keras.models.load_model(
        str(model_directory / "saved_models" / APARENT_MODEL_NAME),
        compile=False,
    )
    encoder = get_aparent_encoder()
    rows = []
    for row in table.itertuples(index=False):
        scores = {}
        for label, pas_pos in (("proximal", row.proximal_pas_pos), ("distal", row.distal_pas_pos)):
            delta_column = f"APARENT_delta_{label}"
            status_column = f"aparent_scoring_status_{label}"
            in_window_column = f"aparent_variant_in_window_{label}"
            if pd.isna(pas_pos):
                scores[delta_column] = np.nan
                scores[status_column] = "no_pas_pos"
                scores[in_window_column] = False
                continue
            relative = transcript_relative_position(int(row.variant_pos), int(pas_pos), row.phenotype_strand)
            if not -APARENT_WINDOW_HALF <= relative <= APARENT_WINDOW_HALF:
                scores[delta_column] = 0.0
                scores[status_column] = "variant_outside_window_delta_zero"
                scores[in_window_column] = False
                continue
            # window of APARENT_WINDOW_HALF nt on each side of the PAS, reverse-complemented on the minus strand
            reference_chrom = row.variant_chr if row.variant_chr in fasta.references else row.variant_chr.removeprefix("chr")
            reference = fasta.fetch(reference_chrom, int(pas_pos) - APARENT_WINDOW_HALF - 1, int(pas_pos) + APARENT_WINDOW_HALF).upper()
            if row.phenotype_strand != "+":
                reference = revcomp(reference)
            index = relative + APARENT_WINDOW_HALF
            ref = transcript_allele(row.variant_ref, row.phenotype_strand)
            alt = transcript_allele(row.variant_alt, row.phenotype_strand)
            if len(reference) != 205 or reference[index:index + len(ref)] != ref:
                scores[delta_column] = np.nan
                scores[status_column] = "reference_conflict"
                scores[in_window_column] = False
                continue
            alternate = reference[:index] + alt + reference[index + len(ref):]
            reference_probability = float(np.clip(model.predict(encoder([reference]), verbose=0)[0][0, 0], 1e-7, 1 - 1e-7))
            alternate_probability = float(np.clip(model.predict(encoder([alternate]), verbose=0)[0][0, 0], 1e-7, 1 - 1e-7))
            scores[delta_column] = math.log(alternate_probability / (1 - alternate_probability)) - math.log(reference_probability / (1 - reference_probability))
            scores[status_column] = "ok"
            scores[in_window_column] = True
        predicted = scores["APARENT_delta_distal"] - scores["APARENT_delta_proximal"]
        rows.append({
            "candidate_id": row.candidate_id,
            **scores,
            "APARENT_predicted_PDUI_delta": predicted,
            "aparent_scoring_status": "ok" if any(scores[f"aparent_variant_in_window_{label}"] for label in ("proximal", "distal")) else "variant_outside_both_windows",
        })
    return pd.DataFrame(rows)


def write_gzip_table(table: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(path.name + ".tmp")
    table.to_csv(temporary, sep="\t", index=False, compression="gzip")
    temporary.replace(path)


def bool_value(value: object) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


# --- load signals and PAS anchors ---

with CONFIG_PATH.open(encoding="utf-8") as handle:
    config = _expand_env(yaml.safe_load(handle))
thresholds = config["candidate_sets"]["evidence_thresholds"]
pip_exploratory = float(thresholds["pip_exploratory"]["value"])
pip_high = float(thresholds["pip_high_confidence"]["value"])

rbp_expression_cfg = config["reference_resources"].get("pfc_rbp_expression") or {}
rbp_expression_path = Path(rbp_expression_cfg["path"]) if rbp_expression_cfg.get("path") else None
rbp_expression_min_mean_tpm = float(
    config["parameters"].get("rbp_expression", {}).get("min_mean_tpm", 1.0)
)
if rbp_expression_path is not None and not rbp_expression_path.exists():
    raise FileNotFoundError(f"Configured pfc_rbp_expression reference is missing: {rbp_expression_path}")
rbp_expression = load_rbp_expression(rbp_expression_path)
if rbp_expression is None:
    print("pfc_rbp_expression reference absent: Level 2 RBP expression support disabled "
          "(sequence-only CISBP evidence)")

associations = pd.read_csv(ASSOCIATIONS, sep="\t")
interaction_signal_map = associations.loc[
    associations["interaction_group_id"].notna(),
    [
        "signal_id", *GXSTAGE_FIELDS, "signal_assignment_method", "signal_lead_variant",
        "ld_r2_to_lead", "cs_id",
    ],
].rename(columns={
    "signal_assignment_method": "interaction_pair_signal_assignment_method",
    "signal_lead_variant": "interaction_pair_signal_lead_variant",
    "ld_r2_to_lead": "interaction_pair_ld_r2_to_signal_lead",
    "cs_id": "interaction_pair_cs_id",
}).drop_duplicates()
if interaction_signal_map["signal_id"].isna().any():
    raise ValueError("Interaction representative pairs must map to a stage-specific QTL signal")
if interaction_signal_map["signal_id"].duplicated().any():
    duplicated = interaction_signal_map.loc[
        interaction_signal_map["signal_id"].duplicated(keep=False), "signal_id"
    ].nunique()
    raise ValueError(f"{duplicated} stage-specific QTL signals map to multiple interaction groups")
representatives = associations.loc[
    associations["independent_signal_assigned"].fillna(False).astype(bool)
    & associations["signal_representative"].fillna(False).astype(bool),
    ["signal_id", "candidate_id", "variant_id", "pip"],
].rename(columns={
    "candidate_id": "representative_candidate_id",
    "variant_id": "representative_variant_id",
    "pip": "representative_variant_pip",
})
signals = associations.loc[
    associations["independent_signal_assigned"].fillna(False).astype(bool)
].copy()

event_info = pd.read_csv(config["inputs"]["apa_phenotype_annotation"]["path"], sep="\t")
event_info = event_info[["phenotype_id", "Predicted_Proximal_APA", "Loci"]].drop_duplicates("phenotype_id")
event_info["proximal_pas_pos"] = pd.to_numeric(event_info["Predicted_Proximal_APA"], errors="coerce").astype("Int64")
signals = signals.merge(event_info, on="phenotype_id", how="left")
signals["distal_pas_pos"] = signals["phenotype_end"]
negative = signals["phenotype_strand"].eq("-")
signals.loc[negative, "distal_pas_pos"] = signals.loc[negative, "phenotype_start"]
signals["dur_start"] = signals[["proximal_pas_pos", "distal_pas_pos"]].min(axis=1).astype("Int64")
signals["dur_end"] = signals[["proximal_pas_pos", "distal_pas_pos"]].max(axis=1).astype("Int64")
signals["dur_length"] = (signals["distal_pas_pos"] - signals["proximal_pas_pos"]).abs().astype("Int64")
signals["pas_anchor_resolution"] = "proximal_pipeline_distal_boundary_approximation"


# --- direct PAS annotation ---

motifs = set(config["parameters"]["pas_hexamers"])
fasta = pysam.FastaFile(config["reference_resources"]["genome_fasta"]["path"])
direct_rows = []
for row in signals.to_dict("records"):
    series = pd.Series(row)
    if pd.isna(series["proximal_pas_pos"]) or pd.isna(series["distal_pas_pos"]) or not bool_value(series["allele_harmonized"]):
        direct_rows.append({
            "candidate_id": series["candidate_id"],
            "direct_pas_evidence": False,
            "direct_pas_sites": "",
            "direct_pas_anchor_high_confidence": False,
            "direct_pas_predicted_pdui_direction": "not_available",
            "direct_pas_reference_conflict": False,
            "proximal_pas_core_change": "not_scored",
            "distal_pas_core_change": "not_scored",
            "proximal_variant_distance_nt": pd.NA,
            "distal_variant_distance_nt": pd.NA,
        })
        continue
    site_results = [
        direct_pas_change(series, int(series["proximal_pas_pos"]), "proximal", fasta, motifs),
        direct_pas_change(series, int(series["distal_pas_pos"]), "distal", fasta, motifs),
    ]
    direct_rows.append({"candidate_id": series["candidate_id"], **collapse_pas_results(site_results)})
signals = signals.merge(pd.DataFrame(direct_rows), on="candidate_id", how="left")


# --- PAS-local RBP and RNA structure ---

signals["local_regulatory_eligible"] = (
    pd.to_numeric(signals["proximal_variant_distance_nt"], errors="coerce").between(-100, -1)
    | pd.to_numeric(signals["proximal_variant_distance_nt"], errors="coerce").between(10, 50)
    | pd.to_numeric(signals["distal_variant_distance_nt"], errors="coerce").between(-100, -1)
    | pd.to_numeric(signals["distal_variant_distance_nt"], errors="coerce").between(10, 50)
)
local_candidates = signals.loc[
    signals["local_regulatory_eligible"] & signals["allele_harmonized"].map(bool_value)
].copy()
cisbp_motifs = parse_cisbp_meme(Path(config["reference_resources"]["cisbp_rna_motifs"]["path"]))
rbp_rows = [
    score_local_rbp(
        pd.Series(row), fasta, cisbp_motifs,
        expression=rbp_expression,
        min_mean_tpm=rbp_expression_min_mean_tpm,
        stage=stage_of(row),
    )
    for row in local_candidates.to_dict("records")
]
if rbp_rows:
    signals = signals.merge(pd.DataFrame(rbp_rows), on="candidate_id", how="left")
structure_threshold = float(config["parameters"]["rna_structure"]["local_delta_mfe_threshold_kcal_mol"])
structure = score_local_structure(local_candidates, fasta, structure_threshold)
signals = signals.merge(structure, on="candidate_id", how="left")


# --- APARENT, conservation and POSTAR3 ---

aparent = score_aparent_variants(
    signals,
    fasta,
    Path(config["reference_resources"]["aparent_model"]["path"]),
)
signals = signals.merge(aparent, on="candidate_id", how="left")
conservation = annotate_conservation(
    signals,
    Path(config["reference_resources"]["phastcons_bigwig"]["path"]),
)
signals = signals.merge(conservation, on="candidate_id", how="left")
postar = stream_postar_overlaps(
    local_candidates,
    Path(config["reference_resources"]["postar3_clip_peaks"]["path"]) / "ALL_RBP_binding_sites_merged.csv",
)
signals = signals.merge(postar, on="candidate_id", how="left")


# --- local regulatory evidence ---

for column in ["v2_top_rbp_gain", "v2_top_rbp_loss", "v2_rbp_variant_delta_support", "local_structure_delta_mfe"]:
    if column not in signals.columns:
        signals[column] = pd.NA
signals["top_rbp_gain"] = signals["v2_top_rbp_gain"]
signals["top_rbp_loss"] = signals["v2_top_rbp_loss"]
signals["postar_variant_overlap_count"] = pd.to_numeric(
    signals["postar_variant_overlap_count"], errors="coerce"
).fillna(0).astype(int)
signals["local_rbp_motif_change"] = signals["local_regulatory_eligible"] & (
    signals["v2_rbp_variant_delta_support"].map(bool_value)
)
signals["RNA_structure_support"] = signals["local_structure_change"].map(bool_value)
signals["local_annotation_status"] = np.where(
    signals["local_regulatory_eligible"],
    "computed_cisbp_rnafold_postar_aparent_phastcons",
    "not_eligible_for_pas_local_model",
)
signals["local_clip_support"] = (
    signals["local_regulatory_eligible"]
    & pd.to_numeric(signals["postar_variant_overlap_count"], errors="coerce").fillna(0).gt(0)
)
signals["local_structure_change"] = signals["local_regulatory_eligible"] & signals["local_structure_change"].map(bool_value)
signals["local_rbp_or_structure_evidence"] = signals["local_rbp_motif_change"] | signals["local_structure_change"]
direction_ready = (
    signals["apa_direction"].isin(["lengthening", "shortening"])
    & pd.to_numeric(signals["beta_alt_PDUI"], errors="coerce").notna()
    & pd.to_numeric(signals["APARENT_predicted_PDUI_delta"], errors="coerce").notna()
)
signals["aparent_beta_direction_consistent"] = pd.NA
signals.loc[direction_ready, "aparent_beta_direction_consistent"] = (
    np.sign(pd.to_numeric(signals.loc[direction_ready, "beta_alt_PDUI"], errors="coerce"))
    == np.sign(pd.to_numeric(signals.loc[direction_ready, "APARENT_predicted_PDUI_delta"], errors="coerce"))
)

# evidence label: first matching rule wins
evidence_labels = []
for _, row in signals.iterrows():
    pip = pd.to_numeric(row.get("pip"), errors="coerce")
    direct = bool_value(row.get("direct_pas_evidence"))
    local = bool_value(row.get("local_rbp_or_structure_evidence"))
    local_supported = local and (
        bool_value(row.get("local_clip_support")) or bool_value(row.get("local_structure_change"))
    )
    predicted = direct or local or pd.notna(row.get("APARENT_predicted_PDUI_delta"))
    if direct and pd.notna(pip) and pip >= pip_high:
        label = "direct_high_confidence"
    elif direct:
        label = "direct_exploratory"
    elif local_supported and pd.notna(pip) and pip > pip_exploratory:
        label = "local_regulatory_supported"
    elif predicted:
        label = "predicted_only"
    else:
        label = "association_mechanism_unresolved"
    evidence_labels.append(label)
signals["snp_to_apa_evidence_label"] = evidence_labels
signals["mechanism_direction_concordant"] = (
    signals["direct_pas_predicted_pdui_direction"].eq(signals["apa_direction"])
).where(signals["direct_pas_predicted_pdui_direction"].isin(["lengthening", "shortening"]), pd.NA)
signals["mechanism_star_eligible"] = (
    pd.to_numeric(signals["pip"], errors="coerce").ge(pip_high)
    & (
        signals["snp_to_apa_evidence_label"].eq("direct_high_confidence")
        | (
            signals["snp_to_apa_evidence_label"].eq("local_regulatory_supported")
            & signals["local_clip_support"].map(bool_value)
        )
    )
    & signals["allele_harmonized"].map(bool_value)
    & ~signals["direct_pas_reference_conflict"].map(bool_value)
)


# --- select mechanism variant per signal ---

mechanism_order = {
    "direct_high_confidence": 0,
    "local_regulatory_supported": 1,
    "direct_exploratory": 2,
    "predicted_only": 3,
    "association_mechanism_unresolved": 4,
}
signals["mechanism_variant_priority"] = signals["snp_to_apa_evidence_label"].map(mechanism_order)
signals["signal_variant_count"] = signals.groupby("signal_id")["candidate_id"].transform("size")
signals["signal_mechanism_annotated_variant_count"] = (
    signals["direct_pas_evidence"].map(bool_value)
    | signals["local_rbp_or_structure_evidence"].map(bool_value)
    | signals["APARENT_predicted_PDUI_delta"].notna()
).groupby(signals["signal_id"]).transform("sum")
signals["_pip_sort"] = pd.to_numeric(signals["pip"], errors="coerce").fillna(-1)
signals["_representative_sort"] = ~signals["signal_representative"].map(bool_value)
signals = (
    signals.sort_values(
        ["mechanism_variant_priority", "_pip_sort", "_representative_sort"],
        ascending=[True, False, True],
        kind="stable",
    )
    .drop_duplicates("signal_id")
    .drop(columns=["_pip_sort", "_representative_sort"])
    .merge(representatives, on="signal_id", how="left")
)
if signals["signal_id"].duplicated().any():
    raise ValueError("Mechanism selection did not produce one row per independent signal")
signals = (
    signals.drop(columns=GXSTAGE_FIELDS, errors="ignore")
    .merge(interaction_signal_map, on="signal_id", how="left", validate="one_to_one")
)
has_interaction_group = signals["interaction_group_id"].notna()
same_variant = (
    has_interaction_group
    & signals["variant_id"].eq(signals["interaction_representative_variant_id"])
)
same_credible_set = (
    has_interaction_group
    & signals["cs_id"].notna()
    & signals["cs_id"].astype(str).str.strip().ne("")
    & signals["cs_id"].astype(str).eq(signals["interaction_pair_cs_id"].astype(str))
)
interaction_ld_to_selected_representative = (
    has_interaction_group
    & signals["variant_id"].eq(signals["interaction_pair_signal_lead_variant"])
    & pd.to_numeric(signals["interaction_pair_ld_r2_to_signal_lead"], errors="coerce").notna()
)
signals["gxstage_mechanism_variant_link"] = np.select(
    [same_variant, same_credible_set, interaction_ld_to_selected_representative, has_interaction_group],
    [
        "same_variant",
        "same_susie_credible_set",
        "interaction_variant_ld_to_selected_signal_representative",
        "same_stage_qtl_signal_pairwise_link_unresolved",
    ],
    default="not_interaction_group_linked",
)
signals["gxstage_mechanism_variant_link_r2"] = np.where(
    same_variant,
    1.0,
    np.where(
        interaction_ld_to_selected_representative,
        pd.to_numeric(signals["interaction_pair_ld_r2_to_signal_lead"], errors="coerce"),
        np.nan,
    ),
)


# --- detail tables for case figures ---

pas_rows = []
for row in signals.to_dict("records"):
    series = pd.Series(row)
    for label in ("proximal", "distal"):
        if pd.isna(row[f"{label}_pas_pos"]):
            continue
        pas_pos = int(row[f"{label}_pas_pos"])
        relative = transcript_relative_position(int(series["variant_pos"]), pas_pos, series["phenotype_strand"])
        window = fetch_pas_window(fasta, series["phenotype_chr"], pas_pos, series["phenotype_strand"])
        core_reference = window[50:100]
        core_alternate = core_reference
        in_core = -50 <= relative <= -1
        mutation_status = "outside_interval"
        pas_core_change = "not_eligible"
        if in_core:
            index = relative + 100
            ref = transcript_allele(series["variant_ref"], series["phenotype_strand"])
            alt = transcript_allele(series["variant_alt"], series["phenotype_strand"])
            if window[index:index + len(ref)] != ref:
                mutation_status = "reference_conflict"
                pas_core_change = "reference_conflict"
            else:
                core_alternate = (window[:index] + alt + window[index + len(ref):])[50:100]
                mutation_status = "applied"
        ref_hits = motif_hits(core_reference, motifs)
        alt_hits = motif_hits(core_alternate, motifs)
        if mutation_status == "applied":
            if sum(alt_hits.values()) > sum(ref_hits.values()):
                pas_core_change = "Gain"
            elif sum(alt_hits.values()) < sum(ref_hits.values()):
                pas_core_change = "Loss"
            elif alt_hits != ref_hits:
                pas_core_change = "Modulation"
            else:
                pas_core_change = "No core change"
        if series["phenotype_strand"] == "+":
            window_start, window_end = pas_pos - 50, pas_pos - 1
        else:
            window_start, window_end = pas_pos + 1, pas_pos + 50
        pas_rows.append({
            "candidate_id": series["candidate_id"],
            "pas_label": label,
            "pas_pos": pas_pos,
            "variant_to_pas_nt": relative,
            "variant_in_pas_core_window": in_core,
            "pas_window_start": window_start,
            "pas_window_end": window_end,
            "pas_window_ref_seq": core_reference,
            "pas_window_alt_seq": core_alternate,
            "pas_window_mutation_status": mutation_status,
            "ref_pas_motif_count": sum(ref_hits.values()),
            "alt_pas_motif_count": sum(alt_hits.values()),
            "pas_motif_count_delta": sum(alt_hits.values()) - sum(ref_hits.values()),
            "ref_pas_motifs": ";".join(sorted(ref_hits.elements())),
            "alt_pas_motifs": ";".join(sorted(alt_hits.elements())),
            "pas_core_change": pas_core_change,
            "predicted_pdui_direction": (
                "shortening"
                if (label == "proximal" and pas_core_change == "Gain")
                or (label == "distal" and pas_core_change == "Loss")
                else "lengthening"
                if pas_core_change in {"Gain", "Loss"}
                else "not_available"
            ),
        })
pas_details = pd.DataFrame(pas_rows)

selected_local = signals.loc[
    signals["local_regulatory_eligible"] & signals["allele_harmonized"].map(bool_value)
]
rbp_rows = []
for row in selected_local.to_dict("records"):
    series = pd.Series(row)
    max_width = max(motif["width"] for motif in cisbp_motifs)
    reference, alternate, offset, status = oriented_variant_sequences(series, fasta, max_width - 1)
    stage = stage_of(row)
    for motif in cisbp_motifs:
        symbol = rbp_gene_symbol(motif["rbp_name"])
        info = expression_info(symbol, rbp_expression)
        prenatal_tpm = info["Prenatal"] if info and pd.notna(info.get("Prenatal")) else pd.NA
        postnatal_tpm = info["Postnatal"] if info and pd.notna(info.get("Postnatal")) else pd.NA
        supported = rbp_expression_supported(symbol, rbp_expression, rbp_expression_min_mean_tpm, stage)
        width = motif["width"]
        start_min = max(0, offset - width + 1)
        start_max = min(offset, len(reference) - width)
        ref_best = 0.0
        alt_best = 0.0
        for start in range(start_min, start_max + 1):
            ref_best = max(ref_best, normalized_pwm_score(reference[start:start + width], motif))
            alt_best = max(alt_best, normalized_pwm_score(alternate[start:start + width], motif))
        delta = alt_best - ref_best
        name = str(motif["rbp_name"])
        if name.startswith("(") and ")" in name:
            match = re.search(r"_\(([^()]*)\)_\(RBD_[0-9.]+\)", name[name.index(")") + 1:])
            motif_source_species = match.group(1) if match else ""
        else:
            motif_source_species = "Homo_sapiens"
        rbp_rows.append({
            "candidate_id": series["candidate_id"],
            "rbp_name": motif["rbp_name"],
            "motif_id": motif["motif_id"],
            "motif_width": width,
            "variant_window_mutation_status": status,
            "ref_best_norm_score": ref_best,
            "alt_best_norm_score": alt_best,
            "delta_norm_score": delta,
            "delta_direction": "gain" if delta > 0 else "loss" if delta < 0 else "unchanged",
            "rbp_gene_symbol": symbol,
            "rbp_motif_source_species": motif_source_species,
            "rbp_expression_gene_id": info["expression_gene_id"] if info else pd.NA,
            "rbp_expression_gene_symbol": info["expression_gene_symbol"] if info else pd.NA,
            "rbp_mean_tpm_prenatal": prenatal_tpm,
            "rbp_mean_tpm_postnatal": postnatal_tpm,
            "rbp_expression_supported": supported,
        })
rbp_details = pd.DataFrame(rbp_rows)
fasta.close()
for table in (pas_details, rbp_details):
    table["analysis_version"] = ANALYSIS_VERSION
    table["pipeline_run_id"] = PIPELINE_RUN_ID
write_gzip_table(pas_details.sort_values(["candidate_id", "pas_label"]), PAS_DETAIL_OUTPUT)
write_gzip_table(rbp_details.sort_values(["candidate_id", "rbp_name", "motif_id"]), RBP_DETAIL_OUTPUT)


# --- signal mechanism labels ---

signals["mechanism_variant_selection"] = "strongest_evidence_then_pip_within_assigned_signal"
signals["aparent_role"] = "basic_APARENT_delta_distal_minus_proximal"
signals["analysis_version"] = ANALYSIS_VERSION
signals["pipeline_run_id"] = PIPELINE_RUN_ID
signals["config_path"] = str(CONFIG_PATH)
signals["genome_fasta"] = config["reference_resources"]["genome_fasta"]["path"]
signals["genome_fasta_release"] = config["reference_resources"]["genome_fasta"]["release"]
signals["cisbp_rna_release"] = config["reference_resources"]["cisbp_rna_motifs"]["release"]
signals["postar_release"] = config["reference_resources"]["postar3_clip_peaks"]["release"]
signals["rbp_expression_release"] = rbp_expression_cfg.get("release", "")
signals["rbp_expression_min_mean_tpm"] = rbp_expression_min_mean_tpm
signals["rbp_expression_stage_matched"] = bool(
    config["parameters"].get("rbp_expression", {}).get("stage_matched", True)
)
signals["pas_core_window_nt"] = "-50:-1"
signals["mechanism_formal_unit"] = "independent_signal"


# --- output ---

front = [
    "signal_id", "candidate_id", "stage", "phenotype_id", "gene_name", "variant_id",
    "representative_candidate_id", "representative_variant_id", "representative_variant_pip",
    "signal_variant_count", "signal_mechanism_annotated_variant_count", "mechanism_variant_selection",
    "pip", "finemap_tier", "credible_set_member", "signal_assignment_method",
    *GXSTAGE_FIELDS, *GXSTAGE_LINK_FIELDS,
    "proximal_pas_pos", "distal_pas_pos", "pas_anchor_resolution", "dur_start", "dur_end", "dur_length",
    "direct_pas_evidence", "proximal_pas_core_change", "distal_pas_core_change",
    "direct_pas_sites", "direct_pas_anchor_high_confidence",
    "direct_pas_predicted_pdui_direction", "mechanism_direction_concordant",
    "local_regulatory_eligible", "local_rbp_motif_change", "local_clip_support", "local_structure_change",
    "local_structure_delta_mfe", "v2_top_rbp_gain_delta", "v2_top_rbp_loss_delta",
    "postar_variant_overlap_count", "postar_variant_rbps", "postar_variant_technologies",
    "variant_phastcons", "proximal_pas_phastcons", "distal_pas_phastcons", "local_annotation_status",
    "APARENT_delta_proximal", "APARENT_delta_distal", "APARENT_predicted_PDUI_delta",
    "aparent_scoring_status", "aparent_beta_direction_consistent", "aparent_role",
    "snp_to_apa_evidence_label", "mechanism_star_eligible", "analysis_version",
]
signals = signals[front + [column for column in signals.columns if column not in front]]
signals.sort_values(["stage", "phenotype_id", "signal_id"], kind="stable").to_csv(
    OUTPUT, sep="\t", index=False, compression="gzip"
)
print(f"Written {len(signals)} assigned independent signals -> {OUTPUT}")
