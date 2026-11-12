#!/usr/bin/env python3
# Scores APA-QTL candidates for mechanism and functional effects.
# Writes the module 07 framework tables and run summary to data/07_coloc.
from __future__ import annotations

import gzip
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pysam

# --- paths ---
ROOT_DIR = Path(os.environ["APA_ROOT"])
REF_DIR = ROOT_DIR / "ref"
DATA_DIR = ROOT_DIR / "data" / "07_coloc"
INPUT_TABLE = DATA_DIR / "functional_support_annotations.tsv.gz"
GENOME_FASTA = REF_DIR / "hg38.fa"
GTF_PATH = REF_DIR / "gencode.v40.annotation.gtf.gz"
MATURE_FA = REF_DIR / "mature.fa"
PAS_MOTIF_LIST = REF_DIR / "motif.list"
CISBP_MEME = REF_DIR / "motif_databases" / "CISBP-RNA" / "Homo_sapiens.meme"
ARE_REFERENCE = REF_DIR / "Intron_ARE.csv"
TARGETSCAN_FAMILY_ZIP = REF_DIR / "Summary_Counts.all_predictions.txt.zip"
BCA_PSEUDOBULK = (ROOT_DIR / ".." / ".." / "effect" / "EXP" / "sc_ref" / "BCA_V1_pseudobulk.csv").resolve()
APARENT_DIR = Path(f"{os.environ['APA_HPC_HOME']}/bin/aparent")
SC_REF_DIR = Path(f"{os.environ['APA_GROUP_SHARE']}/.ref/sc")

# --- parameters ---
# GU-rich element (GRE) core motifs: GUUUG pentamer (GTTTG in DNA) and extended forms
# recognized by CUG-BP/CELF family RBPs; associated with mRNA instability and APA modulation
GRE_MOTIFS = {"GTTTG", "TTTGTT", "GTTTGT", "TGTTTG"}
# U-rich elements (URE): poly-U tracts that function as downstream sequence elements (DSE)
# enhancing PAS cleavage efficiency (AATAAA + DSE model)
URE_MOTIFS = {"TTTTT", "TTTTTT", "TTTTTTT", "TATTTATT", "TATTTAT"}
PAS_AUXILIARY_MOTIF_CLASSES = {
    "ARE": {"ATTTA", "TATTTAT", "AATTTA", "TTATTTATT"},
    "GRE": GRE_MOTIFS,
    "URE": URE_MOTIFS,
}
APARENT_MODEL_NAME = "aparent_large_lessdropout_all_libs_no_sampleweights.h5"
APARENT_WINDOW_HALF = 102  # symmetric 205nt window around proximal PAS (205 = 2*102 + 1)

PAS_RESOLVED_STATUSES = {
    "resolved_region_extrema",
    "resolved_single_pas_plus_distal_flank",
}
CANONICAL_PAS_HEXAMERS = {
    "AATAAA",
    "ATTAAA",
    "AGTAAA",
    "TATAAA",
    "CATAAA",
    "GATAAA",
    "AATATA",
    "AATACA",
    "AATAGA",
    "ACTAAA",
    "AAGAAA",
    "AATGAA",
}
RNAFOLD_MAX_LEN = 4000
RNAFOLD_TIMEOUT_SEC = 30
RNAFOLD_WORKERS = int(os.environ.get("APA05_RNAFOLD_WORKERS", "4"))
RBP_DELTA_THRESHOLD = 0.20
PAS_DELTA_THRESHOLD = 1
PRIORITY_CASE_LIMIT = 50
DENSITY_SUPPORT_QUANTILE = 0.75
DEFAULT_PIP_THRESHOLD = 0.5
DEFAULT_CORE_THRESHOLD = 0.8
FRAMEWORK_TIER1_PIP_THRESHOLD = DEFAULT_PIP_THRESHOLD
FRAMEWORK_HIGH_CONFIDENCE_PIP_THRESHOLD = DEFAULT_CORE_THRESHOLD
FRAMEWORK_TIER2_PIP_THRESHOLD = 0.1
TOTAL_SCORE_WEIGHTS = {
    "pip_norm": 0.2,
    "aparent_delta_norm": 0.2,
    "mirna_density_norm": 0.2,
    "rbp_density_norm": 0.2,
    "phastcons_norm": 0.2,
}
EXPLICIT_DIRECTION_COLUMNS = (
    "utr_effect_direction",
    "effect_direction",
    "direction",
)
SIGNED_EFFECT_COLUMNS = (
    "beta",
    "slope",
    "effect",
    "qtl_beta",
    "qtl_slope",
)


# --- helpers ---
@dataclass
class TranscriptRecord:
    transcript_id: str
    gene_id: str | None = None
    gene_name: str | None = None
    chrom: str | None = None
    strand: str | None = None
    transcript_start: int | None = None
    transcript_end: int | None = None
    exons: list[tuple[int, int]] = field(default_factory=list)
    cds: list[tuple[int, int]] = field(default_factory=list)
    stop_codons: list[tuple[int, int]] = field(default_factory=list)


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    else:
        df.to_csv(path, sep="\t", index=False)


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for record_id, sequence in records:
            handle.write(f">{record_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
    return seq.translate(table)[::-1]


def canonical_pas_hits(text: str) -> set[str]:
    return {
        motif
        for motif in str(text).split(";")
        if motif in CANONICAL_PAS_HEXAMERS
    }


def normalize_effect_direction(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if not text:
        return None
    mapping = {
        "lengthening": "lengthening",
        "lengthen": "lengthening",
        "longer_utr": "lengthening",
        "longer-utr": "lengthening",
        "utr_lengthening": "lengthening",
        "extension": "lengthening",
        "extended": "lengthening",
        "shortening": "shortening",
        "shorten": "shortening",
        "shorter_utr": "shortening",
        "shorter-utr": "shortening",
        "utr_shortening": "shortening",
        "contraction": "shortening",
        "shortened": "shortening",
        "unknown": "unknown",
    }
    return mapping.get(text)


def annotate_effect_directionality(df: pd.DataFrame) -> pd.DataFrame:
    annotated = df.copy()
    original_columns = set(annotated.columns)
    annotated["utr_effect_direction"] = "unknown"
    annotated["effect_direction_status"] = "missing_signed_qtl_effect"
    annotated["effect_direction_source_column"] = ""

    explicit_source: str | None = None
    for col in EXPLICIT_DIRECTION_COLUMNS:
        if col not in original_columns:
            continue
        normalized = annotated[col].map(normalize_effect_direction)
        mask = normalized.notna()
        if mask.any():
            annotated.loc[mask, "utr_effect_direction"] = normalized.loc[mask]
            annotated.loc[mask, "effect_direction_status"] = "provided_by_input_direction_column"
            annotated.loc[mask, "effect_direction_source_column"] = col
            explicit_source = col
            break

    signed_cols = [col for col in SIGNED_EFFECT_COLUMNS if col in annotated.columns]
    if signed_cols:
        signed_available = annotated[signed_cols].notna().any(axis=1)
        unresolved_mask = signed_available & annotated["effect_direction_status"].ne("provided_by_input_direction_column")
        annotated.loc[unresolved_mask, "effect_direction_status"] = (
            "signed_effect_available_but_apa_mapping_unresolved"
        )
        annotated.loc[unresolved_mask, "effect_direction_source_column"] = ",".join(signed_cols)

    if explicit_source is None and not signed_cols:
        annotated["effect_direction_source_column"] = ""
    return annotated


def nonzero_quantile(series: pd.Series, quantile: float) -> float:
    nonzero = series.loc[series.fillna(0).gt(0)]
    if nonzero.empty:
        return 0.0
    return float(nonzero.quantile(quantile))


def minmax_scale(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    lo = numeric.min()
    hi = numeric.max()
    if pd.isna(lo) or pd.isna(hi):
        return pd.Series(0.0, index=series.index, dtype=float)
    if hi == lo:
        return pd.Series(0.5, index=series.index, dtype=float)
    return ((numeric - lo) / (hi - lo)).fillna(0.0)


def resolve_optional_column(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    return next((col for col in candidates if col in df.columns), None)


def compute_reweighted_total_score(
    merged: pd.DataFrame,
    component_weights: dict[str, float],
) -> tuple[pd.Series, pd.Series]:
    component_df = merged[list(component_weights)].copy()
    availability_df = component_df.notna()
    weight_series = pd.Series(component_weights, dtype=float)
    weighted_sum = component_df.fillna(0.0).mul(weight_series, axis=1).sum(axis=1)
    available_weight = availability_df.mul(weight_series, axis=1).sum(axis=1)
    total_score = (weighted_sum / available_weight.where(available_weight.gt(0))).fillna(0.0)
    status = available_weight.map(
        lambda value: (
            "equal_weight_score_reweighted_over_available_components"
            if value > 0 and value < weight_series.sum()
            else "equal_weight_score_all_components_available"
            if value == weight_series.sum()
            else "no_scoreable_components"
        )
    )
    return total_score, status


def summarize_supported_modules(row: pd.Series, columns: list[tuple[str, str]]) -> str:
    labels = [label for label, col in columns if bool(row.get(col, False))]
    return ";".join(labels)


def fetch_genomic_seq(fasta: pysam.FastaFile, chrom: str, start: int, end: int) -> str:
    if start > end or start < 1:
        return ""
    return fasta.fetch(chrom, start - 1, end).upper()


def fetch_oriented_seq(
    fasta: pysam.FastaFile, chrom: str, start: int, end: int, strand: str
) -> str:
    seq = fetch_genomic_seq(fasta, chrom, start, end)
    return seq if strand == "+" else revcomp(seq)


def apply_variant_to_oriented_seq(
    sequence: str,
    interval_start: int,
    interval_end: int,
    strand: str,
    variant_pos: int,
    ref: str,
    alt: str,
) -> tuple[str, str]:
    if len(ref) != 1 or len(alt) != 1:
        return sequence, "non_snv"
    if variant_pos < interval_start or variant_pos > interval_end:
        return sequence, "outside_interval"
    seq_list = list(sequence)
    offset = variant_pos - interval_start if strand == "+" else interval_end - variant_pos
    ref_tx = ref.upper() if strand == "+" else revcomp(ref.upper())
    alt_tx = alt.upper() if strand == "+" else revcomp(alt.upper())
    if offset < 0 or offset >= len(seq_list):
        return sequence, "offset_error"
    if seq_list[offset] != ref_tx:
        return sequence, "ref_mismatch"
    seq_list[offset] = alt_tx
    return "".join(seq_list), "applied"


def parse_cisbp_meme(path: Path) -> list[dict]:
    motifs: list[dict] = []
    current_name: str | None = None
    current_id: str | None = None
    current_width: int | None = None
    rows: list[list[float]] = []

    def flush_current() -> None:
        if current_name is None or current_id is None or current_width is None or not rows:
            return
        pwm = rows[:]
        max_score = sum(math.log2(max(max(row), 1e-6) / 0.25) for row in pwm)
        min_score = sum(math.log2(max(min(row), 1e-6) / 0.25) for row in pwm)
        consensus = "".join("ACGU"[max(range(4), key=lambda idx: row[idx])] for row in pwm)
        motifs.append(
            {
                "motif_id": current_id,
                "rbp_name": current_name,
                "width": current_width,
                "pwm": pwm,
                "max_score": max_score,
                "min_score": min_score,
                "consensus_rna": consensus,
                "consensus_dna": consensus.replace("U", "T"),
            }
        )

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        in_matrix = False
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                if in_matrix and current_width is not None and len(rows) >= current_width:
                    flush_current()
                    current_name = None
                    current_id = None
                    current_width = None
                    rows = []
                    in_matrix = False
                continue
            if line.startswith("MOTIF "):
                if rows:
                    flush_current()
                    rows = []
                    current_width = None
                    in_matrix = False
                parts = line.split()
                current_id = parts[1]
                current_name = parts[2]
                continue
            if line.startswith("letter-probability matrix"):
                match = re.search(r"w=\s*(\d+)", line)
                if match is None:
                    raise ValueError(f"Could not parse motif width from: {line}")
                current_width = int(match.group(1))
                rows = []
                in_matrix = True
                continue
            if in_matrix and current_width is not None:
                parts = line.split()
                if len(parts) == 4:
                    rows.append([float(value) for value in parts])
                    if len(rows) == current_width:
                        flush_current()
                        current_name = None
                        current_id = None
                        current_width = None
                        rows = []
                        in_matrix = False
    return motifs


def normalized_pwm_score(seq_dna: str, motif: dict) -> float:
    if len(seq_dna) != motif["width"]:
        return 0.0
    seq_rna = seq_dna.replace("T", "U")
    index = {"A": 0, "C": 1, "G": 2, "U": 3}
    score = 0.0
    for base, row in zip(seq_rna, motif["pwm"], strict=True):
        if base not in index:
            return 0.0
        score += math.log2(max(row[index[base]], 1e-6) / 0.25)
    denom = motif["max_score"] - motif["min_score"]
    if denom <= 0:
        return 0.0
    return (score - motif["min_score"]) / denom


def overlapping_count(sequence: str, pattern: str) -> int:
    count = 0
    start = 0
    while True:
        idx = sequence.find(pattern, start)
        if idx == -1:
            return count
        count += 1
        start = idx + 1


def parse_human_mirnas(path: Path) -> list[dict]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None and chunks:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.upper().replace("U", "T"))
        if header is not None and chunks:
            records.append((header, "".join(chunks)))

    families: dict[str, dict] = {}
    for header_text, seq in records:
        fields = header_text.split()
        mirna_name = fields[0]
        if not mirna_name.startswith("hsa-") or len(seq) < 8:
            continue
        seed = seq[1:8]
        entry = families.setdefault(
            seed,
            {
                "seed": seed,
                "family_id": f"seed_{seed}",
                "representative_mirna": mirna_name,
                "members": [],
            },
        )
        entry["members"].append(mirna_name)
    for entry in families.values():
        seed_rc = revcomp(entry["seed"])
        sixmer_rc = revcomp(entry["seed"][:6])
        entry["site_8mer"] = "A" + seed_rc
        entry["site_7mer_m8"] = seed_rc
        entry["site_7mer_a1"] = "A" + sixmer_rc
        entry["members"] = sorted(set(entry["members"]))
    return sorted(families.values(), key=lambda item: item["representative_mirna"])


def intersect_segments(segments: list[tuple[int, int]], start: int, end: int) -> list[tuple[int, int]]:
    if start > end:
        return []
    out: list[tuple[int, int]] = []
    for seg_start, seg_end in segments:
        left = max(seg_start, start)
        right = min(seg_end, end)
        if left <= right:
            out.append((left, right))
    return out


def transcript_segments_to_seq(
    fasta: pysam.FastaFile,
    chrom: str,
    strand: str,
    segments: list[tuple[int, int]],
) -> str:
    if not segments:
        return ""
    pieces = [fetch_genomic_seq(fasta, chrom, start, end) for start, end in sorted(segments)]
    sequence = "".join(pieces)
    return sequence if strand == "+" else revcomp(sequence)


def scan_pas_sequence_hits(sequence: str, motifs: list[str]) -> tuple[int, str]:
    found = [motif for motif in motifs if motif in sequence]
    return len(found), ";".join(found[:10])


def load_optional_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    sep = "\t" if path.suffix in {".tsv", ".gz"} else ","
    return pd.read_csv(path, sep=sep)


def run_one_rnafold(record: tuple[str, str]) -> dict:
    record_id, seq = record
    with NamedTemporaryFile("w", delete=False, suffix=".fa", encoding="utf-8") as handle:
        temp_path = Path(handle.name)
        handle.write(f">{record_id}\n{seq}\n")
    try:
        completed = subprocess.run(
            ["RNAfold", "--noPS", str(temp_path)],
            check=True,
            capture_output=True,
            text=True,
            timeout=RNAFOLD_TIMEOUT_SEC,
        )
    except subprocess.TimeoutExpired:
        return {
            "record_id": record_id,
            "seq_length": len(seq),
            "mfe_kcal_per_mol": math.nan,
            "mfe_per_100nt": math.nan,
            "rnafold_status": "timeout",
        }
    except subprocess.CalledProcessError:
        return {
            "record_id": record_id,
            "seq_length": len(seq),
            "mfe_kcal_per_mol": math.nan,
            "mfe_per_100nt": math.nan,
            "rnafold_status": "failed",
        }
    finally:
        temp_path.unlink(missing_ok=True)

    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 3:
        return {
            "record_id": record_id,
            "seq_length": len(seq),
            "mfe_kcal_per_mol": math.nan,
            "mfe_per_100nt": math.nan,
            "rnafold_status": "parse_failed",
        }
    folded_record_id = lines[0].lstrip(">")
    folded_seq = lines[1]
    struct_line = lines[2]
    match = re.search(r"\(\s*([-+]?[0-9]*\.?[0-9]+)\)$", struct_line)
    mfe = float(match.group(1)) if match else math.nan
    return {
        "record_id": folded_record_id,
        "seq_length": len(folded_seq),
        "mfe_kcal_per_mol": mfe,
        "mfe_per_100nt": (mfe / len(folded_seq) * 100.0) if folded_seq else math.nan,
        "rnafold_status": "ok" if match else "parse_failed",
    }


if __name__ == "__main__":
    # --- candidates ---
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    base_df = pd.read_csv(INPUT_TABLE, sep="\t")
    mechanism_df = base_df.copy()
    fate_df = base_df.loc[base_df["dur_length"].notna() & base_df["pas_resolution_status"].isin(PAS_RESOLVED_STATUSES)].copy()
    mechanism_df["analysis_set"] = "mechanism"
    fate_df["analysis_set"] = "functional_effect"
    mechanism_sets_df = pd.concat([mechanism_df, fate_df], ignore_index=True, sort=False)[
        ["candidate_id", "gene_name", "stage", "tier", "pas_resolution_status", "dur_length", "analysis_set"]
    ].drop_duplicates()
    write_table(mechanism_sets_df, DATA_DIR / "candidate_analysis_sets.tsv")

    # --- references ---
    targetscan_rows = []
    for path in [REF_DIR / "Summary_Counts.all_predictions.txt.zip", REF_DIR / "Conserved_Family_Info.txt.zip"]:
        with ZipFile(path) as zf:
            member = zf.namelist()[0]
            targetscan_rows.append({"reference_file": path.name, "member": member})
    targetscan_manifest = pd.DataFrame(targetscan_rows)
    write_table(targetscan_manifest, DATA_DIR / "targetscan_reference_manifest.tsv")

    cache_dir = Path("/tmp/apa_07_coloc_fasta_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    local_fasta = cache_dir / "hg38.fa"
    if local_fasta.exists() or local_fasta.is_symlink():
        target = local_fasta.resolve(strict=False)
        if target != GENOME_FASTA.resolve():
            local_fasta.unlink()
    if not local_fasta.exists():
        local_fasta.symlink_to(GENOME_FASTA)
    fai_path = local_fasta.with_suffix(local_fasta.suffix + ".fai")
    if not fai_path.exists():
        subprocess.run(["samtools", "faidx", str(local_fasta)], check=True)
    fasta = pysam.FastaFile(str(local_fasta))

    pas_motifs = re.findall(r'"([ACGT]+)"', PAS_MOTIF_LIST.read_text(encoding="utf-8"))
    cisbp_motifs = parse_cisbp_meme(CISBP_MEME)
    targetscan_families = []
    if TARGETSCAN_FAMILY_ZIP.exists():
        with ZipFile(TARGETSCAN_FAMILY_ZIP) as zf:
            member = zf.namelist()[0]
            with zf.open(member) as handle:
                header = pd.read_csv(handle, sep="\t", nrows=0)
        columns = list(header.columns)
        seed_col = next(
            (
                col
                for col in columns
                if col.lower() in {"mirna family", "mirna_family", "seed", "seed sequence"}
                or "seed" in col.lower()
            ),
            None,
        )
        family_col = next(
            (col for col in columns if col.lower() in {"representative mirna", "representative_mirna"}),
            None,
        )
        species_col = next((col for col in columns if "species" in col.lower()), None)
        if seed_col is not None:
            usecols = [seed_col]
            if family_col is not None:
                usecols.append(family_col)
            if species_col is not None:
                usecols.append(species_col)
            families: dict[str, dict] = {}
            with ZipFile(TARGETSCAN_FAMILY_ZIP) as zf:
                member = zf.namelist()[0]
                with zf.open(member) as handle:
                    for family_df in pd.read_csv(handle, sep="\t", usecols=usecols, chunksize=200000):
                        if species_col is not None:
                            species_text = family_df[species_col].astype(str).str.lower()
                            human_mask = species_text.str.contains("human|hsa|9606", regex=True, na=False)
                            if human_mask.any():
                                family_df = family_df.loc[human_mask].copy()
                        for _, row in family_df.iterrows():
                            seed = str(row[seed_col]).upper().replace("U", "T").strip()
                            if len(seed) < 6 or set(seed) - {"A", "C", "G", "T"}:
                                continue
                            seed7 = seed[:7] if len(seed) >= 7 else seed
                            family_id = str(row[family_col]).strip() if family_col is not None else f"TargetScan_{seed7}"
                            entry = families.setdefault(
                                seed7,
                                {
                                    "seed": seed7,
                                    "family_id": family_id,
                                    "representative_mirna": family_id,
                                    "members": [],
                                },
                            )
                            entry["members"].append(family_id)
            for entry in families.values():
                seed_rc = revcomp(entry["seed"])
                sixmer_rc = revcomp(entry["seed"][:6])
                entry["site_8mer"] = "A" + seed_rc
                entry["site_7mer_m8"] = seed_rc
                entry["site_7mer_a1"] = "A" + sixmer_rc
                entry["members"] = sorted(set(entry["members"]))
            targetscan_families = sorted(families.values(), key=lambda item: item["representative_mirna"])
    human_mirna_families = targetscan_families or parse_human_mirnas(MATURE_FA)

    # --- ARE reference ---
    df = pd.read_csv(ARE_REFERENCE)
    patterns = {
        str(pattern).upper()
        for pattern in df["Pattern"].dropna().astype(str).tolist()
        if set(str(pattern).upper()) <= {"A", "C", "G", "T"}
    }
    patterns.update({"ATTTA", "TATTTAT", "AATTTA", "TTATTTATT"})
    are_patterns = patterns

    # --- transcripts ---
    transcript_ids = set(fate_df["transcript_id"].unique())
    transcripts: dict[str, TranscriptRecord] = {}
    with gzip.open(GTF_PATH, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            parts = raw_line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr_text = parts
            attrs: dict[str, str] = {}
            for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr_text):
                attrs[match.group(1)] = match.group(2)
            transcript_id = attrs.get("transcript_id")
            if transcript_id not in transcript_ids:
                continue
            record = transcripts.setdefault(transcript_id, TranscriptRecord(transcript_id=transcript_id))
            record.chrom = chrom
            record.strand = strand
            record.gene_id = attrs.get("gene_id", record.gene_id)
            record.gene_name = attrs.get("gene_name", record.gene_name)
            start_i = int(start)
            end_i = int(end)
            if feature == "transcript":
                record.transcript_start = start_i
                record.transcript_end = end_i
            elif feature == "exon":
                record.exons.append((start_i, end_i))
            elif feature == "CDS":
                record.cds.append((start_i, end_i))
            elif feature == "stop_codon":
                record.stop_codons.append((start_i, end_i))
    for record in transcripts.values():
        record.exons.sort()
        record.cds.sort()
        record.stop_codons.sort()
        if record.transcript_start is None and record.exons:
            record.transcript_start = min(start for start, _ in record.exons)
            record.transcript_end = max(end for _, end in record.exons)

    # --- PAS and RBP motif scans ---
    rows: list[dict] = []
    for row in mechanism_df.itertuples(index=False):
        for pas_label, pas_pos in (("proximal", row.proximal_pas_pos), ("distal", row.distal_pas_pos)):
            if pd.isna(pas_pos):
                continue
            pas_pos = int(pas_pos)
            if row.phenotype_strand == "+":
                window_start = max(1, pas_pos - 40)
                window_end = pas_pos - 1
            else:
                window_start = pas_pos + 1
                window_end = pas_pos + 40
            if window_start > window_end:
                continue
            ref_seq = fetch_oriented_seq(fasta, row.phenotype_chr, window_start, window_end, row.phenotype_strand)
            alt_seq, mutate_status = apply_variant_to_oriented_seq(
                ref_seq,
                window_start,
                window_end,
                row.phenotype_strand,
                int(row.variant_pos),
                row.ref_allele,
                row.alt_allele,
            )
            ref_hit_count, ref_hits = scan_pas_sequence_hits(ref_seq, pas_motifs)
            alt_hit_count, alt_hits = scan_pas_sequence_hits(alt_seq, pas_motifs)
            rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "pas_label": pas_label,
                    "pas_pos": pas_pos,
                    "variant_to_pas_nt": int(abs(int(row.variant_pos) - pas_pos)),
                    "variant_in_pas_core_window": window_start <= int(row.variant_pos) <= window_end,
                    "pas_window_start": window_start,
                    "pas_window_end": window_end,
                    "pas_window_ref_seq": ref_seq,
                    "pas_window_alt_seq": alt_seq,
                    "pas_window_mutation_status": mutate_status,
                    "ref_pas_motif_count": ref_hit_count,
                    "alt_pas_motif_count": alt_hit_count,
                    "pas_motif_count_delta": alt_hit_count - ref_hit_count,
                    "ref_pas_motifs": ref_hits,
                    "alt_pas_motifs": alt_hits,
                }
            )
    pas_delta_df = pd.DataFrame(rows)
    write_table(pas_delta_df, DATA_DIR / "pas_motif_delta.tsv.gz")

    rows: list[dict] = []
    for candidate_id, group in pas_delta_df.groupby("candidate_id", sort=False):
        annotated = group.copy()
        annotated["ref_pas_canonical_hits"] = annotated["ref_pas_motifs"].map(
            lambda value: ";".join(sorted(canonical_pas_hits(value)))
        )
        annotated["alt_pas_canonical_hits"] = annotated["alt_pas_motifs"].map(
            lambda value: ";".join(sorted(canonical_pas_hits(value)))
        )
        annotated["pas_canonical_identity_changed"] = (
            annotated["ref_pas_canonical_hits"] != annotated["alt_pas_canonical_hits"]
        )
        annotated["variant_in_pas_core_window"] = (
            annotated["variant_in_pas_core_window"].astype("boolean").fillna(False).astype(bool)
        )
        annotated["pas_count_change_support"] = annotated["pas_motif_count_delta"].abs().ge(PAS_DELTA_THRESHOLD)
        annotated["pas_canonical_rescue"] = (
            annotated["variant_in_pas_core_window"]
            & annotated["pas_canonical_identity_changed"]
        )
        annotated["pas_direct_support"] = (
            annotated["pas_count_change_support"]
            | annotated["pas_canonical_rescue"]
        )
        annotated["pas_direct_support_source"] = annotated.apply(
            lambda row: (
                "motif_count_delta+canonical_identity_rescue"
                if row["pas_count_change_support"] and row["pas_canonical_rescue"]
                else "motif_count_delta"
                if row["pas_count_change_support"]
                else "canonical_identity_rescue"
                if row["pas_canonical_rescue"]
                else "not_supported"
            ),
            axis=1,
        )
        annotated["pas_rescue_status"] = annotated.apply(
            lambda row: (
                "rescued_by_canonical_identity_change"
                if row["pas_canonical_rescue"]
                else "supported_by_motif_count_delta"
                if row["pas_count_change_support"]
                else "identity_changed_outside_pas_core"
                if row["pas_canonical_identity_changed"]
                else "not_supported"
            ),
            axis=1,
        )
        annotated["pas_delta_abs"] = annotated["pas_motif_count_delta"].abs()
        best = annotated.sort_values(
            [
                "pas_direct_support",
                "pas_delta_abs",
                "pas_canonical_identity_changed",
                "variant_in_pas_core_window",
            ],
            ascending=[False, False, False, False],
            kind="stable",
        ).iloc[0]
        rows.append(
            {
                "candidate_id": candidate_id,
                "pas_direct_support": bool(best["pas_direct_support"]),
                "top_pas_support_label": best["pas_label"],
                "top_pas_motif_delta": int(best["pas_motif_count_delta"]),
                "top_pas_variant_distance_nt": int(best["variant_to_pas_nt"]),
                "top_pas_mutation_status": best["pas_window_mutation_status"],
                "ref_pas_canonical_hits": best["ref_pas_canonical_hits"],
                "alt_pas_canonical_hits": best["alt_pas_canonical_hits"],
                "pas_canonical_identity_changed": bool(best["pas_canonical_identity_changed"]),
                "pas_direct_support_source": best["pas_direct_support_source"],
                "pas_rescue_status": best["pas_rescue_status"],
            }
        )
    pas_summary_df = pd.DataFrame(rows)

    rows: list[dict] = []
    for row in mechanism_df.itertuples(index=False):
        for pas_label, pas_pos in (("proximal", row.proximal_pas_pos), ("distal", row.distal_pas_pos)):
            if pd.isna(pas_pos):
                continue
            pas_pos = int(pas_pos)
            if row.phenotype_strand == "+":
                windows = [
                    ("USE", max(1, pas_pos - 100), max(1, pas_pos - 40)),
                    ("DSE", pas_pos + 10, pas_pos + 50),
                ]
            else:
                windows = [
                    ("USE", pas_pos + 40, pas_pos + 100),
                    ("DSE", max(1, pas_pos - 50), max(1, pas_pos - 10)),
                ]
            for auxiliary_region, window_start, window_end in windows:
                if window_start > window_end:
                    window_start, window_end = window_end, window_start
                ref_seq = fetch_oriented_seq(
                    fasta,
                    row.phenotype_chr,
                    window_start,
                    window_end,
                    row.phenotype_strand,
                )
                alt_seq, mutate_status = apply_variant_to_oriented_seq(
                    ref_seq,
                    window_start,
                    window_end,
                    row.phenotype_strand,
                    int(row.variant_pos),
                    row.ref_allele,
                    row.alt_allele,
                )
                for motif_class, motifs in PAS_AUXILIARY_MOTIF_CLASSES.items():
                    ref_count = int(sum(overlapping_count(ref_seq, motif) for motif in motifs))
                    alt_count = int(sum(overlapping_count(alt_seq, motif) for motif in motifs))
                    rows.append(
                        {
                            "candidate_id": row.candidate_id,
                            "pas_label": pas_label,
                            "auxiliary_region": auxiliary_region,
                            "motif_class": motif_class,
                            "pas_pos": pas_pos,
                            "variant_to_pas_nt": int(abs(int(row.variant_pos) - pas_pos)),
                            "variant_in_auxiliary_window": window_start <= int(row.variant_pos) <= window_end,
                            "window_start": window_start,
                            "window_end": window_end,
                            "window_mutation_status": mutate_status,
                            "ref_auxiliary_motif_count": ref_count,
                            "alt_auxiliary_motif_count": alt_count,
                            "auxiliary_motif_count_delta": alt_count - ref_count,
                        }
                    )
    pas_auxiliary_delta_df = pd.DataFrame(rows)
    write_table(pas_auxiliary_delta_df, DATA_DIR / "pas_auxiliary_motif_delta.tsv.gz")

    if pas_auxiliary_delta_df.empty:
        pas_auxiliary_summary_df = pd.DataFrame(columns=[
            "candidate_id",
            "pas_auxiliary_motif_support",
            "top_pas_auxiliary_motif_class",
            "top_pas_auxiliary_motif_delta",
        ])
    else:
        rows: list[dict] = []
        for candidate_id, group in pas_auxiliary_delta_df.groupby("candidate_id", sort=False):
            annotated = group.copy()
            annotated["delta_abs"] = annotated["auxiliary_motif_count_delta"].abs()
            annotated["auxiliary_support"] = (
                annotated["variant_in_auxiliary_window"].astype("boolean").fillna(False).astype(bool)
                & annotated["auxiliary_motif_count_delta"].ne(0)
            )
            best = annotated.sort_values(
                ["auxiliary_support", "delta_abs", "variant_to_pas_nt"],
                ascending=[False, False, True],
                kind="stable",
            ).iloc[0]
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "pas_auxiliary_motif_support": bool(best["auxiliary_support"]),
                    "top_pas_auxiliary_motif_class": best["motif_class"],
                    "top_pas_auxiliary_motif_delta": int(best["auxiliary_motif_count_delta"]),
                }
            )
        pas_auxiliary_summary_df = pd.DataFrame(rows)

    rbp_delta_df = pd.DataFrame()
    if cisbp_motifs:
        max_width = max(motif["width"] for motif in cisbp_motifs)
        rows: list[dict] = []
        for row in mechanism_df.itertuples(index=False):
            flank = max_width - 1
            start = max(1, int(row.variant_pos) - flank)
            end = int(row.variant_pos) + flank
            ref_seq = fetch_oriented_seq(fasta, row.variant_chr, start, end, row.phenotype_strand)
            alt_seq, mutate_status = apply_variant_to_oriented_seq(
                ref_seq,
                start,
                end,
                row.phenotype_strand,
                int(row.variant_pos),
                row.ref_allele,
                row.alt_allele,
            )
            variant_idx = int(row.variant_pos) - start if row.phenotype_strand == "+" else end - int(row.variant_pos)
            for motif in cisbp_motifs:
                width = motif["width"]
                start_idx_min = max(0, variant_idx - width + 1)
                start_idx_max = min(variant_idx, len(ref_seq) - width)
                best_ref = 0.0
                best_alt = 0.0
                for start_idx in range(start_idx_min, start_idx_max + 1):
                    ref_sub = ref_seq[start_idx:start_idx + width]
                    alt_sub = alt_seq[start_idx:start_idx + width]
                    best_ref = max(best_ref, normalized_pwm_score(ref_sub, motif))
                    best_alt = max(best_alt, normalized_pwm_score(alt_sub, motif))
                delta = best_alt - best_ref
                if abs(delta) < 0.01 and max(best_ref, best_alt) < 0.60:
                    continue
                rows.append(
                    {
                        "candidate_id": row.candidate_id,
                        "rbp_name": motif["rbp_name"],
                        "motif_id": motif["motif_id"],
                        "motif_width": width,
                        "variant_window_mutation_status": mutate_status,
                        "ref_best_norm_score": round(best_ref, 6),
                        "alt_best_norm_score": round(best_alt, 6),
                        "delta_norm_score": round(delta, 6),
                        "delta_direction": "gain" if delta > 0 else "loss" if delta < 0 else "no_change",
                    }
                )
        if not rows:
            rbp_delta_df = pd.DataFrame(columns=[
                "candidate_id",
                "rbp_name",
                "motif_id",
                "motif_width",
                "variant_window_mutation_status",
                "ref_best_norm_score",
                "alt_best_norm_score",
                "delta_norm_score",
                "delta_direction",
            ])
        else:
            rbp_delta_df = pd.DataFrame(rows)
    write_table(rbp_delta_df, DATA_DIR / "variant_rbp_motif_delta.tsv.gz")

    rows: list[dict] = []
    if rbp_delta_df.empty:
        rbp_summary_df = pd.DataFrame(columns=[
            "candidate_id",
            "top_rbp_gain",
            "top_rbp_gain_delta",
            "top_rbp_loss",
            "top_rbp_loss_delta",
            "rbp_variant_delta_support",
        ])
    else:
        for candidate_id, group in rbp_delta_df.groupby("candidate_id", sort=False):
            gain_row = group.sort_values("delta_norm_score", ascending=False, kind="stable").iloc[0]
            loss_row = group.sort_values("delta_norm_score", ascending=True, kind="stable").iloc[0]
            support = max(abs(float(gain_row["delta_norm_score"])), abs(float(loss_row["delta_norm_score"]))) >= RBP_DELTA_THRESHOLD
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "top_rbp_gain": gain_row["rbp_name"],
                    "top_rbp_gain_delta": gain_row["delta_norm_score"],
                    "top_rbp_loss": loss_row["rbp_name"],
                    "top_rbp_loss_delta": loss_row["delta_norm_score"],
                    "rbp_variant_delta_support": support,
                }
            )
        rbp_summary_df = pd.DataFrame(rows)
    # --- UTR sequences ---
    utr_records = []
    for _, row in fate_df.iterrows():
        transcript = transcripts.get(row["transcript_id"])
        record = {
            "candidate_id": row["candidate_id"],
            "utr_status": "ok",
            "utr_detail": "coding transcript with PAS-resolved 3UTR segments",
            "gene_id": pd.NA,
            "transcript_chr": pd.NA,
            "transcript_strand": pd.NA,
            "utr_anchor_pos": pd.NA,
            "short_utr_length": pd.NA,
            "long_utr_length": pd.NA,
            "dur_seq_length": pd.NA,
            "short_utr_seq": pd.NA,
            "long_utr_seq": pd.NA,
            "dur_seq": pd.NA,
        }
        if transcript is None:
            record["utr_status"] = "missing_transcript"
            record["utr_detail"] = "Transcript ID was not found in the supplied GTF."
            utr_records.append(record)
            continue
        record["gene_id"] = transcript.gene_id
        record["transcript_chr"] = transcript.chrom
        record["transcript_strand"] = transcript.strand
        if transcript.chrom != row["phenotype_chr"] or transcript.strand != row["phenotype_strand"]:
            record["utr_status"] = "transcript_mismatch"
            record["utr_detail"] = "Transcript chromosome or strand did not match the candidate table."
            utr_records.append(record)
            continue
        if not transcript.exons:
            record["utr_status"] = "missing_exons"
            record["utr_detail"] = "Transcript record had no exon features."
            utr_records.append(record)
            continue
        if not transcript.cds and not transcript.stop_codons:
            record["utr_status"] = "noncoding_or_missing_cds"
            record["utr_detail"] = "Transcript lacks CDS/stop_codon, so transcript-aware 3UTR cannot be defined."
            utr_records.append(record)
            continue

        if transcript.strand == "+":
            utr_anchor = max(end for _, end in transcript.stop_codons) if transcript.stop_codons else max(end for _, end in transcript.cds)
            utr_start = utr_anchor + 1
            utr_end = transcript.transcript_end or max(end for _, end in transcript.exons)
            short_region = (utr_start, int(row["proximal_pas_pos"]))
            long_region = (utr_start, int(row["distal_pas_pos"]))
            dur_region = (int(row["proximal_pas_pos"]) + 1, int(row["distal_pas_pos"]))
        else:
            utr_anchor = min(start for start, _ in transcript.stop_codons) if transcript.stop_codons else min(start for start, _ in transcript.cds)
            utr_start = transcript.transcript_start or min(start for start, _ in transcript.exons)
            utr_end = utr_anchor - 1
            short_region = (int(row["proximal_pas_pos"]), utr_end)
            long_region = (int(row["distal_pas_pos"]), utr_end)
            dur_region = (int(row["distal_pas_pos"]), int(row["proximal_pas_pos"]) - 1)

        record["utr_anchor_pos"] = int(utr_anchor)
        if utr_start > utr_end:
            record["utr_status"] = "empty_utr"
            record["utr_detail"] = "Transcript-aware 3UTR interval is empty after CDS/stop_codon anchoring."
            utr_records.append(record)
            continue

        short_segments = intersect_segments(transcript.exons, short_region[0], short_region[1])
        long_segments = intersect_segments(transcript.exons, long_region[0], long_region[1])
        dur_segments = intersect_segments(transcript.exons, dur_region[0], dur_region[1])
        if not short_segments or not long_segments:
            # Genomic fallback: PAS may be just outside annotated exon boundaries (PolyA_DB vs GTF discrepancy).
            # Use full genomic interval without exon intersection and mark accordingly.
            try:
                short_seq = fetch_genomic_seq(fasta, transcript.chrom, short_region[0], short_region[1])
                long_seq = fetch_genomic_seq(fasta, transcript.chrom, long_region[0], long_region[1])
                dur_seq = fetch_genomic_seq(fasta, transcript.chrom, dur_region[0], dur_region[1]) if dur_region[0] <= dur_region[1] else ""
                if transcript.strand == "-":
                    short_seq, long_seq, dur_seq = revcomp(short_seq), revcomp(long_seq), revcomp(dur_seq)
            except Exception:
                short_seq = long_seq = dur_seq = ""
            if short_seq and long_seq and len(long_seq) >= len(short_seq):
                record["utr_status"] = "ok_genomic_fallback"
                record["utr_detail"] = "PAS outside exonic UTR; full genomic interval used as fallback."
                record["short_utr_length"] = len(short_seq)
                record["long_utr_length"] = len(long_seq)
                record["dur_seq_length"] = len(dur_seq)
                record["short_utr_seq"] = short_seq
                record["long_utr_seq"] = long_seq
                record["dur_seq"] = dur_seq
                utr_records.append(record)
                continue
            record["utr_status"] = "pas_outside_exonic_utr"
            record["utr_detail"] = "Resolved PAS positions did not map cleanly onto transcript 3UTR exons."
            utr_records.append(record)
            continue

        short_seq = transcript_segments_to_seq(fasta, transcript.chrom, transcript.strand, short_segments)
        long_seq = transcript_segments_to_seq(fasta, transcript.chrom, transcript.strand, long_segments)
        dur_seq = transcript_segments_to_seq(fasta, transcript.chrom, transcript.strand, dur_segments)

        record["short_utr_length"] = len(short_seq)
        record["long_utr_length"] = len(long_seq)
        record["dur_seq_length"] = len(dur_seq)
        record["short_utr_seq"] = short_seq
        record["long_utr_seq"] = long_seq
        record["dur_seq"] = dur_seq
        if not short_seq or not long_seq:
            record["utr_status"] = "empty_sequence"
            record["utr_detail"] = "Sequence extraction returned empty short or long 3UTR sequences."
        elif len(long_seq) < len(short_seq):
            record["utr_status"] = "length_inversion"
            record["utr_detail"] = "Long 3UTR sequence was shorter than the short 3UTR sequence."
        utr_records.append(record)
    utr_df = pd.DataFrame(utr_records)
    write_table(
        utr_df.drop(columns=["short_utr_seq", "long_utr_seq", "dur_seq"]),
        DATA_DIR / "utr_sequence_annotations.tsv.gz",
    )

    short_fasta_records = []
    long_fasta_records = []
    dur_fasta_records = []
    for row in utr_df.itertuples(index=False):
        if isinstance(row.short_utr_seq, str) and row.short_utr_seq:
            short_fasta_records.append((row.candidate_id, row.short_utr_seq))
        if isinstance(row.long_utr_seq, str) and row.long_utr_seq:
            long_fasta_records.append((row.candidate_id, row.long_utr_seq))
        if isinstance(row.dur_seq, str) and row.dur_seq:
            dur_fasta_records.append((row.candidate_id, row.dur_seq))
    write_fasta(short_fasta_records, DATA_DIR / "short_utr.fa")
    write_fasta(long_fasta_records, DATA_DIR / "long_utr.fa")
    write_fasta(dur_fasta_records, DATA_DIR / "dur.fa")
    # --- miRNA, dur RBP and ARE scans ---
    mirna_hit_tables = []
    mirna_summary_rows = []
    dur_rbp_rows = []
    are_rows = []
    for row in utr_df.itertuples(index=False):
        if not isinstance(row.dur_seq, str) or not row.dur_seq:
            mirna_summary_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "mirna_total_sites": 0,
                    "mirna_unique_families": 0,
                    "mirna_8mer_sites": 0,
                    "mirna_7mer_m8_sites": 0,
                    "mirna_7mer_a1_sites": 0,
                    "top_mirna_families": "",
                }
            )
            dur_rbp_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "dur_rbp_consensus_total_hits": 0,
                    "dur_rbp_consensus_unique_rbps": 0,
                    "top_dur_rbp_consensus_hits": "",
                }
            )
            are_rows.append(
                {
                    "candidate_id": row.candidate_id,
                    "are_total_hits": 0,
                    "are_attTA_hits": 0,
                    "are_tatttat_hits": 0,
                    "top_are_patterns": "",
                }
            )
            continue
        sequence = row.dur_seq
        candidate_id = row.candidate_id
        hit_rows: list[dict] = []
        family_counter: Counter[str] = Counter()
        site_counter: Counter[str] = Counter()
        for family in human_mirna_families:
            hits_8mer = overlapping_count(sequence, family["site_8mer"])
            hits_7m8 = overlapping_count(sequence, family["site_7mer_m8"])
            hits_7a1 = overlapping_count(sequence, family["site_7mer_a1"])
            total_hits = hits_8mer + hits_7m8 + hits_7a1
            if total_hits == 0:
                continue
            if hits_8mer:
                hit_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "family_id": family["family_id"],
                        "representative_mirna": family["representative_mirna"],
                        "site_type": "8mer",
                        "site_count": hits_8mer,
                    }
                )
            if hits_7m8:
                hit_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "family_id": family["family_id"],
                        "representative_mirna": family["representative_mirna"],
                        "site_type": "7mer-m8",
                        "site_count": hits_7m8,
                    }
                )
            if hits_7a1:
                hit_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "family_id": family["family_id"],
                        "representative_mirna": family["representative_mirna"],
                        "site_type": "7mer-a1",
                        "site_count": hits_7a1,
                    }
                )
            family_counter[family["representative_mirna"]] += total_hits
            site_counter["8mer"] += hits_8mer
            site_counter["7mer-m8"] += hits_7m8
            site_counter["7mer-a1"] += hits_7a1

        top_families = ";".join(
            f"{name}:{count}" for name, count in family_counter.most_common(10)
        )
        mirna_summary = {
            "candidate_id": candidate_id,
            "mirna_total_sites": int(sum(family_counter.values())),
            "mirna_unique_families": int(len(family_counter)),
            "mirna_8mer_sites": int(site_counter["8mer"]),
            "mirna_7mer_m8_sites": int(site_counter["7mer-m8"]),
            "mirna_7mer_a1_sites": int(site_counter["7mer-a1"]),
            "top_mirna_families": top_families,
        }
        hit_df = pd.DataFrame(hit_rows)
        if not hit_df.empty:
            mirna_hit_tables.append(hit_df)
        mirna_summary_rows.append(mirna_summary)
        motif_counter: Counter[str] = Counter()
        for motif in cisbp_motifs:
            count = overlapping_count(sequence, motif["consensus_dna"])
            if count > 0:
                motif_counter[motif["rbp_name"]] += count
        dur_rbp_rows.append({
            "candidate_id": candidate_id,
            "dur_rbp_consensus_total_hits": int(sum(motif_counter.values())),
            "dur_rbp_consensus_unique_rbps": int(len(motif_counter)),
            "top_dur_rbp_consensus_hits": ";".join(
                f"{name}:{count}" for name, count in motif_counter.most_common(10)
            ),
        })
        counts = {pattern: overlapping_count(sequence, pattern) for pattern in are_patterns if len(pattern) <= 9}
        are_rows.append({
            "candidate_id": candidate_id,
            "are_total_hits": int(sum(counts.values())),
            "are_attTA_hits": int(overlapping_count(sequence, "ATTTA")),
            "are_tatttat_hits": int(overlapping_count(sequence, "TATTTAT")),
            "top_are_patterns": ";".join(
                f"{pattern}:{count}" for pattern, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:10] if count > 0
            ),
            "gre_total_hits": int(sum(overlapping_count(sequence, m) for m in GRE_MOTIFS)),
            "ure_total_hits": int(sum(overlapping_count(sequence, m) for m in URE_MOTIFS)),
        })

    mirna_hits_df = pd.concat(mirna_hit_tables, ignore_index=True) if mirna_hit_tables else pd.DataFrame(
        columns=["candidate_id", "family_id", "representative_mirna", "site_type", "site_count"]
    )
    mirna_summary_df = pd.DataFrame(mirna_summary_rows)
    dur_rbp_df = pd.DataFrame(dur_rbp_rows)
    are_df = pd.DataFrame(are_rows)
    write_table(mirna_hits_df, DATA_DIR / "mirna_dur_site_hits.tsv.gz")
    write_table(mirna_summary_df, DATA_DIR / "mirna_dur_site_summary.tsv.gz")
    write_table(dur_rbp_df, DATA_DIR / "dur_rbp_consensus_summary.tsv.gz")
    write_table(are_df, DATA_DIR / "dur_are_summary.tsv.gz")

    # --- RNAfold ---
    rnafold_records = []
    for row in utr_df.itertuples(index=False):
        if isinstance(row.short_utr_seq, str):
            rnafold_records.append((f"{row.candidate_id}__short", row.short_utr_seq))
        if isinstance(row.long_utr_seq, str):
            rnafold_records.append((f"{row.candidate_id}__long", row.long_utr_seq))
        if isinstance(row.dur_seq, str):
            rnafold_records.append((f"{row.candidate_id}__dur", row.dur_seq))
    rows: list[dict] = []
    run_records = [(record_id, seq) for record_id, seq in rnafold_records if seq and len(seq) <= RNAFOLD_MAX_LEN]
    skipped = [(record_id, seq) for record_id, seq in rnafold_records if not seq or len(seq) > RNAFOLD_MAX_LEN]
    if run_records:
        workers = max(1, RNAFOLD_WORKERS)
        with ProcessPoolExecutor(max_workers=workers) as executor:
            future_to_index = {
                executor.submit(run_one_rnafold, record): idx
                for idx, record in enumerate(run_records)
            }
            indexed_rows = []
            for future in as_completed(future_to_index):
                indexed_rows.append((future_to_index[future], future.result()))
        rows.extend(row for _, row in sorted(indexed_rows, key=lambda item: item[0]))
    for record_id, seq in skipped:
        rows.append(
            {
                "record_id": record_id,
                "seq_length": len(seq),
                "mfe_kcal_per_mol": math.nan,
                "mfe_per_100nt": math.nan,
                "rnafold_status": "too_long" if len(seq) > RNAFOLD_MAX_LEN else "empty_sequence",
            }
        )
    rnafold_df = pd.DataFrame(rows)
    write_table(rnafold_df, DATA_DIR / "rnafold_raw.tsv.gz")
    if not rnafold_df.empty:
        rnafold_df[["candidate_id", "sequence_label"]] = rnafold_df["record_id"].str.rsplit("__", n=1, expand=True)
        rnafold_wide = (
            rnafold_df.pivot(index="candidate_id", columns="sequence_label", values="mfe_per_100nt")
            .rename(columns={"dur": "dur_mfe_per_100nt", "long": "long_mfe_per_100nt", "short": "short_mfe_per_100nt"})
            .reset_index()
        )
        rnafold_wide["abs_long_short_mfe_delta_per_100nt"] = (
            rnafold_wide["long_mfe_per_100nt"] - rnafold_wide["short_mfe_per_100nt"]
        ).abs()
    else:
        rnafold_wide = pd.DataFrame(
            columns=[
                "candidate_id",
                "dur_mfe_per_100nt",
                "long_mfe_per_100nt",
                "short_mfe_per_100nt",
                "abs_long_short_mfe_delta_per_100nt",
            ]
        )
    write_table(rnafold_wide, DATA_DIR / "rnafold_candidate_summary.tsv.gz")
    # --- single-cell expression ---
    pseudobulk_df = pd.DataFrame(columns=["gene_id", "gene_id_stripped"])
    try:
        import anndata as ad  # type: ignore
    except ImportError:
        ad = None
    h5ad_path = None
    for fname in ["Velmeshev_2023_pre_postnatal.h5ad", "Fullard_2025_Aging_cohort.h5ad"]:
        if (SC_REF_DIR / fname).exists():
            h5ad_path = SC_REF_DIR / fname
            break
    adata = None
    if ad is not None and h5ad_path is not None:
        try:
            adata = ad.read_h5ad(h5ad_path)
        except Exception:  # noqa: BLE001
            adata = None
    cell_type_col = None
    if adata is not None:
        cell_type_col = next(
            (col for col in ("cell_type", "Cell_type", "celltype", "CellType", "cluster", "subtype", "Subtype", "broad_cell_type") if col in adata.obs.columns),
            None,
        )
    if cell_type_col is not None:
        cell_types = adata.obs[cell_type_col].astype(str).unique().tolist()
        X = adata.X
        gene_names = list(adata.var_names)

        mean_rows: dict[str, list[float]] = {ct: [] for ct in cell_types}
        for ct in cell_types:
            mask = (adata.obs[cell_type_col].astype(str) == ct).values
            sub = X[mask]
            if hasattr(sub, "toarray"):
                sub = sub.toarray()
            mean_rows[ct] = np.asarray(sub.mean(axis=0)).ravel().tolist()

        df = pd.DataFrame(mean_rows, index=gene_names)
        df.index.name = "gene_id"
        df = df.reset_index()
        df["gene_id_stripped"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        pseudobulk_df = df
    if pseudobulk_df.empty or len(pseudobulk_df.columns) <= 2:
        if BCA_PSEUDOBULK.exists():
            pseudobulk_df = pd.read_csv(BCA_PSEUDOBULK)
            gene_col = pseudobulk_df.columns[0]
            pseudobulk_df = pseudobulk_df.rename(columns={gene_col: "gene_id"})
            pseudobulk_df["gene_id_stripped"] = pseudobulk_df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        else:
            pseudobulk_df = pd.DataFrame(columns=["gene_id", "gene_id_stripped"])

    # --- single-cell expression summary ---
    candidate_cols = ["candidate_id", "gene_name", "stage", "tier"]
    merged = base_df[candidate_cols].drop_duplicates().merge(
        utr_df[["candidate_id", "gene_id"]],
        on="candidate_id",
        how="left",
    )
    if pseudobulk_df.empty or "gene_id_stripped" not in pseudobulk_df.columns:
        sc_expression_df = merged.assign(
            sc_expression_status="not_available",
            top_cell_types="",
            top_cell_type_expr=math.nan,
        )
    else:
        pseudobulk_lookup = pseudobulk_df.set_index("gene_id_stripped")
        value_cols = [col for col in pseudobulk_df.columns if col not in {"gene_id", "gene_id_stripped"}]
        if not value_cols:
            sc_expression_df = merged.assign(
                sc_expression_status="not_available",
                top_cell_types="",
                top_cell_type_expr=math.nan,
            )
        else:
            rows: list[dict] = []
            for row in merged.itertuples(index=False):
                gene_id = None if pd.isna(row.gene_id) else str(row.gene_id).split(".")[0]
                lookup_key: str | None = gene_id if (gene_id and gene_id in pseudobulk_lookup.index) else None
                if lookup_key is None:
                    # Fallback: match by gene symbol (handles h5ad files with gene-name var_names)
                    gene_name = str(row.gene_name) if not pd.isna(row.gene_name) else None
                    if gene_name and gene_name in pseudobulk_lookup.index:
                        lookup_key = gene_name
                if lookup_key is None:
                    rows.append(
                        {
                            "candidate_id": row.candidate_id,
                            "gene_name": row.gene_name,
                            "gene_id": row.gene_id,
                            "sc_expression_status": "missing_gene_in_pseudobulk",
                            "top_cell_types": "",
                            "top_cell_type_expr": math.nan,
                        }
                    )
                    continue
                values = pseudobulk_lookup.loc[lookup_key, value_cols]
                if isinstance(values, pd.DataFrame):
                    values = values.iloc[0]
                ranked = values.sort_values(ascending=False, kind="stable")
                rows.append(
                    {
                        "candidate_id": row.candidate_id,
                        "gene_name": row.gene_name,
                        "gene_id": row.gene_id,
                        "sc_expression_status": "ok",
                        "top_cell_types": ";".join(f"{cell}:{ranked[cell]:.4f}" for cell in ranked.index[:5]),
                        "top_cell_type_expr": float(ranked.iloc[0]),
                    }
                )
            sc_expression_df = pd.DataFrame(rows)

    write_table(sc_expression_df, DATA_DIR / "single_cell_expression_mapping.tsv.gz")

    # --- optional candidate evidence ---
    merged = base_df.copy()
    eqtl_path = DATA_DIR / "eqtl_coloc.tsv"
    use_eqtl_coloc = os.environ.get("APA05_USE_EQTL_COLOC", "0") == "1"
    eqtl = load_optional_table(eqtl_path) if use_eqtl_coloc else pd.DataFrame()
    if use_eqtl_coloc and not eqtl.empty and {"candidate_id", "PP.H4", "eQTL_coloc"}.issubset(eqtl.columns):
        eqtl["PP.H4"] = pd.to_numeric(eqtl["PP.H4"], errors="coerce")
        eqtl_summary = (
            eqtl.sort_values("PP.H4", ascending=False, kind="stable")
            .groupby("candidate_id", as_index=False)
            .agg(
                eQTL_coloc=("eQTL_coloc", "max"),
                eqtl_coloc_best_pp4=("PP.H4", "max"),
                eqtl_coloc_best_egene=("egene", "first"),
            )
        )
        merged = merged.merge(eqtl_summary, on="candidate_id", how="left")
        merged["eqtl_coloc_status"] = merged["eQTL_coloc"].map(
            lambda value: "available" if pd.notna(value) else "not_available"
        )
    else:
        merged["eQTL_coloc"] = pd.NA
        merged["eqtl_coloc_best_pp4"] = pd.NA
        merged["eqtl_coloc_best_egene"] = pd.NA
        merged["eqtl_coloc_status"] = "pending_extension"
    optional_specs = [
        (
            DATA_DIR / "protein_evidence.tsv",
            ["candidate_id"],
            ["protein_evidence", "ms_log2_fold_per_allele", "ms_pvalue", "ms_spearman_rho"],
        ),
        (
            DATA_DIR / "gwas_coloc.tsv",
            ["candidate_id"],
            ["GWAS_coloc", "gwas_coloc_pp4", "gwas_trait"],
        ),
        (
            DATA_DIR / "cross_stage_direction.tsv",
            ["candidate_id"],
            [
                "interaction_apa_qtl_status",
                "cross_stage_direction",
                "cross_stage_silenced_substatus",
                "slope_prenatal",
                "slope_postnatal",
            ],
        ),
    ]
    for path, key_cols, value_cols in optional_specs:
        evidence = load_optional_table(path)
        if evidence.empty or not set(key_cols).issubset(evidence.columns):
            continue
        available_cols = [col for col in value_cols if col in evidence.columns and col not in merged.columns]
        if not available_cols:
            continue
        evidence = evidence[key_cols + available_cols].drop_duplicates(subset=key_cols, keep="first")
        merged = merged.merge(evidence, on=key_cols, how="left")
    base_df = merged

    # --- mechanism classes ---
    merged = base_df.merge(pas_summary_df, on="candidate_id", how="left").merge(
        rbp_summary_df, on="candidate_id", how="left"
    )
    if not pas_auxiliary_summary_df.empty:
        merged = merged.merge(pas_auxiliary_summary_df, on="candidate_id", how="left")
    if "pas_auxiliary_motif_support" not in merged.columns:
        merged["pas_auxiliary_motif_support"] = False
        merged["top_pas_auxiliary_motif_class"] = ""
        merged["top_pas_auxiliary_motif_delta"] = 0
    for col, default in {
        "pas_direct_support_source": "not_supported",
        "pas_rescue_status": "not_supported",
    }.items():
        if col not in merged.columns:
            merged[col] = default
    merged["pas_direct_support"] = merged["pas_direct_support"].astype("boolean").fillna(False).astype(bool)
    merged["pas_auxiliary_motif_support"] = (
        merged["pas_auxiliary_motif_support"].astype("boolean").fillna(False).astype(bool)
    )
    merged["pas_direct_support"] = merged["pas_direct_support"] | merged["pas_auxiliary_motif_support"]
    merged["pas_direct_support_source"] = merged.apply(
        lambda row: (
            row["pas_direct_support_source"]
            if not row["pas_auxiliary_motif_support"]
            else "pas_auxiliary_motif_delta"
            if row["pas_direct_support_source"] == "not_supported"
            else f"{row['pas_direct_support_source']}+pas_auxiliary_motif_delta"
        ),
        axis=1,
    )
    merged["rbp_variant_delta_support"] = merged["rbp_variant_delta_support"].astype("boolean").fillna(False).astype(bool)
    merged["postar_rbp_support"] = (
        merged["postar_variant_overlap_count"].fillna(0).gt(0)
        | merged["postar_phenotype_overlap_count"].fillna(0).gt(0)
    )
    prox_distance = (merged["variant_pos"] - merged["proximal_pas_pos"]).abs()
    dist_distance = (merged["variant_pos"] - merged["distal_pas_pos"]).abs()
    merged["variant_near_resolved_pas_250"] = (
        prox_distance.le(250).fillna(False) | dist_distance.le(250).fillna(False)
    )
    merged["rbp_motif_near_pas_support"] = (
        merged["rbp_variant_delta_support"] & merged["variant_near_resolved_pas_250"]
    )
    merged["rbp_support"] = merged["postar_rbp_support"] | merged["rbp_motif_near_pas_support"]
    merged = annotate_effect_directionality(merged)

    merged["mechanism_class"] = merged.apply(
        lambda row: "mixed-support" if row["pas_direct_support"] and row["rbp_support"]
        else "PAS-direct" if row["pas_direct_support"]
        else "RBP-mediated" if row["rbp_support"]
        else "mechanistically unresolved",
        axis=1,
    )
    merged["mechanism_evidence_summary"] = merged.apply(
        lambda row: "; ".join(
            filter(
                None,
                [
                    f"PAS:{row['top_pas_support_label']}:{int(row['top_pas_motif_delta'])}" if pd.notna(row.get("top_pas_support_label")) else "",
                    f"PAS_source:{row['pas_direct_support_source']}" if pd.notna(row.get("pas_direct_support_source")) else "",
                    f"PAS_aux:{row['top_pas_auxiliary_motif_class']}:{row['top_pas_auxiliary_motif_delta']}" if bool(row.get("pas_auxiliary_motif_support", False)) else "",
                    f"RBP_gain:{row['top_rbp_gain']}:{row['top_rbp_gain_delta']}" if pd.notna(row.get("top_rbp_gain")) else "",
                    f"RBP_loss:{row['top_rbp_loss']}:{row['top_rbp_loss_delta']}" if pd.notna(row.get("top_rbp_loss")) else "",
                    f"POSTAR_var={int(row['postar_variant_overlap_count'])}",
                    f"POSTAR_pheno={int(row['postar_phenotype_overlap_count'])}",
                    f"Direction={row['utr_effect_direction']}:{row['effect_direction_status']}",
                ],
            )
        ),
        axis=1,
    )
    mechanism_main_df = merged[
        [
            "candidate_id",
            "gene_name",
            "stage",
            "tier",
            "pas_resolution_status",
            "pip",
            "pas_direct_support",
            "rbp_support",
            "variant_near_resolved_pas_250",
            "mechanism_class",
            "top_pas_support_label",
            "top_pas_motif_delta",
            "top_pas_variant_distance_nt",
            "ref_pas_canonical_hits",
            "alt_pas_canonical_hits",
            "pas_canonical_identity_changed",
            "pas_direct_support_source",
            "pas_rescue_status",
            "pas_auxiliary_motif_support",
            "top_pas_auxiliary_motif_class",
            "top_pas_auxiliary_motif_delta",
            "top_rbp_gain",
            "top_rbp_gain_delta",
            "top_rbp_loss",
            "top_rbp_loss_delta",
            "postar_variant_overlap_count",
            "postar_phenotype_overlap_count",
            "utr_effect_direction",
            "effect_direction_status",
            "effect_direction_source_column",
            "mechanism_evidence_summary",
        ]
    ]

    # --- functional effects ---
    merged = (
        fate_df.merge(utr_df.drop(columns=["short_utr_seq", "long_utr_seq", "dur_seq"]), on="candidate_id", how="left")
        .merge(mirna_summary_df, on="candidate_id", how="left")
        .merge(dur_rbp_df, on="candidate_id", how="left")
        .merge(are_df, on="candidate_id", how="left")
        .merge(rnafold_wide, on="candidate_id", how="left")
    )
    fill_zero_cols = [
        "mirna_total_sites",
        "mirna_unique_families",
        "mirna_8mer_sites",
        "mirna_7mer_m8_sites",
        "mirna_7mer_a1_sites",
        "dur_rbp_consensus_total_hits",
        "dur_rbp_consensus_unique_rbps",
        "are_total_hits",
        "are_attTA_hits",
        "are_tatttat_hits",
        "gre_total_hits",
        "ure_total_hits",
    ]
    for col in fill_zero_cols:
        if col in merged.columns:
            merged[col] = merged[col].fillna(0)
    dur_length_nonzero = merged["dur_length"].where(merged["dur_length"].fillna(0).gt(0))
    merged["mirna_sites_per_100nt"] = (
        merged["mirna_total_sites"] / dur_length_nonzero * 100.0
    ).fillna(0.0)
    merged["dur_rbp_hits_per_100nt"] = (
        merged["dur_rbp_consensus_total_hits"] / dur_length_nonzero * 100.0
    ).fillna(0.0)
    merged["are_hits_per_100nt"] = (
        merged["are_total_hits"] / dur_length_nonzero * 100.0
    ).fillna(0.0)
    merged["gre_hits_per_100nt"] = (
        merged["gre_total_hits"] / dur_length_nonzero * 100.0
    ).fillna(0.0)
    merged["ure_hits_per_100nt"] = (
        merged["ure_total_hits"] / dur_length_nonzero * 100.0
    ).fillna(0.0)
    density_thresholds = {
        "mirna": nonzero_quantile(merged["mirna_sites_per_100nt"], DENSITY_SUPPORT_QUANTILE),
        "rbp": nonzero_quantile(merged["dur_rbp_hits_per_100nt"], DENSITY_SUPPORT_QUANTILE),
        "are": nonzero_quantile(merged["are_hits_per_100nt"], DENSITY_SUPPORT_QUANTILE),
        "gre": nonzero_quantile(merged["gre_hits_per_100nt"], DENSITY_SUPPORT_QUANTILE),
        "ure": nonzero_quantile(merged["ure_hits_per_100nt"], DENSITY_SUPPORT_QUANTILE),
    }
    merged["mirna_density_threshold_per_100nt"] = density_thresholds["mirna"]
    merged["dur_rbp_density_threshold_per_100nt"] = density_thresholds["rbp"]
    merged["are_density_threshold_per_100nt"] = density_thresholds["are"]
    merged["gre_density_threshold_per_100nt"] = density_thresholds["gre"]
    merged["ure_density_threshold_per_100nt"] = density_thresholds["ure"]
    merged["miRNA_support"] = merged["mirna_total_sites"].gt(0)
    merged["RBP_support"] = merged["dur_rbp_consensus_total_hits"].gt(0) | merged["postar_dur_overlap_count"].fillna(0).gt(0)
    merged["ARE_support"] = merged["are_total_hits"].gt(0)
    merged["GRE_support"] = merged["gre_total_hits"].gt(0)
    merged["URE_support"] = merged["ure_total_hits"].gt(0)
    merged["RNA_structure_support"] = merged["abs_long_short_mfe_delta_per_100nt"].fillna(0).ge(1.0)
    merged["miRNA_density_supported"] = (
        merged["miRNA_support"]
        & merged["mirna_sites_per_100nt"].ge(density_thresholds["mirna"])
    )
    merged["RBP_density_supported"] = (
        merged["postar_dur_overlap_count"].fillna(0).gt(0)
        | (
            merged["dur_rbp_consensus_total_hits"].gt(0)
            & merged["dur_rbp_hits_per_100nt"].ge(density_thresholds["rbp"])
        )
    )
    merged["ARE_density_supported"] = (
        merged["ARE_support"]
        & merged["are_hits_per_100nt"].ge(density_thresholds["are"])
    )
    merged["GRE_density_supported"] = (
        merged["GRE_support"]
        & merged["gre_hits_per_100nt"].ge(density_thresholds["gre"])
    )
    merged["URE_density_supported"] = (
        merged["URE_support"]
        & merged["ure_hits_per_100nt"].ge(density_thresholds["ure"])
    )
    merged = annotate_effect_directionality(merged)
    merged["raw_functional_module_count"] = (
        merged["miRNA_support"].astype(int)
        + merged["RBP_support"].astype(int)
        + merged["ARE_support"].astype(int)
        + merged["GRE_support"].astype(int)
        + merged["URE_support"].astype(int)
        + merged["RNA_structure_support"].astype(int)
    )
    merged["functional_module_count"] = (
        merged["miRNA_density_supported"].astype(int)
        + merged["RBP_density_supported"].astype(int)
        + merged["ARE_density_supported"].astype(int)
        + merged["GRE_density_supported"].astype(int)
        + merged["URE_density_supported"].astype(int)
        + merged["RNA_structure_support"].astype(int)
    )
    merged["density_bias_flag"] = merged["raw_functional_module_count"] > merged["functional_module_count"]
    merged["density_supported_modules"] = merged.apply(
        lambda row: summarize_supported_modules(
            row,
            [
                ("miRNA", "miRNA_density_supported"),
                ("RBP", "RBP_density_supported"),
                ("ARE", "ARE_density_supported"),
                ("GRE", "GRE_density_supported"),
                ("URE", "URE_density_supported"),
                ("RNA_structure", "RNA_structure_support"),
            ],
        ),
        axis=1,
    )
    merged["raw_only_load_modules"] = merged.apply(
        lambda row: summarize_supported_modules(
            row,
            [
                ("miRNA", "miRNA_support"),
                ("RBP", "RBP_support"),
                ("ARE", "ARE_support"),
                ("GRE", "GRE_support"),
                ("URE", "URE_support"),
            ],
        ),
        axis=1,
    )

    merged["functional_support_category"] = merged.apply(
        lambda row: "density-supported multi-module hits" if row["functional_module_count"] >= 2
        else "density-supported single-module hit" if row["functional_module_count"] == 1
        else "raw-count-only heuristic load" if row["raw_functional_module_count"] >= 1
        else "limited heuristic support",
        axis=1,
    )
    functional_main_df = merged[
        [
            "candidate_id",
            "gene_name",
            "stage",
            "tier",
            "pas_resolution_status",
            "dur_length",
            "utr_status",
            "mirna_total_sites",
            "mirna_sites_per_100nt",
            "mirna_density_threshold_per_100nt",
            "mirna_unique_families",
            "top_mirna_families",
            "dur_rbp_consensus_total_hits",
            "dur_rbp_hits_per_100nt",
            "dur_rbp_density_threshold_per_100nt",
            "dur_rbp_consensus_unique_rbps",
            "top_dur_rbp_consensus_hits",
            "are_total_hits",
            "are_hits_per_100nt",
            "are_density_threshold_per_100nt",
            "are_attTA_hits",
            "are_tatttat_hits",
            "gre_total_hits",
            "gre_hits_per_100nt",
            "gre_density_threshold_per_100nt",
            "ure_total_hits",
            "ure_hits_per_100nt",
            "ure_density_threshold_per_100nt",
            "long_mfe_per_100nt",
            "short_mfe_per_100nt",
            "abs_long_short_mfe_delta_per_100nt",
            "postar_dur_overlap_count",
            "utr_effect_direction",
            "effect_direction_status",
            "effect_direction_source_column",
            "raw_functional_module_count",
            "functional_module_count",
            "density_bias_flag",
            "miRNA_density_supported",
            "RBP_density_supported",
            "ARE_density_supported",
            "GRE_density_supported",
            "URE_density_supported",
            "RNA_structure_support",
            "density_supported_modules",
            "raw_only_load_modules",
            "functional_support_category",
        ]
    ]
    # --- APARENT ---
    model_path = APARENT_DIR / "saved_models" / APARENT_MODEL_NAME
    if not model_path.exists():
        n = len(mechanism_df)
        aparent_scores_df = pd.DataFrame({
            "candidate_id": mechanism_df["candidate_id"].tolist(),
            "APARENT_delta": [float("nan")] * n,
            "aparent_scoring_status": ["model_not_found"] * n,
        })
    else:
        try:
            if str(APARENT_DIR) not in sys.path:
                sys.path.insert(0, str(APARENT_DIR))
            import keras  # type: ignore
            from aparent.predictor import get_aparent_encoder  # type: ignore
            model = keras.models.load_model(str(model_path))
            encoder = get_aparent_encoder()
        except Exception as exc:  # noqa: BLE001
            n = len(mechanism_df)
            aparent_scores_df = pd.DataFrame({
                "candidate_id": mechanism_df["candidate_id"].tolist(),
                "APARENT_delta": [float("nan")] * n,
                "aparent_scoring_status": [f"import_error:{type(exc).__name__}"] * n,
            })
        else:
            rows: list[dict] = []
            for row in mechanism_df.itertuples(index=False):
                cid = row.candidate_id
                chrom = str(row.variant_chr) if hasattr(row, "variant_chr") else str(row.phenotype_chr)
                strand = str(row.phenotype_strand)
                variant_pos = int(row.variant_pos)
                ref_allele = str(row.ref_allele)
                alt_allele = str(row.alt_allele)

                try:
                    pas_pos = int(row.proximal_pas_pos)
                except (AttributeError, ValueError, TypeError):
                    rows.append({"candidate_id": cid, "APARENT_delta": float("nan"), "aparent_scoring_status": "no_pas_pos"})
                    continue

                win_start = pas_pos - APARENT_WINDOW_HALF
                win_end = pas_pos + APARENT_WINDOW_HALF
                if win_start < 1:
                    rows.append({"candidate_id": cid, "APARENT_delta": float("nan"), "aparent_scoring_status": "window_out_of_bounds"})
                    continue

                ref_seq = fetch_oriented_seq(fasta, chrom, win_start, win_end, strand)
                if len(ref_seq) != 205:
                    rows.append({"candidate_id": cid, "APARENT_delta": float("nan"), "aparent_scoring_status": f"seq_len_error:{len(ref_seq)}"})
                    continue

                alt_seq, apply_status = apply_variant_to_oriented_seq(
                    ref_seq, win_start, win_end, strand, variant_pos, ref_allele, alt_allele
                )
                if apply_status != "applied":
                    rows.append({"candidate_id": cid, "APARENT_delta": float("nan"), "aparent_scoring_status": f"variant_{apply_status}"})
                    continue

                try:
                    iso_ref, _ = model.predict(encoder([ref_seq]), verbose=0)
                    iso_alt, _ = model.predict(encoder([alt_seq]), verbose=0)
                    p_ref = float(np.clip(iso_ref[0, 0], 1e-7, 1 - 1e-7))
                    p_alt = float(np.clip(iso_alt[0, 0], 1e-7, 1 - 1e-7))
                    delta = math.log(p_alt / (1.0 - p_alt)) - math.log(p_ref / (1.0 - p_ref))
                    rows.append({"candidate_id": cid, "APARENT_delta": round(delta, 4), "aparent_scoring_status": "ok"})
                except Exception as exc:  # noqa: BLE001
                    rows.append({"candidate_id": cid, "APARENT_delta": float("nan"), "aparent_scoring_status": f"predict_error:{type(exc).__name__}"})
            aparent_scores_df = pd.DataFrame(rows)
    write_table(aparent_scores_df, DATA_DIR / "aparent_variant_scores.tsv")
    if "APARENT_delta" not in base_df.columns:
        base_df = base_df.merge(
            aparent_scores_df[["candidate_id", "APARENT_delta", "aparent_scoring_status"]],
            on="candidate_id",
            how="left",
        )
    else:
        computed = aparent_scores_df.set_index("candidate_id")["APARENT_delta"]
        status = aparent_scores_df.set_index("candidate_id")["aparent_scoring_status"]
        base_df["APARENT_delta"] = base_df.apply(
            lambda r: computed.get(r["candidate_id"], r["APARENT_delta"])
            if pd.isna(r["APARENT_delta"])
            else r["APARENT_delta"],
            axis=1,
        )
        if "aparent_scoring_status" not in base_df.columns:
            base_df["aparent_scoring_status"] = base_df["candidate_id"].map(status)

    write_table(mechanism_main_df, DATA_DIR / "mechanism_main_table.tsv.gz")
    write_table(functional_main_df, DATA_DIR / "functional_effect_main_table.tsv.gz")
    # --- framework outputs ---
    base_lookup = base_df.drop_duplicates(subset=["candidate_id"], keep="first").set_index("candidate_id")
    merged = (
        mechanism_main_df.merge(functional_main_df, on="candidate_id", how="left", suffixes=("", "_functional"))
        .merge(
            base_df[
                [
                    "candidate_id",
                    "phenotype_id",
                    "phenotype_chr",
                    "phenotype_start",
                    "phenotype_end",
                    "phenotype_strand",
                    "variant_id",
                    "variant_chr",
                    "variant_pos",
                    "ref_allele",
                    "alt_allele",
                    "variant_phastcons",
                    "phenotype_mean_phastcons",
                    "phenotype_max_phastcons",
                    "proximal_pas_phastcons",
                    "distal_pas_phastcons",
                    "dur_mean_phastcons",
                    "dur_max_phastcons",
                ]
            ],
            on="candidate_id",
            how="left",
        )
        .merge(
            sc_expression_df[["candidate_id", "top_cell_types", "top_cell_type_expr"]],
            on="candidate_id",
            how="left",
        )
    )

    merged["framework_tier"] = merged["pip"].map(
        lambda pip: "interaction_context" if pd.isna(pip)
        else "tier1_core" if float(pip) >= FRAMEWORK_TIER1_PIP_THRESHOLD
        else "tier2_extended" if float(pip) >= FRAMEWORK_TIER2_PIP_THRESHOLD
        else "below_framework_threshold"
    )
    merged["framework_tier_rank"] = merged["framework_tier"].map(
        lambda label: {"tier1_core": 0, "tier2_extended": 1, "interaction_context": 2, "below_framework_threshold": 3}.get(str(label), 99)
    )
    merged["high_confidence_subset"] = merged["pip"].fillna(0).ge(FRAMEWORK_HIGH_CONFIDENCE_PIP_THRESHOLD)
    merged["pas_disrupted"] = merged.apply(
        lambda row: "not_supported" if not bool(row.get("pas_direct_support", False))
        else "Gain" if not canonical_pas_hits(row.get("ref_pas_canonical_hits", "")) and canonical_pas_hits(row.get("alt_pas_canonical_hits", ""))
        else "Loss" if canonical_pas_hits(row.get("ref_pas_canonical_hits", "")) and not canonical_pas_hits(row.get("alt_pas_canonical_hits", ""))
        else "Modulation",
        axis=1,
    )
    aparent_col = resolve_optional_column(base_df, ("APARENT_delta", "aparent_delta", "apparent_delta"))
    if aparent_col is None:
        merged["apparent_delta"] = pd.NA
        merged["aparent_status"] = "not_available_in_input"
    else:
        merged["apparent_delta"] = pd.to_numeric(
            base_lookup.reindex(merged["candidate_id"])[aparent_col].values, errors="coerce"
        )
        merged["aparent_status"] = merged["apparent_delta"].map(
            lambda value: "available_from_input" if pd.notna(value) else "missing_value_in_input"
        )

    conservation_columns = [
        "variant_phastcons",
        "phenotype_mean_phastcons",
        "phenotype_max_phastcons",
        "proximal_pas_phastcons",
        "distal_pas_phastcons",
        "dur_mean_phastcons",
        "dur_max_phastcons",
    ]
    merged["conservation_score"] = merged[conservation_columns].max(axis=1)
    conservation_threshold = nonzero_quantile(merged["conservation_score"], DENSITY_SUPPORT_QUANTILE)
    merged["conservation_support_threshold"] = conservation_threshold
    merged["conservation_supported"] = (
        merged["conservation_score"].fillna(0).gt(0)
        & merged["conservation_score"].fillna(0).ge(conservation_threshold)
    )

    protein_col = resolve_optional_column(base_df, ("protein_evidence",))
    if protein_col is None:
        merged["protein_evidence"] = "not_evaluated"
    else:
        merged["protein_evidence"] = base_lookup.reindex(merged["candidate_id"])[protein_col].fillna("not_evaluated").values

    eqtl_col = resolve_optional_column(base_df, ("eQTL_coloc", "eqtl_coloc"))
    if eqtl_col is None:
        merged["eqtl_coloc"] = pd.NA
        merged["eqtl_coloc_status"] = "pending_extension"
    else:
        merged["eqtl_coloc"] = base_lookup.reindex(merged["candidate_id"])[eqtl_col].astype("boolean").values
        status_col = resolve_optional_column(base_df, ("eqtl_coloc_status", "eQTL_coloc_status"))
        if status_col is None:
            merged["eqtl_coloc_status"] = merged["eqtl_coloc"].map(
                lambda value: "available" if pd.notna(value) else "not_available"
            )
        else:
            merged["eqtl_coloc_status"] = (
                base_lookup.reindex(merged["candidate_id"])[status_col].fillna("not_available").values
            )

    gwas_col = resolve_optional_column(base_df, ("GWAS_coloc", "gwas_coloc"))
    if gwas_col is None:
        merged["gwas_coloc"] = pd.NA
    else:
        merged["gwas_coloc"] = base_lookup.reindex(merged["candidate_id"])[gwas_col].astype("boolean").values
    for col, default in {
        "interaction_apa_qtl_status": "not_tested",
        "cross_stage_direction": "",
        "cross_stage_silenced_substatus": "",
    }.items():
        if col not in merged.columns:
            merged[col] = default

    merged["dur_regulatory_effect"] = merged["utr_effect_direction"].map(
        lambda direction: {"lengthening": "gained_in_expressed_region", "shortening": "lost_from_expressed_region"}.get(normalize_effect_direction(direction), "unknown")
    )
    merged["tier1_mechanism_priority"] = merged["framework_tier"].eq("tier1_core")
    merged["tier2_enrichment_only"] = merged["framework_tier"].eq("tier2_extended")
    merged["module_a_eligible"] = merged["framework_tier"].isin(["tier1_core", "tier2_extended", "interaction_context"])
    merged["module_b_eligible"] = merged["dur_length"].fillna(0).gt(0)
    merged["module_c_eligible"] = merged["variant_phastcons"].notna() | merged["dur_mean_phastcons"].notna()
    merged["module_d_status"] = merged["protein_evidence"].map(
        lambda value: "available_from_input" if value not in {"", "not_evaluated", "MS_not_detected"} else "not_evaluated"
    )
    merged["module_e_export_ready"] = True

    merged["pip_norm"] = minmax_scale(merged["pip"])
    merged["aparent_delta_norm"] = minmax_scale(merged["apparent_delta"].abs())
    merged["mirna_density_norm"] = minmax_scale(merged["mirna_sites_per_100nt"].fillna(0.0))
    merged["rbp_density_norm"] = minmax_scale(merged["dur_rbp_hits_per_100nt"].fillna(0.0))
    merged["phastcons_norm"] = minmax_scale(merged["conservation_score"])
    merged["TotalScore"], merged["TotalScore_status"] = compute_reweighted_total_score(
        merged,
        TOTAL_SCORE_WEIGHTS,
    )

    merged["PAS_disrupted"] = merged["pas_disrupted"].where(
        merged["pas_disrupted"].ne("not_supported"),
        "not_supported",
    )
    if "pas_auxiliary_motif_support" in merged.columns:
        merged["PAS_disrupted"] = merged["PAS_disrupted"].where(
            ~(
                merged["PAS_disrupted"].eq("not_supported")
                & merged["pas_auxiliary_motif_support"].fillna(False).astype(bool)
            ),
            "Auxiliary",
        )
    merged["APARENT_delta"] = merged["apparent_delta"]
    merged["DUR_miRNA_affected"] = merged["miRNA_density_supported"].fillna(False).astype(bool)
    merged["DUR_RBP_affected"] = merged["RBP_density_supported"].fillna(False).astype(bool)
    chain_components = pd.DataFrame(
        {
            "PAS_disrupted": merged["PAS_disrupted"].ne("not_supported"),
            "APARENT_delta": merged["APARENT_delta"].notna(),
            "DUR_miRNA_affected": merged["DUR_miRNA_affected"].fillna(False).astype(bool),
            "DUR_RBP_affected": merged["DUR_RBP_affected"].fillna(False).astype(bool),
            "conservation_supported": merged["conservation_supported"].fillna(False).astype(bool),
            "eQTL_coloc": merged["eqtl_coloc"].fillna(False).astype(bool),
            "GWAS_coloc": merged["gwas_coloc"].fillna(False).astype(bool),
            "protein_evidence": merged["protein_evidence"].fillna("").map(
                lambda value: value not in {"", "not_evaluated", "MS_not_detected"}
            ),
        }
    )
    merged["chain_supported_count"] = chain_components.sum(axis=1)
    merged["chain_completeness"] = merged["chain_supported_count"].astype(int).astype(str) + "/8"

    candidates_prioritized = (
        merged[
            [
                "candidate_id",
                "phenotype_id",
                "variant_id",
                "gene_name",
                "stage",
                "tier",
                "framework_tier",
                "framework_tier_rank",
                "high_confidence_subset",
                "pip",
                "utr_effect_direction",
                "pas_disrupted",
                "mechanism_class",
                "conservation_supported",
                "conservation_score",
                "conservation_support_threshold",
                "module_a_eligible",
                "module_b_eligible",
                "module_c_eligible",
                "module_d_status",
                "module_e_export_ready",
                "interaction_apa_qtl_status",
                "cross_stage_direction",
                "cross_stage_silenced_substatus",
                "top_cell_types",
                "TotalScore",
                "TotalScore_status",
            ]
        ]
        .sort_values(["framework_tier_rank", "high_confidence_subset", "pip", "TotalScore"], ascending=[True, False, False, False], kind="stable")
        .drop(columns=["framework_tier_rank"])
    )

    pas_disruption_scores = merged[
        [
            "candidate_id",
            "phenotype_id",
            "variant_id",
            "gene_name",
            "stage",
            "tier",
            "pip",
            "PAS_disrupted",
            "top_pas_support_label",
            "top_pas_motif_delta",
            "top_pas_variant_distance_nt",
            "ref_pas_canonical_hits",
            "alt_pas_canonical_hits",
            "pas_direct_support_source",
            "pas_rescue_status",
            "pas_auxiliary_motif_support",
            "top_pas_auxiliary_motif_class",
            "top_pas_auxiliary_motif_delta",
            "APARENT_delta",
            "aparent_status",
        ]
    ].rename(columns={"PAS_disrupted": "pas_disrupted"})

    dur_regulatory_load = merged[
        [
            "candidate_id",
            "gene_name",
            "stage",
            "tier",
            "phenotype_strand",
            "utr_effect_direction",
            "dur_regulatory_effect",
            "dur_length",
            "mirna_total_sites",
            "mirna_sites_per_100nt",
            "miRNA_density_supported",
            "dur_rbp_consensus_total_hits",
            "dur_rbp_hits_per_100nt",
            "RBP_density_supported",
            "are_total_hits",
            "are_hits_per_100nt",
            "ARE_density_supported",
            "conservation_score",
            "long_mfe_per_100nt",
            "short_mfe_per_100nt",
            "abs_long_short_mfe_delta_per_100nt",
            "functional_support_category",
            "density_supported_modules",
            "raw_only_load_modules",
        ]
    ]

    functional_support = merged[
        [
            "candidate_id",
            "gene_name",
            "stage",
            "tier",
            "framework_tier",
            "pip",
            "high_confidence_subset",
            "mechanism_class",
            "pas_disrupted",
            "utr_effect_direction",
            "dur_regulatory_effect",
            "functional_support_category",
            "density_supported_modules",
            "variant_phastcons",
            "dur_mean_phastcons",
            "phenotype_mean_phastcons",
            "conservation_score",
            "conservation_support_threshold",
            "conservation_supported",
            "mirna_sites_per_100nt",
            "dur_rbp_hits_per_100nt",
            "are_hits_per_100nt",
            "APARENT_delta",
            "TotalScore",
            "TotalScore_status",
            "top_cell_types",
        ]
    ].sort_values(["TotalScore", "pip", "gene_name"], ascending=[False, False, True], kind="stable")
    functional_support["dur_are_density_supported"] = merged.loc[
        functional_support.index,
        "ARE_density_supported",
    ].values

    evidence_chain_summary = merged[
        [
            "candidate_id",
            "PAS_disrupted",
            "APARENT_delta",
            "DUR_miRNA_affected",
            "DUR_RBP_affected",
            "conservation_supported",
            "eqtl_coloc",
            "eqtl_coloc_status",
            "gwas_coloc",
            "protein_evidence",
            "chain_completeness",
            "aparent_status",
            "module_d_status",
            "cross_stage_direction",
            "mechanism_class",
            "functional_support_category",
            "TotalScore",
        ]
    ].rename(columns={"eqtl_coloc": "eQTL_coloc", "gwas_coloc": "GWAS_coloc"})
    # --- total score sensitivity ---
    components = list(TOTAL_SCORE_WEIGHTS)
    base_top = set(
        merged.sort_values("TotalScore", ascending=False, kind="stable")
        .head(50)["candidate_id"]
        .astype(str)
    )
    base_rank = merged.set_index("candidate_id")["TotalScore"].rank(ascending=False, method="average")
    rows: list[dict] = []
    for component in components:
        for component_weight in (0.1, 0.2, 0.3):
            other_weight = (1.0 - component_weight) / (len(components) - 1)
            weights = {col: other_weight for col in components}
            weights[component] = component_weight
            score, _ = compute_reweighted_total_score(merged, weights)
            ranked = merged[["candidate_id"]].copy()
            ranked["sensitivity_score"] = score
            top = set(
                ranked.sort_values("sensitivity_score", ascending=False, kind="stable")
                .head(50)["candidate_id"]
                .astype(str)
            )
            rank = ranked.set_index("candidate_id")["sensitivity_score"].rank(ascending=False, method="average")
            shared = base_rank.index.intersection(rank.index)
            rows.append(
                {
                    "perturbed_component": component,
                    "perturbed_weight": component_weight,
                    "top50_jaccard": len(base_top & top) / len(base_top | top) if base_top or top else pd.NA,
                    "rank_spearman_rho": base_rank.loc[shared].corr(rank.loc[shared], method="spearman")
                    if len(shared) >= 2
                    else pd.NA,
                    "n_candidates": int(merged.shape[0]),
                }
            )
    total_score_sensitivity = pd.DataFrame(rows)

    framework_outputs = {
        "candidates_prioritized": candidates_prioritized,
        "pas_disruption_scores": pas_disruption_scores,
        "dur_regulatory_load": dur_regulatory_load,
        "functional_support": functional_support,
        "evidence_chain_summary": evidence_chain_summary,
        "total_score_sensitivity": total_score_sensitivity,
    }
    # --- APARENT internal concordance ---
    qtl_beta_col = next((col for col in SIGNED_EFFECT_COLUMNS if col in base_df.columns), None)
    aparent_col = next((col for col in ("APARENT_delta", "aparent_delta", "apparent_delta") if col in base_df.columns), None)
    if qtl_beta_col is None or aparent_col is None:
        aparent_internal_concordance_df = pd.DataFrame(
            [
                {
                    "analysis_status": "not_run",
                    "reason": "missing_aparent_or_qtl_beta_columns",
                    "aparent_delta_column": aparent_col or "",
                    "qtl_beta_column": qtl_beta_col or "",
                    "n_tested_candidates": 0,
                    "spearman_rho": pd.NA,
                    "spearman_pvalue": pd.NA,
                }
            ]
        )
    else:
        concordance_df = base_df[[aparent_col, qtl_beta_col]].dropna()
        if concordance_df.empty:
            aparent_internal_concordance_df = pd.DataFrame(
                [
                    {
                        "analysis_status": "not_run",
                        "reason": "no_complete_aparent_qtl_beta_pairs",
                        "aparent_delta_column": aparent_col,
                        "qtl_beta_column": qtl_beta_col,
                        "n_tested_candidates": 0,
                        "spearman_rho": pd.NA,
                        "spearman_pvalue": pd.NA,
                    }
                ]
            )
        else:
            try:
                from scipy.stats import spearmanr  # type: ignore

                rho, pval = spearmanr(
                    pd.to_numeric(concordance_df[aparent_col], errors="coerce"),
                    pd.to_numeric(concordance_df[qtl_beta_col], errors="coerce"),
                    nan_policy="omit",
                )
            except Exception:  # noqa: BLE001
                rho = concordance_df[aparent_col].corr(concordance_df[qtl_beta_col], method="spearman")
                pval = pd.NA
            aparent_internal_concordance_df = pd.DataFrame(
                [
                    {
                        "analysis_status": "ok",
                        "reason": "",
                        "aparent_delta_column": aparent_col,
                        "qtl_beta_column": qtl_beta_col,
                        "n_tested_candidates": int(concordance_df.shape[0]),
                        "spearman_rho": rho,
                        "spearman_pvalue": pval,
                    }
                ]
            )

    write_table(framework_outputs["candidates_prioritized"], DATA_DIR / "candidates_prioritized.tsv")
    write_table(framework_outputs["pas_disruption_scores"], DATA_DIR / "pas_disruption_scores.tsv")
    write_table(aparent_internal_concordance_df, DATA_DIR / "aparent_internal_concordance.tsv")
    write_table(framework_outputs["dur_regulatory_load"], DATA_DIR / "dur_regulatory_load.tsv")
    write_table(framework_outputs["functional_support"], DATA_DIR / "functional_support.tsv.gz")
    write_table(framework_outputs["evidence_chain_summary"], DATA_DIR / "evidence_chain_summary.tsv")
    write_table(framework_outputs["total_score_sensitivity"], DATA_DIR / "total_score_sensitivity.tsv")
    stratified_mechanism = (
        mechanism_main_df.groupby(["stage", "tier", "mechanism_class"], as_index=False)
        .size()
        .rename(columns={"size": "n_candidates"})
    )
    stratified_mechanism["summary_type"] = "mechanism"
    stratified_function = (
        functional_main_df.groupby(["stage", "tier", "functional_support_category"], as_index=False)
        .size()
        .rename(columns={"functional_support_category": "label", "size": "n_candidates"})
    )
    stratified_function["summary_type"] = "functional_effect"
    stratified_function = stratified_function.rename(columns={"label": "mechanism_class"})
    stratified_summary = pd.concat([stratified_mechanism, stratified_function], ignore_index=True, sort=False)
    write_table(stratified_summary, DATA_DIR / "stratified_summary.tsv")
    priority_df = (
        mechanism_main_df.merge(
            functional_main_df[
                [
                    "candidate_id",
                    "raw_functional_module_count",
                    "functional_module_count",
                    "functional_support_category",
                    "density_bias_flag",
                    "density_supported_modules",
                    "raw_only_load_modules",
                    "mirna_total_sites",
                    "mirna_sites_per_100nt",
                    "dur_rbp_consensus_total_hits",
                    "dur_rbp_hits_per_100nt",
                    "are_total_hits",
                    "are_hits_per_100nt",
                ]
            ],
            on="candidate_id",
            how="left",
        )
        .merge(
            sc_expression_df[["candidate_id", "top_cell_types", "top_cell_type_expr"]],
            on="candidate_id",
            how="left",
        )
    )
    priority_df["functional_module_count"] = priority_df["functional_module_count"].fillna(0)
    priority_df["priority_score"] = (
        priority_df["pip"].fillna(0) * 10
        + priority_df["pas_direct_support"].astype(int) * 2
        + priority_df["rbp_support"].astype(int) * 2
        + priority_df["functional_module_count"]
    )
    priority_df = priority_df.sort_values(
        ["priority_score", "pip", "functional_module_count", "gene_name"],
        ascending=[False, False, False, True],
        kind="stable",
    ).head(PRIORITY_CASE_LIMIT)
    write_table(priority_df, DATA_DIR / "priority_cases.tsv")
    run_summary = {
        "input_table": str(INPUT_TABLE),
        "output_dir": str(DATA_DIR),
        "genome_fasta": str(GENOME_FASTA),
        "gtf_path": str(GTF_PATH),
        "bca_pseudobulk": str(BCA_PSEUDOBULK),
        "argv": sys.argv[1:],
        "n_candidates_mechanism": int(mechanism_df.shape[0]),
        "n_candidates_functional_effect": int(fate_df.shape[0]),
        "n_utr_success": int((utr_df["utr_status"] == "ok").sum()),
        "n_pas_direct_supported": int(mechanism_main_df["pas_direct_support"].sum()),
        "n_rbp_supported": int(mechanism_main_df["rbp_support"].sum()),
        "mechanism_class_counts": mechanism_main_df["mechanism_class"].value_counts().to_dict(),
        "functional_support_counts": functional_main_df["functional_support_category"].value_counts().to_dict(),
        "n_priority_cases": int(priority_df.shape[0]),
        "n_framework_candidates_prioritized": int(framework_outputs["candidates_prioritized"].shape[0]),
        "density_support_quantile": DENSITY_SUPPORT_QUANTILE,
        "density_thresholds_per_100nt": {
            "mirna_sites_per_100nt": float(functional_main_df["mirna_density_threshold_per_100nt"].iloc[0]) if not functional_main_df.empty else 0.0,
            "dur_rbp_hits_per_100nt": float(functional_main_df["dur_rbp_density_threshold_per_100nt"].iloc[0]) if not functional_main_df.empty else 0.0,
            "are_hits_per_100nt": float(functional_main_df["are_density_threshold_per_100nt"].iloc[0]) if not functional_main_df.empty else 0.0,
        },
        "effect_direction_status_counts": mechanism_main_df["effect_direction_status"].value_counts().to_dict(),
        "pas_rescue_status_counts": mechanism_main_df["pas_rescue_status"].value_counts().to_dict(),
        "optional_extensions": {
            "eqtl_coloc": "pending_extension_enable_with_APA05_RUN_EQTL_COLOC_or_APA05_USE_EQTL_COLOC",
            "protein_sample_metadata": "not_required_direct_MS_columns_only",
            "pfc_rbp_expression": "not_available_level_2_RBP_expression_support_disabled",
            "enrichment_gene_sets": "external_GMT_optional_use_package_resources_when_available",
            "tested_apa_gene_universe": "data/07_coloc/tested_apa_gene_universe.tsv",
        },
        "framework_outputs": [
            "candidates_prioritized.tsv",
            "pas_disruption_scores.tsv",
            "aparent_internal_concordance.tsv",
            "dur_regulatory_load.tsv",
            "functional_support.tsv.gz",
            "evidence_chain_summary.tsv",
        ],
    }
    with (DATA_DIR / "run_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(run_summary, handle, indent=2, ensure_ascii=False)
    fasta.close()
