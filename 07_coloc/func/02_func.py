#!/usr/bin/env python3
# Functional annotation of APA-QTL candidates: PAS sites, phastCons conservation, POSTAR3 CLIP overlaps.
# Writes tables and run summary to data/07_coloc/.
from __future__ import annotations

import json
import os
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyBigWig

# --- config ---

ROOT_DIR = Path(os.environ["APA_ROOT"])
REF_DIR = ROOT_DIR / "ref"
RESULTS_DIR = ROOT_DIR / "data/07_coloc"
SUMMARY_FILENAME = "functional_interpretation_run_summary.json"
EXTERNAL_SC_REF_DIR = (ROOT_DIR / ".." / ".." / "effect" / "EXP" / "sc_ref").resolve()
STANDARDIZED_QTL = REF_DIR / "all_stages.susie_qtl.standardized.tsv.gz"
POLYA_DB = REF_DIR / "PolyA_DB_hg38.txt.gz"
PHASTCONS_BW = REF_DIR / "hg38.phastCons100way.bw"
POSTAR_MERGED = REF_DIR / "postar3_all_RBPs_data" / "ALL_RBP_binding_sites_merged.csv"

PIP_THRESHOLD = 0.5
CORE_THRESHOLD = 0.8
TIER2_PIP_THRESHOLD = 0.1
SINGLE_PAS_FLANK_BP = 5000
POSTAR_RBP_LIMIT = 20
POSTAR_CHUNK_SIZE = 250000

REFERENCE_SPECS = [
    {
        "reference_id": "qtl_standardized_table",
        "required": True,
        "category": "core_input",
        "analysis_modules": "candidate selection;PAS event standardization",
        "patterns": ["all_stages.susie_qtl.standardized.tsv.gz"],
        "download_source": "generated locally",
        "source_url": "",
        "notes": "Current standardized APA-QTL table.",
    },
    {
        "reference_id": "polyadb_hg38",
        "required": True,
        "category": "core_input",
        "analysis_modules": "PAS event standardization",
        "patterns": ["PolyA_DB_hg38.txt.gz"],
        "download_source": "existing local reference",
        "source_url": "",
        "notes": "PolyA_DB site catalog used for exploratory PAS pairing.",
    },
    {
        "reference_id": "phastcons_100way_bw",
        "required": True,
        "category": "core_input",
        "analysis_modules": "conservation scoring",
        "patterns": ["hg38.phastCons100way.bw"],
        "download_source": "existing local reference",
        "source_url": "",
        "notes": "Used for point and interval conservation summaries.",
    },
    {
        "reference_id": "postar3_merged_clip",
        "required": True,
        "category": "core_input",
        "analysis_modules": "RBP CLIP overlap",
        "patterns": ["postar3_all_RBPs_data/ALL_RBP_binding_sites_merged.csv"],
        "download_source": "generated locally",
        "source_url": "",
        "notes": "Merged POSTAR3 CLIP peaks.",
    },
    {
        "reference_id": "pas_motif_library",
        "required": True,
        "category": "core_input",
        "analysis_modules": "PAS motif scanning",
        "patterns": ["motif.list", "motif_databases/motif_db.csv"],
        "download_source": "existing local reference",
        "source_url": "",
        "notes": "Current local motif list and unpacked MEME motif bundle.",
    },
    {
        "reference_id": "grch38_genome_fasta",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "sequence extraction;PAS motif delta;miRNA;RNAfold",
        "patterns": [
            "hg38.fa",
            "hg38.fa.gz",
            "GRCh38*.fa",
            "GRCh38*.fa.gz",
            "*genome*.fa",
            "*genome*.fa.gz",
            "hg38*.2bit",
            "GRCh38*.2bit",
        ],
        "download_source": "GENCODE or UCSC hg38/GRCh38 genome FASTA",
        "source_url": "https://www.gencodegenes.org/human/ or https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/",
        "notes": "Expect a genome FASTA or 2bit matching hg38/GRCh38, not small motif or miRNA FASTA files.",
    },
    {
        "reference_id": "grch38_gene_annotation_gtf",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "3UTR/DUR definition;long-short isoform sequence extraction",
        "patterns": ["*.gtf", "*.gtf.gz", "*.gff3", "*.gff3.gz"],
        "download_source": "GENCODE comprehensive gene annotation",
        "source_url": "https://www.gencodegenes.org/human/",
        "notes": "Needed to define transcript-aware long and short 3UTR segments.",
    },
    {
        "reference_id": "mirna_reference",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "miRNA site gain-loss analysis",
        "patterns": ["*miRBase*", "mature.fa", "*mature*.fa*", "*mature*.fasta*", "*miRNA*"],
        "download_source": "miRBase mature miRNA sequences",
        "source_url": "https://www.mirbase.org/download/",
        "notes": "Local mature.fa is available, but it includes multiple species and should be filtered to human before analysis.",
    },
    {
        "reference_id": "targetscan_predictions",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "miRNA site gain-loss analysis",
        "patterns": [
            "*TargetScan*",
            "*targetscan*",
            "Summary_Counts.all_predictions.txt.zip",
            "Conserved_Family_Info.txt.zip",
        ],
        "download_source": "TargetScanHuman predictions/family tables",
        "source_url": "https://www.targetscan.org/vert_72/",
        "notes": "Current local TargetScan bundles are present as zipped prediction tables and require human-specific filtering before analysis.",
    },
    {
        "reference_id": "cisbp_rna_motifs",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "RBP motif-level affinity delta",
        "patterns": [
            "*CISBP*",
            "*cisbp*",
            "motif_databases/CISBP-RNA/Homo_sapiens.meme",
            "motif_databases/CISBP-RNA/Homo_sapiens.dna_encoded.meme",
        ],
        "download_source": "CISBP-RNA motif database",
        "source_url": "https://cisbp.org/rna/",
        "notes": "Human CISBP-RNA motif files are present locally and support motif-level RBP delta scanning.",
    },
    {
        "reference_id": "single_cell_brain_reference",
        "required": True,
        "category": "missing_blocker",
        "analysis_modules": "cell-type expression mapping",
        "patterns": ["*.h5ad", "*.loom", "*.mtx", "*.rds", "*.rda", "*Polioudakis*"],
        "download_source": "human fetal cortex scRNA-seq expression matrix plus metadata",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6831089/",
        "notes": "Use ../../effect/EXP/sc_ref as the default external reference root.",
        "search_roots": ["external_sc_ref"],
    },
    {
        "reference_id": "deepripe_models",
        "required": False,
        "category": "optional_enhancement",
        "analysis_modules": "model-based RBP binding prediction",
        "patterns": ["*DeepRiPe*", "*deepripe*", "DeepRiPe/*.h5"],
        "download_source": "DeepRiPe model weights",
        "source_url": "",
        "notes": "DeepRiPe model weight files are present locally as optional enhancement inputs.",
    },
    {
        "reference_id": "zipcode_motif_set",
        "required": False,
        "category": "optional_enhancement",
        "analysis_modules": "neuronal zipcode motif search",
        "patterns": ["*zipcode*", "*Zipcode*"],
        "download_source": "curated neuronal zipcode motif set",
        "source_url": "",
        "notes": "Optional exploratory analysis only.",
    },
]

# --- helpers ---


@dataclass(frozen=True)
class IntervalRecord:
    start: int
    end: int
    candidate_id: str


def write_table(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        df.to_csv(path, sep="\t", index=False, compression="gzip")
    else:
        df.to_csv(path, sep="\t", index=False)


def bw_point_value(bw: pyBigWig.pyBigWig, chrom: str, pos_1based: int) -> float | pd.NA:
    try:
        values = bw.values(chrom, pos_1based - 1, pos_1based)
    except RuntimeError:
        return pd.NA
    if not values:
        return pd.NA
    value = values[0]
    if value is None:
        return pd.NA
    return float(value)


def bw_interval_stat(
    bw: pyBigWig.pyBigWig, chrom: str, start_1based: int, end_1based: int, stat: str
) -> float | pd.NA:
    try:
        values = bw.stats(chrom, start_1based - 1, end_1based, type=stat, exact=True)
    except RuntimeError:
        return pd.NA
    if not values:
        return pd.NA
    value = values[0]
    if value is None:
        return pd.NA
    return float(value)


def build_interval_index(
    df: pd.DataFrame, start_col: str, end_col: str
) -> dict[str, tuple[list[int], list[IntervalRecord]]]:
    index: dict[str, tuple[list[int], list[IntervalRecord]]] = {}
    valid = df.loc[df[start_col].notna() & df[end_col].notna(), ["candidate_id", "phenotype_chr", start_col, end_col]]
    for chrom, group in valid.groupby("phenotype_chr", sort=False):
        intervals = [
            IntervalRecord(int(item[start_col]), int(item[end_col]), item["candidate_id"])
            for item in group.to_dict("records")
        ]
        intervals.sort(key=lambda item: (item.start, item.end, item.candidate_id))
        starts = [item.start for item in intervals]
        index[chrom] = (starts, intervals)
    return index


def overlapping_candidate_ids(
    index: dict[str, tuple[list[int], list[IntervalRecord]]],
    chrom: str,
    start: int,
    end: int,
) -> list[str]:
    payload = index.get(chrom)
    if payload is None:
        return []
    starts, intervals = payload
    stop = bisect_right(starts, end)
    hits = []
    for interval in intervals[:stop]:
        if interval.end >= start:
            hits.append(interval.candidate_id)
    return hits


# --- reference audit ---

reference_rows = []
for spec in REFERENCE_SPECS:
    roots = [REF_DIR]
    for root_name in spec.get("search_roots", []):
        if root_name == "external_sc_ref":
            roots.append(EXTERNAL_SC_REF_DIR)
    matches: list[Path] = []
    for root in roots:
        if not root.exists():
            continue
        for pattern in spec["patterns"]:
            matches.extend(path for path in root.rglob(pattern) if path.is_file())
    matches = [Path(path) for path in sorted({path.resolve() for path in matches})]
    reference_rows.append(
        {
            "reference_id": spec["reference_id"],
            "required": spec["required"],
            "category": spec["category"],
            "analysis_modules": spec["analysis_modules"],
            "status": "present" if matches else "missing",
            "n_matches": len(matches),
            "matched_paths": ";".join(
                str(path.relative_to(ROOT_DIR)) if path.is_relative_to(ROOT_DIR) else str(path)
                for path in matches[:10]
            ),
            "download_source": spec["download_source"],
            "source_url": spec["source_url"],
            "notes": spec["notes"],
        }
    )
reference_audit = pd.DataFrame(reference_rows).sort_values(
    ["required", "status", "reference_id"],
    ascending=[False, True, True],
    kind="stable",
)
write_table(reference_audit, RESULTS_DIR / "reference_audit.tsv")
write_table(
    reference_audit.loc[reference_audit["status"].eq("missing")],
    RESULTS_DIR / "reference_download_manifest.tsv",
)

# --- candidate set ---

qtl = pd.read_csv(STANDARDIZED_QTL, sep="\t")

candidates = qtl.loc[qtl["pip"].ge(PIP_THRESHOLD) | qtl["stage"].eq("interaction_GxStage")].copy()
candidates["tier"] = candidates["pip"].ge(CORE_THRESHOLD).map(
    {True: "core_high_pip", False: "expanded_mid_pip"}
)
candidates.loc[candidates["stage"].eq("interaction_GxStage"), "tier"] = "interaction_GxStage"
candidates["candidate_id"] = (
    candidates["stage"] + "::" + candidates["phenotype_id"] + "::" + candidates["variant_id"]
)
candidates["phenotype_length"] = candidates["phenotype_end"] - candidates["phenotype_start"] + 1
candidates = candidates.sort_values(
    ["tier", "stage", "pip", "gene_name", "phenotype_id", "variant_pos"],
    ascending=[True, True, False, True, True, True],
    kind="stable",
)
write_table(candidates, RESULTS_DIR / "candidate_set_pip_0.5.tsv.gz")

tier2_candidates = qtl.loc[qtl["pip"].ge(TIER2_PIP_THRESHOLD) & qtl["pip"].lt(PIP_THRESHOLD)].copy()
tier2_candidates["tier"] = "tier2_extended"
tier2_candidates["candidate_id"] = (
    tier2_candidates["stage"] + "::" + tier2_candidates["phenotype_id"] + "::" + tier2_candidates["variant_id"]
)
tier2_candidates = tier2_candidates.sort_values(
    ["stage", "pip", "gene_name", "phenotype_id", "variant_pos"],
    ascending=[True, False, True, True, True],
    kind="stable",
)
write_table(tier2_candidates, RESULTS_DIR / "tier2_candidates.tsv")

# --- PolyA_DB PAS events ---

poly = pd.read_csv(POLYA_DB, sep="\t").rename(
    columns={
        "chr": "pas_chr",
        "pA_site": "pas_pos",
        "strand": "pas_strand",
        "Gene_Symbol": "gene_name",
        "PAS_ID": "pas_id",
        "PAS_Signal": "pas_signal",
    }
)
poly["pas_pos"] = poly["pas_pos"].astype(int)
poly = poly.sort_values(
    ["gene_name", "pas_chr", "pas_strand", "pas_pos"],
    kind="stable",
)
poly_lookup: dict[tuple[str, str, str], pd.DataFrame] = {}
for key, group in poly.groupby(["gene_name", "pas_chr", "pas_strand"], sort=False):
    poly_lookup[key] = group.reset_index(drop=True)

pas_records = []
for row in candidates.itertuples(index=False):
    key = (row.gene_name, row.phenotype_chr, row.phenotype_strand)
    pas_table = poly_lookup.get(key)
    if pas_table is None:
        pas_records.append(
            {
                "candidate_id": row.candidate_id,
                "pas_resolution_status": "no_gene_chr_strand_match",
                "pas_resolution_detail": "No PolyA_DB entry for gene/chromosome/strand.",
                "pas_resolution_confidence": "none",
                "n_pas_gene_chr_strand": 0,
                "n_pas_in_region": 0,
                "n_pas_nearby_25kb": 0,
                "n_pas_distal_flank_5kb": 0,
                "proximal_pas_id": pd.NA,
                "proximal_pas_pos": pd.NA,
                "proximal_pas_signal": pd.NA,
                "distal_pas_id": pd.NA,
                "distal_pas_pos": pd.NA,
                "distal_pas_signal": pd.NA,
                "dur_start": pd.NA,
                "dur_end": pd.NA,
                "dur_length": pd.NA,
                "pas_ids_in_region": "",
                "pas_signals_in_region": "",
            }
        )
        continue

    in_region = pas_table.loc[
        (pas_table["pas_pos"] >= row.phenotype_start)
        & (pas_table["pas_pos"] <= row.phenotype_end)
    ].copy()
    nearby = pas_table.loc[
        (pas_table["pas_pos"] >= row.phenotype_start - 25000)
        & (pas_table["pas_pos"] <= row.phenotype_end + 25000)
    ].copy()

    status = "unresolved_no_pas_in_region"
    detail = "No PolyA_DB PAS fell inside the phenotype interval."
    confidence = "none"
    proximal = None
    distal = None
    distal_flank = pd.DataFrame()

    if len(in_region) >= 2:
        if row.phenotype_strand == "+":
            ordered = in_region.sort_values("pas_pos", ascending=True, kind="stable")
        else:
            ordered = in_region.sort_values("pas_pos", ascending=False, kind="stable")
        proximal = ordered.iloc[0]
        distal = ordered.iloc[-1]
        status = "resolved_region_extrema"
        detail = "Using the transcript-ordered extrema among in-region PolyA_DB sites."
        confidence = "moderate"
    elif len(in_region) == 1:
        anchor = in_region.iloc[0]["pas_pos"]
        if row.phenotype_strand == "+":
            distal_mask = pas_table["pas_pos"] > anchor
        else:
            distal_mask = pas_table["pas_pos"] < anchor
        distal_flank = pas_table.loc[
            distal_mask & (pas_table["pas_pos"].sub(anchor).abs() <= SINGLE_PAS_FLANK_BP)
        ].copy()
        if len(distal_flank) == 1:
            chosen = distal_flank.sort_values("pas_pos", kind="stable").iloc[0]
            proximal, distal = in_region.iloc[0], chosen
            status = "resolved_single_pas_plus_distal_flank"
            detail = (
                "One in-region PAS plus one transcript-distal flank PAS within 5 kb."
            )
            confidence = "low"
        elif len(distal_flank) > 1:
            status = "ambiguous_single_pas_multiple_distal_flanks"
            detail = (
                "Exactly one in-region PAS, but multiple distal flank candidates within 5 kb."
            )
            confidence = "none"
        else:
            status = "unresolved_single_pas_in_region"
            detail = "Exactly one in-region PAS and no usable distal flank PAS within 5 kb."
            confidence = "none"

    dur_start = pd.NA
    dur_end = pd.NA
    dur_length = pd.NA
    proximal_pas_id = pd.NA
    proximal_pas_pos = pd.NA
    proximal_pas_signal = pd.NA
    distal_pas_id = pd.NA
    distal_pas_pos = pd.NA
    distal_pas_signal = pd.NA

    if proximal is not None and distal is not None:
        proximal_pas_id = proximal["pas_id"]
        proximal_pas_pos = int(proximal["pas_pos"])
        proximal_pas_signal = proximal["pas_signal"]
        distal_pas_id = distal["pas_id"]
        distal_pas_pos = int(distal["pas_pos"])
        distal_pas_signal = distal["pas_signal"]
        dur_start = int(min(proximal_pas_pos, distal_pas_pos))
        dur_end = int(max(proximal_pas_pos, distal_pas_pos))
        dur_length = int(abs(distal_pas_pos - proximal_pas_pos))

    pas_records.append(
        {
            "candidate_id": row.candidate_id,
            "pas_resolution_status": status,
            "pas_resolution_detail": detail,
            "pas_resolution_confidence": confidence,
            "n_pas_gene_chr_strand": int(len(pas_table)),
            "n_pas_in_region": int(len(in_region)),
            "n_pas_nearby_25kb": int(len(nearby)),
            "n_pas_distal_flank_5kb": int(len(distal_flank)),
            "proximal_pas_id": proximal_pas_id,
            "proximal_pas_pos": proximal_pas_pos,
            "proximal_pas_signal": proximal_pas_signal,
            "distal_pas_id": distal_pas_id,
            "distal_pas_pos": distal_pas_pos,
            "distal_pas_signal": distal_pas_signal,
            "dur_start": dur_start,
            "dur_end": dur_end,
            "dur_length": dur_length,
            "pas_ids_in_region": ";".join(in_region["pas_id"].astype(str).tolist()[:20]),
            "pas_signals_in_region": ";".join(
                in_region["pas_signal"].astype(str).tolist()[:20]
            ),
        }
    )

pas_annotations = pd.DataFrame.from_records(pas_records)
candidate_with_pas = candidates.merge(pas_annotations, on="candidate_id", how="left")
write_table(candidate_with_pas, RESULTS_DIR / "pas_event_annotations.tsv.gz")

# --- phastCons conservation ---

bw = pyBigWig.open(str(PHASTCONS_BW))
conservation_rows = []
try:
    for row in candidate_with_pas.itertuples(index=False):
        record = {
            "candidate_id": row.candidate_id,
            "variant_phastcons": bw_point_value(bw, row.variant_chr, int(row.variant_pos)),
            "phenotype_mean_phastcons": bw_interval_stat(
                bw,
                row.phenotype_chr,
                int(row.phenotype_start),
                int(row.phenotype_end),
                "mean",
            ),
            "phenotype_max_phastcons": bw_interval_stat(
                bw,
                row.phenotype_chr,
                int(row.phenotype_start),
                int(row.phenotype_end),
                "max",
            ),
            "proximal_pas_phastcons": pd.NA,
            "distal_pas_phastcons": pd.NA,
            "dur_mean_phastcons": pd.NA,
            "dur_max_phastcons": pd.NA,
        }
        if pd.notna(row.proximal_pas_pos):
            record["proximal_pas_phastcons"] = bw_point_value(
                bw, row.phenotype_chr, int(row.proximal_pas_pos)
            )
        if pd.notna(row.distal_pas_pos):
            record["distal_pas_phastcons"] = bw_point_value(
                bw, row.phenotype_chr, int(row.distal_pas_pos)
            )
        if pd.notna(row.dur_start) and pd.notna(row.dur_end) and int(row.dur_end) > int(row.dur_start):
            record["dur_mean_phastcons"] = bw_interval_stat(
                bw, row.phenotype_chr, int(row.dur_start), int(row.dur_end), "mean"
            )
            record["dur_max_phastcons"] = bw_interval_stat(
                bw, row.phenotype_chr, int(row.dur_start), int(row.dur_end), "max"
            )
        conservation_rows.append(record)
finally:
    bw.close()
conservation = pd.DataFrame.from_records(conservation_rows)
candidate_with_conservation = candidate_with_pas.merge(conservation, on="candidate_id", how="left")

# --- POSTAR3 CLIP overlaps ---

phenotype_index = build_interval_index(candidate_with_conservation, "phenotype_start", "phenotype_end")
variant_df = candidate_with_conservation.assign(
    variant_start=candidate_with_conservation["variant_pos"],
    variant_end=candidate_with_conservation["variant_pos"],
)
variant_index = build_interval_index(variant_df, "variant_start", "variant_end")
dur_index = build_interval_index(candidate_with_conservation, "dur_start", "dur_end")
candidate_chroms = set(candidate_with_conservation["phenotype_chr"].unique())

candidate_ids = candidate_with_conservation["candidate_id"]
phenotype_rbps = {candidate_id: set() for candidate_id in candidate_ids}
phenotype_techs = {candidate_id: set() for candidate_id in candidate_ids}
phenotype_counts = {candidate_id: 0 for candidate_id in candidate_ids}
variant_rbps = {candidate_id: set() for candidate_id in candidate_ids}
variant_techs = {candidate_id: set() for candidate_id in candidate_ids}
variant_counts = {candidate_id: 0 for candidate_id in candidate_ids}
dur_rbps = {candidate_id: set() for candidate_id in candidate_ids}
dur_techs = {candidate_id: set() for candidate_id in candidate_ids}
dur_counts = {candidate_id: 0 for candidate_id in candidate_ids}

chunk_iter = pd.read_csv(
    POSTAR_MERGED,
    usecols=[0, 4, 20, 21],
    chunksize=POSTAR_CHUNK_SIZE,
    encoding="utf-8-sig",
)

total_rows = 0
for chunk in chunk_iter:
    total_rows += len(chunk)
    chunk = chunk.dropna(subset=["RBP", "l1"]).copy()
    if chunk.empty:
        continue
    coords = chunk["l1"].astype(str).str.extract(
        r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$"
    )
    chunk = pd.concat([chunk.reset_index(drop=True), coords], axis=1)
    chunk = chunk.dropna(subset=["chrom", "start", "end"])
    if chunk.empty:
        continue
    chunk["chrom"] = chunk["chrom"].astype(str)
    chunk = chunk.loc[chunk["chrom"].isin(candidate_chroms)].copy()
    if chunk.empty:
        continue
    chunk["start"] = chunk["start"].astype(int)
    chunk["end"] = chunk["end"].astype(int)

    for row in chunk.itertuples(index=False):
        # Some merged POSTAR3 rows retain coordinates but lose the RBP label fields.
        if pd.notna(row.RBP) and str(row.RBP).strip():
            rbp = str(row.RBP).strip()
        elif pd.notna(row.RBP_Query) and str(row.RBP_Query).strip():
            rbp = str(row.RBP_Query).strip()
        else:
            rbp = "unknown_rbp"
        tech = (
            str(row.CLIP_Technology).strip()
            if pd.notna(row.CLIP_Technology) and str(row.CLIP_Technology).strip()
            else "unknown_technology"
        )
        phenotype_hits = overlapping_candidate_ids(
            phenotype_index, row.chrom, int(row.start), int(row.end)
        )
        for candidate_id in phenotype_hits:
            phenotype_counts[candidate_id] += 1
            phenotype_rbps[candidate_id].add(rbp)
            phenotype_techs[candidate_id].add(tech)

        variant_hits = overlapping_candidate_ids(
            variant_index, row.chrom, int(row.start), int(row.end)
        )
        for candidate_id in variant_hits:
            variant_counts[candidate_id] += 1
            variant_rbps[candidate_id].add(rbp)
            variant_techs[candidate_id].add(tech)

        dur_hits = overlapping_candidate_ids(
            dur_index, row.chrom, int(row.start), int(row.end)
        )
        for candidate_id in dur_hits:
            dur_counts[candidate_id] += 1
            dur_rbps[candidate_id].add(rbp)
            dur_techs[candidate_id].add(tech)

postar_rows = []
for candidate_id in candidate_ids:
    postar_rows.append(
        {
            "candidate_id": candidate_id,
            "postar_variant_overlap_count": variant_counts[candidate_id],
            "postar_variant_rbps": ";".join(
                sorted(variant_rbps[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
            "postar_variant_technologies": ";".join(
                sorted(variant_techs[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
            "postar_phenotype_overlap_count": phenotype_counts[candidate_id],
            "postar_phenotype_rbps": ";".join(
                sorted(phenotype_rbps[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
            "postar_phenotype_technologies": ";".join(
                sorted(phenotype_techs[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
            "postar_dur_overlap_count": dur_counts[candidate_id],
            "postar_dur_rbps": ";".join(
                sorted(dur_rbps[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
            "postar_dur_technologies": ";".join(
                sorted(dur_techs[candidate_id])[:POSTAR_RBP_LIMIT]
            ),
        }
    )
postar = pd.DataFrame.from_records(postar_rows)
final_df = candidate_with_conservation.merge(postar, on="candidate_id", how="left")
write_table(final_df, RESULTS_DIR / "functional_support_annotations.tsv.gz")

# --- summaries ---

tier_summary = (
    final_df.groupby(["tier", "stage"])
    .size()
    .reset_index(name="n_candidates")
    .sort_values(["tier", "stage"], kind="stable")
)
write_table(tier_summary, RESULTS_DIR / "candidate_counts_by_tier.tsv")

pas_summary = (
    final_df.groupby(["tier", "pas_resolution_status"])
    .size()
    .reset_index(name="n_candidates")
    .sort_values(["tier", "n_candidates", "pas_resolution_status"], ascending=[True, False, True])
)
write_table(pas_summary, RESULTS_DIR / "pas_resolution_summary.tsv")

tier_counts = candidates.groupby(["tier", "stage"]).size().reset_index(name="n")
pas_counts = (
    final_df.groupby(["tier", "pas_resolution_status"])
    .size()
    .reset_index(name="n")
    .sort_values(["tier", "n", "pas_resolution_status"], ascending=[True, False, True])
)
overlap_counts = {
    "variant_overlap_positive": int((final_df["postar_variant_overlap_count"] > 0).sum()),
    "phenotype_overlap_positive": int(
        (final_df["postar_phenotype_overlap_count"] > 0).sum()
    ),
    "dur_overlap_positive": int((final_df["postar_dur_overlap_count"] > 0).sum()),
}
summary_payload = {
    "pip_threshold": PIP_THRESHOLD,
    "core_threshold": CORE_THRESHOLD,
    "n_candidates": int(len(candidates)),
    "tier_counts": tier_counts.to_dict(orient="records"),
    "pas_resolution_counts": pas_counts.to_dict(orient="records"),
    "overlap_counts": overlap_counts,
    "missing_required_references": reference_audit.loc[
        reference_audit["required"] & reference_audit["status"].eq("missing"),
        "reference_id",
    ].tolist(),
}
with (RESULTS_DIR / SUMMARY_FILENAME).open("w", encoding="utf-8") as handle:
    json.dump(summary_payload, handle, indent=2, ensure_ascii=True)
