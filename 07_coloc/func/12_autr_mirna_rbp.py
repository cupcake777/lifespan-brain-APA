#!/usr/bin/env python3
# Annotate candidate aUTR sites with miRNA seed matches, CISBP RBP motifs and scAPAatlas Hippocampus overlap.
# Writes the site-level annotation table and a per-site count summary.

# --- packages ---
import gzip
import os
import sys
import tarfile
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pandas as pd


# --- config ---
ROOT = Path(os.environ["APA_ROOT"])
DATA = ROOT / "data/07_coloc"
GROUP_TABLE = ROOT / "data/05_gxs/intx_hc3_core_groups.tsv"
FOCAL_TABLE = ROOT / "data/03_omics/04_cv_ranked_pairs.tsv"
INTERVAL_FASTA = ROOT / "ref/hg38_autr_candidate_intervals.fa"
MATURE_FA = ROOT / "ref/mature.fa"
CISBP_MEME = ROOT / "ref/motif_databases/CISBP-RNA/Homo_sapiens.meme"
MIRBASE_GFF = ROOT / "ref/hsa.mirbase.v23.1.gff3"
GENCODE_GTF = ROOT / "ref/gencode.v40.annotation.gtf.gz"
PFC_EXPRESSION = DATA / "pfc_rbp_expression.tsv.gz"
SCAPA_DIR = Path(f"{os.environ['APA_HPC_HOME']}/project/.ref/scAPAatlas")
SCAPA_APA_TAR = SCAPA_DIR / "specificAPA_Human.tar.gz"
SCAPA_MIRNA_TAR = SCAPA_DIR / "miRNA_Human.tar.gz"
SCAPA_RBP_TAR = SCAPA_DIR / "RBP_Human.tar.gz"
OUT_ANNOTATION = DATA / "12_autr_mirna_rbp_annotation.tsv"
OUT_SUMMARY = DATA / "12_autr_annotation_summary.tsv"
FOCAL_GENES = {"INPP4A", "KPNA4", "PIK3CA", "ARIH1", "RHOT1", "NAPRT"}
MIN_RBP_TPM = 1.0
MAX_POSITIONS_SHOWN = 12
ATP2B1_SITE = "ENST00000551310.1|ATP2B1|chr12:89655583-89656107|-"


# --- helpers ---
def find_positions(sequence, pattern):
    positions = []
    start = 0
    while True:
        index = sequence.find(pattern, start)
        if index == -1:
            return positions
        positions.append(index + 1)
        start = index + 1


def compact_positions(positions):
    shown = ",".join(str(value) for value in positions[:MAX_POSITIONS_SHOWN])
    if len(positions) > MAX_POSITIONS_SHOWN:
        shown += f",...(+{len(positions) - MAX_POSITIONS_SHOWN})"
    return shown


def read_tar_csv(tar_path, member):
    with tarfile.open(tar_path, "r:gz") as archive:
        member_obj = next(item for item in archive.getmembers() if item.name == member)
        return pd.read_csv(archive.extractfile(member_obj))


def base_annotation(record, annotation_type, annotation_name):
    return {
        "site_id": record["site_id"],
        "gene_symbol": record["gene_symbol"],
        "candidate_set": record["candidate_set"],
        "interaction_group_id": record["interaction_group_id"],
        "stage_group_status": record["stage_group_status"],
        "transcript_id": record["transcript_id"],
        "chrom": record["chrom"],
        "autr_start": record["autr_start"],
        "autr_end": record["autr_end"],
        "strand": record["strand"],
        "autr_length": record["autr_length"],
        "annotation_type": annotation_type,
        "annotation_name": annotation_name,
    }


# --- load run plan helpers ---
spec = spec_from_file_location("autr_plan_03", Path(__file__).resolve().parent / "03_run_plan.py")
plan = module_from_spec(spec)
sys.modules["autr_plan_03"] = plan
spec.loader.exec_module(plan)

# --- candidate events ---
groups = pd.read_csv(
    GROUP_TABLE,
    sep="\t",
    usecols=["group_id", "gene_id", "stage_group_status", "lead_phenotype_id"],
).rename(columns={"gene_id": "gene_symbol", "lead_phenotype_id": "site_id"})
groups["candidate_set"] = "interaction_group"

focal = pd.read_csv(FOCAL_TABLE, sep="\t")
focal = focal.loc[
    focal["gene_symbol"].isin(FOCAL_GENES)
    & focal["frozen_primary_candidate"].astype(str).str.upper().eq("TRUE"),
    ["site_id", "gene_symbol"],
].copy()
focal["candidate_set"] = "protein_focal"
focal["group_id"] = ""
focal["stage_group_status"] = ""

events = pd.concat(
    [groups[["site_id", "gene_symbol", "candidate_set", "group_id", "stage_group_status"]], focal],
    ignore_index=True,
)
events = (
    events.groupby("site_id", as_index=False, sort=False)
    .agg(
        gene_symbol=("gene_symbol", "first"),
        candidate_set=("candidate_set", lambda values: ";".join(sorted(set(values)))),
        interaction_group_id=("group_id", "first"),
        stage_group_status=("stage_group_status", "first"),
    )
)

# parse site_id = transcript|gene|chrom:start-end|strand
parsed_rows = []
for site_id in events["site_id"]:
    transcript_id, gene_symbol, region, strand = site_id.split("|")
    chrom, span = region.split(":")
    start, end = span.split("-")
    parsed_rows.append(
        {
            "transcript_id": transcript_id,
            "chrom": chrom,
            "autr_start": int(start),
            "autr_end": int(end),
            "strand": strand,
        }
    )
events = pd.concat([events, pd.DataFrame(parsed_rows)], axis=1)
events["autr_length"] = events["autr_end"] - events["autr_start"] + 1

assert len(events) == 201 and events["site_id"].is_unique
assert ATP2B1_SITE in set(events["site_id"])

# --- aUTR sequences ---
genome = {}
header = None
with INTERVAL_FASTA.open(encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if line.startswith(">"):
            header = line[1:]
            genome[header] = ""
        elif header is not None:
            genome[header] += line.upper()

sequences = {}
for record in events.to_dict("records"):
    key = f"{record['chrom']}:{record['autr_start']}-{record['autr_end']}"
    sequence = genome[key]
    assert len(sequence) == record["autr_length"]
    sequences[record["site_id"]] = sequence if record["strand"] == "+" else plan.revcomp(sequence)

# --- miRNA host gene map ---
primary_loci = {}
with MIRBASE_GFF.open(encoding="utf-8") as handle:
    for raw_line in handle:
        if raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "miRNA_primary_transcript":
            continue
        attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
        primary_loci[attrs["ID"]] = (fields[0], int(fields[3]), int(fields[4]))

genes_by_chrom = {}
with gzip.open(GENCODE_GTF, "rt", encoding="utf-8") as handle:
    for raw_line in handle:
        if raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        attrs = dict(
            item.strip().split(" ", 1) for item in fields[8].strip().split(";") if item.strip()
        )
        gene_name = attrs.get("gene_name", "").strip('"')
        if gene_name:
            genes_by_chrom.setdefault(fields[0], []).append((int(fields[3]), int(fields[4]), gene_name))
for intervals in genes_by_chrom.values():
    intervals.sort()

host_by_primary = {}
for primary_id, (chrom, start, end) in primary_loci.items():
    for gene_start, gene_end, gene_name in genes_by_chrom.get(chrom, []):
        if gene_start <= start and end <= gene_end:
            host_by_primary[primary_id] = gene_name
            break

mirna_host_map = {}
with MIRBASE_GFF.open(encoding="utf-8") as handle:
    for raw_line in handle:
        if raw_line.startswith("#"):
            continue
        fields = raw_line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "miRNA":
            continue
        attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
        host = host_by_primary.get(attrs.get("Derives_from", ""))
        if host:
            mirna_host_map[attrs["Name"]] = host

# --- miRNA seed scan ---
mirna_families = plan.parse_human_mirnas(MATURE_FA)
expression = (
    pd.read_csv(PFC_EXPRESSION, sep="\t")
    .drop_duplicates("gene_symbol")
    .set_index("gene_symbol")[["prenatal_mean_tpm", "postnatal_mean_tpm"]]
)
annotation_rows = []

for record in events.to_dict("records"):
    sequence = sequences[record["site_id"]]
    for family in mirna_families:
        hit_types = []
        hit_positions = []
        total_sites = 0
        for label, pattern in (
            ("8mer", family["site_8mer"]),
            ("7mer-m8", family["site_7mer_m8"]),
            ("7mer-a1", family["site_7mer_a1"]),
        ):
            positions = find_positions(sequence, pattern)
            if positions:
                hit_types.append(f"{label}:{len(positions)}")
                hit_positions.append(f"{label}:{compact_positions(positions)}")
                total_sites += len(positions)
        if total_sites == 0:
            continue
        host_genes = sorted({mirna_host_map[member] for member in family["members"] if member in mirna_host_map})
        host_prenatal = max(
            (float(expression.loc[gene, "prenatal_mean_tpm"]) for gene in host_genes if gene in expression.index),
            default=float("nan"),
        )
        host_postnatal = max(
            (float(expression.loc[gene, "postnatal_mean_tpm"]) for gene in host_genes if gene in expression.index),
            default=float("nan"),
        )
        host_expressed = host_prenatal >= MIN_RBP_TPM or host_postnatal >= MIN_RBP_TPM
        row = base_annotation(record, "mirna_seed_match", family["representative_mirna"])
        row.update(
            {
                "annotation_id": family["family_id"],
                "n_sites": total_sites,
                "site_types": ";".join(hit_types),
                "positions": ";".join(hit_positions),
                "cell_type": "",
                "score": "",
                "eclip_support": "",
                "prenatal_mean_tpm": "",
                "postnatal_mean_tpm": "",
                "host_genes": ";".join(host_genes[:4]) + (f";+{len(host_genes) - 4}" if len(host_genes) > 4 else ""),
                "host_prenatal_max_tpm": host_prenatal,
                "host_postnatal_max_tpm": host_postnatal,
                "expression_context": "PFC_miRNA_expression_reference_not_available",
                "evidence_tier": "host_gene_expression_proxy" if host_expressed else "sequence_only",
                "strict_mirna_candidate": bool(hit_types and hit_types[0].startswith("8mer")) and host_expressed,
                "source": "miRBase_mature_hsa_seed_scan_hg38_aUTR",
            }
        )
        annotation_rows.append(row)

# --- RBP motif scan ---
cisbp_motifs = plan.parse_cisbp_meme(CISBP_MEME)

for record in events.to_dict("records"):
    sequence = sequences[record["site_id"]]
    for motif in cisbp_motifs:
        positions = find_positions(sequence, motif["consensus_dna"])
        if not positions:
            continue
        name = str(motif["rbp_name"])
        symbol = name[1:name.index(")")] if name.startswith("(") and ")" in name else name
        if symbol in expression.index:
            prenatal_tpm = float(expression.loc[symbol, "prenatal_mean_tpm"])
            postnatal_tpm = float(expression.loc[symbol, "postnatal_mean_tpm"])
        else:
            prenatal_tpm = postnatal_tpm = float("nan")
        expression_supported = prenatal_tpm >= MIN_RBP_TPM or postnatal_tpm >= MIN_RBP_TPM
        row = base_annotation(record, "rbp_motif", symbol)
        row.update(
            {
                "annotation_id": motif["motif_id"],
                "n_sites": len(positions),
                "site_types": f"consensus:{motif['consensus_dna']}",
                "positions": compact_positions(positions),
                "cell_type": "",
                "score": "",
                "eclip_support": "",
                "prenatal_mean_tpm": prenatal_tpm,
                "postnatal_mean_tpm": postnatal_tpm,
                "host_genes": "",
                "host_prenatal_max_tpm": "",
                "host_postnatal_max_tpm": "",
                "expression_context": "own_cohort_PFC_RNAseq_mean_TPM",
                "evidence_tier": "expression_supported_rbp" if expression_supported else "sequence_only",
                "strict_mirna_candidate": "",
                "source": "CISBP_RNA_consensus_scan_with_pfc_expression",
            }
        )
        annotation_rows.append(row)

# --- scAPAatlas Hippocampus overlap ---
atlas_apa = read_tar_csv(SCAPA_APA_TAR, "specificAPA_Human/Human_APA_Hippocampus.csv")
gene_events = {gene: [] for gene in events["gene_symbol"].unique()}
for record in events.to_dict("records"):
    gene_events[record["gene_symbol"]].append(record)

gene_matched = atlas_apa.loc[atlas_apa["Gene_Symbol"].isin(gene_events)].copy()
matched_apa_rows = []
for atlas_record in gene_matched.to_dict("records"):
    coordinate = str(atlas_record["Specific_APA_Coordinates"])
    chrom, position_text, atlas_strand = coordinate.split(":")
    atlas_position = int(position_text)
    for record in gene_events[atlas_record["Gene_Symbol"]]:
        if record["chrom"] != chrom or record["strand"] != atlas_strand:
            continue
        inside = record["autr_start"] <= atlas_position <= record["autr_end"]
        if not inside:
            continue
        row = base_annotation(record, "atlas_specific_apa", str(atlas_record["Specific_APA"]))
        row.update(
            {
                "annotation_id": coordinate,
                "n_sites": 1,
                "site_types": "atlas_cell_type_specific_APA_in_aUTR",
                "positions": coordinate,
                "cell_type": atlas_record["Cell_Type"],
                "score": atlas_record["PAU"],
                "eclip_support": "",
                "prenatal_mean_tpm": "",
                "postnatal_mean_tpm": "",
                "expression_context": "adult_Hippocampus_cross_tissue_cross_developmentmental",
                "evidence_tier": "atlas_overlap",
                "source": "scAPAatlas_Hippocampus_specificAPA",
            }
        )
        annotation_rows.append(row)
        matched_apa_rows.append({**atlas_record, "matched_site_id": record["site_id"]})

if matched_apa_rows:
    matched_apa = pd.DataFrame(matched_apa_rows)
    join_keys = ["Gene_Symbol", "Specific_APA", "Specific_APA_Coordinates", "Cell_Type"]
    atlas_mirna = read_tar_csv(SCAPA_MIRNA_TAR, "miRNA_Human/Human_miRNA_Hippocampus.csv")
    atlas_rbp = read_tar_csv(SCAPA_RBP_TAR, "RBP_Human/Human_RBP_Hippocampus.csv")
    record_by_site = events.set_index("site_id", drop=False).to_dict("index")
    for annotation_type, atlas_table, name_column, support_column in (
        ("atlas_mirna", atlas_mirna, "miRNA", "Score"),
        ("atlas_rbp", atlas_rbp, "RBP", "eCLIP_Support"),
    ):
        joined = matched_apa.merge(atlas_table, on=join_keys, how="inner")
        for hit in joined.to_dict("records"):
            record = record_by_site[hit["matched_site_id"]]
            row = base_annotation(record, annotation_type, str(hit[name_column]))
            tier = (
                "atlas_clip_supported"
                if annotation_type == "atlas_rbp" and int(hit["eCLIP_Support"]) == 1
                else "atlas_sequence_only"
            )
            row.update(
                {
                    "annotation_id": str(hit["Specific_APA"]),
                    "n_sites": 1,
                    "site_types": "atlas_annotated_binding_site",
                    "positions": str(hit[f"{name_column}_Binding_Sites"]),
                    "cell_type": hit["Cell_Type"],
                    "score": hit[support_column] if annotation_type == "atlas_mirna" else "",
                    "eclip_support": hit[support_column] if annotation_type == "atlas_rbp" else "",
                    "prenatal_mean_tpm": "",
                    "postnatal_mean_tpm": "",
                    "expression_context": "adult_Hippocampus_cross_tissue_cross_developmentmental",
                    "evidence_tier": tier,
                    "source": "scAPAatlas_Hippocampus_" + ("miRNA" if annotation_type == "atlas_mirna" else "RBP"),
                }
            )
            annotation_rows.append(row)

# --- summary and outputs ---
annotations = pd.DataFrame(annotation_rows)

summary = (
    annotations.pivot_table(
        index=["site_id", "gene_symbol", "candidate_set", "autr_length"],
        columns="annotation_type",
        values="annotation_name",
        aggfunc="count",
    )
    .reset_index()
)
summary.columns = [
    f"n_{column}" if column not in {"site_id", "gene_symbol", "candidate_set", "autr_length"} else column
    for column in summary.columns
]
events_base = events[["site_id", "gene_symbol", "candidate_set", "autr_length"]]
summary = events_base.merge(summary, on=["site_id", "gene_symbol", "candidate_set", "autr_length"], how="left")
count_columns = [column for column in summary.columns if column.startswith("n_")]
summary[count_columns] = summary[count_columns].fillna(0).astype(int)
expression_supported_counts = (
    annotations.loc[lambda table: table["evidence_tier"].eq("expression_supported_rbp")]
    .groupby("site_id")["annotation_name"].nunique().rename("n_expression_supported_rbps")
)
clip_counts = (
    annotations.loc[lambda table: table["evidence_tier"].eq("atlas_clip_supported")]
    .groupby("site_id")["annotation_name"].nunique().rename("n_atlas_clip_supported_rbps")
)
strict_counts = (
    annotations.loc[lambda table: table["strict_mirna_candidate"].astype(str).str.lower().eq("true")]
    .groupby("site_id")["annotation_name"].nunique().rename("n_strict_mirna_candidates")
)
summary = (
    summary.merge(expression_supported_counts, on="site_id", how="left")
    .merge(clip_counts, on="site_id", how="left")
    .merge(strict_counts, on="site_id", how="left")
)
summary[["n_expression_supported_rbps", "n_atlas_clip_supported_rbps", "n_strict_mirna_candidates"]] = summary[
    ["n_expression_supported_rbps", "n_atlas_clip_supported_rbps", "n_strict_mirna_candidates"]
].fillna(0).astype(int)

annotations.to_csv(OUT_ANNOTATION, sep="\t", index=False)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
