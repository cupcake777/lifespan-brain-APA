#!/usr/bin/env python3
# miRNA seed-match feature track for one APA candidate (DUR, DnaFeaturesViewer).
# Plotted-site and layout tables go to data/07_coloc; PNG/PDF to figures/07_coloc.

from __future__ import annotations

import os
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

import matplotlib.pyplot as plt
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
RESULTS_DIR = ROOT / "data" / "07_coloc"
FIG_DIR = ROOT / "figures" / "07_coloc"
TARGETSCAN_FAMILY_ZIP = ROOT / "ref" / "Conserved_Family_Info.txt.zip"
MATURE_FA = ROOT / "ref" / "mature.fa"

# "auto" chooses the highest-ranked candidate with density-supported miRNA.
# Set to a full candidate_id string for a specific locus.
CANDIDATE_ID = "auto"

# Prefer a visually miRNA-focused candidate rather than the absolute top hit
# when possible.
PREFER_MIRNA_SINGLE_MODULE = True

TOP_MIRNA_FAMILIES = 10
MAX_SITES_PER_FAMILY = 12
OUT_PREFIX = FIG_DIR / "figure_mirna_mechanism_track"
TABLE_PREFIX = RESULTS_DIR / "figure_mirna_mechanism_track"
LAYOUT_VALIDATION_PATH = RESULTS_DIR / "figure_mirna_mechanism_track.layout_validation.tsv"

SITE_STYLE = {
    "8mer": {"color": "#c0392b", "height": 18},
    "7mer-m8": {"color": "#d9822b", "height": 14},
    "7mer-a1": {"color": "#2f80ed", "height": 12},
}


# --- helpers ---

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks).upper().replace("U", "T")
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            records[header] = "".join(chunks).upper().replace("U", "T")
    return records


# --- candidate selection ---

ranked = pd.read_csv(RESULTS_DIR / "ranked_candidates.tsv", sep="\t")
func = pd.read_csv(RESULTS_DIR / "functional_effect_main_table.tsv.gz", sep="\t")
mech = pd.read_csv(RESULTS_DIR / "mechanism_main_table.tsv.gz", sep="\t")
ann = pd.read_csv(RESULTS_DIR / "utr_sequence_annotations.tsv.gz", sep="\t")

meta = ranked.merge(
    func[
        [
            "candidate_id",
            "mirna_total_sites",
            "mirna_unique_families",
            "mirna_density_threshold_per_100nt",
            "top_mirna_families",
        ]
    ],
    on="candidate_id",
    how="left",
)
meta = meta.merge(
    mech[["candidate_id", "effect_direction_status", "utr_effect_direction"]],
    on="candidate_id",
    how="left",
)
meta = meta.merge(
    ann[["candidate_id", "short_utr_length", "long_utr_length", "dur_seq_length", "transcript_strand"]],
    on="candidate_id",
    how="left",
)

if CANDIDATE_ID != "auto":
    selected = meta.loc[meta["candidate_id"].eq(CANDIDATE_ID)]
    if selected.empty:
        raise SystemExit(f"Configured CANDIDATE_ID not found: {CANDIDATE_ID}")
    candidate_id = CANDIDATE_ID
    row = selected.iloc[0]
else:
    candidates = meta[
        meta["density_supported_modules"].fillna("").str.contains("miRNA", regex=False)
        & meta["short_utr_length"].notna()
        & meta["dur_seq_length"].notna()
        & meta["mirna_total_sites"].fillna(0).gt(0)
    ].copy()
    if candidates.empty:
        raise SystemExit("No density-supported miRNA candidate found in existing results.")
    if PREFER_MIRNA_SINGLE_MODULE:
        single = candidates[candidates["density_supported_modules"].eq("miRNA")].copy()
        if not single.empty:
            candidates = single
    row = candidates.sort_values(["rank", "TotalScore"], ascending=[True, False]).iloc[0]
    candidate_id = str(row["candidate_id"])

dur_records = parse_fasta(RESULTS_DIR / "dur.fa")
if candidate_id not in dur_records:
    raise SystemExit(f"Candidate missing from results/dur.fa: {candidate_id}")


# --- miRNA families: TargetScan seeds, else mature miRNA seeds ---

families: list[dict] = []
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

        targetscan: dict[str, dict] = {}
        with ZipFile(TARGETSCAN_FAMILY_ZIP) as zf:
            member = zf.namelist()[0]
            with zf.open(member) as handle:
                for family_df in pd.read_csv(handle, sep="\t", usecols=usecols, chunksize=200000):
                    if species_col is not None:
                        species_text = family_df[species_col].astype(str).str.lower()
                        human_mask = species_text.str.contains("human|hsa|9606", regex=True, na=False)
                        if human_mask.any():
                            family_df = family_df.loc[human_mask].copy()
                    for _, fam in family_df.iterrows():
                        seed = str(fam[seed_col]).upper().replace("U", "T").strip()
                        if len(seed) < 6 or set(seed) - {"A", "C", "G", "T"}:
                            continue
                        seed7 = seed[:7]
                        family_id = str(fam[family_col]).strip() if family_col is not None else f"TargetScan_{seed7}"
                        entry = targetscan.setdefault(
                            seed7,
                            {
                                "seed": seed7,
                                "family_id": family_id,
                                "representative_mirna": family_id,
                                "members": [],
                            },
                        )
                        entry["members"].append(family_id)
        for entry in targetscan.values():
            seed_rc = revcomp(entry["seed"])
            sixmer_rc = revcomp(entry["seed"][:6])
            entry["site_8mer"] = "A" + seed_rc
            entry["site_7mer_m8"] = seed_rc
            entry["site_7mer_a1"] = "A" + sixmer_rc
            entry["members"] = sorted(set(entry["members"]))
        families = sorted(targetscan.values(), key=lambda item: item["representative_mirna"])

if not families:
    mature_records = parse_fasta(MATURE_FA)
    mirna_families: dict[str, dict] = {}
    for header_text, seq in mature_records.items():
        mirna_name = header_text.split()[0]
        if not mirna_name.startswith("hsa-") or len(seq) < 8:
            continue
        seed = seq[1:8]
        entry = mirna_families.setdefault(
            seed,
            {
                "seed": seed,
                "family_id": f"seed_{seed}",
                "representative_mirna": mirna_name,
                "members": [],
            },
        )
        entry["members"].append(mirna_name)
    for entry in mirna_families.values():
        seed_rc = revcomp(entry["seed"])
        sixmer_rc = revcomp(entry["seed"][:6])
        entry["site_8mer"] = "A" + seed_rc
        entry["site_7mer_m8"] = seed_rc
        entry["site_7mer_a1"] = "A" + sixmer_rc
        entry["members"] = sorted(set(entry["members"]))
    families = sorted(mirna_families.values(), key=lambda item: item["representative_mirna"])

selected_names = [
    item.split(":", 1)[0]
    for item in str(row.get("top_mirna_families", "") or "").split(";")
    if item
][:TOP_MIRNA_FAMILIES]


# --- seed-match sites on the DUR ---

sequence = dur_records[candidate_id]
short_len = int(row["short_utr_length"])
wanted = set(selected_names)
site_rows: list[dict] = []
for family in families:
    name = family["representative_mirna"]
    if name not in wanted:
        continue
    family_rows: list[dict] = []
    for site_type, motif_key in [
        ("8mer", "site_8mer"),
        ("7mer-m8", "site_7mer_m8"),
        ("7mer-a1", "site_7mer_a1"),
    ]:
        motif = family[motif_key]
        starts: list[int] = []
        if motif:
            start = sequence.find(motif)
            while start != -1:
                starts.append(start)
                start = sequence.find(motif, start + 1)
        for pos in starts:
            family_rows.append(
                {
                    "representative_mirna": name,
                    "family_id": family["family_id"],
                    "site_type": site_type,
                    "motif": motif,
                    "dur_start": pos,
                    "dur_end": pos + len(motif),
                    "plot_start": short_len + pos,
                    "plot_end": short_len + pos + len(motif),
                }
            )
    family_rows = sorted(
        family_rows,
        key=lambda item: (
            {"8mer": 0, "7mer-m8": 1, "7mer-a1": 2}[item["site_type"]],
            item["dur_start"],
        ),
    )[:MAX_SITES_PER_FAMILY]
    site_rows.extend(family_rows)
site_df = pd.DataFrame(site_rows)

feature_df = site_df.copy()
for col in [
    "candidate_id",
    "gene_name",
    "stage",
    "rank",
    "mirna_total_sites",
    "mirna_sites_per_100nt",
    "effect_direction_status",
]:
    feature_df[col] = row[col] if col in row.index else candidate_id
feature_df.to_csv(f"{TABLE_PREFIX}.plotted_sites.tsv", sep="\t", index=False)


# --- figure ---

dur_len = int(row["dur_seq_length"])
features: list[GraphicFeature] = [
    GraphicFeature(
        start=0,
        end=short_len,
        strand=0,
        label="shared short 3'UTR",
        color="#9aa4b2",
        thickness=16,
    ),
    GraphicFeature(
        start=short_len,
        end=short_len + dur_len,
        strand=0,
        label="distal UTR region (DUR)",
        color="#55a868",
        thickness=18,
    ),
    GraphicFeature(
        start=max(0, short_len - 1),
        end=short_len + 1,
        strand=0,
        label="proximal PAS / APA boundary",
        color="#111827",
        thickness=28,
    ),
]
if not site_df.empty:
    rank_by_family = {name: idx for idx, name in enumerate(site_df["representative_mirna"].drop_duplicates())}
    labeled_family: set[str] = set()
    for site in site_df.itertuples(index=False):
        style = SITE_STYLE[site.site_type]
        label = None
        if site.representative_mirna not in labeled_family:
            label = site.representative_mirna
            labeled_family.add(site.representative_mirna)
        features.append(
            GraphicFeature(
                start=int(site.plot_start),
                end=int(site.plot_end),
                strand=0,
                label=label,
                color=style["color"],
                thickness=style["height"],
                linecolor="#2d3748",
                feature_level=rank_by_family[site.representative_mirna] % 4,
                legend_text=site.site_type,
            )
        )
sequence_length = int(row["short_utr_length"]) + int(row["dur_seq_length"])
record = GraphicRecord(sequence_length=sequence_length, features=features)

fig = plt.figure(figsize=(12, 4.9))
gs = fig.add_gridspec(2, 1, height_ratios=[3.3, 1.1], hspace=0.20)
ax_track = fig.add_subplot(gs[0])
record.plot(
    ax=ax_track,
    figure_width=12,
    figure_height=3.4,
    with_ruler=True,
    annotate_inline=False,
    max_label_length=24,
)
ax_track.set_title("APA-linked miRNA site exposure", fontsize=13, pad=12)

ax_summary = fig.add_subplot(gs[1])
ax_summary.axis("off")
site_counts = Counter(site_df["site_type"]) if not site_df.empty else Counter()
summary = (
    f"{row['gene_name']} ({row['stage']}) | rank {int(row['rank'])} | "
    f"PIP {float(row['pip']):.3g} | DUR {int(row['dur_seq_length'])} nt | "
    f"miRNA density {float(row['mirna_sites_per_100nt']):.1f}/100 nt "
    f"(threshold {float(row['mirna_density_threshold_per_100nt']):.1f})"
)
direction_note = (
    "Current APA direction is unresolved; this is an evidence map of miRNA site exposure, "
    "not a validated causal direction model."
)
counts = " | ".join(f"{k}: {site_counts.get(k, 0)} plotted" for k in ["8mer", "7mer-m8", "7mer-a1"])
ax_summary.text(0.01, 0.78, summary, fontsize=9.5, ha="left", va="top", transform=ax_summary.transAxes)
ax_summary.text(0.01, 0.42, counts, fontsize=8.8, ha="left", va="top", transform=ax_summary.transAxes)
ax_summary.text(0.01, 0.10, direction_note, fontsize=8.6, color="#4b5563", ha="left", va="bottom", transform=ax_summary.transAxes)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
tight_bbox = fig.get_tightbbox(renderer)
layout_df = pd.DataFrame(
    [
        {
            "figure_id": "mirna_mechanism_track",
            "width_in": round(float(fig.get_size_inches()[0]), 3),
            "height_in": round(float(fig.get_size_inches()[1]), 3),
            "n_axes": len(fig.axes),
            "n_plotted_sites": len(site_df),
            "tight_bbox_width_in": round(float(tight_bbox.width), 3),
            "tight_bbox_height_in": round(float(tight_bbox.height), 3),
            "layout_within_canvas": bool(
                tight_bbox.width <= fig.get_size_inches()[0] * 1.25
                and tight_bbox.height <= fig.get_size_inches()[1] * 1.25
            ),
            "outputs": f"{OUT_PREFIX.name}.png;{OUT_PREFIX.name}.pdf",
        }
    ]
)
layout_df.to_csv(LAYOUT_VALIDATION_PATH, sep="\t", index=False)

for ext in ["png", "pdf"]:
    fig.savefig(f"{OUT_PREFIX}.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
