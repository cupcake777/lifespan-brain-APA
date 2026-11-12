# Plot paired developmental TWAS effect figures with representative-variant labels.
# Inputs: ranked TWAS tables, GWAS manifest, stage model databases.
from __future__ import annotations

import csv
import os
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


# --- config ---
APA_ROOT = Path(os.environ["APA_ROOT"])
DATA_DIR = APA_ROOT / "data" / "06_disease"
FIGURE_DIR = APA_ROOT / "figures" / "06_disease"
RANKED_ALL_PATH = DATA_DIR / "twas_ranked_all.tsv.gz"
BY_RUN_TABLE_DIR = DATA_DIR / "by_run"
MANIFEST_PATH = DATA_DIR / "gwas_manifest.tsv"
MODEL_DIR = Path(os.environ["APA_HPC_HOME"]) / "project" / "work" / "apa" / "08_twas" / "models"
STAGES = ("prenatal", "postnatal")
STAGE_COLORS = {"prenatal": "#2F6F9F", "postnatal": "#D06945"}
SIGN_COLORS = {True: "#D55E00", False: "#3B5F8A"}
MAX_LABELS_PER_STAGE = 3
BACKGROUND_QUANTILE = 0.985
CHROM_GAP = 5_000_000.0
CHROMOSOME_COLORS = ("#0072B2", "#56B4E9")
POINT_SIZE_MIN = 4.0
POINT_SIZE_MAX = 24.0
PVALUE_SIZE_CAP = 25.0
TRIANGLE_SIZE_ADD = 3.0

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)]
CHROM_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468,
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
    }
)


# --- helpers ---
def bbox_overlap_count(boxes: list[Bbox]) -> int:
    return sum(left.overlaps(right) for i, left in enumerate(boxes) for right in boxes[i + 1 :])


def pvalue_marker_sizes(values: pd.Series | np.ndarray) -> np.ndarray:
    minus_log10_p = np.asarray(values, dtype=float)
    scaled = np.clip(minus_log10_p, 0.0, PVALUE_SIZE_CAP) / PVALUE_SIZE_CAP
    return POINT_SIZE_MIN + scaled * (POINT_SIZE_MAX - POINT_SIZE_MIN)


# --- chromosome layout ---
offsets, centers, ends = {}, {}, {}
current = 0.0
for chrom in CHROM_ORDER:
    offsets[chrom] = current
    centers[chrom] = current + CHROM_LENGTHS[chrom] / 2.0
    ends[chrom] = current + CHROM_LENGTHS[chrom]
    current = ends[chrom] + CHROM_GAP

# --- load ranked TWAS rows ---
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
gwas_ids = ["EUR_SCZ"]
usecols = [
    "gwas_id", "stage", "model_id", "gene_name", "stable_feature_id", "3utr_region",
    "chr", "group_left", "group_right", "effect_size", "pvalue",
    "zscore", "bonferroni_significant",
]
wanted = set(gwas_ids)
chunks = []
for chunk in pd.read_csv(RANKED_ALL_PATH, sep="\t", usecols=usecols, chunksize=50_000):
    selected = chunk.loc[chunk["gwas_id"].astype(str).isin(wanted)].copy()
    if not selected.empty:
        selected["source_table"] = str(RANKED_ALL_PATH.relative_to(APA_ROOT))
        chunks.append(selected)

present_pairs = {
    (str(gwas_id), str(stage))
    for chunk in chunks
    for gwas_id, stage in zip(chunk["gwas_id"], chunk["stage"])
}
for gwas_id in gwas_ids:
    for stage in STAGES:
        if (gwas_id, stage) in present_pairs:
            continue
        path = BY_RUN_TABLE_DIR / stage / f"{gwas_id}.association_ranked.tsv"
        if path.exists():
            by_run = pd.read_csv(path, sep="\t", usecols=usecols)
            by_run["source_table"] = str(path.relative_to(APA_ROOT))
            chunks.append(by_run)

data = pd.concat(chunks, ignore_index=True)

# --- coordinates ---
coord_chr, coord_mid = [], []
for _, row in data.iterrows():
    chrom = str(row.get("chr", ""))
    if chrom and not chrom.startswith("chr"):
        chrom = f"chr{chrom}"
    start = pd.to_numeric(row.get("group_left"), errors="coerce")
    end = pd.to_numeric(row.get("group_right"), errors="coerce")
    if (pd.isna(start) or pd.isna(end)) and pd.notna(row.get("3utr_region")):
        match = re.search(r"(chr\w+):(\d+)-(\d+)", str(row["3utr_region"]))
        if match:
            chrom = match.group(1)
            start = float(match.group(2))
            end = float(match.group(3))
    if pd.notna(start) and pd.notna(end):
        coord_chr.append(chrom)
        coord_mid.append((float(start) + float(end)) / 2.0)
    else:
        coord_chr.append(None)
        coord_mid.append(np.nan)
data["coord_chr"] = coord_chr
data["coord_mid"] = coord_mid
data = data.loc[data["coord_chr"].isin(CHROM_ORDER) & data["coord_mid"].notna()].copy()
data["plot_x"] = data["coord_chr"].map(offsets) + data["coord_mid"]
data["effect_size"] = pd.to_numeric(data["effect_size"], errors="coerce")
data["pvalue"] = pd.to_numeric(data["pvalue"], errors="coerce")
data = data.loc[data["effect_size"].notna() & data["pvalue"].between(np.nextafter(0, 1), 1)].copy()
data["minus_log10_p"] = -np.log10(data["pvalue"].clip(lower=np.nextafter(0, 1)))
data["is_significant"] = data["bonferroni_significant"].astype(str).str.lower().eq("true")
data["feature_label"] = data["gene_name"].fillna("").astype(str)
missing_gene = data["feature_label"].isin(["", "nan", "None"])
data.loc[missing_gene, "feature_label"] = (
    data.loc[missing_gene, "stable_feature_id"].astype(str).str.split("|").str[0]
)

# --- per-GWAS figure ---
output_count = 0
for gwas_id in gwas_ids:
    meta = manifest.loc[manifest["gwas_id"].astype(str).eq(gwas_id)].iloc[0]
    trait_label = str(meta["trait_label"])
    ancestry = str(meta["ancestry"])
    gwas_source = Path(str(meta["gwas_file"]))

    gwas_data = data.loc[data["gwas_id"].astype(str).eq(gwas_id)].copy()
    significant = gwas_data.loc[gwas_data["is_significant"]].copy()
    if significant.empty:
        continue

    background_abs = gwas_data.loc[~gwas_data["is_significant"], "effect_size"].abs()
    background_limit = float(background_abs.quantile(BACKGROUND_QUANTILE)) if not background_abs.empty else 0.0
    significant_limit = float(significant["effect_size"].abs().max())
    display_limit = max(0.6, background_limit * 1.20, significant_limit * 1.18)

    chromosome_colors = {
        chrom: CHROMOSOME_COLORS[index % 2]
        for index, chrom in enumerate(CHROM_ORDER)
    }

    # --- label selection ---
    selected_by_stage = {}
    for stage in STAGES:
        selected_by_stage[stage] = (
            significant.loc[significant["stage"].astype(str).eq(stage)]
            .sort_values(["pvalue", "stable_feature_id"])
            .drop_duplicates("feature_label")
            .head(MAX_LABELS_PER_STAGE)
            .copy()
        )
    selected_all = pd.concat(selected_by_stage.values(), ignore_index=True)

    # --- representative GWAS variants ---
    labeled_hits = selected_all
    targets: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    twas_lookup = labeled_hits.set_index(["stage", "model_id"])["zscore"].to_dict()

    for stage in STAGES:
        model_ids = sorted(
            labeled_hits.loc[labeled_hits["stage"].astype(str).eq(stage), "model_id"]
            .dropna().astype(str).unique()
        )
        if not model_ids:
            continue
        db_path = MODEL_DIR / stage / f"{stage}.full_twas.db"
        placeholders = ",".join("?" for _ in model_ids)
        query = (
            "SELECT rsid, gene AS model_id, weight, ref_allele, eff_allele "
            f"FROM weights WHERE gene IN ({placeholders})"
        )
        with sqlite3.connect(db_path) as connection:
            weights = pd.read_sql_query(query, connection, params=model_ids)
        for row in weights.itertuples(index=False):
            chrom, pos, _, _ = str(row.rsid).split(":")
            ref = str(row.ref_allele).upper()
            eff = str(row.eff_allele).upper()
            targets[(chrom.removeprefix("chr"), pos)].append(
                {
                    "stage": stage,
                    "model_id": str(row.model_id),
                    "rsid": str(row.rsid),
                    "weight": float(row.weight),
                    "ref_allele": ref,
                    "eff_allele": eff,
                    "allele_pair": frozenset((ref, eff)),
                }
            )

    matched: list[dict[str, object]] = []
    with gwas_source.open() as source:
        reader = csv.DictReader(source, delimiter="\t")
        for row in reader:
            candidates = targets.get((row["CHR"].removeprefix("chr"), row["POS"]))
            if not candidates:
                continue
            ea = row["EA"].upper()
            nea = row["NEA"].upper()
            beta = float(row["BETA"])
            se = float(row["SE"])
            pvalue = float(row["P"])
            if not np.isfinite(beta) or not np.isfinite(se) or se <= 0 or not (0 < pvalue <= 1):
                continue
            for candidate in candidates:
                if frozenset((ea, nea)) != candidate["allele_pair"]:
                    continue
                if ea == candidate["eff_allele"] and nea == candidate["ref_allele"]:
                    aligned_beta = beta
                    alignment = "exact"
                elif ea == candidate["ref_allele"] and nea == candidate["eff_allele"]:
                    aligned_beta = -beta
                    alignment = "flipped"
                else:
                    continue
                aligned_z = aligned_beta / se
                signed_term = float(candidate["weight"]) * aligned_z
                twas_z = float(twas_lookup[(candidate["stage"], candidate["model_id"])])
                matched.append(
                    {
                        **{key: value for key, value in candidate.items() if key != "allele_pair"},
                        "gwas_pvalue": pvalue,
                        "gwas_snpid": str(row.get("SNPID", "")),
                        "gwas_variant_id": str(row.get("variant_id", "")),
                        "gwas_beta_aligned": aligned_beta,
                        "gwas_se": se,
                        "gwas_z_aligned": aligned_z,
                        "weight_x_gwas_z": signed_term,
                        "twas_zscore": twas_z,
                        "direction_concordant": bool(np.sign(signed_term) == np.sign(twas_z)),
                        "allele_alignment": alignment,
                    }
                )

    matched_table = pd.DataFrame(matched)
    if matched_table.empty:
        raise ValueError(f"No selected model variants matched {gwas_source}")
    representative_variants = (
        matched_table.sort_values(["stage", "model_id", "gwas_pvalue", "rsid"])
        .drop_duplicates(["stage", "model_id"], keep="first")
        .reset_index(drop=True)
    )
    selected_all = labeled_hits.merge(
        representative_variants,
        on=["stage", "model_id"],
        how="left",
        validate="1:1",
    )
    if selected_all["rsid"].isna().any():
        missing_models = selected_all.loc[selected_all["rsid"].isna(), "model_id"].astype(str).tolist()
        raise ValueError(f"No representative GWAS variant for selected models: {missing_models}")
    selected_all["variant_label"] = selected_all["gwas_snpid"].where(
        selected_all["gwas_snpid"].astype(str).str.strip().ne(""), selected_all["rsid"]
    )
    selected_all["display_label"] = selected_all["feature_label"] + "\n" + selected_all["variant_label"].astype(str)
    selected_by_stage = {
        stage: selected_all.loc[selected_all["stage"].astype(str).eq(stage)].copy()
        for stage in STAGES
    }

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.75), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.105, right=0.89, top=0.875, bottom=0.115, hspace=0.10)
    all_labels: list[plt.Annotation] = []
    clipped_counts: dict[str, int] = {}
    selected_hits: list[pd.DataFrame] = []

    for ax, stage in zip(axes, STAGES):
        stage_data = gwas_data.loc[gwas_data["stage"].astype(str).eq(stage)].copy()
        stage_labels = selected_by_stage[stage]
        selected_hits.append(stage_labels)
        labeled_ids = set(stage_labels["stable_feature_id"].astype(str))
        stage_points = stage_data.loc[~stage_data["stable_feature_id"].astype(str).isin(labeled_ids)].copy()

        within = stage_points["effect_size"].abs().le(display_limit)
        ax.scatter(
            stage_points.loc[within, "plot_x"],
            stage_points.loc[within, "effect_size"],
            s=pvalue_marker_sizes(stage_points.loc[within, "minus_log10_p"]),
            c=stage_points.loc[within, "coord_chr"].map(chromosome_colors),
            alpha=0.48,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
        clipped = stage_points.loc[~within].copy()
        clipped_counts[stage] = len(clipped)
        if not clipped.empty:
            ax.scatter(
                clipped["plot_x"],
                np.sign(clipped["effect_size"]) * display_limit * 0.985,
                s=pvalue_marker_sizes(clipped["minus_log10_p"]) + 4.0,
                marker="_",
                c=clipped["coord_chr"].map(chromosome_colors),
                alpha=0.78,
                linewidths=0.8,
                rasterized=True,
                zorder=3,
            )

        for positive, marker in ((True, "^"), (False, "v")):
            subset = stage_labels.loc[stage_labels["effect_size"].ge(0).eq(positive)]
            if subset.empty:
                continue
            ax.scatter(
                subset["plot_x"], subset["effect_size"],
                s=pvalue_marker_sizes(subset["minus_log10_p"]) + TRIANGLE_SIZE_ADD,
                color=SIGN_COLORS[positive],
                marker=marker,
                edgecolor="#2A2A2A",
                linewidth=0.45,
                alpha=0.98,
                zorder=6,
            )

        ax.axhline(0, color="#59636A", linewidth=0.65, zorder=4)
        ax.set_ylim(-display_limit, display_limit)
        ax.set_xlim(0, max(ends.values()))
        ax.grid(axis="y", color="#E3E7EA", linewidth=0.42)
        ax.set_axisbelow(True)
        ax.text(
            0.012, 0.91, stage,
            transform=ax.transAxes,
            fontsize=8.0,
            fontweight="bold",
            color=STAGE_COLORS[stage],
            ha="left", va="top",
        )

        # --- label placement ---
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        axes_box = ax.get_window_extent(renderer=renderer)
        placed_boxes: list[Bbox] = []
        placed_labels: list[plt.Annotation] = []
        label_offsets = [(0, 22), (18, 19), (-18, 19), (28, -20), (-28, -20), (0, -24)]
        for row in stage_labels.itertuples(index=False):
            sign_color = SIGN_COLORS[float(row.effect_size) >= 0]
            for x_offset, y_offset in label_offsets:
                label = ax.annotate(
                    str(row.display_label),
                    xy=(row.plot_x, row.effect_size),
                    xytext=(x_offset, y_offset),
                    textcoords="offset points",
                    fontsize=6.8,
                    fontweight="bold",
                    ha="center" if x_offset == 0 else ("left" if x_offset > 0 else "right"),
                    va="bottom" if y_offset > 0 else "top",
                    color="#252525",
                    arrowprops={"arrowstyle": "-", "color": sign_color, "linewidth": 0.55, "shrinkB": 4},
                    bbox={"boxstyle": "square,pad=0.10", "facecolor": "white", "edgecolor": "none", "alpha": 0.94},
                    zorder=7,
                )
                fig.canvas.draw()
                box = label.get_window_extent(renderer=renderer).expanded(1.04, 1.12)
                inside = axes_box.contains(box.x0, box.y0) and axes_box.contains(box.x1, box.y1)
                if inside and not any(box.overlaps(other) for other in placed_boxes):
                    placed_labels.append(label)
                    placed_boxes.append(box)
                    break
                label.remove()
        all_labels.extend(placed_labels)

    labeled_chromosomes = CHROM_ORDER[::2]
    axes[-1].set_xticks([centers[chrom] for chrom in labeled_chromosomes])
    axes[-1].set_xticklabels([chrom.replace("chr", "") for chrom in labeled_chromosomes])
    axes[-1].tick_params(axis="x", length=0, pad=3)
    axes[-1].set_xlabel("Chromosome")
    fig.supylabel("SPrediXcan effect estimate", x=0.025, fontsize=8.2)
    title = f"{trait_label} | {ancestry}" if ancestry else trait_label
    fig.suptitle(title, x=0.105, y=0.955, ha="left", fontsize=9.5, fontweight="bold")

    legend = fig.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", markerfacecolor="#78A9C5", markeredgecolor="none", markersize=3.1, label="APA feature"),
            Line2D([], [], marker="^", linestyle="none", markerfacecolor=SIGN_COLORS[True], markeredgecolor="#2A2A2A", markeredgewidth=0.45, markersize=3.8, label="positive"),
            Line2D([], [], marker="v", linestyle="none", markerfacecolor=SIGN_COLORS[False], markeredgecolor="#2A2A2A", markeredgewidth=0.45, markersize=3.8, label="negative"),
        ],
        loc="upper right", bbox_to_anchor=(0.89, 0.955), ncol=3,
        handletextpad=0.30, columnspacing=0.75,
    )
    size_values = np.array([5.0, 15.0, 25.0])
    size_legend = fig.legend(
        handles=[
            Line2D(
                [], [], marker="o", linestyle="none",
                markerfacecolor="#56B4E9", markeredgecolor="none",
                markersize=float(np.sqrt(size)), label=f"{value:.0f}",
            )
            for value, size in zip(size_values, pvalue_marker_sizes(size_values))
        ],
        title=r"$-\log_{10}(P)$",
        loc="center right", bbox_to_anchor=(0.985, 0.49), ncol=1,
        handletextpad=0.35, labelspacing=0.55,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_boxes = [label.get_window_extent(renderer=renderer).expanded(1.03, 1.08) for label in all_labels]
    axes_boxes = [ax.get_window_extent(renderer=renderer) for ax in axes]
    out_of_axes = sum(
        not any(box.x0 >= ax_box.x0 and box.x1 <= ax_box.x1 and box.y0 >= ax_box.y0 and box.y1 <= ax_box.y1 for ax_box in axes_boxes)
        for box in label_boxes
    )
    label_overlap = bbox_overlap_count(label_boxes)
    legend_overlap = any(legend.get_window_extent(renderer=renderer).overlaps(ax_box) for ax_box in axes_boxes)
    size_legend_overlap = any(size_legend.get_window_extent(renderer=renderer).overlaps(ax_box) for ax_box in axes_boxes)
    tick_boxes = [
        tick.get_window_extent(renderer=renderer).expanded(1.04, 1.08)
        for tick in axes[-1].get_xticklabels()
        if tick.get_visible() and tick.get_text()
    ]
    tick_overlap = bbox_overlap_count(tick_boxes)

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    stem_name = f"{gwas_id}.twas_effect_pair"
    stem = FIGURE_DIR / stem_name
    outputs = [Path(f"{stem}.svg"), Path(f"{stem}.pdf"), Path(f"{stem}.png")]
    fig.savefig(outputs[0], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[2], dpi=400, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)

    # --- QC tables ---
    labeled_hits = pd.concat(selected_hits, ignore_index=True)
    labeled_hits.to_csv(
        DATA_DIR / f"{stem_name}.labeled_hits.tsv",
        sep="\t", index=False,
        columns=[
            "gwas_id", "stage", "feature_label", "stable_feature_id", "model_id",
            "effect_size", "zscore", "pvalue", "minus_log10_p", "rsid", "weight",
            "eff_allele", "ref_allele", "variant_label", "gwas_snpid", "gwas_variant_id",
            "gwas_beta_aligned", "gwas_se", "gwas_pvalue",
            "gwas_z_aligned", "weight_x_gwas_z", "direction_concordant",
            "allele_alignment", "source_table",
        ],
    )
    status = "PASS" if label_overlap + out_of_axes + tick_overlap == 0 and not legend_overlap and not size_legend_overlap else "REVIEW"
    pd.DataFrame(
        [{
            "figure_name": stem_name,
            "gwas_id": gwas_id,
            "n_plotted_features": len(gwas_data),
            "n_significant_features": len(significant),
            "pvalue_size_cap": PVALUE_SIZE_CAP,
            "point_size_min": POINT_SIZE_MIN,
            "point_size_max": POINT_SIZE_MAX,
            "triangle_size_add": TRIANGLE_SIZE_ADD,
            "display_limit": display_limit,
            "prenatal_background_clipped": clipped_counts.get("prenatal", 0),
            "postnatal_background_clipped": clipped_counts.get("postnatal", 0),
            "n_placed_labels": len(all_labels),
            "n_representative_variants": len(representative_variants),
            "n_direction_concordant": int(representative_variants["direction_concordant"].sum()),
            "n_direction_discordant": int((~representative_variants["direction_concordant"]).sum()),
            "point_label_overlap_count": label_overlap,
            "out_of_axes_text_count": out_of_axes,
            "tick_label_overlap_count": tick_overlap,
            "legend_overlaps_data_axes": legend_overlap,
            "size_legend_overlaps_data_axes": size_legend_overlap,
            "layout_status": status,
        }]
    ).to_csv(DATA_DIR / f"{stem_name}.layout_check.tsv", sep="\t", index=False)
    output_count += len(outputs)

print(f"Wrote {output_count} effect figure files")
