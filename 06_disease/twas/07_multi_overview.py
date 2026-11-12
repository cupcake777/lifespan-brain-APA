# Plot standalone panels of stage-wise APA-TWAS signals.
# Writes the input audit, the shared-feature table, three panels, and a layout check.
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd


# --- config ---
APA_ROOT = Path(os.environ["APA_ROOT"])
DATA_DIR = APA_ROOT / "data" / "06_disease"
FIGURE_DIR = APA_ROOT / "figures" / "06_disease"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
MANIFEST_PATH = DATA_DIR / "gwas_manifest.tsv"
RANKED_ALL_PATH = DATA_DIR / "twas_ranked_all.tsv.gz"
BY_RUN_TABLE_DIR = DATA_DIR / "by_run"
OVERVIEW_AUDIT_PATH = DATA_DIR / "overview_input_audit.tsv"
LAYOUT_CHECK_PATH = DATA_DIR / "overview_layout_check.tsv"
SHARED_TABLE_PATH = DATA_DIR / "twas_shared_signal_bubble.tsv"
STAGE_ORDER = ("prenatal", "postnatal")
STAGE_COLORS = {"prenatal": "#356E9F", "postnatal": "#C65D35"}
MAX_SHARED_FEATURES = 10

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 12,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 9,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
    }
)


# --- helpers ---
def text_boxes(texts, renderer) -> list[Bbox]:
    return [
        text.get_window_extent(renderer=renderer).expanded(1.02, 1.06)
        for text in texts
        if text.get_visible() and text.get_text()
    ]


def save_figure(fig: plt.Figure, stem: Path) -> list[Path]:
    outputs = [Path(f"{stem}.svg"), Path(f"{stem}.pdf"), Path(f"{stem}.png")]
    fig.savefig(outputs[0], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[2], dpi=400, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return outputs


def validate_layout(
    fig: plt.Figure,
    data_axes: list[plt.Axes],
    legends: list[object] | None = None,
    colorbars: list[object] | None = None,
) -> dict[str, object]:
    legends = legends or []
    colorbars = colorbars or []
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tick_boxes = text_boxes(
        [text for ax in data_axes for text in ax.get_xticklabels() + ax.get_yticklabels()], renderer
    )
    figure_text_boxes = text_boxes(fig.texts, renderer)
    canvas = fig.bbox
    out_of_canvas = sum(
        not (canvas.contains(box.x0, box.y0) and canvas.contains(box.x1, box.y1))
        for box in figure_text_boxes
    )
    axes_boxes = [ax.get_window_extent(renderer=renderer) for ax in data_axes]
    legend_boxes = [legend.get_window_extent(renderer=renderer) for legend in legends]
    colorbar_boxes = [colorbar.ax.get_window_extent(renderer=renderer) for colorbar in colorbars]
    legend_overlap = any(legend_box.overlaps(axis_box) for legend_box in legend_boxes for axis_box in axes_boxes)
    colorbar_overlap = any(colorbar_box.overlaps(axis_box) for colorbar_box in colorbar_boxes for axis_box in axes_boxes)
    colorbar_legend_overlap = any(
        colorbar_box.overlaps(legend_box) for colorbar_box in colorbar_boxes for legend_box in legend_boxes
    )
    tick_overlap = sum(left.overlaps(right) for i, left in enumerate(tick_boxes) for right in tick_boxes[i + 1 :])
    return {
        "tick_label_overlap_count": tick_overlap,
        "out_of_canvas_text_count": out_of_canvas,
        "legend_overlaps_data_axes": legend_overlap,
        "colorbar_overlaps_data_axes": colorbar_overlap,
        "colorbar_overlaps_size_legend": colorbar_legend_overlap,
        "layout_status": (
            "PASS"
            if tick_overlap + out_of_canvas == 0
            and not legend_overlap
            and not colorbar_overlap
            and not colorbar_legend_overlap
            else "REVIEW"
        ),
    }


# --- load input contract ---
manifest = pd.read_csv(MANIFEST_PATH, sep="\t")
usecols = ["gwas_id", "stage", "zscore", "bonferroni_significant"]
ranked = pd.read_csv(RANKED_ALL_PATH, sep="\t", usecols=usecols)
ranked["source_table"] = str(RANKED_ALL_PATH.relative_to(APA_ROOT))
present_pairs = set(zip(ranked["gwas_id"].astype(str), ranked["stage"].astype(str)))

additions = []
for gwas_id in manifest["gwas_id"].astype(str):
    for stage in STAGE_ORDER:
        if (gwas_id, stage) in present_pairs:
            continue
        path = BY_RUN_TABLE_DIR / stage / f"{gwas_id}.association_ranked.tsv"
        if path.exists():
            by_run = pd.read_csv(path, sep="\t", usecols=usecols)
            by_run["source_table"] = str(path.relative_to(APA_ROOT))
            additions.append(by_run)
if additions:
    ranked = pd.concat([ranked, *additions], ignore_index=True)

ranked["gwas_id"] = ranked["gwas_id"].astype(str)
ranked["stage"] = ranked["stage"].astype(str)
ranked["zscore"] = pd.to_numeric(ranked["zscore"], errors="coerce")
ranked["is_significant"] = ranked["bonferroni_significant"].astype(str).str.lower().eq("true")
ranked["positive_hit"] = ranked["is_significant"] & ranked["zscore"].gt(0)
ranked["negative_hit"] = ranked["is_significant"] & ranked["zscore"].lt(0)
rows = (
    ranked.groupby(["gwas_id", "stage"], as_index=False)
    .agg(
        n_tested_features=("gwas_id", "size"),
        n_bonferroni_hits=("is_significant", "sum"),
        n_positive_z_hits=("positive_hit", "sum"),
        n_negative_z_hits=("negative_hit", "sum"),
        source_table=("source_table", lambda values: ";".join(sorted(set(values)))),
    )
    .merge(manifest[["gwas_id", "trait_label", "ancestry"]], on="gwas_id", how="left")
)

# --- input audit ---
records = []
for item in manifest.sort_values("gwas_id").itertuples(index=False):
    record: dict[str, object] = {
        "gwas_id": str(item.gwas_id),
        "trait_label": getattr(item, "trait_label", ""),
        "ancestry": getattr(item, "ancestry", ""),
    }
    for stage in STAGE_ORDER:
        source = rows.loc[
            rows["gwas_id"].astype(str).eq(str(item.gwas_id)) & rows["stage"].astype(str).eq(stage),
            "source_table",
        ]
        record[f"{stage}_complete"] = not source.empty
        record[f"{stage}_source"] = "" if source.empty else str(source.iloc[0])
    record["completed_trait"] = all(record[f"{stage}_complete"] for stage in STAGE_ORDER)
    records.append(record)
audit = pd.DataFrame(records)
audit.to_csv(OVERVIEW_AUDIT_PATH, sep="\t", index=False)

# --- shared significant features ---
wanted_pairs = set(zip(rows["gwas_id"].astype(str), rows["stage"].astype(str)))
usecols = ["gwas_id", "stage", "gene_name", "stable_feature_id", "zscore", "pvalue", "bonferroni_significant"]
chunks = []
consolidated_pairs: set[tuple[str, str]] = set()
for chunk in pd.read_csv(RANKED_ALL_PATH, sep="\t", usecols=usecols, chunksize=50_000):
    chunk["gwas_id"] = chunk["gwas_id"].astype(str)
    chunk["stage"] = chunk["stage"].astype(str)
    consolidated_pairs.update(zip(chunk["gwas_id"], chunk["stage"]))
    keep_pair = chunk.apply(lambda row: (row["gwas_id"], row["stage"]) in wanted_pairs, axis=1)
    keep_sig = chunk["bonferroni_significant"].astype(str).str.lower().eq("true")
    selected = chunk.loc[keep_pair & keep_sig].copy()
    if not selected.empty:
        selected["source_table"] = str(RANKED_ALL_PATH.relative_to(APA_ROOT))
        chunks.append(selected)

for gwas_id, stage in sorted(wanted_pairs - consolidated_pairs):
    path = BY_RUN_TABLE_DIR / stage / f"{gwas_id}.association_ranked.tsv"
    if path.exists():
        by_run = pd.read_csv(path, sep="\t", usecols=usecols)
        by_run = by_run.loc[by_run["bonferroni_significant"].astype(str).str.lower().eq("true")].copy()
        if not by_run.empty:
            by_run["source_table"] = str(path.relative_to(APA_ROOT))
            chunks.append(by_run)

data = pd.concat(chunks, ignore_index=True)
data["pvalue"] = pd.to_numeric(data["pvalue"], errors="coerce")
data["zscore"] = pd.to_numeric(data["zscore"], errors="coerce")
data = data.loc[data["pvalue"].between(np.nextafter(0, 1), 1) & data["zscore"].notna()].copy()
data["minus_log10_p"] = -np.log10(data["pvalue"].clip(lower=np.nextafter(0, 1)))
data["feature_label"] = data["gene_name"].fillna("").astype(str)
missing_gene = data["feature_label"].isin(["", "nan", "None"])
data.loc[missing_gene, "feature_label"] = (
    data.loc[missing_gene, "stable_feature_id"].astype(str).str.split("|").str[0]
)
score = (
    data.groupby("feature_label")
    .agg(n_gwas=("gwas_id", "nunique"), n_cells=("gwas_id", "size"), best_p=("pvalue", "min"))
    .reset_index()
    .sort_values(["n_gwas", "n_cells", "best_p", "feature_label"], ascending=[False, False, True, True])
)
selected_features = score.head(MAX_SHARED_FEATURES)["feature_label"]
bubble = (
    data.loc[data["feature_label"].isin(selected_features)]
    .sort_values(["stage", "gwas_id", "feature_label", "pvalue", "stable_feature_id"])
    .drop_duplicates(["stage", "gwas_id", "feature_label"], keep="first")
    .copy()
)
bubble.to_csv(SHARED_TABLE_PATH, sep="\t", index=False)

# --- panel A: stage signal counts ---
outputs: list[Path] = []
layout_rows: list[dict[str, object]] = []

stage_counts = rows.pivot_table(
    index="gwas_id", columns="stage", values="n_bonferroni_hits", aggfunc="first", fill_value=0
)
stage_counts = stage_counts.reindex(columns=STAGE_ORDER, fill_value=0)
totals = stage_counts.sum(axis=1)
trait_order = totals.sort_values(ascending=False, kind="stable").index.tolist()
y = np.arange(len(trait_order))

fig, ax_counts = plt.subplots(figsize=(6.3, 6.2))
prenatal = stage_counts.loc[trait_order, "prenatal"].to_numpy(dtype=float)
postnatal = stage_counts.loc[trait_order, "postnatal"].to_numpy(dtype=float)
y_prenatal = y - 0.12
y_postnatal = y + 0.12
for index in range(len(trait_order)):
    ax_counts.plot(
        [prenatal[index], postnatal[index]], [index, index],
        color="#B7BDC3", linewidth=1.0, zorder=1,
    )
ax_counts.scatter(
    prenatal, y_prenatal, s=46, color=STAGE_COLORS["prenatal"], edgecolor="white", linewidth=0.55,
    label="prenatal", zorder=3,
)
ax_counts.scatter(
    postnatal, y_postnatal, s=46, color=STAGE_COLORS["postnatal"], edgecolor="white", linewidth=0.55,
    label="postnatal", zorder=3,
)
max_count = max(float(np.max([prenatal.max(), postnatal.max()])), 1.0)
count_ticks = [value for value in [0, 1, 2, 5, 10, 20, 50, 100] if value <= max_count * 1.25]
ax_counts.set_xscale("symlog", linthresh=1, linscale=0.8)
ax_counts.set_xticks(count_ticks, labels=[str(value) for value in count_ticks])
ax_counts.set_xlim(-0.1, max_count * 1.25)
ax_counts.set_yticks(y)
ax_counts.set_yticklabels(trait_order)
ax_counts.set_ylim(len(trait_order) - 0.45, -0.55)
ax_counts.set_xlabel("Bonferroni-significant TWAS features")
ax_counts.set_title("Stage signal counts", loc="center", fontweight="bold")
ax_counts.grid(axis="x", color="#DCE1E5", linewidth=0.5)
ax_counts.set_axisbelow(True)
count_legend = ax_counts.legend(
    loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
    ncol=1, handletextpad=0.4, labelspacing=0.5,
)
fig.subplots_adjust(left=0.22, right=0.80, top=0.91, bottom=0.12)
layout = validate_layout(fig, [ax_counts], legends=[count_legend])
outputs.extend(save_figure(fig, FIGURE_DIR / "stage_signal_counts"))
layout.update(
    {
        "figure_name": "stage_signal_counts",
        "n_plotted_gwas_traits": len(trait_order),
        "n_shared_features": "",
        "n_active_shared_traits": "",
    }
)
layout_rows.append(layout)

# --- panel B: direction of significant signals ---
direction = rows.pivot_table(
    index="gwas_id", columns="stage", values=["n_positive_z_hits", "n_negative_z_hits"],
    aggfunc="first", fill_value=0,
)
direction = direction.reindex(columns=[
    ("n_positive_z_hits", "prenatal"), ("n_negative_z_hits", "prenatal"),
    ("n_positive_z_hits", "postnatal"), ("n_negative_z_hits", "postnatal"),
], fill_value=0)
totals = direction.sum(axis=1)
trait_order = totals.sort_values(ascending=False, kind="stable").index.tolist()
direction_values = direction.loc[trait_order].to_numpy(dtype=float)
max_log = max(float(np.log1p(direction_values.max())), 1.0)
intensity = np.log1p(direction_values) / max_log
stage_columns = ["prenatal", "prenatal", "postnatal", "postnatal"]
rgba = np.ones((*direction_values.shape, 4), dtype=float)
for column_index, stage in enumerate(stage_columns):
    base = np.asarray(to_rgba(STAGE_COLORS[stage]))
    rgba[:, column_index, :] = (1.0 - intensity[:, column_index, None]) + intensity[:, column_index, None] * base

fig, ax_direction = plt.subplots(figsize=(6.5, 6.2))
ax_direction.imshow(rgba, aspect="auto", interpolation="none")
ax_direction.set_xticks(range(4))
ax_direction.set_xticklabels(["Prenatal\n+", "Prenatal\n-", "Postnatal\n+", "Postnatal\n-"])
ax_direction.set_yticks(range(len(trait_order)))
ax_direction.set_yticklabels(trait_order)
ax_direction.tick_params(axis="x", length=0)
ax_direction.set_title("Direction of significant TWAS signals", loc="center", fontweight="bold")
for spine in ax_direction.spines.values():
    spine.set_visible(False)
for position in np.arange(-0.5, len(trait_order), 1):
    ax_direction.axhline(position, color="white", linewidth=0.7)
ax_direction.axvline(1.5, color="white", linewidth=2.4)
direction_legend = ax_direction.legend(
    handles=[
        Patch(facecolor=STAGE_COLORS["prenatal"], edgecolor="none", label="prenatal"),
        Patch(facecolor=STAGE_COLORS["postnatal"], edgecolor="none", label="postnatal"),
    ],
    title="Stage", loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0,
    frameon=False, handlelength=1.1,
)
fig.subplots_adjust(left=0.22, right=0.80, top=0.91, bottom=0.12)
layout = validate_layout(fig, [ax_direction], legends=[direction_legend])
outputs.extend(save_figure(fig, FIGURE_DIR / "stage_signal_direction"))
layout.update(
    {
        "figure_name": "stage_signal_direction",
        "n_plotted_gwas_traits": len(trait_order),
        "n_shared_features": "",
        "n_active_shared_traits": "",
    }
)
layout_rows.append(layout)

# --- panel C: shared significant features ---
active_traits = sorted(bubble["gwas_id"].astype(str).unique())
feature_order = (
    bubble.groupby("feature_label")
    .agg(n_gwas=("gwas_id", "nunique"), best_p=("pvalue", "min"))
    .reset_index()
    .sort_values(["n_gwas", "best_p", "feature_label"], ascending=[False, True, True])["feature_label"]
    .tolist()
)
trait_to_x = {trait: index for index, trait in enumerate(active_traits)}
feature_to_y = {feature: index for index, feature in enumerate(feature_order)}

fig, bubble_axes = plt.subplots(1, 2, figsize=(10.2, 4.8), sharey=True)
max_abs_z = max(2.5, float(np.nanmax(np.abs(bubble["zscore"]))))
norm = TwoSlopeNorm(vmin=-max_abs_z, vcenter=0, vmax=max_abs_z)
bubble_scatter = None
for index, (ax, stage) in enumerate(zip(bubble_axes, STAGE_ORDER)):
    stage_data = bubble.loc[bubble["stage"].astype(str).eq(stage)].copy()
    x_values = stage_data["gwas_id"].map(trait_to_x).to_numpy(dtype=float)
    y_values = stage_data["feature_label"].map(feature_to_y).to_numpy(dtype=float)
    sizes = 32 + np.clip(stage_data["minus_log10_p"].to_numpy(dtype=float), 0, 20) * 6.5
    bubble_scatter = ax.scatter(
        x_values, y_values, s=sizes, c=stage_data["zscore"], cmap="RdBu_r", norm=norm,
        edgecolor="#303030", linewidth=0.45, alpha=0.95,
    )
    ax.set_xticks(range(len(active_traits)))
    ax.set_xticklabels(active_traits, rotation=65, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(feature_order)))
    if index == 0:
        ax.set_yticklabels(feature_order)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.set_xlim(-0.55, len(active_traits) - 0.45)
    ax.set_ylim(len(feature_order) - 0.45, -0.55)
    ax.tick_params(length=0)
    ax.grid(color="#E0E4E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_title(stage, color=STAGE_COLORS[stage], fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)

colorbar_ax = fig.add_axes([0.925, 0.24, 0.016, 0.50])
bubble_colorbar = fig.colorbar(bubble_scatter, cax=colorbar_ax, orientation="vertical")
bubble_colorbar.set_label("TWAS z-score", labelpad=4)
bubble_colorbar.outline.set_linewidth(0.6)
size_handles = [
    Line2D(
        [], [], marker="o", linestyle="none", markerfacecolor="#A9A9A9", markeredgecolor="#303030",
        markersize=np.sqrt(32 + value * 6.5), label=str(value),
    )
    for value in [5, 10, 20]
]
size_legend = fig.legend(
    handles=size_handles, title="$-\\log_{10}(P)$ (capped at 20)", loc="center left",
    bbox_to_anchor=(0.015, 0.52), ncol=1, borderaxespad=0, handletextpad=0.3,
    columnspacing=0.8, labelspacing=0.45,
)
fig.suptitle("Shared significant features across traits", x=0.5, y=0.94, ha="center", fontsize=12, fontweight="bold")
fig.subplots_adjust(left=0.23, right=0.89, top=0.84, bottom=0.27, wspace=0.16)
layout = validate_layout(
    fig, list(bubble_axes), legends=[size_legend], colorbars=[bubble_colorbar]
)
outputs.extend(save_figure(fig, FIGURE_DIR / "shared_signal_bubble"))
layout.update(
    {
        "figure_name": "shared_signal_bubble",
        "n_plotted_gwas_traits": len(active_traits),
        "n_shared_features": len(feature_order),
        "n_active_shared_traits": len(active_traits),
    }
)
layout_rows.append(layout)

# --- write outputs ---
layout = pd.DataFrame(layout_rows)
layout[
    [
        "figure_name", "n_plotted_gwas_traits", "n_shared_features", "n_active_shared_traits",
        "tick_label_overlap_count", "out_of_canvas_text_count", "legend_overlaps_data_axes",
        "colorbar_overlaps_data_axes", "colorbar_overlaps_size_legend", "layout_status",
    ]
].to_csv(LAYOUT_CHECK_PATH, sep="\t", index=False)
print(f"Completed GWAS traits: {int(audit['completed_trait'].sum())}")
print(f"Wrote {len(outputs)} standalone overview panel files")
print(f"Wrote layout validation: {LAYOUT_CHECK_PATH}")
