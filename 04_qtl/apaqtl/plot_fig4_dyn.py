# Figure 4 dynamic APA-QTL panel (3'UTR model, PDUI boxes, WIG group tracks).
# Source tables and layout checks go to data/04_qtl, plots to figures/04_qtl.
import os
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd


# --- helpers ---

def run_text(command):
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout


def add_panel_label(axis, index):
    axis.text(
        -0.11,
        1.06,
        chr(ord("a") + index),
        transform=axis.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
    )


# --- config ---

APA_ROOT = Path(os.environ["APA_ROOT"])
data_dir = APA_ROOT / "data" / "04_qtl"
figure_dir = APA_ROOT / "figures" / "04_qtl"
vcf_file = Path(f"{os.environ['APA_HPC_HOME']}/project/raw_data/CHB_all_anno.vcf.gz")
screen_file = data_dir / "figure4_dynamic_input_screen.tsv"
phenotype_files = {
    "Prenatal": data_dir / "pre_input/phenotype.Prenatal.bed.gz",
    "Postnatal": data_dir / "pre_input/phenotype.Postnatal.bed.gz",
}
wig_profile_file = data_dir / "figure4_dynamic_input_wig_profiles.tsv.gz"
wig_qc_file = data_dir / "figure4_dynamic_input_wig_profile_qc.tsv"
wig_selection_file = data_dir / "figure4_dynamic_input_wig_profile_samples.tsv"

candidate_keys = [
    {"gene": "PACC1", "variant_id": "chr1:212366000:C:T"},
    {"gene": "SLCO1A2", "variant_id": "chr12:21275437:T:C"},
]

png_file = figure_dir / "png/figure4_dynamic.png"
pdf_file = figure_dir / "pdf/figure4_dynamic.pdf"
svg_file = figure_dir / "svg/figure4_dynamic.svg"
phenotype_source_file = data_dir / "figure4_dynamic_phenotype_source.tsv"
coverage_source_file = data_dir / "figure4_dynamic_coverage_source.tsv"
sample_qc_file = data_dir / "figure4_dynamic_coverage_sample_qc.tsv"
summary_file = data_dir / "figure4_dynamic_summary.tsv"
layout_file = data_dir / "figure4_dynamic_layout_validation.tsv"

stage_order = ["Prenatal", "Postnatal"]
genotype_order = ["0|0", "0|1", "1|1"]
stage_colors = {"Prenatal": "#D55E00", "Postnatal": "#0072B2"}
genotype_colors = {"0|0": "#2A9D8F", "0|1": "#E69F00", "1|1": "#3B6FB6"}
neutral = "#666666"
dpi = 360


# --- style ---

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"
mpl.rcParams.update(
    {
        "font.size": 10.2,
        "axes.labelsize": 10.6,
        "axes.titlesize": 11.4,
        "xtick.labelsize": 9.6,
        "ytick.labelsize": 9.6,
        "axes.linewidth": 0.75,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,
    }
)


# --- load audited candidate rows ---

screen_columns = [
    "target_stage",
    "gene",
    "phenotype_id",
    "variant_id",
    "target_p",
    "target_beta",
    "target_q",
    "other_p",
    "other_beta",
    "other_q",
    "p_effect_difference_fdr",
    "region",
    "strand",
    "proximal_pas",
    "common_utr_bp",
    "long_extension_bp",
    "annotation_span_bp",
    "annotation_distal_extension_bp",
    "short_utr_bp",
    "distal_pas",
    "coverage_status",
    "target_coverage_endpoint_delta_bp",
    "target_00_coverage_n",
    "target_00_endpoint_fraction_median",
    "target_01_coverage_n",
    "target_01_endpoint_fraction_median",
    "target_11_coverage_n",
    "target_11_endpoint_fraction_median",
    "other_coverage_endpoint_delta_bp",
    "other_00_coverage_n",
    "other_00_endpoint_fraction_median",
    "other_01_coverage_n",
    "other_01_endpoint_fraction_median",
    "other_11_coverage_n",
    "other_11_endpoint_fraction_median",
]
screen = pd.read_csv(screen_file, sep="\t", usecols=screen_columns)

candidate_rows = []
for key in candidate_keys:
    selected = screen.loc[
        (screen["gene"] == key["gene"])
        & (screen["variant_id"] == key["variant_id"])
    ]
    if len(selected) != 1:
        raise ValueError(f"Expected one screen row for {key}, found {len(selected)}")
    candidate = selected.iloc[0].to_dict()
    candidate["other_stage"] = next(
        stage for stage in stage_order if stage != candidate["target_stage"]
    )
    if int(candidate["common_utr_bp"]) != int(candidate["short_utr_bp"]):
        raise ValueError(f"Short 3'UTR mismatch for {candidate['gene']}")
    if int(candidate["annotation_span_bp"]) != int(candidate["short_utr_bp"] + candidate["long_extension_bp"]):
        raise ValueError(f"Long 3'UTR mismatch for {candidate['gene']}")
    if int(candidate["annotation_distal_extension_bp"]) != int(candidate["long_extension_bp"]):
        raise ValueError(f"Distal extension mismatch for {candidate['gene']}")
    candidate_rows.append(candidate)


# --- load formal APA phenotype values and genotype groups ---

phenotype_rows = []
genotype_cache = {}
for candidate in candidate_rows:
    variant_id = candidate["variant_id"]
    phenotype_id = candidate["phenotype_id"]
    if variant_id not in genotype_cache:
        chrom, position, ref, alt = variant_id.split(":")
        contig = chrom.removeprefix("chr")
        output = run_text(
            [
                "bcftools",
                "query",
                "-r",
                f"{contig}:{position}-{position}",
                "-f",
                "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE=%GT]\\n",
                str(vcf_file),
            ]
        )
        matching_fields = None
        for line in output.rstrip().splitlines():
            fields = line.split("\t")
            if fields[:4] == [contig, position, ref, alt]:
                matching_fields = fields
                break
        if matching_fields is None:
            raise ValueError(f"Variant not found in VCF: {variant_id}")

        rows = []
        dosage_map = {"0|0": 0, "0|1": 1, "1|1": 2}
        for sample_field in matching_fields[4:]:
            sample, genotype = sample_field.rsplit("=", 1)
            alleles = genotype.replace("/", "|").split("|")
            if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
                continue
            genotype = "|".join(sorted(alleles))
            rows.append(
                {
                    "sample": sample,
                    "genotype": genotype,
                    "dosage": dosage_map[genotype],
                }
            )
        genotype_cache[variant_id] = pd.DataFrame(rows)

    for stage in stage_order:
        _, _, region, _ = phenotype_id.split("|")
        chrom, coordinates = region.split(":")
        lines = run_text(
            [
                "tabix",
                "-h",
                str(phenotype_files[stage]),
                f"{chrom.removeprefix('chr')}:{coordinates}",
            ]
        ).rstrip().splitlines()
        header = next(line for line in lines if line.lower().startswith("#chr")).split("\t")
        target = next(
            line.split("\t")
            for line in lines
            if not line.startswith("#") and line.split("\t", 4)[3] == phenotype_id
        )
        phenotype = pd.DataFrame(
            {
                "sample": header[4:],
                "apa_score": pd.to_numeric(target[4:], errors="coerce"),
            }
        )
        merged = phenotype.merge(genotype_cache[variant_id], on="sample", how="inner")
        merged = merged.dropna(subset=["apa_score"]).copy()
        merged["gene"] = candidate["gene"]
        merged["variant_id"] = variant_id
        merged["phenotype_id"] = phenotype_id
        merged["stage"] = stage
        phenotype_rows.append(merged)

phenotype_source = pd.concat(phenotype_rows, ignore_index=True)
phenotype_source.to_csv(phenotype_source_file, sep="\t", index=False)


# --- load representative full-3'UTR WIG coverage ---

wig_profiles = pd.read_csv(wig_profile_file, sep="\t")
wig_selection = pd.read_csv(wig_selection_file, sep="\t")
coverage_source = wig_profiles.merge(
    wig_selection,
    on=["gene", "phenotype_id", "variant_id", "stage", "sample", "genotype"],
    how="inner",
)
coverage_source["position_along_3utr_bp"] = (
    coverage_source["bin_start_bp"] + coverage_source["bin_end_bp"]
) / 2
coverage_source["source_format"] = "WIG"
coverage_source = coverage_source.sort_values(
    ["gene", "stage", "genotype", "position_along_3utr_bp"]
)
coverage_source.to_csv(coverage_source_file, sep="\t", index=False)

wig_qc = pd.read_csv(wig_qc_file, sep="\t")
sample_qc = wig_qc.merge(
    wig_selection[
        ["gene", "phenotype_id", "variant_id", "stage", "sample", "genotype"]
    ],
    on=["gene", "phenotype_id", "variant_id", "stage", "sample", "genotype"],
    how="inner",
)
sample_qc = sample_qc.sort_values(["gene", "stage", "genotype", "sample"])
sample_qc.to_csv(sample_qc_file, sep="\t", index=False)


# --- write candidate summary ---

summary_rows = []
for candidate in candidate_rows:
    gene_selection = wig_selection.loc[
        wig_selection["gene"] == candidate["gene"]
    ].sort_values(["stage", "genotype"])
    sample_summary = ";".join(
        f"{row.stage}:{row.genotype}={row.sample}"
        for row in gene_selection.itertuples(index=False)
    )
    summary_rows.append(
        {
            "gene": candidate["gene"],
            "variant_id": candidate["variant_id"],
            "phenotype_id": candidate["phenotype_id"],
            "target_stage": candidate["target_stage"],
            "comparator_stage": candidate["other_stage"],
            "short_utr_bp": int(candidate["short_utr_bp"]),
            "long_utr_bp": int(candidate["short_utr_bp"] + candidate["long_extension_bp"]),
            "distal_extension_bp": int(candidate["long_extension_bp"]),
            "target_beta": float(candidate["target_beta"]),
            "target_q": float(candidate["target_q"]),
            "comparator_beta": float(candidate["other_beta"]),
            "comparator_q": float(candidate["other_q"]),
            "stage_effect_difference_fdr": float(candidate["p_effect_difference_fdr"]),
            "target_wig_endpoint_delta_bp": float(candidate["target_coverage_endpoint_delta_bp"]),
            "comparator_wig_endpoint_delta_bp": float(candidate["other_coverage_endpoint_delta_bp"]),
            "representative_wig_samples": sample_summary,
            "representative_track_n": int(len(gene_selection)),
            "source_screen_wig_n": int(gene_selection["group_valid_wig_n"].sum()),
            "coverage_normalization": "per-sample mean coverage across the common 3'UTR",
            "coverage_aggregation": "none; one representative WIG per stage-genotype group",
            "representative_selection": (
                "nearest stage-genotype group median endpoint fraction in both stages"
            ),
            "inferential_evidence": "full-cohort APA-QTL distributions and formal q/FDR",
            "evidence_boundary": (
                "Short-read coverage supports distal continuation, not an exact cleavage or PAS coordinate."
            ),
        }
    )

summary = pd.DataFrame(summary_rows)
summary.to_csv(summary_file, sep="\t", index=False)


# --- draw short/long 3'UTR evidence grid ---

png_file.parent.mkdir(parents=True, exist_ok=True)
pdf_file.parent.mkdir(parents=True, exist_ok=True)
svg_file.parent.mkdir(parents=True, exist_ok=True)

fig = plt.figure(figsize=(13.2, 7.4), facecolor="white")
grid = fig.add_gridspec(
    2,
    3,
    left=0.055,
    right=0.985,
    bottom=0.095,
    top=0.955,
    width_ratios=[2.30, 2.55, 3.25],
    hspace=0.64,
    wspace=0.34,
)

global_y_min = float(np.floor(phenotype_source["apa_score"].min() * 2) / 2 - 0.15)
global_y_max = float(np.ceil(phenotype_source["apa_score"].max() * 2) / 2 + 0.15)
layout_rows = []

for row_index, candidate in enumerate(candidate_rows):
    gene = candidate["gene"]
    target_stage = candidate["target_stage"]
    comparator_stage = candidate["other_stage"]
    target_color = stage_colors[target_stage]
    short_bp = int(candidate["short_utr_bp"])
    extension_bp = int(candidate["long_extension_bp"])
    long_bp = short_bp + extension_bp
    gene_phenotype = phenotype_source.loc[phenotype_source["gene"] == gene]
    gene_coverage = coverage_source.loc[coverage_source["gene"] == gene]
    _, _, ref, alt = candidate["variant_id"].split(":")
    allele_map = {"0": ref, "1": alt}
    genotype_labels = {
        genotype: "/".join(allele_map[allele] for allele in genotype.split("|"))
        for genotype in genotype_order
    }

    # Primary model: annotation-defined short/long 3'UTRs and WIG continuation shift.
    model_ax = fig.add_subplot(grid[row_index, 0])
    axis_left = -0.06 * long_bp
    axis_right = 1.06 * long_bp
    model_ax.set_xlim(axis_left, axis_right)
    model_ax.set_ylim(0.02, 1.30)
    model_ax.set_yticks([])
    model_ax.spines[["left", "right", "top"]].set_visible(False)
    model_ax.spines["bottom"].set_color("#888888")
    model_ax.plot([0, short_bp], [0.40, 0.40], color="#202020", linewidth=6, solid_capstyle="butt")
    model_ax.plot([0, short_bp], [0.82, 0.82], color="#202020", linewidth=5, solid_capstyle="butt")
    model_ax.plot(
        [short_bp, long_bp],
        [0.82, 0.82],
        color=target_color,
        linewidth=5,
        solid_capstyle="butt",
    )
    model_ax.scatter([0, short_bp], [0.40, 0.40], s=28, facecolor="white", edgecolor="#202020", zorder=4)
    model_ax.scatter([0, long_bp], [0.82, 0.82], s=28, facecolor="white", edgecolor=target_color, zorder=4)
    model_ax.axvline(short_bp, ymin=0.22, ymax=0.76, color=neutral, linestyle=(0, (3, 3)), linewidth=0.8)

    g0_fraction = float(candidate["target_00_endpoint_fraction_median"])
    g2_fraction = float(candidate["target_11_endpoint_fraction_median"])
    arrow_start = short_bp + extension_bp * g0_fraction
    arrow_end = short_bp + extension_bp * g2_fraction
    model_ax.annotate(
        "",
        xy=(arrow_end, 1.12),
        xytext=(arrow_start, 1.12),
        arrowprops={"arrowstyle": "-|>", "color": target_color, "linewidth": 2.3},
    )
    model_ax.text(arrow_start, 1.18, genotype_labels["0|0"], ha="center", va="bottom", fontsize=8.8, color=neutral)
    model_ax.text(arrow_end, 1.18, genotype_labels["1|1"], ha="center", va="bottom", fontsize=8.8, color=target_color, fontweight="bold")
    model_ax.text(
        (arrow_start + arrow_end) / 2,
        1.04,
        f"+{float(candidate['target_coverage_endpoint_delta_bp']):.0f} bp",
        ha="center",
        va="top",
        color=target_color,
        fontsize=9.4,
        fontweight="bold",
    )
    model_ax.text(0, 0.49, "short 3'UTR", ha="left", va="bottom", fontsize=9.3, fontweight="bold")
    model_ax.text(0, 0.91, "long 3'UTR", ha="left", va="bottom", fontsize=9.3, fontweight="bold")
    model_ax.set_xticks([0, short_bp, long_bp], ["0", str(short_bp), str(long_bp)])
    model_ax.set_xlabel("3'UTR length (bp)")
    model_ax.set_title(
        f"{gene}  |  {target_stage}",
        loc="left",
        color=target_color,
        fontweight="bold",
        pad=9,
    )
    add_panel_label(model_ax, row_index * 3)

    # Formal APA-QTL distributions: statistics come from the QTL screen.
    effect_ax = fig.add_subplot(grid[row_index, 1])
    offsets = {"Prenatal": -0.18, "Postnatal": 0.18}
    rng = np.random.default_rng(20260724 + row_index)
    legend_handles = []
    for stage in stage_order:
        stage_color = stage_colors[stage]
        stage_frame = gene_phenotype.loc[gene_phenotype["stage"] == stage]
        positions = np.arange(3) + offsets[stage]
        values_by_group = [
            stage_frame.loc[stage_frame["genotype"] == genotype, "apa_score"].to_numpy()
            for genotype in genotype_order
        ]
        for position, values in zip(positions, values_by_group):
            jitter = rng.uniform(-0.07, 0.07, size=len(values))
            effect_ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                s=4,
                color=stage_color,
                alpha=0.12 if stage == target_stage else 0.07,
                linewidths=0,
                rasterized=True,
                zorder=1,
            )
        boxes = effect_ax.boxplot(
            values_by_group,
            positions=positions,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": stage_color, "linewidth": 1.4},
            boxprops={"edgecolor": stage_color, "linewidth": 1.1},
            whiskerprops={"color": stage_color, "linewidth": 0.9},
            capprops={"color": stage_color, "linewidth": 0.9},
        )
        for box in boxes["boxes"]:
            box.set_facecolor(stage_color)
            box.set_alpha(0.18 if stage == target_stage else 0.07)
        means = [float(np.mean(values)) for values in values_by_group]
        effect_ax.plot(
            positions,
            means,
            color=stage_color,
            linewidth=2.1 if stage == target_stage else 1.0,
            alpha=1.0 if stage == target_stage else 0.58,
            marker="o",
            markersize=4.2 if stage == target_stage else 3.0,
            zorder=4,
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=stage_color,
                linewidth=2.1 if stage == target_stage else 1.0,
                alpha=1.0 if stage == target_stage else 0.58,
                marker="o",
                markersize=3.5,
                label=stage,
            )
        )

    effect_ax.set_xlim(-0.55, 2.55)
    effect_ax.set_ylim(global_y_min, global_y_max)
    effect_ax.set_xticks(range(3), [genotype_labels[genotype] for genotype in genotype_order])
    effect_ax.set_xlabel("Genotype")
    effect_ax.set_ylabel("Normalized PDUI")
    effect_ax.grid(axis="y", color="#E8E8E8", linewidth=0.45)
    effect_ax.set_axisbelow(True)
    effect_ax.legend(handles=legend_handles, loc="lower right", fontsize=8.6, handlelength=1.5)
    target_q = float(candidate["target_q"])
    contrast_fdr = float(candidate["p_effect_difference_fdr"])
    q_text = f"{target_q:.1e}" if target_q < 1e-4 else f"{target_q:.3f}"
    effect_ax.text(
        0.02,
        0.98,
        f"q={q_text}",
        transform=effect_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.1,
        color=target_color,
        fontweight="bold",
    )
    add_panel_label(effect_ax, row_index * 3 + 1)

    # Dual-stage WIG tracks: one group-median sample per stage and genotype.
    coverage_grid = grid[row_index, 2].subgridspec(2, 3, hspace=0.22, wspace=0.18)
    coverage_y_max = max(
        1.25,
        np.ceil(gene_coverage["normalized_coverage"].max() * 4) / 4 + 0.05,
    )
    gene_selection = wig_selection.loc[wig_selection["gene"] == gene].set_index(
        ["stage", "genotype"]
    )
    distal_means = {stage: {} for stage in stage_order}
    endpoint_fractions = {stage: [] for stage in stage_order}
    selected_samples = []

    for stage_index, stage in enumerate(stage_order):
        for genotype_index, genotype in enumerate(genotype_order):
            selection_row = gene_selection.loc[(stage, genotype)]
            sample = selection_row["sample"]
            curve = gene_coverage.loc[
                (gene_coverage["stage"] == stage)
                & (gene_coverage["genotype"] == genotype)
                & (gene_coverage["sample"] == sample)
            ].sort_values("position_along_3utr_bp")
            x = curve["position_along_3utr_bp"].to_numpy()
            coverage = curve["normalized_coverage"].to_numpy()
            distal_means[stage][genotype] = float(
                curve.loc[
                    curve["position_along_3utr_bp"] >= short_bp,
                    "normalized_coverage",
                ].mean()
            )
            endpoint_fraction = float(selection_row["endpoint_fraction"])
            endpoint_fractions[stage].append(endpoint_fraction)
            selected_samples.append(sample)

            coverage_ax = fig.add_subplot(coverage_grid[stage_index, genotype_index])
            color = genotype_colors[genotype]
            coverage_ax.fill_between(x, 0, coverage, color=color, alpha=0.22, linewidth=0)
            coverage_ax.plot(x, coverage, color=color, linewidth=1.15)
            coverage_ax.axvspan(short_bp, long_bp, color=stage_colors[stage], alpha=0.045, linewidth=0)
            coverage_ax.axvline(
                short_bp,
                color="#505050",
                linestyle=(0, (4, 3)),
                linewidth=0.75,
            )
            coverage_ax.set_xlim(0, long_bp)
            coverage_ax.set_ylim(0, coverage_y_max)
            coverage_ax.grid(axis="y", color="#ECECEC", linewidth=0.35)
            coverage_ax.set_axisbelow(True)
            coverage_ax.text(
                0.97,
                0.90,
                genotype_labels[genotype],
                transform=coverage_ax.transAxes,
                ha="right",
                va="top",
                fontsize=8.3,
                color=color,
                fontweight="bold",
            )
            if genotype_index == 0:
                coverage_ax.set_yticks([0, coverage_y_max], ["0", f"{coverage_y_max:.1f}"])
                coverage_ax.tick_params(axis="y", labelsize=8.0, length=2)
                coverage_ax.set_ylabel(stage, color=stage_colors[stage], fontweight="bold")
            else:
                coverage_ax.set_yticks([])
                coverage_ax.spines["left"].set_visible(False)
            if stage_index == 0:
                coverage_ax.tick_params(axis="x", labelbottom=False, length=0)
            else:
                coverage_ax.set_xticks([0, short_bp, long_bp])
                if genotype_index == 0:
                    coverage_ax.set_xticklabels(["0", str(short_bp), str(long_bp)])
                else:
                    coverage_ax.set_xticklabels(["", str(short_bp), str(long_bp)])
                coverage_ax.tick_params(axis="x", labelsize=7.8)
                if genotype_index == 1:
                    coverage_ax.set_xlabel("Position")
            if stage_index == 0 and genotype_index == 0:
                coverage_ax.set_title(
                    "Group Tracks",
                    loc="left",
                    color=target_color,
                    fontsize=9.0,
                    fontweight="bold",
                    pad=3,
                )
                coverage_ax.text(
                    -0.48,
                    -0.12,
                    "Normalized coverage",
                    transform=coverage_ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=8.2,
                    color=neutral,
                )
                add_panel_label(coverage_ax, row_index * 3 + 2)

    target_coverage_delta = (
        distal_means[target_stage]["1|1"] - distal_means[target_stage]["0|0"]
    )
    comparator_coverage_delta = (
        distal_means[comparator_stage]["1|1"]
        - distal_means[comparator_stage]["0|0"]
    )

    layout_rows.append(
        {
            "figure": "figure4_dynamic",
            "gene": gene,
            "variant_id": candidate["variant_id"],
            "target_stage": target_stage,
            "comparator_stage": comparator_stage,
            "phenotype_samples": int(len(gene_phenotype)),
            "short_utr_bp": short_bp,
            "long_utr_bp": long_bp,
            "distal_extension_bp": extension_bp,
            "target_q": target_q,
            "stage_effect_difference_fdr": contrast_fdr,
            "target_wig_endpoint_delta_bp": float(candidate["target_coverage_endpoint_delta_bp"]),
            "comparator_wig_endpoint_delta_bp": float(candidate["other_coverage_endpoint_delta_bp"]),
            "target_distal_coverage_delta": target_coverage_delta,
            "comparator_distal_coverage_delta": comparator_coverage_delta,
            "representative_track_count": int(sum(len(values) for values in endpoint_fractions.values())),
            "representative_unique_sample_count": int(len(set(selected_samples))),
            "representative_qc_pass_count": int(
                ((sample_qc["gene"] == gene) & (sample_qc["status"] == "PASS")).sum()
            ),
            "selection_max_abs_endpoint_deviation": float(
                (
                    gene_selection["endpoint_fraction"]
                    - gene_selection["group_median_endpoint_fraction"]
                ).abs().max()
            ),
            "representative_endpoint_monotonic": bool(
                np.all(np.diff(endpoint_fractions[target_stage]) > 0)
            ),
            "representative_distal_mean_monotonic": bool(
                distal_means[target_stage]["0|0"]
                < distal_means[target_stage]["0|1"]
                < distal_means[target_stage]["1|1"]
            ),
            "both_stages_present": set(gene_coverage["stage"]) == set(stage_order),
            "stage_genotype_cell_count": int(
                gene_coverage[["stage", "genotype"]].drop_duplicates().shape[0]
            ),
            "coverage_curve_clipping_count": int(
                (
                    gene_coverage["normalized_coverage"]
                    > coverage_y_max
                ).sum()
            ),
            "formal_target_significant": target_q < 0.05,
            "formal_stage_contrast_significant": contrast_fdr < 0.05,
            "target_longer_than_comparator": (
                float(candidate["target_coverage_endpoint_delta_bp"])
                > abs(float(candidate["other_coverage_endpoint_delta_bp"]))
            ),
            "target_coverage_direction_supported": target_coverage_delta > 0,
            "representative_tracks_are_inferential": False,
            "coverage_status": candidate["coverage_status"],
        }
    )


# --- layout checks ---

required_tokens = [
    "PACC1",
    "SLCO1A2",
    "short 3'UTR",
    "long 3'UTR",
    "3'UTR length (bp)",
    "q=",
    "Normalized coverage",
    "Normalized PDUI",
    "Position",
    "Group Tracks",
    "Prenatal",
    "Postnatal",
    "C/C",
    "C/T",
    "T/T",
    "T/C",
]
forbidden_tokens = [
    "exact cleavage",
    "exact PAS",
    "distal PAS",
    "APA-score high",
    "APA-score low",
    "bigWig",
    "BAM",
    "endpoint=",
    "APA score",
    "bp distal",
] + sorted(wig_selection["sample"].unique())
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
figure_box = fig.bbox
tick_text = {
    id(label): axis.bbox
    for axis in fig.axes
    for label in [*axis.get_xticklabels(), *axis.get_yticklabels()]
}
visible_text = []
overflow_labels = []
for text_object in fig.findobj(match=Text):
    label = text_object.get_text().strip()
    if not text_object.get_visible() or not label:
        continue
    box = text_object.get_window_extent(renderer=renderer)
    if box.width <= 0 or box.height <= 0:
        continue
    if id(text_object) in tick_text:
        axis_box = tick_text[id(text_object)]
        center_x = (box.x0 + box.x1) / 2
        center_y = (box.y0 + box.y1) / 2
        if not (
            axis_box.x0 - 4 <= center_x <= axis_box.x1 + 4
            and axis_box.y0 - 4 <= center_y <= axis_box.y1 + 4
        ):
            continue
    visible_text.append((label, box, text_object))
    if (
        box.x0 < figure_box.x0 - 3
        or box.y0 < figure_box.y0 - 3
        or box.x1 > figure_box.x1 + 3
        or box.y1 > figure_box.y1 + 3
    ):
        overflow_labels.append(label.replace("\n", " | "))

overlap_labels = []
for index, (first_label, first_box, first_object) in enumerate(visible_text):
    for second_label, second_box, second_object in visible_text[index + 1 :]:
        if first_object.axes is second_object.axes and first_object.axes is not None:
            if first_object in first_object.axes.get_xticklabels() or second_object in second_object.axes.get_xticklabels():
                continue
            if first_object in first_object.axes.get_yticklabels() or second_object in second_object.axes.get_yticklabels():
                continue
        width = max(0.0, min(first_box.x1, second_box.x1) - max(first_box.x0, second_box.x0))
        height = max(0.0, min(first_box.y1, second_box.y1) - max(first_box.y0, second_box.y0))
        if width * height > 16:
            overlap_labels.append(
                f"{first_label.replace(chr(10), ' | ')} <> "
                f"{second_label.replace(chr(10), ' | ')}"
            )

all_text = " ".join(label for label, _, _ in visible_text)
forbidden_hits = [token for token in forbidden_tokens if token.lower() in all_text.lower()]
layout_validation = {
    "text_overflow_count": len(overflow_labels),
    "overflow_labels": "; ".join(overflow_labels),
    "text_overlap_count": len(overlap_labels),
    "overlap_labels": "; ".join(overlap_labels),
    "required_text_present": all(token in all_text for token in required_tokens),
    "prohibited_text_count": len(forbidden_hits),
    "prohibited_text": "; ".join(forbidden_hits),
}


# --- save figure ---

fig.savefig(png_file, dpi=dpi, facecolor="white")
fig.savefig(pdf_file, facecolor="white")
fig.savefig(svg_file, facecolor="white")
png_pixels = plt.imread(png_file)
nonwhite_fraction = float(np.mean(np.any(png_pixels[..., :3] < 0.99, axis=2)))
plt.close(fig)


# --- layout status ---

layout = pd.DataFrame(layout_rows)
for key, value in layout_validation.items():
    layout[key] = value
layout["canvas_nonwhite_fraction"] = nonwhite_fraction
layout["canvas_nonblank"] = nonwhite_fraction > 0.01
layout["png_width_px"] = int(png_pixels.shape[1])
layout["png_height_px"] = int(png_pixels.shape[0])
layout["status"] = np.where(
    (layout["text_overflow_count"] == 0)
    & (layout["text_overlap_count"] == 0)
    & layout["required_text_present"]
    & (layout["prohibited_text_count"] == 0)
    & layout["formal_target_significant"]
    & layout["formal_stage_contrast_significant"]
    & layout["target_longer_than_comparator"]
    & layout["target_coverage_direction_supported"]
    & layout["representative_endpoint_monotonic"]
    & layout["both_stages_present"]
    & (layout["stage_genotype_cell_count"] == 6)
    & (layout["representative_track_count"] == 6)
    & (layout["representative_unique_sample_count"] == 6)
    & (layout["representative_qc_pass_count"] == 6)
    & (layout["selection_max_abs_endpoint_deviation"] < 0.01)
    & (layout["coverage_curve_clipping_count"] == 0)
    & (layout["coverage_status"] == "PASS_DIRECT_ENDPOINT")
    & layout["canvas_nonblank"],
    "PASS",
    "CHECK",
)
layout.to_csv(layout_file, sep="\t", index=False)

print(layout.to_string(index=False))
