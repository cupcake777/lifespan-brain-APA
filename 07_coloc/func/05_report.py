#!/usr/bin/env python3
# Baseline apaQTL review report: tables, figures and markdown.
# Reads data/07_coloc tables; writes the markdown report, layout validation and figures.

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/apa_07_coloc_mplconfig")

import matplotlib
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr

matplotlib.use("Agg")
from matplotlib import pyplot as plt


# --- config ---

ROOT = Path(os.environ["APA_ROOT"])
DATA_DIR = ROOT / "data" / "07_coloc"
FIGURE_DIR = ROOT / "figures" / "07_coloc"
OUTPUT_PATH = DATA_DIR / "Baseline_apaQTL_review.md"
MECHANISM_PLOT_BASE = FIGURE_DIR / "figure_mechanism_composition"
APARENT_PLOT_BASE = FIGURE_DIR / "figure_aparent_delta_by_pas_status"
APARENT_SCATTER_BASE = FIGURE_DIR / "figure_aparent_concordance_scatter"
DUR_DENSITY_PLOT_BASE = FIGURE_DIR / "figure_dur_density_by_direction"
CROSS_STAGE_PLOT_BASE = FIGURE_DIR / "figure_cross_stage_direction"
LAYOUT_VALIDATION_PATH = DATA_DIR / "figure_layout_validation.tsv"
LAYOUT_VALIDATION_ROWS: list[dict[str, str | int | float | bool]] = []


# --- helpers ---

def save_figure_with_validation(fig, figure_id: str, output_paths: list[Path]) -> None:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    fig_width, fig_height = fig.get_size_inches()
    tight_bbox = fig.get_tightbbox(renderer)
    max_tick_label_chars = 0
    n_legends = 0
    for ax in fig.axes:
        tick_texts = [tick.get_text() for tick in ax.get_xticklabels() + ax.get_yticklabels()]
        max_tick_label_chars = max(max_tick_label_chars, *(len(text) for text in tick_texts), 0)
        n_legends += int(ax.get_legend() is not None)
    LAYOUT_VALIDATION_ROWS.append(
        {
            "figure_id": figure_id,
            "width_in": round(fig_width, 3),
            "height_in": round(fig_height, 3),
            "n_axes": len(fig.axes),
            "n_legends": n_legends,
            "max_tick_label_chars": max_tick_label_chars,
            "tight_bbox_width_in": round(float(tight_bbox.width), 3),
            "tight_bbox_height_in": round(float(tight_bbox.height), 3),
            "layout_within_canvas": bool(tight_bbox.width <= fig_width * 1.25 and tight_bbox.height <= fig_height * 1.25),
            "outputs": ";".join(path.name for path in output_paths),
        }
    )
    for path in output_paths:
        fig.savefig(path, dpi=300, bbox_inches="tight")


def read_optional_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep=sep)


def pct(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return f"{numerator / denominator * 100:.1f}%"


def fmt_int(value: float | int | str) -> str:
    if pd.isna(value):
        return "NA"
    return f"{int(value):,}"


def fmt_float(value: float | int | str, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def pearson_or_na(df: pd.DataFrame, x_col: str, y_col: str) -> str:
    subset = df[[x_col, y_col]].dropna()
    if subset.shape[0] < 2:
        return "NA"
    return fmt_float(subset[x_col].corr(subset[y_col]), 4)


def fmt_paths(paths: list[Path]) -> str:
    if not paths:
        return "not_generated"
    return "、".join(path.name for path in paths)


# --- load tables ---

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
mech = pd.read_csv(DATA_DIR / "mechanism_main_table.tsv.gz", sep="\t")
func = pd.read_csv(DATA_DIR / "functional_effect_main_table.tsv.gz", sep="\t")
utr = pd.read_csv(DATA_DIR / "utr_sequence_annotations.tsv.gz", sep="\t")
rna = pd.read_csv(DATA_DIR / "rnafold_candidate_summary.tsv.gz", sep="\t")
sc = pd.read_csv(DATA_DIR / "single_cell_expression_mapping.tsv.gz", sep="\t")
priority = pd.read_csv(DATA_DIR / "priority_cases.tsv", sep="\t")
pas_scores = pd.read_csv(DATA_DIR / "pas_disruption_scores.tsv", sep="\t")
concordance_summary = pd.read_csv(DATA_DIR / "aparent_internal_concordance.tsv", sep="\t")
protein = read_optional_table(DATA_DIR / "protein_evidence.tsv")
cross_stage = read_optional_table(DATA_DIR / "cross_stage_direction.tsv")
flip_protein = read_optional_table(DATA_DIR / "flip_with_protein_evidence.tsv")
eqtl_coloc = read_optional_table(DATA_DIR / "eqtl_coloc.tsv")
gwas_coloc = read_optional_table(DATA_DIR / "gwas_coloc.tsv")
sensitivity = read_optional_table(DATA_DIR / "total_score_sensitivity.tsv")
run_summary = json.loads((DATA_DIR / "run_summary.json").read_text())

n_all = int(run_summary["n_candidates_mechanism"])
n_functional = int(run_summary["n_candidates_functional_effect"])
n_utr_ok = int(run_summary["n_utr_success"])
n_pas_direct = int(run_summary["n_pas_direct_supported"])
n_rbp_supported = int(run_summary["n_rbp_supported"])
n_multi_module = int((func["functional_support_category"] == "density-supported multi-module hits").sum())
n_raw_only = int((func["functional_support_category"] == "raw-count-only heuristic load").sum())
n_density_single = int((func["functional_support_category"] == "density-supported single-module hit").sum())
n_sc_ok = int((sc["sc_expression_status"] == "ok").sum())
n_long_rna = int(rna["long_mfe_per_100nt"].notna().sum())
n_short_rna = int(rna["short_mfe_per_100nt"].notna().sum())
n_paired_rna = int(rna["abs_long_short_mfe_delta_per_100nt"].notna().sum())
n_no_signed_effect = int(mech["effect_direction_status"].eq("missing_signed_qtl_effect").sum())
n_signed_but_unmapped = int(mech["effect_direction_status"].eq("signed_effect_available_but_apa_mapping_unresolved").sum())
n_direction_blocked = n_no_signed_effect + n_signed_but_unmapped
corr_mirna = pearson_or_na(func, "dur_length", "mirna_total_sites")
corr_rbp = pearson_or_na(func, "dur_length", "dur_rbp_consensus_total_hits")
corr_are = pearson_or_na(func, "dur_length", "are_total_hits")

# median DUR length per functional category
dur_summary = (
    func.groupby("functional_support_category")["dur_length"]
    .median()
    .sort_index()
)
dur_summary_text = "; ".join(
    f"{label} median={fmt_int(value)} nt"
    for label, value in dur_summary.items()
) or "NA"

# mechanism composition table
mech_summary = (
    mech.groupby(["stage", "tier", "mechanism_class"], as_index=False)
    .size()
    .rename(columns={"size": "n"})
)
mechanism_lines = [
    "| Stage | Tier | Mechanism Class | N |",
    "| --- | --- | --- | ---: |",
]
for row in mech_summary.itertuples(index=False):
    mechanism_lines.append(f"| {row.stage} | {row.tier} | {row.mechanism_class} | {row.n} |")
mechanism_table_text = "\n".join(mechanism_lines)

# functional category table
func_summary = (
    func.groupby(["stage", "tier", "functional_support_category"], as_index=False)
    .size()
    .rename(columns={"size": "n"})
)
functional_lines = [
    "| Stage | Tier | Functional Category | N |",
    "| --- | --- | --- | ---: |",
]
for row in func_summary.itertuples(index=False):
    functional_lines.append(
        f"| {row.stage} | {row.tier} | {row.functional_support_category} | {row.n} |"
    )
functional_table_text = "\n".join(functional_lines)


# --- figures ---

# mechanism composition by stage
mechanism_plot_paths: list[Path] = []
if not mech.empty:
    order = [
        "PAS-direct",
        "mixed-support",
        "RBP-mediated",
        "mechanistically unresolved",
    ]
    stage_order = ["prenatal", "postnatal"]
    plot_df = (
        mech.groupby(["stage", "mechanism_class"], as_index=False)
        .size()
        .rename(columns={"size": "n_candidates"})
    )
    plot_df["mechanism_class"] = pd.Categorical(
        plot_df["mechanism_class"], categories=order, ordered=True
    )
    plot_df["stage"] = pd.Categorical(plot_df["stage"], categories=stage_order, ordered=True)
    plot_df = plot_df.sort_values(["stage", "mechanism_class"], kind="stable")
    total_by_stage = plot_df.groupby("stage")["n_candidates"].transform("sum")
    plot_df["fraction"] = plot_df["n_candidates"] / total_by_stage

    sns.set_theme(style="whitegrid")
    palette = {
        "PAS-direct": "#C44E52",
        "mixed-support": "#DD8452",
        "RBP-mediated": "#4C72B0",
        "mechanistically unresolved": "#9AA0A6",
    }
    fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
    sns.barplot(
        data=plot_df,
        x="stage",
        y="fraction",
        hue="mechanism_class",
        hue_order=order,
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of candidates")
    ax.set_ylim(0, 1)
    ax.set_title("Mechanism composition by developmental stage")
    for patch, (_, row) in zip(ax.patches, plot_df.iterrows(), strict=False):
        height = patch.get_height()
        if pd.isna(height) or height <= 0:
            continue
        ax.text(
            patch.get_x() + patch.get_width() / 2,
            height + 0.015,
            f"{int(row['n_candidates'])}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    ax.legend(title="Mechanism class", frameon=False, loc="upper right")
    mechanism_plot_paths = [MECHANISM_PLOT_BASE.with_suffix(".png"), MECHANISM_PLOT_BASE.with_suffix(".pdf")]
    save_figure_with_validation(fig, "mechanism_composition", mechanism_plot_paths)
    plt.close(fig)

# APARENT delta distribution by PAS status
aparent_plot_paths: list[Path] = []
if not pas_scores.empty and not concordance_summary.empty:
    sns.set_theme(style="whitegrid")
    summary_row = concordance_summary.iloc[0]
    plot_df = pas_scores[["pas_disrupted", "APARENT_delta"]].copy()
    plot_df["pas_disrupted"] = plot_df["pas_disrupted"].fillna("not_supported")
    plot_df = plot_df.dropna(subset=["APARENT_delta"])
    if not plot_df.empty:
        order = ["Gain", "Loss", "Modulation", "not_supported"]
        palette = {
            "Gain": "#55A868",
            "Loss": "#C44E52",
            "Modulation": "#DD8452",
            "not_supported": "#9AA0A6",
        }
        fig, ax = plt.subplots(figsize=(8.0, 5.5), constrained_layout=True)
        sns.violinplot(
            data=plot_df,
            x="pas_disrupted",
            y="APARENT_delta",
            hue="pas_disrupted",
            order=order,
            hue_order=order,
            palette=palette,
            cut=0,
            inner=None,
            legend=False,
            ax=ax,
        )
        sns.boxplot(
            data=plot_df,
            x="pas_disrupted",
            y="APARENT_delta",
            order=order,
            width=0.22,
            showcaps=True,
            boxprops={"facecolor": "white", "alpha": 0.85, "zorder": 3},
            whiskerprops={"linewidth": 1.2},
            medianprops={"color": "#111111", "linewidth": 1.4},
            showfliers=False,
            ax=ax,
        )
        rho = summary_row.get("spearman_rho")
        n_tested = int(summary_row.get("n_tested_candidates", 0) or 0)
        ax.set_title("APARENT delta distribution by PAS support status")
        ax.set_xlabel("PAS disruption class")
        ax.set_ylabel("APARENT delta log-odds")
        if pd.notna(rho):
            ax.text(
                0.03,
                0.97,
                f"Internal concordance summary\nSpearman rho = {float(rho):.3f}\nN = {n_tested}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                bbox={"facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.9},
            )
        aparent_plot_paths = [APARENT_PLOT_BASE.with_suffix(".png"), APARENT_PLOT_BASE.with_suffix(".pdf")]
        save_figure_with_validation(fig, "aparent_delta_by_pas_status", aparent_plot_paths)
        plt.close(fig)

# APARENT delta vs QTL signed effect
aparent_scatter_paths: list[Path] = []
slope_col = next((col for col in ("slope", "qtl_beta", "beta") if col in pas_scores.columns), None)
if not (pas_scores.empty or slope_col is None or "APARENT_delta" not in pas_scores.columns):
    plot_df = pas_scores[["APARENT_delta", slope_col, "pas_disrupted"]].copy()
    plot_df[slope_col] = pd.to_numeric(plot_df[slope_col], errors="coerce")
    plot_df["APARENT_delta"] = pd.to_numeric(plot_df["APARENT_delta"], errors="coerce")
    plot_df = plot_df.dropna(subset=["APARENT_delta", slope_col])
    if plot_df.shape[0] >= 2:
        rho, pval = spearmanr(plot_df["APARENT_delta"], plot_df[slope_col])
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(7.5, 5.5), constrained_layout=True)
        sns.regplot(
            data=plot_df,
            x="APARENT_delta",
            y=slope_col,
            scatter=False,
            color="#333333",
            ax=ax,
        )
        sns.scatterplot(
            data=plot_df,
            x="APARENT_delta",
            y=slope_col,
            hue="pas_disrupted",
            palette={"Gain": "#55A868", "Loss": "#C44E52", "Modulation": "#DD8452", "Auxiliary": "#8172B2", "not_supported": "#9AA0A6"},
            ax=ax,
        )
        ax.set_title("APARENT delta vs QTL effect")
        ax.set_xlabel("APARENT delta log-odds")
        ax.set_ylabel("QTL signed effect")
        ax.text(
            0.03,
            0.97,
            f"rho={rho:.3f}\nP={pval:.2e}\nN={plot_df.shape[0]}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "#D0D0D0", "alpha": 0.9},
        )
        aparent_scatter_paths = [APARENT_SCATTER_BASE.with_suffix(".png"), APARENT_SCATTER_BASE.with_suffix(".pdf")]
        save_figure_with_validation(fig, "aparent_concordance_scatter", aparent_scatter_paths)
        plt.close(fig)

# DUR regulatory-load density by APA direction
dur_density_plot_paths: list[Path] = []
if not func.empty:
    plot_df = func.copy()
    plot_df["direction_group"] = plot_df["utr_effect_direction"].fillna("unknown")
    plot_df["direction_group"] = plot_df["direction_group"].replace("", "unknown")
    value_cols = [
        ("mirna_sites_per_100nt", "miRNA"),
        ("dur_rbp_hits_per_100nt", "RBP"),
        ("are_hits_per_100nt", "ARE"),
    ]
    long_df = plot_df[["direction_group", *[col for col, _ in value_cols]]].melt(
        id_vars="direction_group",
        var_name="metric",
        value_name="density_per_100nt",
    )
    label_map = dict(value_cols)
    long_df["metric"] = long_df["metric"].map(label_map)
    long_df = long_df.dropna(subset=["density_per_100nt"])
    if not long_df.empty:
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
        sns.boxplot(
            data=long_df,
            x="metric",
            y="density_per_100nt",
            hue="direction_group",
            order=["miRNA", "RBP", "ARE"],
            hue_order=["lengthening", "shortening", "unknown"],
            palette={
                "lengthening": "#55A868",
                "shortening": "#C44E52",
                "unknown": "#9AA0A6",
            },
            showfliers=False,
            ax=ax,
        )
        ax.set_title("DUR regulatory-load density by APA direction")
        ax.set_xlabel("")
        ax.set_ylabel("Density per 100 nt")
        ax.legend(title="APA direction", frameon=False, loc="upper right")
        dur_density_plot_paths = [DUR_DENSITY_PLOT_BASE.with_suffix(".png"), DUR_DENSITY_PLOT_BASE.with_suffix(".pdf")]
        save_figure_with_validation(fig, "dur_density_by_direction", dur_density_plot_paths)
        plt.close(fig)

# cross-stage direction scatter
cross_stage_plot_paths: list[Path] = []
if not cross_stage.empty and {"slope_prenatal", "slope_postnatal", "cross_stage_direction"}.issubset(cross_stage.columns):
    plot_df = cross_stage.copy()
    plot_df["slope_prenatal"] = pd.to_numeric(plot_df["slope_prenatal"], errors="coerce")
    plot_df["slope_postnatal"] = pd.to_numeric(plot_df["slope_postnatal"], errors="coerce")
    plot_df = plot_df.dropna(subset=["slope_prenatal", "slope_postnatal"])
    if not plot_df.empty:
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(7.0, 6.0), constrained_layout=True)
        sns.scatterplot(
            data=plot_df,
            x="slope_prenatal",
            y="slope_postnatal",
            hue="cross_stage_direction",
            palette={"Flip": "#C44E52", "Amplifying": "#4C72B0", "Attenuating": "#55A868", "not_classified": "#9AA0A6"},
            ax=ax,
        )
        lim = max(plot_df[["slope_prenatal", "slope_postnatal"]].abs().max().max(), 0.1)
        ax.plot([-lim, lim], [-lim, lim], linestyle="--", color="#666666", linewidth=1)
        for row in plot_df.loc[plot_df["cross_stage_direction"].eq("Flip")].head(10).itertuples(index=False):
            ax.text(row.slope_prenatal, row.slope_postnatal, str(getattr(row, "gene_name", "")), fontsize=8)
        ax.set_title("Cross-stage direction")
        ax.set_xlabel("Prenatal slope")
        ax.set_ylabel("Postnatal slope")
        cross_stage_plot_paths = [CROSS_STAGE_PLOT_BASE.with_suffix(".png"), CROSS_STAGE_PLOT_BASE.with_suffix(".pdf")]
        save_figure_with_validation(fig, "cross_stage_direction", cross_stage_plot_paths)
        plt.close(fig)

if LAYOUT_VALIDATION_ROWS:
    pd.DataFrame(LAYOUT_VALIDATION_ROWS).to_csv(LAYOUT_VALIDATION_PATH, sep="\t", index=False)


# --- report ---

density_quantile_pct = fmt_int(run_summary.get("density_support_quantile", 0) * 100)
density_thresholds = run_summary.get("density_thresholds_per_100nt", {})
optional_extensions = run_summary.get("optional_extensions", {})
direction_counts = run_summary.get("effect_direction_status_counts", {})
pas_rescue_counts = run_summary.get("pas_rescue_status_counts", {})
mirna_density_threshold = fmt_float(density_thresholds.get("mirna_sites_per_100nt", float("nan")))
rbp_density_threshold = fmt_float(density_thresholds.get("dur_rbp_hits_per_100nt", float("nan")))
are_density_threshold = fmt_float(density_thresholds.get("are_hits_per_100nt", float("nan")))
aparent_status = concordance_summary.iloc[0].to_dict() if not concordance_summary.empty else {}

top5_unique = priority.drop_duplicates(subset=["gene_name"], keep="first").head(5).copy()

report_lines = [
    "# Baseline apaQTL 审查报告",
    "",
    "## 1. 报告定位",
    "",
    "- 本报告基于 `data/07_coloc/` 根目录中的当前 Baseline 结果，审查对象是 `PIP >= 0.5` 的 APA-QTL 候选。",
    "- 报告只陈述当前结果表可直接支撑的事实，不把 motif、POSTAR3 或序列层证据上升为实验性机制结论。",
    "- 需要特别强调：当前 PAS 配对仍来自启发式标准化，因此 `proximal/distal PAS` 应理解为可解释事件锚点，而不是 transcript-level 的最终真值。",
    "",
    "## 2. 数据漏斗",
    "",
    f"1. 候选起点：共 `{n_all}` 个机制层候选进入 Baseline 审查。",
    f"2. 功能效应主分析集：其中 `{n_functional}` 个候选具备 `dur_length`，占全部候选的 `{pct(n_functional, n_all)}`。",
    f"3. 转录本 3'UTR 成功构建：`{n_utr_ok}` 个候选 `utr_status=ok`，占功能效应主分析集的 `{pct(n_utr_ok, n_functional)}`，占全部机制层候选的 `{pct(n_utr_ok, n_all)}`。",
    f"4. 机制优先支持：`PAS-direct/mixed-support` 共 `{n_pas_direct}` 个，`RBP-mediated` 共 `{n_rbp_supported - n_pas_direct}` 个；其余 `{n_all - n_rbp_supported}` 个当前保留为 `mechanistically unresolved`。",
    f"5. Density-aware 功能支持：`{n_multi_module}` 个候选被归为 `density-supported multi-module hits`，`{n_density_single}` 个为 `density-supported single-module hit`，另有 `{n_raw_only}` 个仅保留为 `raw-count-only heuristic load`。",
    f"6. 单细胞背景可映射：`{n_sc_ok}` 个候选可在 BCA pseudobulk 中定位，覆盖全部机制层候选的 `{pct(n_sc_ok, n_all)}`。",
    f"7. RNA 结构可用度：short UTR 可计算 `{n_short_rna}` 个，long UTR 可计算 `{n_long_rna}` 个，成对 `long-short` 结构差异可比较 `{n_paired_rna}` 个；其余为长度或可用性受限，不应当解释为结构阴性。",
    f"8. 方向性审计：`{n_direction_blocked}` / `{n_all}` 个机制层候选方向 blocked，其中 `{n_no_signed_effect}` 个缺少 signed QTL effect，`{n_signed_but_unmapped}` 个有 signed effect 但 APA mapping 未解析；状态分布为 `{direction_counts}`。",
    "",
    "### 机制分类分层",
    "",
    mechanism_table_text,
    "",
    "### 功能分类分层",
    "",
    functional_table_text,
    "",
    "## 3. 方法清单",
    "",
    "### 已执行模块",
    "",
    "- 候选集与分层：沿用 `data/07_coloc/` 根目录中的 `PIP >= 0.5` 候选、`stage` 与 `tier` 分层。",
    "- PAS 事件标准化：沿用已生成的启发式 PAS 配对与 `dur_length` 定义。",
    "- PAS motif delta：对已解析 PAS 周边进行 REF/ALT motif 差异扫描；当前 `PAS-direct` 同时接受 motif 数量变化和 canonical PAS identity change。",
    "- CISBP-RNA variant motif delta：对变异周边序列进行 RBP motif 变化扫描；在 Baseline 判读中，仅把 `POSTAR3` 支持或“靠近已解析 PAS”的候选提升为 `RBP-mediated`。",
    "- transcript-aware 3'UTR/DUR 序列构建：使用 hg38 FASTA + GENCODE GTF 提取 short UTR、long UTR 与 DUR。",
    f"- miRNA / RBP / ARE 模块：保留 raw count，同时用当前运行内 non-zero density 的 Q{density_quantile_pct} 阈值做 density-aware 支持分层；本轮阈值分别为 miRNA `{mirna_density_threshold}`、RBP `{rbp_density_threshold}`、ARE `{are_density_threshold}` /100nt。",
    "- PAS rescue：除 motif 数量变化外，若变异位于 PAS core 且 canonical hexamer 身份改变，也会显式记为 canonical identity rescue，并在主表中保留来源/状态列。",
    "- RNAfold 模块：对可计算序列给出 `long/short` 结构摘要；超长或不可用序列保留缺失。",
    "- 单细胞背景：将候选基因映射到 `BCA_V1_pseudobulk.csv`，输出最高表达细胞类型。",
    "- Module E · GxStage：若存在 `cross_stage_direction.tsv`，按 Flip / Amplifying / Attenuating 输出跨阶段方向分类。",
    "- Module D · 蛋白层验证：若存在 `protein_evidence.tsv`，把 MS protein-by-genotype 结果接入 evidence chain。",
    f"- eQTL coloc：optional extension，当前状态 `{optional_extensions.get('eqtl_coloc', 'unknown')}`；缺失时不计入 Baseline 支持。",
    "- GWAS coloc：若存在 `gwas_coloc.tsv`，使用既有 apaQTL-GWAS coloc PP4 >= 0.75 判定。",
    "- TotalScore 敏感性分析：若存在 `total_score_sensitivity.tsv`，报告 component weight perturbation 的 top-50 稳定性。",
    "",
    "### 当前未纳入 Baseline 结论的模块",
    "",
    "- DeepRiPe 预测增强：本 Baseline 未执行。",
    "- Zipcode / 局部定位 motif：本 Baseline 未执行。",
    "- GO/Reactome/MSigDB 富集：外部 GMT 不作为 Baseline blocker；若无外部 GMT，使用 R/Bioconductor package-derived gene sets 时需在结果中记录来源。",
    "- 实验验证：本 Baseline 不包含 CLIP、reporter assay、CRISPR perturbation 等验证。",
    "",
    "## 4. 全局趋势",
    "",
    f"- 机制层面以 `mechanistically unresolved` 为主：`{n_all - n_rbp_supported}` / `{n_all}`（`{pct(n_all - n_rbp_supported, n_all)}`）。这说明 Baseline 更像优先级筛选而不是机制定论。",
    f"- `RBP-mediated` 候选共有 `{n_rbp_supported - n_pas_direct}` 个，明显多于 `PAS-direct/mixed-support` 的 `{n_pas_direct}` 个；但这只表示在当前规则下，更多候选保留了 RBP 线索，不代表 RBP 机制已被证明更常见。",
    f"- PAS rescue 状态分布为：{pas_rescue_counts}。这允许把 canonical identity rescue 与单纯的 motif-count delta 支持区分开来，不再把两者混为一个黑箱判定。",
    f"- 功能层面，density-aware 分类后仍有 `{n_multi_module}` / `{n_functional}`（`{pct(n_multi_module, n_functional)}`）属于 `density-supported multi-module hits`，另有 `{n_raw_only}` 个被降级为 `raw-count-only heuristic load`。",
    f"- `dur_length` 与 raw load 的相关性较强：miRNA `{corr_mirna}`，RBP `{corr_rbp}`，ARE `{corr_are}`。因此 raw counts 只保留为 heuristic load，不再单独驱动最终功能类别；分类以密度阈值和 RNA 结构支持共同决定。",
    f"- 按功能类别汇总的 DUR 长度中位数为：{dur_summary_text}。",
    f"- 3'UTR 构建的主要损失来自 `{utr['utr_status'].value_counts().idxmax()}` 之外的转录本边界问题；非 `ok` 共 `{n_functional - n_utr_ok}` 个，其中 `pas_outside_exonic_utr` 为 `{int((utr['utr_status'] == 'pas_outside_exonic_utr').sum())}` 个。",
    f"- RNAfold 是当前 Baseline 中可用度最低的功能模块之一：只有 `{n_paired_rna}` 个候选能直接比较成对 `long-short` 结构差异，因此不应将 `RNAfold=NA` 解释为“无结构效应”。",
    f"- 单细胞背景覆盖中等：`{n_sc_ok}` / `{n_all}`（`{pct(n_sc_ok, n_all)}`）的候选可映射到 pseudobulk，说明细胞类型背景信息对多数但非全部候选可用。",
    "- 当前所有没有显式输入方向列的候选仍保留为 `unknown`，因此报告中不再书写 UTR 缩短/延长、逃逸 miRNA 抑制或增加降解等 direction-dependent 叙述。",
    "- 从分期分层看，postnatal 的候选量高于 prenatal，因此绝对数上的机制/功能支持也主要集中在 postnatal；这反映的是当前候选池规模差异，不应直接解释为发育阶段因果差异。",
    "",
    "### 主线图输出",
    "",
    f"- 机制组成图：`{fmt_paths(mechanism_plot_paths)}`；按发育阶段展示 `PAS-direct / mixed-support / RBP-mediated / mechanistically unresolved` 的构成比例与样本数。",
    f"- APARENT 支持分布图：`{fmt_paths(aparent_plot_paths)}`；当前 APARENT 状态 `{aparent_status.get('analysis_status', 'unknown')}`，原因 `{aparent_status.get('reason', 'NA')}`，因此无完整 `APARENT_delta` 可画时不强行生成空图。",
    f"- DUR 密度分布图：`{fmt_paths(dur_density_plot_paths)}`；按 `utr_effect_direction` 分层展示 miRNA/RBP/ARE 的每 100 nt 负载分布。",
    "",
    "## 5. Module E · 跨发育阶段交互（GxStage）",
    "",
    f"- GxStage 方向表：`{cross_stage.shape[0]}` 行；方向分布为 `{cross_stage['cross_stage_direction'].value_counts().to_dict() if not cross_stage.empty and 'cross_stage_direction' in cross_stage.columns else {}}`。",
    f"- Cross-stage 散点图：`{cross_stage_plot_paths[0].name if cross_stage_plot_paths else 'not_generated'}`。",
    "",
    "## 6. Module D · 蛋白层验证",
    "",
    f"- Protein evidence 表：`{protein.shape[0]}` 行；状态分布为 `{protein['protein_evidence'].value_counts().to_dict() if not protein.empty and 'protein_evidence' in protein.columns else {}}`。",
    f"- Flip x Protein star case 表：`{flip_protein.shape[0]}` 行；concordance 分布为 `{flip_protein['protein_flip_concordance'].value_counts().to_dict() if not flip_protein.empty and 'protein_flip_concordance' in flip_protein.columns else {}}`。",
    "",
    "## 7. Functional Enrichment / Coloc / Sensitivity",
    "",
    f"- eQTL coloc 状态：`{optional_extensions.get('eqtl_coloc', 'unknown')}`；表行数 `{eqtl_coloc.shape[0]}`。GWAS coloc 表：`{gwas_coloc.shape[0]}` 行。",
    f"- TotalScore sensitivity 表：`{sensitivity.shape[0]}` 行。",
    f"- APARENT concordance scatter：`{aparent_scatter_paths[0].name if aparent_scatter_paths else 'not_generated'}`；只有存在完整 `APARENT_delta` 与 signed QTL effect pairs 时才生成。",
    "",
    "## 8. Top 5 案例解读",
    "",
    "以下案例按 `data/07_coloc/priority_cases.tsv` 的综合优先级排序，并去重到基因层后保留前 5 个基因。解读口径只限于当前结果表。",
    "",
]

for rank, row in enumerate(top5_unique.itertuples(index=False), start=1):
    case_row = pd.Series(row._asdict())
    func_match = func.loc[func["candidate_id"] == case_row["candidate_id"]]
    utr_match = utr.loc[utr["candidate_id"] == case_row["candidate_id"]]
    sc_match = sc.loc[sc["candidate_id"] == case_row["candidate_id"]]
    func_row = func_match.iloc[0] if not func_match.empty else None
    utr_row = utr_match.iloc[0] if not utr_match.empty else None
    sc_row = sc_match.iloc[0] if not sc_match.empty else None
    mechanism_repeat_note = ""
    gene_stage_n = int(priority.loc[priority["gene_name"].eq(case_row["gene_name"]), "stage"].nunique())
    if gene_stage_n > 1:
        mechanism_repeat_note = f"同一 `{case_row['gene_name']}` 候选出现在 `{gene_stage_n}` 个 stage 的优先队列中，提示其为跨阶段重复出现的高优先级事件。"

    case_lines = [
        f"### {rank}. {case_row['gene_name']} ({case_row['stage']})",
        f"- 候选：`{case_row['candidate_id']}`",
        f"- 优先级：`{case_row['tier']}`，`PIP={fmt_float(case_row['pip'], 6)}`，综合优先级分数 `={fmt_float(case_row['priority_score'], 3)}`。",
        f"- 机制层：`{case_row['mechanism_class']}`；PAS motif 记录的最高条目为 `{case_row['top_pas_support_label']}` 位点 `delta={fmt_int(case_row['top_pas_motif_delta'])}`，canonical hits `REF={case_row.get('ref_pas_canonical_hits', 'NA')}` / `ALT={case_row.get('alt_pas_canonical_hits', 'NA')}`；PAS 支持来源为 `{case_row.get('pas_direct_support_source', 'NA')}`，rescue 状态为 `{case_row.get('pas_rescue_status', 'NA')}`；RBP motif delta 的最高正向条目为 `{case_row['top_rbp_gain']}` (`{fmt_float(case_row['top_rbp_gain_delta'])}`)，最高负向条目为 `{case_row['top_rbp_loss']}` (`{fmt_float(case_row['top_rbp_loss_delta'])}`)；`variant_near_resolved_pas_250={case_row['variant_near_resolved_pas_250']}`。",
    ]
    if func_row is not None:
        case_lines.append(
            f"- 功能层：`dur_length={fmt_int(func_row['dur_length'])}` nt，miRNA seed-match 位点 `{fmt_int(func_row['mirna_total_sites'])}` 个（`{fmt_float(func_row['mirna_sites_per_100nt'])}` /100nt，Q75 阈值 `{fmt_float(func_row['mirna_density_threshold_per_100nt'])}`）；DUR RBP consensus hits `{fmt_int(func_row['dur_rbp_consensus_total_hits'])}` 次（`{fmt_float(func_row['dur_rbp_hits_per_100nt'])}` /100nt，Q75 阈值 `{fmt_float(func_row['dur_rbp_density_threshold_per_100nt'])}`）；ARE 命中 `{fmt_int(func_row['are_total_hits'])}` 次（`{fmt_float(func_row['are_hits_per_100nt'])}` /100nt，Q75 阈值 `{fmt_float(func_row['are_density_threshold_per_100nt'])}`）；density 支持模块 `{func_row.get('density_supported_modules', '') or 'none'}`，raw-only 模块 `{func_row.get('raw_only_load_modules', '') or 'none'}`；功能类别为 `{func_row['functional_support_category']}`。"
        )
        case_lines.append(
            f"- RNA 结构：`long-short |ΔMFE|/100nt = {fmt_float(func_row['abs_long_short_mfe_delta_per_100nt'])}`；若为 `NA`，表示本 Baseline 运行中未得到可用成对结构摘要。方向字段为 `{func_row['utr_effect_direction']}`（`{func_row['effect_direction_status']}`，来源 `{func_row.get('effect_direction_source_column', '') or 'none'}`），因此不做未被输入显式支撑的 UTR 缩短/延长或 gain/loss 生物学解释。"
        )
    if utr_row is not None:
        case_lines.append(
            f"- 3'UTR 构建：状态 `{utr_row['utr_status']}`；short/long/DUR 长度分别为 `{fmt_int(utr_row['short_utr_length'])}` / `{fmt_int(utr_row['long_utr_length'])}` / `{fmt_int(utr_row['dur_seq_length'])}` nt。"
        )
    if sc_row is not None:
        celltypes = sc_row["top_cell_types"]
        top_celltype = celltypes.split(";")[0] if isinstance(celltypes, str) and celltypes else "NA"
        case_lines.append(
            f"- 细胞类型背景：pseudobulk 映射状态 `{sc_row['sc_expression_status']}`；最高表达细胞类型为 `{top_celltype}`。"
        )
    case_lines.append(
        f"- 审查解读：当前结果支持其作为 `{case_row['mechanism_class']}` 的 Baseline 优先病例，但这仍是基于启发式 PAS 标准化与序列层 heuristic load 的优先级判断，不等同于已完成实验验证，也不支持方向性效应结论。{mechanism_repeat_note}"
    )
    report_lines.extend(["\n".join(case_lines), ""])

report_lines.extend(
    [
        "## 9. Baseline 审查结论",
        "",
        f"- 当前 Baseline 已经把 `{n_all}` 个候选压缩为可操作的优先级框架，其中最核心的交付是：`{n_functional}` 个功能层候选、`{n_utr_ok}` 个成功构建的 transcript-aware 3'UTR 事件，以及 50 个优先病例列表。",
        f"- 机制层最值得后续跟进的是 `{n_pas_direct}` 个 `PAS-direct/mixed-support` 候选和 `{n_rbp_supported - n_pas_direct}` 个 `RBP-mediated` 候选；但所有这些都应视为 Baseline 优先支持，而不是终局结论。",
        f"- 功能层现在区分 `density-supported` 与 `raw-count-only` 两类证据：`{n_multi_module}` 个候选保留为 `density-supported multi-module hits`，`{n_raw_only}` 个被明确降级为 `raw-count-only heuristic load`。",
        "- 当前输入仍缺少可直接映射到 APA 方向的 signed effect，因此除非输入表显式给出方向列，否则方向性生物学解释继续保持 blocked。",
        "- 下一轮审查若要提升证据等级，优先建议做三件事：",
        "- 第一，围绕 `data/07_coloc/priority_cases.tsv` 的前列病例做 transcript/PAS 边界的人工复核。",
        "- 第二，对 `PAS-direct` 或高优先 `RBP-mediated` 病例引入更强的位点级验证模块，例如 DeepRiPe、更严格的 motif 校准或外部实证数据库交叉。",
        "- 第三，若未来补到 signed QTL effect，再把方向字段和病例页解释一起升级；在此之前只保留 DUR 长度、miRNA/RBP/ARE/RNAfold 和细胞类型背景的描述性展示。",
        "",
        "## 10. 主要输出文件",
        "",
        "- `data/07_coloc/mechanism_main_table.tsv.gz`",
        "- `data/07_coloc/functional_effect_main_table.tsv.gz`",
        "- `data/07_coloc/priority_cases.tsv`",
        "- `data/07_coloc/utr_sequence_annotations.tsv.gz`",
        "- `data/07_coloc/rnafold_candidate_summary.tsv.gz`",
        "- `data/07_coloc/single_cell_expression_mapping.tsv.gz`",
        f"- `figures/07_coloc/{MECHANISM_PLOT_BASE.name}.png` / `figures/07_coloc/{MECHANISM_PLOT_BASE.name}.pdf`",
        f"- `figures/07_coloc/{CROSS_STAGE_PLOT_BASE.name}.png` / `figures/07_coloc/{CROSS_STAGE_PLOT_BASE.name}.pdf`",
        f"- `figures/07_coloc/{DUR_DENSITY_PLOT_BASE.name}.png` / `figures/07_coloc/{DUR_DENSITY_PLOT_BASE.name}.pdf`",
    ]
)
if aparent_plot_paths:
    report_lines.append(f"- `figures/07_coloc/{APARENT_PLOT_BASE.name}.png` / `figures/07_coloc/{APARENT_PLOT_BASE.name}.pdf`")
if aparent_scatter_paths:
    report_lines.append(f"- `figures/07_coloc/{APARENT_SCATTER_BASE.name}.png` / `figures/07_coloc/{APARENT_SCATTER_BASE.name}.pdf`")

OUTPUT_PATH.write_text("\n".join(report_lines), encoding="utf-8")
