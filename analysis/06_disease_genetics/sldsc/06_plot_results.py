#!/usr/bin/env python3
# Enrichment figures 1-8 and per-annotation bars for the sLDSC enrichment table.
# Recomputes FDR within each population x trait family and rewrites the summary table in place, as before.
# Writes figures/06_disease/*.png|pdf and data/06_disease/sldsc_enrichment_processed.tsv.
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --- paths ---

APA = os.environ['APA_ROOT']
SUMMARY_PATH = Path(APA, 'data/06_disease/sldsc_enrichment.tsv')
LD_SCORES_DIR = Path(APA, 'data/06_disease/sldsc_run/ldscore')
DATA_DIR = Path(APA, 'data/06_disease')
OUTPUT_DIR = Path(APA, 'figures/06_disease')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --- style and constants ---

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
})

OKABE_ITO = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
             '#0072B2', '#D55E00', '#CC79A7', '#000000']

ANNOT_NAMES = {
    'prenatal_apa_qtl': 'Prenatal APA-QTL',
    'postnatal_apa_qtl': 'Postnatal APA-QTL',
    'lifespan_apa_qtl': 'Lifespan APA-QTL',
    'interaction_apa_qtl': 'Interaction APA-QTL',
}

ANNOT_SHORT_NAMES = {
    'prenatal_apa_qtl': 'Prenatal',
    'postnatal_apa_qtl': 'Postnatal',
    'lifespan_apa_qtl': 'Lifespan',
    'interaction_apa_qtl': 'Interaction',
}

ANNOT_ORDER = [
    'prenatal_apa_qtl',
    'postnatal_apa_qtl',
    'lifespan_apa_qtl',
    'interaction_apa_qtl',
]

MAIN_ANNOT_ORDER = [
    'prenatal_apa_qtl',
    'postnatal_apa_qtl',
    'lifespan_apa_qtl',
]

# centred diverging scale: enrichment = 1 (no enrichment) is white; capped at -5..5
ENRICHMENT_VMIN = -5
ENRICHMENT_VCENTER = 1
ENRICHMENT_VMAX = 5

TRAIT_NAMES = {
    'AD': 'AD', 'ADHD': 'ADHD', 'AIS': 'AIS', 'ALS': 'ALS', 'ANX': 'ANX',
    'ASD': 'ASD', 'BIP': 'BIP', 'INS': 'Insomnia', 'INT': 'INT', 'MDD': 'MDD',
    'PD': 'PD', 'PTSD': 'PTSD', 'SCZ': 'SCZ', 'EP': 'Epilepsy',
}

POP_NAMES = {
    'EUR': 'European (EUR)',
    'EAS': 'East Asian (EAS)',
}

POP_COLORS = {'EUR': OKABE_ITO[4], 'EAS': OKABE_ITO[1]}

ANNOT_COLORS = {
    'prenatal_apa_qtl': OKABE_ITO[0],
    'postnatal_apa_qtl': OKABE_ITO[4],
    'lifespan_apa_qtl': OKABE_ITO[2],
    'interaction_apa_qtl': OKABE_ITO[7],
}

# --- helpers ---

def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    out = np.ones(n)
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n + 1)
    sorted_p = p[order]
    prev = 1.0
    for i in range(n - 1, -1, -1):
        v = sorted_p[i] * n / ranks[order[i]]
        v = min(v, prev)
        out[order[i]] = v
        prev = v
    return out


def fdr_stars(p_value):
    if p_value is None or pd.isna(p_value):
        return ''
    if p_value < 0.001:
        return '***'
    if p_value < 0.01:
        return '**'
    if p_value < 0.05:
        return '*'
    return ''


def add_panel_label(ax, label):
    ax.text(-0.12, 1.05, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top', ha='left',
            clip_on=False)


def clean_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def save_figure(fig, output_path):
    output_path = Path(output_path)
    fig.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path.with_suffix('.pdf'), dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


def pop_legend_handles():
    return [
        plt.Line2D([0], [0], marker='o', linestyle='none', markersize=5,
                   markerfacecolor=POP_COLORS['EUR'], markeredgecolor='black',
                   label='EUR'),
        plt.Line2D([0], [0], marker='o', linestyle='none', markersize=5,
                   markerfacecolor=POP_COLORS['EAS'], markeredgecolor='black',
                   label='EAS'),
    ]


def sorted_traits_by_signal(pop_df):
    trait_rank = (pop_df.groupby('Trait')['FDR_p'].min().reset_index()
                  .sort_values(['FDR_p', 'Trait']))
    return trait_rank['Trait'].tolist()


def signal_row_order(sub):
    rank = sub.sort_values(['FDR_p', 'Population', 'Trait'])
    return rank['PopTrait'].tolist()


def forest_panel_limits(sub, value_col, se_col, mult):
    lo = (sub[value_col] - mult * sub[se_col]).min()
    hi = (sub[value_col] + mult * sub[se_col]).max()
    span = hi - lo
    return lo - 0.05 * span, hi + 0.05 * span

# --- load summary: recompute FDR within population x trait families ---

df = pd.read_csv(SUMMARY_PATH, sep='\t')
df.columns = df.columns.str.strip()

numeric_cols = [
    'Coefficient', 'Coefficient_std_error', 'Coefficient_z-score',
    'Enrichment', 'Enrichment_std_error', 'Enrichment_p', 'FDR_p',
    'FDR_p_global',
]
for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

if 'FDR_p_global' not in df.columns:
    df['FDR_p_global'] = df['FDR_p']
df['FDR_p'] = df.groupby(['Population', 'Trait'], sort=False)[
    'Enrichment_p'].transform(bh_fdr)
df['n_tests_family'] = df.groupby(['Population', 'Trait'])[
    'Annotation'].transform('size')
if (len(df) != 76 or not (df['n_tests_family'] == 4).all() or
        df.duplicated(['Population', 'Trait', 'Annotation']).any()):
    raise ValueError('Expected 76 unique rows in 19 four-annotation families')

summary_columns = list(df.columns)
df['NegLog10P'] = -np.log10(df['Enrichment_p'].clip(lower=1e-300))
df['NegLog10FDR'] = -np.log10(df['FDR_p'].clip(lower=1e-300))
df['Significance'] = df['FDR_p'].apply(fdr_stars).replace('', 'ns')
df['PopTrait'] = df['Population'] + ' ' + df['Trait']

df[summary_columns].to_csv(SUMMARY_PATH, sep='\t', index=False)

# --- annotation sizes: sum per-chromosome M_5_50 ---

m5_50_rows = []
for pop in ['EUR', 'EAS']:
    for annot in ANNOT_ORDER:
        total = 0.0
        for chrom in range(1, 23):
            f = LD_SCORES_DIR / pop / f'{annot}.{chrom}.l2.M_5_50'
            if f.exists():
                total += float(f.read_text().strip())
        m5_50_rows.append({'Population': pop, 'Annotation': annot, 'M_5_50': total})
m5_50 = pd.DataFrame(m5_50_rows)

# --- Figure 1: developmental annotation heatmaps (interaction excluded) ---

fig, axes = plt.subplots(
    1, 2, figsize=(9.6, 5.4),
    gridspec_kw={'left': 0.09, 'right': 0.86, 'top': 0.90,
                 'bottom': 0.14, 'wspace': 0.32})
im = None
for ax, pop in zip(axes, ['EUR', 'EAS']):
    pop_df = df[(df['Population'] == pop) &
                (df['Annotation'].isin(MAIN_ANNOT_ORDER))].copy()
    pivot = pop_df.pivot_table(index='Trait', columns='Annotation',
                               values='Enrichment', aggfunc='first')
    col_order = [c for c in MAIN_ANNOT_ORDER if c in pivot.columns]
    pivot = pivot[col_order]
    fdr_pivot = pop_df.pivot_table(index='Trait', columns='Annotation',
                                   values='FDR_p', aggfunc='first'
                                   ).reindex(index=pivot.index, columns=col_order)
    trait_order = sorted_traits_by_signal(pop_df)
    pivot = pivot.loc[trait_order]
    fdr_pivot = fdr_pivot.loc[trait_order]
    display_columns = [ANNOT_SHORT_NAMES.get(c, c) for c in pivot.columns]
    pivot_clipped = pivot.clip(ENRICHMENT_VMIN, ENRICHMENT_VMAX)
    norm = mcolors.TwoSlopeNorm(vmin=ENRICHMENT_VMIN,
                                vcenter=ENRICHMENT_VCENTER,
                                vmax=ENRICHMENT_VMAX)
    im = ax.imshow(pivot_clipped.values, cmap='RdBu_r', aspect='auto', norm=norm)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            if not np.isnan(val):
                stars = fdr_stars(fdr_pivot.iloc[i, j])
                text_color = 'white' if abs(val - ENRICHMENT_VCENTER) > 2.2 else 'black'
                if val > ENRICHMENT_VMAX:
                    cell = f'>{ENRICHMENT_VMAX:g}'
                elif val < ENRICHMENT_VMIN:
                    cell = f'<{ENRICHMENT_VMIN:g}'
                else:
                    cell = f'{val:.1f}'
                ax.text(j, i, f'{cell}{stars}', ha='center', va='center',
                        fontsize=6, color=text_color)
    ax.set_xticks(range(len(display_columns)))
    ax.set_xticklabels(display_columns, rotation=30, ha='right', fontsize=7)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([TRAIT_NAMES.get(t, t) for t in pivot.index], fontsize=7)
    ax.set_title(POP_NAMES.get(pop, pop), fontsize=9, fontweight='bold')
    add_panel_label(ax, 'A' if pop == 'EUR' else 'B')

cbar_ax = fig.add_axes([0.895, 0.24, 0.014, 0.55])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_ticks([-4, -2, 0, 2, 4])
cbar.set_label('Enrichment (capped; 1 = no enrichment)', fontsize=8)
cbar.ax.tick_params(labelsize=7)
fig.suptitle('APA-QTL heritability enrichment - developmental annotations',
             fontsize=11, fontweight='bold')
save_figure(fig, OUTPUT_DIR / 'sldsc_figure1_development_heatmap.png')

# --- Figure 2: full enrichment forest, one panel per annotation ---

fig, axes = plt.subplots(
    2, 2, figsize=(10.5, 9.2),
    gridspec_kw={'left': 0.20, 'right': 0.97, 'top': 0.93,
                 'bottom': 0.06, 'hspace': 0.42, 'wspace': 0.30})
for ax, annot, plabel in zip(axes.flat, ANNOT_ORDER, ['A', 'B', 'C', 'D']):
    sub = df[df['Annotation'] == annot].copy()
    row_order = signal_row_order(sub)
    sub = sub.set_index('PopTrait').loc[row_order].reset_index()

    y = np.arange(len(sub))
    colors = sub['Population'].map(POP_COLORS)

    ax.errorbar(sub['Enrichment'], y, xerr=sub['Enrichment_std_error'],
                fmt='none', ecolor='0.35', elinewidth=0.9, capsize=2, zorder=1)
    ax.scatter(sub['Enrichment'], y, s=30, c=colors,
               edgecolor='black', linewidth=0.5, zorder=2)

    ax.axvline(1, color='0.5', linestyle='--', linewidth=0.8, zorder=0)
    lo, hi = forest_panel_limits(sub, 'Enrichment', 'Enrichment_std_error', 1.0)
    linthresh = 10 if annot == 'interaction_apa_qtl' else 1
    ax.set_xlim(lo, hi)
    ax.set_xscale('symlog', linthresh=linthresh)

    ax.set_yticks(y)
    ax.set_yticklabels(sub['PopTrait'], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel('Enrichment +/- SE (symlog)')
    ax.set_title(ANNOT_NAMES.get(annot, annot), fontsize=9, fontweight='bold')
    ax.grid(axis='x', color='0.92', linewidth=0.5)
    clean_axes(ax)
    add_panel_label(ax, plabel)

fig.legend(handles=pop_legend_handles(), loc='upper center',
           bbox_to_anchor=(0.5, 0.985), ncol=2, frameon=False)
fig.suptitle('APA-QTL heritability enrichment - all tests',
             fontsize=11, fontweight='bold', y=0.995)
save_figure(fig, OUTPUT_DIR / 'sldsc_figure2_enrichment_forest.png')

# --- Figure 3: full tau-coefficient forest (x 1e-8, 95% CI) ---

fig, axes = plt.subplots(
    2, 2, figsize=(10.5, 9.2),
    gridspec_kw={'left': 0.20, 'right': 0.97, 'top': 0.93,
                 'bottom': 0.06, 'hspace': 0.42, 'wspace': 0.30})
for ax, annot, plabel in zip(axes.flat, ANNOT_ORDER, ['A', 'B', 'C', 'D']):
    sub = df[df['Annotation'] == annot].copy()
    sub['Tau'] = sub['Coefficient'] * 1e8
    sub['Tau_SE'] = sub['Coefficient_std_error'] * 1e8

    row_order = signal_row_order(sub)
    sub = sub.set_index('PopTrait').loc[row_order].reset_index()

    y = np.arange(len(sub))
    colors = sub['Population'].map(POP_COLORS)

    ax.errorbar(sub['Tau'], y, xerr=1.96 * sub['Tau_SE'],
                fmt='none', ecolor='0.35', elinewidth=0.9, capsize=2, zorder=1)
    ax.scatter(sub['Tau'], y, s=30, c=colors,
               edgecolor='black', linewidth=0.5, zorder=2)

    ax.axvline(0, color='0.5', linestyle='--', linewidth=0.8, zorder=0)
    lo, hi = forest_panel_limits(sub, 'Tau', 'Tau_SE', 1.96)
    linthresh = 100 if annot == 'interaction_apa_qtl' else 1
    ax.set_xlim(lo, hi)
    ax.set_xscale('symlog', linthresh=linthresh)

    ax.set_yticks(y)
    ax.set_yticklabels(sub['PopTrait'], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel('Coefficient (tau) x 10-8 (symlog)')
    ax.set_title(ANNOT_NAMES.get(annot, annot), fontsize=9, fontweight='bold')
    ax.grid(axis='x', color='0.92', linewidth=0.5)
    clean_axes(ax)
    add_panel_label(ax, plabel)

fig.legend(handles=pop_legend_handles(), loc='upper center',
           bbox_to_anchor=(0.5, 0.985), ncol=2, frameon=False)
fig.suptitle('APA-QTL stratified coefficient (tau) - all tests',
             fontsize=11, fontweight='bold', y=0.995)
save_figure(fig, OUTPUT_DIR / 'sldsc_figure3_tau_forest.png')

# --- Figure 4: interaction scan, signed -log10(FDR) ---

inter = df[df['Annotation'] == 'interaction_apa_qtl'].copy()
inter = inter.sort_values(['Population', 'FDR_p', 'Trait'])
inter['label'] = inter['Population'] + ' ' + \
    inter['Trait'].map(lambda x: TRAIT_NAMES.get(x, x))
inter['signed_fdr'] = np.where(
    inter['Enrichment'] >= 0, inter['NegLog10FDR'], -inter['NegLog10FDR'])

fig, ax = plt.subplots(figsize=(6.4, 6.8),
                       gridspec_kw={'left': 0.30, 'right': 0.97,
                                    'top': 0.93, 'bottom': 0.07})
y = np.arange(len(inter))
colors = inter['Population'].map(POP_COLORS)

ax.axvline(0, color='0.55', linewidth=0.8)
threshold = -np.log10(0.05)
ax.axvline(threshold, color='0.75', linestyle='--', linewidth=0.8)
ax.axvline(-threshold, color='0.75', linestyle='--', linewidth=0.8)

ax.scatter(inter['signed_fdr'], y, s=32, c=colors,
           edgecolor='black', linewidth=0.4, zorder=2)

for xpos, ypos, fdr in zip(inter['signed_fdr'], y, inter['FDR_p']):
    stars = fdr_stars(fdr)
    if stars:
        dx = 0.015
        ax.text(xpos + dx, ypos, stars, ha='left', va='center',
                fontsize=7, fontweight='bold', zorder=3)

ax.set_yticks(y)
ax.set_yticklabels(inter['label'], fontsize=6.5)
ax.invert_yaxis()
ax.set_xlim(-1.6, 1.6)
ax.set_xticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
ax.set_xlabel('Signed -log10(FDR)')
ax.set_title('Interaction APA-QTL scan', fontsize=10, fontweight='bold')
ax.grid(axis='x', color='0.92', linewidth=0.5)
clean_axes(ax)
ax.legend(handles=pop_legend_handles(), frameon=False, fontsize=7,
          loc='lower right', ncol=1)

save_figure(fig, OUTPUT_DIR / 'sldsc_figure4_interaction_scan.png')

# --- Figure 5: EUR-EAS consistency scatter (shared traits) ---

eur = df[df['Population'] == 'EUR'].copy()
eas = df[df['Population'] == 'EAS'].copy()
shared = sorted(set(eur['Trait']) & set(eas['Trait']))
pair = eur[eur['Trait'].isin(shared)].merge(
    eas[['Trait', 'Annotation', 'Enrichment', 'Enrichment_std_error',
         'FDR_p']],
    on=['Trait', 'Annotation'], suffixes=('_EUR', '_EAS'))

fig, ax = plt.subplots(figsize=(7.0, 6.4),
                       gridspec_kw={'left': 0.13, 'right': 0.95,
                                    'top': 0.93, 'bottom': 0.10})

linthresh = 10
ax.set_xscale('symlog', linthresh=linthresh)
ax.set_yscale('symlog', linthresh=linthresh)

for annot in ANNOT_ORDER:
    sub = pair[pair['Annotation'] == annot]
    ax.errorbar(
        sub['Enrichment_EUR'], sub['Enrichment_EAS'],
        xerr=sub['Enrichment_std_error_EUR'],
        yerr=sub['Enrichment_std_error_EAS'],
        fmt='o', ms=7, color=ANNOT_COLORS[annot], ecolor='0.55',
        elinewidth=0.8, capsize=2, capthick=0.7,
        markeredgecolor='black', markeredgewidth=0.5,
        label=ANNOT_SHORT_NAMES[annot], zorder=2)

lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
ax.plot([0, lim], [0, lim], color='0.5', linestyle='--',
        linewidth=0.8, zorder=1, label='Identity')
ax.axvline(1, color='0.85', linewidth=0.6, zorder=0)
ax.axhline(1, color='0.85', linewidth=0.6, zorder=0)

ax.set_xlabel('EUR enrichment')
ax.set_ylabel('EAS enrichment')
ax.set_title(f'EUR-EAS consistency ({len(shared)} shared traits)',
             fontsize=10, fontweight='bold')
ax.grid(color='0.93', linewidth=0.5, zorder=0)
clean_axes(ax)
ax.legend(frameon=False, fontsize=7, loc='upper left')

save_figure(fig, OUTPUT_DIR / 'sldsc_figure5_population_consistency.png')

# --- Figure 6: stage comparison, postnatal - prenatal (exploratory z-test, BH over tests) ---

post = df[df['Annotation'] == 'postnatal_apa_qtl'].copy()
pre = df[df['Annotation'] == 'prenatal_apa_qtl'].copy()
comp = post[['Population', 'Trait', 'Enrichment', 'Enrichment_std_error']].merge(
    pre[['Population', 'Trait', 'Enrichment', 'Enrichment_std_error']],
    on=['Population', 'Trait'], suffixes=('_post', '_pre'))

comp['Delta'] = comp['Enrichment_post'] - comp['Enrichment_pre']
comp['Delta_SE'] = np.sqrt(comp['Enrichment_std_error_post'] ** 2 +
                           comp['Enrichment_std_error_pre'] ** 2)
comp['Z'] = comp['Delta'] / comp['Delta_SE']
comp['P'] = [math.erfc(abs(z) / math.sqrt(2.0)) for z in comp['Z']]
comp['FDR'] = bh_fdr(comp['P'])

comp = comp.sort_values(['Population', 'Trait']).reset_index(drop=True)
comp['label'] = comp['Population'] + ' ' + comp['Trait']
y = np.arange(len(comp))
colors = comp['Population'].map(POP_COLORS)

fig, ax = plt.subplots(figsize=(6.6, 6.8),
                       gridspec_kw={'left': 0.26, 'right': 0.97,
                                    'top': 0.92, 'bottom': 0.07})

ax.errorbar(comp['Delta'], y, xerr=comp['Delta_SE'],
            fmt='none', ecolor='0.35', elinewidth=0.9, capsize=2, zorder=1)
ax.scatter(comp['Delta'], y, s=34, c=colors,
           edgecolor='black', linewidth=0.5, zorder=2)

for d, yy, fdr in zip(comp['Delta'], y, comp['FDR']):
    stars = fdr_stars(fdr)
    if stars:
        dx = 0.015
        ax.text(d + dx, yy, stars, ha='left', va='center',
                fontsize=7, fontweight='bold', zorder=3)

ax.axvline(0, color='0.5', linestyle='--', linewidth=0.8, zorder=0)
lo, hi = forest_panel_limits(comp, 'Delta', 'Delta_SE', 1.0)
ax.set_xlim(lo, hi)

ax.set_yticks(y)
ax.set_yticklabels(comp['label'], fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel('Postnatal - Prenatal enrichment')
ax.set_title('Stage comparison (postnatal vs prenatal)',
             fontsize=10, fontweight='bold')
ax.grid(axis='x', color='0.92', linewidth=0.5)
clean_axes(ax)
ax.legend(handles=pop_legend_handles(), frameon=False, fontsize=7,
          loc='lower right', ncol=1)

save_figure(fig, OUTPUT_DIR / 'sldsc_figure6_stage_comparison.png')

# --- Figure 7: -log10(FDR) heatmap, EUR and EAS ---

fig, axes = plt.subplots(
    1, 2, figsize=(9.6, 5.6),
    gridspec_kw={'left': 0.09, 'right': 0.86, 'top': 0.90,
                 'bottom': 0.15, 'wspace': 0.32})

for ax, pop in zip(axes, ['EUR', 'EAS']):
    pop_df = df[df['Population'] == pop]
    traits = sorted_traits_by_signal(pop_df)
    matrix = np.full((len(traits), len(ANNOT_ORDER)), np.nan)
    for i, trait in enumerate(traits):
        for j, annot in enumerate(ANNOT_ORDER):
            val = pop_df[(pop_df['Trait'] == trait) &
                         (pop_df['Annotation'] == annot)]['NegLog10FDR']
            if len(val) > 0 and not np.isnan(val.values[0]):
                matrix[i, j] = val.values[0]

    im = ax.imshow(matrix, cmap='viridis', aspect='auto', vmin=0, vmax=3)

    for i in range(len(traits)):
        for j in range(len(ANNOT_ORDER)):
            val = matrix[i, j]
            if np.isnan(val):
                continue
            stars = ''
            if val > -np.log10(0.05):
                stars = '*' if val <= -np.log10(0.01) else \
                    ('**' if val <= -np.log10(0.001) else '***')
            text_color = 'white' if val > 1.5 else 'black'
            ax.text(j, i, f'{val:.1f}{stars}', ha='center', va='center',
                    fontsize=6, color=text_color)

    ax.set_xticks(range(len(ANNOT_ORDER)))
    ax.set_xticklabels([ANNOT_SHORT_NAMES.get(a, a) for a in ANNOT_ORDER],
                       rotation=30, ha='right', fontsize=7)
    ax.set_yticks(range(len(traits)))
    ax.set_yticklabels([TRAIT_NAMES.get(t, t) for t in traits], fontsize=7)
    ax.set_title(POP_NAMES.get(pop, pop), fontsize=9, fontweight='bold')
    add_panel_label(ax, 'A' if pop == 'EUR' else 'B')

cbar_ax = fig.add_axes([0.895, 0.22, 0.014, 0.56])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('-log10(FDR)', fontsize=8)
cbar.ax.tick_params(labelsize=7)
for cutoff in (0.05, 0.01, 0.001):
    cbar.ax.axhline(y=-np.log10(cutoff), color='white',
                    linestyle='--', linewidth=1)

fig.suptitle('APA-QTL enrichment significance (-log10 FDR, * <0.05, '
             '** <0.01, *** <0.001)', fontsize=11, fontweight='bold')
save_figure(fig, OUTPUT_DIR / 'sldsc_figure7_significance_heatmap.png')

# --- Figure 8: annotation size vs mean enrichment SE ---

se_summary = (df.groupby(['Population', 'Annotation'])
              ['Enrichment_std_error'].mean().reset_index()
              .rename(columns={'Enrichment_std_error': 'Mean_SE'}))
qc = se_summary.merge(m5_50, on=['Population', 'Annotation'])

fig, ax = plt.subplots(figsize=(7.0, 5.6),
                       gridspec_kw={'left': 0.14, 'right': 0.96,
                                    'top': 0.90, 'bottom': 0.12})

markers = {'EUR': 'o', 'EAS': 's'}
for annot in ANNOT_ORDER:
    for pop in ['EUR', 'EAS']:
        row = qc[(qc['Annotation'] == annot) & (qc['Population'] == pop)]
        if row.empty:
            continue
        ax.scatter(row['M_5_50'], row['Mean_SE'],
                   s=70, marker=markers[pop], color=ANNOT_COLORS[annot],
                   edgecolor='black', linewidth=0.5, zorder=3,
                   label=None if pop == 'EUR' else f'{pop}')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Annotation size M_5_50 (log scale)')
ax.set_ylabel('Mean enrichment SE across tests (log scale)')
ax.set_title('Annotation size QC', fontsize=10, fontweight='bold')
ax.grid(which='both', color='0.93', linewidth=0.5, zorder=0)
clean_axes(ax)

h1 = [plt.Line2D([0], [0], marker='o', linestyle='none', markersize=6,
                 markerfacecolor=ANNOT_COLORS[a], markeredgecolor='black',
                 label=ANNOT_SHORT_NAMES[a]) for a in ANNOT_ORDER]
h2 = [plt.Line2D([0], [0], marker=m, linestyle='none', markersize=6,
                 markerfacecolor='white', markeredgecolor='black',
                 label=p) for p, m in markers.items()]
legend1 = ax.legend(handles=h1, frameon=False, fontsize=7,
                    loc='lower left', title='Annotation')
legend1.get_title().set_fontsize(7)
ax.add_artist(legend1)
ax.legend(handles=h2, frameon=False, fontsize=7, loc='lower right',
          title='Population')
ax.get_legend().get_title().set_fontsize(7)

save_figure(fig, OUTPUT_DIR / 'sldsc_figure8_annotation_size_qc.png')

# --- per-annotation enrichment bars, EUR and EAS side by side ---

for annot in ANNOT_ORDER:
    annot_df = df[df['Annotation'] == annot].copy()
    if annot_df.empty:
        continue

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax, pop in zip(axes, ['EUR', 'EAS']):
        pop_df = annot_df[annot_df['Population'] == pop]
        traits = sorted(pop_df['Trait'].unique())
        values = [pop_df[pop_df['Trait'] == t]['Enrichment'].values[0]
                  for t in traits]
        errors = [pop_df[pop_df['Trait'] == t]['Enrichment_std_error'].values[0]
                  for t in traits]
        fdr_vals = [pop_df[pop_df['Trait'] == t]['FDR_p'].values[0]
                    for t in traits]

        ax.barh(range(len(traits)), values, xerr=errors,
                color=OKABE_ITO[0] if pop == 'EUR' else OKABE_ITO[1],
                edgecolor='black', capsize=3, linewidth=0.6)

        for yy, v, e, fdr in zip(range(len(traits)), values, errors, fdr_vals):
            if np.isnan(v):
                continue
            stars = fdr_stars(fdr)
            if stars:
                xpos = v + e + 0.03 * max(abs(v), 1.0)
                ax.text(xpos, yy, stars, va='center', fontsize=8,
                        fontweight='bold')

        ax.set_yticks(range(len(traits)))
        ax.set_yticklabels([TRAIT_NAMES.get(t, t) for t in traits], fontsize=7)
        ax.set_xlabel('Enrichment')
        ax.set_title(POP_NAMES.get(pop, pop), fontweight='bold', fontsize=9)
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=0.8)
        ax.invert_yaxis()
        clean_axes(ax)

        max_abs = max(abs(v) for v in values if not np.isnan(v)) if values else 1
        if max_abs > 10:
            ax.set_xscale('symlog', linthresh=1)
            ax.set_xlabel('Enrichment (symlog scale)')

    fig.suptitle(f'Enrichment: {ANNOT_NAMES.get(annot, annot)}',
                 fontsize=11, fontweight='bold')
    save_figure(fig, OUTPUT_DIR / f'enrichment_{annot}.png')

# --- processed table ---

df.to_csv(DATA_DIR / 'sldsc_enrichment_processed.tsv', sep='\t', index=False)
