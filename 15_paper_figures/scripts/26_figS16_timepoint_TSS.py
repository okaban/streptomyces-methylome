# SUPERSEDED 2026-10-15: このスクリプトは figure_manifest.tsv に SuppFigure16.png の産出元として
# 記録されていないのに SuppFigure16.png を書き出す（現行の出荷図を上書きする）。
# 同型の事故を 2026-10-15 に 1 件起こしている（候補図が既存の補足図を上書きし、
# md5 照合のうえバックアップから復元した）。実行を止めてある。
# 内容を再利用する場合は figure_manifest.tsv に産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 26_figS16_timepoint_TSS.py は産出元として登録されていないため実行を止めています。"
    "SuppFigure16.png の現行の産出元は figure_manifest.tsv を参照してください。"
)

"""SUPERSEDED 2026-09-23 — DO NOT SHIP THIS FIGURE.

This script plots the regulatory-gene protection ratio over time from
64_timepoint_TSS/tables/C2_protection_zone_metrics.tsv, which was built on the
FIRST-APPEARANCE site set and gives 0.839 / 0.844 / 1.046 — the series the
manuscript retired on 2026-09-23. On the canonical per-timepoint site set the
ratio is 0.743 / 0.709 / 0.684 (protection maintained at T3);
see 90_per_timepoint_census_audit/tables/protection_zone_ratio_by_timepoint.tsv.

Its output is referenced by nothing: not the manuscript, not figure_manifest.tsv.
The shipped SuppFigure16.png is the C4-corrected crosstalk figure from
78_reviewer_robustness/scripts/P1_figures.py, which is unrelated to this script
despite the similar name. Rebuild from the 90_ table before reusing any of this.
"""
"""
FigS16: Timepoint-resolved TSS methylation protection zone (C-2)

Panels:
  a. 4mC metagene profiles: All TSS, T1/T2/T3 overlay
  b. 6mA metagene profiles: All TSS, T1/T2/T3 overlay
  c. Protection ratio across timepoints (bar chart, mod_type × regulatory/all)
  d. 4mC Regulatory vs Non-regulatory by timepoint (T1/T2/T3 panels)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from shared_utils_local import apply_style, add_panel_label, save_figure, FIG_SUP_DIR
from importlib import import_module as _im
mm_to_inch = _im('00_shared_utils').mm_to_inch

TSS_DIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/64_timepoint_TSS')
TABLES = TSS_DIR / 'tables'

apply_style()

# ── Load data ────────────────────────────────────────────────────────────────
profiles = pd.read_csv(TABLES / 'C2_all_profiles.tsv', sep='\t')
metrics = pd.read_csv(TABLES / 'C2_protection_zone_metrics.tsv', sep='\t')

TP_COLORS = {'T1': '#4C72B0', 'T2': '#DD8452', 'T3': '#55A868'}
TP_LABELS = {'T1': 'T1 (12 h)', 'T2': 'T2 (24 h)', 'T3': 'T3 (50 h)'}
GENE_CAT_COLORS = {'Regulatory': '#d62728', 'Non-regulatory': '#1f77b4'}

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(200)))  # NAR full-width
gs = fig.add_gridspec(3, 3, hspace=0.60, wspace=0.32,
                      height_ratios=[2.2, 2.2, 2.5])

def plot_profile_overlay(ax, mod, gene_cat, panel_label):
    """Overlay T1/T2/T3 metagene profiles for one mod × gene_cat combination."""
    for tp in ['T1', 'T2', 'T3']:
        sub = profiles[(profiles['mod_type'] == mod) &
                       (profiles['gene_cat'] == gene_cat) &
                       (profiles['timepoint'] == tp)]
        if len(sub) == 0:
            continue
        x = sub['bin_center'].values / 1000
        y = sub['density'].values
        ax.plot(x, y, color=TP_COLORS[tp], linewidth=1.8, label=TP_LABELS[tp])
        ax.fill_between(x, y, alpha=0.10, color=TP_COLORS[tp])
    ax.axvline(0, color='red', linestyle='--', linewidth=0.9, alpha=0.7)
    ax.set_xlim(-5, 5)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=6)
    ax.set_ylabel('Density (sites/kb/TSS)', fontsize=6)
    ax.set_title(f'{mod} — {gene_cat}', fontsize=7)
    ax.legend(fontsize=6, loc='upper right')
    add_panel_label(ax, panel_label)

# Panels a, b: All TSS profiles
ax_a = fig.add_subplot(gs[0, :2])
plot_profile_overlay(ax_a, '4mC', 'All', 'a')
ax_a.set_title('4mC methylation profile ±5 kb from TSS (all primary TSS, n=2,771)', fontsize=7)

ax_b = fig.add_subplot(gs[1, :2])
plot_profile_overlay(ax_b, '6mA', 'All', 'b')
ax_b.set_title('6mA methylation profile ±5 kb from TSS (all primary TSS, n=2,771)', fontsize=7)

# Panel c: Protection ratio bar chart
ax_c = fig.add_subplot(gs[0, 2])
tp_order = ['T1', 'T2', 'T3']
x = np.arange(len(tp_order))
width = 0.2
gene_cats = ['All', 'Regulatory', 'Non-regulatory']
gc_colors = {'All': '#78909C', 'Regulatory': '#d62728', 'Non-regulatory': '#1f77b4'}
gc_offsets = [-0.22, 0, 0.22]

for mod, linestyle in [('4mC', '-'), ('6mA', '--')]:
    for i, gc in enumerate(gene_cats):
        sub = metrics[(metrics['mod_type'] == mod) & (metrics['gene_cat'] == gc)]
        sub = sub.set_index('timepoint')
        vals = [sub.loc[tp, 'protection_ratio'] if tp in sub.index else np.nan
                for tp in tp_order]
        ax_c.plot(x + gc_offsets[i], vals, marker='o', linewidth=1.5,
                  linestyle=linestyle,
                  color=gc_colors[gc],
                  label=f'{mod} {gc}' if mod == '4mC' else None,
                  alpha=0.9, markersize=5)

ax_c.axhline(1.0, color='gray', linestyle=':', linewidth=0.8)
ax_c.set_xticks(x)
ax_c.set_xticklabels(['T1', 'T2', 'T3'], fontsize=6)
ax_c.set_ylabel('Protection ratio', fontsize=7)
ax_c.set_title('Protection ratio\nby gene category', fontsize=7)
ax_c.set_ylim(0.65, 1.12)
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color=gc_colors['All'], linewidth=1.5, label='All'),
    Line2D([0], [0], color=gc_colors['Regulatory'], linewidth=1.5, label='Regulatory'),
    Line2D([0], [0], color=gc_colors['Non-regulatory'], linewidth=1.5, label='Non-reg.'),
    Line2D([0], [0], color='gray', linewidth=1.5, linestyle='-', label='4mC'),
    Line2D([0], [0], color='gray', linewidth=1.5, linestyle='--', label='6mA'),
]
ax_c.legend(handles=legend_elements, fontsize=6, loc='lower right', framealpha=0.8)
ax_c.text(1.8, 1.03, '↑ Collapse', color='gray', fontsize=6, ha='right')
ax_c.text(1.8, 0.97, '↓ Protected', color='gray', fontsize=6, ha='right')
add_panel_label(ax_c, 'c')

# Panel d: 4mC Reg vs Non-reg by timepoint
for col, tp in enumerate(['T1', 'T2', 'T3']):
    ax = fig.add_subplot(gs[2, col])
    for gc, color in [('Regulatory', '#d62728'), ('Non-regulatory', '#1f77b4')]:
        sub = profiles[(profiles['mod_type'] == '4mC') &
                       (profiles['gene_cat'] == gc) &
                       (profiles['timepoint'] == tp)]
        if len(sub) == 0:
            continue
        x_vals = sub['bin_center'].values / 1000
        y_vals = sub['density'].values
        ax.plot(x_vals, y_vals, color=color, linewidth=1.8, label=gc)
        ax.fill_between(x_vals, y_vals, alpha=0.10, color=color)
    ax.axvline(0, color='red', linestyle='--', linewidth=0.9, alpha=0.7)
    ax.set_xlim(-5, 5)
    ax.set_xlabel('Distance from TSS (kb)', fontsize=6)
    ax.set_ylabel('Density (sites/kb/TSS)' if col == 0 else '', fontsize=6)
    ax.set_title(f'4mC — {TP_LABELS[tp]}', fontsize=7)
    if col == 0:
        ax.legend(fontsize=6)
    # Annotate protection ratio
    m_sub = metrics[(metrics['mod_type'] == '4mC') &
                    (metrics['timepoint'] == tp) &
                    (metrics['gene_cat'] == 'Regulatory')]
    if len(m_sub) > 0:
        ratio = m_sub.iloc[0]['protection_ratio']
        color = '#d62728' if ratio > 1.0 else '#555555'
        ax.text(0.98, 0.97, f'Reg. ratio = {ratio:.3f}', transform=ax.transAxes,
                ha='right', va='top', fontsize=6, color=color,
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8))
    if col == 0:
        add_panel_label(ax, 'd')

fig.suptitle('Figure S16: Timepoint-resolved TSS methylation protection dynamics',
             fontsize=9, fontweight='bold', y=1.01)

save_figure(fig, FIG_SUP_DIR / 'FigS16_timepoint_TSS_protection.pdf')
