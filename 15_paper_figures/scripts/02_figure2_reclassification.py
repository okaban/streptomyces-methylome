#!/usr/bin/env python3
"""
Figure 2: 4mC Reclassification
  Panel A: Motif attribution pie charts (4mC, 6mA)
  Panel B: Pisciotta 2023 comparison table
  Panel C: 4mC vs 5mC at CCGG sites (genome-wide)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

import importlib
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def classify_4mC_motifs(df):
    """Group 4mC motifs: CCGG family (TGGCCGGC) + AAGCCCG."""
    df_uniq = df.drop_duplicates(subset=['position']).copy()
    ccgg_family = ['TGGCCGGC', 'GGCCGG', 'CCGG (other)', 'CCGG']
    df_uniq['motif_group'] = df_uniq['final_motif'].apply(
        lambda x: 'CCGG family' if x in ccgg_family else x)
    return df_uniq['motif_group'].value_counts()


def classify_6mA_motifs(df):
    """Group 6mA motifs: main motifs + other + unassigned."""
    df_uniq = df.drop_duplicates(subset=['position']).copy()
    main_motifs = ['AAGCCCG', 'GATC', 'CCGKCA']
    def group(x):
        if x in main_motifs:
            return x
        elif x == 'unassigned':
            return 'Unassigned'
        else:
            return 'Other motifs'
    df_uniq['motif_group'] = df_uniq['final_motif'].apply(group)
    return df_uniq['motif_group'].value_counts()


def make_motif_pie(ax, counts, title, color_map):
    """Create motif attribution donut chart."""
    total = counts.sum()
    labels = []
    sizes = []
    colors = []

    for motif in counts.index:
        n = counts[motif]
        pct = n / total * 100
        labels.append(f'{motif}\n{n:,} ({pct:.1f}%)')
        sizes.append(n)
        colors.append(color_map.get(motif, COL_GRAY))

    wedges, _ = ax.pie(
        sizes, labels=None, colors=colors,
        startangle=90, counterclock=False,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.5, 'width': 0.65},
        pctdistance=0.75,
    )

    ax.legend(wedges, labels, loc='center left', bbox_to_anchor=(0.95, 0.5),
              fontsize=7.5, frameon=False, handlelength=1.0, handletextpad=0.4,
              labelspacing=1.2)
    ax.set_title(title, fontsize=11, fontweight='bold', pad=8)


def make_comparison_table(ax):
    """Create Pisciotta 2023 comparison table as figure element."""
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)

    # Table structure
    rows = [
        ('CCGG',     '5mC\n(BS-seq)',     '4mC\n(Nanopore)',     'Reclassified'),
        ('GGCCGG',   '5mC\n(BS-seq)',     '4mC\n(Nanopore)',     'Reclassified'),
        ('TGGCCGGC', '5mC\n(BS-seq)',     '4mC\n(Nanopore)',     'Reclassified'),
        ('AAGCCCG',  'Not reported',       '4mC + 6mA\n(Nanopore)', 'Novel'),
        ('GATC',     'Not reported',       '6mA\n(Nanopore)',     'Novel'),
    ]
    headers = ['Motif', 'Pisciotta et al.\n2023', 'This study', 'Status']

    # Column positions (x coords)
    col_x = [0, 2.0, 4.8, 7.5]
    col_w = [2.0, 2.8, 2.7, 2.5]
    row_h = 0.95
    header_y = 6.2

    # Draw header
    for j, (hdr, cx, cw) in enumerate(zip(headers, col_x, col_w)):
        rect = plt.Rectangle((cx, header_y), cw, row_h,
                              facecolor='#37474F', edgecolor='none')
        ax.add_patch(rect)
        ax.text(cx + cw/2, header_y + row_h/2, hdr,
                ha='center', va='center', fontsize=8, fontweight='bold',
                color='white', linespacing=1.3)

    # Draw data rows
    for i, (motif, pisciotta, this_study, status) in enumerate(rows):
        y = header_y - (i + 1) * row_h
        bg = '#F5F5F5' if i % 2 == 0 else 'white'

        for j, (val, cx, cw) in enumerate(
                zip([motif, pisciotta, this_study, status], col_x, col_w)):
            rect = plt.Rectangle((cx, y), cw, row_h,
                                  facecolor=bg, edgecolor='#E0E0E0', linewidth=0.5)
            ax.add_patch(rect)

            kwargs = dict(ha='center', va='center', fontsize=7.5, linespacing=1.3)
            if j == 0:  # Motif column: monospace bold
                kwargs['fontfamily'] = 'monospace'
                kwargs['fontweight'] = 'bold'
                kwargs['fontsize'] = 8
            elif j == 3:  # Status column: colored
                if status == 'Reclassified':
                    kwargs['color'] = COL_4mC
                    kwargs['fontweight'] = 'bold'
                elif status == 'Novel':
                    kwargs['color'] = COL_6mA
                    kwargs['fontweight'] = 'bold'

            ax.text(cx + cw/2, y + row_h/2, val, **kwargs)

    # Bracket for CCGG family
    bx = col_x[0] - 0.15
    y_top = header_y - 0.1 * row_h
    y_bot = header_y - 3 * row_h + 0.1 * row_h
    ax.annotate('', xy=(bx, y_bot), xytext=(bx, y_top),
                arrowprops=dict(arrowstyle='-', lw=1.5, color=COL_4mC))
    ax.plot([bx, bx + 0.1], [y_top, y_top], color=COL_4mC, lw=1.5)
    ax.plot([bx, bx + 0.1], [y_bot, y_bot], color=COL_4mC, lw=1.5)
    mid_y = (y_top + y_bot) / 2
    ax.text(bx - 0.1, mid_y, 'CCGG\nfamily', ha='right', va='center',
            fontsize=7, color=COL_4mC, fontweight='bold', linespacing=1.2)


def make_4mc_5mc_comparison(ax, df_summary):
    """Panel C: 4mC vs 5mC genome-wide HC sites per timepoint."""
    df = df_summary.copy()

    grouped = df.groupby('timepoint').agg(
        hc_4mc_mean=('4mC_HC_sites', 'mean'),
        hc_4mc_se=('4mC_HC_sites', 'sem'),
        hc_5mc_mean=('5mC_HC_sites', 'mean'),
        hc_5mc_se=('5mC_HC_sites', 'sem'),
    ).reindex(['T1', 'T2', 'T3'])

    x = np.arange(len(grouped))
    width = 0.35

    ax.bar(x - width/2, grouped['hc_4mc_mean'], width,
           yerr=grouped['hc_4mc_se'], capsize=4, error_kw={'lw': 1.2},
           color=COL_4mC, edgecolor='white', linewidth=0.5,
           label='4mC', alpha=0.9, zorder=3)
    bars_5mc = ax.bar(x + width/2, grouped['hc_5mc_mean'], width,
                       color=COL_GRAY, edgecolor='#9E9E9E', linewidth=0.8,
                       label='5mC', alpha=0.7, zorder=3)

    # Mark "0" on 5mC bars
    for xi in x:
        ax.text(xi + width/2, 30, '0', ha='center', va='bottom',
                fontsize=11, fontweight='bold', color=COL_DARK, zorder=4)

    ax.set_xticks(x)
    ax.set_xticklabels(['T1 (24 h)', 'T2 (48 h)', 'T3 (72 h)'])
    ax.set_ylabel('HC sites (freq. > 50%)')
    ax.set_title('Genome-wide: 4mC vs 5mC', fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, frameon=False, loc='upper left')

    # Annotation box
    ax.text(0.97, 0.92,
            'HC 5mC = 0\nacross all 9 samples',
            transform=ax.transAxes, fontsize=8,
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0',
                      edgecolor=COL_4mC, linewidth=1.2, alpha=0.95))

    ax.set_ylim(bottom=0)
    ax.grid(axis='y', alpha=0.3, lw=0.5)


def main():
    apply_style()
    print('=== Figure 2: 4mC Reclassification ===')

    # Load data
    print('Loading data...')
    df_4mc, df_6ma = load_methylation_census()
    _, df_summary = load_5mc_4mc_data()

    # Classify motifs
    counts_4mc = classify_4mC_motifs(df_4mc)
    counts_6ma = classify_6mA_motifs(df_6ma)

    print(f'  4mC motif groups: {dict(counts_4mc)}')
    print(f'  6mA motif groups: {dict(counts_6ma)}')

    # Color maps
    color_4mc = {
        'CCGG family': '#C62828',
        'AAGCCCG': COL_BOTH,
    }
    color_6ma = {
        'AAGCCCG': COL_BOTH,
        'GATC': '#1565C0',
        'CCGKCA': '#0D47A1',
        'Other motifs': '#78909C',
        'Unassigned': '#CFD8DC',
    }

    # Order for pie charts
    order_4mc = ['CCGG family', 'AAGCCCG']
    counts_4mc = counts_4mc.reindex([m for m in order_4mc if m in counts_4mc.index])
    order_6ma = ['AAGCCCG', 'GATC', 'CCGKCA', 'Other motifs', 'Unassigned']
    counts_6ma = counts_6ma.reindex([m for m in order_6ma if m in counts_6ma.index])

    # Create figure — 2 rows × 3 columns layout
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(160)))
    gs = GridSpec(2, 3, figure=fig, hspace=0.5, wspace=0.35,
                  height_ratios=[1, 1.1], width_ratios=[1, 1, 1])

    # Panel A: Pie charts (top row, left and right)
    ax_pie4 = fig.add_subplot(gs[0, 0:2])  # span 2 columns for 4mC
    ax_pie6 = fig.add_subplot(gs[0, 2])

    # Create sub-axes for two pies within the top row
    # Actually, let's use a different approach: manual positioning
    fig.clf()
    # Top row: two pie charts
    ax_pie4 = fig.add_axes([0.02, 0.55, 0.35, 0.40])
    ax_pie6 = fig.add_axes([0.45, 0.55, 0.35, 0.40])
    # Bottom left: table
    ax_table = fig.add_axes([0.02, 0.03, 0.50, 0.45])
    # Bottom right: bar chart
    ax_bar = fig.add_axes([0.60, 0.08, 0.36, 0.38])

    # Panel A: Pie charts
    print('Panel A: Motif attribution pie charts...')
    make_motif_pie(ax_pie4, counts_4mc, '4mC motif attribution', color_4mc)
    make_motif_pie(ax_pie6, counts_6ma, '6mA motif attribution', color_6ma)
    add_panel_label(ax_pie4, 'a', x=-0.05, y=1.10)

    # Panel B: Comparison table
    print('Panel B: Pisciotta comparison table...')
    make_comparison_table(ax_table)
    add_panel_label(ax_table, 'b', x=-0.02, y=1.02)

    # Panel C: 4mC vs 5mC
    print('Panel C: 4mC vs 5mC at CCGG...')
    make_4mc_5mc_comparison(ax_bar, df_summary)
    add_panel_label(ax_bar, 'c', x=-0.18, y=1.10)

    # Save
    out_path = FIG_DIR / 'Figure2_4mC_reclassification'
    save_figure(fig, out_path)

    print('\n=== Figure 2 complete ===')


if __name__ == '__main__':
    main()
