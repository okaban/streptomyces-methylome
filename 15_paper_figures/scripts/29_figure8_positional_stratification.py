#!/usr/bin/env python3
"""
Figure 8: Positional stratification of methylation sites and KEGG/COG signal.

Reviewer A2 + A4 fixes:
  - Drop the "Unassigned 6mA" category in panel A (no biological motif).
  - Show only the two assigned motifs: GCCGGC (4mC) and AAGCCCG (4mC/6mA).
  - Panel A: stacked bar of genomic position categories per motif (two motifs).
  - Panel B: position-stratified KEGG/COG enrichment heatmap (significant only).
  - Publication-quality formatting (Arial 9 pt, dual PDF + SVG).
"""

from pathlib import Path
import importlib
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


KEGG_DIR = EPIGENOME / '62_GO_KEGG_enrichment' / 'tables'

CATEGORY_ORDER = ['promoter', '5UTR_approx', 'CDS_internal', 'intergenic']
# Short legend labels (full definitions are in the caption); the previous
# spelled-out labels with parentheticals overran the figure's right edge and
# were clipped, because save_figure() honours the declared NAR width and does
# not tight-crop (figure-legibility-qc §3).
CATEGORY_LABELS = ['Promoter', "5′ UTR", 'CDS', 'Intergenic']
# Genomic-category qualitative set — muted tones threaded from the unified palette
CATEGORY_COLORS = [COL_ARM, COL_CORE, COL_6mA, COL_GRAY]  # amber, green, blue, grey

MOTIFS = [
    ('GCCGGC', 'GCCGGC (4mC)', COL_4mC),
    ('AAGCCCG', 'AAGCCCG (4mC/6mA)', COL_BOTH),
]


def panel_a(ax):
    """Stacked bar: position categories per motif (Unassigned 6mA dropped)."""
    # Unified region definition (experimental TSS, 300 bp promoter, precedence
    # promoter > 5'-UTR > CDS > intergenic, painted per base so the categories
    # are mutually exclusive). See region_scheme_reconciliation.csv.
    df = pd.read_csv(KEGG_DIR / 'F1_classified_meth_sites_unified.tsv', sep='\t')
    print(f'  classified sites: {len(df)} '
          f'(motifs: {df["motif"].value_counts().to_dict()})')

    counts = (df.groupby(['motif', 'category']).size()
              .unstack(fill_value=0)
              .reindex(index=[m[0] for m in MOTIFS],
                       columns=CATEGORY_ORDER, fill_value=0))
    totals = counts.sum(axis=1)
    pct = counts.div(totals, axis=0) * 100

    y_pos = np.arange(len(MOTIFS))
    height = 0.6
    left = np.zeros(len(MOTIFS))
    for c_idx, cat in enumerate(CATEGORY_ORDER):
        widths = pct[cat].values
        ax.barh(y_pos, widths, height=height, left=left,
                color=CATEGORY_COLORS[c_idx], edgecolor='white',
                linewidth=0.6, label=CATEGORY_LABELS[c_idx],
                zorder=3)
        # In-bar % labels (only if segment > 5 %)
        for i, w in enumerate(widths):
            if w >= 5:
                ax.text(left[i] + w / 2, y_pos[i],
                        f'{w:.0f}%', ha='center', va='center',
                        fontsize=7, color='white', fontweight='bold')
        left += widths

    # Y labels: motif name + sample size
    y_labels = [f'{lbl}\nn = {totals.loc[m]:,}' for m, lbl, _ in MOTIFS]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel('Methylation sites by genomic category (%)', fontsize=9)
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Vertical legend to the right of the bars — clear of the title.
    ax.legend(loc='center left', bbox_to_anchor=(1.01, 0.5), ncol=1,
              fontsize=7, frameon=False, handlelength=1.0, handletextpad=0.4,
              labelspacing=0.7)
    ax.set_title('Site distribution by genomic category (T1)',
                 loc='left', fontsize=10, fontweight='bold', pad=10)


def panel_b(ax):
    """Position-stratified KEGG enrichment heatmap (significant entries only)."""
    df = pd.read_csv(KEGG_DIR / 'F1_stratified_KEGG_by_position_v2_unified.tsv',
                     sep='\t')

    # Numeric-coerce padj/OR
    for c in ('OR', 'padj'):
        df[c] = pd.to_numeric(df[c], errors='coerce')

    # Filter to motifs we keep, drop placeholder rows
    keep = (df['motif'].isin(['GCCGGC 4mC', 'AAGCCCG 6mA']) &
            df['padj'].notna() &
            df['pathway'].notna() &
            ~df['pathway'].astype(str).str.startswith('('))
    df = df[keep].copy()
    df['pathway_short'] = df['pathway'].apply(_short_kegg_name)

    # Pick top pathways by minimum padj (per category) for visual clarity
    top = (df.sort_values('padj')
             .drop_duplicates(['motif', 'category', 'pathway_short'])
             .groupby(['motif', 'category'])
             .head(3))
    pathways = (top.groupby('pathway_short')['padj'].min()
                .sort_values().index.tolist()[:14])
    df = df[df['pathway_short'].isin(pathways)].copy()
    print(f'  KEGG pathways shown: {len(pathways)}')

    # X: (motif, category) tuples; Y: pathway
    motifs_short = ['GCCGGC 4mC', 'AAGCCCG 6mA']
    cats = ['promoter', '5UTR_approx', 'CDS_internal']
    cat_lbl = ['Prom.', "5′UTR", 'CDS']
    x_keys = [(m, c) for m in motifs_short for c in cats]
    x_to_idx = {k: i for i, k in enumerate(x_keys)}
    y_to_idx = {p: i for i, p in enumerate(pathways)}

    # Bubble: size = -log10(padj), color = log2(OR)
    ors = df['OR'].clip(lower=0.05, upper=20)
    log_or = np.log2(ors)

    cmap = LinearSegmentedColormap.from_list(
        'or_cmap', [COL_6mA, '#FFFFFF', COL_4mC])  # unified blue->white->red
    vmax = 4  # log2 OR clipped

    for _, row in df.iterrows():
        if (row['motif'], row['category']) not in x_to_idx:
            continue
        x = x_to_idx[(row['motif'], row['category'])]
        y = y_to_idx[row['pathway_short']]
        size = max(20, -np.log10(max(row['padj'], 1e-6)) * 110)
        c_val = np.clip(np.log2(max(row['OR'], 0.05)), -vmax, vmax)
        sig_edge = '#222' if row['padj'] < 0.10 else '#bbb'
        sig_lw = 1.0 if row['padj'] < 0.10 else 0.3
        ax.scatter(x, y, s=size,
                   c=[cmap((c_val + vmax) / (2 * vmax))],
                   edgecolors=sig_edge, linewidths=sig_lw, zorder=3)
        if row['padj'] < 0.05:
            ax.scatter(x, y, s=10, marker='*', color='black',
                       edgecolors='none', zorder=4)

    # Vertical separator between motifs
    ax.axvline(len(cats) - 0.5, color='#333', linewidth=0.6,
               linestyle=':', zorder=1)

    ax.set_xticks(range(len(x_keys)))
    ax.set_xticklabels(cat_lbl * len(motifs_short), fontsize=7.5)
    # Motif group labels above x ticks (display: 4mC → 4mC; data keys unchanged)
    for i, m in enumerate(motifs_short):
        m_disp = m.replace('GCCGGC 4mC', 'GCCGGC 4mC').replace('AAGCCCG 6mA', 'AAGCCCG 6mA')
        ax.text(i * len(cats) + (len(cats) - 1) / 2, len(pathways) + 0.2,
                m_disp, ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    ax.set_yticks(range(len(pathways)))
    ax.set_yticklabels(pathways, fontsize=7.5)
    ax.set_xlim(-0.6, len(x_keys) - 0.4)
    ax.set_ylim(-0.6, len(pathways) - 0.4 + 0.5)
    ax.invert_yaxis()
    ax.grid(True, axis='both', linestyle=':', linewidth=0.4,
            color=COL_GRID, alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_title('Position-stratified KEGG enrichment',
                 loc='left', fontsize=10, fontweight='bold', pad=18)

    # Color bar (log2 OR)
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.04, pad=0.08, aspect=30)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label(r'log$_2$ Odds Ratio', fontsize=8)


def _short_kegg_name(name: str) -> str:
    import textwrap
    name = str(name).split(' - ')[0]
    name = name.replace(' biosynthesis', ' biosynth.')
    name = name.replace(' metabolism', ' metab.')
    # Wrap long names onto two lines instead of truncating with '…', so no
    # pathway name is lost (e.g. 'Biosynthesis of various nucleotide sugars').
    if len(name) > 30:
        name = '\n'.join(textwrap.wrap(name, width=30, break_long_words=False))
    return name


def main():
    apply_style()
    print('=== Figure 8: Positional stratification (publication) ===')

    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(165)))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.35, 1.0],
                          hspace=0.55, left=0.30, right=0.86,
                          top=0.93, bottom=0.07)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])

    print('Drawing Panel A...')
    panel_a(ax_a)
    print('Drawing Panel B...')
    panel_b(ax_b)

    # lowercase panel letters — uniform offset-points placement (shared helper)
    for ax, L in ((ax_a, 'a'), (ax_b, 'b')):
        _utils.add_panel_label(ax, L)

    # Figure caption note
    fig.text(0.30, 0.02,
             'Bubble size: −log10(FDR);  *: FDR < 0.05;  '
             'thick border: FDR < 0.10.',
             fontsize=7, color='#444')

    out = FIG_DIR / 'Figure8_positional_stratification'
    save_figure(fig, out, formats=('pdf', 'svg', 'png'))
    import shutil
    from pathlib import Path as _P
    slot = _P.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure8.png'
    if slot.parent.is_dir():
        shutil.copyfile(out.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
