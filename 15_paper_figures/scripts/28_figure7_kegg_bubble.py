#!/usr/bin/env python3
"""
Figure 7: KEGG pathway enrichment by methylation motif (publication quality).

Reviewer A3 fix:
  - Combine GCCGGC (4mC), AAGCCCG (4mC/6mA), and Dual-targeted gene sets in a
    single bubble plot.
  - Sort pathways by mean significance, color by motif, size by gene count,
    border by FDR threshold.
  - Use Arial 9 pt, light-grid panel, dual-format (PDF + SVG) output.
"""

from pathlib import Path
import importlib
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)


KEGG_DIR = EPIGENOME / '62_GO_KEGG_enrichment' / 'tables'

# Muted publication palette (consistent with COL_4mC/COL_6mA used elsewhere);
# newline-wrapped labels avoid x-axis tick collisions.
MOTIF_CFG = [
    ('GCCGGC\n(4mC)', 'GCCGGC-proximal_4mC',
     COL_4mC, 'E1_KEGG_enrichment_GCCGGC-proximal_4mC.tsv'),
    ('AAGCCCG\n(4mC/6mA)', 'AAGCCCG-proximal_6mA',
     COL_BOTH, 'E1_KEGG_enrichment_AAGCCCG-proximal_6mA.tsv'),
    # 2026-09-22: was COL_6mA, but blue denotes the 6mA mark elsewhere and this
    # column is a GENE SET (near both motifs), not a modification. Neutral dark.
    ('Dual-\ntargeted', 'Dual-targeted',
     COL_DUAL_SET, 'E1_KEGG_enrichment_Dual-targeted.tsv'),
]

FDR_CUT = 0.10
TOP_PATHWAYS = 12


def load_kegg_long():
    rows = []
    for motif_label, key, color, fname in MOTIF_CFG:
        df = pd.read_csv(KEGG_DIR / fname, sep='\t')
        df = df.assign(motif=motif_label, motif_color=color)
        rows.append(df)
    return pd.concat(rows, ignore_index=True)


def short_pathway_name(name: str) -> str:
    """Shorten KEGG names for plot labels.

    Names are wrapped to at most two lines (no truncation) so every
    pathway label is fully visible on the y-axis.
    """
    import textwrap
    name = name.split(' - ')[0]
    name = name.replace(' biosynthesis', ' biosynth.')
    name = name.replace(' metabolism', ' metab.')
    name = name.replace('Streptomyces coelicolor', 'S. coelicolor')
    # Wrap long names onto two lines instead of truncating with an ellipsis.
    if len(name) > 24:
        name = '\n'.join(textwrap.wrap(name, width=24, break_long_words=False))
    return name


def select_top_pathways(df_long: pd.DataFrame, k: int) -> list:
    """Pick the top-k pathways by best minimum padj across motifs."""
    grp = df_long.groupby('pathway_id').agg(
        best_padj=('padj', 'min'),
        max_n=('n_query', 'max'),
        name=('pathway_name', 'first'),
    )
    sig = grp[grp['best_padj'] < 0.5].sort_values(
        ['best_padj', 'max_n'], ascending=[True, False])
    return sig.index[:k].tolist()


def draw_bubble(ax, df_long: pd.DataFrame, pids: list):
    motif_order = [m[0] for m in MOTIF_CFG]
    motif_colors = {m[0]: m[2] for m in MOTIF_CFG}

    # Build matrix: pathways × motif
    df = df_long[df_long['pathway_id'].isin(pids)].copy()
    df['neg_log_padj'] = -np.log10(df['padj'].clip(1e-6))

    # Sort y-axis by best padj (top significant first → bottom of axis)
    pid_order = pids[::-1]
    pid_to_y = {pid: i for i, pid in enumerate(pid_order)}
    motif_to_x = {m: i for i, m in enumerate(motif_order)}

    # Size scale (gene count): area proportional
    size_min, size_max = 60, 500
    n_min, n_max = max(1, df['n_query'].min()), df['n_query'].max()

    def size(n):
        if n_max == n_min:
            return (size_min + size_max) / 2
        return size_min + (size_max - size_min) * (n - n_min) / (n_max - n_min)

    for _, row in df.iterrows():
        x = motif_to_x[row['motif']]
        y = pid_to_y[row['pathway_id']]
        sig = row['padj'] < FDR_CUT
        ec = '#222' if sig else '#bbb'
        lw = 1.2 if sig else 0.4
        ax.scatter(x, y, s=size(row['n_query']),
                   c=motif_colors[row['motif']],
                   alpha=0.85, edgecolors=ec, linewidths=lw,
                   zorder=3)
        # Star marker for FDR < 0.05
        if row['padj'] < 0.05:
            ax.scatter(x, y, s=12, marker='*', color='white',
                       edgecolors='none', zorder=4)

    # Axes
    ax.set_xticks(list(motif_to_x.values()))
    # Two-line wrapped labels with extra top padding; widened x-limits below
    # spread the three columns so the category labels never collide.
    ax.set_xticklabels(motif_order, fontsize=9, linespacing=1.3)
    ax.set_yticks(list(pid_to_y.values()))
    name_lookup = (df_long.drop_duplicates('pathway_id')
                   .set_index('pathway_id')['pathway_name'])
    ax.set_yticklabels([short_pathway_name(name_lookup[pid]) for pid in pid_order],
                       fontsize=8)
    ax.set_xlim(-0.7, len(motif_order) - 0.3)
    ax.set_ylim(-0.7, len(pid_order) - 0.3)
    ax.grid(True, axis='both', linestyle=':', linewidth=0.5,
            color=COL_GRID, alpha=0.9, zorder=0)
    ax.set_axisbelow(True)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='both', length=3, pad=2)
    ax.set_xlabel('Methylation motif', fontsize=9)


def draw_legends(fig, df_long: pd.DataFrame):
    n_min, n_max = max(1, df_long['n_query'].min()), df_long['n_query'].max()
    sample_ns = sorted({n_min, int((n_min + n_max) / 2), n_max})
    if len(sample_ns) == 1:
        sample_ns = [n_min]

    # Size legend (bottom right)
    size_handles = []
    for n in sample_ns:
        if n_max == n_min:
            s = 200
        else:
            s = 60 + (500 - 60) * (n - n_min) / (n_max - n_min)
        h = Line2D([0], [0], marker='o', linestyle='none',
                   markerfacecolor='#888', markeredgecolor='#222',
                   markersize=np.sqrt(s), label=f'n = {n}')
        size_handles.append(h)
    # Size legend (lower right). Placed in the bottom third of the right band;
    # generous labelspacing keeps the large n bubbles from touching.
    # labelspacing alone measures text rows; the large bubbles (n up to ~500 pt²,
    # ~22 pt across) overflow a text-sized row and touch. handleheight forces each
    # legend row to be tall enough (in font-size units) to hold the biggest marker,
    # so the n = 33 and n = 65 bubbles no longer overlap.
    fig.legend(size_handles, [h.get_label() for h in size_handles],
               loc='upper left', bbox_to_anchor=(0.78, 0.44),
               title='Genes in pathway', title_fontsize=8, fontsize=7,
               frameon=False, labelspacing=2.0, handleheight=3.0,
               handletextpad=1.2, borderaxespad=0)

    # Significance legend
    # ordered by decreasing significance (most significant first)
    fdr_handles = [
        Line2D([0], [0], marker='*', linestyle='none',
               markerfacecolor='white', markeredgecolor='#222',
               markeredgewidth=0.6, markersize=10,
               label='FDR < 0.05 (*)'),
        Line2D([0], [0], marker='o', linestyle='none',
               markerfacecolor='#888', markeredgecolor='#222',
               markeredgewidth=1.2, markersize=10,
               label=f'FDR < {FDR_CUT}'),
        Line2D([0], [0], marker='o', linestyle='none',
               markerfacecolor='#888', markeredgecolor='#bbb',
               markeredgewidth=0.4, markersize=10,
               label=f'FDR ≥ {FDR_CUT}'),
    ]
    # Significance legend (upper right), kept in the top portion of the right
    # band so it never overlaps the 'Genes in pathway' legend below it.
    fig.legend(fdr_handles, [h.get_label() for h in fdr_handles],
               loc='upper left', bbox_to_anchor=(0.78, 0.88),
               title='Significance', title_fontsize=8, fontsize=7,
               frameon=False, labelspacing=0.9, borderaxespad=0)


def main():
    apply_style()
    print('=== Figure 7: KEGG pathway enrichment bubble (publication) ===')
    df_long = load_kegg_long()
    print(f'  Rows loaded: {len(df_long)}')

    pids = select_top_pathways(df_long, TOP_PATHWAYS)
    print(f'  Selected top-{len(pids)} pathways:')
    for pid in pids:
        sub = df_long[df_long['pathway_id'] == pid]
        best = sub['padj'].min()
        print(f'    {pid:10s} padj_min={best:.3g}  '
              f'name={sub["pathway_name"].iloc[0][:60]}')

    # Layout: bubble + 2 separate legend boxes on the right margin.
    # Wider canvas + larger left margin so full (2-line) pathway names fit and
    # the three motif columns are well separated; right band holds the legends.
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(120)))  # NAR full-width
    ax = fig.add_axes([0.30, 0.20, 0.46, 0.70])

    draw_bubble(ax, df_long, pids)
    ax.set_title('KEGG pathway enrichment by methylation motif',
                 fontsize=10, fontweight='bold', pad=8)

    draw_legends(fig, df_long)

    # Save (PNG added alongside PDF/SVG for the manuscript image slot)
    out = FIG_DIR / 'Figure7_kegg_bubble'
    assert_no_text_collisions(fig, 'Figure7_kegg')
    save_figure(fig, out, formats=('pdf', 'svg', 'png'))

    # sync into the Obsidian manuscript slot (Figure7.png) — script previously
    # wrote FIG_DIR only, leaving the slot stale.
    import shutil
    from pathlib import Path as _P
    slot = _P.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure7.png'
    if slot.parent.is_dir():
        shutil.copyfile(out.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')
    print('=== Done ===')


if __name__ == '__main__':
    main()
