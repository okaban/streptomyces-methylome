# SUPERSEDED 2026-10-15: 旧世代の図生成スクリプト。現行の出荷図は
# panels_261013/ と supp_261014/ が産出し、figure_manifest.tsv に登録してある。
# 本スクリプトは撤回済みの入力（Exposed/Shielded 区分、退避版テーブル）を参照するため
# 実行を止めてある。再利用する場合は現行サイト集合で再計算し、産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 21_figS_bgc_methylation_track.py は旧世代のため実行を止めています。"
)

#!/usr/bin/env python3
"""
FigS_BGC_Track: BGC Cluster Methylation Coverage Tracks

Generates methylation site maps for all 4 major BGC clusters:
  - ACT (Actinorhodin)
  - CDA (Calcium-dependent antibiotic)
  - CPK (Coelimycin P1)
  - RED (Undecylprodigiosin)

Style: genomic coverage track (cf. coverage_track_act.pdf reference)
  - T1/T2/T3 panels: methylation frequency as lollipop/stem plots
  - Gene annotation track: arrows colored by functional role
"""

import sys, importlib
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
_utils = importlib.import_module('00_shared_utils')
for _attr in dir(_utils):
    if not _attr.startswith('_'):
        globals()[_attr] = getattr(_utils, _attr)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BGC_DIR    = EPIGENOME / '47_BGC_methylation_geographic_test' / 'tables'
CENSUS_4mC = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census_v1.csv'
CENSUS_6mA = EPIGENOME / '23_expanded_motif_search' / '6mA_final_census_v1.csv'

FLANK = 2000   # flanking bp to include beyond cluster boundaries

TP_COLORS = {'T1': '#43A047', 'T2': '#FB8C00', 'T3': '#1E88E5'}

MOD_COLORS = {
    'TGGCCGGC':      COL_4mC,
    'GGCCGG':        COL_4mC,
    'CCGG (other)':  COL_4mC,
    'AAGCCCG':       COL_BOTH,
    '6mA':           COL_6mA,
}

FUNC_COLORS = {
    'biosynthesis': '#42A5F5',
    'regulator':    '#E53935',
    'resistance':   '#FB8C00',
    'transport':    '#43A047',
}

BGC_META = {
    'act': 'ACT (Actinorhodin)',
    'cda': 'CDA (Calcium-dependent antibiotic)',
    'cpk': 'CPK (Coelimycin P1)',
    'red': 'RED (Undecylprodigiosin)',
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_all_methylation():
    c4 = pd.read_csv(CENSUS_4mC)
    c6 = pd.read_csv(CENSUS_6mA)
    c6['final_motif'] = '6mA'
    return pd.concat([c4, c6], ignore_index=True)


def load_all_genes():
    return pd.read_csv(BGC_DIR / 'BGC_gene_geography.tsv', sep='\t')


def get_cluster_sites(df_all, start, end):
    return df_all[(df_all['position'] >= start) & (df_all['position'] <= end)].copy()


def gene_short_label(product):
    skip = {'protein', 'family', 'domain', 'containing', 'class', 'related',
            'dependent', 'binding', 'superfamily', 'subunit', 'component',
            'putative', 'predicted', 'hypothetical'}
    words = product.split()
    tokens = [w for w in words if w.lower() not in skip and len(w) > 2]
    return ' '.join(tokens[:3]) if tokens else product[:15]


# ── Drawing functions ─────────────────────────────────────────────────────────

def draw_methylation_track(ax, df_sites, timepoint, tp_color, xlim):
    sub = df_sites[df_sites['timepoint'] == timepoint]

    if len(sub) == 0:
        ax.text(0.5, 0.5, f'No HC sites at {timepoint}',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=8, color='#888', style='italic')
    else:
        for _, row in sub.iterrows():
            color = MOD_COLORS.get(row['final_motif'], '#9E9E9E')
            ax.plot([row['position'], row['position']], [0, row['frequency']],
                    color=color, linewidth=1.4, alpha=0.85)
            ax.scatter(row['position'], row['frequency'],
                       color=color, s=22, zorder=4, alpha=0.95,
                       edgecolors='white', linewidths=0.4)

    ax.set_xlim(*xlim)
    ax.set_ylim(-5, 108)
    ax.set_yticks([0, 50, 100])
    ax.set_yticklabels(['0', '50', '100'], fontsize=6.5)
    ax.set_ylabel('Mod. freq. (%)', fontsize=7, labelpad=2)
    ax.axhline(0, color='#bbb', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xticklabels([])
    ax.tick_params(bottom=False)

    ax.text(0.012, 0.90, timepoint, transform=ax.transAxes,
            fontsize=9, fontweight='bold', color='white', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.25', facecolor=tp_color,
                      edgecolor='none', alpha=0.9))

    return len(sub)


def draw_gene_track(ax, df_genes, xlim):
    span = xlim[1] - xlim[0]

    ax.set_xlim(*xlim)
    ax.set_ylim(-3.2, 2.8)
    ax.axis('off')

    ARROW_Y = {'+': 0.9, '-': -0.9}
    ARROW_H = 0.70
    LABEL_Y_PLUS  =  1.85
    LABEL_Y_MINUS = -1.85

    ax.text(xlim[0] - span * 0.01, ARROW_Y['+'], '(+)',
            fontsize=6.5, color='#555', va='center', ha='right')
    ax.text(xlim[0] - span * 0.01, ARROW_Y['-'], '(−)',
            fontsize=6.5, color='#555', va='center', ha='right')
    ax.axhline(0, color='#ccc', linewidth=0.7, linestyle='--')

    prev_end_plus  = xlim[0]
    prev_end_minus = xlim[0]
    min_gap = span * 0.05   # minimum gap between labels (5% of span)

    for _, row in df_genes.iterrows():
        strand = row['strand']
        func   = row['bgc_role']
        color  = FUNC_COLORS.get(func, '#BDBDBD')
        s, e   = row['start'], row['end']
        w      = e - s
        yc     = ARROW_Y[strand]
        hl     = min(w * 0.30, span * 0.025)

        arrow = mpatches.FancyArrow(
            x=s if strand == '+' else e,
            y=yc - ARROW_H / 2,
            dx=w if strand == '+' else -w, dy=0,
            width=ARROW_H, length_includes_head=True,
            head_width=ARROW_H, head_length=hl,
            fc=color, ec='white', linewidth=0.4,
        )
        ax.add_patch(arrow)

        # Label only non-biosynthesis genes
        if func == 'biosynthesis':
            continue
        mid = (s + e) / 2
        short = gene_short_label(row['product'])

        if strand == '+':
            if mid > prev_end_plus + min_gap:
                ax.plot([mid, mid], [yc + ARROW_H / 2, LABEL_Y_PLUS - 0.05],
                        color='#888', linewidth=0.5, linestyle=':')
                ax.text(mid, LABEL_Y_PLUS, short,
                        ha='center', va='bottom', fontsize=4.8,
                        color='#222', clip_on=True)
                prev_end_plus = mid + len(short) * span * 0.012
        else:
            if mid > prev_end_minus + min_gap:
                ax.plot([mid, mid], [yc - ARROW_H / 2, LABEL_Y_MINUS + 0.05],
                        color='#888', linewidth=0.5, linestyle=':')
                ax.text(mid, LABEL_Y_MINUS, short,
                        ha='center', va='top', fontsize=4.8,
                        color='#222', clip_on=True)
                prev_end_minus = mid + len(short) * span * 0.012

    # X-axis ticks — adaptive spacing
    step = max(1000, round(span / 12 / 1000) * 1000)
    start_tick = (xlim[0] // step + 1) * step
    xticks = np.arange(start_tick, xlim[1], step)
    for xt in xticks:
        ax.plot([xt, xt], [-2.55, -2.40], color='#aaa', linewidth=0.5)
        ax.text(xt, -2.72, f'{xt/1e6:.3f} Mb',
                ha='center', fontsize=5.5, color='#555')


# ── Per-cluster figure ────────────────────────────────────────────────────────

def make_cluster_figure(bgc_name, df_all_sites, df_all_genes):
    label = BGC_META[bgc_name]
    genes = df_all_genes[df_all_genes['bgc_name'] == bgc_name].sort_values('start')
    cl_start = genes['start'].min() - FLANK
    cl_end   = genes['end'].max()   + FLANK
    xlim = (cl_start, cl_end)

    sites = get_cluster_sites(df_all_sites, cl_start, cl_end)
    print(f'  {bgc_name}: {cl_start:,}–{cl_end:,} bp | '
          f'4mC={(sites["mod_type"]=="4mC").sum()}, '
          f'6mA={(sites["mod_type"]=="6mA").sum()}')
    for tp in ['T1', 'T2', 'T3']:
        n = (sites['timepoint'] == tp).sum()
        print(f'    {tp}: {n} sites')

    fig, axes = plt.subplots(
        4, 1,
        figsize=(mm_to_inch(180), mm_to_inch(120)),
        gridspec_kw={'height_ratios': [1, 1, 1, 1.3], 'hspace': 0.06},
    )
    fig.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.07)
    fig.suptitle(f'{label} cluster — methylation site map',
                 fontsize=10, fontweight='bold')

    for i, tp in enumerate(['T1', 'T2', 'T3']):
        n = draw_methylation_track(axes[i], sites, tp, TP_COLORS[tp], xlim)

    draw_gene_track(axes[3], genes, xlim)

    # Methylation type legend (T1 panel)
    legend_handles = [
        mpatches.Patch(facecolor=COL_4mC,  edgecolor='none', label='GCCGGC (4mC)'),
        mpatches.Patch(facecolor=COL_BOTH, edgecolor='none', label='AAGCCCG (4mC/6mA dual)'),
        mpatches.Patch(facecolor=COL_6mA,  edgecolor='none', label='6mA'),
    ]
    axes[0].legend(handles=legend_handles, loc='upper right', fontsize=7,
                   frameon=False, ncol=3, handlelength=1.0)

    # Gene function legend — figure level
    func_patches = [
        mpatches.Patch(facecolor=v, edgecolor='none', label=k.capitalize())
        for k, v in FUNC_COLORS.items()
    ]
    fig.legend(handles=func_patches, loc='lower center', ncol=4,
               fontsize=7, frameon=False, handlelength=1.2,
               bbox_to_anchor=(0.5, 0.00))

    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS: BGC Cluster Methylation Tracks ===')

    print('Loading data...')
    df_all_sites = load_all_methylation()
    df_all_genes = load_all_genes()

    for bgc_name in ['act', 'cda', 'cpk', 'red']:
        print(f'\nProcessing {bgc_name.upper()}...')
        fig = make_cluster_figure(bgc_name, df_all_sites, df_all_genes)
        out_path = FIG_SUP_DIR / f'FigS_{bgc_name}_methylation_track'
        save_figure(fig, out_path)
        plt.close(fig)
        print(f'  Saved: {out_path}.pdf / .svg')

    print('\n=== Done ===')


if __name__ == '__main__':
    main()
