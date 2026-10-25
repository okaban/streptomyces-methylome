# SUPERSEDED 2026-10-15: 旧世代の図生成スクリプト。現行の出荷図は
# panels_261013/ と supp_261014/ が産出し、figure_manifest.tsv に登録してある。
# 本スクリプトは撤回済みの入力（Exposed/Shielded 区分、退避版テーブル）を参照するため
# 実行を止めてある。再利用する場合は現行サイト集合で再計算し、産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 23_figS_bgc_integrated_track.py は旧世代のため実行を止めています。"
)

#!/usr/bin/env python3
"""
FigS_BGC_Integrated_Track: BGC Cluster Transcriptome + Methylation Coverage Tracks

Combines RNA-seq coverage (T1/T2/T3 mean±SE) with methylation site maps
for all 4 major BGC clusters: ACT, CDA, CPK, RED

Layout per cluster (6 rows):
  - Row 1: RNA-seq coverage — T1/T2/T3 as colored lines with SE shading
  - Row 2: T1 methylation lollipops
  - Row 3: T2 methylation lollipops
  - Row 4: T3 methylation lollipops
  - Row 5: Gene annotation arrows
"""

import sys, importlib, subprocess
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

# ── Paths ─────────────────────────────────────────────────────────────────────

BAM_DIR = BASE / '02_alignment' / 'analysis' / '02_alignment_260127_v1' / 'bam'
BGC_DIR = EPIGENOME / '47_BGC_methylation_geographic_test' / 'tables'
CENSUS_4mC = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census_v1.csv'
CENSUS_6mA = EPIGENOME / '23_expanded_motif_search' / '6mA_final_census_v1.csv'
CHROM = 'NC_003888.3'

FLANK = 2000
BIN_SIZE = 200   # bp per coverage bin

# ── Sample mapping ────────────────────────────────────────────────────────────

TP_BAMS = {
    'T1': ['M145_1_1', 'M145_1_2', 'M145_1_3'],
    'T2': ['M145_2_1', 'M145_2_3', 'M145_2_4'],
    'T3': ['M145_3_2', 'M145_3_3', 'M145_3_4'],
}

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


# ── Coverage computation ───────────────────────────────────────────────────────

def compute_coverage(sample_names, chrom, start, end):
    """Run samtools depth for given samples and region; return position→depth array."""
    bam_paths = [
        str(BAM_DIR / f'{s}.Aligned.sortedByCoord.out.bam')
        for s in sample_names
    ]
    region = f'{chrom}:{start}-{end}'
    cmd = ['samtools', 'depth', '-a', '-r', region] + bam_paths
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    positions, depths = [], []
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        pos = int(parts[1])
        # Average depth across samples (columns 3+)
        dp_vals = [float(x) for x in parts[2:] if x.strip()]
        positions.append(pos)
        depths.append(np.mean(dp_vals))

    if not positions:
        return np.array([]), np.array([])
    return np.array(positions), np.array(depths)


def compute_coverage_se(sample_names, chrom, start, end):
    """Return positions, mean depth, and SE across replicates."""
    bam_paths = [
        str(BAM_DIR / f'{s}.Aligned.sortedByCoord.out.bam')
        for s in sample_names
    ]
    region = f'{chrom}:{start}-{end}'
    cmd = ['samtools', 'depth', '-a', '-r', region] + bam_paths
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    positions, means, ses = [], [], []
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        pos = int(parts[1])
        dp_vals = np.array([float(x) for x in parts[2:] if x.strip()])
        positions.append(pos)
        means.append(dp_vals.mean())
        ses.append(dp_vals.std(ddof=1) / np.sqrt(len(dp_vals)))

    return np.array(positions), np.array(means), np.array(ses)


def bin_coverage(positions, values, start, end, bin_size=BIN_SIZE):
    """Average values into bins of bin_size bp."""
    bins = np.arange(start, end + bin_size, bin_size)
    bin_centers = bins[:-1] + bin_size / 2
    bin_means = np.zeros(len(bin_centers))
    for i in range(len(bin_centers)):
        mask = (positions >= bins[i]) & (positions < bins[i + 1])
        if mask.sum() > 0:
            bin_means[i] = values[mask].mean()
    return bin_centers, bin_means


# ── Data loading ───────────────────────────────────────────────────────────────

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


# ── Drawing functions ──────────────────────────────────────────────────────────

def draw_coverage_track(ax, xlim, label='RNA-seq coverage'):
    """Draw RNA-seq coverage for T1/T2/T3 with SE shading."""
    start, end = int(xlim[0]), int(xlim[1])

    max_cov = 0
    for tp, tp_color in TP_COLORS.items():
        samples = TP_BAMS[tp]
        positions, means, ses = compute_coverage_se(samples, CHROM, start, end)
        if len(positions) == 0:
            continue
        # Bin for smoother display
        centers, bin_means = bin_coverage(positions, means, start, end)
        _, bin_ses = bin_coverage(positions, ses, start, end)

        ax.fill_between(centers, bin_means - bin_ses, bin_means + bin_ses,
                        color=tp_color, alpha=0.18)
        ax.plot(centers, bin_means, color=tp_color, linewidth=1.4,
                label=tp, alpha=0.90)
        max_cov = max(max_cov, (bin_means + bin_ses).max())

    ax.set_xlim(*xlim)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Read depth', fontsize=7, labelpad=2)
    ax.set_xticklabels([])
    ax.tick_params(bottom=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0, color='#bbb', linewidth=0.5)

    # Legend
    tp_handles = [
        mpatches.Patch(facecolor=TP_COLORS[tp], edgecolor='none', label=tp)
        for tp in ['T1', 'T2', 'T3']
    ]
    ax.legend(handles=tp_handles, loc='upper right', fontsize=7,
              frameon=False, ncol=3, handlelength=0.8)
    ax.text(0.012, 0.88, 'RNA-seq', transform=ax.transAxes,
            fontsize=8, fontweight='bold', color='#37474F', va='top')


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
    min_gap = span * 0.05

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

    step = max(1000, round(span / 12 / 1000) * 1000)
    start_tick = (xlim[0] // step + 1) * step
    xticks = np.arange(start_tick, xlim[1], step)
    for xt in xticks:
        ax.plot([xt, xt], [-2.55, -2.40], color='#aaa', linewidth=0.5)
        ax.text(xt, -2.72, f'{xt/1e6:.3f} Mb',
                ha='center', fontsize=5.5, color='#555')


# ── Per-cluster figure ─────────────────────────────────────────────────────────

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

    # 6 rows: coverage + T1 + T2 + T3 meth + genes
    fig, axes = plt.subplots(
        5, 1,
        figsize=(mm_to_inch(180), mm_to_inch(155)),
        gridspec_kw={'height_ratios': [1.8, 1, 1, 1, 1.3], 'hspace': 0.06},
    )
    fig.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.07)
    fig.suptitle(f'{label} cluster — RNA-seq & methylation track',
                 fontsize=10, fontweight='bold')

    # RNA-seq coverage
    print(f'  Computing RNA-seq coverage...')
    draw_coverage_track(axes[0], xlim)

    # Methylation tracks
    legend_handles = [
        mpatches.Patch(facecolor=COL_4mC,  edgecolor='none', label='GCCGGC (4mC)'),
        mpatches.Patch(facecolor=COL_BOTH, edgecolor='none', label='AAGCCCG (4mC/6mA dual)'),
        mpatches.Patch(facecolor=COL_6mA,  edgecolor='none', label='6mA'),
    ]
    for i, tp in enumerate(['T1', 'T2', 'T3']):
        draw_methylation_track(axes[i + 1], sites, tp, TP_COLORS[tp], xlim)

    axes[1].legend(handles=legend_handles, loc='upper right', fontsize=7,
                   frameon=False, ncol=3, handlelength=1.0)

    # Gene annotation
    draw_gene_track(axes[4], genes, xlim)

    # Gene function legend at figure level
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
    print('=== FigS: BGC Integrated Transcriptome + Methylation Tracks ===')

    print('Loading methylation data...')
    df_all_sites = load_all_methylation()
    df_all_genes = load_all_genes()

    for bgc_name in ['act', 'cda', 'cpk', 'red']:
        print(f'\nProcessing {bgc_name.upper()}...')
        fig = make_cluster_figure(bgc_name, df_all_sites, df_all_genes)
        out_path = FIG_SUP_DIR / f'FigS_{bgc_name}_integrated_track'
        save_figure(fig, out_path)
        plt.close(fig)
        print(f'  Saved: {out_path}.pdf / .svg')

    print('\n=== Done ===')


if __name__ == '__main__':
    main()
