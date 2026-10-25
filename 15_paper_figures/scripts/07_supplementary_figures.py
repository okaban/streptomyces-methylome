# SUPERSEDED 2026-10-15: 旧世代の図生成スクリプト。現行の出荷図は
# panels_261013/ と supp_261014/ が産出し、figure_manifest.tsv に登録してある。
# 本スクリプトは撤回済みの入力（Exposed/Shielded 区分、退避版テーブル）を参照するため
# 実行を止めてある。再利用する場合は現行サイト集合で再計算し、産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 07_supplementary_figures.py は旧世代のため実行を止めています。"
)

#!/usr/bin/env python3
"""
Supplementary Figures S1-S11
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
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
import logomaker

PILEUP_DIR = METHYL / 'analysis' / 'pileup'
SAMPLES = {
    'T1': ['1-1', '1-2', '1-3'],
    'T2': ['2-1', '2-3', '2-4'],
    'T3': ['3-2', '3-3', '3-4'],
}
# modkit pileup BED: col3 = mod code, col9 = N_valid_cov, col10 = pct_mod
MOD_CODES = {'21839': '4mC', 'a': '6mA', 'm': '5mC'}


def _load_pileup(sample, mod_code=None):
    """Load a pileup BED file for a sample."""
    path = PILEUP_DIR / f'{sample}_pileup.bed'
    cols = ['chrom', 'start', 'end', 'mod', 'score', 'strand',
            'thick_s', 'thick_e', 'color', 'n_valid', 'pct_mod',
            'n_mod', 'n_canon', 'n_other', 'n_del', 'n_fail', 'n_diff', 'n_nocall']
    df = pd.read_csv(path, sep='\t', header=None, names=cols,
                     dtype={'mod': str})
    if mod_code is not None:
        df = df[df['mod'] == mod_code]
    return df


def fig_s1_qc():
    """S1: QC — coverage distribution and replicate reproducibility."""
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(mm_to_inch(174), mm_to_inch(160)))

    # Panel A: Coverage distribution per sample (4mC sites)
    ax = axes[0, 0]
    all_cov = []
    labels = []
    for tp, samples in SAMPLES.items():
        for samp in samples:
            df = _load_pileup(samp, mod_code='21839')
            cov = df['n_valid'].values
            all_cov.append(cov)
            labels.append(f'{samp} ({tp})')

    bp = ax.boxplot(all_cov, tick_labels=labels, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', lw=1.5))
    colors_tp = {'T1': '#FFCDD2', 'T2': '#C5CAE9', 'T3': '#C8E6C9'}
    for i, (tp, samples) in enumerate(SAMPLES.items()):
        for j in range(len(samples)):
            idx = i * len(samples) + j
            bp['boxes'][idx].set_facecolor(colors_tp[tp])
    ax.set_ylabel('Coverage (reads)')
    ax.set_title('4mC site coverage per sample', fontsize=10, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3, lw=0.5)
    add_panel_label(ax, 'a')

    # Panel B: Coverage distribution (6mA sites)
    ax = axes[0, 1]
    all_cov = []
    labels = []
    for tp, samples in SAMPLES.items():
        for samp in samples:
            df = _load_pileup(samp, mod_code='a')
            cov = df['n_valid'].values
            all_cov.append(cov)
            labels.append(f'{samp} ({tp})')

    bp = ax.boxplot(all_cov, tick_labels=labels, patch_artist=True, showfliers=False,
                    medianprops=dict(color='black', lw=1.5))
    for i, (tp, samples) in enumerate(SAMPLES.items()):
        for j in range(len(samples)):
            idx = i * len(samples) + j
            bp['boxes'][idx].set_facecolor(colors_tp[tp])
    ax.set_ylabel('Coverage (reads)')
    ax.set_title('6mA site coverage per sample', fontsize=10, fontweight='bold')
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='y', alpha=0.3, lw=0.5)
    add_panel_label(ax, 'b')

    # Panel C: Replicate reproducibility — T1 rep1 vs rep2 (4mC pct_mod)
    ax = axes[1, 0]
    r1 = _load_pileup('1-1', '21839').set_index(['start', 'strand'])['pct_mod']
    r2 = _load_pileup('1-2', '21839').set_index(['start', 'strand'])['pct_mod']
    common = r1.index.intersection(r2.index)
    ax.scatter(r1.loc[common], r2.loc[common], s=1, alpha=0.1, color=COL_4mC,
               rasterized=True)
    corr = np.corrcoef(r1.loc[common], r2.loc[common])[0, 1]
    ax.text(0.05, 0.95, f'r = {corr:.3f}\nn = {len(common):,}',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel('4mC % (rep 1-1)')
    ax.set_ylabel('4mC % (rep 1-2)')
    ax.set_title('T1 replicate reproducibility (4mC)', fontsize=10, fontweight='bold')
    ax.plot([0, 100], [0, 100], 'k--', lw=0.5, alpha=0.3)
    add_panel_label(ax, 'c')

    # Panel D: Replicate reproducibility — T1 rep1 vs rep2 (6mA pct_mod)
    ax = axes[1, 1]
    r1 = _load_pileup('1-1', 'a').set_index(['start', 'strand'])['pct_mod']
    r2 = _load_pileup('1-2', 'a').set_index(['start', 'strand'])['pct_mod']
    common = r1.index.intersection(r2.index)
    ax.scatter(r1.loc[common], r2.loc[common], s=1, alpha=0.1, color=COL_6mA,
               rasterized=True)
    corr = np.corrcoef(r1.loc[common], r2.loc[common])[0, 1]
    ax.text(0.05, 0.95, f'r = {corr:.3f}\nn = {len(common):,}',
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_xlabel('6mA % (rep 1-1)')
    ax.set_ylabel('6mA % (rep 1-2)')
    ax.set_title('T1 replicate reproducibility (6mA)', fontsize=10, fontweight='bold')
    ax.plot([0, 100], [0, 100], 'k--', lw=0.5, alpha=0.3)
    add_panel_label(ax, 'd')

    fig.tight_layout(h_pad=3, w_pad=2)
    save_figure(fig, FIG_SUP_DIR / 'FigS1_QC')


def fig_s2_all_meme_motifs():
    """S2: All MEME motifs (ranks 1-5)."""
    apply_style()
    meme_4mc = parse_meme_pwm(
        EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
        '07_motif_analysis' / 'meme_4mC' / 'meme.txt')
    meme_6ma = parse_meme_pwm(
        EPIGENOME / 'archive' / 'v1_weighted_minreps2' /
        '07_motif_analysis' / 'meme_6mA' / 'meme.txt')

    n_4mc = len(meme_4mc)
    n_6ma = len(meme_6ma)
    n_rows = max(n_4mc, n_6ma)

    fig, axes = plt.subplots(n_rows, 2, figsize=(mm_to_inch(174), mm_to_inch(40 * n_rows)))
    nuc_colors = {'A': '#43A047', 'C': '#1565C0', 'G': '#FFA000', 'T': '#E53935'}

    for i in range(n_rows):
        # 4mC column
        ax = axes[i, 0] if n_rows > 1 else axes[0]
        key = f'MEME-{i + 1}'
        if key in meme_4mc:
            m = meme_4mc[key]
            info = logomaker.transform_matrix(m['pwm'], from_type='probability',
                                              to_type='information')
            logomaker.Logo(info, ax=ax, color_scheme=nuc_colors)
            ax.set_title(f'4mC {key}: {m["name"]}  (n={m["nsites"]:,}, E={m["evalue"]})',
                         fontsize=8, fontweight='bold', color=COL_4mC)
            ax.set_ylim(0, 2.2)
            ax.set_ylabel('Bits', fontsize=7)
        else:
            ax.axis('off')

        # 6mA column
        ax = axes[i, 1] if n_rows > 1 else axes[1]
        if key in meme_6ma:
            m = meme_6ma[key]
            info = logomaker.transform_matrix(m['pwm'], from_type='probability',
                                              to_type='information')
            logomaker.Logo(info, ax=ax, color_scheme=nuc_colors)
            ax.set_title(f'6mA {key}: {m["name"]}  (n={m["nsites"]:,}, E={m["evalue"]})',
                         fontsize=8, fontweight='bold', color=COL_6mA)
            ax.set_ylim(0, 2.2)
            ax.set_ylabel('Bits', fontsize=7)
        else:
            ax.axis('off')

    if n_rows > 1:
        add_panel_label(axes[0, 0], 'a')
        add_panel_label(axes[0, 1], 'b')

    fig.tight_layout(h_pad=2)
    save_figure(fig, FIG_SUP_DIR / 'FigS2_all_MEME_motifs')


def fig_s3_rebase_heatmap():
    """S3: REBASE species × motif heatmap."""
    apply_style()
    df = pd.read_csv(
        EPIGENOME / '21_genuswide_motif_conservation' /
        'rebase_species_motif_matrix.csv')

    motif_cols = [c for c in df.columns if c not in ('species', 'has_any_rm')]
    mat = df.set_index('species')[motif_cols].astype(float)

    # Sort species by total motif count
    mat['total'] = mat.sum(axis=1)
    mat = mat.sort_values('total', ascending=False).drop(columns='total')

    fig, ax = plt.subplots(figsize=(mm_to_inch(174), mm_to_inch(280)))

    cmap = LinearSegmentedColormap.from_list('binary', ['#FFFFFF', '#1565C0'], N=2)
    sns.heatmap(mat, cmap=cmap, cbar=False, linewidths=0.3, linecolor='#E0E0E0',
                ax=ax, xticklabels=True, yticklabels=True)

    ax.set_xlabel('Recognition motif')
    ax.set_ylabel('Streptomyces species (REBASE)')
    ax.set_title('Motif conservation across Streptomyces (REBASE)',
                 fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelsize=5)
    ax.tick_params(axis='x', labelsize=8, rotation=45)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#1565C0', label='Present'),
        mpatches.Patch(facecolor='#FFFFFF', edgecolor='#BDBDBD', label='Absent'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS3_REBASE_heatmap')


def fig_s4_oe_ratio_detail():
    """S4: O/E ratio detail for all motifs."""
    apply_style()
    df = load_motif_summary()

    fig, ax = plt.subplots(figsize=(mm_to_inch(120), mm_to_inch(100)))

    motifs = df['motif'].values
    oe_vals = []
    for _, row in df.iterrows():
        oe_str = str(row['m145_oe'])
        try:
            oe_vals.append(float(oe_str.split('(')[0].strip().split()[0]))
        except (ValueError, IndexError):
            oe_vals.append(float(oe_str))
    oe_vals = np.array(oe_vals)

    order = np.argsort(oe_vals)
    motifs = motifs[order]
    oe_vals = oe_vals[order]

    y = np.arange(len(motifs))
    colors = [COL_BOTH if m == 'AAGCCCG' else
              COL_4mC if 'CCG' in m or 'CCGG' in m.upper() else
              COL_6mA for m in motifs]

    ax.barh(y, oe_vals, height=0.6, color=colors, edgecolor='white',
            linewidth=0.5, alpha=0.9, zorder=3)
    ax.axvline(1.0, color='black', linestyle='--', lw=0.8, alpha=0.5, zorder=2)

    for i in range(len(motifs)):
        ax.text(oe_vals[i] + 0.02, y[i], f'{oe_vals[i]:.2f}',
                va='center', fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels(motifs, fontfamily='monospace', fontsize=9)
    ax.set_xlabel('Observed / Expected ratio')
    ax.set_title('Motif O/E ratios in M145 genome', fontsize=11, fontweight='bold')
    ax.grid(axis='x', alpha=0.3, lw=0.5)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS4_OE_ratios')


def fig_s5_sc_rs17645_homology():
    """S5: SC_RS17645 BLAST and Foldseek homologs."""
    apply_style()
    df_blast = pd.read_csv(
        EPIGENOME / '13_sc_rs17645_analysis' / 'blast_top_hits.csv')
    df_fold = pd.read_csv(
        EPIGENOME / '13_sc_rs17645_analysis' / 'foldseek_structural_homologs.csv')

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(mm_to_inch(174), mm_to_inch(200)))

    # Panel A: BLAST top hits
    df_b = df_blast.head(15).copy()
    df_b = df_b.sort_values('pct_identity')
    y = np.arange(len(df_b))
    labels = [f"{row['organism'][:40]}" for _, row in df_b.iterrows()]

    ax1.barh(y, df_b['pct_identity'].values, height=0.6,
             color=COL_6mA, alpha=0.8, edgecolor='white', zorder=3)
    for i, (_, row) in enumerate(df_b.iterrows()):
        ax1.text(row['pct_identity'] + 0.3, y[i],
                 f"{row['pct_identity']:.1f}%  E={row['evalue']:.0e}",
                 va='center', fontsize=6.5)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=7)
    ax1.set_xlabel('Sequence identity (%)')
    ax1.set_title('SC_RS17645 BLASTp top hits', fontsize=11, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3, lw=0.5)
    add_panel_label(ax1, 'a')

    # Panel B: Foldseek structural homologs
    df_f = df_fold.head(10).copy()
    df_f = df_f.sort_values('seq_identity_pct')
    y = np.arange(len(df_f))
    labels = [f"{row['pdb_id']}:{row['chain']} {str(row['description'])[:35]}"
              for _, row in df_f.iterrows()]

    ax2.barh(y, df_f['seq_identity_pct'].values, height=0.6,
             color=COL_BOTH, alpha=0.8, edgecolor='white', zorder=3)
    for i, (_, row) in enumerate(df_f.iterrows()):
        ax2.text(max(row['seq_identity_pct'], 0) + 0.3, y[i],
                 f"{row['seq_identity_pct']:.1f}%  prob={row['prob']:.0f}",
                 va='center', fontsize=6.5)
    ax2.set_yticks(y)
    ax2.set_yticklabels(labels, fontsize=7)
    ax2.set_xlabel('Sequence identity (%)')
    ax2.set_title('SC_RS17645 Foldseek structural homologs', fontsize=11,
                  fontweight='bold')
    ax2.grid(axis='x', alpha=0.3, lw=0.5)
    add_panel_label(ax2, 'b')

    fig.tight_layout(h_pad=3)
    save_figure(fig, FIG_SUP_DIR / 'FigS5_SC_RS17645_homology')


def fig_s6_mtase_expression():
    """S6: All 22 MTase genes expression heatmap."""
    apply_style()
    df = load_mtase_expression()

    fc_cols = ['log2FC_T2vsT1', 'log2FC_T3vsT1', 'log2FC_T3vsT2']
    padj_cols = ['padj_T2vsT1', 'padj_T3vsT1', 'padj_T3vsT2']

    mat = df.set_index('locus_tag')[fc_cols].copy()
    mat.columns = ['T2 vs T1', 'T3 vs T1', 'T3 vs T2']

    # Fill NaN with 0
    mat = mat.fillna(0)

    # Sort by T2vsT1
    mat = mat.sort_values('T2 vs T1')

    fig, ax = plt.subplots(figsize=(mm_to_inch(120), mm_to_inch(180)))

    vmax = max(abs(mat.values.min()), abs(mat.values.max()), 3)
    cmap = LinearSegmentedColormap.from_list('rg', ['#1565C0', '#FFFFFF', '#E53935'])

    im = sns.heatmap(mat, cmap=cmap, center=0, vmin=-vmax, vmax=vmax,
                     ax=ax, linewidths=0.5, linecolor='white',
                     annot=True, fmt='.2f', annot_kws={'fontsize': 7},
                     cbar_kws={'label': r'$\log_2$FC', 'shrink': 0.6})

    # Mark significant
    padj_mat = df.set_index('locus_tag')[padj_cols].copy()
    padj_mat.columns = mat.columns
    padj_mat = padj_mat.loc[mat.index]
    for i in range(len(mat)):
        for j in range(len(mat.columns)):
            p = padj_mat.iloc[i, j]
            if pd.notna(p) and p < 0.05:
                ax.text(j + 0.85, i + 0.15, '*', ha='center', va='center',
                        fontsize=10, fontweight='bold', color='black')

    # Highlight SC_RS17645
    target_idx = list(mat.index).index('SC_RS17645') if 'SC_RS17645' in mat.index else -1
    if target_idx >= 0:
        ax.add_patch(plt.Rectangle((0, target_idx), len(mat.columns), 1,
                                   fill=False, edgecolor=COL_BOTH, lw=2))

    ax.set_title('MTase gene expression changes', fontsize=12, fontweight='bold')
    ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=8)
    ax.text(-0.15, 1.01, '* padj < 0.05', transform=ax.transAxes,
            fontsize=8, fontstyle='italic')

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS6_MTase_expression')


def fig_s7_tier2_sensitivity():
    """S7: TF BS Tier 2 sensitivity analysis."""
    apply_style()
    df_bs = load_tf_binding_sites()
    df_unique = load_methylation_unique_positions()
    meth_pos = df_unique['position'].values

    fig, ax = plt.subplots(figsize=(mm_to_inch(120), mm_to_inch(100)))

    results = {}
    for tier in [1, 2]:
        bs_tier = df_bs[df_bs['tier'] == tier]
        # Count sites with absolute coordinates (FIMO only has valid abs coords)
        fimo = bs_tier[bs_tier['BS_source'] == 'FIMO_curated_motif']
        n_bs = len(fimo)
        if n_bs == 0:
            continue
        n_overlap = 0
        for _, row in fimo.iterrows():
            if pd.notna(row['BS_start']) and pd.notna(row['BS_end']):
                mask = (meth_pos >= row['BS_start']) & (meth_pos <= row['BS_end'])
                n_overlap += mask.sum()
        results[f'Tier {tier}\n(FIMO, n={n_bs})'] = n_overlap

    # Also show all Tier 1 resolved
    fig4 = importlib.import_module('04_figure4_tf_bs_depletion')
    bs_res = fig4.resolve_bs_jeong2016()
    n_t1 = len(bs_res)
    n_ov_t1 = 0
    for _, row in bs_res.iterrows():
        mask = (meth_pos >= row['abs_start']) & (meth_pos <= row['abs_end'])
        n_ov_t1 += mask.sum()
    results[f'Tier 1 resolved\n(all sources, n={n_t1})'] = n_ov_t1

    labels = list(results.keys())
    values = list(results.values())
    bars = ax.bar(labels, values, color=[COL_4mC, COL_6mA, COL_BOTH],
                  edgecolor='white', alpha=0.9, zorder=3)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 1,
                str(v), ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_ylabel('Methylation sites in BS')
    ax.set_title('BS methylation overlaps by tier/source', fontsize=11,
                 fontweight='bold')
    ax.grid(axis='y', alpha=0.3, lw=0.5)

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS7_Tier2_sensitivity')


def fig_s8_composition_correction():
    """S8: Composition correction details."""
    apply_style()
    df = pd.read_csv(
        BASE / '13_TF_binding-site' / 'analysis' /
        '02_TF_BS_methylation_260207_v1' / 'tables' /
        'BS_methylation_composition_correction.csv')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(mm_to_inch(174), mm_to_inch(100)))

    # Panel A: Fold enrichment by method and mod type
    methods = df['Method'].values
    x = np.arange(len(methods))
    width = 0.25

    ax1.bar(x - width, df['Fold_4mC'].values, width, color=COL_4mC, label='4mC',
            edgecolor='white', alpha=0.9, zorder=3)
    ax1.bar(x, df['Fold_6mA'].values, width, color=COL_6mA, label='6mA',
            edgecolor='white', alpha=0.9, zorder=3)
    ax1.bar(x + width, df['Fold_any'].values, width, color=COL_BOTH, label='Any',
            edgecolor='white', alpha=0.9, zorder=3)

    ax1.axhline(1.0, color='black', linestyle='--', lw=0.8, alpha=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods, fontsize=7, rotation=30, ha='right')
    ax1.set_ylabel('Fold enrichment (Obs / Exp)')
    ax1.set_title('Fold by correction method', fontsize=10, fontweight='bold')
    ax1.legend(fontsize=8, frameon=False)
    ax1.grid(axis='y', alpha=0.3, lw=0.5)
    add_panel_label(ax1, 'a')

    # Panel B: Observed vs Expected counts
    ax2.scatter(df['Exp_any'], df['Obs_any'], s=80, c=COL_BOTH, zorder=3,
                edgecolor='white')
    for _, row in df.iterrows():
        ax2.annotate(row['Method'], (row['Exp_any'], row['Obs_any']),
                     fontsize=7, ha='left', va='bottom',
                     xytext=(5, 3), textcoords='offset points')
    lim = max(df['Exp_any'].max(), df['Obs_any'].max()) * 1.2
    ax2.plot([0, lim], [0, lim], 'k--', lw=0.5, alpha=0.3)
    ax2.set_xlabel('Expected methylation sites')
    ax2.set_ylabel('Observed methylation sites')
    ax2.set_title('Observed vs Expected', fontsize=10, fontweight='bold')
    ax2.grid(alpha=0.3, lw=0.5)
    add_panel_label(ax2, 'b')

    fig.tight_layout(w_pad=3)
    save_figure(fig, FIG_SUP_DIR / 'FigS8_composition_correction')


def fig_s9_promoter_methylation():
    """S9: Promoter methylation density profiles by mod type."""
    apply_style()
    df_tss = load_tss_jeong2016()
    df_unique = load_methylation_unique_positions()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(mm_to_inch(174), mm_to_inch(90)))

    for ax, mod, color, label in [(ax1, '4mC', COL_4mC, '4mC'),
                                   (ax2, '6mA', COL_6mA, '6mA')]:
        positions = df_unique[df_unique['mod_type'] == mod]['position'].values
        pos_sorted = np.sort(positions)
        tss_vals = df_tss['tss'].values.astype(int)
        strands = df_tss['strand'].values

        window = 200
        bin_size = 5
        bins = np.arange(-window, window + 1, bin_size)
        counts = np.zeros(len(bins) - 1)

        for tss, strand in zip(tss_vals, strands):
            if strand == '+':
                local = pos_sorted - tss
            else:
                local = tss - pos_sorted
            in_range = local[(local >= -window) & (local < window)]
            c, _ = np.histogram(in_range, bins=bins)
            counts += c

        # Density: sites per kb per gene
        density = counts / (bin_size / 1000) / len(tss_vals)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        ax.fill_between(bin_centers, density, alpha=0.3, color=color)
        ax.plot(bin_centers, density, color=color, lw=1.5)
        ax.axvline(0, color='gray', linestyle='--', lw=0.8, alpha=0.5)
        ax.axvspan(-12, -7, alpha=0.15, color='#FFC107', label='σ-10 box')

        ax.set_xlabel('Position relative to TSS (bp)')
        ax.set_ylabel(f'{label} density\n(sites/kb/gene)')
        ax.set_title(f'{label} promoter profile', fontsize=10, fontweight='bold')
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.3, lw=0.5)

    add_panel_label(ax1, 'a')
    add_panel_label(ax2, 'b')

    fig.tight_layout(w_pad=3)
    save_figure(fig, FIG_SUP_DIR / 'FigS9_promoter_methylation')


def fig_s10_methylation_vs_expression():
    """S10: Genome-wide methylation vs expression scatter."""
    apply_style()
    df = load_integrated_expression()

    fig, axes = plt.subplots(1, 3, figsize=(mm_to_inch(174), mm_to_inch(80)))

    comparisons = [
        ('log2FC_T2vsT1', 'T2 vs T1'),
        ('log2FC_T3vsT1', 'T3 vs T1'),
        ('log2FC_T3vsT2', 'T3 vs T2'),
    ]

    # Use total methylation sites as x-axis
    for ax, (col, title) in zip(axes, comparisons):
        # Sum all methylation across timepoints
        meth_cols = [c for c in df.columns if 'n_sites' in c.lower() or
                     c.startswith('n_4mC') or c.startswith('n_6mA')]
        if not meth_cols:
            # Try different column naming
            meth_cols = [c for c in df.columns if '4mC' in c or '6mA' in c]
            meth_cols = [c for c in meth_cols if 'n_' in c]

        if meth_cols:
            df['total_meth'] = df[meth_cols].sum(axis=1)
        else:
            # Fallback: use any methylation indicator columns
            df['total_meth'] = 0

        if col not in df.columns:
            ax.text(0.5, 0.5, f'{col}\nnot found', transform=ax.transAxes,
                    ha='center', va='center')
            continue

        has_meth = df['total_meth'] > 0
        no_meth = ~has_meth

        ax.scatter(df.loc[no_meth, 'total_meth'] + np.random.normal(0, 0.1, no_meth.sum()),
                   df.loc[no_meth, col], s=1, alpha=0.05, color=COL_GRAY,
                   rasterized=True, label=f'No meth (n={no_meth.sum():,})')
        ax.scatter(df.loc[has_meth, 'total_meth'],
                   df.loc[has_meth, col], s=5, alpha=0.3, color=COL_BOTH,
                   rasterized=True, label=f'Methylated (n={has_meth.sum():,})')

        ax.axhline(0, color='gray', linestyle=':', lw=0.5)
        ax.set_xlabel('Total methylation sites')
        ax.set_ylabel(f'$\\log_2$FC ({title})')
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.legend(fontsize=6, frameon=False, loc='upper right', markerscale=3)

    add_panel_label(axes[0], 'a')
    add_panel_label(axes[1], 'b')
    add_panel_label(axes[2], 'c')

    fig.tight_layout(w_pad=2)
    save_figure(fig, FIG_SUP_DIR / 'FigS10_methylation_vs_expression')


def fig_s11_bs_overlap_detail():
    """S11: Detailed BS-methylation overlap pairs."""
    apply_style()
    # Load the overlap table we already generated
    st4_path = TABLE_SUP_DIR / 'ST4_BS_methylation_overlaps.tsv'
    if not st4_path.exists():
        print('  ST4 table not found, skipping S11')
        return

    df = pd.read_csv(st4_path, sep='\t')
    if len(df) == 0:
        print('  No overlaps found, skipping S11')
        return

    fig, ax = plt.subplots(figsize=(mm_to_inch(174), mm_to_inch(100)))

    # Summary: overlaps per TF
    tf_counts = df.groupby('TF_name').size().sort_values(ascending=True)
    y = np.arange(len(tf_counts))

    colors = [COL_BOTH if c > 2 else COL_6mA for c in tf_counts.values]
    ax.barh(y, tf_counts.values, height=0.6, color=colors,
            edgecolor='white', alpha=0.9, zorder=3)
    for i, v in enumerate(tf_counts.values):
        ax.text(v + 0.2, y[i], str(v), va='center', fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels(tf_counts.index, fontsize=8)
    ax.set_xlabel('Methylation sites within binding sites')
    ax.set_title('TF binding site methylation overlaps', fontsize=11,
                 fontweight='bold')
    ax.grid(axis='x', alpha=0.3, lw=0.5)

    # Annotation
    ax.text(0.97, 0.05,
            f'Total: {len(df)} overlaps\nacross {df["TF_name"].nunique()} TFs',
            transform=ax.transAxes, fontsize=9, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='#F3E5F5', edgecolor=COL_BOTH,
                      alpha=0.9))

    fig.tight_layout()
    save_figure(fig, FIG_SUP_DIR / 'FigS11_BS_overlap_detail')


def main():
    print('=== Supplementary Figures S1-S11 ===')
    FIG_SUP_DIR.mkdir(parents=True, exist_ok=True)

    funcs = [
        ('S1', fig_s1_qc),
        ('S2', fig_s2_all_meme_motifs),
        ('S3', fig_s3_rebase_heatmap),
        ('S4', fig_s4_oe_ratio_detail),
        ('S5', fig_s5_sc_rs17645_homology),
        ('S6', fig_s6_mtase_expression),
        ('S7', fig_s7_tier2_sensitivity),
        ('S8', fig_s8_composition_correction),
        ('S9', fig_s9_promoter_methylation),
        ('S10', fig_s10_methylation_vs_expression),
        ('S11', fig_s11_bs_overlap_detail),
    ]

    for name, func in funcs:
        print(f'\n--- {name} ---')
        try:
            func()
        except Exception as e:
            print(f'  ERROR in {name}: {e}')
            import traceback
            traceback.print_exc()

    print('\n=== All supplementary figures complete ===')


if __name__ == '__main__':
    main()
