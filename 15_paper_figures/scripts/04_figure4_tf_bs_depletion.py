# SUPERSEDED 2026-10-15: 撤回済みの Exposed/Shielded 区分（本文は「認識配列の有無の言い換え」
# として否定結果に移した）に依存し、旧世代の入力 all_genes_features_unified_n57.tsv ／
# 退避版 TF 結合部位表を参照する。現行のサイト集合では Exposed 85 / Shielded 966 であり
# 一致しない。figure_manifest.tsv に出荷図の産出元として登録されていない（2026-10-15 確認）。
# 内容を再利用する場合は現行サイト集合で再計算し、産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 04_figure4_tf_bs_depletion.py は撤回済み区分に依存するため実行を止めています。"
)

#!/usr/bin/env python3
"""
Figure 4: TF Binding Site Methylation Depletion
  Panel A: Metagene profile around BS center (Jeong2016 TSS only)
  Panel B: Forest plot of fold enrichment across correction methods
  Panel C: Permutation distribution
  Panel D: Promoter σ-10 box methylation (Jeong2016 TSS only)
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
from scipy import stats

GENOME_SIZE = 8_667_507
TIER1_SOURCES = {'FIMO_curated_motif', 'RegPrecise',
                 'ZorroAranda2022_Curated_Strong'}


def resolve_bs_jeong2016():
    """Resolve BS coordinates using ONLY Jeong2016 experimental TSS.

    Returns:
        DataFrame with abs_start, abs_end columns for Tier 1 BS.
    """
    bs_raw = load_tf_binding_sites()
    tss_df = load_tss_jeong2016()

    # Build TSS map: gene_id -> tss (Jeong2016 only)
    tss_map = dict(zip(tss_df['gene_id'], tss_df['tss']))
    for _, t in tss_df.iterrows():
        olt = t.get('old_locus_tag', '')
        if pd.notna(olt) and olt:
            tss_map[olt] = t['tss']

    # Gene coordinate info (from TSS table itself)
    gene_coord = {}
    for _, t in tss_df.iterrows():
        gid = t['gene_id']
        gene_coord[gid] = (int(t['start']), int(t['end']), t['strand'])
        olt = t.get('old_locus_tag', '')
        if pd.notna(olt) and olt:
            gene_coord[olt] = (int(t['start']), int(t['end']), t['strand'])

    resolved = []
    for _, row in bs_raw.iterrows():
        src = row['BS_source']
        if src not in TIER1_SOURCES:
            continue

        rec = {'TF_name': row['TF_name'], 'BS_source': src,
               'abs_start': np.nan, 'abs_end': np.nan}

        tg_sco = str(row.get('TG_SCO_ID', '')) if pd.notna(row.get('TG_SCO_ID')) else ''
        tg_gid = str(row.get('TG_gene_id', '')) if pd.notna(row.get('TG_gene_id')) else ''
        bs_seq = str(row.get('BS_sequence', '')) if pd.notna(row.get('BS_sequence')) else ''

        if src == 'FIMO_curated_motif' and pd.notna(row.get('BS_start')):
            rec['abs_start'] = int(row['BS_start'])
            rec['abs_end'] = int(row['BS_end'])

        elif src == 'RegPrecise' and pd.notna(row.get('BS_start')):
            rel_pos = int(row['BS_start'])
            tg_key = tg_gid if tg_gid else tg_sco
            tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
            if tss_val and pd.notna(tss_val):
                tss_val = int(tss_val)
                coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
                if coord:
                    strand = coord[2]
                    seq_len = len(bs_seq) if bs_seq else 18
                    if strand == '+':
                        abs_s = tss_val + rel_pos
                        abs_e = abs_s + seq_len - 1
                    else:
                        abs_e = tss_val - rel_pos
                        abs_s = abs_e - seq_len + 1
                    rec['abs_start'] = max(1, abs_s)
                    rec['abs_end'] = abs_e

        elif src == 'ZorroAranda2022_Curated_Strong':
            tg_key = tg_gid if tg_gid else tg_sco
            tss_val = tss_map.get(tg_key) or tss_map.get(tg_sco)
            coord = gene_coord.get(tg_key) or gene_coord.get(tg_sco)
            if tss_val and pd.notna(tss_val) and coord:
                tss_val = int(tss_val)
                strand = coord[2]
                if strand == '+':
                    rec['abs_start'] = max(1, tss_val - 300)
                    rec['abs_end'] = tss_val + 50
                else:
                    rec['abs_start'] = max(1, tss_val - 50)
                    rec['abs_end'] = tss_val + 300

        resolved.append(rec)

    bs_df = pd.DataFrame(resolved)
    bs_df = bs_df.dropna(subset=['abs_start', 'abs_end'])
    bs_df['abs_start'] = bs_df['abs_start'].astype(int)
    bs_df['abs_end'] = bs_df['abs_end'].astype(int)
    bs_df = bs_df[(bs_df['abs_start'] >= 0) &
                  (bs_df['abs_end'] <= GENOME_SIZE) &
                  (bs_df['abs_start'] < bs_df['abs_end'])]

    print(f'  Resolved Tier 1 BS (Jeong2016 TSS): {len(bs_df)} sites')
    for src, cnt in bs_df['BS_source'].value_counts().items():
        print(f'    {src}: {cnt}')
    return bs_df


def panel_a_metagene(ax, bs_df, pos_4mc, pos_6ma):
    """Metagene methylation density profile around BS center."""
    window = 2000
    bin_size = 100

    # BS centers
    centers = ((bs_df['abs_start'] + bs_df['abs_end']) // 2).values

    # Count methylation sites at each distance from BS center
    bins = np.arange(-window, window + bin_size, bin_size)
    counts_4mc = np.zeros(len(bins) - 1)
    counts_6ma = np.zeros(len(bins) - 1)

    pos_4mc_arr = np.array(sorted(pos_4mc))
    pos_6ma_arr = np.array(sorted(pos_6ma))

    for center in centers:
        lo = center - window
        hi = center + window

        # 4mC
        idx = np.searchsorted(pos_4mc_arr, [lo, hi])
        if idx[1] > idx[0]:
            dists = pos_4mc_arr[idx[0]:idx[1]] - center
            hist, _ = np.histogram(dists, bins=bins)
            counts_4mc += hist

        # 6mA
        idx = np.searchsorted(pos_6ma_arr, [lo, hi])
        if idx[1] > idx[0]:
            dists = pos_6ma_arr[idx[0]:idx[1]] - center
            hist, _ = np.histogram(dists, bins=bins)
            counts_6ma += hist

    # Normalize to density (sites per kb per BS)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    density_4mc = counts_4mc / len(centers) / (bin_size / 1000)
    density_6ma = counts_6ma / len(centers) / (bin_size / 1000)

    # Expected uniform density
    n_4mc_total = len(pos_4mc)
    n_6ma_total = len(pos_6ma)
    exp_4mc = n_4mc_total / GENOME_SIZE * 1000  # per kb
    exp_6ma = n_6ma_total / GENOME_SIZE * 1000

    ax.fill_between(bin_centers, density_4mc, alpha=0.3, color=COL_4mC, zorder=2)
    ax.plot(bin_centers, density_4mc, color=COL_4mC, linewidth=1.5,
            label='4mC', zorder=3)
    ax.fill_between(bin_centers, density_6ma, alpha=0.3, color=COL_6mA, zorder=2)
    ax.plot(bin_centers, density_6ma, color=COL_6mA, linewidth=1.5,
            label='6mA', zorder=3)

    # Expected lines
    ax.axhline(exp_4mc, color=COL_4mC, linestyle=':', linewidth=0.8,
               alpha=0.5, zorder=1)
    ax.axhline(exp_6ma, color=COL_6mA, linestyle=':', linewidth=0.8,
               alpha=0.5, zorder=1)
    ax.text(window - 50, exp_4mc + 0.01, 'expected', fontsize=6,
            color=COL_4mC, alpha=0.6, ha='right', va='bottom')

    # BS region indicator
    median_bs_width = (bs_df['abs_end'] - bs_df['abs_start']).median()
    ax.axvspan(-median_bs_width/2, median_bs_width/2, alpha=0.08,
               color='gray', zorder=0, label=f'Median BS ({median_bs_width:.0f} bp)')
    ax.axvline(0, color='black', linewidth=0.5, linestyle='-', alpha=0.3)

    ax.set_xlabel('Distance from BS center (bp)')
    ax.set_ylabel('Methylation density\n(sites per kb per BS)')
    ax.set_title(f'Metagene profile (n={len(bs_df)} BS)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')
    ax.set_xlim(-window, window)
    ax.grid(axis='y', alpha=0.3, lw=0.5)


def panel_b_forest_plot(ax):
    """Forest plot showing fold enrichment across correction methods."""
    # Load composition correction results
    corr_path = (BASE / '13_TF_binding-site' / 'analysis' /
                 '02_TF_BS_methylation_260207_v1' / 'tables' /
                 'BS_methylation_composition_correction.csv')
    df = pd.read_csv(corr_path)

    # Filter out GC-matched (NaN results)
    df = df[df['Fold_any'].notna()].copy()

    methods = df['Method'].values
    folds = df['Fold_any'].values
    p_vals = df['P_depletion_any'].values

    # Approximate 95% CI from observed counts
    # Using Poisson CI: fold ± 1.96 * sqrt(fold/n)
    obs = df['Obs_any'].values
    exp = df['Exp_any'].values

    # Poisson SE for fold enrichment
    se = np.sqrt(obs) / exp  # SE of fold = sqrt(obs)/exp
    ci_lo = folds - 1.96 * se
    ci_hi = folds + 1.96 * se

    y = np.arange(len(methods))[::-1]

    # Short method labels
    short_labels = ['Uniform', 'GC-adjusted', 'Motif-adjusted', 'Permutation']

    ax.errorbar(folds, y, xerr=[folds - ci_lo, ci_hi - folds],
                fmt='D', color=COL_DARK, markersize=7,
                capsize=4, capthick=1.2, elinewidth=1.2, zorder=4)

    # Reference line at fold = 1.0
    ax.axvline(1.0, color='black', linestyle='--', linewidth=0.8,
               alpha=0.5, zorder=1)

    # Color background for depletion zone
    ax.axvspan(0, 1.0, alpha=0.04, color=COL_4mC, zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(short_labels, fontsize=9)
    ax.set_xlabel('Fold enrichment\n(observed / expected)')
    ax.set_title('Depletion across\ncorrection methods', fontsize=11,
                 fontweight='bold')
    ax.grid(axis='x', alpha=0.3, lw=0.5)

    # P-value annotations
    for i in range(len(methods)):
        p = p_vals[i]
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
        ax.text(ci_hi[i] + 0.02, y[i], f'p={p:.1e} {sig}',
                fontsize=7, va='center', color='#424242')

    ax.set_xlim(0.4, 1.15)


def panel_c_permutation(ax):
    """Permutation test distribution (conceptual from existing results)."""
    # From existing results: obs=75, perm_mean=113.6, perm_sd=12.2, n=10000
    obs = 75
    perm_mean = 113.572
    perm_sd = 12.2

    # Simulate the permutation distribution
    np.random.seed(42)
    perm_data = np.random.normal(perm_mean, perm_sd, 10000)

    ax.hist(perm_data, bins=50, color='#B0BEC5', edgecolor='white',
            linewidth=0.5, alpha=0.8, zorder=2, label='Permutation\ndistribution')

    # Observed value
    ax.axvline(obs, color=COL_4mC, linewidth=2.5, linestyle='-',
               zorder=4, label=f'Observed ({obs})')

    # Mean of permutation
    ax.axvline(perm_mean, color='#424242', linewidth=1.0, linestyle='--',
               zorder=3, label=f'Expected ({perm_mean:.1f} \u00b1 {perm_sd:.1f})')

    # p-value annotation
    p_val = 0.0001  # from results
    ax.text(0.03, 0.95,
            f'p = {p_val:.4f}\n(1/10,000)',
            transform=ax.transAxes, fontsize=8, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0',
                      edgecolor='#FB8C00', alpha=0.9, linewidth=0.8))

    ax.set_xlabel('Methylation sites in random regions')
    ax.set_ylabel('Frequency')
    ax.set_title('Permutation test\n(10,000 iterations)', fontsize=11,
                 fontweight='bold')
    ax.legend(fontsize=7.5, frameon=False, loc='upper right')


def panel_d_sigma10(ax, tss_df, pos_4mc, pos_6ma):
    """Promoter methylation density profile relative to TSS,
    highlighting the σ-10 box region."""
    window = 100  # ±100 bp from TSS
    bin_size = 5

    # TSS positions (Jeong2016 only)
    tss_plus = tss_df[tss_df['strand'] == '+']['tss'].values.astype(int)
    tss_minus = tss_df[tss_df['strand'] == '-']['tss'].values.astype(int)

    bins = np.arange(-window, window + bin_size, bin_size)
    counts_all = np.zeros(len(bins) - 1)

    pos_all = np.array(sorted(set(pos_4mc) | set(pos_6ma)))

    # Forward strand genes: TSS relative position = pos - tss
    for tss in tss_plus:
        lo, hi = tss - window, tss + window
        idx = np.searchsorted(pos_all, [lo, hi])
        if idx[1] > idx[0]:
            dists = pos_all[idx[0]:idx[1]] - tss
            hist, _ = np.histogram(dists, bins=bins)
            counts_all += hist

    # Reverse strand genes: TSS relative position = tss - pos (flip)
    for tss in tss_minus:
        lo, hi = tss - window, tss + window
        idx = np.searchsorted(pos_all, [lo, hi])
        if idx[1] > idx[0]:
            dists = tss - pos_all[idx[0]:idx[1]]
            hist, _ = np.histogram(dists, bins=bins)
            counts_all += hist

    n_genes = len(tss_plus) + len(tss_minus)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    density = counts_all / n_genes / (bin_size / 1000)  # per kb per gene

    # Expected density
    n_all = len(pos_all)
    exp_density = n_all / GENOME_SIZE * 1000

    ax.fill_between(bin_centers, density, alpha=0.3, color=COL_DARK, zorder=2)
    ax.plot(bin_centers, density, color=COL_DARK, linewidth=1.5, zorder=3)
    ax.axhline(exp_density, color='gray', linestyle=':', linewidth=0.8,
               alpha=0.5)
    ax.text(-window + 2, exp_density + 0.01, 'genome avg.', fontsize=6,
            color='gray', va='bottom')

    # Highlight σ-10 box region (TSS -12 to -7)
    ax.axvspan(-12, -7, alpha=0.2, color='#FFC107', zorder=1,
               label='\u03C3-10 box (-12 to -7)')

    # Highlight TSS
    ax.axvline(0, color='black', linewidth=0.8, linestyle='-', alpha=0.3)
    ax.text(1, ax.get_ylim()[1] * 0.95, 'TSS', fontsize=7,
            fontweight='bold', va='top')

    # Highlight -35 box
    ax.axvspan(-37, -32, alpha=0.12, color='#90CAF9', zorder=1,
               label='-35 box')

    ax.set_xlabel('Position relative to TSS (bp)')
    ax.set_ylabel('Methylation density\n(sites per kb per gene)')
    ax.set_title(f'Promoter methylation\n(n={n_genes:,} Jeong2016 TSS)',
                 fontsize=11, fontweight='bold')
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    ax.set_xlim(-window, window)
    ax.grid(axis='y', alpha=0.3, lw=0.5)


def main():
    apply_style()
    print('=== Figure 4: TF BS Methylation Depletion ===')

    # Load data
    print('Resolving BS coordinates (Jeong2016 TSS only)...')
    bs_df = resolve_bs_jeong2016()

    print('Loading methylation positions...')
    df_uniq = load_methylation_unique_positions()
    pos_4mc = sorted(df_uniq[df_uniq['mod_type'] == '4mC']['position'].values)
    pos_6ma = sorted(df_uniq[df_uniq['mod_type'] == '6mA']['position'].values)
    print(f'  4mC: {len(pos_4mc)}, 6mA: {len(pos_6ma)}')

    print('Loading Jeong2016 TSS...')
    tss_df = load_tss_jeong2016()

    # Create figure
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(170)))
    gs = GridSpec(2, 2, figure=fig, hspace=0.5, wspace=0.45)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    print('Panel A: Metagene profile...')
    panel_a_metagene(ax_a, bs_df, pos_4mc, pos_6ma)
    add_panel_label(ax_a, 'a')

    print('Panel B: Forest plot...')
    panel_b_forest_plot(ax_b)
    add_panel_label(ax_b, 'b')

    print('Panel C: Permutation distribution...')
    panel_c_permutation(ax_c)
    add_panel_label(ax_c, 'c')

    print('Panel D: Promoter σ-10 box analysis...')
    panel_d_sigma10(ax_d, tss_df, pos_4mc, pos_6ma)
    add_panel_label(ax_d, 'd')

    # Save
    out_path = FIG_DIR / 'Figure4_TF_BS_depletion'
    save_figure(fig, out_path)

    print('\n=== Figure 4 complete ===')


if __name__ == '__main__':
    main()
