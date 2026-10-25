# SUPERSEDED 2026-10-15: このスクリプトは figure_manifest.tsv に Figure1.png の産出元として
# 記録されていないのに Figure1.png を書き出す（現行の出荷図を上書きする）。
# 同型の事故を 2026-10-15 に 1 件起こしている（候補図が既存の補足図を上書きし、
# md5 照合のうえバックアップから復元した）。実行を止めてある。
# 内容を再利用する場合は figure_manifest.tsv に産出元として登録してから解除すること。
raise SystemExit(
    "SUPERSEDED 2026-10-15: 01_figure1_landscape.py は産出元として登録されていないため実行を止めています。"
    "Figure1.png の現行の産出元は figure_manifest.tsv を参照してください。"
)

#!/usr/bin/env python3
"""
New Figure 1: Methylation Landscape of S. coelicolor A3(2) M145
  Panel A: Linear genome ideogram with 4mC/6mA density per timepoint
  Panel B: HC site counts by timepoint (true counts, not deduplicated)
  Panel C: Genomic region distribution
  Panel D: Sequence logos (logomaker)
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
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
import logomaker

# ── Unified colour-blind-safe palette (Okabe-Ito) ────────────────────────
# Override shared-utils shades locally for this figure only, so that
# 4mC / 6mA / genome are consistent across panels a, b, c.
COL_4mC = '#A64B44'   # unified rich muted red — 4mC (matches all figs + shared_utils)
COL_6mA = '#3A6B8C'   # unified rich muted blue — 6mA
COL_GRAY = '#9AA7B0'  # unified neutral grey — genome/baseline series

# ── Genome landmarks ─────────────────────────────────────────────────────
GENOME_LEN = 8_667_507
ARM_LEFT = 1_500_000
ARM_RIGHT = 7_167_507
ORIC_POS = 4_270_777       # dnaA (SC_RS21540)
TIR_LEN = 21_653           # Terminal Inverted Repeat

COL_CORE = '#3E7256'       # unified deep muted green — chromosomal core (bands drawn translucent)
COL_ARM = '#C0803A'        # unified calm amber — chromosomal arm


def panel_a_linear(ax):
    """Linear genome ideogram with 4mC and 6mA density per timepoint."""
    # Load HC sites (all timepoints, including overlaps)
    df_hc = load_methylation_hc_all()
    df_4mc = df_hc[df_hc['mod_type'] == '4mC']
    df_6ma = df_hc[df_hc['mod_type'] == '6mA']

    timepoints = ['T1', 'T2', 'T3']
    tp_labels = TP_LABELS
    y_positions = [3.5, 2.0, 0.5]  # Y center for each timepoint row
    chrom_h = 0.18
    bar_h = 0.45  # height for density bars

    genome_mb = GENOME_LEN / 1e6
    bin_size_bp = 50_000
    bin_size_mb = bin_size_bp / 1e6
    bins_bp = np.arange(0, GENOME_LEN + bin_size_bp, bin_size_bp)
    bin_centers_mb = (bins_bp[:-1] + bins_bp[1:]) / 2 / 1e6

    # Compute global max density for consistent scaling
    max_dens = 0
    for tp in timepoints:
        for df_mod in [df_4mc, df_6ma]:
            sub = df_mod[df_mod['timepoint'] == tp]
            if len(sub) > 0:
                h, _ = np.histogram(sub['position'].values, bins=bins_bp)
                d = h / (bin_size_bp / 1000)  # sites per kb
                max_dens = max(max_dens, d.max())
    if max_dens == 0:
        max_dens = 1

    for tp, tp_label, y in zip(timepoints, tp_labels, y_positions):
        # ── Chromosome rectangle ──
        # Left arm
        ax.add_patch(Rectangle((0, y - chrom_h / 2), ARM_LEFT / 1e6, chrom_h,
                               facecolor=COL_ARM, alpha=0.25, edgecolor='none'))
        # Core
        ax.add_patch(Rectangle((ARM_LEFT / 1e6, y - chrom_h / 2),
                               (ARM_RIGHT - ARM_LEFT) / 1e6, chrom_h,
                               facecolor=COL_CORE, alpha=0.12, edgecolor='none'))
        # Right arm
        ax.add_patch(Rectangle((ARM_RIGHT / 1e6, y - chrom_h / 2),
                               (GENOME_LEN - ARM_RIGHT) / 1e6, chrom_h,
                               facecolor=COL_ARM, alpha=0.25, edgecolor='none'))
        # Chromosome outline (linear ends, not circular)
        ax.plot([0, genome_mb], [y - chrom_h / 2, y - chrom_h / 2],
                color='#37474F', linewidth=0.8)
        ax.plot([0, genome_mb], [y + chrom_h / 2, y + chrom_h / 2],
                color='#37474F', linewidth=0.8)
        # Squared ends (linear genome)
        ax.plot([0, 0], [y - chrom_h / 2, y + chrom_h / 2],
                color='#37474F', linewidth=1.2)
        ax.plot([genome_mb, genome_mb], [y - chrom_h / 2, y + chrom_h / 2],
                color='#37474F', linewidth=1.2)

        # ── 4mC density (above chromosome) ──
        sub_4mc = df_4mc[df_4mc['timepoint'] == tp]
        n_4mc = len(sub_4mc)
        if n_4mc > 0:
            h4, _ = np.histogram(sub_4mc['position'].values, bins=bins_bp)
            d4 = h4 / (bin_size_bp / 1000)
            heights = d4 / max_dens * bar_h
            ax.bar(bin_centers_mb, heights, width=bin_size_mb,
                   bottom=y + chrom_h / 2,
                   color=COL_4mC, alpha=0.7, edgecolor='none')

        # ── 6mA density (below chromosome, inverted) ──
        sub_6ma = df_6ma[df_6ma['timepoint'] == tp]
        n_6ma = len(sub_6ma)
        if n_6ma > 0:
            h6, _ = np.histogram(sub_6ma['position'].values, bins=bins_bp)
            d6 = h6 / (bin_size_bp / 1000)
            heights = d6 / max_dens * bar_h
            ax.bar(bin_centers_mb, -heights, width=bin_size_mb,
                   bottom=y - chrom_h / 2,
                   color=COL_6mA, alpha=0.7, edgecolor='none')

        # ── Labels ──
        ax.text(-0.4, y + 0.15, tp_label, ha='right', va='center',
                fontsize=8, fontweight='bold')
        ax.text(-0.4, y - 0.15,
                f'4mC: {n_4mc:,}  6mA: {n_6ma:,}',
                ha='right', va='center', fontsize=6.5, color='#555B61')

    # ── Genome landmarks ──
    # Arm/core boundaries
    for bnd in [ARM_LEFT, ARM_RIGHT]:
        ax.axvline(bnd / 1e6, color='gray', linestyle=':', linewidth=0.5,
                   ymin=0.02, ymax=0.98)

    # oriC marker (recoloured grey/black per reviewer C4 — reduce palette;
    # colour is reserved for the 4mC/6mA data series only)
    ax.axvline(ORIC_POS / 1e6, color='#333333', linestyle='-', linewidth=1.2,
               ymin=0.02, ymax=0.98, alpha=0.8)
    # White bbox so the oriC line does not strike through the letters.
    ax.text(ORIC_POS / 1e6, y_positions[0] + bar_h + chrom_h / 2 + 0.25,
            'oriC', ha='center', va='bottom', fontsize=7,
            fontweight='bold', color='#333333',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor='none', alpha=0.85))

    # TIR markers
    for tir_start, tir_end, label in [
        (0, TIR_LEN, 'TIR'),
        (GENOME_LEN - TIR_LEN, GENOME_LEN, 'TIR'),
    ]:
        for _, _, y in zip(timepoints, tp_labels, y_positions):
            ax.add_patch(Rectangle(
                (tir_start / 1e6, y - chrom_h / 2),
                (tir_end - tir_start) / 1e6, chrom_h,
                facecolor='#E65100', alpha=0.35, edgecolor='none', zorder=5))
    # TIR labels (once, at top)
    ax.text(TIR_LEN / 2 / 1e6, y_positions[0] + bar_h + chrom_h / 2 + 0.25,
            'TIR', ha='center', va='bottom', fontsize=6, color='#E65100',
            fontstyle='italic')
    ax.text((GENOME_LEN - TIR_LEN / 2) / 1e6,
            y_positions[0] + bar_h + chrom_h / 2 + 0.25,
            'TIR', ha='center', va='bottom', fontsize=6, color='#E65100',
            fontstyle='italic')

    # Region labels at bottom (reviewer R2-08). "Core" is centred on the TRUE
    # core-span midpoint (~4.33 Mb) so it sits symmetrically between the two
    # arm/core boundaries. Because that midpoint is ~under the oriC line
    # (4.27 Mb), a white bbox (as on the oriC label) keeps the line from
    # striking through the word rather than shifting the label off-centre.
    CORE_MID = (ARM_LEFT + ARM_RIGHT) / 2 / 1e6
    ax.text(ARM_LEFT / 2 / 1e6, -0.3, 'Left arm', ha='center', va='center',
            fontsize=7, color=COL_ARM, fontstyle='italic')
    ax.text(CORE_MID, -0.3, 'Core',
            ha='center', va='center', fontsize=7, color=COL_CORE,
            fontstyle='italic',
            bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.9))
    ax.text((ARM_RIGHT + GENOME_LEN) / 2 / 1e6, -0.3, 'Right arm',
            ha='center', va='center', fontsize=7, color=COL_ARM,
            fontstyle='italic')

    # Numeric core/arm boundary definition on-panel (reviewer C4): the core
    # spans the two dotted arm/core boundaries; state the coordinates so the
    # definition is legible without the legend. Placed under the (left-shifted)
    # "Core" label so it clears the oriC line at ~4.27 Mb.
    ax.text(CORE_MID, -0.46,
            f'core {ARM_LEFT/1e6:.2f}\u2013{ARM_RIGHT/1e6:.2f} Mb',
            ha='center', va='center', fontsize=6, color=COL_CORE,
            fontstyle='italic',
            bbox=dict(boxstyle='square,pad=0.1', facecolor='white',
                      edgecolor='none', alpha=0.9))

    # Coordinate axis (Mb)
    ax.set_xlim(-0.5, genome_mb + 0.5)
    ax.set_ylim(-0.6, y_positions[0] + bar_h + chrom_h / 2 + 0.5)
    mb_ticks = np.arange(0, genome_mb + 0.5, 1)
    ax.set_xticks(mb_ticks)
    ax.set_xticklabels([f'{v:.0f}' for v in mb_ticks], fontsize=7)
    ax.set_xlabel('Chromosome position (Mb)', fontsize=9)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)

    # Legend
    legend_elements = [
        Line2D([0], [0], color=COL_4mC, lw=6, alpha=0.7, label='4mC'),
        Line2D([0], [0], color=COL_6mA, lw=6, alpha=0.7, label='6mA'),
        mpatches.Patch(facecolor=COL_CORE, alpha=0.15, label='Core'),
        mpatches.Patch(facecolor=COL_ARM, alpha=0.3, label='Arm'),
    ]
    ax.legend(handles=legend_elements, fontsize=7, loc='upper right',
              frameon=True, framealpha=0.9, edgecolor='#CCC', ncol=4)

    ax.set_title('Linear chromosome methylation landscape '
                 '(S. coelicolor A3(2) M145, 8.67 Mb)',
                 fontsize=10, fontweight='bold')


def panel_b_site_counts(ax):
    """HC methylation site counts by timepoint (true counts, not deduplicated)."""
    df_hc = load_methylation_hc_all()

    counts_4mc = (df_hc[df_hc['mod_type'] == '4mC']
                  .groupby('timepoint').size()
                  .reindex(['T1', 'T2', 'T3'], fill_value=0))
    counts_6ma = (df_hc[df_hc['mod_type'] == '6mA']
                  .groupby('timepoint').size()
                  .reindex(['T1', 'T2', 'T3'], fill_value=0))

    x = np.arange(3)
    width = 0.35

    bars1 = ax.bar(x - width / 2, counts_4mc.values, width, color=COL_4mC,
                   edgecolor='white', linewidth=0.5, label='4mC', alpha=0.9,
                   zorder=3)
    bars2 = ax.bar(x + width / 2, counts_6ma.values, width, color=COL_6mA,
                   edgecolor='white', linewidth=0.5, label='6mA', alpha=0.9,
                   zorder=3)

    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, h + 30,
                        f'{int(h):,}', ha='center', va='bottom', fontsize=6.5,
                        fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(TP_LABELS_NL)
    ax.set_ylabel('HC methylation sites')
    ax.set_title('Site counts by timepoint', fontsize=11, fontweight='bold')
    # Legend at upper-left: the tall T2/T3 bars and their value labels occupy the
    # upper-right, so an upper-right legend collided with the '2,446'/'2,295' labels.
    ax.legend(fontsize=8, frameon=False, loc='upper left')
    ax.set_ylim(0, max(counts_4mc.max(), counts_6ma.max()) * 1.28)
    ax.grid(axis='y', alpha=0.3, lw=0.5)

    # E9: "Unique position" annotation removed from panel B
    # (unique position counts reported in text/legend of figure legend instead)


def panel_c_genomic_distribution(ax, ax_inset=None):
    """Genomic region distribution of methylation sites.

    Reviewer A1 fix: CDS dominates (~80%), so non-CDS differences are
    invisible on a linear scale. Added a zoomed inset (0-15% range)
    that highlights Promoter / 5'UTR / Intergenic comparisons.
    """
    df_tss = load_tss_jeong2016()
    gbk = load_reference_gbk()
    genome = load_reference_genome()
    genome_len = len(genome.seq)
    df_unique = load_methylation_unique_positions()

    # Assign each genome position a category
    # 0=Promoter, 1=5'UTR, 2=CDS, 3=Intergenic
    region = np.full(genome_len, 3, dtype=np.int8)

    # CDS from GenBank
    for feat in gbk.features:
        if feat.type == 'CDS':
            s = int(feat.location.start)
            e = int(feat.location.end)
            region[s:e] = 2

    # Promoter and 5'UTR from Jeong2016 TSS (overrides CDS)
    for _, row in df_tss.iterrows():
        tss = int(row['tss'])
        strand = row['strand']
        if strand == '+':
            p_start = max(0, tss - 300)
            p_end = max(0, tss)
            region[p_start:p_end] = 0
            u_end = min(genome_len, tss + 50)
            region[tss:u_end] = 1
        else:
            p_start = min(genome_len, tss + 1)
            p_end = min(genome_len, tss + 301)
            region[p_start:p_end] = 0
            u_start = max(0, tss - 49)
            region[u_start:tss + 1] = 1

    categories = ['Promoter\n(−300 to TSS)', "5′ UTR\n(TSS to +50)",
                  'CDS', 'Intergenic']

    from scipy.stats import chi2_contingency
    from statsmodels.stats.multitest import multipletests

    raw_counts = {}
    results = {}
    for mod_type in ['4mC', '6mA']:
        positions = df_unique[df_unique['mod_type'] == mod_type]['position'].values
        positions = positions[(positions >= 0) & (positions < genome_len)]
        site_regions = region[positions]
        counts = np.bincount(site_regions, minlength=4)
        raw_counts[mod_type] = counts
        results[mod_type] = counts / counts.sum() * 100

    # Genome background
    genome_counts = np.bincount(region, minlength=4)
    results['Genome'] = genome_counts / genome_counts.sum() * 100
    genome_prop = genome_counts / genome_counts.sum()

    # ── Statistical testing (pre-specified, anti-p-hacking protocol) ──────────
    # Step 1: Global chi-square goodness-of-fit for each mod type vs genome
    #   H0: methylation sites distributed proportionally to genome composition
    #   Two global tests → Bonferroni α = 0.05/2 = 0.025
    sig_markers = {}  # {(mod_idx, region_idx): marker}
    for mod_type in ['4mC', '6mA']:
        obs = raw_counts[mod_type]
        exp = genome_prop * obs.sum()
        chi2, p_global, dof, _ = chi2_contingency(
            np.array([obs, exp.astype(int)]))
        print(f'{mod_type} global chi2={chi2:.1f}, p={p_global:.2e}, dof={dof}')

        if p_global < 0.025:   # Bonferroni threshold for 2 global tests
            # Step 2: per-region proportion z-test (post-hoc, only if global sig)
            # Compare each region's observed proportion vs genome proportion
            # Bonferroni correction for 4 regions within each mod type
            pvals = []
            for j in range(4):
                n_total = obs.sum()
                k = obs[j]
                p0 = genome_prop[j]
                # Two-proportion z-test (vs genome background)
                from statsmodels.stats.proportion import proportions_ztest
                stat, p = proportions_ztest(k, n_total, p0)
                pvals.append(p)
            reject, pvals_adj, _, _ = multipletests(pvals, method='bonferroni')
            for j, (rej, p_adj) in enumerate(zip(reject, pvals_adj)):
                if rej:
                    obs_pct = obs[j] / obs.sum() * 100
                    exp_pct = genome_prop[j] * 100
                    direction = '▲' if obs_pct > exp_pct else '▼'
                    marker = ('***' if p_adj < 0.001 else
                              '**'  if p_adj < 0.01  else '*')
                    sig_markers[(mod_type, j)] = (direction, marker)
                    print(f'  {categories[j]}: {direction} p_adj={p_adj:.3e} {marker}')

    # Plot
    x = np.arange(4)
    width = 0.25

    ax.bar(x - width, results['4mC'], width, color=COL_4mC,
           edgecolor='white', linewidth=0.5, alpha=0.9, label='4mC', zorder=3)
    ax.bar(x, results['6mA'], width, color=COL_6mA,
           edgecolor='white', linewidth=0.5, alpha=0.9, label='6mA', zorder=3)
    ax.bar(x + width, results['Genome'], width, color=COL_GRAY,
           edgecolor='white', linewidth=0.5, alpha=0.9, label='Genome', zorder=3)

    # Value labels
    offsets = [-width, 0, width]
    labels_data = [('4mC', results['4mC']), ('6mA', results['6mA']),
                   ('Genome', results['Genome'])]
    for i, (name, cat_data) in enumerate(labels_data):
        for j in range(4):
            val = cat_data[j]
            if val > 1.5 and j != 2:
                ax.text(x[j] + offsets[i], val + 0.4,
                        f'{val:.1f}', ha='center', va='bottom',
                        fontsize=5.5)

    # Significance markers (C5 fix): anchor every marker in a region to a common
    # ceiling above the TALLEST bar of that region, so a marker over a short bar
    # (e.g. 4mC Intergenic = 5.8) no longer lands at the height of a taller
    # neighbour's value label (e.g. 6mA Intergenic = 10.7). 4mC and 6mA markers
    # stay at their own x but share the region ceiling height.
    bar_x_map = {'4mC': x - width, '6mA': x}
    region_ceiling = {
        j: max(results['4mC'][j], results['6mA'][j], results['Genome'][j])
        for j in range(4)
    }
    for (mod_type, j), (direction, marker) in sig_markers.items():
        bx = bar_x_map[mod_type][j]
        color = COL_4mC if mod_type == '4mC' else COL_6mA
        ax.text(bx, region_ceiling[j] + 2.5, f'{direction}{marker}',
                ha='center', va='bottom', fontsize=7,
                color=color, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=7.5)
    ax.set_ylabel('Proportion (%)')
    # Explicit headroom so ▲▼ markers stay inside the axes (no spill below).
    ax.set_ylim(0, max(max(results['4mC']), max(results['6mA']),
                       max(results['Genome'])) * 1.15)
    ax.set_title('Genomic region distribution', fontsize=11, fontweight='bold',
                 pad=8)
    # Only the CDS bar is tall (~85%); promoter/5′UTR/intergenic bars are short,
    # so the upper-left interior is empty. Put the colour key there — clear of
    # the bars, the title, and the ▲▼ significance markers.
    ax.legend(fontsize=7, frameon=False, loc='upper left',
              bbox_to_anchor=(0.0, 0.98), handletextpad=0.4,
              labelspacing=0.3)
    ax.grid(axis='y', alpha=0.3, lw=0.5)
    # (▲/▼ significance-test explanation moved to the figure legend, 2026-06-16,
    #  per reviewer: figures should not carry explanatory text directly.)

    # Non-CDS inset REMOVED 2026-06-16 (redundant with Supplementary Fig 9 occupancy
    #  pie charts, and it overlapped panel c). Inset code below disabled.
    return
    # ── Reviewer A1: Inset zoom for non-CDS categories ─────────────────────
    if ax_inset is None:
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        ax_inset = inset_axes(ax, width='38%', height='42%',
                              bbox_to_anchor=(0.18, 0.45, 1.0, 1.0),
                              bbox_transform=ax.transAxes,
                              loc='upper left', borderpad=0.4)

    nonCDS = [0, 1, 3]      # Promoter, 5'UTR, Intergenic indices
    x_nc = np.arange(len(nonCDS))
    width_nc = 0.25
    ax_inset.bar(x_nc - width_nc, results['4mC'][nonCDS], width_nc,
                 color=COL_4mC, alpha=0.9, edgecolor='white', linewidth=0.4)
    ax_inset.bar(x_nc, results['6mA'][nonCDS], width_nc,
                 color=COL_6mA, alpha=0.9, edgecolor='white', linewidth=0.4)
    ax_inset.bar(x_nc + width_nc, results['Genome'][nonCDS], width_nc,
                 color=COL_GRAY, alpha=0.9, edgecolor='white', linewidth=0.4)

    nc_max = max(max(results['4mC'][nonCDS]),
                 max(results['6mA'][nonCDS]),
                 max(results['Genome'][nonCDS]))
    ax_inset.set_ylim(0, nc_max * 1.30)
    ax_inset.set_xticks(x_nc)
    ax_inset.set_xticklabels(['Prom.', "5′UTR", 'Interg.'], fontsize=6)
    ax_inset.tick_params(axis='y', labelsize=6, length=2, pad=1)
    ax_inset.set_ylabel('%', fontsize=6, labelpad=1)
    ax_inset.set_title('Non-CDS detail', fontsize=6.5, pad=2)
    ax_inset.grid(axis='y', alpha=0.25, lw=0.4)
    ax_inset.spines['top'].set_visible(False)
    ax_inset.spines['right'].set_visible(False)

    # Re-mark significance on inset for non-CDS bars
    for (mod_type, j), (direction, marker) in sig_markers.items():
        if j in nonCDS:
            local_idx = nonCDS.index(j)
            offset = -width_nc if mod_type == '4mC' else 0
            bar_top = results[mod_type][j] + 0.6
            color = COL_4mC if mod_type == '4mC' else COL_6mA
            ax_inset.text(local_idx + offset, bar_top,
                          f'{direction}{marker}',
                          ha='center', va='bottom', fontsize=5,
                          color=color, fontweight='bold')

    # Value labels on inset bars
    for ix, ridx in enumerate(nonCDS):
        for off, name in zip([-width_nc, 0, width_nc],
                             ['4mC', '6mA', 'Genome']):
            v = results[name][ridx]
            if v > 0:
                ax_inset.text(ix + off, v + nc_max * 0.025,
                              f'{v:.1f}', ha='center', va='bottom',
                              fontsize=5, color='#444')


def panel_d_logos(ax_top, ax_bot):
    """Sequence logos for top MEME motifs."""
    meme_4mc = EPIGENOME / 'archive' / 'v1_weighted_minreps2' / \
        '07_motif_analysis' / 'meme_4mC' / 'meme.txt'
    meme_6ma = EPIGENOME / 'archive' / 'v1_weighted_minreps2' / \
        '07_motif_analysis' / 'meme_6mA' / 'meme.txt'

    motifs_4mc = parse_meme_pwm(meme_4mc)
    motifs_6ma = parse_meme_pwm(meme_6ma)

    nuc_colors = {'A': '#43A047', 'C': '#1565C0',
                  'G': '#FFA000', 'T': '#E53935'}

    # ── 4mC logo ──
    pwm_4mc = motifs_4mc['MEME-1']['pwm']
    info_4mc = logomaker.transform_matrix(
        pwm_4mc, from_type='probability', to_type='information')
    logomaker.Logo(info_4mc, ax=ax_top, color_scheme=nuc_colors)
    ax_top.set_ylabel('Bits', fontsize=8)
    m4 = motifs_4mc['MEME-1']
    ax_top.set_title(
        f'4mC: {m4["name"]}  (n={m4["nsites"]:,}, E={m4["evalue"]})',
        fontsize=8, fontweight='bold', color='black')
    ax_top.set_ylim(0, 2.2)
    ax_top.set_xlim(-0.5, m4['width'] - 0.5)
    ax_top.set_xticklabels([])  # Hide x-axis labels on top logo

    # ── 6mA logo ──
    pwm_6ma = motifs_6ma['MEME-1']['pwm']
    info_6ma = logomaker.transform_matrix(
        pwm_6ma, from_type='probability', to_type='information')
    logomaker.Logo(info_6ma, ax=ax_bot, color_scheme=nuc_colors)
    ax_bot.set_ylabel('Bits', fontsize=8)
    m6 = motifs_6ma['MEME-1']
    ax_bot.set_title(
        f'6mA: {m6["name"]}  (n={m6["nsites"]:,}, E={m6["evalue"]})',
        fontsize=8, fontweight='bold', color='black')
    ax_bot.set_ylim(0, 2.2)
    ax_bot.set_xlim(-0.5, m6['width'] - 0.5)


def main():
    apply_style()
    print('=== New Figure 1: Methylation Landscape ===')

    # Layout: 3 rows
    #   Row 1: Panel A (linear genome, full width)
    #   Row 2: Panel B (site counts, left) + Panel C (genomic dist, right)
    #   Row 3: Panel D (logos, full width — two sub-axes)
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(250)))

    # Panel A: Linear genome (top, full width)
    ax_a = fig.add_axes([0.17, 0.65, 0.80, 0.31])

    # Panel B: Site counts (middle left). Lowered so panel A's x-axis label
    # ("Chromosome position (Mb)") clears the panel-B/C titles and the "c" label.
    ax_b = fig.add_axes([0.10, 0.35, 0.37, 0.19])

    # Panel C: Genomic distribution (middle right)
    ax_c = fig.add_axes([0.58, 0.35, 0.38, 0.19])

    # Panel D: Logos (bottom, two sub-axes).
    # Extra vertical gap below panel B/C (0.36 → 0.205) keeps panel C's
    # x-tick labels and ▲▼ markers clear of panel D's title / "d" label.
    # 2026-09-13 (FIG-16): the lower logo sat at y=0.015, so its x-tick labels
    # fell off the canvas bottom; both logo axes lifted by 0.025.
    ax_d1 = fig.add_axes([0.10, 0.150, 0.85, 0.08])
    ax_d2 = fig.add_axes([0.10, 0.040, 0.85, 0.08])

    print('Panel A: Linear genome ideogram...')
    panel_a_linear(ax_a)
    add_panel_label(ax_a, 'a', x=-0.08, y=1.08)

    print('Panel B: Site counts by timepoint...')
    panel_b_site_counts(ax_b)
    add_panel_label(ax_b, 'b', x=-0.15, y=1.10)

    print('Panel C: Genomic region distribution...')
    panel_c_genomic_distribution(ax_c)
    add_panel_label(ax_c, 'c', x=-0.12, y=1.10)

    print('Panel D: Sequence logos...')
    panel_d_logos(ax_d1, ax_d2)
    add_panel_label(ax_d1, 'd', x=-0.06, y=1.25)

    # Save (PNG added for manuscript slot / proofing alongside PDF+SVG)
    out_path = FIG_DIR / 'new_Figure1_methylation_landscape'
    save_figure(fig, out_path, formats=('pdf', 'svg', 'png'))

    # sync into the Obsidian manuscript slot (Figure1.png) — the script previously
    # left this slot stale (7.14" old version), so add an explicit copy.
    import shutil
    # 2026-09-13: sync is opt-in (SYNC_FIG_SLOTS=1); back up the old PNG to
    # fig_images/archive/ first.
    import os
    slot = Path.home() / 'obsidian' / 'Research' / 'rna-seq' / 'Writing' / 'fig_images' / 'Figure1.png'
    if os.environ.get('SYNC_FIG_SLOTS') == '1' and slot.parent.is_dir():
        shutil.copyfile(out_path.with_suffix('.png'), slot)
        print(f'  Synced → {slot}')

    print('\n=== New Figure 1 complete ===')


if __name__ == '__main__':
    main()
