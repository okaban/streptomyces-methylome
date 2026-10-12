#!/usr/bin/env python3
"""
FigS13: Comprehensive Per-motif Analysis

Panel A: Per-motif methylation frequency × downstream gene expression correlation
         (median LFC by motif group: AAGCCCG-only / GCCGGC-only / Dual / Background)
Panel B: Timepoint distribution variation — site counts per motif × timepoint
Panel C: Genomic region distribution — proportion of sites in Promoter / 5'-UTR / CDS / Intergenic
         per motif candidate group (stacked 100% bar chart)
Panel D: TSS protection zone by motif type (spatial profiles for regulatory genes)
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
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy import stats

DIR_40 = EPIGENOME / '40_motif_division_of_labor' / 'tables'
DIR_42 = EPIGENOME / '42_GCCGGC_temporal_derepression' / 'tables'
DIR_44 = EPIGENOME / '44_AAGCCCG_temporal_causality' / 'tables'
DIR_48 = EPIGENOME / '48_TSS_methylation_gradient' / 'tables'
CENSUS_4mC = EPIGENOME / '23_expanded_motif_search' / '4mC_final_census_v1.csv'
CENSUS_6mA = EPIGENOME / '23_expanded_motif_search' / '6mA_final_census_v1.csv'
GFF_PATH = Path('/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff')
GENE_ANN_PATH = BASE / '05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv'
TSS_PATH = EPIGENOME / '18_tss_analyses' / 'comprehensive_tss_table.csv'
JEONG_TSS_PATH = EPIGENOME / '18_tss_analyses' / 'jeong2016_all_tss.csv'
CHROM = 'NC_003888.3'
CHROM_LEN = 8_667_507

# ── Motif group definitions ───────────────────────────────────────────────────
MOTIF_GROUPS = [
    ('GCCGGC\n(4mC)',   '4mC', ['TGGCCGGC', 'GGCCGG', 'CCGG (other)']),
    ('AAGCCCG\n(4mC)',  '4mC', ['AAGCCCG']),
    ('AAGCCCG\n(6mA)',  '6mA', ['AAGCCCG']),
    ('CCGKCA\n(6mA)',   '6mA', ['CCGKCA', 'GCCG', 'CCGC', 'CCGG', 'CCSGG']),
    ('GATC\n(6mA)',     '6mA', ['GATC']),
]
REGION_ORDER  = ['Promoter', "5'-UTR", 'CDS', 'Intergenic']
REGION_COLORS = {
    'Promoter':    '#E53935',
    "5'-UTR":      '#FB8C00',
    'CDS':         '#43A047',
    'Intergenic':  '#B0BEC5',
}
REGION_CODE = {r: i for i, r in enumerate(REGION_ORDER)}
PROMOTER_UP = 300  # bp upstream of an experimental TSS (unified definition)

GROUP_COLORS = {
    'AAGCCCG-only':  COL_BOTH,
    'GCCGGC-only':   COL_4mC,
    'Dual-targeted': '#8E24AA',
    'Background':    '#9E9E9E',
}

def _style(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── Panel A: Expression correlation by motif group ────────────────────────────

def panel_a(ax, df_expr):
    """Grouped bars: median LFC at T2 and T3 by motif group."""
    groups = ['AAGCCCG-only', 'GCCGGC-only', 'Dual-targeted', 'Background']
    lfc_t2 = [df_expr.loc[df_expr['group'] == g, 'median_LFC'].values[0]
               if 'median_LFC' in df_expr.columns
               else df_expr.loc[df_expr['group'] == g, 'median_LFC_T2vsT1'].values[0]
               for g in groups]
    lfc_t3 = [df_expr.loc[df_expr['group'] == g, 'median_LFC_T3vsT1'].values[0]
               for g in groups]

    x = np.arange(len(groups))
    w = 0.38
    colors = [GROUP_COLORS[g] for g in groups]

    bars2 = ax.bar(x - w/2, lfc_t2, width=w, color=colors, alpha=0.85,
                   edgecolor='white', linewidth=0.5, label='T2 vs T1')
    bars3 = ax.bar(x + w/2, lfc_t3, width=w, color=colors, alpha=0.45,
                   edgecolor='white', linewidth=0.5, hatch='///',
                   label='T3 vs T1')

    ax.axhline(0, color='gray', linewidth=0.8, linestyle='--')
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=7)
    ax.set_ylabel('Median log₂FC', fontsize=8)
    ax.set_title('Methylation group × expression change', fontsize=9, fontweight='bold')

    # Custom legend (solid=T2, hatched=T3)
    leg_patches = [
        mpatches.Patch(facecolor='#888', edgecolor='none', label='T2 vs T1 (median LFC)'),
        mpatches.Patch(facecolor='#888', edgecolor='#888', hatch='///',
                       label='T3 vs T1 (median LFC)'),
    ]
    ax.legend(handles=leg_patches, fontsize=6.5, frameon=False, loc='lower right')
    _style(ax)

    # Color patches for group coding
    for xi, g in zip(x, groups):
        ax.text(xi, ax.get_ylim()[0] - 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0]),
                '●', ha='center', va='top', color=GROUP_COLORS[g], fontsize=9)


# ── Panel B: Timepoint distribution ──────────────────────────────────────────

def panel_b(ax, c4, c6):
    """Stacked bars per timepoint × motif."""
    # 4mC motifs grouped
    motif_groups_4mc = {
        'GCCGGC\n(4mC)':  ['TGGCCGGC', 'GGCCGG', 'CCGG (other)'],
        'AAGCCCG\n(4mC)': ['AAGCCCG'],
    }
    motif_groups_6ma = {
        'AAGCCCG\n(6mA)':  ['AAGCCCG'],
        'CCGKCA\n(6mA)':   ['CCGKCA', 'GCCG', 'CCGC', 'CCGG', 'CCSGG'],
        'GATC\n(6mA)':     ['GATC'],
        'Unassigned':      ['GAACCGG', 'CGGCAACC', 'CTGCTCGCCG', 'unassigned'],
    }
    colors_4mc = {
        'GCCGGC\n(4mC)': COL_4mC,
        'AAGCCCG\n(4mC)': COL_BOTH,
    }
    colors_6ma = {
        'AAGCCCG\n(6mA)':  COL_BOTH,
        'CCGKCA\n(6mA)':   '#42A5F5',
        'GATC\n(6mA)':     '#0D47A1',
        'Unassigned':      '#CFD8DC',
    }

    timepoints = ['T1', 'T2', 'T3']
    x = np.arange(len(timepoints))
    width = 0.35

    def count_group(census, motif_vals, tp):
        return ((census['timepoint'] == tp) &
                (census['final_motif'].isin(motif_vals))).sum()

    # 4mC stacked bars (left)
    bottom = np.zeros(3)
    for label, motifs in motif_groups_4mc.items():
        counts = np.array([count_group(c4, motifs, tp) for tp in timepoints])
        ax.bar(x - width/2, counts, width=width, bottom=bottom,
               color=colors_4mc[label], alpha=0.9, edgecolor='white',
               linewidth=0.4, label=label)
        bottom += counts

    # 6mA stacked bars (right)
    bottom = np.zeros(3)
    for label, motifs in motif_groups_6ma.items():
        counts = np.array([count_group(c6, motifs, tp) for tp in timepoints])
        ax.bar(x + width/2, counts, width=width, bottom=bottom,
               color=colors_6ma[label], alpha=0.9, edgecolor='white',
               linewidth=0.4, label=label)
        bottom += counts

    ax.set_xticks(x)
    ax.set_xticklabels(timepoints, fontsize=8)
    ax.set_ylabel('Number of HC sites', fontsize=8)
    ax.set_title('Methylation site count by timepoint\n(left: 4mC, right: 6mA)',
                 fontsize=9, fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=6, frameon=False, loc='upper right',
              ncol=2, handlelength=1.0)
    _style(ax)

    # Text labels above bars
    for i, tp in enumerate(timepoints):
        n4 = len(c4[c4['timepoint'] == tp])
        n6 = len(c6[c6['timepoint'] == tp])
        ax.text(i - width/2, n4 + 20, f'{n4}', ha='center', va='bottom',
                fontsize=6.5, color='#333')
        ax.text(i + width/2, n6 + 20, f'{n6}', ha='center', va='bottom',
                fontsize=6.5, color='#333')


# ── Region annotation helpers ─────────────────────────────────────────────────

def build_region_array(gene_ann, tss_df):
    """Build a per-position region code array for CHROM.

    Unified region definition, shared with Figure 8 and Supp Fig 21
    (region_scheme_reconciliation.csv):

      TSS      : experimental only (Jeong et al. 2016 dRNA-seq). Genes without
                 an experimental TSS contribute a CDS interval and nothing else
                 — no annotated-boundary proxy, because the boundary proxy was
                 shown to behave as a negative control for TSS-proximal signal.
      Priority (highest overwrites): Promoter > 5'-UTR > CDS > Intergenic.
                 The promoter-first direction is the claim-adverse one: it
                 maximises the promoter share, so the paper's "methylation is
                 not concentrated at promoters" reading cannot be an artefact
                 of the category precedence.

    Painting a per-base array makes the four categories mutually exclusive by
    construction; sites that are both promoter-proximal and inside a CDS are
    resolved to Promoter and their count is reported in the legend.
    Returns int8 array of length CHROM_LEN+2 where value = REGION_CODE[region_name].
    """
    arr = np.full(CHROM_LEN + 2, REGION_CODE['Intergenic'], dtype=np.int8)

    g = gene_ann[gene_ann['contig'] == CHROM].copy().reset_index(drop=True)
    # Only PRIMARY ('P') and SECONDARY ('S') TSSs define a gene's promoter.
    # Antisense ('A') TSSs lie on the opposite strand (246/246) and inside the
    # gene body (98.8%); internal ('I') TSSs are inside the gene body (98.7%),
    # so neither one's upstream window is that gene's promoter. Genes with
    # alternative promoters contribute one window per TSS.
    jeong = pd.read_csv(JEONG_TSS_PATH)
    jeong = jeong[jeong['gene_id'].notna() & jeong['category'].isin(('P', 'S'))]
    gidx = {gid: i for i, gid in enumerate(g['gene_id'])}
    starts, ends, strands = g['start'].values, g['end'].values, g['strand'].values
    tss_rows = [(gidx[r.gene_id], int(r.tss_position), strands[gidx[r.gene_id]])
                for r in jeong.itertuples() if r.gene_id in gidx]
    g['has_exp_tss'] = g['gene_id'].isin({r.gene_id for r in jeong.itertuples()})

    # 1. CDS (lowest priority — filled first, overwritten by higher)
    for _, r in g.iterrows():
        arr[int(r['start']):int(r['end']) + 1] = REGION_CODE['CDS']

    # 2. 5'-UTR (between experimental TSS and CDS start)
    for i, tss, strand in tss_rows:
        if strand == '+' and tss < starts[i]:
            arr[tss:int(starts[i])] = REGION_CODE["5'-UTR"]
        elif strand == '-' and tss > ends[i]:
            arr[int(ends[i]) + 1:tss + 1] = REGION_CODE["5'-UTR"]

    # 3. Promoter (highest priority)
    for i, tss, strand in tss_rows:
        if strand == '+':
            a = max(0, tss - PROMOTER_UP)
            b = tss          # exclusive: [a, tss)
        else:
            a = tss + 1
            b = min(CHROM_LEN + 1, tss + PROMOTER_UP + 1)
        arr[a:b] = REGION_CODE['Promoter']

    return arr


def assign_regions(census_df, region_arr):
    """Add 'region' column to census DataFrame using pre-built region_arr."""
    pos = census_df['position'].clip(0, CHROM_LEN).astype(int).values
    codes = region_arr[pos]
    return REGION_ORDER[codes] if False else pd.Categorical.from_codes(
        codes, categories=REGION_ORDER
    )


# ── Panel C: Genomic region distribution ─────────────────────────────────────

def panel_c(ax, c4, c6, region_arr):
    """Stacked 100% bar chart of genomic region distribution per motif group."""
    # Annotate all sites (all timepoints pooled, deduplicated by position+mod_type)
    c4_ann = c4[c4['chrom'] == CHROM].drop_duplicates(subset=['position', 'mod_type']).copy()
    c6_ann = c6[c6['chrom'] == CHROM].drop_duplicates(subset=['position', 'mod_type']).copy()
    c4_ann['region'] = assign_regions(c4_ann, region_arr)
    c6_ann['region'] = assign_regions(c6_ann, region_arr)

    census_by_type = {'4mC': c4_ann, '6mA': c6_ann}

    # Build count table per motif group
    rows = []
    for label, mod_type, motifs in MOTIF_GROUPS:
        df = census_by_type[mod_type]
        sub = df[df['final_motif'].isin(motifs)]
        counts = sub['region'].value_counts().reindex(REGION_ORDER, fill_value=0)
        total = counts.sum()
        rows.append({'label': label, 'total': total, **counts.to_dict()})
    df_plot = pd.DataFrame(rows).set_index('label')

    # Save region distribution table
    out_tsv = DIR_40 / 'region_distribution_by_motif.tsv'
    df_plot.reset_index().to_csv(out_tsv, sep='\t', index=False)
    print(f'  Saved region distribution: {out_tsv}')

    # 100% stacked bar
    x = np.arange(len(df_plot))
    bottoms = np.zeros(len(df_plot))
    for region in REGION_ORDER:
        fracs = df_plot[region].values / df_plot['total'].values * 100
        ax.bar(x, fracs, bottom=bottoms, width=0.65,
               color=REGION_COLORS[region], label=region,
               edgecolor='white', linewidth=0.4)
        # Annotate fraction if large enough
        for xi, (f, b) in enumerate(zip(fracs, bottoms)):
            if f >= 8:
                ax.text(xi, b + f / 2, f'{f:.0f}%', ha='center', va='center',
                        fontsize=6, color='white', fontweight='bold')
        bottoms += fracs

    # n= annotations above bars
    for xi, (lbl, row) in enumerate(df_plot.iterrows()):
        ax.text(xi, 102, f'n={int(row["total"]):,}', ha='center', va='bottom',
                fontsize=6, color='#444')

    ax.set_xticks(x)
    ax.set_xticklabels(df_plot.index, fontsize=7)
    ax.set_ylabel('% of sites', fontsize=8)
    ax.set_ylim(0, 115)
    ax.set_title('Genomic region distribution\nby methylation motif group',
                 fontsize=9, fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=6.5, frameon=False, loc='lower right',
              ncol=2, handlelength=1.0)
    _style(ax)


# ── Panel D: TSS protection zone by motif type ────────────────────────────────

def panel_d(ax, df_spatial):
    """TSS methylation gradient for regulatory genes, by site type."""
    site_configs = [
        ('GCCGGC_4mC',   'GCCGGC (4mC)',  COL_4mC,  '-',  2.0),
        ('All_4mC',      'All 4mC',        COL_4mC,  '--', 1.0),
        ('AAGCCCG_6mA',  'AAGCCCG (6mA)', COL_6mA,  '-',  2.0),
        ('All_6mA',      'All 6mA',        COL_6mA,  '--', 1.0),
    ]

    for site_type, label, color, ls, lw in site_configs:
        sub = df_spatial[
            (df_spatial['site_type'] == site_type) &
            (df_spatial['timepoint'] == 'T1') &
            (df_spatial['gene_category'] == 'regulatory')
        ].sort_values('bin_center_bp')
        if len(sub) == 0:
            continue
        x = sub['bin_center_bp'].values / 1000
        y = sub['density_per_kb_per_gene'].values
        ax.plot(x, y, color=color, linewidth=lw, linestyle=ls,
                alpha=0.9, label=label)

    # Non-regulatory reference (All methylation)
    ref = df_spatial[
        (df_spatial['site_type'] == 'All_methylation') &
        (df_spatial['timepoint'] == 'T1') &
        (df_spatial['gene_category'] == 'non_regulatory')
    ].sort_values('bin_center_bp')
    if len(ref) > 0:
        ax.plot(ref['bin_center_bp'].values / 1000,
                ref['density_per_kb_per_gene'].values,
                color='#BDBDBD', linewidth=1.0, linestyle=':',
                label='All methylation\n(non-regulatory)')

    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.axvspan(-0.3, 0.5, alpha=0.07, color='#1565C0', zorder=0,
               label='Protection zone')

    ax.set_xlabel('Distance from TSS (kb)', fontsize=8)
    ax.set_ylabel('Density (sites/kb/gene)', fontsize=8)
    ax.set_title('TSS protection zone by motif type\n(regulatory genes, T1)',
                 fontsize=9, fontweight='bold', linespacing=1.3)
    ax.legend(fontsize=6.5, frameon=True, fancybox=False, edgecolor='#ccc',
              loc='upper right', ncol=1)
    ax.text(-4.5, ax.get_ylim()[1] * 0.90 if ax.get_ylim()[1] > 0 else 0.3,
            'GCCGGC: depleted near TSS\nAAGCCCG: no depletion',
            fontsize=6.5, color='#444', style='italic', va='top')
    _style(ax)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    apply_style()
    print('=== FigS13: Comprehensive Per-motif Analysis ===')

    # Load data
    df_expr    = pd.read_csv(DIR_40 / 'expression_summary_stats.tsv', sep='\t')
    # df_func: Median LFC (T3 vs T1) by functional category × motif group.
    # No longer plotted in Panel C, but preserved as data for future use.
    df_func    = pd.read_csv(DIR_40 / 'functional_category_by_motif.tsv', sep='\t')
    c4         = pd.read_csv(CENSUS_4mC)
    c6         = pd.read_csv(CENSUS_6mA)
    df_spatial = pd.read_csv(DIR_48 / 'spatial_profile_data.tsv', sep='\t')
    gene_ann   = pd.read_csv(GENE_ANN_PATH, sep='\t')
    tss_df     = pd.read_csv(TSS_PATH)
    print(f'Expression summary: {len(df_expr)} groups')
    print(f'Functional categories (preserved): {df_func["category"].nunique()} categories × {df_func["motif_group"].nunique()} motif groups')
    print(f'4mC sites: {len(c4)}, 6mA sites: {len(c6)}')
    print(f'Genes: {len(gene_ann)}, TSS entries: {len(tss_df)}')

    print('Building region annotation array...')
    region_arr = build_region_array(gene_ann, tss_df)
    counts = {r: int((region_arr == REGION_CODE[r]).sum()) for r in REGION_ORDER}
    print(f'  Region bp counts: {counts}')

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(mm_to_inch(174), mm_to_inch(200)))
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           hspace=0.50, wspace=0.38,
                           left=0.09, right=0.97, top=0.96, bottom=0.07)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_a, 'a', x=-0.14, y=1.08)
    add_panel_label(ax_b, 'b', x=-0.14, y=1.08)
    add_panel_label(ax_c, 'c', x=-0.14, y=1.08)
    add_panel_label(ax_d, 'd', x=-0.14, y=1.08)

    print('Panel A: Expression correlation...')
    panel_a(ax_a, df_expr)

    print('Panel B: Timepoint distribution...')
    panel_b(ax_b, c4, c6)

    print('Panel C: Genomic region distribution...')
    panel_c(ax_c, c4, c6, region_arr)

    print('Panel D: TSS protection zone...')
    panel_d(ax_d, df_spatial)

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = FIG_SUP_DIR / 'FigS13_comprehensive_motif_analysis'
    save_figure(fig, out_path)
    print(f'\nSaved: {out_path}.pdf / .svg')
    print('=== Done ===')


if __name__ == '__main__':
    main()
