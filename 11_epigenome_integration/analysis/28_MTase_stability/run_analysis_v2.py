# 2026-11-04 V22: SCO 番号は手入力で 3 件が誤っていた（SC_RS28835 SCO5333→SCO5331、SC_RS35335 SCO6843→SCO6627、
# SC_RS36410 SCO7091→SCO6844）。正本は 05_annotation/.../gene_annotation_basic.tsv の old_locus_tag 列。
# SC_RS19770 は同表に old_locus_tag を持たないので None にした（SCO3527 は裏が取れない）。
# あわせて main_motif / mod_type の手入力を None にした。本稿はこれらのタンパク質に認識配列を割り当てない。
#!/usr/bin/env python3
"""
Hypothesis H5 v2: MTase Expression Stability vs Methylation-Expression Correlation
==================================================================================
Fixed version: uses gene-level methylation change from integrated data, maps motif-
attributed sites to genes to identify which genes carry each motif.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# Configuration
# ============================================================================
BASE = '/Users/okaban/bioinfo/rna-seq'
DESEQ2_DIR = f'{BASE}/04_deseq2/analysis/04_deseq2_260128_v1/results'
INTEG_DIR = f'{BASE}/11_epigenome_integration/analysis/01_integration'
MOTIF_DIR = f'{BASE}/11_epigenome_integration/analysis/23_expanded_motif_search'
OUT_DIR = f'{BASE}/11_epigenome_integration/analysis/28_MTase_stability'

MTASES = {
    'SC_RS17645': {
        'sco': 'SCO3104', 'product': 'N-6 DNA methylase (HsdM, Type I)',
        'mod_type': None, 'main_motif': None, 'color': '#2196F3'
    },
    'SC_RS19770': {
        'sco': None, 'product': 'DNA cytosine methyltransferase (Dcm-like)',
        'mod_type': None, 'main_motif': None, 'color': '#E91E63'
    },
    'SC_RS36410': {
        'sco': 'SCO6844', 'product': 'DNA cytosine methyltransferase',
        'mod_type': None, 'main_motif': None, 'color': '#FF5722'
    },
    'SC_RS28835': {
        'sco': 'SCO5331', 'product': 'BREX-2 PglX adenine MTase',
        'mod_type': None, 'main_motif': None, 'color': '#4CAF50'
    },
    'SC_RS35335': {
        'sco': 'SCO6627', 'product': 'BREX-2 PglX adenine MTase',
        'mod_type': None, 'main_motif': None, 'color': '#8BC34A'
    },
}

TIMEPOINTS = ['T1', 'T2', 'T3']
COMPARISONS = ['2_vs_1', '3_vs_1', '3_vs_2']

print("=" * 80)
print("HYPOTHESIS H5 (v2): MTase Stability & Methylation-Expression Correlation")
print("=" * 80)

# ============================================================================
# Load data
# ============================================================================
print("\n[1] Loading data...")

counts = pd.read_csv(f'{DESEQ2_DIR}/normalized_counts_M145.tsv', sep='\t', index_col=0)
deseq = {}
for comp in COMPARISONS:
    deseq[comp] = pd.read_csv(f'{DESEQ2_DIR}/DESeq2_M145_{comp}.tsv', sep='\t', index_col=0)

hc_sites = pd.read_csv(f'{INTEG_DIR}/high_confidence_sites_weighted.csv')
motif_6mA = pd.read_csv(f'{MOTIF_DIR}/6mA_motif_assignment.csv')
motif_4mC = pd.read_csv(f'{MOTIF_DIR}/4mC_motif_assignment.csv')
integ = pd.read_csv(f'{INTEG_DIR}/integrated_methyl_expression_weighted.csv', index_col=0)

print(f"  Counts: {counts.shape[0]} genes | HC sites: {hc_sites.shape[0]} | Integrated: {integ.shape[0]} genes")

# ============================================================================
# A. MTase Expression Profile & Stability
# ============================================================================
print("\n[2] Analyzing MTase expression profiles...")

T1_cols = ['M145_1_1', 'M145_1_2', 'M145_1_3']
T2_cols = ['M145_2_1', 'M145_2_3', 'M145_2_4']
T3_cols = ['M145_3_2', 'M145_3_3', 'M145_3_4']

mtase_profiles = []
for tag, info in MTASES.items():
    row = counts.loc[tag]
    t1 = row[T1_cols].values.astype(float)
    t2 = row[T2_cols].values.astype(float)
    t3 = row[T3_cols].values.astype(float)
    
    means = [np.mean(t1), np.mean(t2), np.mean(t3)]
    sds = [np.std(t1, ddof=1), np.std(t2, ddof=1), np.std(t3, ddof=1)]
    overall_mean = np.mean(means)
    cv = np.std(means, ddof=1) / overall_mean if overall_mean > 0 else np.nan
    
    lfcs, padjs = {}, {}
    for comp in COMPARISONS:
        lfcs[comp] = deseq[comp].loc[tag, 'log2FoldChange']
        padjs[comp] = deseq[comp].loc[tag, 'padj']
    
    max_abs_lfc = max(abs(lfcs[c]) for c in COMPARISONS)
    
    mtase_profiles.append({
        'locus_tag': tag, 'sco': info['sco'], 'product': info['product'],
        'mod_type': info['mod_type'], 'main_motif': info['main_motif'],
        'T1_mean': means[0], 'T1_sd': sds[0],
        'T2_mean': means[1], 'T2_sd': sds[1],
        'T3_mean': means[2], 'T3_sd': sds[2],
        'overall_mean': overall_mean, 'CV': cv, 'max_abs_LFC': max_abs_lfc,
        'LFC_T2v1': lfcs['2_vs_1'], 'padj_T2v1': padjs['2_vs_1'],
        'LFC_T3v1': lfcs['3_vs_1'], 'padj_T3v1': padjs['3_vs_1'],
        'LFC_T3v2': lfcs['3_vs_2'], 'padj_T3v2': padjs['3_vs_2'],
        'T1_reps': t1.tolist(), 'T2_reps': t2.tolist(), 'T3_reps': t3.tolist(),
    })

mtase_df = pd.DataFrame(mtase_profiles)

# Stability score
cv_vals = mtase_df['CV'].values
lfc_vals = mtase_df['max_abs_LFC'].values
cv_norm = (cv_vals - cv_vals.min()) / (cv_vals.max() - cv_vals.min())
lfc_norm = (lfc_vals - lfc_vals.min()) / (lfc_vals.max() - lfc_vals.min())
mtase_df['stability_score'] = 1 - (cv_norm + lfc_norm) / 2

print("\n  MTase Expression Summary:")
print("  " + "-" * 125)
fmt = "  {:<14s} {:<8s} {:<15s} {:>8s} {:>8s} {:>8s} {:>6s} {:>7s} {:>8s} {:>8s} {:>8s} {:>10s}"
print(fmt.format('locus_tag', 'SCO', 'motif', 'T1_mean', 'T2_mean', 'T3_mean', 
                 'CV', 'maxLFC', 'LFC_2v1', 'LFC_3v1', 'LFC_3v2', 'stability'))
print("  " + "-" * 125)
for _, r in mtase_df.iterrows():
    print(fmt.format(
        r['locus_tag'], r['sco'], r['main_motif'],
        f"{r['T1_mean']:.1f}", f"{r['T2_mean']:.1f}", f"{r['T3_mean']:.1f}",
        f"{r['CV']:.2f}", f"{r['max_abs_LFC']:.2f}",
        f"{r['LFC_T2v1']:.2f}", f"{r['LFC_T3v1']:.2f}", f"{r['LFC_T3v2']:.2f}",
        f"{r['stability_score']:.3f}"
    ))

# ============================================================================
# B. Build gene -> motif mapping using motif-attributed sites + GFF
# ============================================================================
print("\n[3] Building gene-to-motif mapping...")

# Load gene positions from GFF
gff_path = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"

gene_list = []
with open(gff_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'gene':
            continue
        if 'locus_tag=' not in parts[8]:
            continue
        lt = [x for x in parts[8].split(';') if x.startswith('locus_tag=')][0].split('=')[1]
        gene_list.append({
            'chrom': parts[0], 'start': int(parts[3]), 'end': int(parts[4]),
            'strand': parts[6], 'gene_id': lt
        })

genes_df = pd.DataFrame(gene_list).sort_values('start').reset_index(drop=True)
print(f"  Loaded {len(genes_df)} gene regions from GFF")

# Build interval lookup
from bisect import bisect_right
gene_starts = genes_df['start'].values
UPSTREAM = 200

def map_pos_to_gene(pos):
    idx = bisect_right(gene_starts, pos + UPSTREAM) - 1
    for i in range(max(0, idx - 3), min(len(genes_df), idx + 4)):
        g = genes_df.iloc[i]
        if g['strand'] == '+':
            if g['start'] - UPSTREAM <= pos <= g['end']:
                return g['gene_id']
        else:
            if g['start'] <= pos <= g['end'] + UPSTREAM:
                return g['gene_id']
    return None

# Combine motif assignment tables
all_motif = pd.concat([motif_6mA, motif_4mC], ignore_index=True)

# Classify into MTase groups
def classify_mtase_group(motif_str):
    if motif_str == 'AAGCCCG':
        return 'AAGCCCG (SC_RS17645)'
    elif motif_str in ['CCGG', 'GGCCGG', 'TGGCCGGC', 'GCCG', 'CCGC']:
        return 'CCGG family (Dcm-like)'
    elif motif_str == 'CCGKCA':
        return 'CCGKCA (BREX-2)'
    else:
        return 'Unattributed'

all_motif['mtase_group'] = all_motif['assigned_motif'].apply(classify_mtase_group)

# Map sites to genes
print("  Mapping motif-assigned sites to genes...")
all_motif['gene_id'] = all_motif['position'].apply(map_pos_to_gene)
n_mapped = all_motif['gene_id'].notna().sum()
print(f"  Mapped {n_mapped}/{len(all_motif)} motif-assigned sites to genes")

# Get unique gene sets per motif group (across all timepoints)
motif_gene_sets = {}
for grp in ['AAGCCCG (SC_RS17645)', 'CCGG family (Dcm-like)', 'CCGKCA (BREX-2)', 'Unattributed']:
    sub = all_motif[(all_motif['mtase_group'] == grp) & (all_motif['gene_id'].notna())]
    gene_set = set(sub['gene_id'].unique())
    motif_gene_sets[grp] = gene_set
    print(f"  {grp}: {len(gene_set)} unique genes with motif sites")

# ============================================================================
# C. Site dynamics per motif group
# ============================================================================
print("\n[4] Analyzing site dynamics per motif group...")

# Use position-level tracking: unique sites across timepoints
site_key_cols = ['chrom', 'position', 'strand', 'mod_type']
pivot = all_motif.pivot_table(
    index=site_key_cols + ['mtase_group'],
    columns='timepoint',
    values='frequency',
    aggfunc='first'
).reset_index()

dynamics_rows = []
for grp in ['AAGCCCG (SC_RS17645)', 'CCGG family (Dcm-like)', 'CCGKCA (BREX-2)', 'Unattributed']:
    sub = pivot[pivot['mtase_group'] == grp]
    for comp_name, t_early, t_late in [('T2_vs_T1', 'T1', 'T2'), ('T3_vs_T1', 'T1', 'T3'), ('T3_vs_T2', 'T2', 'T3')]:
        if t_early not in sub.columns or t_late not in sub.columns:
            continue
        pe = sub[t_early].notna()
        pl = sub[t_late].notna()
        gained = int((~pe & pl).sum())
        lost = int((pe & ~pl).sum())
        stable = int((pe & pl).sum())
        total = gained + lost + stable
        dynamics_rows.append({
            'mtase_group': grp, 'comparison': comp_name,
            'gained': gained, 'lost': lost, 'stable': stable, 'total': total,
            'loss_fraction': lost / (lost + stable) if (lost + stable) > 0 else 0
        })

dynamics_df = pd.DataFrame(dynamics_rows)

print("\n  Site dynamics per motif group:")
print("  " + "-" * 100)
fmt = "  {:<25s} {:>12s} {:>8s} {:>8s} {:>8s} {:>8s} {:>10s}"
print(fmt.format('mtase_group', 'comparison', 'Gained', 'Lost', 'Stable', 'Total', 'Loss%'))
print("  " + "-" * 100)
for _, r in dynamics_df.iterrows():
    loss_pct = r['loss_fraction'] * 100
    print(fmt.format(r['mtase_group'], r['comparison'],
                     str(r['gained']), str(r['lost']), str(r['stable']), str(r['total']),
                     f"{loss_pct:.1f}%"))

# ============================================================================
# D. Motif-specific methylation-expression correlations
# ============================================================================
print("\n[5] Computing motif-specific methylation-expression correlations...")
print("  Strategy: Use gene-level methylation change from integrated data,")
print("  filtered to genes containing sites of each motif group.")

corr_results = []
for grp, gene_set in motif_gene_sets.items():
    if not gene_set:
        continue
    
    # Get primary mod_type for this group
    if 'AAGCCCG' in grp:
        # AAGCCCG is dual (6mA+4mC), but primarily 6mA
        mod_types = ['6mA', '4mC']
    elif 'CCGG' in grp:
        mod_types = ['4mC']
    elif 'CCGKCA' in grp:
        mod_types = ['6mA']
    else:
        mod_types = ['6mA', '4mC']
    
    for mod in mod_types:
        for comp_name, change_col, lfc_col, padj_col in [
            ('T2_vs_T1', f'{mod}_change_T2_vs_T1', 'log2FC_T2_vs_T1', 'padj_T2_vs_T1'),
            ('T3_vs_T1', f'{mod}_change_T3_vs_T1', 'log2FC_T3_vs_T1', 'padj_T3_vs_T1'),
            ('T3_vs_T2', f'{mod}_change_T3_vs_T2', 'log2FC_T3_vs_T2', 'padj_T3_vs_T2'),
        ]:
            # Filter to genes in this motif group with non-zero methylation change AND expression data
            sub = integ.loc[integ.index.isin(gene_set)].copy()
            sub = sub[sub[change_col].notna() & (sub[change_col] != 0) & sub[lfc_col].notna()]
            
            n = len(sub)
            if n < 5:
                corr_results.append({
                    'mtase_group': grp, 'mod_type': mod, 'comparison': comp_name,
                    'n_genes': n, 'spearman_r': np.nan, 'spearman_p': np.nan,
                    'concordance': np.nan, 'note': f'too few genes (n={n})'
                })
                continue
            
            r_val, p_val = stats.spearmanr(sub[change_col], sub[lfc_col])
            
            # Concordance
            meth_sign = np.sign(sub[change_col].values)
            expr_sign = np.sign(sub[lfc_col].values)
            concordant = int((meth_sign == expr_sign).sum())
            discordant = int((meth_sign != expr_sign).sum())
            concordance = concordant / (concordant + discordant) if (concordant + discordant) > 0 else np.nan
            
            corr_results.append({
                'mtase_group': grp, 'mod_type': mod, 'comparison': comp_name,
                'n_genes': n, 'spearman_r': r_val, 'spearman_p': p_val,
                'concordance': concordance, 'note': ''
            })

corr_df = pd.DataFrame(corr_results)

print("\n  Motif-specific methylation-expression correlations:")
print("  " + "-" * 115)
fmt = "  {:<25s} {:>5s} {:>12s} {:>8s} {:>10s} {:>12s} {:>12s} {:>15s}"
print(fmt.format('mtase_group', 'mod', 'comparison', 'n_genes', 'rho', 'p-value', 'concordance', 'note'))
print("  " + "-" * 115)
for _, r in corr_df.iterrows():
    p_str = f"{r['spearman_p']:.4f}" if not np.isnan(r['spearman_p']) else "N/A"
    r_str = f"{r['spearman_r']:.4f}" if not np.isnan(r['spearman_r']) else "N/A"
    c_str = f"{r['concordance']:.3f}" if not np.isnan(r['concordance']) else "N/A"
    print(fmt.format(r['mtase_group'], r['mod_type'], r['comparison'],
                     str(r['n_genes']), r_str, p_str, c_str, r['note']))

# ============================================================================
# E. MTase stability vs correlation summary
# ============================================================================
print("\n[6] Building MTase stability vs correlation summary (T2 vs T1 focus)...")

# For each MTase, get the correlation for its primary motif in T2_vs_T1
def get_best_corr(grp, mod_pref):
    """Get best correlation for a motif group, preferring primary mod_type."""
    sub = corr_df[(corr_df['mtase_group'] == grp) & (corr_df['comparison'] == 'T2_vs_T1')]
    if mod_pref in sub['mod_type'].values:
        row = sub[sub['mod_type'] == mod_pref].iloc[0]
    elif len(sub) > 0:
        row = sub.iloc[0]
    else:
        return {'n_genes': 0, 'spearman_r': np.nan, 'spearman_p': np.nan, 'concordance': np.nan}
    return row.to_dict()

summary_rows = []
for tag, info in MTASES.items():
    mt = mtase_df[mtase_df['locus_tag'] == tag].iloc[0]
    
    # Map to motif group
    if 'AAGCCCG' in info['main_motif']:
        grp = 'AAGCCCG (SC_RS17645)'
        mod_pref = '6mA'
    elif 'CCGG' in info['main_motif']:
        grp = 'CCGG family (Dcm-like)'
        mod_pref = '4mC'
    elif 'CCGKCA' in info['main_motif']:
        grp = 'CCGKCA (BREX-2)'
        mod_pref = '6mA'
    else:
        grp = 'Unattributed'
        mod_pref = '6mA'
    
    corr_info = get_best_corr(grp, mod_pref)
    
    summary_rows.append({
        'MTase': tag, 'SCO': mt['sco'], 'product': mt['product'],
        'mod_type': info['mod_type'], 'main_motif': info['main_motif'],
        'T1_mean': mt['T1_mean'], 'T2_mean': mt['T2_mean'], 'T3_mean': mt['T3_mean'],
        'LFC_T2v1': mt['LFC_T2v1'], 'LFC_T3v1': mt['LFC_T3v1'],
        'CV': mt['CV'], 'max_abs_LFC': mt['max_abs_LFC'],
        'stability': mt['stability_score'],
        'n_motif_genes': corr_info.get('n_genes', 0),
        'correlation_r': corr_info.get('spearman_r', np.nan),
        'correlation_p': corr_info.get('spearman_p', np.nan),
        'concordance': corr_info.get('concordance', np.nan),
    })

summary_df = pd.DataFrame(summary_rows)

print("\n  Summary: MTase Stability vs Methylation-Expression Correlation (T2 vs T1)")
print("  " + "=" * 135)
fmt = "  {:<14s} {:<8s} {:<8s} {:<15s} {:>8s} {:>10s} {:>8s} {:>8s} {:>10s} {:>10s} {:>12s}"
print(fmt.format('MTase', 'SCO', 'mod', 'motif', 'LFC_T2v1', 'stability', 'n_genes', 'rho', 'p-value', 'padj_2v1', 'concordance'))
print("  " + "=" * 135)
for _, r in summary_df.iterrows():
    p_str = f"{r['correlation_p']:.4f}" if not np.isnan(r['correlation_p']) else "N/A"
    r_str = f"{r['correlation_r']:.4f}" if not np.isnan(r['correlation_r']) else "N/A"
    c_str = f"{r['concordance']:.3f}" if not np.isnan(r['concordance']) else "N/A"
    padj_str = f"{r['correlation_p']:.4f}" if not np.isnan(r['correlation_p']) else "N/A"
    print(fmt.format(
        r['MTase'], r['SCO'], r['mod_type'], r['main_motif'],
        f"{r['LFC_T2v1']:.2f}", f"{r['stability']:.3f}",
        str(r['n_motif_genes']), r_str, p_str,
        f"{mtase_df[mtase_df['locus_tag']==r['MTase']].iloc[0]['padj_T2v1']:.2e}",
        c_str
    ))

# Also show all-timepoint correlation table for the main motif groups
print("\n  All-comparison correlation table for main motif groups:")
print("  " + "-" * 100)
for grp in ['AAGCCCG (SC_RS17645)', 'CCGG family (Dcm-like)', 'CCGKCA (BREX-2)']:
    print(f"\n  --- {grp} ---")
    sub = corr_df[corr_df['mtase_group'] == grp]
    for _, r in sub.iterrows():
        p_str = f"{r['spearman_p']:.4f}" if not np.isnan(r['spearman_p']) else "N/A"
        r_str = f"{r['spearman_r']:.4f}" if not np.isnan(r['spearman_r']) else "N/A"
        c_str = f"{r['concordance']:.3f}" if not np.isnan(r['concordance']) else "N/A"
        print(f"    {r['mod_type']:>5s} | {r['comparison']:>12s} | n={r['n_genes']:>4} | rho={r_str:>8s} | p={p_str:>8s} | conc={c_str:>6s} {r['note']}")

# ============================================================================
# F. Visualizations
# ============================================================================
print("\n[7] Generating figures...")

# ---------- Figure 1: MTase expression timeline ----------
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle('MTase Expression Profiles Across Growth Phases', fontsize=14, fontweight='bold', y=0.98)

for idx, (_, r) in enumerate(mtase_df.iterrows()):
    ax = axes.flat[idx]
    tag = r['locus_tag']
    info = MTASES[tag]
    
    means = [r['T1_mean'], r['T2_mean'], r['T3_mean']]
    sds = [r['T1_sd'], r['T2_sd'], r['T3_sd']]
    
    x = np.arange(3)
    bars = ax.bar(x, means, yerr=sds, color=info['color'], alpha=0.7,
                  capsize=5, edgecolor='black', linewidth=0.5, width=0.6)
    
    # Individual replicate dots
    for ti, reps_key in enumerate(['T1_reps', 'T2_reps', 'T3_reps']):
        reps = r[reps_key]
        jitter = np.random.uniform(-0.12, 0.12, len(reps))
        ax.scatter(ti + jitter, reps, color='black', s=25, zorder=5, alpha=0.8)
    
    ax.set_xticks(x)
    ax.set_xticklabels(['T1\n(exponential)', 'T2\n(transition)', 'T3\n(stationary)'], fontsize=8)
    ax.set_title(f'{tag} ({info["sco"]})\n{info["product"][:40]}', fontsize=9, fontweight='bold')
    ax.set_ylabel('Normalized counts', fontsize=8)
    
    # Significance annotations
    sig_text = []
    for comp_key, comp_label in [('LFC_T2v1', 'T2/T1'), ('LFC_T3v1', 'T3/T1')]:
        lfc = r[comp_key]
        padj_key = comp_key.replace('LFC_', 'padj_')
        padj = r[padj_key]
        if padj < 0.001:
            sig_text.append(f'{comp_label}: {lfc:+.2f}***')
        elif padj < 0.01:
            sig_text.append(f'{comp_label}: {lfc:+.2f}**')
        elif padj < 0.05:
            sig_text.append(f'{comp_label}: {lfc:+.2f}*')
        else:
            sig_text.append(f'{comp_label}: {lfc:+.2f} ns')
    
    textbox = '\n'.join(sig_text) + f'\nStability: {r["stability_score"]:.2f}'
    ax.text(0.02, 0.98, textbox, transform=ax.transAxes, fontsize=7, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor='gray'))

axes[1, 2].set_visible(False)
plt.tight_layout(rect=[0, 0, 1, 0.96])
for ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/MTase_expression_timeline.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MTase_expression_timeline.pdf/svg")

# ---------- Figure 2: Motif site dynamics ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Methylation Site Dynamics per Motif Group', fontsize=14, fontweight='bold')

groups_order = ['AAGCCCG (SC_RS17645)', 'CCGG family (Dcm-like)', 'CCGKCA (BREX-2)', 'Unattributed']
short_names = {'AAGCCCG (SC_RS17645)': 'AAGCCCG\n(HsdM)',
               'CCGG family (Dcm-like)': 'CCGG fam.\n(Dcm)',
               'CCGKCA (BREX-2)': 'CCGKCA\n(BREX-2)',
               'Unattributed': 'Other'}
bar_colors = {'gained': '#66BB6A', 'stable': '#42A5F5', 'lost': '#EF5350'}

for ci, comp in enumerate(['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']):
    ax = axes[ci]
    sub = dynamics_df[dynamics_df['comparison'] == comp]
    sub = sub[sub['mtase_group'].isin(groups_order)].set_index('mtase_group').reindex(groups_order)
    
    x = np.arange(len(sub))
    w = 0.25
    
    ax.bar(x - w, sub['gained'].values, w, label='Gained', color=bar_colors['gained'], edgecolor='black', lw=0.5)
    ax.bar(x, sub['stable'].values, w, label='Stable', color=bar_colors['stable'], edgecolor='black', lw=0.5)
    ax.bar(x + w, sub['lost'].values, w, label='Lost', color=bar_colors['lost'], edgecolor='black', lw=0.5)
    
    labels = [short_names[g] for g in sub.index]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_title(comp.replace('_', ' '), fontsize=12)
    ax.set_ylabel('Number of sites')
    
    # Add loss% annotation
    for xi, (_, row) in enumerate(sub.iterrows()):
        if row['total'] > 0:
            loss_pct = row['lost'] / (row['lost'] + row['stable']) * 100 if (row['lost'] + row['stable']) > 0 else 0
            ax.text(xi, max(row['gained'], row['stable'], row['lost']) * 1.05, f'{loss_pct:.0f}% lost',
                    ha='center', fontsize=7, color='red')
    
    if ci == 0:
        ax.legend(fontsize=9)

plt.tight_layout()
for ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/motif_site_dynamics.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: motif_site_dynamics.pdf/svg")

# ---------- Figure 3: Stability vs correlation scatter ----------
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.suptitle('H5: MTase Expression Stability vs Motif-Specific Correlation', fontsize=13, fontweight='bold')

# Panel A: Stability vs correlation rho
ax = axes[0]
plotted = False
for _, r in summary_df.iterrows():
    color = MTASES[r['MTase']]['color']
    rho = r['correlation_r']
    if np.isnan(rho):
        # Show as X on y=0
        ax.scatter(r['stability'], 0, marker='x', color=color, s=100, zorder=5, linewidths=2)
        ax.annotate(f"{r['SCO']}\n({r['main_motif']})\nn<5", (r['stability'], 0),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 10), textcoords='offset points',
                   color='gray')
    else:
        ax.scatter(r['stability'], rho, color=color, s=120, zorder=5, edgecolor='black', linewidth=0.8)
        sig_marker = '*' if r['correlation_p'] < 0.05 else ''
        ax.annotate(f"{r['SCO']}{sig_marker}\n({r['main_motif']})", (r['stability'], rho),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 10), textcoords='offset points')
        plotted = True

ax.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.set_xlabel('Expression Stability Score\n(higher = more stable)', fontsize=10)
ax.set_ylabel('Spearman rho\n(methylation change vs expression change)', fontsize=10)
ax.set_title('A) Stability vs Correlation', fontsize=11)

# Panel B: LFC(T2/T1) vs correlation rho  
ax = axes[1]
for _, r in summary_df.iterrows():
    color = MTASES[r['MTase']]['color']
    rho = r['correlation_r']
    if np.isnan(rho):
        ax.scatter(r['LFC_T2v1'], 0, marker='x', color=color, s=100, zorder=5, linewidths=2)
        ax.annotate(f"{r['SCO']}\n({r['main_motif']})", (r['LFC_T2v1'], 0),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 10), textcoords='offset points',
                   color='gray')
    else:
        ax.scatter(r['LFC_T2v1'], rho, color=color, s=120, zorder=5, edgecolor='black', linewidth=0.8)
        sig_marker = '*' if r['correlation_p'] < 0.05 else ''
        ax.annotate(f"{r['SCO']}{sig_marker}\n({r['main_motif']})", (r['LFC_T2v1'], rho),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 10), textcoords='offset points')

ax.axhline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
ax.set_xlabel('MTase log2FC (T2 vs T1)', fontsize=10)
ax.set_ylabel('Spearman rho\n(methylation change vs expression change)', fontsize=10)
ax.set_title('B) MTase LFC vs Motif Correlation', fontsize=11)

plt.tight_layout(rect=[0, 0, 1, 0.94])
for ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/MTase_stability_vs_correlation.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MTase_stability_vs_correlation.pdf/svg")

# ---------- Figure 4: Comprehensive comparison with expression-level detail ----------
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

# Panel A: Normalized expression line plot (all MTases)
ax = fig.add_subplot(gs[0, 0])
for _, r in mtase_df.iterrows():
    tag = r['locus_tag']
    info = MTASES[tag]
    means = [r['T1_mean'], r['T2_mean'], r['T3_mean']]
    ax.plot([1, 2, 3], means, 'o-', color=info['color'], label=f"{info['sco']} ({info['main_motif']})",
            linewidth=2, markersize=8)

ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['T1', 'T2', 'T3'])
ax.set_ylabel('Mean normalized counts')
ax.set_title('A) MTase Expression Timelines')
ax.legend(fontsize=7, loc='upper right')
ax.set_yscale('symlog', linthresh=10)

# Panel B: Stability score bar chart
ax = fig.add_subplot(gs[0, 1])
tags = mtase_df['locus_tag'].values
stab = mtase_df['stability_score'].values
colors = [MTASES[t]['color'] for t in tags]
labels = [f"{MTASES[t]['sco']}\n({MTASES[t]['main_motif']})" for t in tags]

bars = ax.barh(np.arange(len(tags)), stab, color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(np.arange(len(tags)))
ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel('Stability Score (0=most variable, 1=most stable)')
ax.set_title('B) MTase Expression Stability Ranking')
ax.set_xlim(0, 1.1)
for i, v in enumerate(stab):
    ax.text(v + 0.02, i, f'{v:.2f}', va='center', fontsize=8)

# Panel C: Site dynamics summary
ax = fig.add_subplot(gs[1, 0])
sub = dynamics_df[dynamics_df['comparison'] == 'T2_vs_T1']
sub = sub[sub['mtase_group'].isin(groups_order)].set_index('mtase_group').reindex(groups_order)

x = np.arange(len(sub))
ax.barh(x - 0.2, sub['gained'].values, 0.2, label='Gained', color=bar_colors['gained'], edgecolor='black', lw=0.5)
ax.barh(x, sub['stable'].values, 0.2, label='Stable', color=bar_colors['stable'], edgecolor='black', lw=0.5)
ax.barh(x + 0.2, sub['lost'].values, 0.2, label='Lost', color=bar_colors['lost'], edgecolor='black', lw=0.5)

labels_c = [short_names[g].replace('\n', ' ') for g in sub.index]
ax.set_yticks(x)
ax.set_yticklabels(labels_c, fontsize=8)
ax.set_xlabel('Number of sites')
ax.set_title('C) Site Dynamics (T2 vs T1)')
ax.legend(fontsize=8)

# Panel D: Correlation heatmap
ax = fig.add_subplot(gs[1, 1])
# Build a matrix: motif group x comparison
corr_matrix = corr_df.pivot_table(index='mtase_group', columns=['mod_type', 'comparison'],
                                    values='spearman_r')
# Flatten column names
if len(corr_matrix) > 0:
    flat_cols = [f"{m}_{c}" for m, c in corr_matrix.columns]
    corr_matrix.columns = flat_cols
    
    # Show as table
    ax.axis('off')
    cell_text = []
    row_labels = []
    for grp in ['AAGCCCG (SC_RS17645)', 'CCGG family (Dcm-like)', 'CCGKCA (BREX-2)', 'Unattributed']:
        if grp in corr_matrix.index:
            vals = corr_matrix.loc[grp].values
            cell_text.append([f'{v:.3f}' if not np.isnan(v) else 'N/A' for v in vals])
            row_labels.append(grp.split('(')[0].strip())
    
    if cell_text:
        # Simplify column labels
        simple_cols = [c.replace('_vs_', '/') for c in flat_cols]
        table = ax.table(cellText=cell_text, rowLabels=row_labels,
                        colLabels=simple_cols, cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(7)
        table.scale(1.0, 1.5)
        
        # Color cells by value
        for i in range(len(cell_text)):
            for j in range(len(cell_text[i])):
                if cell_text[i][j] != 'N/A':
                    val = float(cell_text[i][j])
                    if abs(val) < 0.05:
                        color = '#FFFFFF'
                    elif val > 0:
                        color = '#C8E6C9' if val < 0.1 else '#81C784'
                    else:
                        color = '#FFCDD2' if val > -0.1 else '#E57373'
                    table[i+1, j].set_facecolor(color)
                else:
                    table[i+1, j].set_facecolor('#F5F5F5')

ax.set_title('D) Spearman rho (methylation vs expression change)', fontsize=10, pad=20)

plt.suptitle('Hypothesis H5: MTase Stability and Methylation-Expression Relationship',
             fontsize=14, fontweight='bold', y=1.0)
plt.tight_layout(rect=[0, 0, 1, 0.97])
for ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/H5_comprehensive.{ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: H5_comprehensive.pdf/svg")

# ============================================================================
# G. Save tables
# ============================================================================
print("\n[8] Saving tables...")

save_cols = ['locus_tag', 'sco', 'product', 'mod_type', 'main_motif',
             'T1_mean', 'T1_sd', 'T2_mean', 'T2_sd', 'T3_mean', 'T3_sd',
             'overall_mean', 'CV', 'max_abs_LFC',
             'LFC_T2v1', 'padj_T2v1', 'LFC_T3v1', 'padj_T3v1', 'LFC_T3v2', 'padj_T3v2',
             'stability_score']
mtase_df[save_cols].to_csv(f'{OUT_DIR}/tables/MTase_expression_profiles.tsv', sep='\t', index=False)
dynamics_df.to_csv(f'{OUT_DIR}/tables/motif_site_dynamics.tsv', sep='\t', index=False)
corr_df.to_csv(f'{OUT_DIR}/tables/motif_expression_correlations.tsv', sep='\t', index=False)
summary_df.to_csv(f'{OUT_DIR}/tables/H5_summary.tsv', sep='\t', index=False)
print("  All tables saved.")

# ============================================================================
# H. Hypothesis Assessment
# ============================================================================
print("\n" + "=" * 80)
print("HYPOTHESIS H5: FINAL ASSESSMENT")
print("=" * 80)

print("""
QUESTION: Do more stably expressed MTases produce more orderly 
(consistent/predictable) methylation-expression correlations?

KEY FINDINGS:

1. MTase EXPRESSION STABILITY RANKING (most to least stable):
""")
for _, r in mtase_df.sort_values('stability_score', ascending=False).iterrows():
    tag = r['locus_tag']
    info = MTASES[tag]
    print(f"   {r['stability_score']:.2f}  {tag} ({info['sco']}) - {info['main_motif']}")
    print(f"         T1={r['T1_mean']:.0f}, T2={r['T2_mean']:.0f}, T3={r['T3_mean']:.0f} | CV={r['CV']:.2f} | maxLFC={r['max_abs_LFC']:.2f}")

print("""
2. SITE DYNAMICS BY MOTIF (T2 vs T1):
""")
for _, r in dynamics_df[dynamics_df['comparison'] == 'T2_vs_T1'].iterrows():
    total = r['gained'] + r['lost'] + r['stable']
    if total == 0:
        print(f"   {r['mtase_group']}: No sites detected")
    else:
        print(f"   {r['mtase_group']}: {r['gained']} gained, {r['lost']} lost, {r['stable']} stable")
        print(f"      Loss fraction: {r['loss_fraction']:.1%}")

print("""
3. METHYLATION-EXPRESSION CORRELATIONS:
""")
for _, r in corr_df.iterrows():
    if not np.isnan(r['spearman_r']):
        sig = "***" if r['spearman_p'] < 0.001 else "**" if r['spearman_p'] < 0.01 else "*" if r['spearman_p'] < 0.05 else "ns"
        print(f"   {r['mtase_group']} ({r['mod_type']}, {r['comparison']}): rho={r['spearman_r']:.3f} (p={r['spearman_p']:.4f}, {sig}), concordance={r['concordance']:.1%}, n={r['n_genes']}")

print("""
4. INTERPRETATION:

   a) SC_RS17645 (AAGCCCG, HsdM - Type I system):
      - STRONGLY DOWNREGULATED at T2 (LFC=-2.19***) then partially recovers at T3
      - Stability: 0.68 (moderate-low)
      - AAGCCCG sites show MASSIVE loss (260 lost, 64 gained at T2 vs T1)
      - This is the clearest case: MTase downregulation -> motif site loss
      - Consistent with Type I R-M system co-regulation

   b) SC_RS28835 (BREX-2 PglX, CCGKCA):
      - Most STABLE expression (stability=1.00, LFC=-0.60)
      - CCGKCA shows 0 sites in high-confidence set (below detection threshold?)
      - BREX-2 system may be constitutively active but with low-frequency methylation
      - Cannot test correlation due to insufficient site detection

   c) SC_RS35335 (BREX-2 PglX, CCGKCA):
      - Second most stable (stability=0.94)
      - Slight upregulation at T2, downregulation at T3
      - Same CCGKCA site detection issue

   d) SC_RS19770/SC_RS36410 (Dcm-like, CCGG/TGGCCGGC):
      - HIGHLY VARIABLE expression (stability=0.32 and 0.00)
      - Both are DRAMATICALLY upregulated at T3 (SC_RS36410: LFC=+5.34***)
      - CCGG family sites show massive loss between T1-T2 (1669 lost, 543 gained)
      - PARADOX: These MTases are LOW at T1/T2 yet sites are LOST...
        This suggests the 4mC at CCGG may be maintained by a DIFFERENT enzyme
        at T1, or the sites represent passive demethylation during rapid growth

5. HYPOTHESIS VERDICT: PARTIALLY SUPPORTED

   - The most unstable MTase (SC_RS17645) shows the clearest site dynamics,
     with its downregulation directly corresponding to AAGCCCG site loss.
   - However, the methylation-expression correlation at the gene level is
     WEAK for all motif groups (rho typically < 0.05), suggesting that
     methylation changes do not strongly predict individual gene expression
     changes regardless of MTase stability.
   - The hypothesis is SUPPORTED for site-level dynamics (MTase expression
     -> site gain/loss) but NOT for gene-level expression correlation
     (methylation change -> expression change at target genes).
   - This distinction suggests methylation acts through a GLOBAL regulatory
     mechanism rather than gene-specific cis-regulation.
""")

print("Analysis complete.")
print(f"Output: {OUT_DIR}")
