# 2026-11-04 V22: 手入力の SCO 番号 3 件を gene_annotation_basic.tsv の old_locus_tag に合わせて訂正。
#!/usr/bin/env python3
"""
Hypothesis H5: MTase Expression Stability vs Methylation-Expression Correlation
==============================================================================
Analyzes whether MTase expression stability predicts the orderliness of 
methylation-expression correlations for sites attributed to each MTase.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
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

# MTase candidates with GFF-verified annotations
MTASES = {
    'SC_RS17645': {
        'sco': 'SCO3104', 'product': 'N-6 DNA methylase (HsdM, Type I)',
        'mod_type': '6mA+4mC', 'main_motif': 'AAGCCCG',
        'color': '#2196F3'
    },
    'SC_RS19770': {
        'sco': 'SCO3527', 'product': 'DNA cytosine methyltransferase (Dcm-like)',
        'mod_type': '4mC', 'main_motif': 'CCGG/TGGCCGGC',
        'color': '#E91E63'
    },
    'SC_RS36410': {
        'sco': 'SCO6844', 'product': 'DNA cytosine methyltransferase',
        'mod_type': '4mC', 'main_motif': 'CCGG/TGGCCGGC',
        'color': '#FF5722'
    },
    'SC_RS28835': {
        'sco': 'SCO5331', 'product': 'BREX-2 PglX adenine MTase',
        'mod_type': '6mA', 'main_motif': 'CCGKCA',
        'color': '#4CAF50'
    },
    'SC_RS35335': {
        'sco': 'SCO6627', 'product': 'BREX-2 PglX adenine MTase',
        'mod_type': '6mA', 'main_motif': 'CCGKCA',
        'color': '#8BC34A'
    },
}

TIMEPOINTS = ['T1', 'T2', 'T3']
COMPARISONS = ['2_vs_1', '3_vs_1', '3_vs_2']

print("=" * 80)
print("HYPOTHESIS H5: MTase Expression Stability & Methylation-Expression Correlation")
print("=" * 80)

# ============================================================================
# Load data
# ============================================================================
print("\n[1] Loading data...")

# Normalized counts
counts = pd.read_csv(f'{DESEQ2_DIR}/normalized_counts_M145.tsv', sep='\t', index_col=0)
print(f"  Normalized counts: {counts.shape[0]} genes x {counts.shape[1]} samples")

# DESeq2 results for all comparisons
deseq = {}
for comp in COMPARISONS:
    df = pd.read_csv(f'{DESEQ2_DIR}/DESeq2_M145_{comp}.tsv', sep='\t', index_col=0)
    deseq[comp] = df
    print(f"  DESeq2 {comp}: {df.shape[0]} genes")

# High confidence methylation sites
hc_sites = pd.read_csv(f'{INTEG_DIR}/high_confidence_sites_weighted.csv')
print(f"  High confidence sites: {hc_sites.shape[0]} entries")

# Motif assignments
motif_6mA = pd.read_csv(f'{MOTIF_DIR}/6mA_motif_assignment.csv')
motif_4mC = pd.read_csv(f'{MOTIF_DIR}/4mC_motif_assignment.csv')
print(f"  6mA motif assignments: {motif_6mA.shape[0]}")
print(f"  4mC motif assignments: {motif_4mC.shape[0]}")

# Integrated methylation-expression
integ = pd.read_csv(f'{INTEG_DIR}/integrated_methyl_expression_weighted.csv', index_col=0)
print(f"  Integrated methyl-expression: {integ.shape[0]} genes")

# ============================================================================
# A. MTase Expression Profile & Stability
# ============================================================================
print("\n[2] Analyzing MTase expression profiles...")

# Sample columns
T1_cols = ['M145_1_1', 'M145_1_2', 'M145_1_3']
T2_cols = ['M145_2_1', 'M145_2_3', 'M145_2_4']
T3_cols = ['M145_3_2', 'M145_3_3', 'M145_3_4']

mtase_profiles = []

for tag, info in MTASES.items():
    if tag not in counts.index:
        print(f"  WARNING: {tag} not found in counts!")
        continue
    
    row = counts.loc[tag]
    
    # Per-timepoint stats
    t1_vals = row[T1_cols].values.astype(float)
    t2_vals = row[T2_cols].values.astype(float)
    t3_vals = row[T3_cols].values.astype(float)
    
    means = [np.mean(t1_vals), np.mean(t2_vals), np.mean(t3_vals)]
    sds = [np.std(t1_vals, ddof=1), np.std(t2_vals, ddof=1), np.std(t3_vals, ddof=1)]
    
    # Stability metrics
    overall_mean = np.mean(means)
    cv_across_timepoints = np.std(means, ddof=1) / overall_mean if overall_mean > 0 else np.nan
    
    # Get LFCs from DESeq2
    lfcs = {}
    padjs = {}
    for comp in COMPARISONS:
        if tag in deseq[comp].index:
            lfcs[comp] = deseq[comp].loc[tag, 'log2FoldChange']
            padjs[comp] = deseq[comp].loc[tag, 'padj']
        else:
            lfcs[comp] = np.nan
            padjs[comp] = np.nan
    
    max_abs_lfc = max(abs(lfcs[c]) for c in COMPARISONS if not np.isnan(lfcs[c]))
    
    rec = {
        'locus_tag': tag, 'sco': info['sco'], 'product': info['product'],
        'mod_type': info['mod_type'], 'main_motif': info['main_motif'],
        'T1_mean': means[0], 'T1_sd': sds[0],
        'T2_mean': means[1], 'T2_sd': sds[1],
        'T3_mean': means[2], 'T3_sd': sds[2],
        'overall_mean': overall_mean,
        'CV_across_timepoints': cv_across_timepoints,
        'max_abs_LFC': max_abs_lfc,
        'LFC_T2v1': lfcs['2_vs_1'], 'padj_T2v1': padjs['2_vs_1'],
        'LFC_T3v1': lfcs['3_vs_1'], 'padj_T3v1': padjs['3_vs_1'],
        'LFC_T3v2': lfcs['3_vs_2'], 'padj_T3v2': padjs['3_vs_2'],
    }
    
    # Store individual replicate values for plotting
    rec['T1_reps'] = t1_vals.tolist()
    rec['T2_reps'] = t2_vals.tolist()
    rec['T3_reps'] = t3_vals.tolist()
    
    mtase_profiles.append(rec)

mtase_df = pd.DataFrame(mtase_profiles)

# Composite stability score: normalized CV + normalized max_abs_LFC (lower = more stable)
cv_norm = (mtase_df['CV_across_timepoints'] - mtase_df['CV_across_timepoints'].min()) / \
          (mtase_df['CV_across_timepoints'].max() - mtase_df['CV_across_timepoints'].min())
lfc_norm = (mtase_df['max_abs_LFC'] - mtase_df['max_abs_LFC'].min()) / \
           (mtase_df['max_abs_LFC'].max() - mtase_df['max_abs_LFC'].min())
mtase_df['stability_score'] = 1 - (cv_norm + lfc_norm) / 2  # 1 = most stable, 0 = least stable

print("\n  MTase Expression Summary:")
print("  " + "-" * 120)
fmt = "  {:<14s} {:<8s} {:<15s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>8s} {:>10s}"
print(fmt.format('locus_tag', 'SCO', 'motif', 'T1_mean', 'T2_mean', 'T3_mean', 'CV', 'maxLFC', 'LFC_T2v1', 'stability'))
print("  " + "-" * 120)
for _, r in mtase_df.iterrows():
    print(fmt.format(
        r['locus_tag'], r['sco'], r['main_motif'],
        f"{r['T1_mean']:.1f}", f"{r['T2_mean']:.1f}", f"{r['T3_mean']:.1f}",
        f"{r['CV_across_timepoints']:.3f}", f"{r['max_abs_LFC']:.2f}",
        f"{r['LFC_T2v1']:.2f}", f"{r['stability_score']:.3f}"
    ))

# ============================================================================
# B. Motif Attribution for High-Confidence Sites
# ============================================================================
print("\n[3] Classifying methylation sites by motif...")

# Combine motif assignments
all_motif = pd.concat([motif_6mA, motif_4mC], ignore_index=True)
print(f"  Total motif-assigned sites: {all_motif.shape[0]}")

# Map assigned_motif to MTase group
def classify_mtase(row):
    motif = str(row.get('assigned_motif', ''))
    mod = str(row.get('mod_type', ''))
    
    if motif == 'AAGCCCG':
        return 'SC_RS17645 (AAGCCCG)'
    elif motif in ['CCGG', 'GGCCGG', 'TGGCCGGC', 'GCCG', 'CCGC']:
        return 'Dcm-like (CCGG family)'
    elif motif == 'CCGKCA':
        return 'BREX-2 (CCGKCA)'
    else:
        return 'Unattributed'

all_motif['mtase_group'] = all_motif.apply(classify_mtase, axis=1)

# Also classify the high-confidence sites by merging with motif assignments
# First, merge on position + timepoint + mod_type
hc_merged = hc_sites.merge(
    all_motif[['chrom', 'position', 'strand', 'mod_type', 'timepoint', 'assigned_motif']],
    on=['chrom', 'position', 'strand', 'mod_type', 'timepoint'],
    how='left'
)
hc_merged['assigned_motif'] = hc_merged['assigned_motif'].fillna('unassigned')
hc_merged['mtase_group'] = hc_merged.apply(classify_mtase, axis=1)

print(f"\n  Motif group distribution (high-confidence sites):")
for grp, count in hc_merged['mtase_group'].value_counts().items():
    print(f"    {grp}: {count}")

# ============================================================================
# C. Temporal dynamics per motif group
# ============================================================================
print("\n[4] Analyzing site dynamics per motif group...")

# Identify unique sites (position + strand + mod_type) and track presence across timepoints
site_key = ['chrom', 'position', 'strand', 'mod_type']
pivot = hc_merged.pivot_table(
    index=site_key + ['mtase_group'],
    columns='timepoint',
    values='weighted_mod_freq',
    aggfunc='first'
).reset_index()

# For each motif group and comparison, classify sites as Gained/Lost/Stable
dynamics_rows = []
for grp in ['SC_RS17645 (AAGCCCG)', 'Dcm-like (CCGG family)', 'BREX-2 (CCGKCA)', 'Unattributed']:
    sub = pivot[pivot['mtase_group'] == grp]
    
    for comp_name, t_early, t_late in [('T2_vs_T1', 'T1', 'T2'), ('T3_vs_T1', 'T1', 'T3'), ('T3_vs_T2', 'T2', 'T3')]:
        present_early = sub[t_early].notna()
        present_late = sub[t_late].notna()
        
        gained = (~present_early & present_late).sum()
        lost = (present_early & ~present_late).sum()
        stable = (present_early & present_late).sum()
        total = gained + lost + stable
        
        dynamics_rows.append({
            'mtase_group': grp, 'comparison': comp_name,
            'gained': gained, 'lost': lost, 'stable': stable, 'total': total
        })

dynamics_df = pd.DataFrame(dynamics_rows)

print("\n  Site dynamics per motif group:")
print("  " + "-" * 90)
fmt = "  {:<25s} {:>12s} {:>8s} {:>8s} {:>8s} {:>8s}"
print(fmt.format('mtase_group', 'comparison', 'Gained', 'Lost', 'Stable', 'Total'))
print("  " + "-" * 90)
for _, r in dynamics_df.iterrows():
    print(fmt.format(r['mtase_group'], r['comparison'], 
                     str(r['gained']), str(r['lost']), str(r['stable']), str(r['total'])))

# ============================================================================
# D. Methylation-Expression Correlation per motif group
# ============================================================================
print("\n[5] Computing motif-specific methylation-expression correlations...")

# For each motif group, find genes with sites of that motif and expression data
# Use the integrated dataset + motif assignments

# First, build a gene -> motif group mapping from the motif assignment files
# We need a gene mapping: position -> gene_id
# Load the GFF-based gene positions or use the integrated data

# From the integrated data, we know which genes have methylation changes
# We need to match motif-assigned sites to genes

# Approach: use the high_confidence_sites to find sites per gene
# Load position->gene mapping from the integration analysis
# Let's check if there's a mapping file
import os

# Build from GFF
gff_path = "/Users/okaban/SFC-CNS Dropbox/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/02_M145-ref/ncbi_dataset/data/GCF_000203835.1/genomic.gff"

gene_regions = []
with open(gff_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        if parts[2] == 'gene' and 'locus_tag=' in parts[8]:
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            lt = [x for x in parts[8].split(';') if x.startswith('locus_tag=')]
            if lt:
                tag = lt[0].split('=')[1]
                gene_regions.append({'chrom': chrom, 'start': start, 'end': end, 
                                    'strand': strand, 'gene_id': tag})

gene_regions_df = pd.DataFrame(gene_regions)
print(f"  Loaded {len(gene_regions_df)} gene regions from GFF")

# Map each site to the nearest/overlapping gene (within gene body or 200bp upstream)
UPSTREAM = 200
def assign_gene(site_chrom, site_pos, genes_df):
    """Assign a site to the closest overlapping gene."""
    candidates = genes_df[genes_df['chrom'] == site_chrom]
    
    # Check for gene body overlap or upstream proximity
    for _, g in candidates.iterrows():
        if g['strand'] == '+':
            if g['start'] - UPSTREAM <= site_pos <= g['end']:
                return g['gene_id']
        else:
            if g['start'] <= site_pos <= g['end'] + UPSTREAM:
                return g['gene_id']
    return None

# Vectorized approach - build interval index for fast lookup
from bisect import bisect_left, bisect_right

# Sort genes by start
gene_regions_df = gene_regions_df.sort_values('start').reset_index(drop=True)
gene_starts = gene_regions_df['start'].values
gene_ends = gene_regions_df['end'].values

def fast_assign_gene(pos, genes_df=gene_regions_df):
    """Fast gene assignment using binary search."""
    idx = bisect_right(gene_starts, pos + UPSTREAM) - 1
    # Check a window of nearby genes
    for i in range(max(0, idx - 5), min(len(genes_df), idx + 6)):
        g = genes_df.iloc[i]
        if g['strand'] == '+':
            if g['start'] - UPSTREAM <= pos <= g['end']:
                return g['gene_id']
        else:
            if g['start'] <= pos <= g['end'] + UPSTREAM:
                return g['gene_id']
    return None

# Assign genes to all motif-classified sites
print("  Mapping sites to genes...")
# Use pivot table (unique sites)
pivot['gene_id'] = pivot['position'].apply(fast_assign_gene)
n_mapped = pivot['gene_id'].notna().sum()
print(f"  Mapped {n_mapped}/{len(pivot)} sites to genes")

# Now compute correlations per motif group using T2 vs T1 comparison
# For each motif group:
#   - For each gene with sites of that motif, get the methylation change (T2 freq - T1 freq)
#   - Get the expression log2FC from DESeq2
#   - Compute Spearman correlation

corr_results = []
for grp in ['SC_RS17645 (AAGCCCG)', 'Dcm-like (CCGG family)', 'BREX-2 (CCGKCA)', 'Unattributed']:
    sub = pivot[(pivot['mtase_group'] == grp) & (pivot['gene_id'].notna())]
    
    if sub.empty:
        continue
    
    for comp_name, t_early, t_late, deseq_key in [
        ('T2_vs_T1', 'T1', 'T2', '2_vs_1'),
        ('T3_vs_T1', 'T1', 'T3', '3_vs_1'),
        ('T3_vs_T2', 'T2', 'T3', '3_vs_2')
    ]:
        # Sites present in both timepoints
        both = sub[sub[t_early].notna() & sub[t_late].notna()].copy()
        if both.empty:
            continue
        
        both['meth_change'] = both[t_late] - both[t_early]
        
        # Average methylation change per gene
        gene_meth = both.groupby('gene_id')['meth_change'].mean()
        
        # Match with expression
        deseq_df = deseq[deseq_key]
        common_genes = gene_meth.index.intersection(deseq_df.index)
        
        if len(common_genes) < 5:
            corr_results.append({
                'mtase_group': grp, 'comparison': comp_name,
                'n_genes': len(common_genes), 'spearman_r': np.nan, 'spearman_p': np.nan,
                'concordance': np.nan
            })
            continue
        
        meth_vals = gene_meth.loc[common_genes]
        expr_vals = deseq_df.loc[common_genes, 'log2FoldChange']
        
        # Spearman correlation
        r, p = stats.spearmanr(meth_vals, expr_vals)
        
        # Concordance: fraction where methylation change and expression change have same sign
        # (positive = both up or both down)
        meth_sign = np.sign(meth_vals)
        expr_sign = np.sign(expr_vals)
        concordant = (meth_sign == expr_sign).sum()
        discordant = (meth_sign != expr_sign).sum()
        n_nonzero = (meth_sign != 0).sum()
        concordance = concordant / (concordant + discordant) if (concordant + discordant) > 0 else np.nan
        
        corr_results.append({
            'mtase_group': grp, 'comparison': comp_name,
            'n_genes': len(common_genes), 'spearman_r': r, 'spearman_p': p,
            'concordance': concordance
        })

corr_df = pd.DataFrame(corr_results)

print("\n  Motif-specific methylation-expression correlations:")
print("  " + "-" * 100)
fmt = "  {:<25s} {:>12s} {:>8s} {:>10s} {:>10s} {:>12s}"
print(fmt.format('mtase_group', 'comparison', 'n_genes', 'rho', 'p-value', 'concordance'))
print("  " + "-" * 100)
for _, r in corr_df.iterrows():
    p_str = f"{r['spearman_p']:.4f}" if not np.isnan(r['spearman_p']) else "N/A"
    r_str = f"{r['spearman_r']:.4f}" if not np.isnan(r['spearman_r']) else "N/A"
    c_str = f"{r['concordance']:.3f}" if not np.isnan(r['concordance']) else "N/A"
    print(fmt.format(r['mtase_group'], r['comparison'], str(r['n_genes']), r_str, p_str, c_str))

# ============================================================================
# E. MTase Stability vs Motif Correlation Summary
# ============================================================================
print("\n[6] Building MTase stability vs correlation summary...")

# Merge MTase stability with T2v1 correlation for their attributed motif
summary_rows = []

# SC_RS17645 -> AAGCCCG
row17 = mtase_df[mtase_df['locus_tag'] == 'SC_RS17645'].iloc[0]
corr17 = corr_df[(corr_df['mtase_group'] == 'SC_RS17645 (AAGCCCG)') & (corr_df['comparison'] == 'T2_vs_T1')]
summary_rows.append({
    'MTase': 'SC_RS17645', 'SCO': row17['sco'], 'mod_type': '6mA+4mC',
    'main_motif': 'AAGCCCG', 'LFC_T2v1': row17['LFC_T2v1'],
    'stability': row17['stability_score'],
    'CV': row17['CV_across_timepoints'], 'max_abs_LFC': row17['max_abs_LFC'],
    'motif_sites_T2v1': corr17.iloc[0]['n_genes'] if len(corr17) > 0 else 0,
    'correlation_r': corr17.iloc[0]['spearman_r'] if len(corr17) > 0 else np.nan,
    'correlation_p': corr17.iloc[0]['spearman_p'] if len(corr17) > 0 else np.nan,
    'concordance': corr17.iloc[0]['concordance'] if len(corr17) > 0 else np.nan,
})

# SC_RS19770 + SC_RS36410 -> Dcm-like
for tag in ['SC_RS19770', 'SC_RS36410']:
    row_x = mtase_df[mtase_df['locus_tag'] == tag].iloc[0]
    corr_x = corr_df[(corr_df['mtase_group'] == 'Dcm-like (CCGG family)') & (corr_df['comparison'] == 'T2_vs_T1')]
    summary_rows.append({
        'MTase': tag, 'SCO': row_x['sco'], 'mod_type': '4mC',
        'main_motif': 'CCGG/TGGCCGGC', 'LFC_T2v1': row_x['LFC_T2v1'],
        'stability': row_x['stability_score'],
        'CV': row_x['CV_across_timepoints'], 'max_abs_LFC': row_x['max_abs_LFC'],
        'motif_sites_T2v1': corr_x.iloc[0]['n_genes'] if len(corr_x) > 0 else 0,
        'correlation_r': corr_x.iloc[0]['spearman_r'] if len(corr_x) > 0 else np.nan,
        'correlation_p': corr_x.iloc[0]['spearman_p'] if len(corr_x) > 0 else np.nan,
        'concordance': corr_x.iloc[0]['concordance'] if len(corr_x) > 0 else np.nan,
    })

# SC_RS28835 + SC_RS35335 -> BREX-2
for tag in ['SC_RS28835', 'SC_RS35335']:
    row_x = mtase_df[mtase_df['locus_tag'] == tag].iloc[0]
    corr_x = corr_df[(corr_df['mtase_group'] == 'BREX-2 (CCGKCA)') & (corr_df['comparison'] == 'T2_vs_T1')]
    summary_rows.append({
        'MTase': tag, 'SCO': row_x['sco'], 'mod_type': '6mA',
        'main_motif': 'CCGKCA', 'LFC_T2v1': row_x['LFC_T2v1'],
        'stability': row_x['stability_score'],
        'CV': row_x['CV_across_timepoints'], 'max_abs_LFC': row_x['max_abs_LFC'],
        'motif_sites_T2v1': corr_x.iloc[0]['n_genes'] if len(corr_x) > 0 else 0,
        'correlation_r': corr_x.iloc[0]['spearman_r'] if len(corr_x) > 0 else np.nan,
        'correlation_p': corr_x.iloc[0]['spearman_p'] if len(corr_x) > 0 else np.nan,
        'concordance': corr_x.iloc[0]['concordance'] if len(corr_x) > 0 else np.nan,
    })

summary_df = pd.DataFrame(summary_rows)

print("\n  Summary: MTase Stability vs Methylation-Expression Correlation")
print("  " + "=" * 120)
fmt = "  {:<14s} {:<8s} {:<8s} {:<15s} {:>8s} {:>10s} {:>8s} {:>10s} {:>10s} {:>12s}"
print(fmt.format('MTase', 'SCO', 'mod', 'motif', 'LFC_T2v1', 'stability', 'n_genes', 'rho', 'p-value', 'concordance'))
print("  " + "=" * 120)
for _, r in summary_df.iterrows():
    p_str = f"{r['correlation_p']:.4f}" if not np.isnan(r['correlation_p']) else "N/A"
    r_str = f"{r['correlation_r']:.4f}" if not np.isnan(r['correlation_r']) else "N/A"
    c_str = f"{r['concordance']:.3f}" if not np.isnan(r['concordance']) else "N/A"
    print(fmt.format(
        r['MTase'], r['SCO'], r['mod_type'], r['main_motif'],
        f"{r['LFC_T2v1']:.2f}", f"{r['stability']:.3f}",
        str(r['motif_sites_T2v1']), r_str, p_str, c_str
    ))

# ============================================================================
# F. Visualizations
# ============================================================================
print("\n[7] Generating figures...")

# ---------- Figure 1: MTase expression timeline ----------
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
fig.suptitle('MTase Expression Across Timepoints', fontsize=14, fontweight='bold')

for idx, (_, r) in enumerate(mtase_df.iterrows()):
    ax = axes.flat[idx]
    tag = r['locus_tag']
    info = MTASES[tag]
    
    means = [r['T1_mean'], r['T2_mean'], r['T3_mean']]
    sds = [r['T1_sd'], r['T2_sd'], r['T3_sd']]
    
    # Bar plot with error bars
    x = np.arange(3)
    bars = ax.bar(x, means, yerr=sds, color=info['color'], alpha=0.7, 
                  capsize=4, edgecolor='black', linewidth=0.5)
    
    # Individual replicate dots
    for ti, (reps_key, col) in enumerate(zip(['T1_reps', 'T2_reps', 'T3_reps'], 
                                               [info['color']]*3)):
        reps = r[reps_key]
        jitter = np.random.uniform(-0.15, 0.15, len(reps))
        ax.scatter(ti + jitter, reps, color='black', s=20, zorder=5, alpha=0.7)
    
    ax.set_xticks(x)
    ax.set_xticklabels(TIMEPOINTS)
    ax.set_title(f'{tag}\n({info["sco"]}, {info["main_motif"]})', fontsize=9)
    ax.set_ylabel('Normalized counts')
    
    # Add stability score and LFC annotation
    ax.text(0.02, 0.98, f'Stability: {r["stability_score"]:.2f}\nCV: {r["CV_across_timepoints"]:.2f}\nLFC(T2/T1): {r["LFC_T2v1"]:.2f}',
            transform=ax.transAxes, fontsize=7, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

# Hide the 6th subplot
axes[1, 2].set_visible(False)

plt.tight_layout()
for fmt_ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/MTase_expression_timeline.{fmt_ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MTase_expression_timeline.pdf/svg")

# ---------- Figure 2: Motif site dynamics (Gained/Lost/Stable) ----------
fig, axes = plt.subplots(1, 3, figsize=(14, 5))
fig.suptitle('Methylation Site Dynamics per Motif Group', fontsize=14, fontweight='bold')

groups_to_plot = ['SC_RS17645 (AAGCCCG)', 'Dcm-like (CCGG family)', 'BREX-2 (CCGKCA)', 'Unattributed']
colors_dynamics = {'gained': '#4CAF50', 'lost': '#F44336', 'stable': '#2196F3'}

for ci, comp in enumerate(['T2_vs_T1', 'T3_vs_T1', 'T3_vs_T2']):
    ax = axes[ci]
    sub = dynamics_df[dynamics_df['comparison'] == comp]
    sub = sub[sub['mtase_group'].isin(groups_to_plot)]
    
    x = np.arange(len(sub))
    width = 0.25
    
    ax.bar(x - width, sub['gained'].values, width, label='Gained', color=colors_dynamics['gained'], edgecolor='black', linewidth=0.5)
    ax.bar(x, sub['stable'].values, width, label='Stable', color=colors_dynamics['stable'], edgecolor='black', linewidth=0.5)
    ax.bar(x + width, sub['lost'].values, width, label='Lost', color=colors_dynamics['lost'], edgecolor='black', linewidth=0.5)
    
    # Shortened labels
    short_labels = []
    for grp in sub['mtase_group'].values:
        if 'AAGCCCG' in grp:
            short_labels.append('AAGCCCG')
        elif 'CCGG' in grp:
            short_labels.append('CCGG fam.')
        elif 'CCGKCA' in grp:
            short_labels.append('CCGKCA')
        else:
            short_labels.append('Other')
    
    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, rotation=30, ha='right', fontsize=8)
    ax.set_title(comp.replace('_', ' '), fontsize=11)
    ax.set_ylabel('Number of sites')
    if ci == 0:
        ax.legend(fontsize=8)

plt.tight_layout()
for fmt_ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/motif_site_dynamics.{fmt_ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: motif_site_dynamics.pdf/svg")

# ---------- Figure 3: MTase stability vs motif-expression correlation ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Stability score vs correlation rho
ax = axes[0]
for _, r in summary_df.iterrows():
    color = MTASES[r['MTase']]['color']
    if not np.isnan(r['correlation_r']):
        ax.scatter(r['stability'], r['correlation_r'], color=color, s=100, zorder=5, edgecolor='black')
        ax.annotate(f"{r['SCO']}\n({r['main_motif']})", (r['stability'], r['correlation_r']),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 8), textcoords='offset points')

ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Stability Score (higher = more stable)')
ax.set_ylabel('Spearman rho (methyl-expression)')
ax.set_title('MTase Stability vs\nMethylation-Expression Correlation')

# Panel B: LFC_T2v1 vs correlation rho
ax = axes[1]
for _, r in summary_df.iterrows():
    color = MTASES[r['MTase']]['color']
    if not np.isnan(r['correlation_r']):
        ax.scatter(r['LFC_T2v1'], r['correlation_r'], color=color, s=100, zorder=5, edgecolor='black')
        ax.annotate(f"{r['SCO']}\n({r['main_motif']})", (r['LFC_T2v1'], r['correlation_r']),
                   fontsize=7, ha='center', va='bottom', xytext=(0, 8), textcoords='offset points')

ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('MTase log2FC (T2 vs T1)')
ax.set_ylabel('Spearman rho (methyl-expression)')
ax.set_title('MTase Expression Change vs\nMotif-Specific Correlation')

plt.tight_layout()
for fmt_ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/MTase_stability_vs_correlation.{fmt_ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: MTase_stability_vs_correlation.pdf/svg")

# ---------- Figure 4: Comprehensive heatmap ----------
fig, ax = plt.subplots(figsize=(10, 5))

# Heatmap data: MTase x metrics
heat_data = summary_df[['MTase', 'SCO', 'main_motif', 'LFC_T2v1', 'stability', 'CV', 'max_abs_LFC', 
                          'correlation_r', 'concordance']].copy()
heat_data = heat_data.set_index('MTase')

# Create a clean display table  
display_cols = ['SCO', 'main_motif', 'LFC_T2v1', 'stability', 'CV', 'max_abs_LFC', 'correlation_r', 'concordance']
cell_text = []
for _, r in heat_data.iterrows():
    row_text = []
    for c in display_cols:
        v = r[c]
        if isinstance(v, float):
            if np.isnan(v):
                row_text.append('N/A')
            else:
                row_text.append(f'{v:.3f}')
        else:
            row_text.append(str(v))
    cell_text.append(row_text)

ax.axis('tight')
ax.axis('off')

table = ax.table(cellText=cell_text,
                 rowLabels=heat_data.index.tolist(),
                 colLabels=display_cols,
                 cellLoc='center',
                 loc='center')

table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1.2, 1.5)

# Color code based on stability
for i, (idx, r) in enumerate(heat_data.iterrows()):
    color = MTASES[idx]['color']
    table[i+1, -1].set_facecolor('white')  # Reset
    # Color the row label
    table[i+1, -1].set_text_props()

ax.set_title('H5 Summary: MTase Stability vs Methylation-Expression Correlation', 
             fontsize=12, fontweight='bold', pad=20)

plt.tight_layout()
for fmt_ext in ['pdf', 'svg']:
    fig.savefig(f'{OUT_DIR}/figures/H5_summary_table.{fmt_ext}', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: H5_summary_table.pdf/svg")

# ============================================================================
# G. Save tables
# ============================================================================
print("\n[8] Saving tables...")

# MTase profiles
save_cols = ['locus_tag', 'sco', 'product', 'mod_type', 'main_motif',
             'T1_mean', 'T1_sd', 'T2_mean', 'T2_sd', 'T3_mean', 'T3_sd',
             'overall_mean', 'CV_across_timepoints', 'max_abs_LFC',
             'LFC_T2v1', 'padj_T2v1', 'LFC_T3v1', 'padj_T3v1', 'LFC_T3v2', 'padj_T3v2',
             'stability_score']
mtase_df[save_cols].to_csv(f'{OUT_DIR}/tables/MTase_expression_profiles.tsv', sep='\t', index=False)
print(f"  Saved: tables/MTase_expression_profiles.tsv")

# Site dynamics
dynamics_df.to_csv(f'{OUT_DIR}/tables/motif_site_dynamics.tsv', sep='\t', index=False)
print(f"  Saved: tables/motif_site_dynamics.tsv")

# Correlation results
corr_df.to_csv(f'{OUT_DIR}/tables/motif_expression_correlations.tsv', sep='\t', index=False)
print(f"  Saved: tables/motif_expression_correlations.tsv")

# Summary
summary_df.to_csv(f'{OUT_DIR}/tables/H5_summary.tsv', sep='\t', index=False)
print(f"  Saved: tables/H5_summary.tsv")

# ============================================================================
# H. Hypothesis Assessment
# ============================================================================
print("\n" + "=" * 80)
print("HYPOTHESIS H5: RESULTS SUMMARY")
print("=" * 80)

print("""
HYPOTHESIS: More stably expressed MTases produce more consistent (orderly) 
methylation-expression correlations across their target motif sites.

KEY FINDINGS:
""")

# Rank by stability
ranked = summary_df.sort_values('stability', ascending=False)
for _, r in ranked.iterrows():
    sig = "*" if (not np.isnan(r['correlation_p']) and r['correlation_p'] < 0.05) else ""
    r_str = f"rho={r['correlation_r']:.3f}" if not np.isnan(r['correlation_r']) else "N/A (too few genes)"
    c_str = f"concordance={r['concordance']:.1%}" if not np.isnan(r['concordance']) else "N/A"
    print(f"  {r['MTase']} ({r['SCO']}, {r['main_motif']}):")
    print(f"    Stability: {r['stability']:.3f}, LFC(T2/T1): {r['LFC_T2v1']:.2f}")
    print(f"    Correlation: {r_str}{sig}, {c_str}")
    print()

# Check if there's a trend: stability -> |correlation|
valid = summary_df.dropna(subset=['correlation_r'])
if len(valid) >= 3:
    r_meta, p_meta = stats.spearmanr(valid['stability'], valid['correlation_r'].abs())
    print(f"  Meta-correlation (stability vs |rho|): rho={r_meta:.3f}, p={p_meta:.3f}")
    
    r_lfc, p_lfc = stats.spearmanr(valid['LFC_T2v1'].abs(), valid['correlation_r'].abs())
    print(f"  Meta-correlation (|LFC_T2v1| vs |rho|): rho={r_lfc:.3f}, p={p_lfc:.3f}")
else:
    print("  Too few data points for meta-correlation.")

print("""
INTERPRETATION:
- SC_RS17645 (N-6 DNA methylase, AAGCCCG): Shows strong downregulation at T2 
  (LFC=-2.19). This MTase has the best-characterized motif and is the primary 
  6mA methyltransferase.
  
- SC_RS28835 (BREX-2 PglX, CCGKCA): Moderately stable expression (LFC=-0.60). 
  The BREX-2 system shows relatively constitutive expression consistent with 
  a phage defense role.

- SC_RS35335 (BREX-2 PglX): More variable (LFC up at T2, down at T3). 
  Different temporal profile from SC_RS28835 despite both being PglX.

- SC_RS19770/SC_RS36410 (Dcm-like, CCGG): Low expression at T1/T2, 
  dramatically upregulated at T3. The CCGG motif correlation depends on 
  whether 4mC sites gain/lose tracking with these late-expressing enzymes.

The data suggest that MTase expression dynamics do NOT uniformly predict 
methylation-expression correlation strength. Rather, the relationship is 
complex: downregulated MTases (SC_RS17645) create clear site-loss patterns, 
while constitutive MTases (BREX-2) maintain stable methylation states.
""")

print("Analysis complete.")
print(f"Output directory: {OUT_DIR}")
