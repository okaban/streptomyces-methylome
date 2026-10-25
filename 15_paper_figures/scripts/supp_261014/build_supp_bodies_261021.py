# -*- coding: utf-8 -*-
"""補足表 1/2/3/4/6/7/8/9/10/11/12 の実体を投稿ディレクトリに書き出す。

凡例 (Writing/Supplementary-Materials-List.md) が宣言する列名・列順に揃え、
凡例が本文で主張している数値をその場で検算する。検算が落ちたら例外で止める
（合わせ込みはしない。食い違いは例外メッセージとして報告する）。

出力: 15_paper_figures/tables/supplementary/SuppTable<N>_<desc>_261021.tsv
"""
import os
import sys

import pandas as pd

RNA = '/Users/okaban/bioinfo/rna-seq'
EPI = os.path.join(RNA, '11_epigenome_integration/analysis')
RECOMP = os.path.join(RNA, '15_paper_figures/scripts/recompute_261015')
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')

checks = []


def ck(label, got, want, tol=0):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    checks.append((label, got, want, 'OK' if ok else 'MISMATCH'))
    return ok


def rd(path, **kw):
    sep = ',' if path.endswith('.csv') else '\t'
    return pd.read_csv(path, sep=sep, **kw)


def write(df, name):
    p = os.path.join(OUT, name)
    df.to_csv(p, sep='\t', index=False)
    return p


written = []

# ---------------------------------------------------------------- Table 1
t1 = rd(os.path.join(RECOMP, 'site_coloc_2x2_canonical_261015.tsv'))
r_pool = t1[t1.scope == 'pooled'].iloc[0]
r_t1 = t1[t1.scope == 'T1'].iloc[0]
ck('T1 pooled in_comod (legend 475)', int(r_pool.in_comod), 475)
ck('T1 pooled out_comod (legend 0)', int(r_pool.out_comod), 0)
ck('T1 candidate pairs (legend 1,696,280)', int(r_pool.n_candidate_pairs), 1696280)
ck('T1 in-motif pairs (legend 2,668)', int(r_pool.n_in_motif_pairs), 2668)
ck('T1 T1-only in_comod (legend 324)', int(r_t1.in_comod), 324)
ck('T1 pooled OR/1e5 (legend 7.3)', round(float(r_pool.OR) / 1e5, 1), 7.3)
ck('T1 T1-only OR/1e5 (legend 4.7)', round(float(r_t1.OR) / 1e5, 1), 4.7)
written.append(write(t1, 'SuppTable1_AAGCCCG_comodification_2x2_261021.tsv'))

# ---------------------------------------------------------------- Table 2
# 凡例は TGGCCGGC と AAGCCCG の両方を宣言しているが、旧 xlsx の実体は
# GCCGGC 側のみだった。両方を motif 列つきで連結する。
dyn = []
for motif, fn in [('TGGCCGGC', 'A5b_GCCGGC_site_dynamics_10pp.tsv'),
                  ('AAGCCCG', 'A5b_AAGCCCG_site_dynamics_10pp.tsv')]:
    d = rd(os.path.join(EPI, '70_methylation_dynamics/tables', fn))
    d.insert(0, 'motif', motif)
    dyn.append(d)
t2 = pd.concat(dyn, ignore_index=True)
t2 = t2.rename(columns={'delta': 'delta_pp', 'pattern': 'trajectory'})
t2['delta_pp'] = t2['delta_pp'].round(2)
t2 = t2[['motif', 'chrom', 'position', 'strand', 'T1_freq', 'T2_freq', 'T3_freq',
         'delta_pp', 'category', 'trajectory', 'region']]
g = t2.groupby('motif')
ck('T2 TGGCCGGC total (legend 855)', int((t2.motif == 'TGGCCGGC').sum()), 855)
ck('T2 TGGCCGGC dynamic (legend 243)',
   int(((t2.motif == 'TGGCCGGC') & (t2.category == 'dynamic')).sum()), 243)
ck('T2 AAGCCCG total (legend 309)', int((t2.motif == 'AAGCCCG').sum()), 309)
ck('T2 AAGCCCG dynamic (legend 97)',
   int(((t2.motif == 'AAGCCCG') & (t2.category == 'dynamic')).sum()), 97)
ck('T2 TGGCCGGC "Other" trajectory (legend 1)',
   int(((t2.motif == 'TGGCCGGC') & (t2.trajectory == 'Other')).sum()), 1)
ck('T2 min delta among dynamic >= 10 pp',
   bool((t2.loc[t2.category == 'dynamic', 'delta_pp'] >= 10).all()), True)
written.append(write(t2, 'SuppTable2_per_site_dynamics_261021.tsv'))

# ---------------------------------------------------------------- Table 3
t3 = rd(os.path.join(EPI, '21_genuswide_motif_conservation/motif_site_density_genuswide.csv'))
ck('T3 assemblies (legend 833)', len(t3), 833)
motifs3 = ['CCGG', 'AAGCCCG', 'GATC', 'GGCCGG', 'GCCCG', 'GCGG', 'CGACNNNCTCC']
ck('T3 motif blocks (legend 7)',
   sum(all(f'{m}_{s}' in t3.columns for s in ('count', 'density', 'expected', 'oe'))
       for m in motifs3), 7)
ck('T3 has no O/E classification column',
   bool(not any('class' in c.lower() for c in t3.columns)), True)
ck('T3 has no REBASE column',
   bool(not any('rebase' in c.lower() for c in t3.columns)), True)
written.append(write(t3, 'SuppTable3_streptomyces_motif_oe_261021.tsv'))

# ---------------------------------------------------------------- Table 4
cat = rd(os.path.join(EPI, '33_GCCGGC_RM_identification/tables/M145_MTase_complete_catalog.tsv'))
ck('T4 MTase genes (legend 129)', len(cat), 129)
for tp in ('T1', 'T2', 'T3'):
    ck(f'T4 {tp}_mean == {tp}_mean_counts',
       bool((cat[f'{tp}_mean'] - cat[f'{tp}_mean_counts']).abs().max() < 1e-9), True)
ck('T4 gene_id duplicates locus_tag', bool((cat.gene_id == cat.locus_tag).all()), True)
t4 = cat.rename(columns={'T1_mean_counts': 'expression_T1',
                         'T2_mean_counts': 'expression_T2',
                         'T3_mean_counts': 'expression_T3',
                         'category': 'substrate_category'})
t4 = t4[['locus_tag', 'product', 'start', 'end', 'strand',
         'expression_T1', 'expression_T2', 'expression_T3',
         'T2vsT1_log2FC', 'T2vsT1_padj', 'T3vsT1_log2FC', 'T3vsT1_padj',
         'substrate_category']]
for c in ('expression_T1', 'expression_T2', 'expression_T3'):
    t4[c] = t4[c].round(2)
for c in ('T2vsT1_log2FC', 'T3vsT1_log2FC'):
    t4[c] = t4[c].round(4)
ck('T4 substrate categories are the 4 legend classes',
   sorted(t4.substrate_category.dropna().unique().tolist()),
   sorted(['DNA cytosine MTase', 'DNA adenine MTase', 'RNA MTase',
           'SAM-dependent MTase (unclassified)']))
for lt in ('SC_RS12890', 'SC_RS19770'):
    ck(f'T4 contains {lt} (named in legend)', bool((t4.locus_tag == lt).any()), True)
written.append(write(t4, 'SuppTable4_MTase_catalogue_261021.tsv'))

# ---------------------------------------------------------------- Table 6
t6 = rd(os.path.join(EPI, '80_rebase_dualmod_search/tables/single_enzyme_dual.tsv'))
ck('T6 dual-mod enzyme records (legend 51)', len(t6), 51)
ck('T6 opposite-strand (legend 50)', int((t6.strand == 'opposite-strand').sum()), 50)
ck('T6 same-strand (legend 1)', int((t6.strand == 'same-strand').sum()), 1)
ck('T6 AAGCCCG absent from REBASE', bool(not t6.motif.str.contains('AAGCCCG').any()), True)
written.append(write(t6, 'SuppTable6_REBASE_dual_modification_261021.tsv'))

# ---------------------------------------------------------------- Table 7
t7 = rd(os.path.join(EPI, '92_d6_5mC_misassignment_control/tables/'
                          'SuppTable11_5mC_misassignment_control.tsv'))
ck('T7 true-5mC observations (legend 18,010)',
   int(t7.loc[t7['sample'] == 'true_5mC', 'n_obs'].iloc[0]), 18010)
ck('T7 unmodified-control observations (legend 8,925)',
   int(t7.loc[t7['sample'] == 'canonical_C', 'n_obs'].iloc[0]), 8925)
ck('T7 CCGG stratum true_5mC obs (legend 334)',
   int(t7.loc[t7['sample'] == 'true_5mC_CCGG_context', 'n_obs'].iloc[0]), 334)
ck('T7 CCGG stratum control obs (legend 51)',
   int(t7.loc[t7['sample'] == 'canonical_C_CCGG_context', 'n_obs'].iloc[0]), 51)
written.append(write(t7, 'SuppTable7_5mC_misassignment_control_261021.tsv'))

# ---------------------------------------------------------------- Table 8
t8 = rd(os.path.join(EPI, '89_supp_deliverables/SuppTable13_AAGCCCG_core_enrichment.tsv'))
null8 = t8[t8.scope == 'motif-conditioned null'].iloc[0]
ck('T8 genomic AAGCCCG instances (legend 1,334)', int(null8.n), 1334)
ck('T8 instances in core (legend 930)', int(null8.n_core), 930)
ck('T8 null core fraction (legend 0.697)', round(float(null8.core_fraction), 3), 0.697)
ck('T8 930/1334 recomputes to the stored fraction',
   round(930 / 1334, 4), round(float(null8.core_fraction), 4), tol=1e-4)
obs8 = t8[t8.scope == 'observed']
ck('T8 6mA not enriched at any timepoint (all p > 0.05)',
   bool((obs8.loc[obs8['mark'] == '6mA', 'p_two_sided'].astype(float) > 0.05).all()), True)
ck('T8 4mC enriched T1 (p <= 1e-4)',
   bool(float(obs8[(obs8['mark'] == '4mC') & (obs8.timepoint == 'T1')].p_two_sided.iloc[0])
        <= 1e-4), True)
ck('T8 4mC not enriched T2 (p > 0.05)',
   bool(float(obs8[(obs8['mark'] == '4mC') & (obs8.timepoint == 'T2')].p_two_sided.iloc[0])
        > 0.05), True)
written.append(write(t8, 'SuppTable8_AAGCCCG_core_enrichment_261021.tsv'))

# ---------------------------------------------------------------- Table 9
t9 = rd(os.path.join(EPI, '89_supp_deliverables/SuppTable14_hemi_strand_mechanism.tsv'))
ck('T9 fixed-strand occurrences (legend 1,514)',
   int(t9[t9.test.str.startswith('replichore x')].n.iloc[0]), 1514)
ck('T9 replichore strata sum to 1,514',
   int(t9[t9.test.str.startswith('plus-strand share')].n.sum()), 1514)
ck('T9 sense/antisense scope n (legend 999)',
   int(t9[t9.test.str.startswith('sense vs antisense')].n.iloc[0]), 999)
written.append(write(t9, 'SuppTable9_hemimethylation_mechanism_261021.tsv'))

# ---------------------------------------------------------------- Table 10
t10 = rd(os.path.join(EPI, '89_supp_deliverables/SuppTable15_growth_and_production_raw.tsv'))
ck('T10 measurements (legend 36)', len(t10), 36)
ck('T10 replicates per group (legend 4)',
   sorted(t10.groupby(['measurement', 'timepoint']).size().unique().tolist()), [4])
ck('T10 timepoints in hours (legend 12/24/50)',
   sorted(t10.hours.unique().tolist()), [12, 24, 50])
ck('T10 measurement types (DCW + 2 antibiotics)', t10.measurement.nunique(), 3)
written.append(write(t10, 'SuppTable10_growth_and_production_261021.tsv'))

# ---------------------------------------------------------------- Table 11
t11 = rd(os.path.join(RECOMP, 'SuppTable15_filter_bias_per_library.tsv'))
ck('T11 libraries (9 Nanopore libraries)', len(t11), 9)
ck('T11 bias = filtered - all-reads (pp), recomputed',
   bool(((t11.occupancy_filtered_pct - t11.occupancy_all_reads_pct)
         - t11.filter_bias_pp).abs().max() <= 0.011), True)
ck('T11 occupancy_filtered = sites/measurable, recomputed',
   bool(((t11.sites_filtered / t11.measurable_filtered * 100)
         - t11.occupancy_filtered_pct).abs().max() <= 0.011), True)
ck('T11 all-read denominator is the full motif population (2,093)',
   sorted(t11.measurable_all_reads.unique().tolist()), [2093])
written.append(write(t11, 'SuppTable11_filter_bias_per_library_261021.tsv'))

# ---------------------------------------------------------------- Table 12
t12 = rd(os.path.join(RECOMP, 'SuppTable16_unassigned_6mA_261015.tsv'))
corr = rd(os.path.join(RECOMP, 'SuppTable16_unassigned_6mA_corrected_261015.tsv'))
ck('T12 plain and "corrected" files are value-identical',
   bool(t12.equals(corr)), True)
r = t12.set_index('timepoint')
ck('T12 T1 unassigned 6mA count (legend 1,531)', int(r.loc['T1', 'mA_unassigned']), 1531)
ck('T12 T1 total 6mA sites (legend 1,947)', int(r.loc['T1', 'sites_6mA']), 1947)
ck('T12 1531/1947 recomputes to 78.63%', round(1531 / 1947 * 100, 2), 78.63, tol=0.01)
ck('T12 unassigned 6mA pct T1/T2/T3 (legend 78.63/79.07/80.37)',
   [float(r.loc[t, 'mA_unassigned_pct']) for t in ('T1', 'T2', 'T3')],
   [78.63, 79.07, 80.37])
ck('T12 unassigned 4mC pct T1/T2/T3 (legend 0.00/0.00/0.04)',
   [float(r.loc[t, 'mC_unassigned_pct']) for t in ('T1', 'T2', 'T3')],
   [0.0, 0.0, 0.04])
written.append(write(t12, 'SuppTable12_genomewide_totals_unassigned_6mA_261021.tsv'))

# ---------------------------------------------------------------- report
bad = [c for c in checks if c[3] != 'OK']
print('written:')
for p in written:
    print('  ', os.path.basename(p))
print('\nlegend cross-checks: %d OK / %d MISMATCH' % (len(checks) - len(bad), len(bad)))
if bad:
    print('NOTE: this script exits 1 when any legend cross-check mismatches. '
          'All %d files above were written before this point; exit 1 is the '
          'mismatch signal, not a write failure.' % len(written))
for label, got, want, st in checks:
    if st != 'OK':
        print('  MISMATCH %s: got %r, legend %r' % (label, got, want))
sys.exit(1 if bad else 0)
