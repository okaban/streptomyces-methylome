# -*- coding: utf-8 -*-
"""補足表 19a: プロモーター窓 × 位置集合ごとの、検定に適格な遺伝子数。

凡例の宣言（Supplementary-Materials-List.md L350-351）:
  「各プロモーター窓・位置集合について、窓内にメチル化可能位置を 1 つ以上持つ
    遺伝子数、すなわちプロモーターメチル化効果が測定可能な遺伝子数。
    列: window, position set, number of genes, number eligible,
    percentage eligible, percentage ineligible」

投稿ディレクトリの旧 19a は上の列を 1 つも持たず、遺伝子集合の構築手順
（step/criterion/n/source）を載せた別物だった。ここで凡例どおりに作り直す。

遺伝子の普遍集合: 制御遺伝子 1,055 のうち TSS が割り当たる 1,051
  （52_shielded_exposed_boundary/tables/all_genes_features_unified_n57.tsv の
    locus_tag / tss / strand のみを使う。同表の is_exposed 系の分類列は
    撤回済みなので読まない）

位置集合: 鎖指定のメチル化可能位置をゲノム配列から直接数え上げる
  （パイルアップの検出結果ではなく、配列上の「マークが載りうる位置」）
"""
import os

import numpy as np
import pandas as pd

RNA = '/Users/okaban/bioinfo/rna-seq'
EPI = os.path.join(RNA, '11_epigenome_integration')
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')
GENOME = os.path.join(EPI, 'data/NC_003888.3.fna')

seq = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
GL = len(seq)
rc = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def motif_offset_positions(motif, offsets):
    """モチーフの指定オフセットの、鎖指定 1-based + 鎖座標。"""
    out = []
    for s, strand in ((seq, '+'), (rc, '-')):
        i = s.find(motif)
        while i != -1:
            for off in offsets:
                j = i + off
                out.append(j + 1 if strand == '+' else GL - j)
            i = s.find(motif, i + 1)
    return np.array(sorted(set(out)))


POSITION_SETS = [
    ('GCCGGC hexamer cytosine (offset 2, both strands)',
     motif_offset_positions('GCCGGC', [2])),
    ('TGGCCGGC 4mC position 5 (recognition sequence)',
     motif_offset_positions('TGGCCGGC', [4])),
    ('AAGCCCG 6mA positions A1/A2 (offsets 0,1)',
     motif_offset_positions('AAGCCCG', [0, 1])),
    ('AAGCCCG 4mC position 5 (offset 4)',
     motif_offset_positions('AAGCCCG', [4])),
    ('Either recognition sequence (TGGCCGGC 4mC or any AAGCCCG position)',
     np.array(sorted(set(motif_offset_positions('TGGCCGGC', [4]).tolist()
                         + motif_offset_positions('AAGCCCG', [0, 1, 4]).tolist())))),
]

genes = pd.read_csv(os.path.join(
    EPI, 'analysis/52_shielded_exposed_boundary/tables/'
         'all_genes_features_unified_n57.tsv'), sep='\t',
    usecols=['locus_tag', 'tss', 'strand'])
n_all = len(genes)
g = genes.dropna(subset=['tss']).copy()
g['tss'] = g['tss'].astype(int)
n_genes = len(g)

WINDOWS = [
    ('TSS +/- 100 bp', lambda t, s: (t - 100, t + 100)),
    ('TSS upstream 293 bp (-293..-1)',
     lambda t, s: (t - 293, t - 1) if s == '+' else (t + 1, t + 293)),
    ('TSS +/- 2 kb', lambda t, s: (t - 2000, t + 2000)),
]

rows = []
for wname, wfun in WINDOWS:
    lo = np.empty(n_genes, dtype=np.int64)
    hi = np.empty(n_genes, dtype=np.int64)
    for k, r in enumerate(g.itertuples()):
        a, b = wfun(r.tss, r.strand)
        lo[k], hi[k] = max(1, a), min(GL, b)
    for pname, pos in POSITION_SETS:
        li = np.searchsorted(pos, lo, side='left')
        ri = np.searchsorted(pos, hi, side='right')
        elig = int(((ri - li) > 0).sum())
        rows.append(dict(
            window=wname,
            position_set=pname,
            n_genes=n_genes,
            n_eligible=elig,
            pct_eligible=round(elig / n_genes * 100, 1),
            pct_ineligible=round((n_genes - elig) / n_genes * 100, 1),
        ))

t19a = pd.DataFrame(rows)
p = os.path.join(OUT, 'SuppTable19a_eligibility_for_promoter_test_261021.tsv')
t19a.to_csv(p, sep='\t', index=False)

print('regulatory genes in table: %d; with an assignable TSS: %d' % (n_all, n_genes))
for pname, pos in POSITION_SETS:
    print('  position set "%s": %d strand-oriented positions' % (pname, len(pos)))
print()
print(t19a.to_string(index=False))
print('\nwritten:', os.path.basename(p))
