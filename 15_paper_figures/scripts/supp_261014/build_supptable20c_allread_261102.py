#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""補足表 20c（再生成・all-read 系統）: クラスタ内部と隣接フランクの per-library 6mA レベル。

    python3 build_supptable20c_allread_261102.py

なぜ作り直したか
----------------
旧 `build_supptable20c_261021.py` は 2 点で正準から外れていた。

1. **入力が 2026-02 の追加フローセルを含まない。** 旧版 :15-16 は
   `methyl/260102_M145/analysis/pileup/<lib>_pileup.bed` を読む。このディレクトリは
   3-2 / 3-4 が追加フローセル前であり、1-1 は不採用の復元版である
   （`Writing/superseded_inputs.tsv` に 2 件とも severity BLOCKER で登録）。
   `Writing/check_data_freshness.py` がこの参照を BLOCKER として検出していた。
2. **信頼度フィルタ版からレベル尺度の量を作っている。** canon F7 により modkit の
   既定フィルタは LEVEL を +15.5 pp 押し上げる。内部／フランクのレベル比較は
   レベル尺度の量なので、フィルタ版を使ってはならない。

2026-11-01 に補足表 20a を all-read で作り直した
`build_supptable20a_allread_261101.py` と**同じ系統**に揃える。

系統の宣言（20a all-read と同一）
--------------------------------
    入力     : analysis/pileup_nofilt_261015/<lib>_nofilt.bed
    由来     : analysis/<lib>_mapped.bam（3-2 / 3-4 は 2026-02 追加フローセルを merge 済み）
    modkit   : v0.6.4、--no-filtering（Nfail == 0 をアサートする）
    6mA 位置 : AAGCCCG の A1 と A2（`a,AAGCCCG,0` / `a,AAGCCCG,1`、鎖指定 2,668 位置）

**旧版からの定義変更（1 点、意図的）**: 旧版はゲノム全体のアデニン（mod code 'a'、
Nvalid > 0）を母集団としていた（1 ライブラリ 35-38 万の内部位置）。all-read
パイルアップはモチーフ限定で作られておりゲノム全体のアデニンを含まないため、
母集団を 20a と同じ AAGCCCG A1/A2 に揃えた。canon の「このゲノムの 6mA は
AAGCCCG に尽きる」（モチーフ内 15.64 % 対モチーフ外 0.0635 %、246 倍）に従えば
こちらが測るべき量だが、**レベルの絶対値は約 5 倍になる**（背景アデニンの希釈が
外れるため）。旧版との差分は `Writing/supptable20c_lineage_diff_261102.tsv` に分解してある。

測定可能性・レベルの定義は旧版のまま: 区間内の測定可能位置（鎖別、Nvalid > 0）
について Nmod/Nvalid を位置ごとに求め、その単純平均。メチル化閾値はかけない。
フランクは領域と同じ長さの上流・下流各 1 本（内部優先）。TSS 除外は
Jeong 2016 実験 TSS から ±300 bp。

出力: 15_paper_figures/tables/supplementary/
        SuppTable20c_cluster_interior_6mA_allread_261102.tsv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

RNA = '/Users/okaban/bioinfo/rna-seq'
EPI = os.path.join(RNA, '11_epigenome_integration')
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')
NOFILT = '/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015'
GENOME = os.path.join(EPI, 'data/NC_003888.3.fna')
LIBS = [('T1', '1-1'), ('T1', '1-2'), ('T1', '1-3'),
        ('T2', '2-1'), ('T2', '2-3'), ('T2', '2-4'),
        ('T3', '3-2'), ('T3', '3-3'), ('T3', '3-4')]
BEDCOL = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'ts', 'te',
          'color', 'Nvalid', 'pct', 'Nmod', 'Ncanon', 'Nother', 'Ndel',
          'Nfail', 'Ndiff', 'Nnocall']
A6_NAMES = ['a,AAGCCCG,0', 'a,AAGCCCG,1']

seq = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
GL = len(seq)
rc = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def motif_positions(motif):
    """両鎖の 1-based 出現開始位置（+ 鎖座標系）。"""
    out = []
    for s, strand in ((seq, '+'), (rc, '-')):
        i = s.find(motif)
        while i != -1:
            out.append(i + 1 if strand == '+' else GL - i - len(motif) + 1)
            i = s.find(motif, i + 1)
    return np.array(sorted(out))


def build_masks(checks):
    t20a = pd.read_csv(os.path.join(
        OUT, 'SuppTable20a_cluster_methylation_expression_allread_261101.tsv'), sep='\t')
    t20b = pd.read_csv(os.path.join(
        OUT, 'SuppTable20b_cluster_motif_occurrence_261021.tsv'), sep='\t')
    assert len(t20a) == 29, len(t20a)

    regions = []
    for r in t20a.itertuples():
        L = int(r.end) - int(r.start) + 1
        regions.append(dict(region=r.region, start=int(r.start), end=int(r.end),
                            up=(max(1, int(r.start) - L), int(r.start) - 1),
                            dn=(int(r.end) + 1, min(GL, int(r.end) + L))))

    for motif in ('TGGCCGGC', 'AAGCCCG'):          # 20b の計数が再現するか検算
        pos = motif_positions(motif)
        ins = [int(((pos >= g['start']) & (pos <= g['end'])).sum()) for g in regions]
        fla = [int((((pos >= g['up'][0]) & (pos <= g['up'][1])) |
                    ((pos >= g['dn'][0]) & (pos <= g['dn'][1]))).sum()) for g in regions]
        checks.append((f'20b {motif} inside counts reproduce',
                       ins == t20b[f'{motif}_inside'].tolist()))
        checks.append((f'20b {motif} flank counts reproduce (equal-length flanks)',
                       fla == t20b[f'{motif}_flank'].tolist()))

    interior = np.zeros(GL + 2, dtype=bool)
    flank = np.zeros(GL + 2, dtype=bool)
    for g in regions:
        interior[g['start']:g['end'] + 1] = True
        flank[g['up'][0]:g['up'][1] + 1] = True
        flank[g['dn'][0]:g['dn'][1] + 1] = True
    flank &= ~interior                              # 内部が優先

    tss = pd.read_csv(os.path.join(EPI, 'analysis/18_tss_analyses/jeong2016_all_tss.csv'))
    near = np.zeros(GL + 2, dtype=bool)
    for p in tss.tss_position.astype(int):
        near[max(1, p - 300):min(GL, p + 300) + 1] = True
    checks += [('experimental TSS loaded', len(tss)),
               ('interior bp', int(interior.sum())),
               ('flank bp (interior-exclusive)', int(flank.sum())),
               ('bp within 300 bp of a TSS', int(near.sum()))]
    return interior, flank, near


def read_allread_6mA(lib):
    d = pd.read_csv(f'{NOFILT}/{lib}_nofilt.bed', sep='\t', header=None, names=BEDCOL)
    assert int(d['Nfail'].sum()) == 0, f'{lib}: not an all-read pileup'
    d = d[d['name'].isin(A6_NAMES)].drop_duplicates(subset=['start', 'strand'])
    # ゲノム上の AAGCCCG A1/A2 は鎖指定 2,668 位置。modkit は被覆 0 の位置を行として
    # 出さないため、ライブラリにより 2,664-2,668 行になる（実測 1-1/1-2 2,664、
    # 1-3 2,666、他 6 本 2,668。Nvalid == 0 の行は 1 本も無い）。
    assert 2600 <= len(d) <= 2668, f'{lib}: AAGCCCG A1/A2 positions out of range ({len(d)})'
    d = d[d['Nvalid'] > 0]
    return d['start'].to_numpy() + 1, (d['Nmod'] / d['Nvalid']).to_numpy()


def main():
    checks = []
    interior, flank, near = build_masks(checks)

    rows = []
    for tp, lib in LIBS:
        pos, lvl = read_allread_6mA(lib)
        m_in, m_fl, m_far = interior[pos], flank[pos], ~near[pos]
        rows.append(dict(
            timepoint=tp, replicate=lib,
            n_positions_interior=int(m_in.sum()),
            level_6mA_interior_all=round(float(lvl[m_in].mean()), 6),
            n_positions_flank=int(m_fl.sum()),
            level_6mA_flank_all=round(float(lvl[m_fl].mean()), 6),
            n_positions_interior_noTSS=int((m_in & m_far).sum()),
            level_6mA_interior_noTSS=round(float(lvl[m_in & m_far].mean()), 6),
            n_positions_flank_noTSS=int((m_fl & m_far).sum()),
            level_6mA_flank_noTSS=round(float(lvl[m_fl & m_far].mean()), 6),
            n_positions_measurable_AAGCCCG_A1A2=int(len(pos)),
        ))

    t = pd.DataFrame(rows)
    t['interior_minus_flank_all'] = (t.level_6mA_interior_all
                                     - t.level_6mA_flank_all).round(6)
    t['interior_minus_flank_noTSS'] = (t.level_6mA_interior_noTSS
                                       - t.level_6mA_flank_noTSS).round(6)
    p = os.path.join(OUT, 'SuppTable20c_cluster_interior_6mA_allread_261102.tsv')
    t.to_csv(p, sep='\t', index=False)

    print(t.to_string(index=False))
    print()
    for k, v in checks:
        print(f'{k}: {v}')
    d_all, d_no = t.interior_minus_flank_all, t.interior_minus_flank_noTSS
    print('\ninterior - flank, all positions : median %+.4f, range %+.4f..%+.4f'
          % (d_all.median(), d_all.min(), d_all.max()))
    print('interior - flank, away from TSS : median %+.4f, range %+.4f..%+.4f'
          % (d_no.median(), d_no.min(), d_no.max()))
    print('Wilcoxon signed-rank (all positions, n=9) p = %.4g'
          % stats.wilcoxon(t.level_6mA_interior_all, t.level_6mA_flank_all).pvalue)
    print('Wilcoxon signed-rank (away from TSS, n=9) p = %.4g'
          % stats.wilcoxon(t.level_6mA_interior_noTSS, t.level_6mA_flank_noTSS).pvalue)
    print('\nwritten:', os.path.basename(p))
    return 0


if __name__ == '__main__':
    sys.exit(main())
