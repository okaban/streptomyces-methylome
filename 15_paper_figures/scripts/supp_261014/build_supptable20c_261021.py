# -*- coding: utf-8 -*-
"""SUPERSEDED 2026-11-02 — 後継 build_supptable20c_allread_261102.py を使うこと。

このスクリプトは残してあるが実行は止めてある（下の raise SystemExit）。理由は 2 点:

  1. 入力 `methyl/260102_M145/analysis/pileup/` が 2026-02 の追加フローセルを含まず
     （3-2 / 3-4）、1-1 は不採用の復元版である。両方とも Writing/superseded_inputs.tsv
     に severity BLOCKER で登録されており、check_data_freshness.py が :16 を検出していた。
  2. 信頼度フィルタ版からレベル尺度の量を作っている（canon F7: LEVEL +15.5 pp）。

実測した入力差（2026-11-02、Writing/supptable20c_lineage_diff_261102.tsv）:
3-2 の level_6mA_interior_all は 0.066395 → 0.092847、3-4 は 0.068067 → 0.077040。
追加フローセルの影響を受けない 6 ライブラリは全 8 セルがビット単位で一致する。

--- 以下、旧版の元の説明 ---

補足表 20c: クラスタ内部と隣接フランクの per-library 6mA レベル。

凡例の宣言:
  「クラスタ内部とフランクの per-library 平均 6mA レベルを、全位置で、および
    実験 TSS から 300 bp 以内の位置を除いて計算する。列: timepoint, replicate,
    interior/flank (all positions), interior/flank (away from TSSs)」

入力:
  - 領域座標: SuppTable20a (29 antiSMASH regions, 1-based inclusive)
  - フランク定義: 領域と同じ長さの上流・下流各 1 本（補足表 20b の
    「two flanks of equal length」）。20b の TGGCCGGC/AAGCCCG フランク計数が
    再現するかをゲノム配列から検算してから使う。
  - 6mA: 全ゲノム信頼度フィルタ版パイルアップ
    methyl/260102_M145/analysis/pileup/<lib>_pileup.bed（1-1 のみ
    1-1_pileup_restored_261009.bed。<lib>_pileup.bed は 5.4M 行で切れている）
    mod code 'a'、Nvalid_cov (col 10) > 0 の行のみ。
  - TSS: 18_tss_analyses/jeong2016_all_tss.csv（Jeong 2016 実験 TSS）

レベルの定義: 区間内の測定可能アデニン位置（鎖別）について Nmod/Nvalid_cov を
位置ごとに求め、その単純平均。分母は「その区間でその library が測定できた位置数」。
"""

raise SystemExit(
    "SUPERSEDED 2026-11-02: build_supptable20c_allread_261102.py "
    "(all-read, canonical pileup lineage) wo tsukau koto."
)

import os
import sys

import numpy as np
import pandas as pd

RNA = '/Users/okaban/bioinfo/rna-seq'
EPI = os.path.join(RNA, '11_epigenome_integration')
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')
A6 = '/tmp/a6mA'
GENOME = os.path.join(EPI, 'data/NC_003888.3.fna')
LIBS = [('T1', '1-1'), ('T1', '1-2'), ('T1', '1-3'),
        ('T2', '2-1'), ('T2', '2-3'), ('T2', '2-4'),
        ('T3', '3-2'), ('T3', '3-3'), ('T3', '3-4')]

# ---------------------------------------------------------------- genome
seq = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
GL = len(seq)
comp = str.maketrans('ACGT', 'TGCA')
rc = seq.translate(comp)[::-1]


def motif_positions(motif):
    """両鎖の 1-based 出現開始位置（+ 鎖座標系）。"""
    out = []
    for s, strand in ((seq, '+'), (rc, '-')):
        i = s.find(motif)
        while i != -1:
            p = i + 1 if strand == '+' else GL - i - len(motif) + 1
            out.append(p)
            i = s.find(motif, i + 1)
    return np.array(sorted(out))


# ---------------------------------------------------------------- regions
t20a = pd.read_csv(os.path.join(OUT, 'SuppTable20a_cluster_methylation_expression_261021.tsv'),
                   sep='\t')
t20b = pd.read_csv(os.path.join(OUT, 'SuppTable20b_cluster_motif_occurrence_261021.tsv'),
                   sep='\t')
assert len(t20a) == 29, len(t20a)

regions = []
for r in t20a.itertuples():
    L = int(r.end) - int(r.start) + 1
    regions.append(dict(region=r.region, start=int(r.start), end=int(r.end), length=L,
                        up=(max(1, int(r.start) - L), int(r.start) - 1),
                        dn=(int(r.end) + 1, min(GL, int(r.end) + L))))

# --- verify the flank definition reproduces Supplementary Table 20b -------
checks = []
for motif in ('TGGCCGGC', 'AAGCCCG'):
    pos = motif_positions(motif)
    ins, fla = [], []
    for g in regions:
        ins.append(int(((pos >= g['start']) & (pos <= g['end'])).sum()))
        fla.append(int((((pos >= g['up'][0]) & (pos <= g['up'][1])) |
                        ((pos >= g['dn'][0]) & (pos <= g['dn'][1]))).sum()))
    pub_in = t20b[f'{motif}_inside'].tolist()
    pub_fl = t20b[f'{motif}_flank'].tolist()
    checks.append((f'20b {motif} inside counts reproduce', ins == pub_in))
    checks.append((f'20b {motif} flank counts reproduce (equal-length flanks)',
                   fla == pub_fl))
    if fla != pub_fl:
        bad = [(g['region'], a, b) for g, a, b in zip(regions, fla, pub_fl) if a != b]
        checks.append((f'   {motif} flank diffs (region, recomputed, published)', bad))

# ---------------------------------------------------------------- masks
interior = np.zeros(GL + 2, dtype=bool)
flank = np.zeros(GL + 2, dtype=bool)
for g in regions:
    interior[g['start']:g['end'] + 1] = True
    flank[g['up'][0]:g['up'][1] + 1] = True
    flank[g['dn'][0]:g['dn'][1] + 1] = True
# 内部が優先（隣接領域のフランクが別領域の内部に掛かる分は内部として数える）
flank &= ~interior

tss = pd.read_csv(os.path.join(EPI, 'analysis/18_tss_analyses/jeong2016_all_tss.csv'))
near = np.zeros(GL + 2, dtype=bool)
for p in tss.tss_position.astype(int):
    near[max(1, p - 300):min(GL, p + 300) + 1] = True

checks.append(('experimental TSS loaded', len(tss)))
checks.append(('interior bp', int(interior.sum())))
checks.append(('flank bp (interior-exclusive)', int(flank.sum())))
checks.append(('bp within 300 bp of a TSS', int(near.sum())))

# ---------------------------------------------------------------- per library
rows = []
for tp, lib in LIBS:
    d = pd.read_csv(os.path.join(A6, f'{lib}.tsv'), sep='\t', header=None,
                    names=['start0', 'strand', 'nvalid', 'nmod'],
                    dtype={'start0': np.int32, 'strand': 'category',
                           'nvalid': np.int32, 'nmod': np.int32})
    pos = d.start0.values + 1                      # bedMethyl は 0-based
    lvl = d.nmod.values / d.nvalid.values
    m_in, m_fl = interior[pos], flank[pos]
    m_far = ~near[pos]
    rows.append(dict(
        timepoint=tp, replicate=lib,
        n_positions_interior=int(m_in.sum()),
        n_positions_flank=int(m_fl.sum()),
        level_6mA_interior_all=round(float(lvl[m_in].mean()), 6),
        level_6mA_flank_all=round(float(lvl[m_fl].mean()), 6),
        n_positions_interior_noTSS=int((m_in & m_far).sum()),
        n_positions_flank_noTSS=int((m_fl & m_far).sum()),
        level_6mA_interior_noTSS=round(float(lvl[m_in & m_far].mean()), 6),
        level_6mA_flank_noTSS=round(float(lvl[m_fl & m_far].mean()), 6),
    ))

t20c = pd.DataFrame(rows)[[
    'timepoint', 'replicate',
    'n_positions_interior', 'level_6mA_interior_all',
    'n_positions_flank', 'level_6mA_flank_all',
    'n_positions_interior_noTSS', 'level_6mA_interior_noTSS',
    'n_positions_flank_noTSS', 'level_6mA_flank_noTSS']]
t20c['interior_minus_flank_all'] = (t20c.level_6mA_interior_all
                                    - t20c.level_6mA_flank_all).round(6)
t20c['interior_minus_flank_noTSS'] = (t20c.level_6mA_interior_noTSS
                                      - t20c.level_6mA_flank_noTSS).round(6)
p = os.path.join(OUT, 'SuppTable20c_cluster_interior_6mA_noTSS_261021.tsv')
t20c.to_csv(p, sep='\t', index=False)

print(t20c.to_string(index=False))
print()
for k, v in checks:
    print(f'{k}: {v}')
print()
# 凡例の主張: 「TSS 近傍位置を除いても interior/flank 差は生じない」
d_all = t20c.interior_minus_flank_all
d_no = t20c.interior_minus_flank_noTSS
from scipy import stats
print('interior - flank, all positions : median %+.4f, range %+.4f..%+.4f'
      % (d_all.median(), d_all.min(), d_all.max()))
print('interior - flank, away from TSS : median %+.4f, range %+.4f..%+.4f'
      % (d_no.median(), d_no.min(), d_no.max()))
print('Wilcoxon signed-rank (all positions, n=9) p = %.4g'
      % stats.wilcoxon(t20c.level_6mA_interior_all, t20c.level_6mA_flank_all).pvalue)
print('Wilcoxon signed-rank (away from TSS, n=9) p = %.4g'
      % stats.wilcoxon(t20c.level_6mA_interior_noTSS, t20c.level_6mA_flank_noTSS).pvalue)
print('\nwritten:', os.path.basename(p))
sys.exit(0)
