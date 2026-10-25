# -*- coding: utf-8 -*-
"""補足表 19b: 各解析集合で 80% 検出力・α = 0.05 で検出可能な最小効果量。

投稿ディレクトリの旧 19b は per-read 共修飾の行の design セルに日本語の内部注記
（「周辺度数が本文に宣言されていないため再計算していない（著者確認）」）が残って
いた。ここで周辺度数を正準抽出から取り直し、最小検出オッズ比を実際に計算し、
design 列を英語の検定設計記述に置き換える。

per-read 2x2 の再抽出:
  入力  11_epigenome_integration/analysis/111_perread_canonical_261015/calls/{1-1,1-2,1-3}_calls.tsv
        （modkit extract calls --no-filtering、AAGCCCG の 5 オフセット両鎖）
  分母  HAAGCCCG インスタンス（5' 隣接塩基が G 以外）を、1 本の read が
        オフセット 0（position 1、6mA）とオフセット 4（position 5、4mC）の
        両方で覆っている read x instance ペア
  判定  単分子の確率閾値 0.5（従来の共修飾解析と同じ。部位コールの
        A 0.6387 / C 0.9375 とは別物）
"""
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.power import NormalIndPower

RNA = '/Users/okaban/bioinfo/rna-seq'
EPI = os.path.join(RNA, '11_epigenome_integration')
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')
CALLS = os.path.join(EPI, 'analysis/111_perread_canonical_261015/calls')
GENOME = os.path.join(EPI, 'data/NC_003888.3.fna')
THR = 0.5

# ---- HAAGCCCG instance set (5' flanking base != G) -----------------------
seq = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
GL = len(seq)
rc = seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

inst = {}   # (strand, offset0 ref_position 0-based) -> instance id
five = []
for s, strand in ((seq, '+'), (rc, '-')):
    i = s.find('AAGCCCG')
    while i != -1:
        b5 = s[i - 1] if i > 0 else 'N'
        p0 = i if strand == '+' else GL - i - 7      # 0-based + strand coord of offset 0
        inst[(strand, p0 if strand == '+' else p0 + 6)] = (strand, p0, b5)
        five.append(b5)
        i = s.find('AAGCCCG', i + 1)
five = pd.Series(five)
n_inst = len(inst)
n_H = int((five != 'G').sum())
n_G = int((five == 'G').sum())

# offset -> ref_position map, keyed on (ref_mod_strand, ref_position)
off_map = {}
for (strand, anchor), (st, p0, b5) in inst.items():
    if b5 == 'G':
        continue
    key = (strand, p0, b5)
    for off in (0, 4):
        rp = p0 + off if strand == '+' else p0 + 6 - off
        off_map[(strand, rp)] = (p0, off)

# ---- per-read 2x2 at T1 ---------------------------------------------------
pairs = {}
for lib in ('1-1', '1-2', '1-3'):
    d = pd.read_csv(os.path.join(CALLS, f'{lib}_calls.tsv'), sep='\t',
                    usecols=['read_id', 'ref_position', 'ref_mod_strand',
                             'call_prob', 'call_code'],
                    dtype={'call_code': str})
    key = list(zip(d.ref_mod_strand, d.ref_position))
    hit = np.array([k in off_map for k in key])
    d = d[hit].copy()
    meta = [off_map[k] for k, h in zip(key, hit) if h]
    d['inst'] = [m[0] for m in meta]
    d['off'] = [m[1] for m in meta]
    d['mod'] = np.where(
        d.off == 0,
        (d.call_code == 'a') & (d.call_prob >= THR),
        (d.call_code == '21839') & (d.call_prob >= THR))
    for rid, instp, off, mod in zip(d.read_id, d.inst, d.off, d['mod']):
        k = (lib, rid, instp)
        e = pairs.setdefault(k, [None, None])
        e[0 if off == 0 else 1] = bool(mod)

both = [(a, c) for a, c in pairs.values() if a is not None and c is not None]
n_both = sum(1 for a, c in both if a and c)
n_A = sum(1 for a, c in both if a and not c)
n_C = sum(1 for a, c in both if not a and c)
n_nei = sum(1 for a, c in both if not a and not c)
N = len(both)
table = np.array([[n_both, n_A], [n_C, n_nei]])
OR, p_f = stats.fisher_exact(table)

# ---- minimum detectable odds ratio at 80% power, alpha = 0.05 -------------
# 周辺度数を固定し、4mC 無しの群の 6mA 率 p0 を観測値に置いて、
# 二群の比率差検定（両側 alpha=0.05, power=0.80）で検出可能な p1 を求め、
# それをオッズ比に直す。群の大きさは観測の 4mC あり / 無しの数。
def mde_odds_ratio(n1, n0, p0):
    """周辺度数を固定したときの最小検出オッズ比（両側 alpha=0.05, power=0.80）。"""
    es = float(np.ravel(NormalIndPower().solve_power(
        effect_size=None, nobs1=n0, alpha=0.05, power=0.80,
        ratio=n1 / n0, alternative='two-sided'))[0])
    p1 = np.sin(np.arcsin(np.sqrt(p0)) + es / 2) ** 2
    return (p1 / (1 - p1)) / (p0 / (1 - p0))


n1 = n_both + n_C          # reads with a 4mC call at position 5
n0 = n_A + n_nei           # reads without
p0 = n_A / n0
mde_recount = mde_odds_ratio(n1, n0, p0)

# 本文・レジストリが引用する分母は COMOD.PAIRS.T1 = 39,799（値レジストリ
# value_registry.py / comod_canonical_recompute_261015.csv）。この分母は上の
# 再計数（80,227 ペア）では再現しない。表に載せる MDE は本文が引用している
# 39,799 の側で与え、再計数との食い違いを note 列に明記する。
N_REG = 39799
scale = N_REG / N
mde_reg = mde_odds_ratio(n1 * scale, n0 * scale, p0)

rows = [
    dict(analysis_set='Promoter-proximal 6mA level (+/-100 bp of a TSS vs elsewhere)',
         unit_of_analysis='culture units', n=9, statistic="Cohen's d",
         detectable_effect=1.067,
         design='one-sample t-test, two-sided, alpha = 0.05, power = 0.80'),
    dict(analysis_set='TSS-level ratio measurable in both windows',
         unit_of_analysis='transcription start sites', n=117, statistic="Cohen's d",
         detectable_effect=0.261,
         design='one-sample t-test, two-sided, alpha = 0.05, power = 0.80'),
    dict(analysis_set='Gene-level methylation-expression correlation',
         unit_of_analysis='regulatory genes with an assignable TSS', n=1051,
         statistic='|rho|', detectable_effect=0.0863,
         design='Fisher z-transformation, two-sided, alpha = 0.05, power = 0.80'),
    dict(analysis_set='Per-read co-modification contingency '
                      '(6mA at position 1 vs 4mC at position 5, HAAGCCCG, T1)',
         unit_of_analysis='read x instance pairs covering both positions', n=N_REG,
         statistic='odds ratio', detectable_effect=round(float(mde_reg), 3),
         design='two-proportion z-test on the contingency marginals '
                '(%.1f%% of reads carry a 4mC call at position 5; 6mA rate in reads '
                'without a 4mC call %.4f), two-sided, alpha = 0.05, power = 0.80'
                % (n1 / N * 100, p0)),
]
t19b = pd.DataFrame(rows)
t19b['note'] = ['', '', '',
                'Marginal proportions recounted from the canonical per-read extraction '
                '(111_perread_canonical_261015, offsets 0 and 4, call probability '
                '>= 0.50); that recount yields %d pairs and an odds ratio of %.4f, '
                'against the %d pairs and odds ratio 1.257 recorded for this analysis. '
                'The bound is given on the recorded denominator; on the recounted '
                'denominator it is %.3f. The discrepancy in pair count is unresolved.'
                % (N, OR, N_REG, mde_recount)]
p = os.path.join(OUT, 'SuppTable19b_detectable_effect_sizes_261021.tsv')
t19b.to_csv(p, sep='\t', index=False)

# ---- verify the recomputed MDE rows that were already in the table --------
from statsmodels.stats.power import TTestPower
chk = []
for nobs, pub in ((9, 1.067), (117, 0.261)):
    got = float(np.ravel(TTestPower().solve_power(
        nobs=nobs, alpha=0.05, power=0.80, alternative='two-sided'))[0])
    chk.append(("n=%d Cohen's d MDE (table %s)" % (nobs, pub), round(got, 3)))
z = stats.norm.isf(0.025) + stats.norm.isf(0.20)
chk.append(('n=1051 |rho| MDE via Fisher z (table 0.0863)',
            round(float(np.tanh(z / np.sqrt(1051 - 3))), 4)))

print('AAGCCCG genomic instances: %d (5\' H = %d, 5\' G = %d)' % (n_inst, n_H, n_G))
print("5' base breakdown:", five.value_counts().to_dict())
print('\nT1 HAAGCCCG per-read 2x2 at probability threshold %.1f' % THR)
print('  n_both = %d, n_A_only = %d, n_C_only = %d, n_neither = %d, total = %d'
      % (n_both, n_A, n_C, n_nei, N))
print('  odds ratio = %.4f, Fisher two-sided p = %.3g' % (OR, p_f))
print('  minimum detectable OR (80%% power, alpha = 0.05): %.4f on the recounted '
      'denominator (n = %d), %.4f on the recorded denominator (n = %d)'
      % (mde_recount, N, mde_reg, N_REG))
print('  recorded for this analysis (value_registry.py): n = 39,799, OR = 1.257, '
      'both-mark share 30.36%%; recounted both-mark share %.2f%%'
      % (n_both / N * 100))
print()
for k, v in chk:
    print('  recomputed %s = %s' % (k, v))
print()
print(t19b.to_string(index=False))
print('\nwritten:', os.path.basename(p))
