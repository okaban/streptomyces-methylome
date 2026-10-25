# -*- coding: utf-8 -*-
"""補足表 16（全検定一覧）の n / correction / p_corrected を埋める。

旧 SuppTable16_statistical_test_summary_261016.tsv は 29 行すべてで
n / correction / p_corrected が空だった。Methods がこの表を「全検定の
多重比較補正の記録」として指名しているため、空のままにできない。

方針:
  n            値レジストリ (value_registry.csv) か正準補足表から直接取れるものだけを
               入れる。取れない行は空白ではなく "not recorded" と明記する（推測しない）。
  correction   原稿が族として宣言しているのは「検定可能な 7 モチーフの領域分布検定に
               対する Bonferroni 閾値 0.05/7 = 0.00714」だけ（原稿 Supplementary
               Figure 12 凡例および補足表 26 凡例で明示）。この族に属する行のみ
               Bonferroni とし、それ以外は "none (single test)" と書く。
  source       原稿の行番号は別エージェントが編集中でずれるため、各行の統計量・p 値の
               文字列を原稿本文から検索して、その箇所を含む直近の見出しを記録する。
               旧 source_line_EN は 2026-10-16 時点の値として残す。
"""
import os
import re
import unicodedata

import pandas as pd

RNA = '/Users/okaban/bioinfo/rna-seq'
W = '/Users/okaban/obsidian/Research/rna-seq/Writing'
OUT = os.path.join(RNA, '15_paper_figures/tables/supplementary')
EN = os.path.join(W, 'Manuscript_EN_NAR_numbered_inline.md')

src = pd.read_csv(os.path.join(OUT, 'superseded',
                                   'SuppTable16_statistical_test_summary_261016.tsv'),
                  sep='\t', dtype=str).fillna('')
assert len(src) == 29, len(src)

# ---- 1) n: 出所のある値だけ -------------------------------------------------
# キー = (test, comparison の先頭 40 字) で照合する
N_BY_COMPARISON = {
    # 3 群 x 4 生物学的反復 = 12 測定。F(2,9) の自由度と整合。出所: 補足表 10
    'dry cell weight': ('12', 'Supplementary Table 10 (4 replicates x 3 timepoints)'),
    'actinorhodin': ('12', 'Supplementary Table 10 (4 replicates x 3 timepoints)'),
    'undecylprodigiosin': ('12', 'Supplementary Table 10 (4 replicates x 3 timepoints)'),
    # 候補ペアの全数。出所: 補足表 1
    'same-instance 6mA/4mC double-mark pairs':
        ('1,696,280 candidate same-strand A/C pairs (10,000 permutation replicates)',
         'Supplementary Table 1'),
    # 出所: value_registry DET.6mA.M3G.CONTROL.GAAGCCCG (universe_n)
    'non-target GAAGCCCG class': ('769', 'value_registry.py DET.6mA.M3G.CONTROL.GAAGCCCG'),
    # 出所: value_registry UNASG.6mA.NEARCOGNATE.OR (universe_n)
    'sequence contexts within one substitution':
        ('2,412,520 eligible non-motif reference adenines',
         'value_registry.py UNASG.6mA.NEARCOGNATE.OR'),
    # 出所: value_registry UNASG.6mA.EXPR.PARTIAL (universe_n)
    'unassigned 6mA site density vs transcript level':
        ('7,863 genes', 'value_registry.py UNASG.6mA.EXPR.PARTIAL'),
    # 出所: value_registry REGION.PROM.HAAGCCCG.ALL / .MOD、補足表 26
    'modified vs measurable positions across promoter':
        ('561 measurable HAAGCCCG positions (259 modified)',
         'Supplementary Table 26; value_registry.py REGION.PROM.HAAGCCCG.*'),
    # 出所: value_registry EPIALL.CHI2.P.T1 (denominator = 410 molecules)
    'observed tail distribution vs binomial expectation':
        ('410 molecules', 'value_registry.py EPIALL.CHI2.P.T1'),
}

# ---- 2) correction family --------------------------------------------------
FAMILY_REGION = ('Bonferroni (region-distribution family, 7 testable motifs, '
                 'alpha = 0.05/7 = 0.00714)')
SINGLE = 'none (single test)'


def superscript_to_plain(s):
    return unicodedata.normalize('NFKC', s)


def p_float(s):
    """'4.2 x 10-29' 等の表記から float を取る。取れなければ None。"""
    t = superscript_to_plain(s).replace('\u2212', '-').replace('\u2009', '')
    t = t.replace('×', 'x').replace('<', '').strip()
    m = re.match(r'^\s*([\d.]+)\s*x\s*10(-?\d+)\s*$', t)
    if m:
        return float(m.group(1)) * 10 ** int(m.group(2))
    m = re.match(r'^\s*10(-?\d+)\s*$', t)
    if m:
        return 10.0 ** int(m.group(1))
    m = re.match(r'^\s*p?\s*=?\s*([\d.]+)\s*$', t)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


# ---- 3) 原稿の該当箇所を内容で同定 -------------------------------------------
lines = open(EN, encoding='utf-8').read().split('\n')
heads = [(i, l.lstrip('#').strip()) for i, l in enumerate(lines)
         if re.match(r'^#{1,4} ', l)]


STOP = {'test', 'tests', 'within', 'across', 'versus', 'against', 'between',
        'instance', 'instances', 'position', 'positions', 'distribution'}


def _head_at(i):
    prev = [h for h in heads if h[0] < i]
    return prev[-1][1] if prev else '(before first heading)'


def section_for(comparison, statistic):
    """検定名と比較内容の語で本文を照合し、最も多く一致した行の見出しを返す。

    行番号は別エージェントの編集でずれるので使わない。p 値の文字列は複数箇所に
    現れて誤爆するため、まず statistic 中の固有数値、次に comparison の語で照合する。
    """
    num = re.findall(r'\d+\.\d{2,}', superscript_to_plain(statistic))
    for nd in num:
        hits = [i for i, l in enumerate(lines) if nd in l]
        if len(hits) == 1:
            return _head_at(hits[0])
    # 語の重なりによる推定は誤った見出しを返すことがあり、投稿表に書く値としては
    # 危険なので採用しない。一意に決まらない行は、test 列と comparison 列そのものを
    # 検索キーとして使うよう明記する。
    return 'identify by the test and comparison columns '\
           '(the 2026-10-16 line numbers are stale)'


rows = []
for r in src.itertuples(index=False):
    comp = r.comparison
    n_val, n_src = r.n, ''
    if not n_val:
        for k, (v, s) in N_BY_COMPARISON.items():
            if comp.startswith(k) or k in comp:
                n_val, n_src = v, s
                break
    if not n_val:
        n_val, n_src = 'not recorded', ''

    is_region_family = comp.startswith('modified vs measurable positions')
    pu = p_float(r.p_uncorrected)
    if is_region_family:
        correction = FAMILY_REGION
        p_corr = ('%.4f' % min(1.0, pu * 7)) if pu is not None else 'not computable'
    else:
        correction = SINGLE
        p_corr = r.p_uncorrected if r.p_uncorrected else 'not recorded'

    sec = section_for(comp, r.statistic)
    rows.append(dict(
        test=r.test, comparison=comp, n=n_val, statistic=r.statistic,
        p_uncorrected=r.p_uncorrected, correction=correction, p_corrected=p_corr,
        n_source=n_src,
        source_section=sec,
        source_line_EN_261016=r.source_line_EN.strip(),
    ))

t16 = pd.DataFrame(rows)
p = os.path.join(OUT, 'SuppTable16_statistical_test_summary_261021.tsv')
t16.to_csv(p, sep='\t', index=False)

print('rows: %d | n filled from a source: %d | n "not recorded": %d'
      % (len(t16), int((t16.n_source != '').sum()), int((t16.n == 'not recorded').sum())))
print('correction = Bonferroni family: %d | none (single test): %d'
      % (int((t16.correction == FAMILY_REGION).sum()),
         int((t16.correction == SINGLE).sum())))
print('blank cells remaining in n/correction/p_corrected: %d'
      % int((t16[['n', 'correction', 'p_corrected']] == '').sum().sum()))
print('rows not located in the manuscript by content: %d'
      % int((t16.source_section == 'not located by content').sum()))
print()
print(t16[['test', 'n', 'correction', 'p_corrected', 'source_section']]
      .to_string(index=False, max_colwidth=46))
print('\nwritten:', os.path.basename(p))
