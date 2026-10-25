#!/usr/bin/env python3
"""補足表 16 — 未帰属 6mA 画分の位置分布（正典部位集合）。

旧記述（本文凡例）は T1 6mA 1,934 件中 1,386 件（71.66%）が未帰属としていたが、
その部位集合は 2026-02 の旧統合センサス由来。正典規則で取り直すと 6mA の総数自体が
大きく変わる。退役させた補足図 5 が持っていた「全配列文脈での部位総数」も
本表に統合する（4mC 列）。これにより 4mC と 6mA の非対称——4mC はほぼ全てが
2 つの認識配列で説明できるのに、6mA は大半が未帰属——が 1 枚で読める。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/build_supptable16.py
"""
import csv
import pathlib
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, 'supp_261014')
from _gccggc import genome  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CACHE = HERE / 'canonical_sites_genomewide_261015.tsv'
TSS_CSV = ('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
           '18_tss_analyses/jeong2016_all_tss.csv')
OUT = HERE / 'SuppTable16_unassigned_6mA_261015.tsv'
TSS_HALF = 300

g = genome()

# 認識配列の鎖指定位置
aag_A, aag_C, tgg_C = set(), set(), set()
for m in re.finditer('AAGCCCG', g):
    i = m.start(); aag_A |= {(i, '+'), (i + 1, '+')}; aag_C.add((i + 4, '+'))
for m in re.finditer('CGGGCTT', g):
    j = m.start(); aag_A |= {(j + 6, '-'), (j + 5, '-')}; aag_C.add((j + 2, '-'))
for m in re.finditer('TGGCCGGC', g):
    tgg_C.add((m.start() + 4, '+'))
for m in re.finditer('GCCGGCCA', g):
    tgg_C.add((m.start() + 3, '-'))

t = pd.read_csv(TSS_CSV)
sel = t[t.category.astype(str).str.upper().isin(['P', 'S']) & t.gene_id.notna()
        & (t.gene_id.astype(str).str.strip().str.len() > 0)]
tss = np.sort(sel.tss_position.astype(int).values)


def near_tss(positions):
    p = np.asarray(sorted(positions))
    if len(p) == 0:
        return 0
    i = np.searchsorted(tss, p)
    lo = np.where(i > 0, p - tss[np.clip(i - 1, 0, None)], 10 ** 9)
    hi = np.where(i < len(tss), tss[np.clip(i, None, len(tss) - 1)] - p, 10 ** 9)
    return int((np.minimum(lo, hi) <= TSS_HALF).sum())


d = pd.read_csv(CACHE, sep='\t')
rows = []
for tp in ['T1', 'T2', 'T3']:
    sub = d[d.timepoint == tp]
    s6 = set(zip(sub[sub['mod'] == '6mA'].pos, sub[sub['mod'] == '6mA'].strand))
    s4 = set(zip(sub[sub['mod'] == '4mC'].pos, sub[sub['mod'] == '4mC'].strand))
    a_assigned = s6 & aag_A
    a_unassigned = s6 - aag_A
    c_assigned = (s4 & aag_C) | (s4 & tgg_C)
    c_unassigned = s4 - aag_C - tgg_C
    rows.append(dict(
        timepoint=tp,
        sites_4mC=len(s4), sites_6mA=len(s6),
        mC_in_recognition_sequence=len(c_assigned),
        mC_unassigned=len(c_unassigned),
        mC_unassigned_pct=round(100 * len(c_unassigned) / len(s4), 2),
        mA_in_AAGCCCG=len(a_assigned),
        mA_unassigned=len(a_unassigned),
        mA_unassigned_pct=round(100 * len(a_unassigned) / len(s6), 2),
        mA_unassigned_near_TSS=near_tss([p for p, _ in a_unassigned]),
        mA_unassigned_near_TSS_pct=round(
            100 * near_tss([p for p, _ in a_unassigned]) / len(a_unassigned), 2)))

with OUT.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)

print(f"{'tp':<4}{'4mC':>7}{'うち認識配列':>13}{'未帰属%':>9}"
      f"{'6mA':>8}{'うちAAGCCCG':>12}{'未帰属%':>9}{'未帰属のTSS近傍%':>17}")
for r in rows:
    print(f"{r['timepoint']:<4}{r['sites_4mC']:>7}{r['mC_in_recognition_sequence']:>13}"
          f"{r['mC_unassigned_pct']:>9}{r['sites_6mA']:>8}{r['mA_in_AAGCCCG']:>12}"
          f"{r['mA_unassigned_pct']:>9}{r['mA_unassigned_near_TSS_pct']:>17}")
print('\n旧記述: T1 6mA 1,934 中 未帰属 1,386（71.66%）')
print('->', OUT)
