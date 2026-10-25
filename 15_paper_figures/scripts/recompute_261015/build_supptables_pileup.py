#!/usr/bin/env python3
"""補足表 13a・13b・18 を正典パイルアップから一度に作る。

3 表とも GCCGGC / AAGCCCG のパイルアップ抽出を入力に持つので、読み込みを共有する。

  13a  GCCGGC 二重鎖のうち両鎖メチル化された数を、各鎖が TGGCCGGC 認識文脈を
       持つかで層別（認識文脈 0 / 1 / 2 個）。
  13b  同じ層別を、信頼度フィルタ版（正典）と all-read 版（--no-filtering）で比較。
  18   AAGCCCG の 6mA（位置 1）と 4mC（位置 5）の占有率を、モチーフ直上流 1 塩基で層別。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/build_supptables_pileup.py
"""
import csv
import pathlib
import sys

import pandas as pd

sys.path.insert(0, 'supp_261014')
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, find_all, genome  # noqa: E402
from _suppstyle import DATA_DIR  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
NOFILT = pathlib.Path('/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015')
g = genome()


def pooled_from_dir(sub, tp):
    acc = {}
    for lib in LIBS[tp]:
        d = pd.read_csv(DATA_DIR / sub / f'{lib}.tsv', sep='\t', header=None,
                        names=['pos', 'strand', 'Nvalid', 'Nmod', 'Nfail'])
        for pos, st, nv, nm, nf in d.itertuples(index=False):
            k = (pos, st); c = acc.get(k)
            if c is None: acc[k] = [nv, nm, nf]
            else: c[0] += nv; c[1] += nm; c[2] += nf
    return acc


def pooled_nofilt(tp, motif_tag):
    acc = {}
    for lib in LIBS[tp]:
        with open(NOFILT / f'{lib}_nofilt.bed') as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if not f[3].startswith(motif_tag): continue
                k = (int(f[1]), f[5]); nv, nm, nf = int(f[9]), int(f[11]), int(f[15])
                c = acc.get(k)
                if c is None: acc[k] = [nv, nm, nf]
                else: c[0] += nv; c[1] += nm; c[2] += nf
    return acc


def called(acc, key):
    v = acc.get(key)
    return v is not None and (v[0] + v[2]) >= DEPTH_MIN and v[0] > 0 and v[1] / v[0] >= FRAC_MIN


def measurable(acc, key):
    v = acc.get(key)
    return v is not None and (v[0] + v[2]) >= DEPTH_MIN and v[0] > 0


# ---- GCCGGC 二重鎖と認識文脈 ------------------------------------------------
duplex = []
for i in find_all(g, 'GCCGGC'):
    n_ctx = (1 if g[max(0, i - 2):i + 6] == 'TGGCCGGC' else 0) + \
            (1 if g[i:i + 8] == 'GCCGGCCA' else 0)
    duplex.append(((i + 2, '+'), (i + 3, '-'), n_ctx))
print(f'GCCGGC 二重鎖 {len(duplex):,}  認識文脈 0/1/2 = '
      + '/'.join(str(sum(1 for d in duplex if d[2] == k)) for k in (0, 1, 2)))

rows13a, rows13b = [], []
for tp in ['T1', 'T2', 'T3']:
    filt = pooled_from_dir('pileup_GCCGGC', tp)
    allr = pooled_nofilt(tp, '21839,GCCGGC')
    for label, acc in (('confidence-filtered', filt), ('all reads', allr)):
        for k in (0, 1, 2):
            sel = [d for d in duplex if d[2] == k]
            meas = [d for d in sel if measurable(acc, d[0]) and measurable(acc, d[1])]
            both = [d for d in meas if called(acc, d[0]) and called(acc, d[1])]
            one = [d for d in meas if (called(acc, d[0]) ^ called(acc, d[1]))]
            row = dict(timepoint=tp, read_set=label, recognition_contexts_on_duplex=k,
                       duplexes=len(sel), duplexes_measurable_both_strands=len(meas),
                       methylated_both_strands=len(both), methylated_one_strand=len(one),
                       both_strand_pct=round(100 * len(both) / len(meas), 2) if meas else None)
            rows13b.append(row)
            if label == 'confidence-filtered':
                rows13a.append(row)

# ---- AAGCCCG を 5' 隣接塩基で層別 -------------------------------------------
inst = []
for i in find_all(g, 'AAGCCCG'):
    inst.append((g[i - 1] if i > 0 else 'N', (i, '+'), (i + 4, '+')))
for j in find_all(g, 'CGGGCTT'):
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    nxt = g[j + 7] if j + 7 < len(g) else 'N'
    inst.append((comp.get(nxt, 'N'), (j + 6, '-'), (j + 2, '-')))
print(f'AAGCCCG インスタンス {len(inst):,}')

rows18 = []
for tp in ['T1', 'T2', 'T3']:
    a6 = pooled_from_dir('pileup_AAGCCCG_6mA', tp)
    c4 = pooled_from_dir('pileup_AAGCCCG_4mC', tp)
    for base in ['A', 'C', 'G', 'T']:
        sel = [x for x in inst if x[0] == base]
        if not sel: continue
        n6 = sum(1 for x in sel if called(a6, x[1]))
        n4 = sum(1 for x in sel if called(c4, x[2]))
        rows18.append(dict(timepoint=tp, five_prime_base=base, motif_instances=len(sel),
                           sites_6mA_position1=n6,
                           occupancy_6mA_pct=round(100 * n6 / len(sel), 2),
                           sites_4mC_position5=n4,
                           occupancy_4mC_pct=round(100 * n4 / len(sel), 2)))

for name, rows in [('SuppTable13a_both_strand_duplex_261015.tsv',
                    [r for r in rows13a]),
                   ('SuppTable13b_both_strand_by_readset_261015.tsv', rows13b),
                   ('SuppTable18_AAGCCCG_by_5prime_base_261015.tsv', rows18)]:
    out = HERE / name
    with out.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'-> {name}  {len(rows)} 行')

print('\n13a（フィルタ版、認識文脈 2 個の二重鎖）')
for r in rows13a:
    if r['recognition_contexts_on_duplex'] == 2:
        print(f"  {r['timepoint']}: 測定可能 {r['duplexes_measurable_both_strands']:>5}  "
              f"両鎖 {r['methylated_both_strands']:>5} ({r['both_strand_pct']}%)  "
              f"片鎖 {r['methylated_one_strand']:>5}")
print('\n18（5′ 塩基別 6mA 占有率）')
for r in rows18:
    if r['timepoint'] == 'T1':
        print(f"  5′{r['five_prime_base']}: n={r['motif_instances']:>4}  "
              f"6mA {r['occupancy_6mA_pct']:>6}%  4mC {r['occupancy_4mC_pct']:>6}%")
