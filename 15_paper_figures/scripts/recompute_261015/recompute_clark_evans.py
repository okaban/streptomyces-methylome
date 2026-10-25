#!/usr/bin/env python3
"""Methods の Clark-Evans 最近接解析を正典部位集合で計算し直す。

旧記述は撤回済みの部位集合（per-replicate depth >= 10）で作られていた。
設計は Methods の記述どおりに写す:

  * 円環染色体上の Clark-Evans 最近接統計量。完全空間ランダム下の期待最近接距離は L/2n。
  * **インスタンス水準**で計算する（1 インスタンスが 1(6mA)-5(4mC) と 2(6mA)-5(4mC) の
    両方を出しうるが、その 2 点は 1 bp しか離れず、イベント水準だと構成上 1 bp の
    隣接点を置いてしまうため）。
  * モチーフ条件付き帰無: ゲノム上の AAGCCCG 1,334 インスタンスから同数を無作為抽出。
    2,000 回、seed 42。
  * 比較のためイベント水準の R も出す。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/recompute_clark_evans.py
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
OUT = HERE / 'clark_evans_canonical_261015.tsv'
NDRAW, SEED = 2000, 42

g = genome()
L = len(g)

# AAGCCCG インスタンス: 中心座標（C4 の位置）と、その A0/A1/C4 の鎖つき座標
inst = []
for m in re.finditer('AAGCCCG', g):
    i = m.start()
    inst.append(dict(centre=i + 4, strand='+', a=[(i, '+'), (i + 1, '+')], c=(i + 4, '+')))
for m in re.finditer('CGGGCTT', g):
    j = m.start()
    inst.append(dict(centre=j + 2, strand='-', a=[(j + 6, '-'), (j + 5, '-')], c=(j + 2, '-')))
assert len(inst) == 1334, len(inst)
all_centres = np.array(sorted(x['centre'] for x in inst))

d = pd.read_csv(CACHE, sep='\t')
sites = {tp: {m: set(zip(s.pos, s.strand)) for m, s in sub.groupby('mod')}
         for tp, sub in d.groupby('timepoint')}


def clark_evans(pos, n_univ_positions=None):
    """円環上の R = 観測平均最近接距離 / (L/2n)。"""
    p = np.sort(np.asarray(pos, dtype=float))
    n = len(p)
    if n < 2:
        return np.nan
    gaps = np.diff(np.concatenate([p, [p[0] + L]]))
    nn = np.minimum(gaps, np.roll(gaps, 1))
    return nn.mean() / (L / (2 * n))


rng = np.random.default_rng(SEED)
rows = []
for tp in ['T1', 'T2', 'T3']:
    m4, m6 = sites[tp]['4mC'], sites[tp]['6mA']
    ev, ins = [], []
    for x in inst:
        hits = [a for a in x['a'] if a in m6] if x['c'] in m4 else []
        if hits:
            ins.append(x['centre'])
            ev.extend(a[0] for a in hits)
    R_i, R_e = clark_evans(ins), clark_evans(ev)
    # モチーフ条件付き帰無: 1,334 インスタンスから同数を無作為抽出
    null = np.array([clark_evans(rng.choice(all_centres, size=len(ins), replace=False))
                     for _ in range(NDRAW)])
    p_two = 2 * min((null <= R_i).mean(), (null >= R_i).mean())
    rows.append(dict(timepoint=tp, n_instances=len(ins), n_events=len(ev),
                     R_instance=round(R_i, 4), R_event=round(R_e, 4),
                     null_mean=round(float(null.mean()), 4),
                     null_sd=round(float(null.std(ddof=1)), 4),
                     p_two_sided=round(min(p_two, 1.0), 4)))
    print(f'{tp}: インスタンス {len(ins):>4}  イベント {len(ev):>4}  '
          f'R(inst) {R_i:.3f}  R(event) {R_e:.3f}  帰無 {null.mean():.3f}±{null.std(ddof=1):.3f}  '
          f'p {min(p_two,1.0):.3f}')

with OUT.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print('\n旧記述: R(inst) 1.010 / 0.990 / 0.963, p 0.855 / 0.847 / 0.444; '
      'R(event) 0.813 / 0.761 / 0.749')
print('->', OUT)
