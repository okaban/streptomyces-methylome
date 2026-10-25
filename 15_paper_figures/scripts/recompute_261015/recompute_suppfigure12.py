#!/usr/bin/env python3
"""補足図 12 の凡例に書かれた数値を再計算する（図は描き直さない）。

目的: 補足図 12 は出荷 PNG のみで生成元スクリプトが存在しない。数値が動いたとき
図は描き直せないが、凡例のどの数値も追跡可能にしておく。本スクリプトは TSV を
1 本出すだけで、PNG には一切触れない。

凡例が主張している値:
    GCCGGC 六塩基    観測 2,539   期待 2,883.8   比 0.8804 (95% CI 0.8478-0.9160)  非整合 0.7525
    AAGCCCG 両鎖     観測   305   期待   322.0   比 0.9471 (0.8424-1.0601)         非整合 1.0898
    GGC または GCC   観測 142,109 期待 146,618.7 比 0.9692 (0.9646-0.9739)         非整合 0.8675
    TSS +/-300 bp の GC 含量 69.54 % 対 ゲノム全体 72.12 %

方法: 実験的 primary / secondary TSS を中心とする 601 bp 窓。ゲノムを同じ幅の
重なりのない窓に分けて GC 含量の五分位に層別し、各 TSS 窓を自身の GC が属する
五分位に割り当て、その五分位の 1 窓あたり平均出現数から期待数を出す（GC 整合）。
非整合の比はゲノム全体の平均密度を使ったもの。区間は観測数の Byar ポアソン近似。

実行: python3 recompute_261015/recompute_suppfigure12.py
"""
import pathlib
import re

import numpy as np
import pandas as pd

BASE = pathlib.Path('/Users/okaban/bioinfo/rna-seq')
GENOME = BASE / '11_epigenome_integration/data/NC_003888.3.fna'
TSS = BASE / '11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv'
OUT = pathlib.Path(__file__).resolve().parent / 'suppfigure12_recomputed.tsv'
HALF, NQ = 300, 5
WIN = 2 * HALF + 1

CLAIMED = {
    'GCCGGC hexamer': (2539, 2883.8, 0.8804, 0.8478, 0.9160, 0.7525),
    'AAGCCCG both strands': (305, 322.0, 0.9471, 0.8424, 1.0601, 1.0898),
    'GGC or GCC trinucleotide': (142109, 146618.7, 0.9692, 0.9646, 0.9739, 0.8675),
}
MOTIFS = {
    'GCCGGC hexamer': ['GCCGGC'],
    'AAGCCCG both strands': ['AAGCCCG', 'CGGGCTT'],
    'GGC or GCC trinucleotide': ['GGC', 'GCC'],
}


def genome():
    with open(GENOME) as fh:
        return ''.join(l.strip() for l in fh if not l.startswith('>')).upper()


def tss_positions():
    d = pd.read_csv(TSS)
    cat = d['category'].astype(str).str.strip().str.upper()
    gid = d['gene_id'].astype(str).str.strip()
    sel = d[cat.isin(['P', 'S']) & gid.ne('') & gid.ne('nan')]
    return sel['tss_position'].astype(int).to_numpy()


def byar_ci(obs, exp):
    lo = obs * (1 - 1 / (9 * obs) - 1.959964 / (3 * np.sqrt(obs))) ** 3
    hi = (obs + 1) * (1 - 1 / (9 * (obs + 1)) + 1.959964 / (3 * np.sqrt(obs + 1))) ** 3
    return lo / exp, hi / exp


def main():
    g = genome()
    glen = len(g)
    b = np.frombuffer(g.encode(), dtype=np.uint8)
    cum = np.concatenate([[0], np.cumsum(((b == ord('G')) | (b == ord('C'))).astype(np.int64))])

    tss = tss_positions()
    tss = tss[(tss - HALF >= 0) & (tss + HALF + 1 <= glen)]
    starts = tss - HALF

    bg = np.arange(0, glen - WIN, WIN)
    gcf = lambda s: (cum[s + WIN] - cum[s]) / WIN            # noqa: E731
    edges = np.quantile(gcf(bg), np.linspace(0, 1, NQ + 1))[1:-1]
    bg_q, tss_q = np.digitize(gcf(bg), edges), np.digitize(gcf(starts), edges)

    rows = []
    for name, pats in MOTIFS.items():
        hits = np.sort(np.array([m.start() for p in pats
                                 for m in re.finditer(f'(?={p})', g)]))
        cnt = lambda s: (np.searchsorted(hits, s + WIN)      # noqa: E731
                         - np.searchsorted(hits, s))
        obs = int(cnt(starts).sum())
        bgc = cnt(bg)
        exp_m = float(np.array([bgc[bg_q == q].mean() for q in range(NQ)])[tss_q].sum())
        exp_p = float(bgc.mean() * len(starts))
        lo, hi = byar_ci(obs, exp_m)
        c = CLAIMED[name]
        rows.append(dict(quantity=name, observed=obs, claimed_observed=c[0],
                         expected_gc_matched=round(exp_m, 1), claimed_expected=c[1],
                         ratio=round(obs / exp_m, 4), claimed_ratio=c[2],
                         ci95_low=round(lo, 4), ci95_high=round(hi, 4),
                         claimed_ci=f'{c[3]}-{c[4]}',
                         ratio_unmatched=round(obs / exp_p, 4), claimed_unmatched=c[5]))
    rows.append(dict(quantity='GC content, TSS +/-300 bp (%)',
                     observed=round(100 * gcf(starts).mean(), 2), claimed_observed=69.54,
                     expected_gc_matched=round(100 * float(cum[-1]) / glen, 2),
                     claimed_expected=72.12))
    rows.append(dict(quantity='TSS windows used', observed=len(starts),
                     claimed_observed=3030))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep='\t', index=False)
    print(df.to_string(index=False))
    print(f'\n-> {OUT}')


if __name__ == '__main__':
    main()
