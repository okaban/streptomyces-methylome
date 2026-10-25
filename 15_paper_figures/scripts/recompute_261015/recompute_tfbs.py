#!/usr/bin/env python3
"""TFBS とメチル化の関係を正典データで再計算する（旧・補足図 11b の置き換え）。

旧版（analysis/49_TF_BS_methylation_protection, 2026-02-26）の問題:
  1. 使用サイト数が撤回済み（GCCGGC T1 = 1,289。正典は TGGCCGGC で 1,886）。
     同ディレクトリの _RETRACTED_VALUES.md に記載済み。
  2. 検定が「TFBS 区間あたり、メチル化サイトが 1 つ以上重なるか」の二項検定で、
     **配列可用性と修飾されやすさを分離していない**。モチーフ配列そのものが
     TFBS 上で偏在するため、fold 1.157（濃縮）はほぼ配列側の偏りだった。

本スクリプトは 3 段階で計算する:
  (S1) 旧検定の再現（実装検証）— 撤回サイト数 1,289 で公称期待値 322.35 を再現する。
  (S2) 配列側とメチル化側を別々に測り、配列可用性で割った修飾固有成分を出す。
  (S3) 配列の有無で条件づけた検定 — モチーフ位置を母集団に、TFBS 内外で
       「呼ばれる率」（Fisher）と「修飾レベル」（Mann-Whitney、局所 GC 五分位整合）を比べる。
  (S4) 交絡確認 — 予測 TFBS がプロモーターに集中していないかを測る。
       集中していれば TFBS 効果とプロモーター効果は分離不能になる。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/recompute_tfbs.py
"""
import csv
import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, 'supp_261014')
from _gccggc import (DEPTH_MIN, LIBS, called, find_all, genome,  # noqa: E402
                     motif_positions, pooled)
from _suppstyle import DATA_DIR  # noqa: E402

GFF = ("/Users/okaban/bioinfo/rna-seq/13_TF_binding-site/analysis/"
       "01_master_TF_list_260206_v1/intermediate/fimo_results/fimo.gff")
TSS_CSV = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
           "18_tss_analyses/jeong2016_all_tss.csv")
OUT = pathlib.Path(__file__).resolve().parent / "tfbs_recompute_261015.tsv"
CHROM, FLANK, TSS_HALF, GC_HALF = "NC_003888.3", 10, 300, 150

g = genome()
L = len(g)
rows = []


def emit(block, quantity, value, unit, note):
    rows.append(dict(block=block, quantity=quantity, value=value, unit=unit, note=note))


# ---- FIMO TFBS ----
raw = []
with open(GFF) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[0] == CHROM:
            raw.append((int(f[3]) - 1, int(f[4]) - 1))
tf = pd.DataFrame(raw, columns=['start', 'end'])
tf['width'] = tf.end - tf.start + 1
emit("input", "FIMO TFBS hits", len(tf), "sites", f"{CHROM}; median width {tf.width.median():.0f} bp")

iv = sorted((s - FLANK, e + FLANK) for s, e in raw)
merged = []
for s, e in iv:
    if merged and s <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], e)
    else:
        merged.append([s, e])
ms = np.array([m[0] for m in merged])
me = np.array([m[1] for m in merged])
emit("input", "merged TFBS intervals", len(merged), "intervals",
     f"+/-{FLANK} bp flank; {100 * sum(e - s + 1 for s, e in merged) / L:.1f}% of the genome")


def inside(pos):
    i = np.searchsorted(ms, pos, 'right') - 1
    return (i >= 0) & (pos <= np.where(i >= 0, me[np.clip(i, 0, None)], -1))


def expected(n_sites):
    """一様配置の帰無で、各 TFBS +/-10 bp にサイトが 1 つ以上入る確率の総和。"""
    return float((1.0 - (1.0 - n_sites / L) ** (tf.width.values + 2 * FLANK)).sum())


def overlap_test(sites):
    arr = np.array(sorted(sites))
    obs = int(((np.searchsorted(arr, tf.end.values + FLANK, 'right')
                - np.searchsorted(arr, tf.start.values - FLANK, 'left')) > 0).sum())
    exp = expected(len(arr))
    p = stats.binomtest(obs, len(tf), exp / len(tf), alternative='two-sided').pvalue
    return len(arr), obs, exp, obs / exp, p


# ---- サイト集合 ----
_, tgg_keys = motif_positions(g)
c_keys = ([(i + 4, '+') for i in find_all(g, "AAGCCCG")]
          + [(j + 2, '-') for j in find_all(g, "CGGGCTT")])
a_keys = ([(i, '+') for i in find_all(g, "AAGCCCG")]
          + [(i + 1, '+') for i in find_all(g, "AAGCCCG")]
          + [(j + 6, '-') for j in find_all(g, "CGGGCTT")]
          + [(j + 5, '-') for j in find_all(g, "CGGGCTT")])


def load(sub, keys):
    acc = {}
    for lib in LIBS['T1']:
        d = pd.read_csv(DATA_DIR / sub / f"{lib}.tsv", sep='\t', header=None,
                        names=['pos', 'strand', 'Nvalid', 'Nmod', 'Nfail'])
        for pos, st, nv, nm, nf in d.itertuples(index=False):
            k = (pos, st)
            cur = acc.get(k)
            if cur is None:
                acc[k] = [nv, nm, nf]
            else:
                cur[0] += nv
                cur[1] += nm
                cur[2] += nf
    return acc


SETS = [("TGGCCGGC 4mC", tgg_keys, pooled('T1')),
        ("AAGCCCG 4mC pos5", c_keys, load('pileup_AAGCCCG_4mC', c_keys)),
        ("AAGCCCG 6mA pos1+2", a_keys, load('pileup_AAGCCCG_6mA', a_keys))]

# ---- S1 実装検証 ----
emit("S1_validation", "expected overlap at n=1289 (retracted site count)",
     round(expected(1289), 2), "TFBS", "published value 322.35 — implementation reproduces it")

# ---- S2 配列側 vs メチル化側 ----
for name, keys, acc in SETS:
    seq = [k[0] for k in keys]
    met = [k[0] for k in called(acc, keys)]
    n_s, o_s, e_s, f_s, p_s = overlap_test(seq)
    n_m, o_m, e_m, f_m, p_m = overlap_test(met)
    emit("S2_sequence", f"{name}: all motif positions, fold", round(f_s, 3), "obs/exp",
         f"n={n_s}, obs={o_s}, exp={e_s:.2f}, p={p_s:.3g}")
    emit("S2_methylated", f"{name}: called sites, fold", round(f_m, 3), "obs/exp",
         f"n={n_m}, obs={o_m}, exp={e_m:.2f}, p={p_m:.3g}")
    emit("S2_specific", f"{name}: methylation-specific component", round(f_m / f_s, 3),
         "ratio of folds", "called-site fold divided by sequence fold")

# ---- S3 配列で条件づけた検定 ----
def local_gc(pos):
    out = np.empty(len(pos))
    for n, p in enumerate(pos):
        w = g[max(0, p - GC_HALF):p + GC_HALF + 1]
        out[n] = (w.count('G') + w.count('C')) / len(w)
    return out


for name, keys, acc in SETS:
    pos = np.array([k[0] for k in keys])
    ins = inside(pos)
    cs = set(called(acc, keys))
    met = np.array([k in cs for k in keys])
    a, b = int((ins & met).sum()), int((ins & ~met).sum())
    c, d = int((~ins & met).sum()), int((~ins & ~met).sum())
    orr, pf = stats.fisher_exact([[a, b], [c, d]])
    emit("S3_called_rate", f"{name}: called rate inside TFBS", round(100 * a / (a + b), 1), "%",
         f"outside {100 * c / (c + d):.1f}%; OR {orr:.3f}; Fisher p={pf:.3g}; "
         f"motif positions in/out {a + b}/{c + d}")

    rec = [(k[0], acc[k][1] / acc[k][0]) for k in keys
           if k in acc and acc[k][0] > 0 and (acc[k][0] + acc[k][2]) >= DEPTH_MIN]
    lp = np.array([r[0] for r in rec])
    lv = np.array([r[1] for r in rec])
    li = inside(lp)
    mw = stats.mannwhitneyu(lv[li], lv[~li]).pvalue
    q = pd.qcut(local_gc(lp), 5, labels=False, duplicates='drop')
    dif = [(lv[li & (q == k)].mean() - lv[~li & (q == k)].mean(), int((li & (q == k)).sum()))
           for k in np.unique(q) if (li & (q == k)).sum() > 0 and (~li & (q == k)).sum() > 0]
    wd = sum(x * n for x, n in dif) / sum(n for _, n in dif)
    emit("S3_level", f"{name}: level inside TFBS", round(float(lv[li].mean()), 4), "fraction",
         f"outside {lv[~li].mean():.4f}; difference {lv[li].mean() - lv[~li].mean():+.4f}; "
         f"Mann-Whitney p={mw:.3g}; measurable positions in/out {li.sum()}/{(~li).sum()}")
    emit("S3_level_gcmatched", f"{name}: level difference, GC-matched", round(float(wd), 4),
         "fraction", f"local GC +/-{GC_HALF} bp in five quantiles, n-weighted")

# ---- S4 プロモーターとの交絡 ----
t = pd.read_csv(TSS_CSV)
sel = t[t.category.astype(str).str.upper().isin(['P', 'S']) & t.gene_id.notna()
        & (t.gene_id.astype(str).str.strip().str.len() > 0)]
tp = np.sort(sel.tss_position.astype(int).values)
emit("S4_input", "TSS used", len(tp), "TSS", "Jeong 2016, category P or S with a non-empty gene_id")


def near_tss(pos):
    i = np.searchsorted(tp, pos)
    lo = np.where(i > 0, pos - tp[np.clip(i - 1, 0, None)], 10 ** 9)
    hi = np.where(i < len(tp), tp[np.clip(i, None, len(tp) - 1)] - pos, 10 ** 9)
    return np.minimum(lo, hi) <= TSS_HALF


sample = np.arange(0, L, 97)
si, sn = inside(sample), near_tss(sample)
emit("S4_collinearity", "TSS-proximal share inside TFBS", round(100 * sn[si].mean(), 1), "%",
     f"outside {100 * sn[~si].mean():.1f}%; enrichment {sn[si].mean() / sn[~si].mean():.2f}x; "
     f"genome sampled every 97 bp (n={len(sample)}); TSS window +/-{TSS_HALF} bp")
for name, keys, _ in SETS[:1] + SETS[2:]:
    p = np.array(sorted(k[0] for k in keys))
    i2, n2 = inside(p), near_tss(p)
    emit("S4_collinearity", f"{name.split(' ')[0]} motif positions: TSS-proximal share inside TFBS",
         round(100 * n2[i2].mean(), 1), "%",
         f"outside {100 * n2[~i2].mean():.1f}%; {i2.sum()}/{len(p)} positions inside TFBS")

with OUT.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['block', 'quantity', 'value', 'unit', 'note'],
                       delimiter='\t')
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(f"{r['block']:<20}{r['quantity']:<58}{str(r['value']):>10}  {r['note']}")
print(f"\n-> {OUT}")
