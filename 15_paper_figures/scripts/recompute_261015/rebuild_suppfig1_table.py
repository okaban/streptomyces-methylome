#!/usr/bin/env python3
"""補足図 1 の 2x2 表を正典部位集合で作り直す。

元実装は 87_site_coloc_stats_C4/site_coloc_stats_C4.py。設計はそのまま写す:
  単位      同一鎖で A が C の 3 bp または 4 bp 上流にある (A, C) 位置対。
            ゲノム全体で 1,696,280 対（配列だけで決まるので規則に依らない）。
  モチーフ内 その C が AAGCCCG の offset 4 である対（2 x 1,334 = 2,668 対）。
  共修飾    C が 4mC 部位 かつ A が 6mA 部位（同一鎖）。
  検定      Fisher 正確確率（片側 greater）、Haldane 補正オッズ比、
            超幾何ラベル置換帰無 10,000 回・seed 42。
  scope     pooled = 部位を T1-T3 で重複排除した和集合 / T1 = T1 の部位のみ。

変更点は部位の定義だけ。元は 07_motif_analysis/methylation_site_sequences.csv の
**撤回済み規則**（per-replicate depth >= 10, freq >= 50%）。本版は正典パイルアップから
3 ライブラリ合算・真深度 (Nvalid+Nfail) >= 30・合算 Nmod/Nvalid >= 0.5 で取り直す。

注意: 「pooled」は時点ごとの件数の和ではない。部位の和集合を取ってから対を数える。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/rebuild_suppfig1_table.py
"""
import pathlib
import re
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, 'supp_261014')
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, genome  # noqa: E402

PILEUP = pathlib.Path('/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015')
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / 'site_coloc_2x2_canonical_261015.tsv'
CACHE = HERE / 'canonical_sites_genomewide_261015.tsv'
CODE = {'a': '6mA', '21839': '4mC'}
# 参照塩基ガード（2026-10-15）。modkit はリードが参照とミスマッチしていても修飾コールを
# 参照座標に計上するため、参照が C の位置に 1 本だけ A を持つリードがあると
# Nvalid=1/Nmod=1 で level 100% になる。真深度を Nvalid+Nfail としていたので
# Nfail 側だけで閾値を越え、T1 6mA 7,297 件のうち 5,350 件が非アデニン位置だった。
REF_OK = {'a': ('A', 'T'), '21839': ('C', 'G')}
_G = None


def ref_base_matches(pos, strand, code):
    global _G
    if _G is None:
        _G = genome()
    plus, minus = REF_OK[code]
    return _G[pos] == (plus if strand == '+' else minus)


def build_sites():
    """時点ごとの (pos, strand, mod) 正典部位集合。キャッシュがあれば読む。"""
    if CACHE.exists():
        d = pd.read_csv(CACHE, sep='\t')
        return {tp: {m: set(zip(g.pos, g.strand)) for m, g in sub.groupby('mod')}
                for tp, sub in d.groupby('timepoint')}
    recs = []
    for tp, libs in LIBS.items():
        acc = {}
        for lib in libs:
            p = subprocess.run(
                ['awk', '-F', '\t',
                 '$4=="a"||$4=="21839"{print $2"\t"$6"\t"$4"\t"$10"\t"$12"\t"$16}',
                 str(PILEUP / f'{lib}_pileup.bed')],
                capture_output=True, text=True, check=True)
            for line in p.stdout.splitlines():
                pos, strand, code, nv, nm, nf = line.split('\t')
                k = (int(pos), strand, code)
                cur = acc.get(k)
                if cur is None:
                    acc[k] = [int(nv), int(nm), int(nf)]
                else:
                    cur[0] += int(nv); cur[1] += int(nm); cur[2] += int(nf)
        for (pos, strand, code), (nv, nm, nf) in acc.items():
            if not ref_base_matches(pos, strand, code):
                continue
            if (nv + nf) >= DEPTH_MIN and nv > 0 and nm / nv >= FRAC_MIN:
                recs.append((tp, pos, strand, CODE[code]))
    pd.DataFrame(recs, columns=['timepoint', 'pos', 'strand', 'mod']).to_csv(
        CACHE, sep='\t', index=False)
    print(f'部位集合をキャッシュ: {CACHE.name}（{len(recs):,} 行）')
    return build_sites()


sites = build_sites()
seq = genome()
G = len(seq)

# AAGCCCG の C4 位置（0 起点、鎖つき）
c4 = {('+', m.start() + 4) for m in re.finditer('AAGCCCG', seq)} | \
     {('-', m.start() + 2) for m in re.finditer('CGGGCTT', seq)}
assert len(c4) == 1334, len(c4)

arr = np.frombuffer(seq.encode(), dtype='S1')
pC, pA, pG, pT = (arr == b'C'), (arr == b'A'), (arr == b'G'), (arr == b'T')
rows = []
for sp in (3, 4):
    i = np.where(pC[sp:] & pA[:-sp])[0] + sp
    rows.append(pd.DataFrame(dict(strand='+', c_pos=i, a_pos=i - sp, spacing=sp)))
    i = np.where(pG[:-sp] & pT[sp:])[0]
    rows.append(pd.DataFrame(dict(strand='-', c_pos=i, a_pos=i + sp, spacing=sp)))
cand = pd.concat(rows, ignore_index=True)
cand['in_motif'] = [(s, c) in c4 for s, c in zip(cand.strand, cand.c_pos)]
assert cand.in_motif.sum() == 2 * 1334, cand.in_motif.sum()
print(f'候補対 {len(cand):,}（うちモチーフ内 {int(cand.in_motif.sum()):,}）')

scopes = {'pooled': {m: set().union(*(sites[tp][m] for tp in sites)) for m in ('4mC', '6mA')},
          'T1': sites['T1']}
out = []
for scope, S in scopes.items():
    m4, m6 = S['4mC'], S['6mA']
    cand['comod'] = [((c, s) in m4) and ((a, s) in m6)
                     for c, a, s in zip(cand.c_pos, cand.a_pos, cand.strand)]
    a = int((cand.in_motif & cand.comod).sum()); b = int((cand.in_motif & ~cand.comod).sum())
    c = int((~cand.in_motif & cand.comod).sum()); d = int((~cand.in_motif & ~cand.comod).sum())
    by_sp = cand[cand.in_motif & cand.comod].spacing.value_counts().to_dict()
    OR_raw, p_f = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
    aa, bb, cc, dd = (np.array([a, b, c, d]) + (0.5 if min(a, b, c, d) == 0 else 0)).astype(float)
    OR = aa * dd / (bb * cc); se = np.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    n_tot, n_in, n_comod = a + b + c + d, a + b, a + c
    exp_in = n_in * n_comod / n_tot
    rng = np.random.default_rng(42); N = 10_000
    ap = rng.hypergeometric(n_comod, n_tot - n_comod, n_in, N)
    n_inst = cand[cand.in_motif & cand.comod].drop_duplicates(['strand', 'c_pos']).shape[0]
    out.append(dict(scope=scope, n_instances_with_comod_pair=n_inst,
                    in_comod=a, in_notcomod=b, out_comod=c, out_notcomod=d,
                    n_candidate_pairs=n_tot, n_in_motif_pairs=n_in,
                    in_comod_spacing4_A0C4=by_sp.get(4, 0), in_comod_spacing3_A1C4=by_sp.get(3, 0),
                    OR=OR, OR_CI95_lower=float(np.exp(np.log(OR) - 1.96 * se)),
                    OR_CI95_upper=float(np.exp(np.log(OR) + 1.96 * se)),
                    haldane_corrected=(min(a, b, c, d) == 0), fisher_p_greater=p_f,
                    expected_in_comod_H0=exp_in, fold_over_expectation=a / exp_in,
                    perm_p=(f'< {1/(N+1):.1e}' if (ap >= a).sum() == 0
                            else f'{((ap >= a).sum() + 1) / (N + 1):.4f}'),
                    perm_null_max_in_comod=int(ap.max()),
                    rate_in=a / n_in, rate_out=c / (c + d),
                    rate_ratio=(a / n_in) / (c / (c + d)) if c else np.inf))
res = pd.DataFrame(out)
res.to_csv(OUT, sep='\t', index=False, float_format='%.6g')
old = {'pooled': (406, 0.638579, 5, 345), 'T1': (214, 0.336591, 4, 194)}
print(f"\n{'scope':<8}{'in_comod':>10}{'(旧)':>7}{'期待':>10}{'帰無最大':>9}{'外側':>6}{'インスタンス':>12}")
for r in out:
    o = old[r['scope']]
    print(f"{r['scope']:<8}{r['in_comod']:>10}{o[0]:>7}{r['expected_in_comod_H0']:>10.4f}"
          f"{r['perm_null_max_in_comod']:>9}{r['out_comod']:>6}"
          f"{r['n_instances_with_comod_pair']:>8}（旧 {o[3]}）")
print('->', OUT)
