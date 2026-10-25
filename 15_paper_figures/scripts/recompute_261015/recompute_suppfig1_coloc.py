#!/usr/bin/env python3
"""補足図 1 — 共修飾ペアのモチーフ特異性を正典部位集合で再計算する。

旧版 (87_site_coloc_stats_C4/site_coloc_stats_C4.py) は
「Co-modified = the C is a canonical 4mC site AND the A a canonical 6mA site
(depth >= 10 ...)」すなわち **撤回済みの per-replicate 深度 10 規則**で部位を定義していた。
正典規則は 3 ライブラリ合算・真深度 (Nvalid+Nfail) >= 30・合算 Nmod/Nvalid >= 0.5。

候補ペアの母集団（同一鎖で 3-4 bp 離れた A/C の組）は配列だけで決まるので規則に依らない。
変わるのは「そのうち両方が部位と呼ばれる数」。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/recompute_suppfig1_coloc.py
"""
import csv
import pathlib
import subprocess
import sys

import numpy as np

sys.path.insert(0, 'supp_261014')
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, find_all, genome  # noqa: E402

PILEUP = pathlib.Path('/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015')
OUT = pathlib.Path(__file__).resolve().parent / 'suppfig1_coloc_recompute_261015.tsv'
CODE = {'a': '6mA', '21839': '4mC'}
REF_OK = {'a': ('A', 'T'), '21839': ('C', 'G')}
_GEN = genome()
SPACING = (3, 4)


def called_sites(libs):
    """合算して正典規則を満たす (pos, strand) を mod 種別ごとに返す。"""
    acc = {}
    for lib in libs:
        # awk で該当 mod code 行だけを抜き、(pos,strand,code) -> [Nvalid, Nmod, Nfail] に積む
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
    out = {'6mA': set(), '4mC': set()}
    for (pos, strand, code), (nv, nm, nf) in acc.items():
        # 参照塩基ガード（2026-10-15）。ミスマッチしたリードの修飾コールが参照座標に
        # 計上されるため、参照が A/T でない位置の 6mA を落とす。rebuild_suppfig1_table.py と同規則。
        plus, minus = REF_OK[code]
        if _GEN[int(pos)] != (plus if strand == '+' else minus):
            continue
        if (nv + nf) >= DEPTH_MIN and nv > 0 and nm / nv >= FRAC_MIN:
            out[CODE[code]].add((pos, strand))
    return out


g = genome()
# AAGCCCG インスタンスの A0/A1 と C4（鎖指定）
motif_A, motif_C = set(), set()
for i in find_all(g, 'AAGCCCG'):
    motif_A |= {(i, '+'), (i + 1, '+')}
    motif_C.add((i + 4, '+'))
for j in find_all(g, 'CGGGCTT'):
    motif_A |= {(j + 6, '-'), (j + 5, '-')}
    motif_C.add((j + 2, '-'))

rows = []
for tp in ['T1', 'T2', 'T3']:
    s = called_sites(LIBS[tp])
    a_sites, c_sites = s['6mA'], s['4mC']
    in_motif = out_motif = 0
    for (cpos, strand) in c_sites:
        for d in SPACING:
            apos = cpos - d if strand == '+' else cpos + d
            if (apos, strand) in a_sites:
                if (cpos, strand) in motif_C and (apos, strand) in motif_A:
                    in_motif += 1
                else:
                    out_motif += 1
    rows.append(dict(timepoint=tp, sites_6mA=len(a_sites), sites_4mC=len(c_sites),
                     comod_pairs_total=in_motif + out_motif,
                     comod_in_AAGCCCG=in_motif, comod_outside_AAGCCCG=out_motif))
    print(f'{tp}: 6mA {len(a_sites):>6}  4mC {len(c_sites):>6}  '
          f'共修飾ペア {in_motif + out_motif:>4}（AAGCCCG 内 {in_motif}, 外 {out_motif}）')

tot_in = sum(r['comod_in_AAGCCCG'] for r in rows)
tot_out = sum(r['comod_outside_AAGCCCG'] for r in rows)
rows.append(dict(timepoint='pooled', sites_6mA='', sites_4mC='',
                 comod_pairs_total=tot_in + tot_out,
                 comod_in_AAGCCCG=tot_in, comod_outside_AAGCCCG=tot_out))
with OUT.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader(); w.writerows(rows)
print(f'\npooled: AAGCCCG 内 {tot_in} / 外 {tot_out}  （旧値 406 / 0）')
print('->', OUT)
