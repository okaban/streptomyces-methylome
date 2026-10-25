# -*- coding: utf-8 -*-
"""region_composition_by_mark_261015.csv を再生成し、第 4 のマーク
`AAGCCCG 4mC + 6mA`（同一インスタンスに両マーク）の行を足す。

2026-10-24。著者コメント #7（補足図 4 を 2 列 2 行の 4 パネルにし、AAGCCCG の
4mC+6mA を加える）に対応するために書いた。元の 9 行を作った実装はリポジトリに
残っていなかったので、README.md の規則どおりに書き直し、**既存 9 行を小数第 2 位
まで完全再現することを assert** している。再現できない変更は 4 行目も信用できない。

区分規則（README.md のとおり、優先順）
  1 promoter      遺伝子体の外、かつ TSS の上流 1-500 bp
  2 5'UTR         TSS の 0 〜 +100 bp
  3 CDS-internal  CDS の内側
  4 intergenic    それ以外

再現に必要だった 2 点（元実装の挙動であり、ここでは意図的に踏襲する）:
  * 遺伝子体マスクは GFF の CDS 8,199 件すべてを座標として union している。
    うち 414 件はプラスミド NC_003903.1 / NC_003904.1 の行で、染色体座標に
    重ねるのは誤り。影響は promoter で最大 0.6 ポイント。
  * パイルアップの pos は 0-based、GFF は 1-based だが、元実装は +1 せずに
    突き合わせている（1 bp のずれ）。
  どちらも図の主張（どのマークも 9 割前後が CDS 内部）を変えない。4 パネルの
  内部整合を優先してそのまま使い、訂正は別途 EPI_ANALYST に回す。

デュアル集合の定義: C5（4mC）が呼ばれ、かつ同一インスタンスの A1 または A2 で
6mA が呼ばれたインスタンス。A1/A2 は鎖の向きに沿って C5 の 4 bp / 3 bp 5' 側
（+ 鎖なら pos-4 / pos-3、- 鎖なら pos+4 / pos+3）。この幾何は参照配列で検証済み。
領域の割り当ては C5 の座標で行う（他の AAGCCCG パネルと同じ基準）。

Run: python3 build_region_composition_261024.py
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent                      # scripts/supp_261014/data
RNA = pathlib.Path('/Users/okaban/bioinfo/rna-seq')
TSS_CSV = RNA / '11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv'
GFF = RNA / ('11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/'
             'data/GCF_000203835.1_ASM20383v1_genomic.gff')
OUT = HERE / 'region_composition_by_mark_261024.csv'

LIBS = {'T1': ['1-1', '1-2', '1-3'], 'T2': ['2-1', '2-3', '2-4'],
        'T3': ['3-2', '3-3', '3-4']}
COLS = ['pos', 'strand', 'Nvalid', 'Nmod', 'Nfail']
ORDER = ['promoter', "5'UTR", 'CDS-internal', 'intergenic']


def load_universe(sub):
    """(pos, strand) の全行と、時点ごとの呼び出し済み部位集合を返す。"""
    universe, called = set(), {}
    for tp, libs in LIBS.items():
        fr = [pd.read_csv(DATA / sub / f'{l}.tsv', sep='\t', header=None, names=COLS)
              for l in libs]
        df = pd.concat(fr)
        universe |= set(zip(df.pos, df.strand))
        df['true'] = df.Nvalid + df.Nfail
        g = df.groupby(['pos', 'strand']).agg(Nvalid=('Nvalid', 'sum'),
                                              Nmod=('Nmod', 'sum'),
                                              true=('true', 'sum'))
        g = g[(g.true >= 30) & (g.Nvalid > 0) & (g.Nmod / g.Nvalid >= 0.5)]
        called[tp] = set(g.index)
    return universe, called


def build_masks():
    cds = []
    for ln in open(GFF):
        if ln.startswith('#'):
            continue
        f = ln.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] == 'CDS':
            cds.append((int(f[3]), int(f[4])))
    glen = max(e for _, e in cds) + 2000
    body = np.zeros(glen, dtype=bool)
    for s, e in cds:
        body[s:min(e + 1, glen)] = True

    tss = pd.read_csv(TSS_CSV)
    sel = tss[tss.category.isin(['P', 'S']) & tss.gene_id.notna()
              & (tss.gene_id.astype(str).str.strip() != '')]
    assert len(sel) == 3030, f'TSS 件数が 3,030 でなく {len(sel)}'
    prom = np.zeros(glen, dtype=bool)
    utr = np.zeros(glen, dtype=bool)
    for p, st in zip(sel.tss_position.values, sel.strand.values):
        if st == '+':
            prom[max(1, p - 500):p] = True
            utr[p:min(glen - 1, p + 100) + 1] = True
        else:
            prom[p + 1:min(glen - 1, p + 500) + 1] = True
            utr[max(1, p - 100):p + 1] = True
    return body, prom, utr


def assign(positions, body, prom, utr):
    p = np.asarray(sorted(positions))
    inb = body[p]
    r = np.full(len(p), 'intergenic', dtype=object)
    r[inb] = 'CDS-internal'
    r[utr[p] & ~(prom[p] & ~inb)] = "5'UTR"
    r[prom[p] & ~inb] = 'promoter'
    return pd.Series(r)


def adenines(pos, strand):
    sg = 1 if strand == '+' else -1
    return [(pos - 4 * sg, strand), (pos - 3 * sg, strand)]


def main():
    body, prom, utr = build_masks()
    uni4C, cal4C = load_universe('pileup_GCCGGC')
    uniA4, calA4 = load_universe('pileup_AAGCCCG_4mC')
    uniA6, calA6 = load_universe('pileup_AAGCCCG_6mA')

    # 幾何の検証: すべての C5 インスタンスが 6mA 側に A1/A2 の行を持つ
    missing = [k for k in uniA4 if not all(a in uniA6 for a in adenines(*k))]
    assert not missing, f'A1/A2 の行が無い C5 インスタンス {len(missing)} 件'

    rows = []
    for key, label in [('pileup_GCCGGC', 'GCCGGC hexamer 4mC'),
                       ('pileup_AAGCCCG_4mC', 'AAGCCCG 4mC'),
                       ('pileup_AAGCCCG_6mA', 'AAGCCCG 6mA')]:
        called = {'pileup_GCCGGC': cal4C, 'pileup_AAGCCCG_4mC': calA4,
                  'pileup_AAGCCCG_6mA': calA6}[key]
        for tp in LIBS:
            pos = [p for p, _ in called[tp]]
            vc = assign(pos, body, prom, utr).value_counts(normalize=True) * 100
            rows.append(dict(mark=label, timepoint=tp, n_sites=len(pos),
                             **{k: round(vc.get(k, 0.0), 2) for k in ORDER}))

    # 第 4 のマーク: 同一インスタンスに 4mC(C5) と 6mA(A1 or A2) の両方
    for tp in LIBS:
        dual = [k for k in uniA4
                if k in calA4[tp] and any(a in calA6[tp] for a in adenines(*k))]
        vc = assign([p for p, _ in dual], body, prom, utr).value_counts(normalize=True) * 100
        rows.append(dict(mark='AAGCCCG 4mC + 6mA', timepoint=tp, n_sites=len(dual),
                         **{k: round(vc.get(k, 0.0), 2) for k in ORDER}))

    new = pd.DataFrame(rows)

    # 既存 9 行の完全再現を要求する
    old = pd.read_csv(HERE / 'region_composition_by_mark_261015.csv')
    chk = new.merge(old, on=['mark', 'timepoint'], suffixes=('_n', '_p'))
    assert len(chk) == 9, f'既存行と突き合わせられたのが {len(chk)} 行'
    assert (chk.n_sites_n == chk.n_sites_p).all(), '部位数が既存表と一致しない'
    diff = np.abs(chk[[k + '_n' for k in ORDER]].to_numpy()
                  - chk[[k + '_p' for k in ORDER]].to_numpy()).max()
    assert diff == 0.0, f'既存 9 行を再現できていない（最大 {diff} ポイント）'

    new.to_csv(OUT, index=False)
    print(f'再現 OK（既存 9 行の最大差 {diff} ポイント）')
    print(new.to_string(index=False))


if __name__ == '__main__':
    main()
