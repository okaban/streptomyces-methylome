#!/usr/bin/env python3
"""Supplementary Table 17（モチーフ別メチル化センサス）の生成器。

## なぜこれが要るのか

2026-10-23 時点で、この表は手 staging されており生成スクリプトが存在しなかった。
本文 L205・補足表 17 凡例 L657・補足図 4 の凡例・value_registry の 6 行は
すべてこの 1 ファイルの写しなので、それらが互いに一致していても独立した
裏づけにはならない。D024（T3 の TGGCCGGC 4mC 部位数の数え直し）で
正準ブロックが 1,877、出荷物が 1,866 と食い違ったとき、どちらが一次データの
値かを決める手段が無かった。

このスクリプトは一次データ（参照 FASTA ＋ 正準パイルアップ 9 本）だけから
表全体を作り直す。D024 の判定もこれで行った。

## 規則（LOCKED CANON の「再現可能規則」）

  - 時点ごとに 3 ライブラリを pos×strand で Nvalid_cov / Nmod / Nfail 加算
  - eligible: プール後の真の深度 (Nvalid_cov + Nfail) >= 30 かつ
    3 ライブラリの少なくとも 1 本がその位置に行を出している
  - methylated: eligible かつ プール Nmod / Nvalid >= 0.5

座標規約は内部対照で検証される（下の EXPECT を参照）。

## 列の定義

  genome_occurrences   モチーフの出現数（＋鎖出現 ＋ 相補配列の出現）。
                       GAACCGG のように 1 出現に修飾位置が 2 つある場合も
                       「出現数」であって位置数ではない。
  expected_by_composition
                       2 * L * Π(単塩基頻度)。IUPAC は該当塩基の頻度和。
                       L はゲノム長、係数 2 は両鎖。
  eligible_positions   上の eligible を満たした鎖指定位置の数。
  methylated_sites     上の methylated を満たした鎖指定位置の数。
  occupancy_pct        methylated_sites / eligible_positions * 100
  occupancy_pct_over_genome_positions
                       methylated_sites / （鎖指定位置の総数）* 100
                       ← 原稿が報告する占有率はこちら（canon: 分母 2,093）。

使い方:  python3 make_supptable17_motif_census.py [--out PATH] [--check]
  --check  既存の出荷ファイルと差分を取り、一致しない行を報告する（書き換えない）
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

REF = Path('/Users/okaban/bioinfo/methyl/260102_M145/data/ref.fa')
PILEUP = Path('/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015')
OUT = Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary'
           '/SuppTable17_motif_census_canonical_261019.tsv')

TIMEPOINTS = {'T1': ['1-1', '1-2', '1-3'],
              'T2': ['2-1', '2-3', '2-4'],
              'T3': ['3-2', '3-3', '3-4']}

DEPTH_MIN = 30      # プール後の真の深度の下限
LEVEL_CUT = 0.5     # プール Nmod/Nvalid の下限

MODCODE = {'4mC': '21839', '6mA': 'a'}   # modkit bedMethyl の 4 列目
# M を落とすと CCGKCA（相補 TGMCGG）の −鎖が丸ごと消える。2026-10-23 に実際に
# 起きた: 出現数が 26,585 ではなく 12,985（＋鎖のみ）になった。COMP が作りうる
# 文字はすべてここに載っていなければならない（下の表明で強制する）。
IUPAC = {'H': 'ACT', 'B': 'CGT', 'D': 'AGT', 'V': 'ACG', 'N': 'ACGT',
         'K': 'GT', 'M': 'AC', 'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT'}
COMP = str.maketrans('ACGTHKMNRYSWBDV', 'TGCADMKNYRSWVHB')
assert set(COMP.values()) <= set(map(ord, 'ACGT')) | set(map(ord, IUPAC)), \
    'COMP の出力に IUPAC 定義の無い文字がある（相補鎖の照合が静かに 0 件になる）'

# motif, mark, modified_position（1 始まり。複数なら '/' 区切り）
MOTIFS = [
    ('AAGCCCG',  '4mC', '5'),
    ('AAGCCCG',  '6mA', '1'),
    ('CCGKCA',   '6mA', '6'),
    ('CGGCAACC', '6mA', '6'),
    ('GAACCGG',  '6mA', '2/3'),
    ('GAAGCCCG', '6mA', '2'),
    ('GCCGGC',   '4mC', '3'),
    ('HAAGCCCG', '6mA', '2'),
    ('TGGCCGGC', '4mC', '5'),
]

# GCCGGC 行だけは独立のセンサスではない。原稿の宣言（value_registry.py
# OCC.4mC.GCCGGC.*）どおり、分子は TGGCCGGC の判定部位数で、分母だけを
# 六塩基の鎖指定位置に置き換えた「配列可用性の分母」である。
# 認識配列は 8 塩基の TGGCCGGC なので、この行を修飾率として読んではならない。
NUMERATOR_FROM = {('GCCGGC', '4mC'): ('TGGCCGGC', '4mC')}

# 内部対照。座標規約・パイルアップの取り違えはここで落ちる。
EXPECT = {
    ('TGGCCGGC', '4mC'): {'genome_occurrences': 2093,
                          'methylated_sites': {'T1': 1886, 'T2': 1892, 'T3': 1866}},
    ('AAGCCCG', '6mA'):  {'genome_occurrences': 1334,
                          'methylated_sites': {'T1': 323, 'T2': 350, 'T3': 335}},
}


def load_genome() -> tuple[str, str]:
    name, chunks = None, []
    with REF.open() as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    raise SystemExit('参照が複数配列: この生成器は単一染色体前提')
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    return name, ''.join(chunks).upper()


def iupac_regex(motif: str) -> str:
    return ''.join(f'[{IUPAC[c]}]' if c in IUPAC else c for c in motif)


def revcomp(motif: str) -> str:
    return motif.translate(COMP)[::-1]


def motif_positions(genome: str, motif: str, offsets: list[int]) -> list[tuple[int, str]]:
    """鎖指定の修飾塩基位置（0 始まり, strand）を返す。

    ＋鎖: 出現開始 + offset
    −鎖: 相補配列の出現開始 + (len - 1 - offset)
    """
    out = []
    fwd = [m.start() for m in re.finditer(f'(?={iupac_regex(motif)})', genome)]
    rev = [m.start() for m in re.finditer(f'(?={iupac_regex(revcomp(motif))})', genome)]
    for s in fwd:
        for off in offsets:
            out.append((s + off, '+'))
    for s in rev:
        for off in offsets:
            out.append((s + len(motif) - 1 - off, '-'))
    return out, len(fwd) + len(rev)


def expected_by_composition(genome: str, motif: str) -> float:
    mono, L = Counter(genome), len(genome)
    p = 1.0
    for ch in motif:
        p *= sum(mono[b] for b in IUPAC.get(ch, ch)) / L
    return p * L * 2


def extract(whitelist: Path, workdir: Path) -> dict[str, dict[tuple[int, str], tuple[int, int, int]]]:
    """各ライブラリから、ホワイトリスト上の 4mC/6mA 行を awk で抜く。

    bedMethyl の列: 2=start(0 始まり) 4=mod_code 6=strand
                    10=Nvalid_cov 12=Nmod 16=Nfail
    """
    procs = {}
    workdir.mkdir(parents=True, exist_ok=True)
    for lib in sorted({l for ls in TIMEPOINTS.values() for l in ls}):
        src = PILEUP / f'{lib}_pileup.bed'
        if not src.exists():
            raise SystemExit(f'正準パイルアップが無い: {src}')
        dst = workdir / f'{lib}.tsv'
        awk = (r'NR==FNR{wl[$1"\t"$2]=1; next} '
               r'($4=="21839" || $4=="a") && (($2"\t"$6) in wl) '
               r'{print $2"\t"$6"\t"$4"\t"$10"\t"$12"\t"$16}')
        procs[lib] = (subprocess.Popen(['awk', awk, str(whitelist), str(src)],
                                       stdout=dst.open('w')), dst)
    tables: dict = {}
    for lib, (proc, dst) in procs.items():
        if proc.wait() != 0:
            raise SystemExit(f'awk 失敗: {lib}')
        t = {}
        with dst.open() as fh:
            for line in fh:
                pos, strand, code, nv, nm, nf = line.rstrip('\n').split('\t')
                t[(int(pos), strand, code)] = (int(nv), int(nm), int(nf))
        tables[lib] = t
    return tables


def census(tables, positions, code) -> dict[str, tuple[int, int]]:
    """時点ごとに (eligible, methylated) を返す。"""
    res = {}
    for tp, libs in TIMEPOINTS.items():
        elig = meth = 0
        for pos, strand in positions:
            nv = nm = nf = 0
            present = 0
            for lib in libs:
                row = tables[lib].get((pos, strand, code))
                if row is None:
                    continue
                present += 1
                nv += row[0]; nm += row[1]; nf += row[2]
            if present >= 1 and (nv + nf) >= DEPTH_MIN:
                elig += 1
                if nv > 0 and nm / nv >= LEVEL_CUT:
                    meth += 1
        res[tp] = (elig, meth)
    return res


def build(workdir: Path) -> list[dict]:
    name, genome = load_genome()
    sys.stderr.write(f'[ref] {name} {len(genome):,} bp\n')

    layout = {}
    for motif, mark, modpos in MOTIFS:
        offs = [int(x) - 1 for x in modpos.split('/')]
        pos, occ = motif_positions(genome, motif, offs)
        if len(set(pos)) != len(pos):
            raise SystemExit(f'{motif} {mark}: 鎖指定位置が重複している')
        layout[(motif, mark)] = (pos, occ, offs)

    wl = workdir / 'whitelist.txt'
    workdir.mkdir(parents=True, exist_ok=True)
    wl.write_text(''.join(f'{p}\t{s}\n' for p, s in
                          sorted({ps for pos, _, _ in layout.values() for ps in pos})))
    tables = extract(wl, workdir)

    # 先に全モチーフのセンサスを取る（GCCGGC 行は TGGCCGGC の結果を借りるので、
    # 行の組み立てより前に全部そろっていなければならない）
    cache = {}
    for motif, mark, _ in MOTIFS:
        pos, _, _ = layout[(motif, mark)]
        cache[(motif, mark)] = census(tables, pos, MODCODE[mark])

    rows = []
    for motif, mark, modpos in MOTIFS:
        pos, occ, _ = layout[(motif, mark)]
        c = cache[(motif, mark)]
        exp = expected_by_composition(genome, motif)
        for tp in TIMEPOINTS:
            elig, meth = c[tp]
            src = NUMERATOR_FROM.get((motif, mark))
            if src is not None:                       # 分子を別モチーフから借りる行
                meth = cache[src][tp][1]
            rows.append({
                'motif': motif, 'mark': mark, 'modified_position': modpos,
                'timepoint': tp,
                'genome_occurrences': occ,
                'expected_by_composition': round(exp, 1),
                'observed_over_expected': round(occ / exp, 3),
                'eligible_positions': elig,
                'methylated_sites': meth,
                'occupancy_pct': round(100 * meth / elig, 2) if elig else '',
                'occupancy_pct_over_genome_positions': round(100 * meth / len(pos), 2),
            })

    for key, exp in EXPECT.items():
        got_occ = layout[key][1]
        if got_occ != exp['genome_occurrences']:
            raise SystemExit(f'内部対照に失敗 {key}: 出現数 {got_occ} '
                             f'≠ {exp["genome_occurrences"]}（座標規約を疑うこと）')
        for tp, n in exp['methylated_sites'].items():
            if cache[key][tp][1] != n:
                raise SystemExit(f'内部対照に失敗 {key} {tp}: 部位数 '
                                 f'{cache[key][tp][1]} ≠ {n}')
    sys.stderr.write('[ok] 内部対照 通過\n')
    return rows


COLUMNS = ['motif', 'mark', 'modified_position', 'timepoint', 'genome_occurrences',
           'expected_by_composition', 'observed_over_expected', 'eligible_positions',
           'methylated_sites', 'occupancy_pct', 'occupancy_pct_over_genome_positions']


# 出荷ファイルの桁数に合わせる（O/E は 3 桁、占有率は末尾 0 を落とした 2 桁、
# 期待値は 1 桁）。桁数が変わると凡例・本文との文字列照合が外れる。
DECIMALS = {'observed_over_expected': 3, 'expected_by_composition': 1}


def render(rows: list[dict]) -> str:
    def fmt(col, v):
        if not isinstance(v, float):
            return str(v)
        if col in DECIMALS:
            return f'{v:.{DECIMALS[col]}f}'
        return f'{v:.2f}'.rstrip('0').rstrip('.') or '0'
    out = ['\t'.join(COLUMNS)]
    for r in sorted(rows, key=lambda r: (r['motif'], r['mark'], r['timepoint'])):
        out.append('\t'.join(fmt(c, r[c]) for c in COLUMNS))
    return '\n'.join(out) + '\n'


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, default=OUT)
    ap.add_argument('--workdir', type=Path,
                    default=Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures'
                                 '/scripts/.supptable17_work'))
    ap.add_argument('--check', action='store_true',
                    help='既存ファイルと差分を取るだけで書き換えない')
    args = ap.parse_args()

    rows = build(args.workdir)
    text = render(rows)

    if args.check:
        if not args.out.exists():
            print(f'既存ファイルが無い: {args.out}')
            return 1
        old = [l.split('\t') for l in args.out.read_text().strip().split('\n')]
        new = [l.split('\t') for l in text.strip().split('\n')]
        oldmap = {tuple(r[:4]): r for r in old[1:]}
        newmap = {tuple(r[:4]): r for r in new[1:]}
        ndiff = 0
        for key in sorted(set(oldmap) | set(newmap)):
            o, n = oldmap.get(key), newmap.get(key)
            if o is None or n is None:
                print(f'  [行の有無] {key}: ' + ('生成器のみ' if o is None else '既存のみ'))
                ndiff += 1
                continue
            for col, a, b in zip(old[0], o, n):
                if a != b:
                    print(f'  [差分] {" ".join(key)} {col}: 既存 {a} → 生成器 {b}')
                    ndiff += 1
        print(f'差分 {ndiff} 件 / 照合 {len(newmap)} 行')
        return 1 if ndiff else 0

    args.out.write_text(text)
    print(f'書き出し {args.out} （{len(rows)} 行）')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
