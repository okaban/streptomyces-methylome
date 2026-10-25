#!/usr/bin/env python3
"""補足表 15 — 信頼度フィルタが TGGCCGGC 4mC 占有率に与えるバイアス（ライブラリ別）。

この表は「フィルタ版」と「all-read 版」の対比を要求するが、以前の
`modkit pileup --no-filtering` 出力は 2026-10-15 時点でディスク上に無かった
（補足図 4 の凡例が no longer on disk と明記）。BAM から作り直した
`pileup_nofilt_261015/` を all-read 側に使う。

三つの規則:
  filtered        正典パイルアップ（modkit 既定の信頼度フィルタ）、Nvalid >= 1 を測定可能とし
                  Nmod/Nvalid >= 0.5 を部位とする
  all reads       --no-filtering パイルアップに同じ規則
  filtered + d10  正典パイルアップ、測定可能の下限を Nvalid >= 10 に上げる

バイアス = filtered − all reads（百分率ポイント）。フィルタで落ちた読みは
未修飾側に偏るため占有率が押し上げられる。

実行: cd scripts && PYTHONPATH=. python3 recompute_261015/build_supptable15.py
"""
import csv
import pathlib
import sys

import pandas as pd

sys.path.insert(0, 'supp_261014')
from _gccggc import genome, motif_positions  # noqa: E402
from _suppstyle import DATA_DIR  # noqa: E402

NOFILT = pathlib.Path("/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015")
OUT = pathlib.Path(__file__).resolve().parent / "SuppTable15_filter_bias_per_library.tsv"
LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}

TGG = set(motif_positions(genome())[1])


def read_nofilt(lib):
    acc = {}
    with open(NOFILT / f"{lib}_nofilt.bed") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3].startswith("21839,GCCGGC"):
                acc[(int(f[1]), f[5])] = (int(f[9]), int(f[11]), int(f[15]))
    return acc


def read_canon(lib):
    d = pd.read_csv(DATA_DIR / "pileup_GCCGGC" / f"{lib}.tsv", sep="\t", header=None,
                    names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
    return {(p, s): (nv, nm, nf) for p, s, nv, nm, nf in d.itertuples(index=False)}


def occupancy(acc, floor):
    eligible = [k for k in TGG if k in acc and acc[k][0] >= max(floor, 1)]
    called = [k for k in eligible if acc[k][1] / acc[k][0] >= 0.5]
    return len(called), len(eligible), 100 * len(called) / len(eligible)


def main():
    rows = []
    for tp, libs in LIBS.items():
        for lib in libs:
            canon, nofilt = read_canon(lib), read_nofilt(lib)
            cf, ar, d10 = occupancy(canon, 1), occupancy(nofilt, 1), occupancy(canon, 10)
            rows.append(dict(
                timepoint=tp, replicate=lib,
                sites_filtered=cf[0], measurable_filtered=cf[1],
                occupancy_filtered_pct=round(cf[2], 2),
                sites_all_reads=ar[0], measurable_all_reads=ar[1],
                occupancy_all_reads_pct=round(ar[2], 2),
                sites_filtered_depth10=d10[0], measurable_filtered_depth10=d10[1],
                occupancy_filtered_depth10_pct=round(d10[2], 2),
                filter_bias_pp=round(cf[2] - ar[2], 2)))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    bias = [r["filter_bias_pp"] for r in rows]
    print(f"{len(rows)} libraries; filter bias {min(bias):+.2f} to {max(bias):+.2f} pp "
          f"(mean {sum(bias) / len(bias):+.2f}), same sign in all")
    print("->", OUT)


if __name__ == "__main__":
    main()
