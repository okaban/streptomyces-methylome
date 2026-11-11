#!/usr/bin/env python3
"""補足表 S2a (SuppTable15_per_library_depth_and_sites.tsv) — ライブラリ別の深度と
TGGCCGGC 4mC 部位数。

2026-11-02 新設。この表はこれまで生成器を持たず手で維持されており、1-1 行だけが
不採用の復元版パイルアップ系統 (superseded_inputs.tsv の 2 行目、severity BLOCKER) の
値 1,182 を保持していた。正準系統 (`pileup_canonical_261015/`、1-1 は
cursor_results/20260108 の無傷の原本) では 1,189。

入力は正準パイルアップの npz キャッシュ `pileup_cache_261015/` のみ
(build_pileup_cache_261015.py が pileup_canonical_261015/ だけを読む)。

規則 (canon・value_registry と同一のライブラリ単位版):
  測定可能   真深度 Nvalid_cov + Nfail >= 30、参照塩基一致 (+鎖 C / −鎖 G)
  部位       Nmod / Nvalid >= 0.5
  深度中央値 GCCGGC 六塩基の 32,118 シトシン位置における真深度の中央値

プール後の時点別 census (T1 1,886 / T2 1,892 / T3 1,866) は同じ規則をプール後に
適用したもので、1-1 の系統選択に影響されない。

実行: cd 15_paper_figures/scripts && PYTHONPATH=. python3 \
        recompute_261015/build_supptable15_depth_and_sites_261102.py
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "supp_261014"))
from _gccggc import genome, motif_positions  # noqa: E402

CACHE = pathlib.Path(__file__).resolve().parent / "pileup_cache_261015"
OUT = (pathlib.Path(__file__).resolve().parents[2] / "tables" / "supplementary"
       / "SuppTable15_per_library_depth_and_sites.tsv")
TP = {"1-1": "T1", "1-2": "T1", "1-3": "T1", "2-1": "T2", "2-3": "T2", "2-4": "T2",
      "3-2": "T3", "3-3": "T3", "3-4": "T3"}
TOPUP = {"3-2": "yes (62.8% of reads from the 2026-02 flow cells)", "3-4": "yes (16.7%)"}
MIN_TRUE_DEPTH = 30


def lookup(key_sorted, order, query):
    """query の各要素について key_sorted 内の一致行インデックスと一致フラグを返す。"""
    at = np.searchsorted(key_sorted, query)
    at_c = np.clip(at, 0, len(key_sorted) - 1)
    hit = (at < len(key_sorted)) & (key_sorted[at_c] == query)
    return order[at_c], hit


def main():
    g = genome()
    duplexes, tgg = motif_positions(g)
    tgg_key = np.array(sorted(p * 2 + (st == "+") for p, st in tgg), dtype=np.int64)
    gcc_pos = np.array(sorted({p for pair in duplexes for p in pair}), dtype=np.int64)
    assert len(tgg_key) == 2093 and len(gcc_pos) == 32118, (len(tgg_key), len(gcc_pos))

    rows = []
    for lib, tp in TP.items():
        z = np.load(CACHE / f"{lib}.npz")
        m = (z["code"] == 1) & z["refok"]                 # 4mC チャネル、参照塩基一致
        pos, plus = z["start"][m].astype(np.int64), z["plus"][m]
        nv, nm, nf = z["Nvalid"][m], z["Nmod"][m], z["Nfail"][m]

        key = pos * 2 + plus.astype(np.int64)
        order = np.argsort(key)
        sel, hit = lookup(key[order], order, tgg_key)
        true_depth = np.where(hit, nv[sel] + nf[sel], -1)
        ratio = np.where(hit & (nv[sel] > 0), nm[sel] / np.maximum(nv[sel], 1), 0.0)
        sites = int(((true_depth >= MIN_TRUE_DEPTH) & (ratio >= 0.5)).sum())

        order_p = np.argsort(pos)
        sel_p, hit_p = lookup(pos[order_p], order_p, gcc_pos)
        depth = (nv[sel_p] + nf[sel_p])[hit_p]
        rows.append((lib, tp, int(np.median(depth)), sites, TOPUP.get(lib, "no")))

    with OUT.open("w", newline="") as fh:
        fh.write("library\ttimepoint\tmedian_true_depth_x\tTGGCCGGC_4mC_sites_called"
                 "\ttop_up_flow_cells\r\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\r\n")
    for r in rows:
        print("\t".join(str(x) for x in r))
    print("->", OUT)


if __name__ == "__main__":
    main()
