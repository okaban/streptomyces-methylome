#!/usr/bin/env python3
"""補足図 10 — 複製に由来する被覆深度（複製数）勾配。

2026-10-15 の再構築（著者指示）。旧版は 3 パネル構成だったが:

  * 旧パネル b（深度帯で整合させたコア／腕比）と旧パネル c（コア対腕のレベル分布）は、
    「コアの 4mC 濃縮は複製数の人工産物ではない」という一つの防御の両輪だった。
    その主張は撤回済み。b は単一レプリケート・≥50% しきい値で計算されており、
    ある深度帯ではコア側の位置が 20 個しかなかった。フィルタ前のパイルアップで
    計算し直すと十分な位置数のある 8 帯でコア／腕比は 0.89-1.39 と 1.0 をまたぐ。
    c は撤回済みの深度フィルタ版サイト集合（n = 1,289）で描かれており、
    排除したかったフィルタ偏りを自分が被る循環になっていた。
    反証された論法の記録は本文の文章が担当する。
  * 旧パネル a は T1 しか描いておらず、凡例の T2/T3 の ori/ter 値に図の裏付けが無かった。
  * 旧版の軸ラベル・凡例は "GCCGGC m4C" 表記で、本文の呼称規約（TGGCCGGC / 4mC）
    と不一致だった。生成元スクリプトが無く画像内文字を直せなかった。

本版は 3 時点を 1 パネルで描く。主張は「T1 に複製由来の勾配があり、T2・T3 で平坦」。
この勾配は T2 で検出されるサイト集合が腕側へ広がることを生物学的変化として
解釈してはならない理由（検出力の変化と分離できない）になっている。

縦軸は各時点のゲノム全体中央値で正規化した被覆深度（marker frequency 表示）。
絶対深度は時点ごとのライブラリ深度に依存するため、正規化しないと勾配の形を
時点間で比較できない。絶対値は図注に記載する。

実行: cd scripts && PYTHONPATH=. python3 supp_261014/suppfigure10_coverage_gradient.py
"""
import pathlib
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, 'supp_261014')
from _suppstyle import (DATA_DIR, META_GREY, NAR_WIDTH_IN, OUT_DIR,  # noqa: E402
                        apply_figure_style,
                        report, save, source_note)

COV = DATA_DIR / "coverage261015" / "coverage_50kb_by_library.csv"
REGISTERED_SRC = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
                  "110_topup_recompute/tables/repl_profile_canonical.csv")

LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}
ORIC = 4.270310e6
CORE = (1.50e6, 7.17e6)
GENOME = 8_667_507
BIN = 50_000
ORITER_W = 200_000
REGISTERED = {"T1": 4.09, "T2": 1.08, "T3": 1.33}   # value_registry REPL.ORITER.* (確定)

# 時点はマークではないので 4mC(赤)/6mA(青) の配色は使わない。
# 勾配を持つ T1 を濃く、平坦な T2/T3 を淡く置いて主張をそのまま符号化する。
STYLE = {"T1": dict(color="#1a1a1a", lw=1.15, ls="-", zorder=5),
         "T2": dict(color="#6a6a6a", lw=0.95, ls=(0, (4, 1.6)), zorder=4),
         "T3": dict(color="#a6a6a6", lw=0.95, ls=(0, (1.3, 1.3)), zorder=3)}

cov = pd.read_csv(COV)
REF = pd.read_csv(REGISTERED_SRC, header=None,
                  names=["timepoint", "start", "size", "mean_true_depth"])


def pooled_profile(libs):
    """3 ライブラリを位置数で重み付けして合算した窓別平均真深度。"""
    sub = cov[cov.library.isin(libs)]
    g = sub.groupby("window_index")
    num = g.apply(lambda d: (d.mean_true_depth * d.n_positions).sum(), include_groups=False)
    den = g.n_positions.sum()
    return (num / den).sort_index()


def pooled_oriter(libs, width):
    """ori/ter は登録確定値 REPL.ORITER.* の出所である repl_profile_canonical.csv から読む。

    図・凡例・value_registry で同じ数値が出ることを保証するため、ここだけは独立に
    計算し直さない。定義はゲノム先頭から等間隔に切った 200 kb グリッドビンで、
    oriC を含むビン / 両終端ビンのうち浅い方。oriC を中心に置いた窓では約 5 % 下がる。
    本スクリプトの 50 kb プロファイルから同じグリッドで計算すると 4.01/1.07/1.33 となり、
    これは既登録の独立再計算値 REPL.ORITER.RECOMPUTE と一致する（深度定義の差。
    向きと桁はどの定義でも不変）。
    """
    tp = [k for k, v in LIBS.items() if v == libs][0]
    sub = REF[REF.timepoint == tp].sort_values("start")
    val = sub.mean_true_depth.values
    return float(val[int(ORIC // width)]), float(min(val[0], val[-1]))


apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))
fig, ax = plt.subplots(figsize=(NAR_WIDTH_IN, 2.55))

# 腕（コア外）を淡く塗る
for x0, x1 in [(0, CORE[0]), (CORE[1], GENOME)]:
    ax.axvspan(x0 / 1e6, x1 / 1e6, color="0.93", lw=0, zorder=0)

ax.axhline(1.0, color="0.70", lw=0.6, zorder=1)

medians, ratios = {}, {}
for tp, libs in LIBS.items():
    prof = pooled_profile(libs)
    med = float(np.median(prof.values))
    medians[tp] = med
    o, t = pooled_oriter(libs, ORITER_W)
    assert abs(o / t - REGISTERED[tp]) < 0.006, (tp, o / t, REGISTERED[tp])
    ratios[tp] = REGISTERED[tp]          # 表示は台帳の丸めに合わせる
    x = (prof.index.values * BIN + BIN / 2) / 1e6
    ax.plot(x, prof.values / med, label=f"{tp}  (ori/ter {ratios[tp]:.2f})", **STYLE[tp])

ax.axvline(ORIC / 1e6, color=META_GREY, lw=0.7, ls=(0, (2.2, 1.8)), zorder=2)
ax.annotate("oriC", xy=(ORIC / 1e6, ax.get_ylim()[0]), xytext=(2.5, 2.0),
            textcoords="offset points", ha="left", va="bottom", fontsize=6.0, color=META_GREY)

ax.set_xlim(0, GENOME / 1e6)
ax.set_xlabel("chromosomal position (Mb)")
ax.set_ylabel("coverage /\ngenome-wide median")
# 凡例はパネル上部に横 3 列で外置きする。腕の上端は「arm」注記のために空けておく。
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=3, frameon=False,
          fontsize=6.0, handlelength=2.0, columnspacing=1.8, borderpad=0.0,
          handletextpad=0.5)

# 腕／コアの注記（データの無い下端に置く）
yhi = ax.get_ylim()[1]
for xc, txt in [((0 + CORE[0]) / 2 / 1e6, "arm"), ((CORE[1] + GENOME) / 2 / 1e6, "arm")]:
    ax.annotate(txt, xy=(xc, yhi), xytext=(0, -2.0), textcoords="offset points",
                ha="center", va="top", fontsize=6.0, color="0.45")

fig.subplots_adjust(left=0.115, right=0.985, top=0.885, bottom=0.325)
source_note(fig,
    "Data source: modkit pileup (canonical, 2026-10-15); three libraries pooled per timepoint "
    "(T1 1-1/1-2/1-3, T2 2-1/2-3/2-4, T3 3-2/3-3/3-4).\n"
    "True depth = Nvalid_cov + Nfail at every C and A position on both strands, averaged in "
    f"{BIN // 1000} kb windows (n = {len(pooled_profile(LIBS['T1']))} windows). Genome-wide median\n"
    f"depth {medians['T1']:.1f}x / {medians['T2']:.1f}x / {medians['T3']:.1f}x at T1 / T2 / T3. "
    f"ori/ter = mean true depth in the {ORITER_W // 1000} kb grid bin containing oriC "
    "divided by that of the shallower terminal bin.", fontsize=5.8)

save(fig, "SuppFigure10")
report(OUT_DIR / "SuppFigure10.png")
print("median depth:", {k: round(v, 2) for k, v in medians.items()})
print("ori/ter 200kb:", {k: round(v, 3) for k, v in ratios.items()})
