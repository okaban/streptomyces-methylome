# 2026-11-05 D188: 補足図 12→9 の削減。出荷名 SuppFigure9.png -> SuppFigure8.png。
#!/usr/bin/env python3
"""Supplementary Figure 18 — two negative controls that separate the mark from the sequence.

Both panels test the same kind of claim: an apparent property of METHYLATION that
turns out to be a property of where the recognition SEQUENCE sits. The control in
each case is the motif's own measurable positions, not a genome-wide background.

Panel a asks whether the modified subset of each motif is distributed across genomic
regions differently from that motif's own positions. Panel b asks whether the
coding-strand expression bias of genes carrying modified sites is larger than the
bias of genes carrying the motif at all.

Scoped to the two negatives that are new in this work. Region composition of the
recognition sequences is Supplementary Figure 2, the promoter-proximal scarcity and
its GC control are Supplementary Figures 12 and 16, and the methylation-expression
null at cluster level is Supplementary Figure 14; none of those is repeated here.

    python suppfigure18_negative_controls.py
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from _suppstyle import (DATA_DIR, META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

A = pd.read_csv(DATA_DIR / "neg261015" / "region_modified_vs_own_T1.csv")
S = pd.read_csv(DATA_DIR / "neg261015" / "strand_orientation_expression_T1.csv")

# 2026-10-22 — 注記の数値は図が描く表そのものから導出する。
#
# それまでは neg_numbers.json という別ファイルから読んでいた。この json は
# 2026-10-15 11:06 のもので、region_modified_vs_own_T1.csv が正準規則で作り直された
# 2026-10-16 07:02（figure_manifest.tsv の「旧表は HAAGCCCG 行が 19 インスタンス
# 不足」）より古く、更新されないまま残っていた。結果として
#
#     注記 : smallest p = 0.068 (AAGCCCG)      <- 古い表に基づく値
#     パネル: AAGCCCG (5′ not G) で p = 0.02   <- 現在の表の値
#
# と、同じ図の中で最小 p が食い違っていた。現在の表の最小 p は HAAGCCCG 行の
# 0.0248（chi2 = 9.36、修飾 259 / 可測 561）で、json の 0.0676 はどの行とも一致せず、
# json の min_p_motif = "AAGCCCG" の行は現在 p = 0.1069 である。したがって古いのは
# json のほうで、パネルが正しい。
#
# 数値を json に書き直すだけでは同じ乖離がまた起きるので、注記の 4 つの数値
# （検定したモチーフ数・Bonferroni 閾値・最小 p・その モチーフ）はすべて A から
# 計算する。これで表を作り直せば注記も必ず追随する。
# 表示名。データキー（motif 列）は照合キーなので変えない。
DISPLAY_NAMES = {"HAAGCCCG": "AAGCCCG (5\u2032 not G)",
                 "GAAGCCCG": "AAGCCCG (5\u2032 G)"}

_imin = A.p.idxmin()
NUM = {
    "n_motifs": int(len(A)),
    "alpha": 0.05 / len(A),
    "min_p": float(A.p.loc[_imin]),
    "min_p_motif": str(A.motif.loc[_imin]),
}
assert not A.passes_bonferroni.any(), A[A.passes_bonferroni]
assert NUM["min_p"] > NUM["alpha"], NUM

C_MOD, C_ALL, C_REF = "#2b2b2b", "#b4b4b4", "#6b6b6b"  # FIG23: マーク色ではないので中立グレー

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 3.35))  # 2026-11-09: 図下の散文を凡例へ移し、空いた 0.80 in を詰めた
    gs = fig.add_gridspec(1, 2, width_ratios=[1.10, 1.00], wspace=0.95,
                          left=0.205, right=0.985, top=0.839, bottom=0.294)

    # (a) 修飾部分集合の領域分布を、そのモチーフ自身の位置分布と比べる
    ax_a = fig.add_subplot(gs[0, 0])
    d = A.sort_values("p")
    y = np.arange(len(d))
    ax_a.hlines(y, d.prom_all, d.prom_mod, color=META_GREY, lw=0.8, zorder=2)
    ax_a.scatter(d.prom_all, y, s=22, color=C_ALL, zorder=3)
    ax_a.scatter(d.prom_mod, y, s=22, color=C_MOD, zorder=4)
    ax_a.set_yticks(y)
    # motif 列は入力表の照合キー。表示名はここで与える（データキーは変えない）。
    # 2026-10-15: HAAGCCCG / GAAGCCCG は AAGCCCG の 5′ 文脈クラスであって別の認識配列では
    # ないので、図では同一モチーフの層として見えるように表示する
    # （DECIDED.6mA.DENOMINATOR.261015）。
    DISPLAY = DISPLAY_NAMES
    ax_a.set_yticklabels([f"{DISPLAY.get(m, m)}\nn = {int(n):,}   p = {pv:.3f}"
                          for m, n, pv in zip(d.motif, d.n_all, d.p)],
                         fontsize=6.0, linespacing=1.25)
    ax_a.set_xlabel("Share of sites in a promoter window (%)\n(\u2212300 bp to the TSS)")
    ax_a.set_xlim(0, 21)
    ax_a.set_ylim(-0.7, len(d) - 0.3)
    ax_a.set_title("", loc="left")

    # (b) 鎖方向と発現の関係は、修飾ではなく配列の性質。パネル a と同じ文法に揃える。
    ax_b = fig.add_subplot(gs[0, 1])
    order = ["AAGCCCG", "TGGCCGGC"]
    mod = S[S.set_en == "modified sites only"].set_index("motif").loc[order]
    alls = S[S.set_en == "all measurable positions"].set_index("motif").loc[order]
    yb = np.arange(len(order))
    ax_b.axvline(1.0, color=C_REF, lw=0.9, ls=(0, (4, 2)), zorder=2)
    ax_b.hlines(yb, alls.ratio, mod.ratio, color=META_GREY, lw=0.8, zorder=3)
    ax_b.scatter(alls.ratio, yb, s=22, color=C_ALL, zorder=4)
    ax_b.scatter(mod.ratio, yb, s=22, color=C_MOD, zorder=5)
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels(
        [f"{m}\nmodified p = {pm:.3f}\nall p = {pa:.3f}"
         for m, pm, pa in zip(order, mod.p, alls.p)],
        fontsize=6.0, linespacing=1.3)
    ax_b.set_xlabel("Median TPM, coding-strand\ndivided by antisense-only")
    ax_b.set_xlim(0.85, 2.05)
    ax_b.set_xticks([1.0, 1.25, 1.5, 1.75, 2.0])
    ax_b.set_ylim(-0.55, len(order) - 0.45)
    ax_b.set_title("", loc="left")

    # 両パネルで系列が共通なので図単位の凡例にする。軸内に置くとデータに乗る。
    fig.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor=C_MOD,
                               markersize=5, label="modified subset"),
                        Line2D([0], [0], marker="o", color="none", markerfacecolor=C_ALL,
                               markersize=5, label="all measurable positions of the same motif"),
                        Line2D([0], [0], color=C_REF, lw=0.9, ls=(0, (4, 2)),
                               label="no strand effect (panel b)")],
               frameon=False, fontsize=6.0, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, 0.077), handlelength=1.2, handletextpad=0.4,
               columnspacing=1.4)

    figure_panel_letters(fig, [ax_a, ax_b], "ab", dx=[0.185, 0.150], dy=0.022)
    return fig


if __name__ == "__main__":
    # 出荷名は SuppFigure9.png（補足図 18 -> 12 の改番、2026-10-24 に出力名を追随）。
    report(save(build(), "SuppFigure9.png"))
