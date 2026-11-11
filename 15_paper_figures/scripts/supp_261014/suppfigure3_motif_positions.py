#!/usr/bin/env python3
"""Supplementary Figure 3 — where inside each recognition sequence the marks sit.

Rebuild of the 321.6 mm original at 174 mm column width (2026-10-14).  Numbers
are read verbatim from data/numbers.json (the canonical reproducible-rule
census); nothing is recomputed here.

    python suppfigure3_motif_positions.py
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from _suppstyle import (BLUE, DATA_DIR, GREY, LIGHT, NAR_WIDTH_IN, RED,
                        apply_figure_style, figure_panel_letters, report,
                        save)

N = json.load(open(DATA_DIR / "numbers.json"))

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))


def schematic(ax, top, bottom, marks, title, groups=()):
    """groups: iterable of (index_tuple, label, colour) drawn as one shared label."""
    ax.set_xlim(-0.95, len(top) + 0.75)
    ax.set_ylim(-0.1, 2.75)
    ax.axis("off")
    ax.text(-0.95, 2.62, title, ha="left", va="bottom", fontsize=7.0)
    for row, (seq, y) in enumerate([(top, 1.35), (bottom, 0.45)]):
        for i, b in enumerate(seq):
            col, tc = LIGHT, "0.15"
            if row == 0 and i in marks:
                col, tc = marks[i][1], "white"
            ax.add_patch(FancyBboxPatch((i + 0.08, y), 0.84, 0.66,
                                        boxstyle="round,pad=0.015,rounding_size=0.08",
                                        fc=col, ec="0.65", lw=0.5))
            ax.text(i + 0.5, y + 0.33, b, ha="center", va="center",
                    fontsize=6.4, color=tc, fontweight="bold")
        if row == 0:
            for i in range(len(seq)):
                ax.plot([i + 0.5, i + 0.5], [y, y - 0.24], color="0.65", lw=0.5)
    grouped = {i for idx, _, _ in groups for i in idx}
    for i, (lab, col) in marks.items():
        if i in grouped:
            continue
        ax.text(i + 0.5, 2.06, lab, ha="center", va="bottom",
                fontsize=6.0, color=col, fontweight="bold")
    for idx, lab, col in groups:
        x0, x1 = min(idx) + 0.12, max(idx) + 0.88
        ax.plot([x0, x0, x1, x1], [2.02, 2.14, 2.14, 2.02], color=col, lw=0.6,
                solid_joinstyle="miter")
        ax.text((x0 + x1) / 2, 2.18, lab, ha="center", va="bottom",
                fontsize=6.0, color=col, fontweight="bold")
    ax.text(-0.52, 1.68, "5'", fontsize=6.0, ha="center")
    ax.text(len(top) + 0.52, 1.68, "3'", fontsize=6.0, ha="center")
    ax.text(-0.52, 0.78, "3'", fontsize=6.0, ha="center")
    ax.text(len(top) + 0.52, 0.78, "5'", fontsize=6.0, ha="center")


def barpanel(ax, labels, values, colours, title, ymax):
    ax.bar(range(len(values)), values, color=colours, width=0.62)
    for i, v in enumerate(values):
        txt = f"{v:.2f}".rstrip("0").rstrip(".") if v < 1 else f"{v:.0f}"
        ax.text(i, v + ymax * 0.035, txt, ha="center", va="bottom", fontsize=6.0,
                color="0.45" if v < 1 else "0.15")
        if v < ymax * 0.01:
            ax.plot([i - 0.31, i + 0.31], [0, 0], color=colours[i], lw=1.4,
                    solid_capstyle="butt")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylim(0, ymax)
    ax.set_title(title, loc="left", fontsize=6.6)
    ax.set_ylabel("Mean modified fraction over\nmeasurable positions (%)")


def countpanel(ax, labels, values, colours, universe, title, nastr=None):
    xs = range(len(values))
    ax.bar(xs, values, color=colours, width=0.52)
    for i, v in enumerate(values):
        if v is None or np.isnan(v):
            ax.text(i, universe * 0.34, nastr, ha="center", va="center",
                    fontsize=6.0, color="0.55", linespacing=1.4)
        else:
            ax.text(i, v + universe * 0.05, f"{v:,}", ha="center", va="bottom",
                    fontsize=6.2)
    # 2026-10-28 著者指摘「点線と数字が重なっている」。母数の点線とその
    # ラベル（all N）を削除した。ラベルを右端へ逃がしても、線が y=universe、
    # 文字の基線が y=universe*1.02 なので、この文字サイズでは線が文字を貫く。
    # 衝突検査は線を文字として見ないので素通りしていた。
    # 母数はパネル題（"1,886 of 2,093 TGGCCGGC"）と凡例に書いてあるので、
    # 図に線で引く必要がない。
    ax.set_xticks(list(xs))
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.65, len(values) - 0.35)
    ax.set_ylim(0, universe * 1.22)
    ax.set_title(title, loc="left", fontsize=6.6)
    ax.set_ylabel("Methylated positions (T1)")


def build(path_name):
    # 2026-10-28 D122: 模式図パネル (a)(d) を削除し 2x3 から 2x2 に詰めた。
    # 削除理由: 平面図では鎖の極性を描けず、「同一鎖に 2 つの化学種」という
    # 本論文の中心主張を支えられない。両配列の二重鎖模式図は補足図 15 にあり、
    # 5'/3' 極性・塩基対・主修飾塩基の塗り円・A2 の白円まで描いた上位互換である。
    # 模式図が担っていた「どちらの配列の話か」は、各棒パネルの題で名指しする。
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 3.90))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.92],
                          hspace=0.80, wspace=0.52,
                          left=0.095, right=0.985, top=0.8829, bottom=0.0687)
    axes = []

    ax_a = fig.add_subplot(gs[0, 0]); axes.append(ax_a)
    t = N["perpos_T1_TGG"]
    # 表示は認識配列 TGGCCGGC の 1 起点位置。格納キー C1/C2/C5 は六塩基 GCCGGC の
    # 0 起点オフセットで、TGGCCGGC の位置 4 / 5 / 8 に対応する（位置 5 が修飾側）。
    barpanel(ax_a, ["C4", "C5", "C8"], [t["C1"], t["C2"], t["C5"]],
             [GREY, RED, GREY], "TGGCCGGC", 95)

    ax_b = fig.add_subplot(gs[0, 1]); axes.append(ax_b)
    countpanel(ax_b, ["4mC\nat C5", "6mA"], [N["TGG_called"]["T1"], np.nan],
               [RED, BLUE], N["universe_TGG"],
               f"{N['TGG_called']['T1']:,} of {N['universe_TGG']:,}\nTGGCCGGC",
               nastr="n.a.\nno A in motif")

    ax_c = fig.add_subplot(gs[1, 0]); axes.append(ax_c)
    a = N["perpos_T1_AAG"]
    barpanel(ax_c, ["A1", "A2", "C4", "C5", "C6"],
             [a["A0"], a["A1"], a["C3"], a["C4"], a["C5"]],
             [BLUE, BLUE, GREY, RED, GREY],
             "AAGCCCG", 95)

    ax_d = fig.add_subplot(gs[1, 1]); axes.append(ax_d)
    countpanel(ax_d, ["4mC\nat C5", "6mA\nat A1/A2"],
               [N["AAG_C4_called"]["T1"], N["AAG_6mA_union"]["T1"]],
               [RED, BLUE], N["universe_AAG"],
               f"{N['AAG_C4_called']['T1']:,} and {N['AAG_6mA_union']['T1']:,}\n"
               f"of {N['universe_AAG']:,} AAGCCCG")

    figure_panel_letters(fig, axes, "abcd",
                         dx=[0.105, 0.098, 0.105, 0.098], dy=0.006)
    return fig


if __name__ == "__main__":
    # 出荷名は SuppFigure1.png（補足図 3 -> 1 の改番、2026-10-24 に出力名を追随）。
    fig = build("SuppFigure3.png")
    # 二本鎖模式図の塩基文字 A/C/G/T は模式図の内容そのもの（マーク上に置くのが正）。
    report(save(fig, "SuppFigure3.png", allow_over_data=("A", "C", "G", "T")))
