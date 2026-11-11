#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 6, single-panel version (2026-11-04).

著者決定（2026-11-04）: 本文が引用していたのは Figure 6a だけで、(b)(c) は補足図 5・12 と
補足表 S5a/S5c/S5d を出典にしており、(d) は Results に対応する文が無かった。Figure 6 を
維持率の 1 パネルに絞り、(b)(c)(d) は補足へ移す。データ・計算・系統は 261102 版から一切
変更していない（panel_a をそのまま読み込んで描く）。
"""
import importlib, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SUPPFIG_OUT", HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014")
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")

import pandas as pd
import matplotlib.pyplot as plt
from _suppstyle import NAR_WIDTH_IN, apply_figure_style, report, save, source_note

src = importlib.import_module("nf_r5_figure_261102")
U = src.U
NAME = "NF-R5_retention_261104.png"


def main():
    apply_figure_style()

    # 2026-11-04 著者決定。4mC / 6mA の別は x 軸ラベルが示しているので、色で重ねて符号化しない。
    # FIG23 の規約どおり中立グレーに統一する（意味色は色だけが情報を担う図で使う）。
    src.RET_COL = ["#4a4a4a"] * 4
    ret = pd.read_csv(f"{src.TAB}/nf_r5_retention_261102.tsv", sep="\t")
    # 2026-11-07: 画像に焼き込んでいた 6 行の注記は、原稿の Figure 6 凡例が同じ内容を
    # そのまま述べており、図内テキストの二重化だった。凡例が唯一の出所になるよう図から外す
    # （figure-submission-conventions: 凡例で言えることを画像に入れない）。分母と
    # 判定規則は凡例と Supplementary Table S5a に残る。
    fig = plt.figure(figsize=(NAR_WIDTH_IN * 0.52, 2.10))
    ax = fig.add_axes([0.185, 0.185, 0.785, 0.780])
    src.panel_a(ax, ret)
    # 2026-11-07: 2 行に折り返していた軸名を 1 行に。分母（「最初の時点で修飾されていた
    # サイトのうち」）は凡例の 1 文目が述べている。
    ax.set_ylabel("Sites retained (%)")
    U.style_axes(ax, frame=True)
    p = save(fig, NAME)
    report(p)


if __name__ == "__main__":
    main()
