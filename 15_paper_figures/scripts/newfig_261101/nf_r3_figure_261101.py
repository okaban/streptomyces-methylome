# 2026-11-04 著者決定: 区別は軸ラベル・パネル内表記・形が担っているので、色で重ねて符号化しない（中立グレー）。
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R3 — where the modifications sit, and the promoter-proximal level ratio.

    python3 nf_r3_figure_261101.py

Inputs: tables/nf_r3_region_261101.tsv and tables/nf_r3_tssratio_261101.tsv,
written by nf_r3_recompute_261101.py. Nothing is hard-coded here.

Panel (a) answers Yuzawa R3 "where in the gene do the modified motifs sit":
the four declared region classes on the author's 2026-11-01 definition
(promoter = -300..TSS, 5'UTR = TSS..+50, CDS, other; TSS = the 3,570 measured
Jeong 2016 dRNA-seq TSSs). Each series is drawn twice - modified sites, and the
measurable positions of the same motif - so the reader can see that the
composition of the modified set is the composition of the motif's occurrences.

Panel (b) puts the one surviving spatial result on the same axis as its 4mC
comparators, computed on the same read set, the same TSS set, the same GC
matching and the same windows. It deliberately does NOT assert a mark-level
contrast: AAGCCCG carries a deficit under BOTH marks, and TGGCCGGC does not,
so the deficit tracks the motif rather than the mark. The retired framing
"the two marks have opposite spatial preferences" must not be reintroduced.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SUPPFIG_OUT", HERE)
SUPP = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014"
sys.path.insert(0, SUPP)
sys.path.insert(0, os.path.dirname(SUPP))

import matplotlib.pyplot as plt  # noqa: E402
from _suppstyle import (NAR_WIDTH_IN, apply_figure_style,  # noqa: E402
                        figure_panel_letters, report, save)

sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
import importlib  # noqa: E402
U = importlib.import_module("00_shared_utils")

TAB = os.path.join(HERE, "tables")
NAME = "NF-R3_local_methylation_261101.png"

REGION_ORDER = ["promoter", "utr5", "cds", "other"]
REGION_LABEL = {"promoter": "Promoter (\u2212300 to TSS)",
                "utr5": "5\u2032 UTR (TSS to +50)",
                "cds": "CDS", "other": "Other"}
# one ordinal ramp, darkest on the class the panel is about
REGION_COLOUR = {"promoter": U.COL_BAR_DARK, "utr5": U.COL_BAR_MID,
                 "cds": U.COL_BAR_LIGHT, "other": U.COL_BAR_PALE}
SERIES_ORDER = ["TGGCCGGC C5", "AAGCCCG C5", "AAGCCCG A1/A2"]
SERIES_MARK = {"TGGCCGGC C5": "4mC", "AAGCCCG C5": "4mC",
               "AAGCCCG A1/A2": "6mA"}


def panel_a(ax, reg):
    d = reg[reg.timepoint == "T1"]
    rows, ylab = [], []
    for s in SERIES_ORDER:
        sub = d[d.series == s].set_index("region")
        for which, col, ncol in (("modified", "pct_of_modified",
                                  "n_modified_total"),
                                 ("measurable", "pct_of_measurable",
                                  "n_measurable_total")):
            rows.append([float(sub.loc[r, col]) for r in REGION_ORDER])
            ylab.append(f"{s} {SERIES_MARK[s]}\n{which}, "
                        f"n = {int(sub.iloc[0][ncol]):,}")
    rows = np.array(rows)
    y = np.array([5.6, 4.85, 3.4, 2.65, 1.2, 0.45])
    left = np.zeros(len(rows))
    for i, r in enumerate(REGION_ORDER):
        ax.barh(y, rows[:, i], left=left, height=0.56,
                color=REGION_COLOUR[r], edgecolor="white", linewidth=0.5,
                label=REGION_LABEL[r])
        left += rows[:, i]
    for yy, v in zip(y, rows[:, 0]):
        ax.text(103, yy, f"{v:.1f}", va="center", ha="left", fontsize=6.0,
                color=U.COL_BAR_DARK)
    ax.text(103, 6.35, "Promoter\nshare (%)", va="center", ha="left",
            fontsize=5.8, color=U.COL_BAR_DARK)
    ax.set_yticks(y)
    ax.set_yticklabels(ylab, fontsize=6.0)
    ax.set_xlim(0, 122)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_ylim(-0.35, 6.9)
    ax.set_xlabel("Share of positions (%)")
    ax.legend(frameon=False, fontsize=5.8, loc="upper center",
              bbox_to_anchor=(0.46, -0.145), ncol=2, handlelength=1.0,
              columnspacing=0.8, borderpad=0.0, handletextpad=0.5)


def panel_b(ax, rat):
    order = [("6mA, AAGCCCG A1/A2", "6mA"), ("4mC, AAGCCCG C5", "4mC"),
             ("4mC, TGGCCGGC C5", "4mC")]
    ypos, ylab = [], []
    y = 0.0
    # 2026-11-07 共著者コメント C08:
    #  (1) 「グラフ内の印の違いに関する記載がない／印に違いを設ける必要はない」
    #      → ウィンドウ幅は y ラベルが明示しているので、記号は 1 種（塗り丸）に統一。
    #  (2) 「1.00 を中央にしていないのが違和感」→ x 軸を 1.00 対称にする。
    span = 0.0
    for uni, mark in order:
        for win in (100, 300):
            r = rat[(rat.universe == uni) & (rat.window_bp == win)].iloc[0]
            col = "#4a4a4a"  # 2026-11-04: y ラベルがマーク名を含むため色では分けない
            ax.plot([r.ci95_lo, r.ci95_hi], [y, y], color=col, lw=1.0,
                    solid_capstyle="butt", zorder=2)
            ax.plot([r.ratio], [y], marker="o", ms=4.0, color=col,
                    mfc=col, mew=0.9, zorder=3)
            span = max(span, abs(1.0 - r.ci95_lo), abs(1.0 - r.ci95_hi))
            ypos.append(y)
            short = uni.split(", ")[1] + " " + uni.split(",")[0]
            ylab.append(f"{short}\n\u00b1{win} bp, n = {int(r.n_near)}")
            y -= 1.0
        y -= 0.3
    ax.axvline(1.0, color="0.55", lw=0.7, ls=(0, (3, 2)), zorder=1)
    half = span * 1.12
    ax.set_xlim(1.0 - half, 1.0 + half)
    ax.set_xticks([round(1.0 + k * 0.05, 2) for k in (-2, -1, 0, 1, 2)])
    ax.set_yticks(ypos)
    ax.set_yticklabels(ylab, fontsize=6.0)
    ax.set_xlabel("Methylation level,\nTSS-proximal / distal (ratio)")
    ax.margins(y=0.10)


def main():
    apply_figure_style()
    reg = pd.read_csv(f"{TAB}/nf_r3_region_261101.tsv", sep="\t")
    rat = pd.read_csv(f"{TAB}/nf_r3_tssratio_261101.tsv", sep="\t")
    rat = rat[rat.universe != "4mC, all GCCGGC C3"]

    fig, axes = plt.subplots(1, 2, figsize=(NAR_WIDTH_IN, 2.97),
                             gridspec_kw={"width_ratios": [1.45, 1.0]})
    fig.subplots_adjust(left=0.150, right=0.985, top=0.9289, bottom=0.1957,
                        wspace=0.62)
    panel_a(axes[0], reg)
    panel_b(axes[1], rat)
    figure_panel_letters(fig, axes, ["a", "b"], dx=(0.143, 0.128), dy=0.004)
    p = save(fig, NAME)
    report(p)


if __name__ == "__main__":
    main()
