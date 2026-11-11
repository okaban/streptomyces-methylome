# SUPERSEDED 2026-10-15: 著者指示により補足図 13 をマーク別（4mC / 6mA）に作り直した。
# 産出元は suppfigure13_bgc_by_mark.py。旧版は認識配列別（TGGCCGGC のみ）で、
# さらに全トラックで ax.set_yticks([]) を実行しており縦軸が読めなかった。
# このまま実行すると現行の SuppFigure13.png を旧版で上書きするため止めてある。
raise SystemExit(
    "SUPERSEDED 2026-10-15: suppfigure13_bgc_vs_density.py は suppfigure13_bgc_by_mark.py に置換済み。"
)

#!/usr/bin/env python3
"""Supplementary Figure 13 — predicted BGCs against the timepoint structure of
TGGCCGGC 4mC density.

Rebuild of the 298.0 mm original at 174 mm column width (2026-10-14).  Site
counts are re-derived from the deposited pileups (data/pileup_GCCGGC) under the
reproducible pooled rule and asserted against the canonical census in
data/numbers.json before anything is drawn.

    python suppfigure13_bgc_vs_density.py
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from _gccggc import called, genome, motif_positions, pooled
from _suppstyle import (DATA_DIR, NAR_WIDTH_IN, RED, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

N = json.load(open(DATA_DIR / "numbers.json"))
bgc = pd.read_csv(DATA_DIR / "bgc_regions_29.csv")

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))

g = genome()
GENOME_LEN = len(g)
_duplex, tgg = motif_positions(g)

BIN = 50_000
nbins = int(np.ceil(GENOME_LEN / BIN))
bins, counts = {}, {}
for tp in ("T1", "T2", "T3"):
    hits = called(pooled(tp), tgg)
    counts[tp] = len(hits)
    assert counts[tp] == N["TGG_called"][tp], f"{tp}: {counts[tp]} != {N['TGG_called'][tp]}"
    bins[tp] = np.bincount(np.array([p // BIN for p, _ in hits]), minlength=nbins)

CORE = (1.50e6, 7.17e6)
ORIC = 4.270310e6
EXPR, LOW, SIL = "#2a6f8e", "#d9a05b", "#9e4a4a"


def shade(ax):
    for x0, x1, col in [(0, CORE[0], "#f4e6dd"), (CORE[0], CORE[1], "#e7eee8"),
                        (CORE[1], GENOME_LEN, "#f4e6dd")]:
        ax.add_patch(Rectangle((x0 / 1e6, 0), (x1 - x0) / 1e6, 1,
                               transform=ax.get_xaxis_transform(), fc=col, ec="none",
                               zorder=0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 5.0))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 2.5], hspace=0.12,
                          left=0.135, right=0.985, top=0.915, bottom=0.115)
    x = (np.arange(nbins) * BIN + BIN / 2) / 1e6
    ymax = max(b.max() for b in bins.values()) * 1.12

    track_axes = []
    for k, tp in enumerate(("T1", "T2", "T3")):
        ax = fig.add_subplot(gs[k, 0])
        track_axes.append(ax)
        shade(ax)
        ax.bar(x, bins[tp], width=BIN / 1e6 * 0.96, color=RED, lw=0, zorder=3)
        ax.axvline(ORIC / 1e6, color="0.45", ls="--", lw=0.8, zorder=4)
        ax.set_xlim(0, GENOME_LEN / 1e6)
        ax.set_ylim(0, ymax)
        ax.set_yticks([])
        ax.set_xticks([])
        for sp in ("top", "right", "left", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.text(-0.012, 0.5, f"{tp}   {counts[tp]:,} sites", transform=ax.transAxes,
                ha="right", va="center", fontsize=6.4)

    a0 = track_axes[0]
    span = GENOME_LEN / 1e6
    a0.text(0.5 / span, 1.10, "left arm", transform=a0.transAxes, ha="center",
            va="bottom", fontsize=6.0, color="#9c7a5e")
    a0.text(3.00 / span, 1.10, "core", transform=a0.transAxes, ha="center",
            va="bottom", fontsize=6.0, color="#5d7a63")
    a0.text(7.95 / span, 1.10, "right arm", transform=a0.transAxes, ha="center",
            va="bottom", fontsize=6.0, color="#9c7a5e")
    a0.text(ORIC / 1e6 + 0.06, 1.10, "oriC", transform=a0.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=6.0, color="0.35")
    fig.text(0.135, 0.995, "TGGCCGGC 4mC sites per 50 kb bin", ha="left", va="top",
             fontsize=6.4, color="0.2")

    axb = fig.add_subplot(gs[3, 0])
    shade(axb)
    lanes = {"expressed": (3, EXPR), "low": (2, LOW), "silent": (1, SIL)}
    nclass = bgc.status.value_counts()
    for st, (y, col) in lanes.items():
        for r in bgc[bgc.status == st].itertuples():
            axb.add_patch(Rectangle((r.start / 1e6, y - 0.16),
                                    max((r.end - r.start) / 1e6, 0.035), 0.32,
                                    fc=col, ec="none", zorder=3))
        axb.text(-0.012, (y + 0.25) / 6.85, f"{st}  n={nclass[st]}",
                 transform=axb.transAxes, ha="right", va="center", fontsize=6.4,
                 color=col)
    axb.axvline(ORIC / 1e6, color="0.45", ls="--", lw=0.8, zorder=4)

    named = {"coelichelin": 4, "desferrioxamine": 10, "CDA": 11, "ACT": 12,
             "RED": 17, "geosmin": 19, "CPK/SCB1": 21, "SCO6431": 22}
    for lab, rn in named.items():
        r = bgc[bgc.region_number == rn].iloc[0]
        mid = (r.start + r.end) / 2 / 1e6
        axb.annotate(lab, xy=(mid, 3.20), xytext=(mid, 3.75), rotation=90,
                     ha="center", va="bottom", fontsize=6.0, color="0.2",
                     arrowprops=dict(arrowstyle="-", lw=0.5, color="0.55"))

    sil = bgc[bgc.status == "silent"].sort_values("start", ascending=False)
    for j, r in enumerate(sil.itertuples()):
        mid = (r.start + r.end) / 2 / 1e6
        axb.annotate(f"{r.known_cluster} (silent)", xy=(mid, 0.84),
                     xytext=(1.55, 0.60 - 0.55 * j), ha="left", va="center",
                     fontsize=6.0, color=SIL,
                     arrowprops=dict(arrowstyle="-", lw=0.5, color=SIL))

    axb.set_ylim(-1.45, 6.6)   # silent ラベル 3 行分の余地（x 目盛りラベルとの衝突回避）
    axb.set_xlim(0, GENOME_LEN / 1e6)
    axb.set_yticks([])
    for sp in ("top", "right", "left"):
        axb.spines[sp].set_visible(False)
    # 自動ロケータは xlim(8.667) の外側にも 9.0 を置き、ラベルだけがキャンバス外に
    # 描かれる。明示指定で上限内に収める。
    axb.set_xticks(list(range(0, 9)))
    axb.set_xlabel("Chromosomal position (Mb)")

    figure_panel_letters(fig, [a0, axb], "ab", dx=0.105, dy=[0.040, 0.004])
    source_note(fig,
                "Site rule: three libraries pooled per timepoint, pre-pooling true depth (Nvalid_cov + Nfail) >= 30, pooled modified fraction\n"
                f">= 50 % (deposited pileups); n = {counts['T1']:,} / {counts['T2']:,} / {counts['T3']:,} "
                "TGGCCGGC 4mC sites and all 29 antiSMASH regions.",
                y=0.008)
    return fig


if __name__ == "__main__":
    report(save(build(), "SuppFigure13.png"))
