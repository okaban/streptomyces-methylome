#!/usr/bin/env python3
"""Supplementary Figure 3 — 4mC calling is motif-exclusive and the hexamer is
marked on one strand rather than both.

ファイル名の 6 は改番前の番号。出荷名は SuppFigure3.png で、凡例
「4mC calling is motif-exclusive, and a methylated GCCGGC hexamer carries the mark
on one strand rather than both」に対応する現役の生成器である（figure_manifest.tsv
も SuppFigure3.png の生成器としてこのファイルを記録している）。2026-10-22 に
出力名が SuppFigure6.png のままだったことから退役図と誤認して superseded/ へ
隔離されたが、2026-10-24 に内容照合のうえ復帰させ、出力名を出荷名に合わせた。

Rebuild of the 387.1 mm original at 174 mm column width (2026-10-14).  Panel a
is read verbatim from data/numbers.json; panel b re-derives the per-strand
levels from the deposited pileups (data/pileup_GCCGGC) under the same site rule.

    python suppfigure6_motif_exclusive_asymmetry.py
"""
import json

import matplotlib.pyplot as plt
import numpy as np

from _gccggc import DEPTH_MIN, eligible, genome, level, motif_positions, pooled
from _suppstyle import (DATA_DIR, GREY, NAR_WIDTH_IN, RED, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

N = json.load(open(DATA_DIR / "numbers.json"))
S = N["s13"]

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))

g = genome()
duplex, _tgg = motif_positions(g)
acc = pooled("T1")
pairs = [(level(acc, (p, "+")), level(acc, (m, "-")))
         for p, m in duplex
         if eligible(acc, (p, "+")) and eligible(acc, (m, "-"))]
fp = np.array([a for a, _ in pairs])
fm = np.array([b for _, b in pairs])
assert len(fp) == S["occ_eligible"], f"{len(fp)} eligible duplexes, expected {S['occ_eligible']}"


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 3.05))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.74, 1.0], wspace=0.40,
                          left=0.105, right=0.985, top=0.845, bottom=0.235)

    ax_a = fig.add_subplot(gs[0, 0])
    pct_nm = 100 * S["nonmotif_ge50"] / S["nonmotif_eligible"]
    ax_a.bar([0, 1], [pct_nm, S["motif_pct"]], color=[GREY, RED], width=0.5)
    ax_a.plot([-0.25, 0.25], [0, 0], color=GREY, lw=1.8, solid_capstyle="butt")
    ax_a.text(0, 0.22, f"{S['nonmotif_ge50']} of {S['nonmotif_eligible'] / 1e6:.2f} M\n({pct_nm:.3f} %)",
              ha="center", va="bottom", fontsize=6.2, linespacing=1.4)
    ax_a.text(0, 1.55,
              f"{S['nonmotif_ge25']} at $\\geq$ 25 %,\n{S['nonmotif_ge10']} at $\\geq$ 10 %",
              ha="center", va="bottom", fontsize=6.0, color="0.45", linespacing=1.4)
    ax_a.text(1, S["motif_pct"] + 0.22,
              f"{S['motif_called']:,} of {S['motif_eligible']:,}\n({S['motif_pct']:.2f} %)",
              ha="center", va="bottom", fontsize=6.2, linespacing=1.4)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["non-motif\ncytosines", "hexamer GCCGGC\nmodified C"])
    ax_a.set_ylabel("Cytosines called 4mC-positive (%)")
    ax_a.set_ylim(0, 7.2)
    ax_a.set_xlim(-0.6, 1.6)
    ax_a.set_title("", loc="left",
                   fontsize=6.6)

    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.scatter(fp * 100, fm * 100, s=1.6, color=RED, alpha=0.30, lw=0, rasterized=True)
    ax_b.axhline(50, color="0.6", ls=":", lw=0.7)
    ax_b.axvline(50, color="0.6", ls=":", lw=0.7)
    ax_b.text(99, 57, f"both strands\n{S['both']} ({S['both_pct']:.1f} %)",
              ha="right", va="bottom", fontsize=6.2)
    ax_b.text(99, 43, f"one strand only\n{S['one']:,} ({S['one_pct']:.1f} %)",
              ha="right", va="top", fontsize=6.2)
    ax_b.set_xlabel("+ strand 4mC frequency (%)")
    ax_b.set_ylabel("$-$ strand 4mC frequency (%)")
    ax_b.set_xlim(-3, 103)
    ax_b.set_ylim(-3, 103)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_yticks([0, 25, 50, 75, 100])
    ax_b.set_title("", loc="left", fontsize=6.6)

    figure_panel_letters(fig, [ax_a, ax_b], "ab", dx=[0.088, 0.072], dy=0.010)
    source_note(fig,
                "Source: deposited modkit pileups of the nine modified-base BAMs (T1 = libraries 1-1/1-2/1-3). Site rule: three libraries pooled\n"
                f"per timepoint, pre-pooling true depth (Nvalid_cov + Nfail) >= {DEPTH_MIN}, pooled modified fraction >= 50 %. "
                f"n = {S['nonmotif_eligible']:,} eligible non-motif\ncytosines (both recognition sequences excluded) and "
                f"{S['occ_eligible']:,} hexamer GCCGGC occurrences with both strands eligible.",
                y=0.010)
    return fig


if __name__ == "__main__":
    # 積み上げ棒の内部に置く構成比ラベル。規則 FIG20 の allow_over_data として明示宣言。
    inbar = (f"both strands\n{S['both']} ({S['both_pct']:.1f} %)",
             f"one strand only\n{S['one']:,} ({S['one_pct']:.1f} %)")
    # 出荷名は SuppFigure3.png（補足図 6 -> 3 の改番、2026-10-24 に出力名を追随）。
    # SuppFigure6.png は現在 GC 統制図（suppfigure6b_gc_control_261022.py）が占める。
    report(save(build(), "SuppFigure3.png", allow_over_data=inbar))
