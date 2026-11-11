#!/usr/bin/env python3
"""Supplementary Figure 15 — single-molecule concordance of the two marks.

Rebuild of the 258.0 mm original at 174 mm column width (2026-10-14).  The
per-read statistics in data/h26/ were regenerated from the nine modified-base
BAMs under the original convention and reproduce the published values exactly
(0.470 / 0.171 / 0.016; 3.0 % and 3.3 %) — see data/h26/README.md.

    python suppfigure15_single_molecule.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _suppstyle import (DATA_DIR, META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

H = pd.read_csv(DATA_DIR / "h26" / "H26.csv")
PCd = pd.read_csv(DATA_DIR / "h26" / "PC26.csv")
D26 = pd.read_csv(DATA_DIR / "h26" / "D26.csv")

six = H[H["mark"].str.contains("6mA")]
items = [("bistable\nswitch\n(prediction)", six.bistable_pred.mean(), "#9aa7b0", None),
         ("same AAGCCCG\nA1\u2013A2\n(1 bp)", PCd.excess.mean(), "#c0552b", PCd.excess.std()),
         ("between\nmolecules\n(>10 bp)", six.excess_conc.mean(), "#e8b4a0",
          six.excess_conc.std())]
assert round(items[0][1], 3) == 0.470 and round(items[1][1], 3) == 0.171 \
    and round(items[2][1], 3) == 0.016, [round(i[1], 3) for i in items]

MC = {"GCCGGC 4mC": "#1b4965", "AAGCCCG 6mA": "#c0552b"}
ORDER = ["GCCGGC 4mC", "AAGCCCG 6mA"]

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 2.85))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.05, 0.95], wspace=0.58,
                          left=0.095, right=0.985, top=0.760, bottom=0.265)

    ax_a = fig.add_subplot(gs[0, 0])
    for i, (lab, v, c, sd) in enumerate(items):
        ax_a.bar(i, v, width=0.6, color=c, edgecolor="none")
        if sd:
            ax_a.errorbar(i, v, yerr=sd, fmt="none", ecolor="black", elinewidth=0.8,
                          capsize=2.0)
        ax_a.text(i, v + 0.020, f"{v:.3f}", ha="center", fontsize=6.0)
    ax_a.set_xticks(range(3))
    ax_a.set_xticklabels([x[0] for x in items], fontsize=6.0, linespacing=1.35)
    ax_a.set_ylabel("Excess concordance\non one molecule")
    ax_a.set_ylim(0, 0.56)
    ax_a.set_title("AAGCCCG 6mA: real coupling is\ndetected within an instance,\nnot between molecules",
                   loc="left", fontsize=6.6)

    ax_b = fig.add_subplot(gs[0, 1])
    for mk in ORDER:
        g = D26[(D26["mark"] == mk) & (D26.dist_max <= 15000)]
        s = g.groupby("dist_max")["excess"].agg(["mean", "std"])
        ax_b.errorbar(s.index, s["mean"], yerr=s["std"], fmt="o-", color=MC[mk],
                      lw=1.0, ms=3.0, capsize=1.8, label=mk)
    ax_b.axhline(0, color=META_GREY, lw=0.7)
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Distance between sites on\none molecule (bp)")
    ax_b.set_ylabel("Excess concordance")
    ax_b.set_ylim(-0.005, 0.048)
    ax_b.set_title("Flat from 200 bp to 15 kb —\nno processive action by\neither enzyme",
                   loc="left", fontsize=6.6)
    ax_b.legend(frameon=False, fontsize=6.0, loc="upper left", handlelength=1.2,
                borderaxespad=0.1, labelspacing=0.25)

    ax_c = fig.add_subplot(gs[0, 2])
    for i, mk in enumerate(ORDER):
        g = H[H["mark"] == mk]
        ax_c.scatter([i] * len(g), g.pct_of_bistable, s=12, color=MC[mk],
                     edgecolor="white", linewidth=0.4, zorder=3)
        ax_c.hlines(g.pct_of_bistable.mean(), i - 0.24, i + 0.24, color=MC[mk],
                    lw=1.4, zorder=2)
        ax_c.text(i + 0.30, g.pct_of_bistable.mean() + 4,
                  f"{g.pct_of_bistable.mean():.1f} %", ha="left", va="center",
                  fontsize=6.0, color=MC[mk])
    ax_c.axhline(100, color=META_GREY, lw=0.8, ls="--")
    ax_c.text(0.5, 103, "a clean ON/OFF switch = 100 %", ha="center", fontsize=6.0,
              color=META_GREY)
    ax_c.set_xticks([0, 1])
    ax_c.set_xticklabels(["GCCGGC\n4mC", "AAGCCCG\n6mA"], fontsize=6.0)
    ax_c.set_ylabel("Observed excess, as % of the\nbistable-switch prediction")
    ax_c.set_ylim(0, 118)
    ax_c.margins(x=0.3)
    ax_c.set_title("Molecule-level coherence is\n3 % of what an epiallele\nwould produce",
                   loc="left", fontsize=6.6)

    figure_panel_letters(fig, [ax_a, ax_b, ax_c], "abc",
                         dx=[0.088, 0.082, 0.082], dy=0.020)
    source_note(fig,
                "Per-read modification calls from the MM/ML tags of the nine modification-tagged alignments (pysam modified_bases, argmax\n"
                "against the canonical base); error bars are the SD over the nine libraries.",
                y=0.008)
    return fig


if __name__ == "__main__":
    raise SystemExit(
        "SUPERSEDED 2026-10-15: Supplementary Figure 15 is now produced by\n"
        "suppfigure15_optionC_rebuilt.py, which leads with the per-molecule\n"
        "distribution and restricts the 6mA panels to HAAGCCCG. build() is kept\n"
        "callable here for reference; it no longer writes the shipped PNG.")
