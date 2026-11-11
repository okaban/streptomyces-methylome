#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R5 — methylation against secondary metabolism (and, provisionally,
primary metabolism).

    python3 nf_r5_figure_261101.py

Reads the tables written by nf_r5_recompute_261101.py. Every level on this
figure is all-read (pileup_nofilt_261015); the confidence-filtered lineage is
not used, because it carries a +15.5 pp upward bias on level-scale statistics
(canon F7) and because its 3-2 / 3-4 libraries in `analysis/pileup/` are
missing the 2026-02 flowcell.

Panels (a) and (b) are null results and are drawn as such: the point of (a) is
that the per-cluster points sit on the identity line, and the point of (b) is
that the cloud has no usable slope. Neither panel asserts an effect, and the
legend must not turn either into one. Power is bounded, not zero: the paired
design detects a difference of about a third of the flank mean, no smaller.

Panel (c) is PROVISIONAL. There is no agreed gene-set definition of "primary
metabolism" in this project, so three candidates are drawn side by side to show
how little the picture moves between them. The author chooses; until then this
panel is not a result.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SUPPFIG_OUT", HERE)
SUPP = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014"
sys.path.insert(0, SUPP)
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from _suppstyle import (NAR_WIDTH_IN, apply_figure_style,  # noqa: E402
                        figure_panel_letters, report, save, source_note)
import importlib  # noqa: E402
U = importlib.import_module("00_shared_utils")

TAB = os.path.join(HERE, "tables")
SUPPTAB = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/"
           "SuppTable20a_cluster_methylation_expression_allread_261101.tsv")
NAME = "NF-R5_metabolism_261101.png"
MIN_POS = 5
COL = {"4mC": U.COL_4mC, "6mA": U.COL_6mA}
MK = {"4mC": "o", "6mA": "^"}
SHORT = {"P1 COG core metabolism (C,E,F,G,H,I)":
         "P1  COG core metabolism\n      C,E,F,G,H,I",
         "P2 COG metabolism + transport/envelope (C,E,F,G,H,I,P,M)":
         "P2  COG metabolism, wider\n      + P, M",
         "P3 all protein-coding genes outside antiSMASH regions":
         "P3  All CDS outside\n      antiSMASH regions",
         "(comparator) genes inside antiSMASH regions":
         "Comparator  CDS inside\n      antiSMASH regions"}


def panel_a(ax, cl, tests):
    for mark in ("4mC", "6mA"):
        d = cl[(cl["mark"] == mark) & (cl.n_interior >= MIN_POS)
               & (cl.n_flank >= MIN_POS)]
        g = d.groupby("region")[["level_interior", "level_flank"]].mean().dropna()
        pv = float(tests[tests.test ==
                         f"{mark} interior vs flank, cluster unit"]
                   .iloc[0].p_value)
        ax.scatter(g.level_flank, g.level_interior, s=15, marker=MK[mark],
                   facecolor="none", edgecolor=COL[mark], linewidth=0.8,
                   label=f"{mark}: {len(g)} clusters, signed-rank p = {pv:.2f}",
                   zorder=3)
    lim = (0.18, 0.80)
    ax.plot(lim, lim, color="0.55", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Mean level, flanks (fraction of reads)")
    ax.set_ylabel("Mean level, cluster interior\n(fraction of reads)")
    ax.legend(frameon=False, fontsize=5.2, loc="lower right",
              bbox_to_anchor=(1.045, -0.035), handlelength=0.8,
              handletextpad=0.35, borderpad=0.0, labelspacing=0.25)


def panel_b(ax, t20a, tests):
    for mark in ("4mC", "6mA"):
        col = "included_4mC" if mark == "4mC" else "included_6mA"
        m = t20a[t20a[col] == "yes"]
        xs, ys = [], []
        for tp in ("T1", "T2", "T3"):
            xs += list(m[f"TPM_{tp}"])
            ys += list(m[f"{mark}_level_{tp}"])
        ax.scatter(xs, ys, s=11, marker=MK[mark], facecolor="none",
                   edgecolor=COL[mark], linewidth=0.7, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 4000)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1k"])
    ax.set_ylim(0.15, 0.85)
    ax.set_xlabel("Cluster mean expression (TPM)")
    ax.set_ylabel("Cluster mean level\n(fraction of reads)")
    txt = []
    for mark in ("4mC", "6mA"):
        r = tests[tests.test == f"{mark} cluster mean level vs mean TPM"].iloc[0]
        txt.append(f"{mark}: \u03c1 = {float(r.statistic):+.2f}, "
                   f"p = {float(r.p_value):.2f}, n = {int(r.n)}")
    ax.text(0.03, 0.97, "\n".join(txt), transform=ax.transAxes, fontsize=5.8,
            va="top", ha="left", color="0.25")


def panel_c(ax, trial):
    order = list(SHORT)
    y = np.arange(len(order))[::-1]
    for mark in ("4mC", "6mA"):
        v = [float(trial[trial.definition == d][f"{mark}_mean_genebody_level"]
                   .iloc[0]) for d in order]
        ax.scatter(v, y, s=22, marker=MK[mark], facecolor=COL[mark],
                   edgecolor=COL[mark], label=mark, zorder=3)
    for yy in y:
        ax.axhline(yy, color="#EAEAEA", lw=0.5, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [SHORT[d] + f"\n      n = {int(trial[trial.definition == d].n_genes.iloc[0]):,}"
         for d in order], fontsize=5.5)
    ax.set_xlim(0.26, 0.76)
    ax.set_ylim(-0.95, len(order) - 0.45)
    ax.set_xlabel("Mean gene-body level at T1\n(fraction of reads)")
    ax.legend(frameon=False, fontsize=5.8, loc="upper center",
              bbox_to_anchor=(0.5, 1.035), ncol=2, handlelength=0.8,
              handletextpad=0.35, borderpad=0.0, columnspacing=1.2)


def main():
    apply_figure_style()
    cl = pd.read_csv(f"{TAB}/nf_r5_cluster_levels_261101.tsv", sep="\t")
    tests = pd.read_csv(f"{TAB}/nf_r5_cluster_tests_261101.tsv", sep="\t")
    trial = pd.read_csv(f"{TAB}/nf_r5_primary_metabolism_trial_261101.tsv",
                        sep="\t")
    t20a = pd.read_csv(SUPPTAB, sep="\t")

    fig = plt.figure(figsize=(NAR_WIDTH_IN, 2.95))
    ax_a = fig.add_axes([0.075, 0.250, 0.225, 0.675])
    ax_b = fig.add_axes([0.400, 0.250, 0.225, 0.675])
    ax_c = fig.add_axes([0.790, 0.250, 0.200, 0.675])
    panel_a(ax_a, cl, tests)
    panel_b(ax_b, t20a, tests)
    panel_c(ax_c, trial)
    figure_panel_letters(fig, [ax_a, ax_b, ax_c], ["a", "b", "c"],
                         dx=(0.068, 0.068, 0.175), dy=0.004)
    source_note(fig, "All levels are all-read (modkit --no-filtering, "
                     "pileup_nofilt_261015); 4mC = TGGCCGGC C5, 6mA = AAGCCCG "
                     "A1/A2.\nFlanks are two intervals of the cluster's own "
                     "length, one each side. The paired design detects a "
                     "difference of about a third of the flank mean, no "
                     "smaller.\n(c) Provisional: the gene-set definition of primary metabolism is an open author decision; the three candidates are shown so the choice can be made on what it changes.",
                y=0.008, fontsize=5.5)
    p = save(fig, NAME)
    report(p)


if __name__ == "__main__":
    main()
