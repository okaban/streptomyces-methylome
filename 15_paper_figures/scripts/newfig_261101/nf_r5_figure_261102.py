#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R5 (Figure 6 under the C04 six-main-figure layout), 2026-11-02.

    python3 nf_r5_figure_261102.py

Supersedes nf_r5_figure_261101.py, which stays on disk unchanged.

Changes against the 261101 version
----------------------------------
(a) NEW.  Per-mark site retention between timepoint pairs, the panel the C04
    reorganisation moves here from the dissolved Figure 4.  Recomputed on the
    all-read lineage by nf_r5_recompute_261102.py, on the published universes
    (2,093 / 1,334 / 1,334 / 565 motif positions) and the published call rule.
    The percentages drawn are the measurement, not the published numbers; they
    sit 3.65 pp below to 2.41 pp above them because the published panel used
    the confidence-filtered pileup, which this lane may not use.
(b), (c)  The 261101 panels (a) and (b), shifted down one letter.  Unchanged
    data, unchanged tables, still null results drawn as null results.
(d) REPLACES the provisional three-definition panel.  Author decision D155
    fixes primary metabolism as the core gene set of the KEGG metabolic
    pathways; the definition, release, retrieval date and gene count are in
    the source note and the legend.  antiSMASH cluster genes are kept as the
    comparator and the remaining protein-coding genes as the background.

Every level on this figure is all-read (modkit --no-filtering,
pileup_nofilt_261015); 4mC = TGGCCGGC C5, 6mA = AAGCCCG A1.  Panels (b) and
(c) assert no effect and the legend must not turn either into one: power is
bounded, not zero (the paired design in (b) detects a difference of about a
third of the flank mean, no smaller).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("SUPPFIG_OUT", HERE)
SUPP = "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts/supp_261014"
sys.path.insert(0, SUPP)
sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from _suppstyle import (NAR_WIDTH_IN, apply_figure_style,  # noqa: E402
                        figure_panel_letters, report, save, source_note)
import importlib  # noqa: E402
U = importlib.import_module("00_shared_utils")

TAB = os.path.join(HERE, "tables")
SUPPTAB = ("/Users/okaban/bioinfo/rna-seq/15_paper_figures/tables/supplementary/"
           "SuppTable20a_cluster_methylation_expression_allread_261101.tsv")
NAME = "NF-R5_metabolism_261102.png"
MIN_POS = 5
COL = {"4mC": U.COL_4mC, "6mA": U.COL_6mA}
MK = {"4mC": "o", "6mA": "^"}

RET_ORDER = ["TGGCCGGC 4mC (position 5)", "AAGCCCG 4mC (position 5)",
             "AAGCCCG 6mA (position 1, all instances)",
             "HAAGCCCG 6mA (position 1, H context)"]
RET_LAB = ["TGGCCGGC\n4mC", "AAGCCCG\n4mC", "AAGCCCG\n6mA", "HAAGCCCG\n6mA"]
RET_COL = [U.COL_4mC, mpl.colors.to_rgba(U.COL_4mC, 0.45), U.COL_6mA,
           mpl.colors.to_rgba(U.COL_6mA, 0.45)]
GROUPS = ["KEGG primary metabolism", "antiSMASH BGC genes",
          "Other protein-coding genes"]
GLAB = {"KEGG primary metabolism": "KEGG primary\nmetabolism",
        "antiSMASH BGC genes": "antiSMASH\ncluster genes",
        "Other protein-coding genes": "All other\nprotein-coding"}


# ------------------------------------------------------------------ panel a
def panel_a(ax, ret):
    for i, (mk, c) in enumerate(zip(RET_ORDER, RET_COL)):
        sub = ret[ret["mark"] == mk]
        assert len(sub) == 3, f"retention table: {mk} has {len(sub)} pairs"
        lo, hi = sub.retention_pct.min(), sub.retention_pct.max()
        ax.plot([i, i], [lo, hi], color=c, lw=1.1, zorder=1,
                solid_capstyle="round")
        for j, (_, r) in enumerate(sub.iterrows()):
            ax.plot(i + (j - 1) * 0.17, r.retention_pct, "o", ms=4.0, mfc=c,
                    mec="white", mew=0.6, zorder=3)
    ax.set_xticks(range(len(RET_ORDER)))
    ax.set_xticklabels(RET_LAB, fontsize=6.0)
    ax.set_yticks([75, 80, 85, 90, 95, 100])
    ax.set_ylim(74, 101)
    ax.set_xlim(-0.55, len(RET_ORDER) - 0.45)
    ax.set_ylabel("Sites retained (% of sites\nmodified at the first timepoint)")


# ------------------------------------------------------------------ panel b
def panel_b(ax, cl, tests):
    for mark in ("4mC", "6mA"):
        d = cl[(cl["mark"] == mark) & (cl.n_interior >= MIN_POS)
               & (cl.n_flank >= MIN_POS)]
        g = d.groupby("region")[["level_interior", "level_flank"]].mean().dropna()
        pv = float(tests[tests.test ==
                         f"{mark} interior vs flank, cluster unit"]
                   .iloc[0].p_value)
        ax.scatter(g.level_flank, g.level_interior, s=15, marker=MK[mark],
                   facecolor="none", edgecolor=COL[mark], linewidth=0.8,
                   label=f"{mark}: n = {len(g)}, p = {pv:.2f}",
                   zorder=3)
    lim = (0.18, 0.80)
    ax.plot(lim, lim, color="0.55", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Mean level, flanks (fraction of reads)")
    ax.set_ylabel("Mean level, cluster interior\n(fraction of reads)")
    ax.legend(frameon=False, fontsize=5.6, loc="lower right",
              bbox_to_anchor=(1.035, 0.022), handlelength=0.8,
              handletextpad=0.35, borderpad=0.0, labelspacing=0.25)


# ------------------------------------------------------------------ panel c
def panel_c(ax, t20a, tests):
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
    ax.set_ylim(0.15, 0.95)
    ax.set_xlabel("Cluster mean expression (TPM)")
    ax.set_ylabel("Cluster mean level\n(fraction of reads)")
    txt = []
    for mark in ("4mC", "6mA"):
        r = tests[tests.test == f"{mark} cluster mean level vs mean TPM"].iloc[0]
        txt.append(f"{mark}: \u03c1 = {float(r.statistic):+.2f}, "
                   f"p = {float(r.p_value):.2f}, n = {int(r.n)}")
    ax.text(0.03, 0.97, "\n".join(txt), transform=ax.transAxes, fontsize=5.8,
            va="top", ha="left", color="0.25")


# ------------------------------------------------------------------ panel d
def panel_d(ax, long, summ):
    ys = np.arange(len(GROUPS))[::-1]
    off = {"4mC": 0.19, "6mA": -0.19}
    for mark in ("4mC", "6mA"):
        for g, y in zip(GROUPS, ys):
            v = long[(long["mark"] == mark) & (long["group"] == g)
                     ]["level_T1"].dropna().to_numpy()
            bp = ax.boxplot([v], positions=[y + off[mark]], vert=False,
                            widths=0.30, showfliers=False, whis=(5, 95),
                            patch_artist=True, zorder=3)
            for b in bp["boxes"]:
                b.set(facecolor=mpl.colors.to_rgba(COL[mark], 0.22),
                      edgecolor=COL[mark], linewidth=0.7)
            for k in ("whiskers", "caps"):
                for e in bp[k]:
                    e.set(color=COL[mark], linewidth=0.7)
            for m_ in bp["medians"]:
                m_.set(color=COL[mark], linewidth=1.3)
    ax.set_yticks(ys)
    lab = []
    for g in GROUPS:
        n4 = int(summ[(summ["mark"] == "4mC") & (summ["group"] == g)]
                 ["n_genes_with_level"].iloc[0])
        n6 = int(summ[(summ["mark"] == "6mA") & (summ["group"] == g)]
                 ["n_genes_with_level"].iloc[0])
        ntot = int(summ[(summ["mark"] == "4mC") & (summ["group"] == g)]
                   ["n_genes_in_group"].iloc[0])
        lab.append(f"{GLAB[g]}\n{ntot:,} genes; measurable\n4mC {n4:,} / 6mA {n6:,}")
    ax.set_yticklabels(lab, fontsize=5.3)
    ax.set_ylim(-0.62, len(GROUPS) - 0.38)
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("Mean gene-body level at T1 (fraction of reads)")
    h = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(COL[m], 0.22),
                           edgecolor=COL[m], linewidth=0.7, label=m)
         for m in ("4mC", "6mA")]
    ax.legend(handles=h, frameon=False, fontsize=5.8, loc="lower right",
              bbox_to_anchor=(1.02, -0.03), ncol=2, handlelength=0.9,
              handletextpad=0.35, borderpad=0.0, columnspacing=1.0)
    p = {m: float(summ[(summ["mark"] == m) &
                       (summ["group"].str.startswith("TEST"))]
                  ["median_level"].iloc[0]) for m in ("4mC", "6mA")}
    ax.text(0.0, 1.035,
            "KEGG vs cluster genes, Mann\u2013Whitney:\n"
            f"4mC p = {p['4mC']:.2f},  6mA p = {p['6mA']:.2f}",
            transform=ax.transAxes, fontsize=5.6, va="bottom", ha="left",
            color="0.25")


# ------------------------------------------------------------------ main
def main():
    apply_figure_style()
    ret = pd.read_csv(f"{TAB}/nf_r5_retention_261102.tsv", sep="\t")
    cl = pd.read_csv(f"{TAB}/nf_r5_cluster_levels_261101.tsv", sep="\t")
    tests = pd.read_csv(f"{TAB}/nf_r5_cluster_tests_261101.tsv", sep="\t")
    long = pd.read_csv(f"{TAB}/nf_r5_genesets_levels_261102.tsv", sep="\t")
    summ = pd.read_csv(f"{TAB}/nf_r5_geneset_summary_261102.tsv", sep="\t")
    t20a = pd.read_csv(SUPPTAB, sep="\t")
    meta = json.load(open(f"{TAB}/nf_r5_kegg_meta_261102.json"))

    fig = plt.figure(figsize=(NAR_WIDTH_IN, 5.05))
    ax_a = fig.add_axes([0.100, 0.620, 0.325, 0.325])
    ax_b = fig.add_axes([0.590, 0.620, 0.325, 0.325])
    ax_c = fig.add_axes([0.100, 0.228, 0.325, 0.265])
    ax_d = fig.add_axes([0.590, 0.228, 0.325, 0.265])
    panel_a(ax_a, ret)
    panel_b(ax_b, cl, tests)
    panel_c(ax_c, t20a, tests)
    panel_d(ax_d, long, summ)
    for ax in (ax_a, ax_b, ax_c, ax_d):
        U.style_axes(ax, frame=True)
    figure_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d], ["a", "b", "c", "d"],
                         dx=(0.085, 0.085, 0.085, 0.085), dy=0.010)
    source_note(
        fig,
        "All levels are all-read (modkit --no-filtering, "
        "pileup_nofilt_261015); 4mC = TGGCCGGC C5, 6mA = AAGCCCG A1;\n"
        "pooled Nvalid \u2265 30 per timepoint. (a) Universes 2,093 / 1,334 / "
        "1,334 / 565 motif positions; a site counts as modified\n"
        "when pooled Nmod/Nvalid \u2265 0.5; points are the three timepoint "
        "pairs and the line spans them. (b) Flanks are two\n"
        "intervals of the cluster's own length, one each side, and p is a "
        "two-sided Wilcoxon signed-rank test with the\n"
        "cluster as the unit; the paired design detects a difference of about "
        "a third of\n"
        "the flank mean, no smaller. (b) and (c) show no effect. (d) Primary "
        "metabolism = KEGG BRITE br:sco00001,\n"
        "A09100 Metabolism, excluding B-categories 09109, 09110, 09111 and "
        f"09112 ({meta['n_primary_core_mapped']:,} genes mapped to locus\n"
        f"tags; KEGG {meta['kegg_release']}, retrieved "
        f"{meta['retrieved_on']}). Boxes are quartiles, whiskers the "
        "5th\u201395th centile.",
        y=0.008, fontsize=5.4)
    p = save(fig, NAME)
    report(p)


if __name__ == "__main__":
    main()
