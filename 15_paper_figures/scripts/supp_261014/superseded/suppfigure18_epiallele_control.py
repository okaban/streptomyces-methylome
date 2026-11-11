#!/usr/bin/env python3
"""Option B — standalone supplementary figure for the sequence-matched epiallele control.

Within-instance coupling of the three AAGCCCG target positions, measured against
non-motif triplets of identical geometry and read depth, and against expression
variability.  Companion to Supplementary Figure 15, which uses a distance-based
control instead.
    python suppfigure18_epiallele_control.py
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from _suppstyle import (DATA_DIR, META_GREY, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

D = DATA_DIR / "epi261015"
ex = pd.read_csv(D / "epi_excess_by_set_261015.csv")
cvt = pd.read_csv(D / "epi_expression_cv_261015.csv")
NUM = json.load(open(D / "epi_numbers.json"))

MOTIF, CTRL = "AAGCCCG instances", "non-motif control"
COL = {MOTIF: "#c0552b", CTRL: "#9aa7b0"}
mot = ex[ex.set == MOTIF]
ctl = ex[ex.set == CTRL]
THR = NUM["control_p95"]

apply_figure_style(frame="open", sizes=(7.0, 6.2, 6.0))


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 3.05))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.10, 0.95], wspace=0.52,
                          left=0.085, right=0.985, top=0.800, bottom=0.315)

    # (a) 超過の分布 — モチーフ対 対照
    ax_a = fig.add_subplot(gs[0, 0])
    rng = np.random.default_rng(7)
    for i, (lab, d) in enumerate(((CTRL, ctl), (MOTIF, mot))):
        x = i + rng.uniform(-0.18, 0.18, len(d))
        ax_a.scatter(x, d.excess, s=3.0, alpha=0.35, color=COL[lab], linewidths=0)
        m = d.excess.median()
        ax_a.plot([i - 0.30, i + 0.30], [m, m], color=COL[lab], lw=1.8, solid_capstyle="butt")
        ax_a.annotate(f"{m:.3f}", (i + 0.34, m), fontsize=6.2, color=COL[lab],
                      va="center", ha="left")
    ax_a.axhline(THR, color=META_GREY, lw=0.7, ls=(0, (3, 2)))
    ax_a.annotate("control 95th percentile", (1.58, THR + 0.012), fontsize=6.0,
                  color=META_GREY, va="bottom", ha="right")
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["non-motif\ntriplets\n(A,A,_,_,C)", "AAGCCCG\ninstances"])
    ax_a.set_xlim(-0.55, 1.62)
    ax_a.set_ylabel("Excess within-molecule structure\n(ΔH above permutation null)")
    ax_a.set_title("Coupling at the motif exceeds a\nsequence-matched technical control",
                   loc="left")

    # (b) 超過とリード深度 — 深度が効いていないこと
    ax_b = fig.add_subplot(gs[0, 1])
    for lab, d in ((CTRL, ctl), (MOTIF, mot)):
        ax_b.scatter(d.n, d.excess, s=3.0, alpha=0.35, color=COL[lab], linewidths=0,
                     label=lab.replace(" instances", "").replace("non-motif control",
                                                                 "non-motif triplets"))
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Molecules scored per instance")
    ax_b.set_ylabel("Excess within-molecule structure")
    ax_b.set_title("The gap is not a read-depth effect\n(means 141 and 147 molecules)",
                   loc="left")
    leg = ax_b.legend(frameon=False, fontsize=6.0, loc="upper left",
                      handlelength=0.9, borderpad=0.1, labelspacing=0.3,
                      handletextpad=0.45)
    for h in leg.legend_handles:
        h.set_sizes([9])

    # (c) 発現のばらつき — 帰無
    ax_c = fig.add_subplot(gs[0, 2])
    grp = [cvt[cvt.above_control_p95]["expression_cv_T1"],
           cvt[~cvt.above_control_p95]["expression_cv_T1"]]
    cols = ["#c0552b", "#9aa7b0"]
    for i, (d, c) in enumerate(zip(grp, cols)):
        x = i + rng.uniform(-0.16, 0.16, len(d))
        ax_c.scatter(x, d, s=3.0, alpha=0.35, color=c, linewidths=0)
        m = d.median()
        ax_c.plot([i - 0.28, i + 0.28], [m, m], color=c, lw=1.8, solid_capstyle="butt")
    ax_c.set_xticks([0, 1])
    ax_c.set_xticklabels([f"with coupled\ninstance\n(n = {NUM['cv_hi_n']})",
                          f"without\n(n = {NUM['cv_lo_n']})"])
    ax_c.set_xlim(-0.5, 1.5)
    ax_c.set_ylim(0, min(0.5, cvt["expression_cv_T1"].quantile(0.98)))
    ax_c.set_ylabel("Expression CV across\nthe three T1 replicates")
    ax_c.set_title(f"Coupling does not reach the\ntranscript ($p$ = {NUM['cv_p']:.2f})",
                   loc="left")

    figure_panel_letters(fig, [ax_a, ax_b, ax_c], "abc",
                         dx=[0.078, 0.075, 0.075], dy=0.020)
    source_note(
        fig,
        "Per-read calls from the canonical nine-library extraction; molecules carrying a call at all three AAGCCCG target positions (6mA position 1, 6mA position 2,\n"
        "4mC position 5), \u2265 20 molecules per instance, call probability \u2265 0.5. Excess = observed (H_independent \u2212 H_observed) minus the mean of 100 permutations\n"
        "that preserve each position's marginal and the read depth. Control triplets carry A, A and C at the same spacing and lie \u2265 60 bp from any motif instance.",
        y=0.012)
    return fig


if __name__ == "__main__":
    raise SystemExit(
        "NOT ADOPTED 2026-10-15: this candidate layout was not selected for\n         Supplementary Figure 15. The shipped figure is produced by\n         suppfigure15_optionC_rebuilt.py. build() is kept callable for reference.")
