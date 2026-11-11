#!/usr/bin/env python3
"""Merged supplementary figure — chromosomal occupancy of both methylation
systems against their own motif-conditioned nulls, and the k-mer cross-talk
control for the AAGCCCG dual modification.

Consolidation of two shipped supplementary figures (author approval 2026-11-03):

    old Supplementary Figure 2  "Both methylation systems occupy the chromosome
                                 in proportion to their recognition sequences"
        generator: supp_261014/suppfigure4_core_share.py
    old Supplementary Figure 4  "The AAGCCCG dual modification is motif-specific,
                                 not a basecaller k-mer artefact"
        generator: 11_epigenome_integration/analysis/78_reviewer_robustness/
                   scripts/P1_figures.py  (its "Fig A" block -> P1_1_crosstalk.png)

The panels are re-laid on ONE canvas; neither original generator is called.
Every panel's computation is the code of its original generator, copied
unchanged, so no plotted number moves.  Panel letters are renumbered:

    old S2 a,b,c,d  ->  a,b,c,d   (chromosomal occupancy block, top)
    old S4 a,b,c    ->  e,f,g     (k-mer cross-talk block, bottom)

Only the axis-label wording changed (units added in parentheses, FIG rule).

    SUPPFIG_OUT=<dir> PYTHONPATH=. python3 supp_261014/suppfigure_merged_occupancy_261103.py
"""
import json
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, genome          # noqa: E402
from _suppstyle import (BLUE, DATA_DIR, NAR_WIDTH_IN, RED,  # noqa: E402
                        apply_figure_style, figure_panel_letters, report, save)

N = json.load(open(DATA_DIR / "numbers.json"))
CORE = (1_500_000, 7_167_507)
ORIC_MB = 4.27
BINW = 50_000
TPC = {"T1": "#8e8e8e", "T2": "#5a5a5a", "T3": "#1f1f1f"}   # FIG23: timepoints are neutral greys
# 2026-11-07 補足図点検 5: パネル (c) の T1 トラックが白地に淡灰色で 6mA の疎な
# ピークを判読できなかった。3 段の下限を #b4b4b4 -> #8e8e8e に上げ、上 2 段も
# 段差を保ったまま濃くした（相対的な明度差は維持）。
XLIM0 = -1.45

P1T = pathlib.Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
                   "78_reviewer_robustness/tables")

G = genome()
GLEN = len(G)

SETS = {
    "tgg": ("pileup_GCCGGC", "21839",
            [m.start() + 4 for m in re.finditer("(?=TGGCCGGC)", G)]
            + [m.start() + 3 for m in re.finditer("(?=GCCGGCCA)", G)]),
    "c4": ("pileup_AAGCCCG_4mC", "21839",
           [m.start() + 4 for m in re.finditer("(?=AAGCCCG)", G)]
           + [m.start() + 2 for m in re.finditer("(?=CGGGCTT)", G)]),
    "a01": ("pileup_AAGCCCG_6mA", "a",
            [m.start() + i for m in re.finditer("(?=AAGCCCG)", G) for i in (0, 1)]
            + [m.start() + 5 + i for m in re.finditer("(?=CGGGCTT)", G) for i in (0, 1)]),
}
TGG_KEYS = {(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", G)} | \
           {(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", G)}


def core_frac(positions):
    a = np.asarray(positions)
    return 100 * ((a >= CORE[0]) & (a <= CORE[1])).mean()


def _pool(setname, timepoint):
    subdir, _code, _univ = SETS[setname]
    acc = {}
    for lib in LIBS[timepoint]:
        d = pd.read_csv(DATA_DIR / subdir / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, strand, nv, nm, nf in d.itertuples(index=False):
            key = (pos, strand)
            if setname == "tgg" and key not in TGG_KEYS:
                continue
            cur = acc.get(key)
            if cur is None:
                acc[key] = [nv, nm, nf]
            else:
                cur[0] += nv; cur[1] += nm; cur[2] += nf
    return acc


def called_positions(setname, timepoint):
    acc = _pool(setname, timepoint)
    return sorted(p for (p, _s), v in acc.items()
                  if (v[0] + v[2]) >= DEPTH_MIN and v[0] > 0 and v[1] / v[0] >= FRAC_MIN)


def measurable_positions(setname, timepoint):
    """Depth condition only (no modified-fraction condition)."""
    acc = _pool(setname, timepoint)
    return sorted(p for (p, _s), v in acc.items() if (v[0] + v[2]) >= DEPTH_MIN)


POS = {(s, tp): called_positions(s, tp) for s in SETS for tp in ("T1", "T2", "T3")}
for tp in ("T1", "T2", "T3"):
    assert len(POS[("tgg", tp)]) == N["TGG_called"][tp], (tp, len(POS[("tgg", tp)]))
    assert len(POS[("c4", tp)]) == N["AAG_C4_called"][tp], (tp, len(POS[("c4", tp)]))
    assert len(POS[("a01", tp)]) == N["AAG_6mA_union"][tp], (tp, len(POS[("a01", tp)]))

NULLT = {(s, tp): core_frac(measurable_positions(s, tp))
         for s in SETS for tp in ("T1", "T2", "T3")}

rows = []
for s, lbl in [("tgg", "GCCGGC hexamer 4mC (TGGCCGGC)"), ("c4", "AAGCCCG 4mC (C5)"),
               ("a01", "AAGCCCG 6mA (A0/A1)")]:
    for tp in ("T1", "T2", "T3"):
        ps = POS[(s, tp)]
        cp = core_frac(ps)
        k = int(round(len(ps) * cp / 100))
        rows.append((s, lbl, tp, len(ps), cp, NULLT[(s, tp)],
                     binomtest(k, len(ps), NULLT[(s, tp)] / 100).pvalue))
CF = pd.DataFrame(rows, columns=["key", "mark", "tp", "n_sites", "core_pct",
                                 "null_pct", "p_binom"])

# ───────────── k-mer cross-talk control (old Supplementary Figure 4) ─────────
rates = pd.read_csv(P1T / "P1_1_crosstalk_rates_C4.tsv", sep="\t")
probs = pd.read_csv(P1T / "P1_1_prob_distributions_C4.tsv", sep="\t")
aag = rates[rates.group == "AAGCCCG_6mA"].iloc[0]
non = rates[rates.group == "nonAAGCCCG_6mA_generic"].iloc[0]

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))
H = 6.23


def landscape(ax, key, title, col, headroom=1.0):
    lane_tops = []
    for i, tp in enumerate(("T1", "T2", "T3")):
        ps = np.array(POS[(key, tp)])
        h, edges = np.histogram(ps, bins=np.arange(0, GLEN + BINW, BINW))
        mx = max(h.max(), 1)
        ax.bar(edges[:-1] / 1e6, h, width=BINW / 1e6 * 0.95, bottom=i * mx * 1.15,
               color=TPC[tp], edgecolor="none")
        ax.text(-0.18, i * mx * 1.15 + mx * 0.45, f"{tp}  n={len(ps):,}", fontsize=6.0,
                ha="right", va="center", color=col)
        lane_tops.append((i + 1) * mx * 1.15)
    ax.axvspan(CORE[0] / 1e6, CORE[1] / 1e6, color="#6b6b6b", alpha=0.07, zorder=0)
    ax.axvline(ORIC_MB, color="#6b6b6b", lw=0.9, ls="--", zorder=1)
    ax.set_xlabel("Chromosome position (Mb)")
    ax.set_yticks([])
    ax.set_ylim(0, max(lane_tops) * headroom)
    ax.set_xlim(XLIM0, 8.9)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_title(title, loc="left", fontsize=6.6)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)


def corepanel(ax, keys, title):
    xs, hs, cs, xl = [], [], [], []
    for ki, k in enumerate(keys):
        for j, tp in enumerate(("T1", "T2", "T3")):
            row = CF[(CF.key == k) & (CF.tp == tp)].iloc[0]
            xs.append(ki * 3.8 + j)
            hs.append(row.core_pct)
            cs.append(TPC[tp])
            xl.append(tp)
    ax.bar(xs, hs, width=0.72, color=cs, edgecolor="none")
    for x, h in zip(xs, hs):
        ax.text(x, h + 3.4, f"{h:.1f}", ha="center", fontsize=6.0)
    nulln = float(np.mean([NULLT[(keys[0], tp)] for tp in ("T1", "T2", "T3")]))
    ax.axhline(nulln, color="#444444", lw=1.0, ls="--")
    ax.set_xticks(xs)
    ax.set_xticklabels(xl)
    ax.set_ylabel("Sites in core (%)")
    ax.set_ylim(0, 100)
    ax.set_title(title, loc="left", fontsize=6.6, pad=12)
    if len(keys) == 2:
        for gx, gl in ((1.0, "4mC (C5)"), (4.8, "6mA (A1/A2)")):
            ax.annotate(gl, xy=(gx, 0), xycoords=ax.get_xaxis_transform(),
                        xytext=(0, -15), textcoords="offset points",
                        ha="center", va="top", fontsize=6.0)


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, H))
    gs_top = fig.add_gridspec(2, 2, width_ratios=[1.95, 1.0], height_ratios=[1.0, 2.30],
                              hspace=0.42, wspace=0.42,
                              left=0.115, right=0.985, top=0.9676, bottom=0.3932)
    # 2026-11-07 補足図点検 3: 旧パネル (f) は (e) と内容が重複していた（(e) が
    # 既に near/far を AAGCCCG と generic の両方について示し、365x と 1.09x として
    # 同じ対比を与えている）。しかも (f) の generic 6mA の棒は軸下限に張り付いて
    # 見えなかった。(f) を削除し、旧 (g) を (f) に繰り上げる。
    gs_bot = fig.add_gridspec(1, 2, wspace=0.42,
                              left=0.105, right=0.985, top=0.2889, bottom=0.0562)

    ax_a = fig.add_subplot(gs_top[0, 0])
    landscape(ax_a, "tgg",
              "TGGCCGGC 4mC sites per 50 kb", RED,
              headroom=1.32)
    ax_b = fig.add_subplot(gs_top[0, 1])
    corepanel(ax_b, ["tgg"], "")

    ax_c = fig.add_subplot(gs_top[1, 0])
    lanes_c = [("c4", "AAGCCCG 4mC (C5)"), ("a01", "AAGCCCG 6mA (A1/A2)")]
    HMAX = max(np.histogram(np.array(POS[(k, tp)]),
                            bins=np.arange(0, GLEN + BINW, BINW))[0].max()
               for k, _ in lanes_c for tp in ("T1", "T2", "T3"))
    LANE = int(HMAX) + 2
    GROUP = 3 * LANE + int(HMAX) + 4
    for i, (key, nm) in enumerate(lanes_c):
        for j, tp in enumerate(("T1", "T2", "T3")):
            h, edges = np.histogram(np.array(POS[(key, tp)]),
                                    bins=np.arange(0, GLEN + BINW, BINW))
            base = i * GROUP + j * LANE
            ax_c.bar(edges[:-1] / 1e6, h, width=BINW / 1e6 * 0.95, bottom=base,
                     color=TPC[tp], edgecolor="none")
            ax_c.text(-0.18, base + LANE * 0.42, f"{tp}  n={len(POS[(key, tp)]):,}",
                      fontsize=6.0, ha="right", va="center", color=BLUE if i else RED)
        ax_c.text(0.05, i * GROUP + 2 * LANE + HMAX + 1, nm, fontsize=6.0,
                  ha="left", va="bottom", fontweight="bold", color=BLUE if i else RED)
    ax_c.axvspan(CORE[0] / 1e6, CORE[1] / 1e6, color="#6b6b6b", alpha=0.07, zorder=0)
    ax_c.axvline(ORIC_MB, color="#6b6b6b", lw=0.9, ls="--", zorder=1)
    ax_c.set_xlabel("Chromosome position (Mb)")
    ax_c.set_ylim(0, GROUP + 2 * LANE + HMAX + 1 + 1.1 * LANE)
    ax_c.set_yticks([])
    ax_c.set_xlim(XLIM0, 8.9)
    ax_c.set_xticks([0, 2, 4, 6, 8])
    ax_c.set_title("AAGCCCG sites per 50 kb",
                   loc="left", fontsize=6.6)
    for sp in ("left", "right", "top"):
        ax_c.spines[sp].set_visible(False)

    ax_d = fig.add_subplot(gs_top[1, 1].subgridspec(2, 1, height_ratios=[1.0, 1.30],
                                                    hspace=0.0)[0])
    corepanel(ax_d, ["c4", "a01"], "")

    # (e) near vs far co-call rate
    ax_e = fig.add_subplot(gs_bot[0, 0])
    x = np.arange(2)
    w = 0.35
    ax_e.bar(x - w / 2, [aag.near_rate, non.near_rate], w, label="near (\u00b13,4 bp)",
             color="#2b2b2b")
    ax_e.bar(x + w / 2, [aag.far_rate, non.far_rate], w, label="far (20\u201350 bp)",
             color="#b4b4b4")
    ax_e.set_yscale("log")
    ax_e.set_xticks(x)
    ax_e.set_xticklabels(["AAGCCCG\n6mA", "non-AAGCCCG\n6mA (generic)"])
    ax_e.set_ylabel("4mC co-call rate\n(per offset, log)")
    ax_e.set_title("4mC near a 6mA call", loc="left", fontsize=6.6)
    ax_e.legend(fontsize=6.0, frameon=False, loc="upper right")
    ax_e.set_ylim(top=ax_e.get_ylim()[1] * 3)
    ax_e.text(0 - w / 2, aag.near_rate * 1.35, f"{aag.near_over_far:.0f}\u00d7",
              ha="center", fontsize=6.5)
    ax_e.text(1 - w / 2, non.near_rate * 1.35, f"{non.near_over_far:.2f}\u00d7",
              ha="center", fontsize=6.5)

    # (f) per-read call-probability distributions（旧 g）
    ax_g = fig.add_subplot(gs_bot[0, 1])
    ax_g.plot(probs.prob_bin_center, probs.AAGCCCG_A_6mA / probs.AAGCCCG_A_6mA.sum(),
              "-o", ms=3, color=BLUE, label="A1/A2 (6mA)")
    ax_g.plot(probs.prob_bin_center, probs.AAGCCCG_C_4mC / probs.AAGCCCG_C_4mC.sum(),
              "-s", ms=3, color=RED, label="C5 (4mC)")
    ax_g.set_xlabel("Per-read probability (0\u20131)")
    ax_g.set_ylabel("Fraction of calls\n(per 0.05 bin)")
    ax_g.set_title("call confidence at AAGCCCG", loc="left", fontsize=6.6)
    ax_g.legend(fontsize=6.0, frameon=False)

    figure_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d], "abcd",
                         dx=[0.105, 0.085, 0.105, 0.085], dy=0.0092)
    figure_panel_letters(fig, [ax_e, ax_g], "ef", dx=0.085, dy=0.0118)
    return fig


if __name__ == "__main__":
    report(save(build(), "SuppFigure_merged_occupancy_261103.png",
                allow_over_data=("365\u00d7", "1.09\u00d7")))
