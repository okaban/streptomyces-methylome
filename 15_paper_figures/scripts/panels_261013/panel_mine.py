# TILES: panel_1a.png panel_1d.png panel_3b.png panel_3e.png panel_4d.png
# [FIG18] 複数タイルを変数経由で書き出すため、生成物を明示宣言する。
# -*- coding: utf-8 -*-
"""Figure panels drawn directly for the 2026-10-13 four-figure set.

Panels 1a, 1d, 3b, 3e and 4d of the submission figures. The remaining
panels have one file each in this directory. Run from the repository with the
project environment `epi-rna`:

2026-10-22: `draw4a` moved out to `panel_F4a.py`, which writes the tile name the
composer actually pastes (`panel_F4a.png`) instead of the pre-relettering name
`panel_4a.png` that the composer then copied. One writer per tile name.

    python panels_261013/panel_mine.py

Every constant comes from 00_shared_utils (project figure rules forbid
hardcoded palette hexes). Output sizes are the composed-grid tile sizes at
174 mm full width, which is the ceiling check_figures.py enforces.
"""
import sys, os, importlib
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data_261013")
U = importlib.import_module("00_shared_utils")
U.apply_unified_style()
TTL = U.FONT["panel_title"]
MK = ["o", "s", "^"]

SIZES = {"1a": (2053, 543), "1d": (2053, 496), "3b": (653, 614),
         "3e": (653, 614), "4d": (1003, 590)}

lib  = pd.read_csv(os.path.join(DATA, "va92d64c2_library_depth_and_sites_261013.csv"))
pos  = pd.read_csv(os.path.join(DATA, "vf7bc4cb9_F3b_positions.csv"))
occ  = pd.read_csv(os.path.join(DATA, "vf8771452_SuppTable_motif_occupancy.csv"))
e4   = pd.read_csv(os.path.join(DATA, "v11136318_F4e_enzyme.csv"))


def _ci(row, tp):
    lo, hi = str(row[f"{tp}_occ_CI95"]).replace("\u2013", "-").split("-")
    return float(lo), float(hi)


# Figure 3e は「同じ 1,334 インスタンス上での 4mC 対 6mA」を示す。
# 2026-10-15 著者決定（DECIDED.6mA.DENOMINATOR.261015）により 6mA の分母は
# HAAGCCCG 565 から AAGCCCG 全インスタンスに変更した。比の注記も同一分母どうしになる。
# motif_mark は CSV の照合キーであって表示ラベルではない（表示は下の文字列で与える）。
rT = occ[occ.motif_mark.str.startswith("TGGCCGGC")].iloc[0]
r4 = occ[occ.motif_mark == "AAGCCCG 4mC (C4)"].iloc[0]
r6 = occ[occ.motif_mark == "AAGCCCG 6mA (A0) — all instances"].iloc[0]
oT = [float(rT[f"{t}_occupancy_pct"]) for t in ("T1", "T2", "T3")]
o4 = [float(r4[f"{t}_occupancy_pct"]) for t in ("T1", "T2", "T3")]
o6 = [float(r6[f"{t}_occupancy_pct"]) for t in ("T1", "T2", "T3")]
eTci = np.array([_ci(rT, t) for t in ("T1", "T2", "T3")]).T
e4ci = np.array([_ci(r4, t) for t in ("T1", "T2", "T3")]).T
e6ci = np.array([_ci(r6, t) for t in ("T1", "T2", "T3")]).T
rat = [round(100 * b / a, 1) for a, b in zip(o4, o6)]   # 6mA / 4mC、どちらも分母 1,334


def tile(px, fn, draw):
    w, h = px
    fig = plt.figure(figsize=(w / 300, h / 300), dpi=300)
    draw(fig)
    U.assert_no_text_collisions(fig, fn)
    fig.savefig(fn, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    return fn


def _x1a(h):
    """Culture hours -> drawing units on the 1a timeline (0 h at 18, 62 h at 97).

    2026-10-24, D058 (author approved): the span was 0-56 h with ticks every 12 h,
    so the last tick (48) sat just left of the T3 node (50 h) with no number above
    it, and T3 could be read as 48 h. The span now runs to 62 h and the ticks to
    60, putting a numbered tick on each side of T3.
    """
    return 18.0 + (h / 62.0) * 79.0


def draw1a(fig):
    """Experimental design.

    2026-10-24, D044 (author comment #25): candidate A of
    figure_readability_261024.md §1.2, approved. Three changes against the
    2026-10-13 panel, each tied to a finding in that table:
      F1  the 27-glyph replicate grid said "3 replicates x 3 timepoints" a
          fourth time -> one branch per timepoint, multiplied by "x 3" (9 glyphs).
      F2  the caption sentence at y=5 was legend prose -> removed.
      F4  nodes were equally spaced (46/64/82) while sampling is 12/24/50 h
          -> nodes at their real times under an hour axis.
      F5  the dead gutter at x=30..36 -> content starts at x=18.
    No value or experimental condition is added; the hours, the replicate count
    and the two library counts are the ones the 2026-10-13 panel already carried.

    The hour unit rides on the row label instead of standing as its own text.
    The candidate render put "h after inoculation" at the right end of the tick
    row, abutting the 48 h numeral; moving it to the left of that row then put
    it under the composer's bold panel letter, which the collision checker
    cannot see because panel a reaches the composer as a raster tile. The row
    label is the one position in this panel that is clear of both. The zero
    point is named on the panel 1b axis ("Time after inoculation (h)").
    Tick marks and numerals go above the arrow: below it is where the branches
    descend, and a numeral there is crossed by a connector.
    """
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.set_axis_off()
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 100); ax.set_ylim(0, 104)
    shades = list(U.TIME_RAMP[:3]); LBL = 1.5
    HOURS = [12, 24, 50]; TAGS = ["T1", "T2", "T3"]
    for name, y in [("Culture time (h)", 86), ("Biological replicates", 63),
                    ("Nanopore sequencing", 41), ("Strand-specific RNA-seq", 18)]:
        ax.text(LBL, y, name, fontsize=U.FONT["tick"], color=U.COL_DARK, ha="left", va="center")
    ax.annotate("", xy=(99.5, 86), xytext=(_x1a(0), 86),
                arrowprops=dict(arrowstyle="-|>", color=U.COL_DARK, lw=0.9))
    for h in (0, 12, 24, 36, 48, 60):
        ax.plot([_x1a(h), _x1a(h)], [86, 87.6], color=U.COL_DARK, lw=0.7, zorder=2)
        ax.text(_x1a(h), 96.0, f"{h}", ha="center", va="center",
                fontsize=U.FONT["annot"], color=U.COL_GRAY)
    for tag, h, c in zip(TAGS, HOURS, shades):
        x = _x1a(h)
        ax.plot([x], [86], "o", ms=5.4, mfc=c, mec=U.COL_DARK, mew=0.7, zorder=4)
        ax.text(x, 90.9, tag, ha="center", va="center",
                fontsize=U.FONT["tick"], color=U.COL_DARK)
        ax.plot([x, x], [83.0, 67.0], color=U.COL_GRAY, lw=0.6, zorder=1)
        ax.plot([x], [63], "o", ms=4.6, mfc=c, mec="white", mew=0.6, zorder=4)
        ax.text(x + 2.4, 63, "\u00d7 3", ha="left", va="center",
                fontsize=U.FONT["annot"], color=U.COL_DARK)
        ax.plot([x, x], [60.5, 44.0], color=U.COL_GRAY, lw=0.6, zorder=1)
        ax.plot([x], [41], "s", ms=4.8, mfc=U.COL_4mC, mec="white", mew=0.6, zorder=4)
        ax.plot([x, x], [38.5, 21.0], color=U.COL_GRAY, lw=0.6, zorder=1)
        ax.plot([x], [18], "s", ms=4.8, mfc=U.COL_6mA, mec="white", mew=0.6, zorder=4)
    ax.text(LBL, 33.0, "9 libraries  →  4mC / 6mA occupancy per site",
            fontsize=U.FONT["annot"], color=U.COL_4mC, ha="left", va="center")
    ax.text(LBL, 10.0, "9 libraries  →  transcript abundance (TPM)",
            fontsize=U.FONT["annot"], color=U.COL_6mA, ha="left", va="center")


def draw1d(fig):
    ax = fig.add_axes([0.085, 0.245, 0.895, 0.545])
    cols = {"T1": U.TIME_RAMP[0], "T2": U.TIME_RAMP[1], "T3": U.TIME_RAMP[2]}
    for i, (_, r) in enumerate(lib.iterrows()):
        out = r.library == "3-2"
        ax.bar(i, r.median_true_depth, width=0.64, color=U.COL_BAR_MID,
               edgecolor=U.COL_4mC if out else "none", lw=1.0 if out else 0, zorder=2)
        ax.text(i, r.median_true_depth + 1.8, f"{int(r.tg_sites):,}", ha="center",
                fontsize=U.FONT["annot"], color=U.COL_4mC if out else U.COL_DARK)
    ax.axhline(30, color=U.COL_GRAY, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.text(-0.45, 31.5, "site-calling minimum 30×", fontsize=U.FONT["annot"],
            color=U.COL_GRAY, ha="left", va="bottom")
    ax.text(6, 24.5, "depth outlier", fontsize=U.FONT["annot"], color=U.COL_4mC, ha="center")
    for j, (tp, x) in enumerate((("T1 (12 h)", 1), ("T2 (24 h)", 4), ("T3 (50 h)", 7))):
        ax.text(x, 64.5, tp, ha="center", fontsize=U.FONT["tick"], color=U.TIME_RAMP[j])
    ax.set_xticks(range(len(lib))); ax.set_xticklabels(lib.library, fontsize=U.FONT["tick"])
    ax.set_ylim(0, 70); ax.set_yticks([0, 20, 40, 60]); ax.set_xlim(-0.6, 8.6)
    ax.set_ylabel("Median true depth (×)")
    ax.set_xlabel("Sequencing library;  number above bar = TGGCCGGC 4mC sites contributed",
                  fontsize=U.FONT["annot"], color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)
    ax.set_title("Library 3-2 is a depth outlier at 15× and yields only 36 sites",
                 fontsize=TTL, loc="left", pad=6)
    ax = fig.axes[0]
    for t_ in list(ax.texts):
        if t_.get_text() == "depth outlier": t_.remove()
        elif t_.get_text().startswith("site-calling"):
            t_.set_position((-0.45, 46.0)); t_.set_ha("left"); t_.set_va("bottom")


def draw3b(fig):
    ax = fig.add_axes([0.215, 0.285, 0.735, 0.525])
    labs = ["A1\n5′ A/C/T", "A1\n5′ G", "A2", "pos. 5"]
    cols = [U.COL_6mA, mpl.colors.to_rgba(U.COL_6mA,0.5), mpl.colors.to_rgba(U.COL_6mA,0.5), U.COL_4mC]
    for i,(_,r) in enumerate(pos.iterrows()):
        vals = [r.T1, r.T2, r.T3]
        for j,v in enumerate(vals):
            ax.plot(i+(j-1)*0.17, v, MK[j], ms=4.0, mfc=cols[i], mec="white", mew=0.5, zorder=3)
        ax.plot([i-0.30,i+0.30],[np.mean(vals)]*2,
                color=cols[i] if isinstance(cols[i],str) else U.COL_6mA, lw=1.4, zorder=2, solid_capstyle="butt")
    ax.text(0.0, 56, "6mA\n46.2", color=U.COL_6mA, fontsize=U.FONT["annot"], ha="center", va="bottom")
    ax.text(3.0, 84, "4mC\n74.3", color=U.COL_4mC, fontsize=U.FONT["annot"], ha="center", va="bottom")
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=U.FONT["tick"])
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(0, 104); ax.set_yticks([0,20,40,60,80,100])
    ax.set_ylabel("Occupancy (%)")
    ax.set_xlabel("Position in AAGCCCG\n● T1  ■ T2  ▲ T3;  bar = mean", fontsize=U.FONT["annot"],
                  color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)


def draw3e(fig):
    ax = fig.add_axes([0.235, 0.26, 0.58, 0.52])
    xx = np.arange(3)
    ax.errorbar(xx-0.20, oT, yerr=[np.array(oT)-eTci[0], eTci[1]-np.array(oT)], fmt="o", ms=5,
                color=U.COL_4mC, capsize=2.5, lw=1.2, zorder=3, mfc="white", mew=1.3)
    ax.errorbar(xx, o4, yerr=[np.array(o4)-e4ci[0], e4ci[1]-np.array(o4)], fmt="o", ms=5,
                color=U.COL_4mC, capsize=2.5, lw=1.2, zorder=3)
    ax.errorbar(xx+0.20, o6, yerr=[np.array(o6)-e6ci[0], e6ci[1]-np.array(o6)], fmt="s", ms=5,
                color=U.COL_6mA, capsize=2.5, lw=1.2, zorder=3)
    ax.text(2.42, oT[2], "4mC TGGCCGGC", color=U.COL_4mC,
            fontsize=U.FONT["annot"], va="center")
    ax.text(2.42, o4[2], "4mC AAGCCCG", color=U.COL_4mC,
            fontsize=U.FONT["annot"], va="center")
    ax.text(2.42, o6[2], "6mA AAGCCCG", color=U.COL_6mA,
            fontsize=U.FONT["annot"], va="center")
    ax.set_xticks(xx); ax.set_xticklabels(["T1","T2","T3"]); ax.set_xlim(-0.5,4.05)
    ax.set_yticks([20,40,60,80,100]); ax.set_ylim(10,104); ax.set_ylabel("Occupancy (%)")
    ax.set_xlabel("n = 2,093 and 1,334; Wilson 95% CI",
                  fontsize=U.FONT["annot"], color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)


def draw4d(fig):
    ax = fig.add_axes([0.135, 0.255, 0.70, 0.52]); ax2 = ax.twinx()
    x = np.arange(3)
    ax.bar(x - 0.17, e4.SC_RS17645_TPM, width=0.30, color=U.COL_GRAY, zorder=2)
    ax2.plot(x + 0.17, e4.occ_6mA_AAGCCCG_pct, "o-", color=U.COL_6mA, lw=1.4, ms=4.6, zorder=3)
    ax.set_ylim(0, 12); ax.set_yticks([0, 4, 8, 12]); ax2.set_ylim(0, 40); ax2.set_yticks([0, 10, 20, 30, 40])
    ax.set_xticks(x); ax.set_xticklabels(["T1", "T2", "T3"]); ax.set_xlim(-0.55, 2.75)
    ax.set_ylabel("Transcript level (TPM)", color=U.COL_GRAY)
    ax2.set_ylabel("6mA occupancy (%)", color=U.COL_6mA)
    ax2.tick_params(axis="y", colors=U.COL_6mA); ax.tick_params(axis="y", colors=U.COL_GRAY)
    ax.text(1.0, 2.9, "transcript falls 6.4×", ha="center", fontsize=U.FONT["annot"], color=U.COL_GRAY)
    ax2.text(2.32, 25.2, "occupancy\nholds at\n24.2–26.2%", fontsize=U.FONT["annot"],
             color=U.COL_6mA, va="center")
    U.style_axes(ax)
    for s in ("top", "left"): ax2.spines[s].set_visible(False)
    ax2.spines["right"].set_color(U.COL_6mA); ax2.spines["right"].set_linewidth(0.8)
    ax.set_xlabel("n = 3 timepoints; no correlation claimed", fontsize=U.FONT["annot"],
                  color=U.COL_GRAY, labelpad=4)
    ax.set_title("$\\it{SC\\_RS17645}$ is expressed at every\ntimepoint; occupancy does not follow it",
                 fontsize=TTL, loc="left", pad=6)
    ax2 = fig.axes[1]
    for t_ in list(ax2.texts):
        if t_.get_text().startswith("occupancy"):
            t_.set_text("occupancy holds\nat 24.2\u201326.2%")
            t_.set_position((0.02, 37.8)); t_.set_va("top"); t_.set_ha("left")


if __name__ == "__main__":
    out = {"1a": draw1a, "1d": draw1d, "3b": draw3b, "3e": draw3e, "4d": draw4d}
    for key, fn in out.items():
        p = tile(SIZES[key], os.path.join(HERE, f"panel_{key}.png"), fn)
        print("wrote", os.path.basename(p))
