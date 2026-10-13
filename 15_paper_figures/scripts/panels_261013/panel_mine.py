# -*- coding: utf-8 -*-
"""Figure panels drawn directly for the 2026-10-13 four-figure set.

Panels 1a, 1d, 3b, 3e, 4a and 4d of the submission figures. The remaining
panels have one file each in this directory. Run from the repository with the
project environment `epi-rna`:

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

SIZES = {"1a": (1003, 614), "1d": (2053, 496), "3b": (653, 614),
         "3e": (653, 614), "4a": (1003, 590), "4d": (1003, 590)}

lib  = pd.read_csv(os.path.join(DATA, "va92d64c2_library_depth_and_sites_261013.csv"))
pos  = pd.read_csv(os.path.join(DATA, "vf7bc4cb9_F3b_positions.csv"))
occ  = pd.read_csv(os.path.join(DATA, "vf8771452_SuppTable_motif_occupancy.csv"))
stab = pd.read_csv(os.path.join(DATA, "v72645d04_STAB_retention_by_mark_261013.csv"))
e4   = pd.read_csv(os.path.join(DATA, "v11136318_F4e_enzyme.csv"))


def _ci(row, tp):
    lo, hi = str(row[f"{tp}_occ_CI95"]).replace("\u2013", "-").split("-")
    return float(lo), float(hi)


r4 = occ[occ.motif_mark.str.startswith("TGGCCGGC")].iloc[0]
r6 = occ[occ.motif_mark.str.startswith("HAAGCCCG")].iloc[0]
o4 = [float(r4[f"{t}_occupancy_pct"]) for t in ("T1", "T2", "T3")]
o6 = [float(r6[f"{t}_occupancy_pct"]) for t in ("T1", "T2", "T3")]
e4ci = np.array([_ci(r4, t) for t in ("T1", "T2", "T3")]).T
e6ci = np.array([_ci(r6, t) for t in ("T1", "T2", "T3")]).T
rat = [round(100 * b / a, 1) for a, b in zip(o4, o6)]


def tile(px, fn, draw):
    w, h = px
    fig = plt.figure(figsize=(w / 300, h / 300), dpi=300)
    draw(fig)
    U.assert_no_text_collisions(fig, fn)
    fig.savefig(fn, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    return fn


def draw1a(fig):
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.set_axis_off()
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 100); ax.set_ylim(0, 104)
    shades = list(U.TIME_RAMP[:3]); XS = [52.0, 70.0, 88.0]; LBL = 2.0
    for name, y in [("Culture time", 82), ("Biological replicates", 62),
                    ("Nanopore sequencing", 42), ("Strand-specific RNA-seq", 22)]:
        ax.text(LBL, y, name, fontsize=U.FONT["tick"], color=U.COL_DARK, ha="left", va="center")
    ax.annotate("", xy=(97, 82), xytext=(46, 82),
                arrowprops=dict(arrowstyle="-|>", color=U.COL_DARK, lw=0.9))
    for lab, x, c in zip(["T1  12 h", "T2  24 h", "T3  50 h"], XS, shades):
        ax.plot([x], [82], "o", ms=5.0, mfc=c, mec=U.COL_DARK, mew=0.7, zorder=3)
        ax.text(x, 88.5, lab, ha="center", fontsize=U.FONT["tick"], color=U.COL_DARK)
    for x, c in zip(XS, shades):
        ax.plot([x, x], [79.5, 65.5], color=U.COL_GRAY, lw=0.6, zorder=1)
        for dx in (-4.0, 0, 4.0):
            ax.plot([x + dx], [62], "o", ms=4.2, mfc=c, mec="white", mew=0.6, zorder=3)
            ax.plot([x, x + dx], [65.5, 63.6], color=U.COL_GRAY, lw=0.5, zorder=1)
            ax.plot([x + dx, x + dx], [60.2, 44.4], color=U.COL_GRAY, lw=0.5, zorder=1)
            ax.plot([x + dx], [42], "s", ms=4.6, mfc=U.COL_4mC, mec="white", mew=0.6, zorder=3)
            ax.plot([x + dx, x + dx], [40.0, 24.0], color=U.COL_GRAY, lw=0.5, zorder=1)
            ax.plot([x + dx], [22], "s", ms=4.6, mfc=U.COL_6mA, mec="white", mew=0.6, zorder=3)
    ax.text(LBL, 36.0, "9 libraries  →  4mC / 6mA occupancy per site",
            fontsize=U.FONT["annot"], color=U.COL_4mC, ha="left", va="center")
    ax.text(LBL, 16.0, "9 libraries  →  transcript abundance (TPM)",
            fontsize=U.FONT["annot"], color=U.COL_6mA, ha="left", va="center")
    ax.text(LBL, 5.0, "three biological replicates per timepoint, each split for both assays",
            fontsize=U.FONT["annot"], color=U.COL_GRAY, ha="left", va="center")
    fig.text(0.012, 0.995, "Nine samples, each sequenced two ways",
             fontsize=TTL, color=U.COL_DARK, ha="left", va="top")


def draw1d(fig):
    ax = fig.add_axes([0.085, 0.245, 0.895, 0.545])
    cols = {"T1": U.TIME_RAMP[0], "T2": U.TIME_RAMP[1], "T3": U.TIME_RAMP[2]}
    for i, (_, r) in enumerate(lib.iterrows()):
        out = r.library == "3-2"
        ax.bar(i, r.median_true_depth, width=0.64, color=cols[r.timepoint],
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
    ax.set_xlabel("sequencing library;  number above bar = TGGCCGGC 4mC sites contributed",
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
    ax.text(0.0, 56, "6mA\n46.6", color=U.COL_6mA, fontsize=U.FONT["annot"], ha="center", va="bottom")
    ax.text(3.0, 84, "4mC\n74.4", color=U.COL_4mC, fontsize=U.FONT["annot"], ha="center", va="bottom")
    ax.set_xticks(range(4)); ax.set_xticklabels(labs, fontsize=U.FONT["tick"])
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(0, 104); ax.set_yticks([0,20,40,60,80,100])
    ax.set_ylabel("Occupancy (%)")
    ax.set_xlabel("position in AAGCCCG\n● T1  ■ T2  ▲ T3;  bar = mean", fontsize=U.FONT["annot"],
                  color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)
    ax.set_title("A 5′ G blocks 6mA at A1", fontsize=TTL, loc="left", pad=6)


def draw3e(fig):
    ax = fig.add_axes([0.235, 0.26, 0.58, 0.52])
    xx = np.arange(3)
    ax.errorbar(xx-0.14, o4, yerr=[np.array(o4)-e4ci[0], e4ci[1]-np.array(o4)], fmt="o", ms=5,
                color=U.COL_4mC, capsize=2.5, lw=1.2, zorder=3)
    ax.errorbar(xx+0.14, o6, yerr=[np.array(o6)-e6ci[0], e6ci[1]-np.array(o6)], fmt="s", ms=5,
                color=U.COL_6mA, capsize=2.5, lw=1.2, zorder=3)
    for j in range(3):
        ax.text(j, 104, f"{rat[j]:.0f}%", ha="center", fontsize=U.FONT["annot"], color=U.COL_DARK)
    ax.text(2.45, o4[2]+2, "4mC\nTGGCCGGC\nn = 2,093", color=U.COL_4mC, fontsize=U.FONT["annot"], va="center")
    ax.text(2.45, o6[2]-2, "6mA\nHAAGCCCG\nn = 565", color=U.COL_6mA, fontsize=U.FONT["annot"], va="center")
    ax.set_xticks(xx); ax.set_xticklabels(["T1","T2","T3"]); ax.set_xlim(-0.5,4.05)
    ax.set_yticks([40,60,80,100]); ax.set_ylim(28,112); ax.set_ylabel("Occupancy (%)")
    ax.set_xlabel("numbers above: 6mA as % of 4mC\nbars are Wilson 95% intervals", fontsize=U.FONT["annot"],
                  color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)
    ax.set_title("6mA runs at about half\nthe 4mC occupancy", fontsize=TTL, loc="left", pad=6)


def draw4a(fig):
    ax = fig.add_axes([0.155, 0.30, 0.815, 0.54])
    order = ["TGGCCGGC 4mC (position 5)","AAGCCCG 4mC (position 5)",
             "AAGCCCG 6mA (A1, all instances)","HAAGCCCG 6mA (A1, H context)"]
    cols  = [U.COL_4mC, mpl.colors.to_rgba(U.COL_4mC,0.45), U.COL_6mA, mpl.colors.to_rgba(U.COL_6mA,0.45)]
    tcols = [U.COL_4mC, U.COL_4mC, U.COL_6mA, U.COL_6mA]
    lab   = ["TGGCCGGC\n4mC","AAGCCCG\n4mC","AAGCCCG\n6mA","HAAGCCCG\n6mA"]
    for i,(mk,c,tc) in enumerate(zip(order, cols, tcols)):
        sub = stab[stab["mark"]==mk]
        lo, hi = sub.retention_pct.min(), sub.retention_pct.max()
        ax.plot([i,i],[lo,hi], color=c, lw=1.1, zorder=1, solid_capstyle="round")
        for j,(_,rr) in enumerate(sub.iterrows()):
            ax.plot(i+(j-1)*0.17, rr.retention_pct, "o", ms=4.2, mfc=c, mec="white", mew=0.6, zorder=3)
        ax.text(i, hi+1.2, f"{lo:.1f}\u2013{hi:.1f}", ha="center", va="bottom", fontsize=U.FONT["annot"], color=tc)
    ax.set_xticks(range(4)); ax.set_xticklabels(lab, fontsize=U.FONT["tick"])
    ax.set_yticks([80,85,90,95,100]); ax.set_ylim(76,102); ax.set_xlim(-0.55,3.55)
    ax.set_ylabel("Sites retained (%)")
    ax.set_xlabel("three timepoint pairs per mark (T1\u2192T2, T2\u2192T3, T1\u2192T3)\nn = 2,093 / 1,334 / 1,334 / 565 motif positions",
                  fontsize=U.FONT["annot"], color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)
    ax.set_title("4mC sites are held across time; 6mA sites less so", fontsize=TTL, loc="left", pad=6)


def draw4d(fig):
    ax = fig.add_axes([0.135, 0.255, 0.70, 0.52]); ax2 = ax.twinx()
    x = np.arange(3)
    ax.bar(x - 0.17, e4.SC_RS17645_TPM, width=0.30, color=U.COL_GRAY, zorder=2)
    ax2.plot(x + 0.17, e4.occ_6mA_HAAGCCCG_pct, "o-", color=U.COL_6mA, lw=1.4, ms=4.6, zorder=3)
    ax.set_ylim(0, 12); ax.set_yticks([0, 4, 8, 12]); ax2.set_ylim(0, 60); ax2.set_yticks([0, 20, 40, 60])
    ax.set_xticks(x); ax.set_xticklabels(["T1", "T2", "T3"]); ax.set_xlim(-0.55, 2.75)
    ax.set_ylabel("Transcript level (TPM)", color=U.COL_GRAY)
    ax2.set_ylabel("6mA occupancy (%)", color=U.COL_6mA)
    ax2.tick_params(axis="y", colors=U.COL_6mA); ax.tick_params(axis="y", colors=U.COL_GRAY)
    ax.text(1.0, 2.9, "transcript falls 6.4×", ha="center", fontsize=U.FONT["annot"], color=U.COL_GRAY)
    ax2.text(2.32, 46.0, "occupancy\nholds at\n45.8–47.8%", fontsize=U.FONT["annot"],
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
            t_.set_text("occupancy holds\nat 45.8\u201347.8%")
            t_.set_position((0.02, 56.5)); t_.set_va("top"); t_.set_ha("left")


if __name__ == "__main__":
    out = {"1a": draw1a, "1d": draw1d, "3b": draw3b, "3e": draw3e, "4a": draw4a, "4d": draw4d}
    for key, fn in out.items():
        p = tile(SIZES[key], os.path.join(HERE, f"panel_{key}.png"), fn)
        print("wrote", os.path.basename(p))
