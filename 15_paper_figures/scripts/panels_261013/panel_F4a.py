# -*- coding: utf-8 -*-
# No `# TILES:` declaration here: that header is for scripts that write several
# tiles through a variable (panel_mine.py). This script writes one tile under a
# literal name, and declaring it twice makes check_figure_conventions FIG18 count
# two writers for one file.
"""Figure 4 panel a — site retention between timepoint pairs, by mark.

    cd 15_paper_figures && python scripts/panels_261013/panel_F4a.py

Provenance (2026-10-22)
-----------------------
This panel's drawing code was `draw4a()` inside `panels_261013/panel_mine.py`,
which wrote the tile under the pre-relettering name `panel_4a.png`; the composer
`40_main_figures_261013.py` then COPIED that file onto `panel_F4a.png`, the name
its TILES table pastes. The copy step made the tile look generator-less: nothing
in panels_261013/ wrote the name that reaches Figure 4, so `check_provenance.py`
reported a BLOCKER even though the drawing code existed. The code is moved here
unchanged (same data, same geometry, same colours) so that exactly one script
writes exactly one tile name, which is the convention panel_4c.py (-> panel_F4b)
and panel_4d.py (-> panel_F4c) already follow. `draw4a` and the composer's copy
step for F4a were removed in the same change: there is one writer, not two.

Quantity drawn
--------------
Retention = sites called modified in BOTH timepoints of a pair, divided by the
union of sites called modified in either (per-mark universe, canonical census
rule: 3 libraries pooled per timepoint, true depth (Nvalid_cov + Nfail) >= 30 in
both timepoints, pooled Nmod/Nvalid >= 0.5). Three pairs per mark
(T1->T2, T2->T3, T1->T3). Unit: per cent of the union.

Source  : data_261013/v72645d04_STAB_retention_by_mark_261013.csv
          (columns `mark`, `retention_pct`; one row per mark x timepoint pair)
Printed : the min-max range per mark, formatted from the CSV. No value typed here.
Output  : panel_F4a.png, EXACTLY 1003 x 590 px, transparent, bbox_inches=None
          (`tight` re-crops the canvas and breaks the composer's grid).
"""
import importlib
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
DATA = os.path.join(os.path.dirname(os.path.dirname(HERE)), "data_261013")
U = importlib.import_module("00_shared_utils")
U.apply_unified_style()

W_PX, H_PX = 1003, 590
OUT = os.path.join(HERE, "panel_F4a.png")

stab = pd.read_csv(os.path.join(DATA, "v72645d04_STAB_retention_by_mark_261013.csv"))

# `order` are matching keys against the input table's `mark` column, not display
# labels (the labels are given by `lab`). Changing a label in the table trips the
# assert below rather than silently dropping a series.
# 2026-10-15: the HAAGCCCG 6mA series was dropped from the main figure, not
# because it is the same quantity — the H-context subset retains 83.33-88.80 %
# against 80.00-87.62 % for all instances, which is slightly higher — but so that
# the main figure does not carry both the recognition sequence and the 5' context
# class as separate names (DECIDED.6mA.DENOMINATOR.261015). The subset values are
# quoted numerically by Results R5. The fourth series remains in the input table.
#
# 2026-10-24 (author decision D036, consistency audit CA09): the fourth series is
# restored. Results L131 attributes FOUR ranges to this panel -- including
# "83.33-88.80% for HAAGCCCG 6mA" -- and no supplementary table carries that
# series, so the figure is where it has to be drawn. This supersedes the
# 261015 omission for Figure 4a only; the 6mA DENOMINATOR decision itself
# (AAGCCCG all instances as the quoted denominator) is unchanged, which is why
# the all-instances series keeps the saturated colour and the H-context subset
# takes the pale one -- exactly as AAGCCCG 4mC is the pale partner of TGGCCGGC.
ORDER = ["TGGCCGGC 4mC (position 5)", "AAGCCCG 4mC (position 5)",
         "AAGCCCG 6mA (position 1, all instances)",
         "HAAGCCCG 6mA (position 1, H context)"]
COLS = [U.COL_4mC, mpl.colors.to_rgba(U.COL_4mC, 0.45), U.COL_6mA,
        mpl.colors.to_rgba(U.COL_6mA, 0.45)]
TCOLS = [U.COL_4mC, U.COL_4mC, U.COL_6mA, U.COL_6mA]
LAB = ["TGGCCGGC\n4mC", "AAGCCCG\n4mC", "AAGCCCG\n6mA", "HAAGCCCG\n6mA"]


def draw(fig):
    ax = fig.add_axes([0.155, 0.30, 0.815, 0.54])
    for i, (mk, c, tc) in enumerate(zip(ORDER, COLS, TCOLS)):
        sub = stab[stab["mark"] == mk]
        assert len(sub), f"保持率表に mark='{mk}' が無い。ラベル変更の取り残し"
        lo, hi = sub.retention_pct.min(), sub.retention_pct.max()
        ax.plot([i, i], [lo, hi], color=c, lw=1.1, zorder=1, solid_capstyle="round")
        for j, (_, rr) in enumerate(sub.iterrows()):
            ax.plot(i + (j - 1) * 0.17, rr.retention_pct, "o", ms=4.2, mfc=c,
                    mec="white", mew=0.6, zorder=3)
        ax.text(i, hi + 1.2, f"{lo:.1f}\u2013{hi:.1f}", ha="center", va="bottom",
                fontsize=U.FONT["annot"], color=tc)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels(LAB, fontsize=U.FONT["tick"])
    ax.set_yticks([80, 85, 90, 95, 100])
    ax.set_ylim(76, 102)
    ax.set_xlim(-0.55, len(ORDER) - 0.45)
    ax.set_ylabel("Sites retained (%)")
    # Denominators are READ FROM THE INPUT TABLE, not typed: the 2026-10-24
    # fourth series (565 H-context instances) would otherwise have left a
    # three-denominator label under a four-series panel.
    ns = " / ".join(f"{int(stab[stab['mark'] == mk].motif_universe_n.iloc[0]):,}"
                    for mk in ORDER)
    ax.set_xlabel("Three timepoint pairs per mark (T1\u2192T2, T2\u2192T3, T1\u2192T3)\n"
                  f"n = {ns} motif positions",
                  fontsize=U.FONT["annot"], color=U.COL_GRAY, labelpad=4)
    U.style_axes(ax)


if __name__ == "__main__":
    fig = plt.figure(figsize=(W_PX / 300, H_PX / 300), dpi=300)
    draw(fig)
    U.assert_no_text_collisions(fig, "Fig4 panel a")
    fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    size = Image.open(OUT).size
    assert size == (W_PX, H_PX), size
    print("wrote", os.path.basename(OUT), size)
