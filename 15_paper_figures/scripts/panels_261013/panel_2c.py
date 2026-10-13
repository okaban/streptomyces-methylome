#!/usr/bin/env python3
"""Figure 2, panel c — site counts passing the calling rule, by series and timepoint.

Message: no 5mC site passes the rule that calls 1,886 4mC sites.

Chart family: horizontal grouped bars, one group per series, three rows per group
(T1 -> T3, top to bottom). The motif names are 8 characters and do not fit as
vertical-bar x tick labels at any size on the project font ladder in a 677 px-wide
tile, so each group is named by direct text above its bars.

Colour encodes the CHEMISTRY only (project figure rule 2): U.COL_4mC for both 4mC
series, U.COL_6mA for the 6mA series, U.COL_GRAY for 5mC / not detected. Timepoint
is a white-blend tint of the same constant, never a new hex.

Output: panel_c.png, exactly 677 x 661 px, transparent, no tight bbox — this is a
TILE for the composer, not a standalone figure.
"""

import importlib
import sys

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

CSV = "{{artifact:22031ca8-89ff-4199-9bcd-c9cf0fb149b1}}"
OUT = "panel_c.png"
PX = (677, 661)          # tile size handed down by the composer
DPI = 300


def shade(col, f):
    """Blend a palette constant toward white by fraction f (1.0 = full colour)."""
    r, g, b = mcolors.to_rgb(col)
    return (r + (1 - r) * (1 - f), g + (1 - g) * (1 - f), b + (1 - b) * (1 - f))


# ---------------------------------------------------------------- data
df = pd.read_csv(CSV)
counts = df.set_index('series')[['T1', 'T2', 'T3']].values

LBL = ["TGGCCGGC 4mC", "AAGCCCG 4mC", "HAAGCCCG 6mA", "5mC, any motif position"]
COLS = [U.COL_4mC, U.COL_4mC, U.COL_6mA, U.COL_GRAY]
F = [0.45, 0.72, 1.0]    # T1 -> T3 tint of the series' chemistry colour
XMAX = 2300

# ---------------------------------------------------------------- layout
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # the style helper sets 'tight', which resizes the canvas

a, p, s = 0.40, 0.26, 0.55            # label-to-first-bar, bar pitch, group gap
stride = a + 2 * p + s
Lrow = [i * stride for i in range(4)]
Brow = [[Lrow[i] + a + j * p for j in range(3)] for i in range(4)]
YTOP, YBOT = -0.22, Brow[3][2] + 0.22

fig = plt.figure(figsize=(PX[0] / DPI, PX[1] / DPI), dpi=DPI)
ax = fig.add_subplot(111)
U.style_axes(ax)
fig.subplots_adjust(left=0.040, right=0.985, bottom=0.215, top=0.860)

# ---------------------------------------------------------------- bars
for i in range(4):
    for j in range(3):
        v, y, c = counts[i, j], Brow[i][j], shade(COLS[i], F[j])
        if v == 0:
            # zero category: a visible baseline dot, never an empty slot
            ax.plot([0], [y], marker='o', ms=2.6, mfc=c, mec=c, clip_on=False, zorder=4)
        else:
            ax.barh(y, v, height=0.24, color=c, linewidth=0, zorder=3)
    ax.text(28, Lrow[i], LBL[i], ha='left', va='center', fontsize=6.5, color=U.COL_DARK)

# timepoint shade key, side by side on the first label row (identity, stated once)
for j, (xv, lab) in enumerate(zip([1520, 1779, 2038], ["T1", "T2", "T3"])):
    ax.plot([xv], [Lrow[0]], marker='s', ms=4.2, mfc=shade(U.COL_4mC, F[j]),
            mec=shade(U.COL_4mC, F[j]), zorder=4)
    ax.annotate(lab, xy=(xv, Lrow[0]), xytext=(5, 0), textcoords='offset points',
                ha='left', va='center', fontsize=6.5, color=U.COL_DARK)

# ---------------------------------------------------------------- axes
ax.set_xlim(0, XMAX)
ax.set_ylim(YBOT, YTOP)
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.set_xticks([0, 500, 1000, 1500, 2000])
ax.set_xlabel("Sites (depth \u2265 30, level \u2265 50%)", labelpad=2)
ax.set_title("No 5mC site passes the rule\nthat calls 1,886 4mC sites",
             loc='left', x=0.153, pad=4)   # x offset keeps the composer's letter corner clear

# ---------------------------------------------------------------- annotations
ax.annotate("1,886", xy=(counts[0, 0], Brow[0][0]), xytext=(3, 0),
            textcoords='offset points', ha='left', va='center',
            fontsize=6.5, color=U.COL_4mC)
ax.annotate("n.d.", xy=(0, Brow[3][1]), xytext=(7, 0), textcoords='offset points',
            ha='left', va='center', fontsize=6.5, color=U.COL_GRAY)
fig.text(0.5, 0.008, "modkit: C 0.9336\u22120.9375, A 0.6230\u22120.6387",
         ha='center', va='bottom', fontsize=6.5, color=U.COL_GRAY)

# ---------------------------------------------------------------- save + verify
U.assert_no_text_collisions(fig, "panel_c")
fig.savefig(OUT, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)

h, w = plt.imread(OUT).shape[:2]
assert (w, h) == PX, f"tile is {(w, h)}, expected {PX}"
print(f"{OUT}: {w}x{h}")
