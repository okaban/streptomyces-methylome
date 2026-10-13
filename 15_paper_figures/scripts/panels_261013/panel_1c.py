#!/usr/bin/env python
"""Figure 1, panel c — replication gradient along the S. coelicolor M145 chromosome.

Message: T1 carries a steep replication gradient peaking at oriC; T2 and T3 are
nearly flat.

Input   : coverage table with columns timepoint, window_start_bp, window_size_bp
          (200 kb), mean_true_depth — 44 windows x 3 timepoints = 132 rows.
Output  : panel_c.png, EXACTLY 2125 x 543 px (a tile for the composer's
          12-column grid, not a standalone figure).

Run in the `epi-rna` environment.
"""

import importlib
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

# --- project style module (name starts with a digit -> import_module) ----------
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
# the tile must keep its exact pixel size: never let savefig re-crop the canvas
mpl.rcParams['savefig.bbox'] = None
# Arial has no U+2082, so the log2 subscript comes from mathtext set in Arial
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial:italic', 'mathtext.default': 'regular'})

# --- data ---------------------------------------------------------------------
CSV = "{{artifact:8e9a3013-8732-40bd-9eb0-ade6bf75e363}}"
df = pd.read_csv(CSV)

piv = df.pivot(index='window_start_bp', columns='timepoint', values='mean_true_depth')
# depth relative to each timepoint's OWN median of its 44 windows, log2
norm = np.log2(piv / piv.median())
WIN = int(df['window_size_bp'].iloc[0])                   # 200 kb
x_mb = (piv.index.values + WIN / 2) / 1e6                 # window midpoints, Mb

# --- constants ----------------------------------------------------------------
ORI = 4.27                                                # oriC, Mb
RATIOS = {'T1': '4.01', 'T2': '1.07', 'T3': '1.35'}       # oriC window / shallower terminal window
TPS = ['T1', 'T2', 'T3']
LW = {'T1': 2.0, 'T2': 1.1, 'T3': 1.1}                    # focal series dominates by weight, not hue
XL, XR = -0.20, 8.95
YB, YT = -1.28, 1.42
W_PX, H_PX = 2125, 543

# --- figure -------------------------------------------------------------------
fig = plt.figure(figsize=(W_PX / 300, H_PX / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.0565, right=0.9975, top=0.945, bottom=80 / H_PX)
U.style_axes(ax)
ax.set_xlim(XL, XR)
ax.set_ylim(YB, YT)

# reference: each timepoint's own median
ax.plot([XL, XR], [0, 0], lw=0.6, ls=(0, (1, 2.2)), color=U.COL_GRID, zorder=1)
# oriC guide, stopped below the "oriC" label
ax.plot([ORI, ORI], [YB, 1.10], lw=0.7, ls=(0, (2, 2)), color=U.COL_GRAY, zorder=1)

for i, tp in enumerate(TPS):
    ax.plot(x_mb, norm[tp].values, color=U.TIME_RAMP[i], ls=U.TIME_STYLE[i],
            lw=LW[tp], solid_capstyle='round', zorder=3 + (tp == 'T1') * 3)

ax.set_yticks([-1, 0, 1])
ax.set_yticklabels(['0.5\u00d7', '1\u00d7', '2\u00d7'])
ax.set_xticks(range(0, 9))
ax.set_ylabel(r'Depth / median ($\log_2$)', labelpad=1.5)
ax.set_xlabel('Chromosome position (Mb), 200 kb windows', labelpad=1.5)
ax.tick_params(axis='both', length=2.0, pad=1.0, width=0.8)

# Direct series labels: a colour/dash swatch carries the identity, the text is
# U.COL_DARK because U.TIME_RAMP[0] is too pale to carry 6.5 pt type.
ANCHORS = {'T1': (1.60, 0.95), 'T2': (6.90, 0.95), 'T3': (6.90, 0.62)}
labels_used = []
for i, tp in enumerate(TPS):
    xa, ya = ANCHORS[tp]
    ax.plot([xa - 0.45, xa - 0.12], [ya, ya], color=U.TIME_RAMP[i], ls=U.TIME_STYLE[i],
            lw=LW[tp], solid_capstyle='butt', zorder=4)
    txt = f'{U.TP_LABELS[i]}   ori/ter {RATIOS[tp]}'
    ax.text(xa, ya, txt, ha='left', va='center', fontsize=6.5, color=U.COL_DARK, zorder=5)
    labels_used.append(txt)

ax.text(ORI + 0.10, 1.22, 'oriC', ha='left', va='center', fontsize=6.5, color=U.COL_DARK)

TITLE = ('T1 carries a steep replication gradient peaking at oriC; '
         'T2 and T3 are nearly flat')
# placed right of the 118 x 71 px top-left area the composer keeps for the panel letter
fig.text(148 / W_PX, 1 - 19 / H_PX, TITLE, ha='left', va='center',
         fontsize=8, color=U.COL_DARK)

# pin the axis labels inside the canvas without shrinking the data box below 75 %
ax.xaxis.set_label_coords(0.5, (36 - 80) / (513.135 - 80))
ax.yaxis.set_label_coords((39.4 - 120.1) / 1999.6, 0.5)

# --- verify, then save --------------------------------------------------------
fig.canvas.draw()
U.assert_no_text_collisions(fig, 'Fig1c')
OUT = 'panel_c.png'
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
assert Image.open(OUT).size == (W_PX, H_PX), Image.open(OUT).size
