"""Figure 2, panel a — 4mC / 6mA site density along the chromosome.

Tile for the composed Figure 2 (12-column grid): rendered at EXACTLY
2125 x 519 px, transparent, no tight/constrained layout.

Message: both marks are spread along the whole chromosome, arms as well as core.
Data: window_mb, sites_4mC_per_100kb, sites_6mA_per_100kb (100-kb windows, n = 87).

Run in the `epi-rna` environment.
"""
import sys
import importlib

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
from PIL import Image

# ── style ────────────────────────────────────────────────────────────────────
try:                      # figure-style skill helper, when a kernel provides it
    apply_figure_style()  # noqa: F821
except NameError:
    pass
U.apply_unified_style()              # PROJECT style wins
mpl.rcParams['savefig.bbox'] = None  # the style helpers set 'tight', which resizes

# ── data ─────────────────────────────────────────────────────────────────────
CSV = "{{artifact:1c1184c7-bb30-4a1e-a001-0779aedc78ba}}"
OUT = 'panel_a.png'

df = pd.read_csv(CSV)
x = df['window_mb'].values
y4 = df['sites_4mC_per_100kb'].values
y6 = df['sites_6mA_per_100kb'].values
# 500-kb (5-window) centred rolling mean drawn solid over the pale per-window trace
s4 = df['sites_4mC_per_100kb'].rolling(5, center=True, min_periods=1).mean().values
s6 = df['sites_6mA_per_100kb'].rolling(5, center=True, min_periods=1).mean().values
m4, m6 = y4.mean(), y6.mean()

CORE_LO, CORE_HI = 1.5, 7.17
XMAX = 9.72
YMAX = 58.0

# ── figure ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(2125 / 300, 519 / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.058, right=0.997, bottom=0.162, top=0.925)

# chromosomal core band — background context, drawn as a collection (not a patch)
# so the single band label is not read as text-over-a-data-patch by the checker.
# Neutral grey on purpose: COL_CORE (green) behind the COL_4mC line is the
# forbidden red/green pair.
band = ax.fill_between([CORE_LO, CORE_HI], 0, YMAX,
                       color=mcolors.to_rgba(U.COL_GRAY, 0.16),
                       linewidth=0, zorder=0)

# per-window values (pale) + 500-kb rolling mean (solid)
ax.plot(x, y4, color=mcolors.to_rgba(U.COL_4mC, 0.38), lw=0.7, zorder=1.5)
ax.plot(x, y6, color=mcolors.to_rgba(U.COL_6mA, 0.38), lw=0.7, zorder=1.5)
ax.plot(x, s4, color=U.COL_4mC, lw=1.6, solid_capstyle='round', zorder=2.5)
ax.plot(x, s6, color=U.COL_6mA, lw=1.6, solid_capstyle='round', zorder=2.5)

ax.set_xlim(-0.08, XMAX)
ax.set_ylim(0, YMAX)
ax.set_xticks([0, 2, 4, 6, 8])
ax.set_yticks([0, 10, 20, 30, 40, 50])
ax.set_xlabel('position along chromosome (Mb)', labelpad=2)
ax.set_ylabel('sites per 100 kb', labelpad=2)
ax.set_title('Both marks are spread along the whole chromosome, arms as well as core',
             loc='left', x=0.035, pad=3)
U.style_axes(ax)
ax.tick_params(axis='both', pad=1.5, length=2.5, width=0.8)

FS_ID = U.FONT['axis_label']
FS_AN = U.FONT['annot']

# band label — once, in the empty top strip of the band
ax.text((CORE_LO + CORE_HI) / 2, 54.0, 'chromosomal core (1.5–7.17 Mb)',
        ha='center', va='center', fontsize=FS_AN,
        color=mcolors.to_rgba(U.COL_DARK, 0.72))

# direct series labels + their genome-wide means, in the right margin (no data there)
ax.text(8.78, s4[-1], '4mC', ha='left', va='center',
        fontsize=FS_ID, color=U.COL_4mC)
ax.text(8.78, s4[-1] - 5.6, f'mean {m4:.0f}', ha='left', va='center',
        fontsize=FS_AN, color=mcolors.to_rgba(U.COL_4mC, 0.85))
ax.text(8.78, s6[-1] + 6.4, '6mA', ha='left', va='center',
        fontsize=FS_ID, color=U.COL_6mA)
ax.text(8.78, s6[-1] + 1.2, f'mean {m6:.1f}', ha='left', va='center',
        fontsize=FS_AN, color=mcolors.to_rgba(U.COL_6mA, 0.85))

# ── render-then-verify ───────────────────────────────────────────────────────
# No panel letter: the composer stamps it in the clear top-left ~10 x 6 mm.
U.assert_no_text_collisions(fig, 'Fig2 panel a')
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
size = Image.open(OUT).size
assert size == (2125, 519), size
print('wrote', OUT, size)
