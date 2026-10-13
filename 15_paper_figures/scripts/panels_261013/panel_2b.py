"""Figure 2 panel b — paired distributions of the 4mC and 5mC basecalling
channels at the same 2,091 TGGCCGGC positions.

Renders a 677 x 661 px transparent tile for the Figure 2 composer:
the 4mC channel (U.COL_4mC) is heavily modified (median 82%), while the
5mC channel (U.COL_GRAY) measured on the same reads at the same positions
is empty (median 0%).

Run as a cell in the `epi-rna` environment: the {{artifact:...}} marker below
is rewritten to the CSV's local path before execution. Outside the harness,
pass the CSV path as argv[1].
"""
import sys, importlib
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # the style helper sets 'tight', which resizes the canvas

# site-level channel levels: level_4mC_channel_pct, level_5mC_channel_pct (n = 2,091)
CSV = "{{artifact:1bc294f3-dd77-4616-8326-9d0e54718441}}"
if len(sys.argv) > 1:
    CSV = sys.argv[1]
OUT = 'panel_b.png'

df = pd.read_csv(CSV)
a = df['level_4mC_channel_pct'].to_numpy(float)
b = df['level_5mC_channel_pct'].to_numpy(float)
n = len(df)
med_a, med_b = np.median(a), np.median(b)

# darker grey for the 5mC outline/text: derived from U.COL_GRAY, never a new hex
GREY_DARK = tuple(np.array(to_rgb(U.COL_GRAY)) * 0.72)

fig = plt.figure(figsize=(677 / 300, 661 / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.035, right=0.985, bottom=0.175, top=0.795)

BINS = np.arange(-1.25, 102.5, 2.5)   # 2.5%-wide bins centred on 0, 2.5, ..., 100
CENT = BINS[:-1] + 1.25
H = 0.30                              # ridge height; each ridge scaled to its own peak


def ridge(vals, base, fill, edge, med):
    """One horizontal density ridge on baseline `base`, with a median tick."""
    h, _ = np.histogram(vals, bins=BINS)
    h = h / h.max() * H
    ax.fill_between(CENT, base, base + h, step='mid', facecolor=fill,
                    edgecolor='none', alpha=0.92, zorder=2)
    nz = np.flatnonzero(h > 0)
    i0, i1 = max(nz[0] - 1, 0), min(nz[-1] + 1, len(h) - 1)
    ax.step(CENT[i0:i1 + 1], base + h[i0:i1 + 1], where='mid',
            color=edge, linewidth=0.7, zorder=3)
    ax.plot([med, med], [base, base + H * 0.98], color=U.COL_DARK,
            linewidth=1.1, solid_capstyle='butt', zorder=4)


ridge(b, 0.00, U.COL_GRAY, GREY_DARK, med_b)
ridge(a, 0.88, U.COL_4mC, U.COL_4mC, med_a)

# identity labels adjacent to each ridge (no legend box); medians + n are narrative
ax.text(0, 1.23, '4mC channel', color=U.COL_4mC, fontsize=7.5,
        ha='left', va='bottom')
ax.text(med_a, 1.23, f'median {med_a:.0f}%', color=U.COL_4mC, fontsize=6.5,
        ha='center', va='bottom')
ax.text(0, 0.35, '5mC channel', color=GREY_DARK, fontsize=7.5,
        ha='left', va='bottom')
ax.text(5.0, 0.10, f'median {med_b:.0f}%', color=GREY_DARK, fontsize=6.5,
        ha='left', va='bottom')
ax.text(103, 0.56, f'n = {n:,} sites', color=U.COL_DARK, fontsize=6.5,
        ha='right', va='bottom')

ax.set_xlim(-4, 103)
ax.set_ylim(-0.05, 1.40)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks([])
ax.set_xlabel('Methylated reads per TGGCCGGC site (%)')
ax.set_title('Only the 4mC channel is methylated', loc='left', pad=5)
U.style_axes(ax)
ax.spines['left'].set_visible(False)   # categorical rows: no y axis
ax.tick_params(axis='y', length=0)

# tile, not a standalone figure: exact pixel size, no tight bbox
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)
U.assert_no_text_collisions(fig, 'Fig2b')

from PIL import Image
assert Image.open(OUT).size == (677, 661), Image.open(OUT).size
print('size', Image.open(OUT).size)
print('median 4mC %.4f  median 5mC %.4f  n %d' % (med_a, med_b, n))
