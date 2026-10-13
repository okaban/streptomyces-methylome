"""Figure 4 panel c — 6mA deficit at transcription start sites, variant ladder."""
import sys, importlib
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial', 'mathtext.bf': 'Arial:bold'})
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pandas as pd

# data artifact marker (resolved to a local path when this source is run as a cell)
DATA = "{{artifact:844c5b6d-22a4-436c-bd32-e6f07e59f11f}}"
df = pd.read_csv(DATA).set_index('variant')

# progression: both biases present -> filter bias removed -> composition matched -> both removed
ORDER = ['filtered, unmatched', 'all reads, unmatched',
         'filtered, GC-matched', 'all reads, GC-matched']
TICK = {'filtered, unmatched':  'filtered reads\nunmatched GC',
        'all reads, unmatched': 'all reads\nunmatched GC',
        'filtered, GC-matched': 'filtered reads\nGC-matched',
        'all reads, GC-matched': 'all reads\nGC-matched'}
d = df.loc[ORDER]
ys = [3, 2, 1, 0]                     # final (quoted) variant at the bottom
focal = ORDER[-1]

fig = plt.figure(figsize=(1039 / 300, 590 / 300), dpi=300)
ax = fig.add_subplot(111)
U.style_axes(ax)
fig.subplots_adjust(left=0.205, right=0.985, top=0.885, bottom=0.245)

REF = 1.0
pale = to_rgba(U.COL_6mA, 0.45)

# reference line: no deficit
ax.axvline(REF, color=U.COL_DARK, lw=0.7, ls=(0, (3, 2)), zorder=1)

for y, (name, row) in zip(ys, d.iterrows()):
    is_focal = name == focal
    # lollipop stem back to the semantic zero (ratio = 1.0)
    ax.plot([row.ratio, REF], [y, y], color=U.COL_GRAY, lw=0.9, zorder=2,
            solid_capstyle='butt')
    if is_focal and pd.notna(row.ci_lo):
        ax.errorbar(row.ratio, y,
                    xerr=[[row.ratio - row.ci_lo], [row.ci_hi - row.ratio]],
                    fmt='none', ecolor=U.COL_6mA, elinewidth=1.3,
                    capsize=2.4, capthick=1.3, zorder=3)
    ax.scatter([row.ratio], [y], s=58 if is_focal else 34,
               facecolor=U.COL_6mA if is_focal else pale,
               edgecolor=U.COL_6mA, linewidth=0.9, zorder=4)

ax.set_yticks(ys)
ax.set_yticklabels([TICK[n] for n in ORDER], fontsize=U.FONT['tick'])
for lab, name in zip(ax.get_yticklabels(), ORDER):
    lab.set_color(U.COL_DARK if name == focal else U.COL_GRAY)
ax.set_ylim(-0.55, 3.85)
ax.set_xlim(0.783, 1.013)
ax.set_xticks([0.80, 0.85, 0.90, 0.95, 1.00])
ax.set_xticklabels(['0.80', '0.85', '0.90', '0.95', '1.00'],
                   fontsize=U.FONT['tick'])
ax.tick_params(axis='y', length=0)

ax.set_xlabel('6mA level near TSS, relative to genome average',
              fontsize=U.FONT['axis_label'], color=U.COL_DARK)
ax.set_title('A 7.9% 6mA deficit at transcription start sites',
             fontsize=U.FONT['panel_title'], color=U.COL_DARK,
             loc='left', pad=5)

# identity label for the reference line
ax.text(0.997, 3.45, 'no deficit', ha='right', va='center',
        fontsize=U.FONT['annot'], color=U.COL_GRAY)
# the one value the text quotes
ax.text(d.loc[focal, 'ratio'], 0.30, '0.921', ha='center', va='bottom',
        fontsize=U.FONT['annot'], color=U.COL_6mA)
# narrative statistics, outside the data area
stats = fig.text(0.985, 0.006,
                 'n = 9 culture units, p = 1.7 \u00d7 10$^{-6}$; '
                 'per-TSS median ratio 0.844 across 117 TSSs',
                 ha='right', va='bottom', fontsize=U.FONT['annot'],
                 color=U.COL_GRAY)

# park the axis label in the middle of the gap between the x tick labels and the
# statistics line; mathtext box heights differ between matplotlib builds, so the
# position is measured rather than hardcoded
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_ticks_bottom = min(t.get_window_extent(_r).y0 for t in ax.get_xticklabels())
_stats_top = stats.get_window_extent(_r).y1
_xl = ax.xaxis.label
_h = _xl.get_window_extent(_r).height
_ax_bb = ax.get_window_extent(_r)
_top = (_ticks_bottom + _stats_top) / 2 + _h / 2
assert _ticks_bottom - _stats_top > _h, 'no room for the axis label'
ax.xaxis.set_label_coords(0.5, (_top - _ax_bb.y0) / _ax_bb.height)

U.assert_no_text_collisions(fig, 'Fig4 panel c')
fig.savefig('panel_c.png', dpi=300, transparent=True,
            bbox_inches=None, pad_inches=0)

from PIL import Image
sz = Image.open('panel_c.png').size
assert sz == (1039, 590), sz
r = fig.canvas.get_renderer()
ax_bb = ax.get_window_extent(r)
print('size', sz, 'axes frac', round(ax_bb.width * ax_bb.height / (1039 * 590), 3))
