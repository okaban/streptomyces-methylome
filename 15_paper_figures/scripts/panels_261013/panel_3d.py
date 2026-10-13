"""Figure 3, panel d — funnel of REBASE dual 4mC/6mA systems.

Message: dual 4mC/6mA systems are not rare, but none places both marks on the
same strand.  Three-step funnel of counts (51 -> 50 -> 0); the terminal zero is
drawn as a visible baseline stub plus an explicit "0" label so the slot is never
empty.

Tile geometry is fixed at 677 x 614 px (2.257 x 2.047 in @ 300 dpi) for the
multi-panel composer, which stamps the panel letter itself — the top-left
~10 x 6 mm is left clear here and no letter is drawn.

Run this as a `python` tool cell (environment: epi-rna) so the {{artifact:...}}
marker on the CSV path resolves to a local file path.
"""

import sys, importlib, csv

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # the style helper sets 'tight',
                                             # which would resize the canvas
import matplotlib.pyplot as plt

# ── data ──────────────────────────────────────────────────────────────────────
CSV = "{{artifact:9ce87d3e-5538-4019-95bf-7501ed3fecf4}}"
rows = list(csv.DictReader(open(CSV)))
stages = [r['stage'] for r in rows]
counts = [int(r['n']) for r in rows]
assert counts == [51, 50, 0], counts

LAB = ['Dual 4mC/6mA systems in REBASE',
       'Methylated positions resolvable',
       'Both marks on the same strand']

# ── figure ────────────────────────────────────────────────────────────────────
F = U.FONT
fig = plt.figure(figsize=(2.257, 2.047), dpi=300)
ax = fig.add_axes([0.045, 0.145, 0.945, 0.735])
U.style_axes(ax)

XMAX = 58.0
STUB = 0.95          # data-unit stub so the terminal zero is a visible mark
BH = 0.40            # bar height
ys = [2, 1, 0]       # top-to-bottom funnel rows

for y, n, lab in zip(ys, counts, LAB):
    col = U.COL_GRAY if n == 0 else U.COL_DUAL_SET
    w = STUB if n == 0 else n
    ax.barh(y, w, height=BH, color=col, edgecolor='none', zorder=3)
    # category identity label, sitting above its own bar (replaces y ticks)
    ax.text(0, y + 0.255, lab, ha='left', va='bottom',
            fontsize=F['tick'], color=U.COL_DARK, zorder=4)
    # count at the open end of the bar
    ax.text(w + 1.1, y, f'{n:d}', ha='left', va='center_baseline',
            fontsize=F['tick'], color=U.COL_DARK, zorder=4)

ax.set_xlim(0, XMAX)
ax.set_ylim(-0.42, 2.78)
ax.set_yticks([])
ax.spines['left'].set_visible(False)
ax.set_xticks([0, 25, 50])
ax.set_xlabel('Number of systems', fontsize=F['axis_label'])
ax.tick_params(axis='x', length=2.5, pad=1.5)
ax.xaxis.labelpad = 1.5

# takeaway title, kept clear of the composer's top-left letter zone
fig.text(0.205, 0.985,
         'No reported dual-mark system\nplaces both marks on one strand',
         ha='left', va='top', fontsize=F['panel_title'],
         color=U.COL_DARK, linespacing=1.22)

# ── save (collision assert first, then an exact-size check) ───────────────────
U.assert_no_text_collisions(fig, 'Fig3d funnel')
fig.savefig('panel_d.png', dpi=300, transparent=True,
            bbox_inches=None, pad_inches=0)

from PIL import Image
sz = Image.open('panel_d.png').size
assert sz == (677, 614), sz
print('size', sz)
