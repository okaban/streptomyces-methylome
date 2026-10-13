"""Figure 3, panel a — schematic: two recognition sequences, the marks they carry.

Tile for the composer's 12-column grid: EXACT 1039 x 685 px, transparent,
no bbox_inches='tight'.  Colour comes only from 00_shared_utils constants
(chemistry-coded: 4mC -> U.COL_4mC, 6mA -> U.COL_6mA).
"""
import sys, importlib
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None

W, H = 1039, 685
FS_BASE = U.FONT['panel_title']   # 8    - title
FS_SEQ  = 12                      # 12   - sequence letters + mark identity (hero glyphs)
FS_ANN  = U.FONT['annot']         # 6.5  - stats annotations, polarity, bracket

# +0.3 px of slack: matplotlib floors the canvas, and the tile must be EXACT
fig = plt.figure(figsize=((W + 0.3) / 300, (H + 0.3) / 300), dpi=300)
ax = fig.add_axes([0, 0, 1, 1])
U.style_axes(ax)
ax.set_xlim(0, W); ax.set_ylim(0, H)
ax.set_xticks([]); ax.set_yticks([])          # schematic: no scale, no tick Texts
for sp in ax.spines.values():
    sp.set_visible(False)

# ---------------------------------------------------------------- layout grid
X0, DX = 80, 60
PRIME_PAD = 42
ANN_X = 610

DUPLEX = [
    dict(top='TGGCCGGC', bot='ACCGGCCG', y_top=518, y_bot=436,
         marks={4: ('4mC', U.COL_4mC)},
         stats='2,093 positions\n90.1% occupied', bracket=None),
    dict(top='AAGCCCG', bot='TTCGGGC', y_top=252, y_bot=170,
         marks={0: ('6mA', U.COL_6mA), 4: ('4mC', U.COL_4mC)},
         stats='1,334 instances\n4mC 76.7%, 6mA 45.8%\n(HAAGCCCG, n = 565)',
         bracket=(0, 4, '4 bp')),
]

# ------------------------------------------------------------------ the title
fig.text(126 / W, 676 / H,
         'Only AAGCCCG carries 4mC and 6mA on the same strand',
         ha='left', va='top', fontsize=FS_BASE, color=U.COL_DARK,
         linespacing=1.4)

# -------------------------------------------------- strands, marks, stat text
letter_art = {}          # (duplex index, 'top'/'bot', column) -> Text
for d, spec in enumerate(DUPLEX):
    for strand, seq, y, col in (('top', spec['top'], spec['y_top'], U.COL_DARK),
                                ('bot', spec['bot'], spec['y_bot'], U.COL_GRAY)):
        for i, ch in enumerate(seq):
            marked = strand == 'top' and i in spec['marks']
            letter_art[(d, strand, i)] = ax.text(
                X0 + i * DX, y, ch, ha='center', va='center',
                fontsize=FS_SEQ,
                color=spec['marks'][i][1] if marked else col,
                fontweight='bold' if marked else 'normal', zorder=4)
    xr = X0 + (len(spec['top']) - 1) * DX
    for lab, x, y in (("5′", X0 - PRIME_PAD, spec['y_top']),
                      ("3′", xr + PRIME_PAD, spec['y_top']),
                      ("3′", X0 - PRIME_PAD, spec['y_bot']),
                      ("5′", xr + PRIME_PAD, spec['y_bot'])):
        ax.text(x, y, lab, ha='center', va='center', fontsize=FS_ANN,
                color=U.COL_GRAY, zorder=4)
    ax.text(ANN_X, (spec['y_top'] + spec['y_bot']) / 2, spec['stats'],
            ha='left', va='center', fontsize=FS_ANN, color=U.COL_DARK,
            linespacing=1.6, zorder=4)

fig.canvas.draw()
r = fig.canvas.get_renderer()
bb = lambda key: letter_art[key].get_window_extent(r)

# base-pair rungs: measured to start below the top letter and stop above the bottom one
for d, spec in enumerate(DUPLEX):
    for i in range(len(spec['top'])):
        y1, y0 = bb((d, 'top', i)).y0 - 4, bb((d, 'bot', i)).y1 + 4
        ax.add_patch(Rectangle((X0 + i * DX - 0.9, y0), 1.8, y1 - y0,
                               facecolor=U.COL_GRID if False else U.COL_GRAY,
                               edgecolor='none', zorder=1))

# methyl marks: filled circle on the modified base, stem to the base, direct label
R = 11
for d, spec in enumerate(DUPLEX):
    for i, (name, colour) in spec['marks'].items():
        x = X0 + i * DX
        ystem = bb((d, 'top', i)).y1 + 3
        ycirc = ystem + 10 + R
        ax.add_patch(Rectangle((x - 0.9, ystem), 1.8, ycirc - R - ystem,
                               facecolor=colour, edgecolor='none', zorder=2))
        ax.add_patch(Circle((x, ycirc), R, facecolor=colour, edgecolor='none',
                            zorder=3))
        ax.text(x, ycirc + R + 4, name, ha='center', va='bottom',
                fontsize=FS_SEQ, color=colour, fontweight='bold', zorder=4)

# measured bracket for the 4 bp spacing between the two marks on AAGCCCG
for d, spec in enumerate(DUPLEX):
    if spec['bracket'] is None:
        continue
    i, j, lab = spec['bracket']
    xa, xb = X0 + i * DX, X0 + j * DX
    ytop = min(bb((d, 'bot', k)).y0 for k in range(len(spec['bot']))) - 14
    ybar = ytop - 11
    ax.add_patch(Rectangle((xa, ybar), xb - xa, 1.6, facecolor=U.COL_GRAY,
                           edgecolor='none', zorder=2))
    for xe in (xa, xb):
        ax.add_patch(Rectangle((xe - 0.8, ybar), 1.6, ytop - ybar,
                               facecolor=U.COL_GRAY, edgecolor='none', zorder=2))
    ax.text((xa + xb) / 2, ybar - 6, lab, ha='center', va='top',
            fontsize=FS_ANN, color=U.COL_DARK, zorder=4)

U.assert_no_text_collisions(fig, 'Fig3 panel a')
out = 'panel_a.png'
fig.savefig(out, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
sz = Image.open(out).size
assert sz == (1039, 685), sz
print('saved', out, sz)
