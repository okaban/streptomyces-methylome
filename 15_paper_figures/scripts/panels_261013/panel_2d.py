"""Figure 2 panel d — composition of contexts previously reported as 5mC.

Single 100%-stacked bar over the 572 cytosine contexts that earlier work reported
as 5mC in this organism, split into the sequence space of the two motifs we
identify (GCCGGC/GGCCGG family, GCCCG-containing) and everything else.

Renders a TILE for the Figure 2 composer: exactly 677 x 661 px, transparent
background, no tight/constrained layout, top-left 10 x 6 mm left clear for the
composer's panel letter.

Run with environment `epi-rna`:
    python panel_2d.py
"""
import sys
import importlib

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import pandas as pd
from PIL import Image

# ── style ───────────────────────────────────────────────────────────────────
U.apply_unified_style()
# apply_unified_style() inherits savefig.bbox='tight', which silently re-crops
# the canvas and would break the composer's grid.
mpl.rcParams['savefig.bbox'] = None

# ── data ────────────────────────────────────────────────────────────────────
DATA_CSV = "{{artifact:45b56046-0596-4fcf-904e-ac5c204efe55}}"
OUT_PNG = 'panel_d.png'

df = pd.read_csv(DATA_CSV)
counts = dict(zip(df['class'], df['n']))
TOT = int(df['n'].sum())
F = counts['GCCGGC / GGCCGG family']
G = counts['GCCCG-containing']
O = counts['other']
assert F + G + O == TOT

pf, pg, po = 100 * F / TOT, 100 * G / TOT, 100 * O / TOT
p_inscope = pf + pg

# Pale tint of the 6mA mark colour, DERIVED from the palette constant (blended
# to white rather than set with alpha, so the tile stays correct when saved on a
# transparent background). Never type a new hex — project figure rule 2.
PALE_6mA = tuple(0.45 * c + 0.55 for c in to_rgb(U.COL_6mA))

# ── figure ──────────────────────────────────────────────────────────────────
# 677/300 x 661/300 in, not the rounded 2.257 x 2.203 — the latter truncates to
# 660 px of height and fails the composer's size check.
fig = plt.figure(figsize=(677 / 300, 661 / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.175, right=0.995, bottom=0.085, top=0.845)

XC, BW = 0.0, 0.76
segs = [
    (pf, U.COL_4mC, 'GCCGGC /\nGGCCGG family', 'white'),
    (pg, PALE_6mA, 'GCCCG-containing', U.COL_DARK),
    (po, U.COL_GRAY, 'other', U.COL_DARK),
]
inside = []
bottom = 0.0
for frac, col, lab, tcol in segs:
    ax.bar(XC, frac, width=BW, bottom=bottom, color=col,
           edgecolor='none', zorder=2)
    t = ax.text(XC, bottom + frac / 2, lab, ha='center', va='center',
                color=tcol, fontsize=U.FONT['annot'], linespacing=1.15,
                zorder=4)
    inside.append(t.get_text()[:28])
    bottom += frac

# Bracket over the two in-scope classes — the one value label of the panel.
xb = BW / 2 + 0.05
ax.plot([xb - 0.045, xb, xb, xb - 0.045],
        [0.6, 0.6, p_inscope - 0.6, p_inscope - 0.6],
        color=U.COL_DARK, lw=0.8, solid_capstyle='butt', zorder=3, clip_on=False)
ax.text(xb + 0.055, p_inscope / 2, f'{p_inscope:.1f}%\nin motif\nspace',
        ha='left', va='center', color=U.COL_DARK,
        fontsize=U.FONT['annot'], linespacing=1.2, zorder=4)

ax.set_xlim(-0.5, 0.80)
ax.set_ylim(0, 100)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xticks([XC])
ax.set_xticklabels([f'n = {TOT}'])
ax.tick_params(axis='x', length=0)
ax.set_ylabel('% of reported 5mC sites', labelpad=2)
ax.set_title('85% of previously reported\n5mC sites lie in our motifs',
             loc='left', x=0.075, linespacing=1.2)
U.style_axes(ax)

# ── save + verify ───────────────────────────────────────────────────────────
fig.savefig(OUT_PNG, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

# The three in-segment texts are identity labels deliberately placed on a bar,
# which is the sanctioned use of allow_over_data.
U.assert_no_text_collisions(fig, 'Fig2d', allow_over_data=tuple(inside))

assert Image.open(OUT_PNG).size == (677, 661), Image.open(OUT_PNG).size
print('size ok (677, 661)')
print('pct', round(pf, 2), round(pg, 2), round(po, 2), round(p_inscope, 3))
