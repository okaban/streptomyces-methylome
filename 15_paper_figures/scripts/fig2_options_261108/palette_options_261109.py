"""Figure 2c logo palette options — comparison sheet only (2026-11-09).

Author feedback: the current logomaker 'colorblind_safe' palette reads "pop",
not rich or refined.  This script draws the two logos of Figure 2c six times,
once per candidate palette, so the palettes can be compared on the real letter
heights.  It writes NOTHING that any other script or the manuscript consumes:
  out  Writing/logo_palette_options_261109.png

Letter heights are read from revA_logo_heights_261108.tsv, the provenance table
the canonical generator (t_revA_logo.py) already wrote, so the comparison uses
exactly the heights of the rendered figure and the canonical generator is not
re-run and not modified.
"""
import csv, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import logomaker

sys.path.insert(0, str(Path(__file__).resolve().parent))
import colour_math as cm   # local copy, see docstring of that file

HERE = Path(__file__).resolve().parent
HEIGHTS = HERE / 'revA_logo_heights_261108.tsv'
OUT = Path('/Users/okaban/obsidian/Research/rna-seq/Writing/logo_palette_options_261109.png')
BASES = ('A', 'C', 'G', 'T')
LO, HI = -6, 6

# --------------------------------------------------------------- candidates
PALETTES = [
 ('Current — logomaker colorblind_safe',
  dict(A='#009980', C='#59B3E6', G='#E69C04', T='#1A1A1A')),
 ('Option 1 — Manuscript harmonic',
  dict(A='#2F5C45', C='#8CA8C0', G='#A5743F', T='#22282E')),
 ('Option 2 — Muted earth',
  dict(A='#4F6B3C', C='#7FA3A0', G='#C08B52', T='#2E2A26')),
 ('Option 3 — Two-hue pairs (purine / pyrimidine)',
  dict(A='#7D4F2A', C='#2B4A60', G='#D7AF86', T='#8FACC0')),
 ('Option 4 — Greyscale with one accent (C)',
  dict(A='#A9B2B8', C='#6B4760', G='#6E7982', T='#1F2225')),
 ('Option 5 — Muted staircase',
  dict(A='#33654D', C='#6C8DAD', G='#D0A97B', T='#23272B')),
]

# ------------------------------------------------------------ letter heights
def matrices():
    rows = {'4mC': {}, '6mA': {}}
    n = {}
    with open(HEIGHTS) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['block'] == 'logo_heights':
                off = int(r['key'].replace('offset_', ''))
                rows[r['mod']][off] = {b: float(r[f'height_{b}']) for b in BASES}
                n[r['mod']] = int(r['n_drawn'])
    out = {}
    for mod, d in rows.items():
        assert sorted(d) == list(range(LO, HI + 1)), sorted(d)
        out[mod] = pd.DataFrame.from_dict(d, orient='index')[list(BASES)]
    return out, n

MAT, NDRAWN = matrices()

# ------------------------------------------------------------------- canvas
mpl.rcParams.update({'font.family': 'Arial', 'savefig.bbox': None,
                     'text.color': '#22282E', 'axes.edgecolor': '#22282E'})
NROW = len(PALETTES)
ROW_H, TOP, BOT = 1.12, 0.46, 0.22          # inches
FIG_W = 10.0
FIG_H = TOP + BOT + NROW * ROW_H
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=200)

LAB_X, LAB_W = 0.012, 0.185
L_X, L_W = 0.215, 0.372
R_X, R_W = 0.608, 0.372
LOGO_H = 0.62                                # fraction of a row

fig.text(L_X, 1 - 0.18 / FIG_H, '4mC — all high-confidence positions '
         f'(n = {NDRAWN["4mC"]:,})', fontsize=8.5, ha='left', va='center')
fig.text(R_X, 1 - 0.18 / FIG_H, '6mA — AAGCCCG context '
         f'(n = {NDRAWN["6mA"]:,})', fontsize=8.5, ha='left', va='center')

for i, (name, pal) in enumerate(PALETTES):
    y0 = (BOT + (NROW - 1 - i) * ROW_H) / FIG_H
    h = ROW_H / FIG_H
    ax_y = y0 + 0.14 * h
    ax_h = LOGO_H * h

    # label + swatches
    fig.text(LAB_X, y0 + 0.86 * h, name, fontsize=8.0, ha='left', va='top',
             weight='bold' if i else 'normal')
    sw = fig.add_axes([LAB_X, ax_y + 0.10 * ax_h, LAB_W * 0.86, ax_h * 0.42])
    sw.set_xlim(0, 4); sw.set_ylim(0, 1); sw.axis('off')
    for k, b in enumerate(BASES):
        sw.add_patch(Rectangle((k + 0.06, 0), 0.88, 1, facecolor=pal[b], lw=0))
        lum = cm.luminance(cm.hex2rgb(pal[b]))
        sw.text(k + 0.5, 0.5, b, ha='center', va='center', fontsize=8.5,
                color='white' if lum < 0.3 else '#15181B')
    for k, b in enumerate(BASES):
        sw.text(k + 0.5, -0.30, pal[b].lstrip('#'), ha='center', va='top',
                fontsize=4.6, color='#555B60', family='monospace')

    for mod, bx, bw in (('4mC', L_X, L_W), ('6mA', R_X, R_W)):
        ax = fig.add_axes([bx, ax_y, bw, ax_h])
        logomaker.Logo(MAT[mod], ax=ax, color_scheme=pal, font_name='Arial',
                       stack_order='big_on_top', vpad=0.02, show_spines=False)
        ax.set_xlim(LO - 0.6, HI + 0.6); ax.set_ylim(0, 2.05)
        ax.set_yticks([0, 1, 2]); ax.set_yticklabels(['0', '1', '2'], fontsize=5.5)
        ax.set_xticks(range(LO, HI + 1))
        ax.set_xticklabels([f'{o:+d}'.replace('+0', '0') for o in range(LO, HI + 1)],
                           fontsize=5.0)
        ax.tick_params(length=1.6, width=0.5, pad=1.2, colors='#22282E')
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        for s in ('left', 'bottom'):
            ax.spines[s].set_visible(True); ax.spines[s].set_linewidth(0.6)
        ax.plot([0], [2.0], marker='v', ms=3.0, color='#606060', clip_on=False, zorder=5)
        if mod == '6mA':
            ax.plot([1], [2.0], marker='v', ms=2.6, mfc='none', mec='#606060',
                    mew=0.6, clip_on=False, zorder=5)

fig.text(LAB_X, 0.06 / FIG_H, 'Offset from the modified base (bp); letter heights '
         'are those of the rendered Figure 2c (revA_logo_heights_261108.tsv). '
         'Comparison sheet only — not a manuscript figure.',
         fontsize=6.0, ha='left', va='bottom', color='#555B60')

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=200, facecolor='white', bbox_inches=None)
print('saved', OUT)
from PIL import Image
print('size', Image.open(OUT).size)
