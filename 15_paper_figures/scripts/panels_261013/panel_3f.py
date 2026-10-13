"""Figure 3, panel f — sequence context of the AAGCCCG 6mA mark.

LEFT  : 6mA occupancy by the base immediately 5' of AAGCCCG.
RIGHT : TGGCCGGC 4mC vs HAAGCCCG 6mA occupancy across local GC quintiles
        (T1 solid, T2/T3 pale behind).

Renders a composer tile at exactly 2125 x 566 px, transparent background.
Colour means the CHEMISTRY (project rule F2): 4mC -> U.COL_4mC, 6mA -> U.COL_6mA.
Run in the `epi-rna` environment.
"""
import sys, importlib
import pandas as pd, numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

# ---- shared project style module (name starts with a digit) -----------------
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None      # the style helper sets 'tight', which resizes

MINUS = '\u2212'
PRIME = '\u2032'
ARROW = '\u2192'
TIMES = '\u00d7'

OUT_PNG = 'panel_f.png'

# ---- data -------------------------------------------------------------------
FLANK_CSV = "{{artifact:2d0505cc-85d5-410b-9ca6-fc214da8966c}}"   # flank_base, occ_6mA_pct
GC_CSV = "{{artifact:10a64b2d-1d05-439d-bb2c-914ac27421c1}}"      # timepoint, mark, gc_quintile, n, gc_median, meth_pct

flank = pd.read_csv(FLANK_CSV)
gc = pd.read_csv(GC_CSV)

fig = plt.figure(figsize=(2125 / 300, 566 / 300), dpi=300)

# ---- geometry (figure fractions); top-left ~10x6 mm left clear for the
#      composer's panel letter --------------------------------------------
BOT, TOP = 0.247, 0.835
H = TOP - BOT
axL = fig.add_axes([0.065, BOT, 0.186, H])
axR = fig.add_axes([0.468, BOT, 0.373, H])

# ============================ LEFT: 5' flank base ===========================
order = ['C', 'T', 'A', 'G']
vals = [float(flank.loc[flank.flank_base == b, 'occ_6mA_pct'].iloc[0]) for b in order]
xs = np.arange(len(order))
axL.bar(xs, vals, width=0.68, color=U.COL_6mA, edgecolor='none', zorder=3)
axL.set_xticks(xs)
axL.set_xticklabels(order)
axL.set_xlim(-0.65, len(order) - 0.35)
axL.set_ylim(0, 58)
axL.set_yticks([0, 20, 40])
axL.set_ylabel('6mA occupancy (%)')
axL.set_xlabel(f'Base 5{PRIME} of AAGCCCG')
axL.set_title(f'A 5{PRIME} G suppresses 6mA', loc='left', pad=4)
U.style_axes(axL)

# headline value on the mark that carries the claim
axL.text(xs[-1], vals[-1] + 2.4, '8.3', ha='center', va='bottom',
         fontsize=6.5, color=U.COL_6mA, zorder=4)

# narrative annotation -> whitespace between the two sub-plots
fig.text(0.262, TOP - 0.035,
         f'5{PRIME}C : 5{PRIME}G ratio\n'
         f'  in motif       5.55{TIMES}\n'
         f'  non-motif      3.00{TIMES}\n'
         f'  residual       1.85{TIMES}',
         ha='left', va='top', fontsize=6.5, color=U.COL_DARK, linespacing=1.8)

# ============================ RIGHT: GC quintiles ===========================
series = [('TGGCCGGC 4mC', U.COL_4mC, f'93.4 {ARROW} 84.9%', f'\u03b2 = {MINUS}6.4'),
          ('HAAGCCCG 6mA', U.COL_6mA, f'25.1 {ARROW} 83.3%', '\u03b2 = +13.4')]
q = np.arange(1, 6)
for name, col, rng, beta in series:
    sub = gc[gc['mark'] == name]
    for tp in ('T2', 'T3'):                      # pale replicates behind
        y = sub[sub.timepoint == tp].sort_values('gc_quintile')['meth_pct'].to_numpy()
        axR.plot(q, y, '-', color=to_rgba(col, 0.30), lw=2.4, zorder=2,
                 solid_capstyle='round')
    y1 = sub[sub.timepoint == 'T1'].sort_values('gc_quintile')['meth_pct'].to_numpy()
    axR.plot(q, y1, '-', color=col, lw=1.4, marker='o', ms=3.4,
             mfc=col, mec='white', mew=0.6, zorder=4)

axR.set_xticks(q)
axR.set_xticklabels([f'Q{i}' for i in q])
axR.set_xlim(0.80, 5.20)
axR.set_ylim(0, 104)
axR.set_yticks([0, 25, 50, 75, 100])
axR.set_ylabel('Occupancy (%)')
axR.set_xlabel(f'Local GC quintile, motif excluded  (low {ARROW} high)')
axR.set_title('Rising GC lifts 6mA toward the 4mC level, which falls',
              loc='left', pad=4)
U.style_axes(axR)

# direct labels in the right margin, one short leader each
LBLX = 5.62
for (name, col, rng, beta), (yend, ylab) in zip(series, [(84.88, 97.0), (83.33, 48.0)]):
    # leader carries no text of its own, so it cannot become a text-collision finding
    axR.annotate('', xy=(5.04, yend), xytext=(LBLX - 0.06, ylab),
                 annotation_clip=False,
                 arrowprops=dict(arrowstyle='-', color=to_rgba(col, 0.55),
                                 lw=0.7, shrinkA=2, shrinkB=2,
                                 connectionstyle='arc3,rad=0'))
    axR.text(LBLX, ylab, f'{name}\n{rng}\n{beta}',
             ha='left', va='center', fontsize=6.5, color=col,
             linespacing=1.6, clip_on=False)

# ---- checks + save ----------------------------------------------------------
U.assert_no_text_collisions(fig, 'panel_3f')
fig.savefig(OUT_PNG, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
sz = Image.open(OUT_PNG).size
assert sz == (2125, 566), sz
print('size', sz)
