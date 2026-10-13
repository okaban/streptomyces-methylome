"""Figure 4 panel d — funnel of counts: 129 SAM-dependent MTases -> 3 candidates.

Standalone generator for the panel-d tile of Figure 4.
Environment: epi-rna.  Output: panel_d.png, exactly 1039 x 590 px, transparent.
"""
import sys, importlib

# project style module (name starts with a digit -> import by string)
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # apply_style() would set 'tight'
import matplotlib.pyplot as plt
import pandas as pd

# stage, n  (129, 128, 21, 12, 5, 3) — resolved from the artifact marker
DATA = "{{artifact:7295f465-86eb-4c1c-9e6d-2874fb07097d}}"
OUT = 'panel_d.png'

df = pd.read_csv(DATA)
stages = df['stage'].tolist()
counts = df['n'].tolist()
final = len(stages) - 1

fig = plt.figure(figsize=(1039 / 300, 590 / 300), dpi=300)
ax = fig.add_axes([0.035, 0.025, 0.955, 0.800])
U.style_axes(ax)

XMAX = 180.0
BARH = 0.34
ax.set_xlim(0, XMAX)
ax.set_ylim(5.78, -0.88)          # inverted: stage 0 on top

for i, (name, n) in enumerate(zip(stages, counts)):
    is_final = (i == final)
    col = U.COL_6mA if is_final else U.COL_GRAY
    ax.barh(i, n, height=BARH, color=col, edgecolor='none', zorder=2)
    # stage identity label: above its own bar, left-aligned at the funnel mouth
    ax.text(0, i - 0.46, name, ha='left', va='center',
            fontsize=U.FONT['tick'],
            color=U.COL_6mA if is_final else U.COL_DARK,
            fontweight='bold' if is_final else 'normal', zorder=3)
    # count label: just past the bar end; the top one carries the unit
    lab = f'{n} proteins' if i == 0 else f'{n}'
    ax.text(n + 2.8, i, lab, ha='left', va='center',
            fontsize=9 if is_final else U.FONT['tick'],
            color=U.COL_6mA if is_final else U.COL_DARK,
            fontweight='bold' if is_final else 'normal', zorder=3)

# every stage is value-labelled, so a count axis would duplicate the encoding
for side in ('left', 'bottom'):
    ax.spines[side].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# narrative block — who survived, and what did not (empty lower-right quadrant)
AX = 82.0
ax.text(AX, 3.95, '6mA candidates: SC_RS17645 (679 aa)', ha='left', va='center',
        fontsize=U.FONT['annot'], color=U.COL_6mA)
ax.text(AX, 4.44, 'plus two BREX-2 PglX proteins', ha='left', va='center',
        fontsize=U.FONT['annot'], color=U.COL_DARK)
ax.text(AX, 4.93, 'No 4mC candidate remains', ha='left', va='center',
        fontsize=U.FONT['annot'], color=U.COL_GRAY)

# title starts at x=0.122 so the composer's panel letter keeps a clear 10x6 mm
fig.text(0.122, 0.955, 'Five filters reduce 129 MTases to three 6mA candidates',
         ha='left', va='top', fontsize=U.FONT['panel_title'], color=U.COL_DARK)

U.assert_no_text_collisions(fig, 'Fig4 panel d')
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
size = Image.open(OUT).size
print('size', size)
assert size == (1039, 590), size
