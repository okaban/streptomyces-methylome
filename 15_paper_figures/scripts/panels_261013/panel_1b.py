"""Figure 1 panel b — growth and antibiotic production time course.

Self-contained producing script for the Figure 1b tile (1039 x 614 px,
transparent PNG). Environment: conda env `epi-rna`.

Source data: long-format CSV, 36 rows = 3 measurements x 3 timepoints x
4 biological replicates (columns: measurement, timepoint, hours,
replicate, value, unit).

The DATA line holds the artifact marker literally. Inside Claude Science it is
resolved to a local path when the source is submitted to a python/bash cell; to
run the file standalone, substitute the marker with the CSV path, e.g.

    DATA = host.artifact_path('054306f1-184b-4cdd-9891-fa1e464dbd6e')
"""
import sys, importlib

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # the style helper sets 'tight', which resizes the canvas
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

DATA = "{{artifact:054306f1-184b-4cdd-9891-fa1e464dbd6e}}"   # act/red/DCW time course, long format
OUT = 'panel_b.png'

df = pd.read_csv(DATA)

# Antibiotics are secondary metabolites, not methylation marks: neutral + green/tan
# hues only, the 4mC/6mA mark palette is deliberately not used in this panel.
SERIES = {
    'dry_cell_weight_g_L':     dict(label='dry cell weight',    color=U.COL_DARK, marker='o', axis='L'),
    'actinorhodin_mg_L':       dict(label='actinorhodin',       color=U.COL_CORE, marker='s', axis='R'),
    'undecylprodigiosin_mg_L': dict(label='undecylprodigiosin', color=U.COL_ARM,  marker='^', axis='R'),
}
HOURS = sorted(df['hours'].unique())
TP = dict(zip(HOURS, ['T1', 'T2', 'T3']))
NREP = df['replicate'].nunique()

fig = plt.figure(figsize=(1039 / 300, 614 / 300), dpi=300)
axL = fig.add_subplot(111)
axR = axL.twinx()
fig.subplots_adjust(left=0.115, right=0.855, bottom=0.155, top=0.885)

XLIM = (9.0, 56.0)
axL.set_xlim(*XLIM)
axL.set_ylim(0, 7.5)
axR.set_yscale('log')
axR.set_ylim(6e-4, 2.5e1)

for key, spec in SERIES.items():
    ax = axL if spec['axis'] == 'L' else axR
    sub = df[df['measurement'] == key]
    # individual biological replicates, pale, behind the summary
    for rep, g in sub.groupby('replicate'):
        g = g.sort_values('hours')
        v = g['value'].to_numpy(dtype=float)
        # a zero cannot be drawn on the log axis; mask it and mark it as n.d.
        v = np.ma.masked_less_equal(v, 0.0) if spec['axis'] == 'R' else v
        ax.plot(g['hours'], v, color=spec['color'], lw=0.5, alpha=0.32,
                marker='o', ms=1.4, mew=0, zorder=2)
    # replicate-level summary: mean +/- s.d. (the zero IS included in the mean)
    m = sub.groupby('hours')['value'].agg(['mean', 'std']).reindex(HOURS)
    ax.errorbar(HOURS, m['mean'], yerr=m['std'], color=spec['color'], lw=1.4,
                marker=spec['marker'], ms=3.4, mew=0.0,
                elinewidth=0.9, capsize=1.8, capthick=0.9, zorder=4)

axL.set_xticks(HOURS)                        # ticks mark the three sampling times
axL.set_xticklabels([f'{TP[h]}\n{h:g} h' for h in HOURS])
axL.set_yticks([0, 2, 4, 6])
axR.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
axR.set_yticklabels(['0.001', '0.01', '0.1', '1', '10'])
axR.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

axL.set_ylabel('Dry cell weight (g/L)')
axR.set_ylabel('Antibiotic titre (mg/L)')

U.style_axes(axL)
U.style_axes(axR)
axR.spines['right'].set_visible(True)        # the twin axis needs its own spine back
axR.spines['right'].set_edgecolor(U.COL_DARK)
axR.spines['right'].set_linewidth(0.8)
axR.spines['left'].set_visible(False)
axR.spines['bottom'].set_visible(False)

axL.set_title('Biomass plateaus by T2; antibiotics only by T3',
              loc='left', x=0.0, y=1.015, fontsize=U.FONT['panel_title'])

# series identity — frameless direct key in the empty lower-right wedge
key_x0, key_x1, key_tx = 0.455, 0.525, 0.545
for i, spec in enumerate(SERIES.values()):
    y = 0.335 - 0.115 * i
    axL.plot([key_x0, key_x1], [y, y], transform=axL.transAxes, color=spec['color'],
             lw=1.4, solid_capstyle='round', clip_on=False, zorder=6)
    axL.plot([(key_x0 + key_x1) / 2], [y], transform=axL.transAxes, color=spec['color'],
             marker=spec['marker'], ms=3.4, mew=0.0, clip_on=False, zorder=6)
    axL.text(key_tx, y, spec['label'], transform=axL.transAxes, ha='left', va='center',
             fontsize=U.FONT['legend'], color=spec['color'])

axL.text(0.995, 0.035, f'n = {NREP}', transform=axL.transAxes, ha='right', va='bottom',
         fontsize=U.FONT['annot'], color=U.COL_DARK)

# one actinorhodin replicate is zero at T1 and cannot sit on a log axis
axR.text(12.9, 1.1e-3, 'n.d.', ha='left', va='center',
         fontsize=U.FONT['annot'], color=SERIES['actinorhodin_mg_L']['color'])

U.assert_no_text_collisions(fig, 'Fig1b growth/production')
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

# --- render-then-verify: exact tile size, no text outside the canvas, no text on a spine
from PIL import Image
sz = Image.open(OUT).size
print('size', sz)
assert sz == (1039, 614), sz

r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
bad = [(t.get_text()[:24]) for t, b in texts
       if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.x1 + 1 or b.y1 > fig.bbox.y1 + 1]
print('outside canvas:', bad)
assert not bad

spines = [(s, s.get_window_extent(r)) for ax in fig.axes for s in ax.spines.values() if s.get_visible()]
tl = {ax: set(ax.get_xticklabels(which='both') + ax.get_yticklabels(which='both')) for ax in fig.axes}
cross = [(t.get_text()[:20], s.axes.get_label() or 'ax') for t, bt in texts for s, bs in spines
         if bt.overlaps(bs) and t not in tl.get(s.axes, set())]
print('text crossing a spine:', cross)
assert not cross

print('means\n', df.groupby(['measurement', 'hours'])['value'].mean().round(4).to_string())
