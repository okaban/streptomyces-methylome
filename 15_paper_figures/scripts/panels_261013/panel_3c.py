"""Figure 3, panel c — single-molecule dual-mark composition at HAAGCCCG instances.

Stacked composition of read x instance observations per timepoint (reads spanning
position 1 and position 5 of one HAAGCCCG instance), with the both-marks percentage
printed inside the bottom segment and the within-molecule odds ratio above each bar.

Environment: epi-rna.  Renders exactly 1039 x 685 px for the Figure 3 composer grid;
no panel letter is drawn (the composer stamps it) and bbox_inches must stay None.
"""
import sys, importlib, csv

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # the style helper sets 'tight', which resizes the canvas
import matplotlib.pyplot as plt

CSV_PATH = "{{artifact:458a1f46-dc23-4a9b-80f1-b237b8db510f}}"

rows = list(csv.DictReader(open(CSV_PATH)))
tp = [r['timepoint'] for r in rows]
N = [int(r['observations']) for r in rows]
cats = ['both', 'only_4mC', 'only_6mA', 'neither']
pct = {c: [100.0 * int(r[c]) / int(r['observations']) for r in rows] for c in cats}
orr = [float(r['odds_ratio']) for r in rows]

# Chemistry-coded mapping (author decision F2): colour means the chemistry, never the motif.
COLS = {'both': U.COL_DUAL_SET, 'only_4mC': U.COL_4mC,
        'only_6mA': U.COL_6mA, 'neither': U.COL_GRAY}
NAMES = {'both': '4mC + 6mA', 'only_4mC': '4mC only',
         'only_6mA': '6mA only', 'neither': 'neither'}

fig = plt.figure(figsize=(1039 / 300, 685 / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.12, right=0.99, top=0.90, bottom=0.175)

x = [0, 1, 2]
BW = 0.74
bottom = [0.0, 0.0, 0.0]
centres = {}
for c in cats:
    ax.bar(x, pct[c], BW, bottom=bottom, color=COLS[c],
           edgecolor='white', linewidth=0.6, zorder=3)
    centres[c] = [b + h / 2 for b, h in zip(bottom, pct[c])]
    bottom = [b + h for b, h in zip(bottom, pct[c])]

# value label (the number carrying the claim) inside the both-marks segment,
# and the within-molecule odds ratio above each bar
value_labels = []
for i in x:
    t = ax.text(i, centres['both'][i], f"{pct['both'][i]:.1f}", ha='center',
                va='center', color='white', fontsize=7.5, zorder=5)
    value_labels.append(t.get_text())
    ax.text(i, 101.5, f"OR {orr[i]:.2f}", ha='center', va='bottom',
            fontsize=6.5, color=U.COL_DUAL_SET, zorder=5)

# identity: direct colour-matched text in the right gutter, not a boxed key
XR = 3.00
for c in cats:
    ax.text(2.40, centres[c][2], NAMES[c], ha='left', va='center',
            fontsize=6.5, color=COLS[c], zorder=5)

U.style_axes(ax)
ax.set_xlim(-0.5, XR)
ax.set_ylim(0, 110)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xticks(x)
ax.set_xticklabels(U.TP_LABELS)
ax.spines['bottom'].set_bounds(-0.42, 2.37)
ax.spines['left'].set_bounds(0, 100)
ax.set_ylabel('Observations (%)')
ax.set_xlabel('n = 84,047 / 79,313 / 60,651 read \u00d7 instance pairs')

# takeaway title, kept clear of the composer's top-left panel-letter reserve (~10 x 6 mm)
fig.text(0.145, 0.99, 'Both marks co-occur on one molecule more than chance',
         ha='left', va='top', fontsize=8, color=U.COL_DARK)

U.assert_no_text_collisions(fig, 'Fig3 panel c', allow_over_data=tuple(value_labels))
fig.savefig('panel_c.png', dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
assert Image.open('panel_c.png').size == (1039, 685)
