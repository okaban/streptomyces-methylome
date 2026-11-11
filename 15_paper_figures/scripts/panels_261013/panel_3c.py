"""Figure 3, panel b — single-molecule dual-mark co-occurrence at HAAGCCCG instances.

Two plots on one tile:
  left   stacked composition of read x instance observations per timepoint (reads
         spanning position 1 and position 5 of one HAAGCCCG instance) — WHAT is
         carried on a molecule;
  right  the instance-stratified Mantel-Haenszel odds ratio with its 95% CI against
         a reference line at 1 — whether that is MORE THAN CHANCE.

2026-11-02 (co-author comment, Yuzawa): the single-plot version printed "OR 1.37"
above a bar whose visible quantity was a 29.7% share, and a reader could not connect
the two. Independence predicts 28.03/27.56/27.23% for the both-marks cell, so the
visible excess is only +1.67/+1.98/+1.87 percentage points (observed/expected
1.059/1.072/1.069) while the odds ratio is 1.37. Separating composition from
association puts each number on the chart that can express it.

Environment: epi-rna.  Renders exactly 1354 x 685 px for the Figure 3 composer grid
(row 0, columns 4-11); no panel letter is drawn (the composer stamps it) and
bbox_inches must stay None.
"""
import csv
import sys, importlib
from pathlib import Path

# [FIG18] タイルを再生成できるよう、artifact マーカーをローカルの入力表に解決しておく。
_DATA_DIR = Path(__file__).resolve().parents[2] / 'data_261013'

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # the style helper sets 'tight', which resizes the canvas
import matplotlib.pyplot as plt
import os as _os
_FD = int(_os.environ.get('FIG_DPI', 300))  # 2026-11-09 解像度上書き

CSV_PATH = str(_DATA_DIR / "v458a1f46_F3c_singlemolecule.csv")

rows = list(csv.DictReader(open(CSV_PATH)))
tp = [r['timepoint'] for r in rows]
N = [int(r['observations']) for r in rows]
cats = ['both', 'only_4mC', 'only_6mA', 'neither']
pct = {c: [100.0 * int(r[c]) / int(r['observations']) for r in rows] for c in cats}
# D104 / N26 (2026-10-27): the bar-top label carries the PRIMARY estimator, which
# Methods declares to be the instance-stratified Mantel-Haenszel OR, not the
# unstratified one. Column odds_ratio is kept in the input table for the Methods
# comparison but is no longer what the figure prints.
orr = [float(r['odds_ratio_mh']) for r in rows]
lo = [float(r['mh_ci_low']) for r in rows]
hi = [float(r['mh_ci_high']) for r in rows]
# independence expectation for the both-marks cell, from this bar's own margins
exp_both = [100.0 * ((int(r['both']) + int(r['only_4mC'])) / int(r['observations']))
           * ((int(r['both']) + int(r['only_6mA'])) / int(r['observations']))
           for r in rows]

# Chemistry-coded mapping (author decision F2): colour means the chemistry, never the motif.
COLS = {'both': U.COL_DUAL_SET, 'only_4mC': U.COL_4mC,
        'only_6mA': U.COL_6mA, 'neither': U.COL_GRAY}
NAMES = {'both': '4mC + 6mA', 'only_4mC': '4mC only',
         'only_6mA': '6mA only', 'neither': 'neither'}

# 2026-11-07 共著者コメント C13「この図いらないかな。そんなに興味ない」（右の
# Mantel-Haenszel オッズ比フォレスト）。右プロットを外し、本タイルは組成だけを示す。
# オッズ比 1.367 (95% CI 1.326-1.409) とその層別の説明は本文 R2 と
# Supplementary Table S3b に残るので、数値は失われない。
W_PX, H_PX = 780, 685
fig = plt.figure(figsize=((W_PX + 0.3) / 300, (H_PX + 0.3) / 300), dpi=_FD)
ax = fig.add_axes([0.160, 0.200, 0.560, 0.700])    # composition

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
    ax.plot([i - BW / 2, i + BW / 2], [exp_both[i]] * 2, lw=0.9, ls=(0, (2.2, 1.4)),
            color='white', zorder=6, solid_capstyle='butt')

# identity: direct colour-matched text in the right gutter, not a boxed key
XR = 3.05
for c in cats:
    ax.text(2.40, centres[c][2], NAMES[c], ha='left', va='center',
            fontsize=6.5, color=COLS[c], zorder=5)

U.style_axes(ax, frame=True)
ax.set_xlim(-0.5, XR)
ax.set_ylim(0, 110)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xticks(x)
ax.set_xticklabels(U.TP_LABELS)
ax.spines['bottom'].set_bounds(-0.42, 2.37)
ax.spines['left'].set_bounds(0, 100)
ax.set_ylabel('Observations (%)')
ax.set_xlabel('n = ' + ' / '.join(f'{v:,}' for v in N)
              + '\nread \u00d7 instance pairs')

# 独立のときに期待される「両方」の割合を破線で示す。1 文に短縮（2026-11-07, C13）。
# 2026-11-09 著者コメント [10]: 破線の意味は凡例文に移した（図中に散文を置かない）。

U.assert_no_text_collisions(fig, 'Fig3 panel b', allow_over_data=tuple(value_labels))
fig.savefig(str(Path(__file__).with_suffix('.png')), dpi=_FD, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
assert all(abs(a-round(b*_FD/300))<=2 for a,b in zip(Image.open(str(Path(__file__).with_suffix('.png'))).size, (W_PX, H_PX)))
