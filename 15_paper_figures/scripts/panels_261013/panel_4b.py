"""Figure 4 panel b — equivalence test on promoter-proximal occupancy vs expression change.

Renders the tile `panel_b.png` at exactly 1039 x 590 px for the Figure 4 composer.
Run in the `epi-rna` environment:  exec(open('panel_4b.py').read())
"""
import sys, importlib

# project style module (name starts with a digit, so import it by string)
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.transforms import blended_transform_factory
import pandas as pd

try:                                    # figure-style floor, if the skill sidecar is loaded
    apply_figure_style()                # noqa: F821
except Exception:
    pass
U.apply_unified_style()                 # project standard wins
mpl.rcParams['savefig.bbox'] = None     # style helpers set 'tight', which resizes the canvas
mpl.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
                     'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold'})

# ── Data ────────────────────────────────────────────────────────────────────
# statistic, value (Cliff's delta), equivalence_margin_low/high, comparison
CSV = '{{artifact:aa6f3d34-d738-44fd-9670-41d5ed57de7e}}'
if CSV.startswith('{{'):                # unresolved marker (file run directly, not as a cell)
    CSV = host.artifact_path('aa6f3d34-d738-44fd-9670-41d5ed57de7e')  # noqa: F821
d = pd.read_csv(CSV).iloc[0]
delta = float(d['value'])
lo, hi = float(d['equivalence_margin_low']), float(d['equivalence_margin_high'])

MINUS = '\u2212'                        # U+2212, never hyphen-minus
def sgn(v, nd=2):
    s = f'{abs(v):.{nd}f}'
    return (MINUS if v < 0 else '') + s

# ── Canvas ──────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(1039 / 300, 590 / 300), dpi=300)
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.045, right=0.995, bottom=0.195, top=0.805)

XL = 0.42
BAND_Y0 = 0.14          # band floor in axes fraction — leaves a clear strip for the value label
DOT_Y = 0.59

# ── Equivalence margin: the zone inside which an effect counts as negligible ─
ax.axvspan(lo, hi, ymin=BAND_Y0, ymax=1.0, facecolor=to_rgba(U.COL_GRAY, 0.20),
           edgecolor='none', zorder=1)
for b in (lo, hi):
    ax.plot([b, b], [BAND_Y0, 1.0], transform=blended_transform_factory(ax.transData, ax.transAxes),
            color=U.COL_GRAY, lw=0.9, ls=(0, (3, 2)), zorder=2, solid_capstyle='butt')
ax.plot([0, 0], [BAND_Y0, 1.0], transform=blended_transform_factory(ax.transData, ax.transAxes),
        color=U.COL_GRAY, lw=0.7, ls=(0, (1, 2)), zorder=2)        # semantic zero

# ── The observed effect ─────────────────────────────────────────────────────
ax.plot([delta], [DOT_Y], transform=blended_transform_factory(ax.transData, ax.transAxes),
        marker='o', ms=7.5, mfc=U.COL_DARK, mec='white', mew=0.9, ls='none', zorder=5,
        clip_on=False)

# identity + headline value for the single mark
ax.text(delta, 0.065, f'observed  \u03b4 = {sgn(delta)}',
        transform=blended_transform_factory(ax.transData, ax.transAxes),
        ha='center', va='center', fontsize=6.5, color=U.COL_DARK, zorder=6)
# identity for the shaded zone, placed outside the data area
ax.text(0.0, 1.035, f'equivalence margin  |\u03b4| \u2264 {hi:.2f}  (negligible)',
        transform=blended_transform_factory(ax.transData, ax.transAxes),
        ha='center', va='bottom', fontsize=6.5, color=U.COL_GRAY, zorder=6)

# ── Axes ────────────────────────────────────────────────────────────────────
ax.set_xlim(-XL, XL)
ax.set_ylim(0, 1)
ax.set_xticks([lo, 0, hi])                       # margin bounds are self-labelled by the ticks
ax.set_xticklabels([sgn(lo), '0', '+' + f'{hi:.2f}'])
ax.set_yticks([])
ax.set_xlabel("Cliff's $\\delta$:  promoter-proximal occupancy vs. log$_2$FC (T2/T1)")
ax.set_title('Promoter occupancy is unrelated to expression change',
             loc='left', x=0.085, pad=11)        # x offset clears the composer's panel letter

U.style_axes(ax)
ax.spines['left'].set_visible(False)             # y carries no meaning here
ax.tick_params(axis='y', left=False)
ax.tick_params(axis='x', length=3, width=0.8)

# ── Save (tile: exact pixel size, no tight bbox) ────────────────────────────
U.assert_no_text_collisions(fig, 'Fig4b equivalence')
fig.savefig('panel_b.png', dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
sz = Image.open('panel_b.png').size
assert sz == (1039, 590), sz

# §9.1 geometric re-check on the saved canvas
r = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
bad = [(a.get_text()[:20], b.get_text()[:20]) for i, (a, ba) in enumerate(texts)
       for b, bb in texts[i + 1:] if ba.overlaps(bb)]
outside = [t.get_text()[:20] for t, b in texts
           if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.x1 + 1 or b.y1 > fig.bbox.y1 + 1]
assert not bad and not outside, (bad, outside)
print('size', sz, '| text-text overlaps', bad, '| outside canvas', outside)
