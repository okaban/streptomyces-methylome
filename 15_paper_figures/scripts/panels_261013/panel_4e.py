"""Figure 4 panel e — candidate 6mA MTase transcript vs HAAGCCCG occupancy.

Standalone generator for the Fig. 4e tile (2125 x 519 px, transparent PNG).
Run in the `epi-rna` environment; the composer stamps the panel letter, so no
U.add_panel_label() call here.
"""
import sys, importlib
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

U = importlib.import_module('00_shared_utils')
try:
    apply_figure_style()          # figure-style skill helper, when a kernel provides it
except NameError:
    pass
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # the style helper sets 'tight', which resizes the canvas

# ── data ────────────────────────────────────────────────────────────────────
CSV = '{{artifact:11136318-70d0-473e-a86f-0d3fae450b6c}}'   # F4e_enzyme.csv
if CSV.startswith('{{'):   # marker left unresolved (plain `python panel_4e.py`)
    CSV = host.artifact_path('11136318-70d0-473e-a86f-0d3fae450b6c')
df = pd.read_csv(CSV)
assert list(df.timepoint) == ['T1', 'T2', 'T3'], df
tpm = df.SC_RS17645_TPM.to_numpy()
occ = df.occ_6mA_HAAGCCCG_pct.to_numpy()
x = [0, 1, 2]
fold = tpm[0] / tpm[1]

FS_T, FS_L, FS_TK, FS_A = (U.FONT['panel_title'], U.FONT['axis_label'],
                           U.FONT['tick'], U.FONT['annot'])
C_EXPR = U.COL_DARK      # transcript: neutral near-black (not a methyl mark)
C_OCC = U.COL_6mA        # occupancy: the 6mA mark colour

fig = plt.figure(figsize=(2125 / 300, 519 / 300), dpi=300)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
fig.subplots_adjust(left=0.059, right=0.952, bottom=0.115, top=0.845)

# ── series ──────────────────────────────────────────────────────────────────
ax.plot(x, tpm, '-o', color=C_EXPR, lw=1.6, ms=4.5, mfc=C_EXPR,
        mec='white', mew=0.6, zorder=4, clip_on=False)
ax2.plot(x, occ, '-s', color=C_OCC, lw=1.6, ms=4.5, mfc=C_OCC,
         mec='white', mew=0.6, zorder=3, clip_on=False)

ax.set_xlim(-0.20, 2.52)
ax.set_ylim(0, 11.5)
ax.set_yticks([0, 2, 4, 6, 8, 10])
# 35-55 % window: the 1.95 pp spread must read as flat (that is the claim);
# a window tight to the data would magnify it into an apparent anti-correlation.
ax2.set_ylim(35, 55)
ax2.set_yticks([35, 40, 45, 50, 55])
ax.set_xticks(x)
ax.set_xticklabels(U.TP_LABELS, fontsize=FS_TK)

ax.set_ylabel('Transcript level (TPM)', fontsize=FS_L, color=C_EXPR)
ax2.set_ylabel('6mA occupancy (%)', fontsize=FS_L, color=C_OCC)
ax2.yaxis.set_label_coords(1.038, 0.5)

U.style_axes(ax)
for side in ('top', 'left', 'bottom'):
    ax2.spines[side].set_visible(False)
ax2.spines['right'].set_visible(True)
ax2.spines['right'].set_edgecolor(C_OCC)
ax2.spines['right'].set_linewidth(0.8)
ax2.tick_params(axis='y', colors=C_OCC, labelcolor=C_OCC, labelsize=FS_TK)
ax.tick_params(axis='y', colors=C_EXPR, labelcolor=C_EXPR, labelsize=FS_TK)

# ── identity labels (direct text at the right end of each line) ─────────────
ax.text(2.16, tpm[2] - 1.30, 'SC_RS17645', fontsize=FS_A, color=C_EXPR,
        style='italic', ha='left', va='center')
ax2.text(2.16, occ[2] + 2.80, 'HAAGCCCG', fontsize=FS_A, color=C_OCC,
         ha='left', va='center')

# ── the one narrative annotation: the headline number ───────────────────────
ax.text(1.0, 0.62, f'{tpm[1]:.2f} TPM ({fold:.1f}\u00d7 below T1)', fontsize=FS_A,
        color=C_EXPR, ha='center', va='center')

# ── title (italic locus tag + regular clause), left-aligned in the top band ─
r = fig.canvas.get_renderer()
gene = fig.text(0.075, 0.945, 'SC_RS17645', fontsize=FS_T, style='italic',
                color=U.COL_DARK, ha='left', va='center')
gw = gene.get_window_extent(r).width
cands = [
    ' is expressed at every timepoint; transcript falls '
    f'{fold:.1f}\u00d7 at T2 while occupancy holds at '
    f'{occ.min():.1f}\u2013{occ.max():.1f}% (n = 3, no correlation claimed)',
    ' is expressed at every timepoint; transcript falls '
    f'{fold:.1f}\u00d7 at T2, occupancy holds at '
    f'{occ.min():.1f}\u2013{occ.max():.1f}% (n = 3, no correlation)',
    f' expressed at every timepoint; transcript falls {fold:.1f}\u00d7 at T2, '
    f'occupancy {occ.min():.1f}\u2013{occ.max():.1f}% (n = 3, no correlation)',
]
x0 = 0.075 + (gw + 2) / fig.bbox.width
avail = (0.985 - x0) * fig.bbox.width
rest = None
for c in cands:
    probe = fig.text(0, -1, c, fontsize=FS_T)
    w = probe.get_window_extent(r).width
    probe.remove()
    if w <= avail:
        rest = c
        break
assert rest is not None, 'no title candidate fits'
fig.text(x0, 0.945, rest, fontsize=FS_T, color=U.COL_DARK, ha='left',
         va='center')

# ── checks, then save at exactly 2125 x 519 px ──────────────────────────────
U.assert_no_text_collisions(fig, 'Fig4e')
OUT = 'panel_e.png'
fig.savefig(OUT, dpi=300, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
size = Image.open(OUT).size
assert size == (2125, 519), size

# §9.1 geometric re-check: no text leaves the canvas
texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
         if t.get_text().strip() and t.get_visible()]
bad = [t.get_text()[:24] for t, b in texts
       if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.x1 + 1 or b.y1 > fig.bbox.y1 + 1]
assert not bad, bad

# top-left 10 x 6 mm must stay clear for the composer's panel letter
zone = [t.get_text()[:24] for t, b in texts if b.x0 < 118 and b.y1 > 519 - 71]
assert not zone, zone

print(f'{OUT} {size[0]}x{size[1]} fold={fold:.2f} occ_span={occ.max() - occ.min():.2f}pp')
