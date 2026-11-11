# 2026-11-05 D188: 補足図 12→9 の削減。出荷名 SuppFigure11.png -> SuppFigure9.png。
#!/usr/bin/env python3
"""Supplementary Figure 15 — duplex schematic of the two recognition sequences.

Relocated from Figure 3a on 2026-10-22: the author approved replacing the main
figure's schematic with sequence logos built from the canonical census, on the
condition that the schematic itself move to the supplementary material rather
than be discarded.  The drawing is the former `panels_261013/panel_3a.py`
re-rendered at the 174 mm supplementary column width.

What this adds over Supplementary Figure 1, which also schematises both
sequences: Supplementary Figure 1 maps mark positions on a single strand and
pairs them with the measured per-position modified fraction.  This figure shows
the DUPLEX — both strands with their 5'/3' polarity and base pairing — and the
4 bp same-strand spacing between the two marks of AAGCCCG, which is the
configuration the paper reports as unprecedented.

DATA INPUT: none.  This is a schematic.  The four occupancy numbers annotated
beside each duplex are literals carried over from the motif/occupancy analysis
and are registered in Writing/value_registry.csv:
    TGGCCGGC : 2,093 positions, 90.1% occupied
    AAGCCCG  : 1,334 instances; 4mC 76.7%, 6mA 24.2%
               (6mA 45.8% on the 565 whose 5' base is not G)

Colour encodes the CHEMISTRY, not the motif (4mC -> RED, 6mA -> BLUE, both
imported from 00_shared_utils via _suppstyle).  FILL encodes which position of a
mark is being shown: a filled circle is the principal modified position, an open
circle the fallback position used when the principal one is blocked.  Figure 3a
carries the same filled/open contract on its logo markers (filled and open
triangles), so a reader meeting one figure first can read the other.

2026-10-24 (consistency audit CA06): the open circle on the second adenine of
AAGCCCG was added.  Until today this figure marked the adenine mark at position 1
only, while the body text, Figure 3b's legend and Supplementary Figure 2c's
legend all give position 2 as a fallback -- the identical defect the author found
in the Figure 3a logo, which this drawing inherited when it was relocated out of
Figure 3a on 2026-10-22.

Run:  PYTHONPATH=. python suppfigure15_motif_duplex_schematic.py
Out:  $SUPPFIG_OUT/SuppFigure11.png  (default Writing/fig_images/)
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

from _suppstyle import (BLUE, DPI, GREY, NAR_WIDTH_IN, RED, apply_figure_style,
                        report, save)

apply_figure_style(frame="boxed", sizes=(7.5, 6.5, 6.0))

# ---------------------------------------------------------------------- canvas
# 1 data unit == 1 px so the mark circles stay round and the measured-extent
# placement below works in the same units the renderer reports.
# 2026-10-24 (CA06): the explanatory note then at the foot of the canvas gained a
# line, and at the design height of 690 px its top ran into the "4 bp, same
# strand" bracket label (assert_no_text_collisions caught it), so the canvas was
# grown by one 6 pt line (36 px) and every drawn y raised by the same amount.
# 2026-11-08 (author rule: no prose inside a figure): that note moved to the
# legend.  The reserved band is reclaimed by the same one-line mechanism run
# backwards -- every drawn y drops by Y_SHIFT_PX and the canvas loses the same
# amount, so the schematic keeps its internal spacing and the foot margin comes
# back to the ~22 px the rest of the supplementary set uses.
Y_SHIFT_PX = -84
# 2026-11-09: この図は版面を「300 dpi の画素」で直接書き、文字の実測幅（表示画素）を
# データ座標として使い回している。両者が一致するのは図の dpi が 300 のときだけなので、
# 作図は 300 で固定し、解像度は保存時のラスタ化（save_figure_raster の dpi）だけに効かせる。
_GEO = 300
# 2026-11-09: 版面を 174 mm 幅で取ると右に 60 mm の白が残る（図は左 2/3 しか使わない）。
# 投稿規定は余白を詰めることを求めるので、内容の実測幅 114 mm に canvas を合わせる。
W, H = 1400, 690 + Y_SHIFT_PX
FS_SEQ = 13                      # sequence letters + mark identity (hero glyphs)
FS_ANN = 7.0                     # stats annotations, polarity, bracket

fig = plt.figure(figsize=(W / _GEO, H / _GEO), dpi=_GEO)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)

# ----------------------------------------------------------------- layout grid
X0, DX = 150, 72       # x of column 1, column pitch (px)
PRIME_PAD = 50         # gap from the terminal base to its 5'/3' label
ANN_X = 900            # left edge of the stats annotations

# Marks are keyed by 0-based column index; positions in prose are 1-based from
# the 5' end (column index 4 == position 5, the cytosine of both motifs).
DUPLEX = [
    dict(top='TGGCCGGC', bot='ACCGGCCG',
         y_top=530 + Y_SHIFT_PX, y_bot=448 + Y_SHIFT_PX,
         marks={4: ('4mC', RED)}, fallback={},
         stats='2,093 strand-oriented positions\n90.1% occupied',
         bracket=None),
    dict(top='AAGCCCG', bot='TTCGGGC',
         y_top=278 + Y_SHIFT_PX, y_bot=196 + Y_SHIFT_PX,
         marks={0: ('6mA', BLUE), 4: ('4mC', RED)},
         # 2026-10-24 (CA06): the adenine mark also occupies position 2 when the
         # 5' base is G blocks position 1, and this figure marked position 1
         # only — the same defect the author found in the Figure 3a logo, which
         # survived here when the schematic was relocated on 2026-10-22.  The
         # body text (Supplementary Results R3) reports position 2 carrying the
         # mark in 6.97%, 7.50% and 8.85% of instances at T1-T3, Figure 3b's
         # legend reads "A1 by 5' base, A2", and Supplementary Figure 2c's reads
         # "6mA at A1/A2".  An OPEN circle marks it, against the FILLED circle of
         # the principal position -- the same filled/open contract Figure 3a now
         # uses for its filled and open triangles.
         fallback={1: ('fallback', BLUE)},
         stats='1,334 instances\n4mC 76.7%, 6mA 24.2%\n'
               '(6mA 45.8% where 5\u2032 is not G)',
         bracket=(0, 4, '4 bp, same strand')),
]

# --------------------------------------------------- strands, polarity, stats
letter_art = {}          # (duplex index, 'top'/'bot', column) -> Text
for d, spec in enumerate(DUPLEX):
    for strand, seq, y, col in (('top', spec['top'], spec['y_top'], '0.15'),
                                ('bot', spec['bot'], spec['y_bot'], GREY)):
        for i, ch in enumerate(seq):
            marked = strand == 'top' and i in spec['marks']
            letter_art[(d, strand, i)] = ax.text(
                X0 + i * DX, y, ch, ha='center', va='center',
                fontsize=FS_SEQ,
                color=spec['marks'][i][1] if marked else col,
                fontweight='bold' if marked else 'normal', zorder=4)
    xr = X0 + (len(spec['top']) - 1) * DX
    for lab, x, y in (("5\u2032", X0 - PRIME_PAD, spec['y_top']),
                      ("3\u2032", xr + PRIME_PAD, spec['y_top']),
                      ("3\u2032", X0 - PRIME_PAD, spec['y_bot']),
                      ("5\u2032", xr + PRIME_PAD, spec['y_bot'])):
        ax.text(x, y, lab, ha='center', va='center', fontsize=FS_ANN,
                color=GREY, zorder=4)
    ax.text(ANN_X, (spec['y_top'] + spec['y_bot']) / 2, spec['stats'],
            ha='left', va='center', fontsize=FS_ANN, color='0.15',
            linespacing=1.6, zorder=4)

# Every mark and rung below is placed from the MEASURED extent of the letter it
# belongs to, so no glyph can ever collide with a drawn mark.
fig.canvas.draw()
r = fig.canvas.get_renderer()


def bb(key):
    return letter_art[key].get_window_extent(r)


# base-pair rungs: start below the top letter, stop above the bottom one
for d, spec in enumerate(DUPLEX):
    for i in range(len(spec['top'])):
        y1, y0 = bb((d, 'top', i)).y0 - 4, bb((d, 'bot', i)).y1 + 4
        ax.add_patch(Rectangle((X0 + i * DX - 0.9, y0), 1.8, y1 - y0,
                               facecolor=GREY, edgecolor='none', zorder=1))

# methyl marks: filled circle on the modified base, stem to the base, direct label
R = 12
for d, spec in enumerate(DUPLEX):
    for i, (name, colour) in spec['marks'].items():
        x = X0 + i * DX
        ystem = bb((d, 'top', i)).y1 + 3
        ycirc = ystem + 11 + R
        ax.add_patch(Rectangle((x - 0.9, ystem), 1.8, ycirc - R - ystem,
                               facecolor=colour, edgecolor='none', zorder=2))
        ax.add_patch(Circle((x, ycirc), R, facecolor=colour, edgecolor='none',
                            zorder=3))
        ax.text(x, ycirc + R + 5, name, ha='center', va='bottom',
                fontsize=FS_SEQ, color=colour, fontweight='bold', zorder=4)

# fallback positions: OPEN circle, same stem, label set beside the circle rather
# than above it.  Above would put a 7 pt word on the same baseline as the 13 pt
# mark-identity label one column to its left and the two overlap; beside the
# circle the word runs into the gap between position 2 and the 4mC at position 5,
# which is empty at this height.  The identity is not repeated on the open mark:
# it is the same chemistry on the same strand as the filled one, and the circle
# it belongs to is unambiguous from the stem.
for d, spec in enumerate(DUPLEX):
    for i, (lab, colour) in spec['fallback'].items():
        x = X0 + i * DX
        ystem = bb((d, 'top', i)).y1 + 3
        ycirc = ystem + 11 + R
        ax.add_patch(Rectangle((x - 0.9, ystem), 1.8, ycirc - R - ystem,
                               facecolor=colour, edgecolor='none', zorder=2))
        ax.add_patch(Circle((x, ycirc), R, facecolor='none', edgecolor=colour,
                            linewidth=1.3, zorder=3))
        ax.text(x + R + 5, ycirc, lab, ha='left', va='center',
                fontsize=FS_ANN, color=colour, zorder=4)

# measured bracket for the 4 bp spacing between the two marks on AAGCCCG
for d, spec in enumerate(DUPLEX):
    if spec['bracket'] is None:
        continue
    i, j, lab = spec['bracket']
    xa, xb = X0 + i * DX, X0 + j * DX
    ytop = min(bb((d, 'bot', k)).y0 for k in range(len(spec['bot']))) - 16
    ybar = ytop - 13
    ax.add_patch(Rectangle((xa, ybar), xb - xa, 1.6, facecolor=GREY,
                           edgecolor='none', zorder=2))
    for xe in (xa, xb):
        ax.add_patch(Rectangle((xe - 0.8, ybar), 1.6, ytop - ybar,
                               facecolor=GREY, edgecolor='none', zorder=2))
    ax.text((xa + xb) / 2, ybar - 7, lab, ha='center', va='top',
            fontsize=FS_ANN, color='0.15', zorder=4)


if __name__ == "__main__":
    # The duplex base letters A/C/G/T are the content of the schematic; placing
    # them on the mark glyph is correct here, as in Supplementary Figure 1.
    report(save(fig, "SuppFigure2.png", allow_over_data=("A", "C", "G", "T")))
