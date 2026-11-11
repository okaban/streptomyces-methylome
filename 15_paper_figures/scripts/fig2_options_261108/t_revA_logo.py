"""Figure 3, panel a — sequence logos of the high-confidence 4mC and 6mA positions.

Message: each recognition sequence, shown as the base composition around its own
mark rather than asserted as a drawing.

The two rows are drawn on DIFFERENT site sets, and the difference is the point:
  4mC  - every high-confidence position.  4mC calling is confined to the two
         recognition sequences (3,067 of 3,068), so the unconditioned set is
         informative; it is a 64:36 MIXTURE of TGGCCGGC and AAGCCCG, which is why
         offsets -4, -3 and +1 carry two letters at about 1 bit.  Those offsets
         are two-state because two sequences are pooled, NOT because either
         sequence is degenerate.
  6mA  - only the AAGCCCG-context positions (rebuilt 2026-10-22 on author
         instruction).  On every high-confidence 6mA position the logo resolves
         nothing: 82% of that set is the unassigned adenine fraction, so the logo
         reports this genome's base composition (0.16-0.53 bits per flanking
         offset) and MEME on the same set returns a degenerate GC-rich 9-mer with
         no AAGCCCG.  Conditioning on the recognition sequence is therefore a
         PRESENTATION of the assignment, not evidence for it; the evidence is the
         rate contrast (15.59% of AAGCCCG adenine positions modified at T1 vs
         0.0634% of motif-external reference adenines, 246-fold), which the legend
         carries.  The conditioning and its denominator are printed in the panel.

DATA INPUT (recomputed at render time, nothing hard-coded):
  recompute_261015/canonical_sites_genomewide_261015.tsv  — the canonical site
      list (timepoint / pos / strand / mod).  `pos` is a 0-based genome index,
      the convention Writing/check_pileup_refbase.py enforces.
  11_epigenome_integration/data/NC_003888.3.fna           — the reference.
The site set is OPTION A as approved by the author on 2026-10-20: every
high-confidence position of each chemistry, pooled over T1-T3 (not split by
timepoint) and not filtered on occupancy.  Flanks are taken strand-aware: a
minus-strand site contributes the reverse complement of its window, so offset 0
is the modified base and offsets increase 5'->3' along the modified strand.

Letter heights are the standard sequence-logo definition with Schneider's
small-sample correction: R_i = 2 - H_i - e_n with e_n = 3 / (2 ln2 n), and the
height of base b at position i is f_b,i * R_i.  Heights are written to
data_261013/fig3a_logo_261022.tsv so the registry can cite them.

Tile for the composer's 12-column grid: EXACT 1039 x 685 px, transparent, and no
bbox_inches='tight' (tight re-crops the canvas and breaks the grid).  From
2026-11-08 the LETTERS encode the BASE, using logomaker's 'colorblind_safe'
palette, which is the field standard for sequence logos (MEME and the rest).
These two panels are therefore the documented exception to the figure-set rule
that colour encodes the chemistry (FIG22 / FIG23): the chemistry is carried by
the identity block, which keeps U.COL_4mC / U.COL_6mA, and by the triangle that
marks the modified position, which is drawn in a neutral grey so that it cannot
be read as a base.

Run:  python panel_3a.py        (conda env: epi-rna)
Out:  panel_3a.png
"""
import csv
import importlib
import math
import sys
from collections import Counter
from pathlib import Path
import os as _os
_FD = int(_os.environ.get('FIG_DPI', 300))  # 2026-11-09 解像度上書き

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')   # module name starts with a digit

import logomaker
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # the style helper sets 'tight', which resizes

# ------------------------------------------------------------------- input paths
ROOT = Path('/Users/okaban/bioinfo/rna-seq')
SITES = ROOT / '15_paper_figures/scripts/recompute_261015/canonical_sites_genomewide_261015.tsv'
GENOME = ROOT / '11_epigenome_integration/data/NC_003888.3.fna'
TABLE_OUT = Path(__file__).resolve().parent / 'revA_logo_heights_261108.tsv'

LO, HI = -6, 6                 # offsets about the modified base, inclusive
BASES = ('A', 'C', 'G', 'T')
COMP = str.maketrans('ACGT', 'TGCA')
# Where the modified base sits inside its recognition sequence (LOCKED CANON).
# 2026-10-24: 6mA 行は offset 0 と +1 の両方にマーカーを置くようになったので、
# パネル内ラベルも「position 1」から両位置を名指す表記へ改める。パネルが 2 位置を
# 示しているのにラベルが 1 位置しか名指さない状態は、著者が指摘した不整合と同型。
MARK_POSITION = {'4mC': 'position 5', '6mA': 'positions 1 and 2'}


def reference():
    with open(GENOME) as fh:
        g = ''.join(l.strip() for l in fh if not l.startswith('>')).upper()
    assert len(g) == 8667507, len(g)
    return g


def windows(genome):
    """Strand-oriented flanks of every high-confidence position, by chemistry.

    The reference base is re-checked here rather than trusted: a 4mC call must
    sit on C (+) or G (-) and a 6mA call on A (+) or T (-).  A pileup can carry
    a call to a position whose reference base is neither (the 2026-10-15 defect),
    and such a window would contribute a bogus centre letter to the logo.
    """
    expect = {'4mC': ('C', 'G'), '6mA': ('A', 'T')}
    seen, out = set(), {'4mC': [], '6mA': []}
    with open(SITES) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            p, s, m = int(r['pos']), r['strand'], r['mod']
            if (p, s, m) in seen:          # pooled over T1-T3: count each position once
                continue
            seen.add((p, s, m))
            assert genome[p] == (expect[m][0] if s == '+' else expect[m][1]), (p, s, m)
            if p + LO < 0 or p + HI >= len(genome):
                continue
            w = genome[p + LO: p + HI + 1]
            if s == '-':
                w = w.translate(COMP)[::-1]
            if 'N' in w:
                continue
            out[m].append(w)
    return out


def heights(wins):
    """Sequence-logo letter heights, small-sample corrected."""
    n = len(wins)
    e_n = 3.0 / (2.0 * math.log(2) * n)
    rows = {}
    for j, off in enumerate(range(LO, HI + 1)):
        c = Counter(w[j] for w in wins)
        f = {b: c.get(b, 0) / n for b in BASES}
        H = -sum(p * math.log2(p) for p in f.values() if p > 0)
        R = max(0.0, 2.0 - H - e_n)
        rows[off] = {b: f[b] * R for b in BASES}
    return pd.DataFrame.from_dict(rows, orient='index')[list(BASES)], n


def context_of(w, mod):
    """Which recognition sequence this window's modified base sits in.

    4mC sits at position 5 of BOTH recognition sequences, so its site set is a
    mixture; 6mA sits at A1 (position 1) or A2 (position 2) of AAGCCCG.
    """
    sl = lambda lo, hi: w[lo - LO: hi - LO + 1]
    if mod == '4mC':
        return ('TGGCCGGC' if sl(-4, 3) == 'TGGCCGGC'
                else 'AAGCCCG' if sl(-4, 2) == 'AAGCCCG' else 'unassigned')
    return ('AAGCCCG' if sl(0, 6) == 'AAGCCCG' or sl(-1, 5) == 'AAGCCCG'
            else 'unassigned')


# ---------------------------------------------------------------------- canvas
W, H = 2125, 400                  # px, full-width tile (C10, 2026-11-02): the two
                                  # logos sit side by side so the recognition
                                  # sequences are shown before any panel uses them
FS_ANN = U.FONT['annot']          # 6.5 - in-panel identity block
FS_TICK = U.FONT['tick']          # 7   - tick labels
FS_LAB = U.FONT['axis_label']     # 7.5 - axis labels
YMAX = 2.45                       # 2 bits of logo plus headroom for the mark symbol

# +0.3 px of slack: matplotlib floors the canvas, and the tile must be EXACT
fig = plt.figure(figsize=((W + 0.3) / 300, (H + 0.3) / 300), dpi=_FD)

GEN = reference()
WINS = windows(GEN)

# Rows: 4mC above, 6mA below.  Row 1 stops at y = 0.885 so the top-left
# 118 x 71 px stays clear for the panel letter the composer places there.
#
# `subset` is the site set each row is drawn on, and the two rows differ:
#   4mC  - every high-confidence position (option A).  4mC calling is confined to
#          the two recognition sequences (3,067 of 3,068 positions), so the whole
#          set is informative even though it is a mixture of the two.
#   6mA  - only the AAGCCCG-context positions.  Drawn on every high-confidence 6mA
#          position the logo resolves nothing, because 82% of that set is the
#          unassigned adenine fraction and the logo then reports this genome's
#          base composition.  Conditioning on the recognition sequence is NOT the
#          evidence that AAGCCCG is the recognition sequence -- that rests on the
#          rate contrast quoted in the legend -- so the conditioning and its
#          denominator are printed inside the panel where a reader cannot miss
#          them.  `anchor` puts each identity block over its own logo's empty
#          flank: the 4mC logo decays to the 3' side, the conditioned 6mA logo
#          to the 5' side.
# 2026-11-08 著者決定: ロゴの文字は塩基ごとに塗る（MEME ほか配列ロゴの分野標準）。
# 2026-11-09 著者決定: 配色は Option 1 "Manuscript harmonic" に確定。
#   出所 Writing/logo_palette_options_261109.md（比較シート .png の Option 1）。
#   採用理由: 本稿の配色トークン族と調和する（A は COL_GREEN から ΔE00 5.1、
#   T は COL_DARK と同一色）。二色覚での最悪 ΔE00 17.0、グレースケール最小 ΔL* 14.8 で、
#   どちらも直前に使っていた logomaker 'colorblind_safe'（12.1 / 0.4）を上回る。
#   colorblind_safe は C と G の ΔL* が 0.4 しかなく白黒印刷で潰れていた。
# FIG23-exempt: 配列ロゴは A/C/G/T の識別そのものが情報で、グレー階調では成立しない。
# FIG22-exempt: この 2 パネルに限り文字色はマーク（4mC/6mA）ではなく塩基を符号化する。
#   マークの区別は (1) パネル内の文字ラベル（COL_4mC / COL_6mA のまま）と
#   (2) 修飾位置の三角で保つ。
# FIG23-exempt: 配列ロゴの 4 塩基配色（下 4 行）。宣言は検査窓が直前 3 行なので
# 各行にも置いてある。T は U.COL_DARK と同一値だが、4 色を 1 つの表として
# 読めるよう出所（Option 1）の表記のまま並べる。
BASE_COLOURS = {'A': '#2F5C45',   # FIG23-exempt: 塩基 A（深緑、COL_GREEN から ΔE00 5.1）
                'C': '#8CA8C0',   # FIG23-exempt: 塩基 C（スレート青）
                'G': '#A5743F',   # FIG23-exempt: 塩基 G（ブロンズ）
                'T': '#22282E'}   # FIG23-exempt: 塩基 T（近黒、COL_DARK と同一）

# 三角は無彩色の中間グレー。黒にすると T（#22282E の近黒）と、マーク色にすると
# 6mA の青が C のスレート青と混同しうるため、4 塩基のどの色相とも重ならない灰を選ぶ。
# Option 1 の 4 色はいずれも彩度を持つので、無彩色であること自体が区別になる。
TRI_COLOUR = '#606060'

ROWS = [('4mC', U.COL_4mC, [0.085, 0.300, 0.380, 0.555], None, 'right'),
        ('6mA', U.COL_6mA, [0.585, 0.300, 0.380, 0.555], 'AAGCCCG', 'left')]

axes, counts, matrices, drawn = [], {}, {}, {}
for mod, colour, box, subset, anchor in ROWS:
    counts[mod] = len(WINS[mod])
    use = ([w for w in WINS[mod] if context_of(w, mod) == subset]
           if subset else WINS[mod])
    mat, n = heights(use)
    drawn[mod], matrices[mod] = n, mat
    ax = fig.add_axes(box)
    axes.append(ax)
    logomaker.Logo(mat, ax=ax, color_scheme=BASE_COLOURS, font_name='Arial',
                   stack_order='big_on_top', vpad=0.02, show_spines=False)
    U.style_axes(ax, frame=True)
    ax.set_xlim(LO - 0.6, HI + 0.6)
    ax.set_ylim(0, YMAX)
    ax.set_yticks([0, 1, 2])
    ax.set_xticks(range(LO, HI + 1))
    ax.set_xticklabels([f'{o:+d}'.replace('+0', '0') for o in range(LO, HI + 1)],
                       fontsize=FS_TICK)
    ax.tick_params(labelsize=FS_TICK, length=2, width=0.6, pad=1.5)

    # The modified base, marked where it is rather than named in prose: a filled
    # triangle over offset 0, in the colour that encodes the chemistry.
    ax.plot([0], [YMAX - 0.12], marker='v', markersize=3.6, color=TRI_COLOUR,
            clip_on=False, zorder=5)

    # 2026-10-24: 6mA の行は offset +1（AAGCCCG の 2 番目の A）にも開いた三角を置く。
    # 著者の指摘で判明した不整合への対応。本文 R3 は「position 2 carrying 6mA in
    # 6.97%, 7.50% and 8.85% as fallback」と報告し、Figure 3b の凡例は
    # 「A1 by 5' base, A2」、補足図 2c の凡例は「6mA at A1/A2」と書いているのに、
    # このロゴだけが offset 0 にしかマーカーを置いていなかった。塗りつぶしの三角が
    # 主たる修飾位置、開いた三角が 5' G に塞がれたときのフォールバック位置を表す。
    if mod == '6mA':
        ax.plot([1], [YMAX - 0.12], marker='v', markersize=3.0, markerfacecolor='none',
                markeredgecolor=TRI_COLOUR, markeredgewidth=0.7, clip_on=False, zorder=5)

    # Identity block over the uninformative flank of this row's own logo, where
    # no letter reaches; placement is checked by assert_no_text_collisions below.
    # For a conditioned row the block states the conditioning AND its denominator,
    # so the selection cannot be read as a result.
    if subset:
        # 2026-10-24: 位置の表記が「positions 1 and 2」に伸びたぶんを行分割で吸収する。
        # 1 行に残すと標準の figsize で assert_no_text_collisions が落ち、パネル単体の
        # 再生成ができなくなる（合成図は幅が広いので通ってしまい、事故が隠れる）。
        label = (f'{mod} at {subset}\n{MARK_POSITION[mod]}\n'
                 f'n = {n:,} of {counts[mod]:,} ({100 * n / counts[mod]:.2f}%)')
    else:
        label = f'{mod}, {MARK_POSITION[mod]}\nn = {n:,} positions'
    ax.text(LO - 0.45 if anchor == 'left' else HI + 0.45, YMAX - 0.04, label,
            ha=anchor, va='top', fontsize=FS_ANN, color=colour,
            linespacing=1.45, zorder=4)

for _ax in axes:
    _ax.set_xlabel('Offset from the modified base (bp)', fontsize=FS_LAB,
                   labelpad=2)


# 2026-11-08 著者ルール「図の中に散文の説明文を焼き込まない。説明は凡例に書く」。
# 旧版がここで描いていた 2 文（6mA ロゴの条件づけの説明）は削除した。
# 同じ内容は Figure 2 凡例の (c) に文として置く。パネル内に残るのは
# 系列名・修飾位置・n と百分率という短い値表示のみ。

# 2026-11-02 (C10, side-by-side layout): one y-label per logo, placed as figure
# text just left of each axes. set_ylabel() puts the rotated label outside the
# 2125 x 430 canvas at this tile height, which the collision guard rejects, and
# 'Information content (bits)' is itself too long to stand rotated on a 430 px
# canvas (392 px of glyph against 217 px of axes), so it is shortened here.
for _box in (ROWS[0][2], ROWS[1][2]):
    fig.text(_box[0] - 0.034, _box[1] + _box[3] / 2,
             'Information (bits)', rotation=90, ha='left', va='center',
             fontsize=FS_LAB, color=U.COL_DARK)

# ---------------------------------------------- recognition-context breakdown
# Written for BOTH chemistries whatever each row is drawn on, because the two
# numbers a reader needs are exactly these: the 4mC row is a MIXTURE of the two
# recognition sequences, and the 6mA row is conditioned because the unconditioned
# set is dominated by the mark's unassigned fraction.
CTX = {mod: Counter(context_of(w, mod) for w in WINS[mod])
       for mod, *_ in ROWS}

# ------------------------------------------------- provenance table, then PNG
TABLE_OUT.parent.mkdir(parents=True, exist_ok=True)
with open(TABLE_OUT, 'w', newline='') as fh:
    wr = csv.writer(fh, delimiter='\t')
    wr.writerow(['block', 'mod', 'n_available', 'n_drawn', 'key', 'value']
                + [f'height_{b}' for b in BASES])
    for mod, _c, _b, subset, _a in ROWS:
        wr.writerow(['site_set', mod, counts[mod], drawn[mod],
                     'subset', subset or 'all high-confidence positions']
                    + [''] * len(BASES))
        for off in range(LO, HI + 1):
            row = matrices[mod].loc[off]
            wr.writerow(['logo_heights', mod, counts[mod], drawn[mod],
                         f'offset_{off:+d}', f'{row.sum():.4f}']
                        + [f'{row[b]:.4f}' for b in BASES])
        for k, v in CTX[mod].most_common():
            wr.writerow(['context_counts', mod, counts[mod], drawn[mod], k,
                         f'{v} ({100 * v / counts[mod]:.2f}%)'] + [''] * len(BASES))

U.assert_no_text_collisions(fig, 'Fig2 panel a (recognition-sequence logos)')
out = str(Path(__file__).with_suffix('.png'))
fig.savefig(out, dpi=_FD, transparent=True, bbox_inches=None, pad_inches=0)

sz = Image.open(out).size
assert all(abs(a-round(b*_FD/300))<=2 for a,b in zip(sz, (W, H))), sz
print('saved', out, sz)
print('site sets (available -> drawn):')
for mod, _c, _b, subset, _a in ROWS:
    print(f'   {mod}: {counts[mod]:,} -> {drawn[mod]:,}'
          f'  [{subset or "all high-confidence positions"}]')
print('consensus per offset (bits on the DRAWN set):')
for mod, *_ in ROWS:
    m = matrices[mod]
    print('  ', mod, ' '.join(f'{o:+d}:{m.loc[o].idxmax()}{m.loc[o].sum():.2f}'
                              for o in range(LO, HI + 1)))
print('recognition context of the AVAILABLE set:')
for mod, *_ in ROWS:
    print('  ', mod, dict(CTX[mod]))
