"""Figure 1 panel b — growth and antibiotic production across the three
sampling timepoints.

Linear axes, one sub-panel per measurement: bars are the mean of four
independent biological replicate cultures, open circles the individual
replicates, error bars the SD. Brackets are Games-Howell after Welch's
one-way ANOVA.

Tukey HSD was used until 2026-11-02 and was withdrawn because it pools the
error variance across all three timepoints. These assays span three orders
of magnitude and their variance scales with the mean (undecylprodigiosin SD
0.0034 / 0.0119 / 0.6649 at T1/T2/T3; Levene p = 0.019), so the pooled
error is set by T3 and swamps the T1-T2 contrast. Tukey therefore returned
p = 0.528 for undecylprodigiosin 12 h vs 24 h despite a 73-fold difference
(0.0042 -> 0.3073 mg/L), and p = 0.994 for actinorhodin despite a 3.4-fold
difference. A Welch/Games-Howell post-hoc, which uses each pair's own
variances, returns p < 0.001 and p = 0.026 for the same two contrasts. A
log transform fixes undecylprodigiosin equally well but is unavailable for
actinorhodin, one of whose T1 replicates is exactly zero.

This panel replaces the log-scale twin-axis version shipped before
2026-10-15. On that version the right-hand axis was logarithmic mg/L while
the left was linear g/L, so the T2 actinorhodin mean of 0.08 mg/L was drawn
at the height of ~3.5 g/L on the left axis and read as substantial
production. On a linear axis the same value is indistinguishable from zero,
which is what the measurement says.

Renders a 2053 x 590 px transparent tile for the Figure 1 composer.
Run as a cell in the `epi-rna` environment.
"""
import sys, importlib
from pathlib import Path

import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import os as _os
_FD = int(_os.environ.get('FIG_DPI', 300))  # 2026-11-09 解像度上書き

sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # the style helper sets 'tight', which resizes the canvas

_DATA_DIR = Path(__file__).resolve().parents[2] / 'data_261013'
# long-format source: 36 rows = 3 measurements x 3 timepoints x 4 biological replicates
CSV = str(_DATA_DIR / 'v054306f1_F1b_growth.csv')
if len(sys.argv) > 1:
    CSV = sys.argv[1]
OUT = str(Path(__file__).with_suffix('.png'))

df = pd.read_csv(CSV)
TPS = ('T1', 'T2', 'T3')
HRS = ('12 h', '24 h', '50 h')
PANELS = [('dry_cell_weight_g_L', 'Dry cell weight (g/L)', 8.0, [0, 2, 4, 6, 8]),
          ('actinorhodin_mg_L', 'Actinorhodin (mg/L)', 15.0, [0, 5, 10, 15]),
          ('undecylprodigiosin_mg_L', 'Undecylprodigiosin (mg/L)', 8.0, [0, 2, 4, 6, 8])]
PAIRS = [('T1', 'T2'), ('T1', 'T3'), ('T2', 'T3')]


AUTHOR = pd.read_csv(
    '/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
    '106_act_red_dcw/tables/act_red_dcw_stats_author_261102.tsv', sep='\t')


def welch_f(groups):
    """Welch's heteroscedastic one-way F (Welch 1951), as GraphPad Prism."""
    k = len(groups)
    n = np.array([len(g) for g in groups], float)
    m = np.array([g.mean() for g in groups])
    v = np.array([g.var(ddof=1) for g in groups])
    w = n / v
    sw = w.sum()
    mw = (w * m).sum() / sw
    A = (w * (m - mw) ** 2).sum() / (k - 1)
    lam = ((1 - w / sw) ** 2 / (n - 1)).sum()
    B = 1 + 2 * (k - 2) / (k ** 2 - 1) * lam
    return A / B, k - 1, (k ** 2 - 1) / (3 * lam)


fig = plt.figure(figsize=(2053 / 300, 590 / 300), dpi=_FD)
axes = [fig.add_subplot(1, 3, i + 1) for i in range(3)]
fig.subplots_adjust(left=0.075, right=0.995, bottom=0.215, top=0.815, wspace=0.42)

report = []
for ax, (meas, ylab, ymax, yticks) in zip(axes, PANELS):
    sub = df[df.measurement == meas]
    assert len(sub) == 12, f'{meas}: {len(sub)} 行（12 行を期待）'
    means = [sub[sub.timepoint == t].value.mean() for t in TPS]
    sds = [sub[sub.timepoint == t].value.std(ddof=1) for t in TPS]

    grp = {t: sub[sub.timepoint == t].value.values for t in TPS}
    rows = AUTHOR[AUTHOR.measurement == meas]
    assert len(rows) == 3, (meas, len(rows))

    # The author's own GraphPad Prism output is the statistical record
    # (Writing/refs/*_author_261102.csv). Reproduce its inputs before using it.
    lev = stats.levene(*[grp[t] for t in TPS]).pvalue
    assert abs(lev - float(rows.levene_p_author.iloc[0])) < 5e-5, (meas, lev)
    test = rows.test_used.iloc[0]
    if test.startswith('one-way'):
        F, _ = stats.f_oneway(*[grp[t] for t in TPS])
    else:
        F, _, _ = welch_f([grp[t] for t in TPS])
    want = float(rows.omnibus.iloc[0].split('=')[1].split(',')[0])
    assert abs(F - want) / want < 1e-4, (meas, F, want)

    star = {tuple(c.split('_vs_')): s for c, s in zip(rows.comparison, rows.stars)}
    report.append((meas, test, F, star))

    for i, (m, s) in enumerate(zip(means, sds)):
        ax.bar(i, m, width=0.62, color=U.COL_BAR_MID, linewidth=0, zorder=2)
        ax.errorbar(i, m, yerr=s, fmt='none', ecolor=U.COL_DARK,
                    elinewidth=0.9, capsize=2.2, capthick=0.9, zorder=3)
    for i, t in enumerate(TPS):
        v = sub[sub.timepoint == t].value.values
        ax.plot(np.full(len(v), i) + np.linspace(-0.13, 0.13, len(v)), v, 'o',
                ms=3.0, mfc='white', mec=U.COL_DARK, mew=0.7, zorder=4)

    # Games-Howell brackets, stacked above the tallest bar+SD
    top = max(m + s for m, s in zip(means, sds))
    step = ymax * 0.105
    base = top + step * 0.55
    for k, (a, b) in enumerate(PAIRS):
        ia, ib = TPS.index(a), TPS.index(b)
        y = base + step * k
        ax.plot([ia, ia, ib, ib], [y, y + step * 0.13, y + step * 0.13, y],
                lw=0.7, color=U.COL_DARK, zorder=3, clip_on=False)
        ax.text((ia + ib) / 2, y + step * 0.16, star[(a, b)], ha='center',
                va='bottom', fontsize=U.FONT['annot'], color=U.COL_DARK, clip_on=False)

    ax.set_xticks(range(3))
    ax.set_xticklabels(HRS)
    ax.set_xlim(-0.65, 2.65)
    ax.set_ylim(0, ymax)
    ax.set_yticks(yticks)
    ax.set_ylabel(ylab)
    ax.set_xlabel('Time after inoculation (h)', fontsize=U.FONT['annot'],
                  color=U.COL_GRAY, labelpad=3)
    U.style_axes(ax, frame=True)

U.assert_no_text_collisions(fig, 'panel_1b')
fig.savefig(OUT, dpi=_FD, transparent=True, bbox_inches=None, pad_inches=0)
plt.close(fig)

for meas, test, F, star in report:
    print(f'{meas:26s} {test:27s} stat={F:9.4f} | '
          + ' '.join(f'{a}-{b} {star[(a, b)]}' for a, b in PAIRS))
print('wrote', OUT)
