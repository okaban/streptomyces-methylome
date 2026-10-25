# 2026-10-24: 出力名を SuppFigure16.png から SuppFigure10.png へ修正。
# renumber_supp.py が PNG を 16 -> 10 に改番したとき、この生成器が追随して
# いなかった。figure_manifest.tsv はこのスクリプトを SuppFigure10.png の
# 生成器として登録しているので、再実行すると出荷されない名前に書き、
# 補足図 10 は作り直せないままだった。check_provenance.py の
# check_figure_name_coverage が検出。
"""Supplementary Figure 16 | Why the 4mC recognition sequence is scarce at transcription start sites.

Partition-free: uses the experimental TSS set of Jeong et al. (2016) and the
reference sequence only. No methylation call, no gene classification.

(a) GC content metaprofile around TSSs (+/-1 kb, 50 bp bins), genome mean overlaid.
(b) Occurrence density of TGGCCGGC (the 4mC recognition sequence, 87.5% GC) and AAGCCCG
    (dual-mark heptamer, 71.4% GC) in the same bins, each scaled to its own flank mean.
(c) Observed/expected occurrence under a GC-matched empirical expectation, at the TSS
    core and in the flanks, with TSS-level bootstrap 95% CIs.

Inputs : data_261013/s16_tss_sequence_profile.csv, s16_tss_sequence_oe.csv
Output : figures_261013/SuppFigure10.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import importlib.util
_u = Path(__file__).resolve().parents[1] / '00_shared_utils.py'
_spec = importlib.util.spec_from_file_location('shared_utils', _u)
U = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(U)

DATA = Path(__file__).resolve().parents[2] / 'data_261013'
OUT = Path(__file__).resolve().parents[2] / 'figures_261013'


def build(out_path=None):
    prof = pd.read_csv(DATA / 's16_tss_sequence_profile.csv')
    oe = pd.read_csv(DATA / 's16_tss_sequence_oe.csv')
    genome_gc = float(prof.attrs.get('genome_gc', 0)) or 72.12
    n_tss = int(oe.loc[oe.series == 'TGGCCGGC_core', 'n_tss'].iloc[0])

    U.apply_unified_style()
    fig = plt.figure(figsize=(U.mm_to_inch(174), U.mm_to_inch(62)))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.88], wspace=0.42,
                          left=0.072, right=0.985, top=0.90, bottom=0.21)
    axa, axb, axc = (fig.add_subplot(gs[0, i]) for i in range(3))

    # ---- (a) GC content -------------------------------------------------
    axa.axhline(genome_gc, color=U.COL_GRAY, lw=0.9, ls='--', zorder=1)
    axa.plot(prof.offset_bp, prof.gc_percent, color=U.COL_DARK, lw=1.6, zorder=3)
    axa.axvline(0, color=U.COL_GRAY, lw=0.8, ls=':', zorder=1)
    axa.set_xlabel('Position relative to TSS (bp)')
    axa.set_ylabel('GC content (%)')
    axa.set_xlim(-1000, 1000)
    axa.set_xticks([-1000, -500, 0, 500, 1000])
    axa.annotate(f'genome mean {genome_gc:.1f}%', xy=(-960, genome_gc), xytext=(-960, genome_gc + 0.9),
                 fontsize=6.4, color=U.COL_GRAY, ha='left', va='bottom')
    lo = prof.gc_percent.min()
    axa.annotate(f'{lo:.1f}%', xy=(0, lo), xytext=(90, lo - 0.1), fontsize=6.6,
                 color=U.COL_DARK, ha='left', va='center',
                 arrowprops=dict(arrowstyle='-', color=U.COL_DARK, lw=0.7))
    U.style_axes(axa, grid=True, grid_axis='y')

    # ---- (b) motif density ----------------------------------------------
    fl = np.abs(prof.offset_bp) > 600
    C4 = U.MOTIF_COLORS.get('TGGCCGGC', U.COL_BAR_DARK)
    for col, colour, lw, ls, lab in (
            ('tggccggc_per_kb', C4, 1.8, '-', 'TGGCCGGC (87.5% GC)'),
            ('aagcccg_per_kb', U.COL_BOTH, 1.5, '-', 'AAGCCCG (71.4% GC)')):
        y = prof[col] / prof.loc[fl, col].mean()
        axb.plot(prof.offset_bp, y, color=colour, lw=lw, ls=ls, label=lab, zorder=3)
    axb.axhline(1.0, color=U.COL_GRAY, lw=0.9, ls='--', zorder=1)
    axb.axvline(0, color=U.COL_GRAY, lw=0.8, ls=':', zorder=1)
    axb.set_xlabel('Position relative to TSS (bp)')
    axb.set_ylabel('Occurrence density\n(relative to flank mean)')
    axb.set_xlim(-1000, 1000)
    axb.set_xticks([-1000, -500, 0, 500, 1000])
    axb.set_ylim(0, 1.6)
    axb.legend(loc='lower center', frameon=False, fontsize=6.3, handlelength=1.5,
               borderpad=0.1, labelspacing=0.25)
    U.style_axes(axb, grid=True, grid_axis='y')

    # ---- (c) GC-matched O/E ---------------------------------------------
    order = ['TGGCCGGC_core', 'TGGCCGGC_flank', 'AAGCCCG_core', 'AAGCCCG_flank']
    d = oe.set_index('series').loc[order]
    xpos = np.array([0, 0.78, 1.95, 2.73])
    err = np.vstack([d.oe - d.ci_lo, d.ci_hi - d.oe])
    axc.bar(xpos, d.oe, width=0.66, color=U.COL_BAR_MID, edgecolor='none', zorder=3)
    axc.errorbar(xpos, d.oe, yerr=err, fmt='none', ecolor=U.COL_DARK, elinewidth=0.9,
                 capsize=2.2, capthick=0.9, zorder=4)
    axc.axhline(1.0, color=U.COL_DARK, lw=0.9, ls='--', zorder=2)
    axc.set_xticks(xpos)
    axc.set_xticklabels(['core', 'flank'] * 2, fontsize=6.6)
    axc.set_ylabel('Observed / expected\n(GC-matched)')
    axc.set_ylim(0, 1.62)
    axc.set_xlim(-0.55, 3.28)
    for xc, lab in ((0.39, 'TGGCCGGC'), (2.34, 'AAGCCCG')):
        axc.annotate(lab, xy=(xc, -0.135), xycoords=('data', 'axes fraction'),
                     ha='center', va='top', fontsize=7.0, color=U.COL_DARK)
    axc.annotate('*', xy=(0, d.loc['TGGCCGGC_core', 'ci_hi'] + 0.09), ha='center', va='bottom',
                 fontsize=10, color=U.COL_DARK)
    U.style_axes(axc, grid=True, grid_axis='y')

    # [FIG03] パネル記号は組版確定後に図座標の一定位置へ置く。軸ボックス左端からの
    # 一定オフセットなので、y 軸ラベルの幅や行数が違っても記号はずれない。
    LETTER_DX = 0.055          # 軸ボックス左端から左へ（図幅に対する比）
    LETTER_Y = 0.965           # 3 パネル共通の上端
    for ax, lab in ((axa, 'a'), (axb, 'b'), (axc, 'c')):
        x0 = ax.get_position().x0
        fig.text(max(0.004, x0 - LETTER_DX), LETTER_Y, lab,
                 fontsize=U.FONT['panel_letter'], fontweight='bold',
                 ha='left', va='top', color='black')

    # [FIG01] 図全体のタイトルは画像に焼き込まない（Supplementary Figure 16 の見出しは凡例側）。
    # n_tss は凡例本文に書くため、ここでは描画しない。

    out = Path(out_path) if out_path else OUT / 'SuppFigure10.png'
    U.save_figure_raster(fig, out, width_class='full', dpi=300)
    plt.close(fig)
    return out


if __name__ == '__main__':
    print(build())
