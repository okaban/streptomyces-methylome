# -*- coding: utf-8 -*-
"""共著者コメント C22 ＋ 著者決定 D193 の判断材料として、Figure 2 の再構成案
A / B / C 用のタイルを本物のデータから描く（2026-11-08、検討用・原稿には未反映）。

入力（すべて既存の正準表。本スクリプトは何も上書きしない）:
  counts  15_paper_figures/data_261013/v22031ca8_F2c_counts.csv
          （現 Figure 2 パネル d の入力そのもの）
  5mC ゲノム全体の分母
          15_paper_figures/scripts/recompute_261015/pileup_cache_261015/*.npz
          正準規則（3 ライブラリ合算・真深度 Nvalid+Nfail>=30・Nvalid>0・
          Nmod/Nvalid>=0.5・refok）で code=2 を数えると 3 時点とも 0 件。
          測定可能シトシンは T1 6,247,708 / T2 6,250,558 / T3 6,250,331。
          （本スクリプトでは fivemc_denominator_261108.csv に書き出した値を読む）
  channels 15_paper_figures/data_261013/v1bc294f3_F2b_channels.csv（現パネル c）
  coverage 11_epigenome_integration/analysis/127_coverage_canonical_261107/
           tables/mece_motif_breakdown_261107.tsv ＋ canonical_coverage_261107.tsv
           （D193 の被覆率図＝artifact 75e1de3f の元データ）

出力: このディレクトリの t_*.png（タイル）。合成は compose_options_261108.py。
実行環境: epi-rna
"""
import importlib
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

SCRIPTS = Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
sys.path.insert(0, str(SCRIPTS))
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None   # 共有スタイルは 'tight' を入れるのでタイルでは外す

HERE = Path(__file__).resolve().parent
DATA = SCRIPTS.parent / 'data_261013'
COVDIR = Path('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis'
              '/127_coverage_canonical_261107/tables')
import os as _os
DPI = int(_os.environ.get('FIG_DPI', 300))


def shade(col, f):
    """パレット定数を白に向けて f だけ混ぜる（新しい hex は作らない）。"""
    r, g, b = mcolors.to_rgb(col)
    return (r + (1 - r) * (1 - f), g + (1 - g) * (1 - f), b + (1 - b) * (1 - f))


PALE_4mC = shade(U.COL_4mC, 0.45)
PALE_6mA = shade(U.COL_6mA, 0.45)
VPALE_6mA = shade(U.COL_6mA, 0.18)
GREY_DARK = tuple(np.array(mcolors.to_rgb(U.COL_GRAY)) * 0.72)
F_TP = [0.45, 0.72, 1.0]          # T1 -> T3 の濃淡

# ── 入力 ────────────────────────────────────────────────────────────────────
counts = pd.read_csv(DATA / 'v22031ca8_F2c_counts.csv')
den = pd.read_csv(HERE / 'fivemc_denominator_261108.csv')   # 自前で再計算した分母
DEN = dict(zip(den['timepoint'], den['measurable_cytosines']))
FIVEMC = dict(zip(den['timepoint'], den['called_5mC']))
assert set(FIVEMC.values()) == {0}, FIVEMC            # 非検出をデータから確認

ch = pd.read_csv(DATA / 'v1bc294f3_F2b_channels.csv')
mece = pd.read_csv(COVDIR / 'mece_motif_breakdown_261107.tsv', sep='\t')
cov = pd.read_csv(COVDIR / 'canonical_coverage_261107.tsv', sep='\t')


# ── 1. 計数パネル（現 Figure 2 d）を全幅で描く ─────────────────────────────
def draw_counts(ax, with_genomewide_5mC):
    """横向きグループ棒。with_genomewide_5mC=True なら 5mC 行をゲノム全体に広げる。"""
    c = counts.set_index('series')[['T1', 'T2', 'T3']].values
    if with_genomewide_5mC:
        lbl = ['TGGCCGGC  4mC', 'AAGCCCG  4mC', 'AAGCCCG  6mA',
               '5mC, genome-wide']
    else:
        lbl = ['TGGCCGGC  4mC', 'AAGCCCG  4mC', 'AAGCCCG  6mA',
               '5mC, any motif position']
    cols = [U.COL_4mC, U.COL_4mC, U.COL_6mA, U.COL_GRAY]

    pitch, gap = 0.30, 0.52
    stride = 2 * pitch + gap
    rows = [[i * stride + j * pitch for j in range(3)] for i in range(4)]
    XMAX = 2250
    annotated = []

    for i in range(4):
        for j in range(3):
            v, y, col = c[i, j], rows[i][j], shade(cols[i], F_TP[j])
            if v == 0:
                ax.plot([0], [y], marker='o', ms=2.8, mfc=col, mec=col,
                        clip_on=False, zorder=4)
            else:
                ax.barh(y, v, height=0.26, color=col, linewidth=0, zorder=3)
                if j == 0:
                    # 現パネル d と同じ約束: 値は T1 の棒にだけ置く（3 本すべてに
                    # 置くと 6.5 pt のラベルが縦に重なる）
                    t = f'{v:,}'
                    ax.annotate(t, xy=(max(c[i]), y), xytext=(5, 0),
                                textcoords='offset points',
                                ha='left', va='center', fontsize=U.FONT['annot'],
                                color=shade(cols[i], 0.85))
                    annotated.append(t)

    # 5mC 行の注記: 非検出は「見えない棒」になるので、分母を文字で置く
    if with_genomewide_5mC:
        note = (f'0 called  (of {DEN["T1"]:,} / {DEN["T2"]:,} / {DEN["T3"]:,}\n'
                'measurable cytosines, T1 / T2 / T3)')
        ax.annotate(note, xy=(0, rows[3][1]), xytext=(14, 0),
                    textcoords='offset points', ha='left', va='center',
                    fontsize=U.FONT['annot'], color=GREY_DARK, linespacing=1.25)
    else:
        ax.annotate('n.d.', xy=(0, rows[3][1]), xytext=(14, 0),
                    textcoords='offset points', ha='left', va='center',
                    fontsize=U.FONT['annot'], color=U.COL_GRAY)

    # 時点の濃淡キー（1 度だけ、右上の空き）
    ykey = rows[0][0] - 0.62
    for j, (xv, lab) in enumerate(zip([1700, 1880, 2060], ['T1', 'T2', 'T3'])):
        ax.plot([xv], [ykey], marker='s', ms=4.4,
                mfc=shade(U.COL_4mC, F_TP[j]), mec=shade(U.COL_4mC, F_TP[j]),
                zorder=4, clip_on=False)
        ax.annotate(lab, xy=(xv, ykey), xytext=(5, 0),
                    textcoords='offset points', ha='left', va='center',
                    fontsize=U.FONT['annot'], color=U.COL_DARK, annotation_clip=False)

    ax.set_xlim(0, XMAX)
    ax.set_ylim(rows[3][2] + 0.30, rows[0][0] - 0.90)
    ax.set_yticks([r[1] for r in rows])
    ax.set_yticklabels(lbl)
    ax.tick_params(axis='y', length=0)
    ax.set_xticks([0, 500, 1000, 1500, 2000])
    ax.set_xlabel('Sites called (pooled depth \u2265 30, level \u2265 50%)', labelpad=2)
    U.style_axes(ax, frame=True)
    ax.spines['left'].set_visible(True)
    return annotated


def tile_counts(name, with_genomewide_5mC, px=(2125, 560)):
    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.168, right=0.995, bottom=0.190, top=0.905)
    ann = draw_counts(ax, with_genomewide_5mC)
    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(ann))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px)


# ── 2. チャネル分布（現 Figure 2 c）を全幅で描く ────────────────────────────
def draw_channels(ax):
    a = ch['level_4mC_channel_pct'].to_numpy(float)
    b = ch['level_5mC_channel_pct'].to_numpy(float)
    n = len(ch)
    med_a, med_b = np.median(a), np.median(b)
    bins = np.arange(-1.25, 102.5, 2.5)
    cent = bins[:-1] + 1.25
    H = 0.30

    def ridge(vals, base, fill, edge, med):
        h, _ = np.histogram(vals, bins=bins)
        h = h / h.max() * H
        ax.fill_between(cent, base, base + h, step='mid', facecolor=fill,
                        edgecolor='none', alpha=0.92, zorder=2)
        nz = np.flatnonzero(h > 0)
        i0, i1 = max(nz[0] - 1, 0), min(nz[-1] + 1, len(h) - 1)
        ax.step(cent[i0:i1 + 1], base + h[i0:i1 + 1], where='mid',
                color=edge, linewidth=0.7, zorder=3)
        ax.plot([med, med], [base, base + H * 0.98], color=U.COL_DARK,
                linewidth=1.1, solid_capstyle='butt', zorder=4)

    ridge(b, 0.00, U.COL_GRAY, GREY_DARK, med_b)
    ridge(a, 0.78, U.COL_4mC, U.COL_4mC, med_a)
    # 系列名は尾根の「左の裾」（密度がほぼ 0 の領域）に置く。右上は 4mC の
    # 尾根本体があり、そこに重ねると contrast 規約に触れる。
    ax.text(2, 1.04, '4mC channel', color=U.COL_4mC,
            fontsize=U.FONT['annot'] + 0.5, ha='left', va='bottom')
    ax.text(8, 0.20, '5mC channel', color=GREY_DARK,
            fontsize=U.FONT['annot'] + 0.5, ha='left', va='bottom')
    # 中央値ラベルは尾根の本体に乗せない（同系色どうしで判読できなくなる）。
    # 尾根の頂（base + H）の上に逃がす。
    ax.annotate(f'median {med_a:.0f}%', xy=(med_a, 0.78 + 0.30), xytext=(4, 2),
                textcoords='offset points', color=U.COL_4mC,
                fontsize=U.FONT['annot'], ha='left', va='bottom')
    ax.text(8, 0.04, f'median {med_b:.0f}%', color=GREY_DARK,
            fontsize=U.FONT['annot'], ha='left', va='bottom')
    ax.text(60, 0.52, f'n = {n:,} TGGCCGGC positions, same reads',
            color=U.COL_DARK, fontsize=U.FONT['annot'], ha='center', va='center')
    ax.set_xlim(-4, 103)
    ax.set_ylim(-0.04, 1.26)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([])
    ax.set_xlabel('Methylated reads per site (%)', labelpad=2)
    U.style_axes(ax, frame=True)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    return med_a, med_b, n


# ── 3. 案 C 用: 計数＋チャネルの 2 段パネル ────────────────────────────────
def tile_merged(name='t_merged.png', px=(2125, 880)):
    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax1 = fig.add_axes([0.168, 0.580, 0.827, 0.330])
    ax2 = fig.add_axes([0.168, 0.115, 0.827, 0.250])
    ann = draw_counts(ax1, with_genomewide_5mC=False)
    draw_channels(ax2)
    # 上段の 5mC 行と下段の 5mC 尾根を「同じゼロ」で結ぶ破線（下段のゼロ基準線）
    ax2.axvline(0, color=GREY_DARK, linewidth=0.6, linestyle=(0, (2, 2)),
                zorder=1)
    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(ann))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px)


def tile_channels_small(name='t_channels_small.png', px=(677, 661)):
    """案 B 用: 現パネル c をそのまま使う（既存タイルをコピーせず再描画）。"""
    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.035, right=0.985, bottom=0.175, top=0.930)
    draw_channels(ax)
    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    plt.close(fig)
    print('wrote', out.name, px)


# ── 4. D193 の被覆率（MECE 内訳）パネル ────────────────────────────────────
CAT = [('TGGCCGGC C5', U.COL_4mC, 'white'),
       ('AAGCCCG C5', PALE_4mC, U.COL_DARK),
       ('AAGCCCG A1', U.COL_6mA, 'white'),
       ('AAGCCCG A2', PALE_6mA, U.COL_DARK),
       ('candidate motifs (60)', VPALE_6mA, U.COL_DARK),
       ('unassigned', U.COL_GRAY, U.COL_DARK)]


def draw_coverage(ax, label_fs, legend_fs, inside_min=9.0, legend_loc='right',
                  ylab='Share of called sites (%)', ylab_fs=None):
    order = [('4mC', 'T1'), ('4mC', 'T2'), ('4mC', 'T3'),
             ('6mA', 'T1'), ('6mA', 'T2'), ('6mA', 'T3')]
    xs = [0, 1, 2, 3.7, 4.7, 5.7]
    totals = cov.set_index(['mark', 'timepoint'])['n_sites'].to_dict()
    inside = []
    for x, (mark, tp) in zip(xs, order):
        sub = mece[(mece['mark'] == mark) & (mece['timepoint'] == tp)]
        bottom = 0.0
        for cname, col, tcol in CAT:
            row = sub[sub['category'] == cname]
            if row.empty:
                continue
            pct = float(row['pct'].iloc[0])
            if pct == 0:
                continue
            ax.bar(x, pct, width=0.76, bottom=bottom, color=col,
                   edgecolor='none', zorder=2)
            if pct >= inside_min:
                t = ax.text(x, bottom + pct / 2, f'{pct:.0f}', ha='center',
                            va='center', color=tcol, fontsize=label_fs, zorder=4)
                inside.append(t.get_text())
            bottom += pct
        ax.annotate(f'{totals[(mark, tp)]:,}', xy=(x, 100), xytext=(0, 3),
                    textcoords='offset points', ha='center', va='bottom',
                    fontsize=label_fs, color=U.COL_GRAY)
    ax.set_xticks(xs)
    ax.set_xticklabels(['T1', 'T2', 'T3', 'T1', 'T2', 'T3'])
    ax.tick_params(axis='x', length=0)
    for xc, lab, col in ((1, '4mC', U.COL_4mC), (4.7, '6mA', U.COL_6mA)):
        ax.annotate(lab, xy=(xc, 0), xytext=(0, -17), textcoords='offset points',
                    ha='center', va='top', fontsize=U.FONT['axis_label'],
                    color=col, annotation_clip=False)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel(ylab, labelpad=2,
                  fontsize=ylab_fs or U.FONT['axis_label'])
    U.style_axes(ax, frame=True)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec='none') for _, c, _ in CAT]
    labels = [n for n, _, _ in CAT]
    if legend_loc == 'right':
        ax.set_xlim(-0.65, 6.3)
        ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.02, 0.5),
                  frameon=False, fontsize=legend_fs, handlelength=1.0,
                  handleheight=0.9, labelspacing=0.45, borderpad=0.0)
    else:
        ax.set_xlim(-0.65, 6.3)
        ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.46),
                  ncol=2, frameon=False, fontsize=legend_fs, handlelength=1.0,
                  handleheight=0.9, labelspacing=0.35, columnspacing=0.6,
                  borderpad=0.0)
    return inside


def tile_coverage_wide(name='t_cover_wide.png', px=(1353, 615)):
    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.105, right=0.655, bottom=0.185, top=0.900)
    inside = draw_coverage(ax, U.FONT['annot'], U.FONT['legend'], legend_loc='right')
    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(inside))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px)


def tile_coverage_narrow(name='t_cover_narrow.png', px=(653, 638)):
    """案 B 用: 4 カラム枠に押し込んだ版（凡例を下 2 段に逃がす）。"""
    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.175, right=0.985, bottom=0.440, top=0.900)
    # 4 カラム枠では y 軸名が枠より長くなり、合成器が打つパネル記号の角に
    # 届いてしまう。この案だけ軸名を短縮する（案 A / C では短縮不要）。
    inside = draw_coverage(ax, 5.8, 5.4, inside_min=12.0, legend_loc='below',
                           ylab='Share of sites (%)', ylab_fs=7.0)
    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(inside))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px)


if __name__ == '__main__':
    tile_counts('t_counts_plain.png', with_genomewide_5mC=False)
    tile_counts('t_counts_5mC.png', with_genomewide_5mC=True)
    tile_merged()
    tile_channels_small()
    tile_coverage_wide()
    tile_coverage_narrow()
