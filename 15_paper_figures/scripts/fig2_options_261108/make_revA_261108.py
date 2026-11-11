# -*- coding: utf-8 -*-
"""案A 改訂版（2026-11-08、著者指示）。Figure 2 を 5 パネルで組み直す。

著者の指摘:
  「案Aのaのグラフを4mC、5mC、6mAにして欲しい。今は4mCを2モチーフ、6mAを1モチーフに
    分けてbarを作っているが、これだと5mCとの比較をする上でフェアではない。実際に
    コールできた数をT1、T2、T3ごとに示すのはどうか？」
  「上記で指示したbarを折れ線で表すのがcです。なので、aの後は現在のcにして、cには
    5mCも入れましょう。まあ、描写はされないと思いますが。」

新しい並び:
  (a) 化学種別の実コール数（4mC / 5mC / 6mA、各 T1-T3）    ← 本スクリプトで新規作成
  (b) 染色体に沿った 100 kb あたり部位密度（5mC 系列を追加）← 本スクリプトで再描画
  (c) 配列ロゴ 2 枚            t_revA_logo.py（panel_3a.py から散文の注記を削除）
  (d) 被覆率 MECE              t_cover_wide.png（既存）
  (e) Pisciotta 572 文脈       本スクリプトで描き直し（縦軸の判読性）

並べ替え（2026-11-08 著者指示）: 似た話題を隣に置く。ロゴ（モチーフの同定）と
MECE 内訳（モチーフで説明できる割合）を隣接させ、先行研究との差分という離れた
話題を最後に置く。

(a) の値は pileup_cache_261015 から正準規則で数え直したもの（panelA_counts_261108.csv）。
4mC 2,909/2,869/2,840 と 6mA 1,947/2,150/2,308 は canonical_coverage_261107.tsv と一致。

出力: Writing/fig2_revA_261108.png（既存の Figure2.png は触らない）
実行環境: epi-rna
"""
import importlib
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

SCRIPTS = Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
sys.path.insert(0, str(SCRIPTS))
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()
mpl.rcParams['savefig.bbox'] = None

HERE = Path(__file__).resolve().parent
DATA = SCRIPTS.parent / 'data_261013'
PAN = SCRIPTS / 'panels_261013'
OUT_DIR = Path('/Users/okaban/obsidian/Research/rna-seq/Writing')
import os as _os
DPI, GUTTER_MM = int(_os.environ.get('FIG_DPI', 300)), 4
F_TP = [0.45, 0.72, 1.0]          # T1 -> T3 の濃淡
TPS = ['T1', 'T2', 'T3']


def shade(col, f):
    r, g, b = mcolors.to_rgb(col)
    return (r + (1 - r) * (1 - f), g + (1 - g) * (1 - f), b + (1 - b) * (1 - f))


# ── (a) 化学種別の実コール数 ─────────────────────────────────────────────────
MARKS = [('4mC', U.COL_4mC, 'cytosines'),
         ('5mC', U.COL_GRAY, 'cytosines'),
         ('6mA', U.COL_6mA, 'adenines')]


def tile_panel_a(name='t_revA_counts.png', px=(2125, 720)):
    df = pd.read_csv(HERE / 'panelA_counts_261108.csv')
    called = {(r.mark, r.timepoint): int(r.called) for r in df.itertuples()}
    meas = {(r.mark, r.timepoint): int(r.measurable) for r in df.itertuples()}
    assert all(called[('5mC', t)] == 0 for t in TPS), called   # 非検出をデータで確認

    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.072, right=0.995, bottom=0.150, top=0.920)

    pitch, note_gap, group_gap = 0.34, 0.32, 0.75
    base = [i * (2 * pitch + note_gap + group_gap) for i in range(3)]
    XMAX = 3500
    annotated = []

    for i, (mark, col, unit) in enumerate(MARKS):
        for j, tp in enumerate(TPS):
            # 2026-11-09 [10]: 時点は灰階調。マークの識別は y 軸ラベルの色が担う。
            v, y, c = called[(mark, tp)], base[i] + j * pitch, U.TIME_RAMP[j]
            if v == 0:
                ax.plot([0], [y], marker='o', ms=3.0, mfc=c, mec=c,
                        clip_on=False, zorder=4)
            else:
                ax.barh(y, v, height=0.28, color=c, linewidth=0, zorder=3)
                t = f'{v:,}'
                ax.annotate(t, xy=(v, y), xytext=(5, 0), textcoords='offset points',
                            ha='left', va='center', fontsize=U.FONT['annot'],
                            color=U.COL_DARK)
                annotated.append(t)
        # 2026-11-09 [10]: 測定可能塩基数の注記は凡例に移した（図中に散文を置かない）。

    # 時点の濃淡キー（1 度だけ、上の空き）
    ykey = base[0] - 0.62
    for j, xv in enumerate([2750, 3000, 3250]):
        ax.plot([xv], [ykey], marker='s', ms=4.4, mfc=U.TIME_RAMP[j],
                mec=U.TIME_RAMP[j], zorder=4, clip_on=False)
        ax.annotate(TPS[j], xy=(xv, ykey), xytext=(5, 0),
                    textcoords='offset points', ha='left', va='center',
                    fontsize=U.FONT['annot'], color=U.COL_DARK,
                    annotation_clip=False)

    ax.set_xlim(0, XMAX)
    ax.set_ylim(base[2] + 2 * pitch + 0.34, base[0] - 0.92)
    ax.set_yticks([b + pitch for b in base])
    ax.set_yticklabels([m[0] for m in MARKS])
    ax.tick_params(axis='y', length=0)
    for lab, (mark, col, _) in zip(ax.get_yticklabels(), MARKS):
        lab.set_color(col if mark != '5mC' else U.COL_DARK)
        lab.set_fontsize(U.FONT['axis_label'])
    ax.set_xticks([0, 1000, 2000, 3000])
    ax.set_xlabel('Sites called (pooled depth \u2265 30, level \u2265 50%)', labelpad=2)
    U.style_axes(ax, frame=True)

    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(annotated))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px, {m[0]: [called[(m[0], t)] for t in TPS] for m in MARKS})


# ── (b) 染色体密度に 5mC 系列を足す ──────────────────────────────────────────
def tile_panel_b(name='t_revA_density.png', px=(2125, 560)):
    """panel_2a.py の作図をそのまま踏襲し、恒常 0 の 5mC 系列を加える。

    5mC は全 87 窓で 0 件（部位が 1 つも無いので密度も 0）。y = 0 の線が軸枠に
    重なると「描いていない」と区別できないので、y 軸の下端を 0 より下に取り、
    ゼロ線が枠から離れて見えるようにする。線そのものは正しく y = 0 に置く。
    """
    df = pd.read_csv(DATA / 'v1c1184c7_F2a_density.csv')
    x = df['window_mb'].values
    y4 = df['sites_4mC_per_100kb'].values
    y6 = df['sites_6mA_per_100kb'].values
    s4 = df['sites_4mC_per_100kb'].rolling(5, center=True, min_periods=1).mean().values
    s6 = df['sites_6mA_per_100kb'].rolling(5, center=True, min_periods=1).mean().values
    m4, m6 = y4.mean(), y6.mean()
    y5 = [0.0] * len(x)           # 5mC: 正準規則で 0 件 -> どの窓でも 0

    CORE_LO, CORE_HI, XMAX, YMAX, YMIN = 1.5, 7.17, 9.72, 58.0, -5.0

    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.058, right=0.997, bottom=0.152, top=0.930)

    ax.fill_between([CORE_LO, CORE_HI], YMIN, YMAX,
                    color=mcolors.to_rgba(U.COL_GRAY, 0.16), linewidth=0, zorder=0)
    ax.plot(x, y4, color=mcolors.to_rgba(U.COL_4mC, 0.38), lw=0.7, zorder=1.5)
    ax.plot(x, y6, color=mcolors.to_rgba(U.COL_6mA, 0.38), lw=0.7, zorder=1.5)
    ax.plot(x, s4, color=U.COL_4mC, lw=1.6, solid_capstyle='round', zorder=2.5)
    ax.plot(x, s6, color=U.COL_6mA, lw=1.6, solid_capstyle='round', zorder=2.5)
    ax.plot(x, y5, color=U.COL_GRAY, lw=1.6, solid_capstyle='round', zorder=2.5)

    ax.set_xlim(-0.08, XMAX)
    ax.set_ylim(YMIN, YMAX)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_yticks([0, 10, 20, 30, 40, 50])
    ax.set_xlabel('Position along chromosome (Mb)', labelpad=2)
    ax.set_ylabel('Sites per 100 kb', labelpad=2)
    U.style_axes(ax, frame=True)
    ax.tick_params(axis='both', pad=1.5, length=2.5, width=0.8)

    FS_ID, FS_AN = U.FONT['axis_label'], U.FONT['annot']
    ax.text((CORE_LO + CORE_HI) / 2, 54.0, 'chromosomal core (1.5\u20137.17 Mb)',
            ha='center', va='center', fontsize=FS_AN,
            color=mcolors.to_rgba(U.COL_DARK, 0.72))
    ax.text(8.78, s4[-1], '4mC', ha='left', va='center', fontsize=FS_ID,
            color=U.COL_4mC)
    ax.text(8.78, s4[-1] - 5.6, f'mean {m4:.0f}', ha='left', va='center',
            fontsize=FS_AN, color=mcolors.to_rgba(U.COL_4mC, 0.85))
    ax.text(8.78, 12.6, '6mA', ha='left', va='center', fontsize=FS_ID,
            color=U.COL_6mA)
    ax.text(8.78, 7.4, f'mean {m6:.1f}', ha='left', va='center',
            fontsize=FS_AN, color=mcolors.to_rgba(U.COL_6mA, 0.85))
    # 5mC の系列は描いてある。他の 2 系列と同じ「名前＋平均」の形で名指しし、
    # 見えない理由（全窓 0）は凡例の文に置く（図に散文を焼き込まない）。
    ax.text(8.78, 1.4, '5mC', ha='left', va='center', fontsize=FS_ID,
            color=U.COL_GRAY)
    ax.text(8.78, -3.6, 'mean 0', ha='left', va='center',
            fontsize=FS_AN, color=U.COL_GRAY)

    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name)
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px)


# ── (e) Pisciotta 572 文脈 — 縦軸が読めるように描き直す ────────────────────
def tile_panel_e(name='t_revA_pisciotta.png', px=(653, 638)):
    """panel_2d.py と同じ入力・同じ分類・同じ色で、左余白と軸名だけを直した版。

    2026-11-08 著者指摘「今の d とグラフの縦軸が重なって縦軸の単位が見えづらい」。
    原因は 2 つ: 軸名 'Share of reported contexts (%)' が 29 字あり、4 カラム枠の
    タイル高いっぱいに伸びて合成器が打つパネル記号の角に届くこと、および左余白が
    軸名＋目盛ラベルに対して狭いこと。軸名を短縮し、左余白を広げ、軸の中心を
    下げて記号の領域から外す。由来・母数は x 軸ラベルに残す（C19 の指示）。
    """
    from matplotlib.colors import to_rgb
    df = pd.read_csv(DATA / 'v45b56046_F2d_pisciotta.csv')
    counts = dict(zip(df['class'], df['n']))
    TOT = int(df['n'].sum())
    F = counts['GCCGGC / GGCCGG family']
    G = counts['GCCCG-containing']
    O = counts['other']
    assert F + G + O == TOT, (F, G, O, TOT)
    pf, pg, po = 100 * F / TOT, 100 * G / TOT, 100 * O / TOT
    p_inscope = pf + pg
    # 2026-11-09 [9][10]: 色相は 4mC の赤だけ。残り 2 区分は灰階調。
    PALE_6mA = U.COL_BAR_MID

    fig = plt.figure(figsize=((px[0] + 0.4) / 300, (px[1] + 0.4) / 300), dpi=DPI)
    ax = fig.add_axes([0.235, 0.232, 0.735, 0.628])

    XC, BW = 0.0, 0.76
    segs = [(pf, U.COL_4mC, 'GCCGGC hexamer\nfamily', 'white'),
            (pg, PALE_6mA, 'GCCCG-containing', 'white'),
            (po, U.COL_BAR_PALE, 'other', U.COL_DARK)]
    inside, bottom = [], 0.0
    for frac, col, lab, tcol in segs:
        ax.bar(XC, frac, width=BW, bottom=bottom, color=col, edgecolor='none',
               zorder=2)
        t = ax.text(XC, bottom + frac / 2, lab, ha='center', va='center',
                    color=tcol, fontsize=U.FONT['annot'], linespacing=1.15,
                    zorder=4)
        inside.append(t.get_text()[:28])
        bottom += frac

    # 2026-11-09 [9]: 括弧注記（85.0% in 4mC motif space）は凡例に移した。

    ax.set_xlim(-0.62, 0.62)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_xticks([XC])
    # 1 行が軸より長いと、棒の真下に中央寄せした途端に y 目盛へ食い込む。
    # 行を短く割り直して軸幅に収める（字詰め・フォント梯子はいじらない）。
    ax.set_xticklabels([f'Pisciotta 2023\nn = {TOT}'])
    ax.tick_params(axis='x', length=0)
    # 29 字 -> 21 字。回転した軸名の長さがタイル高を超えなくなる。
    ax.set_ylabel('Share of contexts (%)', labelpad=3)
    U.style_axes(ax, frame=True)

    out = HERE / name
    fig.savefig(out, dpi=DPI, transparent=True, bbox_inches=None, pad_inches=0)
    U.assert_no_text_collisions(fig, name, allow_over_data=tuple(inside))
    plt.close(fig)
    assert all(abs(a-round(b*DPI/300))<=2 for a,b in zip(Image.open(out).size, px)), Image.open(out).size
    print('wrote', out.name, px, 'pct', round(pf, 2), round(pg, 2), round(po, 2))


# ── 合成 ────────────────────────────────────────────────────────────────────
LAYOUT = dict(
    row_heights_mm=[61, 44, 34, 54],
    panels=[('a', 0, 0, 12, HERE / 't_revA_counts.png'),
            ('b', 1, 0, 12, HERE / 't_revA_density.png'),
            ('c', 2, 0, 12, HERE / 't_revA_logo.png'),
            ('d', 3, 0, 8, HERE / 't_cover_wide.png'),
            ('e', 3, 8, 4, HERE / 't_revA_pisciotta.png')])


def compose(out_name='fig2_revA_261108.png'):
    mm = DPI / 25.4
    W, g, ncol = int(174 * mm), int(GUTTER_MM * mm), 12
    colw = (W - g * (ncol - 1)) // ncol
    rowh = [int(h * mm) for h in LAYOUT['row_heights_mm']]
    row_y = [sum(rowh[:i]) + g * i for i in range(len(rowh))]
    H = row_y[-1] + rowh[-1]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    for L, row, col, span, tile in LAYOUT['panels']:
        w = colw * span + g * (span - 1)
        h = rowh[row]
        x, y = col * (colw + g), row_y[row]
        ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
        ax.imshow(Image.open(tile).convert('RGBA'))
        ax.set_anchor('W')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_axis_off()
        fig.text((x + 0.012 * W) / W, 1 - (y + 0.010 * H) / H, L,
                 fontsize=U.FONT['panel_letter'], fontweight='bold',
                 ha='left', va='top', color='black')
    U.assert_no_text_collisions(fig, 'fig2_revA')
    path = OUT_DIR / out_name
    fig.savefig(path, dpi=DPI, bbox_inches=None, pad_inches=0, facecolor='white')
    plt.close(fig)
    px = Image.open(path).size
    assert px[0] / DPI * 25.4 <= U.NAR_WIDTH_MM['full'] + 0.5, px
    print(f'{px[0]}x{px[1]} px = {px[0]/DPI*25.4:.1f} x {px[1]/DPI*25.4:.1f} mm '
          f'-> {path}')


if __name__ == '__main__':
    tile_panel_a()
    tile_panel_b()
    tile_panel_e()
    compose()
