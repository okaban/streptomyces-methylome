# -*- coding: utf-8 -*-
"""C22 + D193 の判断材料: Figure 2 再構成案 A / B / C を合成する（検討用）。

既存の合成器 40_main_figures_261013.py と同じ 12 カラム・174 mm のグリッド、
同じパネル記号の打ち方を使う。既存のタイル（panel_3a / panel_2a / panel_2d）は
読むだけで、書き換えない。出力は Writing/fig2_option_{A,B,C}_261108.png。

実行環境: epi-rna
"""
import importlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image

SCRIPTS = Path('/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
sys.path.insert(0, str(SCRIPTS))
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

HERE = Path(__file__).resolve().parent
PAN = SCRIPTS / 'panels_261013'
OUT = Path('/Users/okaban/obsidian/Research/rna-seq/Writing')
DPI, GUTTER_MM = 300, 4

LOGO = PAN / 'panel_3a.png'        # 現 Figure 2 a  配列ロゴ
DENS = PAN / 'panel_2a.png'        # 現 Figure 2 b  染色体分布
PISC = PAN / 'panel_2d.png'        # 現 Figure 2 e  572 文脈

OPTIONS = {
    'A': dict(
        row_heights_mm=[48, 36, 43, 54],
        panels=[('a', 0, 0, 12, HERE / 't_counts_5mC.png'),
                ('b', 1, 0, 12, LOGO),
                ('c', 2, 0, 12, DENS),
                ('d', 3, 0, 4, PISC),
                ('e', 3, 4, 8, HERE / 't_cover_wide.png')]),
    'B': dict(
        row_heights_mm=[48, 36, 43, 54],
        panels=[('a', 0, 0, 12, HERE / 't_counts_plain.png'),
                ('b', 1, 0, 12, LOGO),
                ('c', 2, 0, 12, DENS),
                ('d', 3, 0, 4, PAN / 'panel_2b.png'),   # 現 Figure 2 c をそのまま
                ('e', 3, 4, 4, PISC),
                ('f', 3, 8, 4, HERE / 't_cover_narrow.png')]),
    'C': dict(
        row_heights_mm=[74, 36, 43, 54],
        panels=[('a', 0, 0, 12, HERE / 't_merged.png'),
                ('b', 1, 0, 12, LOGO),
                ('c', 2, 0, 12, DENS),
                ('d', 3, 0, 4, PISC),
                ('e', 3, 4, 8, HERE / 't_cover_wide.png')]),
}


def geom(row_heights_mm, ncol=12, width_mm=174):
    mm = DPI / 25.4
    W = int(width_mm * mm)
    g = int(GUTTER_MM * mm)
    colw = (W - g * (ncol - 1)) // ncol
    rowh = [int(h * mm) for h in row_heights_mm]
    row_y = [sum(rowh[:i]) + g * i for i in range(len(rowh))]
    return W, colw, rowh, row_y, g


def compose(tag):
    s = OPTIONS[tag]
    W, colw, rowh, row_y, g = geom(s['row_heights_mm'])
    H = row_y[-1] + rowh[-1]
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    for L, row, col, span, tile in s['panels']:
        w = colw * span + g * (span - 1)
        h = rowh[row]
        x, y = col * (colw + g), row_y[row]
        ax = fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])
        ax.imshow(Image.open(tile).convert('RGBA'))
        ax.set_anchor('W')          # 既存合成器と同じ: 余白は右に寄せ左端を揃える
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_axis_off()
        fig.text((x + 0.012 * W) / W, 1 - (y + 0.010 * H) / H, L,
                 fontsize=U.FONT['panel_letter'], fontweight='bold',
                 ha='left', va='top', color='black')
    U.assert_no_text_collisions(fig, f'fig2_option_{tag}')
    path = OUT / f'fig2_option_{tag}_261108.png'
    fig.savefig(path, dpi=DPI, bbox_inches=None, pad_inches=0,
                facecolor='white')
    plt.close(fig)
    px = Image.open(path).size
    assert px[0] / DPI * 25.4 <= U.NAR_WIDTH_MM['full'] + 0.5, (tag, px)
    print(f'{tag}: {px[0]}x{px[1]} px = {px[0]/DPI*25.4:.1f} x '
          f'{px[1]/DPI*25.4:.1f} mm, {len(s["panels"])} panels -> {path}')
    return path, px


if __name__ == '__main__':
    for t in ('A', 'B', 'C'):
        compose(t)
