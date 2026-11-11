"""Figure 4 panel d — funnel of counts: 129 SAM-dependent MTases -> 3 candidates.

Standalone generator for the panel-d tile of Figure 4.
Environment: epi-rna.  Output: panel_d.png, exactly 1039 x 590 px, transparent.
"""
import sys, importlib
from pathlib import Path

# [FIG18] タイルを再生成できるよう、artifact マーカーをローカルの入力表に解決しておく。
_DATA_DIR = Path(__file__).resolve().parents[2] / 'data_261013'

# project style module (name starts with a digit -> import by string)
sys.path.insert(0, '/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts')
U = importlib.import_module('00_shared_utils')
U.apply_unified_style()

import matplotlib as mpl
mpl.rcParams['savefig.bbox'] = None          # apply_style() would set 'tight'
import matplotlib.pyplot as plt
import pandas as pd
import os as _os
_FD = int(_os.environ.get('FIG_DPI', 300))  # 2026-11-09 解像度上書き

# stage, n  (129, 128, 21, 12, 5, 3) — resolved from the artifact marker
DATA = str(_DATA_DIR / "v7295f465_F4d_funnel.csv")
OUT = str(Path(__file__).with_name('panel_F4c.png'))

df = pd.read_csv(DATA)
stages = df['stage'].tolist()
counts = df['n'].tolist()
final = len(stages) - 1

fig = plt.figure(figsize=(2231 / 300, 590 / 300), dpi=_FD)  # 2026-10-28 行スロットと同比
ax = fig.add_axes([0.0760, 0.025, 0.9140, 0.800])  # 左端をパネル a の軸枠に合わせた  # 2026-10-28 左端をパネル a に合わせた
U.style_axes(ax, frame=True)

XMAX = 180.0
BARH = 0.34
ax.set_xlim(0, XMAX)
ax.set_ylim(5.78, -0.88)          # inverted: stage 0 on top

for i, (name, n) in enumerate(zip(stages, counts)):
    is_final = (i == final)
    col = U.COL_6mA if is_final else U.COL_GRAY
    ax.barh(i, n, height=BARH, color=col, edgecolor='none', zorder=2)
    # stage identity label: above its own bar, left-aligned at the funnel mouth
    # 2026-10-28 著者指示: 図内の文字は 1 色にする（色と字の大きさで強調を
    # 重ねない）。最終段の強調は太字だけに残す。
    ax.text(0, i - 0.46, name, ha='left', va='center',
            fontsize=U.FONT['tick'], color=U.COL_DARK,
            fontweight='bold' if is_final else 'normal', zorder=3)
    # count label: just past the bar end; the top one carries the unit
    lab = f'{n} proteins' if i == 0 else f'{n}'
    ax.text(n + 2.8, i, lab, ha='left', va='center',
            fontsize=U.FONT['tick'], color=U.COL_DARK,
            fontweight='bold' if is_final else 'normal', zorder=3)

# every stage is value-labelled, so a count axis would duplicate the encoding
for side in ('left', 'bottom'):
    ax.spines[side].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# 2026-10-28 著者指示:「グラフ内や周辺にやたらと文字を入れていたり、グラフに
# グレーや黒などの様々なカラーのテキストを入れたりしている点。グラフは本来もっと
# 見やすくシンプルにして、説明は凡例に回すべきです」。説明的なテキストを図から外した。
# 外した文言は凡例に移す（値は value_registry.csv に登録済みのものだけを使う）。
# 削除した説明文（凡例へ移す）:
#   「Three 6mA candidates survive: one 679 aa protein and two BREX-2 PglX」
#   「The two DPPY proteins the floor removes are below the all-gene median TPM」
#   「No 4mC candidate remains」
# 3 色（COL_6mA / COL_DARK / COL_GRAY）を使い分けており、色が実体を指さない使い方だった。

# title starts at x=0.122 so the composer's panel letter keeps a clear 10x6 mm
fig.text(0.122, 0.955, '',
ha='left', va='top', fontsize=U.FONT['panel_title'], color=U.COL_DARK)


# 2026-10-28 著者指示: 四辺の枠を出す。style_axes(frame=True) だけでは
# 後続処理で消える経路があったので、保存直前に明示的に立てる。
for _s in ('top', 'right', 'left', 'bottom'):
    ax.spines[_s].set_visible(True)
    ax.spines[_s].set_edgecolor(U.COL_DARK)
    ax.spines[_s].set_linewidth(0.8)
U.assert_no_text_collisions(fig, 'Fig4 panel d')
fig.savefig(OUT, dpi=_FD, transparent=True, bbox_inches=None, pad_inches=0)

from PIL import Image
size = Image.open(OUT).size
print('size', size)
assert size == (2231, 590), size  # 2026-10-28 行スロットと同比に変更
