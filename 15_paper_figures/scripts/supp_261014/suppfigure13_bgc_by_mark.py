# 2026-11-04 D175: この生成器は補足図の統合で役目を終えた。出荷名に書き込むと別の図を潰すため、
# 出力名を SUPERSEDED_ 前置きに変更した。現行の図を作るのは suppfigure_merged_*_261103.py。
#!/usr/bin/env python3
"""補足図 13 — 二次代謝産物遺伝子クラスタの染色体上の位置と、マーク別のメチル化密度。

2026-10-15 の設計変更（著者指示）:
  1. 旧版はトラックを認識配列別（TGGCCGGC のみ）に出していたが、この図の主張は
     「BGC との関連は特に無い」という否定である。読者が見るべきは認識配列ごとの
     内訳ではなく 4mC と 6mA の 2 本なので、マーク別に置き換えた。
     認識配列別の分布は補足図 4 が担当する。
  2. 旧版は全トラックで ax.set_yticks([]) を実行しており、縦軸が読めなかった。
     目盛と軸ラベルを入れた（単位は 50 kb 窓あたりの高信頼サイト数）。
  3. クラスタ位置のパネルを図の主役にし、29 領域すべてに名称を注記した（旧版は 4 件のみ）。
  4. 2026-10-15 第二版: 3 時点の折れ線が重なって読めないとの指摘を受け、
     時点平均の実線＋T1-T3 の最小最大帯に変更した。帯の薄さがそのまま
     「時点間でほとんど動かない」ことを示す。

4mC = TGGCCGGC の 4mC（位置 5）＋ AAGCCCG の 4mC（位置 5）。件数の足し算なので
合算できる。占有率は母数が異なるため合算してはならない（補足図 4 の凡例を参照）。
6mA = AAGCCCG の A1/A2 和集合。

実行: cd scripts && PYTHONPATH=. python3 supp_261014/suppfigure13_bgc_by_mark.py
"""
import json
import pathlib
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

sys.path.insert(0, 'supp_261014')
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, genome          # noqa: E402
from _suppstyle import (BLUE, DATA_DIR, META_GREY, NAR_WIDTH_IN, RED,  # noqa: E402
                        apply_figure_style, figure_panel_letters, report, save,
                        source_note)

N = json.load(open(DATA_DIR / "numbers.json"))
bgc = pd.read_csv(DATA_DIR / "bgc_regions_29.csv")

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))

BIN = 50_000
CORE = (1.50e6, 7.17e6)
ORIC = 4.270310e6
# 2026-11-03 FIG23: 時点の 3 色は使用箇所が無くなっていたので削除（意味色は 4mC/6mA/デュアルのみ）
# 2026-10-22: 発現状態の 3 色のうち silent は "#9e4a4a" で、これは 2026-10-22 まで
# この同じパイプラインで 4mC を指していた赤そのものだった。パネル a が正準の 4mC 赤に
# なった以上、同じ図の中で赤が「4mC」と「silent クラスター」の二義になる。発現状態は
# 順序尺度（silent < low < expressed）なので、マーク色と交わらない単一色相の濃淡に移す。
# 行の位置と y 目盛ラベルの色が識別子を兼ねるので、色は冗長な第 2 チャネルである。
import importlib, sys as _sys  # noqa: E402
_sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
_U = importlib.import_module("00_shared_utils")
EXPR, LOW, SIL = _U.COL_BAR_DARK, _U.COL_BAR_MID, _U.COL_BAR_LIGHT
SUB = {"tgg": "pileup_GCCGGC", "c4": "pileup_AAGCCCG_4mC", "a01": "pileup_AAGCCCG_6mA"}


def motif_keys(g):
    return {
        "tgg": {(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", g)}
               | {(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", g)},
        "c4": {(m.start() + 4, "+") for m in re.finditer("(?=AAGCCCG)", g)}
              | {(m.start() + 2, "-") for m in re.finditer("(?=CGGGCTT)", g)},
        "a01": {(m.start() + i, "+") for m in re.finditer("(?=AAGCCCG)", g) for i in (0, 1)}
               | {(m.start() + 5 + i, "-") for m in re.finditer("(?=CGGGCTT)", g) for i in (0, 1)},
    }


def called(setname, tp, keys):
    """プール規則: 時点内 3 ライブラリを位置・鎖で合算、プール後真深度 >= 30、
    プール後 Nmod/Nvalid >= 0.5。補足図 4・7 と同一の規則。"""
    acc = {}
    for lib in LIBS[tp]:
        d = pd.read_csv(DATA_DIR / SUB[setname] / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, st, nv, nm, nf in d.itertuples(index=False):
            k = (pos, st)
            if k not in keys[setname]:
                continue
            a = acc.setdefault(k, [0, 0, 0])
            a[0] += nv
            a[1] += nm
            a[2] += nf
    return [k[0] for k, v in acc.items()
            if v[0] + v[2] >= DEPTH_MIN and v[0] and v[1] / v[0] >= FRAC_MIN]


def collect():
    g = genome()
    keys = motif_keys(g)
    nb = int(np.ceil(len(g) / BIN))
    series = {}
    for tp in ("T1", "T2", "T3"):
        tgg = called("tgg", tp, keys)
        c4 = called("c4", tp, keys)
        a01 = called("a01", tp, keys)
        assert len(tgg) == N["TGG_called"][tp], f"{tp} TGGCCGGC {len(tgg)}"
        assert len(c4) == N["AAG_C4_called"][tp], f"{tp} AAGCCCG 4mC {len(c4)}"
        assert len(a01) == N["AAG_6mA_union"][tp], f"{tp} AAGCCCG 6mA {len(a01)}"
        series[tp] = {
            "m4": np.bincount(np.asarray(tgg + c4) // BIN, minlength=nb),
            "m6": np.bincount(np.asarray(a01) // BIN, minlength=nb),
            "n4": len(tgg) + len(c4), "n6": len(a01),
        }
    return len(g), nb, series


SHORT = {
    1: "leinamycin", 2: "isorenieratene", 3: "lobophorin CR4", 4: "coelichelin",
    5: "region 5 (unnamed)", 6: "informatipeptin", 7: "flaviolin", 8: "ectoine",
    9: "istamycin", 10: "desferrioxamine", 11: "CDA", 12: "actinorhodin",
    13: "albaflavenone", 14: "prejadomycin", 15: "spore pigment", 16: "kinamycin",
    17: "undecylprodigiosin", 18: "14-hydroxyisochainin", 19: "geosmin",
    20: "enduracidin", 21: "coelimycin P1", 22: "nogalamycin", 23: "SapB",
    24: "hopene", 25: "Bisenarsan", 26: "region 26 (unnamed)", 27: "germicidin",
    28: "dimethylallylindole", 29: "coelibactin",
}


def spread(xs, sep, lo, hi):
    """ラベルの x を最小間隔 sep まで押し広げる（前方・後方の二段階）。"""
    x = np.array(xs, dtype=float)
    for i in range(1, len(x)):
        x[i] = max(x[i], x[i - 1] + sep)
    if x[-1] > hi:
        x[-1] = hi
        for i in range(len(x) - 2, -1, -1):
            x[i] = min(x[i], x[i + 1] - sep)
    x[0] = max(x[0], lo)
    for i in range(1, len(x)):
        x[i] = max(x[i], x[i - 1] + sep)
    return x


def shade(ax, glen):
    for x0, x1, col in [(0, CORE[0], "#f2f2f2"), (CORE[0], CORE[1], "#e6e6e6"),
                        (CORE[1], glen, "#f2f2f2")]:
        ax.add_patch(Rectangle((x0 / 1e6, 0), (x1 - x0) / 1e6, 1,
                               transform=ax.get_xaxis_transform(),
                               fc=col, ec="none", zorder=0))
    ax.axvline(ORIC / 1e6, color="#6b6b6b", lw=0.8, ls=(0, (4, 2)), zorder=1)


def track(ax, x, series, key, nkey, glen, ylab, colour):
    """時点平均の実線と T1-T3 の最小最大帯。折れ線 3 本の重なりを避ける。"""
    shade(ax, glen)
    stack = np.vstack([series[t][key] for t in ("T1", "T2", "T3")])
    lo, hi, mean = stack.min(0), stack.max(0), stack.mean(0)
    ax.fill_between(x, lo, hi, step="mid", color=colour, alpha=0.30, lw=0,
                    zorder=2, label="T1-T3 range")
    ax.step(x, mean, where="mid", lw=0.85, color=colour, zorder=3,
            label="mean of T1, T2, T3")
    ax.set_xlim(0, glen / 1e6)
    ax.set_ylim(0, hi.max() * 1.30)
    ax.set_ylabel(ylab)
    ns = "  ".join(f"{t} {series[t][nkey]:,}" for t in ("T1", "T2", "T3"))
    # [D119 で削除: sites called の件数注記]
    ax.legend(frameon=False, loc="upper left", ncol=2, fontsize=6.0,
              handlelength=1.3, columnspacing=1.2, borderpad=0.0,
              handletextpad=0.5, labelspacing=0.2)


def build():
    glen, nb, series = collect()
    x = (np.arange(nb) * BIN + BIN / 2) / 1e6
    # 2026-10-22: FIG09 の 6 pt 床に合わせて、この図の 5.0 / 5.8 pt のリテラル 3 件を
    # 6.0 pt に引き上げた（領域名の ax_c.text は 5.0 pt だった。引き上げ前の姿は
    # .bak_supptitles_261015/suppfigure13_bgc_by_mark.py L207 に残っている）。
    # その結果ここの幾何が破綻した: 縦書きラベルの伸びは字数 x 文字幅 x フォント
    # サイズなので、6 pt では下段ラベルの先端が上段の基線 5.55 を越えて衝突する
    # （最長 "14-hydroxyisochainin" で約 1.0 inch）。以下の figsize / BASE / set_ylim は
    # そのための作り直しで、フォントサイズ自体はこの箇所では変えていない。
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 7.05))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.0, 2.75], hspace=0.26,
                          left=0.118, right=0.985, top=0.968, bottom=0.145)

    ax_a = fig.add_subplot(gs[0])
    track(ax_a, x, series, "m4", "n4", glen, "4mC sites\nper 50 kb", RED)
    ax_a.set_xticklabels([])

    ax_b = fig.add_subplot(gs[1])
    track(ax_b, x, series, "m6", "n6", glen, "6mA sites\nper 50 kb", BLUE)
    ax_b.set_xticklabels([])

    ax_c = fig.add_subplot(gs[2])
    shade(ax_c, glen)
    rows = [("expressed", EXPR, 2), ("low", LOW, 1), ("silent", SIL, 0)]
    colour_of = {st: c for st, c, _ in rows}
    y_of = {st: y for st, _, y in rows}
    for st, col, y in rows:
        for r in bgc[bgc.status == st].itertuples():
            ax_c.add_patch(Rectangle((r.start / 1e6, y - 0.22),
                                     max((r.end - r.start) / 1e6, 0.030), 0.44,
                                     fc=col, ec="none", zorder=4))
    ax_c.set_yticks([y for _, _, y in rows])
    ax_c.set_yticklabels([f"{st} ({len(bgc[bgc.status == st])})" for st, _, _ in rows])
    for lab, (_, col, _) in zip(ax_c.get_yticklabels(), rows):
        lab.set_color(col)
    ax_c.tick_params(axis="y", length=0)

    # 29 領域すべてに名称を注記する。縦書きにし、最小間隔まで押し広げて
    # 真の位置へは引き出し線で戻す（最小間隔 29 kb = 約 0.6 mm のため）。
    d = bgc.sort_values("start").reset_index(drop=True)
    mid = ((d.start + d.end) / 2 / 1e6).to_numpy()
    # 2 段の千鳥に振る。1 段あたりの隣接間隔が倍になるので、押し広げ量が半分で
    # 済み、引き出し線が短く縦向きのままになる（core 側は 1.6 Mb に 11 領域ある）。
    tier = np.arange(len(mid)) % 2
    lx = np.empty_like(mid)
    for t in (0, 1):
        m = tier == t
        lx[m] = spread(mid[m], 0.300, 0.10, glen / 1e6 - 0.10)
    BASE = {0: 2.62, 1: 7.45}
    for xl, xm, rn, st, t in zip(lx, mid, d.region_number, d.status, tier):
        col = colour_of[st]
        y0 = y_of[st] + 0.24
        yb = BASE[t]
        ax_c.plot([xm, xm, xl], [y0, yb - 0.26, yb], lw=0.4, color=col,
                  alpha=0.8, zorder=3, solid_joinstyle="miter")
        ax_c.text(xl, yb + 0.06, SHORT[rn], rotation=90, ha="center", va="bottom",
                  fontsize=6.0, color=col)

    ax_c.set_xlim(0, glen / 1e6)
    ax_c.set_ylim(-0.55, 12.6)
    for sp in ("left", "right", "top"):
        ax_c.spines[sp].set_visible(False)
    ax_c.set_xlabel("Chromosome position (Mb)")
    for xt, lab, col in [(0.75, "left arm", "#777777"), (3.25, "core", "#3a3a3a"),
                         (7.92, "right arm", "#777777")]:
        ax_c.text(xt, -0.44, lab, ha="center", va="bottom", fontsize=6, color=col)

    figure_panel_letters(fig, [ax_a, ax_b, ax_c], "abc",
                         dx=[0.082, 0.082, 0.082], dy=0.006)
    source_note(fig,
                "Source: deposited modkit pileups (canonical, 2026-10-15) and antiSMASH 7 regions. Site rule: three libraries pooled per\n"
                "timepoint, pooled true depth (Nvalid_cov + Nfail) >= 30, pooled modified fraction >= 50 %. 4mC pools TGGCCGGC position 5\n"
                "with AAGCCCG position 5 (counts, not occupancies); 6mA is the AAGCCCG A1/A2 union. Per-recognition-sequence landscapes: Supplementary Figure 4.",
                y=0.006)
    return fig, series


if __name__ == "__main__":
    fig, series = build()
    # 出荷名は SUPERSEDED_SuppFigure6_pre261104.png（補足図 13 -> 7 の改番、2026-10-24 に出力名を追随）。
    # 旧版は拡張子なしの "SuppFigure13" に書いたうえ、report() が出荷側の
    # SuppFigure13.png（別の図＝AAGCCCG の 5' 延長）を読んでいた。
    out = save(fig, "SUPERSEDED_SuppFigure6_pre261104.png")
    for t in series:
        print(f"  {t}: 4mC {series[t]['n4']:,} sites, 6mA {series[t]['n6']:,} sites")
    report(out)
