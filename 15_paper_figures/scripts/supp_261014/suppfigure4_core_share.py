# 2026-11-04 D175: この生成器は補足図の統合で役目を終えた。出荷名に書き込むと別の図を潰すため、
# 出力名を SUPERSEDED_ 前置きに変更した。現行の図を作るのは suppfigure_merged_*_261103.py。
#!/usr/bin/env python3
"""Supplementary Figure 4 — chromosomal occupancy of both methylation systems
against their own motif-conditioned nulls.

Rebuild of the 240.1 mm original at 174 mm column width (2026-10-14).  Site
sets are re-derived from the deposited pileups under the reproducible pooled
rule and asserted against the canonical census before anything is drawn.

    python suppfigure4_core_share.py
"""
import json
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, genome
from _suppstyle import (BLUE, RED,  # noqa: E402
                        DATA_DIR, NAR_WIDTH_IN, apply_figure_style,
                        figure_panel_letters, report, save, source_note)

N = json.load(open(DATA_DIR / "numbers.json"))
CORE = (1_500_000, 7_167_507)
ORIC_MB = 4.27
BINW = 50_000
TPC = {"T1": "#b4b4b4", "T2": "#767676", "T3": "#2b2b2b"}  # FIG23: 時点は中立グレー段
XLIM0 = -1.45          # room for the per-timepoint n labels

G = genome()
GLEN = len(G)

SETS = {
    "tgg": ("pileup_GCCGGC", "21839",
            [m.start() + 4 for m in re.finditer("(?=TGGCCGGC)", G)]
            + [m.start() + 3 for m in re.finditer("(?=GCCGGCCA)", G)]),
    "c4": ("pileup_AAGCCCG_4mC", "21839",
           [m.start() + 4 for m in re.finditer("(?=AAGCCCG)", G)]
           + [m.start() + 2 for m in re.finditer("(?=CGGGCTT)", G)]),
    "a01": ("pileup_AAGCCCG_6mA", "a",
            [m.start() + i for m in re.finditer("(?=AAGCCCG)", G) for i in (0, 1)]
            + [m.start() + 5 + i for m in re.finditer("(?=CGGGCTT)", G) for i in (0, 1)]),
}
# the GCCGGC extract covers every GCCGGC cytosine; restrict it to TGGCCGGC
TGG_KEYS = {(m.start() + 4, "+") for m in re.finditer("(?=TGGCCGGC)", G)} | \
           {(m.start() + 3, "-") for m in re.finditer("(?=GCCGGCCA)", G)}


def core_frac(positions):
    a = np.asarray(positions)
    return 100 * ((a >= CORE[0]) & (a <= CORE[1])).mean()


def called_positions(setname, timepoint):
    subdir, _code, _univ = SETS[setname]
    acc = {}
    for lib in LIBS[timepoint]:
        d = pd.read_csv(DATA_DIR / subdir / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, strand, nv, nm, nf in d.itertuples(index=False):
            key = (pos, strand)
            if setname == "tgg" and key not in TGG_KEYS:
                continue
            cur = acc.get(key)
            if cur is None:
                acc[key] = [nv, nm, nf]
            else:
                cur[0] += nv
                cur[1] += nm
                cur[2] += nf
    return sorted(p for (p, _s), v in acc.items()
                  if (v[0] + v[2]) >= DEPTH_MIN and v[0] > 0 and v[1] / v[0] >= FRAC_MIN)


POS = {(s, tp): called_positions(s, tp) for s in SETS for tp in ("T1", "T2", "T3")}
for tp in ("T1", "T2", "T3"):
    assert len(POS[("tgg", tp)]) == N["TGG_called"][tp], (tp, len(POS[("tgg", tp)]))
    assert len(POS[("c4", tp)]) == N["AAG_C4_called"][tp], (tp, len(POS[("c4", tp)]))
    assert len(POS[("a01", tp)]) == N["AAG_6mA_union"][tp], (tp, len(POS[("a01", tp)]))

NULL = {s: core_frac(SETS[s][2]) for s in SETS}

# 2026-10-28: 帰無を「時点ごとの可測位置」に直した。
# 旧実装は全モチーフ位置（時点非依存）を帰無にしており、T1 でのみ原稿・value_registry と
# 食い違っていた（帰無 63.4496 % 対 63.5103 %、binomial p 0.8111 対 0.8483）。
# T2/T3 は全 2,093 位置が可測なので差は出ない。原稿側の値が実データと一致することは
# Writing/core_share_verification_261028.tsv で再計算して確認済み。
# 破線は 1 本のままにする: 3 時点の帰無は 63.45-63.51 % と 0.06 ポイントしか離れておらず、
# 3 本描いても判別できないため。算出法と範囲は source_note に書く。
def measurable_positions(setname, timepoint):
    """深度条件だけを課した位置集合（修飾レベルの条件は課さない）。"""
    subdir, _code, _univ = SETS[setname]
    acc = {}
    for lib in LIBS[timepoint]:
        d = pd.read_csv(DATA_DIR / subdir / f"{lib}.tsv", sep="\t", header=None,
                        names=["pos", "strand", "Nvalid", "Nmod", "Nfail"])
        for pos, strand, nv, nm, nf in d.itertuples(index=False):
            key = (pos, strand)
            if setname == "tgg" and key not in TGG_KEYS:
                continue
            cur = acc.get(key)
            if cur is None:
                acc[key] = [nv, nm, nf]
            else:
                cur[0] += nv
                cur[1] += nm
                cur[2] += nf
    return sorted(p for (p, _s), v in acc.items() if (v[0] + v[2]) >= DEPTH_MIN)


NULLT = {(s, tp): core_frac(measurable_positions(s, tp))
         for s in SETS for tp in ("T1", "T2", "T3")}

rows = []
for s, lbl in [("tgg", "GCCGGC hexamer 4mC (TGGCCGGC)"), ("c4", "AAGCCCG 4mC (C5)"),
               ("a01", "AAGCCCG 6mA (A0/A1)")]:
    for tp in ("T1", "T2", "T3"):
        ps = POS[(s, tp)]
        cp = core_frac(ps)
        k = int(round(len(ps) * cp / 100))
        rows.append((s, lbl, tp, len(ps), cp, NULLT[(s, tp)],
                     binomtest(k, len(ps), NULLT[(s, tp)] / 100).pvalue))
CF = pd.DataFrame(rows, columns=["key", "mark", "tp", "n_sites", "core_pct",
                                 "null_pct", "p_binom"])

apply_figure_style(frame="boxed", sizes=(7.0, 6.2, 6.0))


def landscape(ax, key, title, col, headroom=1.0):
    lane_tops = []
    for i, tp in enumerate(("T1", "T2", "T3")):
        ps = np.array(POS[(key, tp)])
        h, edges = np.histogram(ps, bins=np.arange(0, GLEN + BINW, BINW))
        mx = max(h.max(), 1)
        ax.bar(edges[:-1] / 1e6, h, width=BINW / 1e6 * 0.95, bottom=i * mx * 1.15,
               color=TPC[tp], edgecolor="none")
        ax.text(-0.18, i * mx * 1.15 + mx * 0.45, f"{tp}  n={len(ps):,}", fontsize=6.0,
                ha="right", va="center", color=col)
        lane_tops.append((i + 1) * mx * 1.15)
    ax.axvspan(CORE[0] / 1e6, CORE[1] / 1e6, color="#6b6b6b", alpha=0.07, zorder=0)
    ax.axvline(ORIC_MB, color="#6b6b6b", lw=0.9, ls="--", zorder=1)
    ax.set_xlabel("Chromosome position (Mb)")
    ax.set_yticks([])
    ax.set_ylim(0, max(lane_tops) * headroom)
    ax.set_xlim(XLIM0, 8.9)
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.set_title(title, loc="left", fontsize=6.6)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)


def corepanel(ax, keys, title):
    xs, hs, cs, xl = [], [], [], []
    for ki, k in enumerate(keys):
        for j, tp in enumerate(("T1", "T2", "T3")):
            row = CF[(CF.key == k) & (CF.tp == tp)].iloc[0]
            xs.append(ki * 3.8 + j)
            hs.append(row.core_pct)
            cs.append(TPC[tp])
            xl.append(tp)
    ax.bar(xs, hs, width=0.72, color=cs, edgecolor="none")
    for x, h in zip(xs, hs):
        ax.text(x, h + 3.4, f"{h:.1f}", ha="center", fontsize=6.0)
    # 3 時点の帰無は 0.06 ポイント差なので 1 本に畳む（範囲は source_note に明記）
    nulln = float(np.mean([NULLT[(keys[0], tp)] for tp in ("T1", "T2", "T3")]))
    ax.axhline(nulln, color="#444444", lw=1.0, ls="--")
    # [D119 で削除: 破線の凡例文]
    ax.set_xticks(xs)
    ax.set_xticklabels(xl)
    ax.set_ylabel("Sites in core (%)")
    ax.set_ylim(0, 100)
    ax.set_title(title, loc="left", fontsize=6.6, pad=12)
    if len(keys) == 2:
        for gx, gl in ((1.0, "4mC (C5)"), (4.8, "6mA (A1/A2)")):
            ax.annotate(gl, xy=(gx, 0), xycoords=ax.get_xaxis_transform(),
                        xytext=(0, -15), textcoords="offset points",
                        ha="center", va="top", fontsize=6.0)


def build():
    fig = plt.figure(figsize=(NAR_WIDTH_IN, 5.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.95, 1.0], height_ratios=[1.0, 2.30],
                          hspace=0.42, wspace=0.42,
                          left=0.115, right=0.985, top=0.905, bottom=0.195)

    ax_a = fig.add_subplot(gs[0, 0])
    # 2026-10-22: この見出し色は "#1b4965"（紺）だった。同じ図のパネル c では
    # 4mC の群見出しを正準の赤、6mA を青にしてあるので、4mC の量を指す見出しが
    # パネル a だけ紺だと図の中で色が一致しない。RED に合わせる。
    landscape(ax_a, "tgg",
              "TGGCCGGC 4mC sites per 50 kb — core shaded, dashed line oriC", RED,
              headroom=1.32)

    ax_b = fig.add_subplot(gs[0, 1])
    corepanel(ax_b, ["tgg"], "")

    ax_c = fig.add_subplot(gs[1, 0])
    # レーン間隔は棒の最大高から決める。固定値 12 / 48 では T3 レーンの背の高い棒が
    # 群の見出しに食い込み、見出しを上げるとパネル表題に当たった（FIG20）。
    lanes_c = [("c4", "AAGCCCG 4mC (C5)"), ("a01", "AAGCCCG 6mA (A1/A2)")]
    HMAX = max(np.histogram(np.array(POS[(k, tp)]),
                            bins=np.arange(0, GLEN + BINW, BINW))[0].max()
               for k, _ in lanes_c for tp in ("T1", "T2", "T3"))
    LANE = int(HMAX) + 2
    GROUP = 3 * LANE + int(HMAX) + 4      # 群の見出し行ぶんの余地
    for i, (key, nm) in enumerate(lanes_c):
        for j, tp in enumerate(("T1", "T2", "T3")):
            h, edges = np.histogram(np.array(POS[(key, tp)]),
                                    bins=np.arange(0, GLEN + BINW, BINW))
            base = i * GROUP + j * LANE
            ax_c.bar(edges[:-1] / 1e6, h, width=BINW / 1e6 * 0.95, bottom=base,
                     color=TPC[tp], edgecolor="none")
            ax_c.text(-0.18, base + LANE * 0.42, f"{tp}  n={len(POS[(key, tp)]):,}", fontsize=6.0,
                      ha="right", va="center", color=BLUE if i else RED)
        ax_c.text(0.05, i * GROUP + 2 * LANE + HMAX + 1, nm, fontsize=6.0,
                  ha="left", va="bottom",
                  fontweight="bold", color=BLUE if i else RED)
    ax_c.axvspan(CORE[0] / 1e6, CORE[1] / 1e6, color="#6b6b6b", alpha=0.07, zorder=0)
    ax_c.axvline(ORIC_MB, color="#6b6b6b", lw=0.9, ls="--", zorder=1)
    ax_c.set_xlabel("Chromosome position (Mb)")
    # 自動スケールはテキストを勘定しないので、最上段の群見出しが軸の外（表題の位置）に
    # 出てしまう。見出し 1 行ぶんの天井を明示する。
    ax_c.set_ylim(0, GROUP + 2 * LANE + HMAX + 1 + 1.1 * LANE)
    ax_c.set_yticks([])
    ax_c.set_xlim(XLIM0, 8.9)
    ax_c.set_xticks([0, 2, 4, 6, 8])
    ax_c.set_title("AAGCCCG sites per 50 kb, the two marks kept separate",
                   loc="left", fontsize=6.6)
    for sp in ("left", "right", "top"):
        ax_c.spines[sp].set_visible(False)

    # 下段の行は c のレーン厚を a に合わせるため 2.30 倍にしてある。d を素直に
    # この行いっぱいに置くと棒だけが縦に伸びるので、行の上側 1.0 ぶんだけを使い
    # b と同じ高さに揃える。
    ax_d = fig.add_subplot(gs[1, 1].subgridspec(2, 1, height_ratios=[1.0, 1.30],
                                               hspace=0.0)[0])
    corepanel(ax_d, ["c4", "a01"], "")

    figure_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d], "abcd",
                         dx=[0.105, 0.085, 0.105, 0.085], dy=0.012)
    source_note(fig,
                "Site rule: three libraries pooled per timepoint, pooled true depth (Nvalid_cov + Nfail) >= 30, pooled modified fraction\n"
                ">= 50 % (deposited pileups). The dashed null is the core share of the positions each mark can occupy that\n"
                "clear the depth condition at that timepoint; one line is drawn because the three timepoints differ by\n"
                f"less than 0.1 percentage point (TGGCCGGC {min(NULLT[('tgg', t)] for t in ('T1','T2','T3')):.2f}-"
                f"{max(NULLT[('tgg', t)] for t in ('T1','T2','T3')):.2f} %, "
                f"AAGCCCG {min(NULLT[('c4', t)] for t in ('T1','T2','T3')):.2f}-"
                f"{max(NULLT[('c4', t)] for t in ('T1','T2','T3')):.2f} %).",
                y=0.008)
    return fig


if __name__ == "__main__":
    # 出荷名は SUPERSEDED_SuppFigure2_pre261104.png（補足図 4 -> 2 の改番、2026-10-24 に出力名を追随）。
    report(save(build(), "SUPERSEDED_SuppFigure2_pre261104.png"))
