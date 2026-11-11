# 2026-11-04 D175: この生成器は補足図の統合で役目を終えた。出荷名に書き込むと別の図を潰すため、
# 出力名を SUPERSEDED_ 前置きに変更した。現行の図を作るのは suppfigure_merged_*_261103.py。
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import re
import pathlib
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Global state
mpl.rcParams["savefig.bbox"] = "standard"

# ── helper functions (from figure-style skill) ─────────────────────────────
MM = 1 / 25.4
META_GREY = "#666666"

def apply_figure_style(sizes=(8, 7, 6)):
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.size": sizes[1],
        "axes.titlesize": sizes[0],
        "axes.labelsize": sizes[1],
        "xtick.labelsize": sizes[2],
        "ytick.labelsize": sizes[2],
        # 2026-10-28 D121: 補足図も四辺の枠に揃える。この生成器は _suppstyle を
        # 使わず自前の style を持っていたので、ここで個別に boxed にする。
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "lines.linewidth": 1.2,
        "legend.frameon": False,
        "savefig.dpi": 300,
        "savefig.bbox": "standard",
        "figure.dpi": 150,
    })

def panel_letter(ax, letter, x=None, y=None):
    # 2026-10-21: 軸座標だと目盛ラベルや軸ラベルと衝突するので図座標の固定位置に置く
    fig_ = ax.get_figure()
    bb = ax.get_position()
    fig_.text(max(0.004, bb.x0 - 0.058), 0.985, letter, fontsize=9, fontweight="bold", va="top", ha="left")

# ── data paths ─────────────────────────────────────────────────────────────
RNA = pathlib.Path("/Users/okaban/bioinfo/rna-seq")
OUT = RNA / "15_paper_figures/scripts/recompute_261015"
S = OUT / "pileup_cache_261015"
GEN = RNA / "11_epigenome_integration/data/NC_003888.3.fna"

seq = "".join(l.strip() for l in GEN.read_text().split("\n")[1:]).upper()
COMP = str.maketrans("ACGT", "TGCA")
LIB = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}

# ── pileup pool ─────────────────────────────────────────────────────────────
POOL = {}

def pooled(tp, code):
    k = (tp, code)
    if k in POOL:
        return POOL[k]
    acc = {}
    for lib in LIB[tp]:
        d = np.load(S / f"{lib}.npz")
        m = (d["code"] == code) & (d["refok"] == 1)
        st, pl, nv, nm, nf = d["start"][m], d["plus"][m], d["Nvalid"][m], d["Nmod"][m], d["Nfail"][m]
        key = st.astype(np.int64) * 2 + pl.astype(np.int64)
        for k2, v, mo, fa in zip(key, nv, nm, nf):
            a = acc.get(k2)
            if a is None:
                acc[k2] = [v, mo, fa]
            else:
                a[0] += v; a[1] += mo; a[2] += fa
    POOL[k] = acc
    return acc

def pooled_c(tp, code):
    return pooled(tp, code)

# ── motif utilities ─────────────────────────────────────────────────────────
def instances(motif, mod_offset):
    out = []
    for m in re.finditer(f"(?={motif})", seq):
        out.append((m.start() + mod_offset, 1, m.start()))
    rcm = motif.translate(COMP)[::-1]
    L = len(motif)
    for m in re.finditer(f"(?={rcm})", seq):
        out.append((m.start() + (L - 1 - mod_offset), 0, m.start()))
    return out

def flank_base(start, plus, L, off):
    if plus:
        i = start + (off if off < 0 else L - 1 + off)
        return seq[i] if 0 <= i < len(seq) else None
    i = start + (L - 1 - off if off < 0 else -off)
    if not (0 <= i < len(seq)):
        return None
    return seq[i].translate(COMP)

def status_fn(acc, pos, plus):
    a = acc.get(pos * 2 + plus)
    if a is None or (a[0] + a[2]) < 30:
        return None
    return 1 if a[0] and a[1] / a[0] >= 0.5 else 0

# ── scan function ───────────────────────────────────────────────────────────
is_H = lambda st, plus, L: flank_base(st, plus, L, -1) in ("A", "C", "T")

def scan_cond(motif, mod_offset, code, label, cond=None, positions=range(-8, 9)):
    L = len(motif)
    rows = []
    for tp in ("T1", "T2", "T3"):
        acc = pooled_c(tp, code)
        dat = []
        for pos, plus, st in instances(motif, mod_offset):
            if cond and not cond(st, plus, L):
                continue
            a = acc.get(pos * 2 + plus)
            if a is None or (a[0] + a[2]) < 30:
                continue
            s = 1 if a[0] and a[1] / a[0] >= 0.5 else 0
            dat.append((st, plus, s))
        for off in positions:
            if off == 0:
                continue
            tab = {}
            for st2, plus2, s2 in dat:
                b = flank_base(st2, plus2, L, off)
                if b in "ACGT":
                    tab.setdefault(b, [0, 0])[s2] += 1
            nm = sum(v[1] for v in tab.values())
            nu = sum(v[0] for v in tab.values())
            for b, (u, mo) in sorted(tab.items()):
                if min(mo + u, nm + nu - mo - u) < 5:
                    continue
                orr, p = fisher_exact([[mo, nm - mo], [u, nu - u]])
                rows.append(dict(motif=label, timepoint=tp, position=off, base=b,
                                 n_mod_with=mo, n_mod_total=nm, n_unmod_with=u, n_unmod_total=nu,
                                 odds_ratio=round(float(orr), 4) if np.isfinite(orr) else np.nan, p=p))
    df = pd.DataFrame(rows)
    for tp2 in df.timepoint.unique():
        m2 = df.timepoint == tp2
        df.loc[m2, "q_BH"] = multipletests(df.loc[m2, "p"], method="fdr_bh")[1]
        df.loc[m2, "n_tests"] = int(m2.sum())
    sig = df[df.q_BH < 0.05].groupby(["position", "base"]).timepoint.nunique()
    df["survives_3tp"] = df.apply(lambda r: sig.get((r.position, r.base), 0) == 3, axis=1)
    return df

# ── run scans ───────────────────────────────────────────────────────────────
d6H = scan_cond("AAGCCCG", 0, 0, "HAAGCCCG 6mA(A1)", cond=is_H)
d4A = scan_cond("AAGCCCG", 4, 1, "AAGCCCG 4mC(pos5)", cond=is_H)

# ── TGGCCGGC detection control ──────────────────────────────────────────────
from scipy.stats import mannwhitneyu

ctl_rows = []
for tp in ("T1", "T2", "T3"):
    acc = pooled_c(tp, 1)
    rec = []
    for pos, plus, st in instances("TGGCCGGC", 4):
        a = acc.get(pos * 2 + plus)
        if a is None or (a[0] + a[2]) < 30:
            continue
        s = 1 if a[0] and a[1] / a[0] >= 0.5 else 0
        up = "".join(flank_base(st, plus, 8, o) or "N" for o in range(-10, 0))
        dn = "".join(flank_base(st, plus, 8, o) or "N" for o in range(1, 11))
        f8 = "".join(flank_base(st, plus, 8, o) or "N" for o in list(range(-8, 0)) + list(range(1, 9)))
        rec.append(dict(status=s,
                        gc8=(f8.count("G") + f8.count("C")) / 16,
                        grun=int(bool(re.search(r"G{3,}", up + dn))),
                        failfrac=a[2] / (a[0] + a[2])))
    d = pd.DataFrame(rec)
    mo = d[d.status == 1]; un = d[d.status == 0]
    ctl_rows.append(dict(timepoint=tp, n_mod=len(mo), n_unmod=len(un),
                         gc8_mod=mo.gc8.mean(), gc8_unmod=un.gc8.mean(),
                         grun_mod=mo.grun.mean(), grun_unmod=un.grun.mean(),
                         fail_mod=mo.failfrac.mean(), fail_unmod=un.failfrac.mean()))
ctl = pd.DataFrame(ctl_rows)

# ── occupancy ladder data ───────────────────────────────────────────────────
# 2026-10-22: 以前は host.artifact_path("4f71131b-...") で artifact を直読みしており、
# Claude Science のカーネル外では NameError になって図を再生成できなかった。同じ
# artifact（SuppTable27_occupancy_ladder_261019.tsv、sha256 de1cac66758f1239...）を
# data/ に取り込み、素のパスから読む。内容は同一。
LADDER = pathlib.Path(__file__).resolve().parent / "data" / "occupancy_ladder_261020.tsv"
lad = pd.read_csv(LADDER, sep="\t")
assert len(lad) == 15 and set(lad.timepoint) == {"T1", "T2", "T3"}, lad.shape

lad2 = lad.copy()
lad2["tier"] = lad2.definition.map({
    "AAGCCCG (all)": 1, "HAAGCCCG": 2, "+ -3=G": 3, "+ -5=C & -3=G": 4,
    "(ref) same-instance 4mC C5": "ref"
})
lad2 = lad2[["timepoint", "tier", "definition", "n", "called", "occupancy_pct"]]

# ── plot ────────────────────────────────────────────────────────────────────
apply_figure_style(sizes=(8, 7, 6))
mpl.rcParams["savefig.bbox"] = "standard"

# 2026-10-22: was #1b6ca8 / #b5651d, a fourth 4mC/6mA pair. Bound to canon.
import importlib, sys as _sys
_sys.path.insert(0, "/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
_U = importlib.import_module("00_shared_utils")
C6, C4 = _U.COL_6mA, _U.COL_4mC
tiers = [1, 2, 3, 4]
ref = lad2[lad2.tier == "ref"].occupancy_pct

fig, axes = plt.subplots(1, 3, figsize=(172.7 * MM, 64 * MM),
                         gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]})

# (a) occupancy ladder
ax = axes[0]
for tp, mk in zip(("T1", "T2", "T3"), ("o", "s", "^")):
    g = lad2[(lad2.timepoint == tp) & (lad2.tier != "ref")].sort_values("tier")
    ax.plot(tiers, g.occupancy_pct, marker=mk, color=C6, lw=1.2, ms=3.6, label=tp)
ax.axhspan(ref.min(), ref.max(), color=C4, alpha=0.18, lw=0)
ax.text(1.05, ref.mean() + 3.0, "4mC, same instances", color=C4, fontsize=6.0, va="bottom", ha="left")
ax.set_xticks(tiers)
ax.set_xticklabels(
    [f"{l}\nn={n}" for l, n in zip(["AAGCCCG", "HAAGCCCG", "+ −3 G", "+ −5 C"],
     lad2[(lad2.timepoint == "T1") & (lad2.tier != "ref")].sort_values("tier").n)],
    fontsize=6.0)
ax.set_ylim(0, 100)
ax.set_xlim(0.7, 4.4)
ax.set_ylabel("6mA occupancy (%)")

ax.legend(frameon=False, fontsize=6, loc="lower right", handlelength=1.2)

# (b) positional scan
ax = axes[1]
for d, col, nm, mk in ((d6H, C6, "6mA (position 1)", "o"), (d4A, C4, "4mC (position 5)", "s")):
    t1 = d[d.timepoint == "T1"]
    best = t1.loc[t1.groupby("position").q_BH.idxmin()].sort_values("position")
    ax.plot(best.position, best.odds_ratio, marker=mk, ms=3.0, lw=1.0, color=col, label=nm, alpha=0.85)
    sp = set(d[d.survives_3tp].position)
    for r in best[best.position.isin(sp)].itertuples():
        ax.plot(r.position, r.odds_ratio, marker=mk, ms=7, mfc="none", mec=col, mew=1.4, lw=0)
        ax.annotate(r.base, (r.position, r.odds_ratio), textcoords="offset points", xytext=(0, 8),
                    ha="center", fontsize=6.5, color=col, weight="bold")
ax.axhline(1, color=META_GREY, lw=0.8, ls=":")
ax.set_yscale("log")
ax.set_ylim(0.16, 7.5)
ax.set_yticks([0.2, 1.0, 5.0]); ax.set_yticklabels(["0.2", "1", "5"])
ax.minorticks_off()
ax.set_xlabel("Position relative to the motif")
ax.set_ylabel("Odds ratio, most significant base")

ax.legend(frameon=False, fontsize=6.0, loc="lower left", handlelength=1.1,
          bbox_to_anchor=(0.0, 0.0), ncol=1, labelspacing=0.3, borderpad=0.2)
# [D119 で削除: ring: survives all three timepoints]

# (c) detection control
ax = axes[2]
met = [("±8 bp GC", "gc8_mod", "gc8_unmod"),
       ("G-run within\n±10 bp", "grun_mod", "grun_unmod"),
       ("call-failure\nfraction", "fail_mod", "fail_unmod")]
x = np.arange(3)
w = 0.36
mv = [ctl[c[1]].mean() for c in met]
uv = [ctl[c[2]].mean() for c in met]
ax.bar(x - w / 2, mv, w, color=C4, label="TGGCCGGC methylated, 1,886")
ax.bar(x + w / 2, uv, w, color=META_GREY, label="TGGCCGGC unmethylated, 205")
# 2026-10-28: 棒の上に焼いていた値ラベル 6 個を削除した。
# 衝突検査が `text overlaps text: [('0.72', '0.77')]` と
# `text over data marks: ['0.33', '0.66', '0.72']` を出しており、
# 著者の報告（「グラフ内にテキストとグラフが重なっている」）と一致した。
# 著者の標準的な指示は「図内のテキストは最小限にし、説明的なものは凡例に回す」なので、
# 値は凡例本文に移した。棒の高さは y 軸の目盛りから読める。
# 値を図に戻す必要が出た場合は、ax.bar_label(..., padding=) を使い
# _suppstyle.save(..., allow_over_data=(...)) で明示宣言すること（黙認は不可）。
_PANEL_C_VALUES = {met[i][0].replace(chr(10), " "): (round(mv[i], 2), round(uv[i], 2))
                   for i in range(3)}
print("panel c の値（凡例に記載する）:", _PANEL_C_VALUES)
ax.set_xticks(x)
ax.set_xticklabels(["±8 bp GC", "G-run,\n±10 bp", "call-failure\nfraction"], fontsize=6.0)
ax.set_ylim(0, 1.30)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
ax.set_ylabel("Fraction of TGGCCGGC positions")  # 2026-10-28 D131: 分母は 1,886+205 の
# 片鎖座標なので、規約 (position=片鎖上の 1 座標 / site=モチーフ 1 個が占める領域) では positions。

ax.legend(frameon=False, fontsize=6, loc="upper center", ncol=1, handlelength=1.2,
          bbox_to_anchor=(0.5, 1.02))

for a_, L in zip(axes, "abc"):
    panel_letter(a_, L)

fig.subplots_adjust(left=0.105, right=0.985, top=0.90, bottom=0.30, wspace=0.56)
# 出荷名は SUPERSEDED_SuppFigure13_pre261104.png。2026-10-24 まではスクリプトの実行ディレクトリに
# SuppFigure14_determinants_261020.png を書いており、出荷物との対応が
# ファイル名からは追えなかった。他の補足図と同じく SUPPFIG_OUT（既定は
# Writing/fig_images）に出荷名で書く。
_OUT_DIR = pathlib.Path(os.environ.get(
    "SUPPFIG_OUT", "/Users/okaban/obsidian/Research/rna-seq/Writing/fig_images"))
_OUT_DIR.mkdir(parents=True, exist_ok=True)
_out = _OUT_DIR / "SUPERSEDED_SuppFigure13_pre261104.png"

# 2026-10-28: 図内テキストの衝突検査を通してから保存する。
# それまでこのスクリプトは fig.savefig で直接保存しており、figure_conventions.tsv の
# FIG20（図内テキストが他のテキスト・データマーク・キャンバス外と重なってはならない）と
# FIG21（保存経路は assert_no_text_collisions を呼ばなければならない）を
# 迂回していた。補足図の生成器 12 本のうち検査を迂回していたのはこの 1 本で、
# 著者が「グラフ内にテキストとグラフが重なっている。この指摘毎回している」と
# 報告した図がまさにこれだった。他の 11 本は _suppstyle.save() を経由しており
# 同じ検査を通っている。
def _assert_no_text_collisions(figure, name):
    import importlib.util
    src = pathlib.Path(__file__).resolve().parent.parent / "00_shared_utils.py"
    spec = importlib.util.spec_from_file_location("_shared_utils_for_sf14", src)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.assert_no_text_collisions(figure, name=name)


_assert_no_text_collisions(fig, "SUPERSEDED_SuppFigure13_pre261104.png")
fig.savefig(_out, dpi=300)
print("saved:", _out)