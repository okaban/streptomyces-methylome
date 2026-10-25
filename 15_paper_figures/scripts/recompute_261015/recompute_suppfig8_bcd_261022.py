#!/usr/bin/env python3
"""補足図 8（旧 SuppFigure14）パネル b/c/d を正準規則で再計算する。

現行の suppfigure14_bgc_relationships.py は MIN_VALID = 10（per-library Nvalid >= 10）を
測定可能の下限にしている。この規則は +15.5 pp の上方バイアスを持つとして撤回済み
（retracted: 深度フィルタ版）。本スクリプトは正準規則だけを使う:

    1. 時点ごとに 3 ライブラリをプール（Nvalid / Nmod / Nfail を (start, strand, code) で合算）
    2. 測定可能 = プール真深度 (Nvalid + Nfail) >= 30
    3. 参照塩基一致 refok（6mA: + 鎖 A / - 鎖 T、4mC: + 鎖 C / - 鎖 G）
    4. メチル化部位 = プール Nmod / Nvalid >= 0.5

パネル b/c が描くのは「水準(level)」であり、4 の >= 0.5 を level の算出前にかけると
低水準側を打ち切って上方バイアスが戻る。したがって level は 1-3 を満たす位置の
プール Nmod/Nvalid とし、>= 0.5 は「部位数」を数えるときだけ使う（パネル d の
unassigned 6mA クラスの定義に必要）。

出力: suppfig8_panelbcd_canonical_261022.tsv（分子/分母を列に残す）

実行: python3 recompute_261015/recompute_suppfig8_bcd_261022.py
"""
import pathlib
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent                                    # 15_paper_figures/scripts
sys.path.insert(0, str(ROOT / "supp_261014"))
from _gccggc import genome                            # noqa: E402

CACHE = HERE / "pileup_cache_261015"
DATA = ROOT / "supp_261014" / "data"
TSS_CSV = ("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
           "18_tss_analyses/jeong2016_all_tss.csv")
FEATURECOUNTS = ("/Users/okaban/bioinfo/rna-seq/03_quantification/data/featureCounts/"
                 "featureCounts_M145_s2.txt")
OUT = HERE / "suppfig8_panelbcd_canonical_261022.tsv"

TP = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
      "T3": ["3-2", "3-3", "3-4"]}
DEPTH_MIN = 30            # プール真深度 Nvalid + Nfail
FRAC_MIN = 0.5            # プール Nmod / Nvalid（部位と数えるとき）
CODE_6mA, CODE_4mC = 0, 1
FLANK_BP = 50_000         # パネル d の flank 幅（領域境界の外側、他領域と重なる分は除く）

G = genome()
GLEN = len(G)


# ---------------------------------------------------------------- motif sets
def motif_keys():
    """(pos, plus) の集合。house の定義（depth_rule_audit_261015.py）と同一。"""
    tgg_C, aag_C, aag_A = set(), set(), set()
    for m in re.finditer("(?=TGGCCGGC)", G):
        tgg_C.add((m.start() + 4, True))
    for m in re.finditer("(?=GCCGGCCA)", G):
        tgg_C.add((m.start() + 3, False))
    for m in re.finditer("(?=AAGCCCG)", G):
        i = m.start()
        aag_C.add((i + 4, True))
        aag_A |= {(i, True), (i + 1, True)}
    for m in re.finditer("(?=CGGGCTT)", G):
        j = m.start()
        aag_C.add((j + 2, False))
        aag_A |= {(j + 6, False), (j + 5, False)}
    # 六塩基コア GCCGGC の修飾シトシン（両鎖）。配列 availability の分母としての
    # 呼称であり、4mC 標的の 8 塩基 TGGCCGGC とは別物。
    hexa_C = set()
    for m in re.finditer("(?=GCCGGC)", G):
        hexa_C |= {(m.start() + 2, True), (m.start() + 3, False)}
    return tgg_C, aag_C, aag_A, hexa_C


TGG_C, AAG_C, AAG_A, HEXA_C = motif_keys()


# ---------------------------------------------------------------- pooling
def pool(tp):
    """時点 tp の 3 ライブラリをプール -> DataFrame(pos, plus, code, Nvalid, Nmod, Nfail, refok)."""
    key_acc = None
    for lib in TP[tp]:
        z = np.load(CACHE / f"{lib}.npz")
        code = z["code"].astype(np.int64)
        keep = code <= 2
        k = (z["start"].astype(np.int64)[keep] * 6
             + code[keep] * 2 + z["plus"][keep].astype(np.int64))
        v = np.stack([z["Nvalid"][keep].astype(np.int64),
                      z["Nmod"][keep].astype(np.int64),
                      z["Nfail"][keep].astype(np.int64)], 1)
        del z
        if key_acc is None:
            key_acc, val_acc = k, v
        else:
            key_acc = np.concatenate([key_acc, k])
            val_acc = np.concatenate([val_acc, v])
    uk, inv = np.unique(key_acc, return_inverse=True)
    agg = np.zeros((len(uk), 3), np.int64)
    for j in range(3):
        agg[:, j] = np.bincount(inv, weights=val_acc[:, j], minlength=len(uk)).astype(np.int64)
    del key_acc, val_acc, inv
    pos, rem = uk // 6, uk % 6
    d = pd.DataFrame(dict(pos=pos, code=(rem // 2).astype(np.int8),
                          plus=(rem % 2).astype(bool),
                          Nvalid=agg[:, 0], Nmod=agg[:, 1], Nfail=agg[:, 2]))
    base = np.frombuffer(G.encode(), np.uint8)[d.pos.to_numpy()]
    exp_p = np.where(d.code == CODE_6mA, ord("A"), ord("C"))
    exp_m = np.where(d.code == CODE_6mA, ord("T"), ord("G"))
    d["refok"] = np.where(d.plus, base == exp_p, base == exp_m)
    d["measurable"] = (d.Nvalid + d.Nfail) >= DEPTH_MIN
    d["level"] = np.where(d.Nvalid > 0, d.Nmod / d.Nvalid.replace(0, np.nan), np.nan)
    return d


# ---------------------------------------------------------------- annotation
def promoter_intervals():
    t = pd.read_csv(TSS_CSV)
    ps = t[t.category.isin(["P", "S"]) & t.gene_id.notna()]
    iv = []
    for r in ps.itertuples():
        p = int(r.tss_position)
        iv.append((p - 300, p - 1) if r.strand in ("+", 1, "1") else (p + 1, p + 300))
    return np.array(sorted(set(iv)))


def in_intervals(pos, iv):
    k = np.searchsorted(iv[:, 0], pos, side="right") - 1
    ok = k >= 0
    out = np.zeros(len(pos), bool)
    out[ok] = pos[ok] <= iv[k[ok], 1]
    return out


def main():
    mb = pd.read_csv(DATA / "MB02_BGC_crosswalk_supplementary.csv")
    regions = mb[["region", "start", "end"]].to_numpy()
    prom_iv = promoter_intervals()

    # flank: 領域境界の外側 FLANK_BP。どの領域の内部とも重ならない部分だけを使う。
    interior = np.zeros(GLEN, bool)
    for _, s, e in regions:
        interior[int(s):int(e) + 1] = True
    flank = np.zeros(GLEN, bool)
    for _, s, e in regions:
        s, e = int(s), int(e)
        flank[max(0, s - FLANK_BP):s] = True
        flank[e + 1:min(GLEN, e + 1 + FLANK_BP)] = True
    flank &= ~interior

    CLASSES = {
        "TGGCCGGC 4mC": (TGG_C, CODE_4mC),
        "AAGCCCG C5 4mC": (AAG_C, CODE_4mC),
        "AAGCCCG 6mA": (AAG_A, CODE_6mA),
        "GCCGGC hexamer 4mC": (HEXA_C, CODE_4mC),
    }

    rows_b, rows_c, rows_d = [], [], []
    for tp in ("T1", "T2", "T3"):
        d = pool(tp)
        print(f"[{tp}] pooled rows {len(d):,}", flush=True)
        eligible = d.measurable & d.refok & (d.Nvalid > 0)

        keyset = set(zip(d.pos.to_numpy(), d.plus.to_numpy()))
        for name, (keys, code) in CLASSES.items():
            kk = np.array(sorted(keys & keyset), dtype=object)
            if not len(kk):
                continue
            kpos = np.array([k[0] for k in kk], np.int64)
            kplus = np.array([k[1] for k in kk], bool)
            sel = d.set_index(["pos", "plus", "code"]).reindex(
                pd.MultiIndex.from_arrays([kpos, kplus, np.full(len(kpos), code)],
                                          names=["pos", "plus", "code"]))
            sel = sel[sel.measurable.fillna(False) & sel.refok.fillna(False)
                      & (sel.Nvalid.fillna(0) > 0)].reset_index()

            # ---- panel d: interior vs flank
            ii = interior[sel.pos.to_numpy()]
            ff = flank[sel.pos.to_numpy()]
            for lab, m in (("interior", ii), ("flank", ff)):
                sub = sel[m]
                rows_d.append(dict(timepoint=tp, mark_class=name, locus=lab,
                                   n_positions=len(sub),
                                   Nmod=int(sub.Nmod.sum()), Nvalid=int(sub.Nvalid.sum()),
                                   level_pooled=round(sub.Nmod.sum() / sub.Nvalid.sum(), 5)
                                   if sub.Nvalid.sum() else np.nan,
                                   level_median=round(float(sub.level.median()), 5)
                                   if len(sub) else np.nan,
                                   n_called=int((sub.level >= FRAC_MIN).sum())))

            # ---- panel b: TGGCCGGC 4mC inside BGC, promoter vs rest
            if name == "TGGCCGGC 4mC":
                inb = ii
                sb = sel[inb]
                isprom = in_intervals(sb.pos.to_numpy(), prom_iv)
                for lab, m in (("promoter (-300..-1)", isprom), ("rest of cluster", ~isprom)):
                    s2 = sb[m]
                    rows_b.append(dict(timepoint=tp, locus=lab, n_positions=len(s2),
                                       Nmod=int(s2.Nmod.sum()), Nvalid=int(s2.Nvalid.sum()),
                                       level_pooled_pct=round(100 * s2.Nmod.sum() / s2.Nvalid.sum(), 2)
                                       if s2.Nvalid.sum() else np.nan,
                                       n_called=int((s2.level >= FRAC_MIN).sum())))

                # ---- panel c: per-region level vs core expression
                for reg, s_, e_ in regions:
                    m = (sb.pos >= int(s_)) & (sb.pos <= int(e_))
                    s3 = sb[m]
                    rows_c.append(dict(timepoint=tp, region=int(reg), n_positions=len(s3),
                                       Nmod=int(s3.Nmod.sum()), Nvalid=int(s3.Nvalid.sum()),
                                       level_pooled_pct=round(100 * s3.Nmod.sum() / s3.Nvalid.sum(), 2)
                                       if s3.Nvalid.sum() else np.nan))
        del d

    B = pd.DataFrame(rows_b).assign(panel="b")
    C = pd.DataFrame(rows_c).assign(panel="c")
    D = pd.DataFrame(rows_d).assign(panel="d")
    allr = pd.concat([B, C, D], ignore_index=True)
    allr.to_csv(OUT, sep="\t", index=False)
    print("\n--- panel b ---"); print(B.to_string(index=False))
    print("\n--- panel d ---"); print(D.to_string(index=False))
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
