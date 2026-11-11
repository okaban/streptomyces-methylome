#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""補足表 20a（再生成・all-read 系統）: antiSMASH 29 領域のメチル化水準と発現。

    python3 build_supptable20a_allread_261101.py

なぜ新規に書いたか
------------------
既存の `SuppTable20a_cluster_methylation_expression_261021.tsv` には生成器が
リポジトリに存在しない。2026-11-01 に本スクリプトの作成時点で
`grep -rl 'SuppTable20a' . --include=*.py` を rna-seq ツリー直下で実行した結果は
参照側の `15_paper_figures/scripts/supp_261014/build_supptable20c_261021.py`
1 件のみ（本スクリプト自身と newfig_261101/ を除く）。同じ不在は
`c04_newfig_data_inventory_261101.tsv`（NF-R5 行）でも先に記録されている。
表のレベル列が信頼度フィルタ版
パイルアップ由来か all-read 由来かを、出力ヘッダからもコードからも確定できず、
`value_registry.csv` の BGC.4mC.TPM.RHO / BGC.6mA.TPM.RHO がその表を出所に
している。生成器の無い表を引き継がない方針（generated-records）に従い、
系統を宣言したうえで書き直す。

系統の宣言
----------
レベルは **all-read**（`pileup_nofilt_261015`、modkit 0.6.4 --no-filtering、
Nfail == 0）。正典パイルアップ `pileup_canonical_261015` は modkit の既定の
信頼度フィルタを通しており、canon の F7 により LEVEL を +15.5 pp、OCCUPANCY を
+5.1 pp 押し上げる。クラスタ内外のレベル比較はレベル尺度の量なので、
フィルタ版を使ってはならない。

`pileup_nofilt_261015` は `analysis/<lib>_mapped.bam` から作られており、
3-2 と 3-4 については 2026-02 の追加フローセルを含む（実測: GCCGGC 位置の
Nvalid 中央値 3-2 = 25、3-4 = 36。正典 README の追加後の値 25x / 36x と一致し、
追加前の 15x / 31x ではない）。`methyl/260102_M145/analysis/pileup/` を読む
`build_supptable20c_261021.py` の取りこぼしはここでは起きない。

修飾位置の宣言（canon）
-----------------------
    4mC = TGGCCGGC の C5（鎖指定、認識配列そのもの。GCCGGC 六塩基ではない）
    6mA = AAGCCCG の A1 と A2
AAGCCCG の C5（4mC）は別系列なので 4mC 列には混ぜない（混ぜると分母が
2 種類の認識配列にまたがる）。

測定可能性: 時点ごとに 3 ライブラリをプールし、プール後 Nvalid >= 30。
all-read なので Nvalid が真深度である（Nfail == 0 をアサートする）。
レベル = プール Nmod / プール Nvalid（位置ごと）、領域値はその単純平均。
メチル化閾値は一切かけない（レベル尺度の量を、測っている当の量で選別しない）。

出力: 15_paper_figures/tables/supplementary/
        SuppTable20a_cluster_methylation_expression_allread_261101.tsv
"""
import os
import re

import numpy as np
import pandas as pd

RNA = "/Users/okaban/bioinfo/rna-seq"
NOFILT = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"
REGIONS = f"{RNA}/15_paper_figures/scripts/supp_261014/data/bgc_regions_29.csv"
FEATURECOUNTS = (f"{RNA}/03_quantification/data/featureCounts/"
                 "featureCounts_M145_s2.txt")
OUTDIR = f"{RNA}/15_paper_figures/tables/supplementary"
OUT = os.path.join(OUTDIR,
                   "SuppTable20a_cluster_methylation_expression_allread_261101.tsv")

LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
        "T3": ["3-2", "3-3", "3-4"]}
BEDCOL = ["chrom", "start", "end", "name", "score", "strand", "ts", "te",
          "color", "Nvalid", "pct", "Nmod", "Ncanon", "Nother", "Ndel",
          "Nfail", "Ndiff", "Nnocall"]
MIN_DEPTH = 30
MIN_POS = 5          # a region needs this many measurable positions to be used


def genome_seq():
    return "".join(l.strip() for l in open(GENOME)
                   if not l.startswith(">")).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def tggccggc_positions(seq):
    """0-based strand-oriented positions of the modified C of TGGCCGGC."""
    L = len(seq)
    out = set(m.start() + 4 for m in re.finditer("(?=TGGCCGGC)", seq))
    out |= set(L - 1 - (m.start() + 4) for m in re.finditer("(?=TGGCCGGC)",
                                                            rc(seq)))
    return out


def pooled_levels(names, keep=None):
    """{timepoint: DataFrame(start, strand, Nvalid, Nmod, level)} on all reads."""
    per = {}
    for lib in sorted({l for v in LIBS.values() for l in v}):
        d = pd.read_csv(f"{NOFILT}/{lib}_nofilt.bed", sep="\t", header=None,
                        names=BEDCOL)
        assert int(d["Nfail"].sum()) == 0, f"{lib}: not an all-read pileup"
        d = d[d["name"].isin(names)]
        if keep is not None:
            d = d[d["start"].isin(keep)]
        per[lib] = (d.drop_duplicates(subset=["start", "strand"])
                    .set_index(["start", "strand"])[["Nvalid", "Nmod"]])
    out = {}
    for tp, libs in LIBS.items():
        idx = per[libs[0]].index
        for l in libs[1:]:
            idx = idx.union(per[l].index)
        s = sum(per[l].reindex(idx).fillna(0) for l in libs)
        s = s[s["Nvalid"] >= MIN_DEPTH].copy()
        s["level"] = s["Nmod"] / s["Nvalid"]
        out[tp] = s.reset_index()
    return out


def tpm_table():
    d = pd.read_csv(FEATURECOUNTS, sep="\t", comment="#")
    cnt = d.set_index("Geneid").iloc[:, 5:]
    cnt.columns = [c.split("/")[-1].split(".")[0] for c in cnt.columns]
    rpk = cnt.div(d.set_index("Geneid")["Length"], axis=0)
    tpm = rpk / rpk.sum() * 1e6
    out = pd.DataFrame(index=tpm.index)
    for tp, libs in LIBS.items():
        cols = [f"M145_{l.replace('-', '_')}" for l in libs]
        out[f"TPM_{tp}"] = tpm[cols].mean(axis=1)
    # featureCounts writes ';'-joined coordinates for multi-interval features
    meta = pd.DataFrame({
        "Geneid": d["Geneid"].values,
        "Start": [min(int(x) for x in str(v).split(";")) for v in d["Start"]],
        "End": [max(int(x) for x in str(v).split(";")) for v in d["End"]],
    }).set_index("Geneid")
    out = out.join(meta)
    assert out["Start"].notna().all(), "gene coordinates failed to join"
    return out


def main():
    seq = genome_seq()
    reg = pd.read_csv(REGIONS)
    assert len(reg) == 29, len(reg)

    lv4 = pooled_levels(["21839,GCCGGC,2"], keep=tggccggc_positions(seq))
    lv6 = pooled_levels(["a,AAGCCCG,0", "a,AAGCCCG,1"])
    tpm = tpm_table()

    rows = []
    for _, r in reg.iterrows():
        s0, e0 = int(r["start"]) - 1, int(r["end"]) - 1      # to 0-based
        row = dict(region=f"NC_003888.3_region{int(r.region_number)}",
                   known_cluster=r["known_cluster"], status=r["status"],
                   region_class=r["region_class"],
                   start=int(r["start"]), end=int(r["end"]),
                   length_bp=int(r["end"]) - int(r["start"]) + 1,
                   n_genes=int(r["n_genes"]))
        for mark, lv in (("4mC", lv4), ("6mA", lv6)):
            npos = None
            for tp in LIBS:
                d = lv[tp]
                inside = d[(d["start"] >= s0) & (d["start"] <= e0)]
                row[f"{mark}_level_{tp}"] = (round(float(inside["level"].mean()), 4)
                                             if len(inside) else np.nan)
                row[f"n_{mark}_positions_{tp}"] = int(len(inside))
                npos = len(inside) if npos is None else min(npos, len(inside))
            row[f"included_{mark}"] = "yes" if npos >= MIN_POS else "no"
        mid = (tpm["Start"] + tpm["End"]) / 2
        g = tpm[(mid >= r["start"]) & (mid <= r["end"])]
        for tp in LIBS:
            row[f"TPM_{tp}"] = round(float(g[f"TPM_{tp}"].mean()), 2) if len(g) else np.nan
        row["n_genes_tpm"] = int(len(g))
        row["coreTPM_T1"] = r["coreTPM_T1"]
        row["coreTPM_T2"] = r["coreTPM_T2"]
        row["coreTPM_T3"] = r["coreTPM_T3"]
        rows.append(row)

    out = pd.DataFrame(rows)
    os.makedirs(OUTDIR, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT}  ({len(out)} regions)")
    print(out[["region", "known_cluster", "n_4mC_positions_T1",
               "n_6mA_positions_T1", "4mC_level_T1", "6mA_level_T1",
               "TPM_T1", "included_4mC", "included_6mA"]].to_string(index=False))
    print("\nincluded regions: 4mC",
          int((out.included_4mC == "yes").sum()),
          "/ 6mA", int((out.included_6mA == "yes").sum()))


if __name__ == "__main__":
    main()
