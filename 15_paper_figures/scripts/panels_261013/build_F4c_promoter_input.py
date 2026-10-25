#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenerate the Figure 4b (panel_4c.py) input table from its source of truth.

    python3 build_F4c_promoter_input.py            # dry: print the diff only
    python3 build_F4c_promoter_input.py --write    # rewrite the CSV

Why this script exists
----------------------
`data_261013/v844c5b6d_F4c_promoter.csv` was a hand-staged file: nothing on disk
generated it, so a corrected value had to be pasted into it by hand and the panel
label had to be kept in step by hand as well. On 2026-10-22 the focal ratio moved
(0.9212 -> 0.9105) when the TSS-distance coordinate convention was fixed, which is
exactly the edit that hand-staging gets wrong. The focal row is now DERIVED from
the table that settled it.

All four rows (2026-10-22, second correction)
---------------------------------------------
Source: 11_epigenome_integration/analysis/48_TSS_methylation_gradient/tables/
        F4b_ladder_corrected_261022.tsv
        (written by recompute_F4b_ladder_261022.py, which drives the same
        functions that produced SuppTable22b)
Fixed by that script: tss_set = jeong_measured (the 3,570 measured dRNA-seq
        TSSs, the set the claim is made on), universe = both_A (both adenines of
        AAGCCCG), window_bp = 100, coord_convention = corrected for ALL FOUR
        rungs. `corrected` is the convention, not `ledger`: `ledger` measures the
        distance between a 0-based BED start and a 1-based TSS coordinate and is
        off by one base.

Until 2026-10-22 only the focal rung was derived; the other three were frozen
constants (0.8052 / 0.8945 / 0.8475) carried in this file, on the `ledger`
convention, so one ladder mixed two conventions that differ by one base and the
four rungs were not comparable. They are now recomputed. The three old constants
are not reproducible under EITHER convention from the pileups: the same mechanism
run with coord_convention=ledger gives 0.8075 / 0.8986 / 0.8391
(`recompute_F4b_ladder_261022.py --conv=ledger`), so they came from the retired
filtered-era table, not from a convention difference. Two of them (0.8052,
0.8945) are retracted values in the canon ledger.

The focal rung is still checked against SuppTable22b, which the value registry
and the manuscript cite.
"""
import argparse
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]                                   # 15_paper_figures/
SRC = Path("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
           "48_TSS_methylation_gradient/tables/F4b_ladder_corrected_261022.tsv")
SUPP22B = ROOT / "tables/supplementary/SuppTable22b_TSS_set_controls_261021.tsv"
OUT = ROOT / "data_261013/v844c5b6d_F4c_promoter.csv"

FOCAL_KEY = dict(tss_set="jeong_measured", universe="both_A",
                 window_bp="100", coord_convention="corrected")
FOCAL_NAME = "all reads, GC-matched"
ORDER = ["filtered, unmatched", "filtered, GC-matched",
         "all reads, unmatched", FOCAL_NAME]


def ladder():
    """The four rungs, keyed by variant name, from the recomputed table."""
    with SRC.open() as fh:
        rows = {r["variant"]: r for r in csv.DictReader(fh, delimiter="\t")}
    missing = [v for v in ORDER if v not in rows]
    if missing:
        raise SystemExit(f"{SRC.name} is missing rung(s): {missing}")
    for v in ORDER:
        if rows[v]["coord_convention"] != FOCAL_KEY["coord_convention"]:
            raise SystemExit(
                f"rung '{v}' is on coord_convention={rows[v]['coord_convention']}; "
                "all four rungs must share one convention")
    return {v: dict(ratio=float(rows[v]["ratio"]), ci_lo=float(rows[v]["ci95_lo"]),
                    ci_hi=float(rows[v]["ci95_hi"]),
                    p_value=float(rows[v]["p_value"]),
                    n_near=int(rows[v]["n_near"]), n_distal=int(rows[v]["n_distal"]),
                    n_cultures=int(rows[v]["n_cultures"])) for v in ORDER}


def check_focal(f):
    """The focal rung must still equal the SuppTable22b row the registry cites."""
    with SUPP22B.open() as fh:
        hits = [r for r in csv.DictReader(fh, delimiter="\t")
                if all(r[k] == v for k, v in FOCAL_KEY.items())]
    if len(hits) != 1:
        raise SystemExit(
            f"expected exactly 1 row matching {FOCAL_KEY} in {SUPP22B.name}, got {len(hits)}")
    r = hits[0]
    for key, col in (("ratio", "ratio"), ("ci_lo", "ci95_lo"), ("ci_hi", "ci95_hi"),
                     ("n_near", "n_near"), ("n_distal", "n_distal")):
        if float(f[key]) != float(r[col]):
            raise SystemExit(f"focal rung {key}: ladder {f[key]} vs "
                             f"{SUPP22B.name} {r[col]} — stopping without writing")


def build(rungs):
    return [{"variant": v, "ratio": f"{rungs[v]['ratio']:.4f}",
             "ci_lo": f"{rungs[v]['ci_lo']:.4f}", "ci_hi": f"{rungs[v]['ci_hi']:.4f}",
             "p_value": f"{rungs[v]['p_value']:.4g}"} for v in ORDER]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    rungs = ladder()
    f = rungs[FOCAL_NAME]
    check_focal(f)
    rows = build(rungs)
    new = ("variant,ratio,ci_lo,ci_hi,p_value\n" +
           "".join(f'"{r["variant"]}",{r["ratio"]},{r["ci_lo"]},{r["ci_hi"]},'
                   f'{r["p_value"]}\n' for r in rows))
    old = OUT.read_text() if OUT.exists() else ""

    print(f"four rungs from {SRC.name} (coord_convention=corrected):")
    for v in ORDER:
        r = rungs[v]
        print(f"  {v:<24} {r['ratio']:.4f}  95% CI {r['ci_lo']:.4f}-{r['ci_hi']:.4f}"
              f"  p {r['p_value']:.3g}")
    print(f"panel label should read {f['ratio']:.4f}; "
          f"stats line p = {f['p_value'] * 1e6:.1f} x 10^-6, n = {f['n_cultures']} culture units")
    if old == new:
        print(f"{OUT.name}: already current")
    elif args.write:
        OUT.write_text(new)
        print(f"{OUT.name}: rewritten")
    else:
        print(f"{OUT.name}: WOULD CHANGE (dry run; pass --write)")
        print("--- current ---\n" + old + "--- new ---\n" + new, end="")
