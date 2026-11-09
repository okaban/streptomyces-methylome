#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""128/04 — Q04 verdict against the pre-declared bound.

Threshold: |rho| = 0.0863, the value in the producing table for Supplementary
Table S6d (SuppTable19b_detectable_effect_sizes_261021.tsv). The task brief
phrased the pre-declared prediction as 0.088, which is the number the manuscript
text carries; 0.0863 is the stricter of the two and is what the Fisher-z
calculation on n = 1,051 actually returns, so the verdict is taken on 0.0863.

Reads only tables written by 128/02; computes nothing new from raw data.
Output: tables/verdict_261107.tsv
"""
import numpy as np, pandas as pd
from pathlib import Path

TAB = Path(__file__).resolve().parent.parent/"tables"
THR = 0.0863

full = pd.read_csv(TAB/"external_corr_full_261107.tsv", sep="\t")
head = pd.read_csv(TAB/"external_corr_headline_261107.tsv", sep="\t")

rows = []
for (ds, gs), s in head.groupby(["dataset", "gene_set"]):
    ext = ds != "THIS_STUDY"
    rows.append(dict(
        scope="pre-specified primary (level per kb, all samples pooled)",
        dataset=ds, external=ext, platform=s.platform.iloc[0], gene_set=gs,
        n_tests=len(s), n_genes=int(s.n.iloc[0]),
        max_abs_rho_raw=s.rho_raw.abs().max(),
        n_raw_over_thr=int((s.rho_raw.abs() > THR).sum()),
        max_abs_rho_partial=s.rho_partial.abs().max(),
        n_partial_over_thr=int((s.rho_partial.abs() > THR).sum()),
        min_p_partial=s.p_partial.min(),
        rho_floor_raw=s.rho_floor_raw.iloc[0],
        rho_floor_partial=s.rho_floor_partial.iloc[0],
        power_at_least_as_good=bool(s.rho_floor_partial.iloc[0] <= THR * 1.05)))
for (ds, gs), s in full.groupby(["dataset", "gene_set"]):
    rows.append(dict(
        scope="full sweep (3 exposures x every sample group)",
        dataset=ds, external=(ds != "THIS_STUDY"), platform=s.platform.iloc[0], gene_set=gs,
        n_tests=len(s), n_genes=int(s.n.median()),
        max_abs_rho_raw=s.rho_raw.abs().max(),
        n_raw_over_thr=int((s.rho_raw.abs() > THR).sum()),
        max_abs_rho_partial=s.rho_partial.abs().max(),
        n_partial_over_thr=int((s.rho_partial.abs() > THR).sum()),
        min_p_partial=s.p_partial.min(),
        rho_floor_raw=s.rho_floor_raw.median(),
        rho_floor_partial=s.rho_floor_partial.median(),
        power_at_least_as_good=bool(s.rho_floor_partial.median() <= THR * 1.05)))

v = pd.DataFrame(rows).sort_values(["scope", "gene_set", "external", "dataset"])
v.to_csv(TAB/"verdict_261107.tsv", sep="\t", index=False, float_format="%.6g")
pd.set_option("display.width", 220)
print("threshold |rho| =", THR)
print(v.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# which individual headline tests cross the threshold
print("\n--- headline tests with |rho| > %.4f ---" % THR)
for which in ("rho_raw", "rho_partial"):
    s = head[head[which].abs() > THR]
    print(f"[{which}] {len(s)} of {len(head)}")
    if len(s):
        print(s[["dataset", "gene_set", "mark", "timepoint", "n", which,
                 which.replace("rho", "p")]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))
