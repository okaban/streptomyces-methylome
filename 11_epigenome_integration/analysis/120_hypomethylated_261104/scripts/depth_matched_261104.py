#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""120_hypomethylated_261104 step 3 - the core/arm and promoter signals under depth matching.

Core sites carry roughly twice the pooled depth of arm sites, so a determinate call
(hypomethylated rather than indeterminate) is easier to make in the core. Two controls:
  (1) Mantel-Haenszel stratified on pooled-depth quintile, all sites retained;
  (2) uniform downsampling of every timepoint to a common depth D, reclassification
      with the same criterion, then the same Fisher test.
Outputs: depth_matched_tests_261104.tsv
"""
import os, numpy as np, pandas as pd, resource
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(HERE), "tables")
TPS = ["T1", "T2", "T3"]
CLASSES = ["TGGCCGGC_C5_4mC", "AAGCCCG_C5_4mC", "AAGCCCG_A1_6mA"]
ALPHA = 0.05
TAU_FRAC = 0.5


def cp(k, n, alpha=ALPHA):
    k = np.asarray(k, float); n = np.asarray(n, float)
    lo = np.where(k > 0, stats.beta.ppf(alpha / 2, np.maximum(k, 1e-9), n - k + 1), 0.0)
    hi = np.where(k < n, stats.beta.isf(alpha / 2, k + 1, np.maximum(n - k, 1e-9)), 1.0)
    return np.nan_to_num(lo, 0.0), np.nan_to_num(hi, 1.0)


def mh(tabs):
    num = den = E = V = S = 0.0
    for a, b, c, d in tabs:
        n = a + b + c + d
        if n == 0:
            continue
        num += a * d / n; den += b * c / n
        r1, r2, c1 = a + b, c + d, a + c
        E += r1 * c1 / n
        if n > 1:
            V += r1 * r2 * c1 * (n - c1) / (n * n * (n - 1))
        S += a
    if den == 0 or V == 0:
        return np.nan, np.nan
    return num / den, float(stats.chi2.sf((abs(S - E) - 0.5) ** 2 / V, 1))


def main():
    ann = pd.read_csv(f"{TAB}/site_annotation_261104.tsv", sep="\t")
    ann["depth_mean"] = ann[[f"Nvalid_{t}" for t in TPS]].mean(axis=1)
    rng = np.random.default_rng(20261104)
    rows = []
    AXES = [("region_class", "core"), ("locus_class", "promoter"),
            ("locus_class", "gene_body"), ("in_island", True)]

    for cls in CLASSES + ["ALL_4mC"]:
        sub = (ann[ann["cls"].isin(["TGGCCGGC_C5_4mC", "AAGCCCG_C5_4mC"])] if cls == "ALL_4mC"
               else ann[ann["cls"] == cls]).copy()
        sub["dq"] = pd.qcut(sub["depth_mean"], 5, labels=False, duplicates="drop")
        for hyp in ["hypo_strict", "hypo_relaxed"]:
            for ax, lv in AXES:
                tabs = []
                for q, g in sub.groupby("dq"):
                    h = g[g[hyp]]; b = g[~g[hyp]]
                    tabs.append(((h[ax] == lv).sum(), (~(h[ax] == lv)).sum(),
                                 (b[ax] == lv).sum(), (~(b[ax] == lv)).sum()))
                o, p = mh(tabs)
                rows.append(dict(control="MH_depth_quintile", cls=cls, definition=hyp,
                                 axis=ax, level=str(lv), D=np.nan,
                                 n_hypo=int(sub[hyp].sum()), OR=o, p=p))

        # ---- uniform downsampling control -------------------------------
        for D in [40, 60, 80]:
            m = (sub[[f"Nvalid_{t}" for t in TPS]] >= D).all(axis=1)
            s2 = sub[m].copy()
            if len(s2) < 50:
                continue
            calls = {}
            for t in TPS:
                k = rng.hypergeometric(s2[f"Nmod_{t}"].values.astype(int),
                                       (s2[f"Nvalid_{t}"] - s2[f"Nmod_{t}"]).values.astype(int),
                                       D)
                tau = TAU_FRAC * k.sum() / (D * len(s2))
                lo, hi = cp(k, np.full(len(s2), D))
                calls[t] = np.where(hi < tau, "H", np.where(lo > tau, "M", "I"))
            C = np.vstack([calls[t] for t in TPS])
            s2["hypo_strict"] = (C == "H").all(0)
            s2["hypo_relaxed"] = (C == "H").any(0) & ~(C == "M").any(0)
            for hyp in ["hypo_strict", "hypo_relaxed"]:
                for ax, lv in AXES:
                    h = s2[s2[hyp]]; b = s2[~s2[hyp]]
                    a1 = int((h[ax] == lv).sum()); a0 = len(h) - a1
                    b1 = int((b[ax] == lv).sum()); b0 = len(b) - b1
                    p = stats.fisher_exact([[a1, a0], [b1, b0]])[1]
                    orr = ((a1 + .5) * (b0 + .5)) / ((a0 + .5) * (b1 + .5))
                    rows.append(dict(control=f"downsample_D{D}", cls=cls, definition=hyp,
                                     axis=ax, level=str(lv), D=D, n_sites_retained=len(s2),
                                     n_hypo=len(h), k_hypo=a1, OR=orr, p=p))
    df = pd.DataFrame(rows)
    df.to_csv(f"{TAB}/depth_matched_tests_261104.tsv", sep="\t", index=False)
    q = df[(df.definition == "hypo_relaxed") & (df.axis == "region_class")]
    print(q[["control", "cls", "n_hypo", "OR", "p"]].round(4).to_string(index=False))
    print("peak RSS GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 3))


if __name__ == "__main__":
    main()
