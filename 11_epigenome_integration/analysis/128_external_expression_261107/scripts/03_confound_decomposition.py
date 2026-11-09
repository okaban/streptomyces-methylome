#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""128/03 — why the raw gene-level correlation is non-zero, and which covariate removes it.

The headline test (128/02) gives a raw Spearman of up to |rho| = 0.090 genome-wide and
|rho| = 0.116 on the regulatory set in external data, both above the 0.088 bound the
manuscript quotes, while the covariate-adjusted value stays below 0.047 everywhere.
This script adds the covariates one at a time so the wording in the memo can say which
one carries the raw signal rather than asserting it.

Output: tables/confound_decomposition_261107.tsv
        figures/q04_external_corr_261107.png
"""
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, norm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent.parent
TAB  = HERE/"tables"; FIG = HERE/"figures"; FIG.mkdir(exist_ok=True)
COV  = ["length_bp", "gc", "dist_oriC", "gene_density"]
ZA, ZB = norm.ppf(0.975), norm.ppf(0.80)
def rho_floor(n, k=0): return float(np.tanh((ZA + ZB) / np.sqrt(n - 3 - k)))
def rank(v): return pd.Series(v).rank().values
def partial_spearman(x, y, Z):
    rx, ry = rank(x), rank(y)
    A = np.column_stack([np.ones(len(x))] + [rank(Z[:, j]) for j in range(Z.shape[1])])
    return pearsonr(rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0],
                    ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0])

# rebuild the two cells that exceed 0.088, plus the in-study counterpart
import importlib.util
spec = importlib.util.spec_from_file_location("m02", HERE/"scripts/02_external_expression_corr.py")
# re-running 02 would redo everything; instead reload its inputs directly
meth = pd.read_csv(TAB/"gene_methylation_covariates_261107.tsv", sep="\t", dtype={"old_locus_tag": str})
meth = meth[meth.old_locus_tag.notna() & (meth.old_locus_tag != "")].copy()
meth["sco"] = meth.old_locus_tag.str.strip().str.upper()
meth = meth.drop_duplicates("sco")

RS = Path("/Users/okaban/bioinfo/rna-seq")
reg = pd.read_csv(RS/"11_epigenome_integration/analysis/52_shielded_exposed_boundary/"
                     "tables/SuppTable_regulatory_gene_classification_n1051.tsv", sep="\t",
                  dtype={"old_locus_tag": str})
REG = set(reg.old_locus_tag.dropna().str.strip().str.upper())

import re
D = HERE/"data"
def tidy(df, c, cols, log=True):
    d = df[[c] + list(cols)].copy()
    d["sco"] = d[c].astype(str).str.strip().str.upper()
    d = d[d.sco.str.match(r"^SCO\d+$")]
    v = d[list(cols)].apply(pd.to_numeric, errors="coerce")
    if log: v = np.log2(v.clip(lower=0) + 1)
    return pd.concat([d[["sco"]], v], axis=1).groupby("sco", as_index=False).mean()

E = {}
d = pd.read_csv(D/"GSE132487_M145-F516-518_quantnorm.tsv.gz", sep="\t")
c = [x for x in d.columns if re.match(r"^F51[678]_\d+h$", x)]
E["GSE132487"] = tidy(d, "Identifier", c).assign(expr=lambda t: t[c].mean(axis=1))[["sco", "expr"]]
d = pd.read_csv(D/"GSE174192_Scoel_Transcription_Combined.csv.gz")
c = [x for x in d.columns if re.match(r"^Expression \d+_\d$", x)]
E["GSE174192"] = tidy(d, "Synonym", c).assign(expr=lambda t: t[c].mean(axis=1))[["sco", "expr"]]
d = pd.read_csv(RS/"05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv",
                sep="\t", dtype={"old_locus_tag": str})
c = [x for x in d.columns if x.startswith("norm_M145_")]
E["THIS_STUDY"] = tidy(d, "old_locus_tag", c).assign(expr=lambda t: t[c].mean(axis=1))[["sco", "expr"]]

CELLS = [("GSE132487", "all_chromosomal_genes", "levelkb_4mC_T1"),
         ("THIS_STUDY", "all_chromosomal_genes", "levelkb_4mC_T1"),
         ("GSE174192", "regulatory_n1051", "levelkb_6mA_T1"),
         ("THIS_STUDY", "regulatory_n1051", "levelkb_6mA_T1")]
LADDER = [[], ["gc"], ["length_bp"], ["dist_oriC"], ["gene_density"], COV]

rows = []
for ds, gs, col in CELLS:
    d = meth.merge(E[ds], on="sco", how="inner")
    if gs == "regulatory_n1051": d = d[d.sco.isin(REG)]
    sub = d[["expr", col] + COV].replace([np.inf, -np.inf], np.nan).dropna()
    x, y = sub[col].values, sub["expr"].values
    for Zc in LADDER:
        if not Zc:
            r, p = spearmanr(x, y)
        else:
            r, p = partial_spearman(x, y, sub[Zc].values)
        rows.append(dict(dataset=ds, gene_set=gs, exposure=col,
                         covariates="none" if not Zc else "+".join(Zc),
                         n=len(sub), rho=r, p=p, rho_floor=rho_floor(len(sub), len(Zc))))
    # how each covariate relates to the two variables
    for cv in COV:
        rows.append(dict(dataset=ds, gene_set=gs, exposure=col,
                         covariates=f"[check] rho({cv}, exposure)", n=len(sub),
                         rho=spearmanr(sub[cv], x)[0], p=spearmanr(sub[cv], x)[1], rho_floor=np.nan))
        rows.append(dict(dataset=ds, gene_set=gs, exposure=col,
                         covariates=f"[check] rho({cv}, expression)", n=len(sub),
                         rho=spearmanr(sub[cv], y)[0], p=spearmanr(sub[cv], y)[1], rho_floor=np.nan))

out = pd.DataFrame(rows)
out.to_csv(TAB/"confound_decomposition_261107.tsv", sep="\t", index=False, float_format="%.6g")
print(out[~out.covariates.str.startswith("[check]")].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
print()
print(out[out.covariates.str.startswith("[check]")].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# ---------------------------------------------------------------- figure
h = pd.read_csv(TAB/"external_corr_headline_261107.tsv", sep="\t")
order = ["GSE132487", "GSE174192", "GSE158810", "GSE58666", "THIS_STUDY"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
for ax, gs, flo in zip(axes, ["all_chromosomal_genes", "regulatory_n1051"], [None, None]):
    s = h[h.gene_set == gs]
    yt, lab = [], []
    k = 0
    for ds in order:
        for mark in ["4mC", "6mA"]:
            for tp in ["T1", "T2", "T3"]:
                r = s[(s.dataset == ds) & (s["mark"] == mark) & (s.timepoint == tp)]
                if r.empty: continue
                col = "#888888" if ds == "THIS_STUDY" else ("#1f78b4" if mark == "4mC" else "#e31a1c")
                ax.plot(r.rho_raw.iloc[0], k, "o", ms=5, mfc="none", color=col)
                ax.plot(r.rho_partial.iloc[0], k, "o", ms=5, color=col)
                yt.append(k); lab.append(f"{ds.replace('THIS_STUDY','this study')} {mark} {tp}")
                k += 1
        k += 0.6
    ax.axvline(0, color="k", lw=0.6)
    for v in (-0.0863, 0.0863):
        ax.axvline(v, color="k", ls="--", lw=0.8)
    fl = s.rho_floor_partial.median()
    for v in (-fl, fl):
        ax.axvline(v, color="0.6", ls=":", lw=0.8)
    ax.set_xlim(-0.16, 0.16); ax.set_xlabel("Spearman rho (methylation level per kb vs transcript level)")
    ax.set_title(gs.replace("_", " "), fontsize=10)
axes[0].set_yticks(yt); axes[0].set_yticklabels(lab, fontsize=6)
axes[0].invert_yaxis()
fig.suptitle("Open circles: raw Spearman.  Filled: partial, controlling length, GC, |d(oriC)|, local gene density.\n"
             "Dashed: +/-0.0863, this study's 80%-power bound (Supplementary Table S6d).  Dotted: the dataset's own bound.",
             fontsize=8)
fig.tight_layout(rect=[0, 0, 1, 0.9])
fig.savefig(FIG/"q04_external_corr_261107.png", dpi=200)
print("wrote", FIG/"q04_external_corr_261107.png")
