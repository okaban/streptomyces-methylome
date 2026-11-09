#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""128/02 — Q04: does the gene-level methylation/expression null reproduce in
independent, externally generated expression data?

    python3 02_external_expression_corr.py

Prediction, declared before the result was seen (task brief, 2026-11-07):
if this study's conclusion holds, no external dataset gives a rank correlation
between gene-body methylation level and transcript level whose |rho| exceeds
0.088 — the smallest effect this study's gene-level design resolves at 80%
power (Supplementary Table S6d; the table's own value is 0.0863 on n = 1,051).

Statistic — the same as this study's gene-level test
----------------------------------------------------
  Spearman rho between gene-body methylation level and transcript level, and
  the partial Spearman controlling the four covariates the manuscript declares
  for gene-level correlations: gene length, gene GC, |distance from oriC| and
  local gene density. Partial Spearman = Pearson r between the residuals of
  rank(x) and rank(y) after OLS on the ranks of the covariates — the identical
  estimator used in 80_partial_corr_covariates/partial_corr_covariates.py.

Power
-----
  Fisher z, two-sided, alpha = 0.05, power = 0.80:
      |rho|_min = tanh((z_{0.975} + z_{0.80}) / sqrt(n - 3 - k))
  k = 0 for the raw correlation, 4 for the partial. Reported per dataset and
  per gene set, so the external floor can be put beside this study's 0.0863.

Inputs
------
  methylation : ../tables/gene_methylation_covariates_261107.tsv  (128/01)
  regulatory  : 52_shielded_exposed_boundary/tables/SuppTable_regulatory_gene_classification_n1051.tsv
  internal    : 05_annotation/.../gene_master_DESeq2.tsv   (this study's own RNA-seq,
                norm_M145_* columns; run only as the in-study anchor, clearly labelled)
  external    : ../data/*  downloaded from GEO on 2026-11-07, see 00_fetch notes
"""
import gzip, re, sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, pearsonr, norm

RS   = Path("/Users/okaban/bioinfo/rna-seq")
HERE = Path(__file__).resolve().parent.parent
DATA = HERE/"data"
TAB  = HERE/"tables"; TAB.mkdir(exist_ok=True)

ZA, ZB = norm.ppf(0.975), norm.ppf(0.80)
def rho_floor(n, k=0):
    """smallest |rho| detectable at 80% power, alpha = 0.05, two-sided"""
    d = n - 3 - k
    return float(np.tanh((ZA + ZB) / np.sqrt(d))) if d > 0 else np.nan

def rank(v): return pd.Series(v).rank().values

def partial_spearman(x, y, Z):
    rx, ry = rank(x), rank(y)
    A = np.column_stack([np.ones(len(x))] + [rank(Z[:, j]) for j in range(Z.shape[1])])
    bx = np.linalg.lstsq(A, rx, rcond=None)[0]; by = np.linalg.lstsq(A, ry, rcond=None)[0]
    return pearsonr(rx - A @ bx, ry - A @ by)

# ---------------------------------------------------------------- methylation
meth = pd.read_csv(TAB/"gene_methylation_covariates_261107.tsv", sep="\t",
                   dtype={"old_locus_tag": str})
meth = meth[meth.old_locus_tag.notna() & (meth.old_locus_tag != "")].copy()
meth["sco"] = meth.old_locus_tag.str.strip().str.upper()
meth = meth.drop_duplicates("sco")
COV = ["length_bp", "gc", "dist_oriC", "gene_density"]
print("methylation table: %d genes with an SCO id" % len(meth))

reg = pd.read_csv(RS/"11_epigenome_integration/analysis/52_shielded_exposed_boundary/"
                     "tables/SuppTable_regulatory_gene_classification_n1051.tsv", sep="\t",
                  dtype={"old_locus_tag": str})
REG = set(reg.old_locus_tag.dropna().str.strip().str.upper())
print("regulatory set:", len(REG))

# ---------------------------------------------------------------- expression loaders
def tidy(df, sco_col, value_cols, log=True):
    """-> DataFrame[sco, <sample columns>] with one row per SCO id"""
    d = df[[sco_col] + list(value_cols)].copy()
    d["sco"] = d[sco_col].astype(str).str.strip().str.upper()
    d = d[d.sco.str.match(r"^SCO\d+$")]
    v = d[list(value_cols)].apply(pd.to_numeric, errors="coerce")
    if log: v = np.log2(v.clip(lower=0) + 1)
    out = pd.concat([d[["sco"]], v], axis=1).groupby("sco", as_index=False).mean()
    return out

DATASETS = {}

# --- GSE132487 : M145, 3 biological replicate cultures x 9 sampling times, RNA-seq,
#     submitter-quantile-normalised gene matrix (Sulheim et al., SYSTERACT).
d = pd.read_csv(DATA/"GSE132487_M145-F516-518_quantnorm.tsv.gz", sep="\t")
cols = [c for c in d.columns if re.match(r"^F51[678]_\d+h$", c)]
e = tidy(d, "Identifier", cols)
groups = {}
for h in sorted({c.split("_")[1] for c in cols}, key=lambda s: int(s[:-1])):
    groups[h] = [c for c in cols if c.endswith("_" + h)]
DATASETS["GSE132487"] = dict(expr=e, groups=groups, platform="RNA-seq",
    strain="M145", design="3 biological replicate cultures (F516/F517/F518) x 9 times (21-57 h)",
    value="submitter quantile-normalised gene-level expression, log2(x+1)")

# --- GSE174192 : M145, ISP2 liquid, 8 days, biological triplicate, RNA-seq.
d = pd.read_csv(DATA/"GSE174192_Scoel_Transcription_Combined.csv.gz")
cols = [c for c in d.columns if re.match(r"^Expression \d+_\d$", c)]
e = tidy(d, "Synonym", cols)
groups = {}
for h in sorted({c.split(" ")[1].split("_")[0] for c in cols}, key=int):
    groups[h + "h"] = [c for c in cols if c.startswith("Expression %s_" % h)]
DATASETS["GSE174192"] = dict(expr=e, groups=groups, platform="RNA-seq",
    strain="M145", design="3 biological replicates x 8 times (24-192 h), ISP2 medium",
    value="submitter gene-level read counts, log2(x+1)")

# --- GSE158810 : M145, 50 h, untreated control arm only (3 biological replicates), RNA-seq.
d = pd.read_csv(DATA/"GSE158810_RNAseq_cAMP_supplemented_M145_processed_data.txt.gz", sep="\t")
cols = ["M145_50h_1", "M145_50h_2", "M145_50h_3"]
e = tidy(d, "Locus_tag", cols)
DATASETS["GSE158810"] = dict(expr=e, groups={"50h": cols}, platform="RNA-seq",
    strain="M145", design="3 biological replicates, 50 h, untreated control arm only "
                          "(the cAMP-supplemented arm is excluded)",
    value="submitter gene-level read counts, log2(x+1)")

# --- GSE58666 : M145 wild-type arm, 5 culture times x 3 biological replicates,
#     two-colour array against a genomic-DNA common reference (VALUE = log2 cDNA/gDNA,
#     i.e. a transcript-level measure, not a between-condition ratio).
def load_series_matrix(path):
    with gzip.open(path, "rt") as fh:
        txt = fh.read()
    hdr = {}
    for ln in txt.splitlines():
        if ln.startswith("!Sample_title"):
            hdr["title"] = [t.strip('"') for t in ln.split("\t")[1:]]
        if ln.startswith("!Sample_geo_accession"):
            hdr["gsm"] = [t.strip('"') for t in ln.split("\t")[1:]]
    body = txt.split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    from io import StringIO
    m = pd.read_csv(StringIO(body), sep="\t")
    m.columns = [c.strip('"') for c in m.columns]
    m[m.columns[0]] = m[m.columns[0]].astype(str).str.strip('"')
    return m, hdr

m, hdr = load_series_matrix(DATA/"GSE58666_series_matrix.txt.gz")
t2g = dict(zip(hdr["gsm"], hdr["title"]))
wt = [g for g in hdr["gsm"] if t2g[g].startswith("SI.")]       # SI.* = M145 wild type
# probe -> ORF from the platform table
with gzip.open(DATA/"GPL16698_family.soft.gz", "rt") as fh:
    txt = fh.read().split("!platform_table_begin\n")[1].split("!platform_table_end")[0]
from io import StringIO
plat = pd.read_csv(StringIO(txt), sep="\t", usecols=["ID", "ORF"])
m = m.rename(columns={m.columns[0]: "ID"}).merge(plat, on="ID", how="inner")
e = tidy(m, "ORF", wt, log=False)       # VALUE already log2
groups = {}
for h in ["32", "42", "49", "56", "66"]:
    groups[h + "h"] = [g for g in wt if t2g[g].startswith("SI.%s." % h)]
DATASETS["GSE58666"] = dict(expr=e, groups=groups, platform="two-colour microarray",
    strain="M145", design="M145 wild-type arm only: 5 culture times (32-66 h) x 3 biological "
                          "replicates; the DargR arm is excluded",
    value="submitter VALUE = log2(cDNA / genomic-DNA common reference), probes averaged per ORF")

# --- in-study anchor (NOT external): this study's own RNA-seq, same pipeline.
d = pd.read_csv(RS/"05_annotation/analysis/05_annotation_260128_v1/tables/gene_master_DESeq2.tsv",
                sep="\t", dtype={"old_locus_tag": str})
cols = [c for c in d.columns if c.startswith("norm_M145_")]
e = tidy(d, "old_locus_tag", cols)
groups = {"T1": [c for c in cols if c.startswith("norm_M145_1")],
          "T2": [c for c in cols if c.startswith("norm_M145_2")],
          "T3": [c for c in cols if c.startswith("norm_M145_3")]}
DATASETS["THIS_STUDY"] = dict(expr=e, groups=groups, platform="RNA-seq",
    strain="M145", design="this study: 3 biological replicates x 3 times (12/24/50 h), R5 medium",
    value="DESeq2 size-factor-normalised counts, log2(x+1)")

# ---------------------------------------------------------------- the tests
MARKS = ["4mC", "6mA"]; TPS = ["T1", "T2", "T3"]
EXPOS = {"levelkb": "sum of canonical-site levels per kb of gene body",
         "n":       "number of canonical sites per gene body",
         "meanlevel": "mean level over the gene's canonical sites (genes with >=1 site only)"}

rows = []
for ds, meta in DATASETS.items():
    e = meta["expr"]
    sample_groups = dict(meta["groups"])
    sample_groups["ALL"] = sorted({c for v in meta["groups"].values() for c in v})
    for gset, gname in ((None, "all_chromosomal_genes"), (REG, "regulatory_n1051")):
        for grp, cols in sample_groups.items():
            y_all = e[["sco"]].copy()
            y_all["expr"] = e[cols].mean(axis=1)
            d = meth.merge(y_all, on="sco", how="inner")
            if gset is not None:
                d = d[d.sco.isin(gset)]
            for mark in MARKS:
                for tp in TPS:
                    for ex, exdesc in EXPOS.items():
                        col = {"levelkb": f"levelkb_{mark}_{tp}", "n": f"n_{mark}_{tp}",
                               "meanlevel": f"meanlevel_{mark}_{tp}"}[ex]
                        sub = d[["expr", col] + COV].replace([np.inf, -np.inf], np.nan).dropna()
                        n = len(sub)
                        if n < 20: continue
                        x, y = sub[col].values, sub["expr"].values
                        if ex == "n":
                            x = x / (d.loc[sub.index, "length_bp"].values / 1000.0)
                        r0, p0 = spearmanr(x, y)
                        rp, pp = partial_spearman(x, y, sub[COV].values)
                        rows.append(dict(dataset=ds, platform=meta["platform"],
                                         gene_set=gname, sample_group=grp,
                                         mark=mark, timepoint=tp, exposure=ex,
                                         n=n, rho_raw=r0, p_raw=p0,
                                         rho_partial=rp, p_partial=pp,
                                         rho_floor_raw=rho_floor(n, 0),
                                         rho_floor_partial=rho_floor(n, len(COV))))
    print("done", ds, flush=True)

res = pd.DataFrame(rows)
res.to_csv(TAB/"external_corr_full_261107.tsv", sep="\t", index=False, float_format="%.6g")
print("wrote", TAB/"external_corr_full_261107.tsv", res.shape)

# headline table: primary exposure, pooled over all samples of the dataset
head = res[(res.exposure == "levelkb") & (res.sample_group == "ALL")].copy()
head = head.sort_values(["gene_set", "dataset", "mark", "timepoint"])
head.to_csv(TAB/"external_corr_headline_261107.tsv", sep="\t", index=False, float_format="%.6g")

# per-dataset summary: worst case across every test run for that dataset
summ = []
for (ds, gs), sub in res.groupby(["dataset", "gene_set"]):
    i = sub.rho_raw.abs().idxmax(); j = sub.rho_partial.abs().idxmax()
    summ.append(dict(dataset=ds, platform=sub.platform.iloc[0], gene_set=gs,
                     n_tests=len(sub), n_genes_median=int(sub.n.median()),
                     max_abs_rho_raw=abs(sub.rho_raw[i]),
                     at_raw=f"{sub.mark[i]} {sub.timepoint[i]} {sub.exposure[i]} {sub.sample_group[i]}",
                     max_abs_rho_partial=abs(sub.rho_partial[j]),
                     at_partial=f"{sub.mark[j]} {sub.timepoint[j]} {sub.exposure[j]} {sub.sample_group[j]}",
                     rho_floor_raw=rho_floor(int(sub.n.median()), 0),
                     rho_floor_partial=rho_floor(int(sub.n.median()), len(COV)),
                     exceeds_0088=bool(max(abs(sub.rho_raw[i]), abs(sub.rho_partial[j])) > 0.088)))
summ = pd.DataFrame(summ).sort_values(["gene_set", "dataset"])
summ.to_csv(TAB/"external_corr_summary_261107.tsv", sep="\t", index=False, float_format="%.6g")
print(summ.to_string(index=False))

# dataset provenance
prov = pd.DataFrame([dict(dataset=k, platform=v["platform"], strain=v["strain"],
                          design=v["design"], value=v["value"],
                          n_genes_matched=int(meth.merge(v["expr"][["sco"]], on="sco").shape[0]))
                     for k, v in DATASETS.items()])
prov.to_csv(TAB/"datasets_used_261107.tsv", sep="\t", index=False)
print(prov.to_string(index=False))
print("\nstudy floor S6d (n=1051, k=0): %.4f" % rho_floor(1051, 0))
