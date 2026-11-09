#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""120_hypomethylated_261104 step 2 - is the consistently hypomethylated subset biased?

Denominator throughout is the MEASURABLE SITES OF THE SAME RECOGNITION SEQUENCE
(not the genome, not all sites). Site-level axes (a)-(c) use Fisher's exact test;
the gene-level axis (d) uses the gene as the unit with a GC-quintile Mantel-Haenszel
companion. Every test is accompanied by its minimum detectable odds ratio, because a
null with n_hypo = 3-8 is uninformative without one.

Outputs (tables/): site_annotation_261104.tsv, bias_axis_abc_261104.tsv,
  bias_gene_kegg_cog_261104.tsv, expression_hypo_promoter_261104.tsv,
  power_floor_261104.tsv
"""
import os, re, json, resource, numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(HERE), "tables")
RNA = "/Users/okaban/bioinfo/rna-seq"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"
ANNOT = f"{RNA}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
COGF = f"{RNA}/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables/gene_COG_classification.tsv"
KEGG = f"{RNA}/15_paper_figures/scripts/newfig_261101/tables/kegg_br_sco00001_261102.keg"
TSSCSV = f"{RNA}/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
BGC = f"{RNA}/15_paper_figures/scripts/supp_261014/data/bgc_regions_29.csv"
ISL = f"{RNA}/16_genomic_islands_261102/islandpath/dimob_NC_003888.3.txt"
FC = f"{RNA}/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
CORE = (1_500_000, 7_167_507)     # project convention, 1-based
PROM = 300
CLASSES = ["TGGCCGGC_C5_4mC", "AAGCCCG_C5_4mC", "AAGCCCG_A1_6mA"]
TPS = ["T1", "T2", "T3"]


def genome_seq():
    return "".join(l.strip() for l in open(GENOME) if not l.startswith(">")).upper()


# ---------------------------------------------------------------- annotation
def build_annotation(w, g):
    ann = pd.read_csv(ANNOT, sep="\t")
    ann = ann[ann["contig"] == "NC_003888.3"].copy()
    tss = pd.read_csv(TSSCSV)
    bgc = pd.read_csv(BGC)
    isl = pd.read_csv(ISL, sep="\t", header=None, names=["id", "start", "end"])

    pos = w["pos1"].values                      # 1-based
    n = len(w)

    # (a) promoter / gene body / intergenic -- promoter wins on overlap
    body = np.zeros(n, bool); gene_of = np.array([None] * n, object)
    order = np.argsort(pos); ps = pos[order]
    for s, e, gid in zip(ann["start"], ann["end"], ann["gene_id"]):
        i, j = np.searchsorted(ps, [s, e], side="left"), None
        lo = np.searchsorted(ps, s, "left"); hi = np.searchsorted(ps, e, "right")
        idx = order[lo:hi]
        body[idx] = True
        for k in idx:
            if gene_of[k] is None:
                gene_of[k] = gid
    prom = np.zeros(n, bool); prom_gene = np.array([None] * n, object)
    for t, st, gid in zip(tss["tss_position"], tss["strand"], tss["gene_id"]):
        s, e = (t - PROM, t - 1) if st == "+" else (t + 1, t + PROM)
        lo = np.searchsorted(ps, s, "left"); hi = np.searchsorted(ps, e, "right")
        idx = order[lo:hi]
        prom[idx] = True
        for k in idx:
            if prom_gene[k] is None and isinstance(gid, str):
                prom_gene[k] = gid
    locus = np.where(prom, "promoter", np.where(body, "gene_body", "intergenic"))

    # (b) core vs arm
    region = np.where((pos >= CORE[0]) & (pos <= CORE[1]), "core", "arm")

    # (c) antiSMASH region / genomic island / rRNA operon
    def in_any(iv):
        m = np.zeros(n, bool)
        for s, e in iv:
            m |= (pos >= s) & (pos <= e)
        return m
    in_bgc = in_any(list(zip(bgc["start"], bgc["end"])))
    in_isl = in_any(list(zip(isl["start"], isl["end"])))
    rr = ann[ann["gene_biotype"] == "rRNA"].sort_values("start")
    ops, cur = [], None
    for s, e in zip(rr["start"], rr["end"]):
        if cur and s - cur[1] <= 5000:
            cur[1] = max(cur[1], e)
        else:
            cur = [s, e]; ops.append(cur)
    ops = [(s - 1000, e + 1000) for s, e in ops]      # operon +/- 1 kb
    in_rrn = in_any(ops)

    out = w.copy()
    out["locus_class"] = locus
    out["region_class"] = region
    out["in_antismash"] = in_bgc
    out["in_island"] = in_isl
    out["in_rrn_operon"] = in_rrn
    out["gene_id"] = gene_of
    out["promoter_gene_id"] = prom_gene
    out["n_rrn_operons"] = len(ops)
    # local GC in a 1 kb window, for the stratified companion
    gc = []
    for p in pos:
        s = g[max(0, p - 500):p + 500]
        gc.append((s.count("G") + s.count("C")) / max(len(s), 1))
    out["gc_1kb"] = gc
    return out


# ---------------------------------------------------------------- statistics
def min_detectable_or(n_hypo, n_bg, k_bg, alpha=0.05):
    """Smallest k among n_hypo giving a one-sided Fisher p < alpha, and its OR."""
    if n_hypo == 0:
        return np.nan, np.nan
    for k in range(0, n_hypo + 1):
        tab = [[k, n_hypo - k], [k_bg, n_bg - k_bg]]
        p = stats.fisher_exact(tab, alternative="greater")[1]
        if p < alpha:
            a, b, c, d = k, n_hypo - k, k_bg, n_bg - k_bg
            orr = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
            return k, orr
    return np.nan, np.inf


def axis_tests(ann, hypo_col, label):
    rows = []
    axes = {"locus_class": sorted(ann["locus_class"].unique()),
            "region_class": sorted(ann["region_class"].unique()),
            "in_antismash": [True], "in_island": [True], "in_rrn_operon": [True]}
    for cls in CLASSES + ["ALL_4mC", "ALL_6mA"]:
        if cls == "ALL_4mC":
            sub = ann[ann["cls"].isin(["TGGCCGGC_C5_4mC", "AAGCCCG_C5_4mC"])]
        elif cls == "ALL_6mA":
            sub = ann[ann["cls"] == "AAGCCCG_A1_6mA"]
        else:
            sub = ann[ann["cls"] == cls]
        h = sub[sub[hypo_col]]; b = sub[~sub[hypo_col]]
        for ax, levels in axes.items():
            for lv in levels:
                a1 = int((h[ax] == lv).sum()); a0 = len(h) - a1
                b1 = int((b[ax] == lv).sum()); b0 = len(b) - b1
                p = stats.fisher_exact([[a1, a0], [b1, b0]])[1]
                orr = ((a1 + .5) * (b0 + .5)) / ((a0 + .5) * (b1 + .5))
                kc, mor = min_detectable_or(len(h), len(b), b1)
                rows.append(dict(definition=label, cls=cls, axis=ax, level=str(lv),
                                 n_hypo=len(h), k_hypo=a1,
                                 pct_hypo=100 * a1 / max(len(h), 1),
                                 n_bg=len(b), k_bg=b1, pct_bg=100 * b1 / max(len(b), 1),
                                 OR=orr, p_fisher=p, k_crit=kc, min_detectable_OR=mor))
    df = pd.DataFrame(rows)
    for (d, c), grp in df.groupby(["definition", "cls"]):
        df.loc[grp.index, "p_BH"] = stats.false_discovery_control(grp["p_fisher"].values)
    return df


def mh_or(tabs):
    """Mantel-Haenszel OR + test over a list of 2x2 [[a,b],[c,d]]."""
    num = den = 0.0; E = V = S = 0.0
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
    chi = (abs(S - E) - 0.5) ** 2 / V
    return num / den, float(stats.chi2.sf(chi, 1))


def main():
    g = genome_seq()
    w = pd.read_csv(f"{TAB}/site_classification_261104.tsv", sep="\t")
    cc = w[[f"call_{t}" for t in TPS]]
    w["n_hypo_tp"] = (cc == "hypomethylated").sum(axis=1)
    w["n_meth_tp"] = (cc == "methylated").sum(axis=1)
    w["hypo_strict"] = w["n_hypo_tp"] == 3
    w["hypo_relaxed"] = (w["n_hypo_tp"] >= 1) & (w["n_meth_tp"] == 0)
    ann = build_annotation(w, g)
    ann.to_csv(f"{TAB}/site_annotation_261104.tsv", sep="\t", index=False)

    abc = pd.concat([axis_tests(ann, "hypo_strict", "strict_3of3"),
                     axis_tests(ann, "hypo_relaxed", "relaxed_never_methylated")],
                    ignore_index=True)
    abc.to_csv(f"{TAB}/bias_axis_abc_261104.tsv", sep="\t", index=False)

    # ----------------------------------------------- (d) gene-level KEGG / COG
    annot = pd.read_csv(ANNOT, sep="\t")
    annot = annot[annot["contig"] == "NC_003888.3"]
    cog = pd.read_csv(COGF, sep="\t")[["gene_id", "COG_category"]]
    # KEGG BRITE: C pathway, D sco:SCOxxxx
    kegg, cur = [], None
    for line in open(KEGG, errors="ignore"):
        if line.startswith("C"):
            cur = line[1:].strip()
        elif line.startswith("D") and cur:
            m = re.search(r"sco:(SCO\d+)", line)
            if m:
                kegg.append((m.group(1), cur))
    kegg = pd.DataFrame(kegg, columns=["old_locus_tag", "pathway"]).drop_duplicates()
    kegg = kegg.merge(annot[["gene_id", "old_locus_tag"]], on="old_locus_tag")

    # gene is in the denominator iff it carries >=1 measurable site (body or promoter)
    ann["any_gene"] = ann["gene_id"].fillna(ann["promoter_gene_id"])
    gsite = ann.dropna(subset=["any_gene"])
    gene = gsite.groupby("any_gene").agg(
        n_sites=("pos1", "size"), n_hypo_strict=("hypo_strict", "sum"),
        n_hypo_relaxed=("hypo_relaxed", "sum"), gc=("gc_1kb", "mean")).reset_index()
    gene = gene.rename(columns={"any_gene": "gene_id"})
    gene["hypo_strict"] = gene["n_hypo_strict"] > 0
    gene["hypo_relaxed"] = gene["n_hypo_relaxed"] > 0
    gene["gc_q"] = pd.qcut(gene["gc"], 5, labels=False, duplicates="drop")
    gene.to_csv(f"{TAB}/gene_level_hypo_261104.tsv", sep="\t", index=False)

    rows = []
    for label, col in [("strict_3of3", "hypo_strict"),
                       ("relaxed_never_methylated", "hypo_relaxed")]:
        for src, tabn in [("COG", cog.assign(term=cog["COG_category"])),
                          ("KEGG", kegg.assign(term=kegg["pathway"]))]:
            mm = gene.merge(tabn[["gene_id", "term"]].drop_duplicates(), on="gene_id")
            hset = set(gene.loc[gene[col], "gene_id"])
            for term, grp in mm.groupby("term"):
                gin = set(grp["gene_id"])
                a1 = len(gin & hset); a0 = len(hset - gin)
                rest = set(mm["gene_id"]) - hset
                b1 = len(gin & rest); b0 = len(rest - gin)
                if a1 + b1 < 5:
                    continue
                p = stats.fisher_exact([[a1, a0], [b1, b0]])[1]
                orr = ((a1 + .5) * (b0 + .5)) / ((a0 + .5) * (b1 + .5))
                tabs = []
                for q, sg in mm.groupby("gc_q"):
                    gi = set(sg["gene_id"])
                    tabs.append((len(gi & gin & hset), len((gi & hset) - gin),
                                 len((gi & gin) - hset), len(gi - gin - hset)))
                mor, mp = mh_or(tabs)
                kc, mdor = min_detectable_or(len(hset), len(rest), b1)
                rows.append(dict(definition=label, source=src, term=term,
                                 n_hypo_genes=len(hset), k_hypo=a1,
                                 n_bg_genes=len(rest), k_bg=b1, OR=orr, p_fisher=p,
                                 OR_MH_gc=mor, p_MH_gc=mp, min_detectable_OR=mdor))
    d = pd.DataFrame(rows)
    if len(d):
        for (lab, src), grp in d.groupby(["definition", "source"]):
            d.loc[grp.index, "p_BH"] = stats.false_discovery_control(grp["p_fisher"].values)
    d.to_csv(f"{TAB}/bias_gene_kegg_cog_261104.tsv", sep="\t", index=False)

    # ----------------------------------------------- 4. expression
    fc = pd.read_csv(FC, sep="\t", comment="#")
    cnt = fc.iloc[:, 6:].values.astype(float)
    rpk = cnt / (fc["Length"].values[:, None] / 1000)
    tpm = pd.DataFrame(rpk / rpk.sum(0) * 1e6, index=fc["Geneid"])
    tpm.columns = [f"{t}_{i}" for t in TPS for i in range(3)]
    tpm["TPM_mean"] = tpm.mean(axis=1)

    prom_sites = ann.dropna(subset=["promoter_gene_id"])
    erows = []
    for label, col in [("strict_3of3", "hypo_strict"),
                       ("relaxed_never_methylated", "hypo_relaxed")]:
        for cls in CLASSES + ["ALL"]:
            s = prom_sites if cls == "ALL" else prom_sites[prom_sites["cls"] == cls]
            gh = set(s.loc[s[col], "promoter_gene_id"])
            gm = set(s.loc[~s[col], "promoter_gene_id"]) - gh
            x = tpm.reindex(list(gh))["TPM_mean"].dropna()
            y = tpm.reindex(list(gm))["TPM_mean"].dropna()
            if len(x) == 0 or len(y) == 0:
                erows.append(dict(definition=label, cls=cls, n_hypo_genes=len(x),
                                  n_meth_genes=len(y), p_mw=np.nan)); continue
            u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
            a = u / (len(x) * len(y))                       # AUC / common-language ES
            # smallest |AUC-0.5| detectable at 80% power for these n (normal approx)
            za, zb = stats.norm.isf(0.025), stats.norm.isf(0.20)
            sd = np.sqrt((len(x) + len(y) + 1) / (12.0 * len(x) * len(y)))
            erows.append(dict(definition=label, cls=cls, n_hypo_genes=len(x),
                              n_meth_genes=len(y), median_TPM_hypo=float(x.median()),
                              median_TPM_meth=float(y.median()), AUC=a, p_mw=float(p),
                              min_detectable_AUC_80pct=float(0.5 + (za + zb) * sd)))
    pd.DataFrame(erows).to_csv(f"{TAB}/expression_hypo_promoter_261104.tsv",
                               sep="\t", index=False)

    pf = abc[["definition", "cls", "axis", "level", "n_hypo", "n_bg", "pct_bg",
              "k_crit", "min_detectable_OR"]]
    pf.to_csv(f"{TAB}/power_floor_261104.tsv", sep="\t", index=False)
    print("rRNA operons:", int(ann["n_rrn_operons"].iloc[0]))
    print(ann.groupby("cls")[["hypo_strict", "hypo_relaxed"]].sum().to_string())
    print("sig strict (p_BH<0.05):", int((abc[(abc.definition == 'strict_3of3')]["p_BH"] < .05).sum()),
          "| sig relaxed:", int((abc[(abc.definition != 'strict_3of3')]["p_BH"] < .05).sum()))
    print("peak RSS GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 3))


if __name__ == "__main__":
    main()
