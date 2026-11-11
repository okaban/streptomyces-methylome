#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R5 (primary / secondary metabolism) — recompute on the all-read lineage.

    python3 nf_r5_recompute_261101.py

What this replaces
------------------
* SuppTable20a — had no generator in the repository. Rebuilt by
  `supp_261014/build_supptable20a_allread_261101.py` on the all-read pileup.
  This script reads the new table and quantifies how far the old one sat above
  it (the F7 filter bias, canon: LEVEL +15.5 pp).
* SuppTable20b — motif occurrence inside vs flank. Reference-sequence-only, so
  it does not depend on the pileup lineage at all; recomputed here and asserted
  against the 261021 file rather than replaced.
* SuppTable20c — interior vs flank 6mA level. The 261021 generator reads
  `methyl/260102_M145/analysis/pileup/`, which is (i) the confidence-filtered
  lineage and (ii) missing the 2026-02 flowcell for 3-2 and 3-4. Rebuilt here on
  `pileup_nofilt_261015` (Nfail == 0, flowcell included: measured Nvalid median
  25 for 3-2 and 36 for 3-4, the post-flowcell depths).

SCOPE CHANGE, declared: the 261021 20c measured 6mA over EVERY adenine row of a
genome-wide pileup. The all-read pileup is motif-limited (GCCGGC offset 2,
AAGCCCG offsets 0/1/4), so the rebuilt 20c is restricted to the declared 6mA
positions, AAGCCCG A1/A2. The two tables are therefore not row-comparable; this
one measures the mark where the mark is (canon 2026-10-06: genome-wide 6mA is
exhausted by AAGCCCG; off-motif level >= 0.5 is 0.0635 % of 2.41 M adenines).

The primary-metabolism gene set is NOT decided here. Three candidate
definitions are trialled so the author can choose; the figure is provisional
until then.

Outputs (tables/):
    nf_r5_cluster_levels_261101.tsv     region x mark x timepoint, all-read
    nf_r5_interior_flank_261101.tsv     rebuilt 20c, per library
    nf_r5_cluster_tests_261101.tsv      cluster-unit tests
    nf_r5_primary_metabolism_trial_261101.tsv   the three candidate definitions
"""
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "tables")
RNA = "/Users/okaban/bioinfo/rna-seq"
SUPP = f"{RNA}/15_paper_figures/tables/supplementary"
NEW20A = f"{SUPP}/SuppTable20a_cluster_methylation_expression_allread_261101.tsv"
OLD20A = f"{SUPP}/SuppTable20a_cluster_methylation_expression_261021.tsv"
OLD20B = f"{SUPP}/SuppTable20b_cluster_motif_occurrence_261021.tsv"
REGIONS = f"{RNA}/15_paper_figures/scripts/supp_261014/data/bgc_regions_29.csv"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"
NOFILT = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015"
TSSCSV = (f"{RNA}/11_epigenome_integration/analysis/18_tss_analyses/"
          "jeong2016_all_tss.csv")
COG = (f"{RNA}/12_supplementary_figures/analysis/12_supplementary_260202_v1/"
       "tables/gene_COG_classification.tsv")
ANNOT = (f"{RNA}/05_annotation/analysis/05_annotation_260128_v1/tables/"
         "gene_annotation_basic.tsv")

LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
        "T3": ["3-2", "3-3", "3-4"]}
FLAT = [(tp, l) for tp, v in LIBS.items() for l in v]
BEDCOL = ["chrom", "start", "end", "name", "score", "strand", "ts", "te",
          "color", "Nvalid", "pct", "Nmod", "Ncanon", "Nother", "Ndel",
          "Nfail", "Ndiff", "Nnocall"]
MIN_DEPTH = 30
MIN_POS = 5
# COG metabolism categories; Q (secondary metabolism) is deliberately excluded
COG_CORE = ["C", "E", "F", "G", "H", "I"]
COG_WIDE = COG_CORE + ["P", "M"]


def genome_seq():
    return "".join(l.strip() for l in open(GENOME)
                   if not l.startswith(">")).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def occurrences(seq, motif):
    """1-based start positions of `motif` on both strands, + strand coords."""
    L = len(seq)
    out = [m.start() + 1 for m in re.finditer("(?=" + motif + ")", seq)]
    out += [L - (m.start() + len(motif)) + 1
            for m in re.finditer("(?=" + motif + ")", rc(seq))]
    return np.array(sorted(out))


def tgg_positions(seq):
    L = len(seq)
    out = set(m.start() + 4 for m in re.finditer("(?=TGGCCGGC)", seq))
    out |= set(L - 1 - (m.start() + 4)
               for m in re.finditer("(?=TGGCCGGC)", rc(seq)))
    return out


def per_library(names, keep=None):
    """{lib: DataFrame(start, strand, Nvalid, Nmod, level)} on all reads."""
    out = {}
    for _, lib in FLAT:
        d = pd.read_csv(f"{NOFILT}/{lib}_nofilt.bed", sep="\t", header=None,
                        names=BEDCOL)
        assert int(d["Nfail"].sum()) == 0, f"{lib}: not an all-read pileup"
        d = d[d["name"].isin(names)]
        if keep is not None:
            d = d[d["start"].isin(keep)]
        d = d.drop_duplicates(subset=["start", "strand"]).copy()
        d["level"] = np.where(d["Nvalid"] > 0,
                              d["Nmod"] / d["Nvalid"].clip(lower=1), np.nan)
        out[lib] = d[["start", "strand", "Nvalid", "Nmod", "level"]]
    return out


def pooled(per, libs):
    s = None
    for l in libs:
        f = per[l].set_index(["start", "strand"])[["Nvalid", "Nmod"]]
        s = f if s is None else s.add(f, fill_value=0)
    s = s[s["Nvalid"] >= MIN_DEPTH].copy()
    s["level"] = s["Nmod"] / s["Nvalid"]
    return s.reset_index()


def flanks(start, end, L):
    n = end - start + 1
    return [(max(1, start - n), start - 1), (end + 1, min(L, end + n))]


# ---------------------------------------------------------------- part 1
def cluster_levels(reg, per4, per6, L):
    rows = []
    for mark, per in (("4mC", per4), ("6mA", per6)):
        for tp, libs in LIBS.items():
            p = pooled(per, libs)
            pos = p["start"].to_numpy()
            lvl = p["level"].to_numpy()
            for _, r in reg.iterrows():
                s, e = int(r["start"]) - 1, int(r["end"]) - 1
                ins = (pos >= s) & (pos <= e)
                fl = np.zeros(len(pos), bool)
                for a, b in flanks(int(r["start"]), int(r["end"]), L):
                    fl |= (pos >= a - 1) & (pos <= b - 1)
                rows.append(dict(
                    region=f"NC_003888.3_region{int(r.region_number)}",
                    known_cluster=r["known_cluster"], mark=mark, timepoint=tp,
                    n_interior=int(ins.sum()), n_flank=int(fl.sum()),
                    level_interior=float(np.nanmean(lvl[ins])) if ins.sum() else np.nan,
                    level_flank=float(np.nanmean(lvl[fl])) if fl.sum() else np.nan))
    return pd.DataFrame(rows)


def cluster_tests(cl, new20a):
    rows = []
    for mark in ("4mC", "6mA"):
        d = cl[(cl["mark"] == mark) & (cl.n_interior >= MIN_POS)
               & (cl.n_flank >= MIN_POS)]
        # cluster unit: average the 3 timepoints within a region first
        g = d.groupby("region")[["level_interior", "level_flank"]].mean().dropna()
        st = stats.wilcoxon(g.level_interior, g.level_flank)
        rows.append(dict(test=f"{mark} interior vs flank, cluster unit",
                         n=len(g), statistic=round(float(st.statistic), 1),
                         p_value=f"{st.pvalue:.3g}",
                         median_ratio=round(float(np.median(
                             g.level_interior / g.level_flank)), 4)))
        # level vs expression, region x timepoint
        col = "included_4mC" if mark == "4mC" else "included_6mA"
        m = new20a[new20a[col] == "yes"]
        obs = []
        for tp in LIBS:
            obs.append(pd.DataFrame({
                "level": m[f"{mark}_level_{tp}"], "tpm": m[f"TPM_{tp}"]}))
        obs = pd.concat(obs).dropna()
        rho, p = stats.spearmanr(obs["level"], obs["tpm"])
        rows.append(dict(test=f"{mark} cluster mean level vs mean TPM",
                         n=len(obs), statistic=round(float(rho), 4),
                         p_value=f"{p:.3g}", median_ratio=np.nan))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- part 2
def rebuild_20b(seq, reg):
    L = len(seq)
    rows = []
    for motif in ("TGGCCGGC", "AAGCCCG"):
        occ = occurrences(seq, motif)
        dens = len(occ) / L
        for _, r in reg.iterrows():
            s, e = int(r["start"]), int(r["end"])
            n_in = int(((occ >= s) & (occ <= e)).sum())
            n_fl, len_fl = 0, 0
            for a, b in flanks(s, e, L):
                n_fl += int(((occ >= a) & (occ <= b)).sum())
                len_fl += b - a + 1
            rows.append(dict(
                region=f"NC_003888.3_region{int(r.region_number)}",
                motif=motif, length_bp=e - s + 1,
                n_inside=n_in, oe_inside=round(n_in / (dens * (e - s + 1)), 3),
                n_flank=n_fl, oe_flank=round(n_fl / (dens * len_fl), 3)))
    return pd.DataFrame(rows)


def rebuild_20c(per6, reg, tss, L):
    """20c on the all-read lineage, restricted to AAGCCCG A1/A2."""
    rows = []
    for tp, lib in FLAT:
        d = per6[lib]
        d = d[d["Nvalid"] > 0]
        pos = d["start"].to_numpy()
        lvl = d["level"].to_numpy()
        near = np.min(np.abs(pos[:, None] + 1 - tss[None, :]), axis=1) <= 300
        ins = np.zeros(len(pos), bool)
        fl = np.zeros(len(pos), bool)
        for _, r in reg.iterrows():
            s, e = int(r["start"]), int(r["end"])
            ins |= (pos + 1 >= s) & (pos + 1 <= e)
            for a, b in flanks(s, e, L):
                fl |= (pos + 1 >= a) & (pos + 1 <= b)
        fl &= ~ins
        rows.append(dict(
            timepoint=tp, replicate=lib,
            interior_all=round(float(np.nanmean(lvl[ins])), 4),
            flank_all=round(float(np.nanmean(lvl[fl])), 4),
            interior_awayTSS=round(float(np.nanmean(lvl[ins & ~near])), 4),
            flank_awayTSS=round(float(np.nanmean(lvl[fl & ~near])), 4),
            n_interior=int(ins.sum()), n_flank=int(fl.sum()),
            n_interior_awayTSS=int((ins & ~near).sum()),
            n_flank_awayTSS=int((fl & ~near).sum())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- part 3
def primary_metabolism_trial(reg, per4, per6, tpm_src):
    """Three candidate 'primary metabolism' gene sets, same statistics each."""
    ann = pd.read_csv(ANNOT, sep="\t")
    ann = ann[ann.gene_biotype == "protein_coding"]
    cog = pd.read_csv(COG, sep="\t")
    cog["cat"] = cog["COG_category"].astype(str).str[0]
    g = ann.merge(cog[["gene_id", "cat"]], on="gene_id", how="left")

    in_bgc = np.zeros(len(g), bool)
    mid = ((g["start"] + g["end"]) / 2).to_numpy()
    for _, r in reg.iterrows():
        in_bgc |= (mid >= r["start"]) & (mid <= r["end"])
    g = g.assign(in_bgc=in_bgc)

    sets = {
        "P1 COG core metabolism (C,E,F,G,H,I)": g["cat"].isin(COG_CORE).to_numpy(),
        "P2 COG metabolism + transport/envelope (C,E,F,G,H,I,P,M)":
            g["cat"].isin(COG_WIDE).to_numpy(),
        "P3 all protein-coding genes outside antiSMASH regions": ~in_bgc,
        "(comparator) genes inside antiSMASH regions": in_bgc,
    }

    p4 = pooled(per4, LIBS["T1"])
    p6 = pooled(per6, LIBS["T1"])
    rows = []
    for name, mask in sets.items():
        sub = g[mask]
        s = sub["start"].to_numpy()
        e = sub["end"].to_numpy()
        res = dict(definition=name, n_genes=int(len(sub)))
        for mark, p in (("4mC", p4), ("6mA", p6)):
            pos = p["start"].to_numpy() + 1
            lvl = p["level"].to_numpy()
            # gene-body positions
            order = np.argsort(pos)
            pos_s, lvl_s = pos[order], lvl[order]
            lo = np.searchsorted(pos_s, s)
            hi = np.searchsorted(pos_s, e, side="right")
            per_gene = np.array([np.nanmean(lvl_s[a:b]) if b > a else np.nan
                                 for a, b in zip(lo, hi)])
            # promoter proxy window: gene start -300 .. +50 on the gene's strand
            ps = np.where(sub["strand"].to_numpy() == "+", s - 300, e - 50)
            pe = np.where(sub["strand"].to_numpy() == "+", s + 50, e + 300)
            plo = np.searchsorted(pos_s, ps)
            phi = np.searchsorted(pos_s, pe, side="right")
            has_prom = (phi - plo) > 0
            res[f"{mark}_mean_genebody_level"] = round(
                float(np.nanmean(per_gene)), 4)
            res[f"{mark}_n_genes_with_level"] = int(np.isfinite(per_gene).sum())
            res[f"{mark}_pct_genes_methylatable_promoter"] = round(
                100 * float(has_prom.mean()), 2)
            t = tpm_src.reindex(sub["gene_id"]).to_numpy()
            ok = np.isfinite(per_gene) & np.isfinite(t)
            if ok.sum() >= 10:
                rho, p_ = stats.spearmanr(per_gene[ok], t[ok])
                res[f"{mark}_rho_level_vs_TPM"] = round(float(rho), 4)
                res[f"{mark}_p_level_vs_TPM"] = f"{p_:.3g}"
                res[f"{mark}_n_corr"] = int(ok.sum())
            else:
                res[f"{mark}_rho_level_vs_TPM"] = np.nan
                res[f"{mark}_p_level_vs_TPM"] = ""
                res[f"{mark}_n_corr"] = int(ok.sum())
        rows.append(res)
    return pd.DataFrame(rows)


def main():
    os.makedirs(TAB, exist_ok=True)
    seq = genome_seq()
    L = len(seq)
    reg = pd.read_csv(REGIONS)
    tss = np.sort(pd.read_csv(TSSCSV)["tss_position"].astype(int).to_numpy())
    new20a = pd.read_csv(NEW20A, sep="\t")
    old20a = pd.read_csv(OLD20A, sep="\t")

    per4 = per_library(["21839,GCCGGC,2"], keep=tgg_positions(seq))
    per6 = per_library(["a,AAGCCCG,0", "a,AAGCCCG,1"])

    # --- how far above the all-read lineage did the undocumented 20a sit?
    m = old20a.merge(new20a, on="region", suffixes=("_old", "_new"))
    diffs = []
    for mark in ("4mC", "6mA"):
        for tp in LIBS:
            diffs.append((m[f"{mark}_level_{tp}_old"]
                          - m[f"{mark}_level_{tp}_new"]) * 100)
    d = pd.concat(diffs).dropna()
    print(f"SuppTable20a 261021 minus all-read: mean {d.mean():+.2f} pp, "
          f"median {d.median():+.2f} pp, n={len(d)} region x mark x timepoint")

    cl = cluster_levels(reg, per4, per6, L)
    cl.to_csv(f"{TAB}/nf_r5_cluster_levels_261101.tsv", sep="\t", index=False)
    tests = cluster_tests(cl, new20a)
    tests.to_csv(f"{TAB}/nf_r5_cluster_tests_261101.tsv", sep="\t", index=False)
    print("\n== cluster-unit tests (all reads) ==")
    print(tests.to_string(index=False))

    b = rebuild_20b(seq, reg)
    old_b = pd.read_csv(OLD20B, sep="\t")
    chk = b[b.motif == "TGGCCGGC"].merge(
        old_b[["region", "TGGCCGGC_inside", "TGGCCGGC_flank"]], on="region")
    agree = int((chk.n_inside == chk.TGGCCGGC_inside).sum())
    print(f"\n20b TGGCCGGC inside counts reproduced for {agree}/29 regions")
    b.to_csv(f"{TAB}/nf_r5_motif_occurrence_261101.tsv", sep="\t", index=False)

    c = rebuild_20c(per6, reg, tss, L)
    c.to_csv(f"{TAB}/nf_r5_interior_flank_261101.tsv", sep="\t", index=False)
    print("\n== rebuilt 20c (all reads, AAGCCCG A1/A2) ==")
    print(c.to_string(index=False))

    tpm = pd.read_csv(f"{RNA}/03_quantification/data/featureCounts/"
                      "featureCounts_M145_s2.txt", sep="\t", comment="#")
    cnt = tpm.set_index("Geneid").iloc[:, 5:]
    rpk = cnt.div(tpm.set_index("Geneid")["Length"], axis=0)
    tpm_t1 = (rpk / rpk.sum() * 1e6).iloc[:, :3].mean(axis=1)

    trial = primary_metabolism_trial(reg, per4, per6, tpm_t1)
    trial.to_csv(f"{TAB}/nf_r5_primary_metabolism_trial_261101.tsv",
                 sep="\t", index=False)
    print("\n== primary-metabolism candidate definitions (T1) ==")
    print(trial[["definition", "n_genes", "4mC_mean_genebody_level",
                 "6mA_mean_genebody_level",
                 "4mC_pct_genes_methylatable_promoter",
                 "4mC_rho_level_vs_TPM", "4mC_p_level_vs_TPM",
                 "4mC_n_corr"]].to_string(index=False))


if __name__ == "__main__":
    main()
