#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""NF-R5 (Figure 6 under the C04 six-main-figure layout) — recompute, 2026-11-02.

    python3 nf_r5_recompute_261102.py

Two things change against nf_r5_recompute_261101.py, and nothing else:

(1) PANEL (a) IS NEW.  The C04 reorganisation dissolves the old Figure 4 and
    folds its panel (a) — per-mark, per-timepoint-pair site retention — into
    this figure as its first panel.  The quantity is recomputed here from the
    pileup rather than copied from
    `data_261013/v72645d04_STAB_retention_by_mark_261013.csv`.

    DELIBERATE LINEAGE CHANGE, declared.  The 261013 table was built on the
    CONFIDENCE-FILTERED pileup (its own `modkit_thresholds` column records
    C:0.9375 / A:0.6386719).  This lane may not use that lineage: a floor on
    the modification state biases level-scale statistics upward (canon F7,
    +15.5 pp), and a call rule of `pooled Nmod/Nvalid >= 0.5` is a level-scale
    statistic.  Everything here is therefore all-read
    (`pileup_nofilt_261015`, Nfail == 0 asserted per library), and the
    recomputed percentages are expected to sit BELOW the published ones.  The
    difference is written to the output table, not hidden.

    Held fixed against the published panel: the four mark universes
    (TGGCCGGC 4mC C5 n = 2,093; AAGCCCG 4mC C5 n = 1,334; AAGCCCG 6mA A1
    n = 1,334; HAAGCCCG 6mA A1 n = 565 — all four reproduced from the
    reference sequence, not read in), the three timepoint pairs, the library
    grouping, the depth floor (pooled Nvalid >= 30 in BOTH timepoints) and the
    call rule (pooled Nmod/Nvalid >= 0.5).

    DENOMINATOR, corrected.  `panel_F4a.py`'s docstring states that retention
    is "sites called modified in BOTH timepoints ... divided by the union".
    That is not the quantity its input table carries.  Re-running this code on
    the filtered lineage the published table was built from reproduces its
    `n_measurable_both`, `n_modified_first`, `n_retained` and `n_union`
    EXACTLY in all 12 rows, and the published `retention_pct` then equals
    n_retained / n_modified_FIRST (1834/1886 = 97.24, 283/323 = 87.62, ...),
    while n_retained / n_union is the table's separate `jaccard` column
    (1834/1944 = 0.9434).  The manuscript legend agrees with the first form
    ("fraction of sites modified at one timepoint that are still modified at
    the next"), so that is what is computed and drawn here; the docstring in
    panel_F4a.py is the error.  Both ratios are written to the output table.

(2) PANEL (d) REPLACES the provisional three-candidate panel (c).  Author
    decision D155 fixes the definition: primary metabolism is the core gene
    set of the KEGG metabolic pathways.  Operational definition used here,
    stated in full because KEGG offers several and the choice is ours:

        A gene of S. coelicolor A3(2) (KEGG organism code `sco`) is a PRIMARY
        METABOLISM CORE gene if BRITE hierarchy br:sco00001 places it under
        A09100 Metabolism in at least one B-level category OTHER than
          09109 Metabolism of terpenoids and polyketides,
          09110 Biosynthesis of other secondary metabolites,
          09111 Xenobiotics biodegradation and metabolism,
          09112 Not included in regular maps.
        09109 and 09110 are secondary metabolism and are the categories D155
        means to exclude; 09111 is catabolism of foreign compounds, not core;
        09112 carries no pathway map and is a residue category.  br:sco00001
        contains no global/overview maps (sco01100, sco01110, ...), so no
        aggregate map had to be removed.  A gene in BOTH a kept and an
        excluded category is KEPT — the exclusion removes genes whose ONLY
        metabolic assignment is secondary.

    The antiSMASH BGC gene set is retained as the comparator, and all
    remaining protein-coding genes are drawn as the background reference so
    that neither group is read against nothing.

Outputs (tables/):
    nf_r5_retention_261102.tsv            panel (a), all-read, 4 marks x 3 pairs
    nf_r5_kegg_primary_metabolism_261102.tsv   the KEGG gene set itself
    nf_r5_genesets_levels_261102.tsv      panel (d) per-gene levels, long form
    nf_r5_geneset_summary_261102.tsv      panel (d) group summary + tests
Panels (b) and (c) reuse the 261101 tables unchanged (same all-read lineage).
"""
import json
import os
import re
import sys
import urllib.request
from datetime import date

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(HERE, "tables")
RNA = "/Users/okaban/bioinfo/rna-seq"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"
NOFILT = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015"
ANNOT = (f"{RNA}/05_annotation/analysis/05_annotation_260128_v1/tables/"
         "gene_annotation_basic.tsv")
REGIONS = f"{RNA}/15_paper_figures/scripts/supp_261014/data/bgc_regions_29.csv"
OLD_RET = (f"{RNA}/15_paper_figures/data_261013/"
           "v72645d04_STAB_retention_by_mark_261013.csv")
KEGG_CACHE = os.path.join(TAB, "kegg_br_sco00001_261102.keg")
KEGG_INFO_CACHE = os.path.join(TAB, "kegg_release_261102.txt")

LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
        "T3": ["3-2", "3-3", "3-4"]}
FLAT = [l for v in LIBS.values() for l in v]
PAIRS = [("T1", "T2"), ("T2", "T3"), ("T1", "T3")]
BEDCOL = ["chrom", "start", "end", "name", "score", "strand", "ts", "te",
          "color", "Nvalid", "pct", "Nmod", "Ncanon", "Nother", "Ndel",
          "Nfail", "Ndiff", "Nnocall"]
MIN_DEPTH = 30
CALL = 0.5

EXCLUDE_B = {"09109", "09110", "09111", "09112"}
KEEP_A = "09100"


# ------------------------------------------------------------------ sequence
def genome_seq():
    return "".join(l.strip() for l in open(GENOME)
                   if not l.startswith(">")).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def universe(seq, pattern, offset):
    """{(0-based genome index, strand)} of the base `offset` into `pattern`.

    `pattern` is a regex over the plus strand; the minus-strand hits are found
    on the reverse complement and mapped back to plus-strand coordinates, which
    is the convention modkit's bed uses (the row's `start` is the plus-strand
    index and `strand` says which strand carries the base).
    """
    L = len(seq)
    out = {(m.start() + offset, "+") for m in re.finditer("(?=" + pattern + ")",
                                                          seq)}
    out |= {(L - 1 - (m.start() + offset), "-")
            for m in re.finditer("(?=" + pattern + ")", rc(seq))}
    return out


# ------------------------------------------------------------------ pileup
def load_bed(lib, names):
    d = pd.read_csv(f"{NOFILT}/{lib}_nofilt.bed", sep="\t", header=None,
                    names=BEDCOL, usecols=["start", "name", "strand",
                                           "Nvalid", "Nmod", "Nfail"])
    assert int(d["Nfail"].sum()) == 0, f"{lib}: not an all-read pileup"
    d = d[d["name"].isin(names)]
    return d.drop_duplicates(subset=["start", "strand"])


def pooled_counts(beds, libs, uni):
    """Sum Nvalid/Nmod over `libs`, restricted to the (start, strand) universe."""
    acc = None
    for l in libs:
        d = beds[l]
        k = d.set_index(["start", "strand"])[["Nvalid", "Nmod"]]
        acc = k if acc is None else acc.add(k, fill_value=0)
    idx = pd.MultiIndex.from_tuples(sorted(uni), names=["start", "strand"])
    return acc.reindex(idx).fillna(0.0)


def retention(seq):
    L = len(seq)
    marks = [
        ("TGGCCGGC 4mC (position 5)", "21839,GCCGGC,2",
         universe(seq, "TGGCCGGC", 4)),
        ("AAGCCCG 4mC (position 5)", "21839,AAGCCCG,4",
         universe(seq, "AAGCCCG", 4)),
        ("AAGCCCG 6mA (position 1, all instances)", "a,AAGCCCG,0",
         universe(seq, "AAGCCCG", 0)),
        ("HAAGCCCG 6mA (position 1, H context)", "a,AAGCCCG,0",
         universe(seq, "[ACT]AAGCCCG", 1)),
    ]
    names = {m[1] for m in marks}
    beds = {l: load_bed(l, names) for l in FLAT}
    old = pd.read_csv(OLD_RET)

    rows = []
    for label, bedname, uni in marks:
        sub = {l: beds[l][beds[l]["name"] == bedname] for l in FLAT}
        pool = {tp: pooled_counts(sub, libs, uni) for tp, libs in LIBS.items()}
        for a, b in PAIRS:
            pa, pb = pool[a], pool[b]
            meas = (pa["Nvalid"] >= MIN_DEPTH) & (pb["Nvalid"] >= MIN_DEPTH)
            la = np.where(pa["Nvalid"] > 0, pa["Nmod"] / pa["Nvalid"].clip(lower=1), 0.0)
            lb = np.where(pb["Nvalid"] > 0, pb["Nmod"] / pb["Nvalid"].clip(lower=1), 0.0)
            ma = meas.to_numpy() & (la >= CALL)
            mb = meas.to_numpy() & (lb >= CALL)
            both = int((ma & mb).sum())
            union = int((ma | mb).sum())
            o = old[(old["mark"] == label) & (old["pair"] == f"{a}->{b}")]
            pub = float(o["retention_pct"].iloc[0]) if len(o) else np.nan
            first = int(ma.sum())
            pct = 100.0 * both / first if first else np.nan
            rows.append(dict(
                mark=label, motif_universe_n=len(uni), pair=f"{a}->{b}",
                n_measurable_both=int(meas.sum()),
                n_modified_first=first, n_modified_second=int(mb.sum()),
                n_retained=both, n_union=union,
                retention_pct=round(pct, 2),
                jaccard=round(both / union, 4) if union else np.nan,
                published_261013_pct=pub,
                delta_vs_published_pp=round(pct - pub, 2),
                lineage="all-read (modkit --no-filtering, pileup_nofilt_261015)",
                rule=("per-timepoint pooling of 3 libraries; pooled Nvalid >= 30 "
                      "in BOTH timepoints; site called modified if pooled "
                      "Nmod/Nvalid >= 0.5; retention = n_retained / "
                      "n_modified_first (jaccard column = n_retained / n_union)")))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ KEGG
def kegg_fetch(url, cache):
    if os.path.exists(cache):
        return open(cache, encoding="utf-8").read()
    with urllib.request.urlopen(url, timeout=120) as r:
        txt = r.read().decode("utf-8")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    open(cache, "w", encoding="utf-8").write(txt)
    return txt


def kegg_primary_metabolism():
    """Parse br:sco00001 into the D155 primary-metabolism core gene set."""
    keg = kegg_fetch("https://rest.kegg.jp/get/br:sco00001", KEGG_CACHE)
    info = kegg_fetch("https://rest.kegg.jp/info/kegg", KEGG_INFO_CACHE)
    rel = {}
    for ln in info.splitlines():
        f = ln.strip().split()
        if len(f) >= 3 and f[0] in ("pathway", "brite", "genes"):
            rel[f[0]] = f"{f[1]} ({f[2]})"
    release = ("brite " + rel.get("brite", "?") + "; pathway "
               + rel.get("pathway", "?"))

    a = b = c = None
    hits = {}
    for line in keg.splitlines():
        if not line or line[0] not in "ABCD":
            continue
        lvl, rest = line[0], line[1:].strip()
        if lvl == "A":
            a = rest.split(None, 1)[0]
        elif lvl == "B":
            b = rest.split(None, 1)[0] if rest else None
        elif lvl == "C":
            c = rest.split(None, 1)[0]
        elif lvl == "D" and a == KEEP_A and b:
            sco = rest.split(None, 1)[0]
            if not sco.startswith("SCO"):
                continue
            rec = hits.setdefault(sco, {"keep": set(), "drop": set(),
                                        "maps": set(), "ko": ""})
            (rec["drop"] if b in EXCLUDE_B else rec["keep"]).add(b)
            if b not in EXCLUDE_B:
                rec["maps"].add(c)
            if "\t" in line:
                ko = line.split("\t")[1].strip().split(None, 1)[0]
                if ko.startswith("K"):
                    rec["ko"] = ko
    rows = []
    for sco, rec in sorted(hits.items()):
        rows.append(dict(
            sco_id=sco, ko=rec["ko"],
            kept_B_categories=";".join(sorted(rec["keep"])),
            excluded_B_categories=";".join(sorted(rec["drop"])),
            n_pathway_maps=len(rec["maps"]),
            pathway_maps=";".join(sorted(rec["maps"])),
            in_primary_core=bool(rec["keep"])))
    return pd.DataFrame(rows), release


# ------------------------------------------------------------------ panel (d)
def gene_levels(seq, genes, mark_name, bedname, uni_pattern, uni_offset):
    uni = universe(seq, uni_pattern, uni_offset)
    beds = {l: load_bed(l, {bedname}) for l in LIBS["T1"]}
    beds = {l: d[d["name"] == bedname] for l, d in beds.items()}
    p = pooled_counts(beds, LIBS["T1"], uni).reset_index()
    p = p[p["Nvalid"] >= MIN_DEPTH]
    pos = np.sort(p["start"].to_numpy() + 1)          # 1-based
    lvl = p.sort_values("start")["Nmod"].to_numpy() / \
        p.sort_values("start")["Nvalid"].to_numpy()
    s = genes["start"].to_numpy()
    e = genes["end"].to_numpy()
    lo = np.searchsorted(pos, s)
    hi = np.searchsorted(pos, e, side="right")
    out = np.array([lvl[x:y].mean() if y > x else np.nan
                    for x, y in zip(lo, hi)])
    n_pos = (hi - lo)
    return out, n_pos


def main():
    os.makedirs(TAB, exist_ok=True)
    seq = genome_seq()

    # ---- panel (a)
    ret = retention(seq)
    ret.to_csv(f"{TAB}/nf_r5_retention_261102.tsv", sep="\t", index=False)
    print("== panel (a) site retention, all-read ==")
    print(ret[["mark", "motif_universe_n", "pair", "n_retained", "n_union",
               "retention_pct", "published_261013_pct",
               "delta_vs_published_pp"]].to_string(index=False))
    for m in ret["mark"].unique():
        d = ret[ret["mark"] == m]
        print(f"  range {m}: {d.retention_pct.min():.2f}-{d.retention_pct.max():.2f}% "
              f"(published {d.published_261013_pct.min():.2f}-"
              f"{d.published_261013_pct.max():.2f}%)")

    # ---- KEGG gene set
    kegg, release = kegg_primary_metabolism()
    ann = pd.read_csv(ANNOT, sep="\t")
    pc = ann[ann["gene_biotype"] == "protein_coding"].copy()
    pc["sco_id"] = pc["old_locus_tag"].astype(str)
    kegg = kegg.merge(pc[["sco_id", "gene_id", "start", "end", "strand",
                          "product"]], on="sco_id", how="left")
    kegg["mapped"] = kegg["gene_id"].notna()
    kegg["kegg_release"] = release
    kegg["retrieved_on"] = date.today().isoformat()
    kegg["definition"] = ("BRITE br:sco00001, A09100 Metabolism, B-categories "
                          "other than 09109/09110/09111/09112")
    kegg.to_csv(f"{TAB}/nf_r5_kegg_primary_metabolism_261102.tsv", sep="\t",
                index=False)
    core = kegg[kegg["in_primary_core"] & kegg["mapped"]]
    print(f"\n== KEGG ({release}) ==")
    print(f"  genes under A09100 Metabolism: {len(kegg)}")
    print(f"  primary-metabolism core (D155 definition): "
          f"{int(kegg['in_primary_core'].sum())}")
    print(f"  of which mapped to a locus_tag: {len(core)}")
    print(f"  dropped (only secondary/xenobiotic/no-map): "
          f"{int((~kegg['in_primary_core']).sum())}")

    # ---- panel (d) groups
    reg = pd.read_csv(REGIONS)
    mid = ((pc["start"] + pc["end"]) / 2).to_numpy()
    in_bgc = np.zeros(len(pc), bool)
    for _, r in reg.iterrows():
        in_bgc |= (mid >= r["start"]) & (mid <= r["end"])
    pc = pc.assign(in_bgc=in_bgc)
    core_ids = set(core["gene_id"])
    grp = np.where(pc["gene_id"].isin(core_ids), "KEGG primary metabolism",
                   np.where(pc["in_bgc"], "antiSMASH BGC genes",
                            "Other protein-coding genes"))
    # a gene that is both KEGG-core and inside a BGC region is counted as
    # KEGG-core; report the overlap rather than letting it pass silently
    overlap = int((pc["gene_id"].isin(core_ids) & pc["in_bgc"]).sum())
    pc = pc.assign(group=grp)

    lv4, n4 = gene_levels(seq, pc, "4mC", "21839,GCCGGC,2", "TGGCCGGC", 4)
    lv6, n6 = gene_levels(seq, pc, "6mA", "a,AAGCCCG,0", "AAGCCCG", 0)
    long = pd.concat([
        pd.DataFrame(dict(gene_id=pc["gene_id"], sco_id=pc["sco_id"],
                          group=pc["group"], mark="4mC",
                          n_positions=n4, level_T1=lv4)),
        pd.DataFrame(dict(gene_id=pc["gene_id"], sco_id=pc["sco_id"],
                          group=pc["group"], mark="6mA",
                          n_positions=n6, level_T1=lv6))])
    long.to_csv(f"{TAB}/nf_r5_genesets_levels_261102.tsv", sep="\t",
                index=False)

    rows = []
    for mark in ("4mC", "6mA"):
        d = long[(long["mark"] == mark) & long["level_T1"].notna()]
        gp = {g: d[d["group"] == g]["level_T1"].to_numpy()
              for g in d["group"].unique()}
        for g, v in gp.items():
            rows.append(dict(mark=mark, group=g,
                             n_genes_in_group=int((pc["group"] == g).sum()),
                             n_genes_with_level=len(v),
                             mean_level=round(float(v.mean()), 4),
                             median_level=round(float(np.median(v)), 4),
                             q1=round(float(np.percentile(v, 25)), 4),
                             q3=round(float(np.percentile(v, 75)), 4)))
        a = gp.get("KEGG primary metabolism")
        b = gp.get("antiSMASH BGC genes")
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(mark=mark, group="TEST KEGG vs BGC (Mann-Whitney U)",
                         n_genes_in_group=len(a), n_genes_with_level=len(b),
                         mean_level=round(float(u.statistic), 1),
                         median_level=float(f"{u.pvalue:.3g}"),
                         q1=np.nan, q3=np.nan))
    summ = pd.DataFrame(rows)
    summ.attrs["overlap"] = overlap
    summ.to_csv(f"{TAB}/nf_r5_geneset_summary_261102.tsv", sep="\t",
                index=False)
    print(f"\n== panel (d) groups (T1, all-read); KEGG-core genes also inside "
          f"an antiSMASH region: {overlap} (assigned to KEGG-core) ==")
    print(summ.to_string(index=False))

    meta = dict(kegg_release=release, retrieved_on=date.today().isoformat(),
                organism="sco (Streptomyces coelicolor A3(2))",
                brite="br:sco00001",
                excluded_B_categories=sorted(EXCLUDE_B),
                n_genes_A09100=int(len(kegg)),
                n_primary_core=int(kegg["in_primary_core"].sum()),
                n_primary_core_mapped=int(len(core)),
                n_bgc_genes=int((pc["group"] == "antiSMASH BGC genes").sum()),
                n_other=int((pc["group"] == "Other protein-coding genes").sum()),
                kegg_core_inside_bgc=overlap)
    json.dump(meta, open(f"{TAB}/nf_r5_kegg_meta_261102.json", "w"), indent=1)
    print("\n" + json.dumps(meta))


if __name__ == "__main__":
    main()
