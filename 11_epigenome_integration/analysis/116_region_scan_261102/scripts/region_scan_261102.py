#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Region-by-region scan of the three methylation series (author request of
2026-11-02: "check every region, not only the promoter").

    python3 region_scan_261102.py

Four parts, four tables, one markdown judgement sheet written by the caller.

Part 1  transcription-unit 3' end (terminator side) - the control for the only
        surviving spatial finding of the manuscript, the 6mA level depletion
        near promoters (PROM.6mA.RATIO = 0.9105, CI 0.8955-0.9257,
        p = 1.12e-06). Same window (+/-100 bp), same read set (all-read
        pileup_nofilt_261015, Nfail == 0), same local-GC quintile matching,
        same unit of replication (culture, n = 9), same estimator - only the
        anchor set changes. The Jeong TSS rung is recomputed in the same loop
        and asserted against 0.9105 so the machinery is proved identical.

Part 2  intergenic regions by the orientation of the flanking gene pair
        (tandem / divergent / convergent). Divergent spacers carry two
        promoters; if the promoter effect were dose-dependent the depletion
        should double there. Region length differs systematically by
        orientation, so every statistic is reported both GC-standardised and
        GC x length-standardised.

Part 3  mobile elements (transposase / integrase / IS / recombinase), gene
        body and gene body +/- 2 kb. H12 (methylation controlling mobile
        element transposition) is the one hypothesis never tested on any of
        the three series.

Part 4  motif placement across promoter / 5'UTR / CDS / other against a
        local-GC-matched expectation. Rebuilds the lost canonical table
        motif_region_GCmatched_261004.csv, which is not on disk, and
        explicitly collates with the numbers the manuscript notes recorded
        from it.

Inputs (canonical only)
-----------------------
level scale   pileup_nofilt_261015 via recompute_TSS_set_controls_261021.py
              (all reads, Nfail == 0) - this is the read set on which 0.9105
              was computed, hence the read set any control for it must use.
census scale  pileup_cache_261015/*.npz (9 libraries): per timepoint pool the
              3 libraries, require pooled TRUE depth Nvalid + Nfail >= 30, call
              a site at pooled Nmod / Nvalid >= 0.5. The floor is on true
              depth, never on the number of valid calls.
TSS           jeong2016_all_tss.csv, 3,570 experimental dRNA-seq TSSs
TU            G10_gene_to_transcription_unit_261008.csv (5,440 units /
              7,593 genes; the published count is 5,432, a known 0.15 % gap)
genes         05_annotation .../gene_annotation_basic.tsv
Pfam          113_mtase_completeness_261024/data/pfam_full_chrom.domtbl

Never read: 01_integration/high_confidence_sites_weighted.csv,
70_methylation_dynamics/, 90_per_timepoint_census_audit/, and the sequence
column of 07_motif_analysis/methylation_site_sequences.csv - all of them put a
floor on the methylation state itself.

Series (canon): TGGCCGGC C5 4mC (offset 4) / AAGCCCG C5 4mC (offset 4) /
AAGCCCG A1,A2 6mA (offsets 0,1).
"""
import importlib.util
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "tables")
ROOT = "/Users/okaban/bioinfo/rna-seq"
MOD = (f"{ROOT}/11_epigenome_integration/analysis/48_TSS_methylation_gradient/"
       "scripts/recompute_TSS_set_controls_261021.py")
NFR3 = (f"{ROOT}/15_paper_figures/scripts/newfig_261101/"
        "nf_r3_recompute_261101.py")
CACHE = (f"{ROOT}/15_paper_figures/scripts/recompute_261015/"
         "pileup_cache_261015")
TSSCSV = (f"{ROOT}/11_epigenome_integration/analysis/18_tss_analyses/"
          "jeong2016_all_tss.csv")
TUCSV = (f"{ROOT}/11_epigenome_integration/analysis/107_G07_G10_audit_261008/"
         "tables/G10_gene_to_transcription_unit_261008.csv")
ANNO = (f"{ROOT}/05_annotation/analysis/05_annotation_260128_v1/tables/"
        "gene_annotation_basic.tsv")
PFAM = (f"{ROOT}/11_epigenome_integration/analysis/113_mtase_completeness_261024/"
        "data/pfam_full_chrom.domtbl")
GFF = (f"{ROOT}/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/"
       "data/GCF_000203835.1_ASM20383v1_genomic.gff")
CHROM = "NC_003888.3"
WIN = 100          # the window of PROM.6mA.RATIO
NQ = 5             # local-GC quintiles, as in the 0.9105 computation
TIMEPOINTS = {"T1": ("1-1", "1-2", "1-3"), "T2": ("2-1", "2-3", "2-4"),
              "T3": ("3-2", "3-3", "3-4")}
CODE = {"6mA": 0, "4mC": 1}
PROM_RATIO_CANON = 0.9105


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load("tssmod", MOD)        # cum_gc/local_gc/min_distance/standardised_ratio
R3 = _load("nfr3", NFR3)        # motif_positions/pooled/classify


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------- estimators
def strata_ratio(level, near, strata):
    """Direct standardisation of the far arm to the near arm's stratum mix.

    Identical arithmetic to M.standardised_ratio, but the strata are supplied
    instead of being qcut from a single covariate, so GC x length composite
    strata can be used. Asserted against M.standardised_ratio in main().
    """
    out = []
    for li in range(level.shape[0]):
        num = den = 0.0
        for s in np.unique(strata):
            m = strata == s
            nn, dd = m & near, m & ~near
            if nn.sum() == 0 or dd.sum() == 0:
                continue
            w = nn.sum()
            num += w * np.nanmean(level[li, nn])
            den += w * np.nanmean(level[li, dd])
        out.append(num / den if den else np.nan)
    return np.array(out)


def gc_matched_expectation(motif_pos0, labels, gc_bin, classes):
    """Expected motif count per class, local-GC-standardised.

    exp(c) = SUM_bins n_motif(bin) * n_positions(c, bin) / n_positions(bin)
    Positions are chromosome bases; `labels` and `gc_bin` are per-base arrays.
    """
    nb = gc_bin.max() + 1
    mb = np.bincount(gc_bin[motif_pos0], minlength=nb).astype(float)
    tot = np.bincount(gc_bin, minlength=nb).astype(float)
    exp = {}
    for c in classes:
        cb = np.bincount(gc_bin[labels == c], minlength=nb).astype(float)
        exp[c] = float(np.sum(mb * cb / np.maximum(tot, 1)))
    return exp


def summarise_plus(r):
    """M.summarise plus the 80 %-power floor of the culture-level test.

    mde = exp((t_0.975,df + t_0.80,df) * SE(log ratio)) - 1, the smallest
    deviation from 1.0 this n = 9 design would detect 80 % of the time. Quoted
    beside every null so the null is readable as a bound, not as absence.
    """
    ratio, lo, hi, p, n = M.summarise(r)
    lg = np.log(r)
    se = lg.std(ddof=1) / np.sqrt(n)
    mde = float(np.exp((stats.t.ppf(0.975, n - 1) + stats.t.ppf(0.80, n - 1)) * se))
    return ratio, lo, hi, p, n, float(lg.std(ddof=1)), mde


def binom_mde(exp, n):
    """Smallest detectable |O/E - 1| at 80 % power for the count test."""
    p0 = min(max(exp / n, 1e-12), 1 - 1e-12)
    z = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)
    return float(z * np.sqrt(p0 * (1 - p0) / n) / p0)


def binom_two_sided(obs, exp, n):
    p = min(max(exp / n, 1e-12), 1 - 1e-12)
    return float(stats.binomtest(int(obs), int(n), p).pvalue)


# ------------------------------------------------------------------ universes
def build_universes(genome, L):
    """The three canonical series on the all-read pileup."""
    tgg = set(int(m.start() + 4) for m in re.finditer("(?=TGGCCGGC)", genome))
    tgg |= set(L - 1 - (m.start() + 4)
               for m in re.finditer("(?=TGGCCGGC)", rc(genome)))
    uni = {}
    uni["6mA, AAGCCCG A1/A2"] = M.load_universe(["a,AAGCCCG,0", "a,AAGCCCG,1"])
    raw = M.load_universe(["21839,GCCGGC,2"])
    k = np.fromiter((p in tgg for p in raw[0]), bool, len(raw[0]))
    uni["4mC, TGGCCGGC C5"] = (raw[0][k], raw[1][k], raw[2][:, k], raw[3][:, k])
    uni["4mC, AAGCCCG C5"] = M.load_universe(["21839,AAGCCCG,4"])
    return uni


SERIES = [("TGGCCGGC", "4mC", [4], "4mC, TGGCCGGC C5"),
          ("AAGCCCG", "4mC", [4], "4mC, AAGCCCG C5"),
          ("AAGCCCG", "6mA", [0, 1], "6mA, AAGCCCG A1/A2")]


def census_store():
    """Pooled per-timepoint census arrays, canonical rule (true-depth floor)."""
    return {(tp, mk): R3.pooled(libs, code)
            for tp, libs in TIMEPOINTS.items()
            for mk, code in CODE.items()}


def census_fraction(store, genome, mask_by_name):
    """Modified-site fraction per region mask, per timepoint x series."""
    rows = []
    for tp in TIMEPOINTS:
        for motif, mark, offs, label in SERIES:
            pos0, plus = R3.motif_positions(genome, motif, offs)
            s = store[(tp, mark)]
            idx = pd.MultiIndex.from_arrays([pos0, plus], names=["start", "plus"])
            d = s.reindex(idx)
            keep = (d["dep"] >= 30).to_numpy(dtype=bool)
            p = pos0[keep]
            called = d["lvl"].to_numpy()[keep] >= 0.5
            for name, mask in mask_by_name.items():
                inr = mask[p]
                rows.append(dict(
                    timepoint=tp, series=label, region=name,
                    n_measurable=int(inr.sum()),
                    n_modified=int((inr & called).sum()),
                    pct_modified=(round(100 * (inr & called).sum() / inr.sum(), 2)
                                  if inr.sum() else np.nan),
                    n_measurable_all=int(len(p)),
                    pct_modified_all=round(100 * called.sum() / max(len(p), 1), 2)))
    return pd.DataFrame(rows)


# -------------------------------------------------------------------- part 1
def tu_termini(tu):
    g = tu[tu.contig == CHROM].groupby("TU")
    assert int((g.strand.nunique() > 1).sum()) == 0, "a TU spans both strands"
    a = g.agg(strand=("strand", "first"), lo=("start", "min"),
              hi=("end", "max"), ngene=("gene_id", "size")).reset_index()
    a["end3"] = np.where(a.strand == "+", a.hi, a.lo)
    a["end5"] = np.where(a.strand == "+", a.lo, a.hi)
    return a


def signed_end_distance(pos1, ends, strands):
    """Signed distance to the nearest terminus, along that unit's own strand.

    Negative = before the terminus (inside the unit for a 3' end), positive =
    beyond it. Mirrors nf_r3's signed_tss_distance.
    """
    out = np.full(len(pos1), np.iinfo(np.int64).max // 4, dtype=np.int64)
    for s, sign in (("+", 1), ("-", -1)):
        t = np.sort(ends[strands == s].astype(np.int64))
        i = np.searchsorted(t, pos1)
        cand = []
        for j in (i - 1, i):
            cand.append(sign * (pos1 - t[np.clip(j, 0, len(t) - 1)]))
        cand = np.vstack(cand)
        best = cand[np.argmin(np.abs(cand), axis=0), np.arange(len(pos1))]
        take = np.abs(best) < np.abs(out)
        out[take] = best[take]
    return out


def part1(uni, agg, tss_all, cGC, cN, L):
    rows = []
    for label, (pos0, st0, Nv, Nm) in uni.items():
        level = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
        gc = M.local_gc(pos0, cGC, cN, L, WIN)
        dT = M.min_distance(pos0 + 1, tss_all)
        d3 = signed_end_distance(pos0 + 1, agg.end3.to_numpy(), agg.strand.to_numpy())
        d5 = signed_end_distance(pos0 + 1, agg.end5.to_numpy(), agg.strand.to_numpy())
        anchors = [
            ("TSS_jeong_reference", dT <= WIN,
             "the rung that defines PROM.6mA.RATIO"),
            ("TU_3prime", np.abs(d3) <= WIN, "terminator side, all positions"),
            ("TU_5prime", np.abs(d5) <= WIN, "initiation side, annotation-derived"),
            ("TU_3prime_minus_TSS", (np.abs(d3) <= WIN) & (dT > WIN),
             "3' end, positions within 100 bp of a TSS removed"),
            ("TU_3prime_minus_TU5", (np.abs(d3) <= WIN) & (np.abs(d5) > WIN),
             "3' end, positions within 100 bp of a TU 5' end removed"),
            ("TU_3prime_isolated", (np.abs(d3) <= WIN) & (dT > WIN) & (np.abs(d5) > WIN),
             "3' end, both leaks removed"),
            ("TU_3prime_inside", (d3 <= 0) & (d3 >= -WIN),
             "the 100 bp before the terminus only"),
            ("TU_3prime_beyond", (d3 > 0) & (d3 <= WIN),
             "the 100 bp after the terminus only"),
        ]
        for name, near, note in anchors:
            if near.sum() < 20 or (~near).sum() < 20:
                rows.append(dict(series=label, anchor=name,
                                 note=note + " [NOT TESTED: n_near < 20]",
                                 n_near=int(near.sum()), ratio=np.nan))
                continue
            r, lo, hi, p, n, sd, mde = summarise_plus(
                M.standardised_ratio(level, Nv, Nm, near, gc, nq=NQ))
            rows.append(dict(
                series=label, anchor=name, note=note, window_bp=WIN,
                sd_log_ratio=round(sd, 5), mde80_ratio=round(mde, 4),
                read_set="all reads (pileup_nofilt_261015)",
                gc_matching=f"local GC +/-{WIN} bp, {NQ} quantiles",
                universe_n=int(len(pos0)), n_near=int(near.sum()),
                n_distal=int((~near).sum()), n_cultures=n,
                mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                ratio=round(r, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4),
                p_value=f"{p:.3g}"))
    out = pd.DataFrame(rows)
    ref = out[(out.series == "6mA, AAGCCCG A1/A2")
              & (out.anchor == "TSS_jeong_reference")].iloc[0]
    assert round(float(ref.ratio), 4) == PROM_RATIO_CANON, (
        f"machinery drift: TSS rung {ref.ratio} != {PROM_RATIO_CANON}")
    return out


# -------------------------------------------------------------------- part 2
def intergenic_regions(anno):
    g = (anno[anno.contig == CHROM]
         .sort_values(["start", "end"]).reset_index(drop=True))
    a, b = g.iloc[:-1].reset_index(drop=True), g.iloc[1:].reset_index(drop=True)
    d = pd.DataFrame({
        "left_gene": a.gene_id, "right_gene": b.gene_id,
        "left_strand": a.strand, "right_strand": b.strand,
        "ig_start": a.end + 1, "ig_end": b.start - 1})
    d["length"] = d.ig_end - d.ig_start + 1
    d = d[d.length >= 1].reset_index(drop=True)
    cls = np.where(d.left_strand == d.right_strand, "tandem",
                   np.where((d.left_strand == "-") & (d.right_strand == "+"),
                            "divergent", "convergent"))
    d["orientation"] = cls
    return d


def labels_from_intervals(L, frame, value, base="none"):
    lab = np.full(L, base, dtype=object)
    for s, e, v in zip(frame.iloc[:, 0], frame.iloc[:, 1], value):
        lab[int(s) - 1:int(e)] = v
    return lab


def part2(uni, ig, genome, cGC, cN, L, gc_bin, store):
    classes = ["tandem", "divergent", "convergent"]
    lab = labels_from_intervals(L, ig[["ig_start", "ig_end"]], ig.orientation)
    # length stratum of the containing spacer, tertiles over all spacers
    qt = pd.qcut(ig.length, 3, labels=["short", "mid", "long"])
    labl = labels_from_intervals(L, ig[["ig_start", "ig_end"]],
                                 [f"{o}|{q}" for o, q in zip(ig.orientation, qt)])
    rows = []
    # (a) motif placement, GC-matched and GC x length-matched
    for motif, mark, offs, label in SERIES:
        pos0, plus = R3.motif_positions(genome, motif, offs)
        exp = gc_matched_expectation(pos0, lab, gc_bin, classes)
        n = len(pos0)
        for c in classes:
            obs = int((lab[pos0] == c).sum())
            rows.append(dict(
                part="2a_motif_placement", series=label, orientation=c,
                length_stratum="all", n_regions=int((ig.orientation == c).sum()),
                region_bp=int(ig.loc[ig.orientation == c, "length"].sum()),
                median_region_bp=int(ig.loc[ig.orientation == c, "length"].median()),
                motif_n_total=n, observed=obs, expected_gcmatched=round(exp[c], 2),
                O_over_E=round(obs / exp[c], 3) if exp[c] else np.nan,
                pct_observed=round(100 * obs / n, 3),
                pct_expected=round(100 * exp[c] / n, 3),
                p_value=f"{binom_two_sided(obs, exp[c], n):.3g}",
                mde80_abs_dev_from_1=round(binom_mde(exp[c], n), 3)))
        sub = [f"{c}|{q}" for c in classes for q in ("short", "mid", "long")]
        expl = gc_matched_expectation(pos0, labl, gc_bin, sub)
        for s in sub:
            c, q = s.split("|")
            obs = int((labl[pos0] == s).sum())
            rows.append(dict(
                part="2a_motif_placement", series=label, orientation=c,
                length_stratum=q,
                n_regions=int(((ig.orientation == c) & (qt == q)).sum()),
                region_bp=int(ig.loc[(ig.orientation == c) & (qt == q), "length"].sum()),
                motif_n_total=len(pos0), observed=obs,
                expected_gcmatched=round(expl[s], 2),
                O_over_E=round(obs / expl[s], 3) if expl[s] else np.nan,
                p_value=f"{binom_two_sided(obs, expl[s], len(pos0)):.3g}"))
    # (b) level ratio, culture as replicate unit
    for label, (pos0, st0, Nv, Nm) in uni.items():
        level = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
        gc = M.local_gc(pos0, cGC, cN, L, WIN)
        cl = lab[pos0]
        gq = pd.qcut(gc, NQ, labels=False, duplicates="drop")
        for c in classes + ["any_intergenic"]:
            near = (cl != "none") if c == "any_intergenic" else (cl == c)
            if near.sum() < 20:
                rows.append(dict(part="2b_level_ratio", series=label,
                                 orientation=c, length_stratum="all",
                                 contrast="NOT TESTED: n_near < 20",
                                 n_near=int(near.sum()), ratio=np.nan))
                continue
            r, lo, hi, p, n, sd, mde = summarise_plus(
                M.standardised_ratio(level, Nv, Nm, near, gc, nq=NQ))
            rows.append(dict(
                part="2b_level_ratio", series=label, orientation=c,
                length_stratum="all", contrast="class vs rest of universe",
                sd_log_ratio=round(sd, 5), mde80_ratio=round(mde, 4),
                universe_n=int(len(pos0)), n_near=int(near.sum()),
                n_distal=int((~near).sum()), n_cultures=n,
                mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                ratio=round(r, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4),
                p_value=f"{p:.3g}"))
        # divergent vs tandem, inside the intergenic universe only, GC x length
        keep = np.isin(cl, ["divergent", "tandem"])
        if keep.sum() >= 40:
            lv = level[:, keep]
            near = cl[keep] == "divergent"
            lq = pd.Series(labl[pos0][keep]).str.split("|").str[1].to_numpy()
            comp = np.array([f"{a}|{b}" for a, b in zip(gq[keep], lq)])
            for tag, strata in (("GC quintile", gq[keep].astype(str)),
                                ("GC quintile x length tertile", comp)):
                rr = strata_ratio(lv, near, strata)
                if np.isnan(rr).any():
                    continue
                r, lo, hi, p, n, sd, mde = summarise_plus(rr)
                rows.append(dict(
                    part="2b_level_ratio", series=label, orientation="divergent",
                    length_stratum="all",
                    sd_log_ratio=round(sd, 5), mde80_ratio=round(mde, 4),
                    contrast=f"divergent vs tandem spacers, {tag} matched",
                    universe_n=int(keep.sum()), n_near=int(near.sum()),
                    n_distal=int((~near).sum()), n_cultures=n,
                    mean_level_near=round(float(np.nanmean(lv[:, near])), 6),
                    mean_level_distal=round(float(np.nanmean(lv[:, ~near])), 6),
                    ratio=round(r, 4), ci95_lo=round(lo, 4),
                    ci95_hi=round(hi, 4), p_value=f"{p:.3g}"))
    # (c) census modified-site fraction
    masks = {c: (lab == c) for c in classes}
    masks["any_intergenic"] = lab != "none"
    cen = census_fraction(store, genome, masks)
    cen.insert(0, "part", "2c_census_fraction")
    cen = cen.rename(columns={"region": "orientation"})
    return pd.concat([pd.DataFrame(rows), cen], ignore_index=True)


# -------------------------------------------------------------------- part 3
MOBILE_PRODUCT = re.compile(
    r"transposase|transposon|\bintegrase\b|recombinase|\bIS\d|insertion sequence"
    r"|\binvertase\b|resolvase|\bTnp\b|IstB", re.I)
# housekeeping recombination/repair proteins the product regex would catch
MOBILE_PRODUCT_EXCLUDE = re.compile(
    r"RecA|RuvX|RuvC|RuvA|RuvB|Holliday|XerC|XerD|RecB|RecD|site-specific "
    r"tyrosine recombinase XerC", re.I)
MOBILE_PFAM = {
    "DDE_Tnp_1", "DDE_Tnp_1_2", "DDE_Tnp_4", "DDE_Tnp_ISL3", "DDE_Tnp_IS240",
    "DEDD_Tnp_IS110", "LZ_Tnp_IS481", "Transposase_20", "HTH_Tnp_1",
    "HTH_Tnp_4", "HTH_OrfB_IS605", "OrfB_IS605", "IstB_IS21",
    "Mu-transpos_C_2", "rve", "rve_3", "Phage_integrase", "Integrase_N_phage",
    "Recombinase", "Resolvase",
}


def mobile_genes(anno):
    a = anno[anno.contig == CHROM].copy()
    prod = a["product"].fillna("")
    byprod = prod.str.contains(MOBILE_PRODUCT) & ~prod.str.contains(MOBILE_PRODUCT_EXCLUDE)
    fam = {}
    with open(PFAM) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) > 4:
                fam.setdefault(f[0], set()).add(f[3])
    hit = {p for p, s in fam.items() if s & MOBILE_PFAM}
    bypfam = a["protein_id"].fillna("").isin(hit)
    a["mobile_by_product"] = byprod.to_numpy()
    a["mobile_by_pfam"] = bypfam.to_numpy()
    a["mobile"] = a.mobile_by_product | a.mobile_by_pfam
    return a


def part3(uni, anno_m, genome, cGC, cN, L, gc_bin, store, flank=2000):
    mob = anno_m[anno_m.mobile]
    lab_body = np.full(L, "background", dtype=object)
    for s, e in zip(mob.start, mob.end):
        lab_body[int(s) - 1:int(e)] = "mobile_body"
    lab_f = np.full(L, "background", dtype=object)
    for s, e in zip(mob.start, mob.end):
        lab_f[max(int(s) - 1 - flank, 0):min(int(e) + flank, L)] = "mobile_body_flank2kb"
    rows = []
    for name, lab in (("mobile_body", lab_body),
                      ("mobile_body_flank2kb", lab_f)):
        n_bp = int((lab == name).sum())
        gcf = float(np.mean(np.frombuffer(genome.encode(), dtype=np.uint8)[lab == name]
                            == ord("G")) +
                    np.mean(np.frombuffer(genome.encode(), dtype=np.uint8)[lab == name]
                            == ord("C")))
        for motif, mark, offs, label in SERIES:
            pos0, plus = R3.motif_positions(genome, motif, offs)
            exp = gc_matched_expectation(pos0, lab, gc_bin, [name])
            obs = int((lab[pos0] == name).sum())
            rows.append(dict(
                part="3a_motif_placement", region=name, series=label,
                n_genes=int(len(mob)), region_bp=n_bp,
                region_gc=round(gcf, 4),
                region_pct_of_genome=round(100 * n_bp / L, 3),
                motif_n_total=len(pos0), observed=obs,
                expected_gcmatched=round(exp[name], 2),
                O_over_E=round(obs / exp[name], 3) if exp[name] else np.nan,
                p_value=f"{binom_two_sided(obs, exp[name], len(pos0)):.3g}",
                mde80_abs_dev_from_1=round(binom_mde(exp[name], len(pos0)), 3)))
        # which gene each base of the region belongs to, for the jackknife
        owner = np.full(L, "", dtype=object)
        for gid, s0, e0 in zip(mob.gene_id, mob.start, mob.end):
            a = int(s0) - 1 if name == "mobile_body" else max(int(s0) - 1 - flank, 0)
            b = int(e0) if name == "mobile_body" else min(int(e0) + flank, L)
            owner[a:b] = gid
        for label, (pos0, st0, Nv, Nm) in uni.items():
            level = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
            gc = M.local_gc(pos0, cGC, cN, L, WIN)
            near = lab[pos0] == name
            if near.sum() < 20:
                rows.append(dict(part="3b_level_ratio", region=name,
                                 series=label, n_near=int(near.sum()),
                                 ratio=np.nan,
                                 note="NOT TESTED: n_near < 20"))
                continue
            r, lo, hi, p, n, sd, mde = summarise_plus(
                M.standardised_ratio(level, Nv, Nm, near, gc, nq=NQ))
            # leave-one-gene-out: IS families are multi-copy, so a handful of
            # near-identical copies could carry the whole effect
            own = owner[pos0]
            genes = sorted(set(own[near]) - {""})
            jk = []
            for gid in genes:
                drop = near & (own == gid)
                nn = near & ~drop
                if nn.sum() < 15:
                    continue
                jk.append(M.summarise(
                    M.standardised_ratio(level, Nv, Nm, nn, gc, nq=NQ))[0])
            rows.append(dict(
                part="3b_level_ratio", region=name, series=label,
                universe_n=int(len(pos0)), n_near=int(near.sum()),
                n_distal=int((~near).sum()), n_cultures=n,
                n_genes_with_positions=len(genes),
                sd_log_ratio=round(sd, 5), mde80_ratio=round(mde, 4),
                mean_level_near=round(float(np.nanmean(level[:, near])), 6),
                mean_level_distal=round(float(np.nanmean(level[:, ~near])), 6),
                ratio=round(r, 4), ci95_lo=round(lo, 4), ci95_hi=round(hi, 4),
                p_value=f"{p:.3g}",
                jackknife_min_ratio=round(min(jk), 4) if jk else np.nan,
                jackknife_max_ratio=round(max(jk), 4) if jk else np.nan))
            # mobile elements sit at 0.686 GC against a 0.721 genome, so most
            # of their positions fall in the lowest local-GC quintile; repeat
            # the ratio at finer matching to show it is not a GC artefact
            gq5 = np.bincount(pd.qcut(gc, NQ, labels=False,
                                      duplicates="drop")[near], minlength=NQ)
            for nq in (10, 20):
                r2, lo2, hi2, p2, _, _, _ = summarise_plus(
                    M.standardised_ratio(level, Nv, Nm, near, gc, nq=nq))
                rows.append(dict(
                    part="3d_gcmatch_sensitivity", region=name, series=label,
                    gc_quantiles=nq, n_near=int(near.sum()),
                    near_positions_per_gc_quintile=",".join(map(str, gq5)),
                    mean_local_gc_near=round(float(gc[near].mean()), 4),
                    mean_local_gc_distal=round(float(gc[~near].mean()), 4),
                    ratio=round(r2, 4), ci95_lo=round(lo2, 4),
                    ci95_hi=round(hi2, 4), p_value=f"{p2:.3g}"))
    cen = census_fraction(store, genome,
                          {"mobile_body": lab_body == "mobile_body",
                           "mobile_body_flank2kb": lab_f == "mobile_body_flank2kb"})
    cen.insert(0, "part", "3c_census_fraction")
    return pd.concat([pd.DataFrame(rows), cen], ignore_index=True)


# -------------------------------------------------------------------- part 4
# what the lost motif_region_GCmatched_261004.csv was recorded as holding
RECORDED_261004 = {
    ("4mC, TGGCCGGC C5", "promoter"): dict(pct=1.43, exp=2.25, ratio=0.64, p="0.02"),
    ("4mC, TGGCCGGC C5", "other"): dict(pct=3.11, exp=5.40, ratio=0.58, p="<0.005"),
    ("4mC, TGGCCGGC C5", "cds"): dict(pct=93.69, exp=89.97, ratio=1.04, p=""),
}


def part4(genome, tss, cds_mask, cGC, cN, L, gc_bin):
    pos_all = np.arange(L, dtype=np.int64)
    d = R3.signed_tss_distance(pos_all + 1, tss["tss_position"].to_numpy(),
                               tss["strand"].to_numpy())
    cls = R3.classify(pos_all, d, cds_mask)
    lab = np.full(L, "other", dtype=object)
    for c in ("other", "cds", "utr5", "promoter"):
        lab[cls[c]] = c
    # alternative precedence, tested only to see whether the lost 2026-10-04
    # table was built with CDS outranking promoter (which would shrink the
    # promoter class to its intergenic part and could explain its 2.25 %)
    prom_raw = (d < 0) & (d >= -300)
    utr_raw = (d >= 0) & (d <= 50)
    lab_cdsfirst = np.full(L, "other", dtype=object)
    lab_cdsfirst[utr_raw & ~cds_mask] = "utr5"
    lab_cdsfirst[prom_raw & ~cds_mask] = "promoter"
    lab_cdsfirst[cds_mask] = "cds"
    classes = ["promoter", "utr5", "cds", "other"]
    rows = []
    for prec, LB in (("canonical_261101 (promoter > 5'UTR > CDS)", lab),
                     ("cds_first (CDS > promoter > 5'UTR)", lab_cdsfirst)):
        for motif, mark, offs, label in SERIES:
            pos0, plus = R3.motif_positions(genome, motif, offs)
            exp = gc_matched_expectation(pos0, LB, gc_bin, classes)
            n = len(pos0)
            for c in classes:
                obs = int((LB[pos0] == c).sum())
                rec = (RECORDED_261004.get((label, c), {})
                       if prec.startswith("canonical") else {})
                rows.append(dict(
                    precedence=prec, series=label, region=c, motif_n_total=n,
                    region_bp=int((LB == c).sum()),
                    region_pct_of_genome=round(100 * (LB == c).sum() / L, 3),
                    observed=obs, pct_observed=round(100 * obs / n, 3),
                    expected_gcmatched=round(exp[c], 2),
                    pct_expected=round(100 * exp[c] / n, 3),
                    O_over_E=round(obs / exp[c], 3) if exp[c] else np.nan,
                    p_value=f"{binom_two_sided(obs, exp[c], n):.3g}",
                    mde80_abs_dev_from_1=round(binom_mde(exp[c], n), 3),
                    recorded_261004_pct=rec.get("pct", ""),
                    recorded_261004_pct_expected=rec.get("exp", ""),
                    recorded_261004_ratio=rec.get("ratio", ""),
                    recorded_261004_p=rec.get("p", ""),
                    matches_record=("" if not rec else
                                    str(abs(round(100 * obs / n, 2) - rec["pct"]) <= 0.05
                                        and abs(round(obs / exp[c], 2) - rec["ratio"]) <= 0.02))))
    return pd.DataFrame(rows), lab


# -------------------------------------------------------------------- part 5
# The region-composition block of value_registry.py (Supp Fig 7, recomputed
# 2026-10-15) states a THIRD region rule, and its CDS-internal percentages for
# TGGCCGGC (93.6 / 93.4 / 93.7) sit on the 93.69 of the lost 2026-10-04 table.
# The lost table was therefore composition of CALLED SITES under that rule, not
# placement of genome motif occurrences under the 2026-11-01 rule. Part 5
# rebuilds exactly that so the record can be collated instead of guessed at.
#   promoter = 1-500 bp upstream of the nearest TSS AND outside any gene body
#   5'UTR    = TSS .. +100
#   cds_internal / intergenic
#   TSS      = Jeong category P or S with a non-null gene_id (3,030)
LEGACY_RECORD = {  # the three rows the manuscript notes carry, TGGCCGGC
    "promoter": dict(pct=1.43, exp=2.25, ratio=0.64),
    "intergenic": dict(pct=3.11, exp=5.40, ratio=0.58),
    "cds_internal": dict(pct=93.69, exp=89.97, ratio=1.04),
}


def part5(genome, tss, cds_mask, gene_mask, cGC, cN, L, gc_bin, store):
    sub = tss[(tss["category"].isin(["P", "S"])) & tss["gene_id"].notna()]
    n_tss = int(sub["tss_position"].nunique())
    pos_all = np.arange(L, dtype=np.int64)
    d = R3.signed_tss_distance(pos_all + 1, sub["tss_position"].to_numpy(),
                               sub["strand"].to_numpy())
    lab = np.full(L, "intergenic", dtype=object)
    lab[cds_mask] = "cds_internal"
    lab[(d >= 0) & (d <= 100)] = "utr5"
    lab[(d <= -1) & (d >= -500) & ~gene_mask] = "promoter"
    classes = ["promoter", "utr5", "cds_internal", "intergenic"]
    rows = []
    for tp in TIMEPOINTS:
        for motif, mark, offs, label in SERIES:
            pos0, plus = R3.motif_positions(genome, motif, offs)
            s = store[(tp, mark)]
            idx = pd.MultiIndex.from_arrays([pos0, plus], names=["start", "plus"])
            q = s.reindex(idx)
            keep = (q["dep"] >= 30).to_numpy(dtype=bool)
            p = pos0[keep]
            called = p[q["lvl"].to_numpy()[keep] >= 0.5]
            exp = gc_matched_expectation(called, lab, gc_bin, classes)
            n = len(called)
            for c in classes:
                obs = int((lab[called] == c).sum())
                rec = (LEGACY_RECORD.get(c, {})
                       if label == "4mC, TGGCCGGC C5" else {})
                rows.append(dict(
                    rule="legacy_suppfig7 (promoter 1-500 bp upstream and "
                         "outside any gene body; 5'UTR TSS..+100; TSS n=%d)" % n_tss,
                    timepoint=tp, series=label, region=c,
                    n_called_sites=n, observed=obs,
                    pct_observed=round(100 * obs / n, 2),
                    expected_gcmatched=round(exp[c], 2),
                    pct_expected=round(100 * exp[c] / n, 2),
                    O_over_E=round(obs / exp[c], 3) if exp[c] else np.nan,
                    p_value=f"{binom_two_sided(obs, exp[c], n):.3g}",
                    mde80_abs_dev_from_1=round(binom_mde(exp[c], n), 3),
                    recorded_261004_pct=rec.get("pct", ""),
                    recorded_261004_pct_expected=rec.get("exp", ""),
                    recorded_261004_ratio=rec.get("ratio", ""),
                    delta_pct_vs_record=(round(100 * obs / n - rec["pct"], 2)
                                         if rec else ""),
                    delta_ratio_vs_record=(round(obs / exp[c] - rec["ratio"], 3)
                                           if rec and exp[c] else "")))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    genome = M.read_reference(M.REF)
    cGC, cN, L = M.cum_gc(genome)
    assert L == 8667507, L
    gc_bin = pd.qcut(M.local_gc(np.arange(L), cGC, cN, L, WIN), NQ,
                     labels=False, duplicates="drop").astype(np.int16)

    tss = pd.read_csv(TSSCSV)
    assert len(tss) == 3570
    tss_all = np.array(sorted(tss["tss_position"].astype(int)))
    tu = pd.read_csv(TUCSV)
    anno = pd.read_csv(ANNO, sep="\t", dtype=str)
    anno["start"] = anno["start"].astype(int)
    anno["end"] = anno["end"].astype(int)
    cds_mask = np.zeros(L, bool)
    gene_mask = np.zeros(L, bool)
    for line in open(GFF):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[0] != CHROM:
            continue
        if f[2] == "CDS":
            cds_mask[int(f[3]) - 1:int(f[4])] = True
        elif f[2] == "gene":
            gene_mask[int(f[3]) - 1:int(f[4])] = True

    uni = build_universes(genome, L)
    store = census_store()

    # estimator equivalence: strata_ratio == M.standardised_ratio on GC quintiles
    pos0, st0, Nv, Nm = uni["6mA, AAGCCCG A1/A2"]
    lv = np.where(Nv > 0, Nm / np.maximum(Nv, 1), np.nan)
    gc = M.local_gc(pos0, cGC, cN, L, WIN)
    near = M.min_distance(pos0 + 1, tss_all) <= WIN
    q = pd.qcut(gc, NQ, labels=False, duplicates="drop").astype(str)
    assert np.allclose(strata_ratio(lv, near, q),
                       M.standardised_ratio(lv, Nv, Nm, near, gc, nq=NQ)), \
        "strata_ratio diverges from the canonical estimator"

    agg = tu_termini(tu)
    t1 = part1(uni, agg, tss_all, cGC, cN, L)
    t1.to_csv(f"{OUT}/region_tu3prime_261102.tsv", sep="\t", index=False)

    ig = intergenic_regions(anno)
    t2 = part2(uni, ig, genome, cGC, cN, L, gc_bin, store)
    t2.to_csv(f"{OUT}/region_intergenic_orientation_261102.tsv", sep="\t", index=False)

    anno_m = mobile_genes(anno)
    t3 = part3(uni, anno_m, genome, cGC, cN, L, gc_bin, store)
    t3.to_csv(f"{OUT}/region_mobile_elements_261102.tsv", sep="\t", index=False)
    anno_m[anno_m.mobile][["gene_id", "start", "end", "strand", "product",
                           "mobile_by_product", "mobile_by_pfam"]].to_csv(
        f"{OUT}/region_mobile_gene_list_261102.tsv", sep="\t", index=False)

    t4, lab4 = part4(genome, tss, cds_mask, cGC, cN, L, gc_bin)
    t4.to_csv(f"{OUT}/region_motif_placement_gcmatched_261102.tsv", sep="\t", index=False)
    t5 = part5(genome, tss, cds_mask, gene_mask, cGC, cN, L, gc_bin, store)
    t5.to_csv(f"{OUT}/region_motif_placement_legacyrule_261102.tsv", sep="\t", index=False)

    print("== part 1  TU 3' end vs TSS ==")
    print(t1[["series", "anchor", "n_near", "ratio", "ci95_lo", "ci95_hi",
              "p_value", "mde80_ratio"]].to_string(index=False))
    print("\n== part 2a  intergenic motif O/E (all lengths) ==")
    print(t2[(t2.part == "2a_motif_placement") & (t2.length_stratum == "all")]
          [["series", "orientation", "observed", "expected_gcmatched",
            "O_over_E", "p_value"]].to_string(index=False))
    print("\n== part 2b  level ratios ==")
    print(t2[t2.part == "2b_level_ratio"][["series", "orientation", "contrast",
          "n_near", "ratio", "ci95_lo", "ci95_hi", "p_value",
          "mde80_ratio"]].to_string(index=False))
    print("\n== part 3  mobile elements ==")
    print(t3[t3.part.isin(["3a_motif_placement", "3b_level_ratio",
                           "3d_gcmatch_sensitivity"])]
          [["part", "region", "series", "n_near", "n_genes_with_positions",
            "observed", "expected_gcmatched", "O_over_E", "ratio", "ci95_lo",
            "ci95_hi", "p_value", "jackknife_min_ratio",
            "jackknife_max_ratio", "gc_quantiles",
            "near_positions_per_gc_quintile", "mean_local_gc_near",
            "mean_local_gc_distal"]].to_string(index=False))
    print("\n== part 4  motif placement vs GC-matched expectation ==")
    print(t4[["precedence", "series", "region", "region_pct_of_genome",
              "pct_observed", "pct_expected", "O_over_E", "p_value",
              "recorded_261004_pct", "recorded_261004_pct_expected",
              "recorded_261004_ratio", "matches_record"]].to_string(index=False))
    print("\n== part 5  legacy-rule collation with the lost 2026-10-04 table ==")
    print(t5[["timepoint", "series", "region", "n_called_sites", "pct_observed",
              "pct_expected", "O_over_E", "p_value", "recorded_261004_pct",
              "recorded_261004_pct_expected", "recorded_261004_ratio",
              "delta_pct_vs_record", "delta_ratio_vs_record"]]
          .to_string(index=False))

    print("\n== census fractions ==")
    cen = pd.concat([t2[t2.part == "2c_census_fraction"],
                     t3[t3.part == "3c_census_fraction"]])
    cen["grp"] = cen.orientation.fillna(cen.region)
    print(cen[["part", "timepoint", "series", "grp", "n_measurable",
               "n_modified", "pct_modified", "n_measurable_all",
               "pct_modified_all"]].to_string(index=False))
    print("\n== intergenic region inventory ==")
    print(ig.groupby("orientation").agg(n=("length", "size"),
                                        bp=("length", "sum"),
                                        median_bp=("length", "median"),
                                        mean_bp=("length", "mean")).to_string())
    print("\n== mobile gene inventory ==")
    print(anno_m.groupby(["mobile_by_product", "mobile_by_pfam"]).size().to_string())
    return t1, t2, t3, t4, t5, ig, anno_m


if __name__ == "__main__":
    main()
