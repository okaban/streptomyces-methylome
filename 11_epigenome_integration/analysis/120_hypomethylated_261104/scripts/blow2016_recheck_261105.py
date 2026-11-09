#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""120_hypomethylated_261104 step 4 - independent re-check under the Blow et al. 2016
PLoS Genet (10.1371/journal.pgen.1005854) Methods design.

Three things change relative to steps 1-3, and the point is to see which one (if any)
moves the conclusion:
  (i)  CRITERION. Blow sets two point thresholds interpolated between the motif mean and
       the non-motif mean, with both means taken over the MIDDLE 60% of ranked values:
         unmethylated threshold = 0.1*mean_motif + 0.9*mean_nonmotif
         methylated   threshold = 0.5*mean_motif + 0.5*mean_nonmotif
       Ported to nanopore occupancy (no ipdR): the non-motif reference is the background
       modified fraction b at non-motif positions of the same base, so
         tau_un = 0.1*mu_trim + 0.9*b,  tau_me = 0.5*mu_trim + 0.5*b.
       The parent's port sets b = 0; b is swept because the all-read pileup is
       motif-limited and carries no non-motif rows to measure b from (declared).
       Blow uses POINT estimates, no confidence interval - implemented as stated.
  (ii) WINDOW. Blow's regulatory region is CDS start -100 .. +50, not -300..TSS.
  (iii)NULL. Blow compares against 1,000 draws of composition-matched random sites,
       and separately uses the methylated motifs as the control for the unmethylated ones.

Reuses tables/site_classification_261104.tsv - no pileup is re-read.
Outputs (tables/): blow_criterion_261105.tsv, blow_window_null_261105.tsv,
  blow_cog_261105.tsv, blow_design_grid_261105.tsv
"""
import os, re, resource, numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(HERE), "tables")
RNA = "/Users/okaban/bioinfo/rna-seq"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"
ANNOT = f"{RNA}/05_annotation/analysis/05_annotation_260128_v1/tables/gene_annotation_basic.tsv"
COGF = f"{RNA}/12_supplementary_figures/analysis/12_supplementary_260202_v1/tables/gene_COG_classification.tsv"
TPS = ["T1", "T2", "T3"]
NPERM = 1000
TRIM = 0.20                      # middle 60% of ranked levels
UP, DOWN = 100, 50               # Blow regulatory window around CDS start
# class -> (motif on the + strand, 0-based offset of the modified base)
MOTIF = {"TGGCCGGC_C5_4mC": ("TGGCCGGC", 4),
         "AAGCCCG_C5_4mC":  ("AAGCCCG", 4),
         "AAGCCCG_A1_6mA":  ("AAGCCCG", 0)}
CLASSES = list(MOTIF)


def genome_seq():
    return "".join(l.strip() for l in open(GENOME) if not l.startswith(">")).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def trimmed_mean(x, trim=TRIM):
    x = np.sort(np.asarray(x, float))
    lo, hi = int(np.floor(len(x) * trim)), int(np.ceil(len(x) * (1 - trim)))
    return float(x[lo:hi].mean())


def min_detectable_or(n_h, n_bg, k_bg, alpha=0.05):
    if n_h == 0:
        return np.nan, np.nan
    for k in range(n_h + 1):
        if stats.fisher_exact([[k, n_h - k], [k_bg, n_bg - k_bg]],
                              alternative="greater")[1] < alpha:
            return k, ((k + .5) * (n_bg - k_bg + .5)) / ((n_h - k + .5) * (k_bg + .5))
    return np.nan, np.inf


# ------------------------------------------------------------------ windows
def regulatory_mask(g, ann, up=UP, down=DOWN):
    """Boolean over genome: within `up` bp upstream .. `down` bp downstream of a CDS start."""
    m = np.zeros(len(g) + 1, bool)
    for s, e, st in zip(ann["start"], ann["end"], ann["strand"]):
        c = s if st == "+" else e
        a, b = (c - up, c + down) if st == "+" else (c - down, c + up)
        m[max(1, a):min(len(g), b) + 1] = True
    return m


def tss_promoter_mask(g, tss, up=300):
    """Boolean over genome: within `up` bp upstream of an experimental TSS.
    Must match step 2's promoter definition exactly (build_annotation in test_bias_261104.py)."""
    m = np.zeros(len(g) + 1, bool)
    for t, st in zip(tss["tss_position"], tss["strand"]):
        a, b = (t - up, t - 1) if st == "+" else (t + 1, t + up)
        m[max(1, a):min(len(g), b) + 1] = True
    return m


def regulatory_gene(g, ann, up=UP, down=DOWN):
    """1-based position -> gene_id of the CDS whose regulatory window it falls in."""
    out = {}
    for s, e, st, gid in zip(ann["start"], ann["end"], ann["strand"], ann["gene_id"]):
        c = s if st == "+" else e
        a, b = (c - up, c + down) if st == "+" else (c - down, c + up)
        for p in range(max(1, a), min(len(g), b) + 1):
            out.setdefault(p, gid)
    return out


# ------------------------------------------------- composition-matched pool
def comp_pool(g, motif):
    """1-based start positions whose L-mer has the same base multiset as `motif`."""
    L = len(motif)
    arr = np.frombuffer(g.encode(), dtype=np.uint8)
    counts = {}
    for b in b"ACGT":
        v = (arr == b).astype(np.int32)
        cs = np.concatenate([[0], np.cumsum(v)])
        counts[b] = cs[L:] - cs[:-L]            # count of b in each L-mer window
    want = {b: motif.count(chr(b)) for b in b"ACGT"}
    ok = np.ones(len(counts[b"A"[0]]), bool)
    for b in b"ACGT":
        ok &= counts[b] == want[b]
    return np.flatnonzero(ok) + 1               # 1-based start


def main():
    g = genome_seq()
    ann = pd.read_csv(ANNOT, sep="\t")
    ann = ann[(ann["contig"] == "NC_003888.3") & (ann["gene_biotype"] == "protein_coding")]
    cog = pd.read_csv(COGF, sep="\t")[["gene_id", "COG_category"]].drop_duplicates()
    w = pd.read_csv(f"{TAB}/site_classification_261104.tsv", sep="\t")
    prev = pd.read_csv(f"{TAB}/site_annotation_261104.tsv", sep="\t")
    w = w.drop(columns=[c for c in ("hypo_strict", "hypo_relaxed") if c in w.columns])
    w = w.merge(prev[["pos1", "strand", "cls", "hypo_strict", "hypo_relaxed",
                      "locus_class"]].rename(columns={"locus_class": "locus_TSS300"}),
                on=["pos1", "strand", "cls"], how="left")
    tss = pd.read_csv(f"{RNA}/11_epigenome_integration/analysis/18_tss_analyses/"
                      "jeong2016_all_tss.csv")
    regmask = regulatory_mask(g, ann)
    tssmask = tss_promoter_mask(g, tss)
    reggene = regulatory_gene(g, ann)
    w["in_cds_window"] = regmask[w["pos1"].values]
    w["cds_window_gene"] = [reggene.get(p) for p in w["pos1"]]

    rng = np.random.default_rng(20261105)

    # --------------------------------------------- (i) Blow criterion
    crit, callcols = [], {}
    for cls in CLASSES:
        sub = w[w["cls"] == cls]
        for b in [0.0, 0.01, 0.02, 0.05]:
            calls = {}
            for t in TPS:
                mu = trimmed_mean(sub[f"level_{t}"])
                tun, tme = 0.1 * mu + 0.9 * b, 0.5 * mu + 0.5 * b
                lv = sub[f"level_{t}"].values
                calls[t] = np.where(lv < tun, "U", np.where(lv > tme, "M", "I"))
                if b == 0.0:
                    crit.append(dict(cls=cls, timepoint=t, b_nonmotif=b,
                                     mu_trimmed60=mu,
                                     mu_readweighted=sub[f"level_{t}"].mean(),
                                     tau_unmethylated=tun, tau_methylated=tme,
                                     n_U=int((calls[t] == "U").sum()),
                                     n_M=int((calls[t] == "M").sum()),
                                     n_I=int((calls[t] == "I").sum())))
            C = np.vstack([calls[t] for t in TPS])
            u3 = (C == "U").all(0)
            if b == 0.0:
                callcols[cls] = (sub.index.values, u3, (C == "M").all(0))
            crit.append(dict(cls=cls, timepoint="consistent_3of3", b_nonmotif=b,
                             n_U=int(u3.sum()), n_M=int((C == "M").all(0).sum()),
                             n_I=int(len(sub) - u3.sum() - (C == "M").all(0).sum())))
    w["blow_unmeth"] = False; w["blow_meth"] = False
    for cls, (idx, u, m) in callcols.items():
        w.loc[idx, "blow_unmeth"] = u; w.loc[idx, "blow_meth"] = m
    pd.DataFrame(crit).to_csv(f"{TAB}/blow_criterion_261105.tsv", sep="\t", index=False)

    # ------------------------- (ii)+(iii) window x null design grid
    rows = []
    # name -> (observed-side column or None for locus_TSS300=='promoter', NULL-side mask).
    # The null mask MUST match the window the observed side is scored on; using one mask
    # for both windows was the defect found in review on 2026-11-05.
    WINDOWS = {"Blow_CDS-100+50": ("in_cds_window", regmask),
               "ours_TSS-300": (None, tssmask)}
    assert np.array_equal(tssmask[w["pos1"].values],
                          (w["locus_TSS300"].values == "promoter")), \
        "TSS mask does not reproduce step 2's promoter assignment"
    DEFS = {"Blow_unmeth_3of3": "blow_unmeth", "ours_hypo_strict": "hypo_strict",
            "ours_hypo_relaxed": "hypo_relaxed"}
    for cls in CLASSES:
        sub = w[w["cls"] == cls].copy()
        motif, off = MOTIF[cls]
        # composition-matched pool, strand-aware, restricted to the two compositions
        pool_p = comp_pool(g, motif) + off
        pool_m = comp_pool(g, rc(motif)) + (len(motif) - 1 - off)
        for wname, (wcol, wmask) in WINDOWS.items():
            inw = (sub[wcol].values if wcol else
                   (sub["locus_TSS300"].values == "promoter"))
            for dname, dcol in DEFS.items():
                h = sub[dcol].values.astype(bool)
                nh, kh = int(h.sum()), int((h & inw).sum())
                # --- null A: composition-matched random sites (Blow)
                npl, nmi = int((sub["strand"] == "+").sum()), int((sub["strand"] == "-").sum())
                frac_p = npl / max(npl + nmi, 1)
                def draw(n, mask=wmask):
                    k1 = rng.binomial(n, frac_p)
                    a = rng.choice(pool_p, k1, replace=True)
                    bb = rng.choice(pool_m, n - k1, replace=True)
                    s = np.concatenate([a, bb])
                    s = s[(s >= 1) & (s <= len(g))]
                    return mask[s].mean() if len(s) else np.nan
                if nh > 0:
                    nulls = np.array([draw(nh) for _ in range(NPERM)])
                    obs = kh / nh
                    p_comp = (1 + (nulls >= obs).sum()) / (NPERM + 1)
                else:
                    nulls, obs, p_comp = np.array([np.nan]), np.nan, np.nan
                # --- null B: the methylated motifs of the same recognition sequence
                ctrl = sub["blow_meth"].values if dname.startswith("Blow") else ~h
                nb, kb = int(ctrl.sum()), int((ctrl & inw).sum())
                p_f = stats.fisher_exact([[kh, nh - kh], [kb, nb - kb]])[1] if nh and nb else np.nan
                orr = (((kh + .5) * (nb - kb + .5)) / ((nh - kh + .5) * (kb + .5))
                       if nh and nb else np.nan)
                kc, mor = min_detectable_or(nh, nb, kb)
                rows.append(dict(cls=cls, window=wname, definition=dname,
                                 n_sites=nh, k_in_window=kh, pct_obs=100 * obs if nh else np.nan,
                                 null_comp_mean_pct=100 * np.nanmean(nulls),
                                 p_composition_null=p_comp,
                                 n_control=nb, pct_control=100 * kb / nb if nb else np.nan,
                                 OR_vs_control=orr, p_fisher_control=p_f,
                                 k_crit=kc, min_detectable_OR=mor))
        # all measurable sites of the class vs composition null (Blow's first test)
        for wname, (wcol, wmask) in WINDOWS.items():
            inw = (sub[wcol].values if wcol else (sub["locus_TSS300"].values == "promoter"))
            n = len(sub)
            k1 = rng.binomial(n, (sub["strand"] == "+").mean(), NPERM)
            nulls = np.empty(NPERM)
            for i in range(NPERM):
                s = np.concatenate([rng.choice(pool_p, k1[i], replace=True),
                                    rng.choice(pool_m, n - k1[i], replace=True)])
                nulls[i] = wmask[s].mean()
            obs = inw.mean()
            rows.append(dict(cls=cls, window=wname, definition="ALL_measurable_motif",
                             n_sites=n, k_in_window=int(inw.sum()), pct_obs=100 * obs,
                             null_comp_mean_pct=100 * nulls.mean(),
                             p_composition_null=(1 + (nulls >= obs).sum()) / (NPERM + 1),
                             n_control=np.nan, pct_control=np.nan, OR_vs_control=np.nan,
                             p_fisher_control=np.nan, k_crit=np.nan, min_detectable_OR=np.nan))
    grid = pd.DataFrame(rows)
    for (c, wn), gg in grid.groupby(["cls", "window"]):
        grid.loc[gg.index, "p_fisher_BH"] = stats.false_discovery_control(
            gg["p_fisher_control"].fillna(1).values)
    grid.to_csv(f"{TAB}/blow_design_grid_261105.tsv", sep="\t", index=False)

    # --------------------------------------------- (iv) COG-stratified repeat
    cw = w.dropna(subset=["cds_window_gene"]).merge(
        cog, left_on="cds_window_gene", right_on="gene_id", how="left")
    crows = []
    for dname, dcol in [("Blow_unmeth_3of3", "blow_unmeth"),
                        ("ours_hypo_strict", "hypo_strict"),
                        ("ours_hypo_relaxed", "hypo_relaxed")]:
        for cat, grp in cw.groupby("COG_category"):
            nh = int(grp[dcol].sum())
            ctrl = grp["blow_meth"] if dname.startswith("Blow") else ~grp[dcol]
            nb = int(ctrl.sum())
            if nb < 5:
                continue
            # enrichment of this COG among window sites of that call, vs the control call
            allh = int(w[dcol].sum()); allb = int((w["blow_meth"] if dname.startswith("Blow")
                                                   else ~w[dcol]).sum())
            p = stats.fisher_exact([[nh, allh - nh], [nb, allb - nb]])[1]
            orr = ((nh + .5) * (allb - nb + .5)) / ((allh - nh + .5) * (nb + .5))
            kc, mor = min_detectable_or(allh, allb, nb)
            crows.append(dict(definition=dname, COG_category=cat, n_window_sites=len(grp),
                              k_call=nh, n_call_total=allh, k_control=nb,
                              n_control_total=allb, OR=orr, p_fisher=p,
                              min_detectable_OR=mor))
    cdf = pd.DataFrame(crows)
    if len(cdf):
        for d, gg in cdf.groupby("definition"):
            cdf.loc[gg.index, "p_BH"] = stats.false_discovery_control(gg["p_fisher"].values)
    cdf.to_csv(f"{TAB}/blow_cog_261105.tsv", sep="\t", index=False)
    w.to_csv(f"{TAB}/blow_site_calls_261105.tsv", sep="\t", index=False)

    print(pd.DataFrame(crit).query("timepoint=='consistent_3of3'")
          [["cls", "b_nonmotif", "n_U", "n_M"]].to_string(index=False))
    print("\nsig composition null (p<0.05):",
          int((grid["p_composition_null"] < .05).sum()), "/", len(grid))
    print("sig control Fisher after BH:", int((grid["p_fisher_BH"] < .05).sum()))
    print("sig COG after BH:", int((cdf["p_BH"] < .05).sum()) if len(cdf) else 0, "/", len(cdf))
    print("peak RSS GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 3))


if __name__ == "__main__":
    main()
