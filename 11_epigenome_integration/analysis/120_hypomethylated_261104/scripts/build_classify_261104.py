#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""120_hypomethylated_261104 step 1 - trinary classification of recognition-sequence
positions into methylated / hypomethylated / indeterminate, and the consistency call.

Lineage: all-read (no-filtering) motif-limited pileup, pileup_nofilt_261015 (Nfail == 0).
NOTE the task brief gave the path as 15_paper_figures/scripts/recompute_261015/
pileup_nofilt_261015/; that directory does not exist. The only pileup_nofilt_261015 on
disk is methyl/260102_M145/analysis/pileup_nofilt_261015 (same filenames <lib>_nofilt.bed),
which is the canon all-read lineage used by nf_r5_recompute_261101.py. Used that.

Outputs (tables/):
  site_classification_261104.tsv   one row per recognition-sequence position
  criterion_reference_261104.tsv   per class x timepoint bulk level, threshold, class counts
  depth_independence_261104.tsv    depth quintile x call, plus downsampling control
  criterion_sensitivity_261104.tsv tau and CI-level sweep
"""
import os, re, resource, numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
TAB = os.path.join(os.path.dirname(HERE), "tables")
os.makedirs(TAB, exist_ok=True)
RNA = "/Users/okaban/bioinfo/rna-seq"
NOFILT = "/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015"
GENOME = f"{RNA}/11_epigenome_integration/data/NC_003888.3.fna"

LIBS = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"],
        "T3": ["3-2", "3-3", "3-4"]}
TPS = ["T1", "T2", "T3"]
BEDCOL = ["chrom", "start", "end", "name", "score", "strand", "ts", "te",
          "color", "Nvalid", "pct", "Nmod", "Ncanon", "Nother", "Ndel",
          "Nfail", "Ndiff", "Nnocall"]
MIN_DEPTH = 30          # pooled Nvalid per timepoint, project convention
TAU_FRAC = 0.5          # hypomethylation threshold = TAU_FRAC x class bulk level
ALPHA = 0.05            # Clopper-Pearson two-sided

# pileup `name` field -> (class label, mark)
CLASSROW = {"21839,GCCGGC,2":   ("TGGCCGGC_C5_4mC", "4mC"),
            "21839,AAGCCCG,4":  ("AAGCCCG_C5_4mC",  "4mC"),
            "a,AAGCCCG,0":      ("AAGCCCG_A1_6mA",  "6mA")}
CLASSES = ["TGGCCGGC_C5_4mC", "AAGCCCG_C5_4mC", "AAGCCCG_A1_6mA"]


def genome_seq():
    return "".join(l.strip() for l in open(GENOME) if not l.startswith(">")).upper()


def find_all(s, sub):
    out, i = [], s.find(sub)
    while i != -1:
        out.append(i); i = s.find(sub, i + 1)
    return out


def site_universe(g):
    """0-based (pos, strand) -> class, for the three declared recognition positions."""
    uni = {}
    # TGGCCGGC, modified C is position 5 of the 8-mer (0-based offset 4)
    for i in find_all(g, "TGGCCGGC"):
        uni[(i + 4, "+")] = "TGGCCGGC_C5_4mC"
    for i in find_all(g, "GCCGGCCA"):          # TGGCCGGC on the minus strand
        uni[(i + 3, "-")] = "TGGCCGGC_C5_4mC"
    # AAGCCCG: C5 (offset 4) 4mC and A1 (offset 0) 6mA
    for i in find_all(g, "AAGCCCG"):
        uni[(i + 4, "+")] = "AAGCCCG_C5_4mC"
        uni[(i + 0, "+")] = "AAGCCCG_A1_6mA"
    for j in find_all(g, "CGGGCTT"):           # AAGCCCG on the minus strand
        uni[(j + 2, "-")] = "AAGCCCG_C5_4mC"
        uni[(j + 6, "-")] = "AAGCCCG_A1_6mA"
    return uni


def load_pooled(uni):
    """pooled (Nvalid, Nmod) per (pos, strand, class) x timepoint."""
    frames = []
    for tp, libs in LIBS.items():
        acc = None
        for lib in libs:
            d = pd.read_csv(f"{NOFILT}/{lib}_nofilt.bed", sep="\t", header=None,
                            names=BEDCOL, usecols=["start", "name", "strand",
                                                   "Nvalid", "Nmod", "Nfail"])
            assert int(d["Nfail"].sum()) == 0, f"{lib}: not an all-read pileup"
            d = d[d["name"].isin(CLASSROW)].copy()
            d["cls"] = d["name"].map(lambda s: CLASSROW[s][0])
            # the pileup emits GCCGGC,2 at every GCCGGC; keep only declared positions
            key = list(zip(d["start"], d["strand"]))
            d["uni"] = [uni.get(k) for k in key]
            d = d[d["uni"] == d["cls"]]
            f = d.groupby(["start", "strand", "cls"])[["Nvalid", "Nmod"]].sum()
            acc = f if acc is None else acc.add(f, fill_value=0)
            del d, f
        acc = acc.reset_index()
        acc["timepoint"] = tp
        frames.append(acc)
    out = pd.concat(frames, ignore_index=True)
    out[["Nvalid", "Nmod"]] = out[["Nvalid", "Nmod"]].astype(int)
    return out


def cp_interval(k, n, alpha=ALPHA):
    """Clopper-Pearson two-sided interval, vectorised."""
    k = np.asarray(k, float); n = np.asarray(n, float)
    lo = np.where(k > 0, stats.beta.ppf(alpha / 2, np.maximum(k, 1e-9), n - k + 1), 0.0)
    hi = np.where(k < n, stats.beta.isf(alpha / 2, k + 1, np.maximum(n - k, 1e-9)), 1.0)
    return np.nan_to_num(lo, nan=0.0), np.nan_to_num(hi, nan=1.0)


def classify(df, tau_frac=TAU_FRAC, alpha=ALPHA):
    """Adds level, ci_lo, ci_hi, call. Threshold is per class x timepoint."""
    df = df.copy()
    df["level"] = df["Nmod"] / df["Nvalid"]
    lo, hi = cp_interval(df["Nmod"].values, df["Nvalid"].values, alpha)
    df["ci_lo"], df["ci_hi"] = lo, hi
    ref = (df.groupby(["cls", "timepoint"])[["Nmod", "Nvalid"]].sum()
             .assign(bulk=lambda x: x["Nmod"] / x["Nvalid"]))
    ref["tau"] = tau_frac * ref["bulk"]
    df = df.merge(ref[["bulk", "tau"]], left_on=["cls", "timepoint"],
                  right_index=True, how="left")
    df["call"] = np.where(df["ci_hi"] < df["tau"], "hypomethylated",
                 np.where(df["ci_lo"] > df["tau"], "methylated", "indeterminate"))
    return df, ref.reset_index()


def main():
    g = genome_seq()
    uni = site_universe(g)
    print("genome bp", len(g), "universe positions", len(uni))
    raw = load_pooled(uni)

    # measurability: pooled depth floor met at EVERY timepoint
    piv = raw.pivot_table(index=["start", "strand", "cls"], columns="timepoint",
                          values="Nvalid", fill_value=0)
    ok = piv.reindex(columns=TPS, fill_value=0).min(axis=1) >= MIN_DEPTH
    keep = set(ok[ok].index)
    raw["measurable"] = [(s, st, c) in keep for s, st, c in
                         zip(raw["start"], raw["strand"], raw["cls"])]
    meas = raw[raw["measurable"]].copy()
    df, ref = classify(meas)

    # ---- per-site wide table with the consistency call -------------------
    w = df.pivot_table(index=["start", "strand", "cls"], columns="timepoint",
                       values=["Nvalid", "Nmod", "level", "ci_lo", "ci_hi"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    calls = df.pivot_table(index=["start", "strand", "cls"], columns="timepoint",
                           values="call", aggfunc="first")
    w = w.join(calls.rename(columns={t: f"call_{t}" for t in TPS}))
    cc = w[[f"call_{t}" for t in TPS]]
    w["n_hypo_tp"] = (cc == "hypomethylated").sum(axis=1)
    w["n_meth_tp"] = (cc == "methylated").sum(axis=1)
    w["consistent"] = np.where(w["n_hypo_tp"] == 3, "hypomethylated",
                      np.where(w["n_meth_tp"] == 3, "methylated", "variable"))
    w = w.reset_index().rename(columns={"start": "pos0"})
    w["pos1"] = w["pos0"] + 1
    w.to_csv(f"{TAB}/site_classification_261104.tsv", sep="\t", index=False)

    # ---- criterion reference table --------------------------------------
    cnt = (df.groupby(["cls", "timepoint"])["call"].value_counts().unstack(fill_value=0)
             .reindex(columns=["methylated", "hypomethylated", "indeterminate"], fill_value=0))
    uni_n = pd.Series({c: sum(1 for v in uni.values() if v == c) for c in CLASSES},
                      name="n_universe")
    ref = ref.merge(cnt.reset_index(), on=["cls", "timepoint"])
    ref["n_measurable"] = ref["cls"].map(df.groupby("cls")["start"].nunique() * 0 +
                                         df[df.timepoint == "T1"].groupby("cls").size())
    ref["n_universe"] = ref["cls"].map(uni_n)
    ref.to_csv(f"{TAB}/criterion_reference_261104.tsv", sep="\t", index=False)

    # ---- depth independence ---------------------------------------------
    rows = []
    for c in CLASSES:
        sub = df[df["cls"] == c]
        for t in TPS:
            s = sub[sub["timepoint"] == t].copy()
            if len(s) < 25:
                continue
            s["q"] = pd.qcut(s["Nvalid"], 5, labels=[1, 2, 3, 4, 5], duplicates="drop")
            for q, grp in s.groupby("q", observed=True):
                rows.append(dict(cls=c, timepoint=t, depth_quintile=int(q), n=len(grp),
                                 depth_median=float(grp["Nvalid"].median()),
                                 pct_hypo=100 * (grp["call"] == "hypomethylated").mean(),
                                 pct_meth=100 * (grp["call"] == "methylated").mean(),
                                 pct_indet=100 * (grp["call"] == "indeterminate").mean(),
                                 level_median=float(grp["level"].median())))
            rho, p = stats.spearmanr(s["Nvalid"], s["level"])
            rows.append(dict(cls=c, timepoint=t, depth_quintile=-1, n=len(s),
                             depth_median=float(s["Nvalid"].median()),
                             pct_hypo=np.nan, pct_meth=np.nan, pct_indet=np.nan,
                             level_median=float(rho), spearman_p=float(p)))
    dep = pd.DataFrame(rows)

    # downsampling control: methylated sites thinned to the hypo median depth
    rng = np.random.default_rng(20261104)
    ds = []
    for c in CLASSES:
        for t in TPS:
            s = df[(df["cls"] == c) & (df["timepoint"] == t)]
            h = s[s["call"] == "hypomethylated"]
            m = s[s["call"] == "methylated"]
            if len(h) == 0 or len(m) == 0:
                continue
            nd = int(np.median(h["Nvalid"]))
            tau = float(s["tau"].iloc[0])
            n_new = np.minimum(m["Nvalid"].values, nd)
            k_new = rng.hypergeometric(m["Nmod"].values,
                                       m["Nvalid"].values - m["Nmod"].values, n_new)
            lo, hi = cp_interval(k_new, n_new)
            ds.append(dict(cls=c, timepoint=t, target_depth=nd, n_methylated=len(m),
                           flip_to_hypo=int((hi < tau).sum()),
                           to_indeterminate=int(((hi >= tau) & (lo <= tau)).sum()),
                           stay_methylated=int((lo > tau).sum())))
    pd.DataFrame(ds).to_csv(f"{TAB}/depth_downsample_control_261104.tsv", sep="\t", index=False)
    dep.to_csv(f"{TAB}/depth_independence_261104.tsv", sep="\t", index=False)

    # ---- sensitivity of the criterion -----------------------------------
    sens = []
    for tf in [0.3, 0.5, 0.7]:
        for al in [0.05, 0.01]:
            d2, _ = classify(meas, tau_frac=tf, alpha=al)
            cw = d2.pivot_table(index=["start", "strand", "cls"], columns="timepoint",
                                values="call", aggfunc="first")
            for c in CLASSES:
                cs = cw[cw.index.get_level_values("cls") == c]
                sens.append(dict(tau_frac=tf, alpha=al, cls=c, n_measurable=len(cs),
                                 n_hypo_consistent=int((cs == "hypomethylated").all(axis=1).sum()),
                                 n_hypo_any_tp=int((cs == "hypomethylated").any(axis=1).sum()),
                                 n_meth_consistent=int((cs == "methylated").all(axis=1).sum())))
    pd.DataFrame(sens).to_csv(f"{TAB}/criterion_sensitivity_261104.tsv", sep="\t", index=False)

    print(w.groupby(["cls", "consistent"]).size().to_string())
    print("peak RSS GB", round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 3))


if __name__ == "__main__":
    main()
