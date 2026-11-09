#!/usr/bin/env python
"""E2-3: is the per-read tract-length distribution multimodal (ON/OFF mixture in the
culture) or is it the unimodal nanopore length-error distribution?
Decision rule: every statistic computed on a candidate tract is ranked against the
same statistic computed on the matched control tracts (same base, same reference
length, outside MTase/BREX genes), which by construction carry a single allele."""
import sys, resource
import numpy as np, pandas as pd
import diptest

OUT = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/123_phasevar_261104"
MINDEPTH = 30
rng = np.random.default_rng(20261104)

def stats_for(v, reflen):
    v = np.asarray(v, float)
    n = len(v)
    vals, cnt = np.unique(v, return_counts=True)
    mode = vals[np.argmax(cnt)]
    f_mode = cnt.max()/n
    # second mode: most frequent length at least 1 bp away from the modal length
    off = vals != mode
    f_second = (cnt[off].max()/n) if off.any() else 0.0
    second = vals[off][np.argmax(cnt[off])] if off.any() else np.nan
    # dip test needs jitter: lengths are integers
    jit = v + rng.normal(0, 0.12, n)
    d, p = diptest.diptest(jit)
    return dict(n=n, mode=mode, mode_eq_ref=int(mode == reflen), f_mode=f_mode,
                second_mode=second, f_second=f_second, mean=v.mean(), sd=v.std(ddof=1),
                f_nonmodal=1-f_mode, dip=d, dip_p=p)

def main():
    df = pd.read_csv(f"{OUT}/tables/e2_2_perread_lengths_261104.tsv.gz", sep="\t")
    loci = pd.read_csv(f"{OUT}/tables/e2_2_loci_261104.tsv", sep="\t")
    info = loci.set_index("locus_id")

    # metric = longest run of the tract base inside the inter-anchor segment
    rows = []
    for lid, sub in df.groupby("locus_id"):
        if len(sub) < MINDEPTH: continue
        rl = int(info.loc[lid, "ref_len"])
        s = stats_for(sub.max_run.values, rl)
        s.update(locus_id=lid, locus_class=info.loc[lid, "locus_class"],
                 base=info.loc[lid, "base"], ref_len=rl,
                 genes=info.loc[lid, "genes"], sco=info.loc[lid, "sco"],
                 gname=info.loc[lid, "gname"], region=info.loc[lid, "region"],
                 f_pure=sub.pure.mean())
        rows.append(s)
    st = pd.DataFrame(rows)

    # empirical rank of each candidate against controls matched on base+ref_len
    out = []
    for _, r in st[st.locus_class == "candidate"].iterrows():
        ctl = st[(st.locus_class == "control") & (st.ref_len == r.ref_len) &
                 (st.base == r.base)]
        if len(ctl) < 5:
            ctl = st[(st.locus_class == "control") & (st.ref_len == r.ref_len)]
        d = r.to_dict()
        d["n_controls"] = len(ctl)
        for m in ["f_second", "f_nonmodal", "sd", "dip"]:
            d[f"p_{m}"] = ((ctl[m].values >= r[m]).sum()+1)/(len(ctl)+1)
            d[f"ctl_{m}_median"] = ctl[m].median()
            d[f"ctl_{m}_max"] = ctl[m].max()
        # minimum detectable minor-allele fraction: the control maximum f_second
        d["min_detectable_minor_fraction"] = ctl.f_second.max()
        out.append(d)
    cand = pd.DataFrame(out)

    st.to_csv(f"{OUT}/tables/e2_3_locus_stats_261104.tsv", sep="\t", index=False)

    # per-length histograms for candidates (for the write-up)
    h = (df[df.locus_id.isin(st[st.locus_class=="candidate"].locus_id)]
           .groupby(["locus_id","max_run"]).size().rename("n_reads").reset_index())
    h.to_csv(f"{OUT}/tables/e2_3_candidate_length_histograms_261104.tsv", sep="\t", index=False)

    # multimodal call: candidate exceeds EVERY matched control on f_second
    # AND dip test rejects unimodality at the same threshold controls do not
    if len(cand):
        cand["multimodal_call"] = ((cand.p_f_second <= 0.05) & (cand.dip_p < 0.05)).astype(int)
        cand.to_csv(f"{OUT}/tables/e2_3_candidate_vs_control_261104.tsv", sep="\t", index=False)
        nmm = int(cand.multimodal_call.sum())
    else:
        nmm = 0
    with open(f"{OUT}/tables/e2_3_summary_261104.tsv","w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"loci_with_depth_ge{MINDEPTH}\t{len(st)}\n")
        fh.write(f"candidates_tested\t{(st.locus_class=='candidate').sum()}\n")
        fh.write(f"controls_tested\t{(st.locus_class=='control').sum()}\n")
        fh.write(f"median_depth_candidates\t{st[st.locus_class=='candidate'].n.median():.0f}\n")
        fh.write(f"median_depth_controls\t{st[st.locus_class=='control'].n.median():.0f}\n")
        c = st[st.locus_class=='control']
        fh.write(f"control_f_mode_median\t{c.f_mode.median():.3f}\n")
        fh.write(f"control_f_second_median\t{c.f_second.median():.3f}\n")
        fh.write(f"control_f_second_max\t{c.f_second.max():.3f}\n")
        fh.write(f"control_dip_p_lt0.05_fraction\t{(c.dip_p<0.05).mean():.3f}\n")
        fh.write(f"control_mode_eq_ref_fraction\t{c.mode_eq_ref.mean():.3f}\n")
        fh.write(f"n_multimodal_candidates\t{nmm}\n")
    print(open(f"{OUT}/tables/e2_3_summary_261104.tsv").read())
    if len(cand):
        cols=["locus_id","sco","gname","region","ref_len","n","mode","mode_eq_ref","f_mode",
              "second_mode","f_second","p_f_second","ctl_f_second_median","ctl_f_second_max",
              "sd","ctl_sd_median","dip_p","min_detectable_minor_fraction","multimodal_call"]
        pd.set_option("display.width",260)
        print(cand[cols].to_string(index=False))
    print("peak_rss_gb", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3)

main()
