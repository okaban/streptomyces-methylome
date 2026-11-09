#!/usr/bin/env python
"""E2-3b: detection bound. How large a second tract-length allele (ON/OFF minor
fraction) would have been detectable at each candidate locus?

Null is control-to-control, not multinomial: each control tract carries a single
allele, so the spread of G/n statistics among controls already contains both
nanopore length error and locus-to-locus sequence-context variation."""
import sys, resource
import numpy as np, pandas as pd

OUT = "."
rng = np.random.default_rng(20261104)
ALPHA, POWER, NSIM = 0.05, 0.80, 400
LMIN, LMAX = 2, 14

def hist(v):
    h = np.zeros(LMAX - LMIN + 1)
    for x in v:
        x = int(min(max(x, LMIN), LMAX)); h[x - LMIN] += 1
    return h

def kernel(hs):
    s = hs.sum(axis=0); return s / s.sum()

def gstat(obs, p):
    p = np.where(p <= 0, 1e-9, p); e = p * obs.sum()
    m = obs > 0
    return 2 * np.sum(obs[m] * np.log(obs[m] / e[m]))

def main():
    df = pd.read_csv(f"{OUT}/tables/e2_2_perread_lengths_261104.tsv.gz", sep="\t")
    loci = pd.read_csv(f"{OUT}/tables/e2_2_loci_261104.tsv", sep="\t")
    info = loci.set_index("locus_id")
    H = {lid: hist(sub.max_run.values) for lid, sub in df.groupby("locus_id")}

    ctl = {}
    for lid in H:
        if info.loc[lid, "locus_class"] != "control": continue
        ctl.setdefault((info.loc[lid, "base"], int(info.loc[lid, "ref_len"])), []).append(lid)

    rows = []
    for lid in [l for l in H if info.loc[l, "locus_class"] == "candidate"]:
        base, L = info.loc[lid, "base"], int(info.loc[lid, "ref_len"])
        key = (base, L)
        if key not in ctl or len(ctl[key]) < 3:
            key = [k for k in ctl if k[1] == L]
            key = key[0] if key else None
        if key is None: continue
        cl = ctl[key]
        P = kernel(np.array([H[c] for c in cl]))

        # control-to-control null: each control vs kernel of the others
        nullg = []
        for c in cl:
            others = [x for x in cl if x != c]
            Pm = kernel(np.array([H[x] for x in others]))
            nullg.append(gstat(H[c], Pm) / H[c].sum())
        nullg = np.array(nullg)
        thr = np.percentile(nullg, 100 * (1 - ALPHA))

        obs = H[lid]; n = int(obs.sum())
        gobs = gstat(obs, P) / n
        p_emp = ((nullg >= gobs).sum() + 1) / (len(nullg) + 1)

        # minimum detectable minor-allele fraction, allele L-1 and L+1
        mdf = {}
        for alt in (L - 1, L + 1):
            akey = [k for k in ctl if k[1] == alt and k[0] == base] or \
                   [k for k in ctl if k[1] == alt]
            if not akey: mdf[alt] = np.nan; continue
            Q = kernel(np.array([H[c] for c in ctl[akey[0]]]))
            lo = np.nan
            for f in np.arange(0.01, 0.61, 0.01):
                mix = (1 - f) * P + f * Q
                draws = rng.multinomial(n, mix, size=NSIM)
                g = np.array([gstat(d.astype(float), P) / n for d in draws])
                if (g > thr).mean() >= POWER: lo = round(float(f), 3); break
            mdf[alt] = lo
        rows.append(dict(locus_id=lid, sco=info.loc[lid,"sco"], gname=info.loc[lid,"gname"],
            region=info.loc[lid,"region"], base=base, ref_len=L, n_reads=n,
            n_controls=len(cl), G_per_read=round(gobs,5),
            ctl_null_median=round(float(np.median(nullg)),5),
            ctl_null_p95=round(float(thr),5), p_empirical=round(p_emp,4),
            min_detectable_f_minus1=mdf.get(L-1), min_detectable_f_plus1=mdf.get(L+1)))
    r = pd.DataFrame(rows)
    r.to_csv(f"{OUT}/tables/e2_3b_detection_bounds_261104.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(r.to_string(index=False))
    print("peak_rss_gb", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**3)

if __name__ == "__main__":
    main()
