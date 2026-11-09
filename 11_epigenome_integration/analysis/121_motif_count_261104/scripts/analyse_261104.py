#!/usr/bin/env python3
"""121 - does per-gene recognition-sequence COUNT relate to expression?

Primary estimand: partial Spearman of count vs expression controlling
log10(length) and GC. Naive Spearman reported as contrast. Secondary:
quintile(length) x quintile(GC) stratified Spearman (Fisher-z pooled),
Kruskal-Wallis over count groups 0/1/2/3+, and the same at three statistical
units (gene -> transcription unit -> 20 kb tile).
"""
import os, sys, resource, itertools
import numpy as np, pandas as pd
from scipy import stats

BASE = "/Users/okaban/bioinfo/rna-seq"
HERE = f"{BASE}/11_epigenome_integration/analysis/121_motif_count_261104"
FC = f"{BASE}/03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
DE = f"{BASE}/04_deseq2/results/results"
GLEN = 8667507
TPS = ["T1", "T2", "T3"]
TILE = 20000
ALPHA = 0.05


# ---------------------------------------------------------------- statistics
def rankz(x):
    return stats.rankdata(np.asarray(x, float))


def partial_spearman(x, y, Z):
    """Spearman partial correlation: rank-transform then linear partialling."""
    n = len(x)
    k = Z.shape[1]
    if n - 3 - k < 2:
        return np.nan, np.nan, n
    X = np.column_stack([rankz(x), rankz(y)] + [rankz(Z[:, j]) for j in range(k)])
    X = X - X.mean(0)
    C = np.corrcoef(X, rowvar=False)
    try:
        P = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        return np.nan, np.nan, n
    r = -P[0, 1] / np.sqrt(P[0, 0] * P[1, 1])
    df = n - 2 - k
    if abs(r) >= 1:
        return r, 0.0, n
    t = r * np.sqrt(df / (1 - r * r))
    return float(r), float(2 * stats.t.sf(abs(t), df)), n


def min_detectable_r(n, k=0, alpha=ALPHA, power=0.8):
    """Smallest |r| detectable at given power (Fisher z)."""
    d = n - 3 - k
    if d <= 0:
        return np.nan
    z = (stats.norm.isf(alpha / 2) + stats.norm.isf(1 - power)) / np.sqrt(d)
    return float(np.tanh(z))


def stratified_spearman(x, y, strat, min_n=20):
    """Spearman within each stratum, pooled by Fisher z weighted by n-3."""
    zs, ws, ns = [], [], []
    for s in pd.unique(strat):
        m = strat == s
        if m.sum() < min_n:
            continue
        xx, yy = np.asarray(x)[m], np.asarray(y)[m]
        if len(np.unique(xx)) < 2 or len(np.unique(yy)) < 2:
            continue
        r = stats.spearmanr(xx, yy).statistic
        if not np.isfinite(r) or abs(r) >= 1:
            continue
        zs.append(np.arctanh(r)); ws.append(m.sum() - 3); ns.append(m.sum())
    if not zs:
        return np.nan, np.nan, 0, 0
    zs, ws = np.array(zs), np.array(ws, float)
    zbar = float((zs * ws).sum() / ws.sum())
    se = 1 / np.sqrt(ws.sum())
    p = float(2 * stats.norm.sf(abs(zbar) / se))
    return float(np.tanh(zbar)), p, len(zs), int(sum(ns))


def kw_groups(count, y):
    g = np.clip(np.asarray(count), 0, 3)
    grps = [np.asarray(y)[g == i] for i in range(4)]
    grps = [a for a in grps if len(a) >= 5]
    if len(grps) < 2:
        return np.nan, np.nan, []
    H, p = stats.kruskal(*grps)
    return float(H), float(p), [len(a) for a in grps]


# ---------------------------------------------------------------- inputs
def load_tpm():
    fc = pd.read_csv(FC, sep="\t", comment="#")
    libcols = [c for c in fc.columns if c.endswith(".bam")]
    lab = [os.path.basename(c).split(".")[0].replace("M145_", "") for c in libcols]
    cnt = fc[libcols].to_numpy(float)
    rpk = cnt / fc["Length"].to_numpy(float)[:, None]
    tpm = rpk / rpk.sum(0) * 1e6
    t = pd.DataFrame(tpm, columns=lab); t["gene_id"] = fc["Geneid"].values
    for tp, pre in zip(TPS, ["1_", "2_", "3_"]):
        cs = [c for c in lab if c.startswith(pre)]
        assert len(cs) == 3, cs
        t[f"tpm_{tp}"] = t[cs].mean(1)
    return t[["gene_id"] + [f"tpm_{tp}" for tp in TPS]]


def load_lfc():
    out = {}
    for c in ["2_vs_1", "3_vs_1", "3_vs_2"]:
        d = pd.read_csv(f"{DE}/DESeq2_M145_{c}.tsv", sep="\t")
        out[c] = d[["gene_id", "log2FoldChange"]].rename(
            columns={"log2FoldChange": f"lfc_{c}"})
    return out


def qbin(x, q=5):
    return pd.qcut(pd.Series(x).rank(method="first"), q, labels=False)


# ---------------------------------------------------------------- units
def build_gene_table():
    g = pd.read_csv(f"{HERE}/tables/per_gene_counts_261104.tsv", sep="\t")
    g = g.merge(load_tpm(), on="gene_id", how="left")
    for c, d in load_lfc().items():
        g = g.merge(d, on="gene_id", how="left")
    for tp in TPS:
        g[f"logtpm_{tp}"] = np.log2(g[f"tpm_{tp}"] + 1)
    g["loglen"] = np.log10(g.glen)
    g["mid"] = (g.b0 + g.b1) // 2
    g["dist_centre"] = (g["mid"] - GLEN / 2).abs()
    return g


def build_tu_table(g):
    """Transcription units: same strand, intergenic gap <= 50 bp."""
    s = g.sort_values(["strand", "start"]).reset_index(drop=True)
    tu, cur, prev_end, prev_strand = [], 0, None, None
    for _, r in s.iterrows():
        if prev_strand is None or r.strand != prev_strand or (r.b0 - prev_end) > 50:
            cur += 1
        tu.append(cur)
        prev_end = r.b1 if prev_strand != r.strand else max(prev_end, r.b1)
        prev_strand = r.strand
    s["tu"] = tu
    agg = {}
    cntcols = [c for c in g.columns if c.startswith(("n_", "nmod_"))]
    for c in cntcols:
        agg[c] = (c, "sum")
    rows = []
    for tuid, d in s.groupby("tu"):
        w = d.glen.to_numpy(float)
        rec = {"unit": tuid, "glen": w.sum(), "n_genes": len(d),
               "gene_gc": float((d.gene_gc * w).sum() / w.sum()),
               "prom_gc": float(d.prom_gc.mean()),
               "mid": float((d["mid"] * w).sum() / w.sum())}
        for c in cntcols:
            rec[c] = int(d[c].sum())
        for tp in TPS:
            rec[f"tpm_{tp}"] = float((d[f"tpm_{tp}"] * w).sum() / w.sum())
        for c in ["2_vs_1", "3_vs_1", "3_vs_2"]:
            v = d[f"lfc_{c}"].to_numpy(float); m = np.isfinite(v)
            rec[f"lfc_{c}"] = float((v[m] * w[m]).sum() / w[m].sum()) if m.any() else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    for tp in TPS:
        t[f"logtpm_{tp}"] = np.log2(t[f"tpm_{tp}"] + 1)
    t["loglen"] = np.log10(t.glen)
    t["dist_centre"] = (t["mid"] - GLEN / 2).abs()
    return t


def build_tile_table(g):
    cntcols = [c for c in g.columns if c.startswith(("n_", "nmod_"))]
    mid = ((g.b0 + g.b1) // 2).to_numpy()
    g2 = g.copy(); g2["tile"] = mid // TILE
    rows = []
    for tid, d in g2.groupby("tile"):
        w = d.glen.to_numpy(float)
        rec = {"unit": int(tid), "glen": w.sum(), "n_genes": len(d),
               "gene_gc": float((d.gene_gc * w).sum() / w.sum()),
               "prom_gc": float(d.prom_gc.mean()),
               "mid": (int(tid) + 0.5) * TILE}
        for c in cntcols:
            rec[c] = int(d[c].sum())
        for tp in TPS:
            rec[f"tpm_{tp}"] = float((d[f"tpm_{tp}"] * w).sum() / w.sum())
        for c in ["2_vs_1", "3_vs_1", "3_vs_2"]:
            v = d[f"lfc_{c}"].to_numpy(float); m = np.isfinite(v)
            rec[f"lfc_{c}"] = float((v[m] * w[m]).sum() / w[m].sum()) if m.any() else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    for tp in TPS:
        t[f"logtpm_{tp}"] = np.log2(t[f"tpm_{tp}"] + 1)
    t["loglen"] = np.log10(t.glen)
    t["dist_centre"] = (t["mid"] - GLEN / 2).abs()
    return t


# ---------------------------------------------------------------- driver
MOTIFS = ["TGGCCGGC_4mC", "AAGCCCG_4mC", "AAGCCCG_6mA", "ANY"]


def predictors(df):
    out = []
    for m in MOTIFS:
        for reg in ("gene", "prom"):
            c = f"n_{m}_{reg}"
            if c in df.columns:
                out.append((c, "possible", m, reg))
            c2 = f"nmod_{m}_{reg}_T1"
            if c2 in df.columns:
                out.append((c2, "modified_T1", m, reg))
    return out


def run_unit(df, unit, responses, rows):
    preds = predictors(df)
    lq = qbin(df.loglen); 
    for pcol, kind, motif, reg in preds:
        x = df[pcol].to_numpy(float)
        gccol = "prom_gc" if reg == "prom" else "gene_gc"
        gq = qbin(df[gccol])
        strat = (lq.astype(str) + "_" + gq.astype(str)).to_numpy()
        Z = np.column_stack([df["loglen"].to_numpy(float), df[gccol].to_numpy(float)])
        Z3 = np.column_stack([Z, df["dist_centre"].to_numpy(float)])
        for rcol in responses:
            y = df[rcol].to_numpy(float)
            ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(1)
            if ok.sum() < 50 or len(np.unique(x[ok])) < 2:
                continue
            nr = stats.spearmanr(x[ok], y[ok])
            ar, ap, n = partial_spearman(x[ok], y[ok], Z[ok])
            a3r, a3p, _ = partial_spearman(x[ok], y[ok], Z3[ok])
            sr, sp, nstr, nused = stratified_spearman(x[ok], y[ok], strat[ok])
            H, kp, sizes = kw_groups(x[ok], y[ok])
            rows.append(dict(
                unit=unit, n=int(ok.sum()), predictor=pcol, count_kind=kind,
                motif=motif, region=reg, response=rcol,
                naive_rho=round(float(nr.statistic), 4), naive_p=float(nr.pvalue),
                adj_rho=round(ar, 4) if np.isfinite(ar) else np.nan, adj_p=ap,
                adj3_rho=round(a3r, 4) if np.isfinite(a3r) else np.nan, adj3_p=a3p,
                strat_rho=round(sr, 4) if np.isfinite(sr) else np.nan, strat_p=sp,
                n_strata=nstr, n_in_strata=nused,
                kw_H=round(H, 3) if np.isfinite(H) else np.nan, kw_p=kp,
                kw_group_sizes=",".join(map(str, sizes)),
                mdr_naive=round(min_detectable_r(int(ok.sum()), 0), 4),
                mdr_adj=round(min_detectable_r(int(ok.sum()), 2), 4),
                mdr_adj3=round(min_detectable_r(int(ok.sum()), 3), 4),
                mdr_strat=round(min_detectable_r(nused, 0), 4) if nused else np.nan))


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main():
    g = build_gene_table()
    sys.stderr.write(f"[gene] {len(g)} genes, TPM finite {int(np.isfinite(g.tpm_T1).sum())}\n")
    tu = build_tu_table(g); tile = build_tile_table(g)
    sys.stderr.write(f"[tu] {len(tu)} units  [tile] {len(tile)} tiles\n")
    tu.to_csv(f"{HERE}/tables/per_tu_261104.tsv", sep="\t", index=False)
    tile.to_csv(f"{HERE}/tables/per_tile_261104.tsv", sep="\t", index=False)

    resp_expr = [f"logtpm_{tp}" for tp in TPS]
    resp_lfc = [f"lfc_{c}" for c in ["2_vs_1", "3_vs_1", "3_vs_2"]]
    for d in (g, tu, tile):
        for c in ["2_vs_1", "3_vs_1", "3_vs_2"]:
            d[f"abslfc_{c}"] = d[f"lfc_{c}"].abs()
    resp_abs = [f"abslfc_{c}" for c in ["2_vs_1", "3_vs_1", "3_vs_2"]]

    rows = []
    for df, unit in ((g, "gene"), (tu, "transcription_unit"), (tile, "tile_20kb")):
        run_unit(df, unit, resp_expr + resp_lfc + resp_abs, rows)
    R = pd.DataFrame(rows)
    R["response_class"] = np.where(R.response.str.startswith("logtpm"), "expression_level",
                         np.where(R.response.str.startswith("abslfc"), "expression_change_abs",
                                  "expression_change_signed"))
    for col, new in (("naive_p", "naive_q"), ("adj_p", "adj_q"), ("adj3_p", "adj3_q"),
                     ("strat_p", "strat_q"), ("kw_p", "kw_q")):
        R[new] = np.nan
        for cls, idx in R.groupby("response_class").groups.items():
            v = R.loc[idx, col].to_numpy(float)
            m = np.isfinite(v)
            qq = np.full(len(v), np.nan); qq[m] = bh(v[m])
            R.loc[idx, new] = qq
    R = R.sort_values(["response_class", "unit", "motif", "region", "response"])
    R.to_csv(f"{HERE}/tables/correlation_results_261104.tsv", sep="\t", index=False)
    sys.stderr.write(f"[write] correlation_results_261104.tsv ({len(R)} tests)\n")

    # power summary
    pw = []
    for unit, n in (("gene", len(g)), ("transcription_unit", len(tu)), ("tile_20kb", len(tile))):
        for k, lab in ((0, "naive"), (2, "partial_len_GC"), (3, "partial_len_GC_pos")):
            pw.append(dict(unit=unit, n=n, k_covariates=k, analysis=lab,
                           min_detectable_rho_power80=round(min_detectable_r(n, k), 4),
                           min_detectable_rho_power50=round(min_detectable_r(n, k, power=0.5), 4)))
    sub = R[(R.unit == "gene") & (R.response == "logtpm_T1")]
    for _, r in sub.iterrows():
        pw.append(dict(unit="gene_stratified", n=int(r.n_in_strata), k_covariates=0,
                       analysis=f"stratified_{r.predictor}",
                       min_detectable_rho_power80=r.mdr_strat,
                       min_detectable_rho_power50=np.nan))
    pd.DataFrame(pw).to_csv(f"{HERE}/tables/power_261104.tsv", sep="\t", index=False)

    # headline pick: strongest naive |rho| among gene-level expression-level tests
    prim = R[(R.unit == "gene") & (R.response_class == "expression_level")].copy()
    prim["abs_naive"] = prim.naive_rho.abs()
    top = prim.sort_values("abs_naive", ascending=False).iloc[0]
    prim["abs_adj"] = prim.adj_rho.abs()
    topadj = prim.sort_values("abs_adj", ascending=False).iloc[0]
    hl = pd.DataFrame([
        dict(which="strongest_naive", **{k: top[k] for k in
             ["predictor", "response", "n", "naive_rho", "naive_p", "naive_q",
              "adj_rho", "adj_p", "adj_q", "adj3_rho", "adj3_p", "adj3_q",
              "strat_rho", "strat_p", "kw_H", "kw_p", "mdr_naive", "mdr_adj"]}),
        dict(which="strongest_adjusted", **{k: topadj[k] for k in
             ["predictor", "response", "n", "naive_rho", "naive_p", "naive_q",
              "adj_rho", "adj_p", "adj_q", "adj3_rho", "adj3_p", "adj3_q",
              "strat_rho", "strat_p", "kw_H", "kw_p", "mdr_naive", "mdr_adj"]})])
    hl.to_csv(f"{HERE}/tables/headline_261104.tsv", sep="\t", index=False)
    print(hl.to_string(index=False))
    # confounder sanity: length/GC vs count and vs expression
    san = []
    for c in ["n_ANY_gene", "n_TGGCCGGC_4mC_gene", "n_AAGCCCG_4mC_gene"]:
        for v in ["loglen", "gene_gc"]:
            r = stats.spearmanr(g[c], g[v])
            san.append(dict(pair=f"{c}~{v}", rho=round(r.statistic, 4), p=r.pvalue))
    for v in ["loglen", "gene_gc", "dist_centre"]:
        r = stats.spearmanr(g[v], g["logtpm_T1"])
        san.append(dict(pair=f"{v}~logtpm_T1", rho=round(r.statistic, 4), p=r.pvalue))
    til = pd.read_csv(f"{HERE}/tables/per_tile_261104.tsv", sep="\t")
    for v in ["n_AAGCCCG_4mC_gene", "n_ANY_gene", "logtpm_T1", "lfc_3_vs_1", "gene_gc"]:
        r = stats.spearmanr(til[v], til["dist_centre"])
        san.append(dict(pair=f"TILE {v}~dist_centre", rho=round(r.statistic, 4), p=r.pvalue))
    pd.DataFrame(san).to_csv(f"{HERE}/tables/confounder_structure_261104.tsv",
                             sep="\t", index=False)
    print(pd.DataFrame(san).to_string(index=False))
    sys.stderr.write("[peak RSS GB] %.2f\n" %
                     (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3))


if __name__ == "__main__":
    main()
