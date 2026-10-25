"""支持された 3 所見（F1 convergent gap / F2 センス鎖 / F3 プロモーター解離）の再計算。
実行: cd scripts && PYTHONPATH=. python3 run_supported_findings.py
出力は ../tables/ に書く。正典規則は core.canonical / core.measurable が実装する。
初回はパイルアップキャッシュの読み直しに約 2 分（9 ライブラリ × 2 修飾種 × 3 時点）。
"""
import warnings; warnings.filterwarnings("ignore")
import core, numpy as np, pandas as pd, pathlib, re
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable
from scipy import stats as st

HERE = pathlib.Path(__file__).resolve().parent
TAB = HERE.parent / "tables"; TAB.mkdir(exist_ok=True)
G = core.genome()
aa = core.aagcccg_instances().copy()
tg = core.tggccggc_positions().copy()
gs = core.predict_tus(max_gap=50).sort_values("gstart").reset_index(drop=True)
G0, G1 = gs.gstart.values - 1, gs.gend.values
STR, TUID, TUPOS, TUSZ = gs.gstrand.values, gs.tu_id.values, gs.tu_pos.values, gs.tu_size.values

# ---- vectorised local GC (±250 bp)
_g = np.frombuffer(G.encode(), dtype="S1")
_cs = np.concatenate([[0], np.cumsum(((_g == b"G") | (_g == b"C")).astype(np.int32))])
def gc_vec(pos, half=250):
    L = len(_g); a = np.clip(pos - half, 0, L); b = np.clip(pos + half + 1, 0, L)
    return (_cs[b] - _cs[a]) / (b - a)

# ---- canonical site / measurable caches
CP = TAB / "canonical_sites_cache.npz"
if not CP.exists():
    store = {}
    for tp in ("T1", "T2", "T3"):
        for mod in ("6mA", "4mC"):
            c = core.canonical(tp, mod); m = core.measurable(tp, mod)
            for tag, d in (("hc", c), ("me", m)):
                store[f"{tag}_{tp}_{mod}"] = np.rec.fromarrays(
                    [d.start.values.astype(np.int32), d.plus.values.astype(bool),
                     d.lvl.values, d.depth.values.astype(np.int32)],
                    names="start,plus,lvl,depth")
            core.canonical.cache_clear(); core.measurable.cache_clear()
    np.savez_compressed(CP, **store)
SITES = dict(np.load(CP, allow_pickle=False))

def annotate(pos, strand_plus=None):
    pos = np.asarray(pos)
    j = np.searchsorted(G0, pos, side="right") - 1
    host = np.full(len(pos), -1, np.int64)
    for k in range(len(pos)):
        jj = j[k]
        while jj >= 0 and pos[k] - G0[jj] < 50_000:
            if G0[jj] <= pos[k] < G1[jj]:
                host[k] = jj; break
            jj -= 1
    out = pd.DataFrame({"pos": pos, "host": host})
    out["intragenic"] = host >= 0
    hc = np.clip(host, 0, None)
    out["gene_id"] = np.where(host >= 0, gs.gene_id.values[hc], "")
    out["gstrand"] = np.where(host >= 0, STR[hc], "")
    out["tu_pos"] = np.where(host >= 0, TUPOS[hc], 0)
    out["tu_size"] = np.where(host >= 0, TUSZ[hc], 0)
    jl = np.searchsorted(G1, pos, side="right"); left = jl - 1
    right = np.searchsorted(G0, pos, side="right")
    lv = np.where((left >= 0) & (~out.intragenic.values), left, -1)
    rv = np.where((right < len(G0)) & (~out.intragenic.values), right, -1)
    ls = np.where(lv >= 0, STR[np.clip(lv, 0, None)], "")
    rs = np.where(rv >= 0, STR[np.clip(rv, 0, None)], "")
    cls = np.where((lv < 0) | (rv < 0), "edge",
          np.where((ls == "+") & (rs == "-"), "convergent",
          np.where((ls == "-") & (rs == "+"), "divergent", "tandem")))
    out["pair_class"] = np.where(out.intragenic, "intragenic", cls)
    tuf = gs[gs.tu_pos == 1]
    tu5 = np.where(tuf.gstrand.values == "+", tuf.gstart.values - 1, tuf.gend.values - 1)
    o = np.argsort(tu5); t5, t5s = tu5[o], tuf.gstrand.values[o]
    i = np.searchsorted(t5, pos)
    cand = np.stack([np.clip(i - 1, 0, len(t5) - 1), np.clip(i, 0, len(t5) - 1)])
    d = pos[None, :] - t5[cand]; pick = np.argmin(np.abs(d), axis=0)
    nd = d[pick, np.arange(len(pos))]; nstr = t5s[cand[pick, np.arange(len(pos))]]
    out["d_tu5"] = np.where(nstr == "+", nd, -nd)
    out["replichore"] = np.where(pos >= core.ORIC, "right", "left")
    if strand_plus is not None:
        out["mod_plus"] = np.asarray(strand_plus)
    return out

def tu_class(an):
    c = np.full(len(an), "other_intergenic", dtype=object)
    intra = an.intragenic.values; first = an.tu_pos.values == 1; multi = an.tu_size.values > 1
    c[intra & first & multi] = "TUfirst_body"
    c[intra & (~first) & multi] = "TUinternal_body"
    c[intra & (an.tu_size.values == 1)] = "single_gene_body"
    up = (~intra) & (an.d_tu5.values >= -300) & (an.d_tu5.values < 0)
    c[up] = "TU_upstream300"
    c[(~intra) & (an.pair_class.values == "tandem") & ~up] = "tandem_gap"
    c[(~intra) & (an.pair_class.values == "divergent") & ~up] = "divergent_gap"
    c[(~intra) & (an.pair_class.values == "convergent")] = "convergent_gap"
    return c

def mh(y, expo, gcv, dep, nq=10, dq_n=5):
    gq = pd.qcut(gcv, nq, labels=False, duplicates="drop")
    dq = pd.qcut(dep.astype(float), dq_n, labels=False, duplicates="drop")
    tabs = []
    for i in range(int(np.nanmax(gq)) + 1):
        for j in range(int(np.nanmax(dq)) + 1):
            s = (gq == i) & (dq == j)
            if s.sum() < 50: continue
            t = np.array([[int((y & expo & s).sum()), int((~y & expo & s).sum())],
                          [int((y & ~expo & s).sum()), int((~y & ~expo & s).sum())]])
            if t[:, 0].sum() == 0 or t[0].sum() == 0 or t[1].sum() == 0: continue
            tabs.append(t)
    stx = StratifiedTable(np.dstack(tabs)); lo, hi = stx.oddsratio_pooled_confint()
    return dict(n_strata=len(tabs), OR_MH=round(stx.oddsratio_pooled, 3),
                ci=f"[{lo:.2f},{hi:.2f}]", p_MH=stx.test_null_odds().pvalue,
                p_homogeneity=stx.test_equal_odds().pvalue)

# =============================== F1: convergent gaps ===============================
CLS = {}
for tp in ("T1", "T2", "T3"):
    for mod in ("6mA", "4mC"):
        me = SITES[f"me_{tp}_{mod}"]
        CLS[(tp, mod)] = tu_class(annotate(me["start"], me["plus"]))
rows = []
cmps = [("TU_upstream300", "tandem_gap"), ("TU_upstream300", "gene_body"),
        ("TUfirst_body", "TUinternal_body"), ("convergent_gap", "gene_body"),
        ("divergent_gap", "gene_body")]
for tp in ("T1", "T2", "T3"):
    for mod in ("6mA", "4mC"):
        me = SITES[f"me_{tp}_{mod}"]; cl = CLS[(tp, mod)]
        grp = np.where(np.isin(cl, ["TUfirst_body", "TUinternal_body", "single_gene_body"]),
                       "gene_body", cl)
        y = me["lvl"] >= 0.5; gcv = gc_vec(me["start"]); dep = me["depth"]
        for A, B in cmps:
            cA = (cl == A) if A.endswith("_body") else (grp == A)
            cB = (cl == B) if B.endswith("_body") and B != "gene_body" else (grp == B)
            s = cA | cB
            if cA.sum() < 30 or cB.sum() < 30: continue
            r = mh(y[s], cA[s], gcv[s], dep[s])
            r.update(tp=tp, mod=mod, exposure=A, reference=B, n_exp=int(cA.sum()),
                     hc_exp=int((y & cA).sum()), n_ref=int(cB.sum()), hc_ref=int((y & cB).sum()))
            rows.append(r)
TU = pd.DataFrame(rows); TU["q_MH"] = core.bh(TU.p_MH)
TU.to_csv(TAB / "IH01_IH07_IH09_tu_regions.tsv", sep="\t", index=False)

# --- GC window sensitivity: 6mA call rate depends strongly on LOCAL (+/-25 bp) GC.
# Stratifying on +/-250 bp GC alone over-states the convergent-gap effect; the
# +/-25 bp version is the one that matches bio261016 (2026-10-16) and is reported.
rows = []
for tp in ("T1", "T2", "T3"):
    me = SITES[f"me_{tp}_6mA"]; cl = CLS[(tp, "6mA")]
    y = me["lvl"] >= 0.5; dep = me["depth"]; g25 = gc_vec(me["start"], half=25)
    grp = np.where(np.isin(cl, ["TUfirst_body", "TUinternal_body", "single_gene_body"]),
                   "gene_body", cl)
    inter = np.isin(cl, ["TU_upstream300", "tandem_gap", "divergent_gap", "convergent_gap"])
    for B, sel in (("gene_body", (grp == "convergent_gap") | (grp == "gene_body")),
                   ("all_other_intergenic", inter)):
        r = mh(y[sel], (cl == "convergent_gap")[sel], g25[sel], dep[sel])
        r.update(tp=tp, reference=B, gc_window="+/-25bp"); rows.append(r)
CG25 = pd.DataFrame(rows); CG25["q_MH"] = core.bh(CG25.p_MH)
CG25.to_csv(TAB / "IH09_convergent_gc25_retest.tsv", sep="\t", index=False)
print("\nF1 convergent gap, +/-25 bp GC stratification:")
print(CG25[["tp", "reference", "OR_MH", "ci", "q_MH"]].to_string(index=False))

rows = []
for tp in ("T1", "T2", "T3"):
    for mod in ("6mA", "4mC"):
        me = SITES[f"me_{tp}_{mod}"]; cl = CLS[(tp, mod)]
        y = me["lvl"] >= 0.5; gcv = gc_vec(me["start"]); dep = me["depth"]
        inter = np.isin(cl, ["TU_upstream300", "tandem_gap", "divergent_gap", "convergent_gap"])
        for B in ["tandem_gap", "TU_upstream300"]:
            s = (cl == "convergent_gap") | (cl == B)
            r = mh(y[s], (cl == "convergent_gap")[s], gcv[s], dep[s], nq=8, dq_n=4)
            r.update(tp=tp, mod=mod, exposure="convergent_gap", reference=B); rows.append(r)
        r = mh(y[inter], (cl == "convergent_gap")[inter], gcv[inter], dep[inter], nq=8, dq_n=4)
        r.update(tp=tp, mod=mod, exposure="convergent_gap", reference="all_other_intergenic")
        rows.append(r)
CGS = pd.DataFrame(rows); CGS["q_MH"] = core.bh(CGS.p_MH)
CGS.to_csv(TAB / "IH09_convergent_specificity.tsv", sep="\t", index=False)

# =============================== F2: sense-strand orientation ======================
an = annotate(aa.anchor.values, aa.plus.values)
aa2 = pd.concat([aa.reset_index(drop=True), an.drop(columns=["pos"]).reset_index(drop=True)], axis=1)
aa2["mod_sense"] = np.where(aa2.intragenic, aa2.plus == aa2.gstrand.eq("+"), np.nan)
aa2["mod_leading"] = aa2.plus == aa2.replichore.eq("right")
aa2["gc500"] = gc_vec(aa2.anchor.values)
aa2["det_H"] = aa2.up6.str[-1].isin(list("ACT"))
aa2["det_tier"] = np.select(
    [aa2.det_H & aa2.m3.eq("G") & aa2.m5.eq("C"), aa2.det_H & aa2.m3.eq("G"), aa2.det_H],
    ["H_m3G_m5C", "H_m3G", "H_only"], default="non_H")
def flank(r, lo, hi):
    return G[r.start + lo:r.start + hi] if r.plus else core.rc(G[r.start + 6 - hi + 1:r.start + 6 - lo + 1])
aa2["ctx_up"] = [flank(r, -6, 0) for r in aa2.itertuples()]
aa2["ctx_dn"] = [flank(r, 7, 10) for r in aa2.itertuples()]
for i, p in enumerate(range(-6, 0)): aa2[f"u{abs(p)}"] = aa2.ctx_up.str[i]
for i, p in enumerate(range(1, 4)):  aa2[f"d{p}"] = aa2.ctx_dn.str[i]
CTX = " + ".join([f"C(u{k})" for k in (6, 5, 4, 3, 2, 1)] + [f"C(d{k})" for k in (1, 2, 3)])
gn_i = gs.set_index("gene_id")
aa2["log_tpm_host"] = aa2.gene_id.map(np.log2(gn_i.tpm_mean + 1))
for tp in ("T1", "T2", "T3"):
    for tag, mod, posc in (("6", "6mA", "pos_A1"), ("4", "4mC", "pos_C4")):
        hc = SITES[f"hc_{tp}_{mod}"]; me = SITES[f"me_{tp}_{mod}"]
        hs = set(zip(hc["start"], hc["plus"])); dm = dict(zip(zip(me["start"], me["plus"]), me["depth"]))
        key = list(zip(aa2[posc], aa2.plus))
        aa2[f"meas{tag}_{tp}"] = [k in dm for k in key]
        aa2[f"mod{tag}_{tp}"] = [k in hs for k in key]
        aa2[f"depth{tag}_{tp}"] = [dm.get(k, np.nan) for k in key]
off = np.where(aa2.gstrand.values == "+",
               aa2.pos_A1.values + 1 - np.where(aa2.host.values >= 0, gs.gstart.values[np.clip(aa2.host.values, 0, None)], np.nan),
               np.where(aa2.host.values >= 0, gs.gend.values[np.clip(aa2.host.values, 0, None)], np.nan) - (aa2.pos_A1.values + 1))
aa2["codon_phase"] = np.where(aa2.host.values >= 0, np.mod(off, 3), np.nan)
aa2.to_csv(TAB / "aagcccg_instances_context.tsv", sep="\t", index=False)

rows = []
for tp in ("T1", "T2", "T3"):
    d = aa2[aa2[f"meas6_{tp}"] & aa2.intragenic].dropna(subset=["log_tpm_host"]).copy()
    d["y"] = d[f"mod6_{tp}"].astype(int); d["ms"] = d.mod_sense.astype(int)
    d["ml"] = d.mod_leading.astype(int); d["pl"] = d.plus.astype(int)
    d["in_core"] = d.anchor.between(*core.CORE).astype(int); d["ldep"] = np.log(d[f"depth6_{tp}"])
    m = smf.logit("y ~ ms + ml + pl + gc500 + ldep + in_core + C(codon_phase) + " + CTX,
                  data=d).fit(disp=0, maxiter=400)
    ci = m.conf_int()
    rows.append(dict(tp=tp, n=len(d), OR_sense=round(np.exp(m.params.ms), 3),
                     ci_sense=f"[{np.exp(ci.loc['ms',0]):.2f},{np.exp(ci.loc['ms',1]):.2f}]",
                     p_sense=m.pvalues.ms, OR_plus=round(np.exp(m.params.pl), 3), p_plus=m.pvalues.pl,
                     OR_leading=round(np.exp(m.params.ml), 3), p_leading=m.pvalues.ml))
pd.DataFrame(rows).to_csv(TAB / "IH10_sense_orientation_adjusted.tsv", sep="\t", index=False)

# controls: 4mC at the same instances, and TGGCCGGC
out = []
for tp in ("T1", "T2", "T3"):
    d = aa2[aa2[f"meas4_{tp}"] & aa2.intragenic].copy()
    d["y"] = d[f"mod4_{tp}"].astype(int); d["ms"] = d.mod_sense.astype(int)
    d["ml"] = d.mod_leading.astype(int); d["in_core"] = d.anchor.between(*core.CORE).astype(int)
    d["ldep"] = np.log(d[f"depth4_{tp}"]); d = d.dropna(subset=["ldep"])
    m = smf.logit("y ~ ms + ml + gc500 + ldep + in_core + " + CTX, data=d).fit(disp=0, maxiter=200)
    out.append(dict(target="AAGCCCG C5 (4mC)", tp=tp, n=len(d),
                    occ_sense=round(100 * d[d.ms == 1].y.mean(), 2),
                    occ_anti=round(100 * d[d.ms == 0].y.mean(), 2),
                    OR_ms=round(np.exp(m.params.ms), 3), p_ms=m.pvalues.ms))
at = annotate(tg.anchor.values, tg.plus.values)
tg2 = pd.concat([tg.reset_index(drop=True), at.drop(columns=["pos"]).reset_index(drop=True)], axis=1)
tg2["mod_sense"] = np.where(tg2.intragenic, tg2.plus == tg2.gstrand.eq("+"), np.nan)
tg2["mod_leading"] = tg2.plus == tg2.replichore.eq("right")
tg2["gc500"] = gc_vec(tg2.anchor.values)
for tp in ("T1", "T2", "T3"):
    hc = SITES[f"hc_{tp}_4mC"]; me = SITES[f"me_{tp}_4mC"]
    hs = set(zip(hc["start"], hc["plus"])); dm = dict(zip(zip(me["start"], me["plus"]), me["depth"]))
    d = tg2[tg2.intragenic].copy()
    d["y"] = [int(k in hs) for k in zip(d.start, d.plus)]
    d["ldep"] = np.log([dm.get(k, np.nan) for k in zip(d.start, d.plus)])
    d["ms"] = d.mod_sense.astype(int); d["ml"] = d.mod_leading.astype(int)
    d["in_core"] = d.anchor.between(*core.CORE).astype(int); d = d.dropna(subset=["ldep"])
    m = smf.logit("y ~ ms + ml + gc500 + ldep + in_core", data=d).fit(disp=0, maxiter=200)
    out.append(dict(target="TGGCCGGC (4mC)", tp=tp, n=len(d),
                    occ_sense=round(100 * d[d.ms == 1].y.mean(), 2),
                    occ_anti=round(100 * d[d.ms == 0].y.mean(), 2),
                    OR_ms=round(np.exp(m.params.ms), 3), p_ms=m.pvalues.ms))
pd.DataFrame(out).to_csv(TAB / "IH10_controls_4mC_and_TGGCCGGC.tsv", sep="\t", index=False)

# =============================== F3: promoter dissociation =========================
t = pd.read_csv(core.TSS_CSV)
sel = t[t.category.astype(str).str.upper().isin(["P", "S"]) & t.gene_id.notna()
        & (t.gene_id.astype(str).str.strip().str.len() > 0)]
o = np.argsort(sel.tss_position.astype(int).values)
tssP = sel.tss_position.astype(int).values[o]; tssS = sel.strand.values[o]
def prox(pos, lo=-300, hi=50):
    i = np.searchsorted(tssP, pos)
    cand = np.stack([np.clip(i - 1, 0, len(tssP) - 1), np.clip(i, 0, len(tssP) - 1)])
    d = pos[None, :] - tssP[cand]; pick = np.argmin(np.abs(d), axis=0)
    nd = d[pick, np.arange(len(pos))]; ns = tssS[cand[pick, np.arange(len(pos))]]
    sd = np.where(ns == "+", nd, -nd)
    return sd, (sd >= lo) & (sd <= hi)
rows = []
for tp in ("T1", "T2", "T3"):
    for mod in ("6mA", "4mC"):
        me = SITES[f"me_{tp}_{mod}"]
        sd, ip = prox(me["start"]); y = me["lvl"] >= 0.5
        a, b = y[ip], y[~ip]
        orr, pv = st.fisher_exact([[int(a.sum()), int((~a).sum())], [int(b.sum()), int((~b).sum())]])
        r = mh(y, ip, gc_vec(me["start"]), me["depth"])
        r.update(tp=tp, mod=mod, n_prox=int(ip.sum()), hc_prox=int(a.sum()),
                 rate_prox_per100k=round(1e5 * a.mean(), 2), rate_far_per100k=round(1e5 * b.mean(), 2),
                 OR_crude=round(orr, 3), p_crude=pv,
                 mean_lvl_prox_pct=round(100 * np.nanmean(me["lvl"][ip]), 3),
                 mean_lvl_far_pct=round(100 * np.nanmean(me["lvl"][~ip]), 3))
        rows.append(r)
P = pd.DataFrame(rows); P["q_crude"] = core.bh(P.p_crude); P["q_MH"] = core.bh(P.p_MH)
P.to_csv(TAB / "IH19_promoter_dissociation.tsv", sep="\t", index=False)
print(P[["tp", "mod", "OR_crude", "OR_MH", "ci", "q_crude", "q_MH"]].to_string(index=False))
