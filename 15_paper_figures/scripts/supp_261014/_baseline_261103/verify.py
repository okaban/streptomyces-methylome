"""Compare every number the three merged figures bake in against the values the
six ORIGINAL generators produce (captured in baseline_*.json by capture.py)."""
import importlib.util, json, os, pathlib, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

SCRIPTS = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
OUT = pathlib.Path(__file__).resolve().parent
os.environ.setdefault("SUPPFIG_OUT", str(OUT / "png_verify"))
sys.path.insert(0, str(SCRIPTS)); sys.path.insert(0, str(SCRIPTS / "supp_261014"))
os.chdir(SCRIPTS)


def load(tag, fname):
    spec = importlib.util.spec_from_file_location("m_" + tag, SCRIPTS / "supp_261014" / fname)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


B = {k: json.load(open(OUT / f"baseline_{k}.json")) for k in ("s2", "s4", "s6", "s7", "s12", "s13")}
rows = []


def cmp(name, base, got, tol=0.0):
    if isinstance(base, (int, float)) and isinstance(got, (int, float)):
        ok = abs(float(base) - float(got)) <= tol
    else:
        ok = base == got
    rows.append((name, repr(base), repr(got), "PASS" if ok else "FAIL"))


# ───────────────── merged figure 1: BGC (old S6 + S7) ───────────────────────
m = load("bgc", "suppfigure_merged_bgc_261103.py")
fig, series = m.build()
for t in ("T1", "T2", "T3"):
    cmp(f"BGC a n4 {t}", B["s6"]["n4"][t], int(series[t]["n4"]))
    cmp(f"BGC b n6 {t}", B["s6"]["n6"][t], int(series[t]["n6"]))
    cmp(f"BGC a bin-max {t}", B["s6"]["m4_bin_max"][t], int(series[t]["m4"].max()))
    cmp(f"BGC b bin-max {t}", B["s6"]["m6_bin_max"][t], int(series[t]["m6"].max()))
    cmp(f"BGC a bin-sum {t}", B["s6"]["m4_bin_sums"][t], int(series[t]["m4"].sum()))
    cmp(f"BGC b bin-sum {t}", B["s6"]["m6_bin_sums"][t], int(series[t]["m6"].sum()))
for k, v in B["s6"]["status_counts"].items():
    cmp(f"BGC c regions {k}", v, int((m.bgc.status == k).sum()))
for key in ("True", "False"):
    kk = key == "True"
    cmp(f"BGC d with_site {key}", B["s7"]["panel_a"][key]["w"], int(m.summ.loc[kk, "with_site"]))
    cmp(f"BGC d n_promoters {key}", B["s7"]["panel_a"][key]["n"], int(m.summ.loc[kk, "n_promoters"]))
    cmp(f"BGC d pct {key}", B["s7"]["panel_a"][key]["pct"], float(m.summ.loc[kk, "pct"]), 1e-9)
for t in ("T1", "T2", "T3"):
    cmp(f"BGC e promoter level {t}", B["s7"]["panel_b"][t]["prom_pct"],
        float(m.panel_b[("level_pooled_pct", m.PROM_KEY)][t]), 1e-9)
    cmp(f"BGC e rest level {t}", B["s7"]["panel_b"][t]["rest_pct"],
        float(m.panel_b[("level_pooled_pct", m.REST_KEY)][t]), 1e-9)
    cmp(f"BGC e n_positions promoter {t}", B["s7"]["panel_b"][t]["n_prom"],
        int(m.panel_b[("n_positions", m.PROM_KEY)][t]))
    cmp(f"BGC e n_positions rest {t}", B["s7"]["panel_b"][t]["n_rest"],
        int(m.panel_b[("n_positions", m.REST_KEY)][t]))
cmp("BGC f n", B["s7"]["panel_c"]["n"], int(len(m.L4)))
cmp("BGC f Spearman rho", B["s7"]["panel_c"]["rho"], float(m.r4), 1e-12)
cmp("BGC f Spearman p", B["s7"]["panel_c"]["p"], float(m.p4), 1e-12)
cmp("BGC f level min", B["s7"]["panel_c"]["ymin"], float(m.L4.level_pooled_pct.min()), 1e-9)
cmp("BGC f level max", B["s7"]["panel_c"]["ymax"], float(m.L4.level_pooled_pct.max()), 1e-9)
for k, v in B["s7"]["panel_c"]["cls_counts"].items():
    cmp(f"BGC f n {k}", v, int((m.L4.cls == k).sum()))
for b, r in zip(B["s7"]["panel_d"], m.dpan.itertuples()):
    cmp(f"BGC g ratio {b['mark']} {b['tp']}", b["ratio"], float(r.ratio), 1e-12)

# ───────── merged figure 2: AAGCCCG 5' (old S12 + S13) ──────────────────────
m = load("aag", "suppfigure_merged_aagcccg5p_261103.py")
for o in m.ORDER:
    cmp(f"AAG a n {o}", B["s12"]["ladder_n"][o],
        int(m.LAD[(m.LAD.definition == o) & (m.LAD.tp == "T1")].n_sites.iloc[0]))
    for t in ("T1", "T2", "T3"):
        cmp(f"AAG a occupancy {o} {t}", B["s12"]["ladder_occ"][o][t],
            float(m.LAD[(m.LAD.definition == o) & (m.LAD.tp == t)].occupancy_pct.iloc[0]), 1e-9)
cmp("AAG a 4mC dashed line", B["s12"]["c4_dashed"], float(m.C4.C4_H_occ.mean()), 1e-9)
for t in ("T1", "T2", "T3"):
    for off, q in B["s12"]["profile_q_at"][t].items():
        got = m.PR[(m.PR.tp == t) & (m.PR.offset == int(off))].q.iloc[0]
        cmp(f"AAG b q offset {off} {t}", q, float(got), 1e-15)
for ch in m.CHANS:
    for b_ in "ACGT":
        cmp(f"AAG c occupancy {ch} {b_}", B["s12"]["ctrl_T1"][ch][b_],
            float(m.CTRL[(m.CTRL.channel == ch) & (m.CTRL.tp == "T1") &
                         (m.CTRL.base_at_minus3 == b_)].occupancy_pct.iloc[0]), 1e-9)
for b in B["s13"]["ladder"]:
    g = m.lad2[(m.lad2.timepoint == b["tp"]) & (m.lad2.tier.astype(str) == b["tier"])].iloc[0]
    cmp(f"AAG d n {b['tp']} tier {b['tier']}", b["n"], int(g.n))
    cmp(f"AAG d called {b['tp']} tier {b['tier']}", b["called"], int(g.called))
    cmp(f"AAG d occupancy {b['tp']} tier {b['tier']}", b["occ"], float(g.occupancy_pct), 1e-9)
for nm, d in (("6mA", m.d6H), ("4mC", m.d4A)):
    key = "d6H" if nm == "6mA" else "d4A"
    t1 = d[d.timepoint == "T1"]
    best = t1.loc[t1.groupby("position").q_BH.idxmin()].sort_values("position")
    got = [{"pos": int(r.position), "base": r.base, "or": float(r.odds_ratio)} for r in best.itertuples()]
    exp = [{"pos": r["pos"], "base": r["base"], "or": r["or"]} for r in B["s13"][key + "_best_T1"]]
    cmp(f"AAG e best base/OR per offset ({nm})", exp, got)
    cmp(f"AAG e offsets surviving 3 tp ({nm})", B["s13"][key + "_survives_pos"],
        sorted(int(p) for p in set(d[d.survives_3tp].position)))
    cmp(f"AAG e tests per timepoint ({nm})", B["s13"][key + "_n_tests"],
        {t: int(d[d.timepoint == t].shape[0]) for t in ("T1", "T2", "T3")})
for k, v in B["s13"]["panel_c_means"].items():
    cmp(f"AAG f mean {k}", v, float(m.ctl[k].mean()), 1e-12)
for b, r in zip(B["s13"]["ctl"], m.ctl.itertuples()):
    cmp(f"AAG f n_mod {b['timepoint']}", b["n_mod"], int(r.n_mod))
    cmp(f"AAG f n_unmod {b['timepoint']}", b["n_unmod"], int(r.n_unmod))

# ───────── merged figure 3: occupancy + k-mer (old S2 + S4) ─────────────────
m = load("occ", "suppfigure_merged_occupancy_261103.py")
for b in B["s2"]["CF"]:
    g = m.CF[(m.CF.key == b["key"]) & (m.CF.tp == b["tp"])].iloc[0]
    cmp(f"OCC n_sites {b['key']} {b['tp']}", b["n_sites"], int(g.n_sites))
    cmp(f"OCC core_pct {b['key']} {b['tp']}", b["core_pct"], float(g.core_pct), 1e-12)
    cmp(f"OCC null_pct {b['key']} {b['tp']}", b["null_pct"], float(g.null_pct), 1e-12)
    cmp(f"OCC binomial p {b['key']} {b['tp']}", b["p_binom"], float(g.p_binom), 1e-12)
for k, v in B["s2"]["hist_max"].items():
    s_, t_ = k.split("|")
    cmp(f"OCC lane max {s_} {t_}", v,
        int(np.histogram(np.array(m.POS[(s_, t_)]),
                         bins=np.arange(0, m.GLEN + 50000, 50000))[0].max()))
for k, v in B["s4"]["aag"].items():
    if k != "group":
        cmp(f"OCC e/f AAGCCCG {k}", v, float(getattr(m.aag, k)), 1e-12)
for k, v in B["s4"]["non"].items():
    if k != "group":
        cmp(f"OCC e/f non-motif {k}", v, float(getattr(m.non, k)), 1e-12)
cmp("OCC f specificity fold", B["s4"]["spec_fold"], float(m.aag.near_rate / m.non.near_rate), 1e-9)
cmp("OCC g 6mA prob histogram", [round(x, 12) for x in B["s4"]["prob_A"]],
    [round(float(x), 12) for x in (m.probs.AAGCCCG_A_6mA / m.probs.AAGCCCG_A_6mA.sum())])
cmp("OCC g 4mC prob histogram", [round(x, 12) for x in B["s4"]["prob_C"]],
    [round(float(x), 12) for x in (m.probs.AAGCCCG_C_4mC / m.probs.AAGCCCG_C_4mC.sum())])

with open(OUT / "verification_261103.tsv", "w") as fh:
    fh.write("quantity\toriginal\tmerged\tresult\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
nfail = sum(1 for r in rows if r[3] == "FAIL")
print(f"checked {len(rows)} quantities, {nfail} FAIL")
for r in rows:
    if r[3] == "FAIL":
        print("FAIL", r[0], r[1], "->", r[2])
