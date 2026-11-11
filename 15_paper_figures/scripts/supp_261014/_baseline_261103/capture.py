"""Run each original generator and dump the numbers it bakes into its panels.
Used only to verify that the merged rebuilds change nothing."""
import json, os, runpy, sys, pathlib

SCRIPTS = pathlib.Path("/Users/okaban/bioinfo/rna-seq/15_paper_figures/scripts")
OUT = pathlib.Path(__file__).resolve().parent
os.environ["SUPPFIG_OUT"] = str(OUT / "png")
(OUT / "png").mkdir(parents=True, exist_ok=True)
import matplotlib
matplotlib.use("Agg")
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "supp_261014"))
os.chdir(SCRIPTS)

which = sys.argv[1]
res = {}

if which == "s6":
    g = runpy.run_path(str(SCRIPTS / "supp_261014/suppfigure13_bgc_by_mark.py"), run_name="__main__")
    s = g["series"]
    res["n4"] = {t: int(s[t]["n4"]) for t in s}
    res["n6"] = {t: int(s[t]["n6"]) for t in s}
    res["m4_bin_sums"] = {t: int(s[t]["m4"].sum()) for t in s}
    res["m6_bin_sums"] = {t: int(s[t]["m6"].sum()) for t in s}
    res["m4_bin_max"] = {t: int(s[t]["m4"].max()) for t in s}
    res["m6_bin_max"] = {t: int(s[t]["m6"].max()) for t in s}
    bgc = g["bgc"]
    res["status_counts"] = {k: int(v) for k, v in bgc.status.value_counts().items()}

elif which == "s7":
    g = runpy.run_path(str(SCRIPTS / "supp_261014/suppfigure14_bgc_relationships.py"), run_name="__main__")
    summ, pb, L4, dpan = g["summ"], g["panel_b"], g["L4"], g["dpan"]
    res["panel_a"] = {str(k): dict(n=int(summ.loc[k, "n_promoters"]), w=int(summ.loc[k, "with_site"]),
                                   pct=float(summ.loc[k, "pct"])) for k in (True, False)}
    PK, RK = g["PROM_KEY"], g["REST_KEY"]
    res["panel_b"] = {t: {"prom_pct": float(pb[("level_pooled_pct", PK)][t]),
                          "rest_pct": float(pb[("level_pooled_pct", RK)][t]),
                          "n_prom": int(pb[("n_positions", PK)][t]),
                          "n_rest": int(pb[("n_positions", RK)][t])} for t in ("T1", "T2", "T3")}
    res["panel_c"] = {"n": int(len(L4)), "rho": float(g["r4"]), "p": float(g["p4"]),
                      "ymin": float(L4.level_pooled_pct.min()), "ymax": float(L4.level_pooled_pct.max()),
                      "cls_counts": {k: int(v) for k, v in L4.cls.value_counts().items()}}
    res["panel_d"] = [{"mark": r.mark_class, "tp": r.timepoint, "ratio": float(r.ratio)}
                      for r in dpan.itertuples()]

elif which == "s12":
    g = runpy.run_path(str(SCRIPTS / "recompute_261015/aagcccg261016/make_fig_261016.py"), run_name="__main__")
    LAD, PR, CTRL, C4 = g["LAD"], g["PR"], g["CTRL"], g["C4"]
    order = g["order"]
    res["ladder_n"] = {o: int(LAD[(LAD.definition == o) & (LAD.tp == "T1")].n_sites.iloc[0]) for o in order}
    res["ladder_occ"] = {o: {t: float(LAD[(LAD.definition == o) & (LAD.tp == t)].occupancy_pct.iloc[0])
                             for t in ("T1", "T2", "T3")} for o in order}
    res["c4_dashed"] = float(g["c4"])
    res["profile_minq"] = {t: float(PR[PR.tp == t].q.min()) for t in ("T1", "T2", "T3")}
    res["profile_q_at"] = {t: {int(r.offset): float(r.q) for r in PR[PR.tp == t].itertuples()}
                           for t in ("T1", "T2", "T3")}
    res["ctrl_T1"] = {ch: {b: float(CTRL[(CTRL.channel == ch) & (CTRL.tp == "T1") &
                                         (CTRL.base_at_minus3 == b)].occupancy_pct.iloc[0]) for b in "ACGT"}
                      for ch in g["chans"]}

elif which == "s13":
    g = runpy.run_path(str(SCRIPTS / "supp_261014/suppfigure14_determinants_261020.py"), run_name="__main__")
    lad2, ctl, d6H, d4A = g["lad2"], g["ctl"], g["d6H"], g["d4A"]
    res["ladder"] = [{"tp": r.timepoint, "tier": str(r.tier), "definition": r.definition,
                      "n": int(r.n), "called": int(r.called), "occ": float(r.occupancy_pct)}
                     for r in lad2.itertuples()]
    res["ctl"] = ctl.to_dict(orient="records")
    res["panel_c_means"] = {k: [float(ctl[k].mean())] for k in
                            ("gc8_mod", "gc8_unmod", "grun_mod", "grun_unmod", "fail_mod", "fail_unmod")}
    res["n_mod_unmod"] = [[int(r.n_mod), int(r.n_unmod)] for r in ctl.itertuples()]
    for nm, d in (("d6H", d6H), ("d4A", d4A)):
        t1 = d[d.timepoint == "T1"]
        best = t1.loc[t1.groupby("position").q_BH.idxmin()].sort_values("position")
        res[nm + "_best_T1"] = [{"pos": int(r.position), "base": r.base, "or": float(r.odds_ratio),
                                 "q": float(r.q_BH)} for r in best.itertuples()]
        res[nm + "_survives_pos"] = sorted(int(p) for p in set(d[d.survives_3tp].position))
        res[nm + "_n_tests"] = {t: int(d[d.timepoint == t].shape[0]) for t in ("T1", "T2", "T3")}

elif which == "s2":
    g = runpy.run_path(str(SCRIPTS / "supp_261014/suppfigure4_core_share.py"), run_name="__main__")
    CF, NULLT, POS = g["CF"], g["NULLT"], g["POS"]
    res["CF"] = CF.to_dict(orient="records")
    res["NULLT"] = {f"{k[0]}|{k[1]}": float(v) for k, v in NULLT.items()}
    res["NULL_allpos"] = {k: float(v) for k, v in g["NULL"].items()}
    res["n_pos"] = {f"{k[0]}|{k[1]}": len(v) for k, v in POS.items()}
    import numpy as np
    res["hist_max"] = {f"{k[0]}|{k[1]}": int(np.histogram(np.array(v), bins=np.arange(0, g["GLEN"] + 50000, 50000))[0].max())
                       for k, v in POS.items()}

elif which == "s4":
    import pandas as pd
    T = "/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/78_reviewer_robustness/tables"
    rates = pd.read_csv(f"{T}/P1_1_crosstalk_rates_C4.tsv", sep="\t")
    probs = pd.read_csv(f"{T}/P1_1_prob_distributions_C4.tsv", sep="\t")
    aag = rates[rates.group == "AAGCCCG_6mA"].iloc[0]
    non = rates[rates.group == "nonAAGCCCG_6mA_generic"].iloc[0]
    res["aag"] = {k: (float(v) if not isinstance(v, str) else v) for k, v in aag.items()}
    res["non"] = {k: (float(v) if not isinstance(v, str) else v) for k, v in non.items()}
    res["spec_fold"] = float(aag.near_rate / non.near_rate)
    res["prob_bins"] = [float(x) for x in probs.prob_bin_center]
    res["prob_A"] = [float(x) for x in (probs.AAGCCCG_A_6mA / probs.AAGCCCG_A_6mA.sum())]
    res["prob_C"] = [float(x) for x in (probs.AAGCCCG_C_4mC / probs.AAGCCCG_C_4mC.sum())]

json.dump(res, open(OUT / f"baseline_{which}.json", "w"), indent=1, default=str)
print("wrote", which, "keys:", list(res))
