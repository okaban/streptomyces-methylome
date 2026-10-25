"""Nvalid+Nfail>=30 (現行 _gccggc.eligible) vs Nvalid>=30 の比較監査。"""
import pandas as pd, numpy as np, re, sys, pathlib
sys.path.insert(0,'supp_261014')
from _gccggc import genome
g=genome(); L=len(g)
PD="/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015"
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
COLS=["start","code","strand","Nvalid","Nmod","Nfail"]
aag_A=set(); tgg_C=set(); aag_C=set()
for m in re.finditer("AAGCCCG",g):
    i=m.start(); aag_A|={(i,"+"),(i+1,"+")}; aag_C.add((i+4,"+"))
for m in re.finditer("CGGGCTT",g):
    j=m.start(); aag_A|={(j+6,"-"),(j+5,"-")}; aag_C.add((j+2,"-"))
for m in re.finditer("TGGCCGGC",g): tgg_C.add((m.start()+4,"+"))
for m in re.finditer("GCCGGCCA",g): tgg_C.add((m.start()+3,"-"))
ORI=4270310; BIN=200000
rows=[]; occ=[]
for tp,libs in TP.items():
    acc=[]
    for lib in libs:
        d=pd.read_csv(f"{PD}/{lib}_pileup.bed",sep=r"\s+",header=None,usecols=[1,3,5,9,11,15],names=COLS)
        acc.append(d.groupby(["start","strand","code"])[["Nvalid","Nmod","Nfail"]].sum())
    idx=acc[0].index
    S=sum(a.reindex(idx).fillna(0) for a in acc).reset_index()
    S["lvl"]=S.Nmod/S.Nvalid.replace(0,np.nan)
    S["refbase"]=[g[i] for i in S.start]
    for code,mod,refp,refm in (("a","6mA","A","T"),("21839","4mC","C","G")):
        s=S[S.code==code].copy()
        s["refok"]=np.where(s.strand=="+",s.refbase==refp,s.refbase==refm)
        for rule,mask in (("Nvalid+Nfail>=30",(s.Nvalid+s.Nfail)>=30),("Nvalid>=30",s.Nvalid>=30)):
            c=s[mask&(s.lvl>=0.5)]
            key=set(zip(c.start,c.strand))
            assigned=len(key&(aag_A if mod=="6mA" else (tgg_C|aag_C)))
            if len(c):
                b=((c.start.values-0)//BIN)
                ori_bin=ORI//BIN
                dist=np.minimum(abs(b-ori_bin),(L//BIN)-abs(b-ori_bin))
                near=(dist<=2).sum(); far=(dist>=18).sum()
            rows.append(dict(timepoint=tp,mod=mod,rule=rule,n_sites=len(c),
                             ref_base_correct=int(c.refok.sum()),
                             pct_ref_correct=round(100*c.refok.mean(),1) if len(c) else np.nan,
                             assigned_to_motif=assigned,
                             pct_unassigned=round(100*(1-assigned/len(c)),2) if len(c) else np.nan,
                             median_Nvalid=float(c.Nvalid.median()) if len(c) else np.nan))
    # motif occupancy under both rules
    for name,keys,code in (("TGGCCGGC_4mC",tgg_C,"21839"),("AAGCCCG_A1A2_6mA",aag_A,"a"),("AAGCCCG_C5_4mC",aag_C,"21839")):
        s=S[S.code==code].set_index(["start","strand"])
        sub=s.reindex(pd.MultiIndex.from_tuples(sorted(keys),names=["start","strand"]))
        for rule,mask in (("Nvalid+Nfail>=30",(sub.Nvalid+sub.Nfail)>=30),("Nvalid>=30",sub.Nvalid>=30)):
            e=sub[mask.fillna(False)]
            lv=e.Nmod/e.Nvalid.replace(0,np.nan)
            occ.append(dict(timepoint=tp,motif=name,rule=rule,n_total=len(keys),n_eligible=len(e),
                            pct_eligible=round(100*len(e)/len(keys),1),
                            occupancy_pct=round(100*(lv>=0.5).mean(),2) if len(e) else np.nan,
                            median_Nvalid=float(e.Nvalid.median()) if len(e) else np.nan))
t=pd.DataFrame(rows); o=pd.DataFrame(occ)
t.to_csv("recompute_261015/depth_rule_audit_sites_261015.tsv",sep="\t",index=False)
o.to_csv("recompute_261015/depth_rule_audit_motifs_261015.tsv",sep="\t",index=False)
print(t.to_string(index=False)); print(); print(o.to_string(index=False))
