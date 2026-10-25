"""正典部位集合の参照塩基監査: 6mA/4mC コールの参照塩基が A/T・C/G かを検定し、
現行値と補正値を並べる。規則は _gccggc と同一（depth=Nvalid+Nfail>=30, level>=0.5, 3 ライブラリ合算・和集合インデックス）。"""
import pandas as pd, numpy as np, re, sys
sys.path.insert(0,'supp_261014'); from _gccggc import genome
g=genome(); L=len(g)
PD="/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015"
COLS=["start","code","strand","Nvalid","Nmod","Nfail"]
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
aag_A=set(); tgg_C=set(); aag_C=set()
for m in re.finditer("AAGCCCG",g):
    i=m.start(); aag_A|={(i,"+"),(i+1,"+")}; aag_C.add((i+4,"+"))
for m in re.finditer("CGGGCTT",g):
    j=m.start(); aag_A|={(j+6,"-"),(j+5,"-")}; aag_C.add((j+2,"-"))
for m in re.finditer("TGGCCGGC",g): tgg_C.add((m.start()+4,"+"))
for m in re.finditer("GCCGGCCA",g): tgg_C.add((m.start()+3,"-"))
ORI=4270310; BIN=200000; nb=L//BIN+1
ob=ORI//BIN; dist=np.minimum(abs(np.arange(nb)-ob),nb-abs(np.arange(nb)-ob))
def oriter(pos):
    c=np.bincount(pos//BIN,minlength=nb)
    return round(c[dist<=2].mean()/c[dist>=18].mean(),2)
rows=[]
for tp,libs in TP.items():
    for code,mod,pb,mb,motif in (("a","6mA","A","T",aag_A),("21839","4mC","C","G",tgg_C|aag_C)):
        acc=[]
        for lib in libs:
            d=pd.read_csv(f"{PD}/{lib}_pileup.bed",sep=r"\s+",header=None,usecols=[1,3,5,9,11,15],names=COLS)
            acc.append(d[d.code==code].groupby(["start","strand"])[["Nvalid","Nmod","Nfail"]].sum())
        uni=acc[0].index
        for a in acc[1:]: uni=uni.union(a.index)
        S=sum(a.reindex(uni).fillna(0) for a in acc)
        S["lvl"]=S.Nmod/S.Nvalid.replace(0,np.nan); S=S.reset_index()
        cur=S[((S.Nvalid+S.Nfail)>=30)&(S.lvl>=0.5)].copy()
        cur["refok"]=np.where(cur.strand=="+",[g[i]==pb for i in cur.start],[g[i]==mb for i in cur.start])
        cor=cur[cur.refok]
        for lab,s in (("current (no ref-base filter)",cur),("corrected (ref base must match)",cor)):
            k=set(zip(s.start,s.strand)); asg=len(k&motif)
            rows.append(dict(timepoint=tp,mod=mod,rule=lab,n_sites=len(s),
                             assigned_to_recognition_seq=asg,
                             pct_unassigned=round(100*(1-asg/len(s)),2),
                             median_Nvalid=float(s.Nvalid.median()),
                             ori_ter_200kb=oriter(s.start.values)))
t=pd.DataFrame(rows)
t.to_csv("recompute_261015/refbase_audit_261015.tsv",sep="\t",index=False)
print(t.to_string(index=False))
