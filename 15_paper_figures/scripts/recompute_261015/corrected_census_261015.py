"""キャッシュから正典部位集合を現行規則・補正規則の両方で作り、補足表 16 を作り直す。
補正規則 = 現行規則（3 ライブラリ合算・和集合インデックス・真深度 Nvalid+Nfail>=30・level>=0.5）
          ＋ 参照塩基が修飾種と一致すること。"""
import numpy as np, pandas as pd, pathlib, re, sys
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path(__file__).resolve().parent
C=HERE/"pileup_cache_261015"
g=genome(); L=len(g)
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
TSS_CSV=("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
         "18_tss_analyses/jeong2016_all_tss.csv")
t=pd.read_csv(TSS_CSV)
sel=t[t.category.astype(str).str.upper().isin(["P","S"])&t.gene_id.notna()
      &(t.gene_id.astype(str).str.strip().str.len()>0)]
tss=np.sort(sel.tss_position.astype(int).values)
def near(pos,half=300):
    if len(pos)==0: return 0
    p=np.sort(pos); i=np.searchsorted(tss,p)
    lo=np.where(i>0,p-tss[np.clip(i-1,0,None)],10**9)
    hi=np.where(i<len(tss),tss[np.clip(i,None,len(tss)-1)]-p,10**9)
    return int((np.minimum(lo,hi)<=half).sum())
aag_A=set(); tgg_C=set(); aag_C=set()
for m in re.finditer("AAGCCCG",g):
    i=m.start(); aag_A|={(i,1),(i+1,1)}; aag_C.add((i+4,1))
for m in re.finditer("CGGGCTT",g):
    j=m.start(); aag_A|={(j+6,0),(j+5,0)}; aag_C.add((j+2,0))
for m in re.finditer("TGGCCGGC",g): tgg_C.add((m.start()+4,1))
for m in re.finditer("GCCGGCCA",g): tgg_C.add((m.start()+3,0))
def pooled(libs,code):
    """和集合インデックスで 3 ライブラリを合算。戻り値 DataFrame(start,plus,Nvalid,Nmod,Nfail,refok)"""
    frames=[]
    for lib in libs:
        z=np.load(C/f"{lib}.npz")
        m=(z["code"]==code)
        frames.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],
            "Nvalid":z["Nvalid"][m],"Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m],
            "refok":z["refok"][m]}).set_index(["start","plus"]))
    uni=frames[0].index
    for f in frames[1:]: uni=uni.union(f.index)
    s=sum(f[["Nvalid","Nmod","Nfail"]].reindex(uni).fillna(0) for f in frames)
    s["refok"]=frames[0]["refok"].reindex(uni).groupby(level=[0,1]).first()
    bad=s.refok.isna()
    if bad.any():
        for f in frames[1:]:
            fill=f["refok"].reindex(uni[bad.values])
            s.loc[bad.values,"refok"]=s.loc[bad.values,"refok"].fillna(fill)
    return s.reset_index()
rows=[]; sitesets={}
for tp,libs in TP.items():
    for code,mod,motif in ((0,"6mA",aag_A),(1,"4mC",tgg_C|aag_C)):
        s=pooled(libs,code)
        s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan)
        cur=s[((s.Nvalid+s.Nfail)>=30)&(s.lvl>=0.5)]
        cor=cur[cur.refok.astype(bool)]
        for lab,d in (("current",cur),("corrected",cor)):
            k=set(zip(d.start,d.plus)); asg=k&motif; un=k-motif
            rows.append(dict(timepoint=tp,mod=mod,rule=lab,n_sites=len(k),
                in_recognition_sequence=len(asg),unassigned=len(un),
                unassigned_pct=round(100*len(un)/len(k),2),
                unassigned_near_TSS=near(np.array(sorted(p for p,_ in un))),
                median_Nvalid=float(d.Nvalid.median())))
            if lab=="corrected": sitesets[(tp,mod)]=d
np.save(HERE/"corrected_unassigned_T1_6mA.npy",
        np.array(sorted(set(zip(sitesets[("T1","6mA")].start,sitesets[("T1","6mA")].plus))-aag_A),
                 dtype=[("start","i4"),("plus","?")]))
r=pd.DataFrame(rows)
r["unassigned_near_TSS_pct"]=(100*r.unassigned_near_TSS/r.unassigned).round(2)
r.to_csv(HERE/"corrected_census_261015.tsv",sep="\t",index=False)
print(r.to_string(index=False))
# 補足表 16（補正版）
w=r[r.rule=="corrected"].pivot(index="timepoint",columns="mod")
st=pd.DataFrame({
 "sites_4mC":w[("n_sites","4mC")],"sites_6mA":w[("n_sites","6mA")],
 "mC_in_recognition_sequence":w[("in_recognition_sequence","4mC")],
 "mC_unassigned":w[("unassigned","4mC")],"mC_unassigned_pct":w[("unassigned_pct","4mC")],
 "mA_in_AAGCCCG":w[("in_recognition_sequence","6mA")],
 "mA_unassigned":w[("unassigned","6mA")],"mA_unassigned_pct":w[("unassigned_pct","6mA")],
 "mA_unassigned_near_TSS":w[("unassigned_near_TSS","6mA")],
 "mA_unassigned_near_TSS_pct":w[("unassigned_near_TSS_pct","6mA")]}).reset_index()
st.to_csv(HERE/"SuppTable16_unassigned_6mA_corrected_261015.tsv",sep="\t",index=False)
print(); print(st.to_string(index=False))
