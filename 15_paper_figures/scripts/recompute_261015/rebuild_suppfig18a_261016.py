"""補足図 18a の領域検定表を正準規則で作り直す。

旧表 `supp_261014/data/neg261015/region_modified_vs_own_T1.csv` は HAAGCCCG 行が
19 インスタンス不足し（n_all 544、正準 561）、行ごとに母数の定義が混在していた（生成元は
リポジトリに無い）。本スクリプトは修飾位置を明示的に宣言したうえで全 7 行を作り直す。

修飾位置の宣言（推測していない根拠）:
  TGGCCGGC  4mC 位置 5        — 正典
  AAGCCCG   4mC 位置 5 (C4)   — 正典
  HAAGCCCG  6mA 位置 1 (A0)   — 正典。H = A/C/T
  GAAGCCCG  6mA 位置 1 (A0)   — 正典
  CGGCAACC  6mA 6 位          — 旧表と同じ（81/902 系列の定義）
  CCGKCA    6mA 6 位          — このモチーフのアデニンは 1 つだけなので一意。K = G/T
  GAACCGG   6mA 2 位と 3 位   — アデニンが 2 つあるので「アデニン位置あたり」の規約を使う
                               （AAGCCCG の A1/A2 と同じ扱い。行に明記する）
"""
import numpy as np, pandas as pd, re, sys, pathlib
from scipy import stats
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path("recompute_261015"); C=HERE/"pileup_cache_261015"
g=genome(); L=len(g); comp=str.maketrans("ACGT","TGCA")
rc=lambda s:s.translate(comp)[::-1]
IUPAC={"A":"A","C":"C","G":"G","T":"T","K":"[GT]","H":"[ACT]","N":"."}
def regex(m): return "".join(IUPAC[c] for c in m)
SPEC=[("CCGKCA","a",[5]),("GAACCGG","a",[1,2]),("TGGCCGGC","21839",[4]),
      ("AAGCCCG","21839",[4]),("CGGCAACC","a",[5]),("GAAGCCCG","a",[1]),("HAAGCCCG","a",[1])]
def positions(motif, offs):
    out=[]
    r=regex(motif)
    for m in re.finditer("(?="+r+")",g):
        for o in offs: out.append((m.start()+o,True))
    for m in re.finditer("(?="+r+")",rc(g)):
        for o in offs: out.append((L-1-(m.start()+o),False))
    return sorted(set(out))
TSS="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
t=pd.read_csv(TSS); sel=t[t.category.astype(str).str.upper().isin(["P","S"])&t.gene_id.notna()
      &(t.gene_id.astype(str).str.strip().str.len()>0)]
tss=np.sort(sel.tss_position.astype(int).values)
GFF="/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"
cds=np.zeros(L,bool)
for line in open(GFF):
    if line.startswith("#"): continue
    f=line.rstrip().split("\t")
    if len(f)<8 or f[2]!="CDS": continue
    cds[int(f[3])-1:int(f[4])]=True
def pooled(code):
    fr=[]
    for lib in ("1-1","1-2","1-3"):
        z=np.load(C/f"{lib}.npz"); m=(z["code"]==code)&z["refok"]
        fr.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],"Nvalid":z["Nvalid"][m],
            "Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m]}).set_index(["start","plus"]))
    uni=fr[0].index
    for f in fr[1:]: uni=uni.union(f.index)
    s=sum(f.reindex(uni).fillna(0) for f in fr); s["dep"]=s.Nvalid+s.Nfail
    s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan); return s
S={"a":pooled(0),"21839":pooled(1)}
rows=[]
for motif,code,offs in SPEC:
    pos=positions(motif,offs)
    d=pd.DataFrame(pos,columns=["start","plus"]).set_index(["start","plus"]).join(S[code],how="left")
    d=d[d.dep>=30].reset_index()
    if not len(d): continue
    d["called"]=d.lvl>=0.5
    p=d.start.to_numpy(); i=np.searchsorted(tss,p)
    lo=np.where(i>0,p-tss[np.clip(i-1,0,None)],10**9)
    hi=np.where(i<len(tss),tss[np.clip(i,None,len(tss)-1)]-p,10**9)
    # 宣言されている 4 クラス: promoter(-300..TSS) / 5'UTR(TSS..+50) / CDS / intergenic
    nearest=np.where(lo<=hi,-lo,hi)          # 負=TSS より上流
    d["prom"]=(nearest<0)&(np.abs(nearest)<=300)
    d["utr5"]=(nearest>=0)&(nearest<=50)
    d["cds"]=cds[p]&(~d.prom)&(~d.utr5)
    d["other"]=(~d.prom)&(~d.utr5)&(~d.cds)
    m=d[d.called]
    tab=np.array([[m.prom.sum(),m.utr5.sum(),m.cds.sum(),m.other.sum()],
                  [d.prom.sum(),d.utr5.sum(),d.cds.sum(),d.other.sum()]],float)
    tab[1]-=tab[0]
    keep=tab.sum(0)>0
    chi2,pv=(np.nan,np.nan)
    if keep.sum()>=2 and tab[:,keep].min()>=0 and m.shape[0]>=5:
        chi2,pv,_,_=stats.chi2_contingency(tab[:,keep])
    rows.append(dict(motif=motif,mark="6mA" if code=="a" else "4mC",
        modified_positions="/".join(str(o+1) for o in offs),
        unit="adenine positions" if len(offs)>1 else "strand-oriented positions",
        n_mod=int(m.shape[0]),n_all=int(d.shape[0]),
        prom_mod=round(100*m.prom.mean(),1),prom_all=round(100*d.prom.mean(),1),
        utr5_mod=round(100*m.utr5.mean(),1),utr5_all=round(100*d.utr5.mean(),1),
        cds_mod=round(100*m.cds.mean(),1),cds_all=round(100*d.cds.mean(),1),
        d_prom=round(100*(m.prom.mean()-d.prom.mean()),1),
        chi2=round(chi2,2) if chi2==chi2 else None,
        p=pv if pv==pv else None,
        passes_bonferroni=bool(pv<0.05/len(SPEC)) if pv==pv else False))
out=pd.DataFrame(rows)
out.to_csv("supp_261014/data/neg261015/region_modified_vs_own_T1.csv",index=False)
out.to_csv(HERE/"suppfig18a_region_rebuilt_261016.tsv",sep="\t",index=False)
print(out.to_string(index=False))
print("\nBonferroni alpha =", round(0.05/len(SPEC),5), " 通過:", int(out.passes_bonferroni.sum()))
