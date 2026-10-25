"""補正後の未帰属 6mA 部位と、弱く浮かんだ候補モチーフの特徴づけ。
位置（プロモーター／CDS／その他）、コア・アーム、BGC、既知モチーフとの距離、
局所 GC、ホモポリマー文脈、遺伝子鎖との一致、時点間再現性を、GC・鎖整合の
ランダムアデニン対照と比べる。"""
import numpy as np, pandas as pd, pathlib, re, sys
from scipy import stats
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path(__file__).resolve().parent
C=HERE/"pileup_cache_261015"
g=genome(); L=len(g); comp=str.maketrans("ACGT","TGCA"); rc=lambda s:s.translate(comp)[::-1]
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
CORE=(1_500_000,7_170_000)
# ---- 既知モチーフ
aag_A=set(); tgg=[]; aag=[]
for m in re.finditer("AAGCCCG",g): aag_A|={(m.start(),True),(m.start()+1,True)}; aag.append(m.start())
for m in re.finditer("CGGGCTT",g): aag_A|={(m.start()+6,False),(m.start()+5,False)}; aag.append(m.start())
for m in re.finditer("TGGCCGGC",g): tgg.append(m.start())
for m in re.finditer("GCCGGCCA",g): tgg.append(m.start())
tgg=np.sort(np.array(tgg)); aag=np.sort(np.array(aag))
# ---- 注釈
TSS_CSV=("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
         "18_tss_analyses/jeong2016_all_tss.csv")
t=pd.read_csv(TSS_CSV)
sel=t[t.category.astype(str).str.upper().isin(["P","S"])&t.gene_id.notna()
      &(t.gene_id.astype(str).str.strip().str.len()>0)]
tss=np.sort(sel.tss_position.astype(int).values)
GFF=("/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/"
     "39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff")
cds=np.zeros(L,np.int8)   # 0 none, 1 plus, 2 minus, 3 both
for line in open(GFF):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")
    if len(f)<8 or f[2]!="CDS": continue
    s,e,st=int(f[3])-1,int(f[4]),f[6]
    cds[s:e] |= (1 if st=="+" else 2)
bgc=pd.read_csv(HERE.parent/"supp_261014/data/bgc_regions_29.csv")
bcol=[c for c in bgc.columns if c.lower() in ("start","region_start","from")][0]
ecol=[c for c in bgc.columns if c.lower() in ("end","region_end","to")][0]
inbgc=np.zeros(L,bool)
for s,e in zip(bgc[bcol].astype(int),bgc[ecol].astype(int)): inbgc[s:e]=True
def mind(pos,arr):
    i=np.searchsorted(arr,pos)
    lo=np.where(i>0,pos-arr[np.clip(i-1,0,None)],10**9)
    hi=np.where(i<len(arr),arr[np.clip(i,None,len(arr)-1)]-pos,10**9)
    return np.minimum(lo,hi)
def homopoly(pos,k=5,win=5):
    out=[]
    for p in pos:
        s=g[max(0,p-win):p+win+1]
        out.append(any(b*k in s for b in "ACGT"))
    return np.array(out)
def annotate(pos,plus):
    pos=np.asarray(pos); plus=np.asarray(plus,bool)
    d=pd.DataFrame({"pos":pos,"plus":plus})
    d["d_tss"]=mind(pos,tss); d["promoter"]=d.d_tss<=300
    c=cds[pos]
    d["in_cds"]=c>0
    d["cds_sense"]=np.where(plus,(c&1)>0,(c&2)>0)
    d["intergenic"]=(c==0)&(~d.promoter)
    d["core"]=(pos>=CORE[0])&(pos<CORE[1])
    d["bgc"]=inbgc[pos]
    d["d_TGGCCGGC"]=mind(pos,tgg); d["near_TGGCCGGC_200"]=d.d_TGGCCGGC<=200
    d["d_AAGCCCG"]=mind(pos,aag)
    d["gc50"]=[ (g[max(0,p-25):p+26].count("G")+g[max(0,p-25):p+26].count("C"))/51 for p in pos]
    d["homopolymer"]=homopoly(pos)
    return d
# ---- 部位集合（補正規則）
def sites(tp,code=0):
    frames=[]
    for lib in TP[tp]:
        z=np.load(C/f"{lib}.npz"); m=(z["code"]==code)
        frames.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],
            "Nvalid":z["Nvalid"][m],"Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m],
            "refok":z["refok"][m]}).set_index(["start","plus"]))
    uni=frames[0].index
    for f in frames[1:]: uni=uni.union(f.index)
    s=sum(f[["Nvalid","Nmod","Nfail"]].reindex(uni).fillna(0) for f in frames)
    ok=frames[0]["refok"].reindex(uni)
    for f in frames[1:]: ok=ok.fillna(f["refok"].reindex(uni))
    s["refok"]=ok.fillna(False).astype(bool)
    s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan)
    s=s.reset_index()
    return s[((s.Nvalid+s.Nfail)>=30)&(s.lvl>=0.5)&s.refok]
U={}
for tp in TP:
    s=sites(tp)
    U[tp]=set(zip(s.start,s.plus))-aag_A
print({k:len(v) for k,v in U.items()},flush=True)
t1=sorted(U["T1"])
fg=annotate([p for p,_ in t1],[b for _,b in t1])
# GC・鎖整合背景
rng=np.random.default_rng(1)
arr=np.frombuffer(g.encode(),np.uint8)
pool={True:np.flatnonzero(arr==ord("A")),False:np.flatnonzero(arr==ord("T"))}
bp,bs=[],[]
for (p,pl),gc in zip(t1,fg.gc50):
    for _ in range(400):
        q=int(rng.choice(pool[bool(pl)]))
        if q-25<0 or q+26>L: continue
        w=g[q-25:q+26]
        if abs((w.count("G")+w.count("C"))/51-gc)<=0.02: bp.append(q); bs.append(pl); break
bg_=annotate(bp,bs)
rows=[]
for feat in ["promoter","in_cds","cds_sense","intergenic","core","bgc","near_TGGCCGGC_200","homopolymer"]:
    a=int(fg[feat].sum()); b=len(fg)-a; c_=int(bg_[feat].sum()); d_=len(bg_)-c_
    orr,p=stats.fisher_exact([[a,b],[c_,d_]])
    rows.append(dict(set="unassigned_6mA_T1",feature=feat,n=len(fg),fg_pct=round(100*a/len(fg),2),
                     bg_pct=round(100*c_/len(bg_),2),odds_ratio=round(orr,3),p=f"{p:.3g}"))
for col in ["d_tss","d_TGGCCGGC","d_AAGCCCG","gc50"]:
    u,p=stats.mannwhitneyu(fg[col],bg_[col])
    rows.append(dict(set="unassigned_6mA_T1",feature=col+" (median)",n=len(fg),
                     fg_pct=round(float(fg[col].median()),4),bg_pct=round(float(bg_[col].median()),4),
                     odds_ratio=np.nan,p=f"{p:.3g}"))
# 再現性
rows.append(dict(set="unassigned_6mA_T1",feature="also unassigned at T2",n=len(U["T1"]),
                 fg_pct=round(100*len(U["T1"]&U["T2"])/len(U["T1"]),2),bg_pct=np.nan,odds_ratio=np.nan,p=""))
rows.append(dict(set="unassigned_6mA_T1",feature="also unassigned at T3",n=len(U["T1"]),
                 fg_pct=round(100*len(U["T1"]&U["T3"])/len(U["T1"]),2),bg_pct=np.nan,odds_ratio=np.nan,p=""))
rows.append(dict(set="unassigned_6mA_T1",feature="unassigned at all three",n=len(U["T1"]),
                 fg_pct=round(100*len(U["T1"]&U["T2"]&U["T3"])/len(U["T1"]),2),bg_pct=np.nan,odds_ratio=np.nan,p=""))
# ---- 候補モチーフ別: ゲノム全出現を メチル化済 / 非メチル化 に割って同じ特徴を検定
def instances(pattern,aoff):
    """pattern は N を含んでよい。戻り値 (pos_of_A, plus)"""
    reg=pattern.replace("N",".")
    out=[]
    for m in re.finditer("(?="+reg+")",g): out.append((m.start()+aoff,True))
    for m in re.finditer("(?="+reg+")",rc(g)): out.append((L-1-(m.start()+aoff),False))
    return out
CAND=[("GGCAACCC",4),("CGGCAACC",5),("GGCAACC",4),("GGCGGNNNNGTCAC",12),("GCGGNNNNGTCAC",11)]
meth=U["T1"]
for pat,off in CAND:
    inst=instances(pat,off)
    if not inst: continue
    d=annotate([p for p,_ in inst],[b for _,b in inst])
    d["methylated"]=[(p,b) in meth for p,b in inst]
    nm=int(d.methylated.sum())
    if nm<8: continue
    for feat in ["promoter","in_cds","cds_sense","bgc","core","homopolymer","near_TGGCCGGC_200"]:
        a=int(d[d.methylated][feat].sum()); b=nm-a
        c_=int(d[~d.methylated][feat].sum()); d_=len(d)-nm-c_
        orr,p=stats.fisher_exact([[a,b],[c_,d_]])
        rows.append(dict(set=f"{pat} (meth {nm}/{len(d)})",feature=feat,n=len(d),
                         fg_pct=round(100*a/nm,2),bg_pct=round(100*c_/max(len(d)-nm,1),2),
                         odds_ratio=round(orr,3),p=f"{p:.3g}"))
    for col in ["gc50","d_TGGCCGGC"]:
        u,p=stats.mannwhitneyu(d[d.methylated][col],d[~d.methylated][col])
        rows.append(dict(set=f"{pat} (meth {nm}/{len(d)})",feature=col+" (median)",n=len(d),
                         fg_pct=round(float(d[d.methylated][col].median()),4),
                         bg_pct=round(float(d[~d.methylated][col].median()),4),
                         odds_ratio=np.nan,p=f"{p:.3g}"))
r=pd.DataFrame(rows)
r.to_csv(HERE/"unassigned_6mA_characterisation_261015.tsv",sep="\t",index=False)
fg.assign(timepoint="T1").to_csv(HERE/"unassigned_6mA_T1_annotated_261015.tsv",sep="\t",index=False)
print(r.to_string(index=False))
