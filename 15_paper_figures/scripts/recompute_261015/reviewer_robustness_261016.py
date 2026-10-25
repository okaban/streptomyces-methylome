"""査読対策の 3 本立て。正典規則（フィルタ済み site-level、真深度>=30、level>=0.5）で統一。
 (1) ラダー最良クラスのブートストラップ CI と、同一インスタンス上の 6mA 対 4mC の McNemar 検定
 (2) TGGCCGGC 4mC の両側決定基走査（±20 bp、事前登録規則と同形）
 (3) 6mA 決定基走査のフィルタ済み版（採用した -5C/-3G が規則を変えても残るか）
"""
import numpy as np, pandas as pd, re, sys, pathlib, json
from scipy import stats
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool
sys.path.insert(0,'supp_261014'); from _gccggc import genome
HERE=pathlib.Path(__file__).resolve().parent; C=HERE/"pileup_cache_261015"
g=genome(); L=len(g); comp=str.maketrans("ACGT","TGCA")
TP={"T1":["1-1","1-2","1-3"],"T2":["2-1","2-3","2-4"],"T3":["3-2","3-3","3-4"]}
NBOOT=10000; W=20

def pooled(tp, code):
    fr=[]
    for lib in TP[tp]:
        z=np.load(C/f"{lib}.npz"); m=(z["code"]==code)&z["refok"]
        fr.append(pd.DataFrame({"start":z["start"][m],"plus":z["plus"][m],"Nvalid":z["Nvalid"][m],
            "Nmod":z["Nmod"][m],"Nfail":z["Nfail"][m]}).set_index(["start","plus"]))
    uni=fr[0].index
    for f in fr[1:]: uni=uni.union(f.index)
    s=sum(f.reindex(uni).fillna(0) for f in fr)
    s["dep"]=s.Nvalid+s.Nfail; s["lvl"]=s.Nmod/s.Nvalid.replace(0,np.nan)
    return s

def base_at(pos, plus, off):
    """off<0 は 5′ 側、off>0 は 3′ 側（鎖方向）。off=0 は修飾塩基そのもの。"""
    p = pos+off if plus else pos-off
    if p<0 or p>=L: return "N"
    b=g[p]; return b if plus else b.translate(comp)

aag=[(m.start(),True) for m in re.finditer("(?=AAGCCCG)",g)]
rcg=g.translate(comp)[::-1]
aag+=[(L-1-m.start(),False) for m in re.finditer("(?=AAGCCCG)",rcg)]
tgg=[(m.start()+4,True) for m in re.finditer("(?=TGGCCGGC)",g)]
tgg+=[(m.start()+3,False) for m in re.finditer("(?=GCCGGCCA)",g)]

def scan(sites, s6, label, inner_lo, inner_hi):
    """修飾されたインスタンスと されなかったインスタンスで、モチーフ外の塩基組成を比べる。
    inner_lo..inner_hi はモチーフ本体の offset 範囲（検定から除く）。"""
    d=pd.DataFrame(sites,columns=["start","plus"]).set_index(["start","plus"]).join(s6,how="left")
    d=d[d.dep>=30].reset_index(); d["called"]=d.lvl>=0.5
    rows=[]
    for off in range(-W,W+1):
        if inner_lo<=off<=inner_hi: continue
        b=np.array([base_at(p,pl,off) for p,pl in zip(d.start,d.plus)])
        for nt in "ACGT":
            a=int(((b==nt)&d.called.values).sum()); bb=int(d.called.sum())-a
            c=int(((b==nt)&~d.called.values).sum()); dd=int((~d.called).sum())-c
            if a+c<10: continue
            orr,p=stats.fisher_exact([[a,bb],[c,dd]])
            rows.append(dict(set=label,offset=off,base=nt,n_called=a,n_notcalled=c,
                             odds_ratio=round(orr,3),p=p))
    t=pd.DataFrame(rows)
    t["q"]=multipletests(t.p,method="fdr_bh")[1]
    t["n_tests"]=len(t)
    return t

def boot_one(seed):
    rng=np.random.default_rng(seed)
    idx=rng.integers(0,len(_B6),len(_B6))
    return float(_B6[idx].mean()), float(_B4[idx].mean())

if __name__=="__main__":
    out={}
    # ---------- (1) ラダー最良クラス
    sel=[(p,b) for p,b in aag
         if base_at(p,b,-1)!="G" and base_at(p,b,-3)=="G" and base_at(p,b,-5)=="C"]
    boot_rows=[]
    for tp in TP:
        s6=pooled(tp,0); s4=pooled(tp,1)
        q6=s6.reindex(pd.MultiIndex.from_tuples(sel,names=["start","plus"]))
        c4=[(p+4,b) if b else (p-4,b) for p,b in sel]
        q4=s4.reindex(pd.MultiIndex.from_tuples(c4,names=["start","plus"]))
        ok=(q6.dep>=30).values & (q4.dep>=30).values
        _B6=(q6.lvl.values>=0.5)[ok]; _B4=(q4.lvl.values>=0.5)[ok]
        globals()["_B6"]=_B6; globals()["_B4"]=_B4
        with Pool(10) as p:
            bs=p.map(boot_one, range(NBOOT), chunksize=200)
        a6=np.array([x[0] for x in bs]); a4=np.array([x[1] for x in bs])
        dif=a6-a4
        b=int(((_B6)&(~_B4)).sum()); c=int(((~_B6)&(_B4)).sum())
        mc=stats.binomtest(b,b+c,0.5) if b+c>0 else None
        boot_rows.append(dict(timepoint=tp,n_paired=int(ok.sum()),
            occ_6mA=round(100*_B6.mean(),2), occ_4mC=round(100*_B4.mean(),2),
            ci6_lo=round(100*np.percentile(a6,2.5),2), ci6_hi=round(100*np.percentile(a6,97.5),2),
            ci4_lo=round(100*np.percentile(a4,2.5),2), ci4_hi=round(100*np.percentile(a4,97.5),2),
            diff_pp=round(100*(_B6.mean()-_B4.mean()),2),
            diff_ci_lo=round(100*np.percentile(dif,2.5),2),
            diff_ci_hi=round(100*np.percentile(dif,97.5),2),
            mcnemar_b=b, mcnemar_c=c,
            mcnemar_p=round(mc.pvalue,4) if mc else None))
        print("boot",tp,boot_rows[-1],flush=True)
    pd.DataFrame(boot_rows).to_csv(HERE/"ladder_bootstrap_261016.tsv",sep="\t",index=False)
    # ---------- (2)(3) 決定基走査
    scans=[]
    for tp in TP:
        s6=pooled(tp,0); s4=pooled(tp,1)
        H=[(p,b) for p,b in aag if base_at(p,b,-1)!="G"]
        t=scan(H,s6,f"6mA_HAAGCCCG_{tp}",0,6); scans.append(t)
        t2=scan(tgg,s4,f"4mC_TGGCCGGC_{tp}",-4,3); scans.append(t2)
        print("scan",tp,"6mA tests",len(t),"4mC tests",len(t2),flush=True)
    S=pd.concat(scans,ignore_index=True)
    S.to_csv(HERE/"determinant_scan_filtered_261016.tsv",sep="\t",index=False)
    sig=S[S.q<0.05].copy()
    sig["tp"]=sig.set.str.extract(r"_(T\d)$")[0]; sig["mark"]=sig.set.str.split("_").str[0]
    rep=(sig.groupby(["mark","offset","base"]).tp.nunique()
           .reset_index().rename(columns={"tp":"n_timepoints"}))
    rep=rep[rep.n_timepoints==3].sort_values(["mark","offset"])
    rep.to_csv(HERE/"determinant_scan_replicated_261016.tsv",sep="\t",index=False)
    print("\n3 時点すべてで q<0.05 の位置:"); print(rep.to_string(index=False))
