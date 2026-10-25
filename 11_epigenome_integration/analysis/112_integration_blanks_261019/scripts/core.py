"""112_integration_blanks_261019 — 共有コア。
正典規則のみを実装する（corrected_census_261015.py と同一の合算・フィルタ）。
- 3 ライブラリ和集合合算 / 真深度 (Nvalid+Nfail) >= 30 / pooled Nmod/Nvalid >= 0.5 / refok==1
出力はすべて 0-based start, plus(bool)。
"""
import numpy as np, pandas as pd, pathlib, re, functools

ROOT = pathlib.Path("/Users/okaban/bioinfo/rna-seq")
CACHE = ROOT / "15_paper_figures/scripts/recompute_261015/pileup_cache_261015"
FNA = ROOT / "11_epigenome_integration/analysis/112_integration_blanks_261019/../../data/NC_003888.3.fna"
GFF = ROOT / "11_epigenome_integration/analysis/39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff"
FC = ROOT / "03_quantification/data/featureCounts/featureCounts_M145_s2.txt"
TSS_CSV = ROOT / "11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv"
DESEQ = ROOT / "04_deseq2/results/results"

LIBS = ["1-1", "1-2", "1-3", "2-1", "2-3", "2-4", "3-2", "3-3", "3-4"]
TP = {"T1": ["1-1", "1-2", "1-3"], "T2": ["2-1", "2-3", "2-4"], "T3": ["3-2", "3-3", "3-4"]}
ORIC = 4_270_310
CORE = (1_500_000, 7_170_000)
CODE = {"6mA": 0, "4mC": 1, "5mC": 2}
COMP = str.maketrans("ACGTN", "TGCAN")


@functools.lru_cache(1)
def genome() -> str:
    s = []
    with open(FNA) as fh:
        for line in fh:
            if not line.startswith(">"):
                s.append(line.strip().upper())
    return "".join(s)


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# ---------------------------------------------------------------- motifs
@functools.lru_cache(1)
def aagcccg_instances() -> pd.DataFrame:
    """AAGCCCG 両鎖インスタンス。1 行 = 1 インスタンス。
    pos_A1/pos_A2 = 第1/第2アデニンの 0-based 座標, pos_C4 = 5 番目の C（4mC 部位）。
    m5/m3 = 修飾鎖上での -5 / -3 位の塩基。"""
    g = genome()
    L = len(g)
    rows = []
    for m in re.finditer("(?=AAGCCCG)", g):
        i = m.start()
        up = g[max(0, i - 6):i]
        rows.append(dict(inst_id=f"P{i}", start=i, plus=True, pos_A1=i, pos_A2=i + 1,
                         pos_C4=i + 4,
                         up6=up, m5=up[-5] if len(up) >= 5 else "N",
                         m3=up[-3] if len(up) >= 3 else "N"))
    for m in re.finditer("(?=CGGGCTT)", g):
        j = m.start()          # minus-strand motif occupies j..j+6 on the plus strand
        i = j + 6              # plus-strand coordinate of the motif's first A (minus strand)
        up = rc(g[j + 7:j + 13])
        rows.append(dict(inst_id=f"M{j}", start=j, plus=False, pos_A1=j + 6, pos_A2=j + 5,
                         pos_C4=j + 2,
                         up6=up, m5=up[-5] if len(up) >= 5 else "N",
                         m3=up[-3] if len(up) >= 3 else "N"))
    d = pd.DataFrame(rows)
    d["mid"] = d.pos_A1.where(d.plus, d.pos_C4) + 0  # approximate anchor
    d["anchor"] = ((d.start + d.start + 6) // 2)
    return d.sort_values("anchor").reset_index(drop=True)


@functools.lru_cache(1)
def tggccggc_positions() -> pd.DataFrame:
    """TGGCCGGC の 4mC 部位（鎖指定）。"""
    g = genome()
    rows = []
    for m in re.finditer("(?=TGGCCGGC)", g):
        rows.append(dict(start=m.start() + 4, plus=True, anchor=m.start()))
    for m in re.finditer("(?=GCCGGCCA)", g):
        rows.append(dict(start=m.start() + 3, plus=False, anchor=m.start()))
    return pd.DataFrame(rows)


# ------------------------------------------------------- canonical sites
def pooled(libs, code) -> pd.DataFrame:
    frames = []
    for lib in libs:
        z = np.load(CACHE / f"{lib}.npz")
        m = z["code"] == code
        frames.append(pd.DataFrame({"start": z["start"][m], "plus": z["plus"][m],
                                    "Nvalid": z["Nvalid"][m], "Nmod": z["Nmod"][m],
                                    "Nfail": z["Nfail"][m], "refok": z["refok"][m]}
                                   ).set_index(["start", "plus"]))
        del z
    uni = frames[0].index
    for f in frames[1:]:
        uni = uni.union(f.index)
    s = sum(f[["Nvalid", "Nmod", "Nfail"]].reindex(uni).fillna(0) for f in frames)
    s["refok"] = frames[0]["refok"].reindex(uni).groupby(level=[0, 1]).first()
    bad = s.refok.isna()
    if bad.any():
        for f in frames[1:]:
            fill = f["refok"].reindex(uni[bad.values])
            s.loc[bad.values, "refok"] = s.loc[bad.values, "refok"].fillna(fill)
    s = s.reset_index()
    s["lvl"] = s.Nmod / s.Nvalid.replace(0, np.nan)
    s["depth"] = s.Nvalid + s.Nfail
    return s


@functools.lru_cache(16)
def canonical(tp: str, mod: str) -> pd.DataFrame:
    """正典フィルタ通過部位（高信頼部位）。"""
    s = pooled(TP[tp], CODE[mod])
    d = s[(s.depth >= 30) & (s.lvl >= 0.5) & (s.refok.astype(bool))]
    return d.reset_index(drop=True)


@functools.lru_cache(16)
def measurable(tp: str, mod: str) -> pd.DataFrame:
    """深度条件のみ通過（占有率の分母に使える可測位置）。"""
    s = pooled(TP[tp], CODE[mod])
    d = s[(s.depth >= 30) & (s.refok.astype(bool))]
    return d.reset_index(drop=True)


# ------------------------------------------------------------ expression
@functools.lru_cache(1)
def genes() -> pd.DataFrame:
    fc = pd.read_csv(FC, sep="\t", comment="#")
    fc.columns = list(fc.columns[:6]) + LIBS
    fc = fc.rename(columns={"Geneid": "gene_id", "Start": "gstart", "End": "gend",
                            "Strand": "gstrand", "Length": "length"})
    fc["gstart"] = fc.gstart.astype(str).str.split(";").str[0].astype(int)
    fc["gend"] = fc.gend.astype(str).str.split(";").str[-1].astype(int)
    fc["gstrand"] = fc.gstrand.astype(str).str.split(";").str[0]
    g = genome()
    fc["gc"] = [ (g[a-1:b].count("G") + g[a-1:b].count("C")) / max(1, b-a+1)
                 for a, b in zip(fc.gstart, fc.gend) ]
    # CPM / TPM
    cnt = fc[LIBS].astype(float)
    cpm = cnt / cnt.sum() * 1e6
    rpk = cnt.div(fc.length / 1000.0, axis=0)
    tpm = rpk / rpk.sum() * 1e6
    for lib in LIBS:
        fc[f"cpm_{lib}"] = cpm[lib]
        fc[f"tpm_{lib}"] = tpm[lib]
    for tp, libs in TP.items():
        fc[f"tpm_{tp}"] = tpm[libs].mean(axis=1)
        fc[f"cv_{tp}"] = tpm[libs].std(axis=1, ddof=1) / tpm[libs].mean(axis=1).replace(0, np.nan)
    fc["tpm_mean"] = tpm.mean(axis=1)
    fc["log_tpm"] = np.log2(fc.tpm_mean + 1)
    # replication strand: right replichore (oriC -> right end) leading = plus strand
    fc["gmid"] = (fc.gstart + fc.gend) / 2
    fc["replichore"] = np.where(fc.gmid >= ORIC, "right", "left")
    lead_plus = fc.replichore.eq("right")
    fc["leading"] = np.where(fc.gstrand.eq("+"), lead_plus, ~lead_plus)
    fc["in_core"] = (fc.gmid >= CORE[0]) & (fc.gmid <= CORE[1])
    return fc.reset_index(drop=True)


@functools.lru_cache(1)
def deseq() -> pd.DataFrame:
    out = {}
    for a, b in (("2", "1"), ("3", "1"), ("3", "2")):
        d = pd.read_csv(DESEQ / f"DESeq2_M145_{a}_vs_{b}.tsv", sep="\t")
        idc = d.columns[0]
        d = d.rename(columns={idc: "gene_id"})
        out[f"lfc_{a}v{b}"] = d.set_index("gene_id").log2FoldChange
        out[f"padj_{a}v{b}"] = d.set_index("gene_id").padj
    return pd.DataFrame(out).reset_index().rename(columns={"index": "gene_id"})


@functools.lru_cache(1)
def gff_products() -> pd.DataFrame:
    rows = []
    with open(GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("CDS", "gene"):
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            lt = at.get("locus_tag")
            if not lt:
                continue
            rows.append(dict(gene_id=lt, feat=f[2], product=at.get("product", ""),
                             gene_name=at.get("gene", "")))
    d = pd.DataFrame(rows)
    cds = d[d.feat == "CDS"].drop_duplicates("gene_id")
    gn = d[d.feat == "gene"].drop_duplicates("gene_id")[["gene_id", "gene_name"]]
    return cds[["gene_id", "product"]].merge(gn, on="gene_id", how="outer").fillna("")


# ------------------------------------------------- transcription units (TU)
def predict_tus(max_gap=50, min_r=None) -> pd.DataFrame:
    """同方向・遺伝子間距離 <= max_gap bp を連結した転写単位の予測（標準的ヒューリスティック）。
    戻り値は genes() に tu_id / tu_size / tu_pos(1-based) / tu_first / tu_last を足した表。"""
    g = genes().sort_values("gstart").reset_index(drop=True)
    tu, cur = [], 0
    prev_end, prev_strand = None, None
    for _, r in g.iterrows():
        if prev_end is None:
            cur = 0
        else:
            gap = r.gstart - prev_end - 1
            if not (r.gstrand == prev_strand and gap <= max_gap):
                cur += 1
        tu.append(cur)
        prev_end, prev_strand = max(prev_end or 0, r.gend), r.gstrand
    g["tu_id"] = tu
    sz = g.groupby("tu_id").gene_id.transform("size")
    g["tu_size"] = sz
    # position within TU, 5'->3'
    g["tu_rank_plus"] = g.groupby("tu_id").cumcount() + 1
    g["tu_pos"] = np.where(g.gstrand == "+", g.tu_rank_plus, g.tu_size - g.tu_rank_plus + 1)
    g["tu_first"] = g.tu_pos == 1
    g["tu_last"] = g.tu_pos == g.tu_size
    return g


# -------------------------------------------------------------- utilities
def assign_to_gene(pos: np.ndarray, gdf: pd.DataFrame) -> np.ndarray:
    """0-based 座標を遺伝子 index に割り当てる（重なりは先頭優先）。-1 = 遺伝子間。"""
    gs = gdf.gstart.values - 1
    ge = gdf.gend.values
    order = np.argsort(gs)
    gs_s, ge_s, idx_s = gs[order], ge[order], gdf.index.values[order]
    out = np.full(len(pos), -1, dtype=np.int64)
    j = np.searchsorted(gs_s, pos, side="right") - 1
    for k in range(len(pos)):
        jj = j[k]
        while jj >= 0 and pos[k] - gs_s[jj] < 300_000:
            if gs_s[jj] <= pos[k] < ge_s[jj]:
                out[k] = idx_s[jj]
                break
            jj -= 1
    return out


def local_gc(pos, half=250):
    g = genome()
    L = len(g)
    out = np.empty(len(pos))
    for i, p in enumerate(pos):
        a, b = max(0, p - half), min(L, p + half + 1)
        w = g[a:b]
        out[i] = (w.count("G") + w.count("C")) / len(w)
    return out


def bh(p):
    p = np.asarray(p, float)
    ok = ~np.isnan(p)
    q = np.full(p.shape, np.nan)
    pv = p[ok]
    n = len(pv)
    o = np.argsort(pv)
    r = np.empty(n)
    r[o] = np.arange(1, n + 1)
    qv = np.minimum.accumulate((pv[o] * n / r[o])[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.clip(qv, 0, 1)
    q[ok] = qq
    return q
