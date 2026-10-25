#!/usr/bin/env python3
"""AAGCCCG 6mA の鎖方向性テスト用インスタンス表を組む（宿主遺伝子つき）。

出力: strand_instances_261021.tsv
  1 行 = 1 モチーフインスタンス × 1 時点
  列: timepoint, pos_6mA(0-based), motif_strand, gene, gene_strand, sense(1=コード鎖),
      Nvalid, Nmod, Nfail, measurable, modified, codon_phase, gc_<w>, L<k> (側方塩基)

占有率の定義は build_supptables_pileup.py と同一:
  measurable = (Nvalid+Nfail) >= 30 かつ Nvalid > 0
  modified   = measurable かつ Nmod/Nvalid >= 0.5
3 ライブラリ/時点をプールしてから判定する。

実行: cd scripts/recompute_261015 && python3 strand_cluster_261021_build.py
"""
import pathlib
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'supp_261014'))
from _gccggc import DEPTH_MIN, FRAC_MIN, LIBS, find_all, genome  # noqa: E402

GFF = ('/Users/okaban/bioinfo/rna-seq/11_epigenome_integration/analysis/'
       '39_GCCGGC_MTase_reverse_ID/data/GCF_000203835.1_ASM20383v1_genomic.gff')
CACHE = HERE / 'pileup_cache_261015'
MOTIF = 'AAGCCCG'
RC = 'CGGGCTT'          # revcomp(AAGCCCG)
MA_OFFSET = 0           # 6mA はモチーフ位置 1 (0-based 0)
LAT_UP, LAT_DOWN = 9, 9  # 感度解析用に広めに保存。主解析は上流 4 + 下流 5 を使う
GC_WINDOWS = (101, 201, 501, 1001)
COMP = str.maketrans('ACGT', 'TGCA')


def revcomp(s):
    return s.translate(COMP)[::-1]


# ---- モチーフインスタンス（両鎖） -------------------------------------------
g = genome()
inst = []
for i in find_all(g, MOTIF):
    inst.append(dict(motif_start=i, motif_strand='+', pos_6mA=i + MA_OFFSET))
for i in find_all(g, RC):
    # マイナス鎖上のモチーフ: ゲノム i..i+6 を覆い、モチーフ位置 1 は i+6
    inst.append(dict(motif_start=i, motif_strand='-', pos_6mA=i + 6 - MA_OFFSET))
inst = pd.DataFrame(inst).sort_values('pos_6mA').reset_index(drop=True)
print(f'AAGCCCG インスタンス 全ゲノム {len(inst):,} '
      f'(+{(inst.motif_strand == "+").sum():,} / -{(inst.motif_strand == "-").sum():,})')


# ---- 側方塩基・局所 GC（モチーフ方向で） -------------------------------------
def context(row):
    s, st = row.motif_start, row.motif_strand
    up = g[s - LAT_UP:s]
    down = g[s + 7:s + 7 + LAT_DOWN]
    seq = up + g[s:s + 7] + down
    if st == '-':
        seq = revcomp(seq)
    # モチーフ方向での並び: 上流 LAT_DOWN? -> revcomp 後は上流 = 元の下流
    # 正しくは: マイナス鎖モチーフの上流は g[s+7:] の revcomp 側
    return seq


lat = []
for row in inst.itertuples(index=False):
    s, st = row.motif_start, row.motif_strand
    if st == '+':
        up = g[s - LAT_UP:s]
        down = g[s + 7:s + 7 + LAT_DOWN]
    else:
        # モチーフ方向の上流 = ゲノム座標で下流側、revcomp する
        up = revcomp(g[s + 7:s + 7 + LAT_UP])
        down = revcomp(g[s - LAT_DOWN:s])
    lat.append(list(up) + list(down))
latcols = [f'L{-k}' for k in range(LAT_UP, 0, -1)] + [f'Lp{k}' for k in range(1, LAT_DOWN + 1)]
inst[latcols] = pd.DataFrame(lat, index=inst.index)

garr = np.frombuffer(g.encode(), dtype=np.uint8)
is_gc = ((garr == ord('G')) | (garr == ord('C'))).astype(np.int32)
cum = np.concatenate([[0], np.cumsum(is_gc)])
for w in GC_WINDOWS:
    half = w // 2
    lo = np.clip(inst.pos_6mA.to_numpy() - half, 0, len(g))
    hi = np.clip(inst.pos_6mA.to_numpy() + half + 1, 0, len(g))
    inst[f'gc{w}'] = (cum[hi] - cum[lo]) / (hi - lo)


# ---- 宿主遺伝子 -------------------------------------------------------------
# 宿主遺伝子の帰属は GFF の `gene` フィーチャで行う（既存の SuppTable27 と一致。
# TGGCCGGC で n_sense=1360/occ 90.0000%、n_anti=585/occ 91.4530% を再現して確認済み）。
# コドン位相だけは重なる CDS の読み枠から取る（RNA 遺伝子は phase=-1）。
def read_feats(kind):
    out = []
    with open(GFF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != 'NC_003888.3' or f[2] != kind:
                continue
            a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
            out.append((int(f[3]) - 1, int(f[4]) - 1, f[6],
                        a.get('locus_tag') or a.get('ID'),
                        a.get('gene_biotype', ''), f[7]))
    return pd.DataFrame(out, columns=['s', 'e', 'st', 'gene', 'biotype', 'phase'])


def assign(pos, tab):
    s, e = tab.s.to_numpy(), tab.e.to_numpy()
    o = np.argsort(s)
    s, e = s[o], e[o]
    t = tab.iloc[o].reset_index(drop=True)
    idx = []
    for p in pos:
        j = np.searchsorted(s, p, side='right') - 1
        f = [k for k in range(max(0, j - 8), min(len(s), j + 2)) if s[k] <= p <= e[k]]
        idx.append(f[0] if f else -1)
    return np.array(idx), t


GENE = read_feats('gene')
CDS = read_feats('CDS')
print(f'gene {len(GENE):,} / CDS {len(CDS):,}')

# 帰属アンカーはモチーフのゲノム左端 (motif_start)。この選択で既存 SuppTable27 の
# AAGCCCG 6mA を全桁再現する（T1 sense 247/708=34.887006%、anti 45/448=10.044643%）。
gi, gt = assign(inst.motif_start.to_numpy(), GENE)
inst['gene'] = np.where(gi >= 0, gt.gene.reindex(np.clip(gi, 0, None)).to_numpy(), None)
inst['gene_strand'] = np.where(gi >= 0, gt.st.reindex(np.clip(gi, 0, None)).to_numpy(), None)
inst['biotype'] = np.where(gi >= 0, gt.biotype.reindex(np.clip(gi, 0, None)).to_numpy(), None)
print(f'遺伝子内 {(gi >= 0).sum():,} / 遺伝子間 {(gi < 0).sum():,}')

ci, ct = assign(inst.motif_start.to_numpy(), CDS)
inst['cds_start'] = np.where(ci >= 0, ct.s.reindex(np.clip(ci, 0, None)).to_numpy(), -1)
inst['cds_end'] = np.where(ci >= 0, ct.e.reindex(np.clip(ci, 0, None)).to_numpy(), -1)
inst['cds_phase'] = np.where(ci >= 0,
                             pd.to_numeric(ct.phase.reindex(np.clip(ci, 0, None)),
                                           errors='coerce').fillna(0).to_numpy(), -1)

ing = inst[gi >= 0].copy()
ing['sense'] = (ing.motif_strand == ing.gene_strand).astype(int)

# コドン位相: 6mA 塩基が宿主 CDS の読み枠で何番目か (0/1/2)、CDS 外は -1
off_plus = ing.pos_6mA - ing.cds_start - ing.cds_phase
off_minus = ing.cds_end - ing.pos_6mA - ing.cds_phase
ph = np.where(ing.gene_strand == '+', off_plus, off_minus) % 3
ing['codon_phase'] = np.where(ing.cds_start >= 0, ph, -1).astype(int)


# ---- パイルアップ（時点プール） ----------------------------------------------
def pooled_6mA(tp):
    acc = {}
    for lib in LIBS[tp]:
        d = np.load(CACHE / f'{lib}.npz')
        m = (d['code'] == 0) & d['refok']
        st, pl, nv, nm, nf = (d['start'][m], d['plus'][m], d['Nvalid'][m],
                              d['Nmod'][m], d['Nfail'][m])
        for a, b, c, e, f in zip(st, pl, nv, nm, nf):
            k = (int(a), bool(b))
            v = acc.get(k)
            if v is None:
                acc[k] = [int(c), int(e), int(f)]
            else:
                v[0] += int(c); v[1] += int(e); v[2] += int(f)
    return acc


rows = []
for tp in ('T1', 'T2', 'T3'):
    acc = pooled_6mA(tp)
    sub = ing.copy()
    keys = list(zip(sub.pos_6mA.astype(int), sub.motif_strand == '+'))
    vals = [acc.get(k, [0, 0, 0]) for k in keys]
    sub['Nvalid'] = [v[0] for v in vals]
    sub['Nmod'] = [v[1] for v in vals]
    sub['Nfail'] = [v[2] for v in vals]
    sub['measurable'] = (((sub.Nvalid + sub.Nfail) >= DEPTH_MIN) & (sub.Nvalid > 0)).astype(int)
    frac = np.where(sub.Nvalid > 0, sub.Nmod / sub.Nvalid.replace(0, np.nan), np.nan)
    sub['frac'] = frac
    sub['modified'] = ((sub.measurable == 1) & (sub.frac >= FRAC_MIN)).astype(int)
    sub['timepoint'] = tp
    rows.append(sub)
    m = sub[sub.measurable == 1]
    for s, lab in ((1, 'sense'), (0, 'anti')):
        x = m[m.sense == s]
        print(f'{tp} {lab}: n={len(x)} modified={x.modified.sum()} '
              f'occ={100 * x.modified.mean():.6f}%')

out = pd.concat(rows, ignore_index=True)
out.to_csv(HERE / 'strand_instances_261021.tsv', sep='\t', index=False)
print('wrote strand_instances_261021.tsv', out.shape)
