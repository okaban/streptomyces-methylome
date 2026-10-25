#!/usr/bin/env python3
"""生成元スクリプトが無い補足図 1・9・10・11・12 について、凡例に書かれた数値を
出所まで辿って照合する。図は一切描き直さない（PNG に触れない）。

経緯: この 5 枚は出荷 PNG のみで、ツリー内に描画スクリプトが存在しない
（deposition/figure_producer_map.tsv にも 4 件が no match として記録済み）。
数値が動いたとき図は直せないが、凡例のどの数値も追跡可能にしておくための措置。
著者指示 2026-10-15。

出力: recompute_261015/producerless_verification.tsv
  figure / quantity / claimed / recomputed / status / source

status の意味
  OK        … 凡例の値が出所から再現した（表示桁で一致）
  METHOD    … 点推定は一致するが区間などの算出法が凡例から復元できない
  STALE     … 出所が撤回済み・旧世代の集合を使っている（要対応）

実行: cd scripts && python3 recompute_261015/recompute_producerless.py
"""
import pathlib
import re

import numpy as np
import pandas as pd

RNA = pathlib.Path('/Users/okaban/bioinfo/rna-seq')
EPI = RNA / '11_epigenome_integration/analysis'
GENOME = RNA / '11_epigenome_integration/data/NC_003888.3.fna'
TSSCSV = EPI / '18_tss_analyses/jeong2016_all_tss.csv'
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / 'producerless_verification.tsv'

rows = []


def add(fig, qty, claimed, got, status, source):
    rows.append(dict(figure=fig, quantity=qty, claimed=claimed,
                     recomputed=got, status=status, source=source))


def close(a, b, tol):
    try:
        return abs(float(a) - float(b)) <= tol
    except (TypeError, ValueError):
        return False


# ---- 補足図 1: AAGCCCG 内の共修飾ペアの置換帰無 ----------------------------
def suppfig1():
    src = EPI / '87_site_coloc_stats_C4/tables/site_coloc_2x2_C4.tsv'
    d = pd.read_csv(src, sep='\t').set_index('scope')
    rel = str(src.relative_to(RNA))
    for scope, claim_obs, claim_exp, claim_max in [('pooled', 406, 0.64, 5),
                                                  ('T1', 214, 0.34, 4)]:
        r = d.loc[scope]
        lab = 'all timepoints pooled' if scope == 'pooled' else 'T1 only'
        add('SuppFig 1', f'observed co-modified pairs, {lab}', claim_obs,
            int(r.in_comod), 'OK' if int(r.in_comod) == claim_obs else 'MISMATCH', rel)
        add('SuppFig 1', f'null expectation, {lab}', claim_exp,
            round(float(r.expected_in_comod_H0), 3),
            'OK' if close(r.expected_in_comod_H0, claim_exp, 0.005) else 'MISMATCH', rel)
        add('SuppFig 1', f'null maximum, {lab}', claim_max,
            int(r.perm_null_max_in_comod),
            'OK' if int(r.perm_null_max_in_comod) == claim_max else 'MISMATCH', rel)
    r = d.loc['pooled']
    add('SuppFig 1', 'candidate same-strand A/C pairs at 3-4 bp', 1696280,
        int(r.n_candidate_pairs),
        'OK' if int(r.n_candidate_pairs) == 1696280 else 'MISMATCH', rel)
    add('SuppFig 1', 'Haldane-corrected odds ratio', '6.1e5',
        f'{float(r.OR):.3g}', 'OK' if close(r.OR, 6.1e5, 1.5e4) else 'MISMATCH', rel)


# ---- 補足図 9: near/far 共呼び出し ------------------------------------------
def suppfig9():
    src = EPI / '78_reviewer_robustness/tables/P1_1_crosstalk_rates.tsv'
    d = pd.read_csv(src, sep='\t').set_index('group')
    rel = str(src.relative_to(RNA))
    m, n = d.loc['AAGCCCG_6mA'], d.loc['nonAAGCCCG_6mA_generic']
    checks = [
        ('non-motif 6mA calls scanned', '2.1e6', f'{int(n.n):,}', close(n.n, 2.1e6, 6e4)),
        ('near/far at non-motif 6mA', 1.09, round(float(n.near_over_far), 2),
         close(n.near_over_far, 1.09, 0.006)),
        ('near-offset co-call rate at AAGCCCG', 0.0474, round(float(m.near_rate), 4),
         close(m.near_rate, 0.0474, 5e-5)),
        ('near-offset co-call rate at non-motif 6mA', 0.000110,
         round(float(n.near_rate), 6), close(n.near_rate, 0.000110, 1e-6)),
        ('within-motif near/far ratio', 365, round(float(m.near_over_far)),
         close(m.near_over_far, 365, 0.5)),
        ('motif / non-motif near-rate ratio', 430,
         round(float(m.near_rate / n.near_rate)),
         close(m.near_rate / n.near_rate, 430, 1.0)),
    ]
    for q, c, g, ok in checks:
        add('SuppFig 9', q, c, g, 'OK' if ok else 'MISMATCH', rel)


# ---- 補足図 11: TFBS 占有と発現五分位 ---------------------------------------
def suppfig11():
    s1 = EPI / '49_TF_BS_methylation_protection/tables/statistical_tests.tsv'
    d = pd.read_csv(s1, sep='\t').set_index('test_id')
    rel1 = str(s1.relative_to(RNA))
    r = d.loc['T1_overlap_GCCGGC_T1']
    add('SuppFig 11', 'methylation enrichment within TFBS (fold)', 1.157,
        round(float(r.fold), 3), 'OK' if close(r.fold, 1.157, 5e-4) else 'MISMATCH', rel1)
    add('SuppFig 11', 'methylation enrichment within TFBS (p)', 0.006,
        round(float(r.pvalue), 4), 'OK' if close(r.pvalue, 0.006, 5e-4) else 'MISMATCH', rel1)

    s2 = EPI / '52_shielded_exposed_boundary/tables/expression_quintile_JT_n62.tsv'
    j = pd.read_csv(s2, sep='\t').iloc[0]
    rel2 = str(s2.relative_to(RNA))
    add('SuppFig 11', 'Jonckheere-Terpstra p across expression quintiles', 0.818,
        round(float(j.p), 4), 'OK' if close(j.p, 0.818, 5e-4) else 'MISMATCH', rel2)
    # 凡例は「認識配列から 293 bp 以内かどうかによる層別」と書いているが、
    # p = 0.818 を出している検定の outcome 列は撤回済みの Exposed 分類である。
    add('SuppFig 11', 'outcome variable behind the quintile test',
        'proximity to a TGGCCGGC site within 293 bp', str(j.outcome), 'STALE', rel2)
    add('SuppFig 11', 'n behind the quintile test', '1,051 (canonical with TSS)',
        f'{int(j.n):,}', 'STALE', rel2)


# ---- 補足図 12: GC 整合後のプロモーター近傍希少性 ---------------------------
def suppfig12(nboot=2000, seed=42):
    HALF, NQ = 300, 5
    WIN = 2 * HALF + 1
    g = ''.join(l.strip() for l in open(GENOME) if not l.startswith('>')).upper()
    glen = len(g)
    b = np.frombuffer(g.encode(), np.uint8)
    cum = np.concatenate([[0], np.cumsum(((b == ord('G')) | (b == ord('C'))).astype(np.int64))])
    t = pd.read_csv(TSSCSV)
    cat = t['category'].astype(str).str.strip().str.upper()
    gid = t['gene_id'].astype(str).str.strip()
    tss = t[cat.isin(['P', 'S']) & gid.ne('') & gid.ne('nan')]['tss_position'].astype(int).to_numpy()
    tss = tss[(tss - HALF >= 0) & (tss + HALF + 1 <= glen)]
    st = tss - HALF
    bg = np.arange(0, glen - WIN, WIN)

    def gcf(s):
        return (cum[s + WIN] - cum[s]) / WIN

    ed = np.quantile(gcf(bg), np.linspace(0, 1, NQ + 1))[1:-1]
    bq, tq = np.digitize(gcf(bg), ed), np.digitize(gcf(st), ed)
    rng = np.random.default_rng(seed)
    rel = 'recomputed from NC_003888.3.fna + jeong2016_all_tss.csv'
    spec = {'GCCGGC hexamer': (['GCCGGC'], 2539, 2883.8, 0.8804, '0.8478-0.9160', 0.7525),
            'AAGCCCG both strands': (['AAGCCCG', 'CGGGCTT'], 305, 322.0, 0.9471,
                                     '0.8424-1.0601', 1.0898),
            'GGC or GCC trinucleotide': (['GGC', 'GCC'], 142109, 146618.7, 0.9692,
                                         '0.9646-0.9739', 0.8675)}
    for name, (pats, c_obs, c_exp, c_r, c_ci, c_un) in spec.items():
        hits = np.sort(np.array([m.start() for p in pats
                                 for m in re.finditer(f'(?={p})', g)]))

        def cnt(s):
            return np.searchsorted(hits, s + WIN) - np.searchsorted(hits, s)

        o, bgc = cnt(st), cnt(bg)
        dens = np.array([bgc[bq == q].mean() for q in range(NQ)])[tq]
        exp_m, exp_p = float(dens.sum()), float(bgc.mean() * len(st))
        reps = np.empty(nboot)
        for i in range(nboot):
            k = rng.integers(0, len(o), len(o))
            reps[i] = o[k].sum() / dens[k].sum()
        lo, hi = np.percentile(reps, [2.5, 97.5])
        add('SuppFig 12', f'{name}: observed', c_obs, int(o.sum()),
            'OK' if int(o.sum()) == c_obs else 'MISMATCH', rel)
        add('SuppFig 12', f'{name}: GC-matched expected', c_exp, round(exp_m, 1),
            'OK' if close(exp_m, c_exp, 0.05) else 'MISMATCH', rel)
        add('SuppFig 12', f'{name}: ratio', c_r, round(o.sum() / exp_m, 4),
            'OK' if close(o.sum() / exp_m, c_r, 5e-5) else 'MISMATCH', rel)
        add('SuppFig 12', f'{name}: unmatched ratio', c_un, round(o.sum() / exp_p, 4),
            'OK' if close(o.sum() / exp_p, c_un, 5e-5) else 'MISMATCH', rel)
        add('SuppFig 12', f'{name}: 95% interval', c_ci, f'{lo:.4f}-{hi:.4f}',
            'METHOD', rel + f' (window bootstrap, {nboot} reps, seed {seed})')
    add('SuppFig 12', 'GC content, TSS +/-300 bp (%)', 69.54,
        round(100 * gcf(st).mean(), 2),
        'OK' if close(100 * gcf(st).mean(), 69.54, 0.005) else 'MISMATCH', rel)
    add('SuppFig 12', 'genome GC content (%)', 72.12,
        round(100 * float(cum[-1]) / glen, 2), 'OK', rel)
    add('SuppFig 12', 'TSS windows used', 3030, len(st),
        'OK' if len(st) == 3030 else 'MISMATCH', rel)


# ---- 補足図 10: core/arm の被覆深度と ori/ter -------------------------------
def suppfig10():
    pil = pathlib.Path('/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015')
    rel = 'pileup_canonical_261015 (per-library bed)'
    reps = {'T1': '1-1', 'T2': '2-1', 'T3': '3-2'}          # 凡例が使う代表ライブラリ
    claim_cov = {'T1': 2.26, 'T2': 1.13, 'T3': 1.19}
    claim_ot = {200: {'T1': 4.09, 'T2': 1.08, 'T3': 1.33},
                100: {'T1': 4.12, 'T2': 0.98, 'T3': 1.18},
                500: {'T1': 3.93, 'T2': 1.16, 'T3': 1.33}}
    CORE, ORIC, GLEN = (1_500_000, 7_170_000), 4_270_310, 8_667_507
    for tp, lib in reps.items():
        f = pil / f'{lib}_pileup.bed'
        if not f.exists():
            add('SuppFig 10', f'core/arm coverage ratio, {tp}', claim_cov[tp],
                'input missing', 'MISSING', str(f))
            continue
        # bedMethyl: 1 start, 9 Nvalid_cov, 15 Nfail。真深度 = Nvalid_cov + Nfail。
        # 同一位置に修飾コードごとの行があるため、位置で畳んで最大値を取る。
        d = pd.read_csv(f, sep='\t', header=None, usecols=[1, 9, 15],
                        names=['pos', 'nvalid', 'nfail'])
        d['depth'] = d.nvalid + d.nfail
        d = d.groupby('pos', as_index=False)['depth'].max()
        pos, cov = d.pos.to_numpy(), d.depth.to_numpy()
        incore = (pos >= CORE[0]) & (pos <= CORE[1])
        ratio = cov[incore].mean() / cov[~incore].mean()
        add('SuppFig 10', f'core/arm coverage ratio, {tp} ({lib})', claim_cov[tp],
            round(float(ratio), 2),
            'OK' if close(ratio, claim_cov[tp], 0.015) else 'MISMATCH', rel)
        for w, cl in claim_ot.items():
            hw = w * 1000 // 2
            o = cov[(pos >= ORIC - hw) & (pos <= ORIC + hw)].mean()
            ter = np.concatenate([cov[pos <= hw], cov[pos >= GLEN - hw]]).mean()
            add('SuppFig 10', f'ori/ter depth ratio, {tp}, {w} kb windows', cl[tp],
                round(float(o / ter), 2),
                'OK' if close(o / ter, cl[tp], 0.015) else 'NOT_REPRODUCIBLE', rel)
    add('SuppFig 10', 'method recovery: variants tried 2026-10-15',
        'legend says "computed from the modkit pileups"',
        'depth as Nvalid_cov+Nfail and as Nvalid_cov alone; canonical pileups '
        '(pileup_canonical_261015) and superseded pileups (analysis/pileup); terminus as the '
        'mean of both telomeres. None reproduces the T2/T3 ori/ter values; T1 core/arm differs '
        'on both pileup generations (canonical 2.31, superseded 1.75, claimed 2.26).',
        'NOT_REPRODUCIBLE', 'see recompute_261015/README.md')


if __name__ == '__main__':
    suppfig1()
    suppfig9()
    suppfig11()
    suppfig10()
    suppfig12()
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep='\t', index=False)
    print(df[['figure', 'quantity', 'claimed', 'recomputed', 'status']].to_string(index=False))
    print('\n' + df.status.value_counts().to_string())
    print(f'\n-> {OUT}')
