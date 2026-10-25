#!/usr/bin/env python3
"""AAGCCCG 6mA 鎖方向性: 宿主遺伝子クラスタリングを考慮した再推定。

BLIND_SPOT_AUDITOR 所見 M5 への応答。検定単位が「インスタンス」だったため、
1 遺伝子に複数インスタンスが乗ることによる非独立性が無視されていた。

出力
  SuppTable27b_strand_orientation_clustered_261021.tsv  時点 × 推定法
  strand_cluster_sensitivity_261021.tsv                 側方窓・GC 窓の感度
  strand_cluster_genecounts_261021.tsv                  1 遺伝子あたりインスタンス数分布

被説明変数   modified = (Nvalid+Nfail)>=30 かつ Nvalid>0 かつ Nmod/Nvalid>=0.5
主効果       sense = 1 (モチーフが宿主遺伝子のコード鎖) / 0 (鋳型鎖)
統制（既存踏襲・本スクリプトで明示的に定義）
  側方 9 塩基  モチーフ方向で上流 4 塩基 (-4..-1) + 下流 5 塩基 (+1..+5)、各位置を
               4 水準のカテゴリ変数として投入（計 27 ダミー）
  コドン位相   6mA 塩基の宿主 CDS 読み枠内位置 0/1/2、CDS 外 (rRNA 遺伝子) は -1 を
               独立水準として保持
  局所 GC      モチーフ左端を中心とする 501 bp 窓の GC 率（連続）
  ※ 原稿既報の調整済み OR（T1 2.641 [1.41-4.95] 等）を出したスクリプトは
    ディスク上に残っていないため、上記は本スクリプトで定義し直した仕様である。
    側方窓・GC 窓の選択に対する感度は sensitivity 表に出す。

実行: cd scripts/recompute_261015 && python3 strand_cluster_261021_analyse.py
"""
import pathlib
import subprocess
import tempfile
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

warnings.filterwarnings('ignore')
HERE = pathlib.Path(__file__).resolve().parent
INST = HERE / 'strand_instances_261021.tsv'
B_BOOT = 2000
SEED = 20261021
LAT = ['L-4', 'L-3', 'L-2', 'L-1', 'Lp1', 'Lp2', 'Lp3', 'Lp4', 'Lp5']
LAT_R = ['Lm4', 'Lm3', 'Lm2', 'Lm1', 'Lp1', 'Lp2', 'Lp3', 'Lp4', 'Lp5']
GC = 'gc501'
TPS = ('T1', 'T2', 'T3')

inst = pd.read_csv(INST, sep='\t')
ADJ = ('modified ~ sense + C(codon_phase) + ' + GC + ' + '
       + ' + '.join(f'Q("{c}")' for c in LAT))
UNADJ = 'modified ~ sense'


def subset(tp):
    d = inst[(inst.timepoint == tp) & (inst.measurable == 1)].copy()
    return d.reset_index(drop=True)


def counts(d):
    s, a = d[d.sense == 1], d[d.sense == 0]
    return dict(n_inst_sense=len(s), n_mod_sense=int(s.modified.sum()),
                n_inst_anti=len(a), n_mod_anti=int(a.modified.sum()),
                n_genes_sense=s.gene.nunique(), n_genes_anti=a.gene.nunique(),
                n_genes_total=d.gene.nunique(),
                occ_sense_pct=round(100 * s.modified.mean(), 4),
                occ_anti_pct=round(100 * a.modified.mean(), 4))


def row(tp, est, orr, lo, hi, p, d, note='', conv=None):
    r = dict(timepoint=tp, estimator=est,
             OR=None if orr is None else round(float(orr), 4),
             CI_low=None if lo is None else round(float(lo), 4),
             CI_high=None if hi is None else round(float(hi), 4),
             p=None if p is None else float(p))
    r.update(counts(d))
    r['optimiser_converged'] = conv
    r['note'] = note
    return r


def conv_of(m):
    try:
        return bool(m.mle_retvals.get('converged', True)) and np.isfinite(m.llf)
    except Exception:
        return None


# 重要: 側方 9 塩基を 4 水準カテゴリで入れると設計行列が 33 列になり、statsmodels
# 既定の newton では収束しないことがある。収束しなかった newton は llf=-inf を返し、
# 真の最適解 (OR≈2.4-3.0) ではなく OR≈1.4-2.1 の偽値を返す。そのため収束した
# optimiser を順に探し、収束フラグと採用 optimiser を必ず記録する。
OPTIMISERS = (('newton', 500), ('bfgs', 5000), ('lbfgs', 5000))


def fit_logit(formula, d, cov=None, groups=None):
    kw = {}
    if cov == 'cluster':
        kw = dict(cov_type='cluster', cov_kwds={'groups': groups})
    best = None
    for meth, mi in OPTIMISERS:
        try:
            m = smf.logit(formula, data=d).fit(disp=0, method=meth, maxiter=mi, **kw)
        except Exception:
            continue
        m.fit_method = meth
        if conv_of(m):
            best = m
            break
        if best is None and np.isfinite(m.llf):
            best = m
    if best is None:
        raise RuntimeError('どの optimiser でも収束せず')
    ci = best.conf_int().loc['sense']
    return (np.exp(best.params['sense']), np.exp(ci[0]), np.exp(ci[1]),
            best.pvalues['sense'], best)


def optimiser_check(formula, d, groups):
    """収束した optimiser 間で sense の OR が一致するかを見る。"""
    kw = dict(cov_type='cluster', cov_kwds={'groups': groups})
    ors = []
    for meth, mi in OPTIMISERS:
        try:
            m = smf.logit(formula, data=d).fit(disp=0, method=meth, maxiter=mi, **kw)
            if conv_of(m):
                ors.append(np.exp(m.params['sense']))
        except Exception:
            pass
    if len(ors) < 2:
        return f'収束 optimiser {len(ors)} 個のみ'
    return (f'収束 optimiser 間 OR 一致 {min(ors):.4f}-{max(ors):.4f} '
            f'(相対差 {100 * (max(ors) - min(ors)) / min(ors):.2f}%)')


# ---------------------------------------------------------------- bootstrap
def gene_bootstrap(d, formula, B=B_BOOT, seed=SEED):
    """遺伝子を復元抽出して調整済み logOR の分布を作る。"""
    mod = smf.logit(formula, data=d)
    y, X, names = mod.endog, mod.exog, list(mod.exog_names)
    j = names.index('sense')
    genes = d.gene.to_numpy()
    uniq, inv = np.unique(genes, return_inverse=True)
    idx_by_gene = [np.where(inv == k)[0] for k in range(len(uniq))]
    rng = np.random.default_rng(seed)
    out, fails = [], 0
    for _ in range(B):
        pick = rng.integers(0, len(uniq), len(uniq))
        rows = np.concatenate([idx_by_gene[k] for k in pick])
        Xb, yb = X[rows], y[rows]
        keep = (Xb.std(axis=0) > 0) | (np.arange(Xb.shape[1]) == names.index('Intercept'))
        jj = int(np.sum(keep[:j]))
        got = None
        for meth, mi in (('newton', 500), ('bfgs', 5000)):
            try:
                r = sm.Logit(yb, Xb[:, keep]).fit(disp=0, method=meth, maxiter=mi)
            except Exception:
                continue
            if (np.isfinite(r.llf) and r.mle_retvals.get('converged', True)
                    and np.isfinite(r.params[jj]) and abs(r.params[jj]) < 10):
                got = r.params[jj]
                break
        if got is None:
            fails += 1
        else:
            out.append(got)
    return np.array(out), fails


def boot_summary(b, fails, B=B_BOOT):
    lo, hi = np.percentile(b, [2.5, 97.5])
    frac = min((b <= 0).mean(), (b >= 0).mean())
    p = max(2 * frac, 1.0 / (len(b) + 1))
    return np.exp(np.median(b)), np.exp(lo), np.exp(hi), p, fails


# ------------------------------------------------------------------- glmer
GLMER_R = r'''
suppressMessages(library(lme4))
a <- commandArgs(trailingOnly=TRUE)
d <- read.delim(a[1])
for (c in c("codon_phase", %s)) d[[c]] <- factor(d[[c]])
f <- as.formula(paste("modified ~ sense + codon_phase + gc501 +",
                      paste(c(%s), collapse=" + "), "+ (1|gene)"))
m <- try(glmer(f, data=d, family=binomial,
               control=glmerControl(optimizer="bobyqa",
                                    optCtrl=list(maxfun=200000))), silent=TRUE)
if (inherits(m, "try-error")) { cat("FAIL\t", as.character(m), "\n"); quit(status=0) }
s <- summary(m)
co <- s$coefficients["sense", ]
vc <- as.data.frame(VarCorr(m))$vcov[1]
sing <- isSingular(m, tol=1e-4)
warn <- paste(unlist(m@optinfo$conv$lme4$messages), collapse="; ")
cat(sprintf("OK\t%%.10f\t%%.10f\t%%.10g\t%%.10f\t%%s\t%%s\n",
            co[1], co[2], co[4], vc, sing, warn))
'''


def glmer(d):
    cols = ['modified', 'sense', 'codon_phase', GC, 'gene']
    e = d[cols + LAT].rename(columns=dict(zip(LAT, LAT_R)))
    with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as fh:
        e.to_csv(fh.name, sep='\t', index=False)
        path = fh.name
    quoted = ', '.join(f'"{c}"' for c in LAT_R)
    src = GLMER_R % (quoted, quoted)
    with tempfile.NamedTemporaryFile('w', suffix='.R', delete=False) as fh:
        fh.write(src)
        rpath = fh.name
    try:
        p = subprocess.run(['Rscript', '--vanilla', rpath, path],
                           capture_output=True, text=True, timeout=3600)
    except Exception as ex:
        return None, f'glmer 実行不可: {ex}'
    line = [l for l in p.stdout.splitlines() if l.startswith(('OK', 'FAIL'))]
    if not line:
        return None, f'glmer 出力なし: {(p.stderr or "")[-200:]}'
    f = line[0].split('\t')
    if f[0] == 'FAIL':
        return None, 'glmer 収束せず: ' + f[1][:160].replace('\n', ' ')
    beta, se, pv, vc, sing, warn = (float(f[1]), float(f[2]), float(f[3]),
                                    float(f[4]), f[5].strip(), f[6].strip())
    note = f'gene 切片分散={vc:.4f}; singular={sing}'
    if warn:
        note += f'; {warn[:120]}'
    return (np.exp(beta), np.exp(beta - 1.96 * se), np.exp(beta + 1.96 * se), pv), note


# ----------------------------------------------------------- main per-timepoint
rows, sens, gcnt = [], [], []
for tp in TPS:
    d = subset(tp)
    g = d.gene.to_numpy()

    # 1) インスタンス単位（既報の再現）
    orr, lo, hi, p, m = fit_logit(UNADJ, d)
    rows.append(row(tp, 'A_unadjusted_instance_level', orr, lo, hi, p, d,
                    '既報 SuppTable27 の再現（素の OR）', conv_of(m)))

    # 2) 未調整 + 遺伝子クラスタ頑健 SE
    orr, lo, hi, p, m = fit_logit(UNADJ, d, 'cluster', g)
    rows.append(row(tp, 'B_unadjusted_cluster_robust_gene', orr, lo, hi, p, d,
                    '遺伝子クラスタ頑健 SE', conv_of(m)))

    # 3) 調整済み・インスタンス単位（素の SE）
    orr, lo, hi, p, m = fit_logit(ADJ, d)
    rows.append(row(tp, 'C_adjusted_instance_level', orr, lo, hi, p, d,
                    '側方9塩基+コドン位相+GC501、クラスタリング未考慮', conv_of(m)))

    # 4) 調整済み + 遺伝子クラスタ頑健 SE
    orr, lo, hi, p, mc = fit_logit(ADJ, d, 'cluster', g)
    rows.append(row(tp, 'D_adjusted_cluster_robust_gene', orr, lo, hi, p, d,
                    f'主解析: 遺伝子クラスタ頑健 SE; optimiser={mc.fit_method}; '
                    + optimiser_check(ADJ, d, g), conv_of(mc)))

    # 5) GEE (exchangeable, gene)
    try:
        gee = smf.gee(ADJ, groups='gene', data=d,
                      family=sm.families.Binomial(),
                      cov_struct=sm.cov_struct.Exchangeable()).fit()
        ci = gee.conf_int().loc['sense']
        alpha = gee.cov_struct.dep_params
        rows.append(row(tp, 'E_adjusted_GEE_exchangeable_gene',
                        np.exp(gee.params['sense']), np.exp(ci[0]), np.exp(ci[1]),
                        gee.pvalues['sense'], d,
                        f'作業相関 exchangeable, alpha={float(alpha):.4f}'))
    except Exception as ex:
        rows.append(row(tp, 'E_adjusted_GEE_exchangeable_gene', None, None, None,
                        None, d, f'収束せず: {str(ex)[:120]}'))

    # 6) 混合ロジット（遺伝子変量切片）
    res, note = glmer(d)
    if res is None:
        rows.append(row(tp, 'F_adjusted_mixed_logit_gene_RI', None, None, None,
                        None, d, note))
    else:
        rows.append(row(tp, 'F_adjusted_mixed_logit_gene_RI', *res, d, note))

    # 7) 遺伝子ブートストラップ
    b, fails = gene_bootstrap(d, ADJ)
    orr, lo, hi, p, fails = boot_summary(b, fails)
    rows.append(row(tp, 'G_adjusted_gene_bootstrap_2000', orr, lo, hi, p, d,
                    f'遺伝子を復元抽出 B={B_BOOT}, 成功 {len(b)}, 失敗 {fails}; '
                    f'パーセンタイル CI'))

    # 8) 遺伝子単位に集約した対応あり検定（両方向を持つ遺伝子のみ）
    agg = (d.groupby(['gene', 'sense']).modified.agg(['sum', 'size'])
             .reset_index())
    piv = agg.pivot(index='gene', columns='sense', values=['sum', 'size'])
    both = piv.dropna()
    rs = (both[('sum', 1)] / both[('size', 1)]).to_numpy()
    ra = (both[('sum', 0)] / both[('size', 0)]).to_numpy()
    diff = rs - ra
    nz = diff[diff != 0]
    n_pos, n_neg = int((nz > 0).sum()), int((nz < 0).sum())
    p_sign = stats.binomtest(n_pos, len(nz), 0.5).pvalue if len(nz) else None
    rows.append(row(tp, 'H_gene_paired_sign_test', None, None, None, p_sign, d,
                    f'両方向を持つ遺伝子 n={len(both)}、非同点 {len(nz)} '
                    f'(コード鎖>鋳型鎖 {n_pos} / 逆 {n_neg})、'
                    f'中央値差={np.median(diff):.4f}'))
    if len(nz):
        try:
            w = stats.wilcoxon(rs, ra, zero_method='wilcox')
            pw = w.pvalue
        except Exception:
            pw = None
    else:
        pw = None
    rows.append(row(tp, 'I_gene_paired_wilcoxon', None, None, None, pw, d,
                    f'対応あり Wilcoxon 符号順位、n={len(both)} 遺伝子、'
                    f'平均差={np.mean(diff):.4f}'))

    # 9) 遺伝子単位の非対応比較（全遺伝子、遺伝子ごとの修飾率）
    gs = (d[d.sense == 1].groupby('gene').modified.mean().to_numpy())
    ga = (d[d.sense == 0].groupby('gene').modified.mean().to_numpy())
    u = stats.mannwhitneyu(gs, ga, alternative='two-sided')
    rows.append(row(tp, 'J_gene_level_mannwhitney', None, None, None, u.pvalue, d,
                    f'遺伝子を 1 単位にした修飾率の非対応比較 '
                    f'(コード鎖 {len(gs)} 遺伝子 平均 {gs.mean():.4f} / '
                    f'鋳型鎖 {len(ga)} 遺伝子 平均 {ga.mean():.4f})'))

    # 10) 遺伝子単位に 1 行へ集約した未調整 OR（遺伝子を単位にした 2x2）
    #     各遺伝子×方向を「1 つでも修飾あり」で 2 値化
    any_s = d[d.sense == 1].groupby('gene').modified.max()
    any_a = d[d.sense == 0].groupby('gene').modified.max()
    tab = np.array([[int(any_s.sum()), int((1 - any_s).sum())],
                    [int(any_a.sum()), int((1 - any_a).sum())]])
    orres = stats.fisher_exact(tab)
    ct = sm.stats.Table2x2(tab)
    cilo, cihi = ct.oddsratio_confint()
    rows.append(row(tp, 'K_gene_unit_2x2_any_modified', ct.oddsratio, cilo, cihi,
                    orres.pvalue, d,
                    f'遺伝子単位 2x2: コード鎖 {tab[0,0]}/{tab[0].sum()} 遺伝子で '
                    f'1 つ以上修飾 vs 鋳型鎖 {tab[1,0]}/{tab[1].sum()}（Fisher）'))

    # ---- 1 遺伝子あたりインスタンス数分布 ----
    for v, lab in ((1, 'coding_strand'), (0, 'template_strand')):
        per = d[d.sense == v].groupby('gene').size()
        gcnt.append(dict(timepoint=tp, orientation=lab, n_instances=int(per.sum()),
                         n_genes=int(per.size), median_per_gene=float(per.median()),
                         mean_per_gene=round(float(per.mean()), 4),
                         max_per_gene=int(per.max()),
                         genes_with_1=int((per == 1).sum()),
                         genes_with_2=int((per == 2).sum()),
                         genes_with_3plus=int((per >= 3).sum())))

    # ---- 感度: 側方窓 × GC 窓 ----
    for up, down in ((4, 5), (5, 4), (3, 6), (6, 3), (2, 7), (7, 2), (9, 0), (0, 9)):
        cols = ([f'L-{k}' for k in range(up, 0, -1)]
                + [f'Lp{k}' for k in range(1, down + 1)])
        missing = [c for c in cols if c not in d.columns]
        if missing:
            continue
        for w in ('gc101', 'gc201', 'gc501', 'gc1001'):
            f = ('modified ~ sense + C(codon_phase) + ' + w
                 + (' + ' + ' + '.join(f'Q("{c}")' for c in cols) if cols else ''))
            try:
                o, l, h, pv, mm = fit_logit(f, d, 'cluster', g)
                sens.append(dict(timepoint=tp, lateral_up=up, lateral_down=down,
                                 gc_window=w, estimator='cluster_robust_gene',
                                 OR=round(o, 4), CI_low=round(l, 4),
                                 CI_high=round(h, 4), p=float(pv),
                                 converged=conv_of(mm),
                                 fit_method=getattr(mm, 'fit_method', None), note=''))
            except Exception as ex:
                sens.append(dict(timepoint=tp, lateral_up=up, lateral_down=down,
                                 gc_window=w, estimator='cluster_robust_gene',
                                 OR=None, CI_low=None, CI_high=None, p=None,
                                 converged=False, fit_method=None,
                                 note=str(ex)[:80]))
    print(f'{tp} done', flush=True)

out = pd.DataFrame(rows)
out.to_csv(HERE / 'SuppTable27b_strand_orientation_clustered_261021.tsv',
           sep='\t', index=False)
pd.DataFrame(sens).to_csv(HERE / 'strand_cluster_sensitivity_261021.tsv',
                          sep='\t', index=False)
pd.DataFrame(gcnt).to_csv(HERE / 'strand_cluster_genecounts_261021.tsv',
                          sep='\t', index=False)
print(out[['timepoint', 'estimator', 'OR', 'CI_low', 'CI_high', 'p']].to_string(index=False))
