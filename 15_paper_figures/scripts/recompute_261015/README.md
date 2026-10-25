# recompute_261015 — 生成元スクリプトが無い補足図の数値照合

補足図 1・9・10・11・12 は出荷 PNG のみで、ツリー内に描画スクリプトが無い
（`deposition/figure_producer_map.tsv` にも 4 件が no match として記録済み）。
図を描き直す代わりに、凡例の数値を出所まで辿って照合する。著者指示 2026-10-15。

    python3 recompute_261015/recompute_producerless.py   # 5 図まとめて → producerless_verification.tsv
    python3 recompute_261015/recompute_suppfigure12.py   # 補足図 12 単独（詳細版）

PNG には一切触れない。

## 結果（2026-10-15）

| 図 | 判定 |
|---|---|
| 補足図 1 | 7 値すべて再現。出所 `87_site_coloc_stats_C4/tables/site_coloc_2x2_C4.tsv` |
| 補足図 9 | 6 値すべて再現。出所 `78_reviewer_robustness/tables/P1_1_crosstalk_rates.tsv`（計算スクリプト `P1_1_crosstalk_test.py` も実在） |
| 補足図 11 | 統計値 3 件は再現するが、五分位検定の outcome が撤回済みの Exposed 分類（`is_exposed_canonical(62)`、n = 1,017）。凡例は「認識配列から 293 bp 以内かどうか」と書いており一致しない |
| 補足図 12 | 点推定 12 値すべて再現（ゲノムと TSS だけで完結）。区間は窓ブートストラップで、最大カウントの GGC/GCC が小数第 4 位まで一致するため手法はこれで確定。小さいカウントの 2 件は Monte-Carlo 誤差の範囲 |
| 補足図 10 | **再現せず。** 下記の変種をすべて試したが T2/T3 の ori/ter が合わない |

## 補足図 10 で試した変種（繰り返さないこと）

- 深度の定義: `Nvalid_cov + Nfail`（真深度）と `Nvalid_cov` のみ → 差は 0.01 未満で原因でない
- パイルアップ世代: `pileup_canonical_261015` と旧 `analysis/pileup` の両方
- terminus の取り方: 両テロメアの平均（線状染色体のため）

core/arm は正典 2.31 / 旧 1.75 に対し凡例 2.26。ori/ter 200 kb は T2 が 0.90 対 1.08、
T3 が 0.94 対 1.33 で、丸めでは説明できない。旧 `1-1_pileup.bed` は先頭 100 kb に
行が無く ori/ter が定義できない。

**結論**: 補足図 10 の被覆深度比は、凡例が書いている方法では復元できない。
再構築するか、凡例の記述を実際の方法に合わせて訂正する必要がある。
