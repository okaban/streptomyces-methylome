# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `260226_H25_TSS_methylation_gradient_report.md` | 22 | `1,987` | T1 全 4mC 部位数。正: 2,758 | - **Methylation sites (T1)**: All methylation (3,918), All 4mC (1,987), All 6mA (1,934), GCCGGC 4mC (1,289), A |
| `260226_H25_TSS_methylation_gradient_report.md` | 22 | `1,934` | T1 全 6mA 部位数。正: 2,003 | - **Methylation sites (T1)**: All methylation (3,918), All 4mC (1,987), All 6mA (1,934), GCCGGC 4mC (1,289), A |
| `260226_H25_TSS_methylation_gradient_report.md` | 22 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | - **Methylation sites (T1)**: All methylation (3,918), All 4mC (1,987), All 6mA (1,934), GCCGGC 4mC (1,289), A |
| `260226_H25_TSS_methylation_gradient_report.md` | 63 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | \| AAGCCCG 6mA \| n.d. \| -- \| 0.000 \| +300 bp \| -0.011 \| |
| `260226_H25_TSS_methylation_gradient_report.md` | 221 | `Gatekeeper` | 四層モデル、撤回 | ## Connection to Gatekeeper Model |
| `260226_H25_TSS_methylation_gradient_report.md` | 223 | `Gatekeeper` | 四層モデル、撤回 | This analysis provides the first spatially resolved view of the Layer 2 (Protection/Depletion) component of th |
