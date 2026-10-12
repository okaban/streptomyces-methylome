# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `260224_H7_CCGG_MTase_paradox_report.md` | 13 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | \| 1 \| **ZERO positional overlap** across all timepoints (Jaccard = 0.000 for T1-T2, T1-T3, T2-T3) \| Complet |
| `260224_H7_CCGG_MTase_paradox_report.md` | 14 | `core-to-arm` | BLOCKER-0 の物語語彙 | \| 2 \| **T1: 1,516 CCGG 4mC sites (core genome), T2: 486 (chromosomal arms), T3: 30 (scattered)** \| Dramatic |
| `260224_H7_CCGG_MTase_paradox_report.md` | 56 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | \| T1 intersect T2 \| **0** \| **0.000** \| |
| `260224_H7_CCGG_MTase_paradox_report.md` | 57 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | \| T1 intersect T3 \| **0** \| **0.000** \| |
| `260224_H7_CCGG_MTase_paradox_report.md` | 58 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | \| T2 intersect T3 \| **0** \| **0.000** \| |
| `260224_H7_CCGG_MTase_paradox_report.md` | 64 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | \| TGGCCGGC \| 1,289 \| 407 \| 21 \| 0 \| 0 \| |
| `260224_H7_CCGG_MTase_paradox_report.md` | 94 | `3.7%` | 同 T3 | TGGCCGGC (dominant motif) follows the same pattern: T1 mean 83.7% --> T2 78.4% --> T3 68.3%. |
