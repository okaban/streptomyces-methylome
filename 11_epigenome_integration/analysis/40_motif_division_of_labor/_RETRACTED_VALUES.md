# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `260226_H17_motif_division_of_labor_report.md` | 26 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | - **GCCGGC-proximal (T1):** Mapped 1,289 T1 GCCGGC sites to protein-coding genes within 2 kb using GFF annotat |
| `260226_H17_motif_division_of_labor_report.md` | 27 | `451` | 同 T2。正: 450 | - **Overlap:** 451 genes proximal to both motifs ("dual-targeted") |
| `260226_H17_motif_division_of_labor_report.md` | 108 | `451` | 同 T2。正: 450 | - 451 genes are proximal to both AAGCCCG and GCCGGC sites (Jaccard index = 0.107) |
| `260226_H17_motif_division_of_labor_report.md` | 111 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | - The dual-targeted expression pattern resembles GCCGGC-only more than AAGCCCG-only, consistent with the much  |
| `260226_H17_motif_division_of_labor_report.md` | 119 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | 1. **Compositional confound:** GCCGGC has 1,289 T1 sites covering 3,695 genes (46% of protein-coding genes). T |
| `260226_H17_motif_division_of_labor_report.md` | 134 | `Gatekeeper` | 四層モデル、撤回 | This is consistent with the Gatekeeper Model v2 (H11), where methylation acts as a constraint on gene expressi |
| `260226_H17_motif_division_of_labor_report.md` | 138 | `Gatekeeper` | 四層モデル、撤回 | **REJECTED** -- The hypothesis that GCCGGC targets undergo developmental upregulation while AAGCCCG targets re |
