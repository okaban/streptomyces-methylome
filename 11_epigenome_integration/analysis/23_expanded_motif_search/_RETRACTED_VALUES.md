# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `expanded_motif_catalog.csv` | 2 | `8bp palindrome` | TGGCCGGC は回文ではない。23_expanded_motif_search/expanded_motif_catalog.csv の誤記 | 1,CCGG/GGCCGG/TGGCCGGC,4mC,2034,175,75.9%,5.4%,YES (MspI/HpaII-type),No,Dominant 4mC motif; true recognition c |
| `integration_comprehensive_table.csv` | 2 | `C3/C5` | 同上 | 1,AAGCCCG,4mC + 6mA (DUAL),973 (4mC) + 360 (6mA),"36.3% of 4mC, 11.2% of 6mA","MEME-1 (6mA), STREME (4mC resid |
