# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `README.md` | 8 | `4.84` | per-read OR、read 方向バグ＋ Simpson。修正版 OR 1.04、MH 1.07 | - pooled OR reproduces 79_: 4.84 / 3.17 / 2.40 at p≥0.5 / 0.75 / 0.9. |
| `README.md` | 9 | `4.84` | per-read OR、read 方向バグ＋ Simpson。修正版 OR 1.04、MH 1.07 | - per-replicate 1-1/1-2/1-3: 4.82/4.84/4.85 (0.5); 3.11/3.21/3.19 (0.75); 2.28/2.51/2.42 (0.9). |
| `README.md` | 14 | `698` | T1 AAGCCCG C₄ 4mC 部位数。正: 962 | - `tables/REFA02a_site_level_within_motif_fisher_T1.tsv`: instance level, N=1,334: 698 4mC(C4) × 366 6mA(A0\|A |
| `README.md` | 14 | `418` | T1 AAGCCCG A₀/A₁ 6mA 部位数。正: 416（差は 2 部位。通過率の高い 6mA は床の影響をほぼ受けず、これが補正の内部対照） | - `tables/REFA02a_site_level_within_motif_fisher_T1.tsv`: instance level, N=1,334: 698 4mC(C4) × 366 6mA(A0\|A |
| `README.md` | 17 | `−0.92` | 鎖間 Pearson r。文脈の有無で決まる量であってメチル化の性質ではない | `tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): |
| `README.md` | 17 | `1,493` | ヘミメチル化の分母（レプリケート 1-3 の occurrence）。情報を持つ単位は両鎖に文脈がある 55 二本鎖のみ | `tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): |
| `README.md` | 17 | `1,493` | ヘミメチル化の分母（レプリケート 1-3 の occurrence）。情報を持つ単位は両鎖に文脈がある 55 二本鎖のみ | `tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): |
| `README.md` | 17 | `−0.92` | 鎖間 Pearson r。文脈の有無で決まる量であってメチル化の性質ではない | `tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): |
| `README.md` | 17 | `99.8` | hemi-methylation の出所不明値 | `tables/REFA06_GCCGGC_strand_symmetry_pileup.tsv` (8 pileups). Among methylated instances (≥50% on ≥1 strand): |
