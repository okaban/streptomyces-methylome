# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `README.md` | 5 | `138,440` | 旧 census の site OR。訂正後は飽和（Haldane 補正 ≈ 6.1e5） | Replaces 68_or_permutation (OR 138,440; 244 in / 1,090 not; universe 1,238,215 hard-coded). Same 2x2 / |
| `README.md` | 21 | `4.84` | per-read OR、read 方向バグ＋ Simpson。修正版 OR 1.04、MH 1.07 | 12,771 expected under independence at ≥ 0.5 (fold 1.97; OR 4.84 [4.73, 4.95]). |
| `README.md` | 23 | `418` | T1 AAGCCCG A₀/A₁ 6mA 部位数。正: 416（差は 2 部位。通過率の高い 6mA は床の影響をほぼ受けず、これが補正の内部対照） | 6mA occupancy at AAGCCCG, T1 (corrected): 418 canonical 6mA sites over 2,668 A0/A1 positions (15.7%); 366 of 1 |
| `README.md` | 24 | `698` | T1 AAGCCCG C₄ 4mC 部位数。正: 962 | instances (27.4%) carry 6mA at A0 or A1; 698 of 1,334 C4 positions (52.3%) are canonical 4mC. The manuscript's |
| `README.md` | 25 | `31.6%` | 分母不明の旧占有率 | '447 6mA sites at T1, 31.6% of genomic AAGCCCG positions' has no traceable denominator (447/0.316 = 1,415). |
| `README.md` | 28 | `418` | T1 AAGCCCG A₀/A₁ 6mA 部位数。正: 416（差は 2 部位。通過率の高い 6mA は床の影響をほぼ受けず、これが補正の内部対照） | The AAGCCCG motif also carried N6-methyladenine (6mA) at its two adenines (A₀ and A₁): at T1, 418 of the 2,668 |
| `README.md` | 30 | `698` | T1 AAGCCCG C₄ 4mC 部位数。正: 962 | at A₀ or A₁, whereas the 4mC call lay at the central cytosine (C₄; 698 of 1,334 instances, 52.3%). At the site |
| `README.md` | 40 | `47% of reads` | 同上 | covering A₀/A₁ and C₄ at T1) showed 6mA on 47% of reads carrying 4mC at C₄ versus 24% of all reads (OR 4.84, 9 |
| `README.md` | 40 | `4.84` | per-read OR、read 方向バグ＋ Simpson。修正版 OR 1.04、MH 1.07 | covering A₀/A₁ and C₄ at T1) showed 6mA on 47% of reads carrying 4mC at C₄ versus 24% of all reads (OR 4.84, 9 |
