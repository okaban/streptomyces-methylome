# 撤回値を含むファイル — 自動生成の注記

- 生成日: 2026-09-29
- 判定元: `/Users/okaban/obsidian/Research/rna-seq/Writing/canon_values.json` の `retracted`（75 トークン）
- 検証記録: `s23_hit_verification.csv`（ヒット 1 件ごとに「散文・表中の参照」か「数値列での偶然一致」かを判定）
- **ファイルは移動も編集もしていません。** スクリプトからのパス参照を壊さないため、注記はこのサイドカーのみです。

以下の行は撤回済みの値を含みます。公開・引用・図の再生成の前に該当行を確認してください。

| ファイル | 行 | 撤回トークン | 撤回理由 / 正しい値 | 文脈 |
|---|---|---|---|---|
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 1 | `57 Exposed` | 旧 Exposed 数 | # H34: Temporal Dynamics and Hierarchical Structure of 57 Exposed Transcription Factors |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 12 | `two antagonistic blocs` | 二プログラム枠組み、撤回 | H32 identified that the 57 exposed TFs split into two antagonistic blocs: an **activation bloc** (~35 genes, c |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 118 | `0.813` | TGGCCGGC のゲノム内マルコフ O/E。次数 1–5 で 1.502 / 0.813 / 0.869 / 0.940 / 1.118 と安定せず、分野の基準 0.78 も全次数で上回る。TGGCCGGC の回避は成 | - **Temporal lag (exposed - shielded):** mean = 0.067, Wilcoxon p = 0.813 -- no systematic ordering |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 136 | `0.000` | first-appearance 表の Jaccard。正準 0.700 | - All pairs (early x late responders): rho = 0.000, p = 1.0 -- no cross-module prediction |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 165 | `0.813` | TGGCCGGC のゲノム内マルコフ O/E。次数 1–5 で 1.502 / 0.813 / 0.869 / 0.940 / 1.118 と安定せず、分野の基準 0.78 も全次数で上回る。TGGCCGGC の回避は成 | \| TCS temporal lag (Wilcoxon) \| W=6.0 \| 0.813 \| No systematic exposed-first ordering \| |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 189 | `Gatekeeper` | 四層モデル、撤回 | ## Model Update: Gatekeeper Model v3 Implications |
| `260227_H34_temporal_dynamics_exposed_TF_report.md` | 191 | `Gatekeeper` | 四層モデル、撤回 | The synchronous activation of both blocs at T1->T2 refines the Gatekeeper Model: |
| `260504_v_defense_MTase_expression_dynamics_report.md` | 37 | `palindromic GCCGGC` | GCCGGC は認識配列ではない。認識配列 TGGCCGGC は非回文（相補鎖 GCCGGCCA） | GCCGGC-m4C density: count of called 4mC sites at GCCGGC motifs (TGGCCGGC / GGCCGG context, `final_motif` in ce |
| `260504_v_defense_MTase_expression_dynamics_report.md` | 83 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | \| T1 \| 1,289 \| 83.7 \| **4.10 %** \| |
| `260504_v_defense_MTase_expression_dynamics_report.md` | 140 | `1,289` | T1 GCCGGC 4mC 部位数。Nvalid_cov ≥10 の床がメチル化状態に対する床として働いていた。正: 1,886（緩和規則） | 3. **Disambiguate the T2 detection-floor effect.** Quantify how much of the T1→T2 site-count drop (1,289 → 407 |
