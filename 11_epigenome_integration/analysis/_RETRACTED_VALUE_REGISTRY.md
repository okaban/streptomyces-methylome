# 撤回値レジストリ — 2026-09-29

`canon_values.json` の撤回トークン 75 件を解析ツリー全体（821 ファイル / 960 MB）に照合し、
ヒット 1 件ごとに検証した結果です。**ファイルの移動・編集は一切していません。**

| 判定 | ヒット | ファイル | 意味 |
|---|---|---|---|
| GENUINE | 325 | 77 | 散文・要約表の中で撤回値を実際に参照している。注記対象 |
| AUDIT_RECORD | 184 | 17 | 監査・レビュー記録。撤回前の値を引用するのが正しい。凍結 |
| COINCIDENCE | 698 | 104 | データ列・外部ツール出力での数値の偶然一致。対応不要 |

## 偶然一致の判定根拠
- 3 桁小数トークン（0.813 など）が、20 種類以上の相異なる数値を持つ数値列のセル値そのものである場合
- 整数トークン（443, 451 など）がデータ表の中にある場合（座標・長さ・カウントと区別できない）
- 指数表記（4.84e-44）、桁区切り（30,451）、長い数値の部分一致（443.3018）
- `antismash_output/` 配下（antiSMASH の外部出力）、5 MB 超の生 pileup 抽出

### 当初台帳からの主な訂正
`21_genuswide_motif_conservation/motif_site_density_genuswide.csv` は当初「隔離」対象（33 ヒット）でしたが、
**33 件すべてが 833 ゲノム × モチーフの O/E 値列における偶然一致**でした。このテーブルは復活した
AAGCCCG 回避主張（中央値 O/E 0.703、610/833 が 0.75 未満）の生証拠であり、隔離も列削除も不要です。
同ディレクトリで実際に撤回値を含むのは `GENUSWIDE_MOTIF_CONSERVATION_REPORT.md`（`19/82`、正: 14 of 82 species）と
`analysis_log.txt` の 3 行のみです。

## 注記対象ディレクトリ（`_RETRACTED_VALUES.md` を設置）

| ディレクトリ | ファイル | ヒット |
|---|---|---|
| `.` | 3 | 34 |
| `92_d6_5mC_misassignment_control` | 2 | 25 |
| `57_temporal_dynamics_exposed_TF/results` | 3 | 18 |
| `52_shielded_exposed_boundary` | 1 | 13 |
| `42_GCCGGC_temporal_derepression` | 1 | 13 |
| `60_coexpression_regulon_prediction` | 1 | 12 |
| `53_TSS_sequence_determinants` | 1 | 11 |
| `57_temporal_dynamics_exposed_TF` | 2 | 10 |
| `59_exposed_TF_conservation` | 1 | 10 |
| `34_gatekeeper_model_v2` | 1 | 10 |
| `87_site_coloc_stats_C4` | 1 | 9 |
| `91_comod_stratified_and_hemi` | 1 | 9 |
| `37_defense_island_GCCGGC` | 1 | 8 |
| `39_GCCGGC_MTase_reverse_ID` | 1 | 8 |
| `38_cross_motif_regulatory_avoidance` | 1 | 8 |
| `50_coordinated_regulators_protection` | 1 | 7 |
| `40_motif_division_of_labor` | 1 | 7 |
| `30_CCGG_MTase_paradox` | 1 | 7 |
| `32_CCGG_5mC_misclassification` | 1 | 6 |
| `79_comod_full_denominator` | 4 | 6 |
| `48_TSS_methylation_gradient` | 1 | 6 |
| `78_reviewer_robustness` | 2 | 5 |
| `41_GCCGGC_dose_response` | 1 | 5 |
| `01_integration` | 2 | 5 |
| `73_g2_correction_qc` | 1 | 4 |
| `56_exposed_TF_neighborhood` | 1 | 4 |
| `46_category_specificity_avoidance` | 1 | 4 |
| `52_shielded_exposed_boundary/tables` | 1 | 3 |
| `47_BGC_methylation_geographic_test` | 1 | 3 |
| `45_sequence_level_motif_depletion` | 1 | 3 |
| `07_motif_analysis` | 1 | 3 |
| `44_AAGCCCG_temporal_causality` | 1 | 3 |
| `21_genuswide_motif_conservation` | 2 | 3 |
| `34_gatekeeper_model_v2/tables` | 3 | 3 |
| `58b_exposed_TF_functional_prediction` | 1 | 2 |
| `62_GO_KEGG_enrichment` | 1 | 2 |
| `70_methylation_dynamics/tables` | 1 | 2 |
| `23_expanded_motif_search` | 2 | 2 |
| `49_TF_BS_methylation_protection` | 1 | 2 |
| `88_occupancy_series_and_CE` | 1 | 2 |
| `63_window_methylation` | 1 | 2 |
| `73_SC_RS17645` | 1 | 2 |
| `32_CCGG_5mC_misclassification/tables` | 2 | 2 |
| `73_layer_overlap` | 1 | 2 |
| `58_AAGCCCG_exposed_TF_causal` | 1 | 2 |
| `43_regulatory_avoidance_geographic_test` | 1 | 2 |
| `55_exposed_regulatory_module` | 1 | 2 |
| `51_exposed_regulators_characteristics` | 1 | 1 |
| `89_supp_deliverables` | 1 | 1 |
| `86_geneset_counts_B10` | 1 | 1 |
| `31_coordinated_regulators_characterization` | 1 | 1 |
| `68_or_permutation` | 1 | 1 |
| `76_FIRE_methylation_crossref` | 1 | 1 |
| `75_tss_downstream_protection` | 1 | 1 |
| `33_GCCGGC_RM_identification` | 1 | 1 |
| `70_methylation_dynamics` | 1 | 1 |
| `64_timepoint_TSS` | 1 | 1 |
| `35_GCCGGC_MTase_BLAST` | 1 | 1 |
| `36_AAGCCCG_distribution` | 1 | 1 |
| `54_exposed_TF_downstream_network` | 1 | 1 |
| `94_figure_restructure_260921` | 1 | 1 |

## 監査記録として凍結（`_AUDIT_RECORD_FROZEN.md` を設置）

| ディレクトリ | ファイル | ヒット |
|---|---|---|
| `100_opportunity_feasibility` | 1 | 24 |
| `101_discussion_insertions` | 1 | 3 |
| `102_post_edit_parity` | 1 | 6 |
| `105_referee_B` | 1 | 32 |
| `90_per_timepoint_census_audit` | 1 | 25 |
| `90_per_timepoint_census_audit/tables` | 2 | 22 |
| `93_enjp_parity_audit_260921/tables` | 5 | 27 |
| `95_enjp_parity_round2` | 1 | 4 |
| `96_supp_legend_vs_data` | 1 | 10 |
| `97_opportunity_scan` | 1 | 8 |
| `99_figure_forensics` | 1 | 4 |
| `out_nar_revision_260716` | 1 | 19 |
