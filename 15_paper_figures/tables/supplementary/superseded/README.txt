2026-10-21 退避。現行の凡例（Supplementary-Materials-List.md）から参照されていないファイル。
- ST1-ST8 系 12 点: Exposed/Shielded 区分 (62/955)、57 遺伝子 bloc、コア濃縮 83%、占有率 75.9%、
  REBASE 18/82 など撤回済みの値を含む。投稿パッケージに混入させないため退避。
- SuppTable27_occupancy_ladder_261019.tsv / SuppTable28_occupancy_ladder_261020.tsv:
  SuppTable28_occupancy_ladder_261021.tsv（4mC 列を追加、生成元 build_supptable28_261021.py）に置換。
- SuppTable25_per_library_depth_and_sites.csv: 中身は Supplementary Table 15 の凡例どおりの内容で、
  番号が誤っていた。SuppTable15_per_library_depth_and_sites.tsv として作り直した。
削除ではない。必要なら戻せる。

2026-10-21 (supplementary body completion):
  SuppTable19a_eligibility_for_promoter_test_261019.tsv
    - held the gene-set construction steps (step/criterion/n/source), not the
      columns the legend declares. Replaced by ..._261021.tsv, which gives
      window x position set eligibility counts.
  SuppTable19b_detectable_effect_sizes_261019.tsv
    - the per-read co-modification row carried a Japanese internal note saying
      the marginals had not been recomputed. Replaced by ..._261021.tsv.
  SuppTable16_statistical_test_summary_261016.tsv
    - n / correction / p_corrected were empty in all 29 rows. Replaced by
      ..._261021.tsv.
