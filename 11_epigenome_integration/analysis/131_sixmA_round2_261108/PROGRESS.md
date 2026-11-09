2026-11-08 着手。skill 2 本・前期報告・open_questions/findings 確認。A の 3 値の仕様差を文書から同定、事前登録を確定。
2026-11-08 A 完了。3 値を同一部位集合(1,156/1,158/1,157)・同一転帰で再現し分解。
  S8 2.8943/2.9597/4.1838 と 129 レーン 1.9853/2.1493/2.2038 を両方 4 桁再現。
  差の主因は codon_phase（leading ではない: 投入で OR 変化 <2.5%）。事前登録 P-A1 は棄却。
  codon_phase は sense のほぼ再符号化（Cramer V 0.811、90.4% 予測）。
  表: tables/A_orientation_model_ladder_261108.tsv, A_orientation_within_codonphase_261108.tsv, A_orientation_MH_261108.tsv
2026-11-08 B 主要部完了。正準部位集合を per-library cache から完全再現 (1,947/2,150/2,308、AAGCCCG 1,334)。
  3 ライブラリ独立コール: 非モチーフ A 2,407,166 位置、1 ライブラリ呼出率 0.00107、
  無作為誤差帰無の 3/3 期待 0.004 件に対し観測 795 件。正準未帰属部位の 91.3% が 2/3 以上で再現 (帰属は 95.3%)。
  時点間: 未帰属 union 2,826、T1∩T2∩T3 964。
  ホモポリマー誤り機構は否定 (run>=4 が 0.49% 対 背景 0.45%)。文脈集中は強い (上位10 5-mer が部位の 40.6%、可測位置の 4.7%)。
  公開無修飾対照: ENA 照会 4 件 (MDA/RANDOM PCR/PCR/Streptomyces 全)、計 16,056 ラン取得。
  タイトルに WGA/unmethylated 等を明示するもの 51 ラン (PRJNA1426731 H. pylori 11 等)。Streptomyces の無修飾対照は 0 件。
  GPU 無し (list_compute = providers なし) のため同一モデル再 basecall は本セッションでは不可 → 実行可能性報告に留める。
  注記: 既存 unassigned_6mA_threshold_sweep_261015.tsv は正準規則と不一致 (T1 level>=0.5 で 4,248 対 正準 1,947)。
2026-11-08 C 完了。COG 20 カテゴリ x 3 時点 = 60 検定、GC501 五分位 x 深度五分位 25 層 MH + BH。
  BH 5% 通過は T3 の COG:S のみ (OR 1.256, q 0.026)。解釈可能カテゴリは全帰無 → P-C1 ほぼ支持。
  可動性因子 (transposase/integrase/recombinase 71 CDS): OR 2.221/2.111/2.189, BH q 0.044/0.066/0.026
  → P-C2 を外した。ただし和集合 22 部位・21 遺伝子、3 時点で 8 部位共通 (独立再現ではない)。
  プロファージ: geNomad が NC_003888.3 に virus 0 / provirus 0 → 検証不能 (帰無ではない)。
2026-11-08 報告書を Writing/open_questions_report2_261108.md に保存。数値監査 18/18 一致。レーン完了。
