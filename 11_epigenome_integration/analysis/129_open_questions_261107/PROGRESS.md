# PROGRESS — open questions lane (261107)

- 21:03 着手。open_questions.tsv (Q01-Q08)、findings_ledger 未閉鎖 4 件 (V08/V17/V21 author_pending, V25 open) を確認。
- 21:06 著者指示で焦点を 6mA に変更。最優先=核様体仮説(i)の公開データ検証、次点=複製鎖 vs 遺伝子向き、5′塩基、未帰属 6mA の分布、遺伝子間過剰の内訳。
- 21:07 事前登録 PREREGISTRATION_261107.md を作成（P1 核様体/P2 複製鎖/P3 未帰属分布/P4 遺伝子間内訳）。
- 公開データ探索 (GEO/ENA, 2026-11-07):
  - GEO db=gds 4 クエリ: Hi-C 4 / ChIP-seq 83 / NAP 36 / S.coelicolor ChIP+HiC 19 件。
  - **NAP ChIP-seq は 124_occupancy_tracks_261104 で既に実施済み**（Gbn GSE165795、HupA PRJEB50438 自前処理、TF 11 本 GSE337352 ほか計 22 トラック）→ 全て陰性。重複しない。
  - **新規: ENA `tax_tree(1883) AND library_strategy="Hi-C"` 54 ラン。うち S. coelicolor A3(2) = PRJNA891940（4 ラン: L期 rep1/2 = SRR23618920/21、M期 rep1/2 = SRR23618922/23、NovaSeq PE、計 ~71 GB）。**同種・同座標系**なのでオルソログ対応づけ不要。
  - 同 BioProject に ssRNA-seq 6 ラン（SRR22030115-20）。
  - 他種 Hi-C: S. ambofaciens PRJEB108791 / S. venezuelae PRJEB41527, PRJEB102522 / S. lividans PRJNA1264485（今回は不使用）。
- 21:50 [結果] 遺伝子間 6mA 過剰の内訳（新規・事前登録 P4 に反する）: 非モチーフ 6mA の部位率を GC×深度 5x5 層別 MH OR で比較。divergent(プロモーター同士) OR 1.05/0.97/1.09（過剰なし）、convergent(ターミネーター同士) OR 2.57/2.83/2.89、tandem 2.04/1.89/1.85。3 時点で再現。表: tables/intergenic_subtype_6mA_261107.tsv
- 22:46 [結果] Hi-C (PRJNA891940, L期 SRR23618920) を自前で処理: 10.0M ペア、10 kb 851 ビン、ICE 平衡後の insulation。6mA(未帰属)密度 vs insulation の Spearman は GC・深度を残差化すると rho=-0.012/-0.006/-0.011 (p>=0.72)。境界ビン(16)±1 の MH OR 1.10/1.19/1.18 (MDE 1.33-1.37)。核様体仮説(i) は支持されず（OR 1.4 以上を除外できる検出力）。
- 22:55 M 期 Hi-C (SRR23618922) 3.17M ペアで同じ検定を実行、L 期と同じく帰無（境界 OR 1.02/1.12/1.10、MDE 1.38-1.44）。L 期と M 期の insulation は Spearman 0.734 で一致。
- 22:58 記録更新: Writing/open_questions.tsv（Q01-Q05/Q08 の状態更新、Q09-Q13 を新規追加、計 13 行）、Writing/open_questions_report_261107.md（総括）。
- 既検証と判定して再実行しなかったもの: NAP/TF 実測占有（124_、22 トラック全陰性）、複製鎖 vs 遺伝子向き（125_b1）、外部発現データ（128_）。
- 事前登録に対し P3・P4 を外した（未帰属 6mA は遺伝子間で一様ではない）。基準は動かしていない。
- 23:12 【自己訂正】記憶 mem_c9fc29f6e58c / mem_b74ff8636240 の警告（6mA の位置偏りは GC の窓幅に依存）に従い GC±25 bp で再検。収束 vs 発散は残る（2.40/2.61/2.16）が、「収束が同方向より強い」は消える（1.06/1.20/1.25）。頑健なのは「発散間隙にだけ過剰が無い」。報告と台帳を訂正済み。IH09 の gc25 再検（1.244/1.487/1.436）は小数第 2 位まで再現。mem_2a3031cf7762（2026-10-16、調整後 0.936）は再現せず＝未解決として記録。
- 23:14 成果物: PROGRESS.md / PREREGISTRATION_261107.md / tables 11 本 / Writing/open_questions.tsv（13 行）/ Writing/open_questions_report_261107.md。
