# PROGRESS — Q04 external expression validation (2026-11-07)

- 00:00 lane created; dirs scripts/ tables/ figures/
- read 127 canonical sites (4mC 2,909/2,869/2,840, 6mA 1,947/2,150/2,308) and S6d (|rho| floor 0.0863 on n=1,051 regulatory genes; manuscript text quotes 0.088)
- 01 built tables/gene_methylation_covariates_261107.tsv (7,863 chromosomal genes; 7,619 with an SCO id). Gene-body assignment keeps 93.1-93.3% of 4mC and 85.0-85.6% of 6mA canonical sites.
- GEO screened: 88 unique GSE from 4 term queries. Selected RNA-seq GSE132487 (M145, 3 cultures x 9 times), GSE174192 (M145 ISP2, 8 times x 3 reps), GSE158810 (M145 50 h control, 3 reps). Array candidate GSE58666 (M145 WT arm, 5 times x 3 reps, log2 cDNA/gDNA) pending probe->SCO mapping.
- Rejected so far: GSE53562 (7 samples, no replicates; VALUE is a log-ratio to the 18 h sample, not a level), GSE155796 (strain M1146, not M145), GSE30569/GSE31068 (submitter states no replicates; raw CEL only).
- CORRECTION (著者指摘, 2026-11-07): 上の PROGRESS 行に書いた正準部位数は wc -l の値（ヘッダ込み）で、4mC T2/T3 と 6mA 全時点が 1 ずつ多かった。正しくは 4mC 2,909/2,869/2,840、6mA 1,947/2,150/2,308。誤りは PROGRESS の散文のみ。scripts/01 は pandas read_csv でヘッダを除いて読んでおり、実行ログの sites 値は最初から正準値と一致。再実行して tables/gene_methylation_covariates_261107.tsv の md5 が不変（0e756153e702c4ef68b59e6cfbfcba27）であることを確認したので、被覆率 4mC 93.1/93.2/93.3%、6mA 85.0/85.6/85.4% は正しい分母で計算されている。
- 判定基準: 事前登録の予測は課題文で 0.088 と書かれていたが、S6d の表の値 0.0863 のほうが厳しいので **0.0863 で判定**する。
- FINDING (原稿の欠陥候補、本レーンでは原稿を触らない): 本文 R5 と S4b 段落は「|rho| = 0.088」と書くが、生成元の SuppTable19b_detectable_effect_sizes_261021.tsv は 0.0863。Fisher z で n=1,051, alpha=0.05, power=0.80 を解くと tanh((1.959964+0.841621)/sqrt(1048)) = 0.08632 なので表の値が正しく、本文の 0.088 は出所不明。
- 02/03/04 完了。外部 4 件（GSE132487/GSE174192/GSE158810 RNA-seq、GSE58666 アレイ）＋本研究自身を同一パイプラインで解析。
- 判定（閾値 0.0863）: 事前登録の主解析（levelkb × 全サンプル平均発現、4 共変量統制の偏 Spearman）で外部 48 検定すべてが閾値未満。最大 |rho| = 0.0468（GSE174192, 6mA T1, 調節遺伝子, p=0.152）。→ 陰性は再現。
- 検出力: 全染色体遺伝子 n=7,439-7,532 で床 0.0323-0.0325（本研究 0.0863 の約 1/2.7）。調節遺伝子 n=999-1,003 で床 0.0885-0.0887（本研究と実質同等、わずかに粗い）。「弱い検出力で陰性」ではない。
- 要著者判断: 調整前の生 Spearman では外部 48 検定中 5 件が 0.0863 超（GSE132487 4mC T1 +0.0900、GSE174192 6mA T1 -0.1163 ほか）。ただし本研究自身も生 +0.0845 とほぼ同値で、交絡は全染色体 4mC が遺伝子長、調節遺伝子 6mA が GC と特定済み。本文 R5 がどちらの統計量の主張かを明示する必要がある。
- メモ: /Users/okaban/obsidian/Research/rna-seq/Writing/q04_external_expression_261107.md
