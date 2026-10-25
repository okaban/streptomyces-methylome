# coverage261015 — 被覆深度プロファイル（補足図 10）

    source       : /Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015/<lib>_pileup.bed
    extracted_on : 2026-10-15
    modkit       : v0.6.4

## coverage_50kb_by_library.csv

9 ライブラリ × 174 窓。各 (位置, 鎖) を 1 回だけ数えるため mod code を
`m`（C 位置）と `a`（A 位置）に限定して抽出した。この 2 コードで両鎖を合わせると
ゲノム上の全位置がちょうど 1 回ずつ現れる（+ 鎖が C なら m、G なら − 鎖の m、
A なら a、T なら − 鎖の a）。

    columns : library, window_index, start_bp, mean_true_depth, n_positions
    真の深度 = Nvalid_cov + Nfail

## oriter_windows_by_library.csv

oriC を**中心に置いた** 100/200/500 kb 窓と両終端窓の平均真深度。
補足図 10 の凡例が報告する ori/ter は**これではなく**、登録確定値
`REPL.ORITER.*` の出所である `110_topup_recompute/tables/repl_profile_canonical.csv`
（ゲノム先頭から等間隔に切ったグリッドビン）から取る。中心窓では約 5 % 低く出るため、
図・凡例・value_registry の三者一致を保つ目的で定義を揃えてある。本ファイルは
その差を検証するための参考値。
