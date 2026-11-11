# 退役した補足図の生成器

ここに置いてあるスクリプトは**走らせないこと**。いずれも出荷中の図版ファイル名に
書き込むため、実行すると現行の図を上書きする。2026-10-22 に
`Writing/check_provenance.py` の出力名衝突検査で見つかり、ここへ隔離した。

| ファイル | 書き込む先 | 隔離理由 |
|---|---|---|
| `suppfigure6_motif_exclusive_asymmetry.py` | `SuppFigure6.png` | 別内容の退役図（「4mC コールはモチーフ排他的／六塩基は片鎖のみ」）。現行の補足図 6 は `suppfigure6b_gc_control_261022.py` が作る GC 統制図で、凡例もそちらと一致する |
| `_recovered_suppfig17.py` | `SuppFigure17.png` | 復元作業の中間物（docstring 無し、素の import から始まる）。補足図 17 は現在どの凡例でも宣言されていない |

同じ事故の型: 補足図を 15 から 9 に改番したとき `suppfigure15_optionC_rebuilt.py` の
出力名が `SuppFigure15.png` のまま残り、のちに Figure 3a から移送した二本鎖模式図が
補足図 15 を占めたため、旧生成器の再実行で模式図を上書きする状態になっていた
（2026-10-22 に出力名を `SuppFigure9.png` へ訂正済み）。

ファイルの存在確認も凡例の grep もこの型を検出できない。どちらも「その名前の
ファイルがある」ことしか見ず、中身がどちらのものかを知らないため。

## 2026-10-24 の経緯（訂正を含む）

並行レーンが 2026-10-24 に 5 本をここへ移した。親エージェントが「出荷中の図の唯一の
生成器を退避したのではないか」と疑い 2 本を live に戻したが、**その判断は誤りだった**。
`figure_manifest.tsv` を見ると当該 2 本は別のスクリプトに割り当てられた名前を書くもので、
戻すと衝突する側だった。同日中に再退避した。

| ファイル | 書き込む先 | manifest 上のその名前の持ち主 | 現在の所在 |
|---|---|---|---|
| `suppfigure10_coverage_gradient.py` | `SuppFigure10.png` | `panels_261013/panel_s16_tss_sequence.py` | **退避（再退避済み）** |
| `suppfigure13_bgc_vs_density.py` | `SuppFigure13.png` | `recompute_261015/aagcccg261016/make_fig_261016.py` | **退避（再退避済み）** |
| `suppfigure15_optionA_4panel.py` | （出荷名に書かない） | — | 退避 |
| `suppfigure15_single_molecule.py` | （出荷名に書かない） | — | 退避 |
| `suppfigure18_epiallele_control.py` | （出荷名に書かない） | — | 退避 |

逆に、先のレーンが「別内容の退役図」として退避していた
`suppfigure6_motif_exclusive_asymmetry.py` は manifest 上 `SuppFigure3.png` の正当な
生成器だったので live に戻した。退避理由に書かれていた「`SuppFigure6.png` に書く」は
誤読で、実際に書くのは `SuppFigure3.png` である。

**教訓: 出力名から所有者を推測してはいけない。`figure_manifest.tsv` が割り当ての正。**
この誤りは `Writing/check_provenance.py` に 2026-10-24 に足した
`check_figure_name_coverage` が manifest と突き合わせて検出した。同検査は
「manifest の割り当てと実際の出力名が食い違う」を BLOCKER として出す。
