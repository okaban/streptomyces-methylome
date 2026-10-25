# pileup_AAGCCCG_6mA

AAGCCCG の 6mA 位置 2,668（位置 1 と位置 2、両鎖）での 6mA 行（コード a）を、**正準パイルアップ**から per-library で抽出したもの。

    source  : /Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_canonical_261015/<lib>_pileup.bed
    columns : pos (0-based), strand, Nvalid_cov, Nmod, Nfail
    extracted_on : 2026-10-15

## 2026-10-15 の再抽出について

以前の抽出（2026-10-14）は `analysis/pileup/` を参照していた。そこは 3-2 と 3-4 が
追加フローセル前のままで、どちらも T3 のライブラリであるため **T3 の部位数だけが
過大**になっていた（TGGCCGGC で 1,877、正しくは 1,866）。T1 と T2 は影響を受けない。
来歴は `methyl/260102_M145/analysis/DATA_PROVENANCE.md` を参照。

判定規則（時点ごとに 3 ライブラリを合算、合算前の真の深度 Nvalid+Nfail >= 30、
合算 Nmod/Nvalid >= 0.5）を適用すると 位置1 323 / 350 / 335、位置2 93 / 100 / 118、和集合 416 / 450 / 453 部位。
