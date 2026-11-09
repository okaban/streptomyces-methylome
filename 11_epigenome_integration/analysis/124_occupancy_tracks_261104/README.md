# 制御因子・核様体結合タンパク質の実測占有とメチル化の突き合わせ

解析 `11_epigenome_integration/analysis/124_occupancy_tracks_261104/`（2026-11-04）

棚卸し台帳の段 **C2-1 / C2-2 / C2-3**（制御因子の実測占有）と **I2-2 / I2-3**（NAP 占有）、
および **E1-3 / E1-4**（株パネル）の実現可能性判定。

> **台帳について**: 指定された `Writing/lit_methylation_function_approaches_261104.tsv` は
> ワークスペースにもアーティファクト store にも存在しなかった（`find /Users/okaban` 全走査と
> `host.artifacts(search=...)` の両方で 0 件）。各段の「実行手順」欄を読めていないため、
> 以下の設計はすべて本レーンで新規に立てたものである。

---

## 0. 要約

| 問い | 答え |
|---|---|
| S. coelicolor の実測占有データはあるか | **ある**。17 データセットを同定。うち 4 件が座標つきの処理済みピーク表を公開しており、**22 トラック**を構築できた（うち HupA 1 本は本レーンで生 FASTQ から自前作成） |
| ピーク内外でメチル化水準に差があるか | **ない**。水準の検定 29 件・密度の検定 54 件のすべてが FDR 5% で陰性（最小 BH: 水準 0.713、ピーク単位順位 0.425、密度 0.0505） |
| 予測（TSS）と実測占有で結論は変わるか | **変わる。ただし本研究に有利な向きに**。TSS±100 bp のメチル化欠損は再現するが、同じ欠損は実測占有では出ない。欠損は **TSS に付随し、タンパク質の占有には付随しない** |
| 株パネル（E1-3 / E1-4）は可能か | **不可能**。メチロームと発現が同一株・同一培養で対になった公開データは Streptomyces 属に 1 件も存在しない |
| 原稿に書けるか | **書ける**（C2-2 / C2-3 / I2-2）。permissive な主張を実測占有で補強する否定的証拠として |


---

## 1. C2-1 公開占有データの棚卸し

### 探したアクセッション空間

| アーカイブ | クエリ | 件数 |
|---|---|---|
| ENA portal | `tax_tree(1883) AND library_strategy="ChIP-Seq"` | 301 ラン / 43 スタディ（属全体）、うち S. coelicolor 15 スタディ |
| ENA portal | `tax_tree(1883) AND instrument_platform="PACBIO_SMRT"` | 895 ラン |
| ENA portal | `tax_tree(1883) AND instrument_platform="OXFORD_NANOPORE"` | 2,055 ラン |
| ENA portal | `tax_tree(1883) AND library_strategy="RNA-Seq"` | 3,142 ラン |
| GEO (db=gds) | `"Streptomyces coelicolor"[Organism] AND gse[ETYP]` ほか 4 クエリ | 121 シリーズ（S. coelicolor 104） |
| ArrayExpress / BioStudies | `Streptomyces ChIP` ほか | ChIP-chip 2 件（GSE23401 / GSE28908 相当）以外に新規なし |

候補として挙がっていた因子の探索結果は `tables/C2-1_regulator_search_coverage_261104.tsv`。
**S. coelicolor では存在しなかった**もの: **GlnR, PhoP, DasR, ArgR, Crp, AdpA, WhiA, Lsr2, HupS**
（GEO・ENA のどちらにも ChIP 系のランが 0 件）。
うち AdpA は S. griseus（PRJNA156443 / PRJNA156445）と S. venezuelae（PRJEB91679）のみ、
Lsr2 は S. venezuelae のみ（PRJNA475288、PRJNA1365143）、HupS / WhiA / BldD / BldC も
S. venezuelae のみで、いずれも NC_003888.3 座標には直接使えない。
DasR・ArgR・PhoP・Crp は該当 GEO シリーズがトランスクリプトーム（変異株の発現比較）であって
占有データではない。

### 見つかったデータセット（17 件）

| GEO | ENA_BioProject | protein | organism_strain | assay | processed_data_available | coordinate_system | status |
|---|---|---|---|---|---|---|---|
| GSE165795 | PRJNA697922 | Gbn (SCO1839), nucleoid-associated protein | S. coelicolor M145 | ChIP-seq (BGISEQ-500) | MACS2 peak table (per sample) | NC_003888.3 (native) | used in C2-2 |
| GSE337352 | PRJNA1483500 | 11 FLAG-tagged TFs: SCO1145, SCO1186, SCO1262, SCO2243, SCO3320, SCO3632, SCO3668, SCO3696, SCO3699, SCO4158, SCO5231 | S. coelicolor M145 | ChIP-seq (DNBSEQ-G400) | narrowPeak (per sample) | AL645882 — md5-identical to NC_003888.3 | used in C2-2 |
| GSE234438 | PRJNA981425 | RokL6 (SCO1447) | S. coelicolor M145 | ChIP-seq (NovaSeq 6000) | MACS2 peak table | NC_003888.3 (native) | used in C2-2 |
| GSE64903 | PRJNA272447 | ScbR, ScbR2 | S. coelicolor M145 | ChIP-seq (HiSeq 2000) | MACS1 peak table w/ gene annotation | NC_003888.3 (gi|32141095 header) | used in C2-2 |
| GSE59010 | PRJNA254102 | NdgR | S. coelicolor | ChIP-seq (GAIIx) | per-base GFF signal profile, no peak calls | NC_003888 (no version) | not used — signal only; peak calling would require our own threshold |
| GSE225370 | PRJNA935355 | CutR (and S. venezuelae CutR ortholog) | S. coelicolor M145 | ChIP-seq (NovaSeq 6000) | bedGraph coverage only | contig renamed 'sco'; no peak calls | not used — no peak calls, contig name not resolvable to a versioned accession |
| GSE136166/GSE136168 | PRJNA561422 | WblC (WhiB7) | S. coelicolor A3(2) | ChIP-seq (HiSeq 4000) | none (raw only) | — | not used — raw FASTQ only |
| GSE23401 | — | BldD | S. coelicolor M600 | ChIP-chip (two-colour array) | per-channel probe intensities | array probe IDs; needs platform mapping; M600 not M145 | not used — probe-level, different strain |
| GSE28908 | — | BldC | S. coelicolor M600 | ChIP-chip | per-channel probe intensities | array probe IDs | not used — probe-level, different strain |
| GSE33593 | — | SigR | S. coelicolor A3(2) | ChIP-chip (NimbleGen .pair) | raw .pair files | array probe IDs | not used — probe-level |
| GSE95760 | — | Zur | S. coelicolor | ChIP-chip (NimbleGen) | vsn-normalised probe ratios (GFF) | NC_003888 (no version), probe-level | not used — probe ratios, no peak calls |
| GSE311144 | PRJNA1368754 | ParB1 (SCP1 partitioning) | S. coelicolor A3(2) | ChIP-seq (NextSeq 500) | per-contig RPBPPM coverage csv | separate chr/SCP1/SCP2 csv | not used — coverage only; SCP1-focused and M145 lacks SCP1 |
| GSE263222 | PRJNA1096301 | ParT / ParABS | S. coelicolor A3(2) | ChIP-seq (NextSeq 500) | coverage csv | — | not used — coverage only |
| — (E-MTAB-11369) | PRJEB50438 | HupA, nucleoid-associated protein | S. coelicolor A3(2) HupA-FLAG | ChIP-seq (HiSeq 3000, PE) | none (raw FASTQ only, 19.5 GB) | — | identified; see I2 note on processing |
| — | PRJEB27299 | HrdB (principal sigma) | S. coelicolor A3(2) | ChIP-seq (HiSeq 4000) | none (raw only) | — | not used — raw FASTQ only |
| — | PRJEB31880 / PRJEB47709 | SigE | S. coelicolor M600 / A3(2) | ChIP-seq | none (raw only) | — | not used — raw FASTQ only; M600 for PRJEB31880 |
| — | PRJNA838444 | AtrA | S. coelicolor A3(2) | ChIP-seq | none (raw only) | — | not used — raw FASTQ only |


### 座標系の検証（重要）

GSE337352 の narrowPeak は contig 名が **AL645882**（INSDC）で NC_003888.3 ではない。
ENA から `AL645882.2` の FASTA を取得して本研究の `NC_003888.3.fna` と照合したところ、
**長さ 8,667,507 bp・md5 `fa44769ad49da012d251d8b1ec8672fc` で完全一致**した。
したがって座標変換は不要で、残る違いは BED 系の 0-based 半開区間と MACS の `peaks.xls`
（1-based 閉区間）の慣習だけである。本解析はすべて 0-based 半開に揃えた
（narrowPeak はそのまま、MACS1 / MACS2 のテキストは `start-1`）。

GSE225370（CutR）だけは contig 名が `sco` に改名されており、バージョンつきアクセッションに
解決できない。ピーク呼び出しも公開されていないため使用しなかった。

### 構築したトラック（22 本）

| track | protein | series | condition | n_peaks | median_width | pct_genome | resolution |
|---|---|---|---|---|---|---|---|
| Gbn_25h | Gbn (NAP) | GSE165795 | 25h culture | 2825 | 245 | 9.832 | MACS2 peak |
| Gbn_48h | Gbn (NAP) | GSE165795 | 48h culture | 2919 | 291 | 12.204 | MACS2 peak |
| TF_SCO3632 | SCO3632 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 57 | 288 | 0.198 | MACS narrowPeak |
| TF_SCO2243 | SCO2243 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 114 | 211 | 0.364 | MACS narrowPeak |
| TF_SCO3320 | SCO3320 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 55 | 328 | 0.272 | MACS narrowPeak |
| TF_SCO3699 | SCO3699 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 2 | 580 | 0.013 | MACS narrowPeak |
| TF_SCO4158 | SCO4158 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 13 | 363 | 0.056 | MACS narrowPeak |
| TF_SCO1145 | SCO1145 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 283 | 402 | 6.788 | MACS narrowPeak |
| TF_SCO1186 | SCO1186 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 349 | 207 | 1.062 | MACS narrowPeak |
| TF_SCO1262 | SCO1262 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 268 | 404 | 3.152 | MACS narrowPeak |
| TF_SCO3668 | SCO3668 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 10 | 254 | 0.028 | MACS narrowPeak |
| TF_SCO3696 | SCO3696 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 26 | 248 | 0.092 | MACS narrowPeak |
| TF_SCO5231 | SCO5231 (FLAG ChIP-seq) | GSE337352 | M145 liquid culture | 113 | 248 | 0.466 | MACS narrowPeak |
| RokL6_M24 | RokL6 (SCO1447) | GSE234438 | M24 | 278 | 812 | 3.033 | MACS2 peak |
| RokL6_M48 | RokL6 (SCO1447) | GSE234438 | M48 | 1 | 420 | 0.005 | MACS2 peak |
| RokL6_G24 | RokL6 (SCO1447) | GSE234438 | G24 | 1 | 436 | 0.005 | MACS2 peak |
| RokL6_G48 | RokL6 (SCO1447) | GSE234438 | G48 | 1 | 455 | 0.005 | MACS2 peak |
| ScbR_30h | ScbR | GSE64903 | 30h | 144 | 2430 | 4.332 | MACS1 peak |
| ScbR2_42h | ScbR2 | GSE64903 | 42h | 491 | 1471 | 9.214 | MACS1 peak |
| TF_union_11 | union of 11 FLAG-ChIP TFs | GSE337352 | M145 liquid culture | 680 | 235 | 2.756 | MACS narrowPeak (merged union) |
| Gbn_union | Gbn (NAP), 25h+48h union | GSE165795 | union | 2974 | 261 | 10.685 | MACS2 peak (merged union) |


`TF_union_11` は GSE337352 の 11 TF ピークの和集合、`Gbn_union` は Gbn 25h/48h の和集合
（いずれも本解析で導出）。ピーク数 10 未満のトラック（TF_SCO3699 n=2、TF_SCO3668 n=10 は
境界、RokL6 の M48 / G24 / G48 は各 n=1）は検定から外した。

---

## 2. C2-2 ピーク内外のメチル化水準

### 設計

- **統計単位はピーク**。部位ではない（部位単位だと 1 ピークに複数部位が入ったときに擬似反復になる）。
- **測定値**: 本研究の正本パイルアップ（`pileup_nofilt_261015/*_nofilt.bed`、9 ライブラリ）から
  修飾可能位置のリードを 9 ライブラリ分プールし、部位水準 = ΣNmod / ΣNvalid。
  位置はゲノム配列から再抽出して本文の正準値と一致することを確認した
  （TGGCCGGC C5 = **2,093**、AAGCCCG C5 / A1 = **1,334**）。
- **ヌル**: 各ピークと**幅が完全に一致**し、**領域 GC が一致**するランダム領域。
  GC 許容幅は 0.02 から始めて、必要なら 0.035 / 0.05 / 0.08 / 0.20 と段階的に緩める
  （実際に達成した |ΔGC| は表の `mean_abs_dGC_*` 列）。1 比較あたり 400 反復。
  水準の検定では対照も「同じ配列文脈の測定可能な部位を 1 つ以上含む」ことを要求し、
  密度の検定ではこの条件を課さない（課すと密度がゼロの側が消える）。
- **交絡**: 遺伝子長に相当するのはピーク幅で、これは完全一致で統制した。GC は真の交絡なので
  上記のとおり一致させた。119 レーンで素の検定の有意が遺伝子長の層別で消えた前例に合わせ、
  素の Fisher 的な比較は一切報告していない。
- **検定**: (a) ピーク平均の平均 vs ヌル分布の経験 p、(b) 各ピークを自分自身の対照 400 個の中で
  順位づけた u 統計量（帰無で一様）に対する 1 標本 Wilcoxon。(b) は分布形の影響を受けず
  (a) より検出力が高い。密度・水準・u のそれぞれで BH 補正。

### GC は実際に交絡していた

実測ピークの平均 GC はゲノム平均（0.721）から大きくずれる。とくに `TF_SCO1145` は **0.632**、
`RokL6_M24` は 0.622、`TF_SCO1262` は 0.622 で、GC を合わせない素の比較では
TGGCCGGC（GC 75%）の密度が 0.17-0.50 倍に「枯渇」して見える。GC 一致ヌルを使うと
この見かけの枯渇は消える（`density_ratio` 列）。

### 結果 — すべて陰性

- **密度**: 54 比較。BH 最小 **0.0505**（`TF_SCO1186` の TGGCCGGC、生 p = 0.0025）。閾値の直上だが
  FDR 5% を通るものは **0 件**。名目で p < 0.05 の 10 比較は向きが揃っておらず（0.20-5.28 倍）、
  一貫した信号を成していない。
- **水準**: 29 比較。BH 最小 **0.713**（平均ベース）、**0.425**（ピーク単位 u）。**0 件**。
  効果量は全 29 比較で **0.931 - 1.139 倍**の範囲に収まる。

水準の検定（HupA を除く 26 件、ピーク単位 u の p 昇順。HupA は 4 節の表）:

| track | motif | n_peaks_with_site | mean_peak_GC | obs_mean_level | null_mean_level | level_ratio | p_level | p_level_BH | mean_u | p_perpeak_u | p_perpeak_u_BH | mde_level_ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TF_SCO1145 | 6mA_AAGCCCG | 105 | 0.6322 | 0.37764 | 0.36103 | 1.046 | 0.13466 | 0.68828 | 0.5507 | 0.024512 | 0.38137666666666664 | 0.9436-1.0543 |
| RokL6_M24 | 4mC_AAGCCCG | 51 | 0.6218 | 0.58281 | 0.56175 | 1.0375 | 0.24439 | 0.68828 | 0.5573 | 0.042896 | 0.38137666666666664 | 0.9392-1.0618 |
| TF_SCO1262 | 4mC_AAGCCCG | 38 | 0.6215 | 0.5784 | 0.55275 | 1.0464 | 0.22444 | 0.68828 | 0.572 | 0.044005 | 0.38137666666666664 | 0.9273-1.0795 |
| TF_union_11 | 6mA_AAGCCCG | 41 | 0.7103 | 0.4259 | 0.38875 | 1.0956 | 0.14713 | 0.68828 | 0.5891 | 0.073822 | 0.47984299999999996 | 0.8771-1.1230 |
| Gbn_48h | 6mA_AAGCCCG | 170 | 0.6823 | 0.36111 | 0.37907 | 0.9526 | 0.13466 | 0.68828 | 0.4638 | 0.095181 | 0.49494119999999997 | 0.9438-1.0600 |
| Gbn_25h | 6mA_AAGCCCG | 144 | 0.6828 | 0.36129 | 0.37867 | 0.9541 | 0.16708 | 0.68828 | 0.4659 | 0.129097 | 0.5594203333333333 | 0.9371-1.0677 |
| Gbn_union | 6mA_AAGCCCG | 152 | 0.6855 | 0.36658 | 0.37839 | 0.9688 | 0.29925 | 0.68828 | 0.4717 | 0.241304 | 0.7752246666666667 | 0.9485-1.0743 |
| TF_SCO1145 | 4mC_AAGCCCG | 105 | 0.6322 | 0.58231 | 0.5692 | 1.023 | 0.24439 | 0.68828 | 0.5238 | 0.262502 | 0.7752246666666667 | 0.9607-1.0355 |
| Gbn_25h | 4mC_TGGCCGGC | 159 | 0.6828 | 0.63359 | 0.6455 | 0.9816 | 0.26434 | 0.68828 | 0.4776 | 0.375807 | 0.7752246666666667 | 0.9698-1.0329 |
| Gbn_union | 4mC_TGGCCGGC | 175 | 0.6855 | 0.63617 | 0.64503 | 0.9863 | 0.32918 | 0.68828 | 0.4819 | 0.39863 | 0.7752246666666667 | 0.9735-1.0283 |
| ScbR2_42h | 4mC_AAGCCCG | 112 | 0.715 | 0.54998 | 0.54395 | 1.0111 | 0.64339 | 0.76037 | 0.5189 | 0.425502 | 0.7752246666666667 | 0.9574-1.0431 |
| TF_SCO1186 | 6mA_AAGCCCG | 14 | 0.762 | 0.46887 | 0.41182 | 1.1386 | 0.13466 | 0.68828 | 0.632 | 0.42627 | 0.7752246666666667 | 0.8370-1.1830 |
| ScbR_30h | 4mC_TGGCCGGC | 56 | 0.7194 | 0.65044 | 0.63911 | 1.0177 | 0.45387 | 0.7208423529411765 | 0.5312 | 0.48298 | 0.7752246666666667 | 0.9509-1.0482 |
| TF_SCO1186 | 4mC_AAGCCCG | 13 | 0.762 | 0.49856 | 0.53548 | 0.9311 | 0.34414 | 0.68828 | 0.451 | 0.497314 | 0.7752246666666667 | 0.8406-1.1208 |
| ScbR2_42h | 6mA_AAGCCCG | 112 | 0.715 | 0.39829 | 0.38906 | 1.0237 | 0.5187 | 0.722943 | 0.5189 | 0.503383 | 0.7752246666666667 | 0.9369-1.0679 |
| ScbR_30h | 6mA_AAGCCCG | 50 | 0.7194 | 0.37255 | 0.39134 | 0.952 | 0.34165 | 0.68828 | 0.4696 | 0.505357 | 0.7752246666666667 | 0.9008-1.0958 |
| Gbn_48h | 4mC_TGGCCGGC | 195 | 0.6823 | 0.63897 | 0.64495 | 0.9907 | 0.47132 | 0.7208423529411765 | 0.487 | 0.513981 | 0.7752246666666667 | 0.9715-1.0258 |
| ScbR_30h | 4mC_AAGCCCG | 50 | 0.7194 | 0.53021 | 0.54376 | 0.9751 | 0.43641 | 0.7208423529411765 | 0.4709 | 0.536694 | 0.7752246666666667 | 0.9337-1.0763 |
| ScbR2_42h | 4mC_TGGCCGGC | 153 | 0.715 | 0.63277 | 0.63982 | 0.989 | 0.52868 | 0.722943 | 0.4882 | 0.598598 | 0.8051277 | 0.9707-1.0323 |
| RokL6_M24 | 6mA_AAGCCCG | 51 | 0.6218 | 0.3807 | 0.36499 | 1.0431 | 0.41147 | 0.7208423529411765 | 0.4844 | 0.619329 | 0.8051277 | 0.9099-1.1070 |
| Gbn_25h | 4mC_AAGCCCG | 145 | 0.6828 | 0.5539 | 0.55035 | 1.0064 | 0.75312 | 0.8513530434782609 | 0.5112 | 0.656256 | 0.8125074285714285 | 0.9593-1.0418 |
| TF_union_11 | 4mC_AAGCCCG | 40 | 0.7103 | 0.55097 | 0.54746 | 1.0064 | 0.88778 | 0.93516 | 0.509 | 0.787868 | 0.8979292173913044 | 0.9223-1.0716 |
| TF_union_11 | 4mC_TGGCCGGC | 21 | 0.7103 | 0.65431 | 0.63987 | 1.0226 | 0.55611 | 0.722943 | 0.5248 | 0.794322 | 0.8979292173913044 | 0.9077-1.0824 |
| Gbn_48h | 4mC_AAGCCCG | 166 | 0.6823 | 0.55091 | 0.54996 | 1.0017 | 0.93516 | 0.93516 | 0.5033 | 0.858589 | 0.90454936 | 0.9613-1.0355 |
| TF_SCO1262 | 6mA_AAGCCCG | 38 | 0.6215 | 0.35993 | 0.37144 | 0.969 | 0.60349 | 0.7471780952380952 | 0.4953 | 0.869759 | 0.90454936 | 0.8809-1.1164 |
| Gbn_union | 4mC_AAGCCCG | 150 | 0.6855 | 0.54896 | 0.55005 | 0.998 | 0.91771 | 0.93516 | 0.4989 | 0.966327 | 0.966327 | 0.9590-1.0406 |


### 「0 件」の検出力（必須併記）

「差がなかった」は検出可能だった最小効果量とセットでのみ意味を持つ。`mde_level_ratio` 列は
ヌル分布の 2.5-97.5 パーセンタイルをヌル中央値で割ったもので、両側 5% で棄却できた比の範囲を表す。
最も検出力の高いトラックでは:

| トラック | モチーフ | 部位を含むピーク数 | 棄却できた比の範囲 | 観測比 |
|---|---|---|---|---|
| Gbn_48h | 4mC TGGCCGGC | 195 | < 0.972 または > 1.026 | 0.991 |
| Gbn_union | 4mC TGGCCGGC | 175 | < 0.974 または > 1.028 | 0.986 |
| Gbn_48h | 6mA AAGCCCG | 170 | < 0.944 または > 1.060 | 0.953 |
| Gbn_48h | 4mC AAGCCCG | 166 | < 0.961 または > 1.036 | 1.002 |

すなわち Gbn 占有領域については **4mC で約 3%、6mA で約 6% の水準差があれば検出できた**。
TF トラック（ピーク数が少ない）では 8-20% が限界で、そこは「差がない」ではなく
「中程度の差は検出できない」と書くべき範囲である。

---

## 3. C2-3 予測サイトと実測サイトで結論が変わるか

本研究はこれまで Jeong 2016 の実験 TSS（3,570）に基づく**予測**制御領域を使ってきた。
同じパイプライン（同じ GC・幅一致ヌル、同じ統計単位）を予測領域に通して比較した。

### 予測領域での結果（`tables/C2-3_predicted_vs_measured_261104.tsv`）

| track | motif | n_regions | pct_genome | n_regions_with_site | level_ratio | p_level | mean_u | p_perpeak_u | mde_level_ratio |
|---|---|---|---|---|---|---|---|---|---|
| PRED_promoter_TSS300 | 4mC_TGGCCGGC | 2736 | 10.872 | 168 | 1.0149 | 0.34913 | 0.53 | 0.1808 | 0.9686-1.0330 |
| PRED_promoter_TSS300 | 4mC_AAGCCCG | 2736 | 10.872 | 151 | 0.9655 | 0.10973 | 0.4642 | 0.128305 | 0.9560-1.0396 |
| PRED_promoter_TSS300 | 6mA_AAGCCCG | 2736 | 10.872 | 153 | 0.9882 | 0.70324 | 0.4919 | 0.706811 | 0.9419-1.0608 |
| PRED_TSSproximal_100bp | 4mC_TGGCCGGC | 2840 | 7.307 | 82 | 1.0057 | 0.77805 | 0.4887 | 0.627386 | 0.9613-1.0405 |
| PRED_TSSproximal_100bp | 4mC_AAGCCCG | 2840 | 7.307 | 113 | 0.9185 | 0.00249 | 0.4201 | 0.00497 | 0.9515-1.0455 |
| PRED_TSSproximal_100bp | 6mA_AAGCCCG | 2840 | 7.307 | 114 | 0.8961 | 0.00748 | 0.4217 | 0.005839 | 0.9333-1.0697 |


- **TSS±100 bp**: 4mC AAGCCCG が **0.919 倍（p = 0.0025）**、6mA AAGCCCG が **0.896 倍（p = 0.0075）**、
  ピーク単位 u でも整合（p = 0.005 / 0.006）。**本文の TSS 近傍欠損が、この独立なヌルでも再現した**。
- **-300..TSS プロモーター**（より広い定義）: 4mC AAGCCCG 0.966 倍（p = 0.11）。欠損は検出されない。
  欠損は TSS の直近に限局しており、慣用的なプロモーター窓に広げると薄まる。
- **TGGCCGGC は TSS±100 bp でも水準差なし**（1.006 倍）。密度だけが 0.735 倍（p = 0.0075）で、
  これは水準の話ではなく配列の配置の話である。

### 実測占有での結果

同じ指標を実測占有の 11 トラック（29 比較）で計算すると、前節のとおり**すべて陰性**で
効果量は 0.93-1.14 倍。TSS で見えた 0.90-0.92 倍の欠損に相当するものは実測占有には無い。

### 両者を直接分ける検定

TSS±100 bp の窓 2,840 個を、**実測占有ピーク（HupA を除く 19 トラックの和集合、ゲノムの 23.6%）と
重なるか否か**で二分した（重なる 1,146、重ならない 1,694）:

| subset | motif | n_regions | n_with_site | obs_level | null_level | level_ratio | p_level | mean_u | p_u | mde |
|---|---|---|---|---|---|---|---|---|---|---|
| TSS+-100 OCCUPIED | 4mC_TGGCCGGC | 1146 | 36 | 0.67626 | 0.64776 | 1.044 | 0.13965 | 0.5512 | 0.220408 | 0.9387-1.0515 |
| TSS+-100 OCCUPIED | 4mC_AAGCCCG | 1146 | 50 | 0.52177 | 0.54958 | 0.9494 | 0.19451 | 0.4723 | 0.496389 | 0.9284-1.0669 |
| TSS+-100 OCCUPIED | 6mA_AAGCCCG | 1146 | 50 | 0.34799 | 0.3741 | 0.9302 | 0.18953 | 0.4468 | 0.218398 | 0.8944-1.1114 |
| TSS+-100 UNOCCUPIED | 4mC_TGGCCGGC | 1694 | 46 | 0.62965 | 0.64519 | 0.9759 | 0.35411 | 0.443 | 0.130232 | 0.9484-1.0549 |
| TSS+-100 UNOCCUPIED | 4mC_AAGCCCG | 1694 | 63 | 0.49081 | 0.54891 | 0.8942 | 0.00249 | 0.3837 | 0.001676 | 0.9362-1.0594 |
| TSS+-100 UNOCCUPIED | 6mA_AAGCCCG | 1694 | 64 | 0.32943 | 0.37748 | 0.8727 | 0.00748 | 0.4 | 0.010033 | 0.9097-1.1000 |


**占有されていない TSS 窓でこそ欠損が出る**（4mC AAGCCCG 0.894 倍 p = 0.0025、6mA 0.873 倍 p = 0.0075）。
占有されている TSS 窓では点推定は同じ向き（0.949 / 0.930 倍）だが有意ではなく、検出可能範囲
（0.928-1.067、0.894-1.111）もほぼ同等なので、これは検出力不足だけでは説明しきれない。

ただし **2 つの部分集合の差そのもの**を GC 層別の順列検定（2,000 回、GC 五分位内で占有ラベルを交換）で
直接検定すると:

| motif | n_occ | n_unocc | mean_level_occupied | mean_level_unoccupied | diff | p_perm_GCstratified | mde_abs_diff |
|---|---|---|---|---|---|---|---|
| 4mC_TGGCCGGC | 36 | 46 | 0.67626 | 0.62965 | 0.04661 | 0.03898 | -0.0408 to 0.0436 |
| 4mC_AAGCCCG | 50 | 63 | 0.52177 | 0.49081 | 0.03096 | 0.54923 | -0.0387 to 0.0653 |
| 6mA_AAGCCCG | 50 | 64 | 0.34799 | 0.32943 | 0.01856 | 0.37031 | -0.0585 to 0.0470 |


3 検定中、名目で有意なのは 4mC TGGCCGGC のみ（占有側が **+0.047 高い**、p = 0.039、BH では 0.117）。
AAGCCCG の 4mC・6mA は差なし。つまり**占有の有無が TSS 近傍の欠損を説明する証拠はない**。
向きは占有側がむしろ「高い」側で、これは「タンパク質が DNA を覆ってメチル化を妨げる」という
予測と逆である（補足図 11b の FIMO 予測サイトで 1.157 倍の増加が見えたのと同じ向き）。

### C2-3 の結論

予測（TSS）から実測（ChIP）に替えると**結論は変わる。向きは本研究に有利**:

1. 本研究が報告している TSS 近傍のメチル化欠損は、**転写開始点に付随する性質**であって
   制御因子・NAP の占有に付随する性質ではない。占有の有無で層別しても欠損は占有されない側に残る。
2. したがって「メチル化は制御ゲノムを空間的に組織するが、転写を指示はしない」という permissive な
   主張は、実測占有を入れても崩れない。**むしろ「タンパク質が結合しているからメチル化が低い」という
   対抗説明を 1 つ潰せる**。
3. 逆に「予測サイトでしか見えていない現象ではないか」という査読コメントに対しては、
   「予測 TSS では見えるが実測占有では見えない。ただし両者は別の対象を測っており、
   欠損は TSS 側に帰属する」と答えられる。

---

## 4. I2-2 / I2-3 核様体結合タンパク質（NAP）

S. coelicolor で**実測占有が公開されている NAP は Gbn（SCO1839、GATC 結合 NAP、GSE165795）のみ**。
本解析で最も検出力の高いトラックがこの Gbn で、25h・48h・和集合の 9 比較すべてが陰性
（効果量 0.953-1.006 倍、4mC では約 3% の差まで検出できる条件）。

- **Lsr2**: S. coelicolor には存在しない。S. venezuelae のみ（PRJNA475288 = Lsr2 ChIP-seq、
  PRJNA1365143 = LsrL/Lsr2、GSE115538 ほか）。Gehrke et al. 2019 eLife 10.7554/eLife.47691 が
  潜在クラスターの抑制を報告しているのも S. venezuelae である。NC_003888.3 座標には直接使えず、
  オーソログ経由の座標移植は本レーンでは行わなかった。
- **HupA**: S. coelicolor A3(2) の HupA-FLAG ChIP-seq が **PRJEB50438 / E-MTAB-11369** に存在する
  （8 ラン PAIRED、正常 TopA `thio0` 3 反復 + TopA 枯渇 `thio2` 3 反復 + untagged 対照 PS04 2 反復）。
  公開されているのは生 FASTQ のみで、処理済みピークは無い。
- **HupS, SMC, ParB**: S. venezuelae のみ（PRJEB41570 / PRJEB41571 / PRJEB81515）。
  S. coelicolor の ParB1 は SCP1 分配に関する GSE311144 で、M145 は SCP1 を保持しないため使えない。


---

## 5. E1-3 / E1-4 株パネルの実現可能性

### 判定: **不可能**

「他の Streptomyces 株でメチロームと発現が対になった公開データ」は、属全体を走査して **0 件**だった。

### 探した空間

ENA portal API の `tax_tree(1883)`（Streptomyces 属）に対し、プラットフォーム別・戦略別に全ランを列挙した:

| 条件 | ラン数 |
|---|---|
| `instrument_platform="PACBIO_SMRT"` | 895 |
| 同上のうち速度論を復元できる submitted_format（`PACBIO_HDF5` または `BAM`） | **123** |
| `instrument_platform="OXFORD_NANOPORE"` | 2,055 |
| 同上のうち生シグナルを含む（`OXFORDNANOPORE_NATIVE`） | **50** |
| `library_strategy="RNA-Seq"` | 3,142 |
| 長鎖（PacBio または ONT）と RNA-seq の**両方**がある種 | 45 |

GEO 側では `"Streptomyces"[Organism]` × メチローム語（methylome / DNA methylation / SMRT /
PacBio / nanopore / bisulfite / 6mA / 4mC）で **3 ヒット**のみ、うち実体は
`Streptomyces sp. WMMB 714` の PacBio RS II サンプル 1 件（GSM1711764 / GPL20529）だけだった。
Blow et al. の原核生物メチローム網羅調査 **GSE69872 には Streptomyces サンプルが 0 件**である。

### 候補ごとの判定（`tables/E1-3_E1-4_strain_panel_feasibility_261104.tsv`）

| candidate_strain_species | methylome_source | methylome_data_type | expression_data_available | verdict_reason |
|---|---|---|---|---|
| S. coelicolor A3(2) (NCTC strain) | PRJEB6403 (NCTC_3000) | PacBio RS, 3 runs, submitted_format PACBIO_HDF5 (bax.h5 -> IPD recoverable) | ENA: 504 RNA-seq runs for S. coelicolor A3(2) (other labs/conditions) | Same species as this study, so adds no between-strain contrast; strain field empty (NCTC isolate, not M145); genome-survey sample with no matched culture. bax.h5 requires legacy SMRT Analysis 2.x ipdSummary, which is unsupported on current platforms. |
| S. albus subsp. albus DSM 41398 | PRJEB11297 | PacBio RS, 5 runs, PACBIO_HDF5 | ENA: 25 RNA-seq runs for S. albus, predominantly strain J1074 | Methylome strain (DSM 41398) and transcriptome strain (J1074) differ; no paired culture; legacy bax.h5 pipeline. |
| S. albidoflavus | PRJEB108104 (Revio Swiss Streptomyces genomes) | PacBio Revio, 11 runs, BAM;BAI | ENA: 245 RNA-seq runs | Genome-assembly project; isolate identity not matched to any transcriptome strain; whether the HiFi BAMs retain per-base kinetics (fi/ri/fp/rp) or carry 6mA/4mC MM/ML calls is not stated in the metadata and would have to be checked file by file. |
| S. griseus (NCTC7807, NCTC13033) | PRJEB6403 | PacBio RS, 3 runs, PACBIO_HDF5 | ENA: 6 RNA-seq runs for S. griseus | Different strains between methylome and transcriptome; 6 RNA-seq runs only; legacy bax.h5. |
| S. tendae / rochei / microflavus / alboniger / tanashiensis / avidinii and 20 further species | PRJEB71438 | Oxford Nanopore, 50 runs with submitted_format OXFORDNANOPORE_NATIVE (raw signal present) | ENA: 3-12 RNA-seq runs per species, different strains | Raw-signal ONT exists so a modification caller could be run, but it is a one-isolate-per-species genome survey; no transcriptome from the same strain or culture. |
| Streptomyces sp. WMMB 714 | GSM1711764 / GPL20529 | PacBio RS II methylome sample in GEO | none | Single isolate, no transcriptome of any kind. |
| S. venezuelae ATCC 10712 | — | no PacBio and no raw-signal ONT run found in ENA for this species | abundant RNA-seq and ChIP-seq | No methylome-capable sequencing exists for the species that otherwise has the richest Streptomyces functional-genomics resource. |


### 不可能である理由（3 点）

1. **対になっていない**。メチローム側のデータはすべてゲノム決定プロジェクトの産物で、
   1 分離株・1 条件・RNA は採られていない。発現側は別の研究室の別株・別条件である。
   本研究は T1/T2/T3 × 3 培養という**同一培養からの対**で成り立っており、この設計を
   公開データで再現できる株は属内に 1 つも無い。
2. **株が違う**。S. albus では メチロームが DSM 41398、トランスクリプトームが J1074。
   S. griseus では NCTC7807/NCTC13033 と別株。株間でメチル化系の有無が変わる以上、
   これは比較にならない。
3. **再処理の経路が切れている**。速度論を持つ PacBio はほぼ `bax.h5`（PacBio RS 世代）で、
   `ipdSummary` を回すには legacy SMRT Analysis 2.x が要る。これは現行環境では動かない。
   新しい Revio の `BAM;BAI`（PRJEB108104、S. albidoflavus 11 ラン）は、速度論タグ
   （fi/ri/fp/rp）や 6mA/4mC の MM/ML 呼び出しを残しているかがメタデータから判断できず、
   ファイル単位で確認しないとわからない。

### それでもやれる代替設計（実行はしていない）

株パネルそのものは諦めるとして、主張の一般性を調べたいなら次の 2 つが現実的:

- **(a) モチーフ配置の比較ゲノミクス**。REBASE / PacBio 由来の既報メチル化モチーフを持つ
  Streptomyces ゲノムを集め、各ゲノムで「そのモチーフが TSS 近傍にどう配置されているか」だけを
  比較する。部位ごとの水準も発現も要らないので、公開ゲノム配列と既報モチーフだけで完結する。
  ただし「水準」ではなく「配置」の主張しか作れない。
- **(b) S. albidoflavus の Revio BAM を 1 本だけ検査**（PRJEB108104、ERR16674745 ほか）。
  `samtools view` でヘッダと最初の数リードを見れば、速度論タグまたは MM/ML タグの有無が
  即座に決まる。もし MM/ML が入っていれば、同種の RNA-seq 245 ランとの**条件非対応**の
  比較は可能になる（同一培養の対ではないので、主張は「配置の保存」までに限られる）。

いずれも本研究の主要結論を支える力は弱く、**E1-3 / E1-4 は needs_data として閉じるのが妥当**と判断する。


### HupA を自前で処理した結果（暫定トラック）

PRJEB50438 の生 FASTQ を一部処理して HupA 占有トラックを 1 本作った。

- 取得: ENA FTP から ChIP（`ERR8135850` = HupAFLAG_thio0_1）と untagged 対照（`ERR8135856` = PS04_1）の
  **各 mate 先頭 400 万リード対**のみ。全 8 ラン 19.5 GB をそのまま落とすと、実測転送速度
  （1 ファイル 386 MB に 12.7 分）から 10 時間超になるため、この縮約を選んだ。
- 処理: bowtie2 2 → samtools sort → MACS2 `callpeak -f BAMPE -g 8.67e6 -q 0.01 --keep-dup auto`。
  整列率 68.5%（ChIP）/ 60.1%（対照）。
- 結果: **1,955 ピーク、計 1.10 Mb（ゲノムの 12.71%）、幅中央値 412 bp・平均 563 bp、fold enrichment 1.4-1.6**
  （NAP らしく広く浅い）。

| モチーフ | 部位を含むピーク | 密度比 | p | 水準比 | p | 検出可能だった比 |
|---|---|---|---|---|---|---|
| 4mC TGGCCGGC | 166 | 0.895 | 0.212 | **1.007** | 0.648 | <0.973 / >1.025 |
| 4mC AAGCCCG | 235 | 1.250 | 0.005 | **1.002** | 0.880 | <0.969 / >1.032 |
| 6mA AAGCCCG | 234 | 1.230 | 0.007 | **0.966** | 0.212 | <0.950 / >1.049 |

**HupA 占有領域でもメチル化水準に差はない**（約 3-5% の差なら検出できた条件で）。
AAGCCCG の密度が 1.23-1.25 倍という名目上の濃縮はあるが、BH 補正後は 0.0505 で閾値を通らない。

このトラックは **暫定（PROVISIONAL）**である: 反復 1 本のみ、ファイル先頭からの縮約サンプル
（先頭リードはフローセル上の位置が偏るので厳密な無作為抽出ではない）、ピーク呼び出しは本レーンの自前。
原稿に載せるなら全 8 ランを落としてから再計算するべきで、現状は「水準差が無いという結論を
HupA でも覆せなかった」以上のことは言えない。


---

## 6. 方法の要点と限界

### 使ったデータ
- メチル化: `/Users/okaban/bioinfo/methyl/260102_M145/analysis/pileup_nofilt_261015/*_nofilt.bed`
  （9 ライブラリ、全リード、モチーフ位置に限定済みで各 6.4 MB）。
- 参照: `11_epigenome_integration/data/NC_003888.3.fna`（8,667,507 bp）。
- TSS: `11_epigenome_integration/analysis/18_tss_analyses/jeong2016_all_tss.csv`（3,570）。
- 占有: 上記 4 シリーズの処理済みピーク表（GEO FTP から直接取得、`raw/tracks/`）。

### 限界
1. **AAGCCCG の 4mC と 6mA は同じ 1,334 位置**なので、密度の検定ではこの 2 行が構造的に重複する。
   BH はその分だけ保守側に振れている（独立な部位集合は TGGCCGGC と AAGCCCG の 2 つ）。
2. **実測占有は条件が本研究と一致しない**。Gbn は 25 h / 48 h の液体培養、GSE337352 は M145 の
   液体培養、ScbR/ScbR2 は 30 h / 42 h。本研究の T1/T2/T3 とは培地も時点も違う。
   メチル化は時点間でほぼ不変（占有率の時点差 ±1.5% 以内）なので致命的ではないが、
   占有側が本研究の培養を反映していないことは明記が要る。
3. **TF トラックは検出力が低い**。GSE337352 の 11 TF のうち、水準の検定に乗ったのは 4 本だけ
   （SCO1145 / SCO1186 / SCO1262 と和集合）。ピークが狭く（中央値 207-404 bp）、
   修飾可能位置の密度が低い（TGGCCGGC で 1/4.1 kb）ため、狭いピークにはそもそも部位が入らない。
4. **GC 一致の達成度**。GC が極端なトラック（TF_SCO1145 の平均 GC 0.632 に対しゲノム 0.721）では
   許容幅を 0.02 から段階的に緩めている。達成した平均 |ΔGC| は表に記録した。
   ゲノム全体の 400 bp 窓 GC は 1 パーセンタイルで 0.610 なので、GC 0.57 前後のピークには
   厳密な一致対照が存在しない。この種のトラックの結果は「GC をおおむね揃えた」程度に読むべき。
5. **ChIP-chip は使えなかった**。BldD（GSE23401）と BldC（GSE28908）は S. coelicolor で
   占有が公開されている数少ない発生制御因子だが、公開されているのはプローブ強度で、
   かつ株が M600（M145 ではない）である。den Hengst 2010 の BldD レギュロンを使いたい場合は
   論文の補足表（転写単位リスト）経由になり、ピーク座標ではないので本解析の枠組みには乗らない。
6. **ピーク呼び出しは各原著のパラメータのまま**使った。再呼び出しはしていないので、
   トラック間で閾値（Gbn は q < 1e-20、GSE337352 は narrowPeak 既定）が揃っていない。
   トラックごとに独立に対照と比べる設計なので致命的ではないが、トラック間の効果量の直接比較はできない。
7. **部位単位の検定は一切していない**。ピーク内の部位をまとめて t 検定すると擬似反復になるため、
   すべてピーク単位（またはピークを対照分布の中で順位づけた u 統計量）で報告した。


---

## 7. 結論

| 段 | 判定 |
|---|---|
| **C2-1** | 完了。S. coelicolor の占有データ **17 件**を同定。座標つきピーク表があるのは 4 件（GSE165795 / GSE337352 / GSE234438 / GSE64903）。加えて PRJEB50438（HupA）を自前処理。計 **22 トラック**（うち 11 トラック・29 比較が水準の検定に乗った） |
| **C2-2** | 完了。**ピーク内外でメチル化水準に差はない**。水準 29 比較・密度 54 比較すべて FDR 5% で陰性（最小 BH: 水準 0.713、ピーク単位 u 0.425、密度 0.0505）。最も検出力の高い Gbn / HupA では 4mC で約 3%、6mA で約 5-6% の差まで検出できた |
| **C2-3** | 完了。**予測と実測で結論が変わり、向きは本研究に有利**。TSS±100 bp の欠損は再現する（4mC AAGCCCG 0.919 倍 p=0.0025、6mA 0.896 倍 p=0.0075）が、実測占有では一切出ない。TSS 窓を占有の有無で二分すると、欠損は**占有されていない側**に残る |
| **I2-2 / I2-3** | 完了。S. coelicolor の NAP は Gbn（公開ピーク）と HupA（自前処理、暫定）の 2 つ。どちらも陰性。Lsr2 は S. coelicolor には存在せず S. venezuelae のみ |
| **E1-3 / E1-4** | **needs_data**。メチロームと発現が同一株・同一培養で対になった Streptomyces の公開データは 0 件。代替設計 2 案を 5 節に記載 |

### 原稿に書けるか — 書ける

C2-2 / C2-3 / I2-2 は原稿に載せる価値がある。理由は **方向性のある否定**を含むからである:

> タンパク質が DNA を覆ってメチル化を妨げるなら、占有領域のメチル化水準は低下するはずである。
> 水準の検定に十分な部位を含んでいた実測占有 11 トラック（TF 5 + NAP 3 + ScbR/ScbR2 + RokL6）の
> 29 比較のいずれでも、そのような低下は見られなかった。
> 一方で TSS±100 bp の欠損は同じ枠組みで再現し、しかも占有されていない TSS 窓でこそ残る。

これは補足図 11b の FIMO 予測サイトの所見（1.157 倍の軽度の**増加**、予測と逆向き）を
**実測占有で裏づける**ものであり、かつ「予測サイトでしか見ていない」という留保を外せる。
`mem_e67830e233e2` が記録していた 2 つの留保のうち、「FIMO は予測であって ChIP 実測占有ではない」
という留保はこれで解消する。残るもう 1 つの留保（パネル b が検定するのは密度で、
生き残る主張は 6mA のレベル欠損であり尺度がずれている）は、本解析が密度と水準を分けて
報告しているので、そのまま尺度を揃えて書き直せる。

書くときの注意: **帰無を発見に膨らませないこと**。「占有はメチル化に影響しない」ではなく
「3-6% の差があれば検出できた条件で、差は検出されなかった」と書く。
TF トラックは検出力が 8-20% しかないので、そこまで含めて「影響しない」とは書けない。

---

## 8. 成果物

| ファイル | 内容 |
|---|---|
| `tables/C2-1_occupancy_dataset_inventory_261104.tsv` | 見つかった 17 データセットのアクセッション・条件・解像度・座標系・採否 |
| `tables/C2-1_track_summary_261104.tsv` | 構築した 22 トラックの要約 |
| `tables/C2-1_regulator_search_coverage_261104.tsv` | 候補 30 因子それぞれについて探した GEO / ENA のアクセッション空間と該当件数 |
| `tables/C2-2_peak_vs_matched_control_261104.tsv` | 密度 54 比較・水準 29 比較の本表（効果量・p・BH・検出可能だった比・GC 一致の達成度） |
| `tables/C2-3_predicted_vs_measured_261104.tsv` | 予測領域（TSS±100 bp、-300..TSS）を同じパイプラインに通した結果 |
| `tables/C2-3_TSSwindow_by_occupancy_261104.tsv` | TSS 窓を占有の有無で二分した結果 |
| `tables/C2-3_occupancy_interaction_perm_261104.tsv` | 占有の有無の差を GC 層別順列で直接検定した結果 |
| `tables/E1-3_E1-4_strain_panel_feasibility_261104.tsv` | 株パネル候補 7 件の判定と理由 |
| `tables/site_levels_261104.tsv` | 4,761 修飾可能位置の 9 ライブラリプール水準と局所 GC |
| `raw/tracks/`, `raw/tracks_hupa/` | 取得・生成した占有トラック原本 |
| `scripts/occupancy_methylation_261104.py` | 設計と正準値を記録した解析スクリプト |

ピーク RSS は 0.686 GB（`resource.getrusage(RUSAGE_SELF)`、上限 4 GiB に対して十分な余裕）。
パイルアップはモチーフ位置に限定済みで 1 ライブラリ 6.4 MB しかなく、9 本すべてを同時に
保持しても問題にならなかったため、chunksize 読みは不要だった。


---

## 9. ±25 bp 局所 GC 層別（著者指示による追記、2026-11-04）

### 問い

親から著者へ報告した際、「TSS 近傍のメチル化**修飾率**低下は、単に TSS 付近で GC 含量が
下がるからではないか」という指摘を受けた。区別すべき量は 2 つある。

- **(A) モチーフの出現頻度**が下がる … TSS 近傍で GC が落ちるため自明。本稿の既報はこちら。
- **(B) 出現したモチーフの修飾率**が下がる … 3 節の `level_ratio`（4mC AAGCCCG 0.919、6mA 0.896）はこちら。

3 節の対照は**領域 GC**（ピーク幅と同じ長さの窓の GC）一致であって、**部位直近の統制になっていない**。
本プロジェクトでは 6mA のコール率が ±20〜30 bp の短距離 GC に強く依存することが実証済みで、
±250 bp だけで層別した解析が「統制済みに見えて残差が大きい」失敗を実際に起こしている
（2026-10-19、convergent gap の OR が 2.0 → 1.24）。したがって部位ごとの局所 GC で層別し直す。

### 窓の定義（win_excl 型に統一）

**win_excl 型** = 修飾塩基の ±25 bp（51 塩基）から、その部位の認識モチーフ自身の offset を除いたもの。
モチーフ鎖方向で取る。

| 部位 | 認識モチーフ | 修飾塩基の index | 除外 offset | 残る窓幅 |
|---|---|---|---|---|
| 6mA AAGCCCG A1 | HAAGCCCG（8-mer） | 1 | −1..+6 | **43 bp** |
| 4mC AAGCCCG C5 | AAGCCCG（7-mer） | 4 | −4..+2 | 44 bp |
| 4mC TGGCCGGC C5 | TGGCCGGC（8-mer） | 4 | −4..+3 | 43 bp |

6mA A1 の 43 bp は本プロジェクトの既存定義（`mem_9b2625792a8d`）と一致する。
flank 型（モチーフ外側に 25 bp ずつ、50 bp）は使っていない。比較のため同じ win_excl の
**±250 bp 版**（501 塩基からモチーフを除外）も全解析で並走させた。

### 統制の必要性は実在した

TSS 近傍の部位は実際に局所 GC が低い（win_excl ±25 bp）:

| モチーフ | 近傍（≤100 bp）の平均 GC | 遠方（>100 bp）の平均 GC |
|---|---|---|
| 6mA AAGCCCG | 0.6665 | 0.7131 |
| 4mC AAGCCCG | 0.6684 | 0.7178 |
| 4mC TGGCCGGC | 0.6869 | 0.7329 |

著者の懸念は根拠のあるものだった。

### 1. 部位単位・GC 五分位層別（Mantel–Haenszel 型）

近傍（最近接 Jeong 実験 TSS から ≤100 bp）対 遠方（>100 bp）の修飾率の比を、局所 GC 五分位の
層内で取り、層サイズ重み `w_k = n1k·n0k/(n1k+n0k)` で統合した。p 値は層内で近傍ラベルを
入れ替える層別順列（2,000 回）、CI は層×群内ブートストラップ（2,000 回）。
| motif | stratifier | n_near | n_far | n_bins | ratio | ci_lo | ci_hi | p_perm |
|---|---|---|---|---|---|---|---|---|
| 4mC_TGGCCGGC | none (crude near/far) | 84 | 2009 | 1 | 1.0145 |  |  |  |
| 4mC_TGGCCGGC | win_excl ±25 bp, GC quintile | 84 | 2009 | 5 | 1.0001 | 0.9704 | 1.0312 | 0.9995 |
| 4mC_TGGCCGGC | win_excl ±250 bp, GC quintile | 84 | 2009 | 5 | 1.0092 | 0.9764 | 1.0421 | 0.6807 |
| 4mC_AAGCCCG | none (crude near/far) | 115 | 1219 | 1 | 0.9247 |  |  |  |
| 4mC_AAGCCCG | win_excl ±25 bp, GC quintile | 115 | 1219 | 5 | 0.9051 | 0.8517 | 0.9606 | 0.0005 |
| 4mC_AAGCCCG | win_excl ±250 bp, GC quintile | 115 | 1219 | 5 | 0.9108 | 0.8585 | 0.9628 | 0.0015 |
| 6mA_AAGCCCG | none (crude near/far) | 117 | 1217 | 1 | 0.8449 |  |  |  |
| 6mA_AAGCCCG | win_excl ±25 bp, GC quintile | 117 | 1217 | 5 | 0.8883 | 0.8168 | 0.962 | 0.003 |
| 6mA_AAGCCCG | win_excl ±250 bp, GC quintile | 117 | 1217 | 5 | 0.88 | 0.8105 | 0.9506 | 0.0005 |

**効果は消えない。** 4mC AAGCCCG は素の 0.925 から **0.905**（CI 0.852–0.961、p = 0.0005）へ、
6mA は素の 0.845 から **0.888**（CI 0.817–0.962、p = 0.0030）へ動くだけである。
TGGCCGGC は層別前から効果が無く（1.015）、層別後 1.000（p = 1.00）で変わらない。

±250 bp 版は 4mC 0.911 / 6mA 0.880 で、±25 bp 版とほぼ同じ。
**この量については、convergent gap で起きたような「±250 bp では統制済みに見えて ±25 bp で崩れる」
現象は再現しなかった**。6mA の短距離 GC 依存性は修飾率の絶対水準を動かすが、
近傍対遠方の比には効いていない。

### 2. 層ごとの内訳（6mA、win_excl ±25 bp）
| bin | n_far | n_near | mean_level_far | mean_level_near | med_gc_far | med_gc_near | ratio |
|---|---|---|---|---|---|---|---|
| 0 | 296 | 53 | 0.3521 | 0.3222 | 0.6279 | 0.6047 | 0.9149 |
| 1 | 277 | 28 | 0.3765 | 0.3384 | 0.6977 | 0.6744 | 0.8987 |
| 2 | 150 | 14 | 0.4005 | 0.3055 | 0.7209 | 0.7209 | 0.7628 |
| 3 | 277 | 15 | 0.4294 | 0.39 | 0.7442 | 0.7442 | 0.9082 |
| 4 | 217 | 7 | 0.4557 | 0.4052 | 0.814 | 0.814 | 0.8892 |

**5 層すべてで比が 1 を下回る**（0.763–0.915）。単一の層が駆動しているのではない。
4mC AAGCCCG も同様に 5 層すべてで 0.860–0.946。

### 3. 局所 GC を完全一致させた最強版

win_excl ±25 bp の GC は 43（または 44）分の整数なので離散値である。五分位ではなく
**GC の値そのもの**で層別すれば、層内残差 GC は厳密にゼロになる。
| motif | stratifier | n_strata_used | n_near_used | ratio | p_perm | mean_residual_GC_diff |
|---|---|---|---|---|---|---|
| 4mC_TGGCCGGC | win_excl ±25 bp, EXACT GC value | 14 | 83 | 0.9979 | 0.91504 | 0.0 |
| 4mC_TGGCCGGC | win_excl ±250 bp, GC decile | 9 | 84 | 1.006 | 0.7911 | -0.00209 |
| 4mC_AAGCCCG | win_excl ±25 bp, EXACT GC value | 18 | 115 | 0.8926 | 0.0005 | 0.0 |
| 4mC_AAGCCCG | win_excl ±250 bp, GC decile | 10 | 115 | 0.9084 | 0.0005 | -0.0004 |
| 6mA_AAGCCCG | win_excl ±25 bp, EXACT GC value | 16 | 117 | 0.9051 | 0.0075 | 0.0 |
| 6mA_AAGCCCG | win_excl ±250 bp, GC decile | 10 | 117 | 0.8805 | 0.001 | -0.0007 |

層内残差 GC を 0.00000 に落としても **4mC 0.893（p = 0.0005）、6mA 0.905（p = 0.0075）**。
局所 GC で説明できる成分を完全に取り除いても欠損は残る。これが (B) が (A) の派生ではないことの
最も強い証拠である。

### 4. 既存正準値との同一性 — **同じ量であり、新規所見ではない**

本レーンの近傍/遠方比は、本プロジェクトの「プロモーター近傍 6mA 欠乏」と**同じ量の再測定**である。
正準と同じユニバース（`both_A` = AAGCCCG の A1 と A2 の両方）で計算すると:

| ユニバース | 近傍部位数 | 遠方部位数 |
|---|---|---|
| 本レーン `both_A` | **235** | 2,433 |
| 正準（SuppTable22b, jeong_measured, corrected） | **235** | 2,429 |

近傍部位数が 235 で完全に一致する。分割は同一で、遠方が 4 部位違うのは被覆の扱いによる。

| 条件（培養単位 n = 9） | 比 | 95% CI | p |
|---|---|---|---|
| **正準（現行）** SuppTable22a: all reads, five local-GC quantiles, ±100 | **0.9105** | 0.8955–0.9257 | 1.12e-06 |
| 本レーン `both_A`, win_excl ±25 bp, 10 分位 | 0.9002 | 0.8869–0.9138 | 2.12e-07 |
| 本レーン `both_A`, win_excl ±25 bp, 5 分位 | 0.8859 | 0.8730–0.8990 | 5.95e-08 |
| 本レーン `both_A`, win_excl ±250 bp, 5 分位 | 0.8788 | 0.8665–0.8912 | 2.60e-08 |

10 分位版の CI（0.8869–0.9138）は正準の CI（0.8955–0.9257）と重なる。
**本レーンの 0.896 / 0.919 は独立なヌルでの再現であって、新規所見ではない。**
原稿では「再現」として扱い、新しい効果として書かないこと。

> **依頼文の前提の訂正**: 依頼では正準値を「比 0.9212、CI 0.9086–0.9337、p = 1.7e-6」としていたが、
> この値は本プロジェクト内で既に差し替えられている。0.9212 が載っているのは
> `15_paper_figures/tables/supplementary/**superseded**/SuppTable16_statistical_test_summary_261016.tsv`
> だけで、現行の `SuppTable22a_promoter_deficit_conditions_261020.tsv`（2026-10-23 更新）は
> **0.9105（CI 0.8955–0.9257、p = 1.12e-06）**を載せ、注記で「0.8052/0.8945/0.8475/0.9212 は
> どちらの座標規約でもパイルアップから再現できず、superseded 世代のフィルタ済み表に由来すると見られる」
> と明記している。以後 0.9212 を正準として引かないこと。

### 5. 占有の有無で分けた版（3 節の 0.894 / 0.873 に同じ層別を適用）

3 節では TSS±100 bp の窓を実測占有（19 トラックの和集合、ゲノムの 23.58%）と重なるかで二分し、
欠損が**占有されていない側に残る**ことを示した。同じ分割に局所 GC 層別を掛ける。
参照は (i) 遠方部位すべて、(ii) 占有状態を揃えた遠方部位（占有そのものを定数化した厳しい版）の 2 通り。
| motif | subset | reference | gc_window | n_near | n_far | n_strata | ratio | p_perm |
|---|---|---|---|---|---|---|---|---|
| 4mC_AAGCCCG | near OCCUPIED | all far sites | ±25 win_excl | 38 | 1219 | 5 | 0.9328 | 0.09945 |
| 4mC_AAGCCCG | near OCCUPIED | all far sites | ±250 win_excl | 38 | 1219 | 4 | 0.9354 | 0.11294 |
| 4mC_AAGCCCG | near OCCUPIED | far sites of same occupancy status | ±25 win_excl | 38 | 325 | 5 | 0.9394 | 0.14393 |
| 4mC_AAGCCCG | near OCCUPIED | far sites of same occupancy status | ±250 win_excl | 38 | 325 | 4 | 0.9438 | 0.18291 |
| 4mC_AAGCCCG | near UNOCCUPIED | all far sites | ±25 win_excl | 77 | 1219 | 5 | 0.8921 | 0.0005 |
| 4mC_AAGCCCG | near UNOCCUPIED | all far sites | ±250 win_excl | 77 | 1219 | 5 | 0.9019 | 0.003 |
| 4mC_AAGCCCG | near UNOCCUPIED | far sites of same occupancy status | ±25 win_excl | 77 | 894 | 5 | 0.8885 | 0.001 |
| 4mC_AAGCCCG | near UNOCCUPIED | far sites of same occupancy status | ±250 win_excl | 77 | 894 | 5 | 0.8972 | 0.0025 |
| 6mA_AAGCCCG | near OCCUPIED | all far sites | ±25 win_excl | 38 | 1217 | 5 | 0.9396 | 0.33633 |
| 6mA_AAGCCCG | near OCCUPIED | all far sites | ±250 win_excl | 38 | 1217 | 4 | 0.9426 | 0.3973 |
| 6mA_AAGCCCG | near OCCUPIED | far sites of same occupancy status | ±25 win_excl | 38 | 327 | 5 | 0.9418 | 0.37031 |
| 6mA_AAGCCCG | near OCCUPIED | far sites of same occupancy status | ±250 win_excl | 38 | 327 | 4 | 0.9237 | 0.25737 |
| 6mA_AAGCCCG | near UNOCCUPIED | all far sites | ±25 win_excl | 79 | 1217 | 5 | 0.8645 | 0.002 |
| 6mA_AAGCCCG | near UNOCCUPIED | all far sites | ±250 win_excl | 79 | 1217 | 5 | 0.8526 | 0.0005 |
| 6mA_AAGCCCG | near UNOCCUPIED | far sites of same occupancy status | ±25 win_excl | 79 | 890 | 5 | 0.8629 | 0.0015 |
| 6mA_AAGCCCG | near UNOCCUPIED | far sites of same occupancy status | ±250 win_excl | 79 | 890 | 5 | 0.8607 | 0.0025 |
| 4mC_TGGCCGGC | near OCCUPIED | all far sites | ±25 win_excl | 25 | 2009 | 5 | 1.053 | 0.18441 |
| 4mC_TGGCCGGC | near OCCUPIED | all far sites | ±250 win_excl | 25 | 2009 | 4 | 1.0617 | 0.14243 |
| 4mC_TGGCCGGC | near OCCUPIED | far sites of same occupancy status | ±25 win_excl | 25 | 392 | 5 | 1.0692 | 0.10445 |
| 4mC_TGGCCGGC | near OCCUPIED | far sites of same occupancy status | ±250 win_excl | 25 | 392 | 4 | 1.0781 | 0.05597 |
| 4mC_TGGCCGGC | near UNOCCUPIED | all far sites | ±25 win_excl | 59 | 2009 | 5 | 0.9776 | 0.35082 |
| 4mC_TGGCCGGC | near UNOCCUPIED | all far sites | ±250 win_excl | 59 | 2009 | 5 | 0.9868 | 0.6062 |
| 4mC_TGGCCGGC | near UNOCCUPIED | far sites of same occupancy status | ±25 win_excl | 59 | 1617 | 5 | 0.9727 | 0.26737 |
| 4mC_TGGCCGGC | near UNOCCUPIED | far sites of same occupancy status | ±250 win_excl | 59 | 1617 | 5 | 0.9787 | 0.3913 |

**3 節の結論は層別後も変わらない**。占有されていない近傍部位では
4mC 0.892（p = 0.0005）/ 6mA 0.865（p = 0.0020）、占有を定数化した参照でも 0.889 / 0.863。
占有されている近傍部位では 4mC 0.933（p = 0.099）/ 6mA 0.940（p = 0.336）で有意に達しない。
TGGCCGGC はどちらの側でも効果なし。

### 6. 判定

1. **効果は消えない。** 局所 GC（win_excl ±25 bp）で層別しても、GC 値を完全一致させても、
   近傍の修飾率低下は残る（4mC AAGCCCG 0.893–0.905、6mA 0.888–0.905、いずれも p ≤ 0.0075）。
   著者の懸念（(B) は (A) の派生ではないか）は、**否定された**。
2. **ただし新規所見ではない。** これは既存の「プロモーター近傍 6mA 欠乏」（現行正準 0.9105）の
   独立なヌルでの再現である。近傍部位数 235 が正準と一致することで同一量であることを確認した。
3. **±250 bp でも ±25 bp でも結論は同じ**で、この量については短距離 GC の残差問題は起きていない。
   convergent gap の失敗（2.0 → 1.24）は再現しなかった。
4. **4mC TGGCCGGC には一貫して効果が無い**（1.000–1.006、p ≥ 0.79）。欠損は AAGCCCG 系に限られる。

### 追加した成果物

| ファイル | 内容 |
|---|---|
| `tables/C2-3b_localGC_stratified_siteunit_261104.tsv` | 部位単位・GC 五分位層別（±25 / ±250、MH 統合・順列 p・ブートストラップ CI） |
| `tables/C2-3b_localGC_stratified_cultureunit_261104.tsv` | 培養単位 n = 9 版（5 分位・10 分位） |
| `tables/C2-3b_exact_localGC_stratification_261104.tsv` | 局所 GC の値そのもので層別した版（層内残差 GC = 0） |
| `tables/C2-3b_bothA_canonical_comparison_261104.tsv` | 正準と同じ `both_A` ユニバースでの再計算 |
| `tables/C2-3b_occupancy_split_stratified_261104.tsv` | 占有の有無で分けた層別版（参照 2 通り） |
| `tables/C2-3b_perstratum_6mA_261104.tsv` | 6mA の層ごとの内訳 |
| `tables/site_localGC_261104.tsv` | 全 4,761 部位の win_excl ±25 / ±250 GC、TSS 距離、占有、水準 |

ピーク RSS 0.592 GB。
