# 著者コメント #12 — アデニン MTase 候補 5 件は網羅できているか

解析日 2026-10-24 / `analysis/113_mtase_completeness_261024/` / conda env `epi-rna`
成果物: `tables/mtase_completeness_261024.tsv`, `scripts/mtase_completeness_261024.sh`

**結論（先に）: 候補 5 件は、アノテーション名に依存しない独立走査でも増減しない。原稿の
「5 遺伝子／うち 3 件」という数は変更不要。ただし入口の 129 件は RefSeq の product 名に対する
文字列一致で定義されており、原稿はそれを書いていない。文字列一致で**実体のある取りこぼしが 1 件
あり（SC_RS19860）**、それは S4 相当で落ちるため最終 5 件は変わらないが、漏斗の中間値 21 / 12 は
22 / 13 になる。既存走査が 129 外として挙げたもう 1 件 SC_RS36415 は境界アーティファクトで
候補ではない（§2.4 で根拠を示す）。** 348 aa の床は導出表そのものは実在して再現できたが、`value_registry.csv` の
出所欄が導出元ではなく適用先を指しており、導出表はローカルディスク上に存在しない。

---

## 1. 129 件がどう定義されたか（追跡完了）

| 層 | 現物 |
|---|---|
| 原稿が引く数 | 「129 SAM-dependent methyltransferases」(EN L301, Supplementary-Materials-List L9) |
| 入力テーブル | `107_MTase_candidate_triage/data/M145_MTase_129.faa`（128 配列。1 件は pseudogene） |
| その上流 | `33_GCCGGC_RM_identification/tables/M145_MTase_complete_catalog.tsv`（129 行） |
| さらに上流 | `30_CCGG_MTase_paradox/tables/all_methyltransferases_GFF.tsv`（129 行） |
| 生成コード | `30_CCGG_MTase_paradox/analyze_CCGG_paradox.py` L383-401 |
| 入口の規則 | RefSeq GFF (`GCF_000203835.1/genomic.gff`) の CDS の `product` を小文字化し、**12 語の部分一致** |

キーワード: `methyltransferase`, `methylase`, `dna methylation`, `modification methylase`, `dcm`,
`dam`, `hsdm`, `n-6 adenine`, `c-5 cytosine`, `n-4 cytosine`, `sam-dependent methyltransferase`,
`restriction-modification`。

本日この規則を GFF に再適用して **129 件を正確に再現した**（8,199 CDS 中 129、全件 NC_003888.3 上、
SCP1/SCP2 からは 0 件）。つまり入口は **product 名の文字列一致であり、ドメイン構成ではない**。
名前が付いていない ORF、`DUF…`、`hypothetical protein` は、MTase ドメインを持っていても
構造的に入口で落ちる。

この規則の品質も測った: 129 件のうち **9 件は Pfam のメチル基転移酵素系ファミリーを
gathering threshold で 1 つも持たない**（偽陽性側）。逆方向は §2。

> 同じ問題を、このプロジェクトは一度解決している。`52_shielded_exposed_boundary/scripts/`
> `H29c_pfam_regulatory_gene_selection.py` の冒頭は、制御遺伝子集合を keyword(n=449) から
> Pfam ドメインへ置き換えた理由を「キーワード表が恣意的で ~1.6% の偽陽性を入れ、
> 非標準アノテーションの遺伝子（PhoP = "PhoX family protein"）を落とすため」と書いている。
> 論文の最も荷重の大きい数（62 Exposed / 989 Shielded / 1,051）には適用され、
> MTase カタログには適用されていない。

---

## 2. 独立な走査（アノテーション名に依存しない）

使えたもの: HMMER 3.3（`hmmsearch`, `hmmfetch`）、ローカル Pfam-A（**25,545 モデル**、
`/Users/okaban/bioinfo/data/pfam/Pfam-A.hmm`）、RefSeq `protein.faa` と `genomic.gff`。
InterProScan・CDD のローカル実行環境は無い（**この 2 つでは確認していない**）。

### 2.1 走査 A — 全染色体タンパク質 × Pfam-A 全ファミリー（gathering threshold）

対象: NC_003888.3 の CDS 7,602 件のうち protein_id を持つ **7,553 タンパク質**（重複配列を除いた数）。
`hmmsearch --cut_ga`、12 スレッド。**12,810 ドメインヒット / 6,305 タンパク質 / 3,058 ファミリー**。

- メチル基転移酵素系ファミリーを持つタンパク質: **135 件**（既存 129 のうち 120 件 + **既存 129 に無い 15 件**）
- **DNA amino-MTase ファミリー（N6_Mtase PF02384 / N6_N4_Mtase PF01555 / Eco57I PF07669 /
  MmeI_Mtase PF20473）を持つのは、ゲノム全体で 5 件だけ** — SC_RS17645, SC_RS36625, SC_RS19670,
  SC_RS28835, SC_RS35335。**全件が既存 129 に含まれる。追加はゼロ。**
- C5 系（DNA_methylase PF00145）は SC_RS19770 と SC_RS36410 の 2 件。いずれも既存 129 内。

### 2.2 走査 B — 寛容な E 値での再走査（床下の分岐した酵素対策）

Pfam の gathering threshold は保守的なので、コア DNA-MTase モデル 10 本を取り出し、
`-E 10 --domE 10` で **アセンブリ全体 7,872 タンパク質（SCP1・SCP2 を含む）** に当て直した。
有意域（E < 1e-5）のヒットは §2.1 と同じ 7 件のみ。E ≤ 10 まで下げて現れる 129 外の遺伝子は
**SC_RS19860 1 件のみ（Eco57I, E = 0.31 — 有意ではない）**。**プラスミド側のヒットは 0 件**
（M145 は SCP1/SCP2 を持たないので元々対象外だが、参照アセンブリには含まれるため確認した）。
`Methyltransf_26` はこの Pfam ビルドに存在せず、取得できなかった。

### 2.3 走査 C — HMM を使わない配列パターン走査

7,553 タンパク質に motif IV と motif I を直接当てた。

| パターン | 保有タンパク質数 | うち既存 129 に無い |
|---|---|---|
| motif IV 正準 `[DNS]PP[YFW]` | 43 | 34 |
| motif IV 寛容 `[A-Z]PP[YFW]`（三分類の S2 が使った形） | 410 | — |
| 正準 motif IV ＋ 緩い motif I | 20 | 16 |

単独では情報量が無い（硝酸還元酵素・アセチル CoA 合成酵素などが並ぶ）。意味を持つのは
**motif IV がメチル基転移酵素系 Pfam エンベロープの内側に落ちるか**という S3 相当の条件で、
これをゲノム全体に当てると **18 タンパク質**（既存 129 内 16、**129 外 2**）。

### 2.4 既出の走査との関係（重複していない）

2026-10-06 に `unannotated_ORF_MTase_scan_261006.csv` が別の角度から同じ問いを見ている:
NC_003888.3 の 6 フレーム翻訳 70,506 ORF、うち注釈 CDS との重複 20% 未満の 2,686 ORF に
**Pfam-A 全体を gathering threshold で当てて、いかなるファミリーのヒットも 0 件**。
つまり「注釈 CDS の外」は既に閉じている。本日の走査が閉じたのは **「注釈 CDS の内側だが
product 名が一致しない」** という、そちらが見ていなかった側である。2 つで入口は両方向とも塞がる。

**その走査で 129 外として挙がった遺伝子は 2 件**（`SC_RS19860` と `SC_RS36415`）で、両方を本日
個別に突き合わせた。`SC_RS19860` は §3.1 で扱う実体のある取りこぼし。**`SC_RS36415`
（adenine nucleotide alpha hydrolase family protein, 320 aa, TPM 0.00/0.10/4.03）は境界
アーティファクトである**: 当該 6 フレーム ORF036516 の `DNA_methylase` ドメインはゲノム座標
7,618,135-7,618,380 に落ち、これは `SC_RS36410`（7,617,112-7,618,383, DNA cytosine
methyltransferase, 既存 129 内）の内側で、`SC_RS36415` とは開始座標 7,618,380 で **3 bp
重なるだけ**である。本日の走査では `SC_RS36415` は **gathering threshold で Pfam ドメインを
1 つも持たず**、コアモデルの E ≤ 10 走査でもヒットせず、S3 相当にも届かない。motif IV は
`APPF@3` が 1 つあるがエンベロープ外（そもそもエンベロープが無い）。**DNA MTase 候補ではない。**

---

## 3. 新走査で出て既存 129 に無いもの

### 3.1 S3 相当まで到達する 2 件

| locus_tag | 長さ | product（RefSeq） | Pfam（methyl 系, i-E） | motif IV | TPM T1/T2/T3 | 判定 |
|---|---|---|---|---|---|---|
| **SC_RS19860** | 505 aa | DUF7059 domain-containing protein | **MTS PF05175, 6.7e-12**, envelope 152-293 | **NPPF@238（正準、エンベロープ内）** | **74.69 / 42.21 / 32.91** | S4 で脱落（MTS は DNA amino-MTase ファミリーではない） |
| SC_RS35585 | 359 aa | aromatic ring-hydroxylating dioxygenase subunit alpha | VanA_C PF19112, 1.6e-38 | CPPF@206 | 0.46 / 1.00 / 5.61 | S4 で脱落。そもそも「触媒ドメイン」をどう割り付けるかで S3 に入るかが変わる |

**SC_RS19860 が本件の実体である。** これは product 名が `DUF7059 domain-containing protein` である
ためだけに入口で落ちている。MTS (PF05175) 保有者はゲノム全体で 19 件あり、**129 に入っていないのは
この 1 件だけ**。しかも T1 で 74.7 TPM（全遺伝子の第 85 百分位）と、候補 5 件のどれよりも高く発現している。

ただし **最終 5 件は変わらない**: 既存三分類で MTS を触媒ドメインとする 3 件（SC_RS07160 377 aa,
SC_RS05730 214 aa, SC_RS30050 195 aa）はいずれも S4 で
「catalytic domain MTS is not a DNA amino-MTase family」として落ちている。SC_RS19860 も同じ段で落ちる。

**影響は漏斗の中間値に出る**: S2 は 21 → **22**、S3 は 12 → **13**。S4 = 5、S5 = 3 は不変。

### 3.2 S3 に届かない 129 外の methyl 系ドメイン保有 13 件

methionine synthase (SC_RS10295)、VOC family ×3、FtsJ rRNA MTase (SC_RS10920)、
precorrin-6A synthase (SC_RS37065)、MnmA (SC_RS29625)、NRPS (SC_RS40605)、MerR (SC_RS40675)、
Gfo/Idh/MocA (SC_RS40615)、DUF1295 (SC_RS31735)、および **MGMT family protein 2 件**
（SC_RS28075 135 aa, SC_RS40505 94 aa, DNA_binding_1 PF01035）。MGMT は DNA からアルキル基を
外す修復酵素で SAM 依存 DNA MTase ではない。**DNA MTase 候補になるものは 1 件も無い。**

これとは別に、§2.4 の 2026-10-06 走査が 129 外として挙げた **SC_RS36415**（320 aa）も
突き合わせ済み。Pfam ドメインを GA 閾値で 1 つも持たず、隣接する SC_RS36410 のドメインとの
3 bp 境界重複で名前が出ていただけなので、候補ではない。TSV にも 1 行として収録した。

---

## 4. 348 aa のサイズ床 — 出所の追跡

**導出元は実在し、再現できた。ただしローカルディスク上には無い。**

- 導出表: `MTase_size_vs_recognition_span_261006.csv`（artifact 24833404-5d3c-4ec5-8855-0f8e1538360f）。
  REBASE × reviewed UniProt の結合で、認識配列と酵素長の両方が分かる **65 酵素**（6mA 45 / 4mC 20）。
- 本日読み直して再計算: 認識配列の総長が 7 bp 以上は **6 酵素**、最小 **348 aa = M.BglI**、
  中央値 **436.5 aa**、最大 912 aa（M.BslI）。canon の記述と完全一致。

| 酵素 | 化学種 | 長さ | 認識配列 | 総長 | 連続最長 |
|---|---|---|---|---|---|
| M.BglI | 4mC | **348** | GCCNNNNNGGC | 11 | 3 |
| M.MamI | 6mA | 362 | GATNNNNATC | 10 | 3 |
| M.SfiI | 4mC | 421 | GGCCNNNNNGGCC | 13 | 4 |
| M.EcaI | 6mA | 452 | GGTNACC | 7 | 3 |
| M.MwoI | 4mC | 668 | GCNNNNNNNGC | 11 | 2 |
| M.BslI | 4mC | 912 | CCNNNNNNNGG | 11 | 2 |

**床の弱点を 3 つ名指しする（いずれも原稿が現在書いていない）:**

1. **6 件の最小値である。** n = 6 の標本最小値であって下界ではない。65 酵素全体の最小は 194 aa。
2. **6 件はすべて分断認識配列である。** 65 酵素の中で **連続スパンの最大は 6 bp**（canon の
   「6 bp を超えるものは存在しない」はこの意味で正しい）。つまり 348 aa は
   **「NNNNN を挟む二分割配列を読む酵素」から取った値**で、AAGCCCG の 7 bp 連続・
   TGGCCGGC の 8 bp 連続という**構造的に別の問題**に外挿している。現行原稿本文はこの点に
   触れておらず、“the floor is a heuristic, not a measurement on this organism” とだけ書く。
   旧版 `.bak_261016_supp/Manuscript_EN_NAR_numbered_inline.md` L142 の
   “the smallest REBASE-reviewed enzyme recognising an **unbroken** span of 7 bp or more” は
   6 件が全て分断配列である点で現物と食い違っていた（M.EcaI でさえ GGT**N**ACC）。
   **現行原稿には無い。復活させないこと。**
3. **出所の記録が循環している。** `value_registry.csv` の `MTASE.SIZE_FLOOR.AA` は
   `source_artifact = 107_MTase_candidate_triage/tables/MTase_candidate_elimination_261008.csv`、
   `source_locator = 列 length_aa の閾値` を指す。これは **床を適用した表**であって
   **床を導出した表**ではない。導出表を指しているのは `results_inventory.csv` R7-11 の
   `MTase_size_vs_recognition_span_261006.csv` のみで、**そのファイルはディスク上に存在しない**
   （artifact ストアにのみある）。

---

## 5. 4mC 側に同じ問いを当てる

原稿の現在の記述:

- L226（Discussion 冒頭）: “…with no 4mC methyltransferase candidate in the genome.”
- L248: “The enzyme responsible for either chemistry remains unidentified. … The 4mC activity is
  unassigned at both TGGCCGGC and AAGCCCG position 5.”
- L305: “HPPF places SC_RS17645 on the adenine side, leaving the 4mC methyltransferase unidentified…”
- L258（Conclusion）: “…the cytosine methyltransferase remains unidentified.”

**網羅性の観点での判定:**

1. **列挙の側は、本日の走査で閉じた。** Pfam には 4mC 専用ファミリーが無く、N4 と N6 は
   同じ PF01555 (`N6_N4_Mtase`) / PF02384 を共有する。したがって **4mC 候補のプールは
   アデニン側と同一の 5 件**であり、その 5 件が網羅的であることは §2 で示した。
   C5 系 2 件（SC_RS19770, SC_RS36410）はフォールドが違うので 4mC 候補にならない。
   既存の `4mC_MTase_candidates_261006.csv`（8 遺伝子、個別の除外理由つき）はこの 5+2+1 を
   すでに個別評価しており、本日の走査と矛盾しない。
2. **したがって 4mC 側に「別の網羅性問題」は無い** — ただしそれは
   「4mC 候補は独立に数えられていない。アデニン候補と同じ 5 件を化学種で評価しているだけ」
   だからである。列挙が閉じたのは今日であり、**原稿にはその走査が書かれていない。**
3. **L226 の言い方は、支えられる主張より強い。** “no 4mC methyltransferase candidate in the
   genome” は「ゲノムに候補が存在しない」と読めるが、解析が示しているのは
   「ゲノムの amino-MTase 5 件のうち 4mC 化学種に割り付けられるものが無い」である。
   5 件は存在し、そのうち 2 件（SC_RS36625 / SC_RS19670）は `4mC_MTase_candidates_261006.csv` で
   “UNDETERMINED (motif IV says 6mA; discriminant above 4mC reference minimum)” と判定されている。

---

## 6. 報告（6 列）

| 対象箇所 | 原文または現物 | 問題の種類 | 理由 | 最小限の修正案 | 著者確認 |
|---|---|---|---|---|---|
| EN L301 / Supplementary-Materials-List L9 / Figure 4c | “out of 129 SAM-dependent methyltransferases filtered by sequence retrievability, …” | 方法の未記載 | 129 の入口は `product` 名への 12 語の文字列一致（`analyze_CCGG_paradox.py` L383-401）。読者は入口がドメインベースだと読む。本日の再現で 129 を正確に再現済み | 「129」に由来を 1 句足す: “…out of 129 genes whose RefSeq product name matches a methyltransferase keyword list (Methods)” | **要** |
| Figure 4c（漏斗図、`panel_4d.py` → `panel_F4c.png`、データ `data_261013/v7295f465_F4d_funnel.csv`） | 棒ラベル `SAM-dependent MTases 129 / protein sequence available 128 / motif IV present 21 / motif IV in catalytic domain 12 / amino-MTase family 5 / ≥ 348 aa size floor 3` | 数の食い違い（最終値は不変） | 名前一致の入口が SC_RS19860（DUF7059, 505 aa, MTS PF05175 i-E 6.7e-12, 正準 motif IV NPPF@238 がエンベロープ内, TPM 74.7/42.2/32.9）を取りこぼす。ドメイン基準の入口では 3 段目 = 22、4 段目 = 13。5 段目 = 5、6 段目 = 3 は不変。中間値は本文には無く、この図にのみ出る | どちらかを選ぶ: (a) 中間値をそのままにし、1 段目のラベルを “SAM-dependent MTases (product-name match)” に替える、(b) 入口を Pfam に置換して 22/13 に直す。いずれも `v7295f465_F4d_funnel.csv` の 1 列の変更で済み、タイルを作り直してから `compose(4)` のみ呼ぶ | **要**（どちらを採るかは著者判断） |
| Methods / Supplementary Methods | 網羅性の検証に関する記述が無い | 欠落 | 査読者の自然な質問（「名前の付いていない ORF は？」）に答える材料が、解析としては 2 件（2026-10-06 の 6 フレーム走査、本日の全 Pfam 走査）揃っているのに原稿に出ていない | Supplementary Methods に 2-3 文: 6 フレーム 70,506 ORF × Pfam-A（GA）で注釈外 0 件、染色体 7,553 タンパク質 × Pfam-A（GA）で DNA amino-MTase ファミリーは 5 件のみ | **要** |
| EN L301（“a size floor”）、Supplementary Figure 11 凡例（“the 348 aa size floor”）、Figure 4c 棒ラベル（“≥ 348 aa size floor”） | 現行の 3 箇所いずれも **348 が何から出た数かを書いていない**。本文は “the floor is a heuristic, not a measurement on this organism” と限定するのみ | 根拠が原稿内で追えない | 導出（REBASE × reviewed UniProt の 65 酵素、認識スパン 7 bp 以上の 6 酵素の最小値）は実在し本日再現したが、原稿のどこにも書かれていないので読者は 348 を評価できない。なお “the smallest REBASE-reviewed enzyme recognising an **unbroken** span of 7 bp or more” という旧表現は **現行の原稿には無く**、`.bak_261016_supp/` の旧版にのみ残る（6 酵素はすべて分断配列で、65 酵素中の連続最長は 6 bp なので、この旧表現は現物と食い違っていた。復活させないこと） | Methods か Supplementary Figure 11 凡例に 1 文: “348 aa is the smallest of six REBASE×reviewed-UniProt enzymes whose recognition site spans ≥7 bp (median 436.5 aa, n = 65 enzymes with both site and length known); all six recognise bipartite sites, so the floor is an extrapolation to a contiguous 7-8 bp target.” | **要** |
| `Writing/value_registry.csv` の `MTASE.SIZE_FLOOR.AA` | `source_artifact = …/MTase_candidate_elimination_261008.csv`, `source_locator = 列 length_aa の閾値` | 出所の循環 | 床を**適用した**表を指していて**導出した**表を指していない。導出表 `MTase_size_vs_recognition_span_261006.csv` は artifact ストアにのみ存在しディスク上に無い | 導出表をディスクに書き出し、`source_artifact` をそれに差し替える（原稿本文の変更ではない） | 不要（記録の修正） |
| Discussion L226 | “…and no 4mC methyltransferase candidate in the genome.” | 断定が強い | 解析が支えるのは「ゲノムの amino-MTase 5 件のうち 4mC 化学種に割り付けられるものが無い」。候補が 0 件なのではなく、割付が付かない。うち 2 件は既存解析で UNDETERMINED | “…and no methyltransferase in the genome assignable to 4mC chemistry.” | **要** |
| 候補 5 件・生存 3 件という数そのもの | EN L301, L248, Figure 4c, Supplementary Figure 11 | **問題なし** | アノテーション非依存の全 Pfam 走査（GA）でも寛容走査（E≤10, プラスミド込み）でも、DNA amino-MTase ファミリー保有は同じ 5 件。追加ゼロ | 変更不要 | 不要 |

---

## 7. 確認できなかったこと

- **InterProScan / CDD による独立確認**: ローカルに実行環境が無く、**この方法では確認していない**。
  本日の判定は Pfam-A（25,545 モデル、HMMER 3.3）単独による。
- **SC_RS19765**（RefSeq pseudogene、102 bp）: 翻訳産物が無く、本日の走査でも評価不能。
  既存三分類の S1 脱落と同じ扱い。
- **Pfam `Methyltransf_26`**: このローカル Pfam ビルドに存在せず取得できなかった
  （他 9 本のコアモデルは取得・実行済み）。
- **SC_RS19860 の機能**: MTase ドメインを持ち高発現だが、本解析は
  「DNA amino-MTase ファミリーではない」としか言っていない。実際に何をしている酵素かは追っていない。

## 8. 再現

```
cd 11_epigenome_integration/analysis/113_mtase_completeness_261024
bash scripts/mtase_completeness_261024.sh     # hmmsearch 2 本（12 スレッド）
```
TPM は `03_quantification/data/featureCounts/featureCounts_M145_s2.txt` から実行時に算出
（length 正規化 → ライブラリ合計で割って 1e6。T1 = 1-1,1-2,1-3 / T2 = 2-1,2-3,2-4 /
T3 = 3-2,3-3,3-4 の平均）。全遺伝子中央値 9.14 TPM を再現。
