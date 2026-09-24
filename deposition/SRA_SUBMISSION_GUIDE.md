# SRA 提出ガイド（2026-09-10 作成）

すべてのメタデータは埋まっている。**アカウント操作だけが残っている**。
所要は入力 30 分 ＋ アップロード（**14.3 GB** = BAM 6.18 GB ＋ FASTQ 8.13 GB、回線次第で数時間）。

## 0. pod5 の可否 — **2026-09-22 に SRA ヘルプデスクから回答済み（照会は不要）**

> "SRA does not accept data in POD5 format. We recommend that Nanopore data be submitted in bam or
> fastq format. Your data in bam format should be sufficient for submission to SRA.
> Unfortunately, we don't have a recommendation for an alternative archive to store the raw signal data."
> — The SRA Team, NCBI（2026-09-22）

決定事項:
- **SRA には修飾タグ付き BAM 9 件（6.18 GB）と Illumina FASTQ 18 件（8.13 GB）＝ 計 14.3 GB を提出する。**
  BAM は MM/ML タグを保持しており、SRA 側も「BAM で十分」と明言している。
- **pod5 72.44 GB は SRA には出せない。** NCBI は代替アーカイブの推奨も持っていない。
  → 生シグナルを公開するかどうかは著者判断（`open_items.md` の判断項目を参照）。
  アーカイブしない場合、本文は「modification-tagged BAM が寄託される一次記録である」と明記する必要がある
  （再 basecall は不可能になる）。
- 以前あった「pod5 の返答を待ってから SRA を進める」という保留は**解消**。
  **以下の 1〜5 は pod5 を含まない構成に書き換えてある**（SRA に出すのは BAM と FASTQ だけ）。

## 1. 提出の構造

```
BioProject（1 件）
 └ BioSample × 9（同一培養から DNA/RNA を分取したので 9 個で両アッセイを共有）
    ├ SRA run: ONT BAM（9）
    ├ SRA run: Illumina FASTQ ペア（9）
    （pod5 の run は作らない — SRA は POD5 を受け付けない）
```

## 2. 手順

1. https://submit.ncbi.nlm.nih.gov/ にログイン（NCBI アカウント。ORCID 連携可）
2. **BioProject** → New submission
   - Project type: Raw sequence reads
   - Target: `Streptomyces coelicolor A3(2)`、Sample scope: Multispecies ではなく **Monoisolate**
   - Title / Description は `.zenodo.json` の title/description を流用してよい
   - Release date: **論文 acceptance まで hold**（"Release on specified date" か "upon publication"）
3. **BioSample** → New submission → Package: **Microbe; version 1.0**
   - `biosample_attributes.tsv` をそのままアップロード（9 行、必須項目に空欄なし）
4. **SRA** → New submission → 上記 BioProject / BioSample を紐付け
   - `sra_metadata.tsv` をアップロード（18 行）
   - **pod5 の行は追加しない。** SRA は POD5 を受け付けないため、`pod5_upload_manifest.tsv` は
     SRA 提出には使わない（生シグナルの公開先は著者判断。`著者判断メモ_260922.md` B 項）。
5. **ファイル転送** — Aspera（推奨、`ascp`）または FTP。SRA 画面に出る
   一時ディレクトリのパスと鍵を使う。転送対象:
   - ONT BAM 9 件: `ont_bam_md5.tsv` のパス列（6.18 GB）
   - Illumina FASTQ 18 件: `raw_checksum_report.tsv` のファイル名。実体は
     `~/Library/CloudStorage/Dropbox-SFC-CNS/Takeda Tomoki/my_projects/2026/M145_RNA-seq/data/RH250715199/01_RawData/`（8.13 GB）
6. 転送後、SRA 側の md5 と手元を照合。**2 つのファイルの性質が違うので注意**:
   - `raw_checksum_report.tsv` は Illumina FASTQ 18 本を**受託先の `md5.txt` と照合した結果**
     （2026-09-23、18/18 一致）。独立した基準がある。
   - `ont_bam_md5.tsv` は BAM 9 本から**手元で計算した値**。ONT 側に照合先の manifest が無いため
     **まだ検証されていない**。転送後に SRA 側の md5 と一致することを確認して初めて検証される。
7. 査読者に見せる場合: 提出完了後 SRA から **reviewer link** を発行できる

## 3. 注意点

- **organism は必ず *Streptomyces coelicolor* A3(2)**。受託先報告書の生物種欄
  （*S. albidoflavus* J1074）は転記ミスで、アラインメント率 98.4–98.7% が根拠
- ONT run の `library_strategy` は OTHER / `library_selection` RANDOM / `library_layout` single
- BAM の Reference assembly 欄は `GCF_000203835.1`（NCBI assembly なので fasta 添付不要）
- replicate 番号は元のまま（T1: 1,2,3 / T2: 1,3,4 / T3: 2,3,4）。振り直さない
- **BioProject は 2026-09-23 に取得済み: `PRJNA1533939`**（`SUB16507616`、公開は論文公開時 /
  2027-12-31 のいずれか早い方）。本文 EN/JP の Data availability と Zenodo レコードには反映済み。
  BioSample / SRA の提出時はこの accession を指定する。
- （旧記載の `apply_manuscript_updates.py` は**存在しなかった**。差し替えは 1 文なので直接編集した。
  次に accession が増えたときも同様に EN/JP 対で直接編集し、`make_bootstrap.py` で PASS を確認する。）
  `.zenodo.json` の related_identifiers に BioProject も追加してリリースを更新
