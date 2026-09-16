# data/ — 数据说明与获取方式

> 对照《附件5·代码提交要求》：二-目录结构（`data/ # 数据说明与获取方式（不上传受限数据，注明来源与许可）`）
> 与 三-可复现性与溯源要求·数据（来源与许可、清洗与预处理记录、划分方式、去重与数据泄漏防控、**预训练数据披露**）。

## 0. 定位与原则

1. 本目录 = **数据说明与获取方式**。随仓库附带的数据均为**公开文献整理成果或自有整理表格**，不含受限数据；
   大体量/有再分发约束的外部资源**不入库**，在 §2–§4 给出来源、许可与获取命令。
2. 与其他文档分工：运行环境与外部工具安装 → 根目录《环境配置与下载清单.md》与 `external/README.md`；
   模型权重获取与 Model Card → `models/README.md`；BLAST 库构建 → `BLAST/README.txt`。
3. 数据路径一律由 `configs/paths.yaml` 解析（支持环境变量覆盖），无硬编码绝对路径。
4. **受限材料不入库**：出版商论文 PDF、出版商图片、出版商补充材料压缩包一律不随包分发，仅在 `文献来源/来源与许可.md` 记录获取方式；入包内容限于「公开文献数值提取物」「CC-BY 材料」「本项目自有整理成果」三类。
5. **OligoFormer 四个数据集已随包**：`OligoFormer官方数据/`为官方仓库发布版（Hu 2361 / Taka 702 / Mix 472 / Simone 322）。工作区整理副本（如 464 行 Mix）与 `deprecated_Taka.csv`（勘误前版本）**不随包**，以官方版为准；两版差异在 §2 如实记录。
6. 包内每个文件在 `_MANIFEST.csv` 中登记「路径 / 大小 / md5 / 许可」，可逐条校验。

## 0.5 数据文件位置速查（先看这里）

| 你要找的数据 | 本包内位置 | 提交包内？ |
|--------------|-----------|:----------:|
| **错配主表 v1.0**（= 工作区 v3_4，146 行 / 123 设计）+ 逐条校验 + 数据卡 | `data/主表/` | ✅ |
| **SFRP1 CDS**（任务5 输入靶序列） | `data/参考序列/SFRP1-mRNA.txt` | ✅ |
| 管道基准 CSV | **本目录根**：`Holen50.csv`、`Birmingham362.csv`、`positive_controls.csv`（代码契约：`scripts/run_experiments.py` 按根路径直读，勿移动） | ✅ |
| **来源数据**（674 合并总表 / 93 条 / 两份新增批次 / 溯源说明） | `data/来源数据/` | ✅ |
| 文献数值提取物 + 来源登记 + 许可说明 | `data/文献来源/` | ✅（受限材料已剔除，见该目录 `来源与许可.md`） |
| 研究证据（构建脚本 + 可行性报告） | `data/研究证据/` | ✅ |
| **迭代版本**（v1 / v2 / v3 / v3_1 / v3_2 / v3_3、旧版数据卡） | **仅在工作区 `数据集\`**，不进包 | ❌（迭代链见 §6） |
| Hu / Mix / Taka / Simone（OligoFormer 官方数据） | `data/OligoFormer官方数据/`（官方发布版） | ✅ |
| 出版商 PDF / 图片 / 补充材料压缩包 | **不入库**（受限数据） | ❌（获取方式见 `data/文献来源/来源与许可.md`） |

### 本目录结构

```
data/
├── README.md                本文件（数据说明与获取方式）
├── _MANIFEST.csv            包内每个文件的 路径/大小/md5/许可，可逐条校验
├── Holen50.csv / Birmingham362.csv / positive_controls.csv
│                            管道基准（代码契约：根路径直读）
├── OligoFormer官方数据/      Hu / Taka / Mix / Simone（官方发布版）
├── 主表/                    最终版：错配主表 v1.0 + 校验逐条 + 数据卡 + 通俗解读版
├── 参考序列/                 SFRP1-mRNA.txt
├── 来源数据/                674 合并总表（含处理后与处理说明）/ 93 条 / 两份新增批次 / 溯源说明
├── 文献来源/                文献数值提取物 + sources 登记表 + 来源与许可.md
└── 研究证据/                构建脚本（含 assert）+ 可行性报告
```

**合计 58 个文件 / 2.52 MB**（`_MANIFEST.csv` 记录 57 条 路径/大小/md5/许可，可逐条校验）。

**只放最终版本**：`_合并.csv`(v1,120行) 与 `_v2/_v3/_v3_1/_v3_2/_v3_3` 六个迭代代次、以及旧版 `错配siRNA_主表v2_数据卡.md` **均不进包**——它们只存在于工作区 `数据集\` 供溯源；包内的最终版即 `主表/错配siRNA_guide链_主表_v1.0.csv`（= 工作区 `…_合并_v3_4.csv`，md5 一致）。迭代链文字记录见 §6。

> OligoFormer 官方数据已随包（`OligoFormer官方数据/`，Mix 以官方 472 行版为准）；工作区整理副本（464 行 Mix 等）与 `deprecated_Taka.csv`（勘误前）不随包。

## 1. 随仓库附带的数据（不受限，可直接使用）

| 文件 | 行数 | 内容 | 用途（阶段） | 来源 | 许可 | 获取/整理时间 |
|------|-----:|------|--------------|------|------|---------------|
| `Birmingham362.csv` | 362 | Birmingham 2006 脱靶微阵列数据（siRNA/off-target 基因、log2 比值双重复），转录自原文补充表，效率原始数值未改动；另含 2026-09-10 文献重建列（mRNA_57nt 实时切片 / target_gene / refseq / seed_pos，`mRNA_57nt_note` 逐行记录重建方式与失败原因）。本文件与 `文献来源/birmingham2006_offtargets*.csv`（原始 4 列提取与 8 列中间版）为同一批 362 行的三个加工阶段，保留全链供溯源 | 任务15 脱靶基准、消融评估 | Birmingham et al., 2006, Nature Methods 3:213-217, DOI:10.1038/nmeth854（每行 `数据来源` 列含完整引文） | 文献数据，注明出处即可 | 2026-06 收集整理，2026-09-10 文献重建 |
| `Holen50.csv` | 50 | 错配 siRNA 50 条（guide 19nt + mRNA_57nt + 错配位置 + 效率 + `label`） | 任务14 完全互补/错配基准 | Holen et al., 2005, Nucleic Acids Research 33:4704-4710, DOI:10.1093/nar/gki785 (PMC1188085) | 文献数据，注明出处即可 | 2026-06 收集整理 |
| `positive_controls.csv` | 7 | SFRP1 阳性对照 siRNA（sense/guide19/target19/CDS 坐标/57nt 上下文） | 任务17 阳性对照校准 | Suzuki 2008, Br J Cancer 98:1147（si1/si2/si3）；Saini 2009, Cancer Res 69:6815；Wang & Gao 2020, Front Oncol 10:532581；Wang 2024, Regen Ther 28:161；Broadley 2026, IJMS 27:1815（每行 `source` 列含引文） | 文献数据，注明出处即可 | 2026-09 整理 |

`Holen50.csv` 的 `label` 列 = 沉默效率 min-max 归一化到 [0,1]（0.2605–1.0000），与 OligoFormer 官方数据集的 label 约定对齐。

## 2. OligoFormer 预训练与基准数据（Hu / Taka / Mix / Simone）——预训练数据披露

**用途披露**：本项目**未对 OligoFormer 做任何训练或微调**，仅以官方发布权重做推理（效率辅助分）并将其基准集用于对照评估。按《附件5》三，以下完整披露 OligoFormer 官方训练/评测数据的构成（转录自其仓库 README 数据集表，原始文献链接以该表为准）：

| OligoFormer 数据集 | 官方条数 | 细胞系 | 原始文献 |
|--------------------|-----:|--------|----------|
| **Huesken（Hu，训练集）** | 2431 | H1299 | Huesken et al., 2005, Nat Biotechnol 23:995-1001（nature.com/articles/nbt1118） |
| Reynolds | 240 | HEK293 | Reynolds et al., 2004, Nat Biotechnol（nbt936） |
| Vickers | 76 | T24 | Vickers et al., JBC |
| Harborth | 44 | HeLa | Harborth et al., 2003（liebertonline 10.1089/108729003321629638） |
| Ui-Tei | 62 | HeLa | Ui-Tei et al., 2004, NAR 32:936 |
| Khvorova | 14 | HEK293 | Khvorova et al., 2003, Nat Biotechnol |
| Hsieh | 108 | HEK293T | Hsieh et al., 2004, NAR 32:893 |
| Amarzguioui | 46 | Cos-1/HaCaT | Amarzguioui et al., 2003（PMID 12527766） |
| **Takayuki（Taka）** | 702 | HeLa | Katoh & Suzuki, 2007, NAR 35(4):e27（其仓库已勘误 siRNA/passenger 链标注） |

**仓库内副本与行数核对**（统一列约定 `siRNA/mRNA/label/y/td`；`label` ∈ [0,1] 效率，`td` 为热力学特征标记）：

| 文件 | OligoFormer 官方（`OligoFormer部分/data/`） | 整理副本（提交包外 `数据集/最终数据/`） | 差异说明 |
|------|------:|------:|----------|
| `Hu.csv` | 2361 | 2361 | 官方 README 表记 2431 条，发布 CSV 为 2361 行（差 70 行，属上游分发口径，如实披露） |
| `Taka.csv` | 702 | 702 | 与官方表 702 一致 |
| `Mix.csv` | **472** | **464** | ⚠️ 两副本不一致（唯一 (siRNA,mRNA,label) 键：整理版少 8 条、官方版多 16 条）。Mix 为多源合并集，疑对应上表 Reynolds/Vickers/Harborth/Ui-Tei/Khvorova/Hsieh/Amarzguioui 合并去重；**须以官方版为准** |
| `Simone.csv` | 322 | 322 | ⚠️ 出处**待溯源**：OligoFormer README 数据集表未单列该子集，需以其论文（Bai et al., 2024, bioRxiv）正文为准 |

- 获取方式：`git clone https://github.com/lulab/OligoFormer.git`（数据随仓库 `data/` 分发）；许可：其仓库 License and Disclaimer——学术/非商业自由使用。
- 获取时间：2026-09-12（工作副本与引用副本 MD5 校验一致）。
- 训练命令（官方流程，供复核其预训练口径）：`python scripts/main.py --datasets Hu Mix ...`（见 `OligoFormer部分/README.md`）。

## 3. 错配数据资产清单（现有错配数据全景）

**终版（评测口径）**：

| 资产 | 行数/设计数 | 内容 | 位置 |
|------|------------|------|------|
| 错配主表 v1.0（工作区代号 v3_4） | 146 行 / 123 独立设计 | 4 篇文献（Holen 2005 ×50、Ohnishi 2008 ×74、Sierant 2011 ×19、Kini 2009 ×3），guide/mRNA_57nt/错配位置/效率（截断与非截断双归一化列）/来源 | **`data/主表/`**：`错配siRNA_guide链_主表_v1.0.csv` ＋ `错配siRNA_guide链_主表v1.0_校验逐条.csv` ＋ `主表数据卡_v3_4.docx` |

**管道内使用（提交包内）**：

| 资产 | 行数 | 说明 | 位置 |
|------|-----:|------|------|
| `data/Holen50.csv`（根目录） | 50 | 主表 Holen 子集，OligoFormer 列约定（任务14/15 基准；`run_experiments.py` 根路径直读） | `siRNA_pipeline/data/` |
| `OligoFormer部分/data/mismatch.csv` | 52 | Holen 错配 50 + WT 参考 2，OligoFormer 推理输入 | `siRNA_pipeline/OligoFormer部分/data/` |
| `OligoFormer部分/data/mismatch_td.csv` | 52 | 同上 + 热力学特征列（dG_total/dG_seed/ΔΔG_ends/MFE_guide/GC） | 同上 |
| `OligoFormer部分/data/fasta/`、`data/RNAFM/` | — | mismatch 序列 FASTA 与 RNA-FM 嵌入输入 | 同上 |
| `OligoFormer部分/*_predictions.csv` | 4 | OligoFormer 在 Hu/Taka/Mix/Simone 上的推理产物 | `siRNA_pipeline/OligoFormer部分/` |

**来源数据与原始来源材料**（来源数据已随包，见 §0.5 与 `data/来源数据/`；受限的原始材料仍不随包分发）：

| 资产 | 行数 | 说明 | 位置 |
|------|-----:|------|------|
| 合并总表（7 源） | 674 | Birmingham 362 + HIVsirDB 179 + Ohnishi 55 + Holen 50 + Amarzguioui 13 + Schwarz 10 + Ui-Tei 5 | `数据集\错配siRNA数据集_合并总表*.csv`（含处理说明 .md） |
| 完整 5 项精筛 | 93 | Holen 50 + Ohnishi 有效 43（0–1 效率口径） | `数据集\错配siRNA_明确实验沉默效率_93条.csv` |
| 基准表 | 120 | v2 的前身（Holen 50 + Ohnishi 70） | `数据集\错配siRNA_guide链_论文实测效率_合并.csv` |
| 二次收集集 | 100/96 | mismatch.csv / mismatch_curated.csv（Holen 18 + Ohnishi 78 + Kini 4） | `生科挑战赛\数据集\最终数据\错配数据集\` |
| Ohnishi 2008 数值表 | 78（74 错配 + 4 完全互补） | 数值经 desiRm 2011 Table S2–S4 转录 | 同上 `Ohnishi_2008_78_mutated_siRNA_normalized.csv` |
| Holen 18 条解析 | 18 | ⚠️ 该解析 guide 序列错位（0/18 互补核验通过），**已弃用**，仅存档 | 同上 `Holen_2005_18_mutated_siRNA_final_flanked.csv` |
| Kini 2009 原始表 | 4（3 错配 + 1 WT） | guide 列多 1 个插入 G，入库前已修复 | 同上 `Kini_2009_terminal_mismatch_siRNA.csv` |
| WT 参考 | 2 | Holen 完全互补对照 | 同上 `mismatch_wt_references.csv` |
| 原始补充材料 | — | Birmingham MOESM PDF×7、HIVsirDB.zip、desiRm 转录记录、Schwarz FASTA 等 | `数据集\阶段五_计算验证\task15_mismatch\sources\` |
| task14/15/16/17 数据副本 | — | 任务级评估输入与输出 | `数据集\阶段五_计算验证\task14_fullmatch\` 等 |

## 4. 仓库外其他必需输入（`configs/paths.yaml` 引用）

| 名称 | 版本 | 来源与许可 | 获取方式 | 获取时间 | 用途 |
|------|------|-----------|----------|----------|------|
| `SFRP1-mRNA.txt`（SFRP1 CDS，315..1259） | RefSeq NM_003012.5 | NCBI RefSeq（公有领域，无使用限制） | 见 §5 命令 ①，或 NCBI 网页下载后按 `configs/paths.yaml: cds_fasta` 路径放置 | 2026-09-12 | 任务5 滑窗设计输入 |
| GENCODE v46 转录本库（BLAST 近全同源脱靶层） | Release 46 | gencodegenes.org（引用 GENCODE/Ensembl；数据许可要求见 `BLAST/README.txt` 与《环境配置与下载清单.md》§4.2） | 见 §5 命令 ③ | 库构建时 | BLAST 脱靶第三层（可选） |
| ViennaRNA / NCBI BLAST+ / PITA+TargetScan / RNA-FM | 见对应文档 | 开源工具 | `external/README.md`、《环境配置与下载清单.md》§4 | — | 结构检测 / BLAST / 种子区脱靶 / 序列嵌入（均可选，缺省 graceful skip） |

## 5. 获取命令汇总

```bash
# ① SFRP1 CDS（RefSeq NM_003012.5，CDS 315..1259）
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_003012.5&rettype=fasta&retmode=text" -o SFRP1_NM_003012.5.fa
# 截取 CDS 315..1259（或直接使用仓库外已有的 数据集/SFRP1-mRNA.txt）

# ② OligoFormer 官方仓库（Hu/Mix/Taka/Simone 数据 + 权重 + 代码）
git clone https://github.com/lulab/OligoFormer.git

# ③ GENCODE v46 转录本 FASTA + BLAST 库
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
gunzip gencode.v46.transcripts.fa.gz
makeblastdb -in gencode.v46.transcripts.fa -dbtype nucl -out gencode_v46_pc_simple
```

## 6. 数据清洗与预处理记录（溯源）

- **错配主表 v3_4**：处理链 = 7 源合并总表（674 行，2026-06）→ 完整 5 项精筛（93 行）→ 基准合并（120 行）→ 补 Ohnishi 4 条 + 新增 Sierant 2011（19 条）、Kini 2009（3 条）→ v3_4（146 行 / 123 独立设计）。
  逐条清洗规则、负值截断（censored）政策、Kini guide 序列修复（去除插入 G）、9 组重复测量裁定、缺陷清单 D1–D5 见**数据卡** `data/主表/主表数据卡_v3_4.docx`；
  构建脚本（含 assert：123 设计 / 9 / 6 / 5 族 / 142 eligible）在 **`data/研究证据/`**（`build_v3_main_table.py` → `build_v3_4.py`），作为**转换逻辑的说明**随包提供。
  ⚠️ **可复现性边界（如实披露）**：按"只放最终版本"的口径，作为脚本输入的 `_合并_v2.csv` 等中间代次**未随包**，故该脚本链**不能仅凭包内文件重跑**；包内提供的是**冻结的最终表**（`主表/主表_v1.0.csv`）＋**逐行校验表**（146 行逐项核对）＋**数据卡**（清洗规则、负值截断政策、缺陷 D1–D5）。如需包内端到端重跑，须补回 `_合并_v2.csv` 一个文件（约 65 KB）。
- **Holen50.csv** = 主表 Holen 子集（OligoFormer 列约定）；**Birmingham362.csv** = 合并总表 Birmingham 行转录；生成脚本见阶段五 `task14_fullmatch/`、`task15_mismatch/`。
- 所有 RNA 序列统一大写、仅保留 A/U/G/C；效率原始值一律未改动，归一化列另存（截断版 `silencing_efficiency_norm` 与不截断版 `_minmax_unclipped`，后者仅用于保序分析，**不得作回归标签**）。

## 7. 数据划分、去重与泄漏防控

- 本项目**未开展任何模型训练/微调**（主分为可解释计算方法，OligoFormer 仅作可选预训练推理），故无 train/val/test 训练划分义务；评估按任务 14（全互补基准）/ 15（错配基准）/ 16（消融）/ 17（阳性对照）组织，划分即任务边界本身。
- **去重与设计折叠**：错配主表按 `design_id` 折叠为 123 个独立设计；9 组同设计重复测量（Holen 组内）保留并在 `provenance` 列标注（单次测量 114 / 同mRNA重复测量 4 / 重复-待核 28）。
- **泄漏防控**：族安全交叉验证的折单位 = `family_coarse`（5 族：Holen 主座 / Holen g1 报告基因块 / Ohnishi PRNP 体系 / Sierant PSEN1 / Kini EGFP）；`family_locus`（6）与 `family_construct`（9）**不得作折单位**——错位对齐显示三个 PRNP 构建的 mRNA 窗口一致度 100%（重叠 48–56 位），拆入不同折会泄漏 mRNA 上下文。
- **与 OligoFormer 训练数据的重叠**：本项目使用 Hu/Taka/Mix/Simone 仅作推理与对照评估，不参与本方法任何参数拟合；主表错配数据（Holen/Ohnishi/Sierant/Kini）与 OligoFormer 训练集（Huesken）无重叠来源，不存在训练-评测泄漏路径。
- **量纲隔离**：`Birmingham362.csv` 为微阵列 log2 比值（方向与 0–1 效率相反），明确不与沉默效率混池，仅用于脱靶评估口径。
- Ohnishi 4 条"等位单错配"行 `n_mismatch_from_guide=0`，`eligible_guide_strand_mismatch=False`（有效建模样本 142 行 / 123 独立设计）。

## 8. 曾核查、未纳入主表的外部数据（如实披露）

HIVsirDB 数据库、Schwarz 2006（IC50）、Amarzguioui 2003（剩余%）、Ui-Tei 2008（仅图无数值）：已核查，因效率量纲与 0–1 沉默效率不兼容（fold/IC50/log2/无数值）**未纳入**错配基准；原始材料留存于阶段五 `task15_mismatch/sources/`，未随提交包分发。

## 9. 引用清单

错配/校准来源：Holen et al. 2005, NAR 33:4704-4710；Ohnishi et al. 2008, PLoS ONE 3(5):e2248；Ahmed & Raghava 2011 (desiRm), PLoS ONE 6(8):e23443；Sierant et al. 2011, Int J Alzheimer's Dis 2011:809218；Kini & Walton 2009, FEBS J 276:6576-6585；Birmingham et al. 2006, Nature Methods 3:213-217。
阳性对照：Suzuki et al. 2008, Br J Cancer 98:1147；Saini et al. 2009, Cancer Res 69:6815；Wang & Gao 2020, Front Oncol 10:532581；Wang et al. 2024, Regen Ther 28:161；Broadley et al. 2026, IJMS 27:1815。
OligoFormer 训练/评测集原始文献：Huesken et al. 2005, Nat Biotechnol 23:995-1001；Katoh & Suzuki 2007, NAR 35(4):e27；Reynolds et al. 2004, Nat Biotechnol；Hsieh et al. 2004, NAR 32:893；Ui-Tei et al. 2004, NAR 32:936；Khvorova et al. 2003, Nat Biotechnol；Harborth et al. 2003；Vickers et al., JBC；Amarzguioui et al. 2003, NAR。
模型与数据库：OligoFormer（Bai et al. 2024, bioRxiv；github.com/lulab/OligoFormer）；RNA-FM（Chen et al. 2022, arXiv:2204.00300；github.com/ml4bio/RNA-FM）；GENCODE Release 46；HIVsirDB（github.com/raghavagps/hivsirdb）；NCBI RefSeq NM_003012.5。
