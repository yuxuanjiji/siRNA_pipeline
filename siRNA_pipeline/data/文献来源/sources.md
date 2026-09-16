# 错配 siRNA + 沉默效率数据集 — 来源清单（sources.md）

任务：搜集"存在错配的 siRNA 序列 + 沉默效率/敲降数据"记录。本文件逐来源记录：名称 / URL / 访问状态 / 记录数 / 字段完整度 / 窗口定义与效率口径。
采集日期：2026-09-08。环境：Windows + PowerShell + Python（urllib/curl），网络可直连。

## 主口径说明（适用于本目录全部主表 CSV）
- **siRNA 列**：统一为 guide/antisense 链（与 mRNA 互补的那条）19nt，5'→3'，RNA 字母（U）。
  - Birmingham 2006：原文给的是 sense 链 DNA（19nt）→ 反向互补并 T→U 得 guide；换算已注明。
  - HIVsirDB：原文给的是靶位点 RNA 序列（target RNA，19/21nt，个别 26/27nt）→ 反向互补得 guide。
- **mRNA 列**：任务要求 57nt 窗口（上19+靶19+下19）。所有来源均**未提供**完整 57nt 窗口（Birmingham 只给 off-target 基因号；HIVsirDB 只给靶位点序列）→ 按纪律标注"未获取（…）"，**未自行拼接编造**。
- **沉默效率列**：数值 + 口径，同来源统一（见各来源）。
- **突变位置列**：相对 guide 链 5' 端 1-based；换算规则见各来源；算不出写"未获取"。

---

## 1. siRecords（siRecords 数据库）
- URL：http://sirecords.umn.edu/siRecords/
- 访问状态：**已关停/无法访问**（连接失败 HTTP 000；曾多次重试）。Wayback Machine 可用性接口在本地网络超时（exit 28），未取得存档。
- 记录数：0（未提取）。历史文献记载约 4162 条 siRNA 记录，但当前无法访问，无法核实其 efficacy 标注与错配子集。
- 字段完整度：—（不可访问）。
- 窗口定义与效率口径：—。

## 2. HuSiDa（Human siRNA Database）
- URL：http://www.human-siRNA-database.net/ ；镜像 http://itb.biologie.hu-berlin.de/~nebulus/sirna/
- 访问状态：**已关停/无法访问**（301 跳转后目标不可达 / 404）。
- 记录数：0（未提取）。
- 字段完整度：—。
- 窗口定义与效率口径：—。

## 3. siDirect 2.0
- URL：http://sidirect2.rnai.jp/
- 访问状态：**可访问（HTTP 200）**，但为在线设计/特异性预测工具，无错配+效率数据集下载。
- 记录数：0（无可提取数据）。
- 字段完整度：—（无数据文件）。
- 窗口定义与效率口径：—。

## 4. RNAiAtlas
- URL：http://www.rnaiatlas.ethz.ch/ ；http://rnaiatlas.org/
- 访问状态：**已关停**（连接失败 HTTP 000）。
- 记录数：0（未提取）。
- 字段完整度：—。
- 窗口定义与效率口径：—。

## 5. NCBI GEO — Birmingham 2006（Nat Methods 2006；PMID 16489337；DOI 10.1038/nmeth854）
- URL：https://www.nature.com/articles/nmeth854 ；GEO 检索（esearch gds 多组关键词 + elink pubmed→gds）**未检出该文的微阵列条目**（其敲降数据以论文补充表发布，未入 GEO）。
- 实际提取来源：Nature 论文补充材料 Table 1（MOESM4，PDF）→ 产物 `birmingham2006_offtargets.csv`。
  - Table 1 URL：https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth854/MediaObjects/41592_2006_BFnmeth854_MOESM4_ESM.pdf
  - Table 2（汇总计数，仅计数值）：https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth854/MediaObjects/41592_2006_BFnmeth854_MOESM5_ESM.pdf
  - 附注（MOESM7）：https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth854/MediaObjects/41592_2006_BFnmeth854_MOESM7_ESM.doc
- 访问状态：**可下载**（PDF/DOC 直链，已下载并解析）。
- 记录数：**362 行**（12 条 siRNA 的已验证 off-target 记录；siRNA 维度：C1=27, C14=46, C2=45, C3=12, C4=73, C52=37, G4=6, G41=36, M1=24, M2=11, M3=6, M4=39）。另 50 行排除留档说明见 hivsirdb_mutants_excluded.csv（不含 Birmingham）。
- 字段完整度：siRNA（sense 19nt ✓，已换算 guide）；off-target mRNA 仅给出 RefSeq/Ensembl 转录本号 + 基因名（**无序列、无 57nt 窗口**）；效率 = 双重复微阵列 log2 ratio（✓ 数值）；错配类型 = seed 区匹配型 off-target（论文表注：以达阈值下调为验证标准），**未给单碱基错配位点**。
- 窗口定义与效率口径：mRNA 列 = "未获取（仅off-target基因号，无57nt窗口）"；沉默效率 = "logratio1;logratio2 (微阵列log2比值,双重复)"（生物重复 ×2 或 100nM/50nM 功能重复；负值 = 下调）；突变位置 = "未获取（seed型off-target，论文未给单碱基错配位点）"。

## 6. HIVsirDB（GitHub 仓库 raghavagps/HIVsirDB）
- URL：https://github.com/raghavagps/HIVsirDB ；ZIP：https://codeload.github.com/raghavagps/HIVsirDB/zip/refs/heads/main
- 访问状态：**可下载**（git clone 连接重置，ZIP 下载成功；解压于 `_research_src\HIVsirDB\HIVsirDB-main\hivsir\`）。原始文件：hivsir1.csv（651 条亲本 siRNA）、test_mut_original.csv（92 行突变）、test2_mut.csv、hivsir_esc_seq.csv（107 行逃逸）。
- 记录数（提取）：**突变 72 行**（`hivsirdb_mutants.csv`）+ **逃逸 107 行**（`hivsirdb_escapes.csv`）= **179 行**；另有 **20 行排除留档**（`hivsirdb_mutants_excluded.csv`，reason=verify_fail，见下）。
- 字段完整度（原始 CSV 字段）：id, hiv_strain, ncbi_ac, target_gene, pos_sirna_target, target_site_seq, hiv_esc_seq, len_sirna, efficacy, gc_content, cell_type, sirna_source, transfec_reagent, test_objec, test_method, test_time, pubmed（亲本表）；mut-id, id, mut-seq, mut-eff, mutation(label), number（突变表）；sirna_id, sirna_seq, esc_target, nucleo_len, nt_no, nucleo_status, day, efficacy, pubmed（逃逸表）。
- 窗口定义：靶位点 RNA 序列（19/21nt，逃逸文件个别 26/27nt，含"-#"尾部后缀与 `*`/`-`/`_` 标注）→ **无 57nt 窗口**，mRNA 列 = "未获取（仅靶位点 Nnt：野生型 … → 突变型 …）"。guide = 靶 RNA 反向互补。
- 效率口径：
  - 突变表：fold 敲降（如 "50 fold"，Das 2004 式病毒复制抑制）或百分比数值（保留原样）；亲本效率另存于 detail 表 parent_efficacy 列。
  - 逃逸表：数值（0–99.7，抑制%，"0"=无抑制）或文本 "log"/"High log"（定性：病毒呈对数/高对数逃逸），原样保留。
- 突变位置：guide 5' 1-based；换算 guide_pos = 靶长 L − 靶坐标 pos + 1；逃逸删除记 "del{guide_pos}"（连续删除合并为 "del{start}-{end}"）；星号标注以"两侧均被 `*` 包裹的碱基"为突变位点，与 nt_no 不符时回退为"野生型 vs 突变型比对差异位点"（均已在 detail 表保留原始 esc_target 供复核）。
- **20 条排除记录**（非 bug，为数据库自身不一致，已留档）：
  1) 标签格式不统一（如 '9AU' 无原始碱基标注，'5AU' 与 wt 不符）；
  2) 原始数据 parent_id 关联错误（mut-seq 与所挂亲本 target 完全不同谱，如 vif 序列挂到 Gag 亲本，diff ≥12 处）；
  3) 原始数据本身错误（mut-seq 与 wt 完全相同但标签称有突变；标签与序列对不齐；突变表标签与逃逸表同一记录的星号标注不一致，如 mTt69 的 1263 记录）。

## 7. siRNAEfficacyDB（cellknowledge.com.cn，2025 新库）
- URL：https://cellknowledge.com.cn/siRNAEfficacy ；下载页 https://cellknowledge.com.cn/siRNAEfficacy/download.html
- 数据文件：https://cellknowledge.com.cn/siRNAEfficacy/download/siRNA_all.txt 、sequence-siRNA.fasta 、feature.txt
- 访问状态：**可访问、可下载**（三个文件均已下载至工作区）。
- 记录数：库内含 3544 条效率记录；**按纪律未提取（0 行）**——数据为全匹配（siRNA 与其设计靶位点完全互补）效率数据，**无错配/突变/off-target 字段**（已下载全文检索确认）。
- 字段完整度：—（无错配子集）。
- 窗口定义与效率口径：—。

## 8. OligoFormer（github.com/lulab/OligoFormer）
- URL：https://github.com/lulab/OligoFormer
- 访问状态：可访问；本地已有副本（参考模型\OligoFormer-main，由其他代理核对，本代理仅在线核对 releases/外部数据链接）。
- 记录数：**0（错配训练数据未公开）**。仓库含 scripts/mismatch.py + model/mismatch_model.pth + model/Train.log（2361/473/473 训练/验证/测试对，2023-12 训练，val_rocauc 0.78）——**错配模型权重公开，但训练数据未随仓库或论文（PMC 全文无 mismatch 数据集描述）发布**；仓库内数据文件（data/Hu.csv/Taka.csv/Mix.csv 等）均为全匹配效率数据（不收）；off-target 目录为 TargetScan/PITA 预测脚本与参考序列，无敲降实测。

## 9. siRNADiscovery（github.com/BertramLoong/siRNADiscovery）
- URL：https://github.com/BertramLoong/siRNADiscovery
- 访问状态：可访问；本地副本（siRNADiscovery-2）由其他代理核对，本代理仅在线核对 releases/数据文件。
- 记录数：**0**。仓库 releases 无数据资产；仓库树无 mismatch/off-target 数据文件（仅全匹配效率训练数据）。

## 10. siRNA-Features（github.com/mrichter0/siRNA-Features，治疗性 off-target 预测）
- URL：https://github.com/mrichter0/siRNA-Features
- 访问状态：可访问。
- 记录数：**0（未提取）**。仓库含 gene_alignments3.csv 等 off-target 比对文件，但其 21-mer 比对口径与"guide 链错配 + 敲降幅度"语义结合本任务不可靠（为单一化学修饰 guide 的比对与 RNA-seq 计数，非错配效率记录）→ 留档备查，不进入主表。

## 11. VIRsiRNAdb / ASPsiRNA（crdd.osdd.net）
- URL：http://crdd.osdd.net/servers/virsirnadb/ ；http://crdd.osdd.net/servers/aspsirna/
- 访问状态：**已关停**（HTTP 000）。
- 记录数：0（未提取）。

## 12. Sonnhammer 2008 补充材料（siRNA specificity 论文，Bioinformatics 24:1316–1317）
- URL：https://sonnhammer.sbc.su.se/download/papers/2008_Bioinformatics_24_1316-1317_S1.pdf （Table S1 非特异性 siRNA in siRNAdb；Table S2 基于 Jackson 2003 数据的 off-target 命中）
- 访问状态：可下载（已核对链接存在）。
- 记录数：0（未提取）——为**计算预测的 off-target 命中列表**（特异性打分），无实验敲降幅度与效率数值，不满足"沉默效率"字段要求。

## 13. SeedMatchR（github.com/tacazares/SeedMatchR）
- URL：https://github.com/tacazares/SeedMatchR
- 访问状态：可访问（CRAN + GitHub）。
- 记录数：0（未提取）——R 分析工具（RNA-seq 中检测 seed 介导 off-target），非序列+效率数据集。

## 14. 其他已探明死链/不可用资源（记录状态，不伪造）
- Si_scale：http://www.med.nagoya-u.ac.jp/neurogenetics/Si_scale/ → 404。
- DKFZ siRNA 数据页：http://www.dkfz.de/signaling/ehem-nih/siRNA/ → 404。
- Anderson 2008（RNAi 脱靶，PMC2327361）微阵列原始数据定位于 ArrayExpress **E-MEXP-1402**（仅确认存在，本任务未下载）——可作为后续可选来源，当前记为"已定位未提取"。

---

## 汇总
| 来源 | 访问状态 | 提取行数 | 主表文件 |
|---|---|---|---|
| Birmingham 2006 (Nat Methods) 补充表1 | 可下载 | 362 | birmingham2006_offtargets.csv |
| HIVsirDB mutants | 可下载 | 72（+20 排除留档） | hivsirdb_mutants.csv |
| HIVsirDB escapes | 可下载 | 107 | hivsirdb_escapes.csv |
| siRecords | 已关停 | 0 | — |
| HuSiDa | 已关停 | 0 | — |
| RNAiAtlas | 已关停 | 0 | — |
| siDirect 2.0 | 在线可用，无数据 | 0 | — |
| VIRsiRNAdb / ASPsiRNA | 已关停 | 0 | — |
| siRNAEfficacyDB | 可下载，全匹配 | 0 | — |
| OligoFormer / siRNADiscovery | 可访问，无错配数据 | 0 | — |
| siRNA-Features / SeedMatchR / Sonnhammer 2008 | 可访问，口径不符或无效率 | 0 | — |

**合计：541 行主表记录（362 + 72 + 107）**，另 20 行排除留档（数据库自身不一致）。
