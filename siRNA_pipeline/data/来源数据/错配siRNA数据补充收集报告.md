# 错配 siRNA 数据补充收集报告

> 目标：搜集现有已发表的**错配 siRNA**（mismatch siRNA）定量沉默数据，**不重复**当前文件夹已有的 `mismatch.csv`。
> 现有数据集：`数据集/mismatch.csv`（52 条，源自 **Holen et al., 2005**，对应 OligoFormer 用于错配 CNN 训练的数据）。
> 编写日期：2026-06（数据自开放全文/补充材料读取）。

---

## 一、结论：新增 1 个独立可用的数据集

已落盘文件：

| 文件 | 来源 | 行数 | 说明 |
|------|------|------|------|
| `数据集/mismatch_Ohnishi2008_PRNP_allele_specific.csv` | Ohnishi et al. 2008 (PLoS ONE, PMC2373929)；数值经 desiRm (Ahmed & Raghava 2011, PLoS ONE e23443) 补充表 S2–S4 转录 | **58** | 人 PRNP 基因三个致病等位基因（P102L、P105L、D178N）的等位基因特异性错配 si/shRNA |

**已确认完全不重复**：`mismatch.csv`（52 行）中没有任何一条 guide 序列与本数据集重复（去重校验通过，重叠=0）。两者基因（PRNP vs Holen 的靶标）、错配类型、测量系统均不同。

---

## 二、为什么选这一个数据集？

在开放获取（PMC / Europe PMC / NCBI）能直接拿到表格/数值的前提下，可用的“错配 siRNA 定量沉默”数据集主要有以下几类，我逐一核对后认为 **Ohnishi 2008 是唯一能直接、稳定、逐条提取“序列 + 错配位置 + 实测沉默效率”的数据集**：

| 候选来源 | 状态 | 可否直接提取数值 |
|----------|------|------------------|
| **Ohnishi et al. 2008 (PLoS ONE, PMC2373929)** | ✅ OA，全文 + 补充材料 | ✅ 数值在 desiRm 补充表 S2–S4 以文本给出（guide、靶标、错配位置、实测效率） |
| **Du et al. 2005 (NAR 33:1671, PMC1069010)** | ✅ OA 全文 | ⚠️ 57 个单错配靶点数据**以图（条形图）呈现**，表格只给引物；数值无法逐条文本提取 |
| **Huang et al. 2009 (NAR 37:7560, PMC2794185)** | ✅ OA | ⚠️ 240+ 报告基因数据在**图**；补充 .doc 仅图（S4/S5），无逐条数值 |
| **Schwarz et al. 2006 (PLoS Genet, PMC1560399)** | ✅ OA | ⚠️ 数据为时间-裂解曲线（补充 PDF 图 S1），无逐条效率数值表 |
| **Ui-Tei et al. 2008 (NAR 36:2136, PMC2367719)** | ✅ OA | ⚠️ 表 1 仅 IC50（DNA 取代修饰），非“错配基序”逐条表 |
| **Pusch et al. 2003 (NAR 31:6444, PMC275570)** | ⚠️ 非 OA | 数值置于图；无表 |
| **Amarzguioui et al. 2003 (NAR 31:589, PMC140512)** | ⚠️ 非 OA | 数值在图中（残存 mRNA%），无法逐条提取 |

> 结论：**Du 2005、Huang 2009、Schwarz 2006** 都是“位置扫描”型的经典数据，但数值在**图**里，本报告未做图像 OCR/数字化，故**没有**纳入 CSV，避免猜测。后续如需扩展，可读取这些图或用补充表格（Du 2005 的 57 个靶点最值得做）。

因此本次交付聚焦在**文本可直接、无歧义提取**的 Ohnishi 2008 / desiRm 数据集。

---

## 三、数据来源与出处

### 3.1 原始实验
- **Ohnishi Y, Tamura Y, Yoshida M, Tokunaga K, Hohjoh H.**
  *Enhancement of allele discrimination by introduction of nucleotide mismatches into siRNA in allele-specific gene silencing by RNAi.*
  **PLoS One. 2008; 3(5):e2248.** DOI: 10.1371/journal.pone.0002248　PMID: 18493311　PMCID: PMC2373929
  - 在 PRNP 三个致病 SNP（P102L、P105L、D178N）上，用含人工错配的 si/shRNA 做等位基因特异性沉默，双荧光素酶报告基因定量。

### 3.2 数值转录来源（本报告读取的表格）
- **Ahmed F, Raghava GP.** *Designing of Highly Effective Complementary and Mismatch siRNAs for Silencing a Gene.*
  **PLoS One. 2011; 6(8):e23443.** DOI: 10.1371/journal.pone.0023443　PMCID: PMC3154470
  - 其补充材料 **Table S2 / S3 / S4** 将 Ohnishi 2008 的实验数据逐条转录为（siRNA 反义链序列、mRNA 靶序列、错配注解、实测效率、desiRm 预测效率）。
  - 本报告解析了这三个 .docx 表，共提取 **58 条**记录（含 4 条 0 错配正对照）。

---

## 四、字段说明（`mismatch_Ohnishi2008_PRNP_allele_specific.csv`）

| 字段 | 含义 |
|------|------|
| `gene` | 靶基因，恒为 `PRNP` |
| `allele_variant` | 突变等位基因：`P102L` / `P105L` / `D178N` |
| `target_allele` | 本行对应测量的报告基因等位基因：`mutant`（突变等位）或 `wt`（野生型等位） |
| `construct` | 原始表格中的 siRNA 名称（如 `siPrnp105(T10)-5A`） |
| `guide_antisense_19nt` | 反义链（guide）19 nt，5'→3'（已统一为大写） |
| `target_sense_19nt` | 靶链（mRNA sense）19 nt，5'→3' |
| `n_mismatch` | 序列比对得到的错配数 |
| `mismatch_positions_guide5p` | 错配位置（1-based，从 **guide 5' 端**计），多错配用逗号分隔 |
| `mismatch_types_guide:target` | 每个错配的 `位置:guide碱基:靶碱基`，如 `15:U:U`；多错配以 `;` 分隔 |
| `source_mismatch_annotation` | desiRm 表中的原始注解（可能有个别笔误，见下） |
| `silencing_efficiency` | 实测沉默效率（Ohnishi 双荧光素酶报告基因，0–1 范围，个别为负值表示无沉默/实验噪声） |
| `y` | 二值标签：`1` 当 `silencing_efficiency ≥ 0.70`，否则 `0`（与现有 `mismatch.csv` 的阈值约定一致） |

> **坐标约定**：`mismatch_positions_guide5p` 统一用 **guide 的 5'→3' 1-based 位置**（与现有数据集的 siRNA 5' 端计数习惯一致）。`target_sense_19nt` 是反义链按 Watson–Crick 互补对齐出的 sense 靶窗，因此 guide 位置 i 对应靶 19-mer 的（20−i）号碱基。

---

## 五、质量核查

1. **错配数量一致性**：58 条中仅 1 条（`siPrnp178(A9)-5C` 的 wt 靶标行）`n_mismatch` 与 desiRm 表头标注（`1 (wt)` vs 计算 2）不一致——应是 desiRm 原表笔误（`1` 应为 `2`），本报告以其**序列比对结果**为准。
2. **位置/类型一致性**：58 条中 45 条自算错配位置与 desiRm 注解完全一致；其余 13 条为 desiRm 原注解的个别笔误（如单条 off-by-one，或 `U:U` 实为 `U:C`），本报告以**实际序列比对**为准。
3. **去重**：与现有 `数据集/mismatch.csv` 及 `参考模型/OligoFormer-main/data/mismatch*.csv` 比对，guide 序列零重叠。
4. **效率分布**：`n_mismatch` = 0 / 1 / 2 的条数分别为 4 / 28 / 26；`y=1`（效率≥0.7）共 10 条，多为 0 错配或靠近 3' 端、长度不敏感的错配，符合“3' 端错配容忍度较高”的已有规律。

---

## 六、如何进一步扩展（建议）

- **Du et al. 2005 (siCD46 + siNPY)**：57 个单错配靶点的**图**可数字化为逐条效率数据，是最值得补的第二个错配扫描数据集（本报告已保存其全文与补充图于 `_research_src`）。
- **Huang et al. 2009**：240+ 报告基因数据在正文图片；如能读取图，可补充任意位置×任意错配类型的通用容错谱。
- 现有文件夹中 `Hu.csv`（Huesken）、`Taka.csv`、`Simone.csv`、`Mix.csv` 属于**完全互补** siRNA 效率数据，非错配数据，故不在本次“错配”范围内。

---

## 七、复现与本地取材

本次所有原始下载与解析脚本保存在 `_research_src/`（临时目录）：
- 全文/补充 XML：`_research_src/PMC2373929/`、`_research_src/PMC3154470/` 等
- 解析脚本：`parse_desirm2.py`、`build_dataset2.py`、`verify.py`、`qc.py`
- 关键接口：Europe PMC 补充材料 ZIP、NCBI PubMed/PMC 全文、desiRm (PLoS ONE) .docx 表

> ⚠️ 需求提示：`_research_src` 为本次研究临时目录，若需保留原始取证材料请自行移动；`数据集/` 下新增的是正式交付文件。
