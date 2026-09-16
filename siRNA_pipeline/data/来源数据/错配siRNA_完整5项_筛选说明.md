# 错配siRNA「5项齐备」筛选说明

> 目标：读取 `数据集\` 文件夹内所有数据集，筛出**同时具备**【siRNA、mRNA 序列、可归一化沉默效率、文献地址、突变位置】的记录，整理为单个文件。
> 输出：`数据集\错配siRNA_完整5项_可归一化.csv`（**50 行**）
> 筛选口径：**mRNA 须为完整 57nt（19 上游 + 19 靶位点 + 19 下游）**；可归一化效率指已是 0–1 沉默效率并能 min-max 到 [0,1] 的值。

---

## 一、扫描的数据集（文件夹内全部数据文件）

| 文件 | 行数 | siRNA | mRNA序列 | 可归一化效率 | 文献地址 | 突变位置 | 是否入选 |
|------|------|:----:|:--------:|:----------:|:--------:|:--------:|----------|
| `错配siRNA数据集_合并总表_处理后.csv` | 674 | ✅ | ✅(50行有57nt) | ✅ | ✅ | ✅ | ✅ **数据来源** |
| `mismatch.csv` (Holen) | 52 | ✅ | ✅ | ✅ | ❌ | ✅ | ⚠️ 并入上面 |
| `mismatch_Ohnishi2008_PRNP_allele_specific.csv` | 58 | ⚠️(guide_antisense) | ⚠️(仅19nt) | ✅ | ❌ | ✅(mismatch_positions) | ❌ 列名非标准、无文献地址列、无57nt |
| `Hu.csv` | 2361 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ 非错配数据集、无突变位置、无文献列 |
| `Taka.csv` | 702 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ 同上 |
| `Mix.csv` | 464 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ 同上 |
| `Simone.csv` | 322 | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ 同上 |
| `mismatch.csv` 与 `错配...合并总表_处理后.csv` | — | — | — | — | — | — | 同一数据源(Holen)的不同阶段文件，以处理后表为准 |
| （非数据文件）`*.md` / `*.json` / `*.txt` / `*.docx` | — | — | — | — | — | — | 非记录级数据集，不参与 |

### 关键判断
- **Hu/Taka/Mix/Simone**：它们是**完全互补** siRNA 效率数据集（Silencing 效率回归用），**不含"突变位置"列**，也不带文献地址；`label` 是 0–1 效率，但**缺突变位置 + 文献地址** → **不符合本筛选条件**。
- **mismatch.csv**：Holen 2005 的有效子集，但与"处理后表"里 Holen 行是同一批，且无独立文献地址列 → 直接以**处理后表**为准。
- **mismatch_Ohnishi2008...**：列名非标准（`guide_antisense_19nt` / `target_sense_19nt`），仅 19nt 靶位点无 57nt 全长，且无独立文献地址列 → 不入选本文件（其可归一化效率在"处理后表"里另有体现，但无 57nt mRNA，按你选的**严格口径**排除）。

---

## 二、最终入选：50 行（全部来自 Holen et al., 2005）

在 `错配siRNA数据集_合并总表_处理后.csv`（674 行）中，按以下 5 条件筛选得到 **50 行**：

1. `siRNA` 为有效核酸序列（≥17nt，仅 A/U/G/C）
2. `mRNA_57nt` 已构建且为完整 **57nt**（19+19+19，仅 A/U/G/C）
3. `silencing_efficiency_norm` 非空（已 min-max 归一化到 [0,1]）
4. `数据来源`（文献地址）非空
5. `突变位置` 非空且不为"未获取"

**结果统计**：
- 行数：**50**
- 来源：全部 `Holen et al., 2005, Nucleic Acids Research 33:4704-4710, DOI:10.1093/nar/gki785 (PMC1188085)`
- `silencing_efficiency_norm` 范围：**0.2605 ~ 1.0000**
- 所有 `mRNA_57nt` 均长度 57、仅 ACGU
- 所有 5 项字段均非空

---

## 三、输出文件列说明

| 列 | 说明 |
|----|------|
| `siRNA` | 反义链(guide)，5'→3'，19nt |
| `mRNA_57nt` | 57nt mRNA 靶序列：19 上游 \| 19 靶位点 \| 19 下游 |
| `mismatch_position` | 错配位置（1-based，计于 siRNA 5' 端，多错配逗号分隔） |
| `silencing_efficiency` | 原始沉默效率（0–1，未改动） |
| `silencing_efficiency_norm` | min-max 归一化到 [0,1]，4 位小数（`x'=(x−0.0140)/(0.9200−0.0140)`） |
| `数据来源` | 文献地址（含 DOI / PMCID） |

---

## 四、回读验证

| 验证项 | 结果 |
|--------|------|
| 行数 | 50（符合） |
| 所有 `mRNA_57nt` 长度==57 且仅 ACGU | ✅ |
| 所有 `silencing_efficiency_norm` ∈ [0,1] | ✅ |
| `mismatch_position` / `数据来源` / `siRNA` 均非空 | ✅ |
| 原始 `silencing_efficiency` 未改动 | ✅ |

**抽样核对**（5 行）：`UGAGCUGGAUUCAUGCUGG`(mut 1, eff 0.62→norm 0.6689)、`UCCAUGGCCAGCAGUGAGG`(mut 1, eff 0.28→norm 0.2936)、`GCGCUCGUCAUUGUGCUGU`(mut 3,7,19, eff 0.70→norm 0.7572)、`GCACUCGUCAUUGUGCUGC`(mut 7, eff 0.78→norm 0.8455) —— 均在 0≤norm≤1 且 mRNA_57nt 结构正确。

---

## 五、说明与局限

1. 本文件**未新增任何合成序列**；`mRNA_57nt` 全部取自源数据已提供的完整 57nt 窗口。
2. 归一化仅对"已是 0–1 且方向一致"的样本（Holen 50 行）用 `min=0.0140, max=0.9200` 做 min-max；**不包含**量纲不兼容的来源（Ohnishi 负值、Birmingham log2、HIVsirDB fold、Schwarz IC50、Amarzguioui 剩余%、Ui-Tei 无数值）——如需更全可另行按"有真实 mRNA 序列即收（含 19nt 靶位点）"的宽松口径扩到 93 行。
