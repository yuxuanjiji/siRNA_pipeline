# 经典错配 siRNA 文献来源登记表（mismatch 数据集扩充）

> 生成日期：2026-09-08
> 适用范围：本文件登记【本次会话负责的经典文献】+ 既有 `数据集\mismatch.csv` 的出处鉴定。
> 既有 `artifacts\sources.md`（另一代理版本，覆盖 gki312 / pgen.0020140 / pone.0002248 / gkp835 / desiRm / OligoFormer / siRNADiscovery 七个来源）仍有效，与本文件互补、不冲突。
> 纪律声明：只收"siRNA 与靶 mRNA 存在错配"的记录；全匹配不收；所有数值均有文本/表格出处，无法确认的一律写"未获取"。

---

## 0. 既有 `数据集\mismatch.csv`（52 行）出处鉴定结论

**结论：全部 52 行数据出自 Holen T, Moe SE, Sørbø JG, Meza TJ, Ottersen OP, Klungland A. "Tolerated wobble mutations in siRNAs decrease specificity, but can enhance activity in vivo." Nucleic Acids Res. 2005;33(15):4704–4710. DOI: 10.1093/nar/gki785；PMC1188085。**

鉴定证据链（五重）：
1. **结构**：表头 `siRNA,mRNA,mismatch_position,y,silencing_efficiency`，52 行。第 1–18 行 = 单一 guide `GCGCUCAUCAUUGUGCUGC`（19nt）在固定 57nt mRNA 窗口（含靶 19nt `AGCAGCACAAUGAUGAGUGC`）上的单/双/三错配位置扫描（效率 0.25–0.8）。第 20–21 行 = guide `CCCCCUACCUGUGGACAUA`（完美匹配，**按纪律排除**）。第 22–53 行 = 10 个 guide × 大鼠/小鼠两套 Aqp4 靶窗口（67nt 原始、成对重复，效率 0.28–0.92）。
2. **block1 命中**：FEN1（人 NM_004111.6 pos523 精确；小鼠 NM_007999.4 同源）。论文正文 FEN1 探针用 NM_004111，siRNA Fe775i + 18 个突变体命名（w3i/w7i/w10i/w3/7i/w7/10i/w3/7/10i/w9i/w16i/w6i/w2i/w19i/w3/7/19i/w2/3/7i/m1Ai/m1Ci/m1Ui + hard3/hard7）与 block1 18 行一一对应。
3. **block2 命中**：guide 与报告基因靶位点恰有 **1 个错配（guide 第 1 位：U:G 或 G:A，即 wobble）**；窗口全长精确命中 Holen 2005 报告基因质粒 **GenBank AY785357（pmuAqp4-Luc，小鼠）** 与 **AY785358（praAqp4-Luc，大鼠）**（验证脚本 `scratch\verify_aqp4_constructs.py`）。论文 Figure 5B 明确：所有 siRNAs 在 antisense 5' 端核苷酸突变（C-to-U）+ 互补核苷酸 (G-to-A)，构成 wt/w/ww/dw 四种退火组合。
4. **row20 归属**：guide `CCCCCUACCUGUGGACAUA` 的 67nt 窗口 BLAST 仅弱命中非模式生物 AQP4（E≈6e-04），在 AY785357/358、rat/mouse Aqp4、FEN1 均无精确命中 → **归属依上下文推断**（该行在序列上与 Holen 报告构建体无完美窗口定位），未直接 BLAST 完美吻合。
5. **排除证据**：gkp835 补充 Data S1 20 条 sense 序列与 mismatch.csv 无一条重合；OligoFormer README 对 mismatch.csv 无来源说明（"仓库自带、源自 Huesken" 为另一代理推测，**已被鉴定为错误**，实际来源为 Holen 2005）。

**注意**：mismatch.csv 的 `silencing_efficiency` 是 0–1 连续值（block1 来自 northern 定量、block2 来自荧光素酶归一化），CSV 内已统一口径；产出 CSV 效率列沿用该数值。第 20–21 行（完美匹配）按纪律剔除，不进入产物 CSV。

---

## 1. Holen 2005（NAR 33:4704，gki785）—— 已产出 50 行

- **标题**：Tolerated wobble mutations in siRNAs decrease specificity, but can enhance activity in vivo
- **DOI**：10.1093/nar/gki785 ｜ **期刊/年份**：Nucleic Acids Research, 2005, 33(15):4704–4710 ｜ **PMC**：1188085
- **靶基因/实验体系**：block1 = 人 FEN1（NM_004111.6，siRNA Fe775i 及其 18 个突变体，HeLa/mouse 细胞 northern+荧光定量）；block2 = 大鼠/小鼠 AQP4 荧光素酶报告基因（pra-Aqp4-Luc=AY785358 / pmu-Aqp4-Luc=AY785357，wobble 双链）
- **窗口定义**：block1 = mismatch.csv 原始 57nt（上19+靶19+下19）；block2 = 从 67nt 原始窗口按 19+19+19 重切为 57nt
- **效率口径**：0–1 连续值（沿用 mismatch.csv，block1 northern 定量 / block2 荧光素酶归一化）
- **突变位置**：block1 沿用 mismatch.csv 的 `mismatch_position`（相对 guide 5' 1-based）；block2 = 1（guide 第 1 位 wobble）
- **提取自**：`数据集\mismatch.csv`（52 行中 50 行，剔除 2 行完美匹配）
- **产物**：`artifacts\holen2005_gki785.csv`（50 行）
- **未获取项及原因**：论文原文补充表 1（40 条 siRNA 全序列）未获（Europe PMC 补充 zip 无该表；OUP 被 Cloudflare 拦截）；block2 的"成对重复"在 CSV 中以两行重复呈现

---

## 2. Amarzguioui 2003（NAR 31:589，gkg147）—— 已产出 13 行

- **标题**：Tolerance for mutations and chemical modifications in a siRNA
- **DOI**：10.1093/nar/gkg147 ｜ **期刊/年份**：Nucleic Acids Research, 2003, 31(2):589–595 ｜ **PMC**：140512
- **靶基因/实验体系**：hTF167i 靶向组织因子 hTF（Acc. M16553，靶位 167–187），HaCaT 细胞，northern（TF/GAPDH 归一化）
- **错配结构**：9 个单突变（s1/s2/s3/s4/s7/s10(=M1)/s11/s13/s16）+ 4 个双突变（ds7/10, ds10/11, ds10/13, ds10/16），均为 GC 颠换，命名按 sense 链 5' 端计数
- **窗口定义**：原文仅给 21nt 靶区（M16553 167–187），mRNA 列即该 21nt（未扩 57nt）
- **效率口径**：剩余 mRNA %（mock=100）。文本给出：wt=10%；s10=20%；s1–s3 与 wt 相当（≈10%）；s4/s7/s11=20–30%；s13/s16 文本仅"severely impaired"无数值；4 个双突变仅给出活性排序无数值 → 无数值者写"未获取"
- **突变位置**：相对 guide 5' 端 1-based（guide pos = 22 − sense pos）：s1=21, s2=20, s3=19, s4=18, s7=15, s10=12, s11=11, s13=9, s16=6
- **提取自**：Europe PMC PDF 全文（`scratch\gkg147b.pdf/txt`）+ 图 1 序列表 OCR（`scratch\gkg147_p2_x13.png`）+ M16553 FASTA（`scratch\M16553.fa`）
- **产物**：`artifacts\amarzguioui2003_gkg147.csv`（13 行）
- **名称↔guide 映射**（CSV 的 siRNA 列为 guide 序列，此处给论文命名）：s1=UAUUUGUAGUGCCUGAAGCGG(21)、s2=UAUUUGUAGUGCCUGAAGCCC(20)、s3=UAUUUGUAGUGCCUGAAGGGC(19)、s4=UAUUUGUAGUGCCUGAACCGC(18)、s7=UAUUUGUAGUGCCUCAAGCGC(15)、s10=UAUUUGUAGUGGCUGAAGCGC(12)、s11=UAUUUGUAGUCCCUGAAGCGC(11)、s13=UAUUUGUACUGCCUGAAGCGC(9)、s16=UAUUUCUAGUGCCUGAAGCGC(6)、ds7/10=UAUUUGUAGUGGCUCAAGCGC(12,15)、ds10/11=UAUUUGUAGUCGCUGAAGCGC(11,12)、ds10/13=UAUUUGUACUGGCUGAAGCGC(9,12)、ds10/16=UAUUUCUAGUGGCUGAAGCGC(6,12)（括号内为突变位置）
- **未获取项及原因**：s13/s16 与 4 个双突变的效率数值（图 2 柱状图无柱顶数值标注）；图 1 OCR 有少量碱基识别噪声（如 s3 多一 C），已按论文"GC 颠换"定义从 wt 构造并经 M16553 靶位核验

---

## 3. Ui-Tei 2008（NAR 36:7100，gkn902）—— 已产出 5 行（效率未获取）

- **标题**：Thermodynamic stability and Watson–Crick base pairing in the seed duplex is a major determinant of the efficiency of the siRNA-based off-target effect
- **DOI**：10.1093/nar/gkn902 ｜ **期刊/年份**：Nucleic Acids Research, 2008, 36(21):7100–7109 ｜ **PMC**：2602766
- **靶基因/实验体系**：siVIM-270（人波形蛋白 NM_003380，靶位点 270 附近，报告基因 psiCHECK-cm 靶序列 VIM-270-cm=CGCCAUCAACACCGAGUUCAAGA），HeLa 双荧光素酶
- **错配结构**：guide 种子区（P2/P3/P5/P6/P7）单点突变 siVIM-270m2/m3/m5/m6/m7（21nt guide，突变为 U→C/G→A/A→G/C→U/U→C），与 WT cm 靶构成单错配；另有 revertant 靶（VIM-270m*-cm）恢复 W:C 配对（未入 CSV）
- **窗口定义**：mRNA 列 = 报告基因 psiCHECK-cm 靶序列 23nt（guide 与第 1–21 位完全配对，RC 验证；22–23 位为旁侧序列）
- **效率口径**：未获取（正文图 6 柱状图，无数字表；正文无可引用数值）
- **突变位置**：guide 5' 端 1-based：m2=2, m3=3, m5=5, m6=6, m7=7（已程序化核验 guide 与 cm 靶第 1–21 位仅在对应位点错配）
- **提取自**：Europe PMC 补充 zip（`scratch\epmc_supp.zip`）→ `gkn902_1.pdf` 补充表 S1/S2（`scratch\gkn902_1.txt`）；效率图 `scratch\gkn902f6.jpg`
- **产物**：`artifacts\uitei2008_gkn902.csv`（5 行）
- **未获取项及原因**：所有效率数值（图 6 柱状图无柱顶标注，OCR 仅读出组标签 P2 U→C/A→G 等，无数值）；sm（seed-match，13nt 非同源）靶系列未入 CSV（非单错配结构）；另一篇 Ui-Tei NAR 36:2136（systematic DNA substitution, OUP article/36/7/2136）未提取（正文为 DNA seed arm 替换研究，定位不同，且补充表未见数字效率表）

---

## 4. Birmingham 2006（Nat Methods 3:199，nmeth854）—— 另一代理已产出

- **标题**：3' UTR seed matches, but not overall identity, are associated with RNAi off-targets
- **DOI**：10.1038/nmeth854 ｜ **期刊/年份**：Nature Methods, 2006, 3(3):199–204 ｜ **PMID**：16489337
- **状态**：本代理不再重复提取。另一代理已产出 `artifacts\birmingham2006_offtargets.csv`（sense_seq、off-target gene、双重复 log2 比值）。性质为 seed 型 off-target（非单一错配位置），按既有产出使用
- **未获取项**：nature.com 原文/补充表被 403/robots 拦截，未重复抓取

---

## 5. Jackson 2003（Nat Biotechnol 21:635，nbt831）—— 无可用错配表

- **标题**：Expression profiling reveals off-target gene regulation by RNAi
- **DOI**：10.1038/nbt831 ｜ **期刊/年份**：Nature Biotechnology, 2003, 21(6):635–637
- **靶基因/实验体系**：16 条 IGF1R + 8 条 MAPK14 siRNA，HeLa 表达谱微阵列；发现仅 11nt 连续同源的 off-target 即被沉默
- **结论**：无成体系"具体错配位置+序列+效率"表。off-target 为部分同源/seed 型（非单一错配位点），补充表 1（siRNA 序列清单）需订阅，未取到
- **提取自**：免费 PDF（`scratch\jackson2003.pdf/txt`，哈佛镜像）
- **产物**：无 CSV（按纪律不产行）
- **未获取项及原因**：补充表 1（Nature 付费）；效率为微阵列表达值而非错配-效率配对

---

## 6. Saxena 2003（JBC 278:44312，M307089200）—— 无数字错配表

- **标题**：Small RNAs with imperfect match to endogenous mRNA repress translation. Implications for off-target activity of small inhibitory RNA in mammalian cells
- **DOI**：10.1074/jbc.M307089200 ｜ **期刊/年份**：Journal of Biological Chemistry, 2003, 278(45):44312–44319
- **靶基因/实验体系**：p21/geminin，3–4 个错配的 siRNA 仍可抑制翻译（luciferase/报告基因）
- **结论**：少样本（错配容忍度方向性结论），数据散在正文图表（翻译抑制 %），无补充数字表、无系统性"错配位置×效率"矩阵
- **产物**：无 CSV
- **未获取项及原因**：无成体系数值表；具体序列+错配位置未在可免费获取的正文/补充中成表

---

## 7. Schwarz 2006（PLoS Genet 2:e140，pgen.0020140）—— 确认 DOI，跳过提取

- **标题**：Designing siRNA that distinguish between genes that differ by a single nucleotide
- **DOI**：10.1371/journal.pgen.0020140 ｜ **期刊/年份**：PLoS Genetics, 2006, 2(9):e140
- **状态**：本地 `_research_src\PMC1560399` 已有（归另一代理），仅确认 DOI，不重复提取

---

## 8. 其他检索过但无可用错配表的经典文献

| 文献 | DOI | 结论 |
|---|---|---|
| Holen 2002（NAR 30:1757，组织因子 tiling） | 10.1093/nar/30.8.1757 | 位置效应/tiling 研究，无错配-效率表 |
| Reynolds 2004（Nat Biotechnol 22:326，理性设计规则） | 10.1038/nbt936 | 设计规则论文，无错配数据集 |
| Anderson 2008（RNA 14:853，seed 互补频率验证） | 10.1261/rna.902708 | seed 频率统计研究，无"错配位置+效率"表 |

---

## 9. 近 10 年/近缘检索（带错配序列+效率候选）—— 均核验为图形数据，0 行

| 文献 | DOI | 结论 |
|---|---|---|
| Wang et al. 2012, PLoS ONE 7(11):e49309「siRNA Has Greatly Elevated Mismatch Tolerance at 3′-UTR Sites」(PMC3493533) | 10.1371/journal.pone.0049309 | 补充表 S1 仅 4 条 guide+完美靶（`scratch\pone49309_s005.txt`）；错配效率在补充图 S1–S4（tif 柱状图，无数值标注）→ 0 行 |
| Sun et al. 2018, NAR 46(13):6806「Differences in silencing of mismatched targets by sliced versus diced siRNAs」(PMC6061797) | 10.1093/nar/gky287 | 补充表 S1 仅 P 值、S4 仅寡核苷酸序列；错配靶沉默效率在正文图（相对荧光素酶）→ 0 行 |
| Aleman et al. 2007, RNA 13:1481「Comparison of siRNA-induced off-target RNA and protein effects」(PMC1800510) | 10.1261/rna.352507 | CXCR4-A 错配在 pos10、CXCR4-B 在 pos11（正文图数据，无数字表）→ 0 行 |
| desiRm（pone.0023443，另一代理已覆盖） | 10.1371/journal.pone.0023443 | 55 行四列完整，见既有 sources.md |

---

## 10. 产物清单（本文件对应）

| 文件 | 行数 | 来源 |
|---|---|---|
| `artifacts\holen2005_gki785.csv` | 50 | Holen 2005 gki785（= 数据集\mismatch.csv 剔除 2 行完美匹配） |
| `artifacts\amarzguioui2003_gkg147.csv` | 13 | Amarzguioui 2003 gkg147 |
| `artifacts\uitei2008_gkn902.csv` | 5 | Ui-Tei 2008 gkn902（效率未获取） |

共 68 行错配记录。所有行均可追溯至上述 DOI/补充表；未获取字段已按纪律标注。
