# 变体层参数的文献依据（λ_position / 类型因子 / 改善通道）

> 日期：2026-09-12 · 状态：**已据此重定标**（`ranker.py::_VARIANT_LAMBDA_BANDS` 等）
> 来源：`参考文献/` 下 7 个 PDF（6 篇独立文献）→ 文本抽取件 `outputs/_tmp/lit/*.txt`
> 抽取命令：`python -c "import fitz; ..."`（PyMuPDF；PDF 抽取把 5′/3′ 渲染成 `50`/`30`，引文已还原）

## 0. 文献清单与它们各自回答什么

| 文件 | 文献 | 对本模型的作用 |
|---|---|---|
| `gki785.pdf` = `5端.pdf`（同一篇） | **Holen et al., NAR 2005**, "Tolerated wobble mutations in siRNAs decrease specificity, but can enhance activity in vivo" | **本项目 Holen 数据的原文**；逐位点 wobble 耐受度 + "错配可提高活性"的体内证据 |
| `nar_33_15_4704__1.pdf` | 同上文 **Supplementary Table 1**（68 条 duplex 序列） | 序列级核验：21nt、19bp 核心 + 两端 2nt 3′ 突出端、双链 5′→3′ 书写 |
| `file.pdf` | **Ohnishi et al., PLoS ONE 2008**, "Enhancement of Allele Discrimination by Introduction of Nucleotide Mismatches into siRNA…" | **本项目 Ohnishi 数据的原文**；错配位置对"区分度 vs 活性"的分离 |
| `AGO2.pdf` | "High-Throughput Analysis Reveals Rules for Target RNA Binding and Cleavage by AGO2" | **逐位置错配规则**（切割 vs 结合）、transition/transversion 轴、错配**加速**切割的位点 |
| `KRNB_20_2217400.pdf` | "Biochemistry-informed design selects potent siRNAs against SARS-CoV-2" | 设计规则的**位置白名单**与链选择/不对称性规则 |
| `MicroRNA.pdf` | Kertesz et al., "The role of site accessibility in microRNA target recognition" | 可及性打分 ΔΔG = ΔG_duplex − ΔG_open（备用特征） |

## 1. 逐位点结论（决定 λ_position）

| guide 位点 | 文献结论 | 原文引文（≤30 词） | 出处 |
|---|---|---|---|
| **g1** | 在 RISC 中**本就不与靶配对**；guide 5′ 端 C→U 摆动（mRNA 侧 G:C→U:G）**提高活性**；t1A 使亲和力 +1 kcal/mol | *"a wobble mutation (C-to-U, resulting in U:G wobble) in the 50 terminal siRNA:mRNA interaction (w-version) improved activity"*；*"The first nucleotide of the siRNA guide strand does not need to be complementary to the facing target nucleotide: it is anyway unpaired in the RISC complex"* | Holen Fig 5C；KRNB |
| g2–g5 | 错配主要**拖慢结合（association）**，对切割影响小 | *"For both guides, mismatches with seed nucleotides g2–g5 most slowed association rates"* | AGO2 Fig 2B |
| g2–g8 (seed) | 多数种子错配对切割影响很小，**部分反而加速**；Ohnishi：种子区单错配"hardly affects silencing" | *"the majority of seed mismatches had small effects on single turnover cleavage rates … some seed mismatches accelerated cleavage"* | AGO2 Fig 4B；Ohnishi |
| **g9–g11（中央）** | **最脆弱**：三联错配使 k_cleave 下降 >500×；破坏中央配对的单插入几乎不可检出切割；KRNB 要求 9–13（有时 9–14）完全配对 | *"For all possible let-7a triple mismatches and 21 of 27 miR-21 triple mismatches at t9–t11, cleavage was undetectable (>500-fold kcleave decrease…)"*；*"Guide nt 9–13 (and 9–14 for some siRNA sequences) need to be perfectly paired"* | AGO2 Fig 4D；KRNB Fig 1d |
| **g12** | **常可省略、偶尔有益**：t12A 错配使切割速率**快 2.5×**；60 个含 t12A 的双错配中 37 个快于全互补；但 KRNB 把 9–13 一起列为"必须配对" | *"The rate of cleavage for a target bearing t12A mismatched with the let-7a g12G (0.14 s−1) was 2.5-fold faster than the fully complementary t12C target (0.055 s−1)"* | AGO2 Fig 4B |
| g13–g16 | **t13 比 t12 更扰动切割**；Holen 第 16 位错配亦"急剧下降" | *"mismatches at t13, a position not usually considered part of the central region, actually perturbed cleavage more than t12 mismatches"* | AGO2 Fig 4B；Holen Fig 2 |
| **g17–g21（3′ 端）** | 单/双/三错配**提高**单周转切割速率与细胞内敲低；建议"引入 g18–g21 错配可能增强敲低"；设计白名单 {1,15,17,18,19,20,21} | *"single, double, and triple mismatches from t15 to t21 for miR-21 and t16 to t21 for let-7a increased the single-turnover cleavage rate"*；*"mismatches may be harmlessly introduced at positions 1, 15, 17, 18, 19, 20 or 21"* | AGO2 Fig 4B–4D；KRNB |
| 计数上限 | 2 个错配可容忍，**3 个几乎耗尽活性** | *"Triple-mutations seemed to exhaust most of the activity of Fe775i"* | Holen Fig 2/4 |

**位置形状 = 两端低、中央高的"U 形谷"**（Holen 原文：*"wobble mutations in the central part of the antisense strand caused a pronounced decrease in activity, while mutations in the 50 and 30 ends were tolerated very well"*）。

## 2. 错配**类型**的轴（决定 λ_type）

* AGO2 的通用切割模型就是 **(transition/transversion) × 21 个位置 + bulges = 81 参数**，且该轴**在 t12 处符号反转**：
  *"transversions perturbed cleavage more than transitions at positions t6–t11"*；*"at positions t12–t19, transversions were often cleaved at faster rates"*（AGO2 Fig 6D）。
  → 本模型据此把原来的碱基化学分类（PP/YY/PY）换成 **WC / GU / transition / transversion** 轴，并对位置 ≥12 段单独给因子（符号反转）。
* **GU 摆动**：三篇中**没有任何可标定的 GU 数值**（AGO2："未见"；Kertesz 只给"7/8-mer 允许一个 G:U 摆动"的硬规则）。
  → GU 因子只能标注为**纯先验**。
* 中央区"结合 vs 切割"要求相反：*"A cleavage-competent conformation requires pairing at central bases t9–t11"*，但 *"Guide pairing to target bases t9 and t10 has been shown to reduce RISC affinity"*（AGO2）
  → 不应把"结合"和"切割"乘成一个数，应分别输出（与 `qguide_variant_signal.md` P3 一致）。

## 3. "错配可提高活性"的机制与证据（决定改善通道）

| 证据 | 原文引文 | 出处 |
|---|---|---|
| 体内系统性证据（9/9 个 siRNA） | *"In all cases, a wobble mutation (C-to-U, resulting in U:G wobble) in the 50 terminal siRNA:mRNA interaction (w-version) improved activity"*；*"This is the first demonstration of such an enhancing effect in mammalian cells in vivo."* | Holen Fig 5C/5D |
| 机制 = 降低 5′ 端配对能 | *"high-activity siRNA duplexes seem on average to have lower pairing energy in the 50 end of the antisense strands"* | Holen |
| 体外切割速率提升（3′ 端） | *"Mismatches from t17 to t21 resulted in target reduction equal to or greater than that observed for the perfectly complementary (PC) target"* | AGO2 Fig 7D |
| 设计处方 | *"incorporating 3′ mismatches, particularly from g18 to g21, will likely enhance target knockdown"*；*"an siRNA fully complementary to its target is unlikely to be optimal"* | AGO2 |
| 链选择/不对称性 | *"differential base-pairing stability on the two ends of the duplex principally determines the identity of the guide strand"*；*"Fraying the 5′-most nucleotide of the intended guide strand … therefore improves siRNA efficiency"* | KRNB |

**与效率不可分离的一面（必须同时披露）**：Holen 标题即 *"…decrease specificity, but can enhance activity in vivo"*，并指出 *"there is no absolute specificity of siRNAs"*。
同一处 wobble 既提升活性、又降低对完全互补靶标的特异性 → 报告时须与特异性指标并列，不可只报活性。

## 4. 据此对本模型做的修改

| 项 | 改前 | 改后 | 依据 |
|---|---|---|---|
| λ(g1) | 0.0500（先验） | **0.0000** | g1 在 RISC 中不配对；Holen 9/9 增益 |
| λ(seed g2–g8) | 0.1786（拟合 n=24） | 0.1786（不变） | 拟合样本最大，且文献只给"影响小/可加速"的方向 |
| λ(cleavage g9–g11) | 0.0942（拟合 n=4） | **0.3082** | Holen 第 10 位剧烈失活 + AGO2 t9–t11 >500×；原弱拟合与强证据冲突 → 用文献先验 |
| λ(g12) | 0.3082（拟合 n=2） | **0.1500** | AGO2 t12 可省/t12A 反而更快 vs KRNB 9–13 须配对 → 折中并做敏感性 |
| λ(mid g13–g16) | 0.1685（拟合 n=9） | 0.1685（不变） | 拟合可靠；且 λ(g13)>λ(g12) 与 AGO2 一致 |
| λ(p3 g17–g19) | 0.0800（先验） | **0.0000** | AGO2 t15–21 增益 + KRNB 白名单 |
| 类型轴 | PP/YY/PY 化学分类 | **WC/GU/transition/transversion**（≥12 段符号反转） | AGO2 的 81 参数模型即此轴 |
| 改善通道 | 全位点可增益 | **按位点门控** `gain_positions={1,8,12,17,18,19,20,21}` | 仅上述文献支持的位点允许增益 |
| 多错配 | 无 | **≥3 错配额外 ×0.5** | Holen：三重突变几乎耗尽活性 |

**效果（真实 run，`outputs/analysis/variant_layer_verification.json`）**：
超过本窗口 WT 的变体数 **389 → 795**（/3954）；Top-20/Top-50 Jaccard vs 冻结 run = 0.290 / 0.449；
WT 行仍**逐位一致**（0 处差异），公式核对 3954/3954 通过。

## 5. 冲突、不确定性与仍未标定之处

1. **KRNB vs AGO2 对 9–13 的表述不一致**：KRNB 概括为"9–13（有时 9–14）必须完全配对"，AGO2 原文是 t9–t11 + t13；本模型取 AGO2 原文（g12 给中等代价、g13–g16 保留拟合）。
2. **两篇 Holen/Ohnishi 的数值全在图里**（柱/条状图），抽取文本无逐位数值表；Ohnishi 补充表未提供 → **λ 的绝对幅度不可从这两篇标定**，本文只用于**排序与符号**校正。
3. **AGO2 数据体系**：小鼠 AGO2、体外单周转、37 ℃；外推到人细胞 siRNA 有物种与体系差 → 幅度当先验。
4. **改善通道的 confound**：Holen 的第 1 位 wobble 同时削弱了 antisense 5′ 端的**双链**配对（G:C→A:U），因此"增益"可能包含链选择效应；本模型把它统一记入热力项，未单独拆分（消融时应留意）。
5. **增益幅度 γ（`mech_gain_weight`）仍无标定数据**：建议报告 γ = 0 / 0.5 / 1.0 的敏感性；γ=0 等价于"错配只降效率"的旧假设。
6. **仍未落地的文献建议**：把 (Δaffinity, Δk_cleave, C_match) 分别输出（AGO2 结合/切割分离）、可及性特征 ΔG_duplex − ΔG_open（Kertesz）、以及 Ohnishi 的"对变体 1 个 seed 错配近乎免费、对 WT 2 个错配陡增"的**非线性 seed 惩罚**（当前 λ(seed) 是线性的）。

## 6. 复现

```bash
# 文本抽取（PyMuPDF）
python -c "import fitz,pathlib; [ ... ]"     # 见本文件顶部说明；产物在 outputs/_tmp/lit/
python -m unittest tests.test_rank           # λ/类型轴/门控/回退 的单元测试
python scripts/verify_variant_layer.py       # 真实 run 双跑验收
```
