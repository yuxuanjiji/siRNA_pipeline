# 错配数据集（最终版）模型优化建议

> 数据位置：`数据集/最终数据/错配数据集/`（3 个来源，共 **100 行**，2026-09-11 更新）
> 影响的两条模型线：
> ①热力学错配修正：`项目搭建/热力学参数计算/calibrate_mismatch.py` → `mismatch_calibration.json` →（应）`configs/stages.yaml thermo.overrides`
> ②OligoFormer 错配头：`OligoFormer部分/model/mismatch_model.pth` + `data/mismatch*.csv`
> 一句话结论：**这版数据"质量好但规模极小、且类型分布与旧版相反"** —— 优化重点不是加模型容量，
> 而是①补齐标注可解析率 ②按数据集/guide 分组验证 ③用文献先验做收缩 ④把校准真正接进管道。

---

## 0. 本版数据实测画像（已逐行核对）

| 指标 | 实测值 | 说明 / 风险 |
|---|---|---|
| 文件与行数 | Ohnishi_2008 **78** · Holen_2005 **18** · Kini_2009 **4** | 3 个实验体系（Northern / EGFP 荧光 / 原始 Activity） |
| 总计 | **100 行**；WT 对照（mismatch_count=0）**5 行**，突变体 **95 行** | 对照过少 → 配对 Δeff 只能覆盖少数行 |
| 唯一 guide / 重复组 | 71 / **29 组** | 有同 guide 多行（可做组内比较） |
| `mismatch` 标注可解析（`X:Y(pos)`） | **52 / 95**（解析率 55%） | 43 行不可解析：多错配（如 `G:C(13)/G:A(17)`）与 Kini 的 `Guide 1` |
| 标注 vs 序列一致性（可解析行） | 一致 **50** / 不一致 **2** | 2 行需修正 |
| 派生类型分布（可解析 52 行） | **PP 24 / YY 15 / PY 5 / GU 4 / WC 4** | ⚠️ 与旧 674 行集（GU 主导）**正好相反**；4 行"标注错配但实际 WC"必须剔除 |
| 位置分布（可解析） | 3–8、10–15、17–19 有；**1、2、9、16 缺** | g1/g2 无数据 → 无法从本数据验证"g1 宽容" |
| 目标/mRNA 长度 | target 19nt×78、**20nt×3**、**21nt×19**；mRNA 57×78、58×3、**59×19** | ⚠️ 需统一裁到 19nt 窗口 + 57nt 上下文 |
| mRNA 靶窗居中（=mRNA[19:38]） | 96 ✅ / **4 ❌** | 4 行需重定位（与上一条的长度异常同源） |
| 小写字母标记 | **74 行**（如 `...GCUaGUUC`） | 用小写标突变位；必须统一大写并保留该位信息 |
| 归一化效率（0–1） | Holen 0.508±0.202 · Kini 0.594±0.245 · Ohnishi 0.687±0.223（min 0–max 1） | 体系间均值差 ~0.18 → **数据集是协变量，不能直接混池** |

> 结论：**有效可用于"类型×位置带"建模的行 ≈ 50**（52 可解析 − 2 不一致 − 4 WC）。
> 按"自由参数 ≤ 有效样本/15"经验，本数据只支持 **≤3–4 个自由参数**（即位置带级），
> 类型级（4 类）只能靠**文献先验收缩**，不能裸拟合。

---

## 1. P0：先把这 100 行变成"可直接建模"的一版（0.5–1 天）

1. **统一几何**：target 一律 19 nt、mRNA 一律 57 nt，且 `mRNA[19:38] == target`；
   - Holen 的 21 nt/59 nt 行：按 mRNA 中靶窗实际位置裁出 `[19:38]`；
   - Kini 的 `...CUGT` 尾部乱码：去掉多余碱基后再校验。
2. **标注补全（解析率 55% → 100%）**
   - 多错配（`G:C(13)/G:A(17)`）：拆成"位置列表 + 类型列表"，与 `calibrate_mismatch.py` 已有的
     **列计数特征结构**直接兼容（该脚本支持多错配 run/弱闭合项）；
   - Kini 的 `Guide 1`：用 guide↔target 逐列比对自动推断真实 (位置, 类型)，再入表（这 3 行是
     **端部/g1 错配**，正好补 §0 中缺失的 g1 数据）；
   - 小写标记位 → 转成 `mutation_positions` 字段并与 `mismatch` 交叉核对（一致率应 100%）。
3. **剔除/修正 6 行问题样本**：4 行 WC（标注错配但实际完全互补）+ 2 行标注-序列不一致；
   处理动作与原因写进清洗报告（保留可追溯的 `excluded_reason`）。
4. **统一标签口径**：保留 `silencing_efficiency_norm`（0–1）作为主标签；同时保留
   `actual_efficacy + efficiency_unit + dataset`，供"数据集作为随机效应/分层"使用。
5. **产出两个标准文件**（供下游共用，避免各脚本各读一份）：
   - `mismatch_curated.csv`（干净 100 行，统一列名）；
   - `mismatch_curated_with_td.csv`（追加本工程热力学特征：`dG_total / dG_seed / delta_deltaG_ends / MFE_guide / GC`，
     即 `OligoFormer部分/data/mismatch_td.csv` 的列结构）。

## 2. 热力学修正模型（可解释线）优化

1. **参数结构：先带级、后类型级（分层收缩）**
   - 第一层（本数据能支撑）：位置带 4 参数 `g1 / seed(g2-8) / central(g9-12) / 3′(g13-19)`；
   - 第二层：类型惩罚（GU/PP/YY/PY）以 **Turner 2004 / Xia 1998 文献值为先验均值**做 partial pooling，
     样本不足的类自动收缩回先验，而不是拟合出"错配提高效率"的伪系数（旧校准就出现了 β_GU=+0.02、β_PP=+0.09）。
2. **目标函数换成与用途一致的排序损失**
   - 最终用途是排序（rank 阶段），建议 `pairwise ranking / Spearman` 或对 Δeff 的 Huber 回归；
   - 有 WT 对照的 5 行优先做**配对 Δeff**，其余用按 dataset 中心化后的效率。
3. **验证必须分组**：`GroupKFold(by guide)` + **leave-one-dataset-out（LOSO）**；
   旧版 50 行的 LOO 未分组，指标偏乐观，需重算并与新版并列报告。
4. **接进管道（当前最大缺口）**
   - 现状：`thermo_calculator` 跑的是 `DEFAULT_CORRECTION`，`mismatch_calibration.json` 的结论
     （含 s=1.75、zone 方向）**从未被 configs 读取**；
   - 做法：输出 `mismatch_calibration_v2.json`（n / 划分 / CV 指标 / 系数±CI / s / 数据 hash / 版本），
     由 `configs/stages.yaml → thermo.overrides` 引用（模块已支持 `overrides/override_params`），
     加单测断言 override 生效；`mismatch_validated.csv` 时代的结果保留为对照基线。
5. **联动复核**：若 s（seed 放大）或 zone 权重变化，会改变 `feat_dG_seed / feat_dG_duplex`，
   必须重跑真实 SFRP1 全链 + 消融（现有 rank 消融与 `meta.json` 溯源可直接用），比较 top-50 变化。

## 3. OligoFormer 错配头（DL 线）优化

1. **先对齐数据**：`OligoFormer部分/data/mismatch.csv`（52 行，9/10 生成）**未包含 9/11 新增的 Kini**，
   建议用 §1 的 `mismatch_curated*.csv` 重新生成，保证"校准集 = 训练集"同源，避免两边结论打架。
2. **规模现实 → 冻结 + 小参数微调**：100 行（有效 ~50）不能整模微调；只训分类/回归头，
   或 LoRA/Adapter；多种子集成 + 早停；**按 guide 分组 CV**。
3. **混合物理特征（强烈建议）**：把 `*_with_td.csv` 的 dG_total/dG_seed/ΔΔGends/MFE/GC
   与"类型×带 one-hot"显式拼进输入 —— 小样本下比纯端到端稳，也便于与热力学线对照。
4. **输出与损失**：连续效率用回归头（MSE/Huber）+ pairwise 排序辅助损失；
   不确定性用深度集成/MC-dropout；概率输出做温度/等渗校准。
5. **路线对比实验（关键）**：默认 thermo / 带级校准 / 带级+类型收缩 / hybrid DL 四者在
   **GroupCV + LOSO** 下的 Spearman、top-k、分类型误差 → 用数据决定谁进排序（或加权）。
6. **与本项目候选库结合**：我们的突变库是 **g1/g12/g17/g18/g19 单错配变体**，正是错配模型的落点；
   可把校准后的错配惩罚作为 `feat_dG_mismatch_calibrated` 加进 rank（或用"同窗口 WT↔变体"效率差验证增益）。
   ⚠️ 注意本数据集 **g1/g2 位置无样本**（g1 仅 Kini 3 行待解析），g1 的系数仍要靠先验。

## 4. 排期与验收

| 优先级 | 动作 | 产出 | 成本 |
|---|---|---|---|
| **P0** | §1 标准化（几何/标注/剔除/统一标签 + with_td 特征） | `mismatch_curated.csv`、`*_with_td.csv`、清洗报告 | 0.5–1 天 |
| **P1** | 带级校准 v2 + 类型先验收缩；GroupCV/LOSO 报告 | `mismatch_calibration_v2.json` | 1 天 |
| **P2** | 接入 `thermo.overrides` + 单测 + 真实 SFRP1 重跑与消融 | 前后 top-50 对比、`docs/design/mismatch_calibration.md` | 0.5 天 |
| **P3** | hybrid DL vs 热力学四路线对照 | 结论 + 是否进 rank 的建议 | 2–3 天 |
| **P4** | 定向补数据：g1/g2 与位置 9/16、GU 类（本数据集仅 4 行）、每来源补 WT 对照 | 采样矩阵补齐 | 依赖采集 |

**验收清单**
- [ ] 100 行全部：target=19nt、mRNA=57nt 且靶窗在 `[19:38]`；标注 100% 可解析且与序列一致；无 WC 残留；
- [ ] 标签统一 0–1 且保留 dataset/unit 字段；
- [ ] CV 一律分组（guide）+ LOSO；报告 n、划分、指标、系数±CI、数据 hash；
- [ ] 自由参数 ≤ 有效样本/15（有效 ~50 → ≤3–4）；
- [ ] 校准结果写入 configs 并被单测覆盖；全链重跑与旧版结果一并留档。
