# 第五步综合排序 · 加权方法建议与落地口径（ranking_weighting）

> 配套技术路线第五步。本文回答两个问题：
> ① 加权求和具体怎么加、各项权重/方向取多少（给"建议值"+给出方法而非拍脑袋）；
> ② 实现落在哪里、如何校准、如何复核。
> 相关代码：`stages/rank/ranker.py`（默认值/逻辑）、`configs/stages.yaml`（配置）、
> `common/records.py`（列）、`tests/test_rank.py`（回归）。

---

## 1. 终分公式（已落地）

```
final_score = α·score_thermo + β·score_oligo − penalty_total

score_thermo = Σ_i (w_i · norm_i) / Σ_i w_i      # 主分：四热力特征归一化加权，∈[0,1]
score_oligo  = minmax(oligo_efficacy)            # 辅助分：OligoFormer 效率，∈[0,1]
penalty_total= Σ 软惩罚命中扣分                  # 脱靶/毒性 flag，单位同上
final_rank   = 按 final_score 降序
```

- 归一化默认 **min-max + direction 统一为"越大越好"**（z-score 为可选项）；
- **优雅降级**：无 `oligo_efficacy` 列（或部分缺失）时 β=0、α=1，排序与"纯热力版本"一致且可复现；
- 四特征方向：`feat_mfe`=+1（自折叠越弱越好）、`feat_ddg_ends`=+1（末端越不反转越好）、
  `feat_dG_duplex`=−1（双链结合越强越好）、`feat_dG_seed`=−1（seed 结合越强越好）。

## 2. seed 特征口径 v2（本次关键修正）

原实现把 seed 特征定义成"seed 位错配代价"（`dG_mismatch_by_position` 的 g2–g8 求和），
但本库突变位点 **{g1,g12,g17,g18,g19} 全在 seed 之外** → 该特征对所有候选恒 0，min-max 后退化为
常数 0.5，是"死特征"（白占权重位、稀释其余特征）。修正后：

```
feat_dG_seed = Σ_{k∈seed 相邻位} STACK_DG[guide dimer]   # canonical NN 堆叠，负值、随序列变化
             + Σ dG_mismatch_by_position[seed 位]          # 叠加项，当前库=0，为未来 seed 突变留口
```

- canonical 部分只依赖 `guide_checked` + 共享表 `common/nn_tables.py::STACK_DG`，在 rank 特征层
  补充计算，**不改任何 legacy**；
- 该特征现在随窗口 GC/序列内容变化，恢复区分度；
- **方向 −1 需消融复核**：seed 结合强对沉默有利，但也抬高 miRNA 样脱靶——脱靶已由惩罚项对冲；
  若消融显示强 seed 排名前移反而偏离已知高效集，可把方向翻转为 +1（config 一处即改）。

## 3. 加权方法建议（你要的建议）

### 3.1 主分特征权重（默认建议值）

| 特征 | 建议权重 | 理由 |
|---|---|---|
| `feat_mfe`（自折叠 MFE） | 1.0 | 引导链自身结构越少越易装载；与文献弱相关但独立（Singh 2012 正自由能偏好） |
| `feat_ddg_ends`（末端 ΔΔG） | 1.0 | 链选择/装载方向的硬指标，结构阶段已硬筛（<-0.5 淘汰），入排后继续细分 |
| `feat_dG_duplex`（整体双链 ΔG） | **1.2** | 最综合、与沉默效率相关性最强的热力代理，给最高权重 |
| `feat_dG_seed`（seed 结合能） | **0.0（已剔除）** | 实测（四基准网格 + task16）显示该特征对效率预测无正贡献（错配集 ρ=−0.13；带 seed 权重时 mean ρ 更低）→ 权重置 0，仅保留列供诊断 |

> 归一化后各特征量纲一致（0..1），权重即"相对重要度"；四者做加权平均避免分数膨胀。
> **权重依据**：`outputs/analysis/weighting_sweep.json`（四基准 Hu/Taka/Mix/Simone 网格）：
> 当前四特征权重 mean ρ=0.287 → 去掉 `feat_dG_seed` 后配合下表 α/β 提升到 **0.569**。
>
> **鲁棒性核查**（`scripts/weight_robustness_scan.py` → `outputs/analysis/weight_robustness.json`，
> 合格 4218 行）：单权重 ±20% 时 Top-20 Jaccard ≥0.667、Top-50 ≥0.818，全体 ρ≥0.998；
> 200 次随机扰动（各权重 ×U(0.8,1.2)）Top-20 Jaccard 均值 0.818 / 最差 0.667、ρ 均值 0.9995。
> → 排序不依赖权重精调；但 **rank1 在 69.5% 的扰动下会变**，故对外报告 Top-N 而非单一第 1 名。
> β 的影响远大于单权重（β=0.4/0.8 时 Top-20 Jaccard 掉到 0.29/0.38），故 α/β 只按基准集证据定。
> 详见 `docs/design/weight_robustness.md`。

### 3.2 主分 vs 辅助分（α/β）

**已按证据采用 α=0.4、β=0.6**（α+β=1）。依据与修正说明：

- 四基准集（Hu 2361 / Taka 702 / Mix 464 / Simone 322）上：
  **纯 DL** mean ρ=0.550、**旧默认 α0.8/β0.2** mean ρ=**0.287**、**新默认 α0.4/β0.6** mean ρ=**0.569**
  （最差数据集由 0.096 → 0.323）→ 热力主分主导会**稀释** DL 信号，故 DL 权重应高于热力；
- 旧版本曾建议 β∈[0.1,0.3]（"DL 宜作修正项"），该假设**已被基准数据否定**，在此更正；
- 场景区分：错配场景反而以位置模型（Ridge）为主、DL 弱（ρ≈0.2），
  故两场景分开建模；本权重针对**全匹配排序**；
- 模型不可用时仍自动回退 α=1（见 §1），不影响可运行性。

### 3.3 校准方法（比"拍一个数"更可靠，任务16 执行）

1. **建立锚点集**：取一段有实测/文献沉默效率的序列集（如 Huesken 2005、Takayuki 2006 等
   基准效率集，或本实验室 qPCR 验证的少量 Top 候选），得到 (候选, 实测效率) 锚点；
2. **秩相关网格搜索**：在权重超参上做小网格（w 各 ∈ {0.5,0.8,1.0,1.2,1.5}、α∈{0.6..0.9}），
   对每个组合计算预测排名 vs 锚点的 **Spearman ρ / top-K 命中率**，选最优且邻域平滑者；
   样本少（<30）时退化为**敏感性分析**：选 ρ 平台区中点而非峰值（防过拟合）；
3. **稳健性检查**：权重 ±20% 扰动下 top-20 的交并比 ≥0.7 才算稳定；
4. **消融对照**：分别把四特征之一权重置 0、β 置 0、惩罚置 0，观察 top-K 变化，
   把结果写回 `outputs/results/run_*_summary.json` 与本文档版本记录。

### 3.4 惩罚项标定

- 默认：`tox_viability_flag=0.05`、`imm_high_flag=0.10`、`offtarget_risk=0.10`
  （命中即扣，可叠加，只扣分不淘汰）；
- 校准：统计命中率后按"命中候选平均应被推后 5–15 名"来整体缩放罚分；建议保留独立可调键，
  勿与特征权重耦合。若将来有连续分数（PITA 分/TargetScan 分），可升级为
  `罚分 = base · saturating(score)` 的连续惩罚（现阶段 flag 足够）。

## 4. 落地与核查

- 默认值：`ranker.DEFAULT_RANK_CFG`（features/alpha/beta/penalties）；
- 配置覆盖：`configs/stages.yaml → ranking`（`weights/directions/alpha/beta/penalties/aux`），
  ranker 会把 yaml 的 `weights`/`directions` 字典**按特征名叠加**到默认 features 上（消融只需改 yaml）；
- 产出核查：每行新增 `score_thermo / score_oligo / penalty_total / final_score / final_rank`；
  `meta.json` 记录 `features_stats / weights / directions / alpha / beta / aux(available,reason) / penalties`
  与排除原因分布，保证"谁排第几、为什么"可回放；
- 测试：`tests/test_rank.py`（归一化边界=1/0、惩罚、缺特征排除、辅助分存在/降级、seed 特征非退化、排名单调、可复现）。

## 5. 与第三步（结构检测）的关系（口径说明）

技术路线文字称第三步为"强过滤，排除强发夹和链选择反转"。现有实现按更稳的分层语义执行并保持不变：

- **硬淘汰（真"直接筛出"）**：链选择反转（`delta_deltaG_ends < -0.5`）、杂交过弱（cofold MFE > -10）、
  靶标不可及（若校准）；`structure_pass=1` 才能进排序（默认 `require_structure_pass: true`）；
- **软提示（默认不筛）**：自折叠 MFE 偏强、稳定茎——因为"阈值未校准不误杀"，
  它们作为连续特征 `feat_mfe` 进入第五步细分，且 `strict_self_fold` 可一键升级为硬淘汰；
- **无 ViennaRNA 环境**：structure 自动 Level3 近似并标记 `structure_reliable=0`；
  此时建议 `require_structure_pass: false`（回归/演示口径），近似 MFE 仅作提示特征。

一句话：结构筛选的"不合格直接筛掉"落在硬阈值项，MFE/强发夹以软+细分方式参与排序，
两者不重复、不矛盾（语义差见 `docs/design/phase2 §2.2` 与 `stages/structure` 模块头注释）。

## 6. 默认值速查（configs/stages.yaml 与代码一致）

```yaml
ranking:
  alpha: 0.4
  beta: 0.6
  weights:    {feat_mfe: 1.0, feat_ddg_ends: 1.0, feat_dG_duplex: 1.2, feat_dG_seed: 0.0}
  directions: {feat_mfe: 1, feat_ddg_ends: 1, feat_dG_duplex: -1, feat_dG_seed: -1}
  penalties:  {tox_viability_flag: 0.05, imm_high_flag: 0.10, offtarget_risk: 0.10}
  aux:        {enabled: true, source: oligo_efficacy, direction: 1}
```
