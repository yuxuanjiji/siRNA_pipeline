# 计算验证结论（task14–16 口径修正版 · 2026-09-11）

> 本文记录对 `scripts/run_experiments.py` 的口径修正（L1）与 A/B 两项实验的结果。
> 所有数字由 `outputs/analysis/task14_benchmark.json / task15_holen.json / task16_ablation.json`
> 直接可查；复现命令见文末。**n 与置信区间必须与点估计一起引用**。

## 1. 为什么做这次修正

原图的"纯热力 vs 纯 DL vs 综合"存在三处口径问题：

1. **热力分是脚本自造的三特征等权 min-max**（`dG_total/dG_seed/ΔΔG_ends`，且一律"越负越好"），
   与管道实际使用的 ranker 权重/方向不一致 → 报告数字 ≠ 项目真实模型；
2. **没有置信区间与显著性**：错配集实际只有 **n=52**，0.18 / 0.245 / 0.285 之间的差异全在噪声内；
3. **NaN 未过滤**：`MFE_guide` 列在现有文件中全为空，NaN 参与排序会产出**伪相关**
   （曾一度误报"MFE ρ=0.43"；真实值见下表）。

## 2. 修正内容（`scripts/run_experiments.py`）

| # | 修正 | 说明 |
|---|---|---|
| 1 | 热力分改用**管道 ranker 口径** | 权重/方向读自 `stages/rank/ranker.py::DEFAULT_RANK_CFG`，并记录 `thermo_method` |
| 2 | 全部指标加 **bootstrap 95% CI + 置换检验 p** | `_stats/_metric` 统一输出 `ci95_low/high`、`p_perm`、唯一值计数 |
| 3 | α 改用 **K 折交叉验证** | `recommended_cv`；全体 argmax 仅作对照并标注过拟合风险 |
| 4 | **自动推导错配位点/位置带** | 最佳窗口偏移搜索（兼容 39/57/67nt 框长）；与 Holen 文献标签 **50/50 一致** |
| 5 | Birmingham 可**自动推导二分类标签** | 由 `沉默效率` 的 log2 双重复取均值（阈值可环境变量覆盖）；有预测列即出 AUC |
| 6 | 新增 `features` 子命令 | 任意 `siRNA+57nt` 表 → 现算 `MFE_guide/dG_total/dG_seed/ΔΔG_ends/GC/位点` |
| 7 | NaN/inf 过滤 | 统计前剔除无效值（这是上述伪相关的根因） |

## 3. 结果 A：task14（完全互补基准，特征为管道现算）

| 数据集 | n | 纯 DL | **纯热力（ranker 口径）** | 热力+DL(0.8/0.2) |
|---|---|---|---|---|
| Hu | 2361 | 0.644 [0.620, 0.666] | **0.132 [0.092, 0.168]** | 0.417 [0.377, 0.449] |
| Taka | 702 | 0.565 [0.506, 0.623] | **−0.192 [−0.246, −0.102]** | 0.096 [0.039, 0.180] |
| Mix | 464 | 0.689 [0.638, 0.735] | **0.150 [0.062, 0.231]** | 0.515 [0.447, 0.574] |

- 相比修正前（等权口径：Hu 0.037 / Taka −0.340 / Mix 0.036），**ranker 口径显著更好**
  （Hu 0.037→0.132、Mix 0.036→0.150、Taka −0.340→−0.192）；
- 但完全互补集上 **DL 仍远强于热力**（0.56–0.69 vs 0.13–0.15），故"热力为主 + DL 为辅"的固定
  权重在这一类候选上不成立；
- **Taka 为异常集**（热力为负）：其 mRNA 带 X 填充（特征计算时 X→A），且效率口径可能与其他集不同，
  建议单独复核或标注为离群基准。

## 4. 结果 B：错配集（n=52）权重与 MFE 权重扫描

### 4.1 逐特征 Spearman（含 CI/p）

| 特征 | ρ | 95% CI | p |
|---|---|---|---|
| **MFE_guide（RNAfold 现算）** | **+0.353** | [0.068, 0.593] | **0.014** |
| S_combo（原 0.7热+0.3奥） | +0.285 | [0.012, 0.515] | 0.044 |
| S_thermo（原三特征等权） | +0.245 | [−0.040, 0.484] | 0.076 |
| 管道 ranker 口径热力分 | +0.181 | [−0.128, 0.441] | 0.186 |
| OligoFormer | +0.180 | [−0.102, 0.451] | 0.199 |
| GC / dG_seed / dG_total / ΔΔG_ends | −0.144 / −0.130 / +0.104 / +0.083 | 均跨 0 | 0.31 / 0.35 / 0.45 / 0.53 |

### 4.2 权重预设对照（特征加权，同一批 52 行）

| 预设 | 权重（MFE / ΔΔG_ends / dG_duplex / dG_seed） | ρ | 95% CI | p |
|---|---|---|---|---|
| **mfe_dominant** | 1.8 / 0.4 / 0.4 / 0.2 | **0.360** | [0.055, 0.626] | **0.0085** |
| mfe_only | 1.0 / — / — / — | 0.353 | [0.068, 0.593] | 0.014 |
| no_seed | 1.4 / 0.6 / 0.6 / — | 0.264 | [−0.052, 0.532] | 0.048 |
| ranker_default（现配置） | 1.0 / 1.0 / 1.2 / 0.8 | 0.181 | [−0.128, 0.441] | 0.186 |
| seed_only | — / — / — / 1.0 | 0.130 | [−0.090, 0.331] | 0.352 |

### 4.3 MFE 权重扫描（a = 非 MFE 热力权重；a=0 即纯 MFE）

CV 均值：a=0.0 → 0.338；0.1 → 0.323；0.2 → 0.314；0.3 → 0.356；**0.4 → 0.374（最优）**；
0.5 → 0.262；0.6 → 0.185；0.7 → 0.165；0.8 → 0.087；0.9 → 0.008；**1.0 → −0.057**。
全样本最优同为 a=0.4（ρ=0.392）。**a=1（去掉 MFE）时相关归零甚至为负** → 错配集上 MFE 是主要信息源。

### 4.4 加 DL 是否更好

以 a=0.4 的 MFE 主导热力分为底，再混入 OligoFormer：w_dl=0.2 → ρ=0.284；w_dl=0.3 → ρ=0.310
（均低于纯 MFE 主导的 0.360）→ **错配候选上 DL 不带来增量**。

### 4.5 位置带分层（位点经文献验证）

| 位置带 | n | 热力 ρ | DL ρ |
|---|---|---|---|
| g1 | 35 | −0.063 | 0.163 |
| seed(2–8) | 8 | **−0.383** | 0.036 |

注：现有 52 行的错配以 g1 为主；seed 带样本仅 8 条，方向异常需扩样复核。

## 5. 结论与建议

1. **错配集预测精度已被提升**：0.285（原综合）→ **0.360**（MFE 主导，p=0.0085）；
   核心是把 **MFE（引导链自折叠）** 从"被 dG_duplex 压住的次要特征"提到主导地位；
2. **建议按候选类型分段权重**（写入 `configs/stages.yaml` 前需再扩样确认）：
   - 完全互补（wt）候选：DL 为主（task14 显示 DL 0.56–0.69）；
   - 含错配（mut）候选：MFE 主导热力分 + DL 权重 ≤0.2（本节 4.2–4.4）；
3. **`dG_seed` 建议降权**：错配集 ρ=−0.13（seed 带内 −0.383），当前权重 0.8 是负贡献；
4. **下一步优先级**：① 扩样 52 → 200+ 后再定权重；② 加入真值结构特征
   （cofold 内环/凸起、可访问性）——MFE 已证明结构方向有效；③ Birmingham 需一个预测列才能算 AUC。

## 6. 复现命令

```bash
# 0) 环境：ViennaRNA CLI（RNAfold）在 PATH 或设 VIENNARNA_BIN
# A) 为完全互补基准现算管道特征
python scripts/run_experiments.py features --input ../数据集/阶段五_计算验证/task14_fullmatch/Hu.csv \
       --out-csv outputs/analysis/features/Hu.csv
#    （Taka / Mix 同理）→ 用真实口径重跑基准
python scripts/run_experiments.py benchmark --input-dir OligoFormer部分 \
       --dataset-dir outputs/analysis/features
# B) 错配集：特征表 + MFE 权重扫描 + 权重预设
python scripts/run_experiments.py ablation --input OligoFormer部分/mismatch_ranked.csv --with-mfe
#    错配集单独算某个预测列（含 MFE 现算、位点分层）
python scripts/run_experiments.py mismatch --input OligoFormer部分/mismatch_ranked.csv \
       --kind holen --prediction pipeline_thermo --with-mfe
# C) 图表（Spearman 为主，含 95% CI 与置换显著性标注）
python scripts/plot_validation_figures.py
```

### 图件清单（`outputs/analysis/figures/`）

| 文件 | 内容 |
|---|---|
| `validation_task14_task15_pair.png/.pdf` | **与原汇报版式一致的 1×2**：左=三基准三口径，右=错配集旧 vs 新 |
| `validation_dashboard.png/.pdf` | 2×2 汇总（上述两图 + 逐特征 + MFE 权重扫描） |
| `validation_task14_spearman.png/.pdf` | 三基准：纯 DL / 纯热力(ranker 口径) / 热力+DL，误差棒=95% CI |
| `validation_task15_scores_spearman.png/.pdf` | 错配集得分对比（旧口径灰、新口径绿、负贡献红） |
| `validation_feature_spearman.png/.pdf` | 逐特征 Spearman（绿色=置换检验 p<0.05） |
| `validation_mfe_weight_sweep.png/.pdf` | MFE 权重扫描：CV 与全样本曲线，标出 CV 最优 a=0.4 |

> 旧图 `task14_benchmark.png` / `task15_mismatch.png` 为**修正前口径**，建议在汇报中替换为上述
> `validation_*` 版本（旧图保留以便对照说明"口径修正前后"）。

## 7. 局限（务必随结论一起引用）

- 错配集 **n=52**，除 MFE 与 mfe_dominant 外差异多不显著；CI 宽度约 ±0.28；
- task14 的 **Simone 集缺预测文件**，未纳入对比；
- Taka 集热力为负相关，疑与 X 填充/效率口径有关，待单独复核；
- Birmingham 的二分类标签已可推导，但缺预测列 → AUC 待补（本项按决定暂不做）。
