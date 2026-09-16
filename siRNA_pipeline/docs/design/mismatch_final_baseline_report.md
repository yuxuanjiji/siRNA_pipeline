# 错配数据最终基线报告

日期：2026-09-11

## 1. 数据与划分

- 标准化记录：102 条。
- 来源：Holen 19、Kini 5、Ohnishi 78。
- 输入契约：siRNA 19 nt、mRNA 57 nt、靶窗 `mRNA[19:38]`。
- 外层评估：Leave-One-Dataset-Out。
- 内层选择：按 `dataset + siRNA family` 分组的 GroupKFold。
- WT 参考：同靶窗 WT 63 条；同 family 的已核实 WT 39 条。
- 标准化数据 hash：`6884e98e262d9d00d68c8f6355b452bd57659c83339e646af6504b89044f6d3a`。
- family 分组数：7。

## 2. 统一评估结果

| 模型 | Holen Spearman | Kini Spearman | Ohnishi Spearman |
|---|---:|---:|---:|
| 位置特征 Ridge | 0.752 | 0.577 | 0.720 |
| 位置特征 ElasticNet | 0.749 | 0.577 | 0.698 |
| 位置 + 热力学 Ridge | 0.677 | 0.300 | 0.710 |
| 位置 + 热力学 ElasticNet | 0.564 | 0.500 | 0.702 |
| MismatchModule | -0.330* | 0.100* | 0.054* |
| WT 配对 delta Ridge | 0.756 | 0.000* | 0.639 |

`*` MismatchModule 为 3 个随机种子的均值近似，单种子结果波动较大；它不超过位置 Ridge。
`*` Kini 的 WT 配对 delta 在 5 条样本上没有可定义的相关性，原始效率 Spearman 也仅作探索性参考。

WT 配对 delta Ridge 的原始效率误差为：Holen MAE=0.143、Ohnishi MAE=0.145；
这明显优于未配对位置 Ridge 的绝对误差，但它依赖待评分 family 的 WT 基准，因此不能直接作为没有 WT 标签时的全局排序器。

## 3. 不确定性

位置 Ridge 的 bootstrap/置换结果：

| 来源 | Spearman | 95% bootstrap CI | permutation p |
|---|---:|---:|---:|
| Holen | 0.752 | [0.406, 0.894] | 0.0006 |
| Kini | 0.577 | [-0.395, 1.000] | 0.406 |
| Ohnishi | 0.720 | [0.604, 0.805] | 0.0002 |

Kini 仅 5 条，不能作为稳定泛化证据。Holen 和 Ohnishi 支持位置错配特征存在排序信号。

## 4. 模型选择结论

当前正式全局排序研究基线固定为**位置/错配特征 Ridge**。加入当前热力学特征没有稳定提高跨来源排序，
`MismatchModule` 在严格 LOSO 下不稳定，因此二者暂不接入生产排序。

新增的 WT 配对 delta Ridge 作为**相对错配效应模型**保留：它更适合回答“相对同一 WT，某个错配造成多大效率变化”，
不替代没有 WT 标签时的全局排序基线。

当前结果不能证明 Ridge 已经是最终生物学模型，只能证明它是小样本条件下更稳健、可解释的基线。

## 5. 下一步准入条件

只有在新增独立数据后满足以下条件，才重新评估 OligoFormer 或接入 rank：

1. 仍使用相同的 dataset/family 分组规则；
2. 至少两个外部来源上的 Spearman 稳定为正；
3. 多随机种子和 bootstrap 结果不劣于 Ridge；
4. 原有 rank 测试和 ViennaRNA freeze 全部通过；
5. 新模型通过配置开关接入，默认关闭。

## 6. 可追溯产物

- `outputs/analysis/mismatch/mismatch_curated.csv`
- `outputs/analysis/mismatch/mismatch_curated_td.csv`
- `outputs/analysis/mismatch/mismatch_qc.json`
- `outputs/analysis/mismatch/mismatch_splits.json`
- `outputs/analysis/mismatch/baseline_loso.json`
- `outputs/analysis/mismatch/baseline_uncertainty.json`
- `outputs/analysis/mismatch/mismatch_module_loso.json`
- `outputs/analysis/mismatch/mismatch_ridge_artifact.json`：使用全部标准化数据拟合的研究参数，仅供离线研究，不连接生产。
- `outputs/analysis/mismatch/optimized_delta_ridge_loso.json`：使用 WT 配对效应和错配类型的研究评估。
- `outputs/analysis/mismatch/oligoformer_head_loso.json`：冻结 OligoFormer 主干、只训练线性 delta 效率头的首次微调结果。
- `outputs/analysis/mismatch/mismatch_siamese_loso.json`：WT-错配 Siamese CNN 的研究评估结果。

## 7. OligoFormer 首次微调结果

本次实验使用 `best_model.pth`，RNA-FM 预训练表征，102 条标准化记录，
4888 维 merge 表征；全部 OligoFormer encoder 冻结，只训练线性回归头，
使用 WT 配对 `delta_label`，三随机种子和 LOSO。

结果不超过 Ridge：

| 来源 | OligoFormer 头部 raw Spearman（3 seeds） | 位置 Ridge |
|---|---:|---:|
| Holen | -0.176 / -0.341 / -0.244 | 0.752 |
| Kini | -0.300 / -0.600 / -0.300 | 0.577 |
| Ohnishi | 0.205 / 0.152 / 0.047 | 0.720 |

因此当前不解冻 encoder、不接入生产排序。该结果说明在现有样本规模和来源异质性下，
预训练 OligoFormer 表征加线性头没有提供超过位置 Ridge 的增益；后续应优先补充独立数据，
而不是增加微调自由度。

## 8. Siamese CNN 首次结果

新增 WT-错配 Siamese CNN，使用共享 1D-CNN、`mutant-WT` 与绝对差值、错配位置/类型、
热力学特征，损失为 Huber + 0.25 pairwise hinge。LOSO 结果为三随机种子均值：

| 来源 | raw Spearman | raw MAE | 位置 Ridge raw Spearman |
|---|---:|---:|---:|
| Holen | 0.577 | 0.166 | 0.752 |
| Kini | -0.100 | 0.303 | 0.577 |
| Ohnishi | 0.595 | 0.155 | 0.720 |

Siamese CNN 明显优于旧版 `MismatchModule`，但目前仍未超过位置 Ridge，尤其 Kini 只有 5 条。
因此它作为后续可继续研究的候选模型保存，暂不替代 Ridge，也不接入生产排序。