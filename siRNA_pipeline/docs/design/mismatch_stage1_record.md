# 错配数据第一阶段记录

日期：2026-09-11

本阶段只新增研究用数据准备脚本，不修改生产流水线、排序逻辑或最终冻结结果。

## 固定规则

- 保留三份来源文件的全部记录。
- 按 OligoFormer 输入契约统一为 `siRNA=19 nt`、`mRNA=57 nt`，并使用 `mRNA[19:38]` 作为 19 nt 靶窗。
- 对文献中的 21 nt 序列，明确将第 20-21 位视为 3' 突出端并截断，只保留前 19 nt；对应 mRNA 保留前 57 nt。
- 该截断是来源数据定义的一部分，不是按模型效果选择的裁剪；原始序列仍保存在来源 CSV，标准化动作写入 `normalization_action`。
- 只有标准化后满足 19 nt siRNA、19 nt target、57 nt mRNA，且 target 位于 mRNA `[19:38]` 的记录进入 `model_eligible`。
- 同一 `dataset + sirna family` 作为一个 `group_id`，跨等位基因的 WT 和错配变体不能跨训练/验证/测试。
- 外层使用 leave-one-dataset-out；内层使用 GroupKFold。
- 无法安全标准化的 20/21 nt、58/59 nt 记录不静默裁剪，写入排除原因。
- 仓库中已找到的 Holen `Fe775i` 和 Kini `396` WT 作为独立参考记录加入，来源写入 `reference_source`；不使用推断效率。

## 产出

运行：

```text
python scripts/prepare_mismatch_dataset.py
```

产出目录：`outputs/analysis/mismatch/`

- `mismatch_curated.csv`：全部源记录及建模资格标记
- `mismatch_qc.json`：行数、排除原因、来源统计、输入 hash
- `mismatch_splits.json`：LOSO 和组内 GroupKFold 清单

下一阶段只有在检查这些产物后，才开始实现 Ridge/ElasticNet 与 MismatchModule 基线。

## P1 基线

`evaluate_mismatch_baseline.py` 使用错配数量、位置带和位置 one-hot 特征，执行
LOSO 外层评估，并在训练集内部按 `group_id` 用 GroupKFold 选择正则化强度。
该基线不使用 dataset one-hot，也不使用测试集标签进行标准化。

`evaluate_mismatch_module.py` 只接受存在明确 `mismatch_count=0` WT guide 的靶窗。
没有 WT 的来源不会用 target 或任意变体伪造基准序列，而是在审计中标记为缺失。

`assess_mismatch_baseline_uncertainty.py` 对固定 LOSO 预测执行 95% bootstrap
区间和双侧标签置换检验；这些统计结果只用于模型选择，不改变生产排序。

## 当前统计结论

- 位置特征 Ridge：Holen Spearman=0.752，95% bootstrap CI=[0.406, 0.894]；
	Ohnishi Spearman=0.720，95% CI=[0.604, 0.805]。
- Kini 只有 5 条，位置 Ridge 的 Spearman=0.577，但 95% CI=[-0.395, 1.000]，
	permutation p=0.406，只能作为探索性结果。
- 位置加热力学 Ridge 在 Holen/Ohnishi 的 Spearman 分别为 0.677/0.710，
	未超过位置 Ridge；热力学特征暂不进入正式错配排序模型。
- 当前正式研究基线固定为位置/错配特征 Ridge。`MismatchModule` 和 OligoFormer
	暂不接入生产排序，除非后续新增数据证明其在相同分组评估下稳定超过该基线。