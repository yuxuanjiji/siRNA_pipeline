# 错配数据集更新后的评估与训练集更新（含划分/泄漏核查）

> 日期：2026-09-11 · 状态：research_only（不改生产排序）
> 脚本：`scripts/evaluate_oligoformer_on_mismatch.py`、`scripts/build_mismatch_updated_train_files.py`
> 产物：`outputs/analysis/mismatch/oligoformer_on_updated_mismatch.json`、
> `OligoFormer部分/data/mismatch_updated.csv`、`mismatch_updated_td.csv`

## 1. 旧 52 条训练集 vs 更新后 102 行：泄漏核查

| 检查 | 结果 | 说明 |
|---|---|---|
| 旧集合规模 | 52 行 / **28 条唯一序列** | `data/mismatch.csv`、`mismatch_td.csv`、`mismatch_predictions.csv`、`项目搭建/…/mismatch_validated.csv` 四份**完全同源** |
| **序列级重叠** | **0** | 更新后的 102 行（73 条唯一序列）不含任何旧训练序列 → "删去旧训练集"在序列层面已满足 |
| 家族级重叠 | 0 | 旧 sirna_id 家族与新集无交集 |
| **靶窗级重叠** | **19 行** | 这 19 行与旧训练集**共用同一 19nt 靶窗**（Holen 家族窗口）→ 仍有窗口级泄漏，评估时应排除或与旧行同折 |

> 结论：更新的确是"新数据"；唯一残留泄漏通道是 **19 行同靶窗行**。

## 2. 更新数据上的效果（组安全 CV，OligoFormer 直接推理 + 位置 Ridge）

| 口径 | n | OligoFormer ρ | 位置 Ridge（LOSO/LOFO pooled ρ） |
|---|---|---|---|
| 全部 | 102 | 0.208 | 0.340 |
| 仅 Ohnishi（几何自洽子集） | 78 | 0.208 | **0.569** |
| **剔除与旧训练集同靶窗的 19 行** | 83 | 0.180 | **0.646** |

*位置 Ridge 口径：外层 LOSO（多数据集）/ LOFO（单数据集按 family group），内层 GroupKFold 选 alpha；
OligoFormer 用 `best_model.pth` + RNA-FM 嵌入，与线上同一口径。*

**结论**
1. **OligoFormer 在错配场景依旧弱**：ρ≈0.18–0.21（与旧 n=52 的 0.18 一致）→ 不能当错配预测器；
2. **位置 Ridge 在干净口径下更强**（0.57–0.65），是错配场景的可靠基线；
3. 剔除同靶窗的 19 行后 Ridge 由 0.34 → 0.65，说明这 19 行既是泄漏源也是噪声源；
4. 场景需分开报告：**全匹配排序用 OligoFormer（基准 0.56–0.69）**，**错配用位置 Ridge**。

## 3. 训练集更新（旧文件保留，不覆盖）

| 文件 | 内容 |
|---|---|
| `OligoFormer部分/data/mismatch_updated.csv` | 102 行；列 `siRNA,mRNA,label,y,td`（与旧 `mismatch.csv` 同列，可直接喂 `train_single.py`） |
| `OligoFormer部分/data/mismatch_updated_td.csv` | 追加 `dG_total,dG_seed,delta_deltaG_ends,MFE_guide,GC_content` |
| y 标签 | 沿用旧文件分界（y=0 最大 0.69 / y=1 最小 0.70）→ 阈值 0.695，y1=49 / y0=53 |
| 旧文件 | `data/mismatch.csv` 等**未改动**，需要时可回退 |

用法示例：
```bash
# 训练/评估用更新集（test 集另给，如按数据集留一）
python scripts/train_single.py --path ./data/ --datasets mismatch_updated <test_set> ...
```

## 4. 复现命令（注意解释器）

```bash
# 必须用 OligoFormer 的 venv（fm 依赖 ptflops；系统 python 无）
OligoFormer部分/.venv/Scripts/python.exe scripts/evaluate_oligoformer_on_mismatch.py
OligoFormer部分/.venv/Scripts/python.exe scripts/build_mismatch_updated_train_files.py
```

## 5. 仍需处理
1. 19 行同靶窗行：评估时排除（已实现 `exclude_old_window` 口径）或与旧训练行同折；
2. Holen/Kini 的 21nt 几何与 Kini 端部编号（见 `mismatch_optimization_split_report.md`）；
3. 若要用更新集训练错配模型：用 `mismatch_updated*.csv` + 分组 CV，并**要求 LOSO/LOFO ρ > 位置 Ridge 才允许进生产**。
