# 错配数据集（更新版）模型优化 · 划分严格化报告

> 数据：`数据集/最终数据/错配数据集/`（Ohnishi 78 + Holen 18+1 + Kini 4+1 + wt 参考）
> 代码：`scripts/prepare_mismatch_v2.py`（数据修复 + 划分）、`scripts/optimize_mismatch_model.py`（划分严格化建模）
> 产物：`outputs/analysis/mismatch/mismatch_optimization_{all102,dropkini,ohnishi_lofo}.json`、`mismatch_curated_v2.csv`、`mismatch_splits_v2.json`、`mismatch_qc_v2.json`
> 状态：**research_only —— 不改动生产排序结果**（生产仍是规则/结构/热力学加权管道）

---

## 1. 数据 QC：更新版数据集里发现的四个问题（必须先修）

| # | 问题 | 规模 | 影响面 | 处理 |
|---|---|---|---|---|
| 1 | **Holen/Kini 几何不一致**：原始 21nt 双链被截为 19nt 后 guide 与 target 不再互补；对齐搜索后仍出现 11–13 个"错配"，与论文声明 1–2 个矛盾 | Holen 15 + Kini 5（另 4 行无法对齐） | **只影响序列类特征**（热力学 dG/MFE、RNA-FM/OligoFormer）；位置类特征用论文注释，不受影响 | 在 v2 中标记 `model_eligible_v2=False`；序列类建模只保留几何自洽的 **Ohnishi 78 行** |
| 2 | **注释与序列实算不一致** | 26 / 102 行（Holen 19、Ohnishi 2、Kini 5 中 10 行注释为空） | 位置特征可能错位（如 Ohnishi 出现 13 vs 12 的 1 位偏移） | v2 以"序列实算位点"为准并保留原注释列审计（`repair_note`） |
| 3 | **Kini 注释 `Guide 1` 无法解析** | 5 行 | 位置特征全 0 → 被当成"0 错配/类 WT"，指标不可信 | v2 由序列补齐位点；仍建议人工核对端部错配编号 |
| 4 | **`row_id` 跨数据集重复**（各源文件均从 1 编号） | 全部 | splits 只存 row_id → 跨数据集 train 列表键歧义、join 会串行 | v2 引入全局唯一 `uid = dataset|row_id`，split 全部以 uid 为键，并断言"任何 group 不同时出现在 train/test" |

> 附带修复的环境问题：系统 Python 的 `site-packages/pandas/io/parsers/readers.py` **被覆盖成了项目的 flanking 脚本**（import 时读 fasta），导致 pandas 直接 import 失败；已用 `pip install --force-reinstall pandas==3.0.1` 修复。**请勿直接改写 site-packages 文件**，补丁应留在项目内。

## 2. 划分方案（本次采用）

* **外层**：`LOSO`（Leave-One-Dataset-Out，跨实验体系泛化）为主；另提供 `LOFO`（Leave-One-Family-Out，组内泛化）用于单数据集场景。
* **内层**：`GroupKFold(≤5)` 按 `group_id = dataset|family` 选超参（训练折内完成，测试折绝不参与）。
* **组完整性**：断言 train/test 无共享 group；同 family 的 WT 与突变体、重复测量均同侧。
* **超参网格**：alpha ∈ {0.01,0.1,1,10,100,**300,1000,3000,10000**} —— 原实现的 100 是网格上界，本次多折实际选到 **300/1000/10000**，证明原上限确实受限。
* **指标**：两层拆分下的 **out-of-fold 预测汇总（pooled Spearman）** + 逐数据集 Spearman + **按 group 重采样的 bootstrap 95% CI** + **组内打乱标签的置换检验 p 值**；Kini（n=5）只作提示不作门槛。

## 3. 结果

### 3.1 全 102 行 · 外层 LOSO（复现并强化其报告）
| 模型 | pooled ρ | 95% CI | perm p | Holen | Kini | Ohnishi |
|---|---|---|---|---|---|---|
| **legacy25 + ElasticNet** | **0.370** | [0.175, 0.711] | <0.001 | 0.75 | 0.58 | 0.70 |
| legacy25 + Ridge（其原配置） | 0.327 | [0.045, 0.732] | 0.006 | 0.75 | 0.58 | 0.72 |
| bands_pos + Ridge | 0.294 | [0.016, 0.688] | <0.001 | 0.55 | NA | 0.68 |
| 位置+类型+上下文 + Ridge | 0.290 | [−0.002, 0.710] | 0.024 | 0.86 | NA | 0.70 |
| 位置+类型+上下文+热力学 + Ridge | 0.182 | [−0.105, 0.648] | 0.548 | 0.83 | −0.50 | 0.63 |

*逐数据集数值与他们既有报告一致（0.75/0.58/0.72），说明口径对齐；本次新增 pooled 指标 + CI + 置换检验。*
*注意 Kini 在"无计数特征"的模型里 Spearman 无法计算（位置全 0 → 无方差），这正暴露问题 #3。*

### 3.2 剔除 Kini（97 行）· LOSO 敏感性
| 模型 | pooled ρ | 95% CI | perm p | Holen | Ohnishi |
|---|---|---|---|---|---|
| **legacy25 + ElasticNet** | **0.386** | [0.140, 0.763] | <0.001 | 0.758 | 0.699 |
| legacy25 + Ridge | 0.358 | [0.052, 0.772] | 0.002 | 0.734 | 0.714 |
| bands_pos + Huber | 0.352 | [0.079, 0.721] | 0.002 | 0.852 | 0.673 |

### 3.3 仅几何自洽子集（Ohnishi 78 行）· 外层 LOFO（4 折）
| 模型 | pooled ρ | 95% CI | perm p | 四折范围 |
|---|---|---|---|---|
| **bands_pos + ElasticNet** | **0.600** | [0.566, 0.698] | <0.001 | 0.577–0.727 |
| bands_pos + Huber | 0.584 | [0.545, 0.732] | <0.001 | 0.581–0.704 |
| **bands_only(4 特征) + Ridge** | **0.573** | [0.516, 0.750] | <0.001 | 0.576–0.844 |
| legacy25 + Ridge | 0.569 | [0.490, 0.731] | <0.001 | 0.543–0.760 |
| 位置+类型+上下文 + Ridge | 0.364 | [0.132, 0.730] | 0.004 | 0.60–0.83 |
| 位置+类型+上下文+热力学 + Ridge | 0.269 | [0.178, 0.682] | 0.006 | 0.63–0.84 |

**两个清晰结论**：
1. **特征越多越差**：在 n≤100 下，加入类型/上下文/热力学特征系统性地降低 pooled ρ（0.58 → 0.36 → 0.27）。
2. **最简模型已接近最优**：4 个位置带特征的 Ridge（0.573）与 25 维 Ridge（0.569）持平；ElasticNet 在 pooled 上略优，但 CI 与逐折波动更需留意。

## 4. 建议（按优先级）

1. **研究基线定为**：`bands_only + Ridge`（4 特征、最可解释）或 `legacy25 + ElasticNet`（pooled 略优）。
   报告时**必须同时给出 pooled ρ + 95% CI + 置换 p + 逐数据集 ρ**，并注明 Kini n=5 不设门槛。
2. **暂不使用热力学/类型/上下文特征**（本数据规模下无增益甚至有害）；热力学特征想用，前提是先修好问题 #1 的几何。
3. **数据修复优先级**：Holen（21nt 编号与真实 19nt 窗口对齐）、Kini（`Guide 1` 端部错配编号）→ 修好后才能把热力学/OligoFormer 特征纳入错配建模；在修复前，序列类特征只在 Ohnishi 子集上评估。
4. **划分纪律固化**：以 v2 的 `uid` + `group_id` 为准；外层 LOSO（多体系）/ LOFO（单体系）双报告；任何超参选择只在内层折完成。
5. **生产不变**：以上模型保持 `research_only`；若要进生产排序，需通过晋级门槛（CI 不跨 0、消融不劣化、可解释性审查），并以默认关闭的 config 开关引入。

## 5. 复现命令

```bash
python scripts/prepare_mismatch_v2.py                      # 数据修复 + v2 splits/QC
python scripts/optimize_mismatch_model.py --legacy25 \
    --out outputs/analysis/mismatch/mismatch_optimization_all102.json
python scripts/optimize_mismatch_model.py --legacy25 --drop-kini \
    --out outputs/analysis/mismatch/mismatch_optimization_dropkini.json
python scripts/optimize_mismatch_model.py --legacy25 --datasets Ohnishi_2008 --outer lofo \
    --out outputs/analysis/mismatch/mismatch_optimization_ohnishi_lofo.json
```
