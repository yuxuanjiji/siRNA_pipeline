# models/ —— 模型与权重说明（Model Card 摘要；对照《附件5》· 一算法模型 / 三可复现性·模型）

> 本目录不存放权重二进制（体积与许可原因），按《附件5》"模型权重或调用说明"以
> **说明 + 获取方式** 提供。如需打包提交，请按下方路径补充。

## 1. 一句话结论

本项目主管道为**传统（非深度学习）可解释计算方案**：滑窗生成 → 规则过滤 →
热力学/结构特征计算 → 线性加权综合排序；主排序**不自训练**。
**自训练部分**：C_match 变体层错配模型（Ridge 错配代价矩阵 + WT/mismatch Siamese CNN），
使用自建错配数据集（文献来源，见 `data/README.md`）训练，统一入口 `train.py --model cmatch|siamese`
（脚本 `scripts/fit_mismatch_ridge_artifact.py`、`scripts/train_mismatch_siamese.py`），
跨数据集 LOSO 验证 Spearman ρ=0.611 / LOFO 0.594（见 `outputs/analysis/mismatch/` 与
`docs/design/variant_layer.md`）。
OligoFormer 作为**开源预训练效率辅助分**（默认开启，资产缺失时优雅降级），
本项目不训练/微调其权重，仅推理调用。

## 2. 本工程"模型/计算引擎"构成

| 名称 | 类型 | 版本/位置 | 是否自训练 |
|---|---|---|---|
| 排序评分引擎（主分） | 规则/热力学线性加权（无权重文件，逻辑即模型） | `src/sirna_pipeline/stages/rank/ranker.py`；参数见 `configs/stages.yaml → ranking` | 否（可解释公式，版本见 `src/sirna_pipeline/__init__.py::__version__`） |
| **C_match 变体层错配模型** | Ridge 错配代价矩阵 + WT/mismatch Siamese CNN | 训练：`scripts/fit_mismatch_ridge_artifact.py`、`scripts/train_mismatch_siamese.py`（统一入口 `train.py`）；评估产物 `outputs/analysis/mismatch/*.json` | **是**（自建错配数据集；LOSO ρ=0.611 / LOFO 0.594） |
| 热力学 NN 计算器 | 前向最近邻能量模型（只前向，不训练） | `stages/thermo/legacy/thermo_calculator.py` | 否 |
| 结构检测 | ViennaRNA（外部）/纯 Python NN 近似 | `stages/structure/legacy/` | 否 |
| **OligoFormer 效率模型** | 开源预训练 Transformer（silencing-efficiency 分类） | 权重：`项目搭建/OligoFormer部分/model/best_model.pth`（引用副本 `参考模型/OligoFormer-main`，MD5 同源）；仓库 LICENSE 以其自身为准（`项目搭建/OligoFormer部分/LICENSE`） | 否（仅推理；版本以仓库 commit/训练说明为准） |
| 化学修饰建议引擎 | 规则（ESC-19） | `stages/chemmod/rules.py`（RULE_VERSION） | 否 |

## 3. 调用方式与创新贡献（合规声明）

- 调用参数：所有阶段参数在 `configs/`（pipeline.yaml 开关 / stages.yaml 权重阈值）；
  OligoFormer 启用方式见 `docs/design/oligoformer_keep_delete.md` §5。
- 训练声明：本项目**自训练了 C_match 变体层错配模型**（自建错配数据集，见 §2 表格与 `train.py`）；
  不声称自行训练了 OligoFormer（其使用符合其开源许可，未修改其权重）。
  未修改其权重。本项目的**实际创新贡献**为：①把多来源物理特征（结构 MFE/末端 ΔΔG、
  双链 ΔG、seed 结合能）与 DL 效率预测组合成可解释的"主分+辅助分+惩罚"综合排序；
  ②修复 seed 特征退化并给出 seed 结合能的 canonical NN 定义（`docs/design/ranking_weighting.md`）；
  ③六步流水线的工程化封装与大赛标准结果文件输出。

## 4. Model Card（摘要）

- 适用输入：19 nt 引导链（RNA/U）+ 57 nt mRNA 上下文（或 CDS 坐标）；
- 输出：综合排序（final_score/final_rank）与关键指标（热力学、毒性、脱靶）；
- 已知局限：seed 脱靶/毒性依赖外部表与工具链（默认关）；无 ViennaRNA 时结构为近似；
  OligoFormer 效率分仅当其环境就绪后并入；权重为建议默认值（任务16 校准）。

## 5. 打包/获取（如需提交完整附件）

1. 从 `参考模型/OligoFormer-main`（或其工作副本）复制 `model/best_model.pth`
   到本目录（体积 ~? MB，按仓库实际）；若不允许再分发，保留本 README 获取说明即可。
2. 结果文件由 `python predict.py` 生成于 `outputs/results/results.csv`。
