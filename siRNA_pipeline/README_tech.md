# SFRP1 靶向 siRNA 计算筛选与综合排序（siRNA\_pipeline）

第一届全球大学生生命科学挑战赛・赛道二（AI 基因编辑与核酸工具设计）

以 SFRP1 mRNA（NCBI RefSeq NM\_003012，雄激素性脱发相关 Wnt 通路抑制因子）为靶点，

按「候选生成 → 规则筛选 → 结构检测 → 热力学 / 脱靶 / 毒性评估 → 综合排序 → 化学修饰」

六步流水线，从 **14,832** 条候选 siRNA 中计算筛选并输出 **Top-50** 候选清单与逐位修饰建议。

全部为计算验证（可复现）；AI 方法核心为热力学主分 + 深度学习（OligoFormer 效率分）的分层组合排序，

以及自建错配专项模型（C\_match，变体层）。

## 环境要求



* 详细清单（系统 / 硬件 / Python / 依赖 / 外部工具及版本）见 `环境要求.txt`；pip 依赖锁定见 `requirements.txt`

* Python >= 3.9（建议 3.10–3.13，已验证 3.13）；无需 GPU；Windows / Linux 均可

* 安装：`python -m pip install -r requirements.txt`

* 完整管线必需外部组件：ViennaRNA 2.7.2、NCBI BLAST+ 2.14.1 + GENCODE v46、RNA-FM、

  OligoFormer、PITA/TargetScan（版本、安装与配置见 `环境要求.txt`）

## 运行（主入口 predict.py）



```
python predict.py                              # 真实 SFRP1 全链 → 结果

python predict.py --run-name myrun             # 独立运行目录

python predict.py --until structure            # 断点续跑：只跑到某阶段

python predict.py --from thermo --until chemmod

python predict.py --fasta <其他CDS.fa>         # 换靶基因

\# 小规模演示（合成数据，无需真实数据/外部工具）：notebooks/demo\_pipeline.ipynb
```



* 配置：`configs/paths.yaml`（数据路径，相对路径、零硬编码绝对路径）、`pipeline.yaml`（阶段开关、

  随机种子 run.seed=42、赛道标识）、`stages.yaml`（阶段参数：排序权重 /α/β/ 阈值）

* 预期耗时（CPU，真实全链）：无 ViennaRNA 约 60–90 秒；有 ViennaRNA 约 30–35 分钟

## 输出结果（outputs/）



* `results/results.csv` —— 大赛标准结果文件（UTF-8，Top-50；候选编号 / 赛道 / 候选序列 / 关键预测指标 / 模型版本 / 备注）

* `results/rank_final.csv`、`rank_top.csv`、`chemmod_top.csv` —— 全量排序、Top 榜、化学修饰建议

* `runs/<run>/<n>_<stage>/` —— 逐阶段候选表 + `manifest.json`（输入哈希 / 参数 / 统计，可溯源）

* `logs/pipeline_<run>.log`

## 数据说明（详见 data/README.md）



* SFRP1 mRNA（NM\_003012）：NCBI RefSeq，公共域；不在仓库内，经 `configs/paths.yaml` 相对引用

* 错配校准与验证数据：Holen 2005、Ui-Tei 2008、Amarzguioui 2003、Birmingham 2006、HIVsirDB、

  Hu/Mix/Taka/Simone（Sciabola 2013）等；派生表入库 `artifacts/` 与 `outputs/analysis/`，原始文献不入库

* 主管道为可解释计算；**C_match 变体层模型为自训练**（自建错配数据集：Ridge 错配代价 + Siamese CNN，

  训练入口 `python train.py --model cmatch|siamese`，LOSO 跨来源验证 ρ=0.611 / LOFO 0.594）；

  排序权重在公开基准上最小二乘校准（LOSO 跨来源验证），文献正对照回收为独立验证；

  固定随机种子 42，逐阶段 manifest 留痕防泄漏

## 模型与第三方声明（详见 models/README.md・Model Card）



* 主管道为可解释计算（热力学 NN 前向 + 线性加权 + 规则），非自训练

* 自训练模型：C_match 变体层错配模型（错配数据集；`python train.py --model cmatch|siamese`）

* 第三方：OligoFormer（开源预训练 Transformer，仅推理、默认开，资产缺失自动降级）、RNA-FM（序列嵌入）、ViennaRNA、

  BLAST+；版本、调用与许可见 `external/README.md`

* 本项目创新贡献：①热力学 + 深度学习分层组合排序（窗口层 α/β 加权 + 变体层 C\_match 错配代价模型）；

  ②错配从设计禁区转为可标定设计变量；③六步流水线工程化与可复现交付

## 目录结构



```
src/sirna\_pipeline/   核心代码（stages: generation/rule\_filter/structure/thermo/

&#x20;                     toxicity/offtarget/oligoformer/rank/chemmod; pipeline; common）

configs/              配置（paths/pipeline/stages）

scripts/              分析、图表与验证脚本

tests/                单元测试（74 项）

notebooks/            demo\_pipeline.ipynb

data/  models/  external/   数据/模型/外部工具说明

outputs/              结果（results/logs/runs）

docs/design/          设计依据文档
```

## 验证与复现



```
python -m unittest discover -s tests -v                          # 单元测试（无需外部工具/真实数据）

SIRNA\_REGRESSION=1 python -m unittest discover -s tests -v      # 真实数据回归

python train.py --model cmatch                                  # 复现 C_match 错配模型（Ridge 校准）

python train.py --model siamese                                 # 复现变体层 Siamese CNN

python train.py --model oligo\_head                              # OligoFormer 冻结表征 + 小头 LOSO 校准

python scripts/update\_positive\_controls.py                       # 复算文献正对照回收

python scripts/make\_pitch\_figures.py                             # 复现答辩图
```

## 大赛《附件 5》对照



| 要求           | 实现                                                       |
| ------------ | -------------------------------------------------------- |
| 主运行入口        | `predict.py`（README 写明完整命令）                              |
| 训练入口         | `train.py`（C\_match 变体层：错配数据集 Ridge/Siamese 训练）         |
| 环境依赖         | `requirements.txt` + `环境要求.txt`（解释器 / OS/CUDA/ 硬件 / 耗时）  |
| 数据说明         | `data/README.md`（来源 / 许可 / 清洗 / 去重 / 泄漏防控）               |
| 模型说明         | `models/README.md`（Model Card；C\_match 训练声明与创新贡献）          |
| 可执行 Notebook | `notebooks/demo_pipeline.ipynb`                          |
| 最终结果文件       | `outputs/results/results.csv`（UTF-8，一键生成）                |
| 可复现溯源        | 每阶段 `manifest.json` + `run_*_summary.json`（种子 / 参数 / 哈希） |