# SFRP1 靶向 siRNA 计算筛选与优化

> 基于功能分区差异化互补（Functional-Zone Differential Complementarity）的 SFRP1 靶向 siRNA 全流程计算筛选管线
> 第一届全球大学生生命科学挑战赛・赛道二（AI 基因编辑与核酸工具设计）



![Python](https://img.shields.io/badge/Python-3.9%2B-blue)



![License](https://img.shields.io/badge/License-MIT-green)



![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-lightgrey)



![No GPU required](https://img.shields.io/badge/GPU-Not%20Required-orange)

一个面向雄激素性脱发（Androgenetic Alopecia, AGA）的 siRNA 计算药物设计项目：以 Wnt 通路抑制因子 **SFRP1**（NCBI RefSeq NM\_003012）为靶点，将「非核心区错配」从设计禁区转化为可控设计变量，通过 **生成 → 筛选 → 排序 → 修饰** 的全计算管线，从 14,832 条候选 siRNA 中筛选并输出高效、低毒、可递送的 Top-50 候选清单。全流程可一键复现、可溯源。



***

## 目录



* [背景与科学问题](#背景与科学问题)

* [核心创新](#核心创新)

* [技术路线](#技术路线)

* [验证与结果](#验证与结果)

* [快速开始](#快速开始)

* [目录结构](#目录结构)

* [技术栈与依赖](#技术栈与依赖)

* [输出结果](#输出结果)

* [模型与第三方声明](#模型与第三方声明)

* [可复现性](#可复现性)

* [License](#license)

* [致谢与参考数据](#致谢与参考数据)

## 背景与科学问题

**疾病背景**：雄激素性脱发是发病率最高的进行性脱发类型。SFRP1 是 Wnt/β-catenin 通路的分泌型抑制因子，在脱发区毛乳头细胞中高表达，抑制毛囊干细胞活化；敲低 SFRP1 可延长毛囊生长期、促进毛囊再生，这一结论已有离体人毛囊实验证实。

**科学问题**：头皮透皮给药后，毛囊内 siRNA 浓度常处于亚治疗水平。此时限制沉默效率的瓶颈不再是「互补程度」，而是低浓度下 **RISC 的催化周转效率**。同时，主流深度学习模型对「非核心区错配」序列缺乏训练覆盖，属于域外预测，不可单独作为排序依据。

**核心思路**：在引导链非核心区（g1、g12、g17、g18、g19）引入可控单点错配，可能提升 RISC 的催化周转；配合纯计算的多目标排序，在湿实验之前以计算替代试错，锁定值得进实验室的少数候选。

## 核心创新



1. **错配工程（Mismatch Engineering）**：将错配从「设计禁区」转为「可标定设计变量」。在 g1/g12/g17/g18/g19 五个非核心位点做单位点替换（1 条 WT + 15 条单点突变体），系统搜索错配空间以提升低浓度下的 RISC 催化周转。

2. **热力学主导 + 深度学习辅助的分层组合排序**：以引导链 MFE、双链末端 ΔΔG、seed 区结合能、整体双链 ΔG 四个热力学特征归一化加权作为主排序依据，OligoFormer 深度学习效率分为辅助信号，脱靶与毒性作为软惩罚项，兼顾可解释性与预测能力。

3. **自建错配专项模型 C\_match**：针对深度学习模型对错配序列「域外预测」不可靠的问题，自建错配代价模型（Ridge 错配代价校准 + Siamese CNN），留一来源（LOSO）跨来源验证 Spearman ρ = 0.611。

4. **工程化与可复现交付**：六步流水线模块化实现、固定随机种子、逐阶段 manifest 留痕防泄漏、74 项单元测试，`predict.py` 一键复现全流程。

## 技术路线



```
SFRP1 mRNA (NM\_003012)

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ① 候选生成  generation    │  1 nt 滑窗生成全部 19 nt 完全互补 siRNA

│                           │  + g1/g12/g17/g18/g19 五非核心位点单位点替换

│                           │  → 14,832 条候选

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ② 规则筛选  rule\_filter    │  GC 含量 / 连续相同碱基 / 连续 G/C /

│                           │  回文序列 四项串联过滤

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ③ 结构检测  structure      │  RNAfold 引导链自发折叠自由能 +

│                           │  双链末端稳定性差（强过滤，排除

│                           │  强发夹与链选择反转）

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ④ 热力学评估  thermo       │  可微最近邻模型逐位计算

│                           │  ΔG / ΔH / ΔS / Tm

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ⑤ 脱靶 & 毒性  offtarget   │  PITA + TargetScan seed 脱靶预测

│    / toxicity             │  cell\_viability seed 毒性 + 免疫刺激

│                           │  motif 扫描 + BLAST 基因组脱靶

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ⑥ 综合排序  rank           │  热力学四特征归一化加权（主）

│                           │  + OligoFormer 效率分（辅）

│                           │  − 脱靶 / 毒性软惩罚 → 最终排名

└──────────────────────────┘

&#x20;       │

&#x20;       ▼

┌──────────────────────────┐

│ ⑦ 化学修饰  chemmod        │  Top-50 候选清单 + 逐位化学修饰

│                           │  与递送建议

└──────────────────────────┘
```

## 验证与结果

全部为计算验证，可复现：



| 验证项    | 方法与结果                                                                       |
| ------ | --------------------------------------------------------------------------- |
| 公开基准回测 | 四套公开 siRNA 效率基准集（共 3,849 条）上，排序与实验沉默效率的平均 Spearman 秩相关由 **0.287 提升至 0.569** |
| 阳性对照回收 | 5 条已发表高效 SFRP1 siRNA 作为阳性对照，**4 条进入前 25%**，平均百分位约 14%                       |
| 错配模型泛化 | 自建 C\_match 错配模型：LOSO 跨来源验证 ρ = 0.611，LOFO ρ = 0.594                        |
| 消融实验   | 对比「错配 vs 无错配」及各位点贡献，量化非核心区错配的信号                                             |
| 权重稳健性  | 排序权重在公开基准上最小二乘校准，并做权重扰动稳健性扫描                                                |
| 工程验证   | 74 项单元测试（规则 / 结构 / 热力学 / 脱靶 / 毒性 / 排序 / 管线 / 化学修饰）                          |

## 快速开始



```
\# 1. 克隆并安装依赖

git clone \<your-repo-url>

cd siRNA\_pipeline

python -m pip install -r requirements.txt

\# 2. 运行完整管线（真实 SFRP1 全链 → 结果）

python predict.py

\# 3. 复现自建错配模型 / 运行单元测试

python train.py --model cmatch          # Ridge 错配代价校准

python train.py --model siamese         # 变体层 Siamese CNN

python -m unittest discover -s tests -v
```



* 小规模演示（合成数据，无需真实数据 / 外部工具）：`notebooks/demo_pipeline.ipynb`

* 换靶基因：`python predict.py --fasta <其他CDS.fa>`

* 断点续跑：`python predict.py --until structure` / `--from thermo --until chemmod`

* 配置：`configs/paths.yaml`（相对路径，零硬编码）、`pipeline.yaml`（阶段开关 / 随机种子）、`stages.yaml`（排序权重 / 阈值）

**环境要求**：Python ≥ 3.9（建议 3.10–3.13，已验证 3.13），无需 GPU，Windows / Linux 均可。完整管线必需外部组件：ViennaRNA 2.7.2、NCBI BLAST+ 2.14.1 + GENCODE v46、RNA-FM、OligoFormer、PITA / TargetScan（版本与安装见 `环境要求.txt`）。CPU 实测：无 ViennaRNA 约 60–90 秒，有 ViennaRNA 约 30–35 分钟。

## 目录结构



```
siRNA\_pipeline/

├── src/sirna\_pipeline/      # 核心代码（generation / rule\_filter / structure /

│                            #  thermo / toxicity / offtarget / oligoformer /

│                            #  rank / chemmod / pipeline / common）

├── configs/                 # 配置（paths / pipeline / stages）

├── scripts/                 # 分析、图表与验证脚本

├── tests/                   # 单元测试（74 项）

├── notebooks/               # demo\_pipeline.ipynb 演示

├── data/                    # 数据说明与派生表（来源 / 许可 / 清洗见 data/README.md）

├── models/                  # Model Card 与自训练模型说明

├── external/                # 第三方工具说明（版本 / 调用 / 许可）

├── outputs/                 # 结果（results / logs / runs / analysis）

├── final\_results/           # 大赛标准结果文件（results.csv / results.xlsx）

├── predict.py               # 主运行入口

├── train.py                 # 训练入口（cmatch / siamese / oligo\_head）

└── requirements.txt         # pip 依赖锁定
```

## 技术栈与依赖



| 类别   | 组件                                | 用途                                    |
| ---- | --------------------------------- | ------------------------------------- |
| 语言   | Python ≥ 3.9                      | 主开发语言                                 |
| 深度学习 | PyTorch ≥ 2.0                     | 热力学可微最近邻模块、Siamese CNN、OligoFormer 推理 |
| 结构预测 | ViennaRNA 2.7.2                   | RNAfold 二级结构 / 自由能                    |
| 序列嵌入 | RNA-FM                            | 预训练 RNA 语言模型特征                        |
| 效率预测 | OligoFormer                       | 开源预训练 Transformer，沉默效率推理（辅助信号）        |
| 脱靶   | PITA / TargetScan / BLAST+ 2.14.1 | seed 脱靶预测与基因组脱靶扫描                     |
| 毒性   | cell\_viability 数据库 + motif 扫描    | seed 毒性查表与免疫刺激检测                      |

## 输出结果



* `outputs/results/results.csv` — 大赛标准结果文件（UTF-8，Top-50：候选编号 / 赛道 / 候选序列 / 关键预测指标 / 模型版本 / 备注）

* `outputs/results/rank_final.csv`、`rank_top.csv`、`chemmod_top.csv` — 全量排序、Top 榜、化学修饰建议

* `outputs/runs/<run>/<n>_<stage>/` — 逐阶段候选表 + `manifest.json`（输入哈希 / 参数 / 统计，可溯源）

* `final_results/` — 汇总交付的 results.csv/results.xlsx

## 模型与第三方声明



* **主管道为可解释计算**（热力学最近邻前向 + 线性加权 + 规则），非黑盒；

* **自训练模型**：C\_match 变体层错配模型（Ridge 错配代价 + Siamese CNN），训练入口 `python train.py --model cmatch|siamese`；

* **第三方**：OligoFormer（开源预训练 Transformer，仅推理、资产缺失自动降级）、RNA-FM、ViennaRNA、BLAST+（版本 / 调用 / 许可见 `external/README.md`）；

* 排序权重在公开基准上最小二乘校准（LOSO 跨来源验证），文献正对照回收为独立验证；固定随机种子 42，逐阶段 manifest 留痕防泄漏。

## 可复现性



```
python -m unittest discover -s tests -v                 # 单元测试（无需外部工具 / 真实数据）

SIRNA\_REGRESSION=1 python -m unittest discover -s tests -v   # 真实数据回归

python train.py --model cmatch                          # 复现 C\_match 错配模型

python train.py --model siamese                         # 复现变体层 Siamese CNN

python scripts/update\_positive\_controls.py              # 复算文献正对照回收

python scripts/make\_pitch\_figures.py                    # 复现论文 / 答辩图
```

## License

本项目基于 [MIT License](LICENSE) 开源。第三方组件（OligoFormer、RNA-FM、ViennaRNA 等）遵循各自许可，详见 `external/README.md`。

## 致谢与参考数据



* 公开效率数据集：Holen 2005、Birmingham 2006、Ui-Tei 2008、Amarzguioui 2003、HIVsirDB、Sciabola 2013（Hu / Mix / Taka / Simone）等，派生表入库 `data/` 与 `outputs/analysis/`，原始文献不入库；

* 靶点序列：SFRP1 mRNA（NM\_003012），NCBI RefSeq 公共域；

* 感谢 OligoFormer、RNA-FM、ViennaRNA、BLAST 等开源社区的工作。



***

**方法学意义**：本流程与具体靶点解耦，可迁移至 Wnt 通路其他靶点、纤维化及肿瘤等更多靶点的核酸药物计算设计，为「以计算替代试错」的核酸工具设计提供可复现的方法学参考。