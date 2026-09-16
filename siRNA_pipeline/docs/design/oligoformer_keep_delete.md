# OligoFormer 保留/删除取舍说明（efficiency 辅助分接入）

> 配套技术路线第五步："将 OligoFormer 模型作为辅助（负责沉默效率）"。
> 原则：**不删除/不改动任何旧文件**（项目铁律：`项目搭建/…` 原地保留），
> 主管道只新增"精简调用层"（适配器 + 归一化 + 配置），并在此说明每部分的取舍与原因。

---

## 1. 结论摘要（先读）

| 类别 | 内容 | 处理 |
|---|---|---|
| **保留（供参考/调用）** | 推理主链：`model/best_model.pth` + `scripts/model.py`(Oligo 模型) + `scripts/loader.py`(infer 侧 `data_process_loader_infer`) + `scripts/RNA-FM.sh`（RNA-FM 嵌入）；`scripts/infer.py`（只读参考其 `efficacy` 口径与 td 计算）；`toxicity/cell_viability.txt` | 原位保留；管道内置 executor 以子进程/import 复用其语义（见 §3） |
| **保留（参数已收敛）** | `infer.py` 顶部 `DeltaG/DeltaH` 最近邻表 | 数值与 `siRNA_pipeline/common/nn_tables.py` **逐字一致**（Xia 1998）；管道新代码只引用 `common/nn_tables.py` |
| **不并入主管道（保留原位）** | 训练/评估：`scripts/train*.py / test*.py / analyze_*.py`、`predict_mismatch.py`、`data/`、`Comparison methods/`（Monopoli-RF、siRNAPred） | 研究/对照用途；与"筛选流水线"解耦，避免强环境耦合与体积拖累 |
| **不并入（适配层接管）** | OligoFormer 自带 off-target：`scripts/pita.sh / targetscan.sh` + `off-target/` | 见 `siRNA_pipeline/stages/offtarget/adapter.py`（默认 disabled，graceful skip） |
| **新写精简调用层（管道内）** | `stages/oligoformer/adapter.py`（适配/探测）+ `executor.py`（内置推理）＋ rank 的 `score_oligo`/α·β 归一化 + `OLIGO_COLS` | 见 §3；**默认 enabled=true**，缺模型/RNA-FM 资产时优雅降级（α=1） |
| **删除建议（仅清理建议，未执行）** | `.venv/`、`data/` 演示/中间产物、`Comparison methods/` 若干大数据 | 如需腾空间由负责人手工执行；按铁律本轮零删除 |

---

## 2. 我们"保留"了什么：推理主链（efficiency 语义）

`scripts/infer.py` 中与本项目第五步"沉默效率辅助分"直接相关的语义只有一处：

```python
pred, _, _ = best_model(siRNA, mRNA, siRNA_FM, mRNA_FM, td)   # 分类 logits
efficacy  = pred[:, 1] * 1.341        # “高效类”概率 ×1.341（对齐 PITA 分数量纲，越大越好）
```

即 OligoFormer 对每条 (guide, mRNA-57nt) 输出一个**越高越好的效率分**。管道保留这一语义：

- 输入侧需要：siRNA 序列 + 57 nt mRNA 上下文（管道统一记录已具备：`guide_checked` + `cds_start`→`mrna57_from_cds`）；
- 模型侧需要：`model/best_model.pth`、`Oligo` 网络定义（`scripts/model.py`）、RNA-FM 预训练嵌入（`scripts/RNA-FM.sh`，bash + 外部 FM 权重）；
- 输出侧：一列 `oligo_efficacy`（越大越好）+ `oligo_available`（0/1 可用标记）。

## 3. 管道内新写的"精简调用层"（保留语义、隔离环境）

`siRNA_pipeline/stages/oligoformer/adapter.py` + `executor.py`（内置执行器，同目录）：
按官方 infer 口径对统一记录逐条出 `oligo_efficacy`（RNA-FM 嵌入 → loader 数据集 →
Oligo 模型 → `pred[:,1]×1.341`）。两种模式：

```
内置 executor（默认）：python executor.py --input <csv> --output <csv>
                        --repo <OligoFormer部分> --model <best_model.pth> --cds <CDS.fa>
外部执行器契约（可选）：external_cmd <cmd> --input <csv> --output <csv>
输出 CSV 含 variant_id 列 + oligo_efficacy/oligo_available（可子集）
```

- `enabled: true`（configs 已默认开）：模型权重 / OligoFormer .venv / **RNA-FM 资产**
  （`<repo>/RNA-FM`，含 pretrained/extract_embedding.yml）齐全才真正推理；任一缺失 →
  写出 `oligo_available=0`、不伪造分数（graceful skip，状态原因在 manifest/summary 可见）；
- rank 阶段消费 `oligo_efficacy`：min-max（direction +1）→ `score_oligo` → `β·score_oligo` 并入终分；
  模型不可用/列缺失时 **β 自动置 0、α 回补为 1**（`stages/rank/ranker.py` 的 `aux_used` 逻辑），
  排序结果与纯热力口径完全一致、可复现（见 `docs/design/ranking_weighting.md`）。

## 4. 我们"删除/不并入"了什么（含原因）

| 部件 | 位置 | 取舍原因 |
|---|---|---|
| 训练/微调与全套评估 | `scripts/train.py, train_single.py, train_logger.py, test.py, test_single.py, analyze_*.py` | 训练与研究用；流水线只做推理，不承担训练环境 |
| 错配预测脚本与演示数据 | `predict_mismatch.py`、`data/mismatch_demo/`、`check_data.py`、`flanking_*.py` | 研究/对拍用（错配校准已在 `siRNA_pipeline/stages/thermo/legacy/calibrate_mismatch.py` 独立承载） |
| 比较方法大目录 | `Comparison methods/Monopoli-RF`、`siRNAPred` | 基准对照（综述/消融报告用），与主管道无关 |
| 自带脱靶工具链 | `scripts/pita.sh / targetscan.sh`、`off-target/ref|targetscan` | Perl/bash + ViennaRNA 强耦合 → 适配层默认关（见 `offtarget/adapter.py`） |
| RNA-FM/bash 链路 | `scripts/RNA-FM.sh` | 效率推理的最重外部依赖（bash+FM 权重）；不在主管道内 `os.system`，走外部执行器 |
| `.venv/`、`requirements/environment.yml` | 目录根 | 属于本机运行环境，不随工程发布 |

> 说明：以上"不并入"均指**不进入 `siRNA_pipeline` 的运行链/导入链**，旧文件一个未删。
> 若负责人后续希望物理清理大文件（`.venv` 等），按 §1 删除建议清单手工执行并留痕即可。

## 5. 当前状态与开启/补资产步骤

**状态（2026-09-09 起）**：`configs/pipeline.yaml → stages.oligoformer.enabled: true`，
内置 executor 已接线、默认即启用；rank 在 `oligo_efficacy` 可用时自动并入辅助分，
否则优雅降级（α=1）。

在本机**让效率分真正生效**只需补齐 RNA-FM 资产（模型权重与 .venv 已在
`工程根/OligoFormer部分/`）：
1. 在 `<OligoFormer部分>/RNA-FM/` 放置 RNA-FM 代码与预训练（官方
   `launch/predict.py` + `pretrained/extract_embedding.yml` 层级；需 fairseq/torch 环境，
   网络下载约 1–2GB，注意其数据许可与引用：RNA-FM 论文 Chen et al. 2022, Nat Mach Intell 4:759）；
   或自行安装后确认 `bash scripts/RNA-FM.sh <data_dir>` 能跑通；
2. （可选）覆盖 `configs/pipeline.yaml → oligoformer`：`python / repo_dir / model_path /
   rnafm_required`；
3. `python predict.py --run-name myrun` 即可；rank 的 meta.json 会记录
   `aux.available=true` 与 α/β 生效值，便于核对与复现。

## 6. 一致性声明

- 保留的 NN 表数值已由 `siRNA_pipeline/tests/` 中 thermo 对拍断言与 `infer.py` 表一致
  （`common/nn_tables.py` 头注释声明同源）；本说明只描述取舍，不替代代码内注释与测试。
