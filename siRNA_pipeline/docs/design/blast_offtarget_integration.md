# BLAST 近全长同源脱靶层 · 集成说明

> 将 `siRNA_pipeline/BLAST/`（原脚本）作为**脱靶检测模块的第三层**接入主管道。
> 原文件零修改：`BLAST/blast_offtarget.py` 与 `BLAST/README.txt` 一字未动；
> 接入代码全部在 `stages/offtarget/blast_adapter.py`（新适配层），沿用既有
> 「配置启用 + 状态探测 + graceful skip」适配模式（同 PITA/TargetScan 层）。

## 1. 模块定位

脱靶验证三件套（与既有适配层平行、可叠加）：

| 层 | 检测内容 | 状态 |
|---|---|---|
| PITA/TargetScan（种子区） | seed 脱靶打分/过滤 | `stages/offtarget/adapter.py`，默认关 |
| （热力学可及性） | —— | 由 03 structure 的 plfold 承担 |
| **BLAST 近全长同源** | 转录组中"含 rc(guide) 近全长位点"的非靶基因 | `stages/offtarget/blast_adapter.py`，默认关 |

BLAST 层判据（tolerant 口径，见 `BLAST/README.txt`）：任一命中同时满足
非靶基因 / 反义方向(sstart>send) / 配对数 (length−mismatch) ≥17 / mismatch ≤1 → 判脱靶。

## 2. 接入方式（默认不影响主管道）

- 配置：`configs/pipeline.yaml → stages.offtarget.blast`（默认 `enabled: false`，
  与 stages.yaml 同键合并）。未启用或环境（blastn/库）未就绪时，适配层
  **graceful skip**：写出 `offtarget_blast_available=0`，不阻断、不伪造判据；
- 编排：orchestrator 的 offtarget 阶段在 `blast.enabled=true` 时，先（可选）跑
  PITA 层、再跑 BLAST 层串接；产物统一为 `05_offtarget/candidates_offtarget.csv`
  （或仅 PITA 时维持原文件名），供 rank 阶段 join；
- 回填列（`common/records.py` 只增列）：
  `offtarget_blast_flag`（1=判脱靶）、`offtarget_blast_n_hits_ge16`（审计计数）、
  `offtarget_blast_available`，并写共享列 `offtarget_risk`（=flag，rank 已默认
  将其作为 0.10 软惩罚项，开启 BLAST 后无需改 rank 配置）；
- 顺序语义：BLAST 层在 PITA 层之后执行；`offtarget_risk` 取“任一命中”并集。

## 3. 启用样例（环境就绪后）

```yaml
# configs/pipeline.yaml → stages.offtarget
offtarget:
  blast:
    enabled: true
    blast_bin: "E:/大创2026/blast/ncbi-blast-2.14.1+/bin"   # 或留空走 PATH
    db: "E:/blastdb/gencode_v46_pc_simple"                  # makeblastdb 前缀
    target_gene: SFRP1
    min_paired: 17
    max_mismatch: 1
    python: null        # null=当前解释器（需装有 pandas）
    script: null        # null=自动定位 工程根/BLAST/blast_offtarget.py
```

环境要求与依赖声明见 `external/README.md` §5 与 `BLAST/README.txt`（NCBI BLAST+ 2.14.1、
GENCODE v46 库许可需自查）。

## 4. 单测与不破坏性保证

- `tests/test_blast_offtarget.py`：默认跳过 / 环境缺失优雅降级 / 桩脚本契约回填
  （flag、n_hits、available、risk）；
- 主管道默认全关 → 既有 6 步流程、排序、化学修饰与标准结果文件行为不变
  （全量单测 + 真实 SFRP1 全链回归验证）。

## 5. 与《附件5》一致性的备注

本层仍遵守：原文件零修改、绝对路径不硬编码（bin/db 走 config）、结果可溯源
（manifest 记录 db/blastn/判据/计数）、可选能力默认关闭且不影响主流程。
