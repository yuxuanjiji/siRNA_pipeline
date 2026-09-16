# 第三阶段：Pipeline 设计（Phase 3 · Pipeline Design）

> 数据流决定：**任何候选数据只以统一记录表（CandidateRecord CSV）形式在阶段间流动**；
> 阶段输出都带 `<variant_id>` 键可 join；记录“列只增不删”。

## 3.1 Pipeline 总览

```
[输入] configs/paths.yaml 数据根
   │
   ▼
1. GENERATION(任务5/6, 纯标准库)     输入: 数据集/SFRP1-mRNA.txt
   │  输出: runs/01_generation/candidates_wt.csv + candidates_mut.csv(合并 candidates.csv)
   ▼
2. RULE_FILTER(任务7, 纯标准库)       输入: 01 输出
   │  输出: runs/02_rules/candidates_rules.csv (含 rules_pass) + candidates_passed.csv(rules_pass=1)
   ▼
3. STRUCTURE(任务8, ViennaRNA 可选)   输入: 02 passed
   │  输出: runs/03_structure/candidates_structure.csv (structure_pass/reject_reason + 特征列)
   ▼
4. THERMODYNAMICS(任务11/12, torch)   输入: 03 输出(或 02 passed 直算)
   │  输出: runs/04_thermo/candidates_thermo.csv (+派生 feat_* 列)
   ▼
5. OFF-TARGET(任务9, Adapter, 默认关闭) 输入: guide(seed)
   │  输出: runs/05_offtarget/candidates_offtarget.csv (offtarget_*; 关闭时整列留空)
   ▼
6. TOXICITY(任务10, 标准库+本地表)    输入: 04(或02) 记录
   │  输出: runs/06_toxicity/candidates_toxicity.csv (tox_*/imm_*)
   ▼
7. RANKING(任务11–13, 新模块)         输入: 03/04/05/06 各表 join(按 variant_id)
   │  输出: outputs/results/rank_final.csv + meta.json(权重/统计/淘汰数)
   ▼
[交付] Top-N 列表 + 报告素材(淘汰原因分布、四特征分布)
```

**去重记录原则**：structure 与 thermo 都只需 19nt guide（+57nt 上下文重建自 CDS），故 03 与 04 可各自独立从 02 passed 出发并行；05/06 同理。Join 只在 RANKING 阶段发生。

## 3.2 每步明细（输入/输出/调用/数据结构/配置/错误处理）

| 步 | 输入 | 输出 | 调用模块/函数 | 输入结构 | 输出结构 | 配置 | 错误处理 |
|---|---|---|---|---|---|---|---|
| 1 generation | FASTA | candidates CSV(wt+mut 16/窗口) | `stages.generation.runner` → 复用 task5/6 核心(原逻辑不变) | FASTA | CandidateRecord 01 | 基因/CDS 路径、滑窗长19、突变位点[1,12,17,18,19] | 输入文件缺失/非法碱基→抛错并写 stage 错误日志；校验失败(窗口数/每窗口15条)即中止该阶段 |
| 2 rules | 01 | 02(+passed) | `stages.rule_filter.runner` → task7 `build_candidates/evaluate` | 记录表 | +rules_* | GC 区间、run 阈值、回文臂长(默认=OligoFormer 口径) | 行级 evaluate 有 ValueError 则记 bad rows 跳过；库级校验失败中止 |
| 3 structure | 02 passed | 03 | `stages.structure.runner` → `StructureDetector.detect_batch` | 记录表(19nt+mRNA57重建) | +structure_*/pass/reject | MFE/-6、端ΔΔG/-0.5、cofold≤-10、可及性 None、strict 开关、ViennaRNA 定位 | ViennaRNA 缺失：自动降级标记 `structure_reliable=False` 并在日志告警(不静默)；单条异常行→该行 NaN 并记原因 |
| 4 thermo | 02 passed | 04 | `stages.thermo.runner` → `siRNAThermoCalculator.calculate_batch` | 记录表 | +dG_*/Tm +派生 feat | 默认修正系数(与 calibrate json 一致)、温度37、seed 口径 | torch 不可用→阶段报错退出；非法序列行记 NaN |
| 5 offtarget(关) | 04/02 | 05 | `stages.offtarget.adapter` | guide 表 | +offtarget_* | `enabled:false`、工具路径、阈值 | 未启用→输出“整列空”并注明；启用且缺 Perl/bash→显式报错并建议关闭 |
| 6 toxicity | 04/02 | 06 | `stages.toxicity.runner` → `ToxicityDetector.detect_batch` | 记录表 | +tox_*/imm_* | viability 表路径、阈值50、motif 规则开/关 | 表缺失→按 config 报错或降级(viability 列空) |
| 7 ranking | 03+04+05+06(join) | rank_final.csv+meta.json | `stages.rank.rank` | 宽记录 | +score/rank | 特征权重、归一化方式、惩罚项、Top-N | join 缺行→告警留空；四特征缺→按策略(跳过/报错) |

## 3.3 Pipeline 主入口与阶段运行器（约定）

```
scripts/run_pipeline.py            # 唯一入口
  --config configs/pipeline.yaml   # 阶段开关(offtarget 默认 false 等)
  --paths   configs/paths.yaml     # 数据/旧目录/输出根
  --until   阶段名(可选, 断点续跑: generation|rules|structure|thermo|offtarget|toxicity|ranking)
  --from    阶段名(可选, 只跑某段)
  --overwrite/--keep  (已有 runs 目录处理策略)
```
- 每个 stage 是一个独立可 import 的运行器模块 + `main(argv)`；阶段通过 `StageResult`(manifest 字段)向编排器报告行数/通过数/耗时/校验结果。
- **可回退/断点**：每阶段写 `manifest.json`(输入文件 hash、参数、版本、状态)；重复运行相同入参时可 `--keep` 跳过已成功阶段。
- **退出码**：0=成功；2=配置错误；3=该阶段数据处理失败；4=环境缺失(ViennaRNA/bash/torch)；日志写 `outputs/logs/pipeline_YYYYmmdd_HHMMSS.log` + 每阶段 `stage_<n>_<name>.log`。

## 3.4 配置管理（不再硬编码）

```
configs/
├── pipeline.yaml   # 阶段开关、每窗口16条断言开关、seed口径、continue_on_error
├── stages/
│   ├── generation.yaml   # window_len=19, sites=[1,12,17,18,19], fasta 引用 paths
│   ├── rules.yaml        # gc:[30,65], run_gc:6, run_same:5, palindrome_arm:4
│   ├── structure.yaml    # mfe_threshold:-6.0, end_diff:-0.5, cofold_max:-10, strict:false, backend:auto, vienna_bin_dir:null
│   ├── thermo.yaml       # overrides:{}（默认=模块 DEFAULT_CORRECTION）, temperature:37
│   ├── offtarget.yaml    # enabled:false, pita_threshold, targetscan_threshold, seed_window:6
│   ├── toxicity.yaml     # viability_threshold:50.0, hard_rules:{}, motifs 开
│   └── ranking.yaml      # features:[mfe,ddg_ends,dG_duplex,dG_seed], weights:{}, norm:minmax,
│                         #   direction:{}, penalty:{offtarget_soft, tox_viability, imm_high}, top_n:20
└── paths.yaml            # data_root: ../数据集, ref_root: ../参考模型, legacy_root: ../项目搭建,
                          # runs_dir/outputs_dir/external_dir(相对本文件或绝对, 支持 env 覆盖)
```

## 3.5 日志与输出布局

```
outputs/
├── logs/     pipeline_*.log · stage_<n>_<name>.log（logging: 时间/级别/阶段/行号）
├── runs/     01_generation/…07_ranking/ 各自 candidates_*.csv + manifest.json + 统计 json
└── results/  rank_final.csv · rank_meta.json · 报告用汇总 csv
```
*每条关键决策/校验结果都写入阶段日志与 manifest，保证“可追踪、可解释”。*

## 3.6 错误与一致性策略

- 数值列写入：`None/NaN→空串`；读回统一 float/None——schema 由 `common/records` 校验器保证。
- 每阶段入口做“上游断言”（例：rules 阶段要求存在 01 且 16 条/窗口 校验通过）。
- 不修改原始任务脚本：全部通过“新 runner 内 import 原核心 + 记录转换”完成（原始文件只读拷贝于 `src/…/legacy/` 或引用原位，见 Phase 4）。
