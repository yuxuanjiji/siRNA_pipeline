# 后续任务提示词（可直接转发给其他 Agent）

> 用途：把本文件中的"公共上下文 + 单条任务提示词"整体复制给其他 agent 执行。
> 每条提示词都自带：目标 / 改动点 / 数据与命令 / 验收标准 / 禁止事项。
> 生成时间：2026-09-11；配套结论文档：`docs/design/experiment_findings_2026-09-11.md`

---

## 0. 公共上下文（每条任务都请一并附上）

```
项目：SFRP1 靶向 siRNA 计算筛选与排序（生科挑战赛）
工作目录：D:\HuaweiMoveData\Users\及宇轩\Desktop\生科挑战赛\siRNA_pipeline
运行环境：Windows 10 · Python 3.13（系统解释器，含 torch 2.14 / pandas / matplotlib）
必设环境变量：VIENNARNA_BIN=C:\Program Files (x86)\ViennaRNA Package   # RNAfold CLI
OligoFormer 推理副本：siRNA_pipeline/OligoFormer部分（含 model/best_model.pth、RNA-FM 资产、.venv）
测试命令：python -m unittest discover -s tests      # 当前 71 项用例、7 项真实回归按开关跳过
断点/单测注意：测试临时目录请用 tests/_util.tmpdir（沙箱下 tempfile 会报权限错）

硬性约束（务必遵守）：
1) 旧模块零修改：src/sirna_pipeline/stages/*/legacy/、项目搭建/、BLAST/ 下的原文件不得改动；
   需要改行为时通过"新 runner / 适配层 / config 注入"实现。
2) 所有可选外部能力保持 graceful skip，缺环境不得让主管道失败。
3) 代码中不得硬编码绝对路径（含 Windows 盘符路径），路径一律走 configs/ 或环境变量。
4) 不得编造数据或指标：所有结论必须来自 outputs/analysis/*.json，并同时给出 n、95% CI 与 p 值。
5) 改动后必须跑 `python -m unittest discover -s tests` 全绿，且不改变既有 CLI 语义。

关键数据与产物：
- 验证脚本：scripts/run_experiments.py（子命令 check/benchmark/features/mismatch/ablation/positive/freeze）
- 绘图脚本：scripts/plot_validation_figures.py → outputs/analysis/figures/validation_*.png
- 结论 JSON：outputs/analysis/{task14_benchmark,task15_holen,task16_ablation}.json
- 阶段五数据：../数据集/阶段五_计算验证/{task14_fullmatch,task15_mismatch,task16_ablation,task17_positive_control}
- 错配特征表：OligoFormer部分/mismatch_ranked.csv（52 行，含 label/dG_total/dG_seed/ΔΔG_ends/oligo_pred/S_thermo/S_combo）
```

---

## P0-1 · 分段权重落地 + CV 复跑验证

```
【任务】按候选类型分段加权：完全互补（wt）候选以 DL 为主；含错配（mut）候选以 MFE 主导热力分为主、DL 权重 ≤0.2。

【依据】docs/design/experiment_findings_2026-09-11.md §4–5：
  - 错配集 n=52：MFE 单特征 ρ=0.353(p=0.014)；MFE 主导权重(1.8/0.4/0.4/0.2) ρ=0.360(p=0.0085)；
    现配置 ranker_default ρ=0.181(ns)；dG_seed 为负贡献(ρ=−0.130)；MFE 权重扫描 CV 最优 a=0.4。
  - 完全互补集：纯 DL 0.565–0.689，远高于纯热力 0.132/−0.192/0.150。

【改动点】
1) src/sirna_pipeline/stages/rank/ranker.py：新增按 kind（wt/mut）选择特征权重与 α/β 的能力，
   通过 cfg 读取（如 cfg["profiles"] = {"wt": {...}, "mut": {...}}），缺省保持现有行为向后兼容；
2) configs/stages.yaml 的 ranking 段：写入两个 profile（wt: DL 为主；mut: MFE 主导 + β≤0.2），
   保留原 weights/directions 字段以便回退；
3) meta.json 记录本行实际使用的 profile 与生效权重（可追溯）。

【验收】
- python -m unittest discover -s tests 全绿（新增 ranker 分段权重单测：wt/mut 各自取到正确权重、缺省回退一致）；
- 用 scripts/run_experiments.py ablation --input OligoFormer部分/mismatch_ranked.csv --with-mfe
  复跑，报告中 mut profile 的 Spearman ≥ 0.35 且给出 CI/p；
- 用真实 SFRP1 数据跑一次 predict.py，确认 results.csv 正常产出且 eligible 数量合理（≈4000 量级）。

【禁止】不要改 legacy 热力学/结构代码；不要把 n=52 的结论写成"已确证"（必须带 CI 与样本量）。
```

## P0-2 · 任务17 阳性对照：跑管道并出排名分位

```
【任务】把 7 条已发表 SFRP1 siRNA（../数据集/阶段五_计算验证/task17_positive_control/positive_controls.csv）
转成管道候选，跑完整流程，产出"是否通过 rules/structure + 综合排名分位"的结果表与图。

【现状】数据已备（含 CDS 定位与 57nt 上下文，7 条中 5 条可入管线，2 条靶 UTR 无法入管线），但尚未跑管道。

【改动点】
1) 新增 scripts 下的小工具（或用 run_experiments.py 现有 positive 子命令扩展）：
   把 positive_controls.csv 转成 pipeline 候选 CSV（沿用统一记录列：window_id/variant_id/kind/cds_start/target_mRNA_19/
   guide_checked/rules_pass 等），注意 guide 方向与 rc 口径一致（与 pipeline 生成规则相同）；
2) 跑 predict.py（可 --until ranking）得到全排序，再用 positive 子命令 join，输出
   outputs/analysis/task17_positive_controls.csv（含 final_rank 与 rank_percentile）与 JSON 摘要；
3) 新增图：outputs/analysis/figures/validation_task17_positive_controls.png —— 显示每条对照的 rank 分位
   （0%=最好），并标注是否通过 rules/structure。

【验收】5 条对照全部 matched；给出每条 final_rank / percentile；若某条排名靠后，需在结论中如实说明并分析原因。
【禁止】不得为了"好看"而挑选子集或改动排序权重来让对照组靠前。
```

## P0-3 · 错配集加入结构真值特征（cofold / 可访问性）

```
【任务】把 MFE 之外的结构特征加入错配评估：用 ViennaRNA 计算双链杂交结构相关量（端部配对概率、
内环数、凸起数）与靶标可访问性，检验能否在 n=52 上把 ρ 从 0.360 继续提升。

【改动点】
1) scripts/run_experiments.py 的 features 子命令扩展：新增 RNAcofold/RNAplfold 计算列
   （cofold_mfe、end5/end3 配对概率、internal_loops、bulges、target_accessibility），
   RNAfold/RNAcofold/RNAplfold 一律通过 VIENNARNA_BIN 或 PATH 定位（不得硬编码路径）；
2) 在 ablation 里把这些列纳入 feature_spearman 与 weight_presets（新增预设如 structure_plus_mfe）；
3) 结论写入 docs/design/experiment_findings_2026-09-11.md（追加一节，含 CI/p 与样本量）。

【验收】feature_spearman 中出现新特征且给出 CI/p；若某特征显著（p<0.05），给出与 MFE 的联合权重与 CV 结果。
【禁止】不得用结构模块 legacy 的私有函数直接改写；只读引用。
```

## P1-1 · 错配样本扩样 52 → 200+

```
【任务】把可用于排序评估的错配样本从 52 扩到 ≥200（口径统一后再计算 Spearman/AUC）。

【可用来源】../数据集/阶段五_计算验证/task15_mismatch/
  - 错配siRNA数据集_合并总表_处理后.csv（674 行；含 沉默效率 原始口径、57nt、来源、seed_pos）
  - desiRm 子集 55（效率为相关系数口径，需转换为序关系/秩）
  - Schwarz 10（IC50 nM → 相对活性秩）
  - Amarzguioui ~9、Ui-Tei 5（效率值需从图取数，建议 WebPlotDigitizer；不得猜数）
  - Holen 50 为金标准；Birmingham 362 为 seed 型脱靶，仅做 ROC/AUC，不并入 Spearman

【要求】
1) 每条样本必须标注：来源、效率口径、mismatch_side（guide/mRNA）、错配位点（可用现有推导器校验）；
2) 只把"同一口径可比"的样本进同一 Spearman 计算；不同口径分别出结果；
3) 扩样后重跑 MFE 权重扫描与权重预设，报告 n、CI、p。

【验收】扩样表写入 ../数据集/阶段五_计算验证/task15_mismatch/（含 README 说明），
并给出"扩样前(n=52) vs 扩样后"的对照表。
【禁止】不得用插值/猜数补齐缺失效率；缺就标注缺失。
```

## P1-2 · Taka 集异常复核

```
【任务】查清 task14 中 Taka 集"热力负相关（ρ=−0.192）"的原因：是 X 填充处理、效率归一化口径，还是数据本身特性。

【步骤】
1) 对比 Taka 与其他集的 mRNA 组成（X 占比）与效率分布；
2) 做 X→A 替换的敏感性分析：改为"仅用无 X 的 19nt 窗口计算 dG"或"跳过 X 行"，各算一次 ρ；
3) 与原始文献（Katoh & Suzuki 2007, NAR 35:e27）的效率口径核对；
4) 结论写入 findings 文档：若确认不可比，则把 Taka 标注为离群基准并在图注中说明。

【验收】给出至少 3 种处理下的 ρ 与 n 对照表，并给出明确结论与建议（保留/剔除/标注）。
```

## P1-3 · 工程健壮性两处修复（合并一条提示词）

```
【任务 A】OligoFormer 执行失败改为 fail-open：
  stages/oligoformer/adapter.py 中 executor/外部命令非零退出时，不再抛错中断管道，
  而是写 oligo_available=0、记录原因到 manifest/summary，继续排序；新增 cfg 开关 on_error: skip|raise（默认 skip）。
【任务 B】长任务日志可见性：把 executor 子进程 stdout/stderr 落盘到 outputs/logs/oligoformer_executor.log（成功也写），
  并在 adapter 日志里输出其尾部若干行。

【验收】
- 新增单测：模拟 executor 返回非零 → 管道不失败、oligo_available=0、manifest 含 error 原因；
- 真实跑一次 predict.py 确认无行为回归（aux 可用时仍正常并入 score_oligo）；
- 全量测试全绿。
```

## P1-4 · PITA 层写出 offtarget_risk（让脱靶惩罚真正生效）

```
【任务】stages/offtarget/adapter.py 目前不写 offtarget_risk，导致 ranker 的脱靶软惩罚在默认配置下恒为 0。
请把 PITA/TargetScan 结果映射为 offtarget_risk（任一 filter/pass 命中 → 1，否则 0），
并在 manifest 记录映射规则；保持默认关闭时不影响任何现有结果。

【验收】单测：构造带 offtarget_pita_filter/targetscan_filter/offtarget_pass 的外部结果 →
offtarget_risk 正确置位；rank 阶段对该列扣分符合预期（penalty_total 增加）。
```

## P2-1 · 任务19 正式版成药性论证卡

```
【任务】用真实 Top 候选替换 ../项目搭建/化学修饰/任务19_成药性论证卡.md 里的 21nt 占位（P1/P2/P3），
产出正式版：化学修饰方案（沿用 stages/chemmod 的 ESC-19 引擎输出，如 results/chemmod_top.csv 的前 3 条）+
递送策略 + 可制造性评估，并保持 19nt 主链口径一致（21nt 化作为后续迭代说明）。

【验收】文档中每条候选的 guide/sense 修饰记法与 chemmod_top.csv 完全一致（逐字符核对），
参考文献沿用既有核实状态标注（[已核实]/[待复核]）。
```

## P2-2 · 任务20 交付件（归档 + PPT + 摘要）

```
【任务】
1) 仓库整理：确保 README/项目说明/环境配置与下载清单/模型卡/设计文档齐备；压缩包排除 .venv、__pycache__、
   outputs/runs、OligoFormer部分/RNA-FM 与 data（给出排除清单与重建说明）；
2) 摘要：按附件2 写 300–500 字 + 3–5 关键词，不得出现校名/姓名/LOGO；
3) PPT：16:9、全中文、≤8 分钟录音；页面与台词中不得出现单位/校徽/姓名/照片；
   内容须与评审表一致；图表使用 outputs/analysis/figures/validation_*.png（新口径）。

【验收】清单化交付物路径 + 一份"匿名合规自查表"（逐条打勾）。
```

## P2-3 · 文档与图件收尾

```
【任务】
1) docs/design/experiment_inputs.md 补充 features 子命令与 --with-mfe 的输入契约；
2) README「大赛对照」章节补"验证结果"小节（引用 findings 文档与 validation_* 图）；
3) 图件：在 plot_validation_figures.py 中新增一版"PCC + Spearman 双指标"图（用 task14 JSON 的 pearson/r²），
   并在 findings 文档图件清单中登记；旧图 task14_benchmark.png / task15_mismatch.png 标注为"旧口径"。

【验收】文档互链可达；重新运行绘图脚本后 validation_* 全部更新且无字体告警。
```

---

## 通用验收提示词（每次改动后追加发送）

```
请对本次改动做自检并给出证据：
1) python -m unittest discover -s tests 输出摘要（Ran/OK/FAILED 行原文）；
2) 受影响 JSON 的关键指标（含 n、95% CI、p 值）前后对照表；
3) 是否触碰 legacy/原文件（应为"否"）、是否出现硬编码绝对路径（应为"无"）；
4) 若涉及结论，明确写出"样本量 + 置信区间"与"不显著项"。
```
