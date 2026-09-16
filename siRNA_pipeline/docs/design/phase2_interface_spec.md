# 第二阶段：模块接口规范（Phase 2 · Module Interface Specification）

> 工程：`siRNA_pipeline`（工作区根 `D:\…\生科挑战赛\siRNA_pipeline`）
> 依据：第一阶段全量代码审查（原代码一律不改、原地保留；本阶段仅定规范）
> 数据根默认引用：`../数据集`（原 SFRP1‑mRNA.txt、mismatch 数据、基准效率集），`../参考模型`（OligoFormer / siRNADiscovery），`../项目搭建`（旧模块只读参考）

---

## 2.1 模块接口注册表

图例：**当前接口** = 该模块现状；**推荐接口** = 整合后统一形态。所有 CSV 均 UTF‑8‑sig(BOM)。

| 模块(任务) | 输入 | 输入格式 | 输出 | 输出格式 | 依赖 | 当前接口 | 推荐接口 | 是否需要修改 |
|---|---|---|---|---|---|---|---|---|
| generation·滑窗(5) `task5_…py` | SFRP1 CDS FASTA(单条) | FASTA(DNA/RNA) | 927 行 19nt 全互补窗口表 | CSV | 仅标准库 | CLI `--fasta/--out/--window`；`sliding_windows()/rc_rna()` 供 import | runner 直接调用库函数，窗口表 → 统一记录 CSV（stage=generate, kind=wt） | 不改算法；包装 |
| generation·定点突变(6) `task6_…py` | 任务5 表(或 FASTA) | CSV 列 `guide_antisense_19…` | 13 905 行错配变体库 | CSV | 仅标准库；import task5/seq_utils | CLI `--input/--fasta/--out/--g1-wobble-only` | 同上 → 统一记录追加 15 条 mut/窗口 | 不改算法；包装 |
| rule_filter(7) `task7_…py` | 任务5+6 表 | CSV(列集含 `guide_mut_19/guide_antisense_19`) | `_annotated.csv`(14832) + `_passed.csv`(6524) | CSV | 仅标准库 | CLI `--library/--wt/--no-wt/--out-dir/--prefix`；默认路径猜 候选序列生成/ | 吃统一记录，输出在记录上追加 `rules_*` 列并落 passed | 不改算法；包装+映射列名 |
| structure(8) `structure_detector.py` | (guide19[, mrna57]) | Python 字符串/列表 | dict(键见 §2.2‑C) | dict（不落盘） | ViennaRNA 绑定 `import RNA` 或 CLI RNAfold/RNAcofold/RNAplfold；无则纯 Python 近似 | 类库 `StructureDetector.detect/detect_batch/filter_candidates`；仅有 `__main__` 自检 | 新增 stage runner：读统一记录 → 重建 57nt → detect_batch → 追加 `structure_*` 列 + `structure_pass/reject_reason` | 类不改；新增批量入口 |
| thermodynamics(11/12) `thermo_calculator.py` | (guide19, mrna57) | 字符串/张量 | dict：ΔG/ΔH/ΔS/Tm + 逐位分解 | dict/Tensor | torch | `siRNAThermoCalculator().calculate/calculate_batch/forward` | 保留为唯一热力学实现；runner 追加 `dG_*` 列与派生排序特征列 | 仅新增派生输出(加法)；doubao 版归档 legacy |
| thermo·校准 `calibrate_mismatch.py` `validate_mismatch_dataset.py` | `数据集/mismatch.csv` + 校验表 | CSV | `mismatch_calibration.json` / `mismatch_validated.csv` | JSON/CSV | numpy+scipy / 标准库 | CLI `python calibrate_mismatch.py` | 移入 analysis/，保持 CLI，只规范化路径；产出入 `data/processed/` | 路径/归档，不改逻辑 |
| toxicity(10) `toxicity_detector.py` | guide(+passenger) 或 TSV | 字符串 / TSV | dict(`tox_*`+`imm_*`) / TSV | dict/TSV | 标准库 + 本地 `data/cell_viability.txt`(4096) | 类库 + CLI `--self-test/--trial/--input/--output` | runner 追加 `tox_*` 列；viability 表路径走 config | 类不改；runner 包装 |
| offtarget(9) 【未实现】 | guide 19nt(seed) | 待定 | PITA/TargetScan 分数与过滤 | 待定 | OligoFormer off‑target(Perl/bash)+ViennaRNA | 无 | **Adapter 接口**：`OffTargetDetector.available()/detect_batch()`；环境不满足→配置级跳过 | 全新（适配层先行，工具链后补） |
| ranking(11–13) 【未实现】 | 四特征 + 惩罚项 | 统一记录列 | 归一化特征、热力粗排、软惩罚、最终排名 | CSV/JSON | numpy(或 pandas 可选) | 无 | `rank/`：config 定义权重/归一化/惩罚；输出 `results/rank_final.csv` + `meta.json` | 全新（按文档算法） |
| efficiency·DL 辅助 | OligoFormer / siRNADiscovery 输入格式 | FASTA / td 串 / candidates.json | 效率分 | 分数列 | torch / TF+stellargraph+RNA‑FM | 参考仓库自带入口（含 bash 调用） | 仅作**外部可选步骤**：config 开关 + Wrapper；默认关闭 | 不并入本轮主线 |
| chemmod(19) | Top 候选 | 记录表 | 修饰方案/论证卡字段 | 待定 | — | 占位脚本(21nt，硬编码绝对路径) | 不并入主线；占位保留原位，待任务19排期 | 本轮不动 |

## 2.2 数据格式逐项核对（检查点 1–9 结论）

1. **上一模块输出能否直喂下一模块**：任务5→6→7 可以（列名一致，7 已自带适配：`guide_checked`=wt 的 `guide_antisense_19`/mut 的 `guide_mut_19`）。**断点在 7→8**：`_passed.csv` 无人消费；structure 无文件 IO。**8/10→ranking**：无接收方。
2. **数据类型**：全管道字符串序列一致（RNA 大写、U）。数字：task5/6 写字符串，task7 追加 0/1 整型与浮点 GC；structure/thermo 返回原生 float/int → **runner 统一按“可空数值列”规范写入**（空用空串，读时转换），记录层做 schema 校验。
3. **序列格式/几何**：统一 19nt 引导链（5′→3′），mRNA 窗口 19，外接 57=19+19+19，配对 `guide[p]↔target[18-p]`；唯一不一致点：toxicity 头部注释假定 passenger=rc(guide)（见 §2.4 待裁定项#3）与 chemmod 21nt 占位。
4. **参数名**：跨模块别名需映射——`guide_antisense_19`（t5）/`guide_wt_19`（t6）/`guide_checked`（t7）是同一物理量；`pair_type`（t6：GU/MM）与 thermo 的 WC/GU/PP/YY/PY 分类粒度不同（t6 不区分 PP/YY/PY）→ 统一记录同时保留两粒度：`pair_type`（粗）与可选 `mismatch_class`（细，由 thermo/calibrate 的 classify 补充）。
5. **文件路径**：旧模块默认路径 = “相对本文件猜邻居”，只读引用没问题，写默认输出全部改为显式 config/参数；新工程**只写 `outputs/runs/<stage>/`**，绝不写回旧目录。
6. **模型输出格式**：structure/thermo/toxicity 均返回 dict → runner 按白名单列展开；分数型列统一 float，标志型列统一 `0/1/空`。
7. **相同功能**：RC/GC/配对分类/CSV(BOM) ≥9 处复制；NN 表 3 份（structure/thermo/doubao）；错配 classify 3 份 → 收敛到 `common/seqio.py`、`common/nn_tables.py`、`common/mismatch_class.py`，各模块只 import 不复制（对旧文件保留原样，只在新 runner 中引用共享实现，供长期收敛）。
8. **同名函数/文件冲突**：无同名文件；函数同名不同实现者（如各 `rc`/`normalize`/`classify`）在新包内统一命名空间，避免 `from x import *`。
9. **强耦合点**：a) task6 import task5（同目录，需成对放置）；b) task7 的默认输入路径硬指 候选序列生成 目录；c) OligoFormer infer.py 内部 `os.system('bash …')` + 仓库相对路径（最强外部耦合，仅适配层面对）；d) `calibrate_mismatch.py` 以 `__file__` 推导父目录读 `数据集/mismatch.csv`（结构脆弱但只读）。

## 2.3 统一数据结构：CandidateRecord（列模型）

原则：**列只增不删；阶段白名单校验**；key = `(window_id, variant_id)`；variant_id 空即 wt。

| 分组 | 字段（统一名） | 必/选 | 生产者 | 消费者 |
|---|---|---|---|---|
| 主键/坐标 | `window_id`(`W0001`) | 必 | t5 | 全部 |
| | `variant_id`(`W0001_g1:A>U`；wt 记 `{window_id}_wt`) | 必 | t5/t6 派生 | 全部 |
| | `kind`(`wt`/`mut`) | 必 | t5/t6 | t7、消融 |
| | `cds_start/cds_end/nm_003012_start/nm_003012_end` | 必 | t5 | 57nt 重建、报告 |
| 序列 | `target_mRNA_19` | 必 | t5 | 配对/结构/热力 |
| | `guide_wt_19` | 必 | t5 | 回比 |
| | `guide_checked`（被筛/被算序列） | 必 | wt=t5,mut=t6 | t7/8/10/11 |
| | `sense_strand_19` | 选 | t5 | 修饰/报告 |
| 突变注释(mut) | `site/position/wt_nt/mut_nt/paired_mRNA_nt/pair_type` | 条件必(mut) | t6 | 消融、热力预标注 |
| 规则筛选 | `gc_pct` | 选 | t7 | 报告 |
| | `fail_gc_range/fail_run6_gc/fail_run5_same/fail_palindrome` | 选 | t7 | 统计 |
| | `hit_code/hit_rule` | 选 | t7 | 对拍、统计 |
| | `rules_pass`(旧 `seq_pass` 别名) | 必 | t7 | 8/9/10/11 过滤语义 |
| 结构 | `mfe/mfe_structure/max_stem_length/max_stem_gc` | 选 | structure | 排序特征/消融 |
| | `delta_G_5end/delta_G_3end/delta_deltaG_ends` | 选 | structure | 排序特征 |
| | `cofold_mfe/target_accessibility/…` | 选(可空) | structure | 仅强过滤时 |
| | `structure_reliable` | 选 | structure | 质量标记 |
| | `structure_pass/reject_reason/warnings` | 选 | structure | 过滤+消融 |
| 热力学 | `dG_total/dH_total/dS_total/Tm` | 选 | thermo | 报告 |
| | `dG_core_canonical/dG_mismatch_total/dG_flank_total` | 选 | thermo | ΔΔG 归因 |
| | `dG_mismatch_by_position`(JSON串,len19) | 选 | thermo | seed 特征/消融 |
| 排序特征(派生) | `feat_mfe/feat_ddg_ends/feat_dG_duplex/feat_dG_seed` | 条件必(进入排序前) | 见 §2.5 | rank |
| 脱靶 | `seed6`/`offtarget_score/offtarget_pass`(`pita/targetscan_*`) | 选 | offtarget(Adapter) | rank 软惩罚 |
| 毒性 | `tox_seed/tox_viability_score/tox_viability_flag/imm_high_flag/imm_flag`(及原 `imm_*` 明细列) | 选 | toxicity | rank 软惩罚、报告 |
| 排名 | `score_thermo/final_score/final_rank` | 条件必(输出) | rank | 交付 |

> 读取兼容：t7 旧列 `seq_pass` 在统一层用**别名**，不做硬改名（避免改 t7 输出约定）；`guide_checked` 一律用于后续所有计算，取代各自读 `guide_mut_19/guide_antisense_19` 的分叉。

## 2.4 待裁定项（不擅自决定，先记录）

1. 每窗口条数：**已定稿 16 条/窗口（1 WT + 15 单点突变体）**；文档旧写“17 条”为历史写法，计划/简介已同步更正（2026-09-10）。另：**双模型一致性检验计划已放弃**（任务18 取消），效率辅助分只用 OligoFormer。
2. seed 定义：toxicity `g2–g7`(6mer,复刻 OligoFormer) vs thermo/任务6注释 `g2–g8` → 排序特征与脱靶用哪个区间需负责人确认；默认分别保留 `seed6`/`seed7` 两列由 config 选。
3. passenger 链：合成双链中突变体 sense 是 `target`（与 guide 有一处错配）还是 `rc(guide_mut)`（完全互补）？影响热力学“双链ΔG”与脱靶语义 → 需确认后写死于 duplex.py。
4. “seed 区结合能”操作化定义（见 §2.5）是否接受。

## 2.5 四项排序特征的推荐操作化定义（新代码，不改原算法）

| 特征 | 来源模块 | 定义（在统一记录列上） | 说明 |
|---|---|---|---|
| MFE | structure | `feat_mfe = mfe`（自折叠，软/硬阈值分离） | 直接取列 |
| 末端 ΔΔG | structure | `feat_ddg_ends = delta_deltaG_ends` | 直接取列 |
| 双链 ΔG | thermo | `feat_dG_duplex = dG_total`（理想+错配修正+侧翼） | 直接取列 |
| seed 区结合能 | thermo(派生) | `feat_dG_seed = Σ dG_mismatch_by_position[1:8](g2–g8)`（seed7 口径；WT 天然=0）＋ `dG_core_canonical` 的 seed 段贡献（可选叠加，权重由 config 决定） | **新增派生输出**：在 runner 中基于已输出的 19 长逐位列求和，不改 `siRNAThermoCalculator` 内部 |

归一化：min‑max 或 z‑score（config 选，方向统一“分数越大越好”：ΔG/MFE 取负处理或取反向）；粗排分 = Σ w_i·norm_i；惩罚：脱靶命中/毒信号按 config 罚分；终分 = 粗排 − 惩罚；权重默认值取文档方向（可消融）并在 meta.json 记录。

---

*下一步：Phase 3（Pipeline 设计）→ Phase 4（目录与映射）→ Phase 5 按 common→generation→rule→structure→thermo→toxicity→ranking→offtarget 逐模块实现。*
