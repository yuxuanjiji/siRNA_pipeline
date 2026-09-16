# 第四阶段：目录重构方案与旧→新文件映射（Phase 4）

> 铁律：**旧目录（项目搭建/…、数据集/…）一律原地不动**，作为“原始代码不丢失”的权威备份；
> 新工程内只放：① 共享层新代码；② 旧模块核心文件的**只读拷贝**(legacy/，含原校验)供 runner import；
> ③ 新增 runner/配置/测试/文档。任何“旧→新”都不得改变算法行。

## 4.1 最终目录树（按本项目实际调整，非机械套用）

```
siRNA_pipeline/
├── README.md                      # 快速开始/依赖/测试/运行
├── requirements.txt               # 运行时依赖(numpy, torch 可选, 其余标准库)
├── pyproject.toml                 # 包化: [tool.setuptools] packages
├── configs/
│   ├── paths.yaml  pipeline.yaml  ranking.yaml
│   └── stages/ generation.yaml rules.yaml structure.yaml thermo.yaml toxicity.yaml offtarget.yaml
├── src/
│   └── sirna_pipeline/
│       ├── __init__.py
│       ├── common/                # 共享层(新)
│       │   ├── seqio.py           # to_rna/rc_rna/gc_pct/normalize/pair_type/load_fasta/sliding_windows/CSV(BOM)
│       │   ├── duplex.py          # 反平行配对/57nt 重建(mrna57_from_cds)/seed6/seed7 常量
│       │   ├── nn_tables.py       # Xia1998 ΔG/ΔH 等唯一参数源(供 structure/thermo 引用)
│       │   ├── mismatch_class.py  # WC/GU/PP/YY/PY 分类(统一 calibrate/validate/thermo 语义)
│       │   ├── records.py         # CandidateRecord 列定义/白名单/校验器/别名表
│       │   └── stage_io.py        # 阶段读写(manifest/统计/路径解析/config 载入最小实现)
│       ├── stages/
│       │   ├── generation/  __init__.py runner.py  legacy/{seq_utils.py, task5_…py, task6_…py}
│       │   ├── rule_filter/  __init__.py runner.py legacy/{task7_…py}
│       │   ├── structure/    __init__.py runner.py legacy/{structure_detector.py}
│       │   ├── thermo/       __init__.py runner.py legacy/{thermo_calculator.py, calibrate_mismatch.py, validate_mismatch_dataset.py} legacy_mismatch_doubao/…
│       │   ├── toxicity/     __init__.py runner.py legacy/{toxicity_detector.py} resources/{cell_viability.txt(拷贝)}
│       │   ├── offtarget/    __init__.py adapter.py   # 接口+可用性探测; 工具链执行体留占位
│       │   └── rank/         __init__.py ranker.py normalize.py   # 新模块(文档算法, 权重配置化)
│       ├── pipeline/
│       │   ├── orchestrator.py   # 阶段调度/断点续跑/manifest
│       │   └── logging_setup.py
│       └── analysis/             # calibrate/validate 运行器(研究向, 引用 legacy)
├── scripts/
│   └── run_pipeline.py           # 唯一入口
├── tests/
│   ├── conftest.py
│   ├── test_common*.py  test_generation*.py  test_rules*.py
│   ├── test_structure*.py  test_thermo*.py  test_toxicity*.py
│   ├── test_offtarget_adapter.py  test_rank*.py  test_pipeline_e2e.py
├── data/processed/               # 由 分析 阶段产出的轻量派生(如 mismatch_validated 副本) + 说明
├── external/README.md            # 说明外部引用(../参考模型/../数据集)与可放入的工具链
├── docs/design/phase2..4 等
└── outputs/{runs,results,logs}
```

## 4.2 旧 → 新 映射表（核心代码；均以拷贝/引用方式迁移，算法不改）

| 旧路径（原位置不动） | 新位置（siRNA_pipeline/） | 迁移方式 | 说明 |
|---|---|---|---|
| `项目搭建\候选序列生成\seq_utils.py` | `src/sirna_pipeline/stages/generation/legacy/seq_utils.py` | copy(只读) | 供 legacy task5/6 import；共享收敛见 common/seqio.py |
| `…\task5_sliding_window_siRNA.py` | `…/generation/legacy/task5_sliding_window_siRNA.py` | copy | 原 CLI/校验完整保留，含 `--fasta/--out` |
| `…\task6_site_mutagenesis.py` | `…/generation/legacy/task6_site_mutagenesis.py` | copy | 依赖同目录 task5/seq_utils → 三文件必须同目录放置 |
| `项目搭建\初级规则筛选\task7_sequence_rules_filter.py` | `…/rule_filter/legacy/task7_sequence_rules_filter.py` | copy | 原实现+对拍校验保留 |
| `项目搭建\结构检测\structure_detector.py` | `…/structure/legacy/structure_detector.py` | copy | 类库原样；自检 `python legacy/structure_detector.py` 可跑 |
| `项目搭建\热力学参数计算\thermo_calculator.py` | `…/thermo/legacy/thermo_calculator.py` | copy | 唯一保留的热力学实现 |
| `…\thermo_calculator_mismatch(doubao).py` | `…/thermo/legacy_mismatch_doubao/thermo_calculator_mismatch_doubao.py` | copy+仅文件名去括号 | 归档研究版，不入主 runner（保留可回退） |
| `…\calibrate_mismatch.py` | `…/thermo/legacy/calibrate_mismatch.py`（引用） | copy | 研究用；新入口在 analysis 层，输出 data/processed |
| `…\validate_mismatch_dataset.py` | `…/thermo/legacy/validate_mismatch_dataset.py` | copy | 同上（其输入/输出路径在 analysis runner 中显式化） |
| `项目搭建\毒性检测\toxicity_detector.py` | `…/toxicity/legacy/toxicity_detector.py` | copy | 类库原样 |
| `项目搭建\毒性检测\data\cell_viability.txt` | `…/toxicity/resources/cell_viability.txt` | copy | 4096 行静态表；runner 经 config 定位（默认资源目录） |
| `项目搭建\化学修饰\_*.py` | 不迁移 | 原地保留 | 21nt 占位与主链(19nt)口径冲突，待任务19 重做；避免误导（已在 phase2 记录） |
| `项目搭建\脱靶检测\(空)` | `…/stages/offtarget/adapter.py` | 新建 | 接口先行 |
| 无(排序/排名) | `…/stages/rank/*` | 新建 | 按 phase2 §2.5 |
| `数据集\SFRP1-mRNA.txt` 等 | 不拷贝 | config 引用 | 单数据源原则，避免体积翻倍；`paths.yaml: data_root: ../数据集` |
| `参考模型\OligoFormer-main`、`siRNADiscovery-2` | 不拷贝 | config 引用 | external/README 说明版本与运行条件；duplicate 问题见 phase2 记录 |

## 4.3 无法确定归属、需保留原处并说明的文件

| 文件/目录 | 原因 |
|---|---|
| `项目搭建\结构检测\ViennaRNA-master`、`ViennaRNA-2.7.2`(≈350MB 源码) | 第三方源码树，无编译产物；新工程应通过 conda 安装 ViennaRNA 而非携带源码 → external/README 给出安装法；如需保留本地构建工作，留在旧目录即可 |
| `_pl3/_plfold_work/_pl_probe*/` 临时目录 | RNAplfold 运行垃圾(ACL 受限)；新工程改在 runs/ 内隔离工作目录，旧垃圾清理与否不阻塞 |
| `artifacts/`、`_research_src/`(35 脚本+PMC 原件) | 数据挖掘“研究现场”，与运行管道无关；建议后续单独归档为 research/，本轮不搬（避免误删/丢失溯源） |
| 根目录/数据集 多份 `*项目计划*.docx`、模块下 7 份结构 docx、热力学 2 份 docx | 历史文档版本；新文档统一进 `docs/`，旧版保留原处备查 |
| `项目搭建\Oligoformer部分\OligoFormer-main` 与 `参考模型\OligoFormer-main` | 双副本问题(已 MD5 验证同源)：本轮以 `参考模型/OligoFormer-main` 为引用基准并写入 external/README；是否删一份由负责人定 |
| `siRNA候选生成与筛选管道流程图.html`、`dot.ps`、根目录各申报书 | 对外材料，不属运行代码；保持原处 |

## 4.4 本阶段自检

- [ ] 每个被迁移核心文件在本轮内逐字节拷贝并 diff 校验（git 提交时展示 hash）
- [ ] legacy 文件不参与新包 import 链（只被 runner 显式 import 或 subprocess 调用）
- [ ] 新代码 import 方向唯一：runner → legacy；common ← runner；禁止 legacy → 新代码
