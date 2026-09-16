# 第六步化学修饰 · 规则引擎与行业金标准依据（chemmod_rules）

> 配套技术路线第六步："根据序列生成化学修饰建议，依据行业金标准"。
> 实现：`siRNA_pipeline/stages/chemmod/rules.py`（规则引擎）＋ `runner.py`（阶段运行器），
> 输出 `08_chemmod/candidates_chemmod.csv`（及 `results/chemmod_top.csv`）。
> 本文件说明每条规则的依据、参数默认值与工程口径。

---

## 1. 总体逻辑：ESC（增强稳定化学）骨架

选型依据：临床获批/在研 siRNA 普遍采用 Alnylam 系"增强稳定化学"配比思路——
**2′-O-methyl（OMe）为主干 + 2′-F 定点增强 + 端部硫代磷酸酯（PS）+ 引导链 5′-磷酸**，
其逻辑是"稳定性和降免疫用 OMe、亲和力与活性兼容用 2′-F、抗核酸外切酶用末端 PS"。
本项目把 21 nt 临床配比按 **19 nt 主链**缩放为 ESC-19 默认骨架：

| 链 | 默认布点（1-based） |
|---|---|
| guide（引导链/反义链，进 RISC） | 全链 2′-OMe；偶数位 {2,4,…,18} 改 2′-F；g1、g19 保持 OMe；端部 PS 键 {1,2,18}；**5′-P 开** |
| passenger（有义链/乘客链） | 全链 2′-OMe；{7,9,11,13,15,17} 改 2′-F；端部 PS 键 {1,2,18}；5′-P 关 |

## 2. 叠加掩蔽规则（优先级 > 默认配比）

1. **免疫掩蔽（默认开）**：引导链/乘客链上 **≥3 连续 U**、`UGUGU`、`GUCCUUCAA`
   覆盖的 **U 碱基强制改 2′-OMe**（哪怕该位默认是偶数 2′-F）。
   依据：序列依赖的免疫刺激由这些 U-rich/富含 GU 的序列基序驱动（Judge 2005）；
   2′-OMe 是公认的 TLR7/8 拮抗/掩蔽修饰（Judge 2006；Robbins 2009），临床配比中
   "免疫相关位点优先 OMe"是标准做法。
   > 阈值为 ≥3 连续 U，避免把普通 UU 全部 OMe 而破坏 F/OMe 交替的经典配比。
2. **seed OMe 旋钮（默认关，`seed_ome: false`）**：开启后 guide g2–g8 全部 2′-OMe，
   以降低 miRNA 样 seed 脱靶（Birmingham 2006 seed-match 语义），代价是部分沉默效率。
   因本管道第五步已用脱靶软惩罚处理风险，默认不开、按需双保险。

## 3. 设计约束（与证据一致性）

- **修饰不改碱基** → seed 脱靶、靶标互补性、免疫基序序列本身不变；修饰只能缓解
  降解/免疫，不能修复脱靶——脱靶应在第五步被惩罚掉，不在本阶段补救（同任务19 论证卡 §1.4）；
- **guide 5′ 端保持可磷酸化/不挡装载**：g1 不进入 2′-F 集合、5′-P 单独标注为"建议"，
  与 RISC/AGO2 装载对 5′ 端的要求一致（5′-P 是装载前提；ESC 端部 OMe 不冲突）；
- **central 切割区**沿用 ESC 的 F/OMe 交替布点（临床 ESC 即覆盖中央区），不额外做空窗。

## 4. 记法（与旧占位 `项目搭建/化学修饰/_mod_notation.py` 一致）

`mX`=2′-OMe-X，`fX`=2′-F-X；token 间 `s`=PS 键；guide 行前缀 `5′-P`。
示例（引擎对实际 Top-1 生成，非示例序列）：`5′-P mC fU mU fG mU fC … s mC`。

## 5. 输出列（每条候选追加，叶阶段）

| 列 | 含义 |
|---|---|
| `chem_guide_mod` / `chem_sense_mod` | 逐位修饰记法（含 5′-P 前缀） |
| `chem_guide_f_positions` / `chem_guide_ome_positions` / `chem_guide_ps_bonds` | 结构化布点（JSON） |
| `chem_sense_*` | 乘客链同上 |
| `chem_5p_phosphate` | guide 5′-P 建议（恒 1） |
| `chem_notes` | 命中说明（免疫 motif/多U/seed 旋钮/口径注记） |

## 6. 工程口径与限制（务必知悉）

- **19 nt 主链口径**：输入为管道候选（guide_checked 19nt）；passenger 按 `rc(guide)` 完全互补
  生成（与"sense=target 的 1-错配"口径差异见 `docs/design/phase2 §2.4#3`，已在 notes 标注）；
- **21 nt 化（2 nt 悬垂、临床全长 ESC）属后续成药化迭代**（旧 `任务19_成药性论证卡.md` 的
  21 nt 占位骨架即该目标形态）：本引擎只给 19 nt 建议并在 notes 提示，避免 19/21 口径冲突；
- **旧目录零改动**：`项目搭建/化学修饰/*`（占位候选/记法/论证卡）原样保留，仅作为口径与
  论证素材的参考来源；
- 本阶段**不淘汰、不改序列**，仅输出建议，供合成下单与任务19 正式论证卡替换使用。

## 7. 主要参考文献（标注核实状态）

1. Judge AD, et al. Sequence-dependent stimulation of the innate immune response by
   synthetic siRNA. *Nat Biotechnol* 2005;23:457-462.（UGUGU/GUCCUUCAA、poly-U 免疫基序）[与毒性阶段任务10 清单同源]
2. Judge AD, et al. Design of noninflammatory synthetic siRNA mediating potent gene
   silencing in vivo. *Mol Ther* 2006;13:494-505.（2′-OMe 掩蔽免疫）[同任务10 清单]
3. Robbins M, et al. 2′-O-methyl-modified RNAs act as TLR7 antagonists. *Mol Ther*
   2009;17:883-892.（2′-OMe 免疫拮抗机制）[同任务10 清单]
4. Hu B, et al. Therapeutic siRNA: state of the art. *Signal Transduct Target Ther*
   2020;5:101.（ESC/修饰汇总综述）
5. Birmingham A, et al. 3′ UTR seed matches… RNAi off-targets. *Nat Methods*
   2006;3:199-204.（seed 脱靶语义，支撑 seed OMe 旋钮）
6. 任务19 成药性论证卡（`项目搭建/化学修饰/任务19_成药性论证卡.md`）——21 nt ESC 配比、
   递送与可制造性论证的既有口径（本引擎为其 19 nt 自动前置步骤）。

> 状态说明：上表 1-3 已在任务10 调研文档标注"已核实"，4-6 为项目内既有引用（同论证卡口径）。
> 正式对外引用前请按论证卡 §5 复核卷页/PMID。
