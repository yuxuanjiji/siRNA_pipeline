# -*- coding: utf-8 -*-
"""综合排序模块（Stage 7 · rank）——新模块（原项目无实现）。

按项目计划文档“阶段四·综合排序”口径实现：
  1) 特征：MFE(自折叠)、末端 ΔΔG、双链 ΔG、seed 区结合能 —— 归一化 + 方向统一；
  2) 热力学综合分（主分，加权归一化粗排）；3) OligoFormer 效率输出（辅助分，可选/优雅降级）；
  4) 脱靶/毒性软惩罚；5) 终分 = α·主分 + β·辅助分 − 惩罚 → final_rank。

seed 特征口径（v2，修复原“seed 错配代价恒 0 退化”问题，见 docs/design/ranking_weighting.md）：
  * 本库突变位点 {g1,g12,g17,g18,g19} 全在 seed(g2–g8) 之外 → dG_mismatch_by_position 的 seed 位恒 0，
    旧定义使该特征对所有候选退化。故改为：
      feat_dG_seed = seed 区 canonical NN 堆叠结合能（负值，越负结合越强，逐窗口/逐序列可变）
                   + seed 位错配代价（当前库恒 0，为未来 seed 突变留口）；
  * canonical 堆叠直接用共享 NN 表 common.nn_tables.STACK_DG 对 guide_checked 的 seed 内部
    相邻二核苷酸求和（guide 为 5'→3' 上链），不改任何 legacy 算法。

设计声明（写入 meta.json）：
  * 本模块不触碰任何 legacy 算法，仅消费统一记录列；
  * 权重/α/β 为“建议默认值”（文档要求经任务16 消融确定），全部可配置；
  * 方向(direction)：归一化后“越大越好”。默认 feat_mfe=+1(自折叠越弱越好)、
    feat_ddg_ends=+1、feat_dG_duplex=-1(ΔG 越负结合越强)、feat_dG_seed=-1(seed 结合越强越好)；
  * 辅助分（score_oligo）：OligoFormer efficacy（越大越好，min-max 后并入）；
    模型不可用(oligo_available=0/列缺失)时 β 自动置 0、α 回补为 1，排序不受影响。

分层口径（v2，2026-09-12；见 docs/design/variant_layer.md）：
  * **窗口层**（选哪个 19nt 窗口）：Q_window = α·score_thermo + β·score_oligo，口径不变
    （四基准集 mean ρ=0.569 已验证）；
  * **变体层**（同一窗口内选哪个变体）：final = Q_window(锚点) × C_match − penalty，
    C_match 由错配位置/类型决定（λ 来自错配集标定，见 _VARIANT_LAMBDA_BANDS）；
  * 原因：窗口级分数**在变体粒度上不成立**——实测 ΔQ_guide 的符号对 Δlabel 无信息
    （置换 p=1.000）、Δthermo 与 Δoligo 方向相反（ρ=−0.33）、DL 占 ΔQ 方差 82.5%；
    而 C_match 在其本职（同窗口内变体排序）上有实测效度（组内 ρ=0.611 / LOFO 0.594）；
  * **向后兼容**：WT 行 C_match≡1 → final_score 与旧口径**逐位一致**；
    `variant_layer.enabled=False` 可整体回退旧口径。
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from ...common import duplex, records
from ...common.nn_tables import STACK_DG
from ...common.stage_io import write_manifest

# ---------------------------------------------------------------------------
# 变体层：C_match（错配位置 × 类型）参数
# λ 来源分三类，逐带标注（元数据写 lambda_source，详见 docs/design/mismatch_literature_basis.md）：
#   fitted               = 错配集 v2 合格集（Ohnishi 78 行 / 4 条 guide）单错配 × 配对 WT 拟合
#                          （文献先验收缩 PN=4；组内 ρ=0.611 / LOFO 0.594）
#   literature_prior     = 无数据（n=0），取文献先验
#   literature_prior_over_fit = 有弱拟合但与强外部证据冲突，改用文献先验并记录冲突
#
# 文献依据（本轮精读，引文见 mismatch_literature_basis.md）：
#   * g1：RISC 中 guide 第 1 位**本就不与靶配对**（KRNB）；Holen 2005 在 9/9 个 siRNA 中
#         证明 guide 第 1 位 C→U 摆动（mRNA 侧 G:C→U:G）**提高活性**；AGO2 测得 t1A 使亲和力 +1 kcal/mol
#         → λ(g1)=0（并列入可增益位点）
#   * seed(g2–g8)：Holen 第 2/3/6/7 位显著耐受；AGO2 多数 seed 错配对切割影响很小、部分反而加速
#         → 保留拟合值（n=24，最大样本），不覆盖
#   * cleavage(g9–g11)：Holen 单个第 10 位错配"dramatic loss"；AGO2 三联 t9–t11 错配使
#         k_cleave 下降 >500×、中央区插入几乎不可检出 → **最脆弱位点**。原拟合 n=4 得 0.0942
#         与外部证据冲突，改用文献先验（取本标定尺度的最大量级）
#   * g12：AGO2 t12 错配常可省、t12A 反而快 2.5×，但 KRNB 主张 9–13 必须完全配对
#         → 取"中等代价"折中，并做敏感性；下调 n=2 拟合的 0.3082
#   * mid(g13–g16)：保留拟合（n=9）；λ(g13)>λ(g12) 与 AGO2（t13 比 t12 更扰动切割）一致
#   * p3(g17–g19)：AGO2 报告 t15–t21/t17–t21 的单/双/三错配**提高**切割速率与细胞内敲低，
#         并建议"引入 g18–g21 错配可能增强敲低"；KRNB 白名单同为 {1,15,17,18,19,20,21}
#         → λ=0（并列入可增益位点）
#   n=0 的带（g1 / p3）λ 完全等于先验，meta 与文档中标注为"纯先验"。
# ---------------------------------------------------------------------------
_VARIANT_LAMBDA_BANDS = {
    "g1": 0.0000,
    "seed": 0.1786,
    "cleavage": 0.3082,
    "g12": 0.1500,
    "mid": 0.1685,
    "p3": 0.0000,
}
_VARIANT_LAMBDA_SOURCE = {
    "g1": "literature_prior(n=0; Holen 9/9 增益, KRNB 该位不配对)",
    "seed": "fitted(n=24)",
    "cleavage": "literature_prior_over_fit(n=4; Holen p10 剧烈失活, AGO2 t9-11 >500x)",
    "g12": "literature_prior_over_fit(n=2; AGO2 t12 可省/t12A 更快 vs KRNB 9-13 须配对)",
    "mid": "fitted(n=9)",
    "p3": "literature_prior(n=0; AGO2 t15-21 增益, KRNB 白名单)",
}
# 类型因子：改用文献的 transition/transversion 轴（AGO2 的 81 参数模型正是此轴 × 位置），
# 并保留 GU 摆动为独立类。AGO2 实测该轴在 t12 处**符号反转**（t6–t11 颠换更伤、t12–t19 颠换反而更快），
# 故分两段给因子。**幅度均为先验**（文献只给方向，不给可标定数值）。
_VARIANT_TYPE_FACTOR = {"WC": 0.0, "GU": 0.5, "transition": 0.85, "transversion": 1.0}
_VARIANT_TYPE_FACTOR_3P = {"transition": 1.0, "transversion": 0.7}   # 位置 ≥12 段
_VARIANT_GAIN_POSITIONS = (1, 8, 12, 17, 18, 19, 20, 21)   # 文献支持"错配可增益"的位点
_VARIANT_MULTI3_FACTOR = 0.5      # ≥3 错配近似失活（Holen：三重突变几乎耗尽活性）→ 额外折半
_VARIANT_ONLY_PRIOR = [1, 17, 18, 19]


def variant_band(pos: int) -> str:
    """错配位置(1-based) → λ 带名（与标定脚本 fit_cmatch_v2.py 分带一致）。"""
    if pos == 1:
        return "g1"
    if 2 <= pos <= 8:
        return "seed"
    if 9 <= pos <= 11:
        return "cleavage"
    if pos == 12:
        return "g12"
    if 13 <= pos <= 16:
        return "mid"
    return "p3"


def mismatch_axis_class(guide_nt: str, mrna_nt: str) -> str:
    """按文献口径分类：WC / GU / transition / transversion。

    transition = 靶碱基相对其互补碱基发生 A↔G 或 C↔U；否则 transversion。
    （AGO2 的逐位置切割模型即用此轴；GU 摆动单独成类。）
    """
    g, m = guide_nt.upper(), mrna_nt.upper()
    comp = {"A": "U", "U": "A", "C": "G", "G": "C"}
    if g == comp.get(m, ""):
        return "WC"
    if {g, m} == {"G", "U"}:
        return "GU"
    c = comp.get(g, "")
    return "transition" if (c, m) in (("A", "G"), ("G", "A"), ("C", "U"), ("U", "C")) \
        else "transversion"


def c_match(row: dict, vcfg: dict) -> tuple[float, list[dict]]:
    """C_match = Π_错配 (1 − λ_band × λ_type)；WT（无 position 列）恒为 1.0。

    读候选列：position(突变 guide 位点 1-based)、mut_nt(突变后碱基)、
    paired_mRNA_nt(该位点配对的 mRNA 碱基)、n_mismatch(错配总数，≥3 额外折半)。
    返回 (C_match, 逐错配明细)。λ 缺省的带按 0 处理（不惩罚，保守）。
    """
    pos = records.num(row.get("position"))
    if pos is None:
        return 1.0, []
    mut = str(row.get("mut_nt", "")).strip().upper()
    paired = str(row.get("paired_mRNA_nt", "")).strip().upper()
    cls = mismatch_axis_class(mut, paired) if (mut and paired) else "transversion"
    band = variant_band(int(pos))
    lam = float((vcfg.get("lambda_bands") or {}).get(band, 0.0))
    tmap = dict(vcfg.get("type_factor") or _VARIANT_TYPE_FACTOR)
    if int(pos) >= 12:
        tmap.update(vcfg.get("type_factor_3p") or _VARIANT_TYPE_FACTOR_3P)
    tf = float(tmap.get(cls, _VARIANT_TYPE_FACTOR["transversion"]))
    factor = max(0.0, 1.0 - lam * tf)
    n_mm = records.num(row.get("n_mismatch"))
    multi3 = bool(n_mm is not None and n_mm >= 3)
    if multi3:
        factor *= float(vcfg.get("multi_mismatch_factor", _VARIANT_MULTI3_FACTOR))
    return factor, [{"pos": int(pos), "band": band, "class": cls,
                     "lambda": lam, "type_factor": tf,
                     "multi_mismatch_ge3": multi3, "factor": factor}]

DEFAULT_RANK_CFG = {
    # 权重/α/β 依据四基准集网格实验（outputs/analysis/weighting_sweep.json，2026-09-12）：
    #   当前默认(1,1,1.2,0.8 @ α0.8/β0.2) mean ρ=0.287 / worst=0.096
    #   推荐(1,1,1.2,0 @ α0.4/β0.6)      mean ρ=0.569 / worst=0.323  ← 采用
    #   （纯 DL β=1.0 mean 0.550；dG_seed 权重置 0 后不再贡献排序，仅保留列供诊断）
    "features": [
        {"name": "feat_mfe", "source": "mfe", "weight": 1.0, "direction": 1},
        {"name": "feat_ddg_ends", "source": "delta_deltaG_ends", "weight": 1.0, "direction": 1},
        {"name": "feat_dG_duplex", "source": "dG_total", "weight": 1.2, "direction": -1},
        {"name": "feat_dG_seed", "source": "guide_checked", "weight": 0.0,
         "direction": -1, "kind": "seed_binding"},
    ],
    "norm": "minmax",
    "seed_kind": "seed7",
    "require_rules_pass": True,
    "require_structure_pass": True,
    "missing_feature_policy": "skip",       # skip=不参与排序(保留记录)
    # 主分/辅助分线性组合：final_score = alpha*主分 + beta*辅助分 − 惩罚。
    # 四基准集证据支持 DL 权重更高（β=0.6 最优、纯 DL β=1.0 次之）；
    # alpha+beta 建议=1（可不等，内部仍不归一）。辅助分不可用时 beta→0、alpha→1。
    "alpha": 0.4,
    "beta": 0.6,
    "aux": {                                 # OligoFormer 效率辅助分（可选）
        "enabled": True,                     # True=列存在且可用则启用；缺列/未跑自动降级
        "source": "oligo_efficacy",          # 统一记录列（越大越好；由 efficiency 适配层写入）
        "direction": 1,
    },
    "penalties": {                           # 单位与归一化主/辅助分同量纲[0,1]
        "tox_viability_flag": 0.05,          # 列命中(值为1=毒/风险)即扣
        "imm_high_flag": 0.10,
        "offtarget_risk": 0.10,              # 由 offtarget 适配层写入(1=有脱靶风险)；缺列自动忽略
    },
    # 位置感知错配增强（默认关闭，保持原行为不变）
    "enhancements": {
        "position_aware_mismatch": False,
        "position_weights": {"seed": 2.0, "center": 1.0, "edges": 0.8},
    },
    # ---- 变体层（窗口内排序）：final = (锚点 + 改善通道) × C_match − penalty ----
    # enabled=False 即完全回退旧口径（全部行用 α·thermo+β·oligo−penalty）。
    "variant_layer": {
        "enabled": True,
        "anchor": "window_wt",          # window_wt=该窗口 WT 行的 Q_window；缺 WT 时退化取窗口内最高
        "lambda_bands": dict(_VARIANT_LAMBDA_BANDS),
        "type_factor": dict(_VARIANT_TYPE_FACTOR),
        "prior_only_positions": list(_VARIANT_ONLY_PRIOR),
        "keep_base_score": True,        # 输出 final_score_base(旧口径分) 便于对照/审计
        # 改善通道：允许变体**超过 WT**。错配并非只降效率——文献精读给出明确位点：
        #   * Holen 2005（即本项目 Holen 数据原文）guide 第 1 位 C→U 摆动在 9/9 个 siRNA 中提升活性；
        #   * AGO2 高通量：t15–t21/t17–t21 错配与 t12A 使单周转切割速率上升，并建议"引入 g18–g21 错配"；
        #   * KRNB 白名单同为 {1,15,17,18,19,20,21}。
        # 实现：窗口内**冻结 DL 项**（其变体级 Δ 无信息、与热力反相关），只让热力项随突变移动，
        # 且**按位点门控**（gain_positions）：非门控位点不享受增益。
        #   q_eff = Q_window(锚点) + γ·α·(score_thermo(变体) − score_thermo(WT))   [仅门控位点]
        # γ=0 即关闭该通道。增益幅度由机制/文献先验驱动，未用数据标定
        # （错配集不含 g1/g2 错配，且标签为 WT 归一化保留率）。
        "mech_gain_weight": 1.0,
        "gain_positions": list(_VARIANT_GAIN_POSITIONS),
        "type_factor_3p": dict(_VARIANT_TYPE_FACTOR_3P),
        "multi_mismatch_factor": _VARIANT_MULTI3_FACTOR,
        "lambda_source": dict(_VARIANT_LAMBDA_SOURCE),
    },
    "top_n": None,                            # None=全部
}

# 附加排序过程列（作为 FEATURE 派生结果也写入记录，便于核查）
_FEAT_RAW = ("feat_mfe", "feat_ddg_ends", "feat_dG_duplex", "feat_dG_seed")


def merge_tables(base_csv: str | Path, extra_csvs: list[str | Path] | None = None) -> list[dict]:
    """把基线表(如 structure 通过表)与各阶段(03/04/06…)表按 (window_id, variant_id) 合并。

    extras 仅用于“按主键补列”，不把额外行带回（被基线过滤掉的候选保持缺席）。
    同名键取非空值；键冲突以更靠后的表为准。
    """
    rows = records.read_records(base_csv)
    idx = {(r["window_id"], r["variant_id"]): r for r in rows}
    for path in (extra_csvs or []):
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"阶段表不存在：{p}")
        for r in records.read_records(p):
            key = (r.get("window_id"), r.get("variant_id"))
            if key in idx:
                idx[key].update({k: v for k, v in r.items() if v not in ("", None)})
    merged = [idx[k] for k in [ (r["window_id"], r["variant_id"]) for r in rows ]]
    records.validate_base(merged, "rank merge")
    return merged


def _seed_mismatch_cost(row: dict, seed_kind: str) -> float:
    """seed 位错配代价：dG_mismatch_by_position 在 seed 位求和；列缺失/解析失败按 0。

    当前库突变位点 {g1,g12,g17,g18,g19} 均在 seed(g2–g8) 之外 → 该项恒为 0；
    保留叠加是为未来 seed 位突变留口（见 docs/design/ranking_weighting.md）。
    """
    raw = row.get("dG_mismatch_by_position", "")
    if not raw:
        return 0.0
    try:
        vals = json.loads(raw)
    except (TypeError, ValueError):
        return 0.0
    if not isinstance(vals, list) or len(vals) != 19:
        return 0.0
    pos = duplex.seed_positions(seed_kind)          # 1-based
    s = 0.0
    for p in pos:
        try:
            s += float(vals[p - 1])
        except (TypeError, ValueError):
            pass
    return s


def _seed_canonical_dg(row: dict, seed_kind: str) -> float | None:
    """seed 区结合能（canonical NN 堆叠）。

    以 guide_checked 为 5'→3' 上链，对 seed 内部相邻二核苷酸查共享 STACK_DG 求和
    （负值=结合越强；seed 完全互补时双链 seed 堆叠能等于单链 NN 求和）。
    guide 缺失/非 19nt → None（该候选不参与排序，维持“缺特征跳过”语义）。
    """
    g = str(row.get("guide_checked", ""))
    if len(g) != 19:
        return None
    pos = list(duplex.seed_positions(seed_kind))    # 1-based
    s = 0.0
    for a, b in zip(pos, pos[1:]):
        if b == a + 1:                              # 仅相邻 seed 位构成 NN 堆叠步
            s += STACK_DG.get(g[a - 1] + g[b - 1], 0.0)
    return s


def _seed_binding(row: dict, seed_kind: str) -> float | None:
    """feat_dG_seed = seed 区 canonical 堆叠结合能 + seed 位错配代价（叠加）。"""
    canon = _seed_canonical_dg(row, seed_kind)
    if canon is None:
        return None
    return canon + _seed_mismatch_cost(row, seed_kind)


def _position_aware_mismatch_score(row: dict, seed_kind: str) -> float | None:
    """位置感知错配增强：在默认关闭时不参与排序；启用时只读取已有 dG_mismatch_by_position。

    设计目标：保留项目原行为，额外提供一个更细粒度的逐位错配摘要，用于实验型筛选器。若列缺失
    或无法解析，则返回 None，确保不会影响正常运行。该增强在默认配置中关闭，不会改变现有输出。
    """
    raw = row.get("dG_mismatch_by_position", "")
    if not raw:
        return None
    try:
        vals = json.loads(raw)
    except (TypeError, ValueError):
        return None
    if not isinstance(vals, list) or len(vals) != 19:
        return None
    pos = list(duplex.seed_positions(seed_kind))
    seed_set = set(pos)
    weights = []
    for i in range(1, 20):
        if i in seed_set:
            w = 2.0
        elif i in (1, 2, 18, 19):
            w = 0.8
        elif 9 <= i <= 11:
            w = 1.0
        else:
            w = 1.2
        weights.append(w)
    total = 0.0
    wsum = 0.0
    for idx, val in enumerate(vals, start=1):
        try:
            v = abs(float(val))
        except (TypeError, ValueError):
            continue
        w = weights[idx - 1]
        total += w * v
        wsum += w
    if wsum == 0.0:
        return None
    return total / wsum


def _feature_raw(row: dict, fcfg: dict, seed_kind: str) -> float | None:
    if fcfg.get("kind") == "seed_binding":
        return _seed_binding(row, seed_kind)
    if fcfg.get("kind") == "mismatch_position":
        return _position_aware_mismatch_score(row, seed_kind)
    v = records.num(row.get(fcfg["source"]))
    return v


def _eligible(row: dict, cfg: dict) -> tuple[bool, str]:
    """候选进入排序的条件。返回 (ok, reason) reason='' 表示通过。"""
    if cfg.get("require_rules_pass") and str(row.get("rules_pass", "")) != "1":
        return False, "rules_pass!=1"
    if cfg.get("require_structure_pass"):
        sp = row.get("structure_pass", "")
        if sp == "0":
            return False, "structure_pass=0"
        if sp == "":                            # 未做结构检测：视为不可排序
            return False, "structure 缺失"
    return True, ""


def _flags(row: dict, penalties_cfg: dict) -> dict:
    return {k: (str(row.get(k, "")) == "1") for k in penalties_cfg}


def run_rank(
    input_csv: str | Path,
    extra_csvs: list[str | Path] | None = None,
    out_dir: str | Path | None = None,
    cfg: dict | None = None,
    log=None,
) -> tuple[Path, dict]:
    """执行综合排序。

    返回 (候选排序 CSV, meta dict)。未满足排序条件的行保留在 CSV（排名列空），原因入 meta。
    """
    import logging
    log = log or logging.getLogger("sirna.rank")
    cfg = {**DEFAULT_RANK_CFG, **(cfg or {})}
    out_dir = Path(out_dir) if out_dir else Path(input_csv).parent
    os.makedirs(out_dir, exist_ok=True)

    rows = merge_tables(input_csv, extra_csvs)
    log.info("rank: 合并后 %d 行", len(rows))

    # ---- 1) 资格与特征提取 ----
    feat_cfgs = list(cfg["features"])
    enh = cfg.get("enhancements") or {}
    if bool(enh.get("position_aware_mismatch", False)):
        feat_cfgs.append({
            "name": "feat_mismatch_positional",
            "source": "dG_mismatch_by_position",
            "weight": 0.4,
            "direction": 1,
            "kind": "mismatch_position",
        })
    # 配置叠加：允许 yaml 用 weights/directions 字典按特征名覆盖默认值（便于消融/调参）
    w_over = cfg.get("weights") or {}
    d_over = cfg.get("directions") or {}
    feat_cfgs = [
        dict(f, weight=w_over.get(f["name"], f.get("weight", 1.0)),
             direction=d_over.get(f["name"], f.get("direction", 1)))
        for f in feat_cfgs
    ]
    seed_kind = str(cfg.get("seed_kind", "seed7"))
    scored = []          # 参与排序的行(附加特征)
    excluded: dict[str, int] = {}
    for r in rows:
        ok, why = _eligible(r, cfg)
        if not ok:
            excluded[why] = excluded.get(why, 0) + 1
            continue
        feats = {}
        for fc in feat_cfgs:
            v = _feature_raw(r, fc, seed_kind)
            if v is None:
                break
            feats[fc["name"]] = v
        if len(feats) != len(feat_cfgs):
            excluded["特征缺失"] = excluded.get("特征缺失", 0) + 1
            continue
        # 记录派生原始特征列（FEATURE_COLS），便于核查
        for name, v in feats.items():
            r[name] = "%.10g" % v
        scored.append({"row": r, "feats": feats})

    # ---- 2) 归一化（min-max；方向已把“越大约好”统一）----
    if not scored:
        raise ValueError("无满足排序条件的行")
    stats = {}
    for fc in feat_cfgs:
        name = fc["name"]
        vals = []
        for s in scored:
            x = s["feats"][name]
            if fc.get("direction", 1) < 0:
                x = -x
            s["feats"][name] = x          # 就地转换为“越大越好”
            vals.append(x)
        lo, hi = min(vals), max(vals)
        span = hi - lo
        for s in scored:
            x = s["feats"][name]
            s["feats"][name] = 0.5 if span == 0 else (x - lo) / span
        stats[name] = {"min": lo, "max": hi, "span": span}

    # ---- 3) 辅助分可用性（OligoFormer efficacy；整列可用才启用，否则优雅降级）----
    aux_cfg = cfg.get("aux") or {}
    aux_col = str(aux_cfg.get("source", "oligo_efficacy"))
    alpha = float(cfg.get("alpha", 0.8))
    beta = float(cfg.get("beta", 0.2))
    avail = [(s, records.num(s["row"].get(aux_col))) for s in scored]
    avail = [(s, v) for s, v in avail if v is not None]
    aux_used = bool(aux_cfg.get("enabled", True)) and len(avail) == len(scored)
    if aux_used:
        vals = [v for _, v in avail]
        if aux_cfg.get("direction", 1) < 0:
            vals = [-v for v in vals]
            avail = [(s, -v) for s, v in avail]
        lo, hi = min(vals), max(vals)
        span = hi - lo
        for s, v in avail:
            s["aux_raw"] = v
            s["score_oligo"] = 0.5 if span == 0 else (v - lo) / span
        aux_stats = {"available": True, "column": aux_col,
                     "min": lo, "max": hi, "span": span}
    else:
        alpha, beta = 1.0, 0.0                      # 降级：主分权重回补，排序不受影响
        reason = ("aux.enabled=false" if not aux_cfg.get("enabled", True)
                  else "辅助列缺失/部分缺失")
        aux_stats = {"available": False, "column": aux_col, "reason": reason}
    log.info("rank: alpha=%.2f beta=%.2f aux=%s", alpha, beta, aux_stats)

    # ---- 4) 主分(α·score_thermo) + 辅助分(β·score_oligo) + 软惩罚 ----
    wsum = sum(float(fc["weight"]) for fc in feat_cfgs) or 1.0
    for s in scored:
        r = s["row"]
        raw = sum(float(fc["weight"]) * s["feats"][fc["name"]] for fc in feat_cfgs)
        score = raw / wsum                                  # 0..1
        aux_v = s.get("score_oligo", 0.0)
        combined = alpha * score + beta * aux_v
        fl = _flags(r, cfg["penalties"])
        pen = sum(cfg["penalties"].get(k, 0.0) for k, hit in fl.items() if hit)
        r["score_thermo"] = "%.10g" % score
        if aux_used:
            r["score_oligo"] = "%.10g" % aux_v
        r["penalty_total"] = "%.10g" % pen
        r["q_window"] = "%.10g" % combined            # 窗口层分数（未扣惩罚）
        r["final_score"] = "%.10g" % (combined - pen)
        if (cfg.get("variant_layer") or {}).get("keep_base_score", True):
            r["final_score_base"] = "%.10g" % (combined - pen)   # 旧口径分（变体层会覆盖 final_score）

    # ---- 4b) 变体层：窗口内改由 C_match + 改善通道排序；窗口间仍用窗口层分数 ----
    vcfg = {**DEFAULT_RANK_CFG["variant_layer"], **(cfg.get("variant_layer") or {})}
    is_wt = lambda r: (str(r.get("kind", "")) == "wt"
                       or str(r.get("variant_id", "")).endswith("_wt"))   # noqa: E731
    gamma = float(vcfg.get("mech_gain_weight", 1.0))
    variant_meta = {"enabled": bool(vcfg.get("enabled", True)),
                    "mode": "anchor_plus_thermo_gain_times_c_match",
                    "anchor": vcfg.get("anchor", "window_wt"),
                    "lambda_bands": vcfg.get("lambda_bands"),
                    "type_factor": vcfg.get("type_factor"),
                    "prior_only_positions": vcfg.get("prior_only_positions"),
                    "mech_gain_weight": gamma,
                    "improvement_channel": ("q_eff = 锚点 + γ·α·Δscore_thermo；"
                                            "窗口内 DL 项冻结在 WT（其为噪声），热力项可升可降"),
                    "channel_note": ("错配非纯损失：5′ 端错配可改善链选择不对称性从而提高效率。"
                                     "该通道由机制/文献先验驱动，未用数据标定"
                                     "（错配集不含 g1/g2 错配，且标签为 WT 归一化保留率）。")}
    if variant_meta["enabled"]:
        by_window: dict[str, list[dict]] = {}
        for s in scored:
            by_window.setdefault(str(s["row"].get("window_id", "")), []).append(s)
        anchors: dict[str, float] = {}         # 锚点 = Q_window
        anchor_thermo: dict[str, float] = {}   # 锚点行的 score_thermo（改善通道基线）
        wt_final: dict[str, float] = {}        # 该窗口 WT 的 final_score（用于"超过 WT"计数）
        fallback_windows: set[str] = set()
        n_fallback = 0
        for wid, group in by_window.items():
            wt_rows = [s for s in group if is_wt(s["row"])]
            if wt_rows:                        # 该窗口的 WT
                ref = wt_rows[0]
                wt_final[wid] = float(ref["row"]["final_score"])
            else:                              # 无 WT 行（被前序硬过滤）→ 退化取窗口内最高 Q
                ref = max(group, key=lambda s: float(s["row"]["q_window"]))
                fallback_windows.add(wid)
                n_fallback += 1
            anchors[wid] = float(ref["row"]["q_window"])
            anchor_thermo[wid] = float(ref["row"]["score_thermo"])
        n_variant = 0
        n_above_wt = 0
        n_gain_on = 0
        for s in scored:
            r = s["row"]
            wid = str(r.get("window_id", ""))
            if is_wt(r):                       # WT: C_match≡1、Δthermo≡0 → 与旧口径逐位一致
                r["c_match"] = "1"
                r["q_window_anchor"] = "%.10g" % anchors[wid]
                r["thermo_delta_vs_wt"] = "0"
                r["q_variant_effective"] = "%.10g" % anchors[wid]
                r["variant_layer"] = "wt"
                continue
            cm, detail = c_match(r, vcfg)
            anchor = anchors[wid]
            d_thermo = float(r["score_thermo"]) - anchor_thermo[wid]
            # 改善通道**按位点门控**：只有文献支持"错配可增益"的位点才允许热力改善计入
            pos_v = records.num(r.get("position"))
            gain_pos = tuple(vcfg.get("gain_positions", _VARIANT_GAIN_POSITIONS))
            gated = bool(pos_v is not None and int(pos_v) in gain_pos)
            q_eff = anchor + (gamma * alpha * d_thermo if gated else 0.0)
            if gated:
                n_gain_on += 1
            r["c_match"] = "%.10g" % cm
            r["q_window_anchor"] = "%.10g" % anchor
            r["thermo_delta_vs_wt"] = "%.10g" % d_thermo
            r["q_variant_effective"] = "%.10g" % q_eff
            r["variant_layer"] = "cmatch"
            if detail:
                d = detail[0]
                r["c_match_band"] = d["band"]
                r["c_match_class"] = d["class"]
                r["c_match_lambda"] = "%.10g" % d["lambda"]
                r["c_match_type_factor"] = "%.10g" % d["type_factor"]
            r["final_score"] = "%.10g" % (q_eff * cm - float(r["penalty_total"]))
            n_variant += 1
            # 只统计"相对该窗口真实 WT 行"的超越（无 WT 的退化锚点窗口不计）
            if wid in wt_final and float(r["final_score"]) > wt_final[wid]:
                n_above_wt += 1
        variant_meta.update({"n_variant_rows": n_variant,
                             "n_variant_above_own_wt": n_above_wt,
                             "n_gain_channel_enabled_rows": n_gain_on,
                             "gain_positions": list(vcfg.get("gain_positions",
                                                             _VARIANT_GAIN_POSITIONS)),
                             "lambda_bands": vcfg.get("lambda_bands"),
                             "lambda_source": dict(_VARIANT_LAMBDA_SOURCE),
                             "type_axis": ("WC/GU/transition/transversion（AGO2 轴）；"
                                           "位置≥12 段符号反转，幅度为先验"),
                             "multi_mismatch_ge3_factor": float(
                                 vcfg.get("multi_mismatch_factor", _VARIANT_MULTI3_FACTOR)),
                             "n_windows": len(by_window),
                             "n_windows_with_wt_anchor": len(by_window) - n_fallback,
                             "n_windows_fallback_anchor": n_fallback,
                             "calibration_ref": "outputs/analysis/mismatch/"
                                                "mismatch_cmatch_calibration_v2.json；"
                                                "文献依据 docs/design/mismatch_literature_basis.md",
                             "evidence": "C_match 组内 Spearman 0.611 / LOFO 0.594"
                                         "（4 组 74 对，Ohnishi v2 合格集）；"
                                         "λ/类型幅度按文献先验校正",
                             "note": "窗口间仍用 Q_window（四基准集已验证）；变体层 = C_match × 改善通道。"
                                     "enabled=False 可完全回退旧口径。"})
        log.info("rank: variant_layer on（变体 %d 行，其中 %d 行超过本窗口 WT；"
                 "改善通道生效 %d 行；%d 窗口有 WT 锚点，%d 窗口退化锚点；γ=%.2f）",
                 n_variant, n_above_wt, n_gain_on,
                 variant_meta["n_windows_with_wt_anchor"],
                 variant_meta["n_windows_fallback_anchor"], gamma)

    # 终分降序；并列按原记录顺序（稳定），保证确定性
    scored.sort(key=lambda s: float(s["row"]["final_score"]), reverse=True)
    for rank, s in enumerate(scored, start=1):
        s["row"]["final_rank"] = str(rank)

    out_csv = out_dir / "candidates_ranked.csv"
    records.write_records(out_csv, rows)

    n_eligible = len(scored)
    n_top = len(scored) if not cfg.get("top_n") else min(int(cfg["top_n"]), len(scored))
    meta = {
        "stage": "07_ranking",
        "input_csv": str(input_csv),
        "extra_csvs": [str(p) for p in (extra_csvs or [])],
        "rows": len(rows),
        "eligible": n_eligible, "excluded": excluded,
        "features_stats": stats,
        "weights": {fc["name"]: fc["weight"] for fc in feat_cfgs},
        "directions": {fc["name"]: fc.get("direction", 1) for fc in feat_cfgs},
        "norm": cfg["norm"], "seed_kind": cfg["seed_kind"],
        "alpha": alpha, "beta": beta, "aux": aux_stats,
        "penalties": cfg["penalties"],
        "enhancements": cfg.get("enhancements", {}),
        "variant_layer": variant_meta,
        "top_n": cfg.get("top_n"), "ranked": n_top,
        "note": ("加权方案与默认值建议见 docs/design/ranking_weighting.md；"
                 "分层口径（窗口层/变体层）见 docs/design/variant_layer.md。"),
    }
    write_manifest(out_dir, meta)
    log.info("rank: eligible=%d excluded=%s top=%d", n_eligible, excluded, n_top)
    return out_csv, meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage7 综合排序：四特征归一化加权粗排 + 脱靶/毒性软惩罚 -> final_rank")
    ap.add_argument("--input", type=Path, required=True, help="基线表(02 passed)")
    ap.add_argument("--extra", type=Path, action="append", default=[],
                    help="阶段表(03 structure/04 thermo/06 toxicity…)，可多次")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.rank")
    try:
        csv, meta = run_rank(args.input, args.extra, args.out_dir, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("rank 失败：%s", e)
        return 3
    log.info("rank 完成：%s（eligible=%d）", csv, meta["eligible"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
