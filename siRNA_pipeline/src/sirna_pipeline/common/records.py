# -*- coding: utf-8 -*-
"""统一候选记录（CandidateRecord）schema 与 CSV 读写约定。

原则（docs/design/phase2 §2.3）：
  * 列只增不删；全工程 CSV 用 UTF-8-sig；
  * key = (window_id, variant_id)；variant_id 为 "{window_id}_wt" 表示完全互补 WT；
  * 数值列允许空串 = 缺失（NaN），由 numeric 转换统一处理；
  * 阶段白名单校验：每阶段只允许写自己负责的追加列 + 全量已有列。

本模块只做“容器/校验”，不做任何算法。
"""
from __future__ import annotations

from collections import OrderedDict
from typing import Iterable

from .seqio import read_csv_rows, write_csv

# ---------------------------------------------------------------------------
# 字段分组（顺序即推荐列顺序）—— 只增不删
# ---------------------------------------------------------------------------
BASE_COLS = [
    "window_id", "variant_id", "kind",                 # 主键/类型(wt|mut)
    "cds_start", "cds_end", "nm_003012_start", "nm_003012_end",
    "target_mRNA_19", "guide_wt_19", "guide_checked", "sense_strand_19",
]
MUT_COLS = ["site", "position", "wt_nt", "mut_nt", "paired_mRNA_nt", "pair_type"]
RULES_COLS = [
    "gc_pct", "fail_gc_range", "fail_run6_gc", "fail_run5_same",
    "fail_palindrome", "hit_code", "hit_rule", "rules_pass",
]
STRUCTURE_COLS = [
    "mfe", "mfe_structure", "max_stem_length", "max_stem_gc", "stem_loop_count",
    "delta_G_5end", "delta_G_3end", "delta_deltaG_ends",
    "cofold_mfe", "cofold_structure", "internal_loops", "bulges", "binding_energy",
    "target_accessibility", "target_5end_access", "target_3end_access",
    "structure_reliable", "structure_pass", "reject_reason", "warnings",
]
THERMO_COLS = [
    "dG_total", "dH_total", "dS_total", "Tm",
    "dG_core_canonical", "dG_mismatch_total", "dG_flank_total",
    "dG_mismatch_base", "dG_mismatch_context", "dG_mismatch_positional_extra",
    "dG_mismatch_by_position",        # JSON 串（长度 19 列表）
    "n_mismatch",
]
FEATURE_COLS = [                       # 派生排序特征（归一化输入）
    "feat_mfe", "feat_ddg_ends", "feat_dG_duplex", "feat_dG_seed",
]
OFFTARGET_COLS = [
    "seed6", "offtarget_pita_score", "offtarget_pita_filter",
    "offtarget_targetscan_score", "offtarget_targetscan_filter",
    "offtarget_pass", "offtarget_available",
    "offtarget_risk",            # 排序软惩罚用：1=有脱靶风险（由 offtarget 各层写入，缺省空=忽略）
]
BLAST_OFFTARGET_COLS = [               # BLAST 近全长同源脱靶层（第三层，BLAST/blast_offtarget.py）
    "offtarget_blast_flag",            # 1=命中 tolerant 判据被判脱靶（非靶基因/反义/配对数>=17/错配<=1）
    "offtarget_blast_n_hits_ge16",     # 该候选全部 >=16nt 同源命中数（审计留痕，不作判据）
    "offtarget_blast_available",       # 0=未跑/不可用；1=已回填
]
TOXICITY_COLS = [
    "tox_seed", "tox_viability_score", "tox_viability_flag", "tox_viability_miss",
    "tox_broadcast",
    "imm_ugugu", "imm_guccuucaa",
    "imm_polyu_max_run_guide", "imm_polyu_max_run_passenger",
    "imm_polyu_high", "imm_polyu_warn",
    "imm_u_content_guide", "imm_u_content_passenger", "imm_u_content_warn",
    "imm_high_flag", "imm_warn_flag", "imm_flag", "imm_hits_detail",
]
OLIGO_COLS = [                       # OligoFormer 效率辅助分（efficiency 适配层；可选/优雅降级）
    "oligo_efficacy",                # 越大越好（pred[:,1]*1.341 口径，见 docs/design/oligoformer_keep_delete.md）
    "oligo_available",               # 0=未跑/不可用；1=已回填
]
RANK_COLS = [
    "score_thermo", "score_oligo", "penalty_total", "final_score", "final_rank",
]
VARIANT_LAYER_COLS = [                 # 变体层（窗口内排序 = 锚点 × C_match；见 docs/design/variant_layer.md）
    "q_window",                        # 窗口层分数 α·score_thermo + β·score_oligo（未扣惩罚）
    "q_window_anchor",                 # 该候选所用锚点（同窗口 WT 的 q_window；无 WT 时取窗口内最高）
    "c_match",                         # 匹配代价系数 ∈(0,1]；WT 恒 1
    "c_match_band",                    # 所用 λ 位置带（g1/seed/cleavage/g12/mid/p3）
    "c_match_class",                   # 错配类型类（WC/GU/PP/YY/PY）
    "c_match_lambda",                  # 所用 λ_position
    "c_match_type_factor",             # 所用 λ_type
    "variant_layer",                   # 口径标记：wt / cmatch（便于审计与回退对照）
    "thermo_delta_vs_wt",              # 变体相对本窗口锚点行的 score_thermo 变化（改善通道输入）
    "q_variant_effective",             # 锚点 + γ·α·Δthermo（改善通道后的有效窗口分）
    "final_score_base",                # 旧口径分（α·thermo+β·oligo−penalty），用于对照/回退
]
CHEM_COLS = [                          # 化学修饰建议（第六步 chemmod 引擎输出，叶阶段）
    "chem_guide_mod", "chem_sense_mod",
    "chem_guide_f_positions", "chem_guide_ome_positions", "chem_guide_ps_bonds",
    "chem_sense_f_positions", "chem_sense_ome_positions", "chem_sense_ps_bonds",
    "chem_5p_phosphate", "chem_notes",
]

# 全量列顺序（拼合）
CANONICAL_COLS: list[str] = (
    BASE_COLS + MUT_COLS + RULES_COLS + STRUCTURE_COLS
    + THERMO_COLS + FEATURE_COLS + OFFTARGET_COLS + BLAST_OFFTARGET_COLS
    + TOXICITY_COLS + OLIGO_COLS + RANK_COLS + VARIANT_LAYER_COLS + CHEM_COLS
)
CANONICAL_SET = frozenset(CANONICAL_COLS)

# 数值型列（空串视为缺失；其余一律字符串原样保留）
NUMERIC_COLS = frozenset({
    "cds_start", "cds_end", "nm_003012_start", "nm_003012_end", "position",
    "gc_pct", "fail_gc_range", "fail_run6_gc", "fail_run5_same",
    "fail_palindrome", "hit_code", "rules_pass",
    "mfe", "max_stem_length", "max_stem_gc", "stem_loop_count",
    "delta_G_5end", "delta_G_3end", "delta_deltaG_ends",
    "cofold_mfe", "internal_loops", "bulges", "binding_energy",
    "target_accessibility", "target_5end_access", "target_3end_access",
    "structure_reliable", "structure_pass",
    "dG_total", "dH_total", "dS_total", "Tm",
    "dG_core_canonical", "dG_mismatch_total", "dG_flank_total",
    "dG_mismatch_base", "dG_mismatch_context", "dG_mismatch_positional_extra",
    "n_mismatch",
    "feat_mfe", "feat_ddg_ends", "feat_dG_duplex", "feat_dG_seed",
    "offtarget_pita_score", "offtarget_pita_filter",
    "offtarget_targetscan_score", "offtarget_targetscan_filter",
    "offtarget_pass", "offtarget_available", "offtarget_risk",
    "offtarget_blast_flag", "offtarget_blast_n_hits_ge16",
    "offtarget_blast_available",
    "tox_viability_score", "tox_viability_flag", "tox_viability_miss",
    "tox_broadcast", "imm_ugugu", "imm_guccuucaa",
    "imm_polyu_max_run_guide", "imm_polyu_max_run_passenger",
    "imm_polyu_high", "imm_polyu_warn",
    "imm_u_content_guide", "imm_u_content_passenger", "imm_u_content_warn",
    "imm_high_flag", "imm_warn_flag", "imm_flag",
    "oligo_efficacy", "oligo_available",
    "score_thermo", "score_oligo", "penalty_total", "final_score", "final_rank",
    "chem_5p_phosphate",
})

# 旧列名别名：读取/统计旧产物时用（写入统一记录一律用新名）
ALIASES = {
    "seq_pass": "rules_pass",              # task7 旧输出列 -> 统一列
    "guide_antisense_19": "guide_checked",  # task5 WT 列名 -> 统一被筛序列列
}

REQUIRED_BASE = frozenset({
    "window_id", "variant_id", "kind", "cds_start", "target_mRNA_19",
    "guide_wt_19", "guide_checked",
})


def num(value) -> float | None:
    """空/None/非法 -> None；否则 float。"""
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def to_int01(value) -> str:
    """布尔/0/1 归一为 '0'/'1'（None -> ''）。"""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    s = str(value).strip()
    if s in ("", "nan"):
        return ""
    return "1" if (float(s) != 0.0) else "0"


def validate_base(rows: Iterable[dict], context: str = "") -> None:
    """基础必填断言：主键非空、wt/mut 结构、guide_checked 19nt 合法。"""
    from collections import Counter
    from .seqio import assert_rna
    keys = Counter()
    for r in rows:
        wid, vid = str(r.get("window_id", "")), str(r.get("variant_id", ""))
        if not wid or not vid:
            raise ValueError(f"[{context}] 存在空 window_id/variant_id 行")
        keys[(wid, vid)] += 1
        kind = r.get("kind")
        if kind not in ("wt", "mut"):
            raise ValueError(f"[{context}] kind 须为 wt/mut：{kind!r}")
        if kind == "wt" and vid != f"{wid}_wt":
            raise ValueError(f"[{context}] wt 行 variant_id 应为 {wid}_wt：{vid}")
        g = str(r.get("guide_checked", ""))
        assert_rna(g, f"[{context}] guide_checked")
        if len(g) != 19:
            raise ValueError(f"[{context}] guide_checked 长度 != 19：{vid}")
    dup = {k: v for k, v in keys.items() if v > 1}
    if dup:
        raise ValueError(f"[{context}] 主键重复：{dup}")


def write_records(path, rows, cols: list[str] | None = None) -> int:
    """写统一记录 CSV（UTF-8-sig）。cols 缺省=全量列序（缺失列自动补空）。"""
    fieldnames = cols or CANONICAL_COLS
    fixed = []
    for r in rows:
        fixed.append({c: (r.get(c) if r.get(c) is not None else "")
                      for c in fieldnames})
    return write_csv(path, fieldnames, fixed)


def read_records(path) -> list[dict]:
    """读统一记录 CSV（返回 OrderedDict 行）。"""
    return read_csv_rows(path)


def load_with_aliases(path, alias_map: dict[str, str] | None = None) -> list[dict]:
    """读取（可能为旧格式的）CSV，并做列名别名映射（就地不改原文件）。"""
    rows = read_csv_rows(path)
    amap = dict(ALIASES)
    if alias_map:
        amap.update(alias_map)
    out = []
    for r in rows:
        nr = OrderedDict((amap.get(k, k), v) for k, v in r.items())
        out.append(nr)
    return out
