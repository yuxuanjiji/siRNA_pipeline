# -*- coding: utf-8 -*-
"""大赛标准结果文件导出（对照《附件5-代码提交要求》第四节）。

主运行入口一键产出 `results/results.csv`（UTF-8，无 BOM），列设计对应必填字段：
  候选编号 / 所属赛道 / 候选序列 / 关键预测指标 / 对应模型与运行版本 / 备注。
字段与内部统一记录列的映射说明见 README「大赛提交对照」表；赛道官方《提交模板》
发布后，只需改本模块 HEADER 映射即可对齐（其余管道不动）。
"""
from __future__ import annotations

import json
from pathlib import Path

from .. import __version__ as _PIPE_VERSION
from . import records

# 赛道标识（占位；正式赛道名以官网为准，可在 configs/pipeline.yaml summary.track 覆盖）
TRACK_DEFAULT = "siRNA（反义寡核苷酸）靶向序列设计 — SFRP1(NM_003012)"

# 结果文件列（顺序即表头顺序；值来源 = 内部统一记录列/派生）
HEADER = [
    "candidate_id",          # 候选编号 = variant_id（如 W0001_wt / W0001_g12:C>G）
    "window_id",             # 靶窗编号
    "kind",                  # wt / mut
    "track",                 # 所属赛道（占位常量）
    "cds_start", "cds_end",  # 靶位在 CDS 上的 1-based 区间
    "target_seq",            # mRNA 靶窗 19 nt（5'→3'）
    "guide_seq",             # 引导链/反义链 19 nt（5'→3'，用于合成）
    "passenger_seq",         # 乘客链（按完全互补 rc(guide) 生成，口径注记见 README）
    "mutation_site",         # 突变位点描述（mut 时形如 g12:C>G；wt 为空）
    "rules_pass",            # 四项规则过滤是否通过
    "structure_pass",        # 结构强过滤是否通过
    "mfe",                   # 引导链自折叠 MFE（kcal/mol）
    "delta_deltaG_ends",     # 双链末端 ΔΔG（链选择方向）
    "dG_duplex_total",       # 整体双链 ΔG（kcal/mol）
    "dG_seed_binding",       # seed 区结合能（kcal/mol，越负越强）
    "score_thermo",          # 主分（四热力特征归一化加权）
    "score_oligo",           # OligoFormer 效率辅助分（未启用时空）
    "penalty_total",         # 脱靶/毒性软惩罚
    "final_score",           # 综合终分
    "final_rank",            # 最终排名（1=最优）
    "tox_viability_score",   # seed 毒性 cell-viability（越低越毒）
    "tox_viability_flag",    # 毒性命中（0/1）
    "imm_high_flag",         # 免疫刺激高危 motif 命中（0/1）
    "offtarget_risk",        # 脱靶风险标记（1=有风险；缺列空）
    "pipeline_version",      # 管道版本
    "ranking_params",        # 排序参数（权重/方向/α/β/归一化/seed 口径，JSON）
    "notes",                 # 备注（结构告警/突变说明等）
]

# 参考化学修饰列（第六步 chemmod；给出 chem_csv 时追加在末尾）
CHEM_HEADER_EXTRA = [
    "guide_mod",             # 引导链逐位修饰记法（mX=2′-OMe、fX=2′-F、s=PS、5′-P 前缀）
    "sense_mod",             # 乘客链逐位修饰记法
    "mod_rule",              # 采用的修饰规则版本（如 esc19-v1）
    "mod_notes",             # 修饰依据/注记（JSON 数组串）
    "mod_seed_ome",          # 是否启用 seed(g2–g8) 全 2′-OMe 旋钮（0/1）
]

_COMP = {"A": "U", "U": "A", "C": "G", "G": "C"}


def _rc(seq: str) -> str:
    """RNA 反向互补（结果文件乘客链用；T 视为 U）。"""
    return "".join(_COMP.get(b, b) for b in reversed(seq.upper().replace("T", "U")))


def _mutation_site(r: dict) -> str:
    if str(r.get("kind", "")) != "mut":
        return ""
    site = str(r.get("site", "") or "").strip() or (
        "g" + str(r.get("position", "")) if str(r.get("position", "")).strip() else "")
    wt = str(r.get("wt_nt", "") or "")
    mt = str(r.get("mut_nt", "") or "")
    if site and wt and mt:
        return f"{site}:{wt}>{mt}"
    return site or (wt + ">" + mt if (wt or mt) else "")


def _notes(r: dict) -> str:
    parts = []
    for key in ("warnings", "reject_reason"):
        v = str(r.get(key, "") or "").strip()
        if v:
            parts.append(v)
    if str(r.get("imm_high_flag", "")) == "1":
        parts.append("免疫刺激高危 motif")
    if str(r.get("tox_viability_flag", "")) == "1":
        parts.append("seed 毒性命中")
    if str(r.get("structure_reliable", "")) != "1":
        parts.append("结构为近似结果(structure_reliable=0)")
    return "；".join(parts)


def _rank_params(meta: dict | None) -> str:
    meta = meta or {}
    try:
        return json.dumps({
            "weights": meta.get("weights"),
            "directions": meta.get("directions"),
            "alpha": meta.get("alpha"),
            "beta": meta.get("beta"),
            "norm": meta.get("norm"),
            "seed_kind": meta.get("seed_kind"),
            "penalties": meta.get("penalties"),
            "aux_available": (meta.get("aux") or {}).get("available"),
        }, ensure_ascii=False, sort_keys=True)
    except (TypeError, ValueError):
        return ""


def export_standard_results(
    rank_csv: str | Path,
    out_path: str | Path,
    top_n: int | None = None,
    track: str | None = None,
    rank_meta: dict | None = None,
    chem_csv: str | Path | None = None,
    mod_rule: str = "ESC-19",
    seed_ome: bool | None = None,
) -> Path:
    """把 07_ranking 排序表导出为 submission 用 results.csv（UTF-8）。

    top_n：只导出前 N 名（None=全部已排）；track：赛道标识（None=TRACK_DEFAULT）。
    chem_csv：给定第六步 chemmod 产物时，按 variant_id 追加**参考化学修饰**列
    （guide_mod/sense_mod/mod_rule/mod_notes/mod_seed_ome）。
    """
    out_path = Path(out_path)
    rows = [r for r in records.read_records(rank_csv)
            if str(r.get("final_rank", "")).strip() not in ("", "nan")]
    rows.sort(key=lambda r: int(float(r["final_rank"])))
    if top_n:
        rows = rows[: int(top_n)]

    chem_by: dict[str, dict] = {}
    if chem_csv:
        cp = Path(chem_csv)
        if cp.exists():
            for cr in records.read_records(cp):
                chem_by[str(cr.get("variant_id", ""))] = cr
    header = list(HEADER) + (list(CHEM_HEADER_EXTRA) if chem_csv else [])

    track_val = track or TRACK_DEFAULT
    out_rows = []
    for r in rows:
        out_rows.append([
            str(r.get("variant_id", "")),
            str(r.get("window_id", "")),
            str(r.get("kind", "")),
            track_val,
            str(r.get("cds_start", "")), str(r.get("cds_end", "")),
            str(r.get("target_mRNA_19", "")),
            str(r.get("guide_checked", "")),
            _rc(str(r.get("guide_checked", ""))),
            _mutation_site(r),
            str(r.get("rules_pass", "")),
            str(r.get("structure_pass", "")),
            str(r.get("mfe", "")),
            str(r.get("delta_deltaG_ends", "")),
            str(r.get("dG_total", "")),
            str(r.get("feat_dG_seed", "")),
            str(r.get("score_thermo", "")),
            str(r.get("score_oligo", "")),
            str(r.get("penalty_total", "")),
            str(r.get("final_score", "")),
            str(r.get("final_rank", "")),
            str(r.get("tox_viability_score", "")),
            str(r.get("tox_viability_flag", "")),
            str(r.get("imm_high_flag", "")),
            str(r.get("offtarget_risk", "")),
            f"sirna_pipeline@{_PIPE_VERSION}",
            _rank_params(rank_meta),
            _notes(r),
        ])
        if chem_csv:
            cr = chem_by.get(str(r.get("variant_id", "")), {})
            out_rows[-1] += [
                str(cr.get("chem_guide_mod", "")),
                str(cr.get("chem_sense_mod", "")),
                str(mod_rule or ""),
                str(cr.get("chem_notes", "") or ""),
                ("1" if seed_ome else "0") if seed_ome is not None else "",
            ]

    def _escape(v: str) -> str:
        """RFC4180 风格转义（字段内可能含逗号/引号，如 JSON 参数列）。"""
        if any(ch in v for ch in ',"\n\r'):
            return '"' + v.replace('"', '""') + '"'
        return v

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write(",".join(_escape(v) for v in header) + "\n")
        for row in out_rows:
            fh.write(",".join(_escape(v) for v in row) + "\n")
    return out_path
