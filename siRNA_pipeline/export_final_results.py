#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
export_final_results.py —— 大赛「最终结果文件」提交层导出器

读取 pipeline 已生成的 outputs/results/results.csv（及其排序/修饰表），
按大赛《四、最终结果文件》要求，一键输出提交用：
    final_results/results.xlsx          （带格式）
    final_results/results.csv           （UTF-8 带 BOM，Excel 直接打开不乱码）
    final_results/results_summary.txt   （字段校验审计报告）

必填字段：候选编号 / 所属赛道 / 候选序列（引导链 5'→3'）/ 关键预测指标 /
          对应模型与运行版本 / 备注；序列型候选结构文件列统一标注「无（序列型）」。

接入 predict.py（末尾追加）：
    from export_final_results import export_final_results
    export_final_results(results_path="outputs/results/results.csv",
                         rank_path="outputs/results/rank_final.csv",
                         chemmod_path="outputs/results/chemmod_top.csv",
                         outdir="final_results", topn=50)

独立命令行：
    python export_final_results.py --results outputs/results/results.csv --outdir final_results
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# 字段别名表（覆盖 pipeline 现有列名；模糊匹配，大小写/下划线不敏感）
# ---------------------------------------------------------------------------
FIELD_ALIASES: Dict[str, List[str]] = {
    "candidate_id": [
        "candidate_id", "variant_id", "candidate", "cand_id", "id", "编号",
        "候选编号", "候选序号", "seq_id", "sequence_id",
    ],
    "guide_seq": [
        "guide_seq", "guide_sequence", "guide_checked", "guide", "guide_strand",
        "引导链", "引导链序列", "sense_sequence", "sense", "sequence", "seq",
        "候选序列", "sirna_sequence",
    ],
    "passenger_seq": [
        "passenger_seq", "passenger_sequence", "passenger", "passenger_strand",
        "过客链", "过客链序列", "antisense_sequence", "antisense",
    ],
    "window": [
        "window_id", "window", "窗口", "窗口编号", "target_region", "region",
    ],
    "mutation_site": [
        "mutation_site", "mutation_sites", "mismatch_sites", "mismatch",
        "mutation", "variants", "错配位点", "突变位点", "site",
    ],
    "kind": ["kind", "type", "类别", "类型"],
    "final_score": [
        "final_score", "score", "total_score", "综合得分", "综合排序分",
        "终分", "总分", "排序分", "rank_score", "composite_score", "final",
    ],
    "final_rank": ["final_rank", "rank", "排名", "最终排名"],
    "eff_score": [
        "score_oligo", "eff_score", "eff_pred", "oligoformer_score",
        "oligo_efficacy", "oligoformer", "效率分", "预测效率",
        "silencing_efficiency", "dl_score", "efficiency_prediction",
    ],
    "thermo_score": [
        "score_thermo", "thermo_score", "thermodynamic_score", "thermo",
        "热力学主分", "热力学", "tm_score", "thermodynamic",
    ],
    "gc": ["gc_pct", "gc", "gc_content", "gc含量", "GC含量"],
    "offtarget": [
        "offtarget_risk", "offtarget", "off_target", "offtarget_penalty",
        "脱靶", "脱靶罚分", "offtarget_score",
    ],
    "toxicity": [
        "tox_viability_flag", "toxicity", "tox", "toxicity_penalty",
        "毒性", "毒性罚分", "tox_score", "tox_viability_score",
    ],
    "immuno": [
        "imm_high_flag", "imm_flag", "immuno", "immunostimulatory", "immune",
        "免疫", "免疫刺激", "immune_motif", "immuno_flag",
    ],
    "percentile": [
        "percentile", "percentile_rank", "百分位", "排名百分位", "pct_rank",
    ],
    "model_version": [
        "pipeline_version", "model_version", "model", "模型", "模型版本",
        "version", "tool_version",
    ],
    "ranking_params": ["ranking_params", "rank_params", "排序参数", "meta"],
    "notes": ["notes", "note", "备注", "comment", "comments", "warnings"],
    "chemmod": [
        "guide_mod", "sense_mod", "mod_notes", "chemmod", "modification",
        "modifications", "修饰建议", "mods", "chemical_modification",
    ],
    "delivery": ["delivery", "delivery_strategy", "递送", "递送策略"],
    "status": ["status", "state", "状态"],
}

# 输出列（固定顺序）
OUTPUT_COLUMNS: List[str] = [
    "候选编号",
    "所属赛道",
    "候选序列（引导链5'→3'）",
    "过客链序列",
    "靶基因",
    "窗口/靶区",
    "错配位点",
    "关键预测指标-综合排序分",
    "关键预测指标-最终排名",
    "关键预测指标-OligoFormer效率预测",
    "关键预测指标-热力学主分",
    "关键预测指标-GC含量(%)",
    "关键预测指标-脱靶风险/罚分",
    "关键预测指标-毒性罚分",
    "关键预测指标-免疫刺激风险",
    "对应模型与运行版本",
    "排序参数",
    "备注（修饰/递送/状态）",
    "结构文件",
]

REQUIRED_FIELDS: List[str] = [
    "候选编号",
    "所属赛道",
    "候选序列（引导链5'→3'）",
    "关键预测指标-综合排序分",
    "对应模型与运行版本",
    "备注（修饰/递送/状态）",
]

DEFAULT_TRACK = "赛道二：AI基因编辑与核酸工具设计"
DEFAULT_TARGET = "SFRP1"
DEFAULT_TOPN = 50


# ---------------------------------------------------------------------------
def _norm(s: str) -> str:
    s = s.lower().strip()
    return re.sub(r"[\s_\-（）()\[\]【】]+", "", s)


def _build_index(cols: List[str]) -> Dict[str, str]:
    return {_norm(c): c for c in cols}


def pick_column(cols: List[str], aliases: List[str]) -> Optional[str]:
    idx = _build_index(cols)
    for alias in aliases:
        hit = idx.get(_norm(alias))
        if hit is not None:
            return hit
    for alias in aliases:
        na = _norm(alias)
        for key, orig in idx.items():
            if len(na) >= 3 and (na in key or key in na):
                return orig
    return None


def _read_csv_any(path: str) -> Optional[object]:
    import pandas as pd
    if not path or not os.path.isfile(path):
        return None
    for enc in ("utf-8-sig", "utf-8", "gbk"):
        try:
            return pd.read_csv(path, encoding=enc, dtype=str)
        except (UnicodeDecodeError, pd.errors.ParserError):
            continue
    raise ValueError(f"无法解析 CSV（utf-8/utf-8-sig/gbk 均失败）: {path}")


def _safe_float(x) -> Optional[float]:
    try:
        return float(str(x).replace(",", "").replace("%", "").strip())
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
def export_final_results(
    results_path: Optional[str] = None,
    rank_path: Optional[str] = None,
    chemmod_path: Optional[str] = None,
    outdir: str = "final_results",
    track: str = DEFAULT_TRACK,
    target_gene: str = DEFAULT_TARGET,
    topn: int = DEFAULT_TOPN,
    output_base: str = "results",
) -> Dict[str, str]:
    """生成大赛标准化最终结果文件。返回 {"xlsx","csv","summary"} 路径。"""
    import pandas as pd

    os.makedirs(outdir, exist_ok=True)
    topn = max(1, int(topn))
    summary: List[str] = []

    df = _read_csv_any(results_path or "")
    if df is None and rank_path:
        df = _read_csv_any(rank_path)
        summary.append(f"[INFO] 未找到 --results，已用 rank 文件作主表: {rank_path}")
    if df is None:
        df = pd.DataFrame()
        summary.append("[WARN] 未找到任何上游结果文件，已生成空结果表（仅表头）。")

    cols = list(df.columns)
    summary.append(f"[INFO] 上游主表列数={len(cols)}")

    def col(aliases: List[str]) -> Optional[str]:
        return pick_column(cols, aliases)

    c_cid = col(FIELD_ALIASES["candidate_id"])
    c_guide = col(FIELD_ALIASES["guide_seq"])
    c_pass = col(FIELD_ALIASES["passenger_seq"])
    c_win = col(FIELD_ALIASES["window"])
    c_mut = col(FIELD_ALIASES["mutation_site"])
    c_kind = col(FIELD_ALIASES["kind"])
    c_score = col(FIELD_ALIASES["final_score"])
    c_rank = col(FIELD_ALIASES["final_rank"])
    c_eff = col(FIELD_ALIASES["eff_score"])
    c_thermo = col(FIELD_ALIASES["thermo_score"])
    c_gc = col(FIELD_ALIASES["gc"])
    c_off = col(FIELD_ALIASES["offtarget"])
    c_tox = col(FIELD_ALIASES["toxicity"])
    c_imm = col(FIELD_ALIASES["immuno"])
    c_pct = col(FIELD_ALIASES["percentile"])
    c_model = col(FIELD_ALIASES["model_version"])
    c_params = col(FIELD_ALIASES["ranking_params"])
    c_notes = col(FIELD_ALIASES["notes"])
    c_chem = col(FIELD_ALIASES["chemmod"])

    # ---- 按候选编号从 rank 补全缺失列（GC 等） -------------------------------
    if df.shape[0] > 0:
        if c_cid is None:
            df["__cid"] = [f"C{i+1:05d}" for i in range(len(df))]
            c_cid = "__cid"
            summary.append("[WARN] 主表无候选编号列，已自动生成 C00001… 编号。")
        else:
            df["__cid"] = df[c_cid].astype(str)

        if rank_path and os.path.isfile(rank_path):
            dfr = _read_csv_any(rank_path)
            if dfr is not None and dfr.shape[0] > 0:
                rc_cid = pick_column(list(dfr.columns), FIELD_ALIASES["candidate_id"])
                if rc_cid:
                    add_cols = {}
                    for dst, aliases in [
                        ("__gc", FIELD_ALIASES["gc"]),
                        ("__pct", FIELD_ALIASES["percentile"]),
                        ("__oligo", FIELD_ALIASES["eff_score"]),
                        ("__off", FIELD_ALIASES["offtarget"]),
                        ("__tox", FIELD_ALIASES["toxicity"]),
                        ("__imm", FIELD_ALIASES["immuno"]),
                        ("__score", FIELD_ALIASES["final_score"]),
                    ]:
                        if (dst == "__score" and c_score is not None):
                            continue
                        src = pick_column(list(dfr.columns), aliases)
                        if src:
                            add_cols[dst] = src
                    if add_cols:
                        sub = dfr[[rc_cid] + list(add_cols.values())].rename(
                            columns={rc_cid: "__cid", **{v: k for k, v in add_cols.items()}})
                        df = df.merge(sub, on="__cid", how="left")
                        summary.append(
                            f"[INFO] 已从 rank 文件合并列: {sorted(set(add_cols.values()))}")
                # 更新已解析列引用
                cols = list(df.columns)
                c_gc = col(FIELD_ALIASES["gc"])
                c_pct = col(FIELD_ALIASES["percentile"])
                c_eff = col(FIELD_ALIASES["eff_score"])
                c_score = col(FIELD_ALIASES["final_score"])
                c_off = col(FIELD_ALIASES["offtarget"])
                c_tox = col(FIELD_ALIASES["toxicity"])
                c_imm = col(FIELD_ALIASES["immuno"])

        if chemmod_path and os.path.isfile(chemmod_path):
            dfc = _read_csv_any(chemmod_path)
            if dfc is not None and dfc.shape[0] > 0:
                cc_cid = pick_column(list(dfc.columns), FIELD_ALIASES["candidate_id"])
                cc_mod = pick_column(list(dfc.columns), FIELD_ALIASES["chemmod"])
                if cc_cid and cc_mod:
                    sub = dfc[[cc_cid, cc_mod]].rename(
                        columns={cc_cid: "__cid", cc_mod: "__chem_merge"})
                    df = df.merge(sub, on="__cid", how="left")
                    c_chem = "__chem_merge"
                    summary.append("[INFO] 已从 chemmod 文件按候选编号合并化学修饰信息。")

    # ---- Top-N 排序截取 ------------------------------------------------------
    n_before = df.shape[0]
    if c_score is not None and df.shape[0] > 0:
        df["__score_num"] = df[c_score].map(_safe_float)
        df = df.sort_values("__score_num", ascending=False, na_position="last")
        summary.append(f"[INFO] 按综合排序分降序排列（{c_score}）。")
    else:
        summary.append("[WARN] 未找到综合排序分列，按原表顺序截取 Top-N。")
    df_top = df.head(topn).copy()
    summary.append(f"[INFO] 上游共 {n_before} 条候选，截取 Top-{topn} 输出。")

    # ---- 组装输出 -------------------------------------------------------------
    n = len(df_top)
    out = pd.DataFrame(index=range(n))

    def row_get(colname: Optional[str], i: int, default=""):
        if colname is None:
            return default
        v = df_top.iloc[i][colname]
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return default
        s = str(v).strip()
        return default if s in ("", "nan", "None") else s

    out["候选编号"] = [row_get(c_cid, i) or f"C{i+1:05d}" for i in range(n)]
    out["所属赛道"] = track
    out["候选序列（引导链5'→3'）"] = [row_get(c_guide, i) for i in range(n)]
    out["过客链序列"] = [row_get(c_pass, i) for i in range(n)]
    out["靶基因"] = target_gene
    out["窗口/靶区"] = [row_get(c_win, i) for i in range(n)]
    out["错配位点"] = [
        (row_get(c_mut, i) or ("野生型" if row_get(c_kind, i) in ("wt", "野生型", "WT") else ""))
        for i in range(n)]
    out["关键预测指标-综合排序分"] = [row_get(c_score, i) for i in range(n)]
    out["关键预测指标-最终排名"] = [row_get(c_rank, i) for i in range(n)]
    out["关键预测指标-OligoFormer效率预测"] = [row_get(c_eff, i) for i in range(n)]
    out["关键预测指标-热力学主分"] = [row_get(c_thermo, i) for i in range(n)]
    out["关键预测指标-GC含量(%)"] = [row_get(c_gc, i) for i in range(n)]
    out["关键预测指标-脱靶风险/罚分"] = [row_get(c_off, i) for i in range(n)]
    out["关键预测指标-毒性罚分"] = [row_get(c_tox, i) for i in range(n)]

    imm_vals = [row_get(c_imm, i) for i in range(n)]
    out["关键预测指标-免疫刺激风险"] = [
        ("有" if v in ("1", "TRUE", "True", "有", "yes") else ("无" if v else "")) for v in imm_vals]

    model_cells = []
    for i in range(n):
        m = row_get(c_model, i)
        if not m:
            m = f"sirna_pipeline@{_safe_float('0') or 0}" if False else ""
        model_cells.append(m or "sirna_pipeline（见 ranking_params/manifest）")
    out["对应模型与运行版本"] = model_cells

    out["排序参数"] = [row_get(c_params, i) for i in range(n)]

    notes_parts = []
    for i in range(n):
        bits = []
        chem = row_get(c_chem, i)
        if chem:
            bits.append(f"修饰:{chem[:200]}")
        note = row_get(c_notes, i)
        if note:
            bits.append(note)
        notes_parts.append("；".join(bits))
    out["备注（修饰/递送/状态）"] = notes_parts

    out["结构文件"] = ["无（序列型候选）" for _ in range(n)]

    out = out.reindex(columns=OUTPUT_COLUMNS)

    # ---- 必填校验 ------------------------------------------------------------
    for req in REQUIRED_FIELDS:
        miss = out[req].astype(str).apply(lambda x: x.strip() == "").sum()
        summary.append(f"[CHECK] {req}: 缺失 {miss}/{n}")

    # ---- 写出 -----------------------------------------------------------------
    xlsx_path = os.path.join(outdir, f"{output_base}.xlsx")
    csv_path = os.path.join(outdir, f"{output_base}.csv")
    summary_path = os.path.join(outdir, f"{output_base}_summary.txt")

    _write_xlsx(out, xlsx_path, track=track, topn=topn)
    out.to_csv(csv_path, index=False, encoding="utf-8-sig", lineterminator="\n")
    summary.append(f"[DONE] 已输出 {xlsx_path}")
    summary.append(f"[DONE] 已输出 {csv_path}（UTF-8 带 BOM）")

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary) + "\n")

    return {"xlsx": xlsx_path, "csv": csv_path, "summary": summary_path}


def _write_xlsx(df, path: str, track: str, topn: int) -> None:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter

    wb = Workbook()
    ws = wb.active
    ws.title = "最终结果清单"
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="2F5597")
    side = Side(style="thin", color="000000")
    border = Border(left=side, right=side, top=side, bottom=side)

    for ci, h in enumerate(df.columns, start=1):
        cell = ws.cell(row=1, column=ci, value=h)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = Alignment(horizontal="center", vertical="center")
        cell.border = border
    for ri, row in enumerate(df.itertuples(index=False), start=2):
        for ci, val in enumerate(row, start=1):
            cell = ws.cell(row=ri, column=ci, value=val)
            cell.alignment = Alignment(vertical="center", wrap_text=True)
            cell.border = border
    note_row = len(df) + 3
    ws.cell(row=note_row, column=1,
            value=f"说明：所属赛道={track}；共 {len(df)} 条（Top-{topn}）；"
                  f"字段规范见 results_template.xlsx「字段规范」；序列型候选无结构文件。").font = \
        Font(italic=True, size=9, color="808080")
    ws.freeze_panes = "A2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(df.columns))}{len(df)+1}"
    for ci, h in enumerate(df.columns, start=1):
        sample = [str(h) * 2] + [str(x) for x in df[h].astype(str).head(20)]
        width = max(12, min(42, int(max(len(s) for s in sample) * 1.1)))
        ws.column_dimensions[get_column_letter(ci)].width = width
    wb.save(path)


def create_template(outdir: str = ".", track: str = DEFAULT_TRACK) -> str:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill

    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "results_template.xlsx")
    wb = Workbook()
    ws = wb.active
    ws.title = "清单示例"
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill("solid", fgColor="2F5597")
    for ci, h in enumerate(OUTPUT_COLUMNS, start=1):
        c = ws.cell(row=1, column=ci, value=h)
        c.font = header_font
        c.fill = header_fill
        c.alignment = Alignment(horizontal="center", vertical="center")
    example = [
        "W0655_g18:U>G", track, "UCUUCUUGUCGCCAUUUGC", "GCAAAUGGCGACAAGAAGA",
        "SFRP1", "W0655", "g18:U>G",
        "0.8523", "1", "0.9389", "0.6607", "47.4", "0", "0", "无",
        "sirna_pipeline@0.1.0", '{"alpha":0.4,"beta":0.6}',
        "修饰:5′-P mU s fC …；ESC-19 骨架；杂交偏强(mfe<-25.0)",
        "无（序列型候选）",
    ]
    for ci, v in enumerate(example, start=1):
        c = ws.cell(row=2, column=ci, value=v)
        c.alignment = Alignment(vertical="center", wrap_text=True)
    ws.cell(row=4, column=1, value="↑ 示例行，提交时删除；字段规范见「字段规范」表。").font = \
        Font(italic=True, size=9, color="808080")
    ws.freeze_panes = "A2"

    ws2 = wb.create_sheet("字段规范")
    spec = [
        ("字段名", "必填", "类型", "单位/取值", "规则说明"),
        ("候选编号", "是", "文本", "如 W0655_g18:U>G", "全库唯一；与 rank_final.csv/chemmod_top.csv 对齐"),
        ("所属赛道", "是", "文本", "赛道二：AI基因编辑与核酸工具设计", "与报名赛道一致"),
        ("候选序列（引导链5'→3'）", "是", "文本", "19 nt RNA", "大写 ACGU；引导链（反义链）用于合成"),
        ("过客链序列", "推荐", "文本", "19 nt RNA", "rc(guide)；未设计填 N/A"),
        ("靶基因", "是", "文本", "SFRP1", "HGNC 官方缩写"),
        ("窗口/靶区", "推荐", "文本", "W0001…", "窗口编号（1-based CDS 区间见 rank_final）"),
        ("错配位点", "推荐", "文本", "g18:U>G", "位点按引导链 5' 端 g1 起算；野生型填'野生型'"),
        ("关键预测指标-综合排序分", "是", "数值", "0–1", "终分 = α·热力学主分 + β·OligoFormer − 惩罚"),
        ("关键预测指标-最终排名", "推荐", "整数", "1=最优", "对应 final_rank"),
        ("关键预测指标-OligoFormer效率预测", "推荐", "数值", "0–1", "辅助效率分"),
        ("关键预测指标-热力学主分", "推荐", "数值", "归一化", "四热力特征归一化加权"),
        ("关键预测指标-GC含量(%)", "推荐", "数值", "0–100", "一位小数"),
        ("关键预测指标-脱靶风险/罚分", "推荐", "数值/标记", "0/1 或罚分", "PITA/TargetScan/BLAST 合并"),
        ("关键预测指标-毒性罚分", "推荐", "数值/标记", "0/1 或罚分", "seed 6mer 活力"),
        ("关键预测指标-免疫刺激风险", "推荐", "文本", "有/无", "UGUGU、GUCCUUCA、CUGAAUU motif"),
        ("对应模型与运行版本", "是", "文本", "sirna_pipeline@x.y.z", "含 OligoFormer/ViennaRNA 等第三方版本"),
        ("排序参数", "推荐", "文本", "JSON", "α/β/weights/penalties 等"),
        ("备注（修饰/递送/状态）", "是", "文本", "自由文本", "2′-OMe/2′-F/PS/5′-P；递送建议；状态"),
        ("结构文件", "按需", "文本", "无（序列型候选）", "涉及三维结构时附 PDB/CIF 并标注文件名"),
    ]
    for r, row in enumerate(spec, start=1):
        for c, v in enumerate(row, start=1):
            cell = ws2.cell(row=r, column=c, value=v)
            if r == 1:
                cell.font = Font(bold=True, color="FFFFFF")
                cell.fill = PatternFill("solid", fgColor="2F5597")
            cell.alignment = Alignment(vertical="center", wrap_text=True)
    ws2.column_dimensions["A"].width = 34
    ws2.column_dimensions["B"].width = 8
    ws2.column_dimensions["C"].width = 12
    ws2.column_dimensions["D"].width = 26
    ws2.column_dimensions["E"].width = 60
    ws2.freeze_panes = "A2"
    wb.save(path)
    return path


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="大赛最终结果文件导出器")
    p.add_argument("--results", default=None, help="outputs/results/results.csv")
    p.add_argument("--rank", default=None, help="outputs/results/rank_final.csv（补全 GC/百分位等）")
    p.add_argument("--chemmod", default=None, help="outputs/results/chemmod_top.csv（补全修饰）")
    p.add_argument("--outdir", default="final_results")
    p.add_argument("--track", default=DEFAULT_TRACK)
    p.add_argument("--target-gene", default=DEFAULT_TARGET)
    p.add_argument("--topn", type=int, default=DEFAULT_TOPN)
    p.add_argument("--output-base", default="results")
    p.add_argument("--make-template", action="store_true")
    args = p.parse_args(argv)

    if args.make_template:
        tpl = create_template(args.outdir, track=args.track)
        print(f"[OK] 已生成模板: {tpl}")
        return 0

    paths = export_final_results(
        results_path=args.results, rank_path=args.rank, chemmod_path=args.chemmod,
        outdir=args.outdir, track=args.track, target_gene=args.target_gene,
        topn=args.topn, output_base=args.output_base)
    with open(paths["summary"], "r", encoding="utf-8") as f:
        print(f.read())
    print(f"[OK] 输出目录: {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
