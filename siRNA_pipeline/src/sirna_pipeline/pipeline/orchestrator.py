# -*- coding: utf-8 -*-
"""Pipeline 编排器（把所有阶段串成一条命令）。

职责：解析配置 -> 顺序执行各阶段 runner -> 每阶段目录产出统一记录 + manifest
     -> 阶段间传递文件路径 -> 汇总 summary.json 与 results/rank_final.csv。

数据流：
  01 generation --rules--> 02 rules(仅通过子表) --structure--> 03 structure
  02 --thermo--> 04 ；02 --offtarget(默认跳过)--> 05；02 --toxicity--> 06
  rank 输入 = 03(structure 输出，含 pass) + extras(04/05/06)，按 variant_id join

错误策略：遇错即停（PipelineError）；summary 记录已成功阶段与失败原因。
"""
from __future__ import annotations

import copy
import json
import shutil
import time
from pathlib import Path

from ..common import records, stage_io
from ..common.seed import seed_everything
from ..common.results_export import export_standard_results
from ..stages.generation.runner import run_generation
from ..stages.rule_filter.runner import run_rules
from ..stages.structure.runner import run_structure
from ..stages.thermo.runner import run_thermo
from ..stages.toxicity.runner import run_toxicity
from ..stages.offtarget.adapter import run_offtarget
from ..stages.offtarget.blast_adapter import run_blast_offtarget
from ..stages.oligoformer.adapter import run_oligoformer
from ..stages.chemmod.runner import run_chemmod
from ..stages.rank.ranker import run_rank

STAGE_ORDER = [
    ("generation", "01_generation"),
    ("rules", "02_rules"),
    ("structure", "03_structure"),
    ("thermo", "04_thermo"),
    ("offtarget", "05_offtarget"),
    ("toxicity", "06_toxicity"),
    ("rank", "07_ranking"),
    ("chemmod", "08_chemmod"),
]


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def build_config(config_dir: str | Path, run_name: str | None = None,
                 fasta_override: str | Path | None = None) -> dict:
    """装载并合并 configs/ 下 paths + pipeline + stages 配置。"""
    config_dir = Path(config_dir).resolve()
    paths = stage_io.load_config(config_dir / "paths.yaml")
    pipeline = stage_io.load_config(config_dir / "pipeline.yaml")
    stages_cfg = stage_io.load_config(config_dir / "stages.yaml")
    paths_r = stage_io.resolve_paths(paths, config_dir)
    stages_cfg.pop("stages", None)                 # 若顶层有 stages 段则并入
    pipeline["stages"] = _deep_merge(pipeline.get("stages", {}), stages_cfg)
    pipeline["paths"] = paths_r
    if run_name:
        pipeline["run"]["run_name"] = run_name
    if fasta_override:
        paths_r["cds_fasta"] = str(Path(fasta_override).resolve())
    return pipeline


class PipelineError(Exception):
    """管道执行失败（数据/环境类）。"""


def _stage_dir(cfg: dict, key: str) -> Path:
    outputs = Path(cfg["paths"]["outputs_root"])
    run_name = cfg["run"].get("run_name", "default")
    base = outputs / "runs" / run_name
    d = {k: n for k, n in STAGE_ORDER}[key]
    p = base / d
    p.mkdir(parents=True, exist_ok=True)
    return p


def run_pipeline(cfg: dict, log=None, stage_from: str | None = None,
                 stage_until: str | None = None) -> dict:
    """执行管道并返回 summary。"""
    import logging
    log = log or logging.getLogger("sirna.pipeline")
    outputs = Path(cfg["paths"]["outputs_root"])
    run_name = cfg["run"].get("run_name", "default")
    log_dir = outputs / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    exe = stage_io.setup_logging("sirna.pipeline.exec",
                                 log_dir / f"pipeline_{run_name}.log")
    exe.info("=== pipeline start run=%s ===", run_name)

    seed = int(cfg.get("run", {}).get("seed", 42))
    seed_everything(seed)
    exe.info("random seed = %d（可复现性，见 common/seed.py）", seed)

    fasta = Path(cfg["paths"]["cds_fasta"])
    if not fasta.exists():
        raise PipelineError(f"CDS FASTA 不存在：{fasta}")
    stages_cfg = cfg["stages"]

    keys = [k for k, _ in STAGE_ORDER]
    if stage_from and stage_from in keys:
        keys = keys[keys.index(stage_from):]
    if stage_until and stage_until in keys:
        keys = keys[:keys.index(stage_until) + 1]

    summary: dict = {"run_name": run_name, "seed": seed, "stages": {}, "timings": {}}
    links: dict[str, Path | None] = {}
    structure_src = None
    thermo_src = None

    # 阶段执行键 -> 配置段键（STAGE_ORDER 名称与 configs 段名可不同）
    _cfg_key = {"rank": "ranking"}      # 阶段名 rank；配置段名为 ranking

    for key in keys:
        scfg = stages_cfg.get(_cfg_key.get(key, key), {}) or {}
        # 脱靶阶段 = PITA/TargetScan 层 + BLAST 近全长层，任一启用即进入
        enabled_stage = bool(scfg.get("enabled", True))
        if key == "offtarget":
            enabled_stage = enabled_stage or bool((scfg.get("blast") or {}).get("enabled"))
        if not enabled_stage:
            exe.info("stage %s disabled，跳过", key)
            summary["stages"][key] = {"ok": True, "skipped": True}
            continue
        out_dir = _stage_dir(cfg, key)
        t0 = time.time()
        try:
            if key == "generation":
                out_csv = run_generation(fasta, out_dir, log=exe)
            elif key == "rules":
                full_csv, passed_csv = run_rules(_req(links, "generation", exe),
                                                 out_dir, log=exe)
                links["rules_full"] = full_csv
                out_csv = passed_csv
            elif key == "structure":
                src = links.get("structure_src") or _req(links, "rules", exe)
                out_csv = run_structure(
                    src, fasta, out_dir,
                    detector_kwargs=scfg.get("detector_kwargs"),
                    only_passed=bool(scfg.get("only_passed", True)), log=exe)
                # 结构硬过滤通过子表（供排序；含规则列）
                structure_src = _write_passed(out_dir, out_csv, "structure_pass", exe)
                if not records.read_records(structure_src):
                    exe.warning("structure 通过子表为空（近似/阈值下无通过者），排序将回退到 rules passed")
                    structure_src = None
                links["structure_src"] = structure_src
            elif key == "thermo":
                src = thermo_src or links.get("rules") or _req(links, "rules", exe)
                out_csv = run_thermo(src, fasta, out_dir,
                                     only_passed=True, log=exe)
            elif key == "offtarget":
                src = links.get("rules") or _req(links, "rules", exe)
                pita_on = bool(scfg.get("enabled"))
                blast_on = bool((scfg.get("blast") or {}).get("enabled"))
                if blast_on:
                    # BLAST 近全长层（第三层）：可选叠加 PITA/TargetScan 层后串接，
                    # 逐层在统一记录上追加列；BLAST 层 graceful skip 时结果无判据列。
                    cur = src
                    if pita_on:
                        pita_dir = out_dir / "_pita"
                        cur, st = run_offtarget(
                            cur, pita_dir,
                            cfg={"enabled": True,
                                 "external_cmd": scfg.get("external_cmd"),
                                 "python": scfg.get("python")}, log=exe)
                        summary["offtarget_status"] = st
                    bcfg = dict(scfg.get("blast") or {})
                    bcfg["enabled"] = True
                    cur, bst = run_blast_offtarget(cur, out_dir, cfg=bcfg, log=exe)
                    summary["blast_offtarget_status"] = bst
                    out_csv = out_dir / "candidates_offtarget.csv"
                    shutil.copyfile(str(cur), out_csv)
                else:  # 仅 PITA/TargetScan 层（原行为）
                    out_csv, st = run_offtarget(
                        src, out_dir,
                        cfg={"enabled": True,
                             "external_cmd": scfg.get("external_cmd"),
                             "python": scfg.get("python")}, log=exe)
                    summary["offtarget_status"] = st
            elif key == "toxicity":
                src = links.get("rules") or _req(links, "rules", exe)
                out_csv = run_toxicity(src, out_dir, only_passed=True,
                                       viability_threshold=float(
                                           scfg.get("threshold", 50.0)), log=exe)
            elif key == "rank":
                base = structure_src or _req(links, "rules", exe)
                extras = [p for p in (links.get("thermo"),
                                      links.get("toxicity")) if p]
                ot_cfg = stages_cfg.get("offtarget", {}) or {}
                if links.get("offtarget") and (
                        ot_cfg.get("enabled") or bool((ot_cfg.get("blast") or {}).get("enabled"))):
                    extras.append(links["offtarget"])
                # OligoFormer 效率辅助分（可选）：启用时跑适配层并把结果并入 rank extras；
                # 未启用 → rank 内自动优雅降级（β=0、α=1），不影响排序。
                og_cfg = stages_cfg.get("oligoformer", {}) or {}
                if og_cfg.get("enabled"):
                    og_src = links.get("rules") or _req(links, "rules", exe)
                    og_dir = out_dir.parent / "06_oligoformer"
                    og_run = dict(og_cfg)
                    og_run["enabled"] = True
                    og_run.setdefault("cds_fasta", cfg["paths"].get("cds_fasta"))
                    og_run.setdefault("seed", cfg.get("run", {}).get("seed", 42))
                    og_csv, og_st = run_oligoformer(og_src, og_dir, cfg=og_run,
                                                    log=exe)
                    extras.append(og_csv)
                    summary["oligoformer_status"] = og_st
                    links["oligoformer"] = og_csv
                rank_cfg = {k: v for k, v in scfg.items() if k != "enabled"}
                out_csv, meta = run_rank(base, extras or None, out_dir,
                                         cfg=rank_cfg or None, log=exe)
                summary["rank_meta"] = meta
                _publish_results(cfg, out_csv, meta, run_name, exe)
                _export_submission(cfg, out_csv, meta, run_name, exe)
            elif key == "chemmod":
                rank_csv = links.get("rank") or _req(links, "rank", exe)
                cm_top = scfg.get("top_n")
                if cm_top is None:
                    cm_top = cfg.get("summary", {}).get("top_n")
                cm_csv = run_chemmod(rank_csv, out_dir, top_n=cm_top,
                                     seed_ome=bool(scfg.get("seed_ome", False)),
                                     log=exe)
                _publish_chemmod(cfg, cm_csv, run_name, exe)
                summary["chemmod_csv"] = str(cm_csv)
                # 用 chemmod 产物重导标准结果文件 → results.csv 内含参考化学修饰列
                _export_submission(cfg, rank_csv,
                                   summary.get("rank_meta") or {}, run_name, exe,
                                   chem_csv=cm_csv)
            links[key] = out_csv
            summary["stages"][key] = {"ok": True, "out": str(out_csv)}
        except Exception as e:  # noqa: BLE001
            exe.error("stage %s 失败：%s", key, e)
            summary["stages"][key] = {"ok": False, "error": str(e)}
            summary["timings"][key] = round(time.time() - t0, 3)
            _write_summary(cfg, summary, run_name)
            raise PipelineError(f"stage {key} 失败：{e}") from e
        summary["timings"][key] = round(time.time() - t0, 3)
        exe.info("stage %s 完成（%.2fs）", key, summary["timings"][key])

    exe.info("=== pipeline finished ===")
    _write_summary(cfg, summary, run_name)
    return summary


def _req(links: dict, key: str, log) -> Path:
    p = links.get(key)
    if p is None:
        raise PipelineError(f"缺少前置阶段输出：{key}")
    return Path(p)


def _write_passed(out_dir: Path, src_csv: Path, flag_col: str, log) -> Path:
    """从阶段输出抽 flag_col=='1' 子表（供下一阶段/排序使用）。"""
    rows = records.read_records(src_csv)
    kept = [r for r in rows if r.get(flag_col) == "1"]
    dst = out_dir / "candidates_passed.csv"
    records.write_records(dst, kept)
    log.info("pass 子表 %s=%d/%d -> %s", flag_col, len(kept), len(rows), dst)
    return dst


def _publish_results(cfg: dict, rank_csv: Path, meta: dict, run_name: str,
                     log) -> None:
    results_dir = Path(cfg["paths"]["outputs_root"]) / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(rank_csv, results_dir / "rank_final.csv")
    top_n = int(meta.get("top_n") or 0)
    if top_n > 0:
        # 修复：rank_csv 行序为输入(窗口)序，取 TopN 前必须先按 final_rank 排序
        rows = sorted((r for r in records.read_records(rank_csv) if r.get("final_rank")),
                      key=lambda r: int(r["final_rank"]))
        records.write_records(results_dir / "rank_top.csv", rows[:top_n])
    log.info("results 已发布 -> %s", results_dir)


def _publish_chemmod(cfg: dict, chem_csv: Path, run_name: str, log) -> None:
    """把化学修饰建议叶产物发布到 outputs/results/chemmod_top.csv。"""
    results_dir = Path(cfg["paths"]["outputs_root"]) / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(chem_csv, results_dir / "chemmod_top.csv")
    log.info("chemmod 已发布 -> %s", results_dir)


def _export_submission(cfg: dict, rank_csv: Path, meta: dict, run_name: str,
                       log, chem_csv: Path | None = None) -> Path:
    """大赛标准结果文件导出（《附件5》：results/results.csv，UTF-8）。

    覆盖范围 = 排序 Top-N（meta.top_n → summary.top_n → 全部已排）。
    传入 chem_csv（08_chemmod 产物）时，追加**参考化学修饰**列
    （guide_mod/sense_mod/mod_rule/mod_notes/mod_seed_ome）。
    """
    from ..common.results_export import TRACK_DEFAULT
    results_dir = Path(cfg["paths"]["outputs_root"]) / "results"
    top_n = meta.get("top_n") or cfg.get("summary", {}).get("top_n")
    track = cfg.get("summary", {}).get("track") or TRACK_DEFAULT
    mod_rule, seed_ome = "ESC-19", None
    if chem_csv is not None:
        mf = stage_io.read_manifest(Path(chem_csv).parent) or {}
        mod_rule = mf.get("rule_version") or mod_rule
        seed_ome = bool(mf.get("seed_ome", False))
    out = export_standard_results(
        rank_csv, results_dir / "results.csv",
        top_n=int(top_n) if top_n else None,
        track=track, rank_meta=meta,
        chem_csv=chem_csv, mod_rule=mod_rule, seed_ome=seed_ome)
    log.info("标准结果文件已导出 -> %s（Top-N=%s，化学修饰列=%s）",
             out, top_n, "含" if chem_csv else "无")
    return out


def _write_summary(cfg: dict, summary: dict, run_name: str) -> Path:
    results_dir = Path(cfg["paths"]["outputs_root"]) / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    p = results_dir / f"run_{run_name}_summary.json"
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=2)
    return p
