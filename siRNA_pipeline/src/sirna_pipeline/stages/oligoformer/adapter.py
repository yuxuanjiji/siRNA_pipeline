# -*- coding: utf-8 -*-
"""OligoFormer 效率预测适配层（第五步 · 效率辅助分，内置执行器模式）。

技术路线第五步把 OligoFormer 沉默效率输出作为综合排序的"辅助分"。
本适配层（配置 enabled=true 推荐）采用两种执行方式：
  * **内置执行器**（默认）：调用 `executor.py`（同目录），按官方 infer 口径
    （RNA-FM 嵌入 → loader 数据集 → Oligo 模型 → pred[:,1]×1.341）对统一记录逐条出效率分；
  * **外部执行器契约**（可选）：`external_cmd <cmd> --input <csv> --output <csv>`。

优雅降级（不阻断主管道）：
  * 模型权重 / OligoFormer .venv / RNA-FM 资产缺失 → status 报告具体原因，
    run() 写 oligo_available=0 并 skip；rank 侧自动 α=1、β=0，排序口径不变；
  * RNA-FM 资产（<repo>/RNA-FM）是真实推理的硬前置（官方模型必须 concat FM 嵌入），
    未放置前不会伪造任何效率分。

保留/删除取舍与开启步骤见：docs/design/oligoformer_keep_delete.md。
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from ...common import records
from ...common.stage_io import write_manifest

_PROJECT_ROOT = Path(__file__).resolve().parents[4]      # siRNA_pipeline/
DEFAULT_REPO = _PROJECT_ROOT / "OligoFormer部分"          # 原生推理仓库（.venv/model/scripts）
DEFAULT_EXECUTOR = Path(__file__).resolve().parent / "executor.py"

DEFAULT_CFG = {
    "enabled": False,             # 由 configs 覆盖（推荐 true；缺资产时自动 graceful skip）
    "external_cmd": None,         # 自定义外部执行器（可选；None=内置 executor.py）
    "python": None,               # None→OligoFormer .venv python（存在时），否则 sys.executable
    "repo_dir": None,             # None→工程根/OligoFormer部分
    "model_path": None,           # None→<repo>/model/best_model.pth
    "cds_fasta": None,            # 靶基因 CDS（executor 重建 57nt 上下文用）
    "rnafm_python": None,         # 跑 RNA-FM 的独立解释器（py3.8/3.10 conda 环境，fairseq 用）
    "rnafm_required": True,       # True=缺 RNA-FM 资产视为未就绪（graceful skip）
    "device": None,               # None=executor 自动（有 CUDA 版 torch+显卡则用 GPU，否则 CPU）
    "seed": 42,
}

_OLIGO_OUT_COLS = records.OLIGO_COLS


def _default_venv_python(repo: Path) -> str | None:
    for cand in (repo / ".venv" / "Scripts" / "python.exe",
                 repo / ".venv" / "bin" / "python"):
        if cand.exists():
            return str(cand)
    return None


def _resolve_defaults(cfg: dict) -> dict:
    repo = Path(str(cfg.get("repo_dir") or DEFAULT_REPO)).resolve()
    python = cfg.get("python") or _default_venv_python(repo) or sys.executable
    model = str(cfg.get("model_path")) if cfg.get("model_path") else \
        str(repo / "model" / "best_model.pth")
    return {
        "repo_dir": str(repo),
        "python": python,
        "model_path": model,
        "executor": str(DEFAULT_EXECUTOR),
    }


def _rnafm_assets_ok(repo: Path) -> bool:
    """RNA-FM 前置：repo/RNA-FM/redevelop + pretrained/extract_embedding.yml。"""
    if not (repo / "RNA-FM").exists():
        return False
    yml = repo / "RNA-FM" / "pretrained" / "extract_embedding.yml"
    return yml.exists() or (repo / "RNA-FM" / "redevelop").exists()


def _detection_status(cfg: dict) -> dict:
    """探测依赖并给出状态/原因（不抛异常）。"""
    d = _resolve_defaults(cfg)
    st = {"enabled": bool(cfg.get("enabled")),
          "mode": "external" if cfg.get("external_cmd") else "builtin",
          **d,
          "model_found": Path(d["model_path"]).exists(),
          "repo_exists": Path(d["repo_dir"]).exists(),
          "rnafm_ok": _rnafm_assets_ok(Path(d["repo_dir"])),
          "ok": False, "reason": ""}
    if not st["enabled"]:
        st["reason"] = "disabled（默认未启用 OligoFormer 效率辅助分）"
        return st
    if st["mode"] == "external":
        exe = cfg["external_cmd"]
        if isinstance(exe, (list, tuple)):
            exe = exe[0]
        if not (shutil.which(str(exe)) or Path(str(exe)).exists()):
            st["reason"] = f"external_cmd 不存在或不可执行：{exe}"
            return st
        st["ok"] = True
        st["reason"] = "ready (external)"
        return st
    # 内置 executor 模式
    if not st["repo_exists"]:
        st["reason"] = f"OligoFormer 仓库不存在：{d['repo_dir']}"
        return st
    if not st["model_found"]:
        st["reason"] = f"模型权重不存在：{d['model_path']}"
        return st
    if not Path(d["python"]).exists():
        st["reason"] = f"python 不存在：{d['python']}"
        return st
    if bool(cfg.get("rnafm_required", True)) and not st["rnafm_ok"]:
        st["reason"] = (f"RNA-FM 资产缺失（需 {Path(d['repo_dir']) / 'RNA-FM'}，含 "
                        "pretrained/extract_embedding.yml），效率推理需先放置 RNA-FM；"
                        "当前 graceful skip（α=1、排序口径不变）")
        return st
    rp = cfg.get("rnafm_python")
    if rp and not Path(str(rp)).exists():
        st["reason"] = f"rnafm_python 不存在：{rp}"
        return st
    st["reason"] = "ready (builtin)"
    st["ok"] = True
    return st


class OligoFormerAdapter:
    """OligoFormer 效率适配器：探测 + 批量执行（内置 executor / 外部契约 / 跳过）。"""

    def __init__(self, cfg: dict | None = None):
        self.cfg = {**DEFAULT_CFG, **(cfg or {})}
        self._d = _resolve_defaults(self.cfg)

    def status(self) -> dict:
        return _detection_status(self.cfg)

    # ------------------------------------------------------------------
    def _run_external(self, input_csv: Path, out_csv: Path, log) -> None:
        cmd_base = self.cfg["external_cmd"]
        python = self.cfg.get("python")
        if python:
            cmd = [python]
            cmd += [cmd_base] if isinstance(cmd_base, str) else list(cmd_base)
        elif isinstance(cmd_base, (list, tuple)):
            cmd = list(cmd_base)
        else:
            cmd = [cmd_base]
        cmd += ["--input", str(input_csv), "--output", str(out_csv)]
        log.info("oligoformer external: %s", " ".join(str(c) for c in cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(
                f"oligoformer external 失败（exit={proc.returncode}）："
                f"{(proc.stdout + proc.stderr)[-1500:]}")

    def _run_builtin(self, input_csv: Path, out_csv: Path, log) -> None:
        cmd = [self._d["python"], self._d["executor"],
               "--input", str(input_csv), "--output", str(out_csv),
               "--repo", self._d["repo_dir"], "--model", self._d["model_path"],
               "--cds", str(self.cfg.get("cds_fasta") or ""),
               "--seed", str(int(self.cfg.get("seed", 42))),
               "--work-dir", str(out_csv.parent)]
        if self.cfg.get("rnafm_python"):
            cmd += ["--rnafm-python", str(self.cfg["rnafm_python"])]
        if self.cfg.get("fm_mode"):
            cmd += ["--fm-mode", str(self.cfg["fm_mode"])]
        if self.cfg.get("rnafm_model"):
            cmd += ["--rnafm-model", str(self.cfg["rnafm_model"])]
        if self.cfg.get("device"):
            cmd += ["--device", str(self.cfg["device"])]
        log.info("oligoformer builtin: %s", " ".join(str(c) for c in cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(
                f"oligoformer executor 失败（exit={proc.returncode}）："
                f"{(proc.stdout + proc.stderr)[-2000:]}")

    # ------------------------------------------------------------------
    def run(self, input_csv: str | Path, out_csv: str | Path, log=None) -> dict:
        """执行（或跳过）效率预测，写出统一记录 CSV，返回状态 dict。"""
        import logging
        log = log or logging.getLogger("sirna.oligoformer")
        input_csv, out_csv = Path(input_csv), Path(out_csv)
        if not input_csv.exists():
            raise FileNotFoundError(f"输入记录不存在：{input_csv}")
        status = self.status()

        rows = records.read_records(input_csv)
        for r in rows:
            r["oligo_available"] = "0"

        if not status["ok"]:
            log.info("oligoformer skip：%s", status["reason"])
            records.write_records(out_csv, rows)
            write_manifest(out_csv.parent, {
                "stage": "efficiency_oligoformer", "rows": len(rows),
                "status": status, "mode": "skip",
            })
            return status

        os.makedirs(out_csv.parent, exist_ok=True)
        work = out_csv.parent / "_executor_result.csv"
        if status["mode"] == "external":
            self._run_external(input_csv, work, log)
        else:
            self._run_builtin(input_csv, work, log)
        ext = records.read_records(work)
        ext_by = {r["variant_id"]: r for r in ext if r.get("variant_id")}
        n_map = 0
        for r in rows:
            e = ext_by.get(r["variant_id"])
            if e is None:
                continue
            n_map += 1
            for col in _OLIGO_OUT_COLS:
                if col in e and e[col] not in ("", None):
                    r[col] = e[col]
            r["oligo_available"] = "1"
        if n_map != len(rows):
            raise ValueError(
                f"oligoformer 结果行数与输入不一致：{n_map}/{len(rows)}")
        records.write_records(out_csv, rows)
        write_manifest(out_csv.parent, {
            "stage": "efficiency_oligoformer", "rows": len(rows), "mapped": n_map,
            "status": status, "mode": status["mode"],
            "model": self._d["model_path"],
            "repo": self._d["repo_dir"],
            "seed": self.cfg.get("seed"),
        })
        return status


def run_oligoformer(
    input_csv: str | Path,
    out_dir: str | Path,
    cfg: dict | None = None,
    log=None,
) -> tuple[Path, dict]:
    """便捷入口：返回 (输出 CSV, 状态)。"""
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    adapter = OligoFormerAdapter(cfg)
    out_csv = out_dir / "candidates_oligoformer.csv"
    status = adapter.run(input_csv, out_csv, log)
    return out_csv, status


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="OligoFormer 效率适配：内置 executor（默认）/ 外部契约；缺失资产 graceful skip")
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--enabled", action="store_true")
    ap.add_argument("--external-cmd", type=str, default=None)
    ap.add_argument("--python", type=str, default=None)
    ap.add_argument("--repo-dir", type=str, default=None)
    ap.add_argument("--model-path", type=str, default=None)
    ap.add_argument("--cds-fasta", type=str, default=None)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.oligoformer")
    cfg = {"enabled": args.enabled, "external_cmd": args.external_cmd,
           "python": args.python, "repo_dir": args.repo_dir,
           "model_path": args.model_path, "cds_fasta": args.cds_fasta,
           "seed": args.seed}
    try:
        out, status = run_oligoformer(args.input, args.out_dir, cfg=cfg, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("oligoformer 失败：%s", e)
        return 3
    log.info("oligoformer 完成：%s（status=%s）", out, status["reason"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
