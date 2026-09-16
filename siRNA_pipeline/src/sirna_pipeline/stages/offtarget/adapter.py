# -*- coding: utf-8 -*-
"""脱靶检测适配层（Stage 5 · offtarget）——新模块。

背景：任务9 原计划用 OligoFormer 内置 PITA + TargetScan 做 seed 脱靶预测，
但工具链为 Perl/bash + 仓库内脚本 + ViennaRNA 依赖，属于强外部环境依赖。
本适配层：
  1) 默认 disabled：输出整列留空 + offtarget_available=0（排序软惩罚自动忽略），
     保证管道在没有工具链的机器上也能完整跑通（graceful skip）；
  2) enabled + external_cmd 时：调用“外部脱靶执行器”子进程，按契约读写 CSV；
  3) 预留 status() 探测：Perl/bash/ViennaRNA/脚本路径逐一报告，便于环境就绪后开启。

外部执行器契约（docs 同步说明）：
    <cmd> --input <统一记录CSV> --output <输出CSV>
    输出 CSV 必须含 variant_id 列，可含 offtarget_pita_score/pita_filter/
    targetscan_score/filter/offtarget_pass 中任意子集（其余留空）；每行对应输入行。
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
from pathlib import Path

from ...common import duplex, records
from ...common.stage_io import write_manifest

_OFFTARGET_OUT_COLS = records.OFFTARGET_COLS

DEFAULT_CFG = {
    "enabled": False,             # 未提供工具链前默认关闭
    "external_cmd": None,         # 外部执行器路径（enabled 时必填）
    "python": None,               # 运行外部执行器的解释器（若为 python 脚本）
    "pita_threshold": None,
    "targetscan_threshold": None,
}


def _detection_status(cfg: dict) -> dict:
    """探测外部依赖并给出状态/原因。不抛异常。"""
    st = {
        "enabled": bool(cfg.get("enabled")),
        "external_cmd": bool(cfg.get("external_cmd")),
        "bash": bool(shutil.which("bash")),
        "perl": bool(shutil.which("perl")),
        "python": cfg.get("python") or "sys",
        "ok": False, "reason": "",
    }
    if not st["enabled"]:
        st["reason"] = "disabled（默认：未配置脱靶工具链，跳过）"
        return st
    if not st["external_cmd"]:
        st["reason"] = "enabled 但未配置 external_cmd（外部执行器契约见模块 docstring）"
        return st
    # enabled：要求可执行文件存在
    exe = cfg["external_cmd"]
    if isinstance(exe, (list, tuple)):
        exe = exe[0]
    if not shutil.which(str(exe)) and not Path(str(exe)).exists():
        st["reason"] = f"external_cmd 不存在或不可执行：{exe}"
        return st
    st["reason"] = "ready"
    st["ok"] = True
    return st


class OffTargetAdapter:
    """脱靶适配器：探测 + 批量执行（子进程或跳过）。"""

    def __init__(self, cfg: dict | None = None):
        self.cfg = {**DEFAULT_CFG, **(cfg or {})}

    def status(self) -> dict:
        return _detection_status(self.cfg)

    # ------------------------------------------------------------------
    def _run_external(self, input_csv: Path, out_csv: Path, log) -> None:
        env = dict(os.environ)
        env["PYTHONIOENCODING"] = "utf-8"
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
        log.info("offtarget external: %s", " ".join(str(c) for c in cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", env=env, timeout=3600)
        if proc.returncode != 0:
            raise RuntimeError(
                f"offtarget external 失败（exit={proc.returncode}）：{(proc.stdout + proc.stderr)[-1500:]}")

    def run(self, input_csv: str | Path, out_csv: str | Path, log=None) -> dict:
        """执行（或跳过）脱靶检测，写出统一记录 CSV，返回状态 dict。"""
        import logging
        log = log or logging.getLogger("sirna.offtarget")
        input_csv, out_csv = Path(input_csv), Path(out_csv)
        if not input_csv.exists():
            raise FileNotFoundError(f"输入记录不存在：{input_csv}")
        status = self.status()

        rows = records.read_records(input_csv)
        # seed6 恒可计算（g2–g7），先统一填充
        for r in rows:
            r["seed6"] = duplex.guide_seed(r["guide_checked"], "seed6")
            r["offtarget_available"] = "0"

        if not status["ok"]:
            # graceful skip：不伪造任何预测列
            log.info("offtarget skip：%s", status["reason"])
            records.write_records(out_csv, rows)
            write_manifest(out_csv.parent, {
                "stage": "05_offtarget", "rows": len(rows),
                "status": status, "mode": "skip",
            })
            return status

        # 外部执行器（就绪）
        os.makedirs(out_csv.parent, exist_ok=True)
        work = out_csv.parent / "_external_result.csv"
        self._run_external(input_csv, work, log)
        ext = records.read_records(work)
        ext_by = {r["variant_id"]: r for r in ext if r.get("variant_id")}
        n_map = 0
        for r in rows:
            e = ext_by.get(r["variant_id"])
            if e is None:
                continue
            n_map += 1
            for col in _OFFTARGET_OUT_COLS:
                if col in e and e[col] not in ("", None):
                    r[col] = e[col]
            r["offtarget_available"] = "1"
        if n_map != len(rows):
            raise ValueError(
                f"offtarget external 结果行数与输入不一致：{n_map}/{len(rows)}")
        records.write_records(out_csv, rows)
        write_manifest(out_csv.parent, {
            "stage": "05_offtarget", "rows": len(rows), "mapped": n_map,
            "status": status, "mode": "external",
        })
        return status


def run_offtarget(
    input_csv: str | Path,
    out_dir: str | Path,
    cfg: dict | None = None,
    log=None,
) -> tuple[Path, dict]:
    """便捷入口：返回 (输出 CSV, 状态)。"""
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    adapter = OffTargetAdapter(cfg)
    out_csv = out_dir / "candidates_offtarget.csv"
    status = adapter.run(input_csv, out_csv, log)
    return out_csv, status


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage5 脱靶检测适配：seed6 + PITA/TargetScan（外部可选）")
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--enabled", action="store_true", help="开启（需 --external-cmd）")
    ap.add_argument("--external-cmd", type=str, default=None,
                    help="外部执行器（见模块 docstring 契约）")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.offtarget")
    cfg = {"enabled": args.enabled, "external_cmd": args.external_cmd}
    try:
        out, status = run_offtarget(args.input, args.out_dir, cfg=cfg, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("offtarget 失败：%s", e)
        return 3
    log.info("offtarget 完成：%s（status=%s）", out, status["reason"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
