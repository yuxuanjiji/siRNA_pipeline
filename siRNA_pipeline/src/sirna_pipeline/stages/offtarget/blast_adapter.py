# -*- coding: utf-8 -*-
"""BLAST 近全长同源脱靶适配层（脱靶检测第三层 · BLAST）。

背景
----
`BLAST/blast_offtarget.py` 是本项目脱靶验证的第三层（与 offtarget/adapter.py 的
PITA/TargetScan 种子区层平行）：用 blastn-short（word_size 7）对人类转录组做
"近全长反向互补"同源搜索，再按 tolerant 判据（非靶基因 / 反义方向 /
配对数>=17 / 错配<=1）标记候选（口径详见 BLAST/README.txt，原文件零修改）。

本适配层沿用例规约：
  1) `BLAST/` 原脚本只读引用（不改一行），以子进程调用（复用其 CLI 契约）；
  2) 默认 disabled → graceful skip：写出 offtarget_blast_available=0、不伪造判据；
  3) enabled 但 BLAST+ / 库未就绪 → status 明确报原因，仍 graceful skip 不阻断管道；
  4) 结果回填统一记录列：offtarget_blast_flag / offtarget_blast_n_hits_ge16 /
     offtarget_blast_available，并（flag=1 时）置共享列 offtarget_risk=1 供排序软惩罚。

外部依赖（enabled 时才需要）：python 子进程环境装有 pandas；NCBI BLAST+（blastn）可执行；
BLAST 库（makeblastdb 产物，如 gencode_v46_pc_simple）。库与 bin 均可通过 config 传入。
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from ...common import records
from ...common.stage_io import write_manifest

# 默认定位工程根 BLAST/blast_offtarget.py（siRNA_pipeline/BLAST）
_PROJECT_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_SCRIPT = _PROJECT_ROOT / "BLAST" / "blast_offtarget.py"

DEFAULT_CFG = {
    "enabled": False,             # 未配置 BLAST+/库前默认关闭
    "script": str(DEFAULT_SCRIPT),
    "python": None,               # None -> 当前解释器 sys.executable
    "blast_bin": None,            # BLAST+ bin 目录（None -> PATH 探测；可用 script 内默认）
    "db": None,                   # BLAST 库前缀（makeblastdb 产物；ready 判据之一）
    "target_gene": "SFRP1",
    "min_paired": 17,
    "max_mismatch": 1,
    "require_no_gaps": False,
    "threads": 8,
}

_DB_EXTS = (".nsq", ".nhr", ".nin")     # makeblastdb（nucl）v4 产物；v5(.ndb) 亦接受


def _blastn_path(cfg: dict) -> str | None:
    """返回可用的 blastn 可执行（bin 目录优先，其次 PATH）。"""
    bin_dir = cfg.get("blast_bin")
    if bin_dir:
        for name in ("blastn.exe", "blastn"):
            p = Path(str(bin_dir)) / name
            if p.exists():
                return str(p)
        return None
    return shutil.which("blastn")


def _db_ready(db: str | None) -> bool:
    if not db:
        return False
    pre = str(db)
    for ext in _DB_EXTS + (".ndb",):
        if Path(pre + ext).exists():
            return True
    return False


def _detection_status(cfg: dict) -> dict:
    """探测外部依赖并给出状态/原因（不抛异常）。"""
    st = {
        "enabled": bool(cfg.get("enabled")),
        "script": cfg.get("script"),
        "script_exists": False,
        "python": cfg.get("python") or sys.executable,
        "blastn": None,
        "db": cfg.get("db"),
        "db_ready": False,
        "ok": False, "reason": "",
    }
    sp = Path(str(cfg.get("script", "")))
    st["script_exists"] = sp.exists()
    if not st["enabled"]:
        st["reason"] = "disabled（默认：未配置 BLAST+/库，近全长脱靶层跳过）"
        return st
    if not st["script_exists"]:
        st["reason"] = f"BLAST 脚本不存在：{sp}"
        return st
    st["blastn"] = _blastn_path(cfg)
    st["db_ready"] = _db_ready(cfg.get("db"))
    if not st["blastn"]:
        st["reason"] = ("blastn 不可执行：blast_bin 目录无 blastn(.exe) 且 PATH 中未找到；"
                        "enabled 需在 config 填 blast_bin/db")
        return st
    if not st["db_ready"]:
        st["reason"] = f"BLAST 库未就绪：db={cfg.get('db')}（需 makeblastdb 产物）"
        return st
    st["reason"] = "ready"
    st["ok"] = True
    return st


class BlastOffTargetAdapter:
    """BLAST 近全长脱靶适配器：探测 + 批量执行（子进程调用 BLAST/blast_offtarget.py）。"""

    def __init__(self, cfg: dict | None = None):
        self.cfg = {**DEFAULT_CFG, **(cfg or {})}

    def status(self) -> dict:
        return _detection_status(self.cfg)

    # ------------------------------------------------------------------
    def _run_script(self, work_csv: Path, out_dir: Path, prefix: str, log) -> None:
        """按 BLAST/blast_offtarget.py 的 CLI 契约调用（原脚本不改）。"""
        cfg = self.cfg
        cmd = [cfg.get("python") or sys.executable,
               cfg.get("script") or str(DEFAULT_SCRIPT),
               "--input", str(work_csv), "--column", "siRNA",
               "--out-dir", str(out_dir), "--prefix", prefix,
               "--target-gene", str(cfg.get("target_gene", "SFRP1")),
               "--min-paired", str(int(cfg.get("min_paired", 17))),
               "--max-mismatch", str(int(cfg.get("max_mismatch", 1))),
               "--threads", str(int(cfg.get("threads", 8)))]
        if cfg.get("blast_bin"):
            cmd += ["--blast-bin", str(cfg["blast_bin"])]
        if cfg.get("db"):
            cmd += ["--db", str(cfg["db"])]
        if cfg.get("require_no_gaps"):
            cmd += ["--require-no-gaps"]
        log.info("blast_offtarget external: %s", " ".join(str(c) for c in cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace",
                              timeout=7200)
        if proc.returncode != 0:
            raise RuntimeError(
                f"blast_offtarget 失败（exit={proc.returncode}）："
                f"{(proc.stdout + proc.stderr)[-2000:]}")

    # ------------------------------------------------------------------
    @staticmethod
    def _ids_from(path: Path, key: str = "variant_id") -> set[str]:
        if not path.exists():
            return set()
        with open(path, encoding="utf-8-sig", newline="") as fh:
            rd = csv.DictReader(fh)
            if rd.fieldnames is None or key not in rd.fieldnames:
                return set()
            return {r[key] for r in rd if r.get(key, "") not in ("", "nan")}

    @staticmethod
    def _hit_counts(hits_path: Path) -> dict[str, int]:
        """hits_all.tsv：按 qseqid 统计全部 >=16nt 命中数（审计留痕）。"""
        counts: dict[str, int] = {}
        if not hits_path.exists():
            return counts
        with open(hits_path, encoding="utf-8-sig", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                q = r.get("qseqid", "")
                counts[q] = counts.get(q, 0) + 1
        return counts

    # ------------------------------------------------------------------
    def run(self, input_csv: str | Path, out_csv: str | Path, log=None) -> dict:
        """执行（或跳过）BLAST 近全长脱靶层，写出统一记录 CSV，返回状态 dict。

        输入记录表行经适配后回填 BLAST 层列；flag=1 的候选同时置共享列
        offtarget_risk=1（供 rank 软惩罚；未命中保持原值/置 0）。
        """
        import logging
        log = log or logging.getLogger("sirna.blast_offtarget")
        input_csv, out_csv = Path(input_csv), Path(out_csv)
        if not input_csv.exists():
            raise FileNotFoundError(f"输入记录不存在：{input_csv}")
        status = self.status()

        rows = records.read_records(input_csv)
        for r in rows:
            r["offtarget_blast_available"] = "0"

        if not status["ok"]:
            log.info("blast_offtarget skip：%s", status["reason"])
            records.write_records(out_csv, rows)
            write_manifest(out_csv.parent, {
                "stage": "05_offtarget_blast", "rows": len(rows),
                "status": status, "mode": "skip",
            })
            return status

        # ---- 就绪：构建送检输入并调用原脚本 ----
        prefix = "offtarget_blast"
        work = out_csv.parent / ("_blast_" + prefix)
        work.mkdir(parents=True, exist_ok=True)
        work_csv = work / "input_candidates.csv"
        with open(work_csv, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["variant_id", "siRNA"])
            for r in rows:
                w.writerow([r.get("variant_id", ""), r.get("guide_checked", "")])

        self._run_script(work_csv, work, prefix, log)

        # ---- 解析：flagged/clean + hits + stats ----
        flagged = self._ids_from(work / f"{prefix}_flagged.csv")
        clean = self._ids_from(work / f"{prefix}_clean.csv")
        hits = self._hit_counts(work / f"{prefix}_hits_all.tsv")
        stats = {}
        sp = work / f"{prefix}_stats.json"
        if sp.exists():
            try:
                stats = json.loads(sp.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                stats = {}

        n_flagged = int(stats.get("n_flagged_offtarget", len(flagged)))
        for i, r in enumerate(rows):
            vid = r.get("variant_id", "")
            is_flag = vid in flagged
            r["offtarget_blast_flag"] = "1" if is_flag else "0"
            r["offtarget_blast_n_hits_ge16"] = str(hits.get(f"q{i}", 0))
            r["offtarget_blast_available"] = "1"
            # 共享风险列：脱靶命中 → 1（供 rank offtarget_risk 软惩罚）
            r["offtarget_risk"] = "1" if is_flag else (
                "1" if str(r.get("offtarget_risk", "")) == "1" else "0")
        records.write_records(out_csv, rows)
        write_manifest(out_csv.parent, {
            "stage": "05_offtarget_blast",
            "input_csv": str(input_csv),
            "rows": len(rows),
            "flagged": n_flagged,
            "hits_all_ge16": int(stats.get("n_hits_all_ge16", len(hits))),
            "db": status.get("db"),
            "blastn": status.get("blastn"),
            "target_gene": self.cfg.get("target_gene"),
            "criterion": stats.get("criterion"),
            "top_offtarget_genes": stats.get("top_offtarget_genes"),
            "mode": "external",
        })
        log.info("blast_offtarget 完成：rows=%d flagged=%d", len(rows), n_flagged)
        return status


def run_blast_offtarget(
    input_csv: str | Path,
    out_dir: str | Path,
    cfg: dict | None = None,
    log=None,
) -> tuple[Path, dict]:
    """便捷入口：返回 (输出 CSV, 状态)。"""
    out_dir = Path(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    adapter = BlastOffTargetAdapter(cfg)
    out_csv = out_dir / "candidates_blast_offtarget.csv"
    status = adapter.run(input_csv, out_csv, log)
    return out_csv, status


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage5 BLAST 近全长同源脱靶适配：子进程调用 BLAST/blast_offtarget.py")
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--enabled", action="store_true")
    ap.add_argument("--script", type=str, default=str(DEFAULT_SCRIPT))
    ap.add_argument("--python", type=str, default=None)
    ap.add_argument("--blast-bin", type=str, default=None)
    ap.add_argument("--db", type=str, default=None)
    ap.add_argument("--target-gene", type=str, default="SFRP1")
    ap.add_argument("--min-paired", type=int, default=17)
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.blast_offtarget")
    cfg = {"enabled": args.enabled, "script": args.script, "python": args.python,
           "blast_bin": args.blast_bin, "db": args.db,
           "target_gene": args.target_gene, "min_paired": args.min_paired,
           "max_mismatch": args.max_mismatch, "threads": args.threads}
    try:
        out, status = run_blast_offtarget(args.input, args.out_dir, cfg=cfg, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("blast_offtarget 失败：%s", e)
        return 3
    log.info("blast_offtarget 完成：%s（status=%s）", out, status["reason"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
