# -*- coding: utf-8 -*-
"""候选生成阶段运行器（Stage 1 · generation）。

职责：调用 legacy 核心（task5/task6 原脚本，通过 subprocess 黑盒执行、算法零改动），
再把两张旧格式 CSV 转换为统一 CandidateRecord 记录，写入 runs/01_generation/。

为什么不直接 import legacy 脚本？
  - task6 依赖同目录 task5/seq_utils 的顶层 import；以 subprocess 运行能 100% 保留
    原行为与内置校验（“不改原算法、不复制拼接代码”原则）。
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from ...common import records
from ...common.stage_io import file_sha256, write_manifest

LEGACY_DIR = Path(__file__).resolve().parent / "legacy"
TASK5_SCRIPT = LEGACY_DIR / "task5_sliding_window_siRNA.py"
TASK6_SCRIPT = LEGACY_DIR / "task6_site_mutagenesis.py"
TASK5_DEFAULT_NAME = "SFRP1_task5_sliding_window_19nt_fullcomp.csv"
TASK6_DEFAULT_NAME = "SFRP1_task6_site_mutagenesis_library.csv"

# legacy 旧列 -> 统一列（write 前先由 dict 名映射）
_OUT_FIELDS = records.BASE_COLS + records.MUT_COLS


def _run_script(script: Path, args: list[str], python: str, log) -> None:
    """运行 legacy 脚本并透传退出/输出；失败抛 RuntimeError（中文信息）。"""
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    cmd = [python, str(script)] + [str(a) for a in args]
    log.info("run: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=env, timeout=600)
    tail = (proc.stdout + proc.stderr)[-2000:]
    if proc.returncode != 0:
        raise RuntimeError(
            f"legacy 脚本运行失败（exit={proc.returncode}）：{script.name}\n{tail}")
    log.info("legacy %s 完成，stdout 尾部:\n%s", script.name, proc.stdout[-800:])


# ---------------------------------------------------------------------------
# 转换器：旧 CSV 行 -> 统一记录行
# ---------------------------------------------------------------------------
def wt_rows_from_task5(t5_rows: list[dict]) -> list[dict]:
    out = []
    for r in t5_rows:
        wid = r["window_id"]
        guide = r["guide_antisense_19"]
        out.append({
            "window_id": wid, "variant_id": f"{wid}_wt", "kind": "wt",
            "cds_start": r.get("cds_start", ""), "cds_end": r.get("cds_end", ""),
            "nm_003012_start": r.get("nm_003012_start", ""),
            "nm_003012_end": r.get("nm_003012_end", ""),
            "target_mRNA_19": r.get("target_mRNA_19", ""),
            "guide_wt_19": guide, "guide_checked": guide,
            "sense_strand_19": r.get("sense_strand_19", ""),
        })
    return out


def mut_rows_from_task6(t6_rows: list[dict], nm_by_wid: dict) -> list[dict]:
    out = []
    for r in t6_rows:
        wid = r["window_id"]
        nm = nm_by_wid.get(wid, {})
        out.append({
            "window_id": wid, "variant_id": r["variant_id"], "kind": "mut",
            "cds_start": r.get("cds_start", ""), "cds_end": r.get("cds_end", ""),
            "nm_003012_start": nm.get("nm_003012_start", ""),
            "nm_003012_end": nm.get("nm_003012_end", ""),
            "target_mRNA_19": r.get("target_mRNA_19", ""),
            "guide_wt_19": r.get("guide_wt_19", ""),
            "guide_checked": r.get("guide_mut_19", ""),
            "sense_strand_19": "",          # t6 不输出 sense；wt 行提供
            "site": r.get("site", ""), "position": r.get("position", ""),
            "wt_nt": r.get("wt_nt", ""), "mut_nt": r.get("mut_nt", ""),
            "paired_mRNA_nt": r.get("paired_mRNA_nt", ""),
            "pair_type": r.get("pair_type", ""),
        })
    return out


def merge_and_sort(wt_rows: list[dict], mut_rows: list[dict]) -> list[dict]:
    """合并排序：按任务5窗口顺序逐窗口输出；窗口内 wt 在前、15 条变体保持任务6文件内 v01..v15 顺序。"""
    wt_by_wid = {r["window_id"]: r for r in wt_rows}
    mut_by_wid: dict[str, list[dict]] = {}
    for r in mut_rows:
        mut_by_wid.setdefault(r["window_id"], []).append(r)
    flat: list[dict] = []
    for wid, wt_r in wt_by_wid.items():
        flat.append(wt_r)
        flat.extend(mut_by_wid.get(wid, []))
    if len(flat) != len(wt_rows) + len(mut_rows):
        raise AssertionError("窗口合并后行数不一致（存在未知窗口的突变体？）")
    return flat


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def run_generation(
    fasta: str | Path,
    out_dir: str | Path,
    work_dir: str | Path | None = None,
    python: str | None = None,
    g1_wobble_only: bool = False,
    log=None,
) -> Path:
    """执行候选生成并写出统一记录。

    返回统一记录 CSV 路径（runs/<out_dir>/candidates.csv）。
    """
    import logging
    log = log or logging.getLogger("sirna.generation")
    fasta = Path(fasta)
    if not fasta.exists():
        raise FileNotFoundError(f"SFRP1 FASTA 不存在：{fasta}")
    python = python or sys.executable
    out_dir = Path(out_dir)
    work_dir = Path(work_dir) if work_dir else out_dir / "_legacy_work"
    for d in (out_dir, work_dir):
        os.makedirs(d, exist_ok=True)

    # 1) 调用 legacy task5 / task6（黑盒，原校验保留）
    t5_out = work_dir / TASK5_DEFAULT_NAME
    t6_out = work_dir / TASK6_DEFAULT_NAME
    _run_script(TASK5_SCRIPT, ["--fasta", fasta, "--out", t5_out],
                python, log)
    task6_args = ["--input", t5_out, "--out", t6_out]
    if g1_wobble_only:
        task6_args.append("--g1-wobble-only")
    _run_script(TASK6_SCRIPT, task6_args, python, log)

    # 2) 转换旧格式 -> 统一记录
    from ...common.seqio import read_csv_rows
    t5 = read_csv_rows(t5_out)
    t6 = read_csv_rows(t6_out)
    nm_by_wid = {r["window_id"]: r for r in t5}
    wt = wt_rows_from_task5(t5)
    mut = mut_rows_from_task6(t6, nm_by_wid)
    rows = merge_and_sort(wt, mut)

    # 3) 结构校验：主键/kind/序列 19nt + 每窗口 1wt+15mut
    records.validate_base(rows, "generation")
    from collections import Counter
    per_win = Counter(r["window_id"] for r in rows)
    n_win = len(per_win)
    bad_wins = {k: v for k, v in per_win.items() if v != (16 if not g1_wobble_only else -1)}
    if not g1_wobble_only and bad_wins:
        raise AssertionError(f"窗口内记录数 != 16：{list(bad_wins.items())[:5]}")
    n_wt = sum(1 for r in rows if r["kind"] == "wt")
    n_mut = len(rows) - n_wt
    log.info("generation: windows=%d wt=%d mut=%d total=%d",
             n_win, n_wt, n_mut, len(rows))

    # 4) 写统一记录 + manifest
    out_csv = out_dir / "candidates.csv"
    records.write_records(out_csv, rows)
    write_manifest(out_dir, {
        "stage": "01_generation", "windows": n_win,
        "wt": n_wt, "mut": n_mut, "total": len(rows),
        "g1_wobble_only": g1_wobble_only,
        "fasta_sha256": file_sha256(fasta),
        "input_fasta": str(fasta),
    })
    return out_csv


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Stage1 候选生成：task5 滑窗 + task6 定点突变 → 统一记录 CSV")
    ap.add_argument("--fasta", type=Path, required=True,
                    help="SFRP1 mRNA FASTA（通常 ../数据集/SFRP1-mRNA.txt）")
    ap.add_argument("--out-dir", type=Path, required=True,
                    help="统一记录输出目录（将写入 candidates.csv + manifest.json）")
    ap.add_argument("--work-dir", type=Path, default=None,
                    help="legacy 临时 CSV 目录（默认 out-dir/_legacy_work）")
    ap.add_argument("--python", type=str, default=None,
                    help="运行 legacy 脚本的解释器（默认当前解释器）")
    ap.add_argument("--g1-wobble-only", action="store_true",
                    help="透传给 task6：g1 只保留 G:U 摆动（默认与任务6一致不过滤）")
    args = ap.parse_args(argv)
    from ...common.stage_io import setup_logging
    log = setup_logging("sirna.generation")
    try:
        out = run_generation(fasta=args.fasta, out_dir=args.out_dir,
                             work_dir=args.work_dir, python=args.python,
                             g1_wobble_only=args.g1_wobble_only, log=log)
    except Exception as e:  # noqa: BLE001
        log.error("generation 失败：%s", e)
        return 3
    log.info("generation 完成：%s", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
