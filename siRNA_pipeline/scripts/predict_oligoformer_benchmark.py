# -*- coding: utf-8 -*-
"""为任意基准集生成 OligoFormer 预测（补齐 task14 所需的 <name>_predictions.csv）。

口径与 pipeline 的 stages/oligoformer/executor.py 完全一致：
  siRNA(19nt) + mRNA(57nt，允许 X 填充) → 官方 infer.calculate_td → RNA-FM 嵌入
  → scripts/loader.data_process_loader_infer → Oligo(best_model.pth) → pred[:,1]×1.341

输出 CSV 列：siRNA,mRNA,label,pred
  * siRNA/mRNA 原样写出（task14 按 (siRNA,mRNA) 精确键 join，不能改字符串）；
  * 嵌入前内部统一大写 / T→U。

用法（必须用 OligoFormer 的 venv python：fm 依赖 ptflops）：
  OligoFormer部分/.venv/Scripts/python.exe scripts/predict_oligoformer_benchmark.py \
      --input ../数据集/最终数据/Simone.csv --out OligoFormer部分/Simone_predictions.csv
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True, help="含 siRNA + mRNA(57nt)[ + label] 的 CSV")
    ap.add_argument("--out", type=Path, required=True, help="输出 <name>_predictions.csv")
    ap.add_argument("--repo", type=Path, default=ROOT / "OligoFormer部分")
    ap.add_argument("--model", type=Path, default=None, help="默认 <repo>/model/best_model.pth")
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--device", default="cpu", help="cpu 或 cuda")
    args = ap.parse_args()

    from sirna_pipeline.stages.oligoformer import executor as ex

    repo = args.repo.resolve()
    model = (args.model or (repo / "model" / "best_model.pth")).resolve()
    rows = list(csv.DictReader(open(args.input, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit("输入为空")
    raw_si = [str(r.get("siRNA", "")) for r in rows]
    raw_mr = [str(r.get("mRNA", "")) for r in rows]
    labels = [str(r.get("label", r.get("silencing_efficiency_norm", ""))) for r in rows]

    def norm(s: str) -> str:
        return "".join(s.split()).upper().replace("T", "U")

    df = pd.DataFrame({"siRNA": [norm(s) for s in raw_si],
                       "mRNA": [norm(m) for m in raw_mr]}).astype(object)
    sys.path.insert(0, str(repo / "scripts"))
    import infer as infer_mod  # type: ignore
    df = infer_mod.calculate_td(df).reset_index(drop=True)
    for c in ("siRNA", "mRNA"):
        df[c] = df[c].astype(str)

    uniq = {"siRNA": {}, "mRNA": {}}
    for s in df["siRNA"]:
        uniq["siRNA"].setdefault(ex._md5seq(s), s)
    for s in df["mRNA"]:
        uniq["mRNA"].setdefault(ex._md5seq(s), s)
    print("rows=%d | 唯一 siRNA=%d | 唯一 mRNA=%d"
          % (len(df), len(uniq["siRNA"]), len(uniq["mRNA"])), flush=True)

    ckpt = ex._find_ckpt(repo, None)
    if ckpt is None:
        raise SystemExit("未找到 RNA-FM 权重")
    cache = repo / "data" / "RNAFM_cache"
    ex._embed_repo_fm(uniq, cache, repo, ckpt, seed=42, batch_size=args.batch,
                      log=lambda *a: print("[fm]", *a, flush=True), device=args.device)

    token = "bench_%d" % os.getpid()
    data_dir = repo / "data" / "infer" / token
    shutil.rmtree(data_dir, ignore_errors=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    ex._materialize(data_dir, cache, uniq, log=lambda *a: None)
    try:
        scores = ex._infer(df, token, repo, model, 42, device=args.device)
    finally:
        shutil.rmtree(data_dir, ignore_errors=True)
    if len(scores) != len(rows):
        raise SystemExit("预测行数不一致：%d != %d" % (len(scores), len(rows)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["siRNA", "mRNA", "label", "pred"])
        for si, mr, lab, p in zip(raw_si, raw_mr, labels, scores):
            w.writerow([si, mr, lab, "%.8f" % p])
    print("written:", args.out, "| pred range %.3f–%.3f" % (min(scores), max(scores)))


if __name__ == "__main__":
    main()
