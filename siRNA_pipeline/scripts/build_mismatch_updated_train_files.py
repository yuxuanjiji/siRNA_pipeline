# -*- coding: utf-8 -*-
"""生成更新后的错配训练文件（不覆盖旧文件）。

输出（OligoFormer部分/data/）：
  mismatch_updated.csv      # siRNA,mRNA,label,y,td（与旧 mismatch.csv 同列，供 train_single 使用）
  mismatch_updated_td.csv   # 追加 dG_total/dG_seed/delta_deltaG_ends/MFE_guide/GC_content
y 阈值：沿用旧 mismatch.csv 中 label↔y 的实际分界（自动推断并打印）。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "OligoFormer部分"
CURATED = ROOT / "outputs" / "analysis" / "mismatch" / "mismatch_curated_td.csv"
OLD = REPO / "data" / "mismatch.csv"


def infer_threshold() -> float:
    if not OLD.exists():
        return 0.7
    d = pd.read_csv(OLD, dtype=str)
    lab = pd.to_numeric(d["label"], errors="coerce")
    y = pd.to_numeric(d["y"], errors="coerce")
    ones = lab[y == 1].min()
    zeros = lab[y == 0].max() if (y == 0).any() else None
    if pd.notna(ones) and zeros is not None and pd.notna(zeros):
        print("旧文件 label/y 分界: y=0 最大 %.4f | y=1 最小 %.4f -> 阈值取 %.2f"
              % (zeros, ones, (zeros + ones) / 2))
        return float((zeros + ones) / 2)
    return 0.7


def main() -> None:
    thr = infer_threshold()
    data = pd.read_csv(CURATED, dtype=str)
    data = data[data["model_eligible"] == "True"].reset_index(drop=True)
    label = pd.to_numeric(data["label"], errors="coerce")

    sys.path.insert(0, str(REPO / "scripts"))
    import infer as infer_mod  # noqa: E402

    df = pd.DataFrame({"siRNA": data["siRNA"].str.upper().str.replace("T", "U"),
                       "mRNA": data["mRNA"].str.upper().str.replace("T", "U")}).astype(object)
    df = infer_mod.calculate_td(df).reset_index(drop=True)
    df["siRNA"] = df["siRNA"].astype(str)
    df["mRNA"] = df["mRNA"].astype(str)

    out = pd.DataFrame({
        "siRNA": df["siRNA"],
        "mRNA": df["mRNA"],
        "label": label.round(6),
        "y": (label >= thr).astype(int),
        "td": [",".join("%.6g" % float(v) for v in row) for row in df["td"]],
    })
    out.to_csv(REPO / "data" / "mismatch_updated.csv", index=False, encoding="utf-8-sig")

    td2 = out.copy()
    for col in ("dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content"):
        if col in data:
            td2[col] = pd.to_numeric(data[col], errors="coerce")
    td2.to_csv(REPO / "data" / "mismatch_updated_td.csv", index=False, encoding="utf-8-sig")

    summary = {
        "rows": len(out), "threshold_y": thr,
        "y1": int(out["y"].sum()), "y0": int((1 - out["y"]).sum()),
        "by_dataset": data["dataset"].value_counts().to_dict(),
        "old_file_kept": str(OLD.name),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print("written:", REPO / "data" / "mismatch_updated.csv")
    print("written:", REPO / "data" / "mismatch_updated_td.csv")


if __name__ == "__main__":
    main()
