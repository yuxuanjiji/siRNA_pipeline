# -*- coding: utf-8 -*-
"""在更新后的错配数据集上重评两个模型（并核查与旧 52 条训练集的窗口级重叠）。

  * OligoFormer（best_model.pth）直接推理：复用 stages/oligoformer/executor.py 的
    FM 嵌入 + 官方 loader 口径（不新增模型、不改仓库）；
  * 位置特征 Ridge：外层 LOSO by dataset / 内层 GroupKFold 选 alpha（组安全）；
  * 泄漏核查：序列级、靶窗级、家族级三档重叠。

输出：outputs/analysis/mismatch/oligoformer_on_updated_mismatch.json
状态：research_only
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
BANDS = ((1, 2), (3, 8), (9, 12), (13, 19))
ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0, 300.0, 1000.0)
REPO = ROOT / "OligoFormer部分"
OLD_FILES = [REPO / "data" / "mismatch.csv",
             REPO / "data" / "mismatch_td.csv",
             REPO / "mismatch_predictions.csv",
             ROOT.parent / "项目搭建" / "热力学参数计算" / "mismatch_validated.csv"]


def norm(s):
    return re.sub(r"\s+", "", str(s).upper()).replace("T", "U")


def spearman(a, b):
    a, b = pd.Series(list(a), dtype=float), pd.Series(list(b), dtype=float)
    if len(a) < 3 or a.nunique() < 2 or b.nunique() < 2:
        return None
    return float(a.corr(b, method="spearman"))


def feats(row):
    pos = sorted(int(x) for x in ast.literal_eval(row["mismatch_positions"]))
    bands = [float(sum(a <= p <= b for p in pos)) for a, b in BANDS]
    count = [float(row["mismatch_count"]), float(len(pos))]
    return count + bands + [float(p in pos) for p in range(1, 20)]


def ridge_loso(data: pd.DataFrame, group_col="group_id") -> dict:
    """外层 LOSO（多数据集）或 LOFO（单数据集，按 family group）；内层 GroupKFold 选 alpha。"""
    data = data.reset_index(drop=True)
    x = np.asarray([feats(r) for _, r in data.iterrows()])
    y = data["label"].to_numpy(dtype=float)
    out, oof = {}, {}
    if data["dataset"].nunique() == 1:                    # 单数据集 → 外层按 group 留一
        gkf = GroupKFold(n_splits=min(5, data[group_col].nunique()))
        outer = [("fold%d" % i, te, np.setdiff1d(np.arange(len(data)), te))
                 for i, (_, te) in enumerate(gkf.split(x, groups=data[group_col]))]
    else:
        outer = [(ds, np.where(data["dataset"].to_numpy() == ds)[0],
                  np.where(data["dataset"].to_numpy() != ds)[0])
                 for ds in sorted(data["dataset"].unique())]
    for name, te_idx, tr_idx in outer:
        if len(tr_idx) < 5 or len(te_idx) < 3:
            continue
        inner = []
        k = min(5, data.loc[tr_idx, group_col].nunique())
        if k >= 2:
            for a, b in GroupKFold(n_splits=k).split(x[tr_idx], groups=data.loc[tr_idx, group_col]):
                inner.append((a, b))
        best, best_s = ALPHAS[0], -np.inf
        for al in ALPHAS:
            sc = []
            for a, b in inner:
                m = make_pipeline(StandardScaler(), Ridge(alpha=al)).fit(x[tr_idx][a], y[tr_idx][a])
                s = spearman(y[tr_idx][b], m.predict(x[tr_idx][b]))
                if s is not None:
                    sc.append(s)
            if sc and float(np.mean(sc)) > best_s:
                best, best_s = al, float(np.mean(sc))
        model = make_pipeline(StandardScaler(), Ridge(alpha=best)).fit(x[tr_idx], y[tr_idx])
        pred = model.predict(x[te_idx])
        for i, p in zip(te_idx, pred):
            oof[int(i)] = float(p)
        out[name] = {"n": int(len(te_idx)), "alpha": best, "spearman": spearman(y[te_idx], pred)}
    idx = sorted(oof)
    pooled = spearman([y[i] for i in idx], [oof[i] for i in idx]) if idx else None
    return {"per_fold_or_dataset": out, "pooled_spearman": pooled, "n_evaluated": len(idx)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", type=Path,
                    default=ROOT / "outputs" / "analysis" / "mismatch" / "mismatch_curated_td.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "outputs" / "analysis" / "mismatch" /
                            "oligoformer_on_updated_mismatch.json")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--batch", type=int, default=32)
    args = ap.parse_args()

    data = pd.read_csv(args.curated, dtype=str)
    data = data[data["model_eligible"] == "True"].reset_index(drop=True)
    for c in ("label", "dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content"):
        if c in data:
            data[c] = pd.to_numeric(data[c], errors="coerce")
    data["uid"] = data["dataset"] + "|" + data["row_id"]

    # ---------- 泄漏核查：与旧 52 条训练集的三档重叠 ----------
    old_seq, old_win, old_fam = set(), set(), set()
    for p in OLD_FILES:
        if not p.exists():
            continue
        rows = pd.read_csv(p, dtype=str, encoding="utf-8-sig").fillna("")
        if "siRNA" in rows:
            old_seq |= {norm(s) for s in rows["siRNA"]}
        if "mRNA" in rows:
            old_win |= {norm(m)[19:38] for m in rows["mRNA"] if len(norm(m)) >= 38}
        if "window" in rows:
            old_win |= {norm(w) for w in rows["window"]}
        if "sirna_id" in rows:
            old_fam |= {re.sub(r"-.*$", "", str(x)) for x in rows["sirna_id"]}
    leak = {
        "old_unique_sequences": len(old_seq),
        "sequence_overlap": int(sum(1 for s in data["siRNA"] if norm(s) in old_seq)),
        "target_window_overlap": int(sum(1 for w in data["target"] if norm(w) in old_win)),
        "family_overlap": int(sum(1 for f in data["family"] if f in old_fam)),
    }
    print("泄漏核查:", leak)

    # ---------- OligoFormer 在更新后数据上的推理 ----------
    from sirna_pipeline.stages.oligoformer import executor as ex
    ckpt = ex._find_ckpt(REPO, None)
    if ckpt is None:
        raise SystemExit("未找到 RNA-FM 权重")
    df = pd.DataFrame({"siRNA": data["siRNA"].map(norm), "mRNA": data["mRNA"].map(norm)})
    df = df.astype(object)
    sys.path.insert(0, str(REPO / "scripts"))
    import infer as infer_mod  # noqa: E402
    df = infer_mod.calculate_td(df).reset_index(drop=True)
    df["siRNA"] = df["siRNA"].astype(str)
    df["mRNA"] = df["mRNA"].astype(str)

    uniq = {"siRNA": {}, "mRNA": {}}
    for s in df["siRNA"]:
        uniq["siRNA"].setdefault(ex._md5seq(s), s)
    for s in df["mRNA"]:
        uniq["mRNA"].setdefault(ex._md5seq(s), s)
    cache = REPO / "data" / "RNAFM_cache"
    ex._embed_repo_fm(uniq, cache, REPO, ckpt, seed=42, batch_size=args.batch,
                      log=lambda *a: print("[fm]", *a, flush=True), device=args.device)
    token = "mm_eval_%d" % os.getpid()
    data_dir = REPO / "data" / "infer" / token
    shutil.rmtree(data_dir, ignore_errors=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    ex._materialize(data_dir, cache, uniq, log=lambda *a: None)
    try:
        scores = ex._infer(df, token, REPO, REPO / "model" / "best_model.pth", 42,
                           device=args.device)
    finally:
        shutil.rmtree(data_dir, ignore_errors=True)
    data["oligo_efficacy"] = [float(v) for v in scores]

    # ---------- 评估 ----------
    report = {"n_rows": len(data), "leakage": leak,
              "oligoformer": {}, "ridge_position": {}, "status": "research_only"}
    for tag, sub in (("all", data),
                     ("ohnishi_only", data[data["dataset"] == "Ohnishi_2008"]),
                     ("exclude_old_window", data[~data["target"].map(norm).isin(old_win)])):
        if len(sub) < 5:
            continue
        report["oligoformer"][tag] = {
            "n": len(sub),
            "spearman": spearman(sub["label"], sub["oligo_efficacy"]),
            "pearson": float(np.corrcoef(sub["label"], sub["oligo_efficacy"])[0, 1])
            if sub["label"].nunique() > 1 else None,
            "rmse": float(np.sqrt(np.mean((sub["label"] - sub["oligo_efficacy"]) ** 2))),
            "per_dataset": {ds: spearman(g["label"], g["oligo_efficacy"])
                            for ds, g in sub.groupby("dataset")},
        }
        report["ridge_position"][tag] = ridge_loso(sub)
        print("\n[%s] n=%d" % (tag, len(sub)))
        print("  OligoFormer: ρ=%.3f | Ridge(位置,LOSO pooled): %s"
              % (report["oligoformer"][tag]["spearman"] or float("nan"),
                 report["ridge_position"][tag]["pooled_spearman"]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nwritten:", args.out)


if __name__ == "__main__":
    main()
