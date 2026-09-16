# -*- coding: utf-8 -*-
"""Split-safe optimization of the mismatch position model (research-only).

在 prepare_mismatch_dataset.py 的 curated 数据与划分方案上做严格化：
  * uid = dataset|row_id（修跨数据集 row_id 冲突）；
  * 外层 Leave-One-Dataset-Out，内层 GroupKFold(by group_id) 选 alpha（复刻其 splits 逻辑）；
  * alpha 网格扩到 1e4（原实现选到边界 100）；
  * pooled out-of-fold 预测 + 逐数据集指标 + 按 group 的 bootstrap CI + 组内置换检验；
  * 可选 --drop-kini（Kini 注释未解析、位置特征不可信时做敏感性）。

用法：
  python scripts/optimize_mismatch_model.py                 # 全部 102 行
  python scripts/optimize_mismatch_model.py --drop-kini     # 敏感性
状态：research_only，不改动生产排序结果。
"""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, HuberRegressor, Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

BANDS = ((1, 2), (3, 8), (9, 12), (13, 19))
ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0)
COMP = {"A": "U", "U": "A", "C": "G", "G": "C"}


def mtype(gb: str, tb: str) -> str:
    if gb == COMP.get(tb, ""):
        return "WC"
    if {gb, tb} == {"G", "U"}:
        return "GU"
    pg, pm = gb in "AG", tb in "AG"
    if pg and pm:
        return "PP"
    if not pg and not pm:
        return "YY"
    return "PY"


def feature_sets(row: pd.Series, ds_codes: dict) -> dict[str, list[float]]:
    """位置/类型/上下文/热力学四组特征（位置取自论文注释，与位置模型口径一致）。"""
    ann_pos = sorted(int(x) for x in ast.literal_eval(row["mismatch_positions"]))
    types = [mtype(s.split(":")[0], s.split(":")[1]) for s in ast.literal_eval(row["mismatch_types"])
             ] if row.get("mismatch_types") not in ("", None) else []
    bands = [float(sum(a <= p <= b for p in ann_pos)) for a, b in BANDS]
    ones = [float(p in ann_pos) for p in range(1, 20)]
    tcount = {k: float(types.count(k)) for k in ("GU", "PP", "YY", "PY")}

    guide, target = str(row["siRNA"]).upper(), str(row["target"]).upper()
    def weak(p: int) -> float:
        if not 1 <= p <= 19:
            return 0.0
        gb, tb = guide[p - 1], target[19 - p]
        return 1.0 if (gb == COMP.get(tb, "") and gb in "AU") else 0.0
    n_weak = float(sum(weak(p - 1) + weak(p + 1) for p in ann_pos))
    run_ge2 = float(any({p, p + 1} <= set(ann_pos) for p in ann_pos))
    thermo = ([float(row["dG_total"]), float(row["dG_seed"]), float(row["delta_deltaG_ends"]),
               float(row["MFE_guide"]), float(row["GC_content"])]
              if "dG_total" in row and pd.notna(row.get("dG_total")) else [])
    return {
        "bands_pos": bands + ones,
        "bands_only": bands,
        "bands_types_ctx": bands + [tcount["GU"], tcount["PP"], tcount["YY"], tcount["PY"], n_weak, run_ge2],
        "bands_types_ctx_thermo": bands + [tcount["GU"], tcount["PP"], tcount["YY"], tcount["PY"],
                                           n_weak, run_ge2] + thermo,
        "bands_types_ctx_dataset": bands + [tcount["GU"], tcount["PP"], tcount["YY"], tcount["PY"],
                                            n_weak, run_ge2] + ds_codes[row["dataset"]],
    }


def spearman(a, b) -> float | None:
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if len(a) < 3 or len(set(a)) < 2 or len(set(b)) < 2:
        return None
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def fit_predict(x_tr, y_tr, x_te, model, alpha):
    est = (Ridge(alpha=alpha) if model == "ridge"
           else HuberRegressor(alpha=alpha, epsilon=1.35, max_iter=1000) if model == "huber"
           else ElasticNet(alpha=alpha, l1_ratio=0.2, max_iter=20000))
    return make_pipeline(StandardScaler(), est).fit(x_tr, y_tr).predict(x_te)


def select_alpha(x_tr, y_tr, model, inner):
    best, best_s = ALPHAS[0], -np.inf
    for a in ALPHAS:
        sc = [spearman(y_tr[v], fit_predict(x_tr[t], y_tr[t], x_tr[v], model, a))
              for t, v in inner]
        sc = [x for x in sc if x is not None]
        m = float(np.mean(sc)) if sc else -np.inf
        if m > best_s:
            best, best_s = a, m
    return best


def bootstrap_ci(y, p, groups, n=1000, seed=0):
    rng = np.random.default_rng(seed)
    uniq = sorted(set(groups))
    vals = []
    for _ in range(n):
        pick = rng.choice(len(uniq), size=len(uniq), replace=True)
        idx = [i for g in pick for i, gg in enumerate(groups) if gg == uniq[g]]
        s = spearman([y[i] for i in idx], [p[i] for i in idx])
        if s is not None:
            vals.append(s)
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))) if vals else (None, None)


def permutation_p(y, p, groups, obs, n=500, seed=0):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    null = []
    for _ in range(n):
        perm = y.copy()
        for g in sorted(set(groups)):
            idx = [i for i, gg in enumerate(groups) if gg == g]
            perm[idx] = rng.permutation(perm[idx])
        s = spearman(perm, p)
        if s is not None:
            null.append(s)
    return float(np.mean([abs(x) >= abs(obs) for x in null])) if null else None


def main() -> None:
    root = Path(__file__).parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", type=Path,
                    default=root / "outputs" / "analysis" / "mismatch" / "mismatch_curated_td.csv")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--drop-kini", action="store_true", help="敏感性：剔除 Kini（注释未解析）")
    ap.add_argument("--datasets", default=None, help="只保留这些数据集（逗号分隔），如 Ohnishi_2008")
    ap.add_argument("--outer", default="loso", choices=("loso", "lofo"),
                    help="外层划分：loso=按数据集留一；lofo=按 family group 留一")
    ap.add_argument("--legacy25", action="store_true", help="附加复刻其 25 维特征集")
    args = ap.parse_args()

    data = pd.read_csv(args.curated, dtype=str)
    data = data[data["model_eligible"] == "True"].reset_index(drop=True)
    data["uid"] = data["dataset"] + "|" + data["row_id"]
    if args.drop_kini:
        data = data[data["dataset"] != "Kini_2009"].reset_index(drop=True)
    if args.datasets:
        keep = {x.strip() for x in args.datasets.split(",")}
        data = data[data["dataset"].isin(keep)].reset_index(drop=True)
    for col in ("dG_total", "dG_seed", "delta_deltaG_ends", "MFE_guide", "GC_content", "label"):
        if col in data:
            data[col] = pd.to_numeric(data[col], errors="coerce")
    ds_codes = {d: [float(d == x) for x in sorted(data["dataset"].unique())]
                for d in sorted(data["dataset"].unique())}

    feats = {name: {uid: feature_sets(row, ds_codes)[name] for uid, row in
                    data.set_index("uid").iterrows()}
             for name in ("bands_pos", "bands_only", "bands_types_ctx",
                          "bands_types_ctx_thermo", "bands_types_ctx_dataset")}
    if args.legacy25:
        legacy = {}
        for uid, row in data.set_index("uid").iterrows():
            ann = sorted(int(x) for x in ast.literal_eval(row["mismatch_positions"]))
            legacy[uid] = ([float(row["mismatch_count"]), float(len(ann))]
                           + [float(sum(a <= p <= b for p in ann)) for a, b in BANDS]
                           + [float(p in ann) for p in range(1, 20)])
        feats["legacy25"] = legacy

    # ---- 外层划分 + 内层 GroupKFold（复刻其 splits 逻辑，但以 uid 为键）----
    splits = {}
    if args.outer == "loso":
        outer_keys = [("dataset", d) for d in sorted(data["dataset"].unique())]
    else:
        outer_keys = [("group", g) for g in sorted(data["group_id"].unique())]
    for kind, val in outer_keys:
        col = "dataset" if kind == "dataset" else "group_id"
        te = data[data[col] == val]
        tr = data[data[col] != val]
        inner = []
        k = min(5, tr["group_id"].nunique())
        if k >= 2 and len(tr):
            gkf = GroupKFold(n_splits=k)
            for a, b in gkf.split(tr, groups=tr["group_id"]):
                inner.append((tr.iloc[a].index.to_numpy(), tr.iloc[b].index.to_numpy()))
        splits[str(val)] = {"train_uids": tr["uid"].tolist(), "test_uids": te["uid"].tolist(),
                            "inner": inner}
    print("rows=%d groups=%d datasets=%s | outer=%s folds=%d"
          % (len(data), data["group_id"].nunique(), data["dataset"].value_counts().to_dict(),
             args.outer, len(splits)))

    rec = data.set_index("uid")
    label_by = rec["label"].to_dict()
    group_by = rec["group_id"].to_dict()

    report = {"n_rows": len(data), "n_groups": int(data["group_id"].nunique()),
              "drop_kini": bool(args.drop_kini), "variants": {}}
    for fname, fmap in feats.items():
        for model in ("ridge", "huber", "elasticnet"):
            oof, per_ds = {}, {}
            for ds, sp in splits.items():
                tr_ids, te_ids = sp["train_uids"], sp["test_uids"]
                if not tr_ids or not te_ids:
                    continue
                idx_tr = data.index[data["uid"].isin(tr_ids)].to_numpy()
                idx_te = data.index[data["uid"].isin(te_ids)].to_numpy()
                pos_in_tr = {data.loc[i, "uid"]: pos for pos, i in enumerate(idx_tr)}
                x_tr = np.asarray([fmap[data.loc[i, "uid"]] for i in idx_tr])
                y_tr = data.loc[idx_tr, "label"].to_numpy(dtype=float)
                x_te = np.asarray([fmap[data.loc[i, "uid"]] for i in idx_te])
                inner = []
                for t_arr, v_arr in sp["inner"]:
                    t_pos = [pos_in_tr[data.loc[i, "uid"]] for i in t_arr]
                    v_pos = [pos_in_tr[data.loc[i, "uid"]] for i in v_arr]
                    if t_pos and v_pos:
                        inner.append((np.asarray(t_pos), np.asarray(v_pos)))
                alpha = select_alpha(x_tr, y_tr, model, inner)
                pred = fit_predict(x_tr, y_tr, x_te, model, alpha)
                for i, p in zip(idx_te, pred):
                    oof[data.loc[i, "uid"]] = float(p)
                y_te = data.loc[idx_te, "label"].to_numpy(dtype=float)
                per_ds[ds] = {"n": len(idx_te), "alpha": alpha, "spearman": spearman(y_te, pred)}
            uids = [u for u in data["uid"] if u in oof]
            y = [label_by[u] for u in uids]
            p = [oof[u] for u in uids]
            g = [group_by[u] for u in uids]
            pooled = spearman(y, p)
            ci = bootstrap_ci(y, p, g)
            pv = permutation_p(y, p, g, pooled) if pooled is not None else None
            key = "%s|%s" % (fname, model)
            report["variants"][key] = {"feature_set": fname, "model": model,
                                       "pooled_spearman": pooled, "pooled_ci95": ci,
                                       "perm_p": pv, "per_dataset": per_ds}
            print("%-34s pooled=%6.3f CI=[%6.3f,%6.3f] p=%.3f | %s"
                  % (key, pooled if pooled is not None else float("nan"),
                     ci[0] if ci[0] is not None else float("nan"),
                     ci[1] if ci[1] is not None else float("nan"),
                     pv if pv is not None else float("nan"),
                     " ".join("%s=%.2f(n%d)" % (d, v["spearman"] if v["spearman"] is not None else float("nan"), v["n"])
                              for d, v in per_ds.items())))

    report["status"] = "research_only_not_connected_to_production"
    report["note"] = ("外层 LOSO by dataset + 内层 GroupKFold 选 alpha；pooled 指标为两层拆分下 "
                      "out-of-fold 预测汇总；CI 按 group 重采样；置换检验按组内打乱标签。"
                      "Kini n=5 且注释未解析，--drop-kini 为敏感性对照。")
    out = args.out or (root / "outputs" / "analysis" / "mismatch" /
                       ("mismatch_optimization_dropkini.json" if args.drop_kini
                        else "mismatch_optimization.json"))
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("written:", out)


if __name__ == "__main__":
    main()
