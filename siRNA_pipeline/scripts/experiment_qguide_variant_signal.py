# -*- coding: utf-8 -*-
"""前置实验：Q_guide 的变体信号是否可信（"双信号乘性分解"的可行性检验）。

拟采用的模型（第四方案）：

    final_score = Q_guide(突变后序列) × C_match(位置, 错配类型)
    Q_guide     = α·score_thermo + β·score_oligo        (生产口径 α/β = 0.4/0.6)
    C_match     = Π_错配 (1 − λ_position · λ_type)      (WT 无错配 → C_match = 1)

该设计有一个隐含前提：**Q_guide 对单点突变的相对变化方向要可信**。
因为提拔一个 g1 GU 变体只需 Q_guide 比 WT 高 2.6%（C=0.975），
若这点差异被 OligoFormer（错配场景 ρ≈0.21）的跳变淹没，则整个设计不成立。

本实验**不拟合任何参数**，只测量三件事：

  1. ΔQ = Q(变体) − Q(同家族 WT) 与 Δlabel = label(变体) − label(WT) 的一致性：
     符号准确率、**提拔精度**（模型说"变体能赢"时真赢的比例）、Spearman；
     家族级 cluster bootstrap CI + 家族内置换检验 p 值。
  2. ΔQ 的来源分解：0.4·Δthermo 与 0.6·Δoligo 的方差占比、两者符号一致率。
  3. **需求对照**：用文献先验 λ 算每个变体所需相对提升 1/C−1，与实测 ΔQ/Q(WT) 比较
     —— 直接检验"g1 GU 可被提拔、g12 PP 不可"这两个算例是否在真实数据上成立。

口径（可切换）：全 102 行（项目默认 `mismatch_curated_td.csv`）/ 仅 Ohnishi 78 行；
并把 (α,β) 取 (1,0) 纯热力、(0,1) 纯 DL、(0.4,0.6) 生产值三档对照。

输出：outputs/analysis/mismatch/qguide_variant_signal.json
运行：OligoFormer部分/.venv/Scripts/python.exe scripts/experiment_qguide_variant_signal.py
"""
from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "OligoFormer部分"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from sirna_pipeline.common import mismatch_class  # noqa: E402

# ---- 文献先验（第四方案给出的锚点，未做任何拟合）----
LAMBDA_POS_ANCHORS = {"g1": 0.05, "g2_8_seed": 0.25, "g9_13_central": 0.30,
                      "g14_18": 0.15, "g19": 0.10}
LAMBDA_TYPE = {"WC": 0.0, "GU": 0.5, "PP": 1.0, "YY": 0.8, "PY": 0.8}
SEED = 42
N_BOOT = 2000
N_PERM = 1000


def lambda_position(pos: int) -> float:
    """按锚点分带取 λ_position（g1≈0.05、g12≈0.30、g17/18≈0.15、g19≈0.10）。"""
    if pos == 1:
        return LAMBDA_POS_ANCHORS["g1"]
    if 2 <= pos <= 8:
        return LAMBDA_POS_ANCHORS["g2_8_seed"]
    if 9 <= pos <= 13:
        return LAMBDA_POS_ANCHORS["g9_13_central"]
    if 14 <= pos <= 18:
        return LAMBDA_POS_ANCHORS["g14_18"]
    return LAMBDA_POS_ANCHORS["g19"]


def minmax(vals):
    """生产 ranker 的归一化口径：min-max；span=0 时全取 0.5。"""
    lo, hi = min(vals), max(vals)
    span = hi - lo
    return [0.5] * len(vals) if span == 0 else [(v - lo) / span for v in vals]


def thermo_scores(rows):
    """复用管道 ranker 默认权重/方向（与生产同一事实来源）。"""
    import run_experiments as rx  # noqa: PLC0415

    scores, method = rx._pipeline_thermo_score(rows)
    if scores is None:
        raise SystemExit("热力分计算失败：%s" % method)
    return scores, method


def oligo_raw(df: pd.DataFrame, device: str, batch: int) -> list[float]:
    """OligoFormer + RNA-FM 直接推理，返回原始 efficacy（与线上同一口径）。"""
    from sirna_pipeline.stages.oligoformer import executor as ex  # noqa: PLC0415

    import re
    norm = lambda s: re.sub(r"\s+", "", str(s).upper()).replace("T", "U")  # noqa: E731
    ckpt = ex._find_ckpt(REPO, None)
    if ckpt is None:
        raise SystemExit("未找到 RNA-FM 权重")
    work = pd.DataFrame({"siRNA": df["siRNA"].map(norm), "mRNA": df["mRNA"].map(norm)})
    work = work.astype(object)
    sys.path.insert(0, str(REPO / "scripts"))
    import infer as infer_mod  # noqa: PLC0415

    work = infer_mod.calculate_td(work).reset_index(drop=True)
    work["siRNA"] = work["siRNA"].astype(str)
    work["mRNA"] = work["mRNA"].astype(str)

    uniq = {"siRNA": {}, "mRNA": {}}
    for s in work["siRNA"]:
        uniq["siRNA"].setdefault(ex._md5seq(s), s)
    for s in work["mRNA"]:
        uniq["mRNA"].setdefault(ex._md5seq(s), s)
    cache = REPO / "data" / "RNAFM_cache"
    ex._embed_repo_fm(uniq, cache, REPO, ckpt, seed=SEED, batch_size=batch,
                      log=lambda *a: print("[fm]", *a, flush=True), device=device)
    token = "qguide_%d" % os.getpid()
    data_dir = REPO / "data" / "infer" / token
    shutil.rmtree(data_dir, ignore_errors=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    ex._materialize(data_dir, cache, uniq, log=lambda *a: None)
    try:
        scores = ex._infer(work, token, REPO, REPO / "model" / "best_model.pth", SEED,
                           device=device)
    finally:
        shutil.rmtree(data_dir, ignore_errors=True)
    return [float(v) for v in scores]


def cmatch(positions: list[int], types: list[str]) -> tuple[float, list[dict]]:
    """C_match = Π(1 − λ_p·λ_t)；返回 (C, 每个错配的明细)。"""
    c, detail = 1.0, []
    for pos, t in zip(positions, types):
        cls = mismatch_class.mismatch_class(*t.split(":")) if ":" in t else "PY"
        lp = lambda_position(int(pos))
        lt = LAMBDA_TYPE.get(cls, LAMBDA_TYPE["PY"])
        factor = max(0.0, 1.0 - lp * lt)
        c *= factor
        detail.append({"pos": int(pos), "type": t, "class": cls,
                       "lambda_pos": lp, "lambda_type": lt, "factor": round(factor, 4)})
    return c, detail


def spearman(a, b) -> float | None:
    sa, sb = pd.Series(list(a), dtype=float), pd.Series(list(b), dtype=float)
    if len(sa) < 3 or sa.nunique() < 2 or sb.nunique() < 2:
        return None
    return float(sa.corr(sb, method="spearman"))


def family_bootstrap(pairs: list[dict], fn, rng: random.Random) -> dict:
    """家族级 cluster bootstrap（变体在同一家族内相关，不能按行重采样）。"""
    fams: dict[str, list[dict]] = {}
    for p in pairs:
        fams.setdefault(p["family"], []).append(p)
    keys = list(fams)
    vals = []
    for _ in range(N_BOOT):
        draw = []
        for _ in range(len(keys)):
            draw.extend(fams[rng.choice(keys)])
        v = fn(draw)
        if v is not None:
            vals.append(v)
    if not vals:
        return {"lo": None, "hi": None, "n_valid": 0}
    vals.sort()
    return {"lo": vals[int(0.025 * len(vals))], "hi": vals[int(0.975 * len(vals)) - 1],
            "n_valid": len(vals), "n_families": len(keys)}


def core_metrics(pairs: list[dict]) -> dict:
    """ΔQ vs Δlabel 的一致性指标。提拔精度是唯一体现"突变体可能更好"的指标。"""
    if not pairs:
        return {}
    dq = [p["dQ"] for p in pairs]
    dl = [p["dl"] for p in pairs]
    sign_ok = [1.0 if (a > 0) == (b > 0) else 0.0 for a, b in zip(dq, dl)]
    promoted = [p for p in pairs if p["dQ"] > 0]
    hit = [p for p in promoted if p["dl"] > 0]
    base = sum(1 for x in dl if x > 0) / len(dl)
    return {
        "n_pairs": len(pairs),
        "sign_accuracy": sum(sign_ok) / len(sign_ok),
        "label_better_base_rate": base,
        "balanced_accuracy": None,  # 见 finalize（需要两侧都非空）
        "promoted_n": len(promoted),
        "promote_precision": (len(hit) / len(promoted)) if promoted else None,
        "spearman_dQ_dlabel": spearman(dq, dl),
    }


def finalize(m: dict, pairs: list[dict]) -> dict:
    """补 balanced accuracy（变体更好/更差两侧的召回均值），避免"全判 WT 更好"虚高。"""
    better = [p for p in pairs if p["dl"] > 0]
    worse = [p for p in pairs if p["dl"] <= 0]
    tp = sum(1 for p in better if p["dQ"] > 0) / len(better) if better else None
    tn = sum(1 for p in worse if p["dQ"] <= 0) / len(worse) if worse else None
    m["recall_variant_better"] = tp
    m["recall_wt_better"] = tn
    m["balanced_accuracy"] = (None if tp is None or tn is None else 0.5 * (tp + tn))
    m["n_variant_better"] = len(better)
    m["n_wt_better"] = len(worse)
    return m


def permutation_p(pairs: list[dict], observed: float | None,
                  rng: random.Random) -> float | None:
    """家族内置换 ΔQ 与变体的对应关系（保留家族结构，破坏 ΔQ–Δlabel 对应）。"""
    if observed is None:
        return None
    fams: dict[str, list[dict]] = {}
    for p in pairs:
        fams.setdefault(p["family"], []).append(p)
    count = 0
    for _ in range(N_PERM):
        dq = []
        for f, group in fams.items():
            vals = [p["dQ"] for p in group]
            rng.shuffle(vals)
            dq.extend(vals)
        dl = [p["dl"] for p in pairs]
        acc = sum(1 for a, b in zip(dq, dl) if (a > 0) == (b > 0)) / len(dl)
        if acc >= observed:
            count += 1
    return (count + 1) / (N_PERM + 1)


def build_pairs(df: pd.DataFrame, q: list[float]) -> list[dict]:
    """每个变体与其家族 WT 配对，返回 ΔQ / Δlabel 及错配明细。"""
    import json as _json

    wt_q: dict[str, float] = {}
    for i, (_, r) in enumerate(df.iterrows()):
        if int(r["mismatch_count"]) == 0:
            wt_q[str(r["family"])] = q[i]
    pairs = []
    for i, (_, r) in enumerate(df.iterrows()):
        fam = str(r["family"])
        if int(r["mismatch_count"]) == 0 or fam not in wt_q:
            continue
        pos = [int(x) for x in _json.loads(r["mismatch_positions"])]
        types = list(_json.loads(r["mismatch_types"]))
        c, detail = (cmatch(pos, types) if pos and types else (None, []))
        qw = wt_q[fam]
        pairs.append({
            "row_id": str(r["row_id"]), "dataset": str(r["dataset"]), "family": fam,
            "sirna_id": str(r["sirna_id"]),
            "n_mismatch": int(r["mismatch_count"]), "positions": pos, "types": types,
            "q_variant": q[i], "q_wt": qw,
            "dQ": q[i] - qw, "rel_gain": ((q[i] - qw) / qw) if qw else None,
            "label_variant": float(r["label"]), "label_wt": None,
            "dl": None, "C_match": c, "required_rel_gain": (1.0 / c - 1.0) if c else None,
            "cmatch_detail": detail,
        })
    # 回填 WT 标签
    wt_label = {str(r["family"]): float(r["label"]) for _, r in df.iterrows()
                if int(r["mismatch_count"]) == 0}
    for p in pairs:
        p["label_wt"] = wt_label[p["family"]]
        p["dl"] = p["label_variant"] - p["label_wt"]
    return pairs


def analyse(df: pd.DataFrame, raw_oligo: list[float], alpha: float, beta: float,
            tag: str, rng: random.Random) -> dict:
    """在给定子集上按 (α,β) 组合 Q_guide，输出全部指标。"""
    rows = df.to_dict(orient="records")
    thermo, method = thermo_scores(rows)
    oligo = minmax(raw_oligo)
    q = [alpha * thermo[i] + beta * oligo[i] for i in range(len(df))]

    # 分量口径：纯热力 / 纯 DL（只比较 Δ 的**符号**，故无需再乘 α/β）
    q_thermo, q_dl = thermo, oligo

    pairs = build_pairs(df, q)
    out = {"tag": tag, "n_rows": len(df), "n_pairs": len(pairs),
           "alpha": alpha, "beta": beta, "thermo_method": method}
    for name, series in (("combined", q), ("thermo_only", q_thermo), ("dl_only", q_dl)):
        pr = build_pairs(df, series)
        m = finalize(core_metrics(pr), pr)
        m["sign_accuracy_perm_p"] = permutation_p(pr, m.get("sign_accuracy"), rng)
        m["sign_accuracy_ci"] = family_bootstrap(
            pr, lambda ps: (sum(1 for p in ps if (p["dQ"] > 0) == (p["dl"] > 0)) / len(ps))
            if ps else None, rng)
        if name == "combined":
            out["combined"] = m
            out["_pairs"] = pairs
        else:
            out[name] = m

    # 方差分解：ΔQ 由哪一部分驱动
    import statistics as st

    pair_th = build_pairs(df, q_thermo)
    pair_dl = build_pairs(df, q_dl)
    vth = st.pvariance([p["dQ"] for p in pair_th]) if pair_th else 0.0
    vdl = st.pvariance([p["dQ"] for p in pair_dl]) if pair_dl else 0.0
    out["variance_share"] = {
        "var_dThermo": vth, "var_dOligo": vdl,
        "weighted_thermo": (alpha ** 2) * vth, "weighted_dl": (beta ** 2) * vdl,
        "dl_share_of_weighted_variance": ((beta ** 2) * vdl /
                                          (((alpha ** 2) * vth) + ((beta ** 2) * vdl)))
        if (((alpha ** 2) * vth) + ((beta ** 2) * vdl)) > 0 else None,
    }
    out["sign_agreement_dThermo_dOligo"] = (
        sum(1 for a, b in zip([p["dQ"] for p in pair_th], [p["dQ"] for p in pair_dl])
            if (a > 0) == (b > 0)) / len(pair_th) if pair_th else None)
    out["spearman_dThermo_dOligo"] = spearman([p["dQ"] for p in pair_th],
                                              [p["dQ"] for p in pair_dl])

    # 需求对照：实测相对提升 vs 1/C − 1
    typed = [p for p in pairs if p["C_match"] is not None and p["rel_gain"] is not None]
    promoted = [p for p in typed if p["dQ"] > 0 and p["rel_gain"] >= p["required_rel_gain"]]
    out["required_gain_test"] = {
        "n_typed_pairs": len(typed),
        "n_untyped": len(pairs) - len(typed),
        "n_pass_requirement": len(promoted),
        "pass_row_ids": [p["row_id"] for p in promoted],
        "median_required_rel_gain": (sorted(p["required_rel_gain"] for p in typed)[len(typed) // 2]
                                     if typed else None),
        "median_rel_gain": (sorted(p["rel_gain"] for p in typed)[len(typed) // 2]
                            if typed else None),
        "by_class": {},
    }
    # 按 (位置带, 类型) 分组统计——直接检验 g1/GU 与 g12/PP 两个算例
    buckets: dict[str, list[dict]] = {}
    for p in typed:
        for d in p["cmatch_detail"]:
            band = ("g1" if d["pos"] == 1 else "g2-8" if d["pos"] <= 8 else
                    "g9-13" if d["pos"] <= 13 else "g14-18" if d["pos"] <= 18 else "g19")
            buckets.setdefault("%s|%s" % (band, d["class"]), []).append(p)
    for k, group in sorted(buckets.items()):
        out["required_gain_test"]["by_class"][k] = {
            "n": len(group),
            "C_median": sorted(p["C_match"] for p in group)[len(group) // 2],
            "required_rel_gain_median": sorted(p["required_rel_gain"] for p in group)[len(group) // 2],
            "observed_rel_gain_median": sorted(p["rel_gain"] for p in group)[len(group) // 2],
            "n_pass": sum(1 for p in group
                          if p["dQ"] > 0 and p["rel_gain"] >= p["required_rel_gain"]),
            "label_better_frac": sum(1 for p in group if p["dl"] > 0) / len(group),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--curated", type=Path,
                    default=ROOT / "outputs" / "analysis" / "mismatch" /
                            "mismatch_curated_td.csv")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "outputs" / "analysis" / "mismatch" /
                            "qguide_variant_signal.json")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--batch", type=int, default=32)
    args = ap.parse_args()

    rng = random.Random(SEED)
    data = pd.read_csv(args.curated, dtype=str, encoding="utf-8-sig").fillna("")
    print("载入 %d 行；家族 %d 个" % (len(data), data["family"].nunique()))

    print("OligoFormer 推理中（与线上同一口径）...")
    raw = oligo_raw(data, args.device, args.batch)
    data = data.reset_index(drop=True)

    report: dict = {
        "experiment": "qguide_variant_signal",
        "purpose": "检验 Q_guide 对单点突变的变化方向是否可信（乘性双信号分解的前提）",
        "lambda_priors": {"position_anchors": LAMBDA_POS_ANCHORS, "type": LAMBDA_TYPE,
                          "note": "文献先验，未拟合；多变体按 Π(1−λp·λt) 连乘"},
        "seed": SEED, "n_boot": N_BOOT, "n_perm": N_PERM,
        "scopes": {},
    }

    scopes = {
        "all_102": data,
        "ohnishi_78": data[data["dataset"] == "Ohnishi_2008"].reset_index(drop=True),
    }
    for scope_name, sub in scopes.items():
        pos = data.index.get_indexer_for(sub.index) if scope_name != "all_102" else None
        sub_raw = [raw[i] for i in pos] if pos is not None else raw
        entry = {}
        for tag, (a, b) in {"alpha1.0_beta0.0": (1.0, 0.0),
                            "alpha0.0_beta1.0": (0.0, 1.0),
                            "alpha0.4_beta0.6": (0.4, 0.6)}.items():
            res = analyse(sub.reset_index(drop=True), sub_raw, a, b,
                          "%s|%s" % (scope_name, tag), rng)
            pairs = res.pop("_pairs")
            entry[tag] = res
            if tag == "alpha0.4_beta0.6":
                report["scopes"][scope_name] = entry
                report["scopes"][scope_name]["_pairs_head"] = [
                    {k: v for k, v in p.items() if k != "cmatch_detail"} for p in pairs]
        # 打印摘要
        comb = entry["alpha0.4_beta0.6"]["combined"]
        print("\n=== %s ===" % scope_name)
        print("  变体对数 %d（其中标签上「变体更好」%d 对，基准率 %.1f%%）"
              % (comb["n_pairs"], comb["n_variant_better"],
                 100 * comb["label_better_base_rate"]))
        for name in ("combined", "thermo_only", "dl_only"):
            m = entry["alpha0.4_beta0.6"][name]
            print("  %-12s 符号准确率 %.3f [%.3f,%.3f] p=%.3f | balanced %.3f "
                  "| 提拔 %s 个，精度 %s | ρ=%.3f"
                  % (name, m["sign_accuracy"], m["sign_accuracy_ci"]["lo"],
                     m["sign_accuracy_ci"]["hi"], m["sign_accuracy_perm_p"] or float("nan"),
                     m["balanced_accuracy"] if m["balanced_accuracy"] is not None else float("nan"),
                     m["promoted_n"],
                     ("%.2f" % m["promote_precision"]) if m["promote_precision"] is not None else "—",
                     m["spearman_dQ_dlabel"] if m["spearman_dQ_dlabel"] is not None else float("nan")))
        vs = entry["alpha0.4_beta0.6"]["variance_share"]
        print("  ΔQ 方差占比：DL %.1f%% / 热力 %.1f%%；两者符号一致率 %.2f，ρ=%.2f"
              % (100 * (vs["dl_share_of_weighted_variance"] or 0),
                 100 * (1 - (vs["dl_share_of_weighted_variance"] or 0)),
                 entry["alpha0.4_beta0.6"]["sign_agreement_dThermo_dOligo"],
                 entry["alpha0.4_beta0.6"]["spearman_dThermo_dOligo"] or float("nan")))
        rg = entry["alpha0.4_beta0.6"]["required_gain_test"]
        print("  需求对照：可定型 %d 对（未定型 %d）；达标 %d 个"
              % (rg["n_typed_pairs"], rg["n_untyped"], rg["n_pass_requirement"]))
        print("    需要提升中位 %.3f vs 实测提升中位 %.3f"
              % (rg["median_required_rel_gain"] or float("nan"),
                 rg["median_rel_gain"] or float("nan")))
        for k, v in rg["by_class"].items():
            print("    %-10s n=%-3d C中位 %.3f 需提升 %.3f 实测 %.3f 达标 %d 标签更优占比 %.2f"
                  % (k, v["n"], v["C_median"], v["required_rel_gain_median"],
                     v["observed_rel_gain_median"], v["n_pass"], v["label_better_frac"]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nwritten:", args.out)


if __name__ == "__main__":
    main()
