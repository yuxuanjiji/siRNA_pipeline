# -*- coding: utf-8 -*-
"""按**现行生产口径**重算阳性对照排名，并重画回收图（含全部 5 条可比对照）。

背景：`outputs/analysis/task17_positive_controls.csv` 生成于 2026-09-11 10:02，
早于重定权重（α0.4/β0.6）与变体层；其 `final_rank` 与现行交付 `outputs/results/rank_final.csv`
不一致（例：17-4b 记为 rank 3，现行口径为 rank 75）。本脚本：

  1. 以 guide 序列在**现行** `outputs/results/rank_final.csv` 中回查排名/分数，
     写入新列 `final_rank_current` / `rank_percentile_current` / `current_score` /
     `current_window_id` / `current_variant_id` / `rank_source`；
  2. **保留**原列并改名为 `final_rank_legacy_20260911` / `rank_percentile_legacy`（审计留痕，不删）；
  3. 重画回收图：**全部 5 条可比对照** + 随机期望基线 + 显著性检验（Irwin–Hall 平均百分位），
     并列出 2 条未匹配对照及原因。

产出：
  outputs/analysis/task17_positive_controls.csv         （原地更新，列只增不改语义）
  outputs/analysis/task17_positive_controls.json        （同步 matched/unmatched 与口径说明）
  outputs/analysis/figures/fig_positive_controls_recovery.png/pdf

用法：python scripts/update_positive_controls.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np              # noqa: E402
import pandas as pd             # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
FIG = AN / "figures"
CSV = AN / "task17_positive_controls.csv"
JS = AN / "task17_positive_controls.json"
CURRENT = ROOT / "outputs" / "results" / "rank_final.csv"
LEGACY_CUTOFF = "2026-09-11"

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 200
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

C_TOP, C_MID, C_BASE = "#2E7D32", "#4C78A8", "#90A4AE"


def irwin_hall_cdf(s: float, n: int = 5) -> float:
    """n 个 U(0,1) 之和 ≤ s 的概率（s ≤ 1 段）：s^n / n!。"""
    return (s ** n) / float(math.factorial(n)) if s <= 1 else float("nan")


def main() -> int:
    ctrl = pd.read_csv(CSV, dtype=str, encoding="utf-8-sig")
    cur = pd.read_csv(CURRENT, dtype=str, encoding="utf-8-sig")
    cur = cur[cur["final_rank"].notna() & (cur["final_rank"] != "")].copy()
    n_total = len(cur)
    by_guide = {str(g): r for g, r in zip(cur["guide_checked"], cur.to_dict("records"))}

    # 旧口径列只搬一次：仅在原列仍存在时改名为审计列，重跑不覆盖已保留的旧值
    leg_rank_col = "final_rank_legacy_%s" % LEGACY_CUTOFF.replace("-", "")
    if "final_rank" in ctrl.columns:
        ctrl[leg_rank_col] = ctrl["final_rank"]
        if "rank_percentile" in ctrl.columns:
            ctrl["rank_percentile_legacy"] = ctrl["rank_percentile"]
        ctrl = ctrl.drop(columns=[c for c in ("final_rank", "rank_percentile") if c in ctrl.columns])

    ranks, pcts, scores, wins, vars_ = [], [], [], [], []
    for _, r in ctrl.iterrows():
        rec = by_guide.get(str(r.get("guide19", ""))) if str(r.get("matched")) == "1" else None
        if rec is None:
            ranks.append(""); pcts.append(""); scores.append(""); wins.append(""); vars_.append("")
            continue
        rk = int(rec["final_rank"])
        ranks.append(str(rk))
        pcts.append("%.6f" % (rk / n_total))
        scores.append(rec.get("final_score", ""))
        wins.append(rec.get("window_id", ""))
        vars_.append(rec.get("variant_id", ""))
    ctrl["final_rank_current"] = ranks
    ctrl["rank_percentile_current"] = pcts
    ctrl["current_score"] = scores
    ctrl["current_window_id"] = wins
    ctrl["current_variant_id"] = vars_
    ctrl["rank_source"] = ("outputs/results/rank_final.csv（现行口径：α0.4/β0.6 + 变体层，"
                           "n=%d）" % n_total)
    ctrl.to_csv(CSV, index=False, encoding="utf-8-sig")

    matched = ctrl[ctrl["matched"] == "1"].copy()
    matched["pct"] = matched["rank_percentile_current"].astype(float)
    matched = matched.sort_values("pct")
    pcts_sorted = matched["pct"].tolist()
    mean_pct = float(np.mean(pcts_sorted))
    p_val = irwin_hall_cdf(sum(pcts_sorted), len(pcts_sorted))
    in10 = int((matched["pct"] <= 0.10).sum())
    in25 = int((matched["pct"] <= 0.25).sum())
    unmatched = ctrl[ctrl["matched"] != "1"]

    # ---- 同步 JSON ----
    js = json.loads(JS.read_text(encoding="utf-8"))
    js.update({
        "controls": int(len(ctrl)), "matched": int(len(matched)),
        "unmatched": unmatched["id"].tolist(),
        "rank_source": "outputs/results/rank_final.csv（现行口径 α0.4/β0.6 + 变体层，n=%d）" % n_total,
        "legacy_note": ("原 final_rank/rank_percentile（%s 口径）已保留为 "
                        "final_rank_legacy_20260911 / rank_percentile_legacy，未删除"
                        % LEGACY_CUTOFF),
        "current_ranks": {r["id"]: {"rank": int(r["final_rank_current"]),
                                    "percentile": round(float(r["rank_percentile_current"]), 6),
                                    "source": r["source"]}
                          for _, r in matched.iterrows()},
        "summary_current": {"n_matched": int(len(matched)), "in_top10pct": in10,
                            "in_top25pct": in25, "mean_percentile": round(mean_pct, 4),
                            "p_mean_percentile_irwin_hall": float("%.2e" % p_val)},
        "caveat": ("Suzuki 2008 的 17-4a/4b/4c 原文按三连混合给药，个体效力未分别刊出；"
                   "17-3/17-5 在 CDS 内无精确匹配（可能靶向 UTR/异构体），无法评分"),
    })
    JS.write_text(json.dumps(js, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- 重画图（全部可比对照 + 随机基线 + 显著性）----
    labels = ["%s  %s" % (r["id"], str(r["source"]).split(",")[0]) for _, r in matched.iterrows()]
    vals = [p * 100 for p in matched["pct"]]
    colors = [C_TOP if v <= 5 else (C_MID if v <= 25 else C_BASE) for v in vals]
    fig, ax = plt.subplots(figsize=(11.6, 4.6))
    ys = np.arange(len(vals))[::-1]
    ax.barh(ys, vals, color=colors, height=0.55)
    for y, v, (_, r) in zip(ys, vals, matched.iterrows()):
        ax.text(v + 0.5, y, "rank %s / %d（前 %.1f%%）" %
                (r["final_rank_current"], n_total, v), va="center", fontsize=9.5)
    ax.axvline(50, color="#C62828", ls="--", lw=1.6)
    ax.text(50.8, float(np.mean(ys)), "随机期望\n（均匀零假设）= 50%",
            color="#C62828", fontsize=9.5, va="center")
    ax.set_yticks(ys, labels, fontsize=9.6)
    ax.set_xlim(0, 62)
    ax.set_xlabel("在 %d 条参与排序的 SFRP1 候选中的排名百分位（越左越好；"
                  "结构过滤后共 4417 行，其中 4218 行可排序）" % n_total)
    ax.set_title("文献已发表有效 siRNA 的回收排名（现行口径 α0.4/β0.6 + 变体层）\n"
                 "%d 条可比对照中 %d 条进入前 25%%、%d 条进入前 10%%；"
                 "平均百分位 %.1f%%，p≈%.1e（Irwin–Hall，n=%d）"
                 % (len(matched), in25, in10, mean_pct * 100, p_val, len(matched)),
                 fontsize=11.5, loc="left")
    if len(unmatched):
        ax.text(0.5, -0.30, "未纳入统计的对照：" +
                "；".join("%s（%s）" % (r["id"], str(r["ctx_note"])[:38])
                          for _, r in unmatched.iterrows()),
                transform=ax.transAxes, fontsize=8.6, color="#546E7A")
    ax.text(0.5, -0.42, "注：Suzuki 2008 的 17-4a/4b/4c 原文按三连混合给药，个体效力未分别刊出，"
                        "仅作「该区域可被命中」的证据。",
            transform=ax.transAxes, fontsize=8.6, color="#546E7A")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG / ("fig_positive_controls_recovery.%s" % ext), bbox_inches="tight")
    plt.close(fig)

    print("CSV 已更新：", CSV)
    print(" 现行排名：", {r["id"]: r["final_rank_current"] for _, r in matched.iterrows()})
    print(" 前10%%: %d | 前25%%: %d | 平均百分位 %.4f | p=%.2e" % (in10, in25, mean_pct, p_val))
    print("图已重画：", FIG / "fig_positive_controls_recovery.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
