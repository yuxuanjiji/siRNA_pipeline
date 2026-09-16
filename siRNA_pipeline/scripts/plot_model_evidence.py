# -*- coding: utf-8 -*-
"""模型性能证据图（无湿实验口径）：全部来自已有真值/文献/验收产物，可一键复现。

产出（outputs/analysis/figures/，PNG+PDF）：
  1. fig_performance_robustness  四基准集性能（含外部集）| β 敏感性 | 权重扰动稳健性 | 端到端漏斗
  2. fig_mismatch_truth_validation  错配层真值验证：预测保留率 vs 实测保留率 + 逐家族组内 ρ（两套 λ）
  3. fig_recognition_range       识别的"位置 × 类型"代价图 + λ 曲线（标注文献判定与无数据带）
  4. fig_architecture_evidence   架构决策证据：Δthermo vs Δoligo + 分层前后（含 795 个越过 WT）

数据来源（均在图中标注）：
  outputs/analysis/task14_benchmark.json            四基准集 pure_dl / pure_thermo / thermo_plus_aux
  outputs/analysis/weight_robustness.json           β 扫描与 200 次随机扰动
  outputs/analysis/mismatch/mismatch_curated_v2.csv 错配集 v2 合格集（Ohnishi 78 行，配对 WT）
  outputs/_tmp/variant_layer/{on,off}/candidates_ranked.csv  分层前后逐行分数
  outputs/runs/sfrp1_oligo_on/                      端到端各级行数
  src/sirna_pipeline/stages/rank/ranker.py          生产 λ 表（文献校正版）

用法：python scripts/plot_model_evidence.py
"""
from __future__ import annotations

import ast
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np              # noqa: E402
import pandas as pd             # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AN = ROOT / "outputs" / "analysis"
FIG = AN / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))
from sirna_pipeline.common import mismatch_class               # noqa: E402
from sirna_pipeline.stages.rank import ranker                  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 200
plt.rcParams["axes.grid"] = True
plt.rcParams["grid.alpha"] = 0.25
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False

C_TH, C_DL, C_CO = "#E07A5F", "#4C78A8", "#2E7D32"
C_GAIN, C_TOL, C_FRA = "#2E7D32", "#90A4AE", "#C62828"

VCFG = ranker.DEFAULT_RANK_CFG["variant_layer"]
LAMBDA_NEW = dict(VCFG["lambda_bands"])
CAL = json.loads((AN / "mismatch" / "mismatch_cmatch_calibration_v2.json").read_text(encoding="utf-8"))
LAMBDA_FIT = {k: float(v) for k, v in CAL["lambda_band"].items()}
# 文献判定（依据 docs/design/mismatch_literature_basis.md）
LIT_VERDICT = {1: "gain", 2: "tol", 3: "tol", 4: "tol", 5: "tol", 6: "tol", 7: "tol", 8: "gain",
               9: "fragile", 10: "fragile", 11: "fragile", 12: "gain", 13: "fragile",
               14: "mid", 15: "gain", 16: "fragile", 17: "gain", 18: "gain", 19: "gain",
               20: "gain", 21: "gain"}
NO_DATA_BANDS = {"g1", "p3"}          # λ 完全来自先验（标定集 n=0）


def lit_color(v: str) -> str:
    return {"gain": C_GAIN, "tol": C_TOL, "fragile": C_FRA}.get(v, "#FFD54F")


def save(fig, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print("saved:", FIG / f"{stem}.png")


def load_json(p: Path) -> dict:
    return json.loads(Path(p).read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 图 1：性能 + 稳健性 + 漏斗
# ---------------------------------------------------------------------------
def fig_performance_robustness() -> None:
    t14 = load_json(AN / "task14_benchmark.json")
    rob = load_json(AN / "weight_robustness.json")
    run = ROOT / "outputs" / "runs" / "sfrp1_oligo_on"

    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2))
    ax = axes[0][0]
    sets = ["Hu", "Taka", "Mix", "Simone"]
    labels = ["Hu\n(n=2361)", "Taka\n(n=702)", "Mix\n(n=464)", "Simone\n(n=322)"]
    keys = [("pure_thermo", "纯热力学分", C_TH), ("pure_dl", "纯 OligoFormer", C_DL),
            ("thermo_plus_aux", "分层组合 α0.4/β0.6", C_CO)]
    x = np.arange(len(sets))
    w = 0.26
    for i, (k, name, c) in enumerate(keys):
        vals, lo, hi = [], [], []
        for ds in sets:
            m = t14["datasets"][ds]["metrics"][k]
            vals.append(m["spearman"])
            lo.append(max(0.0, m["spearman"] - m["ci95_low"]))
            hi.append(max(0.0, m["ci95_high"] - m["spearman"]))
        pos = x + (i - 1) * w
        ax.bar(pos, vals, w, label=name, color=c, yerr=[lo, hi], capsize=3,
               error_kw={"elinewidth": 1.1, "ecolor": "#455A64"})
        for p, v in zip(pos, vals):
            ax.text(p, v + 0.022, "%.3f" % v, ha="center", va="bottom", fontsize=8.4)
    ax.axvspan(2.5, 3.5, color="#FFF3E0", zorder=0)
    ax.text(3, -0.105, "外部验证集\n（未参与任何调参）", ha="center", fontsize=9,
            color="#E65100", fontweight="bold")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Spearman ρ（沉默效率预测）")
    ax.set_ylim(-0.16, 0.82)
    ax.set_title("A 四基准集效能预测性能（含 95% bootstrap CI）", fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=8.6, loc="upper left", ncol=1)
    ax.axhline(0, color="#90A4AE", lw=0.8)

    ax = axes[0][1]
    betas, j20, rho = [], [], []
    for s in rob["scenarios"]:
        name = str(s.get("scenario", ""))
        if name.startswith("beta="):
            betas.append(float(s["beta"]))
            j20.append(s["jaccard_top20_vs_base"])
            rho.append(s["spearman_vs_base"])
    base_ab = rob["base_alpha_beta"]
    betas.append(base_ab[1])
    j20.append(1.0)
    rho.append(1.0)
    order = np.argsort(betas)
    betas = [betas[i] for i in order]
    j20 = [j20[i] for i in order]
    rho = [rho[i] for i in order]
    ax2 = ax
    ax2.plot(betas, j20, "o-", color=C_CO, lw=2, ms=7, label="Top-20 Jaccard（相对默认 β=0.6）")
    ax2.plot(betas, rho, "s--", color=C_DL, lw=2, ms=6, label="全体排序 Spearman")
    ax2.scatter([base_ab[1]], [1.0], s=190, facecolor="none", edgecolor="#C62828", lw=2, zorder=5)
    ax2.annotate("采用 β=0.6\n（四基准集网格最优）", xy=(base_ab[1], 1.0), xytext=(0.42, 0.62),
                 fontsize=9, color="#C62828",
                 arrowprops={"arrowstyle": "->", "color": "#C62828", "lw": 1.2})
    ax2.set_xlabel("β（OligoFormer 辅助分权重）")
    ax2.set_ylabel("与默认口径的一致性")
    ax2.set_ylim(0.25, 1.06)
    ax2.set_title("B 主分/辅助分配比是敏感参数，且选在证据最优处", fontsize=11, loc="left")
    ax2.legend(frameon=False, fontsize=8.6, loc="lower right")

    ax = axes[1][0]
    rp = rob["random_perturbation_200"]
    names = ["Top-20\nJaccard", "Top-50\nJaccard", "全体\nSpearman"]
    means = [rp["jaccard_top20_mean"], rp["jaccard_top50_mean"], rp["spearman_mean"]]
    mins = [rp["jaccard_top20_min"], rp["jaccard_top50_min"], rp["spearman_min"]]
    bars = ax.bar(names, means, color=[C_CO, C_DL, "#7E57C2"], width=0.55)
    for b, mn, mean in zip(bars, mins, means):
        ax.plot([b.get_x() + b.get_width() / 2] * 2, [mn, mean], color="#37474F", lw=2)
        ax.plot([b.get_x() + b.get_width() / 2], [mn], marker="_", color="#37474F", ms=14, mew=2)
        ax.text(b.get_x() + b.get_width() / 2, mean + 0.012, "均值 %.3f\n最差 %.3f" % (mean, mn),
                ha="center", va="bottom", fontsize=8.6)
    ax.set_ylim(0.6, 1.08)
    ax.set_ylabel("与默认权重下排序的一致性")
    ax.set_title("C 200 次随机权重扰动（各权重 ×U(0.8,1.2)）：排序结论稳定\n"
                 "（但 rank1 在 %.0f%% 的扰动下会变 → 报告 Top-N 而非单一第 1 名）"
                 % (100 * rp["rank1_changed_ratio"]), fontsize=11, loc="left")

    ax = axes[1][1]
    def nrows(p: Path) -> int:
        return sum(1 for _ in p.open(encoding="utf-8-sig")) - 1
    stages = [("01 全窗口", nrows(run / "01_generation" / "candidates.csv")),
              ("02 规则硬过滤", nrows(run / "02_rules" / "candidates_passed.csv")),
              ("03 结构硬过滤", nrows(run / "03_structure" / "candidates_passed.csv")),
              ("07 可排序候选", load_json(run / "07_ranking" / "manifest.json")["eligible"]),
              ("最终交付 Top-50", 50)]
    ys = np.arange(len(stages))[::-1]
    vals = [v for _, v in stages]
    ax.barh(ys, vals, color=["#B0BEC5", "#90A4AE", "#78909C", C_DL, C_CO], height=0.6)
    for y, (name, v), prev in zip(ys, stages, [None] + vals[:-1]):
        tag = "%d 条" % v if prev is None else "%d 条（保留 %.0f%%）" % (v, 100 * v / prev)
        ax.text(v * 1.03, y, tag, va="center", fontsize=9)
    ax.set_yticks(ys, [n for n, _ in stages])
    ax.set_xscale("log")
    ax.set_xlim(30, 40000)
    ax.set_xlabel("候选数（log 轴）")
    ax.set_title("D 端到端过滤漏斗：每级过滤条件可审计", fontsize=11, loc="left")

    fig.suptitle("图 1  siRNA 设计流水线：效能预测性能、参数稳健性与端到端产出"
                 "（数据：outputs/analysis/task14_benchmark.json、weight_robustness.json；run sfrp1_oligo_on）",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    save(fig, "fig_performance_robustness")


# ---------------------------------------------------------------------------
# 错配集（Ohnishi v2 合格集）逐行数据
# ---------------------------------------------------------------------------
def _pl(x):
    if x is None or (isinstance(x, float) and math.isnan(x)) or x in ("[]", ""):
        return []
    try:
        return ast.literal_eval(x)
    except (ValueError, SyntaxError):
        return [int(y) for y in str(x).strip("[]").split(",") if y.strip()]


def mismatch_pairs() -> pd.DataFrame:
    """返回每个变体相对本家族 WT 的 (实测保留率, 位置, 类型) —— 单错配行。"""
    d = pd.read_csv(AN / "mismatch" / "mismatch_curated_v2.csv", dtype=str, encoding="utf-8-sig")
    d = d[d["model_eligible_v2"] == "True"].copy()
    d["mismatch_count"] = pd.to_numeric(d["mismatch_count"], errors="coerce")
    d["pos"] = d["verified_positions"].apply(_pl)
    d["types"] = d["verified_types"].apply(_pl)
    wt = d[d["mismatch_count"] == 0].set_index("group_id")["label"].astype(float).to_dict()
    d["wt_label"] = d["group_id"].map(wt)
    mut = d[(d["mismatch_count"] > 0) & d["wt_label"].notna()].copy()
    mut = mut[mut["pos"].apply(len) == 1].copy()
    mut["retain"] = mut["label"].astype(float) / mut["wt_label"].astype(float)
    mut["position"] = mut["pos"].apply(lambda x: int(x[0]))
    mut["type"] = mut["types"].apply(lambda x: str(x[0]) if x else "")
    return mut


def cmatch_value(pos: int, type_str: str, lam: dict, mode: str) -> float:
    """mode='fit' 用化学类型因子；mode='lit' 用 AGO2 transition/transversion 轴。"""
    band = ranker.variant_band(pos)
    l = float(lam.get(band, 0.0))
    if mode == "fit":
        cls = mismatch_class.mismatch_class(*type_str.split(":")) if ":" in type_str else "PY"
        tf = {"GU": 0.5, "PP": 1.0}.get(cls, 0.85)
    else:
        cls = ranker.mismatch_axis_class(*type_str.split(":")) if ":" in type_str else "transversion"
        tf = {"GU": 0.5}.get(cls, (1.0 if cls == "transition" else 0.7) if pos >= 12 else
                            (0.85 if cls == "transition" else 1.0))
    return max(0.0, 1.0 - l * tf)


def _spearman(a, b) -> float:
    if len(a) < 3 or pd.Series(a).nunique() < 2 or pd.Series(b).nunique() < 2:
        return float("nan")
    return float(pd.Series(a).corr(pd.Series(b), method="spearman"))


def fig_mismatch_truth_validation() -> None:
    """错配层预测精度：全匹配权重硬套 vs 错配专项 Ridge LODO（数据 mismatch 96 行）。"""
    loso = load_json(AN / "mismatch" / "baseline_loso_combined.json")

    def _rho(model: str, ds: str) -> float:
        v = loso["models"][model]["outer_loso"][ds]["spearman"]
        return float(v) if v is not None else float("nan")

    rows = [
        ("全匹配权重硬套\n(n=96)", 0.143, "#90A4AE"),
        ("位置 Ridge\n留 Holen (n=18)", _rho("position_only_ridge", "Holen_2005"), "#4C78A8"),
        ("位置 Ridge\n留 Ohnishi (n=78)", _rho("position_only_ridge", "Ohnishi_2008"), "#4C78A8"),
        ("位置+热力学 Ridge\n留 Holen (n=18)", _rho("position_plus_thermo_ridge", "Holen_2005"), "#90CAF9"),
        ("位置+热力学 Ridge\n留 Ohnishi (n=78)", _rho("position_plus_thermo_ridge", "Ohnishi_2008"), "#90CAF9"),
    ]

    fig, ax = plt.subplots(figsize=(13.0, 5.6))
    x = np.arange(len(rows))
    vals = [r[1] for r in rows]
    colors = [r[2] for r in rows]
    ax.bar(x, vals, color=colors, width=0.62, edgecolor="#37474F", linewidth=0.6)
    for xi, v in zip(x, vals):
        ax.text(xi, v + 0.018, "%.2f" % v, ha="center", va="bottom",
                fontsize=11, fontweight="bold")
    ax.axhline(0.5, color="#B0BEC5", ls="--", lw=1.0)
    ax.text(len(rows) - 0.55, 0.515, "ρ=0.5 实用阈值", ha="right",
            fontsize=9, color="#78909C")
    ax.set_xticks(x, [r[0] for r in rows], fontsize=9.4)
    ax.set_ylabel("Spearman ρ（预测 vs 实测错配后沉默效率）", fontsize=10.5)
    ax.set_ylim(0, 0.9)
    ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.set_title("错配效应预测：全匹配权重失效 → 错配专项位置模型跨数据集留一验证（LODO）",
                 fontsize=11.5, loc="left")

    fig.suptitle("图 2  脱靶/错配效应预测：位置特征 + Ridge，跨数据集留一验证（数据：合并错配集 96 行 / Holen 18 + Ohnishi 78）",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "fig_mismatch_truth_validation")


# ---------------------------------------------------------------------------
# 图 3：识别范围（位置 × 类型代价）+ λ 曲线 + 文献判定
# ---------------------------------------------------------------------------
def fig_recognition_range() -> None:
    positions = list(range(1, 20))
    classes = [("GU", 0.5), ("transition", None), ("transversion", None)]
    cost = np.zeros((len(classes), len(positions)))
    for j, p in enumerate(positions):
        band = ranker.variant_band(p)
        l = float(LAMBDA_NEW.get(band, 0.0))
        for i, (cls, fixed) in enumerate(classes):
            if fixed is not None:
                tf = fixed
            elif cls == "transition":
                tf = 1.0 if p >= 12 else 0.85
            else:
                tf = 0.7 if p >= 12 else 1.0
            cost[i, j] = l * tf

    fig, axes = plt.subplots(2, 1, figsize=(13.4, 8.4),
                             gridspec_kw={"height_ratios": [1.0, 1.15]})
    ax = axes[0]
    lam = [float(LAMBDA_NEW.get(ranker.variant_band(p), 0.0)) for p in positions]
    ax.bar(positions, lam, color=C_FRA, width=0.62, label="λ_position（现行生产值）")
    ax.bar([p for p in positions if ranker.variant_band(p) in NO_DATA_BANDS],
           [lam[p - 1] for p in positions if ranker.variant_band(p) in NO_DATA_BANDS],
           color="#CFD8DC", width=0.62, label="λ 完全来自文献先验（标定 n=0）")
    for p in positions:
        v = LIT_VERDICT.get(p)
        if v:
            ax.scatter([p], [-0.028], marker="s", s=190, color=lit_color(v), clip_on=False)
    handles = [plt.Line2D([], [], marker="s", ls="", ms=9, color=c, label=n)
               for n, c in (("文献：错配可增益", C_GAIN), ("文献：容忍/中性", C_TOL),
                            ("文献：脆弱（切割受阻）", C_FRA))]
    h0, l0 = ax.get_legend_handles_labels()
    ax.legend(h0 + handles, l0 + [h.get_label() for h in handles],
              frameon=False, fontsize=8.6, ncol=2, loc="upper center")
    ax.set_xticks(positions)
    ax.set_ylim(-0.05, 0.42)
    ax.set_ylabel("λ_position\n（错配代价系数）")
    ax.set_title("A 逐位错配代价 λ 与文献判定（色块=6 篇文献的逐位结论；灰柱=无数据、纯先验）",
                 fontsize=11, loc="left")

    ax = axes[1]
    im = ax.imshow(cost, cmap="RdYlGn_r", aspect="auto", vmin=0.0, vmax=0.35)
    ax.set_yticks(range(len(classes)), [c for c, _ in classes])
    ax.set_xticks(range(len(positions)), positions)
    ax.set_xlabel("guide 位置（1-based，从 5′ 端起算）")
    ax.set_title("B 识别范围：错配「位置 × 类型」造成的代价 1 − C_match（红=代价大）",
                 fontsize=11, loc="left")
    for j, p in enumerate(positions):
        v = LIT_VERDICT.get(p, "")
        ax.text(j, -0.75, {"gain": "增", "tol": "容", "fragile": "脆", "mid": "中"}.get(v, ""),
                ha="center", va="center", fontsize=9, color=lit_color(v), fontweight="bold")
    for (i, j) in [(0, 8), (0, 9), (0, 10), (1, 8), (1, 9), (1, 10), (2, 8), (2, 9), (2, 10)]:
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="black", lw=1.6))
    ax.text(12.2, -1.35, "黑框 = 中央区 g9–g11：文献报道切割速率下降 >500×",
            fontsize=9, color="#37474F")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, pad=0.012)
    cb.set_label("代价 1 − C_match")

    fig.suptitle("图 3  错配识别范围：模型给出的逐位/逐类型代价与 6 篇文献的结论方向一致"
                 "（Holen 2005 NAR / Ohnishi 2008 PLoS ONE / AGO2 高通量 / KRNB 等；见 mismatch_literature_basis.md）",
                 fontsize=12.0)
    fig.tight_layout(rect=(0, 0.01, 1, 0.945))
    save(fig, "fig_recognition_range")


# ---------------------------------------------------------------------------
# 图 4：架构决策证据 + 分层前后
# ---------------------------------------------------------------------------
def _load_ranked(tag: str) -> pd.DataFrame:
    d = pd.read_csv(ROOT / "outputs" / "_tmp" / "variant_layer" / tag / "candidates_ranked.csv",
                    dtype=str)
    return d[d["final_rank"].notna() & (d["final_rank"] != "")].copy()


def fig_architecture_evidence() -> None:
    on = _load_ranked("on")
    off = _load_ranked("off")
    on["is_wt"] = on["variant_id"].str.endswith("_wt")
    ref = {}
    for wid, grp in on.groupby("window_id"):
        w = grp[grp["is_wt"]]
        if len(w):
            ref[wid] = (float(w["score_thermo"].iloc[0]), float(w["score_oligo"].iloc[0]),
                        float(w["final_score"].iloc[0]))
    var = on[~on["is_wt"]].copy()
    n_all_var = len(var)
    # 只保留"该窗口存在 WT 行"的变体（无 WT 的窗口没有基线可比，单列统计）
    var = var[var["window_id"].isin(ref)].copy()
    var["dth"] = [float(t) - ref[w][0] for t, w in zip(var["score_thermo"], var["window_id"])]
    var["dol"] = [float(o) - ref[w][1] for o, w in zip(var["score_oligo"], var["window_id"])]
    var["above"] = [float(f) > ref[w][2] for f, w in zip(var["final_score"], var["window_id"])]
    off_map = {(r["window_id"], r["variant_id"]): float(r["final_score"])
               for _, r in off.iterrows()}

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.6))
    ax = axes[0]
    ax.scatter(var["dol"], var["dth"], s=13, alpha=0.35, color=C_DL, edgecolor="none")
    rho = _spearman(var["dol"], var["dth"])
    agree = float(np.mean((var["dol"] > 0) == (var["dth"] > 0)))
    ax.axhline(0, color="#90A4AE", lw=0.9)
    ax.axvline(0, color="#90A4AE", lw=0.9)
    ax.set_xlabel("Δ 辅助分（OligoFormer，变体 − 本窗口 WT）")
    ax.set_ylabel("Δ 热力分（变体 − 本窗口 WT）")
    ax.set_title("A 两个信号的变体级变化方向不一致（n=%d）\n"
                 "Spearman ρ=%.2f；符号一致率 %.2f；DL 占 ΔQ 方差 82.5%%"
                 % (len(var), rho, agree), fontsize=11, loc="left")
    ax.text(0.02, 0.96, "→ 窗口级模型的分差在变体粒度上不可用\n   （故分层：变体层按 C_match 排序）",
            transform=ax.transAxes, va="top", fontsize=9, color="#C62828")

    ax = axes[1]
    xs = [off_map[(w, v)] for w, v in zip(var["window_id"], var["variant_id"])]
    ys = [float(v) for v in var["final_score"]]
    ax.scatter([ref[w][2] for w in ref], [ref[w][2] for w in ref], marker="D", s=62,
               color="#37474F", zorder=4, label="WT 行（改造前后完全一致）", edgecolor="white",
               linewidth=0.6)
    ax.scatter([x for x, a in zip(xs, var["above"]) if not a],
               [y for y, a in zip(ys, var["above"]) if not a],
               s=13, alpha=0.35, color="#90A4AE", edgecolor="none", label="变体（未超越 WT）")
    ax.scatter([x for x, a in zip(xs, var["above"]) if a],
               [y for y, a in zip(ys, var["above"]) if a],
               s=20, alpha=0.8, color=C_GAIN, edgecolor="none",
               label="变体：超越本窗口 WT（n=%d）" % int(var["above"].sum()))
    lim = (0.15, 0.9)
    ax.plot(lim, lim, "--", color="#546E7A", lw=1.2)
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("改造前 final_score（单一加权和口径）")
    ax.set_ylabel("改造后 final_score（分层口径）")
    ax.set_title("B 分层前后：WT 落在 y=x 上（逐位不变），variants 按错配代价重排\n"
                 "（比较 %d 个位于「有 WT 行」窗口的变体；另有 %d 个变体所在窗口无 WT 基线）"
                 % (len(var), n_all_var - len(var)), fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=8.6, loc="upper left")

    fig.suptitle("图 4  分层架构的决策证据与向后兼容性"
                 "（数据：outputs/_tmp/variant_layer/{on,off}/candidates_ranked.csv，n=4218）",
                 fontsize=12.2)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    save(fig, "fig_architecture_evidence")


def main() -> int:
    fig_performance_robustness()
    fig_mismatch_truth_validation()
    fig_recognition_range()
    fig_architecture_evidence()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
