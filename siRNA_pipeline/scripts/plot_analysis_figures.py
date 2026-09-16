"""Create publication-style figures from the generated task 14-17 reports."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "outputs" / "analysis"
FIGURES = ANALYSIS / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_json(name: str) -> dict:
    return json.loads((ANALYSIS / name).read_text(encoding="utf-8"))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def annotate_bars(ax, bars, fmt="{:.3f}"):
    for bar in bars:
        value = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.025,
                fmt.format(value), ha="center", va="bottom", fontsize=8)


def plot_task14() -> None:
    data = load_json("task14_benchmark.json")
    names = ["Hu", "Taka", "Mix"]
    metrics = [("spearman", "Spearman $\\rho$"),
               ("pearson", "Pearson $r$"),
               ("r2", "$R^2$")]
    methods = [("pure_dl", "Pure DL", "#276FBF"),
               ("pure_thermo", "Pure thermo", "#E07A5F"),
               ("thermo_plus_aux", "Thermo + aux", "#4C956C")]
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.4), sharey=True)
    for ax, (key, label) in zip(axes, metrics):
        x = np.arange(len(names))
        width = 0.24
        for offset, (method, method_label, color) in enumerate(methods):
            values = []
            for name in names:
                metric = data.get("datasets", {}).get(name, {}).get("metrics", {}).get(method)
                values.append(metric.get(key) if metric else np.nan)
            plot_values = np.nan_to_num(values, nan=0.0)
            bars = ax.bar(x + (offset - 1) * width, plot_values, width,
                          label=method_label if key == "spearman" else None,
                          color=color, edgecolor="#30343B", linewidth=0.5)
            for bar, value in zip(bars, values):
                if np.isfinite(value):
                    ax.text(bar.get_x() + bar.get_width() / 2, value + 0.018,
                            f"{value:.3f}", ha="center", fontsize=6.5)
        ax.set_title(label)
        ax.set_ylim(-0.1, 0.8)
        ax.set_xticks(x, names)
        ax.grid(axis="y", color="#D9DEE5", linewidth=0.6)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Metric value")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 0.99))
    fig.suptitle("Task 14. OligoFormer benchmark performance", y=1.06, fontweight="bold")
    fig.text(0.5, -0.02, "Three methods evaluated on the same sequence-keyed records; Simone omitted.",
             ha="center", fontsize=8, color="#555B64")
    save(fig, "task14_benchmark")


def plot_task15() -> None:
    holen = load_json("task15_holen.json")
    birm = load_json("task15_birmingham.json")
    metrics = holen["metrics"]["S_combo"]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4), gridspec_kw={"width_ratios": [1, 1.25]})
    labels = ["Spearman\n$\\rho$", "Pearson\n$r$", "$R^2$"]
    vals = [metrics["spearman"], metrics["pearson"], metrics["r2"]]
    bars = axes[0].bar(labels, vals, color="#3C8DAD", edgecolor="#263238", linewidth=0.7)
    annotate_bars(axes[0], bars)
    axes[0].axhline(0, color="#333", linewidth=0.7)
    axes[0].set_ylim(-0.05, 0.38)
    axes[0].set_ylabel("Metric value")
    axes[0].set_title("Holen 50 (n = 50)")
    axes[0].grid(axis="y", color="#D9DEE5", linewidth=0.6)
    axes[0].set_axisbelow(True)

    axes[1].barh(["Rows", "Binary gold label", "Prediction score", "ROC/AUC"],
                 [362, 0, 0, 0], color=["#3C8DAD", "#D0D4D8", "#D0D4D8", "#D0D4D8"],
                 edgecolor="#30343B", linewidth=0.6)
    axes[1].set_xlim(0, 400)
    axes[1].set_xlabel("Available records / estimable items")
    axes[1].set_title("Birmingham 362: AUC not estimable")
    axes[1].grid(axis="x", color="#D9DEE5", linewidth=0.6)
    axes[1].set_axisbelow(True)
    axes[1].text(7, 1, "Not available", va="center", fontsize=8)
    axes[1].text(7, 2, "Not available", va="center", fontsize=8)
    axes[1].text(7, 3, "Not estimable", va="center", fontsize=8)
    fig.suptitle("Task 15. Mismatch back-testing", y=1.03, fontweight="bold")
    fig.text(0.5, -0.02, "Guide-side/mRNA-side stratification was not available in the supplied records.",
             ha="center", fontsize=8, color="#555B64")
    save(fig, "task15_mismatch")


def plot_task16() -> None:
    data = load_json("task16_ablation.json")
    grid = data["grid"]
    alpha = np.array([row["alpha_thermo"] for row in grid])
    recommended = data["recommended"]["alpha_thermo"]
    specs = [("spearman", "Spearman $\\rho$"), ("pearson", "Pearson $r$"), ("r2", "$R^2$")]
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.2), sharex=True)
    for ax, (key, label) in zip(axes, specs):
        values = np.array([row[key] for row in grid])
        ax.plot(alpha, values, marker="o", color="#C44E52", linewidth=1.8, markersize=5)
        best = int(np.argmax(values))
        ax.scatter([alpha[best]], [values[best]], s=70, facecolor="#F2C14E", edgecolor="#333", zorder=3)
        ax.axvline(recommended, color="#555B64", linestyle="--", linewidth=0.9)
        ax.set_title(label)
        ax.set_xlabel("Thermodynamic weight $\\alpha$")
        ax.set_xticks(alpha)
        ax.grid(color="#D9DEE5", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.text(recommended + 0.02, ax.get_ylim()[0], "recommended", rotation=90,
                va="bottom", fontsize=7, color="#555B64")
    axes[0].set_ylabel("Metric value")
    fig.suptitle("Task 16. Weight ablation on the available mismatch set (n = 52)",
                 y=1.03, fontweight="bold")
    fig.text(0.5, -0.02, "$\\beta = 1-\\alpha$; highlighted point marks the maximum of each metric.",
             ha="center", fontsize=8, color="#555B64")
    save(fig, "task16_ablation")


def plot_task17() -> None:
    rows = list(csv.DictReader((ANALYSIS / "task17_positive_controls.csv").open(
        encoding="utf-8-sig", newline="")))
    labels = [row["id"] for row in rows]
    values = [float(row["rank_percentile"]) * 100 if row.get("rank_percentile") else np.nan for row in rows]
    colors = ["#4C956C" if np.isfinite(value) else "#D0D4D8" for value in values]
    plot_values = [value if np.isfinite(value) else 0 for value in values]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    bars = ax.bar(labels, plot_values, color=colors, edgecolor="#30343B", linewidth=0.6)
    for bar, value in zip(bars, values):
        if np.isfinite(value):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.5, f"{value:.1f}%",
                    ha="center", fontsize=8)
        else:
            ax.text(bar.get_x() + bar.get_width() / 2, 3, "NA", ha="center", fontsize=8)
    ax.invert_yaxis()
    ax.set_ylim(100, -8)
    ax.set_ylabel("Final-rank percentile (lower is better)")
    ax.set_xlabel("Published SFRP1 positive control")
    ax.set_title("Task 17. Positive-control placement in the ranked pipeline", fontweight="bold")
    ax.grid(axis="y", color="#D9DEE5", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.text(0.99, 0.03, "5/7 matched; gray = no exact CDS match", transform=ax.transAxes,
            ha="right", fontsize=8, color="#555B64")
    save(fig, "task17_positive_controls")


def plot_task17_weight_stability() -> None:
    rank_path = ROOT.parent / "数据集" / "阶段五_计算验证" / "task16_ablation" / "rank_final.csv"
    if not rank_path.exists():
        return
    ranked = list(csv.DictReader(rank_path.open(encoding="utf-8-sig", newline="")))
    controls = list(csv.DictReader((ANALYSIS / "task17_positive_controls.csv").open(
        encoding="utf-8-sig", newline="")))
    usable = []
    for row in ranked:
        try:
            usable.append({
                "target": row.get("target_mRNA_19", "").upper(),
                "thermo": float(row["score_thermo"]),
                "oligo": float(row["score_oligo"]),
                "penalty": float(row.get("penalty_total", 0) or 0),
            })
        except (KeyError, TypeError, ValueError):
            continue
    alphas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    output = []
    for control in controls:
        target = control.get("target19", "").upper()
        if not target or not control.get("matched") == "1":
            continue
        matched = [row for row in usable if row["target"] == target]
        if not matched:
            continue
        ranks = []
        for alpha in alphas:
            scores = sorted((alpha * row["thermo"] + (1 - alpha) * row["oligo"] - row["penalty"]
                             for row in usable), reverse=True)
            score = alpha * matched[0]["thermo"] + (1 - alpha) * matched[0]["oligo"] - matched[0]["penalty"]
            ranks.append(scores.index(score) + 1)
        output.append({"id": control["id"], **{f"alpha_{alpha:.1f}": rank for alpha, rank in zip(alphas, ranks)}})
    write_csv(ANALYSIS / "task17_weight_stability.csv", output)
    if not output:
        return
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    for row in output:
        values = [int(row[f"alpha_{alpha:.1f}"]) for alpha in alphas]
        ax.plot(alphas, values, marker="o", linewidth=1.4, markersize=4, label=row["id"])
    ax.invert_yaxis()
    ax.set_xticks(alphas)
    ax.set_xlabel("Thermodynamic weight $\\alpha$")
    ax.set_ylabel("Rank (lower is better)")
    ax.set_title("Task 17. Positive-control rank stability across weights", fontweight="bold")
    ax.grid(color="#D9DEE5", linewidth=0.6)
    ax.legend(frameon=False, ncol=3, fontsize=8)
    fig.text(0.5, -0.02, "This is ranking stability; no experimental efficacy labels were available for these controls.",
             ha="center", fontsize=8, color="#555B64")
    save(fig, "task17_weight_stability")


if __name__ == "__main__":
    plot_task14()
    plot_task15()
    plot_task16()
    plot_task17()
    plot_task17_weight_stability()
    print(f"Wrote figures to {FIGURES}")
