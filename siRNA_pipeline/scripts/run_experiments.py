"""Reproducible analysis entry point for tasks 14-17 and final freezing.

The script deliberately does not download or invent benchmark/control data.  It
consumes CSVs supplied by the user and writes auditable JSON/CSV artifacts.
Examples (from the repository root)::

    python scripts/run_experiments.py check
    python scripts/run_experiments.py benchmark --input-dir "OligoFormer部分"
    python scripts/run_experiments.py mismatch --input Holen50.csv --kind holen
    python scripts/run_experiments.py mismatch --input Birmingham362.csv --kind birmingham
    python scripts/run_experiments.py ablation --input mismatch_features.csv
    python scripts/run_experiments.py positive --input data/positive_controls.csv --rank outputs/results/rank_final.csv

Input contracts are described in docs/design/experiment_inputs.md.

口径修正（2026-09-11，L1）
-------------------------
1. `pure_thermo` 改用**管道 ranker 的默认权重/方向**（`stages/rank/ranker.py::
   DEFAULT_RANK_CFG`）计算，而不是脚本自造的"三特征等权 min-max"，并把所用口径
   记录为 `thermo_method`；
2. 所有指标附加 **bootstrap 95% CI** 与 **置换检验 p 值**，并给出分数/标签的唯一值
   个数（n=52 这类小样本必须看区间，不能只看点估计）；
3. α（热力权重）推荐值改为 **K 折交叉验证**选择（`recommended_cv`），保留全体样本
   argmax 作为 `recommended` 仅作对照，避免过拟合选参；
4. `mismatch` 自动从 siRNA/mRNA 序列**推导错配位点、位置带、guide/mRNA 侧**，
   据此分层输出；Birmingham 行可由"沉默效率"里的 log2 双重复**推导二分类标签**；
5. Birmingham 的 AUC 只要存在可用预测列就会计算（缺预测列时明确报缺，不伪造）。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "outputs" / "analysis"

_COMP = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A", "X": "X"}


def _read(path: Path) -> list[dict[str, str]]:
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            with path.open(encoding=encoding, newline="") as fh:
                return list(csv.DictReader(fh))
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("utf-8", b"", 0, 1, f"无法解码CSV: {path}")


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def _numbers(rows, column):
    values = []
    for row in rows:
        try:
            values.append(float(row[column]))
        except (KeyError, TypeError, ValueError):
            pass
    return values


def _rank(values):
    order = sorted(range(len(values)), key=lambda i: (values[i], i))
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end - 1) / 2.0 + 1
        for i in order[start:end]:
            result[i] = rank
        start = end
    return result


def _pearson(x, y):
    if len(x) != len(y) or len(x) < 2:
        return None
    mx, my = sum(x) / len(x), sum(y) / len(y)
    xx = sum((v - mx) ** 2 for v in x)
    yy = sum((v - my) ** 2 for v in y)
    if xx == 0 or yy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(xx * yy)


# ---------------------------------------------------------------------------
# 显著性：bootstrap CI 与置换检验（纯 Python，规模自适应）
# ---------------------------------------------------------------------------
def _spearman(x, y):
    return _pearson(_rank(x), _rank(y))


def _bootstrap_ci(x, y, n_boot: int, seed: int = 0):
    n = len(x)
    if n < 8 or n_boot <= 0:
        return None, None
    rnd = random.Random(seed)
    out = []
    for _ in range(n_boot):
        idx = [rnd.randrange(n) for _ in range(n)]
        s = _spearman([x[i] for i in idx], [y[i] for i in idx])
        if s is not None:
            out.append(s)
    if not out:
        return None, None
    out.sort()
    return out[int(0.025 * len(out))], out[min(len(out) - 1, int(0.975 * len(out)))]


def _perm_pvalue(x, y, n_perm: int, seed: int = 0):
    obs = _spearman(x, y)
    if obs is None or len(x) < 8 or n_perm <= 0:
        return None
    rnd = random.Random(seed)
    hits = 0
    yy = list(y)
    for _ in range(n_perm):
        rnd.shuffle(yy)
        s = _spearman(x, yy)
        if s is not None and abs(s) >= abs(obs):
            hits += 1
    return (hits + 1) / (n_perm + 1)


def _stats(x, y, detail: bool = True):
    pairs = [(a, b) for a, b in zip(x, y)
             if isinstance(a, (int, float)) and isinstance(b, (int, float))
             and math.isfinite(a) and math.isfinite(b)]
    n_raw = len(list(zip(x, y)))
    if len(pairs) < 2:
        return {"n": len(pairs), "n_raw": n_raw, "spearman": None, "pearson": None, "r2": None,
                "note": "有效数值对不足（NaN/inf 已剔除）"}
    xs = [a for a, _ in pairs]
    ys = [b for _, b in pairs]
    p = _pearson(xs, ys)
    s = _pearson(_rank(xs), _rank(ys))
    out = {"n": len(xs), "n_raw": n_raw, "spearman": s, "pearson": p,
           "r2": None if p is None else p * p}
    if detail and s is not None and len(xs) >= 8:
        n = len(xs)
        n_boot = 200 if n > 500 else 1000
        n_perm = 500 if n > 500 else 2000
        lo, hi = _bootstrap_ci(xs, ys, n_boot)
        out.update({
            "ci95_low": lo, "ci95_high": hi,
            "p_perm": _perm_pvalue(xs, ys, n_perm),
            "n_boot": n_boot, "n_perm": n_perm,
            "n_unique_scores": len(set(xs)), "n_unique_labels": len(set(ys)),
        })
    return out


def _minmax(values, reverse=False):
    lo, hi = min(values), max(values)
    if lo == hi:
        return [0.5] * len(values)
    scaled = [(v - lo) / (hi - lo) for v in values]
    return [1 - v for v in scaled] if reverse else scaled


def _metric(rows, prediction, label="label", detail: bool = True):
    pairs = []
    for row in rows:
        try:
            yv = float(row[label])
            pv = float(row[prediction])
        except (KeyError, TypeError, ValueError):
            continue
        if not (math.isfinite(yv) and math.isfinite(pv)):
            continue
        pairs.append((yv, pv))
    return _stats([p[1] for p in pairs], [p[0] for p in pairs], detail) if pairs else {"n": 0}


def _first_value(row, *names):
    for name in names:
        value = row.get(name, "")
        if str(value).strip() != "":
            return value
    return ""


# ---------------------------------------------------------------------------
# 错配注释：从序列推导位点/位置带（无外部标签时使用）
# ---------------------------------------------------------------------------
def _rc(seq: str) -> str:
    return "".join(_COMP.get(b, "N") for b in reversed(seq.upper().replace("T", "U")))


def _band_of(pos: int) -> str:
    if pos <= 1:
        return "g1"
    if pos <= 8:
        return "seed(2-8)"
    if pos <= 12:
        return "central(9-12)"
    return "3'(13-19)"


def _diff_positions(guide: str, mrna: str, prefer_offsets=(19, 20, 10, 0)):
    """返回 (offset, 差异位点列表) —— 在 mRNA 上搜索最佳 19nt 靶窗。

    guide[p]（0-based）与靶窗 window[18-p] 反平行配对；窗口偏移在
    prefer_offsets 优先、随后全序列搜索，取**差异最少**的窗口（容错 39/57/67nt
    等不同框长与 X 填充）。缺列/长度不足返回 None。
    """
    guide = (guide or "").upper().replace("T", "U")
    mrna = (mrna or "").upper().replace("T", "U")
    if len(guide) != 19 or len(mrna) < 19:
        return None
    offsets = list(dict.fromkeys(list(prefer_offsets) + list(range(0, len(mrna) - 18))))
    best = None
    for off in offsets:
        window = mrna[off:off + 19]
        if len(window) != 19:
            continue
        expect = _rc(window)
        pos = [p + 1 for p in range(19) if guide[p] != expect[p]]
        if best is None or len(pos) < len(best[2]):
            best = (off, len(pos), pos)
    return (best[0], best[2]) if best else None


def _annotate_mismatches(rows) -> int:
    """就地补 diff_positions / diff_count / diff_band / diff_offset / mismatch_side，返回成功行数。"""
    guide_cols = ("guide_checked", "siRNA", "guide_seq", "guide_wt_19")
    mrna_cols = ("mRNA_57nt", "mRNA", "mRNA_window")
    n_ok = 0
    for row in rows:
        g = _first_value(row, *guide_cols)
        m = _first_value(row, *mrna_cols)
        res = _diff_positions(g, m)
        if res is None:
            continue
        off, pos = res
        row["diff_positions"] = ",".join(str(p) for p in pos)
        row["diff_count"] = str(len(pos))
        row["diff_band"] = "|".join(sorted({_band_of(p) for p in pos})) if pos else "none"
        row["diff_offset"] = str(off)
        row.setdefault("mismatch_side", "guide")
        n_ok += 1
    return n_ok


def _fill_mfe(rows, guide_candidates=("guide_checked", "siRNA", "guide_seq", "guide_wt_19")):
    """用 RNAfold CLI 为每行补 MFE_guide（kcal/mol）；不可用时返回 available=False。

    只依赖 PATH 或环境变量 VIENNARNA_BIN，不硬编码机器路径。
    """
    exe = None
    env_bin = os.environ.get("VIENNARNA_BIN")
    if env_bin:
        for name in ("RNAfold.exe", "RNAfold"):
            cand = Path(env_bin) / name
            if cand.exists():
                exe = str(cand)
                break
    if exe is None:
        exe = shutil.which("RNAfold") or shutil.which("RNAfold.exe")
    if exe is None:
        return {"available": False, "reason": "RNAfold 不在 PATH，且未设置 VIENNARNA_BIN"}
    cache, filled = {}, 0
    for row in rows:
        guide = _first_value(row, *guide_candidates).upper().replace("T", "U")
        if len(guide) != 19:
            continue
        if guide not in cache:
            try:
                proc = subprocess.run([exe, "--noPS", "-T", "37"], input=guide + "\n",
                                      capture_output=True, text=True, timeout=20)
                hit = re.search(r"\(\s*([-+]?\d+\.\d+)\s*\)", proc.stdout or "")
                cache[guide] = float(hit.group(1)) if hit else None
            except Exception:  # noqa: BLE001
                cache[guide] = None
        if cache[guide] is not None:
            row["MFE_guide"] = "%.3f" % cache[guide]
            filled += 1
    return {"available": True, "exe": exe, "filled": filled, "unique_guides": len(cache)}


def _parse_log2(text):
    """从 '−0.33;-0.12 (微阵列log2比值,双重复)' 解析双重复均值；无数字返回 None。"""
    vals = [float(v) for v in re.findall(r"[-+]?\d+\.?\d*", str(text or ""))]
    vals = [v for v in vals if abs(v) < 10]          # 排除类似 "2006" 的年份噪声
    return sum(vals) / len(vals) if vals else None


# ---------------------------------------------------------------------------
# 管道特征计算（features 子命令）：MFE + 双链 ΔG + seed 结合能 + 末端 ΔΔG
# ---------------------------------------------------------------------------
def _nn_tables():
    sys.path.insert(0, str(ROOT / "src"))
    from sirna_pipeline.common import duplex as _duplex        # type: ignore
    from sirna_pipeline.common.nn_tables import STACK_DG       # type: ignore
    return STACK_DG, _duplex


def _end_diff(guide: str) -> float:
    """末端 ΔΔG（与 structure 模块同口径，不依赖 ViennaRNA）。"""
    stack, _ = _nn_tables()
    dg5 = stack.get(guide[0:2], 0.0)
    dg3 = stack.get(guide[17:19], 0.0)
    if guide[0] in "AU":
        dg5 += 0.45
    if guide[18] in "AU":
        dg3 += 0.45
    return dg5 - dg3


def _seed_canonical(guide: str, seed_kind: str = "seed7") -> float:
    """seed 区 canonical NN 堆叠能（负值=越强）。"""
    stack, duplex_mod = _nn_tables()
    pos = list(duplex_mod.seed_positions(seed_kind))
    total = 0.0
    for a, b in zip(pos, pos[1:]):
        if b == a + 1:
            total += stack.get(guide[a - 1] + guide[b - 1], 0.0)
    return total


def _load_thermo_calculator():
    """加载 legacy 热力学计算器（只读引用，算法零改动）。"""
    path = ROOT / "src" / "sirna_pipeline" / "stages" / "thermo" / "legacy" / "thermo_calculator.py"
    if not path.exists():
        return None, f"missing:{path}"
    try:
        sys.path.insert(0, str(ROOT / "src"))
        spec = importlib.util.spec_from_file_location("_thermo_legacy_exp", path)
        mod = importlib.util.module_from_spec(spec)          # type: ignore[arg-type]
        spec.loader.exec_module(mod)                        # type: ignore[union-attr]
        return mod.siRNAThermoCalculator(), None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}:{exc}"


def _features(args):
    """为任意 (siRNA, 57nt mRNA[, label]) 表计算管道特征，供 task14/15 用真实口径评估。"""
    rows = _read(Path(args.input))
    out_dir = Path(args.out or DEFAULT_OUT)
    out_csv = Path(args.out_csv) if args.out_csv else out_dir / ("features_" + Path(args.input).stem + ".csv")
    guide_cols = ("siRNA", "guide_checked", "guide_seq")
    mrna_cols = ("mRNA_57nt", "mRNA")
    mfe_info = {"available": False, "reason": "disabled"} if args.no_mfe else _fill_mfe(rows)
    calc, calc_err = _load_thermo_calculator()
    seed_kind = args.seed_kind
    notes, n_thermo = {}, 0
    for row in rows:
        guide = _first_value(row, *guide_cols).upper().replace("T", "U")
        mrna = _first_value(row, *mrna_cols).upper().replace("T", "U")
        if len(guide) != 19:
            notes["guide_len_ne_19"] = notes.get("guide_len_ne_19", 0) + 1
            continue
        row["GC_content"] = "%.1f" % (100.0 * (guide.count("G") + guide.count("C")) / 19.0)
        row["delta_deltaG_ends"] = "%.4f" % _end_diff(guide)
        # 取 57nt（不足则按最佳窗口 + X 补足；X 在热力学计算中替换为 A 并标注）
        if len(mrna) == 57:
            m57 = mrna
        elif len(mrna) >= 19:
            res = _diff_positions(guide, mrna)
            off = res[0] if res else 19
            m57 = "X" * 19 + mrna[off:off + 19] + "X" * 19
            notes["mrna_rebuilt_with_X"] = notes.get("mrna_rebuilt_with_X", 0) + 1
        else:
            notes["mrna_too_short"] = notes.get("mrna_too_short", 0) + 1
            continue
        if any(c not in "ACGU" for c in m57):
            notes["X_replaced_by_A"] = notes.get("X_replaced_by_A", 0) + 1
        if calc is None:
            continue
        try:
            res = calc.calculate(guide, m57.replace("X", "A"))
            row["dG_total"] = "%.4f" % res["dG_total"]
            row["dG_seed_only_canonical"] = "%.4f" % _seed_canonical(guide, seed_kind)
            by_pos = res.get("dG_mismatch_by_position") or []
            _stack, duplex_mod = _nn_tables()
            row["dG_seed"] = "%.4f" % (_seed_canonical(guide, seed_kind) + sum(
                by_pos[p - 1] for p in duplex_mod.seed_positions(seed_kind) if p - 1 < len(by_pos)))
            row["n_mismatch"] = str(res.get("n_mismatch", ""))
            n_thermo += 1
        except Exception:  # noqa: BLE001
            notes["thermo_calc_error"] = notes.get("thermo_calc_error", 0) + 1
    n_annot = _annotate_mismatches(rows)
    _write_csv(out_csv, rows)
    result = {"task": "features", "input": str(args.input), "output": str(out_csv),
              "rows": len(rows), "thermo_rows": n_thermo,
              "mfe": mfe_info, "annotation": {"rows_with_derived_positions": n_annot},
              "notes": notes, "thermo_error": calc_err,
              "columns_added": ["GC_content", "delta_deltaG_ends", "MFE_guide", "dG_total",
                                "dG_seed", "dG_seed_only_canonical", "n_mismatch",
                                "diff_positions", "diff_band"]}
    if args.out_csv:
        result.pop("output", None)
    else:
        _write_json(out_dir / ("features_" + Path(args.input).stem + ".json"), result)
    return result


# ---------------------------------------------------------------------------
# 与管道一致的"热力分"（避免自造等权口径）
# ---------------------------------------------------------------------------
_PIPE_FEATURE_SOURCES = {
    "feat_mfe": ("feat_mfe", "mfe", "MFE_guide"),
    "feat_ddg_ends": ("feat_ddg_ends", "delta_deltaG_ends", "ddg_ends"),
    "feat_dG_duplex": ("feat_dG_duplex", "dG_total", "dG_duplex"),
    "feat_dG_seed": ("feat_dG_seed", "dG_seed", "seed_dG"),
}


def _ranker_alpha_beta() -> tuple[float, float]:
    """从 ranker 默认配置读取主分/辅助分权重（与生产排序保持同一事实来源）。"""
    try:
        sys.path.insert(0, str(ROOT / "src"))
        from sirna_pipeline.stages.rank.ranker import DEFAULT_RANK_CFG  # type: ignore
        return (float(DEFAULT_RANK_CFG.get("alpha", 0.4)),
                float(DEFAULT_RANK_CFG.get("beta", 0.6)))
    except Exception:  # noqa: BLE001
        return 0.4, 0.6


def _pipeline_thermo_score(rows):
    """用管道 ranker 默认配置（权重/方向）计算热力分；返回 (scores, method) 或 (None, 原因)。"""
    try:
        sys.path.insert(0, str(ROOT / "src"))
        from sirna_pipeline.stages.rank.ranker import DEFAULT_RANK_CFG  # type: ignore
    except Exception as exc:  # noqa: BLE001
        return None, f"ranker_unavailable:{type(exc).__name__}"
    feats = DEFAULT_RANK_CFG["features"]
    values, used_cols = {}, {}
    for feat in feats:
        name = feat["name"]
        col = None
        for cand in _PIPE_FEATURE_SOURCES.get(name, (name,)):
            vals = []
            for row in rows:
                try:
                    v = float(row.get(cand, ""))
                except (TypeError, ValueError):
                    continue
                if math.isfinite(v):
                    vals.append(v)
            if len(vals) == len(rows) and vals:      # 该列必须**全行有效**才可用
                col = cand
                break
        if col is None:
            return None, f"missing_feature:{name}"
        values[name], used_cols[name] = vals, col
    norm = {}
    for feat in feats:
        vals = values[feat["name"]]
        if int(feat.get("direction", 1)) < 0:
            vals = [-v for v in vals]
        norm[feat["name"]] = _minmax(vals)
    wsum = sum(float(f["weight"]) for f in feats) or 1.0
    scores = [sum(float(f["weight"]) * norm[f["name"]][i] for f in feats) / wsum
              for i in range(len(rows))]
    return scores, "ranker_default_config:" + ",".join(
        f"{k}<-{v}" for k, v in sorted(used_cols.items()))


# ---------------------------------------------------------------------------
# K 折交叉验证选择 α（避免在全体样本上 argmax）
# ---------------------------------------------------------------------------
def _weighted_score(rows, spec):
    """spec=[(列名, 权重, 方向)]；方向>0 越大越好、<0 越负越好。返回加权归一化分或 None。"""
    cols, weights, norms = [], [], {}
    for col, weight, direction in spec:
        vals = []
        for row in rows:
            try:
                v = float(row.get(col, ""))
            except (TypeError, ValueError):
                return None
            if not math.isfinite(v):
                return None
            vals.append(v if direction > 0 else -v)
        cols.append(col)
        weights.append(float(weight))
        norms[col] = _minmax(vals)
    wsum = sum(weights) or 1.0
    return [sum(w * norms[c][i] for c, w in zip(cols, weights)) / wsum
            for i in range(len(rows))]


def _cv_mix(pairs, alphas, k: int = 5, repeats: int = 50, seed: int = 0):
    """pairs=[(score_a, score_b, label)]；alpha 为 score_a 的权重，返回 (均值表, 训练折众数)。"""
    n = len(pairs)
    if n < 10:
        return {a: None for a in alphas}, None
    k = max(2, min(k, n))
    rnd = random.Random(seed)
    acc = {a: [] for a in alphas}
    picked = []
    for _ in range(repeats):
        idx = list(range(n))
        rnd.shuffle(idx)
        folds = [idx[i::k] for i in range(k)]
        for a in alphas:
            rhos = []
            for fold in folds:
                test = [pairs[i] for i in fold]
                if len(test) < 3:
                    continue
                s = _spearman([a * x + (1 - a) * y for x, y, _ in test],
                              [lab for _, _, lab in test])
                if s is not None:
                    rhos.append(s)
            if rhos:
                acc[a].append(sum(rhos) / len(rhos))
        for fold in folds:
            test = [pairs[i] for i in fold]
            train = [pairs[i] for i in idx if i not in set(fold)]
            if len(train) < 6 or len(test) < 3:
                continue
            best = max(alphas, key=lambda a: (_spearman(
                [a * x + (1 - a) * y for x, y, _ in train],
                [lab for _, _, lab in train]) or -2))
            picked.append(best)
    means = {a: (sum(v) / len(v) if v else None) for a, v in acc.items()}
    return means, (max(set(picked), key=picked.count) if picked else None)


def _cv_alpha(rows, alphas, k: int = 5, repeats: int = 50, seed: int = 0):
    pairs = []
    for row in rows:
        try:
            pairs.append((float(row["thermo_score"]), float(row["dl_score"]), float(row["label"])))
        except (KeyError, TypeError, ValueError):
            continue
    return _cv_mix(pairs, alphas, k, repeats, seed)


# ---------------------------------------------------------------------------
# 子命令
# ---------------------------------------------------------------------------
def _benchmark(args):
    source = Path(args.input_dir)
    dataset_source = Path(args.dataset_dir) if args.dataset_dir else source / "data"
    out = Path(args.out or DEFAULT_OUT)
    datasets = ["Hu", "Taka", "Mix", "Simone"]
    result = {"task": "14", "datasets": {}, "missing": []}
    for name in datasets:
        dataset = dataset_source / f"{name}.csv"
        prediction = source / f"{name}_predictions.csv"
        if not dataset.exists() or not prediction.exists():
            result["missing"].append({"dataset": name, "dataset_file": str(dataset),
                                      "prediction_file": str(prediction)})
            continue
        rows, preds = _read(dataset), _read(prediction)
        prediction_by_key = {}
        has_prediction_keys = all(any(name in pred for pred in preds)
                                  for name in ("siRNA", "mRNA")) if preds else False
        if has_prediction_keys:
            for pred in preds:
                prediction_by_key[(pred.get("siRNA", ""), pred.get("mRNA", ""))] = pred
        joined = []
        unmatched = 0
        for row in rows:
            pred = (prediction_by_key.get((row.get("siRNA", ""), row.get("mRNA", "")))
                    if has_prediction_keys else None)
            if pred is None and not has_prediction_keys:
                index = len(joined)
                pred = preds[index] if index < len(preds) else {}
            if not pred:
                unmatched += 1
                continue
            merged = dict(row)
            merged["dl_score"] = _first_value(pred, "pred", "pred_efficiency")
            merged["label"] = _first_value(row, "label", "true_efficiency",
                                           "silencing_efficiency_norm", "沉默效率")
            joined.append(merged)
        n = len(joined)
        metrics = {"pure_dl": _metric(joined, "dl_score")}
        dataset_result = {"rows": n, "metrics": metrics}
        if unmatched or len(rows) != len(preds) or not has_prediction_keys:
            dataset_result["prediction_count_warning"] = {
                "dataset_rows": len(rows), "prediction_rows": len(preds),
                "used_rows": n, "unmatched_dataset_rows": unmatched,
                "join_mode": "key" if has_prediction_keys else "row_order_fallback",
            }
        # ---- 口径修正：热力分 = 管道 ranker 默认权重/方向（不再自造等权 min-max）----
        thermo_scores, thermo_method = _pipeline_thermo_score(joined)
        dataset_result["thermo_method"] = thermo_method
        labels = _numbers(joined, "label")
        dl = _numbers(joined, "dl_score")
        if thermo_scores is not None and len(labels) == len(thermo_scores) == len(dl):
            dl_n = _minmax(dl)
            combo_alpha, combo_beta = _ranker_alpha_beta()
            scored = [{"label": str(y), "thermo_score": str(t),
                       "combo_score": str(combo_alpha * t + combo_beta * d)}
                      for y, t, d in zip(labels, thermo_scores, dl_n)]
            dataset_result["combo_alpha_beta"] = [combo_alpha, combo_beta]
            metrics["pure_thermo"] = _metric(scored, "thermo_score")
            metrics["thermo_plus_aux"] = _metric(scored, "combo_score")
        else:
            dataset_result.setdefault("missing_methods", []).extend(["pure_thermo", "thermo_plus_aux"])
        result["datasets"][name] = dataset_result
    _write_json(out / "task14_benchmark.json", result)
    return result


def _auc(labels, scores):
    positives = [s for y, s in zip(labels, scores) if y == 1]
    negatives = [s for y, s in zip(labels, scores) if y == 0]
    if not positives or not negatives:
        return None
    wins = sum(1 if p > n else 0.5 if p == n else 0 for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


def _mismatch(args):
    rows = _read(Path(args.input))
    if args.prediction_file:
        predictions = _read(Path(args.prediction_file))
        if args.reference:
            rows = _read(Path(args.reference))
        keyed = {}
        keyed_by_siRNA = {}
        for prediction in predictions:
            key = (prediction.get("siRNA", ""), prediction.get("mRNA", prediction.get("mRNA_57nt", "")))
            keyed.setdefault(key, []).append(prediction)
            keyed_by_siRNA.setdefault(prediction.get("siRNA", ""), []).append(prediction)
        for row in rows:
            key = (row.get("siRNA", ""), row.get("mRNA", row.get("mRNA_57nt", "")))
            prediction = keyed.get(key, [{}]).pop(0) if keyed.get(key) else (
                keyed_by_siRNA.get(row.get("siRNA", ""), [{}]).pop(0)
                if keyed_by_siRNA.get(row.get("siRNA", "")) else {})
            for name in (args.prediction, "pred_efficiency", "true_efficiency", "S_combo"):
                if name in prediction:
                    row[name] = prediction[name]
    out = Path(args.out or DEFAULT_OUT)
    # ---- 先补 MFE（RNAfold CLI），再算口径分；否则"缺 MFE → 热力分算不出" ----
    mfe_info = _fill_mfe(rows) if getattr(args, "with_mfe", False) else None
    n_annot = _annotate_mismatches(rows)
    label_col = args.label
    if args.kind != "birmingham" and rows and not any(label_col in row for row in rows):
        label_col = next((c for c in ("label", "true_efficiency", "silencing_efficiency_norm", "沉默效率")
                          if any(c in row for row in rows)), label_col)
    prediction_cols = [c for c in (args.prediction, "final_score", "pred_efficiency")
                       if c and any(c in row for row in rows)]
    if not prediction_cols and rows:
        prediction_cols = [c for c in ("pred_efficiency", "S_combo", "silencing_efficiency_norm")
                           if any(c in row for row in rows)]
    # ---- 新增：可直接用"管道 ranker 口径"的热力分作为预测列（推荐）----
    thermo_method = None
    if str(args.prediction) in ("pipeline_thermo", "pipeline_thermo_combo"):
        scores, thermo_method = _pipeline_thermo_score(rows)
        if scores is not None:
            dl_vals = _numbers(rows, "dl_score") or _numbers(rows, "oligo_pred") or _numbers(rows, "pred_efficiency")
            dl_norm = _minmax(dl_vals) if len(dl_vals) == len(rows) else [0.5] * len(rows)
            target = "pipeline_thermo" if str(args.prediction) == "pipeline_thermo" else "pipeline_thermo_combo"
            for row, s, d in zip(rows, scores, dl_norm):
                row[target] = str(s if target == "pipeline_thermo" else 0.8 * s + 0.2 * d)
            prediction_cols = [target]
    expected = {"holen": 50, "birmingham": 362}[args.kind]
    result = {"task": "15", "kind": args.kind, "rows": len(rows), "expected_rows": expected,
              "metrics": {},
              "strata": {}, "missing": []}
    if thermo_method:
        result["thermo_method"] = thermo_method
    if len(rows) != expected:
        result["row_count_warning"] = f"{args.kind} 期望 {expected} 行，实际 {len(rows)} 行；请勿作为正式基准结论。"

    # ---- 错配注释与 MFE 已在前面完成，这里只登记 ----
    result["annotation"] = {"rows_with_derived_positions": n_annot}
    if mfe_info is not None:
        result["mfe_fill"] = mfe_info

    # ---- 新增：Birmingham 由 log2 双重复推导二分类标签 ----
    if args.kind == "birmingham":
        thr = float(os.environ.get("SIRNA_EXP_LOG2_THRESHOLD", "-0.3"))
        derived = 0
        for row in rows:
            if str(row.get(label_col, "")).strip() not in ("", "nan"):
                continue
            value = _parse_log2(row.get("沉默效率", ""))
            if value is None:
                continue
            row[label_col] = "1" if value <= thr else "0"
            row["log2_mean"] = str(value)
            derived += 1
        if derived:
            result["derived_labels"] = {
                "count": derived, "rule": f"mean(log2 ratio) <= {thr}",
                "note": "阈值可用环境变量 SIRNA_EXP_LOG2_THRESHOLD 覆盖；论文口径不同请自行调整",
            }

    for col in prediction_cols[:1]:
        result["metrics"][col] = _metric(rows, col, label_col)
    # AUC：只要标签是二分类且有预测列就计算（不再要求 label_col == args.label）
    if prediction_cols:
        labels, scores = [], []
        for row in rows:
            try:
                labels.append(int(float(row[label_col])))
                scores.append(float(row[prediction_cols[0]]))
            except (IndexError, KeyError, TypeError, ValueError):
                continue
        if len(set(labels)) == 2:
            result["metrics"][prediction_cols[0] + "_roc"] = {
                "n": len(labels), "auc": _auc(labels, scores),
                "positives": sum(labels), "negatives": len(labels) - sum(labels),
            }
    if args.kind == "birmingham" and not prediction_cols:
        result["missing"].append(
            "预测列缺失：Birmingham 为 seed 型脱靶数据，需要先为这些行生成预测分"
            "（如 PITA/TargetScan 或 seed-ΔG 打分）后才能算 AUC")
    if args.kind == "birmingham" and label_col == args.label and not any(
            str(row.get(label_col, "")).strip() in {"0", "1", "0.0", "1.0"} for row in rows):
        result["missing"].append("Birmingham binary gold label")

    # ---- 分层：位置带（最有用）+ 显式 side 列，两者都给 ----
    strata_cols = []
    if any("diff_band" in row for row in rows):
        strata_cols.append("diff_band")
    if any(args.side in row for row in rows):
        strata_cols.append(args.side)
    if strata_cols and prediction_cols:
        for col in strata_cols:
            for value in sorted({row.get(col, "") for row in rows}):
                subset = [row for row in rows if row.get(col, "") == value]
                if len(subset) >= 3:
                    result["strata"][f"{col}={value}"] = {
                        "n": len(subset),
                        "metrics": _metric(subset, prediction_cols[0], label_col, detail=False),
                    }
        result["strata_column"] = strata_cols
    else:
        result["strata_note"] = f"未找到分层列 {args.side!r}/diff_band，无法分层。"
    _write_json(out / f"task15_{args.kind}.json", result)
    return result


def _ablation(args):
    rows = _read(Path(args.input))
    out = Path(args.out or DEFAULT_OUT)
    aliases = {"label": ("label", "true_efficiency", "silencing_efficiency_norm", "沉默效率"),
               "dl_score": ("dl_score", "oligo_pred", "pred_efficiency"),
               "thermo_score": ("thermo_score", "S_thermo", "score_thermo")}
    for target, candidates in aliases.items():
        if not any(target in row for row in rows):
            source = next((c for c in candidates if any(c in row for row in rows)), None)
            if source:
                for row in rows:
                    row[target] = row.get(source, "")
    required = ["label", "dl_score", "thermo_score"]
    missing = [c for c in required if not any(c in r for r in rows)]
    result = {"task": "16", "rows": len(rows), "grid": [], "missing": missing}
    mfe_info = None
    if getattr(args, "with_mfe", False):
        mfe_info = _fill_mfe(rows)
    n_annot = _annotate_mismatches(rows)
    result["annotation"] = {"rows_with_derived_positions": n_annot}
    if mfe_info is not None:
        result["mfe_fill"] = mfe_info
    # ---- 新增：优先采用"管道 ranker 口径"热力分（含 MFE 权重/方向）----
    thermo_scores, thermo_method = _pipeline_thermo_score(rows)
    if thermo_scores is not None:
        for row, s in zip(rows, thermo_scores):
            row["thermo_score_ranker"] = "%.10g" % s
            row["thermo_score"] = row["thermo_score_ranker"]      # 覆盖队友等权口径
        result["thermo_method"] = thermo_method
        result["note_thermo_score"] = ("thermo_score 已替换为管道 ranker 口径"
                                       "（含 MFE 权重/方向）；队友等权口径保留在 S_thermo 列供对照")
    if not missing:
        labels = [float(r["label"]) for r in rows]
        alphas = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)      # 保持 6 点（历史图/测试口径）
        for alpha in alphas:
            scores = [alpha * float(r["thermo_score"]) + (1 - alpha) * float(r["dl_score"]) for r in rows]
            result["grid"].append({"alpha_thermo": alpha, "beta_dl": 1 - alpha,
                                   **_stats(scores, labels)})
        result["recommended"] = max(result["grid"], key=lambda r: (r.get("spearman") or -float("inf")))
        result["note_recommended"] = ("在全体样本上取 max Spearman，n 小时易过拟合；"
                                      "请以 recommended_cv 为准")
        # ---- 新增：K 折交叉验证选 α ----
        means, cv_pick = _cv_alpha(rows, list(alphas))
        result["cv_alpha_mean_spearman"] = means
        result["recommended_cv"] = {
            "alpha_thermo": cv_pick,
            "alpha_thermo_from_mean": max(
                (a for a, m in means.items() if m is not None),
                key=lambda a: means[a], default=None),
            "seed": 0,
            "note": "K=5 折、重复 50 次；alpha_thermo=按训练折 argmax 选出的众数",
        }
        # ---- 新增：逐特征 Spearman（含 CI/p）——"哪个特征真正有预测力"的直接证据 ----
        feat_cols = [c for c in ("thermo_score", "dl_score", "MFE_guide", "dG_total", "dG_seed",
                                 "delta_deltaG_ends", "GC_content", "S_thermo", "S_combo")
                     if any(str(r.get(c, "")).strip() not in ("", "nan") for r in rows)]
        result["feature_spearman"] = {c: _metric(rows, c, "label") for c in feat_cols}
        # ---- 新增：权重预设对照（含 MFE 主导）与 MFE 权重扫描 ----
        presets = {
            "ranker_default": [("MFE_guide", 1.0, 1), ("delta_deltaG_ends", 1.0, 1),
                               ("dG_total", 1.2, -1), ("dG_seed", 0.8, -1)],
            "mfe_dominant": [("MFE_guide", 1.8, 1), ("delta_deltaG_ends", 0.4, 1),
                             ("dG_total", 0.4, -1), ("dG_seed", 0.2, -1)],
            "mfe_only": [("MFE_guide", 1.0, 1)],
            "no_seed": [("MFE_guide", 1.4, 1), ("delta_deltaG_ends", 0.6, 1),
                        ("dG_total", 0.6, -1)],
            "seed_only": [("dG_seed", 1.0, -1)],
        }
        result["weight_presets"] = {}
        for name, spec in presets.items():
            sc = _weighted_score(rows, spec)
            result["weight_presets"][name] = (
                {"spec": spec, **_stats(sc, labels)} if sc is not None
                else {"spec": spec, "note": "所需特征列缺失，无法计算"})
        mfe_only = _weighted_score(rows, [("MFE_guide", 1.0, 1)])
        rest = _weighted_score(rows, [("delta_deltaG_ends", 1.0, 1), ("dG_total", 1.2, -1),
                                      ("dG_seed", 0.8, -1)])
        if mfe_only is not None and rest is not None:
            alphas = tuple(round(0.1 * i, 1) for i in range(11))
            full = {a: _stats([a * r + (1 - a) * m for r, m in zip(rest, mfe_only)], labels)
                    for a in alphas}
            means, picked = _cv_mix(list(zip(rest, mfe_only, labels)), alphas)
            result["mfe_weight_sweep"] = {
                "meaning": "score = a*rest_thermo_norm + (1-a)*MFE_norm（a=非MFE热力权重）",
                "full_sample": {str(a): v for a, v in full.items()},
                "cv_mean_spearman": {str(a): v for a, v in means.items()},
                "recommended_a_cv": picked,
            }
            dl_ok = all(str(r.get("dl_score", "")).strip() not in ("", "nan") for r in rows)
            if dl_ok:
                dl_n = _minmax([float(r["dl_score"]) for r in rows])
                a_star = picked if picked is not None else 0.0
                thermo_mix = [a_star * r + (1 - a_star) * m for r, m in zip(rest, mfe_only)]
                for w_dl in (0.2, 0.3):
                    combo = [(1 - w_dl) * t + w_dl * d for t, d in zip(thermo_mix, dl_n)]
                    result.setdefault("mfe_dominant_plus_dl", {})[f"w_dl={w_dl}"] = {
                        "a_rest": a_star, **_stats(combo, labels)}
        # ---- 新增：按位置带分层（列缺失时用序列推导）----
        band_col = "diff_band" if any("diff_band" in r for r in rows) else "mutation_site"
        if any(band_col in r for r in rows):
            result["mutation_sites"] = {}
            for band in sorted({r.get(band_col, "") for r in rows if r.get(band_col, "")}):
                subset = [r for r in rows if r.get(band_col, "") == band]
                if len(subset) < 3:
                    continue
                result["mutation_sites"][band] = {
                    "n": len(subset),
                    "dl": _metric(subset, "dl_score", detail=False),
                    "thermo": _metric(subset, "thermo_score", detail=False),
                }
    _write_json(out / "task16_ablation.json", result)
    return result


def _positive(args):
    controls = _read(Path(args.input))
    ranked = _read(Path(args.rank))
    by_guide = {}
    for ranked_row in ranked:
        for name in ("guide_seq", "target_seq", "target_mRNA_19", "guide_checked", "siRNA"):
            guide = ranked_row.get(name, "").upper().replace("T", "U")
            if guide:
                by_guide[guide] = ranked_row
    result_rows = []
    for row in controls:
        guide = _first_value(row, "guide_seq", "siRNA", "guide19").upper().replace("T", "U")
        hit = by_guide.get(guide, {})
        result_rows.append({**row, "matched": int(bool(hit)),
                            "rules_pass": hit.get("rules_pass", ""),
                            "structure_pass": hit.get("structure_pass", ""),
                            "final_rank": hit.get("final_rank", ""),
                            "rank_percentile": (float(hit["final_rank"]) / len(ranked)) if hit.get("final_rank") else ""})
    out = Path(args.out or DEFAULT_OUT)
    _write_csv(out / "task17_positive_controls.csv", result_rows)
    summary = {"task": "17", "controls": len(controls), "matched": sum(r["matched"] for r in result_rows),
               "unmatched": [r.get("control_id", r.get("id", "")) for r in result_rows if not r["matched"]],
               "note": "该命令只做序列与已跑管道结果 join；若需重跑，请先将阳性对照 FASTA/CSV 纳入输入契约。"}
    _write_json(out / "task17_positive_controls.json", summary)
    return summary


def _check(args):
    root = Path(args.root).resolve()
    paths_file = root / "configs" / "paths.yaml"
    configured_fasta = _configured_fasta(root)
    checks = {
        "sfrp1_fasta": configured_fasta,
        "positive_controls": root / "data" / "positive_controls.csv",
        "simone": root / "OligoFormer部分" / "data" / "Simone.csv",
        "holen50": root / "data" / "Holen50.csv",
        "birmingham362": root / "data" / "Birmingham362.csv",
        "vienna_rna": Path(__import__("os").environ["VIENNARNA_BIN"]) if __import__("os").environ.get("VIENNARNA_BIN") else Path("__VIENNARNA_BIN_UNSET__"),
    }
    result = {"root": str(root), "configured_sfrp1_fasta": str(configured_fasta),
              "available": {k: p.exists() for k, p in checks.items()},
              "missing": [k for k, p in checks.items() if not p.exists()]}
    _write_json(Path(args.out or DEFAULT_OUT) / "input_check.json", result)
    return result


def _freeze(args):
    root = Path(args.root).resolve()
    cfg = root / "configs" / "paths.yaml"
    result = {"task": "freeze", "status": "blocked", "missing": []}
    if not cfg.exists():
        result["missing"].append(str(cfg))
    if not _configured_fasta(root).exists():
        result["missing"].append(str(_configured_fasta(root)))
    vienna_binding = _vienna_binding_available()
    if not __import__("os").environ.get("VIENNARNA_BIN") and not vienna_binding:
        result["missing"].append("VIENNARNA_BIN")
    if result["missing"]:
        _write_json(Path(args.out or DEFAULT_OUT) / "freeze_status.json", result)
        return result
    command = [sys.executable, str(root / "predict.py"), "--run-name", args.run_name]
    completed = subprocess.run(command, cwd=root, text=True, capture_output=True, check=False)
    result.update({"status": "ok" if completed.returncode == 0 else "failed", "returncode": completed.returncode,
                   "vienna_backend": "python_binding" if vienna_binding else "cli",
                   "results": str(root / "outputs" / "results" / "results.csv")})
    _write_json(Path(args.out or DEFAULT_OUT) / "freeze_status.json", result)
    return result


def _vienna_binding_available() -> bool:
    """Probe the active interpreter, including environments with RNA bindings."""
    probe = subprocess.run(
        [sys.executable, "-c", "import RNA; print(RNA.__file__)"],
        capture_output=True, text=True, check=False,
    )
    return probe.returncode == 0 and bool(probe.stdout.strip())


def _configured_fasta(root: Path) -> Path:
    """Resolve the simple relative cds_fasta entry used by paths.yaml."""
    paths_file = root / "configs" / "paths.yaml"
    configured = root / "data" / "SFRP1-mRNA.txt"
    if not paths_file.exists():
        return configured
    for line in paths_file.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("cds_fasta:"):
            raw = line.split(":", 1)[1].split("#", 1)[0].strip()
            return (paths_file.parent / raw).resolve()
    return configured


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check"); check.add_argument("--root", type=Path, default=ROOT); check.add_argument("--out", type=Path)
    bench = sub.add_parser("benchmark"); bench.add_argument("--input-dir", type=Path, default=ROOT / "OligoFormer部分"); bench.add_argument("--dataset-dir", type=Path); bench.add_argument("--out", type=Path)
    feat = sub.add_parser("features"); feat.add_argument("--input", type=Path, required=True, help="含 siRNA + 57nt mRNA（可含 label）的 CSV"); feat.add_argument("--out-csv", type=Path, help="输出特征表路径（默认 outputs/analysis/features_<stem>.csv）"); feat.add_argument("--seed-kind", default="seed7", choices=("seed6", "seed7")); feat.add_argument("--no-mfe", action="store_true", help="跳过 RNAfold MFE 计算"); feat.add_argument("--out", type=Path)
    mm = sub.add_parser("mismatch"); mm.add_argument("--input", type=Path, required=True); mm.add_argument("--kind", choices=("holen", "birmingham"), required=True); mm.add_argument("--prediction", default="final_score"); mm.add_argument("--prediction-file", type=Path); mm.add_argument("--reference", type=Path); mm.add_argument("--label", default="label"); mm.add_argument("--side", default="mismatch_side"); mm.add_argument("--with-mfe", action="store_true", help="用 RNAfold CLI 补 MFE_guide（需 PATH/VIENNARNA_BIN）"); mm.add_argument("--out", type=Path)
    ab = sub.add_parser("ablation"); ab.add_argument("--input", type=Path, required=True); ab.add_argument("--with-mfe", action="store_true", help="用 RNAfold CLI 补 MFE_guide"); ab.add_argument("--out", type=Path)
    pos = sub.add_parser("positive"); pos.add_argument("--input", type=Path, required=True); pos.add_argument("--rank", type=Path, required=True); pos.add_argument("--out", type=Path)
    freeze = sub.add_parser("freeze"); freeze.add_argument("--root", type=Path, default=ROOT); freeze.add_argument("--run-name", default="final_frozen"); freeze.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    result = {"check": _check, "benchmark": _benchmark, "features": _features,
              "mismatch": _mismatch, "ablation": _ablation, "positive": _positive,
              "freeze": _freeze}[args.command](args)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if not result.get("missing") and result.get("status", "ok") not in {"blocked", "failed"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
