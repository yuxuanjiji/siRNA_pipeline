# -*- coding: utf-8 -*-
"""任务 7：序列规则筛选（初级规则 ①，四项串联过滤）。

对候选库中每条 19 nt 引导链（guide / 反义链，5'->3'）依次施加下述四项规则，
命中**任一**规则即淘汰（四项串联过滤 = 逻辑 AND 全部通过才保留）：

  规则1  GC 含量：G+C 占全长的百分比必须落在 [30%, 65%]（含边界）——
         越界即淘汰；
  规则2  连续 ≥6 个相邻 G/C（6 nt 连续串中每个碱基都是 G 或 C）——
         命中即淘汰；
  规则3  连续 ≥5 个相同碱基（AAAAA / UUUUU / CCCCCC / GGGGG 型同聚物
         及其更长延伸）——命中即淘汰；
  规则4  存在回文序列（OligoFormer 口径：19 nt 内存在任意 4 nt 臂，其反向
         互补序列出现在该臂下游 ≥4 nt 处，即 ≥8 nt 的反向重复/发夹骨架）——
         命中即淘汰。

参考实现（逐条对齐，代码内嵌对拍校验）：
    OligoFormer scripts/infer.py 的 func_filter()
    (github.com/lulab/OligoFormer; Bioinformatics 2024;40(10):btae577)

与 func_filter 标签的对应关系（用于直接对拍 / 复现原模型输出）：
    func_filter = 0  通过（seq_pass = 1）
    func_filter = 1  规则1 命中：GC 含量越界
    func_filter = 2  规则3 命中：连续 ≥5 个相同碱基
    func_filter = 3  规则2 命中：连续 ≥6 个相邻 G/C
    func_filter = 4  规则4 命中：回文序列
  命中多条时只记 func_filter 源码中**先触发**的一条（hit_code / hit_rule），
  与源码逐条 continue 的判定顺序一致；seq_pass 只取决于四项是否全部通过，
  与该记分顺序无关（同一序列被淘汰的结论两种口径恒等）。

输入（默认自动定位 项目搭建\\候选序列生成 下的成品，也可 --library/--wt 显式指定）：
  * 任务 6 输出 SFRP1_task6_site_mutagenesis_library.csv：13 905 条突变体
    （每窗口 15 条），被筛序列取字段 guide_mut_19；
  * 任务 5 输出 SFRP1_task5_sliding_window_19nt_fullcomp.csv（默认并入，
    --no-wt 关闭）：927 条完全互补引导链，被筛序列取 guide_antisense_19。
  任务 5 + 任务 6 全表 = 927 × 16 条/窗口口径下的完整候选池（见
  候选序列生成/任务5_6_候选序列生成说明.md §5.2），筛后通过者进入任务 8。

输出（UTF-8-sig BOM，Excel 可直接打开；文件放 初级规则筛选 文件夹）：
  * SFRP1_task7_sequence_rules_annotated.csv  全量注释表（保留输入全部字段，
    追加 kind / guide_checked / gc_pct / fail_* / seq_pass / hit_code / hit_rule）；
  * SFRP1_task7_sequence_rules_passed.csv     通过子表（seq_pass = 1，供任务 8 接力）。

用法:
    python task7_sequence_rules_filter.py
        [--library 任务6CSV]     默认自动定位 候选序列生成/SFRP1_task6_...csv
        [--wt 任务5CSV]          默认自动定位 候选序列生成/SFRP1_task5_...csv
        [--no-wt]                只筛任务6突变库，不并入完全互补引导链
        [--out-dir 目录]         输出目录（默认 初级规则筛选 文件夹）
        [--prefix 前缀]          输出文件名前缀（默认 SFRP1_task7_sequence_rules）
        [--selftest]             只跑内置单元锚点 + func_filter 对拍后退出

运行环境：仅 Python 3.9+ 标准库，无第三方依赖；内置校验在每次运行时对全部
真实候选逐条与 OligoFormer func_filter 对拍，不通过即抛 AssertionError。
"""
from __future__ import annotations

import argparse
import csv
import itertools
import random
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# 一、常量与阈值（与 OligoFormer func_filter 一致；修改前先核对参考源码）
# ---------------------------------------------------------------------------

GC_MIN = 30.0          # GC 含量下限（%），含边界
GC_MAX = 65.0          # GC 含量上限（%），含边界
RUN_GC = 6             # 规则2：连续 ≥6 个相邻 G/C
RUN_SAME = 5           # 规则3：连续 ≥5 个相同碱基
PAL_ARM = 4            # 规则4：OligoFormer 回文判定的 4 nt 臂长（勿改）

_RNA_COMP = str.maketrans("AUGC", "UACG")          # A<->U, G<->C
_HOMO_PATTERNS = tuple(re.compile(p) for p in
                       ("AAAAA", "UUUUU", "CCCCC", "GGGGG"))   # 规则3（OligoFormer 字面）
_GC6_PATTERNS = tuple(re.compile("".join(t)) for t in
                      itertools.product(("G", "C"), repeat=RUN_GC))  # 规则2（字面枚举 64 种）

# hit_code / hit_rule 记分顺序 = func_filter 源码 continue 顺序（GC->5同->6GC->回文）。
# 注意：仅影响“先命中哪条”的标注；淘汰结论（seq_pass）与记分顺序无关。
HIT_ORDER = [
    ("fail_gc_range", 1, "gc_range_out"),
    ("fail_run5_same", 2, "run5_same"),
    ("fail_run6_gc", 3, "run6_gc"),
    ("fail_palindrome", 4, "palindrome"),
]
HIT_CODE_NAMES = {c: name for _, c, name in HIT_ORDER}

# 输出 CSV 列顺序（输入字段保留 + 任务 7 追加字段）
OUT_FIELDS = [
    "window_id", "kind", "variant_id", "site", "position", "wt_nt", "mut_nt",
    "cds_start", "cds_end", "nm_003012_start", "nm_003012_end",
    "target_mRNA_19", "guide_wt_19", "guide_mut_19", "guide_checked",
    "paired_mRNA_nt", "pair_type", "sense_strand_19", "window_gc_pct",
    "gc_pct", "fail_gc_range", "fail_run6_gc", "fail_run5_same",
    "fail_palindrome", "seq_pass", "hit_code", "hit_rule",
]

_HERE = Path(__file__).resolve().parent                 # .../项目搭建/初级规则筛选
_PROJ_DIR = _HERE.parent                                # .../项目搭建
_GEN_DIR = _PROJ_DIR / "候选序列生成"

TASK6_DEFAULT_NAME = "SFRP1_task6_site_mutagenesis_library.csv"
TASK5_DEFAULT_NAME = "SFRP1_task5_sliding_window_19nt_fullcomp.csv"
PREFIX_DEFAULT = "SFRP1_task7_sequence_rules"
ANNOTATED_NAME = PREFIX_DEFAULT + "_annotated.csv"
PASSED_NAME = PREFIX_DEFAULT + "_passed.csv"


# ---------------------------------------------------------------------------
# 二、基础序列工具
# ---------------------------------------------------------------------------

def normalize(seq: str) -> str:
    """清洗输入序列：去空白、大写、T->U；非法字符/长度非 19 时抛错。"""
    s = "".join(seq.split()).upper().replace("T", "U")
    if any(c not in "AUGC" for c in s):
        raise ValueError(f"序列含非 A/U/G/C 字符：{seq!r}")
    if len(s) != 19:
        raise ValueError(f"序列长度 {len(s)} != 19：{s}")
    return s


def rc(seq: str) -> str:
    """RNA 反向互补（5'->3' 转换）。"""
    return seq.translate(_RNA_COMP)[::-1]


def gc_pct(seq: str) -> float:
    """GC 含量（%）。"""
    return (seq.count("G") + seq.count("C")) / len(seq) * 100.0


# ---------------------------------------------------------------------------
# 三、四项规则判定（各自独立，命中返回 True）
# ---------------------------------------------------------------------------

def fail_gc_range(seq: str) -> bool:
    """规则1：GC 含量越界（<30% 或 >65%）即淘汰。"""
    gc = gc_pct(seq)
    return gc < GC_MIN or gc > GC_MAX


def fail_run6_gc(seq: str) -> bool:
    """规则2：存在连续 6 个相邻 G/C（滑窗扫描实现）。"""
    n = len(seq)
    for i in range(n - RUN_GC + 1):
        if all(b in "GC" for b in seq[i:i + RUN_GC]):
            return True
    return False


def fail_run5_same(seq: str) -> bool:
    """规则3：存在连续 5 个相同碱基（滑窗扫描实现）。"""
    n = len(seq)
    for i in range(n - RUN_SAME + 1):
        if len(set(seq[i:i + RUN_SAME])) == 1:
            return True
    return False


def fail_palindrome(seq: str) -> bool:
    """规则4：OligoFormer 回文判定。

    语义（与 scripts/infer.py palindromic_sequence 逐字符等价）：
    对起点 i ∈ [0, n-8]，取 4 nt 臂 seq[i:i+4]，若其反向互补序列出现在
    seq[i+4:]（下游 ≥4 nt，允许任意间隔/紧邻）即判为回文。等价于 19 nt 内
    存在 ≥8 nt 的反向重复（可形成发夹/自互补茎），臂长固定 4 nt。
    """
    for i in range(len(seq) - 2 * PAL_ARM + 1):
        pat = rc(seq[i:i + PAL_ARM])
        if pat in seq[i + PAL_ARM:]:
            return True
    return False


# ---------------------------------------------------------------------------
# 四、单序列注释
# ---------------------------------------------------------------------------

def evaluate(seq: str) -> dict:
    """对一条 19 nt 引导链做四项规则注释。

    返回字段：gc_pct / fail_gc_range / fail_run6_gc / fail_run5_same /
    fail_palindrome / seq_pass / hit_code / hit_rule。
    hit_code 与 hit_rule 按 func_filter 源码的判定顺序（先命中者优先）标注。
    """
    s = normalize(seq)
    gc = gc_pct(s)
    flags = {
        "fail_gc_range": 1 if gc < GC_MIN or gc > GC_MAX else 0,
        "fail_run6_gc": 1 if fail_run6_gc(s) else 0,
        "fail_run5_same": 1 if fail_run5_same(s) else 0,
        "fail_palindrome": 1 if fail_palindrome(s) else 0,
    }
    code, rule = 0, "pass"
    for flag, c, name in HIT_ORDER:              # func_filter 判定顺序
        if flags[flag]:
            code, rule = c, name
            break
    return {
        "gc_pct": round(gc, 2),
        **flags,
        "seq_pass": 1 if code == 0 else 0,
        "hit_code": code,
        "hit_rule": rule,
    }


# ---------------------------------------------------------------------------
# 五、输入读取与候选行归一化
# ---------------------------------------------------------------------------

def read_csv_rows(path: Path) -> list[dict]:
    """读 UTF-8-sig CSV 为行字典列表。"""
    if not path.exists():
        raise FileNotFoundError(f"找不到输入文件：{path}")
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh)]


def build_candidates(mut_rows: list[dict], wt_rows: list[dict]) -> list[dict]:
    """把任务 6（突变）与任务 5（完全互补）两表归一化为统一注释行骨架。

    统一列见 OUT_FIELDS；kind = 'mut'/'wt'；guide_checked 为被筛序列
    （mut -> guide_mut_19，wt -> guide_antisense_19），后续 evaluate 用它。
    """
    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()

    for r in mut_rows:                                   # ---- 突变体（任务6）----
        out = {f: ("" if r.get(f) is None else r[f]) for f in OUT_FIELDS}
        guide = normalize(r.get("guide_mut_19", ""))
        out.update({
            "window_id": r.get("window_id", ""),
            "kind": "mut",
            "cds_start": r.get("cds_start", ""),
            "cds_end": r.get("cds_end", ""),
            "guide_checked": guide,
        })
        key = ("mut", out["variant_id"] or f"{out['window_id']}|{guide}")
        if key in seen:
            raise ValueError(f"任务6表存在重复记录：{key[1]}")
        seen.add(key)
        rows.append(out)

    for r in wt_rows:                                    # ---- 完全互补（任务5）----
        out = {f: ("" if r.get(f) is None else r[f]) for f in OUT_FIELDS}
        guide = normalize(r.get("guide_antisense_19", ""))
        out.update({
            "window_id": r.get("window_id", ""),
            "kind": "wt",
            "cds_start": r.get("cds_start", ""),
            "cds_end": r.get("cds_end", ""),
            "nm_003012_start": r.get("nm_003012_start", ""),
            "nm_003012_end": r.get("nm_003012_end", ""),
            "target_mRNA_19": r.get("target_mRNA_19", ""),
            "guide_wt_19": guide,
            "guide_checked": guide,
            "sense_strand_19": r.get("sense_strand_19", ""),
            "window_gc_pct": r.get("window_gc_pct", ""),
        })
        key = ("wt", out["window_id"] or guide)
        if key in seen:
            raise ValueError(f"任务5表存在重复窗口：{key[1]}")
        seen.add(key)
        rows.append(out)
    return rows


def _sort_key(row: dict) -> tuple:
    """输出排序：按窗口编号；同窗口 WT 在前，随后 v01..v15。"""
    win = row["window_id"]
    try:
        win_no = int(re.sub(r"\D", "", win) or 0)
    except ValueError:
        win_no = 0
    if row["kind"] == "wt":
        return (win_no, 0, 0)
    var_no = int(re.sub(r"\D", "", row.get("variant_no", "v00")) or 0)
    return (win_no, 1, var_no)


# ---------------------------------------------------------------------------
# 六、输出与统计
# ---------------------------------------------------------------------------

def write_csv(path: Path, rows: list[dict]) -> int:
    """UTF-8-sig（带 BOM）写 CSV，字段顺序固定为 OUT_FIELDS。"""
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: ("" if v is None else v) for k, v in r.items()})
    return len(rows)


# ---------------------------------------------------------------------------
# 七、校验：单元锚点 + OligoFormer func_filter 逐条对拍（独立实现）
# ---------------------------------------------------------------------------

def _ref_func_filter_label(seq: str) -> int:
    """OligoFormer scripts/infer.py func_filter 的**逐字符忠实复刻**（单条）。

    返回 0=通过 / 1=GC越界 / 2=连续5相同 / 3=连续6个G·C / 4=回文。
    """
    s = seq
    # 1) GC 含量
    gc = (s.count("G") + s.count("C")) / len(s) * 100
    if gc < 30 or gc > 65:
        return 1
    # 2) 连续 ≥5 个相同碱基
    for p in _HOMO_PATTERNS:
        if re.search(p, s):
            return 2
    # 3) 连续 ≥6 个相邻 G/C（原文 itertools.product 枚举 64 种六联串）
    for p in _GC6_PATTERNS:
        if re.search(p, s):
            return 3
    # 4) 回文（原文窗口 i ∈ range(len-8+1)，臂 4 nt）
    for i in range(len(s) - 8 + 1):
        pat = s[i:i + 4][::-1].translate(str.maketrans("AUCG", "UAGC"))
        if re.search(pat, s[i + 4:]):
            return 4
    return 0


def _alt_flags(seq: str) -> tuple:
    """四项规则的**第二套独立表述**（正则/枚举），用于交叉校验主实现。"""
    s = normalize(seq)
    f_gc = gc_pct(s) < GC_MIN or gc_pct(s) > GC_MAX
    f_run6 = bool(re.search(r"[GC]{%d}" % RUN_GC, s))
    f_homo = bool(re.search(r"(.)\1{%d}" % (RUN_SAME - 1), s))
    # 回文：枚举所有“左臂 i、右臂 j（j ≥ i+4）”配对，任一满足即真
    f_pal = any(
        rc(s[i:i + PAL_ARM]) == s[j:j + PAL_ARM]
        for i in range(len(s) - PAL_ARM + 1)
        for j in range(i + PAL_ARM, len(s) - PAL_ARM + 1)
    )
    return f_gc, f_run6, f_homo, f_pal


def _parity_seq_pool(real_seqs: list[str], seed: int = 20250107) -> list[str]:
    """对拍样本池：全部真实序列 + 固定种子随机 19 nt（覆盖各规则组合）。"""
    rng = random.Random(seed)
    pool = list(real_seqs)
    alphabet = "AUGC"
    for _ in range(12000):                      # 均匀随机
        pool.append("".join(rng.choice(alphabet) for _ in range(19)))
    for _ in range(4000):                       # GC 富集随机（多触发规则1/2）
        pool.append("".join(rng.choices(alphabet, weights=[1, 1, 4, 4], k=19)))
    # 定向构造：保证每种规则/组合都被覆盖（全部 19 nt）
    pool += [
        "GCGCGCGCGCGCGCGCGCG",                  # 连续6个G·C(规则2)，无5同(规则3不触发)
        "ACGUAAAAAUCGCGCAUAG",                  # 连续5个A(规则3)
        "ACGUACGU" + "AAGCAUGC" + "CGC",        # 臂ACGU·回文：rc(ACGU)=ACGU 出现在下游4nt(规则4)
        "GC" * 7 + "AUGCA",                     # GC≈84%>65%(规则1)，兼规则2
        "AUCGAUCGAUCGAUCGAUC",                  # 周期4 自互补(规则4)，人工对照
        "AAUUCCGGAACCUUAACCU",                  # 区块回文对照（rc(UUCC)=GGAA 下游命中）
    ]
    return pool


def run_selftests(real_rows: list[dict] | None = None, seed: int = 20250107) -> None:
    """内置校验：
    1) 单元锚点（逐规则语义）；2) 与 func_filter 忠实复刻对拍 hit_code；
    3) 与第二套独立实现交叉校验四个 fail 标志；4) 真实行结构自检。
    """
    # ---- 单元锚点（构造串长度均为 19；规则 1..4 各自至少一条正向命中）----
    assert fail_run6_gc("GC" * 9 + "G"), "规则2 锚点失败"                       # 全 G/C，六联 G·C
    assert not fail_run6_gc("ACGUACGUACGUACGUAUC"), "规则2 误报"               # 含 A/U，无六联 G·C
    assert fail_run5_same("ACGUAAAAAUCGCGCAUAG"), "规则3 锚点失败"             # 连续 5 个 A
    assert not fail_run5_same("ACGUACGUACGUACGUAUC"), "规则3 误报"             # 无相邻重复
    assert fail_palindrome("ACGUACGU" + "AAGCAUGC" + "CGC"), "规则4 锚点失败"   # rc(ACGU)=ACGU 邻接
    assert fail_gc_range("GC" * 7 + "AUGCA"), "规则1 锚点失败"                 # GC≈84%>65%
    assert not fail_gc_range("ACGUACGUACGUACGUAUC"), "规则1 误报"              # GC≈47%∈[30,65]
    # 说明：回文“不命中”不做手工锚点（反向重复极易被漏看），交由下方与
    # func_filter 忠实复刻 + 第二套独立实现的万级随机序列对拍覆盖。

    # ---- 与 func_filter 对拍（样本池含全部真实候选 + 随机 + 定向构造）----
    pool = []
    if real_rows:
        pool.extend(r["guide_checked"] for r in real_rows)
    pool = _parity_seq_pool(pool, seed=seed)
    bad = []
    for s in pool:
        ev = evaluate(s)
        ref = _ref_func_filter_label(normalize(s))
        if ref != ev["hit_code"]:
            bad.append((s, ref, ev["hit_code"]))
            if len(bad) > 10:
                break
    if bad:
        raise AssertionError(f"与 OligoFormer func_filter 对拍不一致 {len(bad)} 例，"
                             f"前3例：{bad[:3]}")
    # seq_pass 与 hit_code 的自洽
    for s in pool[:5000]:
        ev = evaluate(s)
        assert (ev["seq_pass"] == 1) == (ev["hit_code"] == 0), s

    # ---- 与第二套独立实现交叉校验 ----
    for s in pool:
        ev = evaluate(s)
        f_gc, f_run6, f_homo, f_pal = _alt_flags(s)
        assert ev["fail_gc_range"] == (1 if f_gc else 0), ("gc", s)
        assert ev["fail_run6_gc"] == (1 if f_run6 else 0), ("run6", s)
        assert ev["fail_run5_same"] == (1 if f_homo else 0), ("homo5", s)
        assert ev["fail_palindrome"] == (1 if f_pal else 0), ("pal", s)

    print(f"[校验] 全部通过：样本池 {len(pool)} 条（真实 {len(real_rows) if real_rows else 0}"
          f" + 随机/定向构造）与 OligoFormer func_filter 逐条对拍一致，"
          f"四项规则与第二套独立实现交叉校验一致。")


def verify_real_rows(rows: list[dict]) -> None:
    """对真实候选行做结构与一致性自检。"""
    n_wt = sum(1 for r in rows if r["kind"] == "wt")
    n_mut = sum(1 for r in rows if r["kind"] == "mut")
    assert n_wt + n_mut == len(rows)
    for r in rows:
        s = r["guide_checked"]
        assert len(s) == 19 and all(c in "AUGC" for c in s), r["window_id"]
        if r["kind"] == "wt":
            assert r["guide_checked"] == r["guide_wt_19"]
            assert rc(r["target_mRNA_19"]) == s, f"{r['window_id']} WT≠rc(target)"
        else:
            assert r["guide_checked"] == r["guide_mut_19"]
            # 突变体仅由任务6在 5 个位点做单点替换：与 WT 恰差 1 nt
            wt = r["guide_wt_19"]
            diff = [i for i, (a, b) in enumerate(zip(wt, s)) if a != b]
            assert len(diff) == 1, f"{r['variant_id']} 非单点突变"
            pos = int(r["position"]) if str(r["position"]).isdigit() else -1
            assert diff[0] == pos - 1, f"{r['variant_id']} 突变位点不符"
            assert rc(r["target_mRNA_19"]) == wt, f"{r['variant_id']} WT≠rc(target)"
    print(f"[校验] 真实候选结构自检通过：WT {n_wt} 条 + 突变 {n_mut} 条 = {len(rows)} 条。")


# ---------------------------------------------------------------------------
# 八、主流程
# ---------------------------------------------------------------------------

def default_task6_path() -> Path:
    return _GEN_DIR / TASK6_DEFAULT_NAME


def default_task5_path() -> Path:
    return _GEN_DIR / TASK5_DEFAULT_NAME


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="任务7：序列规则筛选——四项规则串联过滤（参考 OligoFormer func_filter）")
    ap.add_argument("--library", type=Path, default=None,
                    help=f"任务6突变库 CSV（默认自动定位 {default_task6_path()}）")
    ap.add_argument("--wt", type=Path, default=None,
                    help=f"任务5完全互补表 CSV（默认自动定位 {default_task5_path()}）")
    ap.add_argument("--no-wt", action="store_true",
                    help="不并入任务5完全互补引导链，只筛任务6突变库")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="输出目录（默认 初级规则筛选 文件夹）")
    ap.add_argument("--prefix", type=str, default=PREFIX_DEFAULT,
                    help="输出文件名前缀（默认 %s）" % PREFIX_DEFAULT)
    ap.add_argument("--selftest", action="store_true",
                    help="只跑内置校验（单元锚点 + func_filter 对拍）后退出，不读写文件")
    return ap


def main(argv=None) -> int:
    args = build_argparser().parse_args(argv)

    # ---- 纯自检模式：不依赖任何输入文件，只验证规则实现与对拍 ----
    if args.selftest:
        run_selftests(None)
        return 0

    # ---- 输入定位 ----
    library = args.library or default_task6_path()
    if not library.exists():
        raise SystemExit(f"错误：找不到任务6突变库 {library}。\n"
                         f"请先在 项目搭建\\候选序列生成 运行 task6_site_mutagenesis.py，"
                         f"或 --library 指定路径。")
    wt_path = args.wt
    if wt_path is not None and not wt_path.exists():
        raise SystemExit(f"错误：--wt 指定的文件不存在：{wt_path}")
    include_wt = not args.no_wt and (args.wt is not None or default_task5_path().exists())

    mut_rows = read_csv_rows(library)
    wt_rows = read_csv_rows(default_task5_path()) if include_wt else []
    n_mut_in, n_wt_in = len(mut_rows), len(wt_rows)
    print(f"[输入] 任务6突变库：{library}  ({n_mut_in} 行)")
    print(f"[输入] 任务5 WT    ：{'并入 ' + str(default_task5_path()) if include_wt else '跳过（--no-wt）'}  "
          f"({'%d 行' % n_wt_in if include_wt else ''})")

    # ---- 归一化 + 四项规则注释 ----
    cands = build_candidates(mut_rows, wt_rows)
    if not cands:
        raise SystemExit("错误：候选池为空。")
    for r in cands:
        r.update(evaluate(r["guide_checked"]))
    rows = sorted(cands, key=_sort_key)

    # ---- 校验 ----
    verify_real_rows(rows)
    run_selftests(rows)

    # ---- 输出 ----
    out_dir = (args.out_dir or _HERE)
    out_dir.mkdir(parents=True, exist_ok=True)
    annotated_path = out_dir / (args.prefix + "_annotated.csv")
    passed_path = out_dir / (args.prefix + "_passed.csv")
    passed = [r for r in rows if r["seq_pass"] == 1]
    write_csv(annotated_path, rows)
    write_csv(passed_path, passed)
    print(f"[输出] 注释全表：{annotated_path} 共 {len(rows)} 行")
    print(f"[输出] 通过子表：{passed_path} 共 {len(passed)} 行")

    # ---- 统计 ----
    from collections import Counter
    n_total = len(rows)
    n_pass = len(passed)
    n_rej = n_total - n_pass
    rule_hits = {flag: sum(1 for r in rows if r[flag] == 1)
                 for flag, _, _ in HIT_ORDER}
    code_c = Counter(r["hit_code"] for r in rows)
    n_wt = sum(1 for r in rows if r["kind"] == "wt")
    n_mut = n_total - n_wt
    print(f"\n[统计] 总候选 {n_total} 条 = WT {n_wt} + 突变 {n_mut}")
    print(f"[统计] 通过 {n_pass} 条（{n_pass / n_total * 100:.1f}%），"
          f"淘汰 {n_rej} 条（{n_rej / n_total * 100:.1f}%）")
    print("[统计] 各规则独立命中数（同一序列可重叠计多次）：")
    names = {"fail_gc_range": "1.GC越界[30,65]之外",
             "fail_run6_gc": "2.连续≥6个相邻G/C",
             "fail_run5_same": "3.连续≥5个相同碱基",
             "fail_palindrome": "4.回文序列"}
    for flag, _, _ in HIT_ORDER:
        print(f"        {names[flag]:<14} {rule_hits[flag]:>6} 条")
    print("[统计] 按 func_filter 首命中原因（互斥划分，0=通过）：")
    for code in range(5):
        print(f"        func_filter={code} ({HIT_CODE_NAMES.get(code, 'pass'):<10}) "
              f"{code_c.get(code, 0):>6} 条")

    # ---- 样例 ----
    print("\n[样例] 前 3 条淘汰候选：")
    for r in [x for x in rows if x["seq_pass"] == 0][:3]:
        print(f"  {r['window_id']} {r.get('variant_id') or '(WT)':<22} "
              f"guide={r['guide_checked']}  gc={r['gc_pct']}%  "
              f"hit={r['hit_rule']}({r['hit_code']})")
    print("[样例] 前 3 条通过候选：")
    for r in passed[:3]:
        print(f"  {r['window_id']} {r.get('variant_id') or '(WT)':<22} "
              f"guide={r['guide_checked']}  gc={r['gc_pct']}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
