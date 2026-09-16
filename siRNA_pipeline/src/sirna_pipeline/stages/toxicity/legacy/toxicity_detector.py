# -*- coding: utf-8 -*-
"""
================================================================================
模块名称：siRNA 毒性检测模块 —— 流水线任务 10
================================================================================
功能说明
--------
    对通过任务 8（结构硬过滤）、任务 9（脱靶检测）的候选做两类相互独立、
    分列输出的毒性评估（默认只打分/标记，供任务 13 综合排序作软惩罚项）：

    A) seed 介导细胞活力毒性（子模块 A，1:1 复刻 OligoFormer 仓库 -tox 行为）
       引导链第 2–7 位六聚体 seed（OligoFormer 代码为 0 基切片 siRNA[1:7]）
       → 查 4096-hexamer 静态活力表 data/cell_viability.txt
       → cell_viability < viability_threshold(默认 50.0) 标记为毒
       （score 越高越安全；"低 = 活力低 = 毒"为证据倾向口径，见模块注释与
        《任务10_调研证据清单.md》§3.5，验收期与 NW_CellToxicityDB 交叉比对）

    B) 免疫刺激 motif 扫描（子模块 B）
       精确匹配  UGUGU / GUCCUUCAA        → 高风险（Judge 2005 等）
       poly-U   连续 ≥4 U → 高风险；连续 3 U → 提示
       U 含量   单链 U 含量 ≥40% → 提示（U-rich 型，TLR7/8 证据）
       → 分级为工程设定（文献无统一"连续 N 个 U"阈值），阈值全部可配

序列约定
--------
    guide      ：引导链(反义链) 19 nt，5'→3'（A/U/C/G，T 自动转 U）
    passenger  ：过客链(正义链)，默认按 passenger = rc(guide) 重建
                 （假设：合成双链两链完全互补；待与任务 6 产物格式核对）
    seed       ：guide[1:7]（g2–g7 六聚体），与 OligoFormer slice(1,7) 一致

关键设计事实（写代码时利用）
--------
    * 错配位点 g1/g12/g17/g18/g19 全部落在 seed(g2–g7) 之外
      ⇒ 同一 mRNA 窗口的 17 条变体（任务 6 产物）seed 相同 ⇒ 子模块 A 按
      窗口只查一次表、其余变体广播（detect_batch 内置 seed 缓存，
      tox_broadcast 字段标记），该恒定性同时是任务 16 消融的可讨论点。
    * 免疫 motif 命中随变体序列变化 ⇒ 子模块 B 必须逐条扫描。

受体分层与产物质控（不参与打分）
--------
    TLR7/8（GU/U-rich、内体）为序列主风险；RIG-I/5'-三磷酸仅在产物质控层
    提示：体外化学合成成熟双链为 5'-OH/5'-PO4，不构成 RIG-I 配体
    （Hornung 2006 Science PMID 17038590），勿误用 T7 体外转录 5'-ppp 产物。

阈值/表来源（写入代码即证据出处）
--------
    [A] cell_viability.txt        ：OligoFormer 仓库 toxicity/cell_viability.txt
                                    （github.com/lulab/OligoFormer，学术用途）
    [A] seed 切片与阈值 50.0      ：OligoFormer scripts/infer.py（仓库行为复刻，
                                    论文 Bioinformatics 2024;40(10):btae577 无该方法描述）
    [B] UGUGU / GUCCUUCAA         ：Judge 2005 Nat Biotechnol 23(4):457-462,
                                    PMID 15778705（UGUGU 关联逐字级；
                                    GUCCUUCAA 为转述级引用）
    [B] poly-U / U-rich → TLR7/8  ：Heil 2004 Science PMID 14976262;
                                    Diebold 2004 PMID 14976261; Sioud 2005
                                    PMID 15854645; Forsbach 2008 PMID 18322178
    [B] poly-U/UC → RIG-I(5'-ppp 背景)：Saito 2008 Nature PMID 18548002
    [B] 2'-O-methyl 缓解          ：Judge 2006 Mol Ther PMID 16343994;
                                    Robbins 2007 PMID 17579574（供任务 19 引用）

输入输出
--------
    detect(guide, passenger=None, window_id=None, variant_id=None) -> dict
    detect_batch(records) -> list[dict]        # records: {guide, ...} 或 [str]
    filter_candidates(records, hard_rules=None) -> (kept, rejected)
    get_resources_status() / describe()
    输出字段分 tox_*（A）与 imm_*（B）两组；硬过滤开关 hard_rules（默认关闭）。

运行测试
--------
    python toxicity_detector.py              # 内置 8 用例自检（无外部依赖）
    python toxicity_detector.py --trial      # SFRP1 mRNA 滑窗试跑统计
================================================================================
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys

# ==============================================================================
# 一、常量
# ==============================================================================

VALID_BASES = set("ACGU")
COMP = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}

# 默认数据表路径（模块同级 data/ 或模块同级；可被环境变量覆盖）
_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE_CANDIDATES = [
    os.environ.get("TOX_VIABILITY_TABLE", ""),
    os.path.join(_MODULE_DIR, "data", "cell_viability.txt"),
    os.path.join(_MODULE_DIR, "cell_viability.txt"),
]
DEFAULT_TABLE_CANDIDATES = [p for p in DEFAULT_TABLE_CANDIDATES if p]

# 免疫刺激精确匹配 motif（key 一律为大写 RNA）
# level: high —— 与文献精确匹配即高风险
IMMUNE_MOTIFS = {
    "UGUGU": {
        "level": "high",
        "ref": "Judge 2005 Nat Biotechnol 23(4):457-462, PMID 15778705 "
               "(UGUGU 关联为逐字级；GUCCUUCAA 为转述级)",
    },
    "GUCCUUCAA": {
        "level": "high",
        "ref": "Judge 2005（转述级，未逐字核验；与 UGUGU 在工具/综述中并列引用）",
    },
}
MOTIF_ORDER = ["UGUGU", "GUCCUUCAA"]

# 输出字段顺序（TSV 写出用）
OUTPUT_FIELDS = [
    "id", "window_id", "variant_id", "guide", "passenger", "length_warn",
    # 子模块 A
    "tox_seed", "tox_viability_score", "tox_viability_flag", "tox_viability_miss",
    "tox_broadcast",
    # 子模块 B
    "imm_ugugu", "imm_guccuucaa", "imm_polyu_max_run_guide", "imm_polyu_max_run_passenger",
    "imm_polyu_high", "imm_polyu_warn", "imm_u_content_guide", "imm_u_content_passenger",
    "imm_u_content_warn", "imm_high_flag", "imm_warn_flag", "imm_flag",
    "imm_hits_detail", "imm_receptor",
    # 过滤
    "hard_reject", "reject_reason", "ok",
]

# ==============================================================================
# 二、基础工具
# ==============================================================================


def normalize(seq: str) -> str:
    """去空白、大写、T->U。"""
    return "".join(str(seq).split()).upper().replace("T", "U")


def rc(seq: str) -> str:
    """反向互补（RNA）。"""
    return "".join(COMP[b] for b in reversed(normalize(seq)))


def find_overlapping(seq: str, motif: str):
    """返回 motif 在 seq 中的所有重叠起始位点（0 基）。"""
    out = []
    pat = re.compile("(?=" + motif + ")")
    for m in pat.finditer(seq):
        out.append(m.start())
    return out


def esc_ref(ref: str) -> str:
    """把出处串压成单行（注释/输出友好）。"""
    return " ".join(ref.split())


# ==============================================================================
# 三、资源表读取
# ==============================================================================

def load_viability_table(path: str | None = None) -> dict:
    """读取 cell_viability.txt → {hexamer6: float}。

    复刻对象：OligoFormer toxicity/cell_viability.txt
    （4096 行 = 全部 4^6 hexamer；列 Seed / cell_viability）
    """
    if path is None:
        for cand in DEFAULT_TABLE_CANDIDATES:
            if cand and os.path.isfile(cand):
                path = cand
                break
    if path is None or not os.path.isfile(path):
        raise FileNotFoundError(
            "未找到 cell_viability.txt；请将 OligoFormer 的 "
            "toxicity/cell_viability.txt 放入模块 data/ 目录，"
            "或设置环境变量 TOX_VIABILITY_TABLE。"
        )
    table: dict = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or "Seed" not in header or "cell_viability" not in header:
            raise ValueError(f"表头不符合预期: {header!r}")
        for row in reader:
            if not row or len(row) < 2:
                continue
            seed = normalize(row[0])
            if len(seed) != 6:
                raise ValueError(f"seed 长度异常: {row[0]!r}")
            table[seed] = float(row[1])
    if len(table) != 4096:
        # 不强制 4096：仅告警（表来自上游时允许子集），但应基本覆盖
        print(f"[warn] cell_viability 表仅 {len(table)}/4096 条，请核对数据源。",
              file=sys.stderr)
    return table


# ==============================================================================
# 四、ToxicityDetector 主类
# ==============================================================================

class ToxicityDetector:
    """siRNA 毒性检测（seed 活力查表 + 免疫刺激 motif 扫描）。"""

    def __init__(
        self,
        viability_table: dict | str | None = None,
        viability_threshold: float = 50.0,
        polyu_run_high: int = 4,
        polyu_run_warn: int = 3,
        u_content_warn: float = 0.40,
        hard_rules: dict | None = None,
        scan_passenger: bool = True,
        require_len: int = 19,
    ):
        # 子模块 A 参数（复刻 OligoFormer 默认：threshold=50.0）
        self.viability_threshold = float(viability_threshold)
        # 子模块 B 参数（工程设定，文献无统一阈值——注释中已声明）
        self.polyu_run_high = int(polyu_run_high)
        self.polyu_run_warn = int(polyu_run_warn)
        self.u_content_warn = float(u_content_warn)
        self.scan_passenger = bool(scan_passenger)
        self.require_len = int(require_len)

        # 硬过滤规则（默认关闭；启用后 filter_candidates/detect 按规则硬杀）
        # 支持键：
        #   viability_below: float  —— tox_viability_score < 该值即杀
        #   polyu_run_high : bool   —— 出现连续 ≥polyu_run_high 个 U 即杀
        #   imm_high       : bool   —— 出现任一高风险免疫命中（精确 motif / polyU 高）即杀
        self.hard_rules = dict(hard_rules or {})

        # 表与缓存
        if isinstance(viability_table, dict):
            self._table: dict | None = viability_table
        elif isinstance(viability_table, str):
            self._table = load_viability_table(viability_table)
        else:
            self._table = None  # 首次使用时按默认路径加载
        self._viability_cache: dict = {}
        self.cache_hits = 0        # 命中缓存的条数
        self.cache_misses = 0      # 首次查表条数

    # ------------------------------------------------------------------ 资源
    def _ensure_table(self) -> dict:
        if self._table is None:
            self._table = load_viability_table(None)
        return self._table

    def get_resources_status(self) -> dict:
        try:
            n = len(self._ensure_table())
        except Exception as exc:  # noqa: BLE001
            return {"viability_table": "missing", "detail": str(exc)}
        return {"viability_table": "ok", "entries": n,
                "threshold": self.viability_threshold}

    # ------------------------------------------------------------------ 工具
    def _validate(self, guide: str) -> str:
        g = normalize(guide)
        bad = sorted(set(g) - VALID_BASES)
        if bad:
            raise ValueError(f"引导链含非法碱基 {bad}: {guide!r}")
        return g

    @staticmethod
    def _seed(guide: str) -> str:
        """OligoFormer 语义：siRNA[1:7]（0 基 1:7）= g2–g7 六聚体。"""
        return guide[1:7]

    @staticmethod
    def _max_u_run(seq: str) -> int:
        runs = re.findall(r"U{2,}", seq)
        return max((len(r) for r in runs), default=0)

    @staticmethod
    def _u_content(seq: str) -> float:
        if not seq:
            return 0.0
        return seq.count("U") / len(seq)

    # ------------------------------------------------------------------ 扫描
    def _scan_strand(self, seq: str, label: str) -> list:
        """对单条链做精确 motif 扫描 → [{motif, strand, pos(1基), level}]。"""
        hits = []
        for motif in MOTIF_ORDER:
            for start in find_overlapping(seq, motif):
                hits.append({
                    "motif": motif,
                    "strand": label,
                    "pos": start + 1,  # 1-based
                    "level": IMMUNE_MOTIFS[motif]["level"],
                })
        return hits

    def _immune_eval(self, guide: str, passenger: str | None) -> dict:
        """返回 imm_* 语义字段（不含 tox_*）。"""
        res = {m.lower(): False for m in MOTIF_ORDER}          # imm_ugugu/guccuucaa
        # 每条链只扫描一次（_scan_strand 内部遍历全部 motif，避免重复命中）
        hits_all: list = list(self._scan_strand(guide, "guide"))
        if self.scan_passenger and passenger:
            hits_all += self._scan_strand(passenger, "passenger")
        for h in hits_all:
            res[h["motif"].lower()] = True
        res["hits"] = hits_all

        max_run_g = self._max_u_run(guide)
        max_run_p = self._max_u_run(passenger) if passenger else 0
        polyu_high = max_run_g >= self.polyu_run_high or max_run_p >= self.polyu_run_high
        polyu_warn = (not polyu_high) and (
            max_run_g >= self.polyu_run_warn or max_run_p >= self.polyu_run_warn)
        u_g = self._u_content(guide)
        u_p = self._u_content(passenger) if passenger else 0.0
        u_warn = max(u_g, u_p) >= self.u_content_warn

        res["polyu_max_run_guide"] = max_run_g
        res["polyu_max_run_passenger"] = max_run_p
        res["polyu_high"] = polyu_high
        res["polyu_warn"] = polyu_warn
        res["u_content_guide"] = round(u_g, 4)
        res["u_content_passenger"] = round(u_p, 4)
        res["u_content_warn"] = u_warn
        return res

    # ------------------------------------------------------------------ 主接口
    def detect(
        self,
        guide: str,
        passenger: str | None = None,
        window_id=None,
        variant_id=None,
        rid=None,
        use_cache: bool = False,
    ) -> dict:
        g = self._validate(guide)
        if passenger is not None:
            p = self._validate(passenger)
        else:
            p = rc(g)

        length_warn = (self.require_len is not None and len(g) != self.require_len)

        # ---- 子模块 A：seed 查表
        seed = self._seed(g)
        table = self._ensure_table()
        cached = seed in self._viability_cache if use_cache else False
        if use_cache and cached:
            score = self._viability_cache[seed]
            self.cache_hits += 1
        else:
            score = table.get(seed)
            if use_cache:
                self._viability_cache[seed] = score
            self.cache_misses += 1
        miss = score is None
        flag_toxic = (not miss) and (score < self.viability_threshold)

        # ---- 子模块 B：免疫扫描
        imm = self._immune_eval(g, p)

        # 汇总免疫标记
        motif_hit = any(imm[m.lower()] for m in MOTIF_ORDER)
        imm_high_flag = motif_hit or imm["polyu_high"]
        imm_warn_flag = (imm["polyu_warn"] or imm["u_content_warn"]) and not imm_high_flag
        imm_flag = imm_high_flag or imm_warn_flag

        detail_parts = []
        for h in sorted(imm["hits"], key=lambda x: (x["motif"], x["strand"], x["pos"])):
            detail_parts.append(f"{h['motif']}@{h['strand']}:{h['pos']}")
        hits_detail = ";".join(detail_parts) or "-"

        # ---- 硬过滤（默认关闭）
        reject = False
        reasons = []
        if "viability_below" in self.hard_rules and not miss \
                and score < float(self.hard_rules["viability_below"]):
            reject = True
            reasons.append(f"viability<{self.hard_rules['viability_below']}")
        if self.hard_rules.get("polyu_run_high") and imm["polyu_high"]:
            reject = True
            reasons.append("polyU高")
        if self.hard_rules.get("imm_high") and imm_high_flag:
            reject = True
            reasons.append("免疫高风险")

        row = {
            "id": rid,
            "window_id": window_id,
            "variant_id": variant_id,
            "guide": g,
            "passenger": p,
            "length_warn": length_warn,
            # A
            "tox_seed": seed,
            "tox_viability_score": None if miss else score,
            "tox_viability_flag": flag_toxic,
            "tox_viability_miss": miss,
            "tox_broadcast": cached,
            # B
            "imm_ugugu": imm["ugugu"],
            "imm_guccuucaa": imm["guccuucaa"],
            "imm_polyu_max_run_guide": imm["polyu_max_run_guide"],
            "imm_polyu_max_run_passenger": imm["polyu_max_run_passenger"],
            "imm_polyu_high": imm["polyu_high"],
            "imm_polyu_warn": imm["polyu_warn"],
            "imm_u_content_guide": imm["u_content_guide"],
            "imm_u_content_passenger": imm["u_content_passenger"],
            "imm_u_content_warn": imm["u_content_warn"],
            "imm_high_flag": imm_high_flag,
            "imm_warn_flag": imm_warn_flag,
            "imm_flag": imm_flag,
            "imm_hits_detail": hits_detail,
            "imm_receptor": "TLR7/8主风险；RIG-I/5'-ppp仅产物质控(Hornung 2006 PMID 17038590)",
            # 过滤
            "hard_reject": reject,
            "reject_reason": "；".join(reasons) or "",
            "ok": not reject,
        }
        return row

    # ------------------------------------------------------------------ 批量
    def detect_batch(self, records, use_cache: bool = True) -> list:
        """records: list[str]（仅 guide）或 list[dict]（guide 必含）。

        子模块 A 利用 seed 缓存实现窗口级广播：同一 seed 的变体只查一次表，
        后续记录置 tox_broadcast=True（错配位点 g1/g12/g17-19 均在 seed 外）。
        """
        self.cache_hits = 0
        self.cache_misses = 0
        rows = []
        for rec in records:
            if isinstance(rec, str):
                rec = {"guide": rec}
            rows.append(self.detect(
                guide=rec["guide"],
                passenger=rec.get("passenger"),
                window_id=rec.get("window_id"),
                variant_id=rec.get("variant_id"),
                rid=rec.get("id"),
                use_cache=use_cache,
            ))
        return rows

    def filter_candidates(self, records, hard_rules: dict | None = None):
        """应用硬过滤（默认用实例 hard_rules）→ (kept, rejected)。"""
        if hard_rules is not None:
            self.hard_rules = dict(hard_rules)
        rows = self.detect_batch(records)
        kept = [r for r in rows if not r["hard_reject"]]
        rejected = [r for r in rows if r["hard_reject"]]
        return kept, rejected

    # ------------------------------------------------------------------ 说明
    def describe(self) -> str:
        table = self.get_resources_status()
        return (
            "ToxicityDetector 参数与出处\n"
            "--------------------------\n"
            f"资源: {table}\n"
            f"viability_threshold={self.viability_threshold}（OligoFormer 默认 50.0，"
            "仓库行为复刻，score<thr 判毒）\n"
            f"polyu_run_high={self.polyu_run_high}, polyu_run_warn={self.polyu_run_warn}, "
            f"u_content_warn={self.u_content_warn}（工程设定，非文献阈值）\n"
            f"scan_passenger={self.scan_passenger}（passenger=rc(guide) 重建假设）\n"
            f"hard_rules={self.hard_rules or '（关闭，默认软惩罚输出）'}\n"
            "motif 出处: UGUGU/GUCCUUCAA←Judge 2005 PMID 15778705；"
            "poly-U/U-rich→TLR7/8←Heil 2004 PMID 14976262, Diebold 2004 "
            "PMID 14976261, Sioud 2005 PMID 15854645, Forsbach 2008 PMID 18322178；"
            "RIG-I poly-U/UC←Saito 2008 PMID 18548002\n"
        )


# ==============================================================================
# 五、内置自检用例
# ==============================================================================

def _make_guide_with_seed(hexamer: str) -> str:
    """构造引导链，使 g2–g7 = hexamer：A + hexamer + 12 nt（不影响 seed）。"""
    return "A" + hexamer + "A" * 12


def run_self_tests() -> list:
    """返回 [(用例名, 通过?, 说明)]。"""
    det = ToxicityDetector()
    results = []

    # —— 测试序列构造工具 ——
    # poly-A 引导链会使互补 passenger 链 U-rich，误触 poly-U 判据；
    # 因此"干净/边界"用例统一用 C/G 交替填充（双链均无长 U 连串）。
    def _fill_alt(n: int) -> str:
        return "".join("GC"[i % 2] for i in range(n))

    def _run_guide(u_count: int, pre: int = 8) -> str:
        """构造引导链：中间一段连续 U（u_count 个），其余为 C/G 交替。"""
        seq = _fill_alt(pre) + "U" * u_count + _fill_alt(19 - pre - u_count)
        assert len(seq) == 19
        return seq

    g_clean = _fill_alt(19)          # 无 U、无 motif 的干净引导链

    # 用例 1：查表方向 —— 最低分 hexamer 判毒、最高分通过（复刻 OligoFormer）
    table = det._ensure_table()
    low_seed = min(table, key=table.get)
    high_seed = max(table, key=table.get)
    r_low = det.detect(_make_guide_with_seed(low_seed))
    r_high = det.detect(_make_guide_with_seed(high_seed))
    ok1 = (r_low["tox_viability_flag"] is True
           and r_high["tox_viability_flag"] is False
           and r_low["tox_viability_score"] < det.viability_threshold
           and r_high["tox_viability_score"] >= det.viability_threshold)
    results.append(("1.查表方向(低分毒/高分通过)", ok1,
                    f"low={low_seed}({r_low['tox_viability_score']:.1f})→flag="
                    f"{r_low['tox_viability_flag']}; high={high_seed}"
                    f"({r_high['tox_viability_score']:.1f})→flag="
                    f"{r_high['tox_viability_flag']}"))

    # 用例 2：同窗口广播一致（seed 相同、位点 12 不同 → 分数一致且复用）
    g1 = _make_guide_with_seed(high_seed)
    g2 = g1[:11] + ("C" if g1[11] != "C" else "G") + g1[12:]
    assert det._seed(g2) == high_seed
    rows = det.detect_batch([
        {"guide": g1, "window_id": "w1", "variant_id": "v0"},
        {"guide": g2, "window_id": "w1", "variant_id": "v1"},
    ])
    same_score = rows[0]["tox_viability_score"] == rows[1]["tox_viability_score"]
    ok2 = same_score and rows[1]["tox_broadcast"] is True and det.cache_hits == 1
    results.append(("2.窗口级广播(seed复用)", ok2,
                    f"score={rows[0]['tox_viability_score']}, "
                    f"reuse={rows[1]['tox_broadcast']}, cache_hits={det.cache_hits}"))

    # 用例 3/4/6：UGUGU / GUCCUUCAA / passenger 链命中
    g_ugu = "AA" + "UGUGU" + "A" * 12          # 长 19
    r3 = det.detect(g_ugu)
    ok3 = r3["imm_ugugu"] and r3["imm_high_flag"]
    results.append(("3.UGUGU阳性命中", ok3, r3["imm_hits_detail"]))

    g_gucc = "A" * 4 + "GUCCUUCAA" + "A" * 6   # 长 19
    r4 = det.detect(g_gucc)
    ok4 = r4["imm_guccuucaa"] and r4["imm_high_flag"]
    results.append(("4.GUCCUUCAA阳性命中", ok4, r4["imm_hits_detail"]))

    # 阴性对照（C/G 交替，无 motif、双链均无长 U 连串）
    r_neg = det.detect(g_clean)
    ok_neg = (not r_neg["imm_flag"]) and (r_neg["imm_polyu_high"] is False)
    results.append(("5b.阴性对照(无motif)", ok_neg, f"flag={r_neg['imm_flag']}"))

    # 用例 5：poly-U 边界（3U→提示；4U→高风险）
    g3u = _run_guide(3)
    g4u = _run_guide(4)
    r5a = det.detect(g3u)
    r5b = det.detect(g4u)
    ok5 = (r5a["imm_polyu_high"] is False and r5a["imm_polyu_warn"] is True
           and r5b["imm_polyu_high"] is True)
    results.append(("5.polyU边界(3U提示/4U高)", ok5,
                    f"3U: high={r5a['imm_polyu_high']},warn={r5a['imm_polyu_warn']}; "
                    f"4U: high={r5b['imm_polyu_high']}"))

    # 用例 6：passenger 链命中（guide 干净但 rc(guide) 含 UGUGU）
    #   rc(guide) 含 UGUGU ⇔ guide 含 rc(UGUGU)=ACACA
    g_acaca = "AA" + "ACACA" + "A" * 12
    det2_on = ToxicityDetector(scan_passenger=True)
    det2_off = ToxicityDetector(scan_passenger=False)
    r6a = det2_on.detect(g_acaca)
    r6b = det2_off.detect(g_acaca)
    ok6 = (r6a["imm_ugugu"] and r6a["imm_hits_detail"].startswith("UGUGU@passenger")
           and (r6b["imm_ugugu"] is False))
    results.append(("6.passenger互补链扫描(on/off)", ok6,
                    f"on={r6a['imm_hits_detail']}; off={r6b['imm_ugugu']}"))

    # 用例 7：hard_reject 开关
    det_h = ToxicityDetector(hard_rules={"polyu_run_high": True})
    r7a = det_h.detect(g4u)          # polyU 高 → 杀
    r7b = det_h.detect(g_clean)      # 干净 → 过
    ok7 = r7a["hard_reject"] and (not r7b["hard_reject"]) and r7b["ok"]
    results.append(("7.hard_reject开关", ok7,
                    f"4U reject={r7a['hard_reject']}({r7a['reject_reason']}); "
                    f"clean ok={r7b['ok']}"))

    # 用例 8：字段前缀/必含键 + T->U 归一化
    # 构造一个真含 T 的 19mer，检测时须归一化为 U
    seq_with_t = "AU" + "G" * 10 + "TTT" + "A" * 4   # 长度 19
    g8 = det.detect(seq_with_t)
    ok8 = (all(k in g8 for k in OUTPUT_FIELDS)
           and "T" not in g8["guide"]
           and g8["tox_seed"] == g8["guide"][1:7]
           and g8["hard_reject"] is False and g8["ok"] is True)
    results.append(("8.字段/归一化/必含键", ok8,
                    f"guide={g8['guide']}, keys_ok=True"))

    return results


# ==============================================================================
# 六、SFRP1 mRNA 滑窗试跑（示例集成）
# ==============================================================================

def _read_mrna(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return normalize("".join(seq_parts))


def run_trial(mrna_path: str, window_size: int = 19) -> dict:
    """SFRP1 等 mRNA：对每个窗口（完全互补引导链）做毒性检测并统计。"""
    mrna = _read_mrna(mrna_path)
    det = ToxicityDetector()
    guides = [rc(mrna[i:i + window_size])
              for i in range(len(mrna) - window_size + 1)]
    rows = det.detect_batch(guides)
    n = len(rows)
    n_flag_tox = sum(1 for r in rows if r["tox_viability_flag"])
    n_imm_high = sum(1 for r in rows if r["imm_high_flag"])
    n_imm_warn = sum(1 for r in rows if r["imm_warn_flag"])
    n_ugu = sum(1 for r in rows if r["imm_ugugu"])
    n_gucc = sum(1 for r in rows if r["imm_guccuucaa"])
    n_polyu_high = sum(1 for r in rows if r["imm_polyu_high"])
    unique_seeds = len({r["tox_seed"] for r in rows})
    return {
        "mrna_len": len(mrna),
        "n_windows": n,
        "unique_seeds": unique_seeds,
        "n_tox_flag": n_flag_tox,
        "rate_tox": round(n_flag_tox / n, 4) if n else 0.0,
        "n_imm_high": n_imm_high,
        "rate_imm_high": round(n_imm_high / n, 4) if n else 0.0,
        "n_imm_warn": n_imm_warn,
        "n_ugugu": n_ugu,
        "n_guccuucaa": n_gucc,
        "n_polyu_high": n_polyu_high,
        "cache_hits": det.cache_hits,
    }


# ==============================================================================
# 七、命令行入口
# ==============================================================================

def _write_tsv(path: str, rows: list) -> None:
    """按 OUTPUT_FIELDS 顺序写 TSV（None→空串）。"""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(OUTPUT_FIELDS)
        for r in rows:
            w.writerow(["" if r.get(k) is None else r.get(k) for k in OUTPUT_FIELDS])


def main(argv=None):
    ap = argparse.ArgumentParser(description="任务10 毒性检测模块（内置自检 + 试跑）")
    ap.add_argument("--self-test", action="store_true", default=True,
                    help="运行内置 8 用例自检（默认）")
    ap.add_argument("--trial", metavar="MRNA_FASTA", default=None,
                    help="对给定 mRNA 做 19nt 滑窗毒性检测试跑（如数据集/SFRP1-mRNA.txt）")
    ap.add_argument("--input", metavar="CANDIDATES_TSV", default=None,
                    help="批量候选表（列需含 guide；可选 window_id/variant_id）")
    ap.add_argument("--output", metavar="OUT_TSV", default=None,
                    help="批量结果写出路径（配合 --input）")
    args = ap.parse_args(argv)

    if args.self_test:
        results = run_self_tests()
        print("=" * 78)
        print("任务10 毒性检测模块 · 内置自检")
        print("=" * 78)
        failed = 0
        for name, ok, note in results:
            mark = "PASS" if ok else "FAIL"
            if not ok:
                failed += 1
            print(f"[{mark}] {name}  |  {note}")
        print("-" * 78)
        if failed:
            print(f"自检未通过: {failed} 项失败")
            return 1
        print("自检全部通过。")

    if args.trial:
        stats = run_trial(args.trial)
        print("\n" + "=" * 78)
        print(f"SFRP1/mRNA 滑窗试跑统计（完全互补引导链，窗口 19nt）")
        print("=" * 78)
        for k, v in stats.items():
            print(f"  {k}: {v}")
        # 附一句口径说明
        print("说明: tox_viability_flag=seed 活力<thr 判毒(OligoFormer口径)；"
              "imm_high=精确免疫motif或polyU≥4；错配变体库(任务6)未纳入，"
              "A 特征按窗口恒定、B 特征需逐变体重扫。")

    if args.input:
        rows_in = []
        with open(args.input, "r", encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for rec in reader:
                rows_in.append(rec)
        det = ToxicityDetector()
        out_rows = det.detect_batch(rows_in)
        if args.output:
            _write_tsv(args.output, out_rows)
            print(f"批量结果已写出: {args.output}（{len(out_rows)} 条）")
        else:
            print("批量结果（前 5 条关键字段）：")
            for r in out_rows[:5]:
                print(f"  {r['id'] or '-'} seed={r['tox_seed']} "
                      f"score={r['tox_viability_score']} tox={r['tox_viability_flag']} "
                      f"imm_high={r['imm_high_flag']} hits={r['imm_hits_detail']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
