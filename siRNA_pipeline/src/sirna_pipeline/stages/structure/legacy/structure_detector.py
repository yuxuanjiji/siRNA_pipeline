# -*- coding: utf-8 -*-
"""
================================================================================
模块名称：siRNA 结构检测模块（ViennaRNA 套件驱动）—— 流水线任务 8
================================================================================
功能说明
--------
    以 ViennaRNA Package（>=2.5）为第一手段，对候选 siRNA 做三类结构检测并
    依硬阈值淘汰候选：
        1) 引导链自折叠（RNAfold）      —— MFE / 集合自由能 / 质心结构 /
                                            碱基配对概率 / 茎环统计
        2) 末端不对称性（最近邻 NN 表）  —— 5'端 vs 3'端堆叠 ΔG 差，判断
                                            RISC 链选择方向（Khvorova 2003,
                                            Schwarz 2003）
        3) 杂交结构（RNAcofold）        —— 引导链+靶标区杂交 MFE、内环/凸起、
                                            两端配对状态
        4) 靶标可及性（RNAplfold）      —— mRNA 靶标区 19 nt 平均未配对概率
    过滤规则（修订版）见 __init__ 阈值参数：硬淘汰仅限"链选择反转 / 杂交过弱 /
    靶标不可及"；自折叠 MFE 与稳定茎默认降为软提示(warnings)，仅 strict 模式下
    才作硬淘汰；阈值未校准或工具缺失的规则置 None 不误杀。淘汰原因写入
    reject_reason，供任务 16 消融分析。

依赖工具
--------
    * ViennaRNA Python 绑定  import RNA （Level 0）
    * ViennaRNA 命令行       RNAfold / RNAcofold / RNAplfold （Level 1/2）
    * 无 ViennaRNA 时纯 Python 最近邻近似（Level 3，仅引导链 MFE/茎估算，
      cofold/plfold 输出 NaN 并给出降级提示）

输入输出
--------
    detect(sirna_19nt, mrna_57nt=None) -> dict   （每条候选，按工具分组输出）
    detect_batch(sirna_list[, mrna_list]) -> list[dict]
    filter_candidates(sirna_list[, mrna_list]) -> list[dict]（仅通过者）

序列约定
--------
    siRNA 引导链：19 nt，5'→3'（A/U/C/G，T 自动转 U）
    mRNA：57 nt = 5'侧翼19 + 靶标区19(mRNA[19:38]) + 3'侧翼19
    杂交输入：引导链 5'→3' 反平行配对 mRNA 靶标区；guide[p] ↔ target[18-p]

降级分层
--------
    Level 0：Python 绑定可用（RNAfold/RNAcofold 走绑定，RNAplfold 走命令行）
    Level 1：无绑定但三个命令行齐全（全部 subprocess）
    Level 2：仅 RNAfold 可用（cofold/plfold 字段 NaN）
    Level 3：ViennaRNA 完全不可用（RNAfold 用纯 Python NN 近似，其余 NaN，
             打印警告）；get_tool_status() 可查看当前状态

参考文献
--------
    [1] Lorenz R, et al. ViennaRNA Package 2.0. Algorithms Mol Biol. 2011;6:26.
    [2] Khvorova A, et al. Functional siRNAs and miRNAs exhibit strand bias.
        Cell. 2003;115(2):209-216.
    [3] Schwarz DS, et al. Asymmetry in the assembly of the RNAi enzyme complex.
        Cell. 2003;115(2):199-208.
    [4] Ding Y, Lawrence CE. A statistical sampling algorithm for RNA secondary
        structure prediction. Nucleic Acids Res. 2003;31(24):7280-7301.
    [5] Heale BSE, et al. siRNA target site secondary structure predictions
        using local stable substructures. Nucleic Acids Res. 2005;33(1):e30.
    [6] Bernhart SH, et al. Local RNA base pairing probabilities in large
        sequences. Bioinformatics. 2006;22(5):614-615.
    [7] Reynolds A, et al. Rational siRNA design for RNA interference.
        Nat Biotechnol. 2004;22(3):326-330.
    [8] Singh S, et al. (2012) —— 正自由能折叠的 siRNA 更易结合靶标。

运行测试
--------
    python structure_detector.py
    未安装 ViennaRNA 时自动降级并在输出中标注"跳过/降级"；解析器等纯逻辑
    部分用合成输入自检，保证无工具环境下代码逻辑仍可验证。
================================================================================
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter

# ==============================================================================
# 一、常量与基础工具
# ==============================================================================

VALID_BASES = set("ACGU")

# RNA 最近邻堆叠 ΔG 表（37 ℃；Xia/SantaLucia RNA 参数，与热力学参数计算模块
# thermo_calculator.py 的 STACK_DG 完全一致 —— 保证任务 8 与任务 11 字段同口径）
STACK_DG = {
    "AA": -0.93, "AC": -2.24, "AG": -2.08, "AU": -1.10,
    "CA": -2.11, "CC": -3.26, "CG": -2.36, "CU": -2.08,
    "GA": -2.35, "GC": -3.42, "GG": -3.26, "GU": -2.24,
    "UA": -1.33, "UC": -2.35, "UG": -2.11, "UU": -0.93,
}
COMP = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}

# 内部点括号字符
_PAIR_OPEN, _PAIR_CLOSE, _UNPAIRED = "(", ")", "."


def rc(seq: str) -> str:
    """反向互补（RNA）。"""
    return "".join(COMP[b] for b in reversed(seq.upper()))


def normalize(seq: str) -> str:
    """去空白、大写、T->U。"""
    return "".join(seq.split()).upper().replace("T", "U")


def is_comp(a: str, b: str, allow_gu: bool = True) -> bool:
    """a 与 b 是否可配对（WC；allow_gu=True 时含 G:U 摆动）。"""
    if a == COMP[b]:
        return True
    if allow_gu and {a, b} == {"G", "U"}:
        return True
    return False


# ==============================================================================
# 二、主类
# ==============================================================================

class StructureDetector:
    """
    siRNA 结构检测器：RNAfold / RNAcofold / RNAplfold + NN 末端不对称 + 硬阈值。

    阈值默认值依据任务 8 规格；均可通过 __init__ 或属性覆盖。
    """

    # ---- 类属性：ViennaRNA 可用性缓存（首次探测后复用）----
    _binding_ok: bool | None = None
    _cli_cache: dict = {}
    _cached_path: str | None = None
    _cached_bins: str | None = None

    # ==================================================================
    # 构造与工具探测
    # ==================================================================
    def __init__(
        self,
        mfe_threshold: float = -6.0,        # 自折叠 MFE 提示阈值（软提示；仅 strict 时硬淘汰）
        end_diff_threshold: float = -0.5,   # 末端 ΔΔG ≥ 该值才通过（双链端稳定度，硬淘汰）
        stem_threshold: int = 4,            # 最长茎 ≥ 4 bp 判为稳定茎（软提示）
        stem_gc_threshold: float = 60.0,    # 且 GC ≥ 60% 才判为强发夹（软提示）
        cofold_min: float = -25.0,          # （弃用）原"杂交过强"下限，不再作为淘汰依据
        cofold_max: float = -10.0,          # 杂交"过弱"阈值：MFE 须 ≤ 此值（Sfold 结合能≤-10）
        access_threshold: float | None = None,  # 可及性阈值；None=未校准，不参与淘汰
        strict_self_fold: bool = False,     # True 时才把自折叠/稳定茎作为硬淘汰
        temperature: float = 37.0,
        use_cofold: bool = True,
        use_plfold: bool = True,
        backend: str = "auto",              # auto|python|cli|approx（测试用）
        vienna_bin_dir: str | None = None,  # 显式指定 ViennaRNA 预编译 bin 目录
        vienna_python_dir: str | None = None,  # 显式指定含 RNA 模块的目录/环境
    ) -> None:
        self.mfe_threshold = float(mfe_threshold)
        self.end_diff_threshold = float(end_diff_threshold)
        self.stem_threshold = int(stem_threshold)
        self.stem_gc_threshold = float(stem_gc_threshold)
        self.cofold_min = float(cofold_min)     # 保留仅作"杂交偏强"提示阈值
        self.cofold_max = float(cofold_max)     # 语义：杂交过弱阈值（须≤-10）
        self.access_threshold = (float(access_threshold)
                                 if access_threshold is not None else None)
        self.strict_self_fold = bool(strict_self_fold)
        self.temperature = float(temperature)
        self.use_cofold = bool(use_cofold)
        self.use_plfold = bool(use_plfold)
        self.backend = backend                 # 强制后端（测试/降级模拟用）
        self.vienna_bin_dir = vienna_bin_dir   # 显式二进制目录（可选）
        self.vienna_python_dir = vienna_python_dir  # 显式 Python 绑定目录（可选）
        # 结果缓存：同一引导链的 RNAfold 结果不重复计算
        self._fold_cache: dict = {}

        # 探测各工具可用性（先做二进制目录自动发现）
        self.tool_status = self._probe_tools()
        self._warned_approx = False

    # ------------------------------------------------------------------
    def _discover_bin_dirs(self) -> list:
        """
        自动发现 ViennaRNA 预编译二进制目录并加入 os.environ['PATH']。
        候选来源（按优先级）：
          1) __init__ 显式传入 vienna_bin_dir
          2) 环境变量 VIENNARNA_BIN
          3) 模块同目录下常见命名目录（或其 bin/ 子目录）：
             viennarna_windows / viennarna / viennarna_bin / ViennaRNA-win64 ...
          4) 系统级常见安装位置（Windows 官方安装包默认落点 / Unix 常见前缀）
          5) Windows 下兜底扫描 Program Files 与 Program Files (x86) 中
             名字含 "Vienna" 的一级目录（安装向导若改了目录名也能命中）
        只要目录内存在 RNAfold.exe（或 RNAfold）即视为有效。
        返回实际并入 PATH 的目录列表。
        """
        found_dirs = []
        module_dir = os.path.dirname(os.path.abspath(__file__))
        candidates = []
        if self.vienna_bin_dir:
            candidates.append(self.vienna_bin_dir)
        env_bin = os.environ.get("VIENNARNA_BIN")
        if env_bin:
            candidates.append(env_bin)
        for name in ("viennarna_windows", "viennarna", "viennarna_bin",
                     "ViennaRNA-win64", "viennarna-2.7.2", "ViennaRNA-2.7.2"):
            candidates.append(os.path.join(module_dir, name))
        # 再扫描模块目录下名字含 viennarna 的一级子目录
        try:
            for child in os.listdir(module_dir):
                if "viennarna" in child.lower():
                    candidates.append(os.path.join(module_dir, child))
        except OSError:
            pass
        # 系统级常见安装位置（ViennaRNA 官方 Windows 安装包默认装到
        # "C:\Program Files (x86)\ViennaRNA Package"，可执行文件直接放根目录，
        # 不一定是 bin/；注意该目录含空格，subprocess 传列表参数不受影响）
        candidates.extend((
            r"C:\Program Files\ViennaRNA Package",
            r"C:\Program Files (x86)\ViennaRNA Package",
            r"C:\Program Files\ViennaRNA",
            r"C:\Program Files (x86)\ViennaRNA",
            "/usr/local/bin", "/usr/bin", "/opt/homebrew/bin",
            "/usr/local/opt/viennarna/bin",
        ))
        # 兜底：扫 Program Files 两个根下名字含 "Vienna" 的一级目录
        for base in (os.environ.get("ProgramFiles"),
                     os.environ.get("ProgramFiles(x86)"),
                     os.environ.get("ProgramW6432")):
            if not base or not os.path.isdir(base):
                continue
            try:
                for child in os.listdir(base):
                    if "vienna" in child.lower():
                        candidates.append(os.path.join(base, child))
            except OSError:
                pass

        seen = set()
        for cand in candidates:
            if not cand or cand in seen:
                continue
            seen.add(cand)
            for root in (cand, os.path.join(cand, "bin")):
                if not os.path.isdir(root):
                    continue
                has_fold = (os.path.exists(os.path.join(root, "RNAfold.exe"))
                            or os.path.exists(os.path.join(root, "RNAfold")))
                if has_fold and root not in found_dirs:
                    found_dirs.append(root)
        if found_dirs:
            extra = os.pathsep.join(found_dirs)
            cur = os.environ.get("PATH", "")
            os.environ["PATH"] = extra + (os.pathsep + cur if cur else "")
        return found_dirs

    # ------------------------------------------------------------------
    def _probe_tools(self) -> dict:
        """探测 Python 绑定与三个命令行工具；PATH 或二进制目录变化时重探。"""
        # 先做二进制目录发现（会把找到的目录并入 PATH）
        self._found_bin_dirs = self._discover_bin_dirs()
        path_now = os.environ.get("PATH", "")
        bins_key = "|".join(self._found_bin_dirs)
        # PATH 或发现的二进制目录与上次不同 → 清除缓存重新探测
        if (path_now != StructureDetector._cached_path
                or bins_key != StructureDetector._cached_bins):
            StructureDetector._cli_cache.clear()
            StructureDetector._binding_ok = None
            StructureDetector._cached_path = path_now
            StructureDetector._cached_bins = bins_key

        # Python 绑定：显式目录优先（conda 环境等），否则常规 import
        if StructureDetector._binding_ok is None:
            ok = False
            try:
                if self.vienna_python_dir:
                    sys.path.insert(0, self.vienna_python_dir)
                import RNA  # type: ignore   # ViennaRNA Python 绑定
                RNA.version()
                ok = True
            except Exception:
                ok = False
            StructureDetector._binding_ok = ok
        st = {"python_binding": StructureDetector._binding_ok}
        for tool in ("RNAfold", "RNAcofold", "RNAplfold"):
            if tool not in StructureDetector._cli_cache:
                # shutil.which 快速；找不到再试 subprocess（可能 PATH 未刷新）
                found = shutil.which(tool)
                if found is None:
                    try:
                        p = subprocess.run([tool, "--version"],
                                           capture_output=True, timeout=10)
                        found = p.returncode == 0
                    except Exception:
                        found = False
                StructureDetector._cli_cache[tool] = found
            st[tool.lower() + "_cli"] = bool(StructureDetector._cli_cache[tool])
        st["bin_dirs"] = self._found_bin_dirs
        # 记录解析出的绝对路径：Windows 下 RNAplfold 等按 PATH 名称调用不稳定，
        # 统一用绝对路径执行（shutil.which 已随 PATH 注入生效）
        self._tool_paths = {}
        for tool in ("RNAfold", "RNAcofold", "RNAplfold"):
            if st[tool.lower() + "_cli"]:
                self._tool_paths[tool] = shutil.which(tool)
            else:
                self._tool_paths[tool] = None
        # 综合层级：0 完整绑定 / 1 命令行全齐 / 2 仅 RNAfold / 3 无
        bind, f, c, p = (st["python_binding"], st["rnafold_cli"],
                         st["rnacofold_cli"], st["rnaplfold_cli"])
        if bind and f and c and p:
            st["level"] = 0
        elif f and c and p:
            st["level"] = 1
        elif f or bind:
            st["level"] = 2
        else:
            st["level"] = 3
        return st

    # ------------------------------------------------------------------
    def get_tool_status(self) -> dict:
        """返回当前各工具可用状态与降级层级。"""
        return dict(self.tool_status)

    # ------------------------------------------------------------------
    # 输入校验（中文报错）
    # ------------------------------------------------------------------
    def _validate_input(self, sirna: str, mrna: str | None = None) -> tuple:
        """校验序列；返回规范化后的 (sirna, mrna)。"""
        if not isinstance(sirna, str):
            raise TypeError(f"siRNA 应为字符串，实际为 {type(sirna).__name__}")
        s = normalize(sirna)
        if len(s) != 19:
            raise ValueError(
                f"siRNA 引导链长度必须为 19 nt（当前 {len(s)} nt）")
        bad = sorted({c for c in s if c not in VALID_BASES})
        if bad:
            raise ValueError(f"siRNA 含非法碱基：{bad}；仅允许 A/C/G/U/T")
        m = None
        if mrna is not None:
            if not isinstance(mrna, str):
                raise TypeError(f"mRNA 应为字符串，实际为 {type(mrna).__name__}")
            m = normalize(mrna)
            if len(m) != 57:
                raise ValueError(
                    f"mRNA 长度必须为 57 nt（5'侧翼19+靶标19+3'侧翼19；"
                    f"当前 {len(m)} nt）")
            bad = sorted({c for c in m if c not in VALID_BASES})
            if bad:
                raise ValueError(f"mRNA 含非法碱基：{bad}；仅允许 A/C/G/U/T")
        return s, m

    # ==================================================================
    # 核心方法
    # ==================================================================
    def detect(self, sirna: str, mrna: str | None = None) -> dict:
        """单条结构检测。mrna 为 None 时只做引导链相关检测（单链模式）。"""
        s, m = self._validate_input(sirna, mrna)
        out = {"sirna": s, "mrna": m}

        # ---- Step 1: 引导链自折叠（RNAfold）----
        fold = self._call_rnafold(s)
        out.update(fold)

        # ---- Step 2: 末端不对称性（NN 表，不依赖 ViennaRNA）----
        dg5, dg3, ddg = self._calculate_end_diff(s)
        out["delta_G_5end"] = dg5
        out["delta_G_3end"] = dg3
        out["delta_deltaG_ends"] = ddg

        # ---- Step 3: 杂交（RNAcofold）----
        cofold = {"cofold_mfe": None, "cofold_structure": None,
                  "cofold_5end_bp_prob": None, "cofold_3end_bp_prob": None,
                  "internal_loops": None, "bulges": None,
                  "binding_energy": None, "cofold_available": False}
        if m is not None and self.use_cofold:
            target = m[19:38]                      # 靶标区 19 nt
            cofold = self._call_rnacofold(s, target,
                                          guide_self_mfe=fold.get("mfe"))
        out.update(cofold)

        # ---- Step 4: 靶标可及性（RNAplfold）----
        plfold = {"target_accessibility": None, "target_5end_access": None,
                  "target_3end_access": None, "plfold_available": False}
        if m is not None and self.use_plfold:
            plfold = self._call_rnaplfold(m)
        out.update(plfold)

        # ---- Step 5: 硬阈值过滤 ----
        filt = self._apply_filters(out)
        out.update(filt)
        return out

    # ------------------------------------------------------------------
    def detect_batch(self, sirna_list: list, mrna_list: list | None = None) -> list:
        """批量检测。mrna_list 为 None → 单链模式（跳过杂交与可及性）。"""
        if not isinstance(sirna_list, list):
            raise TypeError("sirna_list 应为 list")
        if mrna_list is not None:
            if not isinstance(mrna_list, list):
                raise TypeError("mrna_list 应为 list")
            if len(sirna_list) != len(mrna_list):
                raise ValueError(
                    f"siRNA 与 mRNA 数量不一致：{len(sirna_list)} vs "
                    f"{len(mrna_list)}")
        out = []
        for i, s in enumerate(sirna_list):
            m = mrna_list[i] if mrna_list is not None else None
            out.append(self.detect(s, m))
        return out

    # ------------------------------------------------------------------
    def filter_candidates(self, sirna_list: list,
                          mrna_list: list | None = None) -> list:
        """硬阈值过滤：返回 pass_all=True 的候选（含完整结构字典）。"""
        results = self.detect_batch(sirna_list, mrna_list)
        return [r for r in results if r.get("pass_all")]

    # ==================================================================
    # 工具调用层（分层降级）
    # ==================================================================
    def _vienna_python_ok(self) -> bool:
        return (self.backend in ("auto", "python") and
                self.tool_status["python_binding"])

    def _cli_ok(self, tool: str) -> bool:
        if self.backend == "approx":
            return False
        if self.backend == "cli":
            return bool(self.tool_status[tool.lower() + "_cli"])
        return bool(self.tool_status[tool.lower() + "_cli"])

    # ------------------------------------------------------------------
    # Step 1：RNAfold 引导链自折叠
    # ------------------------------------------------------------------
    def _call_rnafold(self, seq: str) -> dict:
        """引导链 RNAfold：绑定 > 命令行 > 纯 Python NN 近似。"""
        if seq in self._fold_cache:
            return dict(self._fold_cache[seq])

        if self._vienna_python_ok():
            try:
                res = self._rnafold_binding(seq)
                self._fold_cache[seq] = res
                return dict(res)
            except Exception as e:                      # 绑定失败 → 降级
                print(f"[降级] RNAfold 绑定失败({e})，转命令行")
        if self._cli_ok("RNAfold"):
            try:
                res = self._rnafold_cli(seq)
                self._fold_cache[seq] = res
                return dict(res)
            except Exception as e:
                print(f"[降级] RNAfold 命令行失败({e})，转纯 Python 近似")
        # Level 3：纯 Python NN 近似（仅估计引导链发夹；标记 approx）
        if not self._warned_approx:
            print("[警告] ViennaRNA 不可用 → 引导链 MFE 采用纯 Python 最近邻近似；"
                  "cofold/plfold 输出 NaN（建议安装：conda install -c bioconda "
                  "viennarna）")
            self._warned_approx = True
        res = self._rnafold_approx_python(seq)
        self._fold_cache[seq] = res
        return dict(res)

    def _rnafold_binding(self, seq: str) -> dict:
        """Level 0：ViennaRNA Python 绑定。"""
        import RNA  # type: ignore
        fc = RNA.fold_compound(seq)
        mfe_ss, mfe = fc.mfe()
        ens = fc.pf()                          # 集合自由能 (kcal/mol)
        try:
            centroid_ss, centroid_dist = fc.centroid()
        except Exception:
            centroid_ss, centroid_dist = None, None
        try:
            bpp = fc.bpp()                     # (L,L) 配对概率矩阵
            max_bp = max(max(row) for row in bpp) if bpp else None
        except Exception:
            max_bp = None
        parsed = self._parse_dot_bracket(mfe_ss, seq)
        res = {
            "mfe": round(float(mfe), 3), "mfe_structure": mfe_ss,
            "ensemble_dg": round(float(ens), 3),
            "centroid_structure": centroid_ss,
            "centroid_dist": (round(float(centroid_dist), 3)
                              if centroid_dist is not None else None),
            "max_bp_prob": (round(float(max_bp), 4) if max_bp is not None
                            else None),
            "stem_loop_count": parsed["stem_loop_count"],
            "max_stem_length": parsed["max_stem_length"],
            "max_stem_gc": parsed["max_stem_gc"],
            "fold_backend": "python",
        }
        return res

    def _rnafold_cli(self, seq: str) -> dict:
        """Level 1/2：RNAfold --noPS -p -T 37 命令行。"""
        exe = self._tool_paths.get("RNAfold") or "RNAfold"
        cmd = [exe, "--noPS", "-p", "-T", str(self.temperature)]
        p = subprocess.run(cmd, input=seq + "\n", capture_output=True,
                           text=True, timeout=10)
        if p.returncode != 0:
            raise RuntimeError(f"RNAfold 退出码 {p.returncode}: {p.stderr}")
        lines = [ln.strip() for ln in p.stdout.splitlines() if ln.strip()]
        # RNAfold 输出中第一行 '(((...))) ( -4.50)' 即 MFE 结构
        # （-p 开启配分函数后其后还有 pf 行与可能的质心行，勿取末行）
        mfe = None
        mfe_ss = None
        ensemble_dg = None
        for ln in lines:
            mm = re.search(r"^([.()]+)\s*\(\s*([-\d.]+)\s*\)", ln)
            if mm:
                mfe_ss, mfe = mm.group(1), float(mm.group(2))
                break
        m_ens = re.search(r"\[\s*([-\d.]+)\s*\]", p.stdout)
        if m_ens:
            ensemble_dg = float(m_ens.group(1))
        if mfe is None:
            raise RuntimeError("RNAfold 输出无法解析")
        parsed = self._parse_dot_bracket(mfe_ss, seq)
        return {
            "mfe": round(mfe, 3), "mfe_structure": mfe_ss,
            "ensemble_dg": (round(ensemble_dg, 3)
                            if ensemble_dg is not None else None),
            "centroid_structure": None, "centroid_dist": None,
            "max_bp_prob": None,
            "stem_loop_count": parsed["stem_loop_count"],
            "max_stem_length": parsed["max_stem_length"],
            "max_stem_gc": parsed["max_stem_gc"],
            "fold_backend": "cli",
        }

    def _rnafold_approx_python(self, seq: str) -> dict:
        """
        Level 3：纯 Python 最近邻近似自折叠。
        只搜索"发夹"结构（siRNA 自折叠的最常见强结构）：扫描所有可能的
        茎-环组合，茎以 WC/G:U 配对，能量 = Σ 茎内 NN 堆叠 + 环惩罚；
        取最负能量作为 MFE 近似值并构造对应点括号，供茎统计复用。
        结果标 approx=True，仅供无 ViennaRNA 时粗筛，不替代真实 RNAfold。
        """
        n = len(seq)
        best_e = 0.0
        best_pairs = []                      # [(左臂guide位, 右臂guide位), ...]
        # 遍历外端对 (i, j)，j-i>=4 保证环长≥3
        for i in range(n - 4):
            for j in range(i + 4, n):
                pairs = []
                k = 0
                while (i + k < j - k - 3 and
                       is_comp(seq[i + k], seq[j - k])):
                    pairs.append((i + k, j - k))
                    k += 1
                if len(pairs) < 1:
                    continue
                # 茎内堆叠能（沿 5'->3' guide 取二核苷酸）
                stem_energy = 0.0
                # 相邻两对之间的堆叠：左臂位 i+k, i+k+1
                for k in range(len(pairs) - 1):
                    stem_energy += STACK_DG.get(seq[i + k] + seq[i + k + 1], 0.0)
                loop_len = (pairs[-1][1] - pairs[-1][0]) - 1
                loop_pen = 0.9 + 0.2 * max(0, loop_len - 3)   # 近似环惩罚
                hairpin_e = stem_energy + loop_pen
                if hairpin_e < best_e:
                    best_e = hairpin_e
                    best_pairs = list(pairs)
        if best_e < 0:
            bracket = [_UNPAIRED] * n
            for (a, b) in best_pairs:
                bracket[a] = _PAIR_OPEN
                bracket[b] = _PAIR_CLOSE
            mfe_ss = "".join(bracket)
        else:
            mfe_ss = _UNPAIRED * n
        parsed = self._parse_dot_bracket(mfe_ss, seq)
        return {
            "mfe": round(best_e, 3),
            "mfe_structure": mfe_ss,
            "ensemble_dg": None, "centroid_structure": None,
            "centroid_dist": None, "max_bp_prob": None,
            "stem_loop_count": parsed["stem_loop_count"],
            "max_stem_length": parsed["max_stem_length"],
            "max_stem_gc": parsed["max_stem_gc"],
            "fold_backend": "approx",
        }

    # ------------------------------------------------------------------
    # Step 3：RNAcofold 杂交
    # ------------------------------------------------------------------
    def _call_rnacofold(self, guide: str, target: str,
                        guide_self_mfe: float | None = None) -> dict:
        """杂交 RNAcofold：绑定 > 命令行；不可用返回 NaN 字段。"""
        out = {"cofold_mfe": None, "cofold_structure": None,
               "cofold_5end_bp_prob": None, "cofold_3end_bp_prob": None,
               "internal_loops": None, "bulges": None,
               "binding_energy": None, "cofold_available": False}
        if self._vienna_python_ok():
            try:
                import RNA  # type: ignore
                ss, mfe = RNA.cofold(guide + "&" + target)
                return self._finalize_cofold(ss, mfe, guide, target,
                                             guide_self_mfe, "python")
            except Exception as e:
                print(f"[降级] RNAcofold 绑定失败({e})，转命令行")
        if self._cli_ok("RNAcofold"):
            try:
                inp = f">{guide}_{target}\n{guide}&{target}\n"
                exe = self._tool_paths.get("RNAcofold") or "RNAcofold"
                p = subprocess.run([exe, "--noPS", "-T", str(self.temperature)],
                                   input=inp, capture_output=True,
                                   text=True, timeout=10)
                if p.returncode != 0:
                    raise RuntimeError(p.stderr)
                lines = [ln.strip() for ln in p.stdout.splitlines() if ln.strip()]
                ss, mfe = None, None
                for ln in reversed(lines):
                    mm = re.search(r"^([.()&]+)\s*\(\s*([-\d.]+)\s*\)", ln)
                    if mm and "&" in mm.group(1):
                        ss, mfe = mm.group(1), float(mm.group(2))
                        break
                if ss is None:
                    raise RuntimeError("RNAcofold 输出无法解析")
                return self._finalize_cofold(ss, mfe, guide, target,
                                             guide_self_mfe, "cli")
            except Exception as e:
                print(f"[降级] RNAcofold 命令行失败({e})")
        # 不可用：NaN 字段 + 提示
        out["cofold_note"] = "cofold 不可用（ViennaRNA 未安装或被禁用）"
        return out

    def _finalize_cofold(self, ss: str, mfe: float, guide: str, target: str,
                         guide_self_mfe: float | None, backend: str) -> dict:
        """RNAcofold 结构解析与字段组装。"""
        parsed = self._parse_cofold_structure(ss)
        # 纯结合能 ≈ cofold MFE − 引导链自折叠 MFE − 靶窗自折叠 MFE
        tgt_self = None
        tgt_fold = self._call_rnafold(target)        # 靶窗 19 nt 自折叠
        tgt_self = tgt_fold.get("mfe")
        binding = None
        if mfe is not None and guide_self_mfe is not None and tgt_self is not None:
            binding = mfe - guide_self_mfe - tgt_self
        return {
            "cofold_mfe": (round(float(mfe), 3) if mfe is not None else None),
            "cofold_structure": ss,
            "cofold_5end_bp_prob": parsed["end5_prob"],
            "cofold_3end_bp_prob": parsed["end3_prob"],
            "internal_loops": parsed["internal_loops"],
            "bulges": parsed["bulges"],
            "binding_energy": (round(binding, 3) if binding is not None
                               else None),
            "cofold_available": True,
            "cofold_backend": backend,
        }

    # ------------------------------------------------------------------
    # Step 4：RNAplfold 靶标可及性
    # ------------------------------------------------------------------
    def _call_rnaplfold(self, mrna: str) -> dict:
        """
        RNAplfold：-W 40 -L 40 -u 19（2.7.x 默认输出 plfold_lunp）。
        兼容性说明：Windows 沙箱/中文路径下直接 subprocess 调用 RNAplfold 不稳定
        （MinGW 版本偶发访问违例/挂起），故通过 cmd /c + 文件重定向在 ASCII 系统
        临时目录运行，python 仅读写文件——实测 4/4 稳定。
        """
        out = {"target_accessibility": None, "target_5end_access": None,
               "target_3end_access": None, "plfold_available": False}
        if not self._cli_ok("RNAplfold"):
            out["plfold_note"] = "plfold 不可用（ViennaRNA 未安装或被禁用）"
            return out
        exe = self._tool_paths.get("RNAplfold") or "RNAplfold"
        # ASCII 工作目录（系统临时目录下唯一子目录），避开中文路径
        import uuid
        workdir = os.path.join(tempfile.gettempdir(),
                               f"vrnapl_{os.getpid()}_{uuid.uuid4().hex[:8]}")
        try:
            os.makedirs(workdir, exist_ok=True)
            seq_f = os.path.join(workdir, "in.seq")
            with open(seq_f, "w", encoding="utf-8") as fh:
                fh.write(mrna + "\n")
            out_f = os.path.join(workdir, "out.txt")
            err_f = os.path.join(workdir, "err.txt")
            # 关键经验：Windows 沙箱下 RNAplfold 的 stdin/stdout 必须用真实文件
            # 句柄（PIPE/DEVNULL 或 cmd/批处理包装均可能触发 MinGW 访问违例或
            # 码页问题），文件句柄 + ASCII cwd 实测稳定。
            with open(seq_f, "rb") as fi, open(out_f, "wb") as fo, \
                    open(err_f, "wb") as fe:
                p = subprocess.run([exe, "-W", "40", "-L", "40", "-u", "19"],
                                   stdin=fi, stdout=fo, stderr=fe,
                                   timeout=30, cwd=workdir)
            err_txt = ""
            if os.path.exists(err_f):
                with open(err_f, encoding="utf-8", errors="replace") as fh:
                    err_txt = fh.read().strip()
            if p.returncode != 0:
                raise RuntimeError(f"RNAplfold 退出码 {p.returncode}: {err_txt}")
            # 2.7.x 默认输出 plfold_lunp / plfold_dp.ps
            lunp = None
            for fn in os.listdir(workdir):
                if fn.endswith("_lunp"):
                    lunp = os.path.join(workdir, fn)
                    break
            if lunp is None:
                raise RuntimeError("RNAplfold 未生成 _lunp 文件")
            mat = self._parse_lunp(lunp)
            return self._accessibility_from_lunp(mrna, mat)
        except Exception as e:
            print(f"[降级] RNAplfold 失败：{e}")
            out["plfold_note"] = f"RNAplfold 执行失败：{e}"
            return out
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    @staticmethod
    def _parse_lunp(path: str) -> dict:
        """
        解析 RNAplfold 的 *_lunp 文件。
        假设格式：每行 = 起始位置(整数) + 各长度(1..W) 的未配对概率。
        返回 {pos(0-based): [p_len1, p_len2, ...]}。
        若行首位置疑似 1-based 则自动平移为 0-based；越界概率做合法性检查。
        """
        rows = []
        with open(path, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                toks = ln.split()
                try:
                    pos = int(toks[0])
                except ValueError:
                    continue                    # 跳过 '#...' 表头行
                # 逐列解析：'NA' 等占位符视为缺失(None)，不整行丢弃
                vals = []
                for x in toks[1:]:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        vals.append(None)
                if vals:
                    rows.append((pos, vals))
        if not rows:
            return {}
        # 判断基准：文件里位置是否出现 0（0-based）还是从 1 开始（1-based）
        min_pos = min(p for p, _ in rows)
        if min_pos == 0:
            base = 0
        else:
            base = 1
        return {p - base: v for p, v in rows}

    def _accessibility_from_lunp(self, mrna: str, mat: dict) -> dict:
        """
        从 lunp 矩阵提取靶标区可及性（0-based 靶标区 [19,38)）。

        口径说明（重要）：RNAplfold 的"长度 L 列"是"区间 [i,i+L-1] 整体同时未配对"
        的概率——对 19 mer 目标段该值几乎恒为 0，无区分度。文献（Heale 2005、
        Sfold 靶位可及性规则）实际用的是**靶位各核苷酸未配对概率的平均值**，
        即"长度 1 列"在靶标区上的均值。故：
          target_accessibility = mean_i∈[19,38) P(碱基i未配对)
          5'/3'端同理取对应 5 nt。
        整窗 19 mer 同时未配对概率另存 target_full_unpaired_19 作参考。
        """
        def _base_up(i: int):
            row = mat.get(i)
            if not row or row[0] is None:
                return None
            v = row[0]
            return v if 0.0 <= v <= 1.0 else None

        def _mean(idx_range):
            vals = [v for i in idx_range if (v := _base_up(i)) is not None]
            return sum(vals) / len(vals) if vals else None

        acc_full = _mean(range(19, 38))
        a5 = _mean(range(19, 24))
        a3 = _mean(range(33, 38))

        # 参考值：整段 19 mer 同时未配对概率（通常≈0，仅存档）
        p19 = None
        row = mat.get(19)
        if row and len(row) >= 19 and row[18] is not None:
            v = row[18]
            if 0.0 <= v <= 1.0:
                p19 = round(float(v), 6)

        out = {"plfold_available": True, "plfold_backend": "cli"}
        out["target_accessibility"] = (round(acc_full, 4)
                                       if acc_full is not None else None)
        out["target_5end_access"] = (round(a5, 4) if a5 is not None else None)
        out["target_3end_access"] = (round(a3, 4) if a3 is not None else None)
        out["target_full_unpaired_19"] = p19
        if acc_full is None:
            out["plfold_note"] = "未配对概率缺失：未能计算靶标可及性"
        return out

    # ==================================================================
    # Step 2：末端不对称性（NN 表；不依赖 ViennaRNA）
    # ==================================================================
    def _calculate_end_diff(self, sirna: str) -> tuple:
        """
        末端不对称（链选择方向，硬淘汰依据；Khvorova 2003 Cell / Schwarz 2003 Cell）。

        说明（修订）：RISC 装载哪条链由 siRNA **双链两端**的相对热力学稳定性决定。
        对 19 nt 完全配对的 siRNA，乘客链(passenger/sense) = rc(引导链)，因此
        "引导链 5' 端前 2 nt 与自身互补配对的堆叠能" 恰好等于该端双链的端部堆叠能
        —— 单链二核苷酸 NN 表在完美配对时即等价于双链端稳定度（本表隐含互补底链）。
        ΔG_5end = ΔG(siRNA[0:2]) + 端部 A·U 惩罚(5' 端碱基 A/U 时 +0.45)
        ΔG_3end = ΔG(siRNA[17:19]) + 端部 A·U 惩罚(3' 端碱基 A/U 时 +0.45)
        delta_deltaG_ends = ΔG_5end − ΔG_3end
        负值过大（< -0.5）→ 5' 端过稳 → RISC 倾向装载错误链 = 链选择反转。
        注：端部 A·U 惩罚同时隐含"g1 为 U/A 时 5' 端更不稳定"的链选择偏好。
        """
        dg5 = STACK_DG.get(sirna[0:2], 0.0)
        dg3 = STACK_DG.get(sirna[17:19], 0.0)
        if sirna[0] in "AU":
            dg5 += 0.45
        if sirna[18] in "AU":
            dg3 += 0.45
        ddg = dg5 - dg3
        return round(dg5, 3), round(dg3, 3), round(ddg, 3)

    # ==================================================================
    # 点括号解析
    # ==================================================================
    def _parse_dot_bracket(self, struct: str, seq: str | None = None) -> dict:
        """
        解析单链点括号：统计发夹(茎环)数、最长茎长(bp)、最长茎 GC%。
        stem = 连续堆叠的配对列块；发夹 = 茎内不再包含其他配对的茎。
        seq：对应序列（计算最长茎 GC 含量用，可为 None）。
        """
        empty = {"stem_loop_count": 0, "max_stem_length": 0,
                 "max_stem_gc": None}
        if not struct or set(struct) - set(".()"):
            return empty
        # 1) 括号配对：close -> open
        stack = []
        pairs = {}
        for i, ch in enumerate(struct):
            if ch == _PAIR_OPEN:
                stack.append(i)
            elif ch == _PAIR_CLOSE:
                if stack:
                    pairs[i] = stack.pop()
        if not pairs:
            return empty
        # 2) 堆叠块（连续配对列：close 递增 & open 递减）
        blocks = []                       # 每块 = [close索引, ...]（递增）
        used = set()
        for c in sorted(pairs):
            if c in used:
                continue
            block = [c]
            k = c
            while (k + 1) in pairs and pairs[k + 1] == pairs[k] - 1:
                block.append(k + 1)
                k += 1
            k = c
            while (k - 1) in pairs and pairs[k - 1] == pairs[k] + 1:
                block.append(k - 1)
                k -= 1
            block.sort()
            for x in block:
                used.add(x)
            blocks.append(block)
        # 3) 发夹数：某块最内层对（close 最小）内部无配对
        hairpins = 0
        for block in blocks:
            inner_c = min(block)
            inner_o = pairs[inner_c]
            interior = struct[inner_o + 1:inner_c]
            if "(" not in interior and ")" not in interior:
                hairpins += 1
        # 4) 最长茎及 GC 含量
        max_len = max(len(b) for b in blocks)
        longest = max(blocks, key=len)
        stem_gc = None
        if seq is not None:
            open_pos = [pairs[c] for c in longest]       # '(' 侧位置
            close_pos = [c for c in longest]
            idxs = open_pos + close_pos
            idxs = [i for i in idxs if 0 <= i < len(seq)]
            if idxs:
                gc = sum(1 for i in idxs if seq[i] in "GC")
                stem_gc = round(100.0 * gc / len(idxs), 1)
        return {"stem_loop_count": hairpins,
                "max_stem_length": max_len,
                "max_stem_gc": stem_gc}

    def _parse_cofold_structure(self, struct: str) -> dict:
        """
        解析 RNAcofold 杂交结构（含 &）。
        假设引导链与靶标区在杂交区内按设计反平行连续配对：
            列 c (0..18)：guide[c] 对 target[18-c]
        按列判定配对/单链，统计：
            internal_loops : 连续多列两侧均未配对的区段数
            bulges         : 连续多列仅一侧未配对的区段数
            end5/end3 prob : 引导链 5'/3' 端各 2 nt 在 MFE 结构中的配对比例
        （若结构与假设不符——如无杂交——返回 NaN 并附提示。）
        """
        empty = {"internal_loops": None, "bulges": None,
                 "end5_prob": None, "end3_prob": None,
                 "parse_note": "结构不可解析（无 & 或长度不符）"}
        if not struct or "&" not in struct:
            return empty
        parts = struct.split("&")
        if len(parts) != 2:
            return empty
        s1, s2 = parts[0], parts[1]
        # 允许目标串与 19 nt 不一致（RNAcofold 可能给出端部 dangling）
        n_col = min(len(s1), 19, 19)
        if n_col < 1 or len(s2) < 1:
            return empty
        # 列配对状态
        col_pair = []                 # True=该列两侧均配对
        for c in range(n_col):
            t = (18 - c) if c < 19 and (18 - c) < len(s2) else -1
            if t < 0:
                col_pair.append(False)
                continue
            g_dot = s1[c] == _UNPAIRED
            t_dot = s2[t] == _UNPAIRED
            col_pair.append((not g_dot) and (not t_dot))
        # 连续区段统计
        internal = bulge = 0
        prev_state = None
        for c in range(n_col):
            g_dot = s1[c] == _UNPAIRED
            t_dot = (s2[18 - c] == _UNPAIRED) if (18 - c) < len(s2) else True
            if g_dot and t_dot:
                state = "I"                       # 内部环（两侧均单链）
            elif g_dot != t_dot:
                state = "B"                       # 凸起（仅一侧单链）
            else:
                state = "P"
            if state != prev_state:
                if state == "I":
                    internal += 1
                elif state == "B":
                    bulge += 1
                prev_state = state
        # 两端配对比例（引导链 5' 端 c=0,1；3' 端 c=17,18）
        def _end_prob(idxs):
            vals = []
            for c in idxs:
                if c < len(s1) and c < 19:
                    vals.append(1.0 if s1[c] in "()" else 0.0)
            return (sum(vals) / len(vals)) if vals else None

        e5 = _end_prob([0, 1])
        e3 = _end_prob([17, 18])
        if internal == 0 and bulge == 0 and e5 is None and e3 is None:
            return empty
        return {"internal_loops": internal, "bulges": bulge,
                "end5_prob": e5, "end3_prob": e3}

    # ==================================================================
    # 硬阈值过滤
    # ==================================================================
    def _apply_filters(self, d: dict) -> dict:
        """
        硬阈值过滤（修订版）。
        硬淘汰只保留"大概率必失败"项：链选择反转(R3)、杂交过弱(R4)、靶标不可及(R5)。
        自折叠 MFE / 稳定茎(R1/R2) 默认降为"软提示"写入 warnings，仅 strict_self_fold
        且结构可靠（真实 ViennaRNA）时才作为硬淘汰。
        工具缺失或阈值未校准的规则置 None（不据此淘汰）。
        文献依据：R1/R2 见 Singh 2012；R3 见 Khvorova 2003 / Schwarz 2003；
        R4 过弱见 Sfold 规则（结合能 ≤ -10）；R5 见 Heale 2005 / Ding 2003 / Bernhart 2006。
        """
        reasons = []       # 硬淘汰原因
        warnings = []      # 软提示（不淘汰）
        # 结构可靠性：仅真实 ViennaRNA（python/cli）产出才可靠；approx 只作提示
        reliable = d.get("fold_backend") in ("python", "cli")

        # a) MFE（自折叠，软提示；strict 且可靠时才硬淘汰）
        mfe = d.get("mfe")
        if mfe is None:
            pass_mfe = None
        else:
            pass_mfe = mfe >= self.mfe_threshold
            if not pass_mfe:
                msg = (f"自折叠偏强(MFE={mfe:.2f}<-{abs(self.mfe_threshold):.1f})")
                if self.strict_self_fold and reliable:
                    reasons.append(msg)
                else:
                    warnings.append("自折叠提示：" + msg)

        # b) 末端不对称（链选择反转，恒为硬淘汰；NN 表不依赖 ViennaRNA）
        ddg = d.get("delta_deltaG_ends")
        pass_end_diff = ddg >= self.end_diff_threshold
        if not pass_end_diff:
            reasons.append(f"链选择反转(5'端过稳ΔΔG={ddg:.2f}<-{abs(self.end_diff_threshold):.1f})")

        # c) 稳定茎（软提示；strict 且可靠时才硬淘汰）
        stem_len = d.get("max_stem_length")
        stem_gc = d.get("max_stem_gc")
        if stem_len is None or stem_gc is None:
            pass_stem = None
        else:
            bad_stem = (stem_len >= self.stem_threshold
                        and stem_gc >= self.stem_gc_threshold)
            pass_stem = not bad_stem
            if bad_stem:
                msg = f"稳定茎(最长{stem_len}bp,GC{stem_gc:.0f}%)"
                if self.strict_self_fold and reliable:
                    reasons.append("强发夹：" + msg)
                else:
                    warnings.append("强发夹提示：" + msg)

        # d) 杂交强度（仅"过弱"淘汰；"过强"只提示，交任务 9 seed 脱靶复核）
        cm = d.get("cofold_mfe")
        if cm is None or not d.get("cofold_available"):
            pass_cofold = None
        else:
            pass_cofold = cm <= self.cofold_max       # 须 ≤ -10（结合足够强）
            if not pass_cofold:
                reasons.append(f"杂交过弱(结合不足 mfe={cm:.2f}>-{abs(self.cofold_max):.1f})")
            elif cm < self.cofold_min:
                warnings.append(f"杂交偏强(mfe={cm:.2f}<{self.cofold_min:.1f}，"
                                f"请交任务9 seed区脱靶复核)")

        # e) 靶标可及性（阈值未校准=None 时不参与淘汰；仅 plfold 可用时判定）
        acc = d.get("target_accessibility")
        if (acc is None or not d.get("plfold_available")
                or self.access_threshold is None):
            pass_accessibility = None
        else:
            pass_accessibility = acc >= self.access_threshold
            if not pass_accessibility:
                reasons.append(f"靶标不可及(可及性{acc:.3f}<{self.access_threshold})")

        # 汇总：硬判定集合 = 端不对称 + 杂交过弱 + 可及性（可选 + 严格自折叠）
        hard = [pass_end_diff, pass_cofold, pass_accessibility]
        if self.strict_self_fold and reliable:
            hard += [pass_mfe, pass_stem]
        hard = [f for f in hard if f is not None]
        pass_all = bool(hard) and all(hard)
        return {
            "pass_mfe": pass_mfe,
            "pass_end_diff": pass_end_diff,
            "pass_stem": pass_stem,
            "pass_cofold": pass_cofold,
            "pass_accessibility": pass_accessibility,
            "pass_all": pass_all,
            "reject_reason": "；".join(reasons),
            "warnings": "；".join(warnings),
            "structure_reliable": reliable,
        }




# ==============================================================================
# 三、内置测试（python structure_detector.py）
# ==============================================================================

def _fmt(v, digits: int = 4) -> str:
    if v is None:
        return "--"
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


def _pass_summary(d: dict) -> str:
    if d["pass_all"]:
        return "通过 ✓"
    return f"淘汰 ✗ 原因: {d['reject_reason']}"


def _run_tests() -> None:
    print("=" * 80)
    print("siRNA 结构检测模块（任务 8）—— 内置测试")
    print("=" * 80)
    det = StructureDetector()
    status = det.get_tool_status()
    print(f"\n[工具状态] level={status['level']}  "
          f"python绑定={status['python_binding']}  "
          f"RNAfold={status['rnafold_cli']}  RNAcofold={status['rnacofold_cli']}  "
          f"RNAplfold={status['rnaplfold_cli']}")
    if status["level"] == 3:
        print("  → ViennaRNA 不可用：cofold/plfold 相关用例标注[跳过]，"
              "RNAfold 用纯 Python NN 近似（approx）。")
    # ---------- 合成 mRNA 57 nt（靶标区随机、两端给侧翼）----------
    def make_mrna(window19: str) -> str:
        f5, f3 = ("ACGU" * 5)[:19], ("UGCA" * 5)[:19]
        return f5 + window19 + f3

    # ====== 用例 1：正常高效 siRNA（无强发夹、5' 端相对不稳定）======
    print("\n" + "-" * 80)
    print("【用例 1】正常候选（5'端AU丰富）：UAUUCCGACAAUAGUACGA")
    print("-" * 80)
    s1 = "UAUUCCGACAAUAGUACGA"
    m1 = make_mrna(rc(s1))          # 靶标 = 反向互补 → 完全配对
    r1 = det.detect(s1, m1)
    for k in ("mfe", "mfe_structure", "delta_G_5end", "delta_G_3end",
              "delta_deltaG_ends", "max_stem_length", "stem_loop_count",
              "cofold_mfe", "target_accessibility"):
        print(f"  {k:<22}: {_fmt(r1.get(k))}")
    print(f"  → {_pass_summary(r1)}")

    # ====== 用例 2：强发夹候选（高 GC 回文）======
    print("\n" + "-" * 80)
    print("【用例 2】强发夹候选（高GC回文）：GCGCGCGCGCGCGCGCGCG")
    print("-" * 80)
    s2 = "GCGCGCGCGCGCGCGCGCG"
    r2 = det.detect(s2)
    for k in ("mfe", "max_stem_length", "max_stem_gc", "pass_mfe",
              "pass_stem", "pass_all"):
        print(f"  {k:<22}: {_fmt(r2.get(k))}")
    print(f"  → {_pass_summary(r2)}")
    print(f"  warnings(软提示): {r2.get('warnings')}")
    # 严格模式：自折叠/稳定茎才作为硬淘汰
    det_strict = StructureDetector(strict_self_fold=True, backend="approx")
    r2s = det_strict.detect(s2)
    print(f"  [严格模式] pass_all={r2s['pass_all']}  淘汰原因={r2s['reject_reason']}")

    # ====== 用例 3：链选择反转候选（5'端 GC 富集、3'端 AU 富集）======
    print("\n" + "-" * 80)
    print("【用例 3】链选择反转候选：GCCGCCGCCAUUAUUAUUA")
    print("-" * 80)
    s3 = "GCCGCCGCCAUUAUUAUUA"
    r3 = det.detect(s3)
    for k in ("delta_G_5end", "delta_G_3end", "delta_deltaG_ends",
              "pass_end_diff", "pass_all"):
        print(f"  {k:<22}: {_fmt(r3.get(k))}")
    print(f"  → {_pass_summary(r3)}")

    # ====== 用例 4：RNAcofold 杂交 ======
    print("\n" + "-" * 80)
    print("【用例 4】RNAcofold 杂交（完全互补 guide + target）")
    print("-" * 80)
    if det.tool_status["level"] >= 3:
        print("  [跳过] ViennaRNA 不可用，cofold 未执行（cofold_mfe=--）")
        # 用合成结构自检 cofold 解析器
        synth = ("(((((((((((((((((((&)))))))))))))))))))")
        parsed = det._parse_cofold_structure(synth)
        print(f"  [解析器自检] 合成完美杂交结构 → internal_loops="
              f"{parsed['internal_loops']}, bulges={parsed['bulges']}, "
              f"end5_prob={parsed['end5_prob']}, end3_prob={parsed['end3_prob']}")
        synth_mm = ("(((((((((.(((((((((&))))))))).)))))))))")  # 双侧同列单链=内环
        parsed2 = det._parse_cofold_structure(synth_mm)
        print(f"  [解析器自检] 含错配杂交结构 → internal_loops="
              f"{parsed2['internal_loops']}, bulges={parsed2['bulges']}")
    else:
        tgt = rc(s1)
        r4 = det.detect(s1, make_mrna(tgt))
        for k in ("cofold_mfe", "cofold_structure", "internal_loops",
                  "bulges", "binding_energy"):
            print(f"  {k:<22}: {_fmt(r4.get(k))}")

    # ====== 用例 5：RNAplfold 靶标可及性 ======
    print("\n" + "-" * 80)
    print("【用例 5】RNAplfold 可及性（高GC靶标 vs 高AU靶标）")
    print("-" * 80)
    if det.tool_status["level"] >= 3:
        print("  [跳过] ViennaRNA 不可用，RNAplfold 未执行")
        # 合成 lunp 解析自检
        import io
        fake = os.path.join(tempfile.gettempdir(), "fake_pl_lunp.txt")
        with open(fake, "w", encoding="utf-8") as f:
            for i in range(1, 58):
                f.write(f"{i} " + " ".join(["0.9"] * 19) + "\n")
        mat = det._parse_lunp(fake)
        acc = det._accessibility_from_lunp("A" * 57, mat)
        print(f"  [解析器自检] 合成 lunp(全0.9) → target_accessibility="
              f"{acc['target_accessibility']}, 5'端={acc['target_5end_access']}, "
              f"3'端={acc['target_3end_access']}")
        os.remove(fake)
    else:
        gc_win = "GCCGCGCCGCGCGCCGCCG"
        au_win = "AUAUAAUAAUAAUAAUAAU"
        rgc = det.detect(s1, make_mrna(gc_win))
        rau = det.detect(s1, make_mrna(au_win))
        print(f"  高GC靶标 accessibility={_fmt(rgc['target_accessibility'])}")
        print(f"  高AU靶标 accessibility={_fmt(rau['target_accessibility'])}")

    # ====== 用例 6：边界值测试（过滤逻辑直接验证）======
    print("\n" + "-" * 80)
    print("【用例 6】边界值：MFE=-4.0 与 ΔΔG=-0.5 恰好通过")
    print("-" * 80)
    boundary = {
        "mfe": -4.0, "mfe_structure": "." * 19, "delta_deltaG_ends": -0.5,
        "max_stem_length": 0, "max_stem_gc": None,
        "cofold_available": False, "cofold_mfe": None,
        "plfold_available": False, "target_accessibility": None,
    }
    rb = det._apply_filters(boundary)
    print(f"  pass_mfe={rb['pass_mfe']}  pass_end_diff={rb['pass_end_diff']}  "
          f"pass_all={rb['pass_all']}  reject='{rb['reject_reason']}'  "
          f"warnings='{rb['warnings']}'")

    # ====== 用例 7：输入校验 ======
    print("\n" + "-" * 80)
    print("【用例 7】输入校验（应抛出中文 ValueError）")
    print("-" * 80)
    for bad_s, bad_m in (("ACGUACGUACGUACGUACGU", None),      # 长度20
                         ("ACGUACGUACGUACGUACX", None),       # 非法碱基X
                         ("ACGUACGUACGUACGUACG", "A" * 56)):  # mRNA长度错
        try:
            det.detect(bad_s, bad_m)
            print(f"  ✗ 未拦截：{bad_s!r} / {bad_m!r}")
        except ValueError as e:
            print(f"  ✓ 已拦截：{e}")

    # ====== 用例 8：批量 10 条混合候选 + 淘汰原因分布 =======
    print("\n" + "-" * 80)
    print("【用例 8】批量 10 条混合候选")
    print("-" * 80)
    pool = [s1, s2, s3, "UUAUACGCAGCAUAUACUC", "UACGCUGAUACUAACGUAC",
            "AUAUGCUAGUCGUACGAUA", "GCACUCAUCAUUGUGCUGC",
            "CCCCCCCCCCCCCCCCCCC", "AAAAAAAAAAAAAAAAAAA", "UAACGAUCCAGUAUCGUAC"]
    results = det.detect_batch(pool)
    passed = [r for r in results if r["pass_all"]]
    from collections import Counter as _C
    why = _C()
    for r in results:
        if not r["pass_all"]:
            for w in r["reject_reason"].split("；"):
                why[w.split("(")[0]] += 1
    for i, r in enumerate(results, 1):
        print(f"  #{i:<2} {r['sirna']}  MFE={_fmt(r['mfe'],1):>7}  "
              f"ΔΔG={_fmt(r['delta_deltaG_ends'],2):>7}  "
              f"{'通过' if r['pass_all'] else '淘汰'}")
    print(f"  通过 {len(passed)}/{len(results)}；淘汰原因分布: {dict(why)}")

    # ====== 用例 9：降级模拟 =======
    print("\n" + "-" * 80)
    print("【用例 9】降级测试：approx 后端 + filter_candidates 单链模式")
    print("-" * 80)
    det_a = StructureDetector(backend="approx")
    st = det_a.get_tool_status()
    print(f"  get_tool_status → level={st['level']}")
    ok = det_a.filter_candidates([s1, s2, s3])
    print(f"  filter_candidates(3条, 单链) → 通过 {len(ok)} 条")
    print("  ✓ Level 3 降级路径无异常")

    # 二进制目录自动发现自检：临时放入 dummy exe，验证 PATH 注入与缓存失效
    tmpbin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "_dummy_vrna_bin")
    os.makedirs(tmpbin, exist_ok=True)
    try:
        for t in ("RNAfold", "RNAcofold", "RNAplfold"):
            with open(os.path.join(tmpbin, t + ".exe"), "w") as f:
                f.write("dummy")
        det_b = StructureDetector(vienna_bin_dir=tmpbin)
        stb = det_b.get_tool_status()
        print(f"  [发现自检] vienna_bin_dir 注入后 → level={stb['level']}  "
              f"RNAfold={stb['rnafold_cli']} RNAcofold={stb['rnacofold_cli']} "
              f"RNAplfold={stb['rnaplfold_cli']} bin_dirs={stb['bin_dirs']}")
    finally:
        shutil.rmtree(tmpbin, ignore_errors=True)

    print("\n" + "=" * 80)
    print("测试完成。ViennaRNA 未安装的用例已按规格标注[跳过/降级]。")
    print("=" * 80)


if __name__ == "__main__":
    _run_tests()
