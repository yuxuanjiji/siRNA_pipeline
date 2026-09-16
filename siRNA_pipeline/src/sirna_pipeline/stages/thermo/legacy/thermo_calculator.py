# -*- coding: utf-8 -*-
"""
================================================================================
模块名称：siRNA-mRNA 杂交热力学参数计算模块
================================================================================
功能说明
--------
    针对 19 nt siRNA（反义链/引导链）与 57 nt mRNA
    （结构 = 19 nt 5'侧翼 + 19 nt 靶标区 + 19 nt 3'侧翼）的杂交双链，
    基于核酸热力学的最近邻模型（Nearest-Neighbor, NN Model）计算：
        * 基础杂交自由能 ΔG°(37 ℃)、焓变 ΔH°、熵变 ΔS°（完全互补配对区）
        * mRNA 5'/3' 侧翼序列对靶标区两端的"悬垂堆叠（dangling end）"贡献
        * 错配位点能量修正：
            - G:U 摆动配对（专用摆动参数）
            - 单碱基错配（A:A / A:C / C:C / C:U / G:A / G:G / U:U 等非 WC、
              非 G:U 配对），按碱基类型 + 两侧相邻闭合配对（上下文）修正
            - 连续多碱基错配：按"错配区段（run）"叠加计算
            - 位置带权重（位置权重）：按 guide 5'→3' 位点分 4 个位置带 ——
              g1(5'端,宽容) / seed(g2-8,放大) / central(g9-12,切割敏感区,最强)
              / 3'端(g13-19,基准)；默认权重 0.5 / 2.0 / 2.5 / 1.0，
              各权重系数可配置、可训练（override_params / trainable）
        * 最终输出：杂交 ΔG (kcal/mol)、ΔH (kcal/mol)、ΔS (cal/(mol·K))、
          解链温度 Tm (℃)
    全部计算基于 PyTorch 张量实现：支持批量输入（batch 维）、全程可微分，
    可作为可学习子模块嵌入深度学习模型进行端到端训练/微调。

输入
----
    seq_sirna : 19 nt 反义链（RNA，U/T 均可；内部统一转为 U）
    seq_mrna  : 57 nt mRNA（5'侧翼19 + 靶标区19 + 3'侧翼19）
                靶标区固定为 mRNA 的第 20–38 位（0-based 索引 19:38）
    批量模式：seq_sirna_list / seq_mrna_list；亦支持碱基索引张量输入
    碱基编码：A=0, C=1, G=2, U=3

输出
----
    dict（单条：Python 标量；批量：同形状 torch.Tensor，维度 (B,) 或 (B,19)）
        分项拆解 + 汇总，全部中间能量项可见：
        dG_core_stack / dH_core_stack      最近邻堆叠项（理想互补参考双链）
        dG_init_term / dH_init_term        双链起始项
        dG_terminal_au / dH_terminal_au    端部 A·U 惩罚
        dG_sym / dH_sym                    对称性校正（回文自互补双链）
        dG_core_canonical / dH_core_canonical   基础杂交能（上述之和）
        dG_mismatch_base / dG_mismatch_context / dG_mismatch_positional_extra
        dG_mismatch_total                  错配修正（逐项/汇总，ΔH 近似为 0）
        dG_flank_5p_dangle / dG_flank_3p_dangle / dG_flank_total
                                            侧翼悬垂堆叠（稳定化，负值）
        dG_mismatch_by_position            逐位错配惩罚（含位置带加权，B×19）
        n_mismatch                         每条序列错配列数
        pair_type_columns                  每列配对类型编码 0=WC 1=G:U 2=MM
        dG_total / dH_total / dS_total / Tm  最终四项热力学参数
        （说明：dG_mismatch_positional_extra 为"位置带放大附加项"=Σ 惩罚×(权重−1)，
          seed 与 central 带为正、g1 带为负(宽容扣减)；旧键名 dG_mismatch_seed_extra
          保留为同值别名，供历史脚本兼容）
    单条模式另附 mismatch_sites 人类可读的逐位描述列表。

设计约定（重要）
----------------
    1. "理想互补参考"：对给定的 mRNA 靶标区 window，其完美反义链
       perfect_guide = rc(window)。基础杂交能（canonical 各项）按此理想双链
       计算，因此对"同一靶窗的所有突变体"是同一个常量 —— 这正是项目里
       "ΔΔG = 突变体 ΔG − 野生型 ΔG" 的语义基础：
            ΔG(实际双链) = ΔG(理想双链, 常量) + 错配修正项 + 侧翼项
       于是 突变体相对野生型的能量差 ΔΔG 完全由"错配修正项"贡献，
       便于逐位归因（对应项目任务 11/12/16 的 ΔΔG_by_site 分析）。
    2. 最近邻堆叠参数（AA/UU、AC/UG、… 等 16 项，以 5'→3' 上链二核苷酸
       为键）取 RNA 双链 NN 参数（Xia et al. 1998；与 OligoFormer 原版
       热力学模块 thermodynamics_calculator.py 数值完全一致）。
    3. 错配/G:U 摆动/侧翼悬垂等修正项的默认系数为本模块内置的"近似可调
       默认值"（ΔH 修正默认取 0，即以熵项近似、修正主要作用于 ΔG/排序），
       结构遵循 Turner 规则 NN 框架（错配区段起始项 + 闭合配对上下文修正）。
       如需采用 Turner 2004 完整 1×1 内环实测表，可在 __init__ 传入
       overrides 或在实例上调用 override_params() 整体替换，不影响接口。
    4. 碱基序号：siRNA 位置均按 5'→3' 从 1 起计（g1 为 5' 端）；
       反平行配对关系：guide[p]（0-based）与 window[18-p] 配对，
       即 guide 5' 端对着 window 3' 端。

参考文献
--------
    [1] Xia T., SantaLucia J., et al. Thermodynamic parameters for an expanded
        nearest-neighbor model for formation of RNA duplexes with Watson-Crick
        base pairs. Biochemistry, 1998, 37: 14719-14735.
    [2] SantaLucia J., Hicks D. The thermodynamics of DNA structural motifs.
        Annu. Rev. Biophys. Biomol. Struct., 2004, 33: 415-440.
        （NN 模型框架与 Tm / 盐浓度校正公式出处）
    [3] Mathews D.H., Sabina J., Zuker M., Turner D.H. Expanded sequence
        dependence of thermodynamic parameters improves prediction of RNA
        secondary structure. J. Mol. Biol., 1999, 288: 911-940.
        （错配/内环 Turner 规则背景）
    [4] OligoFormer 仓库 scripts/infer.py 中 calculate_td 与
        thermodynamics_calculator.py（接口与典型 NN 参数对照来源）

运行测试
--------
    python thermo_calculator.py
    末尾自带 3 组测试：
      * 用例 1：完全互补 19 nt siRNA + 57 nt mRNA（全套参数输出）
      * 用例 2：1 个 G:U 摆动错配于中央切割敏感区（g10，位置带 ×2.5）
      * 用例 3：1 个种子区单碱基错配（seed 带 ×2.0），对照 3' 端基准带（内位点16）
      * 用例 4：同类型 YY 错配在 g1(×0.5) 与 central(×2.5) 的宽容/敏感对比
      * 用例 5：批量计算 + trainable 端到端反向传播（4 个位置带权重均有梯度）
================================================================================
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Union

import torch
import torch.nn as nn
import torch.nn.functional as F

# ==============================================================================
# 一、全局常量与基础工具
# ==============================================================================

# 碱基字符顺序：A=0, C=1, G=2, U=3（内部统一使用 U，T 视为 U）
BASE_ORDER: str = "ACGU"
BASE2IDX: Dict[str, int] = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3}

# 互补碱基索引：A<->U, C<->G
RC_IDX: List[int] = [3, 2, 1, 0]

# 嘌呤 / 嘧啶：A、G 为嘌呤；C、U 为嘧啶
PURINES = (0, 2)
PYRIMIDINES = (1, 3)

# RNA 双链最近邻堆叠参数（37 ℃, 1 M NaCl；单位 kcal/mol）
# 键 = 5'->3' 上链二核苷酸；数值与 OligoFormer thermodynamics_calculator.py 一致
STACK_DG: Dict[str, float] = {
    "AA": -0.93, "AC": -2.24, "AG": -2.08, "AU": -1.10,
    "CA": -2.11, "CC": -3.26, "CG": -2.36, "CU": -2.08,
    "GA": -2.35, "GC": -3.42, "GG": -3.26, "GU": -2.24,
    "UA": -1.33, "UC": -2.35, "UG": -2.11, "UU": -0.93,
}
STACK_DH: Dict[str, float] = {
    "AA": -6.82, "AC": -11.40, "AG": -10.48, "AU": -9.38,
    "CA": -10.44, "CC": -13.39, "CG": -10.64, "CU": -10.48,
    "GA": -12.44, "GC": -14.88, "GG": -13.39, "GU": -11.40,
    "UA": -7.69, "UC": -12.44, "UG": -10.44, "UU": -6.82,
}

# 双链总能量中的常数项（kcal/mol；沿用 OligoFormer 原模块口径）
INIT_DG: float = 4.09        # 双链起始项（成核惩罚，正值）
INIT_DH: float = 3.61
END_AU_DG: float = 0.45      # 端部 A·U 对惩罚（每端）
END_AU_DH: float = 3.72
SYM_DG: float = 0.43         # 自互补（回文）双链对称性校正
SYM_DH: float = 0.0

# 默认修正系数（kcal/mol）—— 结构遵循 Turner NN 规则，数值为可调默认值：
# 错配列基础惩罚按"非经典配对类型"分类；G:U 摆动单独给较小系数；
# 位置带权重：错配惩罚按 guide 位点分区加权（见类内 3.4 节）
DEFAULT_CORRECTION = {
    "pen_wobble": 0.30,      # G:U 摆动（基础惩罚，远小于真错配）
    "pen_yy": 0.90,          # 嘧啶-嘧啶错配（C:C, C:U, U:U）
    "pen_pp": 1.35,          # 嘌呤-嘌呤错配（A:A, A:G, G:G）
    "pen_py": 1.55,          # 嘌呤-嘧啶非 WC 错配（A:C, C:A）
    "pen_weak_close": 0.30,  # 错配区段由弱闭合对（A·U/U·A）封闭时的附加惩罚
    "pen_run_init": 0.50,    # 连续错配区段（长度≥2）起始项（每区段一次）
    "seed_weight": 2.0,      # 位置带权重：seed 区（guide 2–8 位）放大系数
    "zone_g1": 0.5,          # 位置带权重：g1（5' 端）—— 宽容（错配容忍高）
    "zone_central": 2.5,     # 位置带权重：central（guide 9–12 切割敏感区）—— 最强
    "zone_z3": 1.0,          # 位置带权重：3' 端（guide 13–19）—— 基准
}

# 侧翼悬垂堆叠默认强度（正值表示稳定化幅度，计入时为负号）：
# 行 = 端部配对类型（0=强闭合 GC/CG, 1=弱闭合 AU/UA, 2=非经典/摆动末端）
# 列 = 悬垂碱基（0=嘧啶 C/U, 1=嘌呤 A/G）
DANGLE_STRENGTH: List[List[float]] = [
    [0.15, 0.35],   # 强 GC 末端：嘌呤悬垂堆叠更强
    [0.10, 0.25],   # 弱 AU 末端
    [0.00, 0.05],   # 摆动/错配末端：几乎无稳定悬垂
]


def rc_seq(seq: str) -> str:
    """返回 RNA 序列的反向互补链（5'->3'）。"""
    comp = {"A": "U", "U": "A", "C": "G", "G": "C", "T": "A"}
    return "".join(comp[b] for b in reversed(seq.upper()))


def rc_indices(idx: torch.Tensor) -> torch.Tensor:
    """张量化反向互补：按碱基索引查 RC_IDX 表。"""
    rc = torch.tensor(RC_IDX, dtype=torch.long, device=idx.device)
    return rc[idx]


def base_char(idx) -> str:
    """碱基索引 -> 字符（解码用）。"""
    return BASE_ORDER[int(idx)]


# ==============================================================================
# 二、主计算类
# ==============================================================================

class siRNAThermoCalculator(nn.Module):
    """
    siRNA–mRNA 杂交热力学参数计算器（PyTorch / 可微分 / 支持批量）。

    核心思路（分解 + 汇总）：
        dG_total = dG_core_canonical        # 理想互补参考双链的 NN 能量
                 + dG_mismatch_total        # 错配/G:U 修正（含种子放大）
                 + dG_flank_total           # 两侧翼悬垂堆叠（稳定化，负值）

    典型用法
    --------
        calc = siRNAThermoCalculator()
        # 单条
        res = calc.calculate(sirna_19, mrna_57)
        # 批量
        res_b = calc.calculate_batch([s1, s2], [m1, m2])
    """

    # ------------------------------------------------------------------
    # 构造：加载内置参数表；支持外部配置与"可训练"开关
    # ------------------------------------------------------------------
    def __init__(
        self,
        device: Union[str, torch.device] = "cpu",
        dtype: torch.dtype = torch.float32,
        trainable: bool = False,
        overrides: Optional[Dict[str, float]] = None,
    ) -> None:
        """
        Parameters
        ----------
        device : 计算设备（'cpu' / 'cuda'）
        dtype   : 张量精度
        trainable : True 时把堆叠表与修正系数注册为 nn.Parameter，
                    可嵌入深度学习模型做端到端梯度训练（需自行配优化器）。
        overrides : 可选，形如 {"pen_wobble": 0.5, "seed_weight": 3.0}，
                    覆盖内置修正系数；堆叠表键示例 "GC": -3.42（同名覆盖）。
        """
        super().__init__()
        self.device = torch.device(device)
        self.dtype = dtype
        self.training_params = bool(trainable)

        # ---- 1. 最近邻堆叠参数表：16 项，按 A,C,G,U 顺序排成 (16,) ----
        stack_dg = torch.zeros(16, dtype=dtype)
        stack_dh = torch.zeros(16, dtype=dtype)
        for i, b1 in enumerate(BASE_ORDER):
            for j, b2 in enumerate(BASE_ORDER):
                key = b1 + b2
                stack_dg[i * 4 + j] = STACK_DG[key]
                stack_dh[i * 4 + j] = STACK_DH[key]
        # 支持外部覆盖单个堆叠参数
        for key, val in (overrides or {}).items():
            if len(key) == 2 and key.upper() in STACK_DG:
                i, j = BASE2IDX[key[0].upper()], BASE2IDX[key[1].upper()]
                stack_dg[i * 4 + j] = val

        # ---- 2. 常数项 ----
        self.init_dg = float(INIT_DG)
        self.init_dh = float(INIT_DH)
        self.end_au_dg = float(END_AU_DG)
        self.end_au_dh = float(END_AU_DH)
        self.sym_dg = float(SYM_DG)
        self.sym_dh = float(SYM_DH)

        # ---- 3. 错配修正系数（可训练时注册为 Parameter）----
        cfg = dict(DEFAULT_CORRECTION)
        cfg.update({k: v for k, v in (overrides or {}).items()
                    if k in DEFAULT_CORRECTION})
        # 错配基础惩罚按列的类型编码查表：0=G:U摆动, 1=YY, 2=PP, 3=PY
        self.pen_by_cat = self._make_param(
            [cfg["pen_wobble"], cfg["pen_yy"], cfg["pen_pp"], cfg["pen_py"]],
            trainable)
        self.pen_weak_close = self._make_scalar(cfg["pen_weak_close"], trainable)
        self.pen_run_init = self._make_scalar(cfg["pen_run_init"], trainable)
        # 位置带权重系数（guide 5'->3' 分区：g1 / seed / central / 3' 端）
        self.seed_weight = self._make_scalar(cfg["seed_weight"], trainable)
        self.zone_g1 = self._make_scalar(cfg["zone_g1"], trainable)
        self.zone_central = self._make_scalar(cfg["zone_central"], trainable)
        self.zone_z3 = self._make_scalar(cfg["zone_z3"], trainable)
        # 侧翼悬垂堆叠强度表 (3 端部类型 × 2 悬垂碱基类)
        self.dangle_strength = self._make_param(
            [DANGLE_STRENGTH[i][j]
             for i in range(3) for j in range(2)], trainable)

        # ---- 4. 堆叠表本身（trainable 时一并训练）----
        self.stack_dg = self._make_param(stack_dg.tolist(), trainable)
        self.stack_dh = self._make_param(stack_dh.tolist(), trainable)

        # ---- 5. 熔解温度相关可调参数（不参与梯度，普通 float）----
        self.tm_strand_conc: float = 1e-7       # 总链浓度 (M)，两链等浓度
        self.tm_sodium: float = 1.0             # Na+ 等效浓度 (M)；1 M 时校正为 0
        self.tm_gas_const: float = 1.987        # R, cal/(mol·K)

    # ------------------------------------------------------------------
    # 参数工具
    # ------------------------------------------------------------------
    def _make_scalar(self, value: float, trainable: bool):
        """单值参数：可训练时作为 Parameter（形状 [1]）。"""
        t = torch.tensor([float(value)], dtype=self.dtype)
        if trainable:
            return nn.Parameter(t.to(self.device))
        return t.to(self.device)

    def _make_param(self, values: Sequence[float], trainable: bool):
        """列表参数：可训练时作为 Parameter。"""
        t = torch.tensor([float(v) for v in values], dtype=self.dtype)
        if trainable:
            return nn.Parameter(t.to(self.device))
        return t.to(self.device)

    def override_params(self, **kwargs: float) -> "siRNAThermoCalculator":
        """
        运行期整体替换参数（不破坏计算图）。
        支持：pen_wobble / pen_yy / pen_pp / pen_py / pen_weak_close /
              pen_run_init / seed_weight / zone_g1 / zone_central / zone_z3 /
              dangle_strength / stack_dg / stack_dh
              （堆叠表传长度为 16 的序列或 {"GC": val, ...} 字典）
        """
        with torch.no_grad():
            def _fill(t, src):
                t.copy_(torch.as_tensor(src, dtype=self.dtype, device=t.device))

            for key, val in kwargs.items():
                if key in ("pen_wobble", "pen_yy", "pen_pp", "pen_py"):
                    order = {"pen_wobble": 0, "pen_yy": 1, "pen_pp": 2, "pen_py": 3}
                    self.pen_by_cat[order[key]] = float(val)
                elif key in ("pen_weak_close", "pen_run_init", "seed_weight",
                             "zone_g1", "zone_central", "zone_z3"):
                    getattr(self, key).fill_(float(val))
                elif key == "dangle_strength":
                    _fill(self.dangle_strength, list(val))
                elif key in ("stack_dg", "stack_dh"):
                    if isinstance(val, dict):
                        t = torch.zeros(16, dtype=self.dtype)
                        for k, v in val.items():
                            t[BASE2IDX[k[0]] * 4 + BASE2IDX[k[1]]] = v
                        val = t
                    _fill(getattr(self, key), list(val))
                else:
                    raise KeyError(f"未知参数名：{key}")
        return self

    # ------------------------------------------------------------------
    # 输入校验与编码（中文报错）
    # ------------------------------------------------------------------
    @staticmethod
    def _normalize_str(seq: str, name: str, length: int) -> str:
        """字符串序列规范化：去空白、转大写、T->U，并校验长度与碱基。"""
        if not isinstance(seq, str):
            raise TypeError(f"{name} 应为字符串，实际为 {type(seq).__name__}")
        s = "".join(seq.split()).upper().replace("T", "U")
        if len(s) != length:
            raise ValueError(
                f"{name} 长度必须为 {length} nt（当前 {len(s)} nt）。"
                f"siRNA 反义链固定 19 nt；mRNA 固定 57 nt"
                f"（5'侧翼19 + 靶标19 + 3'侧翼19）。")
        bad = sorted({c for c in s if c not in BASE_ORDER})
        if bad:
            raise ValueError(f"{name} 含非法碱基：{bad}；仅允许 A/C/G/U/T。")
        return s

    def _encode_strs(self, sirnas: List[str], mrnas: List[str]):
        """批量字符串 -> (guide_idx, mrna_idx) 长整型张量 (B,19)/(B,57)。"""
        if len(sirnas) != len(mrnas):
            raise ValueError(
                f"siRNA 与 mRNA 数量不一致：{len(sirnas)} vs {len(mrnas)}")
        g = torch.tensor(
            [[BASE2IDX[c] for c in self._normalize_str(s, "siRNA", 19)]
             for s in sirnas], dtype=torch.long)
        m = torch.tensor(
            [[BASE2IDX[c] for c in self._normalize_str(s, "mRNA", 57)]
             for s in mrnas], dtype=torch.long)
        return g.to(self.device), m.to(self.device)

    @staticmethod
    def _check_idx_tensor(x: torch.Tensor, name: str, length: int) -> torch.Tensor:
        """碱基索引张量校验：末维必须为 length，值必须在 0..3。"""
        if not isinstance(x, torch.Tensor):
            raise TypeError(f"{name} 应为 torch.Tensor（碱基索引），"
                            f"实际为 {type(x).__name__}")
        if x.dim() == 1:
            x = x.unsqueeze(0)               # 单条 -> 批维 1
        if x.shape[-1] != length:
            raise ValueError(f"{name} 末维长度必须为 {length}（当前 {x.shape[-1]}）")
        if (x < 0).any() or (x > 3).any():
            raise ValueError(f"{name} 索引值必须在 0..3（A=0,C=1,G=2,U=3）")
        return x.long()

    # ------------------------------------------------------------------
    # 内部核心计算（纯张量 / 全程可微分）
    # ------------------------------------------------------------------
    def _forward_tensors(self, guide: torch.Tensor, mrna: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        核心前向：输入碱基索引张量 (B,19) 与 (B,57)，返回全部中间与最终
        能量张量。guide 第 p 列（0-based）与靶标区 window[18-p] 反平行配对，
        即 guide 5' 端对着 window 3' 端。
        """
        B = guide.shape[0]
        dev = guide.device

        # ---- 0. 取出靶标区并建立"配对列"坐标系 ----
        # window: mRNA[19:38]（19 nt，5'->3'）
        window = mrna[:, 19:38]
        # b[p] = window[18-p]：guide[p] 的配对碱基（随列序号排列的 mRNA 侧碱基）
        b = torch.flip(window, dims=[1])
        # 理想互补引导链 perfect_guide[p] = rc(b[p])
        gp = rc_indices(b)                       # (B,19)

        # ---- 1. 逐列配对类型：0=WC, 1=G:U 摆动, 2=真错配 ----
        wc = guide == gp                         # Watson-Crick 完全互补
        # G:U / U:G 摆动（索引 G=2, U=3）
        gu = ((guide == 2) & (b == 3)) | ((guide == 3) & (b == 2))
        mm = ~(wc | gu)                          # 其余全部为真错配
        # 错配类型编码用于查基础惩罚表：0=G:U, 1=YY, 2=PP, 3=PY
        pur_g = (guide == 0) | (guide == 2)      # guide 碱基为嘌呤
        pur_b = (b == 0) | (b == 2)              # mRNA 侧碱基为嘌呤
        same_class = pur_g == pur_b              # 同为嘌呤 或 同为嘧啶
        cat = torch.zeros_like(guide, dtype=torch.long)
        cat = torch.where(gu, torch.zeros_like(cat), cat)        # 摆动 -> 0
        cat = torch.where(mm & same_class & pur_g, torch.full_like(cat, 2), cat)  # PP
        cat = torch.where(mm & same_class & (~pur_g), torch.full_like(cat, 1), cat)  # YY
        cat = torch.where(mm & (~same_class), torch.full_like(cat, 3), cat)       # PY

        # 配对类型编码（输出用）：0=WC 1=G:U 2=MM
        pair_code = torch.zeros_like(guide, dtype=torch.long)
        pair_code = torch.where(gu, torch.ones_like(pair_code), pair_code)
        pair_code = torch.where(mm, torch.full_like(pair_code, 2), pair_code)

        # ---- 2. 基础杂交能（canonical）：以理想互补链 gp 为参考 ----
        # 相邻两列均配对即构成一个最近邻堆叠步（共 18 步）
        dimer = gp[:, :-1] * 4 + gp[:, 1:]       # (B,18) 堆叠键索引
        stack_dg_sum = self.stack_dg[dimer].sum(dim=1)
        stack_dh_sum = self.stack_dh[dimer].sum(dim=1)

        # 端部 A·U 惩罚：窗口两端的碱基属于 {A,U} 即构成 A·U 末端（含 U·A）
        ends = torch.cat([window[:, :1], window[:, -1:]], dim=1)   # (B,2)
        is_au_end = (ends == 0) | (ends == 3)                      # A 或 U
        n_au_end = is_au_end.sum(dim=1).to(self.dtype)             # (B,)

        # 对称性校正：window 为回文（window == rc(window)）时双链自互补
        pal = (window == rc_indices(torch.flip(window, dims=[1]))).all(dim=1)

        dG_core = (self.init_dg + stack_dg_sum
                   + self.end_au_dg * n_au_end
                   + torch.where(pal, torch.tensor(self.sym_dg, device=dev,
                                                   dtype=self.dtype),
                                 torch.zeros(B, dtype=self.dtype, device=dev)))
        dH_core = (self.init_dh + stack_dh_sum
                   + self.end_au_dh * n_au_end
                   + torch.where(pal, torch.tensor(self.sym_dh, device=dev,
                                                   dtype=self.dtype),
                                 torch.zeros(B, dtype=self.dtype, device=dev)))

        # ---- 3. 错配能量修正（逐位分解 + 汇总）----
        # 3.1 每列基础惩罚（按类型编码查表 pen_by_cat）
        base_pen_col = self.pen_by_cat[cat].clone()
        base_pen_col = torch.where(wc, torch.zeros_like(base_pen_col),
                                   base_pen_col)   # WC 列惩罚为 0

        # 3.2 弱闭合配对上下文：错配列两侧相邻的 WC 列若是弱 A·U/U·A，
        #     则附加惩罚（Turner 内环"弱闭合对"概念）
        weak_wc = wc & ((guide == 0) | (guide == 3))     # WC 且为弱 A·U/U·A
        left_weak = F.pad(weak_wc[:, :-1], (1, 0))       # 左侧列
        right_weak = F.pad(weak_wc[:, 1:], (0, 1))       # 右侧列
        n_weak_close = ((left_weak & ~wc).to(self.dtype)
                        + (right_weak & ~wc).to(self.dtype))  # 错配列两侧计数
        ctx_pen_col = self.pen_weak_close * n_weak_close

        # 3.3 连续错配区段起始项：仅对长度≥2 的区段加一次
        noncan = ~wc
        starts = noncan & ~F.pad(noncan[:, :-1], (1, 0))      # 区段起点
        has_next = F.pad(noncan[:, 1:], (0, 1))               # 下一列是否错配
        run_start_ge2 = starts & has_next                     # 长度≥2 的起点
        run_pen_col = self.pen_run_init * run_start_ge2.to(self.dtype)

        # 3.4 位置带权重（位置权重）：guide 5'->3' 位点分区，错配惩罚按带加权
        #     g1(位1, 宽容) | seed(位2-8, 放大) | central(位9-12, 切割敏感区, 最强)
        #     | 3'端(位13-19, 基准)
        # 因子向量用"1 + 掩码×(权重−1)"外推法构造，保持与 Parameter 的计算图
        # 连接（trainable 时 zone 权重可获得非零梯度）
        pos = torch.arange(19, device=dev, dtype=self.dtype)   # 0..18 = 位1..19
        sel_g1 = (pos == 0).to(self.dtype)                     # 位 1
        sel_seed = ((pos >= 1) & (pos <= 7)).to(self.dtype)    # 位 2-8
        sel_central = ((pos >= 8) & (pos <= 11)).to(self.dtype)  # 位 9-12
        sel_z3 = (pos >= 12).to(self.dtype)                    # 位 13-19
        factor_vec = (torch.ones(19, device=dev, dtype=self.dtype)
                      + sel_g1 * (self.zone_g1 - 1.0)
                      + sel_seed * (self.seed_weight - 1.0)
                      + sel_central * (self.zone_central - 1.0)
                      + sel_z3 * (self.zone_z3 - 1.0))          # (19,)

        # 每列未加权的惩罚 = 基础 + 上下文 + 区段起始
        term_col_unweighted = base_pen_col + ctx_pen_col + run_pen_col
        # 逐位加权后的惩罚（含位置带加权；输出 by_position 用）
        term_col_weighted = term_col_unweighted * factor_vec    # (B,19)*(19,)广播
        # 分项汇总（输出用）
        dG_mismatch_base = base_pen_col.sum(dim=1)
        dG_mismatch_context = (ctx_pen_col + run_pen_col).sum(dim=1)
        # 位置带放大附加项 = Σ 惩罚×(权重−1)：seed/central 为正，g1(权重<1)为负
        dG_mismatch_positional_extra = (
            term_col_unweighted * (factor_vec - 1.0)).sum(dim=1)
        dG_mismatch_total = term_col_weighted.sum(dim=1)
        n_mismatch = noncan.sum(dim=1)

        # 错配修正的焓近似为 0（修正主要作用于自由能/排序，见模块头注释）

        # ---- 4. 侧翼悬垂堆叠（dangling end）----
        # 末端配对类型：0=强 GC/CG, 1=弱 AU/UA, 2=摆动/错配末端
        def _end_class(pair_wc: torch.Tensor, pair_gu: torch.Tensor,
                       g0: torch.Tensor) -> torch.Tensor:
            """末端列配对类型 -> 0/1/2（查悬垂强度表行）。"""
            cls = torch.where(pair_wc, torch.zeros_like(pair_wc, dtype=torch.long),
                              torch.full_like(pair_wc, 2, dtype=torch.long))
            # WC 末端：A/U 末端 -> 1（弱），GC 末端 -> 0（强）
            weak = (g0 == 0) | (g0 == 3)
            cls = torch.where(pair_wc & weak,
                              torch.ones_like(cls), cls)
            return cls

        # 5' 端（guide[0]，窗口 3' 端 window[18]），相邻侧翼碱基 mRNA[38]
        g5, w5 = guide[:, 0], window[:, 18]
        wc5, gu5 = wc[:, 0], gu[:, 0]
        cls5 = _end_class(wc5, gu5, g5)
        # 悬垂碱基类：0=嘧啶, 1=嘌呤
        dg5_pur = ((mrna[:, 38] == 0) | (mrna[:, 38] == 2)).to(torch.long)
        # 3' 端（guide[18]，窗口 5' 端 window[0]），相邻侧翼碱基 mRNA[18]
        g3, w3 = guide[:, 18], window[:, 0]
        wc3, gu3 = wc[:, 18], gu[:, 0]
        cls3 = _end_class(wc3, gu3, g3)
        dg3_pur = ((mrna[:, 18] == 0) | (mrna[:, 18] == 2)).to(torch.long)

        # 悬垂强度表 (3*2)：索引 = 末端类型 * 2 + 悬垂嘌呤标记
        idx5 = cls5 * 2 + dg5_pur
        idx3 = cls3 * 2 + dg3_pur
        dangle5 = self.dangle_strength[idx5]       # 正值（稳定化幅度）
        dangle3 = self.dangle_strength[idx3]
        # 计入时取负号（悬垂堆叠稳定双链）；焓近似为 0
        dG_flank_5p = -dangle5
        dG_flank_3p = -dangle3
        dG_flank_total = dG_flank_5p + dG_flank_3p

        # ---- 5. 汇总四项热力学参数 ----
        dG_total = dG_core + dG_mismatch_total + dG_flank_total
        dH_total = dH_core                            # 修正项 ΔH≈0（见头注释）
        # ΔS = (ΔH − ΔG) / T，37 ℃=310.15 K；单位 cal/(mol·K)
        dS_total = (dH_total - dG_total) * 1000.0 / 310.15

        # ---- 6. 熔解温度 Tm（SantaLucia & Hicks, 2004）----
        # 盐浓度校正（Na=1 M 时 ln 1 = 0，不生效）：
        #   ΔS(Na) = ΔS + 0.368 × (N/2) × ln[Na+]，N= 双链磷酸基数（19 mer 取 36）
        n_phos = 36
        salt_corr = 0.368 * (n_phos / 2.0) * math.log(self.tm_sodium)
        dS_eff = dS_total + salt_corr
        # 双分子非自互补双链：Tm = ΔH / (ΔS + R·ln(C_T/4)) − 273.15
        x_conc = 1.0 if False else 4.0
        r = self.tm_gas_const
        ct = self.tm_strand_conc
        denom = dS_eff + r * math.log(ct / x_conc)
        tm_k = torch.where(
            denom.abs() > 1e-9,
            dH_total * 1000.0 / denom,
            torch.full_like(denom, float("nan")))
        tm = tm_k - 273.15

        # ---- 7. 组装返回字典 ----
        return {
            # 分项：基础杂交能（canonical）
            "dG_core_stack": stack_dg_sum,
            "dH_core_stack": stack_dh_sum,
            "dG_init_term": torch.full((B,), self.init_dg, device=dev, dtype=self.dtype),
            "dH_init_term": torch.full((B,), self.init_dh, device=dev, dtype=self.dtype),
            "dG_terminal_au": self.end_au_dg * n_au_end,
            "dH_terminal_au": self.end_au_dh * n_au_end,
            "dG_sym": torch.where(pal, torch.tensor(self.sym_dg, device=dev,
                                                    dtype=self.dtype),
                                  torch.zeros(B, dtype=self.dtype, device=dev)),
            "dH_sym": torch.zeros(B, dtype=self.dtype, device=dev),
            "dG_core_canonical": dG_core,
            "dH_core_canonical": dH_core,
            # 分项：错配修正
            "dG_mismatch_base": dG_mismatch_base,
            "dG_mismatch_context": dG_mismatch_context,
            "dG_mismatch_positional_extra": dG_mismatch_positional_extra,
            "dG_mismatch_seed_extra": dG_mismatch_positional_extra,  # 兼容别名
            "dG_mismatch_total": dG_mismatch_total,
            "dG_mismatch_by_position": term_col_weighted,
            "n_mismatch": n_mismatch,
            "pair_type_columns": pair_code,
            # 分项：侧翼悬垂
            "dG_flank_5p_dangle": dG_flank_5p,
            "dG_flank_3p_dangle": dG_flank_3p,
            "dG_flank_total": dG_flank_total,
            # 汇总：四项最终参数
            "dG_total": dG_total,
            "dH_total": dH_total,
            "dS_total": dS_total,
            "Tm": tm,
        }

    # ------------------------------------------------------------------
    # 对外接口 1：单条计算 -> Python 标量 dict（附逐位描述）
    # ------------------------------------------------------------------
    def calculate(self, seq_sirna: str, seq_mrna: str) -> Dict:
        """单条序列计算：返回全部热力学参数（Python 标量）与逐位错配描述。"""
        g, m = self._encode_strs([seq_sirna], [seq_mrna])
        raw = self._forward_tensors(g, m)
        out: Dict = {}
        for k, v in raw.items():
            v = v.detach().cpu()
            if v.numel() == 1:
                # (B=1,) / 0 维 -> Python 标量
                out[k] = float(v.reshape(-1)[0])
            else:
                # (1,19) -> 展开成 19 元素列表
                lst = v.tolist()
                out[k] = lst[0] if (len(lst) == 1 and isinstance(lst[0], list)) else lst

        # ---- 人类可读的逐位错配描述（供打印/报告）----
        sirna = self._normalize_str(seq_sirna, "siRNA", 19)
        mrna = self._normalize_str(seq_mrna, "mRNA", 57)
        window = mrna[19:38]
        pair_types = {0: "WC", 1: "G:U摆动", 2: "MM错配"}
        # 位置带：guide 1-based 位点 -> 名称
        def _zone(p1):
            if p1 == 1:
                return "g1"
            if p1 <= 8:
                return "seed(g2-8)"
            if p1 <= 12:
                return "central(g9-12)"
            return "3'(g13-19)"
        sites = []
        for p in range(19):
            code = out["pair_type_columns"][p]
            if code == 0:
                continue
            g_b = sirna[p]
            m_b = window[18 - p]
            in_seed = 1 <= p <= 7
            pen = out["dG_mismatch_by_position"][p]
            sites.append({
                "位置(1-based)": p + 1,
                "siRNA碱基": g_b,
                "mRNA配对碱基": m_b,
                "配对类型": pair_types[int(code)],
                "位置带": _zone(p + 1),
                "是否种子区(g2-g8)": in_seed,
                "惩罚ΔG(kcal/mol)": round(float(pen), 4),
            })
        out["mismatch_sites"] = sites
        return out

    # ------------------------------------------------------------------
    # 对外接口 2：批量计算 -> torch.Tensor dict（可接梯度）
    # ------------------------------------------------------------------
    def calculate_batch(
        self,
        seq_sirna: Union[str, List[str], torch.Tensor],
        seq_mrna: Union[str, List[str], torch.Tensor],
    ) -> Dict[str, torch.Tensor]:
        """
        批量计算。输入：
          * 均为字符串 / 字符串列表；
          * 或均为碱基索引张量（(B,19)/(B,57)，单条自动补批维）。
        返回 dict：键同 calculate()，值为 (B,) 或 (B,19) 张量（未 detach，
        需要梯度时可直接 backward）。
        """
        if isinstance(seq_sirna, torch.Tensor) or isinstance(seq_mrna, torch.Tensor):
            if not (isinstance(seq_sirna, torch.Tensor)
                    and isinstance(seq_mrna, torch.Tensor)):
                raise TypeError("张量模式下 seq_sirna 与 seq_mrna 必须同时为张量")
            g = self._check_idx_tensor(seq_sirna, "seq_sirna", 19)
            m = self._check_idx_tensor(seq_mrna, "seq_mrna", 57)
            if g.shape[0] != m.shape[0]:
                raise ValueError(f"批次数量不一致：{g.shape[0]} vs {m.shape[0]}")
            g, m = g.to(self.device), m.to(self.device)
        else:
            if isinstance(seq_sirna, str):
                seq_sirna = [seq_sirna]
            if isinstance(seq_mrna, str):
                seq_mrna = [seq_mrna]
            g, m = self._encode_strs(list(seq_sirna), list(seq_mrna))
        return self._forward_tensors(g, m)

    # ------------------------------------------------------------------
    # 便捷工具：dG 张量前向（供外部嵌入模型直接调用）
    # ------------------------------------------------------------------
    def forward(self, guide: torch.Tensor, mrna: torch.Tensor) -> torch.Tensor:
        """nn.Module 前向：返回 (B,) 的 dG_total（最常用标量）。"""
        g = self._check_idx_tensor(guide, "guide", 19)
        m = self._check_idx_tensor(mrna, "mrna", 57)
        return self._forward_tensors(g.to(self.device), m.to(self.device))["dG_total"]


# ==============================================================================
# 三、内置自测 / 演示（python thermo_calculator.py 直接运行）
# ==============================================================================

def _fmt(x, digits: int = 4) -> str:
    """数值格式化：NaN 显示为 '--'。"""
    if isinstance(x, torch.Tensor):
        x = float(x.item()) if x.numel() == 1 else x.detach().cpu().tolist()
    if isinstance(x, float):
        if math.isnan(x):
            return "--"
        return f"{x:.{digits}f}"
    return str(x)


def _run_tests() -> None:
    """三组内置测试用例。"""
    print("=" * 78)
    print("siRNA-mRNA 杂交热力学计算模块 —— 内置测试")
    print("=" * 78)

    # ---------- 公共序列 ----------
    # mRNA 靶标区（19 nt，位于 57 nt 正中 mRNA[19:38]）
    window = "GCAGCACAAUGAUGAGUGC"
    # 5'/3' 侧翼（各 19 nt；两端各只有紧邻的 1 个碱基参与悬垂堆叠）
    flank5 = ("ACGU" * 5)[:19]
    flank3 = ("UGCA" * 5)[:19]
    mrna = flank5 + window + flank3
    assert len(mrna) == 57, "mRNA 测试序列长度异常"

    # 完全互补的反义链（理想引导链）
    guide_wt = rc_seq(window)
    assert len(guide_wt) == 19

    calc = siRNAThermoCalculator()
    print(f"\nmRNA 序列(57 nt): {mrna}")
    print(f"靶标区 window  : {window}")
    print(f"理想引导链(WT) : {guide_wt}")
    print(f"位置带权重     : g1={float(calc.zone_g1):.1f} | "
          f"seed(g2-8)={float(calc.seed_weight):.1f} | "
          f"central(g9-12)={float(calc.zone_central):.1f} | "
          f"3'(g13-19)={float(calc.zone_z3):.1f}")
    print(f"说明: 侧翼只提供悬垂堆叠; 完全互补时错配修正=0\n")

    # ---------- 用例 1：完全互补 ----------
    print("-" * 78)
    print("【用例 1】完全互补配对（WT 引导链）")
    print("-" * 78)
    res1 = calc.calculate(guide_wt, mrna)
    for k in ("dG_core_stack", "dG_init_term", "dG_terminal_au", "dG_sym",
              "dG_core_canonical", "dG_mismatch_total", "dG_flank_total",
              "dG_total", "dH_total", "dS_total", "Tm", "n_mismatch"):
        print(f"  {k:<22s}: {_fmt(res1[k])}")
    print(f"  mismatch_sites : {res1['mismatch_sites']}")

    # ---------- 用例 2：1 个 G:U 摆动错配（位点 10 = 中央切割敏感区 g9-12）----------
    print("\n" + "-" * 78)
    print("【用例 2】1 个 G:U 摆动错配（guide 第 10 位 A->G, 中央区 g9-12）")
    print("-" * 78)
    g2 = list(guide_wt)
    g2[9] = "G"                       # window[9]=U，guide G<->U 构成 G:U 摆动
    guide_gu = "".join(g2)
    res2 = calc.calculate(guide_gu, mrna)
    print(f"  突变引导链  : {guide_gu}")
    for k in ("dG_mismatch_base", "dG_mismatch_context",
              "dG_mismatch_positional_extra", "dG_mismatch_total",
              "dG_total", "Tm", "n_mismatch"):
        print(f"  {k:<26s}: {_fmt(res2[k])}")
    print(f"  mismatch_sites: {res2['mismatch_sites']}")
    ddg2 = res2["dG_total"] - res1["dG_total"]
    print(f"\n  → ΔΔG(突变体 − WT) = {ddg2:+.4f} kcal/mol"
          f"（摆动 0.30 + 弱闭合上下文 0.30 = 0.60，中央区位置带 ×2.5 → "
          f"{_fmt(res2['dG_mismatch_total'])}；悬垂差 0.00）")

    # ---------- 用例 3：1 个种子区单碱基错配（位点 4），对照在 3' 端基准带 ----------
    print("\n" + "-" * 78)
    print("【用例 3】1 个种子区单碱基错配（guide 第 4 位 C->A, seed 区 g2-g8）")
    print("-" * 78)
    g3 = list(guide_wt)
    g3[3] = "A"                       # window[15]=G，A·G = 嘌呤-嘌呤错配，位于种子区
    guide_seed = "".join(g3)
    res3 = calc.calculate(guide_seed, mrna)
    print(f"  突变引导链  : {guide_seed}")
    for k in ("dG_mismatch_base", "dG_mismatch_context",
              "dG_mismatch_positional_extra", "dG_mismatch_total",
              "dG_total", "Tm", "n_mismatch"):
        print(f"  {k:<26s}: {_fmt(res3[k])}")
    print(f"  mismatch_sites: {res3['mismatch_sites']}")
    ddg3 = res3["dG_total"] - res1["dG_total"]
    print(f"\n  → ΔΔG(突变体 − WT) = {ddg3:+.4f} kcal/mol")
    # 对照组：同类型错配移到 3' 端基准带（内部位点 16，避免端部悬垂干扰）
    g3b = list(guide_wt)
    g3b[15] = "A"                     # window[3]=G，A·G 嘌呤-嘌呤错配，3' 端基准带(内位点)
    guide_ctrl = "".join(g3b)
    res_ctrl = calc.calculate(guide_ctrl, mrna)
    ddg_ctrl = res_ctrl["dG_total"] - res1["dG_total"]
    print(f"  对照组(同型PP、3'端基准带, 内位点16 A·G): ΔΔG = {ddg_ctrl:+.4f} kcal/mol")
    print(f"\n  → seed区 vs 3'端 位置带放大倍数 ≈ {ddg3 / ddg_ctrl:.2f}×"
          f"（seed 权重 {float(calc.seed_weight):.1f} vs 3'端 {float(calc.zone_z3):.1f}；"
          f"两侧弱闭合上下文数量不同, 故倍数约等于而非精确 2.0）")

    # ---------- 用例 4：位置带宽容/敏感对比（同类型 YY，不同位置带）----------
    print("\n" + "-" * 78)
    print("【用例 4】位置带对比：同类型 YY 错配 在 g1(×0.5) vs central(×2.5)")
    print("-" * 78)
    g4a = list(guide_wt)
    g4a[0] = "U"                      # 位1: G->U, window[18]=C -> U·C = YY 错配 (g1 带 ×0.5)
    guide_g1yy = "".join(g4a)
    res_g1yy = calc.calculate(guide_g1yy, mrna)
    ddg_g1yy = res_g1yy["dG_total"] - res1["dG_total"]
    g4b = list(guide_wt)
    g4b[9] = "C"                      # 位10: G->C, window[9]=U -> C·U = YY 错配 (central 带 ×2.5)
    guide_cenyy = "".join(g4b)
    res_cenyy = calc.calculate(guide_cenyy, mrna)
    ddg_cenyy = res_cenyy["dG_total"] - res1["dG_total"]
    print(f"  位1  (g1带 ×0.5,      YY: U·C): ΔΔG = {ddg_g1yy:+.4f} kcal/mol"
          f"  sites={res_g1yy['mismatch_sites']}")
    print(f"  位10 (central带 ×2.5, YY: C·U): ΔΔG = {ddg_cenyy:+.4f} kcal/mol"
          f"  sites={res_cenyy['mismatch_sites']}")
    print(f"\n  → central/g1 同类型错配的惩罚比 ≈ {ddg_cenyy / ddg_g1yy:.2f}×"
          f"（位置带权重比 = {float(calc.zone_central):.1f}/{float(calc.zone_g1):.1f} = "
          f"{float(calc.zone_central) / float(calc.zone_g1):.0f}×，"
          f"上下文差异导致倍数略偏）")

    # ---------- 用例 5（附加）：批量模式 + 端到端梯度 ----------
    print("\n" + "-" * 78)
    print("【用例 5】批量计算 & 可微分性（trainable 模式反向传播）")
    print("-" * 78)
    demo = [guide_wt, guide_gu, guide_seed, guide_ctrl, guide_g1yy]
    batch = calc.calculate_batch(demo, [mrna] * len(demo))
    print(f"  批量 dG_total : {[round(float(x), 4) for x in batch['dG_total']]}")
    print(f"  批量 n_mismatch: {batch['n_mismatch'].tolist()}")

    # trainable 模式：把 4 个位置带权重等当作可学习参数，反传一次看梯度是否存在
    calc_t = siRNAThermoCalculator(trainable=True)
    g_t = calc_t.calculate_batch(demo, [mrna] * len(demo))
    loss = g_t["dG_total"].sum()
    loss.backward()
    grads = {
        name: (p.grad.abs().sum().item() if p.grad is not None else 0.0)
        for name, p in calc_t.named_parameters()
    }
    print(f"  trainable 参数梯度(非零即证明可端到端训练): {grads}")
    print("  （zone_g1/seed_weight/zone_central/zone_z3 均非零 => 位置带权重可端到端学习）")

    print("\n" + "=" * 78)
    print("全部测试运行完毕，无异常。")
    print("=" * 78)


if __name__ == "__main__":
    _run_tests()
