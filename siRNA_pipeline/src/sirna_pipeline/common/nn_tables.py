# -*- coding: utf-8 -*-
"""核酸最近邻（NN）热力学参数 —— 全工程唯一参数源（只读常量）。

来源与数值逐字一致：
  - 项目搭建/热力学参数计算/thermo_calculator.py（STACK_DG/STACK_DH/INIT/END_AU/SYM，
    其值又与 OligoFormer thermodynamics_calculator.py 一致）
  - 项目搭建/结构检测/structure_detector.py::STACK_DG（仅 ΔG，注释声明与 thermo 完全一致）

用途：
  - 新代码（runner / rank / 派生特征）只引用本表；
  - 结构/热力学 legacy 模块内部仍用各自常量表 —— 整合测试中以本表对 legacy 数值做等价断言，
    未来如统一替换 Turner 2004 实测表只需改本文件 + 阶段性同步。
"""
from __future__ import annotations

# RNA 双链最近邻堆叠参数（37 ℃, 1 M NaCl；kcal/mol）
# 键 = 5'->3' 上链二核苷酸
STACK_DG: dict[str, float] = {
    "AA": -0.93, "AC": -2.24, "AG": -2.08, "AU": -1.10,
    "CA": -2.11, "CC": -3.26, "CG": -2.36, "CU": -2.08,
    "GA": -2.35, "GC": -3.42, "GG": -3.26, "GU": -2.24,
    "UA": -1.33, "UC": -2.35, "UG": -2.11, "UU": -0.93,
}
STACK_DH: dict[str, float] = {
    "AA": -6.82, "AC": -11.40, "AG": -10.48, "AU": -9.38,
    "CA": -10.44, "CC": -13.39, "CG": -10.64, "CU": -10.48,
    "GA": -12.44, "GC": -14.88, "GG": -13.39, "GU": -11.40,
    "UA": -7.69, "UC": -12.44, "UG": -10.44, "UU": -6.82,
}

INIT_DG: float = 4.09       # 双链起始项（成核惩罚）
INIT_DH: float = 3.61
END_AU_DG: float = 0.45     # 端部 A·U 惩罚（每端）
END_AU_DH: float = 3.72
SYM_DG: float = 0.43        # 自互补（回文）双链对称性校正
SYM_DH: float = 0.0
