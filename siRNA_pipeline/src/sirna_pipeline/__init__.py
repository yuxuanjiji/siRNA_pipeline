# -*- coding: utf-8 -*-
"""sirna_pipeline —— SFRP1 靶向 siRNA 计算筛选与排序集成工程。

包内 import 方向（防循环依赖）：
    common <── stages <── pipeline
    common 内部：seqio / duplex / nn_tables / mismatch_class / records / stage_io 相互独立或单向引用。

原始任务脚本一律不修改；本包通过 stages/<module>/runner 调用 legacy 核心，见 docs/design/。
"""

__version__ = "0.1.0"
