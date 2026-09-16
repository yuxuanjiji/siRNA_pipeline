# -*- coding: utf-8 -*-
"""脱靶检测阶段（Stage 5）适配层。

两层并行、可分别开关（configs: stages.offtarget）：
  * PITA/TargetScan 种子区脱靶（offtarget/adapter.py）；
  * BLAST 近全长同源脱靶（offtarget/blast_adapter.py，默认关）。
"""
from .adapter import OffTargetAdapter, run_offtarget  # noqa: F401
from .blast_adapter import (  # noqa: F401
    BlastOffTargetAdapter,
    run_blast_offtarget,
)
