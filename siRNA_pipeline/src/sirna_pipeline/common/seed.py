# -*- coding: utf-8 -*-
"""固定随机种子（可复现性；《附件5-代码提交要求》六·规范与稳健）。

主管道本身是确定性计算（滑窗/规则/结构近似/NN 热力学/排序均无随机抽样），
此处统一给 numpy/torch/python 打种子，防止任何可复现性问题由环境噪声引入，
并在 manifest / run_*_summary.json 记录所用种子。
"""
from __future__ import annotations

import os
import random


def seed_everything(seed: int = 42) -> None:
    """把 python / numpy / torch 的随机源统一置为 seed。

    numpy / torch 缺失时静默跳过（主管道不依赖它们产生随机性）。
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    try:
        import numpy as np  # type: ignore
        np.random.seed(seed)
    except Exception:  # noqa: BLE001
        pass
    try:
        import torch  # type: ignore
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:  # noqa: BLE001
        pass
