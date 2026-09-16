# -*- coding: utf-8 -*-
"""测试公共工具（临时目录 + 路径引导）。"""
import shutil
import sys
import uuid
from contextlib import contextmanager
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@contextmanager
def tmpdir(prefix="t"):
    """测试临时目录：建于工程 outputs/runs/ 下（沙箱允许 makedirs、拒绝 tempfile 0o700 属性操作）。"""
    base = Path(__file__).resolve().parent.parent / "outputs" / "runs"
    d = base / ("ut_%s_%s" % (prefix, uuid.uuid4().hex[:10]))
    d.mkdir(parents=True, exist_ok=False)
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)
