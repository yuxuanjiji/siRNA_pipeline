# -*- coding: utf-8 -*-
"""Pipeline 唯一入口（薄封装，逻辑在 sirna_pipeline.cli）。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from sirna_pipeline.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
