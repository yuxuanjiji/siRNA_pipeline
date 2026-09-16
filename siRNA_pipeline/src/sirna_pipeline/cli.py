# -*- coding: utf-8 -*-
"""包内 CLI（console script 与 scripts/run_pipeline.py 共用）。"""
import argparse
import logging
from pathlib import Path

from .common.stage_io import setup_logging
from .pipeline.orchestrator import PipelineError, build_config, run_pipeline

_STAGES = ["generation", "rules", "structure", "thermo",
           "offtarget", "toxicity", "rank", "chemmod"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="siRNA 筛选与排序 Pipeline")
    ap.add_argument("--config-dir", type=Path, default=None,
                    help="配置文件目录（默认：包旁 configs/，即工程根 configs）")
    ap.add_argument("--run-name", default=None)
    ap.add_argument("--fasta", type=Path, default=None)
    ap.add_argument("--from", dest="stage_from", choices=_STAGES, default=None)
    ap.add_argument("--until", dest="stage_until", choices=_STAGES, default=None)
    ap.add_argument("--log-level", default="INFO")
    args = ap.parse_args(argv)

    # 默认 configs 位于工程根（src/sirna_pipeline/cli.py 向上三级 = 工程根）
    config_dir = args.config_dir
    if config_dir is None:
        here = Path(__file__).resolve()
        cand = here.parents[2] / "configs"
        if (cand / "pipeline.yaml").exists():
            config_dir = cand
    if config_dir is None:
        raise SystemExit("找不到 configs/（可用 --config-dir 显式指定）")

    log = setup_logging("sirna.cli", level=getattr(logging, args.log_level.upper(),
                                                    logging.INFO))
    try:
        cfg = build_config(config_dir, run_name=args.run_name,
                           fasta_override=args.fasta)
    except Exception as e:  # noqa: BLE001
        log.error("配置加载失败：%s", e)
        return 2
    try:
        summary = run_pipeline(cfg, log=log, stage_from=args.stage_from,
                               stage_until=args.stage_until)
    except PipelineError as e:
        log.error("pipeline 失败：%s", e)
        return 3
    ok = all(s.get("ok") for s in summary["stages"].values())
    log.info("pipeline 结束：ok=%s", ok)
    return 0 if ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
