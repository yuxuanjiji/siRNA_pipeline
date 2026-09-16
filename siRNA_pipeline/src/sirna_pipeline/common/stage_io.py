# -*- coding: utf-8 -*-
"""阶段 IO 公共设施：配置加载 / 路径解析 / manifest / 日志。

- 配置文件用 YAML（pyyaml）；若未安装 pyyaml，报错并提示安装（见 requirements.txt）。
- 所有相对路径以“配置文件所在目录”为基准解析，杜绝硬编码绝对路径。
"""
from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover - 环境缺失提示
    yaml = None


class ConfigError(Exception):
    """配置错误（退出码 2 语义）。"""


def load_yaml(path: str | Path) -> dict:
    if yaml is None:
        raise ConfigError(
            "缺少 pyyaml：请先安装（python -m pip install pyyaml）或改用 JSON 配置。")
    with open(path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ConfigError(f"配置顶层应为映射：{path}")
    return data


def load_config(path: str | Path) -> dict:
    """统一配置加载：.json -> json；.yaml/.yml -> yaml；其余抛 ConfigError。"""
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"配置文件不存在：{p}")
    if p.suffix.lower() in (".yaml", ".yml"):
        return load_yaml(p)
    if p.suffix.lower() == ".json":
        with open(p, "r", encoding="utf-8") as fh:
            return json.load(fh)
    raise ConfigError(f"不支持的配置后缀：{p.suffix}（支持 .yaml/.yml/.json）")


def resolve_paths(paths_cfg: dict, base_dir: str | Path) -> dict:
    """把 paths.yaml 里相对路径解析为绝对路径（相对于配置文件所在目录）。

    支持 ${ENV_VAR} 占位与环境变量覆盖（同名大写键优先，如 SIRNA_DATA_ROOT）。
    """
    import re as _re
    base = Path(base_dir).resolve()
    out = {}
    for key, raw in paths_cfg.items():
        if not isinstance(raw, str):
            out[key] = raw
            continue
        s = raw.strip()
        # 环境变量覆盖（若存在同名大写环境变量）
        env_val = os.environ.get(key.upper())
        if env_val:
            s = env_val
        # ${VAR} 占位
        s = _re.sub(r"\$\{(\w+)\}", lambda m: os.environ.get(m.group(1), ""), s)
        p = Path(s)
        out[key] = p.resolve() if p.is_absolute() else (base / p).resolve()
    return out


# ---------------------------------------------------------------------------
# manifest（每阶段产物附带 json：输入哈希/参数/统计）
# ---------------------------------------------------------------------------
def write_manifest(out_dir: str | Path, meta: dict) -> Path:
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    mp = Path(out_dir) / "manifest.json"
    with open(mp, "w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)
    return mp


def read_manifest(out_dir: str | Path) -> dict | None:
    mp = Path(out_dir) / "manifest.json"
    if not mp.exists():
        return None
    with open(mp, "r", encoding="utf-8") as fh:
        return json.load(fh)


def file_sha256(path: str | Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# 日志：统一输出 outputs/logs/（也可由调用方传入文件路径）
# ---------------------------------------------------------------------------
def setup_logging(name: str = "sirna_pipeline",
                  log_file: str | Path | None = None,
                  level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)
    for h in list(logger.handlers):          # 先关闭旧句柄，避免文件句柄泄漏
        try:
            h.close()
        except Exception:  # noqa: BLE001
            pass
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if log_file is not None:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    logger.propagate = False
    return logger


def stage_logger(stage: str, log_dir: str | Path | None = None) -> logging.Logger:
    """阶段日志：默认到 logs/<stage>.log（若给 log_dir）。"""
    if log_dir is not None:
        return setup_logging(f"sirna.{stage}", Path(log_dir) / f"stage_{stage}.log")
    return setup_logging(f"sirna.{stage}")
