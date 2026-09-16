"""统一训练入口：复现 C_match 变体层模型（自建错配数据集）与 OligoFormer head 校准。

用法：
  python train.py --model cmatch            # C_match 错配代价矩阵（Ridge 校准）
  python train.py --model siamese           # 变体层 WT/mismatch Siamese CNN
  python train.py --model oligo_head        # OligoFormer 冻结表征 + 小头 LOSO 校准（需 OligoFormer 资产就绪）

各脚本的额外参数（--input / --output 等）直接透传给对应 scripts/ 脚本。
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"

MODELS = {
    "cmatch": {
        "script": "fit_mismatch_ridge_artifact.py",
        "help": "C_match 错配代价矩阵（Ridge 校准；错配数据集，跨数据集 LOSO）",
        "needs_assets": False,
    },
    "siamese": {
        "script": "train_mismatch_siamese.py",
        "help": "变体层 WT/mismatch Siamese CNN（错配数据集，跨数据集 LOSO）",
        "needs_assets": False,
    },
    "oligo_head": {
        "script": "train_oligoformer_head_loso.py",
        "help": "OligoFormer 冻结表征 + 小头 LOSO 校准（需要 OligoFormer 资产：权重 / RNA-FM 数据）",
        "needs_assets": True,
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=sorted(MODELS), required=True, help="要复现的训练入口")
    args, unknown = parser.parse_known_args()

    entry = MODELS[args.model]
    script = SCRIPTS / entry["script"]
    if not script.exists():
        raise SystemExit(f"训练脚本缺失：{script}")

    command = [sys.executable, str(script), *unknown]
    print(f"[train] {entry['help']}")
    print(f"[train] 运行：{' '.join(command)}")
    if entry["needs_assets"]:
        print("[train] 提示：该入口需要 OligoFormer 资产就绪（权重 / RNA-FM 数据），缺失时脚本会自行报错。")
    subprocess.run(command, cwd=SCRIPTS, check=True)


if __name__ == "__main__":
    main()
