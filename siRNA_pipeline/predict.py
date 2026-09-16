#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""一键主运行入口（大赛提交规范；与 scripts/run_pipeline.py 功能等效）。

读取 configs/ 配置 → 依序运行 6 步技术路线各阶段
（生成→规则→结构→热力/脱靶/毒性→综合排序→化学修饰），
自动产出大赛标准结果文件 outputs/results/results.csv（UTF-8）。

用法（在工程根目录 siRNA_pipeline/ 下）：
    python predict.py                                   # 真实 SFRP1 全链（默认）
    python predict.py --run-name myrun --until chemmod  # 自定义运行名
    python predict.py --fasta <其他CDS.fa>              # 换靶基因
更多参数见 `python predict.py --help`；完整命令与说明见 README.md。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from sirna_pipeline.cli import main  # noqa: E402

if __name__ == "__main__":
    # 主流程
    code = main()
    # 主流程成功 → 一键生成大赛标准化最终结果文件（results.xlsx + UTF-8 results.csv）
    if code == 0:
        root = Path(__file__).resolve().parent
        sys.path.insert(0, str(root))
        try:
            from export_final_results import export_final_results
            export_final_results(
                results_path=str(root / "outputs" / "results" / "results.csv"),
                rank_path=str(root / "outputs" / "results" / "rank_final.csv"),
                chemmod_path=str(root / "outputs" / "results" / "chemmod_top.csv"),
                outdir=str(root / "final_results"),
                topn=50,
            )
            print("[predict] 已生成大赛最终结果文件 → final_results/results.xlsx + results.csv")
        except Exception as e:  # noqa: BLE001
            print(f"[predict] 最终结果文件导出失败（不影响主流程产物）: {e}", file=sys.stderr)
    raise SystemExit(code)
