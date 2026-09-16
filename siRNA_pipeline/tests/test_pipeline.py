# -*- coding: utf-8 -*-
"""Pipeline 端到端集成测试。

- 合成小 CDS：一条命令跑通 gen→rules→structure→thermo→offtarget(skip)→toxicity→rank；
- 真实 SFRP1 全链（SIRNA_REGRESSION=1 时；ViennaRNA 缺失环境自动降级并关闭结构硬过滤）。

运行：python -m unittest tests.test_pipeline -v
"""
import os
import random
import uuid
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

PROJ = Path(__file__).resolve().parent.parent
CFG_DIR = PROJ / "configs"

from sirna_pipeline.pipeline.orchestrator import build_config, run_pipeline  # noqa: E402
from sirna_pipeline.common import records  # noqa: E402


def _syn_fasta(td, length: int = 120) -> Path:
    fa = td / "syn.fa"
    rng = random.Random(19)
    seq = "".join(rng.choice("ATGC") for _ in range(length))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    return fa


class TestPipelineE2E(unittest.TestCase):
    def test_full_chain_synthetic(self):
        with tmpdir("pipe_syn") as td:
            fa = _syn_fasta(td)
            run_name = "ut_pipe_" + uuid.uuid4().hex[:8]
            cfg = build_config(CFG_DIR, run_name=run_name, fasta_override=fa)
            # 本机无 ViennaRNA：structure 用近似且不把近似结果当硬过滤
            cfg["stages"]["structure"]["detector_kwargs"] = {"backend": "approx"}
            cfg["stages"]["ranking"]["require_structure_pass"] = False
            cfg["run"]["run_name"] = run_name
            # 隔离：e2e 会把 rank/chemmod 产物**发布**到 outputs/results/（大赛交付目录），
            # 不隔离就会用合成序列覆盖真实交付物 → 这里把 outputs_root 指到临时目录。
            cfg["paths"]["outputs_root"] = str(td / "outputs")
            summary = run_pipeline(cfg)
            for k, s in summary["stages"].items():
                self.assertTrue(s["ok"], f"{k} 失败: {s}")
            outs = Path(cfg["paths"]["outputs_root"])
            # 每阶段产物在 runs/<run_name>/xx_*/candidates*.csv
            runs_dir = outs / "runs" / run_name
            self.assertTrue((runs_dir / "01_generation" / "candidates.csv").exists())
            self.assertTrue((runs_dir / "02_rules" / "candidates_passed.csv").exists())
            self.assertTrue((runs_dir / "03_structure" / "candidates_structure.csv").exists())
            self.assertTrue((runs_dir / "04_thermo" / "candidates_thermo.csv").exists())
            self.assertTrue((runs_dir / "07_ranking" / "candidates_ranked.csv").exists())
            # 结果发布
            res = outs / "results"
            self.assertTrue((res / "rank_final.csv").exists())
            self.assertTrue((res / f"run_{run_name}_summary.json").exists())
            rows = records.read_records(res / "rank_final.csv")
            ranked = [r for r in rows if r["final_rank"]]
            self.assertGreater(len(ranked), 0)
            self.assertGreaterEqual(summary["rank_meta"]["eligible"], 1)


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实全链回归需 SIRNA_REGRESSION=1")
class TestPipelineRealE2E(unittest.TestCase):
    def test_real_sfrp1_end_to_end(self):
        run_name = "e2e_real_" + uuid.uuid4().hex[:6]
        cfg = build_config(CFG_DIR, run_name=run_name)
        cfg["stages"]["ranking"]["require_structure_pass"] = False   # 无 ViennaRNA 环境
        # 隔离：同 test_full_chain_synthetic —— 不得覆盖 outputs/results/ 交付目录
        with tmpdir("pipe_real") as td:
            cfg["paths"]["outputs_root"] = str(td / "outputs")
            summary = run_pipeline(cfg)
            outs = Path(cfg["paths"]["outputs_root"]) / "results"
            rows = records.read_records(outs / "rank_final.csv")
            ranked = [r for r in rows if r["final_rank"]]
            self.assertGreater(len(ranked), 100)  # 真实数据应有较多可排序候选
        for k, s in summary["stages"].items():
            self.assertTrue(s["ok"], f"{k} 失败: {s}")
        self.assertEqual(summary["stages"]["generation"]["ok"], True)


if __name__ == "__main__":
    unittest.main(verbosity=2)
