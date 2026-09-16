# -*- coding: utf-8 -*-
"""Stage6 toxicity 整合测试。

- 合成数据端到端：tox_*/imm_* 列契约、seed 缓存广播计数；
- cell_viability 表加载（4096 行）；
- motif 判定与 legacy 语义一致（UGUGU 阳性等）；
- 真实 SFRP1 抽样回归（SIRNA_REGRESSION=1）。

运行：python -m unittest tests.test_toxicity -v
"""
import os
import random
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.generation.runner import run_generation
from sirna_pipeline.stages.rule_filter.runner import run_rules
from sirna_pipeline.stages.toxicity.runner import run_toxicity


def _chain_syn(td) -> tuple[Path, Path]:
    fa = td / "syn.fa"
    rng = random.Random(6)
    seq = "".join(rng.choice("ATGC") for _ in range(60))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    gen_dir, rules_dir = td / "gen", td / "rules"
    run_generation(fasta=fa, out_dir=gen_dir)
    run_rules(gen_dir / "candidates.csv", rules_dir)
    return fa, rules_dir / "candidates_passed.csv"


class TestToxicitySyn(unittest.TestCase):
    def test_columns_and_seed_broadcast(self):
        with tmpdir("tox_syn") as td:
            _, passed = _chain_syn(td)
            out = run_toxicity(passed, td / "tox")
            rows = records.read_records(out)
            passed_rows = records.read_records(passed)
            self.assertEqual(len(rows), len(passed_rows))
            for r in rows:
                for c in records.TOXICITY_COLS:
                    self.assertIn(c, r, "缺列 " + c)
                # seed = guide[1:7]（g2-g7 六聚体）
                self.assertEqual(r["tox_seed"], r["guide_checked"][1:7])
                self.assertIn(r["tox_viability_flag"], ("0", "1", ""))
                self.assertIn(r["imm_high_flag"], ("0", "1"))
            n_bc = sum(1 for r in rows if r["tox_broadcast"] == "1")
            n_win = len({r["window_id"] for r in rows})
            # 同窗口 16 条共享 seed：仅每窗口首条 miss，其余广播
            self.assertEqual(n_bc, len(rows) - n_win)

    def test_viability_table_loaded_4096(self):
        from sirna_pipeline.stages.toxicity.runner import _load_legacy
        legacy = _load_legacy()
        det = legacy.ToxicityDetector()
        table = det._ensure_table()
        self.assertEqual(len(table), 4096)

    def test_motif_semantics(self):
        from sirna_pipeline.stages.toxicity.runner import _load_legacy
        legacy = _load_legacy()
        det = legacy.ToxicityDetector()
        # guide 含 UGUGU -> imm_ugugu True, imm_high_flag True
        g = "UGUGU" + "A" * 14          # 19nt
        r = det.detect(g)
        self.assertTrue(r["imm_ugugu"])
        self.assertTrue(r["imm_high_flag"])
        # 阴性对照（无 motif/polyU）
        g2 = "C" * 19
        r2 = det.detect(g2)
        self.assertFalse(r2["imm_high_flag"])


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实回归测试需 SIRNA_REGRESSION=1")
class TestToxicityRegressionRealData(unittest.TestCase):
    def test_real_passed_sample(self):
        ws = Path(__file__).resolve().parents[2]
        fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        with tmpdir("tox_real") as td:
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fasta, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            passed = records.read_records(rules_dir / "candidates_passed.csv")
            sample = passed[:200]
            sp = td / "passed_200.csv"
            records.write_records(sp, sample)
            out = run_toxicity(sp, td / "tox")
            rows = records.read_records(out)
            self.assertEqual(len(rows), 200)
            # 4096 hexamer 全覆盖 => 200 条都应命中表（无 miss）
            miss = sum(1 for r in rows if r["tox_viability_miss"] == "1")
            self.assertEqual(miss, 0)
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "tox")
            self.assertEqual(m["rows"], 200)


if __name__ == "__main__":
    unittest.main(verbosity=2)
