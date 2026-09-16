# -*- coding: utf-8 -*-
"""Stage4 thermo 整合测试。

- 合成数据：generation->rules->thermo 全链（需 torch，缺失则跳过）；
- 不变式：WT 无错配(n_mismatch=0/惩罚≈0)，mut 单错配(n_mismatch=1)，ΔG_total<0，Tm>0；
- dG_mismatch_by_position 为 len-19 JSON 串；边界无上下文行全空；
- 真实 SFRP1 抽样回归（SIRNA_REGRESSION=1）。

运行：python -m unittest tests.test_thermo -v
"""
import json
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

try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except Exception:  # pragma: no cover
    _HAS_TORCH = False

from sirna_pipeline.stages.thermo.runner import run_thermo


def _chain_syn(td) -> tuple[Path, Path]:
    fa = td / "syn.fa"
    rng = random.Random(5)
    seq = "".join(rng.choice("ATGC") for _ in range(60))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    gen_dir, rules_dir = td / "gen", td / "rules"
    run_generation(fasta=fa, out_dir=gen_dir)
    run_rules(gen_dir / "candidates.csv", rules_dir)
    return fa, rules_dir / "candidates_passed.csv"


@unittest.skipUnless(_HAS_TORCH, "需要 torch")
class TestThermoSyn(unittest.TestCase):
    def test_full_chain_invariants(self):
        with tmpdir("thermo_syn") as td:
            fa, passed = _chain_syn(td)
            out = run_thermo(passed, fa, td / "thermo")
            rows = records.read_records(out)
            passed_rows = records.read_records(passed)
            self.assertEqual(len(rows), len(passed_rows))
            n_mut = 0
            for r in rows:
                for c in records.THERMO_COLS:
                    self.assertIn(c, r, "缺列 " + c)
                ctx = r.get("dG_mismatch_by_position") != ""
                if not ctx:
                    continue                     # 边界行：特征空，已在 manifest 计数
                pos = json.loads(r["dG_mismatch_by_position"])
                self.assertEqual(len(pos), 19)
                dg = records.num(r["dG_total"])
                self.assertIsNotNone(dg)
                self.assertLess(dg, 0.0)
                tm = records.num(r["Tm"])
                self.assertIsNotNone(tm)
                self.assertGreater(tm, 0.0)
                if r["kind"] == "wt":
                    self.assertEqual(records.num(r["n_mismatch"]), 0.0)
                else:
                    n_mut += 1
                    self.assertEqual(records.num(r["n_mismatch"]), 1.0)
                    mm = records.num(r["dG_mismatch_total"])
                    self.assertGreater(mm, 0.0)
            self.assertGreater(n_mut, 0)
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "thermo")
            self.assertEqual(m["rows"], len(rows))

    def test_thermo_equivalence_with_legacy_direct(self):
        """runner 的逐行结果 == 直接调用 legacy calculate()（同上下文）。"""
        with tmpdir("thermo_eq") as td:
            fa, passed = _chain_syn(td)
            run_thermo(passed, fa, td / "thermo")
            from sirna_pipeline.common import seqio, duplex
            from sirna_pipeline.stages.thermo.runner import _load_legacy
            legacy = _load_legacy()
            calc = legacy.siRNAThermoCalculator()
            _, cds = seqio.load_fasta_rna(fa)
            rows = records.read_records(td / "thermo" / "candidates_thermo.csv")
            for r in rows[:30]:
                m57 = duplex.mrna57_from_cds(cds, int(r["cds_start"]))
                if m57 is None:
                    continue
                ref = calc.calculate(r["guide_checked"], m57)
                self.assertAlmostEqual(float(records.num(r["dG_total"])),
                                       float(ref["dG_total"]), places=6)


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1"
                     and _HAS_TORCH, "真实回归需 SIRNA_REGRESSION=1 且 torch")
class TestThermoRegressionRealData(unittest.TestCase):
    def test_real_passed_sample(self):
        ws = Path(__file__).resolve().parents[2]
        fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        with tmpdir("thermo_real") as td:
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fasta, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            passed = records.read_records(rules_dir / "candidates_passed.csv")
            sample_path = td / "passed_200.csv"
            records.write_records(sample_path, passed[:200])
            out = run_thermo(sample_path, fasta, td / "thermo")
            rows = records.read_records(out)
            self.assertEqual(len(rows), 200)
            n_dg = sum(1 for r in rows if records.num(r["dG_total"]) is not None)
            # 200 条中绝大多数应有 57nt 上下文
            self.assertGreater(n_dg, 180)


if __name__ == "__main__":
    unittest.main(verbosity=2)
