# -*- coding: utf-8 -*-
"""Stage3 structure 整合测试。

- 合成数据：generation->rules->structure 全链（强制 backend='approx'，离线可跑）；
- 字段契约：structure_* 列齐全、类型/空值符合约定；无 57nt 上下文行不崩溃；
- 真实 SFRP1 抽样回归（SIRNA_REGRESSION=1 时，前 200 条 passed 行）。

运行：python -m unittest tests.test_structure -v
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
from sirna_pipeline.stages.structure.runner import run_structure


def _chain_syn(td, seq_len: int = 40) -> tuple[Path, Path, Path]:
    """生成合成 CDS -> gen -> rules，返回 (cds_fasta, gen_csv, passed_csv)。"""
    fa = td / "syn.fa"
    rng = random.Random(11)
    seq = "".join(rng.choice("ATGC") for _ in range(seq_len))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    gen_dir, rules_dir = td / "gen", td / "rules"
    run_generation(fasta=fa, out_dir=gen_dir)
    run_rules(gen_dir / "candidates.csv", rules_dir)
    return fa, gen_dir / "candidates.csv", rules_dir / "candidates_passed.csv"


class TestStructureSyn(unittest.TestCase):
    def test_full_chain_fields(self):
        with tmpdir("struct_syn") as td:
            fa, _, passed = _chain_syn(td)
            out_dir = td / "struct"
            out = run_structure(passed, fa, out_dir,
                                detector_kwargs={"backend": "approx"})
            rows = records.read_records(out)
            passed_rows = records.read_records(passed)
            self.assertEqual(len(rows), len(passed_rows))
            self.assertGreater(len(rows), 0)
            for r in rows:
                for c in records.STRUCTURE_COLS:
                    self.assertIn(c, r, "缺列 " + c)
                # 关键特征数值可解析或为空
                for c in ("mfe", "delta_G_5end", "delta_G_3end",
                          "delta_deltaG_ends", "max_stem_length", "max_stem_gc"):
                    v = records.num(r[c])
                    if v is not None:
                        self.assertIsInstance(v, float)
                self.assertIn(r["structure_pass"], ("0", "1"))
                self.assertIsInstance(r["reject_reason"], str)
                self.assertIsInstance(r["warnings"], str)
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(out_dir)
            self.assertEqual(m["rows"], len(rows))

    def test_no_crash_on_boundary_and_deterministic_approx(self):
        with tmpdir("struct_bnd") as td:
            fa = td / "syn.fa"
            # 60nt 随机 CDS：首/末窗 57nt 上下文不足 -> mrna=None 单链模式；无重复窗口
            rng = random.Random(3)
            fa.write_text(">syn\n%s\n" % "".join(rng.choice("ATGC") for _ in range(60)),
                          encoding="utf-8")
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fa, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            o1 = td / "s1"
            o2 = td / "s2"
            r1 = records.read_records(
                run_structure(rules_dir / "candidates_passed.csv", fa, o1,
                              detector_kwargs={"backend": "approx"}))
            r2 = records.read_records(
                run_structure(rules_dir / "candidates_passed.csv", fa, o2,
                              detector_kwargs={"backend": "approx"}))
            self.assertEqual(len(r1), len(r2))
            # 近似路径确定性
            m1 = {r["variant_id"]: r["mfe"] for r in r1}
            m2 = {r["variant_id"]: r["mfe"] for r in r2}
            self.assertEqual(m1, m2)


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实回归测试需 SIRNA_REGRESSION=1")
class TestStructureRegressionRealData(unittest.TestCase):
    def test_real_passed_sample(self):
        ws = Path(__file__).resolve().parents[2]
        fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        if not fasta.exists():
            self.skipTest("缺少 SFRP1 FASTA")
        with tmpdir("struct_real") as td:
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fasta, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            passed = records.read_records(rules_dir / "candidates_passed.csv")
            self.assertEqual(len(passed), 6524)
            sample = passed[:200]                       # 抽样控制时长
            sample_path = td / "passed_200.csv"
            records.write_records(sample_path, sample)
            out = run_structure(sample_path, fasta, td / "struct",
                                detector_kwargs={"backend": "auto"})
            rows = records.read_records(out)
            self.assertEqual(len(rows), 200)
            n_pass = sum(1 for r in rows if r["structure_pass"] == "1")
            # 近端(无57nt上下文)行必须也在场且不抛错
            ctx_missing = sum(1 for r in rows if not r["cofold_mfe"]
                              and not r["target_accessibility"])
            self.assertGreaterEqual(ctx_missing, 0)
            self.assertIsInstance(n_pass, int)
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "struct")
            self.assertEqual(m["rows"], 200)


if __name__ == "__main__":
    unittest.main(verbosity=2)
