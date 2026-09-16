# -*- coding: utf-8 -*-
"""Stage2 rule_filter 整合测试。

- 合成数据端到端（generation -> rules）不变式；
- 与 legacy evaluate() 等价性（构造串覆盖四项规则）；
- 真实 SFRP1 回归（SIRNA_REGRESSION=1 时与旧 task7 成品对拍）。

运行：python -m unittest tests.test_rules -v
"""
import os
import unittest
from collections import Counter
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.generation.runner import run_generation
from sirna_pipeline.stages.rule_filter.runner import run_rules

# 与 legacy 内部模块等价的构造串（见 task7 selftests 注释）
SEQ_PASS = "CCCUUCUACGUAACCUAAA"          # 种子固定随机 19mer，四项全过
SEQ_GC_BAD = "GC" * 7 + "AUGCA"           # GC≈84% > 65% -> 规则1
SEQ_HOMO = "ACGUAAAAAUCGCGCAUAG"          # 连续5个A -> 规则3
SEQ_RUNGC = "GCGCGCGCGCGCGCGCGCG"         # 连续6个G/C -> 规则2
SEQ_PAL = "ACGUACGU" + "AAGCAUGC" + "CGC"  # 回文(4nt臂) -> 规则4


def _syn_generation(td) -> records:
    fa = td / "syn.fa"
    import random
    rng = random.Random(7)
    seq = "".join(rng.choice("ATGC") for _ in range(40))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    out = td / "gen"
    run_generation(fasta=fa, out_dir=out)
    return out / "candidates.csv"


class TestRulesSyn(unittest.TestCase):
    def test_end_to_end_columns_and_passed_subset(self):
        with tmpdir("rules_syn") as td:
            gen_csv = _syn_generation(td)
            full_csv, passed_csv = run_rules(gen_csv, td / "rules")
            full = records.read_records(full_csv)
            passed = records.read_records(passed_csv)
            self.assertEqual(len(full), 22 * 16)
            for r in full:
                for c in records.RULES_COLS:
                    self.assertIn(c, r)
            self.assertTrue(all(r["rules_pass"] == "1" for r in passed))
            self.assertEqual(len(passed), sum(1 for r in full if r["rules_pass"] == "1"))
            self.assertLess(len(passed), len(full))   # 四项规则必淘汰一部分
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "rules")
            self.assertEqual(m["passed"], len(passed))
            self.assertEqual(m["total"], len(full))

    def test_matches_legacy_evaluate_on_crafted(self):
        # 直接调用 legacy evaluate，与 runner._ev_row 语义一致（不重复造算法）
        from sirna_pipeline.stages.rule_filter.runner import _load_legacy
        t7 = _load_legacy()
        cases = [
            (SEQ_PASS, {"fail_gc_range": 0, "fail_run6_gc": 0,
                        "fail_run5_same": 0, "fail_palindrome": 0}, True),
            (SEQ_GC_BAD, {"fail_gc_range": 1}, False),
            (SEQ_HOMO, {"fail_run5_same": 1}, False),
            (SEQ_RUNGC, {"fail_run6_gc": 1}, False),
            (SEQ_PAL, {"fail_palindrome": 1}, False),
        ]
        for seq, expect_flags, expect_pass in cases:
            ev = t7.evaluate(seq)
            self.assertEqual(ev["seq_pass"], 1 if expect_pass else 0, seq)
            for k, v in expect_flags.items():
                self.assertEqual(ev[k], v, (seq, k))


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实回归测试需 SIRNA_REGRESSION=1")
class TestRulesRegressionRealData(unittest.TestCase):
    """旧 task7 成品（annotated/passed）对拍。"""

    def setUp(self):
        ws = Path(__file__).resolve().parents[2]
        self.old_dir = ws / "项目搭建" / "初级规则筛选"
        self.ann = self.old_dir / "SFRP1_task7_sequence_rules_annotated.csv"
        self.pas = self.old_dir / "SFRP1_task7_sequence_rules_passed.csv"
        if not (self.ann.exists() and self.pas.exists()):
            self.skipTest("缺少旧 task7 成品，跳过回归")

    def test_equals_original(self):
        from sirna_pipeline.common.seqio import read_csv_rows
        ws = Path(__file__).resolve().parents[2]
        fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        with tmpdir("rules_real") as td:
            gen_csv = td / "gen" / "candidates.csv"
            run_generation(fasta=fasta, out_dir=gen_csv.parent)
            full_csv, passed_csv = run_rules(gen_csv, td / "rules")
            full = records.read_records(full_csv)
            old = read_csv_rows(self.ann)
            old_pass = read_csv_rows(self.pas)
            # 全量行数一致
            self.assertEqual(len(full), len(old))
            # 通过集合一致（以 guide_checked 多重集比较，语义等价于旧 passed 表）
            got_pass = Counter(r["guide_checked"] for r in full if r["rules_pass"] == "1")
            old_pass_g = Counter(r["guide_checked"] for r in old_pass)
            self.assertEqual(got_pass, old_pass_g)
            self.assertEqual(sum(got_pass.values()), 6524)
            # 抽样逐行对拍判定（mut 行按 variant_id，wt 行按 window_id）
            old_by_key = {}
            for o in old:
                if o["variant_id"]:
                    old_by_key[("mut", o["variant_id"])] = o
                else:
                    old_by_key[("wt", o["window_id"])] = o
            for i in (0, 500, 4000, 8000, 12000, 14831):
                r = full[i]
                key = ("mut", r["variant_id"]) if r["kind"] == "mut" else ("wt", r["window_id"])
                o = old_by_key[key]
                self.assertEqual(r["rules_pass"], records.to_int01(o["seq_pass"]),
                                 (key, r["guide_checked"]))
                self.assertEqual(r["hit_code"], str(o["hit_code"]), key)
                self.assertEqual(r["hit_rule"], o["hit_rule"], key)


if __name__ == "__main__":
    unittest.main(verbosity=2)
