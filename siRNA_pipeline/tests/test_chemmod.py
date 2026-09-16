# -*- coding: utf-8 -*-
"""Stage8 chemmod（第六步 化学修饰）整合测试。

- 规则引擎：ESC-19 默认布点、免疫掩蔽把 U 位（含偶数 F 位）转 2′-OMe、长度校验；
- runner：取 Top-N、chem_* 列齐全、叶产物可复现。

运行：python -m unittest tests.test_chemmod -v
"""
import json
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.chemmod import rules as cm_rules
from sirna_pipeline.stages.chemmod.runner import run_chemmod


def _ranked_rows():
    rows = []
    for i, (wid, guide, rank) in enumerate([
            ("W0001", "CCCUUCUACGUAACCUAAA", 2),
            ("W0002", "AUAUAUAUAUAUAUAUAUA", 1),
            ("W0003", "UGAAACAGCAAAUUACCCC", 3)]):
        rows.append({
            "window_id": wid, "variant_id": wid + "_wt", "kind": "wt",
            "cds_start": "20", "cds_end": "38",
            "target_mRNA_19": "ACGUACGUACGUACGUACG",
            "guide_wt_19": guide, "guide_checked": guide, "sense_strand_19": "",
            "rules_pass": "1", "final_rank": str(rank),
        })
    return rows


class TestChemRules(unittest.TestCase):
    def test_default_esc19_plan(self):
        guide = "CCCUUCUACGUAACCUAAA"
        p = cm_rules.build_chem_plan(guide)
        self.assertEqual(len(p["guide"]), 19)
        self.assertEqual(len(p["passenger"]), 19)
        self.assertEqual(p["passenger"], cm_rules.rc(guide))
        # guide 偶数位 2′-F（2,4,...,18），且 g1/g19 不在 F 集
        self.assertEqual(p["guide_f_positions"],
                         [2, 4, 6, 8, 10, 12, 14, 16, 18])
        self.assertNotIn(1, p["guide_f_positions"])
        self.assertNotIn(19, p["guide_f_positions"])
        self.assertTrue(p["guide_5p_phosphate"])
        self.assertTrue(p["guide_mod"].startswith("5′-P "))

    def test_immune_mask_overrides_even_f(self):
        # 全 U：偶数位本应 2′-F，但整链为 poly-U → 全部 U 位掩蔽为 2′-OMe
        guide = "U" * 19
        p = cm_rules.build_chem_plan(guide)
        self.assertEqual(p["guide_f_positions"], [])
        joined = " ".join(p["notes"])
        self.assertIn("连续 U 串", joined)
        # 记法中不应出现 fU
        self.assertNotIn("fU", p["guide_mod"])

    def test_motif_positions(self):
        self.assertEqual(cm_rules.find_motif_positions("UGUGUAAA"),
                         {1, 2, 3, 4, 5})
        # poly-U：≥3 连续 U 才算（避免普通 UU 过度掩蔽）
        self.assertEqual(cm_rules.find_polyu_positions("UAUUUAU"), {3, 4, 5})

    def test_length_guard(self):
        with self.assertRaises(ValueError):
            cm_rules.build_chem_plan("ACGU")


class TestChemmodRunner(unittest.TestCase):
    def test_top_n_and_columns(self):
        with tmpdir("chem_run") as td:
            base = td / "ranked.csv"
            records.write_records(base, _ranked_rows())

            out_all = run_chemmod(base, td / "all", top_n=None)
            rows_all = records.read_records(out_all)
            self.assertEqual(len(rows_all), 3)
            # runner 按 final_rank 升序输出
            self.assertEqual([r["window_id"] for r in rows_all],
                             ["W0002", "W0001", "W0003"])
            for r in rows_all:
                self.assertTrue(r["chem_guide_mod"].startswith("5′-P "))
                self.assertEqual(r["chem_5p_phosphate"], "1")
                self.assertTrue(r["chem_sense_mod"])
                notes = json.loads(r["chem_notes"])
                self.assertIsInstance(notes, list)

            out_1 = run_chemmod(base, td / "top1", top_n=1)
            rows_1 = records.read_records(out_1)
            self.assertEqual(len(rows_1), 1)
            self.assertEqual(rows_1[0]["window_id"], "W0002")

            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "all")
            self.assertEqual(m["stage"], "08_chemmod")
            self.assertEqual(m["rows"], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
