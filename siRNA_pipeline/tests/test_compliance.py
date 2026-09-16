# -*- coding: utf-8 -*-
"""大赛提交合规测试（对照《附件5-代码提交要求》）。

- 标准结果文件：UTF-8 无 BOM、必填字段齐全、按 final_rank 升序、JSON 参数列可被正常解析；
- 固定随机种子：同种子可复现、不同种子可区分。

运行：python -m unittest tests.test_compliance -v
"""
import csv
import json
import random
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.common.results_export import (
    HEADER,
    TRACK_DEFAULT,
    export_standard_results,
)
from sirna_pipeline.common.seed import seed_everything


def _ranked_csv(td: Path) -> Path:
    rows = []
    guides = ["CCCUUCUACGUAACCUAAA", "AUAUAUAUAUAUAUAUAUA",
              "UGAAACAGCAAAUUACCCC"]
    for i, (wid, g, rank) in enumerate(zip(
            ("W0001", "W0002", "W0003"), guides, ("3", "1", "2"))):
        rows.append({
            "window_id": wid, "variant_id": wid + "_wt", "kind": "wt",
            "cds_start": "20", "cds_end": "38",
            "target_mRNA_19": "ACGUACGUACGUACGUACG",
            "guide_wt_19": g, "guide_checked": g, "sense_strand_19": "",
            "rules_pass": "1", "structure_pass": "1",
            "mfe": "-3.0", "delta_deltaG_ends": "0.2", "dG_total": "-41.0",
            "feat_dG_seed": "-9.5",
            "score_thermo": "0.55", "penalty_total": "0.0",
            "final_score": "0.55", "final_rank": rank,
            "tox_viability_score": "90", "tox_viability_flag": "0",
            "imm_high_flag": "0", "warnings": "自折叠提示：某告警",
        })
    p = td / "ranked.csv"
    records.write_records(p, rows)
    return p


class TestStandardResultsFile(unittest.TestCase):
    def test_header_fields_present(self):
        for col in ("candidate_id", "track", "guide_seq", "final_score",
                    "final_rank", "pipeline_version", "ranking_params", "notes"):
            self.assertIn(col, HEADER)

    def test_export_utf8_no_bom_sorted_and_parseable(self):
        with tmpdir("sub_res") as td:
            rank_csv = _ranked_csv(td)
            out = export_standard_results(rank_csv, td / "results.csv", top_n=2)
            raw = out.read_bytes()
            self.assertFalse(raw.startswith(b"\xef\xbb\xbf"), "结果文件应为纯 UTF-8 无 BOM")
            text = raw.decode("utf-8")
            with open(out, encoding="utf-8", newline="") as fh:
                rd = list(csv.DictReader(fh))
            self.assertEqual(len(rd), 2)                     # top_n=2
            self.assertEqual([r["candidate_id"] for r in rd],
                             ["W0002_wt", "W0003_wt"])       # final_rank 升序 1,2
            self.assertEqual(rd[0]["track"], TRACK_DEFAULT)
            self.assertEqual(rd[0]["pipeline_version"], "sirna_pipeline@0.1.0")
            self.assertTrue(rd[0]["ranking_params"].startswith("{"))
            json.loads(rd[0]["ranking_params"])              # JSON 参数列可解析
            self.assertIn("自折叠提示", rd[0]["notes"])

    def test_export_all_ranked_when_no_top_n(self):
        with tmpdir("sub_all") as td:
            out = export_standard_results(_ranked_csv(td), td / "r.csv", top_n=None)
            with open(out, encoding="utf-8", newline="") as fh:
                self.assertEqual(sum(1 for _ in fh) - 1, 3)


    def test_export_with_chem_reference_columns(self):
        """给定 chemmod 产物时，结果文件追加参考化学修饰列并按候选 join。"""
        with tmpdir("sub_chem") as td:
            rank_csv = _ranked_csv(td)
            chem = td / "chem.csv"
            records.write_records(chem, [{
                "variant_id": "W0002_wt",          # _ranked_csv 中 final_rank=1 的候选
                "chem_guide_mod": "5′-P mA s fU s mC",
                "chem_sense_mod": "mU s mG s mC",
                "chem_notes": json.dumps(["ESC-19 骨架", "命中连续 U 串"], ensure_ascii=False),
            }])
            out = export_standard_results(
                rank_csv, td / "r.csv", top_n=3, chem_csv=chem,
                mod_rule="esc19-v1", seed_ome=False)
            with open(out, encoding="utf-8", newline="") as fh:
                header = next(csv.reader(fh))
            for col in ("guide_mod", "sense_mod", "mod_rule",
                        "mod_notes", "mod_seed_ome"):
                self.assertIn(col, header)
            with open(out, encoding="utf-8", newline="") as fh:
                rd = {r["candidate_id"]: r for r in csv.DictReader(fh)}
            top1 = rd["W0002_wt"]
            self.assertTrue(top1["guide_mod"].startswith("5′-P"))
            self.assertEqual(top1["sense_mod"], "mU s mG s mC")
            self.assertEqual(top1["mod_rule"], "esc19-v1")
            self.assertEqual(top1["mod_seed_ome"], "0")
            self.assertIn("ESC-19", top1["mod_notes"])
            # 未提供修饰的候选留空，不影响解析
            self.assertEqual(rd["W0003_wt"]["guide_mod"], "")


class TestSeed(unittest.TestCase):
    def test_same_seed_reproducible(self):
        seed_everything(7)
        a = [random.random() for _ in range(5)]
        seed_everything(7)
        b = [random.random() for _ in range(5)]
        self.assertEqual(a, b)

    def test_different_seed_differs(self):
        seed_everything(1)
        a = [random.random() for _ in range(5)]
        seed_everything(2)
        b = [random.random() for _ in range(5)]
        self.assertNotEqual(a, b)


if __name__ == "__main__":
    unittest.main(verbosity=2)
