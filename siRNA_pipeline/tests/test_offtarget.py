# -*- coding: utf-8 -*-
"""Stage5 offtarget adapter 整合测试。

- disabled 默认路径：seed6 填充、offtarget_available=0、预测列留空、manifest skip；
- enabled 缺 external_cmd：状态原因明确、graceful skip 不抛；
- 外部执行器契约：mock 脚本子进程 -> 列正确合并（mock 值仅用于验证管道，无科学含义）。

运行：python -m unittest tests.test_offtarget -v
"""
import argparse
import os
import sys
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.offtarget.adapter import run_offtarget, OffTargetAdapter


def _craft(td) -> Path:
    rows = []
    for i, wid in enumerate(("W0001", "W0002", "W0003")):
        g = ["CCCUUCUACGUAACCUAAA", "AUAUAUAUAUAUAUAUAUA",
             "UGAAACAGCAAAUUACCCC"][i]
        rows.append({
            "window_id": wid, "variant_id": wid + "_wt", "kind": "wt",
            "cds_start": "20", "cds_end": "38",
            "target_mRNA_19": "ACGUACGUACGUACGUACG",
            "guide_wt_19": g, "guide_checked": g, "sense_strand_19": "",
            "rules_pass": "1",
        })
    p = td / "base.csv"
    records.write_records(p, rows)
    return p


MOCK_EXECUTOR = r'''
import argparse, csv, sys
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    with open(a.input, encoding="utf-8-sig", newline="") as fh:
        rows = [dict(r) for r in csv.DictReader(fh)]
    out = []
    for r in rows:
        # 仅供管道契约测试；真实打分须由 PITA/TargetScan 工具链提供
        out.append({
            "variant_id": r["variant_id"],
            "offtarget_pita_score": "0.02",
            "offtarget_pita_filter": "0",
            "offtarget_targetscan_score": "-0.5",
            "offtarget_targetscan_filter": "0",
            "offtarget_pass": "1",
        })
    with open(a.output, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
main()
'''


class TestOfftargetAdapter(unittest.TestCase):
    def test_disabled_graceful_skip(self):
        with tmpdir("oft_dis") as td:
            base = _craft(td)
            out_csv, status = run_offtarget(base, td / "out")
            rows = records.read_records(out_csv)
            self.assertEqual(len(rows), 3)
            for r in rows:
                self.assertEqual(r["seed6"], r["guide_checked"][1:7])
                self.assertEqual(r["offtarget_available"], "0")
                self.assertEqual(r.get("offtarget_pass", ""), "")
            self.assertFalse(status["ok"])
            self.assertIn("disabled", status["reason"])
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "out")
            self.assertEqual(m["mode"], "skip")

    def test_enabled_without_executor_reports_reason(self):
        ad = OffTargetAdapter({"enabled": True, "external_cmd": None})
        st = ad.status()
        self.assertFalse(st["ok"])
        self.assertIn("external_cmd", st["reason"])

    def test_external_executor_contract(self):
        with tmpdir("oft_ext") as td:
            script = td / "mock_executor.py"
            script.write_text(MOCK_EXECUTOR, encoding="utf-8")
            base = _craft(td)
            cfg = {"enabled": True, "external_cmd": str(script),
                   "python": sys.executable}
            out_csv, status = run_offtarget(base, td / "out", cfg=cfg)
            self.assertTrue(status["ok"])
            rows = records.read_records(out_csv)
            self.assertEqual(len(rows), 3)
            for r in rows:
                self.assertEqual(r["offtarget_available"], "1")
                self.assertEqual(r["offtarget_pass"], "1")
                self.assertEqual(r["offtarget_pita_score"], "0.02")
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "out")
            self.assertEqual(m["mode"], "external")
            self.assertEqual(m["mapped"], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
