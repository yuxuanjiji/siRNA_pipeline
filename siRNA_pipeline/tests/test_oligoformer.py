# -*- coding: utf-8 -*-
"""OligoFormer 效率适配层整合测试（可选阶段）。

- disabled 默认路径：oligo_available=0、efficacy 留空、manifest skip；
- enabled 缺 external_cmd：状态原因明确、graceful skip 不抛；
- 外部执行器契约：mock 脚本子进程 -> oligo_efficacy 正确回填（mock 值仅验证管道，无科学含义）。

运行：python -m unittest tests.test_oligoformer -v
"""
import sys
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.oligoformer.adapter import (
    OligoFormerAdapter,
    run_oligoformer,
)


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
import argparse, csv
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    with open(a.input, encoding="utf-8-sig", newline="") as fh:
        rows = [dict(r) for r in csv.DictReader(fh)]
    out = []
    for r in rows:
        # 仅供管道契约测试；真实打分须由 OligoFormer 推理环境提供
        out.append({"variant_id": r["variant_id"], "oligo_efficacy": "0.66"})
    with open(a.output, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["variant_id", "oligo_efficacy"])
        w.writeheader()
        w.writerows(out)
main()
'''


class TestOligoFormerAdapter(unittest.TestCase):
    def test_disabled_graceful_skip(self):
        with tmpdir("of_dis") as td:
            base = _craft(td)
            out_csv, status = run_oligoformer(base, td / "out")
            rows = records.read_records(out_csv)
            self.assertEqual(len(rows), 3)
            for r in rows:
                self.assertEqual(r["oligo_available"], "0")
                self.assertEqual(r.get("oligo_efficacy", ""), "")
            self.assertFalse(status["ok"])
            self.assertIn("disabled", status["reason"])
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "out")
            self.assertEqual(m["mode"], "skip")

    def test_enabled_missing_assets_reports_reason(self):
        # 用一个"只有模型、没有 RNA-FM"的临时仓库，确定性地验证缺资产时报因并 graceful skip
        with tmpdir("of_noassets") as td:
            repo = td / "repo"
            (repo / "model").mkdir(parents=True)
            (repo / "model" / "best_model.pth").write_text("", encoding="utf-8")
            ad = OligoFormerAdapter({
                "enabled": True, "repo_dir": str(repo),
                "model_path": str(repo / "model" / "best_model.pth"),
                "python": sys.executable,
            })
            st = ad.status()
            self.assertFalse(st["ok"])
            self.assertIn("RNA-FM", st["reason"])

    def test_enabled_missing_external_reports_reason(self):
        ad = OligoFormerAdapter({"enabled": True,
                                 "external_cmd": "C:/no_such_executor_oligo.exe"})
        st = ad.status()
        self.assertFalse(st["ok"])
        self.assertIn("external_cmd", st["reason"])

    def test_external_executor_contract(self):
        with tmpdir("of_ext") as td:
            script = td / "mock_executor.py"
            script.write_text(MOCK_EXECUTOR, encoding="utf-8")
            base = _craft(td)
            cfg = {"enabled": True, "external_cmd": str(script),
                   "python": sys.executable}
            out_csv, status = run_oligoformer(base, td / "out", cfg=cfg)
            self.assertTrue(status["ok"])
            rows = records.read_records(out_csv)
            self.assertEqual(len(rows), 3)
            for r in rows:
                self.assertEqual(r["oligo_available"], "1")
                self.assertEqual(r["oligo_efficacy"], "0.66")
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "out")
            self.assertEqual(m["mode"], "external")
            self.assertEqual(m["mapped"], 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
