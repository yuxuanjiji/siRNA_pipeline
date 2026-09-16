# -*- coding: utf-8 -*-
"""Stage5 BLAST 近全长同源脱靶适配层整合测试。

- disabled 默认路径：offtarget_blast_available=0、判据列留空、manifest skip；
- enabled 但 blastn/库未就绪：状态原因明确、graceful skip 不抛；
- 桩脚本契约：模拟 BLAST/blast_offtarget.py 的 CLI 输出（flagged/clean/hits/stats）
  -> offtarget_blast_flag/n_hits/available 与共享 offtarget_risk 正确回填。

运行：python -m unittest tests.test_blast_offtarget -v
"""
import sys
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from sirna_pipeline.common import records
from sirna_pipeline.stages.offtarget.blast_adapter import (
    BlastOffTargetAdapter,
    run_blast_offtarget,
)


def _craft(td) -> Path:
    rows = []
    for wid in ("W0001", "W0002", "W0003"):
        g = {"W0001": "CCCUUCUACGUAACCUAAA",
             "W0002": "AUAUAUAUAUAUAUAUAUA",
             "W0003": "UGAAACAGCAAAUUACCCC"}[wid]
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


STUB_SCRIPT = r'''
# -*- coding: utf-8 -*-
# 桩脚本：复刻 blast_offtarget.py 的 CLI 输出契约（不真正跑 BLAST）。
import argparse, csv, json, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input"); ap.add_argument("--column", default="siRNA")
    ap.add_argument("--out-dir"); ap.add_argument("--prefix", default="blast")
    ap.add_argument("--target-gene", default="SFRP1")
    ap.add_argument("--min-paired", type=int, default=17)
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    a, _ = ap.parse_known_args()

    with open(a.input, encoding="utf-8-sig", newline="") as fh:
        rows = [dict(r) for r in csv.DictReader(fh)]

    flagged, clean = [], []
    for r in rows:
        # 预置判定：W0001_wt 判脱靶，其余通过
        (flagged if r["variant_id"] == "W0001_wt" else clean).append(r)

    def _write(fn, items):
        with open(os.path.join(a.out_dir, fn), "w", encoding="utf-8", newline="") as fh:
            if not items:
                fh.write("variant_id,siRNA\n")
                return
            w = csv.DictWriter(fh, fieldnames=["variant_id", "siRNA"])
            w.writeheader(); w.writerows(items)

    _write(a.prefix + "_flagged.csv", flagged)
    _write(a.prefix + "_clean.csv", clean)
    with open(os.path.join(a.out_dir, a.prefix + "_hits_all.tsv"), "w",
              encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["qseqid"])
        for i in (0, 2):          # W0001(0) 有 2 条命中、W0003(2) 有 1 条
            w.writerow(["q%d" % i])
            if i == 0:
                w.writerow(["q0"])
    stats = {"n_flagged_offtarget": len(flagged), "n_hits_all_ge16": 3,
             "criterion": "stub tolerant",
             "top_offtarget_genes": {"SFRP2": 2}}
    with open(os.path.join(a.out_dir, a.prefix + "_stats.json"), "w",
              encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False)
    print("stub done", len(rows))

main()
'''


def _ready_cfg(td: Path, script: Path) -> dict:
    """就绪配置：桩脚本 + 假 blastn + 假库（只过 status 探测，不真正执行 BLAST）。"""
    (td / "blastn.exe").write_text("", encoding="utf-8")
    (td / "db.nsq").write_text("", encoding="utf-8")
    return {
        "enabled": True, "script": str(script),
        "blast_bin": str(td), "db": str(td / "db"),
        "target_gene": "SFRP1", "min_paired": 17, "max_mismatch": 1,
        "threads": 1,
    }


class TestBlastOfftargetAdapter(unittest.TestCase):
    def test_disabled_graceful_skip(self):
        with tmpdir("bl_disable") as td:
            base = _craft(td)
            out_csv, status = run_blast_offtarget(base, td / "out")
            rows = records.read_records(out_csv)
            self.assertEqual(len(rows), 3)
            for r in rows:
                self.assertEqual(r["offtarget_blast_available"], "0")
                self.assertEqual(r.get("offtarget_blast_flag", ""), "")
            self.assertFalse(status["ok"])
            self.assertIn("disabled", status["reason"])
            from sirna_pipeline.common.stage_io import read_manifest
            self.assertEqual(read_manifest(td / "out")["mode"], "skip")

    def test_enabled_env_missing_graceful_skip(self):
        with tmpdir("bl_env") as td:
            ad = BlastOffTargetAdapter({"enabled": True})   # script 默认存在但无 bin/db
            st = ad.status()
            self.assertFalse(st["ok"])
            self.assertTrue(ad.cfg["script"] and Path(ad.cfg["script"]).exists())

    def test_stub_executor_contract(self):
        with tmpdir("bl_ext") as td:
            script = td / "stub_blast.py"
            script.write_text(STUB_SCRIPT, encoding="utf-8")
            base = _craft(td)
            out_csv, status = run_blast_offtarget(
                base, td / "out", cfg=_ready_cfg(td, script))
            self.assertTrue(status["ok"])
            rows = {r["variant_id"]: r for r in records.read_records(out_csv)}
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows["W0001_wt"]["offtarget_blast_flag"], "1")
            self.assertEqual(rows["W0002_wt"]["offtarget_blast_flag"], "0")
            self.assertEqual(rows["W0003_wt"]["offtarget_blast_flag"], "0")
            for r in rows.values():
                self.assertEqual(r["offtarget_blast_available"], "1")
            # 共享风险列：命中 → 1，通过 → 0
            self.assertEqual(rows["W0001_wt"]["offtarget_risk"], "1")
            self.assertEqual(rows["W0002_wt"]["offtarget_risk"], "0")
            # 命中数按 q 序回填（W0001 2 条、W0003 1 条）
            self.assertEqual(rows["W0001_wt"]["offtarget_blast_n_hits_ge16"], "2")
            self.assertEqual(rows["W0003_wt"]["offtarget_blast_n_hits_ge16"], "1")
            self.assertEqual(rows["W0002_wt"]["offtarget_blast_n_hits_ge16"], "0")
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(td / "out")
            self.assertEqual(m["mode"], "external")
            self.assertEqual(m["flagged"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
