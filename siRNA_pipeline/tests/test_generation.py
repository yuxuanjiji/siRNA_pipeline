# -*- coding: utf-8 -*-
"""Stage1 generation 整合测试。

- 合成小 CDS 上的确定性功能/不变式测试（快，总是运行）；
- 真实 SFRP1 回归测试（与旧目录成品 CSV 对拍），需设置环境变量 SIRNA_REGRESSION=1。

运行：python -m unittest tests.test_generation -v
"""
import os
import re
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover - 直接以文件方式运行时
    from _util import tmpdir

from sirna_pipeline.common import seqio, records
from sirna_pipeline.stages.generation.runner import run_generation

RE_VARIANT = re.compile(r"^(W\d+)_(g(?:1|12|17|18|19)):([AUCG])>([AUCG])$")

SITES = {"g1": 1, "g12": 12, "g17": 17, "g18": 18, "g19": 19}


def _syn_fasta(tmp) -> str:
    fa = tmp / "syn_SFRP1.fa"
    import random
    rng = random.Random(7)                 # 固定种子：确定性且无周期性重复窗口
    seq = "".join(rng.choice("ATGC") for _ in range(40))
    fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
    return str(fa)


class TestGenerationSyn(unittest.TestCase):
    def test_counts_and_invariants(self):
        with tmpdir("gen_syn") as td:
            out = td / "out1"
            run_generation(fasta=_syn_fasta(td), out_dir=out)
            rows = records.read_records(out / "candidates.csv")
            self.assertEqual(len(rows), 22 * 16)
            from collections import Counter
            kinds = Counter(r["kind"] for r in rows)
            self.assertEqual(kinds["wt"], 22)
            self.assertEqual(kinds["mut"], 22 * 15)
            records.validate_base(rows, "generation")
            # 每窗口 1wt+15mut
            per_win = Counter(r["window_id"] for r in rows)
            self.assertTrue(all(v == 16 for v in per_win.values()))
            # 变体属性
            muts = [r for r in rows if r["kind"] == "mut"]
            self.assertEqual(len(muts), 22 * 15)
            for r in muts:
                g = r["guide_checked"]
                self.assertEqual(len(g), 19)
                wt = r["guide_wt_19"]
                diff = [i for i, (a, b) in enumerate(zip(wt, g)) if a != b]
                self.assertEqual(len(diff), 1, r["variant_id"])
                m = RE_VARIANT.match(r["variant_id"])
                self.assertIsNotNone(m, r["variant_id"])
                self.assertEqual(SITES[m.group(2)], int(r["position"]))
                self.assertEqual(r["site"], m.group(2))
                self.assertEqual(wt[diff[0]], m.group(3))
                self.assertEqual(g[diff[0]], m.group(4))
                self.assertIn(r["pair_type"], ("GU", "MM"))
                # 配对标示
                target = r["target_mRNA_19"]
                pos = int(r["position"])
                self.assertEqual(r["paired_mRNA_nt"], target[19 - pos])
            # wt 校验：guide == rc(target)
            for r in rows:
                if r["kind"] == "wt":
                    self.assertEqual(seqio.rc_rna(r["target_mRNA_19"]),
                                     r["guide_checked"])

    def test_manifest_written(self):
        with tmpdir("gen_man") as td:
            out = td / "out2"
            run_generation(fasta=_syn_fasta(td), out_dir=out)
            from sirna_pipeline.common.stage_io import read_manifest
            m = read_manifest(out)
            self.assertEqual(m["stage"], "01_generation")
            self.assertEqual(m["windows"], 22)
            self.assertEqual(m["wt"], 22)
            self.assertEqual(m["mut"], 22 * 15)


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实回归测试需 SIRNA_REGRESSION=1")
class TestGenerationRegressionRealData(unittest.TestCase):
    """与旧目录成品（任务5/6 输出 CSV）对拍，证明新管道产物与旧逻辑一致。"""

    def setUp(self):
        ws = Path(__file__).resolve().parents[2]  # 工作区根(生科挑战赛)
        # _util.SRC = siRNA_pipeline/src -> parent.parent = 生科挑战赛
        self.fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        self.old_dir = ws / "项目搭建" / "候选序列生成"
        if not (self.fasta.exists() and (self.old_dir / "SFRP1_task5_sliding_window_19nt_fullcomp.csv").exists()):
            self.skipTest("缺少原始 SFRP1 数据/旧成品，跳过回归")

    def test_match_original_outputs(self):
        from sirna_pipeline.common.seqio import read_csv_rows
        with tmpdir("gen_real") as td:
            out = td / "out"
            run_generation(fasta=self.fasta, out_dir=out)
            rows = records.read_records(out / "candidates.csv")
            old5 = read_csv_rows(self.old_dir / "SFRP1_task5_sliding_window_19nt_fullcomp.csv")
            old6 = read_csv_rows(self.old_dir / "SFRP1_task6_site_mutagenesis_library.csv")
            self.assertEqual(len(old5), 927)
            self.assertEqual(len(old6), 927 * 15)
            self.assertEqual(len(rows), 927 + 927 * 15)
            # 抽样逐字段对拍 wt 行
            wt_map = {r["window_id"]: r for r in rows if r["kind"] == "wt"}
            for i in (0, 99, 500, 926):
                o = old5[i]
                r = wt_map[o["window_id"]]
                self.assertEqual(r["guide_checked"], o["guide_antisense_19"])
                self.assertEqual(r["target_mRNA_19"], o["target_mRNA_19"])
                self.assertEqual(r["cds_start"], o["cds_start"])
                self.assertEqual(r["nm_003012_start"], o["nm_003012_start"])
            # 抽样对拍 mut 行（按 variant_id 索引）
            old_map = {r["variant_id"]: r for r in old6}
            mut_rows = [r for r in rows if r["kind"] == "mut"]
            for i in (0, 1234, 6789, 13904):
                o = old6[i]
                r = mut_rows[i]
                self.assertEqual(r["variant_id"], o["variant_id"])
                self.assertEqual(r["guide_checked"], o["guide_mut_19"])
                self.assertEqual(r["guide_wt_19"], o["guide_wt_19"])
                self.assertEqual(r["site"], o["site"])
                self.assertEqual(r["position"], o["position"])
                self.assertEqual(r["pair_type"], o["pair_type"])
                self.assertEqual(r["paired_mRNA_nt"], o["paired_mRNA_nt"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
