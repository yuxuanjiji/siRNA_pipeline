# -*- coding: utf-8 -*-
"""common 共享层单元测试（unittest，无第三方依赖）。

运行（工程根目录）：
    python -m unittest discover -s tests -p "test_*.py" -v
"""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

from contextlib import contextmanager


@contextmanager
def tmpdir(prefix="t"):
    """测试临时目录：建于工程 outputs/runs/ 下（沙箱允许 makedirs，拒绝 tempfile 的 0o700）。"""
    import shutil
    import uuid
    base = Path(__file__).resolve().parent.parent / "outputs" / "runs"
    d = base / ("ut_%s_%s" % (prefix, uuid.uuid4().hex[:10]))
    d.mkdir(parents=True, exist_ok=False)
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from sirna_pipeline.common import seqio, duplex, mismatch_class, nn_tables, records, stage_io


class TestSeqio(unittest.TestCase):
    def test_to_rna_normalize(self):
        self.assertEqual(seqio.to_rna("AATTGGCC"), "AAUUGGCC")
        self.assertEqual(seqio.normalize(" aa tt gc cc "), "AAUUGCCC")

    def test_rc_rna(self):
        self.assertEqual(seqio.rc_rna("AUGC"), "GCAU")
        self.assertEqual(seqio.rc_rna(seqio.rc_rna("AUGC")), "AUGC")

    def test_gc_pct(self):
        self.assertEqual(seqio.gc_pct("AAAAAAAAAA"), 0.0)
        self.assertEqual(seqio.gc_pct("GGGGGGGGGG"), 100.0)
        s = "ACGUACGUACGUACGUAUC"      # 19nt, GC = 9
        self.assertAlmostEqual(seqio.gc_pct(s), 9 / 19 * 100, places=6)

    def test_pair_type(self):
        self.assertEqual(seqio.pair_type("A", "U"), "WC")
        self.assertEqual(seqio.pair_type("U", "A"), "WC")
        self.assertEqual(seqio.pair_type("G", "C"), "WC")
        self.assertEqual(seqio.pair_type("G", "U"), "GU")
        self.assertEqual(seqio.pair_type("U", "G"), "GU")
        self.assertEqual(seqio.pair_type("A", "A"), "MM")

    def test_assert_rna(self):
        with self.assertRaises(ValueError):
            seqio.assert_rna("ACGN")
        self.assertEqual(seqio.assert_rna("AUGC"), "AUGC")

    def test_load_fasta_and_windows(self):
        with tmpdir() as td:
            fa = Path(td) / "m.fa"
            fa.write_text(">h1\nATGCATGCATGCATGCATGCATGCATGC\n", encoding="utf-8")
            header, rna = seqio.load_fasta_rna(fa)
            self.assertTrue(header.startswith(">h1"))
            self.assertNotIn("T", rna)
            wins = seqio.sliding_windows(rna, window_len=19)
            self.assertEqual(len(wins), len(rna) - 19 + 1)
            for w in wins:
                self.assertEqual(seqio.rc_rna(w["target"]), w["guide"])
                self.assertEqual(w["sense"], w["target"])
            bad = Path(td) / "bad.fa"
            bad.write_text(">h\nATGCN\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                seqio.load_fasta_rna(bad)

    def test_csv_roundtrip(self):
        with tmpdir() as td:
            p = Path(td) / "t.csv"
            rows = [{"a": "x", "b": None}, {"a": "y", "b": 3}]
            n = seqio.write_csv(p, ["a", "b"], rows)
            self.assertEqual(n, 2)
            # BOM 存在
            raw = p.read_bytes()
            self.assertTrue(raw.startswith(b"\xef\xbb\xbf"))
            back = seqio.read_csv_rows(p)
            self.assertEqual(len(back), 2)


class TestDuplex(unittest.TestCase):
    TARGET = "ACGUACGUACGUACGUACG"  # 19

    def test_paired_register(self):
        # g1(5'端) ↔ target 末位；g19 ↔ target 首位
        self.assertEqual(duplex.paired_mrna_nt(self.TARGET, 1), self.TARGET[18])
        self.assertEqual(duplex.paired_mrna_nt(self.TARGET, 19), self.TARGET[0])
        self.assertEqual(duplex.paired_mrna_nt(self.TARGET, 10), self.TARGET[9])

    def test_mrna57(self):
        cds = "A" * 10 + self.TARGET + "U" * 10          # 仅 39nt，两侧不足 19
        self.assertIsNone(duplex.mrna57_from_cds(cds, 11))
        cds57 = "A" * 19 + self.TARGET + "U" * 19        # 57nt
        got = duplex.mrna57_from_cds(cds57, 20)          # 窗口起点 1-based 20 才双侧满
        self.assertEqual(got, cds57)
        self.assertEqual(got[19:38], self.TARGET)
        # 首窗(1)与末窗(57-19+1=39)侧翼不足 -> None
        self.assertIsNone(duplex.mrna57_from_cds(cds57, 1))
        self.assertIsNone(duplex.mrna57_from_cds(cds57, 39))

    def test_seed(self):
        g = "ACGUACGUACGUACGUACG"
        self.assertEqual(duplex.guide_seed(g, "seed6"), g[1:7])
        self.assertEqual(duplex.guide_seed(g, "seed7"), g[1:8])
        with self.assertRaises(ValueError):
            duplex.guide_seed(g, "seed8")
        self.assertEqual(duplex.MUTATION_SITES, (1, 12, 17, 18, 19))

    def test_window_constants(self):
        self.assertEqual(duplex.WINDOW_LEN, 19)
        self.assertEqual(duplex.FLANK, 19)
        self.assertEqual(duplex.MRNA57_LEN, 57)


class TestMismatchClass(unittest.TestCase):
    def test_classes(self):
        self.assertEqual(mismatch_class.mismatch_class("A", "U"), "WC")
        self.assertEqual(mismatch_class.mismatch_class("G", "U"), "GU")
        self.assertEqual(mismatch_class.mismatch_class("A", "G"), "PP")
        self.assertEqual(mismatch_class.mismatch_class("C", "U"), "YY")
        self.assertEqual(mismatch_class.mismatch_class("A", "C"), "PY")
        self.assertEqual(mismatch_class.mismatch_class("C", "A"), "PY")

    def test_column_classes_length(self):
        g = "ACGUACGUACGUACGUACG"
        cols = mismatch_class.column_classes(g, seqio.rc_rna(g))
        self.assertEqual(len(cols), 19)
        self.assertTrue(all(c == "WC" for c in cols))


class TestNNTables(unittest.TestCase):
    def test_values(self):
        self.assertEqual(nn_tables.STACK_DG["GC"], -3.42)
        self.assertEqual(nn_tables.STACK_DG["CG"], -2.36)
        self.assertEqual(nn_tables.STACK_DG["AA"], nn_tables.STACK_DG["UU"])
        self.assertEqual(nn_tables.INIT_DG, 4.09)
        self.assertEqual(nn_tables.END_AU_DG, 0.45)
        self.assertEqual(nn_tables.SYM_DG, 0.43)
        # 与旧模块 thermo_calculator STACK_DG/DH 抽值一致
        self.assertEqual(nn_tables.STACK_DH["GC"], -14.88)


class TestRecords(unittest.TestCase):
    def test_num(self):
        self.assertIsNone(records.num(None))
        self.assertIsNone(records.num(""))
        self.assertIsNone(records.num("nan"))
        self.assertEqual(records.num("3.5"), 3.5)

    def test_to_int01(self):
        self.assertEqual(records.to_int01(True), "1")
        self.assertEqual(records.to_int01(0), "0")
        self.assertEqual(records.to_int01(None), "")

    def test_validate_base(self):
        good = [{
            "window_id": "W0001", "variant_id": "W0001_wt", "kind": "wt",
            "cds_start": "1", "target_mRNA_19": "ACGUACGUACGUACGUACG",
            "guide_wt_19": "CGUACGUACGUACGUACGU", "guide_checked": "CGUACGUACGUACGUACGU",
        }]
        records.validate_base(good, "test")          # 不抛
        bad = dict(good[0]); bad["guide_checked"] = "ACGN"
        with self.assertRaises(ValueError):
            records.validate_base([bad], "test")
        dup = [good[0], dict(good[0])]
        with self.assertRaises(ValueError):
            records.validate_base(dup, "test")

    def test_aliases_and_roundtrip(self):
        with tmpdir() as td:
            legacy = Path(td) / "legacy.csv"
            legacy.write_text(
                "window_id,guide_antisense_19,seq_pass\nW0001,ACGU,1\n",
                encoding="utf-8-sig")
            rows = records.load_with_aliases(legacy)
            self.assertIn("guide_checked", rows[0])
            self.assertIn("rules_pass", rows[0])
            self.assertNotIn("seq_pass", rows[0])
            # 统一记录写读往返
            out = Path(td) / "rec.csv"
            records.write_records(out, rows)
            back = records.read_records(out)
            self.assertEqual(back[0]["rules_pass"], "1")

    def test_canonical_superset(self):
        # 分组拼接后无重复、无遗漏定义冲突
        self.assertEqual(len(records.CANONICAL_COLS), len(records.CANONICAL_SET))
        for c in records.BASE_COLS:
            self.assertIn(c, records.CANONICAL_SET)


class TestStageIO(unittest.TestCase):
    def test_json_and_yaml_config(self):
        with tmpdir() as td:
            jp = Path(td) / "c.json"
            jp.write_text('{"a": 1}', encoding="utf-8")
            self.assertEqual(stage_io.load_config(jp)["a"], 1)
            yp = Path(td) / "c.yaml"
            yp.write_text("a: 1\nb:\n  - x\n  - y\n", encoding="utf-8")
            cfg = stage_io.load_config(yp)
            self.assertEqual(cfg["a"], 1)
            self.assertEqual(cfg["b"], ["x", "y"])

    def test_resolve_paths(self):
        cfg = {"data_root": "data", "abs": "C:/abs/win", "rel": "runs"}
        out = stage_io.resolve_paths(cfg, "C:/base")
        self.assertTrue(str(out["data_root"]).replace("\\", "/").endswith("base/data"))
        self.assertEqual(str(out["abs"]).replace("\\", "/"), "C:/abs/win")
        self.assertTrue(str(out["rel"]).replace("\\", "/").endswith("base/runs"))

    def test_manifest_and_sha(self):
        with tmpdir() as td:
            p = Path(td)
            m = {"n": 3, "k": "中"}
            stage_io.write_manifest(p, m)
            self.assertEqual(stage_io.read_manifest(p)["n"], 3)
            f = p / "f.bin"
            f.write_bytes(b"abc")
            h1 = stage_io.file_sha256(f)
            self.assertEqual(len(h1), 64)
            self.assertEqual(h1, stage_io.file_sha256(f))

    def test_logging(self):
        with tmpdir() as td:
            logf = Path(td) / "l.log"
            lg = stage_io.setup_logging("t.stage", logf)
            lg.info("hello 中文")
            for h in list(lg.handlers):
                h.close()
                lg.removeHandler(h)
            self.assertTrue(logf.exists())
            self.assertIn("hello", logf.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
