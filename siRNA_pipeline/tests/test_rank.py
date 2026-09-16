# -*- coding: utf-8 -*-
"""Stage7 rank 整合测试（新模块）。

- 手工构造数据：归一化边界(最好=1/最差=0)、软惩罚、缺特征排除、rank 顺序；
- 合成全链 smoke：gen->rules->structure(approx,不要求pass)->thermo->toxicity->rank；
- 真实 SFRP1 抽样回归（SIRNA_REGRESSION=1）：可复现、非空、排名单调。

运行：python -m unittest tests.test_rank -v
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
from sirna_pipeline.stages.rank.ranker import run_rank, merge_tables
from sirna_pipeline.stages.generation.runner import run_generation
from sirna_pipeline.stages.rule_filter.runner import run_rules
from sirna_pipeline.stages.structure.runner import run_structure
from sirna_pipeline.stages.thermo.runner import run_thermo
from sirna_pipeline.stages.toxicity.runner import run_toxicity

_J19_ZERO = json.dumps([0.0] * 19)
_J19_SEED14 = json.dumps([0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0] + [0.0] * 11)


def _base_row(wid: str, kind: str, guide: str) -> dict:
    return {
        "window_id": wid, "variant_id": (wid + "_wt" if kind == "wt"
                                         else wid + "_g1:A>U"),
        "kind": kind, "cds_start": "20", "cds_end": "38",
        "target_mRNA_19": "ACGUACGUACGUACGUACG",
        "guide_wt_19": guide, "guide_checked": guide, "sense_strand_19": "",
        "rules_pass": "1",
    }


class TestRankHandcrafted(unittest.TestCase):
    def test_normalization_penalty_and_exclusion(self):
        with tmpdir("rank_craft") as td:
            g_best = "CCCUUCUACGUAACCUAAA"
            g_worst = "AUAUAUAUAUAUAUAUAUA"
            rows = [
                {**_base_row("W0001", "wt", g_best),
                 "mfe": "-1", "delta_deltaG_ends": "0.5", "dG_total": "-45",
                 "dG_mismatch_by_position": _J19_ZERO,
                 "tox_viability_flag": "0", "imm_high_flag": "0"},
                {**_base_row("W0002", "wt", g_worst),
                 "mfe": "-9", "delta_deltaG_ends": "-5", "dG_total": "-20",
                 "dG_mismatch_by_position": _J19_SEED14,
                 "tox_viability_flag": "0", "imm_high_flag": "0"},
                {**_base_row("W0003", "wt", g_worst),
                 "mfe": "-9", "delta_deltaG_ends": "-5", "dG_total": "-20",
                 "dG_mismatch_by_position": _J19_SEED14,
                 "tox_viability_flag": "1", "imm_high_flag": "0"},   # 仅被惩罚
                {**_base_row("W0004", "wt", g_worst),
                 "mfe": "-9", "delta_deltaG_ends": "-5", "dG_total": "",   # 缺特征
                 "dG_mismatch_by_position": _J19_SEED14,
                 "tox_viability_flag": "0", "imm_high_flag": "0"},
            ]
            base = td / "base.csv"
            records.write_records(base, rows)
            out_csv, meta = run_rank(base, out_dir=td / "rank",
                                     cfg={"require_structure_pass": False})
            out = records.read_records(out_csv)
            by = {r["window_id"]: r for r in out}
            # W0001 四特征全最优 -> score=1；W0002 全最差 -> score=0
            self.assertEqual(records.num(by["W0001"]["score_thermo"]), 1.0)
            self.assertEqual(records.num(by["W0002"]["score_thermo"]), 0.0)
            # 惩罚：W0003 与 W0002 特征相同但毒性命中 -> 终分更低
            f2 = records.num(by["W0002"]["final_score"])
            f3 = records.num(by["W0003"]["final_score"])
            self.assertAlmostEqual(f2, 0.0)
            self.assertAlmostEqual(f3, -0.05, places=9)
            self.assertEqual(records.num(by["W0003"]["penalty_total"]), 0.05)
            # 缺特征行不参与排序
            self.assertEqual(by["W0004"]["final_rank"], "")
            self.assertEqual(meta["excluded"].get("特征缺失"), 1)
            self.assertEqual(meta["eligible"], 3)
            # 排名降序 & 命中 1..3
            scored = sorted((r for r in out if r["final_rank"]),
                            key=lambda r: int(r["final_rank"]))
            self.assertEqual([r["window_id"] for r in scored],
                             ["W0001", "W0002", "W0003"])
            # 派生原始特征列已写入
            self.assertEqual(records.num(by["W0001"]["feat_dG_duplex"]), -45.0)
            # 特征统计
            st = meta["features_stats"]["feat_mfe"]
            self.assertEqual(st["min"], -9.0)
            self.assertEqual(st["max"], -1.0)

    def test_merge_tables_updates_by_key(self):
        with tmpdir("rank_merge") as td:
            base = td / "b.csv"
            rows = [_base_row("W0001", "wt", "CCCUUCUACGUAACCUAAA"),
                    _base_row("W0002", "wt", "AUAUAUAUAUAUAUAUAUA")]
            records.write_records(base, rows)
            extra = td / "s.csv"
            extra_rows = [
                {**records.read_records(base)[0], "mfe": "-3"},
                {**records.read_records(base)[1], "mfe": "-8"},
            ]
            records.write_records(extra, extra_rows)
            merged = merge_tables(base, [extra])
            m = {r["window_id"]: r for r in merged}
            self.assertEqual(m["W0001"]["mfe"], "-3")
            self.assertEqual(len(merged), 2)


class TestRankChainSmoke(unittest.TestCase):
    def test_full_chain_smoke(self):
        with tmpdir("rank_smoke") as td:
            fa = td / "syn.fa"
            rng = random.Random(13)
            seq = "".join(rng.choice("ATGC") for _ in range(120))  # 保证多数窗口有 57nt 上下文
            fa.write_text(">syn\n%s\n" % seq, encoding="utf-8")
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fa, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            passed = rules_dir / "candidates_passed.csv"
            struct_csv = run_structure(passed, fa, td / "struct",
                                       detector_kwargs={"backend": "approx"})
            thermo_csv = run_thermo(passed, fa, td / "thermo")
            tox_csv = run_toxicity(passed, td / "tox")
            out_csv, meta = run_rank(
                passed, [struct_csv, thermo_csv, tox_csv], td / "rank",
                cfg={"require_structure_pass": False})
            out = records.read_records(out_csv)
            ranked = [r for r in out if r["final_rank"]]
            self.assertGreater(len(ranked), 0)
            # 按 final_rank 顺序的 final_score 非增（排名单调）
            ordered = sorted(ranked, key=lambda r: int(r["final_rank"]))
            scores = [records.num(r["final_score"]) for r in ordered]
            self.assertEqual(scores, sorted(scores, reverse=True))
            # 未排序列结构完整
            self.assertEqual(len(out), len(records.read_records(passed)))


@unittest.skipUnless(os.environ.get("SIRNA_REGRESSION") == "1",
                     "真实回归测试需 SIRNA_REGRESSION=1")
class TestRankRegressionRealData(unittest.TestCase):
    def test_real_sample_deterministic(self):
        ws = Path(__file__).resolve().parents[2]
        fasta = ws / "数据集" / "SFRP1-mRNA.txt"
        with tmpdir("rank_real") as td:
            gen_dir, rules_dir = td / "gen", td / "rules"
            run_generation(fasta=fasta, out_dir=gen_dir)
            run_rules(gen_dir / "candidates.csv", rules_dir)
            passed_rows = records.read_records(rules_dir / "candidates_passed.csv")
            sample = passed_rows[:300]
            sp = td / "passed_300.csv"
            records.write_records(sp, sample)
            st = run_structure(sp, fasta, td / "struct",
                               detector_kwargs={"backend": "auto"})
            th = run_thermo(sp, fasta, td / "thermo")
            tx = run_toxicity(sp, td / "tox")
            out1, meta1 = run_rank(sp, [st, th, tx], td / "r1",
                                   cfg={"require_structure_pass": False})
            out2, meta2 = run_rank(sp, [st, th, tx], td / "r2",
                                   cfg={"require_structure_pass": False})
            r1 = {r["variant_id"]: r["final_rank"]
                  for r in records.read_records(out1) if r["final_rank"]}
            r2 = {r["variant_id"]: r["final_rank"]
                  for r in records.read_records(out2) if r["final_rank"]}
            self.assertEqual(r1, r2)                # 完全可复现
            self.assertEqual(meta1["eligible"], meta2["eligible"])
            self.assertGreater(meta1["eligible"], 0)


class TestRankAuxAndSeed(unittest.TestCase):
    """第五步增强：OligoFormer 辅助分（存在/优雅降级）与 seed 结合能特征非退化。"""

    def _rows(self, with_aux: bool) -> list[dict]:
        g_best = "CCCUUCUACGUAACCUAAA"
        g_worst = "AUAUAUAUAUAUAUAUAUA"
        rows = [
            {**_base_row("W0001", "wt", g_best),
             "mfe": "-1", "delta_deltaG_ends": "0.5", "dG_total": "-45",
             "dG_mismatch_by_position": _J19_ZERO,
             "tox_viability_flag": "0", "imm_high_flag": "0"},
            {**_base_row("W0002", "wt", g_worst),
             "mfe": "-9", "delta_deltaG_ends": "-5", "dG_total": "-20",
             "dG_mismatch_by_position": _J19_SEED14,
             "tox_viability_flag": "0", "imm_high_flag": "0"},
            {**_base_row("W0003", "wt", g_worst),
             "mfe": "-9", "delta_deltaG_ends": "-5", "dG_total": "-20",
             "dG_mismatch_by_position": _J19_SEED14,
             "tox_viability_flag": "1", "imm_high_flag": "0"},
        ]
        if with_aux:
            for r, v in zip(rows, ("0.9", "0.2", "0.8")):
                r["oligo_efficacy"] = v
        return rows

    def test_seed_feature_non_degenerate(self):
        # 直接验证 seed 结合能特征：GC 富集 seed 比 A/U 富集 seed 更负（更强），且都 <0
        from sirna_pipeline.stages.rank import ranker
        r_gc = {"guide_checked": "GCCGCCGCCGCCGCCGCCG",
                "dG_mismatch_by_position": _J19_ZERO}
        r_at = {"guide_checked": "AUAUAUAUAUAUAUAUAUA",
                "dG_mismatch_by_position": _J19_ZERO}
        v_gc = ranker._seed_binding(r_gc, "seed7")
        v_at = ranker._seed_binding(r_at, "seed7")
        self.assertIsNotNone(v_gc)
        self.assertIsNotNone(v_at)
        self.assertLess(v_gc, 0.0)
        self.assertLess(v_at, 0.0)
        self.assertLess(v_gc, v_at)

    def test_aux_present_affects_final(self):
        with tmpdir("rank_aux") as td:
            base = td / "b.csv"
            records.write_records(base, self._rows(with_aux=True))
            out_csv, meta = run_rank(
                base, out_dir=td / "r",
                cfg={"require_structure_pass": False})
            self.assertTrue(meta["aux"]["available"])
            self.assertAlmostEqual(float(meta["alpha"]), 0.4)
            self.assertAlmostEqual(float(meta["beta"]), 0.6)
            out = {r["window_id"]: r for r in records.read_records(out_csv)}
            # W0001：主分=1、辅助分=1 -> final=0.4+0.6=1
            self.assertAlmostEqual(records.num(out["W0001"]["final_score"]), 1.0)
            # W0003：主分=0、辅助分=(0.8-0.2)/0.7、毒性惩罚 0.05
            exp3 = 0.4 * 0.0 + 0.6 * (0.6 / 0.7) - 0.05
            self.assertAlmostEqual(records.num(out["W0003"]["final_score"]),
                                   exp3, places=9)
            self.assertAlmostEqual(records.num(out["W0002"]["final_score"]), 0.0)
            # 辅助分影响排名：W0003(辅0.8) 反超 W0002(辅0.2)，尽管热力特征相同且 W0003 有毒
            scored = sorted((r for r in out.values() if r["final_rank"]),
                            key=lambda r: int(r["final_rank"]))
            self.assertEqual([r["window_id"] for r in scored],
                             ["W0001", "W0003", "W0002"])

    def test_aux_absent_fallback(self):
        with tmpdir("rank_noaux") as td:
            base = td / "b.csv"
            records.write_records(base, self._rows(with_aux=False))
            out_csv, meta = run_rank(
                base, out_dir=td / "r",
                cfg={"require_structure_pass": False})
            self.assertFalse(meta["aux"]["available"])
            self.assertAlmostEqual(float(meta["alpha"]), 1.0)
            self.assertAlmostEqual(float(meta["beta"]), 0.0)
            out = records.read_records(out_csv)
            for r in out:
                pen = records.num(r.get("penalty_total", "")) or 0.0
                self.assertAlmostEqual(records.num(r["final_score"]),
                                       (records.num(r["score_thermo"]) or 0.0) - pen,
                                       places=9)


class TestRankVariantLayer(unittest.TestCase):
    """变体层（窗口内排序 = 锚点 × C_match）：行为、向后兼容与一键回退。

    关键不变式：
      1) WT 行 C_match≡1 → 分数与旧口径**逐位一致**（这是"不破坏冻结结果"的验收条件）；
      2) 同窗口内变体只按 C_match 排序 —— 即便某变体自身 Q_guide 更差也必须排在前面
         （旧口径恰好相反，故两种口径的排序会反转，可据此验证开关有效）；
      3) `enabled: False` 完全回退旧口径。
    """

    GUIDE = "ACGUACGUACGUACGUACG"

    def _row(self, wid: str, tag: str, feats: dict,
             pos: int | None = None, mut: str = "", paired: str = "") -> dict:
        is_wt = tag == "wt"
        row = {
            "window_id": wid,
            "variant_id": "%s_%s" % (wid, tag),
            "kind": "wt" if is_wt else "mut",
            "guide_wt_19": self.GUIDE,
            "guide_checked": self.GUIDE,
            "rules_pass": "1",
            "mfe": feats["mfe"], "delta_deltaG_ends": feats["ddg"],
            "dG_total": feats["dG"],
            "tox_viability_flag": "0", "imm_high_flag": "0",
        }
        if not is_wt:
            row.update({"position": str(pos), "mut_nt": mut,
                        "paired_mRNA_nt": paired, "n_mismatch": "1"})
        return row

    _BEST = {"mfe": "-1", "ddg": "0.5", "dG": "-45"}
    _WORST = {"mfe": "-9", "ddg": "-5", "dG": "-20"}
    _MID = {"mfe": "-5", "ddg": "0", "dG": "-30"}

    def _table(self) -> list[dict]:
        return [
            # W0001：WT(MID) + 4 变体
            #   g12 PP / MID → 只体现 C_match 代价（λ=0.15）
            #   g19 GU / MID → λ=0（文献：3′ 端可增益位点）→ C_match=1
            #   g17 GU / BEST → 门控位点 + 热力改善 → 应可超过 WT
            #   g10 transversion / BEST → **非门控位点**（中央区）→ 有代价且不得增益
            self._row("W0001", "wt", self._MID),
            self._row("W0001", "g12:C>A", self._MID, pos=12, mut="A", paired="G"),
            self._row("W0001", "g19:A>U", self._MID, pos=19, mut="U", paired="G"),
            self._row("W0001", "g17:C>U", self._BEST, pos=17, mut="U", paired="G"),
            self._row("W0001", "g10:C>A", self._BEST, pos=10, mut="A", paired="G"),
            # W0002：无 WT 行（模拟 WT 被前序硬过滤）→ 锚点退化取窗口内最高
            self._row("W0002", "g17:C>U", self._BEST, pos=17, mut="U", paired="G"),
            self._row("W0002", "g12:C>A", self._WORST, pos=12, mut="A", paired="G"),
        ]

    def _run(self, td, enabled: bool):
        base = td / ("base_%s.csv" % enabled)
        records.write_records(base, self._table())
        return run_rank(base, out_dir=td / ("r_%s" % enabled),
                        cfg={"require_structure_pass": False,
                             "variant_layer": {"enabled": enabled}})

    def test_cmatch_values(self):
        """C_match 取值：λ 逐带（文献校正后）+ AGO2 的 transition/transversion 轴。"""
        from sirna_pipeline.stages.rank import ranker
        vcfg = dict(ranker.DEFAULT_RANK_CFG["variant_layer"])
        # λ：g1=0（文献：该位在 RISC 中不配对且摆动可增益）、seed=0.1786、cleavage=0.3082、
        #     g12=0.15、mid=0.1685、p3=0（文献：3′ 端错配可增益）
        cases = {
            # (position, mut_nt, paired_mRNA_nt) -> C_match
            (12, "A", "G"): 1 - 0.1500 * 0.7,      # pos≥12 段 transversion 因子 0.7
            (12, "A", "C"): 1 - 0.1500 * 1.0,      # pos≥12 段 transition 因子 1.0（AGO2 符号反转）
            (12, "U", "G"): 1 - 0.1500 * 0.5,      # GU 摆动（U:G）
            (10, "A", "G"): 1 - 0.3082 * 1.0,      # 中央区 transversion（pos<12 因子 1.0）
            (10, "A", "C"): 1 - 0.3082 * 0.85,     # 中央区 transition
            (5, "A", "G"): 1 - 0.1786 * 1.0,       # seed
            (1, "U", "G"): 1.0,                    # g1：λ=0（文献支持增益位点）
            (17, "U", "G"): 1.0,                   # p3：λ=0
            (18, "A", "G"): 1.0,
            (19, "U", "U"): 1.0,
        }
        for (pos, mut, paired), exp in cases.items():
            got, detail = ranker.c_match({"position": str(pos), "mut_nt": mut,
                                          "paired_mRNA_nt": paired}, vcfg)
            self.assertAlmostEqual(got, exp, places=9,
                                   msg="pos=%d %s:%s" % (pos, mut, paired))
            self.assertEqual(detail[0]["pos"], pos)
        # 类型轴分类（AGO2 口径：以"靶碱基相对其互补碱基"的变化判定）
        self.assertEqual(ranker.mismatch_axis_class("G", "C"), "WC")
        self.assertEqual(ranker.mismatch_axis_class("G", "U"), "GU")
        self.assertEqual(ranker.mismatch_axis_class("A", "C"), "transition")     # U→C
        self.assertEqual(ranker.mismatch_axis_class("A", "G"), "transversion")   # U→G
        # ≥3 错配额外折半（Holen：三重突变几乎耗尽活性）
        got3, det3 = ranker.c_match({"position": "12", "mut_nt": "A",
                                     "paired_mRNA_nt": "G", "n_mismatch": "3"}, vcfg)
        self.assertAlmostEqual(got3, (1 - 0.1500 * 0.7) * 0.5, places=9)
        self.assertTrue(det3[0]["multi_mismatch_ge3"])
        # WT（无 position 列）→ 恒 1.0，不惩罚
        got, detail = ranker.c_match({"variant_id": "W1_wt"}, vcfg)
        self.assertEqual((got, detail), (1.0, []))

    def test_wt_rows_bitwise_identical_and_variant_order_flips(self):
        with tmpdir("rank_vlayer") as td:
            on_csv, on_meta = self._run(td, True)
            off_csv, off_meta = self._run(td, False)
            on = {r["variant_id"]: r for r in records.read_records(on_csv)}
            off = {r["variant_id"]: r for r in records.read_records(off_csv)}

            # (1) WT 行逐位一致（验收条件：不破坏冻结结果；WT 的 Δthermo≡0、C_match≡1）
            for vid in ("W0001_wt",):
                for col in ("score_thermo", "score_oligo", "penalty_total", "final_score"):
                    self.assertEqual(on[vid][col], off[vid][col],
                                     "%s.%s 在两种口径下不一致" % (vid, col))
            self.assertEqual(on["W0001_wt"]["c_match"], "1")
            self.assertEqual(on["W0001_wt"]["thermo_delta_vs_wt"], "0")

            # (2) 变体：final = q_eff × C_match − penalty；
            #     q_eff = 锚点 + α·Δthermo（**仅门控位点**，见 gain_positions）
            anchor = records.num(on["W0001_wt"]["q_window"])
            alpha = float(on_meta["alpha"])
            gain_pos = on_meta["variant_layer"]["gain_positions"]
            for vid, exp_cm in (("W0001_g12:C>A", 1 - 0.1500 * 0.7),   # pos≥12 transversion
                                ("W0001_g19:A>U", 1.0),                # p3 λ=0
                                ("W0001_g17:C>U", 1.0),                # p3 λ=0
                                ("W0001_g10:C>A", 1 - 0.3082 * 1.0)):  # 中央区 transversion
                r = on[vid]
                d_th = records.num(r["thermo_delta_vs_wt"])
                gated = int(records.num(r["position"])) in gain_pos
                q_eff = anchor + (alpha * d_th if gated else 0.0)
                self.assertAlmostEqual(records.num(r["q_window_anchor"]), anchor, places=9)
                self.assertAlmostEqual(records.num(r["c_match"]), exp_cm, places=9)
                self.assertAlmostEqual(records.num(r["q_variant_effective"]), q_eff, places=9)
                self.assertAlmostEqual(
                    records.num(r["final_score"]),
                    q_eff * exp_cm - (records.num(r["penalty_total"]) or 0.0), places=9)

            # (3) C_match 在"同热力"变体之间决定次序：无代价者胜有代价者
            self.assertGreater(records.num(on["W0001_g19:A>U"]["final_score"]),
                               records.num(on["W0001_g12:C>A"]["final_score"]))
            self.assertLess(records.num(on["W0001_g12:C>A"]["final_score"]), anchor)
            self.assertAlmostEqual(records.num(off["W0001_g12:C>A"]["final_score"]), anchor,
                                   places=9)

            # (4) **改善通道**：门控位点(g17)的热力改善必须能超过 WT；非门控位点(g10)不得增益
            self.assertGreater(records.num(on["W0001_g17:C>U"]["thermo_delta_vs_wt"]), 0.0)
            self.assertGreater(records.num(on["W0001_g17:C>U"]["final_score"]), anchor)
            # g10 同为 BEST 热力但不在 gain_positions → q_eff 不含增益，且带中央区代价
            self.assertAlmostEqual(
                records.num(on["W0001_g10:C>A"]["q_variant_effective"]), anchor, places=9)
            self.assertLess(records.num(on["W0001_g10:C>A"]["final_score"]), anchor)
            self.assertGreater(on_meta["variant_layer"]["n_variant_above_own_wt"], 0)
            self.assertIn(17, on_meta["variant_layer"]["gain_positions"])
            self.assertNotIn(10, on_meta["variant_layer"]["gain_positions"])
            # γ=0 关闭改善通道后退化为"突变只降效率"（变体不再可能超过 WT）
            off2 = td / "base_g0.csv"
            records.write_records(off2, self._table())
            csv_g0, meta_g0 = run_rank(
                off2, out_dir=td / "r_g0",
                cfg={"require_structure_pass": False,
                     "variant_layer": {"mech_gain_weight": 0.0}})
            g0 = {r["variant_id"]: r for r in records.read_records(csv_g0)}
            self.assertEqual(meta_g0["variant_layer"]["n_variant_above_own_wt"], 0)
            # γ=0（无增益）时，p3 位点 λ=0 → C_match=1 → 该变体与 WT **同分**（不会超过）
            self.assertAlmostEqual(
                records.num(g0["W0001_g17:C>U"]["final_score"]),
                records.num(g0["W0001_wt"]["final_score"]), places=9)

            # (5) 回退等价性：旧口径下 final_score == final_score_base
            for vid in ("W0001_g12:C>A", "W0001_g19:A>U", "W0001_g17:C>U", "W0001_g10:C>A"):
                self.assertEqual(off[vid]["final_score"], off[vid]["final_score_base"])

            # (6) 锚点退化：W0002 无 WT → 计数入 meta，锚点 = 窗口内最高 q_window
            self.assertEqual(on_meta["variant_layer"]["n_windows_fallback_anchor"], 1)
            self.assertEqual(on_meta["variant_layer"]["n_windows_with_wt_anchor"], 1)
            self.assertEqual(on_meta["variant_layer"]["n_variant_rows"], 6)
            self.assertFalse(off_meta["variant_layer"]["enabled"])
            w2 = [r for r in on.values() if r["window_id"] == "W0002"]
            self.assertAlmostEqual(records.num(w2[0]["q_window_anchor"]),
                                   max(records.num(r["q_window"]) for r in w2), places=9)

    def test_improving_variant_rule_is_reserved(self):
        with tmpdir("rank_vlayer_gamma") as td:
            base = td / "b.csv"
            records.write_records(base, self._table())
            # γ=0 → 改善通道关闭，等价于"突变只降效率"；γ 越大越容易超过 WT
            for gamma, expect_above in ((0.0, False), (1.0, True)):
                csv, meta = run_rank(
                    base, out_dir=td / ("g%s" % gamma),
                    cfg={"require_structure_pass": False,
                         "variant_layer": {"mech_gain_weight": gamma}})
                above = meta["variant_layer"]["n_variant_above_own_wt"]
                self.assertEqual(bool(above), expect_above, "γ=%s 时 above=%d" % (gamma, above))
                self.assertAlmostEqual(float(meta["variant_layer"]["mech_gain_weight"]), gamma)


if __name__ == "__main__":
    unittest.main(verbosity=2)
