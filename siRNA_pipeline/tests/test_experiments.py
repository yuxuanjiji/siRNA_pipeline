import json
import unittest
from pathlib import Path

try:
    from tests._util import tmpdir
except ImportError:  # pragma: no cover
    from _util import tmpdir

from scripts.run_experiments import _ablation, _mismatch


class TestExperimentReports(unittest.TestCase):
    def test_mismatch_count_is_reported(self):
        with tmpdir("exp_mm") as td:
            source = Path(td) / "mismatch.csv"
            source.write_text("label,score\n1,0.8\n0,0.2\n", encoding="utf-8")
            out = Path(td) / "out"
            args = type("Args", (), {"input": source, "prediction_file": None,
                                     "prediction": "score", "label": "label",
                                     "side": "mismatch_side", "kind": "holen",
                                     "with_mfe": False, "out": out})()
            result = _mismatch(args)
            self.assertEqual(result["rows"], 2)
            self.assertIn("row_count_warning", result)

    def test_ablation_grid_has_six_points(self):
        with tmpdir("exp_ab") as td:
            source = Path(td) / "features.csv"
            source.write_text("label,dl_score,thermo_score\n1,0.9,0.8\n0,0.1,0.2\n", encoding="utf-8")
            out = Path(td) / "out"
            args = type("Args", (), {"input": source, "with_mfe": False, "out": out})()
            result = _ablation(args)
            self.assertEqual(len(result["grid"]), 6)
            self.assertTrue((out / "task16_ablation.json").exists())
            json.loads((out / "task16_ablation.json").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
