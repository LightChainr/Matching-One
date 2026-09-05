import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from standard_constant_pairwise_search import DEFAULT_OUTPUT, build_result, primitive_relation_count, validate_result


class StandardConstantPairwiseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_small_count(self):
        self.assertEqual(primitive_relation_count(1), 6)

    def test_committed_result(self):
        checked = json.loads(DEFAULT_OUTPUT.read_text())
        self.assertEqual(checked, self.result)
        self.assertTrue(validate_result(checked)["all_excluded"])

    def test_all_24_cells_excluded(self):
        self.assertEqual(len(self.result["results"]), 24)
        self.assertTrue(all(row["zero_containing_residuals"] == 0 for row in self.result["results"]))


if __name__ == "__main__":
    unittest.main()
