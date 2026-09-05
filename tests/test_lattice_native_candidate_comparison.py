import copy, json, sys, unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lattice_native_candidate_comparison import DEFAULT_OUTPUT, build_result, validate_result


class LatticeNativeCandidateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.result = build_result()
    def test_committed_result(self):
        checked=json.loads(DEFAULT_OUTPUT.read_text()); self.assertEqual(checked,self.result); validate_result(checked)
    def test_four_candidates_are_exactly_excluded(self):
        self.assertEqual(len(self.result["candidates"]),4)
        self.assertTrue(all(row["excluded_by_all_method_intervals"] for row in self.result["candidates"]))
        self.assertTrue(all(len(row["method_comparisons"])==4 for row in self.result["candidates"]))


if __name__ == "__main__": unittest.main()
