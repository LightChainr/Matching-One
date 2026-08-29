import json
import math
import unittest
from pathlib import Path

from scripts.score_matching_odd_existence_meta import combine


class MatchingOddExistenceMetaTest(unittest.TestCase):
    def test_canonical_ledger_result(self) -> None:
        payload = json.loads(
            Path("results/evidence-ledger/latest.json").read_text(encoding="utf-8")
        )
        result = combine(payload)
        zero = result["models"]["zero_effect"]
        h4 = result["models"]["H4_fixed_predictions"]

        self.assertEqual(zero["degrees_of_freedom"], 4)
        self.assertAlmostEqual(zero["chi_square"], 31.18573555150965, places=12)
        self.assertTrue(math.isclose(zero["chi_square_survival"], 2.805595267905808e-6, rel_tol=1e-12))

        self.assertEqual(h4["degrees_of_freedom"], 4)
        self.assertAlmostEqual(h4["chi_square"], 3.4622795373044295, places=12)
        self.assertAlmostEqual(h4["chi_square_survival"], 0.48363695393249573, places=14)
        self.assertAlmostEqual(
            result["comparison"]["delta_nlpd_H4_minus_zero"],
            -13.806418789729094,
            places=12,
        )

        groups = {row["raw_data_group"] for row in result["blocks"]}
        self.assertEqual(len(groups), 2)


if __name__ == "__main__":
    unittest.main()
