
from __future__ import annotations
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_self_matching_beta_targets import (  # noqa: E402
    beta_matching_bernstein_counts,
    freeze_payload,
    shortest_nontrivial_cycle,
    vector_sha256,
)


class C4SelfMatchingBetaTargetTests(unittest.TestCase):
    def test_n10_geometry_and_beta33_regression(self) -> None:
        cycle = shortest_nontrivial_cycle(3, 1)
        self.assertEqual(cycle.support, 3)
        self.assertNotEqual(cycle.winding, (0, 0))
        self.assertEqual(
            beta_matching_bernstein_counts(10, 3),
            [-1, -10, -45, -100, -100, 0, 100, 100, 45, 10, 1],
        )

    def test_n26_targets_are_geometry_only_and_frozen(self) -> None:
        payload = freeze_payload(5, 1)
        self.assertEqual(payload["status"], "FROZEN_BEFORE_N26_ENUMERATION")
        self.assertEqual(payload["geometry"]["N"], 26)
        self.assertEqual(payload["geometry"]["wrapping_channel"], "either")
        self.assertEqual(payload["geometry_shortest_cycle_certificate"]["support"], 5)
        hypotheses = {item["name"]: item for item in payload["hypotheses"]}
        self.assertEqual(hypotheses["geometry_shortest_support"]["beta_parameters"], [5, 5])
        self.assertEqual(hypotheses["antipodal_orbit_majority"]["beta_parameters"], [7, 7])
        for hypothesis in hypotheses.values():
            vector = hypothesis["bernstein_integer_coefficients"]
            self.assertEqual(len(vector), 27)
            self.assertEqual(vector, [-value for value in reversed(vector)])
            self.assertEqual(hypothesis["bernstein_vector_sha256"], vector_sha256(vector))


if __name__ == "__main__":
    unittest.main()
