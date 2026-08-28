from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_self_matching_cyclic_geometry import build, validate  # noqa: E402
from c4_self_matching_exact import (  # noqa: E402
    CHANNELS,
    c4_self_matching_torus,
    enumerate_exact,
)


class C4SelfMatchingExactTests(unittest.TestCase):
    def test_lifted_graph_matches_cyclic_graph(self) -> None:
        cyclic = build(3, 1)
        validate(cyclic)
        lifted = c4_self_matching_torus(3, 1)
        labels = [cyclic.a * x + cyclic.b * y for x, y in lifted.coordinates]
        edge_labels = {
            tuple(sorted((labels[edge.i] % cyclic.n, labels[edge.j] % cyclic.n)))
            for edge in lifted.primal_edges
        }
        self.assertEqual(edge_labels, set(cyclic.edges))
        self.assertEqual(lifted.primal_edges, lifted.matching_edges)

    def test_exhaustive_central_identity(self) -> None:
        result = enumerate_exact(3, 1)
        self.assertTrue(result["passed"])
        self.assertEqual(result["geometry"]["N"], 10)
        self.assertEqual(result["geometry"]["configurations"], 1024)
        for channel in CHANNELS:
            row = result["channels"][channel]
            self.assertTrue(row["M_coefficients_anti_palindromic"])
            self.assertTrue(
                row["coefficient_identity_Mk_equals_Rk_minus_RNminusk"]
            )
            self.assertEqual(row["configuration_complement_pair_failures"], 0)
            self.assertEqual(row["M_at_p_one_half"]["fraction"], "0")

    def test_central_zero_does_not_claim_polynomial_zero(self) -> None:
        result = enumerate_exact(3, 1)
        # The exact finite polynomial is odd about 1/2 and nonzero away from
        # the center; self-matching proves the central zero, not M(p) == 0.
        for channel in CHANNELS:
            row = result["channels"][channel]
            self.assertFalse(row["M_polynomial_identically_zero"])
            self.assertNotEqual(row["M_at_p_one_third"]["fraction"], "0")


if __name__ == "__main__":
    unittest.main()
