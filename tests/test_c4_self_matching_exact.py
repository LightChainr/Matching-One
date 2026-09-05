
from __future__ import annotations
import math
import sys
import unittest
from fractions import Fraction
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

    def test_n10_matching_cdf_is_exact_beta33(self) -> None:
        """The N=10 control gives M(p)=2 I_p(3,3)-1 exactly.

        This supplies a compact independent regression for polynomial
        conversion, the threshold-rank CDF interpretation, and derivative
        invariants used elsewhere in the project.
        """

        result = enumerate_exact(3, 1)
        bernstein = result["channels"]["either"]["M_bernstein_integer_coefficients"]
        self.assertEqual(
            bernstein,
            [-1, -10, -45, -100, -100, 0, 100, 100, 45, 10, 1],
        )

        # Convert sum_k b_k p^k (1-p)^(N-k) to ascending power coefficients.
        n = len(bernstein) - 1
        power = [0] * (n + 1)
        for k, coefficient in enumerate(bernstein):
            for j in range(n - k + 1):
                power[k + j] += coefficient * ((-1) ** j) * math.comb(n - k, j)
        while power and power[-1] == 0:
            power.pop()

        # 2*BetaCDF(3,3)-1 = 12 p^5 - 30 p^4 + 20 p^3 - 1.
        self.assertEqual(power, [-1, 0, 0, 20, -30, 12])

        half = Fraction(1, 2)

        def derivative(order: int) -> Fraction:
            total = Fraction(0)
            for degree, coefficient in enumerate(power):
                if degree < order:
                    continue
                falling = 1
                for shift in range(order):
                    falling *= degree - shift
                total += coefficient * falling * half ** (degree - order)
            return total

        d1 = derivative(1)
        d3 = derivative(3)
        d5 = derivative(5)
        self.assertEqual(d1, Fraction(15, 4))
        self.assertEqual(d3, Fraction(-60, 1))
        self.assertEqual(d5, Fraction(1440, 1))
        self.assertEqual(d3 / d1**3, Fraction(-256, 225))
        self.assertEqual(d5 / d1**5, Fraction(32768, 16875))


if __name__ == "__main__":
    unittest.main()
