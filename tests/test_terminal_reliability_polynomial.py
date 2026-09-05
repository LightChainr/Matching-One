import copy
from fractions import Fraction
import json
import unittest

from scripts.terminal_reliability_polynomial import (
    DEFAULT_OUTPUT,
    STAR4,
    build_star4_result,
    enumerate_reliability,
    evaluate_at,
    validate_gadget,
    validate_result,
)


class TerminalReliabilityPolynomialTests(unittest.TestCase):
    def test_star_partition_counts_and_normalization(self):
        result = build_star4_result()
        self.assertEqual(result["enumeration"]["configurations"], 16)
        self.assertEqual(result["enumeration"]["nonzero_terminal_partitions"], 12)
        self.assertEqual(result["normalization_counts"], [1, 4, 6, 4, 1])
        self.assertEqual(sum(Fraction(row["probability_at_one_half"]) for row in result["terminal_partition_polynomials"]), 1)

    def test_star_disconnected_and_connected_endpoints(self):
        counts = enumerate_reliability(STAR4)
        disconnected = counts[(0, 1, 2, 3)]
        connected = counts[(0, 0, 0, 0)]
        self.assertEqual(disconnected, (1, 4, 0, 0, 0))
        self.assertEqual(connected, (0, 0, 0, 0, 1))
        self.assertEqual(evaluate_at(disconnected, Fraction(0)), 1)
        self.assertEqual(evaluate_at(connected, Fraction(1)), 1)

    def test_two_terminal_series_control(self):
        gadget = {
            "vertex_count": 3,
            "terminal_count": 2,
            "edges": [[0, 2], [1, 2]],
        }
        counts = enumerate_reliability(gadget)
        self.assertEqual(counts[(0, 0)], (0, 0, 1))
        self.assertEqual(counts[(0, 1)], (1, 2, 0))

    def test_checked_in_result_reproduces_exactly(self):
        checked_in = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked_in, build_star4_result())
        summary = validate_result(checked_in)
        self.assertEqual(summary["normalization_counts"], [1, 4, 6, 4, 1])

    def test_perturbed_polynomial_and_claim_fail_closed(self):
        result = build_star4_result()
        polynomial_tamper = copy.deepcopy(result)
        polynomial_tamper["terminal_partition_polynomials"][0]["bernstein_counts_by_open_edge_count"][0] += 1
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(polynomial_tamper)
        claim_tamper = copy.deepcopy(result)
        claim_tamper["conclusion"]["new_percolation_bound"] = True
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(claim_tamper)


if __name__ == "__main__":
    unittest.main()
