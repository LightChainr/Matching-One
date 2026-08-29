from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rigorous_pc_confidence_gate as gate  # noqa: E402


class RigorousPcConfidenceGateTest(unittest.TestCase):
    def test_binomial_tail_boundary_values(self) -> None:
        p = Fraction(3, 10)
        self.assertEqual(gate.binomial_tail(4, 0, p), 1)
        self.assertEqual(gate.binomial_tail(4, 5, p), 0)
        self.assertEqual(gate.binomial_tail(4, 4, p), p**4)
        self.assertEqual(gate.binomial_tail(4, 1, p), 1 - (1 - p) ** 4)

    def test_tail_is_monotone_in_cutoff(self) -> None:
        p = Fraction(7, 10)
        tails = [gate.binomial_tail(20, cutoff, p) for cutoff in range(22)]
        self.assertTrue(all(left >= right for left, right in zip(tails, tails[1:])))

    def test_legacy_paper_cutoff_is_reproduced_exactly(self) -> None:
        cutoff, tail = gate.minimal_successes(400, gate.LEGACY_P0, gate.PER_RUN_ALPHA)
        self.assertEqual(cutoff, 378)
        self.assertLessEqual(tail, gate.PER_RUN_ALPHA)
        self.assertGreater(
            gate.binomial_tail(400, cutoff - 1, gate.LEGACY_P0),
            gate.PER_RUN_ALPHA,
        )
        self.assertAlmostEqual(float(tail), 1.1489903528940095e-7, places=20)

    def test_modern_cutoff_is_373(self) -> None:
        cutoff, tail = gate.minimal_successes(400, gate.MODERN_P0, gate.PER_RUN_ALPHA)
        self.assertEqual(cutoff, 373)
        self.assertLessEqual(tail, gate.PER_RUN_ALPHA)
        self.assertGreater(
            gate.binomial_tail(400, cutoff - 1, gate.MODERN_P0),
            gate.PER_RUN_ALPHA,
        )
        self.assertAlmostEqual(float(tail), 9.514515135926062e-8, places=20)

    def test_modern_power_table(self) -> None:
        cutoff = 373
        expected = {
            Fraction(90, 100): 0.014941716467950962,
            Fraction(92, 100): 0.2056606089348357,
            Fraction(94, 100): 0.7737583835552556,
            Fraction(95, 100): 0.9520076475965301,
        }
        for probability, target in expected.items():
            power = gate.binomial_tail(400, cutoff, probability)
            self.assertAlmostEqual(float(power), target, places=14)

    def test_error_budget_is_familywise(self) -> None:
        self.assertEqual(
            gate.PER_RUN_ALPHA * gate.SIDES * gate.ATTEMPTS_PER_SIDE,
            gate.FAMILYWISE_ALPHA,
        )

    def test_heuristic_cost_scaling_is_separate(self) -> None:
        row = gate.heuristic_cost_row(5e-4, 1e-4)
        self.assertAlmostEqual(row["linear_scale_multiplier"], 5 ** (4 / 3), places=12)
        self.assertAlmostEqual(row["area_work_multiplier"], 5 ** (8 / 3), places=12)

    def test_checked_in_results_are_reproducible(self) -> None:
        artifact = gate.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/rigorous-pc-confidence-gate/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/rigorous-pc-confidence-gate/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, gate.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
