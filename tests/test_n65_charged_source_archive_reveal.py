from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from reveal_n65_charged_source_archive import (  # noqa: E402
    BirthCell,
    LINE_CHARGES,
    evaluate_batch,
    projectivize_f3,
    quadratic_n,
)


def transform(line: tuple[int, int], matrix: tuple[tuple[int, int], tuple[int, int]]) -> tuple[int, int]:
    x, y = line
    return projectivize_f3(
        matrix[0][0] * x + matrix[0][1] * y,
        matrix[1][0] * x + matrix[1][1] * y,
    )


class N65ChargedSourceRevealTests(unittest.TestCase):
    def test_R_charge_and_D4_reflection_blocks(self) -> None:
        quarter_turn = ((0, -1), (1, 0))
        reflection = ((1, 0), (0, -1))
        for line, charges in LINE_CHARGES.items():
            rotated = LINE_CHARGES[transform(line, quarter_turn)]
            self.assertEqual(rotated["q_A"], -charges["q_A"])
            self.assertEqual(rotated["q_D"], -charges["q_D"])
            reflected = LINE_CHARGES[transform(line, reflection)]
            self.assertEqual(reflected["q_A"], charges["q_A"])
            self.assertEqual(reflected["q_D"], -charges["q_D"])

    def test_fixed_sources_and_current_continuity(self) -> None:
        cells = [
            BirthCell("first", 0, 5, 1, 4, "LINE", 1, 0, 1),
            BirthCell("first", 0, 5, 1, 4, "LINE", 0, 1, 1),
            BirthCell("first", 0, 5, 1, 4, "LINE", 1, 1, 1),
            BirthCell("first", 0, 5, 1, 4, "LINE", 1, -1, 1),
            BirthCell("first", 0, 5, 3, 3, "DIRECT_RANK2", 0, 0, 1),
        ]
        metrics, gates = evaluate_batch(cells, 5, 0.4)
        self.assertGreater(metrics["W_A"], 0.0)
        self.assertGreater(metrics["W_D"], 0.0)
        self.assertAlmostEqual(metrics["W_A"], metrics["W_D"])
        self.assertAlmostEqual(metrics["A_unit"], 0.0)
        self.assertAlmostEqual(metrics["D_unit"], 0.0)
        self.assertEqual(metrics["response_A_from_D"], 0.0)
        self.assertEqual(metrics["response_D_from_A"], 0.0)
        self.assertLess(max(abs(value) for value in gates.values()), 1e-12)

    def test_small_dense_quadratic(self) -> None:
        self.assertAlmostEqual(quadratic_n([1.0, 2.0], [[2.0, 0.0], [0.0, 4.0]]), 1.5)


if __name__ == "__main__":
    unittest.main()
