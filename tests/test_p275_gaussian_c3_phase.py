from __future__ import annotations

import importlib.util
from math import atan2, cos, sin
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "p275-gaussian-c3-phase-20260901"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


producer = load("p275_gaussian_c3_producer", EXPERIMENT / "producer.py")
scorer = load("p275_gaussian_c3_score", EXPERIMENT / "score.py")


class P275GaussianC3GeometryTests(unittest.TestCase):
    def test_equal_norm_geometries_and_rotation(self) -> None:
        producer.validate_geometries()
        self.assertAlmostEqual(cos(scorer.DELTA), 12 / 13, places=15)
        self.assertAlmostEqual(sin(scorer.DELTA), 5 / 13, places=15)
        self.assertGreater(sin(12 * scorer.DELTA) ** 2, 0.999)

    def test_counter_field_is_deterministic_and_width_limited(self) -> None:
        first = producer._counter_mask(17, 23, 130)
        self.assertEqual(first, producer._counter_mask(17, 23, 130))
        self.assertNotEqual(first, producer._counter_mask(17, 24, 130))
        self.assertLess(first, 1 << 130)


class P275GaussianC3ScoreTests(unittest.TestCase):
    def synthetic_rows(self, phase: float) -> list[list[float]]:
        a = [0.7, -0.4]
        rot = scorer.rotation(phase)
        z2 = [sum(rot[i][j] * a[j] for j in range(2)) for i in range(2)]
        offsets = [
            [-2, 1, 1, -1],
            [-1, -2, 2, 1],
            [0, 2, -2, 2],
            [1, -1, 1, -2],
            [2, 0, -1, 1],
            [-2, -1, 2, -2],
            [-1, 2, -1, 2],
            [0, -2, 2, -1],
            [1, 1, -2, 1],
            [2, -1, 1, -2],
        ]
        center = a + z2
        return [
            [center[j] + 1e-3 * offset[j] for j in range(4)]
            for offset in offsets
        ]

    def test_h4_phase_selects_h4(self) -> None:
        result = scorer.score(self.synthetic_rows(+4 * scorer.DELTA))
        self.assertEqual(result["decision"], "H4_SELECTED_H8_STOP")

    def test_h8_phase_selects_h8(self) -> None:
        result = scorer.score(self.synthetic_rows(-8 * scorer.DELTA))
        self.assertEqual(result["decision"], "H8_SELECTED_H4_STOP")


if __name__ == "__main__":
    unittest.main()
