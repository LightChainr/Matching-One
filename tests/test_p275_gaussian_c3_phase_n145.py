from __future__ import annotations

import importlib.util
from math import cos, sin
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT = ROOT / "experiments" / "p275-gaussian-c3-phase-n145-20260901"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


producer = load("p275_gaussian_c3_n145_producer", EXPERIMENT / "producer.py")
scorer = load("p275_gaussian_c3_n145_score", EXPERIMENT / "score.py")


class N145GeometryTests(unittest.TestCase):
    def test_exact_new_angle(self) -> None:
        producer.validate_geometries()
        self.assertAlmostEqual(cos(scorer.DELTA), 4 / 5, places=15)
        self.assertAlmostEqual(sin(scorer.DELTA), 3 / 5, places=15)
        self.assertEqual(producer.PRODUCTION_SAMPLES, 5_000_000)
        self.assertEqual(producer.PRODUCTION_BATCHES, 100)
        self.assertEqual(producer.PRODUCTION_SEED, 20_260_901_277)

    def test_invariant_predictions_are_distinct(self) -> None:
        self.assertAlmostEqual(scorer.INVARIANT_PREDICTIONS["H0"], 1.0)
        self.assertAlmostEqual(scorer.INVARIANT_PREDICTIONS["H4"], cos(8 * scorer.DELTA))
        self.assertAlmostEqual(scorer.INVARIANT_PREDICTIONS["H8"], cos(16 * scorer.DELTA))


class N145DecisionTests(unittest.TestCase):
    def synthetic_rows(self, phase: float) -> list[list[float]]:
        a = [0.00120, -0.000894]
        rot = scorer.BASE.rotation(phase)
        z2 = [0.7366 * sum(rot[i][j] * a[j] for j in range(2)) for i in range(2)]
        offsets = [
            [-2, 1, 1, -1], [-1, -2, 2, 1], [0, 2, -2, 2], [1, -1, 1, -2],
            [2, 0, -1, 1], [-2, -1, 2, -2], [-1, 2, -1, 2], [0, -2, 2, -1],
            [1, 1, -2, 1], [2, -1, 1, -2],
        ]
        center = a + z2
        return [[center[j] + 1e-5 * offset[j] for j in range(4)] for offset in offsets]

    def test_each_unique_model_has_a_distinct_decision(self) -> None:
        expectations = {
            "H0": (0.0, "H0_EVEN_CHARACTER_SELECTED"),
            "H4": (+4 * scorer.DELTA, "H4_SELECTED"),
            "H8": (-8 * scorer.DELTA, "H8_OBSERVER_HARMONIC_SELECTED"),
        }
        for model, (phase, expected_decision) in expectations.items():
            with self.subTest(model=model):
                result = scorer.score(self.synthetic_rows(phase))
                self.assertEqual(result["decision"], expected_decision)
                self.assertEqual(result["surviving_models"], [model])
                self.assertAlmostEqual(
                    result["projective_invariant"]["value"],
                    scorer.INVARIANT_PREDICTIONS[model],
                    delta=0.02,
                )


if __name__ == "__main__":
    unittest.main()
