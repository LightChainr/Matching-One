from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p48_new_geometry_channels import score  # noqa: E402


class P48NewGeometryChannelTests(unittest.TestCase):
    def test_committed_score_recomputes_without_target_fit(self) -> None:
        source = json.loads(
            (ROOT / "results/server-20260828/P48-retrospective/summary.json").read_text(
                encoding="utf-8"
            )
        )
        targets = {
            n: json.loads(
                (
                    ROOT
                    / f"results/server-20260828/P43-heldout-fullcurve-500m/analysis/n{n}.p48.json"
                ).read_text(encoding="utf-8")
            )
            for n in (185, 265)
        }
        actual = score(source, targets)
        committed = json.loads(
            (
                ROOT / "results/server-20260828/P48-new-geometry-score/score.json"
            ).read_text(encoding="utf-8")
        )

        self.assertEqual(actual["target_refit_parameters"], 0)
        self.assertEqual(actual["target_sizes"], [185, 265])
        self.assertEqual(actual["important_observable_boundary"], committed["important_observable_boundary"])

        expected = {
            "P4_S": 1.1387820575181482,
            "P4_D": 0.2808531455704667,
            "P4_D_prime": 0.08760655046918807,
            "P4_S_prime": 52.71633588357711,
        }
        for channel, value in expected.items():
            self.assertAlmostEqual(actual["channels"][channel]["chi_square"], value, places=10)
            self.assertAlmostEqual(
                actual["channels"][channel]["chi_square"],
                committed["channels"][channel]["chi_square"],
                places=12,
            )

    def test_only_sprime_pure_law_has_large_new_geometry_residual(self) -> None:
        committed = json.loads(
            (
                ROOT / "results/server-20260828/P48-new-geometry-score/score.json"
            ).read_text(encoding="utf-8")
        )
        self.assertLess(committed["channels"]["P4_S"]["chi_square"], 5.0)
        self.assertLess(committed["channels"]["P4_D"]["chi_square"], 5.0)
        self.assertLess(committed["channels"]["P4_D_prime"]["chi_square"], 5.0)
        self.assertGreater(committed["channels"]["P4_S_prime"]["chi_square"], 25.0)


if __name__ == "__main__":
    unittest.main()
