from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from square_bond_kappa3 import (  # noqa: E402
    exact_estimate,
    monte_carlo_estimate,
    square_bond_pairs,
)


class SquareBondGeometryTests(unittest.TestCase):
    def test_bond_and_dual_edge_counts(self) -> None:
        pairs = square_bond_pairs(4)
        self.assertEqual(len(pairs), 32)
        self.assertEqual(
            {pair.primal[2:] for pair in pairs}, {(1, 0), (0, 1)}
        )
        self.assertEqual(
            {pair.dual[2:] for pair in pairs}, {(1, 0), (0, 1)}
        )

    def test_l2_exact_regression(self) -> None:
        result = exact_estimate(2)
        self.assertEqual(result["exact"]["mean_observable"], "0/1")
        self.assertEqual(result["exact"]["first_derivative"], "27/8")
        self.assertEqual(result["exact"]["third_derivative"], "-45/1")
        self.assertEqual(result["exact"]["kappa3"], "-2560/2187")
        self.assertAlmostEqual(result["kappa3"], float(Fraction(-2560, 2187)))

    def test_l3_exact_regression(self) -> None:
        result = exact_estimate(3)
        self.assertEqual(result["exact"]["mean_observable"], "0/1")
        self.assertAlmostEqual(result["kappa3"], -1.4555871991242173, places=14)


class SquareBondMonteCarloTests(unittest.TestCase):
    def test_reproducible_across_worker_counts(self) -> None:
        sequential = monte_carlo_estimate(2, 400, 4, 17, 1)
        parallel = monte_carlo_estimate(2, 400, 4, 17, 2)
        self.assertEqual(sequential, parallel)
        self.assertGreater(sequential["first_derivative_variance"], 0.0)
        self.assertGreater(sequential["third_derivative_variance"], 0.0)
        self.assertGreater(sequential["jackknife"]["standard_error"], 0.0)

    def test_cli_writes_reproducible_batch_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "control.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "square_bond_kappa3.py"),
                    "--mode",
                    "monte-carlo",
                    "--sizes",
                    "2",
                    "--samples",
                    "200",
                    "--blocks",
                    "4",
                    "--workers",
                    "1",
                    "--seed",
                    "23",
                    "--output",
                    str(output),
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(json.loads(completed.stdout), payload)
            self.assertEqual(payload["implemented_models"], ["square_bond_square_torus"])
            self.assertEqual(len(payload["results"][0]["block_seeds"]), 4)
            self.assertIn("triangular_site_rhombic_torus", payload["unimplemented_models"])


if __name__ == "__main__":
    unittest.main()
