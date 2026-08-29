from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_triangular_energy_logpair_mc import (  # noqa: E402
    analyze_archives,
    delta_radius,
    read_archives,
    run_batches,
    write_archives,
)


class TriangularEnergyLogpairMonteCarloTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.serial = run_batches(
            length=32,
            denominator=8,
            samples=8,
            batches=4,
            seed=2026082934,
            workers=1,
        )

    def test_frozen_delta_radii(self) -> None:
        self.assertEqual(delta_radius(64, 8), 6)
        self.assertEqual(delta_radius(64, 12), 4)
        self.assertEqual(delta_radius(64, 16), 3)
        with self.assertRaisesRegex(ValueError, "even"):
            delta_radius(63, 8)

    def test_worker_count_does_not_change_streams(self) -> None:
        parallel = run_batches(
            length=32,
            denominator=8,
            samples=8,
            batches=4,
            seed=2026082934,
            workers=2,
        )
        self.assertEqual(self.serial, parallel)

    def test_archive_round_trip_is_exact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "batches.csv"
            write_archives(path, self.serial)
            self.assertEqual(read_archives(path), self.serial)

    def test_analysis_exposes_full_covariance_and_no_sign_draws(self) -> None:
        payload = analyze_archives(self.serial, length=32, denominator=8)
        self.assertEqual(payload["moment_order"], ["LL", "LD", "DD"])
        self.assertEqual(len(payload["covariance_of_mean"]), 3)
        self.assertTrue(all(len(row) == 3 for row in payload["covariance_of_mean"]))
        self.assertEqual(payload["monte_carlo"]["random_cluster_sign_draws"], 0)
        self.assertEqual(
            payload["monte_carlo"]["centering"],
            "unbiased cross-configuration order-2 U-statistic",
        )

    def test_committed_smoke_is_finite_and_labeled(self) -> None:
        path = ROOT / "results/local-20260829/P234-triangular-energy-logpair-smoke/L32-d8-40.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(payload["status"], "raw_parent_pair_phaseA_no_universal_coefficient")
        self.assertEqual(payload["monte_carlo"]["total_samples"], 40)
        self.assertEqual(payload["geometry"]["bilocal_radius_lattice_units"], 3)
        self.assertIn("rescaling_invariant_J", payload["jordan_diagnostics"])
        self.assertTrue(
            all(value == value for value in payload["estimate"])
        )


if __name__ == "__main__":
    unittest.main()
