from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_threshold_rank_covariance import (  # noqa: E402
    _orientation_batches,
    bootstrap_heldout_calibration,
    chi_square_survival,
    constant_heldout_audit,
    covariance_of_mean,
    hotelling_calibration,
    jackknife_pseudovalues,
    stable_inverse,
    validate_coupling_metadata,
)
from validate_threshold_rank_covariance_archive import validate_archive  # noqa: E402


ARCHIVE = ROOT / "results" / "server-20260828" / "P33-cross-size-covariance"


def _record(n: int, orientation: str, batch: int, samples: int = 10) -> dict:
    return {
        "n": n,
        "orientation": orientation,
        "batch": batch,
        "samples": samples,
    }


def _aligned_records(sizes=(65, 85), batches=(0, 1), samples=10) -> dict:
    records = {}
    for n in sizes:
        for orientation in ("first", "second"):
            for batch in batches:
                records[(n, orientation, batch)] = _record(
                    n, orientation, batch, samples=samples
                )
    return records


class ThresholdRankCrossSizeCovarianceTests(unittest.TestCase):
    def test_covariance_is_for_the_mean_not_raw_batches(self) -> None:
        means, covariance = covariance_of_mean(
            [
                [1.0, 2.0],
                [2.0, 4.0],
                [3.0, 6.0],
                [4.0, 8.0],
            ]
        )
        self.assertEqual(means, [2.5, 5.0])
        self.assertAlmostEqual(covariance[0][0], 5.0 / 12.0)
        self.assertAlmostEqual(covariance[0][1], 5.0 / 6.0)
        self.assertAlmostEqual(covariance[1][1], 5.0 / 3.0)

    def test_jackknife_pseudovalues_preserve_full_center(self) -> None:
        pseudo = jackknife_pseudovalues(10.0, [9.0, 10.0, 11.0])
        self.assertEqual(pseudo, [12.0, 10.0, 8.0])
        self.assertAlmostEqual(sum(pseudo) / len(pseudo), 10.0)

    def test_constant_heldout_score_propagates_training_covariance(self) -> None:
        sizes = [65, 85, 130, 145, 170]
        values = [0.80, 0.82, 0.79, 0.81, 0.77]
        variances = [0.01**2, 0.012**2, 0.015**2, 0.02**2, 0.025**2]
        covariance = [
            [
                variances[i]
                if i == j
                else 0.1 * math.sqrt(variances[i] * variances[j])
                for j in range(len(sizes))
            ]
            for i in range(len(sizes))
        ]
        result = constant_heldout_audit(
            values,
            covariance,
            sizes,
            [65, 85, 130],
            [145, 170],
            batch_count=100,
        )
        self.assertEqual(result["heldout_dof"], 2)
        self.assertEqual(result["heldout_chi_square_kind"], "plugin_asymptotic")
        self.assertAlmostEqual(sum(result["training_weights"]), 1.0)
        self.assertGreater(result["amplitude_se"], 0.0)
        self.assertTrue(math.isfinite(result["heldout_chi_square"]))
        residual_covariance = result["heldout_residual_covariance"]
        self.assertNotEqual(residual_covariance[0][1], 0.0)
        self.assertEqual(result["score_solver"], "cholesky")
        self.assertIn("eigenvalues", result["heldout_residual_eigenstructure"])
        self.assertIn("hotelling", result["finite_batch_calibration"])

    def test_alignment_requires_same_batches_for_every_size(self) -> None:
        records = _aligned_records()
        sizes, batches, _grouped = _orientation_batches(records)
        self.assertEqual(sizes, [65, 85])
        self.assertEqual(batches, [0, 1])

        del records[(85, "second", 1)]
        with self.assertRaisesRegex(ValueError, "aligned batch ids"):
            _orientation_batches(records)

    def test_equal_batch_weight_rejects_unequal_cross_size_samples(self) -> None:
        records = _aligned_records()
        records[(85, "first", 0)]["samples"] = 9
        records[(85, "second", 0)]["samples"] = 9
        with self.assertRaisesRegex(ValueError, "one sample count"):
            _orientation_batches(records)

    def test_paired_orientations_still_require_equal_samples(self) -> None:
        records = _aligned_records()
        records[(65, "second", 1)]["samples"] = 11
        with self.assertRaisesRegex(ValueError, "paired orientations"):
            _orientation_batches(records)

    def test_diagonal_covariance_full_and_diagonal_scores_match(self) -> None:
        sizes = [65, 85, 130, 145, 170]
        values = [0.0, 0.0, 0.0, 1.0, 1.0]
        covariance = [
            [1.0 if i == j else 0.0 for j in range(len(sizes))]
            for i in range(len(sizes))
        ]
        full = constant_heldout_audit(
            values, covariance, sizes, [65, 85, 130], [145, 170]
        )
        diagonal = constant_heldout_audit(
            values,
            [
                [covariance[i][i] if i == j else 0.0 for j in range(len(sizes))]
                for i in range(len(sizes))
            ],
            sizes,
            [65, 85, 130],
            [145, 170],
        )
        self.assertAlmostEqual(full["amplitude"], 0.0)
        self.assertAlmostEqual(full["heldout_chi_square"], 1.2)
        self.assertEqual(full["heldout_chi_square"], diagonal["heldout_chi_square"])
        self.assertEqual(full["amplitude"], diagonal["amplitude"])
        self.assertEqual(
            full["heldout_residual_covariance"],
            diagonal["heldout_residual_covariance"],
        )

    def test_analytic_equicorrelated_off_diagonal_scores(self) -> None:
        """GLS constant fit with analytic full vs diagonal answers.

        Sizes 65,85 train and 130,145 hold out; observations [0,0,1,1];
        equicorrelation 1/2. Training GLS weights are 1/2,1/2 in both modes.

        Full residual covariance is [[3/4, 1/4], [1/4, 3/4]], chi-square = 2.
        Diagonal residual covariance is [[3/2, 1/2], [1/2, 3/2]], chi-square = 1.
        """

        sizes = [65, 85, 130, 145]
        values = [0.0, 0.0, 1.0, 1.0]
        covariance = [
            [1.0 if i == j else 0.5 for j in range(len(sizes))]
            for i in range(len(sizes))
        ]
        full = constant_heldout_audit(
            values, covariance, sizes, [65, 85], [130, 145]
        )
        diagonal = constant_heldout_audit(
            values,
            [
                [covariance[i][i] if i == j else 0.0 for j in range(len(sizes))]
                for i in range(len(sizes))
            ],
            sizes,
            [65, 85],
            [130, 145],
        )
        self.assertEqual(full["training_weights"], [0.5, 0.5])
        self.assertEqual(diagonal["training_weights"], [0.5, 0.5])
        self.assertAlmostEqual(full["amplitude"], 0.0)
        self.assertAlmostEqual(diagonal["amplitude"], 0.0)
        self.assertAlmostEqual(full["amplitude_se"] ** 2, 0.75)
        self.assertAlmostEqual(diagonal["amplitude_se"] ** 2, 0.5)
        self.assertAlmostEqual(full["heldout_residual_covariance"][0][0], 0.75)
        self.assertAlmostEqual(full["heldout_residual_covariance"][0][1], 0.25)
        self.assertAlmostEqual(diagonal["heldout_residual_covariance"][0][0], 1.5)
        self.assertAlmostEqual(diagonal["heldout_residual_covariance"][0][1], 0.5)
        self.assertAlmostEqual(full["heldout_chi_square"], 2.0)
        self.assertAlmostEqual(diagonal["heldout_chi_square"], 1.0)
        self.assertNotEqual(full["heldout_chi_square"], diagonal["heldout_chi_square"])

    def test_analytic_three_size_off_diagonal_scores(self) -> None:
        sizes = [65, 85, 130]
        values = [1.0, 1.0, 2.0]
        covariance = [
            [1.0 if i == j else 0.5 for j in range(len(sizes))]
            for i in range(len(sizes))
        ]
        full = constant_heldout_audit(
            values, covariance, sizes, [65, 85], [130]
        )
        diagonal = constant_heldout_audit(
            values,
            [
                [1.0 if i == j else 0.0 for j in range(len(sizes))]
                for i in range(len(sizes))
            ],
            sizes,
            [65, 85],
            [130],
        )
        self.assertAlmostEqual(full["heldout_chi_square"], 4.0 / 3.0)
        self.assertAlmostEqual(diagonal["heldout_chi_square"], 2.0 / 3.0)

    def test_stable_inverse_reports_rank_for_singular_covariance(self) -> None:
        inverse, diagnostics = stable_inverse([[1.0, 1.0], [1.0, 1.0]])
        self.assertEqual(diagnostics["solver"], "svd_pseudoinverse")
        self.assertEqual(diagnostics["effective_rank"], 1)
        self.assertAlmostEqual(inverse[0][0], 0.25)
        self.assertAlmostEqual(inverse[0][1], 0.25)
        self.assertAlmostEqual(inverse[1][1], 0.25)
        self.assertEqual(len(diagnostics["eigenvalues"]), 2)
        self.assertAlmostEqual(diagnostics["eigenvalues"][0], 0.0, places=10)
        self.assertAlmostEqual(diagnostics["eigenvalues"][1], 2.0, places=10)

    def test_near_singular_training_uses_svd_pseudoinverse(self) -> None:
        sizes = [65, 85, 130]
        values = [1.0, 1.0, 1.0]
        covariance = [
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
        result = constant_heldout_audit(
            values, covariance, sizes, [65, 85], [130]
        )
        self.assertEqual(
            result["training_covariance_eigenstructure"]["solver"],
            "svd_pseudoinverse",
        )
        self.assertEqual(
            result["training_covariance_eigenstructure"]["effective_rank"], 1
        )
        self.assertAlmostEqual(result["amplitude"], 1.0)
        self.assertAlmostEqual(result["heldout_chi_square"], 0.0)
        self.assertIsNotNone(result["eigenvalue_truncation_sensitivity"])
        self.assertGreaterEqual(
            len(result["eigenvalue_truncation_sensitivity"]), 4
        )

    def test_hotelling_and_chi_square_survival_match_known_values(self) -> None:
        self.assertAlmostEqual(chi_square_survival(0.0, 2), 1.0)
        self.assertAlmostEqual(
            chi_square_survival(5.53, 2), math.exp(-5.53 / 2.0), places=12
        )
        calibrated = hotelling_calibration(2.0, 2, 100)
        self.assertTrue(calibrated["applicable"])
        self.assertAlmostEqual(calibrated["f_statistic"], 98.0 / 99.0)
        self.assertEqual(calibrated["df"], [2, 98])
        self.assertTrue(0.0 < calibrated["f_survival"] < 1.0)

    def test_bootstrap_calibration_returns_ordered_quantiles(self) -> None:
        vectors = [
            [0.0, 0.0, 0.0, 0.1],
            [0.1, -0.1, 0.0, 0.0],
            [-0.1, 0.2, 0.1, -0.1],
            [0.0, 0.0, -0.1, 0.2],
            [0.2, 0.1, 0.0, 0.0],
            [-0.2, -0.1, 0.1, 0.1],
        ]
        observed = constant_heldout_audit(
            *covariance_of_mean(vectors),
            [65, 85, 130, 145],
            [65, 85],
            [130, 145],
            batch_count=len(vectors),
        )
        result = bootstrap_heldout_calibration(
            vectors,
            [65, 85, 130, 145],
            [65, 85],
            [130, 145],
            float(observed["heldout_chi_square"]),
            replicates=40,
            seed=20260828,
        )
        quantiles = result["quantiles"]
        self.assertLessEqual(quantiles["0.5"], quantiles["0.9"])
        self.assertLessEqual(quantiles["0.9"], quantiles["0.95"])
        self.assertLessEqual(quantiles["0.95"], quantiles["0.99"])
        self.assertEqual(result["replicates_used"] + result["failures"], 40)
        self.assertGreaterEqual(result["observed_survival"], 0.0)
        self.assertLessEqual(result["observed_survival"], 1.0)

    def test_bootstrap_diagonal_flag_changes_the_resampled_scores(self) -> None:
        vectors = [
            [0.0, 0.4, 0.0, 1.2],
            [0.3, 0.0, 0.4, 0.8],
            [-0.2, 0.5, -0.1, 1.5],
            [0.1, -0.3, 0.2, 0.9],
            [0.4, 0.2, 0.3, 1.1],
            [-0.3, 0.1, -0.2, 1.3],
            [0.2, 0.3, 0.1, 0.7],
            [-0.1, -0.2, 0.5, 1.4],
        ]
        sizes = [65, 85, 130, 145]
        kwargs = dict(
            batch_vectors=vectors,
            sizes=sizes,
            training_sizes=[65, 85],
            heldout_sizes=[130, 145],
            observed_chi_square=1.0,
            replicates=40,
            seed=20260828,
        )
        full = bootstrap_heldout_calibration(diagonal_only=False, **kwargs)
        diagonal = bootstrap_heldout_calibration(diagonal_only=True, **kwargs)
        self.assertFalse(full["diagonal_only"])
        self.assertTrue(diagonal["diagonal_only"])
        self.assertNotEqual(full["quantiles"]["0.5"], diagonal["quantiles"]["0.5"])

    def test_coupling_metadata_requires_shared_rng_seed_and_counters(self) -> None:
        records = _aligned_records(samples=100)
        metadata = {
            "engine": "same-N Gaussian threshold-rank Newman-Ziff",
            "rng": "counter-derived SplitMix64 stream plus unbiased Fisher-Yates",
            "seed": 2026093303,
            "replica_counter_first": 4000000000,
            "replica_counter_last_exclusive": 4000000200,
            "batches": 2,
            "samples_per_pair": 200,
            "coupling": "same cyclic permutation shared by same-N orientations",
            "designs": [{"N": 65}, {"N": 85}],
        }
        contract = validate_coupling_metadata(
            metadata, records, "p33-working-tree-10m", expected_seed=2026093303
        )
        self.assertTrue(contract["validated"])
        self.assertEqual(contract["samples_per_batch"], 100)
        self.assertEqual(
            contract["per_batch_counter_ranges"][1]["replica_counter_first"],
            4000000100,
        )

        bad_seed = dict(metadata)
        bad_seed["seed"] = 1
        with self.assertRaisesRegex(ValueError, "disagrees with declared seed"):
            validate_coupling_metadata(
                bad_seed, records, "p33-working-tree-10m", expected_seed=2026093303
            )

        bad_interval = dict(metadata)
        bad_interval["replica_counter_last_exclusive"] = 4000000199
        with self.assertRaisesRegex(ValueError, "counter interval"):
            validate_coupling_metadata(bad_interval, records, "label")

        override = dict(metadata)
        override["designs"] = [{"N": 65, "seed": 1}, {"N": 85}]
        with self.assertRaisesRegex(ValueError, "per-size RNG overrides"):
            validate_coupling_metadata(override, records, "label")

        with self.assertRaisesRegex(ValueError, "seed-label"):
            validate_coupling_metadata(metadata, records, "   ")

    def test_committed_archive_has_stable_design_and_positive_covariance(self) -> None:
        result = validate_archive(
            ARCHIVE / "batch_metrics.csv",
            ARCHIVE / "summary.json",
            1e12,
        )
        batch = result["batch_contract"]
        self.assertEqual(batch["sizes"], [65, 85, 130, 145, 170])
        self.assertEqual(batch["batch_count"], 100)
        self.assertEqual(batch["samples_per_size_batch"], 100000)
        root = result["summary_contract"]["metric_covariance_diagnostics"][
            "root_gap"
        ]
        self.assertLess(root["infinity_norm_condition"], 10.0)


if __name__ == "__main__":
    unittest.main()
