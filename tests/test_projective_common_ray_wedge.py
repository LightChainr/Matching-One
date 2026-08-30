import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from projective_common_ray_wedge import (  # noqa: E402
    analyze,
    build_artifact,
    determinant,
    matrix_rank,
    synthetic_batches,
    validate_artifact,
)


class ProjectiveCommonRayWedgeTest(unittest.TestCase):
    SIZES = (85, 170, 340, 680)

    def test_committed_certificate_reproduces_exactly(self):
        path = ROOT / "analysis" / "projective_common_ray_wedge_certificate.json"
        artifact = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact, build_artifact())
        self.assertEqual(validate_artifact(artifact)["common_wedge_count"], 3)

    def test_common_ray_has_zero_wedges_in_every_replicate(self):
        report = analyze(synthetic_batches((2, 2, 2, 2)), self.SIZES)
        self.assertEqual(report["full_wedges"], {"85": "0", "170": "0", "340": "0"})
        self.assertTrue(
            all(
                all(value == "0" for value in replicate["wedges"].values())
                for replicate in report["delete_one_replicates"]
            )
        )
        self.assertEqual(report["covariance_rank"], 0)

    def test_loading_drift_has_exact_nonzero_wedges(self):
        report = analyze(synthetic_batches((2, "5/2", 3, "7/2")), self.SIZES)
        self.assertEqual(
            report["full_wedges"],
            {"85": "-25/16", "170": "-25/64", "340": "-25/256"},
        )
        self.assertEqual(report["covariance_rank"], 1)
        self.assertTrue(all(report["exact_checks"].values()))

    def test_delete_one_wedges_are_recomputed_not_linearly_adjusted(self):
        report = analyze(synthetic_batches((2, "5/2", 3, "7/2")), self.SIZES)
        first = report["delete_one_replicates"][0]
        last = report["delete_one_replicates"][-1]
        self.assertEqual(first["wedges"]["85"], "-9/4")
        self.assertEqual(last["wedges"]["85"], "-1")
        self.assertNotEqual(first["wedges"], last["wedges"])

    def test_exact_matrix_helpers(self):
        matrix = [[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]]
        self.assertEqual(determinant(matrix), 0)
        self.assertEqual(matrix_rank(matrix), 1)

    def test_float_values_fail_closed(self):
        rows = synthetic_batches((2, 2, 2, 2))
        rows[0]["values"]["85"]["A_M"] = 1.0
        with self.assertRaises(TypeError):
            analyze(rows, self.SIZES)

    def test_missing_generation_and_extra_fields_fail_closed(self):
        rows = synthetic_batches((2, 2, 2, 2))
        del rows[0]["values"]["680"]
        with self.assertRaises(ValueError):
            analyze(rows, self.SIZES)
        rows = synthetic_batches((2, 2, 2, 2))
        rows[0]["values"]["85"]["extra"] = "0"
        with self.assertRaises(ValueError):
            analyze(rows, self.SIZES)

    def test_duplicate_or_empty_batches_fail_closed(self):
        rows = synthetic_batches((2, 2, 2, 2))
        rows[1]["batch"] = rows[0]["batch"]
        with self.assertRaises(ValueError):
            analyze(rows, self.SIZES)
        with self.assertRaises(ValueError):
            analyze(synthetic_batches((2, 2, 2, 2))[:1], self.SIZES)

    def test_non_doubling_lineage_fails_closed(self):
        with self.assertRaises(ValueError):
            analyze(synthetic_batches((2, 2, 2, 2)), (85, 170, 341, 682))

    def test_tampered_certificate_fails_validation(self):
        artifact = build_artifact()
        artifact["loading_drift_control"]["full_wedges"]["85"] = "0"
        with self.assertRaises(ValueError):
            validate_artifact(artifact)


if __name__ == "__main__":
    unittest.main()
