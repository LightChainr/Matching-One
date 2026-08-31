import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from projective_birth_crosswalk import (  # noqa: E402
    BirthRow,
    build_artifact,
    reconstruct_at_threshold,
    synthetic_rows,
    validate_artifact,
)


class ProjectiveBirthCrosswalkTest(unittest.TestCase):
    def test_committed_certificate_reproduces_exactly(self):
        path = ROOT / "analysis" / "projective_birth_crosswalk_certificate.json"
        artifact = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact, build_artifact())
        self.assertEqual(
            validate_artifact(artifact),
            {
                "schema": "matching-one/projective-birth-crosswalk/v1",
                "status": "valid_exact_same_stream_birth_crosswalk_control",
                "row_count": 6,
                "threshold_count": 5,
            },
        )

    def test_all_frozen_threshold_reconstructions(self):
        expected = {
            "0": ((5, 1, 0), (1, 0), {"L0": 1}, "-5/6"),
            "1/4": ((4, 2, 0), (2, 0), {"L0": 2}, "-2/3"),
            "1/2": ((2, 2, 2), (4, 2), {"L0": 1, "L1": 1}, "0"),
            "3/4": ((0, 2, 4), (6, 4), {"L1": 1, "L2": 1}, "2/3"),
            "1": ((0, 0, 6), (6, 6), {}, "1"),
        }
        for threshold, (rank, cdf, lines, matching) in expected.items():
            with self.subTest(threshold=threshold):
                result = reconstruct_at_threshold(synthetic_rows(), threshold)
                self.assertEqual(tuple(result["rank_counts"].values()), rank)
                self.assertEqual(tuple(result["birth_cdf_counts"].values()), cdf)
                self.assertEqual(result["plateau_line_counts"], lines)
                self.assertEqual(result["probabilities"]["M"], matching)
                self.assertTrue(all(result["exact_identities"].values()))

    def test_birth_boundary_is_inclusive(self):
        direct = BirthRow(Fraction(1, 2), Fraction(1, 2), "direct_rank2", None)
        plateau = BirthRow(Fraction(1, 4), Fraction(3, 4), "plateau", "L")
        self.assertEqual(direct.rank_at(Fraction(1, 2)), 2)
        self.assertEqual(plateau.rank_at(Fraction(1, 4)), 1)
        self.assertEqual(plateau.rank_at(Fraction(3, 4)), 2)

    def test_kind_inventory_is_independent_of_threshold(self):
        for threshold in ("0", "1/2", "1"):
            result = reconstruct_at_threshold(synthetic_rows(), threshold)
            self.assertEqual(
                result["input_kind_counts"],
                {"DIRECT_RANK2": 2, "plateau": 4},
            )

    def test_record_schema_fails_closed(self):
        row = synthetic_rows()[0]
        for malformed in (
            {key: value for key, value in row.items() if key != "line"},
            {**row, "extra": 1},
        ):
            with self.subTest(malformed=malformed):
                with self.assertRaises(ValueError):
                    reconstruct_at_threshold([malformed], "1/2")

    def test_inexact_numeric_values_fail_closed(self):
        with self.assertRaises(TypeError):
            reconstruct_at_threshold(
                [{"tau1": 0.25, "tau2": "3/4", "kind": "plateau", "line": "L"}],
                "1/2",
            )
        with self.assertRaises(TypeError):
            reconstruct_at_threshold(synthetic_rows(), 0.5)

    def test_invalid_birth_order_and_range_fail_closed(self):
        for tau1, tau2 in (("3/4", "1/4"), ("-1/4", "1/2"), ("1/2", "5/4")):
            with self.subTest(tau1=tau1, tau2=tau2):
                with self.assertRaises(ValueError):
                    reconstruct_at_threshold(
                        [{"tau1": tau1, "tau2": tau2, "kind": "plateau", "line": "L"}],
                        "1/2",
                    )

    def test_direct_and_plateau_semantics_fail_closed(self):
        malformed = (
            {"tau1": "1/4", "tau2": "1/2", "kind": "direct_rank2", "line": None},
            {"tau1": "1/2", "tau2": "1/2", "kind": "direct_rank2", "line": "L"},
            {"tau1": "1/2", "tau2": "1/2", "kind": "plateau", "line": "L"},
            {"tau1": "1/4", "tau2": "1/2", "kind": "plateau", "line": None},
        )
        for row in malformed:
            with self.subTest(row=row):
                with self.assertRaises(ValueError):
                    reconstruct_at_threshold([row], "1/2")

    def test_empty_input_and_out_of_range_threshold_fail_closed(self):
        with self.assertRaises(ValueError):
            reconstruct_at_threshold([], "1/2")
        for threshold in ("-1/4", "5/4"):
            with self.assertRaises(ValueError):
                reconstruct_at_threshold(synthetic_rows(), threshold)

    def test_tampered_certificate_fails_validation(self):
        artifact = build_artifact()
        artifact["threshold_reconstructions"][2]["rank_counts"]["P1"] = 3
        with self.assertRaises(ValueError):
            validate_artifact(artifact)


if __name__ == "__main__":
    unittest.main()
