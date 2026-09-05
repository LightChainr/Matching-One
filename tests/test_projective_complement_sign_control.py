import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from projective_complement_sign_control import (  # noqa: E402
    analyze_orientation_pair,
    build_artifact,
    canonical_rows,
    certify_complement,
    dualize_rows,
    synthetic_control,
    validate_artifact,
)


class ProjectiveComplementSignControlTest(unittest.TestCase):
    def test_committed_certificate_reproduces_exactly(self):
        path = ROOT / "analysis" / "projective_complement_sign_certificate.json"
        artifact = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(artifact, build_artifact())
        self.assertEqual(validate_artifact(artifact)["orientation_count"], 2)

    def test_row_duality_is_involutive(self):
        control = synthetic_control()
        rows = control["orientations"][0]["rows"]
        dual = dualize_rows(rows, control["line_duality"])
        self.assertEqual(
            canonical_rows(dualize_rows(dual, control["line_duality"])),
            canonical_rows(rows),
        )

    def test_rank_sectors_reverse_and_M_changes_sign(self):
        control = synthetic_control()
        certificate = certify_complement(
            control["orientations"][0]["rows"],
            control["threshold"],
            control["line_duality"],
        )
        self.assertEqual(certificate["original_rank_counts"], {"P0": 1, "P1": 3, "P2": 2})
        self.assertEqual(certificate["dual_rank_counts"], {"P0": 2, "P1": 3, "P2": 1})
        self.assertEqual(certificate["original_M"], "1/6")
        self.assertEqual(certificate["dual_M"], "-1/6")

    def test_ordered_H4_contrast_reverses_sign(self):
        report = analyze_orientation_pair(synthetic_control())
        self.assertEqual(report["H4_contrasts"], {"original": "-1/6", "dual": "1/6"})
        self.assertTrue(all(report["exact_checks"].values()))

    def test_complement_closed_midpoint_is_exactly_odd_zero(self):
        midpoint = analyze_orientation_pair(synthetic_control())["self_matching_midpoint"]
        self.assertEqual(midpoint, {"first": "0", "second": "0", "odd_H4_contrast": "0"})

    def test_float_threshold_and_invalid_birth_fail_closed(self):
        control = synthetic_control()
        control["threshold"] = 0.375
        with self.assertRaises(TypeError):
            analyze_orientation_pair(control)
        control = synthetic_control()
        control["orientations"][0]["rows"][0]["tau2"] = "1/8"
        with self.assertRaises(ValueError):
            analyze_orientation_pair(control)


if __name__ == "__main__":
    unittest.main()
