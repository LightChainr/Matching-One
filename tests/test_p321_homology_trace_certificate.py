from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p321_homology_trace_certificate import (  # noqa: E402
    DEFAULT_CERTIFICATE,
    build_certificate,
    join_adjacent,
    rotate_state,
)


class P321HomologyTraceCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.certificate = build_certificate()
        self.rows = {row["width"]: row for row in self.certificate["widths"]}

    def test_frozen_machine_certificate_matches(self) -> None:
        frozen = json.loads(DEFAULT_CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(self.certificate, frozen)

    def test_periodic_state_actions_are_canonical(self) -> None:
        self.assertEqual(rotate_state((0, 0, 1), 1), (0, 1, 1))
        self.assertEqual(join_adjacent((0, 1, 0), 0), (0, 0, 0))
        self.assertEqual(join_adjacent((0, 1, 1, 2), 3), (0, 1, 1, 0))

    def test_catalan_dimensions_and_translation_order(self) -> None:
        self.assertEqual([self.rows[width]["dimension"] for width in (2, 3, 4)], [2, 5, 14])
        for row in self.rows.values():
            self.assertEqual(row["translation_order_residual"]["nonzero_entries"], 0)

    def test_crossed_seam_pull_through_is_exact(self) -> None:
        for row in self.rows.values():
            residuals = row["local_generator_residuals"]
            self.assertTrue(residuals["all_sites_equal"])
            self.assertEqual(residuals["seam_covariance_residual"]["rank_over_Q"], 0)
            self.assertEqual(residuals["seam_pull_through_residual"]["rank_over_Q"], 0)
            self.assertEqual(
                row["trace_laws"]["crossed_trace_sigma_twisted_max_abs_residual_on_generators"],
                0,
            )

    def test_candidate_is_noncentral_from_width_three(self) -> None:
        self.assertEqual(
            [
                self.rows[width]["local_generator_residuals"]["D_hom_commutator"]["rank_over_Q"]
                for width in (2, 3, 4)
            ],
            [0, 2, 6],
        )
        self.assertEqual(
            [
                self.rows[width]["local_generator_residuals"]["D_hom_pull_through_residual"][
                    "frobenius_norm_squared"
                ]
                for width in (2, 3, 4)
            ],
            [0, 6, 20],
        )

    def test_trace_difference_obeys_neither_single_trace_law(self) -> None:
        self.assertFalse(
            self.rows[2]["trace_laws"]["D_hom_ordinary_cyclicity_witness"]["available"]
        )
        for width in (3, 4):
            witness = self.rows[width]["trace_laws"]["D_hom_ordinary_cyclicity_witness"]
            self.assertTrue(witness["available"])
            self.assertEqual((witness["constant_coefficient"], witness["Q_coefficient"]), (-1, 0))
        self.assertEqual(
            self.rows[3]["trace_laws"]["D_hom_sigma_twisted_witness_A_equals_B_equals_e0"][
                "Q_coefficient"
            ],
            -1,
        )
        self.assertEqual(
            self.rows[4]["trace_laws"]["D_hom_sigma_twisted_witness_A_equals_B_equals_e0"][
                "Q_coefficient"
            ],
            -3,
        )


if __name__ == "__main__":
    unittest.main()
