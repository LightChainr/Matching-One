#!/usr/bin/env python3
"""Tests for typed generic-Q lift descriptors and tangent transport."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import generic_q_lift_semantics as semantics  # noqa: E402


PHASE_A = json.loads(semantics.DEFAULT_PHASE_A.read_text(encoding="utf-8"))


class GenericQLiftSemanticsTests(unittest.TestCase):
    def test_descriptors_expose_every_required_coordinate(self) -> None:
        for lift in ("L_hom", "L_CP"):
            value = semantics.descriptor(lift).to_dict()
            for key in (
                "endpoint_observable_id", "lift_id", "sector_weights_in_Q",
                "normalization", "Q_v_path", "explicit_insertion_Q_dependence",
                "projector_convention", "counterterm_convention",
                "field_normalization_convention",
            ):
                self.assertIn(key, value)
        hom = semantics.descriptor("L_hom").to_dict()
        cp = semantics.descriptor("L_CP").to_dict()
        self.assertEqual(hom["sector_weights_in_Q"]["0D"], "-1")
        self.assertEqual(cp["sector_weights_in_Q"]["0D"], "-Q")
        self.assertEqual(
            semantics.GenericQLiftDescriptor.from_dict(cp),
            semantics.descriptor("L_CP"),
        )
        cp["sector_weights_in_Q"]["0D"] = "-1"
        with self.assertRaisesRegex(ValueError, "descriptor drift"):
            semantics.GenericQLiftDescriptor.from_dict(cp)

    def test_same_endpoint_different_lift_is_not_raw_comparable(self) -> None:
        hom = semantics.descriptor("L_hom")
        cp = semantics.descriptor("L_CP")
        self.assertFalse(semantics.raw_tangents_directly_comparable(hom, cp))
        with self.assertRaisesRegex(ValueError, "raw Q tangents are not directly comparable"):
            semantics.require_raw_tangent_identity(hom, cp)
        semantics.require_raw_tangent_identity(hom, hom)

    def test_exact_normalized_transport_is_minus_pi0(self) -> None:
        endpoint = {"pi_2D": Fraction(69, 256), "pi_1D": Fraction(118, 256), "pi_0D": Fraction(69, 256)}
        result = semantics.transport_tangent(
            semantics.descriptor("L_hom", "fixed_v_1"),
            semantics.descriptor("L_CP", "fixed_v_1"),
            endpoint_sectors=endpoint,
        )
        self.assertEqual(result["tangent_shift_target_minus_source"], "-pi_0D")
        self.assertEqual(result["numeric_shift"], "-69/256")
        reverse = semantics.transport_tangent(
            semantics.descriptor("L_CP", "fixed_v_1"),
            semantics.descriptor("L_hom", "fixed_v_1"),
            endpoint_sectors=endpoint,
        )
        self.assertEqual(reverse["numeric_shift"], "69/256")

    def test_transport_reproduces_both_phase_A_paths_and_sizes(self) -> None:
        manifest = semantics.build_manifest(PHASE_A)
        self.assertTrue(manifest["all_checks_passed"])
        self.assertEqual(len(manifest["phase_A_transport_checks"]), 4)
        for row in manifest["phase_A_transport_checks"]:
            self.assertEqual(row["transport"]["numeric_shift"], "-" + row["pi_0D"])
            self.assertEqual(row["transport"]["transported_tangent"], row["raw_d_L_CP"])

    def test_cp_horizontal_critical_transport_annihilates_homology_tangent(self) -> None:
        row = PHASE_A["finite_tori"][1]
        pi0 = Fraction(row["Q1_v1_sector_counts"]["pi_0D"]["fraction"])
        d_h = Fraction(row["paths"]["critical_square_bond_v_sqrt_Q"]["normalized_dh"]["fraction"])
        result = semantics.to_cp_horizontal(
            semantics.descriptor("L_hom"),
            endpoint_sectors={"pi_2D": pi0, "pi_1D": 1 - 2 * pi0, "pi_0D": pi0},
            tangent=d_h,
        )
        self.assertEqual(result["transported_tangent"], "0")
        self.assertEqual(result["connection"], "CP-horizontal")

    def test_path_or_normalization_drift_fails_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "different .* paths"):
            semantics.transport_tangent(
                semantics.descriptor("L_hom", "fixed_v_1"),
                semantics.descriptor("L_CP", "critical_square_bond_v_sqrt_Q"),
            )
        with self.assertRaisesRegex(ValueError, "different normalizations"):
            semantics.transport_tangent(
                semantics.descriptor("L_hom", normalization="normalized_probability"),
                semantics.descriptor("L_CP", normalization="restricted_state_sum"),
            )


if __name__ == "__main__":
    unittest.main()
