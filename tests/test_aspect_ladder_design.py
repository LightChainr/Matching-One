#!/usr/bin/env python3
"""Lock the r = 1, 2, 4 aspect ladder at N = 1300.

Every assertion here names a wrong number that would otherwise reach a frozen
prediction file and then a claim: a torus that does not have the site count it
claims, a family whose three orientations cannot separate spin 4 from spin 8, a
competitor list that quietly drops the one law the N=290 point happens to sit
on, or an r=4 prediction that is not actually where the competitors separate.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from aspect_ladder_design import (  # noqa: E402
    ASPECTS,
    DEFAULT_OUTPUT,
    SITE_COUNT,
    build_result,
    cos4,
    cos8,
    representations,
    validate_result,
)


class AspectLadderDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()
        cls.committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_committed_artifact_reproduces_exactly(self) -> None:
        self.assertEqual(self.committed, self.result)
        validate_result(self.committed)

    def test_every_torus_really_has_the_site_count(self) -> None:
        """The whole design is 'same N, three moduli'.

        A period matrix whose determinant is not 1300 would compare tori of
        different sizes, and the finite-size correction would then be inside the
        score with nothing to say so.
        """
        for family in self.result["families"]:
            for member in family["members"]:
                a, b, c, d = member["period_matrix_row_major"]
                with self.subTest(member=member["gaussian_integer"]):
                    self.assertEqual(a * d - b * c, SITE_COUNT)
                    self.assertEqual(member["sites"], SITE_COUNT)

    def test_each_family_determines_all_three_of_C_A4_and_A8(self) -> None:
        """Two orientations leave spin 8 as an unremovable systematic.

        That is the caveat the N=290 run had to carry. A singular or
        near-singular design matrix here would put it back while the artifact
        claimed it was fitted out.
        """
        for family in self.result["families"]:
            with self.subTest(aspect=family["aspect_ratio"]):
                self.assertEqual(len(family["members"]), 3)
                determinant = Fraction(family["design_determinant"])
                self.assertNotEqual(determinant, 0)
                self.assertGreater(abs(determinant), Fraction(1, 2))
                values = {member["cos4theta"] for member in family["members"]}
                self.assertEqual(len(values), 3)

    def test_the_three_families_are_equally_conditioned(self) -> None:
        """Otherwise one rung is the noisy one and the ladder is not a ladder.

        This is a fact about 1300, not a design choice, so it is checked rather
        than assumed: it is what makes a single sample budget per family right.
        """
        self.assertTrue(self.result["families_share_one_design_matrix"])
        self.assertTrue(self.result["families_share_one_variance_amplification"])
        amplifications = {f["variance_amplification"]["A4"] for f in self.result["families"]}
        self.assertEqual(len(amplifications), 1)
        # 3-orientation A4 costs 1.64x the 2-orientation variance of the N=290 design
        amplification = Fraction(next(iter(amplifications)))
        self.assertAlmostEqual(float(amplification), 0.8933, places=4)

    def test_the_gaussian_norms_are_the_only_ones_that_work(self) -> None:
        """N = 1300 is claimed to be the smallest site count carrying the ladder.

        If a smaller one existed, this design would be spending compute for
        nothing, so the claim is checked by search rather than by assertion.
        """
        for aspect in ASPECTS:
            self.assertGreaterEqual(len(representations(SITE_COUNT // aspect)), 3, aspect)
        for candidate in range(4, SITE_COUNT, 4):
            with self.subTest(candidate=candidate):
                self.assertFalse(
                    all(len(representations(candidate // aspect)) >= 3 for aspect in ASPECTS),
                    f"{candidate} also carries the ladder and is smaller than {SITE_COUNT}",
                )

    def test_cos8_is_the_chebyshev_image_of_cos4(self) -> None:
        for a, b in representations(325):
            with self.subTest(rep=(a, b)):
                self.assertEqual(cos8(a, b), 2 * cos4(a, b) ** 2 - 1)
                self.assertLessEqual(abs(cos4(a, b)), 1)

    def test_the_post_hoc_law_is_named_as_a_competitor(self) -> None:
        """The reason this file exists.

        `A(r)/A(1) = r` was read off the N=290 point after it came back. Naming
        it in the frozen design is what converts it from a rationalization into
        something that can lose. Dropping it here would let the ladder be
        reported as confirming a law it never risked.
        """
        laws = {row["law"]: row for row in self.result["non_modular_competitors"]}
        self.assertIn("bare_aspect_ratio", laws)
        self.assertEqual(laws["bare_aspect_ratio"]["predicted_by_aspect"]["4"], "4.0")
        self.assertIn("post-hoc", laws["bare_aspect_ratio"]["standing"])
        self.assertIn("prospective", laws["bare_aspect_ratio"]["standing"])

    def test_r4_separates_the_two_live_hypotheses_and_r2_does_not(self) -> None:
        """The reason the ladder goes to 4 and not just to 2.

        If these two numbers were close, the run would be a null by
        construction, and the design would be spending compute to learn nothing.
        """
        shapes = {row["shape"]: row for row in self.result["modular_shape_predictions"]["rows"]}
        weight4 = shapes["E4"]["amplitude_ratio_by_aspect"]
        linear = {
            row["law"]: row["predicted_by_aspect"]
            for row in self.result["non_modular_competitors"]
        }["bare_aspect_ratio"]
        at_two = abs(float(weight4["2"]) - float(linear["2"]))
        at_four = abs(float(weight4["4"]) - float(linear["4"]))
        self.assertLess(at_two, 1.0)
        self.assertGreater(at_four, 6.0)
        self.assertGreater(at_four / at_two, 8.0)

    def test_the_frozen_prediction_file_agrees_with_the_design(self) -> None:
        """The prediction file is what a later run is scored against.

        If it drifts from the design artifact -- a competitor value retyped, a
        period matrix off by a digit -- the run would be scored against numbers
        nothing generated, which is the failure mode the whole generated-tables
        discipline exists to prevent.
        """
        import yaml

        frozen = yaml.safe_load(
            (ROOT / "predictions" / "aspect_ladder_n1300_20260905.yaml").read_text(encoding="utf-8")
        )
        self.assertEqual(frozen["design"]["site_count"], SITE_COUNT)
        self.assertEqual(frozen["design_artifact"], "results/aspect-ladder-design/latest.json")

        by_aspect = {family["aspect_ratio"]: family for family in self.result["families"]}
        for key, block in frozen["design"]["families"].items():
            aspect = int(key[1:])
            with self.subTest(aspect=aspect):
                family = by_aspect[aspect]
                self.assertEqual(block["gaussian_norm"], family["gaussian_norm"])
                self.assertEqual(
                    [list(row) for row in block["period_matrices_row_major"]],
                    [member["period_matrix_row_major"] for member in family["members"]],
                )

        predictions = frozen["competing_predictions"]
        shapes = {row["shape"]: row for row in self.result["modular_shape_predictions"]["rows"]}
        laws = {row["law"]: row for row in self.result["non_modular_competitors"]}
        pairs = [
            ("q4_jordan_weight4", shapes["E4"]["amplitude_ratio_by_aspect"]),
            ("weight8_E8", shapes["E8"]["amplitude_ratio_by_aspect"]),
            ("weight12_E12", shapes["E12"]["amplitude_ratio_by_aspect"]),
            ("weight12_E4_cubed", shapes["E4^3"]["amplitude_ratio_by_aspect"]),
            ("weight12_delta", shapes["Delta"]["amplitude_ratio_by_aspect"]),
            ("bare_aspect_ratio", laws["bare_aspect_ratio"]["predicted_by_aspect"]),
            ("plain_area_scaling", laws["plain_area_scaling"]["predicted_by_aspect"]),
            ("no_modulus_dependence", laws["no_modulus_dependence"]["predicted_by_aspect"]),
        ]
        self.assertEqual(set(predictions), {name for name, _ in pairs})
        for name, generated in pairs:
            with self.subTest(prediction=name):
                self.assertEqual(len(predictions[name]), 2)
                for value, aspect in zip(predictions[name], ("2", "4")):
                    self.assertAlmostEqual(
                        float(value), float(generated[aspect]),
                        delta=abs(float(generated[aspect])) * 1e-9 + 1e-15,
                    )

    def test_the_artifact_does_not_claim_a_measurement(self) -> None:
        self.assertEqual(self.result["status"], "design_only_no_measurement")
        self.assertIn("no block has been run", self.result["not_established_by_this_design"][0])


if __name__ == "__main__":
    unittest.main()
