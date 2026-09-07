#!/usr/bin/env python3
"""Lock the r = 1, 2, 4 aspect ladder at N = 580.

Every assertion names a wrong number that would otherwise reach a frozen
prediction file and then a claim: a torus that does not have the site count it
claims, a rung whose leverage silently differs from the others so one leg is the
noisy one, a spin-8 leakage that compounds in the score instead of cancelling, a
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
    ladder,
    rank_candidates,
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

        A period matrix whose determinant is not the site count would compare
        tori of different sizes, and the finite-size correction would then be
        inside the score with nothing to say so.
        """
        for rung in self.result["rungs"]:
            for key in ("first_matrix_row_major", "second_matrix_row_major"):
                a, b, c, d = rung[key]
                with self.subTest(aspect=rung["aspect_ratio"], matrix=key):
                    self.assertEqual(a * d - b * c, SITE_COUNT)

    def test_all_three_rungs_share_one_leverage(self) -> None:
        """Otherwise one rung is the noisy one and the ladder is not a ladder.

        The leverage divides into every amplitude, so an unequal one would make
        a single sample budget per rung the wrong allocation while the artifact
        reported three comparable numbers.
        """
        self.assertTrue(self.result["leverage_is_shared_across_rungs"])
        self.assertEqual(self.result["shared_leverage"], "8064/4205")
        self.assertEqual({rung["delta_cos4"] for rung in self.result["rungs"]}, {"8064/4205"})

    def test_the_spin8_leakage_cancels_in_the_discriminating_ratio(self) -> None:
        """The systematic that does not shrink with samples.

        Two orientations give A4 + A8 * leakage, not A4. The N=290 pair had
        equal and *opposite* leakage on its two families, so the bias entered
        the score twice. Here the r=1 and r=4 rungs carry the same leakage, so
        it cancels to leading order in the ratio that separates the hypotheses.
        Losing that property would put a systematic back into the one number
        this design exists to measure.
        """
        by_aspect = {rung["aspect_ratio"]: rung for rung in self.result["rungs"]}
        self.assertTrue(self.result["spin8_cancels_in_the_discriminating_ratio"])
        self.assertEqual(
            Fraction(by_aspect[1]["spin8_leakage"]),
            Fraction(by_aspect[4]["spin8_leakage"]),
        )
        self.assertEqual(
            Fraction(by_aspect[2]["spin8_leakage"]),
            -Fraction(by_aspect[1]["spin8_leakage"]),
        )
        self.assertLess(abs(Fraction(by_aspect[1]["spin8_leakage"])), Fraction(1, 15))

    def test_the_search_prefers_this_site_count_over_every_smaller_one(self) -> None:
        """N=580 is not the smallest ladder, so the choice has to be justified.

        100, 200, 260, 340, 400 and 500 all carry it and are cheaper. They lose
        on leverage or on leakage, and if that stopped being true this design
        would be paying for sites it does not need.
        """
        ranked = rank_candidates()
        self.assertEqual(ranked[0]["site_count"], SITE_COUNT)
        cheaper = [row for row in ranked if row["site_count"] < SITE_COUNT]
        self.assertTrue(cheaper, "the search found no cheaper candidate to compare against")
        best_leverage = Fraction(ranked[0]["leverage"])
        best_leakage = Fraction(ranked[0]["max_abs_leakage"])
        for row in cheaper:
            with self.subTest(site_count=row["site_count"]):
                self.assertTrue(
                    Fraction(row["leverage"]) < best_leverage
                    or Fraction(row["max_abs_leakage"]) > best_leakage
                )

    def test_a_rung_needs_two_orientations_of_the_right_norm(self) -> None:
        for aspect in ASPECTS:
            with self.subTest(aspect=aspect):
                self.assertGreaterEqual(len(representations(SITE_COUNT // aspect)), 2)
        self.assertIsNone(ladder(SITE_COUNT + 4), "a site count without the ladder was accepted")

    def test_cos8_is_the_chebyshev_image_of_cos4(self) -> None:
        for a, b in representations(145):
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
            (ROOT / "predictions" / "aspect_ladder_n580_20260905.yaml").read_text(encoding="utf-8")
        )
        self.assertEqual(frozen["design"]["site_count"], SITE_COUNT)
        self.assertEqual(frozen["design_artifact"], "results/aspect-ladder-design/latest.json")

        by_aspect = {rung["aspect_ratio"]: rung for rung in self.result["rungs"]}
        for key, block in frozen["design"]["rungs"].items():
            aspect = int(key[1:])
            with self.subTest(aspect=aspect):
                generated = by_aspect[aspect]
                self.assertEqual(block["gaussian_norm"], generated["gaussian_norm"])
                self.assertEqual(list(block["first_rep"]), generated["first_rep"])
                self.assertEqual(list(block["second_rep"]), generated["second_rep"])
                self.assertEqual(list(block["first_matrix"]), generated["first_matrix_row_major"])
                self.assertEqual(list(block["second_matrix"]), generated["second_matrix_row_major"])
                self.assertEqual(block["delta_cos4"], generated["delta_cos4"])
                self.assertEqual(block["spin8_leakage"], generated["spin8_leakage"])

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
        self.assertIn("no scoring block has been run",
                      self.result["not_established_by_this_design"][0])


if __name__ == "__main__":
    unittest.main()
