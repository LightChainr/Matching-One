#!/usr/bin/env python3
"""Lock the five verdicts of the claimed-threshold filter.

Each assertion names a wrong answer the filter would otherwise give about a
claim someone hands us: calling a refuted claim survivable, calling a known
width artifact a discovery, calling a correct-to-ten-digits value a
contradiction of the literature, or matching a polynomial only in the exact form
it happens to be written in.
"""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from threshold_claim_intake import (  # noqa: E402
    EXCLUSIONS,
    decimal_window,
    exclusion_status,
    normalize,
    render,
    verdict,
)


class Arguments:
    def __init__(self, polynomial=None, decimal=None, expression=None):
        self.polynomial = polynomial
        self.decimal = decimal
        self.expression = expression


class NormalizationTests(unittest.TestCase):
    def test_scaling_and_sign_do_not_change_the_polynomial(self) -> None:
        """Otherwise a claim written as -3P would miss its own census entry."""
        self.assertEqual(
            normalize([-102, 54, 42, 267, -3]), [34, -18, -14, -89, 1]
        )


class DecimalWindowTests(unittest.TestCase):
    def test_a_ten_digit_claim_is_a_ten_digit_window(self) -> None:
        low, high = decimal_window("0.5927460508")
        self.assertEqual(high - low, Fraction(1, 10**10))

    def test_a_correct_ten_digit_value_does_not_contradict_the_literature(self) -> None:
        report = render(Arguments(decimal="0.5927460508"))
        self.assertIn(
            "jacobsen-2015-eigenvalue", report["value"]["intervals_containing_it"]
        )

    def test_a_wrong_value_contradicts_all_four(self) -> None:
        report = render(Arguments(decimal="0.5927460509"))
        self.assertEqual(report["value"]["intervals_containing_it"], [])
        self.assertEqual(
            report["verdict"]["outcome"], "contradicts_every_published_interval"
        )


class CensusReadingTests(unittest.TestCase):
    def test_the_quartic_census_is_not_an_exclusion_everywhere(self) -> None:
        """Treating it as one would refute true claims on the wider intervals."""
        quartics = {
            "id": "degree4-height100",
            "artifact": "results/pslq-degree4-{interval}/latest.json",
        }
        self.assertIs(
            exclusion_status(quartics, "jacobsen-2015-eigenvalue")["excluded"], True
        )
        self.assertIs(
            exclusion_status(quartics, "mertens-2022-p-cell")["excluded"], False
        )

    def test_a_recorded_survivor_is_reported_as_one(self) -> None:
        report = render(Arguments(polynomial="34,-18,-14,-89,1"))
        self.assertEqual(
            report["verdict"]["outcome"], "already_catalogued_as_a_width_artifact"
        )


class NoRootTests(unittest.TestCase):
    def test_a_polynomial_with_no_root_in_any_interval_is_refuted(self) -> None:
        """The Ziff A-lattice threshold is exact -- for a different lattice.
        Reporting it as surviving would be the filter's worst failure mode."""
        report = render(Arguments(polynomial="-1,0,2,3,-4,1"))
        self.assertEqual(
            report["verdict"]["outcome"],
            "the_polynomial_has_no_root_in_any_published_interval",
        )


class RefutationBranchTests(unittest.TestCase):
    def test_a_root_inside_an_excluded_class_is_refuted(self) -> None:
        """No such polynomial exists -- that is the theorem -- so this is the
        one branch that can only be reached by a false claim."""
        report = {
            "polynomial": {
                "per_interval": [
                    {
                        "interval": "jacobsen-2015-eigenvalue",
                        "has_root_in_interval": True,
                        "censused_classes_covering_this_claim": [
                            {
                                "class": "degree6-height3",
                                "artifact": "results/x.json",
                                "class_is_excluded_on_this_interval": True,
                                "claim_is_a_recorded_survivor": False,
                                "recorded_survivors_in_this_class": 0,
                            }
                        ],
                    }
                ]
            }
        }
        self.assertEqual(
            verdict(report)["outcome"], "refuted_by_a_committed_certificate"
        )

    def test_the_height_4_census_is_wired_in_and_covers_the_a_lattice_quintic(self) -> None:
        """The wrong answer here is "we have not censused that".

        A height-4 quintic is exactly what the corrected historical range now
        covers.  If the class were missing from EXCLUSIONS the filter would
        report such a claim as uncensused, which is a weaker -- and false --
        statement about what this repository can already say.
        """
        classes = {entry["id"]: entry for entry in EXCLUSIONS}
        self.assertIn("degree6-height4", classes)
        covering = classes["degree6-height4"]
        self.assertEqual((covering["degree_min"], covering["degree_max"], covering["height_max"]), (1, 6, 4))
        for interval in ("jacobsen-2015-eigenvalue", "mertens-2022-p-cell"):
            with self.subTest(interval=interval):
                self.assertTrue(exclusion_status(covering, interval)["excluded"])

        a_lattice = [-1, 0, 2, 3, -4, 1]  # the Ziff A-lattice quintic, height 4
        report = render(Arguments(polynomial=",".join(str(v) for v in a_lattice)))
        self.assertEqual(report["polynomial"]["height"], 4)
        for row in report["polynomial"]["per_interval"]:
            covering_ids = {entry["class"] for entry in row["censused_classes_covering_this_claim"]}
            self.assertIn("degree6-height4", covering_ids)
            self.assertNotIn("degree4-height100", covering_ids)  # degree 5 is outside it

    def test_the_filter_never_confirms(self) -> None:
        report = render(Arguments(decimal="0.59274605079210"))
        self.assertEqual(report["verdict"]["outcome"], "survives_our_checks")
        self.assertIn(
            "not a confirmation", report["verdict"]["what_this_means"]
        )


if __name__ == "__main__":
    unittest.main()
