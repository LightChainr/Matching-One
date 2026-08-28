from __future__ import annotations

from collections import Counter
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_semigroup_design import (  # noqa: E402
    Gaussian,
    default_catalog,
    lineage_payload,
)


def _ratio(prediction):
    payload = prediction["angular_ratio"]
    return Fraction(payload["numerator"], payload["denominator"])


class GaussianSemigroupDesignTests(unittest.TestCase):
    def test_exact_harmonics(self) -> None:
        value = Gaussian(8, 1)
        self.assertEqual(value.cos4m(1), Fraction(3713, 4225))
        self.assertEqual(value.cos4m(2), 2 * value.cos4() ** 2 - 1)
        self.assertEqual(
            value.cos4m(3), 4 * value.cos4() ** 3 - 3 * value.cos4()
        )

    def test_smith_invariants_distinguish_equal_area_translation_groups(self) -> None:
        self.assertEqual(Gaussian(13, 0).smith_invariants(), (13, 13))
        self.assertFalse(Gaussian(13, 0).cyclic_translation_group)
        self.assertEqual(Gaussian(12, 5).smith_invariants(), (1, 169))
        self.assertTrue(Gaussian(12, 5).cyclic_translation_group)
        with self.assertRaisesRegex(ValueError, "different finite translation groups"):
            lineage_payload(Gaussian(13, 0), Gaussian(12, 5), Gaussian(1, 1))

    def test_doubling_reverses_H4_delta(self) -> None:
        payload = lineage_payload(
            Gaussian(8, 1), Gaussian(7, 4), Gaussian(1, 1)
        )
        prediction = payload["harmonic_predictions"]["H4"]
        self.assertEqual(prediction["parent_delta"]["numerator"], 1152)
        self.assertEqual(prediction["parent_delta"]["denominator"], 845)
        self.assertEqual(prediction["child_delta"]["numerator"], -1152)
        self.assertEqual(prediction["child_delta"]["denominator"], 845)
        self.assertEqual(prediction["angular_ratio"]["numerator"], -1)
        self.assertEqual(prediction["angular_ratio"]["denominator"], 1)
        self.assertEqual(prediction["target_expression"], "(-1/1)*2^(-13/8)")
        self.assertEqual(payload["radial_factor_expression"], "2^(-13/8)")
        self.assertEqual(
            payload["child"]["first_canonical"]["pair"], [9, 7]
        )
        self.assertEqual(
            payload["child"]["second_canonical"]["pair"], [11, 3]
        )
        self.assertEqual(
            payload["pair_translation_group_contract"]["parent_smith_invariants"],
            [1, 65],
        )
        self.assertEqual(
            payload["pair_translation_group_contract"]["child_smith_invariants"],
            [1, 130],
        )

    def test_third_doubling_lineage(self) -> None:
        payload = lineage_payload(
            Gaussian(12, 1), Gaussian(9, 8), Gaussian(1, 1)
        )
        self.assertEqual(payload["child"]["first_canonical"]["pair"], [13, 11])
        self.assertEqual(payload["child"]["second_canonical"]["pair"], [17, 1])
        prediction = payload["harmonic_predictions"]["H4"]
        self.assertEqual(prediction["angular_ratio"]["numerator"], -1)
        self.assertEqual(prediction["angular_ratio"]["denominator"], 1)

    def test_norm5_H4_H8_H12_ratios(self) -> None:
        for first, second, multiplier in (
            (Gaussian(8, 1), Gaussian(7, 4), Gaussian(2, -1)),
            (Gaussian(9, 2), Gaussian(7, 6), Gaussian(2, 1)),
        ):
            payload = lineage_payload(first, second, multiplier)
            predictions = payload["harmonic_predictions"]
            self.assertEqual(
                (
                    predictions["H4"]["angular_ratio"]["numerator"],
                    predictions["H4"]["angular_ratio"]["denominator"],
                ),
                (-14, 25),
            )
            self.assertEqual(
                predictions["H4"]["target_expression"],
                "(-14/25)*5^(-13/8)",
            )
            self.assertEqual(
                (
                    predictions["H8"]["angular_ratio"]["numerator"],
                    predictions["H8"]["angular_ratio"]["denominator"],
                ),
                (-1054, 625),
            )
            self.assertEqual(
                (
                    predictions["H12"]["angular_ratio"]["numerator"],
                    predictions["H12"]["angular_ratio"]["denominator"],
                ),
                (23506, 15625),
            )

    def test_angular_ratios_compose_exactly_under_gaussian_multiplication(self) -> None:
        first = Gaussian(8, 1)
        second = Gaussian(7, 4)
        multiplier_1 = Gaussian(1, 1)
        multiplier_2 = Gaussian(2, -1)

        first_step = lineage_payload(first, second, multiplier_1)
        intermediate_first = first.multiply(multiplier_1)
        intermediate_second = second.multiply(multiplier_1)
        second_step = lineage_payload(
            intermediate_first, intermediate_second, multiplier_2
        )
        direct = lineage_payload(
            first, second, multiplier_1.multiply(multiplier_2)
        )

        self.assertEqual(
            direct["norm_ratio"],
            first_step["norm_ratio"] * second_step["norm_ratio"],
        )
        for harmonic in ("H4", "H8", "H12"):
            self.assertEqual(
                _ratio(direct["harmonic_predictions"][harmonic]),
                _ratio(first_step["harmonic_predictions"][harmonic])
                * _ratio(second_step["harmonic_predictions"][harmonic]),
            )

    def test_catalog_designs_preserve_translation_group_within_every_pair(self) -> None:
        catalog = default_catalog()
        self.assertEqual(catalog["schema_version"], 3)
        for section in (
            "doubling_lineages",
            "norm5_harmonic_discrimination",
            "N1105_edges",
        ):
            for edge in catalog[section].values():
                contract = edge["pair_translation_group_contract"]
                self.assertTrue(contract["parent_pair_matches"])
                self.assertTrue(contract["child_pair_matches"])
                self.assertEqual(contract["parent_smith_invariants"][0], 1)
                self.assertEqual(contract["child_smith_invariants"][0], 1)

    def test_N1105_catalog_has_four_orientations(self) -> None:
        catalog = default_catalog()
        orientations = set()
        for edge in catalog["N1105_edges"].values():
            orientations.add(tuple(edge["child"]["first_canonical"]["pair"]))
            orientations.add(tuple(edge["child"]["second_canonical"]["pair"]))
            self.assertEqual(edge["child"]["first_canonical"]["N"], 1105)
            self.assertEqual(edge["child"]["second_canonical"]["N"], 1105)
        self.assertEqual(
            orientations,
            {(33, 4), (32, 9), (31, 12), (24, 23)},
        )

    def test_N1105_commuting_diagram_has_three_genealogies_per_orientation(self) -> None:
        catalog = default_catalog()
        multiplicity = Counter()
        for edge in catalog["N1105_edges"].values():
            multiplicity[tuple(edge["child"]["first_canonical"]["pair"])] += 1
            multiplicity[tuple(edge["child"]["second_canonical"]["pair"])] += 1

        self.assertEqual(
            multiplicity,
            Counter(
                {
                    (33, 4): 3,
                    (32, 9): 3,
                    (31, 12): 3,
                    (24, 23): 3,
                }
            ),
        )


if __name__ == "__main__":
    unittest.main()
