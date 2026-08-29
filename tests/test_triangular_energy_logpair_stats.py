from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from triangular_energy_logpair_stats import (  # noqa: E402
    PAIR_ORDER,
    PRODUCT_ORDER,
    SufficientSums,
    black_cluster_roots,
    brute_force_sign_average,
    configuration_values,
    covariance_of_mean,
    exact_oracle,
    sign_integrated_four_spin,
    tiny_oracle_pairs,
    translation_averaged_configuration_sums,
    triangular_edges,
    unbiased_moment_vector,
)


class ClusterSignIntegrationTests(unittest.TestCase):
    def test_parity_rule_matches_explicit_touched_cluster_signs(self) -> None:
        width, height = 6, 3
        edges = triangular_edges(width, height)
        pairs = tiny_oracle_pairs(width, height)
        mask_limit = (1 << (width * height)) - 1
        for index in range(64):
            mask = (index * 104729 + 8191) & mask_limit
            roots = black_cluster_roots(width, height, mask, edges)
            for first, second in (
                ("L1", "L2"),
                ("L1", "D2"),
                ("D1", "L2"),
                ("D1", "D2"),
            ):
                analytic = sign_integrated_four_spin(
                    roots, pairs[first], pairs[second]
                )
                explicit = brute_force_sign_average(
                    roots, pairs[first], pairs[second]
                )
                self.assertEqual(Fraction(analytic), explicit)

    def test_cross_configuration_u_statistic_is_unbiased(self) -> None:
        width, height = 6, 3
        edges = triangular_edges(width, height)
        pairs = tiny_oracle_pairs(width, height)
        mask_limit = (1 << (width * height)) - 1
        observations = []
        for index in range(16):
            mask = (index * 104729 + 8191) & mask_limit
            roots = black_cluster_roots(width, height, mask, edges)
            observations.append(configuration_values(roots, pairs))

        count = len(observations)
        pair_means = {
            name: Fraction(sum(row[0][name] for row in observations), count)
            for name in PAIR_ORDER
        }
        four_means = {
            name: Fraction(sum(row[1][name] for row in observations), count)
            for name in PRODUCT_ORDER
        }
        target = (
            four_means["LL"] - pair_means["L1"] * pair_means["L2"],
            (
                four_means["L1D2"] - pair_means["L1"] * pair_means["D2"]
                + four_means["D1L2"] - pair_means["D1"] * pair_means["L2"]
            )
            / 2,
            four_means["DD"] - pair_means["D1"] * pair_means["D2"],
        )
        average = [Fraction(0), Fraction(0), Fraction(0)]
        for first in observations:
            for second in observations:
                sums = SufficientSums.empty()
                sums.add(*first)
                sums.add(*second)
                value = unbiased_moment_vector(sums)
                for coordinate in range(3):
                    average[coordinate] += value[coordinate]
        average = [value / (count * count) for value in average]
        self.assertEqual(tuple(average), target)

    def test_translation_sums_use_declared_placement_normalization(self) -> None:
        width = height = 8
        mask = int("1011010010110100101101001011010010110100101101001011010010110100", 2)
        roots = black_cluster_roots(width, height, mask)
        pair_sums, four_sums = translation_averaged_configuration_sums(
            roots, width, height, delta_radius=2, center_displacement=(4, 4)
        )
        sums = SufficientSums.empty()
        sums.add(pair_sums, four_sums)
        sums.add(pair_sums, four_sums)
        observed = unbiased_moment_vector(sums, placements=width * height)
        placement_count = width * height
        expected = (
            Fraction(four_sums["LL"], placement_count)
            - Fraction(pair_sums["L1"] * pair_sums["L2"], placement_count**2),
            (
                Fraction(four_sums["L1D2"], placement_count)
                - Fraction(pair_sums["L1"] * pair_sums["D2"], placement_count**2)
                + Fraction(four_sums["D1L2"], placement_count)
                - Fraction(pair_sums["D1"] * pair_sums["L2"], placement_count**2)
            )
            / 2,
            Fraction(four_sums["DD"], placement_count)
            - Fraction(pair_sums["D1"] * pair_sums["D2"], placement_count**2),
        )
        self.assertEqual(observed, expected)

    def test_three_coordinate_block_covariance_retains_off_diagonals(self) -> None:
        covariance = covariance_of_mean(
            ((1.0, 2.0, 4.0), (2.0, 4.0, 8.0), (4.0, 8.0, 16.0))
        )
        self.assertEqual(len(covariance), 3)
        self.assertGreater(covariance[0][1], 0.0)
        self.assertEqual(covariance[1][0], covariance[0][1])


class EnergyLogpairProtocolTests(unittest.TestCase):
    def test_committed_tiny_oracle_is_exactly_reproducible(self) -> None:
        path = ROOT / "predictions/p234_triangular_energy_logpair_tiny_oracle_20260829.json"
        self.assertEqual(json.loads(path.read_text(encoding="utf-8")), exact_oracle())

    def test_oracle_rationals_and_boundary_are_frozen(self) -> None:
        payload = exact_oracle()
        self.assertEqual(payload["geometry"]["configurations"], 262144)
        self.assertEqual(payload["explicit_cluster_sign_checks"], 1048576)
        self.assertEqual(
            [row["fraction"] for row in payload["exact_centered_moments"]],
            [
                "584500095/68719476736",
                "714261375/68719476736",
                "584500095/68719476736",
            ],
        )
        self.assertIn("not", payload["boundary"])

    def test_protocol_freezes_ordered_limits_but_not_kappa(self) -> None:
        protocol = yaml.safe_load(
            (ROOT / "predictions/p234_triangular_energy_logpair_protocol_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(protocol["status"], "protocol_freeze_no_universal_coefficient")
        self.assertEqual(
            protocol["unbiased_block_estimator"]["output_order"],
            ["E_a_E_a", "E_a_E_a_delta_symmetric", "E_a_delta_E_a_delta"],
        )
        self.assertEqual(
            protocol["ordered_limits"]["forbidden"],
            "simultaneous_free_fit_in_a_and_delta",
        )
        self.assertFalse(protocol["normalization_boundary"]["universal_coefficient_frozen"])


if __name__ == "__main__":
    unittest.main()
