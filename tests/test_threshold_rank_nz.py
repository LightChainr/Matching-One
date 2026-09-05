
from __future__ import annotations
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_threshold_ranks import (  # noqa: E402
    matching_derivative,
    matching_root,
    matching_value,
    read_counts,
)
from exact_matching_polynomial import (  # noqa: E402
    bernstein_to_power,
    evaluate_power,
)
from integer_period_torus import (  # noqa: E402
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
)
from threshold_rank_nz import (  # noqa: E402
    counter_permutation,
    enumerate_exact,
    simulate,
    threshold_ranks,
    write_counts,
)


def direct_cross_bernstein_counts(geometry):
    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        black, _ = classify_configuration(geometry, active)
        white, _ = classify_configuration(
            geometry, [not value for value in active], matching=True
        )
        counts[sum(active)] += int(black.cross) - int(white.cross)
    return counts


class ThresholdConventionTests(unittest.TestCase):
    def test_reverse_off_by_one_matches_direct_scan(self) -> None:
        geometry = axis_integer_torus(3)
        for counter in range(50):
            permutation = counter_permutation(geometry.n, 17, counter)
            k_minus, k_plus = threshold_ranks(geometry, permutation)
            for k in range(geometry.n + 1):
                active = [False] * geometry.n
                for vertex in permutation[:k]:
                    active[vertex] = True
                black, _ = classify_configuration(geometry, active)
                white, _ = classify_configuration(
                    geometry, [not value for value in active], matching=True
                )
                self.assertEqual(black.cross, k >= k_plus)
                self.assertEqual(white.cross, k < k_minus)

    def test_kminus_not_above_kplus_for_tiny_exact_permutations(self) -> None:
        for geometry in (
            axis_integer_torus(2),
            diamond_integer_torus(2),
            gaussian_integer_torus(2, 1),
        ):
            counts = enumerate_exact(geometry)
            self.assertTrue(all(k_minus <= k_plus for k_minus, k_plus in counts.joint))
            self.assertEqual(counts.sample_count, math.factorial(geometry.n))


class ExactReconstructionTests(unittest.TestCase):
    def test_histograms_reconstruct_direct_M_and_derivative(self) -> None:
        mp.mp.dps = 70
        for geometry in (axis_integer_torus(2), gaussian_integer_torus(2, 1)):
            rank_counts = enumerate_exact(geometry)
            power = bernstein_to_power(direct_cross_bernstein_counts(geometry))
            derivative = [degree * value for degree, value in enumerate(power)][1:]
            for p_text in ("0.17", "0.5", "0.731"):
                p = mp.mpf(p_text)
                reconstructed = matching_value(
                    geometry.n,
                    rank_counts.sample_count,
                    rank_counts.kminus,
                    rank_counts.kplus,
                    p,
                )
                reconstructed_prime = matching_derivative(
                    geometry.n,
                    rank_counts.sample_count,
                    rank_counts.kminus,
                    rank_counts.kplus,
                    p,
                )
                self.assertLess(abs(reconstructed - evaluate_power(power, p)), mp.mpf("1e-60"))
                self.assertLess(
                    abs(reconstructed_prime - evaluate_power(derivative, p)),
                    mp.mpf("1e-60"),
                )

            root = matching_root(
                geometry.n,
                rank_counts.sample_count,
                rank_counts.kminus,
                rank_counts.kplus,
            )
            self.assertLess(abs(evaluate_power(power, root)), mp.mpf("1e-60"))


class PersistenceAndCLITests(unittest.TestCase):
    def test_counts_round_trip_and_counter_chunking(self) -> None:
        geometry = axis_integer_torus(3)
        complete = simulate(geometry, 40, seed=99, counter_start=10)
        first = simulate(geometry, 15, seed=99, counter_start=10)
        second = simulate(geometry, 25, seed=99, counter_start=25)
        self.assertEqual(
            complete.joint,
            {
                key: first.joint.get(key, 0) + second.joint.get(key, 0)
                for key in set(first.joint) | set(second.joint)
            },
        )

        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_counts(
                output,
                geometry,
                complete,
                seed=99,
                counter_start=10,
                elapsed_seconds=1.25,
            )
            n, samples, minus, plus, metadata = read_counts(output)
            self.assertEqual((n, samples), (geometry.n, 40))
            self.assertEqual(minus, complete.kminus)
            self.assertEqual(plus, complete.kplus)
            self.assertEqual(
                metadata["first_second_joint_integer_moments"], complete.moments()
            )
            self.assertEqual(
                metadata["rng"]["sample_counter_stop_exclusive"], 50
            )
            self.assertTrue((output / "checksums.sha256").is_file())

    def test_analysis_cli_writes_recomputable_summary(self) -> None:
        geometry = axis_integer_torus(2)
        counts = enumerate_exact(geometry)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            write_counts(
                output,
                geometry,
                counts,
                seed=0,
                counter_start=0,
                elapsed_seconds=0.0,
            )
            summary = output / "derived_summary.json"
            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "analyze_threshold_ranks.py"),
                    "--input-dir", str(output),
                    "--p", "0.317",
                    "--json", str(summary),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            payload = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(payload["sample_count"], math.factorial(geometry.n))
            root = mp.mpf(payload["root"])
            power = bernstein_to_power(direct_cross_bernstein_counts(geometry))
            self.assertLess(abs(evaluate_power(power, root)), mp.mpf("1e-45"))


if __name__ == "__main__":
    unittest.main()
