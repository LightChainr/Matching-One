"""Tests for the #595 v2 / #596 Gate-2 projective channel design.

Each test names the wrong number it would stop us believing.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p205_projective_channel_design import (  # noqa: E402
    BANDS,
    BASIS_CHANNELS,
    MAX_MULTIPLIER,
    Z_LEVEL,
    decide,
    fieller_set,
    folds,
    poles,
    sample_multiplier,
    window,
)
from score_p205_derivative_smith_loading import induced_a8_over_a4, rho_window  # noqa: E402


def estimator(a4, delta, var_a, var_d, cov):
    return {
        "A4": mp.mpf(a4), "delta": mp.mpf(delta),
        "var_A4": mp.mpf(var_a), "var_delta": mp.mpf(var_d),
        "cov_A4_delta": mp.mpf(cov),
    }


class FiellerGeometry(unittest.TestCase):
    def test_the_set_is_unbounded_exactly_when_A4_is_unresolved(self) -> None:
        """Stops us believing a bounded interval where the truth is unbounded.

        This is the whole correction to #595 v1: a delta-method error on
        ``delta/A4`` always returns a finite interval, including when ``A4`` is
        consistent with zero and no finite bound on the ratio exists.  If this
        returned an interval there, every reported sample cost would be a
        fiction.
        """

        mp.mp.dps = 40
        strong = estimator("1e-3", "0", "1e-9", "1e-9", "0")   # A4 at ~31 sigma
        weak = estimator("1e-5", "0", "1e-9", "1e-9", "0")     # A4 at ~0.3 sigma
        self.assertIsNotNone(fieller_set(strong, mp.mpf(1), True))
        self.assertIsNone(fieller_set(weak, mp.mpf(1), True))
        # And a large enough sample resolves the weak case.
        self.assertIsNotNone(fieller_set(weak, mp.mpf("1e5"), True))

    def test_the_set_shrinks_monotonically_with_sample(self) -> None:
        """Stops us believing a bisection run on a non-monotone target.

        ``sample_multiplier`` bisects on "does the set fit the window".  If more
        sample could widen the set the bisection would return an arbitrary point
        rather than the smallest sufficient multiplier.
        """

        mp.mp.dps = 40
        row = estimator("1e-3", "0", "4e-9", "9e-9", "1.5e-9")
        previous = None
        for multiplier in ("1", "4", "16", "64", "256"):
            interval = fieller_set(row, mp.mpf(multiplier), True)
            self.assertIsNotNone(interval)
            width = interval[1] - interval[0]
            if previous is not None:
                self.assertLess(width, previous)
            previous = width

    def test_the_cross_term_actually_changes_the_set(self) -> None:
        """Stops us believing a design score that ignores the A4/delta geometry.

        The objective #595 v1 got wrong was precisely that it never looked at the
        joint geometry.  If the cross-covariance had no effect here, this script
        would be v1 in different clothes.
        """

        mp.mp.dps = 40
        base = estimator("1e-3", "0", "4e-9", "9e-9", "0")
        tilted = estimator("1e-3", "0", "4e-9", "9e-9", "5e-9")
        a = fieller_set(base, mp.mpf(1), True)
        b = fieller_set(tilted, mp.mpf(1), True)
        self.assertNotAlmostEqual(float(a[0]), float(b[0]), places=6)


class WindowAndPoles(unittest.TestCase):
    def test_the_window_edges_map_exactly_onto_the_band(self) -> None:
        """Stops us believing a sample cost aimed at the wrong contamination level.

        The whole calculation is "shrink the Fieller set until it fits this
        window".  A window that did not correspond to the declared band would
        make every multiplier wrong by an unknown factor.
        """

        mp.mp.dps = 40
        for band in BANDS:
            low, high = window(band)
            for family in ("square", "rectangle"):
                family_low, family_high = rho_window(mp.mpf(band), family)
                self.assertGreaterEqual(low, family_low - mp.mpf("1e-30"))
                self.assertLessEqual(high, family_high + mp.mpf("1e-30"))
            for edge in (low, high):
                worst = max(
                    abs(induced_a8_over_a4(edge, family)["induced_A8_over_A4"])
                    for family in ("square", "rectangle")
                )
                self.assertLessEqual(worst, mp.mpf(band) + mp.mpf("1e-25"))

    def test_the_bands_lie_strictly_inside_the_poles(self) -> None:
        """Stops us believing pole exclusion is an extra requirement when it is implied.

        If a declared band reached past a pole, fitting the band would no longer
        guarantee the fitted A4 survives, and the verdict text saying pole
        exclusion is automatic would be wrong.
        """

        mp.mp.dps = 40
        low_pole, high_pole = poles()
        for band in BANDS:
            low, high = window(band)
            self.assertGreater(low, low_pole)
            self.assertLess(high, high_pole)


class SampleCost(unittest.TestCase):
    def test_cost_scales_as_the_inverse_of_the_offset_variance(self) -> None:
        """Stops us believing v1's figure of merit survived into v2.

        The Fieller half-width at large n is ``z sqrt(var(delta)/n)/|A4|``, so in
        a well-resolved regime doubling ``var(delta)`` must exactly double the
        required sample while doubling ``var(A4)`` must barely move it.  If the
        cost tracked ``var(A4)`` instead, this would be the old objective wearing
        a new name.
        """

        mp.mp.dps = 40
        # var(A4) deliberately negligible, so the asymptotic identity is exact.
        base = estimator("1e-3", "0", "1e-14", "9e-9", "0")
        more_delta = estimator("1e-3", "0", "1e-14", "1.8e-8", "0")
        more_a4 = estimator("1e-3", "0", "2e-14", "9e-9", "0")
        cost = sample_multiplier(base, "0.20")
        self.assertLess(abs(sample_multiplier(more_delta, "0.20") / cost - 2), 0.001)
        self.assertLess(abs(sample_multiplier(more_a4, "0.20") / cost - 1), 0.001)

    def test_the_exact_fieller_cost_exceeds_the_asymptotic_shortcut(self) -> None:
        """Stops us believing an optimistic cost taken from the leading-order formula.

        ``A4^2/var(delta)`` is the asymptotic figure of merit and it is what the
        docstring quotes, but it *understates* the requirement whenever ``A4`` is
        weakly resolved -- which is the entire regime this block sits in.  If the
        code had used the shortcut, every quoted sample cost would be too cheap
        in exactly the cases that matter.
        """

        mp.mp.dps = 40
        reach = min(abs(edge) for edge in window("0.20"))
        for var_a4 in ("1e-14", "4e-9"):
            row = estimator("1e-3", "0", var_a4, "9e-9", "0")
            asymptotic = (
                Z_LEVEL**2 * row["var_delta"] / (row["A4"] ** 2 * reach**2)
            )
            exact = sample_multiplier(row, "0.20")
            self.assertGreaterEqual(exact, asymptotic * mp.mpf("0.999"))
            if var_a4 == "4e-9":  # weakly resolved: the shortcut is visibly cheap
                self.assertGreater(exact / asymptotic, mp.mpf("1.02"))

    def test_a_hopeless_measured_point_is_reported_as_unreachable(self) -> None:
        """Stops us believing more sampling can rescue a real offset.

        At the measured point the Fieller set shrinks toward rho-hat, so if
        rho-hat is outside the band no multiplier works.  Returning a large
        finite number instead of ``inf`` would read as an expensive but buyable
        experiment.
        """

        mp.mp.dps = 40
        far = estimator("1e-3", "1.3e-3", "4e-9", "9e-9", "0")  # rho-hat = 1.3
        self.assertEqual(sample_multiplier(far, "0.20", design_null=False), mp.mpf("inf"))
        self.assertLess(sample_multiplier(far, "0.20", design_null=True), MAX_MULTIPLIER)

    def test_a_tighter_band_never_costs_less(self) -> None:
        """Stops us believing a band ordering that got inverted somewhere.

        0.05 is strictly inside 0.20, so its multiplier must be at least as
        large.  An inversion would mean the window construction or the band
        strings were mismatched, and the cheap-looking column would be the
        stricter one.
        """

        mp.mp.dps = 40
        row = estimator("1e-3", "0", "4e-9", "9e-9", "1e-9")
        self.assertGreaterEqual(
            sample_multiplier(row, "0.05"), sample_multiplier(row, "0.20")
        )


class Folds(unittest.TestCase):
    def test_the_folds_partition_the_batches_without_overlap(self) -> None:
        """Stops us believing a held-out fold that the discovery folds already saw.

        Any overlap would leak the selection into its own validation and turn the
        transport check into a restatement of the discovery result.
        """

        for count in (99, 100, 101):
            partition = folds(count)
            flat = [b for fold in partition for b in fold]
            self.assertEqual(sorted(flat), list(range(count)))
            self.assertEqual(len(flat), len(set(flat)))
            self.assertEqual(len(partition), 3)


class Verdict(unittest.TestCase):
    def _report(self, raw_beats, stable_beats, raw_cost, stable_cost):
        def rounds(beats, cost):
            return {
                "rounds": [{
                    "validation_multiplier": mp.mpf(cost),
                    "validation_saving_vs_baseline": mp.mpf(3),
                    "beats_baseline_held_out": beats,
                    "chosen_readout": "Sp", "chosen_p": "x",
                }],
                "all_rounds_beat_baseline": beats,
                "selection_is_stable": True,
                "held_out_multiplier_spread": 1.0,
            }
        return {
            "cross_fit": {
                band: {"raw": rounds(raw_beats, raw_cost),
                       "stabilized": rounds(stable_beats, stable_cost)}
                for band in BANDS
            }
        }

    def test_disagreeing_variants_produce_no_branch(self) -> None:
        """Stops us believing the friendlier of two variants that contradict.

        The stabilized cross-fit was introduced after seeing the raw one behave
        badly.  If the code silently reported whichever variant gave a branch,
        that post-hoc choice would decide a production purchase.
        """

        got = decide([self._report(False, True, 10, 10)])
        self.assertIn("NOT_ESTABLISHED", got["verdict"])
        self.assertEqual(got["raw_branch"], "B")
        self.assertEqual(got["stabilized_branch"], "A")

    def test_agreeing_variants_give_the_expected_branch(self) -> None:
        """Stops us believing a Gate-2 branch letter assigned by the wrong rule.

        Branch C -- "even the best honest channel is unaffordable" -- is the one
        that stops N=650 being bought as a pure A8 experiment, so it must not be
        reachable from a cheap cost or unreachable from an expensive one.
        """

        self.assertIn("GATE2_A", decide([self._report(True, True, 10, 10)])["verdict"])
        self.assertIn("GATE2_B", decide([self._report(False, False, 10, 10)])["verdict"])
        self.assertIn("GATE2_C", decide([self._report(True, True, 5000, 5000)])["verdict"])


class DeclaredManifest(unittest.TestCase):
    def test_the_declared_constants_are_what_the_write_up_will_quote(self) -> None:
        """Stops us believing a cost computed at a confidence level nobody declared.

        The multiplier scales as z^2, so a quietly changed level rescales every
        number in the artifact.  The bands and the excluded dependent channel are
        pinned for the same reason.
        """

        self.assertEqual(BANDS, ("0.20", "0.05"))
        self.assertNotIn("D", BASIS_CHANNELS)
        self.assertLess(abs(Z_LEVEL - mp.mpf("1.96")), mp.mpf("0.001"))
        low, high = poles()
        self.assertLess(abs(abs(low) - mp.mpf("1.2958379")), mp.mpf("1e-6"))
        self.assertLess(abs(high + low), mp.mpf("1e-30"))


if __name__ == "__main__":
    unittest.main()
