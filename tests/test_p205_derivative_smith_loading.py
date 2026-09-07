"""Tests for the #589 derivative-channel Smith-loading rescore.

Each test names the wrong number it would stop us believing.
"""

from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p205_derivative_smith_loading import (  # noqa: E402
    BANDED_BAND,
    CHANNELS,
    CLEAN_BAND,
    CONTROL_CHANNEL,
    EXPECTED_LEAKAGE,
    GEOMETRY_GENERATORS,
    GEOMETRY_ORDER,
    RESCALING_PAIRS,
    SIGMA_MULTIPLE,
    VERDICT_CHANNELS,
    angular_functional,
    check_h4_null_is_the_cos4_interpolant,
    cos_harmonic,
    h4_offset_functional,
    induced_a8_over_a4,
    leakage,
    quadratic_form,
    rescaling_control,
    rho_window,
    sample_multiple,
    smith_loading,
    verify_leakage,
)

H4_WEIGHTS = {325: (5, -11, 6), 425: (20, 13, -33)}


class AngularCoordinate(unittest.TestCase):
    def test_cos_harmonic_matches_the_transcendental_angle(self) -> None:
        """Stops us believing a cos(4 theta) that is not the lattice angle.

        Every weight in this script -- the frozen null, the A4 functional, the
        N=650 leakage -- is built from ``cos_harmonic``.  If it disagreed with
        ``cos(4 atan2(b, a))`` the exact rational arithmetic would be exactly
        the wrong quantity, and nothing downstream would notice.
        """

        mp.mp.dps = 50
        for a, b in ((15, 10), (17, 6), (18, 1), (25, 5), (23, 11), (19, 17)):
            theta = mp.atan2(mp.mpf(b), mp.mpf(a))
            for spin in (4, 8):
                exact = cos_harmonic(a, b, spin)
                value = mp.mpf(exact.numerator) / exact.denominator
                self.assertLess(abs(value - mp.cos(spin * theta)), mp.mpf("1e-40"))

    def test_frozen_H4_null_is_exactly_the_two_cyclic_row_interpolant(self) -> None:
        """Stops us believing ``residual`` is a pure noncyclic offset when it is not.

        The whole reading of ``delta`` depends on the frozen #205 integer null
        being the exact ``cos4`` interpolation through the two cyclic rows.  If
        it were merely close, ``delta`` would mix an offset with an angular
        model error and ``rho`` would not be the quantity #583 needs.
        """

        for n, weights in H4_WEIGHTS.items():
            self.assertEqual(check_h4_null_is_the_cos4_interpolant(n, weights), 0)


class Functionals(unittest.TestCase):
    def _planted(self, n: int, constant: mp.mpf, amplitude: mp.mpf, offset: mp.mpf):
        generators = GEOMETRY_GENERATORS[n]
        values = []
        for index, (a, b) in enumerate(generators):
            c4 = cos_harmonic(a, b, 4)
            value = constant + amplitude * (mp.mpf(c4.numerator) / c4.denominator)
            if index == 0:  # the noncyclic C row
                value += offset
            values.append(value)
        return values

    def test_offset_functional_reads_the_planted_offset_and_nothing_else(self) -> None:
        """Stops us believing a delta that is contaminated by C or A4.

        A wrong normalization here would rescale every reported ``delta`` by a
        constant, which is invisible in a z-score but changes ``rho`` -- and
        ``rho`` is the only number this script exports to #583.
        """

        mp.mp.dps = 40
        for n in H4_WEIGHTS:
            weights = h4_offset_functional(n, H4_WEIGHTS[n])
            for offset in (mp.mpf(0), mp.mpf("0.25"), mp.mpf("-3")):
                point = self._planted(n, mp.mpf("7"), mp.mpf("-2"), offset)
                got = mp.fsum(weights[i] * point[i] for i in range(3))
                self.assertLess(abs(got - offset), mp.mpf("1e-30"))

    def test_angular_functional_reads_A4_and_ignores_the_noncyclic_row(self) -> None:
        """Stops us believing an A4 denominator that the offset has leaked into.

        If the C row had nonzero weight, the same quotient response would appear
        in the numerator and the denominator of ``rho`` and partially cancel --
        which would understate the contamination this script exists to bound.
        """

        mp.mp.dps = 40
        for n in H4_WEIGHTS:
            weights = angular_functional(n)
            self.assertEqual(weights[0], 0)
            for amplitude in (mp.mpf("1.5"), mp.mpf("-0.125")):
                point = self._planted(n, mp.mpf("7"), amplitude, mp.mpf("11"))
                got = mp.fsum(weights[i] * point[i] for i in range(3))
                self.assertLess(abs(got - amplitude), mp.mpf("1e-28"))


class Leakage(unittest.TestCase):
    def test_recomputed_leakage_reproduces_issue_591(self) -> None:
        """Stops us believing an N=650 leakage that is not #591's.

        The bridge from this measurement to #583 is a multiplication by these
        three rationals.  A wrong A8 leakage would rescale the induced ratio
        directly; a wrong A4 leakage would put the fitted-A4 pole in the wrong
        place and change which confidence intervals are reported as unusable.
        """

        report = verify_leakage()
        self.assertTrue(report["matches_issue_591"])
        square = leakage("square")
        rectangle = leakage("rectangle")
        self.assertEqual(square["A8"], EXPECTED_LEAKAGE["A8"])
        self.assertEqual(square["A8"], rectangle["A8"])
        self.assertEqual(square["A4"], -rectangle["A4"])
        self.assertEqual(square["A4"], EXPECTED_LEAKAGE["A4"])

    def test_leakage_is_the_saturated_fit_of_a_unit_noncyclic_offset(self) -> None:
        """Stops us believing a leakage column read out of the wrong matrix slot.

        Checked by independent means: instead of inverting, solve the 3x3
        harmonic system directly for the response vector that is one on the
        noncyclic row and zero elsewhere, and confirm the fitted coefficients
        are the reported leakage.
        """

        from score_p205_derivative_smith_loading import (
            N650_NONCYCLIC_INDEX,
            N650_RECTANGLE,
            N650_SQUARE,
            invert3,
        )

        for family, rows in (("square", N650_SQUARE), ("rectangle", N650_RECTANGLE)):
            design = [
                [Fraction(1), cos_harmonic(a, b, 4), cos_harmonic(a, b, 8)]
                for a, b in rows
            ]
            response = [Fraction(0)] * 3
            response[N650_NONCYCLIC_INDEX[family]] = Fraction(1)
            inverse = invert3(design)
            fitted = [
                sum(inverse[i][j] * response[j] for j in range(3)) for i in range(3)
            ]
            # The fit must reproduce the response exactly (saturated design).
            for i in range(3):
                self.assertEqual(
                    sum(design[i][k] * fitted[k] for k in range(3)), response[i]
                )
            got = leakage(family)
            self.assertEqual([got["C"], got["A4"], got["A8"]], fitted)


class DeltaMethod(unittest.TestCase):
    def _scored(self, n: int, point, covariance):
        return {
            "N": n,
            "channel": "Sp",
            "point": dict(zip(GEOMETRY_ORDER, point)),
            "covariance": covariance,
        }

    def test_ratio_variance_matches_the_gradient_form(self) -> None:
        """Stops us believing a rho error that drops the delta/A4 cross term.

        The code expands ``var(delta)/A4^2 - 2 delta cov/A4^3 + delta^2
        var(A4)/A4^4``.  Independent means: build the gradient of ``delta/A4``
        in the three-point space and evaluate ``g^T Sigma g`` directly.  A
        dropped or sign-flipped cross term would silently shrink or inflate
        every reported bound, and the verdict is read off those bounds.
        """

        mp.mp.dps = 40
        covariance = [
            [mp.mpf("4e-8"), mp.mpf("1.1e-8"), mp.mpf("-0.7e-8")],
            [mp.mpf("1.1e-8"), mp.mpf("3e-8"), mp.mpf("0.9e-8")],
            [mp.mpf("-0.7e-8"), mp.mpf("0.9e-8"), mp.mpf("5e-8")],
        ]
        for n in H4_WEIGHTS:
            point = [mp.mpf("0.0031"), mp.mpf("0.0019"), mp.mpf("0.0044")]
            loading = smith_loading(self._scored(n, point, covariance), H4_WEIGHTS[n])
            weights = h4_offset_functional(n, H4_WEIGHTS[n])
            angular = angular_functional(n)
            delta, amplitude = loading["delta"], loading["A4"]
            gradient = [
                (weights[i] * amplitude - delta * angular[i]) / amplitude**2
                for i in range(3)
            ]
            expected = mp.sqrt(quadratic_form(gradient, covariance, gradient))
            self.assertLess(
                abs(loading["rho_standard_error"] - expected),
                mp.mpf("1e-25") * expected,
            )

    def test_rho_is_invariant_under_rescaling_the_channel(self) -> None:
        """Stops us believing a rho that still carries the channel's units.

        ``rho`` is exported to #583 across a size change and a channel change,
        so it has to be dimensionless.  Scaling the point and its covariance
        consistently must leave it untouched -- this is the property that makes
        the ``D == M/2`` control in the artifact meaningful rather than lucky.
        """

        mp.mp.dps = 40
        covariance = [[mp.mpf("2e-8") if i == j else mp.mpf("0.3e-8") for j in range(3)] for i in range(3)]
        point = [mp.mpf("0.0031"), mp.mpf("0.0019"), mp.mpf("0.0044")]
        base = smith_loading(self._scored(325, point, covariance), H4_WEIGHTS[325])
        factor = mp.mpf("17.5")
        scaled = smith_loading(
            self._scored(
                325,
                [factor * value for value in point],
                [[factor**2 * value for value in row] for row in covariance],
            ),
            H4_WEIGHTS[325],
        )
        self.assertLess(abs(base["rho"] - scaled["rho"]), mp.mpf("1e-28"))
        self.assertLess(
            abs(base["rho_standard_error"] - scaled["rho_standard_error"]),
            mp.mpf("1e-28"),
        )


class InducedRatio(unittest.TestCase):
    def test_window_edges_hit_the_band_exactly(self) -> None:
        """Stops us believing a sample-size plan aimed at the wrong rho target.

        ``rho_window`` inverts the induced-ratio map; if the inversion were
        approximate the reported ``x samples`` would be off by the square of the
        error, and the number quoted to #583 as the cost of a clean answer would
        be wrong in the direction that makes the experiment look cheap.
        """

        mp.mp.dps = 40
        for band in (CLEAN_BAND, BANDED_BAND, mp.mpf("0.4")):
            for family in ("square", "rectangle"):
                low, high = rho_window(band, family)
                for edge in (low, high):
                    got = induced_a8_over_a4(edge, family)["induced_A8_over_A4"]
                    self.assertLess(abs(abs(got) - band), mp.mpf("1e-30"))

    def test_the_pole_annihilates_the_fitted_A4(self) -> None:
        """Stops us believing a finite contamination ratio where none exists.

        At ``rho = -1/leak_A4`` a quotient response cancels the true angular
        amplitude exactly, so N=650's fitted A4 is pure nuisance and no
        ``A8/A4`` is defined.  If the pole were misplaced, intervals that
        contain this catastrophic case would be reported as merely large.
        """

        mp.mp.dps = 40
        for family in ("square", "rectangle"):
            info = induced_a8_over_a4(mp.mpf(0), family)
            pole = info["pole"]
            self.assertLess(abs(1 + info["leak_A4"] * pole), mp.mpf("1e-30"))
            near = induced_a8_over_a4(pole * (1 - mp.mpf("1e-6")), family)
            self.assertGreater(abs(near["induced_A8_over_A4"]), mp.mpf("1e5"))

    def test_sample_multiple_is_the_square_of_the_error_ratio(self) -> None:
        """Stops us believing a cost estimate that scales like 1/n instead of 1/sqrt(n).

        Getting this exponent wrong would understate a 26x sample requirement as
        5x, which is exactly the difference between "rescore an archive" and
        "commission a new run".
        """

        mp.mp.dps = 40
        error = mp.mpf("0.35")
        self.assertLess(
            abs(sample_multiple(error / 2, BANDED_BAND) - sample_multiple(error, BANDED_BAND) / 4),
            mp.mpf("1e-25"),
        )
        reach = min(
            min(abs(edge) for edge in rho_window(BANDED_BAND, family))
            for family in ("square", "rectangle")
        )
        self.assertLess(
            abs(sample_multiple(reach / SIGMA_MULTIPLE, BANDED_BAND) - 1), mp.mpf("1e-25")
        )
        self.assertEqual(sample_multiple(mp.mpf("nan"), BANDED_BAND), mp.mpf("inf"))


class RescalingControl(unittest.TestCase):
    def test_control_fails_when_two_rescaled_channels_disagree(self) -> None:
        """Stops us believing a passing control that cannot fail.

        The artifact reports ``rho(D) == rho(M)`` as a PASS.  A control that
        agrees no matter what is decoration; this plants a disagreement and
        requires it to be caught.
        """

        blocks = {
            "M": {"325": {"rho": "1.5"}, "425": {"rho": "-1.2"}},
            "D": {"325": {"rho": "1.5"}, "425": {"rho": "-1.2"}},
        }
        self.assertTrue(rescaling_control(blocks)["all_agree"])
        blocks["D"]["425"] = {"rho": "-1.2000000001"}
        self.assertFalse(rescaling_control(blocks)["all_agree"])


class DeclaredManifest(unittest.TestCase):
    def test_channels_exist_and_the_verdict_excludes_the_published_one(self) -> None:
        """Stops us believing a verdict scored on a number that was already seen.

        #205 published the ``M`` channel.  Folding it into the decision would
        break GOVERNANCE 2C -- and because ``D`` is exactly ``M/2``, letting
        ``D`` into the verdict would smuggle the same seen number back in.
        """

        from analyze_matching_parity_derivatives_fast import obs, H

        histogram = H(
            n=4, a=1, b=0, orientation="first", batch=0, samples=2,
            minus=[0, 1, 1, 0, 0], plus=[1, 1, 0, 0, 0],
        )
        available = set(obs(histogram, mp.mpf("0.5")))
        self.assertTrue(set(CHANNELS) <= available)
        self.assertIn(CONTROL_CHANNEL, CHANNELS)
        self.assertNotIn(CONTROL_CHANNEL, VERDICT_CHANNELS)
        for left, right in RESCALING_PAIRS:
            self.assertFalse({left, right} & set(VERDICT_CHANNELS))
        self.assertLess(CLEAN_BAND, BANDED_BAND)


if __name__ == "__main__":
    unittest.main()
