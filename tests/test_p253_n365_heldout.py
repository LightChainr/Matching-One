import math

from scripts.score_p253_n365_heldout import (
    covariance_of_mean,
    dyadic_residual,
    fractional_prediction,
)


def test_dyadic_residual_vanishes_for_shared_recurrence():
    # Both old rows and the held-out row obey T=1.1, D=0.18.
    values = [1.0, 0.7, 0.59, -0.4, 0.2, 0.292, 0.3, -0.1, -0.164]
    assert math.isclose(dyadic_residual(values)[0], 0.0, abs_tol=1e-12)


def test_fractional_real_two_mode_matches_direct_value():
    first, second = 0.8, 0.3
    y0 = 1.2 - 0.4
    y1 = 1.2 * first - 0.4 * second
    expected = 1.2 * first ** math.log(3.5, 2) - 0.4 * second ** math.log(3.5, 2)
    assert math.isclose(
        fractional_prediction(first + second, first * second, y0, y1, "R2"), expected
    )


def test_fractional_complex_pair_matches_direct_value():
    modulus, theta = 0.75, 0.6
    y0 = 1.1
    y1 = modulus * (1.1 * math.cos(theta) - 0.2 * math.sin(theta))
    step = math.log(3.5, 2)
    expected = modulus**step * (1.1 * math.cos(theta * step) - 0.2 * math.sin(theta * step))
    assert math.isclose(
        fractional_prediction(2 * modulus * math.cos(theta), modulus**2, y0, y1, "C2"),
        expected,
    )


def test_covariance_of_mean_keeps_scale():
    covariance = covariance_of_mean([[0.0, 1.0], [2.0, 3.0]])
    assert covariance == [[1.0, 1.0], [1.0, 1.0]]
