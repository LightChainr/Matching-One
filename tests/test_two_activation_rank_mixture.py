from fractions import Fraction

from scripts.two_activation_rank_mixture import (
    beta_density,
    evaluate_power,
    exact_axis_L2_oracle,
    power_from_bernstein,
)


def test_binomial_tail_derivative_is_beta_density_at_exact_probe():
    N, k, p = 7, 3, Fraction(2, 5)
    tail_power = power_from_bernstein(
        [Fraction(int(n >= k)) for n in range(N + 1)]
    )
    derivative_power = [
        degree * tail_power[degree] for degree in range(1, N + 1)
    ]
    assert evaluate_power(derivative_power, p) == beta_density(N, k, p)


def test_existing_axis_L2_histogram_has_exact_two_activation_law():
    payload = exact_axis_L2_oracle()
    assert payload["joint_K1_K2_counts"] == [
        {"K1": 2, "K2": 3, "count": 16},
        {"K1": 3, "K2": 3, "count": 8},
    ]
    assert payload["microcanonical_q_by_n"] == ["-1", "-1", "-1/3", "1", "1"]
    assert payload["matching_power_coefficients_ascending"] == ["-1", "0", "4", "0", "-2"]
    assert payload["density_power_coefficients_ascending"] == ["0", "4", "0", "-4"]


def test_midpoint_gap_semantics_and_neutral_area_are_exact():
    payload = exact_axis_L2_oracle()
    assert payload["joint_midpoint_gap_counts"] == [
        {"C": "5/2", "G": 1, "count": 16},
        {"C": "3", "G": 0, "count": 8},
    ]
    semantics = payload["midpoint_gap_semantics"]
    assert semantics["mean_C"] == "8/3"
    assert semantics["mixture_mean_EC_over_N_plus_1"] == "8/15"
    assert semantics["mean_G"] == "2/3"
    assert semantics["integrated_rank1_probability_EG_over_N_plus_1"] == "2/15"


def test_root_is_two_CDF_balance():
    payload = exact_axis_L2_oracle()
    root = payload["root_two_CDF_balance"]
    assert abs(root["balance_residual"]) < 1e-14
    assert 0.54 < root["numerical_root"] < 0.55
