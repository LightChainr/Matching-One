import math

from scripts.physical_root_uniqueness import (
    activation_bracket,
    activation_density,
    beta_density,
    exact_axis_L2_oracle,
)


def test_every_interior_beta_component_is_strictly_positive():
    for N in range(1, 9):
        for k in range(1, N + 1):
            for p in (0.01, 0.2, 0.5, 0.83, 0.99):
                assert beta_density(N, k, p) > 0.0


def test_endpoint_atoms_are_exactly_the_zero_derivative_degeneracies():
    N = 5
    assert beta_density(N, 0, 0.4) == 0.0
    assert beta_density(N, N + 1, 0.4) == 0.0
    assert activation_density(N, {0: 3, N + 1: 2}, 0.4) == 0.0
    assert activation_density(N, {0: 3, 2: 1, N + 1: 2}, 0.4) > 0.0


def test_stochastically_ordered_activation_histograms_bracket_root():
    bracket = activation_bracket(4, {2: 16, 3: 8}, {3: 24})
    assert bracket["onset_mixture_median"] < bracket["physical_matching_root"]
    assert bracket["physical_matching_root"] < bracket["completion_mixture_median"]
    assert math.isclose(
        bracket["physical_matching_root"], math.sqrt(1.0 - 1.0 / math.sqrt(2.0))
    )
    assert bracket["M_derivative_at_root"] > 0.0


def test_equal_activations_collapse_the_bracket():
    bracket = activation_bracket(6, {3: 10}, {3: 10})
    assert math.isclose(bracket["onset_mixture_median"], bracket["physical_matching_root"])
    assert math.isclose(bracket["physical_matching_root"], bracket["completion_mixture_median"])


def test_committed_tiny_oracle_reports_unique_simple_root():
    payload = exact_axis_L2_oracle()
    bracket = payload["axis_L2_exact_histogram"]["bracket"]
    assert bracket["M_derivative_at_root"] > 0.0
    assert abs(bracket["root_balance_residual"]) < 1e-14
    assert payload["axis_L2_exact_histogram"]["strict_positive_derivative_grid_check"]
