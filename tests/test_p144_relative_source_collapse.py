from fractions import Fraction

from scripts.p144_relative_source_collapse import axis_L3_relative_collapse


def test_axis_L3_terminal_support_is_exactly_three_alexander_monomials():
    payload = axis_L3_relative_collapse()
    factorization = payload["exact_factorization"]
    assert factorization["terminal_rank_pairs"] == [[0, 2], [1, 1], [2, 0]]
    assert factorization["terminal_monomials"] == ["y^2", "x*y", "x^2"]
    assert payload["checks"]["Phi_equals_xy_Zrel_symbolically_per_monomial"]


def test_exact_sector_probabilities_normalize_and_match_relative_derivative():
    payload = axis_L3_relative_collapse(Fraction(2, 5))
    sectors = payload["sector_probabilities"]
    total = sum(
        Fraction(row["numerator"], row["denominator"]) for row in sectors.values()
    )
    assert total == 1
    assert payload["relative_source"]["mean_q_equals_M"] == payload["relative_source"]["matching_from_bivariate_derivative"]


def test_rank_covariance_is_strictly_rank_one_in_open_interval():
    for p in (Fraction(1, 5), Fraction(2, 5), Fraction(1, 2), Fraction(4, 5)):
        covariance = axis_L3_relative_collapse(p)["rank_covariance"]
        assert covariance["rank"] == 1
        assert covariance["determinant"]["numerator"] == 0
        assert covariance["strict_rank_one_for_0_lt_p_lt_1"]


def test_diagonal_source_is_deterministic_and_q_cumulant_closes():
    payload = axis_L3_relative_collapse()
    assert payload["diagonal_source"]["rank_sum"] == 2
    assert payload["diagonal_source"]["stochastic_diagonal_cumulants_order_ge_2"] == 0
    closure = payload["three_state_closure"]
    assert closure["third_cumulant_direct"] == closure["third_cumulant_from_closure"]
    assert payload["checks"]["third_cumulant_closure_exact"]
