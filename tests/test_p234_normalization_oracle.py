import math

from scripts.p234_normalization_oracle import (
    gauge_exponents,
    normalization_dictionary,
    render,
)


def test_paper_C1_enters_endpoint_normalization_through_its_square_root():
    row = normalization_dictionary(C1=9.0, C2=5.0, CL=7.0, G2=11.0, F3=13.0)
    assert row["spin_pair_amplitude_sqrt_C1"] == 3.0
    assert math.isclose(row["LD_continuum"], 55.0 / 3.0)
    assert row["DD_log_2delta_slope"] == -154.0
    assert math.isclose(row["kappa_proxy_sqrt_C1_CL_over_C2"], 21.0 / 5.0)
    assert math.isclose(row["paper_kappa_C1_CL_over_C2"], 63.0 / 5.0)


def test_extra_top_spin_spin_slope_recovers_CL_without_bottom_field_gauge():
    row = normalization_dictionary(C1=4.0, C2=3.0, CL=8.0, G2=5.0, F3=7.0)
    assert row["top_spin_spin_log_2delta_slope"] == -56.0
    assert math.isclose(row["CL_from_gauge_invariant_extra_observable"], 8.0)
    assert math.isclose(row["same_shear_gate"], 1.0)


def test_gauge_bookkeeping_marks_proxy_noninvariant_and_extra_ratio_invariant():
    powers = gauge_exponents()
    assert powers["kappa_proxy"] == [-1, 1]
    assert powers["CL_invariant_minus_2_t3_squared_over_s2"] == [0, 0]


def test_partial_8_over_3_comparison_is_labelled_nonfinal_conjecture():
    payload = render(2.653, 0.448)
    diagnostic = payload["partial_8_over_3_diagnostic"]
    assert math.isclose(diagnostic["target"], 8.0 / 3.0)
    assert abs(diagnostic["z_score_estimate_minus_target"]) < 0.04
    assert diagnostic["classification"] == "high_risk_amplitude_conjecture_not_exponent_theorem"
    assert payload["status"] == "theory_audit_partial_not_final"
