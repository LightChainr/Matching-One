from fractions import Fraction

from scripts.p14_w5_terminal_duality import build_certificate, evaluate


def result():
    return build_certificate()


def test_exact_configuration_and_bernstein_totals():
    payload = result()
    distribution = payload["primal_partition_distribution"]
    assert distribution["configuration_count"] == 256
    assert distribution["partition_count_including_zero"] == 15
    assert distribution["identically_zero_partitions"] == ["01|23"]
    assert sum(
        sum(map(sum, table)) for table in distribution["bivariate_bernstein_counts"].values()
    ) == 256


def test_spherical_duality_is_orbit_swapping_involution_but_not_partition_map():
    dual = result()["spherical_dual"]
    assert dual["configuration_map_is_involution"] is True
    assert dual["parameter_map"] == "(r,s)->(1-s,1-r)"
    assert dual["terminal_partition_is_function_of_primal_partition"] is False
    witness = dual["smallest_nonclosure_witness"]
    assert witness["max_open_edge_count"] == 3
    assert witness["left"]["primal_partition"] == witness["right"]["primal_partition"]
    assert witness["left"]["spherical_dual_partition"] != witness["right"]["spherical_dual_partition"]


def test_disk_relative_dual_closes_partition_complement_but_is_not_w5():
    dual = result()["disk_relative_dual"]
    assert dual["relative_partition_is_function_of_primal_partition"] is True
    assert dual["exact_planar_complement_map"]["0123"] == "0|1|2|3"
    assert dual["exact_planar_complement_map"]["0|1|2|3"] == "0123"
    assert dual["boundary_terminal_preserving_isomorphic_to_primal"] is False
    assert dual["primal_degree_multiset"] == [3, 3, 3, 3, 4]
    assert dual["dual_degree_multiset"] == [1, 1, 1, 1, 3, 3, 3, 3]
    assert dual["primal_terminal_degrees"] == [3, 3, 3, 3]
    assert dual["dual_terminal_degrees"] == [1, 1, 1, 1]


def test_all_none_is_not_a_duality_identity_on_natural_line():
    diagnostic = result()["natural_line_scalar_diagnostic"]
    assert diagnostic["identically_zero"] is False
    assert diagnostic["value_at_r_equals_s_equals_one_half"] == "67/128"
    polynomial = diagnostic["P_all_minus_P_none_power_coefficients_low_to_high"]
    assert evaluate(polynomial, Fraction(1, 2)) == Fraction(67, 128)


def test_configuration_certificate_hash_is_stable():
    assert result()["configuration_rows_sha256"] == "a5064e695651527e2ccf196986e36d3aae67fb99a058ebe4b815dbe63a1f927f"
