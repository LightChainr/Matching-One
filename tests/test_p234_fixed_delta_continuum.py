import math

from scripts.score_p234_fixed_delta_continuum import (
    gls_linear_continuum,
    gls_parent_pair,
    transform,
    transform_jacobian,
)


def test_transform_matches_field_normalization():
    result = transform([2.0, 3.0, 5.0, 0.25], 4.0, 2.0)
    assert result == [32.0, 96.0, 320.0]


def test_transform_jacobian_connection_terms():
    jacobian = transform_jacobian([2.0, 3.0, 5.0, 0.25], 4.0, 2.0)
    assert jacobian[0] == [16.0, 0.0, 0.0, 0.0]
    assert jacobian[1] == [0.0, 32.0, 0.0, -384.0]
    assert jacobian[2] == [0.0, 0.0, 64.0, -2560.0]


def test_exact_parent_pair_collapse_has_zero_chi_square():
    points = [[0.0, 2.0, 5.0], [0.0, 2.0, 5.0], [0.0, 2.0, 5.0]]
    covariances = [
        [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]]
        for _ in points
    ]
    score = gls_parent_pair(points, covariances)
    assert math.isclose(score["coefficients"][0], 2.0)
    assert math.isclose(score["coefficients"][1], 5.0)
    assert math.isclose(score["chi_square"], 0.0, abs_tol=1e-12)
    assert score["degrees_of_freedom"] == 7


def test_bottom_null_contributes_to_joint_score():
    covariance = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]]
    score = gls_parent_pair([[1.0, 2.0, 5.0], [1.0, 2.0, 5.0]], [covariance, covariance])
    assert math.isclose(score["chi_square"], 20.0)
    assert score["degrees_of_freedom"] == 4


def test_linear_continuum_recovers_exact_intercepts():
    sizes = [10.0, 20.0, 40.0]
    points = [[3.0 / L, 2.0 + 4.0 / L, 5.0 - 2.0 / L] for L in sizes]
    covariance = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.3]]
    score = gls_linear_continuum(sizes, points, [covariance] * 3)
    assert math.isclose(score["coefficients"][1], 2.0)
    assert math.isclose(score["coefficients"][3], 5.0)
    assert math.isclose(score["chi_square"], 0.0, abs_tol=1e-11)
    assert score["degrees_of_freedom"] == 4
