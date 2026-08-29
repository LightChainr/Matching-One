import math

from scripts.score_p234_cross_cutoff_shear import cross_cutoff_score


def test_exact_cross_cutoff_log_shear_is_recovered():
    rows = []
    mixed, kappa = 2.0, 1.5
    for delta in (0.09, 0.06, 0.045):
        for size in (64, 96, 128, 192):
            radius_inverse = 1.0 / (size * delta)
            rows.append(
                {
                    "L": size,
                    "realized_delta": delta,
                    "natural_realized_cutoff_point": [
                        0.4 / size,
                        mixed + 0.3 * radius_inverse,
                        5.0 - 2.0 * kappa * mixed * math.log(2.0 * delta) - 0.7 * radius_inverse,
                    ],
                    "natural_realized_cutoff_covariance": [
                        [0.1, 0.0, 0.0],
                        [0.0, 0.2, 0.0],
                        [0.0, 0.0, 0.3],
                    ],
                }
            )
    score = cross_cutoff_score(rows)
    assert math.isclose(score["gls"]["chi_square"], 0.0, abs_tol=1e-10)
    assert math.isclose(score["kappa_proxy"]["estimate"], kappa)
    assert score["gls"]["degrees_of_freedom"] == 30
