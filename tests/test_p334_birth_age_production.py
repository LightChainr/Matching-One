import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p334_birth_age_production as birth


def table(rows):
    output = {}
    for line, age, count, event in rows:
        birth.add_slope_sufficient(output, (line,), age, count, event)
    return output


def test_line_fixed_effect_removes_line_intercepts():
    # Both lines have different hazards but no within-line age dependence.
    rows = []
    for line, exits in ((0, 2), (1, 8)):
        rows.extend(((line, 1, 10, False), (line, 3, 10, False)))
        # Add the same exit mass at each age within a line.
        output = table(rows)
        output[(line,)][3] += 2 * exits
        output[(line,)][4] += 4 * exits
        rows = []
    combined = {
        (0,): [20, 40, 100, 4, 8],
        (1,): [20, 40, 100, 16, 32],
    }
    assert abs(birth.age_slope(combined, 100)["beta_age_per_density"]) < 1e-12


def test_age_slope_detects_older_exit_concentration():
    data = {
        (0,): [20, 40, 100, 8, 24],  # all eight exits concentrated at age 3
        (1,): [20, 40, 100, 4, 12],
    }
    assert birth.age_slope(data, 100)["beta_age_per_density"] > 0


def test_complement_birth_mapping_preserves_risk_and_direct_mass():
    n, k0 = 25, 15
    k0c = n - k0
    for k1 in range(1, n + 1):
        for k2 in range(k1, n + 1):
            k1c, k2c = n + 1 - k2, n + 1 - k1
            assert (k1 <= k0 < k2) == (k1c <= k0c < k2c)
            assert (k1 == k2) == (k1c == k2c)


def test_paired_jackknife_covariance_retains_cross_coordinate_signal():
    replicates = [[float(index), 2.0 * index] for index in range(1, 11)]
    covariance = birth.jackknife_covariance(replicates)
    assert covariance.shape == (2, 2)
    assert np.isclose(covariance[0, 1], 2.0 * covariance[0, 0])


def test_six_arm_ratio_is_frozen_without_free_exponent():
    ratio = (425 / 325) ** (-5 / 6)
    assert abs(ratio - 0.7996722504772675) < 1e-15
