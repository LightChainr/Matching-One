import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "score_geometry", ROOT / "scripts" / "score_p334_current_k0_geometry_pilot.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def row(line, age, geometry, outcome, batch=0, orientation="first"):
    return {
        "ell_u": line,
        "ell_v": 1,
        "age": age,
        "g_size": geometry,
        "g_carriers": 0.0,
        "g_occupied_frontier": 0.0,
        "g_vacant_frontier": 0.0,
        "h2_theta_rate": geometry,
        "h2_figure8_rate": 0.0,
        "h2_separate_rate": 0.0,
        "h2_rate": geometry,
        "y": outcome,
        "batch": batch,
        "orientation": orientation,
    }


def test_line_centering_recovers_age_and_geometry_coefficients():
    rows = []
    for line in (0, 1):
        for index in range(1, 8):
            age = index / 10
            geometry = ((index * 3 + line) % 7) / 10
            outcome = 4.0 + line + 2.0 * age - 3.0 * geometry
            rows.append(row(line, age, geometry, outcome))
    fit = MODULE.centered_fit(rows, ("age", "g_size", "g_carriers"))
    assert abs(fit["coefficients"]["age"] - 2.0) < 1e-10
    assert abs(fit["coefficients"]["g_size"] + 3.0) < 1e-10
    assert fit["dropped_zero_information"] == ["g_carriers"]


def test_jackknife_covariance_keeps_paired_direction():
    deleted = [[1.0, -2.0], [2.0, -4.0], [3.0, -6.0]]
    cov = MODULE.covariance(deleted)
    assert cov[0][0] > 0
    assert abs(cov[0][1] + 2 * cov[0][0]) < 1e-12
    assert abs(cov[1][1] - 4 * cov[0][0]) < 1e-12
