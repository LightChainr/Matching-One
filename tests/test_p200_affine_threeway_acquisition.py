import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p200", ROOT / "scripts" / "p200_affine_threeway_acquisition.py"
)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_q2_two_path_difference_is_scaled_cocycle():
    x1 = [2.0, -1.0]
    x2 = [3.0, 4.0]
    x5 = [7.0, 5.0]
    left = MOD.forecast_from_norm5_pair(x1, x5, "q2_analytic")
    right = MOD.forecast_from_norm2_pair(x1, x2, "q2_analytic")
    c = MOD.cocycle_multiplier("q2_analytic")
    residual = [a - c * b + (c - 1) * d for d, b, a in zip(x1, x2, x5)]
    assert MOD.path_difference_factor("q2_analytic") == 9 / 8
    assert all(abs((a - b) - 9 / 8 * r) < 1e-14 for a, b, r in zip(left, right, residual))


def test_jordan_two_path_difference_is_scaled_cocycle():
    x1 = [2.0, -1.0]
    x2 = [3.0, 4.0]
    x5 = [7.0, 5.0]
    left = MOD.forecast_from_norm5_pair(x1, x5, "rank2_Jordan")
    right = MOD.forecast_from_norm2_pair(x1, x2, "rank2_Jordan")
    c = MOD.cocycle_multiplier("rank2_Jordan")
    factor = 1 + 1 / c
    residual = [a - c * b + (c - 1) * d for d, b, a in zip(x1, x2, x5)]
    assert all(abs((a - b) - factor * r) < 1e-14 for a, b, r in zip(left, right, residual))


def test_width_corrected_jet():
    state = {
        "canonical_width": "0.5",
        "finite_thermal_jet": ["1", "2", "4", "8", "16", "32", "64"],
    }
    assert MOD.width_corrected_jet(state) == [1.0, 1.0, 1.0, 1.0, 1.0]


def test_exact_endpoint_products():
    def multiply(z, w):
        return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]

    assert multiply(multiply((8, 1), (2, -1)), (1, 1)) == (23, 11)
    assert multiply(multiply((8, 1), (1, 1)), (2, -1)) == (23, 11)
    assert multiply(multiply((7, 4), (2, -1)), (1, 1)) == (17, 19)
    assert multiply(multiply((7, 4), (1, 1)), (2, -1)) == (17, 19)


def test_render_consumes_frozen_covariances():
    payload = MOD.render(
        ROOT / "results/full-curve-transfer/p180_n145_n290_affine_clock.json",
        ROOT / "results/server-20260829/P57-norm5-500m/thermal_jet_score.json",
        ROOT / "results/server-20260829/P57-norm5-500m/conjugation_parity_diagnostic.json",
    )
    assert payload["N580_radial_clock"]["designs"] == [[24, 2], [18, 16]]
    assert payload["N650_mixed_clock"]["endpoint_designs_N650"] == [[23, 11], [17, 19]]
    models = payload["N650_mixed_clock"]["width_corrected_jet_models"]
    assert set(models) == {"q2_analytic", "rank2_Jordan"}
    assert len(models["rank2_Jordan"]["path_difference_covariance"]) == 5
    assert len(payload["N650_mixed_clock"]["marked_commutator"]["frozen_P57_odd_template"]) == 5
