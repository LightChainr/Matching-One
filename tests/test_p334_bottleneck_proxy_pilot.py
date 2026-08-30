import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "bottleneck_score", ROOT / "scripts" / "score_p334_bottleneck_proxy_pilot.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_paired_contrast_covariance_cancels_common_noise():
    order = [["first", "M0"], ["first", "Mshape"],
             ["second", "M0"], ["second", "Mshape"]]
    covariance = [[4, 3, 1, 1], [3, 4, 1, 1], [1, 1, 9, 8], [1, 1, 8, 9]]
    score = {
        "vector_order": order,
        "jackknife_covariance": covariance,
        "fits": {
            "first": {"M0": {"beta_age": -2}, "Mshape": {"beta_age": -1}},
            "second": {"M0": {"beta_age": -3}, "Mshape": {"beta_age": -1}},
        },
    }
    delta, cov = MODULE.contrast_covariance(score, "Mshape")
    assert delta == [1, 2]
    assert cov == [[2, 0], [0, 2]]


def test_collinear_frozen_proxy_span_keeps_age_identifiable():
    rows = []
    for line in (0, 1):
        for index in range(1, 8):
            age = index / 10
            proxy = ((index * 3 + line) % 7) / 10
            rows.append({
                "ell_u": line, "ell_v": 1, "age": age,
                "g_size": proxy, "duplicate": 2 * proxy,
                "y": 3 + line - 4 * age + 5 * proxy,
            })
    fit = MODULE.centered_fit(rows, ("age", "g_size", "duplicate"))
    assert abs(fit["coefficients"]["age"] + 4) < 1e-10
    assert fit["collinear_deficiency"] == 1
