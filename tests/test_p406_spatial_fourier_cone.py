import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cone", ROOT / "scripts" / "score_p406_spatial_fourier_cone.py")
cone = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cone)


def test_coordinate_parser_and_alias():
    assert cone.decode("am3_bp2_r1_plus_im") == (-3, 2, 1, "plus", "im")
    assert (-3 - 20) % 101 == 78


def test_positive_fourier_data_fit_exactly():
    coordinates = [(0, 0), (1, 0), (0, 1), (2, 1)]
    matrix = cone.design(coordinates)
    weights = np.zeros(101)
    weights[[1, 7, 23]] = [2.0, 0.4, 1.2]
    fitted, distance = cone.fit_nonnegative(matrix, matrix @ weights)
    assert distance < 1e-18
    assert np.linalg.norm(matrix @ fitted - matrix @ weights) < 1e-9


def test_design_respects_group_aliases():
    first = cone.design([(1, 0)])
    second = cone.design([(-9, -1)])
    assert np.max(np.abs(first - second)) < 1e-14
