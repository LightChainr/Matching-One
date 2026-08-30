import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("mode_analysis", ROOT / "scripts" / "analyze_p334_two_time_modes.py")
mode_analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mode_analysis)


def test_matrix_roundtrip():
    matrix = np.arange(49, dtype=float).reshape(7, 7)
    matrix = matrix + matrix.T
    upper = [matrix[i, j] for i in range(7) for j in range(i, 7)]
    assert np.array_equal(mode_analysis.matrix_from_upper(upper), matrix)


def test_subspace_angles_ignore_internal_rotation():
    first = np.eye(4)
    theta = 0.37
    second = np.eye(4)
    second[:, 1:3] = first[:, 1:3] @ np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)],
    ])
    assert max(mode_analysis.angles(first, second, [1, 2])) < 1e-6
