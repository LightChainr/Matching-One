import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "describe_geometry", ROOT / "scripts" / "describe_p334_geometry_temporal_modes.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_upper_matrix_roundtrip():
    matrix = MODULE.upper_matrix([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    assert np.allclose(matrix, [[1, 2, 3], [2, 4, 5], [3, 5, 6]])


def test_line_centered_geometry_explains_constructed_mode():
    rows, values = [], []
    for line in (0, 1):
        for index in range(1, 8):
            rows.append({
                "ell_u": line, "ell_v": 1,
                "g_size": index / 10,
                "g_carriers": 0.0,
                "g_occupied_frontier": ((3 * index + line) % 7) / 10,
                "g_vacant_frontier": 0.0,
            })
            values.append(4 + line + 2 * rows[-1]["g_size"] - rows[-1]["g_occupied_frontier"])
    result = MODULE.explain(rows, np.asarray(values))
    assert abs(result["cheap_geometry_R2"] - 1.0) < 1e-12
