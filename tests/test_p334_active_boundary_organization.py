import csv
import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "active_boundary", ROOT / "scripts" / "score_p334_active_boundary_organization.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_nested_organization_R2_cannot_decrease():
    rng = np.random.default_rng(91)
    rows = []
    for index in range(100):
        row = {"line": (1, 0), "batch": index % 10, "age": rng.normal(),
               "y": rng.normal(), "k1": 1, "k2": 2}
        for name in MODULE.FEATURES: row[name] = rng.normal()
        for name in MODULE.ORG: row[name] = rng.normal()
        rows.append(row)
    arrays = {("N425", orientation): MODULE.centered_arrays(rows)
              for orientation in MODULE.ORIENTATIONS}
    subspace = MODULE.train_subspace(arrays, list(arrays), 1)
    fit = MODULE.fit_nested(rows, subspace, lambda row: row["y"])
    assert fit["incremental_R2"] >= -1e-12


def test_organization_coordinate_names_are_frozen():
    assert MODULE.ORG == ("axis_anisotropy", "corner_balance",
                          "frontier_components_L", "largest_frontier_arc_L",
                          "frontier_concentration")
