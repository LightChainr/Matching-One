import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SCRIPT = ROOT / "scripts" / "score_p334_morphology_state_transport.py"
SPEC = importlib.util.spec_from_file_location("transport", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_physical_feature_groups_partition_columns():
    flattened = [name for names in MODULE.GROUPS.values() for name in names]
    assert sorted(flattened) == sorted(MODULE.FEATURES)


def test_rank_one_subspace_recovers_shared_direction():
    rng = np.random.default_rng(1234)
    direction = rng.normal(size=len(MODULE.FEATURES))
    rows = {}
    for index, key in enumerate((("N325", "first"), ("N325", "second"))):
        x = rng.normal(size=(400, len(MODULE.FEATURES)))
        age = rng.normal(size=400)
        y = 0.4 * age + (1.0 + index) * (x @ direction)
        rows[key] = (age, y, x)
    trained = MODULE.train_subspace(rows, list(rows), 1)
    assert trained["rank_energy"] > 0.999999999
    assert trained["morphology_rank"] == len(MODULE.FEATURES)


def test_chi_square_uses_covariance_rank():
    statistic, rank, p = MODULE.chi_square([1.0, 1.0], [[1.0, 1.0], [1.0, 1.0]])
    assert rank == 1
    assert abs(statistic - 1.0) < 1e-10
    assert 0.3 < p < 0.4
