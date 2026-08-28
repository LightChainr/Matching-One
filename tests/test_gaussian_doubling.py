from __future__ import annotations

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "analyze_gaussian_doubling", ROOT / "scripts" / "analyze_gaussian_doubling.py"
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_p31_retrospective_doubling_score() -> None:
    batches = ROOT / "results/server-20260828/P31/p31_confirmation_seed2026093001.batches.csv"
    result = MODULE.score(MODULE.load_cross_batches([batches]))
    rows = result["lineages"]
    assert len(rows) == 2
    assert math.isclose(result["target_ratio"], -2.0 ** (-13.0 / 8.0), rel_tol=1e-15)
    assert math.isclose(rows[0]["ratio"], -0.3737554822806287, rel_tol=1e-13)
    assert math.isclose(rows[1]["ratio"], -0.23484766130705906, rel_tol=1e-13)
    assert math.isclose(result["joint_residual_chi2"], 1.468190724967981, rel_tol=1e-12)
