import importlib.util
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


summary = load("p334_summary", ROOT / "scripts" / "summarize_p334_two_time_kernel.py")
score = load("p334_score", ROOT / "scripts" / "score_p334_two_time_kernel.py")


def test_frozen_layers_are_strictly_increasing():
    assert summary.layers(325, 193) == sorted(set(summary.layers(325, 193)))
    assert summary.layers(425, 252) == sorted(set(summary.layers(425, 252)))


def test_two_time_identity_by_enumeration():
    grid = [1, 2, 3]
    for k1 in range(1, 4):
        for k2 in range(k1, 4):
            rank = [int(k1 <= k) + int(k2 <= k) for k in grid]
            for i in range(3):
                for j in range(i, 3):
                    joint = int(k1 <= grid[i] and k2 <= grid[j])
                    rhs = rank[i] * rank[j] - int(k1 <= grid[i]) - 2 * int(k2 <= grid[i])
                    assert joint == rhs


def test_adjacent_defect_detects_rank_one():
    vector = np.arange(1.0, 8.0)
    assert np.max(np.abs(score.adjacent_defects(np.outer(vector, vector)))) < 1e-12
    full = np.eye(7) + np.outer(vector, vector)
    assert np.min(score.adjacent_defects(full)) > 0
