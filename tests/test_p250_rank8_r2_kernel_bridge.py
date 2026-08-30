import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p250_rank8_r2_kernel_bridge as bridge


def synthetic_rank8(seed: int = 7) -> np.ndarray:
    rng = np.random.default_rng(seed)
    left, _ = np.linalg.qr(rng.normal(size=(20, 10)) + 1j * rng.normal(size=(20, 10)))
    right, _ = np.linalg.qr(rng.normal(size=(10, 10)) + 1j * rng.normal(size=(10, 10)))
    values = np.diag([10.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 0.0, 0.0])
    return left @ values @ right.conjugate().T


def test_projector_is_basis_free_and_rank_two():
    matrix = synthetic_rank8()
    projector = bridge.kernel_projector(matrix)
    assert np.linalg.matrix_rank(projector, tol=1e-10) == 2
    assert np.allclose(projector, projector.conjugate().T, atol=1e-12)
    assert np.allclose(projector @ projector, projector, atol=1e-12)


def test_conjugate_bridge_closes_on_exact_control():
    plus = synthetic_rank8()
    minus = plus.conjugate()
    assert np.linalg.norm(bridge.bridge_matrix(plus, minus)) < 1e-12
    geometry = bridge.descriptive_geometry(plus, minus)
    assert np.allclose(geometry["kernel_plane_principal_cosines"], [1.0, 1.0])


def test_nonmatching_kernel_plane_is_detected():
    plus = synthetic_rank8(7)
    minus = synthetic_rank8(8)
    assert np.linalg.norm(bridge.bridge_matrix(plus, minus)) > 0.1


def test_frozen_hashes_and_upstream_gate_are_live():
    manifest = json.loads((ROOT / "analysis/p250_rank8_r2_kernel_bridge_freeze.json").read_text())
    paths = bridge.checked_inputs(manifest)
    gate = bridge.check_upstream_gate(paths["upstream_score"], manifest["statistic"]["alpha"])
    assert gate["plus"]["rank_lower_bound_at_alpha"] == 8
    assert gate["minus"]["rank_lower_bound_at_alpha"] == 8
