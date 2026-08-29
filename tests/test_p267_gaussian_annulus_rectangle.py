from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p267_rectangle", ROOT / "scripts/design_p267_gaussian_annulus_rectangle.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_radial_holdout_rows_reproduce_calibration_basis() -> None:
    coordinate = np.log2(np.asarray(MODULE.RADII, dtype=float) / 2)
    for lam in MODULE.LAMBDAS:
        row, interpolation = MODULE.radial_holdout_row(lam)
        calibration = np.stack([
            MODULE.radial_basis(lam, value) for value in coordinate[:3]
        ])
        assert np.allclose(interpolation @ calibration, MODULE.radial_basis(lam, coordinate[3]))
        assert np.allclose(row[:3], -interpolation)
        assert row[3] == 1


def test_crosswalk_has_no_pseudo_matched_rows() -> None:
    verdict = MODULE.semantic_crosswalk()["exact_pairwise_verdict"]
    assert verdict["same_source_rows"] == []
    assert verdict["same_observer_rows"] == []
    assert verdict["numerically_eligible_cross_context_pairs"] == []


def test_manifest_has_exactly_two_missing_cells() -> None:
    manifest = yaml.safe_load(
        (ROOT / "analysis/p267_gaussian_annulus_missing_cells_20260829.yaml").read_text()
    )
    cells = manifest["rectangle"]["cells"]
    statuses = [entry["status"] for context in cells.values() for entry in context.values()]
    assert statuses.count("existing") == 2
    assert statuses.count("missing") == 2
    assert manifest["production_authorized"] is True


def test_every_frozen_period_matrix_has_declared_determinant() -> None:
    manifest = yaml.safe_load(
        (ROOT / "analysis/p267_gaussian_annulus_missing_cells_20260829.yaml").read_text()
    )
    designs = [
        design
        for lineage in manifest["missing_Gaussian_acquisition"]["cover_lineages"]
        for design in lineage["designs"]
    ]
    assert len(designs) == 16
    for design in designs:
        assert round(np.linalg.det(np.asarray(design["matrix"], dtype=int))) == design["N"]


def test_committed_selection_is_covariance_based_and_weak() -> None:
    payload = __import__("json").loads(
        (ROOT / "results/p267-gaussian-annulus-crosswalk/latest.json").read_text()
    )
    selection = payload["annulus_row_selection"]
    assert selection["selected_existing_annulus_context"] == "N425"
    assert selection["candidates"]["N425"]["joint_minimum_adjacent_mahalanobis_squared"] < 0.05
    assert payload["rectangle"]["numerical_score_now"].startswith("forbidden")
