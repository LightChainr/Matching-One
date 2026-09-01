#!/usr/bin/env python3
"""Build and score the honest current-asset column families for Issue #275.

No samples are generated.  The script separates two questions:

1. What part of the six-coordinate restricted-trace thermal jet is fixed by
   the current vacuum/Ward and thermal-Q4/Jordan theory assets?
2. How well separated are the already-frozen semisimple and Jordan *proxy*
   transfer spaces in the existing K1/K2 covariance?

The proxy score is deliberately not promoted to a physical source-column
score: neither candidate currently supplies a same-source six-coordinate
numeric insertion through the original-U normalizer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "p275_restricted_trace_columns_manifest.yaml"
OUTPUT_SCHEMA = "matching-one.p275-restricted-trace-columns.v1"
MANIFEST_SCHEMA = "matching-one.p275-restricted-trace-columns.manifest.v1"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_rank(matrix: np.ndarray, tolerance: float = 1e-11) -> int:
    singular = np.linalg.svd(np.asarray(matrix, dtype=float), compute_uv=False)
    if not singular.size:
        return 0
    return int(np.count_nonzero(singular > tolerance * singular[0]))


def image_intersection_dimension(left: np.ndarray, right: np.ndarray) -> int:
    return matrix_rank(left) + matrix_rank(right) - matrix_rank(np.column_stack([left, right]))


def whiten(covariance: np.ndarray) -> np.ndarray:
    covariance = (covariance + covariance.T) / 2.0
    values, vectors = np.linalg.eigh(covariance)
    cutoff = max(float(np.max(values)) * 1e-11, 1e-30)
    if float(np.min(values)) <= cutoff:
        raise ValueError("expected positive-definite existing covariance")
    return (vectors / np.sqrt(values)).T


def block_design(kappa: float, block_count: int = 4) -> np.ndarray:
    block = np.asarray([[1.0, 1.0], [1.0, kappa], [1.0, kappa * kappa]])
    design = np.zeros((3 * block_count, 2 * block_count))
    for index in range(block_count):
        design[3 * index : 3 * index + 3, 2 * index : 2 * index + 2] = block
    return design


def jordan_design(block_count: int = 4) -> np.ndarray:
    block = np.asarray([[1.0, 0.0], [1.0, 1.0], [1.0, 2.0]])
    design = np.zeros((3 * block_count, 2 * block_count))
    for index in range(block_count):
        design[3 * index : 3 * index + 3, 2 * index : 2 * index + 2] = block
    return design


def orthonormal_image(matrix: np.ndarray) -> tuple[np.ndarray, list[float]]:
    left, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    rank = matrix_rank(matrix)
    return left[:, :rank], [float(value) for value in singular[:rank]]


def gls_fit(estimate: np.ndarray, covariance: np.ndarray, design: np.ndarray) -> dict[str, Any]:
    whitening = whiten(covariance)
    white_y = whitening @ estimate
    white_x = whitening @ design
    coefficients, _, rank, singular = np.linalg.lstsq(white_x, white_y, rcond=None)
    residual = white_y - white_x @ coefficients
    return {
        "design_rank": int(rank),
        "observation_dimension": int(len(estimate)),
        "residual_degrees_of_freedom": int(len(estimate) - rank),
        "mahalanobis_chi_square": float(residual @ residual),
        "coefficients": [float(value) for value in coefficients],
        "whitened_design_singular_values": [float(value) for value in singular],
    }


def restricted_jet_audit() -> dict[str, Any]:
    # raw order: s0,s1,s2,dp_s0,dp_s1,dp_s2
    transport = np.asarray(
        [
            [-0.5, 0.0, 0.5, 0.0, 0.0, 0.0],
            [-0.5, 1.0, -0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, -0.5, 0.0, 0.5],
            [0.0, 0.0, 0.0, -0.5, 1.0, -0.5],
        ]
    )
    common_normalizers = np.asarray(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 1.0],
            [0.0, 1.0],
        ]
    )
    # Critical-surface Ward identity fixes s2=s0, but not its transverse p jet.
    vacuum_critical = np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0],
        ]
    )
    # Stronger neighbourhood identity also fixes dp_s2=dp_s0.
    vacuum_neighbourhood = np.asarray(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
        ]
    )
    thermal_unconstrained = np.eye(6)

    images = {
        "vacuum_Ward_critical_surface": transport @ vacuum_critical,
        "vacuum_Ward_neighbourhood": transport @ vacuum_neighbourhood,
        "thermal_Q4_Jordan_current_envelope": transport @ thermal_unconstrained,
    }
    thermal = images["thermal_Q4_Jordan_current_envelope"]
    return {
        "raw_order": ["s0", "s1", "s2", "partial_p_s0", "partial_p_s1", "partial_p_s2"],
        "normalized_order": ["eta_B", "xi_B", "partial_p_eta_B", "partial_p_xi_B"],
        "transport_matrix": transport.tolist(),
        "transport_rank": matrix_rank(transport),
        "common_normalizer_kernel_dimension": int(common_normalizers.shape[1]),
        "common_normalizer_is_annihilated": bool(np.allclose(transport @ common_normalizers, 0.0)),
        "candidate_families": {
            "vacuum_Ward_critical_surface": {
                "raw_family_dimension": int(vacuum_critical.shape[1]),
                "normalized_image_rank": matrix_rank(images["vacuum_Ward_critical_surface"]),
                "fixed_relation": "eta_B=0_at_p0",
                "unfixed_coordinates": ["xi_B", "partial_p_eta_B", "partial_p_xi_B"],
            },
            "vacuum_Ward_neighbourhood_stronger_case": {
                "raw_family_dimension": int(vacuum_neighbourhood.shape[1]),
                "normalized_image_rank": matrix_rank(images["vacuum_Ward_neighbourhood"]),
                "fixed_relations": ["eta_B=0", "partial_p_eta_B=0"],
                "unfixed_coordinates": ["xi_B", "partial_p_xi_B"],
            },
            "thermal_Q4_Jordan_current_envelope": {
                "raw_family_dimension": int(thermal_unconstrained.shape[1]),
                "normalized_image_rank": matrix_rank(thermal),
                "fixed_restricted_sector_relations": [],
                "available_only_after_transport": ["E4hat_modulus_shape", "Jordan_logN_slope_shape"],
                "unfixed_coordinates": ["eta_B", "xi_B", "partial_p_eta_B", "partial_p_xi_B"],
            },
        },
        "column_space_comparison": {
            "critical_Ward_intersection_with_thermal_dimension": image_intersection_dimension(
                images["vacuum_Ward_critical_surface"], thermal
            ),
            "critical_Ward_image_is_contained_in_thermal_envelope": bool(
                matrix_rank(np.column_stack([thermal, images["vacuum_Ward_critical_surface"]]))
                == matrix_rank(thermal)
            ),
            "neighbourhood_Ward_intersection_with_thermal_dimension": image_intersection_dimension(
                images["vacuum_Ward_neighbourhood"], thermal
            ),
            "thermal_extra_directions_beyond_critical_Ward": int(
                matrix_rank(thermal) - matrix_rank(images["vacuum_Ward_critical_surface"])
            ),
            "thermal_extra_directions_beyond_neighbourhood_Ward": int(
                matrix_rank(thermal) - matrix_rank(images["vacuum_Ward_neighbourhood"])
            ),
        },
    }


def proxy_audit(forward: dict[str, Any]) -> dict[str, Any]:
    estimate = np.asarray(
        [
            row["N_power_13_over_8_scaled_u"]
            for row in forward["observable"]["estimate_by_coordinate"]
        ],
        dtype=float,
    )
    covariance = np.asarray(forward["observable"]["scaled_covariance"], dtype=float)
    semisimple = block_design(0.5)
    jordan = jordan_design()
    whitening = whiten(covariance)
    white_s = whitening @ semisimple
    white_j = whitening @ jordan
    image_s, singular_s = orthonormal_image(white_s)
    image_j, singular_j = orthonormal_image(white_j)
    cosines = np.linalg.svd(image_s.T @ image_j, compute_uv=False)
    intersection = image_intersection_dimension(semisimple, jordan)
    return {
        "label": "conditional_transfer_proxy_not_physical_candidate_columns",
        "coordinate_count": int(len(estimate)),
        "covariance_rank": matrix_rank(covariance),
        "semisimple_kappa_0p5": gls_fit(estimate, covariance, semisimple),
        "Jordan_kappa_1": gls_fit(estimate, covariance, jordan),
        "column_spaces": {
            "semisimple_rank": matrix_rank(semisimple),
            "Jordan_rank": matrix_rank(jordan),
            "ordinary_intersection_dimension": intersection,
            "ordinary_combined_rank": matrix_rank(np.column_stack([semisimple, jordan])),
            "covariance_whitened_principal_cosines": [float(value) for value in cosines],
            "covariance_whitened_semisimple_singular_values": singular_s,
            "covariance_whitened_Jordan_singular_values": singular_j,
        },
        "archive_crosscheck": {
            "semisimple_chi_square": float(
                forward["scores"]["fixed_semisimple_q2_kappa_0p5"]["mahalanobis_chi_square"]
            ),
            "Jordan_chi_square": float(forward["scores"]["Jordan_kappa_1"]["mahalanobis_chi_square"]),
        },
    }


def build(manifest_path: Path) -> dict[str, Any]:
    started = time.perf_counter()
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unexpected manifest schema")

    inputs: dict[str, dict[str, str]] = {}
    for name, relative in manifest["inputs"].items():
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        inputs[name] = {"path": relative, "sha256": sha256_file(path)}

    forward = json.loads((ROOT / manifest["inputs"]["global_k1_k2_covariance"]).read_text())
    rho = json.loads((ROOT / manifest["inputs"]["rho_child_covariance"]).read_text())
    rho_covariance = np.asarray(rho["reconstruction"]["covariance_of_mean"], dtype=float)
    result = {
        "schema": OUTPUT_SCHEMA,
        "issue": 275,
        "status": "PARTIAL_COLUMNS_COMPLETE_EXISTING_COVARIANCES_NOT_DIRECTLY_SCOREABLE",
        "execution": manifest["execution"],
        "inputs": inputs,
        "restricted_trace_jet": restricted_jet_audit(),
        "existing_covariance_coverage": {
            "rho_child": {
                "observation_dimension": int(rho_covariance.shape[0]),
                "covariance_rank": matrix_rank(rho_covariance),
                "same_source_restricted_trace_coordinates_present": 0,
                "candidate_forward_map_rank": 0,
                "direct_candidate_score_status": "UNAVAILABLE",
            },
            "global_K1_K2": {
                "observation_dimension": len(forward["observable"]["estimate_by_coordinate"]),
                "covariance_rank": matrix_rank(np.asarray(forward["observable"]["scaled_covariance"])),
                "same_source_restricted_trace_coordinates_present": 0,
                "candidate_forward_map_rank": 0,
                "direct_candidate_score_status": "UNAVAILABLE",
            },
        },
        "existing_K1_K2_proxy_score": proxy_audit(forward),
        "decision": {
            "literal": "UNIDENTIFIABLE_WITH_CURRENT_ASSETS",
            "stronger_new_detail": "the_critical_surface_Ward_family_has_rank3_inside_the_rank4_thermal_envelope_after_common_normalizers_are_removed",
            "smallest_current_gap": "one_candidate_specific_restricted_sector_direction_at_the_critical_surface_plus_the_numeric_source_scale_and_original_U_map",
            "why_covariance_does_not_close_it": "both_existing_covariances_are_full_rank_but_neither_archive_contains_the_same_B_source_six_coordinate_jet",
            "next_theory_target": "fix_partial_p_eta_B_and_xi_B_partial_p_xi_B_for_each_named_source_then_apply_the_existing_covariance_once",
        },
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "wall_seconds": time.perf_counter() - started,
            "new_random_samples": 0,
        },
        "claim_boundary": manifest["claim_boundary"],
    }
    proxy = result["existing_K1_K2_proxy_score"]
    for key, archive_key in [
        ("semisimple_kappa_0p5", "semisimple_chi_square"),
        ("Jordan_kappa_1", "Jordan_chi_square"),
    ]:
        if not np.isclose(
            proxy[key]["mahalanobis_chi_square"],
            proxy["archive_crosscheck"][archive_key],
            rtol=1e-9,
            atol=1e-9,
        ):
            raise AssertionError(f"proxy GLS replay disagrees for {key}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "status": result["status"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
