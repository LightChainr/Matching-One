#!/usr/bin/env python3
"""Audit the P537 N25/N65 positive-exposure unit mismatch.

The audit reads immutable Git blobs named in a frozen manifest.  It does not
modify the historical scale fingerprint.  The only numerical transformation
is to put the N25 and N65 positive exposures into the same P/N**a convention
before forming conditional signed density and two-point decay coordinates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Mapping

import yaml


SCHEMA = "matching-one.p537-exposure-unit-audit.v1"
MANIFEST_SCHEMA = "matching-one.p537-exposure-unit-audit.manifest.v1"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def git_blob(repo: Path, commit: str, path: str) -> bytes:
    spec = f"{commit}:{path}"
    completed = subprocess.run(
        ["git", "-C", str(repo), "show", spec],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def checked_inputs(repo: Path, manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    checked: dict[str, dict[str, Any]] = {}
    for name, spec in manifest["git_inputs"].items():
        commit = str(spec["commit"])
        path = str(spec["path"])
        data = git_blob(repo, commit, path)
        observed_blob = git_blob_sha1(data)
        expected_blob = str(spec["git_blob_sha1"])
        if observed_blob != expected_blob:
            raise ValueError(
                f"{name}: Git blob mismatch {observed_blob}; expected {expected_blob}"
            )
        text = data.decode("utf-8")
        marker_counts: dict[str, int] = {}
        for marker in spec.get("required_markers", []):
            count = text.count(str(marker))
            if count != 1:
                raise ValueError(f"{name}: required marker count is {count}: {marker!r}")
            marker_counts[str(marker)] = count
        checked[name] = {
            "commit": commit,
            "path": path,
            "git_blob_sha1": observed_blob,
            "sha256": sha256_bytes(data),
            "required_marker_counts": marker_counts,
            "data": data,
        }
    return checked


def as_json(source: Mapping[str, Any]) -> dict[str, Any]:
    value = json.loads(source["data"])
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object from {source['path']}")
    return value


def close(a: float, b: float, *, rel: float = 2e-12, abs_: float = 1e-18) -> bool:
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


def decay(value25: float, value65: float, ratio: float) -> float:
    if value25 == 0 or value65 == 0 or ratio <= 1:
        raise ValueError("Decay coordinate requires nonzero values and size ratio > 1")
    return -math.log(abs(value65 / value25)) / math.log(ratio)


def convert_exposure(value: float, n: int, stored_power: int, target_power: int) -> float:
    """Convert stored P/N**stored_power to P/N**target_power."""
    if value <= 0 or n <= 0:
        raise ValueError("Positive exposure and system size are required")
    return value * n ** (stored_power - target_power)


def convention_block(
    *,
    name: str,
    target_power: int,
    n25: int,
    n65: int,
    exposure25_stored: float,
    exposure65_stored: float,
    stored_power25: int,
    stored_power65: int,
    signed25: float,
    signed65: float,
) -> dict[str, Any]:
    ratio = n65 / n25
    exposure25 = convert_exposure(exposure25_stored, n25, stored_power25, target_power)
    exposure65 = convert_exposure(exposure65_stored, n65, stored_power65, target_power)
    density25 = signed25 / exposure25
    density65 = signed65 / exposure65
    exposure_decay = decay(exposure25, exposure65, ratio)
    density_decay = decay(density25, density65, ratio)
    signed_decay = decay(signed25, signed65, ratio)
    closure = exposure_decay + density_decay - signed_decay
    if abs(closure) > 2e-12:
        raise AssertionError(f"{name}: decay closure failed: {closure}")
    return {
        "id": name,
        "exposure_coordinate": f"P/N^{target_power}",
        "target_N_denominator_power": target_power,
        "N25": {
            "exposure": exposure25,
            "conditional_signed_density": density25,
            "signed_mass": signed25,
        },
        "N65": {
            "exposure": exposure65,
            "conditional_signed_density": density65,
            "signed_mass": signed65,
        },
        "two_point_decay": {
            "exposure": exposure_decay,
            "conditional_signed_density": density_decay,
            "signed_mass": signed_decay,
            "additive_closure_residual": closure,
        },
    }


def public_provenance(checked: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    return {
        name: {key: value for key, value in source.items() if key != "data"}
        for name, source in checked.items()
    }


def build_result(repo: Path, manifest_path: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    checked = checked_inputs(repo, manifest)
    scale = as_json(checked["historical_scale_result"])
    n65_result = as_json(checked["n65_result"])
    full_t = as_json(checked["full_t_result"])

    if scale.get("schema") != "matching-one/p537-contact-stage-scale/v1":
        raise ValueError("Unexpected historical scale schema")
    if scale.get("status") != "POST_N65_TWO_SCALE_MECHANISM_FINGERPRINT":
        raise ValueError("Historical scale result is not the pinned post-reveal diagnostic")
    if n65_result.get("schema") != "matching-one/p537-contact-stage-n65-score/v1":
        raise ValueError("Unexpected N65 result schema")
    if full_t.get("schema") != "matching-one/p537-full-t-transport/v1":
        raise ValueError("Unexpected full-T result schema")

    contract = manifest["coordinate_contract"]
    n25 = int(contract["sizes"]["n25"])
    n65 = int(contract["sizes"]["n65"])
    row_index = int(contract["selected_cell"]["row_index"])
    column_index = int(contract["selected_cell"]["column_index"])
    stored25 = int(contract["stored_exposure_denominator_power"]["n25"])
    stored65 = int(contract["stored_exposure_denominator_power"]["n65"])
    ratio = n65 / n25

    signed25 = float(scale["N25_matrix"][row_index][column_index])
    signed65 = float(scale["N65_matrix"][row_index][column_index])
    old = scale["entry_double_decomposition"]
    exposure25_stored = float(old["N25_exposure"])
    exposure65_stored = float(old["N65_exposure"])

    mixed_density25 = signed25 / exposure25_stored
    mixed_density65 = signed65 / exposure65_stored
    mixed_exposure_decay = decay(exposure25_stored, exposure65_stored, ratio)
    mixed_density_decay = decay(mixed_density25, mixed_density65, ratio)
    signed_decay = decay(signed25, signed65, ratio)
    if not close(mixed_exposure_decay, float(old["exposure_decay_power"])):
        raise AssertionError("Historical exposure exponent was not reproduced")
    if not close(-mixed_density_decay, float(old["conditional_density_growth_power"])):
        raise AssertionError("Historical conditional-density growth was not reproduced")
    if not close(signed_decay, float(old["signed_mass_decay_power"])):
        raise AssertionError("Historical signed-mass exponent was not reproduced")

    conventions = []
    for name, target_power in contract["conventions"].items():
        conventions.append(
            convention_block(
                name=str(name),
                target_power=int(target_power),
                n25=n25,
                n65=n65,
                exposure25_stored=exposure25_stored,
                exposure65_stored=exposure65_stored,
                stored_power25=stored25,
                stored_power65=stored65,
                signed25=signed25,
                signed65=signed65,
            )
        )

    n65_matrix = [[float(value) for value in row] for row in n65_result["primary"]["matrix"]]
    if any(
        not close(n65_matrix[i][j], float(scale["N65_matrix"][i][j]))
        for i in range(2)
        for j in range(2)
    ):
        raise AssertionError("Scale fingerprint and N65 primary matrix differ")
    determinant = n65_matrix[0][0] * n65_matrix[1][1] - n65_matrix[0][1] * n65_matrix[1][0]
    reported_determinant = float(n65_result["primary"]["Delta"])
    if not close(determinant, reported_determinant, rel=2e-14, abs_=1e-28):
        raise AssertionError("N65 determinant did not close")

    full = full_t["N65_pooled_root"]
    full_t_value = float(full["T_t"]["value"])
    full_j = float(full["J"]["value"])
    selected_t = sum(sum(row) for row in n65_matrix)
    selected_j = selected_t * full_j / full_t_value
    complement_j = full_j - selected_j

    try:
        manifest_label = str(manifest_path.relative_to(repo))
    except ValueError:
        manifest_label = str(manifest_path)

    return {
        "schema": SCHEMA,
        "status": "COMPLETED_FIXED_BLOB_EXISTING_DATA_CORRECTION",
        "issue": 537,
        "manifest": {
            "path": manifest_label,
            "sha256": sha256_file(manifest_path),
            "schema": manifest["schema"],
        },
        "provenance": {
            "git_inputs": public_provenance(checked),
            "acquisition": "fixed_git_blobs_only_zero_new_MC",
        },
        "semantic_audit": {
            "n25_positive_exposure": "P",
            "n25_positive_exposure_N_denominator_power": stored25,
            "n65_positive_exposure": "P/N",
            "n65_positive_exposure_N_denominator_power": stored65,
            "historical_operation": "direct_N25_P_vs_N65_P_over_N_comparison",
            "mismatch_detected": True,
            "N65_factor_to_unweighted_P": n65,
        },
        "selected_cell": {
            "row": str(contract["selected_cell"]["row"]),
            "column": str(contract["selected_cell"]["column"]),
            "historical_mixed_units": {
                "N25_exposure_P": exposure25_stored,
                "N65_exposure_P_over_N": exposure65_stored,
                "N25_conditional_signed_density": mixed_density25,
                "N65_conditional_signed_density": mixed_density65,
                "two_point_exposure_decay": mixed_exposure_decay,
                "two_point_conditional_density_decay": mixed_density_decay,
                "historical_reported_conditional_density_growth": -mixed_density_decay,
                "two_point_signed_mass_decay": signed_decay,
            },
            "consistent_conventions": conventions,
            "scope": "finite_N25_to_N65_comparison_not_a_universal_exponent",
        },
        "frozen_unchanged": {
            "reason": "unit_repair_changes_only_the_positive_exposure_coordinate_and_derived_density",
            "N25_signed_matrix": scale["N25_matrix"],
            "N65_signed_matrix": n65_matrix,
            "N65_determinant": determinant,
            "N65_reported_determinant": reported_determinant,
            "N65_transmission_decision": n65_result["primary"]["decision"],
            "full_J65": {
                "value": full_j,
                "se": float(full["J"]["se"]),
                "T_t": full_t_value,
            },
        },
        "coverage_accounting": {
            "selected_T": selected_t,
            "selected_J_point": selected_j,
            "selected_over_full_J_point": selected_j / full_j,
            "complement_of_selected_J_point": complement_j,
            "complement_of_selected_over_full_J_point": complement_j / full_j,
            "semantic_name": "complement_of_selected",
            "not_identified_as": "spatially_nonlocal_mechanism_or_causal_operator_share",
            "joint_uncertainty": "not_scored_without_shared_full_selected_delete_one_factors",
        },
        "interpretation_corrections": {
            "five_eighth": {
                "historical_split_power": scale["split_power_model"]["powers_in_N_row_major"],
                "historical_Q": float(scale["split_power_model"]["Q"]),
                "historical_nominal_p": float(scale["split_power_model"]["p_value"]),
                "status": "POST_REVEAL_SHAPE_FINGERPRINT_ONLY",
                "withdraw": "additional_five_eighth_assigned_to_conditional_signed_strength",
                "does_not_identify": "field_exponent_or_occurrence_frequency_exponent",
            },
            "commutator": {
                "status": "ALGEBRAIC_ENCODING_CANDIDATE",
                "supported": "signed_stage_by_contact_table_is_nonseparable",
                "not_identified": "physical_noncommuting_birth_and_contact_operations_without_independent_F_B_definitions",
            },
            "remainder": {
                "status": "COMPLEMENT_OF_SELECTED_ONLY",
                "not_identified": "spatial_nonlocality_or_unexplained_causal_fraction",
            },
        },
        "boundary": [
            "no_historical_result_overwritten",
            "no_raw_TSV_or_delete_one_reconstruction",
            "no_new_random_samples_GPU_or_cloud",
            "all_N65_derivatives_share_the_existing_dependency_block",
            "priority_is_attention_not_permission_or_lock",
        ],
    }


def format_markdown(result: Mapping[str, Any]) -> str:
    selected = result["selected_cell"]
    historical = selected["historical_mixed_units"]
    blocks = {block["id"]: block for block in selected["consistent_conventions"]}
    p_block = blocks["unweighted_positive_mass"]
    pn_block = blocks["source_weighted_positive_mass"]
    frozen = result["frozen_unchanged"]
    coverage = result["coverage_accounting"]
    five = result["interpretation_corrections"]["five_eighth"]
    return f"""# P537 exposure-unit correction sidecar

## Decision

The N25 scorer stores the selected positive exposure as `P`; the N65 scorer
stores it as `P/N`.  The historical scale diagnostic divided these coordinates
directly.  This sidecar puts both sizes into one declared convention before
forming conditional signed density.  It does not overwrite the frozen result.

## Entry / double-contact cell

| convention | N25 exposure | N65 exposure | exposure decay | density decay | signed decay |
|:--|--:|--:|--:|--:|--:|
| historical mixed `P` vs `P/N` | {historical['N25_exposure_P']:.12g} | {historical['N65_exposure_P_over_N']:.12g} | {historical['two_point_exposure_decay']:.9f} | {historical['two_point_conditional_density_decay']:.9f} | {historical['two_point_signed_mass_decay']:.9f} |
| common unweighted `P` | {p_block['N25']['exposure']:.12g} | {p_block['N65']['exposure']:.12g} | {p_block['two_point_decay']['exposure']:.9f} | {p_block['two_point_decay']['conditional_signed_density']:.9f} | {p_block['two_point_decay']['signed_mass']:.9f} |
| common source-weighted `P/N` | {pn_block['N25']['exposure']:.12g} | {pn_block['N65']['exposure']:.12g} | {pn_block['two_point_decay']['exposure']:.9f} | {pn_block['two_point_decay']['conditional_signed_density']:.9f} | {pn_block['two_point_decay']['signed_mass']:.9f} |

For common unweighted `P`, the conditional signed densities are
`{p_block['N25']['conditional_signed_density']:.12g}` at N25 and
`{p_block['N65']['conditional_signed_density']:.12g}` at N65.  The repaired
finite-pair decomposition is therefore approximately
`{p_block['two_point_decay']['exposure']:.6f} + {p_block['two_point_decay']['conditional_signed_density']:.6f} = {p_block['two_point_decay']['signed_mass']:.6f}`.

## What remains frozen

- N65 transmission decision: `{frozen['N65_transmission_decision']}`.
- N65 determinant: `{frozen['N65_determinant']:.15g}`.
- Full `J65`: `{frozen['full_J65']['value']:.12g} +/- {frozen['full_J65']['se']:.12g}`.
- The signed N25/N65 matrices and the post-reveal shape score are numerically unchanged.

The repair changes the exposure/density attribution, not signed mass.  The old
`[3, 29/8, 3, 3]` score (`Q={five['historical_Q']:.9f}`, nominal
`p={five['historical_nominal_p']:.9f}`) remains a post-reveal shape fingerprint.
Its additional `5/8` is no longer attributed to conditional signed strength and
does not identify a field or occurrence-frequency exponent.

## Mechanism language

- A nonzero stage-by-contact determinant establishes a nonseparable signed table.
  Calling it a physical commutator still requires independently defined operations
  `F`, `B`, `FB`, and `BF`; the present construction is an algebraic encoding candidate.
- The selected cells account for `{coverage['selected_over_full_J_point']:.6%}` of
  the full `J65` point estimate.  The rest is named only
  `complement_of_selected`, not a spatially nonlocal mechanism or causal share.
- Exact joint uncertainty for that share awaits the shared full/selected delete-one
  factors; none is invented here.

## Provenance and boundary

- Manifest: `{result['manifest']['path']}` (`{result['manifest']['sha256']}`).
- Every input is read from a pinned Git blob and its Git blob SHA-1 is verified.
- Existing fixed data only; no raw TSV replay, new random samples, GPU, cloud job,
  full test suite, or historical-result overwrite.
- These are two-size finite coordinates, not universal exponents.  Priority is
  attention allocation, not permission or a task lock.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()

    repo = Path(
        subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )
    manifest_path = args.manifest if args.manifest.is_absolute() else repo / args.manifest
    manifest = yaml.safe_load(manifest_path.read_text())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"Unexpected manifest schema: {manifest.get('schema')!r}")

    result = build_result(repo, manifest_path, manifest)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.output_md.write_text(format_markdown(result))
    print(
        json.dumps(
            {
                "status": result["status"],
                "mismatch_detected": result["semantic_audit"]["mismatch_detected"],
                "unweighted_decay": result["selected_cell"]["consistent_conventions"][0]["two_point_decay"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
