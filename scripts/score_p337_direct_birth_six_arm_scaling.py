#!/usr/bin/env python3
"""Score the frozen P337 four-generation direct-birth scaling line."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, norm


BETA_FIXED = 5.0 / 6.0
ORIENTATIONS = ("first", "second")
EXTERNAL_COMMIT = "2e99533"
EXTERNAL_PATH = "results/local-20260830/P334-birth-age-production/score.json"
EXTERNAL_SHA256 = "4ab6ddf989b8cce5dad3365d5d440b387a931239c8b88a0d5c62f3e508c3ab29"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def uncompressed_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with gzip.open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path, compression: str | None):
    if compression == "gzip":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def extract_batches(path: Path, spec: Mapping[str, object]) -> dict:
    """Reduce a sparse birth table to paired direct-event batch counts."""
    compression = spec.get("compression")
    batches = int(spec["batches"])
    expected_samples = int(spec["samples_per_orientation"])
    direct = np.zeros((batches, 2), dtype=np.int64)
    totals = np.zeros((batches, 2), dtype=np.int64)
    declared = np.zeros((batches, 2), dtype=np.int64)
    kinds: set[str] = set()
    rows_seen = 0
    with open_text(path, str(compression) if compression else None) as handle:
        reader = csv.DictReader(handle)
        expected_fields = {
            "n", "orientation", "batch", "samples", "kind", "count",
        }
        if not expected_fields.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: sparse birth schema changed")
        for row in reader:
            rows_seen += 1
            if int(row["n"]) != int(spec["N"]):
                raise ValueError(f"{path}: wrong N")
            orientation = row["orientation"]
            if orientation not in ORIENTATIONS:
                raise ValueError(f"{path}: unknown orientation {orientation}")
            oi = ORIENTATIONS.index(orientation)
            batch = int(row["batch"])
            if not 0 <= batch < batches:
                raise ValueError(f"{path}: batch out of range")
            count = int(row["count"])
            samples = int(row["samples"])
            if declared[batch, oi] not in (0, samples):
                raise ValueError(f"{path}: inconsistent samples field")
            declared[batch, oi] = samples
            totals[batch, oi] += count
            kinds.add(row["kind"])
            if row["kind"] == "DIRECT_RANK2":
                direct[batch, oi] += count
    if np.any(declared <= 0):
        raise ValueError(f"{path}: missing orientation/batch cell")
    if not np.array_equal(totals, declared):
        maximum = int(np.max(np.abs(totals - declared)))
        raise ValueError(f"{path}: sparse rows do not partition paths, max error {maximum}")
    if int(declared[:, 0].sum()) != expected_samples or int(declared[:, 1].sum()) != expected_samples:
        raise ValueError(f"{path}: samples per orientation changed")
    if "DIRECT_RANK2" not in kinds or "LINE" not in kinds:
        raise ValueError(f"{path}: required birth types absent")
    return {
        "direct": direct,
        "samples": declared,
        "rows_seen": rows_seen,
        "kinds": sorted(kinds),
    }


def jackknife_covariance(deleted: Sequence[Sequence[float]]) -> np.ndarray:
    values = np.asarray(deleted, dtype=float)
    center = values.mean(axis=0)
    centered = values - center
    return (len(values) - 1) / len(values) * centered.T @ centered


def summarize_size(reduced: Mapping[str, object], spec: Mapping[str, object]) -> dict:
    direct = np.asarray(reduced["direct"], dtype=float)
    samples = np.asarray(reduced["samples"], dtype=float)
    full = direct.sum(axis=0) / samples.sum(axis=0)
    deleted = np.asarray([
        (direct.sum(axis=0) - direct[index]) / (samples.sum(axis=0) - samples[index])
        for index in range(len(direct))
    ])
    covariance = jackknife_covariance(deleted)
    weight = np.asarray([0.5, 0.5])
    dbar = float(weight @ full)
    dbar_deleted = deleted @ weight
    dbar_variance = float(weight @ covariance @ weight)
    dbar_log_variance = dbar_variance / (dbar * dbar)
    return {
        "N": int(spec["N"]),
        "samples_per_orientation": int(spec["samples_per_orientation"]),
        "batches": int(spec["batches"]),
        "direct_counts": {
            orientation: int(np.asarray(reduced["direct"])[:, index].sum())
            for index, orientation in enumerate(ORIENTATIONS)
        },
        "D_by_orientation": {
            orientation: float(full[index]) for index, orientation in enumerate(ORIENTATIONS)
        },
        "orientation_delete_one_covariance": covariance.tolist(),
        "orientation_correlation": float(
            covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])
        ),
        "Dbar": dbar,
        "Dbar_standard_error": math.sqrt(dbar_variance),
        "log_Dbar_variance_delta": dbar_log_variance,
        "Dbar_delete_one": dbar_deleted.tolist(),
        "input_audit": {
            "sparse_rows": int(reduced["rows_seen"]),
            "kinds": list(reduced["kinds"]),
            "partition_gate": True,
        },
    }


def gls(y: np.ndarray, covariance: np.ndarray, design: np.ndarray) -> dict:
    inverse = np.linalg.pinv(covariance, rcond=1e-13)
    normal = design.T @ inverse @ design
    parameter_covariance = np.linalg.inv(normal)
    coefficients = parameter_covariance @ design.T @ inverse @ y
    fitted = design @ coefficients
    residual = y - fitted
    statistic = float(residual @ inverse @ residual)
    degrees = len(y) - design.shape[1]
    return {
        "coefficients": coefficients,
        "parameter_covariance": parameter_covariance,
        "fitted": fitted,
        "residual": residual,
        "chi_square": statistic,
        "degrees_of_freedom": degrees,
        "survival_p": float(chi2.sf(statistic, degrees)) if degrees else None,
    }


def model_scores(sizes: Sequence[dict], alpha: float) -> dict:
    n = np.asarray([row["N"] for row in sizes], dtype=float)
    y = np.log(np.asarray([row["Dbar"] for row in sizes], dtype=float))
    variance = np.diag([row["log_Dbar_variance_delta"] for row in sizes])
    logn = np.log(n)
    generation = np.log2(n / n[0])

    fixed_y = y + BETA_FIXED * logn
    fixed = gls(fixed_y, variance, np.ones((len(n), 1)))
    fixed_a = float(fixed["coefficients"][0])
    fixed_var_a = float(fixed["parameter_covariance"][0, 0])
    fixed_prediction = fixed_a - BETA_FIXED * logn
    fixed["amplitude"] = math.exp(fixed_a)
    fixed["amplitude_standard_error_delta"] = math.exp(fixed_a) * math.sqrt(fixed_var_a)
    fixed["log_prediction"] = fixed_prediction
    fixed["prediction"] = np.exp(fixed_prediction)
    fixed["standardized_residuals_marginal"] = (
        (y - fixed_prediction) / np.sqrt(np.diag(variance))
    )
    fixed["decision"] = "survives" if fixed["survival_p"] >= alpha else "rejected"

    free = gls(y, variance, np.column_stack((np.ones(len(n)), -logn)))
    free_a, free_beta = (float(value) for value in free["coefficients"])
    free["amplitude"] = math.exp(free_a)
    free["beta"] = free_beta
    free["beta_standard_error"] = math.sqrt(float(free["parameter_covariance"][1, 1]))
    free["beta_difference_from_5_6"] = free_beta - BETA_FIXED
    free["beta_difference_z"] = (
        (free_beta - BETA_FIXED) / free["beta_standard_error"]
    )
    free["beta_difference_two_sided_p"] = float(
        2 * norm.sf(abs(free["beta_difference_z"]))
    )

    curvature_column = (generation - 1.5) ** 2
    curved = gls(
        y, variance, np.column_stack((np.ones(len(n)), -logn, curvature_column)),
    )
    curved_a, curved_beta, kappa = (float(value) for value in curved["coefficients"])
    kappa_se = math.sqrt(float(curved["parameter_covariance"][2, 2]))
    curved["amplitude"] = math.exp(curved_a)
    curved["beta"] = curved_beta
    curved["kappa"] = kappa
    curved["kappa_standard_error"] = kappa_se
    curved["kappa_z"] = kappa / kappa_se
    curved["kappa_two_sided_p"] = float(2 * norm.sf(abs(kappa / kappa_se)))
    curved["curvature_resolved_at_alpha"] = curved["kappa_two_sided_p"] < alpha

    fixed_vs_free = fixed["chi_square"] - free["chi_square"]
    free_vs_curved = free["chi_square"] - curved["chi_square"]
    nested = {
        "fixed_5_6_vs_free_beta": {
            "delta_chi_square": fixed_vs_free,
            "degrees_of_freedom": 1,
            "survival_p": float(chi2.sf(fixed_vs_free, 1)),
        },
        "free_beta_vs_log_quadratic": {
            "delta_chi_square": free_vs_curved,
            "degrees_of_freedom": 1,
            "survival_p": float(chi2.sf(free_vs_curved, 1)),
        },
    }

    contrast = np.zeros((3, 4))
    for index in range(3):
        contrast[index, index] = -1.0
        contrast[index, index + 1] = 1.0
    doubling_residual = contrast @ y + BETA_FIXED * math.log(2.0)
    doubling_covariance = contrast @ variance @ contrast.T
    inverse = np.linalg.inv(doubling_covariance)
    doubling_chi = float(doubling_residual @ inverse @ doubling_residual)
    ratios = []
    for index in range(3):
        variance_one = float(doubling_covariance[index, index])
        log_residual = float(doubling_residual[index])
        z = log_residual / math.sqrt(variance_one)
        ratios.append({
            "from_N": int(n[index]),
            "to_N": int(n[index + 1]),
            "observed_ratio": float(math.exp(y[index + 1] - y[index])),
            "fixed_ratio": float(2 ** (-BETA_FIXED)),
            "log_ratio_residual": log_residual,
            "standard_error": math.sqrt(variance_one),
            "z": z,
            "two_sided_p": float(2 * norm.sf(abs(z))),
        })
    doubling = {
        "ratios": ratios,
        "contrast_covariance": doubling_covariance.tolist(),
        "joint_chi_square": doubling_chi,
        "degrees_of_freedom": 3,
        "survival_p": float(chi2.sf(doubling_chi, 3)),
    }

    source = gls(
        y[:3] + BETA_FIXED * logn[:3], variance[:3, :3], np.ones((3, 1)),
    )
    heldout_prediction = float(source["coefficients"][0] - BETA_FIXED * logn[3])
    heldout_variance = float(source["parameter_covariance"][0, 0] + variance[3, 3])
    heldout_residual = float(y[3] - heldout_prediction)
    heldout_z = heldout_residual / math.sqrt(heldout_variance)
    heldout = {
        "source_sizes": [int(value) for value in n[:3]],
        "heldout_N": int(n[3]),
        "source_amplitude": math.exp(float(source["coefficients"][0])),
        "predicted_Dbar": math.exp(heldout_prediction),
        "observed_Dbar": math.exp(float(y[3])),
        "log_residual": heldout_residual,
        "standard_error": math.sqrt(heldout_variance),
        "z": heldout_z,
        "two_sided_p": float(2 * norm.sf(abs(heldout_z))),
        "retrospective_heldout_by_analysis": True,
    }
    return {
        "log_Dbar_order": [int(value) for value in n],
        "log_Dbar_covariance": variance.tolist(),
        "fixed_5_6": json_safe(fixed),
        "free_power": json_safe(free),
        "minimal_log_curvature": json_safe(curved),
        "nested_diagnostics": nested,
        "doubling_contrasts": doubling,
        "N680_heldout": heldout,
    }


def json_safe(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def read_external_lineage() -> dict:
    payload = subprocess.run(
        ["git", "show", f"{EXTERNAL_COMMIT}:{EXTERNAL_PATH}"],
        check=True, stdout=subprocess.PIPE,
    ).stdout
    if hashlib.sha256(payload).hexdigest() != EXTERNAL_SHA256:
        raise ValueError("external N325/N425 score object changed")
    source = json.loads(payload)
    output = {}
    for name in ("N325", "N425"):
        row = source["sizes"][name]
        first = row["orientations"]["first"]["D_N"]
        second = row["orientations"]["second"]["D_N"]
        covariance = np.asarray(row["delete_one_covariance"])
        indices = (1, 4)
        subcovariance = covariance[np.ix_(indices, indices)]
        weight = np.asarray([0.5, 0.5])
        variance = float(weight @ subcovariance @ weight)
        output[name] = {
            "N": row["N"],
            "D_first": first,
            "D_second": second,
            "Dbar": (first + second) / 2,
            "Dbar_standard_error": math.sqrt(variance),
        }
    return {
        "source_commit": EXTERNAL_COMMIT,
        "source_path": EXTERNAL_PATH,
        "source_sha256": EXTERNAL_SHA256,
        "sizes": output,
        "role": "external geometry-lineage comparison only; excluded from all P337 fits",
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P337 direct-birth four-generation scaling score", "",
        "| N | D_first | D_second | Dbar +/- SE |", "|---:|---:|---:|---:|",
    ]
    for name in result["size_order"]:
        row = result["sizes"][name]
        lines.append(
            f"| {row['N']} | {row['D_by_orientation']['first']:.8g} | "
            f"{row['D_by_orientation']['second']:.8g} | "
            f"{row['Dbar']:.8g} +/- {row['Dbar_standard_error']:.3g} |"
        )
    fixed = result["models"]["fixed_5_6"]
    free = result["models"]["free_power"]
    curved = result["models"]["minimal_log_curvature"]
    heldout = result["models"]["N680_heldout"]
    lines += [
        "",
        f"Fixed beta=5/6: A={fixed['amplitude']:.8g}, chi2={fixed['chi_square']:.6g}/3, p={fixed['survival_p']:.6g}, decision `{fixed['decision']}`.",
        f"Free power: beta={free['beta']:.8g} +/- {free['beta_standard_error']:.3g}, chi2={free['chi_square']:.6g}/2, p={free['survival_p']:.6g}.",
        f"Minimal curvature: kappa={curved['kappa']:.8g} +/- {curved['kappa_standard_error']:.3g}, p={curved['kappa_two_sided_p']:.6g}.",
        f"N680 heldout-by-analysis: observed={heldout['observed_Dbar']:.8g}, predicted={heldout['predicted_Dbar']:.8g}, z={heldout['z']:.4g}, p={heldout['two_sided_p']:.6g}.",
        "",
        f"Decision: `{result['decision']}`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--batch-output", type=Path, required=True)
    args = parser.parse_args()
    for path in (args.output, args.markdown, args.batch_output):
        if path.exists():
            raise ValueError(f"refusing to overwrite one-shot artifact: {path}")
    manifest = json.loads(args.manifest.read_text())
    size_rows = []
    batch_payload = []
    for name, spec in manifest["inputs"].items():
        births = Path(spec["births"])
        metadata = Path(spec["metadata"])
        if sha256(births) != spec["births_sha256"] or sha256(metadata) != spec["metadata_sha256"]:
            raise ValueError(f"{name}: input hash changed")
        if spec.get("compression") == "gzip":
            if uncompressed_sha256(births) != spec["births_uncompressed_sha256"]:
                raise ValueError(f"{name}: uncompressed input hash changed")
        metadata_payload = json.loads(metadata.read_text())
        if (
            int(metadata_payload["samples_per_pair"]) != int(spec["samples_per_orientation"])
            or int(metadata_payload["batches"]) != int(spec["batches"])
            or int(metadata_payload["seed"]) != int(spec["seed"])
        ):
            raise ValueError(f"{name}: metadata sampling contract changed")
        reduced = extract_batches(births, spec)
        summarized = summarize_size(reduced, spec)
        size_rows.append(summarized)
        direct = np.asarray(reduced["direct"])
        samples = np.asarray(reduced["samples"])
        for batch in range(int(spec["batches"])):
            batch_payload.append({
                "size": name,
                "N": int(spec["N"]),
                "batch": batch,
                "first_direct_count": int(direct[batch, 0]),
                "first_samples": int(samples[batch, 0]),
                "first_D": float(direct[batch, 0] / samples[batch, 0]),
                "second_direct_count": int(direct[batch, 1]),
                "second_samples": int(samples[batch, 1]),
                "second_D": float(direct[batch, 1] / samples[batch, 1]),
            })
    size_rows.sort(key=lambda row: row["N"])
    alpha = float(manifest["decision_alpha"])
    models = model_scores(size_rows, alpha)
    fixed = models["fixed_5_6"]
    result = {
        "schema": "matching-one/p337-direct-birth-six-arm-scaling-score/v1",
        "status": "existing_same_lineage_archives_scored_once",
        "freeze_commit": "9e73388",
        "scorer_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, text=True, stdout=subprocess.PIPE,
        ).stdout.strip(),
        "size_order": [f"N{row['N']}" for row in size_rows],
        "sizes": {f"N{row['N']}": row for row in size_rows},
        "models": models,
        "external_N325_N425": read_external_lineage(),
        "decision_alpha": alpha,
        "decision": (
            "conditional_fixed_5_6_line_survives"
            if fixed["survival_p"] >= alpha
            else "conditional_fixed_5_6_line_rejected"
        ),
        "claim_boundary": list(manifest["claim_boundary"]),
    }
    args.batch_output.parent.mkdir(parents=True, exist_ok=True)
    with args.batch_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(batch_payload[0]))
        writer.writeheader()
        writer.writerows(batch_payload)
    result["batch_sufficient_statistics"] = {
        "path": args.batch_output.as_posix(),
        "sha256": sha256(args.batch_output),
        "rows": len(batch_payload),
        "scope": "paired-orientation per-batch direct counts and denominators; sufficient to reconstruct every P337 covariance",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
