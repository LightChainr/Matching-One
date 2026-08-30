#!/usr/bin/env python3
"""Score exact iota and cyclic birth-site controls for P334 birth age."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, t


ORIENTATIONS = ("first", "second")
CONTROL_NAMES = ("primary", "iota_pair", "smith_site_pair")
PARENT_SCORE = Path("results/local-20260830/P334-birth-age-production/score.json")
PARENT_SCORE_SHA256 = "4ab6ddf989b8cce5dad3365d5d440b387a931239c8b88a0d5c62f3e508c3ab29"
HEADER = (
    "n,a,b,orientation,batch,samples,k1,k2,direct_0_to_2,site01,site12,"
    "line_null,ell_u,ell_v,iota01,iota12,physical_x,physical_y,chi4_re,chi4_im,"
    "mark01_valid,mark01_axis,mark01_diagonal,mark01_landed,mark01_h4,"
    "mark12_valid,mark12_axis,mark12_diagonal,mark12_landed,mark12_h4,count"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class Batch:
    samples: int = 0
    tables: dict[str, dict[tuple[int, ...], list[int]]] = field(
        default_factory=lambda: {name: {} for name in CONTROL_NAMES}
    )


def add_sufficient(
    table: dict[tuple[int, ...], list[int]], stratum: tuple[int, ...],
    age: int, count: int, event: bool,
) -> None:
    # count, sum(count*age), sum(count*age^2), exits, sum(exits*age)
    target = table.setdefault(stratum, [0, 0, 0, 0, 0])
    target[0] += count
    target[1] += count * age
    target[2] += count * age * age
    if event:
        target[3] += count
        target[4] += count * age


def parse(path: Path, spec: Mapping[str, object]) -> dict:
    n = int(spec["N"])
    k0 = int(spec["k0"])
    batches = int(spec["batches"])
    per_batch = int(spec["samples_per_batch"])
    output = {orientation: [Batch() for _ in range(batches)] for orientation in ORIENTATIONS}
    iota_support: Counter[tuple[int, int]] = Counter()
    site_support: Counter[int] = Counter()
    rows_seen = 1
    with path.open(encoding="utf-8") as handle:
        if handle.readline().rstrip("\n") != HEADER:
            raise ValueError(f"{path}: marked-birth header changed")
        for line_number, line in enumerate(handle, start=2):
            rows_seen += 1
            fields = line.rstrip("\n").split(",")
            if len(fields) != 31:
                raise ValueError(f"{path}:{line_number}: wrong column count")
            if int(fields[0]) != n:
                raise ValueError(f"{path}:{line_number}: wrong N")
            orientation = fields[3]
            batch_index = int(fields[4])
            samples = int(fields[5])
            k1, k2 = int(fields[6]), int(fields[7])
            site01, site12 = int(fields[9]), int(fields[10])
            line_null = int(fields[11])
            ell = (int(fields[12]), int(fields[13]))
            iota = (int(fields[14]), int(fields[15]))
            count = int(fields[30])
            if (
                orientation not in output or not 0 <= batch_index < batches
                or samples != per_batch or count <= 0 or not 1 <= k1 <= k2 <= n
            ):
                raise ValueError(f"{path}:{line_number}: sampling/support contract changed")
            row = output[orientation][batch_index]
            row.samples += count
            if k1 <= k0 < k2:
                if line_null or ell == (0, 0) or not 0 <= site01 < n or not 0 <= site12 < n:
                    raise ValueError(f"{path}:{line_number}: rank-one/site contract changed")
                delta = (site12 - site01) % n
                site_class = math.gcd(delta, n)
                age = k0 - k1
                event = k2 == k0 + 1
                add_sufficient(row.tables["primary"], ell, age, count, event)
                add_sufficient(row.tables["iota_pair"], ell + iota, age, count, event)
                add_sufficient(
                    row.tables["smith_site_pair"], ell + iota + (site_class,), age, count, event,
                )
                iota_support[iota] += count
                site_support[site_class] += count
    if rows_seen != int(spec["lines_including_header"]):
        raise ValueError(f"{path}: line count changed")
    for orientation in ORIENTATIONS:
        for batch_index, row in enumerate(output[orientation]):
            if row.samples != per_batch:
                raise ValueError(f"{path}: {orientation} batch {batch_index} does not partition")
    return {
        "batches": output,
        "iota_support": iota_support,
        "site_support": site_support,
        "rows_seen": rows_seen,
    }


def merge(rows: Sequence[Batch], control: str) -> dict[tuple[int, ...], list[int]]:
    output: dict[tuple[int, ...], list[int]] = {}
    for row in rows:
        for key, values in row.tables[control].items():
            target = output.setdefault(key, [0, 0, 0, 0, 0])
            for index, value in enumerate(values):
                target[index] += value
    return output


def slope(
    full: Mapping[tuple[int, ...], Sequence[int]], n: int,
    omitted: Mapping[tuple[int, ...], Sequence[int]] | None = None,
) -> dict:
    numerator = 0.0
    denominator = 0.0
    survivors = 0
    identifying_survivors = 0
    exits = 0
    contributing = 0
    for key, values in full.items():
        removed = omitted.get(key, (0, 0, 0, 0, 0)) if omitted else (0, 0, 0, 0, 0)
        count, sum_age, sum_age2, events, sum_event_age = (
            values[index] - removed[index] for index in range(5)
        )
        if count <= 0:
            continue
        survivors += count
        exits += events
        local_denominator = sum_age2 - sum_age * sum_age / count
        if local_denominator > 0:
            numerator += sum_event_age - sum_age * events / count
            denominator += local_denominator
            identifying_survivors += count
            contributing += 1
    if denominator <= 0 or survivors <= 0:
        raise ValueError("controlled age slope is not identifiable")
    return {
        "beta_age_per_density": n * numerator / denominator,
        "within_stratum_age_denominator_steps2": denominator,
        "survivors": survivors,
        "identifying_survivors": identifying_survivors,
        "next_step_exits": exits,
        "next_step_hazard": exits / survivors,
        "strata": len(full),
        "contributing_strata": contributing,
    }


def jackknife_covariance(values: Sequence[Sequence[float]]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    centered = array - array.mean(axis=0)
    return (len(array) - 1) / len(array) * centered.T @ centered


def quadratic_score(point: np.ndarray, covariance: np.ndarray) -> dict:
    statistic = float(point @ np.linalg.pinv(covariance, rcond=1e-12) @ point)
    rank = int(np.linalg.matrix_rank(covariance, tol=1e-12 * np.max(np.diag(covariance))))
    return {
        "chi_square": statistic,
        "degrees_of_freedom": rank,
        "survival_p": float(chi2.sf(statistic, rank)) if rank else 1.0,
    }


def score_size(parsed: Mapping[str, object], spec: Mapping[str, object]) -> tuple[dict, list[dict]]:
    n = int(spec["N"])
    batch_count = int(spec["batches"])
    archive = parsed["batches"]
    full_tables = {
        (orientation, control): merge(archive[orientation], control)
        for orientation in ORIENTATIONS for control in CONTROL_NAMES
    }
    details = {
        orientation: {
            control: slope(full_tables[orientation, control], n)
            for control in CONTROL_NAMES
        }
        for orientation in ORIENTATIONS
    }
    point = np.asarray([
        details[orientation][control]["beta_age_per_density"]
        for orientation in ORIENTATIONS for control in CONTROL_NAMES
    ])
    deleted = []
    batch_rows = []
    for omitted_batch in range(batch_count):
        vector = []
        for orientation in ORIENTATIONS:
            omitted = archive[orientation][omitted_batch]
            for control in CONTROL_NAMES:
                vector.append(slope(
                    full_tables[orientation, control], n, omitted.tables[control],
                )["beta_age_per_density"])
        deleted.append(vector)
        batch_rows.append({
            "N": n,
            "batch_deleted": omitted_batch,
            **{
                f"{orientation}_{control}": vector[3 * oi + ci]
                for oi, orientation in enumerate(ORIENTATIONS)
                for ci, control in enumerate(CONTROL_NAMES)
            },
        })
    covariance = jackknife_covariance(deleted)
    df = batch_count - 1
    alpha = 0.01
    control_results = {}
    for ci, control in enumerate(CONTROL_NAMES):
        indices = [ci, 3 + ci]
        joint = quadratic_score(point[indices], covariance[np.ix_(indices, indices)])
        orientation_rows = {}
        denominator_retentions = []
        for oi, orientation in enumerate(ORIENTATIONS):
            index = 3 * oi + ci
            primary_index = 3 * oi
            estimate = float(point[index])
            standard_error = math.sqrt(float(covariance[index, index]))
            difference = estimate - float(point[primary_index])
            difference_variance = float(
                covariance[index, index] + covariance[primary_index, primary_index]
                - 2 * covariance[index, primary_index]
            )
            denominator_retention = (
                details[orientation][control]["within_stratum_age_denominator_steps2"]
                / details[orientation]["primary"]["within_stratum_age_denominator_steps2"]
            )
            denominator_retentions.append(denominator_retention)
            exact_noop = abs(difference_variance) < 1e-24 and abs(difference) < 1e-14
            difference_se = math.sqrt(max(difference_variance, 0.0))
            orientation_rows[orientation] = {
                **details[orientation][control],
                "delete_one_standard_error": standard_error,
                "student_t": estimate / standard_error,
                "two_sided_p_df99": float(2 * t.sf(abs(estimate / standard_error), df)),
                "difference_from_primary": difference,
                "difference_standard_error": difference_se,
                "difference_two_sided_p_df99": (
                    1.0 if exact_noop else float(2 * t.sf(abs(difference / difference_se), df))
                ),
                "exact_noop_vs_primary": exact_noop,
                "absolute_magnitude_retention": (
                    abs(estimate) / abs(float(point[primary_index]))
                ),
                "age_denominator_retention": denominator_retention,
            }
        difference_transform = np.zeros((2, 6))
        difference_transform[0, ci] = 1.0
        difference_transform[0, 0] -= 1.0
        difference_transform[1, 3 + ci] = 1.0
        difference_transform[1, 3] -= 1.0
        difference_point = difference_transform @ point
        difference_covariance = difference_transform @ covariance @ difference_transform.T
        difference_joint = quadratic_score(difference_point, difference_covariance)
        identifiable = min(denominator_retentions) >= 0.25
        control_results[control] = {
            "orientations": orientation_rows,
            "two_orientation_joint_slope": joint,
            "two_orientation_joint_difference_from_primary": difference_joint,
            "minimum_age_denominator_retention": min(denominator_retentions),
            "identifiability_gate_passed": identifiable,
            "association_survives": identifiable and joint["survival_p"] < alpha,
        }
    return {
        "N": n,
        "point_vector_order": [
            f"{orientation}_{control}"
            for orientation in ORIENTATIONS for control in CONTROL_NAMES
        ],
        "point_vector": point.tolist(),
        "delete_one_covariance": covariance.tolist(),
        "controls": control_results,
        "iota_pair_support": [
            {"iota01": key[0], "iota12": key[1], "count": count}
            for key, count in sorted(parsed["iota_support"].items())
        ],
        "smith_site_pair_support": [
            {"gcd_delta_site_N": key, "count": count}
            for key, count in sorted(parsed["site_support"].items())
        ],
        "iota_exactly_saturated": set(parsed["iota_support"]) == {(1, 1)},
        "input_audit": {"lines_including_header": int(parsed["rows_seen"])},
    }, batch_rows


def replay_gate(result: Mapping[str, object], parent: Mapping[str, object]) -> dict:
    maximum_point_difference = 0.0
    maximum_covariance_difference = 0.0
    for name in ("N325", "N425"):
        observed = result[name]
        expected = parent["sizes"][name]
        observed_indices = [0, 3]
        expected_indices = [0, 3]
        point_difference = np.max(np.abs(
            np.asarray(observed["point_vector"])[observed_indices]
            - np.asarray(expected["point_vector"])[expected_indices]
        ))
        observed_covariance = np.asarray(observed["delete_one_covariance"])[np.ix_(observed_indices, observed_indices)]
        expected_covariance = np.asarray(expected["delete_one_covariance"])[np.ix_(expected_indices, expected_indices)]
        covariance_difference = np.max(np.abs(observed_covariance - expected_covariance))
        maximum_point_difference = max(maximum_point_difference, float(point_difference))
        maximum_covariance_difference = max(maximum_covariance_difference, float(covariance_difference))
    passed = maximum_point_difference <= 1e-13 and maximum_covariance_difference <= 1e-13
    return {
        "passed": passed,
        "maximum_point_abs_difference": maximum_point_difference,
        "maximum_covariance_abs_difference": maximum_covariance_difference,
        "tolerances": {"point": 1e-13, "covariance": 1e-13},
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P334 exact iota and birth-site control score", "",
        "| N | control | first beta | second beta | joint p | min age-info retained | survives |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for name in ("N325", "N425"):
        row = result["sizes"][name]
        for control in CONTROL_NAMES:
            item = row["controls"][control]
            lines.append(
                f"| {row['N']} | {control} | "
                f"{item['orientations']['first']['beta_age_per_density']:.7g} | "
                f"{item['orientations']['second']['beta_age_per_density']:.7g} | "
                f"{item['two_orientation_joint_slope']['survival_p']:.6g} | "
                f"{item['minimum_age_denominator_retention']:.4g} | "
                f"{item['association_survives']} |"
            )
    lines += [
        "",
        f"Primary replay: `{result['primary_replay']['passed']}`.",
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
    for target in (args.output, args.markdown, args.batch_output):
        if target.exists():
            raise ValueError(f"refusing to overwrite one-shot artifact: {target}")
    manifest = json.loads(args.manifest.read_text())
    if sha256(PARENT_SCORE) != PARENT_SCORE_SHA256:
        raise ValueError("parent P334 score hash changed")
    parent = json.loads(PARENT_SCORE.read_text())
    sizes = {}
    batch_rows = []
    for name, spec in manifest["inputs"].items():
        path = Path(spec["path"])
        if sha256(path) != spec["sha256"]:
            raise ValueError(f"{name}: raw hash changed")
        size, rows = score_size(parse(path, spec), spec)
        sizes[name] = size
        batch_rows.extend(rows)
    replay = replay_gate(sizes, parent)
    if not replay["passed"]:
        raise ValueError(f"primary replay failed: {replay}")
    site_survives = all(
        sizes[name]["controls"]["smith_site_pair"]["association_survives"]
        for name in sizes
    )
    result = {
        "schema": "matching-one/p334-age-iota-site-control-score/v1",
        "status": "existing_archives_nested_control_scored_once",
        "freeze_commit": "cccc10f",
        "scorer_commit": __import__("subprocess").run(
            ["git", "rev-parse", "HEAD"], check=True, text=True,
            stdout=__import__("subprocess").PIPE,
        ).stdout.strip(),
        "primary_replay": replay,
        "sizes": sizes,
        "decision_alpha": float(manifest["decision_alpha"]),
        "decision": (
            "age_association_survives_exact_iota_and_site_pair_control"
            if site_survives else "age_association_not_resolved_after_exact_site_pair_control"
        ),
        "temporal_boundary": manifest["temporal_boundary"],
        "claim_boundary": manifest["claim_boundary"],
    }
    args.batch_output.parent.mkdir(parents=True, exist_ok=True)
    with args.batch_output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(batch_rows[0]))
        writer.writeheader()
        writer.writerows(batch_rows)
    result["batch_delete_one_slopes"] = {
        "path": args.batch_output.as_posix(),
        "sha256": sha256(args.batch_output),
        "rows": len(batch_rows),
        "contract": "one common batch deletion across both orientations within each size; sufficient to reconstruct both 6x6 covariance matrices",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
