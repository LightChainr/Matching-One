#!/usr/bin/env python3
"""Score the frozen P334 production birth-age and collision diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, t


SCHEMA = "matching-one/p334-birth-age-production-score/v1"
ORIENTATIONS = ("first", "second")
HEADER = (
    "n,a,b,orientation,batch,samples,k1,k2,direct_0_to_2,site01,site12,"
    "line_null,ell_u,ell_v,iota01,iota12,physical_x,physical_y,chi4_re,chi4_im,"
    "mark01_valid,mark01_axis,mark01_diagonal,mark01_landed,mark01_h4,"
    "mark12_valid,mark12_axis,mark12_diagonal,mark12_landed,mark12_h4,count"
)
RAW_ROOT = Path(
    "/Volumes/Mac Data/Research/\u8bba\u6587\u9879\u76ee\u603b\u5e93/Matching-One-large-artifacts/"
    "P267-two-observer-source-rank-2m"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class BatchSummary:
    samples: int = 0
    direct: int = 0
    primary: dict[tuple[int, ...], list[int]] = field(default_factory=dict)
    complement: dict[tuple[int, ...], list[int]] = field(default_factory=dict)
    nested: dict[tuple[int, ...], list[int]] = field(default_factory=dict)
    complement_nested: dict[tuple[int, ...], list[int]] = field(default_factory=dict)


def add_slope_sufficient(
    table: dict[tuple[int, ...], list[int]],
    stratum: tuple[int, ...],
    age_steps: int,
    count: int,
    event: bool,
) -> None:
    # n, sum(n*a), sum(n*a^2), events, sum(events*a)
    row = table.setdefault(stratum, [0, 0, 0, 0, 0])
    row[0] += count
    row[1] += count * age_steps
    row[2] += count * age_steps * age_steps
    if event:
        row[3] += count
        row[4] += count * age_steps


def parse_archive(path: Path, spec: Mapping[str, object]) -> dict[str, list[BatchSummary]]:
    n_expected = int(spec["N"])
    k0 = int(spec["k0"])
    k0c = int(spec["complement_k0"])
    batches = int(spec["batches"])
    per_batch = int(spec["samples_per_batch"])
    output = {orientation: [BatchSummary() for _ in range(batches)] for orientation in ORIENTATIONS}
    rows_seen = 1
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n")
        if header != HEADER:
            raise ValueError(f"unexpected marked-birth header in {path}")
        for line_number, line in enumerate(handle, start=2):
            rows_seen += 1
            fields = line.rstrip("\n").split(",")
            if len(fields) != 31:
                raise ValueError(f"{path}:{line_number}: expected 31 columns")
            n = int(fields[0])
            orientation = fields[3]
            batch = int(fields[4])
            samples = int(fields[5])
            k1, k2 = int(fields[6]), int(fields[7])
            direct = int(fields[8])
            line_null = int(fields[11])
            ell = (int(fields[12]), int(fields[13]))
            mark01_valid = int(fields[20])
            mark01 = tuple(int(fields[index]) for index in (21, 22, 23, 24))
            mark12_valid = int(fields[25])
            mark12 = tuple(int(fields[index]) for index in (26, 27, 28, 29))
            count = int(fields[30])
            if n != n_expected or orientation not in output or not 0 <= batch < batches:
                raise ValueError(f"{path}:{line_number}: size/orientation/batch contract changed")
            if samples != per_batch or count <= 0 or not 1 <= k1 <= k2 <= n:
                raise ValueError(f"{path}:{line_number}: sample or birth support changed")
            if bool(direct) != (k1 == k2) or bool(line_null) != (k1 == k2):
                raise ValueError(f"{path}:{line_number}: direct/null-line contract changed")

            row = output[orientation][batch]
            row.samples += count
            if direct:
                row.direct += count

            k1c, k2c = n + 1 - k2, n + 1 - k1
            risk = k1 <= k0 < k2
            risk_c = k1c <= k0c < k2c
            if risk != risk_c:
                raise ValueError(f"{path}:{line_number}: complement risk-set identity failed")
            if risk:
                if line_null or ell == (0, 0) or not mark01_valid or not mark12_valid:
                    raise ValueError(f"{path}:{line_number}: rank-one survivor mark contract changed")
                add_slope_sufficient(row.primary, ell, k0 - k1, count, k2 == k0 + 1)
                add_slope_sufficient(
                    row.nested, ell + mark01, k0 - k1, count, k2 == k0 + 1
                )
                add_slope_sufficient(
                    row.complement, ell, k0c - k1c, count, k2c == k0c + 1
                )
                add_slope_sufficient(
                    row.complement_nested,
                    ell + mark12,
                    k0c - k1c,
                    count,
                    k2c == k0c + 1,
                )
    if rows_seen != int(spec["lines_including_header"]):
        raise ValueError(f"{path}: line count changed")
    for orientation in ORIENTATIONS:
        for batch, row in enumerate(output[orientation]):
            if row.samples != per_batch:
                raise ValueError(f"{path}: {orientation} batch {batch} sums to {row.samples}")
    return output


def merge_tables(rows: Sequence[BatchSummary], attribute: str) -> dict[tuple[int, ...], list[int]]:
    full: dict[tuple[int, ...], list[int]] = {}
    for row in rows:
        for key, values in getattr(row, attribute).items():
            target = full.setdefault(key, [0, 0, 0, 0, 0])
            for index, value in enumerate(values):
                target[index] += value
    return full


def age_slope(
    full: Mapping[tuple[int, ...], Sequence[int]],
    n: int,
    omitted: Mapping[tuple[int, ...], Sequence[int]] | None = None,
) -> dict[str, float | int]:
    numerator = 0.0
    denominator = 0.0
    survivors = 0
    events = 0
    contributing = 0
    for key, values in full.items():
        removed = omitted.get(key, (0, 0, 0, 0, 0)) if omitted else (0, 0, 0, 0, 0)
        count, sum_a, sum_a2, exits, sum_ay = (
            values[index] - removed[index] for index in range(5)
        )
        if count <= 0:
            continue
        survivors += count
        events += exits
        line_denominator = sum_a2 - sum_a * sum_a / count
        if line_denominator > 0.0:
            numerator += sum_ay - sum_a * exits / count
            denominator += line_denominator
            contributing += 1
    if denominator <= 0.0 or survivors <= 0:
        raise ValueError("birth-age slope is not identifiable")
    return {
        "beta_age_per_density": n * numerator / denominator,
        "next_step_hazard": events / survivors,
        "survivors": survivors,
        "next_step_exits": events,
        "strata": len(full),
        "contributing_strata": contributing,
    }


def jackknife_covariance(replicates: Sequence[Sequence[float]]) -> np.ndarray:
    values = np.asarray(replicates, dtype=float)
    centered = values - values.mean(axis=0)
    return (len(values) - 1) / len(values) * centered.T @ centered


def check_complement_audit(path: Path, expected_sha256: str) -> dict[str, object]:
    if sha256(path) != expected_sha256:
        raise ValueError(f"complement audit hash changed: {path}")
    totals = [0] * 6
    rows = 0
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split(",")
        expected = [
            "endpoint_failures", "site_failures", "line_failures",
            "local_mark_failures", "index_mismatches", "separated_mark_failures",
        ]
        indices = [header.index(name) for name in expected]
        for line in handle:
            fields = line.rstrip("\n").split(",")
            rows += 1
            for position, index in enumerate(indices):
                totals[position] += int(fields[index])
    if any(totals):
        raise ValueError(f"nonzero complement audit: {totals}")
    return {"rows": rows, "failure_sums": dict(zip(expected, totals))}


def score_size(
    archive: dict[str, list[BatchSummary]], spec: Mapping[str, object]
) -> dict[str, object]:
    n = int(spec["N"])
    batches = int(spec["batches"])
    point: list[float] = []
    point_details: dict[str, dict[str, object]] = {}
    full_tables: dict[tuple[str, str], dict[tuple[int, ...], list[int]]] = {}
    for orientation in ORIENTATIONS:
        rows = archive[orientation]
        for attribute in ("primary", "complement", "nested", "complement_nested"):
            full_tables[orientation, attribute] = merge_tables(rows, attribute)
        primary = age_slope(full_tables[orientation, "primary"], n)
        complement = age_slope(full_tables[orientation, "complement"], n)
        nested = age_slope(full_tables[orientation, "nested"], n)
        complement_nested = age_slope(full_tables[orientation, "complement_nested"], n)
        direct = sum(row.direct for row in rows)
        samples = sum(row.samples for row in rows)
        collision = direct / samples
        point.extend((primary["beta_age_per_density"], collision,
                      complement["beta_age_per_density"]))
        point_details[orientation] = {
            "primary": primary,
            "collision_count": direct,
            "samples": samples,
            "D_N": collision,
            "complement": complement,
            "nested_birth_mark_control": nested,
            "complement_nested_birth_mark_control": complement_nested,
        }

    # Append diagnostic coordinates after the six frozen primary coordinates.
    point.extend((
        point_details["first"]["nested_birth_mark_control"]["beta_age_per_density"],
        point_details["first"]["complement_nested_birth_mark_control"]["beta_age_per_density"],
        point_details["second"]["nested_birth_mark_control"]["beta_age_per_density"],
        point_details["second"]["complement_nested_birth_mark_control"]["beta_age_per_density"],
    ))

    deleted: list[list[float]] = []
    for omitted_batch in range(batches):
        vector: list[float] = []
        nested_vector: list[float] = []
        for orientation in ORIENTATIONS:
            omitted = archive[orientation][omitted_batch]
            primary = age_slope(full_tables[orientation, "primary"], n, omitted.primary)
            complement = age_slope(
                full_tables[orientation, "complement"], n, omitted.complement
            )
            nested = age_slope(full_tables[orientation, "nested"], n, omitted.nested)
            complement_nested = age_slope(
                full_tables[orientation, "complement_nested"],
                n,
                omitted.complement_nested,
            )
            direct = sum(row.direct for row in archive[orientation]) - omitted.direct
            samples = sum(row.samples for row in archive[orientation]) - omitted.samples
            vector.extend((primary["beta_age_per_density"], direct / samples,
                           complement["beta_age_per_density"]))
            nested_vector.extend((nested["beta_age_per_density"],
                                  complement_nested["beta_age_per_density"]))
        deleted.append(vector + nested_vector)
    covariance = jackknife_covariance(deleted)
    df = batches - 1
    critical = float(t.ppf(0.995, df))

    for orientation_index, orientation in enumerate(ORIENTATIONS):
        beta_index = 3 * orientation_index
        d_index = beta_index + 1
        complement_index = beta_index + 2
        nested_index = 6 + 2 * orientation_index
        complement_nested_index = nested_index + 1
        row = point_details[orientation]
        for label, index in (("primary", beta_index), ("complement", complement_index)):
            estimate = point[index]
            se = math.sqrt(covariance[index, index])
            statistic = estimate / se
            row[label].update({
                "delete_one_se": se,
                "student_t": statistic,
                "two_sided_p_df99": float(2 * t.sf(abs(statistic), df)),
                "resolved_at_alpha_0_01": bool(2 * t.sf(abs(statistic), df) < 0.01),
            })
        d_se = math.sqrt(covariance[d_index, d_index])
        row["D_N_delete_one_se"] = d_se
        row["D_N_99pct_interval"] = [point[d_index] - critical * d_se,
                                      point[d_index] + critical * d_se]
        row["beta_D_correlation"] = float(
            covariance[beta_index, d_index]
            / math.sqrt(covariance[beta_index, beta_index] * covariance[d_index, d_index])
        )
        for label, nested_i, primary_i in (
            ("nested_birth_mark_control", nested_index, beta_index),
            ("complement_nested_birth_mark_control", complement_nested_index, complement_index),
        ):
            nested_row = row[label]
            nested_estimate = point[nested_i]
            nested_se = math.sqrt(covariance[nested_i, nested_i])
            difference = nested_estimate - point[primary_i]
            difference_variance = (
                covariance[nested_i, nested_i] + covariance[primary_i, primary_i]
                - 2 * covariance[nested_i, primary_i]
            )
            nested_row.update({
                "delete_one_se": nested_se,
                "two_sided_p_df99": float(2 * t.sf(abs(nested_estimate / nested_se), df)),
                "difference_from_primary": difference,
                "difference_delete_one_se": math.sqrt(max(difference_variance, 0.0)),
            })

    primary_indices = [0, 3]
    primary_vector = np.asarray([point[index] for index in primary_indices])
    primary_cov = covariance[np.ix_(primary_indices, primary_indices)]
    joint_statistic = float(primary_vector @ np.linalg.pinv(primary_cov) @ primary_vector)
    complement_indices = [2, 5]
    complement_vector = np.asarray([point[index] for index in complement_indices])
    complement_cov = covariance[np.ix_(complement_indices, complement_indices)]
    complement_statistic = float(
        complement_vector @ np.linalg.pinv(complement_cov) @ complement_vector
    )
    dbar = 0.5 * (point[1] + point[4])
    dbar_vector = np.zeros(len(point))
    dbar_vector[1] = dbar_vector[4] = 0.5
    dbar_variance = float(dbar_vector @ covariance @ dbar_vector)
    return {
        "N": n,
        "k0": int(spec["k0"]),
        "complement_k0": int(spec["complement_k0"]),
        "orientations": point_details,
        "point_vector_order": [
            "first_beta_age", "first_D_N", "first_complement_beta_age",
            "second_beta_age", "second_D_N", "second_complement_beta_age",
            "first_nested_beta_age", "first_complement_nested_beta_age",
            "second_nested_beta_age", "second_complement_nested_beta_age",
        ],
        "point_vector": point,
        "delete_one_covariance": covariance.tolist(),
        "primary_two_orientation_joint": {
            "chi_square": joint_statistic,
            "df": 2,
            "survival_p": float(chi2.sf(joint_statistic, 2)),
        },
        "complement_two_orientation_joint": {
            "chi_square": complement_statistic,
            "df": 2,
            "survival_p": float(chi2.sf(complement_statistic, 2)),
        },
        "Dbar_orientation_mean": dbar,
        "Dbar_delete_one_se": math.sqrt(dbar_variance),
        "Dbar_variance": dbar_variance,
    }


def six_arm_score(sizes: Mapping[str, Mapping[str, object]], ratio: float) -> dict[str, object]:
    low, high = sizes["N325"], sizes["N425"]
    contrast = high["Dbar_orientation_mean"] - ratio * low["Dbar_orientation_mean"]
    variance = high["Dbar_variance"] + ratio * ratio * low["Dbar_variance"]
    statistic = contrast * contrast / variance
    return {
        "fixed_ratio": ratio,
        "observed_ratio": high["Dbar_orientation_mean"] / low["Dbar_orientation_mean"],
        "contrast_Dbar425_minus_ratio_Dbar325": contrast,
        "contrast_se": math.sqrt(variance),
        "chi_square": statistic,
        "df": 1,
        "survival_p": float(chi2.sf(statistic, 1)),
        "status": "conditional_high_risk_secondary",
    }


def score(
    manifest: Mapping[str, object], raw_paths: Mapping[str, Path], audit_paths: Mapping[str, Path]
) -> dict[str, object]:
    sizes: dict[str, object] = {}
    input_audits: dict[str, object] = {}
    for name in ("N325", "N425"):
        spec = manifest["inputs"][name]
        if sha256(raw_paths[name]) != spec["sha256"]:
            raise ValueError(f"{name} raw hash changed")
        complement = check_complement_audit(
            audit_paths[name], spec["complement_audit_sha256"]
        )
        archive = parse_archive(raw_paths[name], spec)
        sizes[name] = score_size(archive, spec)
        input_audits[name] = {
            "raw_path": str(raw_paths[name]),
            "raw_sha256": spec["sha256"],
            "complement_audit": complement,
        }
    ratio = float(manifest["six_arm_collision_adversary"]["fixed_ratio"])
    return {
        "schema": SCHEMA,
        "status": "retrospective_existing_production_reveal",
        "freeze_manifest_sha256": sha256(Path("analysis/p334_birth_age_production_freeze.json")),
        "dependency_groups": {
            "N325": "one 2M paired-orientation archive",
            "N425": "one independent 2M paired-orientation archive",
        },
        "input_audits": input_audits,
        "sizes": sizes,
        "six_arm_collision_adversary": six_arm_score(sizes, ratio),
        "decision_alpha": manifest["inference"]["alpha"],
        "claim_boundary": manifest["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=Path("analysis/p334_birth_age_production_freeze.json"))
    parser.add_argument("--n325", type=Path,
                        default=RAW_ROOT / "N325_2m/N325_2m.marked_births.csv")
    parser.add_argument("--n425", type=Path,
                        default=RAW_ROOT / "N425_2m/N425_2m.marked_births.csv")
    parser.add_argument("--n325-audit", type=Path,
                        default=RAW_ROOT / "N325_2m/N325_2m.complement_audit.csv")
    parser.add_argument("--n425-audit", type=Path,
                        default=RAW_ROOT / "N425_2m/N425_2m.complement_audit.csv")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(
        manifest,
        {"N325": args.n325, "N425": args.n425},
        {"N325": args.n325_audit, "N425": args.n425_audit},
    )
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
