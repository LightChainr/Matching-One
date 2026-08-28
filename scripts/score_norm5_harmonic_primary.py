#!/usr/bin/env python3
"""Score the frozen N=325/425 norm-5 H4/H12 experiment.

The primary statistics are the two signed residuals

    DeltaM_child - r_model * DeltaM_parent

at the preregistered fixed probability.  Covariance is reconstructed from
delete-one batches. Runs with the same seed and exact counter interval form a
common-random-number group; disjoint intervals are independent even when their
batch labels happen to agree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Sequence

import mpmath as mp
import yaml

from analyze_matching_parity_derivatives_fast import H, combine, obs, read, remove


SIZES = (65, 85, 325, 425)
PARENTS = (65, 85)
CHILDREN = (325, 425)
LINKS = ((65, 325), (85, 425))
MODEL_ORDER = ("H4", "H12", "H8", "zero_effect")
EXPECTED_REPRESENTATIONS = {
    65: ((8, 1), (7, 4)),
    85: ((9, 2), (7, 6)),
    325: ((17, 6), (18, 1)),
    425: ((16, 13), (19, 8)),
}


@dataclass(frozen=True)
class Run:
    n: int
    histogram_path: Path
    metadata_path: Path
    metadata: Mapping[str, object]
    data: Mapping[tuple[int, str, int], H]

    @property
    def group_key(self) -> tuple[int, int, int]:
        return (
            int(self.metadata["seed"]),
            int(self.metadata["replica_counter_first"]),
            int(self.metadata["replica_counter_last_exclusive"]),
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_artifact(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("prediction artifact must be a mapping")
    if payload.get("status") != "preregistered_before_child_runs":
        raise ValueError("prediction artifact is not preregistered_before_child_runs")
    lineages = payload.get("lineages")
    if not isinstance(lineages, list) or len(lineages) != 2:
        raise ValueError("prediction artifact must contain two frozen lineages")
    observed_links = tuple(
        (int(row["parent"]["N"]), int(row["child_canonical_lineage"]["N"]))
        for row in lineages
    )
    if observed_links != LINKS:
        raise ValueError("frozen lineage order changed")
    order = tuple(payload.get("discrimination", {}).get("pure_model_scoring_order", ()))
    if order != MODEL_ORDER:
        raise ValueError("frozen pure-model scoring order changed")
    p_ref = mp.mpf(str(payload["common"]["fixed_probability"]))
    predictions = payload["exact_harmonic_predictions"]
    ratios: dict[str, mp.mpf] = {}
    for name in MODEL_ORDER[:-1]:
        angular = predictions[name]["angular_ratio"]
        ratios[name] = (
            mp.mpf(int(angular["numerator"]))
            / int(angular["denominator"])
            * mp.power(5, -mp.mpf(13) / 8)
        )
    ratios["zero_effect"] = mp.mpf(0)
    return {"payload": payload, "p_ref": p_ref, "ratios": ratios}


def load_run(n: int, histogram_path: Path, metadata_path: Path) -> Run:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not isinstance(metadata, dict):
        raise ValueError(f"N={n}: metadata must be an object")
    data = read(histogram_path)
    if {key[0] for key in data} != {n}:
        raise ValueError(f"N={n}: histogram contains a different size")
    for field in (
        "git_commit", "seed", "replica_counter_first",
        "replica_counter_last_exclusive", "samples_per_pair", "batches",
    ):
        if field not in metadata:
            raise ValueError(f"N={n}: metadata lacks {field}")
    commit = str(metadata["git_commit"])
    if len(commit) != 40 or any(ch not in "0123456789abcdef" for ch in commit.lower()):
        raise ValueError(f"N={n}: metadata requires a full hexadecimal commit")
    batches = int(metadata["batches"])
    if batches < 100:
        raise ValueError(f"N={n}: primary score requires at least 100 batches")
    batch_ids = sorted({key[2] for key in data})
    if batch_ids != list(range(batches)):
        raise ValueError(f"N={n}: incomplete zero-based batch grid")
    for orientation, expected in zip(("first", "second"), EXPECTED_REPRESENTATIONS[n]):
        rows = [row for key, row in data.items() if key[1] == orientation]
        if len(rows) != batches or {(row.a, row.b) for row in rows} != {expected}:
            raise ValueError(f"N={n}: frozen {orientation} representation changed")
        if sum(row.samples for row in rows) != int(metadata["samples_per_pair"]):
            raise ValueError(f"N={n}: histogram samples disagree with metadata")
    first = int(metadata["replica_counter_first"])
    last = int(metadata["replica_counter_last_exclusive"])
    if last - first != int(metadata["samples_per_pair"]):
        raise ValueError(f"N={n}: counter interval length disagrees with samples")
    return Run(n, histogram_path, metadata_path, metadata, data)


def validate_groups(runs: Sequence[Run]) -> dict[tuple[int, int, int], list[Run]]:
    if tuple(run.n for run in runs) != SIZES:
        raise ValueError(f"runs must be supplied in frozen size order {SIZES}")
    groups: dict[tuple[int, int, int], list[Run]] = {}
    for run in runs:
        groups.setdefault(run.group_key, []).append(run)
    for key, members in groups.items():
        batch_counts = {int(run.metadata["batches"]) for run in members}
        sample_counts = {int(run.metadata["samples_per_pair"]) for run in members}
        if len(batch_counts) != 1 or len(sample_counts) != 1:
            raise ValueError(f"CRN group {key} does not have aligned batches")
    for i, left in enumerate(runs):
        for right in runs[i + 1 :]:
            if int(left.metadata["seed"]) != int(right.metadata["seed"]):
                continue
            l0, l1 = left.group_key[1:]
            r0, r1 = right.group_key[1:]
            overlaps = max(l0, r0) < min(l1, r1)
            if overlaps and left.group_key != right.group_key:
                raise ValueError(
                    f"N={left.n} and N={right.n} have partial counter overlap; "
                    "the frozen scorer cannot infer that covariance"
                )
    return groups


def estimate_run(run: Run, p_ref: mp.mpf) -> tuple[mp.mpf, list[mp.mpf]]:
    by_orientation: dict[str, list[H]] = {}
    for orientation in ("first", "second"):
        by_orientation[orientation] = [
            run.data[key] for key in sorted(run.data)
            if key[1] == orientation
        ]
    first_total = combine(by_orientation["first"])
    second_total = combine(by_orientation["second"])

    def delta(first: H, second: H) -> mp.mpf:
        return obs(first, p_ref)["M"] - obs(second, p_ref)["M"]

    point = delta(first_total, second_total)
    deleted = [
        delta(remove(first_total, first), remove(second_total, second))
        for first, second in zip(by_orientation["first"], by_orientation["second"])
    ]
    return point, deleted


def jackknife_covariance(left: Sequence[mp.mpf], right: Sequence[mp.mpf]) -> mp.mpf:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("jackknife vectors are not aligned")
    count = len(left)
    mean_left = mp.fsum(left) / count
    mean_right = mp.fsum(right) / count
    return mp.mpf(count - 1) / count * mp.fsum(
        (x - mean_left) * (y - mean_right) for x, y in zip(left, right)
    )


def observation_covariance(
    runs: Sequence[Run], deleted: Mapping[int, Sequence[mp.mpf]]
) -> list[list[mp.mpf]]:
    groups = validate_groups(runs)
    index = {n: i for i, n in enumerate(SIZES)}
    covariance = [[mp.mpf(0) for _ in SIZES] for _ in SIZES]
    for members in groups.values():
        for left in members:
            for right in members:
                covariance[index[left.n]][index[right.n]] = jackknife_covariance(
                    deleted[left.n], deleted[right.n]
                )
    return covariance


def residual_covariance(
    covariance: Sequence[Sequence[mp.mpf]], ratio: mp.mpf
) -> list[list[mp.mpf]]:
    # Rows map [P65,P85,C325,C425] to [C325-rP65,C425-rP85].
    transform = [(-ratio, 0, 1, 0), (0, -ratio, 0, 1)]
    return [[
        mp.fsum(
            transform[i][a] * covariance[a][b] * transform[j][b]
            for a in range(4) for b in range(4)
        )
        for j in range(2)] for i in range(2)]


def quadratic2(vector: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]]) -> mp.mpf:
    a, b = covariance[0]
    c, d = covariance[1]
    determinant = a * d - b * c
    if determinant <= 0:
        raise ValueError("residual covariance is not positive definite")
    x, y = vector
    return (d * x * x - (b + c) * x * y + a * y * y) / determinant


def score_from_summary(
    point: Mapping[int, mp.mpf], covariance: Sequence[Sequence[mp.mpf]],
    ratios: Mapping[str, mp.mpf]
) -> list[dict[str, object]]:
    rows = []
    for name in MODEL_ORDER:
        ratio = ratios[name]
        residual = [point[child] - ratio * point[parent] for parent, child in LINKS]
        residual_cov = residual_covariance(covariance, ratio)
        diagonal = [[residual_cov[0][0], mp.mpf(0)], [mp.mpf(0), residual_cov[1][1]]]
        chi2 = quadratic2(residual, residual_cov)
        diagonal_chi2 = quadratic2(residual, diagonal)
        rows.append({
            "name": name,
            "child_over_parent_ratio": mp.nstr(ratio, 40),
            "residual": [mp.nstr(value, 30) for value in residual],
            "residual_covariance": [
                [mp.nstr(value, 20) for value in row] for row in residual_cov
            ],
            "marginal_signed_z": [
                mp.nstr(residual[i] / mp.sqrt(residual_cov[i][i]), 15)
                for i in range(2)
            ],
            "chi_square": mp.nstr(chi2, 20),
            "degrees_of_freedom": 2,
            "chi_square_survival_df2": mp.nstr(mp.exp(-chi2 / 2), 15),
            "diagonal_sensitivity_chi_square": mp.nstr(diagonal_chi2, 20),
        })
    best = min(mp.mpf(row["chi_square"]) for row in rows)
    for row in rows:
        row["delta_chi_square_from_best"] = mp.nstr(mp.mpf(row["chi_square"]) - best, 20)
    return rows


def render(
    runs: Sequence[Run], artifact_path: Path, artifact: Mapping[str, object]
) -> dict[str, object]:
    p_ref = artifact["p_ref"]
    estimated = {run.n: estimate_run(run, p_ref) for run in runs}
    points = {n: value[0] for n, value in estimated.items()}
    deleted = {n: value[1] for n, value in estimated.items()}
    covariance = observation_covariance(runs, deleted)
    groups = validate_groups(runs)
    return {
        "schema": "norm5 frozen harmonic primary score v1",
        "status": "prospective frozen-model score; no target refit",
        "fixed_probability": mp.nstr(p_ref, 30),
        "size_order": list(SIZES),
        "lineages": [list(link) for link in LINKS],
        "observation_delta_M": {str(n): mp.nstr(points[n], 30) for n in SIZES},
        "observation_covariance": [
            [mp.nstr(value, 20) for value in row] for row in covariance
        ],
        "counter_groups": [
            {
                "seed": key[0], "first": key[1], "last_exclusive": key[2],
                "sizes": [run.n for run in members],
                "covariance_rule": "aligned_delete_one" if len(members) > 1 else "independent_group",
            }
            for key, members in groups.items()
        ],
        "model_order": list(MODEL_ORDER),
        "models": score_from_summary(points, covariance, artifact["ratios"]),
        "provenance": {
            "prediction_artifact": str(artifact_path),
            "prediction_artifact_sha256": sha256(artifact_path),
            "inputs": [
                {
                    "N": run.n,
                    "histogram": str(run.histogram_path),
                    "histogram_sha256": sha256(run.histogram_path),
                    "metadata": str(run.metadata_path),
                    "metadata_sha256": sha256(run.metadata_path),
                    "git_commit": run.metadata["git_commit"],
                }
                for run in runs
            ],
        },
        "interpretation_guard": (
            "Report all four fixed models in the frozen order. Do not fit an exponent, "
            "change a lineage sign, omit a target, or treat matching batch labels from "
            "disjoint counter intervals as common randomness."
        ),
    }


def parse_run(specification: str) -> tuple[int, Path, Path]:
    fields = specification.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("run must be N:HISTOGRAM:METADATA")
    try:
        n = int(fields[0])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("run N must be an integer") from exc
    return n, Path(fields[1]), Path(fields[2])


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", action="append", type=parse_run, required=True,
        help="repeat in frozen order as N:HISTOGRAM:METADATA",
    )
    parser.add_argument(
        "--artifact", type=Path,
        default=root / "predictions/gaussian_norm5_harmonic_discrimination_20260828.yaml",
    )
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    artifact = load_artifact(args.artifact)
    runs = [load_run(n, hist, meta) for n, hist, meta in args.run]
    payload = render(runs, args.artifact, artifact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
