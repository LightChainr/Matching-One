#!/usr/bin/env python3
"""Score the frozen N=145 -> 290 full-curve lineage.

The N145 parent and N290 child production blocks use independent random
streams. Each size is jackknifed internally; cross-size residual covariance is
therefore the sum of the two size-local covariance contributions. No matching
batch labels are interpreted as common randomness.

Primary order follows predictions/p49_slope_two_sector_145_290_20260828.yaml:
  1. intrinsic/thermal-even DeltaM transfer;
  2. raw 2^(3/8) slope baseline;
  3. frozen scalar+H4 finite-size slope correction;
  4. raw -1/4 root-ratio baseline;
  5. frozen induced root-ratio correction.

Derivative P4 channels are reported afterward as correlated diagnostics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import yaml
import mpmath as mp

from analyze_p48_retrospective import (
    Histogram,
    covariance_of_mean,
    project_size,
    read_histograms,
)
from covariance_nullspace import covariance_spectral_diagnostics, serialize_diagnostics
from score_p49_fullcurve_doubling import (
    LEVELS,
    P48_POWERS,
    aggregate,
    orientation_values,
    solve_target,
)

PARENT_N = 145
CHILD_N = 290
EXPECTED_REPRESENTATIONS = {
    PARENT_N: {"first": (12, 1), "second": (9, 8)},
    CHILD_N: {"first": (13, 11), "second": (17, 1)},
}
# Both engines emit first-minus-second in the frozen genealogy order.  The
# negative norm-2 transfer character is already carried by the target ratio.
LINEAGE_SIGN = {PARENT_N: +1.0, CHILD_N: +1.0}
FEATURE_ORDER = (
    "X_even_0.0",
    "X_even_0.025",
    "X_even_0.05",
    "mean_slope",
    "root_gap_lineage",
    "P4_S",
    "P4_D",
    "P4_S_prime",
    "P4_D_prime",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_metadata(path: Path) -> Mapping[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: metadata must be a JSON object")
    for key in (
        "seed",
        "replica_counter_first",
        "replica_counter_last_exclusive",
        "samples_per_pair",
        "batches",
        "git_commit",
    ):
        if key not in payload:
            raise ValueError(f"{path}: metadata lacks {key}")
    return payload


def rng_group(metadata: Mapping[str, object]) -> tuple[int, int, int]:
    return (
        int(metadata["seed"]),
        int(metadata["replica_counter_first"]),
        int(metadata["replica_counter_last_exclusive"]),
    )


def read_one_size(path: Path, n: int) -> Dict[Tuple[int, str, int], Histogram]:
    data = read_histograms(path)
    sizes = {key[0] for key in data}
    if sizes != {n}:
        raise ValueError(f"{path}: expected only N={n}, got {sorted(sizes)}")
    expected = EXPECTED_REPRESENTATIONS[n]
    for orientation in ("first", "second"):
        rows = [row for key, row in data.items() if key[1] == orientation]
        if not rows:
            raise ValueError(f"N={n}: missing {orientation} rows")
        reps = {(row.a, row.b) for row in rows}
        if reps != {expected[orientation]}:
            raise ValueError(
                f"N={n}: {orientation} representation {reps} != {expected[orientation]}"
            )
    return data


def grouped(data: Mapping[Tuple[int, str, int], Histogram], n: int):
    output = {}
    for orientation in ("first", "second"):
        rows = sorted(
            (row for key, row in data.items() if key[:2] == (n, orientation)),
            key=lambda row: row.batch,
        )
        ids = [row.batch for row in rows]
        if ids != list(range(len(rows))):
            raise ValueError(f"N={n}: {orientation} batches are not zero-based contiguous")
        output[orientation] = rows
    if len(output["first"]) != len(output["second"]):
        raise ValueError(f"N={n}: orientation batch counts differ")
    if len(output["first"]) < 2:
        raise ValueError(f"N={n}: at least two batches are required")
    return output


def size_statistics(by_orientation, *, lineage_sign: float, omitted: int = -1):
    n = by_orientation["first"][0].n
    rows = {
        name: aggregate(by_orientation[name], omitted)
        for name in ("first", "second")
    }

    def value(name: str, p: float):
        return orientation_values(n, rows[name], p)

    def mean_matching(p: float) -> float:
        return (value("first", p)["M"] + value("second", p)["M"]) / 2.0

    x_even = {}
    for level in LEVELS:
        if level == 0.0:
            p_minus = p_plus = solve_target(mean_matching, 0.0)
        else:
            p_minus = solve_target(mean_matching, -level)
            p_plus = solve_target(mean_matching, level)
        minus_first = value("first", p_minus)
        minus_second = value("second", p_minus)
        plus_first = value("first", p_plus)
        plus_second = value("second", p_plus)
        x_minus = lineage_sign * (minus_first["M"] - minus_second["M"])
        x_plus = lineage_sign * (plus_first["M"] - plus_second["M"])
        x_even[str(level)] = (x_plus + x_minus) / 2.0

    p0 = solve_target(mean_matching, 0.0)
    center = {name: value(name, p0) for name in ("first", "second")}
    mean_slope = (center["first"]["M_prime"] + center["second"]["M_prime"]) / 2.0
    roots = {
        name: solve_target(lambda p, name=name: value(name, p)["M"], 0.0)
        for name in ("first", "second")
    }
    root_gap = lineage_sign * (roots["first"] - roots["second"])
    p48 = project_size(by_orientation, omitted)
    return {
        "X_even_0.0": x_even["0.0"],
        "X_even_0.025": x_even["0.025"],
        "X_even_0.05": x_even["0.05"],
        "mean_slope": mean_slope,
        "root_gap_lineage": root_gap,
        "P4_S": p48["P4_S"],
        "P4_D": p48["P4_D"],
        "P4_S_prime": p48["P4_S_prime"],
        "P4_D_prime": p48["P4_D_prime"],
        "p0": p0,
    }


def vector(stat: Mapping[str, float]) -> list[float]:
    return [float(stat[name]) for name in FEATURE_ORDER]


def pseudovalue_vectors(
    full: Sequence[float], deleted: Sequence[Sequence[float]]
) -> list[list[float]]:
    batches = len(deleted)
    if any(len(row) != len(full) for row in deleted):
        raise ValueError("ragged delete-one vectors")
    return [
        [batches * full[j] - (batches - 1) * row[j] for j in range(len(full))]
        for row in deleted
    ]


def estimate(by_orientation, *, lineage_sign: float):
    point = size_statistics(by_orientation, lineage_sign=lineage_sign)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        vector(size_statistics(by_orientation, lineage_sign=lineage_sign, omitted=batch))
        for batch in batch_ids
    ]
    point_vector = vector(point)
    pseudo = pseudovalue_vectors(point_vector, deleted)
    covariance = covariance_of_mean(pseudo)
    return point, covariance


def submatrix(matrix: Sequence[Sequence[float]], indices: Sequence[int]) -> list[list[float]]:
    return [[float(matrix[i][j]) for j in indices] for i in indices]


def independent_residual_covariance(
    parent_cov: Sequence[Sequence[float]],
    child_cov: Sequence[Sequence[float]],
    ratio: float,
) -> list[list[float]]:
    if len(parent_cov) != len(child_cov):
        raise ValueError("parent/child covariance dimensions differ")
    return [
        [
            float(child_cov[i][j]) + ratio * ratio * float(parent_cov[i][j])
            for j in range(len(parent_cov))
        ]
        for i in range(len(parent_cov))
    ]


def scalar_score(
    parent: float,
    child: float,
    var_parent: float,
    var_child: float,
    ratio: float,
    *,
    ratio_se: float = 0.0,
) -> dict[str, float]:
    residual = child - ratio * parent
    variance = var_child + ratio * ratio * var_parent + parent * parent * ratio_se * ratio_se
    if variance <= 0.0:
        raise ValueError("nonpositive scalar residual variance")
    z = residual / math.sqrt(variance)
    return {
        "ratio": ratio,
        "ratio_source_se": ratio_se,
        "residual": residual,
        "variance": variance,
        "signed_z": z,
        "chi_square": z * z,
        "degrees_of_freedom": 1,
    }


def generalized_covariance_score(
    residual: Sequence[float],
    covariance: Sequence[Sequence[float]],
    relative_cutoff: float = 1e-10,
) -> dict[str, object]:
    """Score a correlated vector and audit every discarded covariance mode."""
    mp.mp.dps = 80
    values = [mp.mpf(str(value)) for value in residual]
    matrix = [[mp.mpf(str(value)) for value in row] for row in covariance]
    diagnostics = covariance_spectral_diagnostics(
        values,
        matrix,
        mp.mpf(str(relative_cutoff)),
        nullspace_policy="estimated",
    )
    return serialize_diagnostics(diagnostics, float)


def render(
    parent_data,
    child_data,
    parent_meta: Mapping[str, object],
    child_meta: Mapping[str, object],
    prediction: Mapping[str, object],
    provenance: Mapping[str, object],
) -> dict[str, object]:
    if rng_group(parent_meta) == rng_group(child_meta):
        raise ValueError(
            "N145/N290 scorer assumes independent streams; metadata declare one shared group"
        )

    parent_groups = grouped(parent_data, PARENT_N)
    child_groups = grouped(child_data, CHILD_N)
    parent, parent_cov = estimate(
        parent_groups, lineage_sign=LINEAGE_SIGN[PARENT_N]
    )
    child, child_cov = estimate(
        child_groups, lineage_sign=LINEAGE_SIGN[CHILD_N]
    )

    target = prediction["target"]
    induced = prediction["induced_root_prediction_if_primary_deltaM_doubling_holds"]
    delta_ratio = float(induced["deltaM_lineage_ratio"])
    slope_baseline = float(target["asymptotic_ratio_2pow3over8"])
    slope_corrected = float(target["frozen_slope_ratio"])
    slope_corrected_se = float(target["linearized_source_se_approx"])
    root_baseline = -0.25
    root_corrected = float(induced["frozen_root_gap_ratio"])
    root_corrected_se = float(induced["linearized_slope_source_se_only_approx"])

    delta_indices = [0, 1, 2]
    parent_delta = [float(vector(parent)[i]) for i in delta_indices]
    child_delta = [float(vector(child)[i]) for i in delta_indices]
    delta_residual = [
        child_delta[i] - delta_ratio * parent_delta[i]
        for i in range(len(delta_indices))
    ]
    delta_cov = independent_residual_covariance(
        submatrix(parent_cov, delta_indices),
        submatrix(child_cov, delta_indices),
        delta_ratio,
    )
    delta_score = generalized_covariance_score(delta_residual, delta_cov)

    slope_index = FEATURE_ORDER.index("mean_slope")
    root_index = FEATURE_ORDER.index("root_gap_lineage")
    slope_raw = scalar_score(
        float(parent["mean_slope"]),
        float(child["mean_slope"]),
        float(parent_cov[slope_index][slope_index]),
        float(child_cov[slope_index][slope_index]),
        slope_baseline,
    )
    slope_frozen = scalar_score(
        float(parent["mean_slope"]),
        float(child["mean_slope"]),
        float(parent_cov[slope_index][slope_index]),
        float(child_cov[slope_index][slope_index]),
        slope_corrected,
        ratio_se=slope_corrected_se,
    )
    root_raw = scalar_score(
        float(parent["root_gap_lineage"]),
        float(child["root_gap_lineage"]),
        float(parent_cov[root_index][root_index]),
        float(child_cov[root_index][root_index]),
        root_baseline,
    )
    root_frozen = scalar_score(
        float(parent["root_gap_lineage"]),
        float(child["root_gap_lineage"]),
        float(parent_cov[root_index][root_index]),
        float(child_cov[root_index][root_index]),
        root_corrected,
        ratio_se=root_corrected_se,
    )

    p48_names = ("P4_S", "P4_D", "P4_S_prime", "P4_D_prime")
    p48_indices = [FEATURE_ORDER.index(name) for name in p48_names]
    p48_rows = []
    for name, index in zip(p48_names, p48_indices):
        ratio = math.pow(2.0, -P48_POWERS[name])
        p48_rows.append(
            {
                "metric": name,
                **scalar_score(
                    float(parent[name]),
                    float(child[name]),
                    float(parent_cov[index][index]),
                    float(child_cov[index][index]),
                    ratio,
                ),
            }
        )

    return {
        "schema": "matching-one/P50-N145-N290-fullcurve-score/v1",
        "status": "frozen primary full-curve score; independent parent/child streams",
        "scoring_order": prediction["scoring_order"],
        "observations": {
            "N145": {key: parent[key] for key in (*FEATURE_ORDER, "p0")},
            "N290": {key: child[key] for key in (*FEATURE_ORDER, "p0")},
        },
        "primary_deltaM_transfer": {
            "levels": list(LEVELS),
            "child_over_parent_ratio": delta_ratio,
            "residual": delta_residual,
            "residual_covariance": delta_cov,
            **delta_score,
        },
        "slope": {
            "raw_asymptotic_baseline": slope_raw,
            "frozen_scalar_plus_H4_correction": slope_frozen,
        },
        "root": {
            "raw_minus_one_quarter_baseline": root_raw,
            "frozen_induced_ratio": {
                **root_frozen,
                "uncertainty_scope": (
                    "target sampling plus frozen slope-source linearized SE; "
                    "does not yet include the separate DeltaM-source uncertainty noted in the prediction artifact"
                ),
            },
        },
        "p48_diagnostics": p48_rows,
        "covariance_rule": (
            "N145 and N290 are jackknifed internally and treated as independent; "
            "residual covariance is Cov_child + ratio^2 Cov_parent"
        ),
        "provenance": provenance,
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-hist", type=Path, required=True)
    parser.add_argument("--parent-metadata", type=Path, required=True)
    parser.add_argument("--child-hist", type=Path, required=True)
    parser.add_argument("--child-metadata", type=Path, required=True)
    parser.add_argument(
        "--prediction",
        type=Path,
        default=root / "predictions/p49_slope_two_sector_145_290_20260828.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    prediction = yaml.safe_load(args.prediction.read_text(encoding="utf-8"))
    if prediction.get("status") != "frozen_before_fullcurve_N290_reveal":
        raise ValueError("prediction artifact is not the frozen pre-reveal contract")
    parent_meta = load_metadata(args.parent_metadata)
    child_meta = load_metadata(args.child_metadata)
    parent_data = read_one_size(args.parent_hist, PARENT_N)
    child_data = read_one_size(args.child_hist, CHILD_N)
    payload = render(
        parent_data,
        child_data,
        parent_meta,
        child_meta,
        prediction,
        {
            "prediction": str(args.prediction),
            "prediction_sha256": sha256(args.prediction),
            "parent_hist": str(args.parent_hist),
            "parent_hist_sha256": sha256(args.parent_hist),
            "parent_metadata": str(args.parent_metadata),
            "parent_metadata_sha256": sha256(args.parent_metadata),
            "child_hist": str(args.child_hist),
            "child_hist_sha256": sha256(args.child_hist),
            "child_metadata": str(args.child_metadata),
            "child_metadata_sha256": sha256(args.child_metadata),
        },
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
