#!/usr/bin/env python3
"""Score the frozen P205 N25/N50/N125 quotient-character prism.

Each equal-N pair is evaluated under its common priority field.  Independent
seeds across N make the three contrast covariance diagonal.  H4, H8 and H12
are scored in the frozen order as one-amplitude lines with the fixed N^-13/8
radial factor; no offset, exponent, or correction coefficient is fitted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp
import yaml

from analyze_matching_parity_derivatives_fast import H, combine, obs, read, remove
from analyze_threshold_rank_orientation import read_histograms, validate_moments


SIZES = (25, 50, 125)
MODEL_ORDER = ("H4", "H8", "H12")
RADIAL_EXPONENT = Fraction(13, 8)


@dataclass(frozen=True)
class PairRun:
    n: int
    prefix: Path
    histogram_path: Path
    moments_path: Path
    metadata_path: Path
    metadata: Mapping[str, object]
    data: Mapping[tuple[int, str, int], H]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_contract(path: Path) -> dict[str, object]:
    design = yaml.safe_load(path.read_text(encoding="utf-8"))
    if design.get("status") != "design_frozen_after_P205_10M_before_new_target_reveal":
        raise ValueError("design is not the frozen prereveal artifact")
    if design["production_design"]["authorization"] != "frozen_not_launched":
        raise ValueError("frozen acquisition marker changed")
    if tuple(design["frozen_model_order"]) != MODEL_ORDER:
        raise ValueError("frozen model order changed")
    if design["exact_sign_code"] != {
        "H4": ["+", "+", "+"],
        "H8": ["+", "-", "+"],
        "H12": ["+", "+", "-"],
    }:
        raise ValueError("exact character code changed")
    targets = {int(row["N"]): row for row in design["targets"]}
    if tuple(sorted(targets)) != SIZES:
        raise ValueError("frozen target set changed")
    characters = {
        model: {
            n: Fraction(str(targets[n]["exact_character_difference"][model]))
            for n in SIZES
        }
        for model in MODEL_ORDER
    }
    return {
        "design": design,
        "targets": targets,
        "characters": characters,
        "p_ref": mp.mpf(str(design["observable"]["p_ref"])),
    }


def load_run(n: int, prefix: Path, contract: Mapping[str, object]) -> PairRun:
    if n not in SIZES:
        raise ValueError(f"N={n} is not a frozen prism target")
    histogram = Path(f"{prefix}.hist.csv")
    moments = Path(f"{prefix}.moments.csv")
    metadata_path = Path(f"{prefix}.metadata.json")
    data = read(histogram)
    validate_moments(moments, read_histograms(histogram))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    target = contract["targets"][n]
    production = contract["design"]["production_design"]
    design = metadata["designs"][0]
    expected = {
        "N": n,
        "first": [int(value) for value in target["first_representation"]],
        "second": [int(value) for value in target["second_representation"]],
        "first_period_matrix": [
            [int(value) for value in row] for row in target["first_matrix"]
        ],
        "second_period_matrix": [
            [int(value) for value in row] for row in target["second_matrix"]
        ],
        "first_smith_invariants": [
            int(value) for value in target["first_smith_invariants"]
        ],
        "second_smith_invariants": [
            int(value) for value in target["second_smith_invariants"]
        ],
    }
    for field, value in expected.items():
        if design[field] != value:
            raise ValueError(f"N={n}: {field} differs from frozen target")
    required_metadata = {
        "git_commit", "seed", "replica_counter_first",
        "replica_counter_last_exclusive", "samples_per_pair", "batches",
        "threads_requested",
    }
    if required_metadata - set(metadata):
        raise ValueError(f"N={n}: incomplete metadata")
    if metadata["git_commit"] != production["executable_source_commit"]:
        raise ValueError(f"N={n}: runner commit changed")
    if int(metadata["seed"]) != int(target["seed"]):
        raise ValueError(f"N={n}: seed changed")
    if [
        int(metadata["replica_counter_first"]),
        int(metadata["replica_counter_last_exclusive"]),
    ] != [int(value) for value in production["replica_counter"]]:
        raise ValueError(f"N={n}: replica counter domain changed")
    for field, design_field in (
        ("samples_per_pair", "samples_per_pair"),
        ("batches", "batches"),
        ("threads_requested", "threads_per_job"),
    ):
        if int(metadata[field]) != int(production[design_field]):
            raise ValueError(f"N={n}: {field} allocation changed")
    batches = int(metadata["batches"])
    if sorted({key[2] for key in data}) != list(range(batches)):
        raise ValueError(f"N={n}: incomplete batch grid")
    for orientation, representation in (
        ("first", tuple(expected["first"])),
        ("second", tuple(expected["second"])),
    ):
        rows = [row for key, row in data.items() if key[1] == orientation]
        if len(rows) != batches or {(row.a, row.b) for row in rows} != {representation}:
            raise ValueError(f"N={n}: {orientation} descriptor changed")
        if sum(row.samples for row in rows) != int(metadata["samples_per_pair"]):
            raise ValueError(f"N={n}: {orientation} sample total changed")
    return PairRun(n, prefix, histogram, moments, metadata_path, metadata, data)


def aligned_rows(run: PairRun, orientation: str) -> list[H]:
    return [run.data[key] for key in sorted(run.data) if key[1] == orientation]


def jackknife_variance(values: Sequence[mp.mpf]) -> mp.mpf:
    if len(values) < 2:
        raise ValueError("jackknife variance needs at least two values")
    mean = mp.fsum(values) / len(values)
    return mp.mpf(len(values) - 1) / len(values) * mp.fsum(
        (value - mean) ** 2 for value in values
    )


def score_run(run: PairRun, p_ref: mp.mpf) -> dict[str, object]:
    first_rows = aligned_rows(run, "first")
    second_rows = aligned_rows(run, "second")
    first = combine(first_rows)
    second = combine(second_rows)
    point = obs(first, p_ref)["M"] - obs(second, p_ref)["M"]
    deleted = [
        obs(remove(first, first_rows[batch]), p_ref)["M"]
        - obs(remove(second, second_rows[batch]), p_ref)["M"]
        for batch in range(len(first_rows))
    ]
    variance = jackknife_variance(deleted)
    if variance <= 0:
        raise ValueError(f"N={run.n}: pair contrast variance is not positive")
    return {"point": point, "variance": variance, "deleted": deleted}


def model_vector(
    model: str, characters: Mapping[str, Mapping[int, Fraction]]
) -> list[mp.mpf]:
    exponent = mp.mpf(RADIAL_EXPONENT.numerator) / RADIAL_EXPONENT.denominator
    return [
        mp.power(n, -exponent)
        * mp.mpf(characters[model][n].numerator)
        / characters[model][n].denominator
        for n in SIZES
    ]


def score_model(
    points: Sequence[mp.mpf], variances: Sequence[mp.mpf],
    vector: Sequence[mp.mpf],
) -> dict[str, object]:
    information = mp.fsum(vector[i] ** 2 / variances[i] for i in range(3))
    if information <= 0:
        raise ValueError("model vector has no information")
    amplitude = mp.fsum(
        vector[i] * points[i] / variances[i] for i in range(3)
    ) / information
    residuals = [points[i] - amplitude * vector[i] for i in range(3)]
    chi_square = mp.fsum(residuals[i] ** 2 / variances[i] for i in range(3))
    return {
        "fitted_shared_amplitude": mp.nstr(amplitude, 30),
        "amplitude_standard_error": mp.nstr(1 / mp.sqrt(information), 20),
        "model_vector": [mp.nstr(value, 30) for value in vector],
        "residuals": [mp.nstr(value, 30) for value in residuals],
        "signed_standardized_residuals": [
            mp.nstr(residuals[i] / mp.sqrt(variances[i]), 18) for i in range(3)
        ],
        "chi_square": mp.nstr(chi_square, 20),
        "degrees_of_freedom": 2,
        "chi_square_survival_df2": mp.nstr(mp.exp(-chi_square / 2), 18),
    }


def render(
    runs: Mapping[int, PairRun], contract: Mapping[str, object], design_path: Path
) -> dict[str, object]:
    if len({int(runs[n].metadata["seed"]) for n in SIZES}) != len(SIZES):
        raise ValueError("frozen sizes must use independent seeds")
    by_size = {n: score_run(runs[n], contract["p_ref"]) for n in SIZES}
    points = [by_size[n]["point"] for n in SIZES]
    variances = [by_size[n]["variance"] for n in SIZES]
    models = []
    for model in MODEL_ORDER:
        row = score_model(points, variances, model_vector(model, contract["characters"]))
        row["name"] = model
        models.append(row)
    best_chi = min(mp.mpf(row["chi_square"]) for row in models)
    for row in models:
        row["delta_chi_square_from_best"] = mp.nstr(
            mp.mpf(row["chi_square"]) - best_chi, 20
        )
    best = min(models, key=lambda row: mp.mpf(row["chi_square"]))["name"]
    return {
        "schema": "matching-one/p205-quotient-character-prism-score/v1",
        "status": "frozen fixed-sample score; one amplitude per harmonic",
        "fixed_probability": mp.nstr(contract["p_ref"], 30),
        "size_order": list(SIZES),
        "model_order": list(MODEL_ORDER),
        "contrast_vector": [mp.nstr(value, 30) for value in points],
        "contrast_covariance": [
            [mp.nstr(variances[i], 25) if i == j else "0" for j in range(3)]
            for i in range(3)
        ],
        "by_size": {
            str(n): {
                "DeltaM_first_minus_second": mp.nstr(by_size[n]["point"], 30),
                "variance": mp.nstr(by_size[n]["variance"], 25),
                "standard_error": mp.nstr(mp.sqrt(by_size[n]["variance"]), 20),
            }
            for n in SIZES
        },
        "models": models,
        "best_frozen_harmonic_by_chi_square": best,
        "decision": f"best_one_amplitude_character_line_is_{best}",
        "provenance": {
            "design": str(design_path),
            "design_sha256": sha256(design_path),
            "runner_commit": contract["design"]["production_design"][
                "executable_source_commit"
            ],
            "inputs": [
                {
                    "N": n,
                    "histogram": str(runs[n].histogram_path),
                    "histogram_sha256": sha256(runs[n].histogram_path),
                    "moments": str(runs[n].moments_path),
                    "moments_sha256": sha256(runs[n].moments_path),
                    "metadata": str(runs[n].metadata_path),
                    "metadata_sha256": sha256(runs[n].metadata_path),
                }
                for n in SIZES
            ],
        },
        "interpretation_guard": (
            "H4, H8 and H12 were scored in frozen order as exact angular "
            "character lines. Exactly one shared amplitude was fitted per line; "
            "the N exponent remained 13/8 and no offset or correction was fitted."
        ),
    }


def run_spec(value: str) -> tuple[int, Path]:
    fields = value.split(":", 1)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("run must be N:OUTPUT_PREFIX")
    try:
        n = int(fields[0])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("run N must be an integer") from exc
    return n, Path(fields[1])


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--design", type=Path,
        default=root / "predictions/p205_quotient_character_prism_20260829.yaml",
    )
    parser.add_argument("--run", action="append", type=run_spec, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = 80
    contract = load_contract(args.design)
    prefixes = dict(args.run)
    if tuple(sorted(prefixes)) != SIZES:
        raise SystemExit("exactly one --run is required for each N=25,50,125")
    runs = {n: load_run(n, prefixes[n], contract) for n in SIZES}
    payload = render(runs, contract, args.design)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
