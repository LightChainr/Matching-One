#!/usr/bin/env python3
"""Run and score the N=290 square-vs-rectangular spin-4 modulus fingerprint.

Four period matrices, all with 290 sites: a square torus (modulus i) and a 2:1
rectangular torus (modulus 2i), each realized at two lattice orientations whose
cos4theta values are the same pair with the roles exchanged.  Within each family
the frozen P48 projector `P4_S_prime` divides the orientation difference by the
identical rational leverage 8064/4205, leaving the spin-4 amplitude.  The score
is the ratio of the two amplitudes.

The prediction under the Q4-Jordan weight-4 shape is 11/4; plain area scaling
gives 4; no modulus dependence gives 1.  The design is frozen in
`predictions/modulus_fingerprint_n290_*.yaml` and this script does not choose
sample counts, seeds or stopping.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p48_retrospective import project_size, read_histograms  # noqa: E402

SCHEMA = "matching-one.modulus-fingerprint-n290.v1"

FAMILIES = {
    "square": {
        "modulus": "i",
        "first_matrix": (17, -1, 1, 17),
        "second_matrix": (13, -11, 11, 13),
        "first_rep": (17, 1),
        "second_rep": (13, 11),
    },
    "rectangular": {
        "modulus": "2i",
        "first_matrix": (12, -2, 1, 24),
        "second_matrix": (9, -16, 8, 18),
        "first_rep": (12, 1),
        "second_rep": (9, 8),
    },
}

COMPETING = {
    "q4_jordan_weight4": 11.0 / 4.0,
    "plain_area_scaling": 4.0,
    "no_modulus_dependence": 1.0,
    "weight12_delta": 0.125,
    "weight12_E4_cubed": 20.796875,
    "weight12_E12": 32.515625,
}

REPORTED = ("P4_S_prime", "P4_S", "P4_D_prime", "Mbar_prime")


def jackknife_se(full: float, deleted: Sequence[float]) -> float:
    batches = len(deleted)
    if batches < 2:
        raise ValueError("delete-one jackknife needs at least two batches")
    pseudo = [batches * full - (batches - 1) * value for value in deleted]
    mean = math.fsum(pseudo) / batches
    return math.sqrt(
        math.fsum((value - mean) ** 2 for value in pseudo) / (batches * (batches - 1))
    )


def run_family(
    binary: Path,
    family: Mapping[str, object],
    prefix: Path,
    samples: int,
    batches: int,
    seed: str,
    replica_offset: int,
) -> tuple[Path, float]:
    command = [
        str(binary),
        "--first-matrix", *(str(v) for v in family["first_matrix"]),
        "--second-matrix", *(str(v) for v in family["second_matrix"]),
        "--first-rep", *(str(v) for v in family["first_rep"]),
        "--second-rep", *(str(v) for v in family["second_rep"]),
        "--samples", str(samples),
        "--batches", str(batches),
        "--seed", seed,
        "--replica-offset", str(replica_offset),
        "--output-prefix", str(prefix),
    ]
    started = time.time()
    subprocess.run(command, check=True, capture_output=True)
    return Path(f"{prefix}.hist.csv"), time.time() - started


def grouped(path: Path):
    histograms = read_histograms(path)
    by_orientation = collections.defaultdict(list)
    for (_, orientation, _), record in sorted(histograms.items()):
        by_orientation[orientation].append(record)
    return by_orientation


def score(
    binary: Path,
    workdir: Path,
    samples: Mapping[str, int],
    batches: int,
    seed: str,
    replica_offset: int,
) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)
    measured = {}
    for name, family in FAMILIES.items():
        path, seconds = run_family(
            binary, family, workdir / f"n290_{name}", samples[name], batches, seed,
            replica_offset,
        )
        by_orientation = grouped(path)
        full = project_size(by_orientation)
        deleted = [project_size(by_orientation, omitted=b) for b in range(batches)]
        measured[name] = {
            "modulus": family["modulus"],
            "samples": samples[name],
            "seconds": round(seconds, 1),
            "delta_cos4": full["delta_cos4"],
            "channels": {
                key: {
                    "value": full[key],
                    "standard_error": jackknife_se(full[key], [d[key] for d in deleted]),
                }
                for key in REPORTED
            },
            "_deleted": [{key: d[key] for key in REPORTED} for d in deleted],
        }

    ratios = {}
    for key in REPORTED:
        numerator = measured["rectangular"]["channels"][key]["value"]
        denominator = measured["square"]["channels"][key]["value"]
        if denominator == 0.0:
            continue
        full = numerator / denominator
        deleted = [
            measured["rectangular"]["_deleted"][b][key]
            / measured["square"]["_deleted"][b][key]
            for b in range(batches)
            if measured["square"]["_deleted"][b][key] != 0.0
        ]
        ratios[key] = {
            "value": full,
            "standard_error": jackknife_se(full, deleted),
        }

    primary = ratios["P4_S_prime"]
    against = {
        name: {
            "predicted": value,
            "z": (primary["value"] - value) / primary["standard_error"],
        }
        for name, value in COMPETING.items()
    }
    compatible = sorted(name for name, row in against.items() if abs(row["z"]) < 3.0)
    excluded = sorted(name for name, row in against.items() if abs(row["z"]) >= 3.0)

    for name in measured:
        del measured[name]["_deleted"]

    return {
        "schema": SCHEMA,
        "site_count": 290,
        "batches": batches,
        "seed": seed,
        "replica_offset": replica_offset,
        "families": measured,
        "ratios_rectangular_over_square": ratios,
        "primary_channel": "P4_S_prime",
        "comparison": against,
        "compatible_at_3_sigma": compatible,
        "excluded_at_3_sigma": excluded,
        "verdict": (
            "underpowered: more than one competing prediction survives"
            if len(compatible) > 1
            else ("decisive" if len(compatible) == 1 else "no competing prediction survives")
        ),
        "not_established": [
            "identification of the Q4 Jordan module",
            "that the measured amplitude is the log slope rather than a leading "
            "amplitude with the same symmetry",
            "anything separating shapes that carry the same E4 factor",
            "the 11/4 itself, which is conditional on the normalization removing "
            "the same block (docs/astra/Q2-additive-shape-ambiguity.md)",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "tr_period")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--square-samples", type=int, required=True)
    parser.add_argument("--rectangular-samples", type=int, required=True)
    parser.add_argument("--batches", type=int, required=True)
    parser.add_argument("--seed", required=True)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    payload = score(
        arguments.binary,
        arguments.workdir,
        {
            "square": arguments.square_samples,
            "rectangular": arguments.rectangular_samples,
        },
        arguments.batches,
        arguments.seed,
        arguments.replica_offset,
    )
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
