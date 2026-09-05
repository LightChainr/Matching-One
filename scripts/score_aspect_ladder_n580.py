#!/usr/bin/env python3
"""Run and score the N=580 aspect ladder r = 1, 2, 4.

Three tori, all with 580 sites, differing only in modulus: a square torus
(modulus i), a 2:1 rectangle (2i) and a 4:1 rectangle (4i).  Each rung is one
paired run at two lattice orientations of the same Gaussian norm; the frozen P48
projector ``P4_S_prime`` divides the orientation difference by the leverage
``8064/4205``, which is identical on all three rungs.  The score is the pair of
amplitude ratios ``(A4(2i)/A4(i), A4(4i)/A4(i))``.

The design is frozen in ``predictions/aspect_ladder_n580_20260905.yaml`` and
generated into ``results/aspect-ladder-design/latest.json``.  **This script reads
the geometry from that artifact rather than carrying its own copy**, so a
scoring run cannot silently disagree with the design it claims to execute.  It
chooses no sample count, no seed and no stopping rule.

What it does not do, deliberately:

* it does not pool runs, extend a run, or accept a second pass;
* it does not decide anything from the spin-8 leakage beyond reporting it --
  the leakage cancels to leading order in the r=4/r=1 entry and does not cancel
  in the r=2/r=1 entry, and both facts are carried into the output;
* it never reports a competing prediction as confirmed.  ``compatible_at_3_sigma``
  with one survivor is the strongest verdict available, and it is called
  "decisive" only in the sense that the others are out.
"""

from __future__ import annotations

import argparse
import collections
import json
import math
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p48_retrospective import project_size, read_histograms  # noqa: E402

SCHEMA = "matching-one.aspect-ladder-n580.v1"
DESIGN = ROOT / "results" / "aspect-ladder-design" / "latest.json"
SITE_COUNT = 580
REPORTED = ("P4_S_prime", "P4_S", "P4_D_prime", "Mbar_prime")

# From the frozen prediction file.  Each entry is (r=2 value, r=4 value).
COMPETING = {
    "q4_jordan_weight4": (2.75, 10.9908008589),
    "bare_aspect_ratio": (2.0, 4.0),
    "plain_area_scaling": (4.0, 16.0),
    "no_modulus_dependence": (1.0, 1.0),
    "weight12_delta": (0.125, 0.0000279010739704),
    "weight12_E4_cubed": (20.796875, 1327.6635036),
    "weight12_E12": (32.515625, 2080.30719731),
    "weight8_E8": (7.5625, 120.79770352),
}


def load_design() -> dict[int, dict[str, Any]]:
    """The three rungs, keyed by aspect ratio, checked before anything runs."""
    payload = json.loads(DESIGN.read_text(encoding="utf-8"))
    if payload["site_count"] != SITE_COUNT:
        raise ValueError(f"design is for {payload['site_count']} sites, not {SITE_COUNT}")
    rungs = {row["aspect_ratio"]: row for row in payload["rungs"]}
    if set(rungs) != {1, 2, 4}:
        raise ValueError(f"expected rungs 1, 2, 4; design has {sorted(rungs)}")
    leverages = {row["delta_cos4"] for row in rungs.values()}
    if len(leverages) != 1:
        raise ValueError(f"rungs do not share one leverage: {leverages}")
    for aspect, row in rungs.items():
        for key in ("first_matrix_row_major", "second_matrix_row_major"):
            a, b, c, d = row[key]
            if a * d - b * c != SITE_COUNT:
                raise ValueError(f"r={aspect} {key} has determinant {a * d - b * c}")
    return rungs


def jackknife_se(full: float, deleted: Sequence[float]) -> float:
    batches = len(deleted)
    if batches < 2:
        raise ValueError("delete-one jackknife needs at least two batches")
    pseudo = [batches * full - (batches - 1) * value for value in deleted]
    mean = math.fsum(pseudo) / batches
    return math.sqrt(
        math.fsum((value - mean) ** 2 for value in pseudo) / (batches * (batches - 1))
    )


def run_rung(binary: Path, rung: Mapping[str, Any], prefix: Path, samples: int,
             batches: int, seed: str, replica_offset: int) -> tuple[Path, float]:
    command = [
        str(binary),
        "--first-matrix", *(str(v) for v in rung["first_matrix_row_major"]),
        "--second-matrix", *(str(v) for v in rung["second_matrix_row_major"]),
        "--first-rep", *(str(v) for v in rung["first_rep"]),
        "--second-rep", *(str(v) for v in rung["second_rep"]),
        "--samples", str(samples),
        "--batches", str(batches),
        "--seed", seed,
        "--replica-offset", str(replica_offset),
        "--output-prefix", str(prefix),
    ]
    started = time.time()
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        # Surface what the engine said. Swallowing it costs an hour of guessing
        # on a run that is otherwise five hours long.
        raise RuntimeError(
            f"engine failed for r={rung['aspect_ratio']} (exit {completed.returncode})\n"
            f"  command: {' '.join(command)}\n"
            f"  stderr: {completed.stderr.strip() or '(empty)'}\n"
            f"  stdout: {completed.stdout.strip() or '(empty)'}"
        )
    return Path(f"{prefix}.hist.csv"), time.time() - started


def grouped(path: Path):
    histograms = read_histograms(path)
    by_orientation = collections.defaultdict(list)
    for (_, orientation, _), record in sorted(histograms.items()):
        by_orientation[orientation].append(record)
    return by_orientation


def score(binary: Path, workdir: Path, samples: int, batches: int, seed: str,
          replica_offset: int) -> dict[str, Any]:
    rungs = load_design()
    workdir.mkdir(parents=True, exist_ok=True)
    measured: dict[str, Any] = {}

    for aspect in (1, 2, 4):
        rung = rungs[aspect]
        path, seconds = run_rung(
            binary, rung, workdir / f"n580_r{aspect}", samples, batches, seed, replica_offset
        )
        by_orientation = grouped(path)
        full = project_size(by_orientation)
        deleted = [project_size(by_orientation, omitted=b) for b in range(batches)]
        expected = float(Fraction(rung["delta_cos4"]))
        if abs(full["delta_cos4"] - expected) > 1e-12:
            raise ValueError(
                f"r={aspect}: engine leverage {full['delta_cos4']} is not the design's {expected}"
            )
        measured[f"r{aspect}"] = {
            "aspect_ratio": aspect,
            "modulus": rung["modulus"],
            "gaussian_norm": rung["gaussian_norm"],
            "samples": samples,
            "seconds": round(seconds, 1),
            "delta_cos4": full["delta_cos4"],
            "spin8_leakage": rung["spin8_leakage"],
            "channels": {
                key: {
                    "value": full[key],
                    "standard_error": jackknife_se(full[key], [d[key] for d in deleted]),
                }
                for key in REPORTED
            },
            "_deleted": [{key: d[key] for key in REPORTED} for d in deleted],
        }

    ratios: dict[str, dict[str, Any]] = {}
    for label, numerator in (("r2_over_r1", "r2"), ("r4_over_r1", "r4")):
        ratios[label] = {}
        for key in REPORTED:
            denominator = measured["r1"]["channels"][key]["value"]
            if denominator == 0.0:
                continue
            value = measured[numerator]["channels"][key]["value"] / denominator
            deleted = [
                measured[numerator]["_deleted"][b][key] / measured["r1"]["_deleted"][b][key]
                for b in range(batches)
                if measured["r1"]["_deleted"][b][key] != 0.0
            ]
            ratios[label][key] = {
                "value": value,
                "standard_error": jackknife_se(value, deleted),
            }

    primary = {
        label: ratios[label]["P4_S_prime"] for label in ("r2_over_r1", "r4_over_r1")
    }
    comparison = {}
    for name, (at_two, at_four) in COMPETING.items():
        comparison[name] = {
            "predicted": {"r2_over_r1": at_two, "r4_over_r1": at_four},
            "z": {
                label: (primary[label]["value"] - predicted) / primary[label]["standard_error"]
                for label, predicted in (("r2_over_r1", at_two), ("r4_over_r1", at_four))
            },
        }

    # The r=4 entry is the discriminator and the clean one; the r=2 entry keeps a
    # spin-8 systematic, so a competitor is only excluded on the strength of r=4.
    excluded = sorted(n for n, row in comparison.items() if abs(row["z"]["r4_over_r1"]) >= 3.0)
    compatible = sorted(n for n in comparison if n not in excluded)

    for row in measured.values():
        del row["_deleted"]

    return {
        "schema": SCHEMA,
        "site_count": SITE_COUNT,
        "design_artifact": "results/aspect-ladder-design/latest.json",
        "frozen_design": "predictions/aspect_ladder_n580_20260905.yaml",
        "batches": batches,
        "seed": seed,
        "samples_per_rung": samples,
        "replica_offset": replica_offset,
        "rungs": measured,
        "ratios": ratios,
        "primary_channel": "P4_S_prime",
        "discriminating_entry": "r4_over_r1",
        "why_that_entry": (
            "the r=1 and r=4 rungs carry the same spin-8 leakage, so it cancels to "
            "leading order in their ratio. It does not cancel in r2_over_r1, where "
            "the r=2 rung carries the opposite sign, so that entry keeps a "
            "systematic of about -0.055 * (A8/A4 + A8'/A4') and is reported but not "
            "used to exclude"
        ),
        "comparison": comparison,
        "compatible_at_3_sigma": compatible,
        "excluded_at_3_sigma": excluded,
        "verdict": (
            "underpowered: more than one competing prediction survives at r4_over_r1"
            if len(compatible) > 1
            else ("decisive: exactly one competing prediction survives"
                  if len(compatible) == 1
                  else "no competing prediction survives")
        ),
        "not_established": [
            "identification of the Q4 Jordan module",
            "that the measured amplitude is the log slope rather than a leading "
            "amplitude with the same symmetry",
            "the C + A4 cos4 + A8 cos8 form itself, which two orientations "
            "determine exactly and therefore cannot check",
            "freedom from spin-8 in the r2_over_r1 entry",
            "any modular interpretation of the bare aspect ratio, should it be the "
            "survivor. A law that fits three points is a law that fits three points",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "tr_period")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--samples", type=int, required=True,
                        help="samples per rung; the same count on all three")
    parser.add_argument("--batches", type=int, required=True)
    parser.add_argument("--seed", required=True,
                        help="integer seed; the engine rejects non-numeric values")
    parser.add_argument("--replica-offset", type=int, required=True,
                        help="must be disjoint from the pilot's 800000000")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    payload = score(arguments.binary, arguments.workdir, arguments.samples,
                    arguments.batches, arguments.seed, arguments.replica_offset)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
