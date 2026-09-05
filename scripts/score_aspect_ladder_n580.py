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
import hashlib
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
             batches: int, seed: str, replica_offset: int,
             engine_threads: int = 0) -> tuple[Path, float]:
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
    if engine_threads > 0:
        # Threads do not change the output: batches get disjoint replica ranges
        # and each batch's row is written once by index, so this only sets how
        # many cores the wall time is spread over.
        command += ["--threads", str(engine_threads)]
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


def measure_rung(binary: Path, rungs: Mapping[int, Mapping[str, Any]], aspect: int,
                 workdir: Path, samples: int, batches: int, seed: str,
                 replica_offset: int, engine_threads: int = 0) -> dict[str, Any]:
    """Run one rung and return its measured block, `_deleted` included."""
    rung = rungs[aspect]
    path, seconds = run_rung(
        binary, rung, workdir / f"n580_r{aspect}", samples, batches, seed, replica_offset,
        engine_threads
    )
    by_orientation = grouped(path)
    full = project_size(by_orientation)
    deleted = [project_size(by_orientation, omitted=b) for b in range(batches)]
    expected = float(Fraction(rung["delta_cos4"]))
    if abs(full["delta_cos4"] - expected) > 1e-12:
        raise ValueError(
            f"r={aspect}: engine leverage {full['delta_cos4']} is not the design's {expected}"
        )
    return {
        "aspect_ratio": aspect,
        "modulus": rung["modulus"],
        "gaussian_norm": rung["gaussian_norm"],
        "samples": samples,
        "seconds": round(seconds, 1),
        "delta_cos4": full["delta_cos4"],
        "spin8_leakage": rung["spin8_leakage"],
        "histogram_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "channels": {
            key: {
                "value": full[key],
                "standard_error": jackknife_se(full[key], [d[key] for d in deleted]),
            }
            for key in REPORTED
        },
        "_deleted": [{key: d[key] for key in REPORTED} for d in deleted],
    }


def finalize(measured: Mapping[str, Any], seed: str, replica_offset: int) -> dict[str, Any]:
    """Score a complete three-rung measured block. The single source of the
    ratio, comparison and verdict logic -- the split-run path and the one-shot
    path both come through here."""
    batches = len(measured["r1"]["_deleted"])

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

    payload_measured = {
        label: {k: v for k, v in row.items() if k != "_deleted"}
        for label, row in measured.items()
    }

    return {
        "schema": SCHEMA,
        "site_count": SITE_COUNT,
        "design_artifact": "results/aspect-ladder-design/latest.json",
        "frozen_design": "predictions/aspect_ladder_n580_20260905.yaml",
        "batches": batches,
        "seed": seed,
        "samples_per_rung": measured["r1"]["samples"],
        "replica_offset": replica_offset,
        "rungs": payload_measured,
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


SHARD_SCHEMA = "matching-one.aspect-ladder-n580-shard.v1"


def write_shard(path: Path, measured: Mapping[str, Any], samples: int, batches: int,
                seed: str, replica_offset: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": SHARD_SCHEMA,
        "aspect": measured["aspect_ratio"],
        "samples": samples,
        "batches": batches,
        "seed": seed,
        "replica_offset": replica_offset,
        "measured": dict(measured),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_shards(shard_dir: Path) -> tuple[dict[str, Any], int, int, str, int]:
    """Reassemble the three rungs' measured blocks from per-box shards, after
    checking the runs are one experiment and not three. Returns the measured
    blocks plus the shared run parameters (samples, batches, seed, offset)."""
    measured: dict[str, Any] = {}
    fixed: tuple[int, int, str, int] | None = None
    for aspect in (1, 2, 4):
        path = shard_dir / f"rung_r{aspect}.json"
        if not path.is_file():
            raise FileNotFoundError(f"missing shard {path}")
        shard = json.loads(path.read_text(encoding="utf-8"))
        if shard["schema"] != SHARD_SCHEMA or shard["aspect"] != aspect:
            raise ValueError(f"{path} is not a shard for r={aspect}")
        key = (shard["samples"], shard["batches"], shard["seed"], shard["replica_offset"])
        if fixed is None:
            fixed = key
        elif key != fixed:
            raise ValueError(
                f"{path} was run under {key}, not the other shards' {fixed}"
            )
        measured[f"r{aspect}"] = shard["measured"]
    assert fixed is not None
    samples, batches, seed, replica_offset = fixed
    return measured, samples, batches, seed, replica_offset


def score(binary: Path, workdir: Path, samples: int, batches: int, seed: str,
          replica_offset: int, engine_threads: int = 0) -> dict[str, Any]:
    rungs = load_design()
    workdir.mkdir(parents=True, exist_ok=True)
    measured = {
        f"r{aspect}": measure_rung(
            binary, rungs, aspect, workdir, samples, batches, seed, replica_offset,
            engine_threads
        )
        for aspect in (1, 2, 4)
    }
    return finalize(measured, seed, replica_offset)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "tr_period")
    parser.add_argument("--workdir", type=Path)
    parser.add_argument("--samples", type=int,
                        help="samples per rung; the same count on all three")
    parser.add_argument("--batches", type=int)
    parser.add_argument("--seed",
                        help="integer seed; the engine rejects non-numeric values")
    parser.add_argument("--replica-offset", type=int,
                        help="must be disjoint from the pilot's 800000000")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rungs",
                        help="comma-separated subset of {1,2,4} to run, e.g. --rungs 4. "
                             "Writes one shard per rung into --rung-dir and scores nothing")
    parser.add_argument("--rung-dir", type=Path,
                        help="directory for per-rung shard files (default: <workdir>/shards)")
    parser.add_argument("--merge-from", type=Path,
                        help="directory holding rung_r{1,2,4}.json shards; scores the "
                             "merged experiment without running the engine")
    parser.add_argument("--engine-threads", type=int, default=0,
                        help="pass --threads T to the engine (0 = engine default). "
                             "Threads do not change the output, only the wall time")
    arguments = parser.parse_args()

    if arguments.merge_from is not None:
        if arguments.rungs or arguments.workdir or arguments.samples:
            parser.error("--merge-from cannot be combined with running options")
        measured, samples, batches, seed, replica_offset = load_shards(arguments.merge_from)
        payload = finalize(measured, seed, replica_offset)
        payload["merged_from_shards"] = str(arguments.merge_from)
    elif arguments.rungs:
        if arguments.workdir is None or arguments.samples is None or \
                arguments.batches is None or arguments.seed is None or \
                arguments.replica_offset is None:
            parser.error("--rungs needs --workdir --samples --batches --seed --replica-offset")
        aspects = parse_rung_subset(arguments.rungs)
        rungs = load_design()
        workdir = arguments.workdir
        workdir.mkdir(parents=True, exist_ok=True)
        shard_dir = arguments.rung_dir or (workdir / "shards")
        for aspect in aspects:
            measured = measure_rung(
                arguments.binary, rungs, aspect, workdir, arguments.samples,
                arguments.batches, arguments.seed, arguments.replica_offset,
                arguments.engine_threads
            )
            write_shard(shard_dir / f"rung_r{aspect}.json", measured,
                        arguments.samples, arguments.batches, arguments.seed,
                        arguments.replica_offset)
            print(f"rung r={aspect} shard written to {shard_dir / f'rung_r{aspect}.json'}",
                  file=sys.stderr)
        return 0
    else:
        if arguments.workdir is None or arguments.samples is None or \
                arguments.batches is None or arguments.seed is None or \
                arguments.replica_offset is None:
            parser.error("needs --workdir --samples --batches --seed --replica-offset")
        payload = score(arguments.binary, arguments.workdir, arguments.samples,
                        arguments.batches, arguments.seed, arguments.replica_offset,
                        arguments.engine_threads)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


def parse_rung_subset(text: str) -> list[int]:
    aspects = sorted({int(part) for part in text.split(",") if part.strip()})
    if not aspects or not set(aspects) <= {1, 2, 4}:
        raise ValueError(f"--rungs must be a subset of 1,2,4; got {text!r}")
    return aspects


if __name__ == "__main__":
    raise SystemExit(main())
