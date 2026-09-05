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

The three rungs may be run on separate machines (one rung per box) and merged,
because they share a seed and a replica offset: batch ``b`` of one rung and
batch ``b`` of another are the same random block seen through two geometries,
which is what makes the cross-rung covariance a measurement rather than a
coincidence.  ``--rungs`` writes one per-rung shard that carries the delete-one
replicates; ``--merge-from`` scores the merged experiment without running the
engine.

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


def _pseudovalues(full: float, deleted: Sequence[float]) -> list[float]:
    batches = len(deleted)
    return [batches * full - (batches - 1) * value for value in deleted]


def jackknife_covariance(full_x: float, deleted_x: Sequence[float],
                         full_y: float, deleted_y: Sequence[float]) -> float:
    """Paired delete-one covariance of two channels sampled on the same batches.

    The rungs share a seed and a replica offset, so batch b of one rung and
    batch b of another are the same random block seen through two geometries.
    Deleting them together is what makes this a covariance rather than a
    coincidence, and ``load_shards`` is what enforces the sharing.
    """
    batches = len(deleted_x)
    if batches != len(deleted_y):
        raise ValueError("paired jackknife needs the same batches on both channels")
    if batches < 2:
        raise ValueError("delete-one jackknife needs at least two batches")
    px = _pseudovalues(full_x, deleted_x)
    py = _pseudovalues(full_y, deleted_y)
    mx = math.fsum(px) / batches
    my = math.fsum(py) / batches
    return math.fsum((a - mx) * (b - my) for a, b in zip(px, py)) / (batches * (batches - 1))


def fieller_z(numerator: float, var_numerator: float,
              denominator: float, var_denominator: float,
              covariance: float, predicted_ratio: float) -> float:
    """Test ratio == predicted_ratio on the contrast Y - R0*X, not on Y/X.

    Dividing first and then forming (Rhat - R0)/SE(Rhat) assumes Rhat is
    roughly normal with a spread that does not depend on R0.  Both assumptions
    fail when the denominator is only a few sigma from zero: the ratio of two
    normals has Cauchy-like tails, and the spread grows as R0 moves out toward
    the region where the denominator could have been small.  The contrast has
    neither problem -- X and Y are means over the batches, so the CLT applies to
    them directly -- and it is the standard (Fieller) treatment.
    """
    residual = numerator - predicted_ratio * denominator
    variance = (var_numerator
                + predicted_ratio * predicted_ratio * var_denominator
                - 2.0 * predicted_ratio * covariance)
    if variance <= 0.0:
        raise ValueError("non-positive variance in the Fieller contrast")
    return residual / math.sqrt(variance)


def fieller_interval(numerator: float, var_numerator: float,
                     denominator: float, var_denominator: float,
                     covariance: float, sigma: float = 3.0):
    """The Fieller confidence set for the ratio, or None when it is unbounded.

    The set is bounded exactly when the denominator itself clears ``sigma``.
    When it does not, every large ratio is compatible with the data and the
    honest report is that there is no upper limit -- not a wide interval.
    """
    k = sigma * sigma
    a = denominator * denominator - k * var_denominator
    if a <= 0.0:
        return None
    b = -2.0 * denominator * numerator + 2.0 * k * covariance
    c = numerator * numerator - k * var_numerator
    discriminant = b * b - 4.0 * a * c
    if discriminant < 0.0:
        return None
    root = math.sqrt(discriminant)
    return ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a))


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
            top = measured[numerator]["channels"][key]["value"]
            deleted_top = [measured[numerator]["_deleted"][b][key] for b in range(batches)]
            deleted_bottom = [measured["r1"]["_deleted"][b][key] for b in range(batches)]
            deleted = [
                a / b for a, b in zip(deleted_top, deleted_bottom) if b != 0.0
            ]
            var_top = jackknife_se(top, deleted_top) ** 2
            var_bottom = jackknife_se(denominator, deleted_bottom) ** 2
            covariance = jackknife_covariance(top, deleted_top, denominator, deleted_bottom)
            ratios[label][key] = {
                "value": top / denominator,
                "standard_error": jackknife_se(top / denominator, deleted),
                # Kept so the deciding statistic can be recomputed, and the run
                # re-tested against a competitor nobody listed, without the
                # histograms.  The first version of this script threw these away
                # and the run could not be re-scored at all.
                "numerator_value": top,
                "numerator_variance": var_top,
                "denominator_value": denominator,
                "denominator_variance": var_bottom,
                "covariance": covariance,
                "denominator_sigma_from_zero": (
                    abs(denominator) / math.sqrt(var_bottom) if var_bottom > 0.0 else math.inf
                ),
            }

    primary = {
        label: ratios[label]["P4_S_prime"] for label in ("r2_over_r1", "r4_over_r1")
    }
    comparison = {}
    for name, (at_two, at_four) in COMPETING.items():
        row = {"predicted": {"r2_over_r1": at_two, "r4_over_r1": at_four}, "z": {}, "ratio_z_not_used": {}}
        for label, predicted in (("r2_over_r1", at_two), ("r4_over_r1", at_four)):
            entry = primary[label]
            row["z"][label] = fieller_z(
                entry["numerator_value"], entry["numerator_variance"],
                entry["denominator_value"], entry["denominator_variance"],
                entry["covariance"], predicted,
            )
            row["ratio_z_not_used"][label] = (
                (entry["value"] - predicted) / entry["standard_error"]
            )
        comparison[name] = row

    # The r=4 entry is the discriminator and the clean one; the r=2 entry keeps a
    # spin-8 systematic, so a competitor is only excluded on the strength of r=4.
    excluded = sorted(n for n, row in comparison.items() if abs(row["z"]["r4_over_r1"]) >= 3.0)
    compatible = sorted(n for n in comparison if n not in excluded)

    intervals = {}
    for label in ("r2_over_r1", "r4_over_r1"):
        entry = primary[label]
        bounds = fieller_interval(
            entry["numerator_value"], entry["numerator_variance"],
            entry["denominator_value"], entry["denominator_variance"],
            entry["covariance"],
        )
        intervals[label] = {"lower": bounds[0], "upper": bounds[1]} if bounds else {
            "unbounded_because": (
                "the denominator does not clear 3 sigma, so arbitrarily large "
                "ratios are compatible with the data and there is no upper limit"
            )
        }

    # The whole response vector and its whole covariance, before any ratio is
    # taken.  The first version of this script stored ratios and pairwise
    # pieces, which threw away cov(r2, r4) and made the three-rung data
    # untestable against anything except the pair the script happened to form.
    # A model that predicts proportions predicts a ray; storing only ratios
    # quotients the amplitude away before anyone can ask a different question.
    rung_order = ["r1", "r2", "r4"]
    response: dict[str, Any] = {"rung_order": rung_order, "channels": {}}
    for key in REPORTED:
        full = [measured[label]["channels"][key]["value"] for label in rung_order]
        deleted = [[measured[label]["_deleted"][b][key] for b in range(batches)]
                   for label in rung_order]
        covariance = [
            [jackknife_covariance(full[i], deleted[i], full[j], deleted[j])
             for j in range(len(rung_order))]
            for i in range(len(rung_order))
        ]
        response["channels"][key] = {"vector": full, "covariance": covariance}

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
        "response": response,
        "ratios": ratios,
        "fieller_interval_3sigma": intervals,
        "primary_channel": "P4_S_prime",
        "discriminating_entry": "r4_over_r1",
        "test_statistic": {
            "used": "fieller_contrast",
            "definition": "z = (Y - R0*X) / sqrt(varY + R0^2 varX - 2 R0 cov(X,Y))",
            "why": (
                "the denominator A4(i) is only a few sigma from zero, so the ratio "
                "is not approximately normal and its standard error at the observed "
                "point understates the spread out where the large predictions sit. "
                "X and Y are means over batches, so the CLT applies to the contrast "
                "directly. The pre-registered statistic was (Rhat - R0)/SE(Rhat), "
                "reported alongside as ratio_z_not_used and NOT used to decide"
            ),
            "deviation_from_the_frozen_runner": (
                "the ratio z-test was committed in the runner before the run and is "
                "therefore effectively part of the freeze. Changing it after seeing "
                "the data is a real deviation and is recorded here rather than "
                "silently applied. It moves verdicts in both directions"
            ),
        },
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
