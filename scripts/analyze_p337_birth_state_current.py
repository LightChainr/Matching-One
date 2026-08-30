#!/usr/bin/env python3
"""Rebuild P337 birth/completion/collision currents from immutable joint archives."""

from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
from dataclasses import dataclass
from fractions import Fraction
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time
from typing import Any, Iterator, Mapping, Sequence

import yaml


SCHEMA = "matching-one.p337-birth-state-current.v1"
ORIENTATIONS = ("first", "second")
LINE_ORDER = ("x", "d_plus", "d_minus", "y")
SQRT2 = math.sqrt(2.0)
LINE_WEIGHTS = {
    "x": (0.5, 1.0 / SQRT2, 0.0),
    "d_plus": (-0.5, 0.0, 1.0 / SQRT2),
    "d_minus": (-0.5, 0.0, -1.0 / SQRT2),
    "y": (0.5, -1.0 / SQRT2, 0.0),
}
CHARACTERS = ("H", "A", "D")
VECTOR_ORDER = (
    "angular_E_top",
    "angular_A_top",
    "angular_E_log_derivative",
    "angular_A_log_derivative",
    "angular_J01_log",
    "angular_J12_log",
    "angular_J02_log",
    "angular_collision_mass",
    "angular_completion_age1_log",
    "angular_completion_age2_log",
    "angular_age_hazard_beta",
    "angular_K_H_linear",
    "angular_K_A_activity",
    "angular_K_D_activity",
    "angular_H_birth_log",
    "angular_H_completion_log",
    "angular_A_birth_log",
    "angular_A_completion_log",
    "angular_D_birth_log",
    "angular_D_completion_log",
)


def zeros(n: int) -> list[float]:
    return [0.0] * (n + 1)


def nested_zeros(rows: int, n: int) -> list[list[float]]:
    return [zeros(n) for _ in range(rows)]


def projective_line(x: int, y: int) -> str:
    x %= 3
    y %= 3
    if x == 0 and y == 0:
        raise ValueError("zero vector has no projective line")
    if y == 0:
        return "x"
    if x == 0:
        return "y"
    if y == x:
        return "d_plus"
    if y == (-x) % 3:
        return "d_minus"
    raise ValueError(f"unrecognized F3 projective line {(x, y)}")


@dataclass
class BatchBuilder:
    n: int

    def __post_init__(self) -> None:
        self.samples: int | None = None
        self.total_count = 0
        self.rows = 0
        self.h1_line = zeros(self.n)
        self.h2_line = zeros(self.n)
        self.direct = zeros(self.n)
        self.char_birth = nested_zeros(3, self.n)
        self.char_exit = nested_zeros(3, self.n)
        self.age_exit1 = zeros(self.n)
        self.age_exit2 = zeros(self.n)
        width = self.n + 2
        self.risk_diff = [[0.0] * width for _ in LINE_ORDER]
        self.x_k_diff = [[0.0] * width for _ in LINE_ORDER]
        self.x_c_diff = [[0.0] * width for _ in LINE_ORDER]
        self.x2_k2_diff = [[0.0] * width for _ in LINE_ORDER]
        self.x2_k_diff = [[0.0] * width for _ in LINE_ORDER]
        self.x2_c_diff = [[0.0] * width for _ in LINE_ORDER]
        self.exit_y = nested_zeros(len(LINE_ORDER), self.n)
        self.exit_xy = nested_zeros(len(LINE_ORDER), self.n)

    def add(
        self,
        *,
        samples: int,
        k1: int,
        k2: int,
        direct: bool,
        ell_x: int,
        ell_y: int,
        count: int,
    ) -> None:
        if self.samples is None:
            self.samples = samples
        elif self.samples != samples:
            raise ValueError("samples changed inside one batch")
        if not 1 <= k1 <= k2 <= self.n:
            raise ValueError(f"invalid activation pair {(k1, k2)} for N={self.n}")
        if count < 0:
            raise ValueError("negative sparse count")
        self.total_count += count
        self.rows += 1
        if direct:
            if k1 != k2 or ell_x or ell_y:
                raise ValueError("invalid DIRECT_RANK2 atom")
            self.direct[k1] += count
            return
        if k1 >= k2 or math.gcd(abs(ell_x), abs(ell_y)) != 1:
            raise ValueError("invalid LINE atom")
        line = projective_line(ell_x, ell_y)
        line_index = LINE_ORDER.index(line)
        self.h1_line[k1] += count
        self.h2_line[k2] += count
        age = (k2 - k1) / self.n
        self.age_exit1[k2] += count * age
        self.age_exit2[k2] += count * age * age
        for j, weight in enumerate(LINE_WEIGHTS[line]):
            self.char_birth[j][k1] += count * weight
            self.char_exit[j][k2] += count * weight

        start = k1
        stop = k2  # exclusive current-layer risk interval [k1,k2)
        inv_n = 1.0 / self.n
        inv_n2 = inv_n * inv_n
        updates = (
            (self.risk_diff[line_index], count),
            (self.x_k_diff[line_index], count * inv_n),
            (self.x_c_diff[line_index], -count * k1 * inv_n),
            (self.x2_k2_diff[line_index], count * inv_n2),
            (self.x2_k_diff[line_index], -2.0 * count * k1 * inv_n2),
            (self.x2_c_diff[line_index], count * k1 * k1 * inv_n2),
        )
        for array, value in updates:
            array[start] += value
            array[stop] -= value
        current = k2 - 1
        current_age = (current - k1) * inv_n
        self.exit_y[line_index][current] += count
        self.exit_xy[line_index][current] += count * current_age

    def finish(self) -> dict[str, Any]:
        if self.samples is None or self.total_count != self.samples:
            raise ValueError(
                f"batch sparse cells sum to {self.total_count}, expected {self.samples}"
            )
        risk = nested_zeros(len(LINE_ORDER), self.n)
        x = nested_zeros(len(LINE_ORDER), self.n)
        x2 = nested_zeros(len(LINE_ORDER), self.n)
        for line in range(len(LINE_ORDER)):
            running = [0.0] * 6
            sources = (
                self.risk_diff[line],
                self.x_k_diff[line],
                self.x_c_diff[line],
                self.x2_k2_diff[line],
                self.x2_k_diff[line],
                self.x2_c_diff[line],
            )
            for k in range(self.n + 1):
                for j, source in enumerate(sources):
                    running[j] += source[k]
                risk[line][k] = running[0]
                x[line][k] = running[1] * k + running[2]
                x2[line][k] = running[3] * k * k + running[4] * k + running[5]
        return {
            "samples": self.samples,
            "rows": self.rows,
            "h1_line": self.h1_line,
            "h2_line": self.h2_line,
            "direct": self.direct,
            "char_birth": self.char_birth,
            "char_exit": self.char_exit,
            "age_exit1": self.age_exit1,
            "age_exit2": self.age_exit2,
            "risk": risk,
            "risk_x": x,
            "risk_x2": x2,
            "exit_y": self.exit_y,
            "exit_xy": self.exit_xy,
        }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def git_blob(
    root: Path, commit: str, relative: str, expected_sha256: str
) -> Iterator[Path]:
    temporary = tempfile.NamedTemporaryFile(prefix="matching-p337-", delete=False)
    temporary_path = Path(temporary.name)
    digest = hashlib.sha256()
    try:
        process = subprocess.Popen(
            ["git", "cat-file", "blob", f"{commit}:{relative}"],
            cwd=root,
            stdout=subprocess.PIPE,
        )
        assert process.stdout is not None
        with temporary:
            for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
                digest.update(chunk)
                temporary.write(chunk)
        if process.wait() != 0:
            raise RuntimeError(f"git cat-file failed for {commit}:{relative}")
        observed = digest.hexdigest()
        if observed != expected_sha256:
            raise ValueError(
                f"SHA256 mismatch for {relative}: {observed} != {expected_sha256}"
            )
        yield temporary_path
    finally:
        temporary_path.unlink(missing_ok=True)


def open_text(path: Path, compression: str):
    if compression == "gzip":
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    if compression != "none":
        raise ValueError(f"unknown compression {compression!r}")
    return path.open(mode="rt", encoding="utf-8", newline="")


def read_run(root: Path, run: Mapping[str, Any]) -> dict[str, Any]:
    n = int(run["N"])
    column_map = run.get("column_map", {})
    builders: dict[tuple[str, int], BatchBuilder] = {}
    with git_blob(root, run["archive_commit"], run["path"], run["sha256"]) as path:
        with open_text(path, run.get("compression", "none")) as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError("birth archive has no header")
            tau1_key = column_map.get("tau1", "tau1")
            tau2_key = column_map.get("tau2", "tau2")
            ell_x_key = column_map.get("ell_x", "ell_x")
            ell_y_key = column_map.get("ell_y", "ell_y")
            direct_key = column_map.get("direct")
            required = {
                "n", "orientation", "batch", "samples", tau1_key, tau2_key,
                ell_x_key, ell_y_key, "count",
            }
            if direct_key:
                required.add(direct_key)
            else:
                required.add("kind")
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"archive missing columns {sorted(missing)}")
            for raw in reader:
                if int(raw["n"]) != n:
                    raise ValueError("archive N disagrees with manifest")
                orientation = raw["orientation"]
                if orientation not in ORIENTATIONS:
                    raise ValueError(f"unexpected orientation {orientation!r}")
                batch = int(raw["batch"])
                key = (orientation, batch)
                builder = builders.get(key)
                if builder is None:
                    builder = BatchBuilder(n)
                    builders[key] = builder
                direct = (
                    bool(int(raw[direct_key]))
                    if direct_key else raw["kind"] == "DIRECT_RANK2"
                )
                builder.add(
                    samples=int(raw["samples"]),
                    k1=int(raw[tau1_key]),
                    k2=int(raw[tau2_key]),
                    direct=direct,
                    ell_x=int(raw[ell_x_key]),
                    ell_y=int(raw[ell_y_key]),
                    count=int(raw["count"]),
                )
    finished = {key: builder.finish() for key, builder in builders.items()}
    batch_sets = {
        orientation: {batch for found, batch in finished if found == orientation}
        for orientation in ORIENTATIONS
    }
    if batch_sets["first"] != batch_sets["second"]:
        raise ValueError("orientation batch IDs are not aligned")
    batch_ids = sorted(batch_sets["first"])
    if len(batch_ids) != int(run["batches"]):
        raise ValueError("batch count disagrees with manifest")
    for key, row in finished.items():
        if row["samples"] * len(batch_ids) != int(run["samples_per_orientation"]):
            raise ValueError(f"samples per orientation disagree at {key}")
    with git_blob(
        root,
        run["archive_commit"],
        run["metadata_path"],
        run["metadata_sha256"],
    ) as metadata_path:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    return {
        "batches": finished,
        "batch_ids": batch_ids,
        "metadata_git_commit": metadata.get("git_commit"),
    }


ARRAY_KEYS = (
    "h1_line", "h2_line", "direct", "age_exit1", "age_exit2",
)
NESTED_KEYS = (
    "char_birth", "char_exit", "risk", "risk_x", "risk_x2", "exit_y", "exit_xy",
)


def add_in_place(target: list[float], source: Sequence[float]) -> None:
    for i, value in enumerate(source):
        target[i] += value


def totals(rows: Sequence[Mapping[str, Any]], n: int) -> dict[str, Any]:
    output: dict[str, Any] = {key: zeros(n) for key in ARRAY_KEYS}
    output.update({key: nested_zeros(len(rows[0][key]), n) for key in NESTED_KEYS})
    output["samples"] = 0
    output["rows"] = 0
    for row in rows:
        output["samples"] += row["samples"]
        output["rows"] += row["rows"]
        for key in ARRAY_KEYS:
            add_in_place(output[key], row[key])
        for key in NESTED_KEYS:
            for target_part, source_part in zip(output[key], row[key]):
                add_in_place(target_part, source_part)
    return output


def combined(
    total: Mapping[str, Any], omitted: Mapping[str, Any] | None, key: str
) -> list[float] | list[list[float]]:
    if key in ARRAY_KEYS:
        if omitted is None:
            return list(total[key])
        return [a - b for a, b in zip(total[key], omitted[key])]
    if omitted is None:
        return [list(row) for row in total[key]]
    return [
        [a - b for a, b in zip(total_row, omitted_row)]
        for total_row, omitted_row in zip(total[key], omitted[key])
    ]


def binomial_pmf(n: int, p: float) -> list[float]:
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie in (0,1)")
    values = [0.0] * (n + 1)
    mode = min(n, max(0, int(math.floor((n + 1) * p))))
    values[mode] = math.exp(
        math.lgamma(n + 1) - math.lgamma(mode + 1) - math.lgamma(n - mode + 1)
        + mode * math.log(p) + (n - mode) * math.log1p(-p)
    )
    for k in range(mode, 0, -1):
        values[k - 1] = values[k] * k / (n - k + 1) * (1.0 - p) / p
    for k in range(mode, n):
        values[k + 1] = values[k] * (n - k) / (k + 1) * p / (1.0 - p)
    norm = math.fsum(values)
    return [value / norm for value in values]


def tails(pmf: Sequence[float]) -> list[float]:
    output = [0.0] * len(pmf)
    running = 0.0
    for k in range(len(pmf) - 1, -1, -1):
        running += pmf[k]
        output[k] = running
    return output


def dot(left: Sequence[float], right: Sequence[float]) -> float:
    return math.fsum(a * b for a, b in zip(left, right))


def evaluate_orientation(
    total: Mapping[str, Any], omitted: Mapping[str, Any] | None, n: int, p: float
) -> dict[str, Any]:
    sample_count = total["samples"] - (omitted["samples"] if omitted else 0)
    h1_line = combined(total, omitted, "h1_line")
    h2_line = combined(total, omitted, "h2_line")
    direct = combined(total, omitted, "direct")
    char_birth = combined(total, omitted, "char_birth")
    char_exit = combined(total, omitted, "char_exit")
    age_exit1 = combined(total, omitted, "age_exit1")
    age_exit2 = combined(total, omitted, "age_exit2")
    assert isinstance(h1_line[0], float) and isinstance(char_birth[0], list)
    pmf_n = binomial_pmf(n, p)
    tail = tails(pmf_n)
    phi = [0.0] + [n * value for value in binomial_pmf(n - 1, p)]
    inv_samples = 1.0 / sample_count
    h1_all = [(h1_line[k] + direct[k]) * inv_samples for k in range(n + 1)]
    h2_all = [(h2_line[k] + direct[k]) * inv_samples for k in range(n + 1)]
    line_birth = [value * inv_samples for value in h1_line]
    line_exit = [value * inv_samples for value in h2_line]
    direct_mass = [value * inv_samples for value in direct]
    f1 = dot(h1_all, tail)
    f2 = dot(h2_all, tail)
    j01 = dot(line_birth, phi)
    j12 = dot(line_exit, phi)
    j02 = dot(direct_mass, phi)
    scale = p * (1.0 - p)
    character: dict[str, Any] = {}
    for j, name in enumerate(CHARACTERS):
        birth_hist = [value * inv_samples for value in char_birth[j]]
        exit_hist = [value * inv_samples for value in char_exit[j]]
        birth = dot(birth_hist, phi)
        exit_ = dot(exit_hist, phi)
        plateau = dot(birth_hist, tail) - dot(exit_hist, tail)
        character[name] = {
            "plateau": plateau,
            "birth": birth,
            "completion": exit_,
            "K": scale * (birth - exit_) / plateau if abs(plateau) > 1e-15 else None,
        }
    h_birth_hist = [value * inv_samples for value in char_birth[0]]
    h_exit_hist = [value * inv_samples for value in char_exit[0]]
    axis_birth_hist = [0.5 * line_birth[k] + h_birth_hist[k] for k in range(n + 1)]
    axis_exit_hist = [0.5 * line_exit[k] + h_exit_hist[k] for k in range(n + 1)]
    diagonal_birth_hist = [0.5 * line_birth[k] - h_birth_hist[k] for k in range(n + 1)]
    diagonal_exit_hist = [0.5 * line_exit[k] - h_exit_hist[k] for k in range(n + 1)]
    activity = {}
    for name, birth_hist, exit_hist in (
        ("A", axis_birth_hist, axis_exit_hist),
        ("D", diagonal_birth_hist, diagonal_exit_hist),
    ):
        birth = dot(birth_hist, phi)
        exit_ = dot(exit_hist, phi)
        plateau = dot(birth_hist, tail) - dot(exit_hist, tail)
        activity[name] = {
            "plateau": plateau,
            "birth": birth,
            "completion": exit_,
            "K": scale * (birth - exit_) / plateau if abs(plateau) > 1e-15 else None,
        }

    risk = combined(total, omitted, "risk")
    risk_x = combined(total, omitted, "risk_x")
    risk_x2 = combined(total, omitted, "risk_x2")
    exit_y = combined(total, omitted, "exit_y")
    exit_xy = combined(total, omitted, "exit_xy")
    assert isinstance(risk[0], list)
    numerator = 0.0
    denominator = 0.0
    for line in range(len(LINE_ORDER)):
        for k, weight in enumerate(pmf_n):
            r = risk[line][k]
            if r <= 0.0:
                continue
            numerator += weight * (
                exit_xy[line][k] - risk_x[line][k] * exit_y[line][k] / r
            )
            denominator += weight * (
                risk_x2[line][k] - risk_x[line][k] * risk_x[line][k] / r
            )
    age_beta = numerator / denominator if denominator > 0.0 else None
    f1_prime = dot(h1_all, phi)
    f2_prime = dot(h2_all, phi)
    closure = {
        "F1_prime_minus_J01_J02": f1_prime - j01 - j02,
        "F2_prime_minus_J12_J02": f2_prime - j12 - j02,
        "E_prime_direct_cancellation": (f2_prime - f1_prime) - (j12 - j01),
        "A_prime_state_current_closure": (f1_prime + f2_prime) - (j01 + j12 + 2.0 * j02),
    }
    return {
        "p": p,
        "F1": f1,
        "F2": f2,
        "P0": 1.0 - f1,
        "P1": f1 - f2,
        "P2": f2,
        "E_top": 1.0 - f1 + f2,
        "A_top": f1 + f2 - 1.0,
        "J01": j01,
        "J12": j12,
        "J02": j02,
        "E_log_derivative": scale * (j12 - j01),
        "A_log_derivative": scale * (j01 + j12 + 2.0 * j02),
        "J01_log": scale * j01,
        "J12_log": scale * j12,
        "J02_log": scale * j02,
        "collision_mass": math.fsum(direct_mass),
        "completion_age1_log": scale * dot(
            [value * inv_samples for value in age_exit1], phi
        ),
        "completion_age2_log": scale * dot(
            [value * inv_samples for value in age_exit2], phi
        ),
        "age_hazard_beta": age_beta,
        "character": character,
        "natural_activity": activity,
        "closure": closure,
    }


def matching_root(
    totals_by_orientation: Mapping[str, Mapping[str, Any]],
    omitted_by_orientation: Mapping[str, Mapping[str, Any] | None],
    n: int,
) -> float:
    def orientation_a(
        total: Mapping[str, Any], omitted: Mapping[str, Any] | None, p: float
    ) -> float:
        sample_count = total["samples"] - (omitted["samples"] if omitted else 0)
        h1_line = combined(total, omitted, "h1_line")
        h2_line = combined(total, omitted, "h2_line")
        direct = combined(total, omitted, "direct")
        assert isinstance(h1_line[0], float)
        tail = tails(binomial_pmf(n, p))
        f1_plus_f2 = math.fsum(
            (h1_line[k] + h2_line[k] + 2.0 * direct[k]) * tail[k]
            for k in range(n + 1)
        ) / sample_count
        return f1_plus_f2 - 1.0

    def pooled_a(p: float) -> float:
        return 0.5 * math.fsum(
            orientation_a(totals_by_orientation[o], omitted_by_orientation[o], p)
            for o in ORIENTATIONS
        )

    lo, hi = 1e-9, 1.0 - 1e-9
    flo, fhi = pooled_a(lo), pooled_a(hi)
    if not flo < 0.0 < fhi:
        raise ValueError(f"pooled matching root is not bracketed: {flo}, {fhi}")
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        if pooled_a(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def angular_vector(
    first: Mapping[str, Any], second: Mapping[str, Any], delta_cos4: float
) -> list[float]:
    def contrast(key: str) -> float:
        return (second[key] - first[key]) / delta_cos4

    values = {
        "angular_E_top": contrast("E_top"),
        "angular_A_top": contrast("A_top"),
        "angular_E_log_derivative": contrast("E_log_derivative"),
        "angular_A_log_derivative": contrast("A_log_derivative"),
        "angular_J01_log": contrast("J01_log"),
        "angular_J12_log": contrast("J12_log"),
        "angular_J02_log": contrast("J02_log"),
        "angular_collision_mass": contrast("collision_mass"),
        "angular_completion_age1_log": contrast("completion_age1_log"),
        "angular_completion_age2_log": contrast("completion_age2_log"),
        "angular_age_hazard_beta": contrast("age_hazard_beta"),
    }
    for name in CHARACTERS:
        first_character = first["character"][name]
        second_character = second["character"][name]
        values[f"angular_{name}_birth_log"] = (
            second["p"] * (1.0 - second["p"])
            * (second_character["birth"] - first_character["birth"])
            / delta_cos4
        )
        values[f"angular_{name}_completion_log"] = (
            second["p"] * (1.0 - second["p"])
            * (second_character["completion"] - first_character["completion"])
            / delta_cos4
        )
    values["angular_K_H_linear"] = (
        second["character"]["H"]["K"] - first["character"]["H"]["K"]
    ) / delta_cos4
    for name in ("A", "D"):
        values[f"angular_K_{name}_activity"] = (
            second["natural_activity"][name]["K"]
            - first["natural_activity"][name]["K"]
        ) / delta_cos4
    return [float(values[key]) for key in VECTOR_ORDER]


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    mean = [math.fsum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    return [
        [
            (count - 1.0) / count
            * math.fsum((row[i] - mean[i]) * (row[j] - mean[j]) for row in rows)
            for j in range(len(mean))
        ]
        for i in range(len(mean))
    ]


def score_context(
    *,
    run: Mapping[str, Any],
    batches: Mapping[tuple[str, int], Mapping[str, Any]],
    batch_ids: Sequence[int],
    context: Mapping[str, Any],
) -> dict[str, Any]:
    n = int(run["N"])
    rows_by_orientation = {
        orientation: [batches[(orientation, batch)] for batch in batch_ids]
        for orientation in ORIENTATIONS
    }
    totals_by_orientation = {
        orientation: totals(rows_by_orientation[orientation], n)
        for orientation in ORIENTATIONS
    }
    delta_cos4 = float(
        Fraction(str(run["cos4"]["second"])) - Fraction(str(run["cos4"]["first"]))
    )

    def evaluate(omitted_index: int | None) -> tuple[float, dict[str, Any], list[float]]:
        omitted = {
            orientation: (
                rows_by_orientation[orientation][omitted_index]
                if omitted_index is not None else None
            )
            for orientation in ORIENTATIONS
        }
        p = (
            matching_root(totals_by_orientation, omitted, n)
            if context["p"] == "pooled_matching_root"
            else float(context["p"])
        )
        orientation_values = {
            orientation: evaluate_orientation(
                totals_by_orientation[orientation], omitted[orientation], n, p
            )
            for orientation in ORIENTATIONS
        }
        vector = angular_vector(
            orientation_values["first"], orientation_values["second"], delta_cos4
        )
        return p, orientation_values, vector

    point_p, orientation_values, point = evaluate(None)
    leave = [evaluate(index) for index in range(len(batch_ids))]
    leave_vectors = [row[2] for row in leave]
    covariance = jackknife_covariance(leave_vectors)
    maximum_closure = max(
        abs(value)
        for orientation in orientation_values.values()
        for value in orientation["closure"].values()
    )
    return {
        "id": context["id"],
        "p": point_p,
        "delta_cos4_exact": f"{run['cos4']['second']}-({run['cos4']['first']})",
        "delta_cos4": delta_cos4,
        "orientation_values": orientation_values,
        "vector_order": list(VECTOR_ORDER),
        "vector": point,
        "standard_error": [math.sqrt(max(0.0, covariance[i][i])) for i in range(len(point))],
        "covariance": covariance,
        "delete_one": {
            "unit": "same_batch_removed_from_first_and_second_orientations",
            "batch_ids": list(batch_ids),
            "p_values": [row[0] for row in leave],
            "vectors": leave_vectors,
        },
        "maximum_state_current_closure_residual": maximum_closure,
    }


def score(manifest: Mapping[str, Any], root: Path, workers: int) -> dict[str, Any]:
    started = time.perf_counter()
    outputs = []
    for run in manifest["runs"]:
        loaded = read_run(root, run)
        contexts = [
            score_context(
                run=run,
                batches=loaded["batches"],
                batch_ids=loaded["batch_ids"],
                context=context,
            )
            for context in manifest["evaluation_contexts"]
        ]
        direct_by_orientation = {
            orientation: int(
                math.fsum(
                    math.fsum(loaded["batches"][(orientation, batch)]["direct"])
                    for batch in loaded["batch_ids"]
                )
            )
            for orientation in ORIENTATIONS
        }
        outputs.append(
            {
                "id": run["id"],
                "role": run["role"],
                "N": run["N"],
                "dependency_group": run["dependency_group"],
                "source": {
                    "archive_commit": run["archive_commit"],
                    "path": run["path"],
                    "sha256_verified": run["sha256"],
                    "metadata_path": run["metadata_path"],
                    "metadata_sha256_verified": run["metadata_sha256"],
                    "metadata_git_commit": loaded["metadata_git_commit"],
                    "batches": len(loaded["batch_ids"]),
                    "samples_per_orientation": run["samples_per_orientation"],
                    "direct_counts": direct_by_orientation,
                },
                "contexts": contexts,
            }
        )
    primary = [row for row in outputs if row["role"] == "four_generation_primary"]
    sequence = {}
    for coordinate in (
        "angular_E_log_derivative", "angular_J01_log", "angular_J12_log",
        "angular_J02_log", "angular_age_hazard_beta", "angular_K_A_activity",
    ):
        index = VECTOR_ORDER.index(coordinate)
        sequence[coordinate] = [
            {"N": row["N"], "value": row["contexts"][0]["vector"][index]}
            for row in primary
        ]
    return {
        "schema": SCHEMA,
        "status": "completed_retrospective_existing_data_reanalysis",
        "new_samples": False,
        "workers_requested": workers,
        "manifest_schema": manifest["schema"],
        "semantic_contract": manifest["semantic_contract"],
        "covariance_contract": manifest["covariance"],
        "not_scoreable": manifest["not_scoreable"],
        "claim_boundary": manifest["claim_boundary"],
        "runs": outputs,
        "four_generation_sequences": sequence,
        "elapsed_seconds": time.perf_counter() - started,
    }


def markdown(result: Mapping[str, Any]) -> str:
    def cell(context: Mapping[str, Any], coordinate: str) -> str:
        index = VECTOR_ORDER.index(coordinate)
        value = context["vector"][index]
        error = context["standard_error"][index]
        return f"{value:.6g} +/- {error:.3g}"

    lines = [
        "# P337 birth-state current prism",
        "",
        "This is a zero-new-sample reconstruction from immutable joint birth archives. "
        "Each size is one dependency block; directions share one delete-one batch unit.",
        "",
        "## Pooled-root state-current decomposition",
        "",
        "| run | N | p_bar | E-dot | J01 | J12 | J02 | age1 completion | age beta | K_A activity |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in result["runs"]:
        context = run["contexts"][0]
        lines.append(
            "| {id} | {N} | {p:.9f} | {e} | {j01} | {j12} | {j02} | "
            "{age1} | {age} | {ka} |".format(
                id=run["id"], N=run["N"], p=context["p"],
                e=cell(context, "angular_E_log_derivative"),
                j01=cell(context, "angular_J01_log"),
                j12=cell(context, "angular_J12_log"),
                j02=cell(context, "angular_J02_log"),
                age1=cell(context, "angular_completion_age1_log"),
                age=cell(context, "angular_age_hazard_beta"),
                ka=cell(context, "angular_K_A_activity"),
            )
        )
    lines.extend(
        [
            "",
            "## Exact interpretation",
            "",
            "For line paths, `J01` is first birth and `J12` is second completion; "
            "`J02` is the line-free direct `0->2` current. The scorer verifies",
            "",
            "```text",
            "F1' = J01 + J02",
            "F2' = J12 + J02",
            "E_top' = J12 - J01",
            "A_top' = J01 + J12 + 2 J02.",
            "```",
            "",
            "Thus collision cancels from fixed-p E exactly. Its mass and current remain "
            "correlated coordinates, not a third projective-line component.",
            "",
            "Across the four primary sizes the directional collision mass and `J02` current "
            "are unresolved. At N170 the resolved E response is assembled by a negative "
            "first-birth contrast and a positive second-completion contrast, so the two line "
            "currents reinforce rather than cancel. At N340 the resolved piece is concentrated "
            "in second completion; N85 and N680 do not individually resolve both pieces.",
            "",
            "The completion-age first moment is directionally nonzero in all four primary "
            "blocks, while the conditional age-hazard coefficient is unresolved in this coarser "
            "archive. This distinguishes a robust age-weighted completion current from evidence "
            "that age itself causes the directional response. The richer N325/N425 P334 analysis "
            "separately rejects coarse-state age independence.",
            "",
            "The age coefficient is a line/current-layer fixed-effect diagnostic. A nonzero "
            "value rejects sufficiency of the retained coarse state; it cannot distinguish "
            "intrinsic time memory from unrecorded current geometry.",
            "",
            "## Dependency and boundary",
            "",
            "N85/N170/N340/N680 are independent seed blocks. N130 is a separate cross-lineage "
            "control. These outputs do not share sample covariance with the ten-size E_top or "
            "Euler archives. Completion winding, complement line, transporter, ambiguity, "
            "microscopic state and path order are not scoreable.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=root / "analysis/p337_birth_state_current_manifest.yaml",
    )
    parser.add_argument(
        "--output-json", type=Path,
        default=root / "results/p337-birth-state-current/latest.json",
    )
    parser.add_argument(
        "--output-md", type=Path,
        default=root / "results/p337-birth-state-current/latest.md",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--git-root", type=Path, default=root)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    result = score(manifest, args.git_root, args.workers)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(markdown(result), encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": result["schema"],
                "runs": len(result["runs"]),
                "elapsed_seconds": result["elapsed_seconds"],
                "output_json": str(args.output_json),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
