#!/usr/bin/env python3
"""Decompose orientation H4 into the two exact homology activations.

The threshold-rank convention is used without changing the archived data:

    K1 = K_minus,  K2 = K_plus,
    F_i(p) = E[Pr(Binomial(N,p) >= K_i)],
    M(p) = -1 + F1(p) + F2(p).

For every same-N orientation pair, the analyzer solves the root of the pooled
matching curve and evaluates the two orientation differences at that common
coordinate.  The pooled root, both orientation roots, all nonlinear
combinations, and all joint-moment coordinates are recomputed inside aligned
delete-one-batch replicates.  Archives sharing the same counter stream are
retained in one cross-size covariance block.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Optional, Sequence

import yaml
import mpmath as mp


SCHEMA = "matching-one.two-activation-h4.v1"
MANIFEST_SCHEMA = "matching-one.two-activation-h4.manifest.v1"
MOMENT_FIELDS = (
    "sum_kminus",
    "sum_kplus",
    "sum_kminus2",
    "sum_kplus2",
    "sum_product",
    "sum_gap",
    "sum_gap2",
)
MOMENT_METRICS = (
    "C_mean",
    "C_variance",
    "G_mean",
    "G_variance",
    "C_G_covariance",
)
BASE_METRICS = (
    "p_bar",
    "first_F1",
    "first_F2",
    "second_F1",
    "second_F2",
    "delta_F1",
    "delta_F2",
    "delta_M",
    "pooled_M_prime",
    "delta_p1",
    "delta_p2",
    "linearized_root_gap",
    "actual_root_gap",
    "nonlinear_closure_residual",
    "angular_delta_F1",
    "angular_delta_F2",
    "angular_delta_M",
    "angular_delta_p1",
    "angular_delta_p2",
    "angular_linearized_root_gap",
    "angular_actual_root_gap",
    "angular_nonlinear_closure_residual",
)
FULL_METRICS = BASE_METRICS + tuple(
    f"{scope}_{metric}"
    for scope in ("first", "second", "pooled")
    for metric in MOMENT_METRICS
)
DECISION_METRICS = (
    "angular_delta_F1",
    "angular_delta_F2",
    "angular_delta_M",
    "angular_delta_p1",
    "angular_delta_p2",
    "angular_actual_root_gap",
    "angular_nonlinear_closure_residual",
)
GIT_INPUT_RE = re.compile(r"^git\+([0-9a-fA-F]{7,64}):(.+)$")


class ArchiveNotScoreable(ValueError):
    """The archive cannot support the frozen two-activation score."""


@dataclass(frozen=True)
class HistogramBatch:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    k1: tuple[int, ...]
    k2: tuple[int, ...]


@dataclass(frozen=True)
class MomentBatch:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    totals: Mapping[str, int]


@dataclass(frozen=True)
class Archive:
    n: int
    dependency_group: str
    histograms: Mapping[str, tuple[HistogramBatch, ...]]
    moments: Mapping[str, tuple[MomentBatch, ...]]
    metadata: Mapping[str, Any]
    paths: Mapping[str, Path]
    input_provenance: Mapping[str, Mapping[str, str]]


@dataclass(frozen=True)
class MaterializedInput:
    path: Path
    provenance: Mapping[str, str]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve(root: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else root / path


def _git_text(root: Path, arguments: Sequence[str]) -> str:
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), *arguments],
            check=True,
            capture_output=True,
            text=True,
            shell=False,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        raise ArchiveNotScoreable(f"git input lookup failed: {detail.strip()}") from exc
    return completed.stdout.strip()


def _materialize_input(root: Path, value: object) -> MaterializedInput:
    """Resolve a local path or materialize one immutable git blob in ignored tmp/."""

    source = str(value)
    if not source.startswith("git+"):
        path = _resolve(root, source)
        display = str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
        return MaterializedInput(path=path, provenance={"path": display})

    match = GIT_INPUT_RE.fullmatch(source)
    if match is None:
        raise ArchiveNotScoreable(
            "git inputs must use git+<hex-commit>:<repository-relative-path>"
        )
    revision, repository_path_text = match.groups()
    repository_path = PurePosixPath(repository_path_text)
    if (
        repository_path.is_absolute()
        or not repository_path.parts
        or any(part in ("", ".", "..") for part in repository_path.parts)
    ):
        raise ArchiveNotScoreable("git input path must be a normalized repository-relative path")

    resolved_commit = _git_text(root, ["rev-parse", "--verify", f"{revision}^{{commit}}"])
    if not re.fullmatch(r"[0-9a-fA-F]{40,64}", resolved_commit):
        raise ArchiveNotScoreable("git input did not resolve to a full commit object id")
    object_spec = f"{resolved_commit}:{repository_path.as_posix()}"
    blob_oid = _git_text(root, ["rev-parse", "--verify", object_spec])
    if not re.fullmatch(r"[0-9a-fA-F]{40,64}", blob_oid):
        raise ArchiveNotScoreable("git input did not resolve to a blob object id")

    cache_directory = root / "tmp" / "git-blobs" / resolved_commit
    cache_directory.mkdir(parents=True, exist_ok=True)
    suffix = Path(repository_path.name).suffix
    cache_path = cache_directory / f"{blob_oid}{suffix}"
    if not cache_path.is_file():
        try:
            completed = subprocess.run(
                ["git", "-C", str(root), "cat-file", "blob", blob_oid],
                check=True,
                capture_output=True,
                shell=False,
            )
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            detail = getattr(exc, "stderr", b"")
            if isinstance(detail, bytes):
                detail = detail.decode("utf-8", errors="replace")
            raise ArchiveNotScoreable(f"git blob materialization failed: {str(detail).strip()}") from exc
        temporary_name: Optional[str] = None
        try:
            with tempfile.NamedTemporaryFile(
                dir=cache_directory, prefix=".materializing-", delete=False
            ) as handle:
                handle.write(completed.stdout)
                temporary_name = handle.name
            os.replace(temporary_name, cache_path)
            temporary_name = None
        finally:
            if temporary_name is not None:
                Path(temporary_name).unlink(missing_ok=True)

    cached_oid = _git_text(root, ["hash-object", str(cache_path)])
    if cached_oid != blob_oid:
        raise ArchiveNotScoreable("cached git blob hash differs from the pinned object id")
    return MaterializedInput(
        path=cache_path,
        provenance={
            "path": source,
            "uri": source,
            "artifact_commit": resolved_commit,
            "git_blob_oid": blob_oid,
            "materialized_path": str(cache_path.relative_to(root)),
            "cache_policy": "ignored_tmp/git-blobs",
        },
    )


def load_manifest(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"manifest must use schema {MANIFEST_SCHEMA}")
    if payload.get("status") != "retrospective_existing_data_reanalysis":
        raise ValueError("manifest status must remain retrospective_existing_data_reanalysis")
    contract = payload.get("semantic_contract")
    if not isinstance(contract, dict):
        raise ValueError("manifest lacks semantic_contract")
    if contract.get("K1") != "K_minus" or contract.get("K2") != "K_plus":
        raise ValueError("the exact activation mapping K1=K_minus, K2=K_plus changed")
    required = [int(value) for value in payload.get("required_sizes", ())]
    runs = payload.get("runs")
    if not required or len(required) != len(set(required)):
        raise ValueError("required_sizes must be a nonempty unique list")
    if not isinstance(runs, list) or [int(row["N"]) for row in runs] != required:
        raise ValueError("runs must follow required_sizes exactly")
    return payload


def _metadata_contract(metadata: Mapping[str, Any], entry: Mapping[str, Any]) -> None:
    n = int(entry["N"])
    designs = metadata.get("designs")
    if not isinstance(designs, list) or len(designs) != 1:
        raise ArchiveNotScoreable("metadata must contain exactly one Gaussian design")
    design = designs[0]
    if int(design.get("N", -1)) != n:
        raise ArchiveNotScoreable("metadata N differs from manifest N")
    expected_first = [int(value) for value in entry["expected_first"]]
    expected_second = [int(value) for value in entry["expected_second"]]
    if design.get("first") != expected_first or design.get("second") != expected_second:
        raise ArchiveNotScoreable("orientation representatives differ from the manifest")
    if metadata.get("per_batch_joint_moments") is not True:
        raise ArchiveNotScoreable("archive lacks per-batch joint moments")
    coupling = str(metadata.get("coupling", "")).lower()
    if "permutation" not in coupling or "shared" not in coupling:
        raise ArchiveNotScoreable("orientations do not declare a shared same-N random stream")
    if "first black primal cross rank" not in str(metadata.get("K_plus", "")):
        raise ArchiveNotScoreable("K_plus semantics are not the frozen threshold convention")
    if "white matching cross is lost" not in str(metadata.get("K_minus", "")):
        raise ArchiveNotScoreable("K_minus semantics are not the frozen threshold convention")
    actual_group = "crn-{}-{}-{}".format(
        int(metadata["seed"]),
        int(metadata["replica_counter_first"]),
        int(metadata["replica_counter_last_exclusive"]),
    )
    if entry.get("dependency_group") != actual_group:
        raise ArchiveNotScoreable(
            f"dependency_group {entry.get('dependency_group')!r} != {actual_group!r}"
        )


def _read_histograms(path: Path, n: int) -> dict[str, tuple[HistogramBatch, ...]]:
    required = {"n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"}
    mutable: dict[tuple[str, int], dict[str, Any]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ArchiveNotScoreable("histogram CSV missing " + ", ".join(sorted(missing)))
        for raw in reader:
            row_n = int(raw["n"])
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            samples = int(raw["samples"])
            kind = raw["kind"]
            rank = int(raw["k"])
            count = int(raw["count"])
            if row_n != n or orientation not in ("first", "second"):
                raise ArchiveNotScoreable("histogram contains an unexpected N/orientation")
            if kind not in ("minus", "plus") or not 1 <= rank <= n or count <= 0:
                raise ArchiveNotScoreable("histogram contains an invalid rank/count")
            key = (orientation, batch)
            item = mutable.setdefault(
                key,
                {
                    "a": int(raw["a"]),
                    "b": int(raw["b"]),
                    "samples": samples,
                    "minus": [0] * (n + 1),
                    "plus": [0] * (n + 1),
                },
            )
            if (item["a"], item["b"], item["samples"]) != (
                int(raw["a"]), int(raw["b"]), samples
            ):
                raise ArchiveNotScoreable("histogram batch metadata are inconsistent")
            item[kind][rank] += count
    output: dict[str, tuple[HistogramBatch, ...]] = {}
    for orientation in ("first", "second"):
        rows = []
        for (row_orientation, batch), item in sorted(mutable.items(), key=lambda pair: pair[0][1]):
            if row_orientation != orientation:
                continue
            if sum(item["minus"]) != item["samples"] or sum(item["plus"]) != item["samples"]:
                raise ArchiveNotScoreable("marginal histogram total differs from batch samples")
            rows.append(
                HistogramBatch(
                    n=n,
                    a=item["a"],
                    b=item["b"],
                    orientation=orientation,
                    batch=batch,
                    samples=item["samples"],
                    k1=tuple(item["minus"]),
                    k2=tuple(item["plus"]),
                )
            )
        output[orientation] = tuple(rows)
    return output


def _read_moments(path: Path, n: int) -> dict[str, tuple[MomentBatch, ...]]:
    required = {"n", "a", "b", "orientation", "batch", "samples", *MOMENT_FIELDS}
    output: dict[str, list[MomentBatch]] = {"first": [], "second": []}
    seen: set[tuple[str, int]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ArchiveNotScoreable("moments CSV missing " + ", ".join(sorted(missing)))
        for raw in reader:
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            key = (orientation, batch)
            if int(raw["n"]) != n or orientation not in output or key in seen:
                raise ArchiveNotScoreable("moments contain an unexpected or duplicate row")
            seen.add(key)
            totals = {name: int(raw[name]) for name in MOMENT_FIELDS}
            if totals["sum_gap"] != totals["sum_kplus"] - totals["sum_kminus"]:
                raise ArchiveNotScoreable("first gap-moment identity failed")
            if totals["sum_gap2"] != (
                totals["sum_kplus2"] + totals["sum_kminus2"] - 2 * totals["sum_product"]
            ):
                raise ArchiveNotScoreable("second gap-moment identity failed")
            output[orientation].append(
                MomentBatch(
                    n=n,
                    a=int(raw["a"]),
                    b=int(raw["b"]),
                    orientation=orientation,
                    batch=batch,
                    samples=int(raw["samples"]),
                    totals=totals,
                )
            )
    return {orientation: tuple(sorted(rows, key=lambda row: row.batch)) for orientation, rows in output.items()}


def read_archive(root: Path, entry: Mapping[str, Any]) -> Archive:
    n = int(entry["N"])
    materialized = {
        kind: _materialize_input(root, entry[kind])
        for kind in ("histogram", "moments", "metadata")
    }
    paths = {kind: item.path for kind, item in materialized.items()}
    for kind, path in paths.items():
        if not path.is_file():
            raise ArchiveNotScoreable(f"{kind} input is missing: {path}")
    metadata = json.loads(paths["metadata"].read_text(encoding="utf-8"))
    if not isinstance(metadata, dict):
        raise ArchiveNotScoreable("metadata must be a JSON object")
    _metadata_contract(metadata, entry)
    histograms = _read_histograms(paths["histogram"], n)
    moments = _read_moments(paths["moments"], n)
    batches = int(metadata["batches"])
    signatures: dict[str, list[tuple[int, int]]] = {}
    for orientation in ("first", "second"):
        h_rows = histograms[orientation]
        m_rows = moments[orientation]
        if [row.batch for row in h_rows] != list(range(batches)):
            raise ArchiveNotScoreable(f"{orientation} histogram batches are incomplete")
        if [row.batch for row in m_rows] != list(range(batches)):
            raise ArchiveNotScoreable(f"{orientation} moment batches are incomplete")
        for hist, moment in zip(h_rows, m_rows):
            if (hist.a, hist.b, hist.samples) != (moment.a, moment.b, moment.samples):
                raise ArchiveNotScoreable("histogram/moment batch metadata disagree")
            calculated = {
                "sum_kminus": sum(rank * count for rank, count in enumerate(hist.k1)),
                "sum_kplus": sum(rank * count for rank, count in enumerate(hist.k2)),
                "sum_kminus2": sum(rank * rank * count for rank, count in enumerate(hist.k1)),
                "sum_kplus2": sum(rank * rank * count for rank, count in enumerate(hist.k2)),
            }
            if any(moment.totals[name] != value for name, value in calculated.items()):
                raise ArchiveNotScoreable("joint moments disagree with marginal histograms")
        if sum(row.samples for row in h_rows) != int(metadata["samples_per_pair"]):
            raise ArchiveNotScoreable("sample total differs from metadata")
        signatures[orientation] = [(row.batch, row.samples) for row in h_rows]
    if signatures["first"] != signatures["second"]:
        raise ArchiveNotScoreable("orientations are not aligned batch by batch")
    counter_count = int(metadata["replica_counter_last_exclusive"]) - int(
        metadata["replica_counter_first"]
    )
    if counter_count != int(metadata["samples_per_pair"]):
        raise ArchiveNotScoreable("counter interval does not equal samples_per_pair")
    return Archive(
        n=n,
        dependency_group=str(entry["dependency_group"]),
        histograms=histograms,
        moments=moments,
        metadata=metadata,
        paths=paths,
        input_provenance={
            kind: item.provenance for kind, item in materialized.items()
        },
    )


def _sum_histograms(
    rows: Sequence[HistogramBatch], omitted_batch: Optional[int] = None
) -> dict[str, Any]:
    selected = [row for row in rows if row.batch != omitted_batch]
    if len(selected) < 2:
        raise ArchiveNotScoreable("fewer than two batches remain after deletion")
    n = selected[0].n
    k1 = [0] * (n + 1)
    k2 = [0] * (n + 1)
    for row in selected:
        for rank in range(1, n + 1):
            k1[rank] += row.k1[rank]
            k2[rank] += row.k2[rank]
    return {
        "a": selected[0].a,
        "b": selected[0].b,
        "samples": sum(row.samples for row in selected),
        "k1": k1,
        "k2": k2,
    }


def _sum_moments(
    rows: Sequence[MomentBatch], omitted_batch: Optional[int] = None
) -> dict[str, int]:
    selected = [row for row in rows if row.batch != omitted_batch]
    return {
        "samples": sum(row.samples for row in selected),
        **{name: sum(row.totals[name] for row in selected) for name in MOMENT_FIELDS},
    }


def activation_components(
    n: int, samples: int, k1: Sequence[int], k2: Sequence[int], p: float
) -> tuple[float, float, float, float]:
    """Return F1,F2,F1',F2' using stable binomial recurrences."""

    if p <= 0.0:
        return 0.0, 0.0, 0.0, 0.0
    if p >= 1.0:
        return 1.0, 1.0, 0.0, 0.0
    q = 1.0 - p
    probability = q**n
    cumulative1 = 0
    cumulative2 = 0
    total1: list[float] = []
    total2: list[float] = []
    for occupied in range(n + 1):
        if occupied:
            cumulative1 += k1[occupied]
            cumulative2 += k2[occupied]
        total1.append(cumulative1 * probability)
        total2.append(cumulative2 * probability)
        if occupied < n:
            probability *= (n - occupied) * p / ((occupied + 1) * q)
    density = n * q ** (n - 1)
    derivative1: list[float] = []
    derivative2: list[float] = []
    for rank in range(1, n + 1):
        derivative1.append(k1[rank] * density)
        derivative2.append(k2[rank] * density)
        if rank < n:
            density *= (n - rank) * p / (rank * q)
    scale = float(samples)
    return (
        math.fsum(total1) / scale,
        math.fsum(total2) / scale,
        math.fsum(derivative1) / scale,
        math.fsum(derivative2) / scale,
    )


def matching_root(
    n: int,
    samples: int,
    k1: Sequence[int],
    k2: Sequence[int],
    initial: float = 0.59274605079210,
) -> float:
    lower = 0.25
    upper = 0.75
    lower_value = sum(activation_components(n, samples, k1, k2, lower)[:2]) - 1.0
    upper_value = sum(activation_components(n, samples, k1, k2, upper)[:2]) - 1.0
    if not lower_value < 0.0 < upper_value:
        raise ArchiveNotScoreable("matching curve does not bracket a unique physical root")
    p = min(max(initial, lower), upper)
    for _ in range(48):
        f1, f2, d1, d2 = activation_components(n, samples, k1, k2, p)
        value = f1 + f2 - 1.0
        derivative = d1 + d2
        if value < 0.0:
            lower = p
        else:
            upper = p
        if abs(value) <= 2.0e-15 or upper - lower <= 2.0e-15:
            break
        candidate = p - value / derivative if derivative > 0.0 else math.nan
        if not math.isfinite(candidate) or not lower < candidate < upper:
            candidate = (lower + upper) / 2.0
        p = candidate
    # The float Newton/bracketing loop is fast, but its last-ULP stopping point
    # can depend on the platform libm used by ``q**n``.  That is immaterial for
    # the point estimate yet can leak into a delete-one covariance through the
    # difference of two roots.  Finish with a fixed number of high-precision
    # Newton steps so the serialized root is the same correctly rounded float
    # on every supported platform.
    with mp.workdps(50):
        refined = mp.mpf(p)
        mp_samples = mp.mpf(samples)
        for _ in range(4):
            q = 1 - refined
            probability = q**n
            density = n * q ** (n - 1)
            cumulative1 = 0
            cumulative2 = 0
            value = mp.mpf(-1)
            derivative = mp.mpf(0)
            for occupied in range(n + 1):
                if occupied:
                    cumulative1 += k1[occupied]
                    cumulative2 += k2[occupied]
                value += (cumulative1 + cumulative2) * probability / mp_samples
                if occupied < n:
                    probability *= (
                        (n - occupied) * refined / ((occupied + 1) * q)
                    )
            for rank in range(1, n + 1):
                derivative += (k1[rank] + k2[rank]) * density / mp_samples
                if rank < n:
                    density *= (n - rank) * refined / (rank * q)
            refined -= value / derivative
    return float(refined)


def _cos4(a: int, b: int) -> float:
    n = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / float(n * n)


def _joint_coordinates(totals: Mapping[str, int]) -> dict[str, float]:
    samples = float(totals["samples"])
    mean1 = totals["sum_kminus"] / samples
    mean2 = totals["sum_kplus"] / samples
    var1 = totals["sum_kminus2"] / samples - mean1 * mean1
    var2 = totals["sum_kplus2"] / samples - mean2 * mean2
    cov12 = totals["sum_product"] / samples - mean1 * mean2
    c_mean = (mean1 + mean2) / 2.0
    g_mean = mean2 - mean1
    c_variance = (var1 + var2 + 2.0 * cov12) / 4.0
    g_variance = var1 + var2 - 2.0 * cov12
    c_g_covariance = (var2 - var1) / 2.0
    if min(var1, var2, c_variance, g_variance) < -1.0e-9:
        raise ArchiveNotScoreable("joint moments reconstruct a negative variance")
    return {
        "K1_mean": mean1,
        "K2_mean": mean2,
        "K1_variance": max(0.0, var1),
        "K2_variance": max(0.0, var2),
        "K1_K2_covariance": cov12,
        "C_mean": c_mean,
        "C_variance": max(0.0, c_variance),
        "G_mean": g_mean,
        "G_variance": max(0.0, g_variance),
        "C_G_covariance": c_g_covariance,
    }


def estimate(archive: Archive, omitted_batch: Optional[int] = None) -> dict[str, Any]:
    first = _sum_histograms(archive.histograms["first"], omitted_batch)
    second = _sum_histograms(archive.histograms["second"], omitted_batch)
    if first["samples"] != second["samples"]:
        raise ArchiveNotScoreable("orientation sample totals differ")
    pooled_k1 = [left + right for left, right in zip(first["k1"], second["k1"])]
    pooled_k2 = [left + right for left, right in zip(first["k2"], second["k2"])]
    pooled_samples = first["samples"] + second["samples"]
    p_bar = matching_root(archive.n, pooled_samples, pooled_k1, pooled_k2)
    first_values = activation_components(
        archive.n, first["samples"], first["k1"], first["k2"], p_bar
    )
    second_values = activation_components(
        archive.n, second["samples"], second["k1"], second["k2"], p_bar
    )
    first_root = matching_root(
        archive.n, first["samples"], first["k1"], first["k2"], p_bar
    )
    second_root = matching_root(
        archive.n, second["samples"], second["k1"], second["k2"], p_bar
    )
    first_f1, first_f2, first_d1, first_d2 = first_values
    second_f1, second_f2, second_d1, second_d2 = second_values
    delta_f1 = first_f1 - second_f1
    delta_f2 = first_f2 - second_f2
    delta_m = delta_f1 + delta_f2
    pooled_m_prime = (first_d1 + first_d2 + second_d1 + second_d2) / 2.0
    delta_p1 = -delta_f1 / pooled_m_prime
    delta_p2 = -delta_f2 / pooled_m_prime
    linearized = delta_p1 + delta_p2
    actual = first_root - second_root
    residual = actual - linearized
    delta_cos4 = _cos4(first["a"], first["b"]) - _cos4(second["a"], second["b"])
    if abs(delta_cos4) < 1.0e-15:
        raise ArchiveNotScoreable("orientation pair has zero Delta cos(4 theta)")
    first_moments = _sum_moments(archive.moments["first"], omitted_batch)
    second_moments = _sum_moments(archive.moments["second"], omitted_batch)
    pooled_moments = {
        "samples": first_moments["samples"] + second_moments["samples"],
        **{
            name: first_moments[name] + second_moments[name]
            for name in MOMENT_FIELDS
        },
    }
    joint = {
        "first": _joint_coordinates(first_moments),
        "second": _joint_coordinates(second_moments),
        "pooled": _joint_coordinates(pooled_moments),
    }
    metrics: dict[str, float] = {
        "p_bar": p_bar,
        "first_F1": first_f1,
        "first_F2": first_f2,
        "second_F1": second_f1,
        "second_F2": second_f2,
        "delta_F1": delta_f1,
        "delta_F2": delta_f2,
        "delta_M": delta_m,
        "pooled_M_prime": pooled_m_prime,
        "delta_p1": delta_p1,
        "delta_p2": delta_p2,
        "linearized_root_gap": linearized,
        "actual_root_gap": actual,
        "nonlinear_closure_residual": residual,
        "angular_delta_F1": delta_f1 / delta_cos4,
        "angular_delta_F2": delta_f2 / delta_cos4,
        "angular_delta_M": delta_m / delta_cos4,
        "angular_delta_p1": delta_p1 / delta_cos4,
        "angular_delta_p2": delta_p2 / delta_cos4,
        "angular_linearized_root_gap": linearized / delta_cos4,
        "angular_actual_root_gap": actual / delta_cos4,
        "angular_nonlinear_closure_residual": residual / delta_cos4,
    }
    for scope, values in joint.items():
        for metric in MOMENT_METRICS:
            metrics[f"{scope}_{metric}"] = values[metric]
    pooled_at_root = sum(activation_components(
        archive.n, pooled_samples, pooled_k1, pooled_k2, p_bar
    )[:2]) - 1.0
    reconstruction_errors = [
        abs((first_f1 + first_f2 - 1.0) - (first_f1 + first_f2 - 1.0)),
        abs((second_f1 + second_f2 - 1.0) - (second_f1 + second_f2 - 1.0)),
        abs(delta_m - (delta_f1 + delta_f2)),
        abs(linearized - (delta_p1 + delta_p2)),
        abs(pooled_at_root),
    ]
    if max(reconstruction_errors) > 2.0e-12:
        raise ArithmeticError("two-activation reconstruction identity failed")
    return {
        "metrics": metrics,
        "delta_cos4": delta_cos4,
        "first_root": first_root,
        "second_root": second_root,
        "joint": joint,
        "audit": {
            "pooled_M_at_p_bar": pooled_at_root,
            "delta_M_minus_delta_F1_plus_delta_F2": delta_m - delta_f1 - delta_f2,
            "linearized_minus_component_sum": linearized - delta_p1 - delta_p2,
        },
    }


def jackknife_covariance(
    left: Sequence[float], right: Sequence[float]
) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("delete-one vectors must be aligned and nontrivial")
    mean_left = math.fsum(left) / len(left)
    mean_right = math.fsum(right) / len(right)
    return (len(left) - 1.0) / len(left) * math.fsum(
        (x - mean_left) * (y - mean_right) for x, y in zip(left, right)
    )


def _qualitative(
    metrics: Mapping[str, float], standard_errors: Mapping[str, float]
) -> dict[str, Any]:
    first = metrics["angular_delta_p1"]
    second = metrics["angular_delta_p2"]
    denominator = abs(first) + abs(second)
    first_fraction = abs(first) / denominator if denominator else 0.5
    if first_fraction >= 2.0 / 3.0:
        dominant = "first_activation_K1"
    elif first_fraction <= 1.0 / 3.0:
        dominant = "second_activation_K2"
    else:
        dominant = "shared"
    if first * second < 0.0:
        interaction = "cancelling"
    elif first * second > 0.0:
        interaction = "reinforcing"
    else:
        interaction = "one_component_zero"
    cancellation = 1.0 - abs(first + second) / denominator if denominator else 0.0
    return {
        "descriptive_dominant_component": dominant,
        "absolute_K1_fraction": first_fraction,
        "component_interaction": interaction,
        "cancellation_fraction": cancellation,
        "angular_delta_p1_z": (
            first / standard_errors["angular_delta_p1"]
            if standard_errors["angular_delta_p1"] > 0.0 else math.nan
        ),
        "angular_delta_p2_z": (
            second / standard_errors["angular_delta_p2"]
            if standard_errors["angular_delta_p2"] > 0.0 else math.nan
        ),
        "guard": "descriptive per-size classification; not a cross-size independent-evidence score",
    }


def analyze_archive(root_text: str, entry: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(root_text)
    n = int(entry["N"])
    try:
        archive = read_archive(root, entry)
        full = estimate(archive)
        batch_ids = [row.batch for row in archive.histograms["first"]]
        deleted = [estimate(archive, batch) for batch in batch_ids]
        covariance = [
            [
                jackknife_covariance(
                    [row["metrics"][left] for row in deleted],
                    [row["metrics"][right] for row in deleted],
                )
                for right in FULL_METRICS
            ]
            for left in FULL_METRICS
        ]
        standard_errors = {
            name: math.sqrt(max(0.0, covariance[index][index]))
            for index, name in enumerate(FULL_METRICS)
        }
        metrics = full["metrics"]
        result = {
            "N": n,
            "status": "scoreable",
            "dependency_group": archive.dependency_group,
            "representatives": {
                orientation: [
                    archive.histograms[orientation][0].a,
                    archive.histograms[orientation][0].b,
                ]
                for orientation in ("first", "second")
            },
            "samples_per_orientation": int(archive.metadata["samples_per_pair"]),
            "batch_count": len(batch_ids),
            "p_bar": metrics["p_bar"],
            "activation_curves_at_p_bar": {
                "first": {
                    "root": full["first_root"],
                    "F1": metrics["first_F1"],
                    "F2": metrics["first_F2"],
                    "M": metrics["first_F1"] + metrics["first_F2"] - 1.0,
                },
                "second": {
                    "root": full["second_root"],
                    "F1": metrics["second_F1"],
                    "F2": metrics["second_F2"],
                    "M": metrics["second_F1"] + metrics["second_F2"] - 1.0,
                },
                "delta_first_minus_second": {
                    "delta_F1": metrics["delta_F1"],
                    "delta_F2": metrics["delta_F2"],
                    "delta_M": metrics["delta_M"],
                    "pooled_M_prime": metrics["pooled_M_prime"],
                },
            },
            "root_shift_decomposition": {
                "delta_p1": metrics["delta_p1"],
                "delta_p2": metrics["delta_p2"],
                "linearized_root_gap": metrics["linearized_root_gap"],
                "actual_root_gap": metrics["actual_root_gap"],
                "nonlinear_closure_residual": metrics["nonlinear_closure_residual"],
            },
            "angular_normalized": {
                "delta_cos4": full["delta_cos4"],
                **{name.removeprefix("angular_"): metrics[name] for name in BASE_METRICS if name.startswith("angular_")},
            },
            "joint_rank_coordinates": full["joint"],
            "qualitative_decomposition": _qualitative(metrics, standard_errors),
            "estimate_vector_order": list(FULL_METRICS),
            "estimate_vector": [metrics[name] for name in FULL_METRICS],
            "standard_errors": standard_errors,
            "jackknife_covariance": covariance,
            "identity_audit": full["audit"],
            "provenance": {
                "source_commit": archive.metadata["git_commit"],
                "seed": archive.metadata["seed"],
                "counter_first": archive.metadata["replica_counter_first"],
                "counter_last_exclusive": archive.metadata["replica_counter_last_exclusive"],
                "inputs": {
                    kind: {
                        **archive.input_provenance[kind],
                        "sha256": sha256(path),
                    }
                    for kind, path in archive.paths.items()
                },
            },
        }
        return {
            "public": result,
            "deleted": [row["metrics"] for row in deleted],
            "batch_signature": [
                [row.batch, row.samples] for row in archive.histograms["first"]
            ],
        }
    except (ArchiveNotScoreable, ArithmeticError, KeyError, TypeError, json.JSONDecodeError) as exc:
        return {
            "public": {
                "N": n,
                "status": "not_scoreable",
                "dependency_group": str(entry.get("dependency_group", "unknown")),
                "reason": str(exc),
                "inputs": {
                    kind: str(entry.get(kind, ""))
                    for kind in ("histogram", "moments", "metadata")
                },
            },
            "deleted": None,
            "batch_signature": None,
        }


def _dependency_summary(results: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for result in results:
        public = result["public"]
        if public["status"] == "scoreable":
            groups.setdefault(public["dependency_group"], []).append(result)
    scored = [
        result for result in results if result["public"]["status"] == "scoreable"
    ]
    for index, left in enumerate(scored):
        left_provenance = left["public"]["provenance"]
        for right in scored[index + 1:]:
            right_provenance = right["public"]["provenance"]
            if int(left_provenance["seed"]) != int(right_provenance["seed"]):
                continue
            overlap = max(
                int(left_provenance["counter_first"]),
                int(right_provenance["counter_first"]),
            ) < min(
                int(left_provenance["counter_last_exclusive"]),
                int(right_provenance["counter_last_exclusive"]),
            )
            same_group = (
                left["public"]["dependency_group"]
                == right["public"]["dependency_group"]
            )
            if overlap and not same_group:
                raise ValueError(
                    "partially overlapping counter streams must share one dependency_group"
                )
    output = []
    for name, members in groups.items():
        signatures = {tuple(tuple(row) for row in member["batch_signature"]) for member in members}
        if len(members) > 1 and len(signatures) != 1:
            raise ValueError(f"dependency group {name} is not batch aligned")
        output.append({
            "id": name,
            "sizes": [member["public"]["N"] for member in members],
            "rule": "aligned_delete_one_full_covariance" if len(members) > 1 else "independent_archive",
            "independent_evidence_units": 1,
        })
    return output


def render(
    manifest_path: Path,
    manifest: Mapping[str, Any],
    results: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    dependency_groups = _dependency_summary(results)
    scoreable = [result for result in results if result["public"]["status"] == "scoreable"]
    observations = [
        {"N": result["public"]["N"], "metric": metric}
        for result in scoreable
        for metric in DECISION_METRICS
    ]
    estimates = [
        result["deleted"] and result["public"]["estimate_vector"][
            result["public"]["estimate_vector_order"].index(metric)
        ]
        for result in scoreable
        for metric in DECISION_METRICS
    ]
    covariance: list[list[float]] = []
    for left_result in scoreable:
        for left_metric in DECISION_METRICS:
            row = []
            for right_result in scoreable:
                for right_metric in DECISION_METRICS:
                    if left_result["public"]["dependency_group"] != right_result["public"]["dependency_group"]:
                        row.append(0.0)
                    else:
                        row.append(jackknife_covariance(
                            [item[left_metric] for item in left_result["deleted"]],
                            [item[right_metric] for item in right_result["deleted"]],
                        ))
            covariance.append(row)
    by_n = {str(result["public"]["N"]): result["public"] for result in results}
    return {
        "schema": SCHEMA,
        "status": "retrospective existing-data decomposition; no new simulation and no exponent fit",
        "exact_mapping": {
            "K1": "K_minus",
            "K2": "K_plus",
            "F1": "E[Pr(Binomial(N,p)>=K1)]",
            "F2": "E[Pr(Binomial(N,p)>=K2)]",
            "M": "-1+F1+F2",
            "C": "(K1+K2)/2",
            "G": "K2-K1",
        },
        "size_order": [int(value) for value in manifest["required_sizes"]],
        "scoreable_sizes": [result["public"]["N"] for result in scoreable],
        "not_scoreable_sizes": [
            result["public"]["N"] for result in results
            if result["public"]["status"] == "not_scoreable"
        ],
        "by_N": by_n,
        "dependency_groups": dependency_groups,
        "decision_covariance": {
            "metric_order_with_N": observations,
            "estimate_vector": estimates,
            "jackknife_covariance": covariance,
            "zero_cross_group_entries_mean_independent_archives": True,
        },
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
        },
        "interpretation_guard": (
            "K1/K2 is an exact finite-threshold reinterpretation.  The magnitude labels "
            "are descriptive per size; shared dependency groups are one evidence unit, and "
            "no common scaling exponent or continuum operator identity is inferred here."
        ),
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# K1/K2 two-activation H4 decomposition",
        "",
        "This is a retrospective reanalysis of archived threshold-rank histograms. "
        "It generates no Monte Carlo samples and fits no exponent.",
        "",
        "Exact convention: `K1=K_minus`, `K2=K_plus`, "
        "`M(p)=-1+F1(p)+F2(p)`.",
        "",
        "| N | status | K1 share | interaction | angular delta p1 (z) | angular delta p2 (z) | angular root gap | closure residual |",
        "|---:|:---|---:|:---|---:|---:|---:|---:|",
    ]
    for n in payload["size_order"]:
        row = payload["by_N"][str(n)]
        if row["status"] != "scoreable":
            lines.append(f"| {n} | not scoreable: {row['reason']} | | | | | | |")
            continue
        qualitative = row["qualitative_decomposition"]
        angular = row["angular_normalized"]
        lines.append(
            "| {n} | scoreable | {share:.3f} | {interaction} | {p1:.6e} ({z1:.2f}) | "
            "{p2:.6e} ({z2:.2f}) | {gap:.6e} | {residual:.3e} |".format(
                n=n,
                share=qualitative["absolute_K1_fraction"],
                interaction=qualitative["component_interaction"],
                p1=angular["delta_p1"],
                p2=angular["delta_p2"],
                z1=qualitative["angular_delta_p1_z"],
                z2=qualitative["angular_delta_p2_z"],
                gap=angular["actual_root_gap"],
                residual=angular["nonlinear_closure_residual"],
            )
        )
    scoreable_rows = [
        payload["by_N"][str(n)] for n in payload["scoreable_sizes"]
    ]
    first_count = sum(
        row["qualitative_decomposition"]["descriptive_dominant_component"]
        == "first_activation_K1"
        for row in scoreable_rows
    )
    second_count = sum(
        row["qualitative_decomposition"]["descriptive_dominant_component"]
        == "second_activation_K2"
        for row in scoreable_rows
    )
    shared_count = len(scoreable_rows) - first_count - second_count
    cancelling = sum(
        row["qualitative_decomposition"]["component_interaction"] == "cancelling"
        for row in scoreable_rows
    )
    resolved_k1 = sum(
        abs(row["qualitative_decomposition"]["angular_delta_p1_z"]) >= 2.0
        for row in scoreable_rows
    )
    resolved_k2 = sum(
        abs(row["qualitative_decomposition"]["angular_delta_p2_z"]) >= 2.0
        for row in scoreable_rows
    )
    lines.extend([
        "",
        "## Descriptive result",
        "",
        f"Across the {len(scoreable_rows)} scoreable sizes, magnitude classification is "
        f"K1-dominant at {first_count}, K2-dominant at {second_count}, and shared at "
        f"{shared_count}; {cancelling} sizes show opposite-sign component cancellation.",
        f"Using the within-size delete-one standard errors only, K1 has |z|>=2 at "
        f"{resolved_k1} sizes and K2 at {resolved_k2}. Unresolved component point "
        "estimates remain next-target clues rather than confirmed transitions.",
        "These counts map the decomposition; they are not independent-evidence votes. "
        "The dependency groups and their member sizes are recorded explicitly in JSON.",
        "",
        "The `nonlinear closure residual` is the observed orientation root gap minus "
        "`delta_p1+delta_p2`.  Its smallness diagnoses the local linearization only; it "
        "does not establish a continuum exponent or operator identity.",
        "",
        "## Joint rank coordinates",
        "",
        "Every scoreable archive also reports means, variances and covariance for "
        "`C=(K1+K2)/2` and `G=K2-K1` separately for both orientations and their pooled "
        "mixture.  These quantities use the archived paired integer moments, not a "
        "product reconstructed from marginal histograms.",
        "",
        "## Covariance and provenance",
        "",
        "All nonlinear coordinates are recomputed after deleting the same batch from "
        "both orientations.  The JSON contains each size's full covariance over the "
        "declared estimate vector and a cross-size covariance over the decision metrics. "
        "Entries from distinct counter intervals are zero by design; shared streams use "
        "aligned delete-one covariance.  Input SHA256 values and source commits are "
        "stored under each size.",
        "",
        "This is a canonical Phase-D state-coordinate decomposition. It does not "
        "construct the Phase-E `J_top` versus `J_bulk` comparison and is not an "
        "outward-rounded interval or SOS certificate.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path, default=root / "analysis/two_activation_h4_manifest.yaml"
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    try:
        manifest = load_manifest(args.manifest)
        worker_args = [(str(root), entry) for entry in manifest["runs"]]
        if args.workers == 1:
            results = [analyze_archive(*item) for item in worker_args]
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(worker_args))) as pool:
                results = list(pool.map(_analyze_archive_star, worker_args))
        payload = render(args.manifest, manifest, results)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        raise SystemExit(str(exc)) from exc
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(render_markdown(payload), encoding="utf-8")
    print(args.output_json)
    print(args.output_md)
    return 0


def _analyze_archive_star(arguments: tuple[str, Mapping[str, Any]]) -> dict[str, Any]:
    return analyze_archive(*arguments)


if __name__ == "__main__":
    raise SystemExit(main())
