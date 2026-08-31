#!/usr/bin/env python3
"""Run one command and emit bounded, versioned resource telemetry.

This is a measurement primitive for Issue 31, not the transfer-matrix resource
study itself. Peak RSS comes from the command-specific rusage returned by
``os.wait4``; descendant aggregation follows the host's wait4 semantics.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, BinaryIO, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "matching-one/transfer-resource-telemetry/v1"
HEX_DIGITS = frozenset("0123456789abcdef")
ENVIRONMENT_ALLOWLIST = (
    "OMP_NUM_THREADS",
    "OMP_PROC_BIND",
    "OMP_PLACES",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _stream_summary(handle: BinaryIO) -> dict[str, Any]:
    handle.flush()
    handle.seek(0)
    digest = hashlib.sha256()
    size = 0
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
        size += len(chunk)
    return {"bytes": size, "sha256": digest.hexdigest(), "content_retained": False}


def _git_text(*args: str) -> Optional[str]:
    completed = subprocess.run(
        ["git", "-C", str(ROOT), *args],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
        text=True,
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def _repository_provenance() -> dict[str, Any]:
    commit = _git_text("rev-parse", "HEAD")
    status = _git_text("status", "--porcelain", "--untracked-files=no")
    return {
        "commit": commit,
        "tracked_files_dirty": None if status is None else bool(status),
        "root": str(ROOT),
    }


def _resolve_executable(argv0: str, cwd: Path) -> Optional[Path]:
    if os.sep in argv0:
        candidate = Path(argv0)
        if not candidate.is_absolute():
            candidate = cwd / candidate
        candidate = candidate.resolve()
        return candidate if candidate.is_file() else None
    resolved = shutil.which(argv0)
    return Path(resolved).resolve() if resolved else None


def _ru_maxrss_kib(raw_value: float) -> int:
    """Normalize wait4 ru_maxrss to KiB on supported Unix hosts."""

    if sys.platform == "darwin":
        return int(raw_value / 1024.0)
    return int(raw_value)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def run_probe(
    argv: Sequence[str],
    *,
    label: str,
    cwd: Path,
    poll_interval_ms: int = 10,
    environment: Optional[Mapping[str, str]] = None,
) -> dict[str, Any]:
    """Run exactly one argv command and return validated telemetry."""

    _require(bool(argv), "command argv must not be empty")
    _require(all(isinstance(value, str) and value for value in argv), "argv entries must be nonempty strings")
    _require(isinstance(label, str) and bool(label.strip()), "label must be nonempty")
    _require(cwd.is_dir(), "working directory does not exist")
    _require(1 <= poll_interval_ms <= 1000, "poll interval must be in [1,1000] ms")
    env = dict(os.environ if environment is None else environment)
    executable = _resolve_executable(argv[0], cwd)
    _require(executable is not None, "command executable cannot be resolved")
    executable_digest = _sha256_file(executable)
    started_utc = _utc_now()
    start = time.monotonic()
    _require(hasattr(os, "wait4"), "os.wait4 is required for command-specific RSS telemetry")

    with tempfile.TemporaryFile() as stdout_handle, tempfile.TemporaryFile() as stderr_handle:
        process = subprocess.Popen(
            list(argv),
            cwd=str(cwd),
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=stdout_handle,
            stderr=stderr_handle,
            shell=False,
        )
        while True:
            waited_pid, wait_status, usage = os.wait4(process.pid, os.WNOHANG)
            if waited_pid == process.pid:
                exit_code = os.waitstatus_to_exitcode(wait_status)
                process.returncode = exit_code
                peak_rss_kib = _ru_maxrss_kib(usage.ru_maxrss)
                break
            time.sleep(poll_interval_ms / 1000.0)
        wall_seconds = time.monotonic() - start
        stdout_summary = _stream_summary(stdout_handle)
        stderr_summary = _stream_summary(stderr_handle)

    result = {
        "schema": SCHEMA,
        "label": label.strip(),
        "command": {
            "argv": list(argv),
            "shell": False,
            "cwd": str(cwd.resolve()),
            "executable": str(executable),
            "executable_sha256": executable_digest,
        },
        "process": {"exit_code": exit_code},
        "timing": {
            "started_utc": started_utc,
            "ended_utc": _utc_now(),
            "wall_seconds": wall_seconds,
            "clock": "time.monotonic",
            "poll_interval_ms": poll_interval_ms,
        },
        "memory": {
            "peak_rss_kib": peak_rss_kib,
            "rusage_observations": 1,
            "source": "os.wait4(pid) rusage.ru_maxrss",
            "scope": "specific waited command process; descendant aggregation follows host wait4 semantics",
        },
        "stdout": stdout_summary,
        "stderr": stderr_summary,
        "host": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "logical_cpu_count": os.cpu_count(),
        },
        "environment": {
            "allowlisted": {key: env[key] for key in ENVIRONMENT_ALLOWLIST if key in env},
            "policy": "only fixed non-secret parallel-runtime variables are retained",
        },
        "repository": _repository_provenance(),
        "claim_boundary": {
            "included": "single-command wall time, direct-process RSS, exit and provenance telemetry",
            "excluded": "pipeline phases, descendant-process attribution, state counts, branching, extrapolation, or hardware decisions",
            "parent_issue": "remain open",
        },
    }
    validate_result(result)
    return result


def validate_result(result: Mapping[str, Any]) -> None:
    """Fail closed on incomplete, ill-typed, or overclaimed telemetry."""

    _require(result.get("schema") == SCHEMA, "unknown telemetry schema")
    _require(isinstance(result.get("label"), str) and bool(result["label"]), "missing label")
    command = result.get("command", {})
    _require(isinstance(command.get("argv"), list) and bool(command["argv"]), "missing argv")
    _require(command.get("shell") is False, "shell execution is forbidden")
    digest = command.get("executable_sha256")
    _require(isinstance(digest, str) and len(digest) == 64 and set(digest) <= HEX_DIGITS, "bad executable digest")
    process = result.get("process", {})
    _require(isinstance(process.get("exit_code"), int), "exit code must be an integer")
    timing = result.get("timing", {})
    _require(isinstance(timing.get("wall_seconds"), (int, float)) and timing["wall_seconds"] >= 0, "bad wall time")
    _require(timing.get("clock") == "time.monotonic", "unexpected timing clock")
    _require(isinstance(timing.get("poll_interval_ms"), int), "missing poll interval")
    memory = result.get("memory", {})
    _require(isinstance(memory.get("peak_rss_kib"), int) and memory["peak_rss_kib"] > 0, "no positive RSS sample")
    _require(memory.get("rusage_observations") == 1, "unexpected rusage observation count")
    _require(memory.get("source") == "os.wait4(pid) rusage.ru_maxrss", "unexpected RSS source")
    _require("host wait4 semantics" in memory.get("scope", ""), "RSS scope is ambiguous")
    for stream_name in ("stdout", "stderr"):
        stream = result.get(stream_name, {})
        _require(isinstance(stream.get("bytes"), int) and stream["bytes"] >= 0, "%s byte count invalid" % stream_name)
        stream_digest = stream.get("sha256")
        _require(
            isinstance(stream_digest, str) and len(stream_digest) == 64 and set(stream_digest) <= HEX_DIGITS,
            "%s digest invalid" % stream_name,
        )
        _require(stream.get("content_retained") is False, "%s content retention changed" % stream_name)
    allowlisted = result.get("environment", {}).get("allowlisted", {})
    _require(set(allowlisted) <= set(ENVIRONMENT_ALLOWLIST), "non-allowlisted environment value retained")
    boundary = result.get("claim_boundary", {})
    _require(boundary.get("parent_issue") == "remain open", "parent issue boundary changed")
    excluded = boundary.get("excluded", "")
    for phrase in ("pipeline phases", "state counts", "extrapolation", "hardware decisions"):
        _require(phrase in excluded, "missing excluded boundary: %s" % phrase)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--cwd", type=Path, default=Path.cwd())
    parser.add_argument("--poll-interval-ms", type=int, default=10)
    parser.add_argument("--output", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    result = run_probe(
        command,
        label=args.label,
        cwd=args.cwd,
        poll_interval_ms=args.poll_interval_ms,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        sys.stdout.write(rendered)
    else:
        args.output.write_text(rendered, encoding="utf-8")
    return 0 if result["process"]["exit_code"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
