#!/usr/bin/env python3
"""Tests for the bounded Issue 31 resource telemetry primitive."""

from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import transfer_resource_probe as probe  # noqa: E402


class TransferResourceProbeTests(unittest.TestCase):
    def test_success_records_direct_rss_and_stream_digests(self) -> None:
        code = "import time; payload=bytearray(2000000); print('probe-ok'); time.sleep(0.05)"
        result = probe.run_probe(
            [sys.executable, "-c", code],
            label="synthetic-success",
            cwd=ROOT,
            poll_interval_ms=5,
            environment={"OMP_NUM_THREADS": "2", "IGNORED_SECRET": "not-recorded"},
        )
        self.assertEqual(result["process"]["exit_code"], 0)
        self.assertGreater(result["memory"]["peak_rss_kib"], 0)
        self.assertEqual(result["memory"]["rusage_observations"], 1)
        self.assertEqual(result["stdout"]["sha256"], hashlib.sha256(b"probe-ok\n").hexdigest())
        self.assertEqual(result["stdout"]["bytes"], len(b"probe-ok\n"))
        self.assertEqual(result["environment"]["allowlisted"], {"OMP_NUM_THREADS": "2"})

    def test_nonzero_command_is_measured_without_becoming_a_success(self) -> None:
        result = probe.run_probe(
            [sys.executable, "-c", "import time; time.sleep(0.03); raise SystemExit(7)"],
            label="synthetic-failure",
            cwd=ROOT,
            poll_interval_ms=5,
            environment={},
        )
        self.assertEqual(result["process"]["exit_code"], 7)
        probe.validate_result(result)

    def test_tampered_wall_time_and_rss_fail_closed(self) -> None:
        result = probe.run_probe(
            [sys.executable, "-c", "import time; time.sleep(0.03)"],
            label="synthetic-tamper",
            cwd=ROOT,
            poll_interval_ms=5,
            environment={},
        )
        changed = copy.deepcopy(result)
        changed["timing"]["wall_seconds"] = -1
        with self.assertRaisesRegex(ValueError, "wall time"):
            probe.validate_result(changed)
        changed = copy.deepcopy(result)
        changed["memory"]["peak_rss_kib"] = 0
        with self.assertRaisesRegex(ValueError, "RSS"):
            probe.validate_result(changed)

    def test_shell_or_environment_scope_expansion_fails(self) -> None:
        result = probe.run_probe(
            [sys.executable, "-c", "import time; time.sleep(0.03)"],
            label="synthetic-policy",
            cwd=ROOT,
            poll_interval_ms=5,
            environment={},
        )
        changed = copy.deepcopy(result)
        changed["command"]["shell"] = True
        with self.assertRaisesRegex(ValueError, "shell"):
            probe.validate_result(changed)
        changed = copy.deepcopy(result)
        changed["environment"]["allowlisted"]["API_TOKEN"] = "forbidden"
        with self.assertRaisesRegex(ValueError, "non-allowlisted"):
            probe.validate_result(changed)


if __name__ == "__main__":
    unittest.main()
