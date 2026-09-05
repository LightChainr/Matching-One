#!/usr/bin/env python3
"""Tests for the bounded Issue 31 resource telemetry primitive."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import transfer_resource_probe as probe  # noqa: E402


class TransferResourceProbeTests(unittest.TestCase):
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
