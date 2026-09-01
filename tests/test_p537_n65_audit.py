#!/usr/bin/env python3
"""Check the retained N65 covariance contract without replaying raw MC."""

from __future__ import annotations

import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results" / "p537-contact-stage-n65" / "AUDIT.json"
FROZEN = ROOT / "results" / "p537-contact-stage-n65" / "result.json"


def matmul(left: list[list[float]], right: list[list[float]]) -> list[list[float]]:
    return [
        [sum(a * b for a, b in zip(row, col)) for col in zip(*right)]
        for row in left
    ]


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    return [list(col) for col in zip(*matrix)]


class P537N65AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.audit = json.loads(AUDIT.read_text(encoding="utf-8"))
        cls.frozen = json.loads(FROZEN.read_text(encoding="utf-8"))

    def assert_close(self, left: float, right: float, scale: float = 1.0) -> None:
        self.assertLessEqual(abs(left - right), 5e-13 * max(scale, abs(left), abs(right)))

    def test_frozen_primary_is_unchanged(self) -> None:
        primary = self.audit["primary"]
        self.assertEqual(primary["decision"], self.frozen["primary"]["decision"])
        self.assert_close(primary["Delta"], self.frozen["primary"]["Delta"], 1e-13)
        self.assertEqual(primary["theta"]["value"], -1.0)
        self.assertTrue(primary["theta"]["saturated_by_sign_pattern"])

    def test_covariance_groups_and_collapse_close(self) -> None:
        cov = self.audit["matrix"]["covariance"]
        total = cov["total"]
        for i in range(6):
            for j in range(6):
                self.assert_close(total[i][j], total[j][i], 1e-15)
                self.assert_close(total[i][j], cov["production"][i][j] + cov["baseline"][i][j], 1e-15)
        transform = self.audit["primary"]["collapse_transform"]
        collapsed = matmul(matmul(transform, total), transpose(transform))
        stored = self.audit["primary"]["covariance"]["total"]
        for i in range(4):
            for j in range(4):
                self.assert_close(collapsed[i][j], stored[i][j], 1e-15)

    def test_exposure_density_and_selected_total_close(self) -> None:
        matrix = sum(self.audit["matrix"]["estimate"], [])
        exposure = sum(self.audit["positive_exposure"]["estimate"], [])
        density = sum(self.audit["conditional_signed_density"]["estimate"], [])
        self.assertTrue(all(value > 0 for value in exposure))
        for value, mass, conditional in zip(matrix, exposure, density):
            self.assert_close(value, mass * conditional, 1e-7)
        selected = self.audit["selected_carrier_total"]
        self.assert_close(selected["estimate"], sum(matrix), 1e-7)
        self.assert_close(selected["se"] ** 2, selected["variance_total"], 1e-14)
        self.assertLess(selected["ci95"][1], 0.0)
        self.assertTrue(math.isfinite(selected["z"]))


if __name__ == "__main__":
    unittest.main()
