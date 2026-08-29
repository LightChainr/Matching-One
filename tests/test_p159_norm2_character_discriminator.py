#!/usr/bin/env python3
"""Exact checks for the prospective Issue #159 norm-2 discriminator."""

from __future__ import annotations

import cmath
from fractions import Fraction
import math
from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "predictions" / "p159_norm2_character_discriminator_20260829.yaml"


def matmul(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [
        [sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def det(matrix: list[list[int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def tau(matrix: list[list[int]]) -> complex:
    v1 = complex(matrix[0][0], matrix[1][0])
    v2 = complex(matrix[0][1], matrix[1][1])
    return v2 / v1


class P159Norm2CharacterDiscriminatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))

    def test_children_are_exact_norm2_embeddings(self) -> None:
        multiplier = self.payload["embedding"]["left_matrix"]
        self.assertEqual(det(multiplier), 2)
        for row in self.payload["source_rows"]:
            parent = row["parent_period_matrix"]
            child = row["child_period_matrix"]
            with self.subTest(row=row["id"]):
                self.assertEqual(matmul(multiplier, parent), child)
                self.assertEqual(abs(det(parent)), row["parent_N"])
                self.assertEqual(abs(det(child)), row["child_N"])
                self.assertEqual(row["child_N"], 2 * row["parent_N"])
                self.assertAlmostEqual(tau(parent).real, tau(child).real, places=14)
                self.assertAlmostEqual(tau(parent).imag, tau(child).imag, places=14)

    def test_character_phases_and_local_spin_bound(self) -> None:
        theta = math.pi / 4
        self.assertAlmostEqual((cmath.exp(-1j * 4 * theta)).real, -1.0, places=14)
        self.assertAlmostEqual((cmath.exp(+1j * 8 * theta)).real, +1.0, places=14)
        alpha_h4 = (4 - 2) / 2
        alpha_h8_min = (8 - 2) / 2
        self.assertEqual(alpha_h4, 1)
        self.assertEqual(alpha_h8_min, 3)
        self.assertGreater(alpha_h8_min, 1)

    def test_frozen_ratios_and_source_targets(self) -> None:
        ratios = {
            "rank4_H4_Nminus1": -0.5,
            "nonlocal_even_Nminus1": 0.5,
            "quadratic_H4xH4_Nminus2": 0.25,
            "local_H8_bound_saturated_Nminus3": 0.125,
        }
        models = {row["id"]: row for row in self.payload["fixed_models_in_score_order"]}
        self.assertEqual(list(models), list(ratios))
        for source in self.payload["source_rows"]:
            parent = float(source["parent_C"])
            parent_se = float(source["parent_C_SE"])
            for model_id, ratio in ratios.items():
                with self.subTest(source=source["id"], model=model_id):
                    target = source["source_based_targets"][model_id]
                    self.assertAlmostEqual(float(Fraction(target["ratio"])), ratio, places=15)
                    self.assertAlmostEqual(float(target["C"]), ratio * parent, places=15)
                    self.assertAlmostEqual(float(target["source_SE"]), abs(ratio) * parent_se, places=15)

    def test_manifest_is_registered(self) -> None:
        registry = yaml.safe_load(
            (ROOT / "analysis" / "artifact_registry.yaml").read_text(encoding="utf-8")
        )
        paths = {row["path"] for row in registry["frozen_predictions"]}
        self.assertIn(str(MANIFEST.relative_to(ROOT)), paths)


if __name__ == "__main__":
    unittest.main()
