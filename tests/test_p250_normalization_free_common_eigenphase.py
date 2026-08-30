from __future__ import annotations

import cmath
import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p250_normalization_free_common_eigenphase import (  # noqa: E402
    common_eigenphase_contrast,
    q_projection,
    validate_manifest,
)


def vector(a_plus: complex, a_minus: complex, b_plus: complex, b_minus: complex) -> list[float]:
    return [
        a_plus.real, a_plus.imag, a_minus.real, a_minus.imag,
        b_plus.real, b_plus.imag, b_minus.real, b_minus.imag,
    ]


class P250NormalizedCommonEigenphaseTests(unittest.TestCase):
    def test_common_unit_eigenphase_has_zero_frozen_contrast(self) -> None:
        q = cmath.exp(0.7j)
        row = vector(q * (2 + 1j), 2 + 1j, q * (-1 + 3j), -1 + 3j)
        for value in common_eigenphase_contrast(row):
            self.assertAlmostEqual(value, 0.0, places=13)
        self.assertAlmostEqual(q_projection(row).real, q.real)
        self.assertAlmostEqual(q_projection(row).imag, q.imag)

    def test_channel_phase_mismatch_breaks_complex_cross_product(self) -> None:
        row = vector(1j, 1 + 0j, -1 + 0j, 1 + 0j)
        contrast = common_eigenphase_contrast(row)
        self.assertNotAlmostEqual(abs(complex(contrast[0], contrast[1])), 0.0)

    def test_common_nonunit_gain_breaks_norm_contrast(self) -> None:
        row = vector(2j, 1 + 0j, 2j, 1 + 0j)
        contrast = common_eigenphase_contrast(row)
        self.assertAlmostEqual(contrast[0], 0.0)
        self.assertAlmostEqual(contrast[1], 0.0)
        self.assertGreater(contrast[2], 0.0)

    def test_frozen_manifest_hashes_validate(self) -> None:
        manifest_path = ROOT / "analysis/p250_normalization_free_common_eigenphase_20260830.json"
        manifest = validate_manifest(ROOT, manifest_path)
        self.assertEqual(manifest["eligibility"]["separation"], 1)
        self.assertEqual(manifest["charged_channels"]["primitive"], ["C113", "C122"])
        self.assertIn("reopen H4/H8/H12 ordinary-channel voting", manifest["forbidden"])

    def test_checked_result_retains_complex_covariance_when_present(self) -> None:
        path = ROOT / "results/huawei-20260830/P250-normalization-free-common-eigenphase/score.json"
        if not path.exists():
            self.skipTest("normalization-free score not yet revealed")
        result = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(len(result["omega_full_shared_batch_covariance_8x8"]), 8)
        self.assertEqual(
            result["primary_common_unit_eigenphase"]["contrast_order"],
            ["cross_re", "cross_im", "norm_difference"],
        )
        self.assertEqual(len(result["primary_common_unit_eigenphase"]["contrast_full_shared_batch_covariance_3x3"]), 3)
        self.assertAlmostEqual(
            result["primary_common_unit_eigenphase"]["survival_p"],
            0.3717090718193832,
        )
        self.assertEqual(
            result["primary_common_unit_eigenphase"]["decision"],
            "common_unit_eigenphase_survives",
        )
        self.assertGreater(
            result["descriptive_eigenphase"]["common_projection"]["phase_standard_error"],
            math.pi,
        )


if __name__ == "__main__":
    unittest.main()
