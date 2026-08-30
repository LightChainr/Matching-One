from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p418_root_translation_semantics as p418  # noqa: E402


RESULT = ROOT / "results/exact/P418-root-translation-semantics/certificate.json"


class P418RootTranslationSemanticsTests(unittest.TestCase):
    def test_cyclotomic_exact_arithmetic(self) -> None:
        self.assertTrue(p418.cyclotomic_zero([7, 7, 7, 7, 7]))
        self.assertFalse(p418.cyclotomic_zero([1, 0, 0, 0, 0]))
        left = [1, 0, -1, 2, 0]
        right = [0, 2, 1, 0, -2]
        observed = p418.cyclotomic_value(p418.convolve_coefficients(left, right), 1)
        expected = p418.cyclotomic_value(left, 1) * p418.cyclotomic_value(right, 1)
        self.assertLess(abs(observed - expected), 1e-12)

    def test_exact_root_translation_certificate(self) -> None:
        observed = p418.build_certificate()
        self.assertTrue(observed["summary"]["passed"])
        self.assertEqual(observed["summary"]["component_signature_failures"], 0)
        self.assertEqual(observed["summary"]["root_scalar_failures"], 0)
        self.assertEqual(observed["summary"]["section_or_gauge_failures"], 0)
        self.assertEqual(observed["summary"]["translation_orbit_factorization_failures"], 0)
        self.assertEqual(observed["summary"]["archive_semantic_provenance_failures"], 0)
        self.assertGreater(observed["summary"]["fixed_configuration_factorization_counterexamples"], 0)

    def test_committed_certificate_recomputes(self) -> None:
        observed = json.loads(RESULT.read_text())
        self.assertEqual(p418.build_certificate(), observed)


if __name__ == "__main__":
    unittest.main()
