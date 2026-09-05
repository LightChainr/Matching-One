import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_matrix_relation_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    validate_result,
    verify_relations,
)


class ExactMatrixRelationCertificateTests(unittest.TestCase):
    def test_frozen_relation_is_verified_exactly(self) -> None:
        result = build_result()
        self.assertEqual(result["verification"]["status"], "exact_matrix_relations_verified")
        self.assertEqual(result["verification"]["nilpotent_rank"], 1)
        self.assertFalse(result["solver_invoked"])

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_exact_matrix_relation_certificate")
        self.assertEqual(summary["nilpotent_rank"], 1)


if __name__ == "__main__":
    unittest.main()
