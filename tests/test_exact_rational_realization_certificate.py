import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_rational_realization_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    matrix_rank,
    validate_result,
    verify_realization,
)


class ExactRationalRealizationCertificateTests(unittest.TestCase):
    def test_generic_diagonal_realization_verifies_exactly(self) -> None:
        descriptor = {
            "generator": [["1", "0"], ["0", "2"]],
            "channels": [{"id": "sum", "source": ["1", "1"], "readout": ["1", "1"], "moments": ["2", "3", "5"]}],
        }
        verified = verify_realization(descriptor)
        self.assertEqual(verified["channels"][0]["moments"], ["2", "3", "5"])
        self.assertEqual(verified["dimension"], 2)
        self.assertTrue(verified["minimal_on_declared_typed_rows"])

    def test_exact_rank_handles_rectangular_matrices(self) -> None:
        self.assertEqual(matrix_rank([[1, 2, 3], [2, 4, 6]]), 1)
        self.assertEqual(matrix_rank([[1, 0], [0, 1], [1, 1]]), 2)

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        summary = validate_result(checked)
        self.assertEqual(summary["dimension"], 3)
        self.assertEqual(summary["channel_count"], 2)
        self.assertEqual(summary["reachability_rank"], 3)
        self.assertEqual(summary["observability_rank"], 3)


if __name__ == "__main__":
    unittest.main()
