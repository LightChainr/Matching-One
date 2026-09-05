
from __future__ import annotations
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p50_sprime_n290_typed as typed  # noqa: E402


class P50SprimeTypedTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / typed.SEMANTIC_GATE
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / typed.SEMANTIC_GATE, destination)
        return directory, root

    def test_gate_registers_exact_identity(self) -> None:
        gate, source, target, transform = typed.load_semantic_gate(ROOT)
        self.assertEqual(gate["target_size"], 290)
        self.assertEqual(source, target)
        self.assertEqual((transform.scale, transform.offset), (1.0, 0.0))

    def test_delegated_result_is_numerically_unchanged(self) -> None:
        frozen = {
            "schema": "matching-one/P50-N290-Sprime-frozen-score/v1",
            "observable": "P4_S_prime",
            "N": 290,
            "observed": 0.125,
            "models": {
                "q2_even_scalar_correction": {"chi_square": 1.25},
                "rank2_jordan_log": {"chi_square": 0.75},
            },
        }
        with patch.object(typed, "validate_prediction_files"), patch.object(
            typed, "_run_frozen_kernel", return_value=frozen.copy()
        ) as runner:
            result = typed.score_typed(
                ROOT,
                Path("child.csv"),
                Path("q2.yaml"),
                Path("jordan.yaml"),
                runner=runner,
            )
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, frozen)
        self.assertEqual(semantics["applied_transform"]["scale"], 1.0)


if __name__ == "__main__":
    unittest.main()
