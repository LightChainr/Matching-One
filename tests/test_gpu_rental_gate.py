
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from validate_gpu_rental_gate import evaluate_measurement, validate_contract  # noqa: E402


class GpuRentalGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "gpu_rental_gate_manifest.json").read_text(encoding="utf-8")
        )
        cls.source = (ROOT / cls.contract["baseline_source"]["path"]).read_bytes()

    def test_checked_in_contract_is_exact_and_not_authorized(self) -> None:
        result = validate_contract(self.contract, source=self.source)
        self.assertEqual(result["cpu_paired_permutations"], 50_000_000)
        self.assertEqual(result["cpu_elapsed_seconds"], "326/5")
        self.assertEqual(result["cpu_paired_permutations_per_second"], "125000000/163")
        self.assertEqual(
            result["minimum_gpu_paired_permutations_per_second"], "625000000/163"
        )
        self.assertFalse(result["gpu_measurement_present"])
        self.assertEqual(result["rental"], "not_authorized")

    def test_exact_five_x_synthetic_measurement_passes_route_one(self) -> None:
        result = evaluate_measurement(
            self.contract,
            gpu_paired_permutations=50_000_000,
            gpu_elapsed_seconds=Fraction(326, 25),
            end_to_end_timing=True,
            output_contract_equal=True,
            deterministic_regression_passed=True,
        )
        self.assertEqual(result["speedup"], "5")
        self.assertTrue(result["throughput_gate_passed"])
        self.assertEqual(result["rental"], "authorized_by_route_1")

    def test_speed_and_all_three_validation_gates_are_required(self) -> None:
        below = evaluate_measurement(
            self.contract,
            gpu_paired_permutations=50_000_000,
            gpu_elapsed_seconds=Fraction(14),
            end_to_end_timing=True,
            output_contract_equal=True,
            deterministic_regression_passed=True,
        )
        self.assertFalse(below["throughput_gate_passed"])
        self.assertEqual(below["rental"], "not_authorized")

        names = ["end_to_end_timing", "output_contract_equal", "deterministic_regression_passed"]
        for missing in names:
            flags = {name: name != missing for name in names}
            with self.subTest(missing=missing):
                result = evaluate_measurement(
                    self.contract,
                    gpu_paired_permutations=50_000_000,
                    gpu_elapsed_seconds=Fraction(10),
                    **flags,
                )
                self.assertTrue(result["throughput_gate_passed"])
                self.assertEqual(result["rental"], "not_authorized")

    def test_source_and_exact_arithmetic_drift_fail_closed(self) -> None:
        contract = deepcopy(self.contract)
        contract["baseline_source"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "baseline source SHA-256 mismatch"):
            validate_contract(contract, source=self.source)

        contract = deepcopy(self.contract)
        contract["cpu_baseline"]["elapsed_seconds"] = "65.2"
        with self.assertRaisesRegex(ValueError, "not canonical"):
            validate_contract(contract, source=self.source)

        contract = deepcopy(self.contract)
        contract["derived_exact"]["minimum_gpu_paired_permutations_per_second"] = "1"
        with self.assertRaisesRegex(ValueError, "derived GPU threshold drifted"):
            validate_contract(contract, source=self.source)

        contract = deepcopy(self.contract)
        contract["end_to_end_gate"]["minimum_speedup"] = "4"
        with self.assertRaisesRegex(ValueError, "minimum speedup must remain 5"):
            evaluate_measurement(
                contract,
                gpu_paired_permutations=50_000_000,
                gpu_elapsed_seconds=Fraction(10),
                end_to_end_timing=True,
                output_contract_equal=True,
                deterministic_regression_passed=True,
            )

    def test_gpu_result_or_authorization_fields_are_forbidden(self) -> None:
        contract = deepcopy(self.contract)
        contract["estimated_speedup"] = "6"
        with self.assertRaisesRegex(ValueError, "forbidden GPU-result keys"):
            validate_contract(contract, source=self.source)

        contract = deepcopy(self.contract)
        contract["gpu_measurement"] = {"status": "measured"}
        with self.assertRaisesRegex(ValueError, "GPU measurement is checked in"):
            validate_contract(contract, source=self.source)


if __name__ == "__main__":
    unittest.main()
