import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_algebraic_model_problem import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    compile_problem,
    evaluate_equation,
    validate_result,
)


class BuildAlgebraicModelProblemTests(unittest.TestCase):
    def test_frozen_problem_is_canonical_and_exact(self) -> None:
        result = build_result()
        problem = result["problem"]
        self.assertEqual(len(problem["variable_order"]), 8)
        self.assertEqual(problem["maximum_total_degree"], 2)
        self.assertTrue(result["supplied_solution_check"]["all_zero"])

    def test_compilation_is_deterministic(self) -> None:
        first = build_result()["problem"]
        second = build_result()["problem"]
        self.assertEqual(first, second)

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_canonical_algebraic_model_problem")
        self.assertEqual(summary["variable_count"], 8)


if __name__ == "__main__":
    unittest.main()
