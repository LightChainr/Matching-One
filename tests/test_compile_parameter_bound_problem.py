import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from compile_parameter_bound_problem import DEFAULT_OUTPUT, compile_bounds, evaluate, frozen_bounds, validate_result  # noqa: E402


class CompileParameterBoundProblemTests(unittest.TestCase):
    def test_frozen_box_compiles_archimedean_bound(self) -> None:
        problem = compile_bounds(frozen_bounds())
        self.assertEqual(problem["inequality_count"], 5)
        self.assertEqual(problem["archimedean_radius_squared"], "17/4")

    def test_box_endpoints_satisfy_compiled_inequalities(self) -> None:
        problem = compile_bounds(frozen_bounds())
        for assignment in (["-1", "-1/2"], ["2", "1/3"]):
            self.assertTrue(all(evaluate(item["terms"], assignment) >= 0 for item in problem["inequalities"]))

    def test_outside_assignment_violates_a_bound(self) -> None:
        problem = compile_bounds(frozen_bounds())
        values = [evaluate(item["terms"], ["3", "0"]) for item in problem["inequalities"]]
        self.assertTrue(any(value < 0 for value in values))

    def test_reversed_or_duplicate_bounds_fail_closed(self) -> None:
        bounds = frozen_bounds()
        bounds[0]["lower"] = "3"
        with self.assertRaisesRegex(ValueError, "reversed"):
            compile_bounds(bounds)
        duplicate = frozen_bounds()
        duplicate[1]["name"] = duplicate[0]["name"]
        with self.assertRaisesRegex(ValueError, "unique"):
            compile_bounds(duplicate)

    def test_checked_in_result_reproduces(self) -> None:
        summary = validate_result(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")))
        self.assertEqual(summary["status"], "valid_canonical_parameter_bound_problem")


if __name__ == "__main__":
    unittest.main()
