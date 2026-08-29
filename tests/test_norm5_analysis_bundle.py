import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from run_norm5_analysis_bundle import (  # noqa: E402
    REQUIRED_SIZES,
    RunSpec,
    build_commands,
    infer_covariance_groups,
    parse_run,
    validate_runs,
)


class Norm5AnalysisBundleTests(unittest.TestCase):
    def fake_runs(self):
        return [
            RunSpec(n, Path(f"n{n}.hist.csv"), Path(f"n{n}.moments.csv"), Path(f"n{n}.metadata.json"))
            for n in REQUIRED_SIZES
        ]

    def test_parse_run(self) -> None:
        run = parse_run("325:a.hist.csv:a.moments.csv:a.metadata.json")
        self.assertEqual(run.n, 325)
        self.assertEqual(run.histogram, Path("a.hist.csv"))
        self.assertEqual(run.moments, Path("a.moments.csv"))
        self.assertEqual(run.metadata, Path("a.metadata.json"))

    def test_validate_requires_exact_bundle_sizes_but_not_extra_process(self) -> None:
        by_n = validate_runs(self.fake_runs(), check_files=False)
        self.assertEqual(tuple(sorted(by_n)), REQUIRED_SIZES)
        with self.assertRaisesRegex(ValueError, "missing"):
            validate_runs(self.fake_runs()[:-1], check_files=False)

    def test_covariance_groups_are_inferred_from_exact_counter_domain(self) -> None:
        runs = self.fake_runs()
        metadata = {
            65: {"seed": 7, "replica_counter_first": 0, "replica_counter_last_exclusive": 100},
            130: {"seed": 7, "replica_counter_first": 0, "replica_counter_last_exclusive": 100},
            325: {"seed": 7, "replica_counter_first": 0, "replica_counter_last_exclusive": 100},
            85: {"seed": 9, "replica_counter_first": 0, "replica_counter_last_exclusive": 100},
            170: {"seed": 9, "replica_counter_first": 0, "replica_counter_last_exclusive": 100},
            425: {"seed": 12, "replica_counter_first": 500, "replica_counter_last_exclusive": 600},
        }
        self.assertEqual(
            infer_covariance_groups(runs, metadata),
            [[65, 130, 325], [85, 170], [425]],
        )

    def test_bundle_uses_existing_scorers_without_redefining_them(self) -> None:
        by_n = validate_runs(self.fake_runs(), check_files=False)
        groups = [[65, 130, 325], [85, 170], [425]]
        commands = build_commands(
            ROOT,
            by_n,
            groups,
            Path("bundle-out"),
            Path("source-rank-gap.json"),
        )
        names = [row[0] for row in commands]
        self.assertEqual(
            names,
            [
                "primary_harmonic",
                "intrinsic_functional_cocycle",
                "krawtchouk_score_modes",
                "rank_gap_boundary",
            ],
        )
        flattened = "\n".join(" ".join(command) for _, command, _ in commands)
        self.assertIn("score_norm5_harmonic_primary_typed.py", flattened)
        self.assertIn("score_intrinsic_functional_cocycle_typed.py", flattened)
        self.assertIn("threshold_score_modes.py", flattened)
        self.assertIn("score_rank_gap_boundary_targets.py", flattened)
        self.assertIn("65,130,325", flattened)
        self.assertIn("85,170", flattened)


if __name__ == "__main__":
    unittest.main()
