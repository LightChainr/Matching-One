import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_intrinsic_functional_cocycle import (  # noqa: E402
    LEVELS,
    LINEAGES,
    SIZES,
    calculate,
    merge_inputs,
    parse_groups,
    residual_function,
)


class IntrinsicFunctionalCocycleTests(unittest.TestCase):
    def synthetic_sample(self, kind):
        sample = {}
        for n in SIZES:
            levels = {}
            for level in LEVELS:
                leading = 1.2 + 0.7 * level
                correction = -0.4 + 0.3 * level
                if kind == "q2":
                    scaled = leading + correction / n
                else:
                    scaled = leading + correction * math.log(n)
                levels[str(level)] = {"T": n ** (-13.0 / 8.0) * scaled}
            sample[n] = {"levels": levels}
        return sample

    def test_fixed_cocycles_cancel_unknown_functions(self):
        q2 = self.synthetic_sample("q2")
        jordan = self.synthetic_sample("jordan")
        for lineage in LINEAGES:
            for level in LEVELS:
                self.assertAlmostEqual(
                    residual_function(*lineage, level, 8.0 / 5.0)(q2), 0.0, places=12
                )
                self.assertAlmostEqual(
                    residual_function(*lineage, level, math.log(5.0) / math.log(2.0))(jordan),
                    0.0,
                    places=12,
                )

    def write_histogram(self, root, n, samples_per_batch):
        path = root / f"n{n}.hist.csv"
        fields = ["n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"]
        designs = {
            65: ((8, 1), (7, 4)), 85: ((9, 2), (7, 6)),
            130: ((9, 7), (11, 3)), 170: ((11, 7), (13, 1)),
            325: ((17, 6), (18, 1)), 425: ((16, 13), (19, 8)),
        }
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            for batch in range(4):
                for orientation, rep in zip(("first", "second"), designs[n]):
                    shift = (1 if orientation == "first" else -1) * (batch + 1)
                    for kind, base in (("minus", 0.54), ("plus", 0.57)):
                        rank = max(1, min(n, round(base * n) + shift))
                        writer.writerow({
                            "n": n, "a": rep[0], "b": rep[1], "orientation": orientation,
                            "batch": batch, "samples": samples_per_batch,
                            "kind": kind, "k": rank, "count": samples_per_batch,
                        })
        return path

    def test_end_to_end_independent_target_groups(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [
                self.write_histogram(root, n, 100 if n < 300 else 200)
                for n in SIZES
            ]
            groups = parse_groups(("65,85,130,170", "325", "425"))
            result = calculate(merge_inputs(paths), groups)
            self.assertEqual(result["covariance_groups"], [[65, 85, 130, 170], [325], [425]])
            self.assertEqual(result["model_order"], ["q2_analytic", "jordan_rank2"])
            self.assertEqual(list(result["models_in_frozen_order"]), ["q2_analytic", "jordan_rank2"])
            for model in result["models_in_frozen_order"].values():
                self.assertEqual(len(model["residual"]), 6)
                self.assertTrue(math.isfinite(model["chi_square"]))

    def test_rejects_nonpartition_groups(self):
        with self.assertRaisesRegex(ValueError, "partition"):
            parse_groups(("65,85,130,170", "325"))


if __name__ == "__main__":
    unittest.main()
