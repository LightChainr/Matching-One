#!/usr/bin/env python3
"""Integration tests for the C01 Gaussian orientation engine."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gaussian_circulant_geometry import GaussianTorus, best_multiplier  # noqa: E402
from torus_homology import HomologyUnionFind, wrapping_channels  # noqa: E402


class CirculantHomologyAgreementTests(unittest.TestCase):
    def test_circulant_n5_matches_c00_counts(self) -> None:
        torus = GaussianTorus(2, 1)
        n = torus.n
        primal = []
        for j in range(n):
            primal.append(type("E", (), {"i": j, "j": (j + torus.a) % n, "dx": 1, "dy": 0})())
            primal.append(type("E", (), {"i": j, "j": (j + torus.b) % n, "dx": 0, "dy": 1})())
        counts = {"rank0": 0, "rank1": 0, "rank2": 0, "d0": 0, "d1": 0}
        for mask in range(1 << n):
            active = [bool((mask >> v) & 1) for v in range(n)]
            uf = HomologyUnionFind(n, ((2, -1), (1, 2)))
            for edge in primal:
                if active[edge.i] and active[edge.j]:
                    uf.add_edge(edge.i, edge.j, edge.dx, edge.dy)
            channels = wrapping_channels(
                [uf.component(v) for v in range(n) if active[v]] or []
            )
            # Isolated empty mask: no active components.
            if not any(active):
                self.assertEqual(channels.max_rank, 0)
            counts[f"rank{channels.max_rank}"] += 1
            counts["d0"] += int(channels.direction_0)
            counts["d1"] += int(channels.direction_1)
        self.assertEqual(counts, {"rank0": 16, "rank1": 10, "rank2": 6, "d0": 11, "d1": 11})

    def test_structural_multipliers_are_units(self) -> None:
        pairs = [((8, 1), (7, 4)), ((9, 2), (7, 6)), ((12, 1), (9, 8))]
        expected = {65: 4, 85: 3, 145: 8}
        for r1, r2 in pairs:
            g1, g2 = GaussianTorus(*r1), GaussianTorus(*r2)
            t, _score, nn, _ma = best_multiplier(g1, g2)
            self.assertEqual(g1.n, g2.n)
            self.assertEqual(t, expected[g1.n])
            self.assertEqual(math_gcd(t, g1.n), 1)
            self.assertGreaterEqual(nn, 2)


def math_gcd(a: int, b: int) -> int:
    import math
    return math.gcd(a, b)


class GaussianOrientationMCTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler found")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / "gaussian_orientation_mc"
        command = [compiler, "-O3", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            command.append("-fopenmp")
        command += [str(ROOT / "src" / "gaussian_orientation_mc.cpp"), "-o", str(cls.binary)]
        subprocess.run(command, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_exact_self_test(self) -> None:
        completed = subprocess.run(
            [str(self.binary), "--self-test"], check=True, text=True, capture_output=True
        )
        self.assertIn("self-test passed", completed.stdout)

    def test_simulation_reproducible_and_analyzable(self) -> None:
        prefix_a = Path(self.temp.name) / "run_a"
        prefix_b = Path(self.temp.name) / "run_b"
        common = [
            "--rep1", "2,1", "--rep2", "2,1", "--t", "1",
            "--mode", "site", "--p", "0.5", "--samples", "400",
            "--batches", "4", "--seed", "17", "--replica-begin", "0",
            "--threads", "1",
        ]
        subprocess.run([str(self.binary), *common, "--output-prefix", str(prefix_a)], check=True)
        subprocess.run([str(self.binary), *common, "--output-prefix", str(prefix_b)], check=True)
        self.assertEqual(
            Path(f"{prefix_a}.t1.batches.csv").read_bytes(),
            Path(f"{prefix_b}.t1.batches.csv").read_bytes(),
        )
        analysis_json = Path(self.temp.name) / "analysis.json"
        subprocess.run(
            [
                sys.executable, str(ROOT / "scripts" / "analyze_gaussian_orientation.py"),
                "--moments", f"{prefix_a}.t1.moments.json",
                "--metadata", f"{prefix_a}.metadata.json",
                "--json", str(analysis_json),
            ],
            check=True,
        )
        result = json.loads(analysis_json.read_text(encoding="utf-8"))
        self.assertEqual(result["N"], 5)
        self.assertEqual(result["batches"], 4)
        self.assertIn("either", result["channels"])
        self.assertIn("M_delta", result["channels"]["either"])
        self.assertIn("S_delta", result["channels"]["either"])
        # Same orientation pair at t=1 must have zero orientation difference.
        for ch in result["channels"].values():
            self.assertAlmostEqual(ch["M_delta"]["mean"], 0.0, places=12)
            self.assertAlmostEqual(ch["S_delta"]["mean"], 0.0, places=12)

    def test_bond_mode_writes_moments(self) -> None:
        prefix = Path(self.temp.name) / "bond"
        subprocess.run(
            [
                str(self.binary),
                "--rep1", "2,1", "--rep2", "2,1", "--t", "1",
                "--mode", "bond", "--p", "0.5", "--samples", "200",
                "--batches", "4", "--seed", "19", "--threads", "1",
                "--output-prefix", str(prefix),
            ],
            check=True,
        )
        self.assertTrue(Path(f"{prefix}.t1.moments.json").exists())
        meta = json.loads(Path(f"{prefix}.metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta["mode"], "bond")
        self.assertEqual(meta["p"], 0.5)


if __name__ == "__main__":
    unittest.main()
