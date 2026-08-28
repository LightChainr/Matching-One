from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class NewmanZiffThresholdTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        compiler = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
        if compiler is None:
            raise unittest.SkipTest("no C++17 compiler")
        cls.temp = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.temp.name) / "newman_ziff_cpu"
        cmd = [compiler, "-O2", "-std=c++17"]
        if sys.platform != "darwin" and "clang" not in Path(compiler).name:
            cmd.append("-fopenmp")
        cmd += [str(ROOT / "src" / "newman_ziff_cpu.cpp"), "-o", str(cls.binary)]
        subprocess.run(cmd, check=True, cwd=ROOT)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def test_exact_tiny_systems_pass(self) -> None:
        out = Path(self.temp.name) / "exact_out"
        completed = subprocess.run(
            [str(self.binary), "--exact-tests", "--outdir", str(out)],
            check=True, text=True, capture_output=True,
        )
        self.assertIn("exact tests PASS", completed.stderr)
        payload = json.loads((out / "exact" / "exact_tests.json").read_text())
        self.assertEqual(payload["overall"], "PASS")
        self.assertEqual(payload["exhaustive"]["axis_L2"]["Kminus_gt_Kplus"], 0)
        self.assertEqual(payload["exhaustive"]["axis_L3"]["Kminus_gt_Kplus"], 0)


if __name__ == "__main__":
    unittest.main()
