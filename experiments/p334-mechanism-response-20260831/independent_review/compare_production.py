#!/usr/bin/env python3
"""Compare production O(Q) moments to independent exact rational algebra."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from analyze_extension import fast_moments
from verify_ustat import det_estimators

rng = np.random.default_rng(20260831)
checks = []
for q in (4, 8, 64, 72):
    for case in range(10):
        integers = rng.integers(-5, 6, size=(q, 4))
        for denominator in (1, 1000, 1000000):
            rows = [tuple(Fraction(int(x), denominator) for x in row) for row in integers]
            exact = np.array([float(x) for x in det_estimators(rows)])
            x = np.array(rows, dtype=float).reshape(1, q, 1, 2, 2)
            actual = fast_moments(x)[0, :2]
            scales = np.array([denominator**-2, denominator**-4])
            scaled_error = np.abs(actual-exact)/scales
            assert np.all(scaled_error < 1e-10), (q, case, denominator, actual, exact)
            checks.append(float(scaled_error.max()))
report = {
    "checks": len(checks), "quartet_counts": [4,8,64,72],
    "entry_denominators": [1,1000,1000000],
    "max_error_after_rescaling_to_integer_input": max(checks),
    "production_sha256": hashlib.sha256((HERE.parent/"analyze_extension.py").read_bytes()).hexdigest(),
    "independent_algebra_sha256": hashlib.sha256((HERE/"verify_ustat.py").read_bytes()).hexdigest(),
    "scope": "Local independent algebra verification only; no new scientific samples or cloud run."
}
(HERE/"production_comparison.json").write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps(report))
