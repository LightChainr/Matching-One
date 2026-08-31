#!/usr/bin/env python3
"""Compare completed tables; this file does not compute either kernel."""
import csv
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
ours = root / "BELL8_DIAGRAM_RESULTS.csv"
other = root.parent / "results" / "kernel.csv"

def read(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    table = {row["partition"]: row for row in rows}
    assert len(table) == len(rows) == 4140
    return table

a, b = read(ours), read(other)
assert set(a) == set(b)
errors = []
for pi in a:
    x, y = a[pi], b[pi]
    for k, v in (("exterior_blocks", "components"),
                 ("shared_components", "shared"),
                 ("derivative_numerator16", "activation16")):
        if int(x[k]) != int(y[v]):
            errors.append(dict(partition=pi, field=k, ours=x[k], other=y[v]))
result = {"all_4140_partitions_exactly_equal": not errors,
          "compared_fields": ["exterior_blocks", "shared_components", "derivative_numerator16"],
          "errors": errors,
          "equality_diagram_table_sha256": hashlib.sha256(ours.read_bytes()).hexdigest(),
          "coarsening_table_sha256": hashlib.sha256(other.read_bytes()).hexdigest(),
          "other_table_path": "../results/kernel.csv",
          "main_algorithm_imported_other_code": False}
(root / "COMPARISON.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
assert not errors
