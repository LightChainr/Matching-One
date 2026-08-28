#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$ROOT/tests/exact_oracle.py" --L 2 3 \
  --compare "$ROOT/results/issue-7/L02_configs.csv" "$ROOT/results/issue-7/L03_configs.csv" \
  --json-out "$ROOT/results/issue-7/oracle_compare.json"
