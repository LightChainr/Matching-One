#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mkdir -p "$ROOT/results/issue-7/logs"
/usr/bin/time -v "$ROOT/build/exact_torus_matching" --L 5 --threads 8 --outdir "$ROOT/results/issue-7" \
  > "$ROOT/results/issue-7/logs/L05.stdout.txt" \
  2> "$ROOT/results/issue-7/logs/L05.time.txt"
