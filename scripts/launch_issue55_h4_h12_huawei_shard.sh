#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 SHARD_INDEX [OUTPUT_DIR]" >&2
  exit 2
fi

readonly SHARD_INDEX="$1"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly OUTPUT_DIR="${2:-${REPO_ROOT}/results/server-20260830/P55-h4-h12-orthogonal-production/shard${SHARD_INDEX}}"
readonly BINARY="${REPO_ROOT}/build/threshold_rank_integer_period_mc_issue55"
readonly SOURCE="${REPO_ROOT}/src/threshold_rank_integer_period_mc.cpp"
readonly FULL_COMMIT="$(git -C "${REPO_ROOT}" rev-parse HEAD)"

if ! git -C "${REPO_ROOT}" diff --quiet || ! git -C "${REPO_ROOT}" diff --cached --quiet; then
  echo "tracked worktree changes present; refusing a production launch" >&2
  exit 1
fi

mkdir -p "${REPO_ROOT}/build" "${OUTPUT_DIR}"
g++ -O3 -std=c++17 -fopenmp "${SOURCE}" -o "${BINARY}"
"${BINARY}" --self-test
sha256sum "${BINARY}" | tee "${OUTPUT_DIR}/binary.sha256"

python3 "${REPO_ROOT}/scripts/run_issue55_h4_h12_acquisition.py" \
  --stage production \
  --shard-index "${SHARD_INDEX}" \
  --binary "${BINARY}" \
  --output-dir "${OUTPUT_DIR}" \
  --threads 8 \
  --git-commit "${FULL_COMMIT}" \
  --parallel-sizes
