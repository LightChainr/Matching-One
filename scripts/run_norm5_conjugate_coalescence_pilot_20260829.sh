#!/usr/bin/env bash
set -euo pipefail

BIN="${BIN:-build/threshold_rank_integer_period_mc}"
THREADS="${THREADS:-8}"
SAMPLES="${SAMPLES:-10000000}"
BATCHES="${BATCHES:-100}"
GIT_COMMIT="${GIT_COMMIT:-$(git rev-parse HEAD)}"
OUT="${OUT:-results/server-20260829/P205-norm5-conjugate-coalescence/raw}"
mkdir -p "$OUT"

run_pair() {
  local first_a="$1" first_b="$2" second_a="$3" second_b="$4"
  local seed="$5" offset="$6" name="$7"
  "$BIN" \
    --first-matrix "$first_a" "$((-first_b))" "$first_b" "$first_a" \
    --second-matrix "$second_a" "$((-second_b))" "$second_b" "$second_a" \
    --first-rep "$first_a" "$first_b" --second-rep "$second_a" "$second_b" \
    --samples "$SAMPLES" --batches "$BATCHES" --seed "$seed" \
    --replica-offset "$offset" --threads "$THREADS" --git-commit "$GIT_COMMIT" \
    --output-prefix "$OUT/$name"
}

# Same seed/counter inside each size is intentional: C is first in both jobs,
# so its threshold stream must reproduce exactly and supplies the covariance bridge.
run_pair 15 10 17 6 2026105501 9300000000 n325_C_A_10m &
pid_325_a=$!
run_pair 15 10 18 1 2026105501 9300000000 n325_C_B_10m &
pid_325_b=$!
wait "$pid_325_a" "$pid_325_b"

run_pair 20 5 16 13 2026105502 9300000000 n425_C_A_10m &
pid_425_a=$!
run_pair 20 5 19 8 2026105502 9300000000 n425_C_B_10m &
pid_425_b=$!
wait "$pid_425_a" "$pid_425_b"
