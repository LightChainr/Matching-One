#!/usr/bin/env bash
set -euo pipefail

output_root=${1:-results/server-20260829/P154-norm4-production}
commit_sha=$(git rev-parse HEAD)
mkdir -p "$output_root/raw"

run_pair() {
    local label=$1
    shift
    local separator=0
    local -a first=()
    local -a second=()
    for argument in "$@"; do
        if [[ $argument == --PAIR-- ]]; then
            separator=1
        elif (( separator == 0 )); then
            first+=("$argument")
        else
            second+=("$argument")
        fi
    done
    "${first[@]}" &
    local first_pid=$!
    "${second[@]}" &
    local second_pid=$!
    local first_status=0
    local second_status=0
    wait "$first_pid" || first_status=$?
    wait "$second_pid" || second_status=$?
    if (( first_status != 0 || second_status != 0 )); then
        echo "$label failed: first=$first_status second=$second_status" >&2
        return 1
    fi
    echo "$label complete"
}

run_pair targets \
    build/threshold_rank_integer_period_mc --n 260 --samples 1000000000 \
    --batches 100 --seed 2026105401 --replica-offset 8200000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n260_1b" \
    --PAIR-- \
    build/threshold_rank_integer_period_mc --n 340 --samples 1000000000 \
    --batches 100 --seed 2026105402 --replica-offset 8200000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n340_1b"

run_pair source_65_130 \
    build/threshold_rank_orientation_mc --n 65 --samples 1900000000 \
    --batches 100 --seed 2026104501 --replica-offset 5100000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n65_1900m" \
    --PAIR-- \
    build/threshold_rank_orientation_mc --n 130 --samples 1900000000 \
    --batches 100 --seed 2026104501 --replica-offset 5100000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n130_1900m"

run_pair source_85_170 \
    build/threshold_rank_orientation_mc --n 85 --samples 1900000000 \
    --batches 100 --seed 2026104501 --replica-offset 5100000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n85_1900m" \
    --PAIR-- \
    build/threshold_rank_orientation_mc --n 170 --samples 1900000000 \
    --batches 100 --seed 2026104501 --replica-offset 5100000000 \
    --threads 8 --git-commit "$commit_sha" \
    --output-prefix "$output_root/raw/n170_1900m"

echo "norm-4 production complete at $commit_sha"
