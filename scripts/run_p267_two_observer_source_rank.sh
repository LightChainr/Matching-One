#!/usr/bin/env bash
set -euo pipefail

# Frozen Huawei launcher for experiments/p267_two_observer_source_rank_20260830.yaml.
# Usage:
#   TARGET=N325 PREREG_COMMIT=<sha> bash scripts/run_p267_two_observer_source_rank.sh
#   TARGET=N425 PREREG_COMMIT=<sha> bash scripts/run_p267_two_observer_source_rank.sh

: "${TARGET:?set TARGET=N325 or TARGET=N425}"
: "${PREREG_COMMIT:?set PREREG_COMMIT to the frozen prereveal commit}"
P267_OUTPUT_ROOT="${P267_OUTPUT_ROOT:-/workspace/p267-two-observer-source-rank}"
P267_THREADS="${P267_THREADS:-16}"

current_commit="$(git rev-parse HEAD)"
if [[ "${current_commit}" != "${PREREG_COMMIT}" ]]; then
    echo "refusing non-preregistered checkout: ${current_commit} != ${PREREG_COMMIT}" >&2
    exit 2
fi

build_dir="${P267_OUTPUT_ROOT}/build-${PREREG_COMMIT}"
mkdir -p "${build_dir}"
binary="${build_dir}/threshold_rank_integer_period_mc"
g++ -O3 -std=c++17 -fopenmp src/threshold_rank_integer_period_mc.cpp -o "${binary}"
"${binary}" --self-test

case "${TARGET}" in
    N325)
        run_id=N325_2m
        seed=202608303252
        offset=19000000000
        matrix_args=(--first-matrix 17 -6 6 17 --second-matrix 18 -1 1 18)
        rep_args=(--first-rep 17 6 --second-rep 18 1)
        ;;
    N425)
        run_id=N425_2m
        seed=202608304252
        offset=21000000000
        matrix_args=(--first-matrix 16 -13 13 16 --second-matrix 19 -8 8 19)
        rep_args=(--first-rep 16 13 --second-rep 19 8)
        ;;
    *)
        echo "TARGET must be N325 or N425" >&2
        exit 2
        ;;
esac

output_dir="${P267_OUTPUT_ROOT}/${run_id}"
mkdir -p "${output_dir}"
sha256sum "${binary}" > "${output_dir}/BINARY_SHA256.txt"
"${binary}" \
    --samples 2000000 \
    --batches 100 \
    --threads "${P267_THREADS}" \
    "${matrix_args[@]}" \
    "${rep_args[@]}" \
    --seed "${seed}" \
    --replica-offset "${offset}" \
    --git-commit "${PREREG_COMMIT}" \
    --marked-births \
    --far-radius 6 \
    --output-prefix "${output_dir}/${run_id}" \
    > "${output_dir}/run.stdout" \
    2> "${output_dir}/run.stderr"

(
    cd "${output_dir}"
    sha256sum "${run_id}".* BINARY_SHA256.txt run.stdout run.stderr > SHA256SUMS.txt
)
