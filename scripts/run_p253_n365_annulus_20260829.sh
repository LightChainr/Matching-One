#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="results/server-20260829/P253-n365-annulus"
RAW_DIR="${OUT_DIR}/raw"
PREFIX="${RAW_DIR}/n365_200k"

mkdir -p build "${RAW_DIR}"

g++ -O3 -DNDEBUG -std=c++17 -fopenmp \
  src/matching_multiradius_pivotal_mc.cpp \
  -o build/matching_multiradius_pivotal_mc

build/matching_multiradius_pivotal_mc --validate-only \
  --radii 2,4,7,8 --cutoff euclidean \
  --design n365_first,14,13 --design n365_second,19,2

{
  date -u
  uname -a
  lscpu
  free -h
  g++ --version
  git rev-parse HEAD
  git status --short
} > "${OUT_DIR}/environment.txt"

/usr/bin/time -v -o "${RAW_DIR}/run.time.txt" \
  build/matching_multiradius_pivotal_mc \
  --samples 200000 --batches 200 --radii 2,4,7,8 \
  --cutoff euclidean --p 0.592746050790 --threads 16 \
  --seed 25336560829 --replica-offset 25336500000 \
  --design n365_first,14,13 --design n365_second,19,2 \
  --git-commit "$(git rev-parse HEAD)" \
  --output-prefix "${PREFIX}" \
  > "${RAW_DIR}/production.stdout.txt"

find "${RAW_DIR}" -maxdepth 1 -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "${OUT_DIR}/CHECKSUMS.sha256"

echo "completed frozen P253 N365 acquisition"
echo "metadata: ${PREFIX}.metadata.json"
echo "checksums: ${OUT_DIR}/CHECKSUMS.sha256"
