#!/usr/bin/env bash
set -euo pipefail
task_size="${1:?N325 or N425}"
task_mode="${2:?smoke or production}"
task_commit="${3:?frozen commit}"
task_samples=20000
task_batches=20
case "$task_size" in
  N325)
    task_seed=20260831430325; task_offset=43032500000; task_k0=193
    task_geometry=(--first-matrix 325 57 0 1 --second-matrix 325 18 0 1 --first-rep 17 6 --second-rep 18 1)
    ;;
  N425)
    task_seed=20260831430425; task_offset=43042500000; task_k0=252
    task_geometry=(--first-matrix 425 132 0 1 --second-matrix 425 268 0 1 --first-rep 16 13 --second-rep 19 8)
    ;;
  *) exit 2 ;;
esac
if [[ "$task_mode" == smoke ]]; then
  task_samples=100; task_batches=4
  task_seed=$((task_seed+1)); task_offset=$((task_offset+100000))
elif [[ "$task_mode" != production ]]; then
  exit 2
fi
task_prefix="raw/${task_size}-${task_mode}"
mkdir -p raw
if [[ -e "${task_prefix}.metadata.json" || -e "${task_prefix}.exit" ]]; then
  echo "Refusing to overwrite completed or attempted run ${task_prefix}" >&2
  exit 3
fi
trap 'task_exit=$?; printf "%s\n" "$task_exit" > "${task_prefix}.exit"' EXIT
/usr/bin/time -v ./threshold_rank_integer_period_mc \
  "${task_geometry[@]}" --samples "$task_samples" --batches "$task_batches" \
  --seed "$task_seed" --replica-offset "$task_offset" --threads 16 \
  --geometry-pilot-k0 "$task_k0" --branching-clones --cooperative-closure \
  --git-commit "$task_commit" --output-prefix "$task_prefix"
