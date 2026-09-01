#!/bin/bash
set -euo pipefail
cd /workspace/p154-prospective-transmission-20260831
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=14 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/workspace/matching-one-p334-finite-20260831/.venv/bin/python run_production.py --n 340 --batch-begin 100 --batch-end 125 --workers 14 --authorization authorization.json
date -u > QUEUE_COMPLETED_UTC.txt
