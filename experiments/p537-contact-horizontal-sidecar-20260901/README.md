# P537 contact-horizontal sidecar

This is a read-only post-processor for the frozen N25 exact population and N65
contact-stage production.  `score.py` retains each N65 displacement's own
Schur coefficient `beta_d` and reports the two independent 100-batch delete-one
factor sets plus their complete cross-observable covariance.  `score_n25.py`
computes the same horizontal representative from the committed exact tables.

```bash
python3 experiments/p537-contact-horizontal-sidecar-20260901/score.py \
  --baseline results/server-20260828/P45-root-amplitude/n65.hist.csv \
  --tables results/p537-contact-stage-n65/shard-0.tsv \
           results/p537-contact-stage-n65/shard-1.tsv \
           results/p537-contact-stage-n65/shard-2.tsv \
           results/p537-contact-stage-n65/shard-3.tsv \
  --output results/p537-contact-horizontal-sidecar/N65.json

python3 experiments/p537-contact-horizontal-sidecar-20260901/score_n25.py \
  --axis results/p537-finite-collar/axis.csv \
  --tilted results/p537-finite-collar/tilted.csv \
  --aggregates results/p537-finite-collar/schur-aggregates.csv \
  --baseline-axis experiments/p537-landing-matrix-preflight-20260901/baseline-axis.csv \
  --baseline-tilted experiments/p537-landing-matrix-preflight-20260901/baseline-tilted.csv \
  --kernel experiments/p537-landing-matrix-preflight-20260901/kernel.tsv \
  --full-result experiments/p537-landing-matrix-preflight-20260901/result.json \
  --output results/p537-contact-horizontal-sidecar/N25.json
```
