# P263 boundary Q-score Phase-D smoke

Command:

```bash
/usr/bin/time -p build/p263/p263_boundary_qscore_pilot \
  --level 1 --samples 20000 --batches 20 --seed 202608290263 \
  --output results/local-20260829/P263-boundary-qscore-smoke/level1_20k.batches.csv

python3 scripts/score_p263_boundary_qscore_pilot.py \
  --batches results/local-20260829/P263-boundary-qscore-smoke/level1_20k.batches.csv \
  --output results/local-20260829/P263-boundary-qscore-smoke/score.json
```

- Wall time: 7.55 seconds on one local core.
- `14|23` event counts in lambda order `(1/4,1/3,2/3,3/4)`:
  `(19,19,133,210)`.
- Active residual: `(-3.4912085710,-3.8724309998,-4.5147204941)`.
- Joint diagnostic: `chi2=6.7421914932`, covariance rank 3.

The first two coordinates contain only 19 events.  This output freezes an
end-to-end deterministic regression and resource calibration; it is not a
continuum-mechanism score.
