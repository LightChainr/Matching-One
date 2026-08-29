# Phase D: runnable boundary Q-score pilot

This phase turns the Phase C sufficient-statistics contract into one deliberately
small square-bond program.  It is a calibration runner, not a generic FK
framework and not evidence for the continuum ODE by itself.

## Exact lattice observable

For each open rectangle at (p=1/2), the runner records the three mutually
exclusive boundary connectivity indicators

\[
 (I_{1234},I_{12|34},I_{14|23})
\]

together with (J=2k+b).  Thus every batch consists only of integer counts,
(\sum J), (\sum J^2), and their three indicator-weighted versions.  At
(Q=1), the exact finite-graph identity

\[
 \partial_Q P_Q(A)|_{Q=1}=\operatorname{Cov}_1(I_A,J/2)
\]

reconstructs the measure-score tangent without simulating a second value of
(Q).  The event classifier is exactly the one validated by the exhaustive
four-cycle oracle in Phase C.

## Exact-rational four-geometry family

Take bottom-boundary vertices (0,x_2,1,2) after translating away a padding of
(2s).  The four integer realizations are

| lambda | x2 | base span s | conformal K |
|---|---:|---:|---:|
| 1/4 | 2/5 | 15 | 10/3 |
| 1/3 | 1/2 | 14 | 3 |
| 2/3 | 4/5 | 15 | 15/4 |
| 3/4 | 6/7 | 14 | 14/3 |

The box has ((6s+1)(4s+1)) vertices.  Level 2 doubles every span.  This keeps
the cross ratios exact; the price is a known 14-versus-15 lattice-scale
mismatch.  Free top and side boundaries are a finite-box approximation to the
upper half-plane, so the two levels must be scored separately.

## Frozen score

For the `14|23` event, each geometry produces

\[
 z_i=\frac{\operatorname{Cov}(I_i,J/2)}{P_i}
     +4h'\log s_i-2h'\log K_i,
 \qquad h'=\frac{\sqrt3}{3\pi}.
\]

Subtract (z_{1/3}) and the Phase C continuum target.  The active residual is
the three-vector at (lambda=(1/4,2/3,3/4)); synchronized delete-one batches
give its full covariance and (r^T\Sigma^+r).  No individual coordinate is
treated as independent evidence.

## Frozen acquisition and commands

The machine-readable contract is
`experiments/p263_boundary_qscore_phaseD_20260829.yaml`.  The two acquisitions
are 200,000 samples per geometry at level 1 and 500,000 at level 2, each in 100
batches.  The 20,000-sample calibration took 7.55 seconds on one local core;
linear scaling predicts about 76 seconds and 12.6 minutes respectively.  The
largest graph has 21,901 vertices and 43,500 edges, so memory is negligible.

```bash
c++ -O2 -std=c++17 src/p263_boundary_qscore_pilot.cpp \
  -o build/p263/p263_boundary_qscore_pilot

build/p263/p263_boundary_qscore_pilot \
  --level 1 --samples 200000 --batches 100 --seed 2026102631 \
  --output results/p263/level1.batches.csv

python3 scripts/score_p263_boundary_qscore_pilot.py \
  --batches results/p263/level1.batches.csv \
  --output results/p263/level1.score.json
```

Repeat with level 2, 500,000 samples, and seed `2026102632`.

## Smoke result

The committed level-1 smoke used 20,000 samples per geometry in 20 batches.
The rare `14|23` counts were `(19,19,133,210)`.  Its active residual was

\[
 (-3.4912,-3.8724,-4.5147),\qquad \chi^2_3=6.74.
\]

With only 19 events in the first two geometries, this is a software and
sensitivity smoke only.  It demonstrates that the full integer-statistics,
jackknife-covariance, amplitude projection, and GLS path runs end to end; it
does not score the continuum mechanism.

## Claim layers

- Exact finite graph: event semantics and the (J/2) covariance tangent.
- Frozen continuum input: the Phase C high-branch ODE target.
- Scaling hypothesis: these finite rectangles approach the four upper-half-
  plane connectivities after the stated field and conformal normalization.
