# Issue #225: executable two-cutoff pivotal prototype

## Outcome

The PR #224 engine admits the multi-radius extension without a new topology
kernel.  Global crossing and fixed-root pivotal status do not depend on the
landing radius.  The new stream therefore evaluates each black/white pivotal
flag once per configuration, removes the root once, and reads every requested
`LocalLanding(R)` mark from that same environment.

The runtime interface is

```text
--radii 1,2,4
```

and the long batch schema stores `black_h4`, `white_h4`, `h4_plus`,
`h4_minus`, their `(S_t,S_lambda)` products, and the common pivotal totals.
All sizes, radii and channels retain the batch/counter alignment needed for a
single full covariance calculation.

## Fixed R is not fixed delta

The physical second cutoff is

```text
delta = R/sqrt(N).
```

Holding `R` fixed sends `delta` to zero together with the lattice spacing and
measures the UV mark used in PRs #211/#224.  Holding `delta` fixed requires a
size-dependent integer radius `R(N)=round(delta sqrt(N))`.  The analyzer
records both `R` and `delta`; it never silently identifies these limits.

For nested radii it also reports

```text
[A(2R)-A(R)] / log(2)
```

for plus and minus channels, recomputed inside every delete-one batch.  These
shells are correlated coordinates of one block, not separate evidence votes.

## Exact and regression controls

- The inherited complete N10/R1 local/global response oracle still passes.
- A single-radius `--radii 3` stream is integer-identical to the frozen PR #224
  R3 stream for the local/global readouts and all score products.
- Runtime radii must be positive, unique and increasing.
- Annulus injectivity is checked before sampling.  `R=8` fails on both present
  tori, so the honest local dyadic window is `R={1,2,4}`.  A `2,4,8` test needs
  a larger quotient rather than periodic aliasing.

## 20k engineering smoke

The committed smoke uses N130/N170, 100 aligned batches and radii 1,2,4.  It
completed in about 0.62 engine seconds locally.  The conditional amplitudes
are:

| N | R | delta | A_plus | A_minus |
|---:|---:|---:|---:|---:|
| 130 | 1 | 0.0877 | -0.1076 +/- 0.0127 | +0.0098 +/- 0.0098 |
| 130 | 2 | 0.1754 | -0.3089 +/- 0.0134 | +0.0008 +/- 0.0124 |
| 130 | 4 | 0.3508 | -0.0281 +/- 0.0126 | +0.0004 +/- 0.0135 |
| 170 | 1 | 0.0767 | -0.1300 +/- 0.0111 | +0.0082 +/- 0.0112 |
| 170 | 2 | 0.1534 | -0.3193 +/- 0.0138 | +0.0194 +/- 0.0148 |
| 170 | 4 | 0.3068 | +0.0067 +/- 0.0163 | +0.0241 +/- 0.0159 |

The plus shells are strongly nonconstant and reverse sign between `1->2` and
`2->4`; the odd/minus shells are unresolved at 20k.  This is useful design
information: the prototype works and the current radii are not evidence for a
constant logarithmic shell flow.  No Jordan, exponent or asymptotic claim is
made from the smoke.

The near agreement of same-R amplitudes across N130/N170 is a UV-control
observation.  Because their deltas differ, it is not a fixed-delta transport
test.

## Next smallest scientific run

Use larger tori that admit `R={2,4,8}` injectively and choose at least one
cross-size pair with approximately matched `R/sqrt(N)`.  Freeze both views:

1. same `R` across sizes for UV collapse;
2. matched `delta` across sizes for the two-cutoff limit.

Then score the joint plus/minus per-log shell vector with the full aligned
covariance.  Increasing samples on N130/N170 alone cannot create the missing
mesoscopic window.

## Reproduction

```bash
c++ -O3 -std=c++17 src/c4_multiradius_pivotal_mc.cpp \
  -o build/c4_multiradius_pivotal_mc

build/c4_multiradius_pivotal_mc \
  --samples 20000 --batches 100 --radii 1,2,4 \
  --seed 22520260829 --replica-offset 0 --threads 1 \
  --git-commit c6d9768 \
  --output-prefix results/local-20260829/P225-multiradius-pivotal/raw/n130_n170_20k

python3 scripts/analyze_c4_multiradius_pivotal.py \
  --batches results/local-20260829/P225-multiradius-pivotal/raw/n130_n170_20k.batches.csv \
  --metadata results/local-20260829/P225-multiradius-pivotal/raw/n130_n170_20k.metadata.json \
  --output results/local-20260829/P225-multiradius-pivotal/analysis.json
```
