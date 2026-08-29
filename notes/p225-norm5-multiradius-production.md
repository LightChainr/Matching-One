# Issue #225: first scientific norm-5 multi-radius job

This branch converts the PR #240 prototype into a production-ready square-site
matching job.  No 200k production samples were generated locally.

## Why the checkerboard engine cannot simply take the norm-5 reps

PR #240 inherits the checkerboard C4 self-matching control, which requires
both Gaussian coordinates to be odd.  Its quotients therefore have even
order.  Every N325/N425 representation has one even and one odd coordinate,
so feeding `(17,6)`, `(18,1)`, `(16,13)` or `(19,8)` into that geometry would
change or invalidate the graph.

The new engine instead uses the PR #211 semantics:

```text
black field       -> NN primal cross pivotal
complemented white -> NN+NNN matching cross pivotal.
```

The global pivotal decision is still made once per graph and reused for every
radius.

## Cutoff gate

The Chebyshev square at R8 aliases `(16,13)` on N425: two patch points can
differ by that exact period.  This is a mathematical obstruction, not an
implementation inconvenience.  The production job therefore freezes a
Euclidean disk

```text
0 < x^2+y^2 <= R^2
```

with the outer unit shell carrying the eight-sector landing registry.  It is
C4/reflection symmetric, and startup construction verifies injectivity for
every design and radius before any counter is consumed.  R2/R4/R8 pass for all
four norm-5 geometries.  This disk observable is new and is not described as
a replay of the PR #211 square R3 mark.

## Frozen orientation order

The contrasts are exactly

```text
N325: (17,6) - (18,1), Delta cos(4 theta) = -16128/21125
N425: (16,13) - (19,8), Delta cos(4 theta) = -32256/36125.
```

Both angular factors are negative, so no H4 registry reversal occurs between
sizes.  The landing mark remains absolute lattice axis minus diagonal for all
four designs.

## Two distinct scale views

The same production block yields:

1. fixed-R UV contrasts at R2/R4/R8;
2. dyadic per-log shells 2->4 and 4->8;
3. a near-fixed-delta comparison, N325/R7 against N425/R8.

The last pair has

```text
delta325 = 7/sqrt(325) = 0.3882901...
delta425 = 8/sqrt(425) = 0.3880570...
relative mismatch = 6.006e-4.
```

The scorer recomputes every conditional amplitude, orientation contrast,
shell and matched-delta difference inside each of 200 delete-one batches.  It
emits the full 16-coordinate contrast covariance and 8-coordinate shell
covariance.  All outputs are correlated views of one raw block.

## Huawei execution

After checking out commit `f8e02c9` or this branch:

```bash
mkdir -p build results/server-20260829/P225-norm5-multiradius/raw

g++ -O3 -DNDEBUG -std=c++17 -fopenmp \
  src/matching_multiradius_pivotal_mc.cpp \
  -o build/matching_multiradius_pivotal_mc

build/matching_multiradius_pivotal_mc \
  --validate-only --radii 2,4,7,8 --cutoff euclidean

/usr/bin/time -v build/matching_multiradius_pivotal_mc \
  --samples 200000 --batches 200 --radii 2,4,7,8 \
  --cutoff euclidean --p 0.592746050790 --threads 16 \
  --seed 22550260829 --replica-offset 15000000000 \
  --git-commit f8e02c9 \
  --output-prefix results/server-20260829/P225-norm5-multiradius/raw/norm5_200k \
  2> results/server-20260829/P225-norm5-multiradius/raw/run.time.txt

python3 scripts/analyze_norm5_multiradius_pivotal.py \
  --batches results/server-20260829/P225-norm5-multiradius/raw/norm5_200k.batches.csv \
  --metadata results/server-20260829/P225-norm5-multiradius/raw/norm5_200k.metadata.json \
  --output results/server-20260829/P225-norm5-multiradius/analysis.json
```

The counter interval `[15000000000,15000200000)` is fresh relative to the
PR #211 12-billion interval and PR #224's zero-based stream.
