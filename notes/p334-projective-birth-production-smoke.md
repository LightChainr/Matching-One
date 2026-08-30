# Projective essential-birth spectroscopy: Phase B/C prototype

Status: production-capable sparse event prototype plus a local N=65
variance/runtime smoke for Issue #334.  This branch starts from the Phase A
oracle `6f54935` and inherits the unrestricted saturation theorem `c1a72e5`.
It does not start a Huawei production run and does not choose a model from the
pilot mean.

## Minimal stream contract

The general integer-period threshold-rank engine now has an opt-in
`--projective-births` output.  For each permutation it extracts

```text
tau1, primitive projective ell1 up to sign, tau2
```

or the typed atom

```text
DIRECT_RANK2, with tau1=tau2 and no ell.
```

The old `ell1` is retained until the second birth.  There is no `iota` field:
by `c1a72e5` the rank-one carrier is saturated and `iota=1`.  The raw archive
is a sparse joint histogram by aligned batch rather than one row per replica;
it is therefore usable at production scale while retaining every sufficient
event coordinate needed for arbitrary fixed-p scoring.

Every sample is checked against the independent reverse matching filtration:

```text
tau1 == K_minus,
tau2 == K_plus.
```

The exact self-test also reproduces all 120 N=5 paths and the eight
`DIRECT_RANK2` paths of the axis L=2 control.

## N=65 smoke design

The smaller established Gaussian orientation pair was chosen:

```text
first  P = [[8,-1],[1,8]],
second P = [[7,-4],[4,7]],
N=65,
20,000 fresh counter permutations per shape,
20 aligned batches,
seed 20260830, counters [0,20000).
```

Both shapes consume the same counter permutation.  This local one-thread run
took 0.160 s inside the engine (1.05 s end-to-end under `/usr/bin/time`) and
about 3.5 MB maximum resident memory.  The event archive is 676 KB.

## Exact crosswalks

All machine gates pass, with maximum double-precision residuals below
`2.2e-14`:

1. summing all line-bearing births plus `DIRECT_RANK2` reconstructs the K1
   histogram in every batch; tau2 reconstructs K2;
2. the fixed-p plateau character
   `A4=E[1{tau1<=K<tau2} chi4(ell1)]` equals the Issue #156 primitive-sector
   character (maximum residual `6.7e-16`);
3. coefficient scoring gives
   `dA4/dp=j4_birth1-j4_exit2` (maximum residual `6.8e-15`);
4. the trivial-character current obeys
   `M'=line activity + 2*j_DIRECT_RANK2` (maximum residual `2.2e-14`).

There were 785 direct births across the two 20k shape streams.  Those births
cannot be assigned a projective character and remain an explicit atom.

## Does ell add non-micro spin-four support?

Yes, at the representation/support level.  Unlike the tiny Phase A controls,
each N=65 shape has two well-populated inequivalent spin-four values: the
axis-period lines and diagonal-period lines have opposite `chi4`.  Four
distinct values were observed in total once rare longer primitive lines are
included.  The empirical within-path character variance is 0.474 on each
shape; birth-flux weighting raises it to 0.531 and 0.518.  Thus the mark is not
a near-deterministic relabelling of `(K1,K2)` on a non-micro quotient.

The source/sink coordinates are also cheaply measurable.  At
`p_ref=0.592746050790`, the real-part orientation contrasts have batch SEs

```text
j4_birth1(second-first):    0.0236
j4_exit2(second-first):     0.0223
j4_activity(second-first):  0.0414
```

The corresponding pilot means are intentionally not used to select a radial
law, continuum field, or phase convention.  Their purpose is to show that the
new source/sink split has finite variance and that one common covariance block
can be produced without a large run.  `A4` remains a duplicate coordinate of
#156; only the birth/exit decomposition is scientifically new.

## Reproduce

```bash
c++ -O3 -std=c++17 src/threshold_rank_integer_period_mc.cpp \
  -o /tmp/projective_birth_mc
/tmp/projective_birth_mc \
  --samples 20000 --batches 20 --seed 20260830 --replica-offset 0 \
  --threads 1 \
  --first-matrix 8 -1 1 8 --second-matrix 7 -4 4 7 \
  --first-rep 8 1 --second-rep 7 4 \
  --projective-births --git-commit 6f54935 \
  --output-prefix results/local-20260830/P334-projective-birth-N65-smoke/n65_20k

python3 scripts/analyze_projective_birth_smoke.py \
  --births results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.births.csv \
  --histogram results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.hist.csv \
  --metadata results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.metadata.json \
  --json results/local-20260830/P334-projective-birth-N65-smoke/score.json \
  --markdown results/local-20260830/P334-projective-birth-N65-smoke/score.md

python3 -m unittest discover -s tests \
  -p 'test_threshold_rank_integer_period_mc.py'
python3 -m unittest discover -s tests \
  -p 'test_projective_essential_birth_oracle.py'
```

