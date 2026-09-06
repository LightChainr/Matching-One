# The kept direction has a law: A(N) ~ N^-0.97, and its own curvature falsifies it

**Date:** 2026-09-06
**Ticket:** #584 (after Gate 3)
**Claim level:** C2 — reanalysis of committed productions, no new samples
**Artifact:** `results/p582-amplitude-law/latest.json`
**Script:** `scripts/p582_amplitude_law.py` (25 tests in `tests/test_p582_amplitude_law.py`)

Gate 3 (`notes/type582-residual-20260906.md`) screened #582's structured
remainder against every pre-existing discrete label, found that none beats the
exact permutation null, and closed on the decision-table branch *keep the
dominant transferable direction as the robust finite object*.

That branch was taken. This note is about what was kept.

## The question Gate 3 did not ask

Gate 3 reports the five amplitudes of the frozen direction — `-2.18, -1.00,
-1.71, -0.78, -1.05` (x10^-3) — as inputs to a label screen. They are also five
numbers spanning a factor of 2.8, each measured to better than 0.3% (`z` from
335 to 659). Nobody had asked whether they obey a law.

They do, and it is a one-parameter one.

## Not the move #584 forbids

#584 says: *do not fit a free exponent to five residual directions.* Nothing
here fits a direction. `g` is frozen first, by the same power iteration Gate 3
uses, and one exponent is then fitted to the five **amplitudes** along that
already-frozen direction. The direction is the object #584 says to keep; its
amplitude law is a property of it, not a second latent coordinate. The residual
that #584 rules out explaining is left exactly where Gate 3 left it.

## The model, and the factor that is not optional

If the shape part of the law carries a single correction-to-scaling term,

```text
Q(s, u) = A(u) + lambda * exp(-omega * s) * g(u) + ...,     s = log N,
```

then what #582 measures is not the derivative but its exact finite-difference
image,

```text
v = [Q(s_t) - Q(s_b)] / h
  = -lambda * omega * exp(-omega * sbar) * sinh(omega*h/2)/(omega*h/2) * g.
```

The `sinh(x)/x` factor is 1.9% at `h = log 2` and 3.3% at `h = log 2.5`. It is
step-size dependent and the same size as the misfit being measured, so dropping
it would manufacture a difference between the `m=2` and `m=2.5` transitions —
which is precisely the `multiplier` label Gate 3 screened. It is carried
exactly.

## What the five amplitudes say

```text
transition     sbar       h      amplitude          se        z
   65->130  4.52096 0.69315  -2.184054e-03  3.317e-06  -658.5
  130->325  5.32568 0.91629  -9.990715e-04  2.278e-06  -438.5
   85->170  4.78922 0.69315  -1.709832e-03  3.162e-06  -540.7
  170->425  5.59394 0.91629  -7.816930e-04  1.992e-06  -392.5
  145->290  5.32331 0.69315  -1.045978e-03  3.127e-06  -334.5
```

One exponent: `omega = 0.9702`, `lambda = 1.7835e-1`, `chi^2 = 277.4` on 3 df.

The fit is **rejected** — and the misfit is at most 3.8% in amplitude, against
amplitudes measured to 0.3%. So a single exponent describes three lineages and
two step sizes to a few percent and is still statistically excluded. Both halves
of that sentence matter.

Leave-one-transition-out, the exponent never sees the transition it predicts:

```text
held out     omega on the other four   predicted vs measured
  65->130            0.9893                    -1.98%
 130->325            0.9585                    -2.91%
  85->170            0.9671                    +1.46%
 170->425            0.9647                    -1.01%
 145->290            0.9805                    +4.65%
```

A new transition's amplitude is predicted to within 5% by an exponent fitted
without it. The worst case is `p50`, the singleton lineage, which is also the
one with no third size.

## The error budget, and why `omega = 1` survives

The fit's own curvature is not the uncertainty on `omega`, because the fit is
rejected. The honest budget is

```text
omega                                   0.970
leave-one-out half spread              +-0.015
orientation weighting shift            +-0.027
combined                               +-0.031
```

The weighting term is real: rerunning everything under the naive equal-orientation
weighting that #582 warns about gives `omega = 0.997` and moves the frozen
direction by 7.8 degrees. Unity is 0.96 combined bars away. **`omega = 1` is not
excluded.**

In the site count, `N^-1` is `L^-2` on a square torus. That is a coincidence of
numbers until something independent fixes the operator; it is not a measurement
of any named percolation correction-to-scaling exponent, and this note does not
claim one.

## The independent check, and it fired

GOVERNANCE section 2 minimum A asks for one check by independent means. A second
divided difference across a three-size lineage is a different functional of the
same productions. Applied to the model it is `sum_k c_k exp(-omega s_k)` with
weights that annihilate constants. The parameters fitted on the five **first**
differences were frozen and asked to predict it — no refit:

```text
lineage       measured        predicted       ratio
gaussian_13  +2.24453e-03   +1.46051e-03      1.5368
gaussian_17  +1.75342e-03   +1.12583e-03      1.5574
p50           two sizes, no second difference
```

The one-exponent law misses the curvature by 55%, and it misses it **by the same
amount in both lineages**: the two ratios agree with each other to 1.3%, far
more closely than either agrees with one. Different parent primes (13 and 17),
different productions, different seeds. This is a reproducible falsification,
not noise.

## What this changes

The three results now line up:

```text
#582   a dominant transferable direction, plus a small structured remainder
Gate 3 no pre-existing discrete label indexes that remainder
here   the dominant direction's amplitude obeys one exponent to a few percent,
       and that same exponent is wrong by 55% on the curvature, twice over
```

Read together, the natural description of the remainder is **a second smooth
scale**, not a discrete fiber — which is consistent with Gate 3 finding nothing
to label, and is a stronger statement than "structured but unindexed". It is not
proof: a step-size-dependent bias in the quantile reconstruction that the first
differences happen to cancel would look the same. Exactly one measurement
separates those.

## The one production that settles it

`N = 725 = 5^2 * 29`, both orientations `(26, 7)` and `(23, 14)`.

Chosen by declared criteria, not preference. It does three things at once:

1. **A third independent curvature.** `725 = 2.5 * 290` and shares the parent
   prime 29 with 145 and 290, so it extends `p50` from two sizes to three rather
   than starting a fourth lineage. The 55% discrepancy is then tested outside
   the two lineages that produced it.
2. **It breaks the degeneracy Gate 3 exposed.** On the committed sizes, a 5-adic
   valuation of 2 and the *absence* of an interpolating spin-0 combination
   coincide exactly — 325 and 425 have both, so `multiplier`, `valuation` and
   `interpolation flag` are one partition and no screen can tell them apart.
   `725` has valuation 2 and *does* admit an interpolating combination, because
   `(26, 7)` and `(23, 14)` straddle zero in `cos 4 theta`.
3. **A longer lever arm** for the exponent than the current maximum of 425.

`N = 338` was the cheaper candidate for (2) and is **unusable**: `(17, 7)` is
its only primitive representative, so it has no second orientation and no spin-0
combination at all.

## Not established

- The exponent is not identified with any named percolation exponent.
- The second-difference discrepancy is measured, not explained. Two exponents
  would produce it; so would a step-size-dependent reconstruction bias.
- These five transitions are one correlated evidence block with #582 and with
  Gate 3. They are the same histograms read three ways.
- Nothing here revises Gate 3's verdict. No label was retested.
