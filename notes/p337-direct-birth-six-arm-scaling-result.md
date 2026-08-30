# P337 result: the no-correction 5/6 line fails while the doubling flow approaches it

## Four-generation collision law

The direct `0 -> 2` essential-birth mass decreases cleanly across the exact
P337 Gaussian-child lineage:

| N | D_first | D_second | Dbar +/- SE | Dbar N^(5/6) |
|---:|---:|---:|---:|---:|
| 85 | `0.0156250` | `0.0152600` | `0.0154425 +/- 0.0002414` | `0.625989` |
| 170 | `0.00850050` | `0.00847463` | `0.00848756 +/- 0.0000232` | `0.613042` |
| 340 | `0.00471867` | `0.00470692` | `0.00471279 +/- 0.0000135` | `0.606518` |
| 680 | `0.00264499` | `0.00264391` | `0.00264445 +/- 0.00000350` | `0.606400` |

The two orientations remain one paired block per size.  Their measured
correlations are retained in the score and range from `-0.159` to `0.206`;
they are not counted as eight independent observations.

## Frozen 5/6 decision

The single-amplitude, no-correction model

```text
Dbar_N = A N^(-5/6)
```

gives

```text
A = 0.607576 +/- 0.000666,
chi2 = 16.886 / 3,
p = 7.46e-4.
```

It is rejected at the frozen `alpha=0.01`.  The retrospective N680 heldout
view reaches the same strict decision narrowly: the N85/N170/N340 amplitude
predicts `0.00266090`, versus `0.00264445` observed (`z=-2.621`, `p=0.00878`).

The free effective power is

```text
beta_eff = 0.840695 +/- 0.002041,
chi2 = 3.879 / 2,
p = 0.144.
```

Its difference from `5/6=0.833333...` is resolved (`p=3.10e-4`).  The frozen
one-coordinate log-curvature diagnostic gives `kappa=0.00522 +/- 0.00265`,
but is not resolved at the project threshold (`p=0.0489 > 0.01`).  Thus the
present data eliminate the exact no-correction line; they do not identify a
specific correction law.

## The scientifically interesting part: flow toward 5/6

The three adjacent doubling ratios are

```text
N85 -> N170:   0.549624
N170 -> N340:  0.555259
N340 -> N680:  0.561122
fixed 5/6:     0.561231
```

They approach the fixed ratio monotonically.  The final ratio differs by only
`1.95e-4` in log units (`z=-0.062`, `p=0.951`).  Likewise the scaled amplitude
`Dbar N^(5/6)` changes by less than `0.02%` between N340 and N680 after falling
by about `3.1%` over the first two generations.

This produces a sharper judgment than a binary model vote:

> A correction-free four-size 5/6 law is false at current precision, but the
> largest-generation flow lands almost exactly on its parameter-free doubling
> ratio.  The data are a strong conditional clue for a correction-bearing
> six-arm scaling line, not a proof of a six-arm event correspondence.

The free exponent `0.8407` should therefore be read as a four-generation
effective exponent, not automatically as a new asymptotic exponent.

## External N325/N425 lineage

The independent P334 production lineage has

```text
N325 Dbar = 0.00491450 +/- 0.0000359
N425 Dbar = 0.00385275 +/- 0.0000315.
```

It remains only an external geometry control.  Its previously frozen two-point
ratio was compatible with `5/6`, but it is not mixed into this fit or counted
as another generation because its shapes and genealogy differ.

## Provenance and boundary

No new Monte Carlo was run.  The scorer verified all four committed raw hashes,
including the compressed and uncompressed N680 hashes, and checked that every
sparse table partitions every batch/orientation sample count exactly.  The
committed batch-sufficient-statistics table contains 260 paired rows and
reconstructs all orientation and lineage covariances.

This is not an exact arm correspondence, universal amplitude, universality
test, or statement about thermal Q4.  A direct-birth event can have the same
scaling exponent without being literally identified with a continuum six-arm
event.  Conversely, failure of the correction-free form does not eliminate an
asymptotic 5/6 mechanism.

## Scientific card

- Mechanism space changed: the strict single-amplitude 5/6 line is removed,
  while the last doubling step lands on its parameter-free ratio.
- Not proved: arm correspondence, asymptotic exponent, universal amplitude,
  or cross-geometry universality.
- Observer/sector/source/geometry: direct `K1=K2` birth mass; paired P337
  first/second Gaussian orientations; exact N85/170/340/680 child lineage.
- Dependency groups: one paired-orientation block per size; four independent
  seed/counter domains.
- Next lift: another same-lineage generation or an independently frozen
  correction amplitude, not refitting a flexible tail law to these four points.
