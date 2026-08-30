# P250 near-zero reveal: support before phase

Status: existing archives reanalyzed, followed by a capped 4,000-replica local
multi-separation smoke.  No production is proposed or run here.

## The quantity

Let `O_r(x)` be the existing local landing-pivotal H4 row after the transported
Z5 fiber DFT.  For the translated right triangle

```text
x0=(0,0), x1=(d,0), x2=(0,d),  d in {1,2,3},
```

measure the neutral Hermitian two-point function

```text
G_r(d) = 1/4 sum over x/y axes and reversed charges
         O_r(x0) O_-r(x0+d axis).
```

It is exactly real configuration by configuration.  The separation-normalized
cubics are

```text
Omega_113(d) = C_113(d) / sqrt(|G_1(d)|^2 |G_2(d)|),
Omega_122(d) = C_122(d) / sqrt(|G_1(d)| |G_2(d)|^2).
```

The positive denominator removes charged-field magnitudes while deliberately
leaving the phase in the repository's exact transported deck basis.  A second
normalization replaces `G_r` by the local positive variance `V_r=E|O_r|^2`;
this remains stable when a separated two-point function crosses zero.

There is no noisy disconnected subtraction.  Charges `113` and `122` have no
neutral proper subset, so their raw third moments are their connected third
cumulants by exact Z5 charge conservation.

## Why the old files cannot supply it

The one-million-row result at `be80f25` retains batch sums only after the three
local charged rows have already been multiplied.  It has no individual
`O_r(x)`, pair product or separation label.  The P226 archive is a different
global marked-row one-point response with a different root schedule.  It
cannot normalize this local cubic, and the two archives do not contain the
cubic/pair cross-covariance needed for a ratio.

A support-first reanalysis of the existing one-million result gives

```text
raw eight-real cubic support: chi2=8.564036420118496/8, p=0.3803958155142504
phase closure:                chi2=0.8396327481811326/2, p=0.6571674817128721
```

The first line does not detect a nonzero cubic vector.  Therefore the second
line is reclassified as `not_interpretable_until_nonzero_support`, rather than
as evidence that a common phase survived.

## New score order

For each separation:

1. require all four `(G1,G2) x (plus,minus)` denominators to have `|z|>=2`;
2. test the eight-real local-variance-normalized cubic vector against zero;
3. evaluate separation-normalized phase closure only if both denominator and
   cubic-support gates pass.

The closure diagnostic is still saved when a gate fails, but it receives no
scientific phase interpretation.  This prevents a near-zero vector from
making every phase model look successful.

## 4k plumbing/variance smoke

The smoke used seed `25011312220260901`, counters `[0,4000)`, 40 batches,
eight local workers, `p=0.59274605079` and radius one.  The runner has a hard
5,000-replica cap.

| d | weakest two-point denominator | normalized cubic support | phase status |
|---:|---:|---:|---|
| 1 | 4.51 sigma | `11.479/8`, p `0.176` | not interpretable |
| 2 | 0.48 sigma | `9.766/8`, p `0.282` | not interpretable |
| 3 | 0.0039 sigma | `8.542/8`, p `0.382` | not interpretable |

All 5,850 exact parent/root labels pass; maximum Hermitian-pair imaginary
residual is `1.39e-17` and maximum DFT conjugacy residual is `1.25e-17`.

This is a sharper negative selector than simply adding cubic samples.  At
`d=1` the charged two-point denominators are already resolved, but the cubic
support is not.  At `d=2,3`, the separated two-point denominators themselves
are mostly unresolved and the formal ratios become unstable.  The compact
local-H4 charged cubic is therefore contact-range/signal-limited in two
distinct ways.

## Consequence

Do not add another large sample block to the existing `d=1` cubic.  If P250 is
continued, change the insertion before the budget: use a charged operator with
a resolved two-point tail (or a larger-radius/leg-defect row), then rerun this
same support-first score at two separations.  Phase closure becomes meaningful
only after a nonzero cubic vector is established.

This smoke does not reject a continuum charged OPE.  It says the current
radius-one landing-H4 lattice insertion does not expose it with a stable
separation-normalized signal.
