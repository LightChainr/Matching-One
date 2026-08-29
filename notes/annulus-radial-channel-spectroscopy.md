# Annulus-channel spectroscopy from the PR #247 covariance block

## Outcome first

The existing radii are enough to reject or retain a one-mode model, and weakly
test a repeated-root/Jordan model. They are **not** enough to choose between two
real modes and a complex pair in one channel. The reason is structural rather
than statistical: two geometry readouts at the three dyadic radii `2,4,8`
give six values, exactly the six parameters of an unrestricted rank-two
realization.

The point estimates are nevertheless unusually suggestive:

```text
plus : T= 8.26148, D=4.81135, Delta=+49.0066
minus: T=-0.62044, D=1.96608, Delta= -7.47939
```

Thus the plus point recurrence has two positive real roots, while the minus
point recurrence is a complex-conjugate pair. But delta-method uncertainties
are about `215.06` and `11.48` for the two discriminants. Neither sign is yet a
resolved class measurement. The correct present claim is:

> the minus block has a rotating two-state **point realization**, not detected
> complex scaling dimensions.

This is the annular analogue of the exact N10 finding in #256: a scalar failure
can arise because a two-dimensional covariance plane rotates. The two results
must not be counted twice; #256 supplied the coordinate mechanism, while this
post-reveal score uses the same PR #247 block to design the next acquisition.

## Basis-independent recurrence

Write `g_n` for the two-component vector of N325 and N425 orientation
contrasts at `R=2^(n+1)`. Any rank-two one-step realization obeys

```text
g_(n+2) - T g_(n+1) + D g_n = 0.
```

The characteristic discriminant has the usual meaning:

```text
Delta > 0  two real roots (R2)
Delta = 0  repeated root; J2 versus scalar degeneration needs a rank test
Delta < 0  conjugate complex roots (C2)
```

With only `n=0,1,2`, the two vector components solve exactly for `T,D`. This is
why a tiny residual is not evidence for R2 or C2. The scorer carries the full
PR #247 covariance through the nonlinear recurrence and through a fractional
step prediction at R7.

## What is identifiable now

For each plus/minus channel:

| model | parameters from six source values | source df | present status |
|---|---:|---:|---|
| R1 | 2 amplitudes + one exponent | 3 | testable |
| J2 | 4 amplitudes + one exponent | 1 | weakly testable |
| R2-gap1 | 4 amplitudes + one common shift | 1 | weakly testable |
| R2 | 4 amplitudes + two exponents | 0 | saturated |
| C2 | 4 amplitudes + decay/phase | 0 | saturated and phase-aliased |

`R2-gap1` is the mandatory ordinary-real adversary supplied by PR #260, not a
new derivation here: the exact Q=1 spin-4 fields `V_(2,2)` and thermal Q4 have
dimensions `17/4` and `21/4`, hence a unit relative gap and Gaussian-area
relative transfer `Q^-1/2`. The common shift appropriate to the normalized
marked-pivotal shell is unresolved, so the scorer fixes the gap and profiles
only that shared shift. Generic R2 remains as a saturation boundary, not a
preferred physical model.

The joint plus/minus common-spectrum fit has positive degrees of freedom, but
adds a new assumption: matching-even and matching-odd readouts share one radial
generator. It is reported as a mechanism null, not silently imposed.

The covariance-aware profiles sharpen the original shell result without
overriding it:

```text
plus  R1: chi2=0.670/3, p=0.880
minus R1: chi2=11.624/3, p=0.00879, optimum hits the allowed growth boundary

plus  J2: chi2=0.341/1, p=0.559
minus J2: chi2=0.765/1, p=0.382
minus R2-gap1: chi2=0.765/1, p=0.382
```

So the plus channel closes at rank one, while minus rejects a positive-real
one-mode propagation over this declared profile window. J2 and the exact
ordinary-real gap-one adversary are numerically indistinguishable here; that is
the main acquisition result, not a tie to be broken by rhetoric. A joint common
spectrum also fails to distinguish J2, R2-gap1, generic R2 and C2, all near
`chi2=1.7` with 2--3 residual degrees of freedom and some optima running toward
profile boundaries.

R7 is valuable because `log2(7/2)` is not an integer. A positive-root R2
continuation for plus and the principal-phase C2 continuation for minus give
correlated residual scores. These are post-reveal design checks inside the same
raw block. Complex phases also have `2*pi*k` aliases; the principal branch is a
frozen convention, not a discovered frequency.

## Minimal next acquisition

Acquire one third primitive orientation pair,

```text
N365: (14,13) versus (19,2),
radii: 2,4,7,8.
```

The third readout makes the dyadic `T,D` recurrence overidentified, while R7
remains an off-grid propagation target. This is more informative than adding a
dense radius ladder to the same two readouts.

The follow-up production freeze now lives in
`experiments/p253_n365_annulus_recurrence_20260829.yaml`. It keeps the existing
200k/200-batch resolution, uses a fresh RNG domain, and costs only about half
of PR #247: an empirical linear estimate of 1.96 seconds wall time on the same
16-vCPU Huawei class, with a 2--4 second planning envelope. The freeze did not
start the job.

## Reproduction

```bash
python3 scripts/analyze_annulus_channel_recurrence.py \
  --input results/server-20260829/P225-norm5-multiradius/analysis.json \
  --output results/annulus-radial-design/latest.json
```

The script contains exact synthetic J2, R2 and C2 oracles. Each oracle recovers
the correct discriminant class and predicts the fractional R7 step to floating
roundoff.
