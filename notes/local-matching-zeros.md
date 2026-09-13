# Local matching-polynomial zeros near the physical critical root

Status: C1 descriptive for issue #113. Complex-zero scaling route **closed**
at available exact sizes.

## Why this slice exists

PR #78/#84 falsified a global imaginary-RMS extrapolation of the exact
matching-polynomial root cloud. That failure does not logically exclude a
*local* critical-zero trend. The committed L≤5 coefficients already exist, so
the check is a zero-extra-enumeration reanalysis.

These are zeros of the finite matching polynomial. They are not automatically
Fisher zeros of a partition function. This note does not use Lee–Yang or CFT
language.

## Frozen metrics

Definitions live in `predictions/local_matching_zero_metrics_20260829.yaml`
and must not be changed after looking at a future exact L=6 result.

- `physical_root`: unique real root in `(0,1)`
- `nearest_nonreal`: nonreal root closest to the physical root (upper-half
  tie-break)
- `matching_partner_of_physical`: `1-p*`, the dual-polynomial partner
- named diagnostics only: `L^{3/4} |Im|` and `L^4 |z-p*|`

## What the committed catalogue shows

Axis L=1,2 and diamond L=1 are entirely real. For the four polynomials that
have nonreal roots, the nearest nonreal root usually has real part outside
`(0,1)` (diamond L=3 is the exception). The named diagnostics are

```text
L^{3/4} |Im| :  0.429 (axis 3), 0.730 (axis 4), 0.685 (diamond 2), 0.704 (diamond 3)
L^4 |z-p*|  :  48.1,           131.7,          12.6,              41.0
```

Neither column is stable. L≤5 trends are descriptive only; there is no
held-out exact size on `main`.

## Closure

Issue #113 said: if no stable local behavior appears, close the complex-zero
route rather than inventing more cloud statistics. That gate is met. A later
exact polynomial may reopen the route only by the frozen metrics above.

## What this does not establish

- a continuum zero-density law;
- any identification of matching-polynomial zeros with Fisher zeros;
- a replacement for the already-falsified global RMS prediction.
