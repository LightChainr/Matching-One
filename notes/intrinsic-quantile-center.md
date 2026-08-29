# Intrinsic quantile-center spectroscopy of the nonlinear thermal field

Status: C0 definition freeze / C1 N=10 oracle for issue #101. Not a P43 or
Issue #57 target.

## Frozen monomials

The threshold-rank pipeline already solves intrinsic levels. This note freezes
the coordinates themselves as a measurement of the even nonlinear thermal
field, *before* looking at P43.

```text
Mbar_N(p_-^u) = -u,     Mbar_N(p_+^u) = +u
c_u = (p_+^u + p_-^u)/2
w_u = (p_+^u - p_-^u)/2
Q_N = c_{0.05} - c_{0.025}
```

The level set is exactly `u={0.025, 0.05}`. Do not add levels after looking at
outcomes. Under a local inversion of `t(delta)` with `z=t L^{3/4}`,

```text
Q_N ~ C N^{-3/4}.
```

Unknown `p_c` and any common center shift cancel in the midpoint difference.
On a true doubling lineage the leading no-fit prediction is

```text
Q_{2N}/Q_N = 2^{-3/4} ≈ 0.5946035575013605.
```

Separately, `w_u N^{3/8}` is the thermal-metric companion and should connect
to the resolved P49 center-slope correction.

## Exact N=10 oracle

The C4 self-matching control has `M_{10}(p)=2 I_p(3,3)-1`, which is odd about
`p=1/2`. Therefore `c_u=1/2` for every `u` and `Q_10=0` exactly. This is a
solver/oddness check, not a continuum amplitude.

## Descriptive P49 N=130/170

The committed clean P49 archive on `main` contains the doubling *children*
N=130 and N=170, not their parents. Those two sizes are not a doubling pair,
so `Q_170/Q_130` is not a test of `2^{-3/4}`.

Reconstructed full-sample values (rank-2 cross, 100M, 100 batches):

```text
N=130  Q = -4.035e-6    Q N^{3/4} = -1.554e-4    w_0.025 N^{3/8} ≈ 0.01431
N=170  Q = -3.355e-6    Q N^{3/4} = -1.580e-4    w_0.025 N^{3/8} ≈ 0.01432
```

The scaled widths are stable; the scaled midpoint difference is a small
negative number of similar magnitude on the two children. Both statements are
development-only. They are not a P43 or Issue #57 target, and they were not
used to choose `u`.

## Governance

- Retrospective existing-size numbers are development only.
- Any amplitude intended for N=185/265 or 145→290 must be the frozen
  monomials above, committed before those coordinates are read by a scorer.
- A claim-bearing score recomputes `p_±^u` inside each jackknife replicate.

## What this does not establish

- that finite-size corrections reduce to bare-to-scaling-coordinate
  nonlinearity (that is the future test, not this freeze);
- a new operator, Jordan/log identification, or H4 uniqueness claim;
- any P43 or norm-5 score.
