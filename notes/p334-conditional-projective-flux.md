# Conditional projective flux: a line-sorting discriminator

Status: exploratory reanalysis of the existing Issue #334 N65 20k smoke.  No
new random counters and no server were used.  The result is one new coordinate
of the same raw dependency block, not an independent evidence row.

## The missing quotient

The first source/sink analysis retained the complex currents

```text
J4_birth, J4_exit.
```

Those currents still mix two effects:

1. how much line-bearing topological traffic occurs at each boundary;
2. which projective lines carry that traffic.

Divide by the corresponding trivial-character line fluxes:

```text
mu_birth = J4_birth / J0_birth,line,
mu_exit  = J4_exit  / J0_exit,line,

Delta_mu4 = mu_birth - mu_exit.
```

`DIRECT_RANK2` has no projective line and is outside both denominators.  It is
still reported as excluded mass.

This quotient asks whether the projective-line composition entering rank one
is the same as the composition leaving rank one.  If `ell` factorizes from
`(tau1,tau2)`, as it did in every tiny Phase A control, then `Delta_mu4=0`
exactly even when total ingress and egress rates differ.  A nonzero value is a
line-dependent lifetime/exit-sorting effect rather than a scalar traffic
effect.

## Same-modulus baseline cancellation

A single orientation can contain a universal square-modulus line-lifetime
baseline, so its nonzero `Delta_mu4` alone is not a lattice H4 diagnostic.  The
N65 pair solves this without new data: both geometries have modulus `i` but
different microscopic orientation.

For Gaussian representation `(a,b)`, rotate the complex character by

```text
conjugate((a+i b)^4 / N^2).
```

The real component is the declared D4 amplitude in the torus frame.  A
continuum baseline common to the square modulus is the same for both shapes;
their aligned difference therefore isolates microscopic orientation-sensitive
line sorting.

## N65 result

At `p_ref=0.592746050790`, delete-one-batch jackknife over the 20 shared
batches gives

| shape | aligned parallel Delta_mu4 | SE | perpendicular | SE |
|---|---:|---:|---:|---:|
| `(8,1)` | -0.01000528 | 0.00281 | -0.0000134 | 0.0000463 |
| `(7,4)` | +0.00171005 | 0.00265 | -0.0000621 | 0.0000415 |

The two-shape aligned vector has quadratic `12.77 / 2 df`.  More directly, the
same-modulus contrast is

```text
Delta_mu4_parallel(7,4) - Delta_mu4_parallel(8,1)
  = 0.01171533 +/- 0.00435093,

quadratic = 7.25 / 1 df.
```

This is a modest but concrete new discriminator: the first N65 orientation
shows projective line-dependent sorting while its same-modulus partner is
compatible with no sorting at this sample size.  It is not explained by a
common continuum primitive-sector baseline or by changing only total
birth/exit traffic.  The natural next use is a fresh larger counter block or a
child size with the same predeclared quotient; no new angular model is needed.

## Connection to the generic-Q lift

This coordinate lies entirely in the fixed-`Q=1` thermal direction.  If two
generic-Q lifts satisfy

```text
O2(Q,p)-O1(Q,p)=(Q-1)X(Q,p),
```

then every fixed-Q `p` derivative agrees at `Q=1`.  The source/sink process and
`Delta_mu4` are therefore intrinsic horizontal observers under the #333 lift
class.  This is precisely where analysis can proceed without choosing
`L_hom` versus `L_CP`.

A future Q-normal derivative of `Delta_mu4`, or a mixed Q/p collision score,
is different: it acquires the transition counterterm `partial_p X` and must be
transported using the #333 descriptor/connection before pooling.

## Claim boundary

- The `7.25 / 1 df` value is exploratory because the quotient was designed
  after inspecting the same 20k smoke block.
- It is one contraction of the #334 covariance block, not additive evidence.
- No exponent or continuum field is fitted.
- The result supports an orientation-sensitive projective line-lifetime
  coordinate; it does not yet distinguish thermal-Q4 from a charged four-leg
  radial family.

Reproduce with:

```bash
python3 scripts/score_conditional_projective_flux.py \
  --births results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.births.csv \
  --metadata results/local-20260830/P334-projective-birth-N65-smoke/n65_20k.metadata.json \
  --json results/local-20260830/P334-projective-birth-N65-smoke/conditional_flux.json \
  --markdown results/local-20260830/P334-projective-birth-N65-smoke/conditional_flux.md

python3 -m unittest discover -s tests \
  -p 'test_conditional_projective_flux.py'
```

