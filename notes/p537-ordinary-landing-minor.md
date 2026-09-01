# P537: an ordinary landing already has rank two

## Result

The proposed rank-one ordinary-landing block fails on the smallest honest
square quotient on which the relevant ports are unaliased.  On the `3 x 3`
torus, take the two occupied sites `(0,0),(1,0)` and flip `(2,0)`.  Before the
flip the occupied component is a branch-free two-site path; after the flip it
is the primitive horizontal three-cycle.  The four incident ports alternate
vacant/occupied/vacant/occupied, the ambient rank changes `0 -> 1`, and no
extra occupied branch is present.

For the canonical vacant pair used by the exact port kernel, the same flip
changes

```text
g_xy : 0 -> 1/4,
q    : -1 -> 0,
E    :  1 -> 0.
```

Rotate this complete marked state by `0, pi/2, pi, 3pi/2` and sum all four
copies.  Spin four has character `+1` under this `C4` action, so this is the
required angular sum, not a post-hoc choice of one orientation.

For `O in {q,E}`, split the exact covariance derivative on this landing orbit
into its two mandatory pieces,

```text
K_O = sum_orbit w_p(omega) (O_mid-E_p O) D_z g_xy,
P_O = sum_orbit w_p(omega) (g_mid-E_p g_xy) D_z O.
```

Here `g=g16/16` is the canonical completed pair kernel and
`w_p=p^|omega|(1-p)^(N-1-|omega|)`.  The finite landing transfer matrix is

```text
             kernel reconnection   readout pivotal
q                    K_q                  P_q
E                    K_E                  P_E .
```

At the rational control `p=1/2` it is exactly

```text
[ -11/16384    11/8192 ]
[ -47/65536   -11/8192 ]
```

and its determinant is

```text
1001 / 536870912 != 0.
```

This value is only a transparent control.  The decisive statement is at the
actual finite matching root.

## Exact finite-root certificate

For this torus the matching mean is

```text
M(p)=E_p q=-4p^9+18p^8-18p^7+6p^3-1.
```

Monotonicity of ambient rank makes its root in `(0,1)` unique; exact rational
evaluation brackets it by

```text
1173/2000 < p_* < 2933/5000,
```

with `p_*=0.5865114551126757...`.  The canonical determinant polynomial
factors as

```text
p^9 (1-p)^12
 (p^4+9p^3-18p^2+9)
 (48p^9-224p^8+404p^7-332p^6+96p^5+16p^4
  -4p^3-4p^2+1).
```

Its polynomial gcd with `M(p)` over `Q[p]` is exactly `1`.  Hence it cannot
vanish at `p_*`; the displayed decimal value `7.4089043146e-7` is merely a
locator, while nonvanishing is algebraic.

Root conditioning does not rescue rank one.  At `p_*`, put

```text
R=(partial_p E_p E)/(partial_p E_p q)
```

and perform the Schur row operation `E -> E-Rq`.  This left multiplication has
determinant one, so the same `2 x 2` minor is unchanged.  Numerically
`R=0.021325579095...`, but no numerical approximation to `R` enters the exact
nonvanishing certificate.

## Consequence and boundary

This is a physical finite counterexample to the statement that every
branch-free ordinary four-arm landing transfer is thermal-only after `C4` and
root Schur projection.  The four-packet route in the P537 audit therefore
cannot proceed by claiming that nonzero projected terms require an extra
branch: an ordinary first-rank landing already supplies two independent
channels.  The next meaningful object is the surviving signed landing
functional, including both reconnection and readout terms.

The result does **not** prove that the full spatial sum is nonsummable, give an
arm exponent, identify a continuum field, or establish a nonzero asymptotic
original-`U` amplitude.  It only closes the proposed rank-one cancellation
gate at its advertised bounded falsifier.

## Reproduction

```bash
python experiments/p537-landing-minors-20260901/landing_minor.py \
  --output /tmp/p537-landing-minor.json
python -m unittest discover \
  -s experiments/p537-landing-minors-20260901 \
  -p 'test_landing_minor.py' -v
```

The script uses the already integrated physical torus and canonical kernel
implementation from `p337-thermal-gate-audit-20260901`; it enumerates only
`2^9` states and uses standard-library exact rational polynomial arithmetic.
The checked-in machine result is
`results/p537-landing-minors/exact-n9.json`.
