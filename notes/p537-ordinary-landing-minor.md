# P537 supporting control: a single-geometry kernel/readout minor

## Correct status

This N9 calculation is an exact **supporting control**, not by itself the
advertised source/thermal `P4` landing certificate.  It proves that the two
terms in the covariance-derivative identity can form a rank-two block inside
one C4-invariant square geometry.  It does not construct the original
axis-versus-tilted `P4` projector, and its two columns are kernel reconnection
and readout pivotal rather than two independent source/thermal deformations.

The separate N16 three-fibre calculation in
`notes/p537-ordinary-four-arm-landing-minor.md` is stronger geometrically but
has the same final boundary.  Its raw jump coordinates
`k=D_z g_xy`, `h_R=chi D_zE-RD_zq` give the formal determinant `-chi/2`, yet
common off-site occupation does not cancel the fibre-dependent bilinear
`H_i B_i`, and `chi` is not produced by a paired axis/tilted normalization.
Thus the full projected minor remains open.

The N9 control remains useful because it is at an actual finite matching root
and retains both exact midpoint terms.  On the `3 x 3` torus, take the two
occupied sites `(0,0),(1,0)` and flip `(2,0)`.  Before the flip the occupied
component is a branch-free two-site path; after the flip it is the primitive
horizontal three-cycle.  The four incident ports alternate
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
copies.  Spin four has character `+1` under this `C4` action, so the sum is a
valid C4-even local block.  This is **not** the original geometric projector
`P4[X]=(X_axis-X_tilted)/DeltaCos4`: quarter-turn-related copies have the same
`cos(4 theta)`, and no tilted second quotient occurs in this calculation.

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

The simple root row shear does not remove this raw minor.  At `p_*`, put

```text
R=(partial_p E_p E)/(partial_p E_p q)
```

and perform the Schur row operation `E -> E-Rq`.  This left multiplication has
determinant one, so the same `2 x 2` minor is unchanged.  Numerically
`R=0.021325579095...`, but no numerical approximation to `R` enters the exact
nonvanishing certificate.  This row operation is only the `E-Rq` part; it is
not the complete bilinear Schur contribution `-beta H_iB_i`.

## Semantic boundary and consequence

The matrix above has rows `q,E` and columns equal to the **kernel** and
**readout** pieces of one derivative `d_p Cov(O,g_xy)`.  Splitting an identity
into these two terms does not manufacture independent physical sources, so
its nonzero determinant must not be called the Issue's source/thermal minor.
Likewise, invariance under the four quarter turns only supplies the correct C4
representation check; it does not establish a nonzero axis-versus-tilted
contrast.

What the N9 result does prove is narrower and still useful: the exact midpoint
decomposition is not internally rank one, and the row shear `E -> E-Rq` does
not make that raw two-term matrix rank one.  The N16 calculation separately
proves that three physical ordinary landing fibres realize two raw jump
symbols.  Neither calculation includes the complete bilinear Schur term and
separately normalized paired-geometry P4 amplitude, so neither yet retires the
thermal-only landing lemma.

Neither result proves that the full spatial sum is nonsummable, gives an arm
exponent, identifies a continuum field, or establishes a nonzero asymptotic
original-`U` amplitude.  The minimum missing certificate must retain, for two
paired geometries and at least two landing boundary states, the complete
midpoint tensor `(q0,q1,E0,E1,a0,a1,K_-z)`, geometry-specific centering,
`R,beta`, separate normalizations and the final axis-minus-tilted P4.  Before
that object is nonzero, no ordinary-block stop claim is justified.

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
