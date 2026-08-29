# What the P234 natural-cutoff scorer measures

Status: normalization/theory audit stacked on PR #246.  This note does not
change the frozen scorer or reinterpret any completed run as final.

Primary source: Federico Camia and Yu Feng, *Correlations of the percolation
energy field and its logarithmic partner*, arXiv:2508.16047v2, especially
equations (1.6), (1.9), and (1.22)--(1.26).  A small but important convention
is that their spin two-point coefficient is **the square root of** `C1`:

```text
<psi(z) psi(w)> = sqrt(C1) |z-w|^(-5/24).              (CF 1.6)
```

## 1. Paper fields and the raw cutoff field

Put `A=sqrt(C1)` and define the unmixed fused field

```text
B_delta = (2 delta)^(-25/24) eta_delta.
```

Camia--Feng define

```text
hat(phi)_delta = B_delta + kappa_paper log(2 delta) phi,
kappa_paper    = C1 CL/C2 = A^2 CL/C2.                (CF 1.23)
```

Hence `B_delta=hat(phi)_delta-kappa_paper log(2delta)phi`.  At two
macroscopically separated points, write `G2` for the common energy two-point
shape.  On the plane, `G2(z,w)=|z-w|^(-5/2)`.  Equations (1.25)--(1.26) give

```text
<phi B_delta> = C2 G2 + o(1),

d/dlog(2delta) <B_delta B_delta>
  = -2 C1 CL G2 + o(1)                                (equal cutoffs).
```

The second line is independent of the additive representative
`hat(phi)->hat(phi)+c phi`, because `<phi phi>=0`.

## 2. The endpoint-connection normalization divides by sqrt(C1), not C1

The scorer does not know `pi_a` separately.  It uses the same-stream endpoint
connection probability

```text
p_conn(a,delta) = <S(z-delta) S(z+delta)>.
```

Equation (1.6) implies

```text
pi_a^(-2) p_conn(a,delta)
  -> A (2delta)^(-5/24),             A=sqrt(C1).
```

After the natural realized-cutoff correction, the scorer's top field is
therefore

```text
T_delta
 = (2delta)^(-5/4) E_delta/p_conn
 -> B_delta/A.
```

Consequently the continuum coefficients in
`scripts/score_p234_cross_cutoff_shear.py` have the exact dictionary

```text
LD_continuum              = (C2/sqrt(C1)) G2,
DD_log_2delta_slope       = -2 CL G2,
kappa_proxy=-s/(2 B)      = sqrt(C1) CL/C2.
```

Here `B` denotes the fitted `LD_continuum`, not the paper field `B_delta`.
For the plane scorer one strips `G2=|z-w|^(-5/2)` from both coefficients.
On the fixed torus, the ratio still cancels the common shape by the cutoff
shear identity; assigning either individual coefficient a plane amplitude
requires a separate torus-to-plane geometry step.

Thus the proxy is neither `CL/C2` in the paper's notation nor the paper's
mixing coefficient `C1 CL/C2`.  It is the paper coefficient divided by
`sqrt(C1)`.  Equivalently, it is the mixing coefficient in the scorer's
connection-normalized top-field gauge.

## 3. Why `8/3` is a sharp conjecture, not an exponent theorem

The current partial value supplied for this audit is

```text
kappa_proxy = 2.653 +/- 0.448                         (partial, not final),
8/3         = 2.666666...
```

The difference is only `-0.0305` standard errors.  This makes `8/3` an
excellent value to freeze as a high-risk target, but does not derive it.

The tempting observation is

```text
8/3 = 1/(3/8),
```

where `3/8` is the thermal RG eigenvalue when area `N=L^2` is the scale
coordinate (`y_t=3/4` in linear size).  Camia--Feng do not prove a Ward or
Russo identity relating this RG eigenvalue to `C1,C2,CL`.  Their proof obtains
`C2` and `CL` from distinct arm-event limits.  A collision mechanism fixes the
presence and transformation law of the logarithm; it does not by itself fix
the relative normalization of the bottom field.

There is also an exact gauge obstruction.  Under `phi -> lambda phi`, with
the connection-normalized top field held fixed,

```text
C2/sqrt(C1) -> lambda C2/sqrt(C1),
kappa_proxy -> kappa_proxy/lambda,
```

while the exponent `3/8` is unchanged.  Therefore no exponent-only argument
can force `kappa_proxy=8/3`.  In the present natural triangular-lattice energy
normalization the conjecture is the concrete amplitude identity

```text
C2 = (3/8) sqrt(C1) CL.                                (H8/3)
```

A future Russo/dilation derivation would have to prove exactly this
normalization-sensitive identity.  Until such a bridge exists, `(H8/3)` is a
new amplitude conjecture, not a consequence of the thermal exponent.

The Vasseur--Jacobsen--Saleur coefficient `2 sqrt(3)/pi` belongs to a
different generic-`Q` observable and field normalization (arXiv:1206.2312,
equations (16)--(17)).  It supplies no conversion to `(H8/3)`.

## 4. Minimal gauge-invariant extra observable

The cheapest additional statistic is one top field with two
two-point-normalized macroscopic spin insertions.  Define

```text
chi = psi/C1^(1/4),       so <chi(z)chi(w)> has unit coefficient,
t3  = d/dlog(2delta) <T_delta(z1) chi(z2) chi(z3)>.
```

Equations (1.9), (1.23), and (1.24) give, after stripping the known plane
three-point shape `F3=F(z1,z2,z3)`,

```text
t3/F3 = -CL,
s2/G2 = -2 CL,
```

where `s2` is the top--top cutoff slope already fitted by P234.  Therefore

```text
CL_invariant = -2 (t3/F3)^2 / (s2/G2) = CL,            (I)
```

and the accompanying one-parameter shear gate is

```text
(s2/G2) / (2 t3/F3) = 1                               (in the p_conn gauge).
```

The squared ratio `(I)` is invariant under an arbitrary multiplicative
rescaling of the top field; it contains no bottom field and is unaffected by
`hat(phi)->hat(phi)+c phi`.  Normalizing `chi` by its measured two-point
function removes the spin-field scale.  It is therefore suitable for a
cross-lattice comparison of the universal candidate `CL`.

Implementation is low cost.  The new correlator is another four-spin
connected function: one pair defines `E_delta`, and the other pair supplies
the two macroscopic spin insertions.  The cluster-sign parity kernel and
same-stream U-statistic already frozen in PR #246 apply without change.

This measurement separates two questions cleanly:

1. Does the gauge-invariant `CL_invariant` agree across triangular site,
   square bond, and square site?
2. In the declared natural local-energy gauge, does
   `kappa_proxy=CL/(C2/sqrt(C1))` additionally equal `8/3`?

Passing (1) but not (2) means that `8/3` was a microscopic energy-gauge value,
not a universal LCFT coupling.  Passing both would isolate the missing
relation as the strong identity `(H8/3)`.
