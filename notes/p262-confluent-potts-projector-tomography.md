# Exact low-leg Potts projector tomography at `Q -> 1`

This note isolates the smallest exact algebra relevant to Issue #262.  It is not a full partition-algebra construction, and it does not identify a lattice observable with a continuum field.

## 1. Connectivity Gram matrices

For set partitions `pi,sigma` of `n` marked points, use the standard connectivity Gram matrix

```text
G(pi,sigma) = Q^{|pi join sigma|}.
```

The zeta/Moebius factorization of a join matrix gives

```text
det G_n = product_{pi} (Q)_{|pi|}
        = product_{k=1}^n (Q)_k^{S(n,k)}.
```

Thus

```text
det G_2 = (Q)_1 (Q)_2,
det G_4 = (Q)_1 (Q)_2^7 (Q)_3^6 (Q)_4.
```

At `Q=1` both raw Gram matrices have rank one.  In the Moebius-exactified basis, every `k`-cluster direction has pivot `(Q)_k`, so a `k`-cluster norm loses order `k-1`.  This is exact rank loss, but by itself says neither which Potts irrep participates nor which scaling dimensions collide.

## 2. Minimal singlet/[2] carrier

The unordered pair of distinct cluster colours is the smallest carrier that can encode a two-cluster/four-leg insertion.  For generic integer `Q>=4`,

```text
C[{a,b}: a!=b] = [] + [1] + [2].
```

Let `I` be the pair identity, `J` the all-ones tensor, and

```text
X_{ab,cd}=delta_ac+delta_ad+delta_bc+delta_bd.
```

The exact projectors are

```text
P_0 = 2J/[Q(Q-1)],
P_1 = (X-4J/Q)/(Q-2),
P_2 = I-P_0-P_1.
```

Their categorical traces are

```text
tr P_0 = 1,
tr P_1 = Q-1,
tr P_2 = Q(Q-3)/2.
```

The script verifies exact orthogonality and idempotence for several integer `Q` realizations.

## 3. The `Q=1` confluence

Put `epsilon=Q-1`.  Coefficientwise in the diagram basis,

```text
epsilon P_0 ->  2J,
epsilon P_2 -> -2J,
P_1 -> -X+4J.
```

The singlet and `[2]` projectors therefore do not separately have finite limits.  They have opposite, rank-one formal residues, while their sum is regular:

```text
P_0+P_2 -> I+X-4J,
d_Q(P_0+P_2)|_{Q=1}=X.
```

The trace collision is equally sharp:

```text
tr(P_0+P_2)=(Q-1)(Q-2)/2 -> 0.
```

By contrast, `P_1` is coefficientwise regular; only its categorical trace vanishes.  Therefore a normalization pole in a `[1]` observable is not forced by this projector algebra.

## 4. What this does and does not prove about a VJS logarithm

The exact opposite residues are the tensor half of a confluent logarithmic pair.  They are not the spectral half.  Suppose singlet and `[2]` contributions have opposite residue `R/(Q-1)` and dimensions `x_0(Q),x_2(Q)`.  Only if

```text
x_0(1)=x_2(1)=x_*
```

does their finite sum contain the conditional logarithm

```text
2 R [x_2'(1)-x_0'(1)] log(r) r^{-2x_*}.
```

Here `R=2J` in the displayed normalization.  A nonzero slope difference is independent input.  In particular, `d_Q P_0` and `d_Q P_2` each have a double pole, so an individual projector derivative is not a finite top field.  The finite collision combination must be formed first.

This gives a low-cost VJS test: recover both generic-`Q` dimension slopes and check that the measured logarithmic tensor lies in the exact `J` residue direction.  A logarithm with a different colour tensor cannot come from this minimal singlet/[2] collision.

## 5. Keep three derivatives separate

For `O_lambda(Q)=P_lambda(Q) O_bare(Q)`, the formal derivative has three logically distinct pieces:

```text
Cov(P O_bare,T_Q)             # Issue #258 measure score
<finite confluent d_Q P term> # Issue #262 projector geometry
<P d_Q O_bare>                # explicit operator definition
```

The first is a change of measure.  The second is singular until the colliding channels are combined.  The third is not determined by either.

Issue #252's `4:-6:3` Ward row is separate again: it distinguishes the thermal `Q4` descendant from a spin-4 primary such as `V_(2,2)`.  It is a Virasoro-module gate, not a `Q`-derivative contribution.  The efficient order is therefore:

1. use the Ward row to decide primary versus thermal descendant;
2. use the Potts colour tensor to decide singlet versus `[2]` confluence;
3. only then decompose a `Q` tangent into measure, projector, and bare-operator pieces.

## Frozen risky predictions

1. Any minimal singlet/[2] confluent residue has exact tensor direction `J` and residue ratio `-1`.
2. `[1]` is regular-null, not pole-colliding, at this level.
3. A VJS logarithm from this collision must factor as the exact tensor residue times a nonzero, independently measured dimension-slope difference.

All three are falsifiable without fitting a large partition-algebra model.
