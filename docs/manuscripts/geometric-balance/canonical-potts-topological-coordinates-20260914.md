# Canonical Potts topological coordinates: loop fugacity, thermal self-duality and rank charge

Date: 2026-09-14

Status: exact algebra for a toroidal FK model with the bounded rank source.  It reorganizes the generic-Q source/tangent programme and makes the graph-polynomial source a symmetry point rather than an ad hoc closure deformation.

## 1. Start from the sourced FK weight

For a cellular torus FK subgraph `A`,

```text
W(Q,v,h;A)
 = Q^{k(A)} v^{|A|} exp[h X(A)],
X(A)=r(A)-1.
```

The Euler identity gives

```text
X(A)=k(A)-k(A*)+|A|-|V|.
```

## 2. Define three canonical coordinates

Introduce

```text
boxed:
x   = v/sqrt(Q),
eta = h + (1/2) log Q.
```

Equivalently,

```text
v=x sqrt(Q),
h=eta-(1/2)log Q.
```

Substitute into the original weight and use the Euler identity.  All powers of `Q` simplify:

```text
W(Q,x,eta;A)
 = Q^{|V|/2}
   (sqrt Q)^{k(A)+k(A*)}
   x^{|A|}
   exp[eta X(A)].                                   (2.1)
```

The prefactor `Q^{|V|/2}` is global.

Using

```text
k+k*=b+1_(r!=1),
```

we obtain the loop/topology form

```text
boxed:
W ∝
 (sqrt Q)^{b(A)+1_(r(A)!=1)}
 x^{|A|}
 e^{eta(r(A)-1)}.                                    (2.2)
```

This is an exact finite configuration identity.

## 3. Interpretation of the three coordinates

Equation (2.2) cleanly separates three roles.

### `Q`: common loop fugacity

Every medial boundary loop carries the usual factor `sqrt(Q)`, together with the fixed torus extreme-sector bonus already derived.

### `x=v/sqrt(Q)`: thermal / self-duality coordinate

The self-dual bond-FK point is simply

```text
x=1.
```

The logarithm `t=log x` is the natural duality-odd thermal coordinate at fixed `Q`.

### `eta`: primal/dual topological charge

The only remaining rank source is

```text
e^{eta(r-1)}.
```

Thus `eta`, not the raw closure source `h`, is the canonical generic-Q topological chemical potential after the intrinsic FK cluster asymmetry has been removed.

## 4. The graph-polynomial source is exactly `eta=0`

The graph-polynomial balance source is

```text
h_Q=-1/2 log Q.
```

Therefore

```text
boxed:
eta(h_Q,Q)=0.
```

At this source the explicit topological charge factor disappears from (2.2).  The only topology dependence left is the universal torus sector bonus already carried by the medial-loop count.

This is why the primal and dual cluster fugacities become equal in the Euler-localized representation.

Hence the graph-polynomial section is not best thought of as “turn on a special topological field”.  In canonical coordinates it is the **zero topological-field section**.

## 5. The square self-dual graph-polynomial point is `(x,eta)=(1,0)`

At the square-lattice self-dual Potts/FK coupling

```text
v=sqrt(Q),
```

we also have

```text
x=1.
```

Therefore the distinguished point is

```text
boxed:
(log x, eta)=(0,0).
```

It is simultaneously

```text
thermal self-dual,
primal/dual topological-source symmetric.
```

This gives a compact structural explanation for the unusually strong cancellation exploited by the graph-polynomial/eigenvalue method.

## 6. Exact duality action

Let `A*` be the complementary dual state.  Then

```text
|A*|=|E|-|A|,
X(A*)=-X(A),
b(A*)=b(A),
1_(r(A*)!=1)=1_(r(A)!=1).
```

Therefore the weight (2.2) obeys, up to the global factor `x^{|E|}` and replacement by the dual lattice,

```text
boxed:
(log x, eta) -> (-log x, -eta).                       (6.1)
```

On a self-dual lattice the partition function has the corresponding two-coordinate reflection relation.

Thus thermal duality and topological charge conjugation are simply the two odd coordinates around the same fixed point.

## 7. Reinterpretation of the generic-Q shifted rank coordinate

The existing canonical rank coordinate is

```text
b = (1/2) log(Z0/Z2).
```

The graph-polynomial balance occurs at

```text
b=-1/2 log Q.
```

Hence the shifted odd variable

```text
b_tilde=b+(1/2)log Q
```

is precisely the rank-law counterpart of the canonical source coordinate `eta`.

Both vanish on the graph-polynomial symmetric section.

This makes the earlier algebraic shift much less mysterious: it removes the built-in primal FK cluster fugacity and centers the theory on equal primal/dual cluster weight.

## 8. Reinterpretation of the graded Q tangent

Differentiate at fixed `eta=0`.

Since

```text
eta=h+(1/2)log Q,
```

keeping `eta=0` forces

```text
dh/d log Q = -1/2.
```

This is exactly the previously identified kinematic topological-source term in the generic-Q tangent.

Therefore the `-1/2` piece has a geometric meaning:

> it is the compensating source motion required to stay on the primal--dual symmetric `eta=0` surface while Q changes.

After this subtraction, the remaining Q tangent probes the common loop fugacity `sqrt(Q)`, any Q-dependence of the thermal coordinate `x`, and representation/module data.  It is not contaminated by a trivial topological imbalance.

This is the natural coordinate system for #746.

## 9. Near-critical/source scaling should use `(t,eta)`

At fixed Q define

```text
t=log x.
```

The exact duality action is

```text
(t,eta)->(-t,-eta).
```

Therefore a near-critical topology-resolved scaling function should be organized around these two odd coordinates, not raw `(v,h)`.

At Q=1,

```text
eta=h,
x=v,
```

so the distinction disappears and the Matching-One charge source is already canonical.

At generic Q it is essential.

## 10. Possible loop/height interpretation

Equation (2.2) suggests a continuum dictionary:

```text
sqrt(Q) : loop fugacity / Coulomb-gas background parameter,
t         : thermal mass perturbation,
eta       : electric/topological face-charge perturbation.
```

At `eta=0` the primal/dual faces are unbiased.  Moving `eta` changes their relative cluster fugacities while leaving the product fixed.

An oriented-loop/height representation should therefore encode `eta` as an electric/background-charge-like deformation rather than as a change of the unoriented bulk loop fugacity.

This is a targeted representation conjecture, not yet a normalized Coulomb-gas charge assignment.

## 11. Consequence for the massive torus programme

The object sought in #782 can be restated as a finite-volume scaling function in

```text
(Q,t,eta,alpha,tau),
```

where `alpha` is the already-identified rank-one neutral seam.

The graph-polynomial / critical UV anchor is the slice

```text
t=0,
eta=0.
```

The Matching-One rank response is the `eta` derivative (at Q=1), while the thermal response is the `t` derivative.

This provides a cleaner target for a massive modified trace than separate ad hoc projectors for rank0, rank1 and rank2.

## 12. Claim boundary

Exact finite algebra:

```text
x=v/sqrtQ,
eta=h+1/2 logQ,
W ∝ (sqrtQ)^{b+1_(r!=1)} x^{|A|} e^{eta(r-1)},
duality: (log x,eta)->(-log x,-eta),
graph-polynomial source: eta=0.
```

Conjectural/interface:

- identification of `eta` with a particular electric/height operator in the continuum;
- local periodic-TL implementation for arbitrary eta;
- massive/TBA scaling functions in these coordinates.

The main structural payoff is that generic-Q thermal and topological asymmetries are now centered on one explicit duality fixed point.