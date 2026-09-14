# Krushkal topological variables give the finite torus rank-source projector

Date: 2026-09-14

Status: exact finite embedded-graph topology for bond/FK configurations, plus a new interface proposal for #782/#768.  This does **not** by itself construct the missing massive affine-TL/TBA modified trace, and it is not yet a direct square-SITE transfer identity.

## 1. The current blocker can be separated into definition versus realization

The current #782 dictionary already has a continuous rank-one neutral-cluster fugacity from an affine-TL seam:

```text
z = alpha^2/Q.
```

What remained was described as a missing way to distinguish rank-0 and rank-2 configurations inside the zero-noncontractible-loop sector.

At the level of a finite graph embedded on a torus, however, the required topological coordinate is already standard in the Krushkal / topological Tutte polynomial.

The remaining problem is therefore narrower:

> realize the already-defined rank source as a transfer/modified trace compatible with the massive Potts / affine-TL continuum engine.

## 2. Krushkal variables on a torus

Let a spanning subgraph `H` of a graph cellularly embedded in `T^2` have image

```text
V(H) = im[ H_1(H;R) -> H_1(T^2;R) ].
```

Let

```text
r(H)=dim V(H) in {0,1,2}.
```

Krushkal's surface polynomial uses two topological exponents `s(H), s_perp(H)`.  In the symplectic formulation,

```text
s(H)      = dim[ V/(V cap V^perp) ],
s_perp(H) = dim[ V^perp/(V cap V^perp) ],
```

where orthogonality is for the intersection form on `H_1(T^2;R)`.

Because `H_1(T^2)` is a two-dimensional symplectic vector space, there are only three cases.

### Rank 0

`V=0`, so `V^perp=H_1(T^2)`:

```text
s=0,
s_perp=2.
```

### Rank 1

Every one-dimensional subspace is Lagrangian on the torus, hence

```text
V=V^perp,
s=s_perp=0.
```

### Rank 2

`V=H_1(T^2)`, `V^perp=0`:

```text
s=2,
s_perp=0.
```

Therefore, exactly,

```text
boxed:
r(H)-1 = [s(H)-s_perp(H)]/2.
```

No scaling limit or percolation input enters this identity.

## 3. The Matching-One topological source is a direct specialization

The bounded rank source used throughout the repository is

```text
exp[h (r-1)].
```

Using the identity above,

```text
exp[h(r-1)]
 = exp[h s/2] exp[-h s_perp/2].
```

Thus in the Krushkal topological monomial

```text
A^(s/2) B^(s_perp/2)
```

the exact rank-source specialization is simply

```text
boxed:
A=e^h,
B=e^-h.
```

It gives

```text
rank 0 -> e^-h,
rank 1 -> 1,
rank 2 -> e^+h,
```

which is exactly the repository's topological charge source.

The derivative at `h=0` is the observable

```text
X=r-1.
```

## 4. Refined FK partition function

For an embedded FK/bond configuration one may therefore define the topology-refined random-cluster state sum

```text
Z_G(Q,v;h)
 = sum_A v^|A| Q^k(A) exp[h(r(A)-1)].
```

The extra topological factor is exactly the Krushkal `A/B` monomial specialization above.  The ordinary random-cluster/Potts weights occupy the usual connectivity/edge variables; the rank source is an independent surface-topology refinement.

This gives a clean finite definition for the continuum object sought in #782 before any TBA machinery is invoked.

## 5. Combining rank charge with the rank-one neutral fugacity

The affine-TL seam result already gives, inside a fixed rank-one slope sector with `K` parallel essential FK clusters,

```text
z^K,
z=alpha^2/Q.
```

The Krushkal source and the seam source are complementary:

```text
Krushkal A/B variables : distinguish rank 0 versus rank 2;
affine-TL alpha seam   : resolves neutral count K inside rank 1.
```

Thus the desired schematic torus source

```text
T(q0,q2,alpha)
 = q0 Z0 + sum_u Z1,u(alpha^2/Q) + q2 Z2
```

has a natural finite-state interpretation with

```text
q0=e^-h,
q2=e^+h,
z=alpha^2/Q.
```

At `Q=1`, this is precisely the `(charge fugacity, neutral-count fugacity)` coordinate system already used on #771.

## 6. What remains genuinely missing for #782

This observation does **not** yet produce the massive torus scaling function.  The remaining tasks are now more precise:

1. construct a transfer / modified trace realization of the Krushkal `A/B` topology variables in the periodic loop/TL language;
2. combine it consistently with the noncontractible-loop seam `alpha`;
3. identify the corresponding sectors in the massive Potts finite-volume theory and continue toward `Q->1`;
4. verify that the UV limit reproduces the Pinson/Arguin homology weights.

So the verdict should be sharpened from

```text
"the rank0/rank2 projector is unknown"
```

to

```text
FINITE_TOPOLOGICAL_PROJECTOR_EXISTS;
MASSIVE_TRACE_REALIZATION_OPEN.
```

## 7. Why this also matters for the sector-odd correction (#768)

The matching-root observable is the response of the **signed topological source** `X=r-1`.  Therefore the first noncommon correction can be defined without first naming a CFT field:

> it is the leading irrelevant contribution to the `h`-odd part of the topology-refined torus free energy at `h=0`.

If `g` is a microscopic anisotropy/source parameter, the object of interest is schematically the mixed response

```text
partial_g partial_h log Z_G(Q,v;h) |_(h=0)
```

or its sector/free-energy analogue after the appropriate thermodynamic projection.

This supplies a precise map/topology label to the spin-four versus scalar `x=21/4` question.  The continuum operator must live in the channel selected by the Krushkal rank source, not merely have a compatible scaling dimension.

In particular, a scalar eight-arm field of the right dimension matters for the matching root only if it has nonzero matrix element in this `h`-odd topological channel.

## 8. Connection to map-resolved torus tomography

The Krushkal refinement gives a concrete combinatorial-map/topology coordinate rather than an arbitrary modular basis vector.  A future map-resolved torus solution space can therefore be organized by

```text
(rank charge h,
 rank-one slope u,
 neutral fugacity alpha,
 local operator/spin label).
```

This is closer to the repository's exact lattice semantics than fitting isolated modular functions without a source dictionary.

## 9. Site-model boundary

The exact derivation above is for spanning subgraphs / bond-FK configurations on a graph embedded in the torus.  Matching One's primary numerical model is square **site** percolation with the matching graph.

The site model already has the exact digital-Alexander rank variable `r` and source `exp[h(r-1)]`, so the observable itself is unambiguous.  What is not asserted here is that the square-site transfer is literally a Krushkal polynomial specialization without a decorated/incidence-graph construction.

For the massive Potts/FK continuum route of #782, however, the bond/FK formulation is exactly the relevant setting.

## 10. Literature interface

Krushkal's graph-on-surface polynomial introduces variables that record the genus of a regular neighborhood of a spanning subgraph and of its complement, with a duality motivated in part by Potts/Tutte statistical mechanics.  Equivalent information is carried by the Bollobas--Riordan polynomial for ribbon graphs.

The torus identity in Section 2 is a direct specialization of those definitions to the two-dimensional symplectic homology of `T^2`; it is not an additional literature theorem.

Primary references to inspect for the trace realization rather than the state-sum definition:

- V. Krushkal, *Graphs, links, and duality on surfaces*, arXiv:0903.5312;
- the Bollobas--Riordan / Krushkal ribbon-graph equivalence and subsequent transfer/algebra realizations;
- toroidal Potts/TQFT work where twisted sectors are required in partition functions.

## 11. Claim boundary

Exact:

```text
r-1=(s-s_perp)/2 on T^2,
A=e^h, B=e^-h realizes exp[h(r-1)] in the Krushkal topological monomial.
```

Strong interface conclusion:

```text
the finite rank0/rank2 topological projector already exists as a standard surface-graph state-sum coordinate.
```

Open:

```text
a local affine-TL/modified-trace realization compatible with massive Potts and Q->1,
and the resulting continuum matrix elements/scaling functions.
```

This turns an undefined-projector problem into a representation/transfer-realization problem.