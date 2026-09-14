# Euler charge tangent versus the common Potts `Q` tangent

Date: 2026-09-14

Status: source-direction correction to the generic-Q/log-collision narrative.  The finite Euler-source decomposition is exact; continuum branch projections remain to be derived.

## 1. The source directions are not automatically the same

For square-site Matching One, the exact Euler source is

```text
X=r4-1=k4-k8-K+E-F0.
```

Exponentiating gives

```text
e^{hX}
 = e^{+h k4}
   e^{-h k8}
   e^{-h K+hE-hF0}.
```

Thus, at the level of black/white cluster fugacities, the `h` tangent contains

```text
boxed:
(delta log Q_black, delta log Q_white)
 = (+1,-1) delta h,
```

plus explicit local site/edge/plaquette counterterms.

This is an **antisymmetric two-colour cluster-fugacity tangent**.

By contrast, the ordinary generic Potts `Q` derivative changes the common cluster fugacity of a declared random-cluster/Potts model.  In a two-colour extension its cluster-fugacity component is naturally common/symmetric rather than the Euler difference direction.

Therefore

```text
Euler h tangent != generic common-Q tangent
```

unless an additional model-specific map proves otherwise.

## 2. Why this matters for V14/W22 collision

The exact generic-Q Kac collision

```text
x_Kac(Q)-x_W(Q)
 = c_Q (Q-1)+...
```

with

```text
c_Q=-15/(2 pi sqrt(3))
```

is a statement about branch splitting under the **standard Potts Q deformation**.

If an observable is a Q derivative of a colliding pair, pole-cancelled amplitudes can generate a logarithmic size term through

```text
partial_Q L^-x(Q)
 = -x_Q' log L * L^-x + ... .
```

But the Matching-One root is an `h` response of the Euler-sourced paired site family.  The Q-splitting velocity only predicts its log coefficient if the actual `h` tangent has a nonzero projection onto the same generic-Q branch-splitting direction.

That projection has not been established.

Therefore the current V14/W22 note should be read as

```text
representation-theory resonance + possible generic-Q log mechanism,
```

not as a direct prediction that the Matching-One Euler source must show the corresponding log.

## 3. A two-colour deformation space is the natural missing object

Introduce formal local coordinates near percolation such as

```text
q_b = log Q_black,
q_w = log Q_white,
t   = thermal coordinate,
... local counterterm coordinates ...
```

and define

```text
q_common = (q_b+q_w)/2,
q_diff   = (q_b-q_w)/2.
```

Then the Euler source has

```text
partial_h q_diff = 1,
partial_h q_common = 0,
```

before the declared local counterterms are included.

The standard Potts Q tangent is primarily a `q_common` direction, with its own thermal/representation compensation.

The continuum RG problem is therefore to determine the tangent map

```text
(q_common,q_diff,t,...)
 -> scaling-field/block couplings.
```

The V14/W22 branch splitting supplies information along `q_common`; Matching One needs the `q_diff` column.

## 4. Bond/FK canonical coordinates provide an analogy, not the site proof

For toroidal FK the earlier exact canonical coordinates were

```text
x=v/sqrt(Q),
eta=h+(1/2)log Q.
```

There, changing Q at fixed `eta=0` requires a compensating `h` shift, while `eta` is an independent topological-charge direction.  This already shows algebraically that “Q” and “topological charge” are distinct tangent coordinates.

The square-site Euler source has a parallel conceptual structure, but its site/edge/face counterterms make the precise generic extension different.  Do not import the FK coordinate formula as an exact site identity.

## 5. A targeted generic-extension test

Before a large generic-Q production, construct the smallest finite two-colour extension in which:

1. black NN clusters receive fugacity `Q_b`;
2. white matching clusters receive fugacity `Q_w`;
3. local site/edge/plaquette weights are explicit;
4. `Q_b=Q_w=1` recovers ordinary site percolation;
5. the Euler source path is exactly reproduced by

```text
Q_b=e^h,
Q_w=e^-h,
```

with the local Euler counterterms.

Then compute a small response matrix of the candidate post-H4 block to

```text
partial_(q_common),
partial_(q_diff),
partial_t.
```

This is much more informative than differentiating only standard Q and assuming the answer transports.

## 6. Possible outcomes for the log collision

### Common-Q only

If the collision/Jordan logarithmic residue couples only to `q_common` while the Euler `q_diff` column is regular/small, the V14/W22 collision is structurally real but mostly irrelevant to Matching-One root response.

### Shared block

If both common and difference tangents project onto the colliding block, the generic-Q splitting velocity can be used after multiplying by the measured tangent-projection coefficient.

### Difference-specific structure

The Euler tangent may excite a different partner/combination than the standard Potts Q tangent.  Then the correct log coefficient must be derived from the two-colour extension, not inherited from `d_Q Delta x`.

## 7. Relation to the source-normalization audit

The local counterterms `-K+E-F0` are not optional decorations.  They are part of the exact Euler source.  Dropping them changes the tangent direction, just as omitting the row normalizer in #802 changed the physical source response.

Therefore any generic extension must match the **whole source path**, not only the cluster-fugacity signs.

## 8. Claim boundary

Exact:

- Euler source decomposition into `(+ black cluster, - white cluster)` plus local counterterms;
- distinction between common and difference fugacity tangent coordinates.

Known external algebra:

- V14/W22 generic-Q dimension collision and Q splitting velocity.

Open:

- RG projection of the Euler difference tangent onto that colliding block;
- whether the Matching-One post-H4 residual therefore carries a compulsory/visible logarithmic partner.

The main correction is simple: **a generic-Q collision is only relevant after the physical source tangent has been projected onto it.**