# Near-critical rank flux and suppression of direct rank-two births

2026-09-14.

Status: exact finite monotone-flux algebra plus a conditional near-critical scaling prediction. The key missing theorem is a macroscopic arm certificate for a one-site `rank 0 -> rank 2` jump in square/triangular percolation. If that certificate is established, the vanishing of the direct-jump rate follows from standard arm exponents/pivotal scaling.

## 1. Finite monotone rank flux

Under the common uniform-label coupling, insert occupied sites monotonically. Let the ambient homology rank be

```text
r in {0,1,2}.
```

At a single insertion, monotonicity permits

```text
0->0, 0->1, 0->2,
1->1, 1->2,
2->2.
```

For a fixed site and parameter, write the transition probabilities/intensities

```text
alpha = 0->1,
beta  = 0->2,
gamma = 1->2.
```

The exact boundary calculus already developed in the repository gives

```text
-P0' = N(alpha+beta),
 P1' = N(alpha-gamma),
 P2' = N(beta+gamma).
```

Only two of these equations are independent because `P0+P1+P2=1`. The direct-jump channel `beta` is therefore hidden from one-time rank probabilities.

This is the finite predecessor of the continuum flux equation below.

## 2. Near-critical scaling of the rank marginal

Put the torus/cylinder on macroscopic scale one and use

```text
p = pc + lambda * const * L^-3/4.
```

Assuming convergence of the rank law,

```text
Pj,L(lambda) -> Pi_j(lambda;tau).
```

Define continuum transition intensities per unit `lambda`

```text
a(lambda): 0->1,
b(lambda): 0->2,
c(lambda): 1->2.
```

Then the exact finite flux identities formally converge to

```text
-Pi0' = a+b,
 Pi1' = a-c,
 Pi2' = b+c.
```

This relation is independent of a massive-Potts representation and follows from the monotone common-label process if the scaling limits exist.

It also proves an information statement:

> the one-time functions `Pi0,Pi1,Pi2` do not determine the direct-jump channel `b`.

Thus a massive finite-volume engine for the rank marginals cannot by itself reconstruct the genealogy.

## 3. Why direct 0->2 jumps should be absent in the continuum pivotal clock

The lattice `0->2` event is the same one-step quantity called `jump2` / `B_k` in the rank-birth atlas. The local attachment analysis identifies a minimal `T3` geometry: closing the pivotal site leaves one exterior rank-zero component touching the site through at least three attachment germs, while opening the site creates two independent ambient homology generators.

To create two macroscopic homology generators at one microscopic insertion, the local branches cannot all merge inside a contractible disk well below the systole. The expected continuum certificate is therefore a polychromatic six-arm event from the pivotal neighbourhood to a fixed fraction of system size. Split/rose geometries require at least as much macroscopic structure and may be eight-arm.

Assume the certificate

```text
{Delta_v r=2} subset {six alternating arms from v to c L}
```

up to finitely many topology-equivalent variants.

For critical percolation the polychromatic arm exponents are

```text
alpha_4 = 5/4,
alpha_6 = 35/12,
alpha_6-alpha_4 = 5/3.
```

Ordinary near-critical time is normalized by four-arm pivotals. Therefore the fraction of rank-jump-two pivotals among macroscopic pivotals scales as

```text
L^[-(alpha_6-alpha_4)] = L^-5/3.
```

Equivalently, over an O(1) lambda interval,

```text
E[# direct 0->2 macroscopic jumps] = O(L^-5/3) -> 0.
```

Hence the continuum prediction is

```text
b(lambda)=0.
```

This recovers exactly the `R_jump2~L^-5/3` candidate that was independently proposed from the finite pivotal atlas, but here it has a process-level meaning.

## 4. Continuum rank process if b=0

The flux equations collapse to

```text
a(lambda) = -Pi0'(lambda),
c(lambda) =  Pi2'(lambda),
Pi1'      = -Pi0'-Pi2'.
```

At the self-dual/matching critical point, if

```text
Pi2(lambda;tau)=Pi0(-lambda;tau),
```

then

```text
a(0)=c(0).
```

So the aggregate rank projection has nearest-neighbour births only:

```text
0 -> 1 -> 2.
```

This is a statement about rank jumps, not about the full component partition.

## 5. Why this does not restore pure splitting

Macroscopic component mergers can occur in the near-critical pivotal process without producing `0->2` in one step. Examples include

- merger of two lineages in the same primitive homology class, leaving rank one;
- a rank-one component absorbing another lineage without changing ambient rank;
- ordinary `1->2` events driven by four-arm pivotals.

Therefore

```text
b=0
```

does **not** imply the merger factorial `J_tau(lambda1,lambda2)` vanishes. The full genealogy remains richer than the three-state rank process.

This distinction prevents a dangerous shortcut: absence of direct rank-two births cannot be used to infer the fixed-subcritical pure-splitting kernel at finite near-critical lambda.

## 6. Interface to the massive rank-source functions

If a massive-Potts/modified-trace construction supplies

```text
Pi0(lambda;tau), Pi1(lambda;tau), Pi2(lambda;tau),
```

then, conditional on `b=0`, it automatically predicts the continuum nearest-neighbour rank fluxes through derivatives.

This gives a new validation target for #782:

```text
-Pi0'(lambda) >= 0,
 Pi2'(lambda) >= 0,
 Pi1'=-Pi0'-Pi2',
```

with duality relating the first two across `lambda=0`.

But the full merger statistic remains an independent two-time observable.

## 7. What would prove or refute the claim

The decisive object is geometric, not another finite-size exponent fit:

1. prove that every one-site `rank 0->2` jump forces six alternating arms to macroscopic distance on the torus/cylinder;
2. identify any exceptional topology in which fewer arms suffice;
3. once the inclusion is proved, import the standard six-arm estimate in the rigorous triangular near-critical control;
4. for square site, treat the same conclusion as a universality target until the corresponding arm input is available.

A counterexample configuration with a genuine macroscopic `0->2` jump and fewer than six arms would invalidate the proposed continuum suppression mechanism immediately.
