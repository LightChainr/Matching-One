# Euler-critical insertion sites: a discrete-Morse view of the two-birth process

Date: 2026-09-14

Status: exact finite consequence of the square-site Euler identity and rank monotonicity.  The “Morse” terminology is descriptive; no smooth Morse-theory theorem is being imported.

## 1. Local Euler increment

For a monotone insertion permutation, just before a white site `v` turns black define

```text
c_b(v) : number of distinct black NN components touched by v,
c_w(v) : number of distinct white matching components adjacent through v
         after v is removed / equivalently merged when v is added in reverse,
d_b(v) : number of occupied NN neighbours of v,
f_b(v) : number of elementary squares completed all-black by adding v.
```

The configuration Euler identity gives exactly

```text
boxed:
Delta_v X
 = 1 - c_b(v) - c_w(v) + d_b(v) - f_b(v),            (1.1)
```

where

```text
X=r-1.
```

Because ambient homology rank is monotone under adding occupied sites,

```text
Delta_v X in {0,1,2}.                                (1.2)
```

Thus (1.1) implies the nontrivial digital-topology inequality

```text
1-c_b-c_w+d_b-f_b >=0
```

for every actual insertion state, although the separate terms have no obvious sign.

## 2. Only one or two sites in the whole permutation change ambient rank

The empty configuration has

```text
X=-1,
```

and the full configuration has

```text
X=+1.
```

Therefore along every strict insertion order

```text
sum_(v in insertion order) Delta_v X = 2.              (2.1)
```

Together with (1.2), there are only two possible topological histories:

### Two simple critical sites

```text
Delta X=1 at J1,
Delta X=1 at J2,
all other insertions Delta X=0.
```

This is the ordinary rank path

```text
0 -> 1 -> 2.
```

### One double critical site

```text
Delta X=2 at one insertion,
all other insertions Delta X=0.
```

This is the direct

```text
0 -> 2
```

birth with zero rank-one persistence gap.

Hence the entire two-birth process is exactly a process of one or two **Euler-critical insertion sites** inside an otherwise Euler-neutral permutation.

## 3. Forward/reverse union-find computes the whole rank path without lifted homology

For a fixed permutation `pi` of the torus sites:

### Forward black sweep

Add sites in the order `pi_1,...,pi_N` on the NN graph and store for every `k`:

```text
k4(k),
E(k),
F0(k),
K=k.
```

`k4` is updated by an ordinary union-find; `E,F0` are local motif counts.

### Reverse white sweep

Start from no white sites at `k=N` and add sites in reverse permutation order on the matching graph.  This gives

```text
k8_white(k)
```

for the complement of the first `k` black sites, again by an ordinary union-find with no homology gains.

Then at every cardinality

```text
boxed:
X_k
 = k4(k)-k8_white(k)-k+E(k)-F0(k),
r_k=X_k+1.                                             (3.1)
```

Thus the **entire rank path and birth times `(J1,J2)` can be recovered using two ordinary connectivity sweeps**, with the Euler identity supplying the homology rank.

This is an independent algorithmic route to threshold-rank trajectories on honest square tori.  It does not replace lifted homology when projective slope or actual winding class is required.

## 4. Direct birth is a local/global Euler imbalance, not a bare arm count

A `0->2` birth occurs exactly when

```text
1-c_b-c_w+d_b-f_b=2.                                  (4.1)
```

This formula clarifies why both theta/T3 and four-germ/rose geometries can realize the same rank jump.

For example, three black attachment germs can belong to **one** external black component, so `c_b` need not equal the number of visible black arms.  The local degree/germs, black-component merging, white-component splitting and completed-face term jointly determine the index.

Therefore

```text
arm count != Euler critical index.
```

Arm geometry is a continuum/local-separation refinement of the critical site, not its finite definition.

## 5. A natural classification for #769

For every critical insertion, store at minimum

```text
(rank_before, rank_after),
d_b,
f_b,
c_b,
c_w,
number/cyclic order of black attachment germs,
number/cyclic order of white separating germs,
black/white ambient lift data before/after.
```

The tuple `(d_b,f_b,c_b,c_w)` fixes the Euler index through (1.1); the germ/lift data distinguish theta, rose, split-component and small-torus degeneracies.

Then report topology classes separately for

```text
index 1 first birth,
index 1 second birth,
index 2 direct birth.
```

This is more informative than classifying only all `Delta_v X=2` configurations by a guessed arm label.

## 6. Relation to `alpha,beta,gamma`

At a fixed site / parameter, let

```text
alpha : 0->1 transition mass,
beta  : 0->2 transition mass,
gamma : 1->2 transition mass.
```

Then the local Euler index gives

```text
E[Delta_v X]
 = alpha + 2 beta + gamma.
```

The rank-one occupation derivative gives the signed combination

```text
P1'/N = alpha-gamma.
```

So the two natural lattice channels are

```text
activity/index channel : alpha+2beta+gamma,
birth-asymmetry channel: alpha-gamma.
```

A large theta/T3 contribution to `beta` is a statement about index-2 activity.  It does not by itself imply a large contribution to the signed birth-asymmetry channel relevant to the noncommon root correction.

This recovers the recent parity correction without assigning a continuum field sign.

## 7. Matching reflection of critical sites

Under complement and reversed insertion order, the first and second births exchange.  The set of Euler-critical sites is transported accordingly:

```text
index-1 first <-> index-1 second,
index-2 direct <-> index-2 direct.
```

This is the pathwise source of the matching-even persistence-gap/direct-birth law and the odd midpoint/birth-asymmetry coordinate.

Again this is a finite path relation, not an OPE automorphism.

## 8. Research consequence

The useful continuum question can now be phrased:

> what is the scaling geometry of the one or two Euler-critical insertion sites, and which signed combination of their local environments survives in the rank0/rank2 free-energy difference?

This is sharper than asking which arm event “is” the root correction.

Potential mechanisms become:

```text
theta/6-arm geometry dominates index-2 activity beta;
spin-four anisotropic correction lives in the difference between the two index-1 critical-site environments;
scalar/post-H4 corrections alter the common critical-site law or its higher connected responses.
```

These are falsifiable decompositions on actual lattice events.

## 9. Claim boundary

Exact:

- local Euler increment (1.1);
- nonnegative increments / total increment 2;
- one-or-two critical-site classification;
- forward/reverse connectivity reconstruction of the rank path.

Programme:

- use critical-site environment classes as the lattice interface to arm/CFT fields;
- identify which signed environment statistics control the observed H4 and post-H4 blocks.

The conceptual compression is substantial: a length-`N` rank process contains only one or two topologically active insertion events.