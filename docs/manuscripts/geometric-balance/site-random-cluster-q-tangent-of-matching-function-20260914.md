# The Matching function is exactly a site-random-cluster `Q` tangent difference

Date: 2026-09-14

Status: exact finite identity plus a conditional continuum interface.  This note partially corrects `euler-source-vs-common-q-tangent-20260914.md`: the **full nonlinear Euler `h` source family** is not the ordinary common-`Q` deformation, but the first charge moment used by Matching One is nevertheless exactly a difference of ordinary cluster-fugacity `Q` derivatives.

## 1. Site random-cluster generating function

For a finite site graph `G` with `N` vertices, define

```text
Z_G(p,Q)
 = sum_omega
   p^{K(omega)} (1-p)^{N-K(omega)}
   Q^{k_G(omega)}.                                    (1.1)
```

This is the standard site random-cluster deformation of Bernoulli site percolation: `Q=1` gives the independent site measure, while `Q!=1` rewards or penalizes occupied clusters.

Because the Bernoulli weights are normalized,

```text
Z_G(p,1)=1                                             (1.2)
```

for every `p`.

Therefore

```text
boxed:
partial_(log Q) log Z_G(p,Q)|_(Q=1)
 = E_p[k_G].                                           (1.3)
```

## 2. Exact representation of the Matching-One observable

The finite square-site Euler identity is

```text
X=r4-1
 = k4(omega)-k8(omega^c)-K+E-F0.
```

Taking expectation under black density `p`, the complementary white matching configuration has the same law as site occupation with density `1-p` on `G8`.  Also

```text
E[K-E+F0]
 = N [p-2p^2+p^4].                                    (2.1)
```

Hence exactly

```text
boxed:
M_L(p)=P2-P0
 = partial_(log Q)
   [ log Z_G4(p,Q)
    -log Z_G8(1-p,Q) ]_(Q=1)
   -N[p-2p^2+p^4].                                    (2.2)
```

Thus the finite matching function is a **difference of two ordinary site-random-cluster cluster-number tangents**, with the known local Euler polynomial subtracting the bulk/local piece.

This is stronger than saying the full Euler `h` source looks like a two-colour cluster gas.  Equation (2.2) uses the standard one-colour common cluster fugacity `Q` on each graph separately.

## 3. Important correction: first `Q` derivative is automatically thermal-nuisance free at `Q=1`

Let `p_c(Q)` be any differentiable critical curve through `p_c(1)`.  Along that curve,

```text
d/d(log Q) log Z_G(p_c(Q),Q)|_1
 = partial_(log Q) log Z_G|_1
   + (dp_c/dlogQ)|_1 partial_p log Z_G(p,1).
```

But from (1.2),

```text
partial_p log Z_G(p,1)=0                               (3.1)
```

identically.

Therefore

```text
boxed:
[d/d(logQ) along any p(Q)] log Z_G|_(Q=1)
 = partial_(logQ) log Z_G|_(Q=1).                     (3.2)
```

So the **first cluster-number Q tangent does not need a thermal nuisance subtraction**.  The normalized `Q=1` partition function kills that chain-rule term automatically.

This is unlike a generic unnormalized microscopic source response, where pressure/thermal normalization must be handled explicitly.

At second and higher `Q` derivatives, mixed thermal terms re-enter; (3.2) is a first-derivative statement.

## 4. Relation to the Euler `h` source

The full topological source satisfies configurationwise

```text
e^{hX}
 = e^{h k4-h k8-hK+hE-hF0}.
```

That nonlinear family is not the same as changing one common `Q` in both site random-cluster models.

However its first derivative at `h=0` is exactly `E[X]=M`, and equation (2.2) gives an alternative representation of the same first moment.

Thus both statements are true:

```text
full h-source family != common-Q family,
first Matching charge M = difference of common-Q tangents + local Euler term.
```

The earlier note was too pessimistic when it used the first statement to weaken the generic-Q route for the root observable itself.

## 5. Continuum significance

A site random-cluster model with cluster fugacity `Q` is an established statistical model.  Existing numerical work (Wang et al., Phys. Rev. E 92, 022127 (2015), arXiv:1411.4408) reports thermal and magnetic critical exponents matching the bond random-cluster/Potts values for `Q=1.5,2,2.5,3,3.5,4`, with first-order behaviour beyond the usual range.

This is evidence—not a rigorous theorem—that the critical site-RC `Q` family realizes the same Potts/random-cluster continuum universality branch.

If that interface is accepted, then the Q-tangent representation (2.2) gives a much more direct route from Matching One to generic-Q CFT data than an abstract parity assignment.

## 6. Why logarithmic collisions become genuinely relevant again

Suppose a continuum contribution has generic-Q form

```text
A(Q) L^{-x(Q)}.
```

Then its cluster-number derivative contains

```text
partial_Q [A(Q)L^{-x(Q)}]_(Q=1)
 = L^{-x(1)}
   [A'(1)-A(1)x'(1) log L].                         (6.1)
```

Therefore a collision/splitting of dimensions as `Q->1` can generate logarithmic finite-size terms in precisely the kind of `Q` derivative that appears in (2.2).

This does **not** prove that the V14/W(2,2) collision controls the post-H4 residual: one still needs the actual amplitudes/map sectors and the difference between the G4 and G8 regularizations.  But generic-Q splitting is no longer merely an unrelated tangent direction.

## 7. The correct object is a difference of Q tangents between two microscopic regularizations

Write the critical finite-size/site-RC free energy schematically as

```text
log Z_G(p_c,Q)
 = N f_G(Q)
   + F_G^torus(Q,tau)
   + sum_j u_j^G(Q) L^{-omega_j(Q)} ... .
```

Equation (2.2) takes

```text
partial_(logQ) [G4-G8] at Q=1
```

and subtracts the local Euler polynomial.

The infinite-volume Sykes--Essam relation explains the cancellation of the extensive cluster-density part at `Q=1`.  The remaining torus/irrelevant difference is exactly where the root information lives.

Thus the useful continuum question is:

> which generic-Q finite-size blocks have different `Q`-tangent amplitudes in the NN and matching site-RC regularizations after the exact local Euler term is removed?

This is a source-defined alternative to calling fields “matching odd.”

## 8. A revised interpretation of the `x=21/4` and `x=33/4` candidates

For the leading H4 block, the observable lattice fact is that the `Q`-tangent difference between the two site regularizations first becomes visibly nonzero at effective total dimension `21/4`.

For the post-H4 scalar candidate, a generic-Q collision at `x=33/4` can generate a pure-power and/or logarithmic contribution to the same cluster-number tangent difference.

The correct hierarchy is therefore

```text
angular irrep first,
generic-Q Q-tangent block second,
regularization-difference amplitude third,
operator/Jordan naming last.
```

## 9. Practical finite test

A tiny finite-torus implementation can introduce `Q` as a symbolic/automatic-differentiation cluster fugacity independently for

```text
black NN cluster count,
white matching cluster count.
```

At `Q=1`, verify

```text
partial_logQ log Z4 = E k4,
partial_logQ log Z8 = E k8,
```

and reproduce the Euler rank observable through (2.2).

At generic `Q` one can then inspect the two cluster models separately without inventing a local continuum parity.  Large generic-Q production should wait until the relevant angular block/source dictionary is specified.

## 10. Claim boundary

Exact finite:

- site-RC generating function (1.1);
- cluster-number derivative (1.3);
- Matching function representation (2.2);
- first-derivative critical-line invariance (3.2).

External empirical interface:

- site random-cluster and bond random-cluster/Potts universality agreement is numerical evidence in the cited work, not proved here.

Open continuum:

- exact mapping of the square-site Q tangent onto the generic-Q interchiral/logarithmic blocks;
- nonzero amplitude of any particular V14/W22 collision in the G4-G8 difference;
- higher-Q derivative/contact structure.

The key correction is positive: **the root observable itself is already an ordinary site-random-cluster Q tangent difference, even though its full bounded `h` source family is not.**