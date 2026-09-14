# The bounded rank source is algebraically exhausted by `(M, chi)`

Date: 2026-09-14

Status: exact finite algebra and research-priority correction.

## 1. Three-state charge implies a closed source function

Let

```text
X=r-1 in {-1,0,+1},
M=E[X]=P2-P0,
chi=E[X^2]=P0+P2.
```

Then for every finite torus and every source `h`,

```text
Z(h)=E exp(hX)
    =P0 e^-h + P1 + P2 e^h
```

can be rewritten exactly as

```text
boxed:
Z(h)=1 + M sinh h + chi (cosh h-1).                  (1.1)
```

Therefore the complete `h` dependence at fixed microscopic parameters is determined by only two scalars `(M,chi)`.

At the balance root `M=0`,

```text
Z(h)=1+chi(cosh h-1).
```

## 2. Higher pure-h cumulants are not new observables

At the root,

```text
kappa_2(X)=chi,
kappa_3(X)=0,
kappa_4(X)=chi-3chi^2,
```

and every higher pure charge moment/cumulant is a polynomial in `chi`.

Away from the root they are polynomials in `(M,chi)`.

Hence a programme that repeatedly computes higher derivatives with respect to the bounded rank source `h` is algebraically redundant unless another parameter/source is varied simultaneously.

## 3. The informative objects are mixed susceptibilities

New information appears in derivatives such as

```text
partial_p M,
partial_g M,
partial_tau M,
partial_z M,
partial_g chi,
...
```

or equivalently mixed source derivatives

```text
partial_h partial_g log Z |_(h=0),
partial_h partial_p log Z |_(h=0),
```

because these probe how the exact Euler/rank defect couples to a physical perturbation.

For a normalized microscopic source score `H`,

```text
partial_g M = Cov(X,H).
```

This is exactly the numerator entering the root-response covariance quotient.

## 4. Correction to the master-source language

The combined finite object

```text
Z(p;h,z,tau)
 = P0 e^-h + P1 H_p(z) + P2 e^h
```

remains a useful bookkeeping object.  But its `h` coordinate is **not** an independently complicated functional direction: once sector weights are known, the `h` dependence is elementary.

Thus the hard continuum/massive tasks are really to determine

```text
P0(lambda,tau),
P2(lambda,tau),
rank-one neutral law H(z),
and mixed responses to physical sources,
```

not to solve an arbitrary `h`-dependent function from scratch.

A periodic modified trace realizing `h` can still be technically valuable because it extracts the sector sums; the novelty lies in computing the sectors, not in the exponential source algebra itself.

## 5. Research-priority consequence

Prefer

```text
one new mixed h-g response
```

over

```text
many higher h cumulants.
```

This is another instance of the research-compass rule: exact solvability of a generated family does not imply that every derivative carries independent information.

## 6. Claim boundary

Exact: all formulas above.

Interpretive: prioritize mixed source/geometry/thermal responses over pure topological-source elaboration.