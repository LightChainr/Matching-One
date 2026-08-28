# Generic-Q singlet filter for low-dimensional spin-4 competitors

Status: representation-theory safeguard for issue #37.

## 1. Formal Kac enumeration is not the physical Potts spectrum

For the critical Potts CFT at percolation (`c=0`), a non-diagonal module `W(r,s)` has conformal weights

\[
(h_{r,s},h_{r,-s}),
\qquad
h_{r,s}=\frac{(3r-2s)^2-1}{24}.
\]

The conformal spin is

\[
h_{r,s}-h_{r,-s}=-rs.
\]

Thus formal spin `-4` primary candidates obey `rs=4`.

The first three with integer `r>=2` are

\[
W(2,2):\quad x=17/4,
\]

\[
W(3,4/3):\quad x=196/27,
\]

\[
W(4,1):\quad x=49/4.
\]

If one looked only at Kac weights, `x=17/4` would seem to beat the thermal-family descendant `x=21/4`.  That conclusion is incomplete because physical Potts fields also carry the global `S_Q` representation.

## 2. Internal-symmetry data from the Potts space of states

Jacobsen--Ribault--Saleur, *Spaces of states of the two-dimensional O(n) and Potts models* (SciPost Phys. 14, 092 (2023), arXiv:2208.14298) determine the `S_Q` representation `Xi(r,s)` multiplying each non-diagonal conformal module.

Their exact formula implies

\[
\Xi(r,s)=\Xi(r,s+1),
\]

and their explicit examples give

\[
\Xi(2,0)=[2],
\]

\[
\Xi(3,1/3)=[21],
\]

\[
\Xi(4,0)=[4]+[2^2]+[21^2]+[3]+[21]+2[2]+[1]+\boxed{[]}. 
\]

Therefore

\[
\Xi(2,2)=\Xi(2,0)=[2],
\]

\[
\Xi(3,4/3)=\Xi(3,1/3)=[21],
\]

while

\[
\Xi(4,1)=\Xi(4,0)
\]

contains the trivial `S_Q` representation.

So, within this non-diagonal-primary spin-4 sequence, the first formal states at `x=17/4` and `x=196/27` are **not ordinary generic-Q singlet perturbations**.  The first one whose internal representation contains the generic-Q singlet is much higher, `x=49/4`.

See `scripts/generic_potts_singlet_spin4_gap.py` for the exact conformal arithmetic.

## 3. Why this strengthens x=21/4

The microscopic Potts/random-cluster action and the color-blind lattice deformation induced by changing local square-lattice couplings preserve the global Potts permutation symmetry.  A linear bulk irrelevant perturbation of the action should therefore lie in the generic-Q trivial representation.

The thermal primary is a singlet.  Its level-4 spin-4 quasiprimary descendant at

\[
x=x_t+4=5/4+4=21/4
\]

remains in the same singlet conformal family.  It therefore appears **below the first non-diagonal-primary singlet spin-4 candidate** identified above.

This removes the most obvious lower-dimensional formal Kac competitor to the ordinary lattice-perturbation interpretation.

It does not by itself prove that the matching observable couples to the thermal descendant: matching parity and the `Q=1` logarithmic module still have to be established.

## 4. The Q -> 1 logarithmic caveat

At `Q=1`, Potts representation theory becomes non-semisimple.  The energy operator is known to collide/mix logarithmically with a two-cluster field.  Therefore a sector that is nontrivial at generic Q cannot be dismissed from every color-blind percolation correlation function merely by quoting its generic-Q Young diagram.

The correct hierarchy is:

1. **ordinary lattice irrelevant perturbation:** require generic-Q singlet;
2. **special Q=1 logarithmic leakage:** allow non-singlet generic-Q sectors only through a demonstrated collision/Jordan mechanism;
3. test every such leaked candidate using its own fixed finite-size exponent/parity prediction.

In this sense `W(2,2)`, `x=17/4`, remains a controlled logarithmic-mixing competitor rather than the default operator.

## 5. Numerical discriminator already frozen before target data

`x=17/4` would give

\[
\Delta M\propto\Delta\cos4\theta\,N^{-9/8}
\]

and an angular root bias `L^-3`, whereas the thermal descendant gives `N^-13/8` and root `L^-4`.

The separate artifact

`predictions/x17_spin4_competitor_20260828.yaml`

freezes the `N^-9/8` predictions for the not-yet-produced issue #43 targets without modifying the original `N^-13/8` preregistration.

Scoring order remains:

1. original #43 `x=21/4` / `13/8` target;
2. fixed `x=17/4` / `9/8` competitor;
3. zero effect;
4. only then mixtures, logarithmic variants or free exponents.

This is the cleanest way to turn the representation-theory ambiguity into a falsifiable distinction.

## 6. Claim boundary

Supported:

- the `x=17/4` formal spin-4 Kac primary exists;
- at generic Q it carries `[2]`, not the trivial representation;
- the thermal-family `x=21/4` descendant is a generic-Q singlet candidate;
- the current data strongly prefer the `13/8` finite-size law retrospectively and prospectively through doubling.

Not proved:

- that no `Q=1` logarithmic collision feeds an `x=17/4` component into the matching observable;
- that the observed `13/8` field is uniquely the thermal quasiprimary rather than another singlet/logarithmic object with the same quantum numbers.

Those are now precise issues rather than hidden assumptions.
