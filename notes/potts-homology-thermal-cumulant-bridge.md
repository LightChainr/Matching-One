# Potts homology sectors and a thermal-cumulant representation of the matching scaling function

Status: long-horizon analytic bridge.  This note does not claim an exact value for `kappa3`; it identifies a concrete continuum object whose computation would determine it.

## 1. Restricted FK partition functions on a torus

For the Q-state Potts model in FK form,

\[
Z(Q,v)=\sum_{A\subseteq E} v^{b(A)}Q^{k(A)}.
\]

On a torus decompose this into homology classes `H`:

\[
Z=\sum_H Z_H.
\]

The classes include trivial/no-wrap, primitive one-dimensional winding subgroups `{a,b}`, and cross topology `Z x Z`.

Pinson computed the critical percolation (`Q=1`) probabilities, and Arguin generalized the restricted critical partition functions to Potts models.  In particular, the critical restricted weights obey a simple relation between trivial and cross sectors (at `Q=1` they have equal critical weights).

This is exactly the topology language already used by the repository's homology union-find.

References:

- H. T. Pinson, *Critical percolation on the torus*, J. Stat. Phys. 75 (1994).
- L.-P. Arguin, *Homology of Fortuin-Kasteleyn clusters of Potts models on the torus*, J. Stat. Phys. 109 (2002), arXiv:hep-th/0111193.

## 2. Cluster number is a Q derivative

Because Q is the cluster fugacity,

\[
\partial_{\log Q}\log Z=\langle k\rangle.
\]

Similarly, derivatives of the restricted weights `Z_H` generate cluster-number moments conditioned/weighted by homology sector.

This is the field-theory origin of the well-known use of Potts partition-function derivatives to obtain universal excess cluster numbers on a torus (Kleban--Ziff and related work).

Thus the two exact representations used in this repository,

- cluster-number/matching-polynomial side;
- wrapping/homology-probability side,

are naturally two derivatives/projections of the same FK partition-function object.

## 3. Thermal deformation

Let `t` be a continuum thermal coupling normalized so that the perturbation is

\[
S=S_{CFT}+t\int_{T^2}\epsilon(x)\,d^2x.
\]

For a restricted sector define formally

\[
Z_H(t)=Z_{CFT}\left\langle P_H\exp\left[-t\int\epsilon\right]\right\rangle,
\]

where `P_H` is the projector/indicator onto the torus homology sector.

Let

\[
E=\int_{T^2}\epsilon(x)\,d^2x.
\]

Then derivatives of the normalized sector probability

\[
P_H(t)=Z_H(t)/Z(t)
\]

are mixed connected cumulants of the topology indicator and integrated energy.  Schematically,

\[
P_H'(0)=-\kappa(P_H,E),
\]

\[
P_H'''(0)=-\kappa(P_H,E,E,E),
\]

with the standard cumulant subtraction implied by normalization.

The exact signs depend on the convention for `t`; ratios below do not depend on the normalization of `epsilon`.

## 4. Matching odd scaling function

In the scaling limit, the primary and matching lattices share the same universal restricted probability function at opposite thermal coordinates. For a cross-wrapping channel write

\[
P_c(z)=P_{cross}(z).
\]

Then

\[
\mathcal M(z)=P_c(z)-P_c(-z).
\]

Therefore

\[
\mathcal M'(0)=2P_c'(0),\qquad
\mathcal M'''(0)=2P_c'''(0).
\]

The dimensionless invariant becomes

\[
\boxed{
\kappa_3=
\frac{\mathcal M'''(0)}{\mathcal M'(0)^3}
=
\frac{P_c'''(0)}{4P_c'(0)^3}
}
\]

and hence, up to the sign convention for `t`,

\[
\boxed{
\kappa_3=
\frac{\kappa(P_c,E,E,E)}{4\,\kappa(P_c,E)^3}
}.
\]

The normalization of the thermal operator cancels exactly.  This is the continuum analogue of the repository's score/threshold-rank derivative estimators.

This formula is a more precise analytic target than asking whether a decimal happens to equal `-5/3`.

## 5. What is known and what is missing

Known ingredients:

1. exact critical torus homology-sector partition functions/characters for FK Potts;
2. the `Q->1` percolation limit;
3. the critical thermal conformal family (`h=hbar=5/8`);
4. integrable massive thermal perturbations of scaling Potts models and finite-size TBA/NLIE descriptions for many Q;
5. lattice transfer matrices already separate topological/connectivity sectors.

Missing ingredient:

> a practical representation of the **homology-projected finite-volume thermal deformation** whose first and third derivatives can be continued to `Q=1`.

Ground-state TBA alone is not enough.  We need a twisted/excited/defect sector corresponding to the cross/trivial homology projector, or a TCSA/transfer-matrix implementation that retains the same projector.

## 6. Three possible attack paths

### A. Conformal perturbation theory on the torus

Start from the exact critical homology partition functions and compute integrated one- and three-point energy insertions with the sector projector.

This is conceptually direct but technically difficult in a logarithmic `c=0` CFT.  Contact terms and logarithmic partners must be controlled.

### B. Scaling Potts TCSA

Construct the thermal perturbation in finite volume for Q values where the CFT/TCSA is clean, implement the topological/twist sector, compute the finite-volume crossover derivatives, and analytically/numerically continue toward Q=1.

This would give a continuum result independent of the square-site microscopic lattice.

### C. Temperley-Lieb / transfer-matrix sector bridge

The lattice FK transfer matrix naturally carries bridge/connectivity/topological sectors.  Identify which eigenvalue/character combination tends to the continuum homology projector, then study its thermal derivative before taking the scaling limit.

This path may connect most directly to the Jacobsen critical-polynomial machinery already relevant to the project.

## 7. Immediate numerical check enabled by current data

Before attempting a full continuum calculation, use the threshold-rank engine to measure the same mixed-cumulant ratios in several microscopic realizations at fixed torus modulus:

- square site + matching partner;
- a self-matching C4 site control;
- square bond self-dual control where the observable mapping is documented.

If the normalized topology/thermal cumulant ratios collapse, that gives a clean target for any CFT/TCSA computation.

## 8. Research boundary

The existence of exact critical homology probabilities does **not** imply that the off-critical scaling function is already known.  The thermal insertions are the hard part.

Likewise, the integrability of the scaling Potts perturbation does not by itself provide the required torus topology projector.  The value of this bridge is that it isolates precisely what must be solved rather than treating `kappa3` as numerology.
