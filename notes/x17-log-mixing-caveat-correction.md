# Correction: generic [2] membership does not imply logarithmic mixing of W(2,2)

Status: correction to the conservative caveat in the spin-4 competitor notes.

## 1. The known percolation energy collision is dimension-specific

Vasseur--Jacobsen--Saleur, *Logarithmic observables in critical percolation* (arXiv:1206.2312), decompose two-spin Potts observables into a scalar energy-like field and a tensor field transforming in the nontrivial two-cluster representation.

The tensor creates two propagating FK clusters (the four-leg watermelon field).  Their Q->1 logarithmic field is constructed by mixing this tensor with the energy operator.

Crucially, the paper states that a finite Q->1 limit requires the two scaling dimensions to collide, and verifies

\[
\Delta_\epsilon(Q=1)=\Delta_{\hat\psi}(Q=1)=5/4.
\]

The logarithmic Jordan cell therefore comes from a **specific same-dimension collision**, not merely from the fact that one field carries a nontrivial `[2]`-type Potts representation.

## 2. Consequence for the formal W(2,2), x=17/4 competitor

The formal spin-4 module

\[
W(2,2):\quad (h,\bar h)=(1/8,33/8),\qquad x=17/4,\ s=-4
\]

has the same generic-Q internal representation class

\[
\Xi(2,2)=\Xi(2,0)=[2],
\]

but it does **not** share the energy dimension or spin.  The known energy/four-leg logarithmic collision at `x=5/4` therefore gives no direct mechanism for `W(2,2)` to leak into an `S_Q`-invariant spin-4 lattice perturbation.

To produce a logarithmic mixing at `x=17/4`, one would need a singlet field with the same scaling dimension and spin (or a more elaborate indecomposable collision with the required quantum numbers).  No such partner has been identified in the present audit.

Thus the correct theoretical prior is stronger than the earlier caveat:

> `W(2,2)` is excluded as an ordinary generic-Q singlet perturbation, and there is currently **no known Q=1 logarithmic collision that restores it** in the spin-4 singlet channel.

## 3. Why we still keep the frozen N^-9/8 numerical competitor

The artifact `predictions/x17_spin4_competitor_20260828.yaml` was frozen before N=185/265 target data and should not be edited retroactively.  It remains useful as a deliberately adversarial radial model:

\[
\Delta M\propto N^{-9/8}.
\]

Keeping and scoring it has methodological value even though its operator prior is now weak.  A numerical win for `9/8` would force us to discover a missing symmetry/logarithmic mechanism rather than assume one in advance.

The scoring order on #43 remains unchanged:

1. original `13/8` prediction;
2. frozen `9/8` adversarial competitor;
3. zero effect;
4. only then mixtures/free exponents.

## 4. Stronger current operator hierarchy

For an ordinary generic-Q `S_Q`-invariant square-lattice perturbation:

- `x=17/4` non-diagonal spin-4 primary: non-singlet `[2]`, no known same-dimension singlet collision;
- `x=196/27` next formal non-diagonal spin-4 primary: non-singlet `[21]`;
- `x=21/4` thermal-family level-4 spin-4 descendant: singlet because the thermal family is singlet;
- first simple non-diagonal spin-4 primary whose `Xi` contains the generic-Q singlet occurs only at `x=49/4`.

This makes the thermal-family `x=21/4` quasiprimary the lowest currently identified **generic-Q singlet spin-4 field above the dimension-4 sector** that naturally matches the observed `N^-13/8` law.

It is still not a proof of coupling or uniqueness: another singlet descendant/logarithmic field with the same quantum numbers could exist.

## References

- R. Vasseur, J. L. Jacobsen, H. Saleur, *Logarithmic observables in critical percolation*, arXiv:1206.2312.  Their Eqs. (7)--(11) explicitly use the collision `Delta_epsilon=Delta_psi=5/4` at Q=1.
- J. L. Jacobsen, S. Ribault, H. Saleur, *Spaces of states of the two-dimensional O(n) and Potts models*, arXiv:2208.14298, for the `Xi(r,s)` internal-symmetry assignments.
