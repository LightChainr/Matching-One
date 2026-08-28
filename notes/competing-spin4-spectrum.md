# Competing spin-4 spectrum: why x=17/4 must be tested but is not a generic-Q singlet perturbation

Status: operator-identification safeguard.  This note weakens any claim that `x=21/4` is automatically the first formal spin-4 Kac state.

## 1. A lower formal Kac spin-4 primary exists

At percolation (`c=0`) use

\[
h_{r,s}=\frac{(2r-3s)^2-1}{24}.
\]

For the non-diagonal pairing

\[
(h_{r,s},h_{r,-s})
\]

the conformal spin is

\[
h_{r,s}-h_{r,-s}=-rs.
\]

Setting `rs=4`, the smallest Potts-space integer `r>=2` candidate is

\[
r=2,\quad s=2,
\]

with

\[
(h_{2,2},h_{2,-2})=(1/8,33/8),
\]

so

\[
\boxed{x=17/4,\qquad s_{conf}=-4.}
\]

If such a field directly controlled the leading angular correction of a dimensionless torus observable, it would predict

\[
\Delta M\sim L^{-9/4}=N^{-9/8},
\]

and, after division by `M'~L^(3/4)`,

\[
\Delta p^*\sim L^{-3}.
\]

Thus formal Kac-table enumeration alone does **not** make `x=21/4` the lowest spin-4 number.

See `scripts/percolation_spin4_competitors.py`.

## 2. Potts internal symmetry changes the conclusion

The physical Potts space of states is not just a list of Virasoro weights.  Non-diagonal modules `W(r,s)` carry representations `Xi(r,s)` of the Potts permutation group `S_Q`.

Jacobsen--Ribault--Saleur (SciPost Phys. 14, 092 (2023), arXiv:2208.14298) show

\[
\Xi(r,s)=\Xi(r,s+1),
\]

because the interchiral field shifting `s` is an `S_Q` singlet, and explicitly give

\[
\Xi(2,0)=[2],\qquad \Xi(2,1/2)=[1^2].
\]

Therefore

\[
\boxed{\Xi(2,2)=\Xi(2,0)=[2].}
\]

The `x=17/4` primary is in the two-cluster/nontrivial internal-symmetry sector at generic Q; it is **not** the trivial representation `[]`.

This matters because the microscopic Potts/random-cluster action and the color-blind matching/wrapping construction preserve the Potts internal symmetry. A linear bulk irrelevant perturbation of the action must be an `S_Q` singlet at generic Q. Under that criterion `W(2,2)` cannot replace the singlet thermal-family level-4 candidate as the ordinary lattice anisotropy perturbation.

## 3. Why it cannot simply be deleted at Q=1

Percolation is a logarithmic `Q->1` limit.  The two-cluster representation is already known to collide with the energy sector at `Q=1`, producing logarithmic observables.  Thus it would be too strong to say that a field carrying `[2]` at generic Q can never enter a color-blind percolation observable after the singular limit.

The conservative statement is:

> `x=17/4` is excluded as an ordinary generic-Q singlet bulk perturbation, but remains a possible **Q=1 logarithmic-mixing competitor** until its amplitude is bounded or a representation-theoretic no-mixing argument is supplied.

Any such contribution has a distinct radial prediction `N^-9/8`, so it is experimentally separable from the `N^-13/8` thermal-descendant candidate.

## 4. Existing data already disfavor x=17/4

A retrospective one-amplitude fit to the pooled P31 same-N data gives, for the fixed H4 law,

```text
x=17/4  -> alpha_N=9/8:   chi-square about 13.78 / 4
x=21/4  -> alpha_N=13/8:  chi-square about  1.53 / 4
```

The exact numbers here are a diagnostic, not a preregistered model-selection result: the `x=17/4` alternative was identified after the P31/P32 targets were already visible.

The fresh prospective `(1+i)` doubling data also prefer the `13/8` ratio, but the strongest clean discrimination should come from targets frozen **after** the competitor was identified.

## 5. Prospective discrimination on N=185 and N=265

Issue #43 target data have not been produced at the time of this note.  Keep its original `13/8` H4 prediction artifact untouched, and add a separate frozen competitor artifact.

Using all existing pooled P31 source sizes to determine one `9/8` amplitude gives

\[
A_{9/8}\approx0.0865034,
\]

with source-only standard error about `0.00391265`.

For the already frozen #43 orientation pairs this implies approximately

```text
N=185: x=17/4 prediction  DeltaM = +2.86842e-4
       x=21/4 prediction  DeltaM = +1.92226e-4

N=265: x=17/4 prediction  DeltaM = +2.79919e-4
       x=21/4 prediction  DeltaM = +1.56735e-4
```

The separation is large compared with the source-amplitude uncertainties.  Final significance will be set by the target Monte Carlo variance.

Primary scoring order on #43 must remain the original frozen `13/8` prediction first.  The `9/8` model is a newly frozen **competitor**, not a reason to alter the original target.

## 6. Stronger theoretical selection question

The physical operator-identification problem is now sharper:

1. classify all spin-4 modules in the Potts/interchiral space below or near `x=21/4`;
2. attach their `S_Q` representations `Xi(r,s)`;
3. retain ordinary bulk lattice perturbations only from the generic-Q singlet sector;
4. separately classify `Q=1` logarithmic collisions that can leak non-singlet sectors into color-blind observables;
5. compare their fixed radial and matching-parity signatures against prospective data.

This is safer than selecting the first visually attractive Kac weight.

## References

- R. Couvreur, J. L. Jacobsen, R. Vasseur, *Non-scalar operators for the Potts model in arbitrary dimension*, arXiv:1704.02186.
- J. L. Jacobsen, S. Ribault, H. Saleur, *Spaces of states of the two-dimensional O(n) and Potts models*, arXiv:2208.14298, especially the `S_Q` decomposition and `Xi(2,0)=[2]`.
- H. Vasseur, J. L. Jacobsen, H. Saleur, *Indecomposability parameters in chiral logarithmic conformal field theory*, and related percolation-energy/two-cluster logarithmic-mixing work.
