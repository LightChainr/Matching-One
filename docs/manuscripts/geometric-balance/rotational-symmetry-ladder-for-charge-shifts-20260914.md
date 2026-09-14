# A rotational-symmetry ladder for semi-infinite charge-root shifts

2026-09-14.  Broad conjectural synthesis motivated by the square/kagome contrast in Jacobsen's eigenvalue method and the sector-odd thermal-family hypothesis developed on this branch.

## 1. General rule

Consider a two-dimensional percolation lattice whose continuum limit is the ordinary `c=0` percolation fixed point.  Let the semi-infinite-cylinder charge/eigenvalue root be defined by equality of the primal and dual magnetic/topological sectors.

Suppose:

1. lower-dimensional sector-even corrections cancel from the sector difference;
2. the leading sector-odd irrelevant perturbation belongs to the thermal conformal family;
3. the lowest nonredundant thermal-family anisotropy allowed by the microscopic point group has spin `s_*`.

A spin-`s` descendant has total dimension

\[
x_{odd}=x_t+s,
\qquad x_t=5/4.                                                \tag{1.1}
\]

The critical sector mismatch in a per-row excitation energy then scales as

\[
\Theta_w(p_c)\asymp w^{1-(x_t+s_*)},                         \tag{1.2}
\]

while the thermal slope scales as

\[
\Theta'_w(p_c)\asymp w^{1-x_t}.                              \tag{1.3}
\]

Their ratio gives the remarkably simple prediction

\[
\boxed{p_w-p_c\asymp w^{-s_*}.}                              \tag{1.4}
\]

Thus the pseudo-critical shift exponent directly reads the lowest allowed **sector-odd rotational spin**, not a generic scalar irrelevant exponent.

## 2. Point-group ladder

For a nonchiral lattice with reflection/inversion symmetry, odd spin is excluded.  If the rotational subgroup is `C_k`, the first candidate is the smallest positive even integer divisible by `k` (or equivalently by the effective bulk rotation order):

\[
s_*=\operatorname{lcm}(2,k).                                 \tag{2.1}
\]

This gives

| effective point symmetry | first allowed sector-odd anisotropy | predicted shift |
|---|---:|---:|
| `C2` | spin 2 | `w^-2` |
| `C4` | spin 4 | `w^-4` |
| `C3` plus inversion/reflection | spin 6 | `w^-6` |
| `C6` | spin 6 | `w^-6` |

Exact self-matching/self-duality can override the table by forcing the whole sector-odd amplitude to zero.

## 3. Existing evidence

### Square site

The square lattice has `C4` symmetry.  Jacobsen finds

\[
\Delta_1=4,                                                    \tag{3.1}
\]

and the safe transfer directly resolves

\[
\Theta_w(p_c)\asymp w^{-17/4},
\qquad
\Theta'_w\asymp w^{-1/4}.                                    \tag{3.2}
\]

This is exactly the `s_*=4` case.

### Kagome bond

The kagome lattice has hexagonal bulk symmetry (Jacobsen discusses the relevant basis contrast as three-fold versus the square lattice's four-fold rotation).  He finds the square-like exponent-four amplitude absent and

\[
\Delta_{lead}=6,                                               \tag{3.3}
\]

and explicitly suggests rotational symmetry as the reason.

This is exactly the `s_*=6` case.

### Triangular site

Triangular-site percolation is self-matching at `p_c=1/2`.  Digital Alexander plus complement symmetry gives the charge root exactly `1/2` at every finite size.  Thus all sector-odd amplitudes vanish, a stronger cancellation than the rotational rule alone.

## 4. The most discriminating new prediction: C2 should give exponent two

A non-self-dual percolation realization with only two-fold rotational symmetry but reflection/inversion retained should permit the spin-two thermal-family anisotropy.  The ladder therefore predicts

\[
\boxed{p_w-p_c\asymp w^{-2}.}                                \tag{4.1}
\]

At the free-energy level,

\[
\boxed{\Theta_w(p_c)\asymp w^{1-(x_t+2)}=w^{-9/4}.}           \tag{4.2}
\]

This is a much sharper falsification test than adding more square widths: exponent two is qualitatively separated from both four and six.

A useful implementation must avoid a fake `C2` produced only by choosing a rectangular simulation box for an otherwise `C4`-symmetric microscopic lattice.  The **local interaction/occupation rule** itself must break `C4` while remaining non-self-dual, and the continuum metric anisotropy must be accounted for rather than mistaken for an irrelevant spin-two coupling.

## 5. Metric anisotropy versus genuine spin-two irrelevant coupling

Two-fold anisotropy contains a potential trap.  A leading deformation of the continuum metric is a redundant/marginal geometric reparameterization, not the irrelevant descendant responsible for (4.1).  Therefore a clean `C2` test should:

1. determine the physical correlation-length metric;
2. express the cylinder circumference in that isotropized continuum metric;
3. only then fit the residual sector-odd charge mismatch.

The prediction `Delta=2` refers to the first **nonredundant** spin-two thermal-family correction after this metric calibration.

Failure to remove metric anisotropy can manufacture lower-order shape errors that have nothing to do with the charge-sector selection rule.

## 6. Controlled symmetry breaking gives amplitude selection rules

Let `a_s` denote a microscopic anisotropy coupling transforming in spin `s`.  Near a high-symmetry lattice, the charge mismatch should have schematic expansion

\[
\Theta_w(p_c)
=\sum_s C_s a_s\,w^{1-(x_t+s)}+\cdots,                        \tag{6.1}
\]

where only point-group-invariant combinations survive.

This gives several targeted experiments.

- Start from a six-fold non-self-dual lattice and add a weak spin-four distortion.  A `w^-17/4` sector mismatch and `w^-4` root shift should turn on linearly in the spin-four component.
- Start from square `C4` and add a weak `C2` distortion after metric calibration.  A new `w^-9/4` mismatch / `w^-2` root shift should eventually dominate the original spin-four term.
- Preserve self-matching while distorting geometry: the whole odd charge difference should still vanish, providing a parity control against pure metric artifacts.

## 7. Why this is not a generic correction-to-scaling statement

Ordinary percolation observables can have lower correction exponents and lattice-dependent analytic terms.  The eigenvalue/matching root is special because it subtracts two topological magnetic sectors that share the same leading CFT content.  The rotational ladder concerns the **first correction surviving that subtraction**, not the leading irrelevant operator of the bulk theory.

This is why a common `x≈4` correction can be plainly visible in each square-sector energy while the root drift begins only at exponent four through a much higher-dimensional sector-odd field.

## 8. Literature boundary

Jacobsen, arXiv:1507.03027, supplies the crucial square/kagome empirical pattern and explicitly raises rotational symmetry as the explanation for the missing kagome exponent-four amplitude.  Older finite-size-scaling work on square versus triangular/honeycomb critical models also documents lattice-symmetry-dependent cancellation of correction amplitudes.  The general idea that lattice rotations select allowed conformal spin is standard.

The specific formula

\[
\boxed{\Delta_{charge}=s_*}                                  \tag{8.1}
\]

for the primal/dual charge-sector eigenvalue criterion, and especially the `C2 -> Delta=2` prediction, are research conjectures here.

## 9. Claim boundary

The point-group selection rules are exact symmetry statements.  The mapping from the leading allowed spin to a thermal-family descendant of dimension `x_t+s` and hence to the pseudo-critical exponent `s` is conjectural.  It should be tested on a genuinely `C2`, non-self-matching model before being elevated beyond a mechanism hypothesis.
