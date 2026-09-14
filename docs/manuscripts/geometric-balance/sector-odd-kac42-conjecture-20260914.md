# Scalar Kac-(4,2) interpretation: demoted after rotational-symmetry check

2026-09-14.  CORRECTION / RETRACTION OF THE PRIMARY IDENTIFICATION IN THE EARLIER VERSION OF THIS NOTE.

The numerical fact that the square-site charge-sector mismatch is consistent with a correction dimension

\[
x_{odd}=21/4                                                    \tag{1}
\]

remains useful.  What is withdrawn is the earlier claim that the **most natural operator identification** is the scalar diagonal Kac field `(4,2)`, for which

\[
h_{4,2}=21/8,\qquad h+\bar h=21/4.                           \tag{2}
\]

That dimension match is real but insufficient.

## Why the scalar identification is now disfavored

Jacobsen's semi-infinite-cylinder data show two different lattice-symmetry patterns in the SAME percolation continuum theory:

- square-site percolation: the leading pseudo-critical correction has exponent `Delta_1=4`;
- kagome bond percolation: the `Delta=4` amplitude is absent and the leading correction is `Delta_2=6`.

Jacobsen explicitly suggests that the three-fold rotational symmetry of the kagome lattice, replacing the square lattice's four-fold symmetry, may force the first amplitude to vanish.

A scalar correction of dimension `21/4` is invariant under both `C4` and `C3/C6`.  Rotational symmetry therefore gives no natural reason for its amplitude to exist on the square lattice but vanish on the kagome lattice.  This is a substantive negative clue, not a minor aesthetic objection.

Accordingly the scalar `(4,2)` field is no longer the leading hypothesis.

## What survives from the old note

The finite-size scaling arithmetic remains:

\[
x_t=5/4,\qquad
\Theta_w'(p_c)\asymp w^{-1/4},                                \tag{3}
\]

and the observed square-site root shift

\[
p_w^{ch}-p_c\asymp w^{-4}                                    \tag{4}
\]

implies that the first correction which is ODD under exchange of the two topological magnetic sectors has total scaling dimension

\[
\boxed{x_{odd}=x_t+4=21/4.}                                   \tag{5}
\]

The direct transfer diagnostic

\[
\Theta_w(p_c)\asymp w^{-17/4}                                \tag{6}
\]

is consistent with the same value.

The open question is therefore the **spin / lattice-symmetry representation** of that `x=21/4` correction, not its dimension alone.

## Replacement conjecture

The current leading hypothesis is recorded in

`sector-odd-spin4-anisotropy-20260914.md`:

> the first square-lattice sector-odd correction has spin four and total dimension `x_t+4=21/4`; `C4` symmetry permits it, while `C3/C6` symmetry forbids it.  The next allowed rotational correction has spin six and produces a pseudo-critical exponent six, matching the kagome pattern.

This is structurally more compatible with Jacobsen's square/kagome comparison.

## Status of the Kac coincidence

Equation (2) should now be treated only as a **dimension coincidence / secondary alternative**.  It could still become relevant if a logarithmic module mixes a scalar `(4,2)` state with the actual spinful lattice correction, but no such mechanism has been shown here.

Any future attempt to restore the scalar interpretation must answer the rotational-selection objection explicitly.

## Literature boundary

Primary comparison: J. L. Jacobsen, arXiv:1507.03027.  He finds the square-site exponent `4`, the kagome leading exponent `6`, and explicitly proposes rotational symmetry as a possible reason the kagome `Delta=4` amplitude vanishes.  He does not identify either correction with a specific CFT field.

This note therefore records a self-correction: numerical exponent matching alone was not enough to justify the scalar Kac-(4,2) assignment.
