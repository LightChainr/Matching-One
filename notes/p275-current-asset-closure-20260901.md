# #275 current-asset closure: what is and is not identifiable

Status: `UNIDENTIFIABLE_WITH_CURRENT_ASSETS` for vacuum/Ward versus
thermal-Q4/Jordan as mechanisms of the original square-site `U` response.
This is an analysis closure, not a statement that the two continuum mechanisms
are physically identical.  Issue #275 may remain the highest theory attention
line, but no further production is implied.

## Three decisions that are now closed

1. **Gaussian primitive-C3 physical-angle spin identification is closed.**
   The held-out N145 three-way gate at branch-only commit `016a9d69` gives
   `p_H0=.968628`, `p_H4=.000135506`, `p_H8=2.67966e-13` and
   `X=.99992606+/-.00766449`.  Commit `74f55006` then proves that multiplication
   by the Gaussian unit `i` exchanges the unoriented homology lines `ell0` and
   `ell1` on every quotient `Z[i]/(g)`.  Hence `P_ell0=P_ell1` at every finite
   norm and every bond probability, and the baseline-subtracted character is

   ```text
   z = omega * [(P_ell2-pi_ell2) - (P_ell0-pi_ell0)] in omega R.
   ```

   The N145 H0 line is a realization of an exact finite symmetry, not a new
   spin-zero local field.  The earlier N65 H8 selection was a two-model angle
   alias.  No additional Gaussian-ideal angle with this observer can identify
   local H4, H8 or another spin.

2. **The current rho-child production excludes three simple probability
   models.**  On the one N112 2M/100-batch dependency block,
   common-normalizer-only, rank-1-mass-only and independent-real-rescaling
   parameterizations are excluded.  The observed correction changes internal
   winding-resolved rank-1 composition.  This does not identify a field.

3. **The #537 third-size route is stopped.**  The independent N145 200M full-T
   gate gives `J145=-.0006943644+/-.0001106636`, with frozen 95% interval
   `[-.0009112651,-.0004774638]` crossing the boundary `-.0004931`.  Its literal
   verdict is `N145_FULL_T_UNRESOLVED`, with no top-up.  Fixed-power comparisons
   remain post-reveal two-point fingerprints and cannot close #537.

## Candidate-to-observable readiness matrix

The original response requires one candidate-specific, same-source aligned
three-sector jet

```text
J_B = {s_r, partial_p s_r}_{r=0,1,2},   s_r=B_r/Z_r,
```

followed by the original `q/E` map, geometry projection, physical partition
normalizer and pooled moving root.  The current assets cover the following
parts:

| asset or candidate | critical sector/shape | candidate `B_r/Z_r` | same-source `partial_p(B_r/Z_r)` | rank-1 denominator | pooled-root original-U column | closure use |
|---|---:|---:|---:|---:|---:|---|
| vacuum/KdV P231 sector response | yes | no | no | no | no | continuum shape and semantic control only |
| thermal-Q4/Jordan Ward/Hecke response | yes, conditional | no | no | no | no | exact modular ratio, not a lattice-U prediction |
| Gaussian primitive-C3 N65/N145 | finite character only | no | no | observer-normalized | no | spin gate closed exactly by unit rotation |
| rho-child N112 covariance | fixed-p rank probabilities | no | no | observed at fixed p | no | can score a supplied restricted-sector column |
| global K1/K2 archive | baseline p-curves and U residual | no candidate insertion | no candidate insertion | baseline only | observed target | can score a supplied source image |

Neither named continuum candidate currently supplies the left-hand physical
source column.  Consequently the candidate images `C_V` and `C_T` in the
original-U data space are not defined.  A covariance-weighted profile-rank
comparison cannot be performed without silently inventing source amplitudes,
thermal derivatives or normalizer terms.

## Literal current verdict

```text
vacuum/Ward versus thermal-Q4/Jordan for original U:
UNIDENTIFIABLE_WITH_CURRENT_ASSETS

reason:
candidate-specific source--thermal restricted-trace columns are absent
```

This verdict is stronger than “underpowered” and narrower than “the mechanisms
are equivalent.”  More samples of the existing coordinates cannot repair an
undefined forward map.

## The only mechanism-changing reopening condition

Theory must deliver two actual candidate columns, each specifying

- `s_0,s_1,s_2` and `partial_p s_0,partial_p s_1,partial_p s_2` for one named
  microscopic source or a theorem fixing them;
- the rank-1 contribution and physical partition normalizer;
- common versus independent cross-geometry amplitudes and phases;
- the pooled-root counterterm and the exact map to original `q/E/U`.

Then compute the two column spaces with all declared nuisance amplitudes:

- if they coincide, retain `UNIDENTIFIABLE` and record the one missing physical
  relation;
- if they differ, perform one existing-covariance score and demote the failed
  parameterization;
- acquire a new coordinate only when the column-space calculation proves that
  exactly that coordinate is missing.

No new Gaussian angle, contact descriptor, generic certificate, size top-up or
free radial exponent is a substitute for these columns.  Priority allocates
attention; it does not lock or close parallel work.
