# Complete F3 flat-twist curves from the N65 birth archive

The P334 sparse tuple is sufficient for the complete degree-65 Bernstein curves `T_alpha(p)`, their analytic derivatives, all three projective characters and arbitrary-p covariance. No new sample or event field is needed.

All microcanonical partition, D4 and complement/Alexander transport gates pass. The maximum reported residual is `2.5e-15`.

## Twist-sector crossing selector

The central window is `[0.45,0.75]`. Candidate roots are scored by whether the same bracket survives every leave-one-batch reconstruction.

| source | character | root | LOO survival |
|---|---|---:|---:|
| `first` | `H4_axis_diag` | -- | 0% |
| `first` | `axis_odd` | 0.485472926 | 60% |
| `first` | `diagonal_odd` | -- | 0% |
| `second` | `H4_axis_diag` | -- | 0% |
| `second` | `axis_odd` | 0.735437526 | 85% |
| `second` | `diagonal_odd` | 0.736887076 | 55% |
| `second_minus_first` | `H4_axis_diag` | 0.573633326 | 100% |
| `second_minus_first` | `axis_odd` | -- | 45% |
| `second_minus_first` | `diagonal_odd` | -- | 0% |

The unique 20/20-stable candidate is the parameter-free equality of the balanced F3 H4 characters between the two physical Gaussian orientations:

`p_cross=0.573633326`, derivative `0.0602462`, jackknife SE `0.02281`.

Exact complement transport gives the parameter-free dual partner `p_cross_dual=0.426366674` with derivative `-0.0602462`.

Its distance from `p_ref` is `0.01911`, smaller than one current root SE. Under simple inverse-sample scaling, approximately `113988` samples/shape would target 2 sigma and `256474` would target 3 sigma. This run does not request or perform that production.

## Interpretation

Zero twist stays exactly one. Nonzero twists are monotone constraint partition curves; zero-sum characters remove the common rank-zero curve. Complement preserves each line character with `p -> 1-p`, while D4 acts as `S:(H,A,D)->(H,-A,-D)` and reflection as `(H,A,D)->(H,A,-D)`.

The H crossing is a geometry-selector output from the reused 20k block, not a resolved physical root split. The coefficient archive and evaluator are the production-ready result; a future independent block can freeze this crossing without changing the curve model.
