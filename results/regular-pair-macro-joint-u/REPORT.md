# Canonical macro-window joint-U pilot

Status: **completed_valid_pilot_score**.

Decision: **stop the fixed-window field-ratio route at this pilot**.  Both
simultaneous intervals include zero, the two point estimates have opposite
signs, and the projected N400 requirement (2,119,100 configurations) exceeds
the frozen 2,000,000 ceiling.  Therefore D17/D21 are not evaluated, no top-up
or replacement window is allowed, and no local, GPU or cloud production is
launched.

The same fixed canonical Kreg, Euclidean 1/4--2/5 window and complete moving-root U functional are used at both sizes.
Anchors and displacement pairs are within-configuration readouts; the 100 paired batches are the inference units.

| N | root p | D | ESS min fraction | T_N total | simultaneous 95% interval | projected configurations |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.592893780109 | 9.8155809 | 0.999980 | +11.877846 ± 6.21 | [-2.2556162, +26.011308] | 1508400 |
| 400 | 0.59277438044 | 16.590459 | 0.999989 | -542.50382 ± 336 | [-1307.6425, +222.63486] | 2119100 |

## Complete original-U decomposition

| N | support | direct | root motion | source slope | root slope | J2_macro | T_N |
|---:|---|---:|---:|---:|---:|---:|---:|
| 100 | total | +0.0011611163 | +1.4302457e-05 | +1.5455835e-05 | -3.0899707e-06 | +0.0011877846 | +11.877846 |
| 100 | s2 | +0.0011199056 | +1.4200159e-05 | +1.547741e-05 | -3.0678699e-06 | +0.0011465153 | +11.465153 |
| 100 | sge3 | +4.1210719e-05 | +1.0229751e-07 | -2.1574382e-08 | -2.2100839e-08 | +4.1269341e-05 | +0.41269341 |
| 400 | total | -0.0046343615 | +0.0011797315 | +6.5666637e-05 | -1.685524e-06 | -0.0033906489 | -542.50382 |
| 400 | s2 | -0.0046208123 | +0.0011788392 | +6.5763406e-05 | -1.6842491e-06 | -0.0033778939 | -540.46303 |
| 400 | sge3 | -1.3549231e-05 | +8.9231858e-07 | -9.6768346e-08 | -1.274887e-09 | -1.2754956e-05 | -2.0407929 |

## Conditional fixed-power contrasts

Status: `not_evaluated_due_to_unresolved_or_sign_inconsistent_thermal_tail`.

D17/D21 are conditional pilot contrasts only. They do not estimate an exponent, accept/reject a field model or authorize production.
The s2 and sge3 rows are correlated additive coordinates of the same source, not independent evidence.
Raw Cbar diagnostics, all delete-one vectors and the complete covariance are stored in the JSON.
