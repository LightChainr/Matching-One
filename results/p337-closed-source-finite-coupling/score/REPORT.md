# Fixed finite-coupling closed-source turnover

Decision: `positive_source_global_U_turnover_resolved`. Negative exact U_t/A enclosures occur at
`m=[2, 4, 8, 16]`. The already published U(0)>0 and U_t(0)>0, together
with the unique analytic positive-coupling root, imply a positive-coupling
maximum before any resolved negative-derivative point.

| m=exp(t) | critical p | U | dU/dt |
|---:|---:|---:|---:|
| 2 | 0.710261724749 | 0.242368571338 | -1.37077822163 |
| 4 | 0.809483564856 | 0.000283785592447 | -0.00340200365326 |
| 8 | 0.890410580736 | 4.58979443475e-08 | -5.97094985292e-07 |
| 16 | 0.94139193809 | 3.97264799852e-12 | -5.70475721983e-11 |

These are four preselected finite-source laws, not a fitted scan.
The freeze is `b70dc4bd2fddd7676e9536b42bf912ee00ad302f`. Every prescribed value is reported; no grid
extension or new source is allowed by this decision. The baseline values
are imported from the prior exact packet without rescoring.

## Exact finite law and root-complete score

The integer histogram records (K,g,q), where
`g=2N+1-K-Sstar=2K-(beta1+beta_null)`.
With m=exp(t) and h=p/((1-p)m), its normalized weights are h^K/m^g.
Each geometry is normalized separately before the pooled matching mean
Q and P4(E). The pooled-root numerator has degree at most50. The root
uses the frozen 128 rational bisections; all reported enclosures are
outward rational arithmetic bounds, not sampling confidence intervals.

Let D=Q_h, B=Y_h, T=Q_hh, H=Y_hh, and jO=Cov(O,-g) within each geometry.
The four-term response is
`U_t/A=jY_h/D-H*jQ/D^2-B*jQ_h/D^2+B*T*jQ/D^3`.
The moving root is h_t=-jQ/D. Its common p-to-h Jacobian cancels in U;
using -g in this chart exactly differentiates the same homogeneous
Sstar curve. It is not a source substitution on an inhomogeneous
checkerboard saturation chart.

## Scope and next consequence

The result concerns the fixed N25 axis/tilted pair and the named source.
It resolves the monotone-amplification alternative if a downturn is
present; it does not identify a thermodynamic transition or continuum
field. The new full histogram enumerates the same finite configuration
populations as the old first-moment packet; it is not an independent
statistical vote. There were no new random samples, cloud jobs or test
campaign. The finite-volume strong-coupling argument predicts U->0;
the table is the frozen finite-coupling consequence, not a fit of its
asymptotic rate.
