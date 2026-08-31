# An Euler-invisible spatial tangent shifts the birth center and lifetime

The same-degree contact covariance has an exact constructive interpretation:
it is the derivative of a next-site selection rule that preserves the entire
immediate rank/Euler distribution. Its measured effect is nonzero for both
birth clocks. Moreover, E's response is **positive near p_ref but negative
after integration** at both N325 and N425.

This turns an attachment correlation into a precisely defined spatial
perturbation of the permutation ensemble. It does not require physically
editing a loop while holding every other graph feature fixed.

## Exact finite-prefix perturbation

Fix an original R0 prefix Z with d vacant labels. Let

```
A_e = {u : u preserves R0 and has occupied contact degree e},
pi_e = |A_e|/d,
g(u) = e-c(u).
```

Preserve the probability pi_e of each safe degree class, but change the
conditional choice within it to

\[
q_t(u\mid Z)=\pi_e\,
 \frac{\exp\{t\pi_e g(u)\}}
      {\sum_{v\in A_e}\exp\{t\pi_e g(v)\}},\qquad u\in A_e.
\]

All labels outside the safe classes retain probability 1/d. Empty classes
are omitted; non-R0 prefixes are left unchanged. After the selected label,
use the original uniform remaining suffix law. At t=0 the label law is
uniform over all d sites.

For every finite t, each safe class has exactly its original mass. Its
sites have the same rank outcome and Euler increment `1-e`; all other
label probabilities are unchanged. Consequently the joint distribution of
the immediate rank and Euler increment is **exactly unchanged**, not just
unchanged to first order. The original prefix and occupied count also stay
fixed.

For any future response Y with conditional mean m_Y(Z,u), differentiation
gives

\[
\left.\frac{d E_t[Y\mid Z]}{dt}\right|_{t=0}
 =\sum_e\pi_e^2\operatorname{Cov}(g,m_Y\mid Z,A_e).
\]

That is precisely the existing equal-degree half-difference numerator H_Y.
The factor pi_e in the exponent is deliberate: it matches the measured
pi_e-squared weighting without estimating a small conditional frequency
from sixteen next-label draws. This defines a concrete perturbation, not
the only possible normalization of the tangent.

## Measured response from saved sufficient statistics

The following derivatives average the two **orientation-specific** rules
equally over the original prefix populations. They are not one common-label
H4 policy; the two geometries have different safe degree classes.

| Response derivative at t=0 | N325 +/- original-batch SE | N425 +/- original-batch SE |
|---|---:|---:|
| E[K1] | .01035703 +/- .00080415 | .01362344 +/- .00128206 |
| E[K2] | .01553516 +/- .00105874 | .01851172 +/- .00191199 |
| E[C] | .01294609 +/- .00079131 | .01606758 +/- .00143064 |
| E[W] | .00517813 +/- .00101514 | .00488828 +/- .00155301 |
| A(p_ref) | -.000410330 +/- .000024542 | -.000434664 +/- .000042981 |
| E(p_ref) | +.0000829860 +/- .0000223681 | +.0001027765 +/- .0000164317 |
| integral A | -.0000794239 +/- .0000048546 | -.0000754346 +/- .0000067166 |
| integral E | -.0000158838 +/- .0000031139 | -.0000114748 +/- .0000036456 |

The integral identities are exact:

```
d integral A / dt = -2 d E[C]/dt /(N+1),
d integral E / dt = -d E[W]/dt /(N+1).
```

Thus the negative integral E response is the positive lifetime response,
not independent confirmation. The positive E response near p_ref has the
opposite sign. The measured tangent is therefore thermally redistributive;
if these population signs hold, continuity forces it to be negative
somewhere away from p_ref. No crossing location is inferred without reading
the full response curve.

This gives a finite response direction invisible to the immediate
rank/Euler distribution but visible to subsequent topology. It excludes
the corresponding next-response closure on those summaries; it does not
identify a continuum field or show that the unperturbed global H4 mean is
caused by this tangent.

## Sources and reproduction

Source is the saved raw contact covariance at
`b9f79bfb6e1ba4177ff245f74f7b2e51c3bd2fdc`, derived from the same original
e32a8593 conditional tails and 959a7fa2 contact marks. This eight-coordinate
projection took about0.06 seconds and did not read the long fork CSVs again.
`scripts/p334_euler_invisible_tangent.py` produces
`results/p334-euler-invisible-tangent/score.json`, including all twenty raw
batch vectors and their common covariance. Every coordinate lies in the
existing87-column contact covariance span; no new independent evidence
block, simulation, DP or validation campaign was created.
