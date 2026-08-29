# Pinson--Arguin primitive sectors under the c=0 vacuum KdV K4 response

Status: deterministic continuum calculation plus a retrospective covariance
diagnostic on the PR #222 samples.  This is not a preregistered evidence row.

## The numerator is fixed at Q=1

For a finite random-cluster graph,

```text
Z_RC(Q=1) = sum_A p^|A| (1-p)^(|E|-|A|) = 1.
```

Consequently, after the common shape-independent bulk normalization, the
Pinson--Arguin primitive-sector probability `P_a(tau)` is the restricted
numerator `Z_a(tau)` in the Q=1 probability normalization.  There is no
unknown modulus-dependent total factor to reconstruct.  In particular,

```text
K4[Z_a]/Z - P_a K4[Z]/Z = K4[P_a],
K4 = D2 D0 = (delta-E2/6) delta,
delta = q d/dq = (2*pi*i)^(-1) partial_tau.
```

This Q=1 identity removes the normalization obstruction that would be present
if only `P_a=Z_a/Z` were known for a nonconstant `Z`.

## Analytic Gaussian-series derivative

With the PR #213 paper convention, put

```text
y = Im(tau),                 u = a-b*tau,
A_k = 2*pi*k^2 |u|^2/(3y),  T_k = 2 C c_k exp(-A_k),
C = sqrt(2/(3y))/|eta|^2,
c_k = cos(2*pi*k/3)-cos(pi*k).
```

Holding `tau_bar` fixed gives

```text
ell_k := delta log(T_k)
       = 1/(8*pi*y) - E2/24 - k^2 (a-b*tau_bar)^2/(6y^2),

delta ell_k
       = 1/(32*pi^2*y^2) + (E4-E2^2)/288
         - k^2 (a-b*tau_bar)^2/(12*pi*y^3).
```

Therefore the chiral response is the explicit convergent series

```text
K4[Z_ab] = sum_(k>=1) T_k
  [ell_k^2 + delta ell_k - (E2/6) ell_k].
```

The independent oracle complexifies the original numerator to
`Z_ab(tau,tau_bar)`, holds `tau_bar` fixed, and evaluates first and second
Wirtinger derivatives numerically.  Across the committed square, hexagonal,
N30 and N56 controls, the two routes differ by at most `9e-96`; an unrelated
complex modulus test is required to agree below `1e-58`.

## Why the real response is 2 Re, rather than a guess

For a fixed sector label, the Pinson--Arguin numerator is real on the physical
slice `tau_bar=conj(tau)`.  Hence its anti-chiral K4 action is the complex
conjugate of the chiral action.  A real reflection-even coupling aligned with
the registered `omega1` frame therefore produces

```text
R_a(tau) = (K4 + K4bar) Z_a = 2 Re K4[Z_a].
```

This statement holds at a fixed frame and fixed label.  A modular or period
basis transformation must transport both the sector labels and the spin-4
coupling phase; it is not licensed to apply `2Re` after silently relabeling
the sectors.

## N30/N56 deterministic vector

The modular calculation sets `omega1=1`, so the dimension-four perturbation
must be restored with `|omega1|^(2-4)=|omega1|^-2`.  Using `N^-1` would mix
the known shape Jacobian `Im(tau)` into the response because
`N=|omega1|^2 Im(tau)`.  In the registered `(C,Q,S)` basis the correct two
unit-coupling vectors are

| design | C/abs(omega1)^2 | Q/abs(omega1)^2 | S/abs(omega1)^2 |
|---|---:|---:|---:|
| N30, tau=1/2+5i/6, abs(omega1)^2=36 | 5.30237293437e-4 | 0 | -5.90897363014e-5 |
| N56, tau=1/2+7i/8, abs(omega1)^2=64 | 2.66358840128e-4 | approximately 0 | 8.03840950864e-6 |

Two parameter-free consequences are immediate:

```text
C30/C56 = 1.99068779989627 and has the same sign,
Q30 = Q56 = 0 by reflection.
```

The observed PR #222 C ratio is `1.99360563608899`, also same-sign.  Thus the
proper normalized sector KdV response does not inherit the naive scalar-E4
sign reversal used in PR #235.  PR #235 excludes that naive bridge, but not
this inhomogeneous sector-valued KdV direction.

## One-amplitude retrospective score

Using the full PR #222 block covariance and fitting a single common `g4` gives

| score | g4 | chi-square / df | survival p |
|---|---:|---:|---:|
| all N30/N56 C,Q,S coordinates | 12.2597 +/- 1.6384 | 52.7494 / 5 | 3.79e-10 |
| non-scalar C coordinates only | 14.2323 +/- 1.6703 | 0.0000258 / 1 | 0.99595 |

The full one-field sector vector fails because both measured S residuals are
large and positive, whereas the KdV subset-S response is small and changes
sign.  The non-scalar C direction, including its N30-to-N56 transport, is
nevertheless highly compatible.  The narrow conclusion is therefore:

```text
PR #222 is not a pure one-amplitude vacuum-KdV C/Q/S response;
its non-scalar C component remains vacuum-KdV compatible in the presence
of an independent scalar contaminant.
```

No new Monte Carlo was run, and neither covariance score is independent of
PR #222.
